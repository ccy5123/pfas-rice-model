"""
PFAS Rice Compartmental Uptake Model -- Plant Module (Method A: loose coupling)
==============================================================================

Solves the four-compartment (root / stem / leaf / fruit-grain) dynamic plant
uptake ODE system for a fully-dissociated PFAS anion, taking the soil pore-water
free concentration C_w^o(t) and transpiration Q_TP(t) from an EXTERNAL soil model
(e.g. HYDRUS-1D via Phydrus).  The plant-side equations follow the report
"Mechanistic Compartmental Model for PFAS Bioaccumulation in Rice".

Coupling (Method A, one-way):
    HYDRUS-1D / Phydrus  --->  C_w^o(t), Q_TP(t)  --->  this module
    (plant growth model) --->  M_k(t)             --->  this module

Equation map (report -> code):
    j_R  (Eq. JR_pfas)        -> root_uptake()
    B_k  (Eq. binding)        -> binding_factors()
    dC/dt (Eqs. root..fruit)  -> rhs()
    Q_Phl, C_Phl (Eqs. Qphl, Cphl) -> inside rhs()

Unit system (internally consistent; swap freely, just stay consistent):
    time            day
    aqueous conc    ug/L        (C_w^o, C_w,k, C_Phl)
    tissue conc     ug/kg       (state variable C_k)
    mass            kg
    volumetric flow L/day       (Q_TP, Q_Phl)
    binding factor  L/kg        (B_k);  C_k = B_k * C_w,k
    membrane cond.  L/(day*kg)  (kappa_d, see note)

Notes
-----
* B_k uses the Briggs-consistent form  B_k = theta_k + sum_i f_i K_i  with mass
  fractions f_i [kg/kg] and partition coeffs K_i [L/kg]; there is NO density
  prefactor (that term in an early draft was dimensionally inconsistent).
* kappa_d := a_R * P_d^eff lumps specific root membrane area and effective ionic
  permeability into one mass-specific conductance, consistent with the
  identifiability result that BAF data constrain only the lumped influx
  conductance g_in, not a_R and P_d^eff separately.
"""

from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

# ----------------------------------------------------------------------------
# physical constants
# ----------------------------------------------------------------------------
F_FARADAY = 96485.33212   # C / mol
R_GAS = 8.314462618       # J / (mol K)

ROOT, STEM, LEAF, FRUIT = 0, 1, 2, 3   # compartment indices

# Trapp cell-model pH anchors used by the weak-electrolyte phloem ion trap.
LEAF_CYTOSOL_PH = 7.2     # leaf cytosol (the phloem loading SOURCE)
PHLOEM_PH = 8.0           # sieve-tube sap: ALKALINE -> traps a weak acid
# Neutral / ionic membrane permeability ratio P_n/P_d (Trapp 2000).  Re-exported by
# literature_params.PN_OVER_PD; defined here because the ODE itself needs it.
P_N_OVER_P_D = 10.0 ** 3.5


def _speciation(pKa: float, pH: float, is_acid: bool = True):
    """(f_n, f_d) by Henderson-Hasselbalch.  Mirrors literature_params.speciation,
    duplicated here only to keep this module free of project imports (that module
    imports Compound from here)."""
    ex = 1.0 if is_acid else -1.0
    fd = 1.0 / (1.0 + 10.0 ** (ex * (pKa - pH)))
    return 1.0 - fd, fd


def _ion_trap(pKa: float, pH_source: float, pH_sink: float, is_acid: bool = True) -> float:
    """Equilibrium ion-trap enrichment between two compartments (see
    literature_params.ion_trap_factor for the full derivation and caveats)."""
    ex = 1.0 if is_acid else -1.0
    return ((1.0 + 10.0 ** (ex * (pH_sink - pKa)))
            / (1.0 + 10.0 ** (ex * (pH_source - pKa))))


# ----------------------------------------------------------------------------
# Parameter containers (organised by the Tier scheme of the report)
# ----------------------------------------------------------------------------
@dataclass
class Environment:
    """Tier 0 -- known / measurable."""
    T: float = 298.15      # temperature [K]
    E: float = -0.120      # plasmalemma membrane potential [V] (inside-negative)
    z: int = -1            # PFAS anion valence

    @property
    def N(self) -> float:
        """Dimensionless electrochemical driving force  N = zEF/(RT)."""
        return self.z * self.E * F_FARADAY / (R_GAS * self.T)

    def N_for(self, cmpd=None) -> float:
        """N for a specific compound: valence is a property of the CHEMICAL, not the
        environment (a weak base is a cation, z=+1). `Compound.z` wins when supplied;
        otherwise this falls back to `Environment.z` and reproduces `N` exactly."""
        z = self.z if (cmpd is None or getattr(cmpd, "z", None) is None) else cmpd.z
        return z * self.E * F_FARADAY / (R_GAS * self.T)


@dataclass
class Compound:
    """Per-PFAS (chain-length-specific) properties.

    Tier 2 (need inhibitor/kinetic data to separate):  kappa_d, Vmax_*, Km_*
    Tier 3 (independent measurement / QSPR):            K_prot, K_PL, K_cw, L_Ph
    """
    name: str
    # binding partition coefficients [L/kg]
    K_prot: float
    K_PL: float
    K_cw: float
    # root membrane uptake
    kappa_d: float          # lumped ionic conductance a_R*P_d^eff [L/(day kg)]
    Vmax_in: float          # carrier influx capacity [ug/(day kg)]
    Km_in: float            # influx half-saturation [ug/L]
    Vmax_out: float         # carrier efflux capacity [ug/(day kg)]
    Km_out: float           # efflux half-saturation [ug/L]
    # phloem loading partition (carrier/channel, NOT pH ion-trap) [-]
    L_Ph: float
    # root -> xylem loading factor (transpiration-stream concentration factor,
    # TSCF analog) [-]: fraction of the root free aqueous concentration actually
    # loaded into the ascending xylem.  <1 for anions because the endodermal
    # Casparian barrier and tissue binding sequester PFAS in the root and limit
    # root-to-shoot translocation.  f_xy = 1 recovers the unrestricted DPU base
    # (assumption A2); the PFAS limit is f_xy << 1, larger for short chains.
    f_xy: float = 0.1
    # root SURFACE/plaque sorption [L/kg]: Fe/Mn-plaque + rhizoplane + apoplast
    # adsorption that is in equilibrium with the EXTERNAL pore water (Cwo) and is a
    # dead-end pool -- it adds to the measured root burden but does NOT translocate.
    # Field/soil-dependent: ~0 for high-carbon Andosol (Yamazaki2023, root sub-equilibrium)
    # but >0 for high-loading fields (Li2025/Tianjin). Add to root BAF only; ODE unchanged.
    K_surf: float = 0.0
    # lipid-facilitated ("bound") loading conductances [L/kg].  The free-anion
    # xylem/phloem loading (f_xy*Cw, L_Ph*Cw) scales with the FREE aqueous conc
    # Cw=C/B, which collapses ~1/B for high-binding long chains and starves the
    # shoot of them -- yet long chains are observed in straw AND grain. These add
    # a B-INDEPENDENT term that loads the membrane/lipid-associated (bound) pool:
    #   xylem  Cw_xyl = f_xy*Cw_root + g_xy*C_root
    #   phloem C_Phl  = L_Ph*Cw_leaf + g_ph*C_leaf
    # so long chains ride the lipid phase across the endodermis/into the phloem.
    # EXPLORATORY / opt-in: default 0 recovers the free-only model exactly. The
    # value is K_PL-gated (off for short chains); see docs/fxy_longchain_lipid_
    # exploration.md and lipid_loading_conductances(). Still in-sample (Yamazaki).
    g_xy: float = 0.0
    g_ph: float = 0.0
    # speciation (PFAS: fully dissociated)
    fd: float = 1.0
    fn: float = 0.0
    # --- neutral / weak-electrolyte extension (DPU base) -----------------------
    # Every field below defaults to the PFAS limit, where its term vanishes
    # identically, so the anion path is numerically unchanged.  The speciation
    # SWITCH is the (fn, fd) pair, NOT the valence: a weak acid carries a neutral
    # molecule (potential-independent) AND an anion (potential-dependent) at the
    # same time, so one global z cannot express it.  See root_uptake().
    #   neutral      fn=1, fd=0  -> only the P_n Fickian term is alive
    #   weak acid    fn+fd=1     -> both terms alive
    #   PFAS         fn=0, fd=1  -> only the GHK term is alive (current behaviour)
    pKa: float | None = None    # acid dissociation constant; None -> use fn/fd as given
    is_acid: bool = True        # False for a weak base (conjugate acid is the CATION)
    z: int | None = None        # ion valence; None -> Environment.z (acids -1, bases +1)
    P_n: float = 0.0            # a_R*P_n, neutral passive conductance [L/(day kg)].
                                # Trapp: the ion is ~10^3.5 less permeable, P_d ~ P_n*10^-3.5.
    K_lip: float = 0.0          # Briggs lipid partition a*K_ow^b [L/kg lipid]; pairs with
                                # Compartment.f_lip (TOTAL lipid), NOT with f_PL.
    K_AW: float = 0.0           # air-water partition [-]; 0 keeps gas exchange off (A3)
    mol_weight: float = 200.0   # molar mass [g/mol]; only the air-exchange QSPRs use it
    logKow: float | None = None  # only the cuticle permeability QSPR uses it


@dataclass
class Compartment:
    """Per-tissue composition & properties (Tier 0/1)."""
    name: str
    theta: float            # aqueous (water) content [L/kg]
    f_prot: float           # protein mass fraction [kg/kg]
    f_PL: float             # phospholipid mass fraction [kg/kg]
    f_cw: float             # cell-wall mass fraction [kg/kg]
    S: float = 0.0          # specific surface area [m^2/kg] (only leaf/fruit ratio used)
    gamma: float = 0.0      # first-order metabolism [1/day] (PFAS ~ 0)
    # TOTAL lipid mass fraction [kg/kg dw] -- the Briggs sorptive phase, paired with
    # Compound.K_lip.  DISTINCT from f_PL (phospholipid): in the leaf the bulk of the
    # membrane lipid is thylakoid galactolipid (MGDG/DGDG), not phospholipid, so the
    # two must not be conflated.  See params/rice_tissue_params.csv `total_lipid`
    # (root 0.020, stem 0.015, leaf 0.055, brown grain 0.030) and its leaf CAUTION note.
    # 0 keeps the term off, so the PFAS binding factor is unchanged.
    f_lip: float = 0.0
    pH: float | None = None  # compartment pH; drives per-compartment speciation of a
                             # weak electrolyte. None -> no speciation calc (PFAS).
    # NOTE: there is deliberately NO density field.  Density does not appear anywhere
    # in this model's transport, not even in the air-exchange block, where the rho in
    # the report's Q_VOL prefactor cancels the rho inside K_PA exactly (see
    # RiceUptakeModel.air_exchange).  Densities for REPORTING live in
    # model_api.DEFAULT_TISSUE_DENSITY.


@dataclass
class PlantInputs:
    """Time-dependent external drivers (from soil model + growth model).

    Provide arrays on a common time grid `t` [day]:
        Cwo : soil pore-water free concentration  [ug/L]   (from HYDRUS/Phydrus)
        Qtp : transpiration stream                 [L/day]  (from HYDRUS/Phydrus)
        M   : tissue fresh mass, shape (len(t), 4) [kg]     (from growth model)
    Cubic/linear interpolants are built for use inside the ODE RHS.
    """
    t: np.ndarray
    Cwo: np.ndarray
    Qtp: np.ndarray
    M: np.ndarray            # shape (len(t), 4)
    leaf_loss: np.ndarray | None = None   # leaf senescence loss RATE [1/day] (opt; 0 if None)

    def __post_init__(self):
        self.M = np.asarray(self.M, dtype=float)
        assert self.M.shape[1] == 4, "M must have 4 columns (root,stem,leaf,fruit)"
        kw = dict(kind="linear", bounds_error=False, fill_value="extrapolate")
        self._Cwo = interp1d(self.t, self.Cwo, **kw)
        self._Qtp = interp1d(self.t, self.Qtp, **kw)
        self._M = [interp1d(self.t, self.M[:, k], **kw) for k in range(4)]
        # dM/dt by finite difference, then interpolate
        dM = np.gradient(self.M, self.t, axis=0)
        self._dM = [interp1d(self.t, dM[:, k], **kw) for k in range(4)]
        ll = np.zeros(len(self.t)) if self.leaf_loss is None else np.asarray(self.leaf_loss, float)
        self._leaf_loss = interp1d(self.t, ll, **kw)
        # grain FORMATION gate gamma(t): 0 while the grain mass is held at its floor (the
        # panicle has not set, ~pre-flowering), ramping to 1 once it forms (>2% of its
        # season max). The grain's PFAS influx is gated by this, so no solute loads a
        # not-yet-formed organ (DPU-consistent: import is tied to organ existence). It is
        # 1 throughout for a constant-mass driver (no floor period -> grain always present).
        gm = self.M[:, FRUIT]
        glo, ghi = float(gm.min()), float(gm.max())
        # ramp 0->1 as the grain mass LEAVES its floor (glo -> 1.5*glo): gam=1 for the whole
        # of grain filling (so loading is normal), 0 only while it sits at the floor (not set).
        gate = (np.clip((gm - glo) / (0.5 * glo + 1e-30), 0.0, 1.0)
                if (ghi > 0.0 and glo < 0.05 * ghi) else np.ones(len(self.t)))
        self._grain_gate = interp1d(self.t, gate, **kw)

    def Cwo_(self, t):  return float(self._Cwo(t))
    def Qtp_(self, t):  return float(self._Qtp(t))
    def M_(self, t):    return np.array([float(f(t)) for f in self._M])
    def dM_(self, t):   return np.array([float(f(t)) for f in self._dM])
    def leaf_loss_(self, t):  return max(float(self._leaf_loss(t)), 0.0)
    def grain_gate_(self, t):  return min(1.0, max(0.0, float(self._grain_gate(t))))


# ----------------------------------------------------------------------------
# Model functions
# ----------------------------------------------------------------------------
def binding_factors(comps: list[Compartment], cmpd: Compound) -> np.ndarray:
    """B_k = theta_fw + (1-theta_fw)*(f_prot K_prot + f_PL K_PL + f_cw K_cw)  [L/kg fw].

    FRESH-WEIGHT basis (NOTES path A): theta is the fresh-weight water fraction and
    f_* are DRY-weight mass fractions, so the dw binding fractions are rescaled by
    (1 - theta_fw) onto a per-kg-fresh basis -- consistent with the model's fresh-mass
    M and fluxes. Compare to dw-reported data via C_dw = C_fw / (1 - theta_fw).
    4-POOL: f_cw is the WHOLE cell wall (polysaccharide + lignin lumped).

    NEUTRAL/DPU: the extra `f_lip*K_lip` term makes this the general form of the
    Briggs plant-water partition K_PW = W + a*K_ow^b*L -- set K_prot=K_cw=0 and it
    reduces to Briggs exactly.  Both `f_lip` and `K_lip` default to 0, so the PFAS
    binding factor is bit-identical to before.
    """
    return np.array([
        c.theta + (1.0 - c.theta) * (
            c.f_prot * cmpd.K_prot + c.f_PL * cmpd.K_PL + c.f_cw * cmpd.K_cw
            + c.f_lip * cmpd.K_lip)
        for c in comps
    ])


def _ghk_factor(N: float) -> float:
    """N / (exp(N) - 1), with the removable singularity at N=0 handled."""
    if abs(N) < 1e-8:
        return 1.0
    return N / np.expm1(N)


def root_uptake(Cwo: float, Cw_root: float, cmpd: Compound, env: Environment) -> float:
    """Mass-specific root membrane uptake j_R [ug/(day kg)]  (Eq. JR_pfas).

    Three PARALLEL pathways, each weighted by the species fraction it carries:

        j_R = P_n*f_n*(Cwo - Cw)              neutral passive (potential-INdependent)
            + kappa_d*g*f_d*(Cwo - e^N*Cw)    ionic electrodiffusion (GHK)
            + carrier (Michaelis-Menten)

    The GHK factor applies to the ION term ONLY -- that is what lets one code path
    cover the whole speciation spectrum (see Compound):
      * neutral    f_d=0 kills the GHK term, so the valence is irrelevant;
      * weak acid  both terms are alive at the same time;
      * PFAS       f_n=0 kills the neutral term -> identical to the previous code.
    """
    N = env.N_for(cmpd)
    eN = np.exp(N)
    g = _ghk_factor(N)
    # neutral passive permeation (Fickian; no membrane-potential term).  P_n defaults
    # to 0, and 0.0 + x == x exactly, so the PFAS result is bit-identical.
    j_n = cmpd.P_n * (cmpd.fn * Cwo - cmpd.fn * Cw_root)
    # ionic electrodiffusion (membrane + anion channel + aquaporin, lumped in kappa_d)
    j_ed = cmpd.kappa_d * g * (cmpd.fd * Cwo - cmpd.fd * eN * Cw_root)
    # carrier-mediated (active/facilitated), net influx - efflux
    j_carr = (cmpd.Vmax_in * Cwo / (cmpd.Km_in + Cwo)
              - cmpd.Vmax_out * Cw_root / (cmpd.Km_out + Cw_root))
    return j_n + j_ed + j_carr


# ----------------------------------------------------------------------------
# Air exchange (DPU base sections 5.4-5.8) -- OFF unless Compound.K_AW > 0
# ----------------------------------------------------------------------------
# The permeability correlations below are stated in SI (metres, seconds, g/mol)
# while this model runs in day / L / kg / ug.  Every air flux therefore ends up as
#     [m2/kg] * [m/s] * [ug/L]  =  m3/(L*kg*s) * ug
# and m3 = 1000 L, so ONE constant converts all three of them:
M_S_TO_L_DAY = 1000.0 * 86400.0     # (m3->L) * (s->day) = 8.64e7
RHO_WATER_SI = 1000.0               # kg/m3


@dataclass
class AirInputs:
    """Atmospheric drivers for the plant-air exchange block.

    Attach to :class:`RiceUptakeModel` as ``air=``; the whole block stays inert
    unless the compound also has ``K_AW > 0`` (assumption A3 for PFAS).

    Two fields are deliberately NOT named as in the DPU report, because both of
    those symbols are already taken in this model:
      ``RH``      is the report's phi (relative humidity); ``RiceUptakeModel.phi``
                  is the phloem recirculation fraction.
      ``z_path``  is the report's z (aqueous diffusion path); ``Environment.z``
                  is the ion valence.
    """
    C_A: float = 0.0        # atmospheric concentration [ug/L air]
    f_p: float = 0.0        # particle-bound fraction of C_A [-]
    v_dep: float = 1e-3     # particle deposition velocity [m/s] (~0.001)
    RH: float = 0.7         # relative humidity [-]  (report: phi)
    z_path: float = 1e-3    # aqueous boundary-layer path [m]  (report: z)
    D_O2: float = 2.1e-9    # O2 diffusivity in water [m2/s] at 25 C


def _c_h2o_sat(T: float) -> float:
    """Saturation water-vapour concentration [kg/m3] at temperature T [K]."""
    T_C = T - 273.15
    p = 610.7 * 10.0 ** (7.5 * T_C / (237.0 + T_C))          # Pa
    return p / (461.9 * T)                                    # R/M_water = 461.9


def permeabilities(cmpd: Compound, env: Environment, air: AirInputs,
                   Qxyl: np.ndarray, A: np.ndarray) -> np.ndarray:
    """Total plant-air permeability P_P per compartment [m/s]  (report Eq. Pp).

    Cuticle, air boundary layer and the aqueous layer act in SERIES; the resulting
    cuticular path and the stomatal path act in PARALLEL.  Per the report's
    modelling assumptions: no exchange from the root (below ground), cuticle only
    for the stem, cuticle + stomata for leaf and fruit.

    ``A`` is the ABSOLUTE compartment surface area [m2] (= S * M), not the specific
    area: the stomatal term is Q_XYL/A, so mixing a per-kg S with an absolute
    Q_XYL would be a silent basis error.  (In the flux the A cancels again, but
    keeping the bases consistent here is what makes that safe to rely on.)
    """
    m = cmpd.mol_weight
    logKow = 0.0 if cmpd.logKow is None else cmpd.logKow
    P_C = 10.0 ** (0.704 * logKow - 11.2)                     # cuticle
    P_air = np.sqrt(300.0) * cmpd.K_AW / (200.0 * m ** 0.5)   # air boundary layer
    P_aqua = (air.D_O2 / air.z_path) * (m / 32.0) ** -0.5     # aqueous layer
    inv = sum(1.0 / p for p in (P_C, P_air, P_aqua) if p > 0.0)
    P_C_tot = 1.0 / inv if inv > 0.0 else 0.0

    # stomatal path, tied to the transpiration stream (report Eq. Ps)
    P_S = np.zeros(4)
    denom = (1.0 - min(air.RH, 0.999)) * _c_h2o_sat(env.T)
    if denom > 0.0 and cmpd.K_AW > 0.0:
        for k in (LEAF, FRUIT):
            if A[k] <= 0.0:
                continue
            P_S[k] = (RHO_WATER_SI / denom) * (m / 18.0) ** -0.5 * cmpd.K_AW * (
                Qxyl[k] / M_S_TO_L_DAY) / A[k]           # Q_XYL [L/day] -> [m3/s]

    P_P = np.zeros(4)
    for k in (STEM, LEAF, FRUIT):                             # ROOT stays 0
        P_P[k] = P_C_tot + P_S[k]
    return P_P


@dataclass
class RiceUptakeModel:
    env: Environment
    cmpd: Compound
    comps: list[Compartment]            # [root, stem, leaf, fruit]
    inputs: PlantInputs
    phi: float = 0.1                    # phloem recirculation fraction to roots [-]
    T_C_Ph: float = 10.0                # phloem flux per unit grain dry mass [L/kg]
    phloem_pH: float = PHLOEM_PH        # sieve-tube sap pH (weak-electrolyte ion trap)
    air: "AirInputs | None" = None      # atmospheric exchange; inert unless K_AW > 0

    def air_exchange(self, C: np.ndarray, B: np.ndarray, M: np.ndarray,
                     Qxyl: np.ndarray) -> np.ndarray:
        """Net atmospheric exchange per compartment [ug/(kg day)]  (Eqs. Qgas/Qvol/Qdep).

            +Q_GAS  (1-f_p) * C_A * (A/M) * P_P            gaseous uptake
            +Q_DEP  v_dep * f_p * C_A * (A/M)              particle deposition
            -Q_VOL  (A/M) * P_P * C * rho / K_PA           volatilisation

        DENSITY DOES NOT ENTER.  The report writes volatilisation as
        (A*rho/M) * P_P * C / K_PA with K_PA = K_PW*rho/K_AW, and the rho in the
        prefactor cancels the rho inside K_PA exactly:

            (A*rho/M) * P_P * C / (K_PW*rho/K_AW)  ==  (A/M) * P_P * C * K_AW / K_PW

        so the air-side concentration in equilibrium with the tissue is simply
        C*K_AW/B_k.  This was worth checking rather than assuming: the design notes
        for this phase predicted that tissue density would finally become a real
        transport quantity here, and it does not.  The model's standing convention
        (no density prefactor anywhere in transport) therefore holds unchanged.

        Returns zeros unless an ``AirInputs`` is attached AND ``K_AW > 0``, which is
        what keeps the PFAS path (assumption A3: non-volatile) untouched.
        """
        out = np.zeros(4)
        if self.air is None or self.cmpd.K_AW <= 0.0:
            return out
        A = np.array([c.S for c in self.comps]) * M               # absolute area [m2]
        P_P = permeabilities(self.cmpd, self.env, self.air, Qxyl, A)
        a_over_m = A / np.maximum(M, 1e-12)                        # [m2/kg]
        air = self.air
        for k in range(4):
            if a_over_m[k] <= 0.0 or (P_P[k] <= 0.0 and air.f_p <= 0.0):
                continue
            c_air_eq = C[k] * self.cmpd.K_AW / B[k]                # [ug/L air]
            gas = (1.0 - air.f_p) * air.C_A * a_over_m[k] * P_P[k]
            dep = air.v_dep * air.f_p * air.C_A * a_over_m[k]
            vol = a_over_m[k] * P_P[k] * c_air_eq
            out[k] = (gas + dep - vol) * M_S_TO_L_DAY
        return out

    def phloem_loading_factor(self) -> float:
        """Effective leaf->phloem loading partition [-]  (``C_Phl = L * Cw_leaf``).

        Two PARALLEL routes, mirroring ``root_uptake``:

          carrier/channel   L_Ph          -- the fitted PFAS route (assumption A5)
          neutral + trap    Lambda        -- the neutral species crosses the sieve-tube
                                             membrane and re-dissociates in the alkaline
                                             sap, so the TOTAL phloem concentration is
                                             Lambda times the leaf's free concentration

        The trap route is weighted by whether the neutral species can actually deliver
        it, w = Pi/(1+Pi) with Pi = (P_n/P_d)*f_n/f_d -- the permeability-weighted
        neutral fraction.  This is the correction that phase 1.5 exists for: multiplying
        L_Ph by Lambda instead would hand a PERMANENT anion a spurious ~6.3x phloem
        enrichment, because Lambda tends to 10**(dpH) rather than to 1 as pKa falls.
        The trap switches off kinetically (f_n -> 0), not thermodynamically.

        REQUIRES an explicit ``Compound.pKa`` AND a leaf ``Compartment.pH``.  Neither is
        set on the PFAS path, so that path returns ``L_Ph`` unchanged BY CONSTRUCTION --
        a guarantee, not a numerical coincidence (see the note in the module docstring
        about supplying a PFAS pKa deliberately).

        The gate w is phenomenological: the model has no phloem TRANSIT compartment, so
        it represents loading at the leaf, not retention during transport.  The classic
        extra phloem mobility of weak acids comes partly from that retention and is
        therefore only partly captured here.

        Recomputed per RHS call rather than cached: ``calibration.py`` fits by
        ``setattr``-ing L_Ph onto the compound of an already-built model, and a cached
        value would silently ignore that.  The PFAS branch is a single attribute test.
        """
        c = self.cmpd
        pH_leaf = self.comps[LEAF].pH
        if c.pKa is None or pH_leaf is None:
            return c.L_Ph
        fn, fd = _speciation(c.pKa, pH_leaf, c.is_acid)
        lam = _ion_trap(c.pKa, pH_leaf, self.phloem_pH, c.is_acid)
        if fd <= 0.0:                       # strictly neutral: nothing to trap
            return c.L_Ph
        pi = P_N_OVER_P_D * fn / fd
        w = pi / (1.0 + pi)
        return (1.0 - w) * c.L_Ph + w * lam

    def rhs(self, t: float, C: np.ndarray) -> np.ndarray:
        """RHS of dC/dt for the 4 compartments (Eqs. root, stem, leaf, fruit)."""
        Cwo = self.inputs.Cwo_(t)
        Qtp = self.inputs.Qtp_(t)
        M = self.inputs.M_(t)
        dM = self.inputs.dM_(t)
        M = np.maximum(M, 1e-12)                     # guard against division by zero
        mu = dM / M                                  # growth-dilution rates [1/day]

        B = binding_factors(self.comps, self.cmpd)   # [L/kg]
        Cw = C / B                                   # free aqueous conc [ug/L]

        # leaf/fruit xylem split by surface-area fraction
        A3 = self.comps[LEAF].S * M[LEAF]
        A4 = self.comps[FRUIT].S * M[FRUIT]
        split = A3 / (A3 + A4) if (A3 + A4) > 0 else 0.5
        f3, f4 = split, 1.0 - split

        # phloem flow and sap concentration (carrier loading at leaf; NOT pH trap)
        Q_Phl = dM[FRUIT] * self.T_C_Ph + self.phi * Qtp     # [L/day]
        Q_Phl = max(Q_Phl, 0.0)
        # L_Ph_eff == cmpd.L_Ph unless a weak-electrolyte pH trap is configured
        C_Phl = self.phloem_loading_factor() * Cw[LEAF] + self.cmpd.g_ph * C[LEAF]   # free + lipid-bound [ug/L]

        g = [c.gamma for c in self.comps]
        dC = np.zeros(4)

        # root -> xylem loading: a fraction f_xy of the root FREE conc (TSCF;
        # endodermal barrier + binding) PLUS a B-independent lipid-bound term
        # g_xy*C_root (default 0). The SAME Cw_xyl is what the root loses and the
        # stem gains, so the root->stem xylem transfer remains mass-conserving.
        Cw_xyl = self.cmpd.f_xy * Cw[ROOT] + self.cmpd.g_xy * C[ROOT]

        # root
        jR = root_uptake(Cwo, Cw[ROOT], self.cmpd, self.env)
        dC[ROOT] = (jR
                    - (Qtp / M[ROOT]) * Cw_xyl
                    + self.phi * (Q_Phl / M[ROOT]) * C_Phl
                    - g[ROOT] * C[ROOT] - mu[ROOT] * C[ROOT])
        # stem (xylem in from root at the loaded conc, xylem out to leaf+fruit)
        dC[STEM] = ((Qtp / M[STEM]) * (Cw_xyl - Cw[STEM])
                    - g[STEM] * C[STEM] - mu[STEM] * C[STEM])
        # leaf (xylem terminal; phloem source). The leaf supplies the WHOLE
        # phloem export: the grain sink (Q_Phl*C_Phl) plus the fraction phi*Q_Phl
        # that recirculates to the root -> total (1+phi)*Q_Phl*C_Phl. This closes
        # the phloem mass balance (leaf loss = grain gain + root recirculation).
        # leaf senescence loss: the dead/shed leaf carries its PFAS away at the leaf
        # death rate (D/M_leaf), cancelling the spurious -mu*C concentration when the
        # leaf shrinks. 0 unless a senescing biomass driver (ORYZA) supplies the rate.
        # grain FORMATION gate: before the panicle sets (gam~0) the grain takes NO xylem/
        # phloem; its share is rerouted to the leaf (xylem f4 -> leaf; phloem export drops
        # to (gam+phi)) so the balance still closes. gam->1 at formation (DPU-consistent).
        gam = self.inputs.grain_gate_(t)
        dC[LEAF] = ((f3 + (1.0 - gam) * f4) * (Qtp / M[LEAF]) * Cw[STEM]
                    - (gam + self.phi) * (Q_Phl / M[LEAF]) * C_Phl
                    - g[LEAF] * C[LEAF] - mu[LEAF] * C[LEAF]
                    - self.inputs.leaf_loss_(t) * C[LEAF])
        # fruit/grain (small xylem in; phloem-dominated; terminal sink) -- gated by gam
        dC[FRUIT] = (gam * f4 * (Qtp / M[FRUIT]) * Cw[STEM]
                     + gam * (Q_Phl / M[FRUIT]) * C_Phl
                     - g[FRUIT] * C[FRUIT] - mu[FRUIT] * C[FRUIT])
        # plant-air exchange (gaseous uptake + deposition - volatilisation).  Returns
        # exactly zeros for a non-volatile compound, so the PFAS path is untouched.
        if self.air is not None and self.cmpd.K_AW > 0.0:
            Qxyl = np.array([Qtp, Qtp, f3 * Qtp, f4 * Qtp])
            dC = dC + self.air_exchange(C, B, M, Qxyl)
        return dC

    def solve(self, t_eval: np.ndarray, C0: np.ndarray | None = None):
        """Integrate the stiff system with BDF. Returns scipy solution object."""
        if C0 is None:
            C0 = np.zeros(4)
        t_span = (float(t_eval[0]), float(t_eval[-1]))
        sol = solve_ivp(self.rhs, t_span, C0, t_eval=t_eval,
                        method="BDF", rtol=1e-6, atol=1e-9, dense_output=True)
        return sol

    def baf(self, C: np.ndarray, t: float) -> np.ndarray:
        """Bioaccumulation factor BAF_k = C_k / C_w^o  [L/kg]."""
        return C / self.inputs.Cwo_(t)

    def root_baf_total(self, C: np.ndarray, t: float) -> float:
        """Total root BAF = internal partitioning (ODE state) + surface sorption.

        BAF_root,total = C[ROOT]/Cwo + K_surf. The surface pool is a parallel
        equilibrium with the external pore water (Cwo), so it leaves the transport
        ODE untouched and only inflates the *measured* root concentration. With
        K_surf=0 this reduces to the internal BAF. Use K_surf>0 only where field
        data show root BAF exceeding the internal ceiling B_root (high-loading soils).
        """
        return float(C[ROOT] / self.inputs.Cwo_(t) + self.cmpd.K_surf)


# ----------------------------------------------------------------------------
# Synthetic demo (replace inputs with Phydrus/HYDRUS output + a growth model)
# ----------------------------------------------------------------------------
def _logistic(t, M0, Mmax, k, t0):
    return Mmax / (1.0 + (Mmax / M0 - 1.0) * np.exp(-k * (t - t0)))


def _demo():
    season = 120.0                       # days
    t = np.linspace(0.0, season, 481)    # 0.25-day grid

    # --- external drivers (PLACEHOLDERS; supply from HYDRUS/Phydrus + growth model)
    Cwo = np.full_like(t, 1.0)           # soil pore-water free conc [ug/L], constant
    # transpiration: rises with canopy, peaks mid/late season [L/day]
    Qtp = 0.05 + 0.35 * np.exp(-((t - 75.0) ** 2) / (2 * 25.0 ** 2))
    # growth: root/stem/leaf early; grain (fruit) fills from ~flowering (day 65)
    M = np.column_stack([
        _logistic(t, 1e-3, 0.030, 0.10, 20.0),    # root  [kg]
        _logistic(t, 1e-3, 0.040, 0.10, 25.0),    # stem
        _logistic(t, 1e-3, 0.050, 0.12, 30.0),    # leaf
        _logistic(t, 1e-5, 0.025, 0.18, 80.0),    # fruit/grain (late)
    ])
    inputs = PlantInputs(t=t, Cwo=Cwo, Qtp=Qtp, M=M)

    env = Environment()                  # N ~ +4.67 for z=-1, E=-120 mV
    # PFOA-like placeholder compound (values are illustrative, NOT calibrated).
    # The three transport parameters below encode the PFAS-specific physics that
    # produces the empirical root > straw > grain ordering:
    #   f_xy = 0.02  -> low transpiration-stream loading: the anion is strongly
    #                   retained in the root and translocates poorly to the shoot
    #                   (endodermal Casparian barrier + binding); larger for short
    #                   chains, smaller for long chains.
    #   L_Ph = 0.005 -> PFAS is poorly phloem-mobile (no weak-acid ion-trap for a
    #                   fully dissociated anion), so the phloem-fed grain is the
    #                   most poorly supplied compartment.
    cmpd = Compound(
        name="PFOA",
        K_prot=50.0, K_PL=100.0, K_cw=7.0,         # [L/kg] PROVISIONAL: K_cw=whole-cw eff(Mel+poly); K_PL pending QC1 down-weight
        kappa_d=0.5,                                # [L/(day kg)]
        Vmax_in=20.0, Km_in=5.0,                    # carrier influx (must overcome anion exclusion)
        Vmax_out=8.0, Km_out=5.0,                   # carrier efflux
        L_Ph=0.005,                                 # low phloem loading (anion: no pH ion-trap)
        f_xy=0.02,                                  # root->xylem loading (TSCF): limited translocation
    )
    # rice_tissue dw composition (recommended anchors); theta = fresh-weight water
    # fraction; f_* = DRY-weight fractions; f_cw = WHOLE cell wall (poly + lignin).
    comps = [
        Compartment("root",  theta=0.90, f_prot=0.07, f_PL=0.015, f_cw=0.50),
        Compartment("stem",  theta=0.83, f_prot=0.05, f_PL=0.005, f_cw=0.72),
        Compartment("leaf",  theta=0.78, f_prot=0.10, f_PL=0.010, f_cw=0.56, S=20.0),
        Compartment("grain", theta=0.14, f_prot=0.09, f_PL=0.003, f_cw=0.035, S=2.0),
    ]

    model = RiceUptakeModel(env=env, cmpd=cmpd, comps=comps, inputs=inputs)
    sol = model.solve(t)

    B = binding_factors(comps, cmpd)
    Cend = sol.y[:, -1]
    print(f"electrochemical number N = {env.N:.3f}  (e^N = {np.exp(env.N):.1f})")
    print(f"binding factors B_k [L/kg]: " +
          ", ".join(f"{c.name}={b:.2f}" for c, b in zip(comps, B)))
    print("\nfinal tissue concentrations & BAFs:")
    baf = model.baf(Cend, t[-1])
    for c, ck, bk in zip(comps, Cend, baf):
        print(f"  {c.name:5s}  C = {ck:8.3f} ug/kg   BAF = {bk:7.3f} L/kg")
    # empirical target ordering root > straw > grain, where "straw" is the bulk
    # shoot = mass-weighted mean of stem + leaf (standard rice-agronomy tissue).
    Mf = inputs.M_(t[-1])
    straw = (Cend[STEM] * Mf[STEM] + Cend[LEAF] * Mf[LEAF]) / (Mf[STEM] + Mf[LEAF])
    straw_baf = straw / inputs.Cwo_(t[-1])
    ordered = baf[ROOT] > straw_baf > baf[FRUIT]
    print(f"\nordering check (target: root > straw > grain):")
    print(f"  root  BAF = {baf[ROOT]:6.2f}")
    print(f"  straw BAF = {straw_baf:6.2f}   (mass-weighted stem+leaf)")
    print(f"  grain BAF = {baf[FRUIT]:6.2f}")
    print(f"  root > straw > grain : {'OK' if ordered else 'VIOLATED'}"
          f"  (root/straw = {baf[ROOT]/straw_baf:.2f}, straw/grain = {straw_baf/baf[FRUIT]:.2f})")

    # optional figure
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7, 4.2))
        for k, c in enumerate(comps):
            ax.plot(t, sol.y[k], label=c.name, lw=2)
        ax.set_xlabel("time [day]"); ax.set_ylabel("tissue conc [ug/kg]")
        ax.set_title("PFAS rice 4-compartment uptake (synthetic demo)")
        ax.legend(); fig.tight_layout()
        fig.savefig("pfas_rice_demo.png", dpi=130)
        print("\nsaved figure: pfas_rice_demo.png")
    except Exception as e:
        print(f"\n(plot skipped: {e})")

    return model, sol


if __name__ == "__main__":
    _demo()
