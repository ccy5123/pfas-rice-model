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

# ----------------------------------------------------------------------------
# Weak-electrolyte speciation (ported from PR #54)
# ----------------------------------------------------------------------------
# This block is what lets ONE code path cover the whole speciation spectrum:
#
#     neutral (fn=1, fd=0)  <->  weak acid/base (0<fn,fd<1)  <->  PFAS (fn=0, fd=1)
#
# WHY IT IS NOT `z = 0`. `src/neutral_dpu.py` reaches a NEUTRAL compound by setting
# z=0, which makes N=0, the GHK factor 1 and the membrane term degenerate exactly to
# passive Fickian transport. That is correct for a strictly neutral molecule and is
# how every published neutral number in this repo was produced. It cannot express a
# WEAK ACID, which is a neutral molecule AND an anion at the same time: one global
# valence has to be either 0 or -1, and the compound is genuinely both. The switch
# therefore has to be the (fn, fd) PAIR, with the potential-dependent GHK factor on
# the ion term only -- see `root_uptake`.
#
# Every field and constant here defaults to the PFAS limit, where its term vanishes
# identically (not approximately), so the anion path is bit-identical to before.
LEAF_CYTOSOL_PH = 7.2     # leaf cytosol -- the phloem loading SOURCE
PHLOEM_PH = 8.0           # sieve-tube sap: ALKALINE, so it traps a weak acid
# Neutral/ionic membrane permeability ratio P_n/P_d (Trapp 2000). Also re-exported by
# literature_params.PN_OVER_PD; it lives here because the ODE itself needs it.
P_N_OVER_P_D = 10.0 ** 3.5


def _speciation(pKa: float, pH: float, is_acid: bool = True):
    """(f_n, f_d) by Henderson-Hasselbalch. Mirrors `literature_params.speciation`;
    duplicated here to keep this module free of project imports (that module imports
    `Compound` from this one)."""
    ex = 1.0 if is_acid else -1.0
    fd = 1.0 / (1.0 + 10.0 ** (ex * (pKa - pH)))
    return 1.0 - fd, fd


def _ion_trap(pKa: float, pH_source: float, pH_sink: float, is_acid: bool = True) -> float:
    """Equilibrium ion-trap enrichment between two compartments (see
    `literature_params.ion_trap_factor` for the derivation and caveats)."""
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
        """N for a SPECIFIC compound: valence is a property of the chemical, not of
        the environment (a weak base is a cation, z=+1, and would feel the
        inside-negative membrane as an ATTRACTION where an anion feels exclusion).
        `Compound.z` wins when supplied; otherwise this falls back to
        `Environment.z` and reproduces `N` exactly, so the PFAS path is unchanged."""
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
    # --- weak-electrolyte extension (ported from PR #54) -----------------------
    # All four default to the PFAS limit, where their terms vanish identically.
    # `fn`/`fd` above remain the switch; these add what a weak acid needs on top.
    pKa: float | None = None    # acid dissociation constant. None -> use fn/fd as given,
                                # and the phloem ion trap stays off BY CONSTRUCTION.
    is_acid: bool = True        # False for a weak base (its conjugate acid is a CATION)
    z: int | None = None        # ion valence; None -> Environment.z (acids -1, bases +1)
    P_n: float = 0.0            # a_R*P_n, NEUTRAL passive conductance [L/(day kg)].
                                # Trapp: the ion is ~10^3.5 less permeable (P_N_OVER_P_D),
                                # so P_d ~ P_n * 10^-3.5. 0 keeps the term off.
    # APOPLASTIC BYPASS conductance [L/(day kg)]: entry that never crosses a
    # membrane -- solute carried in the water stream through the cell walls and
    # past the endodermis where the Casparian band is absent or broken (root tips,
    # lateral-root emergence points). Rice is the plant this is best documented in;
    # "bypass flow" is the standard explanation for its sodium uptake.
    #
    # It is defined by what it does NOT feel, which is exactly why it is a separate
    # term rather than a bigger kappa_d: no (fn, fd) weighting and no GHK factor,
    # because a route around the membrane cannot be gated by speciation or by the
    # membrane potential. That makes it the one structural lever able to raise a
    # strongly-ionised compound's uptake WITHOUT flattening the speciation ordering
    # the data support -- see docs/neutral_dpu_validation.md section 4l.
    #
    # EXPLORATORY / opt-in: 0 (the default) removes the term identically, so every
    # PFAS and strictly-neutral number is bit-identical. NOT set for PFAS, whose
    # own anion-uptake deficit is carried by the fitted CARRIER (Vmax_in) instead;
    # whether one mechanism should serve both is open.
    g_apo: float = 0.0


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
    pH: float | None = None  # compartment pH, for the weak-electrolyte phloem ion trap
                             # (LEAF_CYTOSOL_PH is the usual leaf value). None -> no
                             # speciation calculation at all, which is the PFAS path.


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
    """
    return np.array([
        c.theta + (1.0 - c.theta) * (
            c.f_prot * cmpd.K_prot + c.f_PL * cmpd.K_PL + c.f_cw * cmpd.K_cw)
        for c in comps
    ])


def _ghk_factor(N: float) -> float:
    """N / (exp(N) - 1), with the removable singularity at N=0 handled."""
    if abs(N) < 1e-8:
        return 1.0
    return N / np.expm1(N)


def root_uptake(Cwo: float, Cw_root: float, cmpd: Compound, env: Environment) -> float:
    """Mass-specific root membrane uptake j_R [ug/(day kg)]  (Eq. JR_pfas).

    FOUR PARALLEL pathways, each weighted by the species fraction that carries it:

        j_R = P_n*f_n*(Cwo - Cw)              neutral passive (potential-INdependent)
            + kappa_d*g*f_d*(Cwo - e^N*Cw)    ionic electrodiffusion (GHK)
            + carrier (Michaelis-Menten)
            + g_apo*(Cwo - Cw)                apoplastic bypass (membrane-INdependent)

    The fourth term is the odd one out and deliberately so: it carries NO (fn, fd)
    weight and NO GHK factor, because a route around the membrane cannot be gated
    by speciation or by the membrane potential. `g_apo` defaults to 0, so it is not
    merely zero but structurally absent for PFAS and for the strict neutral path.

    The GHK factor multiplies the ION term ONLY, and that is precisely what lets one
    code path span the whole speciation spectrum:
      * neutral     f_d=0 kills the GHK term, so the membrane potential is irrelevant
      * weak acid   both terms are alive simultaneously -- the case a single global
                    valence cannot represent (see the speciation block above)
      * PFAS        f_n=0 kills the neutral term -> identical to the previous code

    `P_n` defaults to 0 and `0.0 + x == x` exactly, so the PFAS result is bit-identical.
    """
    N = env.N_for(cmpd)
    eN = np.exp(N)
    g = _ghk_factor(N)
    # neutral passive permeation (Fickian; carries no membrane-potential term)
    j_n = cmpd.P_n * (cmpd.fn * Cwo - cmpd.fn * Cw_root)
    # ionic electrodiffusion (membrane + anion channel + aquaporin, lumped in kappa_d)
    j_ed = cmpd.kappa_d * g * (cmpd.fd * Cwo - cmpd.fd * eN * Cw_root)
    # carrier-mediated (active/facilitated), net influx - efflux
    j_carr = (cmpd.Vmax_in * Cwo / (cmpd.Km_in + Cwo)
              - cmpd.Vmax_out * Cw_root / (cmpd.Km_out + Cw_root))
    # apoplastic bypass: around the membrane, so neither speciation- nor
    # potential-gated. `g_apo` defaults to 0 and `x + 0.0 == x` exactly.
    j_apo = cmpd.g_apo * (Cwo - Cw_root)
    return j_n + j_ed + j_carr + j_apo


@dataclass
class RiceUptakeModel:
    env: Environment
    cmpd: Compound
    comps: list[Compartment]            # [root, stem, leaf, fruit]
    inputs: PlantInputs
    phi: float = 0.1                    # phloem recirculation fraction to roots [-]
    T_C_Ph: float = 10.0                # phloem flux per unit grain dry mass [L/kg]
    # OPTIONAL plant-air exchange (volatilisation + gaseous uptake), see
    # `plant_air.AirExchange` and docs/dpu_model_summary_corrected.tex
    # sec:permeability. None (the default) skips the terms entirely -- they are
    # not merely zero but never evaluated, so every PFAS result is bit-identical.
    # PFAS have K_AW ~ 0 and no air pathway at all (CLAUDE.md section 2); this is
    # for the NEUTRAL path, where volatilisation is a load-bearing leaf sink.
    air: "object | None" = None
    # sieve-tube sap pH for the weak-electrolyte phloem ion trap; inert unless the
    # compound carries a pKa AND the leaf compartment carries a pH.
    phloem_pH: float = PHLOEM_PH

    def phloem_loading_factor(self) -> float:
        """Effective leaf->phloem loading partition [-]  (`C_Phl = L * Cw_leaf`).

        Two PARALLEL routes, mirroring `root_uptake`:

          carrier/channel   L_Ph     -- the fitted PFAS route (assumption A5)
          neutral + trap    Lambda   -- the neutral species crosses the sieve-tube
                                        membrane and re-dissociates in the alkaline
                                        sap, so the TOTAL phloem concentration is
                                        Lambda x the leaf's free concentration

        The trap route is weighted by whether the neutral species can actually
        DELIVER it: w = Pi/(1+Pi) with Pi = (P_n/P_d)*f_n/f_d, the permeability-
        weighted neutral fraction. This weighting is the whole point. Multiplying
        `L_Ph` by Lambda instead would hand a PERMANENT anion a spurious ~6.3x phloem
        enrichment, because Lambda tends to 10^(dpH), NOT to 1, as pKa falls: Lambda
        is an equilibrium ratio derived assuming the neutral species carries
        transport, so what vanishes for a strong acid is the FLUX that would
        establish it, not the ratio. The trap switches off kinetically (f_n -> 0),
        never thermodynamically.

        REQUIRES an explicit `Compound.pKa` AND a leaf `Compartment.pH`. Neither is
        set on the PFAS path, so that path returns `L_Ph` unchanged BY CONSTRUCTION
        -- a guarantee from the branch, not a numerical coincidence.

        The gate w is phenomenological: this model has no phloem TRANSIT
        compartment, so w represents loading at the leaf, not retention during
        transport. The classic extra phloem mobility of weak acids comes partly from
        that retention and is therefore only partly captured here.

        Recomputed per RHS call rather than cached: `calibration.py` fits by
        `setattr`-ing `L_Ph` onto the compound of an already-built model, and a
        cached value would silently ignore that. The PFAS branch is one attribute
        test, so the cost is negligible.
        """
        c = self.cmpd
        pH_leaf = self.comps[LEAF].pH
        if c.pKa is None or pH_leaf is None:
            return c.L_Ph
        fn, fd = _speciation(c.pKa, pH_leaf, c.is_acid)
        if fd <= 0.0:                       # strictly neutral: nothing to trap
            return c.L_Ph
        lam = _ion_trap(c.pKa, pH_leaf, self.phloem_pH, c.is_acid)
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
        C_Phl = (self.phloem_loading_factor() * Cw[LEAF]
                 + self.cmpd.g_ph * C[LEAF])          # free + lipid-bound [ug/L]

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
        # plant-air exchange (gaseous uptake - volatilisation), opt-in and OFF by
        # default: for `air=None` this branch is skipped, so the PFAS path is
        # untouched. The stomatal pathway needs each shoot organ's share of the
        # transpiration stream, which is the same (f3, f4) split used above, with
        # the grain's share gated by its formation.
        if self.air is not None:
            dC += self.air.flux(C=C, M=M, K_PW=B, comps=self.comps, Qtp=Qtp,
                                xyl_share=(f3, gam * f4))
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
