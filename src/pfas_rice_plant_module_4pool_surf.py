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
    Q_g, Q_r, C_Phl (Eqs. Qphl, Cphl) -> phloem_flows(), rhs()

Unit system (internally consistent; swap freely, just stay consistent):
    time            day
    aqueous conc    ug/L        (C_w^o, C_w,k, C_Phl)
    tissue conc     ug/kg       (state variable C_k)
    mass            kg
    volumetric flow L/day       (Q_TP, Q_g, Q_r)
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
    # LINEAR EQUILIBRIUM root entry: j = lam_eq * (r_eq * Cwo - Cw_root). The root's free
    # concentration relaxes toward r_eq x the pore water at rate lam_eq. It is the form the
    # data identify when membrane exchange is fast against xylem export and growth: only
    # the ratio r_eq is then constrained by a root BAF, and lam_eq is numerical (results
    # must not depend on it once it is large -- the convergence rule). It replaces the
    # GHK + Michaelis-Menten pair, whose Km/Vmax have no source and whose uptake the
    # literature finds linear at environmental levels (docs/literature_db/sourced_inputs.csv).
    # 0 (the default) removes the term identically.
    r_eq: float = 0.0
    lam_eq: float = 0.0


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
    K_cw: float | None = None  # this organ's cell-wall binding [L/kg cw]; None -> the
                             # compound's K_cw (one value for every compartment)


M_BASES = (None, "fresh")
from organ_mass import MASS_FLOOR   # the one organ-mass floor (numerical, converged)
# the integrator: "expm" (exponential stepping of the affine system, `linear_expm`) or "bdf";
# SUBSTEPS = exponential steps per driver interval (validation/numerics_convergence.py)
SOLVER = "expm"
SUBSTEPS = 16


def grain_formed(gm, floor=MASS_FLOOR):
    """1 where the grain mass is off the organ-mass floor, 0 where it sits on it; all 1 if
    the series never touches the floor (the grain exists from the start)."""
    gm = np.asarray(gm, dtype=float)
    on_floor = gm <= floor * (1.0 + 1e-12)
    return np.ones(len(gm)) if not on_floor.any() else (~on_floor).astype(float)


def panicle_share(dvs, at_heading, at_maturity):
    """The panicle's share of whole-plant transpiration at ORYZA development stage `dvs`
    (1 = heading/flowering, 2 = maturity): 0 before heading (the panicle is still inside the
    boot), then linear from `at_heading` to `at_maturity`. The two anchors are sourced
    (Sugimoto 1973, 1975; Roberts & Batchelor 1983; Takahashi & Kurata 2007); the linear
    form between them is an assumption (docs/PREREG_sourced_values.md 5.1)."""
    d = float(dvs)
    if not np.isfinite(d) or d < 1.0:
        return 0.0
    return float(at_heading + (at_maturity - at_heading) * min(d - 1.0, 1.0))

# Transpiring surface per kg DRY matter [m^2/kg], for the leaf/grain xylem split when M is
# fresh weight. Leaf 20 is the SLA reading of the original S. Grain 0.5116 = 2 x (1-0.78)/
# (1-0.14): the value at which splitting on DRY matter gives the split the original S
# gave on fresh mass at the shipped theta -- on the site-year ~6% of the leaf + grain
# stream to the grain from panicle formation to harvest, transpiration-weighted (8% from
# heading; nstem_leaf assumes 5%). NOT measured, and a FIXED literal on purpose:
# computing it from theta would tie the split to the organs' water again. Reading the
# original grain S=2 as per kg DRY instead sends ~14% to the grain and leaves 8 of 11
# Yamazaki grain BAFs out of reach at the L_Ph floor (docs/HANDOFF_mass_basis_refit.md §2).
S_XYLEM_DRY = {LEAF: 20.0, FRUIT: 0.5116}


def xylem_leaf_share(comps, M, M_basis=None, s_grain_dry=None):
    """Leaf share of the xylem stream that leaves the stem (the grain takes the rest).

    Split by transpiring surface x mass. For a FRESH M (M_basis='fresh') the split runs
    on DRY matter M*(1-theta), with the binding factor's own theta, and S_XYLEM_DRY --
    the crop-model SLA convention (leaf area = SLA x leaf dry weight), so an organ's
    water content does not set its transpiration share. M_basis None splits on M as
    given with Compartment.S (the basis-free placeholder logistic; Compartment.S is also
    the air-exchange surface).
    """
    if M_basis == "fresh":
        s4 = S_XYLEM_DRY[FRUIT] if s_grain_dry is None else float(s_grain_dry)
        A3 = S_XYLEM_DRY[LEAF] * M[LEAF] * (1.0 - comps[LEAF].theta)
        A4 = s4 * M[FRUIT] * (1.0 - comps[FRUIT].theta)
    else:
        A3 = comps[LEAF].S * M[LEAF]
        A4 = comps[FRUIT].S * M[FRUIT]
    return A3 / (A3 + A4) if (A3 + A4) > 0 else 0.5


def grain_dry_growth(comps, dM, M_basis=None):
    """Grain DRY-matter growth rate [kg/day] for the grain phloem stream Q_g = growth x T_C_Ph,
    whose coefficient is per kg of grain DRY matter. A fresh M (M_basis 'fresh') carries the
    grain's water, so its growth is converted with the binding factor's theta; None (the
    basis-free placeholder) is used as given."""
    if M_basis == "fresh":
        return dM[FRUIT] * (1.0 - comps[FRUIT].theta)
    return dM[FRUIT]


def phloem_flows(comps, dM, Qtp, T_C_Ph, phi, M_basis=None):
    """The two phloem streams from the leaf [L/day]: to the grain, the sap its DRY-matter
    growth draws (`grain_dry_growth` x T_C_Ph); to the root, the recirculating fraction phi
    of transpiration. Separate streams: the recirculation does not pass through the grain."""
    Q_g = max(grain_dry_growth(comps, dM, M_basis) * T_C_Ph, 0.0)
    Q_r = max(phi * Qtp, 0.0)
    return Q_g, Q_r


@dataclass
class PlantInputs:
    """Time-dependent external drivers (from soil model + growth model).

    Provide arrays on a common time grid `t` [day]:
        Cwo : soil pore-water free concentration  [ug/L]   (from HYDRUS/Phydrus)
        Qtp : transpiration stream                 [L/day]  (from HYDRUS/Phydrus)
        M   : tissue fresh mass, shape (len(t), 4) [kg]     (from growth model)
    Cubic/linear interpolants are built for use inside the ODE RHS.

    M_basis says what M is, for the one place that needs DRY matter -- the leaf/grain
    xylem split (`xylem_leaf_share`): 'fresh' = M is fresh weight, so the split takes its
    dry matter M*(1-theta) first; None (default) = split on M as given (the basis-free
    placeholder logistic).
    """
    t: np.ndarray
    Cwo: np.ndarray
    Qtp: np.ndarray
    M: np.ndarray            # shape (len(t), 4)
    leaf_loss: np.ndarray | None = None   # leaf senescence loss RATE [1/day] (opt; 0 if None)
    M_basis: str | None = None            # 'fresh' | None -- see the docstring
    dvs: np.ndarray | None = None         # crop development stage (ORYZA DVS); needed when the
                                          # panicle's xylem share follows the stage
    theta: np.ndarray | None = None       # organ water theta_fw(t) (len(t), 4) that made M fresh
                                          # (PREREG 7.1); None: the compartments' constant theta
    M_husk: np.ndarray | None = None      # panicle remainder (husk + rachis) DRY mass (PREREG 7.2)

    def __post_init__(self):
        if self.M_basis not in M_BASES:
            raise ValueError(f"unknown M_basis {self.M_basis!r}; expected one of {M_BASES}")
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
        # grain FORMATION gate gamma(t): 0 while the grain sits on the organ-mass floor (the
        # crop model has not formed the panicle), 1 from the day it leaves it. The grain's
        # PFAS influx is gated by this, so no solute loads a not-yet-formed organ
        # (DPU-consistent: import is tied to organ existence). A driver whose grain is never
        # on the floor (e.g. constant mass) has the grain present throughout (gate 1). A step,
        # not a ramp: the ramp that used to soften it was not converged
        # (validation/numerics_convergence.py); between daily driver nodes the gate is
        # linearly interpolated like every other driver.
        gate = grain_formed(self.M[:, FRUIT])
        self._grain_gate = interp1d(self.t, gate, **kw)
        self._dvs = (None if self.dvs is None
                     else interp1d(self.t, np.asarray(self.dvs, float), **kw))
        # one table of every driver column for the RHS: a single linear interpolation per
        # call instead of twelve scalar interp1d calls (identical values: same linear
        # interpolation, same linear extrapolation beyond the ends)
        dv = (np.full(len(self.t), np.nan) if self.dvs is None
              else np.asarray(self.dvs, float))
        self._tg = np.asarray(self.t, dtype=float)
        self._tab = np.column_stack([np.asarray(self.Cwo, float), np.asarray(self.Qtp, float),
                                     self.M, dM, ll, gate, dv])
        # the second table (PREREG 7.1/7.2): theta(t), the grain's DRY growth, the husk's dry
        # mass, growth and formation gate
        th = (np.full(self.M.shape, np.nan) if self.theta is None
              else np.asarray(self.theta, float).reshape(self.M.shape))
        gdry = (self.M[:, FRUIT] * (1.0 - th[:, FRUIT]) if self.theta is not None
                else np.full(len(self.t), np.nan))
        mh = (np.full(len(self.t), MASS_FLOOR) if self.M_husk is None
              else np.maximum(np.asarray(self.M_husk, float), MASS_FLOOR))
        self._tab2 = np.column_stack([th, np.gradient(gdry, self.t), mh,
                                      np.gradient(mh, self.t), grain_formed(mh)])

    def extras_at(self, t):
        """(theta[4], grain dry growth, husk dry mass, husk growth, husk gate) at time t."""
        tg = self._tg
        i = int(np.searchsorted(tg, t, side="right")) - 1
        i = min(max(i, 0), len(tg) - 2)
        w = (t - tg[i]) / (tg[i + 1] - tg[i])
        row = self._tab2[i] + w * (self._tab2[i + 1] - self._tab2[i])
        return row[0:4], row[4], row[5], row[6], min(1.0, max(0.0, row[7]))

    def drivers_at(self, t):
        """(Cwo, Qtp, M[4], dM[4], leaf_loss, grain_gate, dvs) at time t, by one linear lerp
        (dvs is NaN for a driver set without a development stage)."""
        tg = self._tg
        i = int(np.searchsorted(tg, t, side="right")) - 1
        i = min(max(i, 0), len(tg) - 2)
        w = (t - tg[i]) / (tg[i + 1] - tg[i])
        row = self._tab[i] + w * (self._tab[i + 1] - self._tab[i])
        return (row[0], row[1], row[2:6], row[6:10], max(row[10], 0.0),
                min(1.0, max(0.0, row[11])), row[12])

    def Cwo_(self, t):  return float(self._Cwo(t))
    def Qtp_(self, t):  return float(self._Qtp(t))
    def M_(self, t):    return np.array([float(f(t)) for f in self._M])
    def dM_(self, t):   return np.array([float(f(t)) for f in self._dM])
    def leaf_loss_(self, t):  return max(float(self._leaf_loss(t)), 0.0)
    def grain_gate_(self, t):  return min(1.0, max(0.0, float(self._grain_gate(t))))
    def dvs_(self, t):
        if self._dvs is None:
            raise ValueError("this driver set carries no 'dvs'; the stage-based panicle "
                             "share needs the crop's development stage")
        return float(self._dvs(t))


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
            c.f_prot * cmpd.K_prot + c.f_PL * cmpd.K_PL
            + c.f_cw * (cmpd.K_cw if c.K_cw is None else c.K_cw))
        for c in comps
    ])


def dry_binding(comps: list[Compartment], cmpd: Compound) -> np.ndarray:
    """X_k = f_prot K_prot + f_PL K_PL + f_cw K_cw [L/kg DRY]: the dry-matter part of B_k, which
    does not depend on the organ's water. B_k(t) = theta_k(t) + (1 - theta_k(t)) X_k when the
    water follows the development stage (PREREG 7.1)."""
    return np.array([c.f_prot * cmpd.K_prot + c.f_PL * cmpd.K_PL
                     + c.f_cw * (cmpd.K_cw if c.K_cw is None else c.K_cw) for c in comps])


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
    # linear equilibrium entry (opt-in; 0 removes it identically)
    j_eq = cmpd.lam_eq * (cmpd.r_eq * Cwo - Cw_root)
    return j_n + j_ed + j_carr + j_apo + j_eq


@dataclass
class RiceUptakeModel:
    env: Environment
    cmpd: Compound
    comps: list[Compartment]            # [root, stem, leaf, fruit]
    inputs: PlantInputs
    # phloem stream to the root as a fraction of transpiration [-] and phloem sap per kg of
    # grain dry-matter growth [L/kg]: no default -- the caller passes the sourced values
    # (model_api.PHI, T_C_PH); a model run without them raises
    phi: float | None = None
    T_C_Ph: float | None = None
    s_grain_dry: float | None = None    # override S_XYLEM_DRY[grain] (an unmeasured constant)
    # (panicle share at heading, at maturity): the leaf/grain xylem split follows the crop's
    # development stage (`panicle_share`); None keeps the transpiring-surface split
    lam_grain_dvs: tuple | None = None
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
    husk_phloem: bool = False           # the husk's dry growth draws phloem sap too (a band end)

    def __post_init__(self):
        if self.phi is None or self.T_C_Ph is None:
            raise ValueError("phi and T_C_Ph have no default: pass the sourced values "
                             "(model_api.PHI, model_api.T_C_PH)")
        # panicle remainder (husk + rachis) as a terminal state 4 (PREREG 7.2)
        self.husk = self.inputs.M_husk is not None
        self.HUSK = 4 if self.husk else None
        self.n_states = 5 if self.husk else 4
        self._staged = self.inputs.theta is not None
        if self.husk_phloem and not self.husk:
            raise ValueError("husk_phloem needs the husk pool (inputs.M_husk)")
        if self._staged and self.air is not None:
            raise ValueError("organ water by stage is not wired into plant-air exchange "
                             "(the neutral path keeps the constant water until amendment 2b)")

    def _binding(self, th):
        """B at this step: the stage's water when the inputs carry theta(t) (PREREG 7.1)."""
        if not self._staged:
            return binding_factors(self.comps, self.cmpd)
        return th + (1.0 - th) * dry_binding(self.comps, self.cmpd)

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
        Cwo, Qtp, M, dM, leaf_loss, gam, dvs = self.inputs.drivers_at(t)
        th, gdry, Mh, dMh, gh = self.inputs.extras_at(t)
        M = np.maximum(M, MASS_FLOOR)                # the drivers are floored already; a guard
        mu = dM / M                                  # growth-dilution rates [1/day]

        B = self._binding(th)                        # [L/kg]
        Cw = C[:4] / B                               # free aqueous conc [ug/L]

        # leaf/fruit xylem split: the panicle's sourced share of transpiration at the crop's
        # development stage (lam_grain_dvs); without it, the transpiring-surface split
        if self.lam_grain_dvs is not None:
            if self.inputs.dvs is None:
                self.inputs.dvs_(t)                 # raises: the stage-based share needs dvs
            f4 = panicle_share(dvs, *self.lam_grain_dvs)
            f3 = 1.0 - f4
        else:
            split = xylem_leaf_share(self.comps, M, self.inputs.M_basis, self.s_grain_dry)
            f3, f4 = split, 1.0 - split

        # phloem flows (carrier loading at leaf; NOT pH trap): the grain draws the sap its
        # dry-matter growth needs (T_C_Ph per kg DM); a fraction phi of transpiration
        # recirculates to the ROOT. The two are separate streams from the same leaf source.
        Q_g, Q_r = phloem_flows(self.comps, dM, Qtp, self.T_C_Ph, self.phi, self.inputs.M_basis)
        if self._staged:                             # the dry driver's own growth (PREREG 7.1)
            Q_g = max(gdry * self.T_C_Ph, 0.0)
        Q_h = max(dMh * self.T_C_Ph, 0.0) * gh if self.husk_phloem else 0.0
        # L_Ph_eff == cmpd.L_Ph unless a weak-electrolyte pH trap is configured
        C_Phl = (self.phloem_loading_factor() * Cw[LEAF]
                 + self.cmpd.g_ph * C[LEAF])          # free + lipid-bound [ug/L]

        g = [c.gamma for c in self.comps]
        dC = np.zeros(self.n_states)
        # the panicle's xylem: to the grain, or with the husk pool to the husk (no xylem flow
        # reaches the intact caryopsis; PREREG 7.2)
        # brown rice (the husk pool on) draws phloem in proportion to its own dry growth, so an
        # organ not yet growing draws nothing: its phloem stream is not gated (the mass-floor
        # gate opened a day early or late on DVS rounding at heading; PREREG section 8)
        wx = gh if self.husk else gam
        gq = 1.0 if self.husk else gam

        # root -> xylem loading: a fraction f_xy of the root FREE conc (TSCF;
        # endodermal barrier + binding) PLUS a B-independent lipid-bound term
        # g_xy*C_root (default 0). The SAME Cw_xyl is what the root loses and the
        # stem gains, so the root->stem xylem transfer remains mass-conserving.
        Cw_xyl = self.cmpd.f_xy * Cw[ROOT] + self.cmpd.g_xy * C[ROOT]

        # root
        jR = root_uptake(Cwo, Cw[ROOT], self.cmpd, self.env)
        dC[ROOT] = (jR
                    - (Qtp / M[ROOT]) * Cw_xyl
                    + (Q_r / M[ROOT]) * C_Phl
                    - g[ROOT] * C[ROOT] - mu[ROOT] * C[ROOT])
        # stem (xylem in from root at the loaded conc, xylem out to leaf+fruit)
        dC[STEM] = ((Qtp / M[STEM]) * (Cw_xyl - Cw[STEM])
                    - g[STEM] * C[STEM] - mu[STEM] * C[STEM])
        # leaf (xylem terminal; phloem source). The leaf supplies the WHOLE phloem
        # export: the grain stream (Q_g*C_Phl) plus the root recirculation (Q_r*C_Phl).
        # This closes the phloem mass balance (leaf loss = grain gain + root gain).
        # leaf senescence loss: the dead/shed leaf carries its PFAS away at the leaf
        # death rate (D/M_leaf), cancelling the spurious -mu*C concentration when the
        # leaf shrinks. 0 unless a senescing biomass driver (ORYZA) supplies the rate.
        # grain FORMATION gate: before the panicle sets (gam~0) the grain takes NO xylem/
        # phloem; its xylem share is rerouted to the leaf and its phloem stream is not
        # exported, so the balance still closes. gam->1 at formation (DPU-consistent).
        dC[LEAF] = ((f3 + (1.0 - wx) * f4) * (Qtp / M[LEAF]) * Cw[STEM]
                    - ((gq * Q_g + Q_r + Q_h) / M[LEAF]) * C_Phl
                    - g[LEAF] * C[LEAF] - mu[LEAF] * C[LEAF]
                    - leaf_loss * C[LEAF])
        # fruit/grain (small xylem in; phloem-dominated; terminal sink) -- gated by gam;
        # brown rice with the husk pool: the phloem only
        dC[FRUIT] = ((0.0 if self.husk else gam * f4 * (Qtp / M[FRUIT]) * Cw[STEM])
                     + gq * (Q_g / M[FRUIT]) * C_Phl
                     - g[FRUIT] * C[FRUIT] - mu[FRUIT] * C[FRUIT])
        if self.husk:                                # terminal, per kg DRY
            dC[self.HUSK] = ((wx * f4 * Qtp * Cw[STEM] + Q_h * C_Phl) / Mh
                             - (dMh / Mh) * C[self.HUSK])
        # plant-air exchange (gaseous uptake - volatilisation), opt-in and OFF by
        # default: for `air=None` this branch is skipped, so the PFAS path is
        # untouched. The stomatal pathway needs each shoot organ's share of the
        # transpiration stream, which is the same (f3, f4) split used above, with
        # the grain's share gated by its formation.
        if self.air is not None:
            dC[:4] += self.air.flux(C=C[:4], M=M, K_PW=B, comps=self.comps, Qtp=Qtp,
                                    xyl_share=(f3, gam * f4))
        return dC

    def affine(self, ts):
        """(A[m,4,4], b[m,4]) with rhs(t, C) = A(t) C + b(t) at the times `ts`, or None when
        the system is not affine (a carrier root entry). The same terms
        as `rhs`, vectorised; `tests/test_linear_expm.py` holds them equal to the rhs."""
        import linear_expm as LE
        kio = LE.root_coefficients(root_uptake, self.cmpd, self.env)
        if kio is None:
            return None
        k_in, k_out = kio
        row = LE.lerp_table(self.inputs._tg, self.inputs._tab, ts)
        row2 = LE.lerp_table(self.inputs._tg, self.inputs._tab2, ts)
        Cwo, Q = row[:, 0], row[:, 1]
        M = np.maximum(row[:, 2:6], MASS_FLOOR)
        dM = row[:, 6:10]
        ll = np.maximum(row[:, 10], 0.0)
        gam = np.clip(row[:, 11], 0.0, 1.0)
        dvs = row[:, 12]
        mu = dM / M
        B = (row2[:, 0:4] + (1.0 - row2[:, 0:4]) * dry_binding(self.comps, self.cmpd)
             if self._staged else np.broadcast_to(binding_factors(self.comps, self.cmpd),
                                                  (len(row), 4)))
        Mh, dMh = row2[:, 5], row2[:, 6]
        gh = np.clip(row2[:, 7], 0.0, 1.0)
        g = np.array([c.gamma for c in self.comps])
        if self.lam_grain_dvs is not None:
            if self.inputs.dvs is None:
                self.inputs.dvs_(0.0)
            f4 = np.array([panicle_share(d, *self.lam_grain_dvs) for d in dvs])
        else:
            f4 = 1.0 - np.array([xylem_leaf_share(self.comps, Mi, self.inputs.M_basis,
                                                  self.s_grain_dry) for Mi in M])
        f3 = 1.0 - f4
        if self._staged:
            dMg = row2[:, 4]
        else:
            dMg = (dM[:, FRUIT] * (1.0 - self.comps[FRUIT].theta)
                   if self.inputs.M_basis == "fresh" else dM[:, FRUIT])
        Q_g = np.maximum(dMg * self.T_C_Ph, 0.0)
        Q_r = np.maximum(self.phi * Q, 0.0)
        Q_h = (np.maximum(dMh * self.T_C_Ph, 0.0) * gh if self.husk_phloem
               else np.zeros(len(row)))
        p = self.phloem_loading_factor() / B[:, LEAF] + self.cmpd.g_ph
        x = self.cmpd.f_xy / B[:, ROOT] + self.cmpd.g_xy
        wx = gh if self.husk else gam
        gq = np.ones_like(gam) if self.husk else gam   # brown rice: phloem not gated (rhs)
        n = self.n_states
        A = np.zeros((len(row), n, n))
        b = np.zeros((len(row), n))
        A[:, ROOT, ROOT] = -k_out / B[:, ROOT] - Q * x / M[:, ROOT] - g[ROOT] - mu[:, ROOT]
        A[:, ROOT, LEAF] = Q_r * p / M[:, ROOT]
        b[:, ROOT] = k_in * Cwo
        A[:, STEM, ROOT] = Q * x / M[:, STEM]
        A[:, STEM, STEM] = -Q / (M[:, STEM] * B[:, STEM]) - g[STEM] - mu[:, STEM]
        A[:, LEAF, STEM] = (f3 + (1.0 - wx) * f4) * Q / (M[:, LEAF] * B[:, STEM])
        A[:, LEAF, LEAF] = (-(gq * Q_g + Q_r + Q_h) * p / M[:, LEAF] - g[LEAF] - mu[:, LEAF] - ll)
        A[:, FRUIT, STEM] = 0.0 if self.husk else gam * f4 * Q / (M[:, FRUIT] * B[:, STEM])
        A[:, FRUIT, LEAF] = gq * Q_g * p / M[:, FRUIT]
        A[:, FRUIT, FRUIT] = -g[FRUIT] - mu[:, FRUIT]
        if self.husk:
            A[:, 4, STEM] = wx * f4 * Q / (Mh * B[:, STEM])
            A[:, 4, LEAF] = Q_h * p / Mh
            A[:, 4, 4] = -dMh / Mh
        if self.air is not None and self.air.K_AW > 0.0:
            # plant-air exchange is affine in C too (volatilisation ~ C, gaseous uptake from a
            # fixed air concentration): read it off `AirExchange.flux` itself, time by time
            for i in range(len(row)):
                kw = dict(M=M[i], K_PW=B[i], comps=self.comps, Qtp=Q[i],
                          xyl_share=(f3[i], gam[i] * f4[i]))
                Aa, ba = LE.affine_parts(lambda _t, C: self.air.flux(C=C, **kw), 0.0, 4)
                A[i, :4, :4] += Aa
                b[i, :4] += ba
        return A, b

    def solve(self, t_eval: np.ndarray, C0: np.ndarray | None = None):
        """Integrate the system. Affine (the shipped PFAS and neutral paths, air exchange
        included): exponential stepping (`linear_expm`, converged in
        validation/numerics_convergence.py); otherwise BDF. Returns a solution object with
        the scipy fields the callers read (t, y, success, sol)."""
        import linear_expm as LE
        if C0 is None:
            C0 = np.zeros(self.n_states)
        if SOLVER == "expm" and self.affine(np.asarray(t_eval[:1], float)) is not None:
            return LE.solve_burden(self.affine, self._burden_scale, t_eval, self.n_states, C0,
                                   substeps=SUBSTEPS)
        t_span = (float(t_eval[0]), float(t_eval[-1]))
        D0, _ = self._burden_scale(np.asarray(t_eval[:1], float))
        sol = solve_ivp(LE.burden_rhs(self.rhs, self._burden_scale), t_span, C0 * D0[0],
                        t_eval=t_eval, method="BDF", rtol=1e-6, atol=1e-9, dense_output=True)
        return LE.BurdenSolution(sol, self._burden_scale)

    def _burden_scale(self, ts):
        """(D, dD) of the burden coordinates `solve` integrates in: the grain's fresh mass and
        the husk's dry mass (and their growth) on those rows, 1 (0) elsewhere
        (`linear_expm.burden_affine`)."""
        import linear_expm as LE
        inp = self.inputs
        ts = np.asarray(ts, float)
        row = LE.lerp_table(inp._tg, inp._tab, ts)
        D, dD = np.ones((len(ts), self.n_states)), np.zeros((len(ts), self.n_states))
        D[:, FRUIT], dD[:, FRUIT] = np.maximum(row[:, 5], MASS_FLOOR), row[:, 9]
        if self.husk:
            row2 = LE.lerp_table(inp._tg, inp._tab2, ts)
            D[:, self.HUSK], dD[:, self.HUSK] = row2[:, 5], row2[:, 6]
        return D, dD

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


def _logistic(t, M0, Mmax, k, t0):
    """Logistic curve (a helper for synthetic test drivers; no model value)."""
    return Mmax / (1.0 + (Mmax / M0 - 1.0) * np.exp(-k * (t - t0)))
