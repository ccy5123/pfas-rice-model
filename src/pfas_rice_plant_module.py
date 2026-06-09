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
    # speciation (PFAS: fully dissociated)
    fd: float = 1.0
    fn: float = 0.0


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

    def Cwo_(self, t):  return float(self._Cwo(t))
    def Qtp_(self, t):  return float(self._Qtp(t))
    def M_(self, t):    return np.array([float(f(t)) for f in self._M])
    def dM_(self, t):   return np.array([float(f(t)) for f in self._dM])


# ----------------------------------------------------------------------------
# Model functions
# ----------------------------------------------------------------------------
def binding_factors(comps: list[Compartment], cmpd: Compound) -> np.ndarray:
    """B_k = theta_k + f_prot K_prot + f_PL K_PL + f_cw K_cw   [L/kg]  (Eq. binding)."""
    return np.array([
        c.theta + c.f_prot * cmpd.K_prot + c.f_PL * cmpd.K_PL + c.f_cw * cmpd.K_cw
        for c in comps
    ])


def _ghk_factor(N: float) -> float:
    """N / (exp(N) - 1), with the removable singularity at N=0 handled."""
    if abs(N) < 1e-8:
        return 1.0
    return N / np.expm1(N)


def root_uptake(Cwo: float, Cw_root: float, cmpd: Compound, env: Environment) -> float:
    """Mass-specific root membrane uptake j_R [ug/(day kg)]  (Eq. JR_pfas).

    Hybrid: ionic electrodiffusion (GHK) + saturable carrier (Michaelis-Menten).
    For PFAS the neutral term is dropped (fn ~ 0).
    """
    N = env.N
    eN = np.exp(N)
    g = _ghk_factor(N)
    # ionic electrodiffusion (membrane + anion channel + aquaporin, lumped in kappa_d)
    j_ed = cmpd.kappa_d * g * (cmpd.fd * Cwo - cmpd.fd * eN * Cw_root)
    # carrier-mediated (active/facilitated), net influx - efflux
    j_carr = (cmpd.Vmax_in * Cwo / (cmpd.Km_in + Cwo)
              - cmpd.Vmax_out * Cw_root / (cmpd.Km_out + Cw_root))
    # optional neutral passive term (negligible for PFAS): cmpd.fn * ...
    return j_ed + j_carr


@dataclass
class RiceUptakeModel:
    env: Environment
    cmpd: Compound
    comps: list[Compartment]            # [root, stem, leaf, fruit]
    inputs: PlantInputs
    phi: float = 0.1                    # phloem recirculation fraction to roots [-]
    T_C_Ph: float = 10.0                # phloem flux per unit grain dry mass [L/kg]

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
        C_Phl = self.cmpd.L_Ph * Cw[LEAF]                    # [ug/L]

        g = [c.gamma for c in self.comps]
        dC = np.zeros(4)

        # root
        jR = root_uptake(Cwo, Cw[ROOT], self.cmpd, self.env)
        dC[ROOT] = (jR
                    - (Qtp / M[ROOT]) * Cw[ROOT]
                    + self.phi * (Q_Phl / M[ROOT]) * C_Phl
                    - g[ROOT] * C[ROOT] - mu[ROOT] * C[ROOT])
        # stem (xylem in from root, xylem out to leaf+fruit)
        dC[STEM] = ((Qtp / M[STEM]) * (Cw[ROOT] - Cw[STEM])
                    - g[STEM] * C[STEM] - mu[STEM] * C[STEM])
        # leaf (xylem terminal; phloem source)
        dC[LEAF] = (f3 * (Qtp / M[LEAF]) * Cw[STEM]
                    - (Q_Phl / M[LEAF]) * C_Phl
                    - g[LEAF] * C[LEAF] - mu[LEAF] * C[LEAF])
        # fruit/grain (small xylem in; phloem-dominated; terminal sink)
        dC[FRUIT] = (f4 * (Qtp / M[FRUIT]) * Cw[STEM]
                     + (Q_Phl / M[FRUIT]) * C_Phl
                     - g[FRUIT] * C[FRUIT] - mu[FRUIT] * C[FRUIT])
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
    # PFOA-like placeholder compound (values are illustrative, NOT calibrated)
    cmpd = Compound(
        name="PFOA",
        K_prot=50.0, K_PL=100.0, K_cw=20.0,        # [L/kg]
        kappa_d=0.5,                                # [L/(day kg)]
        Vmax_in=20.0, Km_in=5.0,                    # carrier influx (must overcome anion exclusion)
        Vmax_out=8.0, Km_out=5.0,                   # carrier efflux
        L_Ph=0.05,
    )
    comps = [
        Compartment("root",  theta=0.70, f_prot=0.02, f_PL=0.010, f_cw=0.30),
        Compartment("stem",  theta=0.80, f_prot=0.01, f_PL=0.005, f_cw=0.08),
        Compartment("leaf",  theta=0.80, f_prot=0.03, f_PL=0.020, f_cw=0.04, S=20.0),
        Compartment("grain", theta=0.15, f_prot=0.08, f_PL=0.010, f_cw=0.10, S=2.0),
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
    # grain accumulation check (terminal sink: should keep rising while it fills)
    print(f"\ngrain BAF ordering vs vegetative (expect root highest, grain trapped by protein):")
    print(f"  root/grain BAF ratio = {baf[ROOT]/baf[FRUIT]:.2f}")

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
