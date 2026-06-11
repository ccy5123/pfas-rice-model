"""
N-segment ("multi-height stem") plant uptake module
====================================================

Resolves the single mass-weighted STEM compartment into N serial stem(+leaf)
segments stacked by height, so the model can represent the *vertical* tissue
gradient that Yamazaki 2023 (Table S18/S19) measured (rice stem 0-20/20-40/
40-60/>60 cm).  This is the GAP-B "open modeling item": with one well-mixed
straw compartment a *monotone* f_xy could not reproduce the data and the W2 fit
had to inflate long-chain f_xy.  The hypothesis under test here is that a serial
stem with transpiration draw-off + radial (tissue) exchange + growth dilution
reproduces BOTH the short-chain upward gradient AND the long-chain flat/down
gradient with a single MONOTONE f_xy.

Mechanism (why the gradient direction flips with chain length)
-------------------------------------------------------------
For each segment the quasi-steady tissue conc is
    C_s ≈ Q_{s-1}·Cw_in / ( M_s·μ_s  +  Q_s/B_s )
i.e. a competition between advective throughput (Q_s/B_s) and growth
accumulation (M_s·μ_s):
  * short chain (low B_s): Q_s/B_s dominates → the segment passes solute on;
    transpiration water loss (Q_s < Q_{s-1}) concentrates Cw upward → UP gradient.
  * long chain (high B_s): M_s·μ_s dominates → the lower segments strip the
    ascending stream, starving the upper ones → DOWN/flat gradient.

Compartments: root(0), stem_1..stem_N(1..N, bottom→top), grain(N+1).
Leaves are folded into each stem segment (Yamazaki lumps "stem incl. leaves").
Xylem is serial root→s1→…→sN→grain with transpiration τ_s drawn at each segment
(water leaves, the non-volatile anion stays → upward concentration). Grain is
phloem-fed from the top segment plus the small residual xylem.

Reuses the basis-A binding, GHK+carrier root influx and Compound/Environment of
``pfas_rice_plant_module_4pool_surf``; only the shoot ODE is generalized.

Mass balance: for γ=0 the sole source is the root membrane flux M_root·j_R; every
internal xylem/phloem transfer telescopes (verified in the test).
"""
from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

from pfas_rice_plant_module_4pool_surf import (
    Environment, Compound, Compartment, binding_factors, root_uptake,
)


@dataclass
class PlantInputsN:
    """Time-dependent drivers for the N-segment model.

    M : tissue fresh mass, shape (len(t), n_comp) with n_comp = N+2
        columns ordered [root, stem_1..stem_N, grain].
    """
    t: np.ndarray
    Cwo: np.ndarray
    Qtp: np.ndarray
    M: np.ndarray

    def __post_init__(self):
        self.M = np.asarray(self.M, dtype=float)
        self.n_comp = self.M.shape[1]
        kw = dict(kind="linear", bounds_error=False, fill_value="extrapolate")
        self._Cwo = interp1d(self.t, self.Cwo, **kw)
        self._Qtp = interp1d(self.t, self.Qtp, **kw)
        self._M = [interp1d(self.t, self.M[:, k], **kw) for k in range(self.n_comp)]
        dM = np.gradient(self.M, self.t, axis=0)
        self._dM = [interp1d(self.t, dM[:, k], **kw) for k in range(self.n_comp)]

    def Cwo_(self, t): return float(self._Cwo(t))
    def Qtp_(self, t): return float(self._Qtp(t))
    def M_(self, t):   return np.array([float(f(t)) for f in self._M])
    def dM_(self, t):  return np.array([float(f(t)) for f in self._dM])


@dataclass
class NStemModel:
    env: Environment
    cmpd: Compound
    comps: list[Compartment]          # [root, stem_1..stem_N, grain]
    inputs: PlantInputsN
    tau: np.ndarray                   # transpiration fraction per stem segment (sum=1)
    phi: float = 0.1                  # phloem recirculation fraction to root
    T_C_Ph: float = 10.0              # phloem flux per unit grain dry-mass gain [L/kg]

    def __post_init__(self):
        self.N = len(self.comps) - 2          # number of stem segments
        self.ROOT = 0
        self.GRAIN = self.N + 1
        self.tau = np.asarray(self.tau, dtype=float)
        assert len(self.tau) == self.N, "tau must have one entry per stem segment"
        # sum(tau) is the fraction of transpiration drawn across the stem; the
        # residual (1 - sum) is the xylem flow that continues to the grain.
        assert 0.0 < self.tau.sum() <= 1.0 + 1e-9, "sum(tau) must be in (0, 1]"

    def rhs(self, t, C):
        Cwo = self.inputs.Cwo_(t)
        Q = self.inputs.Qtp_(t)
        M = np.maximum(self.inputs.M_(t), 1e-12)
        dM = self.inputs.dM_(t)
        mu = dM / M
        B = binding_factors(self.comps, self.cmpd)
        Cw = C / B
        g = np.array([c.gamma for c in self.comps])
        dC = np.zeros_like(C)

        top = self.N                                   # index of top stem segment
        Q_Phl = max(dM[self.GRAIN] * self.T_C_Ph + self.phi * Q, 0.0)
        C_Phl = self.cmpd.L_Ph * Cw[top]               # phloem loads at the top segment

        # --- root: uptake, xylem export (f_xy), phloem return ---
        jR = root_uptake(Cwo, Cw[self.ROOT], self.cmpd, self.env)
        dC[self.ROOT] = (jR
                         - (Q / M[self.ROOT]) * self.cmpd.f_xy * Cw[self.ROOT]
                         + self.phi * (Q_Phl / M[self.ROOT]) * C_Phl
                         - g[self.ROOT] * C[self.ROOT] - mu[self.ROOT] * C[self.ROOT])

        # --- serial stem segments with transpiration draw-off ---
        Cw_in = self.cmpd.f_xy * Cw[self.ROOT]         # stream entering segment 1
        Qin = Q
        for s in range(1, self.N + 1):
            Qout = Qin - self.tau[s - 1] * Q           # water transpired here (solute stays)
            dC[s] = ((Qin * Cw_in - Qout * Cw[s]) / M[s]
                     - g[s] * C[s] - mu[s] * C[s])
            Cw_in = Cw[s]                              # passes its free conc upward
            Qin = Qout
        # top segment additionally exports the full phloem (grain + root recirc.)
        dC[top] -= (1.0 + self.phi) * (Q_Phl / M[top]) * C_Phl

        # --- grain: residual xylem (Qin*Cw_in = Q_out[N]*Cw[top]) + phloem load ---
        dC[self.GRAIN] = ((Qin * Cw_in) / M[self.GRAIN]
                          + (Q_Phl / M[self.GRAIN]) * C_Phl
                          - g[self.GRAIN] * C[self.GRAIN] - mu[self.GRAIN] * C[self.GRAIN])
        return dC

    def solve(self, t_eval, C0=None):
        if C0 is None:
            C0 = np.zeros(len(self.comps))
        sol = solve_ivp(self.rhs, (float(t_eval[0]), float(t_eval[-1])), C0,
                        t_eval=t_eval, method="BDF", rtol=1e-6, atol=1e-9,
                        dense_output=True)
        return sol


def make_stem_compartments(N, stem_kw, root_kw, grain_kw):
    """Build [root, stem_1..stem_N, grain] with identical stem composition."""
    comps = [Compartment("root", **root_kw)]
    comps += [Compartment(f"stem{s}", **stem_kw) for s in range(1, N + 1)]
    comps += [Compartment("grain", **grain_kw)]
    return comps


@dataclass
class NStemKineticModel(NStemModel):
    """N-segment stem with a KINETIC radial xylem<->tissue exchange (sorbing column).

    The equilibrium :class:`NStemModel` lets each segment export the ascending
    xylem at its own free conc (instantaneous radial equilibrium); at realistic
    biomass that gives a ~chain-length-independent upward gradient.  Here the
    xylem and tissue are DECOUPLED: within each segment the (quasi-steady) xylem
    pool exchanges radially with the tissue at a finite mass-specific conductance
    ``k_rad`` [L/(day kg)]:

        quasi-steady xylem:  Q_in*Cw_in = Q_out*Cw_xyl + k_rad*M_s*(Cw_xyl - Cw_s)
            =>  Cw_xyl = (Q_in*Cw_in + k_rad*M_s*Cw_s) / (Q_out + k_rad*M_s)
        tissue:              dC_s/dt = k_rad*(Cw_xyl - Cw_s) - mu_s*C_s

    ``k_rad -> inf`` recovers :class:`NStemModel`.

    HONEST STATUS (validation/nstem_gradient_check.py): finite ``k_rad`` lowers
    the long-chain (high-B) top/bottom ratio somewhat, but across all ``k_rad``
    the predicted gradient range stays ~2.4-4.9, far short of the observed range
    (PFBA top/bot 7.4 down to PFUnDA 0.66).  The reason: here the tissue free
    conc stays in (reversible) balance with the upward-concentrating xylem, so
    the BOUND fraction also concentrates upward.  Reproducing the long-chain
    REVERSAL needs the high-B bound fraction to be delivery-limited and slow to
    re-release -- i.e. IRREVERSIBLE/hysteretic sorption (a sink that strips the
    stream low down), an open elaboration.  This class is therefore a modest,
    mass-conserving improvement, not a full gradient fit.
    """
    k_rad: float = 0.5                # radial xylem<->tissue conductance [L/(day kg)]

    def rhs(self, t, C):
        Cwo = self.inputs.Cwo_(t)
        Q = self.inputs.Qtp_(t)
        M = np.maximum(self.inputs.M_(t), 1e-12)
        dM = self.inputs.dM_(t)
        mu = dM / M
        B = binding_factors(self.comps, self.cmpd)
        Cw = C / B
        g = np.array([c.gamma for c in self.comps])
        dC = np.zeros_like(C)

        top = self.N
        Q_Phl = max(dM[self.GRAIN] * self.T_C_Ph + self.phi * Q, 0.0)
        C_Phl = self.cmpd.L_Ph * Cw[top]

        jR = root_uptake(Cwo, Cw[self.ROOT], self.cmpd, self.env)
        dC[self.ROOT] = (jR
                         - (Q / M[self.ROOT]) * self.cmpd.f_xy * Cw[self.ROOT]
                         + self.phi * (Q_Phl / M[self.ROOT]) * C_Phl
                         - g[self.ROOT] * C[self.ROOT] - mu[self.ROOT] * C[self.ROOT])

        Cw_in = self.cmpd.f_xy * Cw[self.ROOT]
        Qin = Q
        for s in range(1, self.N + 1):
            Qout = Qin - self.tau[s - 1] * Q
            krm = self.k_rad * M[s]                       # k_rad * M_s  [L/day]
            Cw_xyl = (Qin * Cw_in + krm * Cw[s]) / (Qout + krm)   # quasi-steady xylem
            dC[s] = (self.k_rad * (Cw_xyl - Cw[s])        # radial influx into tissue
                     - g[s] * C[s] - mu[s] * C[s])
            Cw_in = Cw_xyl                                 # the (stripped) stream goes up
            Qin = Qout
        dC[top] -= (1.0 + self.phi) * (Q_Phl / M[top]) * C_Phl

        dC[self.GRAIN] = ((Qin * Cw_in) / M[self.GRAIN]
                          + (Q_Phl / M[self.GRAIN]) * C_Phl
                          - g[self.GRAIN] * C[self.GRAIN] - mu[self.GRAIN] * C[self.GRAIN])
        return dC
