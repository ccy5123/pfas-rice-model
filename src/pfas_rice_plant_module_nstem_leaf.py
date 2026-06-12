"""
N-segment stem + explicit leaf ("redistributed shoot") plant uptake module
==========================================================================

Fixes the above-ground *over-translocation* the Tang 2026 out-of-sample test
exposed in the single-straw 4-compartment core (``pfas_rice_plant_module_4pool_surf``):

  * the **stem was an empty pass-through** -- ``dC_stem = (Q/M_stem)(Cw_xyl - Cw_stem)``
    equilibrates the stem to the (f_xy-discounted) xylem free conc and then
    *re-exports* it up the ascending xylem, so the stem never accumulates:
    ``TF_stem = f_xy*B_stem/B_root << 1`` (Tang stalk TF ~ 0.6-1.5);
  * the **leaf was the sole xylem terminal**, so it integrated the WHOLE
    transpiration stream and ran away (it held ~81% of the plant burden; Tang
    leaf TF ~ 0.7-1.7, the single-straw model gave 3-13, worst for the mobile
    ether GenX).

Both are the same structural defect: with one straw compartment ALL the shoot's
xylem-delivered solute can only pile into the leaf, because a serial *mixer*
stem re-exports whatever its transpiration concentrates.  This module cures it by

  (1) **resolving the stem into N serial segments** stacked by height, and
  (2) **redistributing the shoot loading by the transpiration-stream-concentration
      mechanism applied to EVERY shoot organ, not just the leaf**: where canopy
      water evaporates (leaf blades AND stem sheaths/culm AND, a little, the
      panicle), the non-volatile anion it carried is deposited and *retained*
      (it is at a dead-end of the flow path -- there is no route back up against
      the stream).  Each stem segment is therefore a partial transpiration
      *terminal* exactly like the leaf, so the shoot burden is reapportioned
      root -> stem -> leaf -> grain instead of piling into the leaf.

Mechanism (transpiration deposition + retention)
------------------------------------------------
The root loads the xylem at ``Cw_xyl = f_xy*Cw_root (+ g_xy*C_root)`` and exports
the whole stream ``Q*Cw_xyl``.  The stream ascends; organ k transpires ``lam_k*Q``
of water and the solute that water carried, ``lam_k*Q*Cw_xyl``, is deposited
there.  A fraction ``retention`` of each deposit stays in the organ (irreversible
cell-wall/apoplast sequestration at the evaporative terminus -- the terminal
accumulator), the rest is carried on as residual xylem to the grain.  The split
``tau_1..tau_N (stem) + lam_leaf + lam_grain = 1`` closes exactly, so

    sum of shoot deposits  ==  root xylem export  ==  Q*Cw_xyl     (mass-conserving).

``retention=1`` makes every organ a full terminal (the cleanest "redistribution");
``retention=0`` sends all of it to the grain (degenerate).  How transpiration
terminates across blade / sheath+culm / panicle is a CROP-ARCHITECTURE quantity,
independent of PFAS; together with ``retention`` it is the new structural lever.
See ``validation/tang2026_nstem_validation.py`` for the Tang re-check and the
split / retention sensitivity sweep.

Note this is a *different and complementary* generalization of the stem from
``pfas_rice_plant_module_nstem`` (the serial advective *mixer*, which keeps the
solute in reversible balance with the upward-concentrating xylem and is used for
the Yamazaki *within-stem vertical gradient* S18/S19).  Here the solute is
RETAINED at each evaporative terminus -- the mechanism the single-straw model was
missing and that Tang's stalk/leaf/grain split needs.

Compartments: root(0), stem_1..stem_N(1..N, bottom->top), leaf(N+1), grain(N+2).
Reuses the basis-A binding, GHK+carrier root influx and Compound/Environment of
``pfas_rice_plant_module_4pool_surf``; only the shoot ODE is generalized.

Mass balance: for gamma=0 the sole source is the root membrane flux M_root*j_R;
every xylem deposit and phloem transfer telescopes (tests/test_nstem_leaf.py).
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

from pfas_rice_plant_module_4pool_surf import (
    Environment, Compound, Compartment, binding_factors, root_uptake,
)


@dataclass
class PlantInputsNL:
    """Time-dependent drivers for the N-stem + leaf model.

    M : tissue fresh mass, shape (len(t), N+3) with columns ordered
        [root, stem_1..stem_N, leaf, grain].
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
class NStemLeafModel:
    """Redistributed-shoot uptake model (transpiration deposition + retention).

    comps : [root, stem_1..stem_N, leaf, grain]
    tau   : per-stem-segment transpiration draw-off (fraction of Q), len N
    lam_leaf, lam_grain : leaf / panicle transpiration draw-off (fraction of Q)
        Constraint: sum(tau) + lam_leaf + lam_grain == 1.
    retention : fraction of each organ's transpiration-deposited solute retained
        (terminal); 1-retention is carried on to the grain as residual xylem.
    """
    env: Environment
    cmpd: Compound
    comps: list
    inputs: PlantInputsNL
    tau: np.ndarray
    lam_leaf: float = 0.50
    lam_grain: float = 0.05
    retention: float = 1.0
    phi: float = 0.1                  # phloem recirculation fraction to root
    T_C_Ph: float = 10.0              # phloem flux per unit grain dry-mass gain [L/kg]

    def __post_init__(self):
        self.N = len(self.comps) - 3          # number of stem segments
        self.ROOT = 0
        self.LEAF = self.N + 1
        self.GRAIN = self.N + 2
        self.tau = np.asarray(self.tau, dtype=float)
        assert len(self.tau) == self.N, "tau must have one entry per stem segment"
        tot = float(self.tau.sum()) + self.lam_leaf + self.lam_grain
        assert abs(tot - 1.0) < 1e-9, (
            f"transpiration split must sum to 1 (sum(tau)+lam_leaf+lam_grain={tot:.6f})")
        assert self.tau.min() >= 0 and self.lam_leaf >= 0 and self.lam_grain >= 0
        assert 0.0 <= self.retention <= 1.0

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

        leaf, grain, r = self.LEAF, self.GRAIN, self.retention
        # phloem: loaded at the leaf (free + optional lipid-bound), feeds grain + root recirc.
        Q_Phl = max(dM[grain] * self.T_C_Ph + self.phi * Q, 0.0)
        C_Phl = self.cmpd.L_Ph * Cw[leaf] + self.cmpd.g_ph * C[leaf]

        # --- root: uptake, full xylem export (f_xy free + lipid g_xy bound), phloem return ---
        jR = root_uptake(Cwo, Cw[self.ROOT], self.cmpd, self.env)
        Cw_xyl = self.cmpd.f_xy * Cw[self.ROOT] + self.cmpd.g_xy * C[self.ROOT]
        xyl_flux = Q * Cw_xyl                          # total solute loaded into the xylem
        dC[self.ROOT] = (jR
                         - xyl_flux / M[self.ROOT]
                         + self.phi * (Q_Phl / M[self.ROOT]) * C_Phl
                         - g[self.ROOT] * C[self.ROOT] - mu[self.ROOT] * C[self.ROOT])

        # --- stem segments: each retains its transpiration-deposited solute ---
        for s in range(1, self.N + 1):
            dC[s] = (r * self.tau[s - 1] * xyl_flux / M[s]
                     - g[s] * C[s] - mu[s] * C[s])

        # --- leaf: retains its transpiration deposit; phloem source (grain + root recirc) ---
        dC[leaf] = (r * self.lam_leaf * xyl_flux / M[leaf]
                    - (1.0 + self.phi) * (Q_Phl / M[leaf]) * C_Phl
                    - g[leaf] * C[leaf] - mu[leaf] * C[leaf])

        # --- grain: own panicle deposit + the non-retained residual xylem + phloem ---
        residual = (self.lam_grain + (1.0 - r) * (self.tau.sum() + self.lam_leaf)) * xyl_flux
        dC[grain] = (residual / M[grain]
                     + (Q_Phl / M[grain]) * C_Phl
                     - g[grain] * C[grain] - mu[grain] * C[grain])
        return dC

    def solve(self, t_eval, C0=None):
        if C0 is None:
            C0 = np.zeros(len(self.comps))
        return solve_ivp(self.rhs, (float(t_eval[0]), float(t_eval[-1])), C0,
                         t_eval=t_eval, method="BDF", rtol=1e-6, atol=1e-9,
                         dense_output=True)

    def stem_aggregate(self, C, t):
        """Mass-weighted mean stem (stalk) concentration over the N segments."""
        M = self.inputs.M_(t)
        seg = slice(1, self.N + 1)
        Ms = M[seg]
        return float(np.sum(np.asarray(C)[seg] * Ms) / np.sum(Ms))


def make_stem_leaf_compartments(N, stem_kw, root_kw, leaf_kw, grain_kw):
    """Build [root, stem_1..stem_N, leaf, grain] with identical stem composition."""
    comps = [Compartment("root", **root_kw)]
    comps += [Compartment(f"stem{s}", **stem_kw) for s in range(1, N + 1)]
    comps += [Compartment("leaf", **leaf_kw), Compartment("grain", **grain_kw)]
    return comps


def default_tau(N, stem_transp_frac):
    """Per-segment stem draw-off shape (decreasing with height) scaled to a total.

    The lower segments carry the larger sheath/culm transpiring surface, so the
    draw-off decreases upward; the entries sum to ``stem_transp_frac`` (the total
    fraction of canopy transpiration that terminates on the stalk)."""
    shape = np.linspace(1.3, 0.7, N)              # decreasing bottom->top
    return shape / shape.sum() * float(stem_transp_frac)


def split_from_stem_frac(N, stem_transp_frac, lam_grain=0.05):
    """Convenience: (tau, lam_leaf, lam_grain) from a single stalk-transpiration
    fraction, with the leaf taking the remainder.  Closes to sum 1."""
    tau = default_tau(N, stem_transp_frac)
    lam_leaf = 1.0 - stem_transp_frac - lam_grain
    return tau, lam_leaf, lam_grain


if __name__ == "__main__":
    # smoke test against the measured forcings (PFOA), reported as TF
    import os, sys, json
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import forcing_rice as fr
    import growth_rice as gr

    ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    PAR = json.load(open(os.path.join(ROOT_DIR, "params", "parameters.json")))
    comp = PAR["tissue_composition_recommended"]
    carr = PAR["carrier_MichaelisMenten"]
    _kw = lambda d: dict(theta=d["theta_fw"], f_prot=d["f_prot"], f_PL=d["f_PL"], f_cw=d["f_cw"])

    N, SEASON = 4, 150.0
    t = np.linspace(0.0, SEASON, 481)
    Qtp = fr.Q_TP(t, SEASON)
    b = gr.organ_biomass(t, SEASON)
    M = np.column_stack(
        [np.maximum(b["root"], 1e-9)]
        + [np.maximum(b["stem"] / N, 1e-9)] * N
        + [np.maximum(b["leaf"], 1e-9), np.maximum(b["grain"], 1e-9)])
    inputs = PlantInputsNL(t=t, Cwo=np.full_like(t, 1.0), Qtp=Qtp, M=M)
    comps = make_stem_leaf_compartments(
        N, _kw(comp["stem"]), _kw(comp["root"]), _kw(comp["leaf"]), _kw(comp["grain_brown"]))

    c = next(x for x in PAR["congeners"] if x["name"] == "PFOA")
    cmpd = Compound(name="PFOA", K_prot=c["K_prot_Lkg"], K_PL=c["K_PL_Lkg"],
                    K_cw=c["K_cw_wholecw_Lkg"]["stem"], kappa_d=2.0,
                    Vmax_in=carr["Vmax_in"], Km_in=carr["Km_in"],
                    Vmax_out=carr["Vmax_out"], Km_out=carr["Km_out"],
                    L_Ph=0.01, f_xy=c["f_xy_recommended"])
    tau, lam_leaf, lam_grain = split_from_stem_frac(N, 0.45, lam_grain=0.05)
    for ret in (1.0, 0.6):
        m = NStemLeafModel(env=Environment(), cmpd=cmpd, comps=comps, inputs=inputs,
                           tau=tau, lam_leaf=lam_leaf, lam_grain=lam_grain, retention=ret)
        Cend = m.solve(t).y[:, -1]
        root = Cend[0]
        stalk = m.stem_aggregate(Cend, t[-1])
        print(f"PFOA (N=4, stem_transp=0.45, retention={ret}): "
              f"TF stalk={stalk/root:.2f} (Tang 1.45)  "
              f"leaf={Cend[m.LEAF]/root:.2f} (1.66)  grain={Cend[m.GRAIN]/root:.2f} (0.95)")
