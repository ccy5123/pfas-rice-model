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

Optional TWO-POOL root (``k_seq`` > 0) -- the structural merge
---------------------------------------------------------------
Setting ``k_seq`` splits the root into a MOBILE pool (state 0: the one that takes
up from the pore water and loads the xylem) and a SEQUESTERED pool (an extra
state APPENDED last: an irreversible apoplast/cell-wall sink that holds root
burden without feeding the shoot), i.e. the root mechanism of
``validation/twopool_root_exploration.py`` grafted onto this redistributed shoot.
That combination is what a per-organ out-of-sample test needs: the two-pool root
alone keeps the basic 4pool's PASS-THROUGH stem, so its stalk TF collapses and a
per-organ Tang comparison is unfair (``validation/twopool_root_oos_tang.py``,
Result 7); this module's shoot fixes exactly that.  ``k_seq=0`` (default) keeps
the single root pool and the state vector unchanged -- the pre-merge model is
recovered EXACTLY (no extra state is added), so every existing default,
validation number and test is untouched.  Reported root = mobile + sequestered
(``root_total``); ``k_rel`` is the slow seq->mobile desorption (0 = irreversible).
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
import pfas_rice_plant_module_4pool_surf as _P4   # SOLVER / SUBSTEPS (one switch)

from pfas_rice_plant_module_4pool_surf import (
    Environment, Compound, Compartment, binding_factors, dry_binding, root_uptake, grain_formed,
    panicle_share, MASS_FLOOR,
)


@dataclass
class PlantInputsNL:
    """Time-dependent drivers for the N-stem + leaf model.

    M : tissue fresh mass, shape (len(t), N+3) with columns ordered
        [root, stem_1..stem_N, leaf, grain].
    theta : optional organ water theta_fw(t), same shape as M, that made M fresh (PREREG 7.1);
        the binding factors then follow it (B = theta + (1 - theta) X). None: the compartments'
        constant theta.
    M_husk : optional DRY mass of the panicle remainder (husk + rachis; PREREG 7.2), len(t).
        Given, the model carries it as a terminal state that keeps the panicle's transpiration
        deposit, and the grain (brown rice) is fed by the phloem only.
    """
    t: np.ndarray
    Cwo: np.ndarray
    Qtp: np.ndarray
    M: np.ndarray
    leaf_loss: np.ndarray | None = None   # leaf senescence loss RATE [1/day] (opt; 0 if None)
    dvs: np.ndarray | None = None         # crop development stage (ORYZA DVS), for the
                                          # stage-based panicle share
    theta: np.ndarray | None = None       # organ water theta_fw(t), shape of M (PREREG 7.1)
    M_husk: np.ndarray | None = None      # panicle remainder DRY mass (PREREG 7.2)

    def __post_init__(self):
        self.M = np.asarray(self.M, dtype=float)
        self.n_comp = self.M.shape[1]
        kw = dict(kind="linear", bounds_error=False, fill_value="extrapolate")
        self._Cwo = interp1d(self.t, self.Cwo, **kw)
        self._Qtp = interp1d(self.t, self.Qtp, **kw)
        # ONE vectorised interpolator per matrix (axis=0) rather than one per column:
        # the RHS is called ~3k times per solve, so 2 scipy calls instead of 2*(N+3)
        # cuts the solve time ~3x. Same linear interpolant, same extrapolation.
        self._M = interp1d(self.t, self.M, axis=0, **kw)
        dM = np.gradient(self.M, self.t, axis=0)
        self._dM = interp1d(self.t, dM, axis=0, **kw)
        ll = np.zeros(len(self.t)) if self.leaf_loss is None else np.asarray(self.leaf_loss, float)
        self._leaf_loss = interp1d(self.t, ll, **kw)
        # grain FORMATION gate gamma(t): 0 while the grain (last column) sits on the
        # organ-mass floor (panicle not formed), 1 from the day it leaves it; 1 throughout
        # for a driver whose grain never touches the floor. Same rule as the 4-pool core
        # (`grain_formed`); a step because the old ramp was not converged.
        gate = grain_formed(self.M[:, -1])
        self._grain_gate = interp1d(self.t, gate, **kw)
        self._dvs = (None if self.dvs is None
                     else interp1d(self.t, np.asarray(self.dvs, float), **kw))
        # one table of every driver column for the RHS: one linear interpolation per call
        # instead of seven scipy calls (identical values: same linear interpolation, same
        # linear extrapolation beyond the ends)
        n = self.n_comp
        dv = (np.full(len(self.t), np.nan) if self.dvs is None
              else np.asarray(self.dvs, float))
        self._tg = np.asarray(self.t, dtype=float)
        self._tab = np.column_stack([np.asarray(self.Cwo, float), np.asarray(self.Qtp, float),
                                     self.M, dM, ll, gate, dv])
        self._sl = (slice(2, 2 + n), slice(2 + n, 2 + 2 * n), 2 + 2 * n, 3 + 2 * n, 4 + 2 * n)
        # the second table (PREREG 7.1/7.2): theta(t) per compartment, the grain's DRY growth,
        # the husk's dry mass, growth and formation gate; NaN/zeros when not given
        th = (np.full(self.M.shape, np.nan) if self.theta is None
              else np.asarray(self.theta, float).reshape(self.M.shape))
        gdry = (self.M[:, -1] * (1.0 - th[:, -1]) if self.theta is not None
                else np.full(len(self.t), np.nan))
        mh = (np.full(len(self.t), MASS_FLOOR) if self.M_husk is None
              else np.maximum(np.asarray(self.M_husk, float), MASS_FLOOR))
        self._tab2 = np.column_stack([th, np.gradient(gdry, self.t), mh,
                                      np.gradient(mh, self.t), grain_formed(mh)])
        self._sl2 = (slice(0, n), n, n + 1, n + 2, n + 3)

    def drivers_at(self, t):
        """(Cwo, Qtp, M, dM, leaf_loss, grain_gate, dvs) at time t, by one linear lerp."""
        tg = self._tg
        i = int(np.searchsorted(tg, t, side="right")) - 1
        i = min(max(i, 0), len(tg) - 2)
        w = (t - tg[i]) / (tg[i + 1] - tg[i])
        row = self._tab[i] + w * (self._tab[i + 1] - self._tab[i])
        sM, sdM, ill, ig, idv = self._sl
        return (row[0], row[1], row[sM], row[sdM], max(row[ill], 0.0),
                min(1.0, max(0.0, row[ig])), row[idv])

    def extras_at(self, t):
        """(theta[n], grain dry growth, husk dry mass, husk growth, husk gate) at time t."""
        tg = self._tg
        i = int(np.searchsorted(tg, t, side="right")) - 1
        i = min(max(i, 0), len(tg) - 2)
        w = (t - tg[i]) / (tg[i + 1] - tg[i])
        row = self._tab2[i] + w * (self._tab2[i + 1] - self._tab2[i])
        sth, igd, imh, idh, igh = self._sl2
        return row[sth], row[igd], row[imh], row[idh], min(1.0, max(0.0, row[igh]))

    def Cwo_(self, t): return float(self._Cwo(t))
    def Qtp_(self, t): return float(self._Qtp(t))
    def M_(self, t):   return np.asarray(self._M(t), dtype=float)
    def dM_(self, t):  return np.asarray(self._dM(t), dtype=float)
    def leaf_loss_(self, t):  return max(float(self._leaf_loss(t)), 0.0)
    def grain_gate_(self, t):  return min(1.0, max(0.0, float(self._grain_gate(t))))
    def dvs_(self, t):
        if self._dvs is None:
            raise ValueError("this driver set carries no 'dvs'; the stage-based panicle "
                             "share needs the crop's development stage")
        return float(self._dvs(t))



@dataclass
class NStemLeafModel:
    """Redistributed-shoot uptake model (transpiration deposition + retention).

    comps : [root, stem_1..stem_N, leaf, grain]
    tau   : per-stem-segment transpiration draw-off (fraction of Q), len N
    lam_leaf, lam_grain : leaf / panicle transpiration draw-off (fraction of Q)
        Constraint: sum(tau) + lam_leaf + lam_grain == 1.
    retention : fraction of each organ's transpiration-deposited solute retained
        (terminal); 1-retention is carried on to the grain as residual xylem. 1 (the default)
        is the module's own physics: solute delivered to an evaporating organ is at a dead
        end; export from the leaf is the phloem term (docs/PREREG_sourced_values.md 5.2).
    lam_grain_dvs : (panicle share at heading, at maturity). When given, the panicle's share
        of transpiration follows the crop's development stage (`panicle_share`: 0 before
        heading) and the leaf takes 1 - sum(tau) - share; `lam_grain`/`lam_leaf` are then
        unused. Needs `inputs.dvs`.
    k_seq : irreversible root sequestration rate [1/day]. 0 (default) = ONE root
        pool (state vector unchanged, pre-merge model recovered exactly); >0 adds
        a sequestered root state (index `SEQ`, appended last) fed at k_seq*C_root.
    k_rel : slow seq->mobile desorption rate [1/day] (0 = irreversible sink).
    """
    env: Environment
    cmpd: Compound
    comps: list
    inputs: PlantInputsNL
    tau: np.ndarray
    lam_leaf: float | None = None     # static split (when lam_grain_dvs is not given)
    lam_grain: float | None = None
    retention: float = 1.0
    phi: float | None = None          # phloem stream to the root / transpiration: no default
    T_C_Ph: float | None = None       # phloem sap per kg grain DM [L/kg]: no default (model_api)
    k_seq: float = 0.0                # root mobile -> sequestered rate [1/day]
    k_rel: float = 0.0                # sequestered -> mobile desorption [1/day]
    lam_grain_dvs: tuple | None = None   # (share at heading, at maturity); see the docstring
    husk_phloem: bool = False         # the husk's dry growth draws phloem sap too (a band end,
                                      # PREREG 7.2); False: only the leaf loads, only brown rice draws

    def __post_init__(self):
        if self.phi is None or self.T_C_Ph is None:
            raise ValueError("phi and T_C_Ph have no default: pass the sourced values "
                             "(model_api.PHI, model_api.T_C_PH)")
        self.N = len(self.comps) - 3          # number of stem segments
        self.ROOT = 0
        self.LEAF = self.N + 1
        self.GRAIN = self.N + 2
        self.n_comp = len(self.comps)
        # two-pool root: the sequestered state is APPENDED so the organ indices
        # (ROOT/LEAF/GRAIN and the stem slice) are identical in both modes.
        self.two_pool = bool(self.k_seq > 0.0 or self.k_rel > 0.0)
        self.SEQ = self.n_comp if self.two_pool else None
        self.n_states = self.n_comp + (1 if self.two_pool else 0)
        # panicle remainder (husk + rachis) as a terminal state, appended last (PREREG 7.2)
        self.husk = self.inputs.M_husk is not None
        self.HUSK = self.n_states if self.husk else None
        self.n_states += 1 if self.husk else 0
        if self.husk_phloem and not self.husk:
            raise ValueError("husk_phloem needs the husk pool (inputs.M_husk)")
        assert self.k_seq >= 0.0 and self.k_rel >= 0.0
        self.tau = np.asarray(self.tau, dtype=float)
        assert len(self.tau) == self.N, "tau must have one entry per stem segment"
        self._tau_sum = float(self.tau.sum())
        if self.lam_grain_dvs is None:
            if self.lam_leaf is None or self.lam_grain is None:
                raise ValueError("give lam_grain_dvs, or both lam_leaf and lam_grain")
            tot = self._tau_sum + self.lam_leaf + self.lam_grain
            assert abs(tot - 1.0) < 1e-9, (
                f"transpiration split must sum to 1 (sum(tau)+lam_leaf+lam_grain={tot:.6f})")
            assert self.lam_leaf >= 0 and self.lam_grain >= 0
        else:
            if self.inputs.dvs is None:
                raise ValueError("lam_grain_dvs needs the drivers' 'dvs'")
            assert 1.0 - self._tau_sum - max(self.lam_grain_dvs) >= 0.0, \
                "stem + panicle shares exceed the transpiration stream"
        assert self.tau.min() >= 0
        assert 0.0 <= self.retention <= 1.0
        # binding + metabolism depend only on (comps, cmpd), both fixed at
        # construction -- hoist them out of the RHS (called ~3k times per solve).
        self._B = binding_factors(self.comps, self.cmpd)
        self._X = dry_binding(self.comps, self.cmpd)
        self._staged = self.inputs.theta is not None       # B follows theta(t) (PREREG 7.1)
        self._gamma = np.array([c.gamma for c in self.comps])

    def rhs(self, t, C):
        Cwo, Q, M, dM, leaf_loss, gam, dvs = self.inputs.drivers_at(t)
        th, gdry, Mh, dMh, gh = self.inputs.extras_at(t)
        M = np.maximum(M, MASS_FLOOR)
        mu = dM / M
        B = th + (1.0 - th) * self._X if self._staged else self._B
        Cw = C[:self.n_comp] / B
        g = self._gamma
        dC = np.zeros_like(C)

        leaf, grain, r = self.LEAF, self.GRAIN, self.retention
        # phloem: loaded at the leaf (free + optional lipid-bound); two streams -- to the grain
        # the sap its DRY growth draws (T_C_Ph per kg DM; M is fresh), to the root phi x Q_TP
        if not self._staged:
            gdry = dM[grain] * (1.0 - self.comps[grain].theta)
        Q_g = max(gdry * self.T_C_Ph, 0.0)
        Q_r = max(self.phi * Q, 0.0)
        Q_h = max(dMh * self.T_C_Ph, 0.0) * gh if self.husk_phloem else 0.0
        C_Phl = self.cmpd.L_Ph * Cw[leaf] + self.cmpd.g_ph * C[leaf]

        # --- root: uptake, full xylem export (f_xy free + lipid g_xy bound), phloem return ---
        jR = root_uptake(Cwo, Cw[self.ROOT], self.cmpd, self.env)
        Cw_xyl = self.cmpd.f_xy * Cw[self.ROOT] + self.cmpd.g_xy * C[self.ROOT]
        xyl_flux = Q * Cw_xyl                          # total solute loaded into the xylem
        dC[self.ROOT] = (jR
                         - xyl_flux / M[self.ROOT]
                         + (Q_r / M[self.ROOT]) * C_Phl
                         - g[self.ROOT] * C[self.ROOT] - mu[self.ROOT] * C[self.ROOT])

        # --- stem segments: each retains its transpiration-deposited solute ---
        for s in range(1, self.N + 1):
            dC[s] = (r * self.tau[s - 1] * xyl_flux / M[s]
                     - g[s] * C[s] - mu[s] * C[s])

        # transpiration split at this stage: the panicle's sourced share (0 before heading)
        if self.lam_grain_dvs is not None:
            lg = panicle_share(dvs, *self.lam_grain_dvs)
            ll = 1.0 - self._tau_sum - lg
        else:
            lg, ll = self.lam_grain, self.lam_leaf
        # grain FORMATION gate: before the panicle sets (gam~0) the grain takes NO xylem/
        # phloem -- its residual xylem is rerouted to the leaf and its phloem stream is not
        # exported, so the balance closes (DPU-consistent).
        residual = (lg + (1.0 - r) * (self._tau_sum + ll)) * xyl_flux
        # the panicle's xylem goes to the grain, or, with the husk pool, to the husk: no xylem
        # flow reaches the intact caryopsis (Zee 1972; Oparka & Gates 1981; PREREG 7.2)
        # brown rice (the husk pool on) draws phloem in proportion to its own dry growth, so an
        # organ not yet growing draws nothing: its phloem stream is not gated (the mass-floor
        # gate opened a day early or late on DVS rounding at heading; PREREG section 8)
        wx = gh if self.husk else gam
        gq = 1.0 if self.husk else gam

        # --- leaf: own deposit (+ the panicle's rerouted residual pre-formation); phloem
        #     source (gated grain share + root recirc); senescence loss carries PFAS away ---
        dC[leaf] = (r * ll * xyl_flux / M[leaf]
                    + (1.0 - wx) * residual / M[leaf]
                    - ((gq * Q_g + Q_r + Q_h) / M[leaf]) * C_Phl
                    - g[leaf] * C[leaf] - mu[leaf] * C[leaf]
                    - leaf_loss * C[leaf])

        # --- grain (formation-gated): panicle deposit + non-retained residual xylem + phloem;
        #     brown rice with the husk pool: the phloem only ---
        dC[grain] = ((0.0 if self.husk else wx * residual / M[grain])
                     + gq * (Q_g / M[grain]) * C_Phl
                     - g[grain] * C[grain] - mu[grain] * C[grain])
        if self.husk:                    # terminal, per kg DRY (no binding: nothing leaves it)
            dC[self.HUSK] = ((wx * residual + Q_h * C_Phl) / Mh
                             - (dMh / Mh) * C[self.HUSK])

        # --- optional SEQUESTERED root pool: an irreversible (k_rel=0) terminal
        #     sink fed from the mobile root, so the root can hold burden without
        #     draining the xylem feed. Both pools share the root mass M[ROOT], so
        #     the seq pool dilutes with root growth exactly like the mobile one.
        if self.two_pool:
            seq = self.k_seq * C[self.ROOT]
            rel = self.k_rel * C[self.SEQ]
            dC[self.ROOT] += rel - seq
            dC[self.SEQ] = seq - rel - mu[self.ROOT] * C[self.SEQ]
        return dC

    def root_total(self, C):
        """Reported root concentration = mobile + sequestered (axis-0 indexing, so
        it works for a state vector or a full (n_states, n_t) solution array)."""
        C = np.asarray(C)
        return C[self.ROOT] + (C[self.SEQ] if self.two_pool else 0.0)

    def affine(self, ts):
        """(A[m,n,n], b[m,n]) with rhs(t, C) = A(t) C + b(t) at the times `ts`, or None when
        the root entry is not linear (a carrier). The same terms as `rhs`, vectorised;
        `tests/test_linear_expm.py` holds them equal to the rhs."""
        import linear_expm as LE
        kio = LE.root_coefficients(root_uptake, self.cmpd, self.env)
        if kio is None:
            return None
        k_in, k_out = kio
        inp = self.inputs
        row = LE.lerp_table(inp._tg, inp._tab, ts)
        row2 = LE.lerp_table(inp._tg, inp._tab2, ts)
        sM, sdM, ill, ig, idv = inp._sl
        sth, igd, imh, idh, igh = inp._sl2
        Cwo, Q = row[:, 0], row[:, 1]
        M = np.maximum(row[:, sM], MASS_FLOOR)
        dM = row[:, sdM]
        lls = np.maximum(row[:, ill], 0.0)
        gam = np.clip(row[:, ig], 0.0, 1.0)
        dvs = row[:, idv]
        mu = dM / M
        g = self._gamma
        if self._staged:                                     # (m, n): theta(t) per compartment
            th = row2[:, sth]
            B = th + (1.0 - th) * self._X
        else:
            B = np.broadcast_to(self._B, (len(row), self.n_comp))
        R0, Lf, Gr, r = self.ROOT, self.LEAF, self.GRAIN, self.retention
        gdry = (row2[:, igd] if self._staged else dM[:, Gr] * (1.0 - self.comps[Gr].theta))
        Q_g = np.maximum(gdry * self.T_C_Ph, 0.0)
        Q_r = np.maximum(self.phi * Q, 0.0)
        Mh, dMh = row2[:, imh], row2[:, idh]
        gh = np.clip(row2[:, igh], 0.0, 1.0)
        Q_h = (np.maximum(dMh * self.T_C_Ph, 0.0) * gh if self.husk_phloem
               else np.zeros(len(row)))
        p = self.cmpd.L_Ph / B[:, Lf] + self.cmpd.g_ph
        x = self.cmpd.f_xy / B[:, R0] + self.cmpd.g_xy
        if self.lam_grain_dvs is not None:
            lg = np.array([panicle_share(d, *self.lam_grain_dvs) for d in dvs])
            lf = 1.0 - self._tau_sum - lg
        else:
            lg = np.full(len(row), float(self.lam_grain))
            lf = np.full(len(row), float(self.lam_leaf))
        Rc = (lg + (1.0 - r) * (self._tau_sum + lf)) * Q * x     # residual = Rc * C_root
        n = self.n_states
        A = np.zeros((len(row), n, n))
        b = np.zeros((len(row), n))
        wx = gh if self.husk else gam
        gq = np.ones_like(gam) if self.husk else gam   # brown rice: phloem not gated (rhs)
        A[:, R0, R0] = -k_out / B[:, R0] - Q * x / M[:, R0] - g[R0] - mu[:, R0]
        A[:, R0, Lf] = Q_r * p / M[:, R0]
        b[:, R0] = k_in * Cwo
        for s_ in range(1, self.N + 1):
            A[:, s_, R0] = r * self.tau[s_ - 1] * Q * x / M[:, s_]
            A[:, s_, s_] = -g[s_] - mu[:, s_]
        A[:, Lf, R0] = r * lf * Q * x / M[:, Lf] + (1.0 - wx) * Rc / M[:, Lf]
        A[:, Lf, Lf] = -(gq * Q_g + Q_r + Q_h) * p / M[:, Lf] - g[Lf] - mu[:, Lf] - lls
        A[:, Gr, R0] = 0.0 if self.husk else gam * Rc / M[:, Gr]
        A[:, Gr, Lf] = gq * Q_g * p / M[:, Gr]
        A[:, Gr, Gr] = -g[Gr] - mu[:, Gr]
        if self.husk:
            Hk = self.HUSK
            A[:, Hk, R0] = wx * Rc / Mh
            A[:, Hk, Lf] = Q_h * p / Mh
            A[:, Hk, Hk] = -dMh / Mh
        if self.two_pool:
            Sq = self.SEQ
            A[:, R0, R0] -= self.k_seq
            A[:, R0, Sq] += self.k_rel
            A[:, Sq, R0] = self.k_seq
            A[:, Sq, Sq] = -self.k_rel - mu[:, R0]
        return A, b

    def solve(self, t_eval, C0=None):
        import linear_expm as LE
        if C0 is None:
            C0 = np.zeros(self.n_states)
        if _P4.SOLVER == "expm" and self.affine(np.asarray(t_eval[:1], float)) is not None:
            return LE.solve_burden(self.affine, self._burden_scale, t_eval, self.n_states, C0,
                                   substeps=_P4.SUBSTEPS)
        D0, _ = self._burden_scale(np.asarray(t_eval[:1], float))
        sol = solve_ivp(LE.burden_rhs(self.rhs, self._burden_scale),
                        (float(t_eval[0]), float(t_eval[-1])), C0 * D0[0],
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
        sM, sdM = inp._sl[0], inp._sl[1]
        Gr = self.GRAIN
        D, dD = np.ones((len(ts), self.n_states)), np.zeros((len(ts), self.n_states))
        D[:, Gr] = np.maximum(row[:, sM][:, Gr], MASS_FLOOR)
        dD[:, Gr] = row[:, sdM][:, Gr]
        if self.husk:
            row2 = LE.lerp_table(inp._tg, inp._tab2, ts)
            _, _, imh, idh, _ = inp._sl2
            D[:, self.HUSK], dD[:, self.HUSK] = row2[:, imh], row2[:, idh]
        return D, dD

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
    """Per-segment stem draw-off: ``stem_transp_frac`` (the sourced sheath + culm share of
    transpiration) split evenly over the N segments.

    No source gives how the stalk's transpiration divides with height, and nothing reported
    depends on it: with retention 1 a segment exports nothing, so its burden is its own deposit
    and the stem (the mass-weighted mean over equal segments), the straw, the root and the grain
    depend only on the sum. (A height gradient 1.3 -> 0.7 stood here without a source.)"""
    return np.full(int(N), float(stem_transp_frac) / int(N))


def split_from_stem_frac(N, stem_transp_frac, lam_grain):
    """Convenience: (tau, lam_leaf, lam_grain) from a single stalk-transpiration
    fraction, with the leaf taking the remainder.  Closes to sum 1."""
    tau = default_tau(N, stem_transp_frac)
    lam_leaf = 1.0 - stem_transp_frac - lam_grain
    return tau, lam_leaf, lam_grain
