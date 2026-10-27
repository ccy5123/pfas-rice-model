# =============================================================================
# linear_expm.py -- exponential stepping for the plant ODE (it is linear in the state)
# -----------------------------------------------------------------------------
# Every plant model here is affine in its state: dC/dt = A(t) C + b(t). The drivers enter A
# and b; the state never does (root entry, xylem loading, phloem, growth dilution and the
# senescence loss are all first order in C). So the solution over a step is exact for a frozen
# A and a linear b:
#
#     [C(t+h)]          ( [A  b0  b1]     )  [C(t)]
#     [  1   ] = expm(  [0   0   1] * h ) [ 0  ]  ... with the augmented 2 rows
#     [  0   ]          ( [0   0   0]     )  [ 1  ]
#
# A and b are read off the model's OWN right-hand side (A e_j = rhs(t, e_j) - rhs(t, 0),
# b = rhs(t, 0)), so the physics stays in one place. Each driver interval is split into
# `substeps` steps with A frozen at the step midpoint (second order in the variation of A);
# b runs linearly between its values at the step ends. Stiff modes (the root's equilibrium
# rate lam_eq) are integrated exactly, so no step-size control is needed.
#
# Dense output (the harvest operator's sub-daily dead-leaf integral) is a cubic Hermite
# interpolant through the step states and their derivatives.
# =============================================================================
from __future__ import annotations

import numpy as np
from scipy.linalg import expm


# Pade(13) scaling-and-squaring (Higham 2005, the algorithm scipy.linalg.expm uses), batched
# with numpy. scipy's batched expm runs its small products through multithreaded BLAS, which
# is ~40x slower for these 10x10 stacks when several worker processes share the CPUs.
_B13 = (64764752532480000., 32382376266240000., 7771770303897600., 1187353796428800.,
        129060195264000., 10559470521600., 670442572800., 33522128640., 1323241920.,
        40840800., 960960., 16380., 182., 1.)
_THETA13 = 5.371920351148152


def expm_batch(A):
    """exp of each matrix in the stack A[m, n, n] (Pade 13 with scaling and squaring)."""
    A = np.asarray(A, float)
    norm = np.abs(A).sum(axis=1).max(axis=1)            # 1-norm of each matrix
    s = np.maximum(0, np.ceil(np.log2(np.maximum(norm, 1e-300) / _THETA13))).astype(int)
    As = A / (2.0 ** s)[:, None, None]
    n = A.shape[1]
    I = np.broadcast_to(np.eye(n), A.shape)
    A2 = As @ As
    A4 = A2 @ A2
    A6 = A4 @ A2
    b = _B13
    U = As @ (A6 @ (b[13] * A6 + b[11] * A4 + b[9] * A2)
              + b[7] * A6 + b[5] * A4 + b[3] * A2 + b[1] * I)
    V = (A6 @ (b[12] * A6 + b[10] * A4 + b[8] * A2)
         + b[6] * A6 + b[4] * A4 + b[2] * A2 + b[0] * I)
    E = np.linalg.solve(V - U, V + U)
    for k in range(int(s.max()) if len(s) else 0):
        sel = s > k
        E[sel] = E[sel] @ E[sel]
    return E


class _Dense:
    def __init__(self, t, y, f):
        self.t, self.y, self.f = t, y, f

    def __call__(self, tt):
        tt = np.atleast_1d(np.asarray(tt, float))
        t, y, f = self.t, self.y, self.f
        i = np.clip(np.searchsorted(t, tt, side="right") - 1, 0, len(t) - 2)
        h = t[i + 1] - t[i]
        s = (tt - t[i]) / h
        h00 = 2 * s ** 3 - 3 * s ** 2 + 1
        h10 = s ** 3 - 2 * s ** 2 + s
        h01 = -2 * s ** 3 + 3 * s ** 2
        h11 = s ** 3 - s ** 2
        out = (y[:, i] * h00 + f[:, i] * (h10 * h) + y[:, i + 1] * h01
               + f[:, i + 1] * (h11 * h))
        return out if out.shape[1] > 1 else out[:, 0]


class LinearSolution:
    """The fields of a scipy OdeResult that the callers read: t, y, success, sol."""
    def __init__(self, t, y, sol):
        self.t, self.y, self.sol = t, y, sol
        self.success = bool(np.all(np.isfinite(y)))
        self.status = 0 if self.success else -1
        self.message = "exponential stepping (linear system)"


def affine_parts(rhs, t, n):
    """(A, b) of rhs(t, C) = A C + b, read off the right-hand side itself."""
    z = np.zeros(n)
    b = np.asarray(rhs(t, z), float)
    A = np.empty((n, n))
    e = np.zeros(n)
    for j in range(n):
        e[j] = 1.0
        A[:, j] = np.asarray(rhs(t, e), float) - b
        e[j] = 0.0
    return A, b


def lerp_table(tg, tab, ts, hold=None):
    """Row-wise linear interpolation of a driver table, exactly as `drivers_at` does it
    (the same interval choice and linear extrapolation beyond the ends), at many times.
    Columns flagged in `hold` keep the interval's left-node value (the slopes of
    `interval_slopes`)."""
    ts = np.asarray(ts, float)
    i = np.clip(np.searchsorted(tg, ts, side="right") - 1, 0, len(tg) - 2)
    w = (ts - tg[i]) / (tg[i + 1] - tg[i])
    out = tab[i] + w[:, None] * (tab[i + 1] - tab[i])
    if hold is not None:
        out[:, hold] = tab[i][:, hold]
    return out


def lerp_row(tg, tab, t, hold=None):
    """One row of `lerp_table` at a scalar time (the rhs's driver read)."""
    i = int(np.searchsorted(tg, t, side="right")) - 1
    i = min(max(i, 0), len(tg) - 2)
    w = (t - tg[i]) / (tg[i + 1] - tg[i])
    row = tab[i] + w * (tab[i + 1] - tab[i])
    if hold is not None:
        row[hold] = tab[i][hold]
    return row


def interval_slopes(tg, y):
    """The slope of y's piecewise-linear interpolant on each driver interval, stored at the
    interval's left node (the last node repeats the last slope): the growth that matches the
    linearly interpolated mass inside the interval. Used for the growth that draws an inflow
    into a terminal pool (the grain's phloem, the husk): a central difference assigns half of
    the first day's growth to the day before, when the organ still sits on the mass floor, and
    the pool's concentration there is burden / floor (PREREG section 8)."""
    y = np.asarray(y, float)
    dt = np.diff(np.asarray(tg, float)).reshape((-1,) + (1,) * (y.ndim - 1))
    s = np.diff(y, axis=0) / dt
    return np.concatenate([s, s[-1:]], axis=0)


def root_coefficients(root_uptake, cmpd, env):
    """(k_in, k_out) of a LINEAR root entry j = k_in*Cwo - k_out*Cw_root, read off the
    model's own `root_uptake`; None when the entry is not linear (a carrier)."""
    if cmpd.Vmax_in != 0.0 or cmpd.Vmax_out != 0.0:
        return None
    k_in = float(root_uptake(1.0, 0.0, cmpd, env))
    k_out = -float(root_uptake(0.0, 1.0, cmpd, env))
    ok = (abs(root_uptake(3.0, 2.0, cmpd, env) - (3.0 * k_in - 2.0 * k_out))
          <= 1e-12 * max(1.0, abs(k_in), abs(k_out)))
    return (k_in, k_out) if ok else None


# graded steps: the midpoint rule's error grows with how much A changes across a step, so an
# interval where a rate moves a lot over its length (an organ growing from a small mass, the
# grain gate opening) gets more steps, in proportion to the change of the exponent per GRADE,
# up to GRADE_MAX x the base count. For a slow mode (|a| H < 1) that change is |delta a| H; a
# fast one follows its quasi-equilibrium b/a, so only its RELATIVE change |delta a|/|a| counts.
GRADE = 0.2
GRADE_MAX = 64


def step_counts(affine, t_eval, k):
    """Exponential steps per driver interval: `k`, raised where A changes fast."""
    A, _ = affine(t_eval)
    d = np.diagonal(A, axis1=1, axis2=2)
    H = np.diff(t_eval)[:, None]
    dbar = 0.5 * np.abs(d[1:] + d[:-1])
    r = (np.abs(np.diff(d, axis=0)) * H / np.maximum(1.0, dbar * H)).max(axis=1)
    return k * np.clip(np.ceil(r / GRADE), 1, GRADE_MAX).astype(int)


def solve_affine(affine, t_eval, n, C0=None, substeps=8):
    """Integrate dC/dt = A(t) C + b(t) with `affine(ts) -> (A[m,n,n], b[m,n])` vectorised:
    each driver interval in `substeps` steps (more where the rates change fast,
    `step_counts`), A frozen at the step midpoint, b linear over the step (exact for that
    frozen A); dense output by cubic Hermite through the steps."""
    t_eval = np.asarray(t_eval, float)
    ks = step_counts(affine, t_eval, int(substeps))
    tg = np.concatenate([np.linspace(a, b, ki + 1)[:-1]
                         for a, b, ki in zip(t_eval[:-1], t_eval[1:], ks)] + [t_eval[-1:]])
    h = np.diff(tg)
    Am, _ = affine(0.5 * (tg[:-1] + tg[1:]))
    An, bn = affine(tg)
    m = len(h)
    G = np.zeros((m, n + 2, n + 2))
    G[:, :n, :n] = Am
    G[:, :n, n] = bn[:-1]
    G[:, :n, n + 1] = (bn[1:] - bn[:-1]) / h[:, None]
    G[:, n, n + 1] = 1.0
    E = expm_batch(G * h[:, None, None])
    Phi, q = E[:, :n, :n], E[:, :n, n]
    ys = np.empty((n, len(tg)))
    C = np.zeros(n) if C0 is None else np.asarray(C0, float).copy()
    ys[:, 0] = C
    for i in range(m):
        C = Phi[i] @ C + q[i]
        ys[:, i + 1] = C
    fs = (np.einsum("mij,jm->im", An, ys) + bn.T)
    idx = np.concatenate([[0], np.cumsum(ks)])
    return LinearSolution(t_eval, ys[:, idx], _Dense(tg, ys, fs))


# --- terminal pools carried as burden ---------------------------------------------------
# A terminal pool (the grain, the husk) has no outflow but growth dilution, so in the
# concentration form its rate dM/M diverges where the organ starts from the mass floor, and the
# relaxation time M_floor/dM there is far below any step a solver can take. In burden
# coordinates y_k = M_k C_k the dilution term cancels exactly and the pool's equation is its
# inflow alone (PREREG section 8). `scale(ts) -> (D, dD)` gives D = M (dD = dM/dt) on those
# rows and 1 (0) elsewhere; rhs and affine stay in the concentration form.
def burden_affine(affine, scale):
    """The affine system of y = D(t) C: A_y = D A D^-1 + diag(dD/D), b_y = D b."""
    def aff(ts):
        A, b = affine(ts)
        D, dD = scale(ts)
        Ay = D[:, :, None] * A / D[:, None, :]
        k = np.arange(A.shape[1])
        Ay[:, k, k] += dD / D
        return Ay, D * b
    return aff


def burden_rhs(rhs, scale):
    """dy/dt of y = D(t) C from the concentration-form rhs."""
    def f(t, y):
        D, dD = scale(np.array([t], float))
        C = y / D[0]
        return D[0] * np.asarray(rhs(t, C), float) + dD[0] * C
    return f


class _BurdenDense:
    def __init__(self, sol, scale):
        self._sol, self._scale = sol, scale

    def __call__(self, tt):
        tt = np.atleast_1d(np.asarray(tt, float))
        D, _ = self._scale(tt)
        y = np.asarray(self._sol(tt), float).reshape(D.shape[1], len(tt)) / D.T
        return y if y.shape[1] > 1 else y[:, 0]          # as the wrapped dense output


class BurdenSolution:
    """A solution in burden coordinates, read back as concentrations (t, y, success, sol)."""
    def __init__(self, sol, scale):
        D, _ = scale(np.asarray(sol.t, float))
        self.t, self.y = sol.t, np.asarray(sol.y, float) / D.T
        self.success, self.status = sol.success, getattr(sol, "status", 0)
        self.message = sol.message
        self.sol = (_BurdenDense(sol.sol, scale) if getattr(sol, "sol", None) is not None
                    else None)


def solve_burden(affine, scale, t_eval, n, C0=None, substeps=8):
    """`solve_affine` with the terminal pools as burden (`scale`), returned as concentrations."""
    t_eval = np.asarray(t_eval, float)
    D0, _ = scale(t_eval[:1])
    y0 = None if C0 is None else np.asarray(C0, float) * D0[0]
    return BurdenSolution(solve_affine(burden_affine(affine, scale), t_eval, n, y0, substeps),
                          scale)


def solve_linear(rhs, t_eval, n, C0=None, substeps=1):
    """Integrate dC/dt = rhs(t, C) (affine in C) on the grid `t_eval`."""
    t_eval = np.asarray(t_eval, float)
    C = np.zeros(n) if C0 is None else np.asarray(C0, float).copy()
    # the step grid: every driver interval split into `substeps`
    k = int(substeps)
    tg = np.concatenate([np.linspace(a, b, k + 1)[:-1] for a, b in zip(t_eval[:-1], t_eval[1:])]
                        + [t_eval[-1:]])
    ys = np.empty((n, len(tg)))
    fs = np.empty((n, len(tg)))
    ys[:, 0] = C
    _, b_prev = affine_parts(rhs, tg[0], n)   # b at the step's left end
    G = np.zeros((n + 2, n + 2))
    G[n, n + 1] = 1.0
    for i in range(len(tg) - 1):
        t0, t1 = tg[i], tg[i + 1]
        h = t1 - t0
        Am, _ = affine_parts(rhs, 0.5 * (t0 + t1), n)
        b1 = np.asarray(rhs(t1, np.zeros(n)), float)
        G[:n, :n] = Am
        G[:n, n] = b_prev                       # b(t0)
        G[:n, n + 1] = (b1 - b_prev) / h        # slope of b over the step
        E = expm(G * h)
        C = E[:n, :n] @ C + E[:n, n]
        ys[:, i + 1] = C
        b_prev = b1
    for i in range(len(tg)):
        fs[:, i] = np.asarray(rhs(tg[i], ys[:, i]), float)
    idx = np.arange(0, len(tg), k)
    return LinearSolution(t_eval, ys[:, idx], _Dense(tg, ys, fs))
