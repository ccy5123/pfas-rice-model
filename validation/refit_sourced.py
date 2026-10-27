#!/usr/bin/env python3
# =============================================================================
# validation/refit_sourced.py
# -----------------------------------------------------------------------------
# RE-FIT the PFAS transport parameters after the no-unsourced-value amendment
# (docs/PREREG_sourced_values.md, amendment 1 with addendum 5.6 -- committed before this
# script was run). Inputs at their sourced midpoints; observations from the SI
# (data_obs/obs_baf_Yamazaki.csv, built by docs/literature_db/tools/yamazaki_obs_from_si.py).
#
# PROTOCOL (PREREG 5.3 and 5.6)
#   Variants   redistributed_lipid (default) | redistributed_free (monotone f_xy_recommended)
#              | redistributed_fit (f_xy fitted to straw; a reproduction) | 4pool_lipid |
#              4pool_free. Model output: the dry-basis harvest BAF of the two harvest runs
#              (Kasalath 09-01, Koshihikari 09-10), combined by geometric mean as the data are.
#   r_eq       closed form: the model is linear in the state (zero initial state, constant
#              exposure), so every tissue scales with r_eq; checked numerically (self_checks).
#   L_Ph >= 0  solves the grain/root ratio (log10 tolerance 1e-4), no upper bound; L_Ph = 0
#              exactly when the model's ratio at 0 already reaches the observed one; a ratio no
#              L_Ph reaches is recorded ("unreachable", L_Ph at the ratio's maximum).
#   *_fit      f_xy in (0, 1] solves the straw/root ratio, alternating with L_Ph.
#   PFBS       no grain value: L_Ph in [0, L_MQL] (L_MQL: model grain = the censoring bound, the
#              observation's `hi`), run at the midpoint; the range goes to the input-range band.
#   Intervals  per congener: each cultivar fitted on its own observation and harvest run (the
#              bootstrap over the two cultivars).
#   Lipid      per shoot: profiled least squares (r_eq, L_Ph re-solved per evaluation) on the
#              residuals the per-congener fits leave -- every congener's straw and the grain of a
#              congener at L_Ph = 0; physical bounds only (0 < a <= 1, g >= 0, K_half > 0; b and
#              pfsa_ln free); PFDoDA included; start = the previous fit's constants.
#   Elimination each term in turn at its bound (g_ph = 0; g_xy = 0; gate "none"; gate "linear";
#              pfsa_ln = 0; b = 0), the rest refitted; removed when decide(full, reduced) on the
#              profiled per-cell errors is not "worse" for the reduced model; the cheapest removal
#              first; until every remaining removal is "worse".
#
#   python validation/refit_sourced.py            # fits, prints, writes refit_sourced.json
#   python validation/refit_sourced.py --write    # ... and params/refit_sourced.csv +
#                                                 #     params/lipid_loading_sourced.json
# =============================================================================
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.optimize import brentq, least_squares, minimize_scalar

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "docs", "literature_db", "tools"))
import model_api as api                 # noqa: E402
import site_scenario as ss              # noqa: E402
import decision_stats as D              # noqa: E402
import yamazaki_obs_from_si as Y        # noqa: E402

TISS = ("root", "straw", "grain")
TOL = 1e-4                               # log10 residual tolerance of every 1-D solve
# a grain/root ratio that no L_Ph in [0, 1e8] reaches is "unreachable"; L_Ph is then stored at
# the maximum of the model's ratio (the closest approach; the ratio is not monotone in L_Ph --
# a large L_Ph drains the leaf back to the root through the phi stream)
L_TINY = 1e-12                           # numerically 0 for the log-space bracket
LOGG_MIN = -12.0                         # log10 g floor during the optimisation (g*C negligible)
HARVESTS = ss.YAMAZAKI_HARVESTS          # run i <-> cultivar i (Kasalath, Koshihikari)
OBS = api.yamazaki_obs()
CONG = [c for c in api.CONGENERS if c in OBS]
NO_GRAIN = [c for c in CONG if not OBS[c]["use"]["grain"]]
VARIANTS = {
    "redistributed_lipid": ("redistributed", "lipid"),
    "redistributed_free": ("redistributed", "free"),
    "redistributed_fit": ("redistributed", "fit"),
    "4pool_lipid": ("4pool", "lipid"),
    "4pool_free": ("4pool", "free"),
}
# the *_free variants are fitted for the record; they are withdrawn from scoring and decisions
# (model_api.WITHDRAWN_VARIANTS, PREREG 9.4: their f_xy_recommended is untraceable)
# starting point only: the constants of the previous fit (refit_site_year.py, on the 4-pool)
START = dict(a=0.1177, b=0.09053, g_xy=0.01577, g_ph=0.005005, K_half=800.1, pfsa_ln=0.6213)
OUT = os.path.join(HERE, "refit_sourced.json")


# ---------------------------------------------------------------- observations
def per_cultivar_obs():
    """{congener: [ {root, straw, grain} for Kasalath, Koshihikari ]} [L/kg dw], midpoints."""
    si = Y.read_si()
    sh = Y.straw_shares()
    out = {}
    for c in CONG:
        w = si[c]["water"]
        cvs = list(si[c]["cv"].values())
        out[c] = [dict(root=cv["root"][1] / w, straw=Y.straw(cv, sh)[1] / w,
                       grain=Y.brown(cv)[1] / w) for cv in cvs]
    return out


OBS_CV = per_cultivar_obs()


def obs_for(c, runs):
    """Observation matching a run set: both runs -> the table (cultivar geomean); one run ->
    that cultivar's own values."""
    if len(runs) == 2:
        return dict(OBS[c]["baf_dw"])
    return dict(OBS_CV[c][runs[0]])


# ---------------------------------------------------------------- model runs (worker side)
_DRV = None
_CACHE = {}


def _drivers():
    global _DRV
    if _DRV is None:
        _DRV = [api.scenario_drivers(Cwo=1.0, harvest=h) for h in HARVESTS]
    return _DRV


def model_logs(c, shoot, f_xy, g_xy, g_ph, L, runs=(0, 1), r_eq=1.0, Cwo=1.0):
    """log10 {root, straw, grain} dry-basis BAF, geometric mean over `runs` (r_eq given)."""
    key = (c, shoot, f_xy, g_xy, g_ph, L, runs, r_eq, Cwo)
    if key in _CACHE:
        return _CACHE[key]
    logs = []
    for i in runs:
        d = _drivers()[i]
        if Cwo != 1.0:
            d = dict(d, Cwo=d["Cwo"] * Cwo)
        r = api.simulate(c, drivers=d, shoot=shoot, lipid_loading=False, f_xy_override=f_xy,
                         g_xy_override=g_xy, g_ph_override=g_ph, L_Ph_override=L,
                         r_eq_override=r_eq)
        b = r["baf_final_dw"]
        logs.append([np.log10(max(b[k], 1e-300)) for k in TISS])
    v = dict(zip(TISS, np.mean(logs, axis=0)))
    _CACHE[key] = v
    return v


def solve_L(c, shoot, f_xy, g_xy, g_ph, target, runs, L0=None, unreach0=False):
    """L_Ph >= 0 with log10(grain/root) = target. Returns (L, status, logs). `L0` warm-starts
    the search; `unreach0` says L0 was the maximum of an unreachable ratio."""
    def g(L):
        m = model_logs(c, shoot, f_xy, g_xy, g_ph, L, runs)
        return m["grain"] - m["root"] - target, m

    def root_in(lo, hi):
        x = brentq(lambda u: g(10.0 ** u)[0], np.log10(lo), np.log10(hi), xtol=1e-7)
        L = float(10.0 ** x)
        r, m = g(L)
        assert abs(r) < TOL, (c, shoot, r)
        return L, "ok", m

    if unreach0 and L0 is not None and L0 > 0:          # warm: re-maximise near the old peak
        opt = minimize_scalar(lambda u: -g(10.0 ** u)[0],
                              bounds=(np.log10(L0) - 1, np.log10(L0) + 1),
                              method="bounded", options=dict(xatol=1e-3))
        Lm = float(10.0 ** opt.x)
        gm, mm = g(Lm)
        if gm < 0 and abs(opt.x - np.log10(L0)) < 0.99:  # still unreachable, peak inside
            return Lm, "unreachable", mm
    if L0 is not None and np.isfinite(L0) and L0 > 0:   # warm start: widening brackets
        for w in (0.02, 0.2, 1.0):
            lo, hi = L0 * 10.0 ** -w, L0 * 10.0 ** w
            if g(lo)[0] < 0 < g(hi)[0]:
                return root_in(lo, hi)
    g0, m0 = g(0.0)
    if g0 >= -TOL:
        return 0.0, ("bound0" if g0 > TOL else "ok"), m0
    if g(L_TINY)[0] >= 0:                                # rises between 0 and L_TINY
        return L_TINY, "ok", g(L_TINY)[1]
    scan = [(0.0, g0)]
    prev = L_TINY
    for k in range(-6, 9):                                # decades 1e-6 .. 1e8
        L = 10.0 ** k
        gv = g(L)[0]
        if gv >= 0:
            return root_in(prev, L)
        scan.append((L, gv))
        prev = L
    # unreachable: the grain/root ratio peaks below the observation -- keep the L_Ph of its
    # maximum (the closest approach), found within the best decade
    Lb, gb = max(scan, key=lambda z: z[1])
    if Lb == 0.0:
        return 0.0, "unreachable", m0
    opt = minimize_scalar(lambda u: -g(10.0 ** u)[0], bounds=(np.log10(Lb) - 1, np.log10(Lb) + 1),
                          method="bounded", options=dict(xatol=1e-3))
    L = float(10.0 ** opt.x)
    return L, "unreachable", g(L)[1]


def solve_joint(c, shoot, f0, L0, t_straw, t_grain, runs):
    """(f_xy, L_Ph) solving the straw and grain cells together (free anion, *_fit variant),
    started from the alternation's last point. Returns (f, L, logs, status)."""
    def r(x):
        m = model_logs(c, shoot, float(10.0 ** x[0]), 0.0, 0.0, float(10.0 ** x[1]), runs)
        return [m["straw"] - m["root"] - t_straw, m["grain"] - m["root"] - t_grain]
    x0 = [np.log10(min(max(f0, 1e-6), 1.0)), np.log10(max(L0, 1e-8))]
    sol = least_squares(r, x0, bounds=([-6.0, -8.0], [0.0, 8.0]), xtol=1e-10, ftol=1e-12,
                        gtol=1e-12, diff_step=1e-4)
    f, L = float(10.0 ** sol.x[0]), float(10.0 ** sol.x[1])
    m = model_logs(c, shoot, f, 0.0, 0.0, L, runs)
    ok = max(abs(v) for v in sol.fun) < TOL
    return f, L, m, ("ok" if ok else ("bound_f1" if sol.x[0] >= -1e-9 else "not_solved"))


def solve_fxy(c, shoot, L, target, runs):
    """f_xy in (0, 1] with log10(straw/root) = target (free anion). (f_xy, status)."""
    def g(f):
        m = model_logs(c, shoot, f, 0.0, 0.0, L, runs)
        return m["straw"] - m["root"] - target
    if g(1.0) <= TOL:
        return 1.0, ("bound1" if g(1.0) < -TOL else "ok")
    lo = 1e-6
    if g(lo) >= 0:
        return lo, "bound_low"
    x = brentq(lambda u: g(10.0 ** u), np.log10(lo), 0.0, xtol=1e-7)
    return float(10.0 ** x), "ok"


def conductances(c, loading, p=None, f_fit=None):
    rec = api._CONG[c]
    if loading == "explicit":                            # p = (f_xy, g_xy, g_ph)
        return tuple(float(x) for x in p)
    if loading == "lipid":
        return api.lipid_loading_conductances(rec["n_C"], rec["K_PL_Lkg"], rec["group"], p=p)
    if loading == "fit":
        return f_fit, 0.0, 0.0
    return float(rec["f_xy_recommended"]), 0.0, 0.0


def fit_congener(args):
    """Per-congener fit for one variant (and the lipid constants `p` if any).
    Returns dict(r_eq, L_Ph, f_xy, status, res (log10 residuals), logs)."""
    c, shoot, loading, p, runs, warm = args
    runs = tuple(runs)
    o = obs_for(c, runs)
    lo_root = np.log10(o["root"])
    has_grain = OBS[c]["use"]["grain"]
    t_straw = np.log10(o["straw"]) - lo_root
    f_fit = (warm or {}).get("f_xy", float(api._CONG[c]["f_xy_recommended"]))
    L0 = (warm or {}).get("L_Ph")
    unreach = (warm or {}).get("status", {}).get("L_Ph") == "unreachable"
    status = {}
    for _ in range(30):                                   # alternation (one pass unless "fit")
        f, gx, gp = conductances(c, loading, p, f_fit)
        if has_grain:
            t_grain = np.log10(o["grain"]) - lo_root
            L, status["L_Ph"], m = solve_L(c, shoot, f, gx, gp, t_grain, runs, L0, unreach)
            unreach = status["L_Ph"] == "unreachable"
        else:                                             # PFBS: [0, L_MQL], midpoint
            t_mql = np.log10(OBS[c]["hi"]["grain"]) - lo_root
            Lm, st, _ = solve_L(c, shoot, f, gx, gp, t_mql, runs, None if L0 is None else 2 * L0)
            L = 0.5 * Lm if np.isfinite(Lm) else float("inf")
            status["L_Ph"] = "mql_mid" if st == "ok" else st
            status["L_MQL"] = Lm
            m = model_logs(c, shoot, f, gx, gp, L, runs)
        if loading != "fit":
            break
        f_new, status["f_xy"] = solve_fxy(c, shoot, L, t_straw, runs)
        m = model_logs(c, shoot, f_new, 0.0, 0.0, L, runs)
        dL = (0.0 if (L0 is not None and L == L0) else
              abs(np.log10(L / L0)) if (L0 is not None and L > 0 and L0 > 0) else np.inf)
        done = abs(np.log10(f_new / f_fit)) < 1e-5 and dL < 1e-3
        f_fit, L0 = f_new, L
        if done:
            break
    else:
        if loading == "fit" and has_grain:
            # the alternation did not settle (f_xy moves the grain, L_Ph the straw): solve the
            # two cells jointly from where it stopped
            f_fit, L, m, status["joint"] = solve_joint(c, shoot, f_fit, L, t_straw,
                                                      np.log10(o["grain"]) - lo_root, runs)
    f, gx, gp = conductances(c, loading, p, f_fit)
    r_eq = float(10.0 ** (lo_root - m["root"]))
    res = {"root": 0.0, "straw": float(m["straw"] - m["root"] - t_straw)}
    res["grain"] = (float(m["grain"] - m["root"] - (np.log10(o["grain"]) - lo_root))
                    if has_grain else None)
    return c, dict(r_eq=r_eq, L_Ph=L, f_xy=float(f), g_xy=float(gx), g_ph=float(gp),
                   status=status, res=res, logs={k: float(m[k] - m["root"] + lo_root)
                                                 for k in TISS})


# ---------------------------------------------------------------- self checks
def _check_job(args):
    c, shoot = args
    p = START
    f, gx, gp = conductances(c, "lipid", p)
    a = model_logs(c, shoot, f, gx, gp, 0.01, (0, 1), r_eq=1.0)
    b = model_logs(c, shoot, f, gx, gp, 0.01, (0, 1), r_eq=2.0)
    w = model_logs(c, shoot, f, gx, gp, 0.01, (0, 1), r_eq=1.0, Cwo=OBS[c]["water_ng_L"] / 1e3)
    return c, shoot, max(abs(b[k] - a[k] - np.log10(2.0)) for k in TISS), \
        max(abs(w[k] - a[k]) for k in TISS)


def self_checks(pool):
    out = list(pool.map(_check_job, [(c, s) for c in ("PFBA", "PFOS", "PFDoDA")
                                     for s in ("redistributed", "4pool")]))
    lin_req = max(o[2] for o in out)
    lin_cwo = max(o[3] for o in out)
    print(f"self-check: r_eq linearity max |dlog10| {lin_req:.2e}; BAF independent of Cwo "
          f"max |dlog10| {lin_cwo:.2e}")
    assert lin_req < 1e-5 and lin_cwo < 1e-5
    return dict(r_eq_linearity=lin_req, cwo_independence=lin_cwo)


# ---------------------------------------------------------------- per-congener fits
def fit_variant(pool, v, p=None, warm=None):
    shoot, loading = VARIANTS[v]
    jobs = [(c, shoot, loading, p, (0, 1), (warm or {}).get(c)) for c in CONG]
    return dict(pool.map(fit_congener, jobs))


def cultivar_intervals(pool, v, p=None, main=None):
    shoot, loading = VARIANTS[v]
    jobs = [(c, shoot, loading, p, (i,), (main or {}).get(c)) for c in CONG for i in (0, 1)]
    res = list(pool.map(fit_congener, jobs))
    out = {}
    for c, r in res:
        out.setdefault(c, []).append(dict(r_eq=r["r_eq"], L_Ph=r["L_Ph"], f_xy=r["f_xy"],
                                          status=r["status"]))
    return out


def cells(fits):
    """Profiled per-cell errors (log10): straw of every congener, grain of every congener with
    a grain value (0 where L_Ph solved it)."""
    e = []
    for c in CONG:
        e.append(fits[c]["res"]["straw"])
        if fits[c]["res"]["grain"] is not None:
            e.append(fits[c]["res"]["grain"])
    return np.array(e)


def score(fits):
    """Yamazaki observable: log10 RMSE over the used cells (root, straw, grain)."""
    by = {k: [] for k in TISS}
    for c in CONG:
        for k in TISS:
            if OBS[c]["use"][k]:
                by[k].append(fits[c]["logs"][k] - np.log10(OBS[c]["baf_dw"][k]))
    allr = np.concatenate([by[k] for k in TISS])
    rm = lambda a: float(np.sqrt(np.mean(np.square(a))))
    return dict(all=rm(allr), n=int(len(allr)), **{k: rm(by[k]) for k in TISS})


# ---------------------------------------------------------------- lipid constants
FORMS = {   # term -> how the reduced model sets it
    "g_ph": dict(g_ph=0.0), "g_xy": dict(g_xy=0.0), "gate_none": dict(gate="none"),
    "gate_linear": dict(gate="linear"), "pfsa_ln": dict(pfsa_ln=0.0), "b": dict(b=0.0),
}


def free_names(model):
    """Optimised parameters of a lipid model form (`model`: the terms set at their bound)."""
    names = ["a"]
    if model.get("b") != 0.0:
        names.append("b")
    has_xy, has_ph = model.get("g_xy") != 0.0, model.get("g_ph") != 0.0
    if has_xy:
        names.append("g_xy")
    if has_ph:
        names.append("g_ph")
    if (has_xy or has_ph) and model.get("gate", "half") == "half":
        names.append("K_half")
    if model.get("pfsa_ln") != 0.0:
        names.append("pfsa_ln")
    return names


LOGP = ("a", "g_xy", "g_ph", "K_half")


def pack(p, names):
    return np.array([np.log10(p[n]) if n in LOGP else p[n] for n in names], float)


def unpack(x, names, model, base):
    p = dict(base)
    p.update({k: v for k, v in model.items()})
    for n, v in zip(names, x):
        p[n] = float(10.0 ** v) if n in LOGP else float(v)
    return p


def bounds_for(names):
    lo = [(-np.inf if n != "g_xy" and n != "g_ph" else LOGG_MIN) if n in LOGP else -np.inf
          for n in names]
    hi = [0.0 if n == "a" else np.inf for n in names]
    return lo, hi


DLN = 0.01                               # relative step of the conductance derivatives


def _dln_dx(c, p, names):
    """d ln(f_xy, g_xy, g_ph) of congener c / d x (x: log10 for a, g, K_half; b, pfsa_ln)."""
    rec = api._CONG[c]
    n, K_PL, psa = rec["n_C"], rec["K_PL_Lkg"], rec["group"] == "PFSA"
    gate = p.get("gate", "half")
    J = np.zeros((3, len(names)))
    ln10 = np.log(10.0)
    for j, nm in enumerate(names):
        if nm == "a":
            J[0, j] = ln10
        elif nm == "b":
            J[0, j] = -(n - 4)
        elif nm == "pfsa_ln":
            J[:, j] = -1.0 if psa else 0.0
        elif nm == "g_xy":
            J[1, j] = ln10
        elif nm == "g_ph":
            J[2, j] = ln10
        elif nm == "K_half" and gate == "half":
            d = -ln10 * p["K_half"] / (K_PL + p["K_half"])
            J[1, j] = d
            J[2, j] = d
    return J


def _cell_rows(fit):
    rows = [fit["res"]["straw"]]
    if fit["res"]["grain"] is not None:
        rows.append(fit["res"]["grain"])
    return np.array(rows)


def fit_lipid(pool, shoot, model, start, warm):
    """Profiled least squares for one lipid model form (Jacobian by the chain rule through
    each congener's three conductances). Returns dict(p, fits, cells, rmse, nfev)."""
    names = free_names(model)
    state = {"warm": dict(warm), "x": None, "fits": None}

    def fits_at(x):
        if state["x"] is not None and np.array_equal(x, state["x"]):
            return state["fits"]
        p = unpack(x, names, model, start)
        fits = dict(pool.map(fit_congener, [(c, shoot, "lipid", p, (0, 1), state["warm"].get(c))
                                             for c in CONG]))
        state.update(warm=fits, x=np.array(x), fits=fits)
        return fits

    def resid(x):
        return cells(fits_at(x))

    def jac(x):
        fits = fits_at(x)
        p = unpack(x, names, model, start)
        jobs, keys = [], []
        for c in CONG:
            e = (fits[c]["f_xy"], fits[c]["g_xy"], fits[c]["g_ph"])
            for k in range(3):
                if e[k] > 0:
                    ek = list(e)
                    ek[k] = e[k] * np.exp(DLN)
                    jobs.append((c, shoot, "explicit", tuple(ek), (0, 1), fits[c]))
                    keys.append((c, k))
        pert = list(pool.map(fit_congener, jobs))
        dr = {c: np.zeros((len(_cell_rows(fits[c])), 3)) for c in CONG}
        for (c, k), (_, fk) in zip(keys, pert):
            dr[c][:, k] = (_cell_rows(fk) - _cell_rows(fits[c])) / DLN
        return np.vstack([dr[c] @ _dln_dx(c, p, names) for c in CONG])

    x0 = pack(start, names)
    lo, hi = bounds_for(names)
    x0 = np.clip(x0, np.array(lo) + 1e-9, np.array(hi) - 1e-9)
    sol = least_squares(resid, x0, jac=jac, bounds=(lo, hi), x_scale=1.0,
                        xtol=1e-6, ftol=1e-8, gtol=1e-8, max_nfev=200)
    p = unpack(sol.x, names, model, start)
    fits = fits_at(sol.x)
    e = cells(fits)
    return dict(p=p, names=names, model=dict(model), fits=fits, cells=e.tolist(),
                rmse=float(np.sqrt(np.mean(e ** 2))), nfev=int(sol.nfev), njev=int(sol.njev),
                status=int(sol.status),
                at_logg_min=[n for n in names if n in ("g_xy", "g_ph")
                             and np.log10(p[n]) <= LOGG_MIN + 1e-6])


def eliminate(pool, shoot, log):
    """Backward elimination of the lipid terms (PREREG 5.6)."""
    model = {}
    cur = fit_lipid(pool, shoot, model, START, {})
    log.append(dict(step="full", model=dict(model), p=cur["p"], rmse=cur["rmse"]))
    print(f"  [{shoot}] full: rmse {cur['rmse']:.4f}  {fmt(cur['p'], cur['names'])}", flush=True)
    while True:
        cand = []
        has_g = model.get("g_xy") != 0.0 or model.get("g_ph") != 0.0
        for term, setting in FORMS.items():
            if term.startswith("gate"):
                if model.get("gate") or not has_g:
                    continue
            elif model.get(term) == 0.0:
                continue
            red = dict(model, **setting)
            st = dict(cur["p"])
            if term == "gate_linear":                    # g per L/kg: g_half / K_half
                st.update(g_xy=st["g_xy"] / st["K_half"], g_ph=st["g_ph"] / st["K_half"])
            r = fit_lipid(pool, shoot, red, st, cur["fits"])
            d = D.decide(np.array(cur["cells"]), np.array(r["cells"]))
            cand.append((term, r, d))
            print(f"  [{shoot}]   - {term:12s} rmse {r['rmse']:.4f}  {d['verdict']:6s} "
                  f"[{d['lo']:+.4f}, {d['hi']:+.4f}]", flush=True)
        removable = [(t, r, d) for t, r, d in cand if d["verdict"] != "worse"]
        log.append(dict(step="candidates", model=dict(model),
                        cand=[dict(term=t, rmse=r["rmse"], p=r["p"], decide=d)
                              for t, r, d in cand]))
        if not removable:
            break
        t, r, d = min(removable, key=lambda z: z[1]["rmse"])
        model = r["model"]
        cur = r
        log.append(dict(step="remove", term=t, model=dict(model), p=cur["p"], rmse=cur["rmse"]))
        print(f"  [{shoot}] removed {t}: rmse {cur['rmse']:.4f}  {fmt(cur['p'], cur['names'])}",
              flush=True)
        if not free_names(model)[1:]:
            break
    return cur


def fmt(p, names):
    return " ".join(f"{n}={p[n]:.4g}" for n in names) + (
        f" gate={p.get('gate', 'half')}" if "gate" in p else "")


def lipid_record(res):
    """The stored constants of a lipid model (removed terms at their bound)."""
    p, m = res["p"], res["model"]
    rec = dict(a=p["a"], b=p.get("b", 0.0) if m.get("b") != 0.0 else 0.0,
               g_xy=p["g_xy"] if m.get("g_xy") != 0.0 else 0.0,
               g_ph=p["g_ph"] if m.get("g_ph") != 0.0 else 0.0,
               pfsa_ln=p["pfsa_ln"] if m.get("pfsa_ln") != 0.0 else 0.0,
               gate=m.get("gate", "half"))
    if not (rec["g_xy"] or rec["g_ph"]):
        rec["gate"] = "none"                   # no bound term left: the gate is void
    elif rec["gate"] == "half":
        rec["K_half"] = p["K_half"]
    rec["terms_removed"] = sorted(k for k in ("b", "g_xy", "g_ph", "pfsa_ln")
                                  if m.get(k) == 0.0) + (
        [f"gate_{m['gate']}"] if m.get("gate") else [])
    return rec


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true")
    ap.add_argument("--write-only", action="store_true",
                    help="write the parameter files from an existing refit_sourced.json")
    ap.add_argument("--variants", nargs="*", default=list(VARIANTS))
    ap.add_argument("--resume", action="store_true",
                    help="keep the variants an interrupted run already wrote to refit_sourced.json")
    a = ap.parse_args()
    if a.write_only:
        with open(OUT) as f:
            rec = json.load(f)
        if not rec.get("complete"):
            raise SystemExit("refit_sourced.json is a partial run (re-run with --resume)")
        write_params(rec)
        return
    t0 = time.time()
    rec = dict(protocol="docs/PREREG_sourced_values.md 5.3 + 5.6", variants={})
    if a.resume and os.path.exists(OUT):
        with open(OUT) as f:
            rec = json.load(f)
        t0 -= rec.get("seconds", 0.0)
        print(f"resuming: kept {list(rec['variants'])}", flush=True)
    with ProcessPoolExecutor(max_workers=os.cpu_count() or 4) as pool:
        if "self_checks" not in rec:
            rec["self_checks"] = self_checks(pool)
        for v in a.variants:
            if v in rec["variants"]:
                continue
            shoot, loading = VARIANTS[v]
            print(f"\n== {v} ==", flush=True)
            out = dict(shoot=shoot, loading=loading)
            p = None
            if loading == "lipid":
                log = []
                best = eliminate(pool, shoot, log)
                p = lipid_record(best)
                out.update(lipid=p, elimination=log)
                fits = best["fits"]
            else:
                fits = fit_variant(pool, v)
            out["fits"] = fits
            out["score"] = score(fits)
            out["cultivar"] = cultivar_intervals(pool, v, p, fits)
            rec["variants"][v] = out
            sc = out["score"]
            print(f"  Yamazaki observable log10 RMSE {sc['all']:.3f} (root {sc['root']:.3f}, "
                  f"straw {sc['straw']:.3f}, grain {sc['grain']:.3f}; n {sc['n']})")
            for c in CONG:
                f = fits[c]
                iv = out["cultivar"][c]
                print(f"  {c:7s} r_eq {f['r_eq']:9.4g} [{min(x['r_eq'] for x in iv):.4g}, "
                      f"{max(x['r_eq'] for x in iv):.4g}]  L_Ph {f['L_Ph']:9.4g} "
                      f"[{min(x['L_Ph'] for x in iv):.4g}, {max(x['L_Ph'] for x in iv):.4g}] "
                      f"({f['status'].get('L_Ph')})  f_xy {f['f_xy']:.4g}  "
                      f"res straw {f['res']['straw']:+.3f} grain "
                      f"{'   -  ' if f['res']['grain'] is None else format(f['res']['grain'], '+.3f')}",
                      flush=True)
            # checkpoint after every variant: a lost run resumes instead of restarting
            rec["seconds"] = time.time() - t0
            rec["complete"] = False
            with open(OUT, "w") as f:
                json.dump(rec, f, indent=1, default=float)
    rec["seconds"] = time.time() - t0
    rec["complete"] = all(v in rec["variants"] for v in VARIANTS)
    with open(OUT, "w") as f:
        json.dump(rec, f, indent=1, default=float)
    print(f"\n[written] {os.path.relpath(OUT, ROOT)}  ({rec['seconds'] / 60:.1f} min)")
    if a.write:
        write_params(rec)


def write_params(rec):
    import csv
    keys = []
    for v, out in rec["variants"].items():
        keys += [f"r_eq__{v}", f"L_Ph__{v}"] + ([f"f_xy__{v}"] if v.endswith("_fit") else [])
        if any(c in NO_GRAIN for c in CONG):
            keys.append(f"L_Ph_MQL__{v}")
        keys.append(f"L_Ph_unreachable__{v}")
    path = os.path.join(ROOT, "params", "refit_sourced.csv")
    with open(path, "w", newline="") as f:
        f.write("# validation/refit_sourced.py (pre-registered: docs/PREREG_sourced_values.md 5.3, 5.6).\n"
                "# Per variant <shoot>_<lipid|free|fit>: r_eq, L_Ph (and f_xy for *_fit) fitted to Yamazaki 2023 root and\n"
                "# grain (fit: also straw) on the 2020 Tsukuba site-year. L_Ph_MQL: PFBS's censoring bound (L_Ph = its half).\n"
                "# L_Ph_unreachable = 1: no L_Ph reaches the observed grain; L_Ph is the maximum of the model's grain.\n")
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["compound"] + keys)
        for c in CONG:
            row = [c]
            for k in keys:
                name, v = k.split("__")
                fit = rec["variants"][v]["fits"][c]
                if name == "L_Ph_unreachable":
                    val = 1.0 if fit["status"].get("L_Ph") == "unreachable" else None
                elif name == "L_Ph_MQL":
                    val = fit["status"].get("L_MQL")
                else:
                    val = fit.get(name)
                row.append("" if val is None or not np.isfinite(val) else f"{val:.8g}")
            w.writerow(row)
    lip = {VARIANTS[v][0]: out["lipid"] for v, out in rec["variants"].items() if "lipid" in out}
    with open(os.path.join(ROOT, "params", "lipid_loading_sourced.json"), "w") as f:
        json.dump(dict(source="validation/refit_sourced.py (PREREG 5.3, 5.6)",
                       lipid_loading=lip), f, indent=1)
    print("[written] params/refit_sourced.csv, params/lipid_loading_sourced.json -- now run "
          "python build_parameters.py")


if __name__ == "__main__":
    main()
