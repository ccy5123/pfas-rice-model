#!/usr/bin/env python3
# =============================================================================
# validation/refit_site_year.py
# -----------------------------------------------------------------------------
# RE-FIT the PFAS transport parameters on the Yamazaki 2020 SITE-YEAR, against the
# Yamazaki tissue data on the DRY basis, with the harvest observation operator.
#
# Why: every earlier fit (W2 on a placeholder logistic, *_oryza on the reduced ORYZA
# driver) ran on drivers that are not Yamazaki's season -- a 120-d climatology,
# transpiration borrowed from other sites (8.0 L/hill vs 16.5 for the site-year),
# and a fixed season length -- compared the model's GREEN-leaf straw with a harvested
# straw that holds dead blades, and compared fresh-weight values made from ASSUMED
# water contents. All three are corrected here:
#   drivers   src/site_scenario.py: JMA Tsukuba 2020, ORYZA2000 port (IR72.D92) with
#             phenology on the paper's heading/harvest dates, ORYZA2000 transpiration;
#   operator  src/harvest.py: attached dead blades kept in the straw; dry-basis
#             concentrations (burden / dry mass);
#   data      data_obs/obs_baf_Yamazaki.csv on the DRY basis (the fresh-weight table's
#             water assumptions inverted exactly); PFBS grain (wholly imputed) unused.
#
# PRE-REGISTERED (this header is committed and pushed before the script is run).
#
# PROTOCOL
#   Model output   per congener and harvest run: the dry-basis harvest BAF (root; straw
#                  = stem + green leaf + attached dead blades; brown-rice grain), with
#                  Cwo = the irrigation-water concentration of that congener (constant).
#   Two runs       harvest 2020-09-01 and 2020-09-10 (the two cultivars), combined by
#                  geometric mean -- as the observations combine the two cultivars.
#   Fit A ("fit")  per congener (f_xy, L_Ph, kappa_d) -> (straw, grain, root): 1-D brentq
#                  in log space, kappa_d<-root, f_xy<-straw, L_Ph<-grain, repeated until
#                  every residual is < 0.002 log10 or 10 passes; a target outside the
#                  bounds clamps at the bound. Bounds as every earlier fit
#                  (structural_adequacy_fit.py): f_xy [1e-3, 1], L_Ph [1e-5, 1],
#                  kappa_d [1e-2, 50]. PFBS has no grain value: its L_Ph is the geometric
#                  mean of the other congeners' fitted L_Ph, and its predicted grain is
#                  checked against the MQL (2 x the imputed value).
#   Fit B          f_xy = f_xy_recommended (the monotone a-priori TSCF, not fitted);
#   ("recommended") (L_Ph, kappa_d) -> (grain, root) as above. Straw is then a PREDICTION.
#   Lipid fit      the six globals of the K_PL-gated lipid model (a, b, g_xy, g_ph, K_half,
#                  pfsa_ln; model_api.lipid_loading_conductances) by log-space least
#                  squares on all used cells, PFDoDA excluded (as before), with per-
#                  congener (L_Ph, kappa_d) from Fit B, starting from the current values.
#   Unknown inputs predictions with Fit A held fixed, one input at a time: nursery SBDUR
#                  14 / 30 d (default 23), covered seedbed TMPSB 9 C (default open, 0),
#                  density 15 / 22 hills m-2 (default 19.05; 30 x 15-22 cm), 5 plants
#                  per hill (default 3), heading-to-harvest 30 / 31 d (default 30.5),
#                  harvest NOT at maturity (IR72's own DVRR kept), CO2 400 / 430 ppm
#                  (default 414), dead blades all lost (attached 0; default 1), and the
#                  single mean-date run instead of the two-run geometric mean.
#   Scale          the data-anchored reference for "consequential" is the model difference
#                  between the two observed harvest runs, |log10 BAF(run 1/run 2)|, per
#                  congener and tissue. An unknown is CONSEQUENTIAL for a tissue if the
#                  median over congeners of its larger |delta log10 BAF| exceeds the median
#                  reference difference for that tissue. Each consequential unknown is then
#                  re-fitted (Fit A) at both ends and the parameter changes reported.
#
# DECISIONS (fixed before the results)
#   D1  Fit A becomes the fitted transport set in params/parameters.json (f_xy_fit,
#       L_Ph_fit, kappa_d_fit), REPLACING *_W2fit and *_oryza(_dry); Fit B supplies L_Ph /
#       kappa_d for f_xy_source="recommended" (L_Ph_rec, kappa_d_rec); the lipid globals
#       replace LIPID_LOADING(_DRY). This holds whatever the fit quality: the old values
#       were fitted on drivers and an operator now known to be wrong.
#   D2  reproduce_demo.py runs Fit A on this scenario; its RMSE is a REPRODUCTION number
#       (saturated, 3 parameters / 3 cells per congener), stated as such.
#   (The default scenario switch to this site-year was decided on the ORYZA2000 port's
#   evidence before this fit -- recorded in CLAUDE.md, not pre-registered here.)
#
# EXPECTATIONS (before running)
#   E1  Fit A reproduces every used cell within 0.01 log10 except where a bound binds;
#       PFDoDA hits f_xy = 1 with its straw under, as on every earlier driver; a short
#       chain may sit at the L_Ph floor with its grain over (xylem alone over-feeds it).
#   E2  fitted f_xy below f_xy_oryza for most congeners: twice the season transpiration
#       and the dead blades both raise straw for a given f_xy (Q_TP x f_xy is a ridge).
#   E3  kappa_d close to the earlier fits (the root equilibrates largely on its own).
#   E4  Fit B's a-priori straw log10 RMSE between 0.5 and 1.2 (earlier a-priori 0.84-0.95).
#   E5  harvest-not-at-maturity and dead-blades-lost are the consequential unknowns for
#       straw and grain; nursery, density, plants per hill and CO2 are not.
#
#   python validation/refit_site_year.py            # fits + sensitivity, prints the verdicts
#   python validation/refit_site_year.py --write    # ... and applies D1 to parameters.json
# =============================================================================
from __future__ import annotations

import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.optimize import brentq, least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
import model_api as api        # noqa: E402
import site_scenario as ss     # noqa: E402

TISS = ("root", "straw", "grain")
FXY, LPH, KD = (1e-3, 1.0), (1e-5, 1.0), (1e-2, 50.0)
TOL, MAX_PASS = 0.002, 10
HARVESTS = ss.YAMAZAKI_HARVESTS
OBS = api.yamazaki_obs()
CONG = [c for c in api.CONGENERS if c in OBS]
LIPID_NAMES = ("a", "b", "g_xy", "g_ph", "K_half", "pfsa_ln")
LIPID_LOG = ("a", "g_xy", "g_ph", "K_half")
LIPID_BOUNDS = dict(a=(1e-3, 1.0), b=(0.0, 1.5), g_xy=(1e-4, 1.0), g_ph=(1e-5, 0.5),
                    K_half=(1e2, 1e6), pfsa_ln=(0.0, 3.0))
LIPID_EXCLUDED = ("PFDoDA",)
UNKNOWNS = {                                      # name -> (low setting, high setting)
    "sbdur": ({"sbdur": 14}, {"sbdur": 30}),
    "tmpsb": ({"tmpsb": 9.0}, None),
    "density": ({"nh": 15.0}, {"nh": 22.0}),
    "plants_per_hill": ({"nplh": 5.0}, None),
    "heading_to_harvest": ({"days_heading_to_harvest": 30.0}, {"days_heading_to_harvest": 31.0}),
    "harvest_not_maturity": ({"harvest_is_maturity": False}, None),
    "co2": ({"co2": 400.0}, {"co2": 430.0}),
    "dead_blades_lost": ("attached0", None),
    "single_mean_run": ("mean_run", None),
}
OUT = os.path.join(HERE, "refit_site_year.json")


def used(c):
    d = OBS[c]
    return {k: v for k, v in d["baf_dw"].items() if d["use"][k]}


def predict(c, f_xy, L_Ph, kappa_d, unknowns=None, attached=1.0, runs="both", lipid=None):
    """Two-run geometric-mean dry-basis harvest BAF {root, straw, grain}."""
    unknowns = dict(unknowns or {})
    cwo = OBS[c]["water_ng_L"] / 1000.0
    hs = ["mean"] if runs == "mean" else list(HARVESTS)
    kw = {}
    if lipid is not None:
        p = lipid
        rec = api._CONG[c]
        sf = np.exp(-p["pfsa_ln"]) if rec["group"] == "PFSA" else 1.0
        phi = rec["K_PL_Lkg"] / (rec["K_PL_Lkg"] + p["K_half"])
        f_xy = p["a"] * np.exp(-p["b"] * (rec["n_C"] - 4)) * sf
        kw = dict(g_xy_override=p["g_xy"] * phi * sf, g_ph_override=p["g_ph"] * phi * sf)
    logs = {k: [] for k in TISS}
    for h in hs:
        d = api.scenario_drivers(Cwo=cwo, harvest=h, **unknowns)
        r = api.simulate(c, drivers=d, f_xy_override=f_xy, L_Ph_override=L_Ph,
                         kappa_d_override=kappa_d, dead_leaf_attached=attached, **kw)
        b = r["baf_final_dw"]
        for k in TISS:
            logs[k].append(np.log10(max(b[k], 1e-30)))
    return {k: float(10 ** np.mean(v)) for k, v in logs.items()}


def _fit_scalar(c, tissue, var, lo, hi, p, target, **kw):
    """1-D brentq in log space for parameter `var` so that the model's `tissue` hits
    `target`; clamped at a bound when the target is outside the reachable range."""
    def g(x):
        q = dict(p, **{var: 10 ** x})
        return np.log10(predict(c, q["f_xy"], q["L_Ph"], q["kappa_d"], **kw)[tissue]) - np.log10(target)
    a, b = np.log10(lo), np.log10(hi)
    ga, gb = g(a), g(b)
    if ga * gb > 0:
        return lo if abs(ga) < abs(gb) else hi
    return float(10 ** brentq(g, a, b, xtol=1e-4))


def fit_one(c, free=("kappa_d", "f_xy", "L_Ph"), L_Ph_fixed=None, unknowns=None):
    """Sequential saturated fit. `free` lists the parameters fitted (others held);
    returns params, prediction, residuals, passes."""
    rec = api._CONG[c]
    o = used(c)
    p = dict(f_xy=float(rec["f_xy_recommended"]), L_Ph=0.003 if L_Ph_fixed is None else L_Ph_fixed,
             kappa_d=2.0)
    which = {"kappa_d": ("root", KD), "f_xy": ("straw", FXY), "L_Ph": ("grain", LPH)}
    kw = dict(unknowns=unknowns)
    for n in range(1, MAX_PASS + 1):
        for var in ("kappa_d", "f_xy", "L_Ph"):
            if var not in free:
                continue
            tissue, (lo, hi) = which[var]
            if tissue not in o:
                continue
            p[var] = _fit_scalar(c, tissue, var, lo, hi, p, o[tissue], **kw)
        pred = predict(c, p["f_xy"], p["L_Ph"], p["kappa_d"], **kw)
        res = {k: float(np.log10(pred[k] / o[k])) for k in o}
        fitted_t = [which[v][0] for v in free if which[v][0] in o]
        if all(abs(res[k]) < TOL for k in fitted_t) or n == MAX_PASS:
            break
        # stop early when a bound binds and nothing moves any more
    at_bound = {v: bool(np.isclose(p[v], which[v][1][0]) or np.isclose(p[v], which[v][1][1]))
                for v in free}
    return dict(params=p, pred=pred, res=res, passes=n, at_bound=at_bound)


def _fit_a(c, L_Ph_fixed=None, unknowns=None):
    free = ("kappa_d", "f_xy") if L_Ph_fixed is not None else ("kappa_d", "f_xy", "L_Ph")
    return c, fit_one(c, free=free, L_Ph_fixed=L_Ph_fixed, unknowns=unknowns)


def _fit_b(c, L_Ph_fixed=None):
    free = ("kappa_d",) if L_Ph_fixed is not None else ("kappa_d", "L_Ph")
    return c, fit_one(c, free=free, L_Ph_fixed=L_Ph_fixed)


def _star(args):
    f, a = args
    return f(*a)


def run_fits(pool, kind, unknowns=None):
    """Fit A or B for every congener; PFBS (no grain) last with the others' geomean L_Ph."""
    first = [c for c in CONG if "grain" in used(c)]
    rest = [c for c in CONG if "grain" not in used(c)]
    f = _fit_a if kind == "A" else _fit_b
    args = [(f, (c, None, unknowns)) if kind == "A" else (f, (c, None)) for c in first]
    out = dict(pool.map(_star, args))
    gm = float(10 ** np.mean([np.log10(out[c]["params"]["L_Ph"]) for c in first]))
    args = [(f, (c, gm, unknowns)) if kind == "A" else (f, (c, gm)) for c in rest]
    out.update(dict(pool.map(_star, args)))
    return out, gm


def _pred_lipid(args):
    c, p, lph, kd = args
    return c, predict(c, None, lph, kd, lipid=p)


def lipid_fit(pool, fit_b, start):
    cong = [c for c in CONG if c not in LIPID_EXCLUDED]
    names = LIPID_NAMES

    def unpack(x):
        return {n: (10 ** v if n in LIPID_LOG else v) for n, v in zip(names, x)}

    def resid(x):
        p = unpack(x)
        jobs = [(c, p, fit_b[c]["params"]["L_Ph"], fit_b[c]["params"]["kappa_d"]) for c in cong]
        preds = dict(pool.map(_pred_lipid, jobs))
        return np.array([np.log10(preds[c][k] / used(c)[k]) for c in cong for k in used(c)])

    x0 = [np.log10(start[n]) if n in LIPID_LOG else start[n] for n in names]
    lo = [np.log10(LIPID_BOUNDS[n][0]) if n in LIPID_LOG else LIPID_BOUNDS[n][0] for n in names]
    hi = [np.log10(LIPID_BOUNDS[n][1]) if n in LIPID_LOG else LIPID_BOUNDS[n][1] for n in names]
    x0 = np.clip(x0, np.array(lo) + 1e-9, np.array(hi) - 1e-9)
    sol = least_squares(resid, x0, bounds=(lo, hi), diff_step=1e-2, x_scale=1.0, max_nfev=300)
    p = unpack(sol.x)
    r = resid(sol.x)
    return dict(params=p, rmse=float(np.sqrt(np.mean(r ** 2))), n=len(r), nfev=int(sol.nfev))


def _pred_unknown(args):
    c, p, setting = args
    if setting == "attached0":
        return c, predict(c, p["f_xy"], p["L_Ph"], p["kappa_d"], attached=0.0)
    if setting == "mean_run":
        return c, predict(c, p["f_xy"], p["L_Ph"], p["kappa_d"], runs="mean")
    return c, predict(c, p["f_xy"], p["L_Ph"], p["kappa_d"], unknowns=setting)


def _per_run(args):
    c, p = args
    out = []
    for h in HARVESTS:
        d = api.scenario_drivers(Cwo=OBS[c]["water_ng_L"] / 1000.0, harvest=h)
        r = api.simulate(c, drivers=d, f_xy_override=p["f_xy"], L_Ph_override=p["L_Ph"],
                         kappa_d_override=p["kappa_d"])
        out.append({k: r["baf_final_dw"][k] for k in TISS})
    return c, out


def sensitivity(pool, fit_a):
    P = {c: fit_a[c]["params"] for c in CONG}
    base = {c: fit_a[c]["pred"] for c in CONG}
    runs = dict(pool.map(_per_run, [(c, P[c]) for c in CONG]))
    ref = {k: float(np.median([abs(np.log10(runs[c][0][k] / runs[c][1][k])) for c in CONG]))
           for k in TISS}
    eff = {}
    for name, ends in UNKNOWNS.items():
        deltas = {k: [] for k in TISS}
        per_end = []
        for setting in ends:
            if setting is None:
                continue
            pr = dict(pool.map(_pred_unknown, [(c, P[c], setting) for c in CONG]))
            per_end.append({c: {k: float(np.log10(pr[c][k] / base[c][k])) for k in TISS} for c in CONG})
        for k in TISS:
            deltas[k] = [max(abs(e[c][k]) for e in per_end) for c in CONG]
        med = {k: float(np.median(deltas[k])) for k in TISS}
        eff[name] = dict(median_abs_delta=med, consequential={k: med[k] > ref[k] for k in TISS},
                         per_end=per_end)
    return dict(reference=ref, effects=eff, per_run=runs)


def rmse(fit, cells=None):
    e = [fit[c]["res"][k] for c in fit for k in fit[c]["res"] if cells is None or k in cells]
    return float(np.sqrt(np.mean(np.square(e)))) if e else float("nan")


def main(write=False, verbose=True):
    say = print if verbose else (lambda *a, **k: None)
    with ProcessPoolExecutor(max_workers=os.cpu_count() or 2) as pool:
        A, gmA = run_fits(pool, "A")
        B, gmB = run_fits(pool, "B")
        L = lipid_fit(pool, B, api.LIPID_LOADING)
        S = sensitivity(pool, A)
        refits = {}
        for name, e in S["effects"].items():
            if not any(e["consequential"].values()):
                continue
            for i, setting in enumerate(UNKNOWNS[name]):
                if setting is None or isinstance(setting, str):
                    continue
                refits[f"{name}[{i}]"] = {c: v["params"] for c, v in run_fits(pool, "A", setting)[0].items()}

    say("Fit A (f_xy, L_Ph, kappa_d) on the Yamazaki 2020 site-year, dry basis, harvest operator")
    say(f"{'PFAS':8}{'f_xy':>9}{'L_Ph':>10}{'kappa_d':>9} | {'root m/o':>16}{'straw m/o':>16}{'grain m/o':>16} | bound")
    for c in CONG:
        r = A[c]; p = r["params"]; o = OBS[c]["baf_dw"]
        b = ",".join(v for v, x in r["at_bound"].items() if x) or "-"
        say(f"{c:8}{p['f_xy']:>9.4f}{p['L_Ph']:>10.5f}{p['kappa_d']:>9.3f} | "
            + "".join(f"{r['pred'][k]:>7.3g}/{o[k]:<8.3g}" for k in TISS) + f" | {b}")
    say(f"  reproduction log10 RMSE {rmse(A):.3f} (saturated); PFBS L_Ph = geomean {gmA:.4g}; "
        f"PFBS grain {A['PFBS']['pred']['grain']:.3g} vs MQL bound {2 * OBS['PFBS']['baf_dw']['grain']:.3g}")
    say("\nFit B (f_xy = monotone a-priori; L_Ph, kappa_d fitted): straw is a prediction")
    for c in CONG:
        r = B[c]; p = r["params"]
        say(f"{c:8}{p['f_xy']:>9.4f}{p['L_Ph']:>10.5f}{p['kappa_d']:>9.3f} | straw "
            f"{r['pred']['straw']:.3g} vs {OBS[c]['baf_dw']['straw']:.3g} "
            f"[{OBS[c]['straw_lo']:.3g}, {OBS[c]['straw_hi']:.3g}]  ({r['res']['straw']:+.2f})")
    inside = sum(OBS[c]["straw_lo"] <= B[c]["pred"]["straw"] <= OBS[c]["straw_hi"] for c in CONG)
    say(f"  a-priori straw log10 RMSE {rmse(B, ['straw']):.3f}; inside the band range {inside}/{len(CONG)}")
    say(f"\nLipid fit (6 globals, Fit B L_Ph/kappa_d, PFDoDA excluded): RMSE {L['rmse']:.3f} over {L['n']} cells")
    say("  " + ", ".join(f"{k}={v:.4g}" for k, v in L["params"].items()))
    say("\nUnknown inputs (Fit A held): median |delta log10 BAF| vs the two-run reference")
    say(f"  reference (run 1 vs run 2): " + ", ".join(f"{k} {v:.3f}" for k, v in S["reference"].items()))
    for name, e in S["effects"].items():
        say(f"  {name:22}" + "".join(f"{k} {e['median_abs_delta'][k]:.3f}{'*' if e['consequential'][k] else ' '}  "
                                    for k in TISS))
    for key, P in refits.items():
        say(f"  re-fit {key}: " + ", ".join(
            f"{c} f_xy x{P[c]['f_xy'] / A[c]['params']['f_xy']:.2f}" for c in ("PFBA", "PFOA", "PFOS", "PFDoDA")))

    # verdicts against the expectations
    e1 = [c for c in CONG if max(abs(v) for v in A[c]["res"].values()) > 0.01]
    say("\nE1 cells off by >0.01:", e1 or "none", "| bounds:",
        {c: [v for v, x in A[c]["at_bound"].items() if x] for c in CONG if any(A[c]["at_bound"].values())})
    lower = [c for c in CONG if A[c]["params"]["f_xy"] < (api._CONG[c].get("f_xy_oryza") or np.inf)]
    say(f"E2 f_xy below f_xy_oryza: {len(lower)}/{len(CONG)}")
    say("E3 kappa_d / kappa_d_oryza: " + ", ".join(
        f"{c} {A[c]['params']['kappa_d'] / api._CONG[c]['kappa_d_oryza']:.2f}" for c in CONG
        if api._CONG[c].get("kappa_d_oryza")))
    say(f"E4 a-priori straw RMSE {rmse(B, ['straw']):.3f} (expected 0.5-1.2)")
    cons = [n for n, e in S["effects"].items() if any(e["consequential"].values())]
    say(f"E5 consequential unknowns: {cons}")

    out = dict(A={c: {k: A[c][k] for k in ("params", "pred", "res", "at_bound", "passes")} for c in CONG},
               B={c: {k: B[c][k] for k in ("params", "pred", "res", "at_bound", "passes")} for c in CONG},
               L_Ph_geomean_A=gmA, L_Ph_geomean_B=gmB, lipid=L,
               sensitivity=dict(reference=S["reference"],
                                effects={n: {k: e[k] for k in ("median_abs_delta", "consequential")}
                                         for n, e in S["effects"].items()}),
               refits=refits, rmse_A=rmse(A), rmse_B_straw=rmse(B, ["straw"]))
    with open(OUT, "w") as f:
        json.dump(out, f, indent=1)
    if write:
        apply(out)
    return out


def apply(out):
    """D1: write the fitted sets into params/parameters.json (and drop the superseded ones)."""
    path = os.path.join(ROOT, "params", "parameters.json")
    par = json.load(open(path))
    drop = ("f_xy_W2fit", "L_Ph_W2fit", "kappa_d_W2fit", "f_xy_oryza", "L_Ph_oryza", "kappa_d_oryza",
            "f_xy_oryza_dry", "L_Ph_oryza_dry", "kappa_d_oryza_dry")
    for c in par["congeners"]:
        for k in drop:
            c.pop(k, None)
        if c["name"] in out["A"]:
            a, b = out["A"][c["name"]]["params"], out["B"][c["name"]]["params"]
            c["f_xy_fit"], c["L_Ph_fit"], c["kappa_d_fit"] = (round(a[k], 6) for k in ("f_xy", "L_Ph", "kappa_d"))
            c["L_Ph_rec"], c["kappa_d_rec"] = round(b["L_Ph"], 6), round(b["kappa_d"], 6)
    meta = par.setdefault("_meta", {})
    for k in ("refit_oryza", "refit_oryza_dry"):
        meta.pop(k, None)
    meta["refit_site_year"] = (
        "Transport parameters fitted to Yamazaki 2023 on the 2020 Tsukuba site-year (src/site_scenario.py), "
        "dry-basis tissue data, harvest operator with attached dead blades (src/harvest.py); "
        "validation/refit_site_year.py. *_fit: per congener (f_xy, L_Ph, kappa_d), saturated "
        f"(reproduction log10 RMSE {out['rmse_A']:.3f}). *_rec: (L_Ph, kappa_d) with the monotone "
        "f_xy_recommended fixed, so straw is a prediction (a-priori straw log10 RMSE "
        f"{out['rmse_B_straw']:.3f}). PFBS has no grain value: its L_Ph is the geometric mean of the others'.")
    json.dump(par, open(path, "w"), indent=2, ensure_ascii=False)


if __name__ == "__main__":
    main(write="--write" in sys.argv)
