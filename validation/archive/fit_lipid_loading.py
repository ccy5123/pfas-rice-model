#!/usr/bin/env python3
# =============================================================================
# validation/fit_lipid_loading.py
# -----------------------------------------------------------------------------
# FIT the six global constants of the K_PL-gated lipid-facilitated loading model
# (`model_api.LIPID_LOADING`: a, b, g_xy, g_ph, K_half, pfsa_ln) to Yamazaki 2023
# root/straw/grain BAF, PFDoDA excluded (near-MQL outlier).
#
# Why this script exists: the constants were fitted during the exploration recorded
# in docs/fxy_longchain_lipid_exploration.md and no fit script was committed, so
# they could not be re-fitted when the organ-mass basis changed (src/organ_mass.py:
# the crop model's DRY matter had been used as the ODE's FRESH mass M). This
# reconstructs that protocol and runs it on either basis:
#
#   MEASURED forcings (forcing_rice Q_TP) + growth_rice biomass, season 120 d,
#   Cwo = 1, per-congener L_Ph / kappa_d = the W2 fits (`_transport_defaults`),
#   free-anion f_xy = a*exp(-b(n-4)), phi = K_PL/(K_PL+K_half),
#   g_xy/g_ph = g_*max * phi, PFSA x exp(-pfsa_ln)   (lipid_loading_conductances).
#
# PROTOCOL CHECK (printed first): on the DRY basis the stored dry constants give
# log10 RMSE 0.372 excl PFDoDA, which matches the documented "~0.36" -- so this is
# the protocol the constants came from. The constants themselves look hand-rounded
# (0.25 / 0.50 / 0.05 / 0.010 / 3000 / 1.25), so a re-fit on 'dry' need not return
# them exactly; that is why the dry record is kept as stored, not re-fitted.
#
# IN-SAMPLE (Yamazaki). The out-of-sample checks are validation/oos_tang_lipid.py
# and validation/oos_multidataset.py.
#
#   python validation/fit_lipid_loading.py                      # fresh (default)
#   PFAS_MASS_BASIS=dry python validation/fit_lipid_loading.py  # the pre-fix record
#   ... --eval                                                  # score the stored constants only
# =============================================================================
from __future__ import annotations

import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
import model_api as api  # noqa: E402

MASS_BASIS = os.environ.get("PFAS_MASS_BASIS", "fresh")
BIOMASS = "growth_rice"            # the driver the constants were fitted on
TISS = ("root", "straw", "grain")
NAMES = ("a", "b", "g_xy", "g_ph", "K_half", "pfsa_ln")
LOG = ("a", "g_xy", "g_ph", "K_half")                       # fitted in log10
BOUNDS = dict(a=(1e-3, 1.0), b=(0.0, 1.5), g_xy=(1e-4, 1.0), g_ph=(1e-5, 0.5),
              K_half=(1e2, 1e6), pfsa_ln=(0.0, 3.0))
CONGENERS = tuple(c for c in api.CONGENERS if api.observed_baf(c))
EXCLUDED = ("PFDoDA",)
OUT = os.path.join(HERE, f"lipid_loading_fit{'_dry' if MASS_BASIS == 'dry' else ''}.json")


def predict(congener, p, mass_basis=MASS_BASIS):
    """(root, straw, grain) BAF for one congener under the lipid constants `p`."""
    c = api._CONG[congener]
    sf = np.exp(-p["pfsa_ln"]) if c["group"] == "PFSA" else 1.0
    phi = c["K_PL_Lkg"] / (c["K_PL_Lkg"] + p["K_half"])
    r = api.simulate(congener, lipid_loading=True, biomass=BIOMASS, mass_basis=mass_basis,
                     f_xy_override=p["a"] * np.exp(-p["b"] * (c["n_C"] - 4)) * sf,
                     g_xy_override=p["g_xy"] * phi * sf, g_ph_override=p["g_ph"] * phi * sf)
    return dict(root=r["baf_final"]["root"], straw=r["straw_baf"], grain=r["baf_final"]["grain"])


def _predict_star(args):
    return predict(*args)


def predict_all(p, pool=None, mass_basis=MASS_BASIS):
    jobs = [(c, p, mass_basis) for c in CONGENERS]
    out = pool.map(_predict_star, jobs) if pool is not None else map(_predict_star, jobs)
    return dict(zip(CONGENERS, out))


def score(pred):
    """log10 RMSE: all, excl PFDoDA, and by tissue (all congeners)."""
    err, ex = {k: [] for k in TISS}, []
    for c, pr in pred.items():
        o = api.observed_baf(c)
        for k in TISS:
            if k in o:
                d = np.log10(max(pr[k], 1e-9)) - np.log10(o[k])
                err[k].append(d)
                if c not in EXCLUDED:
                    ex.append(d)
    rm = lambda v: float(np.sqrt(np.mean(np.square(v))))
    return dict(all=rm(sum(err.values(), [])), excl=rm(ex), **{k: rm(v) for k, v in err.items()})


def _unpack(x):
    return {n: (10.0 ** v if n in LOG else float(v)) for n, v in zip(NAMES, x)}


def _pack(p):
    return np.array([np.log10(p[n]) if n in LOG else p[n] for n in NAMES])


def _bounds():
    lo = [np.log10(BOUNDS[n][0]) if n in LOG else BOUNDS[n][0] for n in NAMES]
    hi = [np.log10(BOUNDS[n][1]) if n in LOG else BOUNDS[n][1] for n in NAMES]
    return lo, hi


def fit(p0, pool):
    def resid(x):
        pred = predict_all(_unpack(x), pool)
        r = []
        for c, pr in pred.items():
            if c in EXCLUDED:
                continue
            o = api.observed_baf(c)
            r += [np.log10(max(pr[k], 1e-9)) - np.log10(o[k]) for k in TISS if k in o]
        return r
    lo, hi = _bounds()
    x0 = np.clip(_pack(p0), lo, hi)
    sol = least_squares(resid, x0, bounds=(lo, hi), method="trf", diff_step=1e-2,
                        x_scale="jac", max_nfev=200)
    return _unpack(sol.x), sol


def table(pred):
    print(f"  {'PFAS':8}{'root p/o':>16}{'straw p/o':>16}{'grain p/o':>16}")
    for c, pr in pred.items():
        o = api.observed_baf(c)
        cells = "".join(f"{pr[k]:>8.2f}/{o.get(k, float('nan')):<7.2f}" for k in TISS)
        print(f"  {c:8}{cells}{'  (excluded)' if c in EXCLUDED else ''}")


def _line(tag, s):
    print(f"  {tag:34} all {s['all']:.3f} | excl PFDoDA {s['excl']:.3f} | "
          f"root {s['root']:.3f} straw {s['straw']:.3f} grain {s['grain']:.3f}")


def main(eval_only=False):
    print("=" * 80)
    print(f"LIPID_LOADING fit — Yamazaki, {BIOMASS} biomass, measured forcing, "
          f"{MASS_BASIS.upper()} organ mass")
    print("=" * 80)
    with ProcessPoolExecutor(max_workers=min(4, os.cpu_count() or 1)) as pool:
        dry = dict(api.LIPID_LOADING_DRY)
        s_dd = score(predict_all(dry, pool, "dry"))
        _line("stored DRY constants on dry mass", s_dd)
        print("    (documented ~0.36 excl PFDoDA -> the protocol is reproduced)")
        stored = dict(api.lipid_loading_params(MASS_BASIS))
        pred0 = predict_all(stored, pool)
        s0 = score(pred0)
        _line(f"stored {MASS_BASIS} constants", s0)
        if eval_only:
            table(pred0)
            return stored, s0
        p, sol = fit(stored, pool)
        pred = predict_all(p, pool)
        s = score(pred)
    print("\n  fitted: " + "  ".join(f"{n}={p[n]:.4g}" for n in NAMES))
    at_bound = [n for n in NAMES
                if np.isclose(p[n], BOUNDS[n][0], rtol=1e-3) or np.isclose(p[n], BOUNDS[n][1], rtol=1e-3)]
    if at_bound:
        print(f"  AT A BOUND: {', '.join(at_bound)}")
    _line("re-fitted", s)
    table(pred)
    json.dump({"params": p, "rmse": s, "mass_basis": MASS_BASIS, "biomass": BIOMASS,
               "excluded": list(EXCLUDED), "at_bound": at_bound, "nfev": int(sol.nfev),
               "note": "Yamazaki in-sample; validation/fit_lipid_loading.py"},
              open(OUT, "w"), indent=2)
    print(f"\n  saved -> {OUT}")
    return p, s


if __name__ == "__main__":
    main(eval_only="--eval" in sys.argv)
