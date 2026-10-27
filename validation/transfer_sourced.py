#!/usr/bin/env python3
# =============================================================================
# validation/transfer_sourced.py
# -----------------------------------------------------------------------------
# TRANSFER SCORES of every re-fitted variant, and the default decision
# (docs/PREREG_sourced_values.md 5.3; decision rule section 2).
#
# Nothing here is fitted: each variant runs with the values validation/refit_sourced.py
# fitted to Yamazaki (parameters.json r_eq__<v>, L_Ph__<v>, lipid_loading[shoot]) on the
# calibration site-year drivers (the transfer sites' drivers are not in the repo).
#
# Endpoints (items scored; log10 errors):
#   Tang 2026     per-organ TF (stalk, leaf, endosperm) of PFOA, PFOS, GenX at 0.1 ug/g  (9)
#   Kim 2019      grain BAF of PFHpA, PFNA, PFDA, PFUnDA, PFDoDA: point (published, non-detects
#                 as zero) and interval ([0, MDL] non-detects; 0 inside)                 (5)
#   Li 2025       stage-geometric-mean TF: shoot (stem, leaf) and grain (from filling)
#   Yamazaki      the calibration observable (in-sample for every variant)
# Decision (5.3): the incumbent default stays unless another variant is "better"
# (decision_stats.decide: paired bootstrap, 95 % interval below 0) on at least one transfer
# endpoint and "not worse" on every endpoint. The honesty band is the median transfer error
# of the resulting default (model_api.APRIORI_OOS_SCORES).
#
#   python validation/transfer_sourced.py      (writes validation/transfer_sourced.json)
# =============================================================================
from __future__ import annotations

import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, HERE)
import model_api as api                 # noqa: E402
import literature_params as lp          # noqa: E402
import decision_stats as D              # noqa: E402
import kim2019_censored as K            # noqa: E402
import li2025_stage_resolved as L       # noqa: E402

VARIANTS = {
    "redistributed_lipid": dict(shoot="redistributed", lipid_loading=True),
    "redistributed_free": dict(shoot="redistributed", lipid_loading=False),
    "redistributed_fit": dict(shoot="redistributed", lipid_loading=False, f_xy_source="fit"),
    "4pool_lipid": dict(shoot="4pool", lipid_loading=True),
    "4pool_free": dict(shoot="4pool", lipid_loading=False),
}
# the incumbent the amendment-2a decision was taken against (the record's "default" is the
# outcome: redistributed_fit, now model_api's default -- PREREG section 8)
DEFAULT = "redistributed_lipid"
TRANSFER = ("tang", "kim_point", "kim_interval", "li_shoot", "li_grain")
OUT = os.path.join(HERE, "transfer_sourced.json")


def tang_errors(kw):
    e = []
    for c in api.TANG_CONGENERS:
        v = api.tang_tf_validation(c, dose="low", shoot=kw["shoot"],
                                   lipid_loading=kw["lipid_loading"],
                                   f_xy_source=kw.get("f_xy_source", "recommended"))
        for org in v["organs"]:
            e.append(np.log10(max(v["model_tf"][org], 1e-12)) - np.log10(v["tang_tf"][org]))
    return e


def kim_errors(kw):
    b = lp.kim2019_grain_baf_bounds()
    ep, ei = [], []
    for nm in K.SCORED:
        g = float(api.simulate(nm, Cwo=1.0, **kw)["baf_final_dw"]["grain"])   # Kim: freeze-dried
        ep.append(K.point_error(g, b[nm]))
        ei.append(K.interval_error(g, b[nm]))
    return ep, ei


def li_errors(kw, grain="grain"):
    """Li 2025 level errors (shoot, grain); `grain` "rough" scores Li's undefined "grains" as
    brown rice + husk + rachis instead of brown rice (PREREG 7.2; reported, not decided on: 7.13)."""
    data = L.load_data()
    days = L.stage_days()
    res = {c: api.simulate(c, Cwo=1.0, **kw) for c in L.CONGENERS}
    if grain == "rough" and "rough" not in res[L.CONGENERS[0]]["conc_dw"]:
        return [], []
    model = {c: L.model_tf(res[c], days, grain) for c in L.CONGENERS}
    lev = {"shoot": [], "grain": []}
    for c in L.CONGENERS:
        for k in L.SHOOT + (L.GRAIN,):
            if (c, k) not in data:
                continue
            tf_o, _ = L.observed_tf(data, c, k)
            idx = [i for i in L.scored_stages(k) if not np.isnan(tf_o[i])]
            if len(idx) < 2:
                continue
            lo = np.log10(tf_o[idx])
            lm = np.log10(np.maximum(model[c][k][idx], 1e-12))
            lev["grain" if k == L.GRAIN else "shoot"].append(float(lm.mean() - lo.mean()))
    return lev["shoot"], lev["grain"]


def yamazaki_errors(kw):
    obs = api.yamazaki_obs()
    e = []
    for c in api.CONGENERS:
        if c not in obs:
            continue
        m = api.yamazaki_harvest_baf(c, **kw)
        for k in ("root", "straw", "grain"):
            if obs[c]["use"][k]:
                e.append(float(np.log10(m[k] / obs[c]["baf_dw"][k])))
    return e


def score_variant(v):
    kw = VARIANTS[v]
    ep, ei = kim_errors(kw)
    ls, lg = li_errors(kw)
    _, lgr = li_errors(kw, grain="rough")
    return v, dict(tang=tang_errors(kw), kim_point=ep, kim_interval=ei, li_shoot=ls,
                   li_grain=lg, li_grain_rough=lgr, yamazaki=yamazaki_errors(kw))


def rmse(e):
    return float(np.sqrt(np.mean(np.square(e))))


def main():
    with ProcessPoolExecutor(max_workers=min(len(VARIANTS), os.cpu_count() or 4)) as pool:
        err = dict(pool.map(score_variant, VARIANTS))
    ends = TRANSFER + ("yamazaki",)
    print(f"{'variant':22}" + "".join(f"{k:>14}" for k in ends))
    for v in VARIANTS:
        print(f"{v:22}" + "".join(f"{rmse(err[v][k]):>14.3f}" for k in ends)
              + (f"   [withdrawn: {api.WITHDRAWN_VARIANTS[v]}]" if v in api.WITHDRAWN_VARIANTS else ""))
    print("reported, not decided on (PREREG 7.13): Li 2025 grain on the rough-grain basis "
          + "  ".join(f"{v} {rmse(err[v]['li_grain_rough']):.3f}" for v in VARIANTS
                      if err[v].get("li_grain_rough")))
    dec = {}
    for v in VARIANTS:
        if v == DEFAULT:
            continue
        dec[v] = {k: D.decide(np.array(err[DEFAULT][k]), np.array(err[v][k])) for k in ends}
        print(f"\n{v} against {DEFAULT}:")
        for k in ends:
            d = dec[v][k]
            print(f"  {k:14} {d['verdict']:6s} delta {d['delta']:+.3f} [{d['lo']:+.3f}, {d['hi']:+.3f}]")
    # a withdrawn variant (PREREG 9.4) is scored for the record but cannot become the default
    switch = [v for v in dec if v not in api.WITHDRAWN_VARIANTS and any(dec[v][k]["verdict"] == "better" for k in TRANSFER)
              and all(dec[v][k]["verdict"] != "worse" for k in ends)]
    default = DEFAULT if not switch else min(
        switch, key=lambda v: np.median([rmse(err[v][k]) for k in TRANSFER]))
    scores = {k: rmse(err[default][k]) for k in TRANSFER}
    band = float(np.median(list(scores.values())))
    print(f"\nDEFAULT: {default}"
          + (f" (switch candidates: {switch})" if switch else " (no variant better and not worse)"))
    print("transfer errors of the default: " + "  ".join(f"{k} {s:.3f}" for k, s in scores.items())
          + f"   median (band) {band:.3f} log10 = x{10 ** band:.2f}")
    with open(OUT, "w") as f:
        json.dump(dict(errors=err, decisions=dec, default=default, switch_candidates=switch,
                       withdrawn=api.WITHDRAWN_VARIANTS,
                       default_scores=scores, band_log10=band), f, indent=1, default=float)
    print(f"[written] {os.path.relpath(OUT, ROOT)}")
    labels = {"tang": "Tang 2026 per-organ TF, 0.1 ug/g", "kim_point": "Kim 2019 grain BAF (point)",
              "kim_interval": "Kim 2019 grain BAF ([0, MDL] interval)",
              "li_shoot": "Li 2025 shoot/root TF by stage", "li_grain": "Li 2025 grain/root TF by stage"}
    with open(os.path.join(ROOT, "params", "transfer_scores.json"), "w") as f:
        json.dump(dict(source="validation/transfer_sourced.py (PREREG 5.3)", variant=default,
                       scores={labels[k]: round(v, 4) for k, v in scores.items()},
                       yamazaki_in_sample=round(rmse(err[default]["yamazaki"]), 4)), f, indent=1)
    print("[written] params/transfer_scores.json")


if __name__ == "__main__":
    main()
