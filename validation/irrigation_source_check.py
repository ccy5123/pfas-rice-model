#!/usr/bin/env python3
# =============================================================================
# validation/irrigation_source_check.py  (PREREG 9.3, 9.9; expectations E30, E31)
# -----------------------------------------------------------------------------
# The Yamazaki field's PFAS arrive with irrigation water and rain onto a soil that already holds
# them. The Tsukuba Andosol paddy (soil_hydrus, PREREG 9.2/9.8/9.10) is run per congener with
#   initial plough-layer pore water = the measured pre-season soil at equilibrium (Table S16),
#   the subsoil's (unmeasured) at 0 / half / the plough layer's (`sub`: zero / mid / same),
#   inflow = irrigation water + rain (Table S16), flux-weighted each day.
# Recorded:
#   * the root-zone pore water over the irrigation grab (season mean, day 15, end of the first
#     flooding period, end) -- E30;
#   * the OUT-OF-SAMPLE check: predicted post-season 0-15 cm total soil vs the measured one
#     (water + methanol-after-water extracts, pg/g dw), log10 ratio -- E31 (single composite
#     samples: recorded, no statistical decision);
#   * the plant: each Yamazaki cell (root, straw, grain; mean-date run, absolute dry harvest
#     concentration) at the scenario's pore water against the observable's constant grab, delta
#     log10 (enters the band under 5.4 when >= 0.0005).
# The water after midsummer drainage (intermittent irrigation, no sourced cycle; PREREG 9.13) runs at
# both ends: "supplied" (the flooded water balance; the base, every subsoil end) and "rain" (no
# irrigation after the drainage; mid subsoil).
# Nothing here is fitted.
#
#   python validation/irrigation_source_check.py           (~5 min; needs the HYDRUS engine)
# =============================================================================
from __future__ import annotations

import json
import os
import sys
import warnings

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
warnings.filterwarnings("ignore")
import model_api as api      # noqa: E402
import soil_hydrus as sh     # noqa: E402

OUT = os.path.join(HERE, "irrigation_source_record.json")
SUBS = ("mid", "zero", "same")
CASES = [(sub, "supplied") for sub in SUBS] + [("mid", "rain")]


def run_congener(c):
    grab = api.yamazaki_obs(c)["water_ng_L"] / 1000.0          # ug/L
    pre, post = sh.post_season_soil(c)
    rec = {"pre_soil": pre, "post_soil_measured": post, "sub": {}}
    base = api.simulate(c, drivers=api.scenario_drivers(Cwo=grab, harvest="mean"))
    for sub, after in CASES:
        key = sub if after == "supplied" else f"{sub}_{after}"
        drv, res = api.hydrus_drivers(c, source="irrigation", sub=sub, Cwo_ref=grab,
                                      after_drainage=after)
        ratio = res.Cw / sh._s16(c)["irrigation_ng_L"][0]
        flood1 = float(res.soil["drain_start"])
        rec["sub"][key] = dict(
            Kd=res.Kd, ks_pan=res.ks_pan,
            pore_over_irrigation=dict(
                mean=float(np.mean(ratio)), day15=float(np.interp(15.0, res.t, ratio)),
                end_flood1=float(np.interp(flood1, res.t, ratio)), end=float(ratio[-1])),
            post_soil_pred=float(res.plough_total[-1]),
            log10_post_pred_over_obs=float(np.log10(res.plough_total[-1] / post)),
            log10_post_pred_over_pre_pred=float(np.log10(res.plough_total[-1]
                                                         / res.plough_total[0])))
        if sub == "mid":
            # absolute harvest concentrations (simulate's BAF is per its own final pore water, so
            # the two runs' BAFs have different denominators); the observable divides both by the grab
            scen = api.simulate(c, drivers=drv)
            rec["plant_delta_log10" if after == "supplied" else f"plant_delta_log10_{after}"] = {
                k: float(np.log10(max(scen["conc_dw"][k][-1], 1e-30)
                                  / max(base["conc_dw"][k][-1], 1e-30)))
                for k in ("root", "straw", "grain")}
    rec["log10_post_obs_over_pre_obs"] = float(np.log10(post / pre))
    return rec


def main():
    if not sh.hydrus_available():
        print("HYDRUS-1D engine not built; nothing to run.")
        return
    congs = list(api.yamazaki_obs().keys())
    out = {}
    print(f"{'':8}{'Kd':>8} {'pore/irr: mean d15 flood1 end':>32} | "
          f"{'post soil log10 pred/obs (mid zero same | mid rain-only)':>52} | obs post/pre | "
          f"plant dlog root straw grain (| rain-only)")
    for c in congs:
        r = run_congener(c)
        out[c] = r
        m = r["sub"]["mid"]
        p = m["pore_over_irrigation"]
        print(f"{c:8}{m['Kd']:8.3f} {p['mean']:7.2f}{p['day15']:7.2f}{p['end_flood1']:8.2f}"
              f"{p['end']:7.2f}   | "
              + "  ".join(f"{r['sub'][s]['log10_post_pred_over_obs']:+.2f}" for s in SUBS)
              + f" | {r['sub']['mid_rain']['log10_post_pred_over_obs']:+.2f}"
              + f"{'':14} | {r['log10_post_obs_over_pre_obs']:+.2f}        | "
              + "  ".join(f"{r['plant_delta_log10'][k]:+.3f}" for k in ("root", "straw", "grain"))
              + " | " + "  ".join(f"{r['plant_delta_log10_rain'][k]:+.3f}"
                                  for k in ("root", "straw", "grain")))
    with open(OUT, "w") as f:
        json.dump(dict(source="validation/irrigation_source_check.py (PREREG 9.3; nothing fitted)",
                       soil=sh.andosol_soil(), congeners=out), f, indent=1)
    print(f"[written] {os.path.relpath(OUT, ROOT)}")


if __name__ == "__main__":
    main()
