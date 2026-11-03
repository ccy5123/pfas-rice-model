#!/usr/bin/env python3
# =============================================================================
# validation/site_kd_check.py  (PREREG 10.1-10.2; expectations E36-E38)
# -----------------------------------------------------------------------------
# The irrigation-source scenario of the Yamazaki field (PREREG 9.3, unchanged: sub "mid", water
# "supplied") run on the site soil's own desorption Kd (soil_hydrus.site_desorption_kd: the
# PRE-season sequential extraction of SI Table S16; pellet water V_r at its midpoint and both ends).
# Out-of-sample check: predicted post-season 0-15 cm soil vs the measured one (S16, never used to set
# anything). Rule (10.1, per congener): within 0.3 log10 at the V_r midpoint -> admitted.
# Cases per congener:
#   mid / vr0 / vrmax  site Kd, Morita plot C roots (10.2)
#   mid_linear30       site Kd (V_r mid) on the old unsourced root profile
#   qspr               the Koc-QSPR Kd on Morita C roots (separates the root change from the Kd change)
#   corner             the most retentive sourced corner: Kd at V_r = 0, percolation 0.41 cm/d
# For admitted congeners the plant cells (root, straw, grain; mean-date, absolute dry harvest
# concentration) at the scenario's pore water vs the constant grab, delta log10, at mid and both ends.
#
#   python validation/site_kd_check.py        (~15 min; needs the HYDRUS engine)
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

OUT = os.path.join(HERE, "site_kd_record.json")
TOL = 0.3                    # log10, PREREG 10.1


def _plant_delta(c, drv, base):
    scen = api.simulate(c, drivers=drv)
    return {k: float(np.log10(max(scen["conc_dw"][k][-1], 1e-30)
                              / max(base["conc_dw"][k][-1], 1e-30)))
            for k in ("root", "straw", "grain")}


def run_congener(c, vr_max):
    grab = api.yamazaki_obs(c)["water_ng_L"] / 1000.0          # ug/L
    pre, post = sh.post_season_soil(c)
    c_irr = sh._s16(c)["irrigation_ng_L"][0]
    kd = {"mid": sh.site_desorption_kd(c, 0.5 * vr_max), "vr0": sh.site_desorption_kd(c, 0.0),
          "vrmax": sh.site_desorption_kd(c, vr_max)}
    n_C, group = api._CONG[c]["n_C"], api._CONG[c]["group"]
    cases = {
        "mid": dict(Kd=kd["mid"]),
        "vr0": dict(Kd=kd["vr0"]),
        "vrmax": dict(Kd=kd["vrmax"]),
        "mid_linear30": dict(Kd=kd["mid"], root_profile="linear30"),
        "qspr": dict(Kd=sh.paddy_kd(n_C, group)),
        "corner": dict(Kd=kd["vr0"], soil=sh.andosol_soil(percolation=0.41)),
    }
    rec = {"pre_soil": pre, "post_soil_measured": post, "Kd_site": kd, "cases": {}}
    drvs = {}
    for key, kw in cases.items():
        drv, res = api.hydrus_drivers(c, source="irrigation", sub="mid", Cwo_ref=grab, **kw)
        drvs[key] = drv
        ratio = res.Cw / c_irr
        rec["cases"][key] = dict(
            Kd=float(res.Kd),
            pore_over_irrigation=dict(mean=float(np.mean(ratio)), end=float(ratio[-1])),
            post_soil_pred=float(res.plough_total[-1]),
            log10_post_pred_over_obs=float(np.log10(res.plough_total[-1] / post)))
    rec["admitted"] = bool(abs(rec["cases"]["mid"]["log10_post_pred_over_obs"]) <= TOL)
    if rec["admitted"]:
        base = api.simulate(c, drivers=api.scenario_drivers(Cwo=grab, harvest="mean"))
        rec["plant_delta_log10"] = {k: _plant_delta(c, drvs[k], base)
                                    for k in ("mid", "vr0", "vrmax")}
    return rec


def main():
    if not sh.hydrus_available():
        print("HYDRUS-1D engine not built; nothing to run.")
        return
    vr_max = sh.extraction_vr_max()
    congs = list(api.yamazaki_obs().keys())
    out = {}
    hdr = ("mid", "vr0", "vrmax", "mid_linear30", "qspr", "corner")
    print(f"V_r max {vr_max:.4f} mL/g; post-season soil log10 pred/obs per case")
    print(f"{'':8}{'Kd_mid':>9} " + " ".join(f"{h:>12}" for h in hdr) + "  admitted")
    for c in congs:
        r = run_congener(c, vr_max)
        out[c] = r
        print(f"{c:8}{r['Kd_site']['mid']:9.3f} "
              + " ".join(f"{r['cases'][h]['log10_post_pred_over_obs']:+12.2f}" for h in hdr)
              + f"  {r['admitted']}"
              + (("  plant mid " + " ".join(f"{v:+.3f}" for v in r["plant_delta_log10"]["mid"].values()))
                 if r["admitted"] else ""), flush=True)
    with open(OUT, "w") as f:
        json.dump(dict(source="validation/site_kd_check.py (PREREG 10.1-10.2; nothing fitted)",
                       vr_max=vr_max, tol=TOL, soil=sh.andosol_soil(),
                       root_profile=sh.ROOT_PROFILE_DEFAULT, congeners=out), f, indent=1)
    print(f"[written] {os.path.relpath(OUT, ROOT)}")


if __name__ == "__main__":
    main()
