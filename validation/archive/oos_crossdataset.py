#!/usr/bin/env python3
"""
Out-of-sample test on Kim 2019 brown-rice grain: do the Yamazaki-fit transport variants
PREDICT an independent dataset without refitting?
=========================================================================

The transport parameters are fit to Yamazaki 2023 only (validation/refit_site_year.py).
This script transfers three variants -- the monotone a-priori f_xy ("mono"), the
saturated per-congener fit ("fit") and the K_PL-gated lipid loading ("lipid") -- to Kim
2019's brown-rice (grain) BAF, pore-water basis, on the default drivers (the calibration
site-year; Kim's 30 fields have no dated crop record). Grain is compared on the model's
fresh basis (theta 0.14, harvest moisture, = the as-sold brown rice Kim measured).
PFOA is left out of the RMSE (its per-site check is validation/kim2019_pfoa_anchor.py);
the long chains have low rice detection frequency (validation/kim2019_censored.py scores
them with their non-detects). The Li 2025 comparison is validation/li2025_stage_resolved.py.

Run: python validation/oos_crossdataset.py
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
import model_api as api          # noqa: E402
import literature_params as lp   # noqa: E402

MODES = ("lipid", "mono", "fit")


def predict(nm, mode):
    kw = dict(lipid=dict(lipid_loading=True, shoot="4pool"),
              mono=dict(f_xy_source="recommended", **api.FOUR_POOL_FREE),
              fit=dict(f_xy_source="fit", **api.FOUR_POOL_FREE))[mode]
    r = api.simulate(nm, **kw)
    return {"root": r["baf_final"]["root"], "straw": r["straw_baf"], "grain": r["baf_final"]["grain"]}


def kim_grain():
    """(1) Kim 2019 grain BAF — out-of-sample.

    PFOA* is left out of the RMSE as recorded, but no variant is fitted to it (only the
    literature_params demo is); its per-site check is validation/kim2019_pfoa_anchor.py."""
    kim = lp.kim2019_grain_baf("porewater")
    DF = {"PFHpA": 13, "PFOA": 57, "PFNA": 20, "PFDA": 6.7, "PFUnDA": 13, "PFDoDA": 3.3}
    print("=== Kim 2019 grain BAF [porewater]  (OUT-OF-SAMPLE; * = not in the RMSE, see kim2019_pfoa_anchor) ===")
    print(f"{'PFAS':7}{'DF%':>5}{'obs':>9}" + "".join(f"{m:>9}" for m in MODES))
    rows = []
    for nm, o in kim.items():
        if nm not in api.CONGENERS:
            continue
        p = {m: predict(nm, m)["grain"] for m in MODES}
        print(f"{nm:6}{'*' if nm=='PFOA' else ' '}{DF.get(nm,0):>5.0f}{o:>9.2f}"
              + "".join(f"{p[m]:>9.2f}" for m in MODES))
        rows.append((nm, o, p, DF.get(nm, 0)))

    def rmse(mode, hiDF=False):
        e = [(np.log10(max(p[mode], 1e-6)) - np.log10(o)) ** 2
             for nm, o, p, df in rows if nm != "PFOA" and (not hiDF or df >= 15)]
        return float(np.sqrt(np.mean(e))) if e else float("nan")
    print("log10 RMSE (excl PFOA)      : " + "  ".join(f"{m}={rmse(m):.2f}" for m in MODES))
    print("log10 RMSE (DF>=15%, no PFOA): " + "  ".join(f"{m}={rmse(m, True):.2f}" for m in MODES)
          + "   [PFNA alone: PFHpA's rice DF is 13%; see validation/kim2019_censored.py]")
    return {m: (rmse(m), rmse(m, True)) for m in MODES}


if __name__ == "__main__":
    kim_grain()
