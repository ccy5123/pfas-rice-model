#!/usr/bin/env python3
# =============================================================================
# validation/li_grain_from_check.py  (PREREG 9.4, `tr_li_grain_from`)
# -----------------------------------------------------------------------------
# Li 2025's grain is scored from FILLING on (li2025_stage_resolved.GRAIN_FROM): at heading the
# brown rice has barely formed. That start has no source; this check scores the grain from
# HEADING as well and records both, per scored variant (the withdrawn *_free are not scored,
# PREREG 9.4). Brown rice (the decision basis, 7.13) and rough grain. Nothing fitted.
#
#   python validation/li_grain_from_check.py    (~3 min; writes validation/li_grain_from_record.json)
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
sys.path.insert(0, HERE)
warnings.filterwarnings("ignore")
import model_api as api                 # noqa: E402
import li2025_stage_resolved as L       # noqa: E402
import transfer_sourced as T            # noqa: E402

OUT = os.path.join(HERE, "li_grain_from_record.json")
STARTS = ("heading", "filling")


def grain_errors(kw, start, grain):
    data, days = L.load_data(), L.stage_days()
    i0 = L.STAGES.index(start)
    out = []
    for c in L.CONGENERS:
        if (c, L.GRAIN) not in data:
            continue
        res = api.simulate(c, Cwo=1.0, **kw)
        if grain == "rough" and "rough" not in res["conc_dw"]:
            return []
        tf_m = L.model_tf(res, days, grain)[L.GRAIN]
        tf_o, _ = L.observed_tf(data, c, L.GRAIN)
        idx = [i for i in range(i0, 5) if not np.isnan(tf_o[i])]
        if len(idx) < 2:
            continue
        out.append(float(np.log10(np.maximum(tf_m[idx], 1e-12)).mean()
                         - np.log10(tf_o[idx]).mean()))
    return out


def main():
    rec = {}
    for v, kw in T.VARIANTS.items():
        if v in api.WITHDRAWN_VARIANTS:
            continue
        rec[v] = {}
        for grain in ("grain", "rough"):
            for s in STARTS:
                e = grain_errors(kw, s, grain)
                rec[v][f"{grain}_from_{s}"] = dict(errors=e, rmse=T.rmse(e) if e else None)
        r = rec[v]
        print(f"{v:22} brown: heading {r['grain_from_heading']['rmse']:.3f}  filling "
              f"{r['grain_from_filling']['rmse']:.3f}   rough: heading "
              f"{r['rough_from_heading']['rmse']:.3f}  filling {r['rough_from_filling']['rmse']:.3f}")
    with open(OUT, "w") as f:
        json.dump(dict(source="validation/li_grain_from_check.py (PREREG 9.4; nothing fitted)",
                       variants=rec), f, indent=1)
    print(f"[written] {os.path.relpath(OUT, ROOT)}")


if __name__ == "__main__":
    main()
