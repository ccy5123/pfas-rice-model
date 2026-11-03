#!/usr/bin/env python3
# =============================================================================
# validation/yamazaki_husk_check.py  (PREREG 10.5; expectation E41)
# -----------------------------------------------------------------------------
# A record, not a scored endpoint: the default's husk / brown-rice ratio (dry; the model's husk +
# rachis pool, `conc_dw["husk"]`, over its brown rice, `conc_dw["grain"]`) against Yamazaki SI Tables
# S18/S19 (husk; brown = 0.88 white + 0.12 bran, the white share of CLAUDE.md 6.1; "<x" at x/2, the
# midpoint of [0, MQL)). Both cultivars combined by geometric mean, as the observable combines its two
# harvest runs. PFBS is not scored (its husk prints 0 / 1 and its bran and white rice are < MQL).
# Choi 2018's PFOA/PFOS ratios are printed beside it. Nothing is fitted.
#
#   python validation/yamazaki_husk_check.py
# =============================================================================
from __future__ import annotations

import csv
import json
import os
import sys
import warnings

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
warnings.filterwarnings("ignore")
import model_api as api                 # noqa: E402
import site_scenario as ss              # noqa: E402
from choi2018_husk_check import observed as choi_observed   # noqa: E402

RAW = os.path.join(ROOT, "docs", "literature_db", "raw_si", "yamazaki2023_si_s16_s18_s19.csv")
OUT = os.path.join(HERE, "yamazaki_husk_record.json")
WHITE = 0.88


def _v(x):
    return 0.5 * float(x[1:]) if x.startswith("<") else float(x)


def observed():
    rows = list(csv.DictReader(l for l in open(RAW) if not l.startswith("#")))
    out = {}
    for r in rows:
        brown = WHITE * _v(r["white_rice"]) + (1.0 - WHITE) * _v(r["bran"])
        out.setdefault(r["compound"], []).append(_v(r["husk"]) / brown)
    return {c: float(np.exp(np.mean(np.log(v)))) for c, v in out.items() if min(v) > 0}


def model_ratio(c):
    logs = []
    for h in ss.YAMAZAKI_HARVESTS:
        r = api.simulate(c, drivers=api.scenario_drivers(Cwo=1.0, harvest=h))
        logs.append(np.log10(r["conc_dw"]["husk"][-1] / r["conc_dw"]["grain"][-1]))
    return float(10 ** np.mean(logs))


def main():
    obs = observed()
    choi = choi_observed()
    res = {}
    print(f"{'':8}{'model':>9}{'Yamazaki':>10}{'log10 m/o':>11}   Choi 2018")
    for c in api.CONGENERS:
        if c not in obs or c == "PFBS":
            continue
        m = model_ratio(c)
        res[c] = dict(model=m, yamazaki=obs[c], log10_model_over_obs=float(np.log10(m / obs[c])),
                      choi=choi.get(c))
        ch = f"{choi[c][0]:.1f}-{choi[c][1]:.1f}" if c in choi else ""
        print(f"{c:8}{m:9.2f}{obs[c]:10.2f}{res[c]['log10_model_over_obs']:+11.2f}   {ch}")
    with open(OUT, "w") as f:
        json.dump(dict(source="validation/yamazaki_husk_check.py (PREREG 10.5; a record)",
                       variant=api.variant_key(), white_share=WHITE, congeners=res), f, indent=1)
    print(f"[written] {os.path.relpath(OUT, ROOT)}")


if __name__ == "__main__":
    main()
