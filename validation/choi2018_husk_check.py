#!/usr/bin/env python3
# =============================================================================
# validation/choi2018_husk_check.py
# -----------------------------------------------------------------------------
# A CHECK, not a scored endpoint (docs/PREREG_sourced_values.md 7.2, expectation E21).
#
# Choi et al. 2018 (J. Appl. Biol. Chem. 61:191, Table 1; raw_si/choi2018_husk_brown.csv) report
# PFOA and PFOS in rice husk and brown rice from pots of spiked soil. With the panicle's
# transpiration deposit kept in the husk + rachis pool and brown rice fed by the phloem only
# (7.2), the model makes a husk / brown-rice ratio. It runs on the calibration site-year drivers
# (Choi's greenhouse pots, cultivar Sindongjin, 138 d, are not modelled); the ratio is a property
# of the plant model at a linear exposure, so the exposure level does not enter. Nothing is fitted
# to Choi; the model's husk + rachis pool is compared with Choi's husk (their rachis is not
# reported separately).
#
#   python validation/choi2018_husk_check.py
# =============================================================================
from __future__ import annotations

import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
import model_api as api                 # noqa: E402

RAW = os.path.join(ROOT, "docs", "literature_db", "raw_si", "choi2018_husk_brown.csv")


def observed():
    with open(RAW) as f:
        rows = [r for r in csv.DictReader(line for line in f if not line.startswith("#"))]
    out = {}
    for r in rows:
        out.setdefault(r["compound"], []).append(float(r["husk_mg_kg"]) / float(r["brown_rice_mg_kg"]))
    return {c: (min(v), max(v)) for c, v in out.items()}


def model_ratio(c, **kw):
    r = api.simulate(c, **kw)["baf_final_dw"]
    return r["husk"] / r["grain"]


def main(**kw):
    """`kw` pass to `simulate` (the variant; None -> the default)."""
    obs = observed()
    res = {}
    print(f"{'':6} {'model husk/brown':>17} {'Choi 2018':>14} {'log10 distance':>15}")
    for c in ("PFOA", "PFOS"):
        m = model_ratio(c, **kw)
        lo, hi = obs[c]
        d = 0.0 if lo <= m <= hi else (np.log10(m / lo) if m < lo else np.log10(m / hi))
        res[c] = dict(model=m, choi=(lo, hi), log10_distance=float(d))
        print(f"{c:6} {m:17.2f} {lo:6.1f}-{hi:<6.1f} {d:+15.3f}")
    e21 = all(res[c]["model"] > 1.0 for c in res) and res["PFOS"]["model"] > res["PFOA"]["model"]
    print(f"\nE21 (husk/brown > 1 for both, PFOS > PFOA): {'HELD' if e21 else 'WRONG'}")
    return res


if __name__ == "__main__":
    main()
