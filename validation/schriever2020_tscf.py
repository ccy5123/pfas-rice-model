#!/usr/bin/env python3
# =============================================================================
# validation/schriever2020_tscf.py
# -----------------------------------------------------------------------------
# Tests the TSCF QSPR **on its own**, against the 97 measured root-to-shoot
# transfer values in the SI of Schriever & Lamshoeft 2020
# (data_obs/tscf_obs_schriever2020.csv).
#
# WHY THIS IS WORTH A SCRIPT. The neutral path has exactly two a-priori inputs:
# the partition K_PW and the xylem loading factor f_xy = TSCF(log Kow). K_PW has
# been tested directly twice (Liu 2023, and now Li 2019's hydroponic RCF table).
# TSCF never has -- it could only be reached THROUGH the plant ODE, in the Ge
# 2017 per-organ comparison, where its error is inseparable from the unknown
# in-planta half-life. That is precisely why docs/neutral_dpu_validation.md
# records the Ge number as an UPPER BOUND on transport error rather than a
# measurement of it. This script removes the model from the middle: measured
# TSCF in, QSPR prediction out, nothing else.
#
# READ THE IN-SAMPLE WARNING BEFORE QUOTING ANY NUMBER HERE.
#   * `schriever_tscf` was FITTED to these 97 rows. Its score is reproduction.
#     It is reported only as the reference a fitted model achieves on its own
#     training set -- the floor the Briggs bell is being measured against.
#   * `briggs_tscf` was fitted to Briggs, Bromilow & Evans 1982; the barley rows
#     here are Briggs et al. 1987, Inoue 1998 and Shone & Wood 1972/74. Different
#     compounds, so out-of-sample -- but same lab lineage, same species, same
#     method. Out-of-sample, NOT independent. Say it that way.
#
# METRIC. TSCF is a bounded fraction on [0,1] and several measured values are 0,
# so a log10 RMSE (the repo's usual metric) is undefined or dominated by the
# smallest values. Linear RMSE, bias and Spearman rank correlation are reported
# instead; rank correlation is the one that answers "does the QSPR order
# compounds correctly", which is what the plant model actually consumes.
#
#   python validation/schriever2020_tscf.py
# =============================================================================
from __future__ import annotations
import csv, os, sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT_DIR, "src"))

import neutral_dpu as ND                                      # noqa: E402

OBS = os.path.join(ROOT_DIR, "data_obs", "tscf_obs_schriever2020.csv")
MISSING = 999.0        # the SI's sentinel for "not applicable / not reported"


def load(path=OBS):
    """Read the transcribed Table A 3. Returns dicts with floats parsed."""
    rows = []
    with open(path, newline="") as f:
        for r in csv.DictReader(x for x in f if not x.lstrip().startswith("#")):
            if not r.get("id"):
                continue
            rows.append(dict(
                id=r["id"], species=r["species"], authors=r["authors"],
                pKa=float(r["pKa"]), logP=float(r["logP"]),
                logD55=float(r["logD_pH5.5"]), logD74=float(r["logD_pH7.4"]),
                logD_test=float(r["logD_test_pH"]), pH_test=float(r["pH_test"]),
                TSCF=float(r["TSCF"]),
                neutral=bool(int(r["neutral_at_test_pH"]))))
    return rows


def descriptor(r, which):
    """Lipophilicity descriptor for a row. 'logP', or 'logD' at the test pH with
    a documented fallback to pH 5.5 when the source did not report a test pH."""
    if which == "logP":
        return r["logP"]
    if which == "logD":
        return r["logD_test"] if r["pH_test"] != MISSING else r["logD55"]
    raise ValueError(which)


def spearman(a, b):
    """Rank correlation without a SciPy dependency (ties averaged)."""
    def rank(x):
        order = np.argsort(np.asarray(x, float), kind="mergesort")
        rr = np.empty(len(x), float)
        rr[order] = np.arange(len(x), dtype=float)
        # average tied ranks so a column of repeated TSCFs cannot fake a signal
        vals = np.asarray(x, float)[order]
        i = 0
        while i < len(vals):
            j = i
            while j + 1 < len(vals) and vals[j + 1] == vals[i]:
                j += 1
            if j > i:
                rr[order[i:j + 1]] = np.mean(np.arange(i, j + 1))
            i = j + 1
        return rr
    ra, rb = rank(a), rank(b)
    if np.std(ra) == 0 or np.std(rb) == 0:
        return float("nan")
    return float(np.corrcoef(ra, rb)[0, 1])


def score(rows, model, which):
    """Predicted vs measured TSCF for one QSPR on one descriptor."""
    obs = np.array([r["TSCF"] for r in rows], float)
    pred = np.array([ND.tscf(descriptor(r, which), model=model) for r in rows], float)
    return dict(n=len(rows), rmse=float(np.sqrt(np.mean((pred - obs) ** 2))),
                mae=float(np.mean(np.abs(pred - obs))),
                bias=float(np.mean(pred - obs)), rho=spearman(pred, obs),
                pred=pred, obs=obs)


def _line(tag, s):
    print(f"   {tag:34s} n={s['n']:3d}  RMSE {s['rmse']:.3f}   MAE {s['mae']:.3f}"
          f"   bias {s['bias']:+.3f}   Spearman {s['rho']:+.3f}")


def main():
    rows = load()
    neutral = [r for r in rows if r["neutral"]]
    barley = [r for r in neutral if r["species"] == "Barley"]
    other = [r for r in neutral if r["species"] != "Barley"]

    print("=" * 84)
    print("TSCF QSPR vs 97 MEASURED VALUES — Schriever & Lamshoeft 2020 SI, Table A 3")
    print("=" * 84)
    print(f"   rows {len(rows)};  neutral at the test pH {len(neutral)};  "
          f"species {len({r['species'] for r in rows})}")
    print("   TSCF is a fraction on [0,1] with zeros present, so the metric is LINEAR")
    print("   RMSE (a log10 RMSE would be undefined), plus rank correlation.\n")

    print("1. THE REPO DEFAULT (briggs_tscf) ON THE NEUTRAL SUBSET")
    print("   Out-of-sample for the 1982 bell (these are Briggs 1987 / Inoue /")
    print("   Shone & Wood compounds) but NOT independent: same lab, same barley.")
    b_logp = score(neutral, "briggs", "logP")
    b_logd = score(neutral, "briggs", "logD")
    _line("briggs, descriptor = logP", b_logp)
    _line("briggs, descriptor = logD(test pH)", b_logd)

    print("\n2. THE FITTED REFERENCE (schriever_tscf) — IN-SAMPLE, NOT A PREDICTION")
    print("   These 97 rows ARE its training set. This is the floor, not a rival.")
    s_logd = score(neutral, "schriever", "logD")
    s_logp = score(neutral, "schriever", "logP")
    _line("schriever, descriptor = logD (as fitted)", s_logd)
    _line("schriever, descriptor = logP", s_logp)

    print("\n3. WHERE THE BRIGGS BELL HOLDS AND WHERE IT DOES NOT")
    _line("briggs / barley only", score(barley, "briggs", "logP"))
    _line("briggs / every other species", score(other, "briggs", "logP"))
    _line("schriever / barley only", score(barley, "schriever", "logD"))
    _line("schriever / every other species", score(other, "schriever", "logD"))

    print("\n4. THE FULL 97, INCLUDING IONISABLE COMPOUNDS")
    print("   Outside this model's stated scope (it assumes an un-ionised solute);")
    print("   shown to size what the neutral restriction is actually buying.")
    _line("briggs, logP, all 97", score(rows, "briggs", "logP"))
    _line("briggs, logD, all 97", score(rows, "briggs", "logD"))
    _line("schriever, logD, all 97", score(rows, "schriever", "logD"))

    print("\n" + "=" * 84)
    print("VERDICT")
    print("=" * 84)
    print(f"   On the {len(neutral)} un-ionised rows the repo's default Briggs bell reaches")
    print(f"   RMSE {b_logp['rmse']:.3f} (bias {b_logp['bias']:+.3f}) against a FITTED model's")
    print(f"   in-sample {s_logd['rmse']:.3f} on the same rows. The gap between those two is the")
    print("   honest size of the TSCF term's a-priori error -- measured directly, for the")
    print("   first time in this repo, with no plant model in between.")
    print("   Carry forward: TSCF is measured over WHOLE plants in hydroponics on 16")
    print("   species, none of them rice, and the barley rows are half the table.")
    return b_logp, s_logd


if __name__ == "__main__":
    main()
