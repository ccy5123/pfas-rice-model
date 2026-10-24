#!/usr/bin/env python3
# =============================================================================
# validation/kim2019_censored.py
# -----------------------------------------------------------------------------
# KIM 2019 GRAIN BAF WITH ITS NON-DETECTS  (handoff step 4, docs/HANDOFF_mass_basis_refit.md)
#
# PRE-REGISTRATION -- written and pushed before any model number below was produced.
#
# Why: every Kim 2019 grain BAF in the repo is a ratio of two field averages in which the
# non-detects were averaged as ZERO (stated in the paper). For the low-detection
# congeners both the numerator (brown rice, n=30) and the denominator (void water, n=27)
# are mostly non-detects -- PFDoDA rice 1/30, water 3/27 -- so the "measured" BAF is one
# corner of a range set by how the non-detects are treated, not a measurement.
#
# Data: SI Tables S5 (MDL) + S6 (min/avg/max/DF) -> docs/literature_db/raw_si/
# kim2019_si_detection.csv (averages identical to kim2019_field_conc.csv). Each non-detect
# lies in [0, MDL], so the true ratio of means lies in
#     [rice_avg / water_hi,  rice_hi / water_avg],   x_hi = x_avg + n_ND * MDL / n.
# Already computed (data only, before this file): PFHpA [0.39, 2.11], PFOA [4.43, 4.92],
# PFNA [2.21, 7.07], PFDA [1.55, 11.3], PFUnDA [19.8, 168], PFDoDA [5.5, 630] L/kg. The
# shipped ND=0 values are the LOWER bound for PFHpA/PFNA/PFDA/PFOA. Also already seen:
# the "reliable DF>=15%" Kim score in oos_crossdataset is ONE congener (PFNA) -- the
# label there says PFHpA+PFNA, but PFHpA's DF is 13%.
#
# Model: grain BAF at harvest, Cwo constant 1, default ORYZA2000 biomass + measured Q_TP,
# fresh organ mass (as oos_crossdataset). Variants: 4-pool monotone f_xy, 4-pool W2,
# 4-pool lipid loading (the three oos_crossdataset scores), plus nstem_leaf + lipid (the
# configuration that won Tang and Li). Scored congeners: PFHpA, PFNA, PFDA, PFUnDA, PFDoDA
# (PFOA stays out, as in oos_crossdataset). Model fresh-weight BAF vs Kim as published
# (brown-rice basis not stated; dw would move Kim by 0.065 log).
#
# Questions and decision rules:
#   Q1 (data) Is the long-chain RISE in Kim's grain BAF robust to the non-detects? A step
#      counts as robust when the lower bound of the longer chain exceeds the upper bound
#      of the shorter one. (Seen already: PFUnDA lo 19.8 > PFDA hi 11.3 and PFNA hi 7.07
#      -> the C11 rise is robust; PFDoDA lo 5.5 < PFDA hi 11.3 -> the C12 level is not.)
#   Q2 (model) log10 RMSE per variant two ways: POINT (vs the ND=0 value, the existing
#      score) and INTERVAL (error 0 inside the bounds, else distance to the nearer bound).
#   VERDICT RULE: a variant's Kim advantage is "robust to the non-detects" only if it beats
#      every other variant by more than 0.10 log10 on BOTH the point and the interval score.
#   EXPECTATION (stated so it can be wrong): lipid's advantage over monotone and W2
#      survives the interval score (monotone sits below even the lower bounds of the long
#      chains); every RMSE shrinks under the interval score; no prediction for nstem_leaf
#      + lipid vs 4-pool lipid.
#
#   python validation/kim2019_censored.py
# =============================================================================
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, HERE)

import literature_params as lp     # noqa: E402
import model_api as api            # noqa: E402

SCORED = ("PFHpA", "PFNA", "PFDA", "PFUnDA", "PFDoDA")        # PFOA out, as in oos_crossdataset
CHAIN = ("PFHpA", "PFOA", "PFNA", "PFDA", "PFUnDA", "PFDoDA")
VARIANTS = ("mono", "W2", "lipid", "nstem_lipid")
MARGIN = 0.10


def predict(nm, v):
    """Harvest grain BAF [L/kg fw], Cwo = 1, default ORYZA2000 biomass + measured Q_TP."""
    if v == "nstem_lipid":
        r = api.simulate_nstem_leaf(nm, Cwo=1.0, season=120.0, lipid_loading=True,
                                    biomass_fn=api._biomass_fn("oryza"))
        return float(r["baf_final"]["grain"])
    kw = dict(mono=dict(f_xy_source="recommended"), W2=dict(f_xy_source="W2fit"),
              lipid=dict(lipid_loading=True))[v]
    return float(api.simulate(nm, **kw)["baf_final"]["grain"])


def point_error(pred, bounds):
    return np.log10(max(pred, 1e-12)) - np.log10(bounds[1])


def interval_error(pred, bounds):
    lo, _, hi = bounds
    lp_ = np.log10(max(pred, 1e-12))
    if lp_ < np.log10(lo):
        return lp_ - np.log10(lo)
    if lp_ > np.log10(hi):
        return lp_ - np.log10(hi)
    return 0.0


def rise_steps(b):
    """Q1: each step along the chain is robust when lo(longer) > hi(shorter)."""
    out = {}
    for a, z in (("PFNA", "PFUnDA"), ("PFDA", "PFUnDA"), ("PFDA", "PFDoDA"), ("PFUnDA", "PFDoDA")):
        out[(a, z)] = bool(b[z][0] > b[a][2])
    return out


def verdict(scores):
    for v in VARIANTS:
        if all(scores[u][k] - scores[v][k] > MARGIN
               for u in VARIANTS if u != v for k in ("point", "interval")):
            return v
    return None


def main():
    b = lp.kim2019_grain_baf_bounds()
    print("=" * 84)
    print("KIM 2019 GRAIN BAF WITH ITS NON-DETECTS — pre-registered (see header)")
    print("=" * 84)
    print(f"\n  {'PFAS':7}{'DF rice':>8}{'DF water':>9}{'lower':>9}{'published':>11}{'upper':>9}{'width':>8}")
    for nm in CHAIN:
        df_r, df_w = lp.KIM2019_FIELD[nm][3], lp.KIM2019_DETECTION[nm][2]
        lo, pt, hi = b[nm]
        print(f"  {nm:7}{df_r:>7.1f}%{df_w:>8.0f}%{lo:>9.2f}{pt:>11.2f}{hi:>9.2f}{hi / lo:>7.0f}x")

    steps = rise_steps(b)
    print("\nQ1 — is the long-chain rise robust to the non-detects? (lo of the longer > hi of the shorter)")
    for (a, z), ok in steps.items():
        print(f"  {a:>6} -> {z:<7} {'ROBUST' if ok else 'not robust'}   ({b[z][0]:.2f} vs {b[a][2]:.2f})")

    preds = {v: {nm: predict(nm, v) for nm in CHAIN} for v in VARIANTS}
    print("\nQ2 — model grain BAF (L/kg) against the interval:")
    print(f"  {'PFAS':7}{'interval':>17}  " + "".join(f"{v:>12}" for v in VARIANTS))
    for nm in CHAIN:
        lo, _, hi = b[nm]
        tag = " (anchor, not scored)" if nm not in SCORED else ""
        print(f"  {nm:7}{f'[{lo:.2f}, {hi:.2f}]':>17}  "
              + "".join(f"{preds[v][nm]:>12.3g}" for v in VARIANTS) + tag)
    scores = {}
    for v in VARIANTS:
        ep = [point_error(preds[v][nm], b[nm]) for nm in SCORED]
        ei = [interval_error(preds[v][nm], b[nm]) for nm in SCORED]
        scores[v] = dict(point=float(np.sqrt(np.mean(np.square(ep)))),
                         interval=float(np.sqrt(np.mean(np.square(ei)))),
                         inside=int(sum(e == 0.0 for e in ei)))
    print(f"\n  {'variant':12}{'point RMSE':>12}{'interval RMSE':>15}{'inside':>8}")
    for v in VARIANTS:
        s = scores[v]
        print(f"  {v:12}{s['point']:>12.2f}{s['interval']:>15.2f}{s['inside']:>6}/{len(SCORED)}")

    win = verdict(scores)
    print("\n" + "=" * 84)
    print("VERDICT (against the pre-registered rule)")
    print("=" * 84)
    print(f"  robust winner (> {MARGIN} on point AND interval): {win or 'none'}")
    print(f"  long-chain rise: C11 {'robust' if steps[('PFDA', 'PFUnDA')] else 'not robust'}, "
          f"C12 {'robust' if steps[('PFDA', 'PFDoDA')] else 'not robust'}")
    print("  the 'reliable DF>=15%' Kim score in oos_crossdataset is PFNA alone (PFHpA DF 13%)")
    return dict(bounds=b, steps=steps, preds=preds, scores=scores, winner=win)


if __name__ == "__main__":
    main()
