#!/usr/bin/env python3
# =============================================================================
# validation/apriori_prediction.py
# -----------------------------------------------------------------------------
# Honest OUT-OF-SAMPLE (a-priori) prediction test, and a bounded MODEL-IMPROVEMENT
# attempt, driven by the sci-adk verdict that REFUTED "RMSE 0.029 = validation".
#
# The saturated W2 fit (reproduce_demo.py) gets log10 RMSE 0.029 by fitting ~3
# transport params per congener to 3 tissue observations. The genuine predictive
# error uses the theory/QSPR MONOTONE f_xy (a-priori, NOT fit): reproduce_demo.py
# --rec gives log10 RMSE 0.837, dominated by STRAW over-prediction (6-40x) -- the
# single mass-weighted straw compartment cannot host the observed shoot gradient.
#
# IMPROVEMENT HYPOTHESIS: the redistributed-shoot model (nstem_leaf: N stem
# segments + explicit leaf + transpiration deposition/retention) should lower the
# a-priori STRAW error vs the single-straw 4-pool core, WITHOUT any per-congener
# fitting. This script runs both with the SAME monotone f_xy + drivers and reports
# the honest out-of-sample log10 RMSE for each. The result (improvement or honest
# negative) is recorded back into the sci-adk run as append-only Evidence.
#
#   python validation/apriori_prediction.py
# =============================================================================
import csv, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
import model_api as M   # noqa: E402

SEASON = 120.0

# observed Yamazaki tissue BAF
obs = {}
with open(os.path.join(ROOT, "data_obs", "obs_baf_Yamazaki.csv"), newline="") as f:
    for r in csv.DictReader(f):
        obs.setdefault(r["compound"], {})[r["tissue"]] = float(r["baf"])


def _rmse(pred, o):
    e = [(np.log10(max(pred[k], 1e-6)) - np.log10(o[k])) ** 2
         for k in ("root", "straw", "grain") if k in o]
    return e


def _single_straw(cong):
    r = M.simulate(cong, f_xy_source="recommended", season=SEASON)
    return {"root": r["baf_final"]["root"], "straw": r["straw_baf"],
            "grain": r["baf_final"]["grain"]}


def _redistributed_shoot(cong):
    r = M.simulate_nstem_leaf(cong, f_xy_source="recommended", season=SEASON)
    N = r["N"]
    Mfin = r["M"][-1]
    Mstem = float(np.sum(Mfin[1:N + 1])); Mleaf = float(Mfin[N + 1])
    bf = r["baf_final"]
    straw = (bf["stem"] * Mstem + bf["leaf"] * Mleaf) / (Mstem + Mleaf)
    return {"root": bf["root"], "straw": straw, "grain": bf["grain"]}


def main():
    print(f"{'PFAS':8}{'nC':>3} | {'root b/i/o':>22}{'straw b/i/o':>26}{'grain b/i/o':>22}")
    eb, ei = [], []
    for cong in M.CONGENERS:
        if cong not in obs:
            continue
        try:
            base = _single_straw(cong)
            impr = _redistributed_shoot(cong)
        except Exception as e:                      # keep going; record the gap
            print(f"{cong:8} ERROR {e}")
            continue
        o = obs[cong]
        eb += _rmse(base, o); ei += _rmse(impr, o)
        nC = M._CONG[cong]["n_C"]
        def cell(k):
            return (f"{base[k]:.2f}/{impr[k]:.2f}/{o.get(k, float('nan')):.2f}")
        print(f"{cong:8}{nC:>3} | {cell('root'):>22}{cell('straw'):>26}{cell('grain'):>22}")
    rb = float(np.sqrt(np.mean(eb))); ri = float(np.sqrt(np.mean(ei)))
    print(f"\na-priori (monotone f_xy) out-of-sample log10 RMSE vs Yamazaki:")
    print(f"  single-straw 4-pool   : {rb:.3f}")
    print(f"  redistributed-shoot   : {ri:.3f}   "
          f"({'IMPROVED' if ri < rb else 'no improvement'}; delta {rb - ri:+.3f})")
    return rb, ri


if __name__ == "__main__":
    main()
