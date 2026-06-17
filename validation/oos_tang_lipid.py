#!/usr/bin/env python3
# =============================================================================
# validation/oos_tang_lipid.py
# -----------------------------------------------------------------------------
# DOES THE LIPID-FACILITATED LOADING MECHANISM GENERALIZE OUT-OF-SAMPLE?
#
# `validation/oos_tang.py` showed the FREE-ANION model (theory monotone f_xy, NOT
# fit to Tang) FAILS to predict the independent Tang 2026 per-organ TF out-of-sample
# (log10 RMSE 1.232 vs in-sample Tang-refit 0.519) -- the dominant error is PFOS,
# under-predicted ~40-200x (the high-K_PL sulfonate the free-anion route starves).
#
# The long-chain investigation (runs/pfas-rice-longchain, LC1/LC2 SUPPORTED) found a
# B-INDEPENDENT lipid-facilitated bound-loading term (g_xy*C, g_ph*C, K_PL-gated) that
# fixes the IN-SAMPLE long-chain (Yamazaki C10-C12) collapse, and it is independently
# corroborated by Chen 2025 (membrane K_MW rises monotone with chain length; the lipid
# pool, not protein, carries the longest/most-sorptive species). CRUCIALLY, the
# LIPID_LOADING constants were fit on YAMAZAKI (excl. PFDoDA), NOT on Tang -- so applying
# the mechanism to Tang is a genuine OUT-OF-SAMPLE test of whether it GENERALIZES.
#
# This script contrasts the OOS Tang error with the mechanism OFF vs ON. If turning it
# on (no parameters touched for Tang) drops the OOS RMSE toward the in-sample level --
# especially by fixing the PFOS under-prediction -- the mechanism is REAL (generalizes
# across datasets), not an in-sample curve-fit.
#
#   python validation/oos_tang_lipid.py
# =============================================================================
import os, sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))
import model_api as api


def main(dose="low"):
    print("Tang 2026 OUT-OF-SAMPLE: does the lipid-loading mechanism (Yamazaki-fit, NOT")
    print("fit to Tang) generalize to predict the independent Tang per-organ TF (dw)?\n")
    print(f"{'cong':6}{'organ':>11}{'free':>9}{'lipid':>9}{'Tang':>9}")
    e_free, e_lip = [], []
    for c in api.TANG_CONGENERS:
        v0 = api.tang_tf_validation(c, f_xy_source="recommended", dose=dose, lipid_loading=False)
        v1 = api.tang_tf_validation(c, f_xy_source="recommended", dose=dose, lipid_loading=True)
        for org in v0["organs"]:
            f, l, tg = v0["model_tf"][org], v1["model_tf"][org], v0["tang_tf"][org]
            e_free.append((np.log10(max(f, 1e-6)) - np.log10(tg)) ** 2)
            e_lip.append((np.log10(max(l, 1e-6)) - np.log10(tg)) ** 2)
            print(f"{c:6}{org:>11}{f:>9.3f}{l:>9.3f}{tg:>9.3f}")
    free = float(np.sqrt(np.mean(e_free)))
    lip = float(np.sqrt(np.mean(e_lip)))
    print(f"\nOOS log10 RMSE  free-anion (current baseline) = {free:.3f}")
    print(f"OOS log10 RMSE  lipid mechanism (Yamazaki-fit) = {lip:.3f}")
    print(f"in-sample reference (Tang-refit f_xy)          ~ 0.519")
    print("=> the lipid mechanism GENERALIZES out-of-sample (recovers OOS prediction)"
          if lip < free - 0.4 else "=> the lipid mechanism does NOT improve OOS")
    return free, lip


if __name__ == "__main__":
    main()
