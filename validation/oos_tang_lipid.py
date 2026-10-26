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
# pool, not protein, carries the longest/most-sorptive species). CRUCIALLY, the lipid
# constants are fitted on YAMAZAKI only (validation/refit_sourced.py, per shoot model), NOT
# on Tang -- so applying the mechanism to Tang is an OUT-OF-SAMPLE test of whether it
# GENERALIZES. The scored endpoint and the decision are in validation/transfer_sourced.py;
# this script prints the per-organ detail.
#
# This script contrasts the OOS Tang error with the mechanism OFF vs ON. If turning it
# on (no parameters touched for Tang) drops the OOS RMSE toward the in-sample level --
# especially by fixing the PFOS under-prediction -- the mechanism is REAL (generalizes
# across datasets), not an in-sample curve-fit.
#
# Runs on the calibration site-year drivers (Tang's site and dates are not in the repo, so
# the drivers are transferred along with the parameters) and compares dry-weight TF built
# from the harvest operator (model_api.tang_tf_validation).
#
#   python validation/oos_tang_lipid.py
# =============================================================================
import os, sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))
import model_api as api

def _rmse(dose, **kw):
    e, rows = [], {}
    for c in api.TANG_CONGENERS:
        v = api.tang_tf_validation(c, f_xy_source="recommended", dose=dose, **kw)
        for org in v["organs"]:
            m, tg = v["model_tf"][org], v["tang_tf"][org]
            e.append((np.log10(max(m, 1e-6)) - np.log10(tg)) ** 2)
            rows[(c, org)] = (m, tg)
    return float(np.sqrt(np.mean(e))), rows


def main(dose="low"):
    print("Tang 2026 OUT-OF-SAMPLE: does the lipid-loading mechanism (Yamazaki-fit, NOT")
    print("fit to Tang) generalize to predict the independent Tang per-organ TF (dw)?")
    print("Both shoot models are scored: the 4-pool core (simulate) and the redistributed")
    print("shoot (simulate_nstem_leaf).\n")
    variants = [("4pool free", dict(shoot="4pool", lipid_loading=False)),
                ("4pool lipid", dict(shoot="4pool", lipid_loading=True)),
                ("redist free", dict(shoot="redistributed", lipid_loading=False)),
                ("redist lipid", dict(shoot="redistributed", lipid_loading=True))]
    res = {name: _rmse(dose, **kw) for name, kw in variants}
    print(f"{'cong':6}{'organ':>11}" + "".join(f"{n:>14}" for n, _ in variants) + f"{'Tang':>9}")
    for key in res["4pool free"][1]:
        tg = res["4pool free"][1][key][1]
        print(f"{key[0]:6}{key[1]:>11}" + "".join(f"{res[n][1][key][0]:>14.3f}" for n, _ in variants)
              + f"{tg:>9.3f}")
    print()
    for n, _ in variants:
        tag = " (in-sample: f_xy fitted to these data)" if "refit" in n else ""
        print(f"OOS log10 RMSE  {n:14s} = {res[n][0]:.3f}{tag}")
    free, lip = res["redist free"][0], res["redist lipid"][0]
    print("=> the lipid mechanism GENERALIZES out-of-sample (recovers OOS prediction)"
          if lip < free - 0.4 else "=> the lipid mechanism does NOT improve OOS")
    return {n: r[0] for n, r in res.items()}


if __name__ == "__main__":
    main()
