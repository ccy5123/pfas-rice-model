#!/usr/bin/env python3
# =============================================================================
# validation/oos_tang.py
# -----------------------------------------------------------------------------
# OUT-OF-SAMPLE cross-dataset prediction test: do the model's transport parameters
# (theory/QSPR monotone f_xy -- NOT fit to Tang) predict the INDEPENDENT Tang 2026
# per-organ transfer factors (stalk/leaf/endosperm, dry weight)? Tang is a different
# soil (flooded paddy pot), cultivar (Nipponbare), and dose set than Yamazaki, so it
# is a genuine out-of-sample test of predictive validity -- the project's central
# claim. Contrasts the OOS error with the in-sample Tang-refit f_xy.
#
#   python validation/oos_tang.py
# =============================================================================
import os, sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))
import model_api as api


def main(dose="low"):
    print("Tang 2026 OUT-OF-SAMPLE: model with theory f_xy (NOT fit to Tang) vs Tang TF (dw)\n")
    print(f"{'cong':6}{'organ':>11}{'model':>9}{'Tang':>9}{'refit(IS)':>11}")
    e_oos, e_is = [], []
    for c in api.TANG_CONGENERS:
        v = api.tang_tf_validation(c, f_xy_source="recommended", dose=dose)        # OOS
        vr = api.tang_tf_validation(c, f_xy_source="recommended", use_refit=True, dose=dose)  # in-sample refit
        for org in v["organs"]:
            m, tg, rf = v["model_tf"][org], v["tang_tf"][org], vr["model_tf"][org]
            e_oos.append((np.log10(max(m, 1e-6)) - np.log10(tg)) ** 2)
            e_is.append((np.log10(max(rf, 1e-6)) - np.log10(tg)) ** 2)
            print(f"{c:6}{org:>11}{m:>9.3f}{tg:>9.3f}{rf:>11.3f}")
    oos = float(np.sqrt(np.mean(e_oos))); insamp = float(np.sqrt(np.mean(e_is)))
    print(f"\nOOS log10 RMSE (theory params, NOT fit to Tang) = {oos:.3f}")
    print(f"in-sample (Tang-refit f_xy)                     = {insamp:.3f}")
    print("=> the model does NOT predict the independent Tang dataset out-of-sample"
          if oos > insamp + 0.4 else "=> the model transfers to Tang")
    return oos, insamp


if __name__ == "__main__":
    main()
