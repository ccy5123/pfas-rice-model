#!/usr/bin/env python3
# =============================================================================
# validation/oos_multidataset.py
# -----------------------------------------------------------------------------
# MULTI-DATASET out-of-sample report: the Yamazaki-fit transport variants transferred
# WITHOUT refitting (and on the calibration site-year drivers) to three independent
# datasets:
#   * Tang 2026  per-organ TF (dw; PFOA/PFOS/GenX, 0.1 ug/g)    validation/oos_tang_lipid.py
#   * Kim 2019   brown-rice grain BAF (pore water; excl PFOA)   validation/oos_crossdataset.py
#   * Li 2025    stage-resolved TF (dw)                         validation/li2025_stage_resolved.py
# Run: python validation/oos_multidataset.py
# =============================================================================
import contextlib
import io
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, HERE)

import oos_tang_lipid as tang          # noqa: E402
import oos_crossdataset as cross       # noqa: E402
import li2025_stage_resolved as li     # noqa: E402


def main(verbose=True):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf if not verbose else sys.stdout):
        free_tang, lip_tang = tang.main(dose="low")
        print()
        kim = cross.kim_grain()
        print()
        lio = li.main(sensitivity=False)
    out = dict(tang=dict(free=free_tang, lipid=lip_tang), kim=kim, li=lio)
    if verbose:
        print("\n" + "=" * 72)
        print("SUMMARY  (OOS log10 RMSE; lower = better)")
        print("=" * 72)
        print(f"{'Tang 2026 per-organ TF (dw)':40} free anion {free_tang:.2f}   lipid {lip_tang:.2f}")
        print(f"{'Kim 2019 grain, excl PFOA':40} " + "   ".join(f"{m} {kim[m][0]:.2f}" for m in cross.MODES))
        sc = lio["scores"]
        print(f"{'Li 2025 stage TF, shoot / grain':40} " + "   ".join(
            f"{v} {sc[v]['rmse']['shoot']:.2f}/{sc[v]['rmse']['grain']:.2f}" for v in li.VARIANTS))
        print(f"{'Li 2025 pre-registered verdict':40} {lio['winner'] or 'inconclusive'}")
    return out


if __name__ == "__main__":
    main()
