#!/usr/bin/env python3
# =============================================================================
# validation/root_profile_check.py  (PREREG 10.2; expectation E39)
# -----------------------------------------------------------------------------
# The stored flooded shapes (cwo_shape_table: contaminated soil, clean water, Koc-QSPR Kd, Yamazaki
# site-year) on each rice root profile: the old unsourced "linear30" and Morita et al. 1988 plots
# A-F (C = the run). Recorded per congener: season-mean c/c0 and its log10 change from linear30;
# the envelope over A-F is the 5.4 band of the root profile.
#
#   python validation/root_profile_check.py       (~6 min; needs the HYDRUS engine)
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
warnings.filterwarnings("ignore")
import model_api as api      # noqa: E402
import soil_hydrus as sh     # noqa: E402

OUT = os.path.join(HERE, "root_profile_record.json")
PROFILES = ["linear30"] + [f"morita_{k}" for k in sh.ROOT_RLD_MORITA1988]


def main():
    if not sh.hydrus_available():
        print("HYDRUS-1D engine not built; nothing to run.")
        return
    out = {}
    print(f"{'':8}" + "".join(f"{p:>11}" for p in PROFILES) + "   (season-mean c/c0; then log10 vs linear30)")
    for cong in api.CONGENERS:
        c = api._CONG[cong]
        means = {}
        for p in PROFILES:
            _, res = sh.inputs_from_hydrus(c["n_C"], c["group"], root_profile=p)
            means[p] = float(np.mean(res.Cw / res.Cw[0]))
        dl = {p: float(np.log10(means[p] / means["linear30"])) for p in PROFILES}
        mor = [dl[p] for p in PROFILES if p.startswith("morita_")]
        out[cong] = dict(mean_c_over_c0=means, log10_vs_linear30=dl,
                         morita_envelope_vs_C=[min(mor) - dl["morita_C"], max(mor) - dl["morita_C"]])
        print(f"{cong:8}" + "".join(f"{means[p]:11.4f}" for p in PROFILES))
        print(f"{'':8}" + "".join(f"{dl[p]:+11.4f}" for p in PROFILES), flush=True)
    with open(OUT, "w") as f:
        json.dump(dict(source="validation/root_profile_check.py (PREREG 10.2; nothing fitted)",
                       profiles=PROFILES, congeners=out), f, indent=1)
    print(f"[written] {os.path.relpath(OUT, ROOT)}")


if __name__ == "__main__":
    main()
