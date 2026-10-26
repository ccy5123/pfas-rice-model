"""Kim 2019 grain BAF with its non-detects (validation/kim2019_censored.py, pre-registered
in bbb3a8a).

Pinned: (1) the bounds the published averages imply when every non-detect lies in
[0, MDL] -- the ND=0 values are the lower corner and PFDoDA spans 114x; (2) the C11
rise survives the non-detects, the C12 level does not; (3) the scores equal the transfer
record (validation/transfer_sourced.json) of the same fitted variants, and the interval
scoring never exceeds the point scoring.
"""
import contextlib
import io
import os
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "validation"))

import literature_params as lp  # noqa: E402
import kim2019_censored as K    # noqa: E402


@pytest.fixture(scope="module")
def result():
    with contextlib.redirect_stdout(io.StringIO()):
        return K.main()


def test_bounds_from_detection_frequency_and_mdl():
    b = lp.kim2019_grain_baf_bounds()
    pub = lp.kim2019_grain_baf("porewater")
    for nm, (lo, pt, hi) in b.items():
        assert pt == pytest.approx(pub[nm])
        assert lo <= pt <= hi
    assert b["PFDoDA"][0] == pytest.approx(5.5, abs=0.05)
    assert b["PFDoDA"][2] == pytest.approx(629.5, abs=1.0)
    assert b["PFUnDA"][0] == pytest.approx(19.83, abs=0.05)
    assert b["PFOA"][2] / b["PFOA"][0] < 1.15                     # well detected: narrow
    for nm in ("PFHpA", "PFNA", "PFOA"):                         # water well detected:
        assert b[nm][0] == pytest.approx(b[nm][1], rel=0.01)      # ND=0 is ~the lower corner


def test_the_c11_rise_survives_the_non_detects_and_c12_does_not(result):
    s = result["steps"]
    assert s[("PFDA", "PFUnDA")] and s[("PFNA", "PFUnDA")]
    assert not s[("PFDA", "PFDoDA")]


def test_scores_match_the_transfer_record(result):
    import json
    import numpy as np
    with open(os.path.join(ROOT, "validation", "transfer_sourced.json")) as f:
        rec = json.load(f)["errors"]
    same = {"nstem_lipid": "redistributed_lipid", "lipid": "4pool_lipid", "mono": "4pool_free",
            "fit": "redistributed_fit"}
    sc = result["scores"]
    for v, key in same.items():
        for k in ("point", "interval"):
            r = float(np.sqrt(np.mean(np.square(rec[key][f"kim_{k}"]))))
            assert sc[v][k] == pytest.approx(r, abs=1e-6), (v, k)
    assert all(sc[v]["interval"] <= sc[v]["point"] + 1e-12 for v in sc)    # intervals only help


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
