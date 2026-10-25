"""Kim 2019 grain BAF with its non-detects (validation/kim2019_censored.py, pre-registered
in bbb3a8a).

Pinned: (1) the bounds the published averages imply when every non-detect lies in
[0, MDL] -- the ND=0 values are the lower corner and PFDoDA spans 114x; (2) the C11
rise survives the non-detects, the C12 level does not; (3) the scores: both lipid
variants beat the free-anion variants by > 1 log10 on either scoring, but which lipid
variant is better depends on the non-detect treatment, so the pre-registered rule
finds no single winner.
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


def test_scores_and_verdict(result):
    sc = result["scores"]
    assert sc["lipid"]["point"] == pytest.approx(0.50, abs=0.02)
    assert sc["nstem_lipid"]["point"] == pytest.approx(0.43, abs=0.02)
    assert sc["mono"]["point"] == pytest.approx(2.06, abs=0.03)
    assert sc["fit"]["point"] == pytest.approx(0.68, abs=0.03)
    assert sc["nstem_lipid"]["interval"] == pytest.approx(0.18, abs=0.02)
    assert sc["lipid"]["interval"] == pytest.approx(0.28, abs=0.02)
    for v in ("lipid", "nstem_lipid"):
        for k in ("point", "interval"):
            assert sc["mono"][k] - sc[v][k] > 1.0                      # the free anion misses
            assert sc["fit"][k] - sc[v][k] > 0.15                      # so does the Yamazaki fit
    assert all(sc[v]["interval"] <= sc[v]["point"] for v in sc)         # intervals only help
    assert result["winner"] is None


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
