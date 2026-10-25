"""Kim 2019 PFOA grain BAF per site (validation/kim2019_pfoa_anchor.py, pre-registered in bedadbf).

Pinned: (1) the main-text Table 1 transcription reproduces SI Table S6; (2) the candidate
estimators and A-2's weight; (3) the admissible censored per-site estimate, its interval and
the constant literature_params carries; (4) the per-site ratio is not a plant constant
(slope < 1) while the model is linear; (5) the decisions, including that D2 rests on the
interval METHOD (the profile likelihood contains the published value).
"""
import contextlib
import io
import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "validation"))

import literature_params as lp  # noqa: E402
import kim2019_pfoa_anchor as A  # noqa: E402


def test_table1_reproduces_si_table_s6():
    rows = A.load_sites()
    assert len(rows) == 30 and "D-5" not in {r["site"] for r in rows}
    for key, n, mean, det in (("rice_PFOA_ngg", 30, 0.349, 17),
                              ("water_PFOA_ngL", 27, 78.7, 26),
                              ("soil_PFOA_ngg_dw", 30, 0.160, 13)):
        vals = [r[key] for r in rows if r[key] != "NS"]
        assert len(vals) == n
        assert sum(v != "ND" for v in vals) == det
        assert np.mean([0.0 if v == "ND" else float(v) for v in vals]) == pytest.approx(mean, rel=2e-3)   # SI: 3 sig. figures
    pw, soil, rice, df = lp.KIM2019_FIELD["PFOA"]
    assert (pw, soil, rice) == (78.7, 0.160, 0.349) and df == pytest.approx(100 * 17 / 30, abs=0.5)


def test_candidate_estimators_and_a2_weight():
    est = A.estimators(A.load_sites())
    assert est["E1"] == pytest.approx(lp.kim2019_grain_baf()["PFOA"], rel=0.002)
    assert est["E2"] == pytest.approx(4.18, abs=0.01)
    assert est["E3"] == pytest.approx(21.45, abs=0.01)
    assert est["E4"] == pytest.approx(26.96, abs=0.01) and est["n_both"] == 14
    assert est["E4_median"] == pytest.approx(28.3, abs=0.05)
    assert est["a2_water_share"] == pytest.approx(0.828, abs=0.001)


def test_censored_per_site_estimate_and_slope():
    d = A.ratio_data(A.load_sites())
    assert len(d) == 26 and sum(r["cens"] for r in d) == 12
    an = A.anchor(d)
    e5, lo, hi = lp.KIM2019_PFOA_SITE_BAF
    assert an["E5"] == pytest.approx(e5, abs=0.01)
    assert an["E5_ci"] == pytest.approx((lo, hi), abs=0.01)
    assert an["km_median"] == pytest.approx(7.77, abs=0.02)
    assert an["slope"] == pytest.approx(0.52, abs=0.01)
    assert an["slope_ci"][1] < 1.0                                   # not a plant constant
    est = A.estimators(A.load_sites())
    assert est["E1"] < lo and lo - est["E1"] < 0.2                   # D2: outside, by a hair
    plo, phi = A.profile_interval(d)
    assert plo < est["E1"] < phi                                     # ... and method-dependent
    assert (plo, phi) == pytest.approx((4.18, 19.78), abs=0.01)


@pytest.fixture(scope="module")
def result():
    with contextlib.redirect_stdout(io.StringIO()):
        return A.main()


def test_model_is_linear_and_below_the_site_interval(result):
    lo, hi = result["anchor"]["E5_ci"]
    for v, p in result["pred"].items():
        assert p < lo, v
        assert result["slope_model"][v] == pytest.approx(0.964, abs=0.005), v
    assert result["pred"]["nstem_lipid"] == pytest.approx(1.237, abs=0.03)
    assert result["a2_sum17"] == pytest.approx(0.673, abs=0.001)
    assert result["d2"] and result["d4"]


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
