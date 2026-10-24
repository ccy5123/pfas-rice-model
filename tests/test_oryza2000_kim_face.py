"""Independent check of the ORYZA2000 port (and the shipped oryza_growth) on the Shizukuishi
FACE experiment -- validation/oryza2000_kim_face.py, pre-registered in 9f90d55, post-hoc
addendum (station temperatures, finite Q2) in 4fca549.

Pinned: (1) the observations are read and placed by STAGE (five samplings = tillering, PI,
anthesis, mid-ripening, maturity; mid-ripening fractions from the sampling DATs); (2) the
station weather is the JMA record, ~1.1-1.45 C warmer than the POWER cell; (3) the
pre-registered verdict -- D1 FAILS on the primary weather and flips only under S1; (4) the
post-hoc reading -- on station temperatures the rule as written passes only because the
shipped driver's senescent share at anthesis is exactly 0 (infinite Q2), and the finite
variant fails on (b) in every treatment; (5) the shipped driver never matures by Nov 30 at
39.6 N, the port does once the weather is not the cool grid cell.
"""
import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "validation"))

import oryza2000_kim_face as V   # noqa: E402


@pytest.fixture(scope="module")
def out():
    return V.main(verbose=False)


def test_observations_are_placed_by_stage():
    obs, root = V.load_obs()
    assert [obs[y]["panicle_initiation"]["dat"] for y in V.YEARS] == [53, 62, 56]
    assert [obs[y]["anthesis"]["dat"] for y in V.YEARS] == [81, 88, 78]
    assert [obs[y]["maturity"]["dat"] for y in V.YEARS] == [131, 123, 120]
    assert obs[1999]["anthesis"]["GL"] == 187.3 and obs[2000]["maturity"]["SL"] == 85.58
    assert obs[1999]["panicle_initiation"]["TOT"] == 439.0                 # Fig. 2 (overlay-verified; first pass 435)
    assert root == {"panicle_initiation": (45.1, 425.94), "maturity": (37.94, 1495.64)}
    fr = [V.ripening_fraction(obs[y]) for y in V.YEARS]
    assert fr == pytest.approx([0.56, 0.486, 0.524], abs=1e-3)


def test_station_weather_is_warmer_than_the_grid_cell():
    for y in V.YEARS:
        on, ox = V.station_offset(y)
        assert 0.8 < on < 1.5 and 0.8 < ox < 1.7
    w_st, w_pw = V.load_weather(1999, source="station"), V.load_weather(1999)
    assert np.array_equal(w_st["rad_mj"], w_pw["rad_mj"])                  # AMeDAS has no radiation
    assert np.all(np.isfinite(w_st["tmin"])) and np.all(np.isfinite(w_st["tmax"]))
    with pytest.raises(ValueError):
        V.load_weather(1999, source="jma")


def test_the_pre_registered_verdict(out):
    P = out["res"]["primary"]
    assert out["verdict"] is False and out["flips"] == ["S1_lapse+1.7C"]
    assert P["d1"] == {"a": True, "b": False, "c": False, "passed": False}
    assert P["A"]["rmse"]["Q1"] == pytest.approx(0.457, abs=0.01)
    assert P["B"]["rmse"]["Q1"] == pytest.approx(0.205, abs=0.01)
    assert P["A"]["rmse"]["Q4"] == pytest.approx(0.102, abs=0.01)
    assert np.isinf(P["B"]["rmse"]["Q2"])                                  # zero senescence at anthesis
    assert out["res"]["S1_lapse+1.7C"]["d1"]["passed"]


def test_the_post_hoc_station_reading(out):
    S4 = out["post"]["S4_station"]
    assert S4["d1"]["passed"]                                              # only because B's Q2 is inf
    assert S4["d1x"] == {"a": True, "b": False, "c": True, "passed": False}
    assert S4["A"]["rmse"]["Q2x"] == pytest.approx(0.328, abs=0.01)
    assert S4["B"]["rmse"]["Q2x"] == pytest.approx(0.172, abs=0.01)
    assert S4["A"]["rmse"]["Q4"] < S4["B"]["rmse"]["Q4"]                   # port times dry matter better
    assert S4["A"]["rmse"]["Q3"] > S4["B"]["rmse"]["Q3"]                   # shipped keeps leaf share better
    assert not any(R["d1x"]["passed"] for R in list(out["res"].values()) + [S4])


def test_who_reaches_maturity(out):
    every = list(out["res"].items()) + list(out["post"].items())
    for name, R in every:
        assert not any(R["B"]["detail"][y]["reached"]["maturity"] for y in V.YEARS), name
    for name in ("S1_lapse+1.7C",):
        assert all(out["res"][name]["A"]["detail"][y]["reached"]["maturity"] for y in V.YEARS)
    assert all(out["post"]["S4_station"]["A"]["detail"][y]["reached"]["maturity"] for y in V.YEARS)
    assert not any(out["res"]["primary"]["A"]["detail"][y]["reached"]["maturity"] for y in V.YEARS)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
