"""Independent check of the ORYZA2000 port (and the shipped oryza_growth) on the Shizukuishi
FACE experiment -- validation/oryza2000_kim_face.py (pre-registered 9f90d55; station addendum
4fca549; AMENDMENT 2 a8d4e0a: corrected observation operator + bootstrap decision).

Pinned: (1) the observations are read and placed by STAGE; (2) the station weather is the JMA
record, ~1.1-1.45 C warmer than the POWER cell; (3) the corrected operator -- the senescent share
is a floor (only a prediction below it counts) and Q3-Q5 are on living biomass; (4) the decisions
under it: on the station weather (IR72 phenology) the port and the shipped driver are
indistinguishable, under S1 the port is better; (5) the shipped driver never matures by Nov 30
at 39.6 N, the port does once the weather is not the cool grid cell.
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


def test_the_corrected_observation_operator():
    v = {s: {"GL": 2.0, "SL": 1.0, "TOT": 10.0} for s in V.SCORED}
    q = V.quantities(v, {"panicle_initiation": (1.0, 10.0)})
    assert q["Q3"]["anthesis"] == pytest.approx(2.0 / 9.0)              # leaf share of the LIVING plant
    assert q["Q4"]["anthesis"] == pytest.approx(1.0)
    assert q["Q5"]["panicle_initiation"] == pytest.approx(1.0 / 9.0)
    assert V.point_error("Q2", 0.8, 0.4) == 0.0                         # above the floor: no error
    assert V.point_error("Q2", 0.2, 0.4) == pytest.approx(-np.log10(2.0))
    assert V.point_error("Q3", 0.8, 0.4) == pytest.approx(np.log10(2.0))  # other groups stay two-sided


def test_the_decisions_under_the_corrected_scoring(out):
    S4 = out["post"]["S4_station"]
    assert S4["rule"]["verdict"] == "indistinguishable" and S4["rule"]["reached"]
    assert S4["A"]["rmse"]["Q2x"] == 0.0                                 # the port is above the floor
    assert S4["A"]["rmse"]["POOLEDx"] == pytest.approx(0.195, abs=0.005)
    assert S4["B"]["rmse"]["POOLEDx"] == pytest.approx(0.183, abs=0.005)
    assert S4["A"]["rmse"]["Q4"] < S4["B"]["rmse"]["Q4"]                 # port times dry matter better
    assert out["res"]["S1_lapse+1.7C"]["rule"]["verdict"] == "A better"
    assert not out["res"]["primary"]["rule"]["reached"]                  # cool cell: the port never matures
    assert out["res"]["primary"]["B"]["rmse"]["Q2x"] == pytest.approx(0.079, abs=0.005)


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
