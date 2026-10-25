"""Option C -- both biomass drivers' development rates calibrated to Akitakomachi, then the
Shizukuishi FACE comparison re-run (validation/oryza2000_kim_face_phenology.py, pre-registered
in 7743bfc; scoring per AMENDMENT 2 of oryza2000_kim_face.py).

Pinned: (1) the DRATES-style calibration hits the observed PI / anthesis / maturity days; (2) the
leave-one-year-out rates are the other two years' means; (3) the decisions under the corrected
scoring -- the port is better (LOYO, station weather) under every sensitivity and with own-year
rates; C3 says LOYO phenology is NOT adequate.
"""
import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "validation"))

import oryza2000_kim_face as K               # noqa: E402
import oryza2000_kim_face_phenology as C     # noqa: E402


@pytest.fixture(scope="module")
def out():
    return C.main(verbose=False)


def test_calibration_hits_the_observed_days():
    obs, _ = K.load_obs()
    tr = C.TREATMENTS["primary"]
    for arm in C.ARMS:
        rates, f = C.calibrate(arm, 1999, tr, obs)
        r = C.run(arm, 1999, rates, tr)
        for stage, lim in C.STAGE_DVS:
            assert C.crossing(r, lim) == pytest.approx(obs[1999][stage]["dat"], abs=0.01)
        base = C.base_rates(arm)
        assert rates["DVRJ"] / base["DVRJ"] == pytest.approx(f)
        assert rates["DVRI"] / base["DVRI"] == pytest.approx(f)


def test_loyo_rates_are_the_other_years_means(out):
    for arm in C.ARMS:
        own, lr = out["res"]["primary"][arm]["own"], out["res"]["primary"][arm]["loyo_rates"]
        assert lr[1999]["DVRR"] == pytest.approx(np.mean([own[1998][0]["DVRR"], own[2000][0]["DVRR"]]))
        assert lr[1998]["DVRP"] == pytest.approx(np.mean([own[1999][0]["DVRP"], own[2000][0]["DVRP"]]))


def test_the_verdicts_under_the_corrected_scoring(out):
    P = out["res"]["primary"]
    assert out["c1"]["verdict"] == "A better" and out["c1"]["reached"]
    assert out["flips"] == []
    assert P["c1_own"]["verdict"] == "A better"
    assert P["A"]["rmse"]["POOLEDx"] == pytest.approx(0.168, abs=0.005)
    assert P["B"]["rmse"]["POOLEDx"] == pytest.approx(0.240, abs=0.005)
    assert P["A"]["rmse"]["Q1"] == pytest.approx(0.174, abs=0.005)
    assert P["A"]["rmse"]["Q2x"] == 0.0 and P["B"]["rmse"]["Q2x"] == pytest.approx(0.207, abs=0.005)
    assert out["c3"] is False                                         # LOYO phenology misses ~10 d
    assert abs(P["A"]["phen"][1998]["anthesis"]) > 5
    for name, R in out["res"].items():
        assert all(R["A"]["detail"][y]["reached"]["maturity"] for y in C.YEARS), name
        assert all(R["B"]["detail"][y]["reached"]["maturity"] for y in C.YEARS), name


def test_under_the_attached_reading_the_port_is_not_better(out):
    """POST-HOC: with Q2 scored symmetrically (no blade had detached) the port's advantage on
    the FACE data is not resolved under any treatment -- it is indistinguishable, never worse."""
    assert {c["verdict"] for c in out["c1_attached"].values()} == {"indistinguishable"}
    assert out["c1_attached"]["primary"]["delta"] == pytest.approx(0.023, abs=0.005)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
