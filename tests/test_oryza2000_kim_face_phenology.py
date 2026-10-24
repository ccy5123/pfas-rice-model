"""Option C -- both biomass drivers' development rates calibrated to Akitakomachi, then the
Shizukuishi FACE comparison re-run (validation/oryza2000_kim_face_phenology.py, pre-registered
in 7743bfc).

Pinned: (1) the DRATES-style calibration hits the observed PI / anthesis / maturity days; (2) the
leave-one-year-out rates are the other two years' means; (3) the pre-registered verdicts -- C1
passes (LOYO, station weather) and does not flip under any sensitivity, C2 says the port's
dead-leaf excess is mostly its leaf-death table, C3 says LOYO phenology is NOT adequate.
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


def test_the_pre_registered_verdicts(out):
    P = out["res"]["primary"]
    assert out["c1"] == {"a": True, "b": True, "c": True, "passed": True}
    assert out["flips"] == []
    assert P["c1_own"]["passed"]
    assert P["A"]["rmse"]["POOLEDx"] == pytest.approx(0.222, abs=0.005)
    assert P["B"]["rmse"]["POOLEDx"] == pytest.approx(0.245, abs=0.005)
    assert P["A"]["rmse"]["Q1"] == pytest.approx(0.174, abs=0.005)
    assert P["A"]["rmse"]["Q2x"] - P["B"]["rmse"]["Q2x"] == pytest.approx(0.035, abs=0.005)  # (b) margin
    assert out["q2x_uncal"] == pytest.approx(C.UNCALIBRATED_PORT_Q2X, abs=0.002)
    assert out["c2_drop"] == pytest.approx(0.072, abs=0.005) and not out["c2_phenological"]
    assert out["c3"] is False                                         # LOYO phenology misses ~10 d
    assert abs(P["A"]["phen"][1998]["anthesis"]) > 5
    for name, R in out["res"].items():
        assert all(R["A"]["detail"][y]["reached"]["maturity"] for y in C.YEARS), name
        assert all(R["B"]["detail"][y]["reached"]["maturity"] for y in C.YEARS), name


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
