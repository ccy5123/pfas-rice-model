"""Option C redone with the recorded nursery and the 1999 event dates
(validation/oryza2000_kim_face_nursery_events.py, pre-registered in b1774b4).

Pinned: (1) the file's own code reproduces option C's primary exactly, and the recorded nursery
leaves the shipped arm untouched; (2) under the 1999 events the calibration hits PI 54 / heading 78
and the sample-DVS placement equals option C's on target days; (3) the pre-registered verdicts --
D1 passes on N and does not flip under (s1)-(s3), NE FAILS on (b) while the placement control NG
passes (so the flip is the event dates, not the placement), D4 says the dead-leaf excess is the
leaf-death table, C3 is NOT adequate under either design.
"""
import os
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "validation"))

import oryza2000_kim_face as K                       # noqa: E402
import oryza2000_kim_face_nursery_events as E        # noqa: E402


@pytest.fixture(scope="module")
def out():
    return E.main(verbose=False)


def test_the_code_reproduces_option_c_and_leaves_the_shipped_arm_alone(out):
    assert out["checks"] == {"repro": True, "shipped_untouched": True}
    # the recorded nursery moves only the port's juvenile factor (1998 shortest nursery -> fastest)
    f = {y: out["res"]["N"]["A"]["own"][y][1] for y in E.YEARS}
    g = {y: out["res"]["C_repro"]["A"]["own"][y][1] for y in E.YEARS}
    assert f[1998] > g[1998] and f[1999] < g[1999] and f[2000] < g[2000]


def test_events_calibration_and_the_sample_dvs_placement(out):
    obs, _ = K.load_obs()
    tr = E.DESIGNS["NE"]
    tgt = E.targets(tr, obs)
    assert tgt[1999] == E.EVENTS_1999 and tgt[1998]["anthesis"] == obs[1998]["anthesis"]["dat"]
    for arm in E.ARMS:
        r = out["res"]["NE"][arm]["own_runs"][1999]
        for stage, lim in E.C.STAGE_DVS:
            assert E.C.crossing(r, lim) == pytest.approx(E.EVENTS_1999[stage], abs=0.01)
        sd = out["res"]["NE"][arm]["sample_dvs"]
        for y in (1998, 2000):                       # target days -> option C's stage DVS (to the calibration tolerance)
            assert sd[y]["panicle_initiation"] == pytest.approx(0.65, abs=1e-4)
            assert sd[y]["anthesis"] == pytest.approx(1.0, abs=1e-4)
            assert sd[y]["maturity"] == pytest.approx(2.0, abs=1e-4)
        assert 0.74 < sd[1999]["panicle_initiation"] < 0.78 and 1.24 < sd[1999]["anthesis"] < 1.30


def test_the_pre_registered_verdicts(out):
    R = out["res"]
    assert out["d1"] == {"a": True, "b": True, "c": True, "passed": True}
    assert out["flips"] == []
    assert R["N"]["A"]["rmse"]["POOLEDx"] == pytest.approx(0.220, abs=0.005)
    assert R["N"]["B"]["rmse"]["POOLEDx"] == pytest.approx(0.244, abs=0.005)
    assert out["d2_NG"]["passed"] and not out["d2_NE"]["passed"]
    assert out["d2_NE"]["a"] and not out["d2_NE"]["b"]                      # NE fails on (b) only
    assert R["NE"]["A"]["rmse"]["Q2x"] - R["NE"]["B"]["rmse"]["Q2x"] == pytest.approx(0.100, abs=0.01)
    assert abs(R["NG"]["A"]["rmse"]["POOLEDx"] - R["N"]["A"]["rmse"]["POOLEDx"]) < 0.01
    assert not R["NE_s1_tmpsb0"]["c1"]["passed"]
    assert out["d4_drop"] == pytest.approx(0.061, abs=0.005) and out["d4_drop"] < 0.10
    assert out["c3_N"] is False and out["c3_NE"] is False
    assert abs(R["NE"]["A"]["phen"][1999]["anthesis"]) < 3                  # vs the true heading


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
