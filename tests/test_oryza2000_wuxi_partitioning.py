"""Organ partitioning at the Wuxi FACE (Yang 2006), a second independent check of the biomass
drivers (validation/oryza2000_wuxi_partitioning.py, pre-registered in 54bf8c3).

Pinned: (1) the data and weather load as transcribed/downloaded (fractions sum to 1; the weather
check against Hasegawa's season means triggers the offset run for 2003 only); (2) the Wp design hits
its phenology targets; (3) the decisions (paired bootstrap, AMENDMENT 2) -- the port is better, under Wp and every
sensitivity; Gao's partitioning alone is marginally better than IR72's, the full Gao set
indistinguishable.
"""
import os
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "validation"))

import oryza2000_kim_face_phenology as C      # noqa: E402
import oryza2000_wuxi_partitioning as W        # noqa: E402


@pytest.fixture(scope="module")
def out():
    return W.main(verbose=False)


def test_data_and_weather():
    obs = W.load_obs()
    for y in W.YEARS:
        assert set(obs[y]) == {"tillering", "panicle_initiation", "heading", "mid_ripening", "maturity"}
        for v in obs[y].values():
            assert v["leaf"] + v["stem"] + v["spike"] == pytest.approx(1.0, abs=0.0015)
    assert obs[2001]["panicle_initiation"]["dat"] == 51 and obs[2003]["maturity"]["AGB"] == 1674.0
    wc = W.weather_check()
    assert [abs(wc[y]["diff"]) > 0.5 for y in W.YEARS] == [False, False, True]
    w = W.load_weather(2002)
    assert len(w["rad_mj"]) == 365 and w["doy0"] == 1


def test_wp_hits_its_targets(out):
    obs = out["obs"]
    for arm in ("A", "B"):
        r = out["res"]["Wp"][arm]["runs"][2002]
        assert C.crossing(r, 0.65) == pytest.approx(obs[2002]["panicle_initiation"]["dat"], abs=0.01)
        assert C.crossing(r, 2.0) == pytest.approx(W.MATURITY_DAT, abs=0.01)


def test_the_verdicts(out):
    R = out["res"]
    assert out["r1"]["verdict"] == "A better" and out["r1"]["reached"]
    assert out["flips"] == [] and "W_Toffset" in R
    assert R["W"]["A"]["rmse"]["POOLED"] == pytest.approx(0.076, abs=0.005)
    assert R["W"]["B"]["rmse"]["POOLED"] == pytest.approx(0.129, abs=0.005)
    assert R["W"]["A"]["rmse"]["W2"] < R["W"]["B"]["rmse"]["W2"] - 0.04          # spike share
    assert R["W"]["B"]["bias"]["W2"] < -0.08                                       # under-filled panicles
    assert out["r3"]["A+PART"]["verdict"] == "A better"                            # Gao partitioning, marginally
    assert out["r3"]["G"]["verdict"] == "indistinguishable"
    d = R["W"]["A"]["detail"]                                                      # IR72 matures early at Wuxi
    assert all(d[y]["reached"]["maturity"] and d[y]["dat"]["maturity"] < 115 for y in W.YEARS)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
