"""The Gao 2021 cold-region japonica crop set in the ORYZA2000 port, on the Shizukuishi FACE data
(validation/oryza2000_kim_face_gao_japonica.py, pre-registered in 5aada97).

Pinned: (1) the Gao values reach the port as transcribed (tables, SLA switch, single KDF) and the
IR72.D92 / shipped arms reproduce the redo exactly; (2) the decisions under the corrected scoring
(AMENDMENT 2 of oryza2000_kim_face.py) -- the Gao set is better than the shipped driver on N,
indistinguishable on NE, and indistinguishable from the IR72 port; (3) the decomposition and Gao's
own phenology.
"""
import os
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "validation"))

import oryza2000_kim_face_gao_japonica as G     # noqa: E402


@pytest.fixture(scope="module")
def out():
    return G.main(verbose=False)


def test_the_gao_crop_set_reaches_the_port():
    c = G.crop("G")
    assert c["SWISLA"] == "TABLE" and list(c["SLATB"][1][:3]) == [0.0010, 0.0015, 0.0020]
    assert list(c["DRLVT"][0]) == [0.0, 0.4, 1.0, 1.6, 2.1] and c["DRLVT"][1][2] == 0.024
    assert list(c["FSOTB"][1]) == [0.0, 0.0, 0.10, 0.40, 1.0, 1.0]
    assert set(c["KDFTB"][1]) == {0.4} and c["FRPAR"] == 0.45 and c["RGRLMX"] == 0.0088
    assert G.base_rates("G") == {"DVRJ": 0.00089, "DVRI": 0.000758, "DVRP": 0.001055, "DVRR": 0.002}
    a = G.crop("A")
    assert a["SWISLA"] == "FUNCTION" and a["FRPAR"] == 0.5          # IR72.D92 untouched
    assert "COLDMIN" not in c                                       # not applied (no cold-stop rule)


def test_the_verdicts_under_the_corrected_scoring(out):
    R = out["res"]
    assert out["checks"] == {"N": True, "NE": True}
    assert out["g2"]["N"]["verdict"] == "A better"                       # Gao set better than shipped (N)
    assert out["g2"]["NE"]["verdict"] == "indistinguishable"
    assert out["g3"]["N"]["verdict"] == "indistinguishable" == out["g3"]["NE"]["verdict"]   # Gao vs IR72
    assert out["flips"] == []
    assert R["N"]["G"]["rmse"]["Q2x"] == 0.0
    assert R["N"]["G"]["rmse"]["Q1"] < R["N"]["A"]["rmse"]["Q1"] - 0.05
    assert R["N"]["G"]["rmse"]["Q4"] < R["N"]["A"]["rmse"]["Q4"] - 0.05


def test_decomposition_and_gaos_own_phenology(out):
    N = out["res"]["N"]
    assert N["A+DRLVT"]["rmse"]["Q1"] < N["A"]["rmse"]["Q1"]             # the leaf-death table carries the Q1 gain
    assert min(N[a]["rmse"]["POOLEDx"] for a in G.MAIN_ARMS + G.DECOMP_ARMS) == N["A+CANOPY"]["rmse"]["POOLEDx"]
    ph = out["unc"]["phen"]                                         # Gao's own rates, no calibration
    assert ph[1998]["anthesis"] == pytest.approx(9.0, abs=0.3) and ph[1999]["anthesis"] == pytest.approx(-8.7, abs=0.3)
    assert all(ph[y]["maturity"] is not None for y in G.YEARS)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
