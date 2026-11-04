"""PREREG §11 (amendment 4): the sourced whole-cell-wall K_cw."""
import math
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import kcw_sourced as K  # noqa: E402


def test_upper_ends_match_the_registered_values():
    u = K.makselon_upper()
    reg = {"PFOA": 9.3, "PFDoDA": 28.8, "PFTeDA": 95.1, "PFBS": 10.2, "PFHxS": 4.4, "PFOS": 22.1}
    for c, v in reg.items():
        assert round(u[c], 1) == v
    assert "PFUnDA" not in u                      # 3 runs at pH >= 5.8: excluded


def test_pfoa_pfos_from_the_11_7_reading():
    pfoa = K.kcw_sourced(7, "carboxylate")
    assert pfoa.case == "domain"
    assert pfoa.low == pytest.approx(10 ** -0.29) and pfoa.high == pytest.approx(9.329, rel=1e-3)
    assert pfoa.rec == pytest.approx(math.sqrt(pfoa.low * pfoa.high))
    pfos = K.kcw_sourced(8, "sulfonate")
    assert pfos.low == pytest.approx(10 ** -0.11)          # composted bark (PREREG 11.7 reading)


def test_unread_registered_compound_is_an_error_not_out_of_domain():
    lower = {c: v for c, v in K.sorengard_lower().items() if c != "PFBS"}
    with pytest.raises(K.PendingReading):
        K.kcw_sourced(4, "sulfonate", lower=lower)


def test_prereg_11_7_cases():
    assert K.kcw_sourced(6, "carboxylate").case == "removed"          # PFHpA not plotted
    assert K.kcw_sourced(3, "carboxylate").case == "lower_only"       # PFBA
    assert K.kcw_sourced(13, "carboxylate").case == "domain"          # PFTeDA plotted
    assert K.kcw_sourced(5, "carboxylate", n_ether_O=1).case == "removed"   # GenX
    assert any("min_may_be_lower" in f for f in K.kcw_sourced(4, "sulfonate").flags)


def _toy():
    upper = {"PFOA": 10.0, "PFDoDA": 1000.0, "PFBS": 10.0, "PFOS": 100.0}
    lower = {c: {"Kd": 1.0, "sorbent": "Peat", "low_conf": False} for c in K.SORENGARD_REGISTERED}
    lower["PFPeA"] = None                                   # read, no quantifiable sorption
    return upper, lower


def test_decision_tree():
    up, lo = _toy()
    r = K.kcw_sourced(9, "carboxylate", upper=up, lower=lo)  # PFDA: log-linear upper, 100
    assert r.case == "domain" and r.high == pytest.approx(100.0) and r.rec == pytest.approx(10.0)
    assert any("interpolated" in f for f in r.flags)
    r = K.kcw_sourced(3, "carboxylate", upper=up, lower=lo)  # PFBA: lower only
    assert r.case == "lower_only" and r.rec == r.low == r.high == 1.0
    r = K.kcw_sourced(4, "carboxylate", upper=up, lower=lo)  # PFPeA: nothing quantifiable
    assert r.case == "removed" and r.rec == 0.0 and r.nearest_rec == pytest.approx(math.sqrt(10.0))
    r = K.kcw_sourced(5, "carboxylate", n_ether_O=1, upper=up, lower=lo)   # GenX-like ether
    assert r.case == "removed" and r.rec == 0.0
    r = K.kcw_sourced(7, "carboxylate", n_ether_O=2, upper=up, lower=lo)   # ether in the domain
    assert r.case == "domain" and r.rec == pytest.approx(math.sqrt(10.0))
    assert r.low == pytest.approx(1.0 * 10 ** (-0.98)) and r.high == pytest.approx(10.0 * 10 ** 0.98)
    r = K.kcw_sourced(8, "sulfonamide", upper=up, lower=lo)
    assert r.case == "removed"
    r = K.kcw_sourced(12, "carboxylate", upper=up, lower=lo)  # PFTrDA: no reading, beyond the upper ends
    assert r.case == "removed"
