"""Tests for the Bayesian inverse demo (validation/bayesian_inverse_demo.py).

Locks in the identifiability result discussed in the modelling notes: from tissue
C(t) the exposure (Q_TP scale, Cwo level) IS jointly identifiable when transport
is fixed, but Q_TP and f_xy collapse onto a product ridge (only their product is
constrained). Run noise-free for determinism.
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "validation"))
import bayesian_inverse_demo as bid  # noqa: E402


def test_exposure_identifiable_when_transport_fixed():
    """(A) Q_TP-scale & Cwo-level recover from tissue C(t) with transport fixed, and
    are far better conditioned than the transport ridges (B/C)."""
    a = bid.run_scenario("A", ["qtp_scale", "cwo_level"], noise=False)
    b = bid.run_scenario("B", ["qtp_scale", "f_xy_mult"], noise=False)
    assert abs(a["corr"]) < 0.95                         # not a hard ridge
    assert a["cond"] < 0.3 * b["cond"]                   # much better conditioned than the ridge
    assert a["recovered"]["qtp_scale"] == pytest.approx(1.0, abs=0.12)
    assert a["recovered"]["cwo_level"] == pytest.approx(1.0, abs=0.12)


def test_qtp_fxy_is_a_product_ridge():
    """(B) Q_TP and f_xy collapse onto a product ridge: strong correlation, ill-
    conditioned, and only the PRODUCT Q_TP·f_xy is recovered."""
    r = bid.run_scenario("B", ["qtp_scale", "f_xy_mult"], noise=False)
    assert abs(r["corr"]) > 0.9                          # near-±1 ridge
    assert r["cond"] > 100                               # ill-conditioned
    assert r["product"] == pytest.approx(1.0, abs=0.15)  # only the PRODUCT is recovered
