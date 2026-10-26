"""Tests for the Bayesian inverse demo (validation/bayesian_inverse_demo.py).

Locks in the identifiability result on the calibration site-year with the harvest
operator: from tissue concentrations the exposure (Q_TP scale, Cwo level) IS jointly
identifiable when transport is fixed; Q_TP and f_xy trade off strongly (they share the
root->shoot flux, so their product is pinned); Cwo and the root-uptake conductance form
a far worse ridge. Run noise-free for determinism; each scenario is solved once.
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "validation"))
import bayesian_inverse_demo as bid  # noqa: E402


@pytest.fixture(scope="module")
def scenarios():
    return dict(
        A=bid.run_scenario("A", ["qtp_scale", "cwo_level"], noise=False),
        B=bid.run_scenario("B", ["qtp_scale", "f_xy_mult"], noise=False),
        C=bid.run_scenario("C", ["cwo_level", "r_eq_mult"], noise=False),
    )


def test_exposure_identifiable_when_transport_fixed(scenarios):
    """(A) Q_TP-scale & Cwo-level recover from the tissue data with transport fixed."""
    a = scenarios["A"]
    assert abs(a["corr"]) < 0.95                         # not a hard ridge
    assert a["recovered"]["qtp_scale"] == pytest.approx(1.0, abs=0.12)
    assert a["recovered"]["cwo_level"] == pytest.approx(1.0, abs=0.12)


def test_qtp_fxy_trade_off(scenarios):
    """(B) Q_TP and f_xy are strongly anti-correlated (they share the root->shoot flux)
    and their product is pinned. Whether the two separate individually sits near the
    0.95 threshold and moves with the phloem allocation, so it is not asserted."""
    b = scenarios["B"]
    assert b["corr"] < -0.9
    assert b["product"] == pytest.approx(1.0, abs=0.15)


def test_cwo_conductance_is_a_ridge(scenarios):
    """(C) Cwo and the root equilibrium ratio enter the linear root entry only as a product:
    far worse conditioned than either other pair, correlation ~ -1, the product pinned."""
    a, b, c = scenarios["A"], scenarios["B"], scenarios["C"]
    assert abs(c["corr"]) > 0.99
    assert c["cond"] > 1e3 * max(a["cond"], b["cond"])
    assert c["product"] == pytest.approx(1.0, abs=0.02)


def test_emcee_posterior_recovers_well_posed():
    """The full-MCMC cross-check (affine-invariant sampler) recovers the well-posed
    exposure, agreeing with the Laplace verdict. Tiny chain for speed; skips when
    emcee is absent (the forward ODE is ~0.7 s/sample, so this is opt-in).

    Seeded (`seed=0`, the default): a 30-step chain is noisy, and while the sampler
    drew from OS entropy this test passed or failed by luck -- on fresh organ mass
    one CI run landed at cwo_level 0.39. Seeds 0-7 all land within the bound."""
    pytest.importorskip("emcee")
    out = bid.emcee_posterior(["qtp_scale", "cwo_level"], nstep=30, nwalk=6, discard=10)
    assert out is not None and out["n"] > 0
    assert set(out["median"]) == {"qtp_scale", "cwo_level"}
    # recovers the well-posed exposure (a short chain → loose bound); the product is
    # also pinned. The exact sample correlation is chain-noisy, so we don't assert it
    # here — the structural identifiability is locked by the Laplace tests above.
    assert out["median"]["qtp_scale"] == pytest.approx(1.0, rel=0.6)   # within ~1.6x
    assert out["median"]["cwo_level"] == pytest.approx(1.0, rel=0.6)
    assert 0.4 < out["product"] < 2.5
