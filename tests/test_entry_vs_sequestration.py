"""Guards for B2 (validation/entry_vs_sequestration.py).

What is pinned is the thing a later session is most likely to get wrong, because
this run's own first reading got it wrong until a control said otherwise:

  1. the CONTROL result — removing the sequestration term while keeping lipid
     loading flattens the entry trend just as far. Sequestration is NOT what
     absorbs the chain-length dependence, even though the full two-pool arm
     appears to show that it is;
  2. the baseline the whole test speaks from — the single pool really does need
     a strongly chain-length-dependent entry conductance;
  3. that the dependence is REDUCED, never eliminated. A later session reading
     "lipid loading explains the chain-length trend" has over-read it.

Plus the API contract: the new uptake knobs on the two-pool must default to the
cached fit exactly.

Deliberately coarse (few congeners, short grid) — these pin directions and
orderings, not the RMSEs, which are in-sample and driver-dependent.
"""
import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "validation"))

import model_api as api                                          # noqa: E402
import entry_vs_sequestration as EV                              # noqa: E402

GRID = [0.5, 2.0, 5.0, 20.0, 50.0, 200.0]
SHORT, LONG = "PFHxA", "PFDoDA"


# ---------------------------------------------------------------------------
# API contract
# ---------------------------------------------------------------------------
def test_twopool_uptake_knobs_default_to_the_cached_fit():
    a = api.simulate_twopool_seq("PFOA")
    b = api.simulate_twopool_seq("PFOA", uptake="carrier")
    for k in api.TISSUES:
        assert a["baf_final"][k] == b["baf_final"][k], k


def test_twopool_rejects_an_unknown_uptake_mode():
    with pytest.raises(ValueError):
        api.simulate_twopool_seq("PFOA", uptake="apoplast")


def test_twopool_bypass_mode_turns_the_carrier_off():
    r = api.simulate_twopool_seq("PFOA", uptake="bypass")
    explicit = api.simulate_twopool_seq(
        "PFOA", vmax_scale=0.0, g_apo=api.UPTAKE_MODES["bypass"]["g_apo"])
    for k in api.TISSUES:
        assert r["baf_final"][k] == explicit["baf_final"][k], k


# ---------------------------------------------------------------------------
# the findings
# ---------------------------------------------------------------------------
def _best(arm, nm, obs):
    return min(((v, EV.rmse(arm, [nm], obs, v)) for v in GRID),
               key=lambda vr: vr[1])[0]


def test_the_single_pool_needs_a_chain_length_dependent_entry_term():
    """Finding 2 — the baseline. Without it there is nothing for anything else to
    remove, and the run declares itself inconclusive (gate G0)."""
    obs = EV.load_obs()
    short, long = _best("S", SHORT, obs), _best("S", LONG, obs)
    assert long > short, (short, long)
    assert long / short >= 10.0, (short, long)


def test_removing_sequestration_does_not_restore_the_trend():
    """Finding 1, the CONTROL, and the reason B2's answer is not the one it was
    expected to be. The full two-pool flattens the entry trend — but so does the
    SAME model with the sequestration term switched off and only lipid loading
    left, which means sequestration is not what did the flattening. If a later
    session finds the control restoring the single pool's large spread, the
    attribution flips back and that is a real finding, not a test to delete."""
    obs = EV.load_obs()
    s_ratio = _best("S", LONG, obs) / _best("S", SHORT, obs)
    c_ratio = _best("C", LONG, obs) / _best("C", SHORT, obs)
    assert c_ratio < s_ratio / 5.0, (s_ratio, c_ratio)


def test_lipid_loading_alone_flattens_it_in_the_plain_single_pool():
    """The same conclusion reached without any second root pool at all: adding
    lipid loading to the ordinary single-pool model already removes most of the
    chain-length requirement."""
    obs = EV.load_obs()
    s_ratio = _best("S", LONG, obs) / _best("S", SHORT, obs)
    sl_ratio = _best("S+L", LONG, obs) / _best("S+L", SHORT, obs)
    assert sl_ratio < s_ratio / 5.0, (s_ratio, sl_ratio)


def test_the_chain_length_freedom_is_reduced_not_eliminated():
    """Finding 3, the anti-over-claim guard. Even in the arms that flatten the
    trend, letting the conductance vary per congener still fits better than one
    global value — so 'lipid loading explains the chain-length dependence' is an
    over-reading of what this measured."""
    obs = EV.load_obs()
    names = [n for n in EV.load_obs() if n in api._CONG]
    for arm in ("S+L", "C"):
        per = EV.fit_per_congener(arm, names, obs, GRID)
        gv, r_global = EV.fit_global(arm, names, obs, GRID)
        r_per = EV.rmse(arm, names, obs, per)
        assert r_per < r_global, (arm, r_per, r_global)


def test_nothing_was_adopted_by_b2():
    """Lipid loading stays opt-in and k_seq stays unpromoted, whatever B2 found."""
    import json
    with open(os.path.join(ROOT, "params", "parameters.json")) as f:
        P = json.load(f)
    assert not any("k_seq" in k for k in P), sorted(P)
    r = api.simulate("PFOA", season=120.0, n_t=121)
    assert r["params"]["g_xy"] == 0.0 and r["params"]["g_ph"] == 0.0
