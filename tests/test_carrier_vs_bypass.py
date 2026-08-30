"""Guards for the carrier-vs-bypass test (docs/HANDOFF_carrier_vs_bypass.md).

What is pinned here is chosen deliberately. The headline RMSEs are NOT asserted:
they are a-priori-limited and depend on drivers, and pinning them would freeze a
number rather than a finding. What is pinned is the three things a later session
could plausibly get wrong by re-deriving them:

  1. the carrier has ONE global parameter, not one per congener -- this was
     asserted before being checked, and it was wrong;
  2. the depolarisation arm cannot substitute for an added term, which is the
     result that makes the repo's extension of Trapp justified at all;
  3. the two entry mechanisms are NOT separable on this dataset, so a claim that
     either wins is over-reading.

Plus the API contract that keeps every published number safe: vmax_scale=1 and
g_apo=0 must be exactly the shipped model.
"""
import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "validation"))

import model_api as api                                        # noqa: E402
from conftest import SOLVE_INVARIANCE_REL                       # noqa: E402


# ---------------------------------------------------------------------------
# the API contract
# ---------------------------------------------------------------------------
def test_new_overrides_default_to_the_shipped_model():
    """vmax_scale=1 and g_apo=0 must reach the same solve as passing neither --
    same Compound, same RHS, so this is case 1 in conftest and exact."""
    a = api.simulate("PFOA", season=120.0, n_t=121)
    b = api.simulate("PFOA", season=120.0, n_t=121, vmax_scale=1.0, g_apo=0.0)
    for k in api.TISSUES:
        assert a["baf_final"][k] == b["baf_final"][k], k


def test_pfas_keeps_no_bypass():
    """g_apo is for the neutral/weak-electrolyte path; PFAS is untouched by it."""
    import literature_params as LP
    assert LP.literature_compound("PFOA").g_apo == 0.0


# ---------------------------------------------------------------------------
# the findings
# ---------------------------------------------------------------------------
def test_the_carrier_has_one_global_parameter_not_one_per_congener():
    """Finding that made the comparison fair, and the claim this test exists to
    stop being re-derived: the bypass does NOT win on parsimony, because the
    carrier is a single global Vmax rather than a per-congener fit."""
    assert set(api._CARR) >= {"Vmax_in", "Km_in", "Vmax_out", "Km_out"}
    assert isinstance(api._CARR["Vmax_in"], (int, float))
    for rec in api._CONG.values():
        assert not any(k.lower().startswith("vmax") for k in rec), rec.get("name")


def test_turning_the_carrier_off_opens_a_large_gap():
    """The size of what the carrier is doing -- Trapp's PFAS limit with nothing
    added. Order of magnitude, not a pinned value."""
    on = api.simulate("PFOA", season=120.0, n_t=121)["baf_final"]["root"]
    off = api.simulate("PFOA", season=120.0, n_t=121, vmax_scale=0.0)["baf_final"]["root"]
    assert off < on / 5.0, (on, off)


def test_depolarisation_alone_cannot_replace_the_carrier():
    """Arm C, REFUTED -- the load-bearing result. E_m is already inside Trapp's
    GHK, so if pushing it to the end of its recorded plausible range could stand
    in for the carrier, this repo would not need to have extended him at all. It
    cannot: at -90 mV with the carrier off the root is still far below the
    carrier's."""
    kw = dict(season=120.0, n_t=121)
    carrier = api.simulate("PFOA", **kw)["baf_final"]["root"]
    depol = api.simulate("PFOA", vmax_scale=0.0, E_m_mV=-90.0, **kw)["baf_final"]["root"]
    assert depol < carrier / 3.0, (carrier, depol)
    # and it must still be an improvement on no lever at all, or the arm is broken
    none = api.simulate("PFOA", vmax_scale=0.0, E_m_mV=-120.0, **kw)["baf_final"]["root"]
    assert depol > none


def test_a_bypass_large_enough_to_replace_the_carrier_exists():
    """Arm B reaches the carrier's ballpark -- which is WHY the two are not
    separable. Deliberately a bracket, not an equality: the finding is that they
    are indistinguishable, so asserting a winner either way would contradict it."""
    kw = dict(season=120.0, n_t=121)
    carrier = api.simulate("PFOA", **kw)["baf_final"]["root"]
    bypass = api.simulate("PFOA", vmax_scale=0.0, g_apo=20.0, **kw)["baf_final"]["root"]
    assert 0.2 * carrier < bypass < 5.0 * carrier, (carrier, bypass)


def test_the_bypass_is_not_chain_length_independent():
    """Pre-registered item 1, REFUTED, and the first data contradiction of
    theory_anchor.tex's claim that eta is 'essentially independent of tail
    length'. Reproduced cheaply: the g_apo needed to reach a fixed root BAF must
    differ substantially between a short and a long congener. If a later session
    finds these equal, the refutation has been undone and the doc's claim is back
    in play -- which would be a real finding, not a test to delete."""
    import carrier_vs_bypass as CB

    obs = CB.load_obs()
    grid = [0.5, 2.0, 5.0, 20.0, 50.0, 200.0]

    def best(nm):
        return min(((v, CB.rmse([nm], obs, vmax_scale=0.0, g_apo=v)) for v in grid),
                   key=lambda vr: vr[1])[0]

    short, long = best("PFHxA"), best("PFDoDA")
    assert long > short, (short, long)
    assert long / short >= 5.0, (short, long)
