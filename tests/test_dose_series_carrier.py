"""Guards for the dose-series test (validation/dose_series_carrier.py, handoff B1).

What is pinned is the STRUCTURE the finding rests on, not the RMSEs and not the
verdict's wording. The three things a later session could plausibly get wrong:

  1. that the carrier is the model's only nonlinearity in exposure -- if it were
     not, "BCF falls with dose" would stop being a carrier signature and the
     whole discriminator would be invalid;
  2. that entry magnitude divides out of TF -- the claim that lets a TF trend be
     read as NOT the entry term;
  3. that only ONE of Tang's three congeners is well-conditioned for this test,
     which is why the pre-registered 2-of-3 rule was not met and why the result
     is "disfavoured", not "refuted". A later session that reads a clean
     refutation out of this file has over-read it.

Plus the API contract: km_scale=1 must be the shipped model.
"""
import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "validation"))

import model_api as api                                          # noqa: E402
import dose_series_carrier as DS                                 # noqa: E402


# ---------------------------------------------------------------------------
# the API contract
# ---------------------------------------------------------------------------
def test_km_scale_defaults_to_the_shipped_carrier():
    a = api.simulate("PFOA", season=120.0, n_t=121)
    b = api.simulate("PFOA", season=120.0, n_t=121, km_scale=1.0)
    for k in api.TISSUES:
        assert a["baf_final"][k] == b["baf_final"][k], k
    assert b["params"]["Km_in"] == api._CARR["Km_in"]


def test_km_scale_makes_the_carrier_more_linear_not_merely_stronger():
    """What section 7's bound actually needs, and NOT what a first pass at this
    test asserted: raising Km does not raise the flux at a fixed exposure (at
    Cwo >> Km it lowers it, since the term sits further below Vmax). It flattens
    the DOSE RESPONSE -- the BAF ratio between two exposures goes to 1 -- which
    is the whole basis for reading a large Km as "linear over the span"."""
    kw = dict(season=120.0, n_t=121, f_xy_source="recommended")
    ratio = lambda ks: (api.simulate("PFOS", Cwo=10.0, km_scale=ks, **kw)["baf_final"]["root"]
                        / api.simulate("PFOS", Cwo=10000.0, km_scale=ks, **kw)["baf_final"]["root"])
    tight, loose = ratio(1.0), ratio(1000.0)
    assert tight > 3.0, tight                     # saturating: BAF falls with dose
    assert loose == pytest.approx(1.0, abs=0.1), loose      # linear over the span


# ---------------------------------------------------------------------------
# the structure the finding rests on
# ---------------------------------------------------------------------------
def test_the_carrier_is_the_only_nonlinearity_in_exposure():
    """Finding 1. With the carrier off the model must be EXACTLY linear in Cwo
    (GHK and the bypass are both linear), so any curvature in a dose series is
    the carrier's. If this fails, the dose series stops being a clean structural
    discriminator and the B1 result must be re-derived, not repaired."""
    kw = dict(uptake="bypass", season=120.0, n_t=121, f_xy_source="recommended")
    lo = api.simulate("PFOA", Cwo=1.0, **kw)["baf_final"]["root"]
    hi = api.simulate("PFOA", Cwo=1000.0, **kw)["baf_final"]["root"]
    assert lo == pytest.approx(hi, rel=1e-3), (lo, hi)
    # and with the carrier ON it must NOT be linear, or there is nothing to test
    con = dict(season=120.0, n_t=121, f_xy_source="recommended")
    clo = api.simulate("PFOA", Cwo=1.0, **con)["baf_final"]["root"]
    chi = api.simulate("PFOA", Cwo=1000.0, **con)["baf_final"]["root"]
    assert chi < 0.9 * clo, (clo, chi)


def test_entry_magnitude_divides_out_of_the_transfer_factor():
    """Finding 2, and the reason a TF trend can be attributed to something other
    than entry: every downstream term is linear in C, so the shoot/root ratio is
    ~invariant to how much came in."""
    lo = DS.model_tf("PFOA", 1.0)
    hi = DS.model_tf("PFOA", 10000.0)
    assert hi == pytest.approx(lo, rel=0.02), (lo, hi)


def test_only_one_congener_is_well_conditioned_for_this_series():
    """Finding 3 -- the guard against over-reading the verdict. GenX sits far
    above Km at even the lowest dose, so carrier and bypass predict nearly the
    same thing for it; PFOS is the one that crosses Km inside the series."""
    km = api._CARR["Km_in"]
    ratio = {nm: DS.pore_water(nm, 0.1) / km for nm in DS.CONGENERS}
    assert ratio["GenX"] > 100.0, ratio           # saturated throughout -> no signal
    assert ratio["PFOS"] < 10.0, ratio            # crosses Km inside the span
    assert ratio["PFOS"] < ratio["PFOA"] < ratio["GenX"], ratio


def test_the_gate_passes_but_not_for_every_congener():
    """The pre-registration's own flaw, kept visible: the gate was written across
    congeners and passes, while GenX individually fails it. A later session
    should see this rather than rediscover it."""
    doses = [0.1, 1.0, 10.0, 50.0, 100.0]
    dec = {}
    for nm in ("PFOS", "GenX"):
        cw = [DS.pore_water(nm, d) for d in doses]
        dec[nm] = DS.decline([DS.whole_plant_baf(nm, c, "carrier") for c in cw])
    assert dec["PFOS"] >= DS.GATE_MIN_DECLINE, dec       # gate passes overall
    assert dec["GenX"] < DS.GATE_MIN_DECLINE, dec        # but not for GenX


def test_the_measured_series_are_far_flatter_than_a_saturating_carrier():
    """The result itself, as a direction rather than a number: Tang's measured
    BCF barely moves over a 1000x dose span. Bracketed loosely so it survives
    driver changes -- the finding is the SIZE of the mismatch, not its value."""
    data = DS.load_doses()
    for nm in DS.CONGENERS:
        obs = DS.decline([v for _, v, _ in data[nm]["BCF"]])
        assert obs < 2.0, (nm, obs)
    cw = [DS.pore_water("PFOS", d) for d in (0.1, 1.0, 10.0, 50.0, 100.0)]
    carrier = DS.decline([DS.whole_plant_baf("PFOS", c, "carrier") for c in cw])
    obs = DS.decline([v for _, v, _ in data["PFOS"]["BCF"]])
    assert carrier > 3.0 * obs, (carrier, obs)


def test_toxicity_signature_translocation_moves_more_than_uptake():
    """Pre-registered item (b) as measured: for PFOA the TF endpoints fall MORE
    than the BCF does. Entry cannot produce that (see the TF test above), so the
    dose response is not located where a saturating carrier would put it."""
    data = DS.load_doses()
    bcf = DS.decline([v for _, v, _ in data["PFOA"]["BCF"]])
    for e in ("TF_stalk", "TF_leaf", "TF_endosperm"):
        tf = DS.decline([v for _, v, _ in data["PFOA"][e]])
        assert tf > bcf, (e, tf, bcf)


def test_nothing_was_adopted_by_the_dose_series():
    """Km is not re-fitted on a bound from three congeners in one soil."""
    assert api._CARR["Km_in"] == 5.0
    assert api.UPTAKE_MODES[api.DEFAULT_UPTAKE]["vmax_scale"] == 1.0
