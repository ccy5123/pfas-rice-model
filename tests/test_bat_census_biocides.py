"""
Tests for the BAT-census biocide application -- validation/bat_census_biocides.py.

Like the Hwang 2017 tests, most of what matters here is what the work REFUSES to
conclude. These substances have no measured rice endpoint, so every number the
script prints is a prediction; the guards below pin (a) the transcription checks
that stand between the report's table and this model's inputs, (b) the scope
screen that decides what may be quoted, and (c) the caveats, so a later session
cannot quietly promote a prediction into a validation result.

The one genuinely scoreable thing -- propiconazole and triclosan, the two
substances in the report for which this repo already holds a measured root row --
is pinned too, because it is a-priori (log Kow in, partition out, nothing fitted)
and a drift in `k_pw` or in the shipped tissue composition would move it.
"""
import csv
import math
import os
import sys

import pytest

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "src"))
sys.path.insert(0, os.path.join(_ROOT, "validation"))

import bat_census_biocides as BC  # noqa: E402


@pytest.fixture(scope="module")
def rows():
    return BC.load_table()


@pytest.fixture(scope="module")
def screened(rows):
    return BC.screen(rows, quiet=True)


# --------------------------------------------------------------------------
# the transcription: the report's numbers, and the one inversion applied to them
# --------------------------------------------------------------------------
def test_report_percent_ionised_is_on_ph7(rows):
    """The load-bearing check. The report tabulates a percent ionised for the
    anticoagulants and a pKa for three other substances; pKa is recovered from
    the former only because the latter round-trip at pH 7. If this failed, every
    derived pKa -- and so the whole ionisation screen -- would be off in the same
    direction and nothing else in the file would notice."""
    st = BC.check_transcription(rows, quiet=True)
    assert st["n_stated"] >= 3
    assert st["worst_stated"] < 0.01          # percentage points
    assert st["worst_derived"] < 1e-3


def test_derived_pka_is_the_inversion_and_nothing_else(rows):
    for r in rows:
        if (r.get("pka_basis") or "").strip() != "derived":
            continue
        pct = BC._f(r, "pct_ionised_pH7")
        assert pct is not None, r["substance"]
        assert BC._f(r, "pKa") == pytest.approx(
            BC.pka_from_fraction(pct / 100.0), abs=1e-3), r["substance"]


def test_hh_helpers_are_inverses():
    for pka in (2.0, 4.5, 7.0, 8.065, 11.0):
        for acid in (True, False):
            a = BC.frac_ionised(pka, BC.REPORT_PH, acid)
            assert BC.pka_from_fraction(a, BC.REPORT_PH, acid) == pytest.approx(pka, abs=1e-9)
    # a base is ionised BELOW its pKa, an acid above it -- the sign must not drift
    assert BC.frac_ionised(9.0, 7.0, is_acid=True) < 0.5
    assert BC.frac_ionised(9.0, 7.0, is_acid=False) > 0.5


def test_no_input_is_invented(rows):
    """The report's own section 8.14 defect -- entering a number no source
    reports -- must not be reproduced here. A row either carries a log Kow with a
    report section behind it, or carries no log Kow at all."""
    for r in rows:
        if r["status"] == "run":
            assert BC._f(r, "log_kow") is not None, r["substance"]
            assert r["section"].strip(), r["substance"]
        elif r["status"] == "no_logkow_in_report":
            assert (r["log_kow"] or "").strip() == "", r["substance"]
    assert len([r for r in rows if r["status"] == "repo_logkow"]) == 1


def test_table_is_not_an_observation_file(rows):
    """It must stay unusable by the `--obs` harness. That harness would run every
    row on the rice drivers and return a number, and there is no measured plant
    value here for it to be a number ABOUT."""
    forbidden = {"tissue", "value", "endpoint", "basis"}
    assert not (set(rows[0]) & forbidden)


# --------------------------------------------------------------------------
# the scope screen
# --------------------------------------------------------------------------
def test_screen_excludes_exactly_the_strongly_ionised(screened):
    """The report's section 3.0a rule -- 'more than 90% ionised at environmental
    pH' -- which its own screen defined and never implemented. It is implemented
    here at the same threshold, and this pins that it fires on the group the
    report says it would have caught."""
    excluded = {o["substance"] for o in screened if o["verdict"].startswith("EXCLUDED")}
    assert {"Brodifacoum", "Difenacoum", "Bromadiolone", "Flocoumafen"} <= excluded
    for o in screened:
        assert (o["verdict"].startswith("EXCLUDED")) == (
            o["pKa"] is not None and o["f_n"] < BC.F_N_FLOOR), o["substance"]


def test_screen_spans_come_from_the_shipped_tables(screened):
    """Not hardcoded: the bands are recomputed from data_obs/ so that adding a
    measured table widens the scope automatically, and removing one narrows it."""
    spans = BC.measured_spans()
    p_lo, p_hi, n_p = spans["partition"]
    t_lo, t_hi, n_t = spans["tscf"]
    assert n_p > 100 and n_t == 30
    assert p_lo < 0.0 < p_hi and p_hi > 8.0          # the soil table reaches log Kow 8.7
    assert t_hi < 6.0                                # TSCF is measured nowhere above ~5.5
    # the substances the opinions call bioaccumulative all sit past the TSCF anchor
    for o in screened:
        if o["substance"] in ("Difethialone", "Flocoumafen", "Hexaflumuron"):
            assert o["log_kow"] > t_hi - 0.5


def test_muscalure_and_cholecalciferol_are_flagged_extrapolation(screened):
    far = {o["substance"] for o in screened if o["verdict"] == "extrapolation"}
    assert far == {"Cholecalciferol", "Cis-tricos-9-ene (muscalure)"}


# --------------------------------------------------------------------------
# the two substances that can actually be scored -- a priori, nothing fitted
# --------------------------------------------------------------------------
def test_overlap_apriori_predictions(tmp_path):
    ov = BC.overlap_check(quiet=True)
    p = ov["propiconazole"]
    assert p["log_kow"] == 3.72
    assert p["obs"] == pytest.approx(9.323, rel=1e-6)
    # a-priori: log Kow in, K_PW out. Within a factor of ~1.1 of the measured rice root.
    assert abs(p["err"]) < 0.10
    t = ov["triclosan"]
    assert t["n"] == 14
    # the model runs HIGH on the Li 2019 soil table (docs section 4h) and higher
    # still on this acid; the SIGN is the finding, not the size
    assert t["err"] > 0.2
    # and the two independent log Kow sourcings agree: Li 2019's 4.8 against the
    # report's audited 4.76
    assert abs(t["log_kow"] - 4.76) < 0.10


# --------------------------------------------------------------------------
# the cross-model claims that must not silently break
# --------------------------------------------------------------------------
def test_mw_solubility_henry_are_structurally_inert():
    """The report needed 217 runs to establish that BAT's answer does not depend
    on molecular weight, solubility or Henry's constant. Here it is a property of
    the equations, and this test is what keeps it one: two substances at the same
    log Kow must return the SAME number, exactly."""
    a = BC.run_one(6.29, "cyphenothrin")
    b = BC.run_one(6.29, "difethialone")
    for k in ("TSCF", "K_PW_root", "root", "stem", "leaf", "grain", "straw"):
        assert a[k] == b[k], k


def test_air_is_the_one_term_that_separates_them():
    """...and the exception, which is the interesting half: Henry's constant DOES
    reach the answer, through the opt-in air term. Empenthrin (K_AW 1.4e-2) must
    lose leaf burden to volatilisation where cyphenothrin (9e-7) does not."""
    RT = 8.314 * 293.15
    volatile = BC.run_one(6.30, "empenthrin", MW=274.4, K_AW=34.65 / RT, air=True)
    quiet = BC.run_one(6.30, "empenthrin", MW=274.4, K_AW=34.65 / RT, air=False)
    assert volatile["leaf"] < 0.75 * quiet["leaf"]
    involatile_on = BC.run_one(6.29, "cyphenothrin", MW=375.5, K_AW=0.0022 / RT, air=True)
    involatile_off = BC.run_one(6.29, "cyphenothrin", MW=375.5, K_AW=0.0022 / RT, air=False)
    assert involatile_on["leaf"] == pytest.approx(involatile_off["leaf"], rel=1e-3)
    # the root carries NO air term (the derivation excludes root volatilisation),
    # so what little it moves is the coupled solve, not a root flux: <0.01%,
    # against the leaf's tens of percent
    assert volatile["root"] == pytest.approx(quiet["root"], rel=1e-4)


def test_ionisation_is_more_damped_here_than_in_bat():
    """The report's section 8.15 sweep, re-run. At 91% ionised the chemistry says
    0.091 and BAT says 0.245; this model's root says ~0.77, because the membrane
    term sets the RATE of root equilibration and not its level. That is the
    stated reason the strongly-ionised group is excluded rather than reported,
    so it has to stay true for the exclusion to keep making sense."""
    base = BC.run_one(4.60, "dcpp-neutral")["root"]
    weak = BC.run_one(4.60, "dcpp", pKa=6.0, is_acid=True)["root"]
    f_n = 1.0 - BC.frac_ionised(6.0)
    assert f_n == pytest.approx(0.0909, abs=1e-3)
    rel = weak / base
    assert rel > 940.0 / 3830.0        # less damped than BAT
    assert rel > f_n                   # and far less than Henderson-Hasselbalch
    assert rel < 1.0                   # but it is not inert either


# --------------------------------------------------------------------------
# the caveats
# --------------------------------------------------------------------------
def test_the_refusals_are_in_the_file():
    """These substances have no measured rice endpoint. The script says so in as
    many words, and a later session must not be able to delete that quietly while
    keeping the numbers."""
    src = open(os.path.join(_ROOT, "validation", "bat_census_biocides.py")).read()
    for phrase in ("PREDICTION, NOT VALIDATION",
                   "NOT a validation",
                   "SENSITIVITY ONLY",
                   "UPPER BOUND"):
        assert phrase in src, phrase
    table = open(os.path.join(_ROOT, "data_obs", "biocides_bat_census.csv")).read()
    assert "WHAT IT IS NOT" in table
    assert "FISH whole-body biotransformation half-life" in table


def test_fish_km_is_never_the_default(rows, screened):
    """`half_life` must stay opt-in: the headline table is gamma = 0. A fish rate
    silently becoming the plant default is the failure this whole framing is
    guarding against."""
    src = open(os.path.join(_ROOT, "validation", "bat_census_biocides.py")).read()
    assert "def run_one(log_kow, name, half_life=None" in src
    r = BC.run_one(4.60, "x")
    assert r["res"]["params"]["gamma"] == 0.0
