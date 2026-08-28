"""Guards for the two measured tables that arrived with the A1/A3 papers.

These lock down three things that a later session could otherwise undo quietly:

  * the transcriptions themselves (row counts, species tallies, ranges), so a
    re-parse that silently drops rows is caught;
  * the `subset` mechanism, whose entire job is to keep a table's CALIBRATION
    rows out of an a-priori score -- and, just as important, to leave the two
    tables WITHOUT that column (liu2023, ge2017) untouched, so the published
    0.281 / 0.783 cannot move because of it;
  * the root-partition anchor discrepancy that validation/li2019_rcf_apriori.py
    diagnoses. That one is deliberately NOT fixed, which is exactly why it needs
    a test: an undocumented drift back toward "consistent" would erase a
    recorded open question.
"""
from __future__ import annotations

import csv
import os
import sys
from collections import Counter

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "validation"))

import neutral_dpu as ND  # noqa: E402

LI2019 = os.path.join(ROOT, "data_obs", "neutral_obs_li2019_rcf.csv")
SCHRIEVER = os.path.join(ROOT, "data_obs", "tscf_obs_schriever2020.csv")


def _rows(path):
    with open(path, newline="") as f:
        return [r for r in csv.DictReader(x for x in f
                                          if not x.lstrip().startswith("#"))
                if any(v for v in r.values() if isinstance(v, str) and v.strip())]


# --------------------------------------------------------------------------
# the transcriptions
# --------------------------------------------------------------------------
def test_li2019_table_is_intact():
    rows = _rows(LI2019)
    # 48 rows in Table S1 minus clotrimazole, which the header excludes as a
    # weak base (pKa 6.02) that is substantially cationic at hydroponic pH.
    assert len(rows) == 47
    subs = Counter(r["subset"] for r in rows)
    assert subs == {"apriori": 29, "calibration": 18}
    # every calibration row is Briggs' barley -- that is what makes it in-sample
    assert {r["species"] for r in rows if r["subset"] == "calibration"} == {"Barley"}
    lk = [float(r["log_kow"]) for r in rows]
    assert min(lk) == pytest.approx(-0.57) and max(lk) == pytest.approx(5.41)
    # value = 10^(log RCF); the header records the log in the note column
    for r in rows:
        log_rcf = float(r["note"].split("10^")[1])
        assert float(r["value"]) == pytest.approx(10.0 ** log_rcf, rel=1e-3)
    assert sum(1 for r in rows if r["species"].startswith("Rice")) == 4


def test_schriever_table_matches_its_own_species_tally():
    rows = _rows(SCHRIEVER)
    assert len(rows) == 97                      # the SI says n = 97
    by_species = Counter(r["species"] for r in rows)
    # Table A 1 of the same SI: barley 51, poplar 15, cattails 8, wheat 6,
    # tomato 4, potato 3, black willow 2. Parsing that reproduces the paper's
    # own tally is the check that no row was dropped or duplicated.
    assert by_species["Barley"] == 51
    assert by_species["Poplar Hybrid"] + by_species["Hybrid poplar"] == 15
    assert by_species["Common cattails"] == 8
    assert by_species["Wheat"] == 6 and by_species["Tomato"] == 4
    assert by_species["Potato"] == 3 and by_species["Black willow"] == 2
    tscf = [float(r["TSCF"]) for r in rows]
    assert 0.0 <= min(tscf) and max(tscf) <= 1.0
    # 999 is the SI's sentinel and must survive as-is, not be blanked
    assert any(float(r["pKa"]) == 999 for r in rows)


def test_schriever_neutral_flag_is_what_it_claims():
    """`neutral_at_test_pH` must be |logD(test pH) - logP| < 0.1, with the
    documented fall-back to pH 5.5 when no test pH was reported."""
    for r in _rows(SCHRIEVER):
        ref = (float(r["logD_test_pH"]) if float(r["pH_test"]) != 999
               else float(r["logD_pH5.5"]))
        assert int(r["neutral_at_test_pH"]) == int(abs(ref - float(r["logP"])) < 0.1)


# --------------------------------------------------------------------------
# the subset mechanism must not touch the existing tables
# --------------------------------------------------------------------------
def test_subset_filter_leaves_the_published_tables_alone():
    import neutral_dpu_validation as V

    drv = V.drivers()
    for name in ("neutral_obs_liu2023.csv", "neutral_obs_ge2017.csv"):
        path = os.path.join(ROOT, "data_obs", name)
        with open(path) as f:
            header = next(x for x in f if x.lstrip().startswith("compound"))
        assert "subset" not in header.split(",")     # no column -> no filtering
        a = V.compare_to_obs(path, drv, quiet=True)
        b = V.compare_to_obs(path, drv, quiet=True, subset=None)
        assert a == b


def test_subset_filter_actually_holds_back_the_calibration_rows():
    import neutral_dpu_validation as V

    drv = V.drivers()
    apriori = V.compare_to_obs(LI2019, drv, quiet=True)              # default
    everything = V.compare_to_obs(LI2019, drv, quiet=True, subset=None)
    insample = V.compare_to_obs(LI2019, drv, quiet=True, subset="calibration")
    assert apriori != everything != insample
    # the held-back rows are Briggs' own barley, so they must score BETTER than
    # the out-of-sample ones -- if they ever do not, the filter is mislabelled
    assert insample < apriori


# --------------------------------------------------------------------------
# the anchor discrepancy: recorded, not fixed
# --------------------------------------------------------------------------
def test_root_lipid_sits_below_the_briggs_anchor_by_the_documented_factor():
    shipped = ND.TRAPP1994_LIPID_FW["root"] * ND.LIPID_OCTANOL_A
    anchor = 10.0 ** ND.RCF_INTERCEPT
    assert anchor / shipped == pytest.approx(2.48, abs=0.02)
    # the opt-in dict restores the product exactly and changes nothing else
    assert (ND.BRIGGS_ANCHORED_LIPID_FW["root"] * ND.LIPID_OCTANOL_A
            == pytest.approx(anchor))
    for organ in ("stem", "leaf", "grain"):
        assert ND.BRIGGS_ANCHORED_LIPID_FW[organ] == ND.TRAPP1994_LIPID_FW[organ]


def test_briggs_anchored_root_reproduces_briggs_rcf_and_the_default_does_not():
    """The point of the discrepancy, as a number: with the anchor the rice root
    IS the Briggs RCF curve above the water floor; with the shipped lipid it is
    2.5x under it. Compared at a fixed log Kow, so nothing else can drift in."""
    log_kow = 4.0
    briggs = ND.briggs_rcf(log_kow)
    W = ND.RICE_WATER["root"]

    def kpw(lipids):
        return ND.k_pw(log_kow, W=W, L=lipids["root"])

    anchored, shipped = kpw(ND.BRIGGS_ANCHORED_LIPID_FW), kpw(ND.TRAPP1994_LIPID_FW)
    # above the floor the anchored root and Briggs' RCF differ only by W vs 0.82
    assert (anchored - W) == pytest.approx(briggs - ND.RCF_FLOOR, rel=1e-9)
    assert (anchored - W) / (shipped - W) == pytest.approx(2.48, abs=0.02)


def test_default_composition_is_untouched_by_the_new_constant():
    """Importing the anchor must not have changed what the model runs."""
    assert ND.TRAPP1994_LIPID_FW == {"root": 0.01, "stem": 0.03,
                                     "leaf": 0.03, "grain": 0.03}
    comps = {c.name: c for c in ND.rice_compartments()}
    assert comps["root"].f_PL == pytest.approx(0.01 / (1.0 - 0.90))


def test_partition_bias_is_flat_in_exposure_and_grows_with_kow():
    """The diagnosis in section 3b, pinned. The competing explanation for the
    Li 2019 offset is non-equilibrium (Li et al.'s own alpha_pt), which predicts
    the bias shrinks with exposure time. It does not, and it grows with log Kow
    instead -- which is what puts the deficit in the sorption term. If a future
    change ever flips either of those, the recorded conclusion is wrong."""
    rows = [r for r in _rows(LI2019) if r["subset"] == "apriori"]
    W = ND.RICE_WATER["root"]
    La = ND.TRAPP1994_LIPID_FW["root"] * ND.LIPID_OCTANOL_A

    def bias(rs):
        return float(np.mean([
            np.log10(W + La * 10.0 ** (ND.RCF_SLOPE * float(r["log_kow"])))
            - np.log10(float(r["value"])) for r in rs]))

    short = bias([r for r in rows if 1 <= float(r["exposure_d"]) < 3])
    long_ = bias([r for r in rows if float(r["exposure_d"]) >= 3])
    assert abs(short - long_) < 0.05          # flat in exposure time

    by_kow = [bias([r for r in rows if lo <= float(r["log_kow"]) < hi])
              for lo, hi in ((-1, 2), (2, 3.5), (3.5, 4.5), (4.5, 9))]
    assert all(b > nxt for b, nxt in zip(by_kow, by_kow[1:]))   # monotone downward
    assert by_kow[0] > -0.1                   # no bias where the water floor rules
    assert by_kow[-1] < -0.5                  # large where the lipid term rules


# --------------------------------------------------------------------------
# the TSCF QSPR test itself
# --------------------------------------------------------------------------
def test_tscf_qspr_underpredicts_the_measured_values():
    """The headline of validation/schriever2020_tscf.py: on the un-ionised rows
    the default Briggs bell is biased LOW -- the same direction as the root
    partition finding, which is why the two are reported together."""
    import schriever2020_tscf as S

    rows = [r for r in S.load(SCHRIEVER) if r["neutral"]]
    assert len(rows) == 30
    briggs = S.score(rows, "briggs", "logP")
    schriever = S.score(rows, "schriever", "logD")
    assert briggs["bias"] < -0.15                    # under-predicts
    assert schriever["bias"] < 0.0                   # so does the fitted one
    # the fitted model does better on its own training set -- if it ever did
    # not, something is wrong with the transcription or the QSPR
    assert schriever["rmse"] < briggs["rmse"]
    assert briggs["rho"] > 0.7                       # but the bell still ranks well


def test_logd_beats_logp_on_the_ionisable_rows():
    """Schriever & Lamshoeft's central claim, checked on their own table: logD at
    the test pH orders the full 97 far better than logP. It is not a claim this
    repo needs -- the neutral path is un-ionised by construction -- but it is a
    free, independent check that the transcription carries real signal."""
    import schriever2020_tscf as S

    rows = S.load(SCHRIEVER)
    assert S.score(rows, "briggs", "logD")["rho"] > \
        S.score(rows, "briggs", "logP")["rho"] + 0.2


def test_spearman_helper_handles_ties():
    import schriever2020_tscf as S

    assert S.spearman([1, 2, 3, 4], [1, 2, 3, 4]) == pytest.approx(1.0)
    assert S.spearman([1, 2, 3, 4], [4, 3, 2, 1]) == pytest.approx(-1.0)
    # ties must be averaged, not left in arbitrary encounter order
    assert S.spearman([1, 1, 2, 2], [1, 1, 2, 2]) == pytest.approx(1.0)
    assert np.isnan(S.spearman([1, 1, 1], [1, 2, 3]))
