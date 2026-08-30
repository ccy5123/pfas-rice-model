"""Guards for Briggs 1983 Table 1 -- the repo's only MEASURED stem test.

Three things need pinning here, and they are pinned separately on purpose:

  * the TRANSCRIPTION. Every (compound x harvest) row of Table 1 sums to 100 %
    across the three sections, which is what caught two digits the PDF's text
    layer got wrong. If a re-parse ever drops or mangles a row, this fails
    before any science does.
  * the RECONSTRUCTION. The Stem Concentration Factor is derived here (% x total
    dpm / section weight / solution concentration), not published as such, so
    the claim that it IS the quantity the article plotted rests on reproducing
    the article's own eq. (3) with it. That check is the load-bearing one.
  * the RESULT. The stem a-priori RMSE, and -- just as important -- the two
    NEGATIVE findings: the leaf is not an equilibrium, and the stem/leaf split
    moves between harvests, which is why this table cannot do the job the
    handoff wanted it for. A later session that quietly starts treating the
    split as a fixed model target is going backwards.
"""
from __future__ import annotations

import csv
import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "validation"))

import neutral_dpu as ND  # noqa: E402
import briggs1983_shoot as B  # noqa: E402

OBS = os.path.join(ROOT, "data_obs", "neutral_obs_briggs1983_shoot.csv")


@pytest.fixture(scope="module")
def rows():
    return B.load(OBS)


def test_table_is_intact(rows):
    assert len(rows) == 96                                  # 16 x 2 harvests x 3
    compounds = {r["compound"] for r in rows}
    assert len(compounds) == 16
    assert {r["hours"] for r in rows} == {24, 48}
    assert {r["tissue"] for r in rows} == {"leaf", "stem_central", "stem_base"}
    lk = [r["log_kow"] for r in rows]
    assert min(lk) == pytest.approx(-0.57) and max(lk) == pytest.approx(3.7)
    # the two series the article separates
    assert {r["series"] for r in rows} == {"oxime", "phenylurea"}
    assert sum(1 for r in rows if r["series"] == "oxime") == 7 * 2 * 3


def test_every_row_sums_to_one_hundred_percent(rows):
    """The transcription's own check -- this is what caught 17.4 (should be
    77.4) and 3.1 (should be 3.7) in the PDF text layer."""
    by = {}
    for r in rows:
        by.setdefault((r["compound"], r["hours"]), []).append(r["pct_of_shoot"])
    assert len(by) == 32
    for key, pcts in by.items():
        assert len(pcts) == 3
        assert sum(pcts) == pytest.approx(100.0, abs=0.15), key


def test_reconstructed_scf_reproduces_the_articles_own_equation(rows):
    """THE load-bearing check. The SCF is derived, so the evidence that it is
    the article's Fig. 4 quantity is that it lands on the article's eq. (3) --
    which the paper says its own points 'fit quite well'."""
    n, bias, rmse = B._score(rows, "stem_central", B.briggs_eq3)
    assert n == 30                                     # 15 equilibrated x 2 harvests
    assert abs(bias) < 0.10
    assert rmse < 0.35


def test_stem_base_runs_high_exactly_as_the_article_says(rows):
    """The article ascribes the base running 'up to twice' the central section
    partly to direct contact with the treating solution. If the reconstruction
    were wrong this ordering would not survive."""
    _, bias_c, _ = B._score(rows, "stem_central", B.briggs_eq3)
    _, bias_b, _ = B._score(rows, "stem_base", B.briggs_eq3)
    assert bias_b > bias_c + 0.10


def test_the_stem_a_priori_result(rows):
    """The headline: the repo's stem, nothing fitted, against measurement -- and
    against the paper's own FITTED equation, which had this data to fit."""
    n, bias, rmse = B._score(rows, "stem_central", B.repo_scf)
    assert n == 30
    assert rmse == pytest.approx(0.299, abs=0.02)
    assert abs(bias) < 0.10
    # indistinguishable from the fitted equation -- section 4j's "the 4.1x
    # coefficient gap cancels in the observable", now measured rather than argued
    _, _, rmse_briggs = B._score(rows, "stem_central", B.briggs_eq3)
    assert abs(rmse - rmse_briggs) < 0.10


def test_the_stem_result_does_not_depend_on_the_open_root_anchor():
    """Worth its own guard: the lipid_source question is about the ROOT, so this
    new stem result is independent of how that is eventually decided. If a later
    change makes the anchor touch the stem, this catches it."""
    assert (ND.LIPID_SOURCES["briggs_anchor"]["stem"]
            == ND.LIPID_SOURCES["measured"]["stem"])
    for lk in (0.0, 1.78, 3.0):
        a = ND.k_pw(lk, W=ND.RICE_WATER["stem"],
                    L=ND.LIPID_SOURCES["briggs_anchor"]["stem"])
        m = ND.k_pw(lk, W=ND.RICE_WATER["stem"],
                    L=ND.LIPID_SOURCES["measured"]["stem"])
        assert a == pytest.approx(m)


def test_the_leaf_is_not_an_equilibrium_and_the_bias_grows(rows):
    """NEGATIVE finding, pinned. Scored as an equilibrium the leaf drifts, and
    the drift GROWS with exposure -- the terminal-accumulator signature. This is
    why the leaf half of this table is not a partition test."""
    def leaf_eq(x):
        return ND.k_pw(x, W=ND.RICE_WATER["leaf"],
                       L=ND.TRAPP1994_LIPID_FW["leaf"],
                       b=ND.LEAF_SORPTION_EXPONENT) * ND.tscf(x)

    _, b24, _ = B._score(rows, "leaf", leaf_eq, hours=24)
    _, b48, _ = B._score(rows, "leaf", leaf_eq, hours=48)
    assert b48 > b24 + 0.20            # climbing away from equilibrium
    assert b24 > 0.0                   # and already above it at 24 h


def test_the_stem_leaf_split_moves_between_harvests(rows):
    """NEGATIVE finding, pinned. The handoff wanted this table to constrain the
    stem/leaf split. It cannot: the split is a function of exposure duration,
    not of the compound, so a single number from it is not a model target."""
    by = {}
    for r in rows:
        by.setdefault((r["compound"], r["hours"]), {})[r["tissue"]] = r
    ratios = []
    for c in {c for c, _ in by}:
        share = {h: by[(c, h)]["stem_central"]["pct_of_shoot"]
                    + by[(c, h)]["stem_base"]["pct_of_shoot"] for h in (24, 48)}
        ratios.append(share[24] / share[48])
    assert np.median(ratios) > 1.3          # stem share drops as the leaf fills
    assert max(ratios) > 1.9                # and by up to ~2x for some compounds


def test_the_flagged_compounds_are_marked_not_averaged_in(rows):
    """The article names compounds that deviate for stated reasons. They must be
    excluded (equilibrated=0) or labelled (metabolised=1), never quietly scored."""
    flags = {r["compound"]: (r["equilibrated"], r["metabolised"]) for r in rows}
    assert flags["4_4_bromophenoxyphenylurea"][0] == 0     # never equilibrated
    assert flags["aldicarb"][1] == 1                       # parent + oxidation products
    assert flags["3_methylthiophenylurea"][1] == 1
    # and the excluded one really is excluded from the scored set
    n_all = sum(1 for r in rows if r["tissue"] == "stem_central")
    n_scored, _, _ = B._score(rows, "stem_central", B.repo_scf)
    assert n_all == 32 and n_scored == 30
