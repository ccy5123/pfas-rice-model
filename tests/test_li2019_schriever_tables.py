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


def test_partition_bias_is_flat_in_exposure_and_absent_at_low_kow():
    """The part of section 3b that is ROBUST, pinned. The competing explanation
    for the Li 2019 offset is non-equilibrium (Li et al.'s own alpha_pt), which
    predicts the bias shrinks with exposure time; it does not. And the bias is
    absent below log Kow 2, where the water floor dominates, but large above --
    which is what puts the deficit in the sorption term.

    NOT pinned here, deliberately: strict monotonicity across all four bins. It
    holds on the raw rows but not once one study's replicates are collapsed
    (section 3c), so asserting it would freeze an over-claim into the suite."""
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
    assert by_kow[0] > -0.1                   # no bias where the water floor rules
    assert all(b < -0.25 for b in by_kow[1:])  # large wherever the lipid term rules

    # 3c: collapse each compound x study to one row -- the top two bins go flat,
    # so the ladder is not monotone under a fair count. Pinned so the weaker,
    # correct reading is the one the suite defends.
    seen, reps = set(), []
    for r in rows:
        key = (r["compound"], r["source_study"])
        if key not in seen:
            seen.add(key)
            reps.append(r)
    top = [bias([r for r in reps if lo <= float(r["log_kow"]) < hi])
           for lo, hi in ((3.5, 4.5), (4.5, 9))]
    assert top[1] > top[0]                    # NOT still rising once collapsed

    # 3c(ii): raising L is Kow-dependent by construction, so it is not the wrong
    # shape of fix -- the earlier claim that it was is what this pins against.
    La = ND.TRAPP1994_LIPID_FW["root"] * ND.LIPID_OCTANOL_A
    Lb = ND.BRIGGS_ANCHORED_LIPID_FW["root"] * ND.LIPID_OCTANOL_A
    shift = [np.log10((W + Lb * 10 ** (ND.RCF_SLOPE * k))
                      / (W + La * 10 ** (ND.RCF_SLOPE * k))) for k in (1.0, 5.0)]
    assert shift[0] < 0.1 and shift[1] > 0.35


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


# --------------------------------------------------------------------------
# Kodesova 2019 — the A4 table, and the anchor vote it casts
# --------------------------------------------------------------------------
KODESOVA = os.path.join(ROOT, "data_obs", "neutral_obs_kodesova2019.csv")


def test_kodesova_rcf_is_reproducible_from_the_published_numbers():
    """Every value must be recomputable from the three measured quantities its
    own note records, so the file stays a derivation and not a transcription of
    someone's arithmetic. RCF = C_root / (C_soil/KF)^n."""
    import re

    KF = {"HCh": 3.86, "HCa": 2.97, "AE": 0.71}
    n = 1.13
    rows = _rows(KODESOVA)
    assert len(rows) == 21          # 4 plants x 3 soils x 2 treatments, minus 3 NA
    for r in rows:
        root = float(re.search(r"root (\d+) ng/g", r["note"]).group(1))
        soil = float(re.search(r"soil (\d+) ng/g", r["note"]).group(1))
        cw = (soil / 1000.0 / KF[r["soil"]]) ** n
        assert float(r["pore_water_mgL"]) == pytest.approx(cw, rel=1e-3)
        assert float(r["value"]) == pytest.approx(root / 1000.0 / cw, rel=1e-3)
        assert r["basis"] == "dw"   # the article states dry weight explicitly
        assert float(r["log_kow"]) == 2.25


def test_kodesova_isotherm_reading_gives_a_literature_koc():
    """The one pivotal assumption, pinned. Reading KF at mg/kg per mg/L makes it
    the distribution ratio at c = 1 mg/L, so KF/Cox is a Koc -- and it must land
    in carbamazepine's literature band across three very different soils. If it
    ever does not, the derived exposure is wrong and every value above with it."""
    import kodesova2019_carbamazepine as K

    kocs = [K.KF[s] / K.COX[s] for s in K.KF]
    assert all(100 <= k <= 500 for k in kocs)
    assert max(kocs) / min(kocs) < 2.0          # consistent despite 3.8x in Cox
    # derived pore water must sit on the scale of the applied solution (~1 mg/L)
    cw = [float(r["pore_water_mgL"]) for r in _rows(KODESOVA)]
    assert 0.05 < min(cw) and max(cw) < 2.0


def test_kodesova_votes_against_restoring_the_anchor():
    """The decision content. Carbamazepine at log Kow 2.25 separates the two
    compositions by ~1.6x, so this is a well-conditioned vote -- and it goes the
    opposite way to Li 2019. Both tables are pinned here so a later default
    change has to confront the disagreement rather than inherit one side."""
    import neutral_dpu_validation as V

    drv = V.drivers()
    orig = ND.TRAPP1994_LIPID_FW["root"]
    try:
        shipped = [V.compare_to_obs(p, drv, quiet=True) for p in (KODESOVA, LI2019)]
        ND.TRAPP1994_LIPID_FW["root"] = ND.BRIGGS_ANCHORED_LIPID_FW["root"]
        anchored = [V.compare_to_obs(p, drv, quiet=True) for p in (KODESOVA, LI2019)]
    finally:
        ND.TRAPP1994_LIPID_FW["root"] = orig
    assert shipped[0] < anchored[0]          # Kodesova prefers the shipped lipid
    assert anchored[1] < shipped[1]          # Li 2019 prefers the anchor
    assert shipped[0] < 0.25                 # and it is the best a-priori root fit


# --------------------------------------------------------------------------
# the ODE-vs-equilibrium scoring artifact
# --------------------------------------------------------------------------
def test_the_rice_ode_discounts_the_root_in_a_kow_dependent_way():
    """Why `mode="equilibrium"` exists. Running a 120-day rice season to score a
    24-hour barley root measurement imposes a Kow-DEPENDENT discount that has
    nothing to do with the partition being tested: near zero at the extremes,
    but ~0.25 log where the xylem drains the root hardest, around the TSCF peak.
    Kodesova sits at log Kow 2.25, close to the worst point -- so its ODE-basis
    score is flattered. Pinned so the artifact cannot be forgotten again."""
    import neutral_dpu_validation as V

    drv = V.drivers()
    disc = {}
    for lk in (-0.5, 1.78, 2.25, 5.0):
        m = ND.simulate_neutral(ND.NeutralCompound("x", lk), drv)
        disc[lk] = m["baf_final"]["root"] / V.equilibrium_rcf(lk)
    assert disc[-0.5] > 0.85 and disc[5.0] > 0.95      # ~no discount at the ends
    assert disc[1.78] < 0.60                            # ~0.25 log at the peak
    assert disc[2.25] < 0.65                            # and Kodesova sits there


def test_equilibrium_mode_changes_the_root_tables_and_not_ge2017():
    """The appropriate-basis numbers, pinned. Removing the artifact IMPROVES the
    two short-exposure hydroponic tables and WORSENS Kodesova -- which is the
    direction that matters, because it means Kodesova was not in fact the best
    a-priori result in the repo once compared on the right basis."""
    import neutral_dpu_validation as V

    drv = V.drivers()
    liu = os.path.join(ROOT, "data_obs", "neutral_obs_liu2023.csv")
    kod = os.path.join(ROOT, "data_obs", "neutral_obs_kodesova2019.csv")
    ge = os.path.join(ROOT, "data_obs", "neutral_obs_ge2017.csv")

    def both(p):
        return (V.compare_to_obs(p, drv, quiet=True),
                V.compare_to_obs(p, drv, quiet=True, mode="equilibrium"))

    liu_ode, liu_eq = both(liu)
    li_ode, li_eq = both(LI2019)
    kod_ode, kod_eq = both(kod)
    assert liu_eq < liu_ode and li_eq < li_ode        # short-exposure: improves
    assert kod_eq > kod_ode                            # Kodesova was flattered
    assert liu_eq < kod_eq                             # Liu is the best, not Kodesova

    # Ge 2017 carries only stem/leaf tf rows, so equilibrium mode cannot touch it
    assert both(ge)[0] == both(ge)[1]

    # the default must stay "ode", so no published number moves silently
    assert V.compare_to_obs(liu, drv, quiet=True) == liu_ode


def test_equilibrium_mode_does_not_change_the_anchor_verdicts():
    """The anchor tally survives removing the artifact -- and the margins widen,
    which is the useful part: Kodesova's preference for the shipped lipid goes
    from marginal to decisive, and Li 2019's for the anchor likewise."""
    import neutral_dpu_validation as V

    drv = V.drivers()
    kod = os.path.join(ROOT, "data_obs", "neutral_obs_kodesova2019.csv")
    orig = ND.TRAPP1994_LIPID_FW["root"]
    try:
        ship = [V.compare_to_obs(p, drv, quiet=True, mode="equilibrium")
                for p in (kod, LI2019)]
        ND.TRAPP1994_LIPID_FW["root"] = ND.BRIGGS_ANCHORED_LIPID_FW["root"]
        anch = [V.compare_to_obs(p, drv, quiet=True, mode="equilibrium")
                for p in (kod, LI2019)]
    finally:
        ND.TRAPP1994_LIPID_FW["root"] = orig
    assert ship[0] < anch[0] and anch[1] < ship[1]     # same opposite verdicts
    assert anch[0] - ship[0] > 0.10                    # but a much wider margin


def test_spearman_helper_handles_ties():
    import schriever2020_tscf as S

    assert S.spearman([1, 2, 3, 4], [1, 2, 3, 4]) == pytest.approx(1.0)
    assert S.spearman([1, 2, 3, 4], [4, 3, 2, 1]) == pytest.approx(-1.0)
    # ties must be averaged, not left in arbitrary encounter order
    assert S.spearman([1, 1, 2, 2], [1, 1, 2, 2]) == pytest.approx(1.0)
    assert np.isnan(S.spearman([1, 1, 1], [1, 2, 3]))
