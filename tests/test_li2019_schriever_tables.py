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
    shipped = [V.compare_to_obs(p, drv, quiet=True, lipid_source="measured")
               for p in (KODESOVA, LI2019)]
    anchored = [V.compare_to_obs(p, drv, quiet=True, lipid_source="briggs_anchor")
                for p in (KODESOVA, LI2019)]
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
    ship = [V.compare_to_obs(p, drv, quiet=True, mode="equilibrium",
                             lipid_source="measured") for p in (kod, LI2019)]
    anch = [V.compare_to_obs(p, drv, quiet=True, mode="equilibrium",
                             lipid_source="briggs_anchor") for p in (kod, LI2019)]
    assert ship[0] < anch[0] and anch[1] < ship[1]     # same opposite verdicts
    assert anch[0] - ship[0] > 0.10                    # but a much wider margin


# --------------------------------------------------------------------------
# lipid_source as a named mode
# --------------------------------------------------------------------------
def test_lipid_source_default_is_measured_and_changes_nothing():
    """The mode is additive: selecting the default explicitly must be the same
    object graph as not selecting anything. If this ever fails, every published
    neutral number silently changed meaning."""
    assert ND.DEFAULT_LIPID_SOURCE == "measured"
    assert ND.LIPID_SOURCES["measured"] is ND.TRAPP1994_LIPID_FW
    assert ND.LIPID_SOURCES["briggs_anchor"] is ND.BRIGGS_ANCHORED_LIPID_FW

    implicit = ND.rice_compartments()
    explicit = ND.rice_compartments(lipid_source="measured")
    for a, b in zip(implicit, explicit):
        assert (a.name, a.theta, a.f_PL) == (b.name, b.theta, b.f_PL)


def test_lipid_source_selects_the_two_readings_and_lipids_still_overrides():
    anchored = {c.name: c for c in ND.rice_compartments(lipid_source="briggs_anchor")}
    measured = {c.name: c for c in ND.rice_compartments(lipid_source="measured")}
    # 2.48x in the lipid term, root only
    assert anchored["root"].f_PL / measured["root"].f_PL == pytest.approx(2.48, abs=0.02)
    for organ in ("stem", "leaf", "grain"):
        assert anchored[organ].f_PL == pytest.approx(measured[organ].f_PL)
    # an explicit `lipids` dict still wins over the selected mode
    forced = {c.name: c for c in ND.rice_compartments(lipid_source="briggs_anchor",
                                                      lipids={"root": 0.01})}
    assert forced["root"].f_PL == pytest.approx(measured["root"].f_PL)
    with pytest.raises(ValueError):
        ND.rice_compartments(lipid_source="nonsense")


def test_lipid_source_reaches_the_forward_run_and_the_public_api():
    """The selector has to survive the whole way to `model_api`, or the mode is
    documentation rather than a switch. Checked on a lipophilic compound, where
    the two readings are ~2.4x apart, and on a hydrophilic one, where they are
    not -- the Kow dependence is the point (see LIPID_SOURCES)."""
    import model_api as api

    hi = [api.simulate_neutral(4.0, season=60.0, n_t=61, lipid_source=s)
          for s in ("measured", "briggs_anchor")]
    assert hi[1]["baf_final"]["root"] / hi[0]["baf_final"]["root"] > 2.0
    assert hi[0]["params"]["lipid_source"] == "measured"
    assert hi[1]["params"]["lipid_source"] == "briggs_anchor"

    lo = [api.simulate_neutral(0.0, season=60.0, n_t=61, lipid_source=s)
          for s in ("measured", "briggs_anchor")]
    assert lo[1]["baf_final"]["root"] / lo[0]["baf_final"]["root"] < 1.1

    # and omitting it is the default, bit-identical
    dflt = api.simulate_neutral(4.0, season=60.0, n_t=61)
    assert dflt["baf_final"]["root"] == hi[0]["baf_final"]["root"]


def test_compare_lipid_sources_reproduces_the_three_to_one_tally():
    """The anchor decision as a command. The docs claim 3 tables against and 1
    for; this asserts the code still says so, so the doc claim cannot go stale
    while the numbers underneath move."""
    import neutral_dpu_validation as V

    rows = V.compare_lipid_sources(drv=V.drivers())
    by = {label: (m, a) for label, m, a in rows}
    assert len(rows) == 5
    # the anchor wins only on the two the docs say it wins on
    won = {label for label, (m, a) in by.items() if a < m}
    assert won == {"Ge 2017 (per-organ)", "Li 2019 hydroponic"}
    # and the table it damages is the SOIL half of the paper it most helps
    m_soil, a_soil = by["Li 2019 soil"]
    assert a_soil - m_soil > 0.10
    # NOTE the basis: this is mode="ode". docs section 4h quotes 0.639 -> 0.873
    # for the same table, which is the EQUILIBRIUM basis that li2019_soil_table.py
    # scores on (k_pw directly, no ODE). Both are right; neither number means
    # anything without its mode, which is why they are labelled in the docs.
    assert m_soil == pytest.approx(0.549, abs=0.02)


def test_spearman_helper_handles_ties():
    import schriever2020_tscf as S

    assert S.spearman([1, 2, 3, 4], [1, 2, 3, 4]) == pytest.approx(1.0)
    assert S.spearman([1, 2, 3, 4], [4, 3, 2, 1]) == pytest.approx(-1.0)
    # ties must be averaged, not left in arbitrary encounter order
    assert S.spearman([1, 1, 2, 2], [1, 1, 2, 2]) == pytest.approx(1.0)
    assert np.isnan(S.spearman([1, 1, 1], [1, 2, 3]))


# --------------------------------------------------------------------------
# Li 2019 Table S2 -- the soil companion, which disagrees in SIGN
# --------------------------------------------------------------------------
LI2019_SOIL = os.path.join(ROOT, "data_obs", "neutral_obs_li2019_soil.csv")


def test_li2019_soil_table_is_intact():
    rows = _rows(LI2019_SOIL)
    assert len(rows) == 376                       # all of the SI's rows
    assert len({r["species"] for r in rows}) == 13
    assert Counter(r["kom_source"] for r in rows) == {
        "estimated": 261, "experimental": 62, "unmatched": 53}
    # the 20 Zhang 2005 rows are a field survey with no stated duration
    assert sum(1 for r in rows if r["exposure_d"] == "") == 20
    # per-crop lipid must span the ~12x range the SI carries. Note the SI gives
    # radish as both 0.09 and 0.10 % in different source studies while the
    # article text quotes 0.10; the SI value is kept per row.
    fl = [float(r["f_lip"]) for r in rows]
    assert min(fl) == pytest.approx(0.0009) and max(fl) == pytest.approx(0.0114)


def test_the_two_halves_of_li2019_disagree_in_sign():
    """The finding this table exists for. Same paper, same authors, same
    fresh-weight convention -- but the hydroponic half says the model is LOW and
    the soil half, twelve times larger, says it is HIGH. That is what makes 'the
    root partition is too low' a property of one exposure route rather than of
    the partition core, and it is why the anchor decision stays open."""
    import li2019_soil_table as S

    soil = S.load()
    n_hy, b_hy = S._hydroponic_bias()
    n_so, b_so, _ = S.bias(soil)
    assert n_hy == 29 and n_so == 376
    assert b_hy < -0.3 and b_so > +0.2            # opposite signs, both large


def test_soil_table_argues_against_the_anchor_and_own_lipid_removes_the_bias():
    import li2019_soil_table as S

    rows = S.load()
    _, b_rice, e_rice = S.bias(rows, "rice")
    _, b_anch, e_anch = S.bias(rows, "anchor")
    _, b_own, e_own = S.bias(rows, "own")
    assert b_anch > b_rice and e_anch > e_rice    # the anchor makes it worse
    assert abs(b_own) < 0.05                      # own lipid -> unbiased
    assert e_own < e_rice                         # and tighter


def test_measured_kom_rows_are_nearly_unbiased():
    """Where the exposure conversion is measured rather than estimated, the model
    is close to unbiased -- which relocates most of this table's apparent error
    from the plant model to the soil side."""
    import li2019_soil_table as S

    rows = S.load()
    n_exp, b_exp, _ = S.bias([r for r in rows if r["kom_source"] == "experimental"])
    n_est, b_est, _ = S.bias([r for r in rows if r["kom_source"] == "estimated"])
    assert n_exp == 62 and n_est == 261
    assert abs(b_exp) < 0.10 < b_est


# --------------------------------------------------------------------------
# Kodesova's leaf half: translocation, and metabolism that was MEASURED
# --------------------------------------------------------------------------
def test_metabolism_is_measured_not_fitted_and_is_species_dependent():
    """Section 7. Every in-planta half-life in this repo has been an inference
    from a fit; Kodesova measured the metabolites, so the parent fraction is a
    direct observation. Two things must hold, and both bear on how the Ge 2017
    half-life scan should be read."""
    import kodesova2019_carbamazepine as K

    root = {k: v[0] / sum(v) for k, v in K.KOD_ROOT.items()}
    leaf = {k: v[0] / sum(v) for k, v in K.KOD_LEAF.items()}
    assert np.mean(list(root.values())) > 0.85     # root stays mostly parent
    assert np.mean(list(leaf.values())) < 0.60     # the shoot transforms it
    # and it is a SPECIES property, not a compound one -- same compound, same
    # soils, same harvest, 4x apart
    lettuce = np.mean([v for k, v in leaf.items() if k[0] == "lamb's lettuce"])
    radish = np.mean([v for k, v in leaf.items() if k[0] == "radish"])
    assert radish / lettuce > 3.0


def test_metabolism_alone_cannot_close_the_leaf_root_gap():
    """Section 6, and the warning it carries for Ge 2017: even a 1.5-day
    half-life leaves the model far above the measured leaf/root, so a half-life
    fitted to leaf data is absorbing the driver mismatch, not just metabolism."""
    import kodesova2019_carbamazepine as K
    import neutral_dpu_validation as V

    drv = V.drivers()
    t = np.asarray(drv["t"])
    i = int(np.argmin(np.abs(t - K.HARVEST_D)))
    obs = np.median([l / r for r, l in K.KOD_LEAF_ROOT])
    assert 2.0 < obs < 5.0

    def model_leaf_root(half_life):
        g = float(np.log(2.0) / half_life)
        comps = ND.rice_compartments(gammas={k: g for k in
                                             ("root", "stem", "leaf", "grain")})
        m = ND.simulate_neutral(ND.NeutralCompound("carbamazepine", K.LOG_KOW),
                                drv, comps=comps)
        return m["conc"]["leaf"][i] / m["conc"]["root"][i]

    assert model_leaf_root(1.5) > 5.0 * obs        # still far above at 1.5 d
    assert model_leaf_root(1.5) < model_leaf_root(7.0)   # shorter does help, but


# --------------------------------------------------------------------------
# Briggs 1983 -- the stem anchor
# --------------------------------------------------------------------------
def test_briggs1983_scf_reproduces_the_papers_own_stated_maximum():
    """The transcription check. Briggs et al. 1983 state their eq. 3 peaks at
    'about 6 ... at about log Kow = 4.5'; computing that from the coefficients
    alone is what makes the implementation a check rather than an assertion."""
    lk = np.linspace(-1.0, 7.0, 8001)
    scf = np.array([ND.briggs_scf(x) for x in lk])
    i = int(np.argmax(scf))
    assert scf[i] == pytest.approx(6.0, abs=0.6)
    assert lk[i] == pytest.approx(4.5, abs=0.15)


def test_stem_runs_above_its_anchor_but_the_observable_largely_cancels():
    """The stem's mismatch, recorded and NOT fixed -- and recorded at the right
    size. The coefficient gap is 4.1x, but against Briggs' steeper exponent the
    observable (SCF) differs by at most ~0.13 log over the range where TSCF
    delivers anything. A later session must not quote the 4.1x as if it were the
    error in a prediction."""
    shipped = ND.TRAPP1994_LIPID_FW["stem"] * ND.LIPID_OCTANOL_A
    anchor = 10.0 ** ND.STEM_INTERCEPT
    assert shipped / anchor == pytest.approx(4.1, abs=0.1)
    assert ND.STEM_SLOPE > ND.RCF_SLOPE          # steeper than the root's

    W = ND.RICE_WATER["stem"]
    L = ND.TRAPP1994_LIPID_FW["stem"]
    diffs = [abs(np.log10(ND.k_pw(lk, W=W, L=L) * ND.tscf(lk) / ND.briggs_scf(lk)))
             for lk in (0.0, 1.0, 1.78, 2.5, 3.5)]
    assert max(diffs) < 0.15                     # the observable nearly cancels

    # the two partitions cross, so this is a shape difference not an offset
    below = ND.k_pw(2.0, W=W, L=L) > ND.briggs_stem_xylem_partition(2.0)
    above = ND.k_pw(5.0, W=W, L=L) < ND.briggs_stem_xylem_partition(5.0)
    assert below and above


def test_root_and_stem_anchors_are_missed_in_opposite_directions():
    """Read together, the two anchors say the neutral composition was assembled
    from two unrelated sources and neither organ was checked against the anchor
    that existed for it. No single correction fixes both."""
    root_ratio = (ND.TRAPP1994_LIPID_FW["root"] * ND.LIPID_OCTANOL_A) / 10 ** ND.RCF_INTERCEPT
    stem_ratio = (ND.TRAPP1994_LIPID_FW["stem"] * ND.LIPID_OCTANOL_A) / 10 ** ND.STEM_INTERCEPT
    assert root_ratio < 1.0 < stem_ratio
