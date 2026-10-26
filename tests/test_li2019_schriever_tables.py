"""Guards for the two measured tables that arrived with the A1/A3 papers.

These lock down three things that a later session could otherwise undo quietly:

  * the transcriptions themselves (row counts, species tallies, ranges), so a
    re-parse that silently drops rows is caught;
  * the `subset` mechanism, whose entire job is to keep a table's CALIBRATION
    rows out of an a-priori score -- and, just as important, to leave the two
    tables WITHOUT that column (liu2023, ge2017) untouched, so the published
    0.281 / 0.783 cannot move because of it;
  * the neutral composition of PREREG 5.10: lipid on a dry basis restated at the
    organ water, each table on its own plant, and the Briggs anchor as the named
    alternative for the root.
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
def test_default_lipid_is_dry_basis_restated_at_the_organ_water():
    """PREREG 5.10: the rice lipid is carried on a dry basis (Li 2019 cereal roots for the root,
    rice straw ether extract for stem and leaf, brown rice crude fat for the grain) and restated
    at the organ's water, so an override of the water moves the fresh lipid with it."""
    import neutral_composition as NC
    assert ND.RICE_LIPID_DW == pytest.approx({"root": 0.0835, "stem": 0.015, "leaf": 0.015,
                                              "grain": 0.026})
    comps = {c.name: c for c in ND.rice_compartments()}
    for organ, c in comps.items():
        assert c.f_PL == pytest.approx(ND.RICE_LIPID_DW[organ])      # f_PL = L_fw/(1-W) = L_dw
        assert c.theta == ND.RICE_WATER[organ]
    wet = {c.name: c for c in ND.rice_compartments(waters={"root": 0.90})}
    assert wet["root"].f_PL * (1 - 0.90) == pytest.approx(0.0835 * 0.10)
    # the rice root's fresh lipid at its sourced water sits near Briggs' fitted term
    assert ND.RICE_LIPID_FW["root"] == pytest.approx(0.0835 * (1 - ND.RICE_WATER["root"]))
    assert ND.RICE_LIPID_FW["root"] / ND.BRIGGS_ANCHORED_LIPID_FW["root"] == pytest.approx(
        0.97, abs=0.02)
    assert NC.MID == ("mid", "mid") and len(NC.CORNERS) == 4


def test_briggs_anchored_root_reproduces_briggs_rcf():
    """The named alternative: with the anchor the root IS the Briggs RCF curve above the
    water floor. Compared at a fixed log Kow, so nothing else can drift in."""
    log_kow = 4.0
    W = ND.RICE_WATER["root"]
    anchored = ND.k_pw(log_kow, W=W, L=ND.BRIGGS_ANCHORED_LIPID_FW["root"])
    assert (anchored - W) == pytest.approx(ND.briggs_rcf(log_kow) - ND.RCF_FLOOR, rel=1e-9)
    for organ in ("stem", "leaf", "grain"):
        assert ND.BRIGGS_ANCHORED_LIPID_FW[organ] == ND.RICE_LIPID_FW[organ]


def test_each_table_runs_on_its_own_plant():
    """PREREG 5.10: a table's composition comes from its own plant (sources in
    params/neutral_composition_sourced.csv); organs it gives no value for keep the rice
    compartments'."""
    import neutral_composition as NC
    import organ_mass as om
    W, L = NC.spec("li2019_rcf", species="Barley")
    assert W["root"] == 0.90 and L["root"] == pytest.approx(0.0100)    # Li's own basis
    W, L = NC.spec("li2019_rcf", species="Parrot feather")             # not in Li's list
    assert L["root"] == pytest.approx(0.5 * (0.010 + 0.114) * 0.10)    # the 13-crop range
    W, L = NC.spec("kodesova2019", species="radish")
    assert W["root"] == pytest.approx(0.9365)
    assert L["root"] == pytest.approx(0.010 * (1 - 0.9365))
    W, L = NC.spec("ge2017")
    assert W["root"] == pytest.approx(0.848) and W["leaf"] == 0.731
    assert W["stem"] == om.THETA_FW["stem"]                            # no vegetative source
    W, L = NC.spec("li2019_soil", f_lip=0.0053)
    assert L["root"] == pytest.approx(0.0053)                          # the row's own f_lip
    W, L = NC.spec("liu2023", corner=("low", "high"))
    assert W["root"] == 0.824 and L["root"] == pytest.approx(0.114 * (1 - 0.824))


def test_ge2017_stem_rows_are_not_scored():
    """No rice stem water at the vegetative stage was found, so Ge 2017's stem TFs stay out of
    the score (PREREG 5.10) and only its leaf rows count."""
    import neutral_dpu_validation as V
    assert ("ge2017", "stem") in V.UNSCORED
    assert V.table_of(os.path.join(ROOT, "data_obs", "neutral_obs_ge2017.csv")) == "ge2017"


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


def test_sourced_root_beats_the_anchor_on_the_root_tables_but_not_li_hydroponic():
    """On the sourced composition (PREREG 5.10) Kodesova, Liu and Li's soil half are closer to
    the data than with Briggs' fitted root term; Li's hydroponic half is not."""
    import neutral_dpu_validation as V

    drv = V.drivers()
    soil = os.path.join(ROOT, "data_obs", "neutral_obs_li2019_soil.csv")
    liu = os.path.join(ROOT, "data_obs", "neutral_obs_liu2023.csv")
    for p in (KODESOVA, liu, soil):
        assert (V.compare_to_obs(p, drv, quiet=True)
                < V.compare_to_obs(p, drv, quiet=True, lipid_source="briggs_anchor"))
    assert (V.compare_to_obs(LI2019, drv, quiet=True, lipid_source="briggs_anchor")
            < V.compare_to_obs(LI2019, drv, quiet=True))


# --------------------------------------------------------------------------
# the ODE-vs-equilibrium scoring artifact
# --------------------------------------------------------------------------
def test_the_rice_ode_discounts_the_root_in_a_kow_dependent_way():
    """Why `mode="equilibrium"` exists. Running a 120-day rice season to score a
    24-hour barley root measurement imposes a Kow-DEPENDENT discount that has
    nothing to do with the partition being tested: near zero at the extremes,
    but ~0.25 log where the xylem drains the root hardest, around the TSCF peak.
    Kodesova sits at log Kow 2.25, close to the worst point -- so its ODE-basis
    score is flattered. Pinned so the artifact cannot be forgotten again.

    MOSTLY A MASS-BASIS ARTIFACT (src/organ_mass.py): the large discount appeared with
    crop dry matter used as the root mass, a root 10x too light to keep up with its own
    xylem export. On fresh organ mass the ODE root sits within ~12% of K_PW everywhere."""
    import neutral_dpu_validation as V

    drv = V.drivers()
    disc = {}
    for lk in (-0.5, 1.78, 2.25, 5.0):
        m = ND.simulate_neutral(ND.NeutralCompound("x", lk), drv)
        disc[lk] = m["baf_final"]["root"] / V.equilibrium_rcf(lk)
    assert disc[-0.5] > 0.85 and disc[5.0] > 0.95      # ~no discount at the ends
    # and at most ~20 % anywhere (the root is lighter on fresh mass at the sourced water 0.713
    # than it was at 0.90, so it keeps up less with its xylem export)
    assert min(disc.values()) > 0.78


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
    assert kod_eq > kod_ode                            # the ODE flatters Kodesova (log Kow 2.25)
    assert kod_eq < liu_eq                             # on its own plants Kodesova is the best

    # Ge 2017 carries only stem/leaf tf rows, so equilibrium mode cannot touch it
    assert both(ge)[0] == both(ge)[1]

    # the default must stay "ode", so no published number moves silently
    assert V.compare_to_obs(liu, drv, quiet=True) == liu_ode


def test_equilibrium_mode_keeps_the_kodesova_verdict():
    import neutral_dpu_validation as V

    drv = V.drivers()
    ship = V.compare_to_obs(KODESOVA, drv, quiet=True, mode="equilibrium")
    anch = V.compare_to_obs(KODESOVA, drv, quiet=True, mode="equilibrium",
                            lipid_source="briggs_anchor")
    assert ship < anch - 0.10


# --------------------------------------------------------------------------
# lipid_source as a named mode
# --------------------------------------------------------------------------
def test_lipid_source_default_is_measured_and_changes_nothing():
    """The mode is additive: selecting the default explicitly must be the same composition
    as not selecting anything."""
    assert ND.DEFAULT_LIPID_SOURCE == "measured"
    assert ND.LIPID_SOURCES["measured"] is ND.RICE_LIPID_FW
    assert ND.LIPID_SOURCES["briggs_anchor"] is ND.BRIGGS_ANCHORED_LIPID_FW

    implicit = ND.rice_compartments()
    explicit = ND.rice_compartments(lipid_source="measured")
    for a, b in zip(implicit, explicit):
        assert (a.name, a.theta, a.f_PL) == (b.name, b.theta, b.f_PL)


def test_lipid_source_selects_the_two_readings_and_lipids_still_overrides():
    anchored = {c.name: c for c in ND.rice_compartments(lipid_source="briggs_anchor")}
    measured = {c.name: c for c in ND.rice_compartments(lipid_source="measured")}
    assert anchored["root"].f_PL / measured["root"].f_PL == pytest.approx(
        ND.BRIGGS_ANCHORED_LIPID_FW["root"] / ND.RICE_LIPID_FW["root"])
    for organ in ("stem", "leaf", "grain"):
        assert anchored[organ].f_PL == pytest.approx(measured[organ].f_PL)
    # an explicit `lipids` dict still wins over the selected mode
    forced = {c.name: c for c in ND.rice_compartments(
        lipid_source="briggs_anchor", lipids={"root": ND.RICE_LIPID_FW["root"]})}
    assert forced["root"].f_PL == pytest.approx(measured["root"].f_PL)
    with pytest.raises(ValueError):
        ND.rice_compartments(lipid_source="nonsense")


def test_lipid_source_reaches_the_forward_run_and_the_public_api():
    """The selector has to survive the whole way to `model_api`, or the mode is documentation
    rather than a switch."""
    import model_api as api

    hi = [api.simulate_neutral(4.0, lipid_source=s) for s in ("measured", "briggs_anchor")]
    assert hi[1]["baf_final"]["root"] != hi[0]["baf_final"]["root"]
    assert hi[0]["params"]["lipid_source"] == "measured"
    assert hi[1]["params"]["lipid_source"] == "briggs_anchor"
    # and omitting it is the default, bit-identical
    dflt = api.simulate_neutral(4.0)
    assert dflt["baf_final"]["root"] == hi[0]["baf_final"]["root"]


def test_compare_lipid_sources_matches_the_per_table_scores():
    """The two readings as a command, consistent with compare_to_obs table by table."""
    import neutral_dpu_validation as V

    drv = V.drivers()
    rows = V.compare_lipid_sources(drv=drv)
    assert len(rows) == 5
    by = {label: (m, a) for label, m, a in rows}
    won = {label for label, (m, a) in by.items() if a < m}
    assert won == {"Ge 2017 (per-organ)", "Li 2019 hydroponic"}


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


def test_the_two_halves_of_li2019_disagree():
    """Same paper, same fresh-weight convention, both on Li's water 0.90 and species lipid
    (PREREG 5.10): the model is LOW on the hydroponic half and unbiased on the soil half,
    twelve times larger."""
    import li2019_soil_table as S

    soil = S.load()
    n_hy, b_hy = S._hydroponic_bias()
    n_so, b_so, _ = S.bias(soil)
    assert n_hy == 29 and n_so == 376
    assert b_hy < -0.3 and abs(b_so) < 0.05


def test_soil_table_prefers_each_crops_own_lipid():
    import li2019_soil_table as S

    rows = S.load()
    _, b_rice, e_rice = S.bias(rows, "rice")
    _, b_anch, e_anch = S.bias(rows, "anchor")
    _, b_own, e_own = S.bias(rows, "own")
    assert b_anch > b_rice and e_anch > e_rice    # the anchor makes it worse
    assert abs(b_own) < 0.05                      # own lipid -> unbiased
    assert e_own < e_rice                         # and tighter


def test_measured_and_estimated_kom_rows_differ():
    """The soil-to-water conversion moves the result: the rows with a measured K_om sit below
    the model's partition, those with an estimated one above."""
    import li2019_soil_table as S

    rows = S.load()
    n_exp, b_exp, _ = S.bias([r for r in rows if r["kom_source"] == "experimental"])
    n_est, b_est, _ = S.bias([r for r in rows if r["kom_source"] == "estimated"])
    assert n_exp == 62 and n_est == 261
    assert b_exp < 0.0 < b_est


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


def test_metabolism_and_the_leaf_root_gap():
    """Kodesova's leaf/root on fresh organ mass: without metabolism the model's leaf/root
    is far above the measurement (a rice-driver mismatch remains), and a short in-planta
    half-life -- which Kodesova MEASURED in effect (leaf parent fraction 0.49) -- closes
    most of it. (With crop dry matter used as fresh mass the gap had looked like 55x and
    metabolism-proof; that was the mass basis.)"""
    import kodesova2019_carbamazepine as K
    import neutral_dpu_validation as V

    drv = V.drivers()
    t = np.asarray(drv["t"])
    i = int(np.argmin(np.abs(t - K.HARVEST_D)))
    obs = np.median([l / r for r, l in K.KOD_LEAF_ROOT])
    assert 2.0 < obs < 5.0

    def model_leaf_root(half_life):
        g = 0.0 if half_life is None else float(np.log(2.0) / half_life)
        comps = ND.rice_compartments(gammas={k: g for k in ("root", "stem", "leaf", "grain")})
        m = ND.simulate_neutral(ND.NeutralCompound("carbamazepine", K.LOG_KOW), drv, comps=comps)
        return m["conc"]["leaf"][i] / m["conc"]["root"][i]

    assert model_leaf_root(1.5) < model_leaf_root(7.0)   # shorter does help
    assert model_leaf_root(None) > 5.0 * obs             # no metabolism: still far above
    # under PREREG 5.8 (root water 0.713: a lighter fresh root, a higher leaf/root) it takes a
    # ~0.5-d half-life to reach the measurement; 1.5 d leaves it ~2.5x above
    assert model_leaf_root(0.5) < 1.5 * obs < model_leaf_root(1.5)


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


def test_barley_shoot_composition_against_the_stem_anchor():
    """The stem anchor (Briggs 1983) against the barley shoot's sourced composition
    (PREREG 5.10): a lower coefficient and a flatter exponent, so the two partitions cross."""
    import neutral_composition as NC
    W, L = NC.spec("briggs1983")
    W, L = W["stem"], L["stem"]
    assert ND.STEM_SLOPE > ND.RCF_SLOPE          # steeper than the root's
    assert L * ND.LIPID_OCTANOL_A < 10.0 ** ND.STEM_INTERCEPT
    below = ND.k_pw(0.5, W=W, L=L) > ND.briggs_stem_xylem_partition(0.5)
    above = ND.k_pw(3.0, W=W, L=L) < ND.briggs_stem_xylem_partition(3.0)
    assert below and above


# --------------------------------------------------------------------------
# Neutral organics through the SOIL side: the Karickhoff Koc route
# --------------------------------------------------------------------------
def test_koc_neutral_is_karickhoff_and_separate_from_the_pfas_qspr():
    """A neutral compound's soil sorption cannot come from the PFAS chain-length
    Koc QSPR (a per-CF2 group contribution with no meaning for a neutral), so the
    neutral path has its own: Karickhoff 1981, log Koc = 0.989*log Kow - 0.346."""
    import literature_params as lp
    assert lp.koc_neutral(0.0) == pytest.approx(10 ** -0.346, rel=1e-9)
    assert lp.koc_neutral(2.45, log10=True) == pytest.approx(-0.346 + 0.989 * 2.45, rel=1e-9)
    # monotone in log Kow, and a decade of log Kow is ~a decade of Koc (slope ~1)
    assert lp.koc_neutral(5.0) / lp.koc_neutral(4.0) == pytest.approx(10 ** 0.989, rel=1e-9)
    # carbamazepine (log Kow 2.45) lands in the literature band for it, ~100-250 L/kg
    assert 100.0 < float(lp.koc_neutral(2.45)) < 250.0
    # and it is NOT the PFAS QSPR: same numeric argument, different meaning/value
    assert not np.isclose(float(lp.koc_neutral(2.45)), float(lp.koc(2.45)))


def test_neutral_soil_side_plumbing_is_koc_driven_and_pfas_defaults_are_untouched():
    """The neutral compound reaches every soil-side entry point through its own Koc:
    the flooded Cwᵒ(t) shape, its k_leach default and the HYDRUS Kd. The PFAS table
    must not answer for a neutral, and the PFAS defaults must not move."""
    import model_api as api
    import literature_params as lp
    Koc = float(lp.koc_neutral(2.45))
    # k_leach: the HYDRUS-calibrated per-congener TABLE is PFAS-only; a neutral takes
    # the log10(Koc) fallback fit, which is a property of the leaching model, not of PFAS
    rows = {r["congener"]: float(r["k_leach"]) for r in csv.DictReader(
        l for l in open(os.path.join(os.path.dirname(__file__), "..", "params", "cwo_kleach.csv"))
        if not l.startswith("#"))}
    assert api.default_k_leach("PFOA") == pytest.approx(rows["PFOA"])  # the PFAS table answers
    kn = api.default_k_leach(Koc=Koc)
    assert 0.0 <= kn <= 0.15 and kn != pytest.approx(0.02)             # not the flat placeholder
    # the flooded shape runs off log_kow with no congener at all, and is mean-normalised
    t = np.linspace(0.0, 120.0, 241)
    flooded = api.cwo_profile_series(t, level=1.0, profile="flooded", log_kow=2.45)
    assert float(np.mean(flooded)) == pytest.approx(1.0, rel=1e-6)
    assert flooded[-1] < flooded[0]                                    # it leaches
    # a measured Koc overrides the QSPR, and a stronger sorber leaches less
    strong = api.cwo_profile_series(t, level=1.0, profile="flooded", Koc=Koc * 100)
    assert strong[-1] / strong[0] > flooded[-1] / flooded[0]
    # and the whole ODE accepts the shape on the neutral path
    r = api.simulate_neutral(2.45, half_life=7.0, cwo_profile="flooded")
    assert np.isfinite(r["baf_final"]["root"]) and r["baf_final"]["root"] > 0


def test_neutral_bayesian_inverse_recovers_a_known_exposure():
    """The inverse must work for a non-PFAS compound too: same Laplace posterior in
    log10(Cwᵒ), with `simulate_neutral` as the forward model instead of `simulate`."""
    import model_api as api
    truth = 3.0
    fwd = api.simulate_neutral(2.45, half_life=7.0, Cwo=truth)
    meas = {"root": fwd["conc_dw"]["root"][-1], "grain": fwd["conc_dw"]["grain"][-1]}
    est = api.estimate_exposure_bayesian(log_kow=2.45, measured_conc=meas,
                                         neutral_kw=dict(half_life=7.0))
    assert est["median"] == pytest.approx(truth, rel=0.10)
    lo, hi = est["ci95"]
    assert lo < truth < hi
    # the PFAS entry point is unchanged (congener positional, no log_kow)
    pf = api.simulate("PFOA", Cwo=2.0)
    epf = api.estimate_exposure_bayesian("PFOA", {"grain": pf["conc_dw"]["grain"][-1]})
    assert epf["median"] == pytest.approx(2.0, rel=0.10)
