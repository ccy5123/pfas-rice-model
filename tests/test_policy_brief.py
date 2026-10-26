"""
Guards for the policy-council briefing -- validation/policy_brief_runs.py and
docs/POLICY_BRIEF_KR.md.

The briefing is written to be handed to someone (or something) that will build
slides from it without re-deriving anything. That makes two things load-bearing
and worth pinning:

  * the NUMBERS, because a slide quoting a stale figure is worse than no slide
    (the first version of this briefing quoted pre-correction numbers for months
    after the model under it had changed -- that is why these pins exist);
  * the REFUSALS, because the grain is the model's weakest compartment and the
    briefing's job is to stop its numbers being presented as settled.

The heavy ODE runs are not repeated here -- the briefing's own CSV is the
artifact, and these tests check it exists, was made on the current default, agrees
with the document, and that the document still carries the caveats it was written
around.
"""
import csv
import os
import re
import sys

import pytest

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "src"))
sys.path.insert(0, os.path.join(_ROOT, "validation"))

BRIEF = os.path.join(_ROOT, "docs", "POLICY_BRIEF_KR.md")
RESULTS = os.path.join(_ROOT, "docs", "policy_brief_results.csv")


@pytest.fixture(scope="module")
def brief():
    return open(BRIEF, encoding="utf-8").read()


@pytest.fixture(scope="module")
def results():
    """Keyed by (section, key): the worked example is run at two pore-water levels
    and the two share key names, so a key-only dict would silently keep one."""
    if not os.path.exists(RESULTS):
        pytest.skip("run validation/policy_brief_runs.py first")
    with open(RESULTS, encoding="utf-8") as fh:
        return {(r["section"], r["key"]): r for r in csv.DictReader(fh)}


def v(results, section, key):
    return float(results[(section, key)]["value"])


# --------------------------------------------------------------------------
# the numbers
# --------------------------------------------------------------------------
def test_results_csv_covers_every_section(results):
    secs = {s for s, _ in results}
    assert {"0_config", "1_congeners", "2_pfoa", "2_pfoa_50ngL", "3_kim2019",
            "4_inverse", "5_smiles"} <= secs


def test_the_csv_was_made_on_the_current_default(results):
    """A CSV left over from an older default would pass every number test below
    against itself -- so check what it was run on."""
    import model_api as api
    assert results[("0_config", "shoot")]["value"] == api.DEFAULT_SHOOT
    assert bool(v(results, "0_config", "lipid_loading")) == api.DEFAULT_LIPID_LOADING
    assert v(results, "0_config", "band_log10") == pytest.approx(api.APRIORI_LOG10_RMSE, abs=1e-4)


def test_worked_example_matches_the_document(brief, results):
    """Slide 4's table is the one a policy audience will read off the screen."""
    assert v(results, "2_pfoa", "cwo") == pytest.approx(0.0787)
    assert v(results, "2_pfoa", "conc.root") == pytest.approx(0.383, abs=5e-4)
    assert v(results, "2_pfoa", "conc.straw") == pytest.approx(0.247, abs=5e-4)
    assert v(results, "2_pfoa", "conc.grain") == pytest.approx(0.0287, abs=5e-5)
    for s in ("**0.383**", "**0.247**", "**0.0287**"):
        assert s in brief
    assert v(results, "2_pfoa", "efsa_twi_percent") == pytest.approx(8.07, abs=0.01)
    assert "**8 %**" in brief
    assert v(results, "2_pfoa", "grain.band_factor") == pytest.approx(3.08, abs=0.01)
    assert "×/÷3.1" in brief
    assert v(results, "2_pfoa", "grain_via_xylem") == pytest.approx(0.49, abs=0.005)
    assert "물관 **49 %** · 체관 **51 %**" in brief
    assert v(results, "2_pfoa", "root_from_soil") == pytest.approx(0.998, abs=5e-4)
    assert "**99.8 %**" in brief
    # the root entry is linear: the 50 ng/L example is the 78.7 ng/L one scaled
    for k in ("conc.root", "conc.straw", "conc.grain"):
        assert v(results, "2_pfoa_50ngL", k) / v(results, "2_pfoa", k) == pytest.approx(
            0.05 / 0.0787, rel=1e-6)
    assert "농도는 노출에 비례합니다" in brief
    # the deck's 50 ng/L example
    assert v(results, "2_pfoa_50ngL", "conc.grain") == pytest.approx(0.018, abs=5e-4)
    assert v(results, "2_pfoa_50ngL", "efsa_twi_percent") == pytest.approx(5.12, abs=0.01)
    assert "EFSA 대비 5 %" in brief
    # Kim's PFOA brown rice against the model's, both dry
    assert 0.349 / v(results, "2_pfoa", "conc.grain_dw") == pytest.approx(8.26, abs=0.01)
    assert "**약 8배 낮게**" in brief


def test_the_field_grain_table_matches(brief, results):
    """Slide 5: the Yamazaki brown-rice comparison. The grain is fitted per congener
    (L_Ph), so the table is in-sample: six congeners sit on the measurement and four
    stay above it at L_Ph = 0 -- the document must say both, with the CSV's ratios."""
    ratios = {}
    for c in ("PFBA", "PFPeA", "PFHxA", "PFHpA", "PFOA", "PFNA", "PFDA",
              "PFUnDA", "PFDoDA", "PFOS"):
        ratios[c] = v(results, "1_congeners", f"{c}.yz_grain") / \
            v(results, "1_congeners", f"{c}.obs_grain")
    assert not [r for r in ratios.values() if r >= 2.0 or r < 0.99]
    on = sorted(c for c, r in ratios.items() if abs(r - 1.0) < 0.005)
    assert on == sorted(("PFHxA", "PFOA", "PFNA", "PFDA", "PFDoDA", "PFOS"))
    assert "10종 중 6종(PFHxA·PFOA·PFNA·PFDA·PFDoDA·PFOS)" in brief
    for c, s in (("PFUnDA", 1.7), ("PFPeA", 1.5), ("PFHpA", 1.3), ("PFBA", 1.2)):
        assert round(ratios[c], 1) == s and f"{c} {s}배" in brief
    assert "표본 내 비교" in brief
    prev_miss = v(results, "1_congeners", "PFDoDA.obs_grain") / \
        v(results, "1_congeners", "PFDoDA.yz_grain_prev")
    assert 600 < prev_miss < 630 and "약 615배" in brief
    assert v(results, "1_congeners", "mean_grain_shortchain_C4_C6_base") == pytest.approx(1.03, abs=0.005)
    assert v(results, "1_congeners", "mean_grain_longchain_C9_C12_base") == pytest.approx(15.9, abs=0.05)


def test_kim2019_scores_are_the_documented_result(brief, results):
    """Slide 6: no parameter is fitted to Kim, but Kim helped choose the default,
    so the document must say both."""
    assert v(results, "3_kim2019", "rmse_point_default") == pytest.approx(0.49, abs=0.005)
    assert v(results, "3_kim2019", "rmse_interval_default") == pytest.approx(0.37, abs=0.005)
    assert v(results, "3_kim2019", "rmse_point_prev") == pytest.approx(1.93, abs=0.005)
    assert v(results, "3_kim2019", "rmse_interval_prev") == pytest.approx(1.67, abs=0.005)
    for s in ("**0.49**", "**0.37**", "1.93", "1.67"):
        assert s in brief
    # Kim weighed freeze-dried grain: the comparison is on the dry basis
    assert v(results, "3_kim2019", "PFOA.model_grain_baf") == pytest.approx(
        v(results, "1_congeners", "PFOA.grain"))
    assert "건물중으로 읽습니다" in brief
    assert "구성 선택 면에서는 이 자료가 이미 쓰였습니다" in brief
    # the table's measured columns are the literature values, carried in the CSV
    assert v(results, "3_kim2019", "PFOA.porewater") == pytest.approx(78.7)
    assert v(results, "3_kim2019", "PFOA.brownrice") == pytest.approx(0.349)
    assert v(results, "3_kim2019", "PFOA.detection_freq") == pytest.approx(57.0)
    # PFOA is reported against the per-site constant, and the model sits below it
    assert v(results, "3_kim2019", "PFOA.model_grain_baf") < v(results, "3_kim2019",
                                                              "pfoa_site_geomean_lo90")
    assert "9.94" in brief


def test_inverse_recovers_a_known_truth_but_misses_the_field(brief, results):
    """Both halves matter: the method works, and it is not yet field-ready. A
    later session must not quote the first without the second. With a linear root
    entry the synthetic recovery is a check of the procedure only, and the field
    interval contains the truth because it is wide -- the slide says both."""
    truth = v(results, "4_inverse", "synthetic.truth")
    assert v(results, "4_inverse", "synthetic.median") == pytest.approx(truth, rel=0.05)
    assert v(results, "4_inverse", "synthetic.ci95_lo") < truth < \
        v(results, "4_inverse", "synthetic.ci95_hi")
    assert "계산 절차의 점검일 뿐 정확도의 증거가 아닙니다" in brief
    f = v(results, "4_inverse", "field.overestimate_factor")
    assert f == pytest.approx(8.26, abs=0.01)
    assert v(results, "4_inverse", "field.ci95_lo") < truth < v(results, "4_inverse", "field.ci95_hi")
    assert "약 8배 과대추정" in brief and "실측을 겨우 포함함" in brief
    assert v(results, "4_inverse", "field.overestimate_factor_prev") == pytest.approx(13.4, abs=0.05)


def test_structure_input_flags_the_unknown_ones(results):
    assert v(results, "5_smiles", "PFOA.provisional") == 0.0
    assert v(results, "5_smiles", "PFPeS.provisional") == 1.0
    # a SMILES-built PFOA reproduces the named congener exactly
    assert v(results, "5_smiles", "PFOA.grain") == pytest.approx(v(results, "1_congeners", "PFOA.grain"))


def test_the_biocide_slide_quotes_the_biocide_record(brief, results):
    """Slide 8 quotes numbers made by another script; section 6 reads them back."""
    if ("6_biocides", "rho_fishBCF_root") not in results:
        pytest.skip("docs/policy_biocides_results.csv was absent when the brief ran")
    assert round(v(results, "6_biocides", "rho_fishBCF_root"), 2) == 0.89
    assert round(v(results, "6_biocides", "rho_fishBCF_straw"), 2) == -0.78
    assert v(results, "6_biocides", "max_tscf_bioaccumulative") < 0.021
    assert "+0.89" in brief and "−0.78" in brief and "0.021" in brief
    assert "**1.03배**" in brief
    assert v(results, "6_biocides", "n_runnable") == 117


# --------------------------------------------------------------------------
# the refusals -- the part a slide agent must not quietly drop
# --------------------------------------------------------------------------
def test_superseded_numbers_are_named_and_not_quoted(brief):
    """The pre-correction figures are named once in the revision note, so a reader
    holding the old deck knows they are wrong -- and never appear as results."""
    head, body = brief.split("---", 1)
    for old in ("±7배", "0.0136", "451배",                        # 2026-10-25 revision
                "0.376", "0.113", "45 %", "80 %", "2.3–7.3배", "2,300배",   # amendment 1
                "2.7배", "1.05배", "Km"):
        assert old in head
        assert old not in body
    assert "2026-10-26 개정(사전등록 개정 1)" in head


def test_the_model_description_is_the_shipped_model(brief):
    """Amendment 1: a linear equilibrium root entry (no GHK exclusion, no carrier as
    the mechanism), retention 1, one loading constant for every chain."""
    head, body = brief.split("---", 1)
    assert "공극수 농도에 비례하는 평형형" in body
    assert "보유율 1" in body
    assert "물관 적재 상수는 모든 물질에 같은 값 하나씩" in body
    for gone in ("약 100배 배제", "서로 밀어냅니다", "단순 비례하지 않습니다", "포화형(비선형)"):
        assert gone not in body


def test_the_grain_ordering_is_never_presented_as_settled(brief):
    for phrase in ("단쇄 PFAS가 쌀에 더 위험하다", "확정된 결과가 아니라",
                   "실험으로 확정되지 않았습니다"):
        assert phrase in brief


def test_every_slide_carries_its_caveat_block(brief):
    """Slides 3-8 each end on a '반드시 남길 것'; slide 9 is the limits list."""
    n_slides = len(re.findall(r"^## 슬라이드 \d+", brief, re.M))
    assert n_slides >= 10
    assert brief.count("⚠ 반드시 남길 것") >= 6
    assert "## 슬라이드 9 — 지금 이 모델로 하면 안 되는 것" in brief


def test_the_selection_caveat_is_kept(brief):
    """Tang, Kim and Li were not fitted but did choose the default; the ×3 band is
    optimistic for that reason and the document has to say so."""
    assert "기본 구성 선택에 쓰였으므로 실제 오차는 더 클 수 있습니다" in brief
    assert "낙관적일 수 있음" in brief


def test_the_forbidden_phrase_table_exists(brief):
    assert "부록 A — 금지 표현과 대체 표현" in brief
    for banned in ("모델이 검증되었다", "기준치를 초과한다", "AI가 예측한다",
                   "규제 판단에 쓸 수 있다"):
        assert banned in brief


def test_no_legal_limit_is_claimed(brief):
    """There is no PFAS maximum level for rice in either jurisdiction; a slide
    implying otherwise would be the single most damaging error this document
    could seed."""
    assert "법정 기준치는 EU에도 한국에도 없습니다" in brief
    assert "법정 기준 초과 판정은 불가" in brief


def test_reproduction_instructions_point_at_real_files(brief):
    for path in ("validation/policy_brief_runs.py",
                 "validation/bat_census_biocides.py",
                 "validation/policy_biocides_figures.py",
                 "docs/policy_brief_results.csv",
                 "validation/figures/policy_grain_by_chain.png",
                 "validation/figures/policy_kim2019_grain.png"):
        assert path in brief
        assert os.path.exists(os.path.join(_ROOT, path)), path


def test_module_imports_and_helpers_are_pure():
    import model_api as api
    import policy_brief_runs as PB
    assert PB.KIM_PFOA_POREWATER == pytest.approx(0.0787)
    assert PB.PREV == dict(api.FOUR_POOL_FREE)
    assert callable(PB.section1) and callable(PB.figures)
