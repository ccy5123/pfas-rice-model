"""
Guards for the policy-council briefing -- validation/policy_brief_runs.py and
docs/POLICY_BRIEF_KR.md.

The briefing is written to be handed to someone (or something) that will build
slides from it without re-deriving anything. That makes two things load-bearing
and worth pinning:

  * the NUMBERS, because a slide quoting a stale figure is worse than no slide;
  * the REFUSALS, because this model's headline result for the grain is
    mechanism-dependent and the briefing's whole job is to stop that being
    presented as settled.

The heavy ODE runs are not repeated here -- the briefing's own CSV is the
artifact, and these tests check it exists, agrees with the document, and that
the document still carries the caveats it was written around.
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
    if not os.path.exists(RESULTS):
        pytest.skip("run validation/policy_brief_runs.py first")
    with open(RESULTS, encoding="utf-8") as fh:
        return {r["key"]: r for r in csv.DictReader(fh)}


# --------------------------------------------------------------------------
# the numbers
# --------------------------------------------------------------------------
def test_results_csv_covers_every_section(results):
    secs = {r["section"] for r in results.values()}
    assert {"1_congeners", "2_pfoa", "3_kim2019", "4_inverse", "5_smiles"} <= secs


def test_worked_example_matches_the_document(brief, results):
    """Slide 4's table is the one a policy audience will read off the screen."""
    grain = float(results["conc.grain"]["value"])
    assert grain == pytest.approx(0.0136, abs=5e-4)
    assert "0.0136" in brief
    twi = float(results["efsa_twi_percent"]["value"])
    assert twi == pytest.approx(5.4, abs=0.1)
    assert "5.4 %" in brief
    band = float(results["grain.band_factor"]["value"])
    assert band == pytest.approx(7.08, abs=0.05)


def test_kim2019_scores_are_the_documented_out_of_sample_result(results):
    """The briefing's central evidence claim: the base model misses the Korean
    field grain badly, and the lipid pathway -- fitted on a DIFFERENT dataset --
    brings it inside a factor of ~6. These are the numbers CLAUDE.md records for
    the same comparison, so a drift here is a drift in the model, not the doc."""
    mono = float(results["rmse_log10_monotone"]["value"])
    lip = float(results["rmse_log10_lipid"]["value"])
    assert mono > 1.5
    assert lip < 0.7
    assert mono - lip > 1.0
    assert float(results["bias_log10_monotone"]["value"]) < -1.0   # under-predicts


def test_inverse_recovers_a_known_truth_but_misses_the_field(results):
    """Both halves matter: the method works, and it is not yet field-ready. A
    later session must not quote the first without the second."""
    truth = float(results["synthetic.truth"]["value"])
    got = float(results["synthetic.median"]["value"])
    assert got == pytest.approx(truth, rel=0.05)
    assert float(results["synthetic.ci95_lo"]["value"]) < truth < \
        float(results["synthetic.ci95_hi"]["value"])
    assert float(results["field.overestimate_factor"]["value"]) > 10


def test_structure_input_flags_the_unknown_ones(results):
    assert float(results["PFOA.provisional"]["value"]) == 0.0
    assert float(results["PFPeS.provisional"]["value"]) == 1.0


# --------------------------------------------------------------------------
# the refusals -- the part a slide agent must not quietly drop
# --------------------------------------------------------------------------
def test_the_grain_ordering_is_never_presented_as_settled(brief):
    """The base model says short chains reach the grain; both measured datasets
    say the long ones do. The briefing must keep saying so."""
    assert "451배 과소" in brief
    for phrase in ("단쇄 PFAS가 쌀에 더 위험하다", "확정된 결과가 아니라"):
        assert phrase in brief


def test_every_slide_carries_its_caveat_block(brief):
    """Slides 3-9 each end on a '반드시 남길 것' or an explicit limits list."""
    n_slides = len(re.findall(r"^## 슬라이드 \d+", brief, re.M))
    assert n_slides >= 10
    assert brief.count("⚠ 반드시 남길 것") >= 6


def test_the_forbidden_phrase_table_exists(brief):
    assert "부록 A — 금지 표현과 대체 표현" in brief
    for banned in ("모델이 검증되었다", "기준치를 초과한다", "AI가 예측한다",
                   "규제 판단에 쓸 수 있다"):
        assert banned in brief


def test_no_legal_limit_is_claimed(brief):
    """There is no PFAS maximum level for rice in either jurisdiction; the
    briefing has to keep saying that, since a slide implying otherwise would be
    the single most damaging error this document could seed."""
    assert "법정 기준치는 EU에도 한국에도 없습니다" in brief
    assert "법정 기준 초과 판정은 불가" in brief


def test_reproduction_instructions_point_at_real_files(brief):
    for path in ("validation/policy_brief_runs.py",
                 "docs/policy_brief_results.csv",
                 "validation/figures/policy_grain_by_chain.png",
                 "validation/figures/policy_kim2019_grain.png"):
        assert path in brief
        assert os.path.exists(os.path.join(_ROOT, path)) or path.endswith(".csv")


def test_module_imports_and_helpers_are_pure():
    import policy_brief_runs as PB
    assert PB.KIM_PFOA_POREWATER == pytest.approx(0.0787)
    assert callable(PB.section1) and callable(PB.figures)
