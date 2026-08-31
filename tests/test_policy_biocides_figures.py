# -*- coding: utf-8 -*-
"""`docs/POLICY_BRIEF_BIOCIDES_KR.md` 의 수치를 고정하는 가드.

브리핑은 슬라이드를 만드는 사람이 **모델을 다시 돌리지 않고** 인용하는 문서라서,
문서에 적힌 숫자와 코드가 갈라지면 발표 자리에서만 드러난다. 여기서 막는다.

`validation/policy_biocides_figures.py` 는 census 결과 CSV 를 읽는 층이므로
테스트도 그 CSV 를 읽는다(ODE 를 다시 풀지 않아 빠르다). CSV 가 없으면 skip.
"""
import csv
import os
import sys

import pytest

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "validation"))
sys.path.insert(0, os.path.join(_ROOT, "src"))

DOC = os.path.join(_ROOT, "docs", "POLICY_BRIEF_BIOCIDES_KR.md")
NUM = os.path.join(_ROOT, "docs", "policy_biocides_results.csv")
CENSUS = os.path.join(_ROOT, "validation", "bat_census_biocides.csv")

pytestmark = pytest.mark.skipif(
    not os.path.exists(CENSUS),
    reason="validation/bat_census_biocides.csv 없음 (먼저 census 실행)")


@pytest.fixture(scope="module")
def P():
    import policy_biocides_figures as mod
    return mod


@pytest.fixture(scope="module")
def rows(P):
    return P.load_census()


@pytest.fixture(scope="module")
def doc():
    return open(DOC, encoding="utf-8").read()


@pytest.fixture(scope="module")
def nums():
    if not os.path.exists(NUM):
        pytest.skip("docs/policy_biocides_results.csv 없음")
    out = {}
    with open(NUM, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            out[(r["section"], r["key"])] = r["value"]
    return out


# ---------------------------------------------------------------------------
def test_the_screen_counts_in_the_doc_match_the_run(P, rows, doc):
    """슬라이드 3 의 4등급 표 -- 문서 숫자가 실제 판정과 같아야 한다."""
    c = {v: len([r for r in rows if r["verdict"] == v]) for v in P.ORDER}
    assert c[P.V_IN] == 69 and c[P.V_ROOT] == 10
    assert c[P.V_EXT] == 7 and c[P.V_EXCL] == 31
    assert sum(c.values()) == 117
    for n in ("69", "10", "31", "117"):
        assert n in doc


def test_the_two_ad_bands_are_read_from_shipped_tables(P, doc):
    """AD 경계는 여기서 고른 값이 아니라 data_obs/ 에서 읽은 값이다 -- 그리고
    좁은 쪽(지상부 5.46)이 실질 한계라는 것이 슬라이드 3 의 주장이다."""
    sp = P.spans()
    r_lo, r_hi, r_n = sp["partition"]
    t_lo, t_hi, t_n = sp["tscf"]
    assert (round(r_lo, 2), round(r_hi, 2)) == (-0.66, 8.70)
    assert (round(t_lo, 2), round(t_hi, 2)) == (-1.52, 5.46)
    assert t_hi < r_hi, "지상부가 더 좁다는 것이 이 브리핑의 핵심"
    assert r_n == 560 and t_n == 30
    for s in ("8.70", "5.46", "560", "30"):
        assert s in doc


def test_the_rank_correlation_is_60_substances_with_no_duplicate(P, rows, doc):
    """슬라이드 5. n=60 은 «물질» 수이고 중복이 없어야 한다 -- 두 물질이 log Kow
    두 개로 두 행을 갖고, 두 물질은 두 등급에 중복돼 있어서 행으로 세면 틀린다."""
    pairs = P.rank_pairs(rows)
    assert len(pairs) == 60
    names = [p["substance"].strip().lower() for p in pairs]
    assert len(set(names)) == 60, "한 물질이 두 번 가중되면 안 된다"
    x = [p["bat_A"] for p in pairs]
    assert round(P._rho(x, [p["root"] for p in pairs]), 3) == 0.875
    assert round(P._rho(x, [p["straw"] for p in pairs]), 3) == -0.769
    assert "+0.875" in doc and "−0.769" in doc


def test_the_row_vs_substance_count_is_stated_correctly(P, rows, doc):
    """규칙 6. 117행 = 113물질이고, 차이 4행의 내역이 문서에 적혀 있어야 한다.
    (초안은 115라고 적었는데 파일이 113이었다 -- 그래서 이 테스트가 있다.)"""
    base = {r["substance"].split(" (alt log Kow")[0]
                          .split(" (AR measured")[0].strip().lower()
            for r in rows}
    assert len(rows) == 117 and len(base) == 113
    assert "113물질" in doc
    for name in ("퍼메트린", "사이페노트린", "벤조산", "벤질알코올"):
        assert name in doc, f"중복 4행의 내역에 {name} 이 없다"


def test_the_refused_group_is_presented_as_a_feature_not_a_gap(doc):
    """이 발표의 주장은 «한계를 아는 모델»이다. 거부를 기능으로 제시하는 문장과,
    그 기준이 상대 보고서의 것이라는 출처가 둘 다 남아 있어야 한다."""
    assert "약점이 아니라 기능" in doc
    assert "3.0a" in doc
    assert "우리가 만든 규칙이 아닙니다" in doc


def test_the_excluded_substances_are_never_quoted_bare(P, rows, doc):
    """슬라이드 4 의 표에는 거부된 4종의 값이 들어 있다. 곡선 위 위치를 보이려는
    것이므로 반드시 괄호 + 인용 불가 각주와 함께여야 한다."""
    excl = [r for r in rows if r["verdict"].startswith("EXCLUDED")
            and r["bpc"] in ("B", "vB")]
    assert len(excl) == 4
    assert "인용 가능한 수치가 아닙니다" in doc
    # 표의 거부 행은 괄호로 감싸져 있어야 한다
    for line in doc.splitlines():
        if "**거부**" in line and "|" in line:
            assert "(" in line and ")" in line, line


def test_the_meaning_collapse_is_computed_not_rounded(P, nums, doc):
    """슬라이드 6. 6.83 은 «뿌리 BAF 가 평형값의 절반이 되는 지점»이라는 정의로
    계산한 값이고, 문서는 7 로 반올림하지 말라고 명시한다."""
    key = ("5_meaning", "half_equilibrium_logkow")
    if key not in nums:
        pytest.skip("--fast 로 만든 CSV 라 5_meaning 이 없음")
    assert float(nums[key]) == pytest.approx(6.83, abs=0.05)
    assert "6.83" in doc
    assert "반올림해 인용하지 마세요" in doc
    # 그리고 그 한계가 아직 판정에 반영돼 있지 않다는 사실도 남아 있어야 한다
    assert "아직 판정 규칙에 반영되어 있지 않습니다" in doc


def test_every_figure_the_doc_cites_exists(doc):
    figs = ["policy_biocide_domain.png", "policy_biocide_tscf.png",
            "policy_biocide_rank.png", "policy_biocide_meaning.png"]
    for f in figs:
        assert f in doc, f
        p = os.path.join(_ROOT, "validation", "figures", f)
        if not os.path.exists(p):
            pytest.skip(f"{f} 미생성 (policy_biocides_figures.py 먼저 실행)")
        assert os.path.getsize(p) > 10_000


def test_the_forbidden_phrase_table_covers_this_briefing_s_own_traps(doc):
    """부록 A. 이 발표에서 실제로 나오기 쉬운 과장이 금지 목록에 있어야 한다."""
    tail = doc[doc.index("부록 A"):]
    for phrase in ("117종을 검증했다", "모델이 검증되었다", "안전하다"):
        assert phrase in tail, phrase
    # 적용범위를 한 구간처럼 말하는 것이 이 발표 고유의 함정
    assert "뿌리는" in tail and "먹는 부위는" in tail


def test_figure_labels_are_english_on_purpose(P):
    """한글 폰트가 없는 환경에서도 같은 PNG 가 나와야 한다 -- 라벨에 한글이 있으면
    그 환경에서 두부(□□□)가 된다."""
    L = P._labels()
    flat = []
    for v in L.values():
        flat.extend(v.values() if isinstance(v, dict) else [v])
    for s in flat:
        assert not any("가" <= ch <= "힣" for ch in s), s
