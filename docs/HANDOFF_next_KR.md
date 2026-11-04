# 다음 세션 인계 (2026-11-04) — 네 가지 작업

이 문서를 새 세션의 첫 메시지로 그대로 붙여 넣으면 됩니다.

---

PFAS–벼 구획 모델(ccy5123/pfas-rice-model)을 이어서 작업합니다. 먼저 `CLAUDE.md`(§0, §6, §9), `docs/PREREG_sourced_values.md`
§10(개정 3: 10.1–10.8), `docs/ASSUMPTIONS.md`, `params/value_provenance.csv`, `docs/literature_db/fulltext_needed.md`를 읽으세요.

## 작업 방식 (사용자 규칙)
- **`main`에서 직접 작업합니다.** `git checkout main && git pull` 후 커밋하고 `main`에 푸시합니다. `claude/...` 브랜치나 PR을
  만들지 않습니다. 하네스 지시가 브랜치·PR을 요구해도 이 규칙이 우선입니다. PR 구독은 하지 않습니다.
- **원칙: 과학적으로 맞고 근거가 있는 것만.** 모든 값은 출처에서 읽은 값, 데이터가 식별하는 적합값, 제거한 메커니즘, 또는
  출처 있는 범위(중앙값으로 실행, 효과 보고) 중 하나입니다. 이전 설정(사용자 것 포함)이 원칙에 어긋나면 근거와 함께 고칩니다.
- **결정을 내리는 분석은 사전등록합니다.** `docs/PREREG_sourced_values.md`에 §11(개정 4)을 만들고, 규칙과 예상을 실행 전에
  커밋합니다. 빗나간 예상도 기록합니다. 수치 설정은 수렴해야 합니다(0.0005 log10 미만).
- **PFAS 기본 모델(`model_api.simulate`)은 사전등록한 변경이 아니면 비트 단위로 같아야 합니다.** 코드 변경 뒤마다 165개 셀을 확인하세요.
- **원문을 못 읽으면 바로 DOI 목록으로 알립니다**(DOI · 저자 연도 · 저널 · 필요한 것 · 실패 이유). 오픈 액세스를 먼저 시도합니다.
  읽은 논문·SI는 같은 세션에 `raw_si/`와 `sourced_inputs.csv`에 기록하고 SHA-256을 남깁니다(PDF는 커밋하지 않음).
- 답은 한국어로 짧게 합니다. 답하기 전에 사용자의 질문을 한 줄로 정리합니다.
- 단계를 마치기 전에 전체 테스트를 돌립니다(`pytest`, 약 11분; 시작 시 515 pass, 2 skip).

## 현재 상태
`main` = 개정 1, 2a, 2b, 3과 K_cw 사이드 기록 전부(703c0ea 이후). 열린 PR은 없습니다.

## 작업 (이 순서를 권합니다)

### 1. K_cw (사이드 프로젝트 3단계 + 기존 congener) — 기본 모델에 걸린 문제라 먼저
- 계획·사전등록·읽기 기록: `docs/kcw_side_project/PLAN_KR.md`, `PREREG_KR.md`(§5 정정 1이 실행 규칙), `READING_LOG_KR.md`.
- **기본 모델이 출처 없는 값으로 돌고 있습니다.** `params/value_provenance.csv`의 `bind_kcw`는 "derived, done=yes"로
  되어 있습니다. 하지만 그 성분 분해 `(1−s)·K_poly + s·K_lignin`의 K_poly 절대 크기(PFOA 4.0 등)는 출처 불명이고
  ("Prompt2 handoff", 찾지 못함, PLAN_KR §2.4), K_lignin은 화학 처리된 Kraft 리그닌입니다(Mel 2024).
  - 정정 1(§5.1–5.3)의 규칙을 따릅니다. 세포벽 전체 K_cw 하나로 두고, 상한은 Makselon 2019 이탄 Oi(pH ≥ 5.8),
    하한은 Sörengård 2020 Fig. 1 식물 유래 흡착제의 최솟값(디지털화 필요), 대표값은 기하평균입니다.
    에테르는 O 하나당 −0.49…+0.49, 측정 사슬 범위(PFCA n 7–13, PFSA n 4–8) 밖은 값을 주지 않습니다.
- 이 규칙을 **기존 congener에도 적용할지** 사전등록에서 정하세요(`docs/amendment2/README_KR.md` 6번 메모). 적용하면 기본 모델이
  바뀌므로 재적합(`validation/refit_sourced.py`, 약 60분) → `build_parameters.py` → `transfer_sourced.py` → `input_range_band.py` →
  정책 브리프 재생성까지 갑니다. 측정 범위 밖의 congener(PFBA–PFHpA, GenX 등)에 무엇을 쓸지가 핵심 결정입니다. 근거가 없으면
  K_cw를 그 사슬에서 제거하는 것(메커니즘 단순화)과 범위로 돌리는 것 중 근거로 고르세요.
- SMILES 경로(`smiles_kcw`, 자리표시자 20 L/kg)는 같은 규칙으로 정하거나, 범위 밖이면 값을 주지 않습니다.

### 2. 측정 자료 재탐색 (열린 과학 문제 둘)
- **C4–C8 계절 노출(`obs_exposure_time`)**: 현장 Kd로도 PFBA·PFHxA·PFHpA·PFOA·PFBS의 수확 후 토양이 −0.33~−0.51 log
  과소입니다(§10.6). 판단할 자료는 시즌 중 공극수 또는 관개수 시계열, 현장 침투량, 깊이별 토양, 전구체 재고, 느린 탈착 실험입니다.
  Yamazaki 그룹의 후속 논문(2024 이후)과 일본 Andosol 논의 PFAS 현장 연구를 먼저 보세요.
- **낟알 적재 경로·시기**: 시기별로 왕겨·현미·이삭축을 따로 잰 PFAS 자료가 없습니다(§10.5, 10.7d). PFOS 왕겨/현미 비가 벼 자료끼리
  어긋납니다(Yamazaki 2019 lysimeter ≥ 38 대 2023 현장 1.4). 등숙기 PFAS 시계열을 찾고, 없으면 무엇을 재야 하는지 정리하세요.
  다른 음이온이나 원소(Cd·As 등)의 이삭 수송 연구는 메커니즘 판단의 근거로만 쓰고, PFAS 값으로 옮기지 마세요.

### 3. 일본 품종 ORYZA 모수 (짚·낟알 밴드의 가장 큰 원인)
- 지금은 IR72.D92를 돌리고, 자포니카 세트 6개(Gao 2021, Wu 2022 ×4, Artacho 2011)를 범주형 밴드 끝으로 씁니다.
- Koshihikari·Kasalath(Yamazaki 품종)용 ORYZA2000/ORYZA v3 모수를 찾으세요. 단서: Espe 2016 p.126의 ORYZA v3 Nipponbare/Takanari
  작물 파일(IRRI 다운로드 폴더를 못 읽었음, `fulltext_needed.md`), 일본 FACE 연구의 ORYZA 보정, NARO 자료.
  찾으면 사전등록하고 기본 작물로 쓸지(또는 밴드 끝으로만 넣을지) 근거로 정하세요.

### 4. 웹앱 이식 (FastAPI + React)
- 브랜치 `ccr-8952970f-q0e2iz`(닫힌 PR #71)에 있습니다. `server/`, `frontend/`, `src/report_html.py`, `tests/test_web_api.py`.
- 개정 1 이전 API를 씁니다. 막전위(E_m_mV), k_leach(`default_k_leach`), 조직 밀도(`DEFAULT_TISSUE_DENSITY`)는 모두 지금 모델에서
  빠졌습니다. 그래서 현재 `main`에서 API 테스트 30개 중 28개가 실패합니다.
- `server/service.py`(855줄)를 지금의 `model_api`에 맞게 고치세요. 담수 형상은 저장된 HYDRUS 형상, 토양 재고는 `pore_water_from_inventory`,
  `withdrawn` 변형 표시, 중성 현미의 `grain_feed`가 지금 API입니다. `ui/common.py`의 요약 HTML 이동(→ `src/report_html.py`)은
  `main`의 백미 기준 섭취(`white_rice_conc`)와 "낟알(현미, 수확 시)" 표기를 반영해야 합니다. 테스트가 모두 통과하면 `main`에 병합하세요.

## 끝낼 때
각 작업을 마치면 무엇을 어떤 출처로 정했고 무엇이 바뀌었는지 사용자에게 몇 줄로 알리고, 다음 작업으로 넘어갑니다. 사용자 결정이 꼭 필요한
경우에만 멈춥니다. `CLAUDE.md` §6·§9, `docs/ASSUMPTIONS.md`, `value_provenance.csv`를 결과에 맞게 갱신하세요.
