# 개정 2 인계 (다음 세션용)

개정 1(출처 있는 값만 사용)과 부록 5.7–5.10은 끝났고 결과는 `docs/PREREG_sourced_values.md` §6에 있습니다.

**진행 상황 (2026-10-26)**: 개정 2는 2a/2b로 나눠 등록했습니다. **2a 완료** — 사전 등록 §7(7.1–7.13), 구현, 재적합,
전이 시험, 밴드, 결과 §8. 기본 구성은 규칙에 따라 `redistributed_fit`으로 바뀌었습니다(사용자 적용). 아래 목록의
2–5, 11(2a 범위의 문서 갱신)은 끝났고, **2b(§7.9)가 남았습니다**: 6(SMILES K_cw), 7(현장 토양·관개수 시나리오),
8(중성 경로), 9(휘발성), 10(앱 노출 기본값). 원문을 못 구한 논문은 `docs/literature_db/fulltext_needed.md`
(10.1139/er-2025-0116 대기 중).

## 2b 시작 방법 (다음 세션: "개정 2b 시작")

- **브랜치**: 2a는 `claude/clever-carson-yktrl0` (PR ccy5123/pfas-rice-model#72, 초안, CI 통과, 마지막 fa57547)이고
  개정 1 PR(#70, `claude/confident-newton-82t0qo`) 위에 쌓여 있습니다. 둘 다 아직 병합 전이므로 2b는 2a 브랜치 위에서
  시작합니다(새 세션에 지정되는 브랜치를 2a 머리에서 만들기).
- **첫 단계 = 사전 등록**: `docs/PREREG_sourced_values.md` §7.9의 2b 항목을 §9(개정 2b)로 풀어 쓰고, 범위·규칙·예상(E25~)을
  실행 전에 커밋합니다. 2a에서 정한 결정 규칙(§2, 5.3)과 기본 변형(`redistributed_fit`)은 그대로입니다.
- **남은 출처 없는 값** (`params/value_provenance.csv`의 `done=no`):
  - 중성 경로: `neu_lipid_a` 1.22, `neu_leafexp` 0.95, `neu_surface` 20/2, `neu_kaw_cut` 1e-4, `mem_ph` 6.5
    (→ 논 공극수 pH 6.5–7.0 출처). 중성 경로는 아직 상수 수분·이삭 전체 기준(`model_api.neutral_forms`)이라
    θ(t)·현미 기준 적용 여부를 함께 정합니다.
  - 토양: `soil_foc`, `soil_rho`, `soil_depth`, `soil_disp`, `soil_perc`, `soil_flood` (Yamazaki 현장 Andosol;
    f_oc 3.7–3.8 %, 용적밀도 0.61, 침투 0.41–0.67 cm/d, 중간 낙수 — 출처 확인 필요), 관개수 오염원 시나리오
    (`obs_exposure_time`).
  - 앱: `app_density`, `smiles_kcw`(프로토타입 0004), KNHANES 섭취량·체중, 백미/현미 기준.
  - 분석: `tr_drivers`, `tr_li_grain_from`, `calib_sigma`, `fxy_recommended`.
- **원문 대기**: 10.1139/er-2025-0116 (Malaviya 외 2025, Environ. Rev.) — 사용자가 올려 주기로 함.
- **작업 규칙**: 원문을 못 읽으면 즉시 DOI 목록으로 알림(CLAUDE.md §0). 재적합은 약 75분
  (`validation/refit_sourced.py --write` → `build_parameters.py` → `transfer_sourced.py` → `input_range_band.py`),
  전체 테스트 약 9분(503개).

아래는 개정 2 시작 전에 남긴 원래 작업 기록입니다.

- `draft.md` — 개정 2 초안(영문). 항목, 찾은 문헌, 외부 검토(P1–P9) 확인 결과, θ(t) 설계, P3(현미 기준 낟알) 문헌.
  **사전 등록문이 아닙니다**: 등록할 때 이 초안에서 문장을 골라 `PREREG_sourced_values.md` §7로 옮기고, 실행 전에 커밋합니다.
- `prototypes/*.patch` — 이전 세션 scratchpad에서 만든 프로토타입 4개(`git format-patch`). 개정 1 완료 **전** 코드
  위에서 만든 것이라 그대로 `git am` 되지 않을 수 있습니다. 설계 참고용입니다.
  1. 기관 수분 θ(DVS) — 단계별 출처 앵커 사이 선형(`params/rice_theta_stage_sourced.csv`), 생체중·결합인자가
     같은 θ(t)를 읽고, 낟알 체관 흐름은 건물 생장, 수확 연산자는 수확 시점 건물로 나눔.
  2. 밴드·θ 테스트를 단계 앵커에 맞춤.
  3. 공개된 자포니카 작물 세트(Gao 2021 헤이룽장, Wu 2022 랴오닝)를 범주형 밴드 입력으로.
  4. SMILES 입력의 K_cw를 GAP A 사슬길이 사다리에서 유도, SMILES Compound에서 출처 없는 수송 상수 제거.

## 개정 2에서 할 일 (순서 제안)

1. **사전 등록** (§7): 아래 항목, 범위, 재적합 규칙, 예상(E-a, E-b, E-c 등)을 실행 전에 커밋.
2. **θ(t) 생육단계별 수분** (외부 검토 P2) — 재적합 필요(뿌리 수확 수분을 `r_eq`가 흡수).
3. **낟알 = 현미 질량** (P3) — WSO × (1 − 왕겨 비율(t) − 이삭축 비율). 왕겨 비율은 출처 있음(Yoshida 1981,
   Mizuno 1987, Ansari 2016), 이삭축 비율은 출처 없음 → 물리적 범위 또는 왕겨와 함께 이삭 풀에 둠.
4. **실측 표 경로의 잎 손실·보간** (P4, P5), **재식밀도 상수 하나 + 불변성 테스트** (P7), 문서화(P8 plaque 비교
   기준, P9 출수 후 줄기 질량 감소).
5. **자포니카 작물 세트를 밴드로** (IR72와 함께).
6. **SMILES K_cw** (프로토타입 4).
7. **현장 토양** — 관개수가 오염원인 시나리오(지금 HYDRUS 시나리오는 오염 토양 + 깨끗한 관개수로 반대),
   f_oc 3.7–3.8 %, 용적밀도 0.61, 침투 0.41–0.67 cm/d, 중간 낙수, 수리 특성(van Genuchten)은 출처 탐색 필요.
8. **중성 경로의 남은 출처 없는 값** (`params/value_provenance.csv` neu_*): 지질–옥탄올 보정 a = 1.22,
   잎 흡착 지수 0.95, 잎·낟알 표면적 20/2, 휘발성 기준 1e-4; Hwang 2017 상추 수분(진단용).
   부록 5.10에서 중성 조성은 자료별 출처 값으로 옮겼으므로, θ(t)를 중성 경로에도 적용할지 함께 정함.
9. **휘발성 화합물** 검증 자료 탐색 (지금은 민감도 분석만).
10. **앱 노출 기본값**: 쌀 섭취량·체중(KNHANES), 백미/현미 기준 차이(draft.md "BASIS BUG").
11. **최종 재적합** 후 문서 일괄 갱신(CLAUDE.md §6–9, README, MANUAL, OVERVIEW, 브리핑).

## 시작할 때 확인할 것

- `params/value_provenance.csv`에서 `done=no` 행 = 남은 출처 없는 값 목록.
- 결정 규칙은 개정 1과 같음: 값은 출처 / 데이터가 식별하는 적합 / 제거 / 출처 있는 범위(중앙값 실행 + 밴드).
- 전체 테스트 약 8분(`pytest`), 재적합 약 30분(`validation/refit_sourced.py --resume` 가능).
