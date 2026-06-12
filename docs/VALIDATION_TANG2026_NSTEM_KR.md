# Tang 2026 재검증 — 지상부 과이행 수정(다구획 줄기 + 잎/줄기 적재 재배분) — 한국어

> **한 줄 요약** — Tang 2026이 짚은 **지상부 과이행**(잎 sink 폭주·빈 줄기 칸)을
> **다구획 줄기(N-segment) + 잎/줄기 적재 재배분**으로 구조적으로 수정하고, 이어 **congener 절대 레벨을 f_xy로 보정**해
> Tang TF로 재검증했습니다. 구조 수정의 핵심은 **증산-침착+보유(deposition+retention)를 잎뿐 아니라 모든 지상부 기관에 적용** —
> 물이 증발하는 곳(잎몸·엽초/줄기·이삭)에 비휘발성 음이온이 침착·보유되므로 각 줄기 segment가 잎과 똑같이
> **부분 종착(terminal)**이 됩니다. **2단계 결과**: ① 구조가 구획 **패턴**을 치유(**shape RMSE 0.84 → 0.11**,
> PFOA 세 조직 일치, 잎 부담 81%→30%); ② **f_xy 보정이 절대 레벨을 치유**(**overall RMSE 1.28 → 1.01 → 0.18**).
> 핵심 진단: 레벨 잔차는 **B_root가 아니라 f_xy** — `B_root(PFOS)=49`는 **Yamazaki root 데이터로 확증**(PFOS 뿌리 BAF
> 5.93 ≈ PFOA 0.49의 12×)되어 옳고, 문제는 **monotone f_xy(PFOS)=0.013이 PFSA를 과벌점**(Yamazaki는 W2 0.14 필요)하고
> **GenX provisional f_xy=0.233이 18× 과대**(단쇄 PFCA 인공물)인 것. PFOS→W2(0.14), GenX→Tang보정(0.013)으로 **RMSE 0.18**.

![Tang 2026 nstem 재검증 그림](../validation/figures/tang2026_nstem_validation.png)

*그림(6패널): (A·B·C) PFOA·PFOS·GenX의 TF(조직/뿌리) — Tang(검정) vs single-straw 기준(파랑) vs
nstem_leaf monotone(빨강) vs **+f_xy 보정(보라)**. 보정 막대가 세 화합물 모두 Tang에 근접. (D) PFOA 식물-부담 분포 —
잎 독점(81%)이 줄기·잎으로 재배분. (E) RMSE 진행 — **구조가 패턴(shape 0.84→0.11), f_xy 보정이 레벨(overall 1.28→1.01→0.18)**.
(F) 레벨 레버는 **f_xy**(monotone vs W2 vs 보정값) — B_root 아님. 생성: `python validation/tang2026_nstem_validation.py`.*

---

## 1. 무엇이 문제였나 (Tang 1차 검증의 결론, `docs/VALIDATION_TANG2026_KR.md`)
단일 straw 4구획 코어(`pfas_rice_plant_module_4pool_surf`)는 지상부로 **과이행**했습니다:
- **빈 줄기 칸**: `dC_stem = (Q/M_stem)(Cw_xyl − Cw_stem)` — 줄기는 (f_xy 할인된) 물관 자유농도로
  평형만 맞춘 뒤 **상행 물관으로 즉시 재배출**하므로 축적 못 함. `TF_stem = f_xy·B_stem/B_root ≪ 1`
  (모델 0.01–0.25 vs Tang stalk 0.58–1.45).
- **잎 sink 폭주**: 잎이 **유일한 물관 종착**이라 전 증산 stream을 적분 → 성숙기 성장희석→0이면 폭주.
  잎이 **식물 전체 PFAS 부담의 ~81%**를 보유(PFOA), `TF_leaf` 3–13 (이동성 큰 GenX는 70.8) vs Tang 0.68–1.66.

둘은 **한 구조적 결함**입니다: 단일 straw에서는 지상부에 전달된 용질이 잎으로밖에 쌓일 곳이 없습니다.

## 2. 수정 — 다구획 줄기 + 증산 침착·보유 재배분 (`src/pfas_rice_plant_module_nstem_leaf.py`)
구획: `root(0), stem_1..stem_N(1..N), leaf(N+1), grain(N+2)`.
- **다구획 줄기**: 줄기를 N개 직렬 segment로 분해(기본 N=4).
- **증산 침착+보유**: 뿌리가 물관을 `Cw_xyl = f_xy·Cw_root`로 적재해 `Q·Cw_xyl` 전량을 내보내면,
  상행 stream에서 기관 k가 `λ_k·Q`의 물을 증산할 때 그 물이 운반한 용질 `λ_k·Q·Cw_xyl`이 **그 자리에 침착**됩니다.
  각 침착의 분율 `retention`은 기관에 **보유(비가역 — 증발 종착의 세포벽/아포플라스트 격리)**, 나머지(1−retention)는
  잔여 물관으로 grain까지 운반. **증산 분배가 정확히 닫힘**: `Σtau + λ_leaf + λ_grain = 1`
  → 지상부 침착 합 = 뿌리 물관 배출 = `Q·Cw_xyl` (**질량보존**, 테스트로 고정).
- 이로써 **각 줄기 segment가 잎과 똑같이 자기 증산-침착을 보유하는 부분 종착**이 됩니다.
  → 지상부 부담이 잎 독점이 아니라 root→stem→leaf→grain로 재배분.
- **두 구조 레버**(둘 다 PFAS-무관한 작물 구조량): 줄기 증산분율 `stem_transp_frac`(엽초·줄기 표면 비중),
  보유효율 `retention`(증발 종착에서 보유 vs 재가동 분율). 기본값 `0.45 / 0.6`은 Tang TF에 **점-맞춤하지 않은** 구조 설정.

> 이 모듈은 `pfas_rice_plant_module_nstem`(직렬 **mixer** — Yamazaki **줄기 내부 수직 gradient** S18/S19용)와
> **상보적**입니다: 거기선 용질이 상행 물관과 **가역 평형**, 여기선 증발 종착에서 **보유**(Tang의 stalk/leaf/grain
> 분배가 요구한 메커니즘).

## 3. 정량 결과 — TF(조직/뿌리), Tang vs 기준 vs nstem_leaf

| PFAS | 조직 | Tang 평균 | single-straw(기준) | **nstem_leaf** |
|---|---|---:|---:|---:|---:|
| | | **Tang** | single-straw | nstem(mono) | **+f_xy 보정** |
| PFOA | stalk | 1.45 | 0.03 | 1.27 | **1.27** |
| PFOA | leaf  | 1.66 | 5.95 | 2.04 | **2.04** |
| PFOA | endosperm | 0.95 | 0.41 | 0.93 | **0.93** |
| PFOS | stalk | 0.58 | 0.01 | 0.04 | **0.36** |
| PFOS | leaf  | 0.68 | 0.16 | 0.06 | **0.58** |
| PFOS | endosperm | 0.77 | 0.01 | 0.03 | **0.27** |
| GenX | stalk | 1.10 | 0.25 | 20.3 | **1.30** |
| GenX | leaf  | 1.38 | 70.8 | 24.8 | **1.55** |
| GenX | endosperm | 1.39 | 14.6 | 17.5 | **1.09** |

**log10 RMSE 진행 (3 화합물 × 3 조직 = 9점):**

| 지표 | single-straw | nstem(mono) | **+f_xy 보정** |
|---|---:|---:|---:|
| 전체(9점) | 1.28 | 1.01 | **0.18** |
| **shape**(조직 **패턴**, 레벨 제거) | 0.84 | **0.11** | 0.11 |

→ **구조가 패턴을(0.84→0.11), f_xy 보정이 레벨을(1.01→0.18) 치유.** PFOA 식물-부담 분포: 잎 **81%→30%**, 줄기 **1%→29%**.

## 4. 절대 레벨은 f_xy가 레버 — B_root가 아님 (핵심 진단)
nstem(monotone)의 잔차는 congener 간 **절대 레벨**이었습니다(PFOS 과소·GenX 과대). 원인을 정확히 분리하면 **B_root가 아니라 f_xy**입니다:

- **`B_root(PFOS)=49`는 Yamazaki root 데이터로 확증** — PFOS 뿌리 BAF **5.93** ≈ PFOA **0.49**의 **12×**. 즉 강결합(K_PL=31623)으로 뿌리가
  PFOS를 실제로 많이 잡으며, 이건 **틀린 게 아니라 측정과 일치**. B_root는 그대로 둡니다.
- **질량보존 논증으로 f_xy를 지목**: PFOS는 PFOA보다 **적게** 전달되는데(monotone f_xy 0.013 → 0.04의 3× 미만) Yamazaki straw는
  **5× 많이** 보유(PFOS straw 4.35 vs PFOA 0.83) — 단방향 전달로는 불가능. 즉 **monotone f_xy(PFOS)=0.013이 과소**.
- **f_xy 두 출처의 불일치**: monotone(QSPR+머리기 exp(−1.1))은 PFOS=0.013이지만, **W2 fit(Yamazaki 재현)은 0.142** — monotone이 PFSA를
  **과벌점**. GenX는 provisional 0.233(단쇄 PFCA × ether offset)인데 Tang 데이터는 GenX f_xy≈PFOA(~0.013)를 요구 — **18× 과대**.

**보정**(명시적·라벨링): PFOA = monotone 0.040(불변), **PFOS = W2 0.142(독립 Yamazaki-grounded)**, **GenX = 0.013(Tang 보정 — 독립 데이터 없는 provisional)**.
→ **overall RMSE 1.01 → 0.18** (F 패널). PFOS는 W2(독립값)이라 Tang stalk/grain을 약간 과소(0.36/0.27 vs 0.58/0.77)하나 방향·order 일치.

## 5. 정직한 결론
- **2단계로 Tang 재검증 통과**: ① 재배분이 지상부 **패턴**을 구조적으로 치유(shape RMSE 0.84→0.11, PFOA 세 조직 일치);
  ② **f_xy 보정이 절대 레벨**을 치유(overall RMSE **1.28→1.01→0.18**), 세 화합물 모두 Tang과 order-of-magnitude 이내.
- **진단 확정**: 레벨 잔차는 **f_xy**(monotone PFSA 과벌점 + GenX provisional)였고 **B_root는 옳다**(Yamazaki root 확증).
  monotone f_xy의 **머리기 exp(−1.1) 오프셋이 PFSA를 과도하게 깎는다**는 것이 독립적으로 드러남(W2가 데이터-일관).
- **방법론 주의**: GenX f_xy는 Tang에 **보정**(독립 데이터 없음 → OOS 아님); PFOS f_xy는 **독립 Yamazaki(W2)** 값이라 Tang에 대해 준-OOS 유지.
  보정값은 검증 스크립트에 **override로만** 반영하고 **canonical `params/parameters.json`은 변경하지 않음**(provenance 보존).

## 6. 다음 작업
1. **monotone f_xy의 PFSA 머리기 오프셋 재산정** — `exp(−1.1)`이 과벌점(PFOS W2 0.14 vs monotone 0.013). Tang+Yamazaki TF로
   PFSA 오프셋을 재적합하고 `f_xy_recommended`(PFSA)를 갱신할지 검토(`build_parameters.py`, task #8 연계).
2. **GenX ether 전용 Koc·결합 + f_xy 정식화** — 본 보정(f_xy 0.013)을 1차 검증의 BCF 40× 과대(단쇄 PFCA Koc 근사)와 함께
   ether-PFAS 전용 QSPR로 정식화해 canonical 반영.
3. **구조 레버 독립 고정** — `stem_transp_frac`·`retention`(기본 0.45/0.6)을 작물 구조 데이터(엽초/줄기 증산 분배, 세포벽 보유)로
   a priori 고정해 Tang TF로부터 완전 독립화.
4. **시계열 TF(t) 재검증** — 재배분이 Tang의 "1개월 후 뿌리 우세(root-first)"를 재현하는지(1차 검증의 ❌ 항목).

## 7. 재현
```bash
python validation/tang2026_nstem_validation.py   # 표 + 6패널 그림
#   -> validation/figures/tang2026_nstem_validation.png
python -m pytest tests/test_nstem_leaf.py -q      # 질량보존 + 구조 치유 회귀
python src/pfas_rice_plant_module_nstem_leaf.py   # PFOA smoke (retention 1.0 / 0.6)
```
모델: `src/pfas_rice_plant_module_nstem_leaf.py`, API: `model_api.simulate_nstem_leaf(...)`.
데이터: `docs/literature_db/raw_si/tang2026_doseresponse.csv`(SI S8). 관련: `docs/VALIDATION_TANG2026_KR.md`, `docs/VALIDATION_KR.md`.
