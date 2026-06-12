# Tang 2026 재검증 — 지상부 과이행 수정(다구획 줄기 + 잎/줄기 적재 재배분) — 한국어

> **한 줄 요약** — Tang 2026이 짚은 **지상부 과이행**(잎 sink 폭주·빈 줄기 칸)을
> **다구획 줄기(N-segment) + 잎/줄기 적재 재배분**으로 구조적으로 수정하고 Tang TF로 재검증했습니다.
> 핵심 메커니즘은 **증산-침착+보유(deposition+retention)를 잎뿐 아니라 모든 지상부 기관에 적용** —
> 물이 증발하는 곳(잎몸·엽초/줄기·이삭)에 비휘발성 음이온이 침착·보유되므로 각 줄기 segment가 잎과 똑같이
> **부분 종착(terminal)**이 됩니다. **결과**: 구획 **패턴**(stalk/leaf/grain 관계)이 결정적으로 개선 —
> **shape RMSE 0.84 → 0.11**, **PFOA 세 조직 모두 Tang과 일치**(stalk 0.03→1.27, leaf 5.95→2.04,
> grain 0.41→0.93; RMSE 1.03→0.06), **PFOA 잎 부담 81%→30%·줄기 1%→29%로 재배분**.
> **남은 것**: congener 간 **절대 레벨** — Yamazaki 보정 `f_xy`(PFOS 0.013, GenX 0.233)와 basis-A
> `B_root`(PFOS 49)의 스프레드 때문에 **PFOS는 뿌리에 과결합(TF 과소)·GenX는 과이동(TF 과대)**.
> 이건 지상부 구조가 아니라 **결합/이행 magnitude 보정** 문제로, 진단이 한 단계 좁혀졌습니다.

![Tang 2026 nstem 재검증 그림](../validation/figures/tang2026_nstem_validation.png)

*그림(6패널): (A·B·C) PFOA·PFOS·GenX의 TF(조직/뿌리) — Tang(검정) vs **single-straw 기준(파랑)** vs
**nstem_leaf 재배분(빨강)**. (D) PFOA 식물-부담 분포 — 잎 독점(81%)이 줄기·잎으로 재배분. (E) RMSE 분해 —
지상부 **패턴(shape)** 0.84→0.11로 치유, 잔차는 congener 간 레벨. (F) 보유효율(retention) 스윕 —
PFOA 세 조직이 retention≈0.6에서 Tang(점선)에 수렴. 생성: `python validation/tang2026_nstem_validation.py`.*

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
|---|---|---:|---:|---:|
| PFOA | stalk | 1.45 | 0.03 | **1.27** |
| PFOA | leaf  | 1.66 | 5.95 | **2.04** |
| PFOA | endosperm | 0.95 | 0.41 | **0.93** |
| PFOS | stalk | 0.58 | 0.01 | 0.04 |
| PFOS | leaf  | 0.68 | 0.16 | 0.06 |
| PFOS | endosperm | 0.77 | 0.01 | 0.03 |
| GenX | stalk | 1.10 | 0.25 | 20.3 |
| GenX | leaf  | 1.38 | 70.8 | 24.8 |
| GenX | endosperm | 1.39 | 14.6 | 17.5 |

**log10 RMSE 분해 (3 화합물 × 3 조직):**

| 지표 | single-straw | **nstem_leaf** |
|---|---:|---:|
| **shape**(congener 내 조직 **패턴**, 레벨 제거) | 0.84 | **0.11** |
| 전체(9점) | 1.28 | 1.01 |
| **PFOA** | 1.03 | **0.06** |
| PFOS | 1.56 | 1.27 |
| GenX | 1.21 | 1.21 |

**PFOA 식물-부담 분포**: 잎 **81% → 30%**, 줄기 **1% → 29%** (D 패널).

## 4. 읽는 법
- ✅ **지상부 패턴이 치유됨** — 재배분이 노린 차원(조직 간 관계)의 **shape RMSE 0.84 → 0.11**(7.6× ↓).
  빈 줄기(stalk≈0)와 잎 폭주(leaf≫1)가 사라지고 **stalk≈leaf≈grain≈O(1)** 라는 Tang의 패턴을 재현.
- ✅ **PFOA(중심 congener)는 세 조직 모두 일치** — 1.27 / 2.04 / 0.93 vs Tang 1.45 / 1.66 / 0.95 (RMSE **1.03→0.06**).
  `retention≈0.6`에서 보유분이 줄기·잎을, 미보유분(40%)이 grain을 채워 셋이 동시에 맞음(F 패널).
- ✅ **질량보존 정확**(모든 retention에서 유일 소스 `M_root·j_R`; `tests/test_nstem_leaf.py`).
- ❌ **congener 간 절대 레벨은 여전히 빗나감** — PFOS 과소(0.04 vs 0.58), GenX 과대(20.3 vs 1.10).
  원인은 **지상부 구조가 아니라** `TF ∝ f_xy/B_root × (증산/질량)`의 **f_xy·B_root 스프레드**:
  - PFOS: `f_xy=0.013`(작음) **+ `B_root=49`(K_PL=31623, 인지질 결합 큼 → 뿌리 과결합)** → TF 14× 과소.
  - GenX: `f_xy=0.233`(큼) + `B_root=1.3`(작음) → TF 18× 과대.
  머리기 **비율**(PFSA/PFCA≈0.4)은 Tang(0.58/1.45=0.40)과 일치하나, **절대 스케일**이 Yamazaki regime이라 Tang보다 넓음.

## 5. 정직한 결론
- **재배분은 Tang이 짚은 지상부 과이행을 구조적으로 치유**합니다 — 빈 줄기·잎 폭주가 사라지고
  조직 **패턴**이 Tang과 일치(shape RMSE 0.84→0.11), **PFOA는 세 조직 정량 일치**.
- **남은 오차는 한 차원으로 좁혀짐**: congener 간 **절대 레벨**(`f_xy`·`B_root` magnitude).
  진단이 "지상부 구획 분배가 틀림"에서 → "**결합/이행 magnitude가 Yamazaki regime이라 Tang보다 스프레드가 넓음**"으로 이동.
- 이는 **구조 작업이 아니라 보정 작업**입니다(§6). nstem_leaf로 지상부 구조는 더 이상 병목이 아닙니다.

## 6. 다음 작업
1. **congener 절대 레벨 보정** — Tang의 좁은 TF 범위는 `f_xy`·`B_root` 스프레드가 과하다고 시사:
   (a) **PFOS `B_root` 재검토** — K_PL(인지질) 기반 뿌리 결합이 뿌리 잔류를 과대평가하는지(Tang에서 PFOS는 지상부 도달 양호);
   (b) **Tang regime용 `f_xy` 절대 스케일** — Tang은 Yamazaki보다 이행이 더 균일(congener-uniform).
   둘 다 Tang/Kim에 **명시적 보정**으로 하되 OOS와 분리해 보고.
2. **GenX ether 전용 Koc·결합** — 1차 검증의 BCF 40× 과대와 연결(단쇄 PFCA 근사 → 토양·식물 결합 과소).
3. **구조 레버 독립 고정** — `stem_transp_frac`·`retention`을 작물 구조 데이터(엽초/줄기 증산 분배,
   세포벽 보유)로 a priori 고정해 Tang TF로부터 완전 독립화(현재 기본값은 점-맞춤 아님이나 미세조정 여지 있음).
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
