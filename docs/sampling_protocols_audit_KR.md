# 실측 데이터 분취·전처리 감사 + 기관 질량 기준(건물중→생중량) 수정

> 질문: "우리 모델은 어떤 성장모델을 쓰나? 기관별 질량 실측 자료는? 우리가 쓴 실측값은
> 실제 실험에서 각 기관을 어떻게 나눠 분석했나(낟알 도정 여부, 뿌리 세척 등)?"
>
> 이 문서는 그 조사의 기록이다. 원문(본문·SI)을 직접 읽어 확인한 내용만 적었고, 원문에 없는
> 것은 "미보고(NR)"로, 원문을 보지 못한 것은 "미확인"으로 구분했다. 데이터셋×기관별 표는
> `docs/literature_db/Sampling_Protocols.csv`, 문헌의 기관 질량 수치는
> `docs/literature_db/raw_si/rice_organ_biomass_literature.csv`에 있다.
> 원문 PDF는 사용자가 제공했으며 저작권 때문에 커밋하지 않았다.

## 0. 요약

1. **성장모델**: 기본 바이오매스 드라이버는 **ORYZA2000 Level-1의 Python 재구현**
   (`src/oryza_growth.py`, IR72 파라미터)이다. 레거시 `growth_rice`(로지스틱 + IR72 분배),
   실측표 입력 `measured_biomass`가 있다.
2. **버그(이번에 수정)**: 세 드라이버 모두 **건물중(dry matter)** 을 내는데, 식물 ODE의
   `PlantInputs.M`과 결합계수 `B_k`는 **생중량(fresh)** 기준이다. 변환 없이 들어가 있었다.
   → `src/organ_mass.py`가 `fresh = dry/(1−θ_fw)`로 변환한다(θ_fw는 `B_k`와 같은 값).
   기본값 `mass_basis="fresh"`, 수정 전 동작은 `mass_basis="dry"`로 재현된다(§1).
3. **보정·검증 데이터의 실제 실험 조건이 레포 서술과 달랐다**: Yamazaki 2023은 온실이 아닌
   **노지 논**, Li 2025의 TF는 저자 정의와 다른 방식으로 계산됨, Kim 2019의 PFOA BAF는
   **한 필지가 공극수 평균의 83%** 를 차지(§2, §4).
4. **"낟알"이 데이터셋마다 다른 부위**다(현미 재구성 / 현미 / 백미 추정 / 미정의). 모델은 현미 조성.
5. **뿌리**: 어떤 데이터셋도 Fe/Mn plaque를 제거하지 않았다. 표면 탈착 단계가 있는 것은 Liu 2023(아세토니트릴 헹굼)뿐.
6. **기관 수분 함량(θ_fw)** 을 보고한 벼 자료는 없다. 모든 dw↔fw 환산은 가정값에 기댄다.

---

## 1. 성장모델과 기관 질량의 기준

### 1.1 드라이버
| 드라이버 | 용도 | 비고 |
|---|---|---|
| `oryza_growth` (ORYZA2000 Level-1, IR72) | **기본** (`biomass="oryza"`) | `*_oryza` 적합값, `TANG_REFIT_FXY`의 기반 |
| `growth_rice` (로지스틱 + IR72 분배) | 레거시 | `LIPID_LOADING` 상수, two-pool 실측forcing 적합의 기반 |
| placeholder 로지스틱 | `measured_forcing=False`, `reproduce_demo` | W2 적합(RMSE 0.029)이 성립하는 유일한 드라이버. **이번 수정의 영향 없음** |
| `measured_biomass` | 실측 기관 **건물중** 표 → M(t) | 이번 수정으로 생중량 변환 |

### 1.2 건물중/생중량 버그와 수정
- 근거: `organ_biomass_oryza`/`growth_rice.organ_biomass`는 "dry biomass [kg/hill]"을 반환하고,
  `model_api._default_drivers`는 그것을 그대로 `M`에 넣었다. `PlantInputs.M`은 "tissue fresh mass",
  `B_k = θ + (1−θ)·Σf·K` [L/kg fw].
- 영향: 말단 기관 농도(=부담량/M)가 **건물 1 kg당 부담량을 생중량 1 kg당으로 보고**한 셈이었다.
  잎 ×4.5, 줄기 ×5.9, 낟알 ×1.16. 뿌리는 흡수량이 질량에 비례하므로 거의 변하지 않는다.
- 수정: `src/organ_mass.py` — `MASS_BASES = ("fresh", "dry")`, 기본 `"fresh"`.
  적용 지점: `model_api`(`simulate`, `simulate_nstem_leaf`, `tang_tf_validation`, `simulate_neutral`,
  `measured_forcing`, `drivers_from_arrays`, `hydrus_drivers`), `soil_hydrus.inputs_from_hydrus`,
  `oryza_growth.oryza_drivers`, `measured_biomass.biomass_drivers`, `growth_rice.M_for_nstem`,
  `pfas_rice_two_pool`, 앱의 성장표 기본값(`ui/common.py`, "생중량" 표인데 건물중을 채우던 문제),
  `validation/neutral_dpu_validation.py`, `validation/structural_adequacy_fit.py`.
  사용자가 직접 넣는 `drivers=`의 `M`과 앱의 성장표는 원래 생중량 계약이므로 그대로 쓴다.
- `"dry"`는 물리적 대안이 아니라 **수정 전 결과 재현용**이다.

### 1.3 수정 전후 주요 수치 (같은 코드, `mass_basis`만 다름)
| 항목 | dry(수정 전, 문서값) | fresh(수정 후) | 해석 |
|---|---|---|---|
| PFOA 기본 BAF 뿌리/잎/낟알 | 0.479 / 2.52 / 0.148 | 0.488 / 0.614 / 0.036 | 지상부 4–6배 하락 |
| PFOS 기본 BAF 뿌리/잎/낟알 | 6.02 / 0.781 / 0.048 | 6.06 / 0.122 / 0.008 | |
| **Ge 2017 (중성, 적합 없음) log10 RMSE** | **0.783** | **0.311** | 잎 과대 6–14배 → 1.6–3.7배 |
| Ge 2017 반감기 스캔 최적 | ~7 d (0.210) | **≥30 d** (30 d: 0.231; 7 d: 0.773) | **"~7일 식물체 내 반감기"는 질량 기준 오류를 흡수한 값** |
| **Liu 2023 (중성, ODE) log10 RMSE** | 0.281 | **0.213** | (평형 모드 0.206은 드라이버 무관, 불변) |
| Tang 2026 OOS, 자유음이온 | 1.232 | 1.483 | PFSA 과소가 커짐 |
| Tang 2026 OOS, 지질 메커니즘 | 0.516 | 0.561 | 일반화 결론은 유지 |
| Tang 2026 in-sample refit | 0.519 | 0.844 | refit f_xy가 dry 기준이므로 예상된 결과 |
| Kim 2019 낟알 OOS (PFOA 제외) 지질/단조/W2 | 0.48 / 2.05 / 1.07 | 0.55 / 2.77 / 1.67 | |
| Kim, 검출빈도≥15% 지질 | 0.20 | 0.81 | 지질 상수가 dry 기준에서 적합됨 |
| `refit_oryza` 포화 재현 RMSE (Yamazaki) | 0.236 | **0.283** | PFOA f_xy 0.083→0.359; **PFBA·PFOS·PFDoDA f_xy=1 상한**; 단쇄 낟알은 L_Ph 하한에서도 과대 |
| 중성 ODE/평형 뿌리 비(“할인”), logKow 1.78/2.25 | 0.55 / 0.57 | **0.92 / 0.93** | “ODE-vs-평형 채점 인공물”은 대부분 질량 기준 오류였음(뿌리가 10배 가벼워 xylem 유출을 못 따라감) |
| Kodešová 잎/뿌리 (실측 중앙값 3.25): 대사 없음 / 7 d / 1.5 d | 181 / 86 / 28.5 | **35 / 16 / 4.5** | “대사만으로는 못 닫는다”는 **dry에서만 참**. fresh에선 1.5 d 반감기로 ~1.4배까지 접근(Kodešová는 잎 모화합물 비율 0.49를 **측정**함) |
| carrier two-pool PFDoDA 폐합 (`close_longchain`) | f_xy 낮음 + carrier 증강 | **f_xy가 상한 1** | “낮은 f_xy + 강화 carrier” 해석은 dry 기준 결과 |

결론: **적합하지 않은 중성 경로가 크게 좋아졌다**(수정이 옳다는 독립 증거). PFAS 쪽은
질량 기준을 맞추면 지상부가 4–6배 덜 쌓이므로 같은 관측을 맞추려면 f_xy가 ~4배 커져야 하고,
PFOS·PFDoDA 짚과 단쇄 낟알은 **단일 짚 4-pool 구조로 맞출 수 없게 된다**. 즉 기존 "포화 적합으로 재현"
일부는 버그 덕분이었다. `*_oryza`는 생중량 기준으로 **재적합함**(수정 전 값은 `*_oryza_dry`로 보존).
`TANG_REFIT_FXY`, `LIPID_LOADING`, two-pool 적합은 아직 dry 기준이다(§5). 수정 전 수치는 모두
`mass_basis="dry"`로 재현된다.

### 1.4 성장모델 앵커 검증
| 앵커 | 레포 | 원문 확인 | 판정 |
|---|---|---|---|
| Bouman & van Laar 2006 DOI | `…2004.07.011` | `10.1016/j.agsy.2004.09.011` | **수정함** |
| 지상부 1740 g/m² | "Bouman 2006" | Bouman 2006에 없음. Kumar 2017(PMC5291415) 생장상 390 ppm 4시즌 평균 1744 | **출처 수정함** |
| IR72 분배표·SLA·DRLV·DVR·CV | "IR72" | Bouman 2006 본문엔 표 없음(Bouman 2001 매뉴얼로 미룸). APSIM의 IR72.D92 사본과 다름(개화 후 줄기 분배 25–30% vs 0) | **미해결**(§5) |
| `growth_rice` 분배표 출처 | "Li 2017 SI Table S1" | 본문상 Table S1은 적합 통계; 값은 Table S2(IR72) | **표기 수정함** |
| HI | 0.46 | Bouman 0.51–0.55, Kim FACE 0.495–0.53 | 낮음(미수정) |
| 뿌리:지상부 | "0.1–0.25, Nada 2018 verified" | Nada는 **생중량**, 담수 논 아님, 0.18–0.80. 논 건물중: Kim 2001 FACE PI 0.12–0.14, 성숙기 0.026–0.033(블록 채취라 하한) | **DB 수정함** |
| 개화 66일·성숙 116일 | 이식 후 | 원출처는 파종 후 일수 → 이식 후 약 62/103일 | 미수정 |
| θ_fw 0.90/0.83/0.78/0.14 | 고정 | 출처 키가 레포에서 해석되지 않음; 단계별 변화 없음 | **공백**(§5) |

---

## 2. 데이터셋별 기관 분취·전처리

### 2.1 낟알 — 데이터셋마다 다른 부위
| 데이터 | 실제 분석 부위 | 처리 |
|---|---|---|
| Yamazaki 2023 | 왕겨 / 쌀겨 / 백미 (**현미 미측정**) | 70 °C 전체 건조 → 탈곡기 탈각 → 도정(백미:현미 86–90%). 레포는 0.88·백미+0.12·쌀겨로 재구성. 논문 "현미"(S20)는 농도 단순합이라 쓰면 안 됨 |
| Kim 2019 | 현미 (필지당 1주) | −30 °C 동결건조 → 도정기 탈각 → 분쇄. 백미 비교는 시중 쌀 1점 |
| Tang 2026 | 왕겨(chaff) + endosperm | 도정법 미보고. endosperm은 백미 추정(SI가 백미 수분 기준 인용), chaff에 쌀겨 포함 가능. 5개월째에만 존재 |
| Li 2025 | "paddy rice" | 정의·도정 미보고. 출수기부터 보고 |
| 모델 | 현미 조성 B_k | ORYZA WSO는 왕겨 포함 이삭 전체 |

### 2.2 뿌리 — plaque 제거한 연구 없음
Yamazaki 수돗물→Milli-Q · Tang 탈이온수 헹굼 · Li 2025 증류수 세척 · Kim 미분석 ·
**Liu 2023 탈이온수 3회 + 아세토니트릴 신속 헹굼(유일한 탈착)** · Ge 미보고 · Kodešová 헹굼(무=저장근+잔뿌리) ·
Briggs 1982 blot만 · Hwang 흐르는 물 + 종이타월.

### 2.3 지상부·건조·반복
| 데이터 | 지상부 | 건조 / 기준 | 반복 |
|---|---|---|---|
| Yamazaki | 높이 4구간(잎 포함) | 70 °C 1주 후 분리, dw, 기관 수분 없음 | 혼합 1점 |
| Kim | — | 동결건조, 기준 미기재(dw 추정) | 필지당 1주 |
| Tang | 줄기/잎(잎집 미보고) | −80 °C 동결건조, dw | 3반복 × 20주 |
| Li 2025 | 방법 "잎(줄기 포함)" vs 표 줄기·잎 분리 | 동결건조, dw | 3 × 부분시료 5 × 5단계 |
| Liu 2023 | 가위로 분리(위경) | fw | n=3 |
| Ge 2017 | 줄기/잎 | 본문 Table 2 dw/dw, S6 기준 미정 | 3 × 20주, 비담수 |
| Kodešová | 잎/뿌리 | 동결건조 dw | 5포트 혼합 → 1점 |
| Briggs 1983 | 잎몸/위경 상·하 | fw | 6주 합 |
| Hwang | 잎/뿌리 | fw(방법상) | 3 × 20주 |

전체 근거(쪽·절 위치)는 `Sampling_Protocols.csv`.

---

## 3. 기관 질량 실측 자료
- **위 데이터셋**: 기관별 시계열은 Hwang 2017(S5 Fig)뿐. Tang은 수확기 전체·이삭 생중량(기관 질량은 측정했으나 미공개; endosperm ≈ 식물 건물중의 10–16%로 역산).
- **문헌(추출: `raw_si/rice_organ_biomass_literature.csv`)**:
  Kim 2001/2003 Shizukuishi FACE(Akitakomachi) — **논 조건 뿌리 포함**: 뿌리 PI 45 g/m², 성숙기 38–45 g/m²(Table 1, 원문 대조), 잎몸 시계열(Table 2);
  Bouman 2006 Fig 2b IR72(디지털화, 뿌리 없음); Nada 2018(생중량 R/S, 논 아님).
- **공개 자료 후보(미수록)**: ORYZA2000 매뉴얼 IR72DS92 관측값(무료, 뿌리 없음), Ye 2014 Table 1(담수 뿌리 비중 37.5→8.3%), Yoshida 1981(낟알 수분 ~58→≤20%), Brunetti 2021(완두, 기관 생중량+건물중 5회).
- **공백**: 벼 기관별 **수분 함량 시계열**은 찾지 못함(常丽英 2008은 측정했으나 미보고).

---

## 4. 레포 데이터 사용 방식의 문제와 조치
| 문제 | 조치 |
|---|---|
| 건물중을 생중량으로 사용 | **수정**(§1) |
| Yamazaki를 온실/포트/주입/공극수로 서술 | **정정**: 앱(`ui/expert.py`), README, CLAUDE.md, docs, DB |
| Kim "필지별 짝 없음", "Table 4", Andosol | **정정**: `literature_params.py` 주석, `raw_si/README.md`, C1 DB(csv+xlsx). 값 자체는 유지(L_Ph 적합 보존) |
| Kim PFOA BAF = 비율의 평균, A-2 83% | 기록(대안: 필지 짝 기하평균 27, n=14). 재적합은 §5 |
| Li 2025 TF 정의 불일치(레포 1.0–5.8배 낮음), 5단계 평균, 물 농도 디지털화 | 기록. 단계별 비교로 재평가 필요(§5) |
| Tang "HI_dw≈0.5" 주석, endosperm=현미 가정 | 주석 정정; 조성 변경은 §5 |
| Nada R/S "verified" | **DB 정정** |
| Bouman DOI, 1740 출처, Li 2017 표 번호 | **정정** |

---

## 5. 남은 일 (우선순위)
1. **생중량 기준 재적합 — ✅ 다음 세션에서 완료**(결과·해석: `docs/HANDOFF_mass_basis_refit.md` §1, CLAUDE.md §6
   "Fresh-basis RE-FITS"). 아래는 원래 목록. `*_oryza`는 완료(RMSE 0.283). 남은 것: `TANG_REFIT_FXY`(`validation/tang2026_fxy_refit.py`),
   `LIPID_LOADING`, two-pool 적합(`twopool_fitted_params_measured.json`, `twopool_nstem_fitted_params.json`), 그리고
   이들에 기댄 OOS 결론(Tang 지질 0.516, Kim 0.48 등)의 재확인. 수정 전 결과는 `mass_basis="dry"`로 재현 가능.
   직접 M을 만드는 validation 스크립트 다수(`twopool_*`, `tang2026_fxy_refit.py` 등)는 아직 건물중을 그대로 쓴다 —
   재실행 전에 `organ_mass.to_model_mass`로 바꿀 것.
2. **θ_fw 출처 확보·단계별 변화**: 기관 수분 실측(특히 등숙기 낟알, 성숙기 짚).
3. **Li 2025 단계별 비교**(SI Table S6, 5단계 × 기관)와 저자 정의 TF.
4. **Kim PFOA**: 필지 짝 BAF 사용 여부 결정; 장쇄는 검출 1–4건이므로 censored 처리.
5. **ORYZA 파라미터를 IR72.D92 원본과 대조·정렬**(숫자가 다시 크게 바뀌므로 별도 작업) + FACE 뿌리/IR72DS92로 검증.
6. **Tang endosperm을 백미 조성으로**; `tang2026_fxy_TF_validation.py`의 120일 → 150일. — ✅ 완료(백미는 선택지, 기본 현미; 백미면 곡물 격차가
   커짐). 2번(θ_fw)은 민감도까지 완료: 뿌리 θ만 민감, 그리고 물관 분배의 질량 기준이라는 새 결정 사항이 나옴(`docs/HANDOFF_mass_basis_refit.md` §2).

## 6. 저자 문의가 필요한 것 (원문에도 없음)
Tang(도정법·chaff 정의·잎집 배정·기관 질량·담수 농도·월별 수치) · Li 2025(단계별 날짜·낟알 정의·물질별 관개수) ·
Yamazaki(혼합 주수) · Kim(현미 dw/fw, 공극수 채취 깊이) · Ge(SI S6 기준) · Liu(72 h 용액 농도, S4 단위).

## 7. 출처
Yamazaki 2023 `10.1021/acs.est.2c08767` · Kim 2019 `10.1016/j.scitotenv.2019.03.240` · Tang 2026 `10.1016/j.jhazmat.2025.141017` ·
Li 2025 `10.1016/j.jhazmat.2025.138256` · Liu 2023 `10.1016/j.scitotenv.2022.159826` · Ge 2017 `10.1016/j.envpol.2017.04.043` ·
Kodešová 2019 `10.1007/s11356-019-04333-9` · Briggs 1982 `10.1002/ps.2780130506` · Briggs 1983 `10.1002/ps.2780140506` ·
Hwang 2017 `10.1371/journal.pone.0172254` · Li 2019 `10.1016/j.envint.2019.02.020` · Schriever 2020 `10.1016/j.scitotenv.2020.136667` ·
Bouman & van Laar 2006 `10.1016/j.agsy.2004.09.011` · Kim 2001 `10.1046/j.1469-8137.2001.00111.x` · Kim 2003 `10.1046/j.1365-2486.2003.00641.x` ·
Nada 2018 `10.1007/s11738-018-2697-5` · Kumar 2017 PMC5291415.
