# 인계: 생중량 기준 후속 작업 (다음 세션용)

> 이전 세션: 실측 데이터 분취·전처리 감사 + **건물중→생중량 기관 질량 버그 수정**
> (PR [ccy5123/pfas-rice-model#69](https://github.com/ccy5123/pfas-rice-model/pull/69), 브랜치
> `claude/sharp-planck-lszv93`). 전체 기록: `docs/sampling_protocols_audit_KR.md`, 요약: CLAUDE.md §6
> "Sampling-protocol audit + the DRY→FRESH organ-mass fix".
> **아래 순서대로** 진행한다. 각 단계는 앞 단계 결과에 기댄다.

## 0. 먼저 확인할 것
- PR #69의 CI(`tests.yml` 전체 스위트, `rigor.yml`)가 초록인지, 머지됐는지. 로컬 결과: 373 collected, **372 pass, 2 skip**. (1단계 이후: 381 collected. 1b·2·6단계 이후: 390 collected, **390 pass, 1 skip**(emcee 설치 시), ~31분.)
  머지 전이면 이 브랜치에서 이어서, 머지 후면 `main`에서 새 브랜치로 시작.
- 핵심 규약: 성장모델 출력은 **건물중**, ODE의 `M`은 **생중량**. 드라이버를 `M`으로 만드는 곳은 반드시
  `src/organ_mass.py`(`to_model_mass` / `organ_matrix`)를 거친다. `mass_basis="dry"`는 수정 전 재현용일 뿐.
- 이 컨테이너의 스크래치 파일·사용자 제공 PDF는 다음 세션에 **없다**. 본문을 다시 읽어야 하면 사용자에게
  `papers.zip`(PFAS 4편 + 성장모델 7편)과 `DPU4OC_add.zip`(중성 논문)을 다시 요청. 공개 SI는 아래 URL로 직접 받을 수 있다.

## 1. 남은 건물중 기준 적합값을 생중량 기준으로 재적합 — ✅ 완료 (다음 세션, PR #69 위에 쌓은 브랜치 `claude/confident-newton-82t0qo`)
모든 재적합은 생중량에서 다시 했고, 건물중 적합은 **지우지 않고 옆에 보존**했다(`mass_basis="dry"` / `PFAS_MASS_BASIS=dry`로 선택).
각 dry 기록은 문서 수치를 그대로 재현하는지 다시 확인했다. 테스트: `tests/test_mass_basis_refits.py`(6), `test_twopool_nstem_merge.py`(+2),
`test_two_pool.py`가 두 기준을 모두 고정. `parameters.json`, `simulate()` 기본값, `reproduce_demo`(0.029)는 변화 없음.

| 대상 | dry 기록 | **fresh 재적합** | 보존 위치 |
|---|---|---|---|
| `LIPID_LOADING` (Yamazaki excl PFDoDA, in-sample) | 0.372 (root 0.488) | **0.177 (root 0.088)** | `LIPID_LOADING_DRY`; 새 스크립트 `validation/fit_lipid_loading.py` |
| 지질 Tang OOS / 자유 음이온 | 0.516 / 1.232 | **0.520** / 1.483 | `oos_tang_lipid.py` |
| 지질 Kim OOS excl PFOA (reliable DF≥15%) | 0.48 (0.20) | **0.50 (0.74)** | `oos_multidataset.py` |
| `TANG_REFIT_FXY` PFOA/PFOS/GenX | 0.097/0.320/0.017 | **0.570/0.999(상한)/0.089** | `TANG_REFIT_FXY_DRY` |
| Tang in-sample refit RMSE | 0.519 | **0.327** | — |
| two-pool seq (측정 forcing) in-sample / Kim | 0.278 / 0.56 | **0.412 / 0.76** | `twopool_fitted_params_measured_dry.json` |
| two-pool + nstem merge in-sample / Tang OOS | 0.301 / 0.801 | **0.515 / 0.602** | `twopool_nstem_fitted_params_dry.json` |
| carrier 2-pool PFDoDA f_xy / carrier× | 0.43 / 3.15× | **1.0(상한) / 1.50×** | `test_two_pool.py` |

**읽는 법(대부분 in-sample)**
- **지질 로딩의 뿌리 trade-off는 상당 부분 건물중 `M` 버그였다.** 생중량에서 단일 풀 지질 모델이 뿌리와 지상부를 함께 맞춘다(PFUnDA 뿌리 14.7 vs 19.5).
  이 trade-off가 two-pool 뿌리 연구의 출발 동기였다. 단, 6개 전역 상수는 congener별 W2 `kappa_d`/`L_Ph` 위에 얹혀 있어 자유도가 6이 아니다.
- **지질 OOS 일반화는 살아남는다**(Tang 0.52 vs 1.48). 그러나 "OOS ≈ in-sample refit(0.516 vs 0.519)"은 더 이상 성립하지 않는다(0.52 vs 0.33).
  Kim의 신뢰 구간(DF≥15%)은 나빠졌다(0.20→0.74).
- **PFOS f_xy는 Tang·Yamazaki 모두 생중량에서 상한 1** → "조건 의존적(0.14 vs 0.32)"이 아니라 "자유 음이온으로는 둘 다 부족"으로 읽힌다.
- **two-pool은 생중량에서 약해진다.** 11개 중 5개 congener는 sink가 필요 없고(k_seq=0에서도 뿌리 과대), PFOS/PFUnDA의 비-K_PL 분리(dry 4.5×)는 **정의 불가**(두 sink 모두 ~0).
  스크립트가 0에 가까운 두 값의 비로 "HELD (30246×)"를 출력하던 판정 논리를 고쳤다.
- **carrier 2-pool의 "낮은 f_xy + 강화 carrier"는 건물중 해석**이었다. 생중량에서 장쇄 carrier 강화는 1.26–1.50×로 거의 평탄. PFDoDA 짚은 f_xy=1에서도 21 vs 49.7(지상부 적재 상한은 남음).
- `TANG_REFIT_FXY_DRY` 자체도 이후 모델 변경 전의 값이다: 현재 코드로 dry 재적합하면 0.114/0.376/0.020.

**변환만 하고 재실행하지 않은 스크립트**(`PFAS_MASS_BASIS=dry`면 수정 전과 동일한 `M`): `absolute_fxy_fit`, `bayesian_inverse_demo`, `hydrus_coupled_run`,
`refit_realistic_biomass`, `root_shoot_biomass_sensitivity`(분배 통계는 건물중 유지, placeholder 행은 기준 없음), `twopool_longchain`과 이를 import하는 `longchain_*`.
문서에 남은 이들의 수치는 dry 기록이다. sci-adk `runs/*` 기록도 dry 기록으로 동결(rigor guard는 검증만 하고 재계산하지 않음).

### 1b. 1단계가 드러낸 후속 작업 (2단계 전에 할지 결정 필요)
- ✅ **carrier-vs-bypass / dose-series 판정을 생중량에서 재확인 — 완료, 판정 유지** (`docs/HANDOFF_carrier_vs_bypass.md` §5).
  두 스크립트에 `PFAS_MASS_BASIS` 추가(`--fast` dry가 기록 2.640/1.035/0.996/2.289, bootstrap 0.749, corr +0.832를 정확히 재현).
  dose-series는 기준과 무관(PFOS 6.27× vs 1.17×, `Km`≥500, 1/3). "추가 항 필요"는 유지(아무것도 없음 3.004, 탈분극 2.634).
  **"어느 쪽인가"는 공정성 수정 후에만 유지**: 그대로 돌리면 bypass만 재적합되고 carrier `Vmax`는 W2 값이라 P=0.996으로 bypass가 이기는 착시.
  새 사후 섹션 6이 carrier 전역 매개변수도 재적합 → fresh 1.161 vs 1.194, P=0.30(dry 0.16). bypass의 사슬 길이 무관 주장은 더 강하게 반박(+0.909).
  **수렴은 같은 설정에서 처음 확인**: carrier `Vmax`도 congener별로 corr +0.915(fresh)/+0.880(dry). f_xy가 자유인 carrier two-pool에서는
  보정이 사라짐(1.26–1.50×) → 보정은 자유로운 항으로 옮겨 간다(B3 "η 잘못 분할" 해석과 일치, `theory_anchor.tex`에 기록).
- ✅ **곡물 "구조적 ~3–8× 과소" 재검토 — 6단계와 함께 완료** (`docs/tang2026_grain_units_exploration.md` §6): 생중량에서 모델 의존적
  (지질 PFOA 2.8×/PFOS 5.8× 과소, refit PFOA 맞음·PFOS ~4× 과소). 백미 해석은 격차를 **키운다**(아래 6단계).
- **two-pool 위상 결정**: 단일 풀 지질(0.177)이 두 two-pool 변형(0.412 / 0.515)보다 낫다. 승격 게이트(습식 실험)는 그대로지만, 탐색 우선순위를 낮출지 사용자와 결정.
- 변환만 한 스크립트의 재실행(특히 `hydrus_coupled_run`: "constant Cwo가 곡물/짚 BAF를 2–4× 과대" 결론).

## 2. θ_fw(기관 수분 함량)의 근거 확보와 시기별 변화
- 현재 root 0.90 / stem 0.83 / leaf 0.78 / grain 0.14 고정(`params/parameters.json` → `tissue_composition_recommended`),
  출처 키(`S_physiol`, `S_grainMC`)가 레포에서 해석되지 않음. B_k와 질량 변환이 **같은 θ**를 쓴다.
- 후보 자료: Yoshida 1981(낟알 수분 ~58%→≤20%; `http://books.irri.org/9711040522_content.pdf`), Brunetti 2021 SI
  Fig S1/S2(완두 기관 생중량+건물중 5회; DPU4OC zip), 常丽英 2008은 측정만 하고 미보고.
- 설계 결정 필요: θ를 시간 함수로 만들면 B_k도 시간에 따라 변한다(ODE 구조 변경). 먼저 영향 크기(민감도)를 보고 결정.
- ✅ **민감도 완료** (`validation/theta_fw_sensitivity.py`, 사전등록 커밋 `4f74a02`를 결과보다 먼저 push). 한 기관씩 범위 끝값, 건물중 쪽 지표(Tang TF, Ge 2017 TF):
  뿌리(0.88–0.92) Tang 0.039 / Ge **0.053 → INTERMEDIATE**, 줄기 0.000/0.015, 잎 0.001/0.002, 곡물(0.14–0.60) 0.002/0.008 → SMALL.
  지상부 θ는 건물중 기준에서 상쇄된다(생중량 BAF는 줄기 θ로 0.13 log 움직임 → **결과는 건물중으로 인용**). 뿌리 θ는 시기별 변화보다 **값 자체의 불확실성**이
  문제다(담수 논에서 뿌리는 내내 수화) → θ(t) ODE 변경은 불필요, **뿌리 수분 실측**이 필요한 입력. ⇒ 판정 규칙상 뿌리는 사용자 결정 사항.
- ⚠️ **사후 발견 — 물관 분배의 질량 기준(사용자 결정 필요)**: 4-pool 코어는 물관을 잎/곡물에 `S·M`으로 나누는데(`4pool_surf.py:410`), PR #69 이후
  `M`이 생중량이 되어 **수분 함량이 증산 몫을 정한다**. `S`(잎 20, 곡물 2)는 비율로만 정의돼 기준이 없다. 수확기 곡물 몫: 수정 전(건물중) 0.343 →
  현재 0.118. 건물중 분배로 되돌리면 PFOA 곡물 BAF 0.036→0.103(2.8×), 지질 적합 0.177→0.271(재적합 전). PR #69의 곡물 0.148→0.036(4.1×) 중
  질량 환산(×1.16)을 넘는 부분 대부분이 이 몫 변화로 보인다. 결정하면 1단계 재적합(지질·oryza)을 다시 해야 한다.
  `tests/test_theta_fw_sensitivity.py::test_open_issue_the_xylem_split_runs_on_fresh_grain_mass`가 현 상태를 고정(결정 시 이 테스트를 바꿀 것).
- **결정(사용자, 다음 세션): A — 건물중 분배.** 근거는 적합도가 아니라 **정의**다: `S`를 작물 모델의 SLA처럼 건물중 1 kg당
  면적으로 읽는다(잎 면적 = SLA × 잎 건물중; ORYZA의 SLA도 m²/kg 건물). 구현: `PlantInputs.M_basis`("fresh"면 코어가 자기 θ로
  `M·(1−θ)`에서 분배, None/"dry"면 M 그대로) + `xylem_leaf_share()`. 건물중 기록·placeholder(`reproduce_demo` 0.029)·two-pool demo는
  비트 단위 동일. 재적합 전 효과: PFOA 곡물 BAF 0.036→0.103, PFOS(growth_rice) 0.0084→0.0245, 중성 carbamazepine 곡물 3.87→12.6.
  **사전등록 (재적합을 돌리기 전에 기록)**:
  1. `refit_oryza`(재현 RMSE 0.283): 단쇄 곡물은 생중량 분배에서도 이미 L_Ph 하한에서 과대였다. 곡물로 가는 물관 몫이 ~3배가 되면
     단쇄 곡물 과대가 커지고 L_Ph로는 더 낮출 수 없으므로 **RMSE가 나빠질 것**으로 예상(개선되면 이 예상이 틀린 것).
  2. `LIPID_LOADING`(0.177; 재적합 없이 0.271): 재적합이 대부분 회복하리라 예상하되 0.177보다는 나쁠 것. **0.25 초과로 끝나면** 건물중
     분배가 Yamazaki 곡물과 실제로 어긋난다는 뜻으로 기록.
  3. 재적합에 쓰지 않은 독립 점검: Tang 기관별 TF(lipid)는 지상부가 `nstem_leaf`(고정 분배 `lam_grain`)라 **지질 상수를 통해서만**
     바뀐다 → 0.520에서 ±0.05 안으로 예상. Kim 곡물(4-pool)은 방향을 미리 말할 수 없다(곡물이 오르는데 단조 모델은 과소, 지질 모델은
     장쇄 과대 위험).
  4. **판정 규칙**: 분배 기준은 이 숫자들로 정하지 않는다(정의 문제). 단, 독립 점검 둘(Tang, Kim의 lipid OOS)이 **모두 0.10 넘게
     나빠지면** 건물중 분배에 불리한 증거로 기록하고 사용자에게 보고한다.
- **결과 — 1번 예상 맞음, 정도는 예상보다 컸다(재적합 2개만 돌리고 멈춤)**: `refit_oryza` 0.283→**0.394**, 11개 중 **8개의 L_Ph가
  하한(1e-5)**에 붙음 — 체관을 사실상 꺼도 물관만으로 곡물이 과대(PFBA 11.65 vs 관측 1.06, PFPeA 3.31 vs 0.26, PFBS 2.52 vs 0.25).
  two-pool measured 0.412→**0.513**. 원인은 곡물의 **증산 몫**: 등숙기 평균으로 생중량 분배 4.2%, 건물중(S 20:2) 13.7%(growth_rice
  4.4% / 14.8%); 레포의 `nstem_leaf`는 5%(`lam_grain`)를 가정한다. 벼 이삭 증산 비율은 빠른 문헌 검색으로 찾지 못했다.
- **핵심 — 두 기준은 같은 식이다**: θ가 상수인 한 "생중량 분배(S·M_fresh)"는 "건물중 분배 + `S_grain` = 2×(1−θ_잎)/(1−θ_곡물) = 0.5116"과
  **정확히 같다**(확인: PFOA BAF 6자리 일치). 즉 정할 것은 기준이 아니라 **곡물 증산 몫(`S_grain`의 값)**이었다.
- **최종 결정(사용자): ① 건물중 정의 + `S_grain` 0.5116/kg 건물** (`S_XYLEM_DRY = {잎 20, 곡물 0.5116}`, 코어 상수). θ에서 계산하지 않는
  **고정 상수**라서 분배가 θ에 좌우되지 않는다(곡물 θ 0.14→0.60에서 수확기 몫 0.118 그대로). 출하 θ에서는 기존 생중량 분배와 상대오차
  5e-5 이내로 같아 **PR #70의 모든 수치가 그대로**이고 재적합이 필요 없다. 위 20:2 재적합 결과는 증거로만 남기고 되돌렸다.
  `Compartment.S`(20/2)는 "M 그대로" 경로(placeholder·건물중 기록, 비트 단위 동일)와 `plant_air` 표면적에 계속 쓰인다.
  **약점**: 0.5116은 측정값이 아니라 기존 적합·5% 가정과 맞는 값이다. 이삭/낟알 증산 몫(또는 낟알 기부 물관 불연속) 측정이 이 값을 대체할 입력.
  θ 민감도 재실행: 판정 동일(뿌리 INTERMEDIATE, 나머지 SMALL), 곡물 θ의 Ge 영향 0.008→0.000.
  테스트: `test_theta_fw_sensitivity.py`(θ 무관 + 기각된 해석의 증산 몫 14%), `test_organ_mass.py::test_the_xylem_split_follows_the_mass_basis_of_M`.
- θ의 세 사본(`model_api._COMP`, `organ_mass.THETA_FW`, `neutral_dpu.RICE_WATER`) 일치를 테스트로 고정(`test_organ_mass.py`).

## 3. Li 2025를 생육단계별로 다시 비교 — ✅ 완료 (사전등록 `4df3047`, 결과보다 먼저 push)
- **전사**: SI Table S6/S5 → `raw_si/li2025_stage_conc.csv`(15종 × 5단계 × 기관, dw), `li2025_stage_water_soil.csv`(물질군 합만).
  5단계 평균이 `data_obs/Li2025_BAF_TF.csv` 20칸을 모두 재현, 단계별 비율 평균이 저자 PFOA 줄기 TF 2.35를 재현.
- **설계**: TF = 조직/뿌리, **건물중**, 단계별(물 농도 불필요). 단계 날짜가 없어 ORYZA DVS로 대응(분얼 0.4·출수 1.0·등숙 1.25·
  유숙 1.5·성숙 2.0 → 28·66·78·90·116일). 양쪽에 **같은 추정식(단계 기하평균)·같은 단계**. 곡물은 등숙기부터 채점.
- **결과** (`validation/li2025_stage_resolved.py`, ~25초):
  - Q1 추정식만 바꿔도 관측 TF가 중앙 3.0배·최대 5.8배 달라진다(저자 "비율의 평균"은 뿌리 농도가 작은 단계 하나에 끌려 올라감,
    PFOA 줄기 2.35 중 등숙기 비율 9.67).
  - Q2 단계 간 산포가 보고된 SD로 설명되는 양의 **4.5배** → 단계별 모양 비교는 잡음 한계. Q4 모양 신호 없음(모든 변형 ρ≈0).
  - Q3 **판정: `nstem_leaf` + lipid가 구분된다** — 지상부/곡물 RMSE 0.43/0.52 vs 4-pool lipid 0.81/1.52, 단조 1.25/2.42,
    W2 1.31/2.70. DVS ±0.15에서도 동일. **사전등록한 예상(결론 못 냄)은 틀렸다.**
  - 사후(판정에 미사용): 다른 두 추정식에서도 같은 승자(0.66/0.66, 0.36/0.47). 분해 — lipid 없는 `nstem_leaf` 1.11/1.20,
    4-pool + lipid 0.81/1.52 → **지상부 구조와 lipid 항이 둘 다 필요**(Tang에서 이긴 조합과 같다). 잡음 하한 0.31/0.32.
    TF는 노출 수준에 무관(Cwo 1 vs 0.05에서 ≤0.001 log).
  - 한계: 노지 단일 필지, 뿌리 표면 흡착 미제거(관측 TF를 낮추는 쪽 → 모델 과소 예측을 설명하지 못함), 단계 날짜 가정,
    ORYZA IR72 vs Jinyu Japonica No.22. `nstem_leaf`의 곡물은 여전히 과소(PFBS 0.30 vs 2.15, PFOA 0.33 vs 0.93).
  - 기존 `oos_crossdataset.li_tf`(수확일 모델 vs 5단계 평균, 생중량)는 기록으로 남기고 docstring에 이 비교를 가리킴.
    `runs/pfas-rice-oos-multidataset`의 "Li inconclusive"는 동결된 기록.
- (원래 계획)
- 레포의 `data_obs/Li2025_BAF_TF.csv`는 5단계 평균 농도의 비(저자 TF보다 1.0–5.8배 낮음). 저자 TF = 단계별 비율의 평균.
- SI Table S6(단계별 농도, dw): `https://ars.els-cdn.com/content/image/1-s2.0-S0304389425011719-mmc1.docx`
  → `docs/literature_db/raw_si/li2025_stage_conc.csv`로 전사 → 모델(동적)을 단계별(5단계 × 기관)로 비교.
- 주의: 물은 **관개수**(공극수 아님), 물질별 농도는 Fig 2a에서 디지털화한 값. 방법 절은 "leaves (including stems)".
- 영향받는 곳: `validation/oos_crossdataset.py`, `validation/S6_Gap4.py`, two-pool Li OOS.

## 4. Kim 2019 PFOA BAF와 장쇄 처리 — ✅ 완료 (4a 사전등록 `bbb3a8a`, 4b 사전등록 `bedadbf`)
- **4a 장쇄(미검출) 처리** — SI(Table S5 MDL·S6 분포)를 출판사 CDN에서 받아 `raw_si/kim2019_si_detection.csv`로 전사
  (평균은 기존 `kim2019_field_conc.csv`와 동일). 미검출을 0~MDL 어디에 두느냐로 BAF 범위:
  PFHpA 0.39–2.11 · PFOA 4.43–4.92 · PFNA 2.21–7.07 · PFDA 1.55–11.3 · **PFUnDA 19.8–168 · PFDoDA 5.5–630**(114배).
  `literature_params.kim2019_grain_baf_bounds()`. `validation/kim2019_censored.py`:
  - Q1 C11 상승은 미검출 처리와 무관하게 성립, C12 수준은 판단 불가.
  - Q2 기존(미검출=0)/구간 RMSE: 단조 2.77/2.51 · W2 1.67/1.52 · **4-pool lipid 0.50/0.38 · nstem_leaf+lipid 0.60/0.14**.
    판정: **단일 승자 없음** — lipid 계열이 lipid 없는 모델보다 두 점수 모두 1 log 이상 좋지만, 두 lipid 모델의 순위는
    미검출 처리에 달렸다(기존 점수 차 정확히 0.10). 예상대로.
  - `oos_crossdataset`의 "신뢰(DF≥15%)" 점수는 **PFNA 한 물질**뿐(라벨이 PFHpA+PFNA로 잘못돼 있었음 → 수정).
- **4b PFOA 기준점** ✅ — 사용자가 본문 Table 1(p.718)을 전사해 줌 → `raw_si/kim2019_table1_sites.csv`(30필지, SI S6 평균·검출 수 재현).
  `validation/kim2019_pfoa_anchor.py` + `tests/test_kim2019_pfoa_anchor.py`(4):
  - 후보: 평균의 비 4.44(벼 30필지/물 27필지) · 같은 27필지 4.18 · A-2 제외 21.5 · 둘 다 검출된 14필지 기하평균 27.0.
    사전등록한 허용 기준(분자·분모 같은 필지 / 검출 여부·값으로 필지를 빼지 않음 / 한 필지가 절반 넘게 좌우하지 않음)을
    통과하는 것은 **검열 필지별 기하평균 9.94 L/kg, 90% 붓스트랩 [4.57, 20.16]**뿐(26필지: 14 정확값 + 12 벼 미검출을 상한으로).
    평균의 비는 A-2(공극수 PFOA의 82.8%, 비율은 최저 0.98)가 정하고, 27은 비율이 낮은 미검출 12필지를 빼서 위로 치우친다.
  - 벼 ∝ 물^0.52 [0.24, 0.87](일정한 BAF면 1). 모델은 Kim 범위에서 선형(기울기 0.964) → 운반체 포화가 아니라
    수확기 2일 고인 물이 시즌 노출의 잡음 섞인 대리값이라는 뜻(필지 간 sd 0.82 log).
  - 모든 변형이 필지별 구간 아래: nstem_leaf+lipid 1.44(−0.84 log) · lipid 0.34 · 단조 0.036 · W2 0.023.
  - 결정: D1 `kim2019_grain_baf()`는 모든 물질에 평균의 비 유지(다른 다섯 물질은 그것밖에 없음), 필지별 값은
    `KIM2019_PFOA_SITE_BAF`로 옆에. D2 데모 L_Ph 목표 4.43→9.94 — **3% 차이로 결정**(붓스트랩 하한 4.57; 시드 10개 모두 같으나,
    사후 프로파일 우도 구간 [4.18, 19.78]은 4.43을 포함 → 구간 방법에 달림). D3 "PFOA* 적합에 사용" 표시 제거(Kim에 적합된
    배포 파라미터는 없음; L_Ph·LIPID_LOADING은 Yamazaki 적합), PFOA는 별도 OOS 점검으로. D4 A-2가 Σ17 공극수의 67.3% →
    Kim 계열의 절대 수준은 A-2 가중(모양 비교는 유효). 사전등록 예상은 모두 맞음.

## 5. ORYZA를 IR72.D92 원본에 맞추기 — ✅ 이식·검증 완료 (사전등록 `288de63`), 기본값 교체는 사용자 결정 대기
- 범위(사용자와 결정, "누가 물어도 문제가 없는 방법"): 파라미터만이 아니라 **ORYZA2000 잠재생산 경로 전체를 매뉴얼의 Fortran에서 줄 단위로 이식**
  (`src/oryza2000.py`). 잠재생산 모드에서 꺼지는 수분·질소 스트레스 모듈만 제외(계수는 1로 두어 식은 원문 그대로).
- 출처: Bouman et al. 2001 매뉴얼(IRRI 무료 PDF) §7.3 IR72.D92 → `raw_si/oryza_ir72_d92_params.csv`, §7.2 IR72DS92 실측 →
  `raw_si/oryza_ir72ds92_obs.csv`. 매뉴얼 원문·APSIM 사본과 스크립트로 전부 대조.
- 애매한 원문 한 줄: GPPARGET의 AMAX — 두 분기가 이어지려면 REDFT·CO2가 첫 항 전체에 곱해져야 함 → 그 형태가 기본(`amax_form`).
- 원본 엔진·기상은 못 구함: WUR 소스(`O2000V2N13gnufortran.zip`)·`ORYZA2000v2n14s1inclWeather.zip`은 로그인 후에도 403, IRRI v3는 폼 대기,
  옛 드라이브 폴더는 권한 필요. → 기상은 NASA POWER를 매뉴얼에 실린 관측소 12일로 보정(Tmin +1.17, Tmax +1.67 °C, 일사 ×1.211),
  온도만 보정·무보정도 함께 보고. **PHIL1.992를 구하면 그 결과가 대체한다(규칙 동일).**
- 성격: IR72DS92는 매뉴얼의 보정 예제 → 이식본에는 **충실도 시험**이지 ORYZA2000의 독립 검증이 아님.
- 결과(`validation/oryza2000_ir72ds92.py`): 이식본 **ADEQUATE(세 기상 처리 모두)** — PI/개화/성숙 56/82/115 vs 실측 58/83/114,
  지상부 log10 RMSE 0.024(~6%). 기존 `oryza_growth`는 성숙 12일 늦음(126), 초반 생장 느림, 개화 후 줄기 계속 증가, 이삭 덜 채움 →
  이식본이 4/4 계열 + 성숙일에서 우세(R2). 탄소수지 ~1e-15로 닫힘.
- 연결: `simulate(biomass="oryza2000")` 옵트인. 앱 라디오에 추가하고 기존 기본을 "ORYZA Level-1, reduced"로 정직하게 재표기.
  `_biomass_fn`은 모르는 이름이면 이제 오류(전에는 조용히 growth_rice로 대체).
- 사후(규칙 아님): PFAS 기후에서 이식본 성숙 t=93일. 120일 시즌이면 짚 ×1.6–1.9·낟알 ×1.9–2.5, 성숙일 수확이면 짚 ×1.1–1.2·낟알 ×0.6–0.75
  → 낟알 수치는 작물 모델보다 **수확 시점(시즌 길이)**이 더 좌우. 120일 기본 시즌은 늦게 성숙하는 축약 모델에 묶여 있던 값.
- **남은 결정(사용자)**: 기본을 `oryza2000`으로 바꿀지. 바꾸면 ORYZA 기반 적합(`refit_oryza`, `TANG_REFIT_FXY`, 기본값으로 도는 OOS 보고)을
  다시 하고, 시즌을 작물 성숙에 묶는 방식(예: `season=None` → 성숙일)을 정해야 함.
- 아직 안 한 것: 뿌리·잎몸 실측(Kim 2001/2003 FACE, Bouman Fig 2b)과의 비교(다른 품종·기후라 독립 자료), 1740 g/m² 기준점 문제(녹색 vs 고사엽 포함)는
  이식본에서는 기준점을 쓰지 않으므로 해소, 기존 모드에는 그대로 남음.

## 6. Tang 2026 세부 — ✅ 대부분 완료
- ✅ endosperm 백미 조성: `tang_tf_validation(endosperm="brown"|"white")`(`TANG_ENDOSPERM`), 기본 brown(도정법 미보고). 출처 키
  (`S_proximate_rg`·`S_USDA_white`·`S_PL_review`·`S_FAO`)는 레포에서 해석되지 않음(θ 키와 같은 문제). **두 행의 인지질 정의가 다르다**: brown 0.003은
  쌀겨·배아 막 인지질만, white 0.003은 전분 결합 lyso-PL(brown 행이 세지 않는 풀) → brown 정의에 맞춰 white 막 PL = 0.2×0.003/0.88 ≈ 0.0007.
  결과: 백미 계수 PFOA 0.33/PFOS 0.25/GenX 0.66, Tang 점수 모두 악화(지질 0.520→0.700, refit 0.327→0.513). 곡물 격차는 조성으로 안 닫히고
  백미면 더 커진다 → 백미인데 TF≥1이면 쌀겨·배유가 같은 자유 농도를 공유하지 않는다(체관 용질의 배유 축적) 쪽. `tests/test_tang_endosperm.py`.
- chaff TF 사용: 보류(기관 질량 필요 → 저자 문의 목록 그대로).
- ✅ `tang2026_fxy_TF_validation.py` `SEASON` 120→150: 4-pool Tang TF RMSE 1.90→1.74.
- ✅ 뿌리 헹굼만(표면 흡착 → 모든 Tang TF를 낮게 치우칠 수 있음) 주의를 `tang_tf_validation` docstring에 기록.

## 저자 문의 목록 (원문에도 없음)
Tang(Baoliang Chen, blchen@zju.edu.cn: 도정법·chaff 정의·잎집 배정·기관 질량·담수 농도·월별 수치) ·
Li 2025(Fen Jin, jinfen@caas.cn: 단계별 날짜·낟알 정의·물질별 관개수) · Yamazaki(혼합 주수) · Kim(현미 dw/fw, 공극수 채취 깊이, PFOA 외 물질의 필지별 값 — D4를 정량화할 유일한 길) ·
Ge 2017(SI S6 기준) · Liu 2023(72 h 용액 농도, S4 단위).

## 다음 세션 시작 프롬프트 (복사용)
```
docs/HANDOFF_mass_basis_refit.md를 읽자. 1단계, 1b-①·②, 2단계(민감도 + 물관 분배 결정), 6단계는 끝났다.
물관 분배는 건물중 정의 + S_grain 0.5116/kg 건물로 결정됐다(수치 불변, 재적합 불필요).
3단계(Li 2025 단계별)도 끝났다: nstem_leaf + lipid가 구분된다(0.43/0.52).
4a(Kim 장쇄 미검출 처리)도 끝났다: lipid 계열 우위는 견고, 두 lipid 모델 사이는 미정.
4b(Kim PFOA 기준점)도 끝났다: 필지별 검열 기하평균 9.94(평균의 비 4.44는 A-2 한 필지 값), 필지별 비율은 상수가 아님(기울기 0.52).
5단계(ORYZA)도 끝났다: 매뉴얼에서 ORYZA2000을 이식했고 IR72DS92를 재현한다(옵트인 biomass="oryza2000").
먼저 PR #69와 #70의 상태를 확인하고, 사용자에게 기본 바이오매스를 oryza2000으로 바꿀지(바꾸면 재적합 + 시즌을 성숙에 묶기) 묻는다.
PHIL1.992나 ORYZA 원본 소스가 오면 validation/oryza2000_ir72ds92.py를 관측소 기상으로 다시 돌린다.
```
