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
- θ의 세 사본(`model_api._COMP`, `organ_mass.THETA_FW`, `neutral_dpu.RICE_WATER`) 일치를 테스트로 고정(`test_organ_mass.py`).

## 3. Li 2025를 생육단계별로 다시 비교
- 레포의 `data_obs/Li2025_BAF_TF.csv`는 5단계 평균 농도의 비(저자 TF보다 1.0–5.8배 낮음). 저자 TF = 단계별 비율의 평균.
- SI Table S6(단계별 농도, dw): `https://ars.els-cdn.com/content/image/1-s2.0-S0304389425011719-mmc1.docx`
  → `docs/literature_db/raw_si/li2025_stage_conc.csv`로 전사 → 모델(동적)을 단계별(5단계 × 기관)로 비교.
- 주의: 물은 **관개수**(공극수 아님), 물질별 농도는 Fig 2a에서 디지털화한 값. 방법 절은 "leaves (including stems)".
- 영향받는 곳: `validation/oos_crossdataset.py`, `validation/S6_Gap4.py`, two-pool Li OOS.

## 4. Kim 2019 PFOA BAF와 장쇄 처리
- 현재 `kim2019_grain_baf` = 평균의 비(PFOA 4.43). 필지 A-2(하수처리장)가 공극수 평균의 83% → 필지 짝 기하평균 27(n=14, 1–597).
- 결정: L_Ph 적합 앵커로 무엇을 쓸지(평균의 비 / 필지 짝 중앙값 / A-2 제외). 장쇄(PFUnDA·PFDoDA)는 검출 1–4건 →
  censored likelihood 또는 OOS 점수에서 제외 표시.
- 본문 Table 1(필지별 PFOA)은 `papers.zip`의 Kim 본문 p.718(PDF 5쪽)에 있음(전사 필요).

## 5. ORYZA 파라미터를 IR72 원본과 정렬 + 실측 검증
- `src/oryza_growth.py`의 분배표·SLA·DRLV·DVRR·CV가 IR72.D92와 다름(개화 후 줄기 분배 25–30% vs 0; 성숙 116일 vs 원본값으로 97일).
  원본: APSIM 사본 `https://raw.githubusercontent.com/APSIMInitiative/APSIMClassic/master/Model/Oryza.xml`.
- 검증 자료: `docs/literature_db/raw_si/rice_organ_biomass_literature.csv`(Kim 2001/2003 FACE 뿌리·잎몸, Bouman Fig 2b),
  ORYZA2000 매뉴얼 IR72DS92 관측값(`http://books.irri.org/9712201716_content.pdf`).
- 주의: 파라미터를 바꾸면 바이오매스가 바뀌어 **1단계 재적합을 다시 해야** 한다 → 1단계 이후, 그러나 결과를 확정하기 전에 범위를 결정.
- 함께 볼 것: HI 0.46(문헌 0.50–0.55), 개화/성숙일의 파종·이식 기준 혼동, 1740 g/m²(Kumar 2017) 스케일링이 고사엽 미포함 녹색 지상부를 기준으로 하는 문제.

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
Li 2025(Fen Jin, jinfen@caas.cn: 단계별 날짜·낟알 정의·물질별 관개수) · Yamazaki(혼합 주수) · Kim(현미 dw/fw, 공극수 채취 깊이) ·
Ge 2017(SI S6 기준) · Liu 2023(72 h 용액 농도, S4 단위).

## 다음 세션 시작 프롬프트 (복사용)
```
docs/HANDOFF_mass_basis_refit.md를 읽자. 1단계, 1b-①·②, 2단계 민감도, 6단계는 끝났다.
먼저 PR #69와 #70의 상태(CI, 머지 여부)를 확인하고, 2단계의 "물관 분배 질량 기준(S)" 결정과
뿌리 θ 판단을 사용자에게 받은 뒤 진행해. 그 결정에 따라 1단계 재적합을 다시 하거나 3단계(Li 2025)로 간다.
```
