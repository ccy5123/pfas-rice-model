# 인계: 생중량 기준 후속 작업 (다음 세션용)

> 이전 세션: 실측 데이터 분취·전처리 감사 + **건물중→생중량 기관 질량 버그 수정**
> (PR [ccy5123/pfas-rice-model#69](https://github.com/ccy5123/pfas-rice-model/pull/69), 브랜치
> `claude/sharp-planck-lszv93`). 전체 기록: `docs/sampling_protocols_audit_KR.md`, 요약: CLAUDE.md §6
> "Sampling-protocol audit + the DRY→FRESH organ-mass fix".
> **아래 순서대로** 진행한다. 각 단계는 앞 단계 결과에 기댄다.

## 0. 먼저 확인할 것
- PR #69의 CI(`tests.yml` 전체 스위트, `rigor.yml`)가 초록인지, 머지됐는지. 로컬 결과: 373 collected, **372 pass, 2 skip**.
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
- **carrier-vs-bypass / dose-series 판정을 생중량에서 재확인**: 두 스크립트는 기본 `simulate()`(PR #69 이후 fresh)를 쓴다. 테스트는 통과하지만
  핵심 수치(1.035 vs 0.996, bootstrap 0.749, `Km`≥500)가 fresh에서 어떻게 움직였는지 문서와 대조하지 않았다. carrier 강화가 평탄해진 결과(위)는
  "두 진입 항이 같은 사슬 길이 보정을 요구한다"는 수렴 결론을 약화시킬 수 있다.
- **곡물 "구조적 ~3–8× 과소"**(`docs/tang2026_grain_units_exploration.md`) 재검토: 생중량 refit에서 PFOA endosperm은 맞고(1.61 vs 1.37) PFOS는 ~4× 과소.
- **two-pool 위상 결정**: 단일 풀 지질(0.177)이 두 two-pool 변형(0.412 / 0.515)보다 낫다. 승격 게이트(습식 실험)는 그대로지만, 탐색 우선순위를 낮출지 사용자와 결정.
- 변환만 한 스크립트의 재실행(특히 `hydrus_coupled_run`: "constant Cwo가 곡물/짚 BAF를 2–4× 과대" 결론).

## 2. θ_fw(기관 수분 함량)의 근거 확보와 시기별 변화
- 현재 root 0.90 / stem 0.83 / leaf 0.78 / grain 0.14 고정(`params/parameters.json` → `tissue_composition_recommended`),
  출처 키(`S_physiol`, `S_grainMC`)가 레포에서 해석되지 않음. B_k와 질량 변환이 **같은 θ**를 쓴다.
- 후보 자료: Yoshida 1981(낟알 수분 ~58%→≤20%; `http://books.irri.org/9711040522_content.pdf`), Brunetti 2021 SI
  Fig S1/S2(완두 기관 생중량+건물중 5회; DPU4OC zip), 常丽英 2008은 측정만 하고 미보고.
- 설계 결정 필요: θ를 시간 함수로 만들면 B_k도 시간에 따라 변한다(ODE 구조 변경). 먼저 영향 크기(민감도)를 보고 결정.

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

## 6. Tang 2026 세부
- endosperm을 백미 조성으로(현재 `grain_brown`). `params/rice_tissue_params.csv`에 백미(white/polished) 행이 있으니 그것을 출처와 함께 확인해 사용.
- chaff TF(PFOA 1.10 / PFOS 0.65 / GenX 1.18) 사용 여부(왕겨+endosperm 질량가중은 기관 질량이 필요 → 저자 문의).
- `validation/tang2026_fxy_TF_validation.py`의 `SEASON=120` → 150(논문·`model_api`와 일치).
- 뿌리는 헹굼만 → Li와 같은 표면 흡착 주의 표시.

## 저자 문의 목록 (원문에도 없음)
Tang(Baoliang Chen, blchen@zju.edu.cn: 도정법·chaff 정의·잎집 배정·기관 질량·담수 농도·월별 수치) ·
Li 2025(Fen Jin, jinfen@caas.cn: 단계별 날짜·낟알 정의·물질별 관개수) · Yamazaki(혼합 주수) · Kim(현미 dw/fw, 공극수 채취 깊이) ·
Ge 2017(SI S6 기준) · Liu 2023(72 h 용액 농도, S4 단위).

## 다음 세션 시작 프롬프트 (복사용)
```
docs/HANDOFF_mass_basis_refit.md를 읽자. 1단계(생중량 재적합)는 끝났다.
먼저 PR #69와 그 위에 쌓은 PR의 상태(CI, 머지 여부)를 확인하고,
1b(1단계가 드러낸 후속 작업) 중 무엇을 2단계(θ_fw) 전에 할지 정한 뒤 진행해.
```
