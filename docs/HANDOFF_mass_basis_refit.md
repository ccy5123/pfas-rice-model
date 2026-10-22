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

## 1. 남은 건물중 기준 적합값을 생중량 기준으로 재적합  ← 최우선
현재 `*_oryza`만 재적합됨(RMSE 0.236→0.283; PFBA·PFOS·PFDoDA f_xy=1). 나머지는 수정 전(dry) 값 그대로다.

| 대상 | 위치 | 할 일 |
|---|---|---|
| `TANG_REFIT_FXY` | `src/model_api.py`, `validation/tang2026_fxy_refit.py` | 스크립트는 `simulate_nstem_leaf(biomass_fn=…)`를 쓰므로 **재실행하면 자동으로 fresh**. 새 값을 기록하고, dry 값은 `*_dry` 상수로 보존(`mass_basis="dry"`일 때 사용) |
| `LIPID_LOADING` (a, b, g_xy, g_ph, K_half, pfsa_ln) | `src/model_api.py` | **적합 스크립트가 레포에 없음**(문서 `docs/fxy_longchain_lipid_exploration.md`, growth_rice 측정 forcing에서 Yamazaki excl. PFDoDA). 스크립트를 새로 만들어 fresh 기준 재적합 → 그다음 `validation/oos_tang_lipid.py`, `validation/oos_multidataset.py`로 OOS 결론(Tang 0.516, Kim 0.48 → 현재 fresh에서 0.561 / 0.55) 재확인 |
| two-pool seq (측정 forcing) | `validation/twopool_root_measured.py` → `twopool_fitted_params_measured.json` | 내부에서 `growth_rice` 건물중으로 `M`을 직접 만든다 → `organ_mass.organ_matrix`로 바꾸고 재적합(~3분) |
| two-pool + nstem merge | `validation/twopool_nstem_merge.py` (`MASS_BASIS="dry"` 고정) | `MASS_BASIS="fresh"`로 재적합(~40분) → Tang OOS(dry 0.801) 재확인. 캐시 json은 `_dry`로 보존 |
| carrier two-pool 폐합 | `src/pfas_rice_two_pool.py::close_longchain` | fresh에서 PFDoDA f_xy가 상한 1 → "낮은 f_xy + 강화 carrier" 해석 재검토 |
| `M`을 직접 만드는 validation 스크립트 | `absolute_fxy_fit`, `bayesian_inverse_demo`, `hydrus_coupled_run`, `refit_realistic_biomass`, `root_shoot_biomass_sensitivity`, `twopool_longchain`, `twopool_root_measured`, `oryza_growth_validation`(BAF 부분) | 재실행 전에 `organ_mass.to_model_mass` 적용. (`biomass_shape_compare`, `mass_drivers_plot`은 건물중 비교가 목적이라 그대로) |

완료 기준: 재적합 값이 parameters.json/상수에 반영되고 dry 값은 보존, 관련 테스트가 두 기준을 모두 고정, CLAUDE.md §6 수치 갱신.

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
docs/HANDOFF_mass_basis_refit.md를 읽고 1단계부터 순서대로 진행하자.
먼저 PR #69 상태(CI, 머지 여부)를 확인하고, 1단계(남은 건물중 기준 적합값의 생중량 재적합)를 시작해.
각 재적합마다 수정 전 값은 _dry로 보존하고, 테스트는 두 기준을 모두 고정해.
```
