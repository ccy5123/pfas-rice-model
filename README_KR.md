# PFAS–벼 4구획 흡수 모델 — 재현 패키지

[English](README.md) | **한국어** &nbsp;·&nbsp; 앱 매뉴얼: [English](docs/MANUAL_EN.md) | [한국어](docs/MANUAL_KR.md)

논(paddy)에서 자라는 벼(*Oryza sativa*)의 **영구 음이온 PFAS** 흡수를 다루는 메커니즘 모델입니다.
DPU / Trapp 프레임워크의 **이온성유기화합물(IOC) 확장**으로, 4구획(뿌리 → 줄기 → 잎 → 낟알)을 다룹니다.
대상은 **보정된 12종 congener**: PFCA C4–C12(PFBA…PFDoDA), PFSA C4/C6/C8(PFBS, PFHxS, PFOS),
**그리고 GenX**(HFPO-DA, 에터 PFAS — provisional, Tang 2026 검증용으로 추가; `docs/VALIDATION_TANG2026_KR.md`).

이 패키지는 종료된 두 개방 파라미터 작업 — **GAP A(세포벽 분배 `K_cw`)**, **GAP B(뿌리→지상부 적재 `f_xy`)** —
와 모델 코드, basis-A 결합인자, S6 검증을 하나의 재현 가능한 묶음으로 통합합니다.

> 🖥️ **앱 사용법**(일반인~전문가)은 별도 매뉴얼 **[`docs/MANUAL_KR.md`](docs/MANUAL_KR.md)** (영어: [`docs/MANUAL_EN.md`](docs/MANUAL_EN.md)) 참조.

---

## 빠른 시작

```bash
pip install -r requirements.txt           # numpy, scipy, matplotlib
python build_parameters.py                # 소스 표에서 params/parameters.json 재조립
python reproduce_demo.py                  # Yamazaki BAF 전체 ODE 재현 (W2 fit; log10 RMSE ≈ 0.029)
python reproduce_demo.py --rec            # monotone 물리적 f_xy (단일-짚 불일치 — 주의 참조)
python src/literature_params.py           # 문헌 QSPR (K_PL/K_prot/Koc/f_d) + Kim2019 L_Ph fit
python validation/nstem_gradient_check.py # 다중높이 줄기: Yamazaki 줄기 구배 재현
pip install pytest && pytest              # 테스트 (구조/질량보존/QSPR/보정/API)
```

### 구조(SMILES)로 임의의 PFAS 파라미터화

`src/pfas_structure.py`는 **화학 구조**를 모델에 매핑해, 선별 13종이 아닌 *임의의* PFAS도 돌릴 수 있게 합니다.
RDKit이 SMILES를 파싱 → 기술자(descriptor) → 측정 read-across(알려진 congener → 보정 파라미터) 또는
조각 QSPR(신규 → provisional)로 `Compound`를 만듭니다.

```bash
pip install -r requirements-structure.txt        # RDKit (선택)
python src/pfas_structure.py                      # SMILES → 기술자 → Compound 데모
```
```python
import model_api as api
api.simulate_from_smiles("OC(=O)" + "C(F)(F)"*11 + "C(F)(F)F")   # 신규 C13 PFCA를 SMILES로 실행
```
참조: `docs/structure_input.md` (결합/화학종은 구조에서; `f_xy`는 순서만 — 절대 스케일은 fit 유지).

### 실제 토양 측 — HYDRUS-1D (Method A, 연결됨)

Method A의 토양 절반은 실제 HYDRUS-1D 엔진(벤더링된 `external/hydrus_source` FORTRAN을 `phydrus`로 구동)을 돌려
congener별 실제 공극수 궤적 `C_w^o(t)`를 만들어 식물 ODE를 구동합니다. 소스는 벤더링되어 있으므로(상위 submodule이
제한적 네트워크 정책 뒤에서 접근 불가라 de-submodule함) submodule init이 필요 없습니다 — 빌드 + `phydrus` 설치만 하면 됩니다:

```bash
cp external/hydrus_source/makefile external/hydrus_source/source/
(cd external/hydrus_source/source && make)     # gfortran 필요
pip install phydrus
python src/soil_hydrus.py                       # congener별 공극수 요약
python validation/hydrus_coupled_run.py         # 전체 토양→식물 + 그림/CSV
```

약하게 흡착하는 단쇄(Kd≈0.01–0.15, `Koc·f_oc`)는 담수 동안 거의 0으로 leach되어, 상수-`Cwo` 가정은
**낟알/짚 BAF를 ~2–4배 과대예측**합니다(PFBA 낟알 2.07→0.43). 강하게 흡착하는 장쇄(Kd≳7)는 완충됩니다.
엔진/`phydrus`가 없으면 HYDRUS 테스트는 자동 skip됩니다. Claude Code on the web에서는 **SessionStart 훅**
(`.claude/hooks/session-start.sh`)이 엔진 빌드 + 의존성 설치를 자동으로 합니다.

**엔진 없이 시간 변화 노출.** `simulate(cwo_profile="flooded")`는 해석적 Freundlich 희석+leaching `C_w^o(t)`
(단쇄 leach, 장쇄 완충; 시즌 평균을 `Cwo`로 정규화)를 주며, `k_leach` 한 노브로 HYDRUS 방향을 재현합니다
(`python validation/cwo_profile_check.py`). 베이지안 역추정 데모(`python validation/bayesian_inverse_demo.py`)는
조직 `C(t)`가 무엇을 식별할 수 있고 없는지 보여 줍니다: 수송을 고정하면 노출(`Q_TP`-scale, `Cwo`-level)은 복원되나
`Q_TP·f_xy`는 곱 ridge — `Q_TP`/`Cwo`를 절대적으로 고정하려면 독립 측정(수액/공극수 프로브)이 필요합니다.

`reproduce_demo.py`는 `params/parameters.json` + `src/`를 불러와 12종 전부에 대해 4구획 ODE를 풀고,
예측 vs 관측 뿌리/짚/낟알 BAF를 출력합니다.

## 대화형 앱 — 시각화 도구

> 📖 **앱 사용 매뉴얼(일반인~전문가): [한국어](docs/MANUAL_KR.md) · [English](docs/MANUAL_EN.md)** — 두 모드, 모든 탭, 데이터 표/CSV 형식, 베이지안 역추정, 결과 해석, FAQ.

흙 + 벼를 실제 비율로 그리고 **각 구획을 PFAS 축적량 색(heat colormap)**으로 칠하는 Streamlit 대시보드입니다
(한 철을 슬라이더로 스크럽). 대화형 Plotly 시계열(hover/zoom/legend-toggle)도 함께 제공합니다:

```bash
pip install -r requirements.txt     # 앱 전체 스택 (numpy/scipy/matplotlib + streamlit/plotly/pandas)
streamlit run app.py
```

**웹 배포(Streamlit Community Cloud — 무료, URL 공유):** 레포를 GitHub에 올린 뒤
[share.streamlit.io](https://share.streamlit.io)에서 레포/브랜치/`app.py`를 선택 → Deploy.
`requirements.txt`가 전체 앱 스택이라 그대로 배포됩니다(RDKit/HYDRUS는 선택이며, 엔진이 없으면 live-HYDRUS 모드는 자동 숨김).
`requirements-app.txt`에는 선택 extras(kaleido PNG 내보내기, phydrus)만 있습니다. 참조: `docs/deploy.md`.

> 처음 화면은 **일반인 모드(한국어)** 입니다. 사이드바 토글 **🔬 전문가/고급 모드 (Expert / advanced)** 로
> 전체 연구용 인터페이스(영어)로 전환합니다. 자세한 사용법은 매뉴얼 참조.

**🗺️ 식물·토양 지도** — 대표 화면: 수염뿌리 벼(아치형 culm, 긴 잎, 늘어진 낟알 이삭), 각 기관을 공통 컬러바로
농도(또는 BAF)에 따라 칠합니다. 날짜 슬라이더(또는 ▶)로 *언제 어디서* PFAS가 쌓이는지 봅니다 — 잎은 물관-종단,
낟알은 phloem-fed, 뿌리는 음이온을 가둠.

**5가지 입력 모드**(사이드바 "Data source")가 노출 공간 전체를 다룹니다:

| 모드 | 공극수 `Cwᵒ(t)` 출처 | 언제 |
|---|---|---|
| **Model (parametric)** | 직접 지정한 상수 | 빠른 what-if / 교육 |
| **HYDRUS / CSV drivers** | HYDRUS-1D/Phydrus 결과(`t,Cwo,Qtp,M_*` CSV) | 보정된 토양 모델이 있을 때 |
| **Run HYDRUS-1D (live)** | 앱에서 실행한 실제 HYDRUS 엔진 | 앱에서 HYDRUS를 돌리고 싶을 때(빌드 필요) |
| **Soil inventory** | 총 토양 적재량 역산(Freundlich) | 토양 PFAS는 알지만 공극수는 모를 때 |
| **Biomonitoring** | 측정 공극수 값(HYDRUS 불필요) | 현장 조직+물 농도가 있을 때 |

(앱에는 위 5개 + **Custom tables(Cwᵒ + 성장)** 모드가 있어 성장·오염 시계열을 직접 표로 넣을 수 있습니다 — 매뉴얼 §6.2.)

**live HYDRUS-1D** 모드는 실제 엔진(`external/hydrus_source`에서 `phydrus`로 빌드)을 한 철 논 모델로 돌려
congener별 `Cwᵒ(t)`를 만듭니다(단쇄 leach, 장쇄 완충). 엔진을 자동 감지하고 없으면 빌드 단계를 안내합니다.
참조: `src/soil_hydrus.py`, `docs/visualization_tool.md`.

기타 탭: 조직 동역학, **토양 & 드라이버**(`Cwᵒ(t)`, `Q_TP(t)`, `M(t)`, Freundlich 등온선, 깊이 프로파일),
BAF vs 관측/측정, 사슬 길이 추세, congener 비교, 그리고 HYDRUS-1D 입출력 매핑·바이오모니터링 경로를 설명하는 **About** 탭.
계산은 `src/model_api.py`(`simulate(...)`, 토양/드라이버/바이오모니터링 헬퍼), Plotly 그림은 `src/plots.py`
(`fig_plant_schematic` …) — 둘 다 UI 비의존이며 테스트로 커버됩니다. 바로 불러올 예시는 `examples/`에 있습니다.
전체 가이드: `docs/visualization_tool.md`.

---

## 구조(layout)

```
pfas_rice_model/
├── README.md / README_KR.md      ← (영문/국문) 진입점
├── build_parameters.py           ← 소스 표에서 params/parameters.json 조립
├── reproduce_demo.py             ← 독립 실행 ODE 재현(엔트리 포인트)
├── app.py                        ← Streamlit 앱(일반인 한국어 / 전문가 영어)
├── src/                          ← 모델 코드
│   ├── pfas_rice_plant_module_4pool.py        4구획 식물 ODE (basis-A)            ← 정본(CANONICAL)
│   ├── pfas_rice_plant_module_4pool_surf.py   + K_surf (Fe/Mn-plaque 표면 풀)
│   ├── pfas_rice_plant_module_5pool.py        + 명시적 리그닌 풀
│   ├── pfas_rice_plant_module_nstem.py        N개 직렬 줄기 세그먼트(다중높이; GAP-B 보정)
│   ├── pfas_rice_plant_module.py              import 별칭 → 4pool_surf (삭제 금지)
│   ├── soil_paddy.py                          토양↔공극수 (Freundlich)            ← legacy redox 부호
│   ├── soil_paddy_redox_corrected.py          W3-보정 redox (이것을 사용)
│   ├── soil_hydrus.py                         실제 HYDRUS-1D → C_w^o(t),Q_TP(t) (Method A)
│   ├── model_api.py                           UI 비의존 래퍼(simulate/드라이버/표/역추정/내보내기)
│   ├── plots.py                               Plotly 그림(한/영 lang 지원)
│   ├── calibration.py                         BAF→파라미터 fit 기계
│   └── literature_params.py                   문헌 QSPR/앵커(인용) + Kim2019 BAF
├── params/                       ← 파라미터 (parameters.json = 정본)
├── data_obs/                     ← 검증용 관측 BAF/TF (Yamazaki, Li2025, …)
├── validation/                   ← 재현 스크립트 + 산출물 + 그림
└── docs/                         ← 문서 (MANUAL_KR/EN, OVERVIEW_KR, VALIDATION_*_KR, GAP_A/B, …)
```

---

## 한 화면 요약(모델)

결합인자(**basis A, 신선중** — 가장 중요한 규약 하나):

```
B_k = θ_fw + (1 − θ_fw) · ( f_prot·K_prot + f_PL·K_PL + f_cw·K_cw )      [L/kg fw]
```

θ_fw = 조직 수분 분율; f_* = **건조중** 질량 분율; K_* = 분배계수 [L/kg pool-dw].
`(1 − θ_fw)` 인자는 필수 — 빼면(legacy `Bk_table_S5.csv`) B_k를 ~3배 과대평가하고 풀 비중을 망칩니다.
`f_cw`는 전체 세포벽(다당류 + 리그닌), 대응 K는 `K_cw_wholecw_<organ>`.

뿌리 유입 = GHK 전기확산(음이온 배제 e^N ≈ 107 @ E_m = −120 mV) + 포화 Michaelis–Menten 운반체(배제 극복).
물관 적재 = `f_xy · Cw_root`; 낟알은 phloem 우세(`L_Ph · Cw_leaf`).

---

## ⚠ 가장 중요한 한 가지: 어떤 `f_xy`인가

`parameters.json`은 **두 개**의 `f_xy(n)`을 갖습니다:

| 필드 | 무엇인가 | 용도 |
|---|---|---|
| **`f_xy_recommended`** | **monotone** 물리적 TSCF(이론 유도, 교차필드 검증). C4 0.79 → C12 0.003. | 파라미터 인용/보고; 물리적으로 옳은 값 |
| `f_xy_W2fit` | Yamazaki 수송 fit; **C10+에서 가짜로 상승**(0.08→0.67) | *현재* ODE 구조로 Yamazaki를 재현할 때만 |

둘은 장쇄에서 갈립니다 — W2 fit(포화, 3 param/3 obs)이 **모델되지 않은 줄기 축적 구배**를 `f_xy`로 흡수하기 때문.
이론(Trapp GHK + Briggs LFER)과 **물-독립 교차필드 TF**는 monotone *방향*을 요구 → `f_xy_W2fit`의 장쇄 상승은
대부분 단일구획 아티팩트. **파라미터로는 `f_xy_recommended`를 사용**하세요. 참조: `docs/DELIVERABLE_GAP_B_fxy.md`.

**업데이트 — 다중높이 줄기 보정 구현됨**(`src/pfas_rice_plant_module_nstem.py`, `validation/nstem_gradient_check.py`).
줄기를 N개 직렬 세그먼트로 분해(증산 인출 + 반경 교환 + 성장 희석)하면 **monotone f_xy가 관측된 Yamazaki 줄기 구배(PFCA)를
재현**합니다(단쇄는 위로 농축, 장쇄는 평탄/하강; 전환은 `B* ~ Q_s/(M_s·μ_s)`로 결정). 리뷰 주의: (i) *절대* 전환점과
f_xy 스케일은 **측정 `Q_TP(t)`/`M_s(t)`**가 필요(placeholder 증산이 ~5배 높음); (ii) **PFOS/PFSA**는 높은 결합에도
위로 이동 — 결합 기반 monotone f_xy로는 못 잡으므로 PFSA 전용 수송항이 여전히 필요합니다.

---

## 상태(status)

> 리뷰 + 다중높이 줄기 작업 이후의 정직한 상태(이전 "CLOSED/validated" 라벨을 여기서 엄격히 조정):

- **결합 `B_k`** — **측정** congener별 `K_PL`(Chen 2025 K_MW, vs Droge 2019)과 `K_prot`(Zhou 2025 dialysis
  `K_prow`; 콩=식물, BSA=동물) 위에 구축 — `docs/literature_db/raw_si/`, `src/literature_params.py`. basis-A 신선중 규약.
- **GAP A (K_cw)** — 값 제공됨, 단 **앵커링(DFT 사다리 + 측정 리그닌), 직접 측정 아님** — 장기적 최약점이며,
  막 우세 장쇄에는 `K_cw`가 부차적 풀. `docs/DELIVERABLE_GAP_A_Kcw.md`.
- **GAP B (f_xy)** — *형상* 해결(monotone; 단쇄 상한 ≈0.8 Felizeter TSCF 앵커)되고, **다중높이 줄기가 monotone f_xy로
  PFCA 줄기 구배를 재현**. **완전 종료 아님**: 절대 스케일/전환점은 측정 `Q_TP(t)`/`M_s(t)`(과제 2), PFSA는 별도 수송항(과제 3) 필요.
- **검증 주의** — `reproduce_demo.py`의 log10 RMSE 0.029는 **포화 W2 fit**(congener당 3 수송 param을 3 관측 BAF에 fit →
  재현 보장, 예측 검증 *아님*). 진짜 out-of-sample 증거는 물-독립 **교차필드 TF**(monotone 방향)와 **nstem 구배 방향**(PFCA).
  `docs/H8_handoff_S6_final.md`. **전체 검증(국문) + 그림: `docs/VALIDATION_KR.md`** — 보정(Yamazaki, 재현) vs
  out-of-sample(Kim 2019 낟알; Li 2025 TF, 비결정적); 그림은 `python validation/validation_summary.py`.
- **독립 rigor 감사(sci-adk)** — 동결된 기준으로 모델을 판정하는 7개 적대적 run(`sci_adk_review/`): 형식적 기반 SUPPORTED
  (질량보존, 음이온 배제, 토양 QSPR, SMILES read-across), 순진한 예측 주장 REFUTED(0.029는 in-sample), 구조적 적합성
  SUPPORTED(ORYZA2000 biomass), 그리고 **두 독립 데이터셋(Tang 2026, Kim 2019)에 걸쳐 out-of-sample로 일반화되는
  지질 매개 적재 메커니즘**. 통합 합성 논문: `sci_adk_review/runs/pfas-rice-consolidation/paper/draft.tex`
  (재빌드: `python sci_adk_review/build_consolidation.py`; 국문 서술: `sci_adk_review/FINDINGS.md`).
- **Tier-1 fit** — `src/literature_params.py`가 Kim 2019 PFOA 낟알 BAF에 `L_Ph`를 fit(4.43 L/kg 일치).
- **시각화 도구** — `app.py`(+ `src/model_api.py`, `src/plots.py`): 식물/토양 축적 colormap + 5(+1)개 노출 모드.
  일반인 한국어 / 전문가 영어. `docs/visualization_tool.md`, `docs/MANUAL_KR.md`.
- **토양 측(Method A) — 실제 HYDRUS-1D**(`src/soil_hydrus.py`): 엔진 컴파일·연결됨. congener별 `C_w^o(t)`가 식물 ODE
  (및 앱 live 모드)를 구동. 단쇄 leach → 상수-`Cwo`가 낟알/짚 BAF ~2–4배 과대예측; 장쇄 완충. 미빌드 시 테스트 skip.
- **선택형 지질 매개 적재** — `simulate(lipid_loading=True)`가 K_PL-게이트, B-독립 물관/체관 항(`g_xy`/`g_ph`; 기본 0)을
  추가해 고결합 장쇄가 지상부에 도달. in-sample/탐색적(`docs/fxy_longchain_lipid_exploration.md`).
- **테스트** — `pytest`. 전체 스택(RDKit + 빌드된 HYDRUS 엔진 + phydrus)에서 모두 통과(앱: model_api/plots/구조/역추정/표입력 포함).
  bare clone에서는 RDKit 없으면 구조 테스트, 엔진 미빌드면 HYDRUS 테스트가 skip.

**미해결(데이터 한계, 모델링 아님):** 벼(밀 아님) congener별 뿌리 아세포분획 → 막 비중/α; 신뢰할 수 있는 congener별
공극수/수경 RCF → 표면 테스트 + f_xy 절대 스케일; 측정 Q_TP(t), M(t) → f_xy 절대; 직접 K_cw_poly + 벼 세포벽 단당 조성; in-situ 논 E_m.

**미해결(모델링, 지금 가능):** (1) *물리적* monotone f_xy가 장쇄 짚을 재현하도록 다중높이 줄기 구획(현재 W2 fit이 보상);
(2) ~~통합 토양→식물 run~~ **완료**(실제 HYDRUS-1D 엔진) — 남은 것: 측정 현장 담수 일정, 무산소/담수 흡착, 사용자 현장 토양;
(3) f_PL(0.01–0.02, 2배 불확실) 불확실성 전파.

**표준화할 설정:** 뿌리 θ = 0.90(측정 0.90–0.92), 뿌리 f_PL = 0.015, 낟알 θ는 단계 의존(0.14 수확 / 0.30 등숙).
결론은 이 값들에 강건합니다.

> [!WARNING]
> 본 도구/모델은 **연구·교육용**입니다. 규제·식품안전·건강 결정의 근거로 사용하지 마세요.
