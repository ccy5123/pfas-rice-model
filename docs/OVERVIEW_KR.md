# PFAS–벼 흡수 모델 종합 문서 (기능·검증·데이터·실험·표기)

> 이 레포 전체를 한눈에 보는 문서입니다. **무엇을 하는가**, **얼마나 믿을 수 있는가**, **무엇이 부족한가**,
> **그래서 어떤 실험이 필요한가**, 그리고 **기호 정리**를 담았습니다. 수치는 2026-10-25 수정 모델 기준이고,
> 각 수치를 낸 스크립트와 자세한 내용은 `CLAUDE.md` §6, 모든 입력의 출처·민감도는 `docs/ASSUMPTIONS.md`에 있습니다.
> 깊은 유도는 `docs/pfas_rice_compartmental_model.tex`입니다.
>
> 🖥️ **앱 사용법**(일반인~전문가, 한국어)은 **[`docs/MANUAL_KR.md`](MANUAL_KR.md)** 를 보세요.

![모델 구조 모식도](figures/model_schematic.png)

---

## 1. 한 줄 요약

벼(*Oryza sativa*)에서 **영구 음이온성 PFAS**의 생물축적을 푸는 **4구획 동역학 흡수 모델**입니다.
Trapp/Brunetti의 **Dynamic Plant Uptake(DPU)** 틀을 **이온성 유기화합물(IOC)** 로 확장했고, 토양 쪽은
**HYDRUS-1D**와 단방향 결합합니다. **결합은 독립 측정값으로 받치고, 이행(transport)은 한 현장(Yamazaki 2023,
쓰쿠바 2020)에 보정**했습니다. 보정에 쓰지 않은 세 데이터셋으로 옮겨 보면 기본 설정은 log10 1.26–2.33 빗나가고,
지질 결합 경로를 넣으면 그 오차가 크게 줄어듭니다. *예측 검증이 끝난 모델은 아닙니다.*

---

## 2. 모델이 하는 일

### 2.1 과학적 핵심
- **PFAS = 영구 해리 음이온**(`pKa` 매우 낮음 → `f_d ≈ 1`). 중성화합물의 **Kow/Briggs 분배 코어는 적용되지 않습니다.**
- **뿌리 흡수 `j_R`는 하이브리드**: 막전위가 안쪽 음(−)이라 음이온을 **밀어내는** 전기확산(GHK; `e^N ≈ 107`)에
  이를 이겨내는 **포화 운반체**(Michaelis–Menten)가 더해집니다.
- **물관(xylem)** 은 뿌리의 자유 음이온을 적재계수 `f_xy`로 싣고 증산류를 따라 줄기·잎·낟알로 나눕니다
  (잎/낟알 몫은 건물 기준 표면적으로 나눔).
- **체관(phloem)은 잎에서 두 갈래**입니다. 낟알로는 낟알 건물 생장에 드는 수액만큼(`T_C,Ph`), 뿌리로는
  증산량의 비율 `φ`만큼 갑니다. 농도는 둘 다 `L_Ph · C_leaf/B_leaf`.
- **잎·낟알은 종착 축적기**: 배출이 성장 희석뿐이고 성숙기에 `μ→0`이라 정상상태가 없습니다.
  최종 농도 = 시간적분/최종질량이므로 **동역학 모델이 필수**입니다.
- 낟알은 체관급이라 약산 pH 이온트랩이 **적용되지 않습니다**(`f_n ≈ 0`). 적재는 운반체/채널(`L_Ph`).
- 대사 `γ_k ≈ 0`(난분해), 공기 교환 off(`K_AW ≈ 0`).

### 2.2 핵심 관계식 (basis-A, 생체중 기준)
```
B_k = θ_fw + (1−θ_fw)·(f_prot·K_prot + f_PL·K_PL + f_cw·K_cw)     [L/kg fw]
C_k = B_k · C_w,k          BAF = C_k / C_w^o          TF = C_k / C_root
```
`(1−θ_fw)`는 **건물→생체중 변환**입니다(빼면 `B_k`를 약 3배 과대평가). 풀이는 BDF stiff ODE이고 질량이 정확히 보존됩니다.

### 2.3 무엇과 비교하나
- **기관 질량 `M_k`는 생체중**입니다. 작물 모델(ORYZA2000 등)은 건물을 내므로 `B_k`와 같은 `θ_fw`로 변환합니다(`src/organ_mass.py`).
- **관측과 같은 건물중 기준**으로 비교합니다. 수확 시 죽은 잎은 작물 모델의 잎 고사로부터 재구성해 잎·짚에 포함합니다(`src/harvest.py`).
- **기본 구동 조건은 보정 현장의 재배 연도**입니다(`src/site_scenario.py`): JMA 쓰쿠바 2020 기상, Yamazaki 이앙일에 맞추고
  출수·수확일로 생육 단계를 보정한 ORYZA2000 잠재생산 이식(IR72.D92), FAO-56 항으로 계산한 ORYZA 증산.
  `scenario="generic"`은 기후값과 고른 생장 드라이버를 씁니다.

### 2.4 무엇을 입력하면 무엇이 나오나
| 입력(화합물) | 입력(시나리오/작물/토양) | 출력 |
|---|---|---|
| congener 이름, SMILES, 또는 중성물질 log Kow → `K_PL/K_prot/K_cw, f_d, f_xy` | 공극수 `C_w^o(t)`, 증산 `Q_TP(t)`, 생장 `M_k(t)`, `E_m` | 조직별 농도 `C_k(t)`, **BAF/TF**(건물중), 경로별 기여 분해 |

결합·해리는 구조에서 예측되지만 **`f_xy`(이행)는 구조에서 나오지 않습니다**(`docs/structure_input.md`).

### 2.5 모듈 맵
| 파일 | 역할 |
|---|---|
| `src/pfas_rice_plant_module_4pool_surf.py` | **정본** 4구획 ODE 코어(+표면 pool), 체관 두 갈래 `phloem_flows` |
| `src/pfas_rice_plant_module_nstem_leaf.py` | 재분배 지상부: N분절 줄기 + 잎, 기관별 증산 침착·보유 |
| `src/site_scenario.py` · `src/oryza2000.py` · `src/fao56.py` | 보정 현장 연도 구동자(기상·작물·증산) |
| `src/organ_mass.py` · `src/harvest.py` | 건물→생체중 변환 · 수확 관측 연산자(고사엽, 건물중) |
| `src/soil_hydrus.py` | **실제 HYDRUS-1D** 엔진(phydrus) → `C_w^o(t), Q_TP(t)` (Method A) |
| `src/literature_params.py` | 문헌 QSPR/anchor(인용): `Koc, K_PL, K_prot, f_d`, Kim 2019 표 |
| `src/pfas_structure.py` | **SMILES → Compound** (RDKit; read-across + QSPR) |
| `src/neutral_dpu.py` | 중성 유기물 경로(같은 ODE, z=0, Briggs `K_PW`·TSCF) |
| `src/model_api.py` | UI 무관 API: `simulate()`, `yamazaki_harvest_baf()`, `simulate_nstem_leaf()`, `simulate_neutral()`, 역추정 |
| `app.py` + `ui/` + `src/plots.py` | Streamlit 앱(일반 모드 한국어, 전문가 모드 영어) |

---

## 3. 검증 — 보정과 검증을 구분해서

### 3.1 보정 현장(Yamazaki 2023) — 재현이지 검증이 아님
| 설정 | log10 RMSE | 의미 |
|---|---|---|
| Fit A: congener마다 `f_xy, L_Ph, κ_d` 적합 | **0.162** (뿌리 0.000 · 짚 0.035 · 낟알 0.287) | 포화 적합(파라미터 3 / 관측 3) — 잘 맞는 것이 보장됨 |
| Fit B(기본): 단조 물리적 `f_xy`, `L_Ph·κ_d`만 적합 | 짚 **1.40** (예측), 낟알 1.56 | 짚은 예측. PFOA(C8) 이상과 PFOS의 낟알은 `L_Ph = 1` 상한에서도 모자람 |
| 지질 결합 전역 상수 6개 | 0.207 (29셀, PFDoDA 제외) | 표본 내 |

Fit A에서도 PFBA·PFPeA는 `L_Ph` 하한에서 낟알이 넘치고(물관만으로 과공급), PFDA·PFUnDA·PFOS·PFDoDA는 `L_Ph = 1`에서 낟알이 모자랍니다.
기본 설정은 PFDoDA 현미 BAF를 측정값(53)보다 약 2,300배 낮게 계산합니다.

### 3.2 독립 데이터로 옮기기 (보정에 쓰지 않은 데이터)
| 데이터 (스크립트) | 기본 | congener별 적합 | 지질 결합 | 재분배 지상부 + 지질 |
|---|---|---|---|---|
| Tang 2026 부위별 TF, 0.1 µg/g (`oos_tang_lipid.py`) | 1.50 | — | **0.64** | — |
| Kim 2019 현미 BAF, PFOA 제외, 점 / [0, MDL] 구간 (`kim2019_censored.py`) | 2.33 / 2.05 | 0.73 / 0.55 | 0.56 / 0.38 | 0.47 / 0.15 |
| Li 2025 생육단계별 TF, 지상부 / 낟알 (`li2025_stage_resolved.py`) | 1.26 / 1.68 | 0.96 / 1.30 | 0.90 / 1.24 | **0.36 / 0.51** |

- **자유 음이온만으로는 옮겨지지 않습니다.** 기본 설정의 오차 중앙값 log10 1.50(×31)을 앱 일반 모드의 불확실성 범위로 씁니다.
- **지질 결합 경로는 옮겨집니다.** Yamazaki에만 맞춘 상수로 Tang·Kim 오차를 절반 이하로 줄였고, 재분배 지상부와 함께 쓰면
  Li 2025의 사전 등록 판정을 두 종말점 모두에서 통과했습니다. Kim에서 두 지질 변형의 순위는 불검출 처리에 따라 바뀝니다.
- **Kim 2019 PFOA**: 사이트별 중도절단 기하평균 9.94 L/kg(90% 4.57–20.2). 발표된 평균비 4.44는 하수처리장 인근 한 필지가
  83%를 차지한 값입니다. 모든 모델 변형이 이 구간보다 낮습니다.
- **Tang 2026**: 이 데이터에 직접 맞춘 `f_xy`(PFOA 0.93, PFOS 0.999 — 상한, GenX 0.059)로도 0.34. 배유를 백미 조성으로 읽으면 0.64 → 0.83으로 나빠집니다.
- 전이 데이터셋은 각 현장의 기상·일정이 레포에 없어 모두 쓰쿠바 2020 조건으로 계산했습니다.

### 3.3 메커니즘 점검
- **운반체 vs 우회로**(`carrier_vs_bypass.py`, 사전 등록): 아무것도 안 넣음 2.301 · 막 탈분극 −90 mV 1.890 ·
  운반체(고정) 1.227 · 운반체 Vmax 재적합(×8) 1.078 · 아포플라스트 우회(`g_apo` 50) 1.093. *추가는 필요*하지만
  *어느 것인지는 미결*(bootstrap 0.368). congener별로 둘 다 같은 사슬길이 보정을 요구합니다(+0.90/+0.92) → 결합·격리 쪽 문제.
- **용량 시리즈**(Tang 5용량): 규칙 1/3 충족 → 운반체는 *불리하지만 기각은 아님*. 측정만큼 평평하려면 `Km ≥ 150–500 µg/L`(고정값의 30–100배).
- **식별성**(`bayesian_inverse_demo.py`): 수송을 고정하면 노출(`Q_TP` 규모, `Cwo`)은 복원됩니다(cond 18). `Q_TP`와 `f_xy`는
  강하게 맞바뀌고(corr −0.97), `Cwo` 대 뿌리 전도도는 능선(cond 1.6e6)입니다 → 공극수·물관 수액 독립 측정이 필요.
- **θ_fw**: 뿌리 중간(INTERMEDIATE), 줄기·잎·낟알 작음.
- **HYDRUS-1D**: 약흡착 단쇄는 담수기에 용탈(PFBA 낟알 BAF 4.42 상수 노출 → 0.54), 장쇄는 완충.

### 3.4 결과를 크게 움직이는 미측정 입력 (`docs/ASSUMPTIONS.md`)
운반체 `Km` ×100(모든 부위 1.26 log ≈ 18배) · 막전위 −90 mV(0.13 log) · 죽은 잎이 수확 짚에 남는지(짚 0.17) ·
증산류 중 낟알 몫(낟알 0.28) · 낟알 생장당 체관 수액량 `T_C,Ph`(낟알 0.16) · 체관 재순환 `φ`(짚 0.05).
육묘·재식밀도·포기당 주수·CO₂·출수일·수확 시점은 크지 않습니다.

### 3.5 중성 유기물(맞춘 것 없음)
Liu 2023 벼 뿌리 분배 0.213(평형 0.206) · Ge 2017 부위별 TF 0.311 · Kodešová 2019 0.216 · Li 2019 수경 0.546 ·
Li 2019 토양(n=376) 0.550 · Briggs 1983 보리 줄기 0.299. 자세히는 `docs/neutral_dpu_validation.md`.

### 3.6 정직한 결론
- 구조(음이온 배제, 질량 보존, 종착 축적)는 검증됐고, 단쇄의 짚 > 뿌리 방향도 현장과 같습니다.
- 절대 수치는 대략치입니다(기본 설정 ×31). **장쇄 낟알은 가장 불확실**합니다. 현장 측정은 C10–C12에서 올라가는데
  자유 음이온 구조로는 닿지 않고, 지질 결합 경로는 가장 잘 맞는 후보이지만 실험으로 확정되지 않았습니다.
- 기본 설정을 지질 결합 + 재분배 지상부로 바꿀지는 사용자 결정 사항입니다(Yamazaki 표본 내, Kim 순위 불안정).

---

## 4. 부족한 데이터 (공백)

| # | 공백 | 영향 | 비고 |
|---|---|---|---|
| 1 | **Yamazaki 2023 SI 원자료** (S16 수질, S18/S19 조직, MQL 표시, 재식밀도, 육묘) | 보정의 관측값 자체 | 현재 그림·표에서 재구성 |
| 2 | **운반체 동역학** (`Vmax/Km`, 억제제 실험) | 결과를 18배 좌우, 운반체 vs 채널 분리 불가 | 용량 시리즈가 `Km`을 크게 요구 |
| 3 | **조직·시간 분해 BAF/TF + 짝지은 공극수**, 특히 장쇄 낟알 | 장쇄 낟알 메커니즘 판별 | 온실 시계열 |
| 4 | **이삭(영과) 증산 몫** | 낟알 물관 몫(0.28 log) | 미측정 기본값 |
| 5 | **벼 체관 수액량 / 낟알 건물** (`T_C,Ph`) | 낟알 체관 유입(0.16 log) | 미측정 기본값 |
| 6 | **생육기 기관 수분 함량** (`θ_fw`) | 건물↔생체중 모든 단계 | 뿌리가 가장 민감 |
| 7 | **벼 뿌리 세포벽·plaque PFAS 흡착** (`K_cw`, `K_surf`) | 결합·격리 항 | 문헌 계수 없음 |
| 8 | **혐기/담수 흡착 등온선** | 담수기 `C_w^o` | 모두 호기 자료 |
| 9 | **전이 데이터셋의 현장 기상·재배 일정** (Tang, Kim, Li) | 전이 시험의 구동자 | 지금은 쓰쿠바 2020으로 대체 |
| 10 | **에테르/설폰아마이드 Koc·K_PL**, GenX 등 대체물질 벼 자료 | 구조 입력 신규 물질 | 단일 anchor |

---

## 5. 따라서 해야 하는 실험

| 우선 | 실험 | 해소하는 공백 | 모델에 들어가는 것 |
|---|---|---|---|
| **HIGH** | **억제제 ± 농도열 뿌리흡수 assay** | 2 | `Vmax/Km`, 운반체 vs 채널 |
| **HIGH** | **온실 논 시계열**(2품종, 시기별 뿌리/짚/왕겨/현미/백미 + 짝지은 공극수, 재배 일정 기록) | 3, 6, 9 | 진짜 표본외 검증, 장쇄 낟알 |
| **HIGH** | **혐기 배치 흡착**(제어 Eh) | 8 | 담수기 `C_w^o` (HYDRUS 입력) |
| MED | **물관 삼출액 `f_xy` + 이삭자루 체관 `L_Ph`·수액량** | 4, 5 | 이행 파라미터 직접 측정 |
| MED | **뿌리 세포벽 분획·plaque 배치 흡착 + 탈착** (사슬길이 × 머리기) | 7 | `K_cw`, 격리 항 |
| MED | **현장 미세전극 `E_m`** | — | GHK 배제항 검증 |

---

## 6. 기호 정리

> 단위: 시간 **day**, 수용액 농도 **µg/L**, 조직 농도 **µg/kg**(모델 내부 생체중, 보고는 건물중), 질량 **kg**, 유량 **L/day**.
> 기호는 `docs/pfas_rice_compartmental_model.tex`와 1:1 대응합니다.

### 6.1 상태·구획
| 기호 | 의미 | 단위 |
|---|---|---|
| `k = 1,2,3,4` | 구획: root, stem, leaf, grain(fruit) | – |
| `C_k(t)` | 구획 k 조직 농도 | µg/kg fw |
| `C_w,k` | 구획 k 자유(수상) 농도 = `C_k/B_k` | µg/L |

### 6.2 Tier-0 — 입력/기지값
| 기호 | 의미 | 단위 | 출처 |
|---|---|---|---|
| `C_w^o(t)` | 토양 공극수 자유 농도(노출) | µg/L | 상수(기본) · HYDRUS · Freundlich |
| `Q_TP(t)` | 증산 = 뿌리 물흡수 | L/day | ORYZA2000 Tp(현장 연도) / `forcing_rice`(generic) |
| `Q_g = (dM_4^dry/dt)·T_C,Ph` | 낟알로 가는 체관 유량 | L/day | 낟알 건물 생장 × 수액량 |
| `Q_r = φ·Q_TP` | 뿌리로 되돌아가는 체관 유량 | L/day | – |
| `M_k(t)` | 구획 생체중 질량 | kg | ORYZA2000 이식 → `organ_mass` |
| `μ_k = (dM_k/dt)/M_k` | 성장 희석률 | 1/day | – |
| `N = zF·E_m/RT` | 무차원 막전위(GHK); `e^N≈107`=배제 | – | – |
| `E_m` | 뿌리 막전위 | mV (~−120) | Wang 1994 |
| `f_d` / `f_n` | 해리(음이온) / 중성 분율 | – | `f_d≈1, f_n≈0` (Goss 2008) |
| `γ_k` | 대사율 | 1/day | `≈0` |
| `θ_fw` | 생체중 수분 분율 | – | `rice_tissue_params` |
| `f_prot, f_PL, f_cw` | 단백질/인지질/세포벽 **건물** 질량분율 | – | `rice_tissue_params` |
| `T_C,Ph` | 낟알 건물 생장당 체관 수액량 | L/kg dry | 미측정 기본값 |

### 6.3 Tier-1 — BAF로 식별(보정 대상)
| 기호 | 의미 | 단위 | 비고 |
|---|---|---|---|
| `B_k` | 결합인자 (`C_k=B_k·C_w,k`) | L/kg fw | basis-A |
| `f_xy` | 뿌리→물관 적재계수(TSCF 유사), ∈(0,1] | – | 단조 물리값(기본) 또는 적합 |
| `L_Ph` | 체관 적재계수 (∈[1e-5, 1]) | – | 적합 |
| `κ_d` (`kappa_d`) | 이온 막전도(lumped) | L/(day·kg) | 적합 |
| `φ` | 체관 재순환 분율 | – | 미측정 기본값 |

### 6.4 Tier-2 — 운반체 동역학(억제제 자료 필요)
| 기호 | 의미 | 단위 |
|---|---|---|
| `V_max,in / K_m,in` | 운반체 유입 Michaelis–Menten | µg/(kg·day) / µg/L |
| `V_max,out / K_m,out` | 운반체 유출 | 〃 |
| `g_apo` | 아포플라스트 우회 전도도(선택; `uptake="bypass"`) | L/(kg·day) |

### 6.5 Tier-3 — 결합/분배 QSPR(사슬길이 분해)
| 기호 | 의미 | 단위 | 출처 |
|---|---|---|---|
| `K_PL` | 인지질/막–물 분배 | L/kg lipid | Chen 2025 K_MW(측정) |
| `K_prot` | 단백질–물 분배 | L/kg prot | Zhou 2025 투석(측정) |
| `K_cw` | 세포벽–물 분배 | L/kg | **placeholder**(문헌 없음) |
| `K_surf` | 뿌리 표면/plaque 흡착 | L/kg | 현장 의존(기본 0) |

### 6.6 토양
| 기호 | 의미 | 단위 |
|---|---|---|
| `Koc` | 유기탄소 정규화 분배 | L/kg |
| `K_F, n` | Freundlich 용량·지수 (`S=K_F·C_w^n`) | – |
| `Kd = Koc·f_oc` | 선형 분배계수(HYDRUS 실행) | L/kg |
| `R = 1+ρKd/θ` | 지연계수 | – |
| `k_leach` | 해석적 담수 형상의 용탈 속도(HYDRUS에 congener별 보정) | 1/day |

### 6.7 출력 지표 + 구조 어댑터
| 기호 | 의미 |
|---|---|
| `BAF = C_k/C_w^o` | 생물축적계수(공극수 기준, 건물중으로 보고) |
| `TF = C_k/C_root` | 이행계수 |
| `n_perfluoroC, head_group, n_ether_O` | 구조 descriptor(SMILES) |
| `KPL_ETHER_LOG_OFFSET = −0.49` | 에테르당 K_PL 보정(GenX anchor) |
| `KOC_ETHER_LOG_OFFSET = 0` | 에테르당 Koc 보정 **(공백)** |
| `FXY_HEADGROUP_LN_OFFSET` | 머리기 f_xy 오프셋(sulfonate −1.1, ether −0.7) |

---

## 7. 재현 방법
```bash
pip install -r requirements.txt
python reproduce_demo.py                          # Yamazaki 현장 연도, Fit A (log10 RMSE 0.162); --rec = Fit B
python validation/refit_site_year.py              # 보정 다시 하기 (~1 h, 사전 등록)
python validation/oos_multidataset.py             # 전이 시험 요약 (Tang / Kim / Li)
python validation/li2025_stage_resolved.py        # Li 2025 생육단계별
python validation/kim2019_censored.py             # Kim 2019 불검출 반영
python validation/carrier_vs_bypass.py            # 운반체 vs 우회로
python validation/hydrus_coupled_run.py           # 실제 HYDRUS-1D 토양→식물 (엔진 빌드 시)
python validation/neutral_dpu_validation.py --lipid-source both   # 중성 유기물
streamlit run app.py                              # 시각화 도구
pytest                                            # 테스트
```

**관련 문서**: 현재 상태 `CLAUDE.md` · 입력 출처 `docs/ASSUMPTIONS.md` · 유도 `pfas_rice_compartmental_model.tex` ·
구조 입력 `structure_input.md` · 시각화 `visualization_tool.md` · 시료 채취 감사 `sampling_protocols_audit_KR.md`.
`docs/VALIDATION_*`, `docs/HANDOFF_*`, `docs/*_exploration.md`는 당시의 기록이며 현재 수치가 아닙니다.
