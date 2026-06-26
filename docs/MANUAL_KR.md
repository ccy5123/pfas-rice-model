# PFAS–Rice Uptake Explorer — 사용 매뉴얼

> 논(paddy)에서 자라는 벼가 오염된 물·흙에서 **PFAS("영원한 화학물질")**를 얼마나 흡수하고
> 어디(뿌리·짚·낟알)에 쌓이는지를 추정하는 대화형 도구의 **종합 사용 설명서**입니다.
> **일반 사용자**(정책·학생·일반)와 **전문가**(환경과학·연구) 모두를 위해 작성했습니다.
>
> - 일반 사용자는 **[1. 빠른 시작](#1-빠른-시작)** → **[4. 일반인 모드](#4-일반인-모드-사용법-한국어)** 만 읽어도 충분합니다.
> - 전문가는 추가로 **[5. 전문가 모드](#5-전문가-모드-사용법-english-ui)** 와 **[7. 입력 데이터 형식](#7-입력-데이터-형식-csv)**, **[9. 과학적 배경](#9-과학적-배경과-가정)** 을 참고하세요.

> [!WARNING]
> **연구·교육용 모델 — 예시 추정치일 뿐입니다.** 규제·식품안전·건강 판단이 **아니며**,
> 실제 노출·안전 결정에 **사용하지 마세요**. 결과는 모델 가정에 따른 illustration입니다.

---

## 목차

1. [빠른 시작](#1-빠른-시작)
2. [도구 개요](#2-도구-개요)
3. [두 가지 모드와 전환](#3-두-가지-모드와-전환)
4. [일반인 모드 사용법 (한국어)](#4-일반인-모드-사용법-한국어)
5. [전문가 모드 사용법 (English UI)](#5-전문가-모드-사용법-english-ui)
6. [공통 기능](#6-공통-기능)
7. [입력 데이터 형식 (CSV)](#7-입력-데이터-형식-csv)
8. [결과 해석 가이드](#8-결과-해석-가이드)
9. [과학적 배경과 가정](#9-과학적-배경과-가정)
10. [용어집](#10-용어집)
11. [문제 해결 (FAQ)](#11-문제-해결-faq)
12. [재현·검증·테스트](#12-재현검증테스트)
13. [인용·라이선스·참고](#13-인용라이선스참고)

---

## 1. 빠른 시작

### A. 배포된 웹 앱으로 (설치 불필요)
브라우저에서 배포 URL에 접속하면 바로 사용할 수 있습니다(Streamlit Community Cloud, `main` 브랜치 추적).
처음 화면은 **일반인 모드(한국어)** 입니다.

### B. 내 컴퓨터에서 실행
Python 3.9+ 환경에서:

```bash
git clone https://github.com/ccy5123/pfas-rice-model.git
cd pfas-rice-model
pip install -r requirements.txt     # 앱 전체 스택 (numpy/scipy/streamlit/plotly/pandas/rdkit/phydrus)
streamlit run app.py                 # 브라우저에서 http://localhost:8501 열림
```

- 5개 노출 모드 중 4개(Model / HYDRUS-CSV / Soil-inventory / Biomonitoring)와 SMILES 입력은 즉시 동작합니다.
- **Run HYDRUS-1D (live)** 모드만 컴파일된 엔진(gfortran)과 `phydrus`가 필요하며, 없으면 자동으로 숨겨지거나 안내됩니다.
- 정적 그림 PNG 내보내기는 선택 패키지 `kaleido`(+Chrome)가 있을 때만 활성화됩니다(없어도 CSV는 됩니다).

### 첫 30초 사용 흐름 (일반인 모드)
1. 왼쪽 사이드바 **① 화학물질 선택**에서 PFAS를 고릅니다(예: PFOA).
2. **② 오염 정도**에서 낮음/중간/높음을 고릅니다.
3. 본문 **🗺️ 어디로 가나** 지도에서 색이 진한(뜨거운) 부위 = PFAS가 많은 곳을 봅니다.
4. 헤드라인 카드(뿌리/짚/낟알 속 µg/kg)와 한 줄 요약으로 결과를 읽습니다.

---

## 2. 도구 개요

이 도구는 벼의 **4구획 동역학 흡수 모델**(뿌리 → 줄기 → 잎 → 낟알)을 화면으로 보여 줍니다.

- 흙(논 토양)과 벼를 실제 비율로 그리고, 각 구획을 **축적량 색(heat colormap)**으로 칠합니다.
- 한 철(이앙→수확) 동안의 시간 변화를 슬라이더/재생으로 볼 수 있습니다.
- **PFAS 한 종 + 오염 수준**을 주면 → 뿌리·짚·낟알 농도와 **축적 배수(BAF)**를 추정합니다.
- **거꾸로**: 측정한 조직 농도에서 **토양수 오염 수준을 역추정**(불확실성 범위 포함)할 수도 있습니다.

대상 화학물질은 **선별된 13종**(PFCA C4–C12: PFBA·PFPeA·PFHxA·PFHpA·PFOA·PFNA·PFDA·PFUnDA·PFDoDA,
PFSA C4/C6/C8: PFBS·PFHxS·PFOS, 에터 PFAS: GenX)이며, 전문가 모드에서는 **임의의 PFAS를 SMILES 구조로** 넣을 수도 있습니다.

> 계산 로직은 UI와 분리되어 있습니다: `src/model_api.py`(계산), `src/plots.py`(그림), `app.py`(화면).
> 모델 수식·파라미터의 자세한 내용은 `docs/`(특히 `docs/OVERVIEW_KR.md`, `docs/visualization_tool.md`)를 참고하세요.

---

## 3. 두 가지 모드와 전환

| | 일반인 모드 (기본) | 전문가 모드 |
|---|---|---|
| 언어 | **한국어** | 영어(English UI) |
| 대상 | 정책·학생·일반 | 환경과학·연구 |
| 입력 | 화학물질 + 오염 수준(낮음/중간/높음) | 5+1 노출 모드, SMILES, 모든 파라미터 |
| 탭 | 5개(쉬운 말) | 9개(전문 용어) |
| 기호 노출 | 없음(BAF/Cwᵒ/f_xy/eᴺ 숨김) | 전부 노출 |

**전환 방법:** 사이드바 맨 위의 토글 **🔬 전문가/고급 모드 (Expert / advanced)**.
- **끄기(기본)** = 일반인 모드(한국어)
- **켜기** = 전문가 모드(영어) — 일반인 모드의 모든 기능 + 연구용 전체 인터페이스가 복원됩니다(아무것도 사라지지 않음).

---

## 4. 일반인 모드 사용법 (한국어)

### 4.1 사이드바
- **① 화학물질 선택**: 13종 중 하나. 이름 옆 설명(예: "PFOA — 대표적 '영원한 화학물질' 카복실산(C8)").
  사슬이 길수록(C 숫자가 큼) 대체로 식물에 더 잘 달라붙습니다.
- **② 오염 정도**: 토양수에 녹아 있는 PFAS 양의 프리셋
  - 낮음 = 0.1 µg/L, 중간 = 1 µg/L, 높음 = 10 µg/L
  - 높을수록 식물로 더 많이 들어갑니다.
- **📋 내 데이터 표 사용**(선택): 직접 만든 **성장 곡선**과 **시간별 토양수 오염**을 표로 입력 → [6.2 참조](#62-내-데이터-표-성장--토양수)
- 토글 아래 안내대로, 전체 연구 인터페이스가 필요하면 전문가 모드를 켜세요.

### 4.2 헤드라인(요약 카드)
- **뿌리 속 / 짚(줄기+잎) 속 / 낟알(먹는 쌀) 속** — 각 부위의 PFAS 농도(µg/kg)
- 그 아래 한 줄 요약: "선택한 오염 수준에서, 낟알에 약 X µg/kg, 토양수의 약 N배, 대부분은 ○○에 남음"

### 4.3 탭별 사용법
- **🗺️ 어디로 가나** — 벼·논 그림. 색이 진할수록 PFAS가 많음.
  - **이앙 후 일수** 슬라이더로 특정 날짜의 분포를 봅니다.
  - **▶ 한 철 재생** 체크 시 한 철 전체를 애니메이션으로 봅니다.
- **📈 시간에 따른 축적** — 부위별 농도의 시간 변화. **낟알**은 개화 무렵 형성된 뒤부터 흡수가 시작됩니다.
- **📊 얼마나 쌓이나** — 수확 시 부위별 최종 농도 막대.
  - 아래 **🔬 실제 측정값과 비교 (Yamazaki 2023)** 펼치면 모델 vs 실측을 비교합니다(막대가 비슷할수록 잘 맞음).
- **🔎 거꾸로 추정** — 실험실 측정값에서 토양수 오염을 역추정 → [6.1 참조](#61-거꾸로-추정-베이지안-역추정)
- **ℹ️ 안내 & 용어** — 도구 설명 + 보는 법 + **쉬운 용어 사전** + 면책.

### 4.4 결과 내려받기
헤드라인 아래 **⬇️ 결과 내려받기** 펼침에서 요약 표(CSV)/전체 시계열(CSV)/식물 지도(PNG, kaleido 있을 때)를 받습니다.

---

## 5. 전문가 모드 사용법 (English UI)

토글을 켜면 영어 UI로 전환되고, 사이드바가 **1 · Data source / 2 · PFAS compound / 3 · Scenario**로 구성됩니다.

### 5.1 데이터 소스 (1 · Data source) — 노출 Cwᵒ(t)를 어떻게 줄 것인가
| 모드 | Cwᵒ(t) 출처 | Q_TP·M(t) | 언제 쓰나 |
|---|---|---|---|
| **Model (parametric)** | 직접 지정한 상수(또는 flooded shape) | 측정 FAO-56 / ORYZA | 빠른 what-if, 교육 |
| **Custom tables (Cwᵒ + growth)** | 직접 입력한 표 | 표의 성장 + 측정 증산 | 내 성장·오염 시계열이 있을 때 → [6.2](#62-내-데이터-표-성장--토양수) |
| **HYDRUS / CSV drivers** | HYDRUS-1D/Phydrus 결과 CSV | CSV 또는 측정 | 보정된 토양-물-용질 모델이 있을 때 |
| **Run HYDRUS-1D (live)** | 실제 HYDRUS-1D 엔진 실행 | HYDRUS 뿌리흡수 + ORYZA | 앱에서 엔진을 직접 돌리고 싶을 때(빌드 필요) |
| **Soil inventory → pore water** | 총 토양 적재량을 Freundlich로 역산 | 측정 | 총 토양 PFAS만 알 때 |
| **Biomonitoring (measured tissue)** | 측정 토양수 값 | 불필요 | 현장 조직+물 농도가 있을 때 |

### 5.2 화합물 지정 (2 · PFAS compound)
- **Curated congener**: 13종 중 선택(보정된 측정·문헌 파라미터 사용).
- **SMILES (structure)**: 임의의 PFAS 구조를 붙여넣기 → RDKit이 구조 기술자 추출 → (1) 선별 congener와 일치하면 **read-across**, (2) 신규 구조면 **QSPR**(provisional)로 파라미터화. 2-D 구조도 미리보기. (`docs/structure_input.md`)

### 5.3 모델 파라미터(사이드바)
- **E_m [mV]** (root membrane potential): GHK 음이온 배제 레버(벼 −116…−140 mV). 음전위가 셀수록 음이온 흡수가 억제됨.
- **f_xy source**: `recommended`(monotone, 물리적 TSCF) / `W2 fit`(Yamazaki 재현 보정).
- **Biomass driver M(t)**: `ORYZA2000`(메커니즘 탄소수지, 기본) / `growth_rice`(IR72 분배 × 로지스틱).

### 5.4 시나리오 (3 · Scenario) — 모드별 컨트롤
- **Model (parametric)**: `Pore-water Cwᵒ [µg/L]`, `Season length [days]`, `Cwᵒ(t) shape`(constant / flooded(dilution+leaching), per-congener `k_leach` HYDRUS-보정 기본), `Measured forcings` 토글.
- **Custom tables**: 본문 패널에서 표 입력(아래 6.2).
- **HYDRUS/CSV**: 드라이버 CSV 업로드 또는 번들 예시. 컬럼은 [7장](#7-입력-데이터-형식-csv).
- **Run HYDRUS-1D (live)**: `f_oc`, `Flooded until [day]`, `Percolation [cm/day]`. 엔진 미빌드 시 "Build the engine" 버튼.
- **Soil inventory**: `Total soil inventory [µg/kg dry]`, Freundlich `K_F`/`n`/`θ_g`, flooded 여부, `k_leach`.
- **Biomonitoring**: 수동 입력(root/straw/grain conc + Cwᵒ) 또는 CSV.

### 5.5 전문가 탭(9개)
1. **🗺️ Plant & soil map** — 축적 지도(concentration/BAF 토글, day 슬라이더/animate).
2. **📈 Tissue dynamics** — 조직 농도 C_k(t) + **PFAS 질량(burden) C_k·M_k**. B_k/f_xy/L_Ph/κ_d 표시. 낟알 formation-gate 설명.
3. **🟫 Soil & drivers** — 실제 사용된 Cwᵒ(t)·Q_TP(t)·M(t) 드라이버, (소일-인벤토리면) 등온선, 토양 프로파일.
4. **📊 BAF vs observed** — 모델 vs Yamazaki 2023 막대. (선택) **two-pool(seq)** 탐색 모델 오버레이. Yamazaki 조건/매칭 설명.
5. **🔗 Chain-length trends** — 사슬 길이 vs 파라미터(K_PL/K_prot/K_cw/f_xy/B_root/B_grain).
6. **⚖️ Compare congeners** — 사이드바에서 고른 여러 congener의 조직별 BAF 비교.
7. **✅ Tang TF (OOS)** — Tang 2026 per-organ TF(out-of-sample) 검증(PFOA/PFOS/GenX, dw 기준, f_xy refit 비교).
8. **🔎 Inverse (Bayesian)** — 조직 농도 → 노출 Cwᵒ 역추정(식별성 caveat 포함) → [6.1](#61-거꾸로-추정-베이지안-역추정).
9. **ℹ️ About / coupling** — 5개 모드 설명, HYDRUS-1D 입출력 매핑, 용어집.

---

## 6. 공통 기능

> 아래 두 기능은 **일반인·전문가 모드 모두**에서 제공됩니다(라벨만 한/영 다름).

### 6.1 거꾸로 추정 (베이지안 역추정)
"벼에서 PFAS를 측정했는데, 토양수는 얼마나 오염됐었나?"에 **불확실성 범위와 함께** 답합니다.

- **입력**: 뿌리/짚/낟알 측정 농도(µg/kg) 중 하나 이상 + 측정 정밀도(보통 ±~40% / 고정밀 ±~20% / 거침 ±~2×).
- **버튼**: 📐 오염 수준 추정하기(Estimate). (계산은 ODE를 몇 번 풀어 몇 초 걸리므로 버튼식)
- **출력**: 가장 가능성 높은 토양수 수준(µg/L) + 95% 신뢰구간 + 사후분포 곡선 + 입력값 재현 확인.
- **원리**: 뿌리 흡수가 **포화형**(GHK + 운반체)이라 조직 농도는 Cwᵒ의 **비선형** 증가 함수 → 단순 나눗셈이 아닌 진짜 역추정.
  log10(Cwᵒ) 공간에서 **2차 적합 Laplace**로 MAP + 곡률(=사후 폭)을 구합니다(`model_api.estimate_exposure_bayesian`).
- **한계(식별성)**: 여기서는 **노출 수준만** 추정합니다. 조직 데이터만으로는 Q_TP·f_xy(곱), Cwᵒ vs 뿌리흡수 전도도가
  ridge라 분리 불가 → 수송을 절대적으로 고정하려면 독립 측정(수액/공극수 프로브)이 필요합니다.

### 6.2 내 데이터 표 (성장 + 토양수)
직접 만든 두 표를 입력해 모델을 구동합니다(편집 그리드 + CSV 업로드).

- **🌱 성장표**: `day, root, stem, leaf, grain` — 기관별 **신선중(fresh weight)** 시계열.
  - 단위 선택: `g/hill`(기본)·`kg/hill`·`g/m2`·`kg/ha`·`t/ha`.
  - 모델의 M은 **포기당(per-hill) 신선중 질량**이라 단위가 자동 환산됩니다.
- **💧 공극수 표**: `day, Cwo` — 토양수의 **절대** PFAS 농도(µg/L) 시계열.
- **구획 밀도 ρ [kg/L, 신선]**: 기본 root 1.0 · stem 0.30 · leaf 0.30 · grain 1.20(편집 가능).
  - 성장표는 **질량**이고, 모델 수송 ODE는 **질량 기반**(밀도 prefactor 없음)이므로 밀도는 **질량↔부피 환산·보고용**입니다.
  - 화면에 "추정 수확기 기관 부피(= 신선중 ÷ 밀도)"가 표시됩니다(통기조직이 많은 잎/줄기는 밀도 < 1, 낟알은 > 1).
- **부분 입력 허용**: 성장표만 주면 Cwᵒ는 사이드바 값으로, Cwᵒ표만 주면 성장은 선택한 biomass 드라이버로 채워집니다.

### 6.3 다운로드 (CSV / PNG)
- **요약 표(CSV)**: 조직별 model BAF / 최종 농도 / 관측 BAF / (있으면) 측정 BAF.
- **전체 시계열(CSV)**: `t, Cwo, Qtp, conc_*, M_*`.
- **식물 지도(PNG)**: `kaleido`(+Chrome) 설치 시. 없으면 안내 문구로 대체(CSV는 항상 가능).

---

## 7. 입력 데이터 형식 (CSV)

모든 표는 첫 줄이 헤더입니다. **단위 규약**: 시간 `day`, 수용액 농도 `µg/L`, 조직 농도 `µg/kg`, 질량 `kg`, 유량 `L/day`, BAF `L/kg`.

### 7.1 드라이버 CSV (HYDRUS / CSV drivers)
```
t,Cwo,Qtp,M_root,M_stem,M_leaf,M_grain
0,1.0,0.005,0.0001,0.0001,0.0003,0.0001
...
```
- 필수: `t`, `Cwo`. 선택: `Qtp`, `M_root/M_stem/M_leaf/M_grain`(생략 시 측정 forcing/biomass 드라이버로 대체).
- HYDRUS-1D 출력 매핑: `Cwo`←루트존 노드 `Conc`, `Qtp`←`vRoot`/`T_act`, `M_*`←식물 성장 서브모델(HYDRUS 아님).

### 7.2 바이오모니터링 CSV (Biomonitoring)
```
tissue,conc,Cwo
root,0.49,1.0
straw,0.83,
grain,0.46,
```
- `tissue`(root/straw/stem/leaf/grain), `conc`(µg/kg), 선택 `Cwo`(µg/L; 한 번만 있으면 됨). BAF = conc / Cwᵒ.

### 7.3 성장 CSV (내 데이터 표)
```
day,root,stem,leaf,grain
0,0.05,0.02,0.03,0
80,1.0,5.5,4.0,3.0
150,1.1,7.0,4.5,12.0
```
- 기관별 **신선중**(단위는 UI에서 선택). `root` 생략 가능(필요 시 내부 처리).

### 7.4 공극수 Cwᵒ CSV (내 데이터 표)
```
day,Cwo
0,2.0
60,1.5
120,0.5
```

---

## 8. 결과 해석 가이드

- **농도 (µg/kg)**: 조직 1 kg(신선중)당 PFAS 마이크로그램. 일반인 모드 카드/그래프의 기본 표시.
- **축적 배수 BAF (L/kg)**: 조직 농도 ÷ 토양수 농도. "토양수보다 몇 배 진한가". 2면 두 배.
  - 일반인 모드에선 "축적 배수"로만 등장(요약 문장). 전문가 모드 막대/지도에서 직접 봄.
- **축적 지도 색**: 진할수록(뜨거울수록) 그 구획 농도가 높음. 색 스케일은 시점·구획 공통이라 비교 가능.
- **부위 관계**: 보통 **뿌리에 가장 많이** 남고, 짚·낟알로의 이동은 화학물질에 따라 다릅니다(congener-dependent).
  단쇄는 짚이 뿌리보다 높을 수 있고, 장쇄는 뿌리 우세 경향.
- **낟알 형성 게이트**: 낟알은 개화 무렵 형성된 뒤부터 흡수→축적. 그 전 구간은 그리지 않습니다(물리적으로 기관이 없음).
- **불확실성(역추정)**: 사후분포 곡선의 퍼짐이 불확실성. 95% 구간이 넓으면 데이터가 노출을 덜 좁힌 것.
- **dry vs fresh**: 모델은 신선중(fw) 기준. 건조중(dw) 보고 데이터와 비교 시 `C_dw = C_fw/(1−θ_fw)` 환산 필요(전문가).

> [!IMPORTANT]
> Yamazaki 등 관측 막대는 **고정된 실측치**입니다 — 사이드바를 바꿔도 움직이지 않습니다.
> 모델이 그 실험을 재현하도록 보정된 **특정 설정**(Model parametric, Cwᵒ=1, f_xy=W2 fit, E_m −120, season ~120d)에서만
> like-for-like 비교가 됩니다. 다른 설정/모드는 참조 추세이지 보정된 매칭이 아닙니다.

---

## 9. 과학적 배경과 가정

자세한 유도는 `docs/pfas_rice_compartmental_model.{tex,pdf}`, 종합 진입점은 `docs/OVERVIEW_KR.md`.

- **모델 골격**: Trapp/Brunetti **DPU(Dynamic Plant Uptake)**의 **이온성유기화합물(IOC) 확장**, 4구획(뿌리·줄기·잎·낟알) 동역학 ODE.
- **PFAS = 영구 음이온**(매우 낮은 pKa, `f_d≈1`) → 중성 화합물의 Briggs/Kow 분배 코어는 적용 안 됨.
- **뿌리 흡수 `j_R`(하이브리드)**: 이온 전기확산(GHK; 막 안쪽 음전위 ⇒ 음이온 **배제**, `e^N≈107`) **+** 포화 운반체(Michaelis–Menten).
- **내부 이동**: 물관 상향(advection) + 체관(grain은 phloem-fed) + **결합인자 `B_k`**(`θ + f_prot·K_prot + f_PL·K_PL + f_cw·K_cw`, 신선중 basis-A, 밀도 prefactor 없음).
- **뿌리→지상부 적재 `f_xy`(TSCF 유사)**: 음이온은 뿌리에 갇혀 이동이 제한됨. ordering은 congener-dependent.
- **낟알·잎 = 종단 축적체**: 성장 희석이 maturity에서 0 → 정상상태 없음 ⇒ **동역학 모델 필수**(최종 농도 = 시간적분/최종질량).
- **대사 `γ≈0`**(PFAS 난분해), 공기 교환 off.
- **토양 결합(Method A, 단방향)**: HYDRUS-1D/Phydrus → Cwᵒ(t), Q_TP(t)를 받아 식물 ODE는 Python에서 품. HYDRUS는 미수정.

**알려진 한계(정직한 명시)**
- demo의 W2 fit RMSE 0.029는 **포화된 재현**(파라미터=관측 수)이라 예측 검증이 아님. a-priori 예측오차는 log10 RMSE ~0.84–0.95(장쇄 collapse).
- **장쇄(C10–C12) 지상부 과소예측**: 지질 매개 적재/2-pool 분리/장쇄 운반체 용량 등으로 일부 개선되나 PFDoDA 잔차 남음(탐색적).
- `K_cw`(세포벽 분배)는 측정 계수가 문헌에 없어 placeholder. 에터/술폰아미드 Koc는 GAP.
- 자세한 검증·반증(sci-adk)은 `docs/VALIDATION_KR.md`, `sci_adk_review/` 참조.

---

## 10. 용어집

| 한국어(영문) | 쉬운 설명 |
|---|---|
| PFAS | 잘 분해되지 않는 인공 '영원한 화학물질' 무리 |
| 토양수 농도 (Cwᵒ) | 뿌리 주변 토양수에 녹아 있는 PFAS 양 [µg/L] |
| 축적 배수 (BAF) | 조직이 토양수보다 PFAS를 몇 배 진하게 모았는지 [L/kg] |
| 뿌리/짚/낟알 | 식물 부위. 짚 = 줄기+잎, 낟알 = 먹는 현미 |
| 농도 | 조직 1 kg당 PFAS [µg/kg] |
| 화학종(congener) | 특정 PFAS 하나(예: PFOA). 탄소 사슬 길수록 대체로 더 달라붙음 |
| 흡수/이동(translocation) | 뿌리로 들어가 줄기·낟알로 올라가는 과정 |
| 베이지안 추정 | 측정값에서 거꾸로 가장 가능성 높은 원인을 **불확실성 범위**와 함께 추정 |
| f_xy (TSCF) | 뿌리→물관 적재 비율(지상부 이동 용이성) |
| eᴺ (anion exclusion) | GHK 음이온 배제 계수(막 음전위로 음이온 흡수 억제) |
| B_k | 결합인자(조직 농도 = B_k × 국소 자유수 농도) [L/kg fw] |
| TF | transfer factor = C_organ / C_root (Tang 검증) |

---

## 11. 문제 해결 (FAQ)

- **SMILES 모드가 안 떠요 / RDKit 오류**: `pip install rdkit`(또는 `-r requirements-structure.txt`). 없으면 Curated congener 사용.
- **"Run HYDRUS-1D (live)"가 안 보이거나 빌드하라고 해요**: 엔진(FORTRAN, gfortran)이 필요합니다.
  `external/hydrus_source/source`에서 `make` 후 `pip install phydrus`. 웹 세션은 SessionStart 훅이 자동 빌드합니다. 없으면 다른 모드를 쓰세요.
- **PNG 다운로드가 비활성/안내만 나와요**: 선택 패키지 `kaleido`(+Chrome)가 필요합니다. CSV 다운로드는 항상 됩니다.
- **결과가 안 바뀌어요(BAF)**: BAF는 농도/Cwᵒ라 선형구간에서 Cwᵒ에 거의 무관합니다(프리셋을 바꿔도 BAF는 비슷, 농도는 비례).
- **Yamazaki 막대가 안 움직여요**: 정상입니다(고정 실측치). [8장 주의](#8-결과-해석-가이드) 참조.
- **표 입력이 기본값으로 되돌아가요**: 행에 빈 칸/비숫자가 있으면 그 행은 무시됩니다. 최소 2개 완전한 행이 필요합니다.
- **한글이 깨져 보여요(□)**: 브라우저/시스템에 한글 폰트(예: Noto Sans KR)가 있으면 정상 표시됩니다.
- **언어를 바꾸고 싶어요**: 일반인=한국어, 전문가=영어. 사이드바 토글로 전환합니다.

---

## 12. 재현·검증·테스트

```bash
python reproduce_demo.py          # Yamazaki BAF 전체 ODE 재현(W2 fit, log10 RMSE ≈ 0.029)
python reproduce_demo.py --rec    # monotone 물리적 f_xy(a-priori; 단일-짚 불일치 — 정직한 예측오차)
python build_parameters.py        # params/parameters.json 재조립
pip install pytest && pytest      # 전체 테스트(구조/질량보존/QSPR/보정/API/플롯/역추정/표입력 …)
```
- 검증 문서: `docs/VALIDATION_KR.md`, `docs/VALIDATION_TANG2026_*_KR.md`, `validation/*.py`.
- 앱 계산/그림은 head-less 테스트됨: `tests/test_model_api.py`, `tests/test_plots.py`.

---

## 13. 인용·라이선스·참고

- **인용**: PFAS–Rice Compartmental Uptake Model (Trapp/Brunetti DPU의 이온성유기화합물 확장), 2026.
- **소스/문서**: <https://github.com/ccy5123/pfas-rice-model> · `docs/`
- **핵심 참고문헌**:
  - Yamazaki et al. 2023, *Environ. Sci. Technol.* **57**, doi:10.1021/acs.est.2c08767 (관측 BAF)
  - Tang 2026, *J. Hazard. Mater.*, doi:10.1016/j.jhazmat.2025.141017 (out-of-sample TF)
  - Brunetti et al. 2019 *WRR* 10.1029/2019WR025432 · 2021 *ES&T* 10.1021/acs.est.0c07420 · 2022 *JHM* 10.1016/j.jhazmat.2021.127008
  - HYDRUS-1D 4.08 (LGPL-3.0; `external/hydrus_source`)
- **라이선스**: 저장소 `LICENSE` 참조(벤더링된 HYDRUS-1D는 LGPL-3.0).

> [!WARNING]
> 다시 강조: 본 도구는 **연구·교육용**입니다. 규제·식품안전·건강 결정의 근거로 사용하지 마세요.
