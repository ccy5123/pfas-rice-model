# PFAS–벼 모델 검증 정리 (한국어)

> **한 줄 요약** — 이 모델의 **이행(transport) 파라미터는 Yamazaki 2023에 보정(fit)** 되어 있고, **결합(binding)은 다른 논문의 측정값으로 독립적으로** 받칩니다. 다른 논문에 의한 **독립 검증은 "부분적"** 입니다: Kim 2019(현미 grain)에서 *지질결합(lipid) 메커니즘*이 장쇄 상승을 예측한 것이 가장 강한 신호이고, Li 2025는 교란요인으로 판별 불가(inconclusive)이며, 나머지는 방향·부호 수준의 증거입니다. **"Yamazaki에 fit해서 RMSE 0.029로 재현"은 검증이 아니라 재현**입니다(포화 fit).

![검증 요약 그림](../validation/figures/validation_summary.png)

*그림: (A) Yamazaki 2023 보정(fit) — 예측 vs 관측, 포화 fit이라 1:1에 붙음(재현). (B) Kim 2019 현미 grain — OOS(보정에 안 쓴 데이터)에서 사슬길이별 관측 vs 세 모델, 장쇄 상승을 lipid만 따라감. (C) Li 2025 물-무관 조직비 TF — 교란으로 판별 불가. 생성: `python validation/validation_summary.py`.*

---

## 1. 무엇이 보정(fit)이고 무엇이 독립인가

모델 파라미터는 출처가 분명히 갈립니다. 이 구분이 "검증 여부"를 이해하는 핵심입니다.

| 구분 | 파라미터 | 출처 | 성격 |
|---|---|---|---|
| **보정(fit)** | `f_xy`(뿌리→줄기 적재), `L_Ph`(체관 적재), `κ_d`(막 전도) | **Yamazaki 2023** root/straw/grain BAF | congener마다 **3 파라미터 = 3 관측** (포화) |
| (보정, 탐색적) | `g_xy`, `g_ph`(지질결합 적재) | Yamazaki 2023 (PFDoDA 제외) | 전역 fit, **기본 off**의 opt-in 변형 |
| **독립 측정** | `K_PL` | Chen 2025 (측정 `K_MW`) | PFAS 데이터 아님 |
| **독립 측정** | `K_prot` | Zhou 2025 (투석 `K_prow`) | 〃 |
| **독립 측정** | 토양 `Koc` | Higgins & Luthy / Milinovic QSPR | 〃 |
| **독립 측정** | `f_d`(해리분율) | Goss 2008 (pKa) | 〃 |
| **독립 측정** | 뿌리 막전위 `E_m` | Wang 1994 | 〃 |
| **독립 측정** | 증산 `Q_TP`, 생장 `M` | Kumari 2022 / Nay Htoon 2018, ORYZA IR72 | 작물생리(PFAS 무관) |

즉 **결합·토양·작물생리는 PFAS 흡수 데이터에 맞춘 게 아니라 독립적으로 측정된 값**이고, **이행 파라미터만 Yamazaki에 맞춰져** 있습니다.

---

## 2. (A) Yamazaki 2023 "보정"은 검증이 아니라 재현

`python reproduce_demo.py` 결과 **log10 RMSE = 0.029** — 예측이 관측에 거의 정확히 붙습니다(그림 A의 1:1 선). 그러나 이는 **포화 fit**입니다: congener마다 이행 파라미터 3개를 관측 3개(root/straw/grain)에 맞추므로 **잘 맞는 것이 수학적으로 보장**됩니다.

> ⚠️ 따라서 0.029라는 숫자는 **모델이 데이터를 재현할 수 있다**는 뜻이지, **새 데이터를 예측한다**는 증거가 **아닙니다.** (Yamazaki: 일본 Andosol, congener별 청정수 노출, 12종.)

---

## 3. (B) 검증(OOS): Kim 2019 현미 grain — 가장 깨끗한 독립 시험

**Kim et al. 2019**(한국 논, 공극수↔현미 grain BAF를 짝지어 측정; `10.1016/j.scitotenv.2019.03.240`)은 보정에 쓰지 않은 독립 데이터입니다. 같은 BAF 정의(공극수 기준)라 가장 깨끗한 OOS 시험입니다. `PFOA`는 `L_Ph` 보정에 썼으므로 제외(`*`).

세 이행 변형으로 grain BAF를 **재보정 없이** 전이(transfer)한 결과:

| PFAS | 검출빈도 DF | 관측 | lipid | monotone | W2 |
|---|---:|---:|---:|---:|---:|
| PFHpA | 13% | 0.39 | 0.96 | 0.74 | 0.11 |
| PFOA* | 57% | 4.43 | 0.75 | 0.15 | 0.10 |
| PFNA | 20% | 2.21 | 1.31 | 0.03 | 0.08 |
| PFDA | 7% | 1.56 | 5.18 | 0.04 | 0.50 |
| PFUnDA | 13% | **33.1** | 6.41 | 0.05 | 0.82 |
| PFDoDA | 3% | **35.2** | 6.85 | 0.11 | 3.36 |

| log10 RMSE | lipid | monotone | W2 |
|---|---:|---:|---:|
| PFOA 제외 전체 | **0.55** | 2.04 | 1.11 |
| **신뢰(DF≥15%, PFOA 제외: PFHpA·PFNA)** | **0.23** | 1.91 | 1.43 |

**핵심 발견** — 장쇄(PFUnDA·PFDoDA)에서 grain BAF가 **크게 상승(33·35)** 하는데, 기본 모델(monotone·W2)은 이를 **구조적으로 못 맞춥니다(≈0.05)**. 반면 **지질결합(lipid) 메커니즘만 그 상승 방향을 재현**합니다. 신뢰 congener에서 RMSE 0.23 vs 1.9는 프로젝트의 **첫 진짜 예측 신호**입니다.

**단서(반드시 함께 읽을 것)**
- grain **한 조직만**입니다(root/straw 시계열 없음).
- 신뢰할 만한 점이 사실상 **2개(PFHpA·PFNA)** — 장쇄는 DF 3~13%로 **관측 자체가 불안정**합니다.
- 이 신호를 내는 `lipid_loading`은 **기본 off의 탐색적 변형**입니다. 기본 모델(monotone/W2)이 검증된 것이 아닙니다.

---

## 4. (C) 검증(OOS): Li 2025 — 교란으로 판별 불가

**Li 2025**(톈진 현장)은 BAF가 **물 수질에 반비례**합니다(PFOS는 "나쁜 물"에서 BAF≈250, PFOA는 "좋은 물"에서 ≈2) — 즉 **분모(group water)가 신뢰 불가**라 BAF로는 검증할 수 없습니다. 그래서 **물-무관 조직비(TF = straw/root, grain/root)** 만 사용했는데, 이마저 **뿌리 표면흡착 교란**이 있습니다.

| log10 RMSE (TF) | lipid | monotone | W2 |
|---|---:|---:|---:|
| straw/root | 0.86 | 1.20 | **0.39** |
| grain/root | **0.73** | 1.14 | 1.47 |

세 모델의 우열이 지표마다 뒤바뀌고(straw/root는 W2가, grain/root는 lipid가 최선) 절대 오차도 큽니다 → **Li 2025는 모델을 판별하지 못함(inconclusive).** 예상된 결과입니다(현장·group-water·표면흡착 교란).

---

## 5. 방향·부호 수준의 추가 증거 (크기 아님)

1. **물-무관 cross-field TF 형태** — 여러 필드에 걸쳐 조직비의 *방향*이 monotone `f_xy` 형태와 일관(크기가 아니라 방향만; `docs/DELIVERABLE_GAP_B_fxy.md`).
2. **머리기(head-group) 오프셋 부호** — PFSA가 CF₂-매칭 PFCA보다 이행이 **적다**는 것이 **Tang 2026(논, PFOS/PFOA TF 0.26)과 Yamazaki 2023(0.43) 두 독립 데이터에서 모두 확인**됨 → `f_xy(PFSA) = f_xy(PFCA)·e^{−1.1}` (부호 확정). 에터(GenX)는 `e^{−0.7}` 잠정.

---

## 6. 정직한 결론

- **Yamazaki 보정(RMSE 0.029)은 재현이지 검증이 아님** — 포화 fit이라 잘 맞는 게 보장됨.
- **독립 검증은 부분적**:
  - Kim 2019 grain에서 *lipid 메커니즘*이 **장쇄 상승을 예측**(RMSE 0.23 vs 1.9) — 가장 강한 신호이나 grain-only·신뢰점 2개·opt-in 변형.
  - Li 2025는 **교란으로 판별 불가**.
  - 그 외엔 **방향·부호 수준**의 일관성(cross-field TF, head-group).
- **결론**: 이 모델은 *완전히 예측 검증된 모델이 아니라*, **결합은 독립 측정으로 받치고 이행은 Yamazaki로 보정한 뒤, 제한된 OOS로 방향성을 확인한 단계**입니다.
- **가장 큰 공백**: 깨끗한 **독립 시계열(특히 compartment-resolved root/straw 시간경과)** 의 부재. 이것이 채워져야 진짜 예측 검증이 가능합니다(예: Tang 2026 JHM 시계열 — repo의 열린 과제로 기록됨).

---

## 7. 재현 방법

```bash
pip install -r requirements.txt
python reproduce_demo.py                     # (A) Yamazaki 보정: log10 RMSE 0.029
python validation/oos_crossdataset.py        # (B)(C) Kim 2019 / Li 2025 OOS 수치
python validation/validation_summary.py      # 위 세 패널 그림 생성
#   -> validation/figures/validation_summary.png
```

관련 문서: `docs/DELIVERABLE_GAP_B_fxy.md`(f_xy 형태·cross-field), `docs/fxy_longchain_lipid_exploration.md`(lipid 메커니즘), `docs/data_inventory_and_gaps.md`(사용/필요 데이터), `docs/H8_handoff_S6_final.md`(검증 핸드오프).
