# Tang 2026 독립 검증 (본격 OOS) — 한국어

> **한 줄 요약** — Tang et al. 2026(통제 용량, 150일 전 생장주기 논-벼 시스템; PFOA/PFOS/GenX)으로 모델을
> **out-of-sample(OOS) 검증**했습니다. 모델 빌드엔 이 논문의 **head-group *부호*만** 썼으므로 TF/BCF *크기*는
> 진짜 OOS입니다. 결과: **이삭(grain) 이행은 order-of-magnitude로 맞고(특히 lipid·PFOS: 0.80 vs 0.77)**,
> **PFSA<PFCA 머리기 순서도 재현**됩니다. 그러나 **모델의 줄기(stem) 칸이 거의 비어 Tang의 stalk를 크게
> 빗나가고, 잎(leaf)은 과대**합니다 — 즉 Tang이 모델의 **줄기 구획화 약점**을 정량적으로 짚어줍니다.

![Tang 2026 검증 그림](../validation/figures/tang2026_validation.png)

*그림: (좌·중) PFOA·PFOS의 TF(조직/뿌리) — Tang(검정, 5용량 평균±범위) vs 모델 3변형(monotone/W2/lipid),
분모-무관이라 변환 없이 직접 비교. (우) 모델 PFOA 150일 궤적 — Tang의 정성적 시계열(월별 증가, 늦은 이삭
충전, 1개월 내 상향 이동)과 비교. 생성: `python validation/tang2026_validation.py`.*

---

## 1. 받은 자료와 그 구조 (Tang 2026, JHM 502:141017)
- **본문 PDF + SI(mmc1.docx)** 전체를 받아 추출했습니다.
- 실험: Nipponbare 벼, **150일(5개월) 전 생장주기, 월별 샘플링**, 연속 담수 논토양에 PFOA/PFOS/GenX를
  **5용량(0.1·1·10·50·100 µg/g, 토양 기준)** 으로 처리.
- 조직: **root / stalk / leaf / chaff / endosperm** (5구획).
- 정의(본문 Eq.1–5): **BCF = C_rice/C_soil**, **TF_x = C_x/C_root**.
- 수치 데이터: SI **Table S7(BCF)·S8(TF)** 에 5용량 endpoint 전체 → `docs/literature_db/raw_si/tang2026_doseresponse.csv`로 전사.
- 시간 차원: SI **Table S6** 가 "시간(월) 효과 p<0.001, 농도×시간 상호작용 p<0.001"을 보고(시계열은 존재).
  **단, 월별 원시 농도는 본문 Fig.4a 그림**에만 있어(조밀한 다패널) 신뢰성 있는 수치 디지타이즈는 하지 않았습니다.

## 2. 왜 OOS 검증이 가능한가 + 노출 기준 우회
- 모델 빌드에 **이 논문의 head-group 부호만** 사용(`f_xy(PFSA)=f_xy(PFCA)·e^{−1.1}`) → **TF/BCF 크기는 fit에 안 씀 = OOS.**
- **노출 기준 문제 우회**: Tang 노출은 **토양 슬러리 µg/g**(공극수 µg/L 아님)이지만, **TF=조직/뿌리는 분모(토양 vs 공극수)와 무관**합니다.
  따라서 토양→공극수 변환 없이 **모델 TF를 Tang TF와 직접 비교**할 수 있습니다(가장 깨끗한 시험). 매칭: Tang stalk↔모델 stem,
  leaf↔leaf, endosperm↔grain(모델은 껍질·배유를 한 grain 칸으로 묶음).

## 3. 정량 결과 — TF(조직/뿌리), PFOA·PFOS (5용량 평균)

| PFAS | 조직 | Tang 평균 | monotone | W2 | lipid |
|---|---|---:|---:|---:|---:|
| PFOA | stalk→stem | 1.45 | 0.03 | 0.02 | 0.09 |
| PFOA | leaf→leaf | 1.66 | 5.95 | 3.83 | 13.0 |
| PFOA | endosperm→grain | 0.95 | 0.41 | 0.26 | **2.31** |
| PFOS | stalk→stem | 0.58 | 0.01 | 0.09 | 0.39 |
| PFOS | leaf→leaf | 0.68 | 0.16 | 1.58 | 5.66 |
| PFOS | endosperm→grain | 0.77 | 0.01 | 0.11 | **0.80** |

| log10 RMSE | monotone | W2 | lipid |
|---|---:|---:|---:|
| 3조직 전체 | 1.32 | 0.95 | **0.74** |
| leaf+grain만(깨끗한 매칭) | 1.03 | **0.56** | 0.67 |

**읽는 법:**
- ✅ **이삭(endosperm→grain)**: lipid가 **PFOS 0.80 ≈ Tang 0.77** 로 거의 일치, PFOA는 2.31(약간 과대). monotone/W2는 과소.
  → 장쇄 lipid 메커니즘이 *독립 데이터*에서도 grain 이행을 order-of-magnitude로 맞춤(Kim 2019 신호와 일관).
- ✅ **머리기 순서**: 모든 변형에서 PFSA(PFOS)<PFCA(PFOA) TF — Tang의 PFOS/PFOA stalk TF 0.57/2.22=0.26과 부호 일치.
- ✅ **용량 무관성**: 모델 TF는 (선형이라) 용량 무관, Tang TF도 대체로 용량 무관(특히 PFOS) — 구조적 일치.
- ❌ **줄기(stalk→stem)**: 모델 stem TF 0.01–0.39 vs Tang 0.58–1.45 — **모든 변형이 크게 빗나감.** 모델의 줄기 칸이
  고-증산 통과 칸이라 정상농도가 거의 0 → 문서에 기록된 **줄기-기울기(stem-gradient) 약점**이 그대로 드러남.
- ❌ **잎(leaf→leaf) 과대**: 모델 leaf가 물관 종착 sink라 뿌리 대비 과축적(특히 lipid).

## 4. 분배·BCF (정성) — Fig.3 + Table S7
- Tang Fig.3: PFAS는 **대부분 담수(floodwater)에 존재**, 토양·벼는 소수 → **BCF(=C_rice/C_soil)가 낮음**(PFOA 0.18–0.24, PFOS 0.22–0.30).
- 모델도 **음이온 배제(eᴺ≈107)+낮은 f_xy**로 **뿌리 흡수가 낮고 PFAS가 물에 남는** 방향과 정성적으로 일치.
  (절대 BCF 비교는 토양→공극수 변환(Kd)이 필요해 보류 — §6.)

## 5. 시계열 (정성) — Fig.4a + Table S6
- Tang: 조직 농도가 **월별로 증가**(Table S6 시간효과 p<0.001), **이삭은 늦게 충전**, **1개월 내 상향 이동**.
- 모델 150일 궤적(그림 우): **이삭이 늦게(개화 후) 차오르고** 잎/뿌리가 먼저 축적 → **시간 형태(direction)는 일치**.
- **정량 시계열 비교는 미완**: Tang 월별 원시값이 Fig.4a 그림에만 있어 디지타이즈를 보류했습니다(§6).

## 6. 정직한 결론 + 남은 작업
- **Tang 2026은 (Li처럼 inconclusive가 아니라) 모델을 진짜로 검증·판별**합니다:
  - **맞는 것**: 이삭 이행의 order-of-magnitude(특히 lipid·PFOS), PFSA<PFCA 머리기 순서, TF 용량-무관성, 시계열 형태.
  - **틀리는 것(=다음 개선점)**: 모델 **줄기 칸이 비어 stalk를 크게 과소**, **잎은 과대** → 지상부 구획화(stem↔leaf 분배)가 부정확.
    이는 `src/pfas_rice_plant_module_nstem.py`(다구획 줄기)로 가는 방향을 정량적으로 뒷받침합니다.
- **남은 작업(자료/추가 추출 필요)**:
  1. **절대 BCF 검증** — Tang BCF(토양 기준)를 맞추려면 토양→공극수 변환(Koc/Kd) 필요. 모델 soil 서브모델로 가능(다음 단계).
  2. **정량 시계열** — Fig.4a(월별×조직×용량)에서 0.1용량 root/endosperm 궤적을 digitize하면 시간축 RMSE 산출 가능.
  3. **GenX** — ether-PFCA로 core 12에 없음. 추가하려면 `parameters.json`에 congener 정의 필요.

원하시면 1·2를 이어서 진행합니다(둘 다 추가 자료 없이 이 논문 범위에서 가능; 2는 그림 디지타이즈 작업).

## 7. 재현
```bash
python validation/tang2026_validation.py     # TF OOS 표 + 그림
#   -> validation/figures/tang2026_validation.png
```
데이터: `docs/literature_db/raw_si/tang2026_doseresponse.csv`(SI S7/S8 전사). 관련: `docs/VALIDATION_KR.md`(전체 검증 정리).
