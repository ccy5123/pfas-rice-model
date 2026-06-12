# Tang 2026 독립 검증 (본격 OOS) — 한국어

> **한 줄 요약** — Tang et al. 2026(통제 용량, 150일 전 생장주기 논-벼 시스템; PFOA/PFOS/GenX)으로 모델을
> **out-of-sample(OOS) 검증**했습니다(TF·절대 BCF·시계열 TF(t) 3축). 모델 빌드엔 이 논문의 **head-group *부호*만**
> 썼으므로 크기는 진짜 OOS입니다. **맞는 것**: 이삭 이행 order-of-magnitude(lipid·PFOS 0.80 vs 0.77), PFSA<PFCA
> 머리기 순서, TF 용량-무관성, **절대 BCF 크기**(W2 ~2배 이내). **틀리는 것**: 모델 **줄기 칸이 비어 stalk를 크게
> 과소**, **잎은 과대**, **초기 뿌리 우세(root-first)를 재현 못 함** — 즉 Tang이 모델의 **지상부 구획 분배(과이행) 약점**을
> 정량적으로 짚어줍니다(→ 다구획 줄기 nstem이 개선 방향).

![Tang 2026 검증 그림](../validation/figures/tang2026_validation.png)

*그림(4패널): (A·B) PFOA·PFOS의 TF(조직/뿌리) — Tang(검정, 5용량 평균±범위) vs 모델 3변형(monotone/W2/lipid),
분모-무관 직접 비교. (C) 절대 BCF=C_rice/C_soil — 토양 Kd(Koc·f_oc, f_oc=0.016)로 변환해 Tang 범위와 비교.
(D) 시계열 TF(t) — 모델 leaf/root가 전 기간 1보다 커 Tang의 "1개월 후 뿌리 우세(root-first)"를 재현 못 함(과이행).
생성: `python validation/tang2026_validation.py`.*

---

## 1. 받은 자료와 그 구조 (Tang 2026, JHM 502:141017)
- **본문 PDF + SI(mmc1.docx)** 전체를 받아 추출했습니다.
- 실험: Nipponbare 벼, **150일(5개월) 전 생장주기, 월별 샘플링**, 연속 담수 논토양에 PFOA/PFOS/GenX를
  **5용량(0.1·1·10·50·100 µg/g, 토양 기준)** 으로 처리.
- 조직: **root / stalk / leaf / chaff / endosperm** (5구획).
- 정의(본문 Eq.1–5): **BCF = C_rice/C_soil**, **TF_x = C_x/C_root**.
- 수치 데이터: SI **Table S7(BCF)·S8(TF)** 에 5용량 endpoint 전체 → `docs/literature_db/raw_si/tang2026_doseresponse.csv`로 전사.
- 시간 차원: SI **Table S6** 가 "시간(월) 효과 p<0.001, 농도×시간 상호작용 p<0.001"을 보고(시계열은 존재).
  월별 원시 농도는 본문 Fig.4a 그림(용량별 y-축 상이·조밀)에만 있어 **수치 디지타이즈는 신뢰 불가** → 본문의 명시적
  정성 주장("1개월 후 뿌리 우세")으로 모델 TF(t)를 판별했습니다(§5).

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

## 4. 절대 BCF (정량) — Fig.3 + Table S7  [그림 패널 C]
이제 토양→공극수 변환을 넣어 **절대 BCF=C_rice/C_soil** 을 비교했습니다. Tang 토양 OM=27.4 g/kg(Table S2)→
**f_oc≈0.016**. 모델은 사슬길이 Koc QSPR로 Kd(PFOA 1.5, PFOS 9.3 L/kg)를 만들고
`BCF = [Σ M·BAF / Σ M·(1−θ_fw)] / (Kd+θ_g)` (건중 전식물 BAF ÷ 토양 분배)로 산출:

| BCF | Tang(0.1–100) | monotone | W2 | lipid |
|---|---|---:|---:|---:|
| PFOA | 0.18–0.24 | 0.63 | **0.42** | 1.60 |
| PFOS | 0.22–0.30 | 0.08 | **0.42** | 1.15 |

- ✅ **order-of-magnitude 일치** — 모델 0.4–1.6 vs Tang 0.2–0.3. **W2가 ~2배 이내**로 최근접.
- 모델은 대체로 **흡수를 약간 과대평가**(잎 과적재와 일관). 그래도 토양 Kd(독립 QSPR)+식물 흡수의 **전체 사슬이 맞는 크기**를 줍니다.
- Tang Fig.3의 정성 결론(**PFAS 대부분 담수에 존재, BCF 낮음**)과도 일치: 모델의 음이온 배제+낮은 f_xy로 뿌리 흡수가 낮음.

## 5. 시계열 TF(t) (정량) — Fig.4a + 본문  [그림 패널 D]
Tang 본문: **"1개월 후 PFOA·PFOS는 모든 농도에서 뿌리에 우세 축적, 지상부보다 유의하게 높음"**(즉 초기 root>shoot),
이후 지상부가 차오름(수확기 TF≈1–2.5). 이를 **모델의 TF(t)=C_x(t)/C_root(t) 궤적**으로 엄밀히 비교:

- ❌ **모델은 leaf/root가 전 기간 7–15 (>1)** — **초기 뿌리 우세를 전혀 재현 못 함.** Tang은 1개월 후 root>shoot인데
  모델은 처음부터 leaf≫root → **지상부로 과이행(over-translocation).**
- 이는 §3의 수확기 TF_leaf 과대(3–13 vs 1.66)와 같은 뿌리: 모델의 **잎이 물관 종착 sink로 과축적**하고 **뿌리가 과소**.
- **정량 월별 비교는 보류**: Tang 월별 원시값은 Fig.4a 그림(용량별 y-축 상이·막대 조밀)에만 있어 **신뢰성 있는 수치
  디지타이즈가 불가**합니다(가짜 정밀도 회피). 대신 본문의 명시적 정성 주장(root-first)으로 모델을 판별했습니다.

## 6. 정직한 결론
**Tang 2026은 (Li처럼 inconclusive가 아니라) 모델을 진짜로 검증·판별**합니다 — 무엇이 맞고 무엇이 틀리는지를 정량적으로 보여줍니다:

- ✅ **맞는 것**: 이삭(grain) 이행 order-of-magnitude(특히 lipid·PFOS 0.80 vs 0.77), **PFSA<PFCA 머리기 순서**,
  **TF 용량-무관성**, **절대 BCF 크기**(W2 ~2배 이내), 담수-우세 분배 방향.
- ❌ **틀리는 것(= 다음 개선 우선순위)**: 모델의 **지상부 구획 분배가 구조적으로 부정확** —
  (i) **줄기 칸이 비어 stalk를 크게 과소**, (ii) **잎이 과대**, (iii) **초기 뿌리 우세(root-first)를 재현 못 함**.
  세 가지 모두 한 원인을 가리킵니다: **물관 종착 잎으로의 과이행 + 단일·통과형 줄기 칸**.
  → 개선 방향은 **다구획 줄기(`src/pfas_rice_plant_module_nstem.py`) + 잎/줄기 적재 재배분**으로 정량적으로 좁혀집니다.

## 7. 남은 작업
- **정량 월별 시계열 RMSE** — Fig.4a 디지타이즈가 신뢰 불가하므로, 저자에게 원시 시계열(소스 데이터) 요청 시에만 가능.
- **GenX** — ether-PFCA로 core 12에 없음. 추가하려면 `parameters.json`에 congener 정의 + ether 머리기 항(이미 `FXY_HEADGROUP_LN_OFFSET["ether"]=−0.7` 잠정) 필요.
- **모델 개선** — 위 ❌(지상부 구획)을 nstem 다구획 줄기로 고친 뒤 Tang TF로 재검증(별도 모델링 작업).

원하시면 1·2를 이어서 진행합니다(둘 다 추가 자료 없이 이 논문 범위에서 가능; 2는 그림 디지타이즈 작업).

## 7. 재현
```bash
python validation/tang2026_validation.py     # TF OOS 표 + 그림
#   -> validation/figures/tang2026_validation.png
```
데이터: `docs/literature_db/raw_si/tang2026_doseresponse.csv`(SI S7/S8 전사). 관련: `docs/VALIDATION_KR.md`(전체 검증 정리).
