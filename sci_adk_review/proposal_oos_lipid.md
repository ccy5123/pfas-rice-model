# 연구 배경

표본외 검증(`runs/pfas-rice-oos-tang`)은 **자유음이온 모델**(이론 monotone f_xy, Tang 미적합)이
독립 데이터셋 Tang 2026의 조직별 전이계수를 **예측하지 못함**을 보였다(표본외 log10 RMSE 1.232
vs in-sample 재적합 0.519). 지배적 오차는 **PFOS**(고-K_PL 술폰산)로, 자유음이온 경로가
~40–200× 과소예측한다. 한편 장쇄 sub-investigation(`runs/pfas-rice-longchain`, LC1/LC2 SUPPORTED)은
**B-비의존 지질-촉진 결합 로딩 항**(g_xy·C, g_ph·C, K_PL-gated)이 자유음이온의 in-sample 장쇄
(Yamazaki C10–C12) 붕괴를 닫음을 찾았고, 이는 Chen 2025(막–물 K_MW가 사슬길이로 단조↑ — 지질 풀이
최장쇄·최고흡착종을 운반)로 독립 문헌 corroborate된다. **결정적으로 `LIPID_LOADING` 상수는
Yamazaki(PFDoDA 제외)에 적합된 것이지 Tang에 적합된 것이 아니므로**, 이 메커니즘을 Tang에 적용하면
진짜 표본외 일반화 시험이 된다.

# 연구 목표

Yamazaki에 적합되었으나 Tang에는 적합되지 않은 지질-촉진 로딩 메커니즘은, 자유음이온 기준선 대비
독립 데이터셋 Tang 2026 조직별 TF의 표본외 예측을 개선한다 (메커니즘이 일반화하면 표본외 RMSE가
in-sample 재적합 수준(~0.52)에 근접하고, 특히 PFOS 과소예측이 교정된다).

# 연구 방법

`model_api.tang_tf_validation`(ORYZA2000 biomass, dry-weight TF, f_xy_source="recommended"=이론
monotone·Tang 미적합)으로 PFOA·PFOS·GenX 3종에 대해 **lipid_loading=False(자유음이온 기준선)** 와
**lipid_loading=True(Yamazaki-적합 g_xy/g_ph)** 두 경우의 모델 TF를 Tang 측정 TF(dw, 0.1 µg/g)와
비교한다(`validation/oos_tang_lipid.py`). 두 경우의 표본외 log10 RMSE를 in-sample 재적합 수준과
대조한다. lipid 상수는 Yamazaki 적합값이므로 Tang에 대해 어떤 파라미터도 조정하지 않는다. 측정
(Tang 2026)에 대한 평가다.

# 기대 산출물

지질-촉진 로딩 메커니즘의 표본외 일반화 판정: lipid ON이 표본외 RMSE를 in-sample 수준으로 떨어뜨리고
(특히 PFOS 과소예측 교정) 자유음이온 기준선(1.232)을 크게 개선하면, 메커니즘이 데이터셋을 가로질러
**일반화**(SUPPORTED)하는 것으로, 프로젝트의 첫 강한 교차데이터셋 표본외 예측 성공이다. 잔여(GenX
ether 과대예측, PFOS endosperm 잔여 과소)는 정직히 기록. sci-adk verdict + 재현 스크립트.
