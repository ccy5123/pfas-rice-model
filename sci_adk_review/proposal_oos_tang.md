# 연구 배경

메인 sci-adk run은 "Yamazaki 보정 = 예측검증"을 REFUTED로 판정했다(포화 in-sample 적합;
사전적 RMSE 0.84). 그러나 그건 Yamazaki 자기 자신에 대한 평가였다. 예측 타당성의 결정적
시험은 **보정에 쓰지 않은 독립 데이터셋**으로의 표본외(out-of-sample) 예측이다. Tang 2026은
Yamazaki와 다른 토양(담수 paddy pot)·품종(Nipponbare)·투여량 집합으로 측정된 조직별 전이계수
(TF, stalk/leaf/endosperm)를 제공하므로, 이론/QSPR(monotone) transport 파라미터를 Tang에
적합하지 않고 그대로 예측에 쓰면 진짜 표본외 검증이 된다.

# 연구 목표

이론/QSPR 기반 transport 파라미터(Tang에 적합하지 않음)로 구동한 모델은 독립 데이터셋 Tang
2026의 조직별 전이계수를 표본외로 예측한다 (예측이면 표본외 log10 RMSE가 in-sample Tang-재적합
수준에 근접하고, 아니면 크게 초과한다).

# 연구 방법

`model_api.tang_tf_validation`(ORYZA2000 biomass, dry-weight TF)로 PFOA·PFOS·GenX 3종에
대해 f_xy_source="recommended"(이론 monotone, Tang 미적합)의 모델 TF를 Tang 측정 TF와 비교하고,
대조로 Tang-재적합 f_xy(in-sample)도 산출한다(`validation/oos_tang.py`, 0.1 µg/g 저용량).
표본외 log10 RMSE와 in-sample RMSE를 비교한다. 측정(Tang 2026)에 대한 평가다.

# 기대 산출물

독립 데이터셋 표본외 예측 타당성 판정: 표본외 RMSE가 in-sample에 근접하면 예측적, 크게
초과하면(예: 계통적 PFSA 과소·ether 과대) 모델이 한 데이터셋의 파라미터로 다른 데이터셋을
예측하지 못한다는 정직한 한계. sci-adk verdict + 재현 스크립트.
