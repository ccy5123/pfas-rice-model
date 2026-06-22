# 연구 배경

장쇄(C10–C12) PFAS의 잔여 과소예측은 LC1–LC5에서 단계적으로 규명됐다: 원인은
자유음이온 로딩 throttle(LC1), 방향은 지질-촉진 결합 로딩(LC2), 단일-풀은 root를
희생(LC3), 2-풀은 C10–C11만 닫음(LC4), PFDoDA의 uptake 한계는 막 전도도가 아니라
능동 운반체 용량(LC5). 그러나 LC5b가 PFDoDA에 쓴 Vmax ~5×는 한 동족체용 illustrative
값이었다. 이 운반체 강화가 사슬길이의 매끄러운 함수라면 예측적(QSPR) 레버이고,
산만하면 ad-hoc 보정에 불과하다.

# 연구 목표

장쇄 PFAS의 능동-운반체 강화 인자(2-풀에서 측정 root를 재현하는 Vmax 배수)는 사슬길이
n_C의 매끄러운 로그선형 함수여서 QSPR로 예측 가능하다 (선형 적합 R^2 > 0.9이면 예측적,
미만이면 ad-hoc).

# 연구 방법

2-풀 프로토타입(`validation/twopool_longchain.py`, ORYZA2000 biomass)에서 각 장쇄
동족체(PFOA·PFNA·PFDA·PFUnDA·PFDoDA)에 대해 g_xy=g_ph=0(uptake→root 격리)으로 두고
Yamazaki 측정 root를 재현하는 Vmax 배수를 brentq로 적합한다. log10(Vmax 배수)를 사슬길이
n_C 및 막 분배계수 log10 K_PL에 선형회귀하여 기울기와 R^2를 산출한다. 측정(Yamazaki root)에
대한 in-sample 평가임을 명시한다.

# 기대 산출물

운반체 강화 인자의 사슬길이 의존성에 대한 판정: 로그선형(R^2>0.9)이면 장쇄 운반체
강화를 QSPR로 예측 가능(코어 배선 가능)하다는 근거, 아니면 동족체별 ad-hoc 보정에
머문다는 정직한 한계. sci-adk verdict + 재현 스크립트.
