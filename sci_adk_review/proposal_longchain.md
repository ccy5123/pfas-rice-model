# 연구 배경

PFAS-벼 모델은 단쇄~중쇄 동족체의 조직 BAF는 보정으로 재현하지만, **장쇄(C10–C12)
는 구조적으로 과소예측**한다. ORYZA2000 biomass에서 per-congener 재적합을 해도
PFDoDA(C12)는 f_xy=1·L_Ph=1·kappa_d 한계값에서조차 straw 14.6 vs 49.8, grain 7.4 vs
45.5로 도달하지 못한다. 원인 가설은 **자유음이온(free-anion) 로딩의 구조적 throttle**
이다: 물관·체관 적재 플럭스는 자유 수상농도 Cw=C/B에 비례하는데, 결합 B는 사슬
길이에 따라 ~10^(0.5–0.6/CF₂)로 커지므로 장쇄에서 Cw가 붕괴하여 f_xy를 아무리 키워도
shoot가 굶는다. 문서가 제시한 대안은 **지질-촉진(lipid-facilitated) 결합 로딩**
(`g_xy·C_root`, `g_ph·C_leaf`; B-비의존, 막/지질 결합 풀이 그대로 이동)이다.

# 연구 목표

sci-adk로 장쇄 메커니즘을 정직하게 판정한다. 독립 검증 가능한 가설:

가설 LC1. 자유음이온 로딩은 장쇄(C10–C12) shoot 축적을 구조적으로 과소예측한다 —
f_xy=1 천장에서도 적재 플럭스가 Cw=C/B 붕괴로 묶여 shoot가 관측에 도달하지 못한다.

가설 LC2. B-비의존 지질-촉진 결합 로딩 항(g_xy·C, g_ph·C)은 자유 로딩이 닫지 못하는
장쇄 shoot(straw+grain) 격차의 대부분을 닫는다.

가설 LC3. 이 지질 메커니즘은 root나 단·중쇄 적합을 악화시키지 않고 장쇄 shoot를
재현한다(단일 메커니즘의 깨끗한 해법이다).

# 연구 방법

ORYZA2000 biomass(`oryza_growth`) + 측정 증산에서, 각 동족체를 free-only(monotone
f_xy, lipid off)와 lipid-loading(K_PL-gated g_xy/g_ph)로 풀어 root/straw/grain BAF를
Yamazaki 측정과 비교한다(`validation/longchain_mechanism.py`). 장쇄(nC≥10)
straw+grain의 log10 RMSE를 두 경우로 산출하고, 자유농도 Cw=C/B의 붕괴를 B_root와 함께
제시한다. 재적합의 f_xy=1 천장 미달(`validation/refit_oryza.py`)을 LC1의 보강 증거로
쓴다. 모든 비교는 측정(Yamazaki)에 대한 in-sample 평가임을 명시한다.

# 기대 산출물

장쇄 메커니즘에 대한 판정: (1) 자유 로딩 throttle이 원인인가(LC1), (2) 지질 결합이
이를 닫는가(LC2), (3) 대가 없이 닫는가(LC3). 각 판정의 정직한 범위 — in-sample(보정
데이터), 지질 항의 단일-풀 한계(root 트레이드오프), PFDoDA 잔여 floor — 를 기록하고,
후속(2-풀 자유+결합 분리, PFDoDA 잔여 메커니즘)을 제시한다.
