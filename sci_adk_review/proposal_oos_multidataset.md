# 연구 배경

`runs/pfas-rice-oos-lipid`은 지질-촉진 로딩 메커니즘(`LIPID_LOADING` 상수가 Yamazaki에만 적합,
타깃 미적합)이 독립 Tang 2026 조직별 TF의 표본외 예측을 복원함을 보였다(자유음이온 1.232 →
지질 0.516). 그러나 이는 **3 동족체(9 관측)뿐**이다. 이 일반화가 **여러 독립 데이터셋에 걸쳐
견고**한지, 아니면 Tang 한정 우연인지가 핵심이다. 저장소에는 두 번째 독립 표본외 데이터셋
(Kim 2019 현미 곡립 BAF, `validation/oos_crossdataset.py`)과, 현장·교란 데이터셋(Li 2025)이
이미 있다.

# 연구 목표

Yamazaki에만 적합된 지질-촉진 로딩 메커니즘은, **깨끗한 독립 데이터셋**(Tang 2026 조직별 TF,
Kim 2019 곡립 BAF)을 자유음이온(monotone f_xy)·포화 W2 기준선보다 표본외로 더 잘 예측하며,
이 일반화는 여러 데이터셋에 걸쳐 견고하다(현장 교란 데이터셋 Li 2025은 사전등록상 group-water·
표면흡착 교란으로 inconclusive로 취급).

# 연구 방법

세 모델 변형(monotone `f_xy_recommended`=자유음이온, 포화 per-congener `W2fit`, K_PL-gated
`lipid_loading=True`)을 **재적합 없이** 세 독립 데이터셋에 전이해 표본외 log10 RMSE를 비교한다
(`validation/oos_multidataset.py` = `oos_tang_lipid.py`(Tang) + `oos_crossdataset.py`(Kim·Li)).
1차 판정은 **깨끗한 데이터셋(Tang, Kim)**에 둔다; Li는 현장/group-water/표면흡착 교란으로
사전등록상 민감도 점검(1차 시험 아님)이다. 측정(Tang 2026·Kim 2019·Li 2025)에 대한 평가다.

# 기대 산출물

다중 데이터셋 표본외 견고성 판정: 지질이 두 깨끗한 데이터셋(Tang·Kim)에서 기준선을 명확히
이기고 Li가 (예상대로) inconclusive면, 일반화는 Tang 한정 우연이 아니라 **여러 독립
데이터셋에 걸쳐 견고**(SUPPORTED, Li 교란 caveat 명시). 정직한 한계(Li straw/root는 W2가
이김, Kim 장쇄 저-DF 신뢰도)는 기록. sci-adk verdict + 재현 스크립트.
