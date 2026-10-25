# 벼 시료 처리 방법: 레포가 쓰는 논문들은 벼를 어떻게 채취·분리·세척·건조·칭량·도정했나

> 대상 독자: 논에서 벼를 수확해 기관을 나누고, 씻고, 말리고, 무게를 달고, 도정해서 데이터를 만들 실험 파트.
>
> 원칙
> - 모든 내용은 논문 원문(본문과, 있는 경우 SI)을 직접 읽고 적었다. 쪽 번호를 붙였고, 핵심 문장은 원문을 그대로(영문, 25단어 이하) 인용했다.
> - 원문에 없는 항목은 **"보고 안 함"** 으로 적었다. 일반 지식으로 채우거나 추측하지 않았다.
> - 쪽 번호는 학술지에 인쇄된 쪽이다. Kim 2003은 826–837쪽, Yang 2006은 12–19쪽, Sasaki 2005는 1704–1712쪽이다. Gao 2021, Tang 2026, Li 2025는 논문 번호 체계라 PDF 쪽(1–12)이 곧 논문 쪽이다. ORYZA2000 매뉴얼은 인쇄 쪽이고 PDF 쪽은 인쇄 쪽 + 9다. Tang·Li의 SI는 Word 파일이라 고정된 쪽 번호가 없어서 Text/Table 번호로 인용했다.
> - 항목별 전체 표(논문 × 항목, 190행)는 `docs/literature_db/Growth_Sampling_Protocols.csv`에 있다. CSV의 인용문은 모두 PDF 텍스트와 대조해 확인했다(다단 편집 때문에 텍스트로 대조할 수 없었던 1건은 페이지 이미지로 확인).
> - Yamazaki 2023, Kim 2019, Liu 2023은 이전 세션의 감사(`docs/sampling_protocols_audit_KR.md`, `docs/literature_db/Sampling_Protocols.csv`)에서 가져왔다. 이번 세션에서 원문을 다시 읽지는 않았다(PDF가 이번에 제공되지 않음).
> - 원문 PDF는 저작권 때문에 레포에 넣지 않았다.

---

## 0. 요약

포장(노지) 생장 논문 4편(Kim 2003, Sasaki 2005, Yang 2006, Gao 2021)은 큰 틀이 같다. 생육단계(분얼기, 유수형성기, 출수·개화기, 등숙 중기, 성숙기)에 맞춰 날짜를 잡았고, 보호구(border) 포기는 빼고 몇 포기를 통째로 파냈다. 기관은 잎몸 / 줄기(잎집 포함) / 이삭으로 나눴고(Kim·Sasaki·Gao는 뿌리도), 80 °C 오븐(Kim, Yang, Gao)이나 동결건조(Sasaki)로 말린 뒤 **건물중만** 보고했다. 차이가 나는 지점은 이렇다.

- 건조: 80 °C에서 1주(Kim)인지, 항량까지(Yang, Gao)인지, 동결건조(Sasaki, Tang, Li)인지. Yamazaki는 이전 감사 기준으로 식물체 전체를 70 °C에서 1주 말린 **뒤에** 기관을 나눴다.
- 고사엽: 녹색 잎과 고사엽을 따로 칭량했는지. 따로 단 논문은 Kim, Yang, Gao이고, 나누지 않은 논문은 Sasaki, Tang, Li다.
- 뿌리: 정해진 흙 부피를 팠는지. Kim은 포기마다 30 × 17.5 × 15 cm 블록을 팠다. 나머지는 채취 부피를 보고하지 않았고, Yang은 뿌리를 아예 다루지 않았다.
- "낟알"의 실체: 임실 영화의 건물중(Kim), 14% 수분 수량(Yang), "dry grain yield"(Gao), 정의 없는 chaff/endosperm(Tang), 정의 없는 grains/paddy rice(Li). Yamazaki는 왕겨/쌀겨/백미, Kim 2019는 현미다.
- PFAS 시료 세척수: 탈이온수(Tang), 증류수(Li), 수돗물 후 Milli-Q(Yamazaki).

**어느 논문도 보고하지 않은 항목**도 있다.
- 기관별 생중량과 수분 함량. Tang은 조직을 분리 직후 생중량으로 달았지만 값은 내지 않았다.
- 채취 시각, 저울 정밀도, 뿌리 세척에 쓴 체 눈금.
- 철·망간 피막(plaque) 제거.
- **떨어진(탈락) 잎몸을 수거했는지 여부.**
- 제현·도정 조건(생장 논문 전부, Tang, Li), 천립중.

이 공백 가운데 모델에 가장 직접 걸리는 것은 기관별 수분 함량(θ_fw)이다. 레포에는 실측한 계절 값이 하나도 없다.

---

## 1. 논문별 처리 방법

### 1.1 Kim et al. 2003, *Global Change Biology* 9:826–837 (Shizukuishi FACE, Akitakomachi)

| 항목 | 내용 | 쪽 | 원문 인용 |
|---|---|---|---|
| 장소·연도 | 일본 이와테현 Shizukuishi(모리오카 인근) Rice FACE, 1998–2000. FACE 논 4개, 대조(ambient) 논 4개. 대조구는 FACE 링에서 90 m 이상 떨어져 있다 | 826, 827 | "a free air CO2 enrichment (FACE) experiment was conducted at Shizukuishi, Iwate, Japan" |
| 품종 | Akitakomachi | 828 | "seed of rice cv. Akitakomachi was obtained from crops grown locally" |
| 재배 형태 | 노지 담수 논. N 처리구(LN/HN)는 30 cm PVC 판을 흙에 10 cm 박아 분리했다 | 829 | "by a 30 cm PVC board barrier pushed 10 cm into the soil" |
| 재식 | 포기 간격 17.5 cm, 줄 간격 30 cm(19.05 hills m⁻²), 포기당 3본, 손 이앙 | 828–829 | "Hills and rows were 17.5 and 30 cm apart, respectively (equivalent to 19.05 hills m−2)" |
| 육묘·이앙 | 흐르는 물에 침종하고 발아 직후 트레이에 파종했다. 1998·1999년은 비료 함유 암면 매트에 시판 육묘 상토를 덮었고, 2000년은 육묘 상토만 썼다. 트레이는 비닐 챔버에 두었다. 출아 후 14일(1998), 23일(1999), 24일(2000)에 이앙했다(5월 21일/20일/22일) | 828–829 | "hand-transplanted in groups of three (referred to as a hill)" |
| 물관리 | 전 기간 담수. 7월 중순에 5일 중간낙수, 수확 약 10일 전부터 낙수 | 829 | "flooded throughout the season except for a five day summer drainage in mid-July" |
| 채취 시기 | 5회: 분얼 중기, 유수형성기(PI), 개화기, 등숙 중기, 성숙기에 맞췄다. 이앙 후 일수는 Table 1에 있다(MN 1999년: 32/62/88/105/123 DAT) | 829, Table 1 828 | "Sample dates were chosen so as to coincide as close as possible to the mid-tillering, PI, anthesis, mid-ripening and grain maturity growth stages" |
| 생육단계 사건 날짜 | **보고 안 함**. 출수·개화·성숙의 실제 날짜가 없고, 채취일만 단계 이름으로 붙어 있다 | – | – |
| 회당 포기 수·반복 | MN: 1998년은 구당 3곳 × 4포기, 1999·2000년은 구당 2곳 × 6포기. HN: 구당 1곳 × 4포기(1998) 또는 6포기. LN: 1999·2000년 구당 1곳 × 4포기. 블록 4개(블록마다 FACE 1구 + 대조 1구) | 829 | "four hills were sampled from three locations in each plot in 1998 and six hills from two locations in each plot" |
| 포기 선정 | 이웃 포기의 교란을 최소화하도록 골랐다. 처음 두 번의 채취 자리는 보호구에서 가져온 예비 포기로 다시 심었고, 그 포기는 이후 채취하지 않았다. 평균 분얼수 포기를 고르는 식의 대표성 기준은 **보고 안 함** | 829 | "Hills were sampled so as to minimize disturbance of the neighbouring hills." |
| 채취 시각 | **보고 안 함** | – | – |
| 뿌리 채취 | 포기마다 흙 블록(폭 30 cm × 길이 17.5 cm × 깊이 15 cm)을 파냈다. 이 블록 크기는 같은 논문의 줄·포기 간격과 같은 숫자다 | 829 | "a block of soil 30 cm wide, 17.5 cm long and 15 cm deep around each individual hill was dug up" |
| 뿌리 세척 | 1998년은 포기의 1/2, 1999년은 1/3, 2000년은 1/6만 흐르는 물로 씻었다. 나머지 포기의 뿌리량은 씻은 포기의 지상부:뿌리 **건물중** 비로 추정했다. 체 눈금과 피막 처리는 **보고 안 함** | 829 | "soil and other extraneous matter was carefully cleaned from the root system under running water" |
| 기관 분리 | 녹색 잎몸, 고사 잎몸, 줄기(잎집 포함), 이삭(있을 때), 뿌리 | 829 | "The plants were separated into living and senesced leaf laminae, stem (including leaf sheath), panicles (when present) and roots" |
| 고사엽 | 따로 칭량했다(Table 2의 SLdw). **떨어진 잎몸을 모았는지 잃었는지는 보고 안 함** | 831 | "senescent leaf dry weight on a ground area basis" |
| 엽면적 | 포기의 절반에서 녹색 잎면적을 쟀다. 측정 기기·방법은 **보고 안 함** | 829 | "green leaf area was measured for half the hills" |
| 건조·칭량 | 강제송풍 오븐에서 **80 °C, 1주**. 기관별로 따로 말리고 달았다. 생중량과 저울 정밀도는 **보고 안 함** | 829 | "The dry weight of the plant parts was determined separately after drying at 80 °C in a forced air oven for 1 week." |
| 분쇄 | 질소 분석용으로 0.5 mm 체 눈금 분쇄(원문 "milling"은 여기서 분쇄 뜻이고 도정이 아니다). micro-Kjeldahl로 분석 | 829 | "After milling (0.5 mm mesh), the N content of all plant parts was determined separately" |
| 수량 | 성숙기에 LN 18포기(1998·1999)/40포기(2000), MN 54포기(1998)/60포기(1999·2000), HN 10/20/40포기를 수확했다. 임실 영화 수를 세고 낟알 질량을 쟀다. **수량 = 임실 영화의 건물중**(14% 환산 아님). 탈곡, 임실/불임 분리 방법, 제현·도정, 천립중은 **보고 안 함** | 829 | "Grain yield was defined as the dry weight of fertile spikelets." |
| 기관 수분 | **보고 안 함** | – | – |

### 1.2 Sasaki et al. 2005, *Plant Cell Physiol.* 46:1704–1712 (Shizukuishi FACE 1999)

| 항목 | 내용 | 쪽 | 원문 인용 |
|---|---|---|---|
| 장소·연도 | FACE 시설(상세 설명은 Okada et al. 2001 인용). 본문에 지명이 직접 나오지는 않고, 잎 질소를 "같은 장소"(Kim et al. 2003)에서 쟀다고 쓴다. 1999년 | 1710, 1709 | "measured at the same site (Kim et al. 2003)" |
| 품종·재배 | Akitakomachi, 담수 논 | 1710 | "All plants were grown in flooded conditions." |
| 육묘·재식 | 흐르는 물에서 발아시키고, 유근이 나오자마자 비료 함유 암면 매트 트레이에 파종했다. 트레이는 비닐 챔버에 두었다. 출아 후 23일째인 5월 20일에 포기당 3본으로 이앙. 17.5 × 30 cm(19.05 hills m⁻²) | 1710 | "At 23 days after emergence, on 20 May, seedlings were transplanted into the paddy field at three per hill." |
| 생육단계 사건 날짜 | 유수형성기는 이앙 후 약 54일(7월 13일). 출수는 78 DAT(8월 6일, 대조구), 75 DAT(8월 3일, FACE). 성숙일은 **보고 안 함** | 1710 | "The plants reached the panicle initiation stage around 54 days after transplanting (DAT) (13 July), and heading days were 78 DAT" |
| 채취 시기 | 생장 분석은 33, 65, 94, 119 DAT(분얼 중기, 출수 전, 등숙기, 성숙기). LAI는 32, 62, 88, 105, 123 DAT | 1710 | "Plants were sampled for growth analysis at 33 (22 June), 65 (24 July), 94 (22 August), and 119 DAT (16 September)." |
| 회당 포기 수 | 생장 분석은 회당 4포기이고, ¹³C 측정에 쓴 포기와 같다. LAI는 회당 16포기 | 1710 | "Four hills were sampled at each sampling and these were the same plants that we used for measurement of 13C." |
| 포기 선정 | 생장 분석 포기는 ¹³C를 공급한 포기다(투명 PVC 봉지를 120–180분 씌움). 대표성 기준은 **보고 안 함** | 1710 | "Each treated hill was covered with a transparent bag made of 0.10-mm-thick polyvinylchloride film" |
| 채취 시각 | **보고 안 함**. 저녁이라고 적은 것은 군락 호흡 측정뿐이다 | – | – |
| 뿌리 | 흐르는 물로 흙과 이물질을 씻은 뒤 동결건조했다. 채취한 흙 부피·깊이, 체 눈금, 피막 처리는 **보고 안 함** | 1710 | "Roots were carefully cleaned under running water to remove soil and other extraneous matter, and then freeze-dried." |
| 기관 분리 | 잎몸, 줄기(잎집 포함), 뿌리, 이삭. 그림 4 설명에서 줄기는 잎집과 稈(culm)을 합친 것이다. 녹색/고사 잎 구분과 떨어진 잎몸 처리는 **보고 안 함** | 1710, 1706 | "Stem includes leaf sheath and culm." |
| 건조·칭량 | 잎몸·줄기·이삭은 바로 −80 °C에 보관한 뒤 동결건조했다. 건조 후 기관별 건물중을 쟀다. 생중량과 저울 정밀도는 **보고 안 함** | 1710 | "Leaf blades, stems, and ears were immediately stored at –80°C, and then freeze-dried." |
| 엽면적 | 녹색 잎면적. 기기·방법은 **보고 안 함** | 1710 | "Sixteen hills were sampled at each sampling." |
| 낟알·수량 | **보고 안 함**. 수량은 Kim et al. 2001을 인용했다 | – | – |

### 1.3 Yang et al. 2006, *Field Crops Res.* 98:12–19 (Wuxi FACE, Wuxiangjing 14)

| 항목 | 내용 | 쪽 | 원문 인용 |
|---|---|---|---|
| 장소·연도 | 중국 장쑤성 Wuxi, Nianyu Farm(31°35′N, 120°30′E, 해발 5 m), stagnic anthrosol, 벼–밀 윤작. 2001–2003 | 12–13 | "where the soil is classified as stagnic anthrosols (local name, huangni soil)" |
| 품종 | 논문 표기는 "Wuxiangging 14"이고 Hasegawa 2017 SI는 "Wuxiangjing 14"로 쓴다. 조생 만생종 자포니카, 이삭당 영화 평균 155개 | 13 | "Early-maturing late japonica cv. Wuxiangging 14" |
| 재배 형태 | 노지 논. FACE 링 3개(팔각형, 지름 12 m), 대조 링 5개 | 13 | "three FACE rings in octagon, 12 m in diameter (hereinafter called FACE plots), and five similar ambient rings" |
| 재식 | 포기당 3본, 손 이앙. 원문은 17.5 × 30 cm 간격이 24 hills m⁻²에 해당한다고 쓰지만, 17.5 × 30 cm는 계산하면 **19.05 hills m⁻²**다(논문 안의 불일치. Hasegawa 2017 SI도 24로 적었다) | 13 | "Spacing of the hills was 17.5 by 30 cm, (equivalent to 24 hills m−2)." |
| 육묘·이앙 | 5월 18일 파종, 대기 조건에서 육묘, 6월 13일 손 이앙. 육묘 방식은 **보고 안 함** | 13 | "Rice seeds were sown on 18 May. The seedlings grown under ambient air were manually transplanted" |
| 물관리 | 6/13–7/10은 약 5 cm 담수, 7/11–8/4은 여러 번 완전 낙수, 8/5부터 수확 10일 전까지는 간단관개 | 14 | "about 5 cm in depth from 13 June to 10 July, and drained dry several times from 11 July to 4 August" |
| 채취 시기 | 분얼 중기, PI, 출수기(50% 출수), 등숙 중기, 성숙기. 2001년: 28·51 DAT, 출수기, 출수 후 20일, 성숙기. 2002년: 27·47 DAT, 출수기, 출수 후 20일, 성숙기. 2003년: 28·45 DAT, 출수기, 출수 후 30일, 성숙기 | 14 | "coincide as much as possible with the mid-tillering, PI, heading (50% of plant headed), mid-ripening and grain maturity stages" |
| 생육단계 사건 날짜 | **보고 안 함**(출수·성숙 날짜 없음) | – | – |
| 포기 선정·수 | 소구 여러 곳에서 100포기의 줄기 수를 세고, **평균 줄기 수를 가진 5포기**를 골랐다. 채취 뒤 빈자리는 보호구 포기로 다시 심고 이후 채취하지 않았다 | 14 | "the number of stems in 100 hills was counted at different places in each subplot" / "five plants with the mean stem number were selected" |
| 채취 시각 | **보고 안 함** | – | – |
| 뿌리 | **보고 안 함**. 뿌리 채취에 대한 기술이 없고, 표는 지상부만 있다 | – | – |
| 기관 분리 | 녹색 잎, 고사엽, 줄기(잎집 포함), 이삭. 그런데 Table 3은 잎/줄기/이삭(spike) 비율 세 가지만 제시하고, "Leaf"에 고사엽이 들어 있는지는 **보고 안 함**. 떨어진 잎몸 처리도 **보고 안 함** | 14, 16 | "The samples were separated into living and dead leaf tissue, stem (including leaf sheath), and panicle (when applicable)." |
| 엽면적 | 5포기 중 2포기에서 녹색 잎면적을 쟀다. 기기·방법은 **보고 안 함** | 14 | "For two of the five hills green leaf areas were measured." |
| 건조·칭량 | 오븐 **80 °C, 항량까지** 말린 뒤 칭량. 생중량과 저울 정밀도는 **보고 안 함** | 14 | "All the plant parts were oven-dried at 80 °C to constant weight and weighed." |
| 수량 | 소구당 2 m²(보호구 제외)의 식물 전체를 수확했고, **수분 0.14 g H₂O g⁻¹ 생중량(14%)** 으로 환산했다. 탈곡, 임실/불임 분리, 제현·도정, 천립중은 **보고 안 함**(고찰에 임실률·낟알 무게가 나오지만 측정법은 없다) | 14 | "after the grains were adjusted to a moisture content of 0.14 g H2O g−1 fresh weight." |
| 수확지수 정의 | 이 논문의 HI는 성숙기 지상부 건물중 가운데 **이삭(spike) 비율**이다. 이삭축과 왕겨가 포함되므로 낟알 기준 HI가 아니다 | 18 | "the spike DM proportion at grain maturity, which is also called HI" |

### 1.4 Gao et al. 2021, *Agric. Water Manage.* 250:106864 (헤이룽장 Qing'an, 한랭지 자포니카)

| 항목 | 내용 | 쪽 | 원문 인용 |
|---|---|---|---|
| 장소·연도 | Qing'an Irrigation Experimental Station(46°58′06″N, 127°39′37″E), 중국 동북부 쑹넌 평원. 2018, 2019 | 2 | "Field experiments were conducted at Qing'an Irrigation Experimental Station" |
| 품종·재배 | Suijing-18. 라이시미터나 소구가 아닌 정비된 농지에서 관행 담수(basin) 관개 | 3 | "Traditional flood irrigation (i.e. basin irrigation) was applied in our experiment." |
| 재식 | 20 hills m⁻², 포기당 5본. 포기·줄 간격은 **보고 안 함** | 3 | "The planting density was 20 hills per m2 and 5 plants per hill." |
| 육묘·이앙·수확 | 두 해 모두 4월 15일 파종, 5월 16일 이앙. 수확은 9월 18일(2018), 9월 21일(2019). 육묘 방식은 **보고 안 함** | 3 | "The crop was sown on April 15th and transplanted on May 16th in both 2018 and 2019" |
| 생육단계 | Table 2에 연도별 단계 기간이 날짜로 있다(예: 2018년 blooming 7/23–8/3, yellow ripe 8/19–9/18). Table A1의 모델 입력값은 출아 DOY 111, 이앙 136, PI 191, 개화 214, 성숙 261이고 연도는 적혀 있지 않다 | Table 2 p.4, Table A1 p.11 | – |
| 채취 시기·반복 | 시험포 안 3지점에서 약 10일 간격으로 채취했고, 중요한 생육기에는 더 자주 했다. 2018년 9회, 2019년 10회. 반복은 2018년 시험포 3곳 × 3반복, 2019년 시험포 1곳 × 3반복 | 3 | "Crop sampling and observation were taken at three sites within the experimental field about every 10 days during the crop seasons" |
| 회당 포기 수, 포기 선정, 채취 시각 | **보고 안 함** | – | – |
| 기관 분리 | 뿌리, 줄기, 녹색 잎, 고사엽, 이삭. 잎집을 어디에 넣었는지는 **보고 안 함** | 3 | "separated into different components: roots, stems, green leaves, dead leaves and storage organ (i.e. panicles)" |
| 뿌리 | 분리했다고만 쓴다. 채취 부피, 깊이, 세척, 피막 처리는 **보고 안 함**. 뿌리 측정값은 논문에 실리지 않았다(그림 3은 D-AGB, DWST, DWGL, DWSO, LAI만) | 3, Fig. 3 p.8 | – |
| 고사엽 | 고사엽 건물중(DWDL)을 쟀지만 그림에는 없다. D-AGB는 "dry aboveground green biomass"로 정의된다. 떨어진 잎몸 처리는 **보고 안 함** | 3, 8 | "including the dry aboveground green biomass (D-AGB) and the dry biomass of stems, green leaves, dead leaves, and panicles" |
| 엽면적 | 잎마다 길이와 너비를 재서 LAI = 0.75 · ρ · (1/3) · ΣΣ l·b로 계산했다(식 1, Han et al. 2009). 0.75는 형상계수, ρ는 재식밀도 | 3 | "the leaf length and width were measured for each leaf, and the LAI was calculated using Eq. (1)" |
| 건조·칭량 | 오븐 **80 °C, 항량까지**. 생중량과 저울 정밀도는 **보고 안 함** | 3 | "The dry biomass was measured by oven drying at 80 ℃ to the constant mass" |
| 수량 | 수확 때 1 m² × 6반복. "Dry grain yield"라고만 써서 오븐 건조 기준인지 14% 환산인지는 **보고 안 함**. 탈곡, 제현·도정, 천립중도 **보고 안 함** | 3 | "Dry grain yield was collected at 1 m2 with six replicates during harvest." |

### 1.5 Bouman et al. 2001, ORYZA2000 매뉴얼: IR72DS92 시험(IRRI, 1992 건기)

매뉴얼에는 IR72DS92 시험의 **실험 자료 파일**(§7.2)만 실려 있다. 시료를 어떻게 처리했는지는 없다.

| 항목 | 내용 | 쪽 | 원문 인용 |
|---|---|---|---|
| 장소·연도·품종 | IRRI 시험포(14.22N, 121.25E, 23 m), 1992 건기, IR72. 질소 3수준 시험이고(§8.1.1), §7.2의 예시 파일은 225 kg N 처리다(이앙·분얼 중기·PI에 각 60, 개화기 45) | 153, 191–192 | "Experimental site: IRRI farm, 14.22N, 121.25E, 23m" |
| 재식(모델 입력값) | 포기당 5본, 25 hills m⁻², 묘판 1000본 m⁻² | 157 | "NH = 25. ! Number of hills (no m-2)" |
| 육묘·이앙 | 1992년 1월 4일 파종, 1월 16일 이앙(묘판 12일) | 157 | "Transplanting date January 16, 1992; sowing date Jan. 4." |
| 생육단계(관측) | 이앙 DOY 16, PI DOY 58, 개화 DOY 83, 성숙 DOY 114 | 161 | "IDOYPI = 58 ! Day of panicle initiation" |
| 관측 항목 | 녹색 잎(WLVG), 고사엽(WLVD), 줄기(WST), 이삭(WSO, panicle dry weight), 지상부 합계(WAGT), LAI, 잎 질소. 단위는 kg 건물 ha⁻¹이고 관측일은 DOY 4, 16, 34, 58, 68, 83, 97, 114. **뿌리 관측값은 없다** | 161–164 | "Observed values can be given for green leaf weight (WLVG; kg dry matter ha−1)" |
| 잎집 | 매뉴얼은 줄기와 잎집을 함께 다룬다("stem (or sheath)"). IR72DS92 채취에서 잎집을 어디에 넣었는지는 **보고 안 함** | 37 | "stem (or sheath) area index (SAI)" |
| 시료 처리 | 회당 포기 수, 반복, 건조 온도, 뿌리, 탈곡은 모두 **보고 안 함**. 매뉴얼은 시험 정보를 Kropff et al. (1994a)로 미루는데, 이번 세션에서는 그 문헌을 보지 못했다 | 191–192 | "More information on this experiment and on the application of ORYZA2000 is given by" |
| 자료 요건(매뉴얼 지침) | 분배계수를 구하려면 무게 관측을 모두 같은 날 해야 하고, LAI도 가능하면 같은 날 잰다. 고사엽을 재지 않으면 개화 후 감소를 반영할 수 없다. 정조 수량은 14% 수분 기준(WRR14 = WRR/0.86)이다 | 164, 172, 74 | "all the weight observations should be made on the same dates and preferably also for the leaf area index" / "this maximum value is used for partitioning calculations when dead leaves are not measured" |

### 1.6 Hasegawa et al. 2017 SI, Table S1 (두 FACE 시험의 조건 요약)

| 항목 | 내용 | 쪽 |
|---|---|---|
| 일본(Shizukuishi) | 39°38′N/140°57′E. Andosol(모래 43.1%, 실트 31.4%, 점토 25.5%, 용적밀도 0.73 g cm⁻³, 경작층 12.3 cm, 전탄소 82.5 g kg⁻¹). Akitakomachi, 평균 생육기간 146일, 이삭당 영화 81개. 파종은 5/1(1998), 4/23(1999), 4/23(2000), 이앙은 5/21, 5/20, 5/22. **19 hills m⁻²** | SI p.3 |
| 중국(Wuxi) | 31°37′N/120°28′E. Stagnic anthrosol(용적밀도 1.2, 경작층 13.0 cm). Wuxiangjing 14, 158일, 이삭당 영화 155개. 파종 5/18, 이앙 6/13. **24 hills m⁻²** | SI p.3 |
| 시료 처리 | **보고 안 함**(조건 요약표일 뿐이다) | – |

이 표의 파종일로 계산한 Shizukuishi 묘판 기간은 20/27/29일이다. Kim 2003의 이앙 시점("출아 후 14/23/24일")과 모순되지 않는다. 기산점이 파종이냐 출아냐의 차이다.

### 1.7 Tang et al. 2026, *J. Hazard. Mater.* 502:141017 (포트, Nipponbare, PFOA/PFOS/GenX)

| 항목 | 내용 | 쪽 | 원문 인용 |
|---|---|---|---|
| 장소·연도 | 스크린하우스 포트 시험. 토양은 저장대학 농업시험장(항저우) 표토 0–20 cm를 5점 채취법으로 모았다. 스크린하우스 위치와 시험 연도는 **보고 안 함** | 2, 3 | "All experiments were systematically maintained within a screenhouse facility" |
| 포트·토양 | HDPE 버킷(윗지름 18 cm, 높이 16 cm)에 흙 2.5 kg, 묘 20본. 처리마다 3반복, 무작위 배치. 흙은 2 mm 체로 거르고 풍건한 뒤 메탄올 1 L로 PFAS를 넣고 2주 안정화했다. 수분은 건중 기준 40% | 2–3 | "Each treatment consisted of 20 rice seedlings, 2.5 kg of soil, and an individual high-density polyethylene bucket" |
| 육묘 | 종자를 10% H₂O₂로 10분 소독하고 탈이온수로 5번 헹궜다. 30 °C 암조건에서 72시간 발아 → 육묘 트레이 수경 5일 → 1/2 Kimura B 양액(pH 5.6, 72시간마다 교체, 26 °C, RH 80%, 2000 lux) → 4주 묘를 시험에 사용 | 3 | "Four-week-old seedlings exhibiting morphological uniformity were selected for experimental treatments" |
| 물·시비 | 분얼 후기와 등숙기를 빼고 2–3 cm 담수. 관개수는 탈이온수, 무시비 | 3 | "this experiment strictly excluded supplemental fertilization and utilized deionized water for irrigation" |
| 채취 시기 | 노출 150일 뒤 수확, 매월 채취(1–5개월). 날짜와 생육단계는 **보고 안 함** | 3 | "Rice was harvested after 150 days of exposure, with samples collected monthly throughout its entire life cycle." |
| 회당 식물 수 | 버킷당 20본, 처리당 3버킷이라는 것만 있다. **매월 몇 본을 뽑았는지는 보고 안 함** | 3 | "Each treatment including PFAS-free controls was conducted in triplicate." |
| 채취 시각 | **보고 안 함**. 묘 단계의 생리 조작을 "표준화된 채취 시간대"에 맞췄다고만 쓰고, 그 시간대가 언제인지는 없다 | 3 | "synchronized to circadian rhythms through standardized sampling windows" |
| 세척 | 붙은 흙 입자를 탈이온수로 조심스럽게 헹궜다. 성숙기 생장 측정 때도 탈이온수로 흙을 씻었다. **헹굼 횟수, 피막 제거는 보고 안 함**. 종자 소독 뒤의 "5번 헹굼"은 조직 세척이 아니다 | 3 | "Each sample was carefully rinsed with deionized water to remove attached soil particles" |
| 기관 분리 | 뿌리, 줄기(2.6절에서는 stalk), 잎, 왕겨(rice chaff), 배젖(endosperm). 잎집 배정, 고사엽 구분, **chaff와 endosperm의 정의, 탈곡·제현·도정 여부는 보고 안 함** | 3 | "divided into five parts: roots, stems, leaves, rice chaff, and endosperm" |
| 생중량 | 분리 직후 조직 무게를 달았다(생중량). 성숙기에는 식물체와 이삭의 생물량을 따로 쟀다(그림 1, FW). 기관별 생중량 **값은 논문에 없다** | 3, Fig. 1 p.4 | "The plant tissues were immediately weighed" |
| 건조·분쇄 | PE 봉지에 밀봉 → 동결건조기에서 −80 °C 동결건조 → 분말. 분쇄 기기는 **보고 안 함** | 3 | "sealed in polyethylene bags, freeze-dried at − 80 °C in a freeze-dryer, and ground into powders" |
| 분석량 | 동결건조 분말 1 g을 메탄올 5 mL로 초음파 추출했다. SI Text S2에 따르면 추출은 3번 반복 | 3, SI Text S2 | "1 g of freeze-dried powder of rice (root, stalk, leaf, rice chaff, endosperm) and soil samples were extracted by ultrasonication in 5 mL of methanol" |
| 오염 관리 | 전 과정에 폴리프로필렌 재질을 썼고, 모든 용기를 Milli-Q 물로 미리 씻었다 | SI Text S2 | "polypropylene materials were used throughout the experimental process and all containers were pre-washed with Milli-Q water" |
| 보관 | 조직은 PE 봉지에 밀봉해 동결건조했고, 토양은 분석 전 −20 °C에 보관했다. 분말 보관 조건은 **보고 안 함** | 2, 3 | "All samples were stored at –20 ◦ C refrigerator before the analyses." |
| 농도 기준 | 건중 기준 | 7 | "The Bioconcentrations (based on dry weight)" |
| 수량·수분 기준 | **보고 안 함**. SI Text S3의 "백미(GB/T 1354-2018) 수분 15.5% 이하"는 섭취 안전 기준을 계산하는 데 쓴 값이고, 분석 시료에 대한 설명이 아니다 | SI Text S3 | "that the moisture content of japonica rice shall not exceed 15.5%" |
| chaff 관련 언급 | 결론부에서 왕겨·배젖 농축을 "rice bran feed" 같은 부산물을 통한 2차 오염과 엮어 말한다. **정의는 아니다** | 9 | "suggesting vigilance against secondary pollution from agricultural and sideline products (such as rice bran feed)" |

### 1.8 Li et al. 2025, *J. Hazard. Mater.* 492:138256 (톈진 논, Jinyu Japonica No.22)

| 항목 | 내용 | 쪽 | 원문 인용 |
|---|---|---|---|
| 장소·연도·품종 | 톈진 교외 논(38°34′N, 116°43′E), 1년 1작, 100 m × 100 m가 넘는 필지. 2023년 7–10월. Jinyu Japonica NO.22 | 2 | "Expansive paddy fields that cultivate one crop of rice annually, which exceed 100 m × 100 m, were chosen for this study." |
| 재식·이앙·물관리·시비 | **보고 안 함** | – | – |
| 채취 시기 | 분얼기, 출수기, 등숙기(filling), 유숙기(milk-ripening), 완숙기(ripening). **날짜는 보고 안 함** | 2 | "tillering, heading, filling, milk-ripening, and ripening" |
| 반복 | 3반복. 반복마다 20 m × 20 m 구역에서 부분시료 5개를 모았다. 부분시료당 포기 수와 포기 선정 기준은 **보고 안 함** | 2 | "with each replicate comprising five sub-samples sourced from an area spanning 20 m × 20 m" |
| 운반 | PP 봉지에 밀봉하고 아이스백에 넣어 바로 실험실로 옮겼다 | 2 | "All the samples were sealed in polypropylene (PP) bags, stored in ice bags, and transported to the laboratory immediately." |
| 세척 | 채취 당일 증류수로 흙 입자를 씻었다. **헹굼 횟수, 뿌리 채취 부피, 피막 제거는 보고 안 함** | 2 | "On the same day, these plant samples were carefully washed with distilled water to remove soil particles." |
| 기관 분리 | 방법 절은 "뿌리, 잎(줄기 포함), 낟알"이라고 쓰는데, Table 1과 SI Table S6은 뿌리/줄기/잎/낟알을 **따로** 보고한다(논문 안의 불일치). 낟알은 출수기부터 있다. 잎집, 고사엽, **낟알이 왕겨를 포함하는지와 제현·도정 여부는 보고 안 함** | 2, Table 1 p.4 | "Plant samples were divided into roots, leaves (including stems), and grains." |
| 건조·분쇄 | 동결건조(온도·시간 **보고 안 함**) 후 전동 분쇄기로 분말. 조직 1 g을 아세토니트릴 6 mL로 초음파 추출, 3번 반복(SI) | 2, SI | "Both soil and plant samples were freeze-dried and ground into powder using an electric grinder" |
| 낟알 단백질 | 자동 켈달 질소분석기(KDY-9830)로 쟀다. 7.84–14.55% | 6 | "determined by an automatic Kjeldahl nitrogen meter (KDY-9830" |
| 생물량·생중량·수분 | **보고 안 함** | – | – |
| 농도 기준 | 건중(ng/g dw) | 3 | "the root, stem, leaf, and grain (ng/g dw)" |

### 1.9 이전 감사에서 가져온 PFAS 데이터셋 (이번 세션에서 원문 재확인 안 함)

출처는 `docs/literature_db/Sampling_Protocols.csv`다. 쪽 번호는 그 CSV에 적힌 그대로다.

- **Yamazaki 2023**(ES&T, NARO 쓰쿠바 Andosol 노지 논, Kasalath + Koshihikari, 2020년)
  - 뿌리는 수돗물로 충분히 씻고 Milli-Q로 헹궜다. 피막 제거와 탈착 단계는 없었다.
  - **식물체 전체를 70 °C에서 1주 말린 뒤** 기관을 나눴다.
  - 지상부는 잎을 포함한 채 높이 4구간(0–20/20–40/40–60/>60 cm)으로 잘랐다. 잎집과 잎몸은 나누지 않았다.
  - 낟알은 제현기로 왕겨와 현미를 나누고, 현미를 도정해 쌀겨와 백미로 나눴다(백미:현미 86–90%). **현미는 분석하지 않았다.**
  - 기관 질량은 보고되지 않았다. (SI p.S2–S3)
- **Kim 2019**(Sci. Total Environ., 한국 논 30필지, 2016년 수확기)
  - 필지당 1주를 채취했다.
  - 낟알은 −30 °C에서 동결건조하고 도정기로 왕겨를 벗긴 뒤 2 g을 분쇄했다(Fritsch). 즉 **현미**다.
  - 뿌리와 지상부는 분석하지 않았다. (본문 p.715–716)
- **Liu 2023**(수경, 유묘)
  - 뿌리를 탈이온수로 3번 씻고 아세토니트릴로 빠르게 헹궜다. 표면 탈착 단계가 있는 유일한 자료다.
  - 기관은 가위로 잘랐고 생중량 기준이다. (본문 p.3)

---

## 2. 비교표

"–"는 해당 없음, "NR"은 보고 안 함이다. †는 이전 감사에서 인용한 값(이번 세션 재확인 안 함)이다.

| 논문 | 재배·재식 | 채취 단위·반복 | 뿌리 채취·세척 | 기관 분리(잎집 / 고사엽) | 건조 | 생중량 | 낟알·수량 | 엽면적 |
|---|---|---|---|---|---|---|---|---|
| Kim 2003 | 논, 19.05 hills m⁻², 3본 | 4–6포기 × 1–3곳/구, 4블록, 보호구 보식 | 30×17.5×15 cm 블록. 1/2–1/6 포기만 흐르는 물로 세척, 나머지는 지상부:뿌리 건물중 비로 추정. 체 NR | 잎집→줄기 / 녹색·고사 잎몸 분리, 탈락 잎몸 NR | 80 °C 강제송풍 오븐, 1주 | NR | 임실 영화 건물중(10–60포기). 탈곡·도정 NR | 포기 절반, 방법 NR |
| Sasaki 2005 | 논, 19.05 hills m⁻², 3본 | 4포기(¹³C 공급 포기), LAI 16포기 | 흐르는 물 세척 후 동결건조. 부피·체 NR | 잎집→줄기(줄기 = 잎집 + 稈) / 구분 NR | −80 °C 보관 후 동결건조 | NR | NR | 방법 NR |
| Yang 2006 | 논, 24 hills m⁻²로 기재(산술 19.05), 3본 | 평균 줄기 수 5포기/소구, 보호구 보식 | 뿌리 NR | 잎집→줄기 / 녹색·고사 분리(표에서 합쳤는지 NR) | 80 °C 오븐, 항량 | NR | 2 m², 14% 수분 생중량 기준. HI = 이삭 비율 | 2/5포기, 방법 NR |
| Gao 2021 | 논, 20 hills m⁻², 5본 | 3지점, 약 10일 간격, 9–10회. 포기 수 NR | 분리했다고만 씀. 방법·값 NR | 잎집 NR / 녹색·고사 분리(고사엽 값은 미제시) | 80 °C 오븐, 항량 | NR | 1 m² × 6, "dry grain yield". 수분 기준 NR | 잎 길이 × 너비 × 0.75 |
| IR72DS92(매뉴얼) | 논, 25 hills m⁻², 5본(입력값) | NR | 뿌리 관측 없음 | NR / WLVG·WLVD 분리 | NR | NR | WSO = 이삭. 모델 WRR14는 14% | NR |
| Tang 2026 | 포트(2.5 kg, 20본), 3반복 | 매월. 회당 본수 NR | 탈이온수 헹굼. 횟수·피막 NR | 잎집 NR / NR. chaff·endosperm 정의 NR | −80 °C 동결건조 → 분말 | 분리 직후 칭량(값 미공개) | 정의·도정 NR | – |
| Li 2025 | 농가 논. 재식 NR | 3반복 × 부분시료 5(20×20 m) | 증류수 세척. 부피·피막 NR | 방법 "잎(줄기 포함)"이지만 표는 분리 / NR | 동결건조 → 전동 분쇄 | NR | grains·paddy rice라는 용어의 정의 NR | – |
| Yamazaki 2023† | Andosol 노지 논 | 혼합 1점/조직/품종 | 수돗물 → Milli-Q, 피막 제거 없음 | 높이 구간(잎 포함), 잎집·잎몸 미분리 | 70 °C 1주(전체 건조 후 분리) | NR | 왕겨/쌀겨/백미(현미 미측정) | – |
| Kim 2019† | 농가 논 30필지 | 필지당 1주 | 미분석 | 미분석 | 낟알 −30 °C 동결건조 | NR | 현미(도정기로 왕겨 제거) | – |
| Liu 2023† | 수경 유묘 | n = 3 | 탈이온수 3회 + 아세토니트릴 헹굼 | 가위로 분리 | 생중량 | 생중량 기준 | – | – |

---

## 3. 실험 파트를 위한 권장 체크리스트

각 항목에는 (a) 그렇게 한 논문, (b) 그 항목이 메우는 레포의 공백을 붙였다. 원문에 근거가 없는 수치는 넣지 않았다. 논문끼리 다른 곳은 선택지와 출처를 나란히 적었다.

### A. 채취 설계

1. **채취일과 생육단계 사건일을 따로 기록한다.** 이앙일, 유수형성기, 출수기(Yang의 정의: 50% 출수), 개화기, 성숙기를 채취일과 별도로 적는다.
   - (a) 사건 날짜를 보고한 것은 Sasaki 2005(PI 약 54 DAT, 출수 78 DAT), Gao 2021(Table 2의 단계 기간), IR72DS92(이앙·PI·개화·성숙 DOY)뿐이다. Kim 2003, Yang 2006, Tang, Li는 채취일에 단계 이름만 붙였다.
   - (b) 레포는 Kim 1999 채취일의 단계 이름이 실제 사건보다 8–10일 늦다는 것을 Sasaki 2005로 뒤늦게 알았다(CLAUDE.md §6 option C 재분석). Yang과 Li에는 날짜가 아예 없어서 모델의 DVS에 시료를 끼워 맞춰야 했다. ORYZA의 발육속도 보정(DRATES)도 이 날짜들을 입력으로 요구한다(매뉴얼 p.161, 164).
2. **포기 선정 기준을 정하고 기록한다.**
   - (a) 논문마다 방식이 다르다. Yang 2006은 소구 여러 곳에서 100포기의 줄기 수를 세어 평균 줄기 수를 가진 5포기를 골랐다. Kim 2003은 구당 1–3곳에서 4–6포기를 골랐고 이웃 포기의 교란을 최소화했다. Sasaki 2005는 4포기다. 보호구 포기를 빼고 채취 자리를 보호구 포기로 보식한 것은 Kim과 Yang 둘 다다.
   - (b) 회당 포기 수가 보고되지 않은 자료(Gao, Tang, Li)는 분산을 해석할 수 없다.
3. **포기 간격, 줄 간격, 포기당 본수를 모두 적고 계산한 hills m⁻²를 함께 적는다.**
   - (a) Yang 2006은 17.5 × 30 cm 간격을 24 hills m⁻²라고 적었는데 계산하면 19.05다. 간격과 밀도를 둘 다 적어야 이런 불일치를 잡을 수 있다.
   - (b) 포기당 무게를 m² 단위로 바꾸려면 필요하다. 모델 드라이버는 kg/hill, ORYZA는 kg ha⁻¹이다.
4. **채취 시각을 기록한다.**
   - (a) 어느 논문도 보고하지 않았다.
   - (b) 아래 7번의 생중량을 θ_fw 계산에 쓸 것이므로, 같은 조건에서 비교할 수 있게 시각은 남겨 둔다. 시간대를 어떻게 정할지는 원문 근거가 없으니 팀이 정한다.
5. **모든 기관의 무게와 LAI를 같은 날, 같은 포기에서 잰다.**
   - (a) ORYZA 매뉴얼은 분배계수를 추정할 때 이를 요구한다(p.164). Kim, Yang, Sasaki는 같은 포기의 일부에서 엽면적을 쟀다.
   - (b) 기관 분배(잎/줄기/이삭)는 레포가 생장 드라이버를 검증하는 핵심 지표다(Wuxi 분배 검증).

### B. 기관 분리

6. **잎몸을 녹색과 고사로 나눠 칭량한다.**
   - (a) Kim 2003, Yang 2006, Gao 2021, IR72DS92(WLVG/WLVD)는 나눴고, Sasaki, Tang, Li는 나누지 않았다.
   - (b) 레포의 Shizukuishi 비교에서 ORYZA2000 포트가 유일하게 지는 지표가 고사엽 비율이다. 매뉴얼도 고사엽을 재지 않으면 개화 후 감소를 반영할 수 없다고 쓴다(p.172). 모델의 잎 노화 손실 항도 이 비율에 기댄다.
7. **떨어진 잎몸을 가능하면 수거해서 별도 항목("탈락 잎몸")으로 칭량한다.**
   - (a) **어느 논문도 하지 않았다.** 탈락 잎몸을 어떻게 처리했는지조차 보고한 논문이 없다.
   - (b) 레포는 Kim 2003의 고사엽 값을 "탈락분을 잃었으므로 하한"으로 해석하지만(CLAUDE.md §6), 원문에는 그런 말이 없다. 따로 수거한 값이 있어야 이 해석을 검증할 수 있다. 수거 방법(그물, 수면 수거 등)과 면적은 원문 근거가 없으므로 팀이 정해서 기록한다.
8. **잎집은 줄기에 넣고, 그렇게 했다고 기록한다.** 여력이 있으면 잎집과 稈을 따로 달아 합산값도 함께 남긴다.
   - (a) Kim 2003, Yang 2006, Sasaki 2005("Stem includes leaf sheath and culm")가 그렇게 했고, ORYZA도 줄기와 잎집을 함께 다룬다(p.37). Gao, Tang, Li, IR72DS92는 잎집을 어디에 넣었는지 보고하지 않았다.
   - (b) 모델의 줄기 구획과 짚(straw = 줄기 + 잎)을 이 정의에 맞출 수 있다.
9. **이삭을 단계적으로 나눈다:** 이삭 → 이삭축(rachis) / 정조(rough rice) → 임실·불임 → 왕겨 / 현미 → 도정 → 백미 / 쌀겨. 단계마다 무게를 단다.
   - (a) 논문마다 다른 부위를 냈다. Kim 2003은 임실 영화 건물중, Yang 2006은 14% 수분 수량, Gao 2021은 "dry grain yield"(정의 없음)다. Yamazaki†는 왕겨/쌀겨/백미로 나눴고(백미:현미 86–90%), Kim 2019†는 현미다. Tang의 chaff/endosperm과 Li의 grains/paddy rice는 정의가 없다. 이삭축을 따로 단 논문은 없다.
   - (b) 레포의 "낟알"은 데이터셋마다 다른 부위다(감사 §2.1). 모델의 낟알 결합계수는 현미 조성을 쓰지만 ORYZA의 WSO는 왕겨를 포함한 이삭 전체다. 현미 질량과 백미 질량이 둘 다 있어야 서로 비교할 수 있다.
10. **기관별 무게를 남겨서 짚 농도를 질량 가중으로 계산할 수 있게 한다.**
    - (a) Yamazaki†는 높이 구간의 무게를 내지 않았다.
    - (b) 레포는 그래서 짚을 4구간 농도의 **단순 평균**으로 계산하고 있다(감사 CSV).

### C. 무게와 건조

11. **기관마다, 채취 때마다 생중량과 건물중을 둘 다 잰다. 생중량은 분리 직후에 단다.**
    - (a) 생장 논문 4편과 IR72DS92는 건물중만 냈다. Tang은 분리 직후 생중량을 달았지만("The plant tissues were immediately weighed") 기관별 값은 공개하지 않았다. Yamazaki†처럼 식물체 전체를 말린 **뒤에** 나누면 기관별 생중량은 원리적으로 얻을 수 없다.
    - (b) **θ_fw 공백.** 모델의 기관 질량 M과 결합계수 B_k는 생중량 기준인데, 레포의 θ_fw(뿌리 0.90, 줄기 0.83, 잎 0.78, 낟알 0.14)는 출처가 확인되지 않은 고정값이다. 벼 기관 수분의 계절 변화를 실측한 자료는 레포에 하나도 없다. 레포의 민감도 분석에서는 뿌리 θ가 가장 크게 영향을 준다.
    - 세척한 시료는 표면 물기를 어떻게 처리했는지(닦았는지, 어떻게)도 함께 적는다. 이번에 읽은 벼 논문 중 이를 보고한 것은 없다. 이전 감사에 따르면 Briggs 1982는 뿌리 물기를 blot했고, Hwang 2017은 종이타월로 가볍게 닦았다.
12. **건조 방식을 하나로 정해서 기록한다.** 논문마다 다르므로 선택지와 출처를 같이 둔다.
    - 오븐 80 °C, 1주(강제송풍): Kim 2003
    - 오븐 80 °C, 항량까지: Yang 2006, Gao 2021
    - 동결건조: Sasaki 2005(−80 °C 보관 후), Tang 2026(동결건조기 −80 °C), Li 2025(조건 보고 안 함), Kim 2019†(낟알, −30 °C)
    - 오븐 70 °C, 1주(전체 건조 후 분리): Yamazaki†

    PFAS 농도 분석용 시료를 동결건조한 것은 Tang, Li, Kim 2019다. 항량을 확인했는지도 적는다. "항량까지"는 Yang과 Gao가 썼다.
13. **저울 정밀도를 기록한다.** 어느 논문도 보고하지 않았다.

### D. 뿌리

14. **포기마다 정해진 흙 부피를 파서 뿌리를 채취하고, 그 크기와 깊이를 기록한다.**
    - (a) Kim 2003은 폭 30 cm × 길이 17.5 cm × 깊이 15 cm 블록을 팠다(줄·포기 간격과 같은 숫자). 다른 논문은 채취 부피를 보고하지 않았고, Yang 2006은 뿌리를 다루지 않았다.
    - (b) 레포의 논 뿌리 질량 기준값은 Kim 자료에서 나왔다. 블록보다 깊은 곳의 뿌리가 빠지는지는 원문이 말하지 않으므로(레포는 "하한"으로 해석) 깊이를 꼭 남긴다.
15. **뿌리 세척 방식(물 종류, 체 눈금, 횟수)을 기록한다.**
    - (a) 흐르는 물: Kim 2003, Sasaki 2005. 탈이온수: Tang. 증류수: Li. 수돗물 후 Milli-Q: Yamazaki†. 탈이온수 3번 + 아세토니트릴: Liu†. 체 눈금은 **아무도 보고하지 않았고**, 헹굼 횟수를 쓴 것은 Liu†뿐이다.
    - (b) 세척 강도에 따라 미세근 손실과 표면 흙 잔류가 달라지고, 뿌리 건물중과 PFAS 농도가 같이 움직인다.
16. **(PFAS) 철·망간 피막 제거 단계를 추가하되, 제거 전 시료도 같이 분석한다.**
    - (a) 이번에 본 논문과 이전 감사의 데이터셋 **어디에도** 피막 제거는 없다. 탈착 단계가 있는 것은 Liu 2023†의 아세토니트릴 헹굼뿐이다.
    - (b) 레포는 뿌리 농도에 표면·피막의 PFAS가 섞여 있다는 문제를 안고 있다. two-pool 뿌리 연구(`docs/twopool_kseq_mechanism.md`)는 세포벽 대 피막 흡착을 구분하는 실험을 승격 조건으로 꼽는다.
    - 피막 제거 시약과 조건은 이 문헌들에 근거가 없으므로 팀이 따로 정한다. 같은 뿌리를 둘로 나눠 제거 전과 제거 후를 모두 분석하면 기존 자료(전부 미제거)와도 비교할 수 있다.
17. **지상부:뿌리 비는 건물중 기준으로 적는다.**
    - (a) Kim 2003이 건물중 비를 썼다.
    - (b) 생중량 비를 쓴 자료(Nada 2018, 감사 §1.4)와 섞이지 않게 하기 위해서다.

### E. 엽면적

18. **녹색 잎만 재고, 측정 방법을 기록한다.**
    - (a) Gao 2021은 잎 길이 × 너비 × 0.75(형상계수)로 계산했다. Kim, Yang, Sasaki는 녹색 잎면적을 쟀지만 기기를 밝히지 않았다.
    - (b) 생장 드라이버의 LAI·SLA를 검증할 때 쓴다.

### F. PFAS 분석용 처리

19. **세척수 종류와 헹굼 횟수를 기록한다.**
    - (a) 탈이온수: Tang. 증류수: Li. 수돗물 후 Milli-Q: Yamazaki†. 조직 헹굼 횟수를 보고한 PFAS–벼 논문은 없다.
    - (b) 데이터셋끼리 뿌리·지상부 표면 오염 수준을 비교할 수 있다.
20. **채취 기구와 용기를 정해 둔다.**
    - (a) Tang은 전 과정에 폴리프로필렌 재질을 쓰고 용기를 Milli-Q 물로 미리 씻었다(SI Text S2). Li는 PP 봉지에 넣어 아이스백으로 바로 운반했다. Tang은 조직을 PE 봉지에 밀봉해 동결건조했다.
21. **건조·분쇄·분석량을 정한다.**
    - (a) 동결건조 후 분말화한 것은 Tang과 Li(Li는 전동 분쇄기), 2 g 분쇄는 Kim 2019†(Fritsch)다. 추출에 쓴 양은 Tang이 분말 1 g + 메탄올 5 mL, Li가 1 g + 아세토니트릴 6 mL(SI)다.
22. **보관 조건을 기록한다.**
    - (a) Sasaki는 분리 즉시 −80 °C에 넣었다가 동결건조했고, Tang은 토양을 −20 °C에 두었다. 분말 보관 조건은 Tang과 Li 모두 보고하지 않았다.
23. **농도는 건중으로 보고하되, 같은 시료의 생중량과 건물중을 남긴다.**
    - (a) Tang과 Li는 건중 기준이다.
    - (b) 모델은 생중량 기준(B_k, M)이므로 건중 → 생중량 환산에는 같은 시료의 θ_fw가 필요하다(11번).

### G. 수량

24. **수량은 수확 면적이나 포기 수, 원시 무게와 함께 오븐 건조 기준과 14% 수분 기준을 둘 다 남긴다.**
    - (a) Kim 2003은 임실 영화 건물중(10–60포기), Yang 2006은 2 m²·14% 생중량 기준, Gao 2021은 1 m² × 6에 수분 기준 미상이다. ORYZA의 WRR14는 낟알 건물중 / 0.86이다(매뉴얼 p.74).
    - (b) 수확지수의 정의도 논문마다 다르다. Yang의 HI는 이삭 비율이다. 그래서 낟알 질량과 이삭 질량을 둘 다 남긴다.
25. **천립중과 임실률은 선택 항목이다.**
    - (a) 이번에 본 논문 중 측정법까지 보고한 것은 없다. Kim 2003은 임실 영화 수를 세었다.

---

## 4. 이전 감사·레포 서술과의 차이 (이번에 원문을 읽고 확인한 것)

> **조치(2026-10-25)**: 1–6번은 레포에서 바로잡았습니다 — 감사 표·`Sampling_Protocols.csv`·`raw_si` 머리말·`model_api` 주석의
> 서술을 원문 기준으로 고치고 추론은 추론으로 표시했습니다. 4번(고사엽 하한)은 FACE 비교 스크립트가 두 해석(하한 / 전부
> 부착)의 판정을 모두 출력하도록 바꿨고(품종 생육 보정 시 이식 ORYZA2000은 하한 해석에서 우세, 전부 부착 해석에서 구별 불가),
> 5번(Yang 밀도)은 Wuxi 비교를 19.05 hills m⁻²로도 돌려 판정이 그대로임을 확인했습니다.

1. **Tang 반복 수.** 감사 표 §2.3은 "3반복 × 20주"라고 적었다. 원문은 버킷당 묘 20본과 처리당 3버킷만 보고하고, **매월 채취 때 몇 본을 뽑았는지는 보고하지 않는다.** 20본은 심은 수이지 분석한 시료의 수가 아닐 수 있다.
2. **Tang endosperm = 백미라는 추정.** 감사는 "SI가 백미 수분 기준을 인용하므로 백미일 가능성"이라고 적었다. 그 인용(SI Text S3, GB/T 1354-2018 수분 15.5%)은 **섭취 안전 기준 계산**에 들어간 값이고 분석 시료를 설명하는 문장이 아니다. 감사도 이를 추정으로 표시했지만, 근거가 약하다는 점을 덧붙인다. 원문 기준으로는 "보고 안 함"이다.
3. **Kim 뿌리 건물중을 "하한"으로 보는 해석.** 감사 §1.4와 `raw_si/rice_organ_biomass_literature.csv` 머리말은 "일부 포기만 세척 → 하한"이라고 쓴다. 그러나 Kim 2003은 씻지 않은 포기의 뿌리량을 씻은 포기의 지상부:뿌리 건물중 비로 **추정했다**고 명시한다(p.829). 따라서 일부만 씻었다는 사실이 하한의 근거가 되지는 않는다. 블록 깊이(15 cm) 아래의 뿌리가 빠졌을 가능성은 있지만, 원문은 그것을 말하지 않는다.
4. **Kim 2003 고사엽을 "하한"으로 보는 해석.** CLAUDE.md §6은 Kim의 고사엽 비율을 "detached blades are lost"이므로 하한이라고 쓴다. Kim 2003은 탈락 잎몸을 모았는지 잃었는지 **보고하지 않는다.** 레포의 추론이지 원문의 진술이 아니다.
5. **Yang 2006 재식밀도.** `raw_si/yang2006_wuxi_face_dm.csv` 머리말의 "17.5 x 30 cm = 24 hills m-2"는 원문을 옮긴 것이지만, 원문 자체가 계산과 맞지 않는다(19.05). Hasegawa 2017 SI도 24로 적었다.
6. **Gao 2021 전사 머리말**(`raw_si/gao2021_fig3_digitised.csv`, `gao2021_qingan_phenology.csv`). 원문은 DWST를 "dry biomass of stems"로 쓰고 "green stems"라고 하지 않는다. Table A1의 생육 DOY에는 연도가 적혀 있지 않다(성숙 DOY 261이 2018년 수확일 9월 18일과 맞으므로 "2018년"은 타당한 추론이지만 원문 진술은 아니다).
7. 그 밖에 Tang과 Li에 관한 감사 기록(탈이온수·증류수 세척, 피막 제거 없음, 동결건조, 건중 기준, Li의 "잎(줄기 포함)"과 표 사이의 불일치, Li 낟알 정의 미보고, Tang 담수 2–3 cm·탈이온수 관개·무시비·150일·매월 채취)은 원문과 **일치한다.**

---

## 5. 확인하지 못한 자료 / 원문에 없는 항목

### 5.1 확인하지 못한 자료

- **Kim et al. 2001**(*New Phytol.* 150:223). 이번 세션의 파일(`scratchpad/kimface/kim2001.pdf`)은 논문이 아니라 JavaScript 페이지다. 레포는 이 논문의 Table 1(뿌리 건물중, 블록 깊이 20 cm 등)을 `raw_si/rice_organ_biomass_literature.csv`에서 쓰지만, 이번에 확인하지 못했다.
- **Kropff et al. 1994a.** ORYZA 매뉴얼이 IR72DS92 시험의 정보를 이 문헌으로 미룬다. 채취와 건조 절차는 이 문헌에 있을 수 있으나 제공되지 않았다.
- **Okada et al. 2001**(Shizukuishi FACE 시스템), **Liu et al. 2002**(Wuxi FACE 시스템), **Han et al. 2009**(Gao의 LAI 식 출처). 원 논문이 인용만 했고 이번에 읽지 않았다.
- **Hasegawa et al. 2017 본문.** SI의 Table S1만 읽었다.
- **Bouman & van Laar 2006.** 이번에 읽지 않았다(같은 IR72DS92 자료는 매뉴얼 §7.2로 확인).
- **Yamazaki 2023, Kim 2019, Liu 2023.** 이번 세션에 PDF가 없어 이전 감사 기록을 그대로 인용했다.
- **입수 요청용 서지** (DOI는 OpenAlex에서 확인, DOI가 없는 것은 인용한 논문의 참고문헌 목록에서 옮김):
  Kim et al. 2001 `10.1046/j.1469-8137.2001.00111.x` · Okada et al. 2001 `10.1046/j.1469-8137.2001.00097.x` ·
  Hasegawa et al. 2017 `10.1038/s41598-017-13582-y` · Bouman & van Laar 2006 `10.1016/j.agsy.2004.09.011` ·
  Yamazaki et al. 2023 `10.1021/acs.est.2c08767` · Kim et al. 2019 `10.1016/j.scitotenv.2019.03.240` ·
  Liu et al. 2023 `10.1016/j.scitotenv.2022.159826`. DOI 없음: Kropff, van Laar & Matthews (eds) 1994a, *ORYZA1: an
  ecophysiological model for irrigated rice production*, SARP Research Proceedings, IRRI/AB-DLO (ORYZA2000 매뉴얼
  참고문헌) · Liu et al. 2002, *Chin. J. Appl. Ecol.* 13:1253–1258 (Yang 2006 참고문헌) · Han, Peng & Wei 2009,
  *Water Sav. Irrig.* 12:5–7 (Gao 2021 참고문헌).
- **Tang·Li SI의 쪽 번호.** Word 파일이라 고정 쪽이 없다(LibreOffice로 변환하면 Tang SI가 19쪽이 되지만, 원문 표지는 "12 pages"라고 쓴다). 그래서 Text/Table 번호로 인용했다.

### 5.2 이번에 읽은 논문 어디에도 없는 항목

| 항목 | 상태 |
|---|---|
| 기관별 생중량, 기관 수분 함량(θ_fw)의 계절 변화 | 모든 논문 보고 안 함. Tang은 생중량을 쟀으나 값 미공개 |
| 채취 시각 | 모든 논문 보고 안 함 |
| 저울 정밀도 | 모든 논문 보고 안 함 |
| 뿌리 세척에 쓴 체 눈금 | 모든 논문 보고 안 함 |
| 철·망간 피막 제거 | 모든 논문 보고 안 함(이전 감사의 데이터셋도 전부 미제거) |
| 탈락 잎몸 처리 | 모든 논문 보고 안 함 |
| 이삭축 분리 | 모든 논문 보고 안 함 |
| 탈곡 방법, 임실/불임 분리 방법 | 모든 논문 보고 안 함(Kim 2003은 임실 영화를 셌다고만 씀) |
| 제현·도정 조건, 현미/백미 질량 | 생장 논문 4편, Tang, Li 보고 안 함(이전 감사 기준 Yamazaki는 백미:현미 비, Kim 2019는 현미) |
| 천립중 | 모든 논문 보고 안 함 |
| 회당 포기 수 | Gao, Tang(월별), Li(부분시료당) 보고 안 함 |
| 생육단계 사건 날짜 | Kim 2003, Yang 2006, Tang, Li 보고 안 함 |
| 육묘 방식 | Yang 2006, Gao 2021 보고 안 함 |
| 엽면적 측정 기기 | Kim 2003, Sasaki 2005, Yang 2006 보고 안 함 |
