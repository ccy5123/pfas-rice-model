# PFAS–벼 4구획 흡수 모델

[English](README.md) | **한국어** &nbsp;·&nbsp; 앱 매뉴얼: [English](docs/MANUAL_EN.md) | [한국어](docs/MANUAL_KR.md)

**논벼(*Oryza sativa*)의 영구 음이온 PFAS** 흡수를 다루는 기계론적 모델입니다. DPU/Trapp 체계를
이온성 유기화합물(IOC)로 확장했고, 뿌리 → 줄기 → 잎 → 낟알의 네 구획으로 이루어집니다. Yamazaki
관측이 있는 12개 congener(PFCA C4–C12, PFSA C4/C6/C8)와 GenX(잠정), SMILES로 주는 임의의 PFAS,
log Kow로 주는 중성 유기물을 다룹니다.

> **정본**: `CLAUDE.md`(현재 상태와 결과), `docs/ASSUMPTIONS.md`(모든 입력의 출처, 또는 측정값이
> 없을 때의 민감도 범위), `params/parameters.json`(파라미터), `docs/MANUAL_KR.md` / `docs/MANUAL_EN.md`
> (앱). `docs/HANDOFF_*.md`, `docs/VALIDATION_*.md`, `docs/*_exploration.md`, `validation/archive/`는
> 기록이지 현재 상태가 아닙니다.

---

## 빠른 시작

```bash
pip install -r requirements.txt
python reproduce_demo.py          # Yamazaki 2020 사이트-연도, 건물중, 적합 수송(재현)
python reproduce_demo.py --rec    # 단조 물리적 f_xy: 짚은 예측
pytest                            # 테스트
```

```python
import model_api as api                       # (src/를 경로에)
r = api.simulate("PFOA")                      # 보정 사이트-연도, 평균 수확일
r["baf_final_dw"]                             # 뿌리/줄기/잎/낟알/짚, 건물중
api.yamazaki_harvest_baf("PFOA")              # 관측 대응값: 두 번의 수확, 기하평균
api.simulate("PFOA", scenario="generic", season=120.0, biomass="oryza2000")   # 다른 드라이버
api.simulate("PFOA", lipid_loading=True)      # K_PL 게이트 지질 적재(선택)
api.simulate_nstem_leaf("PFOA", lipid_loading=True)   # 재분배 지상부 + 지질
api.simulate_neutral(2.45, name="carbamazepine")      # log Kow로 주는 중성 유기물
```

## 모델이 쓰는 입력

- **드라이버** — 기본 시나리오는 보정 사이트-연도(`src/site_scenario.py`)입니다. JMA 쓰쿠바 2020 기상,
  Yamazaki 이앙일에 맞추고 출수·수확일로 생육 단계를 보정한 ORYZA2000 잠재생산 포트(`src/oryza2000.py`,
  IR72.D92), FAO-56 항으로 계산한 ORYZA2000 증산을 씁니다. `scenario="generic"`은 기후값(`forcing_rice`)과
  고른 바이오매스 드라이버를 고른 기간 동안 씁니다.
- **기관 질량은 생체중(fresh)** 입니다. 작물 모델은 건물을 내므로 변환합니다(`src/organ_mass.py`).
- **수확 관측 연산자**(`src/harvest.py`): 작물 모델의 잎 고사로부터 죽은 잎을 재구성해 잎과 짚에
  포함하고, 결과를 건물중으로 비교합니다.
- **수송 파라미터**는 사이트-연도에서 적합합니다(`validation/refit_site_year.py`, 사전 등록).
  `f_xy_source="recommended"`(기본; 단조 물리적 f_xy, L_Ph와 κ_d는 뿌리·낟알에 맞춤) 또는 `"fit"`(세
  파라미터 모두 congener별 적합; 포화된 재현).

## 결과 요약 (자세한 내용과 스크립트: `CLAUDE.md` §6)

{RESULTS_KR}

일반인 화면은 모든 절대 수치에 ×/÷ 밴드를 함께 보여주며, 그 폭은 기본 설정의 전이 오차
(`model_api.APRIORI_OOS_SCORES`)에서 정해집니다.

## 토양 측 — HYDRUS-1D (Method A)

내장된 HYDRUS-1D 4.08 소스를 gfortran으로 빌드하고 `phydrus`로 구동합니다. 사이트-연도 증산이
토양 계산을 구동하고, 공극수 `C_w^o(t)`와 뿌리 수분 흡수를 돌려줍니다.

```bash
cp external/hydrus_source/makefile external/hydrus_source/source/
(cd external/hydrus_source/source && make)
pip install phydrus
python validation/hydrus_coupled_run.py
```

흡착이 약한 단쇄는 담수기에 용탈되고({HYDRUS_KR}), 장쇄는 완충됩니다. `simulate(cwo_profile="flooded")`는
엔진 없이 쓰는 해석적 형상이며, `k_leach`는 congener별로 HYDRUS 결과에 맞춰 보정돼 있습니다. 웹의
Claude Code에서는 SessionStart 훅이 엔진을 빌드합니다.

## 대화형 앱

```bash
pip install -r requirements.txt
streamlit run app.py
```

농도로 색칠한 식물·토양 지도, 조직 동역학, 드라이버, Yamazaki 비교, Tang 점검, 베이지안 역추정(측정한
조직 농도로부터 노출 추정), 6가지 노출 모드(파라메트릭, 사용자 표, HYDRUS CSV, 실시간 HYDRUS-1D, 토양
재고량, 바이오모니터링), 메커니즘 스위치(carrier/bypass, 지질 적재), 세 가지 화합물 종류(congener,
SMILES, 중성 log Kow — EPA CompTox 조회 선택)를 제공합니다. 일반 모드는 한국어, 전문가 모드는 영어입니다.
`docs/visualization_tool.md`와 매뉴얼을 참고하세요.

## 구조

```
├── reproduce_demo.py / build_parameters.py / app.py
├── src/            모델 코드(핵심 ODE, 수확 연산자, 사이트 시나리오, ORYZA2000 포트, API, 그림, 토양)
├── params/         parameters.json + 원천 표 + refit_site_year.csv + cwo_kleach.csv
├── data_obs/       관측(Yamazaki 건물중, Li 2025, 중성 물질 표)
├── validation/     유지되는 분석; archive/ = 대체된 분석
├── docs/           ASSUMPTIONS.md, 모델 보고서(.tex), 문헌 DB, 매뉴얼
├── external/hydrus_source/   내장 HYDRUS-1D 소스
└── tests/
```

## 한 화면 요약(모델)

```
B_k = θ_fw + (1 − θ_fw) · ( f_prot·K_prot + f_PL·K_PL + f_cw·K_cw )      [L/kg fw]
```

뿌리 유입 = GHK 전기확산(E_m = −120 mV에서 음이온 배제 e^N ≈ 107) + 포화 Michaelis–Menten carrier.
물관 적재 `f_xy · C_root/B_root`(+ 선택적 지질항). 체관은 두 흐름입니다: 낟알로는 건물 생장에 필요한
수액, 뿌리로는 증산의 φ 비율. 낟알과 잎은 종단 축적체라 모델은 동역학이어야 합니다.

## 남은 과제

가정을 대체할 측정은 `docs/ASSUMPTIONS.md` §7에 있습니다(Yamazaki SI 원자료, carrier 동역학, 이삭의
증산 분담, 벼 체관 수액량, 생육기 기관 수분, 벼 뿌리 세포벽·plaque 흡착). 전이 데이터셋은 자체 장소·
일정이 레포에 없어 보정 사이트-연도 드라이버로 돌고 있습니다.
