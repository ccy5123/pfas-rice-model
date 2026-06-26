# HANDOFF — 리뷰어 평가 기반 후속 작업 (next sessions)

> 일반인/전문가 UI 개편 + 한국어화 + 베이지안 역추정 + 데이터 표 + 이중언어 문서 작업
> 직후, **리뷰어 관점 평가**에서 도출한 후속 작업 목록입니다. 종합 평결은
> **"B+ / 연구·교육·시각화 도구로는 릴리스 가능(조건부), 예측·의사결정 도구로는 부적합"**.
> 본 문서는 그 리뷰의 **P1/P2/P3 항목**을 실행 가능한 작업으로 정리합니다.
>
> 현재 상태(기준점): `main` 기준 전체 테스트 **178 passed, 2 skipped**. UI/문서 작업 동안
> `parameters.json`·모델 수식 **불변**. `app.py` 1208줄(단일 파일), 문서 20개, lint nit 일부.

## 작업 우선순위 한눈에
- **P1 (공개 전 권장, 윤리적 핵심)**: 일반인 모드 절대 수치에 **불확실성/등급 가시화**.
- **P2 (사용성 다듬기)**: 역추정 지연, 헤드라인↔지도 불일치, 일반인 모드 용어 누수, 낟알 BAF<1 문구, PNG 내보내기 폴백.
- **P3 (유지보수)**: `app.py` 모듈화, i18n 정리, 문서 정본 표시, lint 정리.
- **별첨**: 기존 과학 백로그(이번 리뷰 범위 밖이나 분실 방지용 포인터).

> 공통 원칙(이번 작업들 전부): **UI/UX 레이어만**. `parameters.json`·모델 수식·`reproduce_demo`(RMSE 0.029)는
> 건드리지 말 것. 새 계산 헬퍼는 `model_api`에 순수 함수로, 그림은 `plots`에 `lang` 인자 유지로. 테스트(`test_model_api`/
> `test_plots`) 추가/유지. 일반인=한국어·전문가=영어 규약 유지.

---

## P1 — 일반인 모드 수치에 불확실성 가시화 (가장 중요)

**문제**: 일반인 모드가 "낟알 0.15 µg/kg" 같은 **소수점 정밀 수치**를 비전문가에게 보여 주는데, 모델의 정직한
a-priori 예측오차는 큼(절대 BAF는 lipid opt-in 외 OOS 예측 실패; 헤드라인 0.029는 in-sample 포화 fit).
면책 배너만으로는 "정밀해 보이는 숫자 → 과신" 인지 간극을 막지 못함. **이 프로젝트에서 유일하게 '윤리적으로' 신경 쓰이는 지점.**

**해야 할 일(택1 또는 조합)**:
1. **수치에 범위/밴드 부착** — 메트릭 카드와 요약 문장의 절대값 옆에 ± 또는 "0.1–0.5" 형태.
   - 근거 소스: (a) 베이지안 forward 불확실성(이미 inverse에서 Laplace 폭을 계산하는 기계가 있음 →
     `model_api.estimate_exposure_bayesian`의 사후폭 산출 로직을 forward 방향으로 재사용), 또는
     (b) 문서화된 a-priori 예측오차(log10 RMSE ~0.84–0.95)를 곱으로 환산한 정성 밴드(예: ×/÷ 약 7배).
   - **권장**: 정량 밴드가 과해 보이면, **정성 등급**("대략적 추정 — 실측과 수배 차이 가능")을 카드 하단 캡션 + 요약 문장에 명시.
2. **절대값 대신 상대/등급을 전면에** — 일반인 헤드라인을 "토양수 대비 ○배(낮음/보통/높음)" 같은 **상대 비교/등급** 중심으로,
   정밀 µg/kg는 보조로.
3. **신뢰 신호 일관화** — 지도/그래프에도 "모델 추정(실측 아님)" 워터마크/캡션을 더 분명히.

**위치**: `app.py` 일반인 헤드라인(`if not expert:` 의 metric 카드 + `st.info` 요약), `📈/📊` 탭 캡션.
**수용 기준**: 일반인 기본 화면에서 **모든 절대 수치 옆에 불확실성/등급 신호가 보임**. 전문가 모드·`parameters.json`·테스트 불변.
새 헬퍼는 `model_api`에 순수 함수 + 테스트.

---

## P2 — 사용성 거친 모서리

### P2-1. 🔎 거꾸로 추정 지연 (8–22초)
- **문제**: 첫 실행이 길고(특히 PFBA ~22초; ODE ~1–2.5s/회 × ~8회) 무반응처럼 보임. 버튼 게이트는 되어 있음.
- **해야 할 일**:
  - 진행 표시 강화(단계별 `st.status`/스피너 텍스트 "ODE 8회 중 n회…"), **또는**
  - 속도 개선: forward 평가 캐시 공유, 또는 추정 그리드/패스 축소(현재 2-pass 6 ODE), 또는 짧은 시즌/저해상 그리드로 초기 추정 후 정밀화.
- **위치**: `model_api.estimate_exposure_bayesian`, `app._render_inverse_estimator`, 캐시 `app._estimate_exposure`.
- **수용 기준**: 일반 사용자가 진행 상태를 인지; 또는 최악 케이스 < ~10초.

### P2-2. 헤드라인 ↔ 지도 불일치
- **문제**: 요약은 "대부분 뿌리에 남음"인데, 지도는 **잎**이 가장 뜨겁게 칠해짐(잎 농도 PFOA 2.52 vs 짚 평균 0.43).
  헤드라인은 root/straw/grain(짚=줄기+잎 평균)으로 판단하기 때문.
- **해야 할 일**: 둘 중 하나로 일관화 — (a) 헤드라인 "대부분 ○○" 판단에 잎을 분리 반영, 또는
  (b) 지도/캡션에 "부위별 **농도**(짚 카드는 줄기+잎 평균)"임을 명시해 혼선 제거.
- **위치**: `app.py` 일반인 헤드라인 `where_most` 계산 + 지도 탭 캡션.

### P2-3. 일반인 모드 용어 누수
- **문제**: "🔬 실제 측정값과 비교" 펼침의 `plots.fig_baf` 축이 영어 `BAF [L/kg]` — "기호 0개" 목표와 충돌(접힌 expander라 경미).
- **해야 할 일**: `fig_baf`에 `lang` 인자 추가(다른 plain 빌더처럼 기본 `"en"`), 일반인 호출 시 `lang="ko"`로 "축적 배수 [L/kg]" 등.
- **위치**: `src/plots.py::fig_baf`, `app.py` 일반인 `📊` 탭. 테스트 영어 기본 유지.

### P2-4. 낟알 BAF<1 문구
- **문제**: 요약의 "약 0.2배"가 어색하게 읽힘(낟알 BAF<1일 때).
- **해야 할 일**: BAF<1이면 "토양수보다 **낮음**(약 1/5 수준)" 식으로 자연스럽게 표현 분기.
- **위치**: `app.py` 일반인 `st.info` 요약 문장.

### P2-5. PNG 내보내기 폴백
- **문제**: `kaleido` 기본 미설치 → 대부분 사용자는 캡션만 봄.
- **해야 할 일**: (a) 폴백 안내를 더 친절하게(설치 명령 1줄), 또는 (b) `requirements-app.txt`에 kaleido 안내 강화, 또는
  (c) matplotlib 기반 정적 대체 이미지 제공(선택). 최소한 안내 문구 개선.
- **위치**: `app._png_bytes` 호출부의 캡션(일반인/전문가 양쪽).

---

## P3 — 유지보수성

### P3-1. `app.py` 모듈화 (1208줄 단일 파일)
- Simple/Expert 분기·헬퍼가 한 파일에 누적. `app/` 하위로 분리 권고: `ui_sidebar.py`, `ui_simple.py`, `ui_expert.py`,
  `ui_common.py`(역추정·표·다운로드·footer). `app.py`는 조립만.
- **주의**: Streamlit 캐시 데코레이터/세션 상태 키 충돌 없게. 동작·스크린샷 동치 확인.

### P3-2. i18n 정리
- 현재 한국어가 **inline 삼항**(`"…" if ko else "…"`)로 흩어져 있어 3번째 언어 추가 시 고통.
- 경량 문자열 테이블(`STRINGS = {"ko": {...}, "en": {...}}` + `t(key, lang)`) 또는 간단 `gettext`로 추출.
- `plots`의 `lang` 인자 패턴은 유지(좋음).

### P3-3. 문서 정본 표시
- 문서 20개 sprawl → 신규 진입자가 "무엇이 정본인지" 헷갈림. CLAUDE.md §6는 사실상 changelog(온보딩 아님).
- README/README_KR 상단에 **"정본 = `CLAUDE.md`(컨텍스트) + `parameters.json`(파라미터) + `docs/MANUAL_*`(사용) + `docs/OVERVIEW_KR.md`(과학)"**,
  나머지 `HANDOFF_*`/`*_exploration.md`는 **이력/탐색**임을 1줄 명시. 또는 `docs/INDEX.md` 추가.

### P3-4. lint 정리
- pyflakes: `app.py:971` 등 f-string-without-placeholder, `calibration.py`/`oryza_growth.py`의 unused `field`,
  `model_api.py` rdkit availability-probe import, `pfas_rice_plant_module.py`의 re-export(별칭이라 의도적 — `# noqa`/`__all__`로 명시).
- 무해하나 정리 가치. CI에 `ruff`/`flake8` 경량 추가 고려(과하지 않게).

---

## 별첨 — 기존 과학 백로그 (이번 리뷰 범위 밖, 분실 방지 포인터)

리뷰는 **UI/UX**에 집중했으나, 모델 자체의 미해결은 이미 잘 기록돼 있음(중복 작성 말 것, 포인터만):
- **장쇄(C10–C12) 지상부 과소예측** 구조적 미해결 — `docs/fxy_longchain_lipid_exploration.md`, `docs/twopool_root_exploration.md`, PR #21 LC1–LC6.
- **`f_xy` 절대 스케일/조건 의존성**(PFOS Yamazaki 0.14 vs Tang 0.32) — `docs/DELIVERABLE_GAP_B_fxy.md`, `docs/VALIDATION_TANG2026_NSTEM_KR.md`.
- **`K_cw` 직접 측정 부재**, 에터/술폰아미드 Koc 공백 — `docs/DELIVERABLE_GAP_A_Kcw.md`, `docs/structure_input.md`.
- **PFSA 전용 수송항** 필요 — README Status / CLAUDE.md §9 task 8.
- **현장 측정 입력**(측정 `Q_TP(t)`/`M_s(t)`, in-situ E_m, 사용자 현장 토양·담수 일정) — CLAUDE.md §9.
- 두-pool seq 모델의 `parameters.json` **승격 여부 결정** — `docs/HANDOFF_BAF_twopool.md`.

이들은 **데이터/실험 의존**이거나 별도 모델링 트랙으로, 본 UI 후속과 분리해 진행.

---

## 권장 진행 순서 (다음 세션들)
1. **P1**(불확실성 가시화) — 공개 전 가장 중요, 단일 세션으로 가능.
2. **P2-2/P2-3/P2-4**(헤드라인·용어·문구) — 가볍고 즉효, 한 세션에 묶기.
3. **P2-1**(역추정 속도/표시) — 측정 후 결정.
4. **P3-1**(app.py 모듈화) — 큰 구조 변경, 별도 세션. 이후 P3-2(i18n)와 함께.
5. **P3-3/P3-4**(문서 정본·lint) — 자투리.

## 각 작업 공통 수용 기준
- `parameters.json`·모델 수식·`reproduce_demo`(RMSE 0.029) **불변**.
- `test_model_api.py`/`test_plots.py` 통과(+신규 커버리지). 전체 suite green 유지.
- 일반인=한국어·전문가=영어 규약 유지(`plots` 빌더 영어 기본값 유지 → 테스트·전문가 영향 0).
- 헤드리스 Streamlit + Playwright 스크린샷으로 일반인/전문가 양쪽 회귀 확인.

## 다음 세션 재개 프롬프트(예시)
> "`docs/HANDOFF_review_followup.md`의 **P1**을 구현: 일반인 모드 절대 수치에 불확실성/등급을 가시화.
> `parameters.json`·모델 수식 불변, UI/UX만. `model_api`에 순수 헬퍼 + 테스트, 일반인=한국어 유지.
> 헤드리스 Streamlit+Playwright로 일반인 화면 회귀 확인."
