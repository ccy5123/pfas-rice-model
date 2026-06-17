# sci-adk 엄밀성 심사 — PFAS-벼 구획모델

> 도구: **sci-adk** (Scientific Agentic Discovery Kit, https://github.com/ccy5123/sci-adk)
> 대상: 이 저장소의 PFAS-벼 4구획 동적흡수 모델
> 재현: `python sci_adk_review/build_review.py` → `sci-adk verify sci_adk_review/runs/pfas-rice`
> 기준 커밋: 이 브랜치(`claude/practical-brahmagupta-q4tlv3`) HEAD

---

## 0. 한 줄 요약

이 모델의 **구조·메커니즘은 통과(SUPPORTED)**하지만, **경험적 예측 주장은
엄밀성 게이트를 통과하지 못한다(REFUTED)**. 모델에 동봉된 데모 BAF를 정직하게
입력하면 sci-adk는 **인증을 거부(HALT)**한다 — 이는 바로 sci-adk가 만들어진
계기인 "rice-failure"의 재현이자 차단이다.

| 가설 | referent | 판정 | 근거 |
|---|---|---|---|
| H1 질량보존 | formal | **SUPPORTED** | 잔차 < 1e-5 (테스트 검증, generated) |
| H2 음이온 배제(e^N≈107) | formal | **SUPPORTED** | E_m·z의 closed-form 귀결 (generated) |
| H3 Yamazaki = 표본외 예측검증 | empirical | **REFUTED** | RMSE 0.029는 **포화 in-sample 적합** (measured) |
| H4 곡립 위해성 예측 정확 | empirical | **REFUTED** | Tang 곡립 3–8× 과소예측 (measured) |
| (trap) 데모 BAF = 현장 예측 | empirical | **HALT** | 데모는 **synthetic_proxy** → 인증 거부 |

`sci-adk verify` 결과: 4개 claim 전부 기록으로부터 재현됨(REPRODUCED), exit 0,
record digest `sha256:8c67be34…`.

---

## 1. sci-adk란 무엇인가 (왜 "배경/목표/방법/예상결과물"을 먼저 쓰게 했나)

sci-adk는 실험을 **수행하는 도구가 아니라**, 연구를 **동결된 기준으로 심판하는
심판/기록자(referee/scorekeeper)**다. 핵심 철학은 **기록(record)과 믿음(belief)의
분리**다.

- **기록(Evidence)**: 단조(monotone)·append-only. null·negative 결과도 1급 시민.
- **믿음(Claim)**: 비단조(non-monotone)·수정 가능. 증거가 바뀌면 강등/철회된다.
- 운영 규칙: *"agents propose; the engine judges by frozen criteria. No
  self-certification."* (에이전트는 제안하고, 엔진이 동결 기준으로 판정하며,
  자기인증은 없다.)

입력 형식이 바로 **4-pane proposal — Background / Goal / Method / Expected Output**
(한글 헤더 `연구 배경 / 연구 목표 / 연구 방법 / 기대 산출물` 지원)이다. 이것을
`sci-adk run`이 **동결된 사전등록 Spec**으로 컴파일한다(anti-HARKing: 결과를 본 뒤
질문·판정기준을 바꾸지 못함). → **그래서 당신이 나에게 그 네 절을 쓰게 한 것이다.**

### 결정적 사실: sci-adk는 이 프로젝트의 실패에서 태어났다

sci-adk 소스가 직접 밝힌다.

- `src/sci_adk/core/validity.py`:
  > *"This is the load-bearing fix for the **rice-failure defect**: a run on an
  > EMPIRICAL proposal used SYNTHETIC data and the harness reported '4/4 SUPPORTED'."*
- `src/sci_adk/core/evidence.py` (Provenance):
  > `synthetic_proxy` — a FABRICATED stand-in for an external referent the data does
  > not contain (**the rice numbers**).

즉, 과거에 이 PFAS-벼 모델(또는 그 전신)을 합성/예시 데이터로 돌렸더니 모든 주장이
"SUPPORTED"로 잘못 인증된 사건이 있었고, sci-adk의 **evidence-validity 게이트**는
바로 그것을 막으려고 만들어졌다. 이번 심사는 그 도구를 **원인이 된 바로 그
프로젝트에 다시 적용**한 것이다.

---

## 2. 무엇을 했나 (방법)

`build_review.py`가 재현 가능한 한 번의 실행으로 다음을 수행한다.

1. **동결 Spec 구성** (`runs/pfas-rice/spec.json`). 휴리스틱 파서는 referent
   클래스·규칙 종류·비순환성을 추론하지 못하므로(README의 "capability가 미리 만든
   Spec을 공급" 경로대로) `proposal.md`의 네 절을 그대로 `raw_proposal`로 싣고,
   4개 가설에 **정직한 referent/규칙/비순환 입증**을 부여했다.
2. **증거 분류** (`runs/pfas-rice/evidence/`). 각 증거를 정직하게 분류:
   - `generated` — 모델의 **형식적 성질**(질량보존, e^N). 모델 내부에서 참인 것.
   - `measured` — **실측 데이터** 비교(Yamazaki/Tang/Kim).
   - `synthetic_proxy` — **예시용 데모 수치**(보정되지 않은 placeholder).
   수치는 모델의 **실제 출력**이다: `reproduce_demo.py` 실행값(RMSE 0.029, PFOA
   root 0.49/0.49), 테스트 통과(질량보존, e^N≈106.8), 문헌 검증(Tang/Kim).
3. **Verdict 작성** (`runs/pfas-rice/verdicts/`). 정성 가설은 in-session 에이전트가
   chief-over-N **판정 trail**을 직접 써야 한다(엔진은 trail 없는 binding 판정을
   F2 게이트로 거부). 각 trail은 동결 rubric을 복사해 자기완결적이다.
4. **루프 + 감사**: `run_checkpoint_loop`이 컴파일→증거 영속화→판정 바인딩을 하고,
   `verify_run` / `sci-adk verify`가 **LLM 없이** 기록만으로 믿음을 재도출한다.

증거는 referent 클래스에 따라 **adequacy 게이트**를 통과해야 한다:
- 경험적(empirical) 가설의 binding 판정은 **≥1개 measured 증거**가 있어야 한다.
- empirical 가설에 **synthetic_proxy** 증거가 닿으면 **무조건 HALT**(범주 오류).
- formal 가설이 generated 증거로 binding되면 **비순환성 입증**을 요구한다(Guard 2).

---

## 3. 결과 (믿음 상태)

### H1 질량보존 — SUPPORTED (formal)
`tests/test_plant_module.py::test_mass_conservation_source_is_root_uptake`가
`dmass == src` (rel 1e-6/abs 1e-9)로 통과. 잔차는 적분기 **밖에서** 독립 회계로
계산되므로(비순환성 입증) 솔버가 보존을 가정한 게 아니다. 임계규칙(< 1e-5) 자동 충족.

### H2 음이온 배제 — SUPPORTED (formal)
`N = z·E_m·F/RT = +4.67` (E_m=−120 mV, z=−1) → `e^N = 106.8`. 이는 Tier-0
입력의 **닫힌 형식 귀결**이지 적합값이 아니며(비순환성 입증), 운반체(Vmax/Km)가
이 배제를 극복해야 한다는 점이 "수동 확산만으로 불충분"이라는 규칙을 충족.

### H3 Yamazaki 표본외 예측검증 — REFUTED (empirical)
`reproduce_demo.py`가 11개 동족체 root/straw/grain BAF를 **log10 RMSE 0.029**로
재현한다. 그러나 W2 적합은 **동족체당 ~3개 transport 파라미터를 3개 관측에**
맞춘 **포화(saturated) 적합**이다 → pred≈obs(0.49/0.49)는 *예측*이 아니라 *재현*의
서명이다. CLAUDE.md 자신이 "a saturated W2 fit … reproduction is guaranteed, NOT
predictive validation"이라 명시한다. **measured 데이터에 근거**하지만, "표본외
예측검증"이라는 가설은 **거부**된다(모델은 Yamazaki와 *일관*될 뿐 *검증*된 게 아님).

### H4 곡립 위해성 예측 — REFUTED (empirical)
measured 비교에서 모델은 곡립을 구조적으로 **3–8× 과소예측**한다(Tang PFOA
endosperm 0.11 vs 0.95, dw; `docs/tang2026_grain_units_exploration.md`의 단위
보정 후에도 닫히지 않음). 유일한 "일치"(Kim)는 **곡립 단일점 in-sample 앵커**로
L_Ph를 강제한 것뿐이다. 위해성 평가가 의존하는 바로 그 구획의 구조적 과소예측은
식이 위해성 평가를 뒷받침할 수 없다 → 거부. **정직한 음성 결과는 1급 시민이다.**

### trap — HALT (synthetic_proxy 거부)
동봉 데모(`_demo()`)의 BAF를 경험적 가설에 입력하자 게이트가 멈췄다:

```
synthetic_proxy Evidence bears on an empirical hypothesis -- a category error:
a fabricated stand-in does not contain the external referent the empirical claim
is about ... See design/evidence-validity.md Guard 3.
```

Claim은 **쓰이지 않았다**(`runs/pfas-rice-trap/VALIDITY_HALT.txt`). 보정되지 않은
예시 수치가 "예측 능력"으로 self-certify되는 길이 구조적으로 막혔다 — rice-failure의
재현이자 차단.

---

## 4. 이것이 프로젝트에 시사하는 것

1. **모델의 기여는 메커니즘에 있다, 예측에 있지 않다 (아직).** 음이온 배제·terminal
   축적·질량보존·basis-A 결합은 견고하다(formal, SUPPORTED). 헤드라인으로 자주
   인용되는 "RMSE 0.029"는 검증이 아니라 **포화 재현**이며, 엄밀성 관점에서 예측
   주장의 근거가 되지 못한다.
2. **곡립은 구조적 한계**이고, 정직한 음성으로 기록되어야 한다(이미 CLAUDE.md가
   인정하는 바를 독립 엔진이 확인).
3. **진짜 검증으로 가는 길**: (a) f_xy/L_Ph를 적합한 데이터와 **분리된** 표본외
   데이터셋으로 평가, (b) 곡립 구조 항(비가역/이력 흡착, 단위·fw/dw) 재정식화,
   (c) 데모를 "예시"로 명확히 라벨링(이미 코드 주석엔 있음) — 결과 보고에서
   `synthetic_proxy`로 취급.

---

## 5. 산출물 (`sci_adk_review/`)

```
proposal.md                     # 동결된 4-pane 사전등록(배경/목표/방법/기대산출물)
build_review.py                 # 재현 드라이버(Spec/Evidence/Verdict + 루프 + verify + trap)
runs/pfas-rice/
  spec.json                     # 동결 Spec v1 (4 가설, referent/규칙/비순환)
  evidence/*.json               # 분류된 증거 (generated/measured) + 문헌(prior-work)
  verdicts/*.json               # chief-over-N 판정 trail (H2/H3/H4)
  claims/*.json                 # 믿음 상태 (SUPPORTED×2, REFUTED×2) + 이력
  paper/draft.tex               # 결정론적 논문 초안 스켈레톤
runs/pfas-rice-trap/
  VALIDITY_HALT.txt             # 거부 증거 (synthetic_proxy → HALT)
```

### 재현
```bash
pip install -e /path/to/sci-adk          # 또는 PYTHONPATH=sci-adk/src
pip install numpy scipy pytest           # 모델 실행/근거 수치용
python sci_adk_review/build_review.py
sci-adk verify sci_adk_review/runs/pfas-rice   # exit 0, 전 claim REPRODUCED
```

> 주의: `runs/`의 증거 수치는 기준 커밋의 모델 출력을 그대로 옮긴 것이고, 각
> 증거는 `provenance.code_ref`(커밋)와 `environment`(실행 경로)를 기록한다.
> 모델 코드가 바뀌어 출력이 달라지면 증거를 **새 항목으로 append**해야 한다
> (Evidence는 append-only). prior-work는 모델이 실제로 인용하는 문헌
> (Yamazaki/Tang/Kim/Brunetti …)을 LITERATURE 증거로 기록했다.
