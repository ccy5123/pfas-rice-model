# sci-adk 엄밀성 심사 — PFAS-벼 구획모델

> 도구: **sci-adk** (Scientific Agentic Discovery Kit, https://github.com/ccy5123/sci-adk)
> 대상: 이 저장소의 PFAS-벼 4구획 동적흡수 모델
> 입력 Spec: `proposal.md` (연구 배경/목표/방법/기대 산출물 **전문** — 동결된 사전등록)
> 재현: `python sci_adk_review/build_review.py` → `sci-adk verify sci_adk_review/runs/pfas-rice`
> 기준 커밋: 이 브랜치(`claude/practical-brahmagupta-q4tlv3`) HEAD

---

## 0. 한 줄 요약

이 모델의 **구조·메커니즘·계산적 기반은 통과(SUPPORTED)**하지만, **경험적 예측
주장은 엄밀성 게이트를 통과하지 못한다(REFUTED)**. 동봉 데모 BAF를 정직하게
입력하면 sci-adk는 **인증을 거부(HALT)**한다 — sci-adk가 만들어진 계기인
"rice-failure"의 재현이자 차단이다.

가설 7개는 `proposal.md`의 **5개 연구 목표**를 독립 검증 가능한 형태로 형식화한 것
(H7은 사용자 재정의: "예측 말고 fitting으로 구조가 실험을 재현 가능한가").

| 가설 (목표) | referent | 판정 | 근거 |
|---|---|---|---|
| H1 질량보존 (목표1) | formal | **SUPPORTED** | 잔차 < 1e-5 (테스트 검증, generated) |
| H2 음이온 배제 e^N≈107 (목표1) | formal | **SUPPORTED** | E_m·z의 closed-form 귀결 (generated) |
| H3 Yamazaki=표본외 예측검증 (목표4) | empirical | **REFUTED** | 포화 적합 0.029 vs **사전적 예측 0.837** (measured) |
| H4 곡립 위해성 예측 (목표4) | empirical | **REFUTED** | Tang 곡립 3–8× 과소예측 (measured) |
| H5 토양 결합 동족체 의존 (목표3) | formal | **SUPPORTED** | Koc 스프레드 4.4 log10 ~25000× (generated) |
| H6 SMILES read-across 재현 (목표5) | formal | **SUPPORTED** | test_pfas_structure 23/23 (generated) |
| **H7 구조적 적합성: shoot 재현(제약 적합)** | empirical | **SUPPORTED** | **straw ~0.18 @ DOF 20 (ORYZA2000)** (measured) |
| (trap) 데모 BAF=현장 예측 | empirical | **HALT** | 데모는 **synthetic_proxy** → 인증 거부 |

**패턴**: formal/계산적 주장(H1·H2·H5·H6) → SUPPORTED, **표본외 예측(H3·H4) →
REFUTED**, 그러나 **"구조가 fitting으로 실험을 재현 가능한가"(H7, 사용자의 핵심
질문)는 — 당신의 ORYZA2000 biomass로 구동 시 — shoot에서 SUPPORTED**(straw ~0.18,
자유도 20; 전 조직 ~factor 2.2 이내 @ DOF 10), grain 장쇄만 잔여 floor.
`sci-adk verify`: 7개 claim 전부 기록으로부터 재현(REPRODUCED), exit 0,
record digest `sha256:50417cd6…`.

---

## 1. sci-adk란 무엇인가 (왜 "배경/목표/방법/예상결과물"을 먼저 쓰게 했나)

sci-adk는 실험을 **수행하는 도구가 아니라**, 연구를 **동결된 기준으로 심판하는
심판/기록자(referee/scorekeeper)**다. 핵심 철학은 **기록(record)과 믿음(belief)의
분리**다.

- **기록(Evidence)**: 단조(monotone)·append-only. null·negative 결과도 1급 시민.
- **믿음(Claim)**: 비단조(non-monotone)·수정 가능. 증거가 바뀌면 강등/철회된다.
- 운영 규칙: *"agents propose; the engine judges by frozen criteria. No
  self-certification."*

입력 형식이 바로 **4-pane proposal — Background / Goal / Method / Expected Output**
(한글 헤더 `연구 배경 / 연구 목표 / 연구 방법 / 기대 산출물` 지원)이다. 이것을
`sci-adk run`이 **동결된 사전등록 Spec**으로 컴파일한다(anti-HARKing: 결과를 본 뒤
질문·판정기준을 바꾸지 못함). → **그래서 그 네 절을 먼저 쓰게 한 것이고, 그 전문이
이번 심사의 동결 Spec(`proposal.md`, `runs/pfas-rice/spec.json`의 `raw_proposal`)이다.**

### 결정적 사실: sci-adk는 이 프로젝트의 실패에서 태어났다

sci-adk 소스가 직접 밝힌다.

- `src/sci_adk/core/validity.py`:
  > *"This is the load-bearing fix for the **rice-failure defect**: a run on an
  > EMPIRICAL proposal used SYNTHETIC data and the harness reported '4/4 SUPPORTED'."*
- `src/sci_adk/core/evidence.py` (Provenance):
  > `synthetic_proxy` — a FABRICATED stand-in for an external referent the data does
  > not contain (**the rice numbers**).

과거에 이 PFAS-벼 모델(또는 전신)을 합성/예시 데이터로 돌렸더니 모든 주장이
"SUPPORTED"로 잘못 인증된 사건이 있었고, sci-adk의 **evidence-validity 게이트**는
바로 그것을 막으려고 만들어졌다. 이번 심사는 그 도구를 **원인이 된 바로 그
프로젝트에 다시 적용**한 것이다.

---

## 2. 무엇을 했나 (방법)

`build_review.py`가 재현 가능한 한 번의 실행으로 다음을 수행한다.

1. **동결 Spec 구성** (`runs/pfas-rice/spec.json`). `proposal.md`의 네 절 전문을
   그대로 `raw_proposal`로 싣고(휴리스틱 파서가 추론 못하는) referent 클래스·규칙
   종류·비순환 입증을 6개 가설에 정직하게 부여했다(README의 "capability가 미리 만든
   Spec 공급" 경로).
2. **증거 분류** (`runs/pfas-rice/evidence/`). 각 증거를 정직하게 분류:
   - `generated` — 모델의 **형식적/계산적 성질**(질량보존, e^N, Koc QSPR 스프레드,
     SMILES read-across). 모델 내부에서 참인 것.
   - `measured` — **실측 데이터** 비교(Yamazaki/Tang/Kim).
   - `synthetic_proxy` — **예시용 데모 수치**(보정되지 않은 placeholder).
   수치는 모델의 **실제 출력**이다: `reproduce_demo.py`(RMSE 0.029, PFOA root
   0.49/0.49), 테스트 통과(질량보존; e^N≈106.8; structure 23/23),
   `literature_params.koc`(라이브 계산), 문헌 검증(Tang/Kim).
3. **Verdict 작성** (`runs/pfas-rice/verdicts/`). 정성 가설은 in-session 에이전트가
   chief-over-N **판정 trail**을 직접 써야 한다(엔진은 trail 없는 binding 판정을
   F2 게이트로 거부). 수치 가설(H1·H5)은 엔진이 자동 판정.
4. **루프 + 감사**: `run_checkpoint_loop`이 컴파일→증거 영속화→판정 바인딩을 하고,
   `verify_run` / `sci-adk verify`가 **LLM 없이** 기록만으로 믿음을 재도출한다.

증거는 referent 클래스에 따라 **adequacy 게이트**를 통과해야 한다:
- 경험적(empirical) 가설의 binding 판정은 **≥1개 measured 증거**가 있어야 한다.
- empirical 가설에 **synthetic_proxy** 증거가 닿으면 **무조건 HALT**(범주 오류).
- formal 가설이 generated 증거로 binding되면 **비순환성 입증**을 요구한다(Guard 2).

---

## 3. 결과 (믿음 상태)

### H1 질량보존 — SUPPORTED (formal, 목표1)
`tests/test_plant_module.py::test_mass_conservation_source_is_root_uptake`가
`dmass == src` (rel 1e-6/abs 1e-9)로 통과. 잔차는 적분기 **밖**의 독립 회계라
(비순환성 입증) 솔버가 보존을 가정한 게 아니다. 임계규칙(< 1e-5) 자동 충족.

### H2 음이온 배제 — SUPPORTED (formal, 목표1)
`N = z·E_m·F/RT = +4.67` (E_m=−120 mV, z=−1) → `e^N = 106.8`. Tier-0 입력의
**닫힌 형식 귀결**이지 적합값이 아니며, 운반체(Vmax/Km)가 이 배제를 극복해야
한다는 점이 "수동 확산만으로 불충분"이라는 규칙을 충족.

### H3 Yamazaki 표본외 예측검증 — REFUTED (empirical, 목표4)
`reproduce_demo.py`가 11개 동족체 BAF를 **log10 RMSE 0.029**로 재현하나, W2 적합은
**동족체당 ~3 파라미터를 3 관측에** 맞춘 **포화 적합**이다 → pred≈obs(0.49/0.49)는
*예측*이 아니라 *재현*의 서명. **루프 이어가기(iteration 2)**: 실제 사전적(a-priori)
예측 — 이론/QSPR monotone f_xy(적합 아님, `reproduce_demo.py --rec`) — 의 오차는
**log10 RMSE 0.837**(포화 0.029의 ~29배, 줄기 6~40배 빗나감)이다. 즉 per-congener로
맞추지 못하게 하면 모델은 조직 BAF를 **예측하지 못한다**. "표본외 예측검증" 가설은
**정량적으로 거부**된다(evi-yamazaki + evi-yamazaki-apriori, 둘 다 measured).

### H4 곡립 위해성 예측 — REFUTED (empirical, 목표4)
measured 비교에서 곡립을 구조적으로 **3–8× 과소예측**(Tang PFOA endosperm 0.11 vs
0.95, dw; 단위 보정 후에도 안 닫힘). 유일한 "일치"(Kim)는 곡립 단일점 in-sample
앵커로 L_Ph를 강제한 것. 위해성 평가가 의존하는 구획의 구조적 과소예측은 식이
위해성 평가를 뒷받침할 수 없다 → 거부. **정직한 음성 결과는 1급 시민이다.**

### H5 토양 결합 동족체 의존 — SUPPORTED (formal, 목표3)
`literature_params.koc`(Higgins–Luthy +0.55/CF2 QSPR) 라이브 계산: Koc(PFBA C4)≈2.1
vs Koc(PFDoDA C12)≈53985 L/kg → **4.4 log10 스프레드(~25000×)**. 토양 지연계수
R=1+ρKd/θ가 동족체에 극도로 의존하므로 단일 상수 Cwo로는 모든 동족체를 표현할 수
없다(단쇄 용탈/장쇄 완충). 스프레드는 플랜트 BAF가 아니라 토양 QSPR에서 나오므로
비순환적(임계규칙 > 2 log10 자동 충족). **주의**: 이는 *coupling rationale*의
formal 확인이지, HYDRUS 결합 자체가 현장 토양 데이터로 검증됐다는 뜻이 아니다.

### H6 SMILES read-across 재현 — SUPPORTED (formal, 목표5)
`tests/test_pfas_structure.py` **23/23 통과**: 곡선 동족체와 일치하는 canonical
SMILES가 `params/parameters.json`의 동일 Compound를 재구성(measured read-across).
구조→파라미터 매핑이 **충실함**을 입증(파라미터가 물리적으로 옳다는 게 아님).
**주의**: 신규(novel) 구조의 f_xy는 provisional(QSPR/보간)이라 이 주장은 **알려진
구조에만** 한정된다.

### H7 구조적 적합성 — SUPPORTED (empirical; 사용자 재정의의 핵심)
"표본외 예측"이 아니라 **"구조가 fitting으로 실험을 재현 가능한가"**를 묻되, 의미가
있으려면 **자유도>0의 제약 적합**이어야 한다(포화 W2는 33/33=DOF0이라 무의미).
`validation/structural_adequacy_fit.py`는 **사용자의 mechanistic ORYZA2000 biomass
(`oryza_growth`) + 측정 증산(`forcing_rice`)**으로 구동(예시 로지스틱 아님)하여 세
제약 시나리오를 적합(11 동족체 × 3 조직 = 33 관측):

| 시나리오 | DOF | root | straw | grain | overall |
|---|---|---|---|---|---|
| A f_xy + 전역 L_Ph + 전역 kappa_d | 20 | 0.45 | 0.18 | 0.52 | 0.41 |
| B + per-cong L_Ph (grain) | 10 | 0.45 | 0.21 | **0.36** | 0.35 |
| C + per-cong kappa_d (root) | 10 | **0.26** | 0.16 | 0.51 | **0.34** |

- **straw(전류) ~0.16–0.18** — 전 동족체(장쇄 포함) **재현**. 포화 아닌 진짜
  goodness-of-fit → **전류 구조(GHK 배제 + f_xy TSCF + 결합) 적합성 입증.**
- **root**: 전역 kappa_d 0.45 → **per-cong kappa_d 0.26**(동족체별 막투과로 개선).
- **grain**: 전역 L_Ph 0.52 → **per-cong L_Ph 0.36**(동족체별 phloem loading으로 개선,
  단 장쇄 상승의 잔여 floor 남음).
- **전체 ~0.34(BAF 평균 ~2배 이내) @ DOF 10.** placeholder biomass의 곡립 파국
  (0.987)은 **현실적 ORYZA2000 biomass에서 사라진다**(biomass 드라이버가 결정적).
- **결론: 당신의 ORYZA2000로 구동 시, 구조는 전 조직을 제약 적합으로 ~2배 이내
  재현(SUPPORTED)** — 사용자의 목표("fitting으로 구조 모사 입증") 달성. 남은 건
  grain 장쇄 상승의 잔여 floor뿐.

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

1. **모델의 기여는 메커니즘·계산 기반에 있다, 경험적 예측에 있지 않다 (아직).**
   음이온 배제·질량보존·동족체 분해 토양 sorption·SMILES read-across는 견고하다
   (formal, SUPPORTED). 그러나 자주 헤드라인으로 인용되는 "RMSE 0.029"는 검증이
   아니라 **포화 재현**이며, 엄밀성 관점에서 예측 주장의 근거가 되지 못한다.
2. **곡립은 구조적 한계**이고 정직한 음성으로 기록되어야 한다(이미 CLAUDE.md가
   인정하는 바를 독립 엔진이 확인).
3. **진짜 검증으로 가는 길**: (a) f_xy/L_Ph를 적합에 쓴 데이터와 **분리된** 표본외
   데이터셋으로 평가, (b) 곡립 구조 항(비가역/이력 흡착, fw/dw 단위) 재정식화,
   (c) HYDRUS 결합을 **현장 토양 데이터**로 검증(현재는 coupling rationale만 formal
   확인), (d) novel-structure f_xy QSPR를 측정으로 검증, (e) 데모를 결과 보고에서
   `synthetic_proxy`로 취급.

---

## 5. 산출물 (`sci_adk_review/`)

```
proposal.md                     # 동결된 4-pane 사전등록 (배경/목표/방법/기대산출물 전문)
build_review.py                 # 재현 드라이버 (Spec/Evidence/Verdict + 루프 + verify + trap)
runs/pfas-rice/
  spec.json                     # 동결 Spec v1 (6 가설; referent/규칙/비순환)
  evidence/*.json               # 분류된 증거 (generated/measured) + 문헌(prior-work)
  verdicts/*.json               # chief-over-N 판정 trail (H2/H3/H4/H6)
  claims/*.json                 # 믿음 상태 (SUPPORTED×4, REFUTED×2) + 이력
  paper/draft.tex               # 결정론적 논문 초안 스켈레톤
runs/pfas-rice-trap/
  VALIDITY_HALT.txt             # 거부 증거 (synthetic_proxy → HALT)
```

### 재현
```bash
pip install -e /path/to/sci-adk          # 또는 PYTHONPATH=sci-adk/src
pip install numpy scipy pytest rdkit     # 모델 실행/근거 수치용
python sci_adk_review/build_review.py
sci-adk verify sci_adk_review/runs/pfas-rice   # exit 0, 전 claim REPRODUCED
```

> 주의: `runs/`의 증거 수치는 기준 커밋의 모델 출력을 그대로 옮긴 것이고, 각 증거는
> `provenance.code_ref`(커밋)와 `environment`(실행 경로)를 기록한다. 모델 코드가
> 바뀌어 출력이 달라지면 증거를 **새 항목으로 append**해야 한다(Evidence는
> append-only). prior-work는 모델이 실제로 인용하는 문헌(Yamazaki/Tang/Kim/Brunetti
> …)을 LITERATURE 증거로 기록했다.

---

## 6. 판정을 받아 실행한 수정 (sci-adk = player가 아니라 referee, 행동은 agent가 이어감)

sci-adk는 판정만 내고 끝이 아니라, **그 판정을 받아 agent가 다음 실험·수정을 돌리고
belief를 갱신하는 루프**다. REFUTED 판정을 받아 다음을 실행했다.

1. **루프 이어가기 (예측 정량화)**: H3 REFUTED가 요구하는 실제 사전적 예측을 돌렸다.
   monotone f_xy(적합 아님, `reproduce_demo.py --rec`)의 **OOS log10 RMSE = 0.837**
   (포화 0.029의 ~29배, straw 6~40배). append-only Evidence `evi-yamazaki-apriori`로 기록.
2. **모델 개선 시도 (bounded, 정직한 음성)**: 오차가 straw 지배라 재배분-shoot
   모델(`nstem_leaf`)을 같은 monotone f_xy로 시도(`validation/apriori_prediction.py`):
   **OOS RMSE 0.987 → 0.951**(소폭 개선). 단쇄 straw는 가까워지나 **장쇄 straw/grain
   붕괴**(PFDoDA straw 0.35 vs 49.75)는 그대로 — hysteretic 고-B sorption 갭. Evidence
   `evi-yamazaki-improve`로 기록. **개선은 미미, REFUTED 유지.**
3. **문서 과대주장 교정**: 저장소는 "0.029=재현"은 이미 정직히 적었으나 **실제 예측오차
   숫자가 없었다.** `reproduce_demo.py`(footer)·`CLAUDE.md`§6·`docs/OVERVIEW_KR.md`에
   **a-priori RMSE ≈0.84/≈0.95**를 명시 추가(0.029은 in-sample이라고 못박음).
4. **엄밀성 루프 상시화**: `tests/test_sci_adk_rigor.py`(커밋된 run 재검증 + **과대주장
   가드**: 경험적 예측 claim이 SUPPORTED가 되면 실패), `.github/workflows/rigor.yml`(push/PR
   마다 가드 실행), `sci_adk_review/run_rigor.sh`(로컬 전체 재생성+verify). sci-adk 미설치
   시 graceful skip. → 앞으로 같은 과대주장은 **자동 차단**된다.
5. **목표 재정의 → 구조적 적합성 입증(H7), ORYZA2000 구동**: "예측 말고 fitting으로 구조가
   실험을 재현 가능함을 보이자"는 재정의에 따라, **사용자의 mechanistic ORYZA2000 biomass
   (`oryza_growth`) + 측정 증산**으로 자유도>0 **제약 적합** 3종(A DOF20 / B,C DOF10)을
   실제로 돌렸다(`validation/structural_adequacy_fit.py`). 결과: **shoot 전류 재현
   (straw ~0.16–0.18)**, root는 per-cong kappa_d로 0.26, grain은 per-cong L_Ph로 0.36,
   **전 조직 overall ~0.34(≈2배 이내) @ DOF 10**. placeholder의 곡립 파국(0.987)은 현실적
   biomass에서 소멸. `evi-adequacy`(measured; hyp-adequacy SUPPORTS + hyp-grain REFUTES)
   + verdict로 기록 → **H7 SUPPORTED.** (직전까지 적합이 예시 로지스틱 biomass를 쓰던
   것을 사용자 지적으로 ORYZA2000으로 교체.)
