# 인수인계 — 중성/약전해질 확장 (Phase 1–3 완료, P4부터)

> 브랜치 `claude/neutral-organic-compound-tutphi` · PR #54 (draft)
> 설계 문서: `docs/NEUTRAL_DPU_EXTENSION_DESIGN_KR.md` (rev.5)
> **Phase 1 / 1.5 / 2 / 3 완료 · 남은 것은 P4(문헌 조사)뿐**

---

## 1. 한 줄 상태

모델이 **PFAS 전용(영구 음이온)에서 화학종 스펙트럼 전체**로 확장됐고(phase 1–2), 이제 그 확장이
**토양(Koc) · 구조(SMILES) · app** 경로에까지 연결됐다(phase 3). 새 모듈을 만들지 않고 중성 DPU
항들을 **다시 켜는** 방식이라, PFAS 경로는 **비트 단위로 불변**이다(20 시나리오 × 172 float, 정확 `==`).

$$\underbrace{f_n=1,\;f_d=0}_{\text{중성}}\;\longleftrightarrow\;
\underbrace{0<f_n<1}_{\text{약산/약염기}}\;\longleftrightarrow\;
\underbrace{f_n=0,\;f_d=1}_{\text{PFAS}}$$

---

## 2. 먼저 알아야 할 것 — PFAS 불변성 검증 방법

이 확장의 **최우선 불변식**은 "PFAS 결과가 한 자리도 변하지 않는다"이다.
`np.allclose`가 아니라 **정확한 `==`** 로 검증했고, 매 Phase마다 재실행했다.

방법 (다음 세션에서도 그대로 재사용할 것):

```bash
mkdir -p /tmp/base && git archive <이전_커밋> | tar -x -C /tmp/base
python validation/pfas_bit_identity.py /tmp/base        # 종료코드 0 = 동일
```

**스크립트는 이제 커밋돼 있다** — `validation/pfas_bit_identity.py` (phase 3에서 추가).
매 세션 재작성하던 것을 레포에 넣었으므로 "방법을 기억"할 필요가 없다. 두 체크아웃에서 각각
`model_api.simulate()`를 subprocess로 돌려 **20개 시나리오 × 172개 float**(4종 congener ×
3개 `f_xy_source` + lipid_loading + flooded profile; BAF·straw·B_k·N·eN)를 `==`로 비교한다.
`np.allclose`가 아닌 이유: 허용오차는 이 가드가 잡아야 할 바로 그 종류의 drift(합 순서 변경,
0을 더했다 빼기, 기본값 변경)를 숨긴다.

**추가 가드**
- `python reproduce_demo.py` → log10 RMSE **0.029** (변하면 안 됨)
- `pytest` 전체 → **239 passed, 2 skipped** (phase 3 기준선 216에서 +23)
- `tests/test_neutral.py`의 `GOLDEN_PFAS` — 확장 이전에 뽑아 고정한 BAF 값

---

## 3. 완료된 것

| Phase | 내용 | 진입점 |
|---|---|---|
| 1 | speciation + 중성 코어 (Briggs 결합/TSCF, $P_n$ Fick 투과) | `model_api.simulate_neutral(logKow, ...)` |
| 1.5 | 약전해질 체관 이온트랩 | `simulate_neutral(pKa=…, ion_trap=…, phloem_pH=…)` |
| 2 | 식물–대기 교환 (흡수/휘발/침적) | `simulate_neutral(K_AW=…, air=AirInputs(...))` |
| 3 | 주변부 통합 (토양 Koc 분기 · SMILES 분기 · app 화합물종류 선택) | `simulate_neutral(cwo_profile=…)` · `simulate_from_smiles()` · Expert 사이드바 |

**핵심 구조 (`src/pfas_rice_plant_module_4pool_surf.py`)**

- `root_uptake()` = `j_n(f_n 가중, Fick)` + `j_ed(f_d 가중, GHK)` + `j_carr`.
  **GHK 인자는 이온 항에만** 걸린다. 중성은 `fd=0`이 GHK를 죽이고, PFAS는 `fn=0`이 중성항을 죽인다.
- `binding_factors()`에 `f_lip*K_lip` 추가 → Briggs $K_{PW}$의 일반형.
- `phloem_loading_factor()` = `(1-w)*L_Ph + w*Λ`, `w = Π/(1+Π)`.
- `air_exchange()` — `K_AW > 0` 게이트.

**신규 API/파일**
- `model_api.simulate_neutral()` → `simulate()`와 **동일한 dict 모양** (phase 3: `cwo_profile=`/`cwo_kw=`)
- `model_api.simulate_from_smiles(..., compound_class=, logKow=, pKa=)` → 클래스 자동 분기
- `pfas_structure`: `Descriptors.compound_class`/`.logKow_crippen`,
  `neutral_compound_from_smiles()`, `compound_from_smiles_auto()`
- `soil_hydrus.paddy_kd(logKow=)` / `inputs_from_hydrus(logKow=)`,
  `model_api.hydrus_drivers(logKow=)` / `default_k_leach(logKow=)` / `cwo_profile_series(logKow=)`
- `validation/pfas_bit_identity.py` — 커밋된 비트 동일성 가드 (§2)
- `tests/test_neutral_phase3.py` — 23개 (dispatch + direction + UI 계약)
- `literature_params`: `speciation`, `ion_trap_factor`, `neutral_pathway_ratio`,
  `briggs_klip/rcf/tscf`, `koc_neutral`, `neutral_compound`, `f_lip_from_fresh_weight`
- `validation/neutral_probe.py` — 구조 검증 6블록 (실행하면 전부 assert)
- `tests/test_neutral.py` — 28개

---

## 4. 구현 중 검출된 설계 오류 5건 (중요 — 다시 되돌리지 말 것)

전부 **설계 문서가 틀렸고 구현이 맞았던** 경우다. 각각 다음 단계가 무엇인지를 바꿨다.

| # | 틀렸던 설계 주장 | 실제 | 위치 |
|---|---|---|---|
| 1 | 중성 스위치 = `Environment.z=0` | 약산은 중성분자+음이온 **동시 존재** → 전역 `z`로 표현 불가. 스위치는 `(f_n,f_d)` | 설계 D7 |
| 2 | PFAS 극한에서 이온트랩 $\Lambda\to1$ | $\Lambda\to10^{\Delta pH}=6.31$. 꺼짐은 **열역학이 아니라 속도론** ($f_n\to0$) | §5.1 |
| 3 | Trapp의 root lipid 0.025 ≈ 우리 0.020 | **기준이 다름**(fw vs dw). 기계론적 $B_k$는 Briggs 기울기(0.770)는 재현하나 **계수는 12.3× 낮음** | §11.4 |
| 4 | straw 농도 피크 = TSCF 피크(1.78) | 농도 $\propto$ TSCF×결합이라 **3.5로 이동**. 원래 수용 기준 `[1.5,3.0]`은 **올바른 구현을 떨어뜨렸을 것** | §3-3 |
| 5 | $\rho_k$가 기체교환에서 실물량이 됨 | **정확히 소거**됨. `CLAUDE.md` §8은 **고칠 필요 없음** | §7.2 |

또한 구현 중 잡힌 실수 2건:
- `neutral_compound`의 `kappa_d=0`이 약산의 이온 경로를 통째로 막고 있었다
  → Trapp 관계식 $P_d=P_n/10^{3.5}$을 기본값으로.
- 체관 계수를 `__post_init__`에 캐시했더니 `calibration.py:89`의
  `setattr(model.cmpd, "L_Ph", …)` 피팅이 무시되어 **기존 테스트 3개가 실패**
  → RHS마다 재계산으로 되돌림. `test_phloem_factor_follows_late_mutation_of_L_Ph`로 고정.

---

## 5. Phase 3 완료 (이 세션) · 남은 작업

설계 문서 §6 "Phase 3 — 주변부 통합" 3항목 **전부 완료**. Phase 3은 **물리를 추가하지 않는다** —
phase 1–2가 만들어 둔 중성 파라미터화를, 아직 PFAS 전용이던 세 경로에 **연결**했을 뿐이다.

### 5.1 토양 $K_{oc}$ 중성 분기 ✅
`koc_neutral`(Karickhoff $0.41K_{ow}$)은 구현돼 있었으나 **아무도 호출하지 않았다**.
`paddy_kd` · `inputs_from_hydrus` · `hydrus_drivers` · `default_k_leach` · `cwo_profile_series`가
`logKow=`를 받아 중성 분기로 가고, `simulate_neutral`은 `simulate`와 **같은 `cwo_profile` 스위치**를 갖는다.

노출 형상이 PFAS와 **같은 물리**에서 나온다: 극성 중성물질은 유출되고(logKow 1 → `k_leach` 0.034,
감소비 0.04) 소수성은 완충된다(logKow 6 → `k_leach` 0, 비 1.00) — PFAS에서 사슬길이가 하는 일 그대로.

> ⚠ **한계**: `params/cwo_kleach.csv`는 **PFAS 13종 전용 HYDRUS 보정**이다. 중성은 그 표가 아니라
> `k_leach(log10 Koc)` **회귀만** 재사용하므로 PFAS의 Koc 범위 밖에서는 **외삽**이다. **방향**은
> 의미가 있고 **속도(rate)는 provisional** — docstring에 명시.

### 5.2 SMILES 중성/약전해질 분기 ✅
`Descriptors.compound_class`가 **강산 머리기 AND 과불소 골격**을 **둘 다** 요구한다.
머리기만으로는 안 된다 — **2,4-D는 head_group='carboxylate'인데 과불소 탄소가 0개**인 평범한
약산이라, 머리기 검사만 했다면 영구 음이온으로 잘못 분류했을 것이다.

그 아래로 `neutral_compound_from_smiles`(Crippen `MolLogP` → Briggs)와,
클래스로 분기하는 `simulate_from_smiles`가 붙었다. 예전에는 중성종을 "가정 위반"으로 **flag만 하고
갈 곳이 없었다**.

> **중성 쪽이 오히려 쉽다**: 중성 파라미터 전체가 **`log K_ow` 단 하나**의 함수이고 RDKit이 그것을
> 구조에서 바로 준다. PFAS에는 대응물이 없다($K_{ow}$가 정의부터 애매) — 그래서 PFAS 쪽은
> read-across + fragment QSPR을 따로 만들어야 했다.
> **대가**: 같은 숫자 하나가 `K_lip` · TSCF · `K_oc`를 **동시에** 몰기 때문에 Crippen 추정오차가
> 전부에 전파된다 → `provisional=True`. **측정 logKow를 주면 플래그가 꺼진다.**
> `pKa`는 **예측하지 않는다**(RDKit에 신뢰할 모델이 없음) → 미지정 = 중성 취급이며, 그 사실을 말한다.

### 5.3 app 화합물 종류 선택 ✅
Expert 사이드바 §2가 "PFAS compound" → **"Compound"**, 맨 위에 클래스 라디오
(PFAS / 중성 / 약전해질). 중성 클래스는 logKow(슬라이더 또는 SMILES), 약전해질은 `pKa`+산/염기,
휘발성은 `K_AW`.

**Simple(일반인·한국어) 모드는 PFAS 전용으로 유지** — 중성 분기는 벼 검증 데이터가 없으므로
일반인용 평문 수치를 낼 근거가 없다.

중성에서 깨지거나 오해를 부르던 4곳을 처리: 헤더의 `p['n_C']`/`p['group']`(중성 params에 없음),
$e^N$ 카드(f_d=0이면 GHK 항이 항등적으로 0 → 중성분율 표시, 약전해질에서는 이온 경로가 실제로
살아 있으므로 $e^N$ 복귀), 그리고 PFAS 전용 탭 3개(chain-length / compare / Bayesian inverse)는
**잘못된 그림을 그리는 대신 이유를 설명**한다.

> 함정 하나: 그 가드를 `st.stop()`으로 쓰면 **탭이 아니라 페이지 전체가 멈춘다**(탭은 eager 렌더).
> `if/else`로 쓸 것.

### 5.4 부수 변경 1건 (되돌리지 말 것)
**퍼플루오로알킬 술폰아미드(FOSA)가 이제 `organic`으로 분류된다.** pKa~6으로 **영구 음이온이 아니고**,
기존 코드도 노트로는 그렇게 말하면서 **PFAS 분기로 돌리고 있었다** — phase 3이 갈 곳을 만들어 줬다.
기존 `test_pfas_structure::test_sulfonamide_speciation_warning`(플래그 검사)은 **수정 없이 통과**하고,
새 라우팅은 `test_sulfonamide_now_has_somewhere_to_go`가 고정한다.

### 5.5 남은 작업 — **P4 문헌 조사** (유일한 잔여 항목)
설계 문서 §8.2. **벼 조직별 중성/약전해질 시계열 데이터**를 찾는 일이며, 코드 작업이 아니다.
탐색 대상: 벼 농약 잔류 시험(조직별 분포), 논 제초제 흡수 실험, **하수 재이용 논의 의약품
(약전해질) 흡수 연구** — 약전해질 쪽이 데이터가 더 있을 수 있다.
**이것이 없으면 중성 분기는 영원히 "구조 검증"에 머문다**(§6).

부수적으로 열린 것(우선순위 낮음): 중성/에터 **soil Koc 측정값 부재**(`KOC_ETHER_LOG_OFFSET=0`은
명시적 GAP), 줄기 비표면적 $S$ 실측, 잎 갈락토지질의 Briggs 적용성(R6).

## 6. 정직한 한계 (보고할 때 반드시 병기)

- **벼 조직별 중성·약전해질 시계열 데이터가 `docs/literature_db/`에 없다.**
  PFAS 쪽 Yamazaki/Tang/Kim/Li에 해당하는 것이 없으므로 지금까지의 결과는
  **이론 대비 구조 검증이지 예측 검증이 아니다.** 레포의 "재현 ≠ 예측" 원칙을 그대로 적용할 것.
- **낟알**: 휘발성 화합물은 이제 보고 가능(휘발이 소실항). 하지만 **비휘발성 중성물질의 낟알은
  여전히 소실항이 대사뿐이라 과대**다.
- **지질 분율**: 현미만 measured, 뿌리/줄기/잎은 문헌 추정치. 잎은 갈락토지질 우세라
  Briggs의 옥탄올 유사성 가정을 벗어날 수 있다(리스크 R6).
- **줄기 $S$=2.7 m²/kg**: 원기둥 기하 유도값이고 실측이 아니다(PROVISIONAL).
- **CI는 전체 테스트를 돌리지 않는다.** `rigor.yml`은 `tests/test_sci_adk_rigor.py`(과잉주장 가드)만
  실행한다. 239 passed는 **로컬 수치**다.
- **(phase 3) 중성 `k_leach`는 외삽이다.** `params/cwo_kleach.csv`는 PFAS 13종 HYDRUS 보정이고
  중성은 `k_leach(log10 Koc)` 회귀만 재사용한다 → **방향은 의미 있고 속도는 provisional**(§5.1).
- **(phase 3) Crippen `MolLogP`는 추정이다.** 원자기여 추정치(보통 ~0.5–1 log 오차, 극성 헤테로고리는
  더 나쁨)이고, 그 하나가 `K_lip`·TSCF·`K_oc`를 동시에 몰기 때문에 오차가 전부에 전파된다 →
  `provisional=True`. **측정 logKow가 있으면 반드시 넣을 것**(§5.2).
- **(phase 3) `pKa`는 예측하지 않는다.** 미지정이면 중성으로 취급하며, 그것은 **가정이지 추론이 아니다**.

---

## 7. 다음 세션 재개 프롬프트 (그대로 붙여넣기)

```
docs/HANDOFF_neutral_extension.md 를 읽고 중성/약전해질 확장의 P4(문헌 조사)를 진행해줘.
Phase 1/1.5/2/3 은 전부 완료됐고 코드 작업은 남아 있지 않아 — 남은 것은 벼 조직별
중성·약전해질 시계열 데이터를 찾아 docs/literature_db/ 에 넣는 일이고, 그게 없으면
중성 분기는 계속 "구조 검증"에 머문다(예측 검증 주장 금지).
PFAS 경로는 비트 단위로 불변이어야 하고, 그 가드는 이제 커밋돼 있다:
  mkdir -p /tmp/base && git archive <이전_커밋> | tar -x -C /tmp/base
  python validation/pfas_bit_identity.py /tmp/base
시작 전 기준선: reproduce_demo RMSE 0.029, pytest 239 passed 2 skipped.
```
