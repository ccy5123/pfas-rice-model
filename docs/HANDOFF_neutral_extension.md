# 인수인계 — 중성/약전해질 확장 (Phase 3부터)

> 브랜치 `claude/neutral-organic-compound-tutphi` · PR #54 (draft)
> 설계 문서: `docs/NEUTRAL_DPU_EXTENSION_DESIGN_KR.md` (rev.5)
> **Phase 1 / 1.5 / 2 완료 · Phase 3 미착수**

---

## 1. 한 줄 상태

모델이 **PFAS 전용(영구 음이온)에서 화학종 스펙트럼 전체**로 확장됐다. 새 모듈을 만들지 않고
중성 DPU 항들을 **다시 켜는** 방식이라, PFAS 경로는 **비트 단위로 불변**이다.

$$\underbrace{f_n=1,\;f_d=0}_{\text{중성}}\;\longleftrightarrow\;
\underbrace{0<f_n<1}_{\text{약산/약염기}}\;\longleftrightarrow\;
\underbrace{f_n=0,\;f_d=1}_{\text{PFAS}}$$

---

## 2. 먼저 알아야 할 것 — PFAS 불변성 검증 방법

이 확장의 **최우선 불변식**은 "PFAS 결과가 한 자리도 변하지 않는다"이다.
`np.allclose`가 아니라 **정확한 `==`** 로 검증했고, 매 Phase마다 재실행했다.

방법 (다음 세션에서도 그대로 재사용할 것):

```bash
SP=/tmp/scratch                     # 아무 임시 경로
mkdir -p $SP/base && git archive <이전_커밋> | tar -x -C $SP/base
python $SP/bitcmp.py $SP/base /home/user/pfas-rice-model
```

`bitcmp.py`는 두 체크아웃에서 각각 `model_api.simulate()`를 돌려
**14개 시나리오 × 140개 float**(4종 congener × 3개 `f_xy_source`, lipid_loading, flooded profile;
BAF·straw·B_k·N·eN)를 `==`로 비교한다. 스크립트는 커밋돼 있지 않으므로 필요하면 재작성한다
(로직은 단순: 두 경로에서 같은 러너를 subprocess로 실행 → JSON → 정확 비교).

**추가 가드**
- `python reproduce_demo.py` → log10 RMSE **0.029** (변하면 안 됨)
- `pytest` 전체 → **216 passed, 2 skipped** (착수 시점 188에서 +28)
- `tests/test_neutral.py`의 `GOLDEN_PFAS` — 확장 이전에 뽑아 고정한 BAF 값

---

## 3. 완료된 것

| Phase | 내용 | 진입점 |
|---|---|---|
| 1 | speciation + 중성 코어 (Briggs 결합/TSCF, $P_n$ Fick 투과) | `model_api.simulate_neutral(logKow, ...)` |
| 1.5 | 약전해질 체관 이온트랩 | `simulate_neutral(pKa=…, ion_trap=…, phloem_pH=…)` |
| 2 | 식물–대기 교환 (흡수/휘발/침적) | `simulate_neutral(K_AW=…, air=AirInputs(...))` |

**핵심 구조 (`src/pfas_rice_plant_module_4pool_surf.py`)**

- `root_uptake()` = `j_n(f_n 가중, Fick)` + `j_ed(f_d 가중, GHK)` + `j_carr`.
  **GHK 인자는 이온 항에만** 걸린다. 중성은 `fd=0`이 GHK를 죽이고, PFAS는 `fn=0`이 중성항을 죽인다.
- `binding_factors()`에 `f_lip*K_lip` 추가 → Briggs $K_{PW}$의 일반형.
- `phloem_loading_factor()` = `(1-w)*L_Ph + w*Λ`, `w = Π/(1+Π)`.
- `air_exchange()` — `K_AW > 0` 게이트.

**신규 API/파일**
- `model_api.simulate_neutral()` → `simulate()`와 **동일한 dict 모양**
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

## 5. 남은 작업 — Phase 3 (다음 세션의 시작점)

설계 문서 §6 "Phase 3 — 주변부 통합"이 그대로 To-Do다.

1. **토양 $K_{oc}$ 중성 분기**
   `literature_params.koc_neutral()` (Karickhoff $K_{oc}=0.41K_{ow}$)은 **이미 구현돼 있으나
   `soil_hydrus`/`cwo_profile`에 배선되지 않았다.** 현재 그 경로들은 PFAS 사슬길이 QSPR 전용이다.
   - `params/cwo_kleach.csv`(HYDRUS 보정)도 PFAS 13종 전용 → 중성은 $K_{oc}$ 회귀 fallback 필요.
   - 검증: `simulate_neutral(cwo_profile="flooded")`이 성립하는가.
2. **SMILES 경로**
   `src/pfas_structure.py:213`이 지금은 중성종을 "가정 위반"으로 **flag만** 한다.
   그 자리에 중성/약전해질 분기를 넣는다. RDKit `Crippen.MolLogP`가 $\log K_{ow}$를 바로 주므로
   **PFAS보다 오히려 쉽다**(PFAS는 $K_{ow}$가 정의부터 애매해 read-across/QSPR을 따로 만들어야 했다).
   `pKa`는 예측이 어려우므로 사용자 입력 또는 미지정(중성 취급)으로.
3. **app 통합** — Expert 모드에 화합물 종류(PFAS/중성/약전해질) 선택.
4. **P4 (병행 가능) — 검증 데이터 문헌 조사**: §6 참조.

---

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
  실행한다. 216 passed는 **로컬 수치**다.

---

## 7. 다음 세션 재개 프롬프트 (그대로 붙여넣기)

```
docs/HANDOFF_neutral_extension.md 와 docs/NEUTRAL_DPU_EXTENSION_DESIGN_KR.md 를 읽고
중성/약전해질 확장의 Phase 3 을 진행해줘. 브랜치는 claude/neutral-organic-compound-tutphi
(PR #54, draft). Phase 1/1.5/2 는 완료됐고 PFAS 경로는 비트 단위로 불변이어야 해 —
시작 전에 reproduce_demo RMSE 0.029 와 pytest 216 passed 를 기준선으로 먼저 확인해줘.
```
