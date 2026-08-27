# 중성 유기화합물(neutral) + 약전해질(weak electrolyte) 확장 — 설계 문서

> 상태: **Phase 1 + 1.5 + 2 구현 완료** · rev.5 · 브랜치 `claude/neutral-organic-compound-tutphi`
> 실행: `python validation/neutral_probe.py` · `pytest tests/test_neutral.py`
> API: `model_api.simulate_neutral(logKow, pKa=None, ...)` → `simulate()`와 동일한 dict
> 선행 문서: `docs/dpu_model_summary_corrected.tex` (중성 DPU 원본),
> `docs/pfas_rice_compartmental_model.tex` (IOC/PFAS 확장), `docs/theory_anchor.tex` (Briggs/Trapp 대조표)
>
> **rev.5 변경점 (Phase 2 — 기체 교환)**: `AirInputs`/`permeabilities()`/`air_exchange()`
> 착지 ($K_{AW}>0$ 게이트). **정정 1건**: gotcha ③($\rho_k$가 실물량이 된다)이 **틀렸음** —
> 정확히 소거되므로 `CLAUDE.md` §8은 손대지 않는다. **부수 효과**: 휘발이 낟알의 소실항이 되어
> §3-4의 비물리적 낟알이 해소된다(휘발성 화합물 한정). 줄기 $S$=2.7 배선. §7.2 참조.
>
> **rev.3 변경점 (Phase 1 구현)**: ① Phase 1 코드 착지 + 회귀 가드(부록 B).
> ② **정정 2건 — 구현 중 단언문이 설계 오류를 검출**: §5.1 이온트랩의 PFAS 극한
> ($\Lambda\to1$이 아니라 $10^{\Delta pH}$; 꺼짐은 속도론적), §11.4 Trapp의 fw 기준 지질
> (우리 dw 값과 같은 양이 아니며 12.3배 차이). ③ D9 `P_n` 기본값 결정.
>
> **rev.2 변경점**: ① 약전해질을 **비목표 → Phase 1 범위**로 편입 (사용자 결정).
> 이에 따라 "중성 스위치 = `z=0`"이 **`(f_n, f_d)` 가중합**으로 바뀜 (D7) — 세 케이스가 한 코드경로.
> ② 열린 항목 **O1(조직별 총 지질 분율) 해소** — 값은 이미 레포에 있었고, 문헌으로 교차검증 (§11).

---

## 0. 한 줄 요약

현재 코어(`src/pfas_rice_plant_module_4pool_surf.py`)는 중성 DPU 위에 IOC 확장을 얹으면서 중성 항들을
**삭제한 것이 아니라 꺼둔** 상태다. 따라서 이 작업은 *새 모델을 만드는 일*이 아니라
**(a) 이미 뚫려 있는 hook 3개를 켜고, (b) 코드에서만 빠진 항 2블록(중성 투과 · 기체 교환)을 복원하고,
(c) PFAS 전용 파라미터 경로를 우회하는 어댑터를 붙이는 일**이다.

약전해질을 함께 여는 결정에 따라 목표는 **하나의 스펙트럼을 전부 덮는 것**이 된다:

$$\underbrace{f_n=1,\;f_d=0}_{\text{중성}}\quad\longleftrightarrow\quad
\underbrace{0<f_n<1}_{\text{약산/약염기}}\quad\longleftrightarrow\quad
\underbrace{f_n=0,\;f_d=1}_{\text{PFAS (현재)}}$$

---

## 1. 목적 / 범위 / 비목표

### 1.1 목적
중성 유기화합물($f_n\approx1$)과 약전해질($0<f_n<1$)에 대해 동일한 4-compartment 동적 ODE로
뿌리/줄기/잎/낟알 농도 시계열을 계산할 수 있게 한다. PFAS 경로는 **수치적으로 완전히 불변**이어야 한다.

### 1.2 범위 (in scope)
- **화학종 분율(speciation)**: Henderson–Hasselbalch, **구획별 pH** → $(f_n, f_d)$
- 중성 분배: Briggs $K_{PW}$ 경로 (현 $B_k$의 특수화)
- 중성 막투과: $P_n$ 항 복원 (이온 GHK 항과 **병렬** 가산)
- 중성 전류(translocation): Briggs TSCF bell → `f_xy`
- **약산 체관 이온트랩**: $L_{Ph}$의 기계론적 대체 (Phase 1.5)
- 기체 교환 3항: $\dot Q_\mathrm{GAS}$, $\dot Q_\mathrm{VOL}$, $\dot Q_\mathrm{DEP}$ (DPU 원본 §5.4–5.8)
- 대사 $\gamma_k > 0$ 활성화
- 토양측 $K_{oc}$ 중성 분기 (Karickhoff)
- SMILES → 파라미터화 (RDKit `MolLogP`, `pKa`)

### 1.3 비목표 (out of scope)
- 다성분($N$종 동시) 텐서 형태 및 대사 변환행렬 $\Gamma$의 비대각 항 — 현 코어는 1종 스칼라.
- Method B(HYDRUS FORTRAN 개조).
- 중성/약전해질에 대한 two-pool / nstem / lipid-loading 등 PFAS 전용 탐색 모델의 이식.
- **양쪽성(zwitterion)·다가 이온** — $f_n+f_d=1$ 단일 해리 평형만 다룬다.
- 액포(vacuole) 격리 — Trapp 세포모델의 액포 구획은 넣지 않는다(4구획 유지).

---

## 2. 왜 fork가 아니라 확장인가 — 이미 존재하는 hook

| hook | 위치 | 현재 상태 | 확장에서의 역할 |
|---|---|---|---|
| `Compound.fd` / `Compound.fn` | `4pool_surf.py:118-120` | `fn`은 **선언만 되고 미사용** | **speciation의 주 스위치** (D7) |
| `Environment.z` (원자가) | `:65-69` | `z=-1` 고정 | 이온 분율에만 적용; 약염기는 $z=+1$ |
| `Compartment.gamma` (1차 대사) | `:132` | 구조 존재, PFAS라 `0` | 중성 유기물의 주 소실경로 |
| `root_uptake()` 중성항 주석 | `:224` | `# optional neutral passive term (negligible for PFAS): cmpd.fn * ...` | 구현 지점이 이미 표시됨 |
| `literature_params.f_d(pKa, pH)` | `literature_params.py:162` | PFAS 정당화용($f_d\ge0.94$)으로만 사용 | **약전해질의 speciation 엔진 그대로 재사용** |

또한 **수식은 이미 문서에 있다**:
- `docs/pfas_rice_compartmental_model.tex:78` — $J_R$ 전체식에 중성 투과항
  $P_n\left(f_n^{o}C_w^{o}-f_n^{i}C_{w,1}\right)$ 이 명시되어 있고, `:90`에서
  "$f_n\to0$이면 소거"라고 적혀 있다. **코드에서만 생략된 것.**
- `docs/theory_anchor.tex:167` — Briggs 중성 관계식(RCF/TSCF) 전사 완료.
- `docs/theory_anchor.tex:253-255` — 약산 이온트랩이 "중성 형태가 있어야 작동"하며 PFAS에서는
  꺼진다는 논증이 이미 서술되어 있다. **약전해질 확장은 이 문장을 코드로 옮기는 일이다.**
- `docs/dpu_model_summary_corrected.tex` — 기체교환 전 항의 완전한 행렬식.

즉 **이론 작업은 이미 끝나 있고, 코드 포팅만 남았다.**

---

## 3. 실현가능성 검증 (완료)

최소 패치 — `z=0` + `Vmax=0` + `K_PL ← a\,K_{ow}^{b}` + `f_xy ← TSCF` — 만으로 **기존 코어를 수정 없이**
돌린 결과 (측정 forcing `forcing_rice.Q_TP` + `growth_rice` biomass, $C_w^o=1$, 120일):

```
 logKow      K_lip   TSCF |     root    straw    grain
   0.00        1.2  0.166 |    0.833    0.517    32.06
   1.00       12.2  0.598 |    0.739    1.388    77.98
   1.78       28.6  0.784 |    0.731    2.087    94.50   <- TSCF 최대
   3.00      249.1  0.426 |    1.199    3.897    59.87
   3.75     1180.5  0.100 |    2.619    4.582    24.11   <- straw 최대
   4.50     3559.3  0.038 |    7.899    3.188     4.64
   6.00    50858.1  0.001 |  102.593    0.108     0.01

z=0  j_R(Cwo=2, Cw=0.5) = 1.5     (정확히 Fick — 기대값과 일치)
z=-1 j_R(Cwo=2, Cw=0.5) = -2.269  (음이온 배제, e^N≈107)
```
(§11 확정 지질값으로 재실행. 0.25 log 간격 전체 스윕 중 발췌.)

**판정 4가지:**

1. **GHK → Fick 축약이 수치적으로 정확하다.** `z=0`에서 $N=0$, `_ghk_factor→1`, $e^N=1$이므로
   $j_{ed}=\kappa_d(C_w^o - C_{w,1})$.
   > ⚠ 단 이 `z=0` 트릭은 **probe 전용 지름길**이다. $P_n$ 항이 생기고 약전해질이 들어오면
   > `z`를 건드리면 안 된다 — D7 참조.
2. **중성 DPU의 특징 거동이 재현된다.** root는 Briggs RCF대로 $K_{ow}$에 단조 증가,
   **straw는 종 모양(bell)** — PFAS 경로에서는 나올 수 없는 형태.
3. **⚠ straw 피크는 TSCF 피크와 위치가 다르다 (설계상 중요).**
   TSCF 자체는 $\log K_{ow}=1.78$에서 최대지만, **straw 농도는 $\log K_{ow}\approx3.75$에서 최대**다.
   조직 농도는 $C_\mathrm{straw}\propto \mathrm{TSCF}\cdot B$인데 $B$가 $K_{ow}$에 **단조 증가**하므로
   가우시안 × 증가함수의 곱이 피크를 오른쪽으로 민다. 해석적으로도 확인된다 —
   $B$의 지질항이 지배적일 때 $\frac{d}{dx}\!\left[-\frac{(x-1.78)^2}{2.44}+0.77\ln\!10\cdot x\right]=0$
   $\Rightarrow x = 1.78 + \frac{0.77\ln 10 \cdot 2.44}{2} \approx 3.94$ (관측 3.75; 저-$K_{ow}$에서
   $B$의 수분항이 지배해 실효 기울기가 0.77보다 작기 때문에 약간 왼쪽).
   ⇒ **"TSCF 피크 = 조직 농도 피크"로 수용 기준을 쓰면 틀린다.** §7 Phase 1 기준에 반영.
4. **grain BAF 32–94는 비물리적이다.** $L_{Ph}=1$로 두면 낟알이 소실항 없는 terminal accumulator가 된다.
   ⇒ **Phase 1만으로는 불충분하고, 대사($\gamma$)·기체교환·기계론적 $L_{Ph}$가 함께 들어가야 한다**는 것이
   이 검증의 핵심 결론이며, Phase 분할의 근거다.

> 검증 스크립트는 Phase 1에서 `validation/neutral_probe.py`로 커밋한다(부록 A에 전문 수록).

---

## 4. 설계 결정 (Design Decisions)

### D1. 모듈을 fork하지 않는다 — `speciation` 스위치로 흡수
`src/pfas_rice_plant_module_*.py`가 이미 6개다. 7번째(`_neutral.py`)를 만들면 향후 모든 수정이
2배가 된다. 대신 **중성 전용 항을 기본값 0이라 PFAS 경로에서 완전히 소거되는 가산항**으로 넣는다.

```python
@dataclass
class Compound:
    ...
    pKa: float | None = None      # None이면 speciation 계산 안 함 (하위호환)
    is_acid: bool = True          # 약염기는 False (z=+1)
    P_n: float = 0.0              # 중성 수동투과 a_R*P_n [L/(day kg)]  (0 → 항 소거)
    K_lip: float = 0.0            # Briggs 지질 분배 a*Kow^b [L/kg lipid] (0 → 항 소거)
    K_AW: float = 0.0             # 공기-물 분배 [-]      (0 → 기체교환 전체 소거)
```

**불변식(invariant): `P_n=0, K_lip=0, K_AW=0, fn=0`이면 RHS가 현재와 부동소수점 수준까지 동일해야 한다.**

### D2. 세 케이스를 하나의 코드경로로 (약전해질 포함)
$j_R$을 `j_neutral + j_electrodiffusion + j_carrier`의 **합**으로 쓰고 각각 $f_n$, $f_d$로 가중한다.

| 화합물 | $f_n$ | $f_d$ | 활성 항 |
|---|---|---|---|
| 중성 | 1 | 0 | $P_n$ (Fick)만 |
| 약산/약염기 | $0<f_n<1$ | $1-f_n$ | $P_n$ + GHK **둘 다** |
| PFAS | 0 | 1 | GHK + carrier (현재와 동일) |

### D7 (신규). **speciation 스위치는 `z`가 아니라 `(f_n, f_d)`다**
rev.1에서는 "중성 = `Environment.z=0`"으로 잡았으나, 약전해질을 열면 **틀린 설계**가 된다:
약산은 *중성 분자*(전위 무관)와 *음이온*(전위 영향)이 **동시에 존재**하므로 전역 `z` 하나로 표현할 수 없다.
올바른 구조는 **GHK 인자를 이온 항에만 적용**하는 것이다:

```python
j_n    = cmpd.P_n * cmpd.fn * (Cwo - Cw_root)                          # 전위 무관
j_ed   = cmpd.kappa_d * g * cmpd.fd * (Cwo - eN * Cw_root)             # 전위 의존 (GHK)
j_carr = Vmax_in*Cwo/(Km_in+Cwo) - Vmax_out*Cw_root/(Km_out+Cw_root)
return j_n + j_ed + j_carr
```
- 중성: `fd=0`이 GHK 항을 곱셈으로 죽인다 → `z` 값과 무관. **`z`를 건드릴 필요가 없다.**
- PFAS: `fn=0`이 중성 항을 죽인다 → 현재와 동일.

**부수 효과**: 원자가는 환경이 아니라 **화합물의 성질**이다(약염기는 $z=+1$).
`Environment.z`는 하위호환을 위해 유지하되, `Compound`가 제공하면 그쪽을 우선한다.
→ Phase 1에서 `Compound.z` 추가, `Environment.z`는 fallback.

### D8 (신규). 구획별 pH — 이온트랩의 전제
약전해질의 $f_n$은 **구획마다 다르다**(토양수 pH ≠ 세포질 pH ≠ 체관액 pH). Trapp 세포모델 표준값:

| 상(phase) | pH | 비고 |
|---|---|---|
| 토양 공극수 | 5.5–6.5 | `literature_params.PADDY_PH` 이미 존재 |
| 아포플라스트 / 물관액 | ~5.5 | 산성 |
| 세포질 | ~7.2 | |
| **체관액** | **~8.0** | **알칼리 → 약산이 갇힘(ion trap)** |

⇒ `Compartment.pH: float | None = None` 신규 필드 + `PHLOEM_PH = 8.0` 상수.
`None`이면 speciation 계산을 건너뛰어 현재 동작 유지.

### D3. 기체 교환은 `K_AW > 0`일 때만 켜지는 옵션 블록
현 `rhs()`에는 $\dot Q_\mathrm{GAS}/\dot Q_\mathrm{VOL}/\dot Q_\mathrm{DEP}$가 **아예 없다**
(가정 A3: `pfas_rice_compartmental_model.tex:63`). `rhs()` 말미에 `if self.cmpd.K_AW > 0:` 가드
하나로 감싼 블록으로 추가하고, 새 입력은 `AirInputs`로 묶어 `RiceUptakeModel`의
**optional 필드**(`air: AirInputs | None = None`)로 둔다.

### D4. `Compartment`에 `f_lip`(총 지질 dw 분율)을 **새 필드로** 추가한다
현 `f_PL`은 *인지질* 분율(PFAS의 막결합 상대)이고, Briggs 지질항은 *총 지질*을 쓴다.
**레포 DB가 이 구분을 이미 명시적으로 경고하고 있다** — `params/rice_tissue_params.csv`의 잎 행:

> `f_PL_membrane … CAUTION: phospholipid is a MINORITY of leaf membrane lipid; thylakoid
> galactolipids (MGDG/DGDG) dominate. Do NOT equate f_PL with total membrane lipid`

문헌도 이를 뒷받침한다(MGDG+DGDG가 틸라코이드 지질의 최대 80%). ⇒ `f_PL` 재활용은 **금지**,
`f_lip: float = 0.0` 신규 필드 + `binding_factors()`에 `c.f_lip * cmpd.K_lip` 항 추가
(기본 0이므로 PFAS 불변). **값은 §11에서 확정.**

### D5. 파라미터는 `parameters.json`에 넣지 않고 **별도 레코드로 주입**한다
`params/parameters.json`의 `congeners`는 PFAS 13종 전용 스키마
(`n_C`, `group`, `f_xy_recommended`, `f_xy_W2fit`, `f_xy_oryza`, `K_cw_wholecw_Lkg`, …)다.
중성 화합물을 여기 섞으면 `_transport_defaults()`(`model_api.py:442`)와 `build_parameters.py`가 오염된다.
⇒ **`model_api.simulate(record=…)` 주입 경로(이미 존재, `:304`)를 재사용**하고,
레코드는 `literature_params.neutral_compound(logKow, pKa=None, …)`가 생성한다.

### D6. 기본 경로 불변 — 회귀 가드
다음이 **한 자리도 변하면 안 된다**: `reproduce_demo.py` log10 RMSE **0.029**,
`tests/` 전체(Phase 1 착수 시점 실측 **188 passed, 2 skipped**), `simulate()`의 기본 인자 동작.
`tests/test_neutral.py`가 골든 BAF 값으로 이를 assert 한다.
**Phase 1 검증 결과**: 변경 전/후 코드를 각각 실행해 14개 시나리오 140개 float를
**정확한 `==`로 비교 → 전부 일치**. RMSE 0.029 유지, 전체 **204 passed, 2 skipped**.
**Phase 1.5 재검증**: 동일 비교 재실행 → 여전히 비트 동일. 전체 **210 passed, 2 skipped**.
**Phase 2 재검증**: 동일 비교 재실행 → 여전히 비트 동일. 전체 **216 passed, 2 skipped**.

### D9 (신규, Phase 1). `P_n` 기본값 = **1000** — 빠른 교환(평형) 극한
DPU base에는 뿌리 막 저항이 **아예 없고**($\dot q_\mathrm{up}$이 외부 BC, $K_{PW}$는 즉시 분배),
Briggs RCF도 마쇄 뿌리의 **평형** 상관식이다. 즉 base 모델은 $P_n\to\infty$ 극한이다.
$P_n=1000$이면 $\log K_{ow}$ 0.5–6 전 구간에서 $C_\mathrm{root}/B_\mathrm{root}\ge0.989$
(≥200이면 이미 수 % 이내)이고 stiffness/런타임 비용도 없다.
$P_n$은 **속도**이므로 낮추면 도달 속도만 느려지고 평형 목표값은 불변 — 테스트로 고정.

---

## 5. 수식 매핑표 — 중성 DPU ↔ 현재 IOC 코어

| 항 | 중성 DPU (원본) | 현재 코어 (PFAS) | 확장에서 할 일 |
|---|---|---|---|
| 화학종 분율 | 없음 ($f_n\equiv1$) | $f_d\equiv1$ | **H–H + 구획별 pH** (D8) |
| 분배 | $K_{PW}=W+a K_{ow}^{b} L$ | $B_k=\theta_{fw}+(1-\theta_{fw})\sum_i f_i K_i$ | **수식 변경 불필요.** `f_lip·K_lip` 항 추가 → 정확히 Briggs |
| 뿌리 흡수 | $\dot q_\mathrm{up}=\dot U_\mathrm{HYD}/M_1$ (외부 BC) | GHK + Michaelis–Menten | $P_n$ 항 복원, **GHK는 $f_d$에만** (D7) |
| 전류(TSCF) | 암묵적 $C/K_{PW}$ (제한 없음) | `f_xy` (PFAS 피팅) | `f_xy ← 0.784\exp[-(\log K_{ow}-1.78)^2/2.44]` |
| 체관 | **없음** (xylem-only) | $C_{Phl}=L_{Ph}C_{w,3}$ | **약산 이온트랩으로 $L_{Ph}$ 기계론화** (Phase 1.5) |
| 기체 교환 | $\dot Q_\mathrm{GAS},\dot Q_\mathrm{VOL},\dot Q_\mathrm{DEP}$ | **없음** (A3: $K_{AW}\approx0$) | **포팅 (Phase 2, 최대 작업)** |
| 대사 | $\dot\Omega=\mathcal{G}C$ | `gamma` (=0) | $\gamma_k>0$ 활성화 |
| 성장 희석 · 낟알 게이트 | 로지스틱 $M(t)$ | ORYZA2000 · `grain_gate_` | 그대로 재사용 |

### 5.1 QSPR / 관계식 목록

| 양 | 식 | 출처 | Phase |
|---|---|---|---|
| $f_n$ (약산) | $1/(1+10^{\,pH-pK_a})$ | Henderson–Hasselbalch | 1 |
| $f_n$ (약염기) | $1/(1+10^{\,pK_a-pH})$ | 〃 | 1 |
| $K_{lip}$ | $a K_{ow}^{b}$, $a=1.22$, $b=0.77$ | Briggs 1982 / Trapp 2004 | 1 |
| TSCF | $0.784\exp\!\left[-\frac{(\log K_{ow}-1.78)^2}{2.44}\right]$ | Briggs 1982 | 1 |
| $P_n$ (막) | $P_d \approx P_n\cdot10^{-3.5}$ | Trapp 2000 | 1 |
| **이온트랩 계수** | $\Lambda=\dfrac{1+10^{\,pH_\mathrm{ph}-pK_a}}{1+10^{\,pH_\mathrm{leaf}-pK_a}}$ | Trapp 2000 세포모델 | **1.5** |
| $K_{oc}$ | $0.41\,K_{ow}$ | Karickhoff 1981 | 3 |
| $P_C$ (큐티클) | $10^{0.704\log K_{ow}-11.2}$ [m/s] | DPU 원본 식 (Pc) | 2 |
| $\log K_{ow}$ | RDKit `Crippen.MolLogP` | — | 3 |

> **이온트랩 sanity**: $pK_a=4$, $pH_\mathrm{leaf}=7.2$, $pH_\mathrm{ph}=8.0$ ⇒ $\Lambda\approx6.3$.
> 즉 약산은 체관에서 ~6배 농축된다 — 2,4-D류 체관이동성 제초제의 고전적 설명.
>
> **⛔ rev.2의 오류 정정 (Phase 1 구현 중 단언문이 검출)**: rev.2는 "PFAS 극한에서 $\Lambda\to1$이라
> 트랩이 스스로 꺼진다"고 적었으나 **틀렸다.** $\Lambda$는 *중성종이 수송을 담당한다고 가정하고* 유도한
> **평형비**이므로, $pK_a\to-\infty$에서 양변이 완전 해리되어
> $\Lambda\to10^{\,pH_\mathrm{ph}-pH_\mathrm{leaf}}=10^{0.8}=6.31$ 로 간다 — 1이 아니다.
> PFAS에서 트랩이 꺼지는 것은 **열역학이 아니라 속도론** 때문이다: 그 평형비를 만들어 줄 flux를
> 중성종이 나르는데 $f_n\to0$이다. 실제 판별량은 **투과도 가중 비율**
> $\Pi = \dfrac{P_n f_n}{P_d f_d}$ (`literature_params.neutral_pathway_ratio`)이고,
> 잎 세포질 pH 7.2에서 $\Pi$ = **2.0 ($pK_a$ 4 제초제)** vs **2.0×10⁻⁷ (PFSA)** — **10⁷배** 붕괴한다.
> ⇒ Phase 1.5는 $L_{Ph}$에 $\Lambda$를 곱하는 것이 **아니라**, 체관 적재 컨덕턴스를 $f_n$에 비례시켜야
> 한다. 그래야 PFAS가 $\Lambda=6.3$의 가짜 농축을 받지 않는다.

---

## 6. 파일별 변경 계획

### Phase 1 — speciation + 중성 코어

| 파일 | 변경 | 위험도 |
|---|---|---|
| `src/pfas_rice_plant_module_4pool_surf.py` | `Compound`에 `pKa`/`is_acid`/`z`/`P_n`/`K_lip` 추가; `Compartment`에 `f_lip`/`pH` 추가; `binding_factors()` 지질항 1줄; `root_uptake()` D7 형태로 재구성 | 낮음 (모두 기본값 0/None) |
| `src/literature_params.py` | `speciation()`, `briggs_klip()`, `briggs_tscf()`, `neutral_compound()` 추가 | 없음 (순수 추가) |
| `src/model_api.py` | `simulate_neutral(logKow, pKa=None, …)` 래퍼 (`record=` 주입) | 낮음 |
| `params/rice_tissue_params.csv` | 변경 없음 — `total_lipid` 행 이미 존재 (§11) | 없음 |
| `validation/neutral_probe.py` | 신규 — §3 검증 재현 | 없음 |
| `tests/test_neutral.py` | 신규 — 회귀 가드 + Briggs bell + speciation 극한 | 없음 |

### Phase 1.5 — 약전해질 완성 (체관 이온트랩)

| 파일 | 변경 |
|---|---|
| `src/pfas_rice_plant_module_4pool_surf.py` | `rhs()`의 `C_Phl` 계산을 $L_{Ph}$ 상수 → $\Lambda(pK_a,pH)$ 기반으로 (PFAS는 $f_n\to0$이라 자동 소거) |
| `src/literature_params.py` | `ion_trap_factor(pKa, pH_leaf, pH_phloem)` |
| `tests/test_neutral.py` | $pK_a\to-\infty$(PFAS 극한)에서 $\Lambda\to1$ 확인 |

### Phase 2 — 기체 교환

| 파일 | 변경 | 위험도 |
|---|---|---|
| `src/pfas_rice_plant_module_4pool_surf.py` | `AirInputs` 신규; `permeabilities()` 신규; `rhs()`에 `if K_AW>0` 가드 블록 | **중** — RHS 수정 |
| `src/model_api.py` | `AirInputs` 노출, 기본 대기 시나리오 | 낮음 |
| `params/` | 조직별 $S$(줄기 미정) · $\rho_k$ 확정 | — |

### Phase 3 — 주변부 통합

| 파일 | 변경 |
|---|---|
| `src/soil_hydrus.py` / `literature_params.koc` | 중성 $K_{oc}=0.41K_{ow}$ 분기; `params/cwo_kleach.csv`(PFAS 전용)는 $K_{oc}$ 회귀 fallback |
| `src/pfas_structure.py` | `:213`의 "중성종 = 가정 위반 flag" 자리에 중성/약전해질 분기; `MolLogP` → Briggs |
| `app.py` | 화합물 종류 선택(PFAS / 중성 / 약전해질) — Expert 모드만 |

---

## 7. 단계별 수용 기준

### Phase 1
- [ ] `reproduce_demo.py` RMSE = 0.029 (변화 없음), 전체 테스트 통과
- [ ] `P_n=K_lip=f_lip=0, fn=0`에서 RHS가 기존과 완전 동일
- [ ] **D7 극한 테스트**: `fd=0` → `z`를 어떤 값으로 줘도 결과 불변 / `fn=0` → 현재 PFAS 결과와 동일
- [ ] **TSCF 함수 자체**의 최대가 $\log K_{ow}=1.78$ (QSPR sanity, ODE와 무관)
- [ ] **straw 농도**의 bell 피크가 $\log K_{ow}\in[3.0,4.5]$ (해석 예측 ~3.9, 관측 3.75 — §3-3).
      **주의: 1.78이 아니다.** 피크가 1.78 근처로 나오면 지질 결합항이 안 걸린 것이다
- [ ] root BAF가 $K_{ow}$에 단조 증가, Briggs RCF와 order-of-magnitude 일치
- [ ] **R6 민감도**: 잎 `f_lip` 0.04/0.055/0.07 스윕에서 straw 피크 위치가 ±0.5 log 이내로 안정
- [ ] **알려진 한계 명시**: 낟알 과대(§3-4) — Phase 1.5/2 전까지 grain 결과는 보고하지 않음

### Phase 1.5 — **DONE**
- [x] ~~$pK_a \le 0$에서 $\Lambda\to1$~~ — **틀린 기준. §5.1 정정 참조** ($\Lambda\to10^{\Delta pH}$)
- [x] 트랩 기여를 $\Pi=\frac{P_nf_n}{P_df_d}$ 게이트로 가중 →
      $L_{Ph}^\mathrm{eff}=(1-w)L_{Ph}+w\Lambda$, $w=\Pi/(1+\Pi)$
- [x] **PFAS는 구조적으로 소거**: `pKa`와 잎 `pH`가 **둘 다** 있어야 트랩이 켜지고,
      `simulate()`는 어느 쪽도 주지 않는다 ⇒ 수치적 우연이 아니라 **보증**
- [x] $pK_a=4$ 약산에서 $\Lambda=6.31$ 재현, 강산($pK_a=-3$)은 $L_{Ph}$ 그대로
- [x] PFAS 결과 불변 (RMSE 0.029, 골든값, 비트 동일 재확인)
- [x] $\kappa_d$ 기본값을 Trapp 관계식 $P_n/10^{3.5}$로 배선 (§7.1 발견 ②)
- [ ] pH 민감도 스윕 — Phase 2로 이월 (grain을 보고 가능하게 된 뒤)

### 7.1 Phase 1.5 구현 중 발견 3건

**① 트랩 효과는 영역 의존적이다 (설계 예측 밖).**
$L_{Ph}^\mathrm{eff}$를 4.2로 올려도 grain이 안 움직이는 구간이 있어 분리 실험을 했다
(`validation/neutral_probe.py` §5, $pK_a$ 4 · $\log K_{ow}$ 2.0):

| 기저 $L_{Ph}$ | $L_{Ph}^\mathrm{eff}$ | leaf(off→on) | grain(off→on) | 트랩 효과 |
|---|---|---|---|---|
| 1e-4 | 4.20 | 65.3 → 0.28 | 3.33 → 24.4 | **7.3×** |
| 1e-3 | 4.20 | 63.2 → 0.28 | 3.99 → 24.4 | **6.1×** |
| 1e-2 | 4.20 | 46.5 → 0.28 | 9.33 → 24.4 | **2.6×** |
| 1e-1 | 4.23 | 9.8 → 0.28 | 22.5 → 24.4 | 1.1× |
| 1.0 | 4.53 | 1.17 → 0.26 | 25.3 → 24.3 | 1.0× |

⇒ **트랩은 체관 *적재*가 병목일 때만 작동**한다. $L_{Ph}\to1$이면 잎이 이미 완전히 비워져
(leaf BAF 65→0.28) 추가 적재 용량이 무의미해지고, grain은 **잎 *공급* 한계**로 넘어간다.
PFAS의 피팅된 $L_{Ph}$는 1e-5~0.06 = **적재 한계 영역**이므로, $\Lambda$를 그냥 곱하는
설계였다면 바로 그 영역에서 최대 피해를 냈을 것이다.

**② `kappa_d=0`이 약산의 이온 경로를 통째로 막고 있었다.**
`neutral_compound`의 초기 기본값 `kappa_d=0`은 중성($f_d=0$)에는 무해하지만, 토양 pH 6.5에서
99.7%가 음이온인 $pK_a$ 4 약산에서는 **뿌리 흡수의 대부분을 조용히 삭제**한다.
문서화만 해두고 배선하지 않았던 Trapp 관계식 $P_d=P_n\cdot10^{-3.5}$을 실제 기본값으로 넣었다.

**③ 캐싱이 `calibration.py`를 깨뜨렸다 (기존 테스트가 검출).**
$L_{Ph}^\mathrm{eff}$를 `__post_init__`에 캐시했더니 `calibration.py:89`의
`setattr(model.cmpd, "L_Ph", …)`(이미 만들어진 model에 값을 주입해 피팅)가 무시되어
캘리브레이션 테스트 3개가 실패했다. RHS마다 재계산하도록 되돌렸다(PFAS 분기는 속성 검사 1회).
`test_phloem_factor_follows_late_mutation_of_L_Ph`로 고정.

### Phase 2
포팅할 식 (DPU 원본 §5.4–5.8):

$$P_{C}=10^{0.704\log K_{ow}-11.2},\quad
P_\mathrm{air}=\frac{\sqrt{300}\,K_{AW}}{200\,m^{0.5}},\quad
P_\mathrm{aqua}=\frac{D_{O_2}}{z_\mathrm{path}}\left(\frac{m}{32}\right)^{-1/2}$$

$$P_{C,\mathrm{tot}}=\left[P_C^{-1}+P_\mathrm{air}^{-1}+P_\mathrm{aqua}^{-1}\right]^{-1},\quad
P_S=\frac{\rho_w}{(1-RH)\,C_{H_2O,\mathrm{sat}}}(m/18)^{-1/2}K_{AW}\frac{\dot Q_\mathrm{XYL}}{A},\quad
P_P=P_S+P_{C,\mathrm{tot}}$$

$$K_{PA}=\frac{B_k\,\rho_k}{K_{AW}},\qquad
\dot Q_\mathrm{VOL}=\frac{A\rho_k}{M}\frac{P_P C}{K_{PA}},\qquad
\dot Q_\mathrm{GAS}=(1-f_p)C_A\frac{A}{M}P_P,\qquad
\dot Q_\mathrm{DEP}=v_\mathrm{DEP}f_p C_A\frac{A}{M}$$

**⚠ 이식 시 반드시 처리할 4가지 (gotcha)**

1. **단위계 충돌.** DPU 투과도 상수들은 **SI(m, s, g/mol)** 기준이고 우리 모델은 **day, L, kg, µg**이다.
   $P\,[\mathrm{m/s}] \times A\,[\mathrm{m^2}] \to \mathrm{m^3/s}$ → $\times 10^3 \times 86400$ 로 L/day 변환.
   **전용 변환 헬퍼 하나로 격리하고 단위 테스트를 붙인다.**
2. **기호 충돌 2건.**
   - `RiceUptakeModel.phi` = *체관 재순환 분율*, DPU의 $\phi$ = *상대습도* → 신규 필드는 **`RH`**.
   - `Environment.z` = *원자가*, DPU $P_\mathrm{aqua}$의 $z$ = *확산 경로 길이* → **`z_path`**.
3. ~~**$\rho_k$의 지위 변경.**~~ **⛔ 이 예측은 틀렸다 — §7.2 발견 ① 참조.**
   $\rho_k$는 기체교환 항에서도 **정확히 소거**된다. `CLAUDE.md` §8은 **정정할 필요가 없다.**
4. **$S$(비표면적) 결측.** 현재 `Compartment.S`는 leaf=20.0, grain=2.0만 있고 **root/stem=0**이다.
   줄기 큐티클 교환을 켜려면 줄기 $S$가 필요 (DPU: 줄기는 큐티클 경로만, 뿌리는 휘발 없음).

**수용 기준 — 전부 충족 (DONE)**
- [x] $K_{AW}=0$에서 Phase 1.5 결과와 **정확히 동일**(`==`), 대기 농도가 붙어 있어도 소거
- [x] 단위 변환을 **손계산과 대조** — 단일 상수 $10^3\times86400=8.64\times10^7$
- [x] 휘발성에서 잎/낟알 감소, **뿌리는 불변**(지하부 교환 없음)
- [x] $C_A>0$ 단독 노출($C_w^o=0$)에서 잎 124 / 낟알 131 vs **뿌리 0.014** µg/kg
- [x] PFAS 비트 동일 재확인, RMSE 0.029, 전체 **216 passed, 2 skipped**

### 7.2 Phase 2 구현 중 발견 2건

**① ⛔ $\rho_k$는 소거된다 — gotcha ③이 틀렸다.**
DPU 원본은 $\dot Q_\mathrm{VOL}=(A\rho/M)P_P C/K_{PA}$, $K_{PA}=K_{PW}\rho/K_{AW}$로 쓰는데,
prefactor의 $\rho$와 $K_{PA}$ 안의 $\rho$가 **정확히 상쇄**된다:

$$\frac{A\rho}{M}P_P\frac{C}{K_{PW}\rho/K_{AW}} = \frac{A}{M}P_P\,\frac{C\,K_{AW}}{K_{PW}}$$

수치 확인: $\rho$를 0.1→5.0 kg/L로 50배 바꿔도 $\dot Q_\mathrm{VOL}$이 **마지막 자리까지 동일**.
⇒ 조직 밀도는 이 모델의 수송에 **어디에도 들어가지 않는다**(기체교환 포함).
`CLAUDE.md` §8의 "밀도는 수송에 안 들어간다"는 **그대로 옳다** — 정정 불필요.
코드에는 `Compartment.rho`를 **두지 않았다**(있으면 쓰이는 것처럼 오해되므로);
밀도는 `model_api.DEFAULT_TISSUE_DENSITY`에 보고용으로만 남는다.

**② 기체교환이 낟알 terminal-sink 폭주를 해결한다 (부수 효과, §3-4의 답).**
낟알이 비물리적이었던 것은 **소실항이 없어서**였다. 휘발이 실제 소실항을 준다:

| $K_{AW}$ | root | stem | leaf | grain |
|---|---|---|---|---|
| 0 (off) | 0.974 | 0.713 | 5.85 | **141.98** |
| 1e-3 | 0.974 | 0.713 | 4.65 | 52.96 |
| 1e-2 | 0.973 | 0.712 | 0.218 | **0.229** |

⇒ **휘발성 중성물질에 대해서는 낟알을 보고할 수 있다.** 단 **비휘발성**($K_{AW}\approx0$)
중성물질의 낟알은 여전히 소실항이 없어 과대 — 대사 $\gamma>0$가 유일한 소실 경로다.

**③ 줄기 $S=0$ (gotcha ④)**: `TISSUE_SPECIFIC_AREA`에 2.7 m²/kg을 넣었다.
자유 파라미터가 아니라 **원기둥 기하** $4/(\rho d)$(직경 5 mm, 0.30 kg/L)에서 나온 값이며,
같은 기하가 기존 앵커를 재현한다(잎 판형 $2/(\rho t)$ 0.2 mm → 33 vs 코드 20;
낟알 2 mm → 1.7 vs 코드 2). 여전히 PROVISIONAL — 벼 줄기 비표면적 실측은 없다.

### Phase 3
- [ ] 중성 $K_{oc}$로 HYDRUS 구동 → `cwo_profile="flooded"` 형상이 중성에도 성립
- [ ] SMILES(중성/약산) → 파라미터 → ODE 완주, `provisional` 플래그 정확

---

## 8. 검증 전략과 데이터 공백

### 8.1 구조 검증 (즉시 가능)
- **Briggs 1982 barley**: RCF/TSCF 관계식 재현 — QSPR 자체의 sanity check.
- **TSCF bell 위치**: $\log K_{ow}=1.78$ 부근 피크.
- **극한 일치**: $f_n\to0$이면 PFAS 결과와 정확히 일치 (D1/D7 불변식).
- **이온트랩 극한**: $pK_a\to$ 매우 낮음 ⇒ $\Lambda\to1$ ⇒ PFAS의 "체관 트랩 없음"이 자동 재현.

### 8.2 데이터 검증 (제약 있음)
- **Trapp 1994 bromacil**: DPU 원본의 검증 세트. 다만 **작물이 벼가 아니다.**
- ⚠ **데이터 공백 (핵심)**: `docs/literature_db/`에 **벼 조직별(뿌리/줄기/잎/현미) 중성·약전해질
  유기물 시계열 데이터가 없다.** PFAS 쪽 Yamazaki/Tang/Kim/Li에 해당하는 것이 없으므로,
  **Phase 1–2는 "구조 검증"까지만 주장 가능하고 예측 검증(OOS)은 주장할 수 없다.**
  이는 PFAS 쪽에서 확립한 "재현 ≠ 예측" 원칙(`CLAUDE.md` §6)을 그대로 따른다.
- 후보 데이터 탐색 대상: 벼 농약 잔류 시험(조직별 분포), 논 제초제 흡수 실험,
  **하수 재이용 논의 의약품(약전해질) 흡수 연구** — 약전해질 쪽이 오히려 데이터가 더 있을 수 있다.
  → **별도 작업 항목 (P4 문헌 조사)**

### 8.3 표기 규칙
결과 보고 시 **반드시** 병기: `speciation`, $\log K_{ow}$, $pK_a$, 사용한 QSPR 출처,
그리고 **"벼 검증 데이터 부재 → 구조 검증 단계"** 라는 단서.

---

## 9. 리스크 · 열린 질문

| ID | 항목 | 영향 | 대응 |
|---|---|---|---|
| R1 | `rhs()` 수정이 PFAS 경로를 미세하게 흔듦 | 높음 (RMSE 0.029 회귀) | D6 회귀 가드; 가산항 기본 0; RHS 동등성 테스트 |
| R2 | 단위계 혼입 (SI ↔ day/L/kg) | 높음 (조용한 오차) | 변환 헬퍼 격리 + 단위 테스트 |
| R3 | 벼 중성/약전해질 검증 데이터 부재 | 중 | 주장 범위를 구조 검증으로 한정 (§8.3) |
| R4 | 낟알 $L_{Ph}$ 근거 부족 | 중 → **완화** | Phase 1.5의 이온트랩이 기계론적 근거 제공 |
| R5 | 구획별 pH가 실측이 아님 | 중 | Trapp 표준값 사용 + 민감도 분석 (Phase 1.5 수용 기준) |
| ~~O1~~ | ~~벼 조직 총 지질 분율~~ | — | **해소 — §11** |
| O2 | 줄기 비표면적 $S$ | — | Phase 2 전 확보 |
| O3 | 대사속도 $\gamma_k$ 출처 | — | 화합물별 문헌; 기본 0 유지 |
| ~~O4~~ | ~~약전해질을 언제 열 것인가~~ | — | **해소 — Phase 1/1.5로 편입 (사용자 결정)** |
| O5 (신규) | 약염기($z=+1$) 검증 사례 | 낮음 | 약산 먼저; 약염기는 부호만 뒤집힘 |

---

## 10. 기호 대조표 (신규 도입분)

| 기호 | 코드 | 단위 | 의미 |
|---|---|---|---|
| $f_n,\;f_d$ | `Compound.fn/.fd` | – | 중성/해리 분율 ($f_n+f_d=1$). **speciation의 주 스위치** |
| $pK_a$ | `Compound.pKa` | – | `None`이면 speciation 미적용(하위호환) |
| $z$ | `Compound.z` | – | 이온 원자가 (약산 −1, 약염기 +1); `Environment.z`는 fallback |
| $pH_k$ | `Compartment.pH` | – | 구획 pH (D8) |
| $\Lambda$ | `ion_trap_factor()` | – | 체관 이온트랩 농축 계수 |
| $P_n$ | `Compound.P_n` | L/(day·kg) | $a_R P_n$ — 중성 수동투과 (질량비) |
| $K_{lip}$ | `Compound.K_lip` | L/kg lipid | $a K_{ow}^{b}$ (Briggs 지질항) |
| $K_{AW}$ | `Compound.K_AW` | – | 공기–물 분배; **0이면 기체교환 전체 소거** |
| $L$ | `Compartment.f_lip` | kg/kg dw | 총 지질 분율 (인지질 `f_PL`과 별개) |
| $RH$ | `AirInputs.RH` | – | 상대습도 (DPU $\phi$; `model.phi`와 충돌 회피) |
| $C_A$ | `AirInputs.C_A` | µg/L | 대기 농도 |
| $f_p$ | `AirInputs.f_p` | – | 입자상 분율 |
| $v_\mathrm{DEP}$ | `AirInputs.v_dep` | m/s | 침적 속도 (~0.001) |
| $z_\mathrm{path}$ | `AirInputs.z_path` | m | 수상층 확산 경로 (원자가 `z`와 충돌 회피) |

---

## 11. O1 해소 — 조직별 총 지질 분율 `f_lip`

### 11.1 결론: 값은 이미 레포에 있었다

`params/rice_tissue_params.csv`에 **organ별 `total_lipid` 행이 이미 존재**한다
(`f_PL_membrane`와 **별도 행**으로 — D4의 설계 근거가 데이터 레벨에서 이미 확립되어 있었다).

| organ | low | **recommended** | high | evidence | source_key |
|---|---|---|---|---|---|
| root | 0.01 | **0.020** | 0.03 | estimate | `S_physiol` |
| stem | 0.01 | **0.015** | 0.02 | estimate | `S_physiol` |
| leaf | 0.04 | **0.055** | 0.07 | estimate | `S_physiol` |
| grain_brown | 0.027 | **0.030** | 0.036 | **measured** | `S_proximate_rg` |
| (husk) | 0.005 | 0.010 | 0.015 | estimate | `S_physiol` |
| (grain_white) | 0.003 | 0.006 | 0.009 | measured | `S_USDA_white` |

> **rev.1의 잠정 anchor는 틀렸다**: 잎을 0.03으로 적었으나 실제 DB 권장값은 **0.055**(~2×),
> 줄기 0.01→0.015, 현미 0.025→0.030. Phase 1은 **위 표의 recommended 값을 쓴다.**

### 11.2 문헌 교차검증

| 조직 | 문헌값 | DB 값 | 판정 |
|---|---|---|---|
| 현미 | 조지방 **3.04–3.59% dw**; 별도 보고 **2.75–4.49% dw**; n-헥산 추출 1.92–2.72% | 0.030 (0.027–0.036) | **일치 (measured 확인)** |
| 볏짚(줄기+잎) | 에테르추출물 **2.04% DM**; 친유성 추출물 **3.4% DM** | — | **아래 정합성 검사** |
| 잎 (유사체) | 라이그래스 EE **4.2% DM**, 알팔파 3.1%, 클로버 2.8%, 옥수수사일리지 2.7% (Palmquist & Jenkins 1980) | 0.055 (0.04–0.07) | 범위 상단이지만 **녹색 활성 엽신으로 타당** |
| 뿌리 | Trapp(2015) DPU 계열 기본 **root lipid 0.025** | 0.020 (0.01–0.03) | **⛔ 정합 아님 — 기준이 다름 (§11.4)** |

**정합성 검사 (독립적 교차검증)** — 볏짚은 줄기(대+엽초)와 엽신의 혼합이다.
질량비 70:30 가정 시 $0.7\times0.015+0.3\times0.055=\mathbf{0.027}$,
50:50 가정 시 $\mathbf{0.035}$. 측정된 볏짚 지질 **2.04–3.4% DM**이 이 구간을 정확히 감싼다.
⇒ **줄기 0.015와 잎 0.055는 서로 독립적인 볏짚 측정치와 모순되지 않는다.**

### 11.3 남은 공백 (정직한 기술)
- **뿌리·줄기·잎은 여전히 `estimate`다** — 벼 조직별 총 지질의 *직접 측정치*는 찾지 못했다.
  문헌은 압도적으로 곡립/미강(상업적 가치)에 편중되어 있고, 영양기관 지질은
  조성(지방산·갈락토지질 종류) 연구는 많으나 **dw 대비 총량**을 보고하지 않는다.
- 따라서 `f_lip`은 **현미만 measured, 나머지는 근연 매트릭스 유추(analogy)**다.
  민감도가 큰 것으로 드러나면(특히 잎, $K_{lip}$이 큰 친유성 화합물에서) 이것이 곧
  **다음 실험 우선순위**가 된다 — PFAS 쪽 `K_cw` 공백과 같은 성격.
- 지질 **조성**의 함의는 별개 리스크: 잎 지질의 대부분은 갈락토지질(MGDG/DGDG)로
  **인지질·중성지질과 극성이 다르다.** Briggs $a K_{ow}^b$는 옥탄올 유사상을 가정하므로
  잎에서 과대평가 가능. → §9 R6으로 추적.

| ID | 항목 | 대응 |
|---|---|---|
| R6 (신규) | 잎 지질이 갈락토지질 우세 → 옥탄올 유사성 가정 위배 가능 | Phase 1에서 잎 `f_lip` 민감도 스윕(0.04/0.055/0.07)을 수용 기준에 포함 |
| R7 (신규) | **측정 지질 vs Briggs "옥탄올 등가"의 12.3배 격차** | §11.4 — 두 경로를 모두 제공, 측정값을 기본으로 |

### 11.4 ⛔ 정정 — Trapp의 0.025는 **fw 기준**이고, 우리 값과 같은 양이 아니다

rev.2의 "Trapp 0.025 vs 우리 0.020 → 정합" 판정은 **틀렸다.** Phase 1 구현 중
`validation/neutral_probe.py`의 Briggs RCF 대조에서 드러났다.

Briggs의 경험식 $\mathrm{RCF}=0.82+0.03\,K_{ow}^{0.77}$을 Trapp은 $W + L\,a\,K_{ow}^{b}$로 읽으며,
여기서 $L\approx0.025$는 **fresh weight 기준**이다. 우리 `total_lipid`는 **dry weight 기준**이고
$B_k$ 식은 $(1-\theta_{fw})$로 dw→fw 변환을 한다. 즉 두 값은 **다른 기준의 다른 양**이다.

| | 지질항 계수 | 비고 |
|---|---|---|
| Briggs 경험식 | $0.03\;K_{ow}^{0.77}$ | fw 기준 "옥탄올 등가" |
| 우리 기계론식 | $(1-\theta)f_\mathrm{lip}^{dw}a = 0.00244\;K_{ow}^{0.77}$ | dw 측정값 |
| **비율** | **12.3×** | |

Briggs를 맞추려면 뿌리 지질이 **dw의 25%** 여야 하는데(= Trapp의 0.025 fw를 $\theta=0.90$에서
dw로 환산한 값), 그런 뿌리는 없다. ⇒ Briggs의 0.03은 **분석적으로 측정되는 지질이 아니라
세포벽·수베린·큐틴 흡착까지 흡수한 실효 흡착용량**이다.

**핵심은 무엇이 맞았는가다**: $K_{ow}$ **기울기는 0.770/log로 정확히 일치**한다(구조는 옳다).
어긋나는 것은 계수(절편)뿐이고, 그것은 정의·기준의 차이다.

**결정**: 두 경로를 모두 제공하되 **측정값을 기본**으로 한다 — 이 모델이 다른 모든 결합 pool을
다루는 방식(basis A, 측정 pool)과 일관되고, 어긋남을 숨기지 않는다.
Briggs 정합 실행이 필요하면 `literature_params.LIPID_OCT_EQUIV_FW` +
`f_lip_from_fresh_weight()`로 fw 기준 앵커를 쓸 수 있다.
이는 PFAS 쪽 `K_cw` 공백과 **같은 성격의 공백**이다 — 실측 흡착계수의 부재.

---

## 부록 A. 검증 스크립트 (§3 재현)

Phase 1에서 `validation/neutral_probe.py`로 커밋 예정.
**주의**: 아래는 $P_n$ 항이 없는 *현재* 코어에 대한 probe라 `z=0` 지름길을 쓴다.
정식 구현에서는 D7에 따라 `z`를 건드리지 않고 `fn=1, fd=0`으로 전환한다.

```python
"""현 코어가 z=0에서 중성 DPU base로 축약되는지 확인하는 probe."""
import sys, numpy as np
sys.path.insert(0, "src")
from pfas_rice_plant_module_4pool_surf import (
    Environment, Compound, Compartment, RiceUptakeModel, PlantInputs,
    binding_factors, root_uptake, ROOT, STEM, LEAF, FRUIT)
import forcing_rice as fr, growth_rice as gr

def briggs(logKow, a=1.22, b=0.77):
    Kow = 10.0 ** logKow
    return a * Kow ** b, 0.784 * np.exp(-(logKow - 1.78) ** 2 / 2.44)   # K_lip, TSCF

season, n = 120.0, 481
t = np.linspace(0, season, n)
Qtp = np.array([fr.Q_TP(x, season) for x in t])
b = gr.organ_biomass(t, season=season)
M = np.column_stack([b["root"], b["stem"], b["leaf"], b["grain"]])
inputs = PlantInputs(t=t, Cwo=np.ones(n), Qtp=Qtp, M=M)

# NOTE: probe 단계에서는 f_PL 슬롯을 총지질로 임시 전용하고 §11의 값을 쓴다.
# 정식 구현에서는 f_lip 신규 필드(D4).
comps = [Compartment("root",  0.90, 0.07, 0.020, 0.50),
         Compartment("stem",  0.83, 0.05, 0.015, 0.72),
         Compartment("leaf",  0.78, 0.10, 0.055, 0.56, S=20.0),
         Compartment("grain", 0.14, 0.09, 0.030, 0.035, S=2.0)]

for logKow in (0.5, 1.78, 3.0, 4.5, 6.0):
    K_lip, tscf = briggs(logKow)
    cmpd = Compound(name=f"neutral{logKow}", K_prot=0.0, K_PL=K_lip, K_cw=0.0,
                    kappa_d=50.0,                      # a_R*P_n (중성 투과, 이온의 ~10^3.5배)
                    Vmax_in=0.0, Km_in=1.0, Vmax_out=0.0, Km_out=1.0,
                    L_Ph=1.0, f_xy=tscf, fd=1.0, fn=0.0)
    m = RiceUptakeModel(env=Environment(z=0), cmpd=cmpd, comps=comps, inputs=inputs)
    Ce = m.solve(t).y[:, -1]
    Mf = inputs.M_(t[-1])
    straw = (Ce[STEM]*Mf[STEM] + Ce[LEAF]*Mf[LEAF]) / (Mf[STEM] + Mf[LEAF])
    print(f"{logKow:5.2f} {K_lip:10.1f} {tscf:6.3f} | "
          f"{Ce[ROOT]:8.2f} {straw:8.2f} {Ce[FRUIT]:8.2f}")

c0 = Compound("x", 0, 0, 0, kappa_d=1.0, Vmax_in=0, Km_in=1, Vmax_out=0, Km_out=1,
              L_Ph=1, f_xy=1)
assert abs(root_uptake(2.0, 0.5, c0, Environment(z=0)) - 1.5) < 1e-12   # 순수 Fick
```

> §3의 표는 위 스크립트(§11 확정 지질값)를 0.25 log 간격 전체 스윕으로 돌려 발췌한 것이다.
> Phase 1에서 정식 파일로 커밋할 때 스윕 범위를 `np.arange(0, 6.01, 0.25)`로 넓히고
> straw 피크 위치를 함께 출력하여 §7의 수용 기준을 자동 검증한다.

---

## 부록 B. 진행 체크리스트

- [x] **P0** 설계 문서 작성
- [x] **O1** 조직 총 지질 분율 조사 → §11 (레포 DB + 문헌 교차검증 + §11.4 정정)
- [x] **O4** 약전해질 범위 결정 → Phase 1/1.5 편입
- [x] **P1** speciation + 중성 코어 + 회귀 가드 — **DONE**
      (`4pool_surf`: `pKa/is_acid/z/P_n/K_lip/K_AW` + `f_lip/pH`, `root_uptake` D7 재구성;
      `literature_params`: `speciation/ion_trap_factor/neutral_pathway_ratio/briggs_*`
      `/koc_neutral/neutral_compound/f_lip_from_fresh_weight`;
      `model_api.simulate_neutral` + `_solve_and_package` 추출;
      `validation/neutral_probe.py`; `tests/test_neutral.py` 16개)
- [x] **P1.5** 약전해질 완성 — **DONE** (체관 트랩을 $\Pi$ 게이트로; `RiceUptakeModel.
      phloem_loading_factor()`, `simulate_neutral(ion_trap=…, phloem_pH=…, tissue_pH=…)`,
      `Compartment.pH`, 테스트 21→26개). 발견 3건은 §7.1.
- [x] **P2** 기체 교환 — **DONE** (`AirInputs`, `permeabilities`, `air_exchange`;
      `simulate_neutral(K_AW=…, air=…)`; 테스트 26→28개; 발견 2건은 §7.2)
- [ ] **P3** 토양 $K_{oc}$ · SMILES · app 통합 ← *다음*
- [ ] **P2** 기체 교환 + 단위 테스트
- [ ] **P3** 토양 $K_{oc}$ · SMILES · app 통합
- [ ] **P4** 검증 데이터 문헌 조사 (§8.2) — 병행 가능

---

## 참고문헌 (§11 조사분)

- FAO, *Rice in human nutrition* — 곡립 구조·조성. https://www.fao.org/4/t0567e/t0567e08.htm
- 현미 조지방 (품종별 proximate, % dw). https://www.researchgate.net/figure/Proximate-composition-of-rice-varieties-dry-weight-basis_tbl1_345127310
- 현미 지질 2.75–4.49% dw (질소·재식밀도 영향). https://pmc.ncbi.nlm.nih.gov/articles/PMC11604421/
- 볏짚 에테르추출물 2.04% DM (FAO, 이집트 반추동물 사료). https://openknowledge.fao.org/server/api/core/bitstreams/e6ccb174-d8f0-454a-9680-3b651c6fa76c/content/x5494e07.htm
- 볏짚 친유성 추출물 3.4% dw. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8981202/
- 조사료 EE (라이그래스 4.2%, 알팔파 3.1% 등; Palmquist & Jenkins 1980). https://www.aafco.org/wp-content/uploads/2023/01/fat_analysis_palmquist.pdf
- Trapp(2015) 계열 DPU 기본 root lipid 0.025 (ionizable PPCP generic model). https://academic.oup.com/etc/article/42/4/793/7730340
- 틸라코이드 갈락토지질 MGDG/DGDG가 총 지질의 최대 80% (D4·R6 근거). https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5897459/
- Brunetti et al. 2019, DPU module for HYDRUS. https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2019WR025432
