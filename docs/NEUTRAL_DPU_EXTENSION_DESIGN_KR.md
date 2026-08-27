# 중성 유기화합물(neutral organic compound) 확장 — 설계 문서

> 상태: **DESIGN (구현 전)** · 대상 브랜치 `claude/neutral-organic-compound-tutphi`
> 선행 문서: `docs/dpu_model_summary_corrected.tex` (중성 DPU 원본),
> `docs/pfas_rice_compartmental_model.tex` (IOC/PFAS 확장), `docs/theory_anchor.tex` (Briggs/Trapp 대조표)

---

## 0. 한 줄 요약

현재 코어(`src/pfas_rice_plant_module_4pool_surf.py`)는 중성 DPU 위에 IOC 확장을 얹으면서 중성 항들을
**삭제한 것이 아니라 꺼둔** 상태다. 따라서 중성 지원은 *새 모델을 만드는 일*이 아니라
**(a) 이미 뚫려 있는 hook 3개를 켜고, (b) 코드에서만 빠진 항 2블록(중성 투과 · 기체 교환)을 복원하고,
(c) PFAS 전용 파라미터 경로를 우회하는 어댑터를 붙이는 일**이다.

---

## 1. 목적 / 범위 / 비목표

### 1.1 목적
중성 유기화합물(농약·의약품·산업용매 등, $f_n \approx 1$)에 대해 동일한 4-compartment 동적 ODE로
뿌리/줄기/잎/낟알 농도 시계열을 계산할 수 있게 한다. PFAS 경로는 **수치적으로 완전히 불변**이어야 한다.

### 1.2 범위 (in scope)
- 중성 분배: Briggs $K_{PW}$ 경로 (현 $B_k$의 특수화)
- 중성 막투과: $P_n$ 항 복원, GHK 전기화학항 비활성화, carrier 비활성화
- 중성 전류(translocation): Briggs TSCF bell → `f_xy`
- 기체 교환 3항: $\dot Q_\mathrm{GAS}$, $\dot Q_\mathrm{VOL}$, $\dot Q_\mathrm{DEP}$ (DPU 원본 §5.4–5.8)
- 대사 $\gamma_k > 0$ 활성화
- 토양측 $K_{oc}$ 중성 분기 (Karickhoff)
- SMILES → 중성 파라미터화 (RDKit `MolLogP`)

### 1.3 비목표 (out of scope, 이번 확장에서 제외)
- **약전해질(weak acid/base)**: $0 < f_n < 1$, pH 이온트랩. 중성($f_n\to1$)과 PFAS($f_n\to0$)의
  *중간* 케이스로, 두 극한이 모두 서면 자연스럽게 얻어지지만 이번 설계에는 포함하지 않는다.
  (구조는 열어 둔다 — §4 D2 참조)
- 다성분($N$종 동시) 텐서 형태 및 대사 변환행렬 $\Gamma$의 비대각 항 — 현 코어는 1종 스칼라.
- Method B(HYDRUS FORTRAN 개조).
- 중성 화합물에 대한 two-pool / nstem / lipid-loading 등 PFAS 전용 탐색 모델의 이식.

---

## 2. 왜 fork가 아니라 확장인가 — 이미 존재하는 hook

| hook | 위치 | 현재 상태 | 중성에서의 역할 |
|---|---|---|---|
| `Compound.fd` / `Compound.fn` | `4pool_surf.py:118-120` | `fn`은 **선언만 되고 미사용** | $f_n$이 중성 투과항의 계수 |
| `Environment.z` (원자가) | `:65-69` | `z=-1` 고정 | `z=0` → $N=0$ → GHK가 **자동으로 Fick 확산으로 축약** |
| `Compartment.gamma` (1차 대사) | `:132` | 구조 존재, PFAS라 `0` | 중성 유기물의 주 소실경로 |
| `root_uptake()` 중성항 주석 | `:224` | `# optional neutral passive term (negligible for PFAS): cmpd.fn * ...` | 구현 지점이 이미 표시됨 |

또한 **수식은 이미 문서에 있다**:
- `docs/pfas_rice_compartmental_model.tex:78` — $J_R$ 전체식에 중성 투과항
  $P_n\left(f_n^{o}C_w^{o}-f_n^{i}C_{w,1}\right)$ 이 명시되어 있고, `:90`에서
  "$f_n\to0$이면 소거"라고 적혀 있다. **코드에서만 생략된 것.**
- `docs/theory_anchor.tex:167` — Briggs 중성 관계식(RCF/TSCF) 전사 완료.
- `docs/dpu_model_summary_corrected.tex` — 기체교환 전 항의 완전한 행렬식.

즉 **이론 작업은 이미 끝나 있고, 코드 포팅만 남았다.**

---

## 3. 실현가능성 검증 (완료)

최소 패치 — `z=0` + `Vmax=0` + `K_PL ← a\,K_{ow}^{b}` + `f_xy ← TSCF` — 만으로 **기존 코어를 수정 없이**
돌린 결과 (측정 forcing `forcing_rice.Q_TP` + `growth_rice` biomass, $C_w^o=1$, 120일):

```
 logKow      K_lip   TSCF   B_root |     root    straw    grain
   0.50        3.0  0.401     0.91 |     0.78     0.91    55.42
   1.78       28.6  0.784     0.96 |     0.73     1.86    94.22
   3.00      249.1  0.426     1.40 |     1.20     2.59    59.73
   4.50     3559.3  0.038     8.02 |     7.90     2.05     5.48
   6.00    50858.1  0.001   102.62 |   102.59     0.10     0.02

z=0  j_R(Cwo=2, Cw=0.5) = 1.5     (정확히 Fick — 기대값과 일치)
z=-1 j_R(Cwo=2, Cw=0.5) = -2.269  (음이온 배제, e^N≈107)
```

**판정 3가지:**

1. **GHK → Fick 축약이 수치적으로 정확하다.** `z=0`에서 $N=0$, `_ghk_factor→1`, $e^N=1$이므로
   $j_{ed}=\kappa_d(C_w^o - C_{w,1})$. 별도 분기 없이 성립.
2. **중성 DPU의 특징 거동이 재현된다.** root는 Briggs RCF대로 $K_{ow}$에 단조 증가,
   **straw는 $\log K_{ow}\approx2$–$3$에서 종 모양(bell) 피크** — PFAS 경로에서는 나올 수 없는 형태.
3. **grain BAF 55–94는 비물리적이다.** $L_{Ph}=1$(중성이라 이온트랩 없이 체관액 ≈ 잎 자유농도)로 두면
   낟알이 소실항 없는 terminal accumulator가 된다. ⇒ **Phase 1만으로는 불충분하고,
   대사($\gamma$)·기체교환·$L_{Ph}(K_{ow})$가 함께 들어가야 한다**는 것이 이 검증의 핵심 결론.

> 검증 스크립트는 Phase 1에서 `validation/neutral_probe.py`로 커밋한다(부록 A에 전문 수록).

---

## 4. 설계 결정 (Design Decisions)

### D1. 모듈을 fork하지 않는다 — `speciation` 스위치로 흡수
`src/pfas_rice_plant_module_*.py`가 이미 6개다. 7번째(`_neutral.py`)를 만들면 향후 모든 수정이
2배가 된다. 대신 **`Compound`에 `speciation` 필드를 추가**하고, 중성 전용 항은
**기본값이 0이라 PFAS 경로에서 완전히 소거되는 가산항**으로 넣는다.

```python
@dataclass
class Compound:
    ...
    speciation: str = "anion"     # "anion" | "neutral"
    P_n: float = 0.0              # 중성 수동투과 a_R*P_n [L/(day kg)]  (0 → 항 소거)
    K_AW: float = 0.0             # 공기-물 분배 [-]      (0 → 기체교환 전체 소거)
```

**불변식(invariant): `P_n=0, K_AW=0, fn=0`이면 RHS가 현재와 부동소수점 수준까지 동일해야 한다.**

### D2. 약전해질을 위한 자리를 남긴다
$j_R$을 `j_neutral + j_electrodiffusion + j_carrier`의 **합**으로 쓰고 각각 $f_n$, $f_d$로 가중한다.
중성은 $(f_n,f_d)=(1,0)$, PFAS는 $(0,1)$. 약전해질은 $f_n+f_d=1$의 임의 값 — 즉
**이번 구현이 자동으로 약전해질의 절반을 완성한다**(나머지 절반은 체관 이온트랩).

### D3. 기체 교환은 `K_AW > 0`일 때만 켜지는 옵션 블록
현 `rhs()`에는 $\dot Q_\mathrm{GAS}/\dot Q_\mathrm{VOL}/\dot Q_\mathrm{DEP}$가 **아예 없다**
(가정 A3: `pfas_rice_compartmental_model.tex:63`). 이를 `rhs()` 말미에
`if self.cmpd.K_AW > 0:` 가드 하나로 감싼 블록으로 추가한다. 새 입력은
`AirInputs` 데이터클래스로 묶어 `RiceUptakeModel`의 **optional 필드**(`air: AirInputs | None = None`)로 둔다.

### D4. `Compartment`에 `f_lip`(총 지질 dw 분율)을 **새 필드로** 추가한다
현 `f_PL`은 *인지질* 분율(PFAS의 막결합 상대)이고, Briggs 지질항은 *총 지질*을 쓴다.
`f_PL`을 재해석해 재사용하면 의미가 충돌하고 `parameters.json`의
`tissue_composition_recommended`(root 0.015 / stem 0.005 / leaf 0.010 / grain 0.003)와도 어긋난다.
⇒ `f_lip: float = 0.0` 신규 필드 + `binding_factors()`에 `c.f_lip * cmpd.K_lip` 항 추가
(기본 0이므로 PFAS 불변).

> **결정 필요(열린 항목 O1)**: 벼 조직의 총 지질 dw 분율 문헌값 확보.
> 잠정 anchor: 뿌리 ~0.02, 줄기 ~0.01, 잎 ~0.03, 현미 ~0.025 (현미 지질 2–3%는 확립된 값).

### D5. 파라미터는 `parameters.json`에 넣지 않고 **별도 레코드로 주입**한다
`params/parameters.json`의 `congeners`는 PFAS 13종 전용 스키마
(`n_C`, `group`, `f_xy_recommended`, `f_xy_W2fit`, `f_xy_oryza`, `K_cw_wholecw_Lkg`, …)다.
중성 화합물을 여기 섞으면 `_transport_defaults()`(`model_api.py:442`)와
`build_parameters.py`가 오염된다.
⇒ **`model_api.simulate(record=…)` 주입 경로(이미 존재, `:304`)를 재사용**하고,
중성 레코드는 `literature_params.neutral_compound(logKow, …)`가 생성한다.
별도 파일 `params/neutral_compounds.json`(선택)에 예시 화합물만 둔다.

### D6. 기본 경로 불변 — 회귀 가드
다음이 **한 자리도 변하면 안 된다**:
- `reproduce_demo.py` log10 RMSE **0.029**
- `tests/` 전체 (현재 174 passed, 2 skipped)
- `simulate()`의 모든 기본 인자 동작
새 테스트 `tests/test_neutral.py`가 이를 명시적으로 assert 한다.

---

## 5. 수식 매핑표 — 중성 DPU ↔ 현재 IOC 코어

| 항 | 중성 DPU (원본) | 현재 코어 (PFAS) | 중성 확장에서 할 일 |
|---|---|---|---|
| 분배 | $K_{PW}=W+a K_{ow}^{b} L$ | $B_k=\theta_{fw}+(1-\theta_{fw})\sum_i f_i K_i$ | **수식 변경 불필요.** $K_{prot}=K_{cw}=0$, `f_lip·K_lip` 항 추가 → 정확히 Briggs |
| 뿌리 흡수 | $\dot q_\mathrm{up}=\dot U_\mathrm{HYD}/M_1$ (외부 BC) | GHK + Michaelis–Menten | $P_n(f_n^o C_w^o - f_n^i C_{w,1})$ 복원, `z=0`, `Vmax=0` |
| 전류(TSCF) | 암묵적 $C/K_{PW}$ (제한 없음) | `f_xy` (TSCF analog, PFAS 피팅) | `f_xy ← 0.784\exp[-(\log K_{ow}-1.78)^2/2.44]` |
| 체관 | **없음** (xylem-only) | $C_{Phl}=L_{Ph}C_{w,3}$ | $L_{Ph}(K_{ow})$ (Kleier 이동성). §7 Phase 3 |
| 기체 교환 | $\dot Q_\mathrm{GAS},\dot Q_\mathrm{VOL},\dot Q_\mathrm{DEP}$ | **없음** (A3: $K_{AW}\approx0$) | **포팅 (Phase 2, 최대 작업)** |
| 대사 | $\dot\Omega=\mathcal{G}C$ | `gamma` (=0) | $\gamma_k>0$ 활성화 |
| 성장 희석 | 로지스틱 $M(t)$ | ORYZA2000 / growth_rice | 그대로 재사용 |
| 낟알 형성 게이트 | 없음 | `grain_gate_` | 그대로 재사용 (DPU 정합) |

### 5.1 중성 QSPR 목록 (Phase 1에서 `literature_params`에 추가)

| 양 | 식 | 출처 | 비고 |
|---|---|---|---|
| $K_{lip}$ | $a K_{ow}^{b}$, $a=1.22$, $b=0.77$ | Briggs 1982 / Trapp 2004 | RCF $=0.82+0.03K_{ow}^{0.77}$에서 역산 |
| TSCF | $0.784\exp\!\left[-\frac{(\log K_{ow}-1.78)^2}{2.44}\right]$ | Briggs 1982 | bell; 최댓값 0.784 @ $\log K_{ow}=1.78$ |
| $P_n$ (막) | $P_d \approx P_n\cdot10^{-3.5}$ | Trapp 2000 | 이온 대비 중성이 $10^{3.5}$배 투과 |
| $K_{oc}$ | $0.41\,K_{ow}$ | Karickhoff 1981 | 토양측; PFAS 사슬길이 QSPR 대체 |
| $P_C$ (큐티클) | $10^{0.704\log K_{ow}-11.2}$ [m/s] | DPU 원본 식 (Pc) | Phase 2 |
| $\log K_{ow}$ | RDKit `Crippen.MolLogP` | — | SMILES 경로 (Phase 3) |

---

## 6. 파일별 변경 계획

### Phase 1 — 중성 코어 (비휘발성 중성물질까지 커버)

| 파일 | 변경 | 위험도 |
|---|---|---|
| `src/pfas_rice_plant_module_4pool_surf.py` | `Compound`에 `speciation`/`P_n`/`K_lip` 추가; `Compartment`에 `f_lip` 추가; `binding_factors()`에 지질항 1줄; `root_uptake()`에 중성항 1줄 | 낮음 (모두 기본값 0) |
| `src/literature_params.py` | `neutral_compound()`, `briggs_klip()`, `briggs_tscf()`, `koc_neutral()` 추가 | 없음 (순수 추가) |
| `src/model_api.py` | `simulate_neutral(logKow, …)` 래퍼 (내부적으로 `record=` 주입) | 낮음 |
| `validation/neutral_probe.py` | 신규 — §3 검증 재현 | 없음 |
| `tests/test_neutral.py` | 신규 — 회귀 가드 + Briggs bell 재현 | 없음 |
| `docs/` (본 문서) | 상태 갱신 | — |

### Phase 2 — 기체 교환 (중성 DPU의 나머지 절반)

| 파일 | 변경 | 위험도 |
|---|---|---|
| `src/pfas_rice_plant_module_4pool_surf.py` | `AirInputs` 신규; `permeabilities()` 신규(P_C/P_air/P_aqua/P_S/P_P/K_PA); `rhs()`에 `if K_AW>0` 가드 블록 | **중** — RHS 수정 |
| `src/model_api.py` | `AirInputs` 노출, 기본 대기 시나리오 | 낮음 |
| `params/` | 조직별 $S$(비표면적, 줄기 미정) · $\rho_k$ 확정 | — |
| `tests/test_neutral.py` | 휘발 질량수지 테스트 | — |

### Phase 3 — 주변부 통합

| 파일 | 변경 |
|---|---|
| `src/soil_hydrus.py` / `literature_params.koc` | 중성 $K_{oc}=0.41K_{ow}$ 분기; `params/cwo_kleach.csv`(PFAS 13종 전용)는 $K_{oc}$ 회귀 fallback 사용 |
| `src/pfas_structure.py` | `:213`의 "중성종 = 가정 위반 flag" 자리에 중성 분기 추가; RDKit `MolLogP` → Briggs 전 파라미터 |
| `src/model_api.py` | $L_{Ph}(K_{ow})$ (Kleier), `simulate_from_smiles` 중성 경로 |
| `app.py` | 화합물 종류 선택(PFAS / 중성) — Expert 모드만 |

---

## 7. 단계별 상세 + 수용 기준

### Phase 1 — 중성 코어

**구현 내용**
```python
# root_uptake() — 세 항의 합으로 재구성 (D2)
j_n    = cmpd.P_n * (cmpd.fn * Cwo - cmpd.fn * Cw_root)          # 신규 (중성 수동)
j_ed   = cmpd.kappa_d * g * (cmpd.fd * Cwo - cmpd.fd * eN * Cw_root)   # 기존
j_carr = Vmax_in*Cwo/(Km_in+Cwo) - Vmax_out*Cw/(Km_out+Cw)             # 기존
return j_n + j_ed + j_carr
```
```python
# binding_factors() — 지질항 추가
c.theta + (1-c.theta)*(c.f_prot*K_prot + c.f_PL*K_PL + c.f_cw*K_cw + c.f_lip*K_lip)
```

**수용 기준**
- [ ] `reproduce_demo.py` RMSE = 0.029 (변화 없음), 전체 테스트 통과
- [ ] `P_n=K_lip=f_lip=0`에서 RHS가 기존과 `np.allclose(rtol=0, atol=0)` 수준 동일
- [ ] $\log K_{ow}$ 스윕에서 straw TSCF bell 피크가 $\log K_{ow}\in[1.5,3.0]$에 존재
- [ ] root BAF가 $K_{ow}$에 단조 증가하고 Briggs RCF와 order-of-magnitude 일치
- [ ] **알려진 한계 명시**: 낟알 과대(§3-3) — Phase 2/3 전까지 grain 결과는 보고하지 않음

### Phase 2 — 기체 교환

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

**⚠ 이식 시 반드시 처리할 3가지 (gotcha)**

1. **단위계 충돌.** DPU 투과도 상수들은 **SI(m, s, g/mol)** 기준이고 우리 모델은 **day, L, kg, µg**이다.
   $P\,[\mathrm{m/s}] \times A\,[\mathrm{m^2}] \to \mathrm{m^3/s}$ → $\times 10^3 \times 86400$ 로 L/day 변환.
   $S$는 $\mathrm{m^2/kg}$, $\rho_k$는 kg/L. **전용 변환 헬퍼 하나로 격리하고 단위 테스트를 붙인다.**
2. **기호 충돌 2건.**
   - `RiceUptakeModel.phi` = *체관 재순환 분율*, DPU의 $\phi$ = *상대습도*. → 신규 필드는 **`RH`**로.
   - `Environment.z` = *원자가*, DPU $P_\mathrm{aqua}$의 $z$ = *확산 경로 길이*. → **`z_path`**로.
3. **$\rho_k$의 지위 변경.** 현재 $\rho_k$(`model_api.DEFAULT_TISSUE_DENSITY`)는 **보고용**이다
   (초기 draft의 차원오류 prefactor와 달리 수송 ODE에 들어가지 않음 — `CLAUDE.md` §8).
   기체교환 항에서는 $K_{PA}=B_k\rho_k/K_{AW}$로 **진짜 물리량이 된다.**
   → 문서/코드 주석에서 "밀도는 수송에 안 들어간다"는 서술을 **"중성 기체교환 항에서만 들어간다"**로 정정 필요.
4. **$S$(비표면적) 결측.** 현재 `Compartment.S`는 leaf=20.0, grain=2.0만 있고 **root/stem=0**이다.
   줄기 큐티클 교환을 켜려면 줄기 $S$가 필요 (DPU: 줄기는 큐티클 경로만, 뿌리는 휘발 없음).

**수용 기준**
- [ ] $K_{AW}=0$에서 Phase 1 결과와 완전 동일
- [ ] 단위 변환 테스트 (알려진 $P$, $A$ → L/day)
- [ ] 대기 청정($C_A=0$) + 휘발성 화합물 → 잎 농도가 휘발로 감소하고 **전체 질량수지가 닫힘**
- [ ] $C_A>0$ 단독 노출(토양 무오염)에서 잎이 오염되고 뿌리는 거의 안 됨

### Phase 3 — 주변부 통합
- [ ] 중성 $K_{oc}$로 HYDRUS 구동 → `cwo_profile="flooded"` 형상이 중성에도 성립
- [ ] SMILES(중성) → 파라미터 → ODE 완주, `provisional` 플래그 정확
- [ ] app Expert 모드에서 중성 화합물 선택 가능

---

## 8. 검증 전략과 데이터 공백

### 8.1 구조 검증 (즉시 가능)
- **Briggs 1982 barley**: RCF/TSCF 관계식 재현 — 모델이 아니라 QSPR 자체의 sanity check.
- **TSCF bell 위치**: $\log K_{ow}=1.78$ 부근 피크 (§7 Phase 1 수용 기준).
- **극한 일치**: $f_n\to0$이면 PFAS 결과와 정확히 일치 (D1 불변식).

### 8.2 데이터 검증 (제약 있음)
- **Trapp 1994 bromacil**: DPU 원본의 검증 세트. 다만 **작물이 벼가 아니다.**
- ⚠ **데이터 공백 (핵심)**: `docs/literature_db/`에 **벼 조직별(뿌리/줄기/잎/현미) 중성 유기물
  시계열 데이터가 없다.** PFAS 쪽 Yamazaki/Tang/Kim/Li에 해당하는 중성 데이터셋이 없으므로,
  **Phase 1–2는 "구조 검증"까지만 주장 가능하고 예측 검증(OOS)은 주장할 수 없다.**
  이는 PFAS 쪽에서 확립한 "재현 ≠ 예측" 원칙(`CLAUDE.md` §6)을 그대로 따른다.
- 후보 데이터 탐색 대상: 벼 농약 잔류 시험(이미다클로프리드·트리사이클라졸 등 조직별 분포),
  논 제초제 흡수 실험. → **별도 작업 항목 (문헌 조사)**

### 8.3 표기 규칙
중성 결과를 보고할 때 **반드시** 다음을 병기한다:
`speciation="neutral"`, $\log K_{ow}$, 사용한 QSPR 출처, 그리고
**"벼 중성물질 데이터 부재 → 구조 검증 단계"** 라는 단서.

---

## 9. 리스크 · 열린 질문

| ID | 항목 | 영향 | 대응 |
|---|---|---|---|
| R1 | `rhs()` 수정이 PFAS 경로를 미세하게 흔듦 | 높음 (RMSE 0.029 회귀) | D6 회귀 가드; 가산항 기본 0; RHS 동등성 테스트 |
| R2 | 단위계 혼입 (SI ↔ day/L/kg) | 높음 (조용한 오차) | 변환 헬퍼 격리 + 단위 테스트 |
| R3 | 벼 중성 검증 데이터 부재 | 중 | 주장 범위를 구조 검증으로 한정 (§8.3) |
| R4 | 낟알 $L_{Ph}(K_{ow})$ 근거 부족 | 중 | Phase 3까지 grain 결과 비보고 |
| O1 | 벼 조직 **총 지질** dw 분율 | — | 문헌 확보 필요 (D4) |
| O2 | 줄기 비표면적 $S$ | — | Phase 2 전 확보 |
| O3 | 대사속도 $\gamma_k$ 출처 | — | 화합물별 문헌; 기본 0 유지 |
| O4 | 약전해질을 언제 열 것인가 | — | 중성이 서면 재검토 (D2) |

---

## 10. 기호 대조표 (신규 도입분)

| 기호 | 코드 | 단위 | 의미 |
|---|---|---|---|
| $f_n$ | `Compound.fn` | – | 중성 분율 (중성=1, PFAS=0) |
| $P_n$ | `Compound.P_n` | L/(day·kg) | $a_R P_n$ — 중성 수동투과 (질량비) |
| $K_{lip}$ | `Compound.K_lip` | L/kg lipid | $a K_{ow}^{b}$ (Briggs 지질항) |
| $K_{AW}$ | `Compound.K_AW` | – | 공기–물 분배; **0이면 기체교환 전체 소거** |
| $L$ | `Compartment.f_lip` | kg/kg dw | 총 지질 분율 (인지질 `f_PL`과 별개) |
| $RH$ | `AirInputs.RH` | – | 상대습도 (DPU $\phi$; `model.phi`와 충돌 회피) |
| $C_A$ | `AirInputs.C_A` | µg/L | 대기 농도 |
| $f_p$ | `AirInputs.f_p` | – | 입자상 분율 |
| $v_\mathrm{DEP}$ | `AirInputs.v_dep` | m/s | 침적 속도 (~0.001) |
| $z_\mathrm{path}$ | `AirInputs.z_path` | m | 수상층 확산 경로 (Environment.`z`와 충돌 회피) |

---

## 부록 A. 검증 스크립트 (§3 재현)

Phase 1에서 `validation/neutral_probe.py`로 커밋 예정.

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

# NOTE: probe 단계에서는 f_PL 슬롯을 총지질로 임시 전용. 정식 구현에서는 f_lip 신규 필드(D4).
comps = [Compartment("root",  0.90, 0.07, 0.020, 0.50),
         Compartment("stem",  0.83, 0.05, 0.010, 0.72),
         Compartment("leaf",  0.78, 0.10, 0.030, 0.56, S=20.0),
         Compartment("grain", 0.14, 0.09, 0.025, 0.035, S=2.0)]

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

---

## 부록 B. 진행 체크리스트

- [ ] **P0** 설계 문서 검토·승인 ← *현재 단계*
- [ ] **P1** 중성 코어 (§7 Phase 1) + 회귀 가드
- [ ] **P1.5** O1(총 지질 분율) 문헌 확보
- [ ] **P2** 기체 교환 (§7 Phase 2) + 단위 테스트
- [ ] **P3** 토양 $K_{oc}$ · SMILES · app 통합
- [ ] **P4** 검증 데이터 문헌 조사 (§8.2) — 병행 가능
