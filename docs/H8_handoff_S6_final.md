# H8 핸드오프 — S6 종결 (α/QC1 · surface cross-field · W3 redox · Gap4)

**선행**: H7 (S6 본작업 W1·W2). 본 세션은 H7 §7 잔여작업을 **전부** 수행 → **S6 본작업 종결**.
**모델**: PFAS–논벼 4-compartment(root/stem/leaf/fruit-grain) mechanistic uptake (IOC 확장, DPU/Trapp). 12 congener: PFCA C4–C12(PFBA…PFDoDA) + PFSA C4/C6/C8(PFBS,PFHxS,PFOS). basis-A `B_k = θ_fw + (1−θ_fw)·Σf_dw·K` 확정.

---

## 0. TL;DR (한 줄 상태)

- **★ α/QC1 — 막 지배는 artifact가 아니라 실재(장쇄)**. basis-A서도 막분율 PFOA 68%→PFDoDA 98% 지속. Prompt3 v3 "f_cw=0.50이 QC1 해결" **basis-A서 붕괴 확정**(잔차 2.3→16.8/48.5/54.8). 막 down-weight α는 **transport ceiling이 금지**(장쇄서 α_F2<α_min) → wheat subcellular anchor는 장쇄 rice에 전이 불가. **장쇄 α≈1(no down-weight) = transport-consistent**. α 점추정 불가(W2 포화로 kappa_d와 교란), ceiling 하한만.
- **★ surface — Li2025 "surface excess"는 water-분모 artifact, surface 증거 아님**. excess가 사슬길이 아닌 **pore-water 품질**을 따라감(good→0, poor→큼). clean 단일점 PFOA(good)는 sub-equilibrium = Yamazaki와 동일. **H7 §4 "field-dependent surface" 정정 → K_surf=0 default 정당**, K_surf-vs-토양 회귀 근거 없음.
- **★ W3 redox — soil_paddy 부호 뒤집힘 정정**. 기존 default는 "침수→흡착 약화→bioavailability↑"(틀림). 정정: flooded → **dilution(θ_g↑)+leaching(용존상 1차손실)로 Cwo↓**. redox→sorption은 부호 불확실 2차효과로 중립 default 강등.
- **★ Gap4 — foundation 검증 + TF cross-field**. basis-A B_k(전기관)+W2 transport이 full-ODE로 Yamazaki 재현(**log10 RMSE 0.029**, PFDoDA만 −25% near-MQL outlier). **TF=tissue/tissue는 water-독립** → BAF로 불가능했던 cross-field 검증 성립: TF_straw 단조감소(TSCF)가 Yamazaki·Li2025 **양 field에서 재현**.

---

## 1. α/QC1 — basis-A 재작업 (Prompt3 무효화 정정 + 신규 식별성)

**구성**: root θ=0.90, f_prot=0.07, f_PL=0.015, f_cw=0.50(poly0.40+lignin0.10), K_cw=rice whole-cw root, K_PL/K_prot=Bk_table_S5.

1. **막 지배 basis-A 지속**: root 막분율 PFHpA 46%·PFOA 68%·PFNA 86%·PFOS 96%·PFDoDA 98%. (1−θ)는 결합 pool 전체를 같이 축소 → 정성 불변. naive 75–98%와 동일.
2. **Prompt3 v3 폐기 확정**: 동일 구성(f_cw0.50+rice K_cw+막→0.30)의 총잔차(|F1|+|F2|+|F3|) = naive **2.3 → basis-A 16.8(θ0.70)/48.5(θ0.90)/54.8(θ0.92)**. (1−θ)가 cw pool 축소 → cw share가 Liu F1 밴드(47.7–59.0%)에 **구조적 도달 불가**(θ0.90서 cw max 30%). H7 §1 재현·확정.
3. **★ NEW — α 식별성**: 막 down-weight α를 두 근거로 경계:
   - **ceiling 하한 α_min** (Yamazaki B_root ≥ obs): C4–C9는 **0**(막 빼도 ceiling>obs), 장쇄 양수 — PFDA .114, PFUnDA .376, **PFDoDA .680**, PFOS .083, PFBS .831.
   - **Liu-F2 anchor α_F2** (막→0.30): PFDA .030, PFUnDA .015, PFDoDA .008, PFOS .018.
   - **C10+·PFOS에서 α_F2 < α_min** → Liu가 요구하는 down-weight 적용 시 B_root가 관측 root BAF 아래로(불가; Yamazaki surface≈0). **장쇄 막 지배는 실재**: 높은 장쇄 root 축적(PFDoDA BAF 69=ceiling 0.69×)이 큰 B_root=큰 막분배를 요구. Prompt3의 "막→wheat F2 down-weight" 프로그램은 **장쇄 rice transport와 모순**.
4. **verdict**: α는 BAF로 점추정 불가(W2 3param/3obs 포화 → α는 kappa_d로 흡수, 예측 BAF 불변). 단·중쇄(C4–C9): 양방향 무제약+막분율 낮음+deep sub-eq → 예측·해석 둘 다 둔감(저우선). 장쇄(C10+): ceiling이 **α≈1 근처로 하한** → **no down-weight 채택**. 결정 데이터 = **rice(wheat 아님) per-congener root subcellular** 또는 직접 in-planta K_PL.
5. **부차**: root `f_PL_membrane`은 estimate(0.01–0.02, **2× 범위**). 장쇄 B_root ∝ f_PL·K_PL → 막분율·B_root 절대값 지배. 직접 측정 우선순위 높음.

산물: `S6_alphaQC1_basisA.py`, `S6_alphaQC1_root_basisA.csv`, `S6_alphaQC1_root.png`.

---

## 2. surface cross-field — Li2025 vs Yamazaki (H7 §4 정정)

`surface_excess = max(0, obs_root_BAF − B_root_ceiling)`.

- **Yamazaki(clean per-congener water)**: 전 11종 obs/ceiling 0.06–0.95, **전부 sub-equilibrium, excess 0**.
- **Li2025(Tianjin, group-water)**: 표면상 5–22× 초과지만 **excess가 사슬 아닌 pore-water 품질을 따라감**:

  | water 품질 | congener | obs/ceiling |
  |---|---|--:|
  | **good** | PFOA(C8) | **0.53** (sub-eq) |
  | med | PFBA(C4) | 6.9 |
  | rough | PFBS(C4)/PFHxA(C6) | 16.4 / 22.1 |
  | poor | PFOS(C8) | 5.1 |

- **진단**: 신뢰 가능한 유일점 PFOA(good)만 sub-equilibrium → Yamazaki와 동일. excess는 강흡착(낮은 aqueous %, PFOS 0.3–1.6%)일수록 water 분모가 과소·불확실 → BAF 인위 팽창. 사슬로 설명 안 됨(C4 PFBA 6.9 > C8 PFOA 0.53).
- **결론**: H7 §4 "surface field 의존(Tianjin 큼)" → **water-분모 artifact로 재해석**. **field-dependent surface 증거 없음. K_surf=0 default 정당**, K_surf-vs-토양(foc/pH) 회귀는 신호 자체가 없어 **근거 없음**. Fe/Mn plaque 흡착 실재(rice_tissue NOTES #5)하므로 K_surf capability는 유지 — Li2025로 calibrate 불가일 뿐. 진짜 검정엔 **신뢰 가능 per-congener pore-water 또는 hydroponic RCF** 필요.

산물: `S6_surface_crossfield.py`, `S6_surface_crossfield.csv`, `S6_surface_crossfield.png`.

---

## 3. W3 — paddy redox 부호 정정 (H7 §7.4)

- **문제**: 원본 `soil_paddy.example_paddy_redox`가 flooded K_F(1.0)<drained(2.0) = "침수→흡착 약화→Cwo↑"(부호 오류; 검증 flooded/drained=**1.76 상승**).
- **정정** `soil_paddy_redox_corrected.py`:
  - **dilution**: flooded θ_g↑(0.35→0.60) → 동일 inventory서 Cw↓ (검증 0.88).
  - **leaching**: 침수기 용존상 1차손실 `k_leach·f_diss·C_T` (이동성 단쇄 빠름) → Cwo 궤적 단조 하락(2.08→1.28).
  - **redox→sorption 중립 default**: 혐기 K_F 변화 부호 불확실(Fe/Mn 환원용해 vs 환원 OM 흡착↑) → "약화" 폐기, 데이터로만 설정하는 hook만.
  - API/import 원본 호환(`PaddyRedoxCorrected`, `FreundlichSoil` 재사용) = drop-in.
- **한계**: θ_g·k_leach 값은 illustrative — 실측 flooding schedule + leaching rate로 calibrate 필요.

산물: `soil_paddy_redox_corrected.py`, `S6_W3_redox_correction.png`.

---

## 4. Gap4 — B_k/TF 독립 검증

1. **Full-ODE 재현**: first-principles basis-A B_k(전기관) + W2 transport(f_xy,L_Ph,kappa_d) → 4-comp ODE → Yamazaki. **log10 RMSE 0.029**(11종 중 10 거의 정확). **PFDoDA −25%**(pred root 54 vs obs 69) = H7 명시 near-MQL water(0.02 ng/L) outlier. → foundation(`S6_Bk_basisA_allorgan.csv`)+W2 파라미터가 관측 재생성 확인.
2. **★ Cross-field TF (water-독립)**: TF=tissue/tissue라 Li2025 water artifact 면역. TF_straw 단조감소(TSCF):
   - Yamazaki C4→C11: 16.2→10.1→4.1→1.4→…→0.42
   - Li2025 C4/C6/C8: 2.73→1.23→0.69
   절대값 차(Yamazaki PFBA 16 vs Li2025 2.7) = Yamazaki stem 높이구배(단쇄 증산농축, 4높이 평균 smear). **경향은 cross-field 일치** → TSCF 두 독립 field 검증. grain phloem-limited(clean TF 0.8–1.35); Li2025 grain TF 팽창(PFBS 19)은 dry-grain(θ0.14) vs wet-root(θ0.90) fw 환산 ~8.6× 일부.

산물: `S6_Gap4.py`, `S6_Gap4_ode_repro.csv`, `S6_Gap4_TF_crossfield.csv`, `S6_Gap4.png`.

---

## 5. Foundation — basis-A B_k(n) all-organ (W1 재생성)

`S6_Bk_basisA_allorgan.csv`: basis-A B_k(n) root/stem/leaf/grain 전 12종(α=1, rice_tissue rec 조성, 기관별 whole-cw K_cw). **이것이 정본 B_k이며 naive `Bk_table_S5.csv`(B_root 22.97 등)를 대체**. 예: PFOA B_root 4.17(θ0.90), PFOS 49.45, PFDoDA 100.99.

---

## 6. config 표준화 결정 필요 (다음 세션 착수 전)

- **root θ**: 0.70(W1 초기) vs **0.90**(rice_tissue rec·Liu 측정 0.90–0.92). B_root 절대값 변동. **권장: 0.90 통일.**
- **root f_PL**: 0.01 vs **0.015**(rec). 장쇄 B_root 지배. 직접측정 전까지 rec 0.015.
- (위 결론 1–4는 두 선택 모두에 robust — θ0.90서 conflict 성립, θ0.70서 막 지배 더 강함.)
- **grain θ**: 수확 0.14(Yamazaki dry grain 정합) vs 등숙 0.30. 단계 의존 — 용도별 선택.

---

## 7. 이월 / 미해결 (전부 data-limited; 모델링 아닌 실험 의존)

1. **rice per-congener root subcellular fractionation**(wheat 아님) → α 점추정 + QC1 종결. [§1]
2. **신뢰 가능 per-congener pore-water 또는 hydroponic RCF**(음이온, per-congener) → surface 진짜 검정 + f_xy↔binding 교란 분리. [§2, H7 §7.2]
3. **실측 Qtp(t)·M(t)**(Yamazaki 증산·기관질량) → f_xy 절대값(현 상대추세만 유효). [H7 §6]
4. **벼 cw 단당 조성**(Guo Table11 벼판) → K_cw_poly grass 재-anchor. **장기 최약점**(직접 K_cw_poly 측정 부재). [Prompt2 backlog #1]
5. **Bk_table_S5 갱신**: §5 정본으로 교체(housekeeping, 신규 데이터 불요).
6. **gap**: PFSA C9–C14 / PFCA C13–14(Kcw_v2는 PFOS/PFTeDA까지) — 필요 시 slope0.10+offset0.40 외삽(명시).
7. **soil_paddy_corrected**: θ_g·k_leach calibration + 실측 flooding schedule.

### 모델링으로 지금 가능한 잔여(데이터 충분)
- (a) Bk_table_S5 → 정본 교체 + 전 코드 참조 갱신.
- (b) **통합 soil+plant run**: `soil_paddy_redox_corrected` → plant ODE 결합, 현실 flooding schedule로 시나리오. (부품 다 있음, 조립만.)
- (c) **불확실성 전파**: f_PL 2× → B_root → 예측 sensitivity.

---

## 8. 파일 manifest

| 파일 | 내용 |
|---|---|
| `H8_handoff_S6_final.md` | 본 핸드오프 |
| **code/** | |
| `S6_alphaQC1_basisA.py` | α/QC1 basis-A 재작업(분해·v3 잔차·ceiling 식별성) |
| `S6_surface_crossfield.py` | surface excess cross-field + water-품질 진단 |
| `soil_paddy_redox_corrected.py` | W3 정정 모듈(dilution+leaching, drop-in) |
| `S6_Gap4.py` | full-ODE 재현 + water-독립 TF cross-field |
| **data/** | |
| `S6_Bk_basisA_allorgan.csv` | **정본 basis-A B_k(n) 전기관** (S5 대체) |
| `S6_alphaQC1_root_basisA.csv` | root 분해·α_min·α_F2·feasibility |
| `S6_surface_crossfield.csv` | obs/ceiling·excess·water 품질, 두 field |
| `S6_Gap4_ode_repro.csv` | ODE pred vs obs(root/straw/grain) |
| `S6_Gap4_TF_crossfield.csv` | Yamazaki vs Li2025 TF |
| **figures/** | |
| `S6_alphaQC1_root.png` | 막분율 & α 경계(ceiling vs anchor) |
| `S6_surface_crossfield.png` | excess vs 사슬 / vs water 품질 |
| `S6_W3_redox_correction.png` | Cwo(t) OLD vs CORRECTED |
| `S6_Gap4.png` | ODE 재현 + TF cross-field |

**입력 원본**(uploads, 본 세션): `pfas_rice_plant_module_4pool(_surf/_5pool).py`, `calibration.py`, `soil_paddy.py`, `Kcw_Klignin_params_v2.csv`+Prompt2 handoff, `rice_tissue` 패키지, `Bk_table_S5.csv`(Prompt3), `obs_baf_{Yamazaki,Li2025}.csv`, `Li2025_BAF_TF.csv`, `W2_transport_fit.csv`, H7.

---

## 9. H7 대비 갱신 요약 (충돌 해소 기록)

| H7 진술 | H8 갱신 |
|---|---|
| §4 surface는 field 의존(Tianjin 큼, Andosol ~0) | Li2025 excess=water-분모 artifact. **field-surface 증거 없음**, K_surf=0 default. |
| §7.1 막 down-weight α는 OPEN | 장쇄 α 점추정 불가지만 **ceiling이 α≈1로 하한 → 막 지배 실재**, no down-weight. Prompt3 program 장쇄서 기각. |
| §7.4 redox 부호 정정 필요 | **완료**(dilution+leaching). |
| §7.5 Gap4 미착수 | **완료**(RMSE 0.029 + TF cross-field). |
| Bk_table_S5 = B_k 입력 | **`S6_Bk_basisA_allorgan.csv`로 대체**(S5는 naive). |

*끝. S6 본작업(W1–W3 + Gap4 + α/QC1 + surface) 종결. 잔여는 전부 실험 데이터 의존 또는 housekeeping.*
