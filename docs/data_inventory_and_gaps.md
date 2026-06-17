# Data inventory & gaps — PFAS bioaccumulation in rice

> What data the model currently rests on, and what data it still needs, organised by the
> model structure. Companion to `DELIVERABLE_GAP_A_Kcw.md`, `DELIVERABLE_GAP_B_fxy.md`,
> `literature_db/Gap_Analysis.csv`, and the two exploration notes
> (`nstem_gradient_exploration.md`, `fxy_longchain_lipid_exploration.md`).

## TL;DR
The model sits on **measured partition coefficients (K_PL, K_prot) + measured crop forcing
(Q_TP, M)**, with the transport parameters (f_xy, L_Ph, root sink) **fitted to a single
greenhouse BAF study (Yamazaki 2023, 33 points)**. Reproduction is therefore *saturated*
(≈3 transport params per 3 observed BAFs per congener), **not validation**. Moving to
prediction is **data-limited, not modelling-limited**: the priorities are (a) a tissue- and
time-resolved greenhouse time-course, (b) xylem/phloem sap (direct translocation & grain
loading), (c) root subcellular fractionation (the sequestration sink), (d) anoxic/flooded
soil sorption, (e) **filling-grain subcellular PFAS + phloem sap** (a *new* structural gap —
on a dry-weight basis the model grain TF is ~3–8× low and unclosable by `L_Ph`/lipid; §2④).

> **Two open issues surfaced this session (both now in §1a/§2④ + Gap_Analysis):**
> (1) **`f_xy` is condition-dependent** — PFOS `f_xy` ≈ 0.14 (Yamazaki clean water) vs ~0.32
> (Tang flooded soil); report the range with conditions, don't pin one value.
> (2) **Grain accumulation is a structural sink gap** — only revealed once the model→Tang
> TF comparison was put on a consistent **dry-weight** basis (`docs/tang2026_grain_units_exploration.md`).

This session's deep-dive resolved the long-chain behaviour into **three distinct mechanisms**,
each needing a different measurement:
1. **free-anion uptake/TSCF** (short chains) — GHK exclusion + carrier;
2. **lipid-facilitated shoot translocation** (long chains) — K_PL-driven; explains the
   U-shaped effective f_xy; wired as the opt-in `g_xy`/`g_ph` term;
3. **root sequestration sink** (very long chains) — chain-length/head-group-specific, **not**
   K_PL-tracking (PFOS & PFUnDA share K_PL=31623 but roots 5.93 vs 19.53), so a two-pool root
   could not close it without per-congener overfitting (negative result).

---

## 1. Data used so far

### 1a. Observational data (calibration / validation targets)

| dataset | scope | role | confidence / caveat |
|---|---|---|---|
| **Yamazaki 2023** (`10.1021/acs.est.2c08767`) | 11 congeners × root/straw/grain = **33 pts**; Andosol, clean per-congener water | **primary calibration** | single greenhouse study; PFDoDA near-MQL outlier |
| Yamazaki S18/S19 | stem by height (0–20…>60 cm) | nstem gradient direction | qualitative |
| **Li 2025** (Tianjin field) | 5 congeners (PFBA, PFHxA, PFOA, PFBS, PFOS) × root/straw/grain | cross-field TF (water-independent direction) | group-water exposure; surface confound |
| **Kim 2019** (`10.1016/j.scitotenv.2019.03.240`) | per-congener brown-rice BAF + paired pore water | L_Ph fit (PFOA, 1 pt) | grain-only |
| **Tang 2026** (`raw_si/tang2026_doseresponse`) | PFOA/PFOS/**GenX** per-organ TF + BCF, **5 soil doses** 0.1–100 µg/g; Nipponbare, flooded paddy, 150 d | f_xy head-group offset sign (0.1 µg/g); dw OOS TF | TF declines with dose (toxicity); soil-spiked, **≠** Yamazaki clean water |

**Meta:** quantitative calibration effectively rests on **Yamazaki alone (33 pts)**. Li has
caveats; Kim is grain-only. Hence the demo's log10 RMSE 0.029 is a *saturated reproduction*,
not a predictive test.

**Condition-dependence of `f_xy` (this session).** Re-fitting `f_xy` to translocation TF
gives **different values across studies because the conditions differ**: PFOS `f_xy` ≈
**0.142** (Yamazaki — Andosol, clean per-congener water, greenhouse) vs **~0.32** (Tang —
flooded paddy-soil pot, Nipponbare, 5 doses). The model build's `f_xy` is therefore
**conditional**: do NOT pin PFOS `f_xy` to one number; carry the **0.14–0.32 range with its
conditions**. GenX, by contrast, agrees (~0.013–0.02) and independently confirms the
provisional 0.233 is ~12–18× too high. (The grain/endosperm miss is a *separate, structural*
gap — see §2④, dry-weight basis.) Closing this needs **paired multi-condition translocation
data** (xylem-sap `f_xy` across soil/hydroponic + cultivar), not a re-fit.

### 1b. Parameter values (measured / QSPR / fitted)

| component | source | status |
|---|---|---|
| **K_PL** membrane/phospholipid–water | Chen 2025 SSLM K_MW, per-congener (vs Droge 2019) | **measured** |
| **K_prot** protein–water | Zhou 2025 dialysis K_prow (soy = plant, BSA = animal) | **measured** |
| **K_cw** cell-wall | GAP A: Guo 2025 DFT ladder + Mel 2024 lignin anchor | **anchored** (not measured) |
| **Koc** soil | Higgins & Luthy slope (+0.55/CF₂, +0.23 sulfonate) + Milinovic 2015 anchor | QSPR |
| **f_d** dissociation | pKa (Goss 2008) → ≈1 | settled |
| **E_m** membrane potential | rice root ~−120 mV (Wang 1994) | literature; in-situ paddy unmeasured |
| carrier **Vmax/Km** | fixed nominal | **unmeasured** (assumed) |
| **f_xy** root→shoot | GAP B theory + W2 fit → revised this session to a **U-shape (lipid-facilitated)** | **fitted** (not measured) |
| **L_Ph** phloem→grain | Kim PFOA, 1-pt fit | **fitted** |
| **g_xy / g_ph** lipid-bound loading | this session, K_PL-gated fit to Yamazaki (excl. PFDoDA) | **fitted**, opt-in, exploratory |

### 1c. Tier-0 forcing + soil engine

| input | source | status |
|---|---|---|
| **Q_TP(t)** transpiration | `forcing_rice`: FAO-56 dual-Kc, Kumari 2022 + Nay Htoon 2018 | measured crop physiology |
| **M_k(t)** organ biomass | `growth_rice`: ORYZA IR72 (partitioning, HI) | crop model |
| **C_w^o(t)** pore water | `soil_hydrus`: compiled HYDRUS-1D, linear Kd (Koc) | real engine, synthetic scenario |

Curated DB: `literature_db/` C1 (rice BAF) · C2 (forcing) · C3 (QSPR/sorption) · C4 (binding)
· C5 (membrane) · C6 (physchem) + Gap_Analysis + `raw_si/` (7 SI extractions).

---

## 2. Data needed — organised by what it resolves

### ⓪ Top priority — validation data (breaks the single-study, saturated-fit problem)
- **Greenhouse time-course**: environmental + elevated doses, 2 cultivars (Indica/Japonica),
  **6–8 destructive harvests**, separating root / straw / husk / brown / polished, with
  **paired pore-water, M(t), Q_TP(t) at every timepoint**. → observations outnumber
  parameters → genuine *out-of-sample prediction* and hold-out. (Gap_Analysis HIGH #1)

### ① Free-anion uptake / membrane crossing (Tier-2; resolves the identifiability gap)
- **Inhibitor ± concentration-series root-uptake assays** → separate carrier vs channel
  (not separable from BAF), pin Vmax/Km.
- **In-situ paddy E_m** (microelectrode) → test the GHK exclusion term directly.

### ② Lipid-facilitated shoot translocation (this session's mechanism — direct test)
- **Xylem-sap f_xy**: root-pressure exudate [PFAS]/[root-water] per congener → directly
  **falsifies/confirms the U-shape and its K_PL gating**.
- **Phloem-sap L_Ph**: peduncle phloem [PFAS] → grain loading (currently a 1-pt fit).

### ③ Root sequestration sink (the new gap exposed by the negative two-pool result)
- **Root subcellular fractionation** (apoplast / symplast / cell-wall / Fe-Mn plaque) **per
  congener, resolved by chain length AND head group** → explain why the long-chain root BAF
  does not track K_PL (PFOS vs PFUnDA).
- **Direct K_cw**: batch sorption on isolated rice-root cell-wall fractions (GAP A; no
  literature coefficient exists).

### ④ Grain accumulation sink (NEW this session — exposed by the dw-consistent Tang TF)
- On a **dry-weight** basis (the only consistent way to compare to Tang's TF; see
  `docs/tang2026_grain_units_exploration.md`), the model's grain/endosperm TF is
  **structurally under-predicted ~3–8×** (PFOA endosperm 0.11 vs Tang 0.95), and **no
  `L_Ph` value — even maxed — nor the lipid `g_ph` term closes it** (ceiling ~0.3–0.7 vs
  Tang 0.9–1.6). The nstem_leaf shoot-distribution fix *starves* the grain relative to the
  4-pool leaf-sink, so shoot split and grain accuracy **trade off** — it is not an `L_Ph`
  calibration issue.
- **Filling-grain subcellular PFAS** (endosperm vs aleurone vs pericarp) **+ peduncle
  phloem-sap [PFAS]** per congener → decide between (i) a higher effective grain binding
  `B_grain` (storage-protein/starch sequestration beyond the current Chen K_PL / Zhou
  K_prot) and (ii) active phloem unloading/enrichment. Until then, grain TF should be
  reported as an **order-of-magnitude** quantity carrying this known low bias.

### ⑤ Soil half
- **Anoxic/flooded (negative Eh) sorption isotherms** on paddy soil (track Fe(II) release)
  → all current isotherms are aerobic; flooded redox governs C_w^o (DB gap, HIGH).
- **Real field flooding schedule + site soil/loading**, and reliable **per-congener** pore
  water (Li is group-water).
- **Ether/sulfonamide soil Koc** — no measured PFECA/sulfonamide paddy Koc exists, so the
  structure adapter's `KOC_ETHER_LOG_OFFSET = 0` (carboxylate approximation; Gap_Analysis C3).
  Needs batch Koc on ether-PFCAs (GenX/ADONA) + FOSA-type sulfonamides.

### ⑥ Scope & reliability
- **GenX / F-53B / 6:2 FTS** rice data (very sparse; authentic standards already procured).
- **Very-long-chain (PFDoDA+)** analytical reliability — extend data above MQL.
- **Structure-adapter binding QSPR** (`src/pfas_structure.py`) — the per-ether-O **K_PL** term
  (`KPL_ETHER_LOG_OFFSET = −0.49 log`) rests on a **single anchor (GenX)**, and the
  **sulfonamide K_PL slope has no data** (uses the carboxylate slope). Needs more measured
  ether/sulfonamide K_MW to multi-anchor (Gap_Analysis C4).

---

## One-line summary
The model currently **reproduces** rice PFAS BAFs by fitting transport parameters to one
greenhouse study on top of measured binding coefficients and crop forcing. To **validate /
predict**, the missing ingredient is data, not model: a tissue- and time-resolved greenhouse
time-course, xylem/phloem sap (translocation & grain loading), root subcellular fractionation
(the sequestration sink), and anoxic paddy-soil sorption.
