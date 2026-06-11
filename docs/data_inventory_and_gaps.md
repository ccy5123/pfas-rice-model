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
soil sorption.

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
| **Tang 2026** (`raw_si/tang2026_tf_bcf`) | paddy PFOS/PFOA TF, BCF | f_xy head-group offset sign | — |

**Meta:** quantitative calibration effectively rests on **Yamazaki alone (33 pts)**. Li has
caveats; Kim is grain-only. Hence the demo's log10 RMSE 0.029 is a *saturated reproduction*,
not a predictive test.

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

### ④ Soil half
- **Anoxic/flooded (negative Eh) sorption isotherms** on paddy soil (track Fe(II) release)
  → all current isotherms are aerobic; flooded redox governs C_w^o (DB gap, HIGH).
- **Real field flooding schedule + site soil/loading**, and reliable **per-congener** pore
  water (Li is group-water).

### ⑤ Scope & reliability
- **GenX / F-53B / 6:2 FTS** rice data (very sparse; authentic standards already procured).
- **Very-long-chain (PFDoDA+)** analytical reliability — extend data above MQL.

---

## One-line summary
The model currently **reproduces** rice PFAS BAFs by fitting transport parameters to one
greenhouse study on top of measured binding coefficients and crop forcing. To **validate /
predict**, the missing ingredient is data, not model: a tissue- and time-resolved greenhouse
time-course, xylem/phloem sap (translocation & grain loading), root subcellular fractionation
(the sequestration sink), and anoxic paddy-soil sorption.
