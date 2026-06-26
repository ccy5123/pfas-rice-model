# PFAS–Rice Uptake Explorer — User Manual

**English** | [한국어](MANUAL_KR.md) &nbsp;·&nbsp; [← README](../README.md)

> The complete user guide for the interactive tool that estimates how much **PFAS
> ("forever chemicals")** a rice plant takes up from contaminated paddy water/soil, and
> **where it ends up** (roots, straw, grain). Written for **both general users**
> (policy, students, the public) **and experts** (environmental science, research).
>
> - General users only need **[1. Quick start](#1-quick-start)** → **[4. Simple mode](#4-simple-mode)**.
> - Experts should also read **[5. Expert mode](#5-expert-mode)**, **[7. Input data formats](#7-input-data-formats-csv)**, and **[9. Scientific background](#9-scientific-background--assumptions)**.

> [!WARNING]
> **Research & educational model — illustrative estimates only.** This is **not** a
> regulatory, food-safety, or health determination. Do **not** use it for real exposure
> or safety decisions. Results are illustrations under the model's assumptions.

> **Language note.** The app is bilingual: the default **general-audience (Simple) view is
> Korean**, and the **Expert view is English**. This manual describes both; the Korean UI
> labels are given in parentheses where helpful, and the Expert section quotes the English
> UI verbatim.

---

## Contents

1. [Quick start](#1-quick-start)
2. [What the tool does](#2-what-the-tool-does)
3. [The two modes and switching](#3-the-two-modes-and-switching)
4. [Simple mode](#4-simple-mode)
5. [Expert mode](#5-expert-mode)
6. [Shared features](#6-shared-features)
7. [Input data formats (CSV)](#7-input-data-formats-csv)
8. [Reading the results](#8-reading-the-results)
9. [Scientific background & assumptions](#9-scientific-background--assumptions)
10. [Glossary](#10-glossary)
11. [Troubleshooting (FAQ)](#11-troubleshooting-faq)
12. [Reproduce / validate / test](#12-reproduce--validate--test)
13. [Cite / license / references](#13-cite--license--references)

---

## 1. Quick start

### A. The deployed web app (no install)
Open the deployment URL in a browser (Streamlit Community Cloud, tracking the `main` branch).
The first screen is the **Simple mode (Korean)**; flip the sidebar toggle for the English Expert UI.

### B. Run it on your computer
With Python 3.9+:

```bash
git clone https://github.com/ccy5123/pfas-rice-model.git
cd pfas-rice-model
pip install -r requirements.txt     # full app stack (numpy/scipy/streamlit/plotly/pandas/rdkit/phydrus)
streamlit run app.py                 # opens http://localhost:8501
```

- 4 of the 5 exposure modes (Model / HYDRUS-CSV / Soil-inventory / Biomonitoring) and SMILES input work immediately.
- Only **Run HYDRUS-1D (live)** needs the compiled engine (gfortran) + `phydrus`; it is auto-hidden/guided when absent.
- Static PNG figure export needs the optional `kaleido` (+ Chrome); without it the CSV downloads still work.

### First 30 seconds (Simple mode)
1. In the sidebar **① Pick a chemical (화학물질 선택)**, choose a PFAS (e.g. PFOA).
2. In **② How contaminated? (오염 정도)**, pick low / medium / high.
3. In the **🗺️ Where it goes (어디로 가나)** map, hotter colour = more PFAS in that part.
4. Read the headline cards (PFAS in roots / straw / grain, µg/kg) and the one-line summary.

---

## 2. What the tool does

It visualises a **4-compartment dynamic uptake model** of rice (root → stem → leaf → grain).

- It draws the paddy soil and the rice plant to scale and fills each compartment with an **accumulation heat colormap**.
- A day slider / ▶ play shows the build-up across one season (transplant → harvest).
- Given **one PFAS + a contamination level**, it estimates root/straw/grain concentrations and the **build-up factor (BAF)**.
- It can also work **backwards**: from measured tissue concentrations it **infers the soil-water contamination level** (with an uncertainty range).

Targets are the **13 curated congeners** (PFCA C4–C12: PFBA·PFPeA·PFHxA·PFHpA·PFOA·PFNA·PFDA·PFUnDA·PFDoDA;
PFSA C4/C6/C8: PFBS·PFHxS·PFOS; ether-PFAS: GenX); in Expert mode **any PFAS can be entered by SMILES structure**.

> Compute is decoupled from the UI: `src/model_api.py` (math), `src/plots.py` (figures), `app.py` (UI).
> Model equations/parameters are detailed under `docs/` (esp. `docs/OVERVIEW_KR.md`, `docs/visualization_tool.md`).

---

## 3. The two modes and switching

| | Simple mode (default) | Expert mode |
|---|---|---|
| Language | **Korean** | **English** |
| Audience | policy / students / public | environmental science / research |
| Input | chemical + level (low/medium/high) | 5+1 exposure modes, SMILES, all parameters |
| Tabs | 5 (plain language) | 9 (technical) |
| Symbols | none (BAF/Cwᵒ/f_xy/eᴺ hidden) | all exposed |

**How to switch:** the toggle at the top of the sidebar — **🔬 전문가/고급 모드 (Expert / advanced)**.
- **Off (default)** = Simple mode (Korean)
- **On** = Expert mode (English) — every Simple feature **plus** the full research interface (nothing is removed).

---

## 4. Simple mode

### 4.1 Sidebar
- **① Pick a chemical**: one of 13, with a plain description (e.g. "PFOA — common 'forever chemical' acid (C8)"). Longer chains generally stick more.
- **② How contaminated?**: a preset for the PFAS dissolved in the soil water — **low = 0.1 µg/L, medium = 1 µg/L, high = 10 µg/L** (higher → more enters the plant).
- **📋 Use my own data tables** (optional): enter your own **growth curve** and **time-varying soil-water contamination** as tables → see [6.2](#62-your-own-data-tables-growth--pore-water).
- A hint points to the Expert toggle for the full research interface.

### 4.2 Headline (summary cards)
- **In the roots / straw (stems+leaves) / grain (edible rice)** — the PFAS concentration (µg/kg) in each part.
- A one-line summary below: at your chosen level, ~X µg/kg in the grain, ~N× the soil water, most stays in ○○.

### 4.3 Tabs
- **🗺️ Where it goes** — the plant + soil picture; hotter colour = more PFAS. Use the **day** slider, or **▶ play the season** for an animation.
- **📈 Build-up over time** — per-part concentration over the season. The **grain** only takes up PFAS once it forms (around flowering).
- **📊 How much builds up** — final per-part concentrations at harvest. Expand **🔬 Compare with real-world measurements (Yamazaki 2023)** to check the model vs data (closer bars = better fit).
- **🔎 Work backwards** — infer the soil-water level from a lab result → see [6.1](#61-work-backwards-bayesian-inverse).
- **ℹ️ About & glossary** — what the tool does, how to read it, a plain-language glossary, and the disclaimer.

### 4.4 Downloads
Under the headline, **⬇️ Download these results** gives a summary table (CSV), the full time series (CSV), and the plant map (PNG, when kaleido is present).

---

## 5. Expert mode

Turning the toggle on switches to the English UI; the sidebar becomes **1 · Data source / 2 · PFAS compound / 3 · Scenario**.

### 5.1 Data source (1 · Data source) — how is Cwᵒ(t) supplied
| Mode | Cwᵒ(t) from | Q_TP·M(t) | When to use |
|---|---|---|---|
| **Model (parametric)** | a constant you set (or a flooded shape) | measured FAO-56 / ORYZA | quick what-ifs, teaching |
| **Custom tables (Cwᵒ + growth)** | tables you enter | the table's growth + measured transpiration | you have your own growth/exposure series → [6.2](#62-your-own-data-tables-growth--pore-water) |
| **HYDRUS / CSV drivers** | a HYDRUS-1D/Phydrus run (CSV) | the CSV, or measured | you have a calibrated soil-water-solute model |
| **Run HYDRUS-1D (live)** | the real HYDRUS-1D engine, executed here | HYDRUS root uptake + ORYZA | you want the engine to run in-app (needs it built) |
| **Soil inventory → pore water** | inverting a total soil load (Freundlich) | measured | you know total soil PFAS, not pore water |
| **Biomonitoring (measured tissue)** | a measured soil-water value | not needed | you have field tissue + water concentrations |

### 5.2 Compound (2 · PFAS compound)
- **Curated congener**: one of 13 (measured/cited calibrated parameters).
- **SMILES (structure)**: paste any PFAS structure → RDKit extracts descriptors → (1) **read-across** if it matches a curated congener, or (2) a **QSPR** (provisional) for a novel structure. A 2-D structure preview is shown. (`docs/structure_input.md`)

### 5.3 Model parameters (sidebar)
- **E_m [mV]** (root membrane potential): the GHK anion-exclusion lever (rice −116…−140 mV). More negative → more anion exclusion.
- **f_xy source**: `recommended` (monotone, physical TSCF) / `W2 fit` (reproduces Yamazaki).
- **Biomass driver M(t)**: `ORYZA2000` (mechanistic carbon balance, default) / `growth_rice` (IR72 partitioning × logistic).

### 5.4 Scenario (3 · Scenario) — per-mode controls
- **Model (parametric)**: `Pore-water Cwᵒ [µg/L]`, `Season length [days]`, `Cwᵒ(t) shape` (constant / flooded(dilution+leaching), per-congener HYDRUS-calibrated `k_leach`), `Measured forcings` toggle.
- **Custom tables**: the tables are entered in the main panel ([6.2](#62-your-own-data-tables-growth--pore-water)).
- **HYDRUS/CSV**: upload a driver CSV or use the bundled example. Columns in [§7](#7-input-data-formats-csv).
- **Run HYDRUS-1D (live)**: `f_oc`, `Flooded until [day]`, `Percolation [cm/day]`. A "Build the engine" button appears if unbuilt.
- **Soil inventory**: `Total soil inventory [µg/kg dry]`, Freundlich `K_F`/`n`/`θ_g`, flooded flag, `k_leach`.
- **Biomonitoring**: manual input (root/straw/grain conc + Cwᵒ) or CSV.

### 5.5 Expert tabs (9)
1. **🗺️ Plant & soil map** — accumulation map (concentration/BAF toggle, day slider/animate).
2. **📈 Tissue dynamics** — tissue concentration C_k(t) + **PFAS mass (burden) C_k·M_k**. Shows B_k/f_xy/L_Ph/κ_d. Explains the grain formation gate.
3. **🟫 Soil & drivers** — the actual Cwᵒ(t)·Q_TP(t)·M(t) drivers, the (soil-inventory) isotherm, and the soil profile.
4. **📊 BAF vs observed** — model vs Yamazaki 2023 bars. Optional **two-pool (seq)** exploratory overlay. Yamazaki conditions/matching explained.
5. **🔗 Chain-length trends** — chain length vs a parameter (K_PL/K_prot/K_cw/f_xy/B_root/B_grain).
6. **⚖️ Compare congeners** — per-tissue BAF across the congeners selected in the sidebar.
7. **✅ Tang TF (OOS)** — Tang 2026 per-organ TF (out-of-sample) check (PFOA/PFOS/GenX, dry-weight, with the f_xy refit).
8. **🔎 Inverse (Bayesian)** — tissue concentrations → exposure Cwᵒ, with the identifiability caveat → [6.1](#61-work-backwards-bayesian-inverse).
9. **ℹ️ About / coupling** — the modes, the HYDRUS-1D input/output mapping, and the glossary.

---

## 6. Shared features

> The two features below are available in **both** modes (only the labels differ KR/EN).

### 6.1 Work backwards (Bayesian inverse)
Answers "I measured PFAS in my rice — how contaminated was the soil water?" **with an uncertainty range**.

- **Input**: measured root/straw/grain concentrations (µg/kg), any subset, + a measurement-precision setting (Typical ±~40% / High precision ±~20% / Rough ±~2×).
- **Button**: 📐 Estimate the contamination level (the estimate solves the ODE a few times, so it is gated behind a button).
- **Output**: the most likely soil-water level (µg/L) + a 95% credible interval + a posterior curve + a check that the model reproduces your inputs.
- **Why it's a real inverse**: root uptake is **saturable** (GHK + carrier), so tissue conc is a **nonlinear** increasing function of Cwᵒ — not a division. It finds the MAP via a **quadratic-fit Laplace** in log10(Cwᵒ) (MAP + curvature = posterior width; `model_api.estimate_exposure_bayesian`).
- **Limit (identifiability)**: only the **exposure level** is estimated. From tissue data alone, Q_TP·f_xy (a product) and Cwᵒ-vs-uptake-conductance are ridges → pinning transport absolutely needs an independent measurement (xylem sap / a pore-water probe).

### 6.2 Your own data tables (growth + pore-water)
Drive the model with two tables you enter (editable grids + CSV upload).

- **🌱 Growth table**: `day, root, stem, leaf, grain` — per-organ **FRESH weight** over time.
  - Units selectable: `g/hill` (default) · `kg/hill` · `g/m2` · `kg/ha` · `t/ha`.
  - The model's M is **per-hill fresh-weight mass**, so units are converted automatically.
- **💧 Pore-water table**: `day, Cwo` — the **absolute** soil-water PFAS concentration (µg/L) over time.
- **Compartment density ρ [kg/L, fresh]**: defaults root 1.0 · stem 0.30 · leaf 0.30 · grain 1.20 (editable).
  - The growth table is **mass** and the transport ODE is **mass-based** (no density prefactor), so density is a **mass↔volume bridge / reporting** quantity.
  - The app shows the "implied end-of-season organ volume (= fresh mass ÷ density)"; airy leaf/culm < 1, dense grain > 1.
- **Partial input allowed**: a growth table alone falls back to the sidebar Cwᵒ; a Cwᵒ table alone runs on the selected biomass driver.

### 6.3 Downloads (CSV / PNG)
- **Summary table (CSV)**: per-tissue model BAF / final concentration / observed BAF / (if any) measured BAF.
- **Full time series (CSV)**: `t, Cwo, Qtp, conc_*, M_*`.
- **Plant map (PNG)**: when `kaleido` (+ Chrome) is installed; otherwise a caption (CSV always works).

---

## 7. Input data formats (CSV)

Every table's first line is the header. **Unit convention**: time `day`, aqueous conc `µg/L`, tissue conc `µg/kg`, mass `kg`, flow `L/day`, BAF `L/kg`.

### 7.1 Driver CSV (HYDRUS / CSV drivers)
```
t,Cwo,Qtp,M_root,M_stem,M_leaf,M_grain
0,1.0,0.005,0.0001,0.0001,0.0003,0.0001
...
```
- Required: `t`, `Cwo`. Optional: `Qtp`, `M_root/M_stem/M_leaf/M_grain` (omitted → measured forcing / biomass driver).
- HYDRUS-1D mapping: `Cwo` ← root-zone node `Conc`, `Qtp` ← `vRoot`/`T_act`, `M_*` ← a plant growth sub-model (not HYDRUS).

### 7.2 Biomonitoring CSV (Biomonitoring)
```
tissue,conc,Cwo
root,0.49,1.0
straw,0.83,
grain,0.46,
```
- `tissue` (root/straw/stem/leaf/grain), `conc` (µg/kg), optional `Cwo` (µg/L; only one row needs it). BAF = conc / Cwᵒ.

### 7.3 Growth CSV (your own data tables)
```
day,root,stem,leaf,grain
0,0.05,0.02,0.03,0
80,1.0,5.5,4.0,3.0
150,1.1,7.0,4.5,12.0
```
- Per-organ **fresh weight** (units selected in the UI). `root` may be omitted (handled internally).

### 7.4 Pore-water Cwᵒ CSV (your own data tables)
```
day,Cwo
0,2.0
60,1.5
120,0.5
```

---

## 8. Reading the results

- **Concentration (µg/kg)**: micrograms of PFAS per kg of tissue (fresh weight). The default display in Simple mode.
- **Build-up factor BAF (L/kg)**: tissue conc ÷ soil-water conc — "how many times more concentrated than the water". 2 = twice.
  - In Simple mode it appears only in words ("build-up factor"); in Expert mode you see it on bars/maps directly.
- **Map colour**: hotter = higher concentration in that compartment. The colour scale is shared across time/compartments, so colours are comparable.
- **Where it goes**: usually **the roots hold the most**; how much reaches straw/grain is **congener-dependent** (short chains can put straw above root; long chains are root-dominated).
- **Grain formation gate**: the grain only takes up → accumulates PFAS after it forms (around flowering); the pre-formation period is not drawn (the organ physically does not exist yet).
- **Uncertainty (inverse)**: the spread of the posterior curve is the uncertainty. A wide 95% interval means the data did not pin down the exposure.
- **fresh vs dry**: the model is fresh-weight (fw). To compare with dry-weight (dw) data, use `C_dw = C_fw/(1−θ_fw)` (expert).

> [!IMPORTANT]
> Observed bars (e.g. Yamazaki) are **fixed measurements** — they do **not** move when you change the sidebar.
> The model was calibrated to reproduce that experiment at one **specific setting** (Model parametric, Cwᵒ=1,
> f_xy=W2 fit, E_m −120, season ~120 d); only there is the overlay a like-for-like match. Other settings/modes
> are a reference trend, not a calibrated match.

---

## 9. Scientific background & assumptions

Full derivation: `docs/pfas_rice_compartmental_model.{tex,pdf}`; the overview entry point is `docs/OVERVIEW_KR.md`.

- **Model skeleton**: an **ionizable-organic-compound (IOC) extension** of the Trapp/Brunetti **DPU (Dynamic Plant Uptake)** framework — a 4-compartment (root/stem/leaf/grain) dynamic ODE.
- **PFAS = permanent anion** (very low pKa, `f_d≈1`) → the neutral-compound Briggs/Kow partition core does not apply.
- **Hybrid root uptake `j_R`**: ionic electrodiffusion (GHK; inside-negative membrane ⇒ anion **exclusion**, `e^N≈107`) **+** a saturable Michaelis–Menten carrier.
- **Internal transport**: xylem advection (up) + phloem (grain is phloem-fed) + the **binding factor `B_k`** (`θ + f_prot·K_prot + f_PL·K_PL + f_cw·K_cw`, basis-A fresh weight, no density prefactor).
- **Root→shoot loading `f_xy`** (TSCF analog): the anion is retained in the root and translocates poorly; the ordering is congener-dependent.
- **Grain & leaf = terminal accumulators**: growth dilution → 0 at maturity ⇒ no steady state ⇒ a **dynamic model is essential** (final conc = time-integral / final mass).
- **Metabolism `γ≈0`** (PFAS recalcitrant), air exchange off.
- **Soil coupling (Method A, one-way)**: HYDRUS-1D/Phydrus → Cwᵒ(t), Q_TP(t); the plant ODE is solved in Python. HYDRUS is unmodified.

**Known limitations (stated honestly)**
- The demo W2-fit RMSE 0.029 is a **saturated reproduction** (params = observations), not predictive validation. The a-priori predictive error is log10 RMSE ~0.84–0.95 (long chains collapse).
- **Long-chain (C10–C12) shoot under-prediction**: lipid-facilitated loading / a 2-pool split / enhanced long-chain carrier capacity each help, but a PFDoDA residual remains (exploratory).
- `K_cw` (cell-wall partition) has no measured coefficient in the literature (placeholder). Ether/sulfonamide Koc is a gap.
- Detailed validation/refutation (sci-adk): see `docs/VALIDATION_KR.md`, `sci_adk_review/`.

---

## 10. Glossary

| Term | Plain meaning |
|---|---|
| PFAS | A family of long-lasting synthetic "forever chemicals" |
| Pore-water level (Cwᵒ) | How much PFAS is dissolved in the soil water around the roots [µg/L] |
| Build-up factor (BAF) | How many times more concentrated the PFAS is in tissue than in the soil water [L/kg] |
| Roots / Straw / Grain | Plant parts. Straw = stems + leaves; Grain = the edible brown rice |
| Concentration | Micrograms of PFAS per kg of tissue [µg/kg] |
| Congener | One specific PFAS (e.g. PFOA). Longer chains generally stick more |
| Uptake / translocation | How the chemical enters the roots and moves up into shoot and grain |
| Bayesian estimate | Working backwards from a measurement to the most likely cause, **with an uncertainty range** |
| f_xy (TSCF) | Root→xylem loading fraction (ease of shoot translocation) |
| eᴺ (anion exclusion) | The GHK anion-exclusion factor (the negative membrane suppresses anion uptake) |
| B_k | The binding factor (tissue conc = B_k × local free-water conc) [L/kg fw] |
| TF | Transfer factor = C_organ / C_root (Tang validation) |

---

## 11. Troubleshooting (FAQ)

- **SMILES mode won't appear / RDKit error**: `pip install rdkit` (or `-r requirements-structure.txt`). Otherwise use a curated congener.
- **"Run HYDRUS-1D (live)" is missing or asks me to build**: the FORTRAN engine (gfortran) is needed. `make` in `external/hydrus_source/source`, then `pip install phydrus`. The web session's SessionStart hook builds it automatically. Use another mode otherwise.
- **PNG download is disabled / shows a note**: the optional `kaleido` (+ Chrome) is required. CSV downloads always work.
- **The BAF doesn't change**: BAF = conc / Cwᵒ, so in the linear regime it is ~independent of Cwᵒ (changing the preset changes the concentration, not the BAF).
- **The Yamazaki bars don't move**: that's correct (fixed measurements). See the [§8 note](#8-reading-the-results).
- **My table reverts to the default**: rows with blank/non-numeric cells are dropped; at least 2 complete rows are needed.
- **Korean text shows as boxes (□)**: it renders correctly when a Korean font (e.g. Noto Sans KR) is available to the browser/system.
- **I want to change the language**: Simple = Korean, Expert = English; switch with the sidebar toggle.

---

## 12. Reproduce / validate / test

```bash
python reproduce_demo.py          # Yamazaki BAF full-ODE reproduction (W2 fit, log10 RMSE ≈ 0.029)
python reproduce_demo.py --rec    # monotone physical f_xy (a-priori; single-straw mismatch — the honest error)
python build_parameters.py        # rebuild params/parameters.json
pip install pytest && pytest      # full test suite (structure/mass-conservation/QSPR/calibration/API/plots/inverse/tables …)
```
- Validation docs: `docs/VALIDATION_KR.md`, `docs/VALIDATION_TANG2026_*_KR.md`, `validation/*.py`.
- App compute/figures are head-less tested: `tests/test_model_api.py`, `tests/test_plots.py`.

---

## 13. Cite / license / references

- **Cite**: PFAS–Rice Compartmental Uptake Model (an IOC extension of the Trapp/Brunetti DPU framework), 2026.
- **Source / docs**: <https://github.com/ccy5123/pfas-rice-model> · `docs/`
- **Key references**:
  - Yamazaki et al. 2023, *Environ. Sci. Technol.* **57**, doi:10.1021/acs.est.2c08767 (observed BAF)
  - Tang 2026, *J. Hazard. Mater.*, doi:10.1016/j.jhazmat.2025.141017 (out-of-sample TF)
  - Brunetti et al. 2019 *WRR* 10.1029/2019WR025432 · 2021 *ES&T* 10.1021/acs.est.0c07420 · 2022 *JHM* 10.1016/j.jhazmat.2021.127008
  - HYDRUS-1D 4.08 (LGPL-3.0; `external/hydrus_source`)
- **License**: see the repository `LICENSE` (the vendored HYDRUS-1D is LGPL-3.0).

> [!WARNING]
> Again: this tool is **for research and education**. Do not use it as the basis for regulatory, food-safety, or health decisions.
