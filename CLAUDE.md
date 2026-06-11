# CLAUDE.md — PFAS Rice Compartmental Uptake Model

> Persistent context for Claude Code. Read this first. Full math lives in `docs/`.

## 1. Goal
Mechanistic **four-compartment dynamic model** for PFAS bioaccumulation in rice
(*Oryza sativa*), built as an **ionizable-organic-compound (IOC) extension** of the
Trapp/Brunetti **Dynamic Plant Uptake (DPU)** framework, designed to couple with
**HYDRUS-1D** for the soil side.

## 2. Scientific summary (see `docs/` for the full derivation)
- PFAS = **permanently dissociated anion** (very low pKa, `f_d ≈ 1`). The neutral-compound
  Briggs/Kow partition core does **not** apply.
- Compartments: `root(1), stem(2), leaf(3), fruit/grain(4)`.
- **Root uptake `j_R` is hybrid**: ionic electrodiffusion (GHK; inside-negative membrane
  ⇒ anion *exclusion*, `e^N ≈ 107`) **+** saturable carrier (Michaelis–Menten). Net uptake
  requires the carrier to overcome electrostatic exclusion.
- Internal compartments exchange by **advection** (xylem up; phloem to grain) plus a
  **binding factor** `B_k = θ_k + f_prot·K_prot + f_PL·K_PL + f_cw·K_cw`
  (Briggs-consistent units, **NO density prefactor**).
- Grain is **phloem-fed**; the weak-acid pH **ion-trap does NOT apply** (`f_n ≈ 0`) ⇒ phloem
  loading is carrier/channel (`L_Ph`), not a pH trap.
- **Grain and leaf are terminal accumulators**: the only sink is growth dilution, which → 0
  at maturity ⇒ no bounded steady state ⇒ final conc = time-integral / final mass. A
  **dynamic** model is therefore essential.
- Metabolism `γ_k ≈ 0` (PFAS recalcitrant). Air exchange off (`K_AW ≈ 0`).

Model report: `docs/pfas_rice_compartmental_model.{tex,pdf}`
Corrected neutral DPU base: `docs/dpu_model_summary_corrected.tex`

## 3. Repo layout
```
.
├── CLAUDE.md  README.md  requirements.txt
├── reproduce_demo.py                 # entry point: Yamazaki BAF via full ODE (W2 fit)
├── build_parameters.py               # (re)assembles params/parameters.json from source tables
├── src/
│   ├── pfas_rice_plant_module_4pool.py       # basis-A 4-compartment ODE (CANONICAL core)
│   ├── pfas_rice_plant_module_4pool_surf.py  #  + K_surf (Fe/Mn-plaque dead-end pool)
│   ├── pfas_rice_plant_module_5pool.py       #  + explicit lignin pool
│   ├── pfas_rice_plant_module_nstem.py       # N serial stem segments (multi-height; GAP-B fix)
│   ├── pfas_rice_plant_module.py             # import alias → 4pool_surf (basis-A); legacy name
│   ├── soil_paddy.py                         # Freundlich soil → C_w^o(t) (legacy redox sign)
│   ├── soil_paddy_redox_corrected.py         # W3-corrected redox (dilution+leaching; USE THIS)
│   ├── soil_hydrus.py                        # REAL HYDRUS-1D soil → C_w^o(t),Q_TP(t) via phydrus (Method A WIRED)
│   ├── calibration.py                        # Tier-1 calibration (scipy)
│   └── literature_params.py                  # literature QSPRs/anchors (cited) + Kim2019 BAF
├── params/                           # parameters.json (CANONICAL) + source CSVs (Bk, f_xy, Kcw, ...)
├── data_obs/                         # observed BAF/TF (Yamazaki, Li2025) + yamazaki_stem_height.csv
├── validation/                       # S6 + nstem + hydrus_coupled_run reproduction scripts + figures
├── docs/
│   ├── pfas_rice_compartmental_model.tex / dpu_model_summary_corrected.tex
│   ├── DELIVERABLE_GAP_A_Kcw.md / DELIVERABLE_GAP_B_fxy.md / theory_anchor.tex / H8_handoff_S6_final.md / sources.csv
│   └── literature_db/                # curated parameter DB (.xlsx + per-sheet .csv) + raw_si/ SI extractions
├── external/hydrus_source/           # git submodule → github.com/phydrus/source_code
├── data/                             # (gitignored)
└── tests/                            # pytest (76): plant, soil, hydrus, calibration, literature params
```

## 4. Coupling strategy
- **Method A — loose, one-way (CURRENT; now WIRED to a real HYDRUS run).** HYDRUS-1D/Phydrus →
  `C_w^o(t)`, `Q_TP(t)`; the plant ODE is solved in Python (`src/pfas_rice_plant_module.py`).
  No FORTRAN edits. Interface = the three arrays in `PlantInputs` (`Cwo`, `Qtp`, `M`).
  `src/soil_hydrus.py` builds & runs the compiled HYDRUS-1D engine (via `phydrus`) for a paddy
  scenario per congener (Kd from the C3 `Koc` QSPR) and returns BOTH the pore-water trajectory
  `C_w^o(t)` and the root water uptake `Q_TP(t)` → `inputs_from_hydrus()` → `PlantInputs`. The
  soil run is driven by the MEASURED transpiration (`forcing_rice.transpiration_mm_d`), so HYDRUS's
  actual uptake `vRoot` carries the measured crop-physiology shape (+ soil-water-stress feedback);
  `qtp_from_hydrus=True` (default) reproduces `forcing_rice.Q_TP` to <1% when unstressed (consistency
  test) and only diverges under water limitation. See `validation/hydrus_coupled_run.py`.
- **Method B — tight (FUTURE).** Modify `external/hydrus_source` (HYDRUS-1D FORTRAN):
  replace/augment the **root solute-uptake routine** with `j_R`, add the plant module,
  rebuild via `makefile`. `external/hydrus_source/source_mcmc/` provides Bayesian
  calibration machinery.

## 5. Parameter tiers (calibration design)
- **Tier 0** inputs/known: `M_k(t), Q_TP(t), C_w^o(t), N(E,z), f_d, γ_k≈0, T_C,Ph`.
- **Tier 1** BAF-identifiable (lumped): `B_k`, `g_in/g_out`, `f_xy` (root→xylem loading/TSCF),
  `Π = Q_Phl·L_Ph/Q_TP`, `φ`.
- **Tier 2** need inhibitor/kinetic data: separate `P_d^eff` (channel) vs `V_max` (carrier);
  influx vs efflux asymmetry.
- **Tier 3** QSPR/measurement (chain-length resolved): `K_prot, K_PL, K_cw, L_Ph, a_R`.
- **Identifiability**: BAF data constrain only the lumped influx conductance
  `g_in = a_R·P_d^eff + carrier clearance` — channel vs carrier are **not separable** from
  BAF alone (need inhibitor experiments).

## 6. Current status
- Derivation + LaTeX docs: complete (TSCF loading factor + mass-conserving phloem added;
  `ρ_k` binding bug fixed). **PDFs are not in the repo — rebuild with pdflatex where available.**
- Python module: runs (BDF stiff solver); reproduces the structural results
  (anion exclusion, terminal-sink accumulation, binding, TSCF-limited translocation).
- Test suite: `tests/test_plant_module.py` (pytest) locks in the structural invariants
  and exact mass conservation (`pip install pytest && pytest`).
- **RESOLVED (was KNOWN ISSUE)**: the terminal-sink runaway is fixed *structurally* by the
  root→xylem loading factor `f_xy` (TSCF, assumption A2): the anion is retained in the root
  and translocates poorly, so leaf/grain no longer out-accumulate the root. The demo now
  reproduces **root > straw > grain** (straw = mass-weighted stem+leaf). Also closed a phloem
  mass-conservation leak (leaf now exports the full `(1+φ)·Q_Phl·C_Phl`). **Demo BAFs remain
  illustrative, NOT calibrated** — real Tier-1 calibration vs data is still task #4.
- **Soil side (task #3)**: `src/soil_paddy.py` adds a Freundlich paddy sorption sub-model
  (`S=K_F·C_w^n`, redox-dependent `K_F`) that inverts a total soil inventory to the
  pore-water `C_w^o(t)`, plus `load_inputs_csv` to drop in real HYDRUS-1D/Phydrus output.
- **REAL HYDRUS-1D run wired (task #3)** — `src/soil_hydrus.py`: the submodule HYDRUS-1D 4.08
  engine is now **compiled** (gfortran; `external/hydrus_source/source/hydrus`) and driven through
  `phydrus` to produce a genuine pore-water `C_w^o(t)` and root water uptake `Q_TP(t)` for a
  one-season paddy (clean-water flooding → drainage), per congener via a **linear Kd** isotherm
  (`Kd = Koc·f_oc`; Freundlich n<1 makes the solute solver diverge at the c→0 clean-water boundary,
  so linear Kd is used — full congener-resolved retardation R=1+ρKd/θ is retained). `inputs_from_hydrus()`
  normalises the series to season-mean exposure and returns `PlantInputs`. **Result** (`validation/
  hydrus_coupled_run.py`): the pore water is strongly **congener-dependent** — weakly-sorbed short
  chains (Kd≈0.01–0.15) leach to near-zero during flooding so the constant-`Cwo` placeholder
  **over-predicts grain/straw BAF ~2–4×** (PFBA grain 2.07→0.43), while strongly-sorbed long chains
  (Kd≳7) stay buffered (BAF≈unchanged). Tests skip when the exe/phydrus is absent. **Remaining**:
  anoxic/flooded sorption + a real field flooding schedule + the user's site soil/loading.
- **Calibration (task #4)**: `src/calibration.py` fits Tier-1 params to observed tissue
  BAFs (log-space weighted least squares, scipy; box bounds; optional global DE). Validated
  by `synthetic_recovery` (recovers known Tier-1 params, incl. under noise). NOTE: tighten
  the finite-diff step (`diff_step≈1e-2`) so the gradient clears the ODE solver's tolerance
  floor. Real fit pending the user's BAF data (`load_baf_csv`).
- **Literature database (task #2 enabler)**: `docs/literature_db/` holds the curated empirical
  parameter database (xlsx + per-sheet CSV; categories C1–C6 + source shortlist + gap analysis)
  plus `raw_si/` (per-congener tables extracted from the cited papers' SI). `src/literature_params.py`
  transcribes the **verified** pieces — soil `Koc(chain length)` QSPR (Higgins & Luthy +0.55/CF₂,
  +0.23 sulfonate; anchored on Milinovic PFOA/PFOS/PFBS), `f_d` from pKa (Goss 2008), rice root
  `E_m` (Wang 1994), and the **MEASURED per-congener `K_PL`** (Chen 2025 K_MW Table S5, L/kg lipid,
  cross-checked vs Droge 2019 SSLM) and **`K_prot`** (Zhou 2025 **Table 1** dialysis `K_prow`: soy
  protein isolate = the plant/grain analog, BSA = animal reference) — into builders
  (`literature_compound`, `literature_environment`, `literature_paddy_soil`). Each value carries a
  citation + `DOI_status`. NOTE: the Chen HSA `K_D` / Zhou BSA `K_A` are *binding constants*; the
  single-site `K_D`→partition route overestimates ~50× vs the dialysis `K_prow`, so dialysis is used
  and `k_prot_albumin` is reference-only. **Still placeholder**: `K_cw` (no coefficient exists in the
  literature); transport params (`f_xy, L_Ph, kappa_d, Vmax/Km`) remain fitted (Tier-1/2).
- **Real Tier-1 calibration (task #4)**: Kim 2019 (`docs/literature_db/raw_si/kim2019_*`) gives
  per-congener brown-rice (grain) BAF paired with paddy pore water. `literature_params.kim2019_grain_baf()`
  exposes it; the demo fits `L_Ph` to the PFOA grain BAF (0.07 → 4.43 L/kg, `L_Ph≈0.84`). The
  measured binding keeps `root > straw > grain` (delivery-limited), but the *grain BAF* now matches
  data. **Limitation**: Kim is grain-only, so `f_xy` (root→shoot) is unconstrained — full
  compartment-resolved TF (root/straw) is a DB gap (greenhouse time-series needed).
- **Integrated advanced workstream (basis-A + GAP A/B + nstem)**: merged the consolidated
  parameter package (`params/parameters.json` + `params/*.csv`, `data_obs/`, `validation/`,
  GAP deliverables in `docs/`, the basis-A plant modules `*_4pool[_surf]/_5pool`, and
  `soil_paddy_redox_corrected`). `pfas_rice_plant_module` is now an **alias to the basis-A
  4pool_surf** core. Key honest-status corrections from the review: (a) `reproduce_demo.py`'s
  log10 RMSE 0.029 is a **saturated W2 fit** (3 transport params/3 obs per congener) — reproduction
  is guaranteed, NOT predictive validation; (b) the empirical ordering is **congener-dependent**
  (Yamazaki: short-chain straw≫root, long-chain root>straw) — `root>straw>grain` is NOT universal
  under basis-A; (c) **GAP B is shape-resolved, not closed** — see task #6.
- **Multi-height stem (task #6)** — `src/pfas_rice_plant_module_nstem.py`: `NStemModel` (equilibrium)
  + `NStemKineticModel` (finite radial `k_rad`); both mass-conserving. Driven by the MEASURED
  forcings (`src/forcing_rice.py` Q_TP from Kumari/NayHtoon; `src/growth_rice.py` M_s from ORYZA
  IR72), the multi-height stem **reproduces the short/mid-chain UPWARD gradient but NOT the
  long-chain reversal** (PFBA top/bot 7.4 → PFUnDA 0.66; model stays ~2.4–4.9 for all k_rad). The
  reversal needs **irreversible/hysteretic high-B sorption** — set aside, fully recorded in
  `docs/nstem_gradient_exploration.md`. NOTE: the earlier "monotone f_xy reproduces the gradient"
  claim was a **placeholder-biomass artifact** (real ORYZA biomass moves the crossover `B* ~
  Q_s/(M_s·μ_s)` above the congener range).
- **f_xy absolute scale (task #7)**: measured `Q_TP(t)` (`forcing_rice`, peak ~0.10 L/d/hill, T/ET=0.42)
  and `M_s(t)` (`growth_rice`, ORYZA IR72, HI~0.53) are built. The absolute f_xy is pinned via the
  **aggregate** root/straw/grain BAF (not the within-stem gradient) — see `validation/`.

## 7. Build & run
- `pip install -r requirements.txt`
- **Main reproduction**: `python reproduce_demo.py` (Yamazaki BAF, W2 fit, RMSE≈0.029);
  `--rec` uses the monotone f_xy. Rebuild params: `python build_parameters.py`.
- Plant demo: `python src/pfas_rice_plant_module_4pool_surf.py` (N, B_k, BAFs; saves `pfas_rice_demo.png`).
- Multi-height stem: `python validation/nstem_gradient_check.py` (stem-gradient direction vs Yamazaki).
- Soil → plant (analytic): `python src/soil_paddy.py` (legacy) / use `soil_paddy_redox_corrected` for redox.
- **Soil → plant (REAL HYDRUS-1D)**: build the engine once, then run the coupling:
  ```
  git submodule update --init external/hydrus_source
  cp external/hydrus_source/makefile external/hydrus_source/source/
  (cd external/hydrus_source/source && make)          # needs gfortran
  pip install phydrus
  python src/soil_hydrus.py                            # per-congener pore-water summary
  python validation/hydrus_coupled_run.py             # full soil→plant + figure/CSV
  ```
- Calibration: `python src/calibration.py`; Literature params: `python src/literature_params.py`.
- Tests: `pip install pytest && pytest` (76 passing, 1 skipped; HYDRUS tests auto-skip if the exe is unbuilt).
- FORTRAN (Method B): init submodule (`git submodule update --init`), then follow
  https://phydrus.readthedocs.io/en/latest/getting_started/compilation.html
  (gfortran + `makefile` / `make.bat`). NOTE: the top-level `makefile` lists the `.FOR` files
  without a path, so build from inside `source/` (copy the makefile in, as above).

## 8. Conventions
- Units: time **day**; aqueous conc **µg/L**; tissue conc **µg/kg**; mass **kg**;
  flow **L/day**; `B_k` in **L/kg fw** (`C_k = B_k · C_w,k`).
- **Binding = basis A (fresh weight)**: `B_k = θ_fw + (1−θ_fw)·(f_prot·K_prot + f_PL·K_PL + f_cw·K_cw)`.
  `θ_fw` = fresh-weight water fraction; `f_*` = **dry-weight** mass fractions; `K_*` in L/kg pool-dw.
  The `(1−θ_fw)` factor is a **dry→fresh conversion** (mandatory; the legacy naive `θ+Σf·K` over-states
  B_k ~3×) — it is NOT the old dimensionally-wrong `ρ_k` density prefactor (still absent). Compare to
  dw-reported data via `C_dw = C_fw/(1−θ_fw)`. `f_cw` = whole cell wall (poly+lignin), K = `K_cw_wholecw`.
- `f_xy` ∈ (0,1] is the root→xylem loading factor (TSCF analog): only `f_xy·C_1/B_1` enters the
  ascending xylem (`f_xy=1` = unrestricted DPU; `f_xy≪1`, monotone-decreasing in chain length, for anions).
  NOTE: it does **not** yield a universal `root>straw>grain` — the ordering is **congener-dependent**
  (short chains: straw>root; long chains: root>straw), matching Yamazaki.
- Symbols map 1:1 to `docs/pfas_rice_compartmental_model.tex` (`j_R, B_k, N, f_xy, L_Ph, ...`).

## 9. Next tasks (prioritized)
1. ~~Physical realism of terminal compartments~~ **DONE** — added the root→xylem loading
   factor `f_xy` (TSCF) + mass-conserving phloem; demo reproduces `root > straw > grain`;
   regression tests in `tests/`. (Calibrating `f_xy`/`L_Ph`/`B_k` to data is task #4.)
2. **Tier-3 QSPR** for `K_prot`, `K_PL` (chain-length descriptors) to populate `B_k`
   **MOSTLY DONE** (`src/literature_params.py` + `docs/literature_db/raw_si/`): **measured
   per-congener** `K_PL` (Chen 2025 K_MW, vs Droge 2019) and `K_prot` (Zhou 2025 Table 1 dialysis
   `K_prow` — soy protein isolate for plant tissues, BSA for animal) are extracted and wired into
   `B_k`. **Remaining**: only a quantitative `K_cw` (no coefficient exists in the literature — batch
   sorption to rice root cell-wall fractions, pectin/hemicellulose).
3. **Freundlich paddy soil sorption** **DONE** (`src/soil_paddy.py`); literature `Koc`→`K_F`
   parametrization in `src/literature_params.py`. **Real HYDRUS-1D run now WIRED** (`src/soil_hydrus.py`,
   `validation/hydrus_coupled_run.py`): the compiled engine produces a genuine per-congener `C_w^o(t)`
   that drives the plant ODE (short chains leach → constant-`Cwo` over-predicts grain BAF ~2–4×).
   **Remaining**: anoxic/flooded sorption (DB gap), a real field flooding schedule, and the user's
   site-specific soil/loading. HYDRUS now also supplies `Q_TP(t)` by default (`qtp_from_hydrus=True`),
   driven by the measured `forcing_rice` transpiration and reproducing it to <1% when unstressed.
4. **Tier-1 calibration machinery** **DONE** + **first real fit done**: **Kim et al. 2019**
   (Korean paddy, paired pore-water/soil/brown-rice, `10.1016/j.scitotenv.2019.03.240`) is wired in
   (`kim2019_grain_baf()`); the demo fits `L_Ph` to the PFOA grain BAF (→ matches 4.43 L/kg).
   **Remaining**: a chain-length series fit and a full compartment-resolved fit — Kim is grain-only,
   so `f_xy` (root→shoot) needs root/straw tissue data (DB gap; greenhouse time-series needed).
5. **Literature parameter DB + measured `B_k`** **DONE** (`docs/literature_db/`, `src/literature_params.py`):
   curated C1–C6 DB + `raw_si/` extractions; measured `K_PL`/`K_prot`/`K_cw` wired into basis-A `B_k`.
6. **Multi-height stem (GAP-B fix)** **DONE (structural)** — `src/pfas_rice_plant_module_nstem.py` +
   `validation/nstem_gradient_check.py`: monotone f_xy reproduces the PFCA stem gradient.
7. **Measured `Q_TP(t)` / `M_s(t)`** → pin the f_xy absolute scale + gradient crossover and run the
   full compartment-resolved fit (currently structural/direction only; placeholder transpiration ~5× high).
   Candidate value source flagged by the user: **Tang 2026 JHM (`10.1016/j.jhazmat.2025.141017`)**.
8. **PFSA-specific transport term** **DONE (sign pinned)** — the headgroup offset on `f_xy` is
   confirmed and quantified: in BOTH Tang 2026 (paddy, PFOS/PFOA TF 0.26) and Yamazaki 2023 (0.43),
   **PFSA translocates LESS** than the CF2-matched PFCA, so `f_xy(PFSA) = f_xy(PFCA)·exp(−1.1)`
   (refines the placeholder `exp(−1.5)`; sign was "uncertain"). Wired as `literature_params.f_xy_headgroup`
   + `FXY_HEADGROUP_LN_OFFSET`; `params/parameters.json` PFSA `f_xy_recommended` rescaled (build via
   `Bk_table_S5.csv`). Ether (GenX) factor `exp(−0.7)` documented (Tang, provisional; not in the core
   12). **Note**: this is distinct from the *long-chain PFCA* shoot mechanism (the f_xy-fit U-shape /
   PFDoDA un-capturable), which remains open (hysteretic sorption — `docs/nstem_gradient_exploration.md`).
9. (Later) **Method B** tight coupling in `external/hydrus_source`.

## 10. Gotchas / external dependencies
- DPU module source is **not public** (author request only). The ionizable extension
  (Brunetti 2022) is **not in the HYDRUS distribution** → reimplement from the papers.
- `phydrus/source_code` is HYDRUS-1D **4.08** (older than official 4.17), **LGPL-3.0**, and
  is the **base soil engine only** (no DPU/PFAS/ionizable modules).
- Key references: Brunetti 2019 *WRR* `10.1029/2019WR025432`; 2021 *ES&T*
  `10.1021/acs.est.0c07420`; 2022 *J. Hazard. Mater.* `10.1016/j.jhazmat.2021.127008`.
