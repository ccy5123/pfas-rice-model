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
├── CLAUDE.md  README.md  requirements.txt  requirements-app.txt
├── reproduce_demo.py                 # entry point: Yamazaki BAF via full ODE (W2 fit)
├── build_parameters.py               # (re)assembles params/parameters.json from source tables
├── app.py                            # Streamlit visualization tool (plant/soil map + 4 input modes)
├── src/
│   ├── pfas_rice_plant_module_4pool.py       # basis-A 4-compartment ODE (CANONICAL core)
│   ├── pfas_rice_plant_module_4pool_surf.py  #  + K_surf (Fe/Mn-plaque dead-end pool)
│   ├── pfas_rice_plant_module_5pool.py       #  + explicit lignin pool
│   ├── pfas_rice_plant_module_nstem.py       # N serial stem segments (multi-height MIXER; Yamazaki gradient)
│   ├── pfas_rice_plant_module_nstem_leaf.py  # N stem segs + explicit leaf (transpiration deposition+RETENTION; Tang over-translocation fix)
│   ├── pfas_rice_plant_module.py             # import alias → 4pool_surf (basis-A); legacy name
│   ├── soil_paddy.py                         # Freundlich soil → C_w^o(t) (legacy redox sign)
│   ├── soil_paddy_redox_corrected.py         # W3-corrected redox (dilution+leaching; USE THIS)
│   ├── soil_hydrus.py                        # REAL HYDRUS-1D run via phydrus → Cwo(t),Qtp(t) (Method A; wired + app live mode)
│   ├── calibration.py                        # Tier-1 calibration (scipy)
│   ├── literature_params.py                  # literature QSPRs/anchors (cited) + Kim2019 BAF
│   ├── model_api.py                          # UI-agnostic wrapper: simulate(), simulate_from_smiles(), driver/soil/biomon helpers
│   ├── pfas_structure.py                      # SMILES (structure) → Compound adapter (RDKit; read-across + QSPR)
│   ├── plots.py                              # Plotly builders: fig_plant_schematic (colormap), drivers, ...
│   ├── forcing_rice.py                       # measured transpiration Q_TP(t) (FAO-56 dual-Kc; Kumari2022 + NayHtoon2018)
│   ├── growth_rice.py                        # ORYZA IR72 partitioning on a logistic → organ M_s(t) (DEFAULT biomass driver)
│   ├── oryza_growth.py                       # MECHANISTIC ORYZA2000 Level-1 carbon balance → weather-responsive M_s(t) (opt-in; drivers=/weather=)
│   └── measured_biomass.py                   # ingest a MEASURED per-organ biomass table → M(t) driver (units→kg/hill; Tang etc.)
├── examples/                         # ready-to-load CSVs for app.py (HYDRUS drivers + biomonitoring)
├── params/                           # parameters.json (CANONICAL) + source CSVs (Bk, f_xy, Kcw, ...)
├── data_obs/                         # observed BAF/TF (Yamazaki, Li2025) + yamazaki_stem_height.csv
├── validation/                       # S6 + nstem + hydrus_coupled_run reproduction scripts + figures
├── docs/
│   ├── OVERVIEW_KR.md                # ★ 종합 진입점: 기능·검증·데이터공백·필요실험·notation 표 (+모식도)
│   ├── pfas_rice_compartmental_model.tex / dpu_model_summary_corrected.tex
│   ├── DELIVERABLE_GAP_A_Kcw.md / DELIVERABLE_GAP_B_fxy.md / theory_anchor.tex / H8_handoff_S6_final.md / sources.csv
│   ├── visualization_tool.md         # app.py guide: plant/soil map, 4 modes, HYDRUS I/O, biomonitoring
│   └── literature_db/                # curated parameter DB (.xlsx + per-sheet .csv) + raw_si/ SI extractions
├── external/hydrus_source/           # VENDORED HYDRUS-1D 4.08 source (de-submoduled from phydrus/source_code; binary gitignored)
├── .claude/                          # SessionStart hook (hooks/session-start.sh): web deps + HYDRUS engine build
├── data/                             # (gitignored)
└── tests/                            # pytest (173 collected): plant, soil, hydrus, calibration, lit params, API (+two-pool, cwo_profile, k_leach), plots, structure(SMILES), oryza, measured-biomass, bayesian-inverse

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
  is guaranteed, NOT predictive validation. The genuine **a-priori predictive error** (theory/QSPR
  monotone f_xy, NOT fit) is **log10 RMSE ≈0.84** (single-straw, `reproduce_demo.py --rec`) /
  **≈0.95** (redistributed-shoot, `validation/apriori_prediction.py`) — straw 6–40× off, long chains
  collapse; i.e. the model does NOT predict out-of-sample. Adjudicated by the sci-adk rigor review
  (`sci_adk_review/FINDINGS.md`: hyp-yamazaki **REFUTED**). (b) the empirical ordering is **congener-dependent**
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
- **Tang over-translocation fix (redistributed shoot)** — `src/pfas_rice_plant_module_nstem_leaf.py`
  (`NStemLeafModel`; `model_api.simulate_nstem_leaf`): the Tang 2026 OOS check flagged the single-straw
  core's **empty stem (pass-through) + leaf-sink runaway** (leaf held ~81% of the plant burden). Fixed by
  resolving the stem into N segments AND **applying transpiration deposition+RETENTION to every shoot organ
  (not just the leaf)** — each organ retains its own transpired solute (a partial terminal), so the shoot
  burden is redistributed root→stem→leaf→grain. Two crop-architecture levers: `stem_transp_frac`,
  `retention` (default 0.45/0.6, NOT point-fit to Tang); mass-conserving (sole source `M_root·j_R`;
  `tests/test_nstem_leaf.py`). **Result** (`validation/tang2026_nstem_validation.py`,
  `docs/VALIDATION_TANG2026_NSTEM_KR.md`): the shoot **tissue PATTERN is improved** (dw-corrected shape RMSE
  0.85→0.39; PFOA stalk 0.02→0.75, leaf 2.70→0.93; leaf burden 81%→30%, stalk 1%→29%) — but the GRAIN stays
  structurally UNDER (PFOA endosperm 0.11 vs Tang 0.95; not closable by L_Ph/lipid — see
  `docs/tang2026_grain_units_exploration.md`, the fresh/dry units fix). NOTE the earlier "shape 0.84→0.11 /
  grain 0.41→0.93 cured" figures were a fresh-vs-dry artifact (model fw TF vs Tang dw TF), now corrected.
  **Then the across-congener absolute LEVEL was calibrated — the lever is `f_xy`, NOT `B_root`**: `B_root`(PFOS)=49
  is CONFIRMED by Yamazaki root data (PFOS root BAF 5.93 ≈ 12× PFOA 0.49) so it is correct; the residual traces to
  (i) the monotone `f_xy`(PFOS)=0.013 OVER-penalizing PFSA (the head-group exp(−1.1) offset) — Yamazaki's own W2 fit
  needs 0.142, and a mass-balance argument confirms 0.013 under-delivers; (ii) the GenX provisional `f_xy`=0.233
  (short-chain-PFCA × ether offset) being ~18× too high. Recalibrating `f_xy` (PFOS → W2 0.142 = independent
  Yamazaki; GenX → 0.013 = Tang, no independent data) drops **overall RMSE (dw) 1.53 → 1.20 (structure) → 0.71
  (f_xy; grain-limited)**, stalk/leaf within order-of-magnitude (grain remains the structural floor). The calibrated f_xy is applied as an **override in the validation
  only** — `params/parameters.json` is UNCHANGED (provenance preserved); follow-up is to re-fit the monotone PFSA
  head-group offset + an ether-PFAS QSPR for GenX (docs §6). COMPLEMENTARY to `nstem` (mixer, Yamazaki within-stem
  gradient): nstem_leaf uses RETENTION for the Tang stalk/leaf/grain split. Default model unchanged (4pool_surf);
  opt-in module.
- **f_xy absolute scale (task #7)**: measured `Q_TP(t)` (`forcing_rice`, peak ~0.10 L/d/hill, T/ET=0.42)
  and `M_s(t)` (`growth_rice`, ORYZA IR72, HI~0.53) are built. The absolute f_xy is pinned via the
  **aggregate** root/straw/grain BAF (not the within-stem gradient) — see `validation/`.
- **Visualization tool (`app.py` + `src/model_api.py` + `src/plots.py`)**: Streamlit dashboard whose
  headline is the **plant + soil accumulation map** — a rice plant (fibrous roots in the paddy soil,
  arching culms, long leaf blades, drooping grain panicles) with each compartment filled by a heat
  **colormap** of its concentration/BAF (`plots.fig_plant_schematic`), a season **day slider / ▶ animate**
  to watch the build-up, plus drivers / soil-profile / isotherm / chain / compare tabs. Covers **four
  exposure modes** via `simulate(..., drivers=…)`: (1) parametric, (2) **HYDRUS/Phydrus CSV** (`t,Cwo,Qtp,M_*`
  → `load_driver_csv`/`drivers_from_arrays`), (3) **soil inventory** (Freundlich inversion,
  `pore_water_from_inventory`), (4) **biomonitoring** (measured tissue conc, no HYDRUS — `baf_from_measurement`).
  `model_api`/`plots` are UI-agnostic + head-less-tested (`tests/test_model_api.py`, `tests/test_plots.py`);
  bundled `examples/` CSVs auto-load. HYDRUS-1D input/output mapping + the biomonitoring path are documented
  in the app's **About** tab and `docs/visualization_tool.md`.
- **Live HYDRUS-1D coupling (`src/soil_hydrus.py`)**: the **real HYDRUS-1D engine** (built from the
  `external/hydrus_source` submodule, gfortran) is driven through **`phydrus`** to run a one-season paddy
  model (Richards + advection-dispersion + **linear Kd** + root uptake) → congener-dependent pore water
  `Cwo(t)` (short chains leach under flooding, long chains buffer; verified: PFBA Cw→0.01, PFOA→0.47,
  PFDoDA→1.00) and actual root uptake `Q_TP(t)`. Per-congener Kd from the C3 Koc(chain-length) QSPR
  (`literature_params.koc`). Wired into the app as the 5th **"Run HYDRUS-1D (live)"** mode via
  `model_api.hydrus_drivers`/`hydrus_available` (graceful fallback when the engine/phydrus are absent);
  `tests/test_soil_hydrus.py` skips the engine tests when unbuilt. Still **Method A** (one-way; HYDRUS
  unmodified). Originally implemented on branch `claude/epic-knuth-npt0cy`; the soil piece is cherry-picked here.

- **Structure (SMILES) input — parameterise ANY PFAS (`src/pfas_structure.py`)**: the "option-3"
  front end that lets a **chemical structure** be the model input, not only the curated 13 congeners.
  RDKit parses the SMILES → structural descriptors (`n_perfluoroC`, `head_group` via SMARTS, `n_ether_O`,
  `n_CF3`, `branched`, MW/formula, `is_linear`) → a `Compound` by **(1) MEASURED read-across** when the
  (canonical) structure matches a curated congener (uses `params/parameters.json` exactly — a SMILES-built
  PFOA reproduces the named PFOA) **or (2) the literature_params QSPR** for a novel structure (per-CF2 slope
  + head-group offset; ether/sulfonamide flagged PROVISIONAL). Binding (`K_PL/K_prot/K_cw`) + speciation
  (`f_d` from head-group pKa) come from structure; **`f_xy` is NOT structure-derivable** — curated monotone
  for knowns, PFCA-series interpolation × head-group offset for novels (provisional). `model_api.simulate_from_smiles()`
  runs the full ODE (delegates to the canonical path for knowns; injects a custom record via the new
  `simulate(..., record=)` arg for novels) and returns the usual dict + `descriptors` + `provisional`.
  Sulfonamides/neutral species are detected and flagged (violate the permanent-anion `f_d≈1` assumption).
  RDKit is **optional** (`requirements-structure.txt`); `tests/test_pfas_structure.py` (23) skips when absent.
  Docs: `docs/structure_input.md`.
- **Ether fragment QSPR term (`literature_params.k_pl`/`koc`)**: `koc`/`k_pl` are now group-contribution —
  `k_pl` adds a per-ether-O term `KPL_ETHER_LOG_OFFSET = -0.49 log` **anchored on the GenX measurement**
  (Chen2025 K_MW 117.5 vs the CF2-only QSPR at nPFC=5 → −0.49; matches "ether REDUCES K_MW"; provisional,
  single anchor). So a novel ether-PFCA (ADONA-type) gets a reduced K_PL, not the carboxylate value. `koc`
  now accepts `ether`/`sulfonamide` head groups (was a ValueError) but `KOC_ETHER_LOG_OFFSET = 0` is an
  explicit **GAP** (no measured ether/sulfonamide soil Koc in the DB; the GenX BCF over-prediction was fixed
  by the f_xy recalibration, not Koc). Wired into `pfas_structure` (novel ethers use the ether term).
  Tests in `test_literature_params.py` (ether term reproduces GenX; koc graceful). Remaining: sulfonamide
  K_PL slope + ether/sulfonamide Koc need data (docs/structure_input.md §Next steps).
- **Mechanistic ORYZA2000 biomass driver (`src/oryza_growth.py`)**: a Python re-implementation of the
  **ORYZA2000 / ORYZA(v3) potential-production (Level-1) carbon balance** — SUCROS astronomy → Gaussian
  day×canopy gross CO₂ assimilation → maintenance+growth respiration → DVS-driven partitioning → SLA-based
  LAI (juvenile RGRL + senescence) → grain fill — so per-organ `M_s(t)` **responds to radiation/temperature**
  instead of the imposed logistic in `growth_rice`. NOT the IRRI binary (Windows exe needing a full weather/
  crop deck; gfortran/`pyoryza` unavailable here) — it is the published Level-1 equation set (Bouman & van Laar
  2006; Goudriaan & van Laar 1994 SUCROS) re-coded, with IR72 standard-set parameters anchored so the potential
  run reproduces the IR72 field anchors (flowering ~day 66, maturity ~116, LAImax 6.6, HI 0.46, shoot scaled to
  1740 g/m²). `oryza_drivers(congener)` returns a `model_api.simulate(drivers=…)` dict (wired via the same
  `drivers=` extension point as HYDRUS); `organ_biomass_oryza(t)` gives kg/hill per organ; `weather=` overrides
  the built-in climatology with a real series. Validation (`validation/oryza_growth_validation.py`) contrasts it
  with `growth_rice` and propagates BOTH biomass drivers through the PFAS ODE: the mechanistic biomass (leaf
  senescence + stem retention) **raises short-chain straw/grain BAF ~40-70%** (e.g. PFBA grain 2.07→3.53) but
  leaves the root-dominated long chains ~unchanged. `tests/test_oryza_growth.py` (6). Opt-in; the canonical path
  (`growth_rice`) is unchanged. Candidate next step: drive it with the measured `M_s(t)`/weather to pin the f_xy
  absolute scale (task #7). **Provenance note**: `oryza_growth.py` + `tests/test_oryza_growth.py` + `validation/
  oryza_growth_validation.py` were *first actually committed* in commit d1f5339 — this §6 description previously
  predated the code (doc-ahead-of-code); they are now in sync (verified by the doc↔code audit below).
- **Measured-biomass ingestion + Tang 2026 TF f_xy re-calibration (this session)**:
  - `src/measured_biomass.py` (+ `examples/measured_biomass_template.csv`, `tests/test_measured_biomass.py`): ingest a
    MEASURED per-organ biomass table → `M(t)` driver (units g/plant·t/ha·g/m²·… → kg/hill; interpolate; optional
    root:shoot reconstruction; pairs with `forcing_rice.Q_TP`). The data-grounded alternative to `growth_rice`/`oryza_growth`.
  - `model_api.simulate_nstem_leaf(biomass_fn=…)`: the redistributed-shoot model can now be driven by the mechanistic
    ORYZA biomass (default still `growth_rice`).
  - **Tang 2026 finding (key, condition-specified)**: Tang reports **NO per-organ biomass time series** — biomass is
    HARVEST-ONLY (whole-plant ~33.5 g FW + ear ~6.8 g FW, control; `raw_si/tang2026_harvest_biomass.csv`, Fig-1
    digitized) so it can anchor final-mass/HI but **cannot drive `M(t)`**. What Tang DOES constrain is the per-organ
    **TF (S8)/BCF (S7)** → `f_xy`. Canonical extraction: `docs/literature_db/raw_si/tang2026_doseresponse.csv` (all 5
    soil doses 0.1–100 µg/g). NOTE the dose CONDITION: TF declines with dose (toxicity) while the linear model gives one
    dose-independent TF, so fits use the **0.1 µg/g** lowest dose (environmentally closest) as PRIMARY, the across-dose
    mean as sensitivity.
  - `validation/tang2026_fxy_refit.py` (ORYZA-driven nstem_leaf; OVERRIDE-only, `parameters.json` UNCHANGED) re-fits
    `f_xy` to Tang TF: overall log10 RMSE 1.23→0.53 (@0.1). **GenX 0.233→0.017–0.020** (independently confirms the
    documented ~12× over-prediction; ≈ the 0.013 recalibration). **PFOS 0.013→~0.32** (current value far too low) — but
    note this **DISAGREES with the Yamazaki-W2 0.142**: PFOS `f_xy` is **dataset/condition-dependent** (Yamazaki = Andosol
    clean per-congener water, greenhouse, Indica+Japonica; Tang = flooded paddy-soil pot, Nipponbare, 5 doses) → do NOT
    pin PFOS `f_xy` to a single value. **PFOA 0.040→0.064–0.097** (dose-condition dependent). This EXTENDS
    `VALIDATION_TANG2026_NSTEM_KR.md` (ORYZA driver + explicit data-file fit), not a re-derivation.
  - `validation/mass_drivers_plot.py`: diagnostic that `M_k(t)` is a time-varying growth curve and the growth-dilution
    sink `μ=(dM/dt)/M → 0` at maturity (terminal leaf/grain ⇒ no steady state).
- **Doc↔code reproducibility audit (this session)**: verified every file referenced in CLAUDE.md/README resolves to a
  real repo file (only the runtime artifact `pfas_rice_demo.png` is "missing" by design), and corrected stale test
  counts (was "111"/"92 passing" → **142 collected, 138 pass, 4 HYDRUS-skip**). The one real doc-ahead-of-code gap
  (`oryza_growth`) was closed by d1f5339.
- **App integration — Tang 2026 validation tab (this session)**: surfaced the Tang TF work in the
  Streamlit app as a new **"✅ Tang TF (OOS)"** tab (`app.py` tabs[6]; About moved to tabs[7]) backed by the
  UI-agnostic `model_api.tang_tf_validation()` / `tang_observed_tf()` + `plots.fig_tang_tf()`
  (`tests/test_model_api.py`, `tests/test_plots.py`). For PFOA/PFOS/GenX it shows the **dry-weight** per-organ
  TF (model vs Tang vs Tang-refit `f_xy`), with the dose toggle (mean / 0.1 µg/g), an optional ORYZA-biomass
  driver, and the three caveats made explicit in-UI (dw basis; `f_xy` condition-dependence incl. PFOS
  0.14–0.32; grain structurally ~3–8× under). Refit `f_xy` is override-only (`parameters.json` unchanged).
- **Selectable biomass driver + Tissue-dynamics mass graph (this session)**: `model_api.simulate(biomass=…)`
  (via `_biomass_fn`) selects the organ-biomass driver — **`"oryza"`** (the mechanistic ORYZA2000 Level-1 carbon
  balance `oryza_growth`; the more first-principles choice, consistent with the model's mechanistic/HYDRUS-coupled
  philosophy) or **`"growth_rice"`** (ORYZA IR72 partitioning on a logistic; the lightweight reconstruction). The
  **app now leads with ORYZA2000** (sidebar "Biomass driver M(t)" radio, default ORYZA2000) so Tissue-dynamics / map /
  BAF run on the mechanistic biomass unless switched; the **Tissue-dynamics tab plots the per-tissue PFAS *mass*
  (burden) C_k·M_k** (`plots.fig_burden`, µg/hill, EXTENSIVE) under the concentration plot — where the chemical
  actually ends up (organ *biomass* M_k(t) is already in the Soil & drivers tab). ORYZA biomass is ~0.01 s (no app-speed cost; `_simulate` is
  cached). **DEFAULT = ORYZA2000 (changed this session, user request "일단 ORYZA2000이 기본")**: `model_api.simulate`,
  `simulate_nstem_leaf`, `_default_drivers`, `_biomass_fn`, and `tang_tf_validation` now default to **`"oryza"`** (the
  mechanistic ORYZA2000), matching the app. **Honest caveat / provenance**: the per-congener `f_xy_W2fit`/`L_Ph_W2fit`
  and the `reproduce_demo.py` RMSE-0.029 reproduction were tuned on a **placeholder/`growth_rice`** driver, so switching
  the live default shifts BAFs (short-chain straw/grain +40–70%) and the W2 fit no longer reproduces Yamazaki under the
  default — **pass `biomass="growth_rice"` to match the legacy artifacts**. `reproduce_demo.py` (placeholder `_logistic`)
  and `calibration.py` (synthetic-recovery demo) use their own drivers and are UNCHANGED. Tests: `test_model_api.py`
  (biomass selectable; **default == oryza**), `test_plots.py`.
- **ORYZA2000 transport re-fit (this session) — `f_xy_source="oryza"`**: since the default biomass is now ORYZA2000,
  the per-congener transport params were RE-FIT on it (`validation/refit_oryza.py`): (f_xy, L_Ph, kappa_d) fit to Yamazaki
  on the mechanistic ORYZA2000 biomass + measured Q_TP, written to `params/parameters.json` as `f_xy_oryza`/`L_Ph_oryza`/
  `kappa_d_oryza` (+ `params/refit_oryza.csv`; the legacy `*_W2fit` are PRESERVED for `reproduce_demo`). `build_parameters.py`
  re-merges `refit_oryza.csv` so a rebuild keeps them. `model_api`'s new **`f_xy_source="oryza"`** (via `_transport_defaults`)
  applies all three; `simulate(f_xy_source="oryza", biomass="oryza")` reproduces Yamazaki at **log10 RMSE 0.236** (saturated
  per congener -> reproduction not prediction; PFDoDA(C12) is a structural long-chain outlier, params at ceilings yet ~4-6x
  under). The default `f_xy_source` stays `"recommended"` (monotone physical TSCF); `"oryza"` is the opt-in reproduction
  calibration (the ORYZA analog of `"W2fit"`). Test: `test_model_api.py::test_oryza_refit_reproduces`. The constrained
  DOF>0 structural-adequacy result (straw ~0.18; `validation/structural_adequacy_fit.py`) is the meaningful goodness-of-fit;
  this saturated re-fit is the operational calibration on the new default driver.
- **Long-chain (C10-C12) mechanism sci-adk sub-investigation (this session)**: `sci_adk_review/proposal_longchain.md`
  + `build_longchain.py` (→ `runs/pfas-rice-longchain`) + `validation/longchain_mechanism.py` adjudicate WHY long chains
  are under-predicted, on the ORYZA2000 biomass. Verdicts: **LC1 SUPPORTED** (free-anion loading structurally starves
  long-chain shoot — free-only long-chain straw+grain log10 RMSE 2.03 ~100×, and the re-fit hits f_xy=1/L_Ph=1 ceilings
  yet PFDoDA straw 14.6 vs 49.8 → the Cw=C/B free-conc collapse throttles loading); **LC2 SUPPORTED** (the B-independent
  lipid bound-loading term `g_xy·C`/`g_ph·C` cuts long-chain straw+grain RMSE 2.03→0.43 ~5×, whole series 1.04→0.39);
  **LC3 REFUTED** (single-pool cost: long-chain root degrades, PFUnDA 20.6→3.9 / PFDoDA 159→4.4, and PFDoDA shoot still
  ~3-4× under). **Conclusion**: lipid-facilitated bound loading is the correct long-chain *direction* but needs a **2-pool
  (free + lipid-bound) split** (so the bound pool feeds the shoot without draining the root) + a PFDoDA residual mechanism
  (irreversible/hysteretic sorption). In-sample. Guard `test_sci_adk_rigor.py::test_longchain_run_reproduces`.
  **LC4 (2-pool root prototype) — CONTESTED**: `validation/twopool_longchain.py` splits the root into a mobile pool
  (water+protein, low binding; feeds the xylem + soil uptake) and a slow-exchanging lipid/cell-wall bound store (holds the
  measured root burden), so lipid-facilitated loading draws from the mobile pool WITHOUT subtracting the large bound store.
  Result: it **closes the LC3 root tradeoff for mid-long chains** (PFDA C10 matches root AND shoot simultaneously
  3.5/4.2·5.0/3.5·4.1/3.4; PFUnDA C11 root within ~2×) — which the single pool could not — **but FAILS for PFDoDA C12**
  (mobile pool rm=0.02 starves → bound root 1.2 vs 69). The PFDoDA residual is an **uptake (jR) mass-balance limit**, not
  internal distribution → needs a different long-chain uptake / irreversible-sorption mechanism. Recorded as hyp-lc-twopool
  (CONTESTED) in `runs/pfas-rice-longchain`. Prototype only (not wired into the core).
  **LC5 (PFDoDA uptake lever)**: scanning the 2-pool, membrane **conductance kappa_d is REFUTED** (LC5a)
  — ×5000 leaves PFDoDA root ~1 vs 69 because GHK anion exclusion caps the internal free conc at Cwo/e^N
  (e^N≈107) regardless of conductance; the **active carrier Vmax is SUPPORTED** (LC5b) — ×5 (20→100)
  overcomes the exclusion and reaches PFDoDA root 62/69 and grain 46/45.5 (straw 102, ~2× over). So the
  longest-chain residual is an **active-carrier-capacity limit**; the complete long-chain resolution =
  2-pool (free+lipid-bound) + lipid-facilitated loading + enhanced long-chain active-carrier uptake
  (consistent with the literature's active carrier-mediated root uptake). `runs/pfas-rice-longchain` now
  holds LC1–LC5 (6 hypotheses); in-sample/prototype, core unchanged.
  **LC6 (carrier-enhancement QSPR) — REFUTED, via the canonical `sci-adk run` CLI**: a separate run
  `runs/pfas-rice-carrier` compiled from `sci_adk_review/proposal_carrier_qspr.md` with the CLI
  (`sci-adk run` → author evidence/verdict → `sci-adk resolve`/`verify`/`prior-work`, not a programmatic
  builder). Tests whether the long-chain carrier enhancement (LC5b's PFDoDA ~5× Vmax) is a smooth
  function of chain length: per-congener Vmax multiplier reproducing the measured root is PFOA 1.2× ·
  PFNA 1.3× · PFDA 1.2× · PFUnDA 2.0× · PFDoDA 5.5×, and log10(multiplier) regresses on n_C at R²=0.70
  (on log K_PL R²=0.62) — NOT log-linear (<0.9): ~no enhancement to C10 then a steep threshold-like
  onset at C11–C12. So the long-chain carrier enhancement is **NOT cleanly QSPR-able** from chain
  length; it stays a longest-chain-specific (ad-hoc) lever. Guard `test_carrier_run_reproduces`.
  **Literature (genuine sci-adk acquisition + source verification)**: `sci-adk prior-work --searched` ran paperforge +
  Unpaywall (contact email `~/.config/sci-adk/config.toml`) over 7 DOIs that corroborate LC1/LC2; ALL 7 are paywalled
  (no OA PDF) → recorded `acquired 0/failed 7` in `evi-lit-*` + a `prior_work_decision` item + `literature/manifest.csv`
  (DOIs still cited in the draft). **5 of 7 were then obtained out-of-band and READ to verify the corroboration at source**
  (`evi-lc-litread`; paywalled PDFs NOT committed — copyright): Chen2025 ES&T 2025,59,82–91 `10.1021/acs.est.4c06734`
  confirms membrane–water K_MW rises **+0.36/CF₂ monotone C4→C16** while protein **HSA affinity peaks at C6–C10** → the
  lipid (membrane) pool, not protein, carries the longest chains (the B-independent lipid-term basis); `newcontam-0025-0007`
  (long-chain root/soil adsorption vs short-chain shoot mobility) + `acsestengg.4c00107` (MW top predictor of TF) +
  `s40726-020-00168-y` + `acs.est.7b06128` corroborate LC1. NOT obtained (2025): `10.1021/acs.est.5c11716`,
  `10.1139/er-2025-0116`. paperforge is the optional `[tools]` extra; the contact email is required for the polite pool (E4).
- **Out-of-sample cross-dataset prediction test (this session) — REFUTED, via the canonical `sci-adk run` CLI**:
  the central predictive-validation result on data NOT used to fit. The main run's H3 ("Yamazaki = predictive
  validation") was REFUTED but that was Yamazaki-on-itself (saturated vs a-priori). The decisive test is out-of-sample
  prediction on an INDEPENDENT dataset. `sci_adk_review/proposal_oos_tang.md` → `runs/pfas-rice-oos-tang` (compiled via
  the `sci-adk run` CLI, then author evidence/verdict → `sci-adk resolve`/`verify`): the model driven by theory/QSPR
  monotone `f_xy` (`f_xy_source="recommended"`, NOT fit to Tang) predicts Tang 2026's per-organ TF (stalk/leaf/endosperm,
  dw; PFOA/PFOS/GenX, 0.1 µg/g) at **OOS log10 RMSE 1.232** vs **in-sample Tang-refit 0.519** (~5× worse; systematic miss —
  PFSA ~40–200× under, GenX ~10× over). **hyp-001 REFUTED**: the structure can REPRODUCE Tang by fitting (0.52, consistent
  with the structural-adequacy result) but does NOT PREDICT an independent dataset with parameters fit elsewhere — confirming
  H3/H4 at the cross-dataset level. The PFSA/GenX directional miss re-confirms the `f_xy` head-group offset / ether QSPR are
  dataset/condition-dependent (Yamazaki Andosol clean water vs Tang flooded paddy), not pinnable to a single value.
  `validation/oos_tang.py`; guard `test_oos_tang_run_reproduces` (`sci-adk verify` exit 0, digest 46d71f24).
- **Does the lipid mechanism GENERALIZE out-of-sample? (this session) — SUPPORTED, via the `sci-adk run` CLI**:
  the positive follow-through to the OOS REFUTED baseline. The OOS failure above was the *free-anion* model. The
  long-chain investigation's **lipid-facilitated loading** (LC2 SUPPORTED; B-independent `g_xy·C`/`g_ph·C`, K_PL-gated)
  had its `LIPID_LOADING` constants **fit on Yamazaki (excl. PFDoDA), NOT on Tang** (`docs/fxy_longchain_lipid_exploration.md`),
  so turning it on for Tang is a genuine out-of-sample generalization test. `sci_adk_review/proposal_oos_lipid.md` →
  `runs/pfas-rice-oos-lipid` (compiled via `sci-adk run` → evi-oos-lipid SUPPORTS → `resolve`/`verify`): with NO parameter
  touched for Tang, `lipid_loading=True` drops the Tang OOS log10 RMSE from **1.232 (free-anion) → 0.516**, matching the
  in-sample Tang-refit (0.519). The dominant free-anion failure (PFOS, the high-K_PL sulfonate, ~40–200× under) is fixed at
  the mechanism level (stalk 0.013→0.620 vs Tang 0.571), exactly as the K_PL-gated lipid term predicts and as Chen2025
  (membrane K_MW monotone) independently corroborates. **hyp-001 SUPPORTED — the project's first strong cross-dataset
  out-of-sample predictive success**: a mechanism fit on one dataset predicts an independent dataset's per-organ pattern
  (the *mechanism* generalizes, not added fitting). Honest residual: GenX (ether) stays over-predicted (provisional ether
  f_xy offset — a separate condition-dependent issue, not lipid loading) and PFOS endosperm ~5× under. `tang_tf_validation`
  gained a `lipid_loading` arg; `validation/oos_tang_lipid.py`; guard `test_oos_lipid_run_reproduces` (verify exit 0,
  digest 684c31e2). This is a genuine OOS success (not in-sample reproduction), so SUPPORTED is justified — distinct from
  the hyp-yamazaki/grain over-claim guard. EXPLORATORY: lipid loading stays opt-in (default off); the core is unchanged.
- **Multi-dataset OOS robustness (this session) — SUPPORTED, via the `sci-adk run` CLI**: the §8.1 lipid OOS success was
  only 3 Tang congeners. `sci_adk_review/proposal_oos_multidataset.md` → `runs/pfas-rice-oos-multidataset` (CLI `sci-adk
  run` → evi-oos-multidataset SUPPORTS → `resolve`/`verify`) transfers three model variants (monotone/free-anion, saturated
  W2, K_PL-gated lipid) WITHOUT refit to three independent datasets via `validation/oos_multidataset.py` (= `oos_tang_lipid.py`
  + `oos_crossdataset.py`). The lipid mechanism wins decisively on BOTH clean datasets: **Tang 2026** per-organ TF 0.52 vs
  free-anion 1.23, and **Kim 2019** brown-rice grain BAF (excl PFOA) 0.48 vs monotone 2.05 vs W2 1.07 (reliable DF≥15%: 0.20
  vs 1.92 vs 1.44) — and lipid uniquely captures the Kim grain long-chain RISE the baselines structurally miss. So the OOS
  generalization is **NOT a Tang artifact** — it holds across two independent datasets (Korean field grain + Chinese pot
  per-organ). Honest limits (pre-registered): **Li 2025** is field/group-water/surface-confounded and inconclusive (W2 wins
  straw/root TF 0.33 vs lipid 0.57, but lipid wins grain/root 0.72 vs 1.15/1.47), and Kim long chains are low-DF (3–13%).
  hyp-001 SUPPORTED (digest 68ebaf39); guard `test_oos_multidataset_run_reproduces`. Core unchanged; lipid stays opt-in.
- **Leaf senescence-loss flux (this session) — fixes the ORYZA leaf-TF artifact**: with the mechanistic
  ORYZA biomass the leaf shrinks (senescence), so the growth-dilution sink `μ=(dM/dt)/M` goes NEGATIVE and the
  `−μ·C` term spuriously CONCENTRATES the leaf — but `oryza_growth` models that loss as leaf DEATH (carbon removed
  from the plant), so the dead/shed leaf should carry its PFAS away. FIX: `oryza_growth` now exposes the leaf death
  rate `drlv(t)` (`organ_biomass_oryza`/`oryza_drivers` extra key `leaf_death_rate`/`leaf_loss`), and the PFAS leaf
  ODE (`4pool_surf` + `nstem_leaf`, via a new optional `PlantInputs.leaf_loss`) subtracts `−leaf_loss·C` with
  `leaf_loss = drlv` (since `D/M_leaf = drlv` EXACTLY), cancelling the death part of `−μ·C` so only the always-diluting
  growth term remains. **Scoped to the ORYZA path** — `growth_rice` has no senescence and supplies no rate (`leaf_loss`
  defaults to 0), so the default/calibration/`reproduce_demo`/tests are UNCHANGED. Effect: PFOA leaf BAF 4.88 (artifact)
  → 2.51 (≈ growth_rice 2.26); the residual small rise is the REAL continued-xylem-input effect (and nudges the Tang
  leaf TF toward the data: 0.93 growth_rice → 1.31 ORYZA vs Tang 1.66). `tests/test_model_api.py::test_oryza_leaf_senescence_loss`.
  NOTE the assumption it encodes: PFAS leaves with the dead leaf at the leaf concentration (uniform); the alternative
  (immobile PFAS retained in situ as mobile dry matter is remobilised) would keep some rise — unmeasured.
- **Grain formation gate (this session) — DPU-consistent; kills the pre-flowering grain spike**: the
  grain/panicle is physically absent until ~flowering, but the ODE floors `M_grain` (1e-4 kg) to avoid 0/0, so a
  trickle of xylem/phloem loaded a tiny burden into the frozen-floor mass → `C=burden/M` ballooned (PFOA grain
  conc spike 2.09 @ d52, **pre**-flowering) then crashed at fruit-set. This is a deviation from the Trapp/Brunetti
  DPU framework, where the grain is a phloem sink whose import is tied to its growth/existence (no loading of a
  not-yet-formed organ). FIX (`4pool_surf` + `nstem_leaf`): a **formation gate** `γ(t)` on `PlantInputs` ramps 0→1
  as `M_grain` LEAVES its floor (glo→1.5·glo); the grain's xylem/phloem influx is scaled by `γ`, and the
  pre-formation share is **rerouted to the leaf (xylem) / not exported (phloem export → (γ+φ))** so the balance
  still closes (mass-conserving). Result: grain rises **monotonically from 0 at flowering** (no spike), terminal
  accumulation intact. **Scoped/robust**: `γ=1` for the whole of grain filling (loading unchanged → `reproduce_demo`
  log10 RMSE stays **0.029**; grain BAF shifts <~5%) and `γ=1` throughout for a **constant-mass driver** (HYDRUS/CSV
  M, no floor → grain always present). `tests/test_model_api.py::test_grain_formation_gate`. The earlier
  display-mask (PR #20) is now backed by the physics gate. NOTE the wrong first cut used a 2%-of-max threshold that
  gated *filling* too (RMSE 0.029→0.34); keying on "mass left the floor" is the correct criterion.
- **Two-pool root — decoupling the root sink from shoot delivery (BAF "고찰" session)**: addresses the central
  mass-balance wall (`docs/fxy_longchain_lipid_exploration.md`): a single root pool cannot reproduce a HIGH long-chain
  root BAF *and* a non-trivial long-chain SHOOT BAF, because the pool whose burden IS the root BAF is the pool that
  feeds the xylem (lipid loading `g·C` fixes long-chain grain but DRAINS the long-chain root). `validation/
  twopool_root_exploration.py` (standalone 5-state ODE `[root_mobile, root_seq, stem, leaf, grain]`; EXPLORATORY,
  in-sample Yamazaki; canonical core + `parameters.json` UNCHANGED) splits the root into a **mobile** pool (binding
  `B_m`; GHK+carrier uptake; loads xylem with the **monotone physical** `f_xy_recommended` + K_PL-gated lipid term) and
  a **sequestered** pool (irreversible apoplast/cell-wall/plaque sink; a TERMINAL accumulator) whose rate `k_seq(n,
  head_group)` is a **NON-K_PL** chain·head-group descriptor. Motivation: PFOS(C8 PFSA) & PFUnDA(C11 PFCA) have
  IDENTICAL K_PL=31623 and near-identical B_k_root (49.4 vs 49.1) yet root BAF 5.93 vs 19.53 (3.3×) — no K_PL-gated sink
  can separate them. **Results**: (1) the structure ties the best prior global model (7 globals, log10 RMSE **0.257**
  excl PFDoDA vs U-shaped-K_PL-f_xy 0.286) **while keeping the monotone physical f_xy** (the straw U-shape emerges from
  lipid loading + root decoupling, NOT a non-physical fitted f_xy). (2) **Root-matched sufficiency test** (back out
  per-congener `k_seq` so model root == obs root): the shoot stays essentially unchanged (straw 0.255, grain 0.307) —
  proving the structure is SUFFICIENT, you CAN hold high long-chain root AND deliver shoot. (3) The empirical `k_seq`
  **separates PFOS (0.047) from PFUnDA (0.210), 4.5× at identical K_PL** — the non-K_PL signature is real & quantified —
  and is **U-shaped** in chain length (PFBA 0.29→PFNA 0.014→PFDoDA 0.49), which is exactly why the LINEAR global `k_seq`
  collapsed (`ks_b→0`). (4) **U-shaped `k_seq(n)` REALIZED (well-posed follow-up, DONE)**: an asymmetric U with the
  RISING arm in **chain length n (NOT K_PL)** — `k_seq=[0.268·e^(−0.52(n−4))+0.615·e^(1.35(n−12))]·{10^+0.18 if PFSA}` —
  fit to the root-matched empirical values then plugged back into the full ODE gives **all-11 (incl PFDoDA) log10 RMSE
  0.251** (root **0.156**, straw 0.260, grain 0.311) and **realizes the separation**: PFOS(C8) k_seq 0.054 vs PFUnDA(C11)
  0.166 (3.1×) → model root PFOS 6.6/PFUnDA 15.9 (was backwards 16.1/9.5 under the linear fit). Root is essentially solved
  incl PFDoDA (82 vs 69); residual is now the **very-long-chain SHOOT** (PFDoDA straw 10.5 vs 49.8 — the C12 carrier-limit
  floor, a shoot problem `k_seq` cannot fix). (5) **OOS transfer (`validation/twopool_root_oos.py`)**: the Yamazaki-fit
  model is transferred WITHOUT re-fitting to independent data, all 4 models on the SAME demo forcings. **Kim 2019 grain
  (clean, PFHpA→PFDoDA): two-pool excl-PFOA log10 RMSE 0.47 = BEST** (mono 1.49, W2 0.57, lipid 1.12) and CAPTURES the
  long-chain grain RISE (2pool PFUnDA 6.1/PFDoDA 8.7 vs monotone 0.19/0.52; obs ~33/35). Honest limits: absolute long-chain
  grain still under (low-DF Kim tail), Kim PFOA grain 4.43 ≫ Yamazaki 0.46 (between-dataset shift), Li 2025 grain/root TF
  root-surface-confounded (inconclusive, as documented). ⇒ OOS SUPPORTS the structure/mechanism but does NOT warrant
  promoting the fitted `k_seq` into `parameters.json` (single clean OOS set; demo forcings). `parameters.json` UNCHANGED
  (exploration-only). Fitted params cached → `validation/twopool_fitted_params.json`. (6) **Long-chain shoot floor diagnosed
  (`validation/twopool_root_seqrelease.py`)**: the residual after the U-shaped k_seq is the very-long-chain SHOOT (PFDoDA
  straw 10.5 vs 49.8). A slow seq→mobile release `k_rel` (added to the ODE, default 0) **cannot** lift it — sweeping k_rel
  the straw barely moves (10.5→13.4) while PFDoDA root COLLAPSES (82→12). The `g_xy` diagnostic localizes the bottleneck to
  the **xylem-LOADING capacity**: reaching PFDoDA straw~50 needs g_xy ×8 (still only 35) and over-feeds PFDA/PFUnDA 3–4×
  (RMSE 0.251→0.665) — **no smooth/QSPR-able loading term selectively lifts C12**. ⇒ the long-chain shoot floor is a
  STRUCTURAL shoot-loading ceiling + near-MQL outlier (obs PFDoDA straw is a 6× jump over PFUnDA for one CF2 vs 3.5× in
  root), outside any ROOT term (k_seq/k_rel) — independently quantifies PR #21 LC5/LC6. The two-pool root (RMSE 0.251) is at
  the achievable floor; residual is NOT a missing root mechanism. (7) **Robust to MEASURED forcings
  (`validation/twopool_root_measured.py`)**: re-fitting the whole model on `forcing_rice.Q_TP` (peak 0.098, ~4× below the
  demo) + `growth_rice` ORYZA-IR72 biomass (HI 0.53) — the forcings the fxy-doc baselines use — gives in-sample RMSE
  **0.278** (root **0.154**), TIES the fxy-doc U-shaped-K_PL-f_xy (0.286) **while keeping monotone physical f_xy**, and the
  **PFOS/PFUnDA separation HOLDS/sharpens to 4.5×** (k_seq 0.031 vs 0.141). Kim grain OOS now apples-to-apples: two-pool
  excl-PFOA **0.56 = ties lipid (0.55)**, crushes mono (2.04)/W2 (1.11) — but keeps the high long-chain root lipid drains.
  ⇒ the structure / monotone f_xy / non-K_PL U-shaped k_seq / separation / OOS all survive realistic biomass+transpiration;
  NOT a placeholder-forcing artifact. Cached → `validation/twopool_fitted_params_measured.json`. Figure
  `validation/figures/twopool_root_exploration.png`; full record `docs/twopool_root_exploration.md`. Still mechanism discovery,
  NOT validation (Yamazaki in-sample fit → OOS transfer; decisive test = per-congener xylem-sap/root-water ratio +
  desorption-resistant root-fraction assay). **Next-session handoff: `docs/HANDOFF_BAF_twopool.md`** (status, open items —
  promotion decision / Tang OOS / opt-in model_api wiring — and a resume prompt).
- **Two-pool wired as a model_api OPT-IN module (this session; handoff item ①)** — `model_api.simulate_twopool_seq(...)`:
  the exploratory two-pool root model is now callable through the UI-agnostic API exactly like `simulate_nstem_leaf`,
  so the app/other validation can use it **without changing any default** (`simulate`/`reproduce_demo`/`parameters.json`
  UNCHANGED). It loads the cached Yamazaki fit (`validation/twopool_fitted_params.json` via the validation module's
  `load_fit()`/`kseq_ushape`/`lipid_g`) and re-implements the 5-state ODE inside `model_api` (driven by the standard
  forcing/`drivers=` machinery) so it returns the **same dict shape** as `simulate()` (root/stem/leaf/grain conc & BAF
  series + finals/`straw_baf`/`tf_final`), plus the root **mobile/seq split** (`conc["root_mobile"|"root_seq"]`,
  `seq_fraction`) and two-pool levers (`k_rel` seq→mobile desorption, `kseq_override`). The reported `root` BAF = mobile
  + sequestered. Defaults (`measured_forcing=False, season=120`) reproduce the documented headline **overall log10 RMSE
  0.251 (root 0.156)** with the **monotone physical `f_xy_recommended`** and the non-K_PL **PFOS/PFUnDA k_seq 3.1×
  separation** at identical K_PL. A drift guard (`tests/test_model_api.py::test_simulate_twopool_seq_matches_validation_and_rmse`)
  pins the wrapper to the standalone validation endpoints (cross-impl RMSE 0.014) so the two implementations cannot
  silently diverge; `test_simulate_twopool_seq_structure_and_keys` / `..._krel_drains_root_to_shoot` lock the I/O contract and
  the Result-5 k_rel behaviour. Still EXPLORATORY / in-sample (the cached fit is on the demo forcings; the measured-forcing
  fit `twopool_fitted_params_measured.json` is not auto-loaded). The §4 promotion decision (handoff item ③) is unchanged —
  **NOT promoted to `parameters.json`** pending the user. **NAMING**: there are now TWO opt-in two-pool root models —
  this SEQUESTRATION one (`simulate_twopool_seq`; irreversible non-K_PL `k_seq` sink, keeps monotone f_xy) and the
  CARRIER one (`simulate_twopool_carrier` / `close_longchain_2pool`, `src/pfas_rice_two_pool.py`; reversible bound store
  tuned by carrier/f_xy levers, the saturated long-chain closure). Different mechanisms; the `_seq`/`_carrier` suffix
  disambiguates (renamed from the colliding `simulate_twopool`/`simulate_two_pool`).
- **App surfacing — two-pool (seq) overlay on the BAF tab (this session)**: the Streamlit **📊 BAF vs observed** tab
  now optionally overlays `simulate_twopool_seq` next to the 4-pool core and the Yamazaki bars (curated congeners
  only; checkbox, EXPLORATORY caveat in-UI) via `plots.fig_baf(res, obs, extra=…)` + the cached `_simulate_twopool_seq`
  in `app.py`. Run at the two-pool's calibrated point (Cwᵒ=1, season≈120) so it is comparable to the fixed observed
  bars and does NOT track the sidebar (the core bar does). The carrier two-pool stays API-only (saturated DOF-0
  closure, ~1 min/congener — too slow to render live). `tests/test_plots.py::test_fig_baf_extra_overlay`. Defaults /
  canonical core / `parameters.json` unchanged.
- **Two-pool seq → Tang 2026 per-organ OOS (this session; handoff item ②) — NEGATIVE/diagnostic**:
  `validation/twopool_root_oos_tang.py` transfers the Yamazaki-fit two-pool (no re-fit) to Tang per-organ TF
  (stalk/leaf/endosperm, dw, 0.1 µg/g). **Result: two-pool OOS log10 RMSE 1.40 — WORSE** than single-pool monotone
  (1.23) and far worse than lipid (0.52, the documented Tang winner). **Why it's informative, not a root-mechanism
  failure**: Tang per-organ is a **SHOOT** test, but the two-pool innovates in the **ROOT** (mobile/seq split); its
  shoot is the unmodified basic 4pool with a **pass-through stem** (PFOA stem 0.008 vs leaf 1.14) → the **stalk TF
  collapses** (the empty-stem defect `nstem_leaf` fixes). The per-organ breakdown isolates it: the two-pool **leaf**
  RMSE 0.38 is the **best of all three models**; only the stalk drags the overall up (the single-pool baselines use
  `nstem_leaf`, so their stalk is populated — an apples-to-oranges SHOOT difference). Tang's congeners are C5–C8, so
  the long-chain root decoupling — the two-pool's whole point — is not even exercised. **Conclusion: Tang is NOT a
  fair OOS of the two-pool root**; a per-organ Tang test needs the two-pool root merged with the `nstem_leaf`
  redistributed shoot (future structural merge). Kim 2019 grain (`twopool_root_oos.py`) stays the informative
  two-pool OOS. Added `simulate_organs(c,p,…)` to `twopool_root_exploration.py` (per-organ stem/leaf split on the
  SAME solve path as `simulate` — root/grain byte-identical). Guard
  `tests/test_model_api.py::test_twopool_simulate_organs_and_tang_passthrough_diagnosis`. Full record:
  `docs/twopool_root_exploration.md` §Result 7. EXPLORATORY; `parameters.json` UNCHANGED (no support for promotion).
- **Time-varying pore-water exposure `cwo_profile` + HYDRUS provisioning + Bayesian identifiability (this session)**:
  the default `simulate(Cwo=…)` holds the pore water CONSTANT (conc==BAF, the BAF-reproduction convention), but a real
  paddy `C_w^o(t)` is time-varying. **`simulate(cwo_profile=…)`** (default `"constant"`, UNCHANGED) makes the time-shape a
  first-class, congener-resolved option: **`"flooded"`** = an analytic Freundlich dilution+leaching shape
  (`pore_water_from_inventory` / `soil_paddy_redox_corrected`; per-congener `K_F = Koc·f_oc`, so short chains LEACH to a
  steep decline and long chains stay BUFFERED — **no HYDRUS engine needed**), **`"hydrus"`** = the real-engine shape. All
  shapes are season-mean-normalised to `Cwo` (the `inputs_from_hydrus` convention) so `Cwo` stays the AVERAGE exposure;
  `cwo_kw` tunes `k_leach` etc. `model_api.cwo_profile_series` is the UI-agnostic builder. **Validated vs the engine**
  (`validation/cwo_profile_check.py`): the analytic `"flooded"` reproduces the HYDRUS DIRECTION (PFBA decline ratio 0.10
  vs HYDRUS 0.08, PFOA 0.64 vs 0.63, corr 0.91–0.95; PFOS/PFDoDA flat in both) with a single `k_leach` knob. **App**: the
  parametric data source has a "Pore-water Cwᵒ(t) shape" toggle + live preview (`plots.fig_cwo_profile`). Tests:
  `test_model_api.py` (constant==default, flooded shape/leach, hydrus-direction guard), `test_plots.py` (`fig_cwo_profile`).
  - **HYDRUS now buildable offline**: the FORTRAN source was **vendored** under `external/hydrus_source/` (de-submoduled —
    the upstream submodule is blocked behind restrictive network policies, binary not in git), and a **SessionStart hook**
    (`.claude/hooks/session-start.sh` + `.claude/settings.json`) auto-installs the Python stack + builds the engine
    (best-effort/non-blocking) on Claude Code on the web. `packages.txt` (gfortran/make) already covers Streamlit Cloud.
  - **Bayesian inverse / identifiability** (`validation/bayesian_inverse_demo.py`): answers "can we infer `Q_TP(t)` &
    `C_w^o(t)` from tissue `C(t)`+`M(t)`?" — YES for the EXPOSURE (`qtp_scale`, `cwo_level`) with transport fixed (Laplace
    posterior from the Fisher Jacobian at truth, cond ~90, recovers), but `Q_TP·f_xy` is a **product ridge** (corr ~−1, cond
    ~500: only the product is constrained — multi-compartment data only PARTIALLY break it, since `Q_TP` also sets intra-shoot
    advection rates `f_xy` does not) and `Cwo` vs root-uptake conductance is even more degenerate (cond ~1e5, no clean product
    invariant — nonlinear GHK+carrier uptake). **Conclusion: pinning `Q_TP`/`Cwo` absolutely needs an independent measurement
    (xylem sap / pore-water probe)**, exactly as §8 notes. `tests/test_bayesian_inverse.py`. Default/`parameters.json` UNCHANGED.
- **Driver-builder biomass bug fix (this session)**: when ORYZA2000 became the `simulate` default, the **driver helpers
  were not updated** — `measured_forcing`, `drivers_from_arrays`, `load_driver_csv`, and `soil_hydrus.inputs_from_hydrus`
  all **hardcoded `growth_rice`** for `M(t)` (their docstrings even claimed "ORYZA"). So the app's **Soil-inventory, CSV-
  driver, and live-HYDRUS** modes (which build `drivers=` and omit `M`) silently ran on growth_rice regardless of the
  sidebar biomass radio — visible as a non-senescing leaf in the Soil & drivers `M(t)` panel. FIX: those helpers now take a
  `biomass="oryza"` arg and build `M` via `_biomass_fn` (and attach the ORYZA `leaf_death_rate` as `drivers["leaf_loss"]` so
  the leaf-senescence correction still applies); `app.py` threads the selected `biomass` into every driver builder
  (`drivers_from_arrays`/`load_driver_csv`/`hydrus_drivers` calls + the cached `_hydrus_drivers_cached`). Now all five
  exposure modes honour the biomass radio (ORYZA2000 by default). The main `simulate`/`_default_drivers` path was already
  correct; `reproduce_demo`/`calibration` use their own drivers and are UNCHANGED.
  `tests/test_model_api.py::test_drivers_from_arrays_respects_biomass_selection`.
- **flooded `k_leach` calibrated to HYDRUS per congener + emcee MCMC cross-check (this session)**:
  - **Per-congener `k_leach` default**: the analytic `cwo_profile="flooded"` shape had a single flat knob (`k_leach=0.02`)
    that under-leached the short chains. `validation/cwo_kleach_calibration.py` now runs the **real HYDRUS-1D engine** for
    all 13 curated congeners, reads each pore-water decline ratio, and fits the `k_leach` that makes the analytic shape
    match → `params/cwo_kleach.csv`. The pattern is **non-monotone** (peaks at PFOA `k_leach`≈0.05, short chains ≈0.025–0.05,
    long chains → 0 since they are buffered), so a per-congener TABLE (not a clean QSPR) is the default; novels/SMILES fall
    back to a `k_leach(log10 Koc)` linear fit (RMSE 0.013). `model_api.default_k_leach(congener|n_C,group)` resolves it;
    `cwo_profile_series(k_leach=None)` (the new default) auto-applies it (explicit `k_leach` still overrides), so PFBA's
    flooded decline now matches HYDRUS (0.072 vs 0.08; was 0.17 at flat 0.02). The app's `k_leach` slider pre-fills the
    calibrated value per congener. `tests/test_model_api.py::test_default_k_leach_is_hydrus_calibrated`. `parameters.json`
    UNCHANGED (the table is a separate artifact loaded directly).
  - **emcee full-MCMC cross-check**: `validation/bayesian_inverse_demo.py` gained `emcee_posterior()` — an affine-invariant
    ensemble sampler that confirms the Laplace verdicts with a sampled posterior (well-posed `(qtp_scale,cwo_level)` recovers;
    `Q_TP·f_xy` is a ridge). It is **OPT-IN** (`python validation/bayesian_inverse_demo.py --emcee`) because the forward ODE
    is ~0.7 s/sample (a chain is minutes); the default run stays Laplace-only. `emcee` is in the new
    `requirements-validation.txt` (optional); `tests/test_bayesian_inverse.py::test_emcee_posterior_recovers_well_posed`
    skips when emcee is absent.
- **General-audience app rework — Simple/Expert split (this session; `docs/HANDOFF_app_general_audience.md` DONE)**:
  `app.py` is re-pitched for a non-expert audience (policy/undergrad/public) via **progressive disclosure** — a sidebar
  **`st.toggle("Expert / advanced controls")`** (default OFF = Simple). **Simple mode** exposes only a friendly congener
  dropdown (`_FRIENDLY_CONG` names) + a low/medium/high **contamination preset** (→ Cwᵒ), a plain-language headline (3
  metric cards + one summary sentence, NO BAF/Cwᵒ/f_xy/eᴺ symbols), and 4 jargon-free tabs (🗺️ Where it goes / 📈 Build-up
  over time / 📊 How much builds up / ℹ️ About & glossary). **Expert mode** restores 100% of the prior UI (5 data-source
  modes, SMILES, E_m/f_xy/biomass, 8 tabs). Added across the page: an **intro card**, a prominent **research/educational
  disclaimer** banner (top + footer, every screen), and a **footer** (version/repo/docs/cite). Plain-language **glossary**
  (`_glossary_md`) in the About tab + Simple. **Biomonitoring fix**: the Tissue-dynamics / Soil & drivers tabs now carry a
  "model reference, not your measured data" warning. **CSV/PNG export** (download buttons): pure helpers
  `model_api.summary_csv()` / `timeseries_csv()` (BAF table + driver/tissue series) + a graceful kaleido PNG of the plant
  map (`_png_bytes`; degrades to a caption when kaleido absent). Two plain-language plot builders
  `plots.fig_buildup_plain` / `fig_where_plain` (friendly tissue names, no symbols). Consistency: the soil-inventory
  `k_leach` slider now uses the per-congener `api.default_k_leach` default (0–0.15), matching the parametric flooded mode.
  **UI-only — `parameters.json`, the model math, `reproduce_demo` (RMSE 0.029), and `simulate()` are UNCHANGED**;
  `model_api` gained only pure string export helpers. Tests: `test_model_api.py::test_export_csv_helpers`,
  `test_plots.py::test_plain_language_figures_build`; full suite **174 passed, 2 skipped**. Verified with headless
  Streamlit + Playwright screenshots of both the Simple landing and the Expert UI.
- **Bayesian inverse exposure estimate wired into the app (this session)**: a user-facing **Bayesian
  parameter estimation** added to BOTH the Simple ("🔎 Work backwards") and Expert ("🔎 Inverse (Bayesian)")
  tabs. `model_api.estimate_exposure_bayesian(congener, measured_conc, sigma_log10=…)` infers the pore-water
  contamination level **Cwᵒ from measured tissue concentrations** (root/straw/grain, any subset) WITH a
  credible interval — the inverse of the forward question. Because root uptake is **saturable** (GHK +
  carrier), tissue conc is a NONLINEAR monotone function of Cwᵒ, so this is a real inverse (not a division):
  it finds the MAP exposure by a **quadratic-fit Laplace** in log10(Cwᵒ) (a coarse then local parabola → MAP +
  curvature = posterior width; ~8 ODE solves, deterministic) and returns median + 68/95% CI + a plotting grid
  + the model's tissue fit at the MAP. Same Laplace idea as `validation/bayesian_inverse_demo.py`, which
  established this is the **well-posed** direction (the EXPOSURE level is identifiable with transport fixed;
  Q_TP·f_xy and Cwᵒ-vs-conductance are ridges — surfaced as the Expert tab's caveat). Plot builder
  `plots.fig_exposure_posterior` (log-x posterior + 95% band + median). App: a cached `_estimate_exposure`
  + a shared `_render_inverse_estimator` (gated behind an Estimate button so the ~8 solves don't run on every
  rerun). Synthetic recovery verified (median recovers the known Cwᵒ within a few %, truth inside the 95% CI).
  Tests: `test_model_api.py::test_estimate_exposure_bayesian_recovers_and_brackets`,
  `test_plots.py::test_fig_exposure_posterior_builds`. UI/inverse only — `parameters.json` and the model math
  are UNCHANGED.
- **Editable data tables (growth curve + Cwᵒ(t)) in the app (this session)**: users can now type/paste their
  own **growth table** (organ FRESH-weight mass over time) and **time-varying pore-water Cwᵒ(t)** (absolute
  µg/L) as editable grids (`st.data_editor`) + CSV upload, in BOTH the Simple ("📋 Use my own data tables"
  checkbox) and Expert ("Custom tables (Cwᵒ + growth)" data-source mode). `model_api.drivers_from_tables(growth,
  cwo, growth_units=…, Cwo_const=…)` sorts+interpolates the rows onto the model grid and builds the standard
  `simulate(drivers=…)` dict (it reuses `measured_biomass.to_kg_per_hill` for the growth units g/hill·kg/hill·
  g/m2·kg/ha·t/ha and `drivers_from_arrays` for the partial-input fallback — a Cwᵒ table alone still runs on the
  biomass driver's M, a growth table alone on a flat `Cwo_const`). **Per-compartment density** `ρ_k` [kg/L fresh]
  is exposed (default `DEFAULT_TISSUE_DENSITY` root1.0/stem0.30/leaf0.30/grain1.20, editable): the growth table
  is FRESH-WEIGHT MASS (the model's M unit), and ρ_k is the **mass↔volume bridge** — used to report the implied
  organ volume (M/ρ) for consistency. The transport ODE is **mass-based** (there is NO density prefactor — the
  early-draft ρ term was dimensionally wrong, see §8/the module headers), so ρ_k does not alter the integration;
  it is a stated per-compartment property + per-volume reporting, exactly as requested. Test:
  `test_model_api.py::test_drivers_from_tables_growth_and_cwo`. UI/driver only — `parameters.json` and the model
  math are UNCHANGED.
- **Korean for the general-audience (Simple) mode (this session)**: the non-expert view now renders in **Korean**,
  while Expert stays English (the sidebar toggle is bilingual `🔬 전문가/고급 모드 (Expert / advanced)`). All Simple-gated
  `app.py` text is Korean (sidebar, intro, disclaimer `_DISCLAIMER_KO`, headline metrics + summary, the 5 tab names +
  captions, glossary `_glossary_md(ko=True)`, friendly congener names `_FRIENDLY_CONG_KO`/`_cong_label_ko`, contamination
  presets `_PRESETS_KO`, the inverse estimator `_render_inverse_estimator(simple=True)`, the custom-tables panel
  `_render_custom_tables(ko=…)`, downloads, footer). The Plotly builders gained a `lang` arg (default `"en"` so Expert +
  tests are unaffected): `fig_buildup_plain`/`fig_where_plain`/`fig_exposure_posterior` and the plant map
  `fig_plant_schematic`/`fig_schematic_from_res`/`fig_schematic_animated` render Korean titles/axes/organ labels
  (뿌리/줄기/잎/낟알/짚) when `lang="ko"`. UI/i18n only — `parameters.json` and the model math are UNCHANGED. Tests:
  `test_plots.py::test_plain_figures_korean_variant` (English defaults still asserted); verified with headless Streamlit
  + Playwright screenshots of the Korean Simple landing and the English Expert UI.

## 7. Build & run
- `pip install -r requirements.txt`
- **Main reproduction**: `python reproduce_demo.py` (Yamazaki BAF, W2 fit, RMSE≈0.029);
  `--rec` uses the monotone f_xy. Rebuild params: `python build_parameters.py`.
- **Visualization tool**: `pip install -r requirements-app.txt && streamlit run app.py`
  (plant/soil accumulation colormap + HYDRUS/soil/biomonitoring modes; see `docs/visualization_tool.md`).
- **Live HYDRUS-1D** (optional, for the "Run HYDRUS-1D (live)" mode): the FORTRAN source is now
  **VENDORED** under `external/hydrus_source/` (de-submoduled — the upstream `phydrus/source_code`
  submodule is unreachable behind restrictive network policies, and the compiled binary is not in
  git), so no submodule init is needed — just build + install phydrus:
  `cp external/hydrus_source/makefile external/hydrus_source/source/ &&
  (cd external/hydrus_source/source && make)` (gfortran); `pip install phydrus`. Demo: `python src/soil_hydrus.py`.
  On **Claude Code on the web** the **SessionStart hook** (`.claude/hooks/session-start.sh`) does all
  of this automatically (installs the Python stack + builds the engine, best-effort/non-blocking).
- Plant demo: `python src/pfas_rice_plant_module_4pool_surf.py` (N, B_k, BAFs; saves `pfas_rice_demo.png`).
- Multi-height stem: `python validation/nstem_gradient_check.py` (stem-gradient direction vs Yamazaki).
- Mechanistic ORYZA biomass: `python src/oryza_growth.py` (IR72 potential sanity);
  `python validation/oryza_growth_validation.py` (vs `growth_rice` + BAF driver-sensitivity + figure).
- Measured-biomass driver: `python src/measured_biomass.py` (template → M(t) drivers demo).
- Mass drivers: `python validation/mass_drivers_plot.py` (M_k(t), dM/dt, growth-dilution μ figure).
- Two-pool root: `python validation/twopool_root_exploration.py` (root sink ↔ shoot decoupling; global fit +
  root-matched sufficiency test + non-K_PL U-shaped k_seq fit; ~3 min, saves `figures/twopool_root_exploration.png` +
  `twopool_fitted_params.json`). OOS transfer: `python validation/twopool_root_oos.py` (Yamazaki-fit → Kim 2019 grain +
  Li 2025 TF, no re-fit; reuses the cached fit, ~5 s). Tang per-organ OOS (NEGATIVE/diagnostic):
  `python validation/twopool_root_oos_tang.py` (Yamazaki-fit → Tang stalk/leaf/endosperm dw TF, no re-fit; ~25 s;
  two-pool 1.40 worse than lipid 0.52 — pass-through stem collapses the stalk; Tang tests the shoot, two-pool fixes
  the root). Long-chain shoot-floor diagnostic:
  `python validation/twopool_root_seqrelease.py` (k_rel seq-release sweep + g_xy xylem-loading diagnostic; ~20 s).
  Measured-forcing robustness re-fit: `python validation/twopool_root_measured.py` (re-fits on forcing_rice + ORYZA
  biomass; in-sample + Kim OOS vs fxy-doc baselines; ~3 min). Opt-in API (no re-fit; reuses the cached fit):
  `model_api.simulate_twopool_seq("PFUnDA")` → the standard `simulate()` dict + root mobile/seq split.
- Tang 2026 f_xy: `python validation/tang2026_fxy_TF_validation.py` (4-pool TF vs Tang, ORYZA-driven);
  `python validation/tang2026_fxy_refit.py` (nstem_leaf + ORYZA f_xy re-calibration; 0.1 µg/g dose primary).
- **Time-varying exposure `cwo_profile`**: `simulate(cwo_profile="flooded")` gives an analytic
  Freundlich dilution+leaching `C_w^o(t)` (short chains leach, long chains buffered; engine-free),
  `"hydrus"` the real-engine shape, `"constant"` the flat default (conc==BAF). Both shapes are
  season-mean-normalised to `Cwo`. The flooded `k_leach` defaults PER CONGENER (calibrated to HYDRUS,
  `params/cwo_kleach.csv`; `model_api.default_k_leach`). Validate the shape vs the engine:
  `python validation/cwo_profile_check.py` (analytic vs HYDRUS direction; saves `figures/cwo_profile_check.png`).
  (Re)calibrate the per-congener `k_leach` table: `python validation/cwo_kleach_calibration.py` (runs the
  engine for all 13 congeners → `params/cwo_kleach.csv`). In the app: the "Model (parametric)" data source
  has a "Pore-water Cwᵒ(t) shape" toggle + live preview (`plots.fig_cwo_profile`), `k_leach` slider pre-filled
  with the calibrated value.
- **Bayesian inverse / identifiability**: `python validation/bayesian_inverse_demo.py` — infers the
  EXPOSURE (`qtp_scale`, `cwo_level`) from tissue `C(t)`+`M(t)`, and a Laplace posterior from the
  Fisher Jacobian at truth shows the ridges: `(qtp_scale, cwo_level)` is identifiable with transport
  fixed (cond ~90), but `Q_TP·f_xy` is a product ridge (corr ~−1, cond ~500) and `Cwo` vs root-uptake
  conductance is even more degenerate (cond ~1e5, no clean product invariant — nonlinear uptake). So
  pinning `Q_TP`/`Cwo` absolutely needs an independent measurement (xylem sap / pore-water probe). Add
  `--emcee` for a full-MCMC cross-check (opt-in; needs `pip install -r requirements-validation.txt`, a
  few minutes — the ODE is ~0.7 s/sample).
- Soil → plant (analytic): `python src/soil_paddy.py` (legacy) / use `soil_paddy_redox_corrected` for redox.
- **Soil → plant (REAL HYDRUS-1D)**: the source is vendored, so build the engine once, then run:
  ```
  cp external/hydrus_source/makefile external/hydrus_source/source/
  (cd external/hydrus_source/source && make)          # needs gfortran
  pip install phydrus
  python src/soil_hydrus.py                            # per-congener pore-water summary
  python validation/hydrus_coupled_run.py             # full soil→plant + figure/CSV
  ```
  (On Claude Code on the web the SessionStart hook builds this automatically.)
- Calibration: `python src/calibration.py`; Literature params: `python src/literature_params.py`.
- **Structure (SMILES) input**: `pip install -r requirements-structure.txt` (RDKit), then
  `python src/pfas_structure.py` (SMILES → descriptors → Compound demo). In code:
  `model_api.simulate_from_smiles("OC(=O)C(F)(F)...")` runs the ODE for any PFAS structure.
- Tests: `pip install pytest && pytest` (173 collected, all pass with the full stack — RDKit + the built
  HYDRUS-1D engine + phydrus, as the SessionStart hook provides on the web; the `test_sci_adk_rigor.py`
  module additionally skips unless `sci-adk` is installed, which CI's `rigor.yml` provides). On a bare
  clone the structure/SMILES tests skip without RDKit and the HYDRUS-engine tests in `test_soil_hydrus.py`
  / the `cwo_profile='hydrus'` guards skip when the engine is unbuilt.
- FORTRAN (Method B): init submodule (`git submodule update --init`), then follow
  https://phydrus.readthedocs.io/en/latest/getting_started/compilation.html
  (gfortran + `makefile` / `make.bat`). NOTE: the top-level `makefile` lists the `.FOR` files
  without a path, so build from inside `source/` (copy the makefile in, as above).

## 8. Conventions
- Units: time **day**; aqueous conc **µg/L**; tissue conc **µg/kg**; mass **kg**;
  flow **L/day**; `B_k` in **L/kg fw** (`C_k = B_k · C_w,k`).
- **Exposure `C_w^o(t)`**: the default scenario holds it CONSTANT (`Cwo`, so conc==BAF). `Q_TP(t)` is
  ALWAYS time-varying (FAO-56 `forcing_rice`). For a time-varying exposure use `simulate(cwo_profile=
  "flooded"|"hydrus")` (mean-normalised to `Cwo`) or supply a `drivers=` series. Architecturally both
  `Cwo` and `Qtp` are time functions (`PlantInputs` interpolants); only the scenario default fills `Cwo` flat.
- **Binding = basis A (fresh weight)**: `B_k = θ_fw + (1−θ_fw)·(f_prot·K_prot + f_PL·K_PL + f_cw·K_cw)`.
  `θ_fw` = fresh-weight water fraction; `f_*` = **dry-weight** mass fractions; `K_*` in L/kg pool-dw.
  The `(1−θ_fw)` factor is a **dry→fresh conversion** (mandatory; the legacy naive `θ+Σf·K` over-states
  B_k ~3×) — it is NOT the old dimensionally-wrong `ρ_k` density prefactor (still absent). Compare to
  dw-reported data via `C_dw = C_fw/(1−θ_fw)`. `f_cw` = whole cell wall (poly+lignin), K = `K_cw_wholecw`.
- `f_xy` ∈ (0,1] is the root→xylem loading factor (TSCF analog): only `f_xy·C_1/B_1` enters the
  ascending xylem (`f_xy=1` = unrestricted DPU). NOTE it does **not** yield a universal
  `root>straw>grain` — the ordering is **congener-dependent** (short: straw>root; long: root>straw).
  **REVISED (`docs/fxy_longchain_lipid_exploration.md`)**: the data require a **non-monotone (U-shaped)**
  effective `f_xy`, not the monotone `f_xy_recommended` — the long-chain rise is REAL (lipid-facilitated
  translocation driven by measured `K_PL`), not the "non-physical W2 artifact" the older framing claimed.
- **Lipid-bound loading (opt-in, default off)**: `Compound.g_xy`/`g_ph` add a B-independent
  `g·C` term to xylem/phloem loading (free anion is `f_xy·Cw`, but `Cw=C/B` starves high-binding long
  chains; the bound pool rides the lipid phase). `model_api.simulate(lipid_loading=True)` uses the
  `K_PL`-gated fit; cuts monotone error 0.98→~0.36 and fixes long-chain grain, but trades off root
  (single-pool limit). EXPLORATORY / in-sample.
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
  is the **base soil engine only** (no DPU/PFAS/ionizable modules). It is now **VENDORED** under
  `external/hydrus_source/` (de-submoduled: the upstream submodule clone is blocked behind restrictive
  network policies, e.g. Claude Code on the web, and the compiled binary is a build artifact, not in
  git). LGPL-3.0 `LICENSE` retained; the 2.7 MB manual PDF and build artifacts (`hydrus`, `*.o/*.mod`)
  are gitignored. Build it with `make` (gfortran); the SessionStart hook does this on the web.
- Key references: Brunetti 2019 *WRR* `10.1029/2019WR025432`; 2021 *ES&T*
  `10.1021/acs.est.0c07420`; 2022 *J. Hazard. Mater.* `10.1016/j.jhazmat.2021.127008`.
