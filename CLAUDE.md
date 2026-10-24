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
│   ├── neutral_dpu.py                        # NEUTRAL-organic (Briggs/Kow) path: same ODE with z=0 (no exclusion/carrier)
│   ├── plant_air.py                          # plant-AIR exchange: volatilisation + gaseous uptake (opt-in; 0 at K_AW=0)
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
│   ├── oryza_growth.py                       # REDUCED re-coding of the ORYZA2000 Level-1 carbon balance → M_s(t) (the biomass="oryza" DEFAULT)
│   ├── oryza2000.py                          # ORYZA2000 potential production ported line by line from the manual, IR72.D92 verbatim (biomass="oryza2000", opt-in)
│   ├── measured_biomass.py                   # ingest a MEASURED per-organ biomass table → M(t) driver (units→kg/hill; Tang etc.)
│   └── organ_mass.py                         # crop-model DRY matter → the FRESH organ mass the ODE runs on (mass_basis)
├── examples/                         # ready-to-load CSVs for app.py (HYDRUS drivers + biomonitoring)
├── params/                           # parameters.json (CANONICAL) + source CSVs (Bk, f_xy, Kcw, ...)
├── data_obs/                         # observed BAF/TF (Yamazaki, Li2025) + yamazaki_stem_height.csv
│                                     #  + NEUTRAL: liu2023/ge2017 (rice) · li2019_rcf (48 hydroponic RCF, 11 spp,
│                                     #    `subset` col holds Briggs' own rows out) · tscf_obs_schriever2020 (97 TSCF)
├── validation/                       # S6 + nstem + hydrus_coupled_run reproduction scripts + figures
├── docs/
│   ├── OVERVIEW_KR.md                # ★ 종합 진입점: 기능·검증·데이터공백·필요실험·notation 표 (+모식도)
│   ├── pfas_rice_compartmental_model.tex / dpu_model_summary_corrected.tex
│   ├── DELIVERABLE_GAP_A_Kcw.md / DELIVERABLE_GAP_B_fxy.md / theory_anchor.tex / H8_handoff_S6_final.md / sources.csv
│   ├── neutral_dpu_validation.md     # NEUTRAL path: anchors, a-priori results (Liu/Ge/Hwang), open gaps
│   ├── HANDOFF_neutral_next.md       # ★ NEXT SESSION: A1 air exchange -> A3 Hwang -> A2 API
│   ├── visualization_tool.md         # app.py guide: plant/soil map, 4 modes, HYDRUS I/O, biomonitoring
│   └── literature_db/                # curated parameter DB (.xlsx + per-sheet .csv) + raw_si/ SI extractions
├── external/hydrus_source/           # VENDORED HYDRUS-1D 4.08 source (de-submoduled from phydrus/source_code; binary gitignored)
├── .claude/                          # SessionStart hook (hooks/session-start.sh): web deps + HYDRUS engine build
├── data/                             # (gitignored)
└── tests/                            # pytest (434 collected): plant, soil, hydrus, calibration, lit params, API (+two-pool, cwo_profile, k_leach), plots, structure(SMILES), oryza, measured-biomass, organ-mass basis (+fresh-basis re-fits), bayesian-inverse, NEUTRAL-organic (Briggs/Kow), twopool-nstem merge

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
  exposes it; the demo fits `L_Ph` to the PFOA grain BAF (now the per-site estimate 9.94 L/kg, 1.27 → 9.94 at `L_Ph`
  0.067; to the old ratio-of-means 4.43 it took 0.022 — the "0.07 → 4.43, `L_Ph≈0.84`" once recorded here was stale). The
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
    note this **DISAGREES with the Yamazaki-W2 0.142**: PFOS `f_xy` is **dataset/condition-dependent** (Yamazaki = an
    Andosol paddy FIELD under ambient irrigation-water exposure, Kasalath+Koshihikari — NOT a greenhouse/clean-water/spiked
    study, corrected by the sampling-protocol audit below; Tang = flooded paddy-soil pot, Nipponbare, 5 doses) → do NOT
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
  dataset/condition-dependent (Yamazaki Andosol paddy field vs Tang flooded paddy pot — both flooded; "clean water" was a
  mis-description, see the sampling-protocol audit), not pinnable to a single value.
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
  **Li re-scored stage by stage later (handoff step 3)** — no longer inconclusive; see the "Li 2025 stage by stage" bullet.
  **NOTE the "reliable DF≥15%" Kim score is ONE congener (PFNA)** — the script's label said PFHpA+PFNA, but PFHpA's
  rice DF is 13%; see the "Kim 2019 with its non-detects" bullet.
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
  grain still under (low-DF Kim tail), Kim PFOA grain 4.43 (per site ~9.9, step 4b) ≫ Yamazaki 0.46 (between-dataset shift), Li 2025 grain/root TF
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
  fit `twopool_fitted_params_measured.json` is not auto-loaded). The §4 promotion decision (handoff item ③) is now
  **DECIDED with the user — DO NOT promote; keep opt-in** (`parameters.json`/defaults/`reproduce_demo` RMSE 0.029 UNCHANGED):
  the Tang per-organ OOS added no support (negative/diagnostic) and the `k_seq` mechanism review
  (`docs/twopool_kseq_mechanism.md`, PARTIALLY SUPPORTED) strengthens the story but provides no direct measured
  root-sorbent dataset; the promotion gate is the §5 rice-root cell-wall/Fe-Mn-plaque batch-sorption + desorption assay
  (chain-length × head-group). **NAMING**: there are now TWO opt-in two-pool root models —
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
- **NEUTRAL-organic (Briggs/Kow) path (this session) — IMPLEMENTED + a-priori VALIDATED (root 0.281 / per-organ 0.783)**:
  `src/neutral_dpu.py` + `validation/neutral_dpu_validation.py` + `tests/test_neutral_dpu.py` +
  `docs/neutral_dpu_validation.md` + `data_obs/neutral_obs_ge2017.csv`. The neutral DPU base was DERIVED
  (`dpu_model_summary_corrected.tex`) but never implemented — `Compound.fn` was pinned at 0.0 and used in NO equation,
  the neutral membrane term existed only as a comment, and `data_obs/` held only PFAS. **Why it matters**: every PFAS
  result is entangled with FITTED PFAS-specific params (`f_xy`, `k_seq`, lipid conductances) — hence "reproduction, not
  prediction". A neutral compound has `K_PW` and `TSCF` both fixed from OUTSIDE by log-Kow QSPRs with nothing to tune, so
  it is the ONLY setting where the DPU **backbone** (4 compartments, xylem advection, growth dilution, terminal
  accumulation) is testable apart from the ionic extension. **No new ODE**: a neutral compound is the SAME 4pool ODE with
  `z=0` (⇒ `N=0`, GHK→1, exclusion `e^N` 107→**1**, so the membrane term degenerates EXACTLY to passive Fickian
  `κ_d(Cwo−Cw)`), `Vmax=0` (no carrier), `f_prot=f_cw=0` + `K_PL=a·Kow^b` (so `binding_factors` returns the Trapp/Briggs
  `K_PW = W + L·a·Kow^b` term for term), `f_xy = TSCF(logKow)` (computed not fitted), and **phloem OFF** (the base
  explicitly excludes it; the phloem is an ADDITION of the ionisable extension).
  **ANCHORS VERIFIED AT SOURCE** (10 papers obtained by the user, `DPU4OC.zip`): Briggs 1982 eqs. 2/3 confirmed
  character-for-character, and Briggs derives the RCF floor 0.82 by attributing it to root WATER content (~90% water →
  ~0.9) — independently confirming the `K_PW` first term maps to `Compartment.theta`. Schriever 2020 reprints Briggs eq. 3
  verbatim (second independent verification) AND supplies a 97-value refit (A=0.746, B=2.160, C=7.230, in log D) — wired
  as `tscf_model="schriever"`. Tissue LIPID contents now cited from Trapp 1994 (root 1%, stem/leaf 3% dw) instead of the
  earlier guesses. **Structural checks**: partition adapter reproduces Briggs RCF to **machine precision** (1.3e-16);
  with zero fitted params the straw/root ratio peaks at **exactly logKow 1.78 = the Briggs TSCF peak**; leaf is an
  unbounded terminal accumulator at γ=0 (leaf BAF 194 → 19 at a 7-d half-life, root unchanged 0.51) ⇒ for neutrals
  metabolism is LOAD-BEARING (volatilisation, the leaf's OTHER sink, is now implemented — see the plant-air bullet).
  **A-PRIORI PREDICTIONS (the repo's first), two independent rice datasets, NOTHING fitted**:
  (1) **Liu 2023** (`10.1016/j.scitotenv.2022.159826`; SI Table S1 log Kow + Tables S3xS4 reconstruction of total
  tissue conc — the fractions sum to exactly 1.00 per tissue, so it is arithmetic on published numbers, not digitising;
  cross-checks against the paper's own text TF_L/S statements) → **root partition RCF for 14 un-ionised pesticides
  (7 neonicotinoid + 7 triazole) spanning logKow −0.66→4.4 at log10 RMSE 0.281**, 11/14 within 1.5x over a 50-fold RCF
  range. Root partition is the CLEANEST possible test: no transpiration/duration/metabolism assumption because it is an
  equilibrium. The 7 sulfonylureas are EXCLUDED (weak acids, ionised at the solution pH 5.6–5.8 → ionic path).
  (2) **Ge 2017** (`10.1016/j.envpol.2017.04.043`, per-organ TF at 60 d, dw) → **log10 RMSE 0.783** vs the PFAS a-priori
  0.84–0.95 (which still has FITTED transport behind it). **The error STRUCTURE is the finding**: stem predicted well
  (0.5x/1.3x/3.9x), leaf over by 6–14x — exactly the γ=0 unbounded-leaf failure mode §3 predicts. Sensitivity scan (NOT a
  fit) has a genuine MINIMUM at ~7 d half-life (0.783→0.210, rising again to 0.515 at 3 d) ⇒ the model predicts a
  SPECIFIC in-planta half-life, a testable prediction. The two datasets CROSS-VALIDATE the diagnosis: half-life
  sensitivity is steep where the endpoint accumulates (leaf) and FLAT where it equilibrates (root, 0.281→0.294).
  Also: the ORIGINAL Briggs bell BEATS the broader modern Schriever refit at every half-life (0.783 vs 1.212) — the
  narrow bell is right that lipophilic compounds do not reach the leaf (difenoconazole obs leaf TF 0.044).
  **CORRECTION the data forced**: the first cut converted Trapp's lipid contents as DRY weight; they are FRESH weight
  (same basis as W, as `K_PW = W + L·a·Kow^b` requires, and as Briggs' own anchor `L·a=10^−1.52` ⇒ L=2.5% next to
  W=0.82 shows). Liu's root data caught it — partition-term RMSE **0.605→0.198**, Ge 1.099→0.783. Mixing the bases
  understates K_PW ~10x for lipophilic compounds.
  **Open**: grain compartment UNTESTED (no dataset exists — same gap as PFAS); rice organ lipid still borrowed from
  soybean; Liu's shoot TFs unused (72-h seedling exposure not comparable to a season run); Inao 2018 samples WHOLE
  SHOOT only (no organ split) so it can only test the HYDRUS coupling against a lumped shoot; Trapp 1994's measured
  concentrations are FIGURE-ONLY; **Brunetti 2021's calibrated pea root `K_RW`=13.3 vs the Briggs `K_PW`≈1.0 for
  carbamazepine is an order-of-magnitude disagreement with the partition core** on the framework's own reference
  implementation — open.
- **PLANT-AIR EXCHANGE — volatilisation + gaseous uptake (this session; handoff A1) — DONE**: `src/plant_air.py`
  + an optional `RiceUptakeModel(air=…)` hook + `simulate_neutral(air=True)` + `tests/test_plant_air.py` (13) +
  §3b of `validation/neutral_dpu_validation.py` / `docs/neutral_dpu_validation.md`. **The gap it closes**: the core
  ODE had NO air terms, so the neutral path's leaf — an unbounded terminal accumulator whose only sinks are
  metabolism and volatilisation — had just ONE of the two, making every volatile compound an upper bound BY
  CONSTRUCTION; `k_aw_warning` merely refused such a compound. Implemented straight from the derivation
  (`dpu_model_summary_corrected.tex` §`sec:permeability`): cuticle `eq:Pc` + air boundary layer `eq:Pair` + aqueous
  layer `eq:Paqua` in **series** (`eq:Pctot`), in **parallel** with the transpiration-linked stomatal conductance
  (`eq:Ps`/`eq:Csat`, `eq:Pp`), driving `eq:Qvol` against `eq:Qgas` through `K_PA` (`eq:Kpa`) — keeping the tex's own
  assumptions (**no volatilisation from roots; stem cuticle-only; leaf+grain cuticle+stomata**). Two pitfalls the
  handoff flagged, both handled + tested: the correlations are **SI (m/s, g/mol)** so everything is converted once to
  m/day and the `m³→L` factor is isolated (pinned by the Henry's-law equilibrium test), and `eq:Ps`'s **`1/(1−φ)`
  pole** is capped (`RH_MAX`). **PFAS-safety is STRUCTURAL, not a special case**: `P_air` and `P_S` are both ∝ `K_AW`
  and `P_air` sits in SERIES, so `K_AW=0` zeroes the whole pathway — and the core's default `air=None` means the term
  is never even evaluated. `reproduce_demo` **RMSE 0.029 re-verified**; enabling air on a `K_AW=0` compound is
  **bit-identical** to air-off. **Result** (K_AW ladder, logKow 2.42): leaf BAF **177 → 0.0025** across K_AW 0 → 0.1
  (leaf t½ ∞ → 0.0008 d) while the **root is invariant**; it also independently CHECKS the old judgement-call warning
  threshold (`K_AW>1e-4`: t½ 787 d at 1e-4 vs 8 d at 1e-3 ⇒ flags early, the safe direction). `k_aw_warning` now names
  the remedy instead of refusing. **Honest limit — the SURFACE AREA, not the equations**: flux ∝ specific surface `S`
  [m²/kg], and `S` entered this repo as a leaf/grain RATIO for the xylem split (never calibrated as an absolute area),
  so absolute volatilisation magnitudes are order-of-magnitude until measured (`AirExchange(S=…)` overrides; the
  shipped **stem `S`=0** leaves the stem term inert, pinned by a test). Particle deposition (`eq:Qdep`) deliberately
  NOT implemented (separate deposition pathway; `f_particle` only excludes the particle-bound share from `eq:Qgas`).
  Opt-in and scoped to the 4pool core + neutral path — `parameters.json`, `simulate()` and the PFAS models unchanged.
- **NEUTRAL path exposed through `model_api` (this session; handoff A2) — DONE**: `model_api.simulate_neutral(log_kow, …)`
  mirrors the `simulate_nstem_leaf`/`simulate_twopool_seq` opt-in pattern — same driver machinery (`drivers=`,
  `biomass=`, `measured_forcing=`) and the same result-dict contract as `simulate()` (`conc`/`baf`/`baf_final`/
  `straw`/`straw_baf`/`tf_final`/`cwo_ref`/`season`/`M`/`params`), plus the neutral-only diagnostics (`K_PW`,
  `TSCF`, `rcf_briggs`, `air_summary`). First arg is a **log Kow**, not a congener (a neutral compound has no
  congener analogue). `half_life=` sets γ; `air=`/`phloem=` stay opt-in. **Drift guards** (`test_model_api.py`):
  bit-identical to `neutral_dpu.simulate_neutral` through BOTH the built-in forcings and `drivers=`, so the published
  a-priori numbers (Liu 0.281 / Ge 0.783) cannot silently stop describing what the API returns; plus `N=0`/`e^N=1`
  and air opt-in ≡ zero at `K_AW=0`. **Streamlit: an EXPERT-ONLY "🧪 Neutral organics" tab** (`ui/expert.py` tabs[8]; About moved to
  tabs[9]) — log Kow + name + half-life + TSCF QSPR (briggs/schriever) + opt-in phloem/air, reporting
  TSCF/`K_PW`/three BAFs/tissue dynamics. Expert-only ON PURPOSE: Simple is congener-driven and
  symbol-free while a neutral compound is a log Kow, and the neutral GRAIN is untested — the opposite
  of what a general-audience screen should show absolute numbers for. No `fig_baf` (there is no observed
  neutral series, so a "predicted vs observed" frame would mislead). Scope limits stated in-UI.
  Defaults / `parameters.json` UNCHANGED.
- **Hwang 2017 lettuce/chlorpyrifos (this session; handoff A3) — DIAGNOSIS, not a score**:
  `validation/hwang2017_lettuce.py` + `tests/test_hwang2017.py` (8) + §4c of `docs/neutral_dpu_validation.md`.
  The only TIME-RESOLVED PER-ORGAN neutral dataset to hand (3 samplings × 2 soil levels), LIPOPHILIC (logKow 4.01,
  Briggs falling limb) and the only one with a MEASURED `Kd` — which is what makes it usable, converting the soil
  residue to the pore-water exposure `Ce(t)=C0·(½)^(t/T)/Kd` the model needs. **Its RMSEs (0.610 fw / 0.726 dw) are
  NOT validation results** — four limits stack (lettuce not rice, 1 compound, **Table 1's fw/dw basis UNSTATED in the
  article** ≈20× lever, growth-curve form `Ig`/`Kg` not in the transcription ⇒ reconstructed log-log, roots grew in
  soil ⇒ contact-confounded like Li 2025). **What it DOES establish**: (i) Table 1 is internally consistent on ONE
  basis — `whole` = mass-weighted mean of leaf+root at root mass fraction **5.4±0.9%**, identical at every sampling —
  and 5.4% is characteristic of FRESH weight (dw would be ~11%); (ii) the modelled root cannot exceed its equilibrium
  partition, so `K_PW`=15.8 L/kg is a STRUCTURAL CEILING and the basis flips the verdict (fw: measurement **2.8–10.4×
  ABOVE** ⇒ unreachable; dw: **under** it); (iii) **the two readings fail on OPPOSITE organs** (fw leaf 0.489/root
  0.711; dw root 0.393/leaf 0.948) ⇒ the basis decides WHERE the model is wrong, not WHETHER — not a units artifact;
  (iv) soil contact can't explain the exceedance (would need 12–49% of washed root mass to be soil). **Payoff**: the
  fw root exceedance is the SAME direction/magnitude as the open **Brunetti 2021** `K_RW`=13.3 vs Briggs ≈1.0
  disagreement ⇒ "Briggs root partition too low for LIPOPHILIC compounds in SOIL-GROWN plants" now has **two
  independent sightings**, making it the best-evidenced open question against the partition core. **Trap named in
  code+docs**: the half-life scan improves dw (0.73→0.30) and worsens fw, so the fit "prefers" dw — using that to pick
  the basis is CIRCULAR and overrides the only non-circular evidence (i). `Tp` scanned, never adopted (the authors say
  it is not their measurement). **Deliberately NOT a `data_obs/` CSV**: the shared `--obs` harness would run it on the
  RICE drivers (120 d, constant Cwo) and return a silently meaningless number. `parameters.json`/model math UNCHANGED.
- **Acquisition-queue papers arrived (this session) — two new measured tables, and the ANCHOR they expose**:
  `DPU4OC_add.zip` delivered rows A1–A4 + the C-row candidates of `docs/literature_db/Acquisition_Queue.csv`, which is
  now a RECORD of what each one actually contained (status column per row) rather than a want-list. Three rows did not
  deliver what they asked for and that is recorded explicitly: **A1 contains NO RICE**, the supplied `Mcfarlane1987.pdf`
  is the **wrong paper** (constructed-wetland N isotopes), and **C1/C2/C3 candidates all screen out** (Deng 2018 samples
  brown rice but every final residue is <LOD; the C2 papers are ROOT morphology while the air term takes no root
  contribution and `a_R` is lumped into `kappa_d`; Honda 2023 is compositional mol%, not a total lipid mass fraction).
  A2/A3 CLOSED; **A4 (Kodešová) is SI-NEEDED** and is the highest-value outstanding request (article gives log Kow 2.25,
  pKa 1.0/13.9, Freundlich `KF` for 3 soils ⇒ with Tables S2/S5 it becomes a per-organ a-priori test with the exposure
  pinned by the paper's own sorption data).
  - **`data_obs/neutral_obs_li2019_rcf.csv`** (A1 Table S1) — 48 hydroponic RCF over 11 species, log Kow −0.57→5.41,
    4 rice rows; the 18 Briggs barley rows ship marked `subset=calibration` and are held out. **A-priori log10 RMSE
    0.598 (n=29)**, and the error is **one offset, not scatter** — all 11 species biased LOW (−0.30…−0.95).
  - **`data_obs/tscf_obs_schriever2020.csv`** (A3 SI Table A3, all 97 rows) + `validation/schriever2020_tscf.py` —
    the first test of **TSCF with no plant model in between** (previously only reachable through the Ge 2017 ODE
    comparison, entangled with the unknown half-life). Default Briggs bell: **RMSE 0.310, bias −0.221** on the 30
    un-ionised rows vs a FITTED model's in-sample 0.234. Free by-product: logD-vs-logP lifts rank corr 0.313→0.653,
    confirming Schriever's own claim on their own table.
  - **THE DIAGNOSIS** (`validation/li2019_rcf_apriori.py`): the offset is INTERNAL. `neutral_dpu` anchors on Briggs'
    RCF (`L·a = 10^−1.52 = 0.0302`) but `rice_compartments` substitutes a MEASURED `L=0.01` and keeps the CONVENTIONAL
    `a=1.22` ⇒ product 0.0122, **2.48× below the anchor it claims** (the header itself says only the product is
    identifiable). Cost on Briggs' OWN barley rows: log10 RMSE **0.266 vs the anchor's 0.111**.
  - **NOT kinetic** (§3b, the PRE-REGISTERED confounder): the bias is **flat in exposure time** (−0.442 @1–3 d vs
    −0.447 @>3 d; within logKow>3.5 alone −0.586 vs −0.517), so Li/Chiou's own `α_pt`(non-equilibrium) story is
    RULED OUT; the deficit is in the **lipophilic sorption term** (absent below logKow 2 where the water floor `W`
    dominates, −0.3…−0.6 above). **TWO CORRECTIONS from an adversarial re-read (§3c)**: (a) the **monotone** ladder is
    NOT robust — Namiki 2015 alone supplies 10/29 rows, all in the top two bins as 2 compounds × 5 species; collapsing
    compound×study gives −0.030/−0.305/−0.576/**−0.501**, i.e. the top two bins go FLAT. What survives is the low-Kow
    zero + all TEN studies biased the same way (−0.20…−0.89). (b) the earlier claim that "a flat lipid rise is the
    wrong instrument for a Kow-dependent deficit" was **WRONG** — `K_PW=W+L·a·Kow^b` is water-floor-dominated at low
    Kow, so scaling `L` IS inherently Kow-dependent (+0.045 log @logKow 1 → +0.391 @5), close to the observed shape and
    worth 60–75% of it ⇒ **Li 2019 is genuine evidence FOR the anchor**.
  - **DEFAULT UNCHANGED — deliberately.** Restoring the anchor improves Li 2019 (0.598→0.331) and Ge 2017
    (0.783→0.651) but DEGRADES **Liu 2023, the only RICE table** (0.281→0.288); an intermediate `L≈0.015` fits all
    three and would make the path's one real claim (nothing is fitted) false. `ND.BRIGGS_ANCHORED_LIPID_FW` makes the
    alternative runnable; `tests/test_li2019_schriever_tables.py` (11) pins the 2.48× so it cannot drift back.
    A1's real payoff is the **DEFINITION**, verified at source: Li/Chiou's `f_lip` is **fresh weight** (dry converted at
    90% root water) inside an RCF expression of this model's form, and their cereals (barley 1.00 / wheat 1.10–1.14 /
    maize 0.53 %) BRACKET the 1% in use ⇒ the value is corroborated, provenance upgraded soybean→cereal root. **Briggs
    1982 measures no lipid at all** (verified: he attributes the 0.82 floor to root WATER), so 1.00% and 2.47% are both
    inferences from a paper that measured neither; the 2.5× excess reads as **non-lipid sorption** (cell wall/lignin =
    the PFAS side's `f_cw·K_cw`, GAP A) which the neutral composition zeroes.
  - **The DATASETS CONTRADICT EACH OTHER — this is the real obstacle (§3d)**: at logKow 3.5–4.5, Li 2019 says the
    model is ~3× LOW (−0.462, n=11) while **Liu 2023 (rice, same hydroponic endpoint) says it is EXACT** (−0.008, n=5);
    on **propiconazole, the one compound BOTH measured** (logKow 3.72) they report RCF **43.65 (lettuce) vs 9.32
    (rice) = 4.7×**, i.e. the spread BETWEEN measurements EXCEEDS the 2.4× the anchor is worth. ⇒ the open question is
    NOT "is the partition too low" but **"which hydroponic dataset describes a RICE root above logKow 3.5"**, and
    nothing in-repo settles it — only a rice measurement in that range. A measured **neutral-organic cell-wall
    coefficient** still serves BOTH paths (PFAS GAP A) and stays the top wet-lab item. `parameters.json`, `simulate()` and `reproduce_demo`
    (RMSE 0.029) UNCHANGED; the `subset` filter is inert on tables without the column so **Liu 0.281 / Ge 0.783 are
    bit-identical** (pinned by a test).
- **Kodešová 2019 SI arrived (this session) — queue A4 CLOSED; the anchor vote flips and Brunetti is superseded**:
  `src/../data_obs/neutral_obs_kodesova2019.csv` + `validation/kodesova2019_carbamazepine.py` +
  `docs/neutral_dpu_validation.md` §4f. Carbamazepine root partition, **4 plants × 3 soils × 2 treatments (n=21)**,
  with the **cleanest exposure in the repo**: root conc (SI Tab S2, ng/g **dw**, basis stated) ÷ pore water derived
  from the SOIL conc measured on the SAME pot at the SAME harvest (Tab S4) via the paper's OWN measured Freundlich
  isotherms (`c=(C_soil/K_F)^n`) — no mass balance, no pot geometry, no dissipation model. CAR is un-ionised
  everywhere (pKa 1.0/13.9) and `DT50>1000 d`, so the authors' own reason for not computing BAFs does not apply.
  - **A-priori log10 RMSE 0.191 (n=21), nothing fitted** (Liu 0.281, Li2019 0.598, Ge 0.783) — **but that 0.191 is
    FLATTERED by a scoring artifact, see the next bullet**; on the appropriate basis it is 0.237 and Liu (0.206) is
    the repo's best a-priori result.
  - **VOTES AGAINST restoring the Briggs anchor**, and is well-conditioned to (at log Kow 2.25 the two compositions
    differ 1.6×): shipped/anchored = Kodešová **0.191**/0.216 · Liu **0.281**/0.288 · Li2019 0.598/**0.331**.
    ⇒ tally **2 tables against the anchor, 1 for it**, and the two against are the soil-grown / rice ones.
    Driver-free cross-check agrees in ORDERING (measured `RCF_fw` median **1.10** vs `K_PW` 1.56 shipped / 2.53
    anchored — the shipped root already runs HIGH here); survives the dw→fw lever (θ 0.85–0.95 keeps the anchor
    0.21 log worse).
  - **SUPERSEDES the Brunetti sighting**: Brunetti's calibrated pea `K_RW`=13.3 is **~12× ABOVE a direct measurement
    of the same compound** (1.10), while Briggs lands within ~1.5×. So that disagreement sits in the calibration,
    not the partition core ⇒ the "four sightings" synthesis is REWRITTEN: the deficit is **NOT global, it is confined
    to log Kow ≳ 3** (Li2019 bias −0.03 below 2 → −0.69 above 4.5; Kodešová at 2.25 is POSITIVE).
  - **Pivotal assumption, flagged not buried**: the Freundlich unit reading. Defended by the `Koc` it implies
    (**222/189/154** across 3 soils spanning 3.8× in organic carbon — inside CAR's literature band) and by the derived
    pore water (0.10–0.70 mg/L) matching both the applied solution (~0.5–1.0 mg/L) and the sorption study's 0.5–10 mg/L
    range; the alternative g/cm³ reading gives `Koc`~1600 and would REVERSE the vote. Tests pin the derivation, the
    `Koc` band, and the opposing votes (`tests/test_li2019_schriever_tables.py`, 15).
  - Limits: 4 leafy/root vegetables, **no rice**; one compound at one lipophilicity (says nothing about the high-Kow
    end where the deficit lives); roots rinsed not exhaustively cleaned (biases obs UP, i.e. against this conclusion).
  `parameters.json`, `simulate()` and `reproduce_demo` (0.029) UNCHANGED.
- **Li 2019 Table S2 (376 SOIL rows) + Kodešová's LEAF half (this session) — the sign FLIPS, and the open
  question is REFRAMED**: `data_obs/neutral_obs_li2019_soil.csv` + `validation/li2019_soil_table.py`, and §6/§7 of
  `validation/kodesova2019_carbamazepine.py`; docs §4h/§4i.
  - **THE SIGN FLIP**: the SAME paper's soil table (n=376, 13 crops) runs the model **HIGH +0.260** where its
    hydroponic half (n=29) runs it **LOW −0.432**. ⇒ "the root partition is too low" is a property of ONE EXPOSURE
    ROUTE, not of the partition core. Restoring the anchor makes the soil table much worse (+0.645, RMSE 0.639→0.873)
    ⇒ tally is now **3 tables against the anchor, 1 for**, and the 1 for is the hydroponic half of a paper whose own
    soil half says the opposite.
  - **The COMPOSITION TERM works**: substituting each crop's OWN measured lipid (11× spread, radish 0.09% → wheat
    1.14%) takes the bias **+0.260 → −0.001** (RMSE 0.639→0.461). So `K_PW = W + L·a·Kow^b`'s FORM is right across a
    large range in `L`. (Caveat: `f_lip` is Li's own model input ⇒ consistency, not independent validation.)
  - **MOST OF THE ERROR IS THE EXPOSURE, NOT THE PLANT**: Li derive `value` as soil conc ÷ `K_om`, and their Table S3
    says which `K_om` was measured. **Experimental `K_om` (n=62): bias +0.033 — essentially unbiased**; QSPR `K_om`
    (n=261): +0.291; and the f_om gradient agrees (+0.499 <1% OM → +0.101 >4%). Combined with Liu (−0.053) and
    Kodešová (+0.162), **on all three tables with a measured/known exposure the bias is −0.05…+0.16** ⇒ the neutral
    path's apparent partition error is largely a **SOIL-SIDE exposure problem**, not a `K_PW` problem. `α_pt` is NOT
    applied (Li median 0.098): what Li put in `α_pt`, Briggs put in the flatter exponent b=0.77 vs Li's `Kow^1.03`
    (20× apart @logKow 5) — dividing would double-count.
  - **Kodešová LEAF (§4i)**: leaf/root cancels the exposure ⇒ a pure TRANSLOCATION test. Measured median **3.25**
    (0.31–9.05) vs model **181** (~55× over) — §3's terminal-leaf runaway, now measured per-organ. **Metabolism CANNOT
    close it**: even a 1.5-d half-life leaves 28.5 ⇒ most of the excess is the **rice-driver mismatch** (a rice season's
    transpiration per unit leaf mass on a 340 cm³ pot of lettuce). ⇒ **WARNING on the Ge 2017 "≈7 d half-life"**: it may
    be absorbing the same mismatch rather than measuring metabolism.
  - **METABOLISM MEASURED, not fitted (first in the repo)**: Kodešová quantified CAR's 4 metabolites, so the parent
    fraction is a direct observation — root **0.919**, leaf **0.489** (lamb's lettuce 0.169, radish 0.810). ⇒ (i) `γ=0`
    is WRONG for carbamazepine despite a soil `DT50>1000 d` — **soil persistence ≠ plant persistence**; (ii) it is
    strongly SPECIES-dependent (4.8×, same compound/soils/harvest) ⇒ an in-planta half-life is **not a compound
    property**, which is what fitting one to a single dataset assumes. The epoxide EXCEEDS the parent in lamb's-lettuce
    leaves (17000 vs 6400 ng/g) ⇒ a parent-only model understates the burden several-fold.
  `tests/test_li2019_schriever_tables.py` (24). `parameters.json`, `simulate()` and `reproduce_demo` (0.029) UNCHANGED.
- **Briggs 1983 — the STEM anchor that was never read (this session)**: `neutral_dpu.briggs_scf` /
  `briggs_stem_xylem_partition` + `validation/briggs1983_stem.py`, docs §4j. The companion to Briggs 1982 sat unread
  in the obtained set; it fits the SAME `K_PW` form to barley shoots (VERIFIED AT SOURCE):
  `log(K_stem/xylem_sap − 0.82) = 0.95·logKow − 2.05` (eq.2), `SCF = K_stem/xylem × TSCF` (eq.3). Self-checking:
  computed from the coefficients the peak is **6.4 @logKow 4.5** vs the paper's "about 6 … at about 4.5". ⇒ the stem
  had **NO anchor** in this repo (Trapp soybean 3% lipid + conventional a=1.22 + the ROOT exponent 0.77):
  **root `L·a` 0.0122 vs anchor 0.0302 = 2.5× BELOW; stem 0.0366 vs 0.0089 = 4.1× ABOVE** — two organs from unrelated
  sources, neither checked, missing in OPPOSITE directions. **BUT the consequences differ and only §3 should be
  quoted**: for the root the coefficient gap IS the disagreement (shared exponent), while for the stem it largely
  **CANCELS** against Briggs' steeper b=0.95 — the observable SCF differs by at most **0.13 log** over logKow 0–3.5
  (where TSCF delivers anything), swinging to −0.20/−0.38 only above 4.5 where TSCF has collapsed and Briggs' own text
  says the decline "was not tested". ⇒ the stem is a **provenance** problem, not a prediction problem; ranks below the
  root question and the exposure-term work. NOTHING CHANGED (one species, a shoot BASE not a true stem, and no measured
  table in-repo would arbitrate it).
- **ODE-vs-EQUILIBRIUM scoring artifact (this session) — `compare_to_obs(mode="equilibrium")`**: all three root-partition
  tables (Liu/Li2019/Kodešová) measure an EQUILIBRIUM over 24 h–26 d, but were scored by running the **120-d rice season**
  and reading `baf_final["root"]`. That imposes a **Kow-DEPENDENT, purely model-side discount** (ODE/K_PW = 0.91 @logKow
  −0.5 · **0.55 @1.78** · **0.57 @2.25** · 0.99 @5.0 — the xylem drains the root hardest near the TSCF peak). Applying a
  rice season's drain to a 24-h barley measurement is an ERROR, not a modelling choice. Rescored: **Liu 0.281→0.206 ·
  Li2019 0.598→0.541 · Kodešová 0.191→0.237** (Ge unchanged — per-organ TF at 60 d, the season IS its endpoint).
  ⇒ **Kodešová was flattered** (it sits at logKow 2.25, near the discount peak) so it is NOT the repo's best a-priori
  result — **Liu 0.206 is**. The §4f anchor verdicts SURVIVE and their **margins widen** (shipped/anchored on the
  equilibrium basis: Liu **0.206**/0.286 · Li2019 0.541/**0.295** · Kodešová **0.237**/0.410) ⇒ removing the artifact
  **SHARPENS the dataset contradiction rather than resolving it**. `mode="ode"` stays the DEFAULT so no published number
  moves silently; which basis should headline is a one-line decision left open with the anchor question.
  `tests/test_li2019_schriever_tables.py` (18) pins the discount shape, the rescored ordering, and that the default is
  unchanged.
- **Structural MERGE — two-pool seq ROOT + `nstem_leaf` redistributed SHOOT (this session) — Result 7 CONFIRMED**:
  the last in-silico item of the two-pool arc (handoff §6: "a fair per-organ Tang test needs the two-pool root merged
  with the redistributed shoot"). `NStemLeafModel` gained optional `k_seq`/`k_rel`: `k_seq>0` APPENDS a sequestered
  root state (irreversible apoplast/cell-wall sink sharing the root mass) so the reported root = mobile+seq;
  **`k_seq=0` (default) adds NO state** → the pre-merge model, every default and every existing number are untouched.
  Wired as `model_api.simulate_twopool_nstem` (root params from the cached Yamazaki fit + monotone physical `f_xy`;
  `twopool_seq_params` exposes that set alone); `simulate_nstem_leaf` also gained `drivers=`, `K_cw_organ`,
  `straw`/`straw_baf`. **Perf**: the RHS rebuilt `binding_factors` and interpolated the driver matrices column-by-column
  on every one of ~3k calls/solve — caching the former + one `axis=0` interpolant per matrix HALVES the solve time,
  Tang TFs identical to 6 dp. **Forcings matter structurally**: the redistributed shoot is transpiration-DEPOSITION fed,
  so on the demo forcings (Q_TP peak 0.40, ~4× measured) the deposition route floods the grain and the re-fit collapses
  onto its bounds (0.658) — the merge is run on the MEASURED forcings (`--demo` reproduces the pathology).
  **Results** (`validation/twopool_nstem_merge.py`, `docs/twopool_root_exploration.md` §Result 8): in-sample Yamazaki
  **0.301** (root 0.153) vs the same root fit behind its own pass-through stem 0.278 — and transferring that fit onto the
  new shoot with **NO re-fit at all** already gives 0.316, i.e. **swapping the entire shoot costs 0.04 log units** and the
  root RMSE is unchanged to 3 dp ⇒ the root mechanism is genuinely separable from the shoot model. Non-K_PL PFOS/PFUnDA
  separation survives and sharpens (k_seq 0.047 vs 0.188 = **4.0×** at identical K_PL). **Tang per-organ OOS
  (no Tang re-fit): 1.398 → 0.801**, recovery carried by the diagnosed organ (**stalk 1.89 → 0.61**; leaf 0.38 → **0.28**,
  best of any model) ⇒ **Result 7's diagnosis is CONFIRMED — the two-pool's Tang failure was a SHOOT artifact.**
  **Honest**: it does NOT beat single-pool lipid loading (0.516, still the Tang winner); the whole residual is the
  **endosperm** (1.21, the documented structural grain under-prediction). The merge also EXPOSES a real tension: Tang's
  TFs are high (stalk TF 2.2 > root) rewarding strong translocation, while Yamazaki's high long-chain root demands a
  retaining `k_seq` — one (f_xy, k_seq) cannot satisfy both (same condition-dependence as `f_xy` PFOS 0.14 vs 0.32).
  Tang is still only 3 C5–C8 congeners, so the long-chain root decoupling remains unexercised (Kim grain stays the
  informative two-pool OOS) ⇒ strengthens the structural case but does NOT move the promotion decision (gate = the §5
  wet-lab assay). Guards `tests/test_twopool_nstem_merge.py` (+ nstem_leaf/model_api). `parameters.json`,
  `simulate()` and `reproduce_demo` (0.029) UNCHANGED.
- **`k_seq` mechanistic provenance — literature synthesis (this session; handoff item ②/§4.5) — PARTIALLY SUPPORTED**:
  `docs/twopool_kseq_mechanism.md` — a fan-out, adversarially-verified deep-research synthesis (17 sources → 66
  candidate claims → 25 verified by 3-vote panels) answering whether the phenomenological U-shaped, non-K_PL,
  head-group-dependent `k_seq` has a physical basis. **Verdict: PARTIALLY SUPPORTED — and it explains the model's
  two-arm form.** Key results: (1) the **U-shape is a SUPERPOSITION of two distinct mechanisms** — a short-chain
  **electrostatic/anion-exchange** arm (biphasic soil Koc; short chains over-sorb the hydrophobic QSPR) + a long-chain
  **hydrophobic, desorption-resistant** arm (irreversibility index TII rising with chain length, ≈0.98 at C10 PFDA) —
  which is exactly why the fitted `k_seq` is a **sum of two exponential arms** and why a single linear `k_seq` collapsed;
  (2) the **irreversible sink** (terminal accumulator, not a reversible K) is justified by the measured desorption
  hysteresis; (3) the model's **PFSA offset `10^+0.18`** ≈ the independently **measured +0.23 log** sulfonate-vs-
  carboxylate offset on lignin/soil (non-lipid, electrostatic) — a quantitative, non-circular corroboration; (4) **Fe-Mn
  root plaque is DEMOTED** as the irreversible-sink candidate (PFAS-ferrihydrite binding is outer-sphere/pH-reversible/
  acidic-only/monotonic; flooded paddies are circumneutral) in favour of **cell-wall entrapment** — though a minor
  PFCA-specific inner-sphere plaque term and the arsenate "molecular-sieve" precedent stay open; (5) **central data gap =
  the decisive experiment**: NO rice-root cell-wall or Fe-Mn-plaque PFAS coefficient resolved by BOTH chain length AND
  head group exists (every anchor is an analog matrix — lignin/ferrihydrite/sediment/soil OC), so `k_seq` is anchored by
  **analogy + superposition**, not a direct root-sorbent dataset. ⇒ strengthens the mechanistic *story* but does NOT by
  itself warrant promoting `k_seq` into `parameters.json` (the §5 cell-wall/plaque batch-sorption + desorption assay is
  the gate). Literature/theory only; `parameters.json`, the model math and `reproduce_demo` (RMSE 0.029) UNCHANGED.
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
- **Root-lipid anchor promoted from a buried constant to a NAMED MODE `lipid_source` (this session)**: the neutral
  path's root lipid was a single constant (`L=0.01`) with the alternative reachable only by importing a dict, even
  though the choice is genuinely unsettled and moves a lipophilic compound's root partition up to **2.5×** (K_PW ratio
  1.02× at logKow 0 → 2.46× at 5; only the PRODUCT `L·a` is identifiable, and **`a=1.22` has no citation in the repo**).
  Now `neutral_dpu.LIPID_SOURCES` = **`"measured"`** (L_root=1% fw, measured cereal roots — Li2019 barley 1.00/wheat
  1.10–1.14/maize 0.53%; DEFAULT) vs **`"briggs_anchor"`** (2.47% fw = `L·a=10^−1.52` exactly, i.e. what Briggs' 1982
  barley regression implies — **Briggs measured no lipid at all**), threaded through `rice_compartments(lipid_source=)`
  → `neutral_dpu.simulate_neutral` → `model_api.simulate_neutral` (reported in `params["lipid_source"]`) →
  `validation.compare_to_obs`/`equilibrium_rcf`. The repo's own idiom (`f_xy_source`, `cwo_profile`, `biomass`,
  `tscf_model`, `mode`); an explicit `lipids=` dict still overrides. **New `compare_lipid_sources()` +
  `--lipid-source both` scores all 5 shipped tables under BOTH readings side by side**, so the anchor evidence is a
  reproducible command instead of a doc claim. **DECISION: default stays `"measured"` — because the evidence no longer
  supports moving, NOT because 1% fw is validated for rice** (it is not; the rice-specific organ lipid gap is open, and
  stem/leaf are still Trapp's soybean values). Evidence is **3 tables against restoring the anchor, 1 for**, and the one
  for is the hydroponic half of a paper whose own 12×-larger soil half says the opposite; the anchor also raises
  predicted root uptake on tables where the model ALREADY runs high. **Nothing moved: Liu 0.281 / Ge 0.783 /
  `reproduce_demo` 0.029 re-verified bit-identical**; the two prior anchor tests were switched off their global
  `TRAPP1994_LIPID_FW` mutation onto the mode. **Docs §5/§6 also made consistent with §4f/§4g** — the "Kodešová 0.191 is
  the best a-priori result" line and the "Brunetti is the best-evidenced open question" synthesis had been retracted in
  the body but SURVIVED in the Honest-summary section (the part most likely to be read alone); both now corrected
  (Liu 0.206 on the equilibrium basis is best; the open anomaly is Li 2019's hydroponic half), and every quoted neutral
  number now names its `lipid_source`/`mode`. Tests: `test_li2019_schriever_tables.py` (default is a no-op; the 2.48×
  and the Kow dependence reach `model_api`; the 3-vs-1 tally pinned).
- **Briggs 1983 Table 1 mined — the repo's FIRST measured STEM test (this session)**: `data_obs/
  neutral_obs_briggs1983_shoot.csv` + `validation/briggs1983_shoot.py` + `tests/test_briggs1983_shoot.py` (9) +
  docs §4k. §4j had compared the repo's stem to Briggs 1983's *fitted equations*; its **Table 1** — the data behind
  them — was never transcribed. 16 non-ionised chemicals (7 O-methylcarbamoyloximes + 9 phenylureas), log Kow
  −0.57→3.7, barley shoots cut into leaf/central-stem/stem-base at 24 and 48 h, from a **nutrient solution of known
  concentration** (the cleanest exposure class, like Liu/Li-hydroponic). Table 1's % distribution × total dpm ÷ the
  section fresh weights (§2.2: 0.54/0.40/1.2 g, six plants) ÷ the solution conc **reconstructs the Stem Concentration
  Factor the article itself defines and plots**. **Reconstruction verified 3×** before scoring: eq.(3) from the
  transcribed coefficients peaks at SCF **6.39 @ logKow 4.43** (paper: "about 6 … at about 4.5"); the reconstructed
  central-stem SCF sits on that curve at bias **+0.049** (the paper says its points "fit quite well"); the stem BASE
  runs high (+0.188) exactly as the article predicts from direct contact with the treating solution. All 32
  (compound×harvest) rows sum to 100±0.1% — which caught **two wrong digits in the PDF text layer** (4-chlorophenylurea
  24h Upper 77.4 not 17.4; phenylurea 48h Middle 3.7 not 3.1), confirmed against the rendered page.
  **RESULT — a-priori, nothing fitted: the repo's stem `K_PW`×TSCF predicts the measured central-stem SCF at log10
  RMSE 0.299 (bias −0.030), statistically indistinguishable from Briggs' OWN FITTED eq.(3) at 0.282** — i.e. at the
  level the root tables reach (Liu 0.281, Li2019 0.598). **This settles §4j empirically**: the 4.1× stem-coefficient
  gap really does cancel in the observable, so it is a **provenance** problem, not a prediction problem (the rice culm
  composition remains unmeasured — it is just not *wrong* where checkable). **NEGATIVE finding, equally important**:
  the table does **NOT** constrain the stem/leaf **split**, which is what handoff item 4 wanted it for — the stem's
  share of the shoot burden moves a **median 1.55× (max 2.11×) between two harvests one day apart**, because the stem
  equilibrates and the leaf does not; scored as an equilibrium the leaf's bias **GROWS +0.229 (24h) → +0.586 (48h)**,
  the terminal-accumulator signature (§3), which the article states independently ("leaf amounts generally increased up
  to 72 or 96 h"). So the leaf half is a **second independent confirmation of the terminal-accumulator structure**
  (different species/exposure from Ge 2017), not a split constraint. The 3 compounds the article flags reproduce its own
  diagnoses (aldicarb 3.6× high = in-planta oxidation, "about three times"; aldoxycarb 3.7× = the logKow<0 TSCF
  underestimate the paper names; 4-(4-bromophenoxy)phenylurea excluded, never equilibrated). Scope: barley seedlings,
  24–48 h, one lab, and §2.2's *typical* section weights applied to every test. **Deliberately NOT wired into `--obs`**
  (the Hwang §4c reason: 48-h barley sections vs a 120-day rice season). Independent of the open `lipid_source`
  question (that anchor touches only the root — pinned by a test). `parameters.json`, `simulate()` and
  `reproduce_demo` (0.029) UNCHANGED.
- **WEAK-ELECTROLYTE speciation ported from PR #54 — the two parallel neutral implementations reconciled
  (this session)**: `src/pfas_rice_plant_module_4pool_surf.py` (speciation block) + `src/literature_params.py`
  (`speciation`/`ion_trap_factor`/`neutral_pathway_ratio`) + `model_api.simulate_neutral(pKa=…, is_acid=…, pH=…)`
  + `tests/test_weak_electrolyte.py` (13). **The problem**: TWO independent neutral-organic implementations existed —
  this repo's `neutral_dpu.py` (separate module, reaches a neutral by `z=0`; carries all 5 measured tables and every
  published a-priori number) and **PR #54/#55** (extends the PFAS core in place with an `(fn, fd)` speciation pair;
  structural verification only, `dirty`/33 commits behind main). They overlapped on air exchange, `simulate_neutral`
  and the Briggs lipid partition — but #54 could do **one thing the merged path cannot: a WEAK ELECTROLYTE**, which is
  a neutral molecule AND an ion simultaneously and so cannot be expressed by one global valence (`z` must be 0 or −1;
  the compound is genuinely both). **ONLY that capability was ported**; the duplicated parts were dropped in favour of
  what already ships. Ported: `(fn,fd)`-weighted 3-pathway `root_uptake` (GHK on the ION term only), `Environment.N_for`
  (valence moves onto the COMPOUND — a weak base is a cation, z=+1), `Compound.pKa/is_acid/z/P_n`, `Compartment.pH`,
  and `RiceUptakeModel.phloem_loading_factor()` (leaf→sieve-tube pH ion trap, weighted by `w=Π/(1+Π)`,
  `Π=(P_n/P_d)·f_n/f_d`). NOT ported: the air block (main has `src/plant_air.py`), `f_lip`/`K_lip` (main's
  `neutral_compartment` carries the Briggs lipid on `f_PL`), the Briggs QSPR duplicates. **Continuity is the safety
  property**: `_weak_electrolyte_kw` sets `P_n = kappa_d` so a weak acid with pKa ≫ pH reproduces the `pKa=None`
  neutral path to <0.1%, i.e. `pKa=` extends this model rather than being a second one; the ION gets
  `kappa_d·10^−3.5` (Trapp `P_N_OVER_P_D`). **Verified bit-identical** (exact `==`, 67 floats over 6 scenarios):
  PFAS `simulate()` for PFBA/PFOA/PFOS/PFDoDA/GenX × recommended/W2fit, AND the neutral path — so Liu 0.281 /
  Ge 0.783 / stem 0.299 / `reproduce_demo` 0.029 all still describe what the code does. **New physics reachable**:
  root BAF collapses 1.51→0.0001 as pKa falls 12→−3, and at the SAME pKa a weak BASE gives 1.51 vs the acid's 0.079
  (~19×) because the cation is ATTRACTED by the inside-negative membrane. The `ion_trap_factor` correction is pinned
  by a test: Λ → `10^ΔpH` = 6.31 as pKa falls, **NOT** 1, so multiplying `L_Ph` by it would hand a permanent anion a
  spurious 6.3× phloem enrichment — the trap switches off **kinetically** (`neutral_pathway_ratio` 2 → 1e-7, a 10⁷
  collapse), never thermodynamically. ~~**NOT VALIDATED**~~ — **superseded: it has now been tested, see the next
  bullet.** `parameters.json` UNCHANGED.
- **Weak electrolyte TESTED — direction SUPPORTED, magnitude REFUTED (this session)**:
  `validation/weak_electrolyte_tscf.py` + `tests/test_weak_electrolyte_tscf.py` (13) + `docs/neutral_dpu_validation.md`
  §4l. The port shipped labelled "no measured weak-electrolyte dataset exists here" — true of **rice** only. §4e scores
  just the **30** rows of Schriever's Table A 3 flagged un-ionised and sets the other **67 ionisable ones** aside as
  "outside this model's stated scope"; the port is exactly what extends the scope to them, so **the test data were
  already in the repo before the capability was**. **`f_n` is derived from the table's OWN logD** (`f_n = 10^(logD−logP)`),
  not its pKa: the two columns disagree by a median **1.41 log** and 20 rows have logD ABOVE logP, which no single-centre
  acid or base can produce. That choice is load-bearing — the pKa route manufactures a **false counterexample** out of the
  8 pKa-1.62 barley rows (TSCF 0.63–0.98, but their own logD says un-ionised and the shipped flag already classes them
  neutral); pinned by a test so it is not rediscovered. **Result**: measured transfer DOES rise with the neutral fraction
  (Spearman **+0.480**, n=67 — the port's first empirical support of any kind, and the sign could have come out flat), and
  speciation ON nearly DOUBLES the model's rank correlation (**+0.284 → +0.520**) — but its influx conductance `Φ` moves
  **~1.6e4-fold** across this table where the measurements move **~3-fold**, so it under-delivers (bias +0.023 → **−0.203**)
  and at `f_n<1e−3` predicts nothing where the measured mean is 0.127. **The two metrics are NOT equally solid**: bootstrapped
  (n=4000) the rank gain survives **94%** of resamples and the RMSE loss only **82%** — a subsample flips the RMSE and one
  did, so the guard test asserts the ordering and the under-delivery but deliberately NOT the RMSE. **Cause is structural and
  already known one directory away**: the model's only entry is transmembrane, while real ions arrive **apoplastically** and
  the PFAS side needed a **fitted carrier** for exactly this reason (§2). ⇒ `pKa=` moves from UNVALIDATED to **BOUNDED** —
  usable for the DIRECTION of a speciation effect, not its size, and not below `f_n≈0.1`. Limits: 16 species, none rice; no
  compound names in Table A 3; acid/base unlabelled so both readings are reported. Opt-in; `parameters.json`, `simulate()`
  and `reproduce_demo` (0.029) UNCHANGED.
- **Apoplastic bypass `g_apo` — the mechanism §4l pointed at; PARTIAL repair (this session)**:
  `Compound.g_apo`/`NeutralCompound.g_apo` + `validation/weak_electrolyte_tscf.py` §5 +
  `docs/neutral_dpu_validation.md` §4m. `root_uptake` gains a FOURTH parallel pathway `g_apo·(Cwo−Cw)`,
  **defined by what it does NOT feel** — no `(f_n,f_d)` weighting, no GHK factor — because a route AROUND
  the membrane cannot be gated by speciation or membrane potential (that is why it is a separate term, not a
  bigger `kappa_d`). Rice is the plant apoplastic "bypass flow" is best documented in. **Self-targeting,
  measured not argued**: conductance sets how FAST the root equilibrates, not the LEVEL, so at `g_apo=2` an
  un-ionised compound moves **+3%** and an `f_n=1e−4` compound **750×** ⇒ one parameter can address the
  ionisable rows without disturbing the 30 un-ionised ones (their RMSE 0.379→0.377). **Result — the two
  optima are in DIFFERENT places, which is the whole answer to the pre-registered question**: a SMALL bypass
  (`g_apo≈0.5`) improves ordering AND scale together (+0.520→**+0.635**, RMSE 0.304→0.291) in **99.6%** of
  bootstrap resamples — the first structural change on this arc to improve both at once — but the
  **RMSE-optimal `g_apo=5`** reaches RMSE **0.245** (beating speciation-OFF's 0.272) by degrading the ordering
  to **+0.450**, worse than the peak in **96.9%** of resamples. That is exactly the **absorption** the handoff
  pre-registered. ⇒ **do NOT fit `g_apo` on RMSE**; the defensible region `g_apo≈0.5–1` rests on a **Pareto**
  argument, not a fit. **NOTHING ADOPTED**: `g_apo` defaults to **0** (structurally absent, not merely zero),
  PFAS untouched (its anion deficit is carried by the fitted CARRIER `Vmax_in` — whether one mechanism should
  serve both paths is the open question this leaves), model still under-delivers at its best point
  (bias −0.184), and one 67-row non-rice table cannot pin a structural conductance. `parameters.json`,
  `simulate()` and `reproduce_demo` (0.029) UNCHANGED.
- **Is the CARRIER necessary? — the repo's one addition to Trapp, tested (this session)**:
  `validation/carrier_vs_bypass.py` + `tests/test_carrier_vs_bypass.py` + `docs/HANDOFF_carrier_vs_bypass.md`.
  `docs/theory_anchor.tex` states the 4-compartment DPU (Rein 2011, Brunetti 2019) is **at the membrane level
  the Trapp (2000, 2004) ionizable-compound cell model**, and marks only `f_xy`, `η` and the **Michaelis–Menten
  carrier** as `— (new)`; `f_xy` is a lumped stand-in for Trapp's own detailed root model and `B_k` is his
  `K_RW` re-parameterised, so **the carrier is the single piece of physics this repo adds** — and it was FITTED
  ("fixed during W2 fit"), never compared with an alternative. Three arms, **one global parameter each**
  (checked, not assumed: `Vmax_in` is a single global 20.0, NOT per-congener), on Yamazaki with the a-priori
  monotone `f_xy` fixed so only the uptake term differs: **nothing 2.640 · A carrier 1.035 · B bypass
  (`g_apo`=20) 0.996 · C depolarisation (`E_m`=−90 mV, NO new term) 2.289** (log10 RMSE).
  **(1) An addition IS necessary — arm C REFUTED**: pushing `E_m` to the far end of its own recorded plausible
  range (`e^N` 107→33) barely moves 2.640→2.289, so the lever already inside Trapp's framework cannot carry
  PFAS and the repo's extension is justified in EXISTENCE. **(2) WHICH addition is UNDECIDED**: A vs B is 0.04
  log over 33 obs, bootstrap **P=0.749** — not a difference. **(3) The bypass's own distinguishing claim FAILS
  (pre-registered item 1 REFUTED)**: per-congener `g_apo` trends with chain length (**corr +0.832, spread 25×**)
  where `theory_anchor.tex` says η (which contains the apoplastic bypass ε) is "essentially independent of tail
  length" — the **first data contradiction of that claim**. **(4) THE CONVERGENCE IS THE FINDING**: LC6's
  per-congener Vmax multiplier (~flat to C10, 2.0× C11, 5.5× C12) and this `g_apo` (2–5 for C4–C8, 20–50 for
  C9–C12) demand the SAME chain-length correction — a requirement common to two different entry terms is not a
  property of entry, it belongs to `B_k`/`φ_free` (sequestration), where the two-pool work put it; pre-registered
  item 3 CONFIRMED from the other side (A and B miss long chain about equally, 1.640 vs 1.439, so no winner read
  there). **NOTHING ADOPTED** — the carrier keeps its place **by default, not by evidence**; `g_apo`=0 for PFAS;
  `parameters.json`, `simulate()` defaults and `reproduce_demo` (0.029) UNCHANGED. `simulate` gained
  `vmax_scale`/`g_apo` in the existing override idiom (defaults bit-identical). Absolute levels are
  a-priori-limited (monotone `f_xy`) and NOT comparable to the fitted 0.029 or the documented a-priori ~0.84.
  **Next (agreed): A1** = expose the bypass as a named mode, default carrier; **then B1** = separate the two on
  Tang's 5-dose series (a carrier saturates, a bypass is linear; data already in `raw_si/`) — see the handoff.
- **A1 DONE — the bypass is a named mode, `simulate(uptake="carrier"|"bypass")` (this session)**:
  `model_api.UPTAKE_MODES` + `DEFAULT_UPTAKE="carrier"`. The repo's idiom for a question that is open rather
  than settled (`lipid_source`, `f_xy_source`, `cwo_profile`, `biomass`, `tscf_model`): **the default does not
  move and the alternative is kept RUNNABLE instead of buried in a validation script**. `"bypass"` = the scored
  arm exactly (carrier off, one global `g_apo`=20); explicit `vmax_scale=`/`g_apo=` still win over the mode so
  the scans keep working; the effective mode is reported in `params`. Recorded in code, not just in docs: **the
  carrier keeps its place by DEFAULT, NOT BY EVIDENCE** — A and B are indistinguishable (1.035 vs 0.996,
  bootstrap 0.749) on the only dataset asked, and the bypass does not win parsimony either. `uptake="carrier"`
  is **bit-identical** to the shipped solve (pinned by a test), so `parameters.json`, `simulate()` defaults and
  `reproduce_demo` (0.029, re-verified) are UNCHANGED.
- **B1 — the DOSE SERIES: can concentration separate carrier from bypass? (this session) — the pre-registered
  rule is NOT met; the carrier is DISFAVOURED, not refuted**: `validation/dose_series_carrier.py` (pre-registration
  `f705de3` committed BEFORE the results) + `tests/test_dose_series_carrier.py` (9). The two entry mechanisms tied
  on Yamazaki because Yamazaki is **one exposure level**; they differ exactly in the response to CONCENTRATION
  (carrier saturates ⇒ BCF falls with dose; bypass is linear ⇒ flat), and the **carrier is the model's ONLY
  nonlinearity in exposure** (pinned by a test), so Tang's 5 soil doses over **1000×** are a structural
  discriminator. **Result**: only **1 of 3** congeners clears the pre-registered bar (needed 2), so **the carrier
  is not refuted on the rule as written** — but the tally is misleading and the run says why. **(a) GenX is
  NON-INFORMATIVE**: `Cwo/Km`=465 at the LOWEST dose ⇒ saturated throughout, carrier and bypass predict within 7%
  (this is failure mode G1 for GenX alone — **the gate should have been PER CONGENER**, a flaw in this file's own
  pre-registration that running it exposed). **(b) PFOA is INCONCLUSIVE**: observed 1.33× IS below the 1.48×
  midpoint but bootstrap 0.671 < 0.90, and it is Kd-limited (at `f_oc` 0.01 × Koc/10 the carrier predicts 1.05×).
  **(c) PFOS is the WELL-CONDITIONED case and it REFUTES**: lowest pore water of the three, the only compound that
  actually crosses `Km` inside the series (`Cwo/Km` 1.7 → 1729) — carrier predicts a **6.28×** decline, measured is
  **1.17×**, bootstrap **1.000**, robust in **8 of 9** Kd combinations. **The endpoint test agrees and is Kd-independent**:
  the model confirms entry magnitude **divides out of TF** (dose-invariant to 3 dp in BOTH arms) so a carrier cannot
  produce a TF trend — yet PFOA's TF falls **2.1–2.3×** while its BCF falls only 1.33×, i.e. the dose response sits in
  **translocation** (the documented toxicity), *larger than the whole uptake signal*, which is why (b) scores as
  unattributable. **THE DURABLE RESULT is a BOUND (§7, POST-HOC and labelled as such)**: to be as flat as measured the
  carrier needs **`Km` ≥ 500 µg/L — 100× the fitted 5** — and PFOA and PFOS land on the SAME bound despite ~6× different
  pore water. A carrier that linear across its whole exposure range **is** the bypass term (`Vmax/Km` as conductance)
  ⇒ the series does not so much choose between the mechanisms as **bound the carrier into the bypass's functional
  form**. Limits: 3 congeners, one soil, harvest only; SHAPES only (normalised to the lowest dose), never levels.
  **NOTHING ADOPTED** — `Km` is not re-fitted; `parameters.json`, `simulate()` defaults and `reproduce_demo` (0.029)
  UNCHANGED. `simulate` gained `km_scale` (default 1.0, bit-identical) as the dose knob.
- **B3 DONE — the η contradiction recorded in `docs/theory_anchor.tex`**: the derivation's parenthetical that η
  (which collects the apoplastic bypass ε and the carrier) is "essentially independent of tail length" is an
  ASSUMPTION of the `f_xy = η·φ_free` factorisation, and both halves of η are now measured to need a chain-length
  term (bypass corr **+0.832**/25× spread; LC6's carrier ~flat to C10 then 2.0×/5.5× at C11/C12). The equation is
  **left as written on purpose** — a requirement common to two DIFFERENT entry terms is unlikely to be a property of
  entry, so the natural reading is that eq. (factor) is **mis-partitioned** and what the fits absorb into η belongs
  in `φ_free`/`B_k` (where the two-pool work independently put it). No replacement factorisation is asserted; that
  is gated on the same wet-lab assay as the `k_seq` promotion decision.
- **CI now runs the WHOLE suite (this session)**: `.github/workflows/tests.yml`. Until now `rigor.yml` was the only
  workflow and it runs a SINGLE module (`test_sci_adk_rigor.py`), which is why every handoff carried "⚠️ CI does not test
  any of this — a green check means nothing here" and why the suite count in this file was maintained by hand from
  whatever a session happened to run locally. The two workflows stay SEPARATE on purpose (the rigor guard is cheap and
  its signal — an empirical claim was over-stated — would be buried in a ~22-min run). RDKit/phydrus come from
  `requirements.txt` and are NOT best-effort, since that file is what Streamlit Cloud installs; only emcee, the gfortran
  HYDRUS-1D build and sci-adk are, and `-rs` prints skip reasons so a silently shrinking suite stays visible.
- **App↔science parity audit + the four PRs the UI had missed (this session)**: an audit against the app found the last
  UI-touching commit was `f5e5593`, so **everything in PR #60–#63 was reachable only from the API** — and one of those,
  the K_PL-gated **lipid loading**, is the repo's single strongest cross-dataset result. Surfaced, all as opt-in named
  modes with the shipped model as the default (`parameters.json`, `simulate()` defaults and `reproduce_demo` 0.029
  UNCHANGED): (a) a sidebar **⚙️ Mechanism** expander carrying `uptake="carrier"|"bypass"` (with the "by default, NOT by
  evidence" verdict in its own help text), `lipid_loading`, and the `Vmax ×`/`Km ×`/`g_apo` entry-constant overrides —
  all threaded through `run_model`'s `sim_kw`, so they reach every tab, the congener comparison AND the Bayesian
  inverse; a **banner next to the headline metrics** whenever one is on, since the expander is collapsed and the run
  would otherwise read as the shipped one. (b) the **Tang tab's `lipid loading (OOS)` series** (`plots.fig_tang_tf`
  gained `val_extra`) — the honest counterpart to the refit bar: the green refit was calibrated ON those measurements,
  the purple lipid run was calibrated on Yamazaki and transferred untouched, and lands in the same place (0.516 vs
  0.519; PFOS stalk 0.013 → 0.620 vs Tang 0.571). (c) **weak electrolytes** in the neutral tab (`pKa`/`is_acid`/`pH`,
  printing `f_n` and which way the ion is pushed) plus **`g_apo`**, both carrying the tested BOUNDS in-UI (direction
  supported, magnitude refuted). (d) the **merged two-pool root + redistributed shoot** as a second BAF overlay
  (`simulate_twopool_nstem`). `estimate_exposure_bayesian` gained a `**sim_kw` pass-through so the inverse can never
  silently run a different model than the forward tabs. **Doc-truth fixes the audit turned up**, in the app text and
  the docs: the neutral tab still called Brunetti 2021's `K_RW`=13.3 an open question against the partition core
  (superseded — Kodešová 2019 measures the same compound at 1.10, so Brunetti is ~12× above a direct measurement and
  the deficit is confined to log Kow ≳ 3), still called the root lipid a soybean value (the ROOT is now corroborated by
  measured cereal roots; stem/leaf remain Trapp's soybean), said "five ways to drive the model" when there are six, and
  told users to `git submodule update --init` a source tree that is now VENDORED. Guards: `test_plots.py::
  test_fig_tang_tf_lipid_extra_series`, `test_model_api.py::test_estimate_exposure_bayesian_runs_under_the_selected_mechanism`.
  Verified in a headless Streamlit + Playwright drive of both views. `ui/` itself remains untested (no UI test layer).

- **Neutral organics promoted to a COMPOUND CLASS in the app (this session; user request)**: the neutral path was
  reachable only as an Expert-only *tab*, which made it a mini-app — no plant map, no exposure modes, no drivers, no
  inverse, no downloads. That was a UI-architecture mistake, not a model limit: "which compound" is a **scenario input**.
  The sidebar's `2 · Compound` now chooses the CLASS — curated congener · SMILES (both PFAS) · **Neutral organic
  (log Kow)** — and `run_model` dispatches to `simulate_neutral`, so every downstream view works on a neutral compound.
  This was cheap because `simulate_neutral` already returns the `simulate()` contract (`t/conc/baf/baf_final/straw/
  straw_baf/Cwo/Qtp/M/season/cwo_ref`); only `B_k` (→ `K_PW`) and the PFAS-only `params` keys needed guarding.
  What was ADDED to reach parity: (a) **`literature_params.koc_neutral`** — Karickhoff 1981 `log Koc = 0.989·log Kow
  − 0.346`, the soil sorption a neutral needs (the PFAS chain-length Koc QSPR is a per-CF2 group contribution and does
  not apply), wired into `cwo_profile_series(log_kow=|Koc=)`, `default_k_leach(Koc=)`, `hydrus_drivers(log_kow=|Kd=)`
  and `soil_hydrus.inputs_from_hydrus(Kd=)`, and exposed as an EDITABLE `soil Koc` field. Flagged **PROVISIONAL** on a
  specific ground: no table here scores a *predicted* Koc (the neutral tables either supply the exposure directly or
  carry the paper's own isotherm), and Li 2019's soil half shows the stakes — bias +0.033 with an experimental `K_om`
  vs +0.291 with an estimated one. (b) `simulate_neutral(cwo_profile=, cwo_kw=, Koc=)`, so the flooded/HYDRUS exposure
  shapes work there. (c) `estimate_exposure_bayesian(log_kow=, neutral_kw=)` — the inverse with the neutral forward
  model (recovers a known Cwᵒ to <1%). What is deliberately ABSENT for a neutral, rather than shown and ignored: `E_m`
  and the `f_xy` source (at z=0 a membrane potential has nothing to act on and `f_xy` IS the computed TSCF), the PFAS
  mechanism switches (the carrier exists to overcome anion exclusion; the lipid term is K_PL-gated on PFAS binding),
  and the four PFAS-only tabs (Yamazaki bars, chain length, congener compare, Tang TF) — dropped, not blanked. The
  headline swaps `eᴺ` (1 by construction) for **TSCF**, and a scope panel states in-UI: nothing fitted, a-priori 0.281
  (Liu)/0.783 (Ge) — 0.206 for Liu on the equilibrium basis — grain UNTESTED, stem/leaf lipid still soybean, γ=0 makes
  the leaf an upper bound (warning), and that this is a RICE model (another crop needs its own composition + drivers).
  The old `🧪 Neutral organics` tab is REMOVED (its content is now the class + the scope panel). PFAS defaults are
  untouched (`simulate("PFOA")` 0.478982/0.432189/0.147641; `default_k_leach("PFOA")` 0.05). Guards in
  `test_li2019_schriever_tables.py` (Karickhoff coefficients + that it is not the PFAS QSPR; the soil-side plumbing +
  a measured-Koc override; the neutral inverse). Verified end-to-end in a headless Streamlit + Playwright drive:
  parametric, flooded, custom tables, soil inventory and a REAL live HYDRUS-1D run for a neutral (Kd = 2.39 L/kg from
  its own Koc), plus the inverse.

- **Compound lookup from a name / CAS-RN / SMILES — EPA CompTox (CTX) (this session; user request)**:
  `src/chem_lookup.py` + `tests/test_chem_lookup.py` (8) + the neutral panel's **🔎 Look up the compound**
  box. The neutral path's inputs (`log Kow`, `MW`, `K_AW`, `pKa`, soil `Koc`) are properties of the COMPOUND,
  not of this model, and typing them by hand is where a user's error enters. A name/CAS/SMILES now resolves to a
  DTXSID (a SMILES via its InChIKey — the CTX chemical search has no structure route) and fills those fields.
  **The design constraint is provenance, not convenience**: CompTox serves EXPERIMENTAL and OPERA-PREDICTED
  values side by side, while every a-priori number this path advertises (Liu 0.206/0.281, Ge 0.783, Briggs stem
  0.299) is on a MEASURED log Kow — so a prediction silently standing in for a measurement would leave those
  numbers describing something the app no longer does. Hence: experimental WINS when both exist, every filled
  field carries an `experimental`/`predicted (OPERA)` badge, a predicted fill raises a warning naming the
  fields, and each field stays editable (a lookup seeds a DEFAULT, it does not commit the run). Two things are
  deliberately never auto-filled: the **in-planta half-life** (no dashboard property corresponds to it, and
  §4i measured the parent fraction varying 4.8× BETWEEN SPECIES for one compound, so it is not a compound
  constant) and any value whose UNIT is unrecognised — `henry_to_kaw` converts atm-m3/mol (`H/RT`) and Pa-m3/mol
  and otherwise returns None, because a wrong `K_AW` silently switches the leaf's volatilisation sink on or off.
  A returned pKa is the one model-PATH change a lookup can make (strictly neutral → weak electrolyte, which is
  TESTED-but-BOUNDED), so it is reported SEPARATELY from the provenance badge — `speciation_note` prints `f_n` at
  the root-zone pH and raises **two INDEPENDENT warnings**, a distinction the first live run made concrete:
  (i) **`f_n` below ≈0.1** is about the MODEL and fires however good the pKa is — **benzoic acid trips it on its own
  MEASURED pKa 4.18** (f_n 0.005), where §4l says the path is direction-supported but magnitude-REFUTED, so the run
  is a lower bound; (ii) **the pKa is PREDICTED** is about the SOURCE, and **carbamazepine is the in-repo
  counterexample**: OPERA returns acidic pKa **5.07** (f_n 0.036, 96% ionised at pH 6.5) while the MEASURED value
  this repo's own best-conditioned table is built on — Kodešová 2019 §4f, whose a-priori 0.191/0.237 assumes CAR is
  **un-ionised everywhere** — is **13.9** (f_n 1.00). Nine log units apart, on the very compound the neutral path is
  validated with, and only the measured one reproduces the published number. Hence the DEFAULT, decided with the user:
  **only a MEASURED pKa switches the model by itself** — a predicted one is filled into the ⚗️ panel and reported but
  the weak-electrolyte checkbox is LEFT OFF (the run stays strictly neutral until the user ticks it), because every
  other looked-up field is a value that may be wrong while a pKa decides WHICH MODEL RUNS. In code `neutral_kwargs()`
  still returns `pKa` (drop it to stay strictly neutral; the CLI note says so), and the tested-floor warning is raised
  in the ⚗️ panel where the path is actually on, so it covers a hand-typed pKa too. **Optional by
  construction**: needs a free EPA key (`CTX_API_KEY` env or the Streamlit secret `ctx_api_key`, never
  committed) and outbound network; with neither the panel says which is missing and everything stays manual.
  Cached a day per query; `CTX_BASE_URL` repoints the client (mirror or stub); 429/502/503/504 retry with
  backoff and anything else is fatal for that record. **CORRECTED against the live API (user-supplied working
  notes) — the first cut had the wrong host and the wrong provenance model**: the base URL is
  `https://comptox.epa.gov/ctx-api` (the older `api-ccte.epa.gov` NO LONGER RESOLVES; the offline tests could
  not catch that, so the constant is now pinned by one), and the property endpoints are SPLIT
  experimental/predicted with `/chemical/fate/…` carrying Koc. Five documented traps, every one of which fails
  SILENTLY (wrong or empty columns, never an exception), are handled and each has a test: (1) **provenance is
  the URL, not the payload** — a property record says nothing about being measured or modelled, so the bucket
  is tagged AT CALL TIME from the endpoint and never inferred from a field (the fate endpoint is the exception:
  it mixes both and carries `propType`); (2) `"NaN"`/`"null"`/`"N/A"` arrive as STRINGS and `float("NaN")`
  SUCCEEDS, so a blank measurement would rank first (experimental beats predicted) and blank out a usable
  prediction — string nulls are rejected before conversion; (3) field naming differs BETWEEN endpoints
  (camelCase `propName`/`propValue` vs snake_case `prop_name`/`prop_value` inside the fate records), so every
  reader takes an alias list; (4) units — Henry's law is atm-m3/mol and reading Pa-m3/mol as atm-m3/mol is a
  clean **5.006 log-unit** offset that still looks plausible, so the conversion is explicit and an unrecognised
  unit returns nothing — and the SHARPEST form of the same trap, **log vs linear**, is the one a LIVE run
  actually sprang: `LogKow` comes as `2.45 [Log10 unitless]` but **`Koc` comes as `549.541 [L/kg]`, LINEAR**,
  though Koc is normally written "log Koc", so assuming log10 is `10**549` (an OverflowError — the lucky case;
  a smaller value passes silently, and benzoic acid did reach `Koc=4.2e+31` before the fix). The scale is now
  read off each row's unit then its name (`is_log_scale`/`to_linear`/`to_log10`) and normalised at parse time,
  with an implausible Koc (>`KOC_MAX_LKG`) refused rather than handed to the soil model; (5) the **OPERA
  applicability domain** is populated only in the `…Global` fields — a
  prediction its own model puts OUTSIDE its domain is out of scope, not merely uncertain, so it loses the tie
  to an in-domain candidate, reaches the badge, and raises its own error in the panel. Parsing stays
  deliberately tolerant (name-substring matching, bare-list or `{"data": …}` envelope) so a schema change
  degrades to a MISSING property rather than a wrong one — and a **SIXTH trap the live run exposed that the
  notes did not list**: a property can come back under a name NONE of the patterns match (OPERA labels its two
  ionisation centres `pKa_a`/`pKa_b`, but a plain `pKa` row exists too — which is why the first live run
  returned no pKa for carbamazepine OR benzoic acid), and that shows as an EMPTY COLUMN, not an error. So a
  `pKa` row with no stated centre is now still used (acid by default, the acid/base radio left to the user),
  the reverse risk is closed too (`"log p"` is a substring of `"log pvap"` — a vapour-pressure row is now
  excluded from the log Kow slot), and `python src/chem_lookup.py <compound> --raw` prints an **INVENTORY of
  every returned property name/value/unit** before the JSON, which is the list to read when a field is blank
  (the JSON dump alone truncates before the interesting rows). The working notes' first rule (probe one
  chemical, read the raw output, then sanity-check on benzoic acid / benzyl alcohol) is in the CLI output.
  Tests are 100% offline fixtures (CI has no key and no network) written against the REAL shapes, pinning the
  base URL, the identifier classification, provenance-from-the-endpoint, that a string `"NaN"` cannot displace
  a prediction, the snake_case fate path, the 5.006-log Henry check, **the linear-Koc reading + the refusal of
  an implausible one + the vapour-pressure exclusion + the unlabelled-pKa fallback**, the AD tie-break, the
  `simulate_neutral` kwarg mapping (and that `half_life` is absent from it), and the no-key/no-network
  degradation. Verified end-to-end against a stub reproducing the LIVE shapes: experimental log Kow 2.45 chosen
  over the OPERA 3.39, the `"NaN"` row dropped rather than displacing it, `Koc 549.541 L/kg` read from the
  nested snake_case fate record as a LINEAR value, and the out-of-domain pKa labelled.

- **Sidebar widget ranges — a live crash, and the class it belongs to (this session)**: a user typed a
  compound into the CompTox box and the **whole app died** with a redacted
  `StreamlitValueBelowMinError`. Cause: **water** (log Kow −1.38) → Karickhoff `Koc` **0.019 L/kg**, below the
  soil-Koc field's `min_value=0.1`. Streamlit RAISES when a widget's `value` is outside its bounds, and
  `app.py` builds the sidebar at the TOP of the script, so one out-of-range value aborts the page **before
  anything renders** — the user cannot even edit the field that broke it. The bug was NOT the one input: the
  field's range has to cover what the panel ITSELF produces, and `koc_neutral` over the log Kow field's own
  span (−2 … 8) gives **0.0047 … 3.7e7**, so [0.1, 1e6] broke at **BOTH** ends — the upper one reachable by
  dragging log Kow to 8 with **no lookup at all**. FIX: the Koc field now spans `[0, KOC_MAX_LKG=1e9]`
  (matching the ceiling the lookup itself will return) with `format="%.4g"` so a small value reads `0.01946`
  instead of `0.0`, and a guarded step; and **every** looked-up seed (log Kow, MW, K_AW, pKa, Koc) goes
  through the new pure `ui.common.clamp_seed`, which clamps into range and **reports it** — never raising
  (which kills the page) and never silent (which would change a scientific input without saying so).
  `tests/test_ui_helpers.py` (4) is the repo's **first `ui/` test layer**, closing part of the documented "no
  UI test layer" gap: it pins the Koc range against `koc_neutral` over the whole log Kow span (the assertion
  that would have caught this), the lookup-ceiling match, `clamp_seed`'s contract, and — as a guard on the
  PATTERN rather than one field — that every seeded widget passes through the clamp. Verified in a headless
  drive: water renders `Koc 0.01946` with no exception, log Kow 8 → `3.681e+07`, −2 → `0.004742`, default 2.45
  → 119.4 unchanged. UI-only; `parameters.json`, the model math and `reproduce_demo` (0.029) UNCHANGED.
- **Sampling-protocol audit + the DRY→FRESH organ-mass fix (this session)**: `docs/sampling_protocols_audit_KR.md`
  (+ `docs/literature_db/Sampling_Protocols.csv`, `raw_si/rice_organ_biomass_literature.csv`). Question: which growth
  model, what measured organ masses exist, and how each dataset actually sampled/fractionated its compartments (grain
  milled or not, roots washed how). Read at source (main texts + SI; PDFs user-supplied, NOT committed).
  - **THE BUG (fixed)**: every biomass driver (`oryza_growth`, `growth_rice`, `measured_biomass` tables) returns DRY
    matter, but `PlantInputs.M` and `B_k` are FRESH weight; the dry matter went into `M` unconverted, so a terminal
    organ's concentration was a burden per kg DRY reported per kg fresh (leaf ×4.5, stem ×5.9, grain ×1.16; the root,
    set by uptake per unit mass, barely moves). **`src/organ_mass.py`** converts `fresh = dry/(1−θ_fw)` with the SAME
    θ_fw as `B_k`, at every point a driver becomes `M` (`model_api` incl. `simulate`/`simulate_nstem_leaf`/
    `simulate_neutral`/`tang_tf_validation`/`measured_forcing`/`drivers_from_arrays`/`hydrus_drivers`, `soil_hydrus`,
    `oryza_drivers`, `measured_biomass.biomass_drivers`, `M_for_nstem`, `pfas_rice_two_pool`, the app's growth-table
    seed, `neutral_dpu_validation.drivers`, `structural_adequacy_fit`). Named mode **`mass_basis="fresh"` (default) |
    `"dry"`** — `"dry"` is NOT a physical alternative, it only reproduces pre-fix artifacts (verified bit-for-bit on
    the documented numbers). The placeholder logistic (`measured_forcing=False`, `reproduce_demo` 0.029) is untouched.
  - **Effect** (dry → fresh): PFOA default BAF root/leaf/grain 0.479/2.52/0.148 → 0.488/0.614/0.036. **The unfitted
    neutral path improves sharply — independent evidence the fix is right: Ge 2017 0.783 → 0.311** (leaf over-prediction
    6–14× → 1.6–3.7×), **Liu 2023 0.281 → 0.213** (ODE mode; equilibrium 0.206 unchanged). **The Ge "~7-d in-planta
    half-life" was absorbing this error**: on the fresh basis shortening the half-life only hurts (7 d 0.773, 30 d 0.231).
    PFAS: Tang OOS free-anion 1.232 → 1.483, lipid 0.516 → 0.561 (the lipid generalisation survives), Tang refit 0.519 →
    0.844 (dry-fit f_xy); Kim grain OOS lipid 0.48 → 0.55 (DF≥15% 0.20 → 0.81). **Two recorded neutral findings were
    mostly this bug**: the "ODE-vs-equilibrium root discount" (ODE/K_PW 0.55 at logKow 1.78 → **0.92** on fresh mass — a
    10×-too-light root could not keep up with its own xylem export), and Kodešová's "metabolism cannot close the leaf/root
    gap" (no-metabolism 181 → **35**, 1.5-d half-life 28.5 → **4.5** vs measured 3.25 — on fresh mass a short half-life
    nearly closes it). The carrier two-pool PFDoDA closure runs its free f_xy to the bound (1) on fresh mass, so its
    "low f_xy + enhanced carrier" reading is dry-basis; `validation/twopool_nstem_merge.py` is pinned to `MASS_BASIS="dry"`
    (fit and Tang transfer consistent). Tests now pin BOTH readings where a recorded finding moved
    (`test_neutral_dpu`, `test_li2019_schriever_tables`, `test_weak_electrolyte`, `test_nstem_leaf`, `test_two_pool`).
  - **`f_xy_source="oryza"` RE-FIT on fresh mass** (`validation/refit_oryza.py`, 1.7 min): reproduction RMSE 0.236 → **0.283**;
    f_xy rises ~4× (PFOA 0.083 → 0.359) and **PFBA/PFOS/PFDoDA hit f_xy=1**, short-chain grain stays high at the L_Ph
    floor ⇒ part of the single-straw 4-pool's "saturated reproduction" was enabled by the bug. `*_oryza` now = fresh fit;
    the pre-fix fit is kept as `*_oryza_dry` / `params/refit_oryza_dry.csv` (used when `mass_basis="dry"`;
    `build_parameters.py` merges both). The remaining dry-basis fits (`TANG_REFIT_FXY`, `LIPID_LOADING`, the two-pool
    caches, scripts that build `M` themselves) were re-fitted in the NEXT session — see the bullet below.
    `tests/test_organ_mass.py` (9).
  - **Dataset descriptions corrected (at source)**: **Yamazaki 2023 is an Andosol paddy FIELD** (NARO Tsukuba,
    Kasalath+Koshihikari, ambient irrigation-water/rain/PM exposure) — NOT greenhouse/pot/spiked/"clean per-congener
    water"; it computes NO BAF (the repo's = tissue ÷ ONE irrigation-water grab), brown rice was NOT analysed (repo
    reconstructs 0.88·white+0.12·bran; the paper's own "brown" sums concentrations), straw = unweighted band mean, all θ
    assumed (whole plants oven-dried before segregation). Fixed in the app (`ui/expert.py`), README, docs, DB. **Kim 2019**:
    one plant per field, non-detects averaged as zero, PFOA IS site-paired and one WWTP field (A-2) is 83% of the pore-water
    mean (ratio of means 4.4 vs site-paired geomean 27), long chains on 1–4 detects, soils NOT Andosol (DB fixed). **Li
    2025**: authors' TF = mean of per-stage ratios, the repo's = ratio of 5-stage means (1.0–5.8× lower); water values are
    digitised irrigation-water shares; Methods say "leaves (including stems)". **Tang 2026**: endosperm most likely milled
    rice (model uses brown), chaff unused, roots rinse-only (surface confound applies), the old "HI_dw≈0.5" note was wrong
    (endosperm ≈10–16% of plant DW). **No dataset removed Fe/Mn plaque**; only Liu 2023 has a solvent desorption rinse.
    Growth anchors: Bouman DOI fixed (`…2004.09.011`), 1740 g/m² is Kumar 2017 not Bouman, `growth_rice` tables are Li 2017
    Table S2, Nada 2018 root:shoot is fresh-weight non-paddy (paddy DW: Kim 2001 FACE 0.026–0.033 at maturity).
    **No rice dataset reports organ water content over a season** — θ_fw remains the open input behind every dw↔fw step.
  - **Next-session handoff: `docs/HANDOFF_mass_basis_refit.md`** — in order: (1) re-fit the remaining dry-basis artifacts
    (`TANG_REFIT_FXY`, `LIPID_LOADING` — no committed fit script —, the two-pool fits) and re-check the OOS conclusions;
    (2) θ_fw provenance/stage dependence; (3) Li 2025 stage-resolved comparison; (4) Kim PFOA anchor + censored long
    chains; (5) align ORYZA parameters with IR72.D92 + validate; (6) Tang endosperm → white-rice composition, season 150.
- **Fresh-basis RE-FITS — handoff step 1 (this session)**: every fit that still ran on the dry-matter `M` was re-made on
  FRESH organ mass, with the dry one KEPT beside it and selected by `mass_basis="dry"` / `PFAS_MASS_BASIS=dry` (the repo's
  idiom; each dry record re-verified to reproduce its documented numbers). `tests/test_mass_basis_refits.py` (6) +
  `test_twopool_nstem_merge.py` (+2) + `test_two_pool.py` pin BOTH bases. `parameters.json`, `simulate()` defaults and
  `reproduce_demo` (0.029) UNCHANGED.
  - **`LIPID_LOADING` — the headline.** No fit script existed, so `validation/fit_lipid_loading.py` reconstructs the
    protocol (Yamazaki excl PFDoDA, measured forcing, `growth_rice`, per-congener W2 `kappa_d`/`L_Ph`); on dry mass the stored
    constants give 0.372, matching the documented ~0.36. **Fresh re-fit: 0.372 → 0.177 (root 0.488 → 0.088; PFUnDA root
    14.7 vs 19.5)** — so the "lipid loading drains the long-chain root" single-pool trade-off, the stated motivation of
    the two-pool root arc, was largely the dry-mass `M` (a root 10× too light). Caveat: the 6 globals sit on per-congener
    W2 `kappa_d`/`L_Ph`, so this is not a 6-DOF fit. `LIPID_LOADING` = fresh fit; `LIPID_LOADING_DRY` = the stored record;
    `lipid_loading_params(mass_basis)`. **OOS re-checked** (`oos_multidataset.py`): Tang per-organ lipid **0.520** (dry 0.516)
    vs free anion 1.483 — the generalisation SURVIVES; Kim grain excl-PFOA lipid **0.50** (dry 0.48) vs mono 2.77 / W2 1.67;
    but Kim's reliable DF≥15% subset worsens (0.20 → 0.74; that subset is PFNA ALONE) and Li straw/root flips to lipid (0.18 vs W2 0.81). The old line
    "lipid OOS ≈ in-sample Tang refit (0.516 vs 0.519)" NO LONGER HOLDS: the fresh in-sample refit reaches 0.327.
  - **`TANG_REFIT_FXY`** (`tang2026_fxy_refit.py`, 0.1 µg/g): fresh **PFOA 0.570 / PFOS 0.999 / GenX 0.089** (dry record
    0.097/0.320/0.017 kept as `TANG_REFIT_FXY_DRY`; `tang_refit_fxy(mass_basis)`). **PFOS sits at the f_xy = 1 ceiling** — the
    free anion cannot deliver Tang's PFOS shoot even unrestricted (Yamazaki's `f_xy_oryza` PFOS is also at 1 on fresh mass),
    so the "PFOS f_xy is condition-dependent, 0.14 vs 0.32" reading becomes "both at the ceiling". GenX's over-prediction by
    the provisional 0.233 shrinks from ~12× to ~2.6×. In-sample Tang refit 0.519 → **0.327**. NOTE the dry record itself
    predates later model changes: re-running the dry fit with the current code gives 0.114/0.376/0.020.
  - **Two-pool seq, measured forcings** (`twopool_root_measured.py`, now through `organ_mass`): in-sample 0.278 → **0.412**,
    Kim excl-PFOA 0.56 → 0.76. Root-matching needs **no sink at all for 5/11 congeners** (root already over-predicted at
    `k_seq`=0, PFOS among them), so the U-shape degenerates and **the non-K_PL PFOS/PFUnDA separation (4.5× dry) is
    UNDEFINED on fresh mass** (both sinks ~0). The script printed "HELD (30246×)" from a ratio of two near-zeros — its verdict
    now requires both sinks in use. Caches: `twopool_fitted_params_measured{,_dry}.json`.
  - **Two-pool + nstem merge** (`twopool_nstem_merge.py`, `PFAS_MASS_BASIS` env; fit ~10 min, not 40): in-sample 0.301 →
    **0.515**, Tang OOS 0.801 → **0.602** (still behind lipid 0.520). Same cause: the staged U-shape is fitted to root-matched
    `k_seq` that floor at 0, so a few points pull it (PFSA offset 10^+0.80 → PFBS root 18 vs 2.3). Caches
    `twopool_nstem_fitted_params{,_dry}.json`; the Tang baselines in its Stage 3 now run on the script's own basis.
  - **Carrier two-pool closure** (`close_longchain`): on fresh mass PFDA/PFUnDA/PFDoDA need f_xy 0.51/0.48/**1.0 (bound)**
    and carrier 1.26/1.27/**1.50×** (dry: 0.18/0.13/0.43, 1.41/1.51/3.15×). "Low f_xy + enhanced carrier" was a dry-mass
    reading; the long-chain carrier enhancement nearly flattens. PFDoDA straw still cannot close (21 vs 49.7 at f_xy=1) — the
    C12 shoot-loading ceiling remains. Pinned in `test_two_pool.py::test_model_api_hook`.
  - **Net reading** (in-sample except where marked OOS): on fresh mass the SINGLE-POOL lipid model (0.177) beats both
    two-pool variants (0.412 / 0.515), and the two-pool's distinguishing evidence (the non-K_PL separation) does not survive.
    This weakens the case for the two-pool root; it does not by itself retire it (B2 / the wet-lab gate are unchanged).
  - Scripts that build `M` themselves now go through `organ_mass` (`absolute_fxy_fit`, `bayesian_inverse_demo`,
    `hydrus_coupled_run`, `refit_realistic_biomass`, `root_shoot_biomass_sensitivity` — partition stats stay on dry matter,
    the placeholder row stays basis-free —, `twopool_longchain` and its `longchain_*` importers). With `PFAS_MASS_BASIS=dry`
    each builds the identical pre-fix `M`; they were converted, NOT re-run, so their recorded numbers are dry-basis records.
    The sci-adk run records (`runs/*`) are frozen dry-basis records too (the rigor guard verifies them, it does not
    recompute). App text (Tang tab, ⚙️ Mechanism help) updated to the fresh numbers.
- **Fresh-basis follow-ups — handoff 1b-①, step 6, step 2 (this session)**: `parameters.json`, `simulate()` defaults and
  `reproduce_demo` (0.029) UNCHANGED; every script below takes `PFAS_MASS_BASIS` and reproduces its dry record.
  - **Carrier vs bypass / dose series re-checked on fresh mass — the verdicts survive.** Both scripts had silently been
    running on fresh mass. Dose series: basis-free (PFOS 6.27× vs 1.17×, `Km` ≥ 500, 1/3). "An addition is necessary":
    holds (nothing 3.004, depolarisation 2.634). **"Which addition" holds only after a fairness fix**: as written the
    script re-fits the bypass's one parameter but keeps the carrier's W2-era `Vmax`, so on fresh mass the bypass "wins" at
    bootstrap 0.996. New POST-HOC section 6 of `carrier_vs_bypass.py` re-fits the carrier's global `Vmax` too: 1.161 vs
    1.194, P 0.30 (dry 0.16) — indistinguishable. Bypass chain-length independence refuted more strongly (+0.909). **The
    convergence is now shown on ONE setup**: per-congener carrier `Vmax` corr(n_C) +0.915 fresh / +0.880 dry vs the bypass's
    +0.909 / +0.832; where `f_xy` is free (carrier two-pool) the carrier correction vanishes ⇒ the correction lands on
    whichever term is free (B3, recorded in `theory_anchor.tex`). The documented record was the `--fast` grid.
  - **Tang endosperm** (`tang_tf_validation(endosperm="brown"|"white")`, default brown — milling unreported): the white-rice
    reading re-partitions the grain's free conc into a white composition with membrane PL counted as the `grain_brown` row
    counts it (the `grain_white` row's 0.003 is starch-complexed lyso-PL, a pool the brown row excludes; the source keys of
    both rows resolve nowhere). Factor PFOA 0.33 / PFOS 0.25 / GenX 0.66, and every Tang score WORSENS (lipid 0.520 → 0.700,
    refit 0.327 → 0.513): composition does not close the grain gap. `tang2026_fxy_TF_validation.py` season 120 → 150 d
    (4-pool Tang TF 1.90 → 1.74). Tang roots rinse-only → surface-sorption caveat recorded. `tests/test_tang_endosperm.py` (4).
  - **θ_fw sensitivity** (`validation/theta_fw_sensitivity.py`, PRE-REGISTERED and pushed before running): one organ at a
    time over `rice_tissue_params.csv` ranges, dry-side metrics. Root (0.88–0.92) Tang 0.039 / Ge **0.053 → INTERMEDIATE**;
    stem/leaf/grain SMALL (≤0.015). Shoot θ cancels on the dry basis (fresh-weight BAFs move up to 0.13 log — quote dry).
    Root θ is a VALUE uncertainty, not a time trend (flooded roots stay hydrated) ⇒ θ(t) is not worth an ODE change; a root
    water measurement is. θ's three copies (`model_api._COMP`, `organ_mass.THETA_FW`, `neutral_dpu.RICE_WATER`) now pinned
    equal. `tests/test_theta_fw_sensitivity.py` (4).
  - **DECIDED (with the user) — the xylem split runs on DRY matter, grain surface 0.5116 per kg dry.** The 4-pool core split
    the xylem between leaf and grain by `S·M`; since the organ-mass fix `M` is FRESH, so an organ's WATER set its transpiration
    share. Now `xylem_leaf_share()` (core) splits a fresh `M` on its dry matter `M·(1−θ)` with **`S_XYLEM_DRY` = leaf 20, grain
    0.5116 m²/kg dry** (`PlantInputs.M_basis`; drivers dicts carry a `mass_basis` tag, untagged = fresh; placeholder = None).
    Reading the original grain `S`=2 as per kg dry was tried first and REJECTED on evidence: it sends ~14% of the transpiration
    stream to the grain during filling (vs ~4% now; `nstem_leaf` assumes 5%), `refit_oryza` 0.283 → 0.394 with **8 of 11 L_Ph at
    the floor** (PFBA grain 11.65 vs 1.06 from the xylem alone), two-pool measured 0.412 → 0.513. **Key point: with constant θ
    the two bases are the same equation up to `S_grain`** — 0.5116 = 2×(1−0.78)/(1−0.14) is exactly the value at which the dry
    split reproduces the old fresh split, so **every number above is unchanged** (≤5e-5 rel.) and no re-fit was needed; what the
    decision buys is that grain θ no longer moves the split (θ 0.14 → 0.60: harvest share 0.118 both; Ge θ-effect 0.008 → 0.000).
    0.5116 is a fixed literal on purpose (computing it from θ would re-couple the split to water) and is NOT measured — a panicle/
    caryopsis transpiration share is the input that would replace it. `Compartment.S` (20/2) still serves the as-given split
    (placeholder, `mass_basis="dry"` records: bit-identical; `reproduce_demo` 0.029) and the `plant_air` surface. Pinned by
    `test_theta_fw_sensitivity.py` + `test_organ_mass.py::test_the_xylem_split_follows_the_mass_basis_of_M`.
  - **App fix found on the way**: the run cache (`ui/common._drivers_tuple`) rebuilt the driver dict from t/Cwo/Qtp/M only, so
    every driver mode silently DROPPED the ORYZA `leaf_loss` (leaf-senescence correction) — PFOA leaf BAF 1.08 instead of 0.55.
    Now every key round-trips (`tests/test_ui_helpers.py::test_the_run_cache_keeps_every_driver_key`).
- **Li 2025 stage by stage (handoff step 3; PRE-REGISTERED in `4df3047`, pushed before any model number) — Li now
  DISCRIMINATES, and the expectation was wrong**: `validation/li2025_stage_resolved.py` + `tests/test_li2025_stage_resolved.py`
  (5) + `docs/literature_db/raw_si/li2025_stage_conc.csv` (SI Table S6, all 15 PFAS × 5 stages × organ, ng/g dw; S5 water/soil
  by class in `li2025_stage_water_soil.csv`). The old `oos_crossdataset.li_tf` set the model's HARVEST-day BAF against the
  data's 5-STAGE MEAN on fresh weight; the authors' TF is a third estimator (mean of per-stage ratios). Transcription checked:
  the 5-stage means reproduce `data_obs/Li2025_BAF_TF.csv` (20/20 cells) and the authors' PFOA stem TF 2.35. Design: TF =
  tissue/root, **dry weight**, per stage (no water needed); stages have no dates, so ORYZA DVS maps them (tillering 0.4 /
  heading 1.0 / filling 1.25 / milk 1.5 / ripening 2.0 → days 28/66/78/90/116); the SAME estimator (stage geometric mean) over
  the SAME stages on both sides; grain scored from filling. **Results**: (Q1) the estimator alone moves the observed TF a
  median 3.0×, max 5.8× (the authors' mean of ratios is pulled up by one low-root stage: PFOA stem 2.35 of which the filling
  ratio is 9.67); (Q2) the stage-to-stage scatter is **4.5×** what the reported SDs imply ⇒ per-stage SHAPE is noise-limited,
  and (Q4) no variant shows a shape signal (ρ≈0); (Q3) **the redistributed shoot + lipid loading (`nstem_leaf`, Tang's
  winning configuration) discriminates: shoot/grain log10 RMSE 0.43/0.52** vs 4-pool lipid 0.81/1.52, monotone 1.25/2.42,
  W2 1.31/2.70 — beyond the pre-registered 0.10 margin on both, unchanged at DVS ±0.15. POST-HOC (not used for the verdict):
  the same winner under the other two estimators (0.66/0.66, 0.36/0.47); **both halves are needed** — `nstem_leaf` without
  lipid 1.11/1.20, 4-pool with lipid 0.81/1.52; the noise floor a perfect model would score is 0.31/0.32; TF is exposure-
  independent (Cwo 1 vs 0.05: ≤0.001 log). ⇒ the second independent per-organ dataset (after Tang) where the SAME combination
  wins out-of-sample. Limits: one field, roots not plaque-cleaned (biases obs TF LOW, so it cannot explain the models' under-
  prediction), assumed stage days, IR72 phenology for Jinyu Japonica No.22; `nstem_leaf` grain stays under (PFBS 0.30 vs 2.15,
  PFOA 0.33 vs 0.93). `oos_crossdataset.li_tf` is kept as the record (docstring points here); `runs/pfas-rice-oos-multidataset`
  stays a frozen record. `parameters.json`, `simulate()` defaults and `reproduce_demo` (0.029) UNCHANGED.
- **Kim 2019 with its non-detects (handoff step 4a; PRE-REGISTERED in `bbb3a8a`) — the lipid advantage is robust, the
  choice between lipid variants is not**: `validation/kim2019_censored.py` + `tests/test_kim2019_censored.py` (3) +
  `docs/literature_db/raw_si/kim2019_si_detection.csv` (SI Tables S5 MDL + S6 min/avg/max/DF; averages identical to
  `kim2019_field_conc.csv`) + `literature_params.KIM2019_DETECTION` / `kim2019_grain_baf_bounds()`. Kim averaged
  non-detects as ZERO, and for the long chains BOTH sides are mostly non-detects (PFDoDA rice 1/30, void water 3/27), so
  each published grain BAF is one corner of a range: with every non-detect in [0, MDL] the ratio of means lies in
  PFHpA [0.39, 2.11] · PFOA [4.43, 4.92] · PFNA [2.21, 7.07] · PFDA [1.55, 11.3] · **PFUnDA [19.8, 168] · PFDoDA [5.5, 630]**
  L/kg (the shipped ND=0 values are the LOWER corner for PFHpA/PFNA/PFOA). (Q1) The **C11 rise survives** (PFUnDA lo
  19.8 > PFDA hi 11.3); **the C12 level does not** (PFDoDA lo 5.5 < PFDA hi). (Q2) point / interval log10 RMSE (PFOA out):
  monotone 2.77/2.51, W2 1.67/1.52, **4-pool lipid 0.50/0.38, `nstem_leaf` + lipid 0.60/0.14** (3/5 inside the intervals
  each). Verdict: **no single winner** — both lipid variants beat the free-anion ones by >1 log on either score, but
  which lipid variant leads depends on the non-detect treatment (point margin exactly 0.10). The expectation held. The
  4-pool lipid's PFUnDA (7.8) sits below the robust lower bound; `nstem_leaf` + lipid (36.6) sits inside it. The
  `oos_crossdataset` label "[reliable: PFHpA, PFNA]" is fixed (PFNA alone). Step 4b (the PFOA anchor) is the next
  bullet. `parameters.json`, `simulate()` defaults and `reproduce_demo` (0.029) UNCHANGED.
- **Kim 2019 PFOA per site (handoff step 4b; PRE-REGISTERED in `bedadbf`) — the published 4.43 is one field's number;
  per site it is ~10, and the per-site ratio is not a plant constant**: `validation/kim2019_pfoa_anchor.py` +
  `tests/test_kim2019_pfoa_anchor.py` (4) + `docs/literature_db/raw_si/kim2019_table1_sites.csv` (main-text Table 1,
  per-site PFOA/PFOS/Σ17 in soil, void water and brown rice for the 30 fields; transcribed by the user from the PDF text
  layer and checked row by row against the rendered page; reproduces SI Table S6: rice 0.3494 ng/g 17/30, water 78.73
  ng/L 26/27, soil 0.1599 13/30) + `literature_params.KIM2019_PFOA_SITE_BAF`. (Q1) Candidates [L/kg]: published ratio
  of means **4.44** (rice over 30 sites, water over 27 — three E sites have no water), same 27 sites 4.18, without A-2
  21.5, per-site geomean of the 14 both-detected sites 27.0 (median 28.3). Pre-registered admissibility (same sites in
  numerator and denominator; no site dropped for detection status or value; no site above half the weight) leaves only
  the **censored per-site geomean: 9.94 L/kg, 90% site bootstrap [4.57, 20.16]** (26 sites: 14 exact + 12 rice<MDL as
  upper limits; site sd 0.82 log10; Kaplan-Meier median 7.77; leave-one-out 8.77–11.65). The ratio of means fails
  because A-2 (WWTP, 1760 ng/L) is **82.8%** of the void-water PFOA and has the LOWEST ratio (0.98); the "27" drops the
  12 rice-ND sites, which are low-ratio sites, so it is biased up. (Q2) **Rice ∝ water^0.52 [0.24, 0.87]** (censored
  slope; a constant BAF means 1) while (Q3b) the model is linear over Kim's range (slope 0.964 in every variant — only
  the carrier is nonlinear) ⇒ the harvest void water (a 2-day snapshot) is a noisy proxy for the season exposure, and the
  per-site ratio carries that noise; it is not the carrier. (Q3) Every variant sits BELOW the per-site interval:
  nstem_leaf+lipid 1.44 (−0.84 log), 4-pool lipid 0.34 (−1.46), monotone 0.036 (−2.44), W2 0.023 (−2.63) — the
  admissible estimator makes the PFOA miss LARGER than 4.43 did. (Q4) A-2 is **67.3%** of the void-water Σ17 too.
  Decisions: **D1** `kim2019_grain_baf()` keeps the ratio of means for every congener (the only estimator the other five
  have; mixed estimators would compare unlike quantities) and the per-site constant sits beside it; **D2** the
  `literature_params` demo's L_Ph target moves 4.43 → 9.94 — but by a 3% margin (E1 4.44 vs bootstrap low 4.57; robust
  over 10 seeds, NOT to the interval method: the POST-HOC profile-likelihood interval [4.18, 19.78] contains 4.43);
  **D3** the "PFOA\* used in the L_Ph fit" labels are removed (no shipped parameter is fitted to Kim — `L_Ph_W2fit`,
  `L_Ph_oryza`, `LIPID_LOADING` are Yamazaki fits) and PFOA is reported as its own OOS check, not folded into the
  pre-registered 4a RMSE; **D4** the absolute LEVEL of every Kim congener is A-2-weighted (the chain-length shape compares
  like with like; the level does not — for PFOA the per-site value is 2.2× the mean). Every pre-registered expectation
  held. Also fixed on the way: `validation_summary` panel B ("PFOA\* (used in fit)", "reliable PFHpA·PFNA") and
  `revalidation_crosscheck` ("DF≥15%" = PFNA alone) labels. `parameters.json`, `simulate()` defaults and
  `reproduce_demo` (0.029) UNCHANGED.
- **ORYZA2000 ported from the manual; the shipped "ORYZA" driver was a reduced re-coding (handoff step 5;
  PRE-REGISTERED in `288de63`)**: `src/oryza2000.py` + `validation/oryza2000_ir72ds92.py` + `tests/test_oryza2000.py` (6) +
  `docs/literature_db/raw_si/oryza_ir72_d92_params.csv` (IR72.D92, manual §7.3) / `oryza_ir72ds92_obs.csv` (the IR72DS92
  experiment, §7.2) / `irri_phil1_1992_printed.csv` / `nasa_power_losbanos_1992.csv`. Source: Bouman et al. 2001,
  *ORYZA2000: modeling lowland rice* (IRRI; free PDF) — every value checked by script against the manual and the
  parameters also against the APSIM `Oryza.xml` copy. **Scope chosen with the user ("누가 물어도 문제가 없는 방법")**:
  port the WHOLE potential-production path from the Fortran the manual prints — ORYZA1, PHENOL/SUBDD (transplanting
  shock), SGPCDT/SASTRO/SSKYC/SGPC1/SGPL/SRDPRF/GPPARGET (sunlit/shaded canopy, leaf-N profile, CO2), respiration, stem
  reserves, partitioning with the grain sink limit, SUBGRN (spikelet sterility), SUBLAI2 (seedbed → dilution → shock →
  SLA), leaf death, the carbon-balance check — and leave out only what potential mode switches off (water/N stress
  modules; their factors kept at 1 so the equations read as printed). One printed line is ambiguous: GPPARGET's AMAX is
  continuous between its branches only if REDFT·CO2 multiply the whole first term; that form is the default
  (`amax_form="as_printed"` keeps the literal text; they differ only below 20 °C daytime or away from 340 ppm). The
  carbon balance closes to ~1e-15. **The original engine/weather could not be obtained**: the WUR ORYZA2000 source
  (`O2000V2N13gnufortran.zip`) and the full PHIL1.992 station file are login-gated (the user's WUR account got 403; the
  IRRI v3 form is pending), so weather = NASA POWER corrected to the station by the 12 days the manual prints (Tmin
  +1.17, Tmax +1.67 °C, radiation ×1.211); temperature-only and uncorrected runs are reported too, and a station run
  will SUPERSEDE them. **Framing**: IR72DS92 is the manual's own calibration example, so for the port this is a
  FIDELITY test, not an independent validation of ORYZA2000. **Result** (scored main-field dates doy 34–114): the port
  is **ADEQUATE under every weather treatment** — PI/flowering/maturity doy 56/82/115 vs observed 58/83/114, log10 RMSE
  aboveground 0.024 (~6%), green leaf 0.085, stem 0.073, panicle 0.052; the shipped `oryza_growth` (anchored, as the
  PFAS default runs it) matures **12 d late** (126), grows too slowly early (WAGT 1373 vs 3498 kg/ha at doy 58), keeps
  growing its stem after flowering (7095 vs 4243 at maturity) and under-fills panicles (6912 vs 9843) — RMSE 0.23/0.29/
  0.17/0.24 ⇒ **the port fits better on 4/4 series and on maturity** (R2). The expectation held (except that the
  shipped maturity miss was 12 d, not 15–20). **Wired as an OPT-IN mode** `simulate(biomass="oryza2000")`
  (`model_api.BIOMASS_DRIVERS`; `_biomass_fn` now RAISES on an unknown name instead of silently falling back to
  `growth_rice`) and in the app's biomass radio, which now labels the default honestly ("ORYZA Level-1, reduced") next to
  "ORYZA2000 IR72.D92 (faithful port)". **POST-HOC (what a switch would move)**: in the PFAS climatology the port
  matures at t = 93 d; at the default 120-d season it raises straw ×1.6–1.9 and grain ×1.9–2.5 (root unchanged), but
  harvested at its own maturity straw is ×1.1–1.2 and grain ×0.6–0.75 — **the harvest date, not the crop model, moves
  the grain most**, and the 120-d default season was tied to the reduced driver's late maturity. **DEFAULT UNCHANGED
  (decision for the user)**: switching needs the ORYZA-based fits redone (`refit_oryza`, `TANG_REFIT_FXY`, the OOS reports
  that run on the default) and a season tied to crop maturity. `parameters.json`, `simulate()` defaults and
  `reproduce_demo` (0.029) UNCHANGED.
- **Independent check of the ORYZA2000 port — Shizukuishi FACE (PRE-REGISTERED in `9f90d55`; post-hoc addendum `4fca549`)
  — the IR72DS92 advantage does NOT transfer cleanly; DEFAULT UNCHANGED**: `validation/oryza2000_kim_face.py` +
  `tests/test_oryza2000_kim_face.py` (5) + `docs/literature_db/raw_si/nasa_power_shizukuishi_1998_2000.csv` /
  `jma_amedas_shizukuishi_1998_2000.csv`. IR72DS92 is the manual's own calibration example, so before any default switch
  the port and the shipped driver were run on an experiment neither was built on: Kim 2003 Table 2 (green/senescent leaf
  blades, exact), Kim 2003 Fig. 2 (totals, digitised) and Kim 2001 Table 1 (roots, a floor) — Akitakomachi, ambient CO2,
  medium N, 1998–2000. **Both drivers carry IR72 phenology**, so every observation is placed by STAGE on each driver's
  own DVS (the five samplings are tillering/PI/anthesis/mid-ripening/maturity — Kim 2003 abstract via OpenAlex) and only
  RATIOS are scored: Q1 green-leaf course, Q2 senescent share, Q3 leaf share of the plant, Q4 dry-matter timing, Q5 root
  share. **Pre-registered verdict: D1 FAILS on the primary weather** (uncorrected NASA POWER cell): the cell sits at 461 m
  vs a ~200 m paddy, so IR72 never reaches maturity by Nov 30 under EITHER driver, and the port's green-leaf course is
  worse (Q1 0.457 vs 0.205); **D2: it flips to PASS only under S1** (+1.7 °C lapse). The pre-registration wrongly said
  no station data exist — the JMA AMeDAS station Shizukuishi (195 m, ~8 km away) does, and runs **+1.1/+1.45/+1.2 °C**
  warmer than the cell (May–Oct), i.e. close to S1. **POST-HOC on station temperatures (S4)**: the port reaches maturity
  in all three years (anthesis DAT 91–102 vs observed 78–88, maturity 140–168 vs 120–131 — IR72 is not Akitakomachi), the
  shipped driver **never matures at 39.6 N in any treatment** (its grain-fill thermal time is ~1.7× IR72.D92's). The rule
  as written PASSES on S4 only because the shipped driver has no leaf death before DVS 1, so its senescent share at
  anthesis is exactly 0 and its Q2 is infinite; scored without that point (Q2x, post-hoc) **D1 fails on (b) in every
  treatment** — the port over-states the dead-leaf share (Q2x 0.328 vs 0.172; port 0.80–0.90 at maturity vs observed
  0.36–0.47, which is a floor since detached blades are lost). By group on S4 (port | shipped): Q1 0.259 | 0.241, Q3
  0.225 | 0.151 (shipped better), **Q4 0.104 | 0.166 (port better)**, Q5 0.403 | 0.358 (both ~2.4× over a root floor —
  not discriminating). Expectation: right that the port times early dry matter better and the shipped driver may not
  mature (it never did); WRONG that D1 passes, that both over-predict Q2, and that Q3 would tie. ⇒ on organ SHAPES the
  port is not better than the shipped driver on independent data; its clear advantages are fidelity to ORYZA2000 and
  finishing the season at a temperate site. Neither driver is a temperate japonica: a cultivar phenology calibration
  would be the next thing needed for either. `parameters.json`, `simulate()` defaults and `reproduce_demo` (0.029)
  UNCHANGED.
- **Option C — both drivers' phenology calibrated to Akitakomachi, then the FACE comparison re-run (PRE-REGISTERED in
  `7743bfc`) — the port now PASSES; the switch can be pre-registered (user decision)**:
  `validation/oryza2000_kim_face_phenology.py` + `tests/test_oryza2000_kim_face_phenology.py` (3). Only the four development
  rates of each driver are fitted, the DRATES way (ORYZA manual §8.2): one common factor on DVRJ/DVRI to the observed PI day
  (the juvenile/photoperiod split is not observable), DVRP to anthesis, DVRR to maturity, by root-finding on each driver's
  own heat units (hits the targets to <0.001 d; port carbon balance still ~1e-16). Partitioning, leaf death, SLA etc. stay
  IR72.D92 / shipped. **Leave-one-year-out** (held-out year runs on the other two years' mean rates), station weather,
  nursery 30 d (an ASSUMPTION — sowing is not in the transcription), Q2 scored without the anthesis point (Q2x).
  **C1 PASSES**: pooled 0.222 (port) vs 0.245 (shipped); Q1 green-leaf course **0.174 vs 0.287**, Q4 dry-matter timing
  **0.210 vs 0.307**, Q3 0.147 vs 0.138 (tie), Q5 0.324 vs 0.274 (both over a root floor), **Q2x 0.256 vs 0.221 — the
  port is still worse on the dead-leaf share, by 0.035 against the 0.05 allowance** (narrow). Holds under every
  sensitivity (nursery 20/40 d, POWER+1.7 °C) and with each year's OWN rates (0.246 vs 0.267). **C2**: the port's Q2x
  falls only 0.328 → 0.256 (−0.072 < 0.10) ⇒ its dead-leaf excess sits mostly in IR72's leaf-death table (DRLVT, from
  DVS 0.6 and steeper), which Kim cannot re-fit (PARAM needs stem/panicle weights). **C3 NOT adequate**: LOYO phenology
  misses anthesis by +10.5 / −10.3 / +0.7 d (1998/1999/2000) — the calibrated rates vary year to year (J/I ×1.09–1.35,
  DVRR ×0.86–1.21) — but the verdict holds with exact (own-year) phenology too. Calibrated rates move the port by
  ×0.86–1.35 and the shipped driver's DVRR by ×1.47–2.08 (its grain fill was far too slow). Expectation: right on C2
  and on no flips; WRONG that C1 fails on (b) and that LOYO phenology would be adequate. ⇒ once both develop like the
  cultivar, the port's organ shapes fit an independent experiment better overall; it still over-states dead leaves. The
  rates are an Akitakomachi-at-Shizukuishi result, not a temperate-japonica set. `parameters.json`, `simulate()` defaults
  and `reproduce_demo` (0.029) UNCHANGED. **Redone with the recorded nursery and the 1999 event dates — see the "Option C
  redone" bullet below: the PASS is robust to the nursery but SENSITIVE to the 1999 event dates.**
- **Second wave of rice-growth sources transcribed (this session; user-supplied PDFs, NOT committed)** — data only, no model
  run yet: `raw_si/yang2006_wuxi_face_dm.csv` (Wuxi FACE 2001–2003, Wuxiangjing 14: above-ground DM + leaf/stem/spike
  fractions at 5 stages — the stem/panicle partitioning Kim lacks), `gao2021_cold_japonica_oryza_params.csv` (an ORYZA2000
  crop set calibrated for a cold-region japonica, incl. DRLV/SLA/partitioning tables) + `gao2021_qingan_phenology.csv` +
  `gao2021_fig3_digitised.csv` (green stem/leaf/panicle + LAI time series, digitised by colour masks, then VERIFIED by
  overlay on the original figure — see the next bullet), `sasaki2005_shizukuishi_1999.csv` and
  `hasegawa2017_face_conditions.csv` (OA Table S1).
  **Two facts that bear on option C**: the Shizukuishi nursery was 20/27/29 d (sowing dates, Hasegawa S1) and covered
  (plastic chambers, Sasaki), not the assumed open 30 d; and Kim 2003's 1999 "PI"/"anthesis" samplings (DAT 62/88) were taken
  ~8/~10 d AFTER the events (PI ~54, heading 78 DAT; Sasaki) — so option C's targets for 1999 were late, and the LOYO anthesis
  "miss" of −10.3 d for 1999 lands on the real heading day. `tests/test_face_transcriptions.py` (6).
- **Gao 2021 Fig. 3 digitisation checked against the original figure (this session; user request "digitize하면 원본
  그림과 겹쳐서 확인")** — `docs/literature_db/tools/gao2021_fig3_check.py <the PDF>` re-fits the axes on every tick
  (≤13 kg/ha, ≤0.006 LAI off), checks marker centres against the error bars (an SD bar is symmetric, so its cap midpoint is
  the plotted value whatever the marker shape: 28 unambiguous pairs, +1.2 ± 0.8 px, max 2.7 ⇒ ±2 px = ±58 kg/ha / ±0.02 LAI
  for clear markers), and draws every value's full nominal marker outline on the original. **The first pass was biased
  wherever a marker was partly hidden** (by another marker, an error-bar cap, or the axis frame, which clips markers): it
  took the centre of the TRUNCATED visible box. Those markers are re-placed on a visible edge + the nominal half-size
  (15 values, e.g. 2019-06-16 DWST 332→199, DWGL 202→35; 2018-06-07 D-AGB 390→244; 2019-05-28 LAI 0→0.06), the two
  "hidden" markers are measured (2018-06-18 DWST 289, 2019-06-26 DWGL 390) and 8 partly visible markers that the colour
  masks missed are added; every row now carries `value_lo`/`value_hi`. Values below ~200 kg/ha (centre within 7 px of the
  axis) are flagged `near_axis` — NOT resolvable (2018-05-30 DWST 133 > D-AGB 71) and not for use in ratios. The organ
  sums now close at the early dates too, and the remaining mismatches (2018-07-26 +14.6%, 2018-08-18/09-17 −4.8/−3.5%) are
  15–29 px, i.e. in the figure, not the reading. Data not yet used by any model run. The earlier digitisations were
  checked next — see the next bullet.
- **The earlier digitisations checked against the originals (this session; user-supplied PDFs, NOT committed)** —
  `docs/literature_db/tools/fig_overlay_checks.py {kim2003|li2025|tang2026} <PDF>` re-fits each figure's axes on its own
  ticks, re-reads the values, prints them beside the stored ones and writes an overlay PNG. `tests/test_digitisation_overlay_checks.py` (3).
  - **Kim 2003 Fig. 2** (FACE total biomass; ambient open circles): corrections −50…+49 g m⁻², the **1999 HN DAT89 point was
    missing** (added), and the three near-zero tillering points were 20–27 g m⁻² too high (50/65/95 → 23/47/75; the 2000 value
    sat ABOVE its FACE disk, which the figure rules out). Hidden open circles are read from the crescent their ring leaves beyond
    the FACE disk. The values run +10…+25 g m⁻² (1–6%) above Kim 2001 Table 1 for the same crop (that table's year is itself
    inferred). **Re-run of 5b/5c on the corrected totals — no verdict changes**: 5b S4 Q3 0.225|0.151 → 0.224|0.148, Q4
    0.104|0.166 → 0.100|0.163; 5c Q3 0.147|0.138 → 0.148|0.133, Q4 0.210|0.307 → 0.207|0.306, pooled 0.222 vs 0.245 → 0.222 vs
    0.244, own-year 0.246 vs 0.267 → 0.246 vs 0.265 (the tillering totals are not scored).
  - **Li 2025 Fig. 2a** (irrigation-water composition): re-read shares close against the TABULATED SI S5 class sums within
    0.4 %-points on 4 of 5 stages and reproduce the printed PFOA 59.3–70.2 % / PFBA 18.6–31.1 % ranges. The water values the
    old `data_obs` Li tables were built with (external package, no script) are off for two congeners: **PFHxA 2.5 → 1.21 ng/L**
    (≈ PFHxA + PFHpA, two adjacent similar colours read as one) and **PFOS 0.4 → 0.71**; PFBA 17.9 → 18.7, PFBS 1.3 → 1.10,
    PFOA 46.7 → 46.5. New `raw_si/li2025_fig2a_water_shares.csv`; the verified water/BAF are added BESIDE the old values in
    `data_obs/Li2025_BAF_TF.csv` / `obs_baf_Li2025.csv` (`*_verified` columns; the old columns stay as the S6 record). The
    figure itself has an anomaly: its ripening bar carries a 4.45 % 6:2 FTS segment that S5 (precursors 0.05 %) rules out.
    TF = tissue/tissue does not use water, so the stage-resolved Li analysis (step 3) and the TF scores are unaffected. **POST-HOC
    consequence for S6** (`docs/H8_handoff_S6_final.md` §2, not re-run in place): root BAF ÷ the internal ceiling becomes PFHxA
    46 / PFBS 19 / PFBA 6.6 / PFOS 2.9 / PFOA 0.53 (was 22 / 16 / 6.9 / 5.1 / 0.53) — the excess orders by chain length, not by
    the reading quality, so the recorded reading "the Li surface excess is a water-denominator (reading) artifact" does not hold;
    what remains is that irrigation water is not the pore water the model needs.
  - **Tang 2026 Fig. 1a/1b**: correct — control biomass 33.47 g FW (stored 33.5), ear 6.86 (6.8), the same control bar in all
    three groups.
  - **Not checked**: Bouman & van Laar 2006 Fig. 2b needs no PDF — it is the IR72DS92 experiment whose exact table is already in
    `raw_si/oryza_ir72ds92_obs.csv`; against it the digitised rows are −12…+11% off (stated ±5%) and no script uses them. Kim 2001
    Fig. 1 (3 root values) is NOT checkable here: the session's `kim2001.pdf` turned out to be a JavaScript stub page, not the paper.
- **Option C redone with the recorded nursery and the 1999 event dates (POST-HOC addendum, PRE-REGISTERED in `b1774b4`) —
  the PASS survives the nursery but NOT the 1999 event dates**: `validation/oryza2000_kim_face_nursery_events.py` +
  `tests/test_oryza2000_kim_face_nursery_events.py` (3). Two option C inputs are contradicted by sources transcribed after it
  ran: the **nursery** (Kim 2003 Crop history, read at source: transplanted **14/23/24 d after EMERGENCE** from trays in plastic
  chambers ⇒ ORYZA SBDUR 14/23/24 + **TMPSB 9 °C**, the manual's Japan value — not an assumed 30-d open seedbed; Hasegawa's
  20/27/29 d count from SOWING) and, for **1999 only**, the **events** (Sasaki 2005: PI ~54, heading 78 DAT, vs the sampling
  targets 62/88). Everything else is option C (DRATES rates, LOYO, station weather, Q2x, C1). The file first re-runs option C
  through its own code (reproduced exactly); the shipped arm is untouched by the nursery (identical). **N (primary, nursery
  only): D1 PASSES**, pooled 0.220 vs 0.244 (option C 0.222); no flip under TMPSB 0 (0.237), nursery-from-sowing (0.206) or
  POWER+1.7 °C (0.220). The recorded nursery does NOT explain the year-to-year development differences: the port's J/I factor
  spreads further (×1.397/1.044/1.075 vs ×1.352/1.092/1.156) and the LOYO PI errors worsen (+14.6/−8.0/−6.4 d). D4: Q2x drop
  0.061 (< 0.10) ⇒ still the leaf-death table. **NE (1999 events; every sample placed at the DVS the driver's own-year run has
  on the sampling day — the 1999 PI/anthesis/mid-ripening samples sit at DVS ~0.76/1.26/1.66, not 0.65/1.0/1.53): FAILS on
  (b)** — port Q2x 0.281 vs shipped 0.181 (worse by 0.100 > 0.05); pooled 0.234 vs 0.214, shipped ahead. The placement control
  NG (same placement, option C targets) PASSES (0.218 vs 0.253) ⇒ the flip is the event dates, not the placement. Mechanism: the
  port OVER-predicts the dead-leaf share in every year (+0.10…+0.37 log), the shipped driver UNDER-predicts it at mid-ripening;
  reading 1999 at its true (later) phase narrows the shipped gap and widens the port's; moving anthesis earlier with maturity
  fixed at the sampling day also lengthens 1999's ripening (port DVRR ×1.209 → ×0.877). NE phenology: 1999 LOYO hits the true
  events (PI +0.0, heading −1.7 d); 1998 +10.6/+8.9/+13.5, 2000 −10.5/−6.6/−4.9 ⇒ C3 still NOT adequate. Expectation: right on
  N (factor spread widens, PI errors no better, pooled < 0.02, D1 PASS), D4, the NE phenology directions, NG ≈ N and no flips;
  **WRONG that NE passes**, and wrong in its reason for C3 ("E1 does not touch DVRR" — it moved 1999's). ⇒ **option C's PASS is
  robust to the nursery but SENSITIVE to the 1999 event dates: correcting the one year that can be corrected removes the margin
  on (b), and what decides it is the port's dead-leaf excess (IR72's DRLVT).** NE mixes event targets (1999) with sampling
  targets (1998/2000, no event dates exist) and keeps 1999's maturity at the sampling day, so it is not the better-informed
  verdict overall — but a default-switch pre-registration can no longer rest on option C's PASS alone. A cultivar leaf-death
  table (needs stem + panicle weights on the same dates) is what would settle (b). `parameters.json`, `simulate()` defaults and
  `reproduce_demo` (0.029) UNCHANGED.
- **A japonica crop set in the port — Gao 2021 (cold-region japonica) transferred to the FACE data (PRE-REGISTERED in
  `5aada97`) — it does NOT fix the dead-leaf excess; it makes it worse**: `validation/oryza2000_kim_face_gao_japonica.py` +
  `tests/test_oryza2000_kim_face_gao_japonica.py` (3). Gao et al. 2021 publish an ORYZA2000 crop set calibrated on Suijing-18
  (Heilongjiang, 47 N): partitioning FLVTB/FSTTB/FSOTB, leaf death DRLVT (from DVS 0.4, 0.024 at DVS 1, 0.020 after 1.6), an
  SLA TABLE (0.0010–0.0023 ha/kg, about half IR72's), a single KDF 0.4, RGRLMX 0.0088, FRPAR 0.45, and DVRs (COLDMIN not
  applied: the port has no cold-stop rule). Its calibration weather is unpublished, so this is a TRANSFER test. Arms A
  (IR72.D92), **G (Gao)**, B (shipped), all with phenology calibrated to Akitakomachi as in the redo (so they differ only in
  non-phenology parameters); A and B reproduce the redo exactly. Design N (primary) / NE (1999 events). **G1 NOT fixed**:
  Q2x G 0.331 vs A 0.267 (worse); **G2 FAILS on (b) under N AND NE** (Q2x 0.331/0.341 vs shipped 0.221/0.181) ⇒ by the
  pre-registered mapping "the Gao set does not transfer". But it is better on everything else: Q1 green-leaf course **0.110**
  vs A 0.173 / B 0.287, Q4 dry-matter timing **0.095** vs 0.183 / 0.306 (Q3 0.173 and Q5 0.415 worse); **G3: no general loss**
  (pooled 0.225 vs A 0.220). No flip under TMPSB 0 / nursery from sowing. **Decomposition** (one group swapped into IR72.D92,
  design N): the leaf-death table alone carries the Q2x rise (A+DRLVT 0.329) and the Q1 gain (0.127); the canopy group (SLA
  table, KDF, RGRLMX, FRPAR) leaves Q2x unchanged (0.267) and gives the best pooled score of any arm (**0.209**: Q3 0.109, Q4
  0.116); partitioning alone does little. **Gao's OWN development rates (no calibration)** miss anthesis by +9.0/−8.7/−2.1 d
  and maturity by −2.3/−10.2/−11.3 (IR72's own rates: +13…+14 d late at anthesis), and score pooled 0.201 (Q2x 0.299).
  Expectation: right on G1 (may rise), A+DRLVT ≈ G on Q2x, Q3 worse, G2 failing NE; WRONG on G3 (expected a general loss)
  and on Gao's own phenology being within 7/10 d. ⇒ **a japonica leaf-death table does not close the port's dead-leaf gap**
  — Gao's kills leaves earlier and the dead pool (WLVD) keeps them. POST-HOC caveat (not tested): Kim's senescent blades are a
  FLOOR (detached blades are lost) while ORYZA keeps every dead leaf, so part of every port arm's "excess" may be blades the
  measurement could not recover — at maturity the port's share (0.70–0.79 vs 0.36–0.47 observed) would need ~3/4 of the dead
  blades lost. Gao's set stays a candidate only; the port's default crop file is IR72.D92. `parameters.json`, `simulate()`
  defaults and `reproduce_demo` (0.029) UNCHANGED.

## 7. Build & run
- `pip install -r requirements.txt`
- **Main reproduction**: `python reproduce_demo.py` (Yamazaki BAF, W2 fit, RMSE≈0.029);
  `--rec` uses the monotone f_xy. Rebuild params: `python build_parameters.py`.
- **Visualization tool**: `pip install -r requirements-app.txt && streamlit run app.py`
  (plant/soil accumulation colormap + HYDRUS/soil/biomonitoring modes; the Expert sidebar's **⚙️ Mechanism**
  expander switches `uptake` carrier/bypass + `lipid_loading`; `2 · Compound` switches the COMPOUND CLASS —
  curated congener / SMILES / **neutral organic (log Kow)**; see `docs/visualization_tool.md`).
- **Compound lookup (EPA CompTox / CTX)**: `python src/chem_lookup.py carbamazepine` (or a CAS-RN,
  or a SMILES) prints log Kow / MW / K_AW / pKa / Koc **with their provenance**; `--raw` dumps every
  endpoint's payload when a field comes back empty. Base URL `https://comptox.epa.gov/ctx-api`
  (`api-ccte.epa.gov` is dead). Needs a free EPA key: `export CTX_API_KEY='…'` (CLI) or the
  Streamlit secret `ctx_api_key` in `.streamlit/secrets.toml` (app; keep it out of git). `CTX_BASE_URL`
  points at a mirror/stub. In the app: the neutral compound panel's **🔎 Look up the compound** box —
  experimental values beat predicted ones, everything is badged and editable, and the in-planta
  half-life is never auto-filled.
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
  **Structural merge** (two-pool root + redistributed shoot; the FAIR per-organ Tang OOS):
  `python validation/twopool_nstem_merge.py` (transfer → Yamazaki re-fit → Tang OOS; ~40 min on the MEASURED forcings;
  `--cached` reuses `twopool_nstem_fitted_params.json`, `--demo` reproduces the demo-forcing pathology). Opt-in API:
  `model_api.simulate_twopool_nstem("PFOA")` → the standard `simulate_nstem_leaf` dict + the root mobile/seq split.
  Measured-forcing robustness re-fit: `python validation/twopool_root_measured.py` (re-fits on forcing_rice + ORYZA
  biomass; in-sample + Kim OOS vs fxy-doc baselines; ~3 min). Opt-in API (no re-fit; reuses the cached fit):
  `model_api.simulate_twopool_seq("PFUnDA")` → the standard `simulate()` dict + root mobile/seq split.
- **Neutral organics (Briggs/Kow base)**: `python src/neutral_dpu.py` (per-compound TSCF/K_PW/BAF demo);
  `python validation/neutral_dpu_validation.py --obs data_obs/neutral_obs_liu2023.csv` (root partition, 14 compounds,
  RMSE 0.281) and `--obs data_obs/neutral_obs_ge2017.csv` (per-organ TF, RMSE 0.783); both also print the partition
  anchor + Kow signature + metabolism scope + switch and the half-life/TSCF sensitivity. Omit `--obs` for the
  structural checks alone — see `docs/neutral_dpu_validation.md`. In code:
  `model_api.simulate_neutral(2.45, name="carbamazepine", half_life=7.0)` → the standard `simulate()` dict
  + `K_PW`/`TSCF`/`rcf_briggs` (first arg is a **log Kow**, not a congener; `drivers=`/`biomass=` as usual).
  Soil-side for a neutral: `cwo_profile="flooded"` + `log_kow=`/`Koc=` (Koc from `literature_params.koc_neutral`,
  Karickhoff 1981, PROVISIONAL — prefer a measured value), `hydrus_drivers(log_kow=…)` or `Kd=` for the live engine,
  and `estimate_exposure_bayesian(log_kow=…, neutral_kw=…)` for the inverse. In the app: sidebar `2 · Compound` →
  **Neutral organic (log Kow)** (all six exposure modes, map/dynamics/drivers, inverse and downloads follow).
- **Weak electrolytes (acids/bases, not just strict neutrals)**: `model_api.simulate_neutral(2.45, pKa=4.5)` —
  `pKa=None` (default) is the strictly-neutral path and is bit-identical to before; a pKa routes through the
  `(f_n, f_d)` speciation pair instead (`is_acid=False` for a base, whose ion is a CATION and is ATTRACTED not
  excluded; `pH=` sets the root-zone split; add `phloem=True` for the leaf→sieve-tube pH ion trap). Helpers:
  `literature_params.speciation` / `ion_trap_factor` / `neutral_pathway_ratio`. BOUNDED by its own test (§6,
  `validation/weak_electrolyte_tscf.py`): direction SUPPORTED, magnitude REFUTED — use it for the direction of a
  speciation effect, not its size, and not below `f_n≈0.1`. Still no measured weak-electrolyte RICE dataset (the
  67-row test table is 16 species, none of them rice). In the app: the Neutral tab's ⚗️ expander (with `g_apo`).
- **Root-lipid anchor, both readings side by side**: `python validation/neutral_dpu_validation.py --lipid-source both`
  (all 5 shipped tables under `"measured"` vs `"briggs_anchor"`; add `--mode equilibrium` for the appropriate basis on
  the root tables, `--obs <table>` to restrict it to one). Single alternative run: `--lipid-source briggs_anchor --obs …`.
  In code: `model_api.simulate_neutral(3.72, lipid_source="briggs_anchor")` / `ND.rice_compartments(lipid_source=…)`.
  Default is `"measured"` — every published neutral number is on it.
- **Briggs 1983 stem anchor**: `python validation/briggs1983_stem.py` (~1 s; transcription self-check · the 4.1×
  coefficient gap · why it largely cancels in the observable). **The DATA behind those equations**:
  `python validation/briggs1983_shoot.py` (~1 s; reconstruction self-checks · the stem a-priori RMSE 0.299 ·
  why the stem/leaf split is NOT extractable · the article's own flagged compounds).
- **Li 2019 SOIL table (376 rows, the sign flip)**: `python validation/li2019_soil_table.py` (~5 s; hydroponic-vs-soil
  sign flip · per-crop-lipid collapse · the measured-vs-estimated `K_om` split · why `α_pt` is not applied).
- **Li 2019 root partition + the anchor diagnosis**: `python validation/li2019_rcf_apriori.py` (a-priori n=29 by
  species → the A1 lipid table → the anchor diagnosis → the full-ODE root-lipid scan; ~3 min, `--fast` skips the scan
  and prints the recorded numbers). Also runnable through the shared harness:
  `python validation/neutral_dpu_validation.py --obs data_obs/neutral_obs_li2019_rcf.csv` (scores `subset=apriori`
  only — the 18 Briggs barley rows are held out; pass `subset=None` in code to score everything).
- **TSCF QSPR on its own (no plant model)**: `python validation/schriever2020_tscf.py` (97 measured TSCF values;
  Briggs bell RMSE 0.310 / bias −0.221 on the 30 un-ionised rows vs the fitted Schriever refit's in-sample 0.234).
- **Weak-electrolyte path, tested (the OTHER 67 rows of the same table)**:
  `python validation/weak_electrolyte_tscf.py` (~35 min — §5 scans `g_apo` over 10 points × 97 rows;
  `--fast` skips every ODE solve and prints §1–2 only, ~5 s).
  Direction SUPPORTED / magnitude REFUTED — Spearman(f_n, TSCF) +0.480, speciation ON lifts the model's
  rank +0.284→+0.520 but under-delivers (bias +0.023→−0.203); bootstrap says the rank gain is robust
  (94%) and the RMSE loss is not (82%). `f_n` comes from the table's own logD, NOT its pKa — the two
  columns disagree by a median 1.41 log, and the pKa route manufactures a false counterexample out of
  the 8 pKa-1.62 rows (pinned by a test).
- **Is the carrier necessary? (the repo's one addition to Trapp)**: `python validation/carrier_vs_bypass.py`
  (~12 min; `--fast` coarser grids, ~4 min). Read the header (the PRE-REGISTRATION) before the VERDICT.
  An addition IS necessary (depolarisation REFUTED) but WHICH is undecided (bootstrap 0.749). The alternative
  is runnable, not buried: **`model_api.simulate(uptake="bypass")`** (default `"carrier"` = the shipped model,
  bit-identical) — see `model_api.UPTAKE_MODES`; explicit `vmax_scale=`/`g_apo=` still override the mode.
- **Dose series — can concentration separate them? (B1)**: `python validation/dose_series_carrier.py`
  (~6 min; `--fast` skips the Kd sweep, ~2 min). Read the header (the PRE-REGISTRATION) before the VERDICT.
  Pre-registered rule NOT met (1/3, needed 2/3) ⇒ carrier **DISFAVOURED, not refuted**: GenX is non-informative
  (saturated at every dose), PFOA inconclusive (bootstrap 0.671, Kd-limited), **PFOS well-conditioned and it
  refutes** (carrier 6.28× vs measured 1.17×, bootstrap 1.000, 8/9 Kd combos). Durable result = the POST-HOC
  bound: the carrier must have **`Km` ≥ 500 µg/L (100× the fitted 5)** to be this flat — i.e. linear over the
  whole span, which is the bypass's own form. `simulate(km_scale=…)` is the knob.
- **Kodešová 2019 carbamazepine (queue A4; the anchor vote)**: `python validation/kodesova2019_carbamazepine.py`
  (exposure from the paper's own measured isotherms + the `Koc` defence of the unit reading; a-priori 0.191;
  the shipped-vs-anchored table across all three root datasets; the Brunetti verdict; the dw→fw sensitivity).
- **Hwang 2017 (lettuce/chlorpyrifos; neutral, handoff A3)**: `python validation/hwang2017_lettuce.py`
  (exposure from the measured `Kd`; internal-consistency check on Table 1; the fw/dw basis span vs the
  `K_PW` ceiling; soil-contact bound; half-life/growth/water sensitivity). Read the VERDICT block — the
  RMSEs are a diagnosis, not a validation score. Not wired into `--obs` on purpose (wrong species/drivers).
- **Plant-air exchange** (volatilisation + gaseous uptake; opt-in, neutral path): `python src/plant_air.py`
  (permeabilities + volatilisation half-life across a volatility ladder). In code:
  `neutral_dpu.simulate_neutral(cmpd, drivers, air=True, air_kw=dict(C_air=…, S=…))` — needs `NeutralCompound(MW=…,
  K_AW=…)`. `K_AW=0` (PFAS) is identically zero and the core default `air=None` skips the term entirely, so no PFAS
  number moves. The §3b ladder in `validation/neutral_dpu_validation.py` reports it against the metabolism ladder.
- Tang 2026 f_xy: `python validation/tang2026_fxy_TF_validation.py` (4-pool TF vs Tang, ORYZA-driven);
  `python validation/tang2026_fxy_refit.py` (nstem_leaf + ORYZA f_xy re-calibration; 0.1 µg/g dose primary).
- **Organ-mass basis of a re-fit**: every re-fit script runs on FRESH organ mass and takes `PFAS_MASS_BASIS=dry` to
  reproduce its pre-fix record — `tang2026_fxy_refit.py`, `fit_lipid_loading.py` (the `LIPID_LOADING` constants;
  `--eval` scores the stored ones; ~2 min), `twopool_root_measured.py`, `twopool_nstem_merge.py`, `refit_oryza.py`,
  `structural_adequacy_fit.py`, the OOS reports `oos_tang_lipid.py` / `oos_crossdataset.py` / `oos_multidataset.py`,
  and `carrier_vs_bypass.py` / `dose_series_carrier.py` (their records are the dry `--fast` / default runs).
- **θ_fw sensitivity**: `python validation/theta_fw_sensitivity.py` (~5 min; pre-registered rule in the header; ends with
  the POST-HOC xylem-split note).
- **Li 2025 stage by stage**: `python validation/li2025_stage_resolved.py` (~25 s; pre-registered header; Q1 estimators ·
  Q2 noise · Q3/Q4 per variant · DVS sensitivity · POST-HOC estimators/decomposition/noise floor).
- **Kim 2019 with its non-detects**: `python validation/kim2019_censored.py` (~20 s; pre-registered header; the
  [0, MDL] bounds per congener · whether the long-chain rise survives · point vs interval score per variant). Tang endosperm as milled rice: `model_api.tang_tf_validation(c, endosperm="white")`.
- **Kim 2019 PFOA per site**: `python validation/kim2019_pfoa_anchor.py` (~1.5 min; pre-registered header; candidate
  estimators · censored per-site geomean + bootstrap · censored slope · model linearity · A-2 weight · decisions ·
  POST-HOC interval-method check).
- **ORYZA2000 port vs the shipped driver on IR72DS92**: `python validation/oryza2000_ir72ds92.py` (~1 min; pre-registered
  header; three weather treatments · phenology + per-organ log10 RMSE · R1/R2 verdict · POST-HOC PFAS BAF change). The port
  itself: `oryza2000.simulate(weather)` (daily ORYZA states from emergence) / `oryza2000.organ_biomass(t)` (kg/hill on days
  after transplanting); in the PFAS model `simulate(biomass="oryza2000")` (opt-in; pass `season=` ≈ the crop's maturity).
- **Independent check on the Shizukuishi FACE data**: `python validation/oryza2000_kim_face.py` (~5 s; pre-registered header
  + post-hoc addendum; stage-matched ratios under the pre-registered weather treatments, then the JMA-station run).
- **Option C (phenology calibrated to the cultivar)**: `python validation/oryza2000_kim_face_phenology.py` (~20 s;
  pre-registered header; DRATES-style rates per year · LOYO phenology error · stage-matched ratios · C1–C3 verdicts).
- **Option C with the recorded nursery + the 1999 event dates**: `python validation/oryza2000_kim_face_nursery_events.py`
  (~40 s; pre-registered header; N primary · NG placement control · NE 1999 events · sensitivities · D1–D5).
- **Gao 2021 japonica crop set on the FACE data**: `python validation/oryza2000_kim_face_gao_japonica.py` (~80 s;
  pre-registered header; arms A/G/B · Gao's own phenology · one-group decomposition · G1–G4).
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
- Tests: `pip install pytest && pytest` (**434 collected; 433 pass, 2 skip** (434 pass, 1 skip with `emcee` installed) — note the two numbers do not add up,
  and that is correct: `test_sci_adk_rigor.py` skips at MODULE level, so it contributes a skip OUTCOME while
  collecting zero tests, and the long-quoted "300 collected" was this same off-by-one. ~27 min with the full stack — RDKit + the built
  HYDRUS-1D engine + phydrus, as the SessionStart hook provides on the web; the `test_sci_adk_rigor.py`
  module additionally skips unless `sci-adk` is installed, which CI's `rigor.yml` provides). On a bare
  clone the structure/SMILES tests skip without RDKit and the HYDRUS-engine tests in `test_soil_hydrus.py`
  / the `cwo_profile='hydrus'` guards skip when the engine is unbuilt.
  **CI runs the WHOLE suite** (`.github/workflows/tests.yml`, added this session) alongside the
  narrow `rigor.yml` over-claim guard — until now `rigor.yml` was the only workflow and it runs a
  SINGLE module, so every handoff carried "a green check means nothing here" and the count above was
  maintained by hand. RDKit/phydrus come from `requirements.txt` and are NOT best-effort (that file is
  what Streamlit Cloud installs); only emcee, the gfortran HYDRUS build and sci-adk are, and `-rs`
  prints the skip reasons so a silently shrinking suite stays visible.
- FORTRAN (Method B): init submodule (`git submodule update --init`), then follow
  https://phydrus.readthedocs.io/en/latest/getting_started/compilation.html
  (gfortran + `makefile` / `make.bat`). NOTE: the top-level `makefile` lists the `.FOR` files
  without a path, so build from inside `source/` (copy the makefile in, as above).

## 8. Conventions
- Units: time **day**; aqueous conc **µg/L**; tissue conc **µg/kg**; mass **kg**;
  flow **L/day**; `B_k` in **L/kg fw** (`C_k = B_k · C_w,k`).
- **Organ mass `M` is FRESH weight** (like `B_k` and every observation). Every crop/biomass driver
  (`oryza_growth`, `growth_rice`, `measured_biomass` tables) produces **DRY** matter, so whatever turns a
  driver into `M` must go through `organ_mass` (`fresh = dry/(1−θ_fw)`, θ_fw = the binding factor's).
  `mass_basis="dry"` exists ONLY to reproduce artifacts fitted before that fix. A `drivers=` `M` (CSV,
  app growth table) is taken as fresh as-is. The ONE place that needs dry matter is the leaf/grain xylem
  split (`xylem_leaf_share`, `S_XYLEM_DRY` per kg dry), so `M` carries its basis: `PlantInputs.M_basis` and a
  `mass_basis` key on drivers dicts (builders set it; untagged = fresh; the placeholder logistic = None).
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
  `K_PL`-gated fit (`LIPID_LOADING`, fresh organ mass; `validation/fit_lipid_loading.py`): Yamazaki excl
  PFDoDA 0.177, root 0.088, and fixes long-chain grain. The root trade-off recorded on the pre-fix dry
  mass (0.36, root 0.49 — "single-pool limit") was largely the dry-matter `M`; `mass_basis="dry"`
  selects that record (`LIPID_LOADING_DRY`). EXPLORATORY / in-sample (OOS: Tang 0.52, Kim 0.50).
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
   (`kim2019_grain_baf()`); the demo fits `L_Ph` to the PFOA grain BAF (per-site 9.94 L/kg since step 4b;
   the ratio of means 4.43 is A-2-weighted — see §6).
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
