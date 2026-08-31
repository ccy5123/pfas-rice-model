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
│   ├── oryza_growth.py                       # MECHANISTIC ORYZA2000 Level-1 carbon balance → weather-responsive M_s(t) (opt-in; drivers=/weather=)
│   └── measured_biomass.py                   # ingest a MEASURED per-organ biomass table → M(t) driver (units→kg/hill; Tang etc.)
├── examples/                         # ready-to-load CSVs for app.py (HYDRUS drivers + biomonitoring)
├── params/                           # parameters.json (CANONICAL) + source CSVs (Bk, f_xy, Kcw, ...)
├── data_obs/                         # observed BAF/TF (Yamazaki, Li2025) + yamazaki_stem_height.csv
│                                     #  + NEUTRAL: liu2023/ge2017 (rice) · li2019_rcf (48 hydroponic RCF, 11 spp,
│                                     #    `subset` col holds Briggs' own rows out) · tscf_obs_schriever2020 (97 TSCF)
│                                     #  + biocides_bat_census (EU biocide INPUTS from the BAT report — NOT observations)
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
└── tests/                            # pytest (365 collected): plant, soil, hydrus, calibration, lit params, API (+two-pool, cwo_profile, k_leach), plots, structure(SMILES), oryza, measured-biomass, bayesian-inverse, NEUTRAL-organic (Briggs/Kow), twopool-nstem merge, BAT-census biocides

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
- **EU biocides from the BAT census run through the neutral path (this session) — PREDICTION, not validation**:
  `data_obs/biocides_bat_census.csv` + `validation/bat_census_biocides.py` + `tests/test_bat_census_biocides.py` (18)
  + `docs/bat_census_biocides.md`. A user-supplied report (`REPORT_bat_census.md`) records a regulatory screening tool
  (BAT) run over every EU biocidal active with a published bioaccumulation opinion, scored on a **fish** BCF. Different
  organism/endpoint/regulation, so exactly ONE thing transfers — an **audited log Kow of the UNCHARGED form**, which is
  the one input the neutral path needs (`K_PW` and `TSCF` are both functions of it, nothing fitted) and whose corruption
  is the report's own headline correction (a distribution ratio D entered where log Kow belongs; §7.5 — an error this
  model would inherit identically). **Census (UPDATED 2026-09-06)**: 44 substances named, 23 with a report-stated log Kow;
  the other 20 were left BLANK rather than filled in from elsewhere, and on request the BAT project then supplied **18 of
  them from its own collection sheet** (`EXPORT_logkow_full_20260906.csv`, 167 rows, with CAS + source + rank) ⇒ then the AUDITED
  provenance for all **61 substances BAT actually entered** (`EXPORT_logkow_entered_into_BAT_20260906.csv`) ⇒ then the
  **52 collection-sheet substances BAT never entered** that survive the file's own quality flags AND carry a pKa
  (`section=COLLECTION_SHEET_not_entered`, a SEPARATE tier — see below) ⇒ **117 rows run**: 22 log Kow printed in the
  report, 19 supplied on request, 21 from substances BAT entered but the report never named individually, 1 from this
  repo's Liu 2023 row, 2 second readings this study will not choose between, 52 collection sheet. **69 inside**
  the measured-data span, 10 root-only, 7 beyond every anchor, **31 EXCLUDED for ionisation** (the jump from 9 is mostly
  the sheet — 22 of its 52 are >90% ionised at pH 7, and the pKa column is exactly what lets the gate fire on them at all;
  without it they would have screened as neutrals, the failure the report records against its own §3.0a). A THIRD file then supplied the
  **BPC class + Scenario A/B fish BCF for all 61**, which is what makes the rank correlation an n=60 statement; the
  transcription this file already held was verified against it at **100 of 101 values** (the one difference: salicylic acid
  Scenario A, report-printed 4 vs export 4.35 — faithful transcription, more precise source). Their four caveats are kept
  per row in a `bat_caveat` column (15 empty Scenario B cells are the FACT — those substances are outside BAT's range so no
  kM was entered; tebuconazole/DCOIT have no opinion; cyphenothrin's nB came from its opinion's METABOLITE table; the
  triamine's BAT output is uninterpretable). **Tebuconazole's acid/base flag is RESOLVED and it caught an error here**: BAT
  entered it as an ACID at pKa 12.6205 and the base reading is a DIFFERENT pKa (3.516), so the two-row treatment that paired
  `is_acid=FALSE` with 12.6205 was the one wrong pairing (reads as 100% ionised); collapsed to one row, both readings >99.9%
  neutral so they cannot diverge. The audited export also replaced the four anticoagulant pKa this file had
  RECOVERED from the report's percentages with SOURCED rank-1 values (bromadiolone 4.5, brodifacoum 4.5, difenacoum 4.84,
  flocoumafen 4.95; three measured) — the recovered values were within 0.03 log, a check on the inversion, but sourced
  supersedes derived and nothing carries `pka_basis=derived` now. **Only creosote stays unrunnable** (no structure); the
  triamine now runs and the ionisation gate refuses it at 99.97% on a COMPUTED criterion instead of this file's assertion.
  Tebuconazole was briefly run BOTH ways (its export row gives `is_acid=TRUE` but `ionisation_class=base` with no
  percentage to arbitrate) until the project resolved it — see above; seven other rows disagree the same way but their
  percentage settles it (`is_acid` is the column that round-trips). Three checks before accepting the export: every stated pKa/%ionised pair round-trips at pH 7 to **<0.005 pp**;
  the pKa this file had RECOVERED from the report's percentages (coumatetralyl 4.781, warfarin 5.183) match the sourced
  values (4.75, 5.19) to 0.03/0.007 log — a check on the inversion, now replaced by the sourced ones; and provenance is
  uneven and recorded per row (IPBC is **rank 3 = a model prediction**; tebuconazole/DCOIT carry no rank). Cyphenothrin is
  deliberately left OPEN — their provenance says the AR's measured 5.79–6.09 supersedes the entered 6.29 and they declined
  to pick one (their §8.14 rule applied to themselves), so both ends are run and neither is called canonical.
  **The export's other 106 substances (never entered into BAT) — 52 of them now run, in a separate tier.** Round 1 read
  the first export as unusable as delivered (82 would fall inside on log Kow alone, but it carried zero pKa and zero
  source rank, so the ionisation gate cannot fire — the §3.0a failure — and a measurement is indistinguishable from a
  model output, the §7.9 trap). **SETTLED over three rounds, and the count that stood was the file's**: the first export genuinely had zero pKa and zero rank for the 106 (accurate about the file), but the inference "those two columns would make them equal to the 43" was WRONG — their sheet has 80 pKa and the real blocker is quality. Their covering note then gave 64 computed medians / 40 clean / 29 clean-with-pKa; counting the delivered `quality_flag` gives **28 / 10 / 68 / 52**. Both are right about different rules and the file shows which: **64 = rows whose `logkow_source_model` contains the STRING "median of N"; 28 = rows the flag marks "not among the candidates", a STRICT SUBSET (36 string matches are medians that ARE among their candidates)**. The string rule reproduces 64/40/29 to the unit, the criterion 28/68/52 — and the criterion is the report's own §8.14 test, the same over-detection that project had already corrected on its entered set (it is what replaced cypermethrin's median 6.175 with the sourced 6.3). They have CONFIRMED all four numbers and found a second defect their side: the flags were joined with `"; "` while the informational flag text contains a semicolon, so the documented split fragments every row — this count survived only because it matched substrings, which is luck and is recorded as such. **⇒ the runnable subset is 52, not 29**, and those 52 are now rows here (`section=COLLECTION_SHEET_not_entered`) marked as a LOWER tier: no source rank exists for any of them, no BPC class or fish BCF so they are excluded from the rank correlation, and their %ionised is computed here. Their pKa is what makes them runnable at all — **22 of the 52 are >90% ionised at pH 7**, so without it they would have been screened as neutrals, exactly the §3.0a failure. Allyl isothiocyanate's log Kow 34.675 is the median of `2.11 · 2.15 · 67.2 · 130.23` — two log values averaged with two LINEAR Kow values — tagged `experimental`; it is in both defect groups and never entered BAT so no result here is touched by it. The exclusion rule is the
  report's OWN §3.0a rule (">90% ionised at environmental pH") which it records as written-down-and-never-implemented;
  it is implemented here at the same threshold, and the pH-7 basis is pinned by round-tripping the six stated pKa to
  **0.002 pp**. **Result**: everything the opinions call bioaccumulative loads the ROOT and does not reach the grain
  (TSCF < 1e-3 for all of them). **The two scoreable substances** (a-priori, equilibrium basis): propiconazole vs the Liu
  2023 **rice** root RCF **+0.023 log (1.05×)**; triclosan vs the Li 2019 soil table (14 radish/carrot rows) +0.625 log —
  same sign as, and larger than, that table's documented +0.260 (§4h), part of it the speciation `K_PW` cannot see. Free
  cross-check: Li 2019's triclosan log Kow 4.8 vs the report's audited 4.76, **0.04 log apart**. **Cross-model**: (a) the
  report needed 217 runs to establish BAT is a one-input model; here it is STRUCTURAL and cyphenothrin/difethialone at
  log Kow 6.29 come back bit-identical; (b) the one property BAT proved inert — Henry — is the one that DOES reach this
  answer, through the opt-in air term (empenthrin K_AW 1.4e-2 loses ~47% of its leaf burden, cyphenothrin 9e-7 none);
  (c) the report's own DCPP pKa sweep re-run on identical inputs shows **this model damps ionisation HARDER than BAT**
  (at 91% ionised: chemistry 0.091, BAT 0.245, rice root **0.768**) because the membrane term sets the RATE of root
  equilibration and not its level — which is WHY the strongly-ionised group is excluded rather than caveated;
  (d) BAT's fish BCF peaks near log Kow 6.3 and turns over, this model's straw peaks at **1.75** and its root never
  turns over but **saturates kinetically** (K_PW/root 0.92 @4 → 0.001 @10.5: above ~7 the root number stops being a
  partition and becomes a rate); (e) Spearman(BAT fish BCF, rice root) **+0.875** vs (·, rice straw) **−0.769** on **60** substances (recomputed as the sample tripled: +0.467/−0.466 at n=20, +0.725/−0.725 at n=35, +0.875/−0.769 at n=60 — every enlargement STRENGTHENED it; two exclusions stated in code, the triamine's uninterpretable BAT output and the two second-log-Kow rows that would weight one substance twice) — a fish
  bioaccumulation class is NOT a statement about grain; (f) the report's largest finding is AMPLIFIED — its metabolic
  input moved fish BCF ≤82×, the same fish kM half-lives move this model's **grain ≤5,573×** (terminal accumulator),
  i.e. the least defensible input is the one the edible compartment is most sensitive to. **Honest limits**: no measured
  rice value exists for any of these substances; of the 7 the opinions call bioaccumulative, 4 are excluded for
  ionisation and only triclosan is inside the TSCF anchor — mirroring the report's §7.3 kM-QSAR Tanimoto 0.19–0.28
  finding from the other side (two independent models weakest exactly where the answer matters). The fish kM are run
  ONCE as a labelled sensitivity, never as a parameterisation; γ=0 makes leaf/grain upper bounds. `parameters.json`,
  `simulate()` defaults, `reproduce_demo` (0.029) and the neutral a-priori numbers (Liu 0.281 / Ge 0.783) UNCHANGED.

- **정책활용협의회 브리핑 (this session) — a slide-ready report, and an honest headline it refuses to make**:
  `docs/POLICY_BRIEF_KR.md` + `validation/policy_brief_runs.py` + `docs/policy_brief_results.csv` +
  `tests/test_policy_brief.py` (11) + two figures. The user is presenting to a Korean policy council and
  wanted the *material*, written so a slide-building agent can build the deck without re-deriving anything:
  every number carries a CSV key, every slide carries a `⚠ 반드시 남길 것` block, and Appendix A is a
  FORBIDDEN-PHRASE table ("모델이 검증되었다" → "한국 현장 자료 1건에 대한 외부 검증"; "기준치를 초과한다" →
  there IS no rice PFAS limit in either jurisdiction). **The headline the run refused to give**: the a-priori
  model says short chains reach the grain and long chains do not (C4–C6 mean 2.15 vs C9–C12 0.052, 41x), but
  BOTH measured datasets say the opposite at the long end — Yamazaki PFDoDA grain BAF **45.5** and Kim 2019
  Korean field **35.2** against the base model's 0.101 (**451x under**), and the base model's grain error runs
  monotonically from ~6x OVER at C5 to 451x UNDER at C12. So the congener ordering is a *mechanism* choice,
  not a result, and the brief says so in the slide, in the caveat block and in the forbidden-phrase table.
  **The evidence slide is Kim 2019** (the only Korean paddy set with paired pore water + brown rice, 6
  congeners): base model log10 RMSE **1.97** (bias −1.69) vs `lipid_loading=True` **0.53** (excl PFOA 2.06 →
  **0.48**) — and the lipid constants were fit on Yamazaki, so it is genuine OOS, reproducing the documented
  §multi-dataset figure exactly. **The inverse demo is deliberately two-sided**: synthetic truth 0.0787 →
  0.0788 µg/L (95% 0.053–0.117, the method works) but the Kim brown-rice measurement 0.349 µg/kg → **3.19
  µg/L against a measured 0.0787 = 40x over**, which is the same grain bias inverted — so the brief states
  the feature is NOT field-ready. Also run: PFOA end-to-end at the measured 78.7 ng/L (grain 0.0136 µg/kg,
  ±7.1x band, **5.4% of the EFSA group TWI**, xylem 78% / phloem 22% via `apportionment`), the saturable-uptake
  non-linearity (effective grain BAF 0.172 at 78.7 ng/L vs 0.148 at 1 µg/L — "농도 × 계수" is wrong), SMILES
  input flagging novels provisional, and the EU biocide census from the earlier PR. `parameters.json`,
  `simulate()` defaults and `reproduce_demo` (0.029) UNCHANGED — this session only reads the model.

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
- **Weak electrolytes (acids/bases, not just strict neutrals)**: `model_api.simulate_neutral(2.45, pKa=4.5)` —
  `pKa=None` (default) is the strictly-neutral path and is bit-identical to before; a pKa routes through the
  `(f_n, f_d)` speciation pair instead (`is_acid=False` for a base, whose ion is a CATION and is ATTRACTED not
  excluded; `pH=` sets the root-zone split; add `phloem=True` for the leaf→sieve-tube pH ion trap). Helpers:
  `literature_params.speciation` / `ion_trap_factor` / `neutral_pathway_ratio`. NOT validated against data —
  no measured weak-electrolyte rice dataset exists here.
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
- **정책 브리핑 수치 재현**: `python validation/policy_brief_runs.py` (~7분; `--fast` 는 그림 생략).
  `docs/POLICY_BRIEF_KR.md` 의 모든 수치 + `docs/policy_brief_results.csv` +
  `validation/figures/policy_grain_by_chain.png` · `policy_kim2019_grain.png`. Read the doc's
  "슬라이드를 만드는 분께" block first — it is written for a downstream slide-building agent, and its
  caveat blocks and forbidden-phrase table are load-bearing (guarded by `tests/test_policy_brief.py`).
- **EU biocides from the BAT census (PREDICTION, not validation)**: `python validation/bat_census_biocides.py`
  (~5 min; `--fast` skips the log Kow sweep, ~2 min). Read the VERDICT block: **117 rows run** (65 from the report + the
  audited entered set, + 52 from the collection sheet BAT never entered — a SEPARATE tier with no fish BCF to compare
  against), 69 inside the measured-data span, 10 root-only, 7 beyond every anchor and **31 excluded for ionisation** (the
  report's own never-implemented §3.0a rule, applied here). Only propiconazole (rice root, +0.023 log) and triclosan (soil radish/carrot, +0.625) can
  be scored at all. Writes `validation/bat_census_biocides.csv`; the input table is
  `data_obs/biocides_bat_census.csv` (NOT an `--obs` observation file — it has no measured plant value in it).
  Full record: `docs/bat_census_biocides.md`.
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
- Tests: `pip install pytest && pytest` (**365 collected; 364 pass, 2 skip** — note the two numbers do not add up,
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
