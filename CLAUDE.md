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
- Internal compartments exchange by **advection** (xylem up; phloem from the leaf in two streams:
  to the grain the sap its dry growth draws, to the root a fraction φ of transpiration) plus a
  **binding factor** `B_k = θ_fw + (1−θ_fw)·(f_prot·K_prot + f_PL·K_PL + f_cw·K_cw)`
  (basis A, fresh weight; **NO density prefactor**).
- Grain is **phloem-fed**; the weak-acid pH **ion-trap does NOT apply** (`f_n ≈ 0`) ⇒ phloem
  loading is carrier/channel (`L_Ph`), not a pH trap.
- **Grain and leaf are terminal accumulators**: the only sink is growth dilution, which → 0
  at maturity ⇒ no bounded steady state ⇒ final conc = time-integral / final mass. A
  **dynamic** model is therefore essential.
- Metabolism `γ_k ≈ 0` (PFAS recalcitrant). Air exchange off (`K_AW ≈ 0`).

Model report: `docs/pfas_rice_compartmental_model.tex` (PDFs are not in the repo — build with pdflatex)
Corrected neutral DPU base: `docs/dpu_model_summary_corrected.tex`

## 3. Repo layout
```
.
├── CLAUDE.md  README.md  requirements.txt  requirements-app.txt
├── reproduce_demo.py                 # entry point: Yamazaki 2020 site-year, dry basis (f_xy_source fit | --rec)
├── build_parameters.py               # (re)assembles params/parameters.json from source tables + refit_site_year.csv
├── app.py                            # Streamlit visualization tool (plant/soil map + exposure modes)
├── src/
│   ├── pfas_rice_plant_module_4pool_surf.py  # basis-A 4-compartment ODE (CANONICAL core; + K_surf)
│   ├── pfas_rice_plant_module_nstem.py       # N serial stem segments (multi-height mixer; Yamazaki gradient)
│   ├── pfas_rice_plant_module_nstem_leaf.py  # N stem segs + explicit leaf (transpiration deposition + retention)
│   ├── pfas_rice_plant_module.py             # import alias -> 4pool_surf
│   ├── harvest.py                            # observation operator: dead blades kept in leaf/straw; dry-basis values
│   ├── site_scenario.py                      # Yamazaki 2020 site-year drivers (JMA Tsukuba, ORYZA2000 port, ORYZA Tp)
│   ├── fao56.py                              # FAO-56 Penman-Monteith terms + ORYZA2000 Tp/Ep (eqs 4.2-4.4)
│   ├── oryza2000.py                          # ORYZA2000 potential production ported from the manual, IR72.D92 verbatim
│   ├── organ_mass.py                         # crop-model DRY matter -> the FRESH organ mass the ODE runs on
│   ├── neutral_dpu.py                        # NEUTRAL-organic (Briggs/Kow) path: same ODE with z=0
│   ├── plant_air.py                          # plant-air exchange (opt-in; 0 at K_AW=0)
│   ├── soil_paddy.py / soil_paddy_redox_corrected.py  # Freundlich soil -> C_w^o(t)
│   ├── soil_hydrus.py                        # REAL HYDRUS-1D run via phydrus -> Cwo(t), Qtp(t) (Method A)
│   ├── calibration.py                        # Tier-1 calibration machinery (synthetic recovery)
│   ├── literature_params.py                  # literature QSPRs/anchors (cited) + Kim 2019 tables
│   ├── model_api.py                          # UI-agnostic API: simulate(), scenario_drivers(), yamazaki_harvest_baf(), ...
│   ├── pfas_structure.py                     # SMILES -> Compound (RDKit; read-across + QSPR)
│   ├── plots.py                              # Plotly builders
│   ├── forcing_rice.py / growth_rice.py / oryza_growth.py  # generic-climatology drivers (scenario="generic")
│   └── measured_biomass.py                   # ingest a MEASURED per-organ biomass table -> M(t)
├── examples/                         # ready-to-load CSVs for app.py
├── params/                           # parameters.json (CANONICAL) + source CSVs + refit_site_year.csv
├── data_obs/                         # observed tables: Yamazaki (DRY basis), Li 2025, neutral tables
├── validation/                       # maintained analyses (refit_site_year, OOS, ORYZA checks, neutral, HYDRUS)
│   └── archive/                      # superseded analyses kept as the record (README says why)
├── docs/
│   ├── ASSUMPTIONS.md                # ★ every input: source, or unmeasured + the range it is run over
│   ├── OVERVIEW_KR.md                # 종합 진입점
│   ├── pfas_rice_compartmental_model.tex / dpu_model_summary_corrected.tex / theory_anchor.tex
│   ├── neutral_dpu_validation.md / visualization_tool.md / sampling_protocols_audit_KR.md
│   └── literature_db/                # curated parameter DB + raw_si/ extractions (weather, SI tables)
├── external/hydrus_source/           # VENDORED HYDRUS-1D 4.08 source (binary gitignored)
├── .claude/                          # SessionStart hook: web deps + HYDRUS engine build
└── tests/                            # pytest
```

## 4. Coupling strategy
- **Method A — loose, one-way (CURRENT; now WIRED to a real HYDRUS run).** HYDRUS-1D/Phydrus →
  `C_w^o(t)`, `Q_TP(t)`; the plant ODE is solved in Python (`src/pfas_rice_plant_module.py`).
  No FORTRAN edits. Interface = the three arrays in `PlantInputs` (`Cwo`, `Qtp`, `M`).
  `src/soil_hydrus.py` builds & runs the compiled HYDRUS-1D engine (via `phydrus`) for a paddy
  scenario per congener (Kd from the C3 `Koc` QSPR) and returns BOTH the pore-water trajectory
  `C_w^o(t)` and the root water uptake `Q_TP(t)` → `inputs_from_hydrus()` → `PlantInputs`. The
  soil run is driven by the scenario's transpiration (the site-year ORYZA2000 Tp by default,
  `forcing_rice` under `scenario="generic"`), so HYDRUS's actual uptake `vRoot` carries that crop
  shape (+ soil-water-stress feedback) and reproduces the driving transpiration when unstressed.
  See `validation/hydrus_coupled_run.py`.
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

Every number below is on the CURRENT model and was produced by the script named next to it.
Superseded analyses are kept in `validation/archive/` (its README says why each one was
retired); git history holds the older states of this file. Every input the model runs on,
with its source or its sensitivity range, is in **`docs/ASSUMPTIONS.md`**.

### 6.1 What the model is
- **Core**: basis-A 4-compartment ODE (`pfas_rice_plant_module_4pool_surf`), hybrid root uptake
  (GHK electrodiffusion + Michaelis–Menten carrier), xylem loading `f_xy`, mass-conserving phloem, grain formation gate, leaf senescence loss flux (ORYZA leaf death carries the
  leaf's PFAS out of the living leaf). `pfas_rice_plant_module` is an alias of it.
- **Organ mass is FRESH weight** (`src/organ_mass.py`): every crop driver produces DRY matter and
  is converted with the binding factor's θ_fw. There is no dry-mass mode.
- **Phloem = two streams from the leaf** (`phloem_flows`): to the grain, the sap its DRY growth
  draws (`grain_dry_growth` × `T_C_Ph`); to the root, φ·Q_TP. (Until amendment 2 of the refit the
  grain also received φ·Q_TP and the root φ times the grain flow — φ applied twice.)
- **Basis-consistent xylem**: the leaf/grain xylem split runs on dry matter
  (`xylem_leaf_share`, `S_XYLEM_DRY` = leaf 20, grain 0.5116 m² per kg dry; the grain value is
  unmeasured, ~6 % of the stream during filling on the site-year).
- **Harvest observation operator** (`src/harvest.py`): the ODE carries GREEN leaves; dead blades
  are reconstructed from the ORYZA leaf death rate (each carries its leaf's concentration at
  death) and kept in the leaf and straw (`dead_leaf_attached=1`; 0 = all lost, a sensitivity).
  Results carry `conc_dw` / `baf_final_dw` (dry basis) next to the fresh-weight values.
- **Default scenario `"yamazaki2020"`** (`src/site_scenario.py`): the calibration site-year —
  JMA Tsukuba 2020 weather, the ORYZA2000 port (IR72.D92) transplanted 2020-05-11, phenology
  calibrated DRATES-style to heading (harvest − 30.5 d) and harvest (= DVS 2), transpiration from
  ORYZA eq. 4.4 on FAO-56 terms. `simulate()` runs the single mean-date harvest (season 117.5 d);
  the Yamazaki observable (`yamazaki_harvest_baf`) combines the two harvest runs (09-01, 09-10)
  by geometric mean. **An explicit `season` other than the site's raises** (`check_season`);
  `scenario="generic"` (forcing_rice + a `biomass` driver over `season` days) is the alternative.
- **Observations**: Yamazaki on DRY basis (`data_obs/obs_baf_Yamazaki.csv`; `use` flags, the
  straw band range, PFBS grain unused — wholly < MQL).
- **Transport parameters** (`params/parameters.json`, from `params/refit_site_year.csv`):
  `f_xy_source="recommended"` (DEFAULT) = monotone physical `f_xy` + `L_Ph_rec`/`kappa_d_rec`
  fitted to root and grain given it; `"fit"` = (f_xy, L_Ph, κ_d) fitted per congener (saturated).
  A record without a fit (GenX, SMILES) takes the fitted congeners' geometric mean.
  `lipid_loading=True` uses `LIPID_LOADING` (six globals, fitted with the `*_rec` values).

### 6.2 Calibration on the site-year — `validation/refit_site_year.py` (pre-registered; amendments 1–2)
- **Fit A** (`"fit"`): log10 RMSE **0.162** over the 32 used cells (root 0.000, straw 0.035,
  grain 0.287) — a REPRODUCTION (3 parameters / 3 cells per congener). PFBA/PFPeA sit at the
  `L_Ph` floor with grain over (1.89 vs 1.23, 0.46 vs 0.31): their xylem alone over-delivers.
  PFDA/PFUnDA/PFOS sit at `L_Ph` = 1 with grain under (PFDA 2.7 vs 3.9, PFOS 1.7 vs 2.3), and
  PFDoDA at f_xy = L_Ph = 1 (straw 190 vs 249, grain 7.6 vs 52.9).
- **Fit B** (`"recommended"`): straw — a prediction — at **1.401** (short chains over, C8 and up under);
  grain 1.560: from PFOA (C8) up, and PFOS, the grain is out of reach at `L_Ph` = 1 (PFOA −0.20, PFDoDA
  −3.36 log), and PFBA/PFPeA sit at the `L_Ph` floor with grain over (+0.63/+0.60).
- **Lipid globals** (`LIPID_LOADING`): 0.207 over 29 cells (PFDoDA excluded).
- **Unmeasured inputs** (fitted parameters held fixed; yardstick = the two harvest runs' model
  difference, root/straw/grain 0.003/0.027/0.067): consequential are carrier Km ×100 (1.26
  everywhere; a re-fit sends κ_d → ~0 and f_xy ×3.9), E_m −90 mV (0.13 everywhere), dead blades lost
  (straw 0.174), the grain xylem surface (grain 0.28, straw 0.044), T_C,Ph (grain 0.160) and φ
  (straw 0.048); nursery, density, plants/hill, CO₂, heading date, harvest-not-at-maturity and the
  single-run shortcut are not. Full table: `docs/ASSUMPTIONS.md`.
- **Amendment 2** (phloem split into a grain stream and a root stream; pre-registered P1–P7): P1–P5
  and P7 held, **P6 wrong** (φ is consequential for straw, not grain); not foreseen: PFDA/PFUnDA/PFOS
  at the L_Ph bound with grain under, and T_C,Ph consequential for grain. Amendment 1: E1–E3 held,
  E4 wrong (a-priori straw above the expected 0.5–1.2), E5 partly wrong.
- `python reproduce_demo.py` prints Fit A (`--rec`: Fit B).

### 6.3 Transfer tests (no parameter fitted to the transfer data)
log10 RMSE; "default" = the configuration the app runs (`recommended`, free anion, 4-pool).

| dataset (script) | default | fit | lipid (4-pool) | redistributed shoot + lipid |
|---|---|---|---|---|
| Tang 2026 per-organ TF, 0.1 µg/g (`oos_tang_lipid.py`) | 1.714 | — | 0.992 | **0.640** |
| Kim 2019 grain BAF excl PFOA, point / [0,MDL] interval (`kim2019_censored.py`) | 2.33 / 2.05 | 0.73 / 0.55 | 0.56 / 0.38 | 0.47 / 0.15 |
| Li 2025 TF by stage, shoot / grain (`li2025_stage_resolved.py`) | 1.26 / 1.68 | 0.96 / 1.30 | 0.90 / 1.24 | **0.36 / 0.51** |

- **Lipid loading transfers; the free anion does not.** Lipid loading (fitted on Yamazaki only)
  cuts every transfer error, and with the redistributed shoot (`simulate_nstem_leaf`) it is best on
  all three datasets; it wins the Li stage test beyond the pre-registered margin on both endpoints,
  robust to the DVS mapping and the estimator. (`tang_tf_validation` runs the redistributed shoot
  unless `shoot="4pool"`; the free anion on the redistributed shoot scores 1.496 on Tang.)
- **Kim**: no robust winner between the two lipid variants (point margin 0.09 < 0.10); the C11
  grain rise survives the non-detects, the C12 level does not. **PFOA per site**
  (`kim2019_pfoa_anchor.py`): censored per-site geomean **9.94 L/kg** [4.57, 20.16] (the published
  ratio of means 4.44 is 83 % one WWTP field); every variant is below the interval (redistributed
  shoot + lipid −0.91 log); rice ∝ water^0.52, so the per-site ratio is not a plant constant.
- **Tang**: the in-sample Tang refit reaches 0.339 (`tang2026_fxy_refit.py`; f_xy PFOA 0.930,
  PFOS 0.999 — at the ceiling —, GenX 0.059; override only). Reading the endosperm as milled rice
  worsens the lipid score (0.640 → 0.831): composition does not close the grain gap.
- **Honesty band**: `model_api.APRIORI_OOS_SCORES` lists the default configuration's transfer
  errors (Yamazaki straw 1.40, Tang 1.71, Kim 2.33, Li 1.26 / 1.68); their median, log10 1.68
  (×48), is the band the Simple view shows. (Until 2026-10-25 the Tang entry was the redistributed
  shoot's 1.496, not the default's, which gave ×31.)
- Every transfer dataset runs on the calibration site-year drivers (Tang/Kim/Li sites and dates
  are not in the repo) — see `docs/ASSUMPTIONS.md` §5.

### 6.4 Mechanism and identifiability checks
- **Carrier vs bypass** (`carrier_vs_bypass.py`; pre-registered): nothing 2.301 · depolarisation
  −90 mV 1.890 · carrier as shipped 1.227 · carrier Vmax re-fitted (×8) 1.078 · bypass
  (`g_apo` 50) 1.093. An entry addition is necessary; which one is undetermined (bypass vs re-fitted
  carrier P = 0.368). Per congener both need the same chain-length correction (g_apo corr +0.898,
  carrier Vmax +0.916) ⇒ it belongs to binding/sequestration, not entry (recorded in
  `docs/theory_anchor.tex`). `simulate(uptake="bypass")` (g_apo 50) keeps the alternative runnable;
  the carrier is the default by convention, not by evidence.
- **Dose series** (`dose_series_carrier.py`; pre-registered): the rule is met by 1/3 congeners —
  PFOS refutes the fitted carrier (predicted 5.30× decline vs measured 1.17×, bootstrap 1.000),
  PFOA inconclusive (bootstrap 0.598), GenX non-informative ⇒ carrier disfavoured, not refuted. To
  be as flat as measured the carrier needs Km ≥ 150 (PFOA) – 500 (PFOS) µg/L, i.e. linear over the
  span.
- **Bayesian identifiability** (`bayesian_inverse_demo.py`): with transport fixed the exposure pair
  (Q_TP scale, Cwo level) recovers (corr +0.79, cond 18); Q_TP and f_xy trade off strongly
  (corr −0.97, cond 150 — only the intra-shoot advection that Q_TP alone sets keeps them apart);
  Cwo vs root conductance is a ridge (cond 1.6e6). `estimate_exposure_bayesian` (the app's inverse)
  estimates Cwo from DRY-weight tissue data with transport fixed.
- **θ_fw** (`theta_fw_sensitivity.py`; pre-registered): root INTERMEDIATE (Tang 0.047, Ge 0.053),
  stem/leaf/grain SMALL.
- **HYDRUS-1D coupling** (`hydrus_coupled_run.py`; Method A, the site-year transpiration drives the
  soil run): short chains leach (PFBA grain BAF 4.42 constant-Cwo → 0.54), long chains buffered.
  The analytic flooded shape's `k_leach` is calibrated per congener to the site-year HYDRUS run
  (`cwo_kleach_calibration.py` → `params/cwo_kleach.csv`).

### 6.5 Crop driver — the ORYZA2000 port
- `src/oryza2000.py` ports ORYZA2000 potential production from the manual (IR72.D92 verbatim;
  carbon balance closes to ~1e-15). On IR72DS92 (the manual's calibration example — a fidelity
  test) it fits all four organ series and the phenology; the reduced re-coding `oryza_growth`
  matures 12 d late and under-fills panicles (`oryza2000_ir72ds92.py`).
- Independent checks: Wuxi FACE organ partitioning (`oryza2000_wuxi_partitioning.py`) — port
  clearly closer on the leaf/stem/spike split; Shizukuishi FACE (`oryza2000_kim_face*.py`) — with
  cultivar phenology the port fits better overall but over-states the dead-leaf share (IR72's
  leaf-death table; a japonica table — Gao 2021 — does not fix it), and the verdict is sensitive
  to the 1999 event dates. The reduced driver never matures at temperate sites.
- The site-year scenario therefore runs the port; `oryza_growth` / `growth_rice` remain as
  `scenario="generic"` alternatives.

### 6.6 Neutral organics (Briggs/Kow path; nothing fitted)
`src/neutral_dpu.py` — the same ODE with z = 0, no carrier, K_PW = W + L·a·Kow^b, f_xy = TSCF(log Kow),
phloem off. The validation tables run on the generic drivers (other species and sites). Current a-priori scores
(log10 RMSE, `lipid_source="measured"`, the default):
Liu 2023 rice root 0.213 (equilibrium basis 0.206) · Ge 2017 per-organ TF 0.311 · Kodešová 2019
carbamazepine root 0.216 (equilibrium 0.237) · Li 2019 hydroponic 0.546 · Li 2019 soil (n = 376)
0.639 as an equilibrium partition (0.550 through the ODE) · Briggs 1983 barley stem 0.299 (Briggs' own fitted equation 0.282) · Schriever TSCF QSPR
alone 0.310. `--lipid-source both`: the Briggs-anchored root lipid is better on 2/5 tables (both
hydroponic) and worse on the rice and soil tables ⇒ default stays "measured". Li 2019's soil half
is unbiased where its exposure used a measured K_om (the apparent partition error is largely a
soil-side exposure error). Hwang 2017 is a diagnosis (fw/dw basis unstated), not a score.
Weak electrolytes (`pKa=`, `weak_electrolyte_tscf.py`, 67 ionisable Schriever rows): the
direction is supported (Spearman(f_n, TSCF) +0.480; speciation ON raises the model's rank
correlation +0.284 → +0.514 and lowers RMSE 0.333 → 0.279) but the scale is not (bias −0.16; the
median prediction is 7× below the median measurement); a small apoplastic bypass (`g_apo` ≈ 0.1)
improves both, the RMSE-optimal one (≈ 1) degrades the ordering — nothing adopted. Plant–air exchange is opt-in and exactly zero at
K_AW = 0. Details: `docs/neutral_dpu_validation.md`.

### 6.7 Data provenance
- `docs/literature_db/` (curated parameter DB + `raw_si/` transcriptions) and
  `docs/literature_db/Sampling_Protocols.csv` / `docs/sampling_protocols_audit_KR.md`: how each
  dataset sampled and fractionated (Yamazaki = an Andosol paddy FIELD under irrigation-water
  exposure; no dataset removed root plaque; brown rice is reconstructed for Yamazaki).
- Digitised values are checked against the original figures
  (`docs/literature_db/tools/fig_overlay_checks.py`, `gao2021_fig3_check.py`).
- Weather: `raw_si/jma_tsukuba_2020.csv` (site-year), NASA POWER / JMA AMeDAS series for the
  FACE checks. PDFs are not committed.

## 7. Build & run
- `pip install -r requirements.txt` (RDKit and phydrus included; `requirements-app.txt` for the app,
  `requirements-validation.txt` for emcee).
- **Reproduction**: `python reproduce_demo.py` (Fit A on the site-year, dry basis; `--rec` = Fit B).
  Re-fit: `python validation/refit_site_year.py` (~1 h; writes `refit_site_year.json`; `--write`
  updates `params/refit_site_year.csv` and `parameters.json`), then `python build_parameters.py`.
- **App**: `pip install -r requirements-app.txt && streamlit run app.py` (Simple = Korean, Expert =
  English; crop drivers site-year | generic; compound class congener | SMILES | neutral log Kow).
- **API**: `model_api.simulate("PFOA")` (site-year, mean harvest) · `yamazaki_harvest_baf(c)` (the
  two-run observable) · `simulate(..., scenario="generic", season=…, biomass=…)` · `simulate_nstem_leaf`
  · `simulate_neutral(log_kow, …)` · `tang_tf_validation(c)` · `estimate_exposure_bayesian(...)` ·
  `apportionment(c)`.
- **Transfer tests**: `validation/oos_multidataset.py` (summary), `oos_tang_lipid.py`,
  `kim2019_censored.py`, `kim2019_pfoa_anchor.py`, `li2025_stage_resolved.py`, `tang2026_fxy_refit.py`.
- **Mechanism checks**: `carrier_vs_bypass.py` (`--fast`), `dose_series_carrier.py` (`--fast`),
  `bayesian_inverse_demo.py` (`--emcee`), `theta_fw_sensitivity.py`.
- **Crop driver checks**: `oryza2000_ir72ds92.py`, `oryza2000_kim_face.py`,
  `oryza2000_kim_face_phenology.py`, `oryza2000_kim_face_nursery_events.py`,
  `oryza2000_kim_face_gao_japonica.py`, `oryza2000_wuxi_partitioning.py`.
- **Neutral path**: `validation/neutral_dpu_validation.py --obs data_obs/neutral_obs_liu2023.csv`
  (or `…ge2017.csv`, `--lipid-source both`, `--mode equilibrium`), `li2019_rcf_apriori.py`,
  `li2019_soil_table.py`, `kodesova2019_carbamazepine.py`, `briggs1983_stem.py`,
  `briggs1983_shoot.py`, `hwang2017_lettuce.py`, `schriever2020_tscf.py`,
  `weak_electrolyte_tscf.py` (`--fast`).
- **HYDRUS-1D** (vendored source): `cp external/hydrus_source/makefile external/hydrus_source/source/ &&
  (cd external/hydrus_source/source && make)`; `python src/soil_hydrus.py`;
  `python validation/hydrus_coupled_run.py`; `python validation/cwo_kleach_calibration.py`. On Claude
  Code on the web the SessionStart hook builds the engine.
- **Compound lookup** (EPA CompTox): `python src/chem_lookup.py carbamazepine` (needs `CTX_API_KEY`;
  base URL `https://comptox.epa.gov/ctx-api`).
- **Tests**: `pytest` (424 collected; all pass with the full stack, plus 1 module-level skip for `test_sci_adk_rigor.py` without sci-adk; ~17 min). CI (`.github/workflows/tests.yml`) runs the whole suite;
  `rigor.yml` runs the sci-adk over-claim guard on the frozen `runs/` records.

## 8. Conventions
- Units: time **day**; aqueous conc **µg/L**; tissue conc **µg/kg**; mass **kg**; flow **L/day**;
  `B_k` in **L/kg fw** (`C_k = B_k · C_w,k`).
- **Organ mass `M` is FRESH weight.** Crop drivers produce DRY matter; whatever turns a driver into
  `M` goes through `organ_mass` (`fresh = dry/(1−θ_fw)`, θ_fw = the binding factor's). A `drivers=`
  `M` is taken as fresh. Drivers dicts carry `mass_basis="fresh"`; the placeholder logistic carries
  None (split on `M` as given).
- **Observations are compared on DRY basis** through the harvest operator (`baf_final_dw`,
  `conc_dw`); fresh-weight values remain in `baf_final`.
- **Scenario**: `"yamazaki2020"` sets the season, the crop and the transpiration; `season=` applies
  only to `"generic"`. Exposure `C_w^o(t)` is constant by default (conc == BAF);
  `cwo_profile="flooded"|"hydrus"` gives a season-mean-normalised time shape.
- **Binding = basis A**: `B_k = θ_fw + (1−θ_fw)·(f_prot·K_prot + f_PL·K_PL + f_cw·K_cw)`; `f_*` dry-weight
  fractions, `K_*` L/kg pool-dw; no density prefactor.
- `f_xy` ∈ (0,1] is the root→xylem loading factor (TSCF analog). The free anion loads `f_xy·C_1/B_1`;
  lipid loading adds a B-independent `g·C` term (opt-in).
- Symbols map 1:1 to `docs/pfas_rice_compartmental_model.tex`.
- A new analysis that decides something is **pre-registered** (rule + expectations committed before
  the run) and its misses are recorded; an unmeasured input gets a stated default and a sensitivity
  range in `docs/ASSUMPTIONS.md`.

## 9. Next tasks / open questions
1. **Measurements that would replace assumptions** (`docs/ASSUMPTIONS.md` §7): Yamazaki 2023 SI
   S16/S18/S19 (raw values, MQL flags, planting density, nursery); a panicle transpiration share
   (grain xylem surface); rice phloem sap volume per grain DM (T_C,Ph); organ water content over a
   season; rice-root cell-wall / plaque PFAS sorption (K_cw, K_surf); carrier kinetics (Km, Vmax).
2. **Transfer datasets on their own drivers**: Tang (site, dates), Kim (field calendars), Li (stage
   dates) — each currently runs on the Yamazaki site-year.
3. **Grain**: the free-anion model under-predicts transferred grain; lipid loading closes part of it.
   On Yamazaki the long-chain grain is out of reach at L_Ph = 1 even with f_xy fitted (PFDoDA ~7×).
4. **Default configuration**: lipid loading + redistributed shoot transfers best on Li and Kim
   (lipid 4-pool on Tang); whether it should become the default is a decision for the user (it is
   in-sample on Yamazaki and the Kim winner is not robust).
5. **Method B** (tight HYDRUS coupling) — later.

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
