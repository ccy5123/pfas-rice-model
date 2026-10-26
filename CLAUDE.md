# CLAUDE.md — PFAS Rice Compartmental Uptake Model

> Persistent context for Claude Code. Read this first. Full math lives in `docs/`.

## 0. Working rules (from the user — apply in every session)
- **Full text not obtained → tell the user, as a DOI list.** Whenever a paper the work needs
  could not be read in full (paywall, download error, page refused, only an abstract or a
  search snippet seen), say so to the user right away and give a list: DOI · first author
  year · journal · what is needed from it (table/figure/value) · why it failed. The user
  obtains papers out of band. Never let an abstract- or snippet-only value pass as read:
  mark it as such in the work and in `sourced_inputs.csv`. Try the open-access copy first
  (OpenAlex `oa_url`, J-STAGE, publisher PDF) and list only what still fails.
  Running list: `docs/literature_db/fulltext_needed.md`.

## 1. Goal
Mechanistic **four-compartment dynamic model** for PFAS bioaccumulation in rice
(*Oryza sativa*), built as an **ionizable-organic-compound (IOC) extension** of the
Trapp/Brunetti **Dynamic Plant Uptake (DPU)** framework, designed to couple with
**HYDRUS-1D** for the soil side.

## 2. Scientific summary (see `docs/` for the full derivation)
- PFAS = **permanently dissociated anion** (very low pKa, `f_d ≈ 1`). The neutral-compound
  Briggs/Kow partition core does **not** apply.
- Compartments: `root(1), stem(2), leaf(3), fruit/grain(4)`.
- **Root uptake**: the theory (`docs/`) derives a hybrid entry, ionic electrodiffusion (GHK;
  inside-negative membrane ⇒ anion *exclusion*) **+** a saturable carrier (Michaelis–Menten). The
  **shipped model** uses a linear equilibrium entry `λ_eq·(r_eq·C_w^o − C_w,root)`: no source gives
  the hybrid's constants (E_m in a paddy root, carrier Km/Vmax), and uptake is linear in exposure
  in every dose series found (amendment 1, `docs/PREREG_sourced_values.md`). `r_eq` is fitted to
  the root; `λ_eq` is numerical (converged).
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
├── reproduce_demo.py                 # entry point: Yamazaki 2020 site-year, dry basis (default | --lipid | --free)
├── build_parameters.py               # (re)assembles params/parameters.json from source tables + refit_sourced.csv
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
├── params/                           # parameters.json (CANONICAL) + source CSVs + refit_sourced.csv + value_provenance.csv
├── data_obs/                         # observed tables: Yamazaki (DRY basis), Li 2025, neutral tables
├── validation/                       # maintained analyses (refit_sourced, transfer, ORYZA checks, neutral, HYDRUS)
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
- **Tier 0** inputs/known: `M_k(t), Q_TP(t), C_w^o(t), f_d, γ_k≈0`; sourced values or sourced
  ranges (run at the midpoint, reported as a band): organ composition, `φ`, `T_C,Ph`, the shoot
  transpiration shares.
- **Tier 1** BAF-identifiable, fitted to Yamazaki (`validation/refit_sourced.py`): `r_eq` (root),
  `L_Ph` (grain), and for the straw `f_xy` per congener (the default, `*_fit`) or the lipid-loading
  constants (`*_lipid`; terms the data do not identify are removed by backward elimination).
- **Tier 2** need inhibitor/kinetic data: channel vs carrier entry, influx vs efflux asymmetry —
  not in the shipped model (the linear entry lumps them).
- **Tier 3** QSPR/measurement (chain-length resolved): `K_prot, K_PL, K_cw`.
- **Identifiability**: BAF data constrain only the lumped root equilibrium `r_eq` — channel vs
  carrier are **not separable** from BAF alone (need inhibitor experiments).

## 6. Current status

Every number below is on the CURRENT model and was produced by the script named next to it.
Superseded analyses are kept in `validation/archive/` (its README says why each one was
retired); git history holds the older states of this file. Every input the model runs on,
with its source or its sensitivity range, is in **`docs/ASSUMPTIONS.md`**.

### 6.1 What the model is
- **No value without a source** (amendments 1 and 2a, `docs/PREREG_sourced_values.md` 5.1-5.10, 7.1-7.13):
  every input is a value read in a source, a fit the data identify, a removed mechanism, or a sourced range
  run at its midpoint with its effect reported (`params/input_range_band.json`). Inventory:
  `params/value_provenance.csv` (guarded by `tests/test_value_provenance.py`; the default-scope rows still
  open are named there); quotes: `docs/literature_db/sourced_inputs.csv`.
- **Default configuration** (re-decided after the amendment-2a re-fit, PREREG §8): `simulate()` runs the
  **redistributed shoot** (`shoot="redistributed"`, `NStemLeafModel`: N = 4 stem segments + leaf; each organ
  keeps the solute its own transpiration deposits — retention 1, derived; the sheath + culm share of
  transpiration 0.29 (0.22-0.36) split evenly over the segments; the panicle's share 0 before heading, 0.21
  (0.09-0.33) at DVS 1 and 0.07 (0-0.14) at DVS 2, linear between) with the **free anion and `f_xy` fitted
  per congener** (`redistributed_fit`: `f_xy_source="fit"`, `lipid_loading=False`) and the **linear
  equilibrium root entry** `λ_eq·(r_eq·C_w^o − C_w,root)`. `model_api.DEFAULT_SHOOT` /
  `DEFAULT_F_XY_SOURCE` / `DEFAULT_LIPID_LOADING`; `resolve_loading` fills an unspecified loading per shoot
  (`shoot="4pool"` alone runs `4pool_lipid`: the single-straw core has no `_fit` variant). On Yamazaki the
  default is a **reproduction** (f_xy fitted to the straw), not a prediction of the straw.
- **Variants, each with its own fitted values** (`model_api.variant_key`; parameters.json `r_eq__<v>`,
  `L_Ph__<v>`, `f_xy__<v>`, `lipid_loading[shoot]`): `redistributed_fit` (default), `redistributed_lipid`
  (amendment 1's default: one `f_xy` and one `g_xy` for every chain), `redistributed_free` (monotone a-priori
  `f_xy`), `4pool_lipid`, `4pool_free`. A switch never runs another variant's values. A record without a fit
  (GenX, SMILES) takes the fitted congeners' median.
- **Root entry**: the GHK + Michaelis–Menten pair left the model (no source for E_m in a paddy root nor
  for any PFAS Km/Vmax; every dose series found is linear). `r_eq` is fitted to Yamazaki's root in closed
  form (the model is linear in it); `λ_eq` = 1e6 is numerical (converged).
- **4-pool core** (`pfas_rice_plant_module_4pool_surf`): basis-A 4-compartment ODE, xylem loading `f_xy`,
  mass-conserving phloem, grain formation gate (a step at the 1e-9 kg organ-mass floor), leaf senescence
  loss flux. `pfas_rice_plant_module` is an alias of it. Phloem constants `phi`, `T_C_Ph` have no class
  default: callers pass the sourced midpoints (`model_api.PHI`, `T_C_PH`).
- **Organ water by stage** (amendment 2a, PREREG 7.1; `organ_water="stage"`): root, stem and leaf θ_fw
  follow DVS linearly between sourced anchors (`params/rice_theta_stage_sourced.csv`: root 0.8245 at DVS ≤ 0.65
  → 0.6365 at DVS 2; stem 0.816 → 0.7095; leaf 0.731 → 0.659); the fresh mass, the binding factor and the
  harvest operator read the same θ(t). Brown rice is not staged (0.21). `organ_water="constant"` is amendment 1.
- **Organ mass is FRESH weight** (`src/organ_mass.py`): crop drivers produce DRY matter, converted with θ_fw.
- **Grain = brown rice** (PREREG 7.2; `grain_basis="brown"`, `src/panicle.py`): brown = (WSO − rachis) ×
  (1 − h(t)), h the husk share of rough grain (1 at heading; Ansari 2017 Nipponbare course scaled to the
  maturity anchor 0.2075 (0.161-0.254)); rachis = 0.5 (0-1) × WSO at DVS 0.95. Brown rice takes **no xylem**
  (phloem only; Zee 1972, Oparka & Gates 1981): the panicle's transpiration deposit stays in a terminal husk +
  rachis pool (per kg dry), reported as `conc_dw["husk"]`, with the rough grain `conc_dw["rough"]`.
  `grain_basis="panicle"` is amendment 1 (ORYZA's whole panicle). A measured growth table's grain is taken as given.
- **Phloem = two streams from the leaf** (`phloem_flows`): to the grain, the sap its DRY growth draws
  (`grain_dry_growth` × `T_C_Ph`, 6.1 L/kg (4.4-7.8)); to the root, φ·Q_TP (φ 0.06 (0.02-0.10)).
- **Harvest observation operator** (`src/harvest.py`): the ODE carries GREEN leaves; dead blades are
  reconstructed from the ORYZA leaf death rate and kept in the leaf and straw (`dead_leaf_attached=1`:
  rice does not shed its dead blades; 0 = all lost, a sensitivity). Results carry `conc_dw` /
  `baf_final_dw` (dry basis) next to the fresh-weight values.
- **Default scenario `"yamazaki2020"`** (`src/site_scenario.py`): JMA Tsukuba 2020 weather, the ORYZA2000
  port (IR72.D92 verbatim) transplanted 2020-05-11, sown 04-13 (nursery 25.5 d after emergence, the 25- and
  26-d drivers interpolated), seedbed cover 2.5 C (0-5.0), 18.5 hills/m², 4.5 plants/hill (4-5), CO₂ 414 ppm;
  phenology calibrated DRATES-style to each cultivar's measured heading (Kasalath 08-03, Koshihikari 08-11) and
  harvest (= DVS 2); transpiration from ORYZA eq. 4.4 on FAO-56 terms. **Japonica crop sets** (Gao 2021; Wu
  2022 Huanren, Donggang, Liaozhong, Dawa; `src/crop_sets.py`) are a categorical band end, IR72 runs.
  `simulate()` runs the mean-date harvest; the Yamazaki observable (`yamazaki_harvest_baf`) combines the two
  cultivar runs by geometric mean. The generic scenario's planting density is defined once
  (`oryza2000.GENERIC_HILLS_PER_M2`, the crop file's 25); its BAFs do not depend on it.
- **Observations** (`data_obs/obs_baf_Yamazaki.csv`, built from the SI by
  `docs/literature_db/tools/yamazaki_obs_from_si.py`): dry basis; water from Table S16; straw = the height
  bands weighted by their share of straw dry mass (Tanaka et al. 2022 Fig. 3, mean of the 10 and 49 d after
  heading profiles, each profile a range end); brown rice at a white share of 0.88 (0.86-0.90); < MQL cells
  at the midpoint of [0, MQL); PFBS grain not scored (its censoring bound sets PFBS's L_Ph range).

### 6.2 Calibration on the site-year — `validation/refit_sourced.py` (pre-registered; amendment 2a)
Every variant fits `r_eq` (root, closed form) and `L_Ph` (grain) per congener; `*_fit` adds `f_xy` per
congener (straw); the lipid variants add two constants shared by all chains (backward elimination left
`f_xy = a` 0.0659 and `g_xy` 0.00712 on `redistributed_lipid`). Yamazaki two-harvest observable, log10 RMSE
over 32 cells (records: `validation/refit_sourced.json`; amendment 1's in `validation/archive/amendment1/`;
outcome: `docs/PREREG_sourced_values.md` §8):

| variant | all | root | straw | grain |
|---|---|---|---|---|
| **redistributed_fit** (default; a reproduction) | **0.116** | 0.000 | 0.156 | 0.127 |
| redistributed_lipid | 0.154 | 0.000 | 0.262 | 0.011 |
| 4pool_lipid | 0.151 | 0.000 | 0.258 | 0.000 |
| redistributed_free | 1.261 | 0.000 | 1.639 | 1.461 |
| 4pool_free | 1.309 | 0.000 | 1.794 | 1.394 |

- **Default**: every cell is reproduced except PFDoDA (f_xy at 1; straw −0.52, grain −0.40 log). Every
  congener with a grain cell fits `L_Ph` > 0 in every variant (the panicle deposit now goes to the husk).
- **Free anion**: from PFOA up, and PFOS, the grain is out of reach at any `L_Ph` and the straw is under
  by 0.8–3.1 log; the stored `L_Ph` of an unreachable chain is merely large.
- **Input-range band** (`validation/input_range_band.py` → `params/input_range_band.json`; PREREG 5.4, 7):
  every sourced range moved to its ends with the fits held, now with the stage-water anchors, the husk
  share, the rachis, the husk phloem draw and the crop sets; median over 11 congeners grain ÷2.97 / ×3.65,
  straw ÷1.40 / ×1.89, root ÷1.29 / ×1.29. The crop sets lead the straw band (11/11) and the grain band
  (7/11; φ 3/11), root composition the root band. The one-at-a-time envelope covers the all-corner run in
  38/66 cases only (inputs interact). The incumbent's band (E19/E23 record):
  `validation/input_range_band__redistributed_lipid.json`.
- **Numerics** (`validation/numerics_convergence.py`): exponential integrator (`src/linear_expm.py`),
  `λ_eq` 1e6; the terminal pools (grain, husk) are integrated as burden y = M·C (their dilution cancels
  exactly), and brown rice's phloem is not gated by the mass floor (the gate stays on the panicle basis).
  This fixed a grain not converged in the floor (PREREG §8; the re-fit was repeated on it).
  The brown rice's pre-harvest water moves no dry-basis number.
- Pre-registered expectations (amendment 2a): E17–E19, E21 (direction; Choi's husk/brown ratio is not
  reached: PFOA −0.19, PFOS −0.89 log on the incumbent), E22, E24 held; E20 (the straw moved 0.086 with the
  shared lipid constants) and E23 (the crop sets, not the panicle share, lead the grain band) wrong.
- `python reproduce_demo.py` prints the default (`--lipid`, `--free`: the other redistributed variants).

### 6.3 Transfer tests (no parameter fitted to the transfer data) — `validation/transfer_sourced.py`
log10 RMSE (records: `validation/transfer_sourced.json`, `params/transfer_scores.json`). Li grain on brown
rice (the decision basis, PREREG 7.13); the rough-grain value in brackets.

| dataset | 4pool_free | redistributed_free | 4pool_lipid | redistributed_lipid | **default (redistributed_fit)** |
|---|---|---|---|---|---|
| Tang 2026 per-organ TF, 0.1 µg/g | 1.876 | 1.738 | 0.916 | 0.687 | **0.686** |
| Kim 2019 grain BAF excl PFOA, point / [0, MDL] interval | 2.037 / 1.777 | 2.118 / 1.866 | 0.559 / 0.445 | 0.554 / 0.444 | **0.562 / 0.443** |
| Li 2025 TF by stage, shoot / grain [rough] | 1.494 / 1.815 | 1.279 / 1.863 | 0.937 / 1.373 | 0.364 / 1.361 [0.864] | **0.467 / 1.334 [1.042]** |

- **Decision (PREREG 5.3, §8; applied by the user)**: `redistributed_fit` is better than the incumbent
  (`redistributed_lipid`) on the Li grain (−0.027, 95 % [−0.045, −0.003]) and ties on every other endpoint
  (Li shoot point 0.364 → 0.467, inside the tie). The margin is small and the chosen variant reproduces
  Yamazaki's straw.
- **Grain against amendment 1** (0.489 Kim, 1.019 Li): worse, since brown rice takes no xylem; Li's grain
  is under-predicted; on the rough basis Li grain is 0.864 (incumbent).
- **Selection caveat**: Tang, Kim and Li chose the default, so its scores on them — and the band built
  from them — may be optimistic.
- **Honesty band**: `model_api.APRIORI_OOS_SCORES` (read from `params/transfer_scores.json`; no number
  typed in code); their median 0.562 log10 (**×3.65**) is the band the Simple view shows. Yamazaki is not
  in the list (in-sample).
- **Kim PFOA per site** (`kim2019_pfoa_anchor.py`): censored per-site geomean 9.94 L/kg [4.57, 20.16];
  the default's 0.538 L/kg is below the interval (−1.27 log).
- **Inverse** (`estimate_exposure_bayesian`, Kim PFOA brown rice 0.349 µg/kg dw, void water
  0.0787 µg/L): the default estimates 0.649 µg/L (8.2× over; 95 % 0.051–8.20, which contains the
  measured value) (`validation/policy_brief_runs.py`).
- Earlier decisions (the amendment-1 records, the pre-amendment default switch, the panicle-size share, the
  gate numerics) are in `validation/archive/` and git history; their numbers do not describe the current model.
- Every transfer dataset runs on the calibration site-year drivers (Tang/Kim/Li sites and dates
  are not in the repo) — see `docs/ASSUMPTIONS.md` §5.

### 6.4 Mechanism and identifiability checks
- **Root entry**: linear (amendment 1). The carrier-vs-bypass and dose-series analyses (archived)
  found an entry addition necessary but no carrier constant sourced, and every dose series linear.
- **Bayesian identifiability** (`bayesian_inverse_demo.py`; dry tissue CONCENTRATIONS, since a BAF hides
  the exposure in a linear model): (A) Q_TP scale and Cwo recover with transport fixed (corr −0.78,
  cond 27); (B) Q_TP and f_xy trade off (corr −0.95, cond 87: only their product is pinned); (C) Cwo and
  `r_eq` are an exact ridge (cond 1.3e16: the linear entry carries only their product).
  `estimate_exposure_bayesian` (the app's inverse) estimates Cwo from DRY-weight tissue data with
  transport fixed.
- **HYDRUS-1D coupling** (`hydrus_coupled_run.py`; Method A, the site-year transpiration drives the
  soil run): weakly sorbed short chains leach during flooding, long chains stay buffered. The paddy
  soil inputs (f_oc, bulk density, percolation, drainage calendar, the irrigation-water source) are
  not yet site-sourced (amendment 2), so no soil-scenario number is reported here.

### 6.5 Crop driver — the ORYZA2000 port
- `src/oryza2000.py` ports ORYZA2000 potential production from the manual (IR72.D92 verbatim;
  carbon balance closes to ~1e-15). On IR72DS92 (the manual's calibration example — a fidelity
  test) it fits all four organ series and the phenology; the reduced re-coding `oryza_growth`
  matures 12 d late and under-fills panicles (`oryza2000_ir72ds92.py`).
- Independent checks: Wuxi FACE organ partitioning (`oryza2000_wuxi_partitioning.py`) — port
  clearly closer on the leaf/stem/spike split in every design. Shizukuishi FACE
  (`oryza2000_kim_face*.py`; living-biomass ratios + paired bootstrap): with cultivar phenology the
  port is better if the observed senescent share is read as a floor (detached blades missing) and
  indistinguishable if it is read as complete (no blade detached); Kim 2003 does not report which
  (`docs/rice_sample_processing_KR.md` §4), so both readings are printed. With IR72 phenology the
  shipped driver is better on the cool grid cell under the complete reading. The reduced driver
  never matures at temperate sites; the port does once the weather is not the cool grid cell.
- The site-year scenario therefore runs the port; `oryza_growth` / `growth_rice` remain as
  `scenario="generic"` alternatives.

### 6.6 Neutral organics (Briggs/Kow path; nothing fitted)
`src/neutral_dpu.py` — the same ODE with z = 0, no carrier, K_PW = W + L·a·Kow^b, f_xy = TSCF(log Kow),
phloem off. Composition (addendum 5.10, `src/neutral_composition.py`, `params/neutral_composition_sourced.csv`):
lipid carried on a dry basis and restated at the organ water (rice root = Li 2019 cereal roots, stem/leaf =
rice straw ether extract, grain = brown rice crude fat); each validation table runs on its own plant's water
and lipid, with its dry-weight observations converted by that water. The tables run on the generic drivers.
Current a-priori scores (log10 RMSE at the midpoints, [corners]; `validation/neutral_scores.py` →
`params/neutral_scores.json`): Liu 2023 rice root 0.212 [0.208, 0.278] (equilibrium 0.204; read as dry
weight 0.821) · Ge 2017 leaf TF 0.631 (stem rows unscored: no vegetative rice stem water found) · Kodešová
2019 0.157 (equilibrium 0.172) · Li 2019 hydroponic 0.677 [0.525, 1.115] · Li 2019 soil 0.418 (equilibrium
0.461) · Briggs 1983 barley stem 0.357. Expectations: N1 wrong (both Li halves moved > 0.05, opposite ways),
N2 and N3 held.
Weak electrolytes (`pKa=`, `weak_electrolyte_tscf.py`, 67 ionisable Schriever rows): the direction is
supported but the scale is not; nothing adopted (`docs/neutral_dpu_validation.md`). Plant–air exchange
is opt-in and exactly zero at K_AW = 0; volatile compounds are not tested against data.

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
- **Reproduction**: `python reproduce_demo.py` (the default on the site-year, dry basis; `--lipid`,
  `--free`). Re-fit: `python validation/refit_sourced.py` (~30 min; every variant, checkpointed,
  `--resume`; writes `validation/refit_sourced.json`, `params/refit_sourced.csv`,
  `params/lipid_loading_sourced.json`), then `python build_parameters.py`.
- **App**: `pip install -r requirements-app.txt && streamlit run app.py` (Simple = Korean, Expert =
  English; crop drivers site-year | generic; compound class congener | SMILES | neutral log Kow).
- **API**: `model_api.simulate("PFOA")` (default configuration, site-year, mean harvest) ·
  `simulate("PFOA", **FOUR_POOL_FREE)` (previous default) · `yamazaki_harvest_baf(c)` (the
  two-run observable) · `simulate(..., scenario="generic", season=…, biomass=…)` · `simulate_nstem_leaf`
  · `simulate_neutral(log_kow, …)` · `tang_tf_validation(c)` · `estimate_exposure_bayesian(...)` ·
  `apportionment(c)`.
- **Amendment 1 records**: `validation/transfer_sourced.py` (transfer scores + decisions),
  `validation/input_range_band.py` (the sourced-range band), `validation/neutral_scores.py`,
  `validation/prereg_a1_outcome.py` (reads E7–E16 off the records),
  `validation/numerics_convergence.py` (the numerical settings).
- **Transfer tests** (per dataset): `oos_tang_lipid.py`, `kim2019_censored.py`,
  `kim2019_pfoa_anchor.py`, `li2025_stage_resolved.py`.
- **Mechanism checks**: `bayesian_inverse_demo.py` (`--emcee`). Superseded analyses (default switch,
  sensitivity, panicle share, gate numerics, carrier vs bypass, dose series, θ_fw) are in
  `validation/archive/`.
- **Crop driver checks**: `oryza2000_ir72ds92.py`, `oryza2000_kim_face.py`,
  `oryza2000_kim_face_phenology.py`, `oryza2000_kim_face_nursery_events.py`,
  `oryza2000_kim_face_gao_japonica.py`, `oryza2000_wuxi_partitioning.py`.
- **Neutral path**: `validation/neutral_dpu_validation.py --obs data_obs/neutral_obs_liu2023.csv`
  (or `…ge2017.csv`, `--lipid-source both`, `--mode equilibrium`), `neutral_scores.py`,
  `li2019_soil_table.py`, `kodesova2019_carbamazepine.py`, `briggs1983_stem.py`,
  `briggs1983_shoot.py`, `hwang2017_lettuce.py`, `schriever2020_tscf.py`,
  `weak_electrolyte_tscf.py` (`--fast`).
- **HYDRUS-1D** (vendored source): `cp external/hydrus_source/makefile external/hydrus_source/source/ &&
  (cd external/hydrus_source/source && make)`; `python src/soil_hydrus.py`;
  `python validation/hydrus_coupled_run.py`; `python validation/cwo_kleach_calibration.py`. On Claude
  Code on the web the SessionStart hook builds the engine.
- **Policy briefs** (Korean slide scripts, every number from a script): `docs/POLICY_BRIEF_KR.md` (PFAS;
  `python validation/policy_brief_runs.py` → `docs/policy_brief_results.csv`, ~3 min) and
  `docs/POLICY_BRIEF_BIOCIDES_KR.md` + `docs/bat_census_biocides.md` (EU biocide census on the neutral path;
  `python validation/bat_census_biocides.py` then `policy_biocides_figures.py`, ~3 min). Run the census before
  the PFAS brief (its slide 8 reads the census record). Guards: `test_policy_brief.py`,
  `test_policy_biocides_figures.py`, `test_bat_census_biocides.py`. Rebuilt on the default configuration; the
  first version (ccy5123/pfas-rice-model#65) carried pre-correction numbers.
- **Compound lookup** (EPA CompTox): `python src/chem_lookup.py carbamazepine` (needs `CTX_API_KEY`;
  base URL `https://comptox.epa.gov/ctx-api`).
- **Tests**: `pytest` (486 pass with the full stack, plus 1 module-level skip for `test_sci_adk_rigor.py` without sci-adk; ~8 min). CI (`.github/workflows/tests.yml`) runs the whole suite;
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
1. **Measurements that would narrow the sourced ranges** (`docs/ASSUMPTIONS.md`,
   `params/value_provenance.csv`): japonica crop parameters for the site's cultivars (the largest straw and
   grain lever since amendment 2a), a rice panicle transpiration share by stage, rice phloem sap per grain DM
   (T_C,Ph), φ, organ composition over a season, the rachis share of the panicle, rice-root cell-wall / plaque
   PFAS sorption (K_cw, K_surf), Yamazaki 2023 SI raw values (MQL flags, planting density, nursery), the
   site's paddy soil.
2. **Transfer datasets on their own drivers**: Tang (site, dates), Kim (field calendars), Li (stage
   dates) — each currently runs on the Yamazaki site-year.
3. **Grain**: brown rice is phloem-fed only, and the default under-predicts Li's grain by stage (1.334;
   rough basis better) and Kim's PFOA; the husk/brown ratio is below Choi 2018's (PFOS by ~1 log). The
   timing and route of grain loading are open; harvest-only data cannot decide them.
4. **Amendment 2b** (registered in PREREG 7.9; draft, prototypes and a Korean handoff in
   `docs/amendment2/`): the neutral path (organ water, dry-basis lipid, the composition basis per dataset),
   site-sourced paddy soil and the irrigation-source exposure shape, app values (KNHANES intake, the
   white/brown basis), SMILES K_cw, a volatile-compound data search; then a docs sweep.
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
