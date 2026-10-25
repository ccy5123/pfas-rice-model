# PFAS–Rice 4-Compartment Uptake Model

**English** | [한국어](README_KR.md) &nbsp;·&nbsp; App manual: [English](docs/MANUAL_EN.md) | [한국어](docs/MANUAL_KR.md)

Mechanistic uptake model for **permanently-anionic PFAS in paddy rice** (*Oryza sativa*): an
ionizable-organic-compound extension of the DPU / Trapp framework over four compartments
(root → stem → leaf → grain). Covers 12 congeners with Yamazaki observations (PFCA C4–C12,
PFSA C4/C6/C8) plus GenX (provisional), any PFAS given as a SMILES, and neutral organics given
by log Kow.

> **Canonical sources**: `CLAUDE.md` (current status and results), `docs/ASSUMPTIONS.md` (every
> input with its source or its sensitivity range), `params/parameters.json` (parameters),
> `docs/MANUAL_EN.md` / `docs/MANUAL_KR.md` (the app). `docs/HANDOFF_*.md`, `docs/VALIDATION_*.md`,
> `docs/*_exploration.md` and `validation/archive/` are history, not the current state.

---

## Quickstart

```bash
pip install -r requirements.txt
python reproduce_demo.py          # Yamazaki 2020 site-year, dry basis, fitted transport (a reproduction)
python reproduce_demo.py --rec    # monotone physical f_xy: straw is a prediction
pytest                            # tests
```

```python
import model_api as api                       # (src/ on the path)
r = api.simulate("PFOA")                      # default: redistributed shoot + lipid loading, site-year
r["baf_final_dw"]                             # root / stem / leaf / grain / straw, dry basis
api.yamazaki_harvest_baf("PFOA")              # the observable: two harvest runs, geometric mean
api.simulate("PFOA", scenario="generic", season=120.0, biomass="oryza2000")   # other drivers
api.simulate("PFOA", **api.FOUR_POOL_FREE)    # the previous default: 4-pool, free anion
api.simulate_neutral(2.45, name="carbamazepine")      # a neutral organic by log Kow
```

## What the model runs on

- **Drivers** — the default scenario is the calibration site-year (`src/site_scenario.py`): JMA
  Tsukuba 2020 weather, the ORYZA2000 potential-production port (`src/oryza2000.py`, IR72.D92)
  transplanted on Yamazaki's date with phenology calibrated to its heading and harvest dates, and
  ORYZA2000 transpiration on FAO-56 terms. `scenario="generic"` uses a climatology
  (`forcing_rice`) and a chosen biomass driver over a chosen season.
- **Organ mass is fresh weight**; crop drivers produce dry matter and are converted
  (`src/organ_mass.py`).
- **Harvest observation operator** (`src/harvest.py`): dead leaf blades are reconstructed from the
  crop model's leaf death and kept in the leaf and straw; results are compared on dry basis.
- **Default configuration** (pre-registered switch, `validation/default_switch.py`): the
  **redistributed shoot** (N = 4 stem segments + leaf; each shoot organ retains the solute its own
  transpiration deposits) with **K_PL-gated lipid loading** and the carrier root entry.
  `shoot="4pool"` is the single-straw core; `api.FOUR_POOL_FREE` the previous default.
- **Transport parameters** are fitted on the site-year with the 4-pool shoot
  (`validation/refit_site_year.py`, pre-registered): `f_xy_source="recommended"` (monotone physical
  f_xy, with L_Ph and κ_d fitted to root and grain) or `"fit"` (all three fitted per congener; a
  saturated reproduction; free anion only). Lipid loading supplies its own f_xy (six globals fitted
  on Yamazaki).

## Results in brief (details and scripts: `CLAUDE.md` §6)

- **Reproduction** (4-pool Fit A, saturated): log10 RMSE 0.162 on Yamazaki's root/straw/grain.
  With the monotone f_xy the straw is a prediction and misses by log10 1.40.
- **Transfer to independent data** (no parameter fitted to them): the default configuration is
  off by log10 0.64 on Tang 2026 per-organ TF, 0.47 on Kim 2019 grain and 0.36 / 0.51 on Li 2025
  TF by stage; the previous default (4-pool free anion) by 1.71, 2.33 and 1.26 / 1.68. On Kim the
  ranking of the two lipid variants depends on the non-detects.
- These three datasets were not fitted, but they chose the default (the switch scored the
  candidates on them), so the default's scores on them — and the ×3.1 band — may be optimistic.
- **On Yamazaki** (in-sample for the lipid constants) the default scores 0.380, with the grain
  over-predicted 2–7× for most congeners.
- **Unmeasured inputs that matter** for the default: the carrier constants (Km above all), E_m,
  the shoot levers `retention` and `lam_grain` (the grain's share of the stream), and whether dead
  blades stay in the harvested straw (`docs/ASSUMPTIONS.md`).
- **Neutral organics** (nothing fitted): rice root partition log10 0.21 (Liu 2023), per-organ TF
  0.31 (Ge 2017), barley stem 0.30 (Briggs 1983).

The general-audience view shows every absolute number with a ×/÷ band set by the default
configuration's transfer error (`model_api.APRIORI_OOS_SCORES`: median log10 0.49, ×3.1).

## Soil side — HYDRUS-1D (Method A)

The vendored HYDRUS-1D 4.08 source is built with gfortran and driven via `phydrus`; the site-year
transpiration drives the soil run, which returns the pore-water `C_w^o(t)` and root water uptake.

```bash
cp external/hydrus_source/makefile external/hydrus_source/source/
(cd external/hydrus_source/source && make)
pip install phydrus
python validation/hydrus_coupled_run.py
```

Weakly sorbed short chains leach during flooding (default configuration: PFBA grain BAF 6.90 at constant exposure →
0.69), long chains stay buffered. `simulate(cwo_profile="flooded")` gives an engine-free analytic
shape whose `k_leach` is calibrated per congener to the HYDRUS run. On Claude Code on the web the
SessionStart hook builds the engine.

## Interactive app

```bash
pip install -r requirements.txt
streamlit run app.py
```

A plant-and-soil map coloured by concentration, tissue dynamics, drivers, the comparison with
Yamazaki, the Tang check, a Bayesian inverse (exposure from measured tissue), six exposure modes
(parametric, custom tables, HYDRUS CSV, live HYDRUS-1D, soil inventory, biomonitoring), the
mechanism switches (shoot model, carrier/bypass, lipid loading) and three compound classes (congener, SMILES,
neutral log Kow with an optional EPA CompTox lookup). Simple mode is in Korean, Expert in English.
See `docs/visualization_tool.md` and the manuals.

## Layout

```
├── reproduce_demo.py / build_parameters.py / app.py
├── src/            model code (core ODE, harvest operator, site scenario, ORYZA2000 port, API, plots, soil)
├── params/         parameters.json + source tables + refit_site_year.csv + cwo_kleach.csv
├── data_obs/       observations (Yamazaki dry basis, Li 2025, neutral tables)
├── validation/     maintained analyses; archive/ = superseded ones
├── docs/           ASSUMPTIONS.md, model report (.tex), literature DB, manuals
├── external/hydrus_source/   vendored HYDRUS-1D source
└── tests/
```

## Model in one screen

```
B_k = θ_fw + (1 − θ_fw) · ( f_prot·K_prot + f_PL·K_PL + f_cw·K_cw )      [L/kg fw]
```

Root influx = GHK electrodiffusion (anion exclusion e^N ≈ 107 at E_m = −120 mV) + a saturable
Michaelis–Menten carrier. Xylem loading `f_xy · C_root/B_root` + the K_PL-gated lipid term
`g_xy · C_root` (on by default). The leaf feeds two phloem streams at concentration
`L_Ph · C_leaf/B_leaf`: to the grain, the sap its dry-matter growth draws; to the root, a fraction
φ of transpiration. In the default (redistributed) shoot the stream passes N stem segments, the
leaf and the grain, and each organ retains part of the solute its own transpiration deposits.
Grain and leaf are terminal accumulators, so the model is dynamic.

## Open

Measurements that would replace assumptions are listed in `docs/ASSUMPTIONS.md` §7 (Yamazaki SI
raw values, carrier kinetics, a panicle transpiration share, rice phloem sap volume, organ water
content over a season, rice-root cell-wall / plaque sorption) — and, for the default shoot, the
`retention` and `lam_grain` levers. The transfer datasets still run on the calibration site-year's
drivers because their own sites and dates are not in the repo. How the growth and PFAS papers
physically processed their rice (sampling, organ separation, washing, drying, milling), and a
checklist for new sampling, is in `docs/rice_sample_processing_KR.md`.
