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
python reproduce_demo.py          # Yamazaki 2020 site-year, dry basis, the default (a reproduction)
python reproduce_demo.py --lipid  # shared lipid constants: the straw is a prediction
python reproduce_demo.py --free   # monotone physical f_xy
pytest                            # tests
```

```python
import model_api as api                       # (src/ on the path)
r = api.simulate("PFOA")                      # default: redistributed shoot, f_xy fitted per congener, site-year
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
  (`src/organ_mass.py`) with each organ's water at the run's development stage (sourced anchors).
- **Grain = brown rice**, fed by the phloem only; the panicle's transpired solute stays in a husk +
  rachis pool (`src/panicle.py`; reported as `husk` and `rough`).
- **Harvest observation operator** (`src/harvest.py`): dead leaf blades are reconstructed from the
  crop model's leaf death and kept in the leaf and straw; results are compared on dry basis.
- **No value without a source** (`docs/PREREG_sourced_values.md`, `params/value_provenance.csv`):
  every input is a value read in a source, a fit the data identify, a mechanism removed for lack
  of either, or a sourced range run at its midpoint whose effect is reported
  (`params/input_range_band.json`).
- **Default configuration**: the **redistributed shoot** (N = 4 stem segments + leaf; each shoot
  organ keeps the solute its own transpiration deposits; sourced stem and panicle shares of
  transpiration) with the **free anion and `f_xy` fitted per congener** (`redistributed_fit`, chosen by the
  pre-registered transfer rule after amendment 2a) and the **linear equilibrium root entry**
  `λ_eq·(r_eq·C_w^o − C_w,root)`. `lipid_loading=True` is amendment 1's default; `shoot="4pool"` is the
  single-straw core.
- **Fitted parameters** (`validation/refit_sourced.py`, pre-registered), per variant: `r_eq` (root)
  and `L_Ph` (grain) per congener; `f_xy` per congener (straw) in the `*_fit` variants, the default — a
  saturated reproduction of Yamazaki; the lipid-loading constants across congeners in `*_lipid` (terms the
  data do not identify removed by backward elimination).

## Results in brief (details and scripts: `CLAUDE.md` §6)

- **On Yamazaki** (the calibration data; `r_eq`, `f_xy` and `L_Ph` fitted per congener) the default
  scores log10 0.116 — a reproduction; only PFDoDA is out of reach. The lipid loading (two shared
  constants, the straw predicted) scores 0.154; the free anion with the a-priori `f_xy` 1.26.
- **Transfer to independent data** (no parameter fitted to them): the default is off by log10 0.69
  on Tang 2026 per-organ TF, 0.56 on Kim 2019 grain (0.44 on the non-detect interval), and 0.47 /
  1.33 on Li 2025 shoot / grain TF by stage (brown rice, phloem-fed only, is under-predicted). It
  replaced the lipid loading by the pre-registered rule (better on the Li grain, tied elsewhere; a
  small margin). The 4-pool free anion: 1.88, 2.04, 1.49 / 1.82.
- These datasets were not fitted, but they chose the default, so its scores on them — and the ×3.65
  band — may be optimistic.
- **Sourced input ranges** move the default's harvest BAF by grain ÷3.0 / ×3.6, straw ÷1.40 / ×1.89,
  root ÷1.29 / ×1.29 (median of 11 congeners); the japonica crop parameter sets lead the straw and grain
  band (`params/input_range_band.json`).
- **Neutral organics** (nothing fitted; each table on its own plant's sourced composition): rice
  root partition log10 0.21 (Liu 2023), carbamazepine roots 0.16 (Kodešová 2019), Li 2019 soil 0.42
  / hydroponic 0.68, barley stem 0.36 (Briggs 1983).

The general-audience view shows every absolute number with a ×/÷ band set by the default
configuration's transfer error (`model_api.APRIORI_OOS_SCORES`: median log10 0.56, ×3.65).

## Soil side — HYDRUS-1D (Method A)

The vendored HYDRUS-1D 4.08 source is built with gfortran and driven via `phydrus`; the site-year
transpiration drives the soil run, which returns the pore-water `C_w^o(t)` and root water uptake.

```bash
cp external/hydrus_source/makefile external/hydrus_source/source/
(cd external/hydrus_source/source && make)
pip install phydrus
python validation/hydrus_coupled_run.py
```

Weakly sorbed short chains leach during flooding, long chains stay buffered (the paddy-soil inputs
are not yet site-sourced, so no soil-scenario number is reported). `simulate(cwo_profile="flooded")` gives an engine-free analytic
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
mechanism switches (shoot model, lipid loading) and three compound classes (congener, SMILES,
neutral log Kow with an optional EPA CompTox lookup). Simple mode is in Korean, Expert in English.
See `docs/visualization_tool.md` and the manuals.

## Layout

```
├── reproduce_demo.py / build_parameters.py / app.py
├── src/            model code (core ODE, harvest operator, site scenario, ORYZA2000 port, API, plots, soil)
├── params/         parameters.json + source tables + refit_sourced.csv + value_provenance.csv + records
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

Root entry = `λ_eq · (r_eq · C_w^o − C_w,root)` (linear equilibrium; the theory's GHK + carrier
pair has no sourced constants). Xylem loading `f_xy · C_root/B_root` + a lipid term `g_xy · C_root`
(on by default). The leaf feeds two phloem streams at concentration `L_Ph · C_leaf/B_leaf`: to the
grain, the sap its dry-matter growth draws; to the root, a fraction φ of transpiration. In the
default (redistributed) shoot each organ draws its own share of the transpiration stream and keeps
the solute it deposits. Grain and leaf are terminal accumulators, so the model is dynamic.

## Open

What is still open is listed in `docs/ASSUMPTIONS.md` (marked **open**) and
`params/value_provenance.csv`: a crop parameter set for Koshihikari/Kasalath, the seasonal course
of the irrigation-water exposure, the soil-side values, the neutral path's composition, the app's
per-volume density. The transfer datasets still run on the calibration site-year's drivers because
their own sites and dates are not in the papers. How the growth and PFAS papers
physically processed their rice (sampling, organ separation, washing, drying, milling), and a
checklist for new sampling, is in `docs/rice_sample_processing_KR.md`. Korean policy-briefing scripts built
from script output: `docs/POLICY_BRIEF_KR.md` (PFAS) and `docs/POLICY_BRIEF_BIOCIDES_KR.md` (EU biocides).
