# PFAS Rice Compartmental Uptake Model

A mechanistic four-compartment (root / stem / leaf / grain) dynamic model for the uptake
and tissue distribution of **PFAS in rice**, treated as an **ionizable-organic-compound
extension** of the Trapp/Brunetti Dynamic Plant Uptake (DPU) framework and designed to
couple with **HYDRUS-1D** for the soil side.

> Agent context lives in [`CLAUDE.md`](CLAUDE.md). Full model derivation is in [`docs/`](docs/).

## What it does
Solves the four-compartment ODE system for a fully-dissociated PFAS anion, with:
- a **hybrid root uptake** flux (ionic electrodiffusion + saturable carrier),
- compartment **binding factors** (protein / phospholipid / cell-wall),
- **xylem** (upward) and **phloem** (grain-loading) transport,
- the grain as a **terminal accumulator**.

Soil pore-water free concentration `C_w^o(t)`, transpiration `Q_TP(t)`, and tissue growth
`M_k(t)` are **external inputs** (Method A: one-way coupling from HYDRUS-1D/Phydrus).

## Quickstart
```bash
pip install -r requirements.txt
python src/pfas_rice_plant_module.py
```
This runs a synthetic demo (placeholder parameters — **not calibrated**) and prints the
electrochemical number, binding factors, and final tissue concentrations/BAFs.

## Layout
- `src/pfas_rice_plant_module.py` — the plant ODE module (Method A).
- `docs/` — LaTeX model report (`pfas_rice_compartmental_model`) and the corrected
  neutral DPU base summary (`dpu_model_summary_corrected`).
- `external/hydrus_source/` — git submodule of `phydrus/source_code` (HYDRUS-1D 4.08
  FORTRAN, LGPL-3.0) for the future tight (Method B) coupling.

## Status
Structure + derivation complete; Python module runs and reproduces the expected structural
behaviour. **Not yet calibrated** — see `CLAUDE.md` §6 and the prioritized next tasks in §9.

## References
Brunetti et al. 2019 (WRR, DOI 10.1029/2019WR025432), 2021 (ES&T, 10.1021/acs.est.0c07420),
2022 (J. Hazard. Mater., 10.1016/j.jhazmat.2021.127008).
