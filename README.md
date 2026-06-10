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
- **xylem** (upward) and **phloem** (grain-loading) transport, mass-conserving,
- **TSCF-limited root→xylem loading** (`f_xy`): the anion is retained in the root and
  translocates poorly, reproducing the empirical **root > straw > grain** ordering,
- the grain as a **terminal accumulator**.

Soil pore-water free concentration `C_w^o(t)`, transpiration `Q_TP(t)`, and tissue growth
`M_k(t)` are **external inputs** (Method A: one-way coupling from HYDRUS-1D/Phydrus). A
**Freundlich paddy soil** sub-model (`src/soil_paddy.py`) can supply `C_w^o(t)` from a soil
inventory (redox-dependent, non-linear sorption), and a **Tier-1 calibration** module
(`src/calibration.py`) fits the identifiable parameters to observed tissue BAFs.

## Quickstart
```bash
pip install -r requirements.txt
python src/pfas_rice_plant_module.py    # plant demo: prints the root > straw > grain check
python src/soil_paddy.py                # Freundlich soil → C_w^o(t) under a flooding schedule
python src/calibration.py               # synthetic parameter recovery + identifiability demo
pip install pytest && pytest            # structural, soil, mass-conservation & calibration tests
```
The plant demo uses placeholder parameters (**not calibrated**) and prints the electrochemical
number, binding factors, final tissue concentrations/BAFs, and the `root > straw > grain` check.

## Layout
- `src/pfas_rice_plant_module.py` — the plant ODE module (Method A).
- `src/soil_paddy.py` — Freundlich paddy soil sorption (`C_w^o(t)`) + CSV / HYDRUS input adapters.
- `src/calibration.py` — Tier-1 calibration to observed BAFs (scipy); synthetic recovery.
- `tests/` — pytest suite (plant, soil, calibration) locking in the structural results,
  exact mass conservation, and parameter recovery.
- `docs/` — LaTeX model report (`pfas_rice_compartmental_model`) and the corrected
  neutral DPU base summary (`dpu_model_summary_corrected`).
- `external/hydrus_source/` — git submodule of `phydrus/source_code` (HYDRUS-1D 4.08
  FORTRAN, LGPL-3.0) for the future tight (Method B) coupling.

## Status
Structure + derivation complete; Python module runs and reproduces the expected structural
behaviour, including the empirical **root > straw > grain** ordering via the TSCF loading
factor. **Not yet calibrated** — see `CLAUDE.md` §6 and the prioritized next tasks in §9.

## References
Brunetti et al. 2019 (WRR, DOI 10.1029/2019WR025432), 2021 (ES&T, 10.1021/acs.est.0c07420),
2022 (J. Hazard. Mater., 10.1016/j.jhazmat.2021.127008).
