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
├── CLAUDE.md                         # this file
├── README.md
├── requirements.txt
├── src/
│   └── pfas_rice_plant_module.py     # Method A plant ODE module (runnable)
├── docs/
│   ├── pfas_rice_compartmental_model.tex / .pdf
│   └── dpu_model_summary_corrected.tex / .pdf
├── external/
│   └── hydrus_source/                # git submodule → github.com/phydrus/source_code
├── data/                             # (gitignored) HYDRUS output, BAF datasets, params
└── tests/
```

## 4. Coupling strategy
- **Method A — loose, one-way (CURRENT).** HYDRUS-1D/Phydrus → `C_w^o(t)`, `Q_TP(t)`;
  the plant ODE is solved in Python (`src/pfas_rice_plant_module.py`). No FORTRAN edits.
  Interface = the three arrays in `PlantInputs` (`Cwo`, `Qtp`, `M`).
- **Method B — tight (FUTURE).** Modify `external/hydrus_source` (HYDRUS-1D FORTRAN):
  replace/augment the **root solute-uptake routine** with `j_R`, add the plant module,
  rebuild via `makefile`. `external/hydrus_source/source_mcmc/` provides Bayesian
  calibration machinery.

## 5. Parameter tiers (calibration design)
- **Tier 0** inputs/known: `M_k(t), Q_TP(t), C_w^o(t), N(E,z), f_d, γ_k≈0, T_C,Ph`.
- **Tier 1** BAF-identifiable (lumped): `B_k`, `g_in/g_out`, `Π = Q_Phl·L_Ph/Q_TP`, `φ`.
- **Tier 2** need inhibitor/kinetic data: separate `P_d^eff` (channel) vs `V_max` (carrier);
  influx vs efflux asymmetry.
- **Tier 3** QSPR/measurement (chain-length resolved): `K_prot, K_PL, K_cw, L_Ph, a_R`.
- **Identifiability**: BAF data constrain only the lumped influx conductance
  `g_in = a_R·P_d^eff + carrier clearance` — channel vs carrier are **not separable** from
  BAF alone (need inhibitor experiments).

## 6. Current status
- Derivation + LaTeX docs: complete.
- Python module: runs (BDF stiff solver); reproduces the structural results
  (anion exclusion, terminal-sink accumulation, binding).
- **KNOWN ISSUE**: with placeholder parameters, leaf & grain over-accumulate
  (terminal-sink runaway). Reproducing the empirical **root > straw > grain** ordering is
  the Tier-1 calibration target; it requires genuinely low-translocation parameters and/or
  realistic terminal-tissue loss terms. **Demo BAFs are NOT calibrated** — ignore absolute
  values until calibration.

## 7. Build & run
- Python: `pip install -r requirements.txt` then `python src/pfas_rice_plant_module.py`
  (prints N, B_k, final tissue concentrations/BAFs; saves `pfas_rice_demo.png`).
- FORTRAN (Method B): init submodule (`git submodule update --init`), then follow
  https://phydrus.readthedocs.io/en/latest/getting_started/compilation.html
  (gfortran + `makefile` / `make.bat`).

## 8. Conventions
- Units: time **day**; aqueous conc **µg/L**; tissue conc **µg/kg**; mass **kg**;
  flow **L/day**; `B_k` in **L/kg** (`C_k = B_k · C_w,k`).
- `B_k` has **no density factor** (the `ρ_k` in the current `.tex` is dimensionally wrong;
  fix it on the next doc edit).
- Symbols map 1:1 to `docs/pfas_rice_compartmental_model.tex` (`j_R, B_k, N, L_Ph, ...`).

## 9. Next tasks (prioritized)
1. **Physical realism of terminal compartments** — add loss/limiting terms or constrain
   translocation so leaf/grain don't run away; verify against `root > straw > grain`.
2. **Tier-3 QSPR** for `K_prot`, `K_PL` (chain-length descriptors) to populate `B_k`.
3. **Plug real Phydrus output** into `PlantInputs`; add Freundlich soil sorption (paddy).
4. **Tier-1 calibration** vs rice BAF / lysimeter data (chain-length trends).
5. (Later) **Method B** tight coupling in `external/hydrus_source`.

## 10. Gotchas / external dependencies
- DPU module source is **not public** (author request only). The ionizable extension
  (Brunetti 2022) is **not in the HYDRUS distribution** → reimplement from the papers.
- `phydrus/source_code` is HYDRUS-1D **4.08** (older than official 4.17), **LGPL-3.0**, and
  is the **base soil engine only** (no DPU/PFAS/ionizable modules).
- Key references: Brunetti 2019 *WRR* `10.1029/2019WR025432`; 2021 *ES&T*
  `10.1021/acs.est.0c07420`; 2022 *J. Hazard. Mater.* `10.1016/j.jhazmat.2021.127008`.
