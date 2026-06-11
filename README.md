# PFAS–Rice 4-Compartment Uptake Model — Reproduction Package

Mechanistic uptake model for **permanently-anionic PFAS in paddy rice** (*Oryza sativa*):
an IOC (ionizable organic compound) extension of the DPU / Trapp framework over four
compartments (root → stem → leaf → grain). Covers **12 congeners**: PFCA C4–C12
(PFBA…PFDoDA) and PFSA C4/C6/C8 (PFBS, PFHxS, PFOS).

This package consolidates the two closed open-parameter workstreams — **GAP A (cell-wall
partition `K_cw`)** and **GAP B (root→shoot loading `f_xy`)** — plus the model code, the
basis-A binding factors, and the S6 validation, into one reproducible bundle.

---

## Quickstart

```bash
pip install -r requirements.txt           # numpy, scipy, matplotlib
python build_parameters.py                # (re)build params/parameters.json from sources
python reproduce_demo.py                  # Yamazaki BAF via full ODE (W2 fit; log10 RMSE ≈ 0.029)
python reproduce_demo.py --rec            # monotone physical f_xy (single-straw mismatch — see note)
python src/literature_params.py           # literature QSPRs (K_PL/K_prot/Koc/f_d) + Kim2019 L_Ph fit
python validation/nstem_gradient_check.py # multi-height stem: reproduces the Yamazaki stem gradient
pip install pytest && pytest              # tests (structure, mass conservation, QSPRs, calibration, API)
```

`reproduce_demo.py` loads `params/parameters.json` + `src/` and runs the 4-compartment ODE
for all 12 congeners, printing predicted vs observed root/straw/grain BAF.

## Interactive app
A Streamlit dashboard to run the model and see the results (tissue concentrations/BAFs,
binding factors, chain-length parameter trends, and the measured forcings):

```bash
pip install -r requirements.txt -r requirements-app.txt
streamlit run app.py
```

Pick a congener and scenario (pore-water `C_wᵒ`, membrane potential `E_m`, `f_xy`
source, measured vs placeholder forcings). Compute lives in `src/model_api.py`
(`simulate(...)`), which is UI-agnostic and used by the app and the tests.

---

## Layout

```
pfas_rice_model/
├── README.md                     ← you are here
├── build_parameters.py           ← assembles params/parameters.json from source tables
├── reproduce_demo.py             ← self-contained ODE reproduction (entry point)
├── src/                          ← model code
│   ├── pfas_rice_plant_module_4pool.py        4-compartment plant ODE (basis-A)   ← CANONICAL
│   ├── pfas_rice_plant_module_4pool_surf.py   + K_surf (Fe/Mn-plaque surface pool)
│   ├── pfas_rice_plant_module_5pool.py        + explicit lignin pool
│   ├── pfas_rice_plant_module_nstem.py        N serial stem segments (multi-height; GAP-B fix)
│   ├── pfas_rice_plant_module.py              import alias → 4pool_surf (do not delete)
│   ├── soil_paddy.py                          soil↔porewater (Freundlich)         ← legacy redox sign
│   ├── soil_paddy_redox_corrected.py          W3-CORRECTED redox (USE THIS)
│   ├── calibration.py                         BAF→parameter fitting machinery
│   └── literature_params.py                   literature QSPRs/anchors (cited) + Kim2019 BAF data
├── params/                       ← parameters
│   ├── parameters.json           ★ CANONICAL consolidated parameter set
│   ├── f_xy_recommended.csv      ★ GAP B f_xy(n) (monotone) vs W2-fit
│   ├── S6_Bk_basisA_allorgan.csv basis-A B_k(n) all organs (supersedes Bk_table_S5)
│   ├── Kcw_Klignin_params_v2.csv GAP A source (K_cw poly/lignin + whole-cw per organ)
│   ├── rice_tissue_params.csv    tissue composition (θ, f_prot, f_PL, f_cw)
│   ├── Bk_table_S5.csv           legacy naive-basis assembly (K_PL/K_prot/f_xy source only)
│   └── W2_transport_fit.csv      S6 transport fit (f_xy/L_Ph/kappa_d per congener)
├── data_obs/                     ← observed BAF/TF for validation
│   ├── obs_baf_Yamazaki.csv      Andosol, clean per-congener water (main calibration)
│   ├── obs_baf_Li2025.csv        Tianjin field (group-water; see surface caveat)
│   ├── Li2025_BAF_TF.csv         Li2025 BAF + TF summary
│   └── yamazaki_stem_height.csv  Yamazaki S18/S19 per-height stem gradient (for nstem)
├── validation/                   ← reproduction scripts + outputs + figures
│   ├── S6_alphaQC1_basisA.py     membrane-share / α identifiability
│   ├── S6_surface_crossfield.py  surface-excess (water-quality confound)
│   ├── S6_Gap4.py                full-ODE reproduction + cross-field TF
│   ├── nstem_gradient_check.py   multi-height stem: stem-gradient direction vs Yamazaki
│   └── figures/*.png
└── docs/
    ├── DELIVERABLE_GAP_A_Kcw.md  GAP A verdict, recommended values, experimental design
    ├── DELIVERABLE_GAP_B_fxy.md  GAP B verdict, f_xy(n), F1/F4, the W2↔theory reconciliation
    ├── theory_anchor.tex         GAP B theory (Trapp+Briggs → monotone f_xy; IOC pH term)
    ├── H8_handoff_S6_final.md     S6 session handoff (validation + open items)
    ├── sources.csv               DOI shortlist
    └── literature_db/            curated empirical parameter DB (xlsx + CSV) + raw_si/ extractions
```

---

## Model in one screen

Binding factor (**basis A, fresh weight** — the single most important convention):

```
B_k = θ_fw + (1 − θ_fw) · ( f_prot·K_prot + f_PL·K_PL + f_cw·K_cw )      [L/kg fw]
```

θ_fw = tissue water fraction; f_* = **dry-weight** mass fractions; K_* = partition coeffs
[L/kg pool-dw]. The `(1 − θ_fw)` factor is mandatory — omitting it (the legacy
`Bk_table_S5.csv`) over-states B_k ~3× and corrupts pool shares. `f_cw` is the whole cell
wall (polysaccharide + lignin); the matching K is `K_cw_wholecw_<organ>`.

Root influx = GHK electrodiffusion (anion exclusion e^N ≈ 107 at E_m = −120 mV) + saturable
Michaelis–Menten carrier (overcomes exclusion). Xylem loading = `f_xy · Cw_root`; grain is
phloem-dominated (`L_Ph · Cw_leaf`).

---

## ⚠ The one thing to get right: which `f_xy`

`parameters.json` carries **two** `f_xy(n)`:

| field | what it is | use it for |
|---|---|---|
| **`f_xy_recommended`** | **monotone** physical TSCF (theory-derived, cross-field-validated). C4 0.79 → C12 0.003. | citing/reporting the parameter; the physically-correct value |
| `f_xy_W2fit` | transport fit to Yamazaki; **rises spuriously for C10+** (0.08→0.67) | reproducing Yamazaki through the *current* ODE structure only |

The two diverge at long chains because the W2 fit (saturated, 3 param / 3 obs) absorbs an
**unmodeled stem accumulation gradient** into `f_xy`. Theory (Trapp GHK + Briggs LFER) and the
**water-independent cross-field TF** require a monotone *direction* → `f_xy_W2fit` long-chain rise
is largely a single-compartment artifact. **Use `f_xy_recommended`** for the parameter. See
`docs/DELIVERABLE_GAP_B_fxy.md`.

**Update — the multi-height stem fix is now implemented** (`src/pfas_rice_plant_module_nstem.py`,
`validation/nstem_gradient_check.py`). Resolving the stem into N serial segments (transpiration
draw-off + radial exchange + growth dilution) lets a **monotone f_xy reproduce the observed
Yamazaki stem gradient for the PFCAs** (short chains concentrate upward, long chains flat/down;
the flip is set by `B* ~ Q_s/(M_s·μ_s)`). Caveats from the review: (i) the *absolute* crossover
and f_xy scale need **measured `Q_TP(t)`/`M_s(t)`** (placeholder transpiration is ~5× too high);
(ii) **PFOS/PFSA** translocate upward despite high binding — a binding-driven monotone f_xy misses
them, so a PFSA-specific transport term is still needed.

---

## Status

> Honest status after review + the multi-height-stem work (the README's earlier
> "CLOSED/validated" labels are tightened here):

- **Binding `B_k`** — built on **measured** per-congener `K_PL` (Chen 2025 K_MW, vs Droge 2019) and
  `K_prot` (Zhou 2025 dialysis `K_prow`; soy = plant, BSA = animal) — see `docs/literature_db/raw_si/`
  and `src/literature_params.py`. basis-A fresh-weight convention.
- **GAP A (K_cw)** — values delivered, but **anchored (DFT ladder + measured lignin), not directly
  measured** — the long-term weakest point; `K_cw` is also a minor pool for the membrane-dominated
  long chains. `docs/DELIVERABLE_GAP_A_Kcw.md`.
- **GAP B (f_xy)** — *shape* resolved (monotone; short-chain ceiling ≈0.8 anchored to Felizeter TSCF)
  and the **multi-height stem reproduces the PFCA stem gradient with a monotone f_xy**. **Not fully
  closed**: the absolute scale/crossover needs measured `Q_TP(t)`/`M_s(t)` (task 2) and PFSA needs a
  separate transport term (task 3). `docs/DELIVERABLE_GAP_B_fxy.md`.
- **Validation caveat** — `reproduce_demo.py`'s log10 RMSE 0.029 uses the **saturated W2 fit**
  (3 transport params per congener fit to 3 observed BAFs → reproduction is guaranteed, *not* a
  predictive test). The genuine out-of-sample evidence is the water-independent **cross-field TF**
  (monotone direction) and the **nstem gradient direction** (PFCAs). `docs/H8_handoff_S6_final.md`.
- **Tier-1 fit** — `src/literature_params.py` fits `L_Ph` to the Kim 2019 PFOA grain BAF (matches 4.43 L/kg).
- **Tests** — 52 passing (`pytest`).

**Open (data-limited, not modeling work):** rice (not wheat) per-congener root subcellular →
membrane-share/α; reliable per-congener pore-water or hydroponic RCF → surface test + f_xy
absolute scale; measured Q_TP(t), M(t) → f_xy absolute; direct K_cw_poly + rice cw monosaccharide
composition; in-situ paddy E_m.

**Open (modeling, doable now):** (1) multi-height stem compartment so the *physical* monotone f_xy
reproduces long-chain straw (currently the W2 fit compensates); (2) integrated soil→plant run with
`soil_paddy_redox_corrected` + a realistic flooding schedule; (3) f_PL (0.01–0.02, 2× uncertain)
uncertainty propagation.

**Config to standardize:** root θ = 0.90 (measured 0.90–0.92), root f_PL = 0.015, grain θ stage-
dependent (0.14 harvest / 0.30 filling). Conclusions are robust to these.
