# PFAS–Rice 4-Compartment Uptake Model — Reproduction Package

Mechanistic uptake model for **permanently-anionic PFAS in paddy rice** (*Oryza sativa*):
an IOC (ionizable organic compound) extension of the DPU / Trapp framework over four
compartments (root → stem → leaf → grain). Covers **12 calibrated congeners**: PFCA C4–C12
(PFBA…PFDoDA) and PFSA C4/C6/C8 (PFBS, PFHxS, PFOS), **plus GenX** (HFPO-DA, an ether-PFAS —
provisional, added for the Tang 2026 validation; `docs/VALIDATION_TANG2026_KR.md`).

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

### Parameterise any PFAS from its structure (SMILES)

`src/pfas_structure.py` maps a **chemical structure** onto the model, so it runs for *any*
PFAS, not only the curated 13. RDKit parses the SMILES → descriptors → a `Compound` by
measured read-across (known congener → curated params) or a fragment QSPR (novel → provisional):

```bash
pip install -r requirements-structure.txt        # RDKit (optional)
python src/pfas_structure.py                      # SMILES -> descriptors -> Compound demo
```
```python
import model_api as api
api.simulate_from_smiles("OC(=O)" + "C(F)(F)"*11 + "C(F)(F)F")   # a novel C13 PFCA, run from SMILES
```
See `docs/structure_input.md` (binding/speciation come from structure; `f_xy` ordering only — absolute scale stays fitted).

### Real soil side — HYDRUS-1D (Method A, now wired)

The soil half of Method A runs the genuine HYDRUS-1D engine (compiled from the
`external/hydrus_source` submodule, driven via `phydrus`) to produce a real
per-congener pore-water trajectory `C_w^o(t)` that drives the plant ODE:

```bash
git submodule update --init external/hydrus_source
cp external/hydrus_source/makefile external/hydrus_source/source/
(cd external/hydrus_source/source && make)     # needs gfortran
pip install phydrus
python src/soil_hydrus.py                       # per-congener pore-water summary
python validation/hydrus_coupled_run.py         # full soil→plant run + figure/CSV
```

Weakly-sorbed short chains (Kd≈0.01–0.15 from `Koc·f_oc`) leach to near-zero during
flooding, so the constant-`Cwo` placeholder **over-predicts grain/straw BAF ~2–4×**
(PFBA grain 2.07→0.43); strongly-sorbed long chains (Kd≳7) stay buffered. The HYDRUS
tests auto-skip where the executable/`phydrus` is unavailable.

`reproduce_demo.py` loads `params/parameters.json` + `src/` and runs the 4-compartment ODE
for all 12 congeners, printing predicted vs observed root/straw/grain BAF.

## Interactive app — the visualization tool
A Streamlit dashboard that **draws the soil + rice plant to scale and colours each
compartment by its PFAS accumulation** (a heat colormap you can scrub through the season),
alongside interactive Plotly time series (hover, zoom, legend-toggle):

```bash
pip install -r requirements.txt     # full app stack (numpy/scipy/matplotlib + streamlit/plotly/pandas)
streamlit run app.py
```

**Deploy to the web (Streamlit Community Cloud — free, share by URL):** push this repo to
GitHub, then on [share.streamlit.io](https://share.streamlit.io) pick the repo, branch, and
`app.py` → Deploy. `requirements.txt` is the full app stack, so it deploys as-is (RDKit/HYDRUS
are optional and the app hides the live-HYDRUS mode when the engine is absent). `requirements-app.txt`
holds only optional extras (kaleido PNG export, phydrus). See `docs/deploy.md`.

**🗺️ Plant & soil map** — the headline view: a fibrous-rooted rice plant with arching
culms, long leaves and drooping grain panicles, each organ filled on a shared colorbar by
its concentration (or BAF). A day slider (or ▶ animate) shows *where and when* PFAS builds
up — leaf is xylem-terminal, grain is phloem-fed, the root retains the anion.

**Five input modes** (sidebar “Data source”) cover the whole exposure space:

| Mode | Pore-water `Cwᵒ(t)` from | When |
|---|---|---|
| **Model (parametric)** | a constant you set | quick what-ifs / teaching |
| **HYDRUS / CSV drivers** | a HYDRUS-1D/Phydrus run (`t,Cwo,Qtp,M_*` CSV) | you have a calibrated soil model |
| **Run HYDRUS-1D (live)** | the real HYDRUS engine, executed in-app | you want HYDRUS to run here (needs it built) |
| **Soil inventory** | inverting a total soil load (Freundlich) | you know soil PFAS, not pore water |
| **Biomonitoring** | a measured pore-water value (no HYDRUS) | you have field tissue + water data |

The **live HYDRUS-1D** mode runs the genuine engine (built from `external/hydrus_source` via
`phydrus`) for a one-season paddy model → congener-dependent `Cwᵒ(t)` (short chains leach, long
chains buffer); it auto-detects the engine and shows build steps if absent. See
`src/soil_hydrus.py` and `docs/visualization_tool.md`.

Other tabs: tissue dynamics, **soil & drivers** (`Cwᵒ(t)`, `Q_TP(t)`, `M(t)`, Freundlich
isotherm, depth profile), BAF vs observed/measured, chain-length trends, compare congeners,
and an **About** tab documenting the HYDRUS-1D input/output mapping and the biomonitoring
path. Compute is in `src/model_api.py` (`simulate(...)`, soil/driver/biomonitoring helpers);
the Plotly figures in `src/plots.py` (`fig_plant_schematic`, …) — both UI-agnostic and
covered by the tests. Ready-to-load examples are in `examples/`. Full guide:
`docs/visualization_tool.md`.

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
│   ├── soil_hydrus.py                         REAL HYDRUS-1D → C_w^o(t),Q_TP(t) via phydrus (Method A)
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
│   ├── hydrus_coupled_run.py     REAL HYDRUS-1D soil → plant coupling (vs constant-Cwo baseline)
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
  **Full validation write-up (KR) with figure: `docs/VALIDATION_KR.md`** — calibration (Yamazaki,
  reproduction) vs out-of-sample (Kim 2019 grain; Li 2025 TF, inconclusive); run
  `python validation/validation_summary.py` for the figure.
- **Independent rigor audit (sci-adk)** — seven adversarial runs adjudicate the model by frozen
  criteria (`sci_adk_review/`): formal foundations SUPPORTED (mass balance, anion exclusion, soil
  QSPR spread, SMILES read-across), naive predictive claims REFUTED (0.029 is in-sample), structural
  adequacy SUPPORTED on the ORYZA2000 biomass, and a **lipid-facilitated loading mechanism that
  generalizes out-of-sample** across two independent datasets (Tang 2026, Kim 2019). Consolidated,
  citable manuscript: **`docs/sci_adk_rigor_review.tex`** (KR narrative: `sci_adk_review/FINDINGS.md`).
- **Tier-1 fit** — `src/literature_params.py` fits `L_Ph` to the Kim 2019 PFOA grain BAF (matches 4.43 L/kg).
- **Visualization tool** — `app.py` (+ `src/model_api.py`, `src/plots.py`): plant/soil accumulation
  colormap + five exposure modes (model / HYDRUS CSV / **live HYDRUS-1D** / soil inventory /
  biomonitoring). `docs/visualization_tool.md`.
- **Soil side (Method A) — real HYDRUS-1D** (`src/soil_hydrus.py`): the engine is compiled and wired
  (built from `external/hydrus_source` via `phydrus`); per-congener `C_w^o(t)` from the soil transport
  solve drives the plant ODE (and the app's live mode). Short chains leach → the constant-`Cwo`
  placeholder over-predicts grain/straw BAF ~2–4×; long chains stay buffered. Tests skip gracefully
  when the engine isn't built.
- **Opt-in lipid-facilitated loading** — `simulate(lipid_loading=True)` adds a K_PL-gated, B-independent
  xylem/phloem term (`g_xy`/`g_ph`; default 0) so high-binding long chains reach the shoot; in-sample/
  exploratory (`docs/fxy_longchain_lipid_exploration.md`).
- **Tests** — 138 passing, 4 skipped (`pytest`; 142 collected — the 4 skips are the HYDRUS engine tests, which run when the engine is built).

**Open (data-limited, not modeling work):** rice (not wheat) per-congener root subcellular →
membrane-share/α; reliable per-congener pore-water or hydroponic RCF → surface test + f_xy
absolute scale; measured Q_TP(t), M(t) → f_xy absolute; direct K_cw_poly + rice cw monosaccharide
composition; in-situ paddy E_m.

**Open (modeling, doable now):** (1) multi-height stem compartment so the *physical* monotone f_xy
reproduces long-chain straw (currently the W2 fit compensates); (2) ~~integrated soil→plant run~~
**DONE** via the real HYDRUS-1D engine (`src/soil_hydrus.py`, `validation/hydrus_coupled_run.py`) —
remaining: a measured field flooding schedule, anoxic/flooded sorption, and the user's site soil;
(3) f_PL (0.01–0.02, 2× uncertain) uncertainty propagation.

**Config to standardize:** root θ = 0.90 (measured 0.90–0.92), root f_PL = 0.015, grain θ stage-
dependent (0.14 harvest / 0.30 filling). Conclusions are robust to these.
