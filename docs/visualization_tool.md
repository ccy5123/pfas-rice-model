# PFAS–Rice visualization tool (`app.py`)

An interactive Streamlit + Plotly dashboard that makes the soil + plant model **visible**:
it draws the paddy soil and a rice plant to scale, colours each compartment by how much
PFAS has accumulated in it (a heat colormap), and lets you scrub the season to watch the
build-up. It also covers the **whole exposure-input space** — from a full HYDRUS-1D soil
run down to bare biomonitoring data with no soil model at all.

```bash
pip install -r requirements.txt -r requirements-app.txt
streamlit run app.py
```

The compute is UI-agnostic and unit-tested head-less:

| Layer | File | Role |
|---|---|---|
| Model API | `src/model_api.py` | `simulate(...)`, driver/soil/biomonitoring helpers, colormap series |
| Plotly figures | `src/plots.py` | `fig_plant_schematic`, `fig_schematic_*`, `fig_soil_profile`, `fig_drivers`, `fig_isotherm`, … |
| UI | `app.py` | Streamlit widgets + tabs (no science) |

---

## The plant & soil accumulation map

`plots.fig_plant_schematic(values, cmin, cmax, label, Cwo)` draws:

- a **paddy soil box** with a water-table line and the pore-water `Cwᵒ` annotated;
- a **fibrous root mass** in the soil, **arching culms/tillers**, **long leaf blades**, and two
  **drooping grain panicles** — the silhouette of ripening *Oryza sativa*;
- every organ filled with the colour sampled from a shared colorbar at that compartment's
  value, so you can read *where* PFAS concentrates at a glance.

`fig_schematic_from_res(res, metric, t_index)` builds it straight from a `simulate(...)`
result; the colour limits span all organs over the whole season (set by
`model_api.metric_series`), so colours stay comparable while you drag the **day slider** or
press **▶ animate** (`fig_schematic_animated`). Switch the colorbar between **concentration**
(µg/kg) and **BAF** (L/kg).

For **biomonitoring** input the same figure is fed the measured tissue values directly; when
only `root/straw/grain` are reported (no separate stem/leaf), the whole shoot takes the straw
colour.

---

## Tang 2026 validation tab (out-of-sample)

The **✅ Tang TF (OOS)** tab checks the root→shoot loading `f_xy` against **Tang et al. 2026**
(flooded paddy, Nipponbare, 150 d; PFOA/PFOS/GenX) — the per-organ transfer factor
**TF = C_organ/C_root** (SI Table S8), shown for the selected congener as grouped bars:
**Tang (measured)** vs **model** vs **model with the Tang-refit `f_xy`** (`model_api.tang_tf_validation`
→ `plots.fig_tang_tf`). Only Tang's head-group *sign* went into the build, so the magnitudes are OOS.

Three things the tab makes explicit (the rigor points from this work):
- **Dry-weight basis.** Tang's TF is dry/dry; the model conc is fresh-weight, and the
  `(1−θ)` factor differs by tissue, so `TF_dw = TF_fw·(1−θ_root)/(1−θ_tissue)` is applied
  (comparing fresh model TF to dry Tang TF flatters the grain ~8×).
- **`f_xy` is condition-dependent.** PFOS `f_xy` ≈ 0.14 (Yamazaki, clean water) vs ~0.32
  (Tang, flooded soil); GenX's provisional 0.233 is ~12× too high (refit ≈ 0.013). The refit
  is **override-only** (`parameters.json` unchanged). Dose toggle (across-dose mean vs 0.1 µg/g)
  and an optional ORYZA-biomass driver are exposed.
- **Grain is structurally under-predicted** ~3–8× and is *not* closable by `L_Ph`/lipid
  (`docs/tang2026_grain_units_exploration.md`).

---

## Biomass driver M(t) & the Tissue-dynamics mass plot

The sidebar **“Biomass driver M(t)”** radio chooses the organ-biomass forcing for the built-in scenarios:
- **ORYZA2000 (mechanistic)** — the Level-1 carbon balance (`oryza_growth`): radiation/temperature →
  assimilation → respiration → DVS partitioning. **The app default** (the first-principles choice for this
  mechanistic, HYDRUS-coupled model).
- **growth_rice (partition + logistic)** — ORYZA IR72 partitioning imposed on a logistic total-biomass curve;
  the lightweight reconstruction and the historical **calibration basis**.

The **Tissue dynamics** tab now shows two plots: tissue **concentration** `C_k(t)` [µg/kg] (top, intensive) and
the per-tissue **PFAS mass / burden** `C_k(t)·M_k(t)` [µg/hill] (bottom, `plots.fig_burden`, extensive) — *where
the chemical actually ends up* (a tissue can be high-concentration yet low-mass). The organ **biomass** `M_k(t)`
itself is in the *Soil & drivers* tab (`fig_drivers`). **Caveat (biomass driver):** the
`f_xy` calibration was done on `growth_rice`, so switching to ORYZA2000 shifts BAFs (short-chain straw/grain
~+40–70%); the code-level `simulate(biomass=)` default stays `growth_rice` for reproducibility, while the app
leads with ORYZA2000.

---

## Four ways to drive the model

Only the pore-water free-anion concentration `Cwᵒ(t)` is PFAS-specific. The transpiration
stream `Q_TP(t)` and organ masses `M(t)` are crop physiology (measured FAO-56 transpiration +
ORYZA IR72 biomass) and are reused across modes unless you supply your own.

| Mode | `Cwᵒ(t)` source | `Q_TP`, `M(t)` | API entry |
|---|---|---|---|
| **Model (parametric)** | a constant you set | measured / placeholder | `simulate(congener, Cwo=…)` |
| **HYDRUS / CSV drivers** | a HYDRUS-1D / Phydrus run | from the CSV, or measured | `simulate(congener, drivers=load_driver_csv(...))` |
| **Run HYDRUS-1D (live)** | a real HYDRUS-1D run executed in-app | HYDRUS root uptake + ORYZA | `hydrus_drivers(congener, …)` → `simulate(drivers=…)` |
| **Soil inventory** | inverting a soil load (Freundlich) | measured | `pore_water_from_inventory(...)` → `drivers_from_arrays` |
| **Biomonitoring** | a measured pore-water value | — (not needed) | `baf_from_measurement(conc, Cwo)` |

### Live HYDRUS-1D run (`src/soil_hydrus.py`)
The "Run HYDRUS-1D (live)" mode executes the **genuine HYDRUS-1D engine** (built from the
`external/hydrus_source` submodule) through `phydrus`: a one-season paddy model — Richards
flow + advection–dispersion + **linear Kd** sorption + root water uptake — returns the
**congener-dependent** pore water `Cwᵒ(t)` and the actual root water uptake `Q_TP(t)`. Kd comes
from the Koc(chain-length) QSPR (`literature_params.koc`), so weakly-sorbed short chains leach
under flooding (Cwᵒ falls and rebounds on drainage) while strongly-sorbed long chains stay
buffered (flat) — structure a constant Cwᵒ cannot represent. `Cwᵒ(t)` is normalised to
season-mean `Cwo_ref` so the average exposure matches a constant-Cwo run; only the temporal
shape and congener-to-congener contrast differ.

`model_api.hydrus_available()` gates the UI; when the engine isn't built the app shows the build
steps and stays usable. To enable it:

```bash
git submodule update --init external/hydrus_source
cp external/hydrus_source/makefile external/hydrus_source/source/
(cd external/hydrus_source/source && make)      # needs gfortran
pip install phydrus
```

This is still **Method A** (one-way): HYDRUS computes the soil water+solute, the plant ODE runs
in Python; HYDRUS itself is not modified.

---

## HYDRUS-1D coupling (Method A, one-way) — inputs & outputs

The plant ODE is solved in Python; HYDRUS is **not** modified. HYDRUS-1D (optionally via
Phydrus) supplies the soil-water-solute side; you hand it off as a CSV.

**HYDRUS *inputs* (you set up, soil side):**

- soil-hydraulic parameters (van Genuchten θ_r, θ_s, α, n, K_s);
- the atmospheric boundary condition — precipitation / irrigation / evaporation (paddy
  ponding / drainage schedule);
- the root water-uptake distribution (and a root-depth / rooting-density profile);
- solute-transport parameters — a linear `K_d` **or** the Freundlich `K_F, n`, plus
  dispersivity and (here) negligible degradation;
- the initial and boundary PFAS concentration.

**HYDRUS *outputs* the tool consumes** (map to the driver CSV columns):

| CSV column | HYDRUS-1D source | meaning |
|---|---|---|
| `t` | output times | day after transplant |
| `Cwo` | `Conc` at the root-zone node (`Obs_Node.out` / `solute1.out`) | pore-water free anion [µg/L] |
| `Qtp` | `vRoot` (root water uptake) / `T_act` (`T_Level.out`) | transpiration stream [L/day] |
| `M_root,M_stem,M_leaf,M_grain` | a plant **growth** sub-model (not HYDRUS) | organ fresh mass [kg] |

`Qtp` and the `M_*` columns are **optional**; if omitted the tool fills them from the measured
crop forcings on the same time grid, so a bare `t,Cwo` series is enough to run the plant model.

A depth-resolved solute field (`Nod_Inf.out`: depth × time × `Conc`) can be passed to
`plots.fig_soil_profile(res, profile=...)` for a soil heatmap (the *tight* Method B —
embedding the root-uptake term `j_R` inside the HYDRUS FORTRAN — is future work).

> **Tight coupling note.** This tool implements Method A only. Method B (modifying
> `external/hydrus_source`) would feed the plant uptake back into HYDRUS's solute mass balance;
> it is out of scope here.

---

## Biomonitoring — when HYDRUS is unnecessary

If you already have **measured tissue concentrations** and a **measured pore-water (or
soil-solution) concentration**, the bioaccumulation factor is simply

```
BAF_tissue = C_tissue / Cwᵒ
```

— no transport simulation is required. The Biomonitoring mode reads BAFs straight off the
data, colours the plant map by the measured concentrations, and overlays the model BAF for a
sanity check. Use it for field-survey data, or when a soil model isn't available/needed.

---

## Bundled examples (`examples/`)

| File | Columns | Use |
|---|---|---|
| `hydrus_drivers_example.csv` | `t,Cwo,Qtp,M_root,M_stem,M_leaf,M_grain` | HYDRUS / CSV mode (synthetic HYDRUS-style run) |
| `biomonitoring_example.csv` | `tissue,conc,Cwo` | Biomonitoring mode (Yamazaki 2023 PFOA brown-rice BAFs) |

Both load automatically when no upload is provided, so every mode is demonstrable out of the box.

---

## Static export (optional)

The Plotly figures render in the browser without extra dependencies. To export PNGs
(e.g. for a report), install kaleido and a headless Chrome:

```bash
pip install kaleido && plotly_get_chrome
python -c "import sys; sys.path.insert(0,'src'); import model_api as api, plots; \
          plots.fig_schematic_from_res(api.simulate('PFOA'),'conc',-1).write_image('map.png', scale=2)"
```
