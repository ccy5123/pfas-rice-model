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

## Two audiences: Simple (default) vs Expert

The app opens in a **Simple mode** for a general audience (policy makers, students, the
public): a plain-language intro, a friendly chemical dropdown + a low/medium/high
**contamination preset**, three plain metric cards + a one-sentence summary, and five
jargon-free tabs (🗺️ *Where it goes* · 📈 *Build-up over time* · 📊 *How much builds up* ·
🔎 *Work backwards* · ℹ️ *About & glossary*) — no `BAF`/`Cwᵒ`/`f_xy`/`eᴺ` symbols in the
default view. A prominent **research/educational disclaimer** sits on every screen (top +
footer), and a **Download** expander exports the results (CSV; PNG when `kaleido` is present).

**🔎 Work backwards (Bayesian inverse).** Already have a lab result for the rice? Enter the
measured tissue concentrations (root/straw/grain) and the app estimates **how contaminated
the soil water most likely was, with a 95% credible interval** — a Bayesian parameter
estimation (`model_api.estimate_exposure_bayesian`; a Laplace posterior in log Cwᵒ that
inverts the model's *saturable* uptake, so it is a real nonlinear inverse, not a division).
The Expert tab adds the identifiability caveat (only the exposure level is well-posed from
tissue data; pinning transport needs an independent sap/pore-water measurement).

**📋 Your own data tables (growth + Cwᵒ).** Both modes can take a user-entered **growth
table** (organ FRESH-weight mass over time) and a **time-varying pore-water Cwᵒ(t)** (absolute
µg/L) as editable grids (`st.data_editor`) or a CSV upload — Simple via the "Use my own data
tables" checkbox, Expert via the "Custom tables (Cwᵒ + growth)" data source
(`model_api.drivers_from_tables`). Each compartment has an editable **density** ρ_k [kg/L
fresh] (default root 1.0 / stem 0.30 / leaf 0.30 / grain 1.20) that links the entered weight
to tissue volume (rice leaf/culm hold aerenchyma ⇒ < 1; grain is denser ⇒ > 1). The growth
table is fresh-weight **mass** (the model's M unit) and the transport ODE is mass-based, so
ρ_k is the explicit mass↔volume bridge (and per-volume reporting), not a transport prefactor.
Either table can be left at its default to fall back to the built-in value.

Flip the sidebar **🔬 Expert / advanced controls** toggle to restore the full research
interface documented below — the six exposure modes, SMILES structure input, the opt-in
mechanism switches, every model parameter, and all ten tabs. Nothing is removed; it is
layered behind the toggle.

The compute is UI-agnostic and unit-tested head-less:

| Layer | File | Role |
|---|---|---|
| Model API | `src/model_api.py` | `simulate(...)`, driver/soil/biomonitoring helpers, colormap series |
| Plotly figures | `src/plots.py` | `fig_plant_schematic`, `fig_schematic_*`, `fig_soil_profile`, `fig_drivers`, `fig_isotherm`, … |
| UI | `app.py` | Streamlit widgets + tabs (no science) |

---

## Mechanism switches (sidebar → ⚙️ Mechanism, Expert only)

Two of the model's open questions are *named modes*, not buried constants — the repo's idiom
(`lipid_source`, `f_xy_source`, `cwo_profile`, `biomass`, `uptake`). Both default to the
shipped model and, when either is on, a banner sits next to the headline metrics, because the
sidebar expander is collapsed and the run would otherwise look like the shipped one.

- **Root entry — carrier (default) vs apoplastic bypass** (`simulate(uptake=…)`,
  `model_api.UPTAKE_MODES`). The Michaelis–Menten carrier is this repo's one addition to the
  Trapp cell model, and it keeps its place **by default, not by evidence**: on Yamazaki the
  carrier (log10 RMSE 1.035) and one global bypass `g_apo`=20 (0.996) are indistinguishable
  (bootstrap 0.749), while adding nothing at all fails (2.640) — so *an* addition is necessary
  and *which* is undecided (`validation/carrier_vs_bypass.py`). An "Override the entry
  constants" checkbox exposes `Vmax ×`, `Km ×` and `g_apo` — the levers of that scan and of
  the dose series (whose durable result is a POST-HOC bound: the carrier must have
  `Km ≥ 500 µg/L`, 100× the fitted 5, to be as flat as Tang's 5 doses measure).
- **Lipid-facilitated loading** (`simulate(lipid_loading=True)`) — the K_PL-gated,
  B-independent `g_xy·C` / `g_ph·C` term. It is the repo's strongest cross-dataset result
  (Tang per-organ 1.232 → 0.516 with **nothing re-fit to Tang**; best of the variants on Kim
  2019 grain), but it supplies its own `f_xy`, so the sidebar's `f_xy` choice is ignored while
  it is on, and against the *Yamazaki* bars it is in-sample. EXPLORATORY, default off.

Both switches propagate to every tab, the congener comparison **and** the Bayesian inverse
(`estimate_exposure_bayesian(**sim_kw)`), so the inverse is never run under a different model
than the forward tabs. Neither touches `parameters.json`.

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

## Neutral organics tab (the Briggs/Kow DPU base)

The **🧪 Neutral organics** tab (Expert only) runs `model_api.simulate_neutral(log_kow, …)` —
the framework's neutral base on the *same* 4-compartment ODE with `z = 0`, so the GHK factor
→ 1, anion exclusion `eᴺ` falls 107 → 1, the membrane term degenerates exactly to passive
diffusion and the carrier is off. Inputs: **log Kow** (the one required input — a neutral
compound has no congener), compound name, in-planta half-life, the TSCF QSPR (Briggs 1982 vs
the broader Schriever 2020 refit), the **root-lipid reading** (`lipid_source`: `measured` vs
`briggs_anchor`), an opt-in phloem toggle, and opt-in **plant–air exchange**
(`src/plant_air.py`; needs MW and `K_AW`, identically zero at `K_AW = 0`). It reports TSCF,
`K_PW`, the three tissue BAFs and the tissue-dynamics curve.

**Weak electrolytes and the apoplastic bypass** (⚗️ expander). A `pKa` turns the compound into
one that is a neutral molecule *and* an ion at once — which the `z = 0` trick cannot express,
since one valence must be either 0 or −1 — so it routes through the `(f_n, f_d)` speciation
pair instead (`simulate_neutral(pKa=, is_acid=, pH=)`); the panel prints the resulting `f_n`
and says which way the ion is pushed (an anion is excluded by the inside-negative membrane, a
cation is *attracted* — the same pKa gives a ~19× higher root BAF as a base than as an acid).
`g_apo` adds the apoplastic bypass, a route *around* the membrane that feels neither
speciation nor GHK. How far this is tested is stated in-UI: on Schriever 2020's 67 ionisable
rows the **direction is supported** (Spearman(f_n, TSCF) +0.480; speciation lifts the model's
rank correlation +0.284 → +0.520) and the **magnitude is refuted** (bias −0.203; ~nothing
predicted below `f_n ≈ 1e−3`), and a small `g_apo ≈ 0.5–1` improves rank *and* scale while the
RMSE-optimal 5 wrecks the ordering — so it must not be fitted on RMSE. Nothing is adopted:
`pKa=None` and `g_apo=0` are the defaults and are bit-identical to the strictly neutral path.

**Expert only, on purpose.** The Simple view is congener-driven and symbol-free; a neutral
compound is described by a log Kow, so it does not belong there — and the neutral **grain
compartment has never been tested against data**, which is the opposite of what a
general-audience screen should show absolute numbers for.

**Nothing in it is fitted**: `K_PW` and `TSCF` both follow from log Kow via published QSPRs, so
this is the one place in the app where the DPU *backbone* is exposed without the fitted PFAS
transport behind it (a-priori vs measured rice: log10 RMSE **0.281** root partition, **0.783**
per-organ TF, on the season-ODE basis; on the equilibrium basis appropriate to a root-partition
measurement Liu is **0.206**, the repo's best a-priori result). The tab states the standing scope
limits in-UI: grain untested; the **root** lipid (1 % fw) corroborated by measured *cereal* roots
(Li 2019) but **stem/leaf still Trapp 1994's soybean 3 %**; the Ge leaf residual confounded with
the missing half-life; and, for **lipophilic** compounds, a root-partition question that is now
narrower than it was — Kodešová 2019 measures carbamazepine's root partition directly (n=21,
median `RCF_fw` 1.10), Briggs lands within ~1.5× of it and **Brunetti 2021's calibrated pea
`K_RW` = 13.3 sits ~12× above the same compound's measurement**, so that disagreement is in the
calibration, not the partition core. What survives is confined to log Kow ≳ 3, where the
measurements contradict each other (propiconazole RCF 43.65 in lettuce vs 9.32 in *rice*).
Records: `docs/neutral_dpu_validation.md`.

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

A **lipid loading (OOS)** checkbox adds a third model series
(`tang_tf_validation(lipid_loading=True)` → `plots.fig_tang_tf(val, val_refit, val_extra)`).
It is the honest counterpart of the refit bar: the green refit was calibrated *on* these
measurements, while the lipid mechanism's constants were fit on **Yamazaki** and transferred
untouched — yet it lands in the same place (0.516 vs 0.519 across Tang's three congeners) and
fixes the dominant free-anion failure at the mechanism level (PFOS stalk 0.013 → 0.620 vs Tang
0.571). Residuals it does not fix: GenX (provisional ether `f_xy`) and the PFOS endosperm
(~5× under). See `validation/oos_tang_lipid.py`.

---

## BAF vs observed — two-pool overlays

The **📊 BAF vs observed** tab compares the model's root/straw/grain BAF to the fixed Yamazaki
2023 bars. For a **curated congener** it can additionally overlay the EXPLORATORY **sequestration
two-pool** model (checkbox, on by default) — `model_api.simulate_twopool_seq` → `plots.fig_baf(res, obs, extra=…)`.
The two-pool seq model adds an irreversible non-K_PL `k_seq` root sink (mobile + sequestered pools)
that captures the long-chain **root** BAF and the **PFOS/PFUnDA** split the single-pool 4-pool core
misses, while keeping the monotone physical `f_xy` (overall log10 RMSE 0.251; `docs/twopool_root_exploration.md`).

A second checkbox overlays the **structural merge** (`model_api.simulate_twopool_nstem`): the same
two-pool sequestration root driving the *redistributed* N-stem+leaf shoot instead of the 4-pool
pass-through stem. Swapping the whole shoot costs only ~0.04 log on Yamazaki (0.301 vs 0.278) with
the root RMSE unchanged — the root mechanism is separable from the shoot model — and it is what makes
a per-organ Tang test fair (that OOS moves 1.398 → 0.801, the recovery carried by the diagnosed
stalk; `docs/twopool_root_exploration.md` Result 8).

Caveats made explicit in the tab:
- **EXPLORATORY / in-sample** (Yamazaki fit); opt-in. `parameters.json` and the canonical 4-pool core
  are unchanged — the overlay does not alter the model bar.
- The overlay runs at the two-pool's **calibrated operating point** (Cwᵒ=1, season≈120, demo forcings),
  so it is comparable to the fixed Yamazaki bars but does **not** track the sidebar `f_xy`/Cwᵒ/biomass
  (the 4-pool core bar does). Curated congeners only (the cached fit has no record for a novel SMILES).
- The **carrier** two-pool (`close_longchain_2pool`, a saturated DOF-0 long-chain closure) is API-only —
  it reproduces the observed bars by construction and is too slow (~1 min/congener) to render live.

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
~+40–70%); pass `biomass="growth_rice"` to reproduce the legacy artifacts (`reproduce_demo.py` and
`calibration.py` use their own drivers and are unaffected). Both the app and the code-level
`simulate(biomass=)` default are **ORYZA2000**.

---

## Six ways to drive the model

Only the pore-water free-anion concentration `Cwᵒ(t)` is PFAS-specific. The transpiration
stream `Q_TP(t)` and organ masses `M(t)` are crop physiology (measured FAO-56 transpiration +
ORYZA IR72 biomass) and are reused across modes unless you supply your own.

| Mode | `Cwᵒ(t)` source | `Q_TP`, `M(t)` | API entry |
|---|---|---|---|
| **Model (parametric)** | a constant you set (flat, or the analytic flooded shape) | measured / placeholder | `simulate(congener, Cwo=…, cwo_profile=…)` |
| **Custom tables** | a `Cwᵒ(t)` table you type or paste | your own growth table, or measured | `drivers_from_tables(growth, cwo)` → `simulate(drivers=…)` |
| **HYDRUS / CSV drivers** | a HYDRUS-1D / Phydrus run | from the CSV, or measured | `simulate(congener, drivers=load_driver_csv(...))` |
| **Run HYDRUS-1D (live)** | a real HYDRUS-1D run executed in-app | HYDRUS root uptake + ORYZA | `hydrus_drivers(congener, …)` → `simulate(drivers=…)` |
| **Soil inventory** | inverting a soil load (Freundlich) | measured | `pore_water_from_inventory(...)` → `drivers_from_arrays` |
| **Biomonitoring** | a measured pore-water value | — (not needed) | `baf_from_measurement(conc, Cwo)` |

### Live HYDRUS-1D run (`src/soil_hydrus.py`)
The "Run HYDRUS-1D (live)" mode executes the **genuine HYDRUS-1D engine** (built from the
`external/hydrus_source`, **vendored in the repo** — no submodule to initialise) through `phydrus`: a one-season paddy model — Richards
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
