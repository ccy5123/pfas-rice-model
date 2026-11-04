# Web app — FastAPI backend + React front end

A second front end for the same model, built because the Streamlit dashboard
(`app.py`) feels slow. It does **not** change the model: every number still comes
from `src/model_api.py`, every chart from `src/plots.py`, and the one-page summary
from `src/report_html.py` (the same builder the Streamlit download calls).

It runs the **current** model (amendments 1–3): the default configuration of
CLAUDE.md §6.1 (`redistributed_fit` on the Yamazaki 2020 site-year, linear
equilibrium root entry, staged organ water, brown rice + husk pool, dry-basis
harvest operator). It was first written on the pre-amendment-1 API (closed PR #71)
and ported in November 2026; what changed in that port is listed at the end.

## Why the Streamlit app feels slow (and what this changes)

Streamlit re-runs the **whole script on the server on every widget touch**, and it
executes the body of **every tab**, visible or not. The web app splits the work:

| | Streamlit (`app.py`) | Web app (`server/` + `frontend/`) |
|---|---|---|
| Day slider on the plant map | full script re-run | **in the browser** — the season is shipped once, the SVG is re-drawn client-side (`frontend/src/plantSvg.ts`, byte-identical to `src/plot_svg.py`) |
| Switching tabs / colour metric | full script re-run (all tabs computed) | in the browser; a tab fetches **its own** figure only when opened |
| Changing a model input | full re-run | one `/api/run` call after a 300 ms debounce; the previous result stays on screen until the new one arrives |
| Repeat requests | `st.cache_data` per function | a shared, de-duplicating LRU memo per solve / run / figure (`server/memo.py`) — identical concurrent requests wait for one solve |
| Tang check, Yamazaki observable, inverse | cached per session | memoised across users; the inverse streams live progress |
| First visit | cold solves | the Simple view's default + preset scenarios are pre-solved at server start (background warm-up) |

A genuinely new solve costs what the model costs (measured on the port's machine, warm
process: a new congener on the site-year ~0.05 s with the exponential integrator, the
Yamazaki two-harvest observable ~1.6 s, a live HYDRUS-1D run ~8 s); the web app removes
the repeated work around it (plot serialisation, the re-run of every tab). The timing
table of the first version (PR #71) was measured on the pre-amendment-1 model and has
not been re-measured; it is not reproduced here.

## Run it

```bash
pip install -r requirements-web.txt          # = requirements.txt + fastapi/uvicorn/httpx
(cd frontend && npm ci && npm run build)     # once; writes frontend/dist (Node 18+)
uvicorn server.app:app --port 8000           # from the repo root
```

Open <http://localhost:8000> (Simple / Korean) or
<http://localhost:8000/?mode=expert> (Expert / English). Interactive API docs:
<http://localhost:8000/api/docs>.

Front-end development with hot reload: keep `uvicorn` running and, in `frontend/`,
`npm run dev` → <http://localhost:5173> (Vite proxies `/api` to port 8000).

Environment variables: `PFAS_WEB_WARMUP=0` skips the start-up warm-up;
`PFAS_WEB_DIST` points at another build directory; `CTX_API_KEY` enables the EPA
CompTox lookup (as in the Streamlit app).

For a shared server run one uvicorn process per CPU you want to give it, behind
any reverse proxy (`uvicorn server.app:app --host 0.0.0.0 --workers 4`). Each
worker keeps its own memo.

## Layout

```
server/
  app.py       FastAPI routes + serves frontend/dist at "/"
  schemas.py   Scenario = the whole sidebar as one JSON document (pydantic)
  service.py   Scenario -> model_api calls (mirrors ui/sidebar.build + ui/common.run_model)
  memo.py      thread-safe LRU memo with in-flight de-duplication
frontend/      React 19 + TypeScript + Vite
  src/views/Simple.tsx        general-audience (Korean) view  = ui/simple.py
  src/views/Expert.tsx        research (English) view         = ui/expert.py
  src/views/ExpertSidebar.tsx Expert sidebar                   = ui/sidebar.py
  src/plantSvg.ts             plant map                         = src/plot_svg.py (byte-identical)
  src/format.ts               Python-compatible "{:.3g}" / "{:.0f}" formatting
src/report_html.py            one-page summary HTML, shared by both front ends
```

The bilingual copy shared by both front ends (disclaimer, glossary, inverse and
custom-table strings, congener labels, presets, scenarios) is served from
`ui/i18n.py` through `/api/meta`, so the two UIs cannot drift apart. Every default the
web UI starts from (the transport variant, the crop, the intake and body weight, the
site f_oc, the inventory seeds) is read from `model_api` / `soil_hydrus` through
`/api/meta` or the seed endpoints, never typed in the front end.

## The scenario (what the sidebar sends)

```jsonc
{
  "compound":  {"kind": "congener" | "smiles" | "neutral", "congener": "PFOA", "smiles": null,
                "neutral": {"log_kow": 2.45, "half_life": 7, "tscf_model": "briggs", "lipid_source": "measured",
                            "phloem": false, "air": false, "weak_electrolyte": false, "g_apo": 0,
                            "Koc": null, "grain_feed": "xylem", ...}},
  "mechanism": {"crop": "site2020" | "oryza2000" | "oryza" | "growth_rice",
                "shoot": "redistributed" | "4pool", "fxy_source": "fit" | "recommended",
                "lipid_loading": false},
  "exposure":  {"mode": "parametric" | "tables" | "csv" | "hydrus" | "inventory" | "biomon",
                "Cwo": 1.0, "season": 120, "cwo_profile": "constant" | "flooded",
                "tables": {...}, "csv_text": null, "f_oc": null,
                "C_total": null, "K_F": null, "n_F": 1.0, "biomon": null, "biomon_csv_text": null},
  "lang": "en"
}
```

* `crop` is the Streamlit "Crop drivers" choice (`ui.common.CROPS` / `crop_kw`): the
  calibration site-year (default; it sets its own season, transplanting to harvest) or a
  generic-climatology biomass driver, for which `season` (90–160 d) applies.
* `shoot` / `fxy_source` / `lipid_loading` pick the fitted variant (`model_api.variant_key`).
  As in the sidebar, f_xy `fit` with lipid loading or on the 4-pool shoot runs
  `recommended` with a warning message; a withdrawn variant (`WITHDRAWN_VARIANTS`,
  PREREG 9.4) is flagged in the payload (`variant.withdrawn`) and in the Expert headline.
* `cwo_profile: "flooded"` is the stored HYDRUS shape of the Tsukuba Andosol paddy
  (`cwo_profile="flooded"`, `params/cwo_shapes_hydrus.csv`); there is no leaching-rate
  setting.
* The soil inventory seeds (`C_total`, `K_F` left `null`) are the Streamlit ones: the
  Yamazaki pre-season soil (SI Table S16) and the site Kd (Koc × f_oc 3.75 %);
  `pore_water_from_inventory` puts the flooded plough layer's water and the paddy's
  HYDRUS shape on them.
* The live HYDRUS run takes only `f_oc` (`null` = the site's 0.0375); the water calendar
  and percolation are sourced (PREREG 9.2–9.14).
* `grain_feed` is `simulate_neutral`'s own switch for the neutral brown rice
  (`"xylem"` = the model default, PREREG 9.15; `null` = not predicted).

Unspecified fields take the model defaults; unknown or removed fields are rejected (422).

## API

| Endpoint | What |
|---|---|
| `GET /api/meta` | congeners, labels, presets, scenarios, i18n strings, crops, the default variant and withdrawn variants, intake defaults, the a-priori band, neutral scores, what is installed (RDKit / HYDRUS / CompTox key) |
| `POST /api/run` | a `Scenario` → the run payload: time series, finals (root, straw, brown rice, white rice), fresh and dry BAFs, bands, input-range factors, EFSA intake on white rice, variant, plant-map colour limits, messages |
| `POST /api/figure` | `{scenario, name, lang, opts}` → `{figure: <Plotly JSON>, data}`; names: `tissue burden drivers soil_profile isotherm buildup_plain where_plain baf baf_observable biomon_compare chain compare congener_compare tang_tf intake_gauge` |
| `POST /api/cwo-preview` | the flooded Cwᵒ(t) preview under the shape radio (`plots.fig_cwo_profile`) |
| `POST /api/inverse` | Bayesian exposure estimate, streamed as NDJSON (`progress` lines, then `result`); `rsd_percent` is the user's analytical RSD (0 = not given), `view` picks the Simple or Expert mechanism set |
| `POST /api/download/{kind}` | `summary_csv timeseries_csv summary_md summary_html map_html` |
| `GET /api/tables/defaults`, `/api/inventory/defaults`, `/api/smiles/svg`, `/api/ctx` | seeds for the editable tables (`crop`, `season`, `Cwo`) and for the soil inventory, the RDKit structure drawing, the CompTox lookup |
| `GET /api/hydrus/status`, `POST /api/hydrus/build` | live HYDRUS-1D engine |

`baf` is this run's dry-basis BAF against the Yamazaki bars (the Simple "얼마나 쌓이나"
expander, and the Expert tab for a SMILES compound); `baf_observable` is the model's
version of the Yamazaki observation (`yamazaki_harvest_baf`: the site-year, two harvests,
the congener's irrigation level), following the Expert mechanism switches
(`opts.view = "expert"`) or the f_xy source only (`"simple"`), as the two Streamlit tabs do.

Example: `curl -s localhost:8000/api/run -H 'content-type: application/json' -d '{"compound":{"congener":"PFOS"},"exposure":{"Cwo":10}}'`.

## Differences from the Streamlit app

* **Driver dicts are passed whole.** The Streamlit cache tuple (`ui/common._drivers_tuple`)
  carries every key the PFAS model reads; the one key it drops is `area` (the leaf's
  air-exchange surface), which matters only for a neutral compound with plant-air
  exchange on. The web run keeps it.
* **Two Streamlit inputs that the current model refuses are not offered**: the
  "Crop drivers on (off = illustrative placeholder logistics)" checkbox
  (`measured_forcing=False` raises: the placeholders were removed) and the live-HYDRUS
  "Flooded until" / "Percolation excess" sliders (`soil_hydrus.run_paddy_hydrus` has no
  such keywords; the Streamlit live mode fails with a `TypeError`). The live run's f_oc
  starts at the site's 0.0375 (`soil_hydrus.FOC`), not the Streamlit slider's 0.02.
* The neutral **brown-rice feed** (`grain_feed`) is a checkbox (default on, the model
  default); the Streamlit sidebar runs the default without offering it.
* Expert **Custom tables** mode shows the fallback constant Cwᵒ it uses where the Cwᵒ
  table is empty (Streamlit uses a hidden 1 µg/L).
* The plant-map **PNG** is rendered in the browser from the SVG (no kaleido/Chrome
  needed on the server); the Expert downloads also offer the one-page summary.

## What the port changed (from the pre-amendment-1 version, PR #71)

Removed from the UI and the API, because the model no longer has them: the membrane
potential `E_m`, the carrier / apoplastic-bypass root entry and its `Vmax ×`, `Km ×`,
`g_apo` overrides (the linear equilibrium entry `r_eq` replaced them), `k_leach` and the
inventory's own flooding/`θ_g` (the flooded shape is the stored HYDRUS shape and the
inventory uses the site's plough-layer water), the per-organ tissue density and implied
volume (per-volume reporting removed, PREREG 9.6), the "W2 fit" f_xy, the measured-forcing
switch, the two-pool overlays and the Tang refit bars, the eᴺ metric, the precision
presets of the inverse.

Added: the crop-driver choice (site-year default), the shoot model, the f_xy `fit` /
`recommended` choice with the sidebar's fall-backs and the withdrawn flag, dry-basis
BAFs in the Expert headline, the Yamazaki two-harvest observable, the white-rice intake
(`api.white_rice_conc`, brown rice as the fallback), the brown-rice labels, the
input-range band caption, the growth-table interpolation spread, the soil-inventory seeds,
the RSD input of the inverse, the current Simple trust table and copy.

## Tests

* `tests/test_web_api.py` (pytest, FastAPI `TestClient`): the default run is
  **bit-identical** to `model_api.simulate()`, the variant switches and their fall-backs
  reach the model, removed parameters are refused, every data source and compound class
  runs, the inventory seeds and the flooded preview are the sidebar's, figures and
  downloads build, the Yamazaki observable follows the view, the inverse streams and
  recovers a known exposure, and the Streamlit `_summary_html` goes through the shared
  `report_html` builder.
* `frontend/src/__tests__/` (vitest): the TypeScript plant map is byte-identical to
  `plot_svg.plant_svg` on Python-generated fixtures, and the number formatting
  matches Python's `format()` including half-to-even ties.
* CI: `.github/workflows/web.yml` type-checks, tests and builds the front end, and runs
  `tests/test_web_api.py` with `requirements-web.txt` (in `tests.yml` that module skips:
  it installs `requirements.txt` only).
