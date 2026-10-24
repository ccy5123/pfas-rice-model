# Web app — FastAPI backend + React front end

A second front end for the same model, built because the Streamlit dashboard
(`app.py`) feels slow. It does **not** change the model: every number still comes
from `src/model_api.py`, and every chart from `src/plots.py`.

## Why the Streamlit app feels slow (and what this changes)

Streamlit re-runs the **whole script on the server on every widget touch**, and it
executes the body of **every tab**, visible or not. In this app that means each
interaction re-checks the sidebar, rebuilds every tab's figures, and — in Expert
mode — re-solves the uncached Tang validation (2–3 solves, ~0.6 s each). A single
ODE solve itself costs ~0.6–1.8 s (`scipy` BDF), so most of the waiting is the
same computation being asked for again.

The web app splits the work differently:

| | Streamlit (`app.py`) | Web app (`server/` + `frontend/`) |
|---|---|---|
| Day slider on the plant map | full script re-run | **in the browser** — the season is shipped once, the SVG is re-drawn client-side (`frontend/src/plantSvg.ts`, byte-identical to `src/plot_svg.py`) |
| Switching tabs / colour metric | full script re-run (all tabs computed) | in the browser; a tab fetches **its own** figure only when opened |
| Changing a model input | full re-run | one `/api/run` call after a 300 ms debounce; the previous result stays on screen until the new one arrives |
| Repeat requests | `st.cache_data` per function | a shared, de-duplicating LRU memo per solve / run / figure (`server/memo.py`) — identical concurrent requests wait for one solve |
| Tang check, two-pool overlays, inverse | Tang uncached; others cached | all memoised; the inverse streams live progress |
| First visit | cold solves | the Simple view's default + preset scenarios are pre-solved at server start (background warm-up) |

What does **not** get faster: a genuinely new solve still takes ~0.6–1.8 s, because
that is the scipy ODE itself. Making *that* faster is a separate job (a compiled
RHS / solver), not a front-end one.

## Measured (same interactions, both front ends)

Headless Chromium, same machine, Expert mode on PFOA unless noted; time from the
action until the screen has settled (Streamlit: its "running" indicator gone; web:
no pending result/figure). The web timings include a fixed ~350 ms settle wait in
the harness, so they are upper bounds; the server had been running (warm memo).

| Interaction | Streamlit | Web |
|---|---:|---:|
| First load (Simple view) | 6.1 s | 0.47 s |
| Switch Simple → Expert | 5.3 s | 1.6 s |
| Move the day slider on the plant map | 2.0 s | 0.03 s (no server call) |
| Open the "Tissue dynamics" tab | 3.3 s | 0.56 s |
| Change congener PFOA → PFOS | 2.7 s | 0.43 s |
| Change back PFOS → PFOA (already computed) | 2.0 s | 0.36 s |

Read the last two rows carefully: the Streamlit app takes ~2 s even for a result
it has cached, because the whole script (every tab, including the uncached Tang
solves) re-runs; that overhead is what the web app removes. A congener the server
has **never** solved still costs one ODE solve (~0.6–1.8 s) in either app.

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
`ui/i18n.py` through `/api/meta`, so the two UIs cannot drift apart.

## API

| Endpoint | What |
|---|---|
| `GET /api/meta` | congeners, labels, presets, scenarios, i18n strings, defaults, what is installed (RDKit / HYDRUS / CompTox key) |
| `POST /api/run` | a `Scenario` → the run payload: time series, finals, BAFs, bands, EFSA intake, plant-map colour limits, messages |
| `POST /api/figure` | `{scenario, name, lang, opts}` → `{figure: <Plotly JSON>, data}`; names: `tissue burden drivers soil_profile isotherm buildup_plain where_plain baf biomon_compare chain compare congener_compare tang_tf intake_gauge` |
| `POST /api/cwo-preview` | the flooded Cwᵒ(t) preview next to the `k_leach` slider |
| `POST /api/inverse` | Bayesian exposure estimate, streamed as NDJSON (`progress` lines, then `result`) |
| `POST /api/download/{kind}` | `summary_csv timeseries_csv summary_md summary_html map_html` |
| `GET /api/tables/defaults`, `/api/k-leach`, `/api/smiles/svg`, `/api/ctx` | seeds for the editable tables, the calibrated `k_leach`, the RDKit structure drawing, the CompTox lookup |
| `GET /api/hydrus/status`, `POST /api/hydrus/build` | live HYDRUS-1D engine |

Example: `curl -s localhost:8000/api/run -H 'content-type: application/json' -d '{"compound":{"congener":"PFOS"},"exposure":{"Cwo":10}}'`.
Unspecified fields take the Streamlit defaults; unknown fields are rejected (422).

## Deliberate differences from the Streamlit app

* **Driver modes keep ORYZA's leaf senescence.** `drivers_from_arrays` attaches the
  ORYZA leaf death rate (`leaf_loss`) so a shrinking leaf does not spuriously
  concentrate PFAS (CLAUDE.md §6, "Leaf senescence-loss flux"). The Streamlit cache
  packs drivers into a `(t, Cwo, Qtp, M)` tuple and **drops it**, so in the
  Soil-inventory and live-HYDRUS modes the Streamlit app reports a leaf BAF ~2×
  too high (PFOA soil inventory: leaf 6.26 vs 3.13, straw 1.05 vs 0.53). The web
  backend passes the driver dict whole. `tests/test_web_api.py` pins this.
* The plant-map **PNG** is rendered in the browser from the SVG (no kaleido/Chrome
  needed on the server).
* The flooded-shape **preview for a neutral compound** uses that compound's own Koc
  (the Streamlit preview falls back to a generic PFCA shape).
* Expert **Custom tables** mode shows the fallback constant Cwᵒ it uses where the Cwᵒ
  table is empty (Streamlit uses a hidden 1 µg/L).
* The Bayesian inverse's precision presets are labelled in Korean in the Simple view.

## Tests

* `tests/test_web_api.py` (pytest, FastAPI `TestClient`): the default run is
  **bit-identical** to `model_api.simulate()`, every data source and compound class
  runs, figures/downloads build, the inverse streams and recovers a known exposure,
  and the leaf-loss difference above is real.
* `frontend/src/__tests__/` (vitest): the TypeScript plant map is byte-identical to
  `plot_svg.plant_svg` on Python-generated fixtures, and the number formatting
  matches Python's `format()` including half-to-even ties.
* CI: `.github/workflows/web.yml` type-checks, tests and builds the front end.
