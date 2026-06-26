"""
PFAS–Rice Uptake Model — interactive dashboard (Plotly + Streamlit)
===================================================================

Run locally:
    pip install -r requirements.txt
    streamlit run app.py

What this tool does
-------------------
Estimates how much PFAS ("forever chemicals") a rice plant takes up from
contaminated paddy water/soil, and where it ends up (roots, straw, grain). It
draws the soil + 4-compartment rice plant to scale and colours each compartment
by its PFAS accumulation (a heat colormap you can scrub through the season).

Two audiences, one app
----------------------
* **Simple mode** (default) — plain language, one chemical + one contamination
  level → a clear picture of where the PFAS goes. No jargon, no expert sliders.
* **Expert / advanced** (sidebar toggle) — restores the full research UI: five
  exposure modes (parametric, HYDRUS/CSV, live HYDRUS-1D, soil inventory,
  biomonitoring), SMILES structure input, and every model parameter.

Compute lives in src/model_api.py; the Plotly builders in src/plots.py (both
UI-agnostic and unit-tested head-less).  See docs/visualization_tool.md.
"""
import streamlit as st

st.set_page_config(page_title="PFAS–Rice Uptake Model", layout="wide")

from ui import common, sidebar, simple, expert     # noqa: E402

cfg = sidebar.build()
common.render_header(cfg)
common.render_custom_tables_panel(cfg)
common.run_model(cfg)                               # may st.stop() on a bad SMILES
if cfg.expert:
    expert.render(cfg)
else:
    simple.render(cfg)
common.render_footer(cfg)
