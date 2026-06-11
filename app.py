"""
PFAS–Rice Uptake Model — interactive dashboard (Plotly)
=======================================================

Run locally:
    pip install -r requirements.txt -r requirements-app.txt
    streamlit run app.py

Pick a PFAS congener and scenario; the charts are interactive (hover for values,
drag to zoom, click legend entries to toggle). Compute is in src/model_api.py;
the Plotly figure builders are in src/plots.py (both UI-agnostic).
"""
import os
import sys

import streamlit as st

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))
import model_api as api          # noqa: E402
import plots                     # noqa: E402

st.set_page_config(page_title="PFAS–Rice Uptake Model", layout="wide")
st.title("PFAS–Rice Compartmental Uptake Model")
st.caption("Mechanistic 4-compartment dynamic model for permanently-anionic PFAS in paddy rice "
           "(IOC extension of the Trapp/Brunetti DPU). Parameters are measured/cited "
           "(docs/literature_db). Charts are interactive — hover, zoom, toggle. Outputs illustrative.")

# ---------------------------------------------------------------- sidebar
with st.sidebar:
    st.header("Scenario")
    congener = st.selectbox("PFAS congener", api.CONGENERS, index=api.CONGENERS.index("PFOA"))
    Cwo = st.number_input("Pore-water C_wᵒ  [µg/L]", min_value=0.0, value=1.0, step=0.1,
                          help="Free anion conc. driving root uptake. 1.0 → tissue conc equals BAF.")
    E_m = st.slider("Root membrane potential E_m  [mV]", -160, -90, -120, 5,
                    help="GHK anion-exclusion lever (rice −116…−140 mV; NH₄⁺ depolarises).")
    fxy_label = st.radio("Root→shoot loading f_xy",
                         ["recommended (monotone, physical)", "W2 fit (reproduces Yamazaki)"])
    fxy_source = "recommended" if fxy_label.startswith("recommended") else "W2fit"
    measured = st.checkbox("Measured forcings (Q_TP, M_s)", value=True,
                           help="On: transpiration from Kumari/NayHtoon, biomass from ORYZA IR72.")
    season = st.slider("Season length  [days]", 90, 160, 120, 5)
    st.divider()
    compare = st.multiselect("Compare congeners (overlay)", api.CONGENERS,
                             default=["PFBA", "PFOA", "PFDA", "PFOS"],
                             help="Shown in the 'Compare' tab.")

res = api.simulate(congener, Cwo=Cwo, E_m_mV=E_m, f_xy_source=fxy_source,
                   season=float(season), measured_forcing=measured)
obs = api.observed_baf(congener)
p = res["params"]

# ---------------------------------------------------------------- headline metrics
st.subheader(f"{congener}  (C{p['n_C']} {p['group']})")
c1, c2, c3, c4 = st.columns(4)
for col, tis in zip((c1, c2, c3), ("root", "straw", "grain")):
    pred = res["straw_baf"] if tis == "straw" else res["baf_final"][tis]
    col.metric(f"{tis} BAF [L/kg]", f"{pred:.2f}",
               f"obs {obs[tis]:.2f}" if tis in obs else "no obs", delta_color="off")
c4.metric("anion exclusion eᴺ", f"{res['eN']:.0f}", f"N={res['N']:.2f}", delta_color="off")

tab1, tab2, tab3, tab4, tab5 = st.tabs(
    ["Tissue dynamics", "BAF vs observed", "Chain-length trends", "Compare congeners", "Forcings"])

with tab1:
    st.plotly_chart(plots.fig_tissue(res), use_container_width=True)
    st.markdown(f"**B_k [L/kg fw]** — " + ", ".join(f"{k}: {v:.2f}" for k, v in res["B_k"].items())
                + f"  ·  f_xy={p['f_xy']:.4g}, L_Ph={p['L_Ph']:.4g}, κ_d={p['kappa_d']:.3g}")

with tab2:
    st.plotly_chart(plots.fig_baf(res, obs), use_container_width=True)
    if not obs:
        st.info("No Yamazaki BAF for this congener (model prediction only).")

with tab3:
    key = st.selectbox("Parameter", ["K_PL", "K_prot", "K_cw_root",
                                     "f_xy_recommended", "B_root", "B_grain"], index=0)
    st.plotly_chart(plots.fig_chain(api.chain_table(), congener, key), use_container_width=True)

with tab4:
    tissue = st.radio("Tissue", ["root", "straw", "grain"], index=1, horizontal=True)
    if compare:
        results = {nm: api.simulate(nm, Cwo=Cwo, E_m_mV=E_m, f_xy_source=fxy_source,
                                    season=float(season), measured_forcing=measured)
                   for nm in compare}
        st.plotly_chart(plots.fig_compare(results, tissue), use_container_width=True)
    else:
        st.info("Select congeners in the sidebar to compare.")

with tab5:
    st.plotly_chart(plots.fig_forcings(res["t"], float(season)), use_container_width=True)
    st.caption("These measured forcings drive the run when 'Measured forcings' is on.")

st.caption("Open item (recorded, not modelled): long-chain PFCA shoot accumulation "
           "(hysteretic sorption) — see docs/nstem_gradient_exploration.md.")
