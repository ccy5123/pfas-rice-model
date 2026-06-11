"""
PFAS–Rice Uptake Model — interactive dashboard
==============================================

Run locally:
    pip install -r requirements.txt -r requirements-app.txt
    streamlit run app.py

Pick a PFAS congener and scenario (pore-water concentration, root membrane
potential, f_xy source, measured vs placeholder forcings) and see the predicted
tissue concentrations/BAFs, binding factors, the chain-length parameter trends,
and the measured crop-physiology forcings. Compute lives in src/model_api.py.
"""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import streamlit as st

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))
import model_api as api          # noqa: E402
import forcing_rice as fr        # noqa: E402
import growth_rice as gr         # noqa: E402

st.set_page_config(page_title="PFAS–Rice Uptake Model", layout="wide")
st.title("PFAS–Rice Compartmental Uptake Model")
st.caption("Mechanistic 4-compartment dynamic model for permanently-anionic PFAS in paddy rice "
           "(IOC extension of the Trapp/Brunetti DPU). Parameters are measured/cited "
           "(see docs/literature_db). Outputs are illustrative, not regulatory.")

# ---------------------------------------------------------------- sidebar
with st.sidebar:
    st.header("Scenario")
    congener = st.selectbox("PFAS congener", api.CONGENERS,
                            index=api.CONGENERS.index("PFOA"))
    Cwo = st.number_input("Pore-water C_wᵒ  [µg/L]", min_value=0.0, value=1.0, step=0.1,
                          help="Free anion conc. driving root uptake. 1.0 → tissue conc equals BAF.")
    E_m = st.slider("Root membrane potential E_m  [mV]", -160, -90, -120, 5,
                    help="GHK anion-exclusion lever (rice −116…−140 mV; NH₄⁺ depolarises).")
    fxy_label = st.radio("Root→shoot loading f_xy",
                         ["recommended (monotone, physical)", "W2 fit (reproduces Yamazaki)"])
    fxy_source = "recommended" if fxy_label.startswith("recommended") else "W2fit"
    measured = st.checkbox("Measured forcings (Q_TP, M_s)", value=True,
                           help="On: transpiration from Kumari/NayHtoon, biomass from ORYZA IR72. "
                                "Off: illustrative logistic placeholders.")
    season = st.slider("Season length  [days]", 90, 160, 120, 5)

res = api.simulate(congener, Cwo=Cwo, E_m_mV=E_m, f_xy_source=fxy_source,
                   season=float(season), measured_forcing=measured)
obs = api.observed_baf(congener)
p = res["params"]

# ---------------------------------------------------------------- headline metrics
st.subheader(f"{congener}  (C{p['n_C']} {p['group']})")
c1, c2, c3, c4 = st.columns(4)
for col, tis in zip((c1, c2, c3), ("root", "straw", "grain")):
    pred = res["straw_baf"] if tis == "straw" else res["baf_final"][tis]
    delta = f"obs {obs[tis]:.2f}" if tis in obs else "no obs"
    col.metric(f"{tis} BAF [L/kg]", f"{pred:.2f}", delta, delta_color="off")
c4.metric("anion exclusion eᴺ", f"{res['eN']:.0f}", f"N={res['N']:.2f}", delta_color="off")

tab1, tab2, tab3, tab4 = st.tabs(
    ["Tissue dynamics", "BAF vs observed", "Chain-length trends", "Forcings"])

# ---- tab 1: concentration-time
with tab1:
    fig, ax = plt.subplots(figsize=(7, 4))
    for tis in api.TISSUES:
        ax.plot(res["t"], res["conc"][tis], lw=2, label=tis)
    ax.plot(res["t"], res["straw"], "k--", lw=1, label="straw (stem+leaf)")
    ax.set_xlabel("days after transplant"); ax.set_ylabel("tissue conc [µg/kg]")
    ax.set_title(f"{congener} tissue concentrations"); ax.legend()
    st.pyplot(fig)
    st.markdown(f"**Binding factor B_k [L/kg fw]** — "
                + ", ".join(f"{k}: {v:.2f}" for k, v in res["B_k"].items())
                + f"  ·  f_xy={p['f_xy']:.4g}, L_Ph={p['L_Ph']:.4g}, κ_d={p['kappa_d']:.3g}")

# ---- tab 2: BAF vs obs
with tab2:
    tissues = ["root", "straw", "grain"]
    pred = [res["baf_final"]["root"], res["straw_baf"], res["baf_final"]["grain"]]
    fig, ax = plt.subplots(figsize=(6, 4))
    x = np.arange(len(tissues)); w = 0.38
    ax.bar(x - w / 2, pred, w, label="model")
    if obs:
        ax.bar(x + w / 2, [obs.get(t_, np.nan) for t_ in tissues], w, label="Yamazaki 2023")
    ax.set_xticks(x); ax.set_xticklabels(tissues); ax.set_ylabel("BAF [L/kg]")
    ax.set_title(f"{congener}: predicted vs observed BAF"); ax.legend()
    st.pyplot(fig)
    if not obs:
        st.info("No Yamazaki BAF for this congener (model prediction only).")

# ---- tab 3: chain-length trends
with tab3:
    rows = api.chain_table()
    nC = {g: [r["n_C"] for r in rows if r["group"] == g] for g in ("PFCA", "PFSA")}
    def series(key, g): return [r[key] for r in rows if r["group"] == g]
    fig, axs = plt.subplots(2, 2, figsize=(10, 7))
    panels = [("K_PL", "K_PL [L/kg] (log)", True), ("f_xy_recommended", "f_xy recommended (log)", True),
              ("B_root", "B_root [L/kg]", False), ("K_cw_root", "K_cw root [L/kg]", False)]
    for ax, (key, lab, logy) in zip(axs.ravel(), panels):
        for g, mk in (("PFCA", "o-"), ("PFSA", "s--")):
            ax.plot(nC[g], series(key, g), mk, label=g)
        sel = next(r for r in rows if r["name"] == congener)
        ax.scatter([sel["n_C"]], [sel[key]], s=120, facecolors="none", edgecolors="red", zorder=5)
        ax.set_xlabel("perfluoro-C"); ax.set_ylabel(lab)
        if logy: ax.set_yscale("log")
        ax.legend(fontsize=8)
    fig.suptitle(f"Chain-length parameter trends (red ring = {congener})")
    fig.tight_layout()
    st.pyplot(fig)

# ---- tab 4: forcings (the measured crop-physiology drivers)
with tab4:
    t = res["t"]
    b = gr.organ_biomass(t, float(season))
    fig, (axa, axb) = plt.subplots(1, 2, figsize=(11, 4))
    axa.plot(t, fr.Q_TP(t, float(season)), lw=2, color="tab:blue")
    axa.set_xlabel("days"); axa.set_ylabel("Q_TP [L/day/hill]")
    axa.set_title("Transpiration stream (Kumari 2022 + Nay Htoon 2018)")
    for k in ("root", "stem", "leaf", "grain"):
        axb.plot(t, b[k], lw=2, label=k)
    axb.set_xlabel("days"); axb.set_ylabel("biomass [kg/hill]")
    axb.set_title("Organ biomass M_s(t) (ORYZA IR72)"); axb.legend()
    st.pyplot(fig)
    st.caption("These measured forcings drive the run when 'Measured forcings' is on "
               "(otherwise illustrative logistic placeholders are used).")

st.caption("Open item (recorded, not modelled): long-chain PFCA shoot accumulation "
           "(hysteretic sorption) — see docs/nstem_gradient_exploration.md.")
