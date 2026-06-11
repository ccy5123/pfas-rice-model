"""
Interactive (Plotly) figure builders for the PFAS-rice dashboard.

Pure functions: take results from `model_api` and return Plotly Figures. No
Streamlit import, so they can be unit-tested head-less.  `app.py` renders them
with `st.plotly_chart`.
"""
from __future__ import annotations
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

import model_api as api
import forcing_rice as fr
import growth_rice as gr

_COL = {"root": "#8c564b", "stem": "#2ca02c", "leaf": "#1f77b4", "grain": "#ff7f0e"}
_LAYOUT = dict(template="plotly_white", hoverlabel=dict(namelength=-1),
               margin=dict(l=60, r=20, t=50, b=50))


def fig_tissue(res):
    """Tissue concentration vs time (hover-unified, legend-toggle, zoom)."""
    fig = go.Figure()
    for tis in api.TISSUES:
        fig.add_scatter(x=res["t"], y=res["conc"][tis], name=tis, mode="lines",
                        line=dict(width=2.5, color=_COL[tis]),
                        hovertemplate=f"{tis}: %{{y:.3g}} µg/kg<extra></extra>")
    fig.add_scatter(x=res["t"], y=res["straw"], name="straw", mode="lines",
                    line=dict(width=1.5, dash="dash", color="black"),
                    hovertemplate="straw: %{y:.3g} µg/kg<extra></extra>")
    fig.update_layout(title=f"{res['congener']} — tissue concentrations",
                      xaxis_title="days after transplant", yaxis_title="conc [µg/kg]",
                      hovermode="x unified", **_LAYOUT)
    return fig


def fig_baf(res, obs):
    """Predicted vs observed (Yamazaki) root/straw/grain BAF, grouped bars."""
    tis = ["root", "straw", "grain"]
    pred = [res["baf_final"]["root"], res["straw_baf"], res["baf_final"]["grain"]]
    fig = go.Figure()
    fig.add_bar(x=tis, y=pred, name="model", marker_color="#1f77b4",
                hovertemplate="model %{x}: %{y:.3g}<extra></extra>")
    if obs:
        fig.add_bar(x=tis, y=[obs.get(t_, None) for t_ in tis], name="Yamazaki 2023",
                    marker_color="#ff7f0e", hovertemplate="obs %{x}: %{y:.3g}<extra></extra>")
    fig.update_layout(barmode="group", title=f"{res['congener']} — predicted vs observed BAF",
                      yaxis_title="BAF [L/kg]", **_LAYOUT)
    return fig


_CHAIN_LOG = {"K_PL", "K_prot", "f_xy_recommended", "f_xy_W2fit", "B_root", "B_grain"}


def fig_chain(rows, congener, key="K_PL"):
    """A chosen per-congener parameter vs chain length; selected congener ringed."""
    fig = go.Figure()
    for g, dash in (("PFCA", "solid"), ("PFSA", "dash")):
        gr_ = [r for r in rows if r["group"] == g]
        fig.add_scatter(x=[r["n_C"] for r in gr_], y=[r[key] for r in gr_],
                        name=g, mode="lines+markers", line=dict(dash=dash),
                        text=[r["name"] for r in gr_],
                        hovertemplate="%{text}<br>nC=%{x}<br>" + key + "=%{y:.4g}<extra></extra>")
    sel = next(r for r in rows if r["name"] == congener)
    fig.add_scatter(x=[sel["n_C"]], y=[sel[key]], mode="markers", showlegend=False,
                    marker=dict(size=16, symbol="circle-open", color="red", line=dict(width=3)),
                    hovertemplate=f"{congener}<extra></extra>")
    fig.update_layout(title=f"{key} vs chain length", xaxis_title="perfluoro-C", yaxis_title=key,
                      yaxis_type="log" if key in _CHAIN_LOG else "linear",
                      hovermode="closest", **_LAYOUT)
    return fig


def fig_forcings(t, season):
    """Measured transpiration Q_TP(t) and ORYZA organ biomass (dual y-axis)."""
    fig = make_subplots(specs=[[{"secondary_y": True}]])
    fig.add_scatter(x=t, y=fr.Q_TP(t, season), name="Q_TP", line=dict(color="royalblue", width=2.5),
                    hovertemplate="Q_TP: %{y:.3f} L/d/hill<extra></extra>", secondary_y=False)
    b = gr.organ_biomass(t, season)
    for k in api.TISSUES:
        fig.add_scatter(x=t, y=b[k], name=k, line=dict(color=_COL[k]),
                        hovertemplate=f"{k}: %{{y:.4f}} kg/hill<extra></extra>", secondary_y=True)
    fig.update_layout(title="Measured forcings (Q_TP, M_s)", xaxis_title="days",
                      hovermode="x unified", **_LAYOUT)
    fig.update_yaxes(title_text="Q_TP [L/day/hill]", secondary_y=False)
    fig.update_yaxes(title_text="biomass [kg/hill]", secondary_y=True)
    return fig


def fig_compare(results, tissue="straw"):
    """Bar of a chosen tissue's BAF across several congeners (comparison view)."""
    names = list(results)
    def val(res):
        return res["straw_baf"] if tissue == "straw" else res["baf_final"][tissue]
    grp = [results[n]["params"]["group"] for n in names]
    fig = go.Figure(go.Bar(x=names, y=[val(results[n]) for n in names],
                           marker_color=["#1f77b4" if g == "PFCA" else "#ff7f0e" for g in grp],
                           hovertemplate="%{x}: %{y:.3g} L/kg<extra></extra>"))
    fig.update_layout(title=f"{tissue} BAF across selected congeners (blue=PFCA, orange=PFSA)",
                      yaxis_title="BAF [L/kg]", yaxis_type="log", **_LAYOUT)
    return fig
