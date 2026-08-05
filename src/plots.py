"""
Interactive (Plotly) figure builders for the PFAS-rice dashboard.

Pure functions: take results from `model_api` and return Plotly Figures. No
Streamlit import, so they can be unit-tested head-less.  `app.py` renders them
with `st.plotly_chart`.
"""
from __future__ import annotations
import re
import numpy as np
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots
from plotly.colors import sample_colorscale

import model_api as api
import forcing_rice as fr
import growth_rice as gr

_COL = {"root": "#8c564b", "stem": "#2ca02c", "leaf": "#1f77b4", "grain": "#ff7f0e"}
# plain-language tissue names for the general-audience (Simple) views
_PLAIN = {"root": "Roots", "stem": "Stems", "leaf": "Leaves", "grain": "Grain", "straw": "Straw"}
_PLAIN_KO = {"root": "뿌리", "stem": "줄기", "leaf": "잎", "grain": "낟알", "straw": "짚"}


def _plain(lang):
    return _PLAIN_KO if lang == "ko" else _PLAIN
# ---- brand Plotly template (applied to every non-schematic chart via _LAYOUT) ----
_CATEGORICAL = ["#0E7A63", "#E69F00", "#56B4E9", "#8E5AA8", "#D55E00", "#7A5230", "#0072B2", "#C7A93B"]
_SEQUENTIAL = [[0.0, "#F1F7F2"], [0.5, "#63BC9C"], [1.0, "#0A5A49"]]


def pfas_template():
    """Brand template — colours/typography only, THEME-NEUTRAL surface. Transparent
    paper/plot backgrounds and no hard-coded font colour, so when a chart is drawn
    with Streamlit's own theme (the default, NOT theme=None) it follows light/dark:
    the app paints the surface + text colour, while our brand colorway + neutral
    grid + margins carry through. (config.toml chartCategoricalColors match, so the
    palette is identical either way.)"""
    grid = "rgba(140,132,116,0.22)"
    return go.layout.Template(layout=dict(
        font=dict(family="Pretendard, 'Malgun Gothic', sans-serif", size=14),
        paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)", colorway=_CATEGORICAL,
        colorscale=dict(sequential=_SEQUENTIAL), title=dict(font=dict(size=17)),
        xaxis=dict(gridcolor=grid, zerolinecolor=grid, linecolor=grid),
        yaxis=dict(gridcolor=grid, zerolinecolor=grid, linecolor=grid),
        margin=dict(l=56, r=24, t=48, b=48)))


pio.templates["pfas"] = pfas_template()
_LAYOUT = dict(template="pfas", hoverlabel=dict(namelength=-1),
               margin=dict(l=60, r=20, t=50, b=50))
# Default "accumulation heat" colour scale for the plant map (more = hotter). A warm
# agricultural ramp cream→wheat→gold→ochre→terracotta: still reads "more = more
# intense", but the top is a soft terracotta rather than an alarm pure-red, so a
# high *relative* organ isn't misread as a danger flag (the colorbar carries the
# real, usually-small µg/kg number). Test suite checks only that figures build.
_HEAT = [[0.0, "#F4F7EE"], [0.22, "#EADF9B"], [0.5, "#E7B24C"],
         [0.75, "#D98A44"], [1.0, "#B4562E"]]


def _formed(res, j, y):
    """Mask the GRAIN's pre-formation transient. The panicle/grain is physically absent
    until ~flowering, but the ODE floors its mass to avoid 0/0, so a tiny xylem/phloem
    influx divides into a SPURIOUS concentration spike before the grain exists (it then
    crashes via growth dilution once the grain bulks up). Blank the grain while its mass is
    still <2% of its season max. Root/stem/leaf form early (no such floor), and a
    constant-mass driver (HYDRUS/CSV) never trips the threshold, so both are untouched."""
    y = np.asarray(y, float).copy()
    if api.TISSUES[j] == "grain":
        M = np.asarray(res["M"], float)[:, j]
        if M.max() > 0.0:
            y[M < 0.02 * M.max()] = np.nan
    return y


def fig_tissue(res):
    """Tissue concentration vs time (hover-unified, legend-toggle, zoom)."""
    fig = go.Figure()
    for j, tis in enumerate(api.TISSUES):
        fig.add_scatter(x=res["t"], y=_formed(res, j, res["conc"][tis]), name=tis, mode="lines",
                        line=dict(width=2.5, color=_COL[tis]),
                        hovertemplate=f"{tis}: %{{y:.3g}} µg/kg<extra></extra>")
    fig.add_scatter(x=res["t"], y=res["straw"], name="straw", mode="lines",
                    line=dict(width=1.5, dash="dash", color="black"),
                    hovertemplate="straw: %{y:.3g} µg/kg<extra></extra>")
    fig.update_layout(title=f"{res['congener']} — tissue concentrations",
                      xaxis_title="days after transplant", yaxis_title="conc [µg/kg]",
                      hovermode="x unified", **_LAYOUT)
    return fig


def fig_burden(res):
    """Per-tissue PFAS *mass* (burden) over the season [µg/hill] = C_k(t)·M_k(t).

    The chemical inventory in each organ (EXTENSIVE), complementing the intensive
    concentration in `fig_tissue`: a tissue can be high-concentration yet low-mass
    (small organ) or vice-versa. Organ biomass M_k(t) itself is in the Soil & drivers
    tab (`fig_drivers`)."""
    M = np.asarray(res["M"], float)
    fig = go.Figure()
    burden = {tis: np.asarray(res["conc"][tis], float) * M[:, j] for j, tis in enumerate(api.TISSUES)}
    for j, tis in enumerate(api.TISSUES):
        fig.add_scatter(x=res["t"], y=_formed(res, j, burden[tis]), name=tis, mode="lines",
                        line=dict(width=2.5, color=_COL[tis]),
                        hovertemplate=f"{tis}: %{{y:.3g}} µg<extra></extra>")
    fig.add_scatter(x=res["t"], y=sum(burden.values()), name="whole plant", mode="lines",
                    line=dict(width=1.5, dash="dash", color="black"),
                    hovertemplate="total: %{y:.3g} µg<extra></extra>")
    fig.update_layout(title=f"{res['congener']} — PFAS mass per tissue (burden)",
                      xaxis_title="days after transplant", yaxis_title="PFAS mass [µg/hill]",
                      **_LAYOUT)
    return fig


_BAF_EXTRA_COLORS = ["#2ca02c", "#9467bd", "#8c564b", "#17becf"]


def fig_baf(res, obs, extra=None, lang="en"):
    """Predicted vs observed (Yamazaki) root/straw/grain BAF, grouped bars.

    `extra` optionally overlays additional model series (e.g. the EXPLORATORY
    two-pool model): a dict {label: {"root","straw","grain"}} (values may be None);
    each is added as its own grouped bar so the canonical core, the overlay model(s)
    and the observed data are compared side by side. `lang="ko"` localises the axis,
    title and the core/observed legend names for the Simple (general-audience) view.
    """
    ko = lang == "ko"
    nm = _plain(lang)
    tis = ["root", "straw", "grain"]
    x = [nm[t_] for t_ in tis]
    pred = [res["baf_final"]["root"], res["straw_baf"], res["baf_final"]["grain"]]
    fig = go.Figure()
    fig.add_bar(x=x, y=pred, name="모델" if ko else "model (4-pool core)", marker_color="#1f77b4",
                hovertemplate="core %{x}: %{y:.3g}<extra></extra>")
    for i, (label, vals) in enumerate(dict(extra or {}).items()):
        fig.add_bar(x=x, y=[vals.get(t_) for t_ in tis], name=label,
                    marker_color=_BAF_EXTRA_COLORS[i % len(_BAF_EXTRA_COLORS)],
                    marker_line=dict(width=0.5, color="#333"),
                    hovertemplate=label + " %{x}: %{y:.3g}<extra></extra>")
    if obs:
        fig.add_bar(x=x, y=[obs.get(t_, None) for t_ in tis],
                    name="Yamazaki 2023 (실측)" if ko else "Yamazaki 2023",
                    marker_color="#ff7f0e", hovertemplate="obs %{x}: %{y:.3g}<extra></extra>")
    fig.update_layout(
        barmode="group",
        title=(f"{res['congener']} — 모델 예측 vs 실측 (Yamazaki)" if ko
               else f"{res['congener']} — predicted vs observed BAF"),
        yaxis_title="축적 배수 [L/kg]" if ko else "BAF [L/kg]", **_LAYOUT)
    return fig


def fig_buildup_plain(res, lang="en"):
    """Plain-language tissue concentration over the season (Simple-mode view).

    Same data as `fig_tissue` but with friendly tissue names (Roots/Stems/Leaves/
    Grain) and jargon-free axis/title text, for a general audience. `lang="ko"`
    renders Korean. Reuses the grain pre-formation mask so the empty pre-flowering
    period is not drawn."""
    nm = _plain(lang)
    ko = lang == "ko"
    fig = go.Figure()
    for j, tis in enumerate(api.TISSUES):
        fig.add_scatter(x=res["t"], y=_formed(res, j, res["conc"][tis]), name=nm[tis],
                        mode="lines", line=dict(width=2.8, color=_COL[tis]),
                        hovertemplate=f"{nm[tis]}: %{{y:.3g}} µg/kg<extra></extra>")
    fig.update_layout(
        title="한 철 동안 식물 속 PFAS 축적" if ko else "How PFAS builds up in the plant over the season",
        xaxis_title="모내기 후 일수" if ko else "days after the rice is transplanted",
        yaxis_title="조직 속 PFAS [µg/kg]" if ko else "PFAS in the plant tissue [µg per kg]",
        hovermode="x unified", **_LAYOUT)
    return fig


def fig_where_plain(res, lang="en", band=False):
    """Plain-language bar of the final PFAS level in roots / straw / grain.

    A jargon-free read of where the chemical ends up at harvest (the same numbers
    as the BAF bars, but labelled as a concentration build-up, no 'BAF' symbol).
    `lang="ko"` renders Korean. `band=True` overlays the coarse a-priori predictive
    uncertainty (`model_api.predictive_band`, a ×/÷ ~7 honesty band) as error bars,
    so the general-audience view never shows an absolute number without a spread."""
    nm = _plain(lang)
    ko = lang == "ko"
    order = ["root", "straw", "grain"]
    vals = {"root": res["conc"]["root"][-1], "straw": res["straw"][-1],
            "grain": res["conc"]["grain"][-1]}
    xs = [nm[t_] for t_ in order]
    ys = [float(vals[t_]) for t_ in order]
    cols = [_COL.get(t_, "#1f77b4") for t_ in order]
    if band:
        # The band is multiplicative (×/÷ ~7). On a LOG y-axis its up/down whiskers
        # are symmetric (log(v·7)−log(v) == log(v)−log(v/7)) instead of the lop-sided
        # look a linear axis gives, so we plot the estimate as a MARKER with error
        # bars (bars read wrong on a log scale — their baseline is arbitrary).
        bands = {t_: api.predictive_band(v) for t_, v in zip(order, ys)}
        err_y = dict(type="data", symmetric=False,
                     array=[bands[t_]["hi"] - v for t_, v in zip(order, ys)],
                     arrayminus=[v - bands[t_]["lo"] for t_, v in zip(order, ys)],
                     color="#9a9384", thickness=1.6, width=9)
        fig = go.Figure(go.Scatter(
            x=xs, y=ys, mode="markers+text", error_y=err_y,
            marker=dict(size=16, color=cols, line=dict(color="#ffffff", width=1.5)),
            text=[f"{v:.2g}" for v in ys], textposition="middle right",
            textfont=dict(size=13), cliponaxis=False,
            hovertemplate="%{x}: %{y:.3g} µg/kg<extra></extra>"))
    else:
        fig = go.Figure(go.Bar(
            x=xs, y=ys, marker_color=cols,
            text=[f"{v:.2g}" for v in ys], textposition="outside",
            hovertemplate="%{x}: %{y:.3g} µg/kg<extra></extra>"))
    fig.update_layout(
        title=(f"{res['congener']}가 식물에서 모이는 곳 (수확 시)" if ko
               else f"Where {res['congener']} ends up in the plant (at harvest)"),
        yaxis_title="조직 속 PFAS [µg/kg]" if ko else "PFAS in the tissue [µg per kg]",
        xaxis_title="벼 부위" if ko else "part of the rice plant", **_LAYOUT)
    if band:
        # decade ticks only (dtick=1 in log space) so the labels don't crowd/overlap
        fig.update_yaxes(type="log", dtick=1, tickformat=".3g", minor=dict(showgrid=False))
        fig.update_xaxes(range=[-0.5, 2.7])                 # room for the value labels
    return fig


def fig_exposure_posterior(est, lang="en"):
    """Posterior over the estimated soil-water contamination level Cwᵒ [µg/L].

    From `model_api.estimate_exposure_bayesian`: a shaded probability curve vs Cwᵒ
    (log x), with the 95% credible interval darker and the median (most likely
    value) marked. Plain labels for a general audience; `lang="ko"` renders Korean."""
    ko = lang == "ko"
    Cwo = np.asarray(est["grid"]["Cwo"], float)
    dens = np.asarray(est["grid"]["density"], float)
    med = est["median"]
    lo, hi = est["ci95"]
    fig = go.Figure()
    fig.add_scatter(x=Cwo, y=dens, mode="lines", line=dict(color="#1f77b4", width=2.5),
                    fill="tozeroy", fillcolor="rgba(31,119,180,0.12)", name="probability",
                    hovertemplate="%{x:.3g} µg/L<extra></extra>")
    if np.isfinite(lo) and np.isfinite(hi):
        m = (Cwo >= lo) & (Cwo <= hi)
        fig.add_scatter(x=Cwo[m], y=dens[m], mode="lines", line=dict(width=0),
                        fill="tozeroy", fillcolor="rgba(31,119,180,0.30)",
                        name="95% range", hoverinfo="skip")
    fig.add_vline(x=med, line=dict(color="#d62728", width=2, dash="dash"),
                  annotation_text=(f"가장 가능성 높음 {med:.3g} µg/L" if ko
                                   else f"most likely {med:.3g} µg/L"),
                  annotation_position="top")
    fig.update_layout(
        title=("토양수 오염 수준 추정 (불확실성 포함)" if ko
               else "Estimated contamination level in the soil water (with uncertainty)"),
        xaxis_title="토양수에 녹아 있는 PFAS [µg/L]" if ko else "PFAS dissolved in the soil water [µg/L]",
        yaxis_title="상대 확률" if ko else "relative probability", xaxis_type="log",
        showlegend=False, **_LAYOUT)
    fig.update_yaxes(rangemode="tozero", showticklabels=False)
    return fig


_SIGNAL_COLOR = {"green": "#2ca02c", "amber": "#e6a817", "red": "#d62728"}


def fig_intake_gauge(info, lang="en"):
    """Half-circle gauge: predicted grain PFAS as a % of the EFSA group TWI.

    `info` is a `model_api.intake_fraction` dict. A traffic-light gauge (green
    <10 %, amber 10-100 %, red >=100 %) with a threshold marker at 100 % of the
    tolerable weekly intake. Plain labels; `lang="ko"` renders Korean. This is an
    intake reference to a HEALTH-BASED guidance value, not a food-standard (MRL)
    check -- the surrounding UI copy carries that caveat."""
    ko = lang == "ko"
    pct = float(info["percent"])
    if not np.isfinite(pct):
        pct = 0.0
    axis_max = max(150.0, pct * 1.15)
    color = _SIGNAL_COLOR.get(info.get("signal"), "#888")
    title = ("주간 안전섭취량(EFSA) 대비 — 참고" if ko
             else "Share of the weekly safe intake (EFSA) — reference")
    fig = go.Figure(go.Indicator(
        mode="gauge+number", value=pct,
        number={"suffix": "%", "font": {"size": 42, "color": color}},
        gauge={
            "axis": {"range": [0, axis_max], "ticksuffix": "%", "tickwidth": 1},
            "bar": {"color": color, "thickness": 0.72},
            "borderwidth": 0,
            "steps": [
                {"range": [0, 10], "color": "rgba(44,160,44,0.16)"},
                {"range": [10, min(100, axis_max)], "color": "rgba(230,168,23,0.16)"},
                {"range": [min(100, axis_max), axis_max], "color": "rgba(214,39,40,0.16)"},
            ],
            "threshold": {"line": {"color": "#d62728", "width": 4},
                          "thickness": 0.85, "value": min(100.0, axis_max)},
        }))
    # transparent bg + neutral tick/title ink so the gauge blends into the page in
    # BOTH light and dark (it is rendered with theme=None to keep the signal colours)
    fig.update_layout(title=dict(text=title, x=0.5, xanchor="center", font=dict(color="#8a8a8a")),
                      margin=dict(l=30, r=30, t=60, b=10), height=280,
                      paper_bgcolor="rgba(0,0,0,0)", font=dict(color="#8a8a8a"))
    return fig


def fig_tang_tf(val, val_refit=None):
    """Tang 2026 per-organ TF (DRY weight): measured vs model, optional refit-f_xy bar.

    `val` (and optional `val_refit`) come from `model_api.tang_tf_validation`. Log-y
    grouped bars over stalk/leaf/endosperm.
    """
    organs = val["organs"]
    fig = go.Figure()
    fig.add_bar(x=organs, y=[val["tang_tf"].get(o) for o in organs], name="Tang 2026 (measured)",
                marker_color="#444444", hovertemplate="Tang %{x}: %{y:.2f}<extra></extra>")
    fig.add_bar(x=organs, y=[val["model_tf"][o] for o in organs],
                name=f"model · f_xy={val['f_xy']:.3g}", marker_color="#bdbdbd",
                marker_line=dict(width=0.5, color="#333"),
                hovertemplate="model %{x}: %{y:.3g}<extra></extra>")
    if val_refit is not None:
        fig.add_bar(x=organs, y=[val_refit["model_tf"][o] for o in organs],
                    name=f"model · Tang-refit f_xy={val_refit['f_xy']:.3g}", marker_color="#2ca02c",
                    marker_line=dict(width=0.5, color="#333"),
                    hovertemplate="refit %{x}: %{y:.3g}<extra></extra>")
    fig.update_layout(barmode="group", yaxis_type="log",
                      title=f"{val['congener']} — per-organ TF (dry weight) vs Tang 2026 (OOS)",
                      yaxis_title="TF = C_organ / C_root  [dry wt]",
                      xaxis_title="tissue (model grain ↔ Tang endosperm)", **_LAYOUT)
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


def fig_congener_compare(results, lang="en", order=None):
    """Grouped bars of the final root / straw / grain concentration across congeners.

    A policy 'at a glance' view of how the accumulation pattern shifts across PFAS
    chemicals: `results` = {congener: simulate() dict}. Ordered by `order` (a
    chain-length list) for the short- vs long-chain message. Friendly tissue names;
    `lang="ko"` renders Korean. Reads the same final concentrations as the where-it-
    ends-up bar, just across several chemicals side by side."""
    ko = lang == "ko"
    nm = _plain(lang)
    names = [c for c in order if c in results] if order else list(results)

    def _final(res, tis):
        return float(res["straw"][-1]) if tis == "straw" else float(res["conc"][tis][-1])
    fig = go.Figure()
    for tis in ("root", "straw", "grain"):
        col = _COL["stem"] if tis == "straw" else _COL[tis]
        fig.add_bar(name=nm[tis], x=names, y=[_final(results[c], tis) for c in names],
                    marker_color=col,
                    hovertemplate="%{x} " + nm[tis] + ": %{y:.3g} µg/kg<extra></extra>")
    # log y: the values span ~3 orders of magnitude (long-chain root ≫ short-chain
    # grain), so a linear axis lets one giant bar hide the grain pattern that is the
    # policy message. Log keeps every bar readable.
    fig.update_layout(
        barmode="group", yaxis_type="log",
        title=("물질에 따라 쌓이는 곳이 다릅니다 (수확 시)" if ko
               else "Where each chemical ends up (at harvest)"),
        yaxis_title="조직 속 PFAS [µg/kg, 로그]" if ko else "PFAS in the tissue [µg/kg, log]",
        xaxis_title="PFAS 물질 (사슬이 길어지는 순)" if ko else "PFAS chemical (increasing chain length)",
        **_LAYOUT)
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


# ===========================================================================
# Plant + soil MAP  --  the model drawn to scale, compartments coloured by the
# accumulation metric (a heat colormap).  This is the "see the model" view.
# ===========================================================================
# A rice plant (Oryza sativa): a fibrous root mass in a paddy-soil cross-section,
# an arching culm with broad leaves at the nodes, and a nodding golden ear of grain.
# The silhouette is the SVG illustration's geometry (coordinate space x:0..620,
# y:110..690 with y DOWN); we flip y into Plotly's y-up convention (yf = _K - y)
# so the paths render upright. Soil/rachis/nodes are fixed structural colours; the
# root/leaf/stem/grain organs are filled/stroked with the metric-scale colour at
# that compartment's value -- that is the accumulation map.
_K = 512.0                       # y-flip constant (yf = _K - y_svg), keeps yf > 0
_NONE = "rgba(0,0,0,0)"
_NUM = re.compile(r"[-+]?\d*\.?\d+")


def _flip_path(path):
    """Flip the y of every coordinate in an SVG path into y-up space (yf=_K-y).

    Every M/L/C/Q command consumes whole (x,y) pairs, so every 2nd number is a y."""
    n = 0
    def repl(m):
        nonlocal n
        n += 1
        v = float(m.group())
        return f"{(_K - v):g}" if n % 2 == 0 else f"{v:g}"
    return _NUM.sub(repl, path)


def _fy(y):
    return _K - y


def _ellipse_path(cx, cy, rx, ry, rot_deg=0.0, npts=16):
    """Closed path for a rotated ellipse (SVG coords) flipped into y-up space."""
    th = np.linspace(0.0, 2 * np.pi, npts, endpoint=False)
    a = np.radians(rot_deg)
    xs = cx + rx * np.cos(th) * np.cos(a) - ry * np.sin(th) * np.sin(a)
    ys = _K - (cy + rx * np.cos(th) * np.sin(a) + ry * np.sin(th) * np.cos(a))
    return " ".join(f"{'M' if i == 0 else 'L'} {x:.1f},{y:.1f}"
                    for i, (x, y) in enumerate(zip(xs, ys))) + " Z"


# --- rice-plant glyph (a clean bold icon, 512-box SOURCE coords, y-down) --------
# Per-organ colouring hooks are preserved: roots + seed → root colour, the four
# blades → leaf colour, the central culm → stem colour, the four top ellipses →
# grain colour. SVG elliptical arcs are SAMPLED to polylines so the y-flip
# (yf = _K - y, _K=512) works on plain coordinate pairs (as _flip_path needs).
import math as _math


def _svg_arc(x1, y1, rx, ry, phi_deg, large, sweep, x2, y2, n=26):
    """Sample an SVG 'A' elliptical-arc (endpoint parameterisation) to (x,y) pts."""
    phi = _math.radians(phi_deg)
    cp, sp = _math.cos(phi), _math.sin(phi)
    dx, dy = (x1 - x2) / 2.0, (y1 - y2) / 2.0
    x1p, y1p = cp * dx + sp * dy, -sp * dx + cp * dy
    rx, ry = abs(rx), abs(ry)
    lam = x1p * x1p / (rx * rx) + y1p * y1p / (ry * ry)
    if lam > 1:
        s = _math.sqrt(lam); rx *= s; ry *= s
    num = rx*rx*ry*ry - rx*rx*y1p*y1p - ry*ry*x1p*x1p
    den = rx*rx*y1p*y1p + ry*ry*x1p*x1p
    co = _math.sqrt(max(num, 0.0) / den) if den else 0.0
    if large == sweep:
        co = -co
    cxp, cyp = co * rx * y1p / ry, -co * ry * x1p / rx
    cx = cp * cxp - sp * cyp + (x1 + x2) / 2.0
    cy = sp * cxp + cp * cyp + (y1 + y2) / 2.0

    def _ang(ux, uy, vx, vy):
        d = (ux*vx + uy*vy) / (_math.hypot(ux, uy) * _math.hypot(vx, vy))
        a = _math.acos(max(-1.0, min(1.0, d)))
        return -a if (ux*vy - uy*vx) < 0 else a
    th1 = _ang(1, 0, (x1p - cxp)/rx, (y1p - cyp)/ry)
    dth = _ang((x1p - cxp)/rx, (y1p - cyp)/ry, (-x1p - cxp)/rx, (-y1p - cyp)/ry)
    if not sweep and dth > 0:
        dth -= 2*_math.pi
    if sweep and dth < 0:
        dth += 2*_math.pi
    pts = []
    for i in range(n + 1):
        t = th1 + dth * i / n
        pts.append((cp*rx*_math.cos(t) - sp*ry*_math.sin(t) + cx,
                    sp*rx*_math.cos(t) + cp*ry*_math.sin(t) + cy))
    return pts


def _pline(pts):
    """(x,y) SOURCE points -> a flipped plotly path 'M x,yf L x,yf …'."""
    return " ".join(f"{'M' if i == 0 else 'L'} {x:.2f},{_K - y:.2f}"
                    for i, (x, y) in enumerate(pts))


def _seg(x1, y1, x2, y2):
    return _pline([(x1, y1), (x2, y2)])


# roots (root colour) — two curls + a taproot; (path, px width)
_ROOTS = [(_pline(_svg_arc(246.0, 394.0, 42, 42, 0, 0, 1, 198.9, 435.7)), 13),
          (_pline(_svg_arc(266.0, 394.0, 42, 42, 0, 0, 0, 313.1, 435.7)), 13),
          (_seg(256, 395, 256, 452), 13)]
_SEED = (256, 378, 24)                                       # cx, cy(src), r → root
# four leaf blades (leaf colour)
_LEAVES = [_seg(178.1, 172.6, 216.9, 211.4),
           _seg(333.3, 215.5, 294.5, 254.3),
           _seg(178.3, 258.2, 217.1, 297.0),
           _seg(333.3, 296.7, 294.5, 335.5)]
# central culm (stem colour): vertical shaft + nodding top arc
_STEM_PATH = _pline([(256, 378), (256, 170)] + _svg_arc(256, 170, 80, 80, 0, 0, 1, 296, 100.7))
# four grain ellipses (grain colour) at the ear
_GRAINS = [(349.3, 63.5, 28.8, 18.7, -22), (398.7, 101.3, 28.8, 18.7, 27),
           (352.3, 123.1, 28.8, 18.7, 52), (409.9, 157.3, 28.8, 18.7, 51)]
_GRAIN_PATHS = [_ellipse_path(*g) for g in _GRAINS]
# soil band at the base (structural) with a gently wavy top; roots dip into it
_SOIL_TOPY = 424.0
_SOIL = _pline([(128, _SOIL_TOPY + 7), (200, _SOIL_TOPY - 3), (256, _SOIL_TOPY + 2),
                (320, _SOIL_TOPY - 3), (452, _SOIL_TOPY + 7), (452, 512), (128, 512)]) + " Z"
_SOIL_TOPLINE = _pline([(128, _SOIL_TOPY + 7), (200, _SOIL_TOPY - 3), (256, _SOIL_TOPY + 2),
                        (320, _SOIL_TOPY - 3), (452, _SOIL_TOPY + 7)])

# marker / label anchors (SOURCE coords -> flipped) shared by static + animated
_MARK = {"root": (256, _fy(415)), "stem": (272, _fy(300)), "leaf": (200, _fy(235)),
         "grain": (378, _fy(110)), "straw": (262, _fy(300))}
# label-box anchors (SOURCE coords -> flipped) + grow-outward (xanchor,yanchor) so
# the boxes stay inside the frame regardless of side.
_LABEL = {"root":  (150, _fy(452), "left",  "middle"),
          "stem":  (150, _fy(300), "left",  "middle"),
          "leaf":  (150, _fy(205), "left",  "middle"),
          "grain": (452, _fy(92),  "right", "top"),
          "straw": (150, _fy(285), "left",  "middle")}


def _frac(v, cmin, cmax):
    if v is None or not np.isfinite(v):
        return None
    if cmax <= cmin:
        return 0.5
    return float(np.clip((v - cmin) / (cmax - cmin), 0.0, 1.0))


def _color(v, cmin, cmax, scale, nan="#dcdcdc"):
    f = _frac(v, cmin, cmax)
    return nan if f is None else sample_colorscale(scale, [f])[0]


def _shapes_for(colors, line="#7a6a52"):
    """Return the rice-plant glyph + soil band shapes coloured by per-organ `colors`.

    The soil band is a fixed structural colour; the roots+seed, leaf blades, central
    culm and grain ellipses are stroked/filled with the metric-scale colour for that
    compartment's value (the accumulation map). Each coloured stroke gets a slightly
    wider dark under-stroke first, for a clean outlined-icon look."""
    S = [dict(type="path", path=_SOIL, fillcolor="#7A5A34", line=dict(width=0), layer="below"),
         dict(type="path", path=_SOIL_TOPLINE, fillcolor=_NONE,
              line=dict(color="#5A3E22", width=1.4), layer="below")]

    def stroke(pth, w, col):                                 # dark edge + colour fill
        S.append(dict(type="path", path=pth, fillcolor=_NONE, line=dict(color=line, width=w + 4)))
        S.append(dict(type="path", path=pth, fillcolor=_NONE, line=dict(color=col, width=w)))
    # central culm (behind), then leaf blades in front
    stroke(_STEM_PATH, 22, colors["stem"])
    for pth in _LEAVES:
        stroke(pth, 22, colors["leaf"])
    # roots (root colour) + seed/crown circle
    for pth, w in _ROOTS:
        stroke(pth, w, colors["root"])
    cx, cy, r = _SEED
    S.append(dict(type="circle", x0=cx - r, x1=cx + r, y0=_fy(cy) - r, y1=_fy(cy) + r,
                  fillcolor=colors["root"], line=dict(color=line, width=2)))
    # ear: nodding grain ellipses (grain colour)
    for pth in _GRAIN_PATHS:
        S.append(dict(type="path", path=pth, fillcolor=colors["grain"], line=dict(color="#B7861A", width=1.0)))
    return S


def _marker_points(values, straw_only):
    """[(name, value, x, y)] for the colourbar markers / labels (skips missing)."""
    order = ["root", "straw", "grain"] if straw_only else ["root", "stem", "leaf", "grain"]
    out = []
    for name in order:
        v = values.get(name, values.get("straw"))
        if v is None or not np.isfinite(v):
            continue
        out.append((name, float(v), _MARK[name][0], _MARK[name][1]))
    return out


def fig_plant_schematic(values, *, cmin, cmax, label="tissue conc [µg/kg]",
                        Cwo=None, colorscale=_HEAT, title=None, t=None,
                        obs=None, lang="en"):
    """Draw the rice plant + paddy soil with each compartment coloured by `values`.

    values : dict with root/stem/leaf/grain (and optional 'straw'). If stem/leaf are
        absent but 'straw' is present (e.g. biomonitoring root/straw/grain), the
        whole shoot is coloured by straw.
    cmin, cmax : shared colour limits (so colours are comparable across time/parts).
    Cwo : pore-water concentration to annotate in the soil [µg/L] (optional).
    obs : optional dict {tissue: value} drawn as a small reference label per organ.
    """
    straw = values.get("straw")
    root_v = values.get("root")
    stem_v = values.get("stem", straw)
    leaf_v = values.get("leaf", straw)
    grain_v = values.get("grain")
    straw_only = ("stem" not in values and "leaf" not in values and straw is not None)

    colors = {k: _color(v, cmin, cmax, colorscale)
              for k, v in (("root", root_v), ("stem", stem_v), ("leaf", leaf_v), ("grain", grain_v))}
    fig = go.Figure()
    fig.update_layout(shapes=_shapes_for(colors))

    # marker points (carry the colorbar + hover); labelled with arrows
    pts = _marker_points(values, straw_only)
    _ink = "#4a463c"                                     # dark ink (readable on the cream bg)
    fig.add_scatter(
        x=[p[2] for p in pts], y=[p[3] for p in pts], mode="markers",
        marker=dict(size=15, color=[p[1] for p in pts], colorscale=colorscale,
                    cmin=cmin, cmax=cmax, line=dict(color="#333", width=1),
                    colorbar=dict(title=dict(text=label, font=dict(color=_ink)),
                                  tickfont=dict(color=_ink), thickness=16, len=0.85,
                                  outlinecolor="#c9c2b2", outlinewidth=1)),
        text=[p[0] for p in pts], customdata=[p[1] for p in pts],
        hovertemplate="%{text}: %{customdata:.3g}<extra></extra>", showlegend=False)

    ko = lang == "ko"
    nmf = _PLAIN_KO if ko else None
    ann = []
    for name, v, mx, my in pts:
        lx, ly, xa, ya = _LABEL.get(name, (mx, my, "center", "middle"))
        disp = nmf.get(name, name) if nmf else name
        txt = f"<b>{disp}</b><br>{v:.3g}"
        if obs and name in obs and obs[name] is not None:
            _ob = ("관측" if ko else "obs")
            txt += f"<br><span style='color:#8a8270'>{_ob} {obs[name]:.3g}</span>"
        ann.append(dict(x=mx, y=my, ax=lx, ay=ly, axref="x", ayref="y",
                        xanchor=xa, yanchor=ya,
                        text=txt, showarrow=True, arrowhead=0, arrowcolor="#9a9384", arrowwidth=1,
                        font=dict(size=12, color=_ink), align="center",
                        bgcolor="rgba(255,253,247,0.92)", bordercolor="#b9b2a2", borderwidth=1))
    if Cwo is not None:
        _pw = (f"토양수<br>PFAS={Cwo:.3g} µg/L" if ko else f"pore water<br>Cwᵒ={Cwo:.3g} µg/L")
        ann.append(dict(x=392, y=_fy(474), text=_pw, showarrow=False,
                        xanchor="center", yanchor="middle",
                        font=dict(size=11, color="#f3ead4"), align="center",
                        bgcolor="rgba(85,58,31,0.85)", bordercolor="#835C36", borderwidth=1))
    if title is None:
        title = "식물·토양 축적 지도" if ko else "Plant + soil accumulation map"
        if t is not None:
            title += (f"  ({t:.0f}일째)" if ko else f"  (day {t:.0f})")
    # Self-contained light styling (explicit cream paper + dark ink) so the figure
    # stays readable under Streamlit dark mode too (render it with theme=None).
    fig.update_layout(
        title=dict(text=title, font=dict(color="#5a554a")), annotations=ann,
        xaxis=dict(visible=False, range=[118, 472], fixedrange=True),
        yaxis=dict(visible=False, range=[-6, 492], fixedrange=True,
                   scaleanchor="x", scaleratio=1.0),
        template="plotly_white", margin=dict(l=14, r=14, t=52, b=14),
        height=580, paper_bgcolor="#FAF7EF", plot_bgcolor="#FAF7EF",
        font=dict(color="#5a554a"))
    return fig


def fig_schematic_from_res(res, metric="conc", t_index=-1, colorscale=_HEAT, obs=None, lang="en"):
    """Convenience: build the plant map from a `model_api` result at one time index."""
    sv = api.schematic_values(res, metric, t_index)
    ko = lang == "ko"
    if ko:
        kind = "농축계수" if metric == "baf" else "농도"
        title = f"{res['congener']} — {kind} 지도  ({sv['t']:.0f}일째)"
        label = "농축계수 [L/kg]" if metric == "baf" else "농도 [µg/kg]"
    else:
        title = f"{res['congener']} — {'BAF' if metric=='baf' else 'concentration'} map  (day {sv['t']:.0f})"
        label = sv["label"]
    return fig_plant_schematic(sv["values"], cmin=sv["cmin"], cmax=sv["cmax"],
                               label=label, Cwo=sv["Cwo"], colorscale=colorscale,
                               title=title, t=sv["t"], obs=obs, lang=lang)


def fig_schematic_animated(res, metric="conc", n_frames=24, colorscale=_HEAT, lang="en"):
    """Autoplay version of the plant map: a play button + slider scrub the season,
    each compartment's colour tracking its accumulation through time."""
    ko = lang == "ko"
    ms = api.metric_series(res, metric)
    cmin, cmax = ms["cmin"], ms["cmax"]
    n_t = len(res["t"])
    idx = np.unique(np.linspace(0, n_t - 1, min(n_frames, n_t)).astype(int))

    base = fig_schematic_from_res(res, metric, int(idx[0]), colorscale, lang=lang)
    _kind = ("농축계수" if metric == "baf" else "농도") if ko else metric.upper()
    frames = []
    for ti in idx:
        sv = api.schematic_values(res, metric, int(ti))
        v = sv["values"]
        straw_only = ("stem" not in v and "leaf" not in v and v.get("straw") is not None)
        colors = {k: _color(v.get(k, v.get("straw")), cmin, cmax, colorscale)
                  for k in ("root", "stem", "leaf", "grain")}
        pts = _marker_points(v, straw_only)
        _ttl = (f"{res['congener']} — {_kind} 지도 ({res['t'][ti]:.0f}일째)" if ko
                else f"{res['congener']} — {_kind} map (day {res['t'][ti]:.0f})")
        frames.append(go.Frame(
            name=f"{res['t'][ti]:.0f}",
            data=[go.Scatter(x=[p[2] for p in pts], y=[p[3] for p in pts],
                             marker=dict(color=[p[1] for p in pts], colorscale=colorscale,
                                         cmin=cmin, cmax=cmax))],
            layout=go.Layout(shapes=_shapes_for(colors), title=_ttl)))
    base.frames = frames
    base.update_layout(
        updatemenus=[dict(type="buttons", showactive=False, x=0.02, y=1.06, xanchor="left",
                          buttons=[dict(label="▶ play", method="animate",
                                        args=[None, dict(frame=dict(duration=120, redraw=True),
                                                         fromcurrent=True, transition=dict(duration=0))]),
                                   dict(label="❚❚ pause", method="animate",
                                        args=[[None], dict(frame=dict(duration=0, redraw=False),
                                                           mode="immediate")])])],
        sliders=[dict(active=0, x=0.12, len=0.8, currentvalue=dict(prefix="day "),
                      steps=[dict(method="animate", label=f.name,
                                  args=[[f.name], dict(mode="immediate",
                                                       frame=dict(duration=0, redraw=True))])
                             for f in frames])])
    return base


# ===========================================================================
# Soil / drivers
# ===========================================================================
def fig_soil_profile(res, profile=None):
    """Soil bioavailability over time.

    profile : optional dict {'depth':(nz,), 'time':(nt,), 'conc':(nz,nt)} for a
        depth-resolved HYDRUS-1D solute field -> heatmap (depth downward). Without
        it, show the root-zone pore-water Cwᵒ(t) line + a single-row heat strip.
    """
    if profile is not None:
        fig = go.Figure(go.Heatmap(
            x=profile["time"], y=profile["depth"], z=profile["conc"],
            colorscale="YlGnBu", colorbar=dict(title="C_w [µg/L]"),
            hovertemplate="day %{x:.0f}, depth %{y:.0f} cm: %{z:.3g} µg/L<extra></extra>"))
        fig.update_layout(title="Soil pore-water profile (HYDRUS-1D)", xaxis_title="days",
                          yaxis_title="depth [cm]", yaxis=dict(autorange="reversed"), **_LAYOUT)
        return fig
    t, Cwo = res["t"], res["Cwo"]
    fig = make_subplots(rows=2, row_heights=[0.78, 0.22], shared_xaxes=True, vertical_spacing=0.04)
    fig.add_scatter(x=t, y=Cwo, name="Cwᵒ", line=dict(color="#1f77b4", width=2.5),
                    hovertemplate="day %{x:.0f}: %{y:.3g} µg/L<extra></extra>", row=1, col=1)
    fig.add_heatmap(x=t, y=[0], z=[Cwo], colorscale="YlGnBu", showscale=False,
                    hovertemplate="day %{x:.0f}: %{z:.3g} µg/L<extra></extra>", row=2, col=1)
    fig.update_layout(title="Root-zone pore-water Cwᵒ(t) (soil → plant driver)", **_LAYOUT)
    fig.update_yaxes(title_text="Cwᵒ [µg/L]", row=1, col=1)
    fig.update_yaxes(visible=False, row=2, col=1)
    fig.update_xaxes(title_text="days after transplant", row=2, col=1)
    return fig


def fig_drivers(res):
    """The three drivers ACTUALLY used by the run: Cwᵒ(t), Q_TP(t), organ mass M(t)."""
    fig = make_subplots(rows=3, cols=1, shared_xaxes=True, vertical_spacing=0.05,
                        subplot_titles=("pore water Cwᵒ [µg/L]", "transpiration Q_TP [L/day]",
                                        "organ fresh mass M [kg]"))
    t = res["t"]
    fig.add_scatter(x=t, y=res["Cwo"], name="Cwᵒ", line=dict(color="#1f77b4", width=2.5),
                    hovertemplate="Cwᵒ: %{y:.3g}<extra></extra>", row=1, col=1)
    fig.add_scatter(x=t, y=res["Qtp"], name="Q_TP", line=dict(color="royalblue", width=2.5),
                    hovertemplate="Q_TP: %{y:.3g}<extra></extra>", row=2, col=1)
    for i, k in enumerate(api.TISSUES):
        fig.add_scatter(x=t, y=res["M"][:, i], name=k, line=dict(color=_COL[k]),
                        hovertemplate=f"{k}: %{{y:.4g}} kg<extra></extra>", row=3, col=1)
    fig.update_layout(title="Drivers used by this run", hovermode="x unified",
                      showlegend=True, **_LAYOUT)
    fig.update_xaxes(title_text="days after transplant", row=3, col=1)
    return fig


def fig_cwo_profile(congener, level=1.0, profile="flooded", season=120.0,
                    n_t=121, k_leach=0.02, height=230):
    """Compact preview of the pore-water exposure shape C_w^o(t) for the chosen
    `cwo_profile`, overlaid on the flat constant baseline (both season-mean ==
    `level`, so only the temporal shape differs). Short chains leach to a decline,
    long chains stay buffered. `congener` may be None (falls back to a generic PFCA
    descriptor) for SMILES-specified compounds."""
    c = api._CONG.get(congener) if congener else None
    n_C = c["n_C"] if c else 8
    group = c["group"] if c else "PFCA"
    t = np.linspace(0.0, float(season), int(n_t))
    flat = api.cwo_profile_series(t, level, "constant")
    shaped = api.cwo_profile_series(t, level, profile, n_C=n_C, group=group,
                                    congener=congener, k_leach=k_leach)
    fig = go.Figure()
    fig.add_scatter(x=t, y=flat, name="constant", line=dict(color="#bbb", width=1.5, dash="dash"),
                    hovertemplate="constant: %{y:.3g}<extra></extra>")
    fig.add_scatter(x=t, y=shaped, name=profile, line=dict(color="#1f77b4", width=2.5),
                    fill="tonexty", fillcolor="rgba(31,119,180,0.10)",
                    hovertemplate=f"{profile}: %{{y:.3g}}<extra></extra>")
    fig.update_layout(title=f"Cwᵒ(t) preview — {profile} (mean={level:g})",
                      hovermode="x unified", height=height, showlegend=True,
                      legend=dict(orientation="h", y=1.12, x=0), **_LAYOUT)
    fig.update_yaxes(title_text="Cwᵒ [µg/L]", rangemode="tozero")
    fig.update_xaxes(title_text="days after transplant")
    return fig


def fig_isotherm(soil, Cwo_now=None):
    """Freundlich sorption isotherm S(C_w)=K_F·C_wⁿ for the paddy soil sub-model.

    `soil` is a FreundlichSoil (has K_F, n, theta_g). Marks the current operating
    point if `Cwo_now` is given. Shows why flooding (dilution) lowers bioavailability.
    """
    K_F, n = getattr(soil, "K_F", None), getattr(soil, "n", None)
    if K_F is None:                                   # redox pair -> use drained leaf
        soil = getattr(soil, "drained", soil)
        K_F, n = soil.K_F, soil.n
    cw = np.linspace(1e-3, max(Cwo_now or 1.0, 1.0) * 2.5, 200)
    fig = go.Figure()
    fig.add_scatter(x=cw, y=soil.sorbed(cw), name=f"S = {K_F:g}·C_w^{n:g}",
                    line=dict(color="#8c564b", width=2.5),
                    hovertemplate="C_w=%{x:.3g} µg/L → S=%{y:.3g} µg/kg<extra></extra>")
    if Cwo_now is not None:
        fig.add_scatter(x=[Cwo_now], y=[float(soil.sorbed(Cwo_now))], mode="markers",
                        name="operating point", marker=dict(size=12, color="#d62728"),
                        hovertemplate="now: C_w=%{x:.3g}, S=%{y:.3g}<extra></extra>")
    fig.update_layout(title="Paddy soil Freundlich isotherm (sorbed vs pore-water)",
                      xaxis_title="pore water C_w [µg/L]", yaxis_title="sorbed S [µg/kg dry]", **_LAYOUT)
    return fig


def fig_biomon_compare(measured_baf, model_baf=None):
    """Biomonitoring BAFs (measured tissue conc / pore water) vs the model, grouped bars."""
    tissues = [t for t in ("root", "stem", "straw", "leaf", "grain") if t in measured_baf]
    fig = go.Figure()
    fig.add_bar(x=tissues, y=[measured_baf[t] for t in tissues], name="measured",
                marker_color="#ff7f0e", hovertemplate="meas %{x}: %{y:.3g}<extra></extra>")
    if model_baf:
        fig.add_bar(x=tissues, y=[model_baf.get(t) for t in tissues], name="model",
                    marker_color="#1f77b4", hovertemplate="model %{x}: %{y:.3g}<extra></extra>")
    fig.update_layout(barmode="group", title="Biomonitoring BAF (measured vs model)",
                      yaxis_title="BAF [L/kg]", **_LAYOUT)
    return fig
