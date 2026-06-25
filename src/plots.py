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
from plotly.subplots import make_subplots
from plotly.colors import sample_colorscale

import model_api as api
import forcing_rice as fr
import growth_rice as gr

_COL = {"root": "#8c564b", "stem": "#2ca02c", "leaf": "#1f77b4", "grain": "#ff7f0e"}
# plain-language tissue names for the general-audience (Simple) views
_PLAIN = {"root": "Roots", "stem": "Stems", "leaf": "Leaves", "grain": "Grain", "straw": "Straw"}
_LAYOUT = dict(template="plotly_white", hoverlabel=dict(namelength=-1),
               margin=dict(l=60, r=20, t=50, b=50))
# default "accumulation heat" colour scale for the plant map (more = hotter)
_HEAT = "YlOrRd"


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


def fig_baf(res, obs, extra=None):
    """Predicted vs observed (Yamazaki) root/straw/grain BAF, grouped bars.

    `extra` optionally overlays additional model series (e.g. the EXPLORATORY
    two-pool model): a dict {label: {"root","straw","grain"}} (values may be None);
    each is added as its own grouped bar so the canonical core, the overlay model(s)
    and the observed data are compared side by side.
    """
    tis = ["root", "straw", "grain"]
    pred = [res["baf_final"]["root"], res["straw_baf"], res["baf_final"]["grain"]]
    fig = go.Figure()
    fig.add_bar(x=tis, y=pred, name="model (4-pool core)", marker_color="#1f77b4",
                hovertemplate="core %{x}: %{y:.3g}<extra></extra>")
    for i, (label, vals) in enumerate(dict(extra or {}).items()):
        fig.add_bar(x=tis, y=[vals.get(t_) for t_ in tis], name=label,
                    marker_color=_BAF_EXTRA_COLORS[i % len(_BAF_EXTRA_COLORS)],
                    marker_line=dict(width=0.5, color="#333"),
                    hovertemplate=label + " %{x}: %{y:.3g}<extra></extra>")
    if obs:
        fig.add_bar(x=tis, y=[obs.get(t_, None) for t_ in tis], name="Yamazaki 2023",
                    marker_color="#ff7f0e", hovertemplate="obs %{x}: %{y:.3g}<extra></extra>")
    fig.update_layout(barmode="group", title=f"{res['congener']} — predicted vs observed BAF",
                      yaxis_title="BAF [L/kg]", **_LAYOUT)
    return fig


def fig_buildup_plain(res):
    """Plain-language tissue concentration over the season (Simple-mode view).

    Same data as `fig_tissue` but with friendly tissue names (Roots/Stems/Leaves/
    Grain) and jargon-free axis/title text, for a general audience. Reuses the
    grain pre-formation mask so the empty pre-flowering period is not drawn."""
    fig = go.Figure()
    for j, tis in enumerate(api.TISSUES):
        fig.add_scatter(x=res["t"], y=_formed(res, j, res["conc"][tis]), name=_PLAIN[tis],
                        mode="lines", line=dict(width=2.8, color=_COL[tis]),
                        hovertemplate=f"{_PLAIN[tis]}: %{{y:.3g}} µg/kg<extra></extra>")
    fig.update_layout(title="How PFAS builds up in the plant over the season",
                      xaxis_title="days after the rice is transplanted",
                      yaxis_title="PFAS in the plant tissue [µg per kg]",
                      hovermode="x unified", **_LAYOUT)
    return fig


def fig_where_plain(res):
    """Plain-language bar of the final PFAS level in roots / straw / grain.

    A jargon-free read of where the chemical ends up at harvest (the same numbers
    as the BAF bars, but labelled as a concentration build-up, no 'BAF' symbol)."""
    order = ["root", "straw", "grain"]
    vals = {"root": res["conc"]["root"][-1], "straw": res["straw"][-1],
            "grain": res["conc"]["grain"][-1]}
    fig = go.Figure(go.Bar(
        x=[_PLAIN[t_] for t_ in order], y=[vals[t_] for t_ in order],
        marker_color=[_COL.get(t_, "#1f77b4") for t_ in order],
        text=[f"{vals[t_]:.2g}" for t_ in order], textposition="outside",
        hovertemplate="%{x}: %{y:.3g} µg/kg<extra></extra>"))
    fig.update_layout(title=f"Where {res['congener']} ends up in the plant (at harvest)",
                      yaxis_title="PFAS in the tissue [µg per kg]",
                      xaxis_title="part of the rice plant", **_LAYOUT)
    return fig


def fig_exposure_posterior(est):
    """Posterior over the estimated soil-water contamination level Cwᵒ [µg/L].

    From `model_api.estimate_exposure_bayesian`: a shaded probability curve vs Cwᵒ
    (log x), with the 95% credible interval darker and the median (most likely
    value) marked. Plain labels for a general audience."""
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
                  annotation_text=f"most likely {med:.3g} µg/L", annotation_position="top")
    fig.update_layout(title="Estimated contamination level in the soil water (with uncertainty)",
                      xaxis_title="PFAS dissolved in the soil water [µg/L]",
                      yaxis_title="relative probability", xaxis_type="log",
                      showlegend=False, **_LAYOUT)
    fig.update_yaxes(rangemode="tozero", showticklabels=False)
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
_K = 800.0                       # y-flip constant (yf = _K - y_svg), keeps yf > 0
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


# --- the SVG silhouette (raw illustration coords), flipped to y-up at import ----
_SOIL_PATH = _flip_path("M62 590 Q 158 576, 250 587 Q 342 597, 432 583 Q 512 573, 558 589 "
                        "L 558 662 L 62 662 Z")
_SOIL_BAND = _flip_path("M62 632 L 558 632 L 558 662 L 62 662 Z")
_SOIL_TOP = _flip_path("M62 590 Q 158 576, 250 587 Q 342 597, 432 583 Q 512 573, 558 589")
_SOIL_TEX = [(140, 612, 3, 2), (210, 640, 2.5, 1.8), (330, 615, 3, 2), (400, 645, 2.5, 1.8),
             (470, 620, 3, 2), (500, 648, 2, 1.5), (100, 642, 2.5, 1.8), (260, 650, 2, 1.5),
             (360, 656, 2.5, 1.8), (440, 610, 2, 1.5)]
_CROWN = (300, 585, 11, 6)                                   # cx, cy, rx, ry
_ROOTS = [(_flip_path(p), w) for p, w in (
    ("M300 586 C 296 606, 284 624, 264 648", 2.6),
    ("M300 586 C 300 610, 301 632, 297 654", 2.6),
    ("M300 586 C 306 606, 320 624, 340 648", 2.6),
    ("M299 588 C 290 604, 274 616, 250 632", 1.8),
    ("M301 588 C 313 602, 331 612, 354 626", 1.8),
    ("M300 590 C 295 612, 289 634, 281 656", 2.2),
    ("M300 590 C 305 612, 313 634, 323 656", 2.2),
    ("M283 624 C 276 630, 270 638, 263 646", 1.1),
    ("M320 624 C 328 630, 334 638, 342 648", 1.1),
    ("M297 632 C 292 640, 290 648, 288 656", 1.1))]
_LEAVES = [_flip_path("M288 548 C 281 487, 183 366, 145 392 C 194 380, 295 489, 294 549 Z"),
           _flip_path("M290 558 C 298 488, 396 392, 452 422 C 405 400, 308 494, 298 559 Z"),
           _flip_path("M286 538 C 294 456, 426 312, 494 345 C 436 322, 304 464, 294 539 Z")]
_LEAF_VEINS = [_flip_path("M288 548 C 281 487, 183 366, 145 392"),
               _flip_path("M290 558 C 298 488, 396 392, 452 422"),
               _flip_path("M286 538 C 294 456, 426 312, 494 345")]
_STEM_PATH = _flip_path("M300 585 C 266 460, 260 290, 318 178")
_NODES = [(293, 560, 6, 4, -8), (279, 471, 5.5, 3.8, 0), (276, 342, 5, 3.6, 20)]
_RACHIS = [_flip_path("M318 178 C 332 150, 356 148, 376 172 C 398 198, 406 250, 396 306"),
           _flip_path("M340 156 C 356 166, 364 184, 360 208")]
_GRAINS = [(324, 172, 5, 9, 30), (332, 160, 5, 9, 25), (343, 152, 5, 9, 15), (356, 150, 5, 9, 5),
           (369, 154, 5, 9, -10), (380, 164, 5, 9, -25), (388, 178, 5, 9, -35), (394, 196, 5, 9, -30),
           (397, 216, 5, 9, -20), (398, 238, 5, 9, -12), (398, 260, 5, 9, -8), (397, 282, 5, 9, -5),
           (395, 302, 5, 9, 0), (318, 162, 4.5, 8, 28), (348, 146, 4.5, 8, 4), (374, 150, 4.5, 8, -18),
           (386, 230, 4.5, 8, -6), (385, 266, 4.5, 8, -4), (384, 296, 4.5, 8, 0)]
_GRAIN_PATHS = [_ellipse_path(*g) for g in _GRAINS]

# marker / label anchors (SVG coords -> flipped) shared by static + animated builders
_MARK = {"root": (300, _fy(612)), "stem": (264, _fy(430)), "leaf": (210, _fy(455)),
         "grain": (395, _fy(240)), "straw": (240, _fy(470))}
# label-box anchors (SVG coords -> flipped) + grow-inward (xanchor,yanchor) so the
# boxes always stay inside the [20,600]x[120,680] frame regardless of side.
_LABEL = {"root":  (150, _fy(645), "center", "bottom"),
          "stem":  (45,  _fy(478), "left",   "middle"),
          "leaf":  (45,  _fy(330), "left",   "middle"),
          "grain": (560, _fy(168), "right",  "top"),
          "straw": (45,  _fy(500), "left",   "middle")}


def _frac(v, cmin, cmax):
    if v is None or not np.isfinite(v):
        return None
    if cmax <= cmin:
        return 0.5
    return float(np.clip((v - cmin) / (cmax - cmin), 0.0, 1.0))


def _color(v, cmin, cmax, scale, nan="#dcdcdc"):
    f = _frac(v, cmin, cmax)
    return nan if f is None else sample_colorscale(scale, [f])[0]


def _shapes_for(colors, line="#4f4f4f"):
    """Return the rice-plant + paddy-soil shapes coloured by the per-organ `colors`.

    Soil cross-section, rachis and nodes are fixed structural colours; the root
    (crown + fibrous roots), leaves, stem and grain ear are filled/stroked with the
    metric-scale colour for that compartment's value (the accumulation map)."""
    S = [dict(type="path", path=_SOIL_PATH, fillcolor="#6B4A2A", line=dict(width=0), layer="below"),
         dict(type="path", path=_SOIL_BAND, fillcolor="#553A1F", opacity=0.85, line=dict(width=0), layer="below"),
         dict(type="path", path=_SOIL_TOP, fillcolor=_NONE, line=dict(color="#835C36", width=1.4), layer="below")]
    for (cx, cy, rx, ry) in _SOIL_TEX:
        S.append(dict(type="circle", x0=cx - rx, x1=cx + rx, y0=_fy(cy) - ry, y1=_fy(cy) + ry,
                      fillcolor="#5A3E22", line=dict(width=0), layer="below"))
    # root crown + fibrous roots (root colour)
    cx, cy, rx, ry = _CROWN
    S.append(dict(type="circle", x0=cx - rx, x1=cx + rx, y0=_fy(cy) - ry, y1=_fy(cy) + ry,
                  fillcolor=colors["root"], line=dict(color=line, width=0.8)))
    for pth, w in _ROOTS:
        S.append(dict(type="path", path=pth, fillcolor=_NONE, line=dict(color=colors["root"], width=w)))
    # stem C-curve (drawn BEFORE the leaves so the leaves sit in front): edge + fill + nodes
    S.append(dict(type="path", path=_STEM_PATH, fillcolor=_NONE, line=dict(color=line, width=8)))
    S.append(dict(type="path", path=_STEM_PATH, fillcolor=_NONE, line=dict(color=colors["stem"], width=7)))
    for (cx, cy, rx, ry, rot) in _NODES:
        S.append(dict(type="path", path=_ellipse_path(cx, cy, rx, ry, rot, 14),
                      fillcolor="#33401f", line=dict(width=0)))
    # leaves (leaf colour) + a subtle midrib -- in front of the stem
    for pth in _LEAVES:
        S.append(dict(type="path", path=pth, fillcolor=colors["leaf"], line=dict(color=line, width=1)))
    for pth in _LEAF_VEINS:
        S.append(dict(type="path", path=pth, fillcolor=_NONE, line=dict(color="rgba(35,35,35,0.22)", width=1)))
    # ear: rachis (structural) + nodding grains (grain colour)
    for pth in _RACHIS:
        S.append(dict(type="path", path=pth, fillcolor=_NONE, line=dict(color="#9c7a2a", width=2.0)))
    for pth in _GRAIN_PATHS:
        S.append(dict(type="path", path=pth, fillcolor=colors["grain"], line=dict(color="#C8881A", width=0.5)))
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
                        obs=None):
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

    ann = []
    for name, v, mx, my in pts:
        lx, ly, xa, ya = _LABEL.get(name, (mx, my, "center", "middle"))
        txt = f"<b>{name}</b><br>{v:.3g}"
        if obs and name in obs and obs[name] is not None:
            txt += f"<br><span style='color:#8a8270'>obs {obs[name]:.3g}</span>"
        ann.append(dict(x=mx, y=my, ax=lx, ay=ly, axref="x", ayref="y",
                        xanchor=xa, yanchor=ya,
                        text=txt, showarrow=True, arrowhead=0, arrowcolor="#9a9384", arrowwidth=1,
                        font=dict(size=12, color=_ink), align="center",
                        bgcolor="rgba(255,253,247,0.92)", bordercolor="#b9b2a2", borderwidth=1))
    if Cwo is not None:
        ann.append(dict(x=476, y=_fy(610), text=f"pore water<br>Cwᵒ={Cwo:.3g} µg/L", showarrow=False,
                        xanchor="center", yanchor="middle",
                        font=dict(size=11, color="#f3ead4"), align="center",
                        bgcolor="rgba(85,58,31,0.85)", bordercolor="#835C36", borderwidth=1))
    if title is None:
        title = "Plant + soil accumulation map"
        if t is not None:
            title += f"  (day {t:.0f})"
    # Self-contained light styling (explicit cream paper + dark ink) so the figure
    # stays readable under Streamlit dark mode too (render it with theme=None).
    fig.update_layout(
        title=dict(text=title, font=dict(color="#5a554a")), annotations=ann,
        xaxis=dict(visible=False, range=[20, 600], fixedrange=True),
        yaxis=dict(visible=False, range=[120, 680], fixedrange=True,
                   scaleanchor="x", scaleratio=1.0),
        template="plotly_white", margin=dict(l=14, r=14, t=52, b=14),
        height=580, paper_bgcolor="#FAF7EF", plot_bgcolor="#FAF7EF",
        font=dict(color="#5a554a"))
    return fig


def fig_schematic_from_res(res, metric="conc", t_index=-1, colorscale=_HEAT, obs=None):
    """Convenience: build the plant map from a `model_api` result at one time index."""
    sv = api.schematic_values(res, metric, t_index)
    title = f"{res['congener']} — {'BAF' if metric=='baf' else 'concentration'} map  (day {sv['t']:.0f})"
    return fig_plant_schematic(sv["values"], cmin=sv["cmin"], cmax=sv["cmax"],
                               label=sv["label"], Cwo=sv["Cwo"], colorscale=colorscale,
                               title=title, t=sv["t"], obs=obs)


def fig_schematic_animated(res, metric="conc", n_frames=24, colorscale=_HEAT):
    """Autoplay version of the plant map: a play button + slider scrub the season,
    each compartment's colour tracking its accumulation through time."""
    ms = api.metric_series(res, metric)
    cmin, cmax = ms["cmin"], ms["cmax"]
    n_t = len(res["t"])
    idx = np.unique(np.linspace(0, n_t - 1, min(n_frames, n_t)).astype(int))

    base = fig_schematic_from_res(res, metric, int(idx[0]), colorscale)
    frames = []
    for ti in idx:
        sv = api.schematic_values(res, metric, int(ti))
        v = sv["values"]
        straw_only = ("stem" not in v and "leaf" not in v and v.get("straw") is not None)
        colors = {k: _color(v.get(k, v.get("straw")), cmin, cmax, colorscale)
                  for k in ("root", "stem", "leaf", "grain")}
        pts = _marker_points(v, straw_only)
        frames.append(go.Frame(
            name=f"{res['t'][ti]:.0f}",
            data=[go.Scatter(x=[p[2] for p in pts], y=[p[3] for p in pts],
                             marker=dict(color=[p[1] for p in pts], colorscale=colorscale,
                                         cmin=cmin, cmax=cmax))],
            layout=go.Layout(shapes=_shapes_for(colors),
                             title=f"{res['congener']} — {metric.upper()} map (day {res['t'][ti]:.0f})")))
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
