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
from plotly.colors import sample_colorscale

import model_api as api
import forcing_rice as fr
import growth_rice as gr

_COL = {"root": "#8c564b", "stem": "#2ca02c", "leaf": "#1f77b4", "grain": "#ff7f0e"}
_LAYOUT = dict(template="plotly_white", hoverlabel=dict(namelength=-1),
               margin=dict(l=60, r=20, t=50, b=50))
# default "accumulation heat" colour scale for the plant map (more = hotter)
_HEAT = "YlOrRd"


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


# ===========================================================================
# Plant + soil MAP  --  the model drawn to scale, compartments coloured by the
# accumulation metric (a heat colormap).  This is the "see the model" view.
# ===========================================================================
# A rice plant (Oryza sativa): a fibrous root mass in the paddy soil, arching
# culms/tillers, long slender leaf blades, and DROOPING panicles heavy with grain
# (the bent-over heads of ripe rice).  Geometry on a 0..10 (x) by -0.6..11 (y)
# canvas; each compartment is stroked/filled by the colour sampled from the
# metric scale at that compartment's value.  Paths are generated below at import.
_SOIL = dict(x0=0.3, y0=-0.4, x1=9.7, y1=3.05)
_WATER_Y = 2.75
_CROWN = (5.0, 3.0)


def _cubic(p0, p1, p2, p3, ts):
    p0, p1, p2, p3 = (np.asarray(p, float) for p in (p0, p1, p2, p3))
    return [(1 - t) ** 3 * p0 + 3 * (1 - t) ** 2 * t * p1
            + 3 * (1 - t) * t * t * p2 + t ** 3 * p3 for t in ts]


def _blade(base, tip, width=0.4, arch=0.0):
    """SVG path for an arching, lens-shaped leaf blade from `base` to `tip`."""
    base, tip = np.asarray(base, float), np.asarray(tip, float)
    d = tip - base
    L = float(np.hypot(*d)) or 1.0
    perp = np.array([-d[1], d[0]]) / L
    mid = (base + tip) / 2 + np.array([0.0, arch])
    a, b = mid + perp * width, mid - perp * width
    return (f"M {base[0]:.2f},{base[1]:.2f} Q {a[0]:.2f},{a[1]:.2f} {tip[0]:.2f},{tip[1]:.2f} "
            f"Q {b[0]:.2f},{b[1]:.2f} {base[0]:.2f},{base[1]:.2f} Z")


def _root_paths(n=11):
    """Open Bézier strokes fanning down from the crown -> a fibrous rice root mass."""
    out = []
    for xe in np.linspace(2.6, 7.4, n):
        s = (xe - 5.0) / 2.4
        ex, ey = xe, 0.05 + 0.55 * abs(s)               # outer roots end shallower
        c1 = (5.0 + s * 0.6, 2.1)
        c2 = (xe + s * 0.3, 1.1)
        out.append(f"M {_CROWN[0]:.2f},{_CROWN[1]+0.1:.2f} "
                   f"C {c1[0]:.2f},{c1[1]:.2f} {c2[0]:.2f},{c2[1]:.2f} {ex:.2f},{ey:.2f}")
    return out


# arching culms (main + two tillers) and the long leaf blades
_CULM_MAIN = "M 5.0,3.0 C 4.55,5.2 5.0,7.1 5.55,8.55"
_TILLERS = ["M 5.0,3.05 C 3.95,4.6 3.25,6.0 3.0,7.0",
            "M 5.0,3.05 C 6.05,4.5 6.85,5.6 7.2,6.6"]
_LEAVES = [_blade((4.15, 4.8), (1.2, 3.0), 0.42, -0.8),     # long lower-left, arching down
           _blade((5.5, 5.5), (8.6, 4.3), 0.40, -0.7),      # long lower-right, arching down
           _blade((4.4, 5.7), (2.2, 5.0), 0.30, -0.35),     # mid-left, medium
           _blade((5.4, 6.0), (7.3, 5.4), 0.30, -0.3)]      # mid-right, medium


def _panicle_beads(rachis, counts, hang=0.36, droop=0.14):
    """Grain ellipse centres hanging off an arching rachis (a drooping head)."""
    pts = _cubic(*rachis, np.linspace(0.06, 1.0, len(counts)))
    beads = []
    for P, nc in zip(pts, counts):
        beads.append((float(P[0]), float(P[1])))                 # node on the rachis
        for j in range(int(nc)):
            beads.append((float(P[0] - droop * j), float(P[1] - hang * (j + 1))))
    return beads


# main (right, large) + secondary (left, smaller) drooping panicles
_RACHIS1 = ((5.55, 8.55), (6.8, 9.55), (8.2, 9.15), (8.75, 7.95))
_RACHIS2 = ((3.0, 7.05), (2.45, 8.05), (1.7, 8.05), (1.35, 7.15))
_GRAIN_RACHISES = [_RACHIS1, _RACHIS2]
_GRAIN_BEADS = (_panicle_beads(_RACHIS1, [1, 1, 2, 2, 3, 3, 3, 2])
                + _panicle_beads(_RACHIS2, [1, 1, 2, 2, 2]))
# marker / label anchors (canvas coords) shared by the static + animated builders
_MARK = {"root": (5.0, 1.25), "stem": (4.55, 4.8), "leaf": (2.6, 4.0),
         "grain": (7.9, 8.25), "straw": (3.25, 5.3)}
_LABEL = {"root": (1.4, 0.6), "stem": (1.1, 6.2), "leaf": (0.9, 2.5),
          "grain": (8.9, 9.7), "straw": (0.9, 6.2)}
_NONE = "rgba(0,0,0,0)"


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
    """Return the rice-plant + soil shapes coloured by the per-part `colors`."""
    S = [dict(type="rect", **_SOIL, fillcolor="#cdb891", line=dict(color="#9c8a63", width=1), layer="below"),
         dict(type="line", x0=_SOIL["x0"], x1=_SOIL["x1"], y0=_WATER_Y, y1=_WATER_Y,
              line=dict(color="#6fa8dc", width=2, dash="dot"), layer="below")]
    # fibrous roots (open strokes, in soil)
    for pth in _root_paths():
        S.append(dict(type="path", path=pth, fillcolor=_NONE, line=dict(color=colors["root"], width=3)))
    # leaf blades (filled, behind the culms)
    for pth in _LEAVES:
        S.append(dict(type="path", path=pth, fillcolor=colors["leaf"], line=dict(color=line, width=1)))
    # culms / tillers (thick strokes)
    S.append(dict(type="path", path=_CULM_MAIN, fillcolor=_NONE, line=dict(color=colors["stem"], width=11)))
    for pth in _TILLERS:
        S.append(dict(type="path", path=pth, fillcolor=_NONE, line=dict(color=colors["stem"], width=7)))
    # panicle rachises (thin) + drooping grain beads
    for rc in _GRAIN_RACHISES:
        S.append(dict(type="path", fillcolor=_NONE, line=dict(color=colors["grain"], width=2.5),
                      path=f"M {rc[0][0]:.2f},{rc[0][1]:.2f} C {rc[1][0]:.2f},{rc[1][1]:.2f} "
                           f"{rc[2][0]:.2f},{rc[2][1]:.2f} {rc[3][0]:.2f},{rc[3][1]:.2f}"))
    for (bx, by) in _GRAIN_BEADS:
        S.append(dict(type="circle", x0=bx - 0.17, y0=by - 0.27, x1=bx + 0.17, y1=by + 0.27,
                      fillcolor=colors["grain"], line=dict(color=line, width=0.7)))
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
    fig.add_scatter(
        x=[p[2] for p in pts], y=[p[3] for p in pts], mode="markers",
        marker=dict(size=15, color=[p[1] for p in pts], colorscale=colorscale,
                    cmin=cmin, cmax=cmax, line=dict(color="#333", width=1),
                    colorbar=dict(title=label, thickness=16, len=0.85)),
        text=[p[0] for p in pts], customdata=[p[1] for p in pts],
        hovertemplate="%{text}: %{customdata:.3g}<extra></extra>", showlegend=False)

    ann = []
    for name, v, mx, my in pts:
        lx, ly = _LABEL.get(name, (mx, my))
        txt = f"<b>{name}</b><br>{v:.3g}"
        if obs and name in obs and obs[name] is not None:
            txt += f"<br><span style='color:#888'>obs {obs[name]:.3g}</span>"
        ann.append(dict(x=mx, y=my, ax=lx, ay=ly, axref="x", ayref="y",
                        text=txt, showarrow=True, arrowhead=0, arrowcolor="#999",
                        font=dict(size=12), align="center",
                        bgcolor="rgba(255,255,255,0.75)", bordercolor="#ccc", borderwidth=1))
    if Cwo is not None:
        ann.append(dict(x=8.6, y=2.5, text=f"pore water<br>Cwᵒ={Cwo:.3g} µg/L", showarrow=False,
                        font=dict(size=11, color="#33597a"), align="center",
                        bgcolor="rgba(214,233,247,0.85)", bordercolor="#9cc3e0", borderwidth=1))
    if title is None:
        title = "Plant + soil accumulation map"
        if t is not None:
            title += f"  (day {t:.0f})"
    fig.update_layout(
        title=title, annotations=ann,
        xaxis=dict(visible=False, range=[0, 10], fixedrange=True),
        yaxis=dict(visible=False, range=[-0.6, 11], fixedrange=True,
                   scaleanchor="x", scaleratio=1.0),
        template="plotly_white", margin=dict(l=10, r=10, t=50, b=10),
        height=560, plot_bgcolor="rgba(235,244,250,0.5)")
    return fig


def fig_schematic_from_res(res, metric="conc", t_index=-1, colorscale=_HEAT, obs=None):
    """Convenience: build the plant map from a `model_api` result at one time index."""
    sv = api.schematic_values(res, metric, t_index)
    title = f"{res['congener']} — {'BAF' if metric=='baf' else 'concentration'} map  (day {sv['t']:.0f})"
    return fig_plant_schematic(sv["values"], cmin=sv["cmin"], cmax=sv["cmax"],
                               label=sv["label"], Cwo=sv["Cwo"], colorscale=colorscale,
                               title=title, t=sv["t"], obs=obs)


def fig_schematic_animated(res, metric="conc", n_frames=36, colorscale=_HEAT):
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
