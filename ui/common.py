"""Shared UI layer for the PFAS-rice dashboard: constants, cached model helpers,
and the render building blocks that app.py assembles. Split out of the monolithic
app.py (HANDOFF P3-1); behaviour is unchanged."""
import os

import numpy as np
import streamlit as st
import streamlit.components.v1 as components

import model_api as api
import plots
from ui import i18n
from ui.i18n import t as _t

# repo root (this file lives in <root>/ui/), and the bundled example CSVs
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_EX = os.path.join(_ROOT, "examples")


APP_VERSION = "1.0 (general-audience UI)"
REPO_URL = "https://github.com/ccy5123/pfas-rice-model"
DOCS_URL = "https://github.com/ccy5123/pfas-rice-model/tree/main/docs"

# Bilingual UI copy now lives in ui/i18n.py (HANDOFF P3-2). These module-level names
# are kept as thin views onto the table so the rest of the UI (sidebar/simple/expert)
# imports them unchanged.
_DISCLAIMER = _t("disclaimer", "en")         # top banner + footer (Expert / English)
_DISCLAIMER_KO = _t("disclaimer", "ko")      # top banner + footer (Simple / Korean)

# Friendlier congener names for the dropdowns (value stays the symbol).
_FRIENDLY_CONG = i18n.CONGENER_LABELS["en"]
_FRIENDLY_CONG_KO = i18n.CONGENER_LABELS["ko"]


def _cong_label(name):
    return _FRIENDLY_CONG.get(name, name)


def _cong_label_ko(name):
    return _FRIENDLY_CONG_KO.get(name, name)


# Plain "how contaminated?" presets → pore-water concentration [µg/L]. Medium =
# 1.0 µg/L keeps tissue conc == build-up factor (the model's reference point). The
# Korean variant carries the short word (낮은/중간/높은) reused in the headline.
_PRESETS = i18n.PRESETS["en"]
_PRESETS_KO = i18n.PRESETS["ko"]

# One-click policy story scenarios (Simple mode) -> (congener, Cwᵒ, word, desc).
_SCENARIOS_KO = i18n.SCENARIOS_KO


# ---------------------------------------------------------------- helpers
def _nearest_index(t, day):
    return int(np.argmin(np.abs(np.asarray(t) - day)))


@st.cache_data(show_spinner=False)
def _simulate(congener, **kw):
    """Cache model runs (drivers passed as a hashable tuple, rebuilt here)."""
    drv = kw.pop("drivers_tuple", None)
    if drv is not None:
        t, Cwo, Qtp, Mflat, ncol = drv
        kw["drivers"] = dict(t=np.array(t), Cwo=np.array(Cwo), Qtp=np.array(Qtp),
                             M=np.array(Mflat).reshape(-1, ncol))
    kl = kw.pop("cwo_k_leach", None)                       # scalar (hashable) -> cwo_kw dict
    if kl is not None and kw.get("cwo_profile", "constant") != "constant":
        kw["cwo_kw"] = {"k_leach": float(kl)}
    return api.simulate(congener, **kw)


def _drivers_tuple(d):
    return (tuple(d["t"]), tuple(d["Cwo"]), tuple(d["Qtp"]),
            tuple(np.asarray(d["M"]).ravel()), int(np.asarray(d["M"]).shape[1]))


@st.cache_data(show_spinner=False)
def _simulate_twopool_seq(congener):
    """Cache the EXPLORATORY sequestration two-pool run (curated congener only).

    Run at the model's calibrated operating point (defaults: Cwo=1, season=120,
    demo forcings) so it reproduces the documented Yamazaki headline -- an
    apples-to-apples reference against the fixed observed bars. Returns the
    root/straw/grain BAF, or None if the run fails."""
    try:
        r = api.simulate_twopool_seq(congener)
        return {"root": r["baf_final"]["root"], "straw": r["straw_baf"],
                "grain": r["baf_final"]["grain"]}
    except Exception:                                        # noqa: BLE001
        return None


@st.cache_data(show_spinner="Parameterising structure (RDKit)…")
def _simulate_smiles(smiles, **kw):
    """Cache a SMILES (structure) run: RDKit → descriptors → Compound → full ODE."""
    drv = kw.pop("drivers_tuple", None)
    if drv is not None:
        t, Cwo, Qtp, Mflat, ncol = drv
        kw["drivers"] = dict(t=np.array(t), Cwo=np.array(Cwo), Qtp=np.array(Qtp),
                             M=np.array(Mflat).reshape(-1, ncol))
    kl = kw.pop("cwo_k_leach", None)                       # scalar (hashable) -> cwo_kw dict
    if kl is not None and kw.get("cwo_profile", "constant") != "constant":
        kw["cwo_kw"] = {"k_leach": float(kl)}
    return api.simulate_from_smiles(smiles, **kw)


@st.cache_data(show_spinner=False)
def _mol_svg(smiles, w=290, h=170):
    """2-D structure as an SVG string (RDKit, Cairo-free → works on Streamlit Cloud).
    Returns (svg, None) on success or (None, reason) so the UI can show why it failed."""
    try:
        from rdkit import Chem
        from rdkit.Chem.Draw import rdMolDraw2D
    except Exception as e:                                   # noqa: BLE001
        return None, f"RDKit import failed: {e}"
    try:
        m = Chem.MolFromSmiles(smiles)
        if m is None:
            return None, "RDKit could not parse this SMILES."
        d = rdMolDraw2D.MolDraw2DSVG(w, h)
        d.drawOptions().padding = 0.12
        d.DrawMolecule(m)
        d.FinishDrawing()
        return d.GetDrawingText(), None
    except Exception as e:                                   # noqa: BLE001
        return None, f"draw error ({type(e).__name__}): {e}"


@st.cache_data(show_spinner="Running HYDRUS-1D…")
def _hydrus_drivers_cached(congener, season, f_oc, flood_until, percolation, biomass="oryza"):
    """Cache a real HYDRUS-1D paddy run (a few seconds) per parameter set."""
    drv, _ = api.hydrus_drivers(congener, season=season, f_oc=f_oc, biomass=biomass,
                                flood_until=float(flood_until), percolation=float(percolation))
    return drv


@st.cache_data(show_spinner=False)
def _hydrus_soil_congener(smiles):
    """Map a SMILES to a curated congener for the HYDRUS soil Kd: the read-across
    match if known, else the nearest curated congener by perfluoro-C in the same
    head-group family. Returns (name, how) or (None, None) if unparseable."""
    try:
        from pfas_structure import descriptors
        d = descriptors(smiles)
    except Exception:                                       # noqa: BLE001
        return None, None
    if d.matched_name:
        return d.matched_name, "match"
    fam = "PFSA" if d.head_group == "sulfonate" else ("ether" if d.n_ether_O else "PFCA")
    cands = [c for c in api._CONG.values() if c["group"] == fam] \
        or [c for c in api._CONG.values() if c["group"] == "PFCA"]
    npfc = lambda c: c["n_C"] - 1 if c["group"] == "PFCA" else c["n_C"]   # noqa: E731
    best = min(cands, key=lambda c: abs(npfc(c) - d.n_perfluoroC))
    return best["name"], "nearest"


def _png_bytes(fig, scale=2):
    """Static PNG of a Plotly figure via kaleido. Returns (bytes, None) or (None, reason)
    so the UI degrades gracefully when kaleido (and its Chrome) are not installed."""
    try:
        return fig.to_image(format="png", scale=scale), None
    except Exception as e:                                   # noqa: BLE001
        return None, f"{type(e).__name__}: {e}"


def _html_bytes(fig):
    """Self-contained interactive HTML of a Plotly figure. Needs NO kaleido/Chrome,
    so it is the always-available export fallback (keeps hover/zoom; loads plotly.js
    from the CDN)."""
    try:
        return fig.to_html(full_html=True, include_plotlyjs="cdn").encode("utf-8"), None
    except Exception as e:                                   # noqa: BLE001
        return None, f"{type(e).__name__}: {e}"


def _summary_html(cfg, lang="ko"):
    """A self-contained one-page HTML summary for a slide / handout (download).

    Combines the plant accumulation map, the EFSA-TWI intake gauge, the key tissue
    numbers (with the coarse a-priori band) and the caveats/sources into a single
    offline HTML file (Plotly embedded, so it works without a network). Returns
    (bytes, None) or (None, reason) so the UI can degrade gracefully."""
    ko = lang != "en"
    try:
        res = cfg.res
        cong = cfg.congener
        scenario = getattr(cfg, "preset_label", None) or getattr(cfg, "scen_label", None)
        grain_c = float(res["conc"]["grain"][-1])
        info = api.intake_fraction(grain_c, congener=cong)
        map_fig = plots.fig_schematic_from_res(res, "conc", -1, lang=lang)
        gauge = plots.fig_intake_gauge(info, lang=lang)
        # embed plotly.js once (offline-safe), then the gauge without re-embedding
        map_div = map_fig.to_html(full_html=False, include_plotlyjs=True,
                                  default_height="480px")
        gauge_div = gauge.to_html(full_html=False, include_plotlyjs=False,
                                  default_height="300px")
        rows = [("뿌리" if ko else "Roots", float(res["conc"]["root"][-1])),
                ("짚(줄기+잎)" if ko else "Straw", float(res["straw"][-1])),
                ("낟알(먹는 쌀)" if ko else "Grain", grain_c)]
        trs = ""
        for name, v in rows:
            b = api.predictive_band(v)
            trs += (f"<tr><td>{name}</td><td>{v:.2g}</td>"
                    f"<td>{b['lo']:.2g}–{b['hi']:.2g}</td></tr>")
        efsa = ""
        if np.isfinite(info.get("percent", float("nan"))):
            grp = ("" if info.get("in_group")
                   else ("(참고 — EFSA 4종 합 미포함)" if ko else "(reference — not in the EFSA four)"))
            efsa = (f"<p><b>EFSA 안전기준(TWI) 대비 약 {info['percent']:.0f}%</b> {grp} "
                    f"— 하루 {info['rice_intake_g_day']:.0f} g 쌀 섭취 가정.</p>" if ko else
                    f"<p><b>~{info['percent']:.0f}% of the EFSA TWI</b> {grp} "
                    f"— assuming {info['rice_intake_g_day']:.0f} g rice/day.</p>")
        disc = _DISCLAIMER_KO if ko else _DISCLAIMER
        disc = disc.replace("**", "")
        title = (f"PFAS–벼 흡수 요약 — {cong}" if ko else f"PFAS–Rice uptake summary — {cong}")
        scen_line = (f"<p style='color:#555'>시나리오: {scenario}</p>" if (ko and scenario)
                     else (f"<p style='color:#555'>Scenario: {scenario}</p>" if scenario else ""))
        th = ("부위,예측 농도 (µg/kg),대략 범위" if ko else "Tissue,Predicted (µg/kg),Rough range").split(",")
        src = ("출처: EFSA 2020 그룹 TWI 4.4 ng/kg 체중/주 (doi:10.2903/j.efsa.2020.6223) · "
               "쌀 섭취량 KOSIS/국민건강영양조사 · 실측 Yamazaki et al. 2023." if ko else
               "Sources: EFSA 2020 TWI 4.4 ng/kg bw/week (doi:10.2903/j.efsa.2020.6223); "
               "rice intake KOSIS/KNHANES; observed Yamazaki et al. 2023.")
        html = (
            f"<!doctype html><html lang='{'ko' if ko else 'en'}'><head><meta charset='utf-8'>"
            f"<meta name='viewport' content='width=device-width,initial-scale=1'>"
            f"<title>{title}</title><style>"
            "body{font-family:'Malgun Gothic',-apple-system,Segoe UI,Arial,sans-serif;"
            "max-width:960px;margin:24px auto;padding:0 16px;color:#222}"
            "h1{font-size:1.5rem;margin:.2rem 0}"
            "table{border-collapse:collapse;margin:8px 0}"
            "td,th{border:1px solid #ccc;padding:5px 12px;text-align:left}"
            "th{background:#f5f5f5}"
            ".row{display:flex;gap:24px;flex-wrap:wrap;align-items:flex-start}"
            ".caveat{background:#fff8e1;border-left:4px solid #f0ad4e;padding:10px 14px;"
            "border-radius:4px;margin:12px 0;font-size:.92rem}"
            ".src{color:#666;font-size:.82rem;margin-top:10px}</style></head><body>"
            f"<h1>🌾 {title}</h1>{scen_line}"
            f"<div class='row'><div style='flex:2;min-width:340px'>{map_div}</div>"
            f"<div style='flex:1;min-width:280px'>{gauge_div}"
            f"<table><tr><th>{th[0]}</th><th>{th[1]}</th><th>{th[2]}</th></tr>{trs}</table>"
            f"{efsa}</div></div>"
            f"<div class='caveat'>⚠ {disc}</div>"
            f"<p class='src'>{src}</p></body></html>")
        return html.encode("utf-8"), None
    except Exception as e:                                   # noqa: BLE001
        return None, f"{type(e).__name__}: {e}"


def _glossary_md(ko=False):
    """Plain-language glossary (rendered in the About tab and a Simple-mode expander)."""
    return _t("glossary", "ko" if ko else "en")


# uncertainty presets for the inverse estimator (measurement + model noise, log10 units)
_UNC = {"Typical (±~40%)": 0.15, "High precision (±~20%)": 0.10, "Rough (±~2×)": 0.30}


def _render_inverse_estimator(congener, *, E_m_mV, f_xy_source, biomass, key, simple=True):
    """Shared 'work backwards' panel: Bayesian estimate of the soil-water contamination
    level Cwᵒ from measured tissue concentrations, with a credible interval. Used by
    both the Simple (Korean) and Expert (English) tabs. `key` namespaces the widgets."""
    lang = "ko" if simple else "en"
    if congener not in api.CONGENERS:
        st.info(_t("inv.not_curated", lang))
        return
    st.markdown(_t("inv.intro", lang))
    c1, c2, c3 = st.columns(3)
    root = c1.number_input(_t("inv.in_root", lang), 0.0, 1e6, 0.0, 0.1, key=f"{key}_root")
    straw = c2.number_input(_t("inv.in_straw", lang), 0.0, 1e6, 0.0, 0.1, key=f"{key}_straw")
    grain = c3.number_input(_t("inv.in_grain", lang), 0.0, 1e6, 0.0, 0.1, key=f"{key}_grain")
    unc_label = st.radio(_t("inv.precision_label", lang), list(_UNC), horizontal=True,
                         key=f"{key}_unc", help=_t("inv.precision_help", lang))
    sigma = _UNC[unc_label]
    have = any(v > 0 for v in (root, straw, grain))
    run = st.button(_t("inv.estimate_btn", lang), key=f"{key}_btn", disabled=not have, type="primary")
    sig = (congener, root, straw, grain, sigma, E_m_mV, f_xy_source, biomass)
    if run:
        st.session_state[f"{key}_sig"] = sig
    if not have:
        st.caption(_t("inv.enter_first", lang))
        return
    if st.session_state.get(f"{key}_sig") != sig:
        st.caption(_t("inv.press_estimate", lang))
        return
    # Result cache in session_state (keyed by sig) so reruns are instant; the slow
    # first compute (~8 ODE solves) shows a live step-by-step progress bar so it
    # never looks frozen (HANDOFF P2-1).
    res_key = f"{key}_result"
    stored = st.session_state.get(res_key)
    if stored is not None and stored[0] == sig:
        est = stored[1]
    else:
        meas = {k: v for k, v in (("root", root), ("straw", straw), ("grain", grain))
                if v is not None and v > 0}
        prog = st.progress(0.0, text=_t("inv.preparing", lang))

        def _cb(done, total):
            prog.progress(min(done / total, 1.0), text=_t("inv.running", lang, done=done, total=total))
        try:
            est = api.estimate_exposure_bayesian(
                congener, meas, sigma_log10=sigma, E_m_mV=E_m_mV,
                f_xy_source=f_xy_source, biomass=biomass, progress=_cb)
        except Exception as e:                               # noqa: BLE001
            prog.empty()
            st.error(_t("inv.error", lang, e=e))
            return
        prog.empty()
        st.session_state[res_key] = (sig, est)
    med = est["median"]
    lo, hi = est["ci95"]
    mc1, mc2 = st.columns([1, 2])
    delta = _t("inv.range95", lang, lo=lo, hi=hi) if np.isfinite(lo) else _t("inv.range_unconstrained", lang)
    mc1.metric(_t("inv.metric_label", lang), f"{med:.3g} µg/L", delta, delta_color="off")
    summary = (_t("inv.summary_lead", lang, congener=congener, med=med)
               + (_t("inv.summary_ci", lang, lo=lo, hi=hi) if np.isfinite(lo)
                  else _t("inv.summary_no_ci", lang))
               + _t("inv.summary_tail", lang))
    mc2.markdown(summary)
    st.plotly_chart(plots.fig_exposure_posterior(est, lang=lang), width="stretch")
    # how well the model reproduces the entered measurements at the best estimate
    names = plots._PLAIN_KO if lang == "ko" else {}
    fit_rows = " · ".join(
        _t("inv.fit_row", lang, name=names.get(t_, t_),
           meas=est["measured"][t_], model=est["model_fit"][t_])
        for t_ in est["used_tissues"])
    st.caption(_t("inv.fit_caption", lang, rows=fit_rows))


_ORGANS4 = ("root", "stem", "leaf", "grain")


def _default_growth_df(season, biomass, n_rows=7):
    """Editable growth table seeded from the selected biomass driver (FRESH g/hill)."""
    import pandas as pd
    t = np.linspace(0.0, float(season), n_rows)
    b = api._biomass_fn(biomass)(t, float(season))
    return pd.DataFrame({"day": np.round(t, 0),
                         **{o: np.round(np.asarray(b[o], float) * 1e3, 2) for o in _ORGANS4}})


def _default_cwo_df(season, level, n_rows=4):
    import pandas as pd
    return pd.DataFrame({"day": np.round(np.linspace(0.0, float(season), n_rows), 0),
                         "Cwo": np.full(n_rows, float(level))})


def _clean_table(df, value_cols):
    """data_editor DataFrame -> {col: array}, dropping incomplete rows, day required."""
    import pandas as pd
    if df is None or "day" not in getattr(df, "columns", []):
        raise ValueError("the table needs a 'day' column")
    d = df.dropna(subset=["day"])
    out = {"day": d["day"].to_numpy(float)}
    mask = np.isfinite(out["day"])
    for c in value_cols:
        if c in d.columns:
            out[c] = pd.to_numeric(d[c], errors="coerce").to_numpy(float)
            mask &= np.isfinite(out[c])
    out = {k: v[mask] for k, v in out.items()}
    if len(out["day"]) < 2:
        raise ValueError("need at least 2 complete rows (day + values)")
    return out


def _render_custom_tables(*, biomass, Cwo_const, season0, key, ko=False):
    """Editable growth + pore-water tables (+ per-compartment density). Returns
    (drivers_dict_or_None, density_dict). Used by both Simple (Korean) and Expert."""
    import pandas as pd
    lang = "ko" if ko else "en"
    organ = i18n.ORGAN_LABELS[lang]
    st.markdown(_t("ct.intro", lang))
    c1, c2 = st.columns(2)
    with c1:
        st.caption(_t("ct.growth_caption", lang))
        up_g = st.file_uploader(_t("ct.growth_upload", lang), type=["csv"], key=f"{key}_gup")
        g_seed = pd.read_csv(up_g) if up_g is not None else _default_growth_df(season0, biomass)
        gdf = st.data_editor(g_seed, num_rows="dynamic", width="stretch",
                             key=f"{key}_growth_{up_g.name if up_g else 'def'}")
        gunit = st.selectbox(_t("ct.growth_units", lang),
                             ["g/hill", "kg/hill", "g/m2", "kg/ha", "t/ha"],
                             index=0, key=f"{key}_gunit")
    with c2:
        st.caption(_t("ct.cwo_caption", lang))
        up_c = st.file_uploader(_t("ct.cwo_upload", lang), type=["csv"], key=f"{key}_cup")
        c_seed = pd.read_csv(up_c) if up_c is not None else _default_cwo_df(season0, Cwo_const)
        cdf = st.data_editor(c_seed, num_rows="dynamic", width="stretch",
                             key=f"{key}_cwo_{up_c.name if up_c else 'def'}")
    st.markdown(_t("ct.density_md", lang))
    dc = st.columns(4)
    density = {o: dc[i].number_input(f"ρ {organ[o]}", 0.05, 2.0,
                                     float(api.DEFAULT_TISSUE_DENSITY[o]), 0.05, key=f"{key}_rho_{o}")
               for i, o in enumerate(_ORGANS4)}
    try:
        growth = _clean_table(gdf, list(_ORGANS4))
        cwo = _clean_table(cdf, ["Cwo"])
        drivers = api.drivers_from_tables(growth, cwo, growth_units=gunit,
                                          Cwo_const=Cwo_const, biomass=biomass)
    except Exception as e:                                   # noqa: BLE001
        st.warning(_t("ct.read_error", lang, e=e))
        return None, density
    Mf = np.asarray(drivers["M"], float)[-1]
    vols = {o: Mf[i] / max(density[o], 1e-6) for i, o in enumerate(_ORGANS4)}
    rows = " · ".join(f"{organ[o]} {vols[o] * 1e3:.0f} mL" for o in _ORGANS4)
    st.caption(_t("ct.implied_volume", lang, rows=rows))
    return drivers, density


# A small, theme-agnostic CSS polish on top of the config.toml design tokens. Uses
# NEUTRAL rgba overlays (not hardcoded light/dark colours) so it reads correctly in
# BOTH the light and dark themes; config.toml carries the palette/font/radius.
_APP_CSS = """
<style>
:root{
  --pfas-safe:#0E6B4F; --pfas-safe-bg:#DCEFE6; --pfas-safe-bd:#A9D6C3;
  --pfas-warn:#9A5A00; --pfas-warn-bg:#FAEBD1; --pfas-warn-bd:#E7C27F;
  --pfas-dang:#B23A2E; --pfas-dang-bg:#FADEDA; --pfas-dang-bd:#EBA99F;
  --pfas-accent:#0E7A63; --pfas-border:#E4DCCE; --pfas-surface:#FFFFFF;
  --pfas-ink:#211E18; --pfas-sub:#6B6456;
  --pfas-shadow:0 1px 2px rgba(40,34,24,.05), 0 6px 20px rgba(40,34,24,.07);
  --pfas-shadow-sm:0 1px 2px rgba(40,34,24,.06);
}
/* Dark tokens key off the ACTUAL Streamlit app theme (stamped on <html> by the
   JS probe in _sync_theme), NOT the OS `prefers-color-scheme` — otherwise a user
   who forces the light theme while their OS is dark gets dark cards/shadows on a
   light app. Default (no attribute yet / probe failed) stays light. */
:root[data-pfas-theme="dark"]{
    --pfas-safe:#3BCB9C; --pfas-safe-bg:#16302A; --pfas-safe-bd:#2E6152;
    --pfas-warn:#E8B24C; --pfas-warn-bg:#322813; --pfas-warn-bd:#6B5726;
    --pfas-dang:#F0796E; --pfas-dang-bg:#351F1C; --pfas-dang-bd:#7A413A;
    --pfas-accent:#35C79E; --pfas-border:#332F26; --pfas-surface:#201D16;
    --pfas-ink:#ECE6D9; --pfas-sub:#B7AE9C;
    --pfas-shadow:0 1px 2px rgba(0,0,0,.4), 0 8px 24px rgba(0,0,0,.35);
    --pfas-shadow-sm:0 1px 2px rgba(0,0,0,.4);
}
/* Toss-style: soft-shadow rounded cards, minimal borders, generous spacing. */
.block-container{ padding-top:2.4rem; padding-bottom:3.2rem; max-width:1120px; }
h1,h2,h3{ letter-spacing:-.02em; }
/* safety-signal badge: colour + SHAPE/ICON + label (colour-blind 4-way encoding) */
.pfas-badge{ display:inline-flex; gap:6px; align-items:center; font-size:13px;
  font-weight:700; padding:5px 13px; border-radius:999px; }
.pfas-badge.safe{ color:var(--pfas-safe); background:var(--pfas-safe-bg); }
.pfas-badge.warn{ color:var(--pfas-warn); background:var(--pfas-warn-bg); }
.pfas-badge.dang{ color:var(--pfas-dang); background:var(--pfas-dang-bg); }
/* caveat pill (uncertainty reminder) */
.pfas-caveat{ display:inline-flex; gap:6px; align-items:center; font-size:12px;
  padding:5px 12px; border-radius:999px; background:color-mix(in oklab, var(--pfas-ink) 5%, transparent);
  color:var(--pfas-sub); }
/* metric cards: floating white surface, soft shadow, strong value */
[data-testid="stMetric"]{ background:var(--pfas-surface); border:1px solid var(--pfas-border);
  border-radius:18px; padding:18px 20px; box-shadow:var(--pfas-shadow);
  transition:transform .12s ease, box-shadow .12s ease; }
[data-testid="stMetric"]:hover{ transform:translateY(-2px);
  box-shadow:0 2px 4px rgba(17,24,39,.05), 0 12px 30px rgba(17,24,39,.09); }
[data-testid="stMetricValue"]{ font-weight:800; letter-spacing:-.03em; }
[data-testid="stMetric"] [data-testid="stMetricLabel"]{ color:var(--pfas-sub); font-weight:600; }
/* tabs: pill-style, active tab filled with the blue accent */
.stTabs [data-baseweb="tab-list"]{ gap:6px; flex-wrap:wrap; border-bottom:none; }
.stTabs [data-baseweb="tab"]{ font-weight:600; padding:9px 24px; border-radius:999px;
  color:var(--pfas-sub); background:color-mix(in oklab, var(--pfas-ink) 4%, transparent); }
.stTabs [aria-selected="true"]{ color:#fff; font-weight:700;
  background:var(--pfas-accent); }
/* primary buttons: rounded, bold, blue */
.stButton>button, .stDownloadButton>button{ border-radius:12px; font-weight:700;
  border:1px solid var(--pfas-border); box-shadow:var(--pfas-shadow-sm); }
.stButton>button[kind="primary"]{ border:none;
  box-shadow:0 4px 14px color-mix(in oklab, var(--pfas-accent) 36%, transparent); }
/* banners + expanders: rounded, shadowed, borderless */
[data-testid="stAlert"]{ border-radius:14px; border:none; box-shadow:var(--pfas-shadow-sm); }
[data-testid="stExpander"]{ border:1px solid var(--pfas-border); border-radius:16px;
  box-shadow:var(--pfas-shadow-sm); overflow:hidden; }
/* inputs / sliders a touch rounder */
[data-baseweb="select"]>div, .stNumberInput input, .stTextInput input{ border-radius:12px; }
/* sidebar scenario radio -> tappable cards, selected one accented */
section[data-testid="stSidebar"] div[role="radiogroup"] label{
  border:1px solid var(--pfas-border); border-radius:14px; background:var(--pfas-surface);
  padding:12px 14px; margin-bottom:8px; width:100%; box-shadow:var(--pfas-shadow-sm);
  transition:border-color .12s ease; }
section[data-testid="stSidebar"] div[role="radiogroup"] label:has(input:checked){
  border:2px solid var(--pfas-accent);
  background:color-mix(in oklab, var(--pfas-accent) 8%, var(--pfas-surface)); }
/* ---- Simple-mode hero + slim disclaimer + result strip ---- */
.pfas-hero{ border:1px solid var(--pfas-border); border-radius:22px;
  padding:26px 28px 22px; margin:2px 0 14px; background:var(--pfas-surface);
  box-shadow:var(--pfas-shadow); }
.pfas-hero-kicker{ display:inline-block; font-size:12.5px; font-weight:700;
  color:var(--pfas-accent); background:color-mix(in oklab, var(--pfas-accent) 10%, transparent);
  padding:4px 11px; border-radius:999px; }
.pfas-hero-title{ font-size:32px; font-weight:800; letter-spacing:-.03em;
  margin:.42em 0 .3em; line-height:1.18; color:var(--pfas-ink); }
.pfas-hero-sub{ font-size:15.5px; line-height:1.6; color:var(--pfas-sub); max-width:62ch; }
.pfas-disc{ display:flex; gap:8px; align-items:flex-start; font-size:12.5px;
  line-height:1.5; color:var(--pfas-warn); background:var(--pfas-warn-bg);
  border-radius:12px; padding:9px 14px; margin:0 0 8px; }
/* the headline result strip: big signal + one-line takeaway on a soft card */
.pfas-result{ display:flex; gap:14px; align-items:center; flex-wrap:wrap;
  border-radius:18px; padding:16px 20px; margin:4px 0 8px; background:var(--pfas-surface);
  box-shadow:var(--pfas-shadow); }
.pfas-result.safe{ background:color-mix(in oklab, var(--pfas-safe-bg) 60%, var(--pfas-surface)); }
.pfas-result.warn{ background:color-mix(in oklab, var(--pfas-warn-bg) 60%, var(--pfas-surface)); }
.pfas-result.dang{ background:color-mix(in oklab, var(--pfas-dang-bg) 60%, var(--pfas-surface)); }
.pfas-result .pfas-badge{ font-size:15px; padding:8px 16px; }
.pfas-result-text{ font-size:16px; line-height:1.5; flex:1; min-width:240px; color:var(--pfas-ink); }
.pfas-result-text b{ font-weight:800; }
</style>
"""


# A 0-height helper iframe whose JS reads the REAL app background (Streamlit does
# not expose its theme as a CSS var or attribute) and stamps data-pfas-theme on the
# parent <html>, so the CSS above follows the app theme, not the OS. Same-origin
# srcdoc iframe → window.parent access is allowed. Re-runs each rerun.
_THEME_PROBE = """
<script>
(function(){
  function apply(){
    try{
      var doc = window.parent.document;
      var app = doc.querySelector('.stApp');
      if(!app) return;
      var m = getComputedStyle(app).backgroundColor.match(/\\d+/g);
      if(!m || m.length < 3) return;
      var lum = 0.2126*(+m[0]) + 0.7152*(+m[1]) + 0.0722*(+m[2]);
      doc.documentElement.setAttribute('data-pfas-theme', lum < 128 ? 'dark' : 'light');
    }catch(e){}
  }
  apply();
  var n = 0, id = setInterval(function(){ apply(); if(++n > 12) clearInterval(id); }, 200);
})();
</script>
"""


def _sync_theme():
    """Detect the actual Streamlit theme and mirror it onto <html data-pfas-theme>."""
    components.html(_THEME_PROBE, height=0)


def inject_css():
    """Inject the small CSS polish once per render (idempotent), then sync the theme."""
    st.markdown(_APP_CSS, unsafe_allow_html=True)
    _sync_theme()


# Safety signal: colour-blind-safe 4-way encoding (colour + shape + KO + EN label
# + fixed order safe->caution->exceed). `signal` is the model_api intake signal.
_SIGNAL = {
    "green": ("safe", "●", "안전", "Safe"),
    "amber": ("warn", "▲", "주의", "Caution"),
    "red":   ("dang", "⬢", "초과", "Exceeds"),
}


def signal_badge_html(signal, *, ko=True, extra=""):
    """An HTML pill badge for a safety signal: colour + shape + label (not colour
    alone), so it reads under colour-blindness. `extra` appends a short suffix
    inside the badge (e.g. a percentage). Returns '' for an unknown signal."""
    if signal not in _SIGNAL:
        return ""
    cls, shape, ko_l, en_l = _SIGNAL[signal]
    label = f"{ko_l} · {en_l}" if ko else en_l
    tail = f" {extra}" if extra else ""
    return f"<span class='pfas-badge {cls}'>{shape} {label}{tail}</span>"


# ---------------------------------------------------------------- render building blocks
def render_header(cfg):
    """Title + disclaimer + intro (both modes)."""
    expert = cfg.expert
    inject_css()
    if expert:
        st.title("🌾 PFAS in Rice — Uptake Explorer")
        st.warning(_DISCLAIMER)
        st.caption(_t("header.expert_caption", "en"))
        return
    # --- Simple (Korean): a compact hero, then a slim (not heavy) disclaimer ---
    st.markdown(
        "<div class='pfas-hero'>"
        "<div class='pfas-hero-kicker'>논 PFAS 벼 축적 시뮬레이터</div>"
        "<h1 class='pfas-hero-title'>🌾 벼의 어디에, 얼마나 쌓일까?</h1>"
        "<div class='pfas-hero-sub'>논의 물·흙에 녹은 <b>'영원한 화학물질' PFAS</b>가 "
        "벼의 <b>뿌리·짚·먹는 낟알</b>에 얼마나 쌓이는지 한눈에 보여줍니다.</div>"
        "</div>", unsafe_allow_html=True)
    st.markdown(
        "<div class='pfas-disc'>⚠️ 연구·교육용 예시 추정치입니다. 규제·식품안전·건강 판단이 "
        "<b>아니며</b>, 실제 노출·안전 결정에 사용하지 마세요.</div>", unsafe_allow_html=True)
    st.caption(_t("header.intro2", "ko"))


def render_custom_tables_panel(cfg):
    """Optional editable growth + Cwᵒ tables (both modes); may override cfg.drivers."""
    expert = cfg.expert
    use_custom_tables = cfg.use_custom_tables
    biomass = cfg.biomass
    Cwo_const = cfg.Cwo_const
    season = cfg.season
    drivers = cfg.drivers
    custom_density = None
    if use_custom_tables:
        _tbl_title = _t("ct.panel_title", "en" if expert else "ko")
        with st.expander(_tbl_title, expanded=True):
            _drv, custom_density = _render_custom_tables(biomass=biomass, Cwo_const=Cwo_const,
                                                         season0=season, key="ctbl", ko=not expert)
            if _drv is not None:
                drivers = _drv
    cfg.drivers = drivers
    cfg.custom_density = custom_density


def run_model(cfg):
    """Run the full ODE (curated congener or SMILES) and attach results to cfg."""
    congener = cfg.congener
    smiles = cfg.smiles
    drivers = cfg.drivers
    Cwo_const = cfg.Cwo_const
    season = cfg.season
    measured = cfg.measured
    cwo_profile = cfg.cwo_profile
    cwo_kleach = cfg.cwo_kleach
    E_m = cfg.E_m
    fxy_source = cfg.fxy_source
    biomass = cfg.biomass
    measured_bio = cfg.measured_bio
    sim_kw = dict(E_m_mV=E_m, f_xy_source=fxy_source, biomass=biomass)
    desc = None
    provisional = False
    if smiles:                                              # compound specified by structure
        if not api.rdkit_available():
            st.error("RDKit not installed — cannot parameterise a SMILES structure. "
                     "`pip install rdkit`, or switch to **Curated congener** in the sidebar.")
            st.stop()
        try:
            if drivers is not None:
                res = _simulate_smiles(smiles, drivers_tuple=_drivers_tuple(drivers), **sim_kw)
            else:
                res = _simulate_smiles(smiles, Cwo=Cwo_const, season=season,
                                       measured_forcing=measured, cwo_profile=cwo_profile,
                                       cwo_k_leach=cwo_kleach, **sim_kw)
        except Exception as e:                              # noqa: BLE001
            st.error(f"Could not build a compound from that SMILES — check the structure.\n\n`{e}`")
            st.stop()
        congener = res["congener"]
        desc = res.get("descriptors")
        provisional = bool(res.get("provisional", False))
    elif drivers is not None:
        with st.spinner("🌾 벼 한 철 축적을 계산하는 중…" if not cfg.expert else "Running the model…"):
            res = _simulate(congener, drivers_tuple=_drivers_tuple(drivers), **sim_kw)
    else:
        with st.spinner("🌾 벼 한 철 축적을 계산하는 중…" if not cfg.expert else "Running the model…"):
            res = _simulate(congener, Cwo=Cwo_const, season=season, measured_forcing=measured,
                            cwo_profile=cwo_profile, cwo_k_leach=cwo_kleach, **sim_kw)
    obs = api.observed_baf(congener)
    p = res["params"]

    # biomonitoring-derived BAFs (measured side)
    bio_baf = None
    if measured_bio and measured_bio.get("Cwo"):
        bio_baf = api.baf_from_measurement(measured_bio["conc"], measured_bio["Cwo"])
    cfg.sim_kw = sim_kw
    cfg.res = res
    cfg.desc = desc
    cfg.provisional = provisional
    cfg.congener = congener
    cfg.obs = obs
    cfg.p = p
    cfg.bio_baf = bio_baf


def render_footer(cfg):
    """Footer shown on every screen."""
    expert = cfg.expert
    lang = "en" if expert else "ko"
    st.divider()
    st.caption(_DISCLAIMER if expert else _DISCLAIMER_KO)
    fc1, fc2, fc3 = st.columns(3)
    fc1.caption(f"**PFAS–Rice Uptake Model** · v{APP_VERSION}")
    fc2.caption(_t("footer.links", lang, repo=REPO_URL, docs=DOCS_URL))
    fc3.caption(_t("footer.cite", lang))
