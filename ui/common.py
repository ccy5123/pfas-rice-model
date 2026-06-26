"""Shared UI layer for the PFAS-rice dashboard: constants, cached model helpers,
and the render building blocks that app.py assembles. Split out of the monolithic
app.py (HANDOFF P3-1); behaviour is unchanged."""
import os

import numpy as np
import streamlit as st

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


# ---------------------------------------------------------------- render building blocks
def render_header(cfg):
    """Title + disclaimer + intro (both modes)."""
    expert = cfg.expert
    st.title("🌾 PFAS in Rice — Uptake Explorer")
    st.warning(_DISCLAIMER if expert else _DISCLAIMER_KO)

    if not expert:
        st.markdown(_t("header.intro1", "ko"))
        st.markdown(_t("header.intro2", "ko"))
    else:
        st.caption(_t("header.expert_caption", "en"))


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
        res = _simulate(congener, drivers_tuple=_drivers_tuple(drivers), **sim_kw)
    else:
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
