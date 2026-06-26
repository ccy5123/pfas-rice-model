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
import os
import sys

import numpy as np
import streamlit as st

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))
import model_api as api          # noqa: E402
import plots                     # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
_EX = os.path.join(_HERE, "examples")

APP_VERSION = "1.0 (general-audience UI)"
REPO_URL = "https://github.com/ccy5123/pfas-rice-model"
DOCS_URL = "https://github.com/ccy5123/pfas-rice-model/tree/main/docs"

# Plain-language disclaimer reused at the top of the page and in the footer.
_DISCLAIMER = (
    "**Research & educational model — illustrative estimates only.** "
    "This is **not** a regulatory, food-safety, or health determination. "
    "Do **not** use it for real exposure or safety decisions."
)
# Korean disclaimer for the general-audience (Simple) mode.
_DISCLAIMER_KO = (
    "**연구·교육용 모델 — 예시 추정치일 뿐입니다.** "
    "규제·식품안전·건강 판단이 **아니며**, 실제 노출·안전 결정에 **사용하지 마세요**."
)

# Friendlier congener names for the Simple-mode dropdown (value stays the symbol).
_FRIENDLY_CONG = {
    "PFBA":   "PFBA — short-chain acid (C4)",
    "PFPeA":  "PFPeA — short-chain acid (C5)",
    "PFHxA":  "PFHxA — short-chain acid (C6)",
    "PFHpA":  "PFHpA — medium-chain acid (C7)",
    "PFOA":   "PFOA — common 'forever chemical' acid (C8)",
    "PFNA":   "PFNA — long-chain acid (C9)",
    "PFDA":   "PFDA — long-chain acid (C10)",
    "PFUnDA": "PFUnDA — long-chain acid (C11)",
    "PFDoDA": "PFDoDA — long-chain acid (C12)",
    "PFBS":   "PFBS — short-chain sulfonate (C4)",
    "PFHxS":  "PFHxS — sulfonate (C6)",
    "PFOS":   "PFOS — common 'forever chemical' sulfonate (C8)",
    "GenX":   "GenX — newer PFOA replacement (ether)",
}


def _cong_label(name):
    return _FRIENDLY_CONG.get(name, name)


# Korean congener labels for the Simple-mode dropdown.
_FRIENDLY_CONG_KO = {
    "PFBA":   "PFBA — 단쇄 카복실산 (C4)",
    "PFPeA":  "PFPeA — 단쇄 카복실산 (C5)",
    "PFHxA":  "PFHxA — 단쇄 카복실산 (C6)",
    "PFHpA":  "PFHpA — 중쇄 카복실산 (C7)",
    "PFOA":   "PFOA — 대표적 '영원한 화학물질' 카복실산 (C8)",
    "PFNA":   "PFNA — 장쇄 카복실산 (C9)",
    "PFDA":   "PFDA — 장쇄 카복실산 (C10)",
    "PFUnDA": "PFUnDA — 장쇄 카복실산 (C11)",
    "PFDoDA": "PFDoDA — 장쇄 카복실산 (C12)",
    "PFBS":   "PFBS — 단쇄 술폰산 (C4)",
    "PFHxS":  "PFHxS — 술폰산 (C6)",
    "PFOS":   "PFOS — 대표적 '영원한 화학물질' 술폰산 (C8)",
    "GenX":   "GenX — PFOA 대체물질 (에터)",
}


def _cong_label_ko(name):
    return _FRIENDLY_CONG_KO.get(name, name)


# Plain "how contaminated?" presets → pore-water concentration [µg/L].
# Medium = 1.0 µg/L keeps tissue conc == build-up factor (the model's reference point).
_PRESETS = {
    "Low — lightly contaminated (0.1 µg/L)": 0.1,
    "Medium — moderately contaminated (1 µg/L)": 1.0,
    "High — heavily contaminated (10 µg/L)": 10.0,
}
# Korean presets for Simple mode (label -> µg/L). The short word (저/중/고) is reused
# in the headline sentence.
_PRESETS_KO = {
    "낮음 — 약하게 오염 (0.1 µg/L)": (0.1, "낮은"),
    "중간 — 보통 오염 (1 µg/L)": (1.0, "중간"),
    "높음 — 심하게 오염 (10 µg/L)": (10.0, "높은"),
}

st.set_page_config(page_title="PFAS–Rice Uptake Model", layout="wide")


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


@st.cache_data(show_spinner="Estimating the contamination level…")
def _estimate_exposure(congener, root, straw, grain, sigma, E_m_mV, f_xy_source, biomass):
    """Cache the Bayesian inverse estimate (a handful of ODE solves, a few seconds)."""
    meas = {k: v for k, v in (("root", root), ("straw", straw), ("grain", grain))
            if v is not None and v > 0}
    return api.estimate_exposure_bayesian(congener, meas, sigma_log10=sigma, E_m_mV=E_m_mV,
                                          f_xy_source=f_xy_source, biomass=biomass)


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


def _glossary_md(ko=False):
    """Plain-language glossary (rendered in the About tab and a Simple-mode expander)."""
    if ko:
        return (
            "| 보게 될 용어 | 쉬운 설명 |\n"
            "|---|---|\n"
            "| **PFAS** | 잘 분해되지 않는 인공 '영원한 화학물질' 무리. |\n"
            "| **토양수 농도** *(Cwᵒ)* | 뿌리 주변 토양수에 녹아 있는 PFAS의 양. |\n"
            "| **축적 배수** *(BAF)* | 식물 조직이 토양수보다 PFAS를 몇 배나 진하게 모았는지. 2배면 토양수의 두 배. |\n"
            "| **뿌리 / 짚 / 낟알** | 식물 부위. *짚* = 줄기 + 잎, *낟알* = 먹는 현미. |\n"
            "| **농도** *(µg/kg)* | 식물 조직 1 kg당 PFAS 마이크로그램. |\n"
            "| **화학종(congener)** | 특정 PFAS 하나(예: PFOA, PFOS). 탄소 사슬이 길수록 대체로 식물에 더 잘 달라붙음. |\n"
            "| **흡수 / 이동** | 화학물질이 뿌리로 들어가 줄기·낟알로 올라가는 과정. |\n"
            "| **베이지안 추정** | 측정값에서 거꾸로 가장 가능성 높은 원인을 찾되, 하나의 숫자가 아니라 **불확실성 범위**까지 제시. |\n")
    return (
        "| Term you'll see | What it means in plain words |\n"
        "|---|---|\n"
        "| **PFAS** | A family of long-lasting synthetic 'forever chemicals'. |\n"
        "| **Pore-water level** *(Cwᵒ)* | How much PFAS is dissolved in the soil water around the roots. |\n"
        "| **Build-up factor** *(BAF)* | How many times more concentrated the PFAS is in the plant tissue than in the soil water. A factor of 2 means the tissue holds twice the water's level. |\n"
        "| **Roots / Straw / Grain** | The plant parts. *Straw* = the stems + leaves together. *Grain* = the edible brown rice. |\n"
        "| **Concentration** *(µg/kg)* | Micrograms of PFAS per kilogram of plant tissue. |\n"
        "| **Congener** | One specific PFAS chemical (e.g. PFOA, PFOS). Longer carbon chains generally stick to the plant more. |\n"
        "| **Uptake / translocation** | How the chemical gets into the roots and then moves up into the shoot and grain. |\n"
        "| **Bayesian estimate** | Working backwards from a measurement to the most likely cause, **with an uncertainty range** instead of a single number. |\n")


# uncertainty presets for the inverse estimator (measurement + model noise, log10 units)
_UNC = {"Typical (±~40%)": 0.15, "High precision (±~20%)": 0.10, "Rough (±~2×)": 0.30}


def _render_inverse_estimator(congener, *, E_m_mV, f_xy_source, biomass, key, simple=True):
    """Shared 'work backwards' panel: Bayesian estimate of the soil-water contamination
    level Cwᵒ from measured tissue concentrations, with a credible interval. Used by
    both the Simple (Korean) and Expert (English) tabs. `key` namespaces the widgets."""
    ko = simple
    if congener not in api.CONGENERS:
        st.info("선별된 13종 화학물질에서만 작동합니다 — 사이드바에서 하나를 고르세요 "
                "(보정된 모델이 필요하며, 직접 입력한 SMILES 구조에는 없습니다)." if ko else
                "This works with the curated chemicals — pick one of them in the sidebar "
                "(it needs the calibrated model, which a custom SMILES structure doesn't have).")
        return
    st.markdown(
        "오염된 땅에서 자란 벼의 **실험실 측정값**이 있으신가요? 식물에서 측정된 값을 입력하면 "
        "모델을 거꾸로 돌려 **토양수가 얼마나 오염됐을지**를 — 하나의 숫자가 아니라 "
        "**불확실성 범위**까지(베이지안 추정) — 추정합니다." if ko else
        "Already have a **lab result** for rice grown on contaminated land? Enter what was "
        "measured in the plant and this estimates **how contaminated the soil water likely "
        "was** — working the model backwards, with an **uncertainty range** (a Bayesian "
        "estimate), not just a single number.")
    c1, c2, c3 = st.columns(3)
    root = c1.number_input("뿌리 측정값 [µg/kg]" if ko else "Measured in roots [µg/kg]",
                           0.0, 1e6, 0.0, 0.1, key=f"{key}_root")
    straw = c2.number_input("짚(줄기+잎) 측정값 [µg/kg]" if ko else "Measured in straw (stems+leaves) [µg/kg]",
                            0.0, 1e6, 0.0, 0.1, key=f"{key}_straw")
    grain = c3.number_input("낟알 측정값 [µg/kg]" if ko else "Measured in grain [µg/kg]",
                            0.0, 1e6, 0.0, 0.1, key=f"{key}_grain")
    unc_label = st.radio("측정값이 얼마나 정밀한가요?" if ko else "How precise are the measurements?",
                         list(_UNC), horizontal=True, key=f"{key}_unc",
                         help="추정에 쓰이는 측정+모델 불확실성을 설정합니다." if ko else
                              "Sets the measurement+model uncertainty used in the estimate.")
    sigma = _UNC[unc_label]
    have = any(v > 0 for v in (root, straw, grain))
    run = st.button("📐 오염 수준 추정하기" if ko else "📐 Estimate the contamination level",
                    key=f"{key}_btn", disabled=not have, type="primary")
    sig = (congener, root, straw, grain, sigma, E_m_mV, f_xy_source, biomass)
    if run:
        st.session_state[f"{key}_sig"] = sig
    if not have:
        st.caption("위에 측정값을 하나 이상 입력한 뒤 버튼을 누르세요." if ko else
                   "Enter at least one measured tissue concentration above, then press the button.")
        return
    if st.session_state.get(f"{key}_sig") != sig:
        st.caption("**추정** 버튼을 누르세요 (값을 바꾼 뒤에는 다시 누르세요)." if ko else
                   "Press **Estimate** to run (or re-run after changing a value).")
        return
    try:
        est = _estimate_exposure(congener, root, straw, grain, sigma, E_m_mV, f_xy_source, biomass)
    except Exception as e:                                   # noqa: BLE001
        st.error((f"추정을 실행할 수 없습니다: {e}") if ko else f"Could not run the estimate: {e}")
        return
    med = est["median"]
    lo, hi = est["ci95"]
    mc1, mc2 = st.columns([1, 2])
    mc1.metric("가장 가능성 높은 토양수 수준" if ko else "Most likely soil-water level", f"{med:.3g} µg/L",
               (f"95% 범위 {lo:.2g}–{hi:.2g}" if np.isfinite(lo) else "범위 미확정") if ko else
               (f"95% range {lo:.2g}–{hi:.2g}" if np.isfinite(lo) else "range unconstrained"),
               delta_color="off")
    if ko:
        mc2.markdown(
            f"입력하신 **{congener}** 측정값으로 보면, 토양수에 녹아 있던 PFAS는 가장 가능성 높게 "
            f"**{med:.3g} µg/L**"
            + (f"이며, **{lo:.2g}~{hi:.2g} µg/L** 사이일 확률이 95%입니다."
               if np.isfinite(lo) else "입니다 (측정값으로는 뚜렷한 범위가 좁혀지지 않았습니다).")
            + " 아래 곡선의 퍼짐이 불확실성입니다.")
    else:
        mc2.markdown(
            f"Given your measurements of **{congener}**, the model estimates the PFAS dissolved "
            f"in the soil water was most likely **{med:.3g} µg/L**"
            + (f", and we're 95% confident it was between **{lo:.2g}** and **{hi:.2g} µg/L**."
               if np.isfinite(lo) else " (the measurements didn't pin down a clear range).")
            + " The spread of the curve below is the uncertainty.")
    st.plotly_chart(plots.fig_exposure_posterior(est, lang="ko" if ko else "en"), width="stretch")
    # how well the model reproduces the entered measurements at the best estimate
    if ko:
        fit_rows = " · ".join(f"{plots._PLAIN_KO.get(t_, t_)}: 입력 {est['measured'][t_]:.3g} vs 모델 {est['model_fit'][t_]:.3g}"
                              for t_ in est["used_tissues"])
        st.caption(f"최적 추정에서 모델이 입력값을 재현합니다 — {fit_rows} (µg/kg). "
                   "이는 **오염 수준**만 추정하며, 모델의 흡수 거동이 옳다고 가정합니다 — 특정 논의 실측이 아닌 예시입니다.")
    else:
        fit_rows = " · ".join(f"{t_}: you {est['measured'][t_]:.3g} vs model {est['model_fit'][t_]:.3g}"
                              for t_ in est["used_tissues"])
        st.caption(f"At the best estimate the model reproduces your inputs — {fit_rows} (µg/kg). "
                   "This estimates only the **contamination level**; it assumes the model's "
                   "plant-uptake is right and can't separately tell apart water uptake vs how "
                   "the chemical moves inside the plant (that needs a sap/soil-water measurement).")


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
    st.markdown(
        "**나만의 한 철**을 두 표로 입력하세요 — **성장 = 기관별 신선중**, **공극수 = 토양수에 "
        "녹아 있는 PFAS 절대 농도(µg/L)**의 시간 변화. 셀을 편집·행 추가/삭제하거나 스프레드시트에서 "
        "붙여넣으세요; 날짜는 모델 타임라인으로 보간됩니다. 어느 한 표를 기본값으로 두면 내장값을 씁니다."
        if ko else
        "Enter your **own season** as two tables — **growth = organ FRESH weight** and "
        "**pore water = the absolute PFAS dissolved in the soil water (µg/L)** over time. "
        "Edit cells, add/remove rows, or paste from a spreadsheet; the days are interpolated "
        "onto the model timeline. Leave either table at its default to use the built-in value.")
    c1, c2 = st.columns(2)
    with c1:
        st.caption("🌱 성장 — 포기당 기관별 신선중" if ko else "🌱 Growth — organ FRESH weight per hill")
        up_g = st.file_uploader("…또는 성장 CSV 업로드 (day,root,stem,leaf,grain)" if ko else
                                "…or upload a growth CSV (day,root,stem,leaf,grain)",
                                type=["csv"], key=f"{key}_gup")
        g_seed = pd.read_csv(up_g) if up_g is not None else _default_growth_df(season0, biomass)
        gdf = st.data_editor(g_seed, num_rows="dynamic", width="stretch",
                             key=f"{key}_growth_{up_g.name if up_g else 'def'}")
        gunit = st.selectbox("성장 무게 단위" if ko else "Growth weight units",
                             ["g/hill", "kg/hill", "g/m2", "kg/ha", "t/ha"],
                             index=0, key=f"{key}_gunit")
    with c2:
        st.caption("💧 공극수 오염 Cwᵒ(t) — 절대 µg/L" if ko else
                   "💧 Pore-water contamination Cwᵒ(t) — absolute µg/L")
        up_c = st.file_uploader("…또는 Cwᵒ CSV 업로드 (day,Cwo)" if ko else
                                "…or upload a Cwᵒ CSV (day,Cwo)", type=["csv"], key=f"{key}_cup")
        c_seed = pd.read_csv(up_c) if up_c is not None else _default_cwo_df(season0, Cwo_const)
        cdf = st.data_editor(c_seed, num_rows="dynamic", width="stretch",
                             key=f"{key}_cwo_{up_c.name if up_c else 'def'}")
    st.markdown("**구획 밀도** ρ [kg/L, 신선] — 입력한 무게를 조직 부피와 연결합니다 "
                "(벼 잎/줄기는 통기조직으로 < 1, 낟알은 더 조밀해 > 1)." if ko else
                "**Compartment density** ρ [kg/L, fresh] — links the entered weight to tissue "
                "volume (rice leaf/culm hold air spaces ⇒ < 1; grain is denser ⇒ > 1).")
    dc = st.columns(4)
    _rho_lbl = {"root": "뿌리", "stem": "줄기", "leaf": "잎", "grain": "낟알"}
    density = {o: dc[i].number_input(f"ρ {_rho_lbl[o] if ko else o}", 0.05, 2.0,
                                     float(api.DEFAULT_TISSUE_DENSITY[o]), 0.05, key=f"{key}_rho_{o}")
               for i, o in enumerate(_ORGANS4)}
    try:
        growth = _clean_table(gdf, list(_ORGANS4))
        cwo = _clean_table(cdf, ["Cwo"])
        drivers = api.drivers_from_tables(growth, cwo, growth_units=gunit,
                                          Cwo_const=Cwo_const, biomass=biomass)
    except Exception as e:                                   # noqa: BLE001
        st.warning((f"표를 읽지 못해 기본 시나리오를 사용합니다: {e}") if ko else
                   f"Using the default scenario — couldn't read the tables: {e}")
        return None, density
    Mf = np.asarray(drivers["M"], float)[-1]
    vols = {o: Mf[i] / max(density[o], 1e-6) for i, o in enumerate(_ORGANS4)}
    if ko:
        st.caption("추정된 수확기 기관 **부피** (신선중 ÷ 밀도): "
                   + " · ".join(f"{_rho_lbl[o]} {vols[o] * 1e3:.0f} mL" for o in _ORGANS4)
                   + ". 수송 모델은 신선중(질량)으로 적분하며, 밀도는 질량↔부피 환산용입니다.")
    else:
        st.caption("Implied end-of-season organ **volume** (fresh mass ÷ density): "
                   + " · ".join(f"{o} {vols[o] * 1e3:.0f} mL" for o in _ORGANS4)
                   + ". The transport model integrates on fresh mass; density sets the mass↔volume scale.")
    return drivers, density


# ---------------------------------------------------------------- sidebar
with st.sidebar:
    expert = st.toggle(
        "🔬 전문가/고급 모드 (Expert / advanced)", value=False,
        help="끄면(기본): 쉬운 한국어 화면 (화학물질 + 오염 수준만 선택). "
             "켜면: 전체 연구용 인터페이스(영어) — 5가지 노출 모드, SMILES 구조 입력, 모든 모델 파라미터.")

    # ---- shared scenario defaults (Simple mode uses these as-is) ----
    drivers = None
    measured = True
    Cwo_const = 1.0
    season = 120.0
    soil_obj = None
    profile = None
    measured_bio = None
    cwo_profile = "constant"
    cwo_kleach = 0.02
    smiles = None
    spec = "Curated congener"
    E_m = -120
    fxy_source = "recommended"
    biomass = "oryza"
    compare = []
    preset_label = None
    preset_word = None
    use_custom_tables = False

    if not expert:
        # ----------------------------- SIMPLE sidebar (한국어) -----------------------------
        mode = "Model (parametric)"
        st.header("① 화학물질 선택")
        congener = st.selectbox("PFAS 화학물질", api.CONGENERS,
                                index=api.CONGENERS.index("PFOA"),
                                format_func=_cong_label_ko,
                                help="특정 '영원한 화학물질' 하나. PFOA·PFOS가 가장 잘 알려져 있고, "
                                     "사슬이 길수록 대체로 식물에 더 잘 달라붙습니다.")

        st.header("② 오염 정도")
        preset_label = st.radio("논의 오염 수준", list(_PRESETS_KO), index=1,
                                help="토양수에 녹아 있는 PFAS의 양. 높을수록 식물로 더 많이 들어갑니다.")
        Cwo_const, preset_word = _PRESETS_KO[preset_label]
        use_custom_tables = st.checkbox(
            "📋 내 데이터 표 사용", value=False,
            help="성장 곡선과 시간에 따른 토양수 오염 수준을 편집 가능한 표로 직접 입력합니다 "
                 "(본문에 나타납니다).")
        st.caption("전체 연구용 인터페이스가 필요하면 위의 **전문가/고급 모드**를 켜세요.")
    else:
        # ----------------------------- EXPERT sidebar -----------------------------
        st.header("1 · Data source")
        mode = st.radio(
            "How is the soil exposure supplied?",
            ["Model (parametric)", "Custom tables (Cwᵒ + growth)", "HYDRUS / CSV drivers",
             "Run HYDRUS-1D (live)", "Soil inventory → pore water",
             "Biomonitoring (measured tissue)"],
            help="Ways to feed the plant model. 'Custom tables' lets you type/paste your own "
                 "growth curve + time-varying Cwᵒ; 'Run HYDRUS-1D (live)' executes the real "
                 "engine (if built); biomonitoring needs no soil model.")
        use_custom_tables = (mode == "Custom tables (Cwᵒ + growth)")

        st.header("2 · PFAS compound")
        spec = st.radio("Specify by", ["Curated congener", "SMILES (structure)"], horizontal=True,
                        help="Pick one of the 13 curated congeners, or paste ANY PFAS structure (SMILES) "
                             "to parameterise it from chemistry (RDKit read-across / QSPR).")
        if spec == "Curated congener":
            congener = st.selectbox("PFAS congener", api.CONGENERS,
                                    index=api.CONGENERS.index("PFOA"), format_func=_cong_label)
        else:
            congener = None
            _EXSMI = {
                "PFOA  (known)": "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
                "PFOS  (known)": "OS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
                "GenX / HFPO-DA  (known ether)": "OC(=O)C(F)(OC(F)(F)C(F)(F)C(F)(F)F)C(F)(F)F",
                "PFTrDA — NOVEL C13 PFCA": "OC(=O)" + "C(F)(F)" * 11 + "C(F)(F)F",
                "ADONA-like — NOVEL ether-PFCA": "OC(=O)C(F)(F)OC(F)(F)C(F)(F)OC(F)(F)C(F)(F)F",
            }
            ex = st.selectbox("Example structure", list(_EXSMI), index=0,
                              help="Known → reads the curated parameters; NOVEL → QSPR (provisional).")
            smiles = st.text_area("SMILES", _EXSMI[ex], height=70).strip()
            if not api.rdkit_available():
                st.warning("RDKit is not installed — the SMILES mode needs it.\n\n"
                           "`pip install rdkit`  (or `-r requirements-structure.txt`). "
                           "Meanwhile use **Curated congener**.")
            elif smiles:                                        # show the 2-D structure
                svg, why = _mol_svg(smiles)
                if svg is not None:
                    import streamlit.components.v1 as components
                    components.html(
                        "<div style='background:#fff;border:1px solid #ddd;border-radius:6px;"
                        f"display:flex;justify-content:center'>{svg}</div>", height=185)
                    st.caption("structure (RDKit)")
                else:
                    st.caption(f"⚠ {why}")
        E_m = st.slider("Root membrane potential E_m  [mV]", -160, -90, -120, 5,
                        help="GHK anion-exclusion lever (rice −116…−140 mV; NH₄⁺ depolarises).")
        fxy_label = st.radio("Root→shoot loading f_xy",
                             ["recommended (monotone, physical)", "W2 fit (reproduces Yamazaki)"])
        fxy_source = "recommended" if fxy_label.startswith("recommended") else "W2fit"
        bm_label = st.radio("Biomass driver M(t)",
                            ["ORYZA2000 (mechanistic)", "growth_rice (partition + logistic)"],
                            help="ORYZA2000 = the Level-1 carbon balance (radiation/temperature-driven; "
                                 "first-principles). growth_rice = ORYZA IR72 partitioning imposed on a "
                                 "logistic (lightweight; the historical calibration basis). Drives M(t) "
                                 "when the scenario uses built-in forcings (ignored if a driver CSV supplies M).")
        biomass = "oryza" if bm_label.startswith("ORYZA2000") else "growth_rice"

        st.header("3 · Scenario")
        if mode == "Model (parametric)":
            Cwo_const = st.number_input("Pore-water Cwᵒ  [µg/L]", min_value=0.0, value=1.0, step=0.1,
                                        help="Free anion conc. driving root uptake. 1.0 → tissue conc equals BAF.")
            season = float(st.slider("Season length  [days]", 90, 160, 120, 5))
            cwo_label = st.radio("Pore-water Cwᵒ(t) shape",
                                 ["constant (flat)", "flooded (dilution + leaching)"],
                                 help="constant → Cwᵒ held flat (tissue conc == BAF). flooded → analytic "
                                      "Freundlich paddy shape: short chains LEACH (decline), long chains stay "
                                      "buffered (~flat); the season-MEAN is held at Cwᵒ, so only the time shape "
                                      "changes. No HYDRUS engine needed (for the real engine use the 'Run "
                                      "HYDRUS-1D (live)' data source).")
            cwo_profile = "constant" if cwo_label.startswith("constant") else "flooded"
            if cwo_profile == "flooded":
                cwo_kleach = st.slider("Leaching rate k_leach  [1/day]", 0.0, 0.15,
                                       float(api.default_k_leach(congener)), 0.0025,
                                       help="Default is CALIBRATED per congener to a HYDRUS-1D run "
                                            "(short chains leach fast, long chains stay buffered ≈0). "
                                            "Higher → faster short-chain pore-water decline.")
                try:                                            # immediate shape feedback
                    st.plotly_chart(plots.fig_cwo_profile(congener, level=Cwo_const,
                                                          profile="flooded", season=season,
                                                          k_leach=cwo_kleach),
                                    width="stretch", theme=None)
                except Exception:                               # noqa: BLE001 (preview is non-essential)
                    pass
            measured = st.checkbox("Measured forcings (Q_TP, M_s)", value=True,
                                   help="On: transpiration from Kumari/NayHtoon, biomass M(t) from the "
                                        "sidebar biomass driver (ORYZA2000 mechanistic, or growth_rice).")

        elif mode == "Custom tables (Cwᵒ + growth)":
            st.caption("Enter the **growth** and **pore-water Cwᵒ(t)** tables in the main panel → "
                       "(editable grids + per-compartment density). Q_TP defaults to the measured "
                       "transpiration; omit either table to fall back to the built-in value.")

        elif mode == "HYDRUS / CSV drivers":
            st.caption("CSV columns: `t, Cwo, Qtp, M_root, M_stem, M_leaf, M_grain` "
                       "(Qtp/M optional → measured forcings). See the **About** tab for the HYDRUS map.")
            up = st.file_uploader("Driver CSV (HYDRUS-1D / Phydrus output)", type=["csv"])
            use_ex = st.checkbox("Use bundled example", value=up is None)
            try:
                if up is not None:
                    import pandas as pd
                    df = pd.read_csv(up)
                    cols = {c.lower(): c for c in df.columns}
                    t = df[cols["t"]].to_numpy(float)
                    Cwo = df[cols["cwo"]].to_numpy(float)
                    Q = df[cols["qtp"]].to_numpy(float) if "qtp" in cols else None
                    M = (df[[cols["m_root"], cols["m_stem"], cols["m_leaf"], cols["m_grain"]]].to_numpy(float)
                         if all(k in cols for k in ("m_root", "m_stem", "m_leaf", "m_grain")) else None)
                    drivers = api.drivers_from_arrays(t, Cwo, Qtp=Q, M=M, biomass=biomass)
                    st.success(f"Loaded {len(t)} rows from {up.name}.")
                elif use_ex:
                    drivers = api.load_driver_csv(os.path.join(_EX, "hydrus_drivers_example.csv"),
                                                  biomass=biomass)
                    st.info("Using examples/hydrus_drivers_example.csv (synthetic HYDRUS-style run).")
            except Exception as e:                                  # noqa: BLE001
                st.error(f"Could not read drivers: {e}")

        elif mode == "Run HYDRUS-1D (live)":
            # SMILES compounds have no congener name -> use a curated congener for the soil
            # Kd (read-across match, else nearest by chain length). The PLANT run still uses
            # the actual SMILES compound (the drivers are just Cwᵒ(t)/Q_TP).
            soil_cong, how = (congener, None)
            if smiles:
                soil_cong, how = _hydrus_soil_congener(smiles)
            if not api.hydrus_available():
                st.warning("HYDRUS-1D engine not built in this environment yet — build it once "
                           "(compiles the FORTRAN solver with gfortran; ~1 min, cached).")
                if st.button("⚙ Build the HYDRUS-1D engine now"):
                    with st.spinner("Fetching source + compiling HYDRUS-1D with gfortran…"):
                        ok, blog = api.build_hydrus_engine()
                    if ok:
                        st.success("✓ Engine built — loading the live mode…")
                        st.rerun()
                    else:
                        st.error("Build failed — details below.")
                        st.code("\n".join(blog))
                st.caption("On **Streamlit Cloud** the build uses `packages.txt` (gfortran/make, bundled) + "
                           "`phydrus` (requirements.txt). Locally: `git submodule update --init "
                           "external/hydrus_source`, `make` in `source/`, `pip install phydrus`. "
                           "Until built, the tool falls back to the parametric model.")
            elif smiles and soil_cong is None:
                st.error("Could not parse the SMILES for the HYDRUS soil run — check the structure "
                         "or switch the compound to a curated congener.")
            else:
                if smiles:
                    st.caption(f"SMILES compound → soil Kd uses **{soil_cong}** "
                               f"({'read-across match' if how == 'match' else 'nearest curated congener by chain length'}); "
                               f"the plant uptake still uses your structure.")
                st.caption("Runs a **real HYDRUS-1D** paddy model (Richards + advection–dispersion + "
                           "linear Kd + root uptake) → Cwᵒ(t), Q_TP(t). Cached per setting.")
                f_oc = st.slider("Soil organic carbon f_oc", 0.005, 0.05, 0.02, 0.005,
                                 help="Kd = Koc(chain length)·f_oc → per-congener retardation R = 1+ρKd/θ.")
                flood_until = st.slider("Flooded until  [day]", 30, 120, 90, 5,
                                        help="Continuous flooding (clean irrigation) until this day, then drainage.")
                percolation = st.slider("Percolation excess  [cm/day]", 0.0, 1.0, 0.30, 0.05,
                                        help="Clean-water through-flow that leaches the dissolved pool.")
                season = 120.0
                try:
                    drivers = _hydrus_drivers_cached(soil_cong, season, f_oc, flood_until, percolation,
                                                     biomass=biomass)
                    st.success(f"HYDRUS-1D run complete — Cwᵒ(t) for {soil_cong} "
                               f"(Kd-retarded; mean-normalised).")
                except Exception as e:                                  # noqa: BLE001
                    st.error(f"HYDRUS run failed: {e}")

        elif mode == "Soil inventory → pore water":
            st.caption("Freundlich paddy soil S = K_F·C_wⁿ inverts a total soil load to pore water Cwᵒ(t).")
            C_total = st.number_input("Total soil inventory  [µg/kg dry]", 0.0, 1e4, 5.0, 0.5)
            K_F = st.slider("Freundlich K_F  [L/kg]", 0.2, 20.0, 2.0, 0.2,
                            help="Sorption capacity; long-chain PFAS sorb harder (higher K_F).")
            n_F = st.slider("Freundlich exponent n", 0.6, 1.0, 0.85, 0.01)
            theta_g = st.slider("Drained water content θ_g  [L/kg]", 0.2, 0.6, 0.35, 0.01)
            flood = st.checkbox("Flooded early season (dilution + leaching)", value=True)
            season = 120.0
            t = np.linspace(0.0, season, 241)
            flooded = t < (0.75 * season) if flood else None
            # Unified with the parametric flooded mode: default to the per-congener
            # HYDRUS-calibrated k_leach, same 0–0.15 range.
            k_leach = st.slider("Leaching rate k_leach  [1/day]", 0.0, 0.15,
                                float(api.default_k_leach(congener)), 0.0025,
                                help="Per-congener HYDRUS-calibrated default (short chains leach "
                                     "fast, long chains stay buffered ≈0).") if flood else 0.0
            Cwo, soil_obj = api.pore_water_from_inventory(
                t, C_total, K_F=K_F, n=n_F, theta_g=theta_g, flooded=flooded, k_leach=k_leach)
            drivers = api.drivers_from_arrays(t, Cwo, season=season, biomass=biomass)

        else:  # Biomonitoring
            st.caption("Enter MEASURED tissue concentrations + the pore-water/soil-solution Cwᵒ. "
                       "BAF is read straight off the data — no HYDRUS run needed.")
            src = st.radio("Input", ["Manual", "Upload CSV (tissue,conc[,Cwo])"], horizontal=True)
            if src == "Manual":
                bw = st.number_input("Pore-water Cwᵒ  [µg/L]", 1e-6, 1e4, 1.0, 0.1, format="%.4f")
                c_root = st.number_input("root conc  [µg/kg]", 0.0, 1e6, 0.49, 0.1)
                c_straw = st.number_input("straw conc  [µg/kg]", 0.0, 1e6, 0.83, 0.1)
                c_grain = st.number_input("grain conc  [µg/kg]", 0.0, 1e6, 0.46, 0.1)
                measured_bio = dict(conc={"root": c_root, "straw": c_straw, "grain": c_grain}, Cwo=bw)
            else:
                up = st.file_uploader("Biomonitoring CSV", type=["csv"])
                if up is not None:
                    try:
                        measured_bio = api.load_biomonitoring_csv(up)
                        st.success(f"Loaded {len(measured_bio['conc'])} tissues; Cwᵒ={measured_bio['Cwo']}.")
                    except Exception as e:                          # noqa: BLE001
                        st.error(f"Could not read CSV: {e}")
                else:
                    measured_bio = api.load_biomonitoring_csv(os.path.join(_EX, "biomonitoring_example.csv"))
                    st.info("Using examples/biomonitoring_example.csv (Yamazaki PFOA).")

        st.divider()
        compare = st.multiselect("Compare congeners (overlay)", api.CONGENERS,
                                 default=["PFBA", "PFOA", "PFDA", "PFOS"],
                                 help="Shown in the 'Compare' tab.")


# ---------------------------------------------------------------- header (title + disclaimer + intro)
st.title("🌾 PFAS in Rice — Uptake Explorer")
st.warning(_DISCLAIMER_KO if not expert else _DISCLAIMER)

if not expert:
    st.markdown(
        "**PFAS는 잘 분해되지 않는 '영원한 화학물질'입니다.** 이 도구는 벼가 오염된 논의 물·흙에서 "
        "PFAS를 얼마나 흡수하고 **어디에 쌓이는지** — 뿌리, 짚(줄기+잎), 먹는 **낟알** — 를 추정합니다.")
    st.markdown(
        "👉 **여기서 시작:** 왼쪽에서 **화학물질**과 **오염 수준**을 고른 뒤, 아래 "
        "**🗺️ 어디로 가나** 지도를 보세요. 화학 배경지식은 필요 없습니다.")
else:
    st.caption("Mechanistic 4-compartment dynamic model for permanently-anionic PFAS in paddy rice "
               "(IOC extension of the Trapp/Brunetti DPU). Parameters are measured/cited "
               "(docs/literature_db). Charts are interactive — hover, zoom, toggle. Outputs illustrative.")


# ---------------------------------------------------------------- custom tables (both modes)
custom_density = None
if use_custom_tables:
    _tbl_title = ("📋 내 데이터 표 — 성장 곡선 + 토양수 오염" if not expert
                  else "📋 Your data tables — growth curve + pore-water contamination")
    with st.expander(_tbl_title, expanded=True):
        _drv, custom_density = _render_custom_tables(biomass=biomass, Cwo_const=Cwo_const,
                                                     season0=season, key="ctbl", ko=not expert)
        if _drv is not None:
            drivers = _drv

# ---------------------------------------------------------------- run the model
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


# ============================================================================
#  SIMPLE MODE
# ============================================================================
if not expert:
    # ---- plain-language headline -------------------------------------------
    grain_c = float(res["conc"]["grain"][-1])
    root_c = float(res["conc"]["root"][-1])
    straw_c = float(res["straw"][-1])
    grain_baf = res["baf_final"]["grain"]
    tops = {"roots": root_c, "straw (stems + leaves)": straw_c, "grain": grain_c}
    where_most = max(tops, key=tops.get)

    st.subheader(f"{_cong_label_ko(congener)}")
    m1, m2, m3 = st.columns(3)
    m1.metric("뿌리 속", f"{root_c:.2g} µg/kg")
    m2.metric("짚(줄기+잎) 속", f"{straw_c:.2g} µg/kg")
    m3.metric("낟알(먹는 쌀) 속", f"{grain_c:.2g} µg/kg")

    _where_ko = {"roots": "뿌리", "straw (stems + leaves)": "짚(줄기+잎)", "grain": "낟알"}[where_most]
    _lead = ("입력하신 **성장 + 오염 표**를 바탕으로, " if use_custom_tables
             else f"선택하신 **{preset_word}** 오염 수준에서, ")
    st.info(
        _lead +
        f"이 모델은 벼 **낟알**에 {congener}가 약 **{grain_c:.2g} µg/kg** 들어 있을 것으로 추정합니다 "
        f"(토양수 농도의 약 **{grain_baf:.1f}배**). 대부분의 화학물질은 **{_where_ko}**에 남습니다.")
    st.caption("예시용 모델 추정치이며 — 식품안전·건강 판단이 아닙니다.")

    s_tabs = st.tabs(["🗺️ 어디로 가나", "📈 시간에 따른 축적",
                      "📊 얼마나 쌓이나", "🔎 거꾸로 추정", "ℹ️ 안내 & 용어"])

    # ---- Simple tab 1: the plant + soil map --------------------------------
    with s_tabs[0]:
        cc1, cc2 = st.columns([1, 1])
        animate = cc1.checkbox("▶ 한 철 재생", value=False,
                               help="벼가 자라는 동안 PFAS가 하루하루 쌓이는 모습을 봅니다.")
        day = cc2.slider("이앙 후 일수", float(res["t"][0]), float(res["t"][-1]),
                         float(res["t"][-1]), 1.0, disabled=animate)
        if animate:
            st.plotly_chart(plots.fig_schematic_animated(res, "conc", lang="ko"),
                            width="stretch", theme=None)
        else:
            ti = _nearest_index(res["t"], day)
            st.plotly_chart(plots.fig_schematic_from_res(res, "conc", ti, obs=None, lang="ko"),
                            width="stretch", theme=None)
        st.caption("벼와 논 흙을 실제 비율로 그렸습니다. **색이 진할수록(뜨거울수록) 그 부위에 PFAS가 더 많습니다.** "
                   "날짜 슬라이더를 끌거나 ▶를 눌러, 이앙부터 수확까지 화학물질이 어디에 쌓이는지 보세요.")

    # ---- Simple tab 2: build-up over time ----------------------------------
    with s_tabs[1]:
        st.plotly_chart(plots.fig_buildup_plain(res, lang="ko"), width="stretch")
        st.caption("한 철 동안 각 식물 부위의 PFAS 농도가 어떻게 변하는지. **낟알**은 형성된 뒤(개화 무렵)부터 "
                   "PFAS를 흡수하기 시작해 수확까지 계속 쌓입니다.")

    # ---- Simple tab 3: how much builds up ----------------------------------
    with s_tabs[2]:
        st.plotly_chart(plots.fig_where_plain(res, lang="ko"), width="stretch")
        st.caption("수확 시 각 부위의 PFAS 농도. 보통 뿌리에 가장 많이 남고, 먹는 낟알까지 "
                   "얼마나 도달하는지는 화학물질에 따라 다릅니다.")
        if obs:
            with st.expander("🔬 실제 측정값과 비교 (Yamazaki 2023)"):
                st.plotly_chart(plots.fig_baf(res, obs), width="stretch")
                st.caption("막대는 모델의 축적 배수를 출판된 온실 벼 연구(Yamazaki et al. 2023)의 측정값과 "
                           "비교한 것입니다. 막대가 비슷할수록 이 화학물질에 대해 모델이 실제 데이터와 잘 맞습니다.")

    # ---- Simple tab 4: work backwards (Bayesian inverse estimate) -----------
    with s_tabs[3]:
        _render_inverse_estimator(congener, E_m_mV=E_m, f_xy_source=fxy_source,
                                  biomass=biomass, key="inv_simple", simple=True)

    # ---- Simple tab 5: about & glossary ------------------------------------
    with s_tabs[4]:
        st.markdown(
            "### 이 도구가 하는 일\n"
            "벼가 흙에서 물과 녹아 있는 화학물질을 빨아들여 줄기 위로 올리고 낟알에 저장하는 과정을 "
            "**메커니즘 모델**로 계산합니다. PFAS 화학물질과 논의 오염 정도를 주면, 한 철 동안 뿌리·짚·"
            "먹는 낟알에 쌓이는 양을 추정합니다.\n\n"
            "### 보는 법\n"
            "- **🗺️ 어디로 가나** — 식물 그림; 색이 뜨거울수록 PFAS가 많음.\n"
            "- **📈 시간에 따른 축적** — 이앙부터 수확까지 농도 변화.\n"
            "- **📊 얼마나 쌓이나** — 최종 농도와 실제 측정값과의 비교.\n"
            "- **🔎 거꾸로 추정** — 실험실 측정값이 있으면 토양수 오염도를 (불확실성 범위와 함께 — "
            "베이지안 추정) 역추정.\n\n"
            "### 쉬운 용어 사전")
        st.markdown(_glossary_md(ko=True))
        st.warning(_DISCLAIMER_KO)
        st.caption("수식·파라미터·토양 결합(HYDRUS-1D)·구조(SMILES) 입력이 필요하면 사이드바의 "
                   "**전문가/고급 모드**를 켜세요.")

    # ---- downloads (Simple) ------------------------------------------------
    with st.expander("⬇️ 결과 내려받기"):
        cda, cdb = st.columns(2)
        cda.download_button("요약 표 (CSV)", api.summary_csv(res, obs, bio_baf),
                            file_name=f"{congener}_summary.csv", mime="text/csv")
        cdb.download_button("전체 시계열 (CSV)", api.timeseries_csv(res),
                            file_name=f"{congener}_timeseries.csv", mime="text/csv")
        png, why = _png_bytes(plots.fig_schematic_from_res(res, "conc", -1, lang="ko"))
        if png is not None:
            st.download_button("식물 지도 (PNG)", png, file_name=f"{congener}_map.png", mime="image/png")
        else:
            st.caption(f"PNG 이미지 내보내기는 선택 패키지 `kaleido`(+Chrome)가 필요합니다. "
                       f"CSV 내려받기는 그것 없이도 됩니다. ({why})")

# ============================================================================
#  EXPERT MODE  (the full research interface)
# ============================================================================
else:
    # ---------------------------------------------------------------- headline metrics
    st.subheader(f"{congener}  (C{p['n_C']} {p['group']})  ·  source: {mode}")
    c1, c2, c3, c4 = st.columns(4)
    for col, tis in zip((c1, c2, c3), ("root", "straw", "grain")):
        pred = res["straw_baf"] if tis == "straw" else res["baf_final"][tis]
        sub = f"obs {obs[tis]:.2f}" if tis in obs else "no obs"
        if bio_baf and tis in bio_baf:
            sub = f"measured {bio_baf[tis]:.2f}"
        col.metric(f"{tis} BAF [L/kg]", f"{pred:.2f}", sub, delta_color="off")
    c4.metric("anion exclusion eᴺ", f"{res['eN']:.0f}", f"N={res['N']:.2f}", delta_color="off")

    # ---- structure → parameters panel (SMILES mode) -------------------------
    if desc is not None:
        kind = "read-across (matches a curated congener)" if desc.matched_name else "QSPR (novel structure)"
        with st.expander(f"🧬 Structure → parameters  ·  {kind}", expanded=True):
            if provisional:
                st.warning("⚠ **PROVISIONAL** — novel / non-calibrated structure. Binding (K_PL/K_prot) "
                           "is from the QSPR and translocation `f_xy` is a head-group estimate, so treat "
                           "the magnitudes as indicative (ordering is more reliable than absolute level).")
            else:
                st.success(f"Matches **{desc.matched_name}** — uses the curated, measured-anchored parameters.")
            dc1, dc2 = st.columns(2)
            dc1.markdown(
                f"**Structure** &nbsp;`{desc.formula}`, MW {desc.mol_weight:.1f}\n\n"
                f"- perfluoro-C: **{desc.n_perfluoroC}** &nbsp;·&nbsp; head group: **{desc.head_group}**\n"
                f"- ether-O: {desc.n_ether_O} &nbsp;·&nbsp; CF₃: {desc.n_CF3} &nbsp;·&nbsp; "
                f"branched: {desc.branched} &nbsp;·&nbsp; linear: {desc.is_linear}\n"
                f"- read-across match: **{desc.matched_name or 'none (novel)'}**")
            dc2.markdown(
                f"**Model parameters**\n\n"
                f"- K_PL = {p['K_PL']:.0f} &nbsp;·&nbsp; K_prot = {p['K_prot']:.0f} &nbsp;·&nbsp; "
                f"K_cw = {p['K_cw']:.0f}  [L/kg]\n"
                f"- f_xy = {p['f_xy']:.4g} &nbsp;·&nbsp; L_Ph = {p['L_Ph']:.3g} &nbsp;·&nbsp; "
                f"κ_d = {p['kappa_d']:.3g}\n"
                f"- B_k: " + ", ".join(f"{k} {v:.1f}" for k, v in res["B_k"].items()))
            for n in getattr(desc, "notes", []):
                st.caption("• " + n)

    is_biomon = (mode == "Biomonitoring (measured tissue)" and measured_bio)
    _biomon_note = ("⚠ This chart shows a **model reference run** (Cwᵒ=1), **not** your measured "
                    "biomonitoring data. Your measured values are on the **map** and **BAF** tabs.")

    tabs = st.tabs(["🗺️ Plant & soil map", "📈 Tissue dynamics", "🟫 Soil & drivers",
                    "📊 BAF vs observed", "🔗 Chain-length trends", "⚖️ Compare congeners",
                    "✅ Tang TF (OOS)", "🔎 Inverse (Bayesian)", "ℹ️ About / coupling"])

    # ---- Tab 1: the plant + soil accumulation map ---------------------------
    with tabs[0]:
        if is_biomon:
            vals = dict(measured_bio["conc"])
            finite = [v for v in vals.values() if v is not None and np.isfinite(v)]
            cmax = max(finite) if finite else 1.0
            fig = plots.fig_plant_schematic(
                vals, cmin=0.0, cmax=cmax, label="measured conc [µg/kg]",
                Cwo=measured_bio.get("Cwo"), title=f"{congener} — measured tissue map")
            st.plotly_chart(fig, width="stretch", theme=None)
            st.caption("Compartments coloured by the MEASURED concentration. Stem/leaf share the "
                       "straw colour when only straw is reported.")
        else:
            cc1, cc2, cc3 = st.columns([1.3, 1, 1])
            metric = cc1.radio("Colour by", ["concentration", "BAF"], horizontal=True)
            metric_key = "baf" if metric == "BAF" else "conc"
            animate = cc2.checkbox("▶ animate season", value=False,
                                   help="Autoplay the accumulation through the season.")
            day = cc3.slider("Day after transplant", float(res["t"][0]), float(res["t"][-1]),
                             float(res["t"][-1]), 1.0, disabled=animate)
            if animate:
                st.plotly_chart(plots.fig_schematic_animated(res, metric_key), width="stretch", theme=None)
            else:
                ti = _nearest_index(res["t"], day)
                st.plotly_chart(plots.fig_schematic_from_res(res, metric_key, ti, obs=None),
                                width="stretch", theme=None)
            st.caption("Each compartment is filled by its accumulation on a shared colorbar — drag the "
                       "day slider (or hit ▶) to watch where PFAS builds up. Leaf is xylem-terminal, "
                       "grain is phloem-fed; the root retains the anion (low f_xy).")

    # ---- Tab 2: tissue concentration dynamics -------------------------------
    with tabs[1]:
        if is_biomon:
            st.warning(_biomon_note)
        st.plotly_chart(plots.fig_tissue(res), width="stretch")
        st.plotly_chart(plots.fig_burden(res), width="stretch")
        st.markdown(f"**B_k [L/kg fw]** — " + ", ".join(f"{k}: {v:.2f}" for k, v in res["B_k"].items())
                    + f"  ·  f_xy={p['f_xy']:.4g}, L_Ph={p['L_Ph']:.4g}, κ_d={p['kappa_d']:.3g}")
        st.caption("Top: tissue **concentration** C_k(t) [µg/kg] (intensive). Bottom: **PFAS mass** "
                   "(burden) = C_k(t)·M_k(t) [µg/hill] (extensive) — where the chemical actually ends up. "
                   "A tissue can be high-concentration yet low-mass (small organ), so the two views differ; "
                   "the terminal leaf/grain keep gaining mass as the organ grows. (Organ *biomass* M_k(t) is "
                   "in the **🟫 Soil & drivers** tab.) The **grain takes up no PFAS until it forms (~flowering)** "
                   "and then accumulates — its loading is *formation-gated* (the panicle is absent before then, "
                   "so no solute enters it; DPU-consistent), and the empty pre-set period is not plotted.")

    # ---- Tab 3: soil & drivers ----------------------------------------------
    with tabs[2]:
        if is_biomon:
            st.warning(_biomon_note)
        st.plotly_chart(plots.fig_drivers(res), width="stretch")
        if soil_obj is not None:
            st.plotly_chart(plots.fig_isotherm(soil_obj, Cwo_now=float(res["Cwo"][-1])),
                            width="stretch")
        st.plotly_chart(plots.fig_soil_profile(res, profile=profile), width="stretch")
        st.caption("Cwᵒ(t) is the only PFAS-specific driver; Q_TP(t) and M(t) are crop physiology "
                   "(measured FAO-56 transpiration + the selected biomass driver) reused across modes.")

    # ---- Tab 4: BAF vs observed ---------------------------------------------
    with tabs[3]:
        if bio_baf:
            model_b = {"root": res["baf_final"]["root"], "straw": res["straw_baf"],
                       "grain": res["baf_final"]["grain"]}
            st.plotly_chart(plots.fig_biomon_compare(bio_baf, model_b), width="stretch")
            st.caption("Measured biomonitoring BAF (tissue conc ÷ pore water) vs the model prediction.")
        else:
            extra = None
            if spec == "Curated congener" and congener is not None:
                show_tp = st.checkbox(
                    "Overlay the two-pool (seq) model — EXPLORATORY", value=True,
                    help="The sequestration two-pool root model (model_api.simulate_twopool_seq; "
                         "docs/twopool_root_exploration.md): a mobile pool + an irreversible non-K_PL "
                         "k_seq sink. Shown at its calibrated operating point (Cwo=1, season≈120, demo "
                         "forcings) so it is comparable to the fixed Yamazaki bars. In-sample / opt-in; "
                         "parameters.json and the canonical core are unchanged.")
                if show_tp:
                    tp = _simulate_twopool_seq(congener)
                    if tp is not None:
                        extra = {"two-pool (seq, exploratory)": tp}
            # The overlay is an optional EXPLORATORY nicety -- never let it crash the BAF
            # tab (e.g. a stale/old plots.py without the `extra` param right after a deploy,
            # before the module cache refreshes). Fall back to the plain core-vs-observed plot.
            try:
                fig = plots.fig_baf(res, obs, extra=extra)
            except Exception:                                       # noqa: BLE001
                fig = plots.fig_baf(res, obs)
                extra = None
            st.plotly_chart(fig, width="stretch")
            if extra is not None:
                st.caption("🧪 **two-pool (seq)** is EXPLORATORY / in-sample (Yamazaki fit): it captures the "
                           "long-chain **root** BAF and the PFOS/PFUnDA split the single-pool core misses, while "
                           "keeping the monotone physical f_xy — overall log10 RMSE 0.251. Run at its calibrated "
                           "point, so it does **not** track the sidebar f_xy/Cwᵒ/biomass (unlike the 4-pool core "
                           "bar). The *carrier* two-pool (`close_longchain_2pool`, saturated DOF-0 closure) is "
                           "API-only — it reproduces the observed bars by construction and is too slow (~1 min) to "
                           "render live.")
            if not obs:
                st.info("No Yamazaki BAF for this congener (model prediction only).")
            else:
                with st.expander("ℹ️ What are the Yamazaki 2023 conditions? (and how to match them)"):
                    st.markdown(
                        "**Yamazaki et al. 2023**, *Environ. Sci. Technol.* **57** "
                        "([doi:10.1021/acs.est.2c08767](https://doi.org/10.1021/acs.est.2c08767)).\n\n"
                        "- **Design**: greenhouse **pot** study, Japanese **Andosol** soil; each congener "
                        "spiked **individually** with clean irrigation water; grown a **full cycle** to maturity.\n"
                        "- **Cultivars**: the plotted observed BAFs are the **geomean of Indica + Japonica** "
                        "(SI tables S16/S18/S19), fresh-weight.\n"
                        "- **BAF definition**: **tissue conc ÷ pore-water (soil-solution) conc** [L/kg], for "
                        "**root / straw (stem+leaf) / grain (brown rice)** — the same definition the model reports.\n\n"
                        "**These observed points are fixed measurements** — they do **not** move when you change "
                        "the sidebar. The model was calibrated to reproduce them at one operating point, so the "
                        "overlay is a like-for-like comparison **only at those settings**:")
                    st.markdown(
                        "| sidebar control | match-Yamazaki value |\n|---|---|\n"
                        "| Data source | **Model (parametric)**, *Measured forcings* on |\n"
                        "| Pore-water Cwᵒ | **1.0 µg/L** (so tissue conc = BAF) |\n"
                        "| Root→shoot f_xy | **W2 fit (reproduces Yamazaki)** |\n"
                        "| E_m | **−120 mV** (default) · Season ~**120 d** |")
                    st.caption("Changing E_m / f_xy / Cwᵒ / season, or a non-parametric source (HYDRUS, soil "
                               "inventory…), moves the model AWAY from the Yamazaki experiment — then the overlay "
                               "is a reference trend, not a calibrated match. (Yamazaki = clean per-congener Andosol "
                               "pot study; your scenario may differ in soil, exposure and crop.)")

    # ---- Tab 5: chain-length trends -----------------------------------------
    with tabs[4]:
        key = st.selectbox("Parameter", ["K_PL", "K_prot", "K_cw_root",
                                         "f_xy_recommended", "B_root", "B_grain"], index=0)
        # a novel SMILES compound is not in the curated chain series -> ring a reference instead
        chain_cong = congener if congener in api.CONGENERS else (
            desc.matched_name if desc and desc.matched_name in api.CONGENERS else "PFOA")
        st.plotly_chart(plots.fig_chain(api.chain_table(), chain_cong, key), width="stretch")
        if congener not in api.CONGENERS:
            st.caption(f"'{congener}' is a novel structure (not in the curated 13); the ring marks "
                       f"**{chain_cong}** for reference. Its own parameters are in the 🧬 panel above.")

    # ---- Tab 6: compare congeners -------------------------------------------
    with tabs[5]:
        tissue = st.radio("Tissue", ["root", "straw", "grain"], index=1, horizontal=True)
        if compare:
            if drivers is not None:
                dt = _drivers_tuple(drivers)
                results = {nm: _simulate(nm, drivers_tuple=dt, **sim_kw) for nm in compare}
            else:
                results = {nm: _simulate(nm, Cwo=Cwo_const, season=season,
                                         measured_forcing=measured, **sim_kw) for nm in compare}
            st.plotly_chart(plots.fig_compare(results, tissue), width="stretch")
        else:
            st.info("Select congeners in the sidebar to compare.")

    # ---- Tab 7: Tang TF (OOS) -----------------------------------------------
    with tabs[6]:
        st.markdown("**Tang 2026** (flooded paddy, Nipponbare, 150 d; PFOA/PFOS/GenX) — per-organ "
                    "transfer factor **TF = C_organ/C_root**, an *out-of-sample* check of the root→shoot "
                    "loading `f_xy` (only Tang's head-group *sign* went into the model build).")
        if smiles or congener not in api.TANG_CONGENERS:
            st.info(f"Tang 2026 covers **{', '.join(api.TANG_CONGENERS)}** only — "
                    "pick one of these curated congeners (sidebar) to see the comparison.")
        else:
            c1, c2 = st.columns(2)
            dose = ("low" if c1.radio("Tang dose", ["across-dose mean", "0.1 µg/g (env-closest)"],
                                      horizontal=True).startswith("0.1") else "mean")
            bm = ("oryza" if c2.checkbox("ORYZA biomass driver", value=False,
                                         help="Drive the shoot model with the mechanistic ORYZA2000 "
                                              "biomass instead of the logistic (slower).") else "growth_rice")
            val = api.tang_tf_validation(congener, dose=dose, biomass=bm)
            valr = api.tang_tf_validation(congener, dose=dose, biomass=bm, use_refit=True)
            st.plotly_chart(plots.fig_tang_tf(val, valr), width="stretch")
            _verdict = {"GenX": "GenX over-predicted ~12× by the provisional 0.233; the refit ≈ the documented 0.013.",
                        "PFOS": "PFOS is **dataset-dependent** — Yamazaki W2 0.142 vs Tang ~0.32; report the **range with conditions**, not one value.",
                        "PFOA": "PFOA f_xy is dose-sensitive (0.064 across-dose mean → 0.097 at 0.1 µg/g)."}[congener]
            st.markdown(
                f"- **f_xy**: recommended **{val['f_xy']:.3g}** → Tang-refit **{valr['f_xy']:.3g}**. {_verdict}\n"
                f"- **Dry-weight basis** — `TF_dw = TF_fw·(1−θ_root)/(1−θ_tissue)`; comparing the model's "
                f"fresh-weight TF to Tang's dry-weight TF (without this factor) is a units error that "
                f"flatters the grain ~8×.\n"
                f"- **Grain/endosperm is structurally under-predicted ~3–8×** and is **not** closable by "
                f"`L_Ph`/lipid — a phloem-delivery + dry-weight limit, not a calibration gap.")
        st.caption("Source: docs/literature_db/raw_si/tang2026_doseresponse.csv (SI S7/S8). Model: "
                   "redistributed-shoot `simulate_nstem_leaf`; the f_xy refit is OVERRIDE-only (parameters.json "
                   "unchanged). Details: docs/VALIDATION_TANG2026_NSTEM_KR.md · docs/tang2026_grain_units_exploration.md.")

    # ---- Tab 8: Inverse (Bayesian exposure estimate) ------------------------
    with tabs[7]:
        st.markdown("**Bayesian inverse** — estimate the pore-water exposure Cwᵒ from measured "
                    "tissue concentrations (Laplace posterior in log Cwᵒ; the well-posed "
                    "direction of `validation/bayesian_inverse_demo.py`). Transport is held at "
                    "the sidebar defaults.")
        _render_inverse_estimator(congener, E_m_mV=E_m, f_xy_source=fxy_source,
                                  biomass=biomass, key="inv_expert", simple=False)
        st.caption("Identifiability: only the EXPOSURE level is estimated here. From tissue data "
                   "alone Q_TP·f_xy is a product ridge and Cwᵒ vs root-uptake conductance is "
                   "degenerate, so pinning transport absolutely needs an independent measurement "
                   "(xylem sap / pore-water probe). See docs + validation/bayesian_inverse_demo.py.")

    # ---- Tab 9: About / coupling --------------------------------------------
    with tabs[8]:
        st.markdown(
            """
### Five ways to drive the plant model

| Mode | PFAS driver Cwᵒ(t) comes from | Q_TP, M(t) | When to use |
|---|---|---|---|
| **Model (parametric)** | a constant you set | measured FAO-56 / ORYZA | quick what-ifs, teaching |
| **HYDRUS / CSV drivers** | a HYDRUS-1D / Phydrus run (CSV) | from the same CSV, or measured | you have a calibrated soil-water-solute model |
| **Run HYDRUS-1D (live)** | a real HYDRUS-1D run executed here | HYDRUS root uptake + ORYZA | you want the engine to run in-app (needs it built) |
| **Soil inventory** | inverting a soil load with a Freundlich isotherm | measured | you know total soil PFAS but not pore water |
| **Biomonitoring** | a measured pore-water / soil-solution value | — (not needed) | you have field tissue + water concentrations |

**Live HYDRUS-1D** runs the genuine engine (built from `external/hydrus_source`) via `phydrus`:
a one-season paddy model (Richards flow + advection–dispersion + linear Kd sorption + root
water uptake) gives the **congener-dependent** pore water — weakly-sorbed short chains leach
under flooding (Cwᵒ falls), strongly-sorbed long chains stay buffered (flat). Kd comes from the
Koc(chain-length) QSPR. Build it with `git submodule update --init external/hydrus_source`,
`make` in `source/` (needs gfortran), and `pip install phydrus`; the mode auto-detects the engine.

### HYDRUS-1D coupling (Method A, one-way) — inputs & outputs

**HYDRUS *inputs* you set up (soil side):** soil-hydraulic parameters (van Genuchten),
the atmospheric boundary condition (precip / irrigation / evaporation), the root
water-uptake distribution, solute-transport parameters (a linear K_d **or** the
Freundlich K_F, n; dispersivity), and the initial / boundary PFAS concentration.

**HYDRUS *outputs* this tool consumes** (map them into the driver CSV columns):

| CSV column | HYDRUS-1D source | meaning |
|---|---|---|
| `t` | output times | day after transplant |
| `Cwo` | `Conc` at the root-zone node (`Obs_Node.out` / `solute1.out`) | pore-water free anion [µg/L] |
| `Qtp` | `vRoot` (root water uptake) / `T_act` (`T_Level.out`) | transpiration stream [L/day] |
| `M_root,M_stem,M_leaf,M_grain` | a plant **growth** sub-model (not HYDRUS) | organ fresh mass [kg] |

The plant ODE is solved in Python; HYDRUS is not modified (the *tight* Method B —
embedding the root-uptake term `j_R` in the HYDRUS FORTRAN — is future work).

### Biomonitoring — when HYDRUS is unnecessary
If you already have **measured tissue concentrations** and a **measured pore-water
(or soil-solution) concentration**, the bioaccumulation factor is just
`BAF = C_tissue / Cwᵒ` — no transport simulation required. Use this mode to read
BAFs straight off field data and to overlay them on the model for a sanity check.

### Plain-language glossary
""")
        st.markdown(_glossary_md())
        st.caption("Open item (recorded, not modelled): long-chain PFCA shoot accumulation "
                   "(hysteretic sorption) — see docs/nstem_gradient_exploration.md.")

    # ---- downloads (Expert) -------------------------------------------------
    with st.expander("⬇️ Download results"):
        cda, cdb = st.columns(2)
        cda.download_button("BAF summary (CSV)", api.summary_csv(res, obs, bio_baf),
                            file_name=f"{congener}_summary.csv", mime="text/csv")
        cdb.download_button("Driver + tissue time series (CSV)", api.timeseries_csv(res),
                            file_name=f"{congener}_timeseries.csv", mime="text/csv")
        png, why = _png_bytes(plots.fig_schematic_from_res(res, "conc", -1))
        if png is not None:
            st.download_button("Plant map (PNG)", png, file_name=f"{congener}_map.png", mime="image/png")
        else:
            st.caption(f"PNG figure export needs the optional `kaleido` package (+ Chrome). CSV works "
                       f"without it. ({why})")


# ---------------------------------------------------------------- footer (every screen)
st.divider()
st.caption(_DISCLAIMER_KO if not expert else _DISCLAIMER)
fc1, fc2, fc3 = st.columns(3)
fc1.caption(f"**PFAS–Rice Uptake Model** · v{APP_VERSION}")
if not expert:
    fc2.caption(f"[소스 코드]({REPO_URL}) · [문서]({DOCS_URL})")
    fc3.caption("인용: PFAS–Rice Compartmental Uptake Model (Trapp/Brunetti DPU의 "
                "이온성유기화합물 확장), 2026.")
else:
    fc2.caption(f"[Source code]({REPO_URL}) · [Documentation]({DOCS_URL})")
    fc3.caption("How to cite: PFAS–Rice Compartmental Uptake Model (IOC extension of the "
                "Trapp/Brunetti DPU framework), 2026.")
