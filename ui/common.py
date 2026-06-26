"""Shared UI layer for the PFAS-rice dashboard: constants, cached model helpers,
and the render building blocks that app.py assembles. Split out of the monolithic
app.py (HANDOFF P3-1); behaviour is unchanged."""
import os

import numpy as np
import streamlit as st

import model_api as api
import plots

# repo root (this file lives in <root>/ui/), and the bundled example CSVs
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_EX = os.path.join(_ROOT, "examples")


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
        prog = st.progress(0.0, text="오염 수준 추정 준비 중…" if ko else "Preparing the estimate…")

        def _cb(done, total):
            prog.progress(min(done / total, 1.0),
                          text=(f"모델을 거꾸로 계산하는 중… ({done}/{total} 단계)" if ko
                                else f"Running the model backwards… (step {done}/{total})"))
        try:
            est = api.estimate_exposure_bayesian(
                congener, meas, sigma_log10=sigma, E_m_mV=E_m_mV,
                f_xy_source=f_xy_source, biomass=biomass, progress=_cb)
        except Exception as e:                               # noqa: BLE001
            prog.empty()
            st.error((f"추정을 실행할 수 없습니다: {e}") if ko else f"Could not run the estimate: {e}")
            return
        prog.empty()
        st.session_state[res_key] = (sig, est)
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


# ---------------------------------------------------------------- render building blocks
def render_header(cfg):
    """Title + disclaimer + intro (both modes)."""
    expert = cfg.expert
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
        _tbl_title = ("📋 내 데이터 표 — 성장 곡선 + 토양수 오염" if not expert
                      else "📋 Your data tables — growth curve + pore-water contamination")
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
