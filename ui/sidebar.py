"""The sidebar (Simple vs Expert) -> a config namespace consumed by app.py.
Split out of app.py (HANDOFF P3-1); behaviour is unchanged."""
import os
from types import SimpleNamespace

import numpy as np
import streamlit as st

import model_api as api
import plots

from ui.common import (_EX, _cong_label, _cong_label_ko, _PRESETS_KO, _SCENARIOS_KO, _mol_svg,
                       _hydrus_drivers_cached, _hydrus_soil_congener)



def _neutral_panel():
    """Sidebar inputs for a NEUTRAL organic -> the kwargs of `model_api.simulate_neutral`.

    A neutral compound has no congener: the ONE required input is a log Kow, which
    drives both `K_PW` (partition) and `TSCF` (translocation) through published QSPRs
    with nothing fitted. Everything else here is a scope switch or an optional
    mechanism, and each carries how far it is actually tested.
    """
    n = {}
    n["log_kow"] = st.number_input("log Kow", value=2.45, min_value=-2.0, max_value=8.0, step=0.05,
                                   help="The one required input. Drives K_PW and TSCF. The "
                                        "measured tables behind this path span about -0.7 to 5.4; "
                                        "outside that you are extrapolating.")
    n["name"] = st.text_input("compound name", value="carbamazepine") or "neutral"
    hl = st.number_input("in-planta half-life [d]  (0 = none)", value=7.0, min_value=0.0,
                         max_value=365.0, step=1.0,
                         help="STRONGLY recommended. With no metabolism the leaf is an "
                              "unbounded terminal accumulator, so a run without it is an UPPER "
                              "BOUND. Do not treat one number as a compound constant: Kodesova "
                              "measured the parent fraction varying 4.8x BETWEEN SPECIES for the "
                              "same compound, and soil persistence is not plant persistence.")
    n["half_life"] = float(hl) or None
    n["tscf_model"] = st.radio("TSCF QSPR", ["briggs", "schriever"], horizontal=True,
                               help="Briggs 1982 (narrow bell, peak 0.784 @ logKow 1.78) vs the "
                                    "Schriever 2020 refit (97 values, ~3x broader). TSCF is an "
                                    "INPUT here, not a fitted parameter, so the gap between them "
                                    "measures how well it is known.")
    with st.expander("🧪 Composition, phloem, air"):
        n["lipid_source"] = st.radio(
            "root lipid reading", ["measured", "briggs_anchor"], horizontal=True,
            help="An OPEN question, so it is a switch. 'measured' = 1% fresh weight from "
                 "measured cereal roots (default, and what every published number here is "
                 "on); 'briggs_anchor' = 2.47%, what Briggs' 1982 regression implies — he "
                 "measured no lipid at all. 2.5x apart for lipophilic compounds, ~equal "
                 "below log Kow 1. Evidence is 3 measured tables to 1 against the anchor.")
        n["phloem"] = st.checkbox(
            "phloem ON (departs from the base)", value=False,
            help="The neutral base excludes phloem transport — it is an addition of the "
                 "ionisable extension. Turning it on drives the small terminal grain hard: "
                 "a statement about assuming unrestricted loading, not a prediction.")
        n["air"] = st.checkbox(
            "plant–air exchange (volatilisation + gaseous uptake)", value=False,
            help="Off by default: needs K_AW and a molar mass, which the strict Kow-only "
                 "a-priori run does not use. Identically zero at K_AW = 0. Leaving it off "
                 "for a volatile compound makes the leaf an upper bound BY CONSTRUCTION.")
        if n["air"]:
            a1, a2 = st.columns(2)
            n["MW"] = a1.number_input("MW [g/mol]", value=236.3, min_value=1.0, step=1.0)
            n["K_AW"] = a2.number_input("K_AW [-]", value=1e-3, min_value=0.0, max_value=10.0,
                                        step=1e-4, format="%.2e",
                                        help="Dimensionless Henry's-law constant. 0 = the air "
                                             "pathway is structurally absent, not just small.")
            n["air_kw"] = dict(C_air=float(st.number_input(
                "ambient C_air [µg/m³]", value=0.0, min_value=0.0, step=0.1,
                help="0 = clean air, i.e. volatilisation only.")))
    with st.expander("⚗️ Weak electrolyte (pKa) + apoplastic bypass"):
        we = st.checkbox("this compound is an acid / base (has a pKa)", value=False,
                         help="OFF = the strictly neutral path. ON = the compound is a neutral "
                              "molecule AND an ion at once, weighted by Henderson-Hasselbalch; "
                              "the ion feels the GHK membrane term, the neutral species does not.")
        if we:
            # two columns, not three: the sidebar is narrow enough that a third
            # squeezes the acid/base radio into one letter per line
            w1, w2 = st.columns(2)
            n["pKa"] = w1.number_input("pKa", -5.0, 14.0, 4.5, 0.1)
            n["pH"] = w2.number_input("root-zone pH", 3.0, 10.0, 6.5, 0.1)
            n["is_acid"] = st.radio("acid or base", ["acid", "base"], horizontal=True) == "acid"
            import literature_params as LP
            f_n, f_d = LP.speciation(float(n["pKa"]), float(n["pH"]), bool(n["is_acid"]))
            st.caption(
                f"→ neutral fraction **f_n = {f_n:.3g}**, ionic **f_d = {f_d:.3g}**. "
                + ("The ion is an ANION: excluded by the inside-negative membrane (the PFAS "
                   "case at f_n→0)." if n["is_acid"] else
                   "The ion is a CATION: ATTRACTED by the inside-negative membrane, not "
                   "excluded — at pKa 4.5 / pH 6.5 / log Kow 2.45 the base's root BAF is "
                   "1.51 against the acid's 0.079 (~19×)."))
        n["g_apo"] = st.number_input(
            "apoplastic bypass g_apo  [L/kg/d]", 0.0, 50.0, 0.0, 0.5,
            help="A route AROUND the membrane, so it feels neither speciation nor GHK. "
                 "0 = structurally absent (the default — NOTHING is adopted).")
        st.caption("Speciation is TESTED and BOUNDED: direction supported (Spearman +0.480 on "
                   "Schriever's 67 ionisable rows; rank +0.284 → +0.520), magnitude refuted "
                   "(bias −0.203, ~nothing predicted below f_n ≈ 1e−3). A small g_apo ≈ 0.5–1 "
                   "improves rank AND scale; the RMSE-optimal 5 wrecks the ordering, so do not "
                   "fit it on RMSE. §4l–4m of docs/neutral_dpu_validation.md.")
    # Soil sorption: needed only by the exposure modes that model the soil. The PFAS
    # chain-length Koc QSPR does not apply, so it comes from log Kow (Karickhoff 1981,
    # PROVISIONAL -- nothing in this repo scores a predicted Koc) and is editable.
    import literature_params as LP
    koc_default = float(LP.koc_neutral(float(n["log_kow"])))
    n["Koc"] = st.number_input(
        "soil Koc  [L/kg]", 0.1, 1e6, koc_default, koc_default / 10.0, format="%.1f",
        key=f"koc_neutral_{n['log_kow']:.2f}",
        help=f"Used by the flooded Cwᵒ(t) shape and the live HYDRUS run (Kd = Koc·f_oc). "
             f"Default {koc_default:.1f} = Karickhoff 1981 (log Koc = 0.989·log Kow − 0.346) — "
             f"PROVISIONAL: no table in this repo scores a predicted Koc, and Li 2019's soil "
             f"half shows the estimated sorption term is where the error collects (bias +0.033 "
             f"with a measured K_om vs +0.291 with an estimated one). Type a MEASURED Kd/f_oc "
             f"here whenever you have one.")
    return n


def _mechanism_panel():
    """The PFAS mechanism switches (sidebar expander) -> (uptake, lipid_loading,
    vmax_scale, g_apo, km_scale). Both switches are OPEN questions, so they are named
    modes with the shipped model as the default, never silent constants.

    PFAS-ONLY, and not by omission: the carrier exists to overcome anion exclusion and
    the lipid term is K_PL-gated on PFAS binding, so neither has a neutral counterpart
    (the neutral path's own entry switch, the apoplastic bypass, lives in its panel).
    """
    vmax_scale = g_apo = None
    km_scale = 1.0
    with st.expander("⚙️ Mechanism (advanced)"):
        up_label = st.radio(
            "Root entry", ["carrier (default — shipped)", "bypass (apoplastic)"],
            help="How the anion gets past the inside-negative membrane. The "
                 "Michaelis–Menten CARRIER is this repo's one addition to the Trapp "
                 "cell model — and it keeps its place by DEFAULT, NOT BY EVIDENCE: "
                 "on Yamazaki the carrier (1.035) and a single global apoplastic "
                 "bypass g_apo=20 (0.996) are indistinguishable (bootstrap 0.749), "
                 "while adding nothing at all fails (2.640). See "
                 "validation/carrier_vs_bypass.py.")
        uptake = "carrier" if up_label.startswith("carrier") else "bypass"
        lipid_loading = st.checkbox(
            "Lipid-facilitated loading (K_PL-gated)", value=False,
            help="Adds the B-independent bound-loading term g_xy·C / g_ph·C: the free "
                 "anion Cw=C/B starves high-binding long chains, but the membrane-lipid-"
                 "bound pool still rides the transpiration stream. Constants were fit on "
                 "YAMAZAKI, so transferring them is genuine out-of-sample — and it is the "
                 "repo's strongest cross-dataset result (Tang per-organ 1.232→0.516; best "
                 "on Kim 2019 grain too). EXPLORATORY, default off. NOTE it supplies its "
                 "own f_xy, so the 'Root→shoot f_xy' choice above is ignored while it is on.")
        st.caption("Both default to the shipped model; turning either on changes every "
                   "tab (map, dynamics, BAF), never `parameters.json`.")
        adv = st.checkbox("Override the entry constants", value=False,
                          help="Scan levers from validation/carrier_vs_bypass.py and "
                               "validation/dose_series_carrier.py. Leave off to use the "
                               "mode's own values.")
        if adv:
            vmax_scale = st.number_input(
                "Vmax ×", 0.0, 100.0, float(api.UPTAKE_MODES[uptake]["vmax_scale"]), 0.5,
                # key includes the mode so switching carrier<->bypass RESEEDS the
                # override with that mode's own value instead of keeping a stale one
                key=f"vmax_scale_{uptake}",
                help="Carrier capacity multiplier (0 = carrier off). The long-chain work "
                     "needed ~5× for PFDoDA and that enhancement is NOT QSPR-able (LC6).")
            km_scale = st.number_input(
                "Km ×", 0.01, 1000.0, 1.0, 1.0,
                help="Half-saturation multiplier — the model's ONLY nonlinearity in "
                     "exposure. Tang's 5-dose series needs Km ≥ 500 µg/L (100× the fitted "
                     "5) to be as flat as measured, i.e. linear over the whole span, which "
                     "IS the bypass's functional form (validation/dose_series_carrier.py).")
            g_apo = st.number_input(
                "g_apo  [L/kg/d]", 0.0, 100.0, float(api.UPTAKE_MODES[uptake]["g_apo"]), 0.5,
                key=f"g_apo_{uptake}",
                help="Apoplastic bypass: a route AROUND the membrane, so it feels neither "
                     "speciation nor the GHK factor. 0 = structurally absent (the default; "
                     "nothing is adopted).")

    return uptake, lipid_loading, vmax_scale, g_apo, km_scale


def build():
    """Render the sidebar and return a SimpleNamespace with the scenario config."""
    cfg = SimpleNamespace()
    with st.sidebar:
        expert = st.toggle(
            "🔬 전문가/고급 모드 (Expert / advanced)", value=False,
            help="끄면(기본): 쉬운 한국어 화면 (화학물질 + 오염 수준만 선택). "
                 "켜면: 전체 연구용 인터페이스(영어) — 6가지 노출 모드, SMILES 구조 입력, 메커니즘 스위치, 모든 모델 파라미터.")

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
        neutral = None            # NEUTRAL-organic kwargs (None = PFAS path)
        spec = "Curated congener"
        E_m = -120
        fxy_source = "recommended"
        biomass = "oryza"
        uptake = api.DEFAULT_UPTAKE
        lipid_loading = False
        vmax_scale = None
        g_apo = None
        km_scale = 1.0
        compare = []
        preset_label = None
        preset_word = None
        use_custom_tables = False

        if not expert:
            # ----------------------------- SIMPLE sidebar (한국어) -----------------------------
            mode = "Model (parametric)"
            # A one-click policy scenario sets BOTH the chemical and the level, so a
            # presenter can flip the whole story instantly; "직접 설정" reveals the
            # manual chemical + level controls for hands-on exploration.
            st.header("① 시나리오")
            scen_label = st.radio(
                "빠른 시나리오", ["✏️ 직접 설정"] + list(_SCENARIOS_KO), index=1,
                help="정책 상황별로 화학물질과 오염 수준을 한 번에 설정합니다. "
                     "'직접 설정'을 고르면 물질·오염도를 각각 지정할 수 있습니다.")
            if scen_label == "✏️ 직접 설정":
                st.header("② 화학물질 선택")
                congener = st.selectbox("PFAS 화학물질", api.CONGENERS,
                                        index=api.CONGENERS.index("PFOA"),
                                        format_func=_cong_label_ko,
                                        help="특정 '영원한 화학물질' 하나. PFOA·PFOS가 가장 잘 알려져 있고, "
                                             "사슬이 길수록 대체로 식물에 더 잘 달라붙습니다.")
                st.header("③ 오염 정도")
                preset_label = st.radio("논의 오염 수준", list(_PRESETS_KO), index=1,
                                        help="토양수에 녹아 있는 PFAS의 양. 높을수록 식물로 더 많이 들어갑니다.")
                Cwo_const, preset_word = _PRESETS_KO[preset_label]
            else:
                congener, Cwo_const, preset_word, _scen_desc = _SCENARIOS_KO[scen_label]
                preset_label = scen_label            # carried into the summary handout
                st.success(f"**{scen_label}**\n\n{_scen_desc}")
                st.caption(f"→ 물질 **{_cong_label_ko(congener)}** · 오염 **{preset_word}** "
                           f"({Cwo_const:g} µg/L). '✏️ 직접 설정'에서 바꿀 수 있습니다.")
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

            st.header("2 · Compound")
            spec = st.radio("Compound class",
                            ["Curated congener", "SMILES (structure)", "Neutral organic (log Kow)"],
                            help="WHICH COMPOUND CLASS the run is about, not just which molecule — the "
                                 "first two are PFAS (permanently dissociated anion: GHK exclusion + "
                                 "carrier), the third is the Briggs/Kow neutral base (z=0, no exclusion, "
                                 "no carrier, nothing fitted). Everything downstream — map, dynamics, "
                                 "drivers, exposure modes, inverse, downloads — follows this choice.")
            if spec == "Curated congener":
                congener = st.selectbox("PFAS congener", api.CONGENERS,
                                        index=api.CONGENERS.index("PFOA"), format_func=_cong_label)
            elif spec == "Neutral organic (log Kow)":
                congener = None
                neutral = _neutral_panel()
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
            # E_m and f_xy are PFAS levers and are MEANINGLESS on the neutral path:
            # at z=0 the GHK term degenerates to passive diffusion (so a membrane
            # potential has nothing to act on) and f_xy is the computed Briggs TSCF,
            # not a fitted loading factor. Showing them would invite a reading the
            # model does not support, so they are absent rather than ignored.
            if not neutral:
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

            # PFAS mechanism switches (carrier/bypass, lipid loading) -- see _mechanism_panel
            if not neutral:
                uptake, lipid_loading, vmax_scale, g_apo, km_scale = _mechanism_panel()

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
                                           float(api.default_k_leach(
                                               congener, Koc=(neutral or {}).get("Koc"))), 0.0025,
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
                               "`phydrus` (requirements.txt). Locally the FORTRAN source is VENDORED in the "
                               "repo (no submodule to init): copy `external/hydrus_source/makefile` into "
                               "`source/`, run `make` there (gfortran), `pip install phydrus`. "
                               "Until built, the tool falls back to the parametric model.")
                elif neutral:
                    # A neutral organic has no congener, so the soil Kd comes from its
                    # own Koc (sidebar; Karickhoff default, editable) -- the PFAS
                    # chain-length QSPR would be meaningless here.
                    st.caption("Runs a **real HYDRUS-1D** paddy model (Richards + advection–dispersion + "
                               "linear Kd + root uptake) → Cwᵒ(t), Q_TP(t). Cached per setting.")
                    f_oc = st.slider("Soil organic carbon f_oc", 0.005, 0.05, 0.02, 0.005,
                                     help="Kd = Koc·f_oc → retardation R = 1+ρKd/θ. Koc is the "
                                          "value in the compound panel above.")
                    flood_until = st.slider("Flooded until  [day]", 30, 120, 90, 5)
                    percolation = st.slider("Percolation excess  [cm/day]", 0.0, 1.0, 0.30, 0.05)
                    season = 120.0
                    kd_n = float(neutral["Koc"]) * float(f_oc)
                    try:
                        drivers = _hydrus_drivers_cached(None, season, f_oc, flood_until, percolation,
                                                         biomass=biomass, Kd=kd_n)
                        st.success(f"HYDRUS-1D run complete — Kd = {kd_n:.3g} L/kg "
                                   f"(Koc {neutral['Koc']:.1f} × f_oc {f_oc:g}; mean-normalised).")
                    except Exception as e:                                  # noqa: BLE001
                        st.error(f"HYDRUS run failed: {e}")

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
                                    float(api.default_k_leach(
                                        congener, Koc=(neutral or {}).get("Koc"))), 0.0025,
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

    cfg.expert = expert
    cfg.mode = mode
    cfg.congener = congener
    cfg.drivers = drivers
    cfg.measured = measured
    cfg.Cwo_const = Cwo_const
    cfg.season = season
    cfg.soil_obj = soil_obj
    cfg.profile = profile
    cfg.measured_bio = measured_bio
    cfg.cwo_profile = cwo_profile
    cfg.cwo_kleach = cwo_kleach
    cfg.smiles = smiles
    cfg.neutral = neutral
    cfg.spec = spec
    cfg.E_m = E_m
    cfg.fxy_source = fxy_source
    cfg.biomass = biomass
    cfg.uptake = uptake
    cfg.lipid_loading = lipid_loading
    cfg.vmax_scale = vmax_scale
    cfg.g_apo = g_apo
    cfg.km_scale = km_scale
    cfg.compare = compare
    cfg.preset_label = preset_label
    cfg.preset_word = preset_word
    cfg.use_custom_tables = use_custom_tables
    return cfg
