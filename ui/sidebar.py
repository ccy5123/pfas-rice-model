"""The sidebar (Simple vs Expert) -> a config namespace consumed by app.py.
Split out of app.py (HANDOFF P3-1); behaviour is unchanged."""
import os
from types import SimpleNamespace

import numpy as np
import streamlit as st

import model_api as api
import plots

from ui.common import (_EX, _cong_label, _cong_label_ko, _PRESETS_KO, _mol_svg,
                       _hydrus_drivers_cached, _hydrus_soil_congener)


def build():
    """Render the sidebar and return a SimpleNamespace with the scenario config."""
    cfg = SimpleNamespace()
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
    cfg.spec = spec
    cfg.E_m = E_m
    cfg.fxy_source = fxy_source
    cfg.biomass = biomass
    cfg.compare = compare
    cfg.preset_label = preset_label
    cfg.preset_word = preset_word
    cfg.use_custom_tables = use_custom_tables
    return cfg
