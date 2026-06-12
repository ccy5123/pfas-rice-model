"""
PFAS–Rice Uptake Model — interactive dashboard (Plotly + Streamlit)
===================================================================

Run locally:
    pip install -r requirements.txt -r requirements-app.txt
    streamlit run app.py

What this tool does
-------------------
Draws the soil + 4-compartment rice plant to scale and colours each compartment
by its PFAS accumulation (a heat colormap you can scrub through the season), and
plots the supporting time series.  It covers the whole input space:

  1. Model (parametric)        — built-in measured/illustrative forcings.
  2. HYDRUS / CSV drivers      — drop in a HYDRUS-1D/Phydrus run (Cwᵒ, Q_TP, M).
  3. Run HYDRUS-1D (live)      — execute the real HYDRUS engine here (if built).
  4. Soil inventory            — invert a soil PFAS load to pore water (Freundlich).
  5. Biomonitoring             — measured tissue concentrations; HYDRUS not needed.

The compound can be one of the 13 curated congeners OR **any PFAS by SMILES
structure** (RDKit → descriptors → read-across/QSPR → Compound; src/pfas_structure.py).

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

st.set_page_config(page_title="PFAS–Rice Uptake Model", layout="wide")
st.title("🌾 PFAS–Rice Compartmental Uptake Model")
st.caption("Mechanistic 4-compartment dynamic model for permanently-anionic PFAS in paddy rice "
           "(IOC extension of the Trapp/Brunetti DPU). Parameters are measured/cited "
           "(docs/literature_db). Charts are interactive — hover, zoom, toggle. Outputs illustrative.")


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
    return api.simulate(congener, **kw)


def _drivers_tuple(d):
    return (tuple(d["t"]), tuple(d["Cwo"]), tuple(d["Qtp"]),
            tuple(np.asarray(d["M"]).ravel()), int(np.asarray(d["M"]).shape[1]))


@st.cache_data(show_spinner="Parameterising structure (RDKit)…")
def _simulate_smiles(smiles, **kw):
    """Cache a SMILES (structure) run: RDKit → descriptors → Compound → full ODE."""
    drv = kw.pop("drivers_tuple", None)
    if drv is not None:
        t, Cwo, Qtp, Mflat, ncol = drv
        kw["drivers"] = dict(t=np.array(t), Cwo=np.array(Cwo), Qtp=np.array(Qtp),
                             M=np.array(Mflat).reshape(-1, ncol))
    return api.simulate_from_smiles(smiles, **kw)


@st.cache_data(show_spinner=False)
def _mol_png(smiles, w=330, h=190):
    """Render a SMILES to a 2-D structure PNG (RDKit); None if it can't be parsed."""
    try:
        from rdkit import Chem
        from rdkit.Chem.Draw import rdMolDraw2D
        m = Chem.MolFromSmiles(smiles)
        if m is None:
            return None
        d = rdMolDraw2D.MolDraw2DCairo(w, h)
        d.drawOptions().padding = 0.12
        d.DrawMolecule(m)
        d.FinishDrawing()
        return d.GetDrawingText()
    except Exception:                                       # noqa: BLE001
        return None


@st.cache_data(show_spinner="Running HYDRUS-1D…")
def _hydrus_drivers_cached(congener, season, f_oc, flood_until, percolation):
    """Cache a real HYDRUS-1D paddy run (a few seconds) per parameter set."""
    drv, _ = api.hydrus_drivers(congener, season=season, f_oc=f_oc,
                                flood_until=float(flood_until), percolation=float(percolation))
    return drv


# ---------------------------------------------------------------- sidebar
with st.sidebar:
    st.header("1 · Data source")
    mode = st.radio(
        "How is the soil exposure supplied?",
        ["Model (parametric)", "HYDRUS / CSV drivers", "Run HYDRUS-1D (live)",
         "Soil inventory → pore water", "Biomonitoring (measured tissue)"],
        help="Five ways to feed the plant model. 'Run HYDRUS-1D (live)' executes the "
             "real engine (if built); biomonitoring needs no soil model.")

    st.header("2 · PFAS compound")
    spec = st.radio("Specify by", ["Curated congener", "SMILES (structure)"], horizontal=True,
                    help="Pick one of the 13 curated congeners, or paste ANY PFAS structure (SMILES) "
                         "to parameterise it from chemistry (RDKit read-across / QSPR).")
    smiles = None
    if spec == "Curated congener":
        congener = st.selectbox("PFAS congener", api.CONGENERS, index=api.CONGENERS.index("PFOA"))
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
            png = _mol_png(smiles)
            if png is not None:
                st.image(png, caption="structure (RDKit)", use_container_width=True)
            else:
                st.caption("⚠ could not parse this SMILES into a structure")
    E_m = st.slider("Root membrane potential E_m  [mV]", -160, -90, -120, 5,
                    help="GHK anion-exclusion lever (rice −116…−140 mV; NH₄⁺ depolarises).")
    fxy_label = st.radio("Root→shoot loading f_xy",
                         ["recommended (monotone, physical)", "W2 fit (reproduces Yamazaki)"])
    fxy_source = "recommended" if fxy_label.startswith("recommended") else "W2fit"

    st.header("3 · Scenario")
    drivers = None
    measured = True
    Cwo_const = 1.0
    season = 120.0
    soil_obj = None
    profile = None
    measured_bio = None

    if mode == "Model (parametric)":
        Cwo_const = st.number_input("Pore-water Cwᵒ  [µg/L]", min_value=0.0, value=1.0, step=0.1,
                                    help="Free anion conc. driving root uptake. 1.0 → tissue conc equals BAF.")
        measured = st.checkbox("Measured forcings (Q_TP, M_s)", value=True,
                               help="On: transpiration from Kumari/NayHtoon, biomass from ORYZA IR72.")
        season = float(st.slider("Season length  [days]", 90, 160, 120, 5))

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
                drivers = api.drivers_from_arrays(t, Cwo, Qtp=Q, M=M)
                st.success(f"Loaded {len(t)} rows from {up.name}.")
            elif use_ex:
                drivers = api.load_driver_csv(os.path.join(_EX, "hydrus_drivers_example.csv"))
                st.info("Using examples/hydrus_drivers_example.csv (synthetic HYDRUS-style run).")
        except Exception as e:                                  # noqa: BLE001
            st.error(f"Could not read drivers: {e}")

    elif mode == "Run HYDRUS-1D (live)":
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
        else:
            st.caption("Runs a **real HYDRUS-1D** paddy model (Richards + advection–dispersion + "
                       "linear Kd + root uptake) for this congener → Cwᵒ(t), Q_TP(t). Cached per setting.")
            f_oc = st.slider("Soil organic carbon f_oc", 0.005, 0.05, 0.02, 0.005,
                             help="Kd = Koc(chain length)·f_oc → per-congener retardation R = 1+ρKd/θ.")
            flood_until = st.slider("Flooded until  [day]", 30, 120, 90, 5,
                                    help="Continuous flooding (clean irrigation) until this day, then drainage.")
            percolation = st.slider("Percolation excess  [cm/day]", 0.0, 1.0, 0.30, 0.05,
                                    help="Clean-water through-flow that leaches the dissolved pool.")
            season = 120.0
            try:
                drivers = _hydrus_drivers_cached(congener, season, f_oc, flood_until, percolation)
                st.success(f"HYDRUS-1D run complete — Cwᵒ(t) for {congener} "
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
        k_leach = st.slider("Leaching rate k_leach  [1/day]", 0.0, 0.1, 0.02, 0.005) if flood else 0.0
        Cwo, soil_obj = api.pore_water_from_inventory(
            t, C_total, K_F=K_F, n=n_F, theta_g=theta_g, flooded=flooded, k_leach=k_leach)
        drivers = api.drivers_from_arrays(t, Cwo, season=season)

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


# ---------------------------------------------------------------- run the model
sim_kw = dict(E_m_mV=E_m, f_xy_source=fxy_source)
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
                                   measured_forcing=measured, **sim_kw)
    except Exception as e:                              # noqa: BLE001
        st.error(f"Could not build a compound from that SMILES — check the structure.\n\n`{e}`")
        st.stop()
    congener = res["congener"]
    desc = res.get("descriptors")
    provisional = bool(res.get("provisional", False))
elif drivers is not None:
    res = _simulate(congener, drivers_tuple=_drivers_tuple(drivers), **sim_kw)
else:
    res = _simulate(congener, Cwo=Cwo_const, season=season, measured_forcing=measured, **sim_kw)
obs = api.observed_baf(congener)
p = res["params"]

# biomonitoring-derived BAFs (measured side)
bio_baf = None
if measured_bio and measured_bio.get("Cwo"):
    bio_baf = api.baf_from_measurement(measured_bio["conc"], measured_bio["Cwo"])

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

# ---- structure → parameters panel (SMILES mode) -----------------------------
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

tabs = st.tabs(["🗺️ Plant & soil map", "📈 Tissue dynamics", "🟫 Soil & drivers",
                "📊 BAF vs observed", "🔗 Chain-length trends", "⚖️ Compare congeners",
                "ℹ️ About / coupling"])

# ---- Tab 1: the plant + soil accumulation map -------------------------------
with tabs[0]:
    if mode == "Biomonitoring (measured tissue)" and measured_bio:
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

# ---- Tab 2: tissue concentration dynamics -----------------------------------
with tabs[1]:
    st.plotly_chart(plots.fig_tissue(res), width="stretch")
    st.markdown(f"**B_k [L/kg fw]** — " + ", ".join(f"{k}: {v:.2f}" for k, v in res["B_k"].items())
                + f"  ·  f_xy={p['f_xy']:.4g}, L_Ph={p['L_Ph']:.4g}, κ_d={p['kappa_d']:.3g}")

# ---- Tab 3: soil & drivers --------------------------------------------------
with tabs[2]:
    st.plotly_chart(plots.fig_drivers(res), width="stretch")
    if soil_obj is not None:
        st.plotly_chart(plots.fig_isotherm(soil_obj, Cwo_now=float(res["Cwo"][-1])),
                        width="stretch")
    st.plotly_chart(plots.fig_soil_profile(res, profile=profile), width="stretch")
    st.caption("Cwᵒ(t) is the only PFAS-specific driver; Q_TP(t) and M(t) are crop physiology "
               "(measured FAO-56 transpiration + ORYZA IR72 biomass) reused across modes.")

# ---- Tab 4: BAF vs observed -------------------------------------------------
with tabs[3]:
    if bio_baf:
        model_b = {"root": res["baf_final"]["root"], "straw": res["straw_baf"],
                   "grain": res["baf_final"]["grain"]}
        st.plotly_chart(plots.fig_biomon_compare(bio_baf, model_b), width="stretch")
        st.caption("Measured biomonitoring BAF (tissue conc ÷ pore water) vs the model prediction.")
    else:
        st.plotly_chart(plots.fig_baf(res, obs), width="stretch")
        if not obs:
            st.info("No Yamazaki BAF for this congener (model prediction only).")

# ---- Tab 5: chain-length trends ---------------------------------------------
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

# ---- Tab 6: compare congeners -----------------------------------------------
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

# ---- Tab 7: about / coupling ------------------------------------------------
with tabs[6]:
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
        """)
    st.caption("Open item (recorded, not modelled): long-chain PFCA shoot accumulation "
               "(hysteretic sorption) — see docs/nstem_gradient_exploration.md.")
