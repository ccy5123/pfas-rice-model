"""Full research (English) Expert-mode view. Split out of app.py (HANDOFF P3-1)."""
import numpy as np
import streamlit as st

import model_api as api
import plots

from ui.common import (_nearest_index, _simulate, _drivers_tuple, _simulate_twopool_seq,
                       _render_inverse_estimator, _glossary_md, _png_bytes, _html_bytes)


def render(cfg):
    """Render the Expert (full research, English) view from a populated cfg."""
    congener = cfg.congener
    res = cfg.res
    obs = cfg.obs
    bio_baf = cfg.bio_baf
    p = cfg.p
    mode = cfg.mode
    spec = cfg.spec
    smiles = cfg.smiles
    desc = cfg.desc
    provisional = cfg.provisional
    soil_obj = cfg.soil_obj
    profile = cfg.profile
    compare = cfg.compare
    drivers = cfg.drivers
    sim_kw = cfg.sim_kw
    Cwo_const = cfg.Cwo_const
    season = cfg.season
    measured = cfg.measured
    E_m = cfg.E_m
    fxy_source = cfg.fxy_source
    biomass = cfg.biomass
    measured_bio = cfg.measured_bio
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
        st.markdown("**B_k [L/kg fw]** — " + ", ".join(f"{k}: {v:.2f}" for k, v in res["B_k"].items())
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
        _map_fig = plots.fig_schematic_from_res(res, "conc", -1)
        png, why = _png_bytes(_map_fig)
        if png is not None:
            st.download_button("Plant map (PNG)", png, file_name=f"{congener}_map.png", mime="image/png")
        else:
            # No kaleido/Chrome -> offer an interactive HTML instead (always works).
            html, _ = _html_bytes(_map_fig)
            if html is not None:
                st.download_button("Plant map (interactive HTML)", html,
                                   file_name=f"{congener}_map.html", mime="text/html")
            st.caption("Static **PNG** export needs the optional `kaleido` package — "
                       "`pip install kaleido && plotly_get_chrome`, then re-run. Meanwhile the "
                       "**interactive HTML** above (zoom/tooltips in a browser) and the CSVs work without it.")
