"""Full research (English) Expert-mode view. Split out of app.py (HANDOFF P3-1)."""
import numpy as np
import streamlit as st

import model_api as api
import plots
import plot_svg

from ui.common import (_nearest_index, _simulate, _drivers_tuple, _yamazaki_observable,
                       _render_inverse_estimator, _glossary_md,
                       _png_bytes, _html_bytes)


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
    fxy_source = cfg.fxy_source
    biomass = cfg.biomass
    measured_bio = cfg.measured_bio
    lipid_loading = bool(getattr(cfg, "lipid_loading", api.DEFAULT_LIPID_LOADING))
    shoot = getattr(cfg, "shoot", api.DEFAULT_SHOOT)
    # ---------------------------------------------------------------- headline metrics
    # A non-default mechanism has to be visible next to the numbers it moved -- the
    # sidebar expander is collapsed, so the run would otherwise look like the shipped one.
    neutral = getattr(cfg, "neutral", None)
    _mech = ([] if lipid_loading == api.DEFAULT_LIPID_LOADING else
             [f"**lipid loading {'ON' if lipid_loading else 'OFF'}**"]) + \
            ([] if shoot == api.DEFAULT_SHOOT else [f"shoot: **{shoot}**"])
    if neutral:
        st.subheader(f"{congener}  (neutral organic · log Kow {p['log_kow']:.2f})  ·  source: {mode}")
    else:
        st.subheader(f"{congener}  (C{p['n_C']} {p['group']})  ·  source: {mode}")
    if _mech and not neutral:
        st.warning("Non-default mechanism — " + " · ".join(_mech) +
                   ". An alternative fitted variant (sidebar → ⚙️ Mechanism); every number "
                   "on this screen moves with it, and it is not the default model.")
    c1, c2, c3, c4 = st.columns(4)
    # DRY basis on both sides: the Yamazaki observations are dry weight, so the model's
    # harvest operator values (`baf_final_dw`) are shown, never the fresh-weight BAF
    dw = res.get("baf_final_dw")
    for col, tis in zip((c1, c2, c3), ("root", "straw", "grain")):
        pred = dw[tis] if dw else (res["straw_baf"] if tis == "straw" else res["baf_final"][tis])
        sub = f"obs {obs[tis]:.2f}" if (tis in obs and dw) else "no obs"
        if bio_baf and tis in bio_baf:
            sub = f"measured {bio_baf[tis]:.2f}"
        col.metric(f"{tis} BAF [L/kg {'dw' if dw else 'fw'}]", f"{pred:.2f}", sub,
                   delta_color="off")
    if neutral:
        # eᴺ is 1 by construction here, so it carries no information; TSCF is the
        # quantity that actually drives this run and it is COMPUTED, not fitted.
        c4.metric("TSCF", f"{res['TSCF']:.3f}", f"K_PW root {res['K_PW']['root']:.2f}",
                  delta_color="off")
    else:
        c4.metric("root entry r_eq", f"{p['r_eq']:.3g}", "fitted to Yamazaki's root",
                  delta_color="off")

    # ---- neutral scope panel (what this path is, and what it is not) --------
    if neutral:
        with st.expander("🧪 Neutral organic — what this run is, and its scope", expanded=True):
            st.markdown(
                "The framework's **Briggs/Kow base**, run on the *same* 4-compartment ODE with "
                "`z = 0`: no charge, so the membrane term is exactly passive diffusion and there "
                "is no carrier. "
                "**Nothing here is fitted** — `K_PW` and `TSCF` both follow from log Kow via "
                "published QSPRs, so this is the one setting where the DPU *backbone* is exposed "
                "without the fitted PFAS transport behind it.")
            _ns, _eq = api.NEUTRAL_APRIORI_OOS_SCORES, api.NEUTRAL_EQUILIBRIUM_SCORES
            _lo, _hi = (min(_ns.values()), max(_ns.values())) if _ns else (float("nan"),) * 2
            st.markdown(
                (f"- a-priori vs measured rice: **log10 RMSE {_ns['Liu 2023 rice root']:.3f}** root "
                 f"partition (Liu 2023) and **{_ns['Ge 2017 per-organ TF']:.3f}** per-organ TF "
                 f"(Ge 2017) — on the season-ODE basis this run uses; on the equilibrium basis "
                 f"proper to a root-partition measurement Liu is "
                 f"**{_eq.get('Liu 2023 rice root', float('nan')):.3f}**. Across the {len(_ns)} "
                 f"tables it has been scored on the error is log10 {_lo:.2f}–{_hi:.2f} "
                 f"(**×{10 ** _lo:.1f}–{10 ** _hi:.1f}**), with the ordering more reliable than "
                 f"the level (validation/neutral_scores.py).\n"
                 if {"Liu 2023 rice root", "Ge 2017 per-organ TF"} <= set(_ns) else
                 "- a-priori scores: run validation/neutral_scores.py.\n")
                + f"- `TSCF` = **{res['TSCF']:.3f}** ({p['tscf_model']} QSPR) · `K_PW` root "
                f"**{res['K_PW']['root']:.2f}** / stem {res['K_PW']['stem']:.2f} / leaf "
                f"{res['K_PW']['leaf']:.2f} L/kg · Briggs barley RCF at this log Kow "
                f"**{res['rcf_briggs']:.2f}** · lipid reading `{p['lipid_source']}`"
                + (f" · half-life **{res['half_life']:.3g} d**" if res.get("half_life")
                   else " · **no metabolism** (γ=0)"))
            if not res.get("half_life"):
                st.warning("γ = 0: the leaf is an **unbounded terminal accumulator**, so the leaf "
                           "and straw numbers are an UPPER BOUND, not a prediction. Set an "
                           "in-planta half-life in the sidebar.")
            st.caption("Scope: the **grain compartment is UNTESTED** for neutrals (no dataset "
                       "reports a neutral Kow series in grain under root-only exposure) — read it "
                       "as a structural output. Stem/leaf lipid are still Trapp 1994's soybean "
                       "values (the root is corroborated by measured cereal roots). This is a "
                       "RICE model: rice composition and rice crop drivers (by default the Tsukuba "
                       "2020 site-year) — another crop needs its own composition and drivers, not "
                       "just a log Kow.")
        if res.get("warning"):
            st.warning(res["warning"])
        if res.get("air_summary"):
            hl = res["air_summary"]["leaf"]["half_life"]
            st.caption(
                f"Air exchange ON — leaf volatilisation half-life "
                f"**{'∞' if not np.isfinite(hl) else f'{hl:.3g} d'}** (snapshot at final "
                f"mass/transpiration), leaf P_plant {res['air_summary']['leaf']['P_plant']:.3g} m/d. "
                "The flux scales with each tissue's specific SURFACE AREA, which this repo has only "
                "ever used as a leaf/grain *ratio* — treat absolute magnitudes as order-of-magnitude.")

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
                f"K_cw (root) = "
                f"{(p['K_cw'].get('root') if isinstance(p['K_cw'], dict) else p['K_cw']):.0f}  [L/kg]\n"
                f"- f_xy = {p['f_xy']:.4g} &nbsp;·&nbsp; L_Ph = {p['L_Ph']:.3g} &nbsp;·&nbsp; "
                f"r_eq = {p['r_eq']:.3g}\n"
                f"- B_k: " + ", ".join(f"{k} {v:.1f}" for k, v in res["B_k"].items()))
            for n in getattr(desc, "notes", []):
                st.caption("• " + n)

    is_biomon = (mode == "Biomonitoring (measured tissue)" and measured_bio)
    _biomon_note = ("⚠ This chart shows a **model reference run** (Cwᵒ=1), **not** your measured "
                    "biomonitoring data. Your measured values are on the **map** and **BAF** tabs.")

    # The PFAS-only tabs are DROPPED rather than shown empty for a neutral compound:
    # chain length, congener comparison, the Tang per-organ check and the Yamazaki bars
    # are all statements about the PFAS series and have no neutral counterpart.
    _PFAS_ONLY = ("📊 BAF vs observed", "🔗 Chain-length trends", "⚖️ Compare congeners",
                  "✅ Tang TF (OOS)")
    _names = ["🗺️ Plant & soil map", "📈 Tissue dynamics", "🟫 Soil & drivers",
              "📊 BAF vs observed", "🔗 Chain-length trends", "⚖️ Compare congeners",
              "✅ Tang TF (OOS)", "🔎 Inverse (Bayesian)", "ℹ️ About / coupling"]
    if neutral:
        _names = [n for n in _names if n not in _PFAS_ONLY]
    _t = dict(zip(_names, st.tabs(_names)))

    tabs = [_t.get(n) for n in
            ["🗺️ Plant & soil map", "📈 Tissue dynamics", "🟫 Soil & drivers",
             "📊 BAF vs observed", "🔗 Chain-length trends", "⚖️ Compare congeners",
             "✅ Tang TF (OOS)", "🔎 Inverse (Bayesian)", "ℹ️ About / coupling"]]
    # index-aligned: a dropped tab is None and its block is skipped below

    # ---- Tab 1: the plant + soil accumulation map ---------------------------
    with tabs[0]:
        # Expert reuses the SAME plant/soil SVG map as the Simple view (rendered
        # inline so the wide layout never clips it), just in English.
        if is_biomon:
            vals = dict(measured_bio["conc"])
            finite = [v for v in vals.values() if v is not None and np.isfinite(v)]
            cmax = max(finite) if finite else 1.0
            st.markdown(plot_svg.plant_svg(vals, cmin=0.0, cmax=cmax,
                                           cwo=measured_bio.get("Cwo"), lang="en"),
                        unsafe_allow_html=True)
            st.caption("Compartments coloured by the MEASURED concentration [µg/kg]; the tag shows "
                       "the pore-water Cwᵒ [µg/L]. Stem/leaf share the straw colour when only straw "
                       "is reported.")
        else:
            metric = st.radio("Colour by", ["concentration", "BAF"], horizontal=True)
            metric_key = "baf" if metric == "BAF" else "conc"
            day = st.slider("Day after transplant", float(res["t"][0]), float(res["t"][-1]),
                            float(res["t"][-1]), 1.0)
            ti = _nearest_index(res["t"], day)
            st.markdown(plot_svg.plant_svg_from_res(res, ti, lang="en", metric=metric_key),
                        unsafe_allow_html=True)
            st.caption("Each compartment is filled by its accumulation on a shared colour scale — drag "
                       "the day slider to watch where PFAS builds up (numbers in µg/kg, or L/kg for "
                       "BAF; the tag is pore-water µg/L). Leaf is xylem-terminal, grain is phloem-fed; "
                       "the root retains the anion (low f_xy).")

    # ---- Tab 2: tissue concentration dynamics -------------------------------
    with tabs[1]:
        if is_biomon:
            st.warning(_biomon_note)
        st.plotly_chart(plots.fig_tissue(res), width="stretch")
        st.plotly_chart(plots.fig_burden(res), width="stretch")
        if neutral:
            # the neutral partition IS K_PW (= W + L·a·Kow^b), the same B_k slot with
            # protein and cell-wall zeroed; f_xy is the computed TSCF, not a fitted one
            st.markdown("**K_PW [L/kg fw]** — "
                        + ", ".join(f"{k}: {v:.2f}" for k, v in res["K_PW"].items())
                        + f"  ·  TSCF={res['TSCF']:.4g} (computed, not fitted), "
                        f"κ_d={p['kappa_d']:.3g}, γ={p['gamma']:.4g} 1/d")
        else:
            st.markdown("**B_k [L/kg fw]** — " + ", ".join(f"{k}: {v:.2f}" for k, v in res["B_k"].items())
                        + f"  ·  f_xy={p['f_xy']:.4g}, L_Ph={p['L_Ph']:.4g}, root entry "
                        f"r_eq={p['r_eq']:.3g}")
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
                   "from the selected crop drivers (default: the Tsukuba 2020 site-year — the ORYZA2000 "
                   "port and its transpiration on JMA weather), reused across modes.")

    # PFAS-only tabs: dropped (not blanked) for a neutral compound -- the Yamazaki
    # bars, the chain-length series, the congener comparison and the Tang per-organ
    # check are all statements about the PFAS series with no neutral counterpart.
    if not neutral:
        # ---- Tab 4: BAF vs observed ---------------------------------------------
        with tabs[3]:
            if bio_baf:
                model_b = {"root": res["baf_final"]["root"], "straw": res["straw_baf"],
                           "grain": res["baf_final"]["grain"]}
                st.plotly_chart(plots.fig_biomon_compare(bio_baf, model_b), width="stretch")
                st.caption("Measured biomonitoring BAF (tissue conc ÷ pore water) vs the model prediction.")
            else:
                if lipid_loading:
                    st.info("**Lipid-facilitated loading is ON** (the default), so against these "
                            "Yamazaki bars the model is *in-sample* — its six lipid constants were fit "
                            "on them. Its out-of-sample record is on Tang, Kim and Li.")
                if obs and spec == "Curated congener":
                    # the model's version of the observation: dry basis, the two harvest runs of
                    # the 2020 site-year, the congener's irrigation-water level -- it tracks the
                    # transport/mechanism settings, not the sidebar exposure scenario
                    ykw = {k: v for k, v in sim_kw.items() if k not in ("biomass", "scenario")}
                    st.plotly_chart(plots.fig_baf(_yamazaki_observable(congener, **ykw), obs),
                                    width="stretch")
                    st.caption("Model bar = the Yamazaki 2020 site-year run the way the data were taken "
                               "(dry weight; straw with attached dead blades; two harvests combined as the "
                               "two cultivars are; exposure = the congener's irrigation-water level). It "
                               "follows the f_xy source and mechanism settings, not the sidebar "
                               "exposure scenario.")
                else:
                    st.plotly_chart(plots.fig_baf(res, obs), width="stretch")
                if not obs:
                    st.info("No Yamazaki BAF for this congener (model prediction only).")
                else:
                    with st.expander("ℹ️ What are the Yamazaki 2023 conditions? (and how to match them)"):
                        st.markdown(
                            "**Yamazaki et al. 2023**, *Environ. Sci. Technol.* **57** "
                            "([doi:10.1021/acs.est.2c08767](https://doi.org/10.1021/acs.est.2c08767)).\n\n"
                            "- **Design**: an open **paddy FIELD** (NARO, Tsukuba; 10 × 50 m, Andosol), continuously "
                            "flooded except the midsummer drainage, **ambient** exposure (irrigation water, rain, "
                            "airborne particles) — **nothing was spiked**; grown a **full cycle**, harvested ~30 d "
                            "after heading.\n"
                            "- **Cultivars**: the plotted observed BAFs are the **geomean of Kasalath (indica) + "
                            "Koshihikari (japonica)** (SI tables S16/S18/S19), one composite sample per tissue.\n"
                            "- **BAF definition**: the paper computes **no BAF**. This repo builds it as tissue conc "
                            "(**dry weight**, as reported) ÷ the **single irrigation-water grab** (no pore water was "
                            "measured), for **root / straw / grain**. Straw = unweighted mean of the four stem+leaf "
                            "height bands; grain = brown rice **reconstructed** as 0.88·white + 0.12·bran (brown rice "
                            "itself was not analysed); PFBS grain is wholly imputed and not shown. See "
                            "`docs/sampling_protocols_audit_KR.md`.\n\n"
                            "**The model bar is that same observation, computed**: the 2020 Tsukuba site-year (JMA "
                            "weather, the ORYZA2000 port on the paper's transplanting / heading / harvest dates), "
                            "two harvest runs combined as the two cultivars are, exposure = the congener's "
                            "irrigation-water level, dry weight with the dead blades still on the straw. It follows "
                            "the f_xy source and mechanism settings; the sidebar's exposure scenario does not "
                            "enter it. With f_xy source **fit** it reproduces the bars by construction (3 fitted "
                            "parameters per congener); otherwise the straw is a prediction.")

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
                c1, c3 = st.columns([2, 1])
                dose = ("low" if c1.radio("Tang dose", ["across-dose mean", "0.1 µg/g (env-closest)"],
                                          horizontal=True).startswith("0.1") else "mean")
                show_free = c3.checkbox(
                    "free anion (no lipid loading)", value=False,
                    help="Add the same shoot without the lipid loading. The grey bars are "
                         "the default model, whose fitted constants come from YAMAZAKI with nothing "
                         "fitted to Tang (a transfer, not a fit).")
                val = api.tang_tf_validation(congener, dose=dose)
                vall = (api.tang_tf_validation(congener, dose=dose, lipid_loading=False)
                        if show_free else None)
                st.plotly_chart(plots.fig_tang_tf(val, vall),
                                width="stretch")
                err = [np.log10(val["model_tf"][o] / val["tang_tf"][o]) for o in val["organs"]
                       if val["tang_tf"].get(o)]
                st.markdown(
                    f"- **{congener}, this dose**: log10 RMSE over the three organs "
                    f"**{float(np.sqrt(np.mean(np.square(err)))):.3f}** (computed now; nothing "
                    f"was fitted to Tang).\n"
                    f"- **Dry-weight basis** — Tang's TF is dry weight; the model's is formed from the "
                    f"harvest operator's dry-basis concentrations (`conc_dw`).")
            st.caption("Source: docs/literature_db/raw_si/tang2026_doseresponse.csv (SI S7/S8). Model: "
                       "the default configuration on the calibration site-year drivers (Tang's site and dates "
                       "are not in the repo). Scoring script: validation/oos_tang_lipid.py.")

    # ---- Tab 8: Inverse (Bayesian exposure estimate) ------------------------
    with tabs[7]:
        st.markdown("**Bayesian inverse** — estimate the pore-water exposure Cwᵒ from measured "
                    "tissue concentrations (Laplace posterior in log Cwᵒ; the well-posed "
                    "direction of `validation/bayesian_inverse_demo.py`). Transport is held at "
                    "the sidebar defaults.")
        if neutral:
            # Same inverse question, neutral forward model: the compound is a log Kow,
            # so `congener` is only a label and every neutral setting rides along.
            _render_inverse_estimator(
                congener, f_xy_source=fxy_source, biomass=biomass,
                key="inv_expert_neutral", simple=False, curated_only=False,
                mech=dict(log_kow=float(neutral["log_kow"]),
                          neutral_kw={k: v for k, v in neutral.items()
                                      if k not in ("Koc", "log_kow", "name")}))
        else:
            _render_inverse_estimator(
                congener, f_xy_source=fxy_source, biomass=biomass,
                key="inv_expert", simple=False,
                # run the inverse under the SAME mechanism as every forward tab
                mech={k: v for k, v in sim_kw.items()
                      if k not in ("f_xy_source", "biomass", "scenario")})
        st.caption("Identifiability: only the EXPOSURE level is estimated here. From tissue data "
                   "alone Q_TP·f_xy is a product ridge and Cwᵒ vs root-uptake conductance is "
                   "degenerate, so pinning transport absolutely needs an independent measurement "
                   "(xylem sap / pore-water probe). See docs + validation/bayesian_inverse_demo.py.")

    # ---- Tab 9: About / coupling --------------------------------------------
    with tabs[8]:
        st.markdown(
            """
### Six ways to drive the plant model

| Mode | PFAS driver Cwᵒ(t) comes from | Q_TP, M(t) | When to use |
|---|---|---|---|
| **Model (parametric)** | a constant you set (flat, or the analytic flooded shape) | measured FAO-56 / ORYZA | quick what-ifs, teaching |
| **Custom tables** | a Cwᵒ(t) table you type or paste | your own growth table, or measured | you have field numbers but no soil model |
| **HYDRUS / CSV drivers** | a HYDRUS-1D / Phydrus run (CSV) | from the same CSV, or measured | you have a calibrated soil-water-solute model |
| **Run HYDRUS-1D (live)** | a real HYDRUS-1D run executed here | HYDRUS root uptake + ORYZA | you want the engine to run in-app (needs it built) |
| **Soil inventory** | inverting a soil load with a Freundlich isotherm | measured | you know total soil PFAS but not pore water |
| **Biomonitoring** | a measured pore-water / soil-solution value | — (not needed) | you have field tissue + water concentrations |

**Live HYDRUS-1D** runs the genuine engine (built from `external/hydrus_source`) via `phydrus`:
a one-season paddy model (Richards flow + advection–dispersion + linear Kd sorption + root
water uptake) gives the **congener-dependent** pore water — weakly-sorbed short chains leach
under flooding (Cwᵒ falls), strongly-sorbed long chains stay buffered (flat). Kd comes from the
Koc(chain-length) QSPR. The FORTRAN source is **vendored in the repo** (no submodule to
initialise): copy `external/hydrus_source/makefile` into `source/`, run `make` there (needs
gfortran), and `pip install phydrus`; the mode auto-detects the engine.

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
