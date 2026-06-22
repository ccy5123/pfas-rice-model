"""
High-level model API for the PFAS-rice uptake model
===================================================

A thin, UI-agnostic wrapper that loads the canonical measured parameters
(`params/parameters.json`) and runs the 4-compartment ODE for a chosen congener
and scenario.  Used by the Streamlit app (`app.py`), the tests, and any script.

Everything here is plain functions returning plain dicts/arrays -- no plotting,
no Streamlit -- so it can be exercised head-less.
"""
from __future__ import annotations
import json, os
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)

from pfas_rice_plant_module_4pool_surf import (
    Environment, Compound, Compartment, RiceUptakeModel, PlantInputs,
    binding_factors, root_uptake, _logistic, ROOT, STEM, LEAF, FRUIT)
import forcing_rice as fr
import growth_rice as gr
from soil_paddy import FreundlichSoil
from soil_paddy_redox_corrected import PaddyRedoxCorrected

with open(os.path.join(_ROOT, "params", "parameters.json")) as _f:
    PARAMS = json.load(_f)
_CARR = PARAMS["carrier_MichaelisMenten"]
_COMP = PARAMS["tissue_composition_recommended"]
_CONG = {c["name"]: c for c in PARAMS["congeners"]}

CONGENERS = [c["name"] for c in PARAMS["congeners"]]            # 12, ordered
TISSUES = ("root", "stem", "leaf", "grain")
# tissue keys that compose "straw" (the bulk shoot reported in agronomy)
STRAW_PARTS = ("stem", "leaf")

# --- EXPLORATORY lipid-facilitated loading model (opt-in; default OFF) ----------
# Two superposed transport mechanisms (docs/fxy_longchain_lipid_exploration.md):
#   free anion TSCF (declining) loads the FREE conc; a K_PL-GATED, B-independent
#   "bound" term loads the membrane/lipid-associated pool into xylem (g_xy) and
#   phloem (g_ph) so high-binding long chains reach the shoot. Global params fit to
#   Yamazaki (excl. near-MQL PFDoDA); IN-SAMPLE shape/mechanism, NOT validated.
LIPID_LOADING = dict(a=0.25, b=0.50, g_xy=0.05, g_ph=0.010, K_half=3000.0, pfsa_ln=1.25)


def lipid_loading_conductances(n_C, K_PL, group="PFCA"):
    """(f_xy_free, g_xy, g_ph) for the K_PL-gated lipid-loading model.

    f_xy_free = a*exp(-b*(n-4))             free-anion TSCF (declining)
    phi       = K_PL/(K_PL+K_half)          lipid takeover (off for short chains)
    g_xy,g_ph = g_*max * phi                B-independent bound loading
    PFSA carries a head-group factor exp(-pfsa_ln). See LIPID_LOADING.
    """
    p = LIPID_LOADING
    sf = np.exp(-p["pfsa_ln"]) if group == "PFSA" else 1.0
    phi = K_PL / (K_PL + p["K_half"])
    return (p["a"] * np.exp(-p["b"] * (n_C - 4)) * sf, p["g_xy"] * phi * sf, p["g_ph"] * phi * sf)


def _compartments():
    g = _COMP
    return [Compartment("root",  g["root"]["theta_fw"], g["root"]["f_prot"], g["root"]["f_PL"], g["root"]["f_cw"]),
            Compartment("stem",  g["stem"]["theta_fw"], g["stem"]["f_prot"], g["stem"]["f_PL"], g["stem"]["f_cw"]),
            Compartment("leaf",  g["leaf"]["theta_fw"], g["leaf"]["f_prot"], g["leaf"]["f_PL"], g["leaf"]["f_cw"], S=20.0),
            Compartment("grain", g["grain_brown"]["theta_fw"], g["grain_brown"]["f_prot"],
                        g["grain_brown"]["f_PL"], g["grain_brown"]["f_cw"], S=2.0)]


def observed_baf(congener):
    """Yamazaki 2023 root/straw/grain BAF for a congener, or {} if not measured."""
    out = {}
    path = os.path.join(_ROOT, "data_obs", "obs_baf_Yamazaki.csv")
    import csv
    with open(path) as f:
        for r in csv.DictReader(f):
            if r["compound"] == congener:
                out[r["tissue"]] = float(r["baf"])
    return out


def simulate_twopool_carrier(congener, *, f_xy, vmax_in, g_xy=0.0, g_ph=0.0, k_off=0.02,
                             seq=1.0, biomass="oryza", season=120.0):
    """Opt-in CARRIER two-pool root model (the long-chain breakthrough; src/pfas_rice_two_pool.py).

    Splits the root into a mobile pool (feeds the xylem / sets the uptake gradient) and a
    slow bound store (holds the measured root burden), so a LOW f_xy (strong root retention)
    and an ENHANCED active carrier (high uptake) are independent levers -- which closes the
    long chains the single-pool core cannot. Returns root/straw/grain BAF (Cwo=1). The
    canonical ``simulate`` (4pool_surf) is unchanged; this is additive.

    NOTE: distinct from `simulate_twopool_seq` -- that one's second pool is an irreversible
    SEQUESTRATION sink (non-K_PL U-shaped k_seq, Yamazaki-fit); this one's is a reversible
    bound STORE tuned by the carrier/f_xy levers (saturated long-chain closure)."""
    import pfas_rice_two_pool as tp
    return tp.simulate(congener, f_xy, vmax_in, g_xy=g_xy, g_ph=g_ph, k_off=k_off,
                       seq=seq, biomass=biomass, season=season)


def close_longchain_2pool(congener, *, obs=None, biomass="oryza", season=120.0):
    """Saturated 3-param structural-adequacy fit of the CARRIER two-pool model
    (`simulate_twopool_carrier`) to the measured BAFs (free f_xy -> straw, active carrier
    -> root, g_ph -> grain). ``obs`` defaults to the Yamazaki BAF for the congener. Returns
    the fitted (f_xy, vmax_in, g_ph) + the reproduced root/straw/grain. Reproduction (DOF 0),
    NOT a-priori prediction."""
    import pfas_rice_two_pool as tp
    obs = obs or observed_baf(congener)
    if not {"root", "straw", "grain"} <= set(obs):
        raise KeyError(f"need root/straw/grain BAF for {congener!r}; got {sorted(obs)}")
    return tp.close_longchain(congener, obs, biomass=biomass, season=season)


# --- Tang 2026 per-organ TF validation (dry-weight basis) -----------------------
TANG_CONGENERS = ("PFOA", "PFOS", "GenX")
# f_xy re-fit to Tang TF at the 0.1 ug/g dose (OVERRIDE-only; validation/tang2026_fxy_refit.py).
# PFOS is dataset-dependent (Yamazaki W2 0.142 vs Tang ~0.32) -> a condition, not one value.
TANG_REFIT_FXY = {"PFOA": 0.097, "PFOS": 0.320, "GenX": 0.017}
_TANG_TF_ENDPOINT = {"TF_stalk": "stalk", "TF_leaf": "leaf", "TF_endosperm": "endosperm"}
# (model grain pool, Tang organ name, theta key) for the fresh->dry conversion
_TANG_ORGANS = (("stem", "stalk", "stem"), ("leaf", "leaf", "leaf"), ("grain", "endosperm", "grain_brown"))


def tang_observed_tf(congener, dose="mean"):
    """Tang 2026 MEASURED per-organ TF (DRY weight) -> {stalk,leaf,endosperm: TF}, or {}.

    dose='mean' averages the 5 soil doses; 'low' uses 0.1 ug/g (environmentally closest,
    least toxicity-confounded). Source: raw_si/tang2026_doseresponse.csv (SI Table S8).
    """
    if congener not in TANG_CONGENERS:
        return {}
    import csv
    path = os.path.join(_ROOT, "docs", "literature_db", "raw_si", "tang2026_doseresponse.csv")
    by = {}
    with open(path) as f:
        for r in csv.DictReader(x for x in f if not x.lstrip().startswith("#")):
            org = _TANG_TF_ENDPOINT.get(r["endpoint"])
            if org and r["compound"] == congener:
                by.setdefault(org, {})[float(r["dose_ugg"])] = float(r["value"])
    if dose == "low":
        return {o: v[min(v)] for o, v in by.items()}
    return {o: float(np.mean(list(v.values()))) for o, v in by.items()}


def tang_tf_validation(congener, f_xy_source="recommended", use_refit=False,
                       dose="mean", biomass="oryza", season=150.0, lipid_loading=False):
    """Model per-organ transfer factor (DRY-weight) vs Tang 2026, for a Tang congener.

    The model conc is fresh-weight (C=B_k*Cw, basis A); Tang TF is dry/dry, and the
    (1-theta_fw) factor differs by tissue (root 0.90 vs grain 0.14) so it does NOT cancel
    in C_tissue/C_root: TF_dw = TF_fw*(1-theta_root)/(1-theta_tissue). Uses the
    redistributed-shoot model (`simulate_nstem_leaf`) for a sensible stem~leaf split;
    biomass='oryza' drives it with the mechanistic ORYZA2000 biomass. use_refit applies
    the Tang-calibrated f_xy override. ``lipid_loading`` turns on the K_PL-gated
    lipid-facilitated loading mechanism (constants fit on Yamazaki, NOT Tang -> applying it
    here is out-of-sample; see validation/oos_tang_lipid.py). Returns None for non-Tang congeners.
    NOTE: grain/endosperm is structurally under-predicted ~3-8x (docs/tang2026_grain_units_exploration.md).
    """
    if congener not in TANG_CONGENERS:
        return None
    fxy = TANG_REFIT_FXY[congener] if use_refit else None
    r = simulate_nstem_leaf(congener, Cwo=1.0, season=season, biomass_fn=_biomass_fn(biomass),
                            f_xy_source=f_xy_source, f_xy_override=fxy, lipid_loading=lipid_loading)
    froot = 1.0 - _COMP["root"]["theta_fw"]
    model = {org: r["tf_final"][mk] * froot / (1.0 - _COMP[tk]["theta_fw"])
             for mk, org, tk in _TANG_ORGANS}
    return dict(congener=congener, organs=[o for _, o, _ in _TANG_ORGANS],
                model_tf=model, tang_tf=tang_observed_tf(congener, dose),
                f_xy=float(r["params"]["f_xy"]), use_refit=bool(use_refit), dose=dose,
                biomass=biomass, refit_fxy=TANG_REFIT_FXY[congener])


def chain_table():
    """Per-congener summary parameters (for the chain-length overview)."""
    rows = []
    for c in PARAMS["congeners"]:
        rows.append(dict(name=c["name"], n_C=c["n_C"], group=c["group"],
                         K_PL=c["K_PL_Lkg"], K_prot=c["K_prot_Lkg"],
                         K_cw_root=c["K_cw_wholecw_Lkg"]["root"],
                         f_xy_recommended=c["f_xy_recommended"],
                         f_xy_W2fit=c.get("f_xy_W2fit"),
                         B_root=c["B_k_basisA_Lkg_fw"]["root"],
                         B_grain=c["B_k_basisA_Lkg_fw"]["grain"]))
    return rows


def _biomass_fn(biomass="oryza"):
    """Organ-biomass callable (t, season) -> dict{root,stem,leaf,grain} [kg/hill].

    'oryza'       -> the mechanistic ORYZA2000 Level-1 carbon balance (`oryza_growth`;
                     radiation/temperature-driven). THE DEFAULT (first-principles).
    'growth_rice' -> ORYZA IR72 DVS-partitioning on a logistic total-biomass curve
                     (`growth_rice`; the lightweight reconstruction). NOTE: the per-congener
                     f_xy_W2fit/L_Ph_W2fit and `reproduce_demo.py` reproduction were tuned on
                     a placeholder/growth_rice driver, so use 'growth_rice' to match those
                     legacy artifacts; the live default is the mechanistic ORYZA2000.
    """
    if biomass == "oryza":
        import oryza_growth as og
        return lambda t, s: og.organ_biomass_oryza(t, p=og.OryzaParams(season=s))
    return gr.organ_biomass


def _default_drivers(t, season, Cwo, measured_forcing, biomass="oryza"):
    """Build the (Cwo, Qtp, M) driver arrays for the built-in scenarios."""
    if measured_forcing:
        Qtp = fr.Q_TP(t, season)
        b = _biomass_fn(biomass)(t, season)
        M = np.maximum(np.column_stack([b["root"], b["stem"], b["leaf"], b["grain"]]), 1e-4)
        leaf_loss = b.get("leaf_death_rate")            # only the senescing ORYZA driver supplies it
    else:
        Qtp = 0.05 + 0.35 * np.exp(-((t - 0.62 * season) ** 2) / (2 * (season / 4.8) ** 2))
        M = np.column_stack([_logistic(t, 1e-3, 0.030, 0.10, season * 0.17),
                             _logistic(t, 1e-3, 0.040, 0.10, season * 0.21),
                             _logistic(t, 1e-3, 0.050, 0.12, season * 0.25),
                             _logistic(t, 1e-5, 0.025, 0.18, season * 0.67)])
        leaf_loss = None
    Cwo_series = np.full_like(t, float(Cwo))
    return Cwo_series, Qtp, M, leaf_loss


def simulate(congener="PFOA", Cwo=1.0, E_m_mV=-120.0, f_xy_source="recommended",
             f_xy_override=None, L_Ph_override=None, kappa_d_override=None,
             lipid_loading=False, g_xy_override=None, g_ph_override=None,
             season=120.0, n_t=241, measured_forcing=True, biomass="oryza",
             drivers=None, K_surf=0.0, record=None):
    """Run the 4-compartment ODE for one congener and scenario.

    Parameters
    ----------
    congener : one of CONGENERS.
    Cwo : constant pore-water free concentration C_w^o [ug/L] (1.0 -> conc == BAF).
    E_m_mV : root membrane potential [mV] (GHK anion-exclusion lever).
    f_xy_source : "recommended" (monotone) or "W2fit" (reproduces Yamazaki).
    f_xy_override : if given, use this f_xy instead.
    L_Ph_override, kappa_d_override : if given, use these phloem-loading / root
        membrane-conductance values instead of the per-congener W2 fits.
    measured_forcing : True -> Q_TP from forcing_rice, organ biomass M from the
        `biomass` driver; False -> the illustrative logistic placeholders.
    biomass : 'oryza' (the mechanistic ORYZA2000 Level-1 carbon balance,
        `oryza_growth`; THE DEFAULT, first-principles) or 'growth_rice' (ORYZA IR72
        partitioning on a logistic; the lightweight reconstruction that matches the
        legacy f_xy_W2fit/`reproduce_demo` calibration). Ignored when `drivers` supply M.
    drivers : optional dict with arrays {'t','Cwo','Qtp','M'} (M shape (n,4)) that
        OVERRIDE the built-in forcings -- the entry point for a HYDRUS-1D / Phydrus
        run or a soil-inventory inversion (see `drivers_from_arrays`,
        `pore_water_from_inventory`, `load_driver_csv`). When given, season/n_t/
        measured_forcing and the scalar Cwo are ignored for the driver series.
    K_surf : root surface/plaque sorption [L/kg] added to the *measured* root BAF
        only (dead-end pool; leaves the ODE untouched). 0 for low-carbon soils.

    Returns a dict with t, per-compartment conc & BAF time series, finals, straw,
    the driver series actually used (Cwo, Qtp, M), B_k, and the effective params.
    """
    if record is not None:
        c = record                       # custom (e.g. SMILES-derived) congener record
    elif congener not in _CONG:
        raise KeyError(f"unknown congener {congener!r}; known: {CONGENERS}")
    else:
        c = _CONG[congener]

    if drivers is not None:
        t = np.asarray(drivers["t"], dtype=float)
        Cwo_series = np.asarray(drivers["Cwo"], dtype=float)
        Qtp = np.asarray(drivers["Qtp"], dtype=float)
        M = np.asarray(drivers["M"], dtype=float)
        leaf_loss = drivers.get("leaf_loss")            # e.g. ORYZA senescence rate
        season = float(t[-1])
    else:
        t = np.linspace(0.0, season, n_t)
        Cwo_series, Qtp, M, leaf_loss = _default_drivers(t, season, Cwo, measured_forcing, biomass)
    inputs = PlantInputs(t=t, Cwo=Cwo_series, Qtp=Qtp, M=M, leaf_loss=leaf_loss)

    f_xy_def, L_Ph_def, kappa_d_def, g_xy_def, g_ph_def = _transport_defaults(c, f_xy_source, lipid_loading)
    f_xy = float(f_xy_override) if f_xy_override is not None else float(f_xy_def)
    g_xy = float(g_xy_override) if g_xy_override is not None else float(g_xy_def)
    g_ph = float(g_ph_override) if g_ph_override is not None else float(g_ph_def)
    kappa_d = float(kappa_d_override) if kappa_d_override is not None else float(kappa_d_def)
    L_Ph = float(L_Ph_override) if L_Ph_override is not None else float(L_Ph_def)
    cmpd = Compound(name=congener, K_prot=c["K_prot_Lkg"], K_PL=c["K_PL_Lkg"],
                    K_cw=c["K_cw_wholecw_Lkg"]["root"], kappa_d=kappa_d,
                    Vmax_in=_CARR["Vmax_in"], Km_in=_CARR["Km_in"],
                    Vmax_out=_CARR["Vmax_out"], Km_out=_CARR["Km_out"],
                    L_Ph=L_Ph, f_xy=f_xy, g_xy=g_xy, g_ph=g_ph, K_surf=float(K_surf))
    comps = _compartments()
    env = Environment(E=E_m_mV / 1000.0)
    model = RiceUptakeModel(env=env, cmpd=cmpd, comps=comps, inputs=inputs)
    sol = model.solve(t)
    C = sol.y                                            # (4, n_t)
    Mf = inputs.M_(t[-1])
    straw = (C[STEM] * Mf[STEM] + C[LEAF] * Mf[LEAF]) / (Mf[STEM] + Mf[LEAF])
    B = binding_factors(comps, cmpd)
    # time-resolved BAF normaliser (per-step pore water); guard zeros
    cw = np.where(Cwo_series > 0, Cwo_series, np.nan)
    cwo_ref = float(Cwo_series[-1]) if Cwo_series[-1] > 0 else (
        float(np.nanmax(Cwo_series)) if np.nanmax(Cwo_series) > 0 else 1.0)
    baf_series = {k: C[i] / cw for i, k in enumerate(TISSUES)}
    root_baf_total = float(C[ROOT, -1] / cwo_ref + cmpd.K_surf)
    return dict(
        t=t, congener=congener, success=bool(sol.success), season=season,
        conc={k: C[i] for i, k in enumerate(TISSUES)},
        straw=straw,
        Cwo=Cwo_series, Qtp=Qtp, M=M,                    # the drivers actually used
        baf=baf_series,
        baf_final={k: float(C[i, -1] / cwo_ref) for i, k in enumerate(TISSUES)},
        straw_baf=float(straw[-1] / cwo_ref),
        root_baf_total=root_baf_total,
        cwo_ref=cwo_ref,
        B_k={k: float(B[i]) for i, k in enumerate(TISSUES)},
        N=float(env.N), eN=float(np.exp(env.N)),
        params=dict(f_xy=f_xy, L_Ph=L_Ph, kappa_d=kappa_d, g_xy=g_xy, g_ph=g_ph,
                    K_PL=c["K_PL_Lkg"], K_prot=c["K_prot_Lkg"],
                    K_cw=c["K_cw_wholecw_Lkg"]["root"], K_surf=float(K_surf),
                    n_C=c["n_C"], group=c["group"]),
    )


def apportionment(congener="PFOA", Cwo=1.0, E_m_mV=-120.0, f_xy_source="recommended",
                  f_xy_override=None, L_Ph_override=None, kappa_d_override=None,
                  lipid_loading=False, g_xy_override=None, g_ph_override=None,
                  season=120.0, n_t=241, measured_forcing=True, biomass="growth_rice"):
    """Source apportionment for one congener -- the PFAS analog of dynamiCROP's Fig. 2.

    Builds the SAME model as ``simulate()`` (identical compound / composition / driver
    resolution) and integrates the inter-compartment MASS fluxes over the season, so the
    cumulative PFAS mass delivered to each compartment is decomposed by transport pathway:
    grain by xylem vs phloem, root by soil-uptake vs phloem recirculation. The headline
    audit number is ``fraction['grain']['phloem_from_leaf']`` -- the model asserts the
    grain is phloem-fed (loading ``L_Ph``), and this quantifies it.

    Unlike dynamiCROP (linear -> split by initial source pool), this model is nonlinear
    (GHK + Michaelis-Menten) and continuously sourced, so the decomposition is flux-based
    and by pathway. Returns the ``apportionment.apportion`` dict plus congener/params.
    Scoped to the canonical 4-compartment core (`simulate`); the redistributed-shoot model
    (`simulate_nstem_leaf`) would need its own flux map.
    """
    import apportionment as ap
    if congener not in _CONG:
        raise KeyError(f"unknown congener {congener!r}; known: {CONGENERS}")
    t = np.linspace(0.0, float(season), int(n_t))
    Cwo_series, Qtp, M, leaf_loss = _default_drivers(t, season, Cwo, measured_forcing, biomass)
    inputs = PlantInputs(t=t, Cwo=Cwo_series, Qtp=Qtp, M=M, leaf_loss=leaf_loss)
    cmpd = _compound_for(congener, f_xy_source, f_xy_override, L_Ph_override,
                         kappa_d_override, lipid_loading, g_xy_override, g_ph_override,
                         K_cw_organ="root")
    comps = _compartments()
    env = Environment(E=E_m_mV / 1000.0)
    model = RiceUptakeModel(env=env, cmpd=cmpd, comps=comps, inputs=inputs)
    sol = model.solve(t)
    res = ap.apportion(model, sol)
    res.update(congener=congener, success=bool(sol.success),
               params=dict(f_xy=float(cmpd.f_xy), L_Ph=float(cmpd.L_Ph),
                           n_C=_CONG[congener]["n_C"], group=_CONG[congener]["group"]))
    return res


def simulate_from_smiles(smiles, *, name=None, f_xy=None, f_xy_override=None, **kw):
    """Run the 4-compartment model for an arbitrary PFAS given by a SMILES string.

    Uses the RDKit structure adapter (``pfas_structure``): a structure that matches
    a known congener is run through the canonical curated parameters; a novel
    structure gets binding from the QSPR and a head-group-offset f_xy (PROVISIONAL).
    Extra keyword args (Cwo, season, drivers, lipid_loading, ...) pass to ``simulate``.

    Returns the usual ``simulate`` dict plus ``descriptors`` (the parsed structure)
    and ``provisional`` (True for novel / non-calibrated structures).  Requires RDKit.
    """
    from pfas_structure import compound_from_smiles
    cmpd, d = compound_from_smiles(smiles, name=name, f_xy=f_xy)
    fxy_ov = f_xy_override if f_xy_override is not None else f_xy
    known = name or d.matched_name
    if known in _CONG:
        res = simulate(known, f_xy_override=fxy_ov, **kw)
    else:
        group = ("PFSA" if d.head_group == "sulfonate"
                 else "ether" if d.n_ether_O > 0 else "PFCA")
        record = {"name": cmpd.name, "n_C": d.n_C, "group": group,
                  "K_prot_Lkg": cmpd.K_prot, "K_PL_Lkg": cmpd.K_PL,
                  "K_cw_wholecw_Lkg": {"root": cmpd.K_cw},
                  "f_xy_recommended": cmpd.f_xy, "f_xy_W2fit": None,
                  "kappa_d_W2fit": cmpd.kappa_d, "L_Ph_W2fit": cmpd.L_Ph}
        res = simulate(cmpd.name, record=record, f_xy_override=fxy_ov, **kw)
    res["descriptors"] = d
    res["provisional"] = (d.matched_name is None) or (not d.is_linear)
    return res


def _transport_defaults(c, f_xy_source, lipid_loading):
    """Resolve the (f_xy, L_Ph, kappa_d, g_xy, g_ph) defaults for a congener record.

    f_xy_source: 'recommended' (monotone physical TSCF), 'W2fit' (per-congener fit on
    the placeholder/growth_rice driver -- reproduces Yamazaki there; reproduce_demo),
    or 'oryza' (per-congener fit on the mechanistic ORYZA2000 biomass -- the
    reproduction calibration for the new default driver; validation/refit_oryza.py).
    """
    if lipid_loading:
        f_xy, g_xy, g_ph = lipid_loading_conductances(c["n_C"], c["K_PL_Lkg"], c["group"])
        return f_xy, (c.get("L_Ph_W2fit") or 0.01), (c.get("kappa_d_W2fit") or 2.0), g_xy, g_ph
    if f_xy_source == "oryza":
        return (c.get("f_xy_oryza") or c["f_xy_recommended"],
                c.get("L_Ph_oryza") or (c.get("L_Ph_W2fit") or 0.01),
                c.get("kappa_d_oryza") or (c.get("kappa_d_W2fit") or 2.0), 0.0, 0.0)
    f_xy = c["f_xy_recommended"] if f_xy_source == "recommended" else (c.get("f_xy_W2fit") or c["f_xy_recommended"])
    return f_xy, (c.get("L_Ph_W2fit") or 0.01), (c.get("kappa_d_W2fit") or 2.0), 0.0, 0.0


def _compound_for(congener, f_xy_source="recommended", f_xy_override=None,
                  L_Ph_override=None, kappa_d_override=None, lipid_loading=False,
                  g_xy_override=None, g_ph_override=None, K_cw_organ="stem"):
    """Build a Compound for `congener` with the same parameter resolution as
    simulate() (f_xy source / lipid loading / overrides). `K_cw_organ` picks the
    cell-wall K used (the nstem stem segments use the stem cell wall)."""
    c = _CONG[congener]
    f_xy_def, L_Ph_def, kappa_d_def, g_xy_def, g_ph_def = _transport_defaults(c, f_xy_source, lipid_loading)
    f_xy = float(f_xy_override) if f_xy_override is not None else float(f_xy_def)
    g_xy = float(g_xy_override) if g_xy_override is not None else float(g_xy_def)
    g_ph = float(g_ph_override) if g_ph_override is not None else float(g_ph_def)
    kappa_d = float(kappa_d_override) if kappa_d_override is not None else float(kappa_d_def)
    L_Ph = float(L_Ph_override) if L_Ph_override is not None else float(L_Ph_def)
    return Compound(name=congener, K_prot=c["K_prot_Lkg"], K_PL=c["K_PL_Lkg"],
                    K_cw=c["K_cw_wholecw_Lkg"][K_cw_organ], kappa_d=kappa_d,
                    Vmax_in=_CARR["Vmax_in"], Km_in=_CARR["Km_in"],
                    Vmax_out=_CARR["Vmax_out"], Km_out=_CARR["Km_out"],
                    L_Ph=L_Ph, f_xy=f_xy, g_xy=g_xy, g_ph=g_ph)


def simulate_nstem_leaf(congener="PFOA", Cwo=1.0, E_m_mV=-120.0,
                        f_xy_source="recommended", f_xy_override=None,
                        L_Ph_override=None, kappa_d_override=None,
                        lipid_loading=False, g_xy_override=None, g_ph_override=None,
                        season=150.0, n_t=361, N=4, stem_transp_frac=0.45,
                        lam_grain=0.05, retention=0.6, biomass_fn=None):
    """Redistributed-shoot uptake run (N stem segments + explicit leaf), the Tang
    2026 over-translocation fix (see `pfas_rice_plant_module_nstem_leaf`).

    Mirrors `simulate()` for the congener/transport parameters but resolves the
    stem into N transpiration terminals and applies the deposition+retention
    redistribution so the leaf no longer monopolizes the shoot burden. Returns a
    dict with the SAME tissue keys as `simulate()` (root/stem/leaf/grain), where
    `stem` is the mass-weighted mean over the N segments, so the Tang validation
    and plots can treat the two models interchangeably.

    stem_transp_frac : fraction of canopy transpiration terminating on the stalk
        (the leaf takes the remainder, minus `lam_grain`). Crop-architecture lever.
    retention : fraction of each organ's transpiration deposit that is retained
        (terminal); 1-retention flows on to the grain as residual xylem.
    biomass_fn : callable (t, season) -> dict{root,stem,leaf,grain} [kg/hill] for the
        organ biomass driver. Default None -> the mechanistic ORYZA2000 (`oryza_growth`).
        Pass `gr.organ_biomass` for the lightweight growth_rice reconstruction.
    """
    from pfas_rice_plant_module_nstem_leaf import (
        NStemLeafModel, PlantInputsNL, make_stem_leaf_compartments, split_from_stem_frac)
    if congener not in _CONG:
        raise KeyError(f"unknown congener {congener!r}; known: {CONGENERS}")

    t = np.linspace(0.0, float(season), int(n_t))
    Qtp = fr.Q_TP(t, season)
    b = (biomass_fn or _biomass_fn("oryza"))(t, season)   # default: mechanistic ORYZA2000
    M = np.column_stack(
        [np.maximum(b["root"], 1e-9)]
        + [np.maximum(b["stem"] / N, 1e-9)] * N
        + [np.maximum(b["leaf"], 1e-9), np.maximum(b["grain"], 1e-9)])
    inputs = PlantInputsNL(t=t, Cwo=np.full_like(t, float(Cwo)), Qtp=Qtp, M=M,
                           leaf_loss=b.get("leaf_death_rate"))

    g = _COMP
    _kw = lambda d: dict(theta=d["theta_fw"], f_prot=d["f_prot"], f_PL=d["f_PL"], f_cw=d["f_cw"])
    comps = make_stem_leaf_compartments(
        N, _kw(g["stem"]), _kw(g["root"]), _kw(g["leaf"]), _kw(g["grain_brown"]))
    cmpd = _compound_for(congener, f_xy_source, f_xy_override, L_Ph_override,
                         kappa_d_override, lipid_loading, g_xy_override, g_ph_override,
                         K_cw_organ="stem")
    tau, lam_leaf, lam_grain = split_from_stem_frac(N, stem_transp_frac, lam_grain)
    env = Environment(E=E_m_mV / 1000.0)
    model = NStemLeafModel(env=env, cmpd=cmpd, comps=comps, inputs=inputs, tau=tau,
                           lam_leaf=lam_leaf, lam_grain=lam_grain, retention=retention)
    sol = model.solve(t)
    Y = sol.y                                            # (N+3, n_t)
    seg = slice(1, N + 1)
    Mseg = M[:, seg]
    stem_conc = np.sum(Y[seg, :] * Mseg.T, axis=0) / np.sum(Mseg.T, axis=0)
    conc = {"root": Y[0], "stem": stem_conc, "leaf": Y[model.LEAF], "grain": Y[model.GRAIN]}
    cwo_ref = float(Cwo) if Cwo else 1.0
    baf_final = {k: float(v[-1] / cwo_ref) for k, v in conc.items()}
    return dict(
        t=t, congener=congener, success=bool(sol.success), season=float(season),
        conc=conc,
        baf={k: v / cwo_ref for k, v in conc.items()},
        baf_final=baf_final,
        tf_final={k: baf_final[k] / baf_final["root"] for k in conc},
        M=M, N=N,
        params=dict(f_xy=float(cmpd.f_xy), retention=float(retention),
                    stem_transp_frac=float(stem_transp_frac), lam_leaf=float(lam_leaf),
                    lam_grain=float(lam_grain), n_C=_CONG[congener]["n_C"],
                    group=_CONG[congener]["group"]),
    )


# ---------------------------------------------------------------------------
# SEQUESTRATION two-pool root model (EXPLORATORY, opt-in) -- mobile + sequestered
#   docs/twopool_root_exploration.md ; validation/twopool_root_exploration.py
# Splits the root into a MOBILE pool (binding B_m; feeds the xylem via the
# monotone-physical f_xy + a K_PL-gated lipid term) and a SEQUESTERED pool (an
# irreversible apoplast/cell-wall/Fe-Mn-plaque terminal sink at rate k_seq, a
# non-K_PL U-shaped chain/head-group descriptor). The decoupling lets the model
# hold a HIGH long-chain root BAF *and* deliver to the shoot while KEEPING the
# monotone physical f_xy. In-sample Yamazaki fit (RMSE 0.251); canonical core and
# parameters.json are UNCHANGED -- the fitted globals + U-shape come from the
# cached fit (validation/twopool_fitted_params.json) loaded lazily below.
# NOTE: distinct from the CARRIER two-pool model (`simulate_twopool_carrier`,
# src/pfas_rice_two_pool.py) -- that one's second pool is a reversible bound store
# tuned by carrier/f_xy levers; this one's is an irreversible k_seq sink.
# ---------------------------------------------------------------------------
_TP_RM, _TP_RS, _TP_ST, _TP_LF, _TP_GR = 0, 1, 2, 3, 4   # 5-state indices
_TWOPOOL_SEQ = None                                      # ((p, q), module) cache


def _twopool_seq():
    """Lazily import the two-pool exploration module and load its cached fit.

    Returns ((p, q), TP) where `p` are the fitted GLOBAL params, `q` the U-shaped
    k_seq coefficients (validation/twopool_fitted_params.json) and `TP` the module
    (for its pure descriptor functions kseq_ushape/lipid_g). Reusing the cached fit
    keeps this wrapper byte-consistent with the validation script."""
    global _TWOPOOL_SEQ
    if _TWOPOOL_SEQ is None:
        import sys
        vdir = os.path.join(_ROOT, "validation")
        if vdir not in sys.path:
            sys.path.insert(0, vdir)
        import twopool_root_exploration as TP
        _TWOPOOL_SEQ = (TP.load_fit(), TP)
    return _TWOPOOL_SEQ


def _twopool_seq_rhs(t_grid, Cwo_s, Qtp_s, M_s, dM_s, cmpd, comps, B, env,
                     gxy, gph, kseq, phi=0.1, T_C_Ph=10.0, k_rel=0.0):
    """5-state RHS (mass-conserving; sole source M_root*j_R into the mobile pool).

    Faithful re-implementation of validation/twopool_root_exploration.make_rhs,
    generalised to arbitrary driver arrays (Cwo, Qtp, M on `t_grid`) so it returns
    a full time series and supports the standard simulate() driver options. Guarded
    against drift by tests/test_model_api.py::test_simulate_twopool_seq_matches_validation_and_rmse."""
    def rhs(t, C):
        Qtp = float(np.interp(t, t_grid, Qtp_s))
        Cwo = float(np.interp(t, t_grid, Cwo_s))
        M = np.array([np.interp(t, t_grid, M_s[:, k]) for k in range(4)])
        dM = np.array([np.interp(t, t_grid, dM_s[:, k]) for k in range(4)])
        M = np.maximum(M, 1e-12)
        mu = dM / M                                       # growth dilution [1/day]
        Mr = M[ROOT]

        Cw = np.empty(5)
        Cw[_TP_RM] = C[_TP_RM] / B[ROOT]                  # mobile-root free conc
        Cw[_TP_ST] = C[_TP_ST] / B[STEM]
        Cw[_TP_LF] = C[_TP_LF] / B[LEAF]
        Cw[_TP_GR] = C[_TP_GR] / B[FRUIT]

        A3 = comps[LEAF].S * M[LEAF]; A4 = comps[FRUIT].S * M[FRUIT]
        split = A3 / (A3 + A4) if (A3 + A4) > 0 else 0.5
        f3, f4 = split, 1.0 - split

        Q_Phl = max(dM[FRUIT] * T_C_Ph + phi * Qtp, 0.0)
        C_Phl = cmpd.L_Ph * Cw[_TP_LF] + gph * C[_TP_LF]
        Cw_xyl = cmpd.f_xy * Cw[_TP_RM] + gxy * C[_TP_RM]

        jR = root_uptake(Cwo, Cw[_TP_RM], cmpd, env)
        seq = kseq * C[_TP_RM]                            # mobile -> seq
        rel = k_rel * C[_TP_RS]                           # seq -> mobile (slow desorption)

        dC = np.zeros(5)
        dC[_TP_RM] = (jR - (Qtp / Mr) * Cw_xyl + phi * (Q_Phl / Mr) * C_Phl
                      - seq + rel - mu[ROOT] * C[_TP_RM])
        dC[_TP_RS] = seq - rel - mu[ROOT] * C[_TP_RS]     # near-terminal accumulator
        dC[_TP_ST] = (Qtp / M[STEM]) * (Cw_xyl - Cw[_TP_ST]) - mu[STEM] * C[_TP_ST]
        dC[_TP_LF] = (f3 * (Qtp / M[LEAF]) * Cw[_TP_ST]
                      - (1.0 + phi) * (Q_Phl / M[LEAF]) * C_Phl - mu[LEAF] * C[_TP_LF])
        dC[_TP_GR] = (f4 * (Qtp / M[FRUIT]) * Cw[_TP_ST]
                      + (Q_Phl / M[FRUIT]) * C_Phl - mu[FRUIT] * C[_TP_GR])
        return dC
    return rhs


def simulate_twopool_seq(congener="PFOA", Cwo=1.0, E_m_mV=-120.0, season=120.0,
                         n_t=481, measured_forcing=False, biomass="growth_rice",
                         drivers=None, k_rel=0.0, kseq_override=None):
    """Sequestration two-pool root run (mobile + sequestered root pools) -- EXPLORATORY, opt-in.

    Mirrors `simulate()` / `simulate_nstem_leaf()` for the I/O so the app and other
    validation can treat the two-pool model interchangeably: returns a dict with the
    SAME tissue keys (root/stem/leaf/grain), where the reported `root` BAF is the
    SUM of the mobile and sequestered root pools. The transport globals (kappa_d,
    L_Ph, the lipid conductances) and the non-K_PL U-shaped `k_seq(n, head_group)`
    come from the cached Yamazaki fit (validation/twopool_fitted_params.json); the
    xylem free-loading uses the monotone physical `f_xy_recommended`.

    NOTE: distinct from `simulate_twopool_carrier` (the long-chain-closure model,
    src/pfas_rice_two_pool.py) -- that one's second pool is a reversible bound store
    tuned by carrier/f_xy levers; here it is an irreversible k_seq sink.

    NOTE the cached fit is on the DEMO forcings, so the defaults
    (`measured_forcing=False`, `season=120`) reproduce the documented in-sample
    result (overall log10 RMSE 0.251; PFOS/PFUnDA root separation at identical K_PL).
    `measured_forcing=True` / `drivers=` / `biomass=` are supported for exploration
    but pair the demo fit with non-demo forcings (Result 6 has a separate measured
    re-fit, not auto-loaded). parameters.json is UNCHANGED.

    k_rel : slow seq->mobile desorption rate [1/day]; 0 (default) = irreversible seq
        sink. kseq_override : bypass the U-shaped descriptor with a fixed k_seq [1/day].
    """
    if congener not in _CONG:
        raise KeyError(f"unknown congener {congener!r}; known: {CONGENERS}")
    (p, q), TP = _twopool_seq()
    c = _CONG[congener]

    if drivers is not None:
        t = np.asarray(drivers["t"], dtype=float)
        Cwo_series = np.asarray(drivers["Cwo"], dtype=float)
        Qtp = np.asarray(drivers["Qtp"], dtype=float)
        M = np.asarray(drivers["M"], dtype=float)
        season = float(t[-1])
    else:
        t = np.linspace(0.0, float(season), int(n_t))
        Cwo_series, Qtp, M, _ = _default_drivers(t, season, Cwo, measured_forcing, biomass)

    comps = _compartments()
    env = Environment(E=E_m_mV / 1000.0)
    cmpd = Compound(name=congener, K_prot=c["K_prot_Lkg"], K_PL=c["K_PL_Lkg"],
                    K_cw=c["K_cw_wholecw_Lkg"]["root"], kappa_d=p["kappa_d"],
                    Vmax_in=_CARR["Vmax_in"], Km_in=_CARR["Km_in"],
                    Vmax_out=_CARR["Vmax_out"], Km_out=_CARR["Km_out"],
                    L_Ph=p["L_Ph"], f_xy=c["f_xy_recommended"])
    B = binding_factors(comps, cmpd)
    gxy, gph = TP.lipid_g(c["K_PL_Lkg"], c["group"], p["gxy"], p["gph"], p["K_half"], p["pfsa_ln"])
    kseq = float(kseq_override) if kseq_override is not None else \
        float(TP.kseq_ushape(c["n_C"], c["group"], q))

    dM = np.gradient(M, t, axis=0)
    rhs = _twopool_seq_rhs(t, Cwo_series, Qtp, M, dM, cmpd, comps, B, env,
                           gxy, gph, kseq, k_rel=k_rel)
    from scipy.integrate import solve_ivp
    sol = solve_ivp(rhs, (t[0], t[-1]), np.zeros(5), t_eval=t,
                    method="BDF", rtol=1e-6, atol=1e-9)
    Y = sol.y                                            # (5, n_t)
    root_total = Y[_TP_RM] + Y[_TP_RS]
    conc = {"root": root_total, "stem": Y[_TP_ST], "leaf": Y[_TP_LF], "grain": Y[_TP_GR],
            "root_mobile": Y[_TP_RM], "root_seq": Y[_TP_RS]}
    Mf = M[-1]
    straw = (Y[_TP_ST] * Mf[STEM] + Y[_TP_LF] * Mf[LEAF]) / (Mf[STEM] + Mf[LEAF])
    cwo_ref = float(Cwo_series[-1]) if Cwo_series[-1] > 0 else (
        float(np.nanmax(Cwo_series)) if np.nanmax(Cwo_series) > 0 else 1.0)
    baf_final = {k: float(v[-1] / cwo_ref) for k, v in conc.items()}
    return dict(
        t=t, congener=congener, success=bool(sol.success), season=float(season),
        conc=conc, straw=straw,
        Cwo=Cwo_series, Qtp=Qtp, M=M,
        baf={k: v / cwo_ref for k, v in conc.items()},
        baf_final=baf_final,
        straw_baf=float(straw[-1] / cwo_ref),
        tf_final={k: baf_final[k] / baf_final["root"] for k in ("root", "stem", "leaf", "grain")},
        seq_fraction=float(Y[_TP_RS, -1] / max(root_total[-1], 1e-12)),
        cwo_ref=cwo_ref,
        B_k={k: float(B[i]) for i, k in enumerate(TISSUES)},
        N=float(env.N), eN=float(np.exp(env.N)),
        params=dict(f_xy=float(cmpd.f_xy), L_Ph=float(p["L_Ph"]), kappa_d=float(p["kappa_d"]),
                    g_xy=float(gxy), g_ph=float(gph), k_seq=kseq, k_rel=float(k_rel),
                    K_PL=c["K_PL_Lkg"], n_C=c["n_C"], group=c["group"]),
    )


# ---------------------------------------------------------------------------
# Schematic / colormap helpers (UI-agnostic; used by plots.fig_plant_schematic)
# ---------------------------------------------------------------------------
def metric_series(res, metric="conc"):
    """Per-tissue series + shared colour limits for the plant-map colormap.

    metric : 'conc' -> tissue concentration [ug/kg]
             'baf'  -> bioaccumulation factor C_k/C_w^o [L/kg]
    Returns dict(data={tissue: array}, label, cmin, cmax). The limits span ALL
    tissues over the WHOLE season so the colour scale is stable while scrubbing
    the time slider (you can compare compartments and times on one colorbar).
    """
    if metric == "baf":
        data = {k: np.asarray(res["baf"][k], float) for k in TISSUES}
        label = "BAF  [L/kg]"
    else:
        data = {k: np.asarray(res["conc"][k], float) for k in TISSUES}
        label = "tissue conc  [µg/kg]"
    finite = [v[np.isfinite(v)] for v in data.values()]
    allv = np.concatenate(finite) if finite else np.array([0.0, 1.0])
    cmin = float(np.nanmin(allv)) if allv.size else 0.0
    cmax = float(np.nanmax(allv)) if allv.size else 1.0
    return dict(data=data, label=label, cmin=cmin, cmax=max(cmax, cmin + 1e-12))


def schematic_values(res, metric="conc", t_index=-1):
    """Per-compartment scalar (+ Cwo) at one time index, for the plant map."""
    ms = metric_series(res, metric)
    ti = int(t_index)
    vals = {k: float(ms["data"][k][ti]) for k in TISSUES}
    vals["straw"] = float(
        (res["conc"]["stem"][ti] * res["M"][ti, STEM]
         + res["conc"]["leaf"][ti] * res["M"][ti, LEAF])
        / (res["M"][ti, STEM] + res["M"][ti, LEAF]))
    if metric == "baf":
        cwv = res["Cwo"][ti]
        vals["straw"] = vals["straw"] / cwv if cwv > 0 else np.nan
    return dict(values=vals, label=ms["label"], cmin=ms["cmin"], cmax=ms["cmax"],
                Cwo=float(res["Cwo"][ti]), t=float(res["t"][ti]))


# ---------------------------------------------------------------------------
# Driver builders -- the three ways to feed the plant ODE
#   (1) MODEL          : built-in measured/placeholder forcings (simulate(...))
#   (2) HYDRUS / CSV   : drivers_from_arrays / load_driver_csv  (one-way coupling)
#   (3) SOIL INVENTORY : pore_water_from_inventory (Freundlich) + measured Q/M
# ---------------------------------------------------------------------------
def measured_forcing(t, season=None):
    """The measured transpiration Q_TP(t) [L/d/hill] and ORYZA organ biomass M(t)
    [kg/hill, (n,4)] on grid `t`. Reused by the HYDRUS/soil modes when the user
    supplies only a concentration series (Q_TP and M are crop physiology, not PFAS)."""
    t = np.asarray(t, float)
    season = float(t[-1]) if season is None else float(season)
    Qtp = fr.Q_TP(t, season)
    b = gr.organ_biomass(t, season)
    M = np.maximum(np.column_stack([b["root"], b["stem"], b["leaf"], b["grain"]]), 1e-4)
    return Qtp, M


def drivers_from_arrays(t, Cwo, Qtp=None, M=None, season=None):
    """Assemble a `drivers` dict for `simulate(drivers=...)`.

    Only `Cwo(t)` is PFAS-specific (from HYDRUS / soil inversion / a probe). If
    `Qtp`/`M` are omitted they default to the measured crop-physiology forcings on
    the same grid, so a bare concentration series is enough to run the plant model.
    """
    t = np.asarray(t, float)
    Cwo = np.asarray(Cwo, float)
    if Qtp is None or M is None:
        Q_def, M_def = measured_forcing(t, season)
        Qtp = Q_def if Qtp is None else np.asarray(Qtp, float)
        M = M_def if M is None else np.asarray(M, float)
    return dict(t=t, Cwo=Cwo, Qtp=np.asarray(Qtp, float), M=np.asarray(M, float))


_DRIVER_COLS = ("t", "Cwo", "Qtp", "M_root", "M_stem", "M_leaf", "M_grain")


def load_driver_csv(path_or_buffer):
    """Load a HYDRUS/external driver table into a `drivers` dict.

    Required columns: t, Cwo, Qtp, M_root, M_stem, M_leaf, M_grain (see
    `docs/visualization_tool.md` for the HYDRUS-1D output mapping). Accepts a path
    or any file-like object. Missing Qtp/M columns fall back to measured forcings.
    """
    data = np.genfromtxt(path_or_buffer, delimiter=",", names=True)
    cols = data.dtype.names or ()
    if "t" not in cols or "Cwo" not in cols:
        raise ValueError(f"driver CSV needs at least 't' and 'Cwo' columns; found {cols}")
    t = np.atleast_1d(data["t"]).astype(float)
    Cwo = np.atleast_1d(data["Cwo"]).astype(float)
    have_M = all(c in cols for c in ("M_root", "M_stem", "M_leaf", "M_grain"))
    M = (np.column_stack([data["M_root"], data["M_stem"], data["M_leaf"], data["M_grain"]]).astype(float)
         if have_M else None)
    Qtp = np.atleast_1d(data["Qtp"]).astype(float) if "Qtp" in cols else None
    return drivers_from_arrays(t, Cwo, Qtp=Qtp, M=M)


def pore_water_from_inventory(t, C_total, K_F=2.0, n=0.85, theta_g=0.35,
                              theta_g_flooded=0.60, flooded=None, k_leach=0.0):
    """Invert a soil PFAS inventory to pore-water C_w^o(t) via a Freundlich paddy soil.

    C_total : total soil inventory [µg/kg dry] (scalar or series on `t`).
    K_F, n  : Freundlich capacity / exponent (S = K_F*C_w^n); long chains sorb harder.
    theta_g : drained gravimetric water content [L/kg dry].
    flooded : optional boolean schedule on `t`; when given, flooded steps use the
        redox-corrected soil (higher water content -> dilution) plus first-order
        `k_leach` on the dissolved pool. None -> single drained isotherm.
    Returns (Cwo_series, soil_object) -- soil_object exposes the isotherm for plots.
    """
    t = np.asarray(t, float)
    C_total = np.full_like(t, float(C_total)) if np.ndim(C_total) == 0 else np.asarray(C_total, float)
    if flooded is None:
        soil = FreundlichSoil(K_F=K_F, n=n, theta_g=theta_g, name="paddy soil")
        Cwo = soil.pore_water_series(C_total)
        return Cwo, soil
    flooded = np.asarray(flooded, bool)
    redox = PaddyRedoxCorrected(
        drained=FreundlichSoil(K_F=K_F, n=n, theta_g=theta_g, name="drained/aerobic"),
        flooded=FreundlichSoil(K_F=K_F, n=n, theta_g=theta_g_flooded, name="flooded (diluted)"),
        k_leach=k_leach)
    C_T = redox.apply_leaching(t, C_total, flooded) if k_leach > 0 else C_total
    Cwo = redox.pore_water_series(C_T, flooded)
    return Cwo, redox


# ---------------------------------------------------------------------------
# Real HYDRUS-1D soil run (live coupling) -- optional, needs the built engine
# ---------------------------------------------------------------------------
def rdkit_available():
    """True if RDKit is importable (needed for SMILES/structure input).

    The app's 'SMILES (structure)' compound input parameterises any PFAS from its
    structure via `src/pfas_structure.py` (RDKit). Gate the UI on this so the rest
    of the tool stays usable when RDKit is absent."""
    try:
        import rdkit  # noqa: F401
    except Exception:
        return False
    return True


def hydrus_available():
    """True if the compiled HYDRUS-1D engine + phydrus are usable here.

    The live 'Run HYDRUS' mode runs a genuine HYDRUS-1D paddy soil model
    (`src/soil_hydrus.py`); it needs the engine built from `external/hydrus_source`
    (gfortran) and `pip install phydrus`. Gate the UI on this so the rest of the
    tool stays usable when the engine is absent (fresh clone / Streamlit Cloud)."""
    try:
        import soil_hydrus as sh
    except Exception:
        return False
    return bool(sh.hydrus_available())


def build_hydrus_engine():
    """Best-effort in-place build of the HYDRUS-1D engine (gfortran). Returns
    (ok, log_lines). Used by the app's 'Run HYDRUS-1D (live)' mode so the engine
    can be compiled on Streamlit Cloud (needs gfortran/make via packages.txt)."""
    try:
        import soil_hydrus as sh
    except Exception as e:                                   # noqa: BLE001
        return False, [f"soil_hydrus import failed: {e}"]
    return sh.build_hydrus_engine()


def hydrus_drivers(congener, season=120.0, Cwo_ref=1.0, f_oc=0.02, n_t=241,
                   qtp_from_hydrus=True, **run_kw):
    """Run a real HYDRUS-1D paddy soil model for `congener` and return a `drivers`
    dict (+ the raw PaddyResult) for `simulate(drivers=…)`.

    HYDRUS supplies the congener-dependent pore-water `Cwᵒ(t)` (short chains leach
    under flooding, long chains stay buffered) and the actual root water uptake
    `Q_TP(t)`; `M(t)` is ORYZA biomass. `Cwᵒ(t)` is normalised to season-mean
    `Cwo_ref` so the average exposure matches a constant-Cwo run. Extra `run_kw`
    (flood_until, percolation, …) pass through to `soil_hydrus.run_paddy_hydrus`."""
    import soil_hydrus as sh
    if congener not in _CONG:
        raise KeyError(f"unknown congener {congener!r}; known: {CONGENERS}")
    c = _CONG[congener]
    pin, res = sh.inputs_from_hydrus(c["n_C"], c["group"], season=season,
                                     Cwo_ref=Cwo_ref, f_oc=f_oc, n_t=n_t,
                                     qtp_from_hydrus=qtp_from_hydrus, **run_kw)
    drivers = dict(t=np.asarray(pin.t, float), Cwo=np.asarray(pin.Cwo, float),
                   Qtp=np.asarray(pin.Qtp, float), M=np.asarray(pin.M, float))
    return drivers, res


# ---------------------------------------------------------------------------
# Biomonitoring -- measured tissue concentrations, no soil model needed
# ---------------------------------------------------------------------------
def baf_from_measurement(conc_by_tissue, Cwo):
    """BAF_k = C_k(measured) / C_w^o for each supplied tissue [L/kg].

    The biomonitoring path: when you have measured tissue concentrations (and a
    measured pore-water / soil-solution concentration) HYDRUS is not needed -- the
    BAF is read straight off the data. Returns {} if Cwo<=0."""
    if not Cwo or Cwo <= 0:
        return {}
    return {k: float(v) / float(Cwo) for k, v in conc_by_tissue.items()
            if v is not None and np.isfinite(v)}


def load_biomonitoring_csv(path_or_buffer):
    """Load measured tissue concentrations: columns `tissue,conc` (+ optional `Cwo`).

    Returns dict(conc={tissue: value}, Cwo=float|None). `tissue` is free text but
    root/straw/stem/leaf/grain map onto the plant schematic."""
    import csv as _csv
    if hasattr(path_or_buffer, "read"):
        text = path_or_buffer.read()
        if isinstance(text, bytes):
            text = text.decode("utf-8")
        lines = text.splitlines()
    else:
        with open(path_or_buffer) as f:
            lines = f.read().splitlines()
    conc, Cwo = {}, None
    for row in _csv.DictReader(lines):
        tis = (row.get("tissue") or "").strip().lower()
        if not tis:
            continue
        try:
            conc[tis] = float(row["conc"])
        except (KeyError, ValueError, TypeError):
            continue
        if Cwo is None and row.get("Cwo"):
            try:
                Cwo = float(row["Cwo"])
            except (ValueError, TypeError):
                pass
    return dict(conc=conc, Cwo=Cwo)


if __name__ == "__main__":
    r = simulate("PFOA")
    print("PFOA:", {k: round(v, 3) for k, v in r["baf_final"].items()},
          "straw_baf", round(r["straw_baf"], 3), "B_k", {k: round(v, 2) for k, v in r["B_k"].items()})
