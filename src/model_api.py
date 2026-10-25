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
    binding_factors, root_uptake, xylem_leaf_share, _logistic, ROOT, STEM, LEAF, FRUIT)
import forcing_rice as fr
import growth_rice as gr
import organ_mass as om
import harvest as hv
import site_scenario as ss
from organ_mass import MASS_BASIS          # crop-model DRY matter -> ODE FRESH mass (tag)
from soil_paddy import FreundlichSoil
from soil_paddy_redox_corrected import PaddyRedoxCorrected

with open(os.path.join(_ROOT, "params", "parameters.json")) as _f:
    PARAMS = json.load(_f)
_CARR = PARAMS["carrier_MichaelisMenten"]
_COMP = PARAMS["tissue_composition_recommended"]
_CONG = {c["name"]: c for c in PARAMS["congeners"]}

CONGENERS = [c["name"] for c in PARAMS["congeners"]]            # 12, ordered
TISSUES = ("root", "stem", "leaf", "grain")
# built-in drivers: the calibration site-year (default) or the generic climatology drivers
SCENARIOS = ("yamazaki2020", "generic")
DEFAULT_SCENARIO = "yamazaki2020"

# ---------------------------------------------------------------------------
# ROOT-ENTRY MECHANISM -- a named mode, because the choice is NOT settled.
#
# docs/theory_anchor.tex marks the Michaelis-Menten CARRIER as one of only three
# things this repo adds to the Trapp (2000, 2004) ionizable-compound cell model,
# and its constants were never measured (parameters.json carrier block). validation/
# carrier_vs_bypass.py compared it, on identical everything else, with a passive
# APOPLASTIC BYPASS of the DPU base's own Fickian form -- one global parameter
# each -- and with membrane depolarisation, which adds no term at all:
#
#     nothing (Trapp's PFAS limit)      log10 RMSE 2.301
#     depolarisation (E_m = -90 mV)               1.890
#     carrier  (fixed Vmax)                       1.227
#     carrier  (global Vmax re-fitted, x8)        1.078
#     bypass   (one global g_apo = 50)            1.093
#
# So an addition IS necessary (depolarisation is REFUTED), but WHICH addition is
# not decided by that data: bootstrap P(bypass beats the re-fitted carrier) = 0.368,
# and the bypass cannot claim parsimony either (Vmax_in is a single global, not
# per-congener -- checked, not assumed).
#
# THE CARRIER IS THE DEFAULT BY DEFAULT, NOT BY EVIDENCE. It is the incumbent and
# every published number runs on it; the bypass is kept RUNNABLE rather than
# buried in a validation script, in the repo's idiom for an open question
# (`lipid_source`, `f_xy_source`, `cwo_profile`, `biomass`, `tscf_model`).
#
# Against the bypass specifically: fitted per congener, g_apo trends with chain
# length (corr +0.898, spread 100x) where theory_anchor.tex says eta -- which
# contains the apoplastic bypass -- is "essentially independent of tail length".
# On its own criterion that is a relabelled carrier. The bypass value below is
# the single global optimum from that scan, NOT a calibrated parameter.
# ---------------------------------------------------------------------------
# SHOOT MODEL -- which shoot structure carries the xylem-delivered solute.
#   "4pool"         the single-straw core (`RiceUptakeModel`): one stem pool, the leaf/grain
#                   xylem split by dry-matter surface (`xylem_leaf_share`).
#   "redistributed" N stem segments + leaf, each organ retaining the solute its own
#                   transpiration deposits (`simulate_nstem_leaf`, the Tang over-translocation
#                   fix; architecture levers stem_transp_frac / retention / lam_grain).
# DEFAULT = the redistributed shoot with K_PL-gated lipid loading, adopted 2026-10-25 by the
# pre-registered rule in validation/default_switch.py: best transfer on all three independent
# datasets (Tang 0.640, Kim 0.469, Li 0.358 / 0.510 vs the 4-pool free anion's 1.714 / 2.33 /
# 1.26 / 1.68), Yamazaki 0.380 (in-sample for the lipid constants). The old default stays
# runnable as FOUR_POOL_FREE.
SHOOTS = ("4pool", "redistributed")
DEFAULT_SHOOT = "redistributed"
DEFAULT_LIPID_LOADING = True
# the configuration every calibration record and 4-pool mechanism analysis was run on
FOUR_POOL_FREE = dict(shoot="4pool", lipid_loading=False)
# How the redistributed shoot hands the residual xylem stream to the grain: "gate" (all of
# it once the panicle forms) or "size" (a share growing with the panicle's mass). Decided by
# the pre-registered rule in validation/panicle_share_amendment.py.
GRAIN_SHARES = ("gate", "size")
DEFAULT_GRAIN_SHARE = "gate"

UPTAKE_MODES = {
    "carrier": dict(vmax_scale=1.0, g_apo=0.0),                 # incumbent = shipped
    "bypass":  dict(vmax_scale=0.0, g_apo=50.0),                # carrier off, one global g_apo
}
DEFAULT_UPTAKE = "carrier"

# Per-congener flooded-profile leaching rate k_leach, calibrated to a real HYDRUS-1D
# run (validation/cwo_kleach_calibration.py -> params/cwo_kleach.csv). Short chains
# leach faster (larger k_leach); long chains stay buffered (k_leach -> 0). Loaded once.
_KLEACH_MAX = 0.15


def _load_cwo_kleach():
    path = os.path.join(_ROOT, "params", "cwo_kleach.csv")
    table, xs, ys = {}, [], []
    try:
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("congener"):
                    continue
                p = line.split(",")
                table[p[0]] = float(p[4])
                xs.append(np.log10(float(p[3])))
                ys.append(float(p[4]))
        b, a = np.polyfit(xs, ys, 1) if len(xs) >= 2 else (0.0, 0.02)
        return table, (float(a), float(b))
    except (OSError, ValueError, IndexError):
        return {}, (0.02, 0.0)                                 # flat fallback if missing


_CWO_KLEACH, _CWO_KLEACH_FALLBACK = _load_cwo_kleach()


def default_k_leach(congener=None, n_C=None, group=None, Koc=None):
    """Per-congener default leaching rate for `cwo_profile='flooded'`, calibrated to
    HYDRUS-1D. Known congeners use the calibrated table; a novel/SMILES compound falls
    back to k_leach(log10 Koc) (short chains leach faster, long chains -> ~0).

    `Koc` [L/kg] takes that same fallback without a congener at all -- the route for a
    NEUTRAL organic, whose Koc comes from log Kow (`literature_params.koc_neutral`) or
    a measurement, not from the PFAS chain-length QSPR. The fit itself is a property of
    the paddy leaching model rather than of PFAS, so it transfers; what does not
    transfer is the per-congener table, which is why a neutral never reads it."""
    if congener and congener in _CWO_KLEACH:
        return _CWO_KLEACH[congener]
    if Koc is None and n_C is not None:
        import literature_params as lp
        Koc = lp.koc(n_C, {"PFCA": "carboxylate", "PFSA": "sulfonate",
                           "ether": "ether"}.get(group, "carboxylate"))
    if Koc is not None:
        a, b = _CWO_KLEACH_FALLBACK
        return float(np.clip(a + b * np.log10(float(Koc)), 0.0, _KLEACH_MAX))
    return 0.02
# tissue keys that compose "straw" (the bulk shoot reported in agronomy)
STRAW_PARTS = ("stem", "leaf")

# --- EXPLORATORY lipid-facilitated loading model (opt-in; default OFF) ----------
# Two superposed transport mechanisms (docs/fxy_longchain_lipid_exploration.md):
#   free anion TSCF (declining) loads the FREE conc; a K_PL-GATED, B-independent
#   "bound" term loads the membrane/lipid-associated pool into xylem (g_xy) and
#   phloem (g_ph) so high-binding long chains reach the shoot. Global params fit to
#   Yamazaki (excl. near-MQL PFDoDA); IN-SAMPLE shape/mechanism, NOT validated.
#   Six globals fitted by validation/refit_site_year.py on the Yamazaki 2020 site-year (dry-basis data,
#   harvest operator), PFDoDA excluded, riding on the per-congener (L_Ph, kappa_d) of the
#   "recommended" set (L_Ph_rec, kappa_d_rec) -- so it is not a 6-parameter model on its own.
#   In-sample log10 RMSE 0.207 over 29 cells. With the grain drawing only its growth sap
#   (separate phloem streams, refit amendment 2) the lipid PHLOEM term is now used
#   (g_ph 0.005): the free anion alone cannot carry the mid/long-chain grain.
LIPID_LOADING = dict(a=0.1177, b=0.09053, g_xy=0.01577, g_ph=0.005005, K_half=800.1, pfsa_ln=0.6213)


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


def yamazaki_obs(congener=None):
    """The Yamazaki 2023 observation table (data_obs/obs_baf_Yamazaki.csv), DRY basis.

    Returns {congener: {"baf_dw": {tissue: value}, "use": {tissue: bool}, "straw_lo",
    "straw_hi", "water_ng_L", "sigma"}} (or one congener's record). BAF is tissue pg/g dw over
    the irrigation-water pg/mL [L/kg dw]; straw is the unweighted band mean, bounded by
    straw_lo/straw_hi; `use` is False for a cell that is wholly imputed (PFBS grain)."""
    import csv
    path = os.path.join(_ROOT, "data_obs", "obs_baf_Yamazaki.csv")
    out = {}
    with open(path) as f:
        for r in csv.DictReader(line for line in f if not line.startswith("#")):
            d = out.setdefault(r["compound"], {"baf_dw": {}, "use": {}})
            d["baf_dw"][r["tissue"]] = float(r["baf_dw"])
            d["use"][r["tissue"]] = r["use"] == "1"
            d["water_ng_L"] = float(r["water_ng_L"])
            d["sigma"] = float(r["sigma"])
            if r["tissue"] == "straw":
                d["straw_lo"], d["straw_hi"] = float(r["straw_lo"]), float(r["straw_hi"])
    return out if congener is None else out.get(congener, {})


def observed_baf(congener):
    """Yamazaki 2023 root/straw/grain BAF [L/kg DRY weight] for a congener -- only the cells
    that are measurements (PFBS grain, wholly imputed, is left out) -- or {} if not measured.
    Compare with the model's dry-basis harvest values (`simulate(...)["baf_final_dw"]`)."""
    d = yamazaki_obs(congener)
    if not d:
        return {}
    return {k: v for k, v in d["baf_dw"].items() if d["use"][k]}


def yamazaki_harvest_baf(congener, runs="both", dead_leaf_attached=1.0, scenario_kw=None,
                         **sim_kw):
    """The model's version of the Yamazaki 2023 observation: dry-basis harvest BAF
    {root, straw, grain} [L/kg dw] on the 2020 site-year, exposure = the congener's
    irrigation-water level (constant), the two harvest runs (2020-09-01, 2020-09-10)
    combined by geometric mean as the data combine the two cultivars (`runs="mean"`: the
    single mean-date run). `sim_kw` pass to `simulate` (f_xy_source, overrides, ...)."""
    import site_scenario as ss
    d = yamazaki_obs(congener)
    cwo = d["water_ng_L"] / 1000.0 if d else 1.0
    hs = ["mean"] if runs == "mean" else list(ss.YAMAZAKI_HARVESTS)
    logs = {k: [] for k in ("root", "straw", "grain")}
    for h in hs:
        drv = scenario_drivers(Cwo=cwo, harvest=h, **(scenario_kw or {}))
        r = simulate(congener, drivers=drv, dead_leaf_attached=dead_leaf_attached, **sim_kw)
        for k in logs:
            logs[k].append(np.log10(max(r["baf_final_dw"][k], 1e-30)))
    return {k: float(10 ** np.mean(v)) for k, v in logs.items()}


# --- Tang 2026 per-organ TF validation (dry-weight basis) -----------------------
TANG_CONGENERS = ("PFOA", "PFOS", "GenX")
# f_xy re-fit to Tang TF at the 0.1 ug/g dose (OVERRIDE-only; validation/tang2026_fxy_refit.py),
# on the same drivers `tang_tf_validation` runs (the calibration site-year).
TANG_REFIT_FXY = {"PFOA": 0.930, "PFOS": 0.999, "GenX": 0.0585}
_TANG_TF_ENDPOINT = {"TF_stalk": "stalk", "TF_leaf": "leaf", "TF_endosperm": "endosperm"}
# (model grain pool, Tang organ name, theta key) for the fresh->dry conversion
_TANG_ORGANS = (("stem", "stalk", "stem"), ("leaf", "leaf", "leaf"), ("grain", "endosperm", "grain_brown"))

# What Tang's "endosperm" is compared against. The model's grain pool is the whole
# caryopsis (BROWN rice). Tang does not report whether its endosperm was dehulled only or
# milled (the SI's white-rice moisture standard, Text S3, enters its intake calculation and does
# not describe the analysed sample; docs/rice_sample_processing_KR.md section 4), so "brown"
# stays the default and "white" is the named alternative.
#   "white" re-partitions the grain's free concentration into an endosperm of white-rice
#   composition (params/rice_tissue_params.csv grain_white: theta 0.11, protein 0.075, cell
#   wall 0.005 + lignin 0.001), assuming bran and endosperm share one free concentration.
#   Membrane phospholipid is counted the way the grain_brown row counts it: that row puts
#   ~80% of brown rice's membrane PL in the bran/embryo, so the endosperm (~88% of brown dw)
#   carries 0.2*f_PL_brown/0.88. The grain_white row's own 0.003 is mostly starch-complexed
#   lyso-PL, which the brown row does NOT count -- using it would compare two different PL pools.
TANG_ENDOSPERM = ("brown", "white")
_ENDO_PL_BRAN_SHARE, _ENDO_MASS_SHARE = 0.80, 0.88


def _grain_white_composition():
    b = _COMP["grain_brown"]
    return dict(theta_fw=0.11, f_prot=0.075, f_cw=0.006,
                f_PL=(1.0 - _ENDO_PL_BRAN_SHARE) * b["f_PL"] / _ENDO_MASS_SHARE)


def _dry_binding(comp, c, K_cw):
    """C_dw / Cw for a tissue composition: theta/(1-theta) + sum f*K (L/kg dw)."""
    th = comp["theta_fw"]
    return (th / (1.0 - th) + comp["f_prot"] * c["K_prot_Lkg"]
            + comp["f_PL"] * c["K_PL_Lkg"] + comp["f_cw"] * K_cw)


def endosperm_white_factor(congener, K_cw_organ="stem"):
    """White-rice/brown-rice DRY concentration ratio at a shared free concentration."""
    c = _CONG[congener]
    K_cw = c["K_cw_wholecw_Lkg"][K_cw_organ]
    return (_dry_binding(_grain_white_composition(), c, K_cw)
            / _dry_binding(_COMP["grain_brown"], c, K_cw))


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


TANG_SHOOTS = ("redistributed", "4pool")


def tang_tf_validation(congener, f_xy_source="recommended", use_refit=False,
                       dose="mean", lipid_loading=None, endosperm="brown", scenario=DEFAULT_SCENARIO,
                       shoot=None, grain_share=None):
    """Model per-organ transfer factor (DRY weight) vs Tang 2026, for a Tang congener.

    Tang's TF is tissue/root on dry weight; the model's is formed the same way from the
    harvest operator (`conc_dw`: stalk = stem, leaf = green + attached dead blades,
    endosperm = the grain pool). Runs the redistributed-shoot model (`simulate_nstem_leaf`)
    on `scenario` -- by default the calibration site-year, since Tang's site and dates are
    not in the repo (a transfer of drivers as well as parameters). use_refit applies the
    Tang-calibrated f_xy override. ``lipid_loading`` turns on the K_PL-gated lipid loading
    (constants fit on Yamazaki, NOT Tang -> out-of-sample; validation/oos_tang_lipid.py).
    ``endosperm`` -- 'brown' (default; the whole-caryopsis grain) or 'white' (milled rice;
    see TANG_ENDOSPERM). NOTE Tang's roots were rinsed only, so its root -- the TF
    denominator -- may carry surface-sorbed PFAS (no Fe/Mn plaque removal; same confound as
    Li 2025), which would bias every Tang TF low. ``shoot`` -- 'redistributed'
    (`simulate_nstem_leaf`) or '4pool' (the single-straw core), None -> DEFAULT_SHOOT, so
    each shoot model can be scored on the same data; ``lipid_loading`` None ->
    DEFAULT_LIPID_LOADING. ``use_refit`` is a FREE-ANION f_xy fitted to Tang, so it runs
    without lipid loading. Returns None for non-Tang congeners.
    """
    if congener not in TANG_CONGENERS:
        return None
    if endosperm not in TANG_ENDOSPERM:
        raise ValueError(f"unknown endosperm {endosperm!r}; expected one of {TANG_ENDOSPERM}")
    shoot = DEFAULT_SHOOT if shoot is None else shoot
    if shoot not in TANG_SHOOTS:
        raise ValueError(f"unknown shoot {shoot!r}; expected one of {TANG_SHOOTS}")
    lipid_loading = (False if use_refit else
                     DEFAULT_LIPID_LOADING if lipid_loading is None else bool(lipid_loading))
    fxy = TANG_REFIT_FXY[congener] if use_refit else None
    if shoot == "4pool":
        r = simulate(congener, Cwo=1.0, f_xy_source=f_xy_source, f_xy_override=fxy,
                     lipid_loading=lipid_loading, scenario=scenario, shoot="4pool")
        kcw_organ = "root"              # simulate() binds every compartment with the root K_cw
    else:
        r = simulate_nstem_leaf(congener, Cwo=1.0, f_xy_source=f_xy_source, f_xy_override=fxy,
                                lipid_loading=lipid_loading, scenario=scenario,
                                grain_share=grain_share)
        kcw_organ = "stem"              # simulate_nstem_leaf's default K_cw_organ
    dw = r["baf_final_dw"]
    model = {org: dw[mk] / dw["root"] for mk, org, _ in _TANG_ORGANS}
    if endosperm == "white":
        model["endosperm"] *= endosperm_white_factor(congener, K_cw_organ=kcw_organ)
    return dict(congener=congener, organs=[o for _, o, _ in _TANG_ORGANS],
                model_tf=model, tang_tf=tang_observed_tf(congener, dose),
                f_xy=float(r["params"]["f_xy"]), use_refit=bool(use_refit), dose=dose,
                scenario=scenario, refit_fxy=TANG_REFIT_FXY[congener], endosperm=endosperm,
                shoot=shoot, lipid_loading=lipid_loading)


def chain_table():
    """Per-congener summary parameters (for the chain-length overview)."""
    rows = []
    for c in PARAMS["congeners"]:
        rows.append(dict(name=c["name"], n_C=c["n_C"], group=c["group"],
                         K_PL=c["K_PL_Lkg"], K_prot=c["K_prot_Lkg"],
                         K_cw_root=c["K_cw_wholecw_Lkg"]["root"],
                         f_xy_recommended=c["f_xy_recommended"],
                         f_xy_fit=c.get("f_xy_fit"),
                         B_root=c["B_k_basisA_Lkg_fw"]["root"],
                         B_grain=c["B_k_basisA_Lkg_fw"]["grain"]))
    return rows


BIOMASS_DRIVERS = ("oryza", "oryza2000", "growth_rice")


def _biomass_fn(biomass="oryza2000"):
    """Organ-biomass callable (t, season) -> dict{root,stem,leaf,grain} [kg DRY/hill], for
    the "generic" scenario only (the default "yamazaki2020" scenario runs the ORYZA2000 port
    on the site's weather and dates -- `scenario_drivers`). Every place that turns the dry
    matter into the ODE's mass goes through `organ_mass` (fresh organ mass).

    'oryza2000'   -> ORYZA2000 potential production ported line by line from the manual's
                     Fortran with the IR72.D92 crop file verbatim (`oryza2000`) on the built-in
                     climatology. It reproduces the IR72DS92 experiment (phenology within 2 d,
                     aboveground biomass within ~6%; validation/oryza2000_ir72ds92.py).
    'oryza'       -> a REDUCED re-coding of the ORYZA2000 Level-1 carbon balance (`oryza_growth`,
                     anchored to a 1740 g/m2 shoot). On IR72DS92 it matures 12 d late and keeps
                     growing its stem after flowering; kept as an alternative for sensitivity.
    'growth_rice' -> ORYZA IR72 DVS-partitioning on a logistic total-biomass curve
                     (`growth_rice`; the lightweight reconstruction; no leaf death).
    """
    if biomass == "oryza":
        import oryza_growth as og
        return lambda t, s: og.organ_biomass_oryza(t, p=og.OryzaParams(season=s))
    if biomass == "oryza2000":
        import oryza2000 as o2k
        return lambda t, s: o2k.organ_biomass(t, season=s)
    if biomass == "growth_rice":
        return gr.organ_biomass
    raise ValueError(f"biomass must be one of {BIOMASS_DRIVERS}, got {biomass!r}")


def _site_run_cached():
    import functools
    import site_scenario as ss

    @functools.lru_cache(maxsize=16)
    def run(harvest, items):
        kw = dict(items)
        if harvest == "mean":
            return ss.yamazaki2020_mean_run(**kw)
        return ss.plant_run(harvest, **kw)
    return run


_SITE_RUN = None


def site_run(harvest="mean", **unknowns):
    """One Yamazaki 2020 plant run (src/site_scenario.py), cached. `harvest` = "mean" (the
    midpoint of the two harvest dates; the default single run) or a date / day of year;
    `unknowns` override the inputs the paper does not report (sbdur, tmpsb, nh, nplh, co2,
    days_heading_to_harvest, harvest_is_maturity)."""
    global _SITE_RUN
    if _SITE_RUN is None:
        _SITE_RUN = _site_run_cached()
    return _SITE_RUN(harvest, tuple(sorted(unknowns.items())))


def scenario_drivers(scenario=DEFAULT_SCENARIO, Cwo=1.0, harvest="mean", **unknowns):
    """A `simulate(drivers=...)` dict for a site-year scenario.

    "yamazaki2020" -- the PFAS calibration site-year (NARO Tsukuba, 2020): JMA weather, the
    ORYZA2000 port (IR72.D92) transplanted on the paper's date with phenology calibrated to its
    heading and harvest dates, ORYZA2000 potential transpiration from FAO-56 terms. Days after
    transplanting to harvest; organ mass FRESH; the ORYZA leaf death rate as `leaf_loss`."""
    if scenario != "yamazaki2020":
        raise ValueError(f"no site-year drivers for scenario {scenario!r}")
    import site_scenario as ss
    run = site_run(harvest, **unknowns)
    d = ss.drivers_from_run(run, Cwo=Cwo)
    return {k: (np.array(v, copy=True) if isinstance(v, np.ndarray) else v) for k, v in d.items()}


def _default_drivers(season, n_t, Cwo, measured_forcing, biomass="oryza2000",
                     scenario=DEFAULT_SCENARIO, scenario_kw=None):
    """The built-in drivers: (t, Cwo, Qtp, M, leaf_loss, M_basis).

    measured_forcing and scenario "yamazaki2020" (DEFAULT): the calibration site-year
        (`scenario_drivers`); the scenario sets the season (transplanting to harvest) and the
        crop model, so `n_t` and `biomass` do not apply and an explicit `season` other than the
        site's is an error (`site_scenario.check_season`).
    measured_forcing and scenario "generic": `forcing_rice` transpiration and the `biomass`
        driver over `season` days (default 120; alternative drivers, for sensitivity).
    not measured_forcing: the illustrative placeholder logistics (no basis of their own)."""
    if scenario not in SCENARIOS:
        raise ValueError(f"scenario must be one of {SCENARIOS}, got {scenario!r}")
    if measured_forcing and scenario == "yamazaki2020":
        d = scenario_drivers(scenario, Cwo=Cwo, **(scenario_kw or {}))
        ss.check_season(season, d["t"])
        return d["t"], d["Cwo"], d["Qtp"], d["M"], d["leaf_loss"], d["mass_basis"]
    season = 120.0 if season is None else float(season)
    t = np.linspace(0.0, season, int(n_t))
    if measured_forcing:
        Qtp = fr.Q_TP(t, season)
        b = _biomass_fn(biomass)(t, season)
        M = om.organ_matrix(b)
        leaf_loss = b.get("leaf_death_rate")            # only the senescing ORYZA driver supplies it
        basis = MASS_BASIS
    else:
        Qtp = 0.05 + 0.35 * np.exp(-((t - 0.62 * season) ** 2) / (2 * (season / 4.8) ** 2))
        M = np.column_stack([_logistic(t, 1e-3, 0.030, 0.10, season * 0.17),
                             _logistic(t, 1e-3, 0.040, 0.10, season * 0.21),
                             _logistic(t, 1e-3, 0.050, 0.12, season * 0.25),
                             _logistic(t, 1e-5, 0.025, 0.18, season * 0.67)])
        leaf_loss = None
        basis = None
    Cwo_series = np.full_like(t, float(Cwo))
    return t, Cwo_series, Qtp, M, leaf_loss, basis


def _harvest(t, conc, M, comps, M_basis, leaf_loss, attached, sol, leaf_index):
    """The harvest observation operator (src/harvest.py) on a 4-tissue solve: dead blades
    kept in the leaf and straw, dry-basis concentrations for dw-reported data."""
    dry = hv.dry_fractions({o: comps[i].theta for i, o in enumerate(TISSUES)}, M_basis)
    return hv.harvest_observables(
        t, conc, M, dry, leaf_loss=leaf_loss, attached=attached,
        C_leaf_fn=(lambda tt: sol.sol(tt)[leaf_index]) if getattr(sol, "sol", None) else None)


def simulate(congener="PFOA", Cwo=1.0, E_m_mV=-120.0, f_xy_source="recommended",
             f_xy_override=None, L_Ph_override=None, kappa_d_override=None,
             lipid_loading=None, g_xy_override=None, g_ph_override=None,
             season=None, n_t=241, measured_forcing=True, biomass="oryza2000",
             drivers=None, K_surf=0.0, record=None,
             cwo_profile="constant", cwo_kw=None,
             uptake=DEFAULT_UPTAKE, vmax_scale=None, g_apo=None, km_scale=1.0,
             dead_leaf_attached=1.0, scenario=DEFAULT_SCENARIO, scenario_kw=None,
             phi=0.1, T_C_Ph=10.0, s_grain_dry=None, shoot=None,
             N=4, stem_transp_frac=0.45, lam_grain=0.05, retention=0.6, grain_share=None):
    """Run the 4-compartment ODE for one congener and scenario.

    Parameters
    ----------
    congener : one of CONGENERS.
    Cwo : constant pore-water free concentration C_w^o [ug/L] (1.0 -> conc == BAF).
    E_m_mV : root membrane potential [mV] (GHK anion-exclusion lever).
    f_xy_source : "recommended" (DEFAULT; the monotone a-priori f_xy, with L_Ph/kappa_d
        fitted to root and grain given it, so straw is a prediction) or "fit" (all three
        fitted per congener -- reproduces Yamazaki by construction). See F_XY_SOURCES.
    f_xy_override : if given, use this f_xy instead.
    L_Ph_override, kappa_d_override : if given, use these phloem-loading / root
        membrane-conductance values instead of the fitted ones.
    scenario : "yamazaki2020" (DEFAULT; the calibration site-year -- JMA Tsukuba 2020, the
        ORYZA2000 port on the paper's dates, ORYZA2000 transpiration; `scenario_drivers`;
        `scenario_kw` overrides its unreported inputs) or "generic" (forcing_rice
        transpiration + the `biomass` driver over `season` days).
    measured_forcing : False -> the illustrative logistic placeholders (any scenario).
    biomass, season, n_t : the "generic" scenario's biomass driver and grid
        ('oryza2000' / 'oryza' / 'growth_rice'; see `_biomass_fn`).
    dead_leaf_attached : fraction of dead blades still on the plant at harvest, for the
        harvest operator (`harvest`): 1 (default) keeps them in the leaf and straw, 0 drops
        them. Not measured anywhere here -- run both ends as a sensitivity.
    phi, T_C_Ph, s_grain_dry : the phloem recirculation fraction, the phloem flow per kg of
        grain dry-matter growth [L/kg] and the grain's transpiring surface per kg dry matter
        (the xylem split; default `S_XYLEM_DRY`). None of the three is measured for rice;
        validation/refit_site_year.py runs them as sensitivities.
    drivers : optional dict with arrays {'t','Cwo','Qtp','M'} (M shape (n,4)) that
        OVERRIDE the built-in forcings -- the entry point for a HYDRUS-1D / Phydrus
        run or a soil-inventory inversion (see `drivers_from_arrays`,
        `pore_water_from_inventory`, `load_driver_csv`). When given, season/n_t/
        measured_forcing and the scalar Cwo are ignored for the driver series.
    K_surf : root surface/plaque sorption [L/kg] added to the *measured* root BAF
        only (dead-end pool; leaves the ODE untouched). 0 for low-carbon soils.
    cwo_profile : pore-water exposure SHAPE when `drivers` is not supplied --
        'constant' (default; flat at Cwo, so conc == BAF), 'flooded' (analytic
        Freundlich dilution+leaching; short chains leach, long chains buffered;
        no HYDRUS engine needed) or 'hydrus' (the shape from a real HYDRUS-1D run;
        needs the built engine). The shape is normalised to season-mean == Cwo, so
        Cwo stays the AVERAGE exposure. `cwo_kw` passes extra options to
        `cwo_profile_series` (flood_fraction, k_leach, C_total, ...). Ignored when
        `drivers` are given (those carry their own Cwo).
    uptake : which ROOT-ENTRY mechanism carries the anion past GHK exclusion --
        'carrier' (default; the fitted Michaelis-Menten carrier, the incumbent
        every published number runs on) or 'bypass' (carrier off, one global
        passive apoplastic conductance g_apo). The two are INDISTINGUISHABLE on
        the only dataset that has been asked (re-fitted carrier 1.078 vs bypass
        1.093 log10 RMSE, bootstrap 0.368), so the carrier is the default by default, not by
        evidence -- see UPTAKE_MODES and validation/carrier_vs_bypass.py.
    vmax_scale, g_apo : explicit overrides of the carrier capacity multiplier /
        apoplastic conductance. Either one given wins over `uptake`, so a scan
        can move one term while the mode sets the other.
    km_scale : multiplier on the carrier half-saturation Km (both influx and
        efflux). 1.0 (default) is the shipped Km. This is the knob a DOSE series
        acts on -- Km is what decides WHERE the carrier stops being linear, and
        `validation/dose_series_carrier.py` uses it to bound Km from below.
    shoot : "4pool" or "redistributed" (None -> `DEFAULT_SHOOT`); see `SHOOTS`. The
        redistributed shoot takes its own architecture levers N, stem_transp_frac,
        lam_grain and retention (ignored by the 4-pool core) and splits the xylem by
        lam_grain, so `s_grain_dry` (the 4-pool's grain surface) does not apply to it.
    lipid_loading : None -> `DEFAULT_LIPID_LOADING`. f_xy_source="fit" is a free-anion fit
        and raises when combined with lipid loading (the lipid model supplies its own f_xy).
    grain_share : the redistributed shoot's routing of the residual xylem stream to the grain,
        "gate" or "size" (None -> `DEFAULT_GRAIN_SHARE`); see `GRAIN_SHARES`. The 4-pool core
        splits by dry matter already and refuses a non-default value.

    Returns a dict with t, per-compartment conc & BAF time series (FRESH weight; the green
    leaf), finals, straw (with attached dead blades), the dry-basis harvest values
    (`conc_dw`, `baf_final_dw` -- compare these with dw-reported data), the dead-blade pool,
    the driver series actually used (Cwo, Qtp, M), B_k, and the effective params.
    """
    if uptake not in UPTAKE_MODES:
        raise ValueError(f"unknown uptake {uptake!r}; expected one of {sorted(UPTAKE_MODES)}")
    shoot = DEFAULT_SHOOT if shoot is None else shoot
    if shoot not in SHOOTS:
        raise ValueError(f"unknown shoot {shoot!r}; expected one of {SHOOTS}")
    lipid_loading = DEFAULT_LIPID_LOADING if lipid_loading is None else bool(lipid_loading)
    if lipid_loading and f_xy_source == "fit":
        raise ValueError("f_xy_source='fit' is a free-anion fit; pass lipid_loading=False "
                         "(the lipid model supplies its own f_xy)")
    grain_share = DEFAULT_GRAIN_SHARE if grain_share is None else grain_share
    if grain_share not in GRAIN_SHARES:
        raise ValueError(f"unknown grain_share {grain_share!r}; expected one of {GRAIN_SHARES}")
    if shoot == "4pool" and grain_share != DEFAULT_GRAIN_SHARE:
        raise ValueError("grain_share applies to the redistributed shoot; the 4-pool core splits "
                         "the xylem by dry matter")
    if shoot == "redistributed" and s_grain_dry is not None:
        raise ValueError("s_grain_dry is the 4-pool's grain xylem surface; the redistributed "
                         "shoot splits the stream by lam_grain")
    _up = UPTAKE_MODES[uptake]
    vmax_scale = _up["vmax_scale"] if vmax_scale is None else float(vmax_scale)
    g_apo = _up["g_apo"] if g_apo is None else float(g_apo)
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
        M_basis = drivers.get("mass_basis", "fresh")
    else:
        t, Cwo_series, Qtp, M, leaf_loss, M_basis = _default_drivers(
            season, n_t, Cwo, measured_forcing, biomass, scenario, scenario_kw)
        season = float(t[-1])
        if cwo_profile != "constant":
            Cwo_series = cwo_profile_series(
                t, level=Cwo, profile=cwo_profile, n_C=c["n_C"], group=c["group"],
                congener=congener,
                scenario=scenario if measured_forcing else "generic", **(cwo_kw or {}))
    inputs = PlantInputs(t=t, Cwo=Cwo_series, Qtp=Qtp, M=M, leaf_loss=leaf_loss,
                         M_basis=M_basis)

    f_xy_def, L_Ph_def, kappa_d_def, g_xy_def, g_ph_def = _transport_defaults(
        c, f_xy_source, lipid_loading)
    f_xy = float(f_xy_override) if f_xy_override is not None else float(f_xy_def)
    g_xy = float(g_xy_override) if g_xy_override is not None else float(g_xy_def)
    g_ph = float(g_ph_override) if g_ph_override is not None else float(g_ph_def)
    kappa_d = float(kappa_d_override) if kappa_d_override is not None else float(kappa_d_def)
    L_Ph = float(L_Ph_override) if L_Ph_override is not None else float(L_Ph_def)
    # the redistributed shoot binds every compartment with the STEM cell-wall K (its own
    # convention, as scored); the 4-pool core with the root's
    kcw_organ = "stem" if shoot == "redistributed" else "root"
    K_cw = c["K_cw_wholecw_Lkg"].get(kcw_organ, c["K_cw_wholecw_Lkg"]["root"])
    cmpd = Compound(name=congener, K_prot=c["K_prot_Lkg"], K_PL=c["K_PL_Lkg"],
                    K_cw=K_cw, kappa_d=kappa_d,
                    Vmax_in=_CARR["Vmax_in"] * vmax_scale,
                    Km_in=_CARR["Km_in"] * float(km_scale),
                    Vmax_out=_CARR["Vmax_out"] * vmax_scale,
                    Km_out=_CARR["Km_out"] * float(km_scale),
                    L_Ph=L_Ph, f_xy=f_xy, g_xy=g_xy, g_ph=g_ph, K_surf=float(K_surf),
                    g_apo=float(g_apo))
    if shoot == "redistributed":
        res = _run_redistributed(
            congener, cmpd, t, Cwo_series, Qtp, M, leaf_loss, season, E_m_mV,
            N=N, stem_transp_frac=stem_transp_frac, lam_grain=lam_grain,
            retention=retention, dead_leaf_attached=dead_leaf_attached,
            phi=phi, T_C_Ph=T_C_Ph, grain_share=grain_share)
        res["params"].update(
            K_PL=c["K_PL_Lkg"], K_prot=c["K_prot_Lkg"], K_cw=K_cw, K_surf=float(K_surf),
            n_C=c["n_C"], group=c["group"], uptake=uptake, vmax_scale=float(vmax_scale),
            g_apo=float(g_apo), Km_in=float(_CARR["Km_in"] * km_scale), mass_basis=M_basis,
            scenario=(None if drivers is not None or not measured_forcing else scenario),
            shoot=shoot, lipid_loading=lipid_loading)
        res["root_baf_total"] = float(res["baf_final"]["root"] + cmpd.K_surf)
        return res
    comps = _compartments()
    env = Environment(E=E_m_mV / 1000.0)
    model = RiceUptakeModel(env=env, cmpd=cmpd, comps=comps, inputs=inputs, phi=float(phi),
                            T_C_Ph=float(T_C_Ph), s_grain_dry=s_grain_dry)
    sol = model.solve(t)
    C = sol.y                                            # (4, n_t)
    conc = {k: C[i] for i, k in enumerate(TISSUES)}
    obs = _harvest(t, conc, M, comps, inputs.M_basis, leaf_loss, dead_leaf_attached, sol, LEAF)
    straw = obs["straw_fw"]
    B = binding_factors(comps, cmpd)
    # time-resolved BAF normaliser (per-step pore water); guard zeros
    cw = np.where(Cwo_series > 0, Cwo_series, np.nan)
    cwo_ref = float(Cwo_series[-1]) if Cwo_series[-1] > 0 else (
        float(np.nanmax(Cwo_series)) if np.nanmax(Cwo_series) > 0 else 1.0)
    baf_series = {k: C[i] / cw for i, k in enumerate(TISSUES)}
    root_baf_total = float(C[ROOT, -1] / cwo_ref + cmpd.K_surf)
    return dict(
        t=t, congener=congener, success=bool(sol.success), season=season,
        conc=conc,
        straw=straw,
        conc_dw=obs["conc_dw"], baf_final_dw={k: float(v[-1] / cwo_ref)
                                              for k, v in obs["conc_dw"].items()},
        dead_leaf=dict(burden=obs["dead_burden"], dry=obs["dead_dry"], attached=obs["attached"]),
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
                    n_C=c["n_C"], group=c["group"],
                    uptake=uptake, vmax_scale=float(vmax_scale), g_apo=float(g_apo),
                    Km_in=float(_CARR["Km_in"] * km_scale),
                    mass_basis=M_basis,
                    scenario=(None if drivers is not None or not measured_forcing else scenario),
                    shoot=shoot, lipid_loading=lipid_loading),
    )


def apportionment(congener="PFOA", Cwo=1.0, E_m_mV=-120.0, f_xy_source="recommended",
                  f_xy_override=None, L_Ph_override=None, kappa_d_override=None,
                  lipid_loading=None, g_xy_override=None, g_ph_override=None,
                  season=None, n_t=241, measured_forcing=True, biomass="growth_rice",
                  scenario=DEFAULT_SCENARIO, shoot=None, grain_share=None):
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
    ``shoot`` (None -> DEFAULT_SHOOT) picks the shoot model; the redistributed shoot uses
    its own flux map (`apportionment.apportion_redistributed`), with the same root/grain keys.
    ``lipid_loading=None`` -> DEFAULT_LIPID_LOADING, as in `simulate`.
    """
    import apportionment as ap
    if congener not in _CONG:
        raise KeyError(f"unknown congener {congener!r}; known: {CONGENERS}")
    shoot = DEFAULT_SHOOT if shoot is None else shoot
    if shoot not in SHOOTS:
        raise ValueError(f"unknown shoot {shoot!r}; expected one of {SHOOTS}")
    lipid_loading = DEFAULT_LIPID_LOADING if lipid_loading is None else bool(lipid_loading)
    t, Cwo_series, Qtp, M, leaf_loss, M_basis = _default_drivers(
        season, n_t, Cwo, measured_forcing, biomass, scenario=scenario)
    if shoot == "redistributed":
        cmpd = _compound_for(congener, f_xy_source, f_xy_override, L_Ph_override,
                             kappa_d_override, lipid_loading, g_xy_override, g_ph_override,
                             K_cw_organ="stem")
        run = _run_redistributed(congener, cmpd, t, Cwo_series, Qtp, M, leaf_loss,
                                 float(t[-1]), E_m_mV, grain_share=grain_share, _keep_model=True)
        res = ap.apportion_redistributed(run["_model"], run["_sol"])
        res.update(congener=congener, success=bool(run["success"]), shoot=shoot,
                   params=dict(f_xy=float(cmpd.f_xy), L_Ph=float(cmpd.L_Ph),
                               n_C=_CONG[congener]["n_C"], group=_CONG[congener]["group"]))
        return res
    inputs = PlantInputs(t=t, Cwo=Cwo_series, Qtp=Qtp, M=M, leaf_loss=leaf_loss,
                         M_basis=M_basis)
    cmpd = _compound_for(congener, f_xy_source, f_xy_override, L_Ph_override,
                         kappa_d_override, lipid_loading, g_xy_override, g_ph_override,
                         K_cw_organ="root")
    comps = _compartments()
    env = Environment(E=E_m_mV / 1000.0)
    model = RiceUptakeModel(env=env, cmpd=cmpd, comps=comps, inputs=inputs)
    sol = model.solve(t)
    res = ap.apportion(model, sol)
    res.update(congener=congener, success=bool(sol.success), shoot=shoot,
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
                  "f_xy_recommended": cmpd.f_xy}   # no fit: L_Ph/kappa_d -> fitted geomean
        res = simulate(cmpd.name, record=record, f_xy_override=fxy_ov, **kw)
    res["descriptors"] = d
    res["provisional"] = (d.matched_name is None) or (not d.is_linear)
    return res


F_XY_SOURCES = ("recommended", "fit")


def _fitted_geomean(key):
    vals = [c[key] for c in PARAMS["congeners"] if c.get(key)]
    return float(10 ** np.mean(np.log10(vals)))


def _transport_defaults(c, f_xy_source, lipid_loading):
    """Resolve the (f_xy, L_Ph, kappa_d, g_xy, g_ph) defaults for a congener record.

    f_xy_source 'recommended': the monotone physical f_xy (`f_xy_recommended`) with L_Ph and
        kappa_d fitted to Yamazaki's grain and root GIVEN it (`L_Ph_rec`, `kappa_d_rec`).
    f_xy_source 'fit': (f_xy, L_Ph, kappa_d) fitted per congener to root, straw and grain
        (`*_fit`; saturated). Both from validation/refit_site_year.py on the 2020 site-year.
    lipid_loading: the K_PL-gated lipid model's f_xy and bound terms on the 'recommended'
        L_Ph / kappa_d (the protocol LIPID_LOADING was fitted with).
    A record without a fit (GenX, a SMILES-derived structure) takes the geometric mean of the
    fitted congeners' values -- the same rule the fit uses for PFBS's missing grain."""
    if f_xy_source not in F_XY_SOURCES:
        raise ValueError(f"f_xy_source must be one of {F_XY_SOURCES}, got {f_xy_source!r}")
    sfx = "_fit" if (f_xy_source == "fit" and not lipid_loading) else "_rec"
    L_Ph = c.get("L_Ph" + sfx) or _fitted_geomean("L_Ph" + sfx)
    kappa_d = c.get("kappa_d" + sfx) or _fitted_geomean("kappa_d" + sfx)
    if lipid_loading:
        f_xy, g_xy, g_ph = lipid_loading_conductances(c["n_C"], c["K_PL_Lkg"], c["group"])
        return f_xy, L_Ph, kappa_d, g_xy, g_ph
    if f_xy_source == "fit":
        return (c.get("f_xy_fit") or c["f_xy_recommended"]), L_Ph, kappa_d, 0.0, 0.0
    return c["f_xy_recommended"], L_Ph, kappa_d, 0.0, 0.0


def _compound_for(congener, f_xy_source="recommended", f_xy_override=None,
                  L_Ph_override=None, kappa_d_override=None, lipid_loading=False,
                  g_xy_override=None, g_ph_override=None, K_cw_organ="stem"):
    """Build a Compound for `congener` with the same parameter resolution as
    simulate() (f_xy source / lipid loading / overrides). `K_cw_organ` picks the
    cell-wall K used (the nstem stem segments use the stem cell wall)."""
    c = _CONG[congener]
    f_xy_def, L_Ph_def, kappa_d_def, g_xy_def, g_ph_def = _transport_defaults(
        c, f_xy_source, lipid_loading)
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
                        season=None, n_t=361, N=4, stem_transp_frac=0.45,
                        lam_grain=0.05, retention=0.6, biomass_fn=None,
                        k_seq=0.0, k_rel=0.0, K_cw_organ="stem", drivers=None,
                        dead_leaf_attached=1.0, scenario=DEFAULT_SCENARIO, scenario_kw=None,
                        grain_share=None):
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
    scenario : "yamazaki2020" (DEFAULT) runs on the calibration site-year drivers, exactly as
        `simulate()`; with `biomass_fn` given (or scenario "generic") it runs on forcing_rice
        transpiration and that biomass driver over `season` days.
    biomass_fn : callable (t, season) -> dict{root,stem,leaf,grain} [kg DRY/hill] for the
        generic drivers (default the ORYZA2000 port, `_biomass_fn("oryza2000")`); its dry
        matter becomes fresh organ mass through `organ_mass`.
    k_seq, k_rel : optional TWO-POOL root (mobile + sequestered) rates [1/day].
        k_seq=0 (default) keeps the single root pool. With k_seq>0 the reported `root`
        is mobile+sequestered and the extra keys `root_mobile`/`root_seq`/`seq_fraction`
        are returned. No fitted two-pool parameter set ships (the two-pool explorations
        are archived: validation/archive/README.md).
    K_cw_organ : which tissue's cell-wall K_cw parameterises the binding ('stem'
        default, as before; 'root' matches the two-pool root fit's convention).
    drivers : optional standard driver dict {t, Cwo, Qtp, M} with M as the FOUR
        organ columns [root, stem, leaf, grain] (the stem column is split evenly
        over the N segments), plus optional 'leaf_loss'. Overrides season/n_t/
        biomass_fn -- used to run this model on another script's forcings.
    """
    if congener not in _CONG:
        raise KeyError(f"unknown congener {congener!r}; known: {CONGENERS}")

    if drivers is not None:
        t = np.asarray(drivers["t"], dtype=float)
        season = float(t[-1])
        Qtp = np.asarray(drivers["Qtp"], dtype=float)
        Cwo_series = np.asarray(drivers["Cwo"], dtype=float)
        M4 = np.asarray(drivers["M"], dtype=float)           # (n_t, 4) organ columns
        leaf_loss = drivers.get("leaf_loss")
    elif biomass_fn is None and scenario == "yamazaki2020":
        d = scenario_drivers(scenario, Cwo=Cwo, **(scenario_kw or {}))
        ss.check_season(season, d["t"], "simulate_nstem_leaf")
        t, Qtp, Cwo_series, leaf_loss = d["t"], d["Qtp"], d["Cwo"], d["leaf_loss"]
        season = float(t[-1])
        M4 = d["M"]
    else:
        season = 150.0 if season is None else float(season)
        t = np.linspace(0.0, season, int(n_t))
        Qtp = fr.Q_TP(t, season)
        b = (biomass_fn or _biomass_fn("oryza2000"))(t, season)
        b = om.to_model_mass(b)                              # crop DRY matter -> FRESH mass
        M4 = np.column_stack([b["root"], b["stem"], b["leaf"], b["grain"]])
        Cwo_series = np.full_like(t, float(Cwo))
        leaf_loss = b.get("leaf_death_rate")
    cmpd = _compound_for(congener, f_xy_source, f_xy_override, L_Ph_override,
                         kappa_d_override, lipid_loading, g_xy_override, g_ph_override,
                         K_cw_organ=K_cw_organ)
    return _run_redistributed(congener, cmpd, t, Cwo_series, Qtp, M4, leaf_loss, season,
                              E_m_mV, N=N, stem_transp_frac=stem_transp_frac,
                              lam_grain=lam_grain, retention=retention, k_seq=k_seq,
                              k_rel=k_rel, dead_leaf_attached=dead_leaf_attached,
                              grain_share=grain_share, segments_in_M=True)


def _run_redistributed(congener, cmpd, t, Cwo_series, Qtp, M4, leaf_loss, season, E_m_mV,
                       N=4, stem_transp_frac=0.45, lam_grain=0.05, retention=0.6,
                       k_seq=0.0, k_rel=0.0, dead_leaf_attached=1.0, phi=0.1, T_C_Ph=10.0,
                       grain_share=None, segments_in_M=False, _keep_model=False):
    """Solve the redistributed-shoot model on prepared drivers and a built Compound.

    M4 is the FRESH organ-mass matrix with the four organ columns [root, stem, leaf, grain];
    the stem column is split evenly over the N segments. Shared by `simulate_nstem_leaf`
    and `simulate(shoot="redistributed")`, so both run the identical solve. The returned
    `M` is the segment matrix when `segments_in_M` (simulate_nstem_leaf's contract) and the
    four-organ matrix otherwise (simulate()'s contract; the segments are in `M_segments`)."""
    from pfas_rice_plant_module_nstem_leaf import (
        NStemLeafModel, PlantInputsNL, make_stem_leaf_compartments, split_from_stem_frac)
    M4 = np.asarray(M4, dtype=float)
    M = np.column_stack(
        [np.maximum(M4[:, 0], 1e-9)]
        + [np.maximum(M4[:, 1] / N, 1e-9)] * N
        + [np.maximum(M4[:, 2], 1e-9), np.maximum(M4[:, 3], 1e-9)])
    inputs = PlantInputsNL(t=t, Cwo=Cwo_series, Qtp=Qtp, M=M, leaf_loss=leaf_loss)
    g = _COMP
    _kw = lambda d: dict(theta=d["theta_fw"], f_prot=d["f_prot"], f_PL=d["f_PL"], f_cw=d["f_cw"])
    comps = make_stem_leaf_compartments(
        N, _kw(g["stem"]), _kw(g["root"]), _kw(g["leaf"]), _kw(g["grain_brown"]))
    tau, lam_leaf, lam_grain = split_from_stem_frac(N, stem_transp_frac, lam_grain)
    env = Environment(E=E_m_mV / 1000.0)
    grain_share = DEFAULT_GRAIN_SHARE if grain_share is None else grain_share
    model = NStemLeafModel(env=env, cmpd=cmpd, comps=comps, inputs=inputs, tau=tau,
                           lam_leaf=lam_leaf, lam_grain=lam_grain, retention=retention,
                           phi=float(phi), T_C_Ph=float(T_C_Ph),
                           k_seq=float(k_seq), k_rel=float(k_rel), grain_share=grain_share)
    sol = model.solve(t)
    Y = sol.y                                            # (N+3 [+1 if two-pool], n_t)
    seg = slice(1, N + 1)
    Mseg = M[:, seg]
    stem_conc = np.sum(Y[seg, :] * Mseg.T, axis=0) / np.sum(Mseg.T, axis=0)
    conc = {"root": model.root_total(Y), "stem": stem_conc,
            "leaf": Y[model.LEAF], "grain": Y[model.GRAIN]}
    if model.two_pool:                                   # expose the root split
        conc["root_mobile"] = Y[model.ROOT]
        conc["root_seq"] = Y[model.SEQ]
    cwo_ref = float(Cwo_series[-1]) if Cwo_series[-1] > 0 else (
        float(np.nanmax(Cwo_series)) if np.nanmax(Cwo_series) > 0 else 1.0)
    baf_final = {k: float(v[-1] / cwo_ref) for k, v in conc.items()}
    # straw = stem segments + green leaf + attached dead blades (the harvest operator)
    Mo = np.column_stack([M[:, 0], np.sum(Mseg, axis=1), M[:, model.LEAF], M[:, model.GRAIN]])
    obs = hv.harvest_observables(
        t, conc, Mo, hv.dry_fractions({o: _COMP["grain_brown" if o == "grain" else o]["theta_fw"]
                                       for o in TISSUES}),
        leaf_loss=leaf_loss, attached=dead_leaf_attached,
        C_leaf_fn=lambda tt: sol.sol(tt)[model.LEAF])
    straw = obs["straw_fw"]
    B = binding_factors(comps, cmpd)
    rec = _CONG.get(congener, {})
    return dict(
        t=t, congener=congener, success=bool(sol.success), season=float(season),
        conc=conc, Cwo=Cwo_series, Qtp=Qtp, leaf_loss=leaf_loss,
        straw=straw, straw_baf=float(straw[-1] / cwo_ref), cwo_ref=cwo_ref,
        conc_dw=obs["conc_dw"], baf_final_dw={k: float(v[-1] / cwo_ref)
                                              for k, v in obs["conc_dw"].items()},
        dead_leaf=dict(burden=obs["dead_burden"], dry=obs["dead_dry"], attached=obs["attached"]),
        baf={k: v / cwo_ref for k, v in conc.items()},
        baf_final=baf_final,
        tf_final={k: baf_final[k] / baf_final["root"]
                  for k in ("root", "stem", "leaf", "grain")},
        seq_fraction=(float(Y[model.SEQ, -1] / max(conc["root"][-1], 1e-12))
                      if model.two_pool else 0.0),
        M=(M if segments_in_M else Mo), M_segments=M, N_segments=N, N=float(env.N),
        eN=float(np.exp(env.N)),
        B_k={"root": float(B[0]), "stem": float(B[1]), "leaf": float(B[model.LEAF]),
             "grain": float(B[model.GRAIN])},
        params=dict(f_xy=float(cmpd.f_xy), retention=float(retention),
                    stem_transp_frac=float(stem_transp_frac), lam_leaf=float(lam_leaf),
                    lam_grain=float(lam_grain), n_C=rec.get("n_C"),
                    group=rec.get("group"), L_Ph=float(cmpd.L_Ph),
                    kappa_d=float(cmpd.kappa_d), g_xy=float(cmpd.g_xy),
                    g_ph=float(cmpd.g_ph), k_seq=float(k_seq), k_rel=float(k_rel),
                    grain_share=grain_share),
        **({"_model": model, "_sol": sol} if _keep_model else {}),
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
# NEUTRAL-organic (Briggs/Kow) path -- the DPU base, exposed like the other models
#   src/neutral_dpu.py ; docs/neutral_dpu_validation.md
# Everything else in this API is PFAS: a permanently dissociated anion for which
# the Briggs/Kow partition core explicitly does NOT apply. A neutral compound is
# the SAME 4-compartment ODE with z=0 (no exclusion, no carrier), K_PW from log
# Kow and f_xy = TSCF(log Kow) -- both fixed from OUTSIDE by published QSPRs with
# nothing to tune, which is what makes it the one setting where the DPU BACKBONE
# is testable apart from the fitted PFAS transport.
# ---------------------------------------------------------------------------
def simulate_neutral(log_kow, name="neutral", Cwo=1.0, season=None, n_t=241,
                     MW=float("nan"), K_AW=0.0, kappa_d=20.0, half_life=None,
                     tscf_model="briggs", phloem=False, L_Ph=1.0,
                     air=False, air_kw=None, waters=None, lipids=None,
                     biomass="oryza", measured_forcing=True, drivers=None,
                     tscf=None, lipid_source=None,
                     pKa=None, is_acid=True, pH=6.5, g_apo=0.0,
                     cwo_profile="constant", cwo_kw=None, Koc=None,
                     dead_leaf_attached=1.0, scenario=DEFAULT_SCENARIO, scenario_kw=None):
    """Run the 4-compartment DPU for a NEUTRAL (non-ionised) organic.

    The neutral analogue of `simulate()`: same driver machinery, same result-dict
    shape (root/stem/leaf/grain `conc`/`baf`/`baf_final`/`straw`/`tf_final`), so
    the app, plots and validation can treat it interchangeably. It adds the
    neutral diagnostics -- `K_PW` per tissue, the `TSCF` actually used, the Briggs
    `rcf_briggs`, and the air-exchange summary when that is on.

    NOTHING here is fitted. `K_PW` and `TSCF` both follow from `log_kow` via
    published QSPRs, which is the whole point: pass a log Kow, get a prediction.

    log_kow : the ONE required input (the PFAS path's congener name has no
        analogue -- a neutral compound is described by its lipophilicity).
    half_life : in-planta dissipation half-life [d] -> per-compartment gamma.
        STRONGLY recommended: with no metabolism the leaf is an unbounded terminal
        accumulator, so a run without it is an UPPER BOUND (docs section 3).
    air : plant-air exchange (volatilisation + gaseous uptake), off by default.
        `True` needs `MW` and a non-zero `K_AW`; `air_kw` passes ambient `C_air`,
        `rh`, the specific areas `S`, ... See `src/plant_air.py`.
    phloem : OFF by default, faithfully to the neutral base, which excludes it --
        the phloem is an addition of the ionisable extension. Turning it on is an
        explicit departure, not a default (it drives the small terminal grain hard).
    tscf : override the QSPR with a MEASURED TSCF; `tscf_model` selects between the
        Briggs 1982 bell (default) and the broader Schriever 2020 refit.
    lipid_source : which tissue-lipid reading the compartments are built from --
        "measured" (default, L_root = 1% fw) or "briggs_anchor" (2.47% fw, the
        value Briggs' own regression implies). The choice is genuinely open and
        moves lipophilic compounds by up to 2.5x, so it is a named mode rather
        than a buried constant; see `neutral_dpu.LIPID_SOURCES` and
        docs/neutral_dpu_validation.md section 5. Every published neutral RMSE is
        on "measured".
    waters / lipids : per-tissue fresh-weight water and lipid fractions, for a
        species other than rice, or to override one tissue of the selected
        `lipid_source`. Rice-specific organ lipid remains an open gap: the root
        value is corroborated by measured CEREAL roots (Li 2019), but stem/leaf
        are still Trapp 1994's soybean figures.
    pKa / is_acid / pH : turn this into a WEAK ELECTROLYTE instead of a strictly
        neutral compound. `pKa=None` (default) is the neutral path and is
        bit-identical to before. With a pKa the compound is part neutral molecule
        and part ion at once -- a case the neutral path's `z=0` trick cannot
        express, since one valence has to be either 0 or -1 -- so it routes through
        the `(f_n, f_d)` speciation pair instead (see the speciation block in
        `pfas_rice_plant_module_4pool_surf`). `pH` is the root-zone pH that sets the
        split; `is_acid=False` makes it a weak base, whose ion is a CATION and so
        is ATTRACTED by the inside-negative membrane rather than excluded. Turning
        `phloem=True` on as well activates the leaf->sieve-tube pH ion trap, the
        textbook mechanism for phloem-mobile acidic herbicides.
        TESTED, and BOUNDED by that test (validation/weak_electrolyte_tscf.py, docs
        section 4l): on Schriever 2020's 67 ionisable TSCF rows the DIRECTION is
        supported -- measured transfer rises with f_n (Spearman +0.480) and turning
        speciation on nearly doubles the model's rank correlation (+0.284 -> +0.514)
        -- but the MAGNITUDE is refuted: the influx conductance moves ~1.6e4-fold
        where the measurements move ~3-fold, so it under-delivers (bias -0.164) and
        predicts ~nothing below f_n ~ 1e-3. Use it for the direction of a speciation
        effect, not its size. Still no measured weak-electrolyte RICE dataset (that
        table is 16 species, none of them rice).
    cwo_profile / cwo_kw / Koc : the pore-water time SHAPE, exactly as `simulate()`
        ("constant" default, "flooded", "hydrus"). Sorption is the one compound-specific
        input the shapes need, and the PFAS chain-length Koc QSPR does not apply to a
        neutral, so it comes from `log_kow` via the PROVISIONAL Karickhoff QSPR unless
        an explicit measured `Koc` [L/kg] is given -- prefer the measurement.
    drivers / biomass / measured_forcing / scenario : exactly as `simulate()`.
    """
    import neutral_dpu as ND

    if drivers is not None:
        t = np.asarray(drivers["t"], dtype=float)
        season = float(t[-1])
        Cwo_series = np.asarray(drivers["Cwo"], dtype=float)
        Qtp = np.asarray(drivers["Qtp"], dtype=float)
        M = np.asarray(drivers["M"], dtype=float)
        leaf_loss = drivers.get("leaf_loss")
        M_basis = drivers.get("mass_basis", "fresh")
    else:
        t, Cwo_series, Qtp, M, leaf_loss, M_basis = _default_drivers(
            season, n_t, Cwo, measured_forcing, biomass, scenario, scenario_kw)
        season = float(t[-1])
        if cwo_profile != "constant":
            Cwo_series = cwo_profile_series(
                t, level=Cwo, profile=cwo_profile, log_kow=float(log_kow), Koc=Koc,
                scenario=scenario if measured_forcing else "generic", **(cwo_kw or {}))
    drv = dict(t=t, Cwo=Cwo_series, Qtp=Qtp, M=M, leaf_loss=leaf_loss, mass_basis=M_basis)

    gam = float(np.log(2.0) / half_life) if half_life else 0.0
    lipid_source = lipid_source or ND.DEFAULT_LIPID_SOURCE
    comps = ND.rice_compartments(lipids=lipids, waters=waters,
                                 gammas=dict.fromkeys(TISSUES, gam),
                                 lipid_source=lipid_source)
    cmpd = ND.NeutralCompound(name=name, log_kow=float(log_kow), MW=float(MW),
                              K_AW=float(K_AW), kappa_d=float(kappa_d),
                              tscf=tscf, tscf_model=tscf_model,
                              pKa=pKa, is_acid=bool(is_acid), pH=float(pH),
                              g_apo=float(g_apo))
    res = ND.simulate_neutral(cmpd, drv, comps=comps, phloem=phloem, L_Ph=L_Ph,
                              air=air, air_kw=air_kw, dead_leaf_attached=dead_leaf_attached)

    # re-shape onto the `simulate()` contract (the neutral module normalises the
    # BAF on the final exposure; keep that reference explicit for the caller)
    cwo_ref = float(Cwo_series[-1]) if Cwo_series[-1] > 0 else (
        float(np.nanmax(Cwo_series)) if np.nanmax(Cwo_series) > 0 else 1.0)
    res.update(
        congener=name, season=float(season), M=M, cwo_ref=cwo_ref,
        leaf_loss=leaf_loss, half_life=(float(half_life) if half_life else None),
        params=dict(log_kow=float(log_kow), TSCF=res["TSCF"], gamma=gam,
                    kappa_d=float(kappa_d), K_AW=float(K_AW), MW=float(MW),
                    tscf_model=tscf_model, phloem=bool(phloem), air=res["air"],
                    lipid_source=lipid_source,
                    mass_basis=M_basis,
                    scenario=(None if drivers is not None or not measured_forcing else scenario)),
    )
    return res


# ---------------------------------------------------------------------------
# Driver builders -- the three ways to feed the plant ODE
#   (1) MODEL          : built-in measured/placeholder forcings (simulate(...))
#   (2) HYDRUS / CSV   : drivers_from_arrays / load_driver_csv  (one-way coupling)
#   (3) SOIL INVENTORY : pore_water_from_inventory (Freundlich) + measured Q/M
# ---------------------------------------------------------------------------
def site_season(scenario_kw=None):
    """The calibration site-year's season [days, transplanting to harvest]."""
    return float(scenario_drivers("yamazaki2020", **(scenario_kw or {}))["t"][-1])


def crop_on_grid(t, scenario=DEFAULT_SCENARIO, biomass="oryza2000", season=None):
    """The crop physiology -- transpiration Q_TP [L/d/hill], FRESH organ mass M [kg/hill,
    (n,4)] and the leaf death rate (or None) -- on an arbitrary grid `t` [days after
    transplanting]. "yamazaki2020" interpolates the calibration site-year (held at its
    harvest values beyond its season); "generic" is forcing_rice + the `biomass` driver."""
    t = np.asarray(t, float)
    if scenario == "yamazaki2020":
        d = scenario_drivers(scenario)
        Qtp = np.interp(t, d["t"], d["Qtp"])
        M = np.column_stack([np.interp(t, d["t"], d["M"][:, k]) for k in range(4)])
        ll = np.interp(t, d["t"], d["leaf_loss"], right=0.0)
        return Qtp, M, ll
    if scenario != "generic":
        raise ValueError(f"scenario must be one of {SCENARIOS}, got {scenario!r}")
    season = float(t[-1]) if season is None else float(season)
    b = _biomass_fn(biomass)(t, season)
    return fr.Q_TP(t, season), om.organ_matrix(b), b.get("leaf_death_rate")


def measured_forcing(t, season=None, biomass="oryza2000", scenario=DEFAULT_SCENARIO):
    """Transpiration Q_TP(t) [L/d/hill] and FRESH organ mass M(t) [kg/hill, (n,4)] on grid
    `t` -- the crop physiology of `scenario` (`crop_on_grid`). Reused by the HYDRUS/soil
    modes when the user supplies only a concentration series."""
    Qtp, M, _ = crop_on_grid(t, scenario, biomass, season)
    return Qtp, M


def drivers_from_arrays(t, Cwo, Qtp=None, M=None, season=None, biomass="oryza2000",
                        scenario=DEFAULT_SCENARIO):
    """Assemble a `drivers` dict for `simulate(drivers=...)`.

    Only `Cwo(t)` is PFAS-specific (from HYDRUS / soil inversion / a probe). If `Qtp`/`M`
    are omitted they default to the crop physiology of `scenario` (`crop_on_grid`; the
    calibration site-year by default), so a bare concentration series is enough to run the
    plant model; a defaulted M brings its leaf death rate as `leaf_loss`. A supplied `M` is
    taken as FRESH organ mass as-is.
    """
    t = np.asarray(t, float)
    Cwo = np.asarray(Cwo, float)
    q, m, ll = crop_on_grid(t, scenario, biomass, season)
    leaf_loss = None
    if Qtp is None:
        Qtp = q
    if M is None:
        M, leaf_loss = m, ll
    drv = dict(t=t, Cwo=Cwo, Qtp=np.asarray(Qtp, float), M=np.asarray(M, float),
               mass_basis=MASS_BASIS)
    if leaf_loss is not None:
        drv["leaf_loss"] = np.asarray(leaf_loss, float)
    return drv


# Per-compartment FRESH-weight density [kg/L] (for the growth-table mass<->volume
# bridge / per-volume reporting). The transport ODE itself is mass-based, so these
# default to ~water for the watery organs; rice leaf/culm hold aerenchyma (air spaces)
# so are < 1, and the dense grain is > 1. Editable in the app.
DEFAULT_TISSUE_DENSITY = {"root": 1.0, "stem": 0.30, "leaf": 0.30, "grain": 1.20}


def drivers_from_tables(growth=None, cwo=None, *, growth_units="g/hill", Cwo_const=1.0,
                        n_t=241, season=None, Qtp=None, biomass="oryza2000",
                        plants_per_hill=None, hills_per_m2=25.0, scenario=DEFAULT_SCENARIO):
    """Build a `simulate(drivers=...)` dict from user-entered TABLES.

    growth : dict with 'day' + any of root/stem/leaf/grain -- organ **fresh-weight**
        mass in `growth_units` (g/hill, kg/hill, g/m2, kg/ha, t/ha, ... per
        `measured_biomass.to_kg_per_hill`). None/empty -> M falls back to the crop
        physiology of `scenario` (`crop_on_grid`; so a Cwᵒ table alone still runs).
    cwo : dict with 'day' + 'Cwo' -- the **absolute** pore-water concentration [µg/L].
        None/empty -> a flat `Cwo_const` series.

    Rows are sorted by day and interpolated onto the model grid. Returns the same dict
    shape as `drivers_from_arrays` (t, Cwo, Qtp, M[, leaf_loss]). The per-compartment
    DENSITY is a display/consistency property handled by the caller (the transport ODE
    is fresh-mass based, so it is not needed to integrate).
    """
    import measured_biomass as mb

    def _sorted(table, keys):
        day = np.asarray(table["day"], float)
        order = np.argsort(day)
        return day[order], {k: np.asarray(table[k], float)[order]
                            for k in keys if table.get(k) is not None}

    g_day = g_cols = None
    if growth and growth.get("day") is not None and len(growth["day"]):
        g_day, g_cols = _sorted(growth, mb.ORGANS)
    c_day = c_vals = None
    if cwo and cwo.get("day") is not None and len(cwo["day"]) and cwo.get("Cwo") is not None:
        c_day, _c = _sorted(cwo, ("Cwo",))
        c_vals = _c["Cwo"]

    days = list(g_day) if g_day is not None else []
    days += list(c_day) if c_day is not None else []
    if not days:
        raise ValueError("at least one of the growth / pore-water tables must have rows")
    season = float(np.max(days)) if season is None else float(season)
    t = np.linspace(0.0, season, int(n_t))

    M = None
    if g_cols:
        organs_kg = {o: (mb.to_kg_per_hill(g_cols[o], growth_units, plants_per_hill=plants_per_hill,
                                           hills_per_m2=hills_per_m2) if o in g_cols else None)
                     for o in mb.ORGANS}
        M = mb.biomass_matrix(t, g_day, organs_kg)
    Cwo_arr = (np.interp(t, c_day, c_vals) if c_day is not None
               else np.full_like(t, float(Cwo_const)))
    return drivers_from_arrays(t, Cwo_arr, Qtp=Qtp, M=M, season=season, biomass=biomass,
                               scenario=scenario)


_DRIVER_COLS = ("t", "Cwo", "Qtp", "M_root", "M_stem", "M_leaf", "M_grain")


def load_driver_csv(path_or_buffer, biomass="oryza2000", scenario=DEFAULT_SCENARIO):
    """Load a HYDRUS/external driver table into a `drivers` dict.

    Required columns: t, Cwo, Qtp, M_root, M_stem, M_leaf, M_grain (see
    `docs/visualization_tool.md` for the HYDRUS-1D output mapping). Accepts a path
    or any file-like object. Missing Qtp/M columns fall back to the crop physiology of
    `scenario` (`crop_on_grid`; the calibration site-year by default).
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
    return drivers_from_arrays(t, Cwo, Qtp=Qtp, M=M, biomass=biomass, scenario=scenario)


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


def hydrus_drivers(congener=None, season=None, Cwo_ref=1.0, f_oc=0.02, n_t=241,
                   qtp_from_hydrus=True, biomass="oryza", Kd=None, log_kow=None,
                   scenario=DEFAULT_SCENARIO, **run_kw):
    """Run a real HYDRUS-1D paddy soil model for `congener` and return a `drivers`
    dict (+ the raw PaddyResult) for `simulate(drivers=…)`.

    HYDRUS supplies the congener-dependent pore-water `Cwᵒ(t)` (short chains leach
    under flooding, long chains stay buffered) and the actual root water uptake
    `Q_TP(t)`; the crop side (transpiration driving the soil run, `M(t)`) is `scenario`'s
    -- the calibration site-year by default (`soil_hydrus.inputs_from_hydrus`).
    `Cwᵒ(t)` is normalised to season-mean `Cwo_ref` so the average exposure matches a
    constant-Cwo run. Extra `run_kw` (flood_until, percolation, …) pass through to
    `soil_hydrus.run_paddy_hydrus`.

    For a NEUTRAL organic there is no congener: pass `log_kow` (Kd = Koc(log Kow)·f_oc
    via the PROVISIONAL Karickhoff QSPR, `literature_params.koc_neutral`) or, better,
    an explicit measured `Kd` [L/kg], which overrides everything else."""
    import soil_hydrus as sh
    n_C, group = None, "PFCA"
    if Kd is None and log_kow is not None:
        import literature_params as lp
        Kd = float(lp.koc_neutral(float(log_kow))) * float(f_oc)
    if Kd is None:
        if congener not in _CONG:
            raise KeyError(f"unknown congener {congener!r}; known: {CONGENERS}. "
                           "For a neutral organic pass log_kow= or Kd= instead.")
        n_C, group = _CONG[congener]["n_C"], _CONG[congener]["group"]
    pin, res = sh.inputs_from_hydrus(n_C, group, season=season,
                                     Cwo_ref=Cwo_ref, f_oc=f_oc, n_t=n_t, Kd=Kd,
                                     qtp_from_hydrus=qtp_from_hydrus, biomass=biomass,
                                     scenario=scenario, **run_kw)
    drivers = dict(t=np.asarray(pin.t, float), Cwo=np.asarray(pin.Cwo, float),
                   Qtp=np.asarray(pin.Qtp, float), M=np.asarray(pin.M, float),
                   mass_basis=pin.M_basis)
    if getattr(pin, "leaf_loss", None) is not None:
        drivers["leaf_loss"] = np.asarray(pin.leaf_loss, float)
    return drivers, res


# ---------------------------------------------------------------------------
# Time-varying pore-water exposure shape C_w^o(t) -- the `cwo_profile` switch
# ---------------------------------------------------------------------------
# By default `simulate(Cwo=...)` holds the pore water CONSTANT (so conc == BAF,
# the BAF-reproduction convention). But a real paddy pore water is time-varying:
# short chains LEACH under flooding while long chains stay buffered (the HYDRUS
# coupling quantifies this -- constant-Cwo over-predicts short-chain grain/straw
# ~2-4x). `cwo_profile` makes that time-shape a first-class, congener-resolved
# option WITHOUT changing the default. All shapes are normalised to season-mean
# == `level` (same convention as `inputs_from_hydrus`'s Cwo_ref) so `Cwo` keeps
# its meaning as the AVERAGE exposure and only the temporal shape changes.
_HEADGROUP = {"PFCA": "carboxylate", "PFSA": "sulfonate", "ether": "ether"}


def cwo_profile_series(t, level=1.0, profile="constant", *, n_C=8, group="PFCA",
                       congener=None, flood_fraction=1.0, k_leach=None, f_oc=0.02,
                       C_total=5.0, theta_g_drained=0.35, theta_g_flooded=0.60,
                       log_kow=None, Koc=None, scenario=DEFAULT_SCENARIO, **hydrus_kw):
    """Time-resolved pore-water C_w^o(t) on grid `t`, normalised to mean == `level`.

    profile:
      'constant' -> flat at `level` (the default; conc == BAF).
      'flooded'  -> analytic Freundlich paddy shape (dilution + first-order
                    leaching of the dissolved pool; `soil_paddy_redox_corrected`).
                    The per-congener capacity K_F = Koc(n_C, head_group)*f_oc comes
                    from the C3 Koc QSPR, so SHORT chains (low K_F -> large dissolved
                    fraction) LEACH to a steep decline while LONG chains stay buffered
                    (~flat). Needs NO HYDRUS engine -- pure Python, works everywhere.
                    `k_leach` defaults (None) to the per-congener value CALIBRATED to a
                    real HYDRUS-1D run (`default_k_leach`); pass a float to override.
      'hydrus'   -> the Cwᵒ(t) shape from a real HYDRUS-1D run (needs the built
                    engine; raises if absent) whose crop side is `scenario`'s; under
                    "yamazaki2020" `t` must end at the site's harvest. Only the Cwᵒ shape
                    is taken (Qtp/M stay the caller's), so the three profiles are
                    directly comparable.

    `flood_fraction` of the season is flooded (remainder drained -> terminal
    concentration). Defaults give continuous flooding (=1.0): a monotone leaching
    decline for short chains, ~flat for long chains.

    NEUTRAL organics take the same shapes through `log_kow` (Koc from the PROVISIONAL
    Karickhoff QSPR) or an explicit measured `Koc` [L/kg], either of which replaces the
    PFAS chain-length QSPR. Sorption is the only compound-specific input here, so the
    shape logic itself is unchanged.
    """
    t = np.asarray(t, float)
    if profile == "constant":
        return np.full_like(t, float(level))
    season = float(t[-1]) if t[-1] > 0 else 1.0
    import literature_params as lp
    if Koc is None and log_kow is not None:
        Koc = float(lp.koc_neutral(float(log_kow)))

    if profile == "hydrus":
        if congener is None and Koc is None:
            raise ValueError("cwo_profile='hydrus' needs the congener name, "
                             "or log_kow=/Koc= for a neutral organic")
        hd, _ = hydrus_drivers(congener, season=season, Cwo_ref=float(level),
                               f_oc=f_oc, n_t=len(t), scenario=scenario,
                               Kd=(None if Koc is None else float(Koc) * float(f_oc)),
                               **hydrus_kw)
        return np.interp(t, hd["t"], hd["Cwo"])

    if profile == "flooded":
        if Koc is None:
            Koc = lp.koc(n_C, _HEADGROUP.get(group, "carboxylate"))
        if k_leach is None:                              # per-congener HYDRUS-calibrated default
            k_leach = default_k_leach(congener, n_C, group, Koc=Koc)
        K_F = lp.koc_to_KF(Koc, f_oc, n=1.0)             # linear Kd (stable; n=1)
        flooded = t < (flood_fraction * season)
        Cwo, _ = pore_water_from_inventory(
            t, C_total, K_F=K_F, n=1.0, theta_g=theta_g_drained,
            theta_g_flooded=theta_g_flooded, flooded=flooded, k_leach=k_leach)
        m = float(np.mean(Cwo))
        return Cwo * (float(level) / m) if m > 0 else np.full_like(t, float(level))

    raise ValueError(f"unknown cwo_profile {profile!r}; "
                     "use 'constant', 'flooded', or 'hydrus'")


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


# ---------------------------------------------------------------------------
# Bayesian inverse -- estimate the exposure (pore-water Cwᵒ) from measured tissue
# concentrations, WITH uncertainty (Laplace posterior in log space).
# ---------------------------------------------------------------------------
# upper-bound ODE-solve count of estimate_exposure_bayesian: ref(1) + two
# 3-point parabola passes(6) + the final model_fit solve(1). A non-convex first
# pass early-exits with fewer; used only to drive the progress bar.
_ESTIMATE_MAX_SOLVES = 8


def estimate_exposure_bayesian(congener=None, measured_conc=None, *, sigma_log10=0.15,
                               season=None, E_m_mV=-120.0, f_xy_source="recommended",
                               biomass="oryza2000", measured_forcing=True, n_plot=121,
                               progress=None, log_kow=None, neutral_kw=None, **sim_kw):
    """Bayesian estimate of the soil-water contamination level Cwᵒ [µg/L] that best
    explains measured rice tissue concentrations.

    The inverse of the forward question (contamination -> tissue): given measured
    `measured_conc` (a dict with any of root/straw/grain in µg/kg DRY weight), infer the
    pore-water Cwᵒ. Because root uptake is saturable (GHK + carrier), tissue conc is
    a NONLINEAR, monotone function of Cwᵒ, so this is a real inverse, not a division.

    Method (matching validation/bayesian_inverse_demo.py): a Laplace posterior in
    log10(Cwᵒ) -- find the MAP exposure (1-D bounded MAP fit) under a lognormal
    measurement model (σ = `sigma_log10` in log10 units, a flat prior in log Cwᵒ),
    then a Gaussian posterior whose width comes from the curvature of the
    log-likelihood at the MAP. Only the EXPOSURE level is estimated; transport is
    held at the model defaults (this is the well-posed direction -- separating
    water-uptake vs root->shoot loading needs an independent measurement).

    Extra keyword args (`uptake`, `lipid_loading`, `km_scale`, ...) pass straight to
    `simulate`, so the inverse is run under the SAME mechanism as the forward model
    rather than silently falling back to the shipped one.

    `log_kow` switches the forward model to the NEUTRAL path (`simulate_neutral`), with
    `neutral_kw` carrying its own arguments (half_life, tscf_model, pKa, air, ...) --
    the same inverse question for a non-PFAS compound. `congener` is then only a label.
    The identifiability caveat is unchanged and arguably sharper there: the neutral
    forward model has no fitted transport at all, so what the posterior width reports is
    the exposure's own conditioning, and for a compound with a short in-planta half-life
    it will be wide unless `neutral_kw["half_life"]` is set to something defensible.

    `progress`, if given, is called as progress(done, total) after each ODE solve
    (done = solves so far, total = `_ESTIMATE_MAX_SOLVES`) so a UI can show a live
    bar -- the run is ~8 forward solves (~1-2.5 s each) and otherwise looks frozen.
    Progress-callback exceptions are swallowed (never break the estimate).

    Returns a dict: median (Cwᵒ), ci68/ci95 (low, high), sd_log10, used_tissues,
    measured, model_fit (predicted tissue conc at the median), grid (Cwᵒ, density
    for plotting), sigma_log10, n_eval, note. Raises ValueError if no usable tissue.
    """
    used = {k: float(v) for k, v in (measured_conc or {}).items()
            if k in ("root", "straw", "grain") and v is not None
            and np.isfinite(v) and float(v) > 0}
    if not used:
        raise ValueError("Provide at least one positive measured tissue "
                         "concentration (root, straw and/or grain) in µg/kg.")

    n_eval = [0]

    def conc_at(Cwo):
        n_eval[0] += 1
        if log_kow is not None:                # NEUTRAL organic: no congener to name
            r = simulate_neutral(float(log_kow), name=(congener or "neutral"),
                                 Cwo=float(Cwo), season=season, biomass=biomass,
                                 measured_forcing=measured_forcing,
                                 **{**(neutral_kw or {}), **sim_kw})
        else:
            r = simulate(congener, Cwo=float(Cwo), season=season, E_m_mV=E_m_mV,
                         f_xy_source=f_xy_source, biomass=biomass,
                         measured_forcing=measured_forcing, **sim_kw)
        if progress is not None:
            try:
                progress(min(n_eval[0], _ESTIMATE_MAX_SOLVES), _ESTIMATE_MAX_SOLVES)
            except Exception:                              # noqa: BLE001 (UI cb never breaks the fit)
                pass
        dw = r["conc_dw"]                      # harvest values, DRY weight (as tissue data are reported)
        return {"root": float(dw["root"][-1]), "straw": float(dw["straw"][-1]),
                "grain": float(dw["grain"][-1])}

    def neg_loglike(L):
        c = conc_at(10.0 ** L)
        return 0.5 * sum(((np.log10(used[k]) - np.log10(max(c[k], 1e-12))) / sigma_log10) ** 2
                         for k in used)

    # rough first guess: tissue conc ~ BAF·Cwᵒ near Cwᵒ=1 -> geomean(meas/BAF)
    ref = conc_at(1.0)
    L0 = float(np.mean([np.log10(used[k]) - np.log10(max(ref[k], 1e-12)) for k in used]))

    # Quadratic-fit Laplace: the log-likelihood is smooth & unimodal in L=log10(Cwᵒ),
    # so a parabola through 3 points gives BOTH the MAP (vertex) and the curvature
    # (=> posterior width). The likelihood is mildly asymmetric (saturable uptake), so
    # iterate -- recentre on the vertex and shrink the span until it settles and the
    # span is local enough for an accurate curvature. Same Laplace idea as
    # bayesian_inverse_demo.py at a fraction of the ODE solves.
    # Two passes: a coarse parabola to locate, then a LOCAL one (small span) centred on
    # it for an accurate MAP + curvature. Deterministic (6 ODE solves), which bounds the
    # interactive cost.
    Lc, curv = L0, 0.0
    for d in (0.35, 0.18):
        xs = np.array([Lc - d, Lc, Lc + d])
        ys = np.array([neg_loglike(x) for x in xs])
        A, B, _C = np.polyfit(xs, ys, 2)
        if A <= 1e-9:                              # not convex -> data uninformative
            curv = 0.0
            Lc = float(xs[int(np.argmin(ys))])
            break
        curv = 2.0 * A
        Lc = float(np.clip(-B / (2.0 * A), Lc - d, Lc + d))      # no extrapolation
    Lmap = float(Lc)
    sdL = float(np.sqrt(1.0 / curv)) if curv > 1e-6 else float("inf")

    median = float(10.0 ** Lmap)
    if np.isfinite(sdL):
        ci68 = (float(10.0 ** (Lmap - sdL)), float(10.0 ** (Lmap + sdL)))
        ci95 = (float(10.0 ** (Lmap - 1.96 * sdL)), float(10.0 ** (Lmap + 1.96 * sdL)))
        span = max(sdL, 1e-3)
        Lg = np.linspace(Lmap - 4 * span, Lmap + 4 * span, int(n_plot))
        dens = np.exp(-0.5 * ((Lg - Lmap) / span) ** 2)
    else:                                          # data did not constrain Cwᵒ
        ci68 = ci95 = (float("nan"), float("nan"))
        Lg = np.linspace(Lmap - 2.0, Lmap + 2.0, int(n_plot))
        dens = np.ones_like(Lg)

    return dict(
        congener=congener, median=median, ci68=ci68, ci95=ci95, sd_log10=sdL,
        used_tissues=list(used), measured=used, model_fit=conc_at(median),
        grid=dict(Cwo=(10.0 ** Lg).tolist(), density=(dens / dens.max()).tolist()),
        sigma_log10=float(sigma_log10), n_eval=int(n_eval[0]),
        note=("Estimates the soil-water PFAS level that best explains your "
              "measurements, assuming the model's plant-uptake behaviour. Transport "
              "is held fixed (the exposure level is the well-posed quantity)."))


# Out-of-sample error of the DEFAULT configuration (the redistributed shoot with K_PL-gated
# lipid loading, monotone f_xy with L_Ph/kappa_d_rec, the calibration site-year): log10 RMSE
# on the three independent datasets, none of which entered any fit (validation/
# default_switch.py, which also records the scripts). They did, however, CHOSE this
# configuration (the switch scored the candidates on them), so the band built from them may
# be optimistic. Yamazaki is NOT in the list: the lipid constants were fitted on it, so its
# 0.380 is in-sample.
APRIORI_OOS_SCORES = {
    "Tang 2026 per-organ TF, 0.1 ug/g": 0.640,       # oos_tang_lipid.py, redistributed + lipid
    "Kim 2019 grain BAF excl PFOA (point)": 0.469,   # kim2019_censored.py, nstem_lipid
    "Li 2025 shoot/root TF by stage": 0.358,         # li2025_stage_resolved.py, nstem_lipid
    "Li 2025 grain/root TF by stage": 0.510,         # li2025_stage_resolved.py, nstem_lipid
}
# The band uses their MEDIAN (no single test speaks for the others). This exists so the
# Simple view never shows a precise-looking absolute number without a reminder of how far
# the default model's predictions sit from measurements. It is an HONESTY reminder, NOT a
# calibrated confidence interval. The previous default (4-pool, free anion; FOUR_POOL_FREE)
# scored Yamazaki straw 1.401, Tang 1.714, Kim 2.33, Li 1.26 / 1.68 (median 1.68, x48).
APRIORI_LOG10_RMSE = float(np.median(list(APRIORI_OOS_SCORES.values())))


def uncertainty_factor(log10_rmse=None):
    """Multiplicative (×/÷) uncertainty factor from the a-priori predictive log10 RMSE
    of the default configuration (`APRIORI_LOG10_RMSE`, 10**that). A coarse honesty
    reminder for the general-audience view, NOT a calibrated CI. Pure."""
    r = APRIORI_LOG10_RMSE if log10_rmse is None else float(log10_rmse)
    return float(10.0 ** max(r, 0.0))


def predictive_band(value, log10_rmse=None):
    """A coarse ×/÷ uncertainty band around a forward model estimate.

    From the a-priori predictive error of the default configuration (see
    `APRIORI_OOS_SCORES` / `APRIORI_LOG10_RMSE`). The model's out-of-sample predictions routinely
    differ from measurements by several fold, so a single decimal figure overstates
    the precision -- this returns a band to show alongside the number for a general
    audience. It is an HONESTY band, NOT a calibrated confidence interval. Pure;
    head-less testable.

    Returns dict(value, lo, hi, factor, log10_rmse); a non-finite/negative value
    gives a NaN lo/hi (factor still reported)."""
    f = uncertainty_factor(log10_rmse)
    r = APRIORI_LOG10_RMSE if log10_rmse is None else float(log10_rmse)
    v = float(value)
    if not np.isfinite(v) or v < 0:
        return dict(value=v, lo=float("nan"), hi=float("nan"), factor=f, log10_rmse=r)
    return dict(value=v, lo=v / f, hi=v * f, factor=f, log10_rmse=r)


# ---------------------------------------------------------------------------
# Dietary-intake reference (policy hook) -- EFSA 2020 group TWI
# ---------------------------------------------------------------------------
# EFSA CONTAM Panel 2020 (doi:10.2903/j.efsa.2020.6223) set a GROUP tolerable
# weekly intake of 4.4 ng/kg body weight per week for the SUM of four PFAS:
# PFOA + PFOS + PFNA + PFHxS (immune effect the critical endpoint). Only these
# four are covered by the TWI; every other congener is compared for REFERENCE
# only (flagged by `in_group=False`). There is NO legal maximum level (MRL) for
# PFAS in rice in either the EU (Reg. 2023/915 covers eggs/fish/meat only) or
# Korea, so an intake-based comparison to the health-based TWI is the defensible
# policy anchor -- NOT a food-standard exceedance.
EFSA_TWI_NG_PER_KG_BW_WEEK = 4.4
EFSA_TWI_CONGENERS = ("PFOA", "PFOS", "PFNA", "PFHxS")
# Korean rice intake: KOSIS per-capita rice supply 2023 ~56.4 kg/yr (~154 g/day);
# KNHANES 2016 ~169.6 g/day. A round 150 g/day default, adjustable in the UI.
DEFAULT_RICE_INTAKE_G_PER_DAY = 150.0
# Korean adult reference body weight (adjustable); EFSA uses 70 kg for adults.
DEFAULT_BODY_WEIGHT_KG = 60.0


def intake_fraction(grain_conc_ug_per_kg, *, congener=None,
                    rice_intake_g_day=DEFAULT_RICE_INTAKE_G_PER_DAY,
                    body_weight_kg=DEFAULT_BODY_WEIGHT_KG):
    """Estimated share of the EFSA group TWI from eating rice at `grain_conc`.

    Converts a predicted grain concentration (µg/kg == ng/g) into an estimated
    weekly PFAS intake -- via a daily rice intake and a body weight -- and
    expresses it as a percentage of the EFSA 2020 group tolerable weekly intake
    (4.4 ng/kg bw/week for PFOA+PFOS+PFNA+PFHxS). Pure; head-less testable.

    HONESTY (carry into the UI): the TWI is for the SUM of the four group PFAS,
    so (a) comparing a SINGLE congener is a conservative "this one chemical
    alone" framing, and (b) it is only defined for the four group members --
    `in_group` flags whether `congener` is one of them (None if unspecified).
    This is an intake reference to a health-based guidance value, NOT a
    food-standard (MRL) exceedance -- no PFAS MRL exists for rice.

    Returns dict(grain_conc, rice_intake_g_day, body_weight_kg,
    weekly_intake_ng_per_kg_bw, twi, fraction, percent, in_group, signal).
    `signal` is a coarse traffic light: "green" <10 %, "amber" 10-100 %,
    "red" >=100 % of the TWI (a NaN/negative conc gives signal=None)."""
    c = float(grain_conc_ug_per_kg)                       # µg/kg == ng/g
    intake = float(rice_intake_g_day)
    bw = float(body_weight_kg)
    in_group = (congener in EFSA_TWI_CONGENERS) if congener else None
    if not np.isfinite(c) or c < 0 or bw <= 0:
        return dict(grain_conc=c, rice_intake_g_day=intake, body_weight_kg=bw,
                    weekly_intake_ng_per_kg_bw=float("nan"),
                    twi=EFSA_TWI_NG_PER_KG_BW_WEEK, fraction=float("nan"),
                    percent=float("nan"), in_group=in_group, signal=None)
    weekly = c * intake * 7.0 / bw                        # ng/kg bw / week
    frac = weekly / EFSA_TWI_NG_PER_KG_BW_WEEK
    pct = frac * 100.0
    signal = "green" if pct < 10.0 else ("amber" if pct < 100.0 else "red")
    return dict(grain_conc=c, rice_intake_g_day=intake, body_weight_kg=bw,
                weekly_intake_ng_per_kg_bw=weekly, twi=EFSA_TWI_NG_PER_KG_BW_WEEK,
                fraction=frac, percent=pct, in_group=in_group, signal=signal)


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


# ---------------------------------------------------------------------------
# Export helpers (UI-agnostic; pure string builders -- NOT model math)
# ---------------------------------------------------------------------------
def summary_csv(res, obs=None, measured=None):
    """One-row-per-tissue BAF summary as a CSV string (for a download button).

    Columns: tissue, model_BAF_L_per_kg, final_conc_ug_per_kg, observed_BAF,
    measured_BAF. `obs` is the Yamazaki dict from `observed_baf`; `measured` an
    optional biomonitoring BAF dict. Pure -- no Streamlit, head-less testable."""
    obs = obs or {}
    measured = measured or {}
    baf = {"root": res["baf_final"]["root"], "straw": res["straw_baf"],
           "grain": res["baf_final"]["grain"]}
    conc = {"root": float(res["conc"]["root"][-1]), "straw": float(res["straw"][-1]),
            "grain": float(res["conc"]["grain"][-1])}

    def _f(x):
        return "" if x is None or not np.isfinite(float(x)) else f"{float(x):.6g}"
    lines = ["tissue,model_BAF_L_per_kg,final_conc_ug_per_kg,observed_BAF,measured_BAF"]
    for tis in ("root", "straw", "grain"):
        lines.append(",".join([tis, _f(baf[tis]), _f(conc[tis]),
                               _f(obs.get(tis)), _f(measured.get(tis))]))
    return "\n".join(lines) + "\n"


def timeseries_csv(res):
    """Full driver + tissue time series as a CSV string (for a download button).

    Columns: t, Cwo, Qtp, conc_<tissue>, M_<tissue> for each of the 4 tissues.
    Pure string builder (no pandas/Streamlit) so it is head-less testable."""
    t = np.asarray(res["t"], float)
    Cwo = np.asarray(res["Cwo"], float)
    Qtp = np.asarray(res["Qtp"], float)
    M = np.asarray(res["M"], float)
    cols = ["t", "Cwo", "Qtp"] + [f"conc_{k}" for k in TISSUES] + [f"M_{k}" for k in TISSUES]
    lines = [",".join(cols)]
    for i in range(len(t)):
        row = [f"{t[i]:.6g}", f"{Cwo[i]:.6g}", f"{Qtp[i]:.6g}"]
        row += [f"{float(res['conc'][k][i]):.6g}" for k in TISSUES]
        row += [f"{M[i, j]:.6g}" for j in range(len(TISSUES))]
        lines.append(",".join(row))
    return "\n".join(lines) + "\n"


def summary_report_md(res, *, congener=None, obs=None, intake=None,
                      scenario=None, lang="ko"):
    """One-page plain-language run summary as a Markdown string (for a handout).

    Pure string builder (no Streamlit/plotly) so it is head-less testable. `intake`
    is an `intake_fraction` dict (optional, adds the EFSA-TWI line); `obs` an
    `observed_baf` dict (optional); `scenario` a short label. `lang="ko"` (default)
    renders Korean, "en" English. The numbers carry the coarse a-priori band."""
    ko = lang != "en"
    cong = congener or res.get("congener", "")
    root_c = float(res["conc"]["root"][-1])
    straw_c = float(res["straw"][-1])
    grain_c = float(res["conc"]["grain"][-1])
    fold = uncertainty_factor()
    rows = [("뿌리" if ko else "roots", root_c), ("짚(줄기+잎)" if ko else "straw", straw_c),
            ("낟알(먹는 쌀)" if ko else "grain", grain_c)]
    tops = {"뿌리" if ko else "roots": root_c, "짚" if ko else "straw": straw_c,
            "낟알" if ko else "grain": grain_c}
    where = max(tops, key=tops.get)
    L = []
    if ko:
        L.append(f"# PFAS–벼 흡수 모델 요약 — {cong}")
        if scenario:
            L.append(f"**시나리오:** {scenario}")
        L.append("")
        L.append("| 부위 | 예측 농도 (µg/kg) | 대략 범위 |")
        L.append("|---|---|---|")
        for name, v in rows:
            b = predictive_band(v)
            L.append(f"| {name} | {v:.2g} | {b['lo']:.2g}–{b['hi']:.2g} |")
        L.append("")
        L.append(f"- 대부분의 PFAS는 **{where}**에 남습니다.")
        if intake is not None and np.isfinite(intake.get("percent", float("nan"))):
            grp = "" if intake.get("in_group") else " (참고 — EFSA 4종 합에 미포함 물질)"
            L.append(f"- 이 쌀을 하루 {intake['rice_intake_g_day']:.0f} g 섭취 시, PFAS 주간 섭취량은 "
                     f"**EFSA 안전기준(TWI)의 약 {intake['percent']:.0f}%**{grp}.")
        L.append("")
        L.append("**한계 (반드시 함께 읽어주세요)**")
        L.append(f"- 절대 수치는 **대략적 예측**으로 실측과 약 ×{fold:.0f}까지 차이날 수 있습니다.")
        L.append("- 쌀에 대한 법적 PFAS 기준(MRL)은 EU·한국 모두 없으며, EFSA 비교는 건강기반 "
                 "섭취기준 환산 **참고치**입니다.")
        L.append("- 정책 스크리닝·교육용 예시이며, 개별 지점의 규제 판단 근거로 단독 사용은 부적절합니다.")
        L.append("")
        L.append("출처: EFSA 2020 그룹 TWI 4.4 ng/kg 체중/주 (doi:10.2903/j.efsa.2020.6223) · "
                 "쌀 섭취량 KOSIS/국민건강영양조사 · 실측 비교 Yamazaki et al. 2023.")
    else:
        L.append(f"# PFAS–Rice uptake summary — {cong}")
        if scenario:
            L.append(f"**Scenario:** {scenario}")
        L.append("")
        L.append("| Tissue | Predicted conc (µg/kg) | Rough range |")
        L.append("|---|---|---|")
        for name, v in rows:
            b = predictive_band(v)
            L.append(f"| {name} | {v:.2g} | {b['lo']:.2g}–{b['hi']:.2g} |")
        L.append("")
        L.append(f"- Most of the PFAS stays in the **{where}**.")
        if intake is not None and np.isfinite(intake.get("percent", float("nan"))):
            grp = "" if intake.get("in_group") else " (reference — not in the EFSA group of four)"
            L.append(f"- Eating {intake['rice_intake_g_day']:.0f} g rice/day, the weekly intake is "
                     f"**~{intake['percent']:.0f}% of the EFSA TWI**{grp}.")
        L.append("")
        L.append("**Limits (read alongside)**")
        L.append(f"- Absolute figures are coarse predictions — up to ~{fold:.0f}× from measurements.")
        L.append("- No legal PFAS limit (MRL) exists for rice (EU/Korea); the EFSA comparison is an "
                 "intake reference to a health-based value.")
        L.append("- Screening/education only; not a single-site regulatory basis.")
        L.append("")
        L.append("Sources: EFSA 2020 group TWI 4.4 ng/kg bw/week (doi:10.2903/j.efsa.2020.6223); "
                 "rice intake KOSIS/KNHANES; observed comparison Yamazaki et al. 2023.")
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    r = simulate("PFOA")
    print("PFOA:", {k: round(v, 3) for k, v in r["baf_final"].items()},
          "straw_baf", round(r["straw_baf"], 3), "B_k", {k: round(v, 2) for k, v in r["B_k"].items()})
