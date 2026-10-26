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
_COMP = PARAMS["tissue_composition_recommended"]
_CONG = {c["name"]: c for c in PARAMS["congeners"]}

CONGENERS = [c["name"] for c in PARAMS["congeners"]]            # 12, ordered
TISSUES = ("root", "stem", "leaf", "grain")
# built-in drivers: the calibration site-year (default) or the generic climatology drivers
SCENARIOS = ("yamazaki2020", "generic")
DEFAULT_SCENARIO = "yamazaki2020"

# ---------------------------------------------------------------------------
# ROOT ENTRY -- the linear equilibrium form the data identify (docs/PREREG_sourced_values.md 5.2).
#
# The GHK + Michaelis-Menten carrier this replaced had Km/Vmax with no source; no rice (or
# any readable) PFAS Km/Vmax exists, and every concentration series found is linear at
# environmental levels (Felizeter 2012/2014, Blaine 2014; Tang's dose series here needed
# Km >= 150-500 ug/L). The root's exchange with the pore water is fast against xylem export
# and growth, so a root BAF constrains only the equilibrium ratio r_eq; the rate LAM_EQ is
# numerical. The carrier-vs-bypass and dose-series analyses that compared root-entry forms
# on the old structure are in validation/archive/.
# ---------------------------------------------------------------------------
# SHOOT MODEL -- which shoot structure carries the xylem-delivered solute.
#   "4pool"         the single-straw core (`RiceUptakeModel`): one stem pool (a mixer), the
#                   leaf/grain xylem split by the panicle's sourced transpiration share.
#   "redistributed" N stem segments + leaf, each organ keeping the solute its own transpiration
#                   deposits (`simulate_nstem_leaf`; retention 1, the sourced stem and panicle
#                   shares).
# DEFAULT = the redistributed shoot with lipid loading; re-decided after the
# amendment-1 re-fit by the rule in docs/PREREG_sourced_values.md 5.3.
SHOOTS = ("4pool", "redistributed")
DEFAULT_SHOOT = "redistributed"
DEFAULT_LIPID_LOADING = True
# the configuration every calibration record and 4-pool mechanism analysis was run on
FOUR_POOL_FREE = dict(shoot="4pool", lipid_loading=False)

# --- inputs known only as a SOURCED RANGE (docs/PREREG_sourced_values.md 5.1, 5.4) ---------
# The model runs at the midpoint; the interval the ends produce is carried with the output
# (`INPUT_RANGE_BAND`). Sources: docs/literature_db/sourced_inputs.csv.
RANGE_INPUTS = {
    "stem_transp_frac": (0.22, 0.36),     # sheath + culm share of transpiration (Sugimoto 1973, 1975)
    "lam_grain_heading": (0.09, 0.33),    # panicle share at heading (Sugimoto; Roberts & Batchelor 1983;
                                          #   Takahashi & Kurata 2007)
    "lam_grain_maturity": (0.0, 0.14),    # panicle share at maturity (Sugimoto 1973 Table 2)
    "phi": (0.02, 0.10),                  # phloem/xylem volume flow, daytime (Windt et al. 2006)
    "T_C_Ph": (4.4, 7.8),                 # phloem sap per kg grain DM [L/kg] (rice sieve-tube sucrose)
    # site-year inputs the paper does not report (src/site_scenario.py has the derivations)
    "sbdur": ss.SBDUR_RANGE,              # nursery days after emergence (ORYZA: whole days)
    "tmpsb": ss.TMPSB_RANGE,              # seedbed cover on the daily maximum [C]
    "nplh": ss.NPLH_RANGE,                # plants per hill
}
# the organ-composition ranges (params/rice_composition_sourced.csv) are in
# parameters.json["tissue_composition_ranges"]; the model runs at their midpoints


def range_mid(name):
    lo, hi = RANGE_INPUTS[name]
    return 0.5 * (lo + hi)


STEM_TRANSP_FRAC = range_mid("stem_transp_frac")
LAM_GRAIN_DVS = (range_mid("lam_grain_heading"), range_mid("lam_grain_maturity"))
PHI = range_mid("phi")
T_C_PH = range_mid("T_C_Ph")
# solute delivered to an evaporating organ stays there (the module's physics); phloem export
# is a separate term -- 1 is derived, not chosen (PREREG 5.2)
RETENTION = 1.0
# amendment 2a (PREREG 7.1, 7.2): organ water follows the development stage; the grain pool is
# brown rice and the panicle's transpired solute stays in the husk + rachis pool. The amendment-1
# forms stay runnable: organ_water="constant", grain_basis="panicle".
ORGAN_WATERS = ("stage", "constant")
GRAIN_BASES = ("brown", "panicle")
DEFAULT_ORGAN_WATER = "stage"
DEFAULT_GRAIN_BASIS = "brown"
_CROP_KEYS = ("organ_water", "grain_basis", "theta_anchors", "panicle_kw")
# LINEAR EQUILIBRIUM root entry lam_eq*(r_eq*Cwo - Cw_root): r_eq is fitted per congener;
# lam_eq is numerical and admissible only because the results do not depend on it once it
# is this large (validation/numerics_convergence.py)
LAM_EQ = 1.0e6

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

# --- lipid-facilitated loading (the default loading form) -------------------------------
# Two superposed xylem/phloem loading terms (docs/fxy_longchain_lipid_exploration.md): the
# free anion loads f_xy x the FREE root concentration, and a term g x the TOTAL concentration
# loads the membrane/lipid-associated pool, so high-binding long chains reach the shoot. The
# full form lets f_xy decline with chain length and gates g by K_PL; the fit keeps only what the
# data identify (amendment 1: one f_xy = a and one g_xy for every chain, both shoots). Its global constants are fitted per shoot model
# with the per-congener (r_eq, L_Ph) by validation/refit_sourced.py and live in
# parameters.json["lipid_loading"][shoot] -- no fitted number is written in this file. A term
# the data do not identify is removed there (stored at its bound: g = 0, b = 0, pfsa_ln = 0;
# gate "none" = every congener loads the bound pool, "linear" = proportional to K_PL).
LIPID_NAMES = ("a", "b", "g_xy", "g_ph", "K_half", "pfsa_ln")


def lipid_globals(shoot=None):
    """The fitted lipid-loading constants of a shoot model (parameters.json)."""
    shoot = DEFAULT_SHOOT if shoot is None else shoot
    try:
        return dict(PARAMS["lipid_loading"][shoot])
    except KeyError:
        raise ValueError(f"no fitted lipid-loading constants for shoot {shoot!r} in "
                         "params/parameters.json -- run validation/refit_sourced.py --write, "
                         "then build_parameters.py") from None


def lipid_loading_conductances(n_C, K_PL, group="PFCA", shoot=None, p=None):
    """(f_xy_free, g_xy, g_ph) of the lipid-loading model (full form; the fitted set may have
    removed terms -- gate "none", b = 0 -- as amendment 1 did for both shoots).

    f_xy_free = a*exp(-b*(n-4))             free-anion loading (declining)
    gate      = K_PL/(K_PL+K_half)          lipid takeover (gate "half", the full form);
                1 (gate "none") or K_PL [L/kg] (gate "linear", g then per L/kg)
    g_xy,g_ph = g_* x gate                  B-independent bound loading
    PFSA carries a head-group factor exp(-pfsa_ln). `p` overrides the shoot's fitted set."""
    p = lipid_globals(shoot) if p is None else p
    sf = np.exp(-p.get("pfsa_ln", 0.0)) if group == "PFSA" else 1.0
    gate = p.get("gate", "half")
    phi = (K_PL / (K_PL + p["K_half"]) if gate == "half" else
           1.0 if gate == "none" else float(K_PL))
    return (p["a"] * np.exp(-p.get("b", 0.0) * (n_C - 4)) * sf,
            p.get("g_xy", 0.0) * phi * sf, p.get("g_ph", 0.0) * phi * sf)


def kcw_by_organ(c, organ=None):
    """{root, stem, leaf, grain} -> K_cw [L/kg cw] for congener record `c`.

    organ=None (the default everywhere) gives each organ its OWN whole-cell-wall K_cw
    (parameters.json `K_cw_wholecw_Lkg`, docs/DELIVERABLE_GAP_A_Kcw.md). A record that
    carries only the root value (a SMILES structure) uses it for every organ. `organ`
    names one organ whose K_cw is used for all (the convention before 2026-10-25: the
    redistributed shoot used the stem's, the 4-pool the root's; kept for reproduction)."""
    kc = c["K_cw_wholecw_Lkg"]
    if organ is not None:
        v = float(kc.get(organ, kc["root"]))
        return {o: v for o in TISSUES}
    return {o: float(kc.get(o, kc["root"])) for o in TISSUES}


def _compartments(kcw=None):
    """[root, stem, leaf, grain] compartments; `kcw` ({organ: K_cw}) binds each organ with
    its own cell-wall K (None -> the compound's single K_cw)."""
    g = _COMP
    k = kcw or {}
    return [Compartment("root",  g["root"]["theta_fw"], g["root"]["f_prot"], g["root"]["f_PL"], g["root"]["f_cw"],
                        K_cw=k.get("root")),
            Compartment("stem",  g["stem"]["theta_fw"], g["stem"]["f_prot"], g["stem"]["f_PL"], g["stem"]["f_cw"],
                        K_cw=k.get("stem")),
            Compartment("leaf",  g["leaf"]["theta_fw"], g["leaf"]["f_prot"], g["leaf"]["f_PL"], g["leaf"]["f_cw"], S=20.0,
                        K_cw=k.get("leaf")),
            Compartment("grain", g["grain_brown"]["theta_fw"], g["grain_brown"]["f_prot"],
                        g["grain_brown"]["f_PL"], g["grain_brown"]["f_cw"], S=2.0, K_cw=k.get("grain"))]


def yamazaki_obs(congener=None):
    """The Yamazaki 2023 observation table (data_obs/obs_baf_Yamazaki.csv, built from the SI by
    docs/literature_db/tools/yamazaki_obs_from_si.py), DRY basis.

    Returns {congener: {"baf_dw": {tissue: value}, "lo": {tissue: value}, "hi": {tissue: value},
    "use": {tissue: bool}, "straw_lo", "straw_hi", "water_ng_L"}} (or one congener's record).
    BAF is tissue pg/g dw over the irrigation-water pg/mL [L/kg dw]; where the SI leaves a value
    open (straw band weights, brown-rice milling share, < MQL cells) baf_dw is the midpoint and
    lo/hi the ends; `use` is False for a cell with no measured value (PFBS grain)."""
    import csv
    path = os.path.join(_ROOT, "data_obs", "obs_baf_Yamazaki.csv")
    out = {}
    with open(path) as f:
        for r in csv.DictReader(line for line in f if not line.startswith("#")):
            d = out.setdefault(r["compound"], {"baf_dw": {}, "lo": {}, "hi": {}, "use": {}})
            t = r["tissue"]
            d["baf_dw"][t] = float(r["baf_dw"])
            d["lo"][t], d["hi"][t] = float(r["lo"]), float(r["hi"])
            d["use"][t] = r["use"] == "1"
            d["water_ng_L"] = float(r["water_ng_L"])
    for d in out.values():
        d["straw_lo"], d["straw_hi"] = d["lo"]["straw"], d["hi"]["straw"]
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
_TANG_TF_ENDPOINT = {"TF_stalk": "stalk", "TF_leaf": "leaf", "TF_endosperm": "endosperm"}
# (model grain pool, Tang organ name, theta key) for the fresh->dry conversion
_TANG_ORGANS = (("stem", "stalk", "stem"), ("leaf", "leaf", "leaf"), ("grain", "endosperm", "grain_brown"))

# What Tang's "endosperm" is compared against. The model's grain pool is the whole
# caryopsis (BROWN rice). Tang does not report whether its endosperm was dehulled only or
# milled (the SI's white-rice moisture standard, Text S3, enters its intake calculation and does
# not describe the analysed sample; docs/rice_sample_processing_KR.md section 4), so "brown"
# stays the default and "white" is the named alternative.
#   "white" re-partitions the grain's free concentration into an endosperm of white-rice
#   composition, assuming bran and endosperm share one free concentration. The composition is
#   the sourced grain_white row of params/rice_composition_sourced.csv (milled-rice protein,
#   NDF and phospholipid, FAO 1993 and Kim et al. 2024; midpoints), carried in parameters.json.
#   Brown and white rice phospholipid are both the measured totals (diacyl + lyso), so the two
#   rows count the same pool.
TANG_ENDOSPERM = ("brown", "white")


def _grain_white_composition():
    return dict(_COMP["grain_white"])


def _dry_binding(comp, c, K_cw):
    """C_dw / Cw for a tissue composition: theta/(1-theta) + sum f*K (L/kg dw)."""
    th = comp["theta_fw"]
    return (th / (1.0 - th) + comp["f_prot"] * c["K_prot_Lkg"]
            + comp["f_PL"] * c["K_PL_Lkg"] + comp["f_cw"] * K_cw)


def endosperm_white_factor(congener, K_cw_organ="grain"):
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


def tang_tf_validation(congener, f_xy_source="recommended",
                       dose="mean", lipid_loading=None, endosperm="brown", scenario=DEFAULT_SCENARIO,
                       shoot=None):
    """Model per-organ transfer factor (DRY weight) vs Tang 2026, for a Tang congener.

    Tang's TF is tissue/root on dry weight; the model's is formed the same way from the
    harvest operator (`conc_dw`: stalk = stem, leaf = green + attached dead blades,
    endosperm = the grain pool). Runs the redistributed-shoot model (`simulate_nstem_leaf`)
    on `scenario` -- by default the calibration site-year, since Tang's site and dates are
    not in the repo (a transfer of drivers as well as parameters). ``lipid_loading`` turns on the lipid loading
    (constants fit on Yamazaki, NOT Tang -> out-of-sample; validation/oos_tang_lipid.py).
    ``endosperm`` -- 'brown' (default; the whole-caryopsis grain) or 'white' (milled rice;
    see TANG_ENDOSPERM). NOTE Tang's roots were rinsed only, so its root -- the TF
    denominator -- may carry surface-sorbed PFAS (no Fe/Mn plaque removal; same confound as
    Li 2025), which would bias every Tang TF low. ``shoot`` -- 'redistributed'
    (`simulate_nstem_leaf`) or '4pool' (the single-straw core), None -> DEFAULT_SHOOT, so
    each shoot model can be scored on the same data; ``lipid_loading`` None ->
    DEFAULT_LIPID_LOADING. Returns None for non-Tang congeners.
    """
    if congener not in TANG_CONGENERS:
        return None
    if endosperm not in TANG_ENDOSPERM:
        raise ValueError(f"unknown endosperm {endosperm!r}; expected one of {TANG_ENDOSPERM}")
    shoot = DEFAULT_SHOOT if shoot is None else shoot
    if shoot not in TANG_SHOOTS:
        raise ValueError(f"unknown shoot {shoot!r}; expected one of {TANG_SHOOTS}")
    lipid_loading = DEFAULT_LIPID_LOADING if lipid_loading is None else bool(lipid_loading)
    if shoot == "4pool":
        r = simulate(congener, Cwo=1.0, f_xy_source=f_xy_source,
                     lipid_loading=lipid_loading, scenario=scenario, shoot="4pool")
    else:
        r = simulate_nstem_leaf(congener, Cwo=1.0, f_xy_source=f_xy_source,
                                lipid_loading=lipid_loading, scenario=scenario)
    dw = r["baf_final_dw"]
    model = {org: dw[mk] / dw["root"] for mk, org, _ in _TANG_ORGANS}
    if endosperm == "white":
        model["endosperm"] *= endosperm_white_factor(congener)     # grain's own K_cw
    return dict(congener=congener, organs=[o for _, o, _ in _TANG_ORGANS],
                model_tf=model, tang_tf=tang_observed_tf(congener, dose),
                f_xy=float(r["params"]["f_xy"]), dose=dose,
                scenario=scenario, endosperm=endosperm,
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
    `unknowns` override the site inputs (sbdur, tmpsb, nh, nplh, co2, heading_doy,
    harvest_is_maturity; defaults in src/site_scenario.py)."""
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
    crop = _crop_kw(unknowns)
    run = site_run(harvest, **unknowns)
    d = ss.drivers_from_run(run, Cwo=Cwo, **crop)
    return {k: (np.array(v, copy=True) if isinstance(v, np.ndarray) else v) for k, v in d.items()}


def _crop_kw(kw):
    """Pop the amendment-2a crop options (organ water, grain basis and their band ends) out of a
    scenario dict, with the defaults filled in."""
    crop = {k: kw.pop(k) for k in _CROP_KEYS if k in kw}
    crop.setdefault("organ_water", DEFAULT_ORGAN_WATER)
    crop.setdefault("grain_basis", DEFAULT_GRAIN_BASIS)
    if crop["organ_water"] not in ORGAN_WATERS:
        raise ValueError(f"organ_water must be one of {ORGAN_WATERS}, got {crop['organ_water']!r}")
    if crop["grain_basis"] not in GRAIN_BASES:
        raise ValueError(f"grain_basis must be one of {GRAIN_BASES}, got {crop['grain_basis']!r}")
    return crop


def _crop_drivers(t, b, crop):
    """(M, extras) from a DRY biomass dict on the grid `t` under the crop options: extras holds
    "theta" (PREREG 7.1) and "husk"/"panicle" (PREREG 7.2) when they apply (need b["dvs"])."""
    import panicle as pn
    has_dvs = b.get("dvs") is not None
    extras = {}
    if crop["grain_basis"] == "brown" and has_dvs:
        b = pn.split_biomass(b, t, **(crop.get("panicle_kw") or {}))
        extras.update(husk=np.asarray(b["husk"], float), panicle=np.asarray(b["panicle"], float))
    stage = crop["organ_water"] == "stage" and has_dvs
    anchors = crop.get("theta_anchors")
    M = om.organ_matrix(b, stage=stage, anchors=anchors)
    if stage:
        extras["theta"] = om.theta_matrix(b, len(t), stage=True, anchors=anchors)
    return M, extras


def _default_drivers(season, n_t, Cwo, measured_forcing, biomass="oryza2000",
                     scenario=DEFAULT_SCENARIO, scenario_kw=None):
    """The built-in drivers as a dict {t, Cwo, Qtp, M, leaf_loss, mass_basis, dvs[, theta, husk,
    panicle]} (the last three under the amendment-2a crop options, `_crop_kw`).

    measured_forcing and scenario "yamazaki2020" (DEFAULT): the calibration site-year
        (`scenario_drivers`); the scenario sets the season (transplanting to harvest) and the
        crop model, so `n_t` and `biomass` do not apply and an explicit `season` other than the
        site's is an error (`site_scenario.check_season`).
    measured_forcing and scenario "generic": `forcing_rice` transpiration and the `biomass`
        driver over `season` days (required; alternative drivers, for sensitivity).
    not measured_forcing: refused (the placeholder logistics were removed)."""
    if scenario not in SCENARIOS:
        raise ValueError(f"scenario must be one of {SCENARIOS}, got {scenario!r}")
    if not measured_forcing:
        raise ValueError("measured_forcing=False selected the illustrative placeholder "
                         "logistics, which had no basis; they were removed (no value without a "
                         "source, docs/PREREG_sourced_values.md)")
    if scenario == "yamazaki2020":
        d = scenario_drivers(scenario, Cwo=Cwo, **(scenario_kw or {}))
        ss.check_season(season, d["t"])
        return d
    if season is None:
        raise ValueError("the generic drivers need an explicit season (days)")
    crop = _crop_kw(dict(scenario_kw or {}))
    season = float(season)
    t = np.linspace(0.0, season, int(n_t))
    Qtp = fr.Q_TP(t, season)
    b = _biomass_fn(biomass)(t, season)
    M, extras = _crop_drivers(t, b, crop)
    leaf_loss = b.get("leaf_death_rate")                # only the senescing ORYZA driver supplies it
    return dict(t=t, Cwo=np.full_like(t, float(Cwo)), Qtp=Qtp, M=M, leaf_loss=leaf_loss,
                mass_basis=MASS_BASIS, dvs=b.get("dvs"), **extras)


def _harvest(t, conc, M, comps, M_basis, leaf_loss, attached, sol, leaf_index, theta=None):
    """The harvest observation operator (src/harvest.py) on a 4-tissue solve: dead blades
    kept in the leaf and straw, dry-basis concentrations for dw-reported data. `theta` (n, 4):
    the organ water that made M fresh, when it follows the stage (PREREG 7.1)."""
    th = ({o: np.asarray(theta, float)[:, i] for i, o in enumerate(TISSUES)} if theta is not None
          else {o: comps[i].theta for i, o in enumerate(TISSUES)})
    dry = hv.dry_fractions(th, M_basis)
    return hv.harvest_observables(
        t, conc, M, dry, leaf_loss=leaf_loss, attached=attached,
        C_leaf_fn=(lambda tt: sol.sol(tt)[leaf_index]) if getattr(sol, "sol", None) else None)


def simulate(congener="PFOA", Cwo=1.0, f_xy_source="recommended",
             f_xy_override=None, L_Ph_override=None, r_eq_override=None,
             lipid_loading=None, g_xy_override=None, g_ph_override=None,
             season=None, n_t=241, measured_forcing=True, biomass="oryza2000",
             drivers=None, K_surf=0.0, record=None,
             cwo_profile="constant", cwo_kw=None,
             dead_leaf_attached=1.0, scenario=DEFAULT_SCENARIO, scenario_kw=None,
             phi=None, T_C_Ph=None, s_grain_dry=None, shoot=None,
             N=4, stem_transp_frac=None, lam_grain=None, retention=RETENTION,
             kcw_organ=None, lam_eq=LAM_EQ, organ_water=None, grain_basis=None,
             husk_phloem=False):
    """Run the 4-compartment ODE for one congener and scenario.

    Parameters
    ----------
    congener : one of CONGENERS.
    Cwo : constant pore-water free concentration C_w^o [ug/L] (1.0 -> conc == BAF).
    f_xy_source : "recommended" (DEFAULT; the monotone a-priori f_xy, with r_eq/L_Ph
        fitted to root and grain given it, so straw is a prediction) or "fit" (f_xy, L_Ph and
        r_eq fitted per congener -- reproduces Yamazaki by construction). See F_XY_SOURCES.
    f_xy_override : if given, use this f_xy instead.
    L_Ph_override, r_eq_override : if given, use these phloem-loading / root equilibrium
        ratio values instead of the fitted ones.
    lam_eq : rate of the linear equilibrium root entry [L/(day kg)]; numerical (`LAM_EQ`).
    scenario : "yamazaki2020" (DEFAULT; the calibration site-year -- JMA Tsukuba 2020, the
        ORYZA2000 port on the paper's dates, ORYZA2000 transpiration; `scenario_drivers`;
        `scenario_kw` overrides its unreported inputs) or "generic" (forcing_rice
        transpiration + the `biomass` driver over `season` days).
    measured_forcing : must be True; False selected logistic placeholders that had no source
        and were removed (it raises).
    biomass, season, n_t : the "generic" scenario's biomass driver and grid
        ('oryza2000' / 'oryza' / 'growth_rice'; see `_biomass_fn`).
    dead_leaf_attached : fraction of dead blades still on the plant at harvest, for the
        harvest operator (`harvest`): 1 (default) keeps them in the leaf and straw -- rice does
        not shed its dead blades (docs/ASSUMPTIONS.md section 2); 0 drops them (a sensitivity).
    phi, T_C_Ph : the phloem stream to the root as a fraction of transpiration and the phloem
        flow per kg of grain dry-matter growth [L/kg]; None -> the midpoints of their sourced
        ranges (`RANGE_INPUTS`).
    s_grain_dry : the legacy transpiring-surface split's grain surface; used only when
        `lam_grain` is a float for the 4-pool core (no source; kept for the neutral path).
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
    shoot : "4pool" or "redistributed" (None -> `DEFAULT_SHOOT`); see `SHOOTS`.
    stem_transp_frac : the sheath + culm share of transpiration (redistributed shoot; None ->
        the sourced midpoint). retention : `RETENTION` (1; derived). N : stem segments
        (results do not depend on it).
    lam_grain : the panicle's share of transpiration. None (default) -> stage-based, 0 before
        heading and the sourced anchors `LAM_GRAIN_DVS` at heading and maturity (needs the
        drivers' 'dvs'); a (heading, maturity) tuple -> other anchors; a float -> a constant
        share without stage (redistributed shoot only).
    lipid_loading : None -> `DEFAULT_LIPID_LOADING`. f_xy_source="fit" is a free-anion fit
        and raises when combined with lipid loading (the lipid model supplies its own f_xy).
    organ_water : "stage" (DEFAULT, PREREG 7.1: root, stem and leaf water follow the drivers'
        DVS, for the fresh mass and the binding alike) or "constant" (amendment 1).
    grain_basis : "brown" (DEFAULT, PREREG 7.2: the grain pool is brown rice, fed by the phloem
        only; the panicle's transpiration deposit stays in a husk + rachis pool, reported as
        `conc_dw["husk"]`, with the rough grain as `conc_dw["rough"]`) or "panicle" (amendment 1:
        the grain is ORYZA's whole panicle). Both need drivers with DVS; without it they do not
        apply. `scenario_kw` may carry their band ends (`theta_anchors`, `panicle_kw`).
    husk_phloem : the husk + rachis dry growth draws phloem sap too (a band end of 7.2).

    Returns a dict with t, per-compartment conc & BAF time series (FRESH weight; the green
    leaf), finals, straw (with attached dead blades), the dry-basis harvest values
    (`conc_dw`, `baf_final_dw` -- compare these with dw-reported data), the dead-blade pool,
    the driver series actually used (Cwo, Qtp, M), B_k, and the effective params.
    """
    shoot = DEFAULT_SHOOT if shoot is None else shoot
    if shoot not in SHOOTS:
        raise ValueError(f"unknown shoot {shoot!r}; expected one of {SHOOTS}")
    lipid_loading = DEFAULT_LIPID_LOADING if lipid_loading is None else bool(lipid_loading)
    if lipid_loading and f_xy_source == "fit":
        raise ValueError("f_xy_source='fit' is a free-anion fit; pass lipid_loading=False "
                         "(the lipid model supplies its own f_xy)")
    if shoot == "redistributed" and s_grain_dry is not None:
        raise ValueError("s_grain_dry is the 4-pool's grain xylem surface; the redistributed "
                         "shoot splits the stream by lam_grain")
    phi = PHI if phi is None else float(phi)
    T_C_Ph = T_C_PH if T_C_Ph is None else float(T_C_Ph)
    stem_transp_frac = STEM_TRANSP_FRAC if stem_transp_frac is None else float(stem_transp_frac)
    if record is not None:
        c = record                       # custom (e.g. SMILES-derived) congener record
    elif congener not in _CONG:
        raise KeyError(f"unknown congener {congener!r}; known: {CONGENERS}")
    else:
        c = _CONG[congener]

    user_drivers = drivers is not None
    skw = dict(scenario_kw or {})
    if organ_water is not None:
        skw["organ_water"] = organ_water
    if grain_basis is not None:
        skw["grain_basis"] = grain_basis
    if drivers is None:
        drivers = _default_drivers(season, n_t, Cwo, measured_forcing, biomass, scenario, skw)
        if cwo_profile != "constant":
            drivers = dict(drivers, Cwo=cwo_profile_series(
                drivers["t"], level=Cwo, profile=cwo_profile, n_C=c["n_C"], group=c["group"],
                congener=congener,
                scenario=scenario if measured_forcing else "generic", **(cwo_kw or {})))
    t = np.asarray(drivers["t"], dtype=float)
    Cwo_series = np.asarray(drivers["Cwo"], dtype=float)
    Qtp = np.asarray(drivers["Qtp"], dtype=float)
    M = np.asarray(drivers["M"], dtype=float)
    leaf_loss = drivers.get("leaf_loss")                # e.g. ORYZA senescence rate
    season = float(t[-1])
    M_basis = drivers.get("mass_basis", "fresh")
    dvs = drivers.get("dvs")
    theta = drivers.get("theta")                        # PREREG 7.1 (None: constant water)
    husk = drivers.get("husk")                          # PREREG 7.2 (None: grain = panicle)
    inputs = PlantInputs(t=t, Cwo=Cwo_series, Qtp=Qtp, M=M, leaf_loss=leaf_loss,
                         M_basis=M_basis, dvs=dvs, theta=theta, M_husk=husk)

    f_xy, L_Ph, r_eq, g_xy, g_ph = _resolve_transport(
        c, f_xy_source, lipid_loading, shoot, f_xy_override, L_Ph_override, r_eq_override,
        g_xy_override, g_ph_override)
    # every compartment binds with its own organ's cell-wall K (before 2026-10-25 the
    # redistributed shoot used the stem's for all and the 4-pool the root's, an unsourced
    # convention; `kcw_organ=` reproduces it)
    kcw = kcw_by_organ(c, kcw_organ)
    K_cw = kcw["root"]
    cmpd = _equilibrium_compound(congener, c, K_cw, f_xy, L_Ph, r_eq, g_xy, g_ph,
                                 K_surf=float(K_surf), lam_eq=lam_eq)
    if shoot == "redistributed":
        res = _run_redistributed(
            congener, cmpd, t, Cwo_series, Qtp, M, leaf_loss, season,
            N=N, stem_transp_frac=stem_transp_frac, lam_grain=lam_grain,
            retention=retention, dead_leaf_attached=dead_leaf_attached,
            phi=phi, T_C_Ph=T_C_Ph, kcw=kcw, dvs=dvs, theta=theta, husk=husk,
            husk_phloem=husk_phloem)
        res["params"].update(
            K_PL=c["K_PL_Lkg"], K_prot=c["K_prot_Lkg"], K_cw=dict(kcw), K_surf=float(K_surf),
            n_C=c["n_C"], group=c["group"], r_eq=r_eq, lam_eq=float(lam_eq), mass_basis=M_basis,
            scenario=(None if user_drivers or not measured_forcing else scenario),
            shoot=shoot, lipid_loading=lipid_loading)
        res["root_baf_total"] = float(res["baf_final"]["root"] + cmpd.K_surf)
        return res
    if isinstance(lam_grain, (int, float)):
        raise ValueError("a constant lam_grain applies to the redistributed shoot; the 4-pool "
                         "core takes None or (heading, maturity)")
    comps = _compartments(kcw)
    env = Environment()          # the linear equilibrium entry carries no membrane potential
    model = RiceUptakeModel(env=env, cmpd=cmpd, comps=comps, inputs=inputs, phi=float(phi),
                            T_C_Ph=float(T_C_Ph), s_grain_dry=s_grain_dry,
                            lam_grain_dvs=(LAM_GRAIN_DVS if lam_grain is None else tuple(lam_grain)),
                            husk_phloem=bool(husk_phloem))
    sol = model.solve(t)
    C = sol.y                                            # (4 [+1 husk], n_t)
    conc = {k: C[i] for i, k in enumerate(TISSUES)}
    obs = _harvest(t, conc, M, comps, inputs.M_basis, leaf_loss, dead_leaf_attached, sol, LEAF,
                   theta=theta)
    if model.husk:
        obs["conc_dw"].update(hv.panicle_observables(
            conc["grain"], M[:, FRUIT], 1.0 - comps[FRUIT].theta, C[model.HUSK], husk))
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
        params=dict(f_xy=f_xy, L_Ph=L_Ph, r_eq=r_eq, lam_eq=float(lam_eq), g_xy=g_xy, g_ph=g_ph,
                    K_PL=c["K_PL_Lkg"], K_prot=c["K_prot_Lkg"],
                    K_cw=dict(kcw), K_surf=float(K_surf),
                    n_C=c["n_C"], group=c["group"],
                    phi=float(phi), T_C_Ph=float(T_C_Ph),
                    mass_basis=M_basis,
                    scenario=(None if user_drivers or not measured_forcing else scenario),
                    shoot=shoot, lipid_loading=lipid_loading,
                    organ_water=("stage" if theta is not None else "constant"),
                    grain_basis=("brown" if husk is not None else "panicle"),
                    husk_phloem=bool(husk_phloem)),
    )


def apportionment(congener="PFOA", Cwo=1.0, f_xy_source="recommended",
                  f_xy_override=None, L_Ph_override=None, r_eq_override=None,
                  lipid_loading=None, g_xy_override=None, g_ph_override=None,
                  season=None, n_t=241, measured_forcing=True, biomass="oryza2000",
                  scenario=DEFAULT_SCENARIO, shoot=None):
    """Source apportionment for one congener -- the PFAS analog of dynamiCROP's Fig. 2.

    Builds the SAME model as ``simulate()`` (identical compound / composition / driver
    resolution) and integrates the inter-compartment MASS fluxes over the season, so the
    cumulative PFAS mass delivered to each compartment is decomposed by transport pathway:
    grain by xylem vs phloem, root by soil-uptake vs phloem recirculation. The headline
    audit number is ``fraction['grain']['phloem_from_leaf']`` -- the model asserts the
    grain is phloem-fed (loading ``L_Ph``), and this quantifies it.

    Unlike dynamiCROP (one application pulse -> split by initial source pool), this model is
    continuously sourced from the pore water, so the decomposition is flux-based and by
    pathway. Returns the ``apportionment.apportion`` dict plus congener/params.
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
    d = _default_drivers(season, n_t, Cwo, measured_forcing, biomass, scenario=scenario)
    t, Cwo_series, Qtp, M, leaf_loss, M_basis, dvs = (d["t"], d["Cwo"], d["Qtp"], d["M"],
                                                      d["leaf_loss"], d["mass_basis"], d.get("dvs"))
    theta, husk = d.get("theta"), d.get("husk")
    if shoot == "redistributed":
        cmpd = _compound_for(congener, f_xy_source, f_xy_override, L_Ph_override,
                             r_eq_override, lipid_loading, g_xy_override, g_ph_override,
                             shoot=shoot)
        run = _run_redistributed(congener, cmpd, t, Cwo_series, Qtp, M, leaf_loss,
                                 float(t[-1]), _keep_model=True, dvs=dvs, theta=theta, husk=husk)
        res = ap.apportion_redistributed(run["_model"], run["_sol"])
        res.update(congener=congener, success=bool(run["success"]), shoot=shoot,
                   params=dict(f_xy=float(cmpd.f_xy), L_Ph=float(cmpd.L_Ph),
                               n_C=_CONG[congener]["n_C"], group=_CONG[congener]["group"]))
        return res
    inputs = PlantInputs(t=t, Cwo=Cwo_series, Qtp=Qtp, M=M, leaf_loss=leaf_loss,
                         M_basis=M_basis, dvs=dvs, theta=theta, M_husk=husk)
    cmpd = _compound_for(congener, f_xy_source, f_xy_override, L_Ph_override,
                         r_eq_override, lipid_loading, g_xy_override, g_ph_override,
                         shoot=shoot)
    comps = _compartments(kcw_by_organ(_CONG[congener]))
    env = Environment()          # the linear equilibrium entry carries no membrane potential
    model = RiceUptakeModel(env=env, cmpd=cmpd, comps=comps, inputs=inputs, phi=PHI,
                            T_C_Ph=T_C_PH, lam_grain_dvs=LAM_GRAIN_DVS)
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


def variant_key(shoot=None, f_xy_source="recommended", lipid_loading=None):
    """The name of a fitted transport variant: '<shoot>_<loading>', loading 'lipid' (the lipid
    model's own f_xy), 'free' (the monotone free-anion f_xy_recommended) or 'fit' (f_xy fitted
    per congener). Each variant's (r_eq, L_Ph[, f_xy]) were fitted with that shoot and loading
    (validation/refit_sourced.py); parameters.json keys are r_eq__<variant> etc."""
    shoot = DEFAULT_SHOOT if shoot is None else shoot
    lipid_loading = DEFAULT_LIPID_LOADING if lipid_loading is None else bool(lipid_loading)
    if f_xy_source not in F_XY_SOURCES:
        raise ValueError(f"f_xy_source must be one of {F_XY_SOURCES}, got {f_xy_source!r}")
    if lipid_loading and f_xy_source == "fit":
        raise ValueError("f_xy_source='fit' is a free-anion fit; pass lipid_loading=False")
    return f"{shoot}_{'lipid' if lipid_loading else ('fit' if f_xy_source == 'fit' else 'free')}"


def _pooled(key):
    """Read-across value for a record without a fit (GenX, a SMILES structure): the median of
    the fitted congeners' values (a fitted L_Ph can be exactly 0, so not a geometric mean)."""
    vals = [c[key] for c in PARAMS["congeners"] if c.get(key) is not None]
    return float(np.median(vals)) if vals else None


def _transport_defaults(c, f_xy_source, lipid_loading, shoot=None):
    """Resolve the (f_xy, L_Ph, r_eq, g_xy, g_ph) defaults of congener record `c` for the
    variant (`variant_key`): r_eq and L_Ph fitted to Yamazaki's root and grain with that
    shoot and loading form ('fit' also fits f_xy to the straw; saturated); f_xy is the lipid
    model's (lipid), the monotone `f_xy_recommended` (free) or the fitted one (fit). A record
    without a fit takes the pooled value (`_pooled`)."""
    v = variant_key(shoot, f_xy_source, lipid_loading)
    L_Ph = c.get(f"L_Ph__{v}")
    r_eq = c.get(f"r_eq__{v}")
    L_Ph = _pooled(f"L_Ph__{v}") if L_Ph is None else L_Ph
    r_eq = _pooled(f"r_eq__{v}") if r_eq is None else r_eq
    if L_Ph is None or r_eq is None:
        raise ValueError(f"no fitted L_Ph__{v} / r_eq__{v} in params/parameters.json -- run "
                         "validation/refit_sourced.py --write, then build_parameters.py")
    if v.endswith("_lipid"):
        f_xy, g_xy, g_ph = lipid_loading_conductances(c["n_C"], c["K_PL_Lkg"], c["group"],
                                                      shoot=v.rsplit("_", 1)[0])
        return f_xy, float(L_Ph), float(r_eq), g_xy, g_ph
    if v.endswith("_fit"):
        f_xy = c.get(f"f_xy__{v}")
        f_xy = _pooled(f"f_xy__{v}") if f_xy is None else f_xy
        return float(f_xy), float(L_Ph), float(r_eq), 0.0, 0.0
    return float(c["f_xy_recommended"]), float(L_Ph), float(r_eq), 0.0, 0.0


def _resolve_transport(c, f_xy_source, lipid_loading, shoot, f_xy=None, L_Ph=None, r_eq=None,
                       g_xy=None, g_ph=None):
    """(f_xy, L_Ph, r_eq, g_xy, g_ph): the overrides given, the variant's fitted values for the
    rest (looked up only when one is missing, so a fit can run before any fit exists)."""
    lipid = DEFAULT_LIPID_LOADING if lipid_loading is None else bool(lipid_loading)
    if not lipid:                         # the free anion has no bound-loading term
        g_xy = 0.0 if g_xy is None else g_xy
        g_ph = 0.0 if g_ph is None else g_ph
    given = (f_xy, L_Ph, r_eq, g_xy, g_ph)
    if any(v is None for v in given):
        dflt = _transport_defaults(c, f_xy_source, lipid_loading, shoot)
        given = tuple(d if v is None else v for v, d in zip(given, dflt))
    return tuple(float(v) for v in given)


def _equilibrium_compound(name, c, K_cw, f_xy, L_Ph, r_eq, g_xy, g_ph, K_surf=0.0,
                          lam_eq=LAM_EQ):
    """The Compound with the linear equilibrium root entry (PREREG 5.2): the GHK and
    Michaelis-Menten terms are off (kappa_d = Vmax = 0) and lam_eq*(r_eq*Cwo - Cw_root) is the
    whole root membrane flux."""
    return Compound(name=name, K_prot=c["K_prot_Lkg"], K_PL=c["K_PL_Lkg"], K_cw=K_cw,
                    kappa_d=0.0, Vmax_in=0.0, Km_in=1.0, Vmax_out=0.0, Km_out=1.0,
                    L_Ph=L_Ph, f_xy=f_xy, g_xy=g_xy, g_ph=g_ph, K_surf=float(K_surf),
                    r_eq=float(r_eq), lam_eq=float(lam_eq))


def _compound_for(congener, f_xy_source="recommended", f_xy_override=None,
                  L_Ph_override=None, r_eq_override=None, lipid_loading=False,
                  g_xy_override=None, g_ph_override=None, K_cw_organ=None, shoot=None):
    """Build a Compound for `congener` with the same parameter resolution as
    simulate() (shoot variant / f_xy source / lipid loading / overrides). The Compound's single K_cw is
    the fallback for a compartment that carries none: the root's by default, or that of
    the organ `K_cw_organ` names. Compartments built by this module carry their own."""
    c = _CONG[congener]
    f_xy, L_Ph, r_eq, g_xy, g_ph = _resolve_transport(
        c, f_xy_source, lipid_loading, shoot, f_xy_override, L_Ph_override, r_eq_override,
        g_xy_override, g_ph_override)
    return _equilibrium_compound(congener, c, kcw_by_organ(c, K_cw_organ)["root"],
                                 f_xy, L_Ph, r_eq, g_xy, g_ph)


def simulate_nstem_leaf(congener="PFOA", Cwo=1.0,
                        f_xy_source="recommended", f_xy_override=None,
                        L_Ph_override=None, r_eq_override=None,
                        lipid_loading=None, g_xy_override=None, g_ph_override=None,
                        season=None, n_t=361, N=4, stem_transp_frac=None,
                        lam_grain=None, retention=RETENTION, biomass_fn=None,
                        k_seq=0.0, k_rel=0.0, K_cw_organ=None, drivers=None,
                        dead_leaf_attached=1.0, scenario=DEFAULT_SCENARIO, scenario_kw=None,
                        phi=None, T_C_Ph=None):
    """Redistributed-shoot uptake run (N stem segments + explicit leaf), the Tang
    2026 over-translocation fix (see `pfas_rice_plant_module_nstem_leaf`).

    Mirrors `simulate()` (same defaults: lipid loading, the sourced shoot shares, the linear
    equilibrium root entry) and returns the SAME tissue keys (root/stem/leaf/grain), where
    `stem` is the mass-weighted mean over the N segments and `M` the segment matrix.

    scenario : "yamazaki2020" (DEFAULT) runs on the calibration site-year drivers, exactly as
        `simulate()`; with `biomass_fn` given (or scenario "generic") it runs on forcing_rice
        transpiration and that biomass driver over `season` days (`season` required).
    biomass_fn : callable (t, season) -> dict{root,stem,leaf,grain[,dvs]} [kg DRY/hill] for the
        generic drivers (default the ORYZA2000 port, `_biomass_fn("oryza2000")`); its dry
        matter becomes fresh organ mass through `organ_mass`.
    k_seq, k_rel : optional TWO-POOL root (mobile + sequestered) rates [1/day]; 0 keeps one pool.
    K_cw_organ : None (default) binds every organ with its own cell-wall K_cw.
    drivers : optional standard driver dict {t, Cwo, Qtp, M[, leaf_loss, dvs]} with M as the
        FOUR organ columns (the stem column is split evenly over the N segments).
    """
    if congener not in _CONG:
        raise KeyError(f"unknown congener {congener!r}; known: {CONGENERS}")
    lipid_loading = DEFAULT_LIPID_LOADING if lipid_loading is None else bool(lipid_loading)
    if drivers is None and biomass_fn is None and scenario == "yamazaki2020":
        drivers = scenario_drivers(scenario, Cwo=Cwo, **(scenario_kw or {}))
        ss.check_season(season, drivers["t"], "simulate_nstem_leaf")
    elif drivers is None:
        if season is None:
            raise ValueError("the generic drivers need an explicit season (days)")
        season = float(season)
        t = np.linspace(0.0, season, int(n_t))
        b = (biomass_fn or _biomass_fn("oryza2000"))(t, season)
        M, extras = _crop_drivers(t, b, _crop_kw(dict(scenario_kw or {})))   # DRY -> FRESH
        drivers = dict(t=t, Qtp=fr.Q_TP(t, season), Cwo=np.full_like(t, float(Cwo)), M=M,
                       leaf_loss=b.get("leaf_death_rate"), dvs=b.get("dvs"), **extras)
    t = np.asarray(drivers["t"], dtype=float)
    season = float(t[-1])
    Qtp = np.asarray(drivers["Qtp"], dtype=float)
    Cwo_series = np.asarray(drivers["Cwo"], dtype=float)
    M4 = np.asarray(drivers["M"], dtype=float)           # (n_t, 4) organ columns
    leaf_loss = drivers.get("leaf_loss")
    dvs = drivers.get("dvs")
    cmpd = _compound_for(congener, f_xy_source, f_xy_override, L_Ph_override,
                         r_eq_override, lipid_loading, g_xy_override, g_ph_override,
                         K_cw_organ=K_cw_organ, shoot="redistributed")
    return _run_redistributed(congener, cmpd, t, Cwo_series, Qtp, M4, leaf_loss, season,
                              N=N, stem_transp_frac=stem_transp_frac,
                              lam_grain=lam_grain, retention=retention, k_seq=k_seq,
                              k_rel=k_rel, dead_leaf_attached=dead_leaf_attached,
                              phi=phi, T_C_Ph=T_C_Ph, segments_in_M=True,
                              kcw=kcw_by_organ(_CONG[congener], K_cw_organ), dvs=dvs,
                              theta=drivers.get("theta"), husk=drivers.get("husk"))


def _run_redistributed(congener, cmpd, t, Cwo_series, Qtp, M4, leaf_loss, season,
                       N=4, stem_transp_frac=None, lam_grain=None, retention=RETENTION,
                       k_seq=0.0, k_rel=0.0, dead_leaf_attached=1.0, phi=None, T_C_Ph=None,
                       segments_in_M=False, _keep_model=False, kcw=None, dvs=None,
                       theta=None, husk=None, husk_phloem=False):
    """Solve the redistributed-shoot model on prepared drivers and a built Compound.

    M4 is the FRESH organ-mass matrix with the four organ columns [root, stem, leaf, grain];
    the stem column is split evenly over the N segments. Shared by `simulate_nstem_leaf`
    and `simulate(shoot="redistributed")`, so both run the identical solve. The returned
    `M` is the segment matrix when `segments_in_M` (simulate_nstem_leaf's contract) and the
    four-organ matrix otherwise (simulate()'s contract; the segments are in `M_segments`)."""
    from pfas_rice_plant_module_nstem_leaf import (
        NStemLeafModel, PlantInputsNL, make_stem_leaf_compartments, split_from_stem_frac,
        default_tau)
    M4 = np.asarray(M4, dtype=float)
    M = np.column_stack(
        [np.maximum(M4[:, 0], om.MASS_FLOOR)]
        + [np.maximum(M4[:, 1] / N, om.MASS_FLOOR)] * N
        + [np.maximum(M4[:, 2], om.MASS_FLOOR), np.maximum(M4[:, 3], om.MASS_FLOOR)])
    th_seg = None
    if theta is not None:                 # organ water by stage, the stem's repeated per segment
        th4 = np.asarray(theta, float)
        th_seg = np.column_stack([th4[:, 0]] + [th4[:, 1]] * N + [th4[:, 2], th4[:, 3]])
    inputs = PlantInputsNL(t=t, Cwo=Cwo_series, Qtp=Qtp, M=M, leaf_loss=leaf_loss, dvs=dvs,
                           theta=th_seg, M_husk=husk)
    g = _COMP
    phi = PHI if phi is None else float(phi)
    T_C_Ph = T_C_PH if T_C_Ph is None else float(T_C_Ph)
    stem_transp_frac = STEM_TRANSP_FRAC if stem_transp_frac is None else float(stem_transp_frac)
    if kcw is None:                       # each organ's own K_cw unless the caller gave a set
        kcw = kcw_by_organ(_CONG[congener]) if congener in _CONG else {}
    _kw = lambda d, o: dict(theta=d["theta_fw"], f_prot=d["f_prot"], f_PL=d["f_PL"],
                            f_cw=d["f_cw"], K_cw=kcw.get(o))
    comps = make_stem_leaf_compartments(
        N, _kw(g["stem"], "stem"), _kw(g["root"], "root"), _kw(g["leaf"], "leaf"),
        _kw(g["grain_brown"], "grain"))
    if isinstance(lam_grain, (int, float)):             # a constant share, no stage
        tau, lam_leaf, lam_g = split_from_stem_frac(N, stem_transp_frac, float(lam_grain))
        split = dict(lam_leaf=lam_leaf, lam_grain=lam_g)
    else:                                                # stage-based (sourced anchors)
        tau = default_tau(N, stem_transp_frac)
        split = dict(lam_grain_dvs=(LAM_GRAIN_DVS if lam_grain is None else tuple(lam_grain)))
    env = Environment()          # the linear equilibrium entry carries no membrane potential
    model = NStemLeafModel(env=env, cmpd=cmpd, comps=comps, inputs=inputs, tau=tau,
                           retention=retention, phi=float(phi), T_C_Ph=float(T_C_Ph),
                           k_seq=float(k_seq), k_rel=float(k_rel),
                           husk_phloem=bool(husk_phloem), **split)
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
    th = ({o: np.asarray(theta, float)[:, i] for i, o in enumerate(TISSUES)} if theta is not None
          else {o: _COMP["grain_brown" if o == "grain" else o]["theta_fw"] for o in TISSUES})
    obs = hv.harvest_observables(
        t, conc, Mo, hv.dry_fractions(th),
        leaf_loss=leaf_loss, attached=dead_leaf_attached,
        C_leaf_fn=lambda tt: sol.sol(tt)[model.LEAF])
    if model.husk:
        obs["conc_dw"].update(hv.panicle_observables(
            conc["grain"], Mo[:, 3], 1.0 - np.asarray(th["grain"], float), Y[model.HUSK], husk))
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
        M=(M if segments_in_M else Mo), M_segments=M, N_segments=N,
        B_k={"root": float(B[0]), "stem": float(B[1]), "leaf": float(B[model.LEAF]),
             "grain": float(B[model.GRAIN])},
        params=dict(f_xy=float(cmpd.f_xy), retention=float(retention),
                    stem_transp_frac=float(stem_transp_frac),
                    lam_grain=(split.get("lam_grain_dvs") or split.get("lam_grain")),
                    n_C=rec.get("n_C"), group=rec.get("group"), L_Ph=float(cmpd.L_Ph),
                    r_eq=float(cmpd.r_eq), lam_eq=float(cmpd.lam_eq), g_xy=float(cmpd.g_xy),
                    g_ph=float(cmpd.g_ph), phi=float(phi), T_C_Ph=float(T_C_Ph),
                    k_seq=float(k_seq), k_rel=float(k_rel),
                    organ_water=("stage" if theta is not None else "constant"),
                    grain_basis=("brown" if husk is not None else "panicle"),
                    husk_phloem=bool(husk_phloem)),
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
def neutral_forms(drivers):
    """A driver dict on the neutral path's forms (PREREG 7.9: constant organ water, whole panicle).
    An amendment-2a dict (fresh M at the staged water, `theta`; brown rice plus `husk`) is turned
    back exactly: dry matter = M (1 - theta), the grain's dry matter + the husk = the panicle, then
    fresh at the constant water. A dict without those keys is returned as it is."""
    if drivers.get("theta") is None and drivers.get("husk") is None:
        return drivers
    M = np.asarray(drivers["M"], float)
    th = (np.asarray(drivers["theta"], float) if drivers.get("theta") is not None
          else np.broadcast_to([om.THETA_FW[o] for o in TISSUES], M.shape))
    dry = M * (1.0 - th)
    if drivers.get("husk") is not None:
        dry = dry.copy()
        dry[:, 3] = dry[:, 3] + np.asarray(drivers["husk"], float)
    Mc = np.maximum(dry / (1.0 - np.array([om.THETA_FW[o] for o in TISSUES])), om.MASS_FLOOR)
    return {k: v for k, v in dict(drivers, M=Mc).items() if k not in ("theta", "husk", "panicle")}


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
        "measured" (default: the sourced dry-basis lipid restated at the organ water,
        PREREG 5.10) or "briggs_anchor" (2.47% fw, the value Briggs' own regression
        implies); see `neutral_dpu.LIPID_SOURCES`.
    waters / lipids : per-tissue fresh-weight water and lipid fractions, for a
        species other than rice, or to override one tissue of the selected
        `lipid_source`. Rice organ lipid is sourced on a dry basis (root: measured
        cereal roots, Li 2019; stem/leaf: rice straw ether extract; grain: brown rice
        crude fat) -- rice total root lipid itself is not measured.
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
        drivers = neutral_forms(drivers)
        t = np.asarray(drivers["t"], dtype=float)
        season = float(t[-1])
        Cwo_series = np.asarray(drivers["Cwo"], dtype=float)
        Qtp = np.asarray(drivers["Qtp"], dtype=float)
        M = np.asarray(drivers["M"], dtype=float)
        leaf_loss = drivers.get("leaf_loss")
        M_basis = drivers.get("mass_basis", "fresh")
    else:
        # the neutral path keeps the constant organ water and the whole panicle until
        # amendment 2b decides them (PREREG 7.9)
        d = _default_drivers(season, n_t, Cwo, measured_forcing, biomass, scenario,
                             dict(scenario_kw or {}, organ_water="constant", grain_basis="panicle"))
        t, Cwo_series, Qtp, M, leaf_loss, M_basis = (d["t"], d["Cwo"], d["Qtp"], d["M"],
                                                     d["leaf_loss"], d["mass_basis"])
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
    (n,4)], the leaf death rate (or None) and the development stage DVS -- on an arbitrary
    grid `t` [days after transplanting]. "yamazaki2020" interpolates the calibration
    site-year (held at its harvest values beyond its season); "generic" is forcing_rice +
    the `biomass` driver."""
    d = _crop_grid(t, scenario, biomass, season, crop=_A1_FORMS)
    return d["Qtp"], d["M"], d["leaf_loss"], d["dvs"]


# a bare (Qtp, M) tuple cannot carry the organ water theta(t) or the husk mass that 2a's fresh M
# is built with, so the tuple helpers return the self-consistent amendment-1 forms (constant water,
# whole panicle); `drivers_from_arrays` returns the full 2a dict
_A1_FORMS = {"organ_water": "constant", "grain_basis": "panicle"}


def _crop_grid(t, scenario=DEFAULT_SCENARIO, biomass="oryza2000", season=None, crop=None):
    """`crop_on_grid` as a dict, with the organ water theta(t) and the husk mass that the fresh
    M was built with (PREREG 7.1/7.2), so a driver set built from it stays self-consistent."""
    t = np.asarray(t, float)
    crop = dict(crop or {})
    if scenario == "yamazaki2020":
        d = scenario_drivers(scenario, **crop)
        lerp = lambda a: (np.column_stack([np.interp(t, d["t"], a[:, k]) for k in range(a.shape[1])])
                          if np.ndim(a) == 2 else np.interp(t, d["t"], a))
        out = dict(Qtp=lerp(d["Qtp"]), M=lerp(d["M"]),
                   leaf_loss=np.interp(t, d["t"], d["leaf_loss"], right=0.0), dvs=lerp(d["dvs"]))
        out.update({k: lerp(d[k]) for k in ("theta", "husk", "panicle") if k in d})
        return out
    if scenario != "generic":
        raise ValueError(f"scenario must be one of {SCENARIOS}, got {scenario!r}")
    season = float(t[-1]) if season is None else float(season)
    b = _biomass_fn(biomass)(t, season)
    M, extras = _crop_drivers(t, b, _crop_kw(crop))
    return dict(Qtp=fr.Q_TP(t, season), M=M, leaf_loss=b.get("leaf_death_rate"),
                dvs=b["dvs"], **extras)


def measured_forcing(t, season=None, biomass="oryza2000", scenario=DEFAULT_SCENARIO):
    """Transpiration Q_TP(t) [L/d/hill] and FRESH organ mass M(t) [kg/hill, (n,4)] on grid
    `t` -- the crop physiology of `scenario` (`crop_on_grid`). Reused by the HYDRUS/soil
    modes when the user supplies only a concentration series."""
    Qtp, M, _, _ = crop_on_grid(t, scenario, biomass, season)
    return Qtp, M


def drivers_from_arrays(t, Cwo, Qtp=None, M=None, season=None, biomass="oryza2000",
                        scenario=DEFAULT_SCENARIO, dvs=None):
    """Assemble a `drivers` dict for `simulate(drivers=...)`.

    Only `Cwo(t)` is PFAS-specific (from HYDRUS / soil inversion / a probe). If `Qtp`/`M`
    are omitted they default to the crop physiology of `scenario` (`crop_on_grid`; the
    calibration site-year by default), so a bare concentration series is enough to run the
    plant model; a defaulted M brings its leaf death rate as `leaf_loss`. A supplied `M` is
    taken as FRESH organ mass as-is. The development stage `dvs` (which sets the panicle's
    share of transpiration) defaults to the scenario crop's.
    """
    t = np.asarray(t, float)
    Cwo = np.asarray(Cwo, float)
    crop = _crop_grid(t, scenario, biomass, season)
    leaf_loss = None
    extras = {}
    if Qtp is None:
        Qtp = crop["Qtp"]
    if M is None:            # the crop's own mass: with its water by stage and its husk pool
        M, leaf_loss = crop["M"], crop["leaf_loss"]
        extras = {k: crop[k] for k in ("theta", "husk") if k in crop}   # what the model reads
    drv = dict(t=t, Cwo=Cwo, Qtp=np.asarray(Qtp, float), M=np.asarray(M, float),
               mass_basis=MASS_BASIS, dvs=np.asarray(crop["dvs"] if dvs is None else dvs, float),
               **extras)
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
                        plants_per_hill=None, hills_per_m2=None, scenario=DEFAULT_SCENARIO,
                        interp="linear"):
    """Build a `simulate(drivers=...)` dict from user-entered TABLES.

    growth : dict with 'day' + any of root/stem/leaf/grain -- organ **fresh-weight**
        mass in `growth_units` (g/hill, kg/hill, g/m2, kg/ha, t/ha, ... per
        `measured_biomass.to_kg_per_hill`). None/empty -> M falls back to the crop
        physiology of `scenario` (`crop_on_grid`; so a Cwᵒ table alone still runs).
    cwo : dict with 'day' + 'Cwo' -- the **absolute** pore-water concentration [µg/L].
        None/empty -> a flat `Cwo_const` series.

    Rows are sorted by day and interpolated onto the model grid (`interp`: "linear" or
    "pchip"; the form between samples is not measured, so both are worth running, PREREG 7.3).
    A falling leaf mass is read as leaf death (`measured_biomass.leaf_loss_rate`; a
    `dead_leaf` column, same units, is used when given), so the green leaf is not concentrated
    by its own senescence. The table's grain column is taken as given and receives the panicle's
    transpiration deposit (no husk column: the amendment-1 grain). Returns the same dict
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
    leaf_loss = None
    if g_cols:
        kg = lambda v: mb.to_kg_per_hill(v, growth_units, plants_per_hill=plants_per_hill,
                                         hills_per_m2=hills_per_m2)
        organs_kg = {o: (kg(g_cols[o]) if o in g_cols else None) for o in mb.ORGANS}
        M = mb.biomass_matrix(t, g_day, organs_kg, interp=interp)
        if "leaf" in g_cols:
            dead = growth.get("dead_leaf")
            dead = (mb.interpolate(t, g_day, kg(np.asarray(dead, float)[np.argsort(
                np.asarray(growth["day"], float))]), interp) if dead is not None else None)
            leaf_loss = mb.leaf_loss_rate(t, M[:, 2], dead)
    Cwo_arr = (np.interp(t, c_day, c_vals) if c_day is not None
               else np.full_like(t, float(Cwo_const)))
    drv = drivers_from_arrays(t, Cwo_arr, Qtp=Qtp, M=M, season=season, biomass=biomass,
                              scenario=scenario)
    if leaf_loss is not None:
        drv["leaf_loss"] = leaf_loss
    return drv


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
    if getattr(pin, "dvs", None) is not None:          # the panicle share follows DVS
        drivers["dvs"] = np.asarray(pin.dvs, float)
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


def estimate_exposure_bayesian(congener=None, measured_conc=None, *, sigma_log10=None,
                               meas_sigma_log10=0.0,
                               season=None, f_xy_source="recommended",
                               biomass="oryza2000", measured_forcing=True, n_plot=121,
                               progress=None, log_kow=None, neutral_kw=None, **sim_kw):
    """Bayesian estimate of the soil-water contamination level Cwᵒ [µg/L] that best
    explains measured rice tissue concentrations.

    The inverse of the forward question (contamination -> tissue): given measured
    `measured_conc` (a dict with any of root/straw/grain in µg/kg DRY weight), infer the
    pore-water Cwᵒ. The root entry is linear, so each tissue concentration is proportional
    to a constant Cwᵒ; with several tissues the estimate weighs their disagreement, and the
    model's own transfer error sets the width.

    Method (matching validation/bayesian_inverse_demo.py): a Laplace posterior in
    log10(Cwᵒ) -- find the MAP exposure (1-D bounded MAP fit) under a lognormal
    measurement model (σ in log10 units, a flat prior in log Cwᵒ),
    then a Gaussian posterior whose width comes from the curvature of the
    log-likelihood at the MAP. Only the EXPOSURE level is estimated; transport is
    held at the model defaults (this is the well-posed direction -- separating
    water-uptake vs root->shoot loading needs an independent measurement).

    Extra keyword args (`shoot`, `lipid_loading`, ...) pass straight to
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

    σ combines the two errors that separate a measurement from the model:
    σ² = meas_sigma_log10² + model², where the model error is the default
    configuration's transfer error (`APRIORI_LOG10_RMSE`, or `NEUTRAL_APRIORI_LOG10_RMSE`
    on the neutral path) and `meas_sigma_log10` is the user's own analytical precision
    (0 -> not given; it is a property of their lab, not of the model). Passing
    `sigma_log10` uses that total instead.
    """
    if sigma_log10 is None:
        model_err = NEUTRAL_APRIORI_LOG10_RMSE if log_kow is not None else APRIORI_LOG10_RMSE
        sigma_log10 = float(np.hypot(float(meas_sigma_log10), model_err))
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
            r = simulate(congener, Cwo=float(Cwo), season=season,
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


# Out-of-sample error of the DEFAULT configuration: log10 RMSE on the independent transfer
# datasets (Tang 2026, Kim 2019, Li 2025), none of which entered any fit, written by
# validation/transfer_sourced.py to params/transfer_scores.json (no score is typed here).
# These datasets also took part in choosing the default (PREREG 5.3), so the band built from
# them may be optimistic. Yamazaki is not in the list: every fit is made on it.
# The band uses their MEDIAN (no single test speaks for the others). It exists so the Simple
# view never shows a precise-looking absolute number without a reminder of how far the
# default model's predictions sit from measurements: an HONESTY reminder, NOT a calibrated
# confidence interval.
def _load_json(name):
    path = os.path.join(_ROOT, "params", name)
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


_TRANSFER = _load_json("transfer_scores.json") or {}
APRIORI_OOS_SCORES = dict(_TRANSFER.get("scores", {}))
APRIORI_LOG10_RMSE = (float(np.median(list(APRIORI_OOS_SCORES.values())))
                      if APRIORI_OOS_SCORES else float("nan"))
# The input-range band (PREREG 5.4): per congener and tissue, the log10 [down, up] envelope
# of the default's Yamazaki observable when each sourced range input moves to its ends with
# the fits held (validation/input_range_band.py -> params/input_range_band.json).
INPUT_RANGE_BAND = (_load_json("input_range_band.json") or {}).get("band", {})

# The NEUTRAL path's a-priori scores (nothing fitted; docs/neutral_dpu_validation.md): the
# model-run tables only (the Schriever TSCF QSPR alone is not a model output), written by
# validation/neutral_scores.py to params/neutral_scores.json (no score is typed here). Their
# median is the neutral path's model error in the inverse below.
_NEUTRAL = _load_json("neutral_scores.json") or {}
NEUTRAL_APRIORI_OOS_SCORES = dict(_NEUTRAL.get("scores", {}))
NEUTRAL_EQUILIBRIUM_SCORES = dict(_NEUTRAL.get("equilibrium_basis", {}))
NEUTRAL_APRIORI_LOG10_RMSE = (float(np.median(list(NEUTRAL_APRIORI_OOS_SCORES.values())))
                              if NEUTRAL_APRIORI_OOS_SCORES else float("nan"))


def input_range_factors(congener, tissue):
    """(down, up) multiplicative factors of the input-range band (`INPUT_RANGE_BAND`, PREREG
    5.4) for a congener and tissue ("root", "straw", "grain") of the default configuration,
    or None when the band has no entry (a congener without Yamazaki data, a SMILES record).
    down <= 1 <= up; the value times each factor is the end of the band."""
    b = INPUT_RANGE_BAND.get(congener, {}).get(tissue)
    if b is None:
        return None
    return float(10.0 ** b[0]), float(10.0 ** b[1])


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
# Dietary exposure defaults (docs/literature_db/sourced_inputs.csv exp_*), adjustable in the UI.
# Rice intake: KNHANES 2023 (KHIDI food-intake statistics, 멥쌀 = milled non-glutinous rice), the
# adult age-group means 19-29 / 30-49 / 50-64 y span 101.68-126.13 g/day; the midpoint runs.
RICE_INTAKE_G_PER_DAY_RANGE = (101.68, 126.13)
DEFAULT_RICE_INTAKE_G_PER_DAY = 0.5 * (RICE_INTAKE_G_PER_DAY_RANGE[0] + RICE_INTAKE_G_PER_DAY_RANGE[1])
# Korean adult mean body weight, KNHANES 2013-2016 (Yoon et al. 2022, JESEE 33:1013)
DEFAULT_BODY_WEIGHT_KG = 64.5
# the intake is milled white rice as bought, at the storage moisture standard (sourced
# tis_theta_grain); the model's grain pool is brown rice on the plant (`white_rice_conc`)
RICE_STORAGE_MOISTURE = 0.14


def white_rice_conc(res):
    """White (milled) rice concentration as eaten [µg/kg at the storage moisture] from a PFAS
    `simulate()` result: the harvest dry brown-rice concentration re-partitioned into white-rice
    composition at a shared free concentration (as `endosperm_white_factor`), then put on the
    14 % storage moisture. None for a result without PFAS binding constants (neutral path)."""
    p = res.get("params", {})
    if "K_PL" not in p or "conc_dw" not in res:
        return None
    kcw = p["K_cw"]["grain"] if isinstance(p.get("K_cw"), dict) else p.get("K_cw")
    c = {"K_prot_Lkg": p["K_prot"], "K_PL_Lkg": p["K_PL"]}
    factor = (_dry_binding(_grain_white_composition(), c, kcw)
              / _dry_binding(_COMP["grain_brown"], c, kcw))
    return float(res["conc_dw"]["grain"][-1] * factor * (1.0 - RICE_STORAGE_MOISTURE))


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
    `signal` is "green" below the TWI and "red" at or above it (a NaN/negative conc gives
    signal=None). The TWI itself is the only sourced threshold; an intermediate band would
    need a source of its own."""
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
    signal = "green" if pct < 100.0 else "red"
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
            ("낟알(현미, 수확 시)" if ko else "grain (brown rice at harvest)", grain_c)]
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
