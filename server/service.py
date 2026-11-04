"""Scenario -> ``model_api`` calls -> JSON-ready payloads.

This is the web backend's engine. It reproduces what the Streamlit sidebar +
``ui/common.run_model`` do (``ui/sidebar.build``: the data source, the crop drivers,
the PFAS transport variant, the neutral path) -- build the exposure drivers for the
selected data source, run the selected compound class, attach observed /
biomonitoring BAFs -- but as plain functions with a shared, de-duplicating memo, so:

* the browser handles every purely visual interaction (tab, day slider, colour
  metric) without calling back here;
* a solve is computed once and then served to every view (and every user) that
  asks for the same inputs;
* the slow extras (congener comparison, Tang check, the Yamazaki two-harvest
  observable, the Bayesian inverse) run only when their view is actually opened.

Driver dicts are passed to the model WHOLE. The Streamlit cache packs them into a
tuple (``ui/common._drivers_tuple``) that carries every key the PFAS model reads; the
only key it leaves out is ``area`` (the leaf's air-exchange surface), which matters
only for a neutral compound with plant-air exchange switched on -- there the web
run uses the crop's own surface, as ``model_api`` intends.

Two Streamlit inputs are not offered because the current model refuses them: the
"crop drivers off" placeholder logistics (``measured_forcing=False`` raises) and the
live-HYDRUS "flooded until" / "percolation" sliders (``soil_hydrus.run_paddy_hydrus``
takes the sourced Ibaraki water calendar and percolation; the keywords raise).
"""
from __future__ import annotations

import io
import inspect
import json
import math
import os
import queue
import threading
from dataclasses import asdict, dataclass, field, is_dataclass

import numpy as np

from . import ROOT
import model_api as api
import plots
import literature_params as LP
from ui import i18n

from .memo import Memo, key_of
from .schemas import (CompoundSpec, CwoPreviewRequest, DownloadRequest, FigureRequest,
                      InverseRequest, MechanismSpec, NeutralSpec, Scenario)

EX = os.path.join(ROOT, "examples")
APP_VERSION = "2.1 (web: FastAPI + React)"
REPO_URL = "https://github.com/ccy5123/pfas-rice-model"
DOCS_URL = "https://github.com/ccy5123/pfas-rice-model/tree/main/docs"
ORGANS4 = ("root", "stem", "leaf", "grain")
GENERIC_SEASON = 120.0          # ui/sidebar: the generic crops' season when no slider is shown

# ui.common.CROPS (labels) / crop_kw (the simulate() pair) -- one table, so the web and the
# Streamlit crop choice cannot drift apart
CROPS = {"site2020": "Yamazaki 2020 site-year (Tsukuba; ORYZA2000 port) — default",
         "oryza2000": "generic climatology: ORYZA2000 port",
         "oryza": "generic climatology: ORYZA Level-1, reduced",
         "growth_rice": "generic climatology: growth_rice (partition + logistic)"}

try:
    from chem_lookup import KOC_MAX_LKG as KOC_MAX, F_N_TESTED_FLOOR
except Exception:                                   # noqa: BLE001 -- optional module
    KOC_MAX, F_N_TESTED_FLOOR = 1e9, 0.1

SIM = Memo(256)     # individual model solves (shared across scenarios and users)
RUN = Memo(64)      # whole scenario runs
DRV = Memo(32)      # driver builds (a live HYDRUS-1D run takes seconds)
FIG = Memo(256)     # serialised figures
AUX = Memo(128)     # Tang check, Yamazaki observable, inverse estimates, lookups


class ApiError(Exception):
    def __init__(self, status, detail):
        super().__init__(detail)
        self.status = int(status)
        self.detail = detail


# ---------------------------------------------------------------- JSON helpers
def clean(o):
    """Recursively convert numpy / dataclass values to JSON-safe Python (NaN -> None)."""
    if isinstance(o, np.ndarray):
        return clean(o.tolist())
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, (bool, np.bool_)):
        return bool(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (float, np.floating)):
        f = float(o)
        return f if math.isfinite(f) else None
    if is_dataclass(o) and not isinstance(o, type):
        return clean(asdict(o))
    return o


def dumps(o):
    return json.dumps(clean(o), ensure_ascii=False, allow_nan=False)


# ---------------------------------------------------------------- cached model calls
def _same(a, b):
    """Scalar equality for default-stripping (arrays/dicts are never 'equal')."""
    if a is None or b is None:
        return a is b
    if isinstance(a, bool) or isinstance(b, bool):
        return type(a) is type(b) and a == b
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return float(a) == float(b)
    if isinstance(a, str) and isinstance(b, str):
        return a == b
    return False


def _canonical(fn, kw):
    """Drop kwargs that equal the function's own defaults, so two calls that mean
    the same solve share one cache entry however they were spelled."""
    sig = inspect.signature(fn).parameters
    out = {}
    for k, v in kw.items():
        p = sig.get(k)
        if p is not None and p.default is not inspect.Parameter.empty and _same(p.default, v):
            continue
        out[k] = v
    return out


def sim(fn_name, *args, **kw):
    fn = getattr(api, fn_name)
    kw = _canonical(fn, kw)
    return SIM.get(key_of(fn_name, *args, **kw), lambda: fn(*args, **kw))


# ---------------------------------------------------------------- scenario pieces
def crop_kw(crop):
    """The crop choice -> simulate()'s (scenario, biomass) pair (`ui.common.crop_kw`)."""
    if crop == "site2020":
        return dict(scenario="yamazaki2020", biomass="oryza2000")
    return dict(scenario="generic", biomass=crop)


def site_season():
    return AUX.get(key_of("site_season"), api.site_season)


def season_of(crop, e=None):
    """The season the Streamlit sidebar runs: the site-year sets its own (transplanting to
    harvest); a generic crop takes the parametric season slider, else 120 d."""
    if crop == "site2020":
        return site_season()
    if e is not None and e.mode == "parametric":
        return float(e.season)
    return GENERIC_SEASON


def karickhoff(log_kow):
    return float(LP.koc_neutral(float(log_kow)))


def neutral_koc(n: NeutralSpec):
    return float(n.Koc) if n.Koc is not None else karickhoff(n.log_kow)


def neutral_kwargs(n: NeutralSpec):
    """``simulate_neutral`` kwargs -- the dict ``ui/sidebar._neutral_panel`` returns
    (without ``Koc``, which only the soil-side exposure shapes use). ``grain_feed`` is
    passed only when it departs from the model default."""
    kw = dict(log_kow=float(n.log_kow), name=(n.name or "neutral"),
              half_life=(float(n.half_life) or None), tscf_model=n.tscf_model,
              lipid_source=n.lipid_source, phloem=bool(n.phloem), air=bool(n.air),
              g_apo=float(n.g_apo))
    if n.air:
        kw.update(MW=float(n.MW), K_AW=float(n.K_AW), air_kw=dict(C_air=float(n.C_air)))
    if n.weak_electrolyte:
        kw.update(pKa=float(n.pKa), pH=float(n.pH), is_acid=bool(n.is_acid))
    if n.grain_feed != "xylem":
        kw["grain_feed"] = n.grain_feed
    return kw


def resolve_mechanism(m: MechanismSpec):
    """(sim kwargs, messages): the PFAS ``simulate`` kwargs ``ui/common.run_model`` builds,
    after the sidebar's own fall-backs -- f_xy 'fit' exists for the free anion on the
    redistributed shoot only, so with lipid loading or on the 4-pool shoot the run uses the
    monotone 'recommended' f_xy, with a message (``ui/sidebar.build``)."""
    fxy, msgs = m.fxy_source, []
    if m.lipid_loading and fxy == "fit":
        msgs.append({"level": "warning", "text": (
            "f_xy 'fit' is a free-anion fit and applies only with lipid loading off (⚙️ "
            "Mechanism); this run uses lipid loading's own f_xy.")})
        fxy = "recommended"
    if m.shoot == "4pool" and fxy == "fit":
        msgs.append({"level": "warning", "text": (
            "f_xy 'fit' was fitted on the redistributed shoot only (no 4-pool fit); this "
            "4-pool run uses the monotone 'recommended' f_xy.")})
        fxy = "recommended"
    kw = dict(f_xy_source=fxy, **crop_kw(m.crop), lipid_loading=bool(m.lipid_loading),
              shoot=m.shoot)
    return kw, msgs


def sim_kw(m: MechanismSpec):
    return resolve_mechanism(m)[0]


def variant_info(kw):
    """The fitted variant a set of sim kwargs runs and its withdrawn flag (PREREG 9.4)."""
    v = api.variant_key(kw.get("shoot"), kw.get("f_xy_source"), kw.get("lipid_loading"))
    return v, api.WITHDRAWN_VARIANTS.get(v)


# ---- soil-inventory seeds (ui/sidebar "Soil inventory -> pore water") ----
def inventory_defaults(comp: CompoundSpec):
    """The values the Streamlit inventory widgets are seeded with: the Yamazaki field's
    measured pre-season soil (SI Table S16, water + methanol extracts) for a curated congener
    (0 when it has none, and for a SMILES or neutral compound), and the compound's linear Kd
    on the site soil (Koc x f_oc; a SMILES compound has no congener and gets 1.0, as there)."""
    neutral = comp.kind == "neutral"
    congener = comp.congener if comp.kind == "congener" else None
    pre = None
    if congener and not neutral:
        import soil_hydrus as sh
        try:
            pre = sh.post_season_soil(api.YAMAZAKI_S16_NAME.get(congener, congener))[0] / 1000.0
        except KeyError:
            pre = None
    if neutral:
        kd = float(neutral_koc(comp.neutral) or 0.0) * api.SOIL_FOC
    else:
        kd = api.site_kd(congener) if congener else 1.0
    return dict(C_total=float(pre) if pre is not None else 0.0,
                C_total_sourced=pre is not None, K_F=max(float(kd), 1e-4))


# ---- editable tables (pure versions of ui/common's DataFrame helpers) ----
def default_growth_table(season, crop, n_rows=7):
    """`ui/common._default_growth_df`: the crop drivers' FRESH organ mass [g/hill]."""
    kw = crop_kw(crop)
    if kw["scenario"] == "yamazaki2020":
        season = site_season()
    t = np.linspace(0.0, float(season), n_rows)
    _, M, _, _ = api.crop_on_grid(t, season=season, **kw)
    out = {"day": np.round(t, 0).tolist()}
    for i, o in enumerate(ORGANS4):
        out[o] = np.round(M[:, i] * 1e3, 2).tolist()
    return out


def default_cwo_table(season, level, n_rows=4):
    return {"day": np.round(np.linspace(0.0, float(season), n_rows), 0).tolist(),
            "Cwo": [float(level)] * n_rows}


def _num(v):
    try:
        f = float(v)
    except (TypeError, ValueError):
        return float("nan")
    return f


def clean_table(tbl, value_cols):
    """{column: [values]} -> {column: array}, dropping incomplete rows; day required."""
    if not tbl or "day" not in tbl:
        raise ValueError("the table needs a 'day' column")
    day = np.array([_num(v) for v in tbl["day"]], float)
    out = {"day": day}
    mask = np.isfinite(day)
    for c in value_cols:
        if c in tbl:
            col = list(tbl[c])[:len(day)] + [None] * max(0, len(day) - len(tbl[c]))
            arr = np.array([_num(v) for v in col], float)
            out[c] = arr
            mask &= np.isfinite(arr)
    out = {k: v[mask] for k, v in out.items()}
    if len(out["day"]) < 2:
        raise ValueError("need at least 2 complete rows (day + values)")
    return out


def hydrus_soil_congener(smiles):
    """SMILES -> curated congener for the HYDRUS soil Kd (read-across match, else the
    nearest by perfluoro-C in the same head-group family). (name, how) or (None, None)."""
    try:
        from pfas_structure import descriptors
        d = descriptors(smiles)
    except Exception:                                   # noqa: BLE001
        return None, None
    if d.matched_name:
        return d.matched_name, "match"
    fam = "PFSA" if d.head_group == "sulfonate" else ("ether" if d.n_ether_O else "PFCA")
    cands = [c for c in api._CONG.values() if c["group"] == fam] \
        or [c for c in api._CONG.values() if c["group"] == "PFCA"]
    npfc = lambda c: c["n_C"] - 1 if c["group"] == "PFCA" else c["n_C"]   # noqa: E731
    best = min(cands, key=lambda c: abs(npfc(c) - d.n_perfluoroC))
    return best["name"], "nearest"


def _hydrus_drivers(congener, season, f_oc, crop, Kd=None):
    def go():
        drv, _ = api.hydrus_drivers(congener, season=season, f_oc=f_oc, Kd=Kd, **crop_kw(crop))
        return drv
    return DRV.get(key_of("hydrus", congener, season, f_oc, crop, Kd), go)


# ---------------------------------------------------------------- drivers per data source
@dataclass
class Built:
    drivers: dict | None = None
    drivers_pchip: dict | None = None
    soil_obj: object = None
    measured_bio: dict | None = None
    messages: list = field(default_factory=list)
    # the effective "built-in forcing" settings (used when drivers is None)
    Cwo: float = 1.0
    season: float | None = None
    cwo_profile: str = "constant"
    extra: dict = field(default_factory=dict)


def build(scn: Scenario) -> Built:
    """Exposure drivers for the selected data source -- ``ui/sidebar.build`` minus widgets."""
    e, comp, lang = scn.exposure, scn.compound, scn.lang
    crop = scn.mechanism.crop
    ckw = crop_kw(crop)
    season = season_of(crop, e)
    b = Built(season=None if crop == "site2020" else season)

    def msg(level, text):
        b.messages.append({"level": level, "text": text})

    # Only the parametric and custom-table sources expose the level and the shape; the
    # others run the Streamlit sidebar defaults (Cwᵒ 1, flat).
    if e.mode in ("parametric", "tables"):
        b.Cwo = float(e.Cwo)
    if e.mode == "parametric":
        b.cwo_profile = e.cwo_profile

    if e.mode == "tables":
        t = e.tables
        try:
            growth = clean_table(t.growth or default_growth_table(season, crop), ORGANS4)
            cwo = clean_table(t.cwo or default_cwo_table(season, e.Cwo), ["Cwo"])
            b.drivers = api.drivers_from_tables(growth, cwo, growth_units=t.growth_units,
                                                Cwo_const=float(e.Cwo), **ckw)
            # the form between the table's samples is not measured (PREREG 7.3): the
            # monotone-cubic run rides along so the result can show the spread
            b.drivers_pchip = api.drivers_from_tables(growth, cwo, growth_units=t.growth_units,
                                                      Cwo_const=float(e.Cwo), interp="pchip",
                                                      **ckw)
        except Exception as ex:                          # noqa: BLE001
            msg("warning", i18n.t("ct.read_error", lang, e=ex))

    elif e.mode == "csv":
        try:
            if e.csv_text:
                import pandas as pd
                df = pd.read_csv(io.StringIO(e.csv_text))
                cols = {c.lower().strip(): c for c in df.columns}
                tt = df[cols["t"]].to_numpy(float)
                Cwo = df[cols["cwo"]].to_numpy(float)
                Q = df[cols["qtp"]].to_numpy(float) if "qtp" in cols else None
                M = (df[[cols["m_root"], cols["m_stem"], cols["m_leaf"], cols["m_grain"]]]
                     .to_numpy(float)
                     if all(k in cols for k in ("m_root", "m_stem", "m_leaf", "m_grain")) else None)
                b.drivers = api.drivers_from_arrays(tt, Cwo, Qtp=Q, M=M, **ckw)
                msg("success", f"Loaded {len(tt)} rows from the uploaded CSV.")
            elif e.csv_text is None:
                b.drivers = DRV.get(key_of("example_csv", crop), lambda: api.load_driver_csv(
                    os.path.join(EX, "hydrus_drivers_example.csv"), **ckw))
                msg("info", "Using examples/hydrus_drivers_example.csv (synthetic HYDRUS-style run).")
            else:
                msg("info", "No driver CSV — running the default parametric scenario.")
        except Exception as ex:                          # noqa: BLE001
            msg("error", f"Could not read drivers: {ex}")

    elif e.mode == "hydrus":
        f_oc = None if e.f_oc is None else float(e.f_oc)
        if not api.hydrus_available():
            msg("warning", "HYDRUS-1D engine not built in this environment yet — build it once "
                           "(compiles the FORTRAN solver with gfortran; ~1 min). Until built, the "
                           "tool falls back to the parametric model.")
        elif comp.kind == "neutral":
            import soil_hydrus as sh
            kd = neutral_koc(comp.neutral) * (sh.FOC if f_oc is None else f_oc)
            try:
                b.drivers = _hydrus_drivers(None, season, f_oc, crop, Kd=kd)
                msg("success", f"HYDRUS-1D run complete — Kd = {kd:.3g} L/kg "
                               f"(Koc {neutral_koc(comp.neutral):.1f} × f_oc "
                               f"{sh.FOC if f_oc is None else f_oc:g}; mean-normalised).")
            except Exception as ex:                      # noqa: BLE001
                msg("error", f"HYDRUS run failed: {ex}")
        else:
            soil_cong, how = (comp.congener, None)
            if comp.kind == "smiles":
                soil_cong, how = hydrus_soil_congener(comp.smiles or "")
            if soil_cong is None:
                msg("error", "Could not parse the SMILES for the HYDRUS soil run — check the "
                             "structure or switch the compound to a curated congener.")
            else:
                if comp.kind == "smiles":
                    msg("info", f"SMILES compound → soil Kd uses **{soil_cong}** "
                                f"({'read-across match' if how == 'match' else 'nearest curated congener by chain length'}); "
                                "the plant uptake still uses your structure.")
                try:
                    b.drivers = _hydrus_drivers(soil_cong, season, f_oc, crop)
                    msg("success", f"HYDRUS-1D run complete — Cwᵒ(t) for {soil_cong} "
                                   "(Kd-retarded; mean-normalised).")
                except Exception as ex:                  # noqa: BLE001
                    msg("error", f"HYDRUS run failed: {ex}")

    elif e.mode == "inventory":
        neutral = comp.kind == "neutral"
        congener = comp.congener if comp.kind == "congener" else None
        seed = inventory_defaults(comp)
        C_total = seed["C_total"] if e.C_total is None else float(e.C_total)
        K_F = seed["K_F"] if e.K_F is None else float(e.K_F)
        tt = np.linspace(0.0, season, 241)
        Cwo, soil = api.pore_water_from_inventory(
            tt, C_total, K_F=K_F, n=float(e.n_F), congener=None if neutral else congener,
            Koc=neutral_koc(comp.neutral) if neutral else None)
        b.drivers = api.drivers_from_arrays(tt, Cwo, season=season, **ckw)
        b.soil_obj = soil
        b.extra["inventory"] = dict(C_total=C_total, K_F=K_F, n_F=float(e.n_F))

    elif e.mode == "biomon":
        try:
            if e.biomon_csv_text:
                b.measured_bio = api.load_biomonitoring_csv(io.StringIO(e.biomon_csv_text))
                msg("success", f"Loaded {len(b.measured_bio['conc'])} tissues; "
                               f"Cwᵒ={b.measured_bio['Cwo']}.")
            elif e.biomon is not None:
                b.measured_bio = dict(conc=dict(e.biomon.conc), Cwo=e.biomon.Cwo)
            else:
                b.measured_bio = api.load_biomonitoring_csv(
                    os.path.join(EX, "biomonitoring_example.csv"))
                msg("info", "Using examples/biomonitoring_example.csv (Yamazaki PFOA).")
        except Exception as ex:                          # noqa: BLE001
            msg("error", f"Could not read CSV: {ex}")
    return b


# ---------------------------------------------------------------- the run
@dataclass
class Run:
    scn: Scenario
    built: Built
    res: dict
    congener: str
    obs: dict
    bio_baf: dict | None
    sim_kw: dict
    messages: list
    desc: object = None
    provisional: bool = False


def _forcing_kw(b: Built):
    """The built-in-forcing kwargs of a run without drivers (`run_model`'s scalar path)."""
    return dict(Cwo=b.Cwo, season=b.season, measured_forcing=True, cwo_profile=b.cwo_profile)


def _run(scn: Scenario) -> Run:
    b = build(scn)
    comp, m = scn.compound, scn.mechanism
    skw, mech_msgs = resolve_mechanism(m)
    messages = list(b.messages)
    desc, provisional = None, False
    if comp.kind == "neutral":
        # the PFAS levers (f_xy source, shoot, lipid loading) have no meaning at z = 0
        nk = neutral_kwargs(comp.neutral)
        log_kow = nk.pop("log_kow")
        if b.drivers is not None:
            res = sim("simulate_neutral", log_kow, drivers=b.drivers, **crop_kw(m.crop), **nk)
        else:
            res = sim("simulate_neutral", log_kow, **_forcing_kw(b), **crop_kw(m.crop),
                      Koc=neutral_koc(comp.neutral), **nk)
    else:
        messages += mech_msgs
        if comp.kind == "smiles":
            smiles = (comp.smiles or "").strip()
            if not smiles:
                raise ApiError(400, "Enter a SMILES structure.")
            if not api.rdkit_available():
                raise ApiError(400, "RDKit not installed — cannot parameterise a SMILES structure. "
                                    "`pip install rdkit`, or switch to a curated congener.")
            try:
                if b.drivers is not None:
                    res = sim("simulate_from_smiles", smiles, drivers=b.drivers, **skw)
                else:
                    res = sim("simulate_from_smiles", smiles, **_forcing_kw(b), **skw)
            except ApiError:
                raise
            except Exception as ex:                      # noqa: BLE001
                raise ApiError(400, "Could not build a compound from that SMILES — check the "
                                    f"structure. ({ex})")
            desc = res.get("descriptors")
            provisional = bool(res.get("provisional", False))
        else:
            if comp.congener not in api.CONGENERS:
                raise ApiError(400, f"unknown congener {comp.congener!r}; known: {api.CONGENERS}")
            if b.drivers is not None:
                res = sim("simulate", comp.congener, drivers=b.drivers, **skw)
                if b.drivers_pchip is not None:          # a measured table: both forms
                    alt = sim("simulate", comp.congener, drivers=b.drivers_pchip,
                              **skw)["baf_final_dw"]
                    res = dict(res, interp_spread={
                        k: float(np.log10(alt[k] / res["baf_final_dw"][k]))
                        for k in ("root", "straw", "grain")})
            else:
                res = sim("simulate", comp.congener, **_forcing_kw(b), **skw)
    congener = res["congener"]
    obs = {} if comp.kind == "neutral" else api.observed_baf(congener)
    mb = b.measured_bio
    bio_baf = (api.baf_from_measurement(mb["conc"], mb["Cwo"])
               if mb and mb.get("Cwo") else None)
    return Run(scn=scn, built=b, res=res, congener=congener, obs=obs, bio_baf=bio_baf,
               sim_kw=skw, messages=messages, desc=desc, provisional=provisional)


def run(scn: Scenario) -> Run:
    return RUN.get(key_of("run", scn), lambda: _run(scn))


def white_conc(res):
    """The intake basis: milled white rice as eaten (`api.white_rice_conc`), or the brown
    rice on the plant when the run has no white-rice value (the neutral path)."""
    w = api.white_rice_conc(res)
    return float(res["conc"]["grain"][-1]) if w is None else float(w)


def payload(r: Run) -> dict:
    """Everything the front end needs to draw the headline, the plant map (scrubbed
    client-side over the whole season) and the series -- one compact JSON object."""
    res, b = r.res, r.built
    root_c = float(res["conc"]["root"][-1])
    straw_c = float(res["straw"][-1])
    grain_c = float(res["conc"]["grain"][-1])
    white_c = white_conc(res)
    M = np.asarray(res["M"], float)
    schematic = {}
    for metric in ("conc", "baf"):
        ms = api.metric_series(res, metric)
        schematic[metric] = {"cmin": ms["cmin"], "cmax": ms["cmax"], "label": ms["label"]}
    p = dict(res.get("params", {}))
    neutral = r.scn.compound.kind == "neutral"
    variant = None
    if not neutral:
        v, why = variant_info(r.sim_kw)
        variant = dict(key=v, withdrawn=why)
    out = dict(
        congener=r.congener, kind=r.scn.compound.kind, mode=r.scn.exposure.mode,
        success=bool(res.get("success", True)),
        t=res["t"], conc=res["conc"], baf=res["baf"], straw=res["straw"],
        Cwo=res["Cwo"], Qtp=res["Qtp"], M={k: M[:, i] for i, k in enumerate(api.TISSUES)},
        baf_final=res["baf_final"], straw_baf=res["straw_baf"], cwo_ref=res.get("cwo_ref"),
        baf_final_dw=res.get("baf_final_dw"),
        conc_dw_final={k: float(np.asarray(v, float)[-1])
                       for k, v in (res.get("conc_dw") or {}).items()},
        B_k=res.get("B_k"), params=p, variant=variant,
        finals=dict(root=root_c, straw=straw_c, grain=grain_c, white=white_c),
        bands={k: api.predictive_band(v) for k, v in
               (("root", root_c), ("straw", straw_c), ("grain", grain_c), ("white", white_c))},
        input_range={k: api.input_range_factors(r.congener, k) for k in ("root", "straw", "grain")},
        uncertainty_factor=api.uncertainty_factor(),
        intake=api.intake_fraction(white_c, congener=r.congener),
        interp_spread=res.get("interp_spread"),
        schematic=schematic, obs=r.obs, bio_baf=r.bio_baf, measured_bio=b.measured_bio,
        descriptors=(asdict(r.desc) if r.desc is not None else None),
        provisional=r.provisional, messages=r.messages,
        has_isotherm=b.soil_obj is not None, inventory=b.extra.get("inventory"),
        effective=dict(Cwo=b.Cwo, season=float(res["t"][-1]), cwo_profile=b.cwo_profile,
                       drivers=b.drivers is not None, sim_kw=(None if neutral else r.sim_kw)),
    )
    if neutral:
        out.update(TSCF=res.get("TSCF"), K_PW=res.get("K_PW"), rcf_briggs=res.get("rcf_briggs"),
                   half_life=res.get("half_life"), warning=res.get("warning"),
                   air_summary=res.get("air_summary"), tf_final=res.get("tf_final"),
                   grain_predicted=res.get("grain_predicted"),
                   koc=neutral_koc(r.scn.compound.neutral))
    return clean(out)


# ---------------------------------------------------------------- figures
def _fig_out(fig, data=None):
    """'{"figure": <plotly JSON>, "data": {...}}' -- figure JSON straight from plotly
    (its encoder already maps NaN to null)."""
    return '{"figure":' + fig.to_json() + ',"data":' + dumps(data or {}) + '}'


def tang(congener, dose, lipid=False):
    """`ui/expert` Tang tab: the default configuration (and optionally the lipid loading)."""
    kw = dict(dose=dose)
    if lipid:
        kw["lipid_loading"] = True
    return AUX.get(key_of("tang", congener, kw),
                   lambda: api.tang_tf_validation(congener, **kw))


def yamazaki_observable(congener, mech):
    """The model's version of the Yamazaki 2023 observation (`ui/common._yamazaki_observable`)."""
    return AUX.get(key_of("yamazaki_obs", congener, mech), lambda: {
        "congener": congener, "baf_final_dw": api.yamazaki_harvest_baf(congener, **mech)})


def observable_mech(view, skw):
    """The transport settings the observable follows: Expert passes the sidebar's f_xy source
    and mechanism switches; Simple passes the f_xy source only (`ui/simple`)."""
    if view == "simple":
        return dict(f_xy_source=skw["f_xy_source"])
    return {k: v for k, v in skw.items() if k not in ("biomass", "scenario")}


SIMPLE_COMPARE = ("PFBA", "PFHxA", "PFOA", "PFDA", "PFDoDA", "PFOS")


def _intake_info(r, rice, bw):
    white_c = white_conc(r.res)
    info = api.intake_fraction(white_c, congener=r.congener, rice_intake_g_day=rice,
                               body_weight_kg=bw)
    band = api.predictive_band(white_c)
    lo = api.intake_fraction(band["lo"], congener=r.congener, rice_intake_g_day=rice,
                             body_weight_kg=bw)
    hi = api.intake_fraction(band["hi"], congener=r.congener, rice_intake_g_day=rice,
                             body_weight_kg=bw)
    return info, lo, hi, white_c


def _figure(req: FigureRequest) -> str:
    name, lang, o = req.name, req.lang, req.opts
    r = run(req.scenario)
    res, obs = r.res, r.obs
    comp = req.scenario.compound
    if name == "tissue":
        return _fig_out(plots.fig_tissue(res))
    if name == "burden":
        return _fig_out(plots.fig_burden(res))
    if name == "drivers":
        return _fig_out(plots.fig_drivers(res))
    if name == "soil_profile":
        return _fig_out(plots.fig_soil_profile(res, profile=None))
    if name == "isotherm":
        if r.built.soil_obj is None:
            raise ApiError(400, "the isotherm needs the soil-inventory data source")
        return _fig_out(plots.fig_isotherm(r.built.soil_obj, Cwo_now=float(res["Cwo"][-1])))
    if name == "buildup_plain":
        return _fig_out(plots.fig_buildup_plain(res, lang=lang))
    if name == "where_plain":
        return _fig_out(plots.fig_where_plain(res, lang=lang, band=True))
    if name == "baf":                                     # this run's dry-basis BAF vs Yamazaki
        if comp.kind == "neutral":
            raise ApiError(400, "no Yamazaki observation for a neutral compound")
        return _fig_out(plots.fig_baf(res, obs, lang=lang))
    if name == "baf_observable":                          # the Yamazaki two-harvest observable
        if comp.kind != "congener" or not obs:
            raise ApiError(400, "the Yamazaki observable needs a curated congener with data")
        mech = observable_mech(str(o.get("view", "expert")), r.sim_kw)
        return _fig_out(plots.fig_baf(yamazaki_observable(r.congener, mech), obs, lang=lang),
                        {"mech": mech})
    if name == "biomon_compare":
        if not r.bio_baf:
            raise ApiError(400, "no biomonitoring BAF for this scenario")
        model_b = {"root": res["baf_final"]["root"], "straw": res["straw_baf"],
                   "grain": res["baf_final"]["grain"]}
        return _fig_out(plots.fig_biomon_compare(r.bio_baf, model_b))
    if name == "chain":
        key = str(o.get("key", "K_PL"))
        if key not in ("K_PL", "K_prot", "K_cw_root", "f_xy_recommended", "B_root", "B_grain"):
            raise ApiError(400, f"unknown chain parameter {key!r}")
        desc = r.desc
        cc = r.congener if r.congener in api.CONGENERS else (
            desc.matched_name if desc is not None and desc.matched_name in api.CONGENERS else "PFOA")
        return _fig_out(plots.fig_chain(api.chain_table(), cc, key), {"ring": cc})
    if name == "compare":                                 # Expert: overlay congeners
        congs = [c for c in (o.get("congeners") or []) if c in api.CONGENERS]
        tissue = str(o.get("tissue", "straw"))
        b = r.built
        results = {}
        for nm in congs:
            if b.drivers is not None:
                results[nm] = sim("simulate", nm, drivers=b.drivers, **r.sim_kw)
            else:      # ui/expert: the sidebar level and season, the flat exposure
                results[nm] = sim("simulate", nm, Cwo=b.Cwo, season=b.season,
                                  measured_forcing=True, **r.sim_kw)
        return _fig_out(plots.fig_compare(results, tissue))
    if name == "congener_compare":                        # Simple: six representative
        b = r.built
        m = req.scenario.mechanism
        kw = dict(f_xy_source=r.sim_kw["f_xy_source"], **crop_kw(m.crop))
        rep = [c for c in SIMPLE_COMPARE if c in api.CONGENERS]
        results = {c: sim("simulate", c, Cwo=b.Cwo, season=b.season, measured_forcing=True, **kw)
                   for c in rep}
        return _fig_out(plots.fig_congener_compare(results, lang=lang, order=rep))
    if name == "tang_tf":
        if comp.kind != "congener" or r.congener not in api.TANG_CONGENERS:
            raise ApiError(400, f"Tang 2026 covers {', '.join(api.TANG_CONGENERS)} only")
        dose = "low" if o.get("dose") == "low" else "mean"
        val = tang(r.congener, dose)
        vall = tang(r.congener, dose, lipid=True) if o.get("lipid") else None
        err = [np.log10(val["model_tf"][g] / val["tang_tf"][g]) for g in val["organs"]
               if val["tang_tf"].get(g)]
        return _fig_out(plots.fig_tang_tf(val, vall, extra_label="model · lipid loading"),
                        {"f_xy": val["f_xy"], "rmse": float(np.sqrt(np.mean(np.square(err))))})
    if name == "intake_gauge":
        rice = float(o.get("rice_g", api.DEFAULT_RICE_INTAKE_G_PER_DAY))
        bw = float(o.get("bw", api.DEFAULT_BODY_WEIGHT_KG))
        info, lo, hi, white_c = _intake_info(r, rice, bw)
        return _fig_out(plots.fig_intake_gauge(info, lang=lang),
                        {"info": info, "lo_percent": lo["percent"], "hi_percent": hi["percent"],
                         "white": white_c})
    raise ApiError(404, f"unknown figure {name!r}")


def figure(req: FigureRequest) -> str:
    return FIG.get(key_of("fig", req), lambda: _figure(req))


def cwo_preview(req: CwoPreviewRequest) -> str:
    """The flooded-shape preview under the shape radio (`plots.fig_cwo_profile`, the stored
    HYDRUS shapes). A neutral compound's shape comes from its own Koc; a SMILES compound
    takes the generic PFCA descriptor, as in the Streamlit preview."""
    comp = req.compound
    congener = comp.congener if comp.kind == "congener" else None
    koc = neutral_koc(comp.neutral) if comp.kind == "neutral" else None
    season = site_season() if req.crop == "site2020" else float(req.season)

    def go():
        return _fig_out(plots.fig_cwo_profile(congener, level=req.level, profile="flooded",
                                              season=season, Koc=koc))
    return FIG.get(key_of("cwo_preview", congener, koc, req.level, season), go)


# ---------------------------------------------------------------- Bayesian inverse
def meas_sigma_log10(rsd_percent):
    """`ui/common.meas_sigma_log10`: relative SD [%] -> log10 SD (lognormal); 0 -> 0."""
    cv = max(float(rsd_percent), 0.0) / 100.0
    return float(np.sqrt(np.log1p(cv * cv)) / np.log(10.0))


def _inverse_args(req: InverseRequest):
    """The `estimate_exposure_bayesian` call `ui/common._render_inverse_estimator` makes."""
    scn, lang = req.scenario, req.lang
    comp = scn.compound
    meas = {k: float(v) for k, v in req.measured.items()
            if k in ("root", "straw", "grain") and v is not None and float(v) > 0}
    if not meas:
        raise ApiError(400, i18n.t("inv.enter_first", lang))
    skw = sim_kw(scn.mechanism)
    base = dict(meas_sigma_log10=meas_sigma_log10(req.rsd_percent),
                f_xy_source=skw["f_xy_source"], **crop_kw(scn.mechanism.crop))
    if comp.kind == "neutral":
        nk = neutral_kwargs(comp.neutral)
        label = nk["name"]
        mech = dict(log_kow=float(nk["log_kow"]),
                    neutral_kw={k: v for k, v in nk.items() if k not in ("log_kow", "name")})
        return label, meas, {**base, **mech}
    if comp.kind != "congener" or comp.congener not in api.CONGENERS:
        raise ApiError(400, i18n.t("inv.not_curated", lang))
    mech = ({} if req.view == "simple" else
            {k: v for k, v in skw.items() if k not in ("f_xy_source", "biomass", "scenario")})
    return comp.congener, meas, {**base, **mech}


def inverse_stream(req: InverseRequest):
    """NDJSON lines: {"type":"progress",done,total} ... then {"type":"result",...} (or
    {"type":"error"}). The estimate (~8 ODE solves) is memoised, so a repeat streams
    the result at once."""
    try:
        label, meas, kw = _inverse_args(req)
    except ApiError as ex:
        yield dumps({"type": "error", "detail": ex.detail}) + "\n"
        return
    key = key_of("inverse", label, meas, kw)
    q = queue.Queue()
    box = {}

    def progress(done, total):
        q.put(("progress", int(done), int(total)))

    def work():
        try:
            box["est"] = AUX.get(key, lambda: api.estimate_exposure_bayesian(
                label, meas, progress=progress, **kw))
        except Exception as ex:                          # noqa: BLE001
            box["err"] = str(ex)
        finally:
            q.put(("done",))

    threading.Thread(target=work, daemon=True).start()
    while True:
        item = q.get()
        if item[0] == "progress":
            yield dumps({"type": "progress", "done": item[1], "total": item[2]}) + "\n"
        else:
            break
    if "err" in box:
        yield dumps({"type": "error", "detail": i18n.t("inv.error", req.lang, e=box["err"])}) + "\n"
        return
    est = box["est"]
    fig = plots.fig_exposure_posterior(est, lang=req.lang)
    yield ('{"type":"result","estimate":' + dumps(est) + ',"figure":' + fig.to_json() + '}\n')


# ---------------------------------------------------------------- downloads
def download(kind: str, req: DownloadRequest):
    """-> (bytes, media_type, filename_ascii, filename_utf8)."""
    r = run(req.scenario)
    res, cong, lang = r.res, r.congener, req.lang
    ko = lang == "ko"
    if kind == "summary_csv":
        return (api.summary_csv(res, r.obs, r.bio_baf).encode("utf-8"), "text/csv",
                f"{cong}_summary.csv", None)
    if kind == "timeseries_csv":
        return (api.timeseries_csv(res).encode("utf-8"), "text/csv",
                f"{cong}_timeseries.csv", None)
    if kind == "summary_md":                 # ui/simple: the intake on white rice as eaten
        intake = api.intake_fraction(white_conc(res), congener=cong)
        md = api.summary_report_md(res, congener=cong, obs=r.obs, intake=intake,
                                   scenario=req.scenario_label, lang=lang)
        return (md.encode("utf-8"), "text/markdown", f"{cong}_summary.md",
                f"{cong}_요약.md" if ko else None)
    if kind == "summary_html":
        import report_html
        html = report_html.summary_html(res, congener=cong, disclaimer=i18n.t("disclaimer", lang),
                                        scenario=req.scenario_label, lang=lang)
        return (html, "text/html", f"{cong}_summary.html", f"{cong}_요약.html" if ko else None)
    if kind == "map_html":
        fig = plots.fig_schematic_from_res(res, "conc", -1, lang=lang)
        return (fig.to_html(full_html=True, include_plotlyjs="cdn").encode("utf-8"),
                "text/html", f"{cong}_map.html", None)
    raise ApiError(404, f"unknown download {kind!r}")


# ---------------------------------------------------------------- small endpoints
def tables_defaults(crop, season, Cwo):
    s = site_season() if crop == "site2020" else float(season)
    return dict(growth=default_growth_table(s, crop), cwo=default_cwo_table(s, Cwo))


def mol_svg(smiles, w=290, h=170):
    """2-D structure SVG (RDKit) -> {"svg": str} or {"error": reason}."""
    def go():
        try:
            from rdkit import Chem
            from rdkit.Chem.Draw import rdMolDraw2D
        except Exception as ex:                          # noqa: BLE001
            return {"error": f"RDKit import failed: {ex}"}
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return {"error": "RDKit could not parse this SMILES."}
            d = rdMolDraw2D.MolDraw2DSVG(w, h)
            d.drawOptions().padding = 0.12
            d.DrawMolecule(mol)
            d.FinishDrawing()
            return {"svg": d.GetDrawingText()}
        except Exception as ex:                          # noqa: BLE001
            return {"error": f"draw error ({type(ex).__name__}): {ex}"}
    return AUX.get(key_of("molsvg", smiles, w, h), go)


def ctx_lookup(q):
    """EPA CompTox lookup -> the dict the Streamlit panel renders. Only successful
    lookups are memoised (a missing key / network error should be retried)."""
    import chem_lookup as cx
    key = key_of("ctx", q)
    hit = AUX.peek(key)
    if hit is not None:
        return hit
    r = cx.lookup(q)
    out = dict(ok=r.ok, note=r.note, dtxsid=r.dtxsid, name=r.preferred_name, casrn=r.casrn,
               smiles=r.smiles, kwargs=clean(r.neutral_kwargs()) if r.ok else {},
               badges={k: v.badge() for k, v in r.props.items()},
               sources={k: v.source for k, v in r.props.items()},
               outside_ad={k: bool(v.outside_ad) for k, v in r.props.items()})
    out["pka_measured"] = any(out["sources"].get(k) == "experimental"
                              for k in ("pka_acidic", "pka_basic", "pka_unlabelled"))
    if r.ok:
        AUX.get(key, lambda: out)
    return out


def ctx_available():
    try:
        import chem_lookup as cx
        return bool(cx.api_key())
    except Exception:                                    # noqa: BLE001
        return False


EXAMPLE_SMILES = [
    ("PFOA  (known)", "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"),
    ("PFOS  (known)", "OS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"),
    ("GenX / HFPO-DA  (known ether)", "OC(=O)C(F)(OC(F)(F)C(F)(F)C(F)(F)F)C(F)(F)F"),
    ("PFTrDA — NOVEL C13 PFCA", "OC(=O)" + "C(F)(F)" * 11 + "C(F)(F)F"),
    ("ADONA-like — NOVEL ether-PFCA", "OC(=O)C(F)(F)OC(F)(F)C(F)(F)OC(F)(F)C(F)(F)F"),
]


def meta():
    import soil_hydrus as sh
    import neutral_dpu as ND
    return clean(dict(
        app_version=APP_VERSION, repo_url=REPO_URL, docs_url=DOCS_URL,
        congeners=api.CONGENERS,
        congener_info={n: {"n_C": c["n_C"], "group": c["group"]} for n, c in api._CONG.items()},
        labels=i18n.CONGENER_LABELS,
        presets={"en": [[k, v] for k, v in i18n.PRESETS["en"].items()],
                 "ko": [[k, v[0], v[1]] for k, v in i18n.PRESETS["ko"].items()]},
        scenarios_ko=[[k, *v] for k, v in i18n.SCENARIOS_KO.items()],
        strings=i18n.STRINGS,
        crops=[[k, v] for k, v in CROPS.items()], default_crop="site2020",
        site_season=site_season(), generic_season=GENERIC_SEASON,
        shoots=list(api.SHOOTS), default_shoot=api.DEFAULT_SHOOT,
        default_f_xy_source=api.DEFAULT_F_XY_SOURCE,
        default_lipid_loading=api.DEFAULT_LIPID_LOADING,
        default_variant=api.variant_key(), withdrawn_variants=api.WITHDRAWN_VARIANTS,
        tang_congeners=list(api.TANG_CONGENERS),
        example_smiles=EXAMPLE_SMILES,
        intake=dict(rice_g=api.DEFAULT_RICE_INTAKE_G_PER_DAY,
                    rice_range=list(api.RICE_INTAKE_G_PER_DAY_RANGE),
                    bw=api.DEFAULT_BODY_WEIGHT_KG, twi=api.EFSA_TWI_NG_PER_KG_BW_WEEK,
                    congeners=list(api.EFSA_TWI_CONGENERS)),
        uncertainty_factor=api.uncertainty_factor(),
        apriori_log10_rmse=api.APRIORI_LOG10_RMSE,
        neutral_scores=api.NEUTRAL_APRIORI_OOS_SCORES,
        neutral_equilibrium_scores=api.NEUTRAL_EQUILIBRIUM_SCORES,
        neutral_lipid_fw=dict(measured=ND.RICE_LIPID_FW["root"],
                              briggs_anchor=ND.BRIGGS_ANCHORED_LIPID_FW["root"]),
        karickhoff=dict(slope=LP.KOC_NEUTRAL_SLOPE, intercept=LP.KOC_NEUTRAL_INTERCEPT),
        koc_max=KOC_MAX, f_n_tested_floor=F_N_TESTED_FLOOR,
        soil_foc=api.SOIL_FOC, hydrus_foc=sh.FOC,
        available=dict(rdkit=api.rdkit_available(), hydrus=api.hydrus_available(),
                       ctx=ctx_available()),
    ))


def warmup():
    """Pre-solve the Simple view's default + preset scenarios in the background so a
    first visitor does not wait for them (the slowest part of the old app's first load)."""
    from .schemas import ExposureSpec
    try:
        cases = [("PFOA", 1.0), ("PFOA", 0.1), ("PFOA", 10.0)] + \
                [(c, cw) for _, c, cw, _w, _d in [[k, *v] for k, v in i18n.SCENARIOS_KO.items()]]
        for cong, cw in cases:
            for lang in ("ko", "en"):
                run(Scenario(compound=CompoundSpec(congener=cong),
                             exposure=ExposureSpec(Cwo=cw), lang=lang))
        kw = dict(f_xy_source=api.DEFAULT_F_XY_SOURCE, **crop_kw("site2020"))
        for c in SIMPLE_COMPARE:
            for cw in (0.1, 1.0, 10.0):
                sim("simulate", c, Cwo=cw, season=None, measured_forcing=True, **kw)
    except Exception:                                    # noqa: BLE001 -- best-effort
        pass
