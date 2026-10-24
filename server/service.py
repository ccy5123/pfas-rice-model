"""Scenario -> ``model_api`` calls -> JSON-ready payloads.

This is the web backend's engine. It reproduces what the Streamlit sidebar +
``ui/common.run_model`` do -- build the exposure drivers for the selected data
source, run the selected compound class, attach observed/biomonitoring BAFs --
but as plain functions with a shared, de-duplicating memo, so:

* the browser handles every purely visual interaction (tab, day slider, colour
  metric) without calling back here;
* a solve is computed once and then served to every view (and every user) that
  asks for the same inputs;
* the slow extras (congener comparison, Tang check, two-pool overlays, the
  Bayesian inverse) run only when their view is actually opened.

Deliberate difference from the Streamlit app: driver dicts are passed through
WHOLE. The Streamlit cache packs drivers into a tuple of (t, Cwo, Qtp, M) and so
drops the ORYZA ``leaf_loss`` senescence rate that ``drivers_from_arrays`` /
``drivers_from_tables`` attach; here it reaches the model, as model_api intends.
"""
from __future__ import annotations

import inspect
import io
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
APP_VERSION = "2.0 (web: FastAPI + React)"
REPO_URL = "https://github.com/ccy5123/pfas-rice-model"
DOCS_URL = "https://github.com/ccy5123/pfas-rice-model/tree/main/docs"
ORGANS4 = ("root", "stem", "leaf", "grain")

try:
    from chem_lookup import KOC_MAX_LKG as KOC_MAX, F_N_TESTED_FLOOR
except Exception:                                   # noqa: BLE001 -- optional module
    KOC_MAX, F_N_TESTED_FLOOR = 1e9, 0.1

SIM = Memo(256)     # individual model solves (shared across scenarios and users)
RUN = Memo(64)      # whole scenario runs
DRV = Memo(32)      # driver builds (a live HYDRUS-1D run takes seconds)
FIG = Memo(256)     # serialised figures
AUX = Memo(128)     # Tang check, two-pool overlays, inverse estimates, lookups


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
def karickhoff(log_kow):
    return float(LP.koc_neutral(float(log_kow)))


def neutral_koc(n: NeutralSpec):
    return float(n.Koc) if n.Koc is not None else karickhoff(n.log_kow)


def neutral_kwargs(n: NeutralSpec):
    """``simulate_neutral`` kwargs -- the dict ``ui/sidebar._neutral_panel`` returns
    (without ``Koc``, which only the soil-side exposure modes use)."""
    kw = dict(log_kow=float(n.log_kow), name=(n.name or "neutral"),
              half_life=(float(n.half_life) or None), tscf_model=n.tscf_model,
              lipid_source=n.lipid_source, phloem=bool(n.phloem), air=bool(n.air),
              g_apo=float(n.g_apo))
    if n.air:
        kw.update(MW=float(n.MW), K_AW=float(n.K_AW), air_kw=dict(C_air=float(n.C_air)))
    if n.weak_electrolyte:
        kw.update(pKa=float(n.pKa), pH=float(n.pH), is_acid=bool(n.is_acid))
    return kw


def sim_kw(m: MechanismSpec):
    """The PFAS ``simulate`` kwargs the Streamlit ``run_model`` builds."""
    kw = dict(E_m_mV=float(m.E_m), f_xy_source=m.fxy_source, biomass=m.biomass,
              uptake=m.uptake, lipid_loading=bool(m.lipid_loading),
              km_scale=float(m.km_scale or 1.0))
    if m.vmax_scale is not None:
        kw["vmax_scale"] = float(m.vmax_scale)
    if m.g_apo is not None:
        kw["g_apo"] = float(m.g_apo)
    return kw


def default_k_leach(comp: CompoundSpec):
    """The flooded-shape leaching default the Streamlit slider is seeded with."""
    cong = comp.congener if comp.kind == "congener" else None
    koc = neutral_koc(comp.neutral) if comp.kind == "neutral" else None
    return float(api.default_k_leach(cong, Koc=koc))


# ---- editable tables (pure versions of ui/common's DataFrame helpers) ----
def default_growth_table(season, biomass, n_rows=7):
    t = np.linspace(0.0, float(season), n_rows)
    b = api._biomass_fn(biomass)(t, float(season))
    out = {"day": np.round(t, 0).tolist()}
    for o in ORGANS4:
        out[o] = np.round(np.asarray(b[o], float) * 1e3, 2).tolist()
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


def _hydrus_drivers(congener, season, f_oc, flood_until, percolation, biomass, Kd=None):
    def go():
        drv, _ = api.hydrus_drivers(congener, season=season, f_oc=f_oc, biomass=biomass, Kd=Kd,
                                    flood_until=float(flood_until),
                                    percolation=float(percolation))
        return drv
    return DRV.get(key_of("hydrus", congener, season, f_oc, flood_until, percolation,
                          biomass, Kd), go)


# ---------------------------------------------------------------- drivers per data source
@dataclass
class Built:
    drivers: dict | None = None
    soil_obj: object = None
    measured_bio: dict | None = None
    messages: list = field(default_factory=list)
    implied_volume_ml: dict | None = None
    # the effective "built-in forcing" settings (used when drivers is None)
    Cwo: float = 1.0
    season: float = 120.0
    measured: bool = True
    cwo_profile: str = "constant"
    k_leach: float | None = None


def build(scn: Scenario) -> Built:
    """Exposure drivers for the selected data source -- ``ui/sidebar.build`` minus widgets."""
    e, comp, lang = scn.exposure, scn.compound, scn.lang
    biomass = scn.mechanism.biomass
    b = Built()

    def msg(level, text):
        b.messages.append({"level": level, "text": text})

    # Only the parametric and custom-table sources expose these; the others run the
    # Streamlit sidebar defaults (Cwᵒ 1, 120 d, measured forcings, flat Cwᵒ).
    if e.mode in ("parametric", "tables"):
        b.Cwo, b.season = float(e.Cwo), float(e.season)
        b.measured, b.cwo_profile = bool(e.measured_forcing), e.cwo_profile
    if b.cwo_profile == "flooded":
        b.k_leach = float(e.k_leach) if e.k_leach is not None else default_k_leach(comp)

    if e.mode == "tables":
        t = e.tables
        try:
            growth = clean_table(t.growth or default_growth_table(e.season, biomass), ORGANS4)
            cwo = clean_table(t.cwo or default_cwo_table(e.season, e.Cwo), ["Cwo"])
            drv = api.drivers_from_tables(growth, cwo, growth_units=t.growth_units,
                                          Cwo_const=float(e.Cwo), biomass=biomass)
            b.drivers = drv
            Mf = np.asarray(drv["M"], float)[-1]
            b.implied_volume_ml = {o: float(Mf[i] / max(float(t.density.get(o, 1.0)), 1e-6) * 1e3)
                                   for i, o in enumerate(ORGANS4)}
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
                b.drivers = api.drivers_from_arrays(tt, Cwo, Qtp=Q, M=M, biomass=biomass)
                msg("success", f"Loaded {len(tt)} rows from the uploaded CSV.")
            elif e.csv_text is None:
                b.drivers = DRV.get(key_of("example_csv", biomass), lambda: api.load_driver_csv(
                    os.path.join(EX, "hydrus_drivers_example.csv"), biomass=biomass))
                msg("info", "Using examples/hydrus_drivers_example.csv (synthetic HYDRUS-style run).")
            else:
                msg("info", "No driver CSV — running the default parametric scenario.")
        except Exception as ex:                          # noqa: BLE001
            msg("error", f"Could not read drivers: {ex}")

    elif e.mode == "hydrus":
        b.season = 120.0
        if not api.hydrus_available():
            msg("warning", "HYDRUS-1D engine not built in this environment yet — build it once "
                           "(compiles the FORTRAN solver with gfortran; ~1 min). Until built, the "
                           "tool falls back to the parametric model.")
        elif comp.kind == "neutral":
            kd = neutral_koc(comp.neutral) * float(e.f_oc)
            try:
                b.drivers = _hydrus_drivers(None, 120.0, e.f_oc, e.flood_until, e.percolation,
                                            biomass, Kd=kd)
                msg("success", f"HYDRUS-1D run complete — Kd = {kd:.3g} L/kg "
                               f"(Koc {neutral_koc(comp.neutral):.1f} × f_oc {e.f_oc:g}; "
                               f"mean-normalised).")
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
                    b.drivers = _hydrus_drivers(soil_cong, 120.0, e.f_oc, e.flood_until,
                                                e.percolation, biomass)
                    msg("success", f"HYDRUS-1D run complete — Cwᵒ(t) for {soil_cong} "
                                   "(Kd-retarded; mean-normalised).")
                except Exception as ex:                  # noqa: BLE001
                    msg("error", f"HYDRUS run failed: {ex}")

    elif e.mode == "inventory":
        season = 120.0
        tt = np.linspace(0.0, season, 241)
        flooded = tt < (0.75 * season) if e.inv_flood else None
        if e.inv_flood:
            kl = float(e.inv_k_leach) if e.inv_k_leach is not None else default_k_leach(comp)
        else:
            kl = 0.0
        Cwo, soil = api.pore_water_from_inventory(tt, e.C_total, K_F=e.K_F, n=e.n_F,
                                                  theta_g=e.theta_g, flooded=flooded, k_leach=kl)
        b.drivers = api.drivers_from_arrays(tt, Cwo, season=season, biomass=biomass)
        b.soil_obj = soil

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
    desc: object = None
    provisional: bool = False


def _cwo_kw(b: Built):
    return {"k_leach": float(b.k_leach)} if b.cwo_profile != "constant" and b.k_leach is not None else None


def _run(scn: Scenario) -> Run:
    b = build(scn)
    comp, m = scn.compound, scn.mechanism
    skw = sim_kw(m)
    desc, provisional = None, False
    if comp.kind == "neutral":
        nk = neutral_kwargs(comp.neutral)
        log_kow = nk.pop("log_kow")
        if b.drivers is not None:
            res = sim("simulate_neutral", log_kow, drivers=b.drivers, biomass=m.biomass, **nk)
        else:
            res = sim("simulate_neutral", log_kow, Cwo=b.Cwo, season=b.season, biomass=m.biomass,
                      measured_forcing=b.measured, cwo_profile=b.cwo_profile, cwo_kw=_cwo_kw(b),
                      Koc=neutral_koc(comp.neutral), **nk)
    elif comp.kind == "smiles":
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
                res = sim("simulate_from_smiles", smiles, Cwo=b.Cwo, season=b.season,
                          measured_forcing=b.measured, cwo_profile=b.cwo_profile,
                          cwo_kw=_cwo_kw(b), **skw)
        except ApiError:
            raise
        except Exception as ex:                          # noqa: BLE001
            raise ApiError(400, "Could not build a compound from that SMILES — check the "
                                f"structure. ({ex})")
        desc = res.get("descriptors")
        provisional = bool(res.get("provisional", False))
    else:
        if comp.congener not in api.CONGENERS:
            raise ApiError(400, f"unknown congener {comp.congener!r}; known: {api.CONGENERS}")
        if b.drivers is not None:
            res = sim("simulate", comp.congener, drivers=b.drivers, **skw)
        else:
            res = sim("simulate", comp.congener, Cwo=b.Cwo, season=b.season,
                      measured_forcing=b.measured, cwo_profile=b.cwo_profile,
                      cwo_kw=_cwo_kw(b), **skw)
    congener = res["congener"]
    obs = {} if comp.kind == "neutral" else api.observed_baf(congener)
    mb = b.measured_bio
    bio_baf = (api.baf_from_measurement(mb["conc"], mb["Cwo"])
               if mb and mb.get("Cwo") else None)
    return Run(scn=scn, built=b, res=res, congener=congener, obs=obs, bio_baf=bio_baf,
               sim_kw=skw, desc=desc, provisional=provisional)


def run(scn: Scenario) -> Run:
    return RUN.get(key_of("run", scn), lambda: _run(scn))


def payload(r: Run) -> dict:
    """Everything the front end needs to draw the headline, the plant map (scrubbed
    client-side over the whole season) and the series -- one compact JSON object."""
    res, b = r.res, r.built
    root_c = float(res["conc"]["root"][-1])
    straw_c = float(res["straw"][-1])
    grain_c = float(res["conc"]["grain"][-1])
    M = np.asarray(res["M"], float)
    schematic = {}
    for metric in ("conc", "baf"):
        ms = api.metric_series(res, metric)
        schematic[metric] = {"cmin": ms["cmin"], "cmax": ms["cmax"], "label": ms["label"]}
    p = dict(res.get("params", {}))
    out = dict(
        congener=r.congener, kind=r.scn.compound.kind, mode=r.scn.exposure.mode,
        success=bool(res.get("success", True)),
        t=res["t"], conc=res["conc"], baf=res["baf"], straw=res["straw"],
        Cwo=res["Cwo"], Qtp=res["Qtp"], M={k: M[:, i] for i, k in enumerate(api.TISSUES)},
        baf_final=res["baf_final"], straw_baf=res["straw_baf"], cwo_ref=res.get("cwo_ref"),
        B_k=res.get("B_k"), N=res.get("N"), eN=res.get("eN"), params=p,
        finals=dict(root=root_c, straw=straw_c, grain=grain_c),
        bands={k: api.predictive_band(v) for k, v in
               (("root", root_c), ("straw", straw_c), ("grain", grain_c))},
        uncertainty_factor=api.uncertainty_factor(),
        intake=api.intake_fraction(grain_c, congener=r.congener),
        schematic=schematic, obs=r.obs, bio_baf=r.bio_baf, measured_bio=b.measured_bio,
        descriptors=(asdict(r.desc) if r.desc is not None else None),
        provisional=r.provisional, messages=b.messages,
        implied_volume_ml=b.implied_volume_ml, has_isotherm=b.soil_obj is not None,
        effective=dict(Cwo=b.Cwo, season=b.season, measured=b.measured,
                       cwo_profile=b.cwo_profile, k_leach=b.k_leach,
                       drivers=b.drivers is not None),
    )
    if r.scn.compound.kind == "neutral":
        out.update(TSCF=res.get("TSCF"), K_PW=res.get("K_PW"), rcf_briggs=res.get("rcf_briggs"),
                   half_life=res.get("half_life"), warning=res.get("warning"),
                   air_summary=res.get("air_summary"), tf_final=res.get("tf_final"),
                   koc=neutral_koc(r.scn.compound.neutral))
    return clean(out)


# ---------------------------------------------------------------- figures
def _fig_out(fig, data=None):
    """'{"figure": <plotly JSON>, "data": {...}}' -- figure JSON straight from plotly
    (its encoder already maps NaN to null)."""
    return '{"figure":' + fig.to_json() + ',"data":' + dumps(data or {}) + '}'


def _twopool(kind, congener):
    def go():
        try:
            r = (api.simulate_twopool_seq(congener) if kind == "seq"
                 else api.simulate_twopool_nstem(congener))
            return {"root": r["baf_final"]["root"], "straw": r["straw_baf"],
                    "grain": r["baf_final"]["grain"]}
        except Exception:                                # noqa: BLE001
            return None
    return AUX.get(key_of("twopool", kind, congener), go)


def tang(congener, dose, biomass, use_refit=False, lipid=False):
    return AUX.get(key_of("tang", congener, dose, biomass, use_refit, lipid),
                   lambda: api.tang_tf_validation(congener, dose=dose, biomass=biomass,
                                                  use_refit=use_refit, lipid_loading=lipid))


SIMPLE_COMPARE = ("PFBA", "PFHxA", "PFOA", "PFDA", "PFDoDA", "PFOS")


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
    if name == "baf":
        extra = None
        if comp.kind == "congener" and lang == "en":
            over = {}
            if o.get("twopool_seq"):
                tp = _twopool("seq", r.congener)
                if tp is not None:
                    over["two-pool (seq, exploratory)"] = tp
            if o.get("twopool_nstem"):
                mg = _twopool("nstem", r.congener)
                if mg is not None:
                    over["two-pool + redistributed shoot (merged)"] = mg
            extra = over or None
        return _fig_out(plots.fig_baf(res, obs, extra=extra, lang=lang),
                        {"overlays": sorted(extra) if extra else []})
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
            else:
                results[nm] = sim("simulate", nm, Cwo=b.Cwo, season=b.season,
                                  measured_forcing=b.measured, **r.sim_kw)
        return _fig_out(plots.fig_compare(results, tissue))
    if name == "congener_compare":                        # Simple: six representative
        b = r.built
        m = req.scenario.mechanism
        kw = dict(E_m_mV=float(m.E_m), f_xy_source=m.fxy_source, biomass=m.biomass)
        rep = [c for c in SIMPLE_COMPARE if c in api.CONGENERS]
        results = {c: sim("simulate", c, Cwo=b.Cwo, season=b.season,
                          measured_forcing=b.measured, **kw) for c in rep}
        return _fig_out(plots.fig_congener_compare(results, lang=lang, order=rep))
    if name == "tang_tf":
        if comp.kind != "congener" or r.congener not in api.TANG_CONGENERS:
            raise ApiError(400, f"Tang 2026 covers {', '.join(api.TANG_CONGENERS)} only")
        dose = "low" if o.get("dose") == "low" else "mean"
        bm = "oryza" if o.get("biomass") == "oryza" else "growth_rice"
        val = tang(r.congener, dose, bm)
        valr = tang(r.congener, dose, bm, use_refit=True)
        vall = tang(r.congener, dose, bm, lipid=True) if o.get("lipid") else None
        return _fig_out(plots.fig_tang_tf(val, valr, vall),
                        {"f_xy": val["f_xy"], "f_xy_refit": valr["f_xy"]})
    if name == "intake_gauge":
        grain_c = float(res["conc"]["grain"][-1])
        rice = float(o.get("rice_g", api.DEFAULT_RICE_INTAKE_G_PER_DAY))
        bw = float(o.get("bw", api.DEFAULT_BODY_WEIGHT_KG))
        info = api.intake_fraction(grain_c, congener=r.congener, rice_intake_g_day=rice,
                                   body_weight_kg=bw)
        band = api.predictive_band(grain_c)
        lo = api.intake_fraction(band["lo"], congener=r.congener, rice_intake_g_day=rice,
                                 body_weight_kg=bw)
        hi = api.intake_fraction(band["hi"], congener=r.congener, rice_intake_g_day=rice,
                                 body_weight_kg=bw)
        return _fig_out(plots.fig_intake_gauge(info, lang=lang),
                        {"info": info, "lo_percent": lo["percent"], "hi_percent": hi["percent"]})
    raise ApiError(404, f"unknown figure {name!r}")


def figure(req: FigureRequest) -> str:
    return FIG.get(key_of("fig", req), lambda: _figure(req))


def cwo_preview(req: CwoPreviewRequest) -> str:
    """The flooded-shape preview next to the k_leach slider. For a NEUTRAL compound the
    shape comes from its own Koc (the PFAS chain-length QSPR does not apply)."""
    comp = req.compound
    kl = float(req.k_leach) if req.k_leach is not None else default_k_leach(comp)

    def go():
        if comp.kind != "neutral":
            cong = comp.congener if comp.kind == "congener" else None
            return _fig_out(plots.fig_cwo_profile(cong, level=req.level, profile="flooded",
                                                  season=req.season, k_leach=kl))
        import plotly.graph_objects as go_
        t = np.linspace(0.0, float(req.season), 121)
        flat = api.cwo_profile_series(t, req.level, "constant")
        shaped = api.cwo_profile_series(t, req.level, "flooded", k_leach=kl,
                                        Koc=neutral_koc(comp.neutral))
        fig = go_.Figure()
        fig.add_scatter(x=t, y=flat, name="constant",
                        line=dict(color="#bbb", width=1.5, dash="dash"))
        fig.add_scatter(x=t, y=shaped, name="flooded", line=dict(color="#1f77b4", width=2.5),
                        fill="tonexty", fillcolor="rgba(31,119,180,0.10)")
        fig.update_layout(title=f"Cwᵒ(t) preview — flooded (mean={req.level:g})",
                          hovermode="x unified", height=230, showlegend=True,
                          legend=dict(orientation="h", y=1.12, x=0), **plots._LAYOUT)
        fig.update_yaxes(title_text="Cwᵒ [µg/L]", rangemode="tozero")
        fig.update_xaxes(title_text="days after transplant")
        return _fig_out(fig)
    return FIG.get(key_of("cwo_preview", req, kl), go)


# ---------------------------------------------------------------- Bayesian inverse
def _inverse_args(req: InverseRequest):
    scn, lang = req.scenario, req.lang
    comp = scn.compound
    meas = {k: float(v) for k, v in req.measured.items()
            if k in ("root", "straw", "grain") and v is not None and float(v) > 0}
    if not meas:
        raise ApiError(400, i18n.t("inv.enter_first", lang))
    m = scn.mechanism
    base = dict(sigma_log10=float(req.sigma_log10), E_m_mV=float(m.E_m),
                f_xy_source=m.fxy_source, biomass=m.biomass)
    if comp.kind == "neutral":
        nk = neutral_kwargs(comp.neutral)
        label = nk["name"]
        mech = dict(log_kow=float(nk["log_kow"]),
                    neutral_kw={k: v for k, v in nk.items() if k not in ("log_kow", "name")})
        return label, meas, {**base, **mech}
    if comp.kind != "congener" or comp.congener not in api.CONGENERS:
        raise ApiError(400, i18n.t("inv.not_curated", lang))
    mech = {k: v for k, v in sim_kw(m).items() if k not in ("E_m_mV", "f_xy_source", "biomass")}
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
    if kind == "summary_md":
        grain_c = float(res["conc"]["grain"][-1])
        intake = api.intake_fraction(grain_c, congener=cong)
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
def tables_defaults(season, biomass, Cwo):
    return dict(growth=default_growth_table(season, biomass),
                cwo=default_cwo_table(season, Cwo),
                density=api.DEFAULT_TISSUE_DENSITY)


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


_UNC = {"en": [("Typical (±~40%)", 0.15), ("High precision (±~20%)", 0.10),
               ("Rough (±~2×)", 0.30)],
        "ko": [("보통 (±~40%)", 0.15), ("고정밀 (±~20%)", 0.10), ("대략 (±~2배)", 0.30)]}

EXAMPLE_SMILES = [
    ("PFOA  (known)", "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"),
    ("PFOS  (known)", "OS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"),
    ("GenX / HFPO-DA  (known ether)", "OC(=O)C(F)(OC(F)(F)C(F)(F)C(F)(F)F)C(F)(F)F"),
    ("PFTrDA — NOVEL C13 PFCA", "OC(=O)" + "C(F)(F)" * 11 + "C(F)(F)F"),
    ("ADONA-like — NOVEL ether-PFCA", "OC(=O)C(F)(F)OC(F)(F)C(F)(F)OC(F)(F)C(F)(F)F"),
]


def meta():
    return clean(dict(
        app_version=APP_VERSION, repo_url=REPO_URL, docs_url=DOCS_URL,
        congeners=api.CONGENERS,
        congener_info={n: {"n_C": c["n_C"], "group": c["group"]} for n, c in api._CONG.items()},
        labels=i18n.CONGENER_LABELS,
        presets={"en": [[k, v] for k, v in i18n.PRESETS["en"].items()],
                 "ko": [[k, v[0], v[1]] for k, v in i18n.PRESETS["ko"].items()]},
        scenarios_ko=[[k, *v] for k, v in i18n.SCENARIOS_KO.items()],
        strings=i18n.STRINGS, organ_labels=i18n.ORGAN_LABELS,
        uptake_modes=api.UPTAKE_MODES, default_uptake=api.DEFAULT_UPTAKE,
        tang_congeners=list(api.TANG_CONGENERS), tang_refit_fxy=api.TANG_REFIT_FXY,
        example_smiles=EXAMPLE_SMILES, tissue_density=api.DEFAULT_TISSUE_DENSITY,
        intake=dict(rice_g=api.DEFAULT_RICE_INTAKE_G_PER_DAY, bw=api.DEFAULT_BODY_WEIGHT_KG,
                    twi=api.EFSA_TWI_NG_PER_KG_BW_WEEK,
                    congeners=list(api.EFSA_TWI_CONGENERS)),
        uncertainty_factor=api.uncertainty_factor(),
        default_k_leach={c: float(api.default_k_leach(c)) for c in api.CONGENERS},
        karickhoff=dict(slope=LP.KOC_NEUTRAL_SLOPE, intercept=LP.KOC_NEUTRAL_INTERCEPT),
        koc_max=KOC_MAX, f_n_tested_floor=F_N_TESTED_FLOOR,
        inverse_precision=_UNC,
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
        for c in SIMPLE_COMPARE:
            for cw in (0.1, 1.0, 10.0):
                sim("simulate", c, Cwo=cw, season=120.0, measured_forcing=True,
                    E_m_mV=-120.0, f_xy_source="recommended", biomass="oryza")
    except Exception:                                    # noqa: BLE001 -- best-effort
        pass
