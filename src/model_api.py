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
    binding_factors, _logistic, ROOT, STEM, LEAF, FRUIT)
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


def _default_drivers(t, season, Cwo, measured_forcing):
    """Build the (Cwo, Qtp, M) driver arrays for the built-in scenarios."""
    if measured_forcing:
        Qtp = fr.Q_TP(t, season)
        b = gr.organ_biomass(t, season)
        M = np.maximum(np.column_stack([b["root"], b["stem"], b["leaf"], b["grain"]]), 1e-4)
    else:
        Qtp = 0.05 + 0.35 * np.exp(-((t - 0.62 * season) ** 2) / (2 * (season / 4.8) ** 2))
        M = np.column_stack([_logistic(t, 1e-3, 0.030, 0.10, season * 0.17),
                             _logistic(t, 1e-3, 0.040, 0.10, season * 0.21),
                             _logistic(t, 1e-3, 0.050, 0.12, season * 0.25),
                             _logistic(t, 1e-5, 0.025, 0.18, season * 0.67)])
    Cwo_series = np.full_like(t, float(Cwo))
    return Cwo_series, Qtp, M


def simulate(congener="PFOA", Cwo=1.0, E_m_mV=-120.0, f_xy_source="recommended",
             f_xy_override=None, season=120.0, n_t=241, measured_forcing=True,
             drivers=None, K_surf=0.0):
    """Run the 4-compartment ODE for one congener and scenario.

    Parameters
    ----------
    congener : one of CONGENERS.
    Cwo : constant pore-water free concentration C_w^o [ug/L] (1.0 -> conc == BAF).
    E_m_mV : root membrane potential [mV] (GHK anion-exclusion lever).
    f_xy_source : "recommended" (monotone) or "W2fit" (reproduces Yamazaki).
    f_xy_override : if given, use this f_xy instead.
    measured_forcing : True -> Q_TP from forcing_rice, M from growth_rice (ORYZA);
        False -> the illustrative logistic placeholders.
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
    if congener not in _CONG:
        raise KeyError(f"unknown congener {congener!r}; known: {CONGENERS}")
    c = _CONG[congener]

    if drivers is not None:
        t = np.asarray(drivers["t"], dtype=float)
        Cwo_series = np.asarray(drivers["Cwo"], dtype=float)
        Qtp = np.asarray(drivers["Qtp"], dtype=float)
        M = np.asarray(drivers["M"], dtype=float)
        season = float(t[-1])
    else:
        t = np.linspace(0.0, season, n_t)
        Cwo_series, Qtp, M = _default_drivers(t, season, Cwo, measured_forcing)
    inputs = PlantInputs(t=t, Cwo=Cwo_series, Qtp=Qtp, M=M)

    if f_xy_override is not None:
        f_xy = float(f_xy_override)
    else:
        f_xy = c["f_xy_recommended"] if f_xy_source == "recommended" else (c.get("f_xy_W2fit") or c["f_xy_recommended"])
    kappa_d = c.get("kappa_d_W2fit") or 2.0
    L_Ph = c.get("L_Ph_W2fit") or 0.01
    cmpd = Compound(name=congener, K_prot=c["K_prot_Lkg"], K_PL=c["K_PL_Lkg"],
                    K_cw=c["K_cw_wholecw_Lkg"]["root"], kappa_d=kappa_d,
                    Vmax_in=_CARR["Vmax_in"], Km_in=_CARR["Km_in"],
                    Vmax_out=_CARR["Vmax_out"], Km_out=_CARR["Km_out"],
                    L_Ph=L_Ph, f_xy=f_xy, K_surf=float(K_surf))
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
        params=dict(f_xy=f_xy, L_Ph=L_Ph, kappa_d=kappa_d, K_PL=c["K_PL_Lkg"],
                    K_prot=c["K_prot_Lkg"], K_cw=c["K_cw_wholecw_Lkg"]["root"],
                    K_surf=float(K_surf), n_C=c["n_C"], group=c["group"]),
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
