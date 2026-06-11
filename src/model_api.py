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

with open(os.path.join(_ROOT, "params", "parameters.json")) as _f:
    PARAMS = json.load(_f)
_CARR = PARAMS["carrier_MichaelisMenten"]
_COMP = PARAMS["tissue_composition_recommended"]
_CONG = {c["name"]: c for c in PARAMS["congeners"]}

CONGENERS = [c["name"] for c in PARAMS["congeners"]]            # 12, ordered
TISSUES = ("root", "stem", "leaf", "grain")


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


def simulate(congener="PFOA", Cwo=1.0, E_m_mV=-120.0, f_xy_source="recommended",
             f_xy_override=None, season=120.0, n_t=241, measured_forcing=True):
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

    Returns a dict with t, per-compartment conc & BAF time series, finals, straw,
    B_k, and the effective parameters used.
    """
    if congener not in _CONG:
        raise KeyError(f"unknown congener {congener!r}; known: {CONGENERS}")
    c = _CONG[congener]
    t = np.linspace(0.0, season, n_t)

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
    inputs = PlantInputs(t=t, Cwo=np.full_like(t, float(Cwo)), Qtp=Qtp, M=M)

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
                    L_Ph=L_Ph, f_xy=f_xy)
    comps = _compartments()
    env = Environment(E=E_m_mV / 1000.0)
    model = RiceUptakeModel(env=env, cmpd=cmpd, comps=comps, inputs=inputs)
    sol = model.solve(t)
    C = sol.y                                            # (4, n_t)
    Mf = inputs.M_(t[-1])
    straw = (C[STEM] * Mf[STEM] + C[LEAF] * Mf[LEAF]) / (Mf[STEM] + Mf[LEAF])
    B = binding_factors(comps, cmpd)
    cwo = float(Cwo) if Cwo else 1.0
    return dict(
        t=t, congener=congener, success=bool(sol.success),
        conc={k: C[i] for i, k in enumerate(TISSUES)},
        straw=straw,
        baf_final={k: float(C[i, -1] / cwo) for i, k in enumerate(TISSUES)},
        straw_baf=float(straw[-1] / cwo),
        B_k={k: float(B[i]) for i, k in enumerate(TISSUES)},
        N=float(env.N), eN=float(np.exp(env.N)),
        params=dict(f_xy=f_xy, L_Ph=L_Ph, kappa_d=kappa_d, K_PL=c["K_PL_Lkg"],
                    K_prot=c["K_prot_Lkg"], K_cw=c["K_cw_wholecw_Lkg"]["root"],
                    n_C=c["n_C"], group=c["group"]),
    )


if __name__ == "__main__":
    r = simulate("PFOA")
    print("PFOA:", {k: round(v, 3) for k, v in r["baf_final"].items()},
          "straw_baf", round(r["straw_baf"], 3), "B_k", {k: round(v, 2) for k, v in r["B_k"].items()})
