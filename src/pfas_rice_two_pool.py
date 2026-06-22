"""Two-pool root model -- the long-chain breakthrough, promoted to a core component.

The canonical single-mobile-pool core (``pfas_rice_plant_module_4pool_surf``) cannot
reproduce the long chains: one root pool both feeds the xylem AND sets the uptake
gradient, so matching the high measured long-chain root forces the shoot over (or, if
the pool is suppressed, the uptake gradient rises and the root over-inflates). The
single-pool refit hit f_xy=1/L_Ph=1 ceilings ~4-6x under (validation/refit_oryza.py).

This module splits the root into two pools:

  rm  MOBILE  (water + protein; low binding)  -- feeds the xylem and sets jR
  rb  BOUND   (lipid + cell wall; high K_PL/K_cw) -- holds the measured root burden,
              exchanges with the mobile pool at rate k_off (k_on = ratio*k_off*seq;
              ratio = (B_full - B_mobile)/B_mobile is the equilibrium capacity ratio).

Two INDEPENDENT physical levers close the long chains (validation/longchain_closure.py,
long-chain log10 RMSE ~0.08, the breakthrough):
  * a LOW per-congener f_xy  -- strong root retention / low TSCF (Casparian, cell wall);
  * an ENHANCED active carrier Vmax_in -- the high uptake that builds the measured root.
Conflating them (forcing a high f_xy) was what broke the earlier fits.

This is OPT-IN. The canonical 4pool_surf core is unchanged. The saturated per-congener
fit (close_longchain) is structural ADEQUACY (reproduction), NOT a-priori prediction.

Units: BAF (tissue conc / pore-water conc), i.e. Cwo=1 normalization, matching model_api.
"""
from __future__ import annotations

import json
import os
from functools import lru_cache

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from scipy.optimize import brentq

from pfas_rice_plant_module_4pool_surf import Environment, Compound, root_uptake

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
_PAR = json.load(open(os.path.join(_ROOT, "params", "parameters.json")))
_CARR = _PAR["carrier_MichaelisMenten"]
_COMP = _PAR["tissue_composition_recommended"]
_CONG = {c["name"]: c for c in _PAR["congeners"]}
_ENV = Environment()

CONGENERS = [c["name"] for c in _PAR["congeners"]]


def _comps():
    g = _COMP
    return dict(root=g["root"], stem=g["stem"], leaf=g["leaf"], grain=g["grain_brown"])


def _B_full(d, c):
    """Fresh-weight binding factor (basis A): water + protein + lipid + cell wall."""
    return d["theta_fw"] + (1 - d["theta_fw"]) * (
        d["f_prot"] * c["K_prot_Lkg"] + d["f_PL"] * c["K_PL_Lkg"]
        + d["f_cw"] * c["K_cw_wholecw_Lkg"]["root"])


def _B_mobile(d, c):
    """Mobile root pool: water + protein only (no lipid / cell wall) -- does not collapse."""
    return d["theta_fw"] + (1 - d["theta_fw"]) * d["f_prot"] * c["K_prot_Lkg"]


@lru_cache(maxsize=8)
def _drivers(season: float, biomass: str):
    """Interpolated forcings: transpiration Q_TP(t) and per-organ biomass M(t)+dM/dt.

    biomass='oryza' = the mechanistic ORYZA2000 carbon balance (default, matches the
    validated breakthrough); 'growth_rice' = the ORYZA IR72 logistic reconstruction.
    """
    import forcing_rice as fr
    t = np.linspace(0.0, season, int(season) + 1)
    if biomass == "growth_rice":
        import growth_rice as gr
        b = gr.organ_biomass(t)
    else:
        import oryza_growth as og
        b = og.organ_biomass_oryza(t, p=og.OryzaParams(season=season))
    Mr = np.maximum(b["root"], 1e-9)
    Ms = np.maximum(b["stem"], 1e-9)
    Ml = np.maximum(b["leaf"], 1e-9)
    Mg = np.maximum(b["grain"], 1e-4)
    Q = fr.Q_TP(t, season)
    kw = dict(kind="linear", bounds_error=False, fill_value="extrapolate")
    fQ = interp1d(t, Q, **kw)
    fMr, fMs, fMl, fMg = (interp1d(t, x, **kw) for x in (Mr, Ms, Ml, Mg))
    dMr, dMs, dMl, dMg = (interp1d(t, np.gradient(x, t), **kw) for x in (Mr, Ms, Ml, Mg))
    return fQ, fMr, fMs, fMl, fMg, dMr, dMs, dMl, dMg


def simulate(name, f_xy, vmax_in, *, g_xy=0.0, g_ph=None, k_off=0.02, seq=1.0,
             kappa_d=None, L_Ph=None, season=120.0, biomass="oryza"):
    """Run the 2-pool root ODE for one congener; return BAFs (Cwo=1).

    f_xy      : free xylem-loading factor (TSCF analog); LOW for long chains.
    vmax_in   : active-carrier influx capacity (the long-chain enhancement lever).
    g_xy/g_ph : optional B-independent lipid-facilitated xylem/phloem loading.
    seq       : irreversible-sequestration factor on the bound store (1 = equilibrium).
    """
    c = _CONG[name]
    cc = _comps()
    L_Ph = (c.get("L_Ph_oryza") or 0.01) if L_Ph is None else L_Ph
    kappa_d = (c.get("kappa_d_oryza") or 2.0) if kappa_d is None else kappa_d
    g_ph = 0.0 if g_ph is None else g_ph   # optional B-independent lipid phloem term (off by default)
    cmpd = Compound(name=name, K_prot=c["K_prot_Lkg"], K_PL=c["K_PL_Lkg"],
                    K_cw=c["K_cw_wholecw_Lkg"]["root"], kappa_d=kappa_d,
                    Vmax_in=vmax_in, Km_in=_CARR["Km_in"], Vmax_out=_CARR["Vmax_out"],
                    Km_out=_CARR["Km_out"], L_Ph=L_Ph, f_xy=f_xy)
    Bm = _B_mobile(cc["root"], c)
    Bt = _B_full(cc["root"], c)
    ratio = max((Bt - Bm) / Bm, 0.0)
    k_on = ratio * k_off * seq
    Bst = _B_full(cc["stem"], c)
    Blf = _B_full(cc["leaf"], c)
    S_leaf, S_grain = 20.0, 2.0
    fQ, fMr, fMs, fMl, fMg, dMr, dMs, dMl, dMg = _drivers(float(season), biomass)

    def rhs(tt, C):
        rm, rb, st, lf, gr = C
        Q = float(fQ(tt))
        mr, ms, ml, mg = float(fMr(tt)), float(fMs(tt)), float(fMl(tt)), float(fMg(tt))
        mur = float(dMr(tt)) / mr if mr > 0 else 0.0
        mus = float(dMs(tt)) / ms if ms > 0 else 0.0
        mul = float(dMl(tt)) / ml if ml > 0 else 0.0
        mug = float(dMg(tt)) / mg if mg > 0 else 0.0
        Cw_m, Cw_st, Cw_lf = rm / Bm, st / Bst, lf / Blf
        jR = root_uptake(1.0, Cw_m, cmpd, _ENV)
        Cw_xyl = f_xy * Cw_m + g_xy * rm
        A3, A4 = S_leaf * ml, S_grain * mg
        sp = A3 / (A3 + A4) if A3 + A4 > 0 else 0.5
        f3, f4 = sp, 1 - sp
        QPh = max(float(dMg(tt)) * 10.0 + 0.1 * Q, 0.0)
        CPh = L_Ph * Cw_lf + g_ph * lf
        drm = jR - (Q / mr) * Cw_xyl - k_on * rm + k_off * rb + 0.1 * (QPh / mr) * CPh - mur * rm
        drb = k_on * rm - k_off * rb - mur * rb
        dst = (Q / ms) * (Cw_xyl - Cw_st) - mus * st
        dlf = f3 * (Q / ml) * Cw_st - 1.1 * (QPh / ml) * CPh - mul * lf
        dgr = f4 * (Q / mg) * Cw_st + (QPh / mg) * CPh - mug * gr
        return [drm, drb, dst, dlf, dgr]

    sol = solve_ivp(rhs, (0.0, season), np.zeros(5), method="BDF",
                    rtol=1e-5, atol=1e-8, t_eval=[season])
    rm, rb, st, lf, gr = sol.y[:, -1]
    mlf, mst = float(fMl(season)), float(fMs(season))
    straw = (st * mst + lf * mlf) / (mst + mlf)
    return {"root": rm + rb, "straw": straw, "grain": gr, "rm": rm, "rb": rb}


def _carrier_for_root(name, f_xy, obs_root, biomass="oryza", season=120.0):
    base = _CARR["Vmax_in"]

    def g(lm):
        r = simulate(name, f_xy, base * 10 ** lm, biomass=biomass, season=season)["root"]
        return np.log10(max(r, 1e-6)) - np.log10(obs_root)

    lo, hi = -2.0, 7.0
    if g(lo) > 0:
        return base * 10 ** lo
    if g(hi) < 0:
        return base * 10 ** hi
    return base * 10 ** brentq(g, lo, hi, xtol=1e-3, maxiter=80)


def close_longchain(name, obs, *, biomass="oryza", season=120.0):
    """Saturated 3-param structural-adequacy fit (the breakthrough recipe):
    free f_xy -> straw, active carrier -> root, g_ph -> grain (g_xy=0).
    ``obs`` = {'root':.., 'straw':.., 'grain':..} measured BAFs. Returns params + sim."""
    o = obs

    def straw_err(lf):
        f_xy = 10 ** lf
        vm = _carrier_for_root(name, f_xy, o["root"], biomass, season)
        s = simulate(name, f_xy, vm, biomass=biomass, season=season)["straw"]
        return np.log10(max(s, 1e-6)) - np.log10(o["straw"])

    a, b = -6.0, 0.0
    if straw_err(a) * straw_err(b) < 0:
        f_xy = 10 ** brentq(straw_err, a, b, xtol=1e-3, maxiter=60)
    else:
        f_xy = 10 ** a if abs(straw_err(a)) < abs(straw_err(b)) else 10 ** b
    vm = _carrier_for_root(name, f_xy, o["root"], biomass, season)

    def gp_err(lg):
        r = simulate(name, f_xy, vm, g_ph=10 ** lg, biomass=biomass, season=season)["grain"]
        return np.log10(max(r, 1e-6)) - np.log10(o["grain"])

    a2, b2 = -10.0, 2.0
    g_ph = 10 ** brentq(gp_err, a2, b2, xtol=1e-2, maxiter=60) if gp_err(a2) * gp_err(b2) < 0 else (
        10 ** a2 if abs(gp_err(a2)) < abs(gp_err(b2)) else 10 ** b2)
    sim = simulate(name, f_xy, vm, g_ph=g_ph, biomass=biomass, season=season)
    return {"f_xy": f_xy, "vmax_in": vm, "carrier_x": vm / _CARR["Vmax_in"],
            "g_ph": g_ph, "sim": sim}


if __name__ == "__main__":
    import csv
    obs = {}
    with open(os.path.join(_ROOT, "data_obs", "obs_baf_Yamazaki.csv"), newline="") as f:
        for r in csv.DictReader(f):
            obs.setdefault(r["compound"], {})[r["tissue"]] = float(r["baf"])
    print("src/pfas_rice_two_pool.py -- long-chain closure (breakthrough), promoted to src\n")
    print(f"{'PFAS':7}{'f_xy':>8}{'carrier':>9} | {'root p/o':>14}{'straw p/o':>14}{'grain p/o':>14}")
    for nm in ("PFDA", "PFUnDA", "PFDoDA"):
        d = close_longchain(nm, obs[nm])
        r, o = d["sim"], obs[nm]
        print(f"{nm:7}{d['f_xy']:>8.3f}{d['carrier_x']:>8.1f}x | "
              f"{r['root']:6.1f}/{o['root']:<7.1f}{r['straw']:6.1f}/{o['straw']:<7.1f}"
              f"{r['grain']:6.1f}/{o['grain']:<7.1f}")
