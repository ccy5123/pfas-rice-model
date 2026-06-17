#!/usr/bin/env python3
# =============================================================================
# validation/twopool_longchain.py
# -----------------------------------------------------------------------------
# Prototype 2-POOL root test for the long-chain (C10-C12) tradeoff (LC4).
#
# LC3 (refuted) showed the SINGLE-pool lipid term reproduces long-chain shoot only
# by DRAINING the root (Cw_xyl = g_xy*C_root subtracts the whole root burden). The
# literature (Chen2025, evi-lc-litread) says the membrane/lipid pool keeps rising
# with chain length while protein peaks at C6-C10 -> the long-chain root burden is a
# large, slowly-exchanging LIPID/CELL-WALL bound store, distinct from the small
# mobile pool that actually feeds the xylem.
#
# This prototype splits the ROOT into two sub-pools:
#   rm  mobile (water + protein; low binding B_m) -- feeds the xylem and exchanges
#       with soil (jR) and with the bound store.
#   rb  bound  (lipid + cell-wall; the K_PL/K_cw part) -- a slow store (rate k_off)
#       that holds most of the MEASURED root burden but is not loaded to the xylem.
# measured root = rm + rb; equilibrium rb/rm = (B_total - B_m)/B_m (the binding split).
# Lipid-facilitated loading draws from the mobile pool (g_xy*rm), so the shoot is fed
# WITHOUT subtracting the large bound store -> the test of whether root AND shoot can
# be matched simultaneously (closing LC3). Mass-conserving (sole source jR + phloem).
#
#   python validation/twopool_longchain.py
# =============================================================================
import csv, os, sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
import json
from pfas_rice_plant_module_4pool_surf import Environment, Compound, Compartment, root_uptake  # noqa
import oryza_growth as og, forcing_rice as fr  # noqa

PAR = json.load(open(os.path.join(ROOT, "params", "parameters.json")))
carr = PAR["carrier_MichaelisMenten"]; comp = PAR["tissue_composition_recommended"]
obs = {}
with open(os.path.join(ROOT, "data_obs", "obs_baf_Yamazaki.csv"), newline="") as f:
    for r in csv.DictReader(f):
        obs.setdefault(r["compound"], {})[r["tissue"]] = float(r["baf"])
CONG = {c["name"]: c for c in PAR["congeners"]}

SEASON = 120.0
t = np.linspace(0.0, SEASON, 121)
_b = og.organ_biomass_oryza(t, p=og.OryzaParams(season=SEASON))
Mr = np.maximum(_b["root"], 1e-9); Ms = np.maximum(_b["stem"], 1e-9)
Ml = np.maximum(_b["leaf"], 1e-9); Mg = np.maximum(_b["grain"], 1e-4)
Qtp = fr.Q_TP(t, SEASON)
from scipy.interpolate import interp1d
kw = dict(kind="linear", bounds_error=False, fill_value="extrapolate")
fQ = interp1d(t, Qtp, **kw)
fMr, fMs, fMl, fMg = (interp1d(t, x, **kw) for x in (Mr, Ms, Ml, Mg))
dMr, dMs, dMl, dMg = (interp1d(t, np.gradient(x, t), **kw) for x in (Mr, Ms, Ml, Mg))
env = Environment()
RM, RB, ST, LF, GR = 0, 1, 2, 3, 4


def _comps():
    g = comp
    return dict(root=g["root"], stem=g["stem"], leaf=g["leaf"], grain=g["grain_brown"])


def _Bfull(d, c):
    return d["theta_fw"] + (1 - d["theta_fw"]) * (
        d["f_prot"] * c["K_prot_Lkg"] + d["f_PL"] * c["K_PL_Lkg"] + d["f_cw"] * c["K_cw_wholecw_Lkg"]["root"])


def _Bmobile(d, c):
    # mobile root pool: water + protein only (no lipid / cell wall) -> not collapsed
    return d["theta_fw"] + (1 - d["theta_fw"]) * d["f_prot"] * c["K_prot_Lkg"]


def simulate2(name, f_xy, g_xy, g_ph, k_off=0.02, L_Ph=None, kappa_d=None):
    c = CONG[name]; cc = _comps()
    L_Ph = c.get("L_Ph_oryza") or 0.01 if L_Ph is None else L_Ph
    kappa_d = c.get("kappa_d_oryza") or 2.0 if kappa_d is None else kappa_d
    cmpd = Compound(name=name, K_prot=c["K_prot_Lkg"], K_PL=c["K_PL_Lkg"],
                    K_cw=c["K_cw_wholecw_Lkg"]["root"], kappa_d=kappa_d,
                    Vmax_in=carr["Vmax_in"], Km_in=carr["Km_in"],
                    Vmax_out=carr["Vmax_out"], Km_out=carr["Km_out"], L_Ph=L_Ph, f_xy=f_xy)
    Bm = _Bmobile(cc["root"], c); Bt = _Bfull(cc["root"], c)
    ratio = max((Bt - Bm) / Bm, 0.0); k_on = ratio * k_off
    Bst = _Bfull(cc["stem"], c); Blf = _Bfull(cc["leaf"], c)
    # stem/leaf use full B (single pool); grain via phloem
    S_leaf, S_grain = 20.0, 2.0

    def rhs(tt, C):
        rm, rb, st, lf, gr = C
        Q = float(fQ(tt)); mr, ms, ml, mg = float(fMr(tt)), float(fMs(tt)), float(fMl(tt)), float(fMg(tt))
        mur = float(dMr(tt))/mr if mr > 0 else 0.0; mus = float(dMs(tt))/ms if ms > 0 else 0.0
        mul = float(dMl(tt))/ml if ml > 0 else 0.0; mug = float(dMg(tt))/mg if mg > 0 else 0.0
        Cw_m = rm / Bm; Cw_st = st / Bst; Cw_lf = lf / Blf
        jR = root_uptake(1.0, Cw_m, cmpd, env)               # Cwo=1
        Cw_xyl = f_xy * Cw_m + g_xy * rm                     # mobile-pool loading
        A3 = S_leaf*ml; A4 = S_grain*mg; sp = A3/(A3+A4) if A3+A4>0 else .5; f3, f4 = sp, 1-sp
        QPh = max(float(dMg(tt))*10.0 + 0.1*Q, 0.0); CPh = L_Ph*Cw_lf + g_ph*lf
        drm = jR - (Q/mr)*Cw_xyl - k_on*rm + k_off*rb + 0.1*(QPh/mr)*CPh - mur*rm
        drb = k_on*rm - k_off*rb - mur*rb
        dst = (Q/ms)*(Cw_xyl - Cw_st) - mus*st
        dlf = f3*(Q/ml)*Cw_st - 1.1*(QPh/ml)*CPh - mul*lf
        dgr = f4*(Q/mg)*Cw_st + (QPh/mg)*CPh - mug*gr
        return [drm, drb, dst, dlf, dgr]

    sol = solve_ivp(rhs, (0.0, SEASON), np.zeros(5), method="BDF", rtol=1e-5, atol=1e-8,
                    t_eval=[SEASON])
    rm, rb, st, lf, gr = sol.y[:, -1]
    mlf, mst = float(fMl(SEASON)), float(fMs(SEASON))
    straw = (st*mst + lf*mlf)/(mst+mlf)
    return {"root": rm+rb, "straw": straw, "grain": gr, "rm": rm, "rb": rb}


def _fit(name, k_off):
    """fit g_xy->straw and g_ph->grain (brentq, monotone); f_xy from oryza refit."""
    c = CONG[name]; o = obs[name]
    f_xy = c.get("f_xy_oryza") or c["f_xy_recommended"]
    def gx(lg):
        return np.log10(max(simulate2(name, f_xy, 10**lg, 0.0, k_off)["straw"], 1e-6)) - np.log10(o["straw"])
    a, b = -6, 1
    g_xy = 10**brentq(gx, a, b, xtol=1e-2, maxiter=40) if gx(a)*gx(b) < 0 else (10**a if abs(gx(a))<abs(gx(b)) else 10**b)
    def gp(lg):
        return np.log10(max(simulate2(name, f_xy, g_xy, 10**lg, k_off)["grain"], 1e-6)) - np.log10(o["grain"])
    g_ph = 10**brentq(gp, a, b, xtol=1e-2, maxiter=40) if gp(a)*gp(b) < 0 else (10**a if abs(gp(a))<abs(gp(b)) else 10**b)
    return f_xy, g_xy, g_ph


def main():
    k_off = 0.02
    print(f"2-POOL root prototype (k_off={k_off}/d slow bound store; lipid loads from mobile pool)\n")
    print(f"{'PFAS':7}{'nC':>3} | {'root p/o':>16}{'straw p/o':>16}{'grain p/o':>16}  rm/rb")
    errs = {"root": [], "shoot": []}
    for name in ("PFOA", "PFDA", "PFUnDA", "PFDoDA"):
        f_xy, g_xy, g_ph = _fit(name, k_off)
        r = simulate2(name, f_xy, g_xy, g_ph, k_off); o = obs[name]
        for k in ("root",):
            errs["root"].append((np.log10(max(r[k],1e-6))-np.log10(o[k]))**2)
        for k in ("straw", "grain"):
            errs["shoot"].append((np.log10(max(r[k],1e-6))-np.log10(o[k]))**2)
        print(f"{name:7}{CONG[name]['n_C']:>3} | "
              f"{r['root']:6.1f}/{o.get('root',float('nan')):<7.1f}"
              f"{r['straw']:6.1f}/{o.get('straw',float('nan')):<7.1f}"
              f"{r['grain']:6.1f}/{o.get('grain',float('nan')):<7.1f}  {r['rm']:.2f}/{r['rb']:.1f}")
    import math
    print(f"\nlong-incl root RMSE {math.sqrt(np.mean(errs['root'])):.3f}  shoot RMSE "
          f"{math.sqrt(np.mean(errs['shoot'])):.3f}")
    print("LC4 question: can the 2-pool match long-chain ROOT and SHOOT simultaneously?")


if __name__ == "__main__":
    main()
