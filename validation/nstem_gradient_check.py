#!/usr/bin/env python3
"""
N-segment stem validation with MEASURED drivers (GAP-B / task 2)
================================================================

Drives the multi-height stem (`pfas_rice_plant_module_nstem`) with the measured
crop-physiology forcings -- transpiration Q_TP(t) from `forcing_rice` (Kumari
2022 Kc/ET0 + Nay Htoon 2018 T/E split) and per-organ biomass M_s(t) from
`growth_rice` (ORYZA IR72 partitioning, Bouman 2006 / Li 2017) -- and compares
the predicted vertical stem gradient to Yamazaki 2023 (S18/S19).

Honest findings (this is a direction/structure test, not a passing validation):
  1. Mass conservation: exact (sole source = M_root*j_R).
  2. With the MEASURED Q_TP and ORYZA biomass, a monotone f_xy gives an UPWARD
     stem gradient for the short/mid-chain PFCAs (correct direction), but the
     gradient magnitude is ~uniform across chain length and the long-chain
     flat/down (near-MQL) reversal is NOT reproduced.
  3. Reason: the equilibrium (well-mixed) segment only discriminates by binding
     when the growth-retention sink mu*M_s is comparable to the advective
     throughput Q_s/B; at realistic (ORYZA) biomass that crossover B* ~
     Q_s/(M_s*mu_s) ~ 100+, far above the congener range (B 1-60), so all
     congeners read "up".  (An earlier apparent flip used inflated placeholder
     biomass.)  Reproducing the chain-length-resolved gradient needs a KINETIC
     radial xylem<->tissue exchange (a sorbing column), not instantaneous
     equilibrium -- the open model-structure item.

Run: python validation/nstem_gradient_check.py
"""
import json, os, sys, csv
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT_DIR, "src"))

from pfas_rice_plant_module_nstem import NStemModel, PlantInputsN, make_stem_compartments
from pfas_rice_plant_module_4pool_surf import (
    Environment, Compound, binding_factors, root_uptake)
import forcing_rice as fr
import growth_rice as gr

PAR = json.load(open(os.path.join(ROOT_DIR, "params", "parameters.json")))
comp = PAR["tissue_composition_recommended"]
carr = PAR["carrier_MichaelisMenten"]

OBS = {}
with open(os.path.join(ROOT_DIR, "data_obs", "yamazaki_stem_height.csv")) as f:
    for r in csv.DictReader(filter(lambda l: not l.startswith("#"), f)):
        OBS[r["pfas"]] = [float(r[k]) for k in
                          ("stem_0_20", "stem_20_40", "stem_40_60", "stem_gt60")]

N, SEASON = 4, 120.0
t = np.linspace(0.0, SEASON, 481)
Qtp = fr.Q_TP(t, SEASON)                 # MEASURED transpiration (forcing_rice)
M = gr.M_for_nstem(t, N, SEASON)         # ORYZA IR72 biomass (growth_rice)
inputs = PlantInputsN(t=t, Cwo=np.full_like(t, 1.0), Qtp=Qtp, M=M)
_kw = lambda d: dict(theta=d["theta_fw"], f_prot=d["f_prot"], f_PL=d["f_PL"], f_cw=d["f_cw"])
comps = make_stem_compartments(N, _kw(comp["stem"]), _kw(comp["root"]), _kw(comp["grain_brown"]))
tau = np.array([0.30, 0.28, 0.24, 0.18]); tau = tau / tau.sum() * 0.85


def build(nm, f_xy):
    c = next(x for x in PAR["congeners"] if x["name"] == nm)
    cmpd = Compound(name=nm, K_prot=c["K_prot_Lkg"], K_PL=c["K_PL_Lkg"],
                    K_cw=c["K_cw_wholecw_Lkg"]["stem"], kappa_d=2.0,
                    Vmax_in=carr["Vmax_in"], Km_in=carr["Km_in"],
                    Vmax_out=carr["Vmax_out"], Km_out=carr["Km_out"], L_Ph=0.01, f_xy=f_xy)
    return NStemModel(env=Environment(), cmpd=cmpd, comps=comps, inputs=inputs, tau=tau), c


def main():
    print(f"MEASURED drivers: peak Q_TP={Qtp.max():.3f} L/d/hill, "
          f"final straw={M[-1,1]*N:.4f} kg/hill, grain={M[-1,-1]:.4f} kg/hill")
    # mass conservation (PFOA)
    pfoa = next(x for x in PAR["congeners"] if x["name"] == "PFOA")
    m, _ = build("PFOA", pfoa["f_xy_recommended"])
    sol = m.solve(t); B = binding_factors(comps, m.cmpd)
    ti = 70.0; C = sol.sol(ti); dC = m.rhs(ti, C); Mv = inputs.M_(ti); dMv = inputs.dM_(ti)
    jR = root_uptake(inputs.Cwo_(ti), (C / B)[0], m.cmpd, m.env)
    dmass = float(np.sum(dMv * C + Mv * dC)); src = float(Mv[0] * jR)
    print(f"[1] mass conservation @t=70: {dmass:.3e} vs {src:.3e} -> "
          f"{'OK' if abs(dmass - src) < 1e-6 * abs(src) + 1e-9 else 'LEAK'}")
    print(f"\n[2] stem top/bottom gradient (monotone f_xy):")
    print(f"    {'PFAS':7}{'B_stem':>8} | {'pred':>6} {'obs':>6}  dir")
    hit = tot = 0
    for cc in PAR["congeners"]:
        nm = cc["name"]
        if nm not in OBS:
            continue
        mm, _ = build(nm, cc["f_xy_recommended"])
        seg = mm.solve(t).y[1:1 + N, -1]
        pr = seg[-1] / seg[0] if seg[0] > 0 else float("inf")
        orr = OBS[nm][-1] / OBS[nm][0]
        fp = "UP" if pr > 1.3 else ("flat" if pr > 0.77 else "DOWN")
        fo = "UP" if orr > 1.3 else ("flat" if orr > 0.77 else "DOWN")
        tot += 1; hit += (fp == fo)
        print(f"    {nm:7}{cc['B_k_basisA_Lkg_fw']['stem']:8.1f} | {pr:6.2f} {orr:6.2f}  {fp}/{fo}")
    print(f"    direction match: {hit}/{tot}.  Short/mid PFCAs UP (correct); long-chain "
          f"(near-MQL) flat/down NOT captured -> equilibrium segment needs kinetic radial sorption.")


if __name__ == "__main__":
    main()
