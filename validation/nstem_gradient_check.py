#!/usr/bin/env python3
"""
N-segment stem validation (GAP-B "open modeling item")
=======================================================

Tests whether the multi-height stem module (`pfas_rice_plant_module_nstem`)
reproduces the *vertical* tissue gradient Yamazaki 2023 measured (Table S18/S19,
stem 0-20/20-40/40-60/>60 cm) using a single MONOTONE f_xy -- the open item the
single mass-weighted "straw" compartment could not handle.

Findings reproduced here:
  1. Mass conservation: for gamma=0 the sole source is M_root*j_R (exact).
  2. Chain dependence: with a serial stem (transpiration draw-off + radial
     exchange + growth dilution) the gradient direction FLIPS with chain length
     -- short chains concentrate UPWARD (transpiration), long chains stay
     flat/down -- exactly as observed, AND with a monotone f_xy.
  3. The flip location is set by B* ~ Q_s/(M_s*mu_s): it requires a REALISTIC
     transpiration/biomass ratio.  With placeholder drivers (peak Q_TP ~0.4 L/d)
     B* is far above the data range so everything reads "up"; lowering Q to a
     plausible per-plant scale moves the flip to ~C11, matching the data.
  4. Open: PFOS (high B but observed UP) is not captured by a binding-driven
     monotone f_xy -> a PFSA-specific transport term is needed (task 3).

Run: python validation/nstem_gradient_check.py
Needs measured Q_TP(t)/M_s(t) for an absolute fit (task 2); this is a structural
direction test only.
"""
import json, os, sys, csv
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT_DIR, "src"))

from pfas_rice_plant_module_nstem import NStemModel, PlantInputsN, make_stem_compartments
from pfas_rice_plant_module_4pool_surf import (
    Environment, Compound, binding_factors, root_uptake, _logistic)

PAR = json.load(open(os.path.join(ROOT_DIR, "params", "parameters.json")))
comp = PAR["tissue_composition_recommended"]
carr = PAR["carrier_MichaelisMenten"]

# observed per-height geomean (Indica+Japonica), pg/g dw
OBS = {}
with open(os.path.join(ROOT_DIR, "data_obs", "yamazaki_stem_height.csv")) as f:
    for r in csv.DictReader(filter(lambda l: not l.startswith("#"), f)):
        OBS[r["pfas"]] = [float(r[k]) for k in
                          ("stem_0_20", "stem_20_40", "stem_40_60", "stem_gt60")]

N = 4
t = np.linspace(0.0, 120.0, 481)
straw = _logistic(t, 1e-3, 0.090, 0.11, 27.0)
M = np.column_stack([_logistic(t, 1e-3, 0.030, 0.10, 20.0)] + [straw / N] * N
                    + [_logistic(t, 1e-5, 0.025, 0.18, 80.0)])
root_kw = dict(theta=comp["root"]["theta_fw"], f_prot=comp["root"]["f_prot"],
               f_PL=comp["root"]["f_PL"], f_cw=comp["root"]["f_cw"])
stem_kw = dict(theta=comp["stem"]["theta_fw"], f_prot=comp["stem"]["f_prot"],
               f_PL=comp["stem"]["f_PL"], f_cw=comp["stem"]["f_cw"])
grain_kw = dict(theta=comp["grain_brown"]["theta_fw"], f_prot=comp["grain_brown"]["f_prot"],
                f_PL=comp["grain_brown"]["f_PL"], f_cw=comp["grain_brown"]["f_cw"])
comps = make_stem_compartments(N, stem_kw, root_kw, grain_kw)
tau = np.array([0.30, 0.28, 0.24, 0.18]); tau = tau / tau.sum() * 0.85   # residual to grain


def build(nm, f_xy, Qtp):
    c = next(x for x in PAR["congeners"] if x["name"] == nm)
    cmpd = Compound(name=nm, K_prot=c["K_prot_Lkg"], K_PL=c["K_PL_Lkg"],
                    K_cw=c["K_cw_wholecw_Lkg"]["stem"], kappa_d=2.0,
                    Vmax_in=carr["Vmax_in"], Km_in=carr["Km_in"],
                    Vmax_out=carr["Vmax_out"], Km_out=carr["Km_out"],
                    L_Ph=0.01, f_xy=f_xy)
    inp = PlantInputsN(t=t, Cwo=np.full_like(t, 1.0), Qtp=Qtp, M=M)
    return NStemModel(env=Environment(), cmpd=cmpd, comps=comps, inputs=inp, tau=tau), c


def main():
    # 1) mass conservation (PFOA)
    Qtp = (0.05 + 0.35 * np.exp(-((t - 75.0) ** 2) / (2 * 25.0 ** 2))) / 0.4 * 0.08
    pfoa_fxy = next(x for x in PAR["congeners"] if x["name"] == "PFOA")["f_xy_recommended"]
    m, _ = build("PFOA", pfoa_fxy, Qtp)
    sol = m.solve(t); B = binding_factors(comps, m.cmpd)
    ti = 70.0; C = sol.sol(ti); dC = m.rhs(ti, C); Mv = m.inputs.M_(ti); dMv = m.inputs.dM_(ti)
    jR = root_uptake(m.inputs.Cwo_(ti), (C / B)[0], m.cmpd, m.env)
    dmass = float(np.sum(dMv * C + Mv * dC)); src = float(Mv[0] * jR)
    print(f"[1] mass conservation @t=70: d/dt sum(M C)={dmass:.3e} vs M_root*jR={src:.3e}"
          f"  -> {'OK' if abs(dmass - src) < 1e-6 * abs(src) + 1e-9 else 'LEAK'}")

    # 2) chain-dependent gradient at a realistic transpiration scale
    print(f"\n[2] stem top/bottom gradient (monotone f_xy, peak Q_TP={Qtp.max():.3f} L/d):")
    print(f"    {'PFAS':7}{'B_stem':>8}{'f_xy':>9} | {'pred':>6} {'obs':>6}  match")
    hit = tot = 0
    for cc in PAR["congeners"]:
        nm = cc["name"]
        if nm not in OBS:
            continue
        mm, _ = build(nm, cc["f_xy_recommended"], Qtp)
        seg = mm.solve(t).y[1:1 + N, -1]
        pr = seg[-1] / seg[0] if seg[0] > 0 else float("inf")
        orr = OBS[nm][-1] / OBS[nm][0]
        fp = "UP" if pr > 1.3 else ("flat" if pr > 0.77 else "DOWN")
        fo = "UP" if orr > 1.3 else ("flat" if orr > 0.77 else "DOWN")
        tot += 1; hit += (fp == fo)
        print(f"    {nm:7}{cc['B_k_basisA_Lkg_fw']['stem']:8.1f}{cc['f_xy_recommended']:9.4f} | "
              f"{pr:6.2f} {orr:6.2f}  {fp}/{fo}{'  OK' if fp == fo else ''}")
    print(f"    direction match: {hit}/{tot}  (PFCA carboxylates match; PFOS is the PFSA outlier -> task 3)")


if __name__ == "__main__":
    main()
