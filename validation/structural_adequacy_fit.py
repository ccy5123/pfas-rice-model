#!/usr/bin/env python3
# =============================================================================
# validation/structural_adequacy_fit.py
# -----------------------------------------------------------------------------
# STRUCTURAL ADEQUACY via a CONSTRAINED (DOF>0) fit.
#
# The question is NOT prediction but: "can the model STRUCTURE reproduce the
# measured experiments when calibrated?" A fit only answers this if it has
# degrees of freedom > 0 -- the W2 fit (3 params/congener vs 3 obs/congener;
# 33 params / 33 obs globally) is SATURATED (DOF 0), so its RMSE 0.029 is
# guaranteed and says nothing about the structure.
#
# This runs a CONSTRAINED global fit: per-congener f_xy (the physical TSCF, the
# one lever expected to vary by chain) but a SINGLE shared L_Ph and a SINGLE
# shared kappa_d across ALL congeners. That is 11 + 2 = 13 parameters for
# 11 congeners x 3 tissues = 33 observations  ->  DOF = 20.
#
# If the structure reproduces the data at DOF 20, structural adequacy is
# demonstrated (not just expressiveness). The result is reported with the RMSE
# AND the degrees of freedom so the claim is honest.
#
#   python validation/structural_adequacy_fit.py
# =============================================================================
import json, csv, os, sys
import numpy as np
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
from pfas_rice_plant_module_4pool_surf import (  # noqa: E402
    Environment, Compound, Compartment, RiceUptakeModel, PlantInputs,
    _logistic, ROOT as R, STEM, LEAF, FRUIT)

PAR = json.load(open(os.path.join(ROOT, "params", "parameters.json")))
obs = {}
with open(os.path.join(ROOT, "data_obs", "obs_baf_Yamazaki.csv"), newline="") as f:
    for r in csv.DictReader(f):
        obs.setdefault(r["compound"], {})[r["tissue"]] = float(r["baf"])

comp = PAR["tissue_composition_recommended"]
carr = PAR["carrier_MichaelisMenten"]
env = Environment()
t = np.linspace(0.0, 120.0, 121)
inputs = PlantInputs(
    t=t, Cwo=np.full_like(t, 1.0),
    Qtp=0.05 + 0.35 * np.exp(-((t - 75.0) ** 2) / (2 * 25.0 ** 2)),
    M=np.column_stack([_logistic(t, 1e-3, 0.030, 0.10, 20.0), _logistic(t, 1e-3, 0.040, 0.10, 25.0),
                       _logistic(t, 1e-3, 0.050, 0.12, 30.0), _logistic(t, 1e-5, 0.025, 0.18, 80.0)]))
Mf = inputs.M_(t[-1])


def _compartments():
    return [Compartment("root", comp["root"]["theta_fw"], comp["root"]["f_prot"], comp["root"]["f_PL"], comp["root"]["f_cw"]),
            Compartment("stem", comp["stem"]["theta_fw"], comp["stem"]["f_prot"], comp["stem"]["f_PL"], comp["stem"]["f_cw"]),
            Compartment("leaf", comp["leaf"]["theta_fw"], comp["leaf"]["f_prot"], comp["leaf"]["f_PL"], comp["leaf"]["f_cw"], S=20.0),
            Compartment("grain", comp["grain_brown"]["theta_fw"], comp["grain_brown"]["f_prot"], comp["grain_brown"]["f_PL"], comp["grain_brown"]["f_cw"], S=2.0)]


# congeners with obs + a (W2) transport entry, as in reproduce_demo
CONG = [c for c in PAR["congeners"] if c["name"] in obs and c.get("f_xy_W2fit") is not None]
NAMES = [c["name"] for c in CONG]
TISS = ("root", "straw", "grain")


def _predict(c, f_xy, L_Ph, kappa_d):
    cmpd = Compound(name=c["name"], K_prot=c["K_prot_Lkg"], K_PL=c["K_PL_Lkg"],
                    K_cw=c["K_cw_wholecw_Lkg"]["root"], kappa_d=kappa_d,
                    Vmax_in=carr["Vmax_in"], Km_in=carr["Km_in"],
                    Vmax_out=carr["Vmax_out"], Km_out=carr["Km_out"], L_Ph=L_Ph, f_xy=f_xy)
    C = RiceUptakeModel(env=env, cmpd=cmpd, comps=_compartments(), inputs=inputs).solve(t).y[:, -1]
    straw = (C[STEM] * Mf[STEM] + C[LEAF] * Mf[LEAF]) / (Mf[STEM] + Mf[LEAF])
    return {"root": C[R], "straw": straw, "grain": C[FRUIT]}


def _unpack(x):
    f_xy = 10.0 ** x[:len(CONG)]
    L_Ph = 10.0 ** x[-2]
    kappa_d = 10.0 ** x[-1]
    return f_xy, L_Ph, kappa_d


def _resid(x):
    f_xy, L_Ph, kappa_d = _unpack(x)
    r = []
    for c, fx in zip(CONG, f_xy):
        try:
            pr = _predict(c, fx, L_Ph, kappa_d)
        except Exception:
            r += [10.0] * 3; continue
        o = obs[c["name"]]
        for k in TISS:
            if k in o:
                r.append(np.log10(max(pr[k], 1e-6)) - np.log10(o[k]))
    return np.asarray(r)


def main():
    n = len(CONG)
    # warm start: recommended monotone f_xy; global L_Ph/kappa_d mid-range
    x0 = np.concatenate([np.log10([max(c["f_xy_recommended"], 1e-3) for c in CONG]),
                         [np.log10(0.1), np.log10(2.0)]])
    lo = np.concatenate([np.full(n, np.log10(1e-3)), [np.log10(1e-4), np.log10(1e-2)]])
    hi = np.concatenate([np.full(n, np.log10(1.0)), [np.log10(1.0), np.log10(50.0)]])
    res = least_squares(_resid, x0, bounds=(lo, hi), method="trf",
                        diff_step=1e-2, xtol=1e-9, ftol=1e-9, max_nfev=400)
    f_xy, L_Ph, kappa_d = _unpack(res.x)
    n_obs = len(res.fun); n_par = n + 2; dof = n_obs - n_par
    rmse = float(np.sqrt(np.mean(res.fun ** 2)))

    print(f"CONSTRAINED global fit: per-congener f_xy ({n}) + global L_Ph + global kappa_d "
          f"= {n_par} params for {n_obs} obs  ->  DOF = {dof}")
    print(f"global L_Ph = {L_Ph:.4g}   global kappa_d = {kappa_d:.4g}\n")
    print(f"{'PFAS':8}{'nC':>3}{'f_xy':>8} | {'root p/o':>16}{'straw p/o':>16}{'grain p/o':>16}")
    for c, fx in zip(CONG, f_xy):
        pr = _predict(c, fx, L_Ph, kappa_d); o = obs[c["name"]]
        print(f"{c['name']:8}{c['n_C']:>3}{fx:>8.4f} | "
              f"{pr['root']:>7.2f}/{o.get('root', float('nan')):<7.2f}"
              f"{pr['straw']:>7.2f}/{o.get('straw', float('nan')):<7.2f}"
              f"{pr['grain']:>7.2f}/{o.get('grain', float('nan')):<7.2f}")
    print(f"\nCONSTRAINED-fit log10 RMSE = {rmse:.3f}   (DOF {dof}; vs saturated W2 0.029 at DOF 0, "
          f"a-priori monotone ~0.84)")
    return rmse, dof


if __name__ == "__main__":
    main()
