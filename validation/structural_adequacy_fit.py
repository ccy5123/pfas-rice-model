#!/usr/bin/env python3
# =============================================================================
# validation/structural_adequacy_fit.py
# -----------------------------------------------------------------------------
# STRUCTURAL ADEQUACY via CONSTRAINED (DOF>0) fits, driven by the MECHANISTIC
# ORYZA2000 biomass (src/oryza_growth.py) + the measured transpiration
# (src/forcing_rice.py) -- NOT the illustrative logistic placeholder.
#
# Question (the user's reframe): not out-of-sample prediction, but "can the model
# STRUCTURE reproduce the measured experiments via fitting?" A fit only answers
# this with degrees of freedom > 0 (the W2 fit is 33 params / 33 obs, DOF 0 --
# saturated, so its RMSE 0.029 is guaranteed and structure-blind).
#
# Three constrained scenarios decide WHICH per-congener lever the structure needs,
# always keeping DOF > 0 (11 congeners x 3 tissues = 33 obs):
#   A  per-congener f_xy + global L_Ph + global kappa_d        (13 params, DOF 20)  [shoot]
#   B  per-congener f_xy + per-congener L_Ph + global kappa_d  (23 params, DOF 10)  [grain]
#   C  per-congener f_xy + global L_Ph + per-congener kappa_d  (23 params, DOF 10)  [root]
#
# Each is fit by a fast staged procedure exploiting the model's near-separability
# (root<-kappa_d, straw<-f_xy, grain<-L_Ph) with 1-D brentq solves on the
# monotone tissue responses, and a loose-tol terminal-only ODE solve. Reports the
# per-tissue log10 RMSE + DOF per scenario so the adequacy claim is honest.
#
#   python validation/structural_adequacy_fit.py
# =============================================================================
import json, csv, os, sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
from pfas_rice_plant_module_4pool_surf import (  # noqa: E402
    Environment, Compound, Compartment, RiceUptakeModel, PlantInputs,
    ROOT as R, STEM, LEAF, FRUIT)
import oryza_growth as og        # noqa: E402  (the user's mechanistic ORYZA2000)
import forcing_rice as fr        # noqa: E402  (measured FAO-56 transpiration)

PAR = json.load(open(os.path.join(ROOT, "params", "parameters.json")))
obs = {}
with open(os.path.join(ROOT, "data_obs", "obs_baf_Yamazaki.csv"), newline="") as f:
    for r in csv.DictReader(f):
        obs.setdefault(r["compound"], {})[r["tissue"]] = float(r["baf"])

comp = PAR["tissue_composition_recommended"]
carr = PAR["carrier_MichaelisMenten"]
env = Environment()

# ---- biomass driver: the MECHANISTIC ORYZA2000 (oryza_growth) + measured Q_TP --
SEASON = 120.0
t = np.linspace(0.0, SEASON, 121)
_b = og.organ_biomass_oryza(t, p=og.OryzaParams(season=SEASON))
M = np.maximum(np.column_stack([_b["root"], _b["stem"], _b["leaf"], _b["grain"]]), 1e-4)
inputs = PlantInputs(t=t, Cwo=np.full_like(t, 1.0), Qtp=fr.Q_TP(t, SEASON), M=M,
                     leaf_loss=_b.get("leaf_death_rate"))
Mf = inputs.M_(t[-1])


def _compartments():
    return [Compartment("root", comp["root"]["theta_fw"], comp["root"]["f_prot"], comp["root"]["f_PL"], comp["root"]["f_cw"]),
            Compartment("stem", comp["stem"]["theta_fw"], comp["stem"]["f_prot"], comp["stem"]["f_PL"], comp["stem"]["f_cw"]),
            Compartment("leaf", comp["leaf"]["theta_fw"], comp["leaf"]["f_prot"], comp["leaf"]["f_PL"], comp["leaf"]["f_cw"], S=20.0),
            Compartment("grain", comp["grain_brown"]["theta_fw"], comp["grain_brown"]["f_prot"], comp["grain_brown"]["f_PL"], comp["grain_brown"]["f_cw"], S=2.0)]


CONG = [c for c in PAR["congeners"] if c["name"] in obs and c.get("f_xy_W2fit") is not None]
TISS = ("root", "straw", "grain")
_cache = {}


def _predict(c, f_xy, L_Ph, kappa_d):
    key = (c["name"], round(f_xy, 5), round(L_Ph, 6), round(kappa_d, 4))
    if key in _cache:
        return _cache[key]
    cmpd = Compound(name=c["name"], K_prot=c["K_prot_Lkg"], K_PL=c["K_PL_Lkg"],
                    K_cw=c["K_cw_wholecw_Lkg"]["root"], kappa_d=kappa_d,
                    Vmax_in=carr["Vmax_in"], Km_in=carr["Km_in"],
                    Vmax_out=carr["Vmax_out"], Km_out=carr["Km_out"], L_Ph=L_Ph, f_xy=f_xy)
    model = RiceUptakeModel(env=env, cmpd=cmpd, comps=_compartments(), inputs=inputs)
    sol = solve_ivp(model.rhs, (float(t[0]), float(t[-1])), np.zeros(4), method="BDF",
                    rtol=1e-4, atol=1e-7, t_eval=[float(t[-1])])
    C = sol.y[:, -1]
    straw = (C[STEM] * Mf[STEM] + C[LEAF] * Mf[LEAF]) / (Mf[STEM] + Mf[LEAF])
    out = {"root": C[R], "straw": straw, "grain": C[FRUIT]}
    _cache[key] = out
    return out


def _le(a, b):
    return (np.log10(max(a, 1e-6)) - np.log10(max(b, 1e-6))) ** 2


def _fit_scalar(c, tissue, var, lo, hi, f_xy, L_Ph, kappa_d):
    """1-D brentq for the scalar `var` ('f_xy'|'L_Ph'|'kappa_d') matching the
    observed `tissue` (monotone response), clamped when the target is unreachable."""
    o = obs[c["name"]]
    if tissue not in o:
        return {"f_xy": f_xy, "L_Ph": L_Ph, "kappa_d": kappa_d}[var]
    def pred(x):
        kw = {"f_xy": f_xy, "L_Ph": L_Ph, "kappa_d": kappa_d}; kw[var] = 10.0 ** x
        return _predict(c, kw["f_xy"], kw["L_Ph"], kw["kappa_d"])[tissue]
    a, b = np.log10(lo), np.log10(hi)
    ga = np.log10(max(pred(a), 1e-6)) - np.log10(o[tissue])
    gb = np.log10(max(pred(b), 1e-6)) - np.log10(o[tissue])
    if ga * gb < 0:
        return 10.0 ** brentq(lambda x: np.log10(max(pred(x), 1e-6)) - np.log10(o[tissue]),
                              a, b, xtol=1e-2, rtol=1e-3, maxiter=40)
    return lo if abs(ga) < abs(gb) else hi


def _global_grid(tissue, grid, fxy_by, L_Ph, kappa_d, var):
    best = (1e9, grid[0])
    for v in grid:
        e = 0.0
        for c in CONG:
            if tissue in obs[c["name"]]:
                kw = {"f_xy": fxy_by[c["name"]], "L_Ph": L_Ph, "kappa_d": kappa_d}; kw[var] = v
                e += _le(_predict(c, kw["f_xy"], kw["L_Ph"], kw["kappa_d"])[tissue], obs[c["name"]][tissue])
        best = min(best, (e, v))
    return best[1]


def _rmse(params_by, glob):
    per = {k: [] for k in TISS}
    for c in CONG:
        f_xy = params_by[c["name"]]["f_xy"]
        L_Ph = params_by[c["name"]].get("L_Ph", glob.get("L_Ph"))
        kappa_d = params_by[c["name"]].get("kappa_d", glob.get("kappa_d"))
        pr = _predict(c, f_xy, L_Ph, kappa_d); o = obs[c["name"]]
        for k in TISS:
            if k in o:
                per[k].append(_le(pr[k], o[k]))
    allp = sum(per.values(), [])
    return (float(np.sqrt(np.mean(allp))),
            {k: float(np.sqrt(np.mean(v))) for k, v in per.items()})


KD_GRID = np.logspace(-2, np.log10(50), 9)
LPH_GRID = np.logspace(-4, 0, 9)
FXY = (1e-3, 1.0); LPH = (1e-5, 1.0); KD = (1e-2, 50.0)


def scenario(name, per_lph, per_kd):
    """Fit with f_xy always per-congener; L_Ph and kappa_d per-congener or global."""
    fxy_by = {c["name"]: c["f_xy_recommended"] for c in CONG}
    glob = {}
    # kappa_d (root)
    if per_kd:
        kd_by = {c["name"]: _fit_scalar(c, "root", "kappa_d", *KD, fxy_by[c["name"]], 0.01, 2.0) for c in CONG}
    else:
        kd_g = _global_grid("root", KD_GRID, fxy_by, 0.01, 2.0, "kappa_d"); glob["kappa_d"] = kd_g
        kd_by = {c["name"]: kd_g for c in CONG}
    # f_xy (straw) given kappa_d, a starting L_Ph
    L0 = 0.003
    for c in CONG:
        fxy_by[c["name"]] = _fit_scalar(c, "straw", "f_xy", *FXY, fxy_by[c["name"]], L0, kd_by[c["name"]])
    # L_Ph (grain)
    if per_lph:
        lph_by = {c["name"]: _fit_scalar(c, "grain", "L_Ph", *LPH, fxy_by[c["name"]], L0, kd_by[c["name"]]) for c in CONG}
    else:
        lph_g = _global_grid("grain", LPH_GRID, fxy_by, L0, kd_by[CONG[0]["name"]], "L_Ph"); glob["L_Ph"] = lph_g
        lph_by = {c["name"]: lph_g for c in CONG}
    # one refine of f_xy with final L_Ph
    for c in CONG:
        fxy_by[c["name"]] = _fit_scalar(c, "straw", "f_xy", *FXY, fxy_by[c["name"]], lph_by[c["name"]], kd_by[c["name"]])

    params_by = {c["name"]: {"f_xy": fxy_by[c["name"]], "L_Ph": lph_by[c["name"]], "kappa_d": kd_by[c["name"]]}
                 for c in CONG}
    rmse, per = _rmse(params_by, glob)
    n_obs = sum(len([k for k in TISS if k in obs[c["name"]]]) for c in CONG)
    n_par = len(CONG) + (len(CONG) if per_lph else 1) + (len(CONG) if per_kd else 1)
    dof = n_obs - n_par
    gtxt = ", ".join(f"{k}={v:.4g}" for k, v in glob.items()) or "(none)"
    print(f"[{name}]  params={n_par}/{n_obs}  DOF={dof}  global: {gtxt}")
    print(f"        per-tissue RMSE: root {per['root']:.3f}  straw {per['straw']:.3f}  grain {per['grain']:.3f}"
          f"   overall {rmse:.3f}")
    return name, dof, per, rmse


def main():
    print("biomass driver: MECHANISTIC ORYZA2000 (oryza_growth) + measured Q_TP (forcing_rice)\n")
    rows = [scenario("A shoot   (f_xy; global L_Ph, global kappa_d)", per_lph=False, per_kd=False),
            scenario("B +grain  (f_xy, per-cong L_Ph; global kappa_d)", per_lph=True, per_kd=False),
            scenario("C +root   (f_xy, per-cong kappa_d; global L_Ph)", per_lph=False, per_kd=True)]
    print("\nSUMMARY (vs saturated W2 0.029 at DOF 0; a-priori monotone ~0.84):")
    for name, dof, per, rmse in rows:
        print(f"  DOF {dof:>2}  overall {rmse:.3f}  | root {per['root']:.3f} straw {per['straw']:.3f} "
              f"grain {per['grain']:.3f}  | {name}")
    return rows


if __name__ == "__main__":
    main()
