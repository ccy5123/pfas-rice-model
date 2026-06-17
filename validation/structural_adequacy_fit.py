#!/usr/bin/env python3
# =============================================================================
# validation/structural_adequacy_fit.py
# -----------------------------------------------------------------------------
# STRUCTURAL ADEQUACY via a CONSTRAINED (DOF>0) fit.
#
# The question is NOT prediction but: "can the model STRUCTURE reproduce the
# measured experiments when calibrated?" A fit only answers this with degrees of
# freedom > 0 -- the W2 fit (3 params/congener vs 3 obs/congener; 33 params / 33
# obs globally) is SATURATED (DOF 0), so its RMSE 0.029 is guaranteed and says
# nothing about the structure.
#
# CONSTRAINED fit: per-congener f_xy (the physical TSCF, the one lever expected to
# vary by chain) but a SINGLE shared L_Ph and a SINGLE shared kappa_d across ALL
# congeners. 11 + 2 = 13 params for 11 congeners x 3 tissues = 33 obs -> DOF = 20.
#
# Fit by a fast STAGED procedure that exploits the model's near-separability
# (root<-kappa_d+binding, straw<-f_xy, grain<-L_Ph), so it is ~500 ODE solves
# instead of the tens of thousands a joint least_squares needs. The result is an
# ACHIEVABLE RMSE at DOF 20 (an upper bound on the joint optimum) -- enough to
# decide structural adequacy. Reported with the per-tissue breakdown + DOF so the
# claim is honest.
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


CONG = [c for c in PAR["congeners"] if c["name"] in obs and c.get("f_xy_W2fit") is not None]
TISS = ("root", "straw", "grain")
_cache = {}


def _predict(c, f_xy, L_Ph, kappa_d):
    """Final-time tissue BAFs. Fast solve: only the terminal value is needed, and
    the structural-adequacy question is robust to a looser tol (rtol 1e-4)."""
    key = (c["name"], round(f_xy, 5), round(L_Ph, 5), round(kappa_d, 4))
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


def _le(a, b):  # squared log10 error, guarded
    return (np.log10(max(a, 1e-6)) - np.log10(max(b, 1e-6))) ** 2


def _best_fxy(c, L_Ph, kappa_d, lo=1e-3, hi=1.0):
    """per-congener f_xy matching the observed straw (straw is monotone in f_xy);
    brentq on log10(straw_pred/straw_obs), clamped when the target is unreachable."""
    o = obs[c["name"]]
    if "straw" not in o:
        return c["f_xy_recommended"]
    def g(lx):
        return np.log10(max(_predict(c, 10.0 ** lx, L_Ph, kappa_d)["straw"], 1e-6)) - np.log10(o["straw"])
    a, b = np.log10(lo), np.log10(hi)
    ga, gb = g(a), g(b)
    if ga * gb < 0:
        return 10.0 ** brentq(g, a, b, xtol=1e-2, rtol=1e-3, maxiter=40)
    return lo if abs(ga) < abs(gb) else hi          # target outside reach -> nearest bound


def _rmse(fxy_by, L_Ph, kappa_d):
    errs = []
    for c in CONG:
        pr = _predict(c, fxy_by[c["name"]], L_Ph, kappa_d); o = obs[c["name"]]
        errs += [_le(pr[k], o[k]) for k in TISS if k in o]
    return float(np.sqrt(np.mean(errs))), errs


def main():
    kd_grid = np.logspace(-2, np.log10(50), 9)
    lph_grid = np.logspace(-4, 0, 9)

    # Stage 1: global kappa_d from ROOT (root ~ uptake/dilution; weak in f_xy/L_Ph).
    L_Ph = 0.1
    best = (1e9, kd_grid[0])
    for kd in kd_grid:
        e = sum(_le(_predict(c, c["f_xy_recommended"], L_Ph, kd)["root"], obs[c["name"]]["root"])
                for c in CONG if "root" in obs[c["name"]])
        best = min(best, (e, kd))
    kappa_d = best[1]

    # Stage 2: per-congener f_xy from STRAW (kappa_d fixed) via brentq.
    fxy_by = {c["name"]: _best_fxy(c, L_Ph, kappa_d) for c in CONG}

    # Stage 3: global L_Ph from GRAIN (kappa_d, f_xy fixed).
    best = (1e9, lph_grid[0])
    for lph in lph_grid:
        e = sum(_le(_predict(c, fxy_by[c["name"]], lph, kappa_d)["grain"], obs[c["name"]]["grain"])
                for c in CONG if "grain" in obs[c["name"]])
        best = min(best, (e, lph))
    L_Ph = best[1]

    # one refinement pass of f_xy with the updated L_Ph
    fxy_by = {c["name"]: _best_fxy(c, L_Ph, kappa_d) for c in CONG}

    rmse, errs = _rmse(fxy_by, L_Ph, kappa_d)
    n_obs = len(errs); n_par = len(CONG) + 2; dof = n_obs - n_par
    # per-tissue RMSE
    per = {}
    for k in TISS:
        ek = [_le(_predict(c, fxy_by[c["name"]], L_Ph, kappa_d)[k], obs[c["name"]][k])
              for c in CONG if k in obs[c["name"]]]
        per[k] = float(np.sqrt(np.mean(ek)))

    print(f"CONSTRAINED global fit: per-congener f_xy ({len(CONG)}) + global L_Ph + global kappa_d "
          f"= {n_par} params for {n_obs} obs  ->  DOF = {dof}")
    print(f"global L_Ph = {L_Ph:.4g}   global kappa_d = {kappa_d:.4g}\n")
    print(f"{'PFAS':8}{'nC':>3}{'f_xy':>8} | {'root p/o':>16}{'straw p/o':>16}{'grain p/o':>16}")
    for c in CONG:
        pr = _predict(c, fxy_by[c["name"]], L_Ph, kappa_d); o = obs[c["name"]]
        print(f"{c['name']:8}{c['n_C']:>3}{fxy_by[c['name']]:>8.4f} | "
              f"{pr['root']:>7.2f}/{o.get('root', float('nan')):<7.2f}"
              f"{pr['straw']:>7.2f}/{o.get('straw', float('nan')):<7.2f}"
              f"{pr['grain']:>7.2f}/{o.get('grain', float('nan')):<7.2f}")
    print(f"\nper-tissue log10 RMSE: root {per['root']:.3f}  straw {per['straw']:.3f}  "
          f"grain {per['grain']:.3f}")
    print(f"CONSTRAINED-fit log10 RMSE = {rmse:.3f}   (DOF {dof}; vs saturated W2 0.029 at DOF 0, "
          f"a-priori monotone ~0.84)")
    return rmse, dof, per


if __name__ == "__main__":
    main()
