"""Numerical choices on the DEFAULT configuration: do the results converge? (descriptive)

The tool's rule is that no value enters without a source. A purely numerical value (a mass
floor that keeps the ODE off zero, a solver tolerance, a sub-daily grid) has no literature
source; it is acceptable only if the results do not depend on it -- i.e. moving it further
toward its limit (floor -> 0, tolerance -> 0, grid -> fine, N -> many) changes no reported
number. The yardstick is the precision the results are reported at: 3 decimals of log10
(|change| < 0.0005) for scores and for log10 BAFs. A value whose results keep moving as it
approaches its limit is not a numerical choice but a structural one, and has to be fixed in
the model, not tuned.

Varied one at a time on the default configuration:
  * the organ-mass floor `organ_mass.MASS_FLOOR` (1e-9 kg; it also decides when the grain
    formation gate opens: the day the grain leaves the floor);
  * the integrator: exponential steps per driver interval (`SUBSTEPS` 16, more where the rates
    change fast, `linear_expm.GRADE` 0.2; the ODE is affine in the state), against more steps,
    a finer grading and BDF at rtol 1e-8;
  * the harvest operator's sub-daily grid (`harvest.REFINE_PER_DAY` 24);
  * the number of stem segments N (4);
  * the rate of the linear equilibrium root entry `model_api.LAM_EQ` (1e6 per day per kg;
    1e4 moved the BAFs by 0.0008 log10 and was raised under this rule, PREREG section 6).
Reported: the Yamazaki two-harvest observable (log10 RMSE), the largest |delta log10 BAF| over
the 11 congeners x 3 tissues against the shipped setting, the PFOA young-panicle peak, and the
Li 2025 stage-resolved shoot / grain TF scores (the time-course-sensitive endpoint).
The daily driver grid is not varied: it is ORYZA2000's own integration step.

    python validation/numerics_convergence.py        (~10 min with 4 workers)
"""
import functools
import multiprocessing as mp
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "src"))
sys.path.insert(0, HERE)

PRECISION = 0.0005          # half the last reported digit of a 3-decimal log10 number
SHIPPED = dict(floor=1e-9, solver="expm", substeps=16, grade=0.2, rtol=1e-6, atol=1e-9,
               refine=24, N=4, lam_eq=1.0e6)
VARIANTS = [("shipped", {})]
for f in (1e-8, 1e-10, 1e-12):
    VARIANTS.append((f"floor {f:g} kg", dict(floor=f)))
VARIANTS.append(("substeps 32", dict(substeps=32)))
VARIANTS.append(("substeps 64", dict(substeps=64)))
VARIANTS.append(("grade 0.05", dict(grade=0.05)))
VARIANTS.append(("BDF rtol 1e-8, atol 1e-11", dict(solver="bdf", rtol=1e-8, atol=1e-11)))
VARIANTS.append(("REFINE_PER_DAY 96", dict(refine=96)))
VARIANTS.append(("N 8", dict(N=8)))
VARIANTS.append(("N 16", dict(N=16)))
VARIANTS.append(("lam_eq 1e5", dict(lam_eq=1.0e5)))
VARIANTS.append(("lam_eq 1e7", dict(lam_eq=1.0e7)))


def _patch(s):
    import harvest
    import model_api as api
    import organ_mass as om
    import pfas_rice_plant_module_4pool_surf as P4
    import pfas_rice_plant_module_nstem_leaf as NL
    om.MASS_FLOOR = s["floor"]                      # model_api's inner floors read it at call time
    P4.MASS_FLOOR = NL.MASS_FLOOR = s["floor"]      # the RHS floors (imported by value)
    orig_matrix = om.__dict__.setdefault("_orig_organ_matrix", om.organ_matrix)   # patch once
    om.organ_matrix = functools.wraps(orig_matrix)(
        lambda b, **kw: orig_matrix(b, **{**kw, "floor": s["floor"]}))
    orig_formed = P4.__dict__.setdefault("_orig_grain_formed", P4.grain_formed)
    formed = lambda gm, floor=s["floor"]: orig_formed(gm, floor=floor)
    P4.grain_formed = NL.grain_formed = formed
    api._SITE_RUN = None
    harvest.REFINE_PER_DAY = s["refine"]
    P4.SOLVER, P4.SUBSTEPS = s["solver"], s["substeps"]
    import linear_expm as LE
    LE.GRADE = s["grade"]
    for mod in (NL, P4):
        orig_ivp = mod.__dict__.setdefault("_orig_solve_ivp", mod.solve_ivp)
        mod.solve_ivp = (lambda o: lambda *a, **k: o(*a, **{**k, "rtol": s["rtol"],
                                                            "atol": s["atol"]}))(orig_ivp)


def _yamazaki(api, kw):
    obs = api.yamazaki_obs()
    geo, res = {}, {k: [] for k in ("root", "straw", "grain")}
    for c in api.CONGENERS:
        if c not in obs:
            continue
        geo[c] = api.yamazaki_harvest_baf(c, **kw)
        for k in res:
            if obs[c]["use"][k]:
                res[k].append(np.log10(geo[c][k] / obs[c]["baf_dw"][k]))
    rm = lambda e: float(np.sqrt(np.mean(np.square(e))))
    yz = {k: rm(v) for k, v in res.items()}
    yz["all"] = rm(np.concatenate([res[k] for k in res]))
    return geo, yz


def run(v):
    name, over = v
    s = {**SHIPPED, **over}
    _patch(s)
    import model_api as api
    import li2025_stage_resolved as L
    kw = dict(N=s["N"], lam_eq=s["lam_eq"])
    r = api.simulate("PFOA", **kw)
    cg = np.asarray(r["conc"]["grain"], float)
    i = int(np.argmax(cg))
    geo, yz = _yamazaki(api, kw)
    days = L.stage_days()
    model = {c: L.model_tf(api.simulate(c, Cwo=1.0, **kw), days) for c in L.CONGENERS}
    li = L.score(L.load_data(), model, np.random.default_rng(0))
    return dict(name=name, peak=float(cg[i] / cg[-1]), peak_day=float(r["t"][i]),
                yamazaki=yz, baf=geo, li_shoot=li["rmse"]["shoot"], li_grain=li["rmse"]["grain"])


def main():
    with mp.get_context("fork").Pool(4, maxtasksperchild=1) as pool:   # fresh, unpatched child per variant
        out = pool.map(run, VARIANTS, chunksize=1)
    b = next(r for r in out if r["name"] == "shipped")
    print("Numerical choices on the default configuration (descriptive)")
    print(f"converged = every reported number moves by < {PRECISION} (3-decimal log10)\n")
    print(f"{'variant':24s} {'Yamazaki all/root/straw/grain':>30s} "
          f"{'max|dlogBAF|':>12s} {'PFOA peak':>13s} {'Li shoot':>9s} {'Li grain':>9s}  verdict")
    for r in out:
        dmax = max(abs(np.log10(r["baf"][c][k] / b["baf"][c][k])) for c in b["baf"] for k in b["baf"][c])
        dscore = max([abs(r["yamazaki"][k] - b["yamazaki"][k]) for k in r["yamazaki"]]
                     + [abs(r["li_shoot"] - b["li_shoot"]), abs(r["li_grain"] - b["li_grain"])])
        verdict = "" if r["name"] == "shipped" else (
            "same" if max(dmax, dscore) < PRECISION else "MOVES")
        y = r["yamazaki"]
        print(f"{r['name']:24s} {y['all']:.3f}/{y['root']:.3f}/{y['straw']:.3f}/{y['grain']:.3f}"
              f"      {dmax:12.4f} {r['peak']:6.2f}x d{r['peak_day']:<4.0f} {r['li_shoot']:9.3f} "
              f"{r['li_grain']:9.3f}  {verdict}")


if __name__ == "__main__":
    main()
