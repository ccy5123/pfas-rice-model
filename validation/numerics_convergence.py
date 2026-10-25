"""Numerical choices on the DEFAULT configuration: do the results converge? (descriptive)

The tool's rule is that no value enters without a source. A purely numerical value (a mass
floor that keeps the ODE off zero, a solver tolerance, a sub-daily grid) has no literature
source; it is acceptable only if the results do not depend on it -- i.e. moving it further
toward its limit (floor -> 0, tolerance -> 0, grid -> fine, N -> many) changes no reported
number. The yardstick is the precision the results are reported at: 3 decimals of log10
(|change| < 0.0005) for scores and for log10 BAFs. A value whose results keep moving as it
approaches its limit is not a numerical choice but a structural one, and has to be fixed in
the model, not tuned.

Varied one at a time, both panicle routings ("gate" = the default, "size"):
  * the organ-mass floor in `organ_mass.organ_matrix` (1e-4 kg as shipped);
  * the grain formation gate's ramp (fully open at 1.5 x the floor);
  * the BDF tolerances (rtol 1e-6, atol 1e-9);
  * the harvest operator's sub-daily grid (`harvest.REFINE_PER_DAY` 24);
  * the number of stem segments N (4).
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
SHIPPED = dict(floor=1e-4, ramp=1.5, rtol=1e-6, atol=1e-9, refine=24, N=4)
VARIANTS = []
for share in ("gate", "size"):
    VARIANTS.append((share, "shipped", {}))
    for f in (1e-3, 1e-5, 1e-6, 1e-7):
        VARIANTS.append((share, f"floor {f:g} kg", dict(floor=f)))
    VARIANTS.append((share, "rtol 1e-8, atol 1e-11", dict(rtol=1e-8, atol=1e-11)))
    VARIANTS.append((share, "REFINE_PER_DAY 96", dict(refine=96)))
    VARIANTS.append((share, "N 8", dict(N=8)))
    VARIANTS.append((share, "N 16", dict(N=16)))
for r in (1.01, 1.1, 3.0):
    VARIANTS.append(("gate", f"ramp {r:g}x", dict(ramp=r)))


def _patch(s):
    from scipy.interpolate import interp1d
    import harvest
    import model_api as api
    import organ_mass as om
    import pfas_rice_plant_module_nstem_leaf as NL
    orig_matrix = om.__dict__.setdefault("_orig_organ_matrix", om.organ_matrix)   # patch once
    om.organ_matrix = functools.wraps(orig_matrix)(lambda b: orig_matrix(b, floor=s["floor"]))
    api._SITE_RUN = None
    harvest.REFINE_PER_DAY = s["refine"]
    orig_ivp = NL.__dict__.setdefault("_orig_solve_ivp", NL.solve_ivp)
    NL.solve_ivp = lambda *a, **k: orig_ivp(*a, **{**k, "rtol": s["rtol"], "atol": s["atol"]})
    orig_post = NL.__dict__.setdefault("_orig_post_init", NL.PlantInputsNL.__post_init__)

    def post(self):
        orig_post(self)
        gm = self.M[:, -1]
        glo, ghi = float(gm.min()), float(gm.max())
        gate = (np.clip((gm - glo) / ((s["ramp"] - 1.0) * glo + 1e-30), 0.0, 1.0)
                if (ghi > 0.0 and glo < 0.05 * ghi) else np.ones(len(self.t)))
        self._grain_gate = interp1d(self.t, gate, kind="linear", bounds_error=False,
                                    fill_value="extrapolate")
    NL.PlantInputsNL.__post_init__ = post


def run(v):
    share, name, over = v
    s = {**SHIPPED, **over}
    _patch(s)
    import model_api as api
    import default_switch as DS
    import li2025_stage_resolved as L
    kw = dict(grain_share=share, N=s["N"])
    r = api.simulate("PFOA", **kw)
    cg = np.asarray(r["conc"]["grain"], float)
    i = int(np.argmax(cg))
    geo = {c: api.yamazaki_harvest_baf(c, **kw) for c in DS.YCONG}
    res = DS.residuals(geo)
    yz = {k: DS.rmse(res, (k,)) for k in DS.TISS}
    yz["all"] = DS.rmse(res)
    days = L.stage_days()
    model = {c: L.model_tf(api.simulate_nstem_leaf(c, Cwo=1.0, lipid_loading=True,
                                                   grain_share=share, N=s["N"]), days)
             for c in L.CONGENERS}
    li = L.score(L.load_data(), model, np.random.default_rng(0))
    return dict(share=share, name=name, peak=float(cg[i] / cg[-1]), peak_day=float(r["t"][i]),
                yamazaki=yz, baf={c: geo[c] for c in DS.YCONG},
                li_shoot=li["rmse"]["shoot"], li_grain=li["rmse"]["grain"])


def main():
    with mp.get_context("fork").Pool(4, maxtasksperchild=1) as pool:   # fresh, unpatched child per variant
        out = pool.map(run, VARIANTS, chunksize=1)
    base = {r["share"]: r for r in out if r["name"] == "shipped"}
    print("Numerical choices on the default configuration (descriptive)")
    print(f"converged = every reported number moves by < {PRECISION} (3-decimal log10)\n")
    print(f"{'share':5s} {'variant':24s} {'Yamazaki all/root/straw/grain':>30s} "
          f"{'max|dlogBAF|':>12s} {'PFOA peak':>13s} {'Li shoot':>9s} {'Li grain':>9s}  verdict")
    for r in out:
        b = base[r["share"]]
        dmax = max(abs(np.log10(r["baf"][c][k] / b["baf"][c][k])) for c in b["baf"] for k in b["baf"][c])
        dscore = max([abs(r["yamazaki"][k] - b["yamazaki"][k]) for k in r["yamazaki"]]
                     + [abs(r["li_shoot"] - b["li_shoot"]), abs(r["li_grain"] - b["li_grain"])])
        verdict = "" if r["name"] == "shipped" else (
            "same" if max(dmax, dscore) < PRECISION else "MOVES")
        y = r["yamazaki"]
        print(f"{r['share']:5s} {r['name']:24s} {y['all']:.3f}/{y['root']:.3f}/{y['straw']:.3f}/{y['grain']:.3f}"
              f"      {dmax:12.4f} {r['peak']:6.2f}x d{r['peak_day']:<4.0f} {r['li_shoot']:9.3f} "
              f"{r['li_grain']:9.3f}  {verdict}")


if __name__ == "__main__":
    main()
