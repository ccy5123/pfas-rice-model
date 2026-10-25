"""The grain formation gate's numerical set-up on the DEFAULT configuration (descriptive;
decides nothing).

Two numbers with no source decide when the panicle starts taking PFAS:
  * `organ_mass.organ_matrix(floor=1e-4)` -- every organ's mass is floored at 1e-4 kg so the
    ODE never divides by zero; before panicle initiation the grain sits on this floor;
  * the gate ramps 0 -> 1 while the grain mass goes from the floor to 1.5 x the floor
    (`PlantInputsNL.__post_init__`, the same rule in the 4-pool core).
On the site-year the panicle leaves the floor on day 70 and passes 1.5 x the floor a day
later, so the gate is in effect a step at panicle appearance -- at a mass the floor sets.
In "gate" mode the formed panicle takes the whole residual xylem stream at once, so the
young panicle's concentration depends on how small it is when that happens.

This script moves the two numbers (ramp 1.1 / 3 / 10 x the floor; floor 1e-5 / 1e-3 kg) with
everything else as shipped and reports what moves:
  * the young-panicle peak (PFOA max grain concentration / harvest value, and its day);
  * the Yamazaki two-harvest observable (log10 RMSE, dry basis, the switch's scoring);
  * the five transfer endpoints of the switch / panicle amendment (Tang, Kim point and
    interval, Li shoot and grain).
The floor also floors the other organs, so the floor rows are not a pure gate change.

    python validation/gate_numerics_sensitivity.py        (~10 min, one worker per variant)
"""
import functools
import multiprocessing as mp
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "src"))
sys.path.insert(0, HERE)

# name -> (organ-mass floor [kg], gate fully open at this multiple of the floor)
VARIANTS = {"shipped (floor 1e-4 kg, ramp 1.5x)": (1e-4, 1.5),
            "ramp 1.1x": (1e-4, 1.1), "ramp 3x": (1e-4, 3.0), "ramp 10x": (1e-4, 10.0),
            "floor 1e-5 kg": (1e-5, 1.5), "floor 1e-3 kg": (1e-3, 1.5)}


def _patch(floor, ramp):
    """Set the floor and the ramp in this (worker) process only."""
    from scipy.interpolate import interp1d
    import model_api as api
    import organ_mass as om
    import pfas_rice_plant_module_nstem_leaf as NL
    orig_matrix = om.__dict__.setdefault("_orig_organ_matrix", om.organ_matrix)   # patch once
    om.organ_matrix = functools.wraps(orig_matrix)(lambda b: orig_matrix(b, floor=floor))
    api._SITE_RUN = None                     # drop the cached site runs built on the old floor
    orig_post = NL.__dict__.setdefault("_orig_post_init", NL.PlantInputsNL.__post_init__)

    def post(self):
        orig_post(self)
        gm = self.M[:, -1]
        glo, ghi = float(gm.min()), float(gm.max())
        gate = (np.clip((gm - glo) / ((ramp - 1.0) * glo + 1e-30), 0.0, 1.0)
                if (ghi > 0.0 and glo < 0.05 * ghi) else np.ones(len(self.t)))
        self._grain_gate = interp1d(self.t, gate, kind="linear", bounds_error=False,
                                    fill_value="extrapolate")
    NL.PlantInputsNL.__post_init__ = post


def run(name):
    _patch(*VARIANTS[name])
    import model_api as api
    import default_switch as DS
    import panicle_share_amendment as PA
    d = api.scenario_drivers()
    g, t = np.asarray(d["M"])[:, 3], np.asarray(d["t"])
    r = api.simulate("PFOA")
    cg = np.asarray(r["conc"]["grain"], float)
    i = int(np.argmax(cg))
    geo = {c: api.yamazaki_harvest_baf(c) for c in DS.YCONG}
    res = DS.residuals(geo)
    yz = {k: DS.rmse(res, (k,)) for k in DS.TISS}
    yz["all"] = DS.rmse(res)
    transfer, _ = PA.transfer_scores("gate")
    return dict(name=name, open_day=float(t[np.argmax(g > g.min())]),
                peak=float(cg[i] / cg[-1]), peak_day=float(r["t"][i]), yamazaki=yz,
                grain={c: geo[c]["grain"] for c in DS.YCONG}, transfer=transfer)


def main():
    with mp.get_context("fork").Pool(len(VARIANTS), maxtasksperchild=1) as pool:   # one unpatched child per variant
        out = pool.map(run, list(VARIANTS), chunksize=1)
    base = out[0]
    short = ("Tang", "Kim pt", "Kim int", "Li shoot", "Li grain")
    print("Grain formation gate numerics (default configuration; descriptive)")
    print(f"{'variant':36s} {'opens':>5s} {'PFOA peak':>13s}  {'Yamazaki all/root/straw/grain':>29s}"
          f"  {'med|dlog grain|':>15s}  " + " ".join(f"{s:>8s}" for s in short))
    for r in out:
        dg = float(np.median([abs(np.log10(r["grain"][c] / base["grain"][c])) for c in base["grain"]]))
        y = r["yamazaki"]
        print(f"{r['name']:36s} d{r['open_day']:<4.0f} {r['peak']:5.2f}x d{r['peak_day']:<4.0f}"
              f"  {y['all']:.3f}/{y['root']:.3f}/{y['straw']:.3f}/{y['grain']:.3f}"
              f"            {dg:15.3f}  " + " ".join(f"{v:8.3f}" for v in r["transfer"].values()))
    peaks = [r["peak"] for r in out]
    dy = max(abs(r["yamazaki"]["all"] - base["yamazaki"]["all"]) for r in out)
    dt = max(abs(r["transfer"][k] - base["transfer"][k]) for r in out for k in r["transfer"])
    print(f"\nPFOA young-panicle peak {min(peaks):.2f}-{max(peaks):.2f}x its harvest value across the "
          f"variants; largest change of the Yamazaki score {dy:.3f}, of a transfer score {dt:.3f}.")


if __name__ == "__main__":
    main()
