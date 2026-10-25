"""Unmeasured inputs on the DEFAULT configuration (descriptive; decides nothing).

`docs/ASSUMPTIONS.md` reports, for every unmeasured input, how far the Yamazaki BAFs move
when the input goes to the end of its range. Those numbers were produced by
`validation/refit_site_year.py` on the 4-pool free-anion fits. Since the pre-registered switch
(`validation/default_switch.py`) the default is the redistributed shoot + lipid loading, so the
register needs the same table on that configuration. This script produces it.

Method (the refit_site_year "Unknown inputs" protocol, applied to the default):
  * the Yamazaki two-harvest observable (two runs, 09-01 and 09-10, geometric mean, dry basis,
    each congener at its irrigation-water Cwo), 11 congeners;
  * one input at a time, to each end of its range; fitted parameters held fixed;
  * effect = median over congeners of the larger |delta log10 BAF| of the two ends;
  * yardstick = the model difference between the two harvest runs (as in the register);
    an effect above it is marked consequential.
Nothing is re-fitted here (the switch decision concerned transfer, not a re-fit), so unlike
the register's 4-pool rows there is no re-fit column.

The driver inputs (nursery, density, plants per hill, heading, harvest convention, CO2) and the
two observation conventions (dead blades, single mean-date run) are the register's; the model
constants are the switch's pre-registered SENS ranges (E_m, Km, phi, T_C_Ph, and the three shoot
architecture levers). The grain xylem surface (`s_grain_dry`) does not exist in the
redistributed shoot (its grain share is `lam_grain`), so it is not run.

    python validation/default_sensitivity.py        (~4 min with 4 workers)
"""
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "src"))
sys.path.insert(0, HERE)
import model_api as api            # noqa: E402
import site_scenario as ss         # noqa: E402

TISS = ("root", "straw", "grain")
OBS = api.yamazaki_obs()
YCONG = [c for c in api.CONGENERS if c in OBS]
DRIVER_KEYS = ("sbdur", "tmpsb", "nh", "nplh", "days_heading_to_harvest", "harvest_is_maturity",
               "co2")

# name -> list of settings (each a dict; driver keys go to scenario_drivers, the rest to simulate)
INPUTS = {
    "nursery (SBDUR 14 / 30 d)": [dict(sbdur=14), dict(sbdur=30)],
    "covered seedbed (TMPSB 9 C)": [dict(tmpsb=9.0)],
    "density (15 / 22 hills m-2)": [dict(nh=15.0), dict(nh=22.0)],
    "plants per hill (5)": [dict(nplh=5.0)],
    "heading (30 / 31 d before harvest)": [dict(days_heading_to_harvest=30.0),
                                           dict(days_heading_to_harvest=31.0)],
    "harvest not at maturity": [dict(harvest_is_maturity=False)],
    "CO2 (400 / 430 ppm)": [dict(co2=400.0), dict(co2=430.0)],
    "dead blades lost": [dict(dead_leaf_attached=0.0)],
    "single mean-date run": ["mean_run"],
    "E_m (-90 mV)": [dict(E_m_mV=-90.0)],
    "carrier Km x100": [dict(km_scale=100.0)],
    "phi (0.05 / 0.2)": [dict(phi=0.05), dict(phi=0.2)],
    "T_C,Ph (4 / 12 L/kg)": [dict(T_C_Ph=4.0), dict(T_C_Ph=12.0)],
    "stem_transp_frac (0.30 / 0.60)": [dict(stem_transp_frac=0.30), dict(stem_transp_frac=0.60)],
    "retention (0.40 / 0.80)": [dict(retention=0.40), dict(retention=0.80)],
    "lam_grain (0.025 / 0.18)": [dict(lam_grain=0.025), dict(lam_grain=0.18)],
}


def _observable(args):
    """Two-run (or single mean-date) geometric-mean dry-basis BAF for one congener."""
    c, setting = args
    cwo = OBS[c]["water_ng_L"] / 1000.0
    if setting == "mean_run":
        harvests, drv_kw, sim_kw = ["mean"], {}, {}
    else:
        setting = dict(setting or {})
        drv_kw = {k: setting.pop(k) for k in list(setting) if k in DRIVER_KEYS}
        harvests, sim_kw = list(ss.YAMAZAKI_HARVESTS), setting
    runs = []
    for h in harvests:
        d = api.scenario_drivers(Cwo=cwo, harvest=h, **drv_kw)
        r = api.simulate(c, drivers=d, **sim_kw)       # the default configuration
        runs.append({k: r["baf_final_dw"][k] for k in TISS})
    geo = {k: float(10 ** np.mean([np.log10(r[k]) for r in runs])) for k in TISS}
    return c, geo, runs


def main(workers=None):
    print("=" * 88)
    print(f"UNMEASURED INPUTS ON THE DEFAULT (shoot={api.DEFAULT_SHOOT}, "
          f"lipid_loading={api.DEFAULT_LIPID_LOADING}, uptake={api.DEFAULT_UPTAKE})")
    print("=" * 88)
    with ProcessPoolExecutor(max_workers=workers or os.cpu_count() or 2) as pool:
        base = {c: (g, r) for c, g, r in pool.map(_observable, [(c, None) for c in YCONG])}
        ref = {k: float(np.median([abs(np.log10(base[c][1][0][k] / base[c][1][1][k]))
                                   for c in YCONG])) for k in TISS}
        print(f"yardstick (two harvest runs): root {ref['root']:.3f}  straw {ref['straw']:.3f}  "
              f"grain {ref['grain']:.3f}\n")
        print(f"  {'input':36s}{'root':>8s}{'straw':>8s}{'grain':>8s}   consequential")
        eff = {}
        for name, settings in INPUTS.items():
            per = []
            for st in settings:
                out = {c: g for c, g, _ in pool.map(_observable, [(c, st) for c in YCONG])}
                per.append({c: {k: float(np.log10(out[c][k] / base[c][0][k])) for k in TISS}
                            for c in YCONG})
            med = {k: float(np.median([max(abs(p[c][k]) for p in per) for c in YCONG]))
                   for k in TISS}
            cons = [k for k in TISS if med[k] > ref[k]]
            eff[name] = dict(median_abs_delta=med, consequential=cons)
            print(f"  {name:36s}{med['root']:8.3f}{med['straw']:8.3f}{med['grain']:8.3f}   "
                  f"{', '.join(cons) or 'none'}")
    return dict(reference=ref, effects=eff)


if __name__ == "__main__":
    main()
