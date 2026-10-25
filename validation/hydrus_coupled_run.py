#!/usr/bin/env python3
"""
HYDRUS-1D -> plant coupling: full Method-A soil->plant run
==========================================================

Drives the plant model (the default configuration: redistributed shoot + lipid loading)
with a REAL HYDRUS-1D pore-water trajectory C_w^o(t) (built in ``src/soil_hydrus.py`` from
the compiled HYDRUS engine via phydrus) and contrasts it with the constant-Cwo baseline.

Both runs share the SAME crop drivers -- the calibration site-year (src/site_scenario.py;
the HYDRUS soil run is driven by its ORYZA2000 transpiration) -- and the HYDRUS pore water
is normalised to the same season-mean exposure (1.0) as the baseline -- so the only
difference is the realistic soil temporal structure: weakly-sorbed short chains leach
during flooding (pore water collapses, so the late-filling grain sees little),
strongly-sorbed long chains stay buffered (BAF ~ unchanged). Values are dry-basis harvest
BAF (src/harvest.py) of the default configuration (`model_api.simulate` defaults).

Run:  python validation/hydrus_coupled_run.py
Saves: validation/figures/hydrus_coupled.png, validation/hydrus_coupled.csv
"""
import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))

import model_api as api  # noqa: E402
import soil_hydrus as sh  # noqa: E402

CONGENERS = ["PFBA", "PFHxA", "PFOA", "PFNA", "PFOS", "PFDoDA"]


def solve_baf(name, drivers=None):
    """Final dry-basis root/straw/grain BAF (relative to season-mean pore water = 1.0)."""
    r = api.simulate(name, drivers=drivers)
    ref = float(np.mean(r["Cwo"]))
    return {k: r["conc_dw"][k][-1] / ref for k in ("root", "straw", "grain")}


def main():
    if not sh.hydrus_available():
        print("HYDRUS executable not built (see src/soil_hydrus.py). Aborting.")
        return 1

    rows, traj = [], {}
    print(f"{'cong':8}{'Kd':>9}  {'root b/h':>16}{'straw b/h':>16}{'grain b/h':>16}")
    for name in CONGENERS:
        c = api._CONG[name]
        drv, res = api.hydrus_drivers(name, Cwo_ref=1.0)
        traj[name] = (drv["t"], drv["Cwo"], res.Kd)
        b = solve_baf(name)
        h = solve_baf(name, drv)
        rows.append(dict(congener=name, n_C=c["n_C"], group=c["group"], Kd=res.Kd,
                         **{f"{k}_baseline": b[k] for k in b},
                         **{f"{k}_hydrus": h[k] for k in h}))
        print(f"{name:8}{res.Kd:>9.3f}  "
              f"{b['root']:>7.2f}/{h['root']:<7.2f}"
              f"{b['straw']:>7.2f}/{h['straw']:<7.2f}"
              f"{b['grain']:>7.2f}/{h['grain']:<7.2f}")

    out_csv = os.path.join(HERE, "hydrus_coupled.csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {out_csv}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.3))
        for name in CONGENERS:
            t, cwo, Kd = traj[name]
            ax1.plot(t, cwo, lw=2, label=f"{name} (Kd={Kd:.2g})")
        ax1.set_xlabel("time [day]"); ax1.set_ylabel("pore-water C_w^o(t)  [norm., mean=1]")
        ax1.set_title("HYDRUS-1D pore water (clean-water flooding → drainage)")
        ax1.legend(fontsize=7); ax1.axvline(90, ls=":", c="grey")   # flood_until default

        x = np.arange(len(CONGENERS)); w = 0.35
        gb = [r["grain_baseline"] for r in rows]
        gh = [r["grain_hydrus"] for r in rows]
        ax2.bar(x - w / 2, gb, w, label="grain BAF (constant Cwo)")
        ax2.bar(x + w / 2, gh, w, label="grain BAF (HYDRUS Cwo)")
        ax2.set_xticks(x); ax2.set_xticklabels(CONGENERS, rotation=45, ha="right")
        ax2.set_ylabel("grain BAF [L/kg dw]"); ax2.set_title("late-filling grain: effect of leaching")
        ax2.legend(fontsize=8)
        fig.tight_layout()
        figpath = os.path.join(HERE, "figures", "hydrus_coupled.png")
        os.makedirs(os.path.dirname(figpath), exist_ok=True)
        fig.savefig(figpath, dpi=130)
        print(f"wrote {figpath}")
    except Exception as e:
        print(f"(plot skipped: {e})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
