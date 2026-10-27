#!/usr/bin/env python3
"""
The paddy scenario's own pore-water shapes, per congener (PREREG 9.11).

`cwo_profile="flooded"` and the soil-inventory mode take the time course of the root-zone pore water
from this table instead of an analytic surrogate: each column is c(t)/c(0) of the Tsukuba Andosol paddy
(`soil_hydrus.run_paddy_hydrus`, the contaminated-soil / clean-water case, source "soil") on the
Yamazaki 2020 site-year (its transpiration, evaporation, rain, heading and season). The well-mixed
first-order surrogate it replaces could not take the leaching front this scenario produces (9 of 13
congeners ended on its numerical ceiling; PREREG 9.11). A structure without a column interpolates in
log10 Kd (`model_api.hydrus_shape`).

Run (needs the built HYDRUS-1D engine + phydrus):
    python validation/cwo_shape_table.py      (~1 min)
"""
from __future__ import annotations

import os
import sys
import warnings

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
warnings.filterwarnings("ignore")
import model_api as api      # noqa: E402
import soil_hydrus as sh     # noqa: E402

OUT = os.path.join(ROOT, "params", "cwo_shapes_hydrus.csv")


def build():
    cols, t_ref = {}, None
    for cong in api.CONGENERS:
        c = api._CONG[cong]
        Kd = sh.paddy_kd(c["n_C"], c["group"])
        _, res = sh.inputs_from_hydrus(c["n_C"], c["group"])
        t = np.asarray(res.t, float)
        t_ref = t if t_ref is None else t_ref
        cols[cong] = (Kd, np.interp(t_ref, t, res.Cw / res.Cw[0]))
        print(f"  {cong:8} Kd {Kd:9.3f}  c/c0 at end {cols[cong][1][-1]:.3f}  "
              f"(mean {np.mean(cols[cong][1]):.3f})")
    names = sorted(cols, key=lambda k: cols[k][0])           # ascending Kd, for interpolation
    with open(OUT, "w") as f:
        f.write("# The Tsukuba Andosol paddy's root-zone pore water c(t)/c(0) per congener (PREREG 9.11),\n"
                "# built by validation/cwo_shape_table.py from soil_hydrus.run_paddy_hydrus (source 'soil':\n"
                "# contaminated soil, clean water) on the Yamazaki 2020 site-year. Row 'Kd' = the column's\n"
                "# linear Kd [L/kg] (Koc QSPR x site f_oc); t_day = days after transplanting.\n")
        f.write("t_day," + ",".join(names) + "\n")
        f.write("Kd," + ",".join(f"{cols[n][0]:.6g}" for n in names) + "\n")
        for i, tt in enumerate(t_ref):
            f.write(f"{tt:g}," + ",".join(f"{cols[n][1][i]:.6g}" for n in names) + "\n")
    print(f"wrote {os.path.relpath(OUT, ROOT)}  ({len(names)} congeners x {len(t_ref)} times)")


if __name__ == "__main__":
    if not sh.hydrus_available():
        sys.exit("HYDRUS-1D engine/phydrus not available")
    build()
