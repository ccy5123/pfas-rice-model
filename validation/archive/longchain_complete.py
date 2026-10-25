#!/usr/bin/env python3
# =============================================================================
# validation/longchain_complete.py
# -----------------------------------------------------------------------------
# Does the COMPLETE long-chain resolution actually close the long chains?
#
# FINDINGS sec.7 proposed the "complete long-chain resolution = 2-pool (free +
# lipid-bound) root + lipid-facilitated loading + enhanced long-chain active
# carrier." LC4/LC5b/LC6 each tested ONE lever in isolation. This combines all
# three into ONE model and asks the decisive question the proposal never tested:
# does the complete recipe reproduce the long-chain ROOT and SHOOT (straw)
# SIMULTANEOUSLY across C10-C12?
#
# Recipe, per congener (reusing the LC prototype, validation/twopool_longchain.py):
#   1. find the active-carrier multiplier (Vmax_in) that reproduces the measured
#      ROOT (the LC6 root-matching multiplier), then
#   2. WITH that carrier, fit the lipid-facilitated loading g_xy (-> straw) and
#      g_ph (-> grain) on the 2-pool root.
# Report root/straw/grain pred/obs vs Yamazaki and the long-chain (nC>=10) RMSE.
#
# RESULT (honest, refining sec.7): the carrier that fixes the ROOT structurally
# OVER-feeds the SHOOT -- the enhanced mobile pool's free loading f_xy*Cw_m alone
# exceeds the observed straw, and lipid g_xy>=0 cannot subtract -- so straw is
# over-predicted for C11-C12 (PFUnDA ~3.3x, PFDoDA ~2.3x). Root and grain close;
# the shoot does NOT close simultaneously. The "complete resolution" is therefore
# NOT a simultaneous closure: the long chains need a root->shoot DECOUPLING
# (irreversible root sequestration that does not translocate), consistent with the
# long-chain "strongly retained in root" literature (LC1).
#
#   python validation/longchain_complete.py
# =============================================================================
import os
import sys
import math

import numpy as np
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "src"))
sys.path.insert(0, HERE)
import twopool_longchain as T  # noqa: E402  (the LC 2-pool prototype + helpers)

BASE = T.carr["Vmax_in"]
K_OFF = 0.02
TISS = ("root", "straw", "grain")
SERIES = ("PFOA", "PFNA", "PFDA", "PFUnDA", "PFDoDA")


def _root_matching_mult(name):
    """The LC6 active-carrier multiplier that reproduces the measured ROOT
    (g_xy=g_ph=0, isolating uptake -> root)."""
    c = T.CONG[name]
    o = T.obs[name]
    f_xy = c.get("f_xy_oryza") or c["f_xy_recommended"]

    def g(lm_):
        r = T.simulate2(name, f_xy, 0.0, 0.0, K_OFF, vmax_in=BASE * 10 ** lm_)["root"]
        return np.log10(max(r, 1e-6)) - np.log10(o["root"])

    a, b = 0.0, 4.0
    if g(a) * g(b) < 0:
        return 10 ** brentq(g, a, b, xtol=1e-2, maxiter=40)
    return 10 ** a if abs(g(a)) < abs(g(b)) else 10 ** b


def complete(name):
    """The complete recipe for one congener: root-matching carrier, then fit lipid."""
    mult = _root_matching_mult(name)
    vm = BASE * mult
    f_xy, g_xy, g_ph = T._fit(name, K_OFF, vmax_in=vm)
    r = T.simulate2(name, f_xy, g_xy, g_ph, K_OFF, vmax_in=vm)
    return mult, r


def run():
    rows = []
    err_lc = {t: [] for t in TISS}
    for name in SERIES:
        c = T.CONG[name]
        o = T.obs[name]
        nC = c["n_C"]
        mult, r = complete(name)
        rows.append((name, nC, mult, r, o))
        if nC >= 10:
            for t in TISS:
                if t in o:
                    err_lc[t].append((np.log10(max(r[t], 1e-6)) - np.log10(o[t])) ** 2)

    def rmse(x):
        return math.sqrt(np.mean(x)) if x else float("nan")

    metrics = {
        "rmse_root_lc": rmse(err_lc["root"]),
        "rmse_straw_lc": rmse(err_lc["straw"]),
        "rmse_grain_lc": rmse(err_lc["grain"]),
        "straw_ratio_PFUnDA": rows[3][3]["straw"] / rows[3][4]["straw"],
        "straw_ratio_PFDoDA": rows[4][3]["straw"] / rows[4][4]["straw"],
    }
    return rows, metrics


def main():
    rows, m = run()
    print("COMPLETE long-chain recipe = 2-pool + LC6 root-matching carrier + fitted lipid\n")
    print(f"{'PFAS':7}{'nC':>3}{'Vmax x':>8} | "
          f"{'root p/o':>16}{'straw p/o':>16}{'grain p/o':>16}")
    for name, nC, mult, r, o in rows:
        print(f"{name:7}{nC:>3}{mult:>7.1f}x | "
              f"{r['root']:6.1f}/{o['root']:<7.1f}"
              f"{r['straw']:6.1f}/{o['straw']:<7.1f}"
              f"{r['grain']:6.1f}/{o['grain']:<7.1f}")
    print(f"\nlong-chain (nC>=10) log10 RMSE  "
          f"root {m['rmse_root_lc']:.3f}  straw {m['rmse_straw_lc']:.3f}  grain {m['rmse_grain_lc']:.3f}")
    print(f"straw OVER-feed ratio (pred/obs)  PFUnDA {m['straw_ratio_PFUnDA']:.2f}x  "
          f"PFDoDA {m['straw_ratio_PFDoDA']:.2f}x")
    print("\nVerdict: the carrier that closes the ROOT over-feeds the SHOOT (lipid g_xy>=0 "
          "cannot\nsubtract), so root+grain close but straw does NOT close simultaneously. "
          "The 'complete\nresolution' is not a simultaneous closure -- the long chains need a "
          "root->shoot decoupling.")
    return rows, m


if __name__ == "__main__":
    main()
