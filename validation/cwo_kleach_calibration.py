#!/usr/bin/env python3
"""
Calibrate the analytic 'flooded' cwo_profile leaching rate k_leach to HYDRUS-1D
==============================================================================

`cwo_profile='flooded'` builds a time-varying pore water C_w^o(t) from an analytic
Freundlich paddy (per-congener capacity K_F = Koc·f_oc + a single first-order
leaching knob `k_leach`). With one flat default (`k_leach=0.02`) the short-chain
decline is too mild vs the real engine (`validation/cwo_profile_check.py`).

This script runs the REAL HYDRUS-1D engine for every curated congener, reads its
pore-water decline ratio (end/start, season-mean-normalised), and finds the
`k_leach` that makes the analytic 'flooded' shape match it. It then fits a simple
descriptor relationship k_leach(Koc) and writes the per-congener table to
`params/cwo_kleach.csv`, which `model_api.cwo_profile_series` loads to default
`k_leach` PER CONGENER (instead of a flat 0.02). Short chains (weakly sorbed,
low Koc) get a larger k_leach; long chains (buffered) collapse toward ~0.

Run (needs the built HYDRUS-1D engine + phydrus):
    python validation/cwo_kleach_calibration.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import model_api as api          # noqa: E402
import literature_params as lp   # noqa: E402

SEASON = 120.0
N_T = 121
K_LEACH_SCAN = np.round(np.arange(0.0, 0.151, 0.0025), 4)
OUT_CSV = os.path.join(os.path.dirname(__file__), "..", "params", "cwo_kleach.csv")


def _ratio(c):
    c = np.asarray(c, float)
    return float(c[-1] / c[0]) if c[0] > 0 else np.nan


def calibrate(congeners=None):
    congeners = congeners or api.CONGENERS
    t = np.linspace(0.0, SEASON, N_T)
    rows = []
    for cong in congeners:
        c = api._CONG[cong]
        hyd = api.cwo_profile_series(t, 1.0, "hydrus", congener=cong)   # real engine
        hyd_r = _ratio(hyd)
        # best k_leach: match the HYDRUS decline ratio in log space
        best = None
        for kl in K_LEACH_SCAN:
            flo = api.cwo_profile_series(t, 1.0, "flooded", n_C=c["n_C"], group=c["group"],
                                         congener=cong, k_leach=float(kl))
            score = abs(np.log10(max(_ratio(flo), 1e-6)) - np.log10(max(hyd_r, 1e-6)))
            if best is None or score < best[0]:
                best = (score, float(kl))
        Koc = lp.koc(c["n_C"], {"PFCA": "carboxylate", "PFSA": "sulfonate",
                                "ether": "ether"}.get(c["group"], "carboxylate"))
        rows.append(dict(congener=cong, n_C=c["n_C"], group=c["group"],
                         Koc=Koc, log10Koc=float(np.log10(Koc)),
                         hyd_ratio=hyd_r, k_leach=best[1]))
        print(f"  {cong:8} C{c['n_C']:<2} {c['group']:5} log10Koc={np.log10(Koc):5.2f}  "
              f"HYDRUS ratio {hyd_r:5.2f}  ->  k_leach {best[1]:.4f}")
    return rows


def _fit_kleach_of_logKoc(rows):
    """k_leach decreases with sorption: fit k_leach = max(0, a - b·log10Koc),
    clipped to the scan range. Reported for novel (SMILES) congeners."""
    x = np.array([r["log10Koc"] for r in rows])
    y = np.array([r["k_leach"] for r in rows])
    b, a = np.polyfit(x, y, 1)                            # y = b*x + a
    yhat = np.clip(a + b * x, 0.0, float(K_LEACH_SCAN.max()))
    rmse = float(np.sqrt(np.mean((yhat - y) ** 2)))
    return float(a), float(b), rmse


def write_csv(rows, fit):
    a, b, rmse = fit
    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    with open(OUT_CSV, "w") as f:
        f.write("# Per-congener flooded-profile k_leach, calibrated to a real HYDRUS-1D run\n")
        f.write("# (validation/cwo_kleach_calibration.py). model_api.cwo_profile_series uses\n")
        f.write("# these as the per-congener default k_leach for cwo_profile='flooded'.\n")
        f.write(f"# Novel/SMILES fallback: k_leach = clip({a:.5f} + ({b:.5f})*log10(Koc), 0, "
                f"{K_LEACH_SCAN.max():.3f})  [RMSE {rmse:.4f}]\n")
        f.write("congener,n_C,group,Koc,k_leach\n")
        for r in rows:
            f.write(f"{r['congener']},{r['n_C']},{r['group']},{r['Koc']:.4g},{r['k_leach']:.4f}\n")
    print(f"\nwrote {os.path.relpath(OUT_CSV)}  ({len(rows)} congeners)")
    print(f"fallback k_leach(log10Koc) = clip({a:.4f} + ({b:.4f})·log10Koc, 0, "
          f"{K_LEACH_SCAN.max():.3f})  RMSE {rmse:.4f}")


def run():
    if not api.hydrus_available():
        print("HYDRUS-1D engine/phydrus not available — cannot calibrate. "
              "Build the engine (see CLAUDE.md §7) and retry.")
        return None
    print("Calibrating flooded k_leach to HYDRUS-1D per congener:")
    rows = calibrate()
    fit = _fit_kleach_of_logKoc(rows)
    write_csv(rows, fit)
    return rows


if __name__ == "__main__":
    run()
