#!/usr/bin/env python3
# =============================================================================
# twopool_root_measured.py
# -----------------------------------------------------------------------------
# ROBUSTNESS check: re-fit the two-pool root + U-shaped k_seq model on the
# MEASURED forcings (forcing_rice Q_TP + growth_rice ORYZA-IR72 organ biomass)
# instead of the demo logistic placeholders. The demo transpiration peaks ~5x
# too high (CLAUDE.md), so this removes the biggest methodological caveat behind
# the twopool_root_exploration / _oos results and makes the Kim-grain OOS
# directly comparable to the fxy-doc baselines (which ARE on measured forcings:
# lipid 0.55 | monotone 2.04 | W2 1.11, Kim grain excl PFOA).
#
# Method: monkpatch the forcing globals in twopool_root_exploration (T, CWO,
# QTP, MMAT, _dM) with the measured series, then re-run its fit. EXPLORATORY /
# in-sample; canonical core + parameters.json UNCHANGED.
#
#   python validation/twopool_root_measured.py   (re-fits ~3 min)
# =============================================================================
from __future__ import annotations
import os, sys, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT_DIR, "src"))
sys.path.insert(0, HERE)

import forcing_rice as fr
from growth_rice import organ_biomass
import twopool_root_exploration as TP
from literature_params import KIM2019_FIELD, kim2019_grain_baf

SEASON = 120.0
LONG = ("PFDA", "PFUnDA", "PFDoDA")
MEAS_CACHE = os.path.join(HERE, "twopool_fitted_params_measured.json")


def install_measured_forcings():
    """Replace TP's demo forcing globals with the measured ones (in place)."""
    t = np.linspace(0.0, SEASON, 241)
    Qtp = fr.Q_TP(t, SEASON)
    b = organ_biomass(t, SEASON)
    M = np.maximum(np.column_stack([b["root"], b["stem"], b["leaf"], b["grain"]]), 1e-4)
    TP.T = t
    TP.CWO = np.full_like(t, 1.0)
    TP.QTP = Qtp
    TP.MMAT = M
    TP._dM = np.gradient(M, t, axis=0)
    return t, Qtp, M


def kim_oos(p, q):
    kim = kim2019_grain_baf("porewater")
    DF = {k: v[3] for k, v in KIM2019_FIELD.items()}
    rows = []
    for nm in kim:
        c = next((x for x in TP.CONGENERS if x["name"] == nm), None)
        if c is None:
            continue
        g = TP.simulate(c, p, kseq_override=TP.kseq_ushape(c["n_C"], c["group"], q))[2]
        rows.append((nm, c["n_C"], DF.get(nm, 0.0), kim[nm], g))
    return rows


def _rmse(rows, names):
    e = [(np.log10(max(g, 1e-6)) - np.log10(o)) ** 2
         for (nm, n, df, o, g) in rows if nm in names]
    return float(np.sqrt(np.mean(e))) if e else float("nan")


def main():
    print("=" * 78)
    print("Two-pool root + U-shaped k_seq — RE-FIT on MEASURED forcings (robustness)")
    print("=" * 78)
    t, Qtp, M = install_measured_forcings()
    print(f"measured forcings: Q_TP peak={Qtp.max():.3f} L/d/hill (demo ~0.40), "
          f"T/ET={fr.seasonal_T_over_ET(SEASON):.2f}")
    print(f"  biomass finals [kg/hill]: root={M[-1,0]:.3f} stem={M[-1,1]:.3f} "
          f"leaf={M[-1,2]:.3f} grain={M[-1,3]:.3f}  (HI={M[-1,3]/M[-1,1:].sum():.2f})")

    # re-fit on measured forcings (do NOT use the demo cache)
    print("\nre-fitting (global + root-match + U-shape) on measured forcings ...")
    p, q = TP.compute_fit_quiet()
    json.dump({"global": p, "ushape_q": list(map(float, q))}, open(MEAS_CACHE, "w"), indent=2)

    # in-sample fit quality
    print(f"\n{'PFAS':8}{'nC':>3}{'grp':>5}{'k_seq':>9} | "
          f"{'root p/o':>14}{'straw p/o':>14}{'grain p/o':>14}")
    err = {"root": [], "straw": [], "grain": []}
    for c in TP.CONGENERS:
        ks = TP.kseq_ushape(c["n_C"], c["group"], q)
        r, s, g = TP.simulate(c, p, kseq_override=ks)
        o = TP.OBS[c["name"]]; pred = {"root": r, "straw": s, "grain": g}
        for k in err:
            if k in o:
                err[k].append((np.log10(max(pred[k], 1e-6)) - np.log10(o[k])) ** 2)
        print(f"{c['name']:8}{c['n_C']:>3}{c['group'][-3:]:>5}{ks:>9.3f} | "
              f"{r:>7.2f}/{o.get('root',float('nan')):<5.2f} "
              f"{s:>7.2f}/{o.get('straw',float('nan')):<5.2f} "
              f"{g:>7.2f}/{o.get('grain',float('nan')):<5.2f}")
    allsq = err["root"] + err["straw"] + err["grain"]
    rt = {k: float(np.sqrt(np.mean(v))) for k, v in err.items()}
    print(f"\nin-sample log10 RMSE (measured forcings) = {np.sqrt(np.mean(allsq)):.3f}"
          f"  root={rt['root']:.3f} straw={rt['straw']:.3f} grain={rt['grain']:.3f}")
    print("   (demo-forcing fit was 0.251; fxy-doc U-shaped-K_PL-f_xy on measured = 0.286)")

    ks_os = TP.kseq_ushape(8, "PFSA", q); ks_un = TP.kseq_ushape(11, "PFCA", q)
    print(f"\nPFOS/PFUnDA separation: k_seq {ks_os:.3f} vs {ks_un:.3f} "
          f"({ks_un/ks_os:.1f}x at identical K_PL)  "
          f"{'HELD' if ks_un/ks_os > 1.5 else 'LOST'}")
    print(f"U-shape q (measured) = [{', '.join(f'{x:.3f}' for x in q)}]")
    print(f"   global: kappa_d={p['kappa_d']:.2f} L_Ph={p['L_Ph']:.2e} "
          f"gxy={p['gxy']:.4f} gph={p['gph']:.4f}")

    # Kim grain OOS on measured forcings (directly comparable to fxy-doc baselines)
    print("\n" + "-" * 78)
    print("Kim 2019 grain OOS (measured forcings) — vs fxy-doc baselines")
    print("-" * 78)
    rows = kim_oos(p, q)
    print(f"{'PFAS':8}{'nC':>3}{'DF%':>5} | {'obs':>8}{'2pool':>9}")
    for (nm, n, df, o, g) in rows:
        print(f"{nm:8}{n:>3}{df:>5.0f} | {o:>8.2f}{g:>9.2f}")
    alln = [r[0] for r in rows]
    excl = [nm for nm in alln if nm != "PFOA"]
    rel = [r[0] for r in rows if r[2] >= 15.0]
    print(f"\n2pool Kim grain RMSE:  all={_rmse(rows, alln):.2f}  "
          f"excl-PFOA={_rmse(rows, excl):.2f}  reliable={_rmse(rows, rel):.2f}")
    print("  fxy-doc baselines (measured forcings, excl PFOA): lipid 0.55 | mono 2.04 | W2 1.11")

    # robustness vs the demo-forcing fit
    print("\n" + "=" * 78)
    print("VERDICT (robustness to measured forcings):")
    held_sep = ks_un / ks_os > 1.5
    rmse_meas = np.sqrt(np.mean(allsq))
    print(f"  in-sample RMSE {rmse_meas:.3f} (demo 0.251); root {rt['root']:.3f}; "
          f"PFOS/PFUnDA separation {'HELD' if held_sep else 'LOST'} ({ks_un/ks_os:.1f}x).")
    print("  The qualitative result (monotone physical f_xy + non-K_PL U-shaped k_seq ->")
    print("  root solved + the separation) is" + (" ROBUST" if held_sep and rmse_meas < 0.45
          else " SENSITIVE") + " to the realistic biomass/transpiration.")
    print(f"  Measured-forcing fit cached -> {MEAS_CACHE}.  EXPLORATORY; parameters.json UNCHANGED.")


if __name__ == "__main__":
    main()
