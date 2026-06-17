#!/usr/bin/env python3
"""
Grain shape with the formation gate (after-fix)
===============================================

Tissue-dynamics concentrations with the DPU-consistent grain formation gate
(4pool_surf / nstem_leaf) in place: the grain takes NO PFAS until it forms (~flowering),
then rises monotonically from 0 to its harvest value -- no pre-flowering spike (the floored,
not-yet-set grain mass is no longer loaded). Left panel zooms the grain to its own range so
the shape is visible; right panel shows all compartments.

Run:  python validation/grain_gate_figure.py
      -> validation/figures/grain_formation_gate.png
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
import model_api as api                                   # noqa: E402

_COL = {"root": "#8c564b", "stem": "#2ca02c", "leaf": "#1f77b4", "grain": "#ff7f0e"}


def main():
    r = api.simulate("PFOA", biomass="oryza", season=120.0)
    t, cg = r["t"], np.asarray(r["conc"]["grain"])

    fig, ax = plt.subplots(1, 2, figsize=(13.5, 5.0))
    # left: grain only, y auto-scaled to its real range (no longer dwarfed by a spike)
    ax[0].plot(t, cg, color="#ff7f0e", lw=2.6)
    ax[0].set_title("grain concentration (formation-gated)")
    ax[0].set_ylabel("conc [µg/kg]")
    i0 = int(np.argmax(cg > 1e-6))
    ax[0].annotate(f"forms ~d{t[i0]:.0f}\n(rises from 0)", (t[i0], cg[i0]), (t[i0] + 7, 0.05),
                   arrowprops=dict(arrowstyle="->", color="grey"), fontsize=8)
    ax[0].annotate(f"harvest {cg[-1]:.3f}", (t[-1], cg[-1]), (96, 0.10), fontsize=8)
    # right: full tissue dynamics
    for k in api.TISSUES:
        ax[1].plot(t, r["conc"][k], color=_COL[k], lw=2.3, label=k)
    ax[1].plot(t, r["straw"], "k--", lw=1.2, label="straw")
    ax[1].set_title("full Tissue dynamics"); ax[1].legend(fontsize=8)
    for a in ax:
        a.axvline(66, color="grey", ls=":", lw=1); a.grid(alpha=0.25)
        a.set_xlabel("days after transplant")
        a.text(66, a.get_ylim()[1] * 0.93, "flowering", rotation=90, fontsize=7,
               color="grey", va="top")
    fig.suptitle("PFOA (ORYZA biomass) — grain rises from 0 at flowering (formation-gated)",
                 fontsize=12)
    fig.tight_layout()
    out = os.path.join(HERE, "figures", "grain_formation_gate.png")
    fig.savefig(out, dpi=135)
    print("grain conc: " + "  ".join(f"d{d}:{cg[np.argmin(abs(t - d))]:.3g}"
                                     for d in (50, 55, 66, 90, 120)))
    print(f"[written] {out}")


if __name__ == "__main__":
    main()
