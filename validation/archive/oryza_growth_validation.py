#!/usr/bin/env python3
"""
ORYZA2000 Level-1 growth core — validation & driver-sensitivity
===============================================================

(1) Sanity-checks the Python ORYZA Level-1 potential run against IR72 field
    anchors (flowering ~DVS1, maturity ~DVS2, LAImax 5-7, shoot ~1740 g/m^2,
    HI ~0.45-0.50, root:shoot ~0.1).
(2) Contrasts the MECHANISTIC ORYZA biomass M_s(t) with the lightweight
    `growth_rice` DVS-partition-on-a-logistic reconstruction.
(3) Propagates both biomass drivers through the full 4-pool PFAS ODE and
    reports how the tissue BAFs shift (driver sensitivity) for a chain-length
    spread of congeners.

Run:  python validation/oryza_growth_validation.py
      -> validation/figures/oryza_growth.png  (+ console table)
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
import growth_rice as gr                       # noqa: E402
import oryza_growth as oz                      # noqa: E402
import model_api as api                        # noqa: E402

SEASON = 120.0
ORGANS = ["root", "stem", "leaf", "grain"]
COL = {"root": "saddlebrown", "stem": "goldenrod", "leaf": "forestgreen", "grain": "darkorange"}
CONGENERS = ["PFBA", "PFHxA", "PFOA", "PFDA", "PFOS"]
F_XY_SOURCE = "recommended"   # held fixed; only the biomass DRIVER changes


def sanity(p=None):
    p = p or oz.OryzaParams()
    sim = oz.simulate_oryza(p)
    t = sim["t"]
    i_fl = int(np.argmax(sim["dvs"] >= 1.0))
    i_mat = int(np.argmax(sim["dvs"] >= 2.0)) or (len(t) - 1)
    shoot = sim["wlv"] + sim["wst"] + sim["wso"]
    hi = sim["wso"][-1] / shoot[-1]
    rs = sim["wrt"][-1] / shoot[-1]
    print("== ORYZA2000 Level-1 (Python) — IR72 potential sanity ==")
    print(f"  flowering DVS=1   day {t[i_fl]:.0f}   (target ~65)")
    print(f"  maturity  DVS=2   day {t[i_mat]:.0f}   (target ~115)")
    print(f"  LAI max           {sim['lai'].max():.2f}    (target 5-7)")
    print(f"  shoot (unscaled)  {shoot[-1]/10:.0f} g/m^2 (anchor 1740; scaled at output)")
    print(f"  HI                {hi:.2f}    (target 0.45-0.50)")
    print(f"  root:shoot        {rs:.2f}    (target ~0.1)")
    checks = {
        "flowering 55-75 d": 55 <= t[i_fl] <= 75,
        "maturity 105-120 d": 105 <= t[i_mat] <= 120,
        "LAImax 4.5-7.5": 4.5 <= sim["lai"].max() <= 7.5,
        "HI 0.40-0.52": 0.40 <= hi <= 0.52,
    }
    for k, v in checks.items():
        print(f"   [{'PASS' if v else 'FAIL'}] {k}")
    return sim, all(checks.values())


def main():
    t = np.linspace(0.0, SEASON, 241)
    sim, ok = sanity()

    b_oz = oz.organ_biomass_oryza(t)                 # mechanistic (anchored)
    b_gr = gr.organ_biomass(t, SEASON)               # logistic-table reconstruction

    print("\n== final biomass [kg/hill] @ day 120 ==")
    print(f"  {'organ':6s} {'ORYZA':>10s} {'growth_rice':>12s}")
    for org in ORGANS:
        print(f"  {org:6s} {b_oz[org][-1]:10.4f} {b_gr[org][-1]:12.4f}")

    # propagate both drivers through the PFAS ODE
    print(f"\n== tissue BAF: ORYZA biomass driver vs growth_rice  (f_xy={F_XY_SOURCE}) ==")
    print(f"  {'cong':6s} | {'root oz/gr':>16s} {'straw oz/gr':>16s} {'grain oz/gr':>16s}")
    rows = {}
    for cong in CONGENERS:
        dr = oz.oryza_drivers(cong, Cwo=1.0, season=SEASON, n_t=241)
        r_oz = api.simulate(cong, drivers=dr, f_xy_source=F_XY_SOURCE)
        r_gr = api.simulate(cong, measured_forcing=True, season=SEASON, f_xy_source=F_XY_SOURCE)
        rows[cong] = (r_oz, r_gr)
        def f(r, k): return r["straw_baf"] if k == "straw" else r["baf_final"][k]
        print(f"  {cong:6s} | {f(r_oz,'root'):7.3f}/{f(r_gr,'root'):<8.3f}"
              f"{f(r_oz,'straw'):7.3f}/{f(r_gr,'straw'):<8.3f}"
              f"{f(r_oz,'grain'):7.3f}/{f(r_gr,'grain'):<8.3f}")

    # ---- figure ----
    fig, ax = plt.subplots(2, 2, figsize=(13, 10))

    for org in ORGANS:
        ax[0, 0].plot(t, b_oz[org] * 1e3, color=COL[org], lw=2, label=org)
        ax[0, 1].plot(t, b_gr[org] * 1e3, color=COL[org], lw=2, label=org)
    ymax = max(max(b_oz[o].max() for o in ORGANS), max(b_gr[o].max() for o in ORGANS)) * 1e3 * 1.05
    for a, ttl in ((ax[0, 0], "ORYZA2000 L1 (mechanistic, carbon-driven)"),
                   (ax[0, 1], "growth_rice (DVS-table on a logistic)")):
        a.set_title(ttl, fontsize=10, fontweight="bold")
        a.set_ylabel("M_k(t)  [g/hill]"); a.set_xlabel("day"); a.set_ylim(0, ymax)
        a.grid(alpha=0.25); a.legend(fontsize=8, ncol=2)

    # LAI & DVS
    td = sim["t"]
    a2 = ax[1, 0]; a2b = a2.twinx()
    a2.plot(td, sim["lai"], color="forestgreen", lw=2, label="LAI")
    a2b.plot(td, sim["dvs"], color="purple", lw=2, ls="--", label="DVS")
    a2.axhline(0, color="k", lw=0.5)
    a2b.axhline(1.0, color="grey", ls=":", lw=1); a2b.axhline(2.0, color="grey", ls=":", lw=1)
    a2.set_ylabel("LAI  [m²/m²]", color="forestgreen"); a2.set_xlabel("day")
    a2b.set_ylabel("DVS", color="purple")
    a2.set_title("ORYZA canopy & phenology (radiation/temperature-driven)", fontsize=10, fontweight="bold")
    a2.grid(alpha=0.25)

    # BAF grouped bars (grain & root)
    a3 = ax[1, 1]
    x = np.arange(len(CONGENERS)); w = 0.2
    g_oz = [rows[c][0]["baf"]["grain"][-1] for c in CONGENERS]
    g_gr = [rows[c][1]["baf"]["grain"][-1] for c in CONGENERS]
    r_oz_ = [rows[c][0]["baf"]["root"][-1] for c in CONGENERS]
    r_gr_ = [rows[c][1]["baf"]["root"][-1] for c in CONGENERS]
    a3.bar(x - 1.5 * w, r_oz_, w, color="saddlebrown", label="root ORYZA")
    a3.bar(x - 0.5 * w, r_gr_, w, color="saddlebrown", alpha=0.45, label="root growth_rice")
    a3.bar(x + 0.5 * w, g_oz, w, color="darkorange", label="grain ORYZA")
    a3.bar(x + 1.5 * w, g_gr, w, color="darkorange", alpha=0.45, label="grain growth_rice")
    a3.set_xticks(x); a3.set_xticklabels(CONGENERS); a3.set_yscale("log")
    a3.set_ylabel("BAF  [L/kg]  (log)")
    a3.set_title(f"PFAS BAF: ORYZA vs growth_rice biomass driver (f_xy={F_XY_SOURCE})",
                 fontsize=10, fontweight="bold")
    a3.grid(alpha=0.25, axis="y"); a3.legend(fontsize=7, ncol=2)

    fig.suptitle("Mechanistic ORYZA2000 biomass driver vs the logistic reconstruction, and its effect on PFAS BAF",
                 fontsize=12.5, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out = os.path.join(HERE, "figures", "oryza_growth.png")
    fig.savefig(out, dpi=130)
    print("\n[written]", out)
    print("sanity:", "PASS" if ok else "FAIL")


if __name__ == "__main__":
    main()
