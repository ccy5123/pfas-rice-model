#!/usr/bin/env python3
"""
Validation summary figure (attached to docs/VALIDATION_KR.md)
=============================================================

Three panels in one figure (English labels; the doc prose is Korean):
  (A) Calibration (fit) to Yamazaki 2023: predicted vs observed BAF. The W2 fit is
      SATURATED (3 transport params = 3 observed BAF per congener) so it sits on the
      1:1 line by construction -> reproduction, NOT validation. Uses the same
      placeholder drivers as reproduce_demo.py (the forcing the W2 fit was made with)
      so the metric matches the documented log10 RMSE ~ 0.029.
  (B) Out-of-sample: Kim 2019 brown-rice grain BAF vs chain length, obs vs the three
      transport variants (lipid / monotone / W2). Same measured-forcing transfer as
      validation/oos_crossdataset.py.
  (C) Out-of-sample: Li 2025 water-independent tissue ratio TF (straw/root). Confounded
      -> inconclusive.

Run:  python validation/validation_summary.py  ->  validation/figures/validation_summary.png
"""
import csv
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
import model_api as api          # noqa: E402
import literature_params as lp   # noqa: E402
from pfas_rice_plant_module_4pool_surf import _logistic   # noqa: E402


def _repro_drivers():
    """The exact placeholder drivers reproduce_demo.py uses (the forcing the W2 fit
    was made with) so panel (A) matches the documented log10 RMSE = 0.029."""
    t = np.linspace(0.0, 120.0, 481)
    Qtp = 0.05 + 0.35 * np.exp(-((t - 75.0) ** 2) / (2 * 25.0 ** 2))
    M = np.column_stack([_logistic(t, 1e-3, 0.030, 0.10, 20.0), _logistic(t, 1e-3, 0.040, 0.10, 25.0),
                         _logistic(t, 1e-3, 0.050, 0.12, 30.0), _logistic(t, 1e-5, 0.025, 0.18, 80.0)])
    return dict(t=t, Cwo=np.full_like(t, 1.0), Qtp=Qtp, M=M)

plt.rcParams["axes.unicode_minus"] = False
_C = {"lipid": "#d62728", "mono": "#1f77b4", "W2": "#2ca02c", "obs": "#222222"}
_LBL = {"lipid": "lipid loading (opt-in)", "mono": "monotone f_xy", "W2": "W2 fit"}


def _sim(nm, mode, measured):
    kw = dict(measured_forcing=measured)
    if mode == "lipid":
        return api.simulate(nm, lipid_loading=True, **kw)
    if mode == "mono":
        return api.simulate(nm, f_xy_source="recommended", **kw)
    return api.simulate(nm, f_xy_source="W2fit", **kw)


def panel_A(ax):
    """Yamazaki 2023 calibration: predicted (W2) vs observed — saturated fit."""
    obs = {}
    with open(os.path.join(ROOT, "data_obs", "obs_baf_Yamazaki.csv")) as f:
        for r in csv.DictReader(f):
            obs.setdefault(r["compound"], {})[r["tissue"]] = float(r["baf"])
    # only the congeners that actually carry a W2 transport fit (others would fall
    # back to the monotone f_xy and are not part of the reproduction, as in reproduce_demo)
    w2 = {row["name"] for row in api.chain_table() if row.get("f_xy_W2fit") is not None}
    drv = _repro_drivers()
    o, p = [], []
    for nm in obs:
        if nm not in api.CONGENERS or nm not in w2:
            continue
        rr = api.simulate(nm, f_xy_source="W2fit", drivers=drv)   # exact reproduce_demo setup
        for tis in ("root", "straw", "grain"):
            if tis not in obs[nm]:
                continue
            pred = rr["straw_baf"] if tis == "straw" else rr["baf_final"][tis]
            o.append(obs[nm][tis]); p.append(pred)
    o, p = np.array(o), np.array(p)
    rmse = np.sqrt(np.mean((np.log10(p) - np.log10(o)) ** 2))
    lim = [min(o.min(), p.min()) * 0.5, max(o.max(), p.max()) * 2]
    ax.plot(lim, lim, "--", color="#999", lw=1, zorder=1, label="1:1")
    ax.scatter(o, p, s=44, c="#1f77b4", edgecolor="white", lw=0.6, zorder=3,
               label="root / straw / grain")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("observed BAF [L/kg]  —  Yamazaki 2023")
    ax.set_ylabel("predicted BAF [L/kg]")
    ax.set_title(f"(A) Calibration (fit): Yamazaki 2023\n"
                 f"log10 RMSE = {rmse:.3f}   ·   saturated fit = reproduction, not validation",
                 fontsize=10.5)
    ax.text(0.045, 0.95, "3 transport params = 3 obs / congener\n→ close fit is guaranteed",
            transform=ax.transAxes, va="top", fontsize=8.5, color="#555",
            bbox=dict(boxstyle="round", fc="#fff7e6", ec="#e0c890"))
    ax.legend(fontsize=8.5, loc="lower right")
    ax.grid(True, which="both", alpha=0.22)


def panel_B(ax):
    """Kim 2019 grain BAF (OOS): obs vs the three transport variants, vs chain length."""
    kim = lp.kim2019_grain_baf("porewater")
    DF = {"PFHpA": 13, "PFOA": 57, "PFNA": 20, "PFDA": 6.7, "PFUnDA": 13, "PFDoDA": 3.3}
    nC = {"PFHpA": 7, "PFOA": 8, "PFNA": 9, "PFDA": 10, "PFUnDA": 11, "PFDoDA": 12}
    names = [n for n in kim if n in api.CONGENERS]
    x = [nC[n] for n in names]
    ax.plot(x, [kim[n] for n in names], "o-", color=_C["obs"], lw=2, ms=7,
            label="observed (Kim 2019)", zorder=5)
    for m in ("lipid", "mono", "W2"):
        ax.plot(x, [_sim(n, m, measured=True)["baf_final"]["grain"] for n in names], "s--",
                color=_C[m], lw=1.6, ms=5, label=_LBL[m], alpha=0.9)
    for n in names:
        if n == "PFOA":
            ax.annotate("PFOA*\n(used in fit)", (nC[n], kim[n]), textcoords="offset points",
                        xytext=(4, 8), fontsize=7.5, color="#b00")
        elif DF.get(n, 0) < 15:
            ax.annotate(f"DF {DF[n]:.0f}%", (nC[n], kim[n]), textcoords="offset points",
                        xytext=(3, -15), fontsize=7, color="#999")
    ax.set_yscale("log")
    ax.set_xlabel("perfluorocarbon chain length")
    ax.set_ylabel("grain BAF [L/kg]")
    ax.set_title("(B) Out-of-sample: Kim 2019 brown-rice grain\n"
                 "only lipid tracks the long-chain rise (reliable PFHpA·PFNA: RMSE 0.23 vs mono 1.91)",
                 fontsize=10.5)
    ax.legend(fontsize=8.5, loc="upper left")
    ax.grid(True, which="both", alpha=0.25)


def panel_C(ax):
    """Li 2025 water-independent tissue ratio TF (straw/root) — inconclusive."""
    obs = {}
    with open(os.path.join(ROOT, "data_obs", "obs_baf_Li2025.csv")) as f:
        for r in csv.DictReader(f):
            obs.setdefault(r["compound"], {})[r["tissue"]] = float(r["baf"])
    order = [n for n in ("PFBA", "PFHxA", "PFOA", "PFBS", "PFOS") if n in obs]
    xo = np.arange(len(order))
    def tf(P): return P["straw_baf"] / P["baf_final"]["root"]
    ax.bar(xo - 0.3, [obs[n]["straw"] / obs[n]["root"] for n in order], 0.2,
           color=_C["obs"], label="observed (Li 2025)")
    for k, m in enumerate(("lipid", "mono", "W2")):
        ax.bar(xo - 0.1 + k * 0.2, [tf(_sim(n, m, measured=True)) for n in order], 0.2,
               color=_C[m], label=_LBL[m], alpha=0.9)
    ax.set_yscale("log")
    ax.set_xticks(xo); ax.set_xticklabels(order, fontsize=9)
    ax.set_ylabel("TF = straw / root   (water-independent ratio)")
    ax.set_title("(C) Out-of-sample: Li 2025 tissue ratio TF\n"
                 "water-quality & surface-sorption confounds → inconclusive", fontsize=10.5)
    ax.legend(fontsize=8, loc="upper right", ncol=2)
    ax.grid(True, axis="y", which="both", alpha=0.25)


def main():
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.4))
    panel_A(axes[0]); panel_B(axes[1]); panel_C(axes[2])
    fig.suptitle("PFAS–rice model validation:  calibration on Yamazaki 2023 (reproduction)  "
                 "+  out-of-sample tests on independent studies (partial)",
                 fontsize=13, fontweight="bold", y=1.02)
    fig.tight_layout()
    out = os.path.join(HERE, "figures", "validation_summary.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=140, bbox_inches="tight")
    print("saved:", out)


if __name__ == "__main__":
    main()
