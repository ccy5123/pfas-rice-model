#!/usr/bin/env python3
"""
Tang 2026 TF — f_xy re-calibration on the redistributed-shoot (nstem_leaf) model
================================================================================

Part A (refit) + Part B (ORYZA driver) of the Tang TF workstream.

The 4-pool model mal-distributes the shoot (stem pass-through ~0, leaf-sink huge),
so fitting f_xy to per-organ TF there is ill-posed (see tang2026_fxy_TF_validation).
The redistributed-shoot model (`simulate_nstem_leaf`) restores a sensible stem~leaf
split, so here we:
  (B) drive nstem_leaf with the MECHANISTIC ORYZA biomass (biomass_fn=ORYZA), and
  (A) re-calibrate f_xy per congener to the MEASURED Tang TF (data_obs/tang2026_TF.csv),
      holding the crop-architecture levers (stem_transp_frac/retention) and L_Ph fixed.

f_xy is fit (1-D, log-space) to the two f_xy-controlled shoot organs (stalk, leaf);
grain/endosperm is phloem-fed (L_Ph) and only reported.  TF is dw-converted to match
Tang: TF_dw = TF_fw * (1-theta_root)/(1-theta_organ).

Provenance: the fitted f_xy is reported as an OVERRIDE only -- params/parameters.json
is NOT modified (same policy as the existing Tang nstem validation).

Run:  python validation/tang2026_fxy_refit.py
      -> table (current vs fitted f_xy, RMSE before/after) + figure
"""
import csv
import json
import os
import sys

import numpy as np
from scipy.optimize import minimize_scalar
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
import model_api as api            # noqa: E402
import oryza_growth as og          # noqa: E402

CONGENERS = ["PFOA", "PFOS", "GenX"]
SEASON = 150.0                                    # Tang growth cycle
# model organ -> (theta key, Tang organ); fit on the shoot pair, report grain
SHOOT = {"stem": ("stem", "stalk"), "leaf": ("leaf", "leaf")}
ALL = {**SHOOT, "grain": ("grain_brown", "endosperm")}
# CLAUDE.md-documented recalibration targets (for context only)
DOC_RECAL = {"PFOS": 0.142, "GenX": 0.013}

_OB = lambda t, s: og.organ_biomass_oryza(t, p=og.OryzaParams(season=s))


def load_theta():
    tc = json.load(open(os.path.join(ROOT, "params", "parameters.json")))["tissue_composition_recommended"]
    return {k: tc[k]["theta_fw"] for k in tc}


def load_tang_mean():
    d = {}
    with open(os.path.join(ROOT, "data_obs", "tang2026_TF.csv"), newline="") as f:
        for r in csv.DictReader(x for x in f if not x.lstrip().startswith("#")):
            d.setdefault((r["congener"], r["organ"]), []).append(float(r["TF_mean"]))
    return {k: float(np.mean(v)) for k, v in d.items()}


def model_tf_dw(cong, theta, f_xy=None):
    """nstem_leaf + ORYZA biomass -> dw TF dict over model organs (stem/leaf/grain)."""
    r = api.simulate_nstem_leaf(cong, Cwo=1.0, season=SEASON, biomass_fn=_OB,
                                f_xy_override=f_xy)
    froot = 1.0 - theta["root"]
    return {mk: r["tf_final"][mk] * froot / (1.0 - theta[tk]) for mk, (tk, _) in ALL.items()}


def _shoot_logresid(cong, theta, tang, f_xy):
    tf = model_tf_dw(cong, theta, f_xy)
    return [np.log10(max(tf[mk], 1e-9)) - np.log10(tang[(cong, to)])
            for mk, (_, to) in SHOOT.items()]


def fit_fxy(cong, theta, tang):
    def obj(log_fxy):
        r = _shoot_logresid(cong, theta, tang, 10.0 ** log_fxy)
        return float(np.mean(np.square(r)))
    res = minimize_scalar(obj, bounds=(-3.0, 0.0), method="bounded",
                          options={"xatol": 1e-3})
    return 10.0 ** res.x


def rmse_all(cong, theta, tang, f_xy):
    tf = model_tf_dw(cong, theta, f_xy)
    r = [np.log10(max(tf[mk], 1e-9)) - np.log10(tang[(cong, to)]) for mk, (_, to) in ALL.items()]
    return float(np.sqrt(np.mean(np.square(r)))), tf


def main():
    theta = load_theta()
    tang = load_tang_mean()
    cur = {c: api._CONG[c]["f_xy_recommended"] for c in CONGENERS}

    print(f"{'cong':6}{'f_xy_cur':>10}{'f_xy_fit':>10}{'doc_recal':>11}"
          f"{'RMSE_before':>13}{'RMSE_after':>12}")
    fitted, before, after = {}, {}, {}
    tf_before, tf_after = {}, {}
    rb_all, ra_all = [], []
    for c in CONGENERS:
        fx = fit_fxy(c, theta, tang)
        fitted[c] = fx
        rb, tfb = rmse_all(c, theta, tang, cur[c])
        ra, tfa = rmse_all(c, theta, tang, fx)
        before[c], after[c], tf_before[c], tf_after[c] = rb, ra, tfb, tfa
        rb_all.append(rb); ra_all.append(ra)
        print(f"{c:6}{cur[c]:>10.4f}{fx:>10.4f}{DOC_RECAL.get(c, float('nan')):>11.4f}"
              f"{rb:>13.2f}{ra:>12.2f}")
    print(f"\noverall log10 RMSE  before={np.sqrt(np.mean(np.square(rb_all))):.2f}"
          f"  ->  after={np.sqrt(np.mean(np.square(ra_all))):.2f}")
    print("(f_xy fit to shoot stalk+leaf; grain/endosperm phloem-fed, reported not fit.")
    print(" fitted f_xy is OVERRIDE-only; params/parameters.json unchanged.)")

    _figure(theta, tang, cur, fitted, tf_before, tf_after)


def _figure(theta, tang, cur, fitted, tf_before, tf_after):
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.6), sharey=True)
    organs = [("stem", "stalk"), ("leaf", "leaf"), ("grain", "endosperm")]
    x = np.arange(len(organs))
    w = 0.27
    for ax, c in zip(axes, CONGENERS):
        obs = [tang[(c, to)] for _, to in organs]
        bef = [tf_before[c][mk] for mk, _ in organs]
        aft = [tf_after[c][mk] for mk, _ in organs]
        ax.bar(x - w, obs, w, label="Tang (obs)", color="#444")
        ax.bar(x, bef, w, label=f"model f_xy={cur[c]:.3f}", color="#bbb", edgecolor="k", lw=0.4)
        ax.bar(x + w, aft, w, label=f"refit f_xy={fitted[c]:.3f}", color="#2ca02c", edgecolor="k", lw=0.4)
        ax.set_yscale("log")
        ax.set_xticks(x); ax.set_xticklabels([to for _, to in organs])
        ax.set_title(c); ax.grid(axis="y", alpha=0.25, which="both")
        ax.legend(fontsize=7.5, loc="upper right")
    axes[0].set_ylabel("transfer factor TF$_{dw}$ = C$_{organ}$/C$_{root}$")
    fig.suptitle("Tang 2026 TF — f_xy re-calibration on nstem_leaf + ORYZA biomass "
                 "(shoot fit; grain phloem-fed)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = os.path.join(HERE, "figures", "tang2026_fxy_refit.png")
    fig.savefig(out, dpi=140)
    print(f"\n[written] {out}")


if __name__ == "__main__":
    main()
