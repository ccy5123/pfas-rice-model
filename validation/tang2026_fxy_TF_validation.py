#!/usr/bin/env python3
"""
Tang 2026 transfer-factor (TF) — out-of-sample f_xy validation
==============================================================

Tests whether the model's root->shoot loading f_xy (the `f_xy_recommended` monotone
TSCF + the PFSA exp(-1.1) / ether exp(-0.7) head-group offsets) reproduces the
per-organ transfer factors MEASURED by Tang 2026 (JHM 502:141017, Table S8).
Reads the CANONICAL extraction docs/literature_db/raw_si/tang2026_doseresponse.csv
(all 5 soil doses 0.1-100 ug/g). Independent of the Yamazaki/Kim calibration that
set f_xy, so this is a genuine out-of-sample check.

DOSE CONDITION: Tang TF decreases with dose (toxicity/saturation). Here we use the
across-dose MEAN (matching validation/tang2026_validation.py); the f_xy head-group
offset itself was derived from the lowest dose 0.1 ug/g (raw_si/tang2026_tf_bcf.csv,
environmentally closest). NOTE this 4-pool check OVERLAPS validation/tang2026_validation.py
and re-confirms its documented finding (stem pass-through under, leaf-sink over);
the *new* contribution is the explicit f_xy verdict + the ORYZA biomass driver. The
structural fix + f_xy re-calibration live in validation/tang2026_fxy_refit.py and
docs/VALIDATION_TANG2026_NSTEM_KR.md.

DRIVER: the mechanistic ORYZA2000 biomass M_s(t) (`oryza_growth.oryza_drivers`),
NOT the logistic growth_rice -- the biomass driver changes the leaf TF ~2x
(ORYZA leaf senescence raises late-season leaf concentration).

UNITS: Tang TF is DRY-WEIGHT (C_organ_dw/C_root_dw); the model concentrations are
FRESH-WEIGHT (C = B_k*Cw, basis A).  A ratio of two tissues does not cancel the
water content, so each organ is converted C_dw = C_fw/(1-theta_fw) before the ratio:
    TF_dw = TF_fw * (1 - theta_root) / (1 - theta_organ).
Organ map: model stem<->Tang stalk, leaf<->leaf, grain(brown)<->endosperm (edible).

Run:  python validation/tang2026_fxy_TF_validation.py
      -> table + f_xy verdict, figure validation/figures/tang2026_fxy_TF.png
"""
import csv
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
import model_api as api            # noqa: E402
import oryza_growth as og          # noqa: E402

SEASON = 120.0
CONGENERS = ["PFOA", "PFOS", "GenX"]
# model organ -> (theta key, Tang organ name)
ORGAN_MAP = {"stem": ("stem", "stalk"), "leaf": ("leaf", "leaf"), "grain": ("grain_brown", "endosperm")}


def load_theta():
    PAR = json.load(open(os.path.join(ROOT, "params", "parameters.json")))
    tc = PAR["tissue_composition_recommended"]
    return {k: tc[k]["theta_fw"] for k in tc}


TANG_DOSERESPONSE = os.path.join(ROOT, "docs", "literature_db", "raw_si", "tang2026_doseresponse.csv")
_TF_ENDPOINT = {"TF_stalk": "stalk", "TF_leaf": "leaf", "TF_endosperm": "endosperm"}


def load_tang_tf():
    """Canonical Tang TF (SI S8) from raw_si/tang2026_doseresponse.csv.
    -> dict[(congener, organ)] = {dose_ugg: TF_mean}."""
    d = {}
    with open(TANG_DOSERESPONSE, newline="") as f:
        for r in csv.DictReader(x for x in f if not x.lstrip().startswith("#")):
            org = _TF_ENDPOINT.get(r["endpoint"])
            if org:
                d.setdefault((r["compound"], org), {})[float(r["dose_ugg"])] = float(r["value"])
    return d


def model_tf_dw(theta):
    """Run the 4-pool ODE with ORYZA biomass; return dict[congener][organ]=TF_dw."""
    out = {}
    for cong in CONGENERS:
        drv = og.oryza_drivers(cong, Cwo=1.0, season=SEASON, p=og.OryzaParams(season=SEASON))
        bf = api.simulate(cong, drivers=drv)["baf_final"]      # fw conc at Cwo=1
        froot = 1.0 - theta["root"]
        tf = {}
        for mk, (tk, _) in ORGAN_MAP.items():
            tf_fw = bf[mk] / bf["root"]
            tf[mk] = tf_fw * froot / (1.0 - theta[tk])
        out[cong] = tf
    return out


def main():
    theta = load_theta()
    tang = load_tang_tf()
    model = model_tf_dw(theta)

    rows = []
    print(f"{'cong':6}{'organ':10}{'model_TF_dw':>12}{'Tang_mean':>11}{'Tang@0.1':>10}{'log10 resid':>13}")
    resid = []
    for cong in CONGENERS:
        for mk, (tk, tang_org) in ORGAN_MAP.items():
            mt = model[cong][mk]
            obs = tang.get((cong, tang_org), {})
            obs_mean = float(np.mean(list(obs.values()))) if obs else np.nan
            obs_low = obs.get(0.1, np.nan)
            r = np.log10(max(mt, 1e-9)) - np.log10(max(obs_mean, 1e-9)) if obs else np.nan
            if np.isfinite(r):
                resid.append(r)
            rows.append((cong, tang_org, mt, obs_mean, obs_low))
            print(f"{cong:6}{tang_org:10}{mt:>12.3f}{obs_mean:>11.3f}{obs_low:>10.3f}{r:>13.2f}")
    rmse = float(np.sqrt(np.mean(np.square(resid))))
    print(f"\nlog10 RMSE (model vs Tang mean, all organs) = {rmse:.2f}")

    # f_xy verdict: congener ordering at the shoot (leaf) + the diagnostic ratios
    print("\n--- f_xy verdict ---")
    for org in ("stalk", "leaf"):
        mk = {"stalk": "stem", "leaf": "leaf"}[org]
        mo = {c: model[c][mk] for c in CONGENERS}
        to = {c: float(np.mean(list(tang[(c, org)].values()))) for c in CONGENERS}
        print(f" {org}: model order {sorted(mo, key=mo.get)}  |  Tang order {sorted(to, key=to.get)}")
        print(f"        PFOS/PFOA  model {mo['PFOS']/mo['PFOA']:.2f} vs Tang {to['PFOS']/to['PFOA']:.2f}"
              f"   GenX/PFOA model {mo['GenX']/mo['PFOA']:.2f} vs Tang {to['GenX']/to['PFOA']:.2f}")

    _figure(model, tang, rmse)


def _figure(model, tang, rmse):
    col = {"PFOA": "#1f77b4", "PFOS": "#d62728", "GenX": "#2ca02c"}
    mark = {"stalk": "o", "leaf": "s", "endosperm": "^"}
    fig, ax = plt.subplots(figsize=(6.4, 6.0))
    lo, hi = 1e-3, 1e2
    ax.plot([lo, hi], [lo, hi], "k-", lw=1, label="1:1")
    for f in (10, 100):
        ax.plot([lo, hi], [lo * f, hi * f], "k:", lw=0.6, alpha=0.5)
        ax.plot([lo, hi], [lo / f, hi / f], "k:", lw=0.6, alpha=0.5)
    for cong in CONGENERS:
        for mk, (_, tang_org) in ORGAN_MAP.items():
            obs = tang.get((cong, tang_org), {})
            if not obs:
                continue
            x = float(np.mean(list(obs.values())))
            y = model[cong][mk]
            ax.scatter(x, y, c=col[cong], marker=mark[tang_org], s=80, edgecolor="k", lw=0.5,
                       label=None, zorder=3)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_xlabel("Tang 2026 measured TF$_{dw}$  (mean over exposures)")
    ax.set_ylabel("model TF$_{dw}$  (4-pool + ORYZA biomass, f_xy_recommended)")
    ax.set_title(f"Tang 2026 per-organ TF — out-of-sample f_xy check\nlog10 RMSE = {rmse:.2f}")
    # legends: congener color + organ marker
    from matplotlib.lines import Line2D
    h1 = [Line2D([], [], marker="o", color=col[c], ls="", mec="k", label=c) for c in CONGENERS]
    h2 = [Line2D([], [], marker=mark[o], color="grey", ls="", mec="k", label=o) for o in mark]
    leg1 = ax.legend(handles=h1, title="congener", loc="upper left", fontsize=8)
    ax.add_artist(leg1)
    ax.legend(handles=h2, title="organ", loc="lower right", fontsize=8)
    ax.grid(alpha=0.25, which="both")
    out = os.path.join(HERE, "figures", "tang2026_fxy_TF.png")
    fig.tight_layout(); fig.savefig(out, dpi=140)
    print(f"\n[written] {out}")


if __name__ == "__main__":
    main()
