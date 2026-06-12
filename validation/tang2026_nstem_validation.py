#!/usr/bin/env python3
"""
Tang 2026 re-validation of the REDISTRIBUTED-SHOOT fix (nstem + leaf/stem loading)
=================================================================================

The first Tang 2026 out-of-sample check (`validation/tang2026_validation.py`,
`docs/VALIDATION_TANG2026_KR.md`) showed the single-straw 4-compartment core
*over-translocates* into the shoot: the stem compartment is an empty pass-through
(stalk TF ~ 0.03 vs Tang 0.6-1.5) while the leaf is the sole xylem terminal and
runs away (leaf TF 3-13 vs Tang 0.7-1.7; the leaf holds ~81% of the plant burden).

This script re-validates the structural fix — `pfas_rice_plant_module_nstem_leaf`:
the stem is resolved into N transpiration terminals and the shoot loading is
redistributed by the transpiration-deposition+retention mechanism applied to
EVERY shoot organ (not just the leaf), so each organ retains its own transpired
solute instead of piling it into the leaf.

What it reports (model TF = C_tissue/C_root at maturity, denominator-free, so
directly comparable to Tang SI Table S8 with no soil->porewater conversion):

  * per-tissue TF (stalk/leaf/grain), Tang vs baseline single-straw vs nstem_leaf;
  * the per-compartment BURDEN redistribution (leaf monopoly -> shared);
  * an RMSE decomposition: the SHAPE RMSE (within-congener stalk/leaf/grain
    pattern, level removed) -- what the redistribution targets -- vs the overall
    RMSE (which still carries the across-congener f_xy/B_root magnitude spread);
  * a sensitivity sweep over the two structural levers (stem transpiration
    fraction, retention efficiency).

Honest result: the redistribution CURES the shoot-pattern error (shape RMSE
0.84 -> 0.11; PFOA matches Tang across all three tissues, RMSE 1.03 -> 0.06). The
residual is the across-congener LEVEL: the Yamazaki-calibrated f_xy (PFOS 0.013,
GenX 0.233) and basis-A B_root (PFOS 49) spread make PFOS too root-bound (TF too
low) and GenX too mobile (TF too high) -- a binding/translocation magnitude issue,
NOT the shoot-structure issue the redistribution fixed.

Run:  python validation/tang2026_nstem_validation.py
      -> validation/figures/tang2026_nstem_validation.png
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

plt.rcParams["axes.unicode_minus"] = False
SEASON = 150.0
KEYS = ("stem", "leaf", "grain")
TISSUE_LBL = {"stem": "stalk", "leaf": "leaf", "grain": "endosperm"}
COMPOUNDS = ("PFOA", "PFOS", "GenX")
GRP = {"PFOA": "PFCA", "PFOS": "PFSA", "GenX": "ether*"}
# nstem_leaf structural defaults (crop-architecture levers; NOT fit to Tang TF)
RETENTION, STEM_FRAC = 0.6, 0.45
_C = {"tang": "#222222", "base": "#1f77b4", "nl": "#d62728"}


def load_tang_tf():
    """Across-dose mean TF {compound: {tissue: mean}} from SI Table S8."""
    raw = {}
    path = os.path.join(ROOT, "docs", "literature_db", "raw_si", "tang2026_doseresponse.csv")
    with open(path) as f:
        for r in csv.DictReader(filter(lambda l: not l.startswith("#"), f)):
            ep = r["endpoint"]
            if not ep.startswith("TF_"):
                continue
            raw.setdefault(r["compound"], {}).setdefault(ep, []).append(float(r["value"]))
    # map Tang tissues -> model keys (model lumps husk+endosperm into one grain pool)
    tmap = {"TF_stalk": "stem", "TF_leaf": "leaf", "TF_endosperm": "grain"}
    out = {}
    for nm, eps in raw.items():
        out[nm] = {tmap[k]: float(np.mean(v)) for k, v in eps.items() if k in tmap}
    return out


def baseline_tf(nm, src="recommended"):
    r = api.simulate(nm, f_xy_source=src, measured_forcing=True, season=SEASON, n_t=361)
    root = r["baf_final"]["root"]
    return {k: r["baf_final"][k] / root for k in KEYS}


def nstem_tf(nm, retention=RETENTION, stem_frac=STEM_FRAC, **kw):
    r = api.simulate_nstem_leaf(nm, retention=retention, stem_transp_frac=stem_frac,
                                season=SEASON, **kw)
    return {k: r["tf_final"][k] for k in KEYS}


def _rmse(model, tang):
    e = [(np.log10(max(model[nm][k], 1e-6)) - np.log10(tang[nm][k])) ** 2
         for nm in tang for k in KEYS]
    return float(np.sqrt(np.mean(e)))


def _shape_rmse(model, tang):
    """Within-congener PATTERN RMSE (each congener's overall log-level removed)."""
    e = []
    for nm in tang:
        m = np.array([np.log10(max(model[nm][k], 1e-6)) for k in KEYS]); m -= m.mean()
        o = np.array([np.log10(tang[nm][k]) for k in KEYS]); o -= o.mean()
        e += list((m - o) ** 2)
    return float(np.sqrt(np.mean(e)))


def _per_cong(model, tang):
    return {nm: float(np.sqrt(np.mean([(np.log10(max(model[nm][k], 1e-6)) - np.log10(tang[nm][k])) ** 2
                                       for k in KEYS]))) for nm in tang}


def burden_fractions(nm, which):
    """Plant-burden fraction per organ (root/stalk/leaf/grain) at maturity."""
    if which == "base":
        r = api.simulate(nm, f_xy_source="recommended", measured_forcing=True, season=SEASON, n_t=361)
        Mf = r["M"][-1]
        b = {"root": r["conc"]["root"][-1] * Mf[0], "stalk": r["conc"]["stem"][-1] * Mf[1],
             "leaf": r["conc"]["leaf"][-1] * Mf[2], "grain": r["conc"]["grain"][-1] * Mf[3]}
    else:
        r = api.simulate_nstem_leaf(nm, retention=RETENTION, stem_transp_frac=STEM_FRAC, season=SEASON)
        Mf = r["M"][-1]; N = r["N"]; seg = slice(1, N + 1)
        b = {"root": r["conc"]["root"][-1] * Mf[0],
             "stalk": float(np.sum(r["conc"]["stem"][-1] * np.sum(Mf[seg]))),
             "leaf": r["conc"]["leaf"][-1] * Mf[N + 1], "grain": r["conc"]["grain"][-1] * Mf[N + 2]}
    tot = sum(b.values())
    return {k: v / tot for k, v in b.items()}


def panel_tf(ax, nm, tang, base, nl, tag):
    x = np.arange(len(KEYS))
    ax.bar(x - 0.26, [tang[nm][k] for k in KEYS], 0.24, color=_C["tang"], label="Tang 2026 (dose mean)")
    ax.bar(x + 0.00, [base[nm][k] for k in KEYS], 0.24, color=_C["base"], alpha=0.9, label="single-straw (baseline)")
    ax.bar(x + 0.26, [nl[nm][k] for k in KEYS], 0.24, color=_C["nl"], alpha=0.9, label="nstem_leaf (redistributed)")
    ax.set_yscale("log"); ax.set_xticks(x); ax.set_xticklabels([TISSUE_LBL[k] for k in KEYS])
    ax.set_ylabel("TF = tissue / root")
    ax.axhline(1.0, color="#bbb", lw=0.8, ls=":")
    ax.set_title(f"({tag}) {nm} ({GRP[nm]}) — TF vs Tang (OOS)", fontsize=10)
    ax.legend(fontsize=7, loc="upper left"); ax.grid(True, axis="y", which="both", alpha=0.25)


def panel_burden(ax):
    organs = ["root", "stalk", "leaf", "grain"]
    bb = burden_fractions("PFOA", "base"); nn = burden_fractions("PFOA", "nl")
    x = np.arange(len(organs))
    ax.bar(x - 0.2, [bb[o] for o in organs], 0.38, color=_C["base"], alpha=0.9, label="single-straw")
    ax.bar(x + 0.2, [nn[o] for o in organs], 0.38, color=_C["nl"], alpha=0.9, label="nstem_leaf")
    ax.set_xticks(x); ax.set_xticklabels(organs); ax.set_ylabel("fraction of plant PFAS burden")
    ax.set_title("(D) PFOA burden distribution\nleaf monopoly (81%) → redistributed across the shoot", fontsize=9.5)
    ax.legend(fontsize=8); ax.grid(True, axis="y", alpha=0.25)


def panel_rmse(ax, tang, base, nl):
    cats = ["shape\n(pattern)", "overall\n(9 pts)", "PFOA", "PFOS", "GenX"]
    pb, pn = _per_cong(base, tang), _per_cong(nl, tang)
    bvals = [_shape_rmse(base, tang), _rmse(base, tang), pb["PFOA"], pb["PFOS"], pb["GenX"]]
    nvals = [_shape_rmse(nl, tang), _rmse(nl, tang), pn["PFOA"], pn["PFOS"], pn["GenX"]]
    x = np.arange(len(cats))
    ax.bar(x - 0.2, bvals, 0.38, color=_C["base"], alpha=0.9, label="single-straw")
    ax.bar(x + 0.2, nvals, 0.38, color=_C["nl"], alpha=0.9, label="nstem_leaf")
    ax.set_xticks(x); ax.set_xticklabels(cats, fontsize=8); ax.set_ylabel("log10 RMSE vs Tang")
    ax.set_title("(E) RMSE decomposition\nshoot PATTERN cured (0.84→0.11); residual = across-congener level", fontsize=9.5)
    ax.legend(fontsize=8); ax.grid(True, axis="y", alpha=0.25)


def panel_sensitivity(ax, tang):
    rets = np.linspace(0.2, 1.0, 9)
    for k, col in zip(KEYS, ("#2ca02c", "#d62728", "#ff7f0e")):
        ax.plot(rets, [nstem_tf("PFOA", retention=r)[k] for r in rets], "-o", ms=3, color=col,
                label=f"model {TISSUE_LBL[k]}")
        ax.axhline(tang["PFOA"][k], color=col, lw=1.0, ls=":", alpha=0.8)
    ax.axvline(RETENTION, color="#888", lw=0.9, ls="--")
    ax.set_xlabel("retention efficiency"); ax.set_ylabel("PFOA TF = tissue/root")
    ax.set_title("(F) Retention sweep (PFOA)\ndotted = Tang; dashed = default 0.6 (matches all 3 tissues)", fontsize=9.5)
    ax.legend(fontsize=7.5, ncol=3, loc="upper center"); ax.grid(True, alpha=0.25)


def print_table(tang, base, nl):
    print("=== Tang 2026 TF (tissue/root): baseline single-straw vs nstem_leaf (OOS) ===")
    print(f"{'PFAS':6}{'tissue':10}{'Tang':>7}{'baseline':>10}{'nstem_leaf':>12}")
    for nm in COMPOUNDS:
        for k in KEYS:
            print(f"{nm:6}{TISSUE_LBL[k]:10}{tang[nm][k]:>7.2f}{base[nm][k]:>10.2f}{nl[nm][k]:>12.2f}")
    print("\nlog10 RMSE vs Tang:")
    print(f"  overall (9 pts) : baseline {_rmse(base,tang):.2f}  ->  nstem_leaf {_rmse(nl,tang):.2f}")
    print(f"  SHAPE (pattern) : baseline {_shape_rmse(base,tang):.2f}  ->  nstem_leaf {_shape_rmse(nl,tang):.2f}")
    pb, pn = _per_cong(base, tang), _per_cong(nl, tang)
    for nm in COMPOUNDS:
        print(f"  {nm:6}        : baseline {pb[nm]:.2f}  ->  nstem_leaf {pn[nm]:.2f}")
    bb = burden_fractions("PFOA", "base"); nn = burden_fractions("PFOA", "nl")
    print(f"\nPFOA leaf burden fraction: baseline {bb['leaf']:.0%}  ->  nstem_leaf {nn['leaf']:.0%}")
    print(f"           stalk fraction: baseline {bb['stalk']:.0%}  ->  nstem_leaf {nn['stalk']:.0%}")


def main():
    tang = load_tang_tf()
    base = {nm: baseline_tf(nm) for nm in COMPOUNDS}
    nl = {nm: nstem_tf(nm) for nm in COMPOUNDS}
    print_table(tang, base, nl)

    fig, axes = plt.subplots(2, 3, figsize=(18.5, 10.0))
    panel_tf(axes[0, 0], "PFOA", tang, base, nl, "A")
    panel_tf(axes[0, 1], "PFOS", tang, base, nl, "B")
    panel_tf(axes[0, 2], "GenX", tang, base, nl, "C")
    panel_burden(axes[1, 0])
    panel_rmse(axes[1, 1], tang, base, nl)
    panel_sensitivity(axes[1, 2], tang)
    fig.suptitle("Tang 2026 re-validation: redistributed-shoot fix (N-segment stem + leaf/stem transpiration "
                 "deposition+retention)\nshoot tissue PATTERN cured; residual = across-congener f_xy/B_root level",
                 fontsize=13, fontweight="bold", y=1.005)
    fig.tight_layout()
    out = os.path.join(HERE, "figures", "tang2026_nstem_validation.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=140, bbox_inches="tight")
    print("\nsaved:", out)


if __name__ == "__main__":
    main()
