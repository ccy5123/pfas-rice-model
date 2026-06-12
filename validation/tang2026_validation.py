#!/usr/bin/env python3
"""
Out-of-sample validation against Tang et al. 2026 (J. Hazard. Mater. 502:141017)
================================================================================

Tang 2026 is a controlled, full-growth-cycle (150 d, monthly sampling) paddy
soil-rice study of PFOA, PFOS, GenX at 5 soil doses (0.1-100 ug/g). Only the
head-group *sign* of this paper was used in the model build (f_xy PFSA offset),
so the TF/BCF *magnitudes* are genuinely OUT-OF-SAMPLE.

The clean, denominator-free target is the transfer factor TF = C_tissue/C_root
(SI Table S8): it is independent of the exposure basis (Tang's soil vs the model's
pore water), so model TF can be compared to Tang TF directly with no soil->porewater
conversion. We compare PFOA and PFOS (the two model congeners; GenX is an ether not
in the core 12) for the three transport variants (monotone / W2 / lipid), against the
across-dose mean and range.

What Tang adds that Yamazaki/Kim/Li could not:
  * controlled dose-response (not a single field exposure);
  * 4 shoot tissues resolved (stalk, leaf, chaff, endosperm) -> tests the model's
    per-compartment translocation, not just a lumped straw;
  * a real time axis (months; Table S6 shows a highly significant time effect) — the
    model's 150-d trajectory is shown for a qualitative temporal check (Tang's monthly
    raw values are only in Fig. 4a and are not digitised here).

Run:  python validation/tang2026_validation.py  ->  validation/figures/tang2026_validation.png
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
SEASON = 150.0                                  # Tang: 150-day cycle
_C = {"lipid": "#d62728", "mono": "#1f77b4", "W2": "#2ca02c", "obs": "#222222"}
_LBL = {"lipid": "model: lipid (opt-in)", "mono": "model: monotone f_xy", "W2": "model: W2 fit"}
# Tang tissue -> model compartment (model has one grain pool for husk+endosperm)
TISSUE_MAP = [("TF_stalk", "stem", "stalk→stem"),
              ("TF_leaf", "leaf", "leaf→leaf"),
              ("TF_endosperm", "grain", "endosperm→grain")]


def load_tang():
    """{compound: {endpoint: {dose: (mean, sd)}}} from the SI transcription."""
    out = {}
    path = os.path.join(ROOT, "docs", "literature_db", "raw_si", "tang2026_doseresponse.csv")
    with open(path) as f:
        for r in csv.DictReader(filter(lambda l: not l.startswith("#"), f)):
            out.setdefault(r["compound"], {}).setdefault(r["endpoint"], {})[float(r["dose_ugg"])] = (
                float(r["value"]), float(r["sd"]))
    return out


def model_tf(nm, mode):
    """Model TF = C_tissue/C_root at maturity (dose-independent: BAF ratio)."""
    kw = dict(measured_forcing=True, season=SEASON)
    if mode == "lipid":
        r = api.simulate(nm, lipid_loading=True, **kw)
    elif mode == "mono":
        r = api.simulate(nm, f_xy_source="recommended", **kw)
    else:
        r = api.simulate(nm, f_xy_source="W2fit", **kw)
    root = r["baf_final"]["root"]
    return {k: r["baf_final"][k] / root for k in ("stem", "leaf", "grain")}


def panel_tf(ax, nm, tang):
    tf = tang[nm]
    mt = {m: model_tf(nm, m) for m in ("mono", "W2", "lipid")}
    labels = [lab for _, _, lab in TISSUE_MAP]
    x = np.arange(len(labels))
    # Tang: across-dose mean with min..max range as the error bar
    means, lo, hi = [], [], []
    for tk, _, _ in TISSUE_MAP:
        vals = np.array([v for v, _ in tf[tk].values()])
        means.append(vals.mean()); lo.append(vals.mean() - vals.min()); hi.append(vals.max() - vals.mean())
    ax.bar(x - 0.30, means, 0.22, yerr=[lo, hi], capsize=3, color=_C["obs"],
           label="Tang 2026 (dose mean ± range)")
    for k, m in enumerate(("mono", "W2", "lipid")):
        ax.bar(x - 0.08 + k * 0.20, [mt[m][mk] for _, mk, _ in TISSUE_MAP], 0.18,
               color=_C[m], alpha=0.9, label=_LBL[m])
    ax.set_yscale("log"); ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("TF = tissue / root")
    ax.set_title(f"{nm} — TF vs Tang 2026 (out-of-sample, denominator-free)", fontsize=11)
    ax.axhline(1.0, color="#bbb", lw=0.8, ls=":")
    ax.legend(fontsize=8, ncol=1, loc="lower left")
    ax.grid(True, axis="y", which="both", alpha=0.25)


def panel_time(ax):
    """Model 150-d tissue trajectory (PFOA) — qualitative temporal check vs Tang."""
    r = api.simulate("PFOA", measured_forcing=True, season=SEASON, lipid_loading=True)
    t = r["t"]
    for k in ("root", "stem", "leaf", "grain"):
        ax.plot(t, r["conc"][k], lw=2, label=k)
    for mo in range(1, 6):
        ax.axvline(mo * 30, color="#ddd", lw=0.7)
    ax.set_xlabel("days after transplant (Tang: 150 d, monthly sampling)")
    ax.set_ylabel("tissue conc [µg/kg]  (Cwᵒ=1)")
    ax.set_title("Model PFOA 150-d trajectory (lipid)\nTang: builds over months (Table S6 time p<0.001),\n"
                 "grain fills late, upward migration within month 1", fontsize=10)
    ax.legend(fontsize=8.5, loc="upper left"); ax.grid(True, alpha=0.25)


def scores(tang):
    """log10 RMSE on the clean mappings (leaf, endosperm→grain) for PFOA+PFOS."""
    print(f"{'PFAS':6}{'tissue':16}{'Tang mean':>10}{'mono':>8}{'W2':>8}{'lipid':>8}")
    err = {m: [] for m in ("mono", "W2", "lipid")}
    err_clean = {m: [] for m in ("mono", "W2", "lipid")}
    for nm in ("PFOA", "PFOS"):
        mt = {m: model_tf(nm, m) for m in ("mono", "W2", "lipid")}
        for tk, mk, lab in TISSUE_MAP:
            to = np.mean([v for v, _ in tang[nm][tk].values()])
            print(f"{nm:6}{lab:16}{to:>10.2f}" + "".join(f"{mt[m][mk]:>8.2f}" for m in ("mono", "W2", "lipid")))
            for m in ("mono", "W2", "lipid"):
                e = (np.log10(max(mt[m][mk], 1e-6)) - np.log10(to)) ** 2
                err[m].append(e)
                if mk != "stem":                       # stalk/stem is the known weak mapping
                    err_clean[m].append(e)
    print("log10 RMSE (all 3 tissues) : " + "  ".join(f"{m}={np.sqrt(np.mean(err[m])):.2f}" for m in ("mono", "W2", "lipid")))
    print("log10 RMSE (leaf+grain only): " + "  ".join(f"{m}={np.sqrt(np.mean(err_clean[m])):.2f}" for m in ("mono", "W2", "lipid")))


def main():
    tang = load_tang()
    print("=== Tang 2026 TF (tissue/root) — model vs observed (OOS) ===")
    scores(tang)
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.2))
    panel_tf(axes[0], "PFOA", tang)
    panel_tf(axes[1], "PFOS", tang)
    panel_time(axes[2])
    fig.suptitle("Out-of-sample validation vs Tang 2026 (paddy soil–rice, 150 d, 5 doses): "
                 "TF tissue/root is denominator-free", fontsize=13, fontweight="bold", y=1.02)
    fig.tight_layout()
    out = os.path.join(HERE, "figures", "tang2026_validation.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=140, bbox_inches="tight")
    print("saved:", out)


if __name__ == "__main__":
    main()
