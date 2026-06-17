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
F_OC = 0.016                                    # Tang soil: OM 27.4 g/kg (Table S2) -> ~1.6% OC
THETA_G = 0.35                                  # gravimetric water content [L/kg]
THETA_FW = np.array([0.90, 0.83, 0.78, 0.14])   # root/stem/leaf/grain fresh-water fraction
_C = {"lipid": "#d62728", "mono": "#1f77b4", "W2": "#2ca02c", "obs": "#222222"}
_LBL = {"lipid": "model: lipid (opt-in)", "mono": "model: monotone f_xy", "W2": "model: W2 fit"}
# Tang tissue -> model compartment (model has one grain pool for husk+endosperm)
TISSUE_MAP = [("TF_stalk", "stem", "stalk→stem"),
              ("TF_leaf", "leaf", "leaf→leaf"),
              ("TF_endosperm", "grain", "endosperm→grain")]
CONG = {"PFOA": (8, "PFCA"), "PFOS": (8, "PFSA"), "GenX": (5, "ether")}
COMPOUNDS = ("PFOA", "PFOS", "GenX")
TANG_BCF_RANGE = {"PFOA": (0.181, 0.240), "PFOS": (0.217, 0.295),   # Table S7 min..max over doses
                  "GenX": (0.358, 0.523)}


def load_tang():
    """{compound: {endpoint: {dose: (mean, sd)}}} from the SI transcription."""
    out = {}
    path = os.path.join(ROOT, "docs", "literature_db", "raw_si", "tang2026_doseresponse.csv")
    with open(path) as f:
        for r in csv.DictReader(filter(lambda l: not l.startswith("#"), f)):
            out.setdefault(r["compound"], {}).setdefault(r["endpoint"], {})[float(r["dose_ugg"])] = (
                float(r["value"]), float(r["sd"]))
    return out


def _run(nm, mode):
    kw = dict(measured_forcing=True, season=SEASON)
    if mode == "lipid":
        return api.simulate(nm, lipid_loading=True, **kw)
    if mode == "mono":
        return api.simulate(nm, f_xy_source="recommended", **kw)
    return api.simulate(nm, f_xy_source="W2fit", **kw)


def model_tf(nm, mode):
    """Model TF = C_tissue(dw)/C_root(dw) at maturity, on Tang's DRY-WEIGHT basis.

    TF is independent of the EXPOSURE basis (Cw cancels: Tang soil vs model pore water),
    but NOT of the tissue-moisture basis: the model concentrations are FRESH-weight
    (C = B_k·Cw, basis A) while Tang's TF is dry/dry, and the fresh->dry factor (1−θ_fw)
    differs between tissues (root θ=0.90 vs grain θ=0.14), so it does NOT cancel. It must
    be applied:  TF_dw = TF_fw · (1−θ_root)/(1−θ_tissue).
    [Correction: earlier this returned the FRESH-weight ratio, which understated what
    Tang's dry-weight TF requires by (1−θ_root)/(1−θ_tissue) — ~0.59 stem, ~0.45 leaf,
    ~0.12 grain. The "grain matches" result was a fw/dw artifact; see
    docs/tang2026_grain_units_exploration.md.]"""
    r = _run(nm, mode)
    th = {"stem": THETA_FW[1], "leaf": THETA_FW[2], "grain": THETA_FW[3]}
    froot = 1.0 - THETA_FW[0]
    root = r["baf_final"]["root"]
    return {k: (r["baf_final"][k] / root) * froot / (1.0 - th[k]) for k in ("stem", "leaf", "grain")}


def model_bcf(nm, mode):
    """Model BCF = C_rice(dw)/C_soil, denominator built from the soil Kd (Koc·f_oc).

    BCF = [Σ M·BAF / Σ M·(1−θ_fw)] / (Kd + θ_g)   (the Cw cancels; dw whole-plant
    BAF over the soil partition Kd). Kd from the chain-length Koc QSPR at Tang's f_oc.
    """
    import soil_hydrus as sh
    r = _run(nm, mode)
    BAF = np.array([r["baf_final"][k] for k in ("root", "stem", "leaf", "grain")])
    Mf = r["M"][-1]
    dw_baf = float(np.sum(Mf * BAF) / np.sum(Mf * (1 - THETA_FW)))
    n_C, grp = CONG[nm]
    Kd = sh.paddy_kd(n_C, grp, f_oc=F_OC)
    return dw_baf / (Kd + THETA_G), Kd


def model_tf_traj(nm, mode):
    """C_leaf/C_root and C_grain/C_root over the season (temporal TF)."""
    r = _run(nm, mode)
    t = r["t"]
    root = np.maximum(r["conc"]["root"], 1e-9)
    return t, r["conc"]["leaf"] / root, r["conc"]["grain"] / root


def panel_tf(ax, nm, tang, tag):
    tf = tang[nm]
    mt = {m: model_tf(nm, m) for m in ("mono", "W2", "lipid")}
    labels = [lab for _, _, lab in TISSUE_MAP]
    x = np.arange(len(labels))
    means, lo, hi = [], [], []
    for tk, _, _ in TISSUE_MAP:
        vals = np.array([v for v, _ in tf[tk].values()])
        means.append(vals.mean()); lo.append(vals.mean() - vals.min()); hi.append(vals.max() - vals.mean())
    ax.bar(x - 0.30, means, 0.22, yerr=[lo, hi], capsize=3, color=_C["obs"],
           label="Tang 2026 (dose mean ± range)")
    for k, m in enumerate(("mono", "W2", "lipid")):
        ax.bar(x - 0.08 + k * 0.20, [mt[m][mk] for _, mk, _ in TISSUE_MAP], 0.18,
               color=_C[m], alpha=0.9, label=_LBL[m])
    ax.set_yscale("log"); ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=8.5)
    ax.set_ylabel("TF = tissue / root")
    grp = {"PFOA": "PFCA", "PFOS": "PFSA", "GenX": "ether*"}[nm]
    ax.set_title(f"({tag}) {nm} ({grp}) — TF vs Tang 2026 (OOS, denominator-free)", fontsize=10)
    ax.axhline(1.0, color="#bbb", lw=0.8, ls=":")
    ax.legend(fontsize=7, ncol=1, loc="lower left")
    ax.grid(True, axis="y", which="both", alpha=0.25)


def panel_bcf(ax):
    """Absolute BCF = C_rice/C_soil: model (soil Kd from Koc·f_oc) vs Tang range."""
    names = list(COMPOUNDS)
    x = np.arange(len(names))
    lo = [TANG_BCF_RANGE[n][0] for n in names]; hi = [TANG_BCF_RANGE[n][1] for n in names]
    mean = [(a + b) / 2 for a, b in zip(lo, hi)]
    ax.bar(x - 0.30, mean, 0.22, yerr=[[m - l for m, l in zip(mean, lo)], [h - m for h, m in zip(hi, mean)]],
           capsize=4, color=_C["obs"], label="Tang 2026 (dose range)")
    for k, m in enumerate(("mono", "W2", "lipid")):
        ax.bar(x - 0.08 + k * 0.20, [model_bcf(n, m)[0] for n in names], 0.18,
               color=_C[m], alpha=0.9, label=_LBL[m])
    ax.set_yscale("log"); ax.set_xticks(x); ax.set_xticklabels(names)
    ax.set_ylabel("BCF = C_rice(dw) / C_soil")
    kd = {n: model_bcf(n, "mono")[1] for n in names}
    ax.set_title(f"(D) Absolute BCF vs Tang  (f_oc={F_OC}; Kd: PFOA {kd['PFOA']:.1f}, "
                 f"PFOS {kd['PFOS']:.1f}, GenX {kd['GenX']:.2f} L/kg)\norder-of-magnitude OK; "
                 "model over-predicts uptake", fontsize=9.5)
    ax.legend(fontsize=7, loc="upper left")
    ax.grid(True, axis="y", which="both", alpha=0.25)


def panel_fxy(ax):
    """Head-group / chain effect on translocation f_xy (why GenX is the most mobile)."""
    fxy = {nm: api._CONG[nm]["f_xy_recommended"] for nm in COMPOUNDS}
    cols = ["#1f77b4", "#2ca02c", "#9467bd"]
    ax.bar(COMPOUNDS, [fxy[n] for n in COMPOUNDS], color=cols, alpha=0.9)
    for i, n in enumerate(COMPOUNDS):
        ax.text(i, fxy[n] * 1.05, f"{fxy[n]:.3f}", ha="center", fontsize=9)
    ax.set_ylabel("f_xy  (root→shoot loading)")
    ax.set_title("(F) Translocation lever f_xy (head-group × chain)\n"
                 "GenX (ether, short) > PFOA > PFOS — matches Tang's strong GenX upward migration",
                 fontsize=9.5)
    ax.grid(True, axis="y", alpha=0.25)


def panel_traj(ax):
    """Temporal TF(t): does the model reproduce Tang's root-dominance early?"""
    for mode, c in (("mono", _C["mono"]), ("lipid", _C["lipid"])):
        t, lr, gr = model_tf_traj("PFOA", mode)
        ax.plot(t, lr, color=c, lw=2, label=f"leaf/root ({mode})")
    t, lr, gr = model_tf_traj("PFOA", "lipid")
    ax.plot(t, np.clip(gr, 0, 30), color="#ff7f0e", lw=1.6, ls="--", label="grain/root (lipid)")
    ax.axhline(1.0, color="#444", lw=1.0, ls=":")
    ax.axhspan(0.01, 1.0, color="#e8f0e8", alpha=0.7)
    ax.text(8, 0.4, "root-dominant\n(Tang: month 1)", fontsize=8, color="#3a7a3a")
    for mo in range(1, 6):
        ax.axvline(mo * 30, color="#eee", lw=0.7)
    ax.set_yscale("log"); ax.set_ylim(0.05, 40)
    ax.set_xlabel("days after transplant (Tang: 150 d, monthly)")
    ax.set_ylabel("model TF = tissue / root")
    ax.set_title("(E) Temporal TF(t), PFOA — model vs Tang's root-first finding\n"
                 "Tang: roots dominate after month 1 (TF<1); model keeps leaf/root≫1 → over-translocation",
                 fontsize=9.5)
    ax.legend(fontsize=8, loc="upper right"); ax.grid(True, which="both", alpha=0.2)


def scores(tang):
    print(f"{'PFAS':6}{'tissue':16}{'Tang mean':>10}{'mono':>8}{'W2':>8}{'lipid':>8}")
    err = {m: [] for m in ("mono", "W2", "lipid")}
    err_clean = {m: [] for m in ("mono", "W2", "lipid")}
    for nm in COMPOUNDS:
        mt = {m: model_tf(nm, m) for m in ("mono", "W2", "lipid")}
        for tk, mk, lab in TISSUE_MAP:
            to = np.mean([v for v, _ in tang[nm][tk].values()])
            print(f"{nm:6}{lab:16}{to:>10.2f}" + "".join(f"{mt[m][mk]:>8.2f}" for m in ("mono", "W2", "lipid")))
            for m in ("mono", "W2", "lipid"):
                e = (np.log10(max(mt[m][mk], 1e-6)) - np.log10(to)) ** 2
                err[m].append(e)
                if mk != "stem":
                    err_clean[m].append(e)
    print("log10 RMSE (all, 3 cmpd × 3 tissue): " + "  ".join(f"{m}={np.sqrt(np.mean(err[m])):.2f}" for m in ("mono", "W2", "lipid")))
    print("log10 RMSE (leaf+grain only)       : " + "  ".join(f"{m}={np.sqrt(np.mean(err_clean[m])):.2f}" for m in ("mono", "W2", "lipid")))
    print("\nAbsolute BCF = C_rice/C_soil (f_oc=%.3f):" % F_OC)
    for nm in COMPOUNDS:
        b = {m: model_bcf(nm, m)[0] for m in ("mono", "W2", "lipid")}
        print(f"  {nm}: Tang {TANG_BCF_RANGE[nm][0]:.2f}-{TANG_BCF_RANGE[nm][1]:.2f}  | "
              + "  ".join(f"{m}={b[m]:.2f}" for m in ("mono", "W2", "lipid")))


def main():
    tang = load_tang()
    print("=== Tang 2026 TF (tissue/root) — model vs observed (OOS) ===")
    scores(tang)
    fig, axes = plt.subplots(2, 3, figsize=(18.5, 9.8))
    panel_tf(axes[0, 0], "PFOA", tang, "A")
    panel_tf(axes[0, 1], "PFOS", tang, "B")
    panel_tf(axes[0, 2], "GenX", tang, "C")
    panel_bcf(axes[1, 0])
    panel_traj(axes[1, 1])
    panel_fxy(axes[1, 2])
    fig.suptitle("Out-of-sample validation vs Tang 2026 (paddy soil–rice, 150 d, 5 doses; PFOA / PFOS / GenX):  "
                 "TF (denominator-free), absolute BCF, temporal TF(t)", fontsize=13, fontweight="bold", y=1.01)
    fig.tight_layout()
    out = os.path.join(HERE, "figures", "tang2026_validation.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=140, bbox_inches="tight")
    print("saved:", out)


if __name__ == "__main__":
    main()
