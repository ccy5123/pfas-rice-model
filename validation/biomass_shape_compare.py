#!/usr/bin/env python3
"""
Per-compartment biomass shape diagnostic: growth_rice vs ORYZA2000
==================================================================

Plots and quantifies the SHAPE of each organ's biomass trajectory M_k(t) over the
season for the two biomass drivers — `growth_rice` (ORYZA IR72 DVS-partitioning on a
logistic total-biomass curve) and `oryza_growth` (the mechanistic ORYZA2000 Level-1
carbon balance) — so the difference that matters for the PFAS model is explicit.

Key finding (season 120 d, dry g/hill):
  * root  — early-saturating plateau (peak ~65 d); ORYZA ~2x larger.
  * stem  — vegetative rise to a maturity plateau (largest vegetative organ).
  * grain — ~0 until flowering, then the post-anthesis S-curve to a maturity plateau
            (a terminal accumulator; no decline).
  * leaf  — THE decisive difference: growth_rice plateaus (NO senescence) while
            ORYZA2000 peaks ~71 d then DECLINES ~57% (leaf senescence + remobilisation).
The leaf senescence is why driving the model with ORYZA biomass RAISES the leaf TF: the
growth-dilution sink mu=(dM/dt)/M goes NEGATIVE on a senescing leaf, so the -mu*C term
concentrates. BUT that rise is likely a PARTIAL ARTIFACT: `oryza_growth` models the loss as
leaf DEATH (dlv=drlv*wlv, carbon REMOVED from the plant), yet the PFAS ODE only sees the net
M(t) and conserves the leaf burden -- there is NO flux removing PFAS with the dead/fallen
leaf. A consistent litterfall term (-drlv*C) would cancel the death part of -mu*C, leaving
only the (always-diluting) growth term -G/M*C. So the concentration rise is defensible ONLY
if the senescing leaf stays attached and sheds mobile dry matter while retaining the
(immobile, bound) PFAS; under the biomass model's own "death = removed" bookkeeping it
OVER-states the leaf TF. Open issue (the M(t) shape is robust; its PFAS coupling is not).

Run:  python validation/biomass_shape_compare.py
      -> table + figure validation/figures/biomass_shape_compare.png
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
import growth_rice as gr            # noqa: E402
import oryza_growth as og           # noqa: E402

SEASON = 120.0
ORGANS = ["root", "stem", "leaf", "grain"]
COL = {"root": "#8c564b", "stem": "#2ca02c", "leaf": "#1f77b4", "grain": "#ff7f0e"}
FLOWER_D, MATURE_D = 66.0, 116.0    # IR72 anchors


def trajectories(t):
    G = {k: np.asarray(v) for k, v in gr.organ_biomass(t, SEASON).items()}
    O = {k: np.asarray(v) for k, v in og.organ_biomass_oryza(t, p=og.OryzaParams(season=SEASON)).items()}
    return G, O


def shape_features(t, M):
    """Per-organ (onset_day, peak_day, peak_g, final_g, peak->final decline %)."""
    out = {}
    for k in ORGANS:
        y = M[k] * 1000.0                       # g/hill
        ip = int(np.argmax(y))
        peak, fin = y[ip], y[-1]
        onset = t[np.argmax(y > 0.05 * peak)] if peak > 1e-9 else np.nan
        decline = 100.0 * (peak - fin) / peak if peak > 1e-9 else 0.0
        out[k] = (onset, t[ip], peak, fin, decline)
    return out


def main():
    t = np.linspace(0.0, SEASON, 241)
    G, O = trajectories(t)
    fG, fO = shape_features(t, G), shape_features(t, O)

    print(f"{'organ':6}{'driver':12}{'onset_d':>8}{'peak_d':>8}{'peak_g':>9}{'final_g':>9}{'declined%':>10}")
    for k in ORGANS:
        for nm, f in (("growth_rice", fG), ("ORYZA2000", fO)):
            on, pd, pk, fn, dc = f[k]
            print(f"{k:6}{nm:12}{on:>8.0f}{pd:>8.0f}{pk:>9.2f}{fn:>9.2f}{dc:>10.0f}")
    print("\n-> root/stem/grain shapes agree (rise->plateau; grain post-anthesis S-curve);")
    print("   leaf is the decisive split: growth_rice plateaus, ORYZA2000 senesces (~57% down).")

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), sharey=True)
    for ax, (M, title) in zip(axes, ((G, "growth_rice (ORYZA partition × logistic)"),
                                     (O, "ORYZA2000 (mechanistic carbon balance)"))):
        for k in ORGANS:
            ax.plot(t, M[k] * 1000.0, color=COL[k], lw=2.4, label=k)
        ax.plot(t, sum(M[k] for k in ORGANS) * 1000.0, "k--", lw=1.3, label="whole plant")
        for d, lab in ((FLOWER_D, "flowering"), (MATURE_D, "maturity")):
            ax.axvline(d, color="grey", ls=":", lw=1)
            ax.text(d, ax.get_ylim()[1] * 0.96, lab, rotation=90, va="top", fontsize=7, color="grey")
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("days after transplant")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8, loc="upper left")
    axes[0].set_ylabel("organ biomass [g/hill, dry]")
    fig.suptitle("Per-compartment biomass M_k(t): growth_rice vs ORYZA2000 (season 120 d)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = os.path.join(HERE, "figures", "biomass_shape_compare.png")
    fig.savefig(out, dpi=140)
    print(f"\n[written] {out}")


if __name__ == "__main__":
    main()
