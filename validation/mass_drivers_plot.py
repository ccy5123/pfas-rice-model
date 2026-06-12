#!/usr/bin/env python3
"""
Compartment-mass drivers M_k(t), dM/dt, growth-dilution mu(t)  (driver diagnostic)
==================================================================================

Visualises the per-compartment mass driver that simulate() feeds the ODE, for the
two built-in sources:
  (left)  measured_forcing=True  -> ORYZA IR72 biomass (growth_rice)         [DEFAULT]
  (right) measured_forcing=False -> illustrative logistic placeholders

Point of the figure: M_k is a TIME-VARYING growth curve, not a constant, and it
enters the model in three ways -- (a) the advective dilution denominator Q/M_k,
(b) the growth-dilution sink mu_k = (dM_k/dt)/M_k (the ONLY sink for the terminal
leaf & grain compartments; mu -> 0 at maturity => no steady state, hence a dynamic
model is required), and (c) the phloem flow Q_Phl ~ dM_grain/dt.  This replicates
exactly how model_api._default_drivers builds M and how PlantInputs/RHS derive
dM (np.gradient) and mu = dM/M (with the 1e-4 mass floor).

Run:  python validation/mass_drivers_plot.py  ->  validation/figures/mass_drivers.png
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
import growth_rice as gr                                          # noqa: E402
from pfas_rice_plant_module_4pool_surf import _logistic           # noqa: E402

SEASON = 120.0
t = np.linspace(0.0, SEASON, 481)
ORGANS = ["root", "stem", "leaf", "grain"]
COL = {"root": "saddlebrown", "stem": "goldenrod", "leaf": "forestgreen", "grain": "darkorange"}
FLOOR = 1e-4  # model_api._default_drivers floors M at 1e-4

# (1) ORYZA default (measured_forcing=True)
b = gr.organ_biomass(t, SEASON)
M_oryza = np.maximum(np.column_stack([b["root"], b["stem"], b["leaf"], b["grain"]]), FLOOR)

# (2) logistic placeholder (measured_forcing=False)
M_plac = np.maximum(np.column_stack([
    _logistic(t, 1e-3, 0.030, 0.10, SEASON * 0.17),
    _logistic(t, 1e-3, 0.040, 0.10, SEASON * 0.21),
    _logistic(t, 1e-3, 0.050, 0.12, SEASON * 0.25),
    _logistic(t, 1e-5, 0.025, 0.18, SEASON * 0.67)]), FLOOR)


def derived(M):
    dM = np.gradient(M, t, axis=0)          # exactly what PlantInputs does
    mu = dM / np.maximum(M, 1e-12)          # growth-dilution rate, the RHS sink
    return dM, mu


def main():
    sources = [("measured_forcing=True  —  ORYZA IR72  (DEFAULT)", M_oryza),
               ("measured_forcing=False  —  logistic placeholder", M_plac)]

    fig, ax = plt.subplots(3, 2, figsize=(12.5, 11.0), sharex=True)
    t_flower = 0.54 * SEASON

    for j, (title, M) in enumerate(sources):
        dM, mu = derived(M)
        for i, org in enumerate(ORGANS):
            ax[0, j].plot(t, M[:, i] * 1e3, color=COL[org], lw=2, label=org)
            ax[1, j].plot(t, dM[:, i] * 1e3, color=COL[org], lw=2, label=org)
            ax[2, j].plot(t, mu[:, i], color=COL[org], lw=2, label=org)
        ax[0, j].set_title(title, fontsize=11, fontweight="bold")
        for r in range(3):
            ax[r, j].axvline(t_flower, color="grey", ls=":", lw=1)
            ax[r, j].grid(alpha=0.25)
        ax[0, j].text(t_flower + 1.5, 0.97 * ax[0, j].get_ylim()[1], "flowering (DVS=1)",
                      fontsize=8, color="grey", va="top")

    ax[0, 0].set_ylabel("M_k(t)   [g / hill]")
    ax[1, 0].set_ylabel("dM_k/dt   [g / hill / day]")
    ax[2, 0].set_ylabel(r"$\mu_k=(dM_k/dt)/M_k$   [1/day]")
    for j in range(2):
        ax[2, j].set_xlabel("day after transplant")
        ax[2, j].axhline(0.0, color="k", lw=0.8)
        ax[2, j].set_ylim(-0.03, 0.82)
    ax[0, 0].legend(loc="upper left", fontsize=9, ncol=2, framealpha=0.9)

    for j in range(2):
        ax[2, j].annotate("growth-dilution sink → 0\n(leaf & grain are terminal\naccumulators; no steady state)",
                          xy=(113, 0.01), xytext=(55, 0.46), fontsize=8.5, color="firebrick",
                          arrowprops=dict(arrowstyle="->", color="firebrick", lw=1.0))

    fig.suptitle("Compartment-mass drivers in the PFAS–rice model:  M is a time-varying growth curve, not a constant",
                 fontsize=12.5, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.975])
    out = os.path.join(HERE, "figures", "mass_drivers.png")
    fig.savefig(out, dpi=130)
    print("[written]", out)

    # compact numeric summary
    print("\nfinal masses [g/hill] @ day 120:")
    for name, M in sources:
        fin = M[-1] * 1e3
        print(f"  {name.split('—')[1].strip():28s} "
              f"root={fin[0]:6.2f} stem={fin[1]:6.2f} leaf={fin[2]:6.2f} grain={fin[3]:6.2f}")
    print("\npeak growth-dilution mu [1/day] (ORYZA), and mu(120):")
    _, mu_o = derived(M_oryza)
    for i, org in enumerate(ORGANS):
        k = np.argmax(mu_o[:, i])
        print(f"  {org:5s}: mu_max={mu_o[k, i]:.3f} @ day {t[k]:.0f}  ->  mu(120)={mu_o[-1, i]:.4f}")


if __name__ == "__main__":
    main()
