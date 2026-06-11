#!/usr/bin/env python3
"""
HYDRUS-1D -> plant coupling: full Method-A soil->plant run
==========================================================

Drives the 4-compartment plant ODE with a REAL HYDRUS-1D pore-water trajectory
C_w^o(t) (built in ``src/soil_hydrus.py`` from the compiled HYDRUS engine via
phydrus) and contrasts it with the constant-Cwo baseline.

Both runs share the SAME measured growth M(t) (ORYZA) and transpiration Q_TP(t)
(forcing_rice), and the HYDRUS pore water is normalised to the same season-mean
exposure (1.0) as the baseline -- so the only difference is the realistic soil
temporal structure: weakly-sorbed short chains leach during flooding (pore water
collapses, so the late-filling grain sees little), strongly-sorbed long chains
stay buffered (BAF ~ unchanged).

Run:  python validation/hydrus_coupled_run.py
Saves: validation/figures/hydrus_coupled.png, validation/hydrus_coupled.csv
"""
import csv
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))

from pfas_rice_plant_module_4pool_surf import (  # noqa: E402
    Environment, Compound, Compartment, RiceUptakeModel, PlantInputs,
    _logistic, ROOT as I_ROOT, STEM, LEAF, FRUIT)
import soil_hydrus as sh  # noqa: E402
import growth_rice as gr  # noqa: E402
import forcing_rice as fr  # noqa: E402

PAR = json.load(open(os.path.join(ROOT, "params", "parameters.json")))
CARR = PAR["carrier_MichaelisMenten"]
COMP = PAR["tissue_composition_recommended"]
CONG = {c["name"]: c for c in PAR["congeners"]}

CONGENERS = ["PFBA", "PFHxA", "PFOA", "PFNA", "PFOS", "PFDoDA"]
SEASON = 120.0
N_T = 241


def compartments():
    g = COMP
    return [Compartment("root", g["root"]["theta_fw"], g["root"]["f_prot"], g["root"]["f_PL"], g["root"]["f_cw"]),
            Compartment("stem", g["stem"]["theta_fw"], g["stem"]["f_prot"], g["stem"]["f_PL"], g["stem"]["f_cw"]),
            Compartment("leaf", g["leaf"]["theta_fw"], g["leaf"]["f_prot"], g["leaf"]["f_PL"], g["leaf"]["f_cw"], S=20.0),
            Compartment("grain", g["grain_brown"]["theta_fw"], g["grain_brown"]["f_prot"],
                        g["grain_brown"]["f_PL"], g["grain_brown"]["f_cw"], S=2.0)]


def compound(name):
    c = CONG[name]
    return Compound(name=name, K_prot=c["K_prot_Lkg"], K_PL=c["K_PL_Lkg"],
                    K_cw=c["K_cw_wholecw_Lkg"]["root"],
                    kappa_d=c.get("kappa_d_W2fit") or 2.0,
                    Vmax_in=CARR["Vmax_in"], Km_in=CARR["Km_in"],
                    Vmax_out=CARR["Vmax_out"], Km_out=CARR["Km_out"],
                    L_Ph=c.get("L_Ph_W2fit") or 0.01,
                    f_xy=c["f_xy_recommended"])


def solve_baf(inputs, name):
    """Final root/straw/grain BAF (relative to season-mean pore water = 1.0)."""
    model = RiceUptakeModel(env=Environment(), cmpd=compound(name),
                            comps=compartments(), inputs=inputs)
    sol = model.solve(inputs.t)
    C = sol.y[:, -1]
    Mf = inputs.M_(inputs.t[-1])
    straw = (C[STEM] * Mf[STEM] + C[LEAF] * Mf[LEAF]) / (Mf[STEM] + Mf[LEAF])
    cwo_ref = float(np.mean(inputs.Cwo))
    return dict(root=C[I_ROOT] / cwo_ref, straw=straw / cwo_ref, grain=C[FRUIT] / cwo_ref)


def baseline_inputs():
    t = np.linspace(0.0, SEASON, N_T)
    b = gr.organ_biomass(t, SEASON)
    M = np.maximum(np.column_stack([b["root"], b["stem"], b["leaf"], b["grain"]]), 1e-4)
    return PlantInputs(t=t, Cwo=np.full_like(t, 1.0), Qtp=fr.Q_TP(t, SEASON), M=M)


def main():
    if not sh.hydrus_available():
        print("HYDRUS executable not built (see src/soil_hydrus.py). Aborting.")
        return 1

    base = baseline_inputs()
    rows, traj = [], {}
    print(f"{'cong':8}{'Kd':>9}  {'root b/h':>16}{'straw b/h':>16}{'grain b/h':>16}")
    for name in CONGENERS:
        c = CONG[name]
        inp_h, res = sh.inputs_from_hydrus(c["n_C"], c["group"], season=SEASON,
                                           Cwo_ref=1.0, n_t=N_T)
        traj[name] = (inp_h.t, inp_h.Cwo, res.Kd)
        b = solve_baf(base, name)
        h = solve_baf(inp_h, name)
        rows.append(dict(congener=name, n_C=c["n_C"], group=c["group"], Kd=res.Kd,
                         **{f"{k}_baseline": b[k] for k in b},
                         **{f"{k}_hydrus": h[k] for k in h}))
        print(f"{name:8}{res.Kd:>9.3f}  "
              f"{b['root']:>7.2f}/{h['root']:<7.2f}"
              f"{b['straw']:>7.2f}/{h['straw']:<7.2f}"
              f"{b['grain']:>7.2f}/{h['grain']:<7.2f}")

    out_csv = os.path.join(HERE, "hydrus_coupled.csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {out_csv}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.3))
        for name in CONGENERS:
            t, cwo, Kd = traj[name]
            ax1.plot(t, cwo, lw=2, label=f"{name} (Kd={Kd:.2g})")
        ax1.set_xlabel("time [day]"); ax1.set_ylabel("pore-water C_w^o(t)  [norm., mean=1]")
        ax1.set_title("HYDRUS-1D pore water (clean-water flooding → drainage)")
        ax1.legend(fontsize=7); ax1.axvline(90, ls=":", c="grey")

        x = np.arange(len(CONGENERS)); w = 0.35
        gb = [r["grain_baseline"] for r in rows]
        gh = [r["grain_hydrus"] for r in rows]
        ax2.bar(x - w / 2, gb, w, label="grain BAF (constant Cwo)")
        ax2.bar(x + w / 2, gh, w, label="grain BAF (HYDRUS Cwo)")
        ax2.set_xticks(x); ax2.set_xticklabels(CONGENERS, rotation=45, ha="right")
        ax2.set_ylabel("grain BAF [L/kg]"); ax2.set_title("late-filling grain: effect of leaching")
        ax2.legend(fontsize=8)
        fig.tight_layout()
        figpath = os.path.join(HERE, "figures", "hydrus_coupled.png")
        os.makedirs(os.path.dirname(figpath), exist_ok=True)
        fig.savefig(figpath, dpi=130)
        print(f"wrote {figpath}")
    except Exception as e:
        print(f"(plot skipped: {e})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
