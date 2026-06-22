#!/usr/bin/env python3
"""
Root:shoot / biomass-partitioning sensitivity of the W2 BAF reproduction
========================================================================

Audits how the (largely unstated) BIOMASS-PARTITIONING assumption controls the
Yamazaki root/straw/grain BAF reproduction. It contrasts three organ-biomass
drivers, holding the W2 transport fit (params/parameters.json) fixed:

  1. reproduce_demo placeholder logistic  (root:shoot ~0.30, HI ~0.07)  <- the
       biomass the W2 f_xy/L_Ph fit was actually calibrated against.
  2. growth_rice (ORYZA IR72, app/simulate default)  (root:shoot ~0.035, HI ~0.51)
  3. growth_rice with the root trajectory rescaled to a LITERATURE maturity
       root:shoot (0.10 / 0.15) while keeping the realistic HI.

Key finding (see docs/biomass_partitioning_rootshoot.md):
  the log10 RMSE 0.029 "reproduction" is only attained with the placeholder, whose
  root:shoot (0.30) and especially HI (0.07) are NON-PHYSICAL (modern lowland rice
  HI ~0.45-0.55, root:shoot ~0.08-0.13). With any realistic HI the SAME W2 params
  give RMSE ~0.25-0.31 regardless of root:shoot -> the transport fit is entangled
  with a non-physical biomass and must be RE-FITTED on a literature-consistent one.

Literature root:shoot anchor (maturity, lowland rice): ~0.08-0.13 (root ~7-12% of
total plant), declining from ~0.2 at seedling; high-yield cultivars lower (grain-
fill dilution). Sources: Frontiers Plant Sci. 2021 (10.3389/fpls.2021.713814,
japonica paddy, tillering->maturity SDW/RDW); Yoshida 1981 "Fundamentals of Rice
Crop Science" (IRRI). Aboveground split validated by Ntanos & Koutroubas 2002
(Field Crops Res. 74:93-101, 10.1016/S0378-4290(01)00203-9; HI 0.47-0.61) and
Amanullah & Inamullah 2016 (Rice Sci. 23(2):78-87; maturity panicle 43-58% /
culm 24-33% / leaf 18-24%).
"""
import json, csv, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT_DIR, "src"))
from pfas_rice_plant_module_4pool_surf import (   # noqa: E402
    Environment, Compound, Compartment, RiceUptakeModel, PlantInputs,
    _logistic, ROOT, STEM, LEAF, FRUIT)
import growth_rice as gr                            # noqa: E402

PAR = json.load(open(os.path.join(ROOT_DIR, "params", "parameters.json")))
OBS = {}
with open(os.path.join(ROOT_DIR, "data_obs", "obs_baf_Yamazaki.csv")) as f:
    for r in csv.DictReader(f):
        OBS.setdefault(r["compound"], {})[r["tissue"]] = float(r["baf"])
COMP = PAR["tissue_composition_recommended"]
CARR = PAR["carrier_MichaelisMenten"]
ENV = Environment()
T = np.linspace(0.0, 120.0, 481)
QTP = 0.05 + 0.35 * np.exp(-((T - 75.0) ** 2) / (2 * 25.0 ** 2))


def _compartments():
    g = COMP
    return [Compartment("root",  g["root"]["theta_fw"],  g["root"]["f_prot"],  g["root"]["f_PL"],  g["root"]["f_cw"]),
            Compartment("stem",  g["stem"]["theta_fw"],  g["stem"]["f_prot"],  g["stem"]["f_PL"],  g["stem"]["f_cw"]),
            Compartment("leaf",  g["leaf"]["theta_fw"],  g["leaf"]["f_prot"],  g["leaf"]["f_PL"],  g["leaf"]["f_cw"], S=20.0),
            Compartment("grain", g["grain_brown"]["theta_fw"], g["grain_brown"]["f_prot"], g["grain_brown"]["f_PL"], g["grain_brown"]["f_cw"], S=2.0)]


def placeholder_M():
    return np.column_stack([_logistic(T, 1e-3, 0.030, 0.10, 20.0), _logistic(T, 1e-3, 0.040, 0.10, 25.0),
                            _logistic(T, 1e-3, 0.050, 0.12, 30.0), _logistic(T, 1e-5, 0.025, 0.18, 80.0)])


def growth_rice_M(root_shoot=None):
    b = gr.organ_biomass(T, 120.0, root_shoot=root_shoot)
    return np.maximum(np.column_stack([b["root"], b["stem"], b["leaf"], b["grain"]]), 1e-6)


def partition_stats(M):
    f = M[-1]; tot = f.sum(); shoot = f[STEM] + f[LEAF] + f[FRUIT]
    return dict(root_pct=100 * f[ROOT] / tot, root_shoot=f[ROOT] / shoot, HI=f[FRUIT] / tot)


def w2_rmse(M):
    """log10 RMSE of root/straw/grain BAF vs Yamazaki, W2 fit fixed, on biomass M."""
    inp = PlantInputs(t=T, Cwo=np.full_like(T, 1.0), Qtp=QTP, M=M)
    errs = []
    for c in PAR["congeners"]:
        p = c["name"]
        if p not in OBS or c["f_xy_W2fit"] is None:
            continue
        cmpd = Compound(name=p, K_prot=c["K_prot_Lkg"], K_PL=c["K_PL_Lkg"], K_cw=c["K_cw_wholecw_Lkg"]["root"],
                        kappa_d=c["kappa_d_W2fit"], Vmax_in=CARR["Vmax_in"], Km_in=CARR["Km_in"],
                        Vmax_out=CARR["Vmax_out"], Km_out=CARR["Km_out"], L_Ph=c["L_Ph_W2fit"], f_xy=c["f_xy_W2fit"])
        sol = RiceUptakeModel(env=ENV, cmpd=cmpd, comps=_compartments(), inputs=inp).solve(T)
        C = sol.y[:, -1]; Mf = inp.M_(T[-1])
        pr = {"root": C[ROOT],
              "straw": (C[STEM] * Mf[STEM] + C[LEAF] * Mf[LEAF]) / (Mf[STEM] + Mf[LEAF]),
              "grain": C[FRUIT]}
        for k in ("root", "straw", "grain"):
            if k in OBS[p]:
                errs.append((np.log10(max(pr[k], 1e-6)) - np.log10(OBS[p][k])) ** 2)
    return float(np.sqrt(np.mean(errs)))


def main():
    def growth_rice_M_B(target_rs):
        b = gr.organ_biomass(T, 120.0, target_root_shoot=target_rs)   # method B (FRT)
        return np.maximum(np.column_stack([b["root"], b["stem"], b["leaf"], b["grain"]]), 1e-6)

    rows = [("reproduce placeholder", placeholder_M()),
            ("growth_rice (default ~0.035)", growth_rice_M(None)),
            ("growth_rice C rescale R/S=0.10", growth_rice_M(0.10)),
            ("growth_rice B FRT R/S=0.10", growth_rice_M_B(0.10)),
            ("growth_rice R/S=0.15 (C)", growth_rice_M(0.15)),
            ("growth_rice R/S=0.30 (~placeholder, C)", growth_rice_M(0.30))]
    print(f"{'biomass driver':36s} {'root%':>6s} {'R/S':>6s} {'HI':>6s} {'W2 RMSE':>8s}")
    print("-" * 66)
    for lab, M in rows:
        s = partition_stats(M)
        print(f"{lab:36s} {s['root_pct']:6.1f} {s['root_shoot']:6.3f} {s['HI']:6.3f} {w2_rmse(M):8.3f}")
    print("\nLiterature maturity root:shoot ~0.08-0.13; HI ~0.45-0.55.")
    print("=> RMSE 0.029 needs the non-physical placeholder (HI 0.07); realistic HI")
    print("   gives ~0.25-0.31 for ANY root:shoot -> the W2 fit must be re-fitted on a")
    print("   literature-consistent biomass (root:shoot ~0.10, HI ~0.5).")


if __name__ == "__main__":
    main()
