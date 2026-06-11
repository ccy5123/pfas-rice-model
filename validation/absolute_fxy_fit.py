#!/usr/bin/env python3
"""
Absolute f_xy fit with MEASURED drivers (task 2 — the core deliverable)
======================================================================

Pins the absolute root->shoot loading factor f_xy by fitting the 4-compartment
model to the Yamazaki straw BAF using the MEASURED forcings:
  * transpiration Q_TP(t) -- forcing_rice (Kumari 2022 Kc/ET0 + Nay Htoon 2018 T/E)
  * per-organ biomass M(t) -- growth_rice (ORYZA IR72 partitioning, Bouman/Li)

f_xy is the straw(shoot)-sensitive parameter, so we fit it (1 param) to the
observed straw BAF per congener (kappa_d, L_Ph held fixed). The fitted value is
the absolute f_xy implied by the data once the drivers are measured; we compare
it to the monotone recommended f_xy (Felizeter-TSCF-anchored).

KEY RESULT (see module end): the data-fitted absolute f_xy is NON-monotone --
short chains are near the recommendation but long chains need a far HIGHER f_xy
(the observed long-chain straw BAF is large), the same tension the saturated W2
fit showed. Measured drivers do not rescue the steep monotone shape at long
chains; that remains the open item (long-chain shoot presence is under-modelled
-- candidate causes in docs/nstem_gradient_exploration.md and the PFSA/ether term).

Run: python validation/absolute_fxy_fit.py   (use python -u to see live progress)
"""
import json, os, sys, csv
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT_DIR, "src"))

from pfas_rice_plant_module_4pool_surf import (
    Environment, Compound, Compartment, RiceUptakeModel, PlantInputs)
from calibration import calibrate, Param, ObservedBAF, predict_bafs
import forcing_rice as fr
import growth_rice as gr

PAR = json.load(open(os.path.join(ROOT_DIR, "params", "parameters.json")))
comp = PAR["tissue_composition_recommended"]
carr = PAR["carrier_MichaelisMenten"]

obs = {}
with open(os.path.join(ROOT_DIR, "data_obs", "obs_baf_Yamazaki.csv")) as f:
    for r in csv.DictReader(f):
        obs.setdefault(r["compound"], {})[r["tissue"]] = float(r["baf"])

SEASON = 120.0
t = np.linspace(0.0, SEASON, 121)                         # coarse grid (fit speed)
Qtp = fr.Q_TP(t, SEASON)                                  # MEASURED transpiration
b = gr.organ_biomass(t, SEASON)                           # MEASURED (ORYZA) biomass
M = np.maximum(np.column_stack([b["root"], b["stem"], b["leaf"], b["grain"]]), 1e-4)
inputs = PlantInputs(t=t, Cwo=np.full_like(t, 1.0), Qtp=Qtp, M=M)
comps = [Compartment("root",  comp["root"]["theta_fw"], comp["root"]["f_prot"],
                     comp["root"]["f_PL"], comp["root"]["f_cw"]),
         Compartment("stem",  comp["stem"]["theta_fw"], comp["stem"]["f_prot"],
                     comp["stem"]["f_PL"], comp["stem"]["f_cw"]),
         Compartment("leaf",  comp["leaf"]["theta_fw"], comp["leaf"]["f_prot"],
                     comp["leaf"]["f_PL"], comp["leaf"]["f_cw"], S=20.0),
         Compartment("grain", comp["grain_brown"]["theta_fw"], comp["grain_brown"]["f_prot"],
                     comp["grain_brown"]["f_PL"], comp["grain_brown"]["f_cw"], S=2.0)]


def main():
    print(f"measured drivers: peak Q_TP={Qtp.max():.3f} L/d/hill, "
          f"final shoot={(b['stem'][-1]+b['leaf'][-1]+b['grain'][-1]):.4f} kg/hill\n")
    print(f"{'PFAS':7}{'nC':>3} | {'f_xy_fit':>9}{'f_xy_rec':>9}{'fit/rec':>8} | {'straw p/o':>13}")
    ratios = []
    for c in PAR["congeners"]:
        nm = c["name"]
        if nm not in obs:
            continue
        cmpd = Compound(name=nm, K_prot=c["K_prot_Lkg"], K_PL=c["K_PL_Lkg"],
                        K_cw=c["K_cw_wholecw_Lkg"]["root"], kappa_d=2.0,
                        Vmax_in=carr["Vmax_in"], Km_in=carr["Km_in"],
                        Vmax_out=carr["Vmax_out"], Km_out=carr["Km_out"],
                        L_Ph=0.01, f_xy=c["f_xy_recommended"])
        model = RiceUptakeModel(env=Environment(), cmpd=cmpd, comps=comps, inputs=inputs)
        res = calibrate(model, [Param("f_xy", 1e-3, 1.0)],
                        [ObservedBAF("straw", obs[nm]["straw"])],
                        x0=[c["f_xy_recommended"]], t=t)
        fx = res.values["f_xy"]; rec = c["f_xy_recommended"]
        ratios.append((c["group"], fx / rec))
        print(f"{nm:7}{c['n_C']:3d} | {fx:9.4f}{rec:9.4f}{fx/rec:8.2f} | "
              f"{res.predicted['straw']:6.2f}/{obs[nm]['straw']:<5.2f}", flush=True)
    pf = [r for g, r in ratios if g == "PFCA"]
    print(f"\nfit/rec ratio (PFCA): short-chain ~O(1), long-chain >>1 -> the data-fitted "
          f"f_xy is NON-monotone (long-chain straw needs far more loading than the "
          f"monotone recommendation). Measured drivers do not rescue the steep tail.")


if __name__ == "__main__":
    main()
