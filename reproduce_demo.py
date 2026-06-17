#!/usr/bin/env python3
# =============================================================================
# reproduce_demo.py — self-contained reproduction entry point
# -----------------------------------------------------------------------------
# Loads params/parameters.json + the basis-A plant module (src/) and reproduces
# the Yamazaki2023 root/straw/grain BAF via the full 4-compartment ODE, using
# the S6 W2 transport fit. Reproduces the Gap4 validation (log10 RMSE ≈ 0.029).
#
#   python reproduce_demo.py            # uses f_xy_W2fit (reproduces obs)
#   python reproduce_demo.py --rec      # uses f_xy_recommended (monotone; see note)
#
# NOTE on the two f_xy:
#   * f_xy_W2fit reproduces Yamazaki but rises spuriously for C10+ (single-straw-
#     compartment entanglement, H7 §7.2) -> NOT the physical TSCF.
#   * f_xy_recommended is the monotone physical TSCF (theory + cross-field TF).
#     With --rec the predicted straw does NOT match obs (over-predicted for short
#     chains, under-predicted for long): the single mass-weighted straw compartment
#     cannot host the observed long-chain stem accumulation gradient, and the W2 fit
#     was absorbing that. Refining the stem (multi-height) is the route to using the
#     physical f_xy directly. See docs/DELIVERABLE_GAP_B_fxy.md.
# =============================================================================
import json, csv, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "src"))
from pfas_rice_plant_module_4pool_surf import (   # noqa: E402
    Environment, Compound, Compartment, RiceUptakeModel, PlantInputs,
    _logistic, ROOT, STEM, LEAF, FRUIT)

USE_REC = "--rec" in sys.argv

with open(os.path.join(HERE, "params", "parameters.json")) as f:
    PAR = json.load(f)

# observed Yamazaki BAF
obs = {}
with open(os.path.join(HERE, "data_obs", "obs_baf_Yamazaki.csv"), newline="") as f:
    for r in csv.DictReader(f):
        obs.setdefault(r["compound"], {})[r["tissue"]] = float(r["baf"])

# tissue compartments (recommended composition; module needs a single K_cw per
# compound, so we use the per-congener ROOT whole-cw — organ spread is ~15%)
comp = PAR["tissue_composition_recommended"]
def compartments():
    return [
        Compartment("root",  comp["root"]["theta_fw"],  comp["root"]["f_prot"],  comp["root"]["f_PL"],  comp["root"]["f_cw"]),
        Compartment("stem",  comp["stem"]["theta_fw"],  comp["stem"]["f_prot"],  comp["stem"]["f_PL"],  comp["stem"]["f_cw"]),
        Compartment("leaf",  comp["leaf"]["theta_fw"],  comp["leaf"]["f_prot"],  comp["leaf"]["f_PL"],  comp["leaf"]["f_cw"], S=20.0),
        Compartment("grain", comp["grain_brown"]["theta_fw"], comp["grain_brown"]["f_prot"], comp["grain_brown"]["f_PL"], comp["grain_brown"]["f_cw"], S=2.0),
    ]

# demo drivers (placeholders; Cwo=1 -> C == BAF). Replace with HYDRUS/growth output.
t = np.linspace(0.0, 120.0, 481)
inputs = PlantInputs(
    t=t, Cwo=np.full_like(t, 1.0),
    Qtp=0.05 + 0.35 * np.exp(-((t - 75.0) ** 2) / (2 * 25.0 ** 2)),
    M=np.column_stack([_logistic(t, 1e-3, 0.030, 0.10, 20.0), _logistic(t, 1e-3, 0.040, 0.10, 25.0),
                       _logistic(t, 1e-3, 0.050, 0.12, 30.0), _logistic(t, 1e-5, 0.025, 0.18, 80.0)]))
carr = PAR["carrier_MichaelisMenten"]; env = Environment()

print(f"f_xy source: {'RECOMMENDED (monotone, physical)' if USE_REC else 'W2 fit (reproduces obs)'}")
print(f"{'PFAS':8}{'nC':>3}{'f_xy':>8} | {'root p/o':>16}{'straw p/o':>16}{'grain p/o':>16}")
errs = []
for c in PAR["congeners"]:
    p = c["name"]
    if p not in obs or c["f_xy_W2fit"] is None:
        continue                                   # need obs + a transport fit
    f_xy = c["f_xy_recommended"] if USE_REC else c["f_xy_W2fit"]
    cmpd = Compound(name=p, K_prot=c["K_prot_Lkg"], K_PL=c["K_PL_Lkg"],
                    K_cw=c["K_cw_wholecw_Lkg"]["root"],
                    kappa_d=c["kappa_d_W2fit"], Vmax_in=carr["Vmax_in"], Km_in=carr["Km_in"],
                    Vmax_out=carr["Vmax_out"], Km_out=carr["Km_out"],
                    L_Ph=c["L_Ph_W2fit"], f_xy=f_xy)
    sol = RiceUptakeModel(env=env, cmpd=cmpd, comps=compartments(), inputs=inputs).solve(t)
    C = sol.y[:, -1]; Mf = inputs.M_(t[-1])
    pr = {"root": C[ROOT],
          "straw": (C[STEM]*Mf[STEM]+C[LEAF]*Mf[LEAF])/(Mf[STEM]+Mf[LEAF]),
          "grain": C[FRUIT]}
    o = obs[p]
    for k in ("root", "straw", "grain"):
        if k in o:
            errs.append((np.log10(max(pr[k], 1e-6)) - np.log10(o[k])) ** 2)
    print(f"{p:8}{c['n_C']:>3}{f_xy:>8.4f} | "
          f"{pr['root']:>7.2f}/{o.get('root', float('nan')):<7.2f}"
          f"{pr['straw']:>7.2f}/{o.get('straw', float('nan')):<7.2f}"
          f"{pr['grain']:>7.2f}/{o.get('grain', float('nan')):<7.2f}")
print(f"\nlog10 RMSE (pred vs obs) = {np.sqrt(np.mean(errs)):.3f}"
      f"   {'(monotone f_xy does NOT reproduce obs via this structure — straw mismatch, single-compartment limit; see note)' if USE_REC else '(W2 fit reproduces; PFDoDA near-MQL outlier)'}")
if not USE_REC:
    print("NOTE: 0.029 is IN-SAMPLE saturated reproduction (3 fitted transport params per "
          "congener vs 3 obs), NOT predictive validation. The genuine a-priori predictive "
          "error (monotone f_xy: `reproduce_demo.py --rec`) is log10 RMSE ~0.84, straw 6-40x "
          "off. See validation/apriori_prediction.py and the sci-adk rigor review "
          "(sci_adk_review/FINDINGS.md: hyp-yamazaki REFUTED).")
