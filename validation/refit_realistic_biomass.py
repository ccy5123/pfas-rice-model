#!/usr/bin/env python3
"""
Re-fit the Tier-1 transport parameters on a LITERATURE-CONSISTENT biomass
=========================================================================

`validation/root_shoot_biomass_sensitivity.py` showed the W2 transport fit
(`params/parameters.json`) only reproduces Yamazaki on the non-physical
`reproduce_demo` placeholder biomass (root:shoot 0.30, HI 0.07). This script
RE-FITS, per congener, the three Tier-1 transport parameters

    f_xy   (root->xylem loading, TSCF analog)
    L_Ph   (phloem loading -> grain)
    kappa_d(lumped root membrane conductance -> root level)

to the SAME observed Yamazaki root/straw/grain BAF but on a literature-consistent
biomass: `growth_rice` with HI ~0.5 and the field root:shoot anchor (~0.10; Japanese
flooded paddy, see docs/biomass_partitioning_rootshoot.md). Drivers are otherwise
identical to reproduce_demo (Cwo=1, the same transpiration shape), so the ONLY change
from the W2 fit is the biomass -> the re-fitted parameters are the W2 fit "translated"
onto realistic biomass.

OVERRIDE-ONLY: params/parameters.json is NOT modified. The re-fit is written to
params/refit_realistic_biomass.csv for inspection / a future promotion decision.

What to look for
----------------
* Whether the saturated 3-param/3-obs fit recovers a low RMSE on realistic biomass
  (it should: the problem is determined) -> realistic biomass is fittable.
* How f_xy shifts vs the W2 fit, and whether the new f_xy is more physical
  (monotone-ish decline with chain length) than the W2 fit (which rose for C10+).
"""
import json, csv, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT_DIR, "src"))
from pfas_rice_plant_module_4pool_surf import (   # noqa: E402
    Environment, Compound, Compartment, RiceUptakeModel, PlantInputs)
import growth_rice as gr                            # noqa: E402
from calibration import Param, ObservedBAF, calibrate   # noqa: E402

ROOT_SHOOT = 0.10                                    # field anchor (Japanese paddy)
PAR = json.load(open(os.path.join(ROOT_DIR, "params", "parameters.json")))
CARR = PAR["carrier_MichaelisMenten"]
COMP = PAR["tissue_composition_recommended"]
OBS = {}
with open(os.path.join(ROOT_DIR, "data_obs", "obs_baf_Yamazaki.csv")) as f:
    for r in csv.DictReader(f):
        OBS.setdefault(r["compound"], {})[r["tissue"]] = float(r["baf"])

T = np.linspace(0.0, 120.0, 241)
QTP = 0.05 + 0.35 * np.exp(-((T - 75.0) ** 2) / (2 * 25.0 ** 2))   # same shape as reproduce_demo
_b = gr.organ_biomass(T, 120.0, root_shoot=ROOT_SHOOT)
M = np.maximum(np.column_stack([_b["root"], _b["stem"], _b["leaf"], _b["grain"]]), 1e-6)
INPUTS = PlantInputs(t=T, Cwo=np.full_like(T, 1.0), Qtp=QTP, M=M)
ENV = Environment()


def _compartments():
    g = COMP
    return [Compartment("root",  g["root"]["theta_fw"],  g["root"]["f_prot"],  g["root"]["f_PL"],  g["root"]["f_cw"]),
            Compartment("stem",  g["stem"]["theta_fw"],  g["stem"]["f_prot"],  g["stem"]["f_PL"],  g["stem"]["f_cw"]),
            Compartment("leaf",  g["leaf"]["theta_fw"],  g["leaf"]["f_prot"],  g["leaf"]["f_PL"],  g["leaf"]["f_cw"], S=20.0),
            Compartment("grain", g["grain_brown"]["theta_fw"], g["grain_brown"]["f_prot"], g["grain_brown"]["f_PL"], g["grain_brown"]["f_cw"], S=2.0)]


def refit_congener(c):
    """Return (new params dict, predicted, observed, rmse) for one congener record."""
    cmpd = Compound(name=c["name"], K_prot=c["K_prot_Lkg"], K_PL=c["K_PL_Lkg"],
                    K_cw=c["K_cw_wholecw_Lkg"]["root"], kappa_d=c["kappa_d_W2fit"],
                    Vmax_in=CARR["Vmax_in"], Km_in=CARR["Km_in"],
                    Vmax_out=CARR["Vmax_out"], Km_out=CARR["Km_out"],
                    L_Ph=c["L_Ph_W2fit"], f_xy=c["f_xy_W2fit"])
    model = RiceUptakeModel(env=ENV, cmpd=cmpd, comps=_compartments(), inputs=INPUTS)
    obs = [ObservedBAF(tis, OBS[c["name"]][tis]) for tis in ("root", "straw", "grain")
           if tis in OBS[c["name"]]]
    params = [Param("f_xy", 1e-3, 1.0), Param("L_Ph", 1e-4, 1.0), Param("kappa_d", 1e-3, 20.0)]
    # start from the W2 values, clipped into the (wider) bounds
    x0 = [min(max(v, p.low), p.high) for v, p in
          zip([c["f_xy_W2fit"], c["L_Ph_W2fit"], c["kappa_d_W2fit"]], params)]
    # local trf from the clipped W2 start; the 3-param/3-obs problem is determined,
    # so a local solve is enough (and fast/deterministic). A residual that stays high
    # is itself informative (that BAF is structurally unreachable on this biomass).
    res = calibrate(model, params, obs, x0=x0)
    errs = [(np.log10(max(res.predicted[o.tissue], 1e-9)) - np.log10(o.value)) ** 2 for o in obs]
    return res.values, res.predicted, {o.tissue: o.value for o in obs}, float(np.sqrt(np.mean(errs)))


def main():
    print(f"Re-fit on literature biomass (root:shoot={ROOT_SHOOT}, HI~{M[-1,3]/M[-1,1:].sum():.2f})\n")
    print(f"{'cong':8s}{'nC':>3s}{'grp':>5s} | {'f_xy old':>9s}{'f_xy new':>9s} | "
          f"{'L_Ph new':>9s}{'kapd new':>9s} | {'RMSE':>6s}")
    rows, all_err = [], []
    for c in PAR["congeners"]:
        p = c["name"]
        if p not in OBS or c.get("f_xy_W2fit") is None:
            continue
        vals, pred, obs, rmse = refit_congener(c)
        all_err.append(rmse)
        rows.append(dict(name=p, n_C=c["n_C"], group=c["group"], f_xy_old=c["f_xy_W2fit"],
                         f_xy_new=vals["f_xy"], L_Ph_new=vals["L_Ph"], kappa_d_new=vals["kappa_d"],
                         rmse=rmse))
        print(f"{p:8s}{c['n_C']:>3d}{c['group']:>5s} | {c['f_xy_W2fit']:>9.4f}{vals['f_xy']:>9.4f} | "
              f"{vals['L_Ph']:>9.4f}{vals['kappa_d']:>9.3f} | {rmse:>6.3f}")
    overall = float(np.sqrt(np.mean([e ** 2 for e in all_err])))
    print(f"\noverall per-congener log10 RMSE = {overall:.3f}  "
          f"(saturated 3-param/3-obs -> reproduction, NOT predictive validation)")

    # monotonicity of the re-fitted f_xy across the PFCA series
    pfca = [r for r in rows if r["group"] == "PFCA"]
    pfca.sort(key=lambda r: r["n_C"])
    drops = sum(1 for a, b in zip(pfca, pfca[1:]) if b["f_xy_new"] <= a["f_xy_new"] + 1e-9)
    print(f"PFCA f_xy_new monotone-decreasing over {drops}/{max(len(pfca)-1,1)} chain-length steps "
          f"(W2 fit rose for C10+; physical TSCF should decline).")

    out = os.path.join(ROOT_DIR, "params", "refit_realistic_biomass.csv")
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["name", "n_C", "group", "f_xy_old",
                                          "f_xy_new", "L_Ph_new", "kappa_d_new", "rmse"])
        w.writeheader()
        for r in rows:
            w.writerow({k: (round(v, 5) if isinstance(v, float) else v) for k, v in r.items()})
    print(f"\nwrote {os.path.relpath(out, ROOT_DIR)} (OVERRIDE-only; parameters.json unchanged)")


if __name__ == "__main__":
    main()
