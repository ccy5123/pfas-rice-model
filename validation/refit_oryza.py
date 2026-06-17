#!/usr/bin/env python3
# =============================================================================
# validation/refit_oryza.py
# -----------------------------------------------------------------------------
# RE-FIT the per-congener transport parameters (f_xy, L_Ph, kappa_d) on the
# MECHANISTIC ORYZA2000 biomass (now the default driver) so the model reproduces
# the Yamazaki tissue BAFs under the new default -- the W2 fit in
# params/parameters.json was tuned on the placeholder/growth_rice driver and no
# longer reproduces once biomass="oryza".
#
# This is the ORYZA analog of the per-congener W2 fit: 3 transport params fit to
# the 3 tissue obs per congener (root<-kappa_d, straw<-f_xy, grain<-L_Ph), solved
# by sequential 1-D brentq on the monotone tissue responses (a few coupling
# passes). It is SATURATED per congener (DOF 0) -> reproduces by construction;
# the *constrained* DOF>0 goodness-of-fit lives in structural_adequacy_fit.py.
#
# Writes the fitted values back into params/parameters.json as new per-congener
# fields f_xy_oryza / L_Ph_oryza / kappa_d_oryza (the legacy *_W2fit are PRESERVED
# for provenance / reproduce_demo) and a params/refit_oryza.csv artifact.
#
#   python validation/refit_oryza.py
# =============================================================================
import json, csv, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)            # reuse the fast ORYZA-driven machinery
sys.path.insert(0, os.path.join(ROOT, "src"))
from structural_adequacy_fit import (  # noqa: E402
    _predict, _fit_scalar, obs, CONG, FXY, LPH, KD, _le)

PARJSON = os.path.join(ROOT, "params", "parameters.json")
TISS = ("root", "straw", "grain")
N_PASS = 3                          # coupling passes (root/straw/grain are near-separable)


def refit_one(c):
    """Saturated per-congener fit on ORYZA biomass: kappa_d<-root, f_xy<-straw,
    L_Ph<-grain, iterated for the weak cross-coupling."""
    f_xy = float(c["f_xy_recommended"]); L_Ph = 0.003; kappa_d = 2.0
    o = obs[c["name"]]
    for _ in range(N_PASS):
        if "root" in o:
            kappa_d = _fit_scalar(c, "root", "kappa_d", *KD, f_xy, L_Ph, kappa_d)
        if "straw" in o:
            f_xy = _fit_scalar(c, "straw", "f_xy", *FXY, f_xy, L_Ph, kappa_d)
        if "grain" in o:
            L_Ph = _fit_scalar(c, "grain", "L_Ph", *LPH, f_xy, L_Ph, kappa_d)
    return f_xy, L_Ph, kappa_d


def main():
    print("RE-FIT per-congener (f_xy, L_Ph, kappa_d) on ORYZA2000 biomass vs Yamazaki\n")
    print(f"{'PFAS':8}{'nC':>3} | {'f_xy':>9}{'L_Ph':>10}{'kappa_d':>9} | "
          f"{'root p/o':>15}{'straw p/o':>15}{'grain p/o':>15}")
    fitted = {}
    errs = []
    for c in CONG:
        f_xy, L_Ph, kappa_d = refit_one(c)
        fitted[c["name"]] = dict(f_xy=f_xy, L_Ph=L_Ph, kappa_d=kappa_d)
        pr = _predict(c, f_xy, L_Ph, kappa_d); o = obs[c["name"]]
        errs += [_le(pr[k], o[k]) for k in TISS if k in o]
        print(f"{c['name']:8}{c['n_C']:>3} | {f_xy:>9.4f}{L_Ph:>10.5f}{kappa_d:>9.3f} | "
              f"{pr['root']:>6.2f}/{o.get('root', float('nan')):<7.2f}"
              f"{pr['straw']:>6.2f}/{o.get('straw', float('nan')):<7.2f}"
              f"{pr['grain']:>6.2f}/{o.get('grain', float('nan')):<7.2f}")
    rmse = float(np.sqrt(np.mean(errs)))
    print(f"\nreproduction log10 RMSE (ORYZA2000 biomass, per-congener fit) = {rmse:.3f}"
          f"   (saturated per congener -> reproduces; cf. placeholder W2 0.029)")

    # ---- write artifacts -------------------------------------------------
    with open(os.path.join(ROOT, "params", "refit_oryza.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["compound", "n_C", "f_xy_oryza", "L_Ph_oryza", "kappa_d_oryza"])
        for c in CONG:
            v = fitted[c["name"]]; w.writerow([c["name"], c["n_C"],
                                               f"{v['f_xy']:.6g}", f"{v['L_Ph']:.6g}", f"{v['kappa_d']:.6g}"])

    par = json.load(open(PARJSON))
    for c in par["congeners"]:
        if c["name"] in fitted:
            v = fitted[c["name"]]
            c["f_xy_oryza"] = round(v["f_xy"], 6)
            c["L_Ph_oryza"] = round(v["L_Ph"], 6)
            c["kappa_d_oryza"] = round(v["kappa_d"], 6)
    par.setdefault("_meta", {})["refit_oryza"] = (
        "Per-congener (f_xy, L_Ph, kappa_d) re-fit to Yamazaki on the MECHANISTIC "
        "ORYZA2000 biomass (oryza_growth) + measured Q_TP, the new default driver "
        f"(validation/refit_oryza.py; reproduction log10 RMSE {rmse:.3f}). Saturated "
        "per congener (3 params / 3 obs) -> reproduces by construction; the *_W2fit "
        "values (placeholder/growth_rice driver) are preserved for reproduce_demo. "
        "Use model_api f_xy_source='oryza' to apply these with biomass='oryza'.")
    json.dump(par, open(PARJSON, "w"), indent=2, ensure_ascii=False)
    print(f"\nwrote f_xy_oryza/L_Ph_oryza/kappa_d_oryza for {len(fitted)} congeners -> "
          f"params/parameters.json (+ params/refit_oryza.csv)")
    return fitted, rmse


if __name__ == "__main__":
    main()
