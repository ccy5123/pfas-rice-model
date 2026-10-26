#!/usr/bin/env python3
# =============================================================================
# reproduce_demo.py -- self-contained reproduction entry point
# -----------------------------------------------------------------------------
# Runs the model on the Yamazaki 2020 site-year (JMA Tsukuba weather, the ORYZA2000 port on
# the paper's transplanting / heading / harvest dates, ORYZA2000 transpiration;
# src/site_scenario.py) and compares the harvest values (dry basis; straw with attached dead
# blades; src/harvest.py) with Yamazaki 2023's root / straw / grain
# (data_obs/obs_baf_Yamazaki.csv, built from the SI).
#
#   python reproduce_demo.py           # the DEFAULT configuration (redistributed_fit; amendment
#                                      #   2a): r_eq, f_xy and L_Ph fitted per congener to root,
#                                      #   straw and grain -> a REPRODUCTION (3 params / 3 cells)
#   python reproduce_demo.py --lipid   # redistributed_lipid (amendment 1's default): the straw
#                                      #   predicted by two lipid constants shared by all chains
#   python reproduce_demo.py --free    # redistributed_free: the monotone a-priori f_xy
#
# The fits are made by validation/refit_sourced.py (pre-registered,
# docs/PREREG_sourced_values.md). PFBS grain is wholly < MQL and is not scored.
# =============================================================================
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "src"))
import model_api as api   # noqa: E402

VARIANTS = {"default": dict(), "lipid": dict(lipid_loading=True),
            "free": dict(lipid_loading=False)}
CHOICE = next((v for v in ("lipid", "free") if f"--{v}" in sys.argv), "default")
TISS = ("root", "straw", "grain")


def main(choice=CHOICE, verbose=True):
    say = print if verbose else (lambda *a, **k: None)
    kw = VARIANTS[choice]
    v = api.variant_key(kw.get("shoot"), kw.get("f_xy_source"), kw.get("lipid_loading"))
    say(f"variant: {v}  (Yamazaki 2020 site-year, dry basis, two harvest runs)")
    say(f"{'PFAS':8}{'nC':>3}{'f_xy':>8} | " + "".join(f"{k + ' m/o':>18}" for k in TISS))
    err = {k: [] for k in TISS}
    out = {}
    for c in api.CONGENERS:
        o = api.observed_baf(c)
        if not o:
            continue
        pr = api.yamazaki_harvest_baf(c, **kw)
        out[c] = pr
        for k in o:
            err[k].append(np.log10(pr[k] / o[k]))
        fxy = api._resolve_transport(api._CONG[c], kw.get("f_xy_source"),
                                     kw.get("lipid_loading"), kw.get("shoot"))[0]
        say(f"{c:8}{api._CONG[c]['n_C']:>3}{fxy:>8.4f} | "
            + "".join(f"{pr[k]:>8.3g}/{o.get(k, float('nan')):<9.3g}" for k in TISS))
    allerr = [e for vals in err.values() for e in vals]
    rmse = float(np.sqrt(np.mean(np.square(allerr))))
    per = {k: float(np.sqrt(np.mean(np.square(vals)))) for k, vals in err.items()}
    say(f"\nlog10 RMSE = {rmse:.3f}  (" + ", ".join(f"{k} {x:.3f}" for k, x in per.items()) + ")")
    if v.endswith("_fit"):
        say("NOTE: a saturated fit (3 parameters per congener vs 3 cells) -- a REPRODUCTION, not a "
            "prediction.")
    else:
        say("NOTE: root and grain are fitted per congener (r_eq, L_Ph; grain where L_Ph reaches it); "
            "the straw is not fitted per congener.")
    return dict(rmse=rmse, per_tissue=per, pred=out, variant=v)


if __name__ == "__main__":
    main()
