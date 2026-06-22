#!/usr/bin/env python3
# =============================================================================
# validation/longchain_mechanism.py
# -----------------------------------------------------------------------------
# Long-chain (C10-C12) mechanism investigation, on the ORYZA2000 biomass.
#
# The re-fit showed PFDoDA (C12) is unreachable even at all-parameter ceilings:
# straw/grain saturate far below the observed values. This probes WHY and whether
# the lipid-facilitated bound-loading term closes it.
#
# DIAGNOSIS (free-anion loading throttle): the xylem/phloem loading uses the FREE
# aqueous conc Cw = C/B. B grows ~10^(0.5-0.6 per CF2) with chain length, so the
# free-loading flux f_xy*Cw collapses ~1/B for long chains -- the shoot starves no
# matter how large f_xy is (PFDoDA hits f_xy=1 and still under-predicts).
#
# FIX (lipid-facilitated bound loading): the opt-in g_xy*C_root / g_ph*C_leaf terms
# are B-INDEPENDENT (use total C, not free Cw), so the membrane/lipid-associated
# pool rides into the xylem/phloem and does not collapse for long chains.
#
# Compares, per congener on ORYZA2000 biomass: free-only (recommended monotone f_xy)
# vs lipid-loading, root/straw/grain BAF vs Yamazaki; reports the long-chain
# (nC>=10) straw+grain RMSE for each and the Cw=C/B collapse.
#
#   python validation/longchain_mechanism.py
# =============================================================================
import csv, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
import model_api as api   # noqa: E402

obs = {}
with open(os.path.join(ROOT, "data_obs", "obs_baf_Yamazaki.csv"), newline="") as f:
    for r in csv.DictReader(f):
        obs.setdefault(r["compound"], {})[r["tissue"]] = float(r["baf"])

TISS = ("root", "straw", "grain")


def _baf(r):
    return {"root": r["baf_final"]["root"], "straw": r["straw_baf"], "grain": r["baf_final"]["grain"]}


def _le(a, b):
    return (np.log10(max(a, 1e-6)) - np.log10(max(b, 1e-6))) ** 2


def _rmse(pairs):
    return float(np.sqrt(np.mean(pairs))) if pairs else float("nan")


def main():
    print("Long-chain mechanism on ORYZA2000 biomass: FREE-only (monotone f_xy) vs LIPID-loading\n")
    print(f"{'PFAS':7}{'nC':>3}{'B_root':>9}{'Cw_root':>9} | "
          f"{'straw free/lip/obs':>22}{'grain free/lip/obs':>22}{'root free/lip/obs':>22}")
    e_free = {"all": [], "long": [], "lc_sg": []}
    e_lip = {"all": [], "long": [], "lc_sg": []}
    for nm in api.CONGENERS:
        if nm not in obs:
            continue
        free = api.simulate(nm, f_xy_source="recommended", lipid_loading=False, biomass="oryza")
        lip = api.simulate(nm, lipid_loading=True, biomass="oryza")
        nC = api._CONG[nm]["n_C"]
        bf, bl, o = _baf(free), _baf(lip), obs[nm]
        B_root = free["B_k"]["root"]
        Cw_root = free["conc"]["root"][-1] / B_root          # free aqueous conc available to load
        for k in TISS:
            if k in o:
                e_free["all"].append(_le(bf[k], o[k])); e_lip["all"].append(_le(bl[k], o[k]))
                if nC >= 10:
                    e_free["long"].append(_le(bf[k], o[k])); e_lip["long"].append(_le(bl[k], o[k]))
                    if k in ("straw", "grain"):
                        e_free["lc_sg"].append(_le(bf[k], o[k])); e_lip["lc_sg"].append(_le(bl[k], o[k]))
        def c3(k):
            return f"{bf[k]:.2f}/{bl[k]:.2f}/{o.get(k, float('nan')):.2f}"
        print(f"{nm:7}{nC:>3}{B_root:>9.1f}{Cw_root:>9.4f} | "
              f"{c3('straw'):>22}{c3('grain'):>22}{c3('root'):>22}")

    print(f"\nlog10 RMSE                         free-only   lipid")
    print(f"  long-chain (nC>=10) straw+grain   {_rmse(e_free['lc_sg']):.3f}      {_rmse(e_lip['lc_sg']):.3f}")
    print(f"  long-chain (nC>=10) all tissues   {_rmse(e_free['long']):.3f}      {_rmse(e_lip['long']):.3f}")
    print(f"  whole series       all tissues    {_rmse(e_free['all']):.3f}      {_rmse(e_lip['all']):.3f}")
    print("\nCw_root = C_root/B_root is the FREE conc the xylem loads (f_xy*Cw). It collapses")
    print("as B_root grows with chain length -> the free-only shoot starves the long chains;")
    print("the B-independent lipid term (g_xy*C, g_ph*C) is what does NOT collapse.")
    return e_free, e_lip


if __name__ == "__main__":
    main()
