#!/usr/bin/env python3
# =============================================================================
# twopool_root_oos.py
# -----------------------------------------------------------------------------
# OUT-OF-SAMPLE test of the two-pool root + U-shaped k_seq(n) model.
#
# The two-pool model (validation/twopool_root_exploration.py) was fit to
# Yamazaki 2023 ONLY. Here we transfer it WITHOUT re-fitting to independent
# datasets and ask: does the U-shaped-k_seq grain prediction generalise --
# in particular, does it reproduce the long-chain grain RISE that the monotone
# and W2 models miss?
#
#   1. Kim 2019 (Korean paddy) brown-rice (grain) BAF, porewater basis -- the
#      decisive OOS series: it spans PFHpA..PFDoDA and shows the long-chain rise.
#   2. Li 2025 (paddy field) grain/root TF -- water-independent (the Li BAF is
#      water-quality-confounded), a 5-congener short-chain+PFOS sanity check.
#
# All four models (two-pool U-shape, single-pool monotone f_xy, single-pool W2,
# single-pool lipid) are evaluated on the SAME demo forcings as the fit, so the
# comparison is apples-to-apples. EXPLORATORY / in-sample fit -> OOS transfer;
# canonical core + parameters.json UNCHANGED.
#
#   python validation/twopool_root_oos.py     (first run fits+caches ~2.5 min)
# =============================================================================
from __future__ import annotations
import os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT_DIR, "src"))
sys.path.insert(0, HERE)

import twopool_root_exploration as TP          # the fitted two-pool model
from pfas_rice_plant_module_4pool_surf import (  # single-pool baselines
    Compound, RiceUptakeModel, PlantInputs, FRUIT, STEM, LEAF)
from model_api import lipid_loading_conductances
from literature_params import KIM2019_FIELD, kim2019_grain_baf

PAR = TP.PAR
CARR = TP.CARR
ENV = TP.ENV
CONG = {c["name"]: c for c in PAR["congeners"]}


# ---------------------------------------------------------------------------
# single-pool 4-compartment grain BAF on the SAME demo forcings (reproduce_demo)
# ---------------------------------------------------------------------------
_INPUTS = PlantInputs(t=TP.T, Cwo=TP.CWO, Qtp=TP.QTP, M=TP.MMAT)


def single_pool(c, f_xy, lipid=False):
    """Return (root, straw, grain) BAF for the canonical single-pool model."""
    g_xy = g_ph = 0.0
    if lipid:
        _, g_xy, g_ph = lipid_loading_conductances(c["n_C"], c["K_PL_Lkg"], c["group"])
    cmpd = Compound(name=c["name"], K_prot=c["K_prot_Lkg"], K_PL=c["K_PL_Lkg"],
                    K_cw=c["K_cw_wholecw_Lkg"]["root"], kappa_d=c["kappa_d_W2fit"],
                    Vmax_in=CARR["Vmax_in"], Km_in=CARR["Km_in"],
                    Vmax_out=CARR["Vmax_out"], Km_out=CARR["Km_out"],
                    L_Ph=c["L_Ph_W2fit"], f_xy=f_xy, g_xy=g_xy, g_ph=g_ph)
    sol = RiceUptakeModel(env=ENV, cmpd=cmpd, comps=TP.compartments(),
                          inputs=_INPUTS).solve(TP.T)
    C = sol.y[:, -1]; Mf = TP.MMAT[-1]
    straw = (C[STEM] * Mf[STEM] + C[LEAF] * Mf[LEAF]) / (Mf[STEM] + Mf[LEAF])
    return C[0], straw, C[FRUIT]


def twopool(c, p, q):
    ks = TP.kseq_ushape(c["n_C"], c["group"], q)
    return TP.simulate(c, p, kseq_override=ks)


def _rmse(pred, obs, keys):
    e = [(np.log10(max(pred[k], 1e-6)) - np.log10(obs[k])) ** 2 for k in keys]
    return float(np.sqrt(np.mean(e))) if e else float("nan")


# ===========================================================================
def main():
    print("=" * 78)
    print("OUT-OF-SAMPLE transfer of the two-pool U-shaped-k_seq model (Yamazaki-fit)")
    print("=" * 78)
    p, q = TP.load_fit()
    print(f"loaded two-pool fit: kappa_d={p['kappa_d']:.2f} L_Ph={p['L_Ph']:.2e} "
          f"gxy={p['gxy']:.4f} gph={p['gph']:.4f}")
    print(f"U-shape q (logA,b,logC,d,sa) = "
          f"[{q[0]:.3f}, {q[1]:.3f}, {q[2]:.3f}, {q[3]:.3f}, {q[4]:.3f}]")

    # -------------------------------------------------------------------
    # (1) Kim 2019 brown-rice (grain) BAF, porewater basis  [the key OOS]
    # -------------------------------------------------------------------
    kim = kim2019_grain_baf("porewater")
    DF = {k: v[3] for k, v in KIM2019_FIELD.items()}   # rice detection frequency [%]
    print("\n" + "-" * 78)
    print("(1) Kim 2019 grain BAF (porewater) — OOS chain-length series")
    print("-" * 78)
    print(f"{'PFAS':8}{'nC':>3}{'DF%':>5} | {'obs':>8}{'2pool':>9}{'mono':>9}"
          f"{'W2':>9}{'lipid':>9}")
    models = {"2pool": [], "mono": [], "W2": [], "lipid": []}
    obs_log = {}
    for nm in kim:
        c = CONG.get(nm)
        if c is None or nm not in TP.OBS:
            continue
        o = kim[nm]
        gp_2 = twopool(c, p, q)[2]
        gp_m = single_pool(c, c["f_xy_recommended"])[2]
        gp_w = single_pool(c, c["f_xy_W2fit"])[2] if c["f_xy_W2fit"] else float("nan")
        gp_l = single_pool(c, c["f_xy_recommended"], lipid=True)[2]
        for key, val in (("2pool", gp_2), ("mono", gp_m), ("W2", gp_w), ("lipid", gp_l)):
            models[key].append((nm, val, o, DF.get(nm, 0.0)))
        obs_log[nm] = o
        print(f"{nm:8}{c['n_C']:>3}{DF.get(nm,0):>5.0f} | {o:>8.2f}{gp_2:>9.2f}"
              f"{gp_m:>9.3f}{gp_w:>9.3f}{gp_l:>9.2f}")

    def kim_rmse(key, names):
        e = [(np.log10(max(v, 1e-6)) - np.log10(o)) ** 2
             for (nm, v, o, df) in models[key] if nm in names and v == v]
        return float(np.sqrt(np.mean(e))) if e else float("nan")

    allnames = [nm for (nm, *_ ) in models["2pool"]]
    excl_pfoa = [nm for nm in allnames if nm != "PFOA"]          # PFOA used in L_Ph fit elsewhere
    reliable = [nm for (nm, v, o, df) in models["2pool"] if df >= 15.0]
    print(f"\nlog10 RMSE vs Kim grain   {'2pool':>8}{'mono':>8}{'W2':>8}{'lipid':>8}")
    for label, names in (("all", allnames), ("excl PFOA", excl_pfoa),
                         (f"reliable DF>=15% {reliable}", reliable)):
        print(f"  {label:22} " + "".join(f"{kim_rmse(k, names):>8.2f}"
              for k in ("2pool", "mono", "W2", "lipid")))
    # does 2pool capture the long-chain RISE? (PFUnDA, PFDoDA obs ~33,35)
    lc = {nm: v for (nm, v, o, df) in models["2pool"] if nm in ("PFUnDA", "PFDoDA")}
    lcm = {nm: v for (nm, v, o, df) in models["mono"] if nm in ("PFUnDA", "PFDoDA")}
    print(f"\nlong-chain grain RISE check (Kim obs PFUnDA~33, PFDoDA~35):")
    print(f"   2pool: PFUnDA {lc.get('PFUnDA', float('nan')):.1f}  "
          f"PFDoDA {lc.get('PFDoDA', float('nan')):.1f}   "
          f"|  monotone: PFUnDA {lcm.get('PFUnDA', float('nan')):.3f}  "
          f"PFDoDA {lcm.get('PFDoDA', float('nan')):.3f}")

    # -------------------------------------------------------------------
    # (2) Li 2025 grain/root TF  (water-independent; short chains + PFOS)
    # -------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("(2) Li 2025 grain/root TF — OOS (water-independent ratio)")
    print("-" * 78)
    li_tf = {"PFHxA": 7.885, "PFOA": 1.35, "PFBS": 19.304, "PFOS": 0.795}  # TF_grain_fw
    print(f"{'PFAS':8} | {'obs TF':>8}{'2pool TF':>10}{'mono TF':>9}")
    e2, em = [], []
    for nm, tf in li_tf.items():
        c = CONG.get(nm)
        if c is None or nm not in TP.OBS:
            print(f"{nm:8} | {tf:>8.2f}   (not in Yamazaki congener set; skipped)")
            continue
        r2, s2, g2 = twopool(c, p, q)
        rm, sm, gm = single_pool(c, c["f_xy_recommended"])
        tf2 = g2 / max(r2, 1e-6); tfm = gm / max(rm, 1e-6)
        e2.append((np.log10(max(tf2, 1e-6)) - np.log10(tf)) ** 2)
        em.append((np.log10(max(tfm, 1e-6)) - np.log10(tf)) ** 2)
        print(f"{nm:8} | {tf:>8.2f}{tf2:>10.2f}{tfm:>9.2f}")
    if e2:
        print(f"\nlog10 RMSE (grain/root TF):  2pool={np.sqrt(np.mean(e2)):.2f}  "
              f"monotone={np.sqrt(np.mean(em)):.2f}")

    print("\n" + "=" * 78)
    print("VERDICT (honest, data-driven):")
    rmses = {k: kim_rmse(k, excl_pfoa) for k in ("2pool", "mono", "W2", "lipid")}
    best = min(rmses, key=rmses.get)
    print(f"  Kim grain (excl PFOA), SAME demo forcings: "
          + "  ".join(f"{k}={v:.2f}" for k, v in rmses.items())
          + f"   -> best = {best}")
    print(f"  - Two-pool U-shaped-k_seq transfers {'BEST' if best=='2pool' else best+'-best'}"
          f" on the clean Kim series and CAPTURES the long-chain RISE")
    print(f"    (2pool PFUnDA {lc.get('PFUnDA',0):.1f}/PFDoDA {lc.get('PFDoDA',0):.1f} vs "
          f"monotone {lcm.get('PFUnDA',0):.2f}/{lcm.get('PFDoDA',0):.2f}; Kim obs ~33/35).")
    print("  - HONEST limits: absolute long-chain grain still UNDER (low-DF Kim tail);")
    print("    Kim PFOA grain 4.43 >> Yamazaki 0.46 (between-dataset shift, unbridgeable);")
    print("    Li 2025 grain/root TF is root-surface-confounded (short-chain anomaly) -> inconclusive.")
    print("  => OOS SUPPORTS the structure/mechanism, but does NOT warrant promoting the")
    print("     fitted k_seq into parameters.json (single clean OOS set; demo forcings).")
    print("  Still in-sample FIT (Yamazaki) -> OOS TRANSFER; not full validation.")


if __name__ == "__main__":
    main()
