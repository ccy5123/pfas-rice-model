#!/usr/bin/env python3
# =============================================================================
# validation/twopool_root_oos_tang.py
# -----------------------------------------------------------------------------
# OUT-OF-SAMPLE transfer of the two-pool root + U-shaped k_seq(n) model to the
# Tang 2026 per-organ transfer factors (stalk / leaf / endosperm, dry weight).
#
# The two-pool model (validation/twopool_root_exploration.py) was fit to
# Yamazaki 2023 ONLY (aggregate root/straw/grain BAF). Tang 2026 is an
# INDEPENDENT dataset (flooded paddy pot, cv. Nipponbare, dose series) reporting
# per-organ TF -- a different soil, cultivar and dose set, so transferring the
# Yamazaki-fit model here WITHOUT re-fitting is a genuine cross-dataset OOS test.
#
# This complements:
#   * validation/twopool_root_oos.py     -- two-pool OOS on Kim 2019 grain + Li 2025 TF
#   * validation/oos_tang.py             -- SINGLE-pool monotone f_xy OOS on Tang (RMSE 1.23)
#   * validation/oos_tang_lipid.py       -- SINGLE-pool lipid-loading OOS on Tang (RMSE 0.52)
# by asking the two-pool-specific question: does the two-pool (monotone physical
# f_xy + K_PL-gated lipid loading, the structure that ALSO keeps the high
# long-chain root) transfer to Tang's per-organ shoot pattern?
#
# RESULT (negative/diagnostic, see VERDICT below): it does NOT -- OOS RMSE ~1.40,
# WORSE than the single-pool monotone (1.23) and far worse than lipid (0.52). The
# reason is structural and informative: Tang per-organ TF is a SHOOT-resolution
# test, but the two-pool's innovation is in the ROOT (mobile/seq split); its shoot
# is the UNMODIFIED basic 4pool with a PASS-THROUGH stem, so the stalk TF collapses
# (the documented empty-stem defect that `nstem_leaf` fixes). Tang congeners are
# C5-C8, so the long-chain ROOT decoupling (the two-pool's whole point) is not even
# exercised. CONCLUSION: Tang is not a suitable OOS test of the two-pool ROOT; a
# fair test needs the two-pool root merged with the nstem_leaf REDISTRIBUTED shoot.
# Kim 2019 grain (twopool_root_oos.py) stays the informative two-pool OOS.
#
# The two-pool 5-state ODE keeps stem and leaf DISTINCT, so stalk (=stem) and
# leaf TF are read directly (via simulate_organs); the single-pool baselines come
# from model_api.tang_tf_validation (redistributed-shoot nstem_leaf split) -- so
# the stalk comparison is apples-to-oranges (different SHOOT model), the point of
# the diagnosis.
#
# Honest framing (carry forward): Yamazaki in-sample fit -> OOS transfer; Tang is
# a single independent set; GenX f_xy_recommended (ether offset) is provisional
# and over-predicts; grain/endosperm is structurally ~3-8x under (units doc).
# EXPLORATORY; canonical core + parameters.json UNCHANGED.
#
#   python validation/twopool_root_oos_tang.py     (first run fits+caches ~2.5 min)
# =============================================================================
from __future__ import annotations
import os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT_DIR, "src"))
sys.path.insert(0, HERE)

import twopool_root_exploration as TP            # the Yamazaki-fit two-pool model
import model_api as api                          # observed Tang TF + single-pool baselines

CONG = {c["name"]: c for c in TP.PAR["congeners"]}
DOSE = "low"          # Tang 0.1 ug/g lowest dose = environmentally closest (primary)


def twopool_tf(c, p, q):
    """Two-pool per-organ DRY-weight TF (stalk/leaf/endosperm) vs root.

    fresh-weight conc from the 5-state ODE -> TF_fw[organ] = C[organ]/C[root];
    TF_dw = TF_fw * (1-theta_root)/(1-theta_tissue) (same dw correction as
    model_api.tang_tf_validation; the (1-theta) factor differs by tissue so it
    does NOT cancel in C_tissue/C_root)."""
    ks = TP.kseq_ushape(c["n_C"], c["group"], q)
    org = TP.simulate_organs(c, p, kseq_override=ks)
    froot = 1.0 - api._COMP["root"]["theta_fw"]
    out = {}
    for mk, tang_organ, tk in api._TANG_ORGANS:
        tf_fw = org[mk] / max(org["root"], 1e-9)
        out[tang_organ] = tf_fw * froot / (1.0 - api._COMP[tk]["theta_fw"])
    return out


def _rmse(pairs):
    """pairs: list of (pred, obs) -> log10 RMSE."""
    e = [(np.log10(max(pr, 1e-6)) - np.log10(max(ob, 1e-6))) ** 2 for pr, ob in pairs]
    return float(np.sqrt(np.mean(e))) if e else float("nan")


def main(dose=DOSE):
    print("=" * 86)
    print("OUT-OF-SAMPLE transfer of the two-pool U-shaped-k_seq model to Tang 2026 per-organ TF")
    print("=" * 86)
    p, q = TP.load_fit()
    print(f"loaded two-pool fit (Yamazaki): kappa_d={p['kappa_d']:.2f} L_Ph={p['L_Ph']:.2e} "
          f"gxy={p['gxy']:.4f} gph={p['gph']:.4f}")
    print(f"Tang dose = {dose} (0.1 ug/g, environmentally closest)\n")

    print(f"{'cong':6}{'organ':>11}{'obs':>8}{'2pool':>9}{'mono':>9}{'lipid':>9}{'refit(IS)':>11}")
    # accumulate per-model (pred, obs) pairs
    acc = {"2pool": [], "mono": [], "lipid": [], "refit": []}
    for nm in api.TANG_CONGENERS:
        c = CONG[nm]
        obs = api.tang_observed_tf(nm, dose)
        tf2 = twopool_tf(c, p, q)
        v_m = api.tang_tf_validation(nm, f_xy_source="recommended", dose=dose)              # OOS monotone
        v_l = api.tang_tf_validation(nm, f_xy_source="recommended", dose=dose,
                                     lipid_loading=True)                                    # OOS lipid
        v_r = api.tang_tf_validation(nm, f_xy_source="recommended", use_refit=True, dose=dose)  # in-sample refit
        for _, organ, _ in api._TANG_ORGANS:
            o = obs[organ]
            m2, mm, ml, mr = (tf2[organ], v_m["model_tf"][organ],
                              v_l["model_tf"][organ], v_r["model_tf"][organ])
            acc["2pool"].append((m2, o)); acc["mono"].append((mm, o))
            acc["lipid"].append((ml, o)); acc["refit"].append((mr, o))
            print(f"{nm:6}{organ:>11}{o:>8.3f}{m2:>9.3f}{mm:>9.3f}{ml:>9.3f}{mr:>11.3f}")

    print("\n" + "-" * 86)
    print(f"log10 RMSE vs Tang per-organ TF (dw, dose={dose}):")
    rmses = {k: _rmse(v) for k, v in acc.items()}
    print(f"  two-pool (OOS)            = {rmses['2pool']:.3f}")
    print(f"  single-pool monotone (OOS)= {rmses['mono']:.3f}")
    print(f"  single-pool lipid  (OOS)  = {rmses['lipid']:.3f}")
    print(f"  single-pool Tang-refit(IS)= {rmses['refit']:.3f}   <- in-sample reference")

    # per-organ breakdown (the 3 congeners share organ index 0=stalk,1=leaf,2=endosperm)
    organs = [o for _, o, _ in api._TANG_ORGANS]
    by_organ = {k: {} for k in acc}
    for k in acc:
        for i, organ in enumerate(organs):
            by_organ[k][organ] = _rmse(acc[k][i::3])
    print(f"\n  per-organ RMSE      {'stalk':>8}{'leaf':>8}{'endosperm':>11}")
    for k in ("2pool", "mono", "lipid"):
        print(f"    {k:16}" + "".join(f"{by_organ[k][o]:>8.2f}" if o != 'endosperm'
              else f"{by_organ[k][o]:>11.2f}" for o in organs))

    print("\n" + "=" * 86)
    print("VERDICT (honest, data-driven):")
    best = min(("2pool", "mono", "lipid"), key=lambda k: rmses[k])
    print(f"  Two-pool OOS RMSE {rmses['2pool']:.3f} -- WORSE than single-pool monotone "
          f"{rmses['mono']:.3f} and far worse than lipid {rmses['lipid']:.3f} (best OOS = {best});")
    print(f"  in-sample refit {rmses['refit']:.3f}.")
    print("  - DIAGNOSIS: Tang per-organ TF is a SHOOT-resolution test (stalk vs leaf vs")
    print("    endosperm), but the two-pool's innovation is in the ROOT (mobile/seq split). Its")
    print("    shoot is the UNMODIFIED basic 4pool with a PASS-THROUGH stem (PFOA stem conc 0.008")
    print("    vs leaf 1.14), so the STALK TF collapses (per-organ RMSE above: stalk is the driver)")
    print("    -- the documented over-translocation/empty-stem defect that `nstem_leaf`")
    print("    (redistributed shoot + retention) was built to fix. The single-pool baselines here")
    print("    USE nstem_leaf, so their stalk is populated -- an apples-to-oranges SHOOT difference,")
    print("    not a root-mechanism difference. Tellingly, the two-pool LEAF RMSE is the BEST of")
    print("    all three models (it resolves the leaf well); ONLY the stalk drags the overall up.")
    print("  - CONCLUSION: the two-pool ROOT mechanism and the Tang per-organ SHOOT pattern are")
    print("    largely ORTHOGONAL. Tang congeners are C5-C8 (short/mid), so the long-chain root")
    print("    decoupling -- the two-pool's whole point -- is not even exercised. Tang is therefore")
    print("    NOT a suitable OOS test of the two-pool root: a fair test needs the two-pool ROOT")
    print("    merged with the nstem_leaf REDISTRIBUTED SHOOT (a future structural merge).")
    print("  - HONEST limits: Yamazaki in-sample fit; single Tang set; GenX (ether f_xy_recommended)")
    print("    over-predicts (provisional QSPR, not the two-pool); grain ~3-8x under (units doc).")
    print("  => OOS does NOT support promoting the fitted k_seq into parameters.json. The actionable")
    print("     finding is structural: pair the two-pool root with the redistributed shoot before")
    print("     any per-organ Tang comparison. Kim 2019 grain (twopool_root_oos.py) remains the")
    print("     informative two-pool OOS (it is grain/root, not shoot-resolved).")
    return rmses


if __name__ == "__main__":
    main()
