#!/usr/bin/env python3
# =============================================================================
# validation/li2019_rcf_apriori.py
# -----------------------------------------------------------------------------
# The second a-priori test of the neutral path's PARTITION term, and the
# diagnosis it forced.
#
# WHAT ARRIVED. docs/literature_db/Acquisition_Queue.csv row A1 -- Li, Chiou, Li
# & Schnoor 2019, Environ. Int. 126:46-53 -- was requested for a rice ROOT LIPID
# value. Its SI does NOT contain one (no rice row anywhere in Tables S2/S4; the
# article's own crop list is wheat, barley, carrot, radish, celery, maize,
# pumpkin, turnip, onion, spinach, Chinese cabbage, ryegrass, amaranth). So the
# top-ranked gap in the queue is NOT closed by the paper that was supposed to
# close it. What it delivered instead is worth more, and this script is it:
#
#   1. Table S1 = 48 hydroponic root concentration factors over 11 species,
#      log Kow -0.57..5.41. Shipped as data_obs/neutral_obs_li2019_rcf.csv, it
#      is a SECOND independent a-priori test of K_PW -- twice the size of the
#      Liu 2023 table, over a wider range, and including four rice rows.
#   2. The operational DEFINITION of the lipid fraction, verified at source,
#      which is what the queue's own DEFINITION NOTE said mattered more than any
#      number: Li et al.'s f_lip is FRESH weight (they convert dry-basis reports
#      at 90% root water), and it enters an RCF expression of exactly this
#      model's form. So their per-crop values ARE commensurable with this
#      model's L, and their cereals -- barley 1.00%, wheat 1.10-1.14%, maize
#      0.53% -- bracket the 1% this repo already runs.
#
# THE DIAGNOSIS. Running the test exposed something the existing tables could
# not: the model's rice root partition is systematically LOW, by a near-constant
# offset in every one of the 11 species (mean log10 bias -0.30 to -0.95). The
# cause is internal, not in the data. `neutral_dpu` anchors on Briggs' RCF, whose
# lipid term is L*a = 10^-1.52 = 0.0302; but `rice_compartments` substitutes a
# MEASURED lipid L = 0.01 while keeping the conventional a = 1.22, which silently
# changes the anchored product to 0.0122 -- 2.5x below Briggs. The module knows L
# and a are identifiable only as their product; the rice compartment nonetheless
# replaces one factor and not the other.
#
# Section 3 shows what that costs on Briggs' OWN data, which is the sharpest form
# of the point: the shipped composition fits the 18 barley rows the RCF QSPR was
# fitted to at log10 RMSE 0.266, where the anchor itself reaches 0.111.
#
# WHY THE DEFAULT IS NOT CHANGED HERE. Restoring the anchor improves two of the
# three tables and DEGRADES the third -- and the third, Liu 2023, is the only one
# measured on rice. An intermediate L would fit all three, but fitting L is
# exactly what this path must not do: the neutral path's whole claim is that
# K_PW and TSCF come from log Kow with nothing tuned. So the anchor is offered as
# an opt-in (`neutral_dpu.BRIGGS_ANCHORED_LIPID_FW`), the discrepancy is pinned
# by a test, and the promotion decision is left to the user.
#
#   python validation/li2019_rcf_apriori.py            # full run (~3 min)
#   python validation/li2019_rcf_apriori.py --fast     # skip the ODE scan
# =============================================================================
from __future__ import annotations
import csv, os, sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT_DIR, "src"))
sys.path.insert(0, HERE)

import neutral_dpu as ND                                       # noqa: E402
import neutral_dpu_validation as V                             # noqa: E402

LI2019 = os.path.join(ROOT_DIR, "data_obs", "neutral_obs_li2019_rcf.csv")
LIU2023 = os.path.join(ROOT_DIR, "data_obs", "neutral_obs_liu2023.csv")
GE2017 = os.path.join(ROOT_DIR, "data_obs", "neutral_obs_ge2017.csv")

# Root lipid contents Li et al. 2019 report for crops, FRESH weight (article
# p.49 and SI Table S2). Rice is absent -- that absence is the finding.
LI2019_ROOT_LIPID_FW = {
    "wheat": 0.0114, "barley": 0.0100, "maize": 0.0053, "pumpkin": 0.0070,
    "chinese cabbage": 0.0068, "spinach": 0.0034, "ryegrass": 0.0032,
    "amaranth": 0.0032, "carrot": 0.0024, "celery": 0.0017, "radish": 0.0010,
    "turnip": 0.0010, "onion": 0.0010,
}


def _root_rows(path, subset=None):
    """(log Kow, measured root BAF) for the root rows of an obs table."""
    lk, v = [], []
    with open(path, newline="") as f:
        for r in csv.DictReader(x for x in f if not x.lstrip().startswith("#")):
            if not r.get("compound") or r.get("tissue") != "root":
                continue
            if subset is not None and r.get("subset") != subset:
                continue
            lk.append(float(r["log_kow"]))
            v.append(float(r["value"]))
    return np.array(lk), np.array(v)


def _kpw_rmse(lk, obs, La, b, W=None):
    """log10 RMSE of the EQUILIBRIUM partition K_PW = W + L*a*Kow^b.

    Used only for the anchor diagnosis. It is not the plant model: the ODE's root
    BAF sits a few percent under K_PW (growth dilution and xylem export), so
    these numbers are close to, but not identical with, section 1's.
    """
    W = ND.RICE_WATER["root"] if W is None else W
    pred = W + La * 10.0 ** (b * lk)
    return float(np.sqrt(np.mean((np.log10(pred) - np.log10(obs)) ** 2)))


def _fit_anchor(lk, obs, W=None):
    """FIT (L*a, b) to a root table. Explicitly not a-priori -- a diagnostic of
    what the data want, to be compared against the two published anchors."""
    from scipy.optimize import least_squares
    W = ND.RICE_WATER["root"] if W is None else W

    def res(p):
        return np.log10(W + 10.0 ** p[0] * 10.0 ** (p[1] * lk)) - np.log10(obs)

    r = least_squares(res, [np.log10(0.0302), ND.RCF_SLOPE])
    return 10.0 ** r.x[0], r.x[1], float(np.sqrt(np.mean(r.fun ** 2)))


def section1_apriori(drv):
    print("=" * 84)
    print("1. A-PRIORI PREDICTION — Li 2019 Table S1 hydroponic RCF, nothing fitted")
    print("=" * 84)
    rmse = V.compare_to_obs(LI2019, drv, group_by="species")
    return rmse


def section2_lipid_table():
    print("\n" + "=" * 84)
    print("2. WHAT A1 ACTUALLY DELIVERED ON THE LIPID QUESTION")
    print("=" * 84)
    print("   Li et al. 2019 root lipid contents, FRESH weight (their conversion:")
    print("   dry-basis reports are rescaled at 90% root water). Rice is ABSENT.")
    for k, v in sorted(LI2019_ROOT_LIPID_FW.items(), key=lambda kv: -kv[1]):
        mark = "  <- cereal" if k in ("wheat", "barley", "maize") else ""
        print(f"      {k:18s} {v * 100:5.2f} %{mark}")
    cereals = [LI2019_ROOT_LIPID_FW[k] for k in ("wheat", "barley", "maize")]
    print(f"\n   cereal mean {np.mean(cereals) * 100:.2f} % fw; this repo runs "
          f"{ND.TRAPP1994_LIPID_FW['root'] * 100:.2f} % (Trapp 1994 soybean).")
    print("   => The VALUE the repo uses is corroborated. Its PROVENANCE improves")
    print("      from a soybean model run to measured cereal roots on a stated,")
    print("      RCF-operational basis. Rice itself remains unmeasured (gap C3).")
    print("\n   Note what this does NOT do: Briggs' fitted lipid term implies")
    print(f"   L = {10 ** ND.RCF_INTERCEPT / ND.LIPID_OCTANOL_A * 100:.2f} % fw for his barley, whereas Li et al. assign")
    print("   that same barley 1.00 %. Briggs 1982 itself reports NO lipid content")
    print("   (verified at source: he attributes the 0.82 floor to root WATER and")
    print("   never measures the lipid), so both figures are inferences, and they")
    print("   disagree 2.5x. Section 3 is what that disagreement costs.")


def section3_anchor():
    print("\n" + "=" * 84)
    print("3. THE ANCHOR DIAGNOSIS — the shipped root sits 2.5x below Briggs")
    print("=" * 84)
    shipped = ND.TRAPP1994_LIPID_FW["root"] * ND.LIPID_OCTANOL_A
    anchor = 10.0 ** ND.RCF_INTERCEPT
    print(f"   Briggs 1982 anchored lipid term  L*a = {anchor:.4f}   (10^{ND.RCF_INTERCEPT})")
    print(f"   what rice_compartments actually runs  = {shipped:.4f}   "
          f"(L={ND.TRAPP1994_LIPID_FW['root']} x a={ND.LIPID_OCTANOL_A})")
    print(f"   ratio {anchor / shipped:.2f}x low\n")
    tables = [
        ("Li2019 hydroponic RCF, a-priori rows", *_root_rows(LI2019, "apriori")),
        ("Li2019 rows Briggs FITTED his QSPR to", *_root_rows(LI2019, "calibration")),
        ("Liu2023 rice root (the only RICE table)", *_root_rows(LIU2023)),
    ]
    cands = [("shipped  L*a=%.4f b=%.2f" % (shipped, ND.RCF_SLOPE), shipped, ND.RCF_SLOPE),
             ("Briggs anchor  L*a=%.4f b=%.2f" % (anchor, ND.RCF_SLOPE), anchor, ND.RCF_SLOPE),
             ("Li/Chiou  f_lip=0.010 K_lip=1.27Kow^1.03", 0.0127, 1.03)]
    print(f"   {'table':42s}" + "".join(f"{n.split()[0]:>12}" for n, _, _ in cands)
          + f"{'fitted':>10}{'L*a':>9}{'b':>7}")
    for name, lk, obs in tables:
        cells = [f"{_kpw_rmse(lk, obs, La, b):>12.3f}" for _, La, b in cands]
        La, b, r = _fit_anchor(lk, obs)
        print(f"   {name:42s}" + "".join(cells) + f"{r:>10.3f}{La:>9.4f}{b:>7.3f}")
    print("\n   equilibrium K_PW log10 RMSE; 'fitted' is 2 free parameters and is")
    print("   NOT a-priori -- it is there to show what each table wants.")
    print("\n   The second row is the sharp one: on the very barley data the RCF")
    print("   QSPR was fitted to, the anchor reaches 0.111 and the shipped")
    print("   composition only 0.266. That is an internal inconsistency, not a")
    print("   disagreement with the world.")
    print("   The third row is why it is not simply fixed: rice prefers the LOW")
    print("   value. The tables disagree, and only one of them is rice.")


def section4_scan(drv, fast=False):
    print("\n" + "=" * 84)
    print("4. WHAT MOVES IF THE ROOT LIPID MOVES — all three tables, full ODE")
    print("=" * 84)
    if fast:
        print("   [--fast: skipped. Recorded result from the full run:]")
        for L, a, b, c in ((0.0050, 0.818, 0.382, 0.902), (0.0100, 0.598, 0.281, 0.783),
                           (0.0150, 0.471, 0.257, 0.719), (0.0200, 0.386, 0.267, 0.678),
                           (0.0247, 0.331, 0.288, 0.651), (0.0300, 0.290, 0.316, 0.629)):
            print(f"      L={L:.4f}   li2019 {a:.3f}   liu2023 {b:.3f}   ge2017 {c:.3f}")
        return
    orig = ND.TRAPP1994_LIPID_FW["root"]
    briggs_L = 10.0 ** ND.RCF_INTERCEPT / ND.LIPID_OCTANOL_A
    print(f"   {'root L':>8}{'li2019 n=29':>14}{'liu2023 n=14':>14}{'ge2017 n=6':>13}")
    try:
        for L in (0.0050, 0.0100, 0.0150, 0.0200, briggs_L, 0.0300):
            ND.TRAPP1994_LIPID_FW["root"] = L
            vals = [V.compare_to_obs(p, drv, quiet=True)
                    for p in (LI2019, LIU2023, GE2017)]
            tag = "  <- shipped" if abs(L - orig) < 1e-9 else (
                "  <- Briggs anchor" if abs(L - briggs_L) < 1e-9 else "")
            print(f"   {L:>8.4f}{vals[0]:>14.3f}{vals[1]:>14.3f}{vals[2]:>13.3f}{tag}")
    finally:
        ND.TRAPP1994_LIPID_FW["root"] = orig
    print("\n   Li 2019 and Ge 2017 improve monotonically toward the Briggs anchor;")
    print("   Liu 2023 -- the rice table -- has its optimum near L = 0.015 and gets")
    print("   worse beyond it. No single value is best for all three.")


def main(fast=False):
    drv = V.drivers()
    print("LI 2019 ROOT-PARTITION A-PRIORI TEST + THE ANCHOR DIAGNOSIS")
    print(f"forcings: measured Q_TP, growth_rice biomass, season {V.SEASON:.0f} d\n")
    rmse = section1_apriori(drv)
    section2_lipid_table()
    section3_anchor()
    section4_scan(drv, fast=fast)
    print("\n" + "=" * 84)
    print("VERDICT")
    print("=" * 84)
    print(f"   a-priori log10 RMSE on 29 out-of-sample rows = {rmse:.3f}, against 0.281")
    print("   for Liu 2023 (n=14, rice, narrower Kow range). The error is NOT scatter:")
    print("   every one of the 11 species is biased the same way, low, so this is one")
    print("   offset rather than eleven disagreements.")
    print("   That offset is traced in section 3 to the model's own composition, and")
    print("   section 4 shows correcting it helps two tables and hurts the rice one.")
    print("   DEFAULT UNCHANGED. Raising L to fit would turn the neutral path's one")
    print("   genuine claim -- that nothing is fitted -- into a fitted result.")
    print("   The physical reading is that the excess sorption Briggs' coefficient")
    print("   carries is NOT lipid: it is the non-lipid solid phase (cell wall,")
    print("   lignin) that this repo's PFAS side models explicitly as f_cw*K_cw and")
    print("   documents as GAP A, and that the neutral composition sets to zero.")
    print("   Closing it needs a measured neutral-organic cell-wall coefficient,")
    print("   not a larger lipid fraction.")
    return rmse


if __name__ == "__main__":
    main(fast="--fast" in sys.argv)
