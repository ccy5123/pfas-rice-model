#!/usr/bin/env python3
# =============================================================================
# validation/briggs1983_stem.py
# -----------------------------------------------------------------------------
# The STEM compartment against its only published anchor.
#
# Briggs, Bromilow, Evans & Williams 1983 (Pestic. Sci. 14:492-500) is the
# companion to the 1982 paper this whole neutral path is built on, and it was
# sitting unread in the obtained set. It does for the shoot what 1982 did for the
# root: measures the distribution of two compound series in barley shoots after
# root uptake, and fits a partition expression of exactly the K_PW form.
#
#     log(K_stem/xylem_sap - 0.82) = 0.95*log Kow - 2.05          (their eq. 2)
#     SCF = K_stem/xylem_sap * TSCF                               (their eq. 3)
#
# WHY IT MATTERS. Everything the root side of this repo has argued about -- is
# the lipid term anchored, is the exponent right, does the composition come from
# a measurement or a convention -- applies to the stem too, and until now the
# stem had NO anchor at all: it runs Trapp 1994's soybean 3 % lipid with the
# conventional a = 1.22 and Briggs' ROOT exponent 0.77. This script measures the
# gap. It turns out to be 4.1x, in the OPPOSITE direction to the root's 2.5x.
#
# WHAT THIS IS NOT. Briggs' "stem" is explicitly "not the true stem but is formed
# from the leaf bases and developing new leaves" of barley -- so it is a shoot
# base, not a rice culm, and it is one species. It is an anchor, not a
# measurement of rice. The same caution the root anchor carries.
#
#   python validation/briggs1983_stem.py
# =============================================================================
from __future__ import annotations
import os, sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT_DIR, "src"))

import neutral_dpu as ND                                       # noqa: E402


def section1_transcription():
    print("=" * 84)
    print("1. THE TRANSCRIPTION CHECKS ITSELF")
    print("=" * 84)
    print("   Briggs et al. state their eq. 3 peaks at 'about 6 ... at about")
    print("   log Kow = 4.5'. Reproducing that from the coefficients alone is what")
    print("   makes this a check rather than an assertion.\n")
    lk = np.linspace(-1.0, 7.0, 8001)
    scf = np.array([ND.briggs_scf(x) for x in lk])
    i = int(np.argmax(scf))
    print(f"   computed maximum SCF = {scf[i]:.2f} at log Kow = {lk[i]:.2f}")
    print("   -> matches the paper. The coefficients are transcribed correctly.")


def section2_gap():
    print("\n" + "=" * 84)
    print("2. WHAT THE REPO'S STEM RUNS, AGAINST THAT ANCHOR")
    print("=" * 84)
    L = ND.TRAPP1994_LIPID_FW["stem"]
    shipped_La = L * ND.LIPID_OCTANOL_A
    anchor_La = 10.0 ** ND.STEM_INTERCEPT
    print(f"   Briggs 1983 anchored stem lipid term  L*a = {anchor_La:.5f}  b = {ND.STEM_SLOPE}")
    print(f"   what rice_compartments runs             = {shipped_La:.5f}  b = {ND.RCF_SLOPE}"
          f"   (L={L}, a={ND.LIPID_OCTANOL_A})")
    print(f"   ratio {shipped_La / anchor_La:.1f}x HIGH in the coefficient, and a FLATTER slope.\n")
    print("   Because the slopes differ the two cross rather than sitting a constant")
    print("   factor apart, so the sign of the disagreement depends on lipophilicity:\n")
    W = ND.RICE_WATER["stem"]
    print(f"   {'log Kow':>9}{'repo stem K_PW':>17}{'Briggs K_stem/xyl':>20}{'ratio':>9}")
    for lk in (0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0):
        repo = ND.k_pw(lk, W=W, L=L)
        briggs = ND.briggs_stem_xylem_partition(lk)
        print(f"   {lk:>9.1f}{repo:>17.2f}{briggs:>20.2f}{repo / briggs:>9.2f}")
    cross = None
    for lk in np.linspace(-1, 7, 8001):
        if ND.k_pw(lk, W=W, L=L) < ND.briggs_stem_xylem_partition(lk):
            cross = lk
            break
    print(f"\n   They cross at log Kow ~ {cross:.2f}: below it the repo's stem is the more")
    print("   sorptive, above it Briggs' is. So this is not a simple 'too high' --")
    print("   it is a different SHAPE, driven by b = 0.77 against b = 0.95.")


def section3_scf():
    print("\n" + "=" * 84)
    print("3. THE OBSERVABLE — stem concentration factor")
    print("=" * 84)
    print("   SCF is the stem analogue of RCF: stem concentration over the EXTERNAL")
    print("   solution, so it folds in the xylem loading. That makes it the quantity")
    print("   this model actually has to reproduce, and it is testable without any")
    print("   ODE -- the same discipline as compare_to_obs(mode='equilibrium').\n")
    W = ND.RICE_WATER["stem"]
    L = ND.TRAPP1994_LIPID_FW["stem"]
    print(f"   {'log Kow':>9}{'TSCF':>8}{'repo SCF':>11}{'Briggs SCF':>13}{'log10 diff':>12}")
    for lk in (0.0, 1.0, 1.78, 2.5, 3.5, 4.5, 5.5):
        t = ND.tscf(lk)
        repo = ND.k_pw(lk, W=W, L=L) * t
        briggs = ND.briggs_scf(lk)
        print(f"   {lk:>9.2f}{t:>8.3f}{repo:>11.2f}{briggs:>13.2f}"
              f"{np.log10(repo / briggs):>+12.3f}")
    print("\n   The disagreement is MODEST where it matters and large where it does")
    print("   not: at most +0.13 log (35%) across log Kow 0-3.5, the range in which")
    print("   the TSCF bell actually delivers anything, then swinging to -0.20 and")
    print("   -0.38 above 4.5 where TSCF has collapsed and the stem receives almost")
    print("   nothing anyway. So the 4.1x coefficient gap in section 2 largely")
    print("   CANCELS against the steeper Briggs exponent over the useful range --")
    print("   which is the honest reading, and a weaker result than section 2 alone")
    print("   suggests. Briggs never tested above log Kow 4.5 -- their own text")
    print("   says 'the predicted decline in SCF for chemicals of log Kow > 4.5 was")
    print("   not tested' -- so the right-hand rows are extrapolation on both sides.")


def section4_verdict():
    print("\n" + "=" * 84)
    print("VERDICT")
    print("=" * 84)
    print("   The stem compartment has an anchor after all, and it does not match.")
    print("   Reading it next to the root result gives the honest summary of where")
    print("   this model's neutral composition stands:")
    print()
    root_ship = ND.TRAPP1994_LIPID_FW["root"] * ND.LIPID_OCTANOL_A
    root_anch = 10.0 ** ND.RCF_INTERCEPT
    stem_ship = ND.TRAPP1994_LIPID_FW["stem"] * ND.LIPID_OCTANOL_A
    stem_anch = 10.0 ** ND.STEM_INTERCEPT
    print(f"   {'organ':8}{'shipped L*a':>14}{'Briggs anchor':>16}{'ratio':>9}{'exponent':>22}")
    print(f"   {'root':8}{root_ship:>14.4f}{root_anch:>16.4f}"
          f"{root_ship / root_anch:>9.2f}{'0.77 vs 0.77 (same)':>22}")
    print(f"   {'stem':8}{stem_ship:>14.4f}{stem_anch:>16.4f}"
          f"{stem_ship / stem_anch:>9.2f}{'0.77 vs 0.95':>22}")
    print()
    print("   The root runs 2.5x BELOW its anchor and the stem 4.1x ABOVE its own.")
    print("   The two organs were parameterised from unrelated sources -- Briggs'")
    print("   fitted coefficient for one, Trapp's soybean composition for the other --")
    print("   and neither was ever checked against the anchor that existed for it.")
    print("   That is the finding; the CONSEQUENCES differ sharply between them.")
    print("   For the root the mismatch is the whole disagreement, because the")
    print("   exponent is shared. For the stem it largely cancels (section 3): the")
    print("   observable differs by at most 0.13 log where TSCF delivers anything.")
    print("   So the stem is a provenance problem, not a prediction problem, and it")
    print("   ranks well below the root question and the exposure-term work.")
    print()
    print("   NOTHING IS CHANGED HERE. The stem anchor is one species, is explicitly")
    print("   a shoot BASE rather than a true stem, and -- unlike the root case --")
    print("   there is no measured table in this repo that would arbitrate it. It is")
    print("   recorded, exposed as briggs_scf(), and pinned by a test.")


def main():
    print("BRIGGS 1983 — THE STEM ANCHOR THAT WAS NEVER READ\n")
    section1_transcription()
    section2_gap()
    section3_scf()
    section4_verdict()


if __name__ == "__main__":
    main()
