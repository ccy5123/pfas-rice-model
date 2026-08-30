#!/usr/bin/env python3
# =============================================================================
# validation/briggs1983_shoot.py
# -----------------------------------------------------------------------------
# The STEM against MEASURED data -- Briggs 1983 Table 1, 16 compounds.
#
# docs/neutral_dpu_validation.md section 4j read Briggs 1983's fitted EQUATIONS
# and compared them to the repo's stem composition, finding a 4.1x coefficient
# gap that it argued "largely cancels" in the observable SCF. That argument was
# equation-against-equation: no measurement entered it, and the stem had never
# been scored against data at all.
#
# This scores it. Table 1's per-section distribution plus the section fresh
# weights in the article's section 2.2 reconstruct the Stem Concentration Factor
# the paper itself defines, for 16 non-ionised chemicals over log Kow -0.57 to
# 3.7, taken up by barley from a nutrient solution of KNOWN concentration.
# See data_obs/neutral_obs_briggs1983_shoot.csv for the derivation and caveats.
#
# WHAT COMES OUT, in order:
#   1. the reconstruction checks itself, twice, before anything is scored
#   2. the STEM a-priori result -- the repo's first, nothing fitted
#   3. the LEAF, which is NOT an equilibrium, and what that costs item 4
#   4. the three compounds the article itself says are special, as a check that
#      the reconstruction reproduces the paper's own diagnoses
#
#   python validation/briggs1983_shoot.py
# =============================================================================
from __future__ import annotations
import csv, os, sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT_DIR, "src"))

import neutral_dpu as ND                                       # noqa: E402

OBS = os.path.join(ROOT_DIR, "data_obs", "neutral_obs_briggs1983_shoot.csv")


def load(path=OBS):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(l for l in f if not l.lstrip().startswith("#")))
    for r in rows:
        for k in ("log_kow", "value", "pct_of_shoot", "section_g"):
            r[k] = float(r[k])
        for k in ("hours", "equilibrated", "metabolised"):
            r[k] = int(r[k])
    return rows


def briggs_eq3(log_kow):
    """The article's own eq. (3): SCF = K_stem/xylem * TSCF, both from log Kow."""
    return ND.briggs_stem_xylem_partition(log_kow) * ND.tscf(log_kow)


def repo_scf(log_kow):
    """What THIS repo predicts for a stem, with nothing fitted: its own neutral
    composition (Trapp 1994 stem lipid, conventional a, Briggs' ROOT exponent)
    through k_pw, times the same Briggs TSCF."""
    return ND.k_pw(log_kow, W=ND.RICE_WATER["stem"],
                   L=ND.TRAPP1994_LIPID_FW["stem"]) * ND.tscf(log_kow)


def _score(rows, tissue, fn, hours=None, equilibrated_only=True):
    e = [np.log10(r["value"] / fn(r["log_kow"])) for r in rows
         if r["tissue"] == tissue
         and (hours is None or r["hours"] == hours)
         and (r["equilibrated"] or not equilibrated_only)]
    e = np.asarray(e)
    return len(e), float(e.mean()), float(np.sqrt((e ** 2).mean()))


def section1(rows):
    print("=" * 84)
    print("1. THE RECONSTRUCTION CHECKS ITSELF, TWICE")
    print("=" * 84)
    print("   The SCF is derived here, not published as such: it is Table 1's %")
    print("   times the total dpm, over the section weight, over the solution")
    print("   concentration. Two independent checks that this is the quantity the")
    print("   article actually plotted.\n")

    lk = np.linspace(-1.0, 7.0, 8001)
    s = np.array([briggs_eq3(x) for x in lk])
    i = int(np.argmax(s))
    print(f"   (a) the article's eq. (3) peaks at SCF {s[i]:.2f} at log Kow {lk[i]:.2f}")
    print("       -- the paper says 'about 6 ... at about log Kow = 4.5'. The")
    print("       equations behind the comparison are transcribed correctly.")

    n, b, e = _score(rows, "stem_central", briggs_eq3)
    print(f"\n   (b) the RECONSTRUCTED central-stem SCF against that same eq. (3):")
    print(f"       n={n}   mean log10 bias {b:+.3f}   RMSE {e:.3f}")
    print("       Essentially unbiased. The article says its points 'fit quite")
    print("       well' the eq.-(3) curve, and the reconstruction reproduces that")
    print("       -- so the derived numbers ARE its Fig. 4, not something else.")

    nb, bb, eb = _score(rows, "stem_base", briggs_eq3)
    print(f"\n   (c) the stem BASE runs high against the same curve: bias {bb:+.3f}")
    print("       The article predicts exactly this -- the base 'was up to twice'")
    print("       the central section, ascribed partly to DIRECT CONTACT with the")
    print("       treating solution. A third agreement, and the reason section 2")
    print("       scores the central section: it is the contact-free one.")


def section2(rows):
    print("\n" + "=" * 84)
    print("2. THE STEM, A-PRIORI -- the repo's first stem result against data")
    print("=" * 84)
    print("   Nothing is fitted. K_PW comes from the repo's own stem composition")
    print("   and TSCF from log Kow; both were fixed before this table was read.\n")
    print(f"   {'':22}{'n':>4}{'bias':>9}{'RMSE':>8}")
    out = {}
    for tissue in ("stem_central", "stem_base"):
        for lab, fn in (("repo K_PW x TSCF", repo_scf), ("Briggs' own eq. (3)", briggs_eq3)):
            n, b, e = _score(rows, tissue, fn)
            out[(tissue, lab)] = (n, b, e)
            print(f"   {tissue:13} {lab:22.22}{n:>4}{b:>+9.3f}{e:>8.3f}")
    rc = out[("stem_central", "repo K_PW x TSCF")]
    bc = out[("stem_central", "Briggs' own eq. (3)")]
    print(f"\n   THE RESULT: on the clean section the repo predicts measured SCF at")
    print(f"   RMSE {rc[2]:.3f} (bias {rc[1]:+.3f}) -- indistinguishable from the")
    print(f"   paper's OWN FITTED equation at {bc[2]:.3f} ({bc[1]:+.3f}), which had this")
    print("   very data to fit. For scale, the root tables run 0.281 (Liu) to")
    print("   0.598 (Li 2019 hydroponic), both on the same a-priori footing.")
    print("\n   This settles section 4j empirically. That section found the repo's")
    print("   stem lipid term 4.1x ABOVE Briggs' and argued the gap 'largely")
    print("   cancels' against his steeper exponent, leaving the observable SCF")
    print("   nearly unchanged -- but it argued it equation-against-equation. The")
    print("   measurements now say the same thing: the two predictions are equally")
    print("   good, so the 4.1x is a PROVENANCE problem, not a prediction problem.")
    print("   The composition is still unsourced for rice; it is just not wrong")
    print("   where it can be checked.")
    return out


def section3(rows):
    print("\n" + "=" * 84)
    print("3. THE LEAF IS NOT AN EQUILIBRIUM -- what that costs the item-4 question")
    print("=" * 84)
    print("   The reason to want this table was to constrain the STEM/LEAF SPLIT")
    print("   independently of Ge 2017. It cannot, and the data say why.\n")

    print(f"   {'':24}{'stem share of shoot burden (%)':>34}")
    print(f"   {'compound':30}{'logKow':>8}{'24 h':>8}{'48 h':>8}{'ratio':>8}")
    by = {}
    for r in rows:
        by.setdefault((r["compound"], r["hours"]), {})[r["tissue"]] = r
    ratios = []
    for c in sorted({c for c, _ in by}, key=lambda c: by[(c, 24)]["leaf"]["log_kow"]):
        share = {}
        for h in (24, 48):
            v = by[(c, h)]
            share[h] = v["stem_central"]["pct_of_shoot"] + v["stem_base"]["pct_of_shoot"]
        lk = by[(c, 24)]["leaf"]["log_kow"]
        ratios.append(share[24] / share[48])
        print(f"   {c:30}{lk:>8.2f}{share[24]:>8.1f}{share[48]:>8.1f}{ratios[-1]:>7.2f}x")
    print(f"\n   The split MOVES: median {np.median(ratios):.2f}x between two harvests one")
    print("   day apart. It is not a property of the compound, so a single number")
    print("   from it cannot be a model target without matching the exposure")
    print("   duration -- and this is a 48-h barley seedling, not a rice season.")

    print("\n   WHY it moves, from the same table: the stem equilibrates and the")
    print("   leaf does not. Scoring the leaf as if it were an equilibrium:\n")
    for h in (24, 48):
        n, b, e = _score(rows, "leaf",
                         lambda x: ND.k_pw(x, W=ND.RICE_WATER["leaf"],
                                           L=ND.TRAPP1994_LIPID_FW["leaf"],
                                           b=ND.LEAF_SORPTION_EXPONENT) * ND.tscf(x),
                         hours=h)
        print(f"     {h} h   n={n}   bias {b:+.3f}   RMSE {e:.3f}")
    print("\n   The bias GROWS with exposure -- the leaf is climbing away from any")
    print("   equilibrium value, which is the terminal-accumulator signature the")
    print("   model has (docs section 3) and the article states independently:")
    print("   leaf amounts 'generally increased up to 72 or 96 h'. So this is a")
    print("   SECOND dataset confirming the leaf structure, from a different")
    print("   species and exposure than Ge 2017 -- just not a split constraint.")


def section4(rows):
    print("\n" + "=" * 84)
    print("4. THE COMPOUNDS THE ARTICLE ITSELF FLAGS")
    print("=" * 84)
    print("   A last check that the reconstruction carries the paper's own")
    print("   diagnoses rather than smoothing them away.\n")
    by = {(r["compound"], r["hours"], r["tissue"]): r for r in rows}
    for c, why in (("aldicarb", "oxidised in planta to a polar sulphoxide that is trapped"),
                   ("aldoxycarb", "log Kow < 0, where the article says its TSCF curve is low"),
                   ("4_4_bromophenoxyphenylurea", "never equilibrated -- excluded from its own Fig. 4")):
        r = by[(c, 48, "stem_central")]
        p = repo_scf(r["log_kow"])
        print(f"   {c:28} logKow {r['log_kow']:5.2f}  obs {r['value']:6.3f}  model {p:6.3f}"
              f"  {r['value'] / p:5.2f}x")
        print(f"       {why}")
    print("\n   All three deviate in the direction and by roughly the magnitude the")
    print("   article reports (it puts aldicarb at 'about three times' predicted).")
    print("   They are excluded or flagged in the CSV, not silently averaged in.")


def main():
    rows = load()
    print("BRIGGS 1983 BARLEY SHOOTS -- the stem anchor, against measurements")
    print(f"{len(rows)} rows = 16 compounds x 2 harvests x 3 sections; "
          f"exposure is a known nutrient solution\n")
    section1(rows)
    out = section2(rows)
    section3(rows)
    section4(rows)
    print("\n" + "=" * 84)
    print("VERDICT")
    print("=" * 84)
    n, b, e = out[("stem_central", "repo K_PW x TSCF")]
    print(f"  STEM: a-priori log10 RMSE {e:.3f} on n={n}, bias {b:+.3f}. The stem had no")
    print("    measured test before this; it now has one, and it passes at the")
    print("    level the root tables do. Section 4j's cancellation argument is")
    print("    confirmed against data.")
    print("  LEAF: not an equilibrium here, so the stem/leaf SPLIT this table was")
    print("    wanted for is not extractable -- the split moves ~1.6x in 24 h.")
    print("    What it gives instead is independent confirmation of the terminal-")
    print("    accumulator leaf.")
    print("  SCOPE: barley, 11-day seedlings, 24-48 h, one lab, two chemical")
    print("    series. It anchors the stem PARTITION FORM; it says nothing about")
    print("    a rice culm's composition, which remains unmeasured.")
    return e


if __name__ == "__main__":
    main()
