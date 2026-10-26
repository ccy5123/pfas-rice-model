#!/usr/bin/env python3
# =============================================================================
# validation/li2019_soil_table.py
# -----------------------------------------------------------------------------
# The 376-row SOIL companion to Li 2019's hydroponic table, and the three things
# it settles that the hydroponic one could not.
#
# WHY IT MATTERS. Li's hydroponic Table S1 and this soil Table S2 -- same authors, same
# fresh-weight convention, twelve times larger, 13 crops -- give the model's root partition
# opposite signs. Composition (PREREG 5.10): both run on Li's own basis, water 0.90, and this
# table on each row's own f_lip (Li's model input, so section 2 is a consistency check inside
# their framework). Sections:
#
#   1. the headline comparison: the default (each crop's own lipid), the rice root's
#      sourced lipid restated at 0.90, and the Briggs anchor;
#   2. the per-crop lipid test;
#   3. where the over-prediction comes from -- splitting on whether Li's K_om was measured
#      or estimated, the step that converts soil concentration into the water concentration
#      the endpoint needs.
#
# WHAT THIS TABLE IS NOT. `value` is DERIVED (soil concentration / K_om), not measured
# against a known solution, so it is weaker per row than a hydroponic RCF.
#
#   python validation/li2019_soil_table.py
# =============================================================================
from __future__ import annotations
import csv, os, sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT_DIR, "src"))

import neutral_dpu as ND
import neutral_composition as NC                                       # noqa: E402

OBS = os.path.join(ROOT_DIR, "data_obs", "neutral_obs_li2019_soil.csv")
HYDRO = os.path.join(ROOT_DIR, "data_obs", "neutral_obs_li2019_rcf.csv")


def load(path=OBS):
    with open(path, newline="") as f:
        rows = [r for r in csv.DictReader(x for x in f
                                          if not x.lstrip().startswith("#"))
                if r.get("compound")]
    for r in rows:
        for k in ("log_kow", "value", "f_lip", "f_om", "alpha_pt", "exposure_d"):
            if k in r:
                # exposure_d is blank for the 20 Zhang 2005 field rows
                r[k] = float(r[k]) if r[k] != "" else None
    return rows


LI_WATER = NC.value("li2019", "root", "water")


def bias(rows, lipid="own"):
    """log10(model K_PW) - log10(measured RCF) at Li's water 0.90, by the root lipid
    named: 'own' = each crop's own f_lip (what the model runs, PREREG 5.10), 'rice' = the
    rice root's sourced dry lipid restated at 0.90, 'anchor' = the Briggs lipid term."""
    if not rows:
        return 0, float("nan"), float("nan")
    W = LI_WATER
    L_rice = ND.RICE_LIPID_DW["root"] * (1.0 - W)
    e = []
    for r in rows:
        L = {"rice": L_rice,
             "anchor": ND.BRIGGS_ANCHORED_LIPID_FW["root"]}.get(lipid, r["f_lip"])
        e.append(np.log10(ND.k_pw(r["log_kow"], W=W, L=L)) - np.log10(r["value"]))
    e = np.asarray(e)
    return len(e), float(e.mean()), float(np.sqrt((e ** 2).mean()))


def _hydroponic_bias():
    with open(HYDRO, newline="") as f:
        rs = [r for r in csv.DictReader(x for x in f
                                        if not x.lstrip().startswith("#"))
              if r.get("subset") == "apriori"]
    e = []
    for r in rs:
        W, L = NC.spec("li2019_rcf", species=r["species"])
        e.append(np.log10(ND.k_pw(float(r["log_kow"]), W=W["root"], L=L["root"]))
                 - np.log10(float(r["value"])))
    return len(e), float(np.mean(e))


def section1(rows):
    print("=" * 84)
    print("1. THE SIGN FLIP — the same paper's soil table disagrees with its hydroponic one")
    print("=" * 84)
    nh, bh = _hydroponic_bias()
    n, b, e = bias(rows)
    print(f"   Li 2019 Table S1, HYDROPONIC   n={nh:3d}   mean log10 bias {bh:+.3f}")
    print(f"   Li 2019 Table S2, SOIL         n={n:3d}   mean log10 bias {b:+.3f}")
    print("\n   Same authors, same fresh-weight convention, same endpoint definition, and")
    print("   both on Li's water 0.90 (PREREG 5.10).\n")
    print(f"   {'root composition':34}{'n':>6}{'bias':>10}{'RMSE':>9}")
    for lab, key in (("each crop's own f_lip (default)", "own"),
                     ("rice root lipid at 0.90 water", "rice"),
                     ("Briggs anchor, L=0.0247", "anchor")):
        n, b, e = bias(rows, key)
        print(f"   {lab:34}{n:>6}{b:>+10.3f}{e:>9.3f}")


def section2(rows):
    print("\n" + "=" * 84)
    print("2. DOES EACH CROP'S OWN LIPID COLLAPSE THE OFFSET?")
    print("=" * 84)
    print("   13 crops whose measured root lipid spans 0.10-1.14 %, an 11x range,")
    print("   against one composition for all of them (the rice root's, at 0.90):\n")
    for lab, key in (("rice composition for all crops", "rice"),
                     ("each crop's OWN measured f_lip", "own")):
        n, b, e = bias(rows, key)
        print(f"   {lab:34}{n:>6}{b:>+10.3f}{e:>9.3f}")
    print("   Caveat, and it is not small: f_lip is Li's own model input, so this is")
    print("   a consistency check inside their framework, not independent validation.")
    print(f"\n   {'crop':18}{'n':>5}{'f_lip %':>9}{'rice bias':>12}{'own bias':>11}")
    by = {}
    for r in rows:
        by.setdefault(r["species"], []).append(r)
    for sp, rs in sorted(by.items(), key=lambda kv: -len(kv[1])):
        n, b1, _ = bias(rs, "rice")
        _, b2, _ = bias(rs, "own")
        print(f"   {sp:18}{n:>5}{rs[0]['f_lip'] * 100:>9.2f}{b1:>+12.3f}{b2:>+11.3f}")


def section3(rows):
    print("\n" + "=" * 84)
    print("3. THE SOIL-TO-WATER CONVERSION — measured vs estimated K_om")
    print("=" * 84)
    print("   `value` is not measured against a known solution: Li derive it from the")
    print("   measured SOIL concentration by dividing by K_om. Their Table S3 records")
    print("   which K_om was experimental and which came from a QSPR, so the cost of")
    print("   that step is directly measurable.\n")
    print(f"   {'K_om source':28}{'n':>6}{'bias':>12}{'RMSE':>9}")
    for src in ("experimental", "estimated", "unmatched"):
        rs = [r for r in rows if r["kom_source"] == src]
        n, b, e = bias(rs)
        if n:
            print(f"   {src:28}{n:>6}{b:>+12.3f}{e:>9.3f}")
    n, b, _ = bias([r for r in rows if r["kom_source"] == "experimental"])
    print(f"\n   On the {n} rows with a MEASURED K_om the mean bias is {b:+.3f}; compare it")
    print("   with the estimated-K_om rows above to see what the exposure conversion costs.")
    print("\n   The soil organic matter trend points the same way (the K_om step scales")
    print("   with f_om, so an error in it should show up as an f_om gradient):")
    print(f"   {'f_om':>14}{'n':>6}{'bias':>10}")
    for lo, hi in ((0.0, 0.01), (0.01, 0.02), (0.02, 0.04), (0.04, 1.0)):
        rs = [r for r in rows if lo <= r["f_om"] < hi]
        n, b, _ = bias(rs)
        if n:
            print(f"   {lo * 100:5.1f}-{hi * 100:5.1f} %{n:>6}{b:>+10.3f}")


def section4(rows):
    print("\n" + "=" * 84)
    print("4. WHY alpha_pt IS NOT APPLIED")
    print("=" * 84)
    apt = np.array([r["alpha_pt"] for r in rows])
    print(f"   Li's own quasi-equilibrium factor over these rows: median {np.median(apt):.3f}")
    print(f"   (range {apt.min():.3f}-{apt.max():.3f}). Their model is")
    print("   RCF_water = alpha_pt * [f_pw + f_ch*K_ch + f_lip*K_lip], so their")
    print("   equilibrium term is ~10x the measured RCF here.")
    print("\n   This model has NO alpha_pt and lands on the measured values anyway,")
    print("   because Briggs' exponent b = 0.77 is far flatter than Li's")
    print("   K_lip ~ Kow^1.03. At log Kow 5 the two lipid terms differ ~20x:")
    for lk in (3.0, 4.0, 5.0):
        briggs = ND.LIPID_OCTANOL_A * 10.0 ** (ND.RCF_SLOPE * lk)
        li = 1.27 * 10.0 ** (1.03 * lk)
        print(f"      log Kow {lk:.1f}   Briggs a*Kow^0.77 = {briggs:10.1f}   "
              f"Li K_lip = {li:12.1f}   ratio {li / briggs:5.1f}x")
    print("\n   So the two frameworks split the same data differently: what Li put in")
    print("   alpha_pt, Briggs put in the exponent. Dividing by alpha_pt here would")
    print("   double-count the attenuation.")


def main():
    rows = load()
    print(f"LI 2019 TABLE S2 — {len(rows)} SOIL-GROWN ROOT RCFs, "
          f"{len({r['species'] for r in rows})} CROPS\n")
    section1(rows)
    section2(rows)
    section3(rows)
    section4(rows)
    print("\n" + "=" * 84)
    print("VERDICT")
    print("=" * 84)
    nh, bh = _hydroponic_bias()
    n, b, _ = bias(rows)
    print(f"   Li's hydroponic half: mean log10 bias {bh:+.3f} (n={nh}); this soil half")
    print(f"   {b:+.3f} (n={n}), both on Li's own water and lipid (PREREG 5.10). The rows")
    print("   with a measured K_om (section 3) isolate the plant side from the soil side.")
    return bias(rows)[1]


if __name__ == "__main__":
    main()
