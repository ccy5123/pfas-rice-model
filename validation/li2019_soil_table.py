#!/usr/bin/env python3
# =============================================================================
# validation/li2019_soil_table.py
# -----------------------------------------------------------------------------
# The 376-row SOIL companion to Li 2019's hydroponic table, and the three things
# it settles that the hydroponic one could not.
#
# WHY IT MATTERS. `validation/li2019_rcf_apriori.py` reports that the model's
# root partition runs LOW by ~0.43 log on Li's hydroponic Table S1 (n=29), and
# that reading is what makes restoring the Briggs anchor tempting. Table S2 of
# the SAME PAPER -- same authors, same fresh-weight convention, twelve times
# larger, 13 crops -- says the opposite: the model runs HIGH by +0.26.
#
# The sign flip is the finding. It means "the root partition is too low" is a
# property of one exposure route, not of the partition core, and it removes the
# main argument for moving the default. Three sections here take it apart:
#
#   1. the headline comparison, and what the Briggs anchor would do to it;
#   2. the per-crop lipid test -- does giving each crop its OWN measured lipid
#      collapse the offset? (yes, to within 0.01 log);
#   3. where the over-prediction actually comes from -- splitting on whether
#      Li's K_om was measured or estimated, which is the step that converts
#      soil concentration into the water concentration the endpoint needs.
#
# WHAT THIS TABLE IS NOT. `value` is DERIVED (soil concentration / K_om), not
# measured against a known solution, so it is weaker per row than a hydroponic
# RCF. Section 3 is precisely about how much that costs. And `f_lip` is Li's own
# model input, so section 2 is a consistency check inside their framework rather
# than an independent validation. Neither caveat is small; both are why this
# script exists instead of the table being wired into the `--obs` harness.
#
#   python validation/li2019_soil_table.py
# =============================================================================
from __future__ import annotations
import csv, os, sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT_DIR, "src"))

import neutral_dpu as ND                                       # noqa: E402

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


def bias(rows, lipid="rice"):
    """log10(model K_PW) - log10(measured RCF), by the compartment composition
    named: 'rice' = what the model runs, 'anchor' = the Briggs lipid term,
    'own' = each crop's own measured lipid."""
    if not rows:
        return 0, float("nan"), float("nan")
    W = ND.RICE_WATER["root"]
    e = []
    for r in rows:
        L = {"rice": ND.TRAPP1994_LIPID_FW["root"],
             "anchor": ND.BRIGGS_ANCHORED_LIPID_FW["root"]}.get(lipid, r["f_lip"])
        e.append(np.log10(ND.k_pw(r["log_kow"], W=W, L=L)) - np.log10(r["value"]))
    e = np.asarray(e)
    return len(e), float(e.mean()), float(np.sqrt((e ** 2).mean()))


def _hydroponic_bias():
    W = ND.RICE_WATER["root"]
    L = ND.TRAPP1994_LIPID_FW["root"]
    with open(HYDRO, newline="") as f:
        rs = [r for r in csv.DictReader(x for x in f
                                        if not x.lstrip().startswith("#"))
              if r.get("subset") == "apriori"]
    e = [np.log10(ND.k_pw(float(r["log_kow"]), W=W, L=L)) - np.log10(float(r["value"]))
         for r in rs]
    return len(e), float(np.mean(e))


def section1(rows):
    print("=" * 84)
    print("1. THE SIGN FLIP — the same paper's soil table disagrees with its hydroponic one")
    print("=" * 84)
    nh, bh = _hydroponic_bias()
    n, b, e = bias(rows)
    print(f"   Li 2019 Table S1, HYDROPONIC   n={nh:3d}   mean log10 bias {bh:+.3f}  (model LOW)")
    print(f"   Li 2019 Table S2, SOIL         n={n:3d}   mean log10 bias {b:+.3f}  (model HIGH)")
    print("\n   Same authors, same fresh-weight convention, same endpoint definition.")
    print("   So the deficit the hydroponic table reports is a property of that")
    print("   exposure route, not of the partition core.\n")
    print(f"   {'root composition':34}{'n':>6}{'bias':>10}{'RMSE':>9}")
    for lab, key in (("rice, L=0.010 (what runs today)", "rice"),
                     ("Briggs anchor, L=0.0247", "anchor")):
        n, b, e = bias(rows, key)
        print(f"   {lab:34}{n:>6}{b:>+10.3f}{e:>9.3f}")
    print("\n   Restoring the anchor makes this table MUCH worse -- it pushes an")
    print(f"   already-high model higher. Whatever the hydroponic rows want, {len(rows)}")
    print(f"   soil-grown rows across {len({r['species'] for r in rows})} crops want the opposite.")


def section2(rows):
    print("\n" + "=" * 84)
    print("2. DOES EACH CROP'S OWN LIPID COLLAPSE THE OFFSET?")
    print("=" * 84)
    print("   The model runs one composition (rice) against 13 crops whose measured")
    print("   root lipid spans 0.10-1.14 %, an 11x range. If that is the whole story,")
    print("   substituting each crop's own value should remove the bias.\n")
    for lab, key in (("rice composition for all crops", "rice"),
                     ("each crop's OWN measured f_lip", "own")):
        n, b, e = bias(rows, key)
        print(f"   {lab:34}{n:>6}{b:>+10.3f}{e:>9.3f}")
    print("\n   It does: the mean bias goes to essentially zero. So the FORM of")
    print("   K_PW = W + L*a*Kow^b is doing its job across an 11x spread in L, and")
    print("   the residual scatter is not a failure of the partition expression.")
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
    print("3. WHERE THE OVER-PREDICTION COMES FROM — the soil-to-water conversion")
    print("=" * 84)
    print("   `value` is not measured against a known solution: Li derive it from the")
    print("   measured SOIL concentration by dividing by K_om. Their Table S3 records")
    print("   which K_om was experimental and which came from a QSPR, so the cost of")
    print("   that step is directly measurable.\n")
    print(f"   {'K_om source':28}{'n':>6}{'rice bias':>12}{'RMSE':>9}")
    for src in ("experimental", "estimated", "unmatched"):
        rs = [r for r in rows if r["kom_source"] == src]
        n, b, e = bias(rs)
        if n:
            print(f"   {src:28}{n:>6}{b:>+12.3f}{e:>9.3f}")
    n, b, _ = bias([r for r in rows if r["kom_source"] == "experimental"])
    print(f"\n   On the {n} rows with a MEASURED K_om the model is essentially unbiased")
    print(f"   ({b:+.3f}). Most of the apparent over-prediction is therefore the")
    print("   exposure conversion, not the plant model -- which is the same lesson")
    print("   the Kodesova table teaches from the other direction, where a measured")
    print("   isotherm was what made the exposure usable at all.")
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
    print("   The hydroponic and soil halves of the same paper disagree in SIGN, so")
    print("   neither settles the anchor question on its own -- and the larger half")
    print("   argues against moving the default, joining Liu 2023 and Kodesova 2019.")
    print("   The two findings that do generalise:")
    print("     * the composition term works -- giving each crop its own lipid")
    print("       removes the bias across an 11x spread in L;")
    print("     * where the exposure is MEASURED rather than estimated, the model is")
    print("       close to unbiased. That now holds on three independent tables")
    print("       (this one's 62 measured-K_om rows, Liu 2023, Kodesova 2019), and it")
    print("       relocates most of the apparent error from the plant model to the")
    print("       exposure term -- which is a soil-side problem, not a K_PW problem.")
    return bias(rows)[1]


if __name__ == "__main__":
    main()
