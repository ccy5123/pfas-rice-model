#!/usr/bin/env python3
# =============================================================================
# validation/kodesova2019_carbamazepine.py
# -----------------------------------------------------------------------------
# Acquisition_Queue.csv row A4, closed. Carbamazepine root uptake from three
# soils by four plants (Kodesova et al. 2019), with the exposure pinned by the
# paper's OWN measured Freundlich isotherms.
#
# WHY THIS ONE MATTERS MORE THAN ITS SIZE SUGGESTS. It was requested to settle a
# specific open question -- Brunetti 2021's calibrated green-pea root partition
# K_RW = 13.3 L/kg fw for carbamazepine, against a Briggs K_PW of ~1.6 -- which
# `docs/neutral_dpu_validation.md` section 5 lists as one of four sightings that
# the root partition core is too LOW. This is the direct measurement that
# question was asking for, on the same compound.
#
# Composition (PREREG 5.10): each plant's own root water and lipid
# (params/neutral_composition_sourced.csv; radish from Ievinsh 2023, the others at Li
# 2019's 0.90 convention), and the dry-weight RCFs are converted with that same water.
#
# Sections 6 and 7 use the LEAF half of the same table, which the root comparison
# leaves untouched. Section 6 tests TRANSLOCATION with the exposure cancelled
# (a leaf/root ratio needs no exposure at all, so it does not inherit the
# isotherm assumption section 1 defends). Section 7 is the more valuable: the
# paper measured carbamazepine's METABOLITES alongside the parent, so in-planta
# transformation stops being the free parameter every half-life statement in
# this repo has had to treat it as.
#
# THE EXPOSURE, which is what makes the dataset usable. Three measured quantities
# and no model in between:
#     C_root  SI Table S2, ng/g dry weight (basis stated in the article)
#     C_soil  SI Table S4, ng/g dry -- same pot, same harvest
#     KF, n   article Table 1, measured Freundlich isotherms, per soil
# then  c = (C_soil/KF)^n  from  S = KF*c^(1/n),  and  RCF = C_root / c.
# No mass balance, no pot geometry, no dissipation model. Section 1 argues the
# isotherm convention rather than assuming it, because that is the one input that
# could flip the conclusion.
#
#   python validation/kodesova2019_carbamazepine.py
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
import neutral_composition as NC                               # noqa: E402

TISSUES = ("root", "stem", "leaf", "grain")

OBS = os.path.join(ROOT_DIR, "data_obs", "neutral_obs_kodesova2019.csv")
LIU = os.path.join(ROOT_DIR, "data_obs", "neutral_obs_liu2023.csv")
LI2019 = os.path.join(ROOT_DIR, "data_obs", "neutral_obs_li2019_rcf.csv")

LOG_KOW = 2.25                    # article Table 1, sourced to Jones et al. 2002
KF = {"HCh": 3.86, "HCa": 2.97, "AE": 0.71}          # article Table 1, per soil
COX = {"HCh": 0.0174, "HCa": 0.0157, "AE": 0.0046}   # article Table 2
BRUNETTI_K_RW = 13.3              # Brunetti 2021 Table 1, calibrated green pea


def rows():
    with open(OBS, newline="") as f:
        return [r for r in csv.DictReader(x for x in f
                                          if not x.lstrip().startswith("#"))
                if r.get("compound")]


def species_comp(sp, corner=NC.MID, lipid_source="measured"):
    """(W, L_fw) of one plant's root (PREREG 5.10); "briggs_anchor" swaps in Briggs'
    fitted lipid term."""
    W, L = NC.spec("kodesova2019", species=sp, corner=corner)
    Lr = ND.LIPID_SOURCES["briggs_anchor"]["root"] if lipid_source != "measured" else L["root"]
    return W["root"], Lr


def rcf_fw(rs=None, theta=None):
    """Measured RCFs on fresh weight: each dry-weight value times (1 - its plant's root
    water), or (1 - theta) when a water is given."""
    rs = rows() if rs is None else rs
    return np.array([float(r["value"]) * (1.0 - (species_comp(r["species"])[0]
                                                if theta is None else theta)) for r in rs])


def kpw_rows(rs=None, lipid_source="measured"):
    rs = rows() if rs is None else rs
    return np.array([ND.k_pw(LOG_KOW, *species_comp(r["species"], lipid_source=lipid_source))
                     for r in rs])


def section1_isotherm():
    print("=" * 84)
    print("1. THE EXPOSURE — and why the isotherm convention is argued, not assumed")
    print("=" * 84)
    print("   At c = 1 mg/L the Freundlich isotherm gives S = KF, so KF is numerically")
    print("   the distribution ratio there. Dividing by each soil's measured organic")
    print("   carbon is therefore a direct test of the unit reading:\n")
    print(f"   {'soil':6}{'KF':>7}{'Cox %':>8}{'Koc':>8}")
    for s in ("HCh", "HCa", "AE"):
        print(f"   {s:6}{KF[s]:>7.2f}{COX[s] * 100:>8.2f}{KF[s] / COX[s]:>8.0f}")
    print("\n   Three soils spanning 3.8x in organic carbon give Koc 154-222, consistent")
    print("   and inside carbamazepine's literature band (~100-500). Reading KF in")
    print("   g/cm3 instead would give Koc ~1600 in every soil, outside it.")
    print("   The EXPONENT is settled separately: Table 1's unit string carries 1/n")
    print("   with n = 1.13, so S = KF*c^0.885 -- concave, the shape soils actually")
    print("   have. The alternative exponent 1.13 would make sorption MORE favourable")
    print("   as concentration rises.")
    cw = np.array([float(r["pore_water_mgL"]) for r in rows()])
    print(f"\n   Derived pore water: {cw.min():.3f}-{cw.max():.3f} mg/L, against an applied")
    print("   irrigation solution of ~0.5-1.0 mg/L and a sorption study calibrated over")
    print("   0.5-10 mg/L. Same scale as both -- which a wrong unit reading would not be.")


def section2_apriori(drv):
    print("\n" + "=" * 84)
    print("2. A-PRIORI PREDICTION — nothing fitted, through the full ODE")
    print("=" * 84)
    rmse = V.compare_to_obs(OBS, drv, group_by="species")
    print("\n   For scale, the repo's other a-priori root tables on the same path:")
    for lab, path, kw in (("Liu 2023, rice, n=14", LIU, {}),
                          ("Li 2019, 11 species, n=29", LI2019, {})):
        print(f"      {lab:28s} {V.compare_to_obs(path, drv, quiet=True, **kw):.3f}")
    print(f"      {'Kodesova 2019, n=21':28s} {rmse:.3f}   <- this table")
    return rmse


def section3_anchor(drv):
    """The sourced composition against Briggs' fitted root lipid term."""
    print("\n" + "=" * 84)
    print("3. SOURCED COMPOSITION vs THE BRIGGS ANCHOR")
    print("=" * 84)
    out = {lab: [V.compare_to_obs(p, drv, quiet=True, lipid_source=src)
                 for p in (OBS, LIU, LI2019)]
           for lab, src in (("sourced (PREREG 5.10)", "measured"),
                            ("Briggs anchor", "briggs_anchor"))}
    print(f"   {'root lipid':22}{'Kodesova n=21':>15}{'Liu rice n=14':>15}{'Li2019 n=29':>14}")
    for lab, v in out.items():
        print(f"   {lab:22}{v[0]:>15.3f}{v[1]:>15.3f}{v[2]:>14.3f}")
    obs = rcf_fw()
    print("\n   Driver-free cross-check (equilibrium K_PW per plant, no ODE):")
    print(f"      measured RCF_fw   median {np.median(obs):.2f}   range "
          f"{obs.min():.2f}-{obs.max():.2f}")
    for lab, src in (("sourced K_PW", "measured"), ("anchored K_PW", "briggs_anchor")):
        k = kpw_rows(lipid_source=src)
        b = float(np.mean(np.log10(k) - np.log10(obs)))
        print(f"      {lab:16} median {np.median(k):.2f}   mean log10 bias {b:+.3f}")


def section4_brunetti():
    print("\n" + "=" * 84)
    print("4. THE BRUNETTI QUESTION, WHICH IS WHAT THIS PAPER WAS FETCHED FOR")
    print("=" * 84)
    obs = rcf_fw()
    kpw = kpw_rows()
    print(f"   Brunetti 2021, calibrated green-pea root K_RW   {BRUNETTI_K_RW:.1f} L/kg fw")
    print(f"   Kodesova 2019, MEASURED, 4 plants x 3 soils     {np.median(obs):.2f} "
          f"(range {obs.min():.2f}-{obs.max():.2f})")
    print(f"   Briggs K_PW on each plant's sourced composition {np.median(kpw):.2f} (median)")
    print(f"\n   Brunetti's value is {BRUNETTI_K_RW / np.median(obs):.0f}x the direct measurement of the")
    print(f"   same compound; the a-priori K_PW is {np.median(kpw) / np.median(obs):.2f}x it. A calibrated")
    print("   posterior can absorb whatever else the model was missing; the partition core")
    print("   here has nothing fitted.")
    print("\n   Note the species caveat both ways: Brunetti is green pea, this is four")
    print("   leafy/root vegetables, and neither is rice.")


def section5_sensitivity():
    print("\n" + "=" * 84)
    print("5. WHAT THE CONCLUSION RESTS ON")
    print("=" * 84)
    rs = rows()
    print("   The root water converts the dry-weight RCF and sets K_PW (with the lipid")
    print("   restated at it), so the sourced range ends are run, not asserted:")
    print(f"   {'corner (water, lipid)':24}{'equilibrium bias':>18}")
    for c in ((NC.MID,) + NC.CORNERS):
        obs = np.array([float(r["value"]) * (1.0 - species_comp(r["species"], c)[0])
                        for r in rs])
        k = np.array([ND.k_pw(LOG_KOW, *species_comp(r["species"], c)) for r in rs])
        print(f"   {str(c):24}{float(np.mean(np.log10(k) - np.log10(obs))):>+18.3f}")
    print("\n   What COULD flip the verdict: the isotherm unit reading in section 1.")
    print("   Reading KF in g/cm3 lowers the derived pore water ~6x and raises every")
    print("   RCF ~6x. That reading is rejected on the Koc evidence, but it is the")
    print("   assumption to attack if this result is ever doubted.")


# SI Table S2, carbamazepine and its four measured metabolites, ng/g dry weight.
# Values reported as "< x" (below quantification) are entered as 0, which makes
# the parent fraction below a lower bound on transformation rather than an
# estimate of it. Spinach and arugula are omitted: their metabolite rows carry
# NAs and censored values in the soils where the roots were not analysed.
#            plant, soil, treatment : [CAR, EPX, OXC, RTC, DHC]
KOD_ROOT = {
    ("lamb's lettuce", "HCh", "M"): [1900, 170, 17, 0, 3.5],
    ("lamb's lettuce", "HCh", "S"): [3000, 320, 36, 0, 4],
    ("lamb's lettuce", "HCa", "M"): [2000, 290, 30, 0, 3.7],
    ("lamb's lettuce", "HCa", "S"): [2400, 440, 46, 0, 4.7],
    ("lamb's lettuce", "AE", "M"): [4400, 390, 40, 19, 8.2],
    ("lamb's lettuce", "AE", "S"): [8600, 1300, 160, 98, 15],
    ("radish", "HCh", "M"): [2100, 80.5, 0.66, 0, 3],
    ("radish", "HCh", "S"): [2600, 67, 1.1, 0, 4.3],
    ("radish", "HCa", "M"): [2100, 84.5, 2.2, 0, 2.55],
    ("radish", "HCa", "S"): [3100, 110, 1.5, 0, 4.2],
    ("radish", "AE", "M"): [5700, 195, 4.2, 0, 9.7],
    ("radish", "AE", "S"): [7900, 340, 9.7, 0, 12],
}
KOD_LEAF = {
    ("lamb's lettuce", "HCh", "M"): [6400, 17000, 2400, 190, 10],
    ("lamb's lettuce", "HCh", "S"): [4600, 20000, 3300, 220, 8.2],
    ("lamb's lettuce", "HCa", "M"): [3300, 17000, 2500, 300, 5.6],
    ("lamb's lettuce", "HCa", "S"): [1400, 10000, 1500, 220, 0],
    ("lamb's lettuce", "AE", "M"): [1600, 5700, 710, 100, 0],
    ("lamb's lettuce", "AE", "S"): [2700, 12000, 1800, 580, 0],
    ("radish", "HCh", "M"): [19000, 4100, 200, 100, 28],
    ("radish", "HCh", "S"): [21000, 5600, 270, 150, 32],
    ("radish", "HCa", "M"): [11000, 1500, 66, 44, 16.5],
    ("radish", "HCa", "S"): [15000, 2000, 100, 68, 22],
    ("radish", "AE", "M"): [31500, 8800, 450, 285, 48],
    ("radish", "AE", "S"): [53000, 15000, 1000, 620, 89],
}
# every (plant, soil, treatment) cell of Table S2 that reports both organs
KOD_LEAF_ROOT = [
    (1900, 6400), (3000, 4600), (2000, 3300), (2400, 1400), (4400, 1600), (8600, 2700),
    (990, 2900), (1500, 2000), (2800, 2300), (2800, 1900), (5200, 4600),
    (1600, 6900), (1800, 8400), (4000, 13000), (3400, 23000),
    (2100, 19000), (2600, 21000), (2100, 11000), (3100, 15000),
    (5700, 31500), (7900, 53000),
]
HARVEST_D = 23.0        # 20 d (lamb's lettuce, radish) to 26 d (spinach, arugula)


def section6_translocation(drv):
    """The leaf data, which the root-only comparison above leaves on the table.

    A leaf/root ratio needs no exposure at all -- it cancels -- so it is the one
    thing in this dataset that tests TRANSLOCATION rather than partition, and it
    does so without inheriting the isotherm assumption section 1 defends.
    """
    print("\n" + "=" * 84)
    print("6. TRANSLOCATION — leaf/root, with the exposure cancelled")
    print("=" * 84)
    obs = np.array([l / r for r, l in KOD_LEAF_ROOT])
    print(f"   measured leaf/root, n={len(obs)}: median {np.median(obs):.2f}, "
          f"range {obs.min():.2f}-{obs.max():.2f}")
    m = ND.simulate_neutral(ND.NeutralCompound("carbamazepine", LOG_KOW), drv)
    t = np.asarray(drv["t"])
    i = int(np.argmin(np.abs(t - HARVEST_D)))
    lr = m["conc"]["leaf"][i] / m["conc"]["root"][i]
    print(f"   model leaf/root at the matched harvest ({HARVEST_D:.0f} d): {lr:.1f}"
          f"   -> {lr / np.median(obs):.0f}x over")
    print("\n   That is the terminal-leaf runaway of section 3 of the validation doc,")
    print("   measured for the first time against per-organ data with the exposure")
    print("   divided out. But read the next block before calling it a model error.")
    print(f"\n   {'in-planta half-life':>22}{'model leaf/root':>18}")
    for hl in (None, 30, 14, 7, 3, 1.5):
        g = 0.0 if hl is None else float(np.log(2.0) / hl)
        comps = ND.rice_compartments(gammas={k: g for k in TISSUES})
        mm = ND.simulate_neutral(ND.NeutralCompound("carbamazepine", LOG_KOW), drv,
                                 comps=comps)
        lab = "none (recalcitrant)" if hl is None else f"{hl} d"
        last = mm['conc']['leaf'][i] / mm['conc']['root'][i]
        print(f"   {lab:>22}{last:>18.1f}")
    print(f"\n   Metabolism CANNOT close it: even a 1.5-day half-life leaves the model")
    print(f"   {last / np.median(obs):.0f}x above the measured median of {np.median(obs):.2f}. So the leaf excess is")
    print("   NOT only the missing gamma -- most of it is the driver mismatch, a rice")
    print("   season's transpiration per unit leaf mass applied to a 340 cm3 pot of")
    print("   lettuce. That is a WARNING about the Ge 2017 result: its half-life")
    print("   minimum at ~7 d may be absorbing the same mismatch rather than")
    print("   measuring metabolism.")


def section7_measured_metabolism():
    """gamma stops being a free parameter, for one compound at least.

    Every half-life statement in this repo so far has been an inference from a
    fit. Kodesova measured the metabolites alongside the parent, so the parent
    fraction is a direct observation of in-planta transformation.
    """
    print("\n" + "=" * 84)
    print("7. IN-PLANTA TRANSFORMATION, MEASURED RATHER THAN FITTED")
    print("=" * 84)

    def parent_fraction(d):
        return {k: v[0] / sum(v) for k, v in d.items()}

    pr, pl = parent_fraction(KOD_ROOT), parent_fraction(KOD_LEAF)
    print("   parent fraction = CAR / (CAR + its four measured metabolites)")
    print("   ('< LOQ' entered as 0, so these are UPPER bounds on the parent share")
    print("   and therefore LOWER bounds on transformation)\n")
    print(f"   {'plant':18}{'root':>16}{'leaf':>16}")
    for plant in ("lamb's lettuce", "radish"):
        r = [v for k, v in pr.items() if k[0] == plant]
        l = [v for k, v in pl.items() if k[0] == plant]
        print(f"   {plant:18}{np.mean(r):>9.3f} ±{np.std(r):.3f}"
              f"{np.mean(l):>9.3f} ±{np.std(l):.3f}")
    ar, al = list(pr.values()), list(pl.values())
    print(f"   {'both':18}{np.mean(ar):>16.3f}{np.mean(al):>16.3f}")
    print("\n   Two things follow, and the second is the more useful.")
    print("   (i) Transformation is a SHOOT process here: the root stays ~92 % parent")
    print("       while the leaf drops to 49 % on average. So gamma = 0 is simply")
    print("       wrong for carbamazepine, and 'recalcitrant' cannot be assumed for")
    print("       a neutral organic the way it is defensible for PFAS.")
    print("   (ii) It is strongly SPECIES-dependent -- lamb's lettuce leaf 0.17,")
    print("       radish leaf 0.81, a 4.8x difference on the same compound, the same")
    print("       soils and the same harvest. An in-planta half-life is therefore")
    print("       not a compound property, which is what fitting one to a single")
    print("       dataset implicitly assumes.")
    print("\n   Note the metabolites are not inert: in lamb's lettuce leaves the")
    print("   epoxide EXCEEDS the parent (17000 vs 6400 ng/g), so a model that")
    print("   tracks only the parent understates the total burden several-fold.")


def main():
    drv = V.drivers()
    print("KODESOVA 2019 — CARBAMAZEPINE ROOT PARTITION FROM THREE SOILS")
    print("acquisition queue row A4, closed by the SI supplied this session\n")
    section1_isotherm()
    rmse = section2_apriori(drv)
    section3_anchor(drv)
    section4_brunetti()
    section5_sensitivity()
    section6_translocation(drv)
    section7_measured_metabolism()
    print("\n" + "=" * 84)
    print("VERDICT")
    print("=" * 84)
    eq = V.compare_to_obs(OBS, drv, quiet=True, mode="equilibrium")
    print(f"   a-priori log10 RMSE {rmse:.3f} on 21 rows through the ODE, {eq:.3f} on the")
    print("   equilibrium basis, nothing fitted, each plant on its own sourced root.")
    print("\n   What this table does settle:")
    print("     * the Brunetti sighting weakens: a direct measurement of the compound")
    print("       Brunetti calibrated 13.3 for (section 4);")
    print("     * section 7 measures in-planta transformation instead of fitting it:")
    print("       carbamazepine is NOT recalcitrant in the shoot even though its soil")
    print("       DT50 exceeds 1000 d. Soil persistence does not imply plant persistence.")
    print("\n   And what it does NOT settle: four leafy/root vegetables, one compound,")
    print("   one lipophilicity, no rice; the leaf side (section 6) runs on rice drivers.")
    return rmse


if __name__ == "__main__":
    main()
