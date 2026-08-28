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
# It also lands on the other side of the decision left open by section 4d: the Li
# 2019 table wants the Briggs anchor restored, the Liu 2023 rice table does not,
# and the tie-break was missing. This is a third vote, and a well-conditioned one:
# carbamazepine sits at log Kow 2.25, where the shipped and anchored compositions
# differ by 1.6x, so the two are cleanly separable.
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
    """The decision this table was fetched to inform."""
    print("\n" + "=" * 84)
    print("3. WHICH WAY DOES IT VOTE ON THE BRIGGS ANCHOR?")
    print("=" * 84)
    orig = ND.TRAPP1994_LIPID_FW["root"]
    out = {}
    try:
        for lab, L in (("shipped  L=0.010", orig),
                       ("anchored L=0.0247", ND.BRIGGS_ANCHORED_LIPID_FW["root"])):
            ND.TRAPP1994_LIPID_FW["root"] = L
            out[lab] = [V.compare_to_obs(p, drv, quiet=True)
                        for p in (OBS, LIU, LI2019)]
    finally:
        ND.TRAPP1994_LIPID_FW["root"] = orig
    print(f"   {'composition':22}{'Kodesova n=21':>15}{'Liu rice n=14':>15}{'Li2019 n=29':>14}")
    for lab, v in out.items():
        print(f"   {lab:22}{v[0]:>15.3f}{v[1]:>15.3f}{v[2]:>14.3f}")
    print("\n   Carbamazepine at log Kow 2.25 separates the two compositions by 1.6x,")
    print("   so this is a well-conditioned vote rather than a marginal one -- and it")
    print("   votes with Liu (the rice table) and against Li 2019.")
    print("   That makes the tally 2 tables against restoring the anchor, 1 for it,")
    print("   and the two against are the ones measured on soil-grown plants at")
    print("   moderate lipophilicity -- the regime this model is actually used in.")

    # driver-free cross-check: the ODE root BAF is not K_PW, so show both framings
    obs = np.array([float(r["value"]) for r in rows()]) * (1.0 - ND.RICE_WATER["root"])
    kpw_ship = ND.k_pw(LOG_KOW, W=ND.RICE_WATER["root"], L=ND.TRAPP1994_LIPID_FW["root"])
    kpw_anch = ND.k_pw(LOG_KOW, W=ND.RICE_WATER["root"],
                       L=ND.BRIGGS_ANCHORED_LIPID_FW["root"])
    print("\n   Driver-free cross-check (equilibrium K_PW, no ODE, so no rice drivers):")
    print(f"      measured RCF_fw   median {np.median(obs):.2f}   range "
          f"{obs.min():.2f}-{obs.max():.2f}")
    for lab, k in (("shipped K_PW", kpw_ship), ("anchored K_PW", kpw_anch)):
        b = float(np.mean(np.log10(k) - np.log10(obs)))
        print(f"      {lab:16} {k:.2f}   mean log10 bias {b:+.3f}")
    print("   The absolute numbers differ from section 2 -- the ODE's root BAF sits")
    print("   below K_PW because the xylem drains it, and log Kow 2.25 is near the")
    print("   TSCF peak where that drain is largest -- but the ORDERING is the same")
    print("   under both framings, which is what the decision rests on.")


def section4_brunetti():
    print("\n" + "=" * 84)
    print("4. THE BRUNETTI QUESTION, WHICH IS WHAT THIS PAPER WAS FETCHED FOR")
    print("=" * 84)
    obs = np.array([float(r["value"]) for r in rows()]) * (1.0 - ND.RICE_WATER["root"])
    print(f"   Brunetti 2021, calibrated green-pea root K_RW   {BRUNETTI_K_RW:.1f} L/kg fw")
    print(f"   Kodesova 2019, MEASURED, 4 plants x 3 soils     {np.median(obs):.2f} "
          f"(range {obs.min():.2f}-{obs.max():.2f})")
    print(f"   Briggs K_PW as this model runs it               "
          f"{ND.k_pw(LOG_KOW, W=ND.RICE_WATER['root'], L=ND.TRAPP1994_LIPID_FW['root']):.2f}")
    print(f"\n   Brunetti's value is ~{BRUNETTI_K_RW / np.median(obs):.0f}x ABOVE a direct measurement of the")
    print("   same compound, while Briggs lands within a factor of ~1.5 of it. So the")
    print("   disagreement is far more likely to be in Brunetti's calibration -- a")
    print("   posterior that can absorb whatever else the model was missing -- than in")
    print("   the partition core. Sighting 1 of the four in section 5 of the validation")
    print("   doc WEAKENS, and this table is a COUNTER-EXAMPLE to the group.")
    print("\n   Note the species caveat both ways: Brunetti is green pea, this is four")
    print("   leafy/root vegetables, and neither is rice.")


def section5_sensitivity():
    print("\n" + "=" * 84)
    print("5. WHAT THE CONCLUSION RESTS ON")
    print("=" * 84)
    raw = np.array([float(r["value"]) for r in rows()])          # L/kg DRY
    kpw_ship = ND.k_pw(LOG_KOW, W=ND.RICE_WATER["root"], L=ND.TRAPP1994_LIPID_FW["root"])
    kpw_anch = ND.k_pw(LOG_KOW, W=ND.RICE_WATER["root"],
                       L=ND.BRIGGS_ANCHORED_LIPID_FW["root"])
    print("   The dw->fw conversion is a real lever, so it is scanned rather than")
    print("   asserted. Root water content assumed:")
    print(f"   {'theta':>8}{'median RCF_fw':>16}{'shipped bias':>15}{'anchored bias':>16}")
    for th in (0.85, 0.90, 0.93, 0.95):
        v = raw * (1.0 - th)
        print(f"   {th:>8.2f}{np.median(v):>16.2f}"
              f"{float(np.mean(np.log10(kpw_ship) - np.log10(v))):>+15.3f}"
              f"{float(np.mean(np.log10(kpw_anch) - np.log10(v))):>+16.3f}")
    print("\n   Across the whole plausible range the anchored composition is further")
    print("   from the data than the shipped one, by a constant 0.21 log. The lever")
    print("   moves both together, so it cannot flip the comparison -- only the")
    print("   absolute verdict on whether the model is high or low.")
    print("\n   What COULD flip it: the isotherm unit reading in section 1. Reading KF")
    print("   in g/cm3 lowers the derived pore water ~6x, raises every RCF ~6x, and")
    print("   would put the data ABOVE both compositions -- reversing the vote. That")
    print("   reading is rejected on the Koc evidence, but it is the assumption to")
    print("   attack if this result is ever doubted.")


def main():
    drv = V.drivers()
    print("KODESOVA 2019 — CARBAMAZEPINE ROOT PARTITION FROM THREE SOILS")
    print("acquisition queue row A4, closed by the SI supplied this session\n")
    section1_isotherm()
    rmse = section2_apriori(drv)
    section3_anchor(drv)
    section4_brunetti()
    section5_sensitivity()
    print("\n" + "=" * 84)
    print("VERDICT")
    print("=" * 84)
    print(f"   a-priori log10 RMSE {rmse:.3f} on 21 rows, nothing fitted -- the best")
    print("   a-priori result in this repo, on the cleanest exposure it has (measured")
    print("   soil concentration + the paper's own measured isotherm, same pot, same")
    print("   harvest, on a compound that is un-ionised everywhere and effectively")
    print("   non-degrading, DT50 > 1000 d).")
    print("   It does TWO things to the open questions. It weakens the Brunetti")
    print("   sighting, by measuring ~1.1 where Brunetti calibrated 13.3. And it votes")
    print("   AGAINST restoring the Briggs anchor, with Liu and against Li 2019.")
    print("   Neither is decisive on its own -- four vegetables, one compound, one")
    print("   lipophilicity, no rice -- but the direction is now 2:1 on soil-grown")
    print("   plants, which is the regime this model is used in.")
    return rmse


if __name__ == "__main__":
    main()
