#!/usr/bin/env python3
# =============================================================================
# validation/oryza2000_kim_face_gao_japonica.py
# -----------------------------------------------------------------------------
# A JAPONICA CROP SET IN THE ORYZA2000 PORT: Gao et al. 2021 (cold-region japonica,
# Suijing-18, Heilongjiang) transferred to the Shizukuishi FACE data (Akitakomachi)
#
# PRE-REGISTRATION -- written and pushed before the port was run with the Gao crop set on
# any weather. Already seen, and the expectation below is written knowing it: every result
# of oryza2000_kim_face.py (5b), oryza2000_kim_face_phenology.py (option C) and
# oryza2000_kim_face_nursery_events.py (the redo): with the recorded nursery the IR72.D92
# port PASSES option C's rule C1 against the shipped driver (N: pooled 0.220 vs 0.244), but
# with the 1999 event dates it FAILS on (b) (NE: dead-leaf share Q2x 0.281 vs 0.181),
# because it over-predicts the senescent share in every year (+0.10...+0.37 log) -- the
# excess option C's C2 put in IR72's leaf-death table DRLVT, which Kim's data cannot
# re-fit (no stem/panicle weights on the same dates).
#
# Why: the port's remaining weakness on independent data is a CROP-FILE question (IR72.D92
# is a tropical indica), not a code question. Gao et al. 2021 (Agric. Water Manag. 250:106864;
# docs/literature_db/raw_si/gao2021_cold_japonica_oryza_params.csv) publish an ORYZA2000
# crop set calibrated on a temperate (cold-region) japonica: development rates, leaf/stem/
# panicle partitioning, a leaf-death table, an SLA table, RGRLMX, FRPAR, a single KDF and
# COLDMIN. Its calibration data (Qing'an 2018) cannot be re-run here -- the daily weather
# was not published -- so this is a TRANSFER test: does a japonica crop set, calibrated
# on another japonica at another site, fix the port's organ shapes on Kim's data?
#
# Arms (all three with phenology calibrated to Akitakomachi, DRATES-style, exactly as in
# the redo's N design; the calibration overwrites each arm's DVRJ..DVRR, so the arms differ
# ONLY in their non-phenology parameters):
#   A  the port with IR72.D92 (= the redo's arm A).
#   G  the port with Gao's calibrated values for FLVTB, FSTTB, FSOTB, DRLVT, SLATB (switch
#      SWISLA = TABLE), KDF (as a constant KDFTB), RGRLMX, FRPAR, and Gao's DVRJ..DVRR as the
#      base the calibration starts from. Everything Gao does not list stays IR72.D92 (as in
#      Gao's own crop file). NOT applied: COLDMIN (the port has no cold-stop rule; Gao lowers
#      it from 12 to 8 oC, i.e. it binds less, and Shizukuishi's season means stay well above
#      8 oC before maturity). Gao's management values (5 plants/hill, 20 hills/m2) are NOT
#      used -- Kim's management stays (3 plants/hill, 19.05 hills/m2).
#   B  the shipped oryza_growth (= the redo's arm B).
#
# Designs (LOYO, station weather, as the redo):
#   N    PRIMARY. Recorded nursery (SBDUR 14/23/24 d, TMPSB 9 oC), sampling-day targets,
#        option C's stage placement.
#   NE   SECONDARY. N + the 1999 event dates (PI 54, heading 78), sample-DVS placement --
#        the design on which arm A failed.
#   Sensitivities of N (rule C1 for G vs B only): (s1) TMPSB 0; (s2) nursery from sowing.
#   DESCRIPTIVE (reported, no rule):
#     * G with Gao's OWN development rates (no calibration; IR72-like 5b scoring on its own
#       DVS) -- does Gao's phenology transfer? Phenology error in days per year and stage.
#     * One group at a time swapped into IR72.D92 (each with calibrated phenology), design N:
#         A+DRLVT   Gao's leaf-death table only;
#         A+PART    Gao's FLVTB/FSTTB/FSOTB only;
#         A+CANOPY  Gao's SLATB (TABLE), KDF, RGRLMX, FRPAR only.
#       Says which part of the crop set moves which score.
#
# Scoring: option C's stage-matched ratios (Q1 green-leaf course, Q2x senescent share
#   without the anthesis point, Q3 leaf share, Q4 dry-matter timing, Q5 root share),
#   log10 RMSE pooled over years, POOLEDx = mean of the five.
#
# Decision rules:
#   G1 (primary, N) -- the Gao leaf-death table FIXES the dead-leaf excess if
#      Q2x(G) <= Q2x(A) - 0.05; otherwise recorded as NOT fixed by it.
#   G2 -- option C's rule C1 with G in place of A, on N and on NE:
#        (a) G POOLEDx <= B POOLEDx + 0.05; (b) G not worse than B by > 0.05 on Q1 or Q2x;
#        (c) G reaches DVS 2 in all three held-out years.
#      Consequence recorded:
#        passes N AND NE  -> a japonica crop set makes the port's pass robust to the 1999
#                            events: a default switch would need a japonica crop file, not
#                            IR72.D92 (Gao's set is a CANDIDATE, not adopted -- G5);
#        passes N, fails NE -> the same fragility as IR72.D92;
#        fails N           -> the Gao set does not transfer to Akitakomachi.
#   G3 -- no general loss: G POOLEDx <= A POOLEDx + 0.02 on N (reported with G1/G2).
#   G4 -- robustness: G2 on N re-evaluated under (s1), (s2); a flip is reported.
#   G5 -- no default, no parameter file and the port's default crop file (IR72.D92) do not
#      change. The Gao values live in this file's crop override only.
#
# Limits stated in advance: Gao calibrated with ORYZA2000's water- and N-limited modules
#   (PADDY, N balance), so their crop values may compensate for stresses this potential-
#   production port does not simulate; their calibration matched GREEN leaf weight, not dead
#   leaves, so their DRLVT is tuned to green-leaf decline, not to the senescent share scored
#   here; Suijing-18 at 47 N is a different (earlier) cultivar than Akitakomachi at 39.6 N.
#
# EXPECTATION (stated so it can be wrong; from the parameter tables only, NOT from any run
#   with the Gao set):
#   * Gao's leaf death starts at DVS 0.4 (0.010 d-1) and runs ABOVE IR72's up to DVS ~1.3
#     (0.024 vs 0.015 at DVS 1), BELOW it after DVS 1.6 (0.020 vs 0.025-0.050). Dead leaves
#     accumulate in WLVD, so the early extra death offsets the late saving at maturity and
#     raises the share at mid-ripening: G1 NOT fixed (Q2x falls by < 0.05, may rise).
#     A+DRLVT alone behaves like G on Q2x.
#   * Gao's lower leaf fraction at DVS 0.5 (0.50 vs 0.60) and the early leaf death lower the
#     leaf share: the port already under-predicts it (Q3 bias -0.076), so Q3 gets worse.
#   * G3 fails (G pooled > A pooled + 0.02); G2 fails on NE (as A) and its N verdict is not
#     better than A's.
#   * Gao's OWN development rates sit inside the range option C calibrated for Akitakomachi
#     (DVRJ x1.15, DVRP x1.35, DVRR x1.12 of IR72's), so uncalibrated they come within 7 d of
#     the observed anthesis and 10 d of maturity in every year (descriptive).
#
#   python validation/oryza2000_kim_face_gao_japonica.py
#
# AMENDMENT 2 (shared scorer): the observation operator and the decision rule were corrected --
# see the AMENDMENT 2 block in validation/oryza2000_kim_face.py (senescent share scored as a
# floor, Q3-Q5 on living biomass, paired-bootstrap decision instead of fixed margins).
# =============================================================================
from __future__ import annotations

import csv
import os
import sys

import numpy as np
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, HERE)

import oryza2000 as O2K                         # noqa: E402
import oryza2000_kim_face as K                  # noqa: E402
import oryza2000_kim_face_phenology as C        # noqa: E402
import oryza2000_kim_face_nursery_events as E   # noqa: E402

YEARS = K.YEARS
GAO_CSV = os.path.join(ROOT, "docs", "literature_db", "raw_si", "gao2021_cold_japonica_oryza_params.csv")
RATE_KEYS = C.RATE_KEYS
GROUPS_GAO = {
    "DRLVT": ("DRLVT",),
    "PART": ("FLVTB", "FSTTB", "FSOTB"),
    "CANOPY": ("SLATB", "SWISLA", "KDFTB", "RGRLMX", "FRPAR"),
    "DVR": RATE_KEYS,
}
SCORE_GROUPS = C.SCORE_GROUPS


def _table(v):
    pairs = [p.split(":") for p in v.split(";")]
    return np.array([float(a) for a, _ in pairs]), np.array([float(b) for _, b in pairs])


def gao_values():
    """Gao 2021 Table 3 calibrated values, in the port's crop-dict form (COLDMIN not applied)."""
    with open(GAO_CSV) as f:
        rows = {r["parameter"]: r["calibrated"] for r in csv.DictReader(l for l in f if not l.startswith("#"))}
    out = {k: _table(rows[k]) for k in ("FLVTB", "FSTTB", "FSOTB", "DRLVT", "SLATB")}
    out["SWISLA"] = "TABLE"
    out["KDFTB"] = (np.array([0.0, 2.5]), np.array([float(rows["KDF"])] * 2))    # Table A1: a single KDF
    for k in ("RGRLMX", "FRPAR") + RATE_KEYS:
        out[k] = float(rows[k])
    return out


GAO = gao_values()
ARM_CROP = {                                     # port arms: overrides on IR72.D92
    "A": {},
    "G": {k: GAO[k] for g in GROUPS_GAO.values() for k in g},
    "A+DRLVT": {k: GAO[k] for k in GROUPS_GAO["DRLVT"]},
    "A+PART": {k: GAO[k] for k in GROUPS_GAO["PART"]},
    "A+CANOPY": {k: GAO[k] for k in GROUPS_GAO["CANOPY"]},
}
MAIN_ARMS = ("A", "G", "B")
DECOMP_ARMS = ("A+DRLVT", "A+PART", "A+CANOPY")


def crop(arm):
    c = O2K.load_crop()
    c.update(ARM_CROP[arm])
    return c


def base_rates(arm):
    if arm == "B":
        return C.base_rates("B")
    c = crop(arm)
    return {k: float(c[k]) for k in RATE_KEYS}


def run(arm, year, rates, tr):
    if arm == "B":
        return K.run_shipped(year, tr["dT"], source=tr["source"], params_kw=dict(rates))
    c = crop(arm)
    c.update(rates)
    return K.run_port(year, tr["dT"], tr["sbdur"][year], tr["source"], crop=c, tmpsb=tr["tmpsb"])


def calibrate(arm, year, tr, tgt):
    """The redo's DRATES-style calibration (common J/I factor -> PI; DVRP; DVRR)."""
    base = base_rates(arm)
    rates = dict(base)

    def miss(lim, day):
        def g(r_):
            x = C.crossing(run(arm, year, r_, tr), lim)
            return (x if x is not None else 1.0e3) - day
        return g

    def solve(fun, lo, hi):
        for _ in range(6):
            if fun(lo) * fun(hi) < 0:
                return brentq(fun, lo, hi, xtol=1e-7, rtol=1e-7)
            lo, hi = lo / 2.0, hi * 2.0
        raise RuntimeError(f"no bracket ({arm}, {year})")

    g1 = miss(0.65, tgt[year]["panicle_initiation"])
    f = solve(lambda x: g1({**rates, "DVRJ": base["DVRJ"] * x, "DVRI": base["DVRI"] * x}), 0.3, 3.0)
    rates.update(DVRJ=base["DVRJ"] * f, DVRI=base["DVRI"] * f)
    g2 = miss(1.0, tgt[year]["anthesis"])
    rates["DVRP"] = solve(lambda x: g2({**rates, "DVRP": x}), base["DVRP"] / 3, base["DVRP"] * 3)
    g3 = miss(2.0, tgt[year]["maturity"])
    rates["DVRR"] = solve(lambda x: g3({**rates, "DVRR": x}), base["DVRR"] / 3, base["DVRR"] * 3)
    return rates, f


def loyo_rates(own, arm):
    base = base_rates(arm)
    out = {}
    for y in YEARS:
        others = [own[o] for o in YEARS if o != y]
        fbar = float(np.mean([o[1] for o in others]))
        out[y] = {"DVRJ": base["DVRJ"] * fbar, "DVRI": base["DVRI"] * fbar,
                  "DVRP": float(np.mean([o[0]["DVRP"] for o in others])),
                  "DVRR": float(np.mean([o[0]["DVRR"] for o in others]))}
    return out


def run_design(tr, obs, root_obs, arms):
    tgt = E.targets(tr, obs)
    out = dict(targets=tgt)
    for arm in arms:
        own = {y: calibrate(arm, y, tr, tgt) for y in YEARS}
        lr = loyo_rates(own, arm)
        runs = {y: run(arm, y, lr[y], tr) for y in YEARS}
        own_runs = {y: run(arm, y, own[y][0], tr) for y in YEARS}
        rl, bl, dl, sd = E.score(runs, obs, root_obs, tr, own_runs)
        ro, _, _, _ = E.score(own_runs, obs, root_obs, tr, own_runs)
        out[arm] = dict(own=own, loyo_rates=lr, runs=runs, own_runs=own_runs, rmse=rl, bias=bl,
                        detail=dl, sample_dvs=sd, own_rmse=ro, phen=E.phenology_errors(runs, tgt))
    return out


def c1(R, arm, ref="B"):
    """AMENDMENT 2 rule: paired bootstrap of POOLEDx(ref) - POOLEDx(arm) + arm reaches DVS 2."""
    return K.switch_rule(R[arm]["detail"], R[ref]["detail"])


def uncalibrated_gao(obs, root_obs, tr):
    """G with Gao's own development rates; option C's stage placement on its own DVS."""
    runs = {y: run("G", y, base_rates("G"), tr) for y in YEARS}
    rmse, bias, detail = K.evaluate(runs, obs, root_obs)
    tgt = E.targets(tr, obs)
    return dict(runs=runs, rmse=rmse, bias=bias, detail=detail, phen=E.phenology_errors(runs, tgt))


def main(verbose=True):
    obs, root_obs = K.load_obs()
    D = E.DESIGNS
    res = {"N": run_design(D["N"], obs, root_obs, MAIN_ARMS + DECOMP_ARMS),
           "NE": run_design(D["NE"], obs, root_obs, MAIN_ARMS)}
    for s in ("N_s1_tmpsb0", "N_s2_sowing"):
        res[s] = run_design(D[s], obs, root_obs, ("G", "B"))
    redo = E.main(verbose=False)["res"]
    checks = {d: all(abs(res[d][a]["rmse"][g] - redo[d][a]["rmse"][g]) < 1e-9
                     for a in ("A", "B") for g in SCORE_GROUPS + ("POOLEDx",)) for d in ("N", "NE")}
    unc = uncalibrated_gao(obs, root_obs, D["N"])
    N, NE = res["N"], res["NE"]
    g2 = {d: c1(res[d], "G") for d in ("N", "NE")}                 # G vs shipped
    g3 = {d: c1(res[d], "G", ref="A") for d in ("N", "NE")}        # G vs IR72 port
    g4 = {s: c1(res[s], "G") for s in ("N_s1_tmpsb0", "N_s2_sowing")}
    flips = [s for s, v in g4.items() if v["verdict"] != g2["N"]["verdict"]]
    out = dict(res=res, checks=checks, unc=unc, g2=g2, g3=g3, g4=g4, flips=flips)
    # POST-HOC: the decisions under the other Q2 reading (see K's Q2 READING note)
    out["attached"] = {(d, lab): K.switch_rule(res[d][arm]["detail"], res[d][ref]["detail"],
                                               reading="attached")
                       for d in ("N", "NE")
                       for lab, arm, ref in (("G vs shipped", "G", "B"), ("G vs IR72 port", "G", "A"),
                                             ("IR72 port vs shipped", "A", "B"))}
    if not verbose:
        return out

    print("=" * 78)
    print("GAO 2021 JAPONICA CROP SET IN THE PORT -- Shizukuishi FACE, LOYO, station weather")
    print("=" * 78)
    print(f"\n0. Checks: arms A and B reproduce the redo exactly: N {checks['N']}, NE {checks['NE']}")

    print("\n1. Gao's OWN development rates (no calibration), design N -- phenology error (days) and scores")
    for y in YEARS:
        e = unc["phen"][y]
        print(f"    {y}: " + "  ".join(f"{s.split('_')[0]} {'n/r' if v is None else f'{v:+.1f}'}" for s, v in e.items()))
    r = unc["rmse"]
    print("    log10 RMSE " + "  ".join(f"{g} {r[g]:.3f}" for g in SCORE_GROUPS) + f"  POOLEDx {r['POOLEDx']:.3f}")

    print("\n2. Calibrated rates per year (own-year, design N), relative to each arm's base")
    for arm in MAIN_ARMS:
        base = base_rates(arm)
        print(f"  {arm:<3} " + " | ".join(
            f"{y} J/I x{N[arm]['own'][y][1]:.3f} P x{N[arm]['own'][y][0]['DVRP'] / base['DVRP']:.3f}"
            f" R x{N[arm]['own'][y][0]['DVRR'] / base['DVRR']:.3f}" for y in YEARS))

    print("\n3. Stage-matched ratios, design N (observed | A | G | B)")
    for g in K.GROUPS:
        print(f"  {g}" + ("   [anthesis point reported, not scored]" if g == "Q2" else ""))
        for y in YEARS:
            qo = N["A"]["detail"][y]["qo"][g]
            for st in qo:
                print(f"    {y} {st:<19} {qo[st]:7.3f} | " + " | ".join(
                    f"{N[a]['detail'][y]['qm'][g][st]:7.3f}" for a in MAIN_ARMS))

    print("\n4. log10 RMSE, LOYO (bias for Q1/Q3/Q4/Q5; AMENDMENT 2 scoring)")
    print("                " + "".join(f"{g:>16}" for g in SCORE_GROUPS) + f"{'POOLED':>9}  own-yr")
    for d in res:
        for arm in [a for a in MAIN_ARMS + DECOMP_ARMS if a in res[d]]:
            rr, b = res[d][arm]["rmse"], res[d][arm]["bias"]
            cells = [f"  {rr[g]:6.3f}" + ("         " if g == "Q2x" else f" ({b[g]:+.3f})") for g in SCORE_GROUPS]
            print(f"  {d:<11} {arm:<9}" + "".join(cells) + f"  {rr['POOLEDx']:6.3f}  {res[d][arm]['own_rmse']['POOLEDx']:6.3f}")

    print("\n5. Decisions (AMENDMENT 2: paired bootstrap, 95 % interval; G1 void -- its premise was a")
    print("   dead-leaf ERROR, removed by the floor)")
    for d in ("N", "NE"):
        for lab, v in (("G vs shipped", g2[d]), ("G vs IR72 port", g3[d])):
            print(f"  {d:<3} {lab:<15} B-A {v['delta']:+.3f} [{v['lo']:+.3f}, {v['hi']:+.3f}] -> "
                  f"{v['verdict'].replace('A better', 'G better')}; G matures every year: {v['reached']}")
    for s_, v in g4.items():
        print(f"  {s_:<12} G vs shipped B-A {v['delta']:+.3f} [{v['lo']:+.3f}, {v['hi']:+.3f}] -> "
              f"{v['verdict'].replace('A better', 'G better')}")
    print("  " + (f"verdict differs under {', '.join(flips)}" if flips else "same verdict under (s1), (s2)"))
    print("  POST-HOC, Q2 read as 'attached' (no blade had detached; Q2 scored symmetrically):")
    for (d, lab), v in out["attached"].items():
        first, second = lab.split(" vs ")
        verdict = v["verdict"].replace("A better", first + " better").replace("B better", second + " better")
        print(f"  {d:<3} {lab:<21} B-A {v['delta']:+.3f} [{v['lo']:+.3f}, {v['hi']:+.3f}] -> {verdict}")
    return out


if __name__ == "__main__":
    main()
