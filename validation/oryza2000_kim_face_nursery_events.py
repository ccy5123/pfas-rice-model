#!/usr/bin/env python3
# =============================================================================
# validation/oryza2000_kim_face_nursery_events.py
# -----------------------------------------------------------------------------
# OPTION C REDONE WITH THE RECORDED NURSERY AND THE 1999 EVENT DATES
# (POST-HOC ADDENDUM to validation/oryza2000_kim_face_phenology.py)
#
# PRE-REGISTRATION -- written and pushed before any driver was run with these settings.
# Already seen, and the expectation below is written knowing it: every result of option C
# (oryza2000_kim_face_phenology.py, re-run after the Kim 2003 Fig. 2 re-read): C1 PASSES,
# pooled 0.222 (port) vs 0.244 (shipped); Q1 0.174/0.287, Q2x 0.256/0.221, Q3 0.148/0.133,
# Q4 0.207/0.306, Q5 0.324/0.274; own-year 0.246/0.265; port J/I factor x1.352/1.092/1.156
# (1998/1999/2000); LOYO port phenology error PI +10.6/-6.3/-2.4, anthesis +10.5/-10.3/+0.7,
# maturity +2.0/-2.1/-3.4 d; C2 drop 0.072; C3 NOT adequate; its SBDUR sensitivity moved
# the port's pooled score 0.233 (20 d) / 0.222 (30 d) / 0.214 (40 d), verdict unchanged.
#
# Why: option C used two inputs that sources transcribed AFTER it ran now contradict.
#   N1. NURSERY. Option C ran the port from emergence = transplanting - 30 d on an OPEN
#       seedbed (SBDUR 30 assumed, TMPSB 0). Kim et al. 2003 state the seedlings were
#       transplanted 14, 23 and 24 d after emergence (1998/1999/2000), raised in trays in
#       PLASTIC CHAMBERS; Sasaki et al. 2005 confirm 23 d and the plastic chambers for 1999
#       (docs/literature_db/raw_si/sasaki2005_shizukuishi_1999.csv). ORYZA's input for a
#       covered seedbed is TMPSB, the rise of the daily maximum temperature under the cover
#       while CROPSTA <= 2; the ORYZA2000 manual gives "In Japan, the maximum temperature is
#       raised by 9 oC (Horie)" (its example file: 0 without cover, 9.5 with). The seedbed
#       matters beyond the emergence date: the port's transplanting shock halts development
#       for SHCKD x TSTR heat units (SHCKD 0.4, TSTR = heat sum at transplanting), and the
#       seedbed heat sets the seedling size at transplanting. The shipped driver starts at
#       transplanting (DVS 0), so N1 changes ONLY the port arm.
#   E1. EVENTS (1999 only). Option C calibrated to the SAMPLING days Kim labels PI /
#       anthesis / maturity (1999: DAT 62/88/123). Sasaki 2005 -- same experiment, same
#       year -- gives the events: PI ~54 DAT, heading 78 DAT (ambient plots), and calls the
#       DAT-62 and DAT-88 samplings "pre-heading" and "early grain filling". So in 1999 the
#       targets were 8 and 10 d late, and the two samples were not at DVS 0.65 / 1.0.
#       No event dates exist for 1998 or 2000 in any source in hand; no maturity date exists
#       for any year (Sasaki's last dry-weight sampling is DAT 119, Kim's DAT 123).
#
# What changes (nothing else does -- rates calibrated the DRATES way, LOYO, station weather,
# Q2 scored without the anthesis point, C1 rule, all exactly as option C):
#   N1  port SBDUR = 14 / 23 / 24 d (Kim 2003), TMPSB = 9 oC (ORYZA manual, Japan).
#   E1  1999 targets PI 54, DVS 1.0 at heading 78 (ORYZA's DVS 1 is flowering; rice flowers
#       ~1-3 d after heading -- taken as equal, stated not corrected), maturity stays the
#       sampling DAT 123. 1998/2000 targets stay the sampling DATs (option C). NOT applied:
#       shifting 1998/2000 by the 1999 lag. The samplings were timed by the observed crop,
#       and 1998 was the coolest year: subtracting 1999's 10 d would put 1998 heading at
#       71 DAT, before the warmer 1999's 78 -- not plausible.
#   PLACEMENT under E1: a 1999 sample is no longer at a target DVS, so every sample is
#       placed at the DVS that the driver's OWN-YEAR-calibrated run has on the sampling day
#       ("sample DVS"), and on any other run (LOYO) it is read where that run reaches the
#       sample DVS. For a sample taken on a target day this is exactly DVS 0.65 / 1.0 / 2.0,
#       so for 1998/2000 PI, anthesis and maturity it is option C's placement. It differs
#       from option C for MID-RIPENING (option C: at the observed day-fraction between
#       anthesis and maturity), so a control isolates the placement from the event dates.
#
# Designs (all LOYO, station weather):
#   N    PRIMARY. Option C with N1 only; option C's placement. The shipped arm is untouched
#        by N1 and must reproduce option C's shipped numbers exactly (a self-check).
#   NG   CONTROL. N with the sample-DVS placement, targets unchanged (= option C's).
#   NE   SECONDARY. N + E1 with the sample-DVS placement.
#   Sensitivities of N: (s1) TMPSB 0 (cover ignored); (s2) SBDUR counted from SOWING,
#        20 / 27 / 29 d (Hasegawa et al. 2017 Table S1 sowing dates; treats sowing as
#        emergence), TMPSB 9; (s3) POWER cell + 1.7 oC (option C's S_power+1.7C), N1.
#   NE is also run under (s1).
#   N is primary because the nursery is recorded for all three years while the events are
#   known for one; NE uses better information for 1999 only and mixes event targets with
#   sampling targets, so it is a check that the 1999 lag does not drive the verdict.
#
# Decision rules:
#   D1 (primary) -- C1 of option C, unchanged, on N: (a) port POOLEDx <= shipped + 0.05;
#      (b) port not worse by > 0.05 on Q1 or Q2x; (c) the port reaches DVS 2 in all three
#      held-out years. If the verdict differs from option C's PASS, option C's verdict is
#      recorded as dependent on the assumed nursery.
#   D2 -- the same rule on NE and on NG. If NE and N disagree, the verdict is recorded as
#      sensitive to the 1999 event dates; NG says whether a disagreement is the placement.
#   D3 -- the same rule under each sensitivity; a flip is reported as "sensitive to <s>".
#   D4 -- option C's C2 on N: does the port's Q2x fall by >= 0.10 from its uncalibrated
#      station value (0.328)?
#   D5 (descriptive) -- LOYO phenology error per year and stage under N and NE; under NE
#      the 1999 errors are against the events (PI 54, heading 78). Option C's C3 threshold
#      (|anthesis| <= 5 d, |maturity| <= 7 d every year) is reported for both.
#   D6 -- no default, no parameter file and no recorded option C number changes; option C
#      stays the record of what it pre-registered, this file sits beside it.
#
# EXPECTATION (stated so it can be wrong; from option C's results and a heat-sum estimate
#   of the seedbeds, NOT from any run with these settings):
#   * N: the calibration absorbs N1 mostly into the port's J/I factor. With TMPSB 9 the
#     1999 and 2000 covered nurseries (23/24 d) accumulate MORE seedbed heat than option C's
#     30 open days, the 14-d 1998 nursery about as much, so the 1999/2000 factors fall and
#     1998 stays the fastest year: the year-to-year spread of the factor does NOT shrink,
#     and the LOYO PI errors do not improve. The port's pooled score moves by < 0.02 (option
#     C's own SBDUR sensitivity moved it 0.019); D1 PASSES; the shipped arm is identical.
#   * D4: the Q2x drop stays < 0.10 (the dead-leaf excess is IR72's leaf-death table).
#   * NE: 1999's targets move 8-10 d earlier, so 1999's J/I factor and DVRP rise toward
#     1998's. LOYO anthesis errors shrink in magnitude: 1999 within 3 d of the true heading,
#     1998 below option C's +10.5, 2000 turns negative. C3 stays NOT adequate (the maturity
#     errors come from DVRR, which E1 does not touch). D2: NE PASSES; NG and N agree within
#     0.01 on the pooled scores.
#   * D3: no flip under (s1)-(s3).
#
#   python validation/oryza2000_kim_face_nursery_events.py
#
# AMENDMENT 2 (shared scorer): the observation operator and the decision rule were corrected --
# see the AMENDMENT 2 block in validation/oryza2000_kim_face.py (senescent share scored as a
# floor, Q3-Q5 on living biomass, paired-bootstrap decision instead of fixed margins).
# =============================================================================
from __future__ import annotations

import os
import sys

import numpy as np
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, HERE)

import oryza2000 as O2K                    # noqa: E402
import oryza2000_kim_face as K             # noqa: E402
import oryza2000_kim_face_phenology as C   # noqa: E402

YEARS = K.YEARS
ARMS = C.ARMS
SBDUR_EMERGENCE = {1998: 14, 1999: 23, 2000: 24}      # Kim 2003 text (1999 also Sasaki 2005)
SBDUR_SOWING = {1998: 20, 1999: 27, 2000: 29}         # Hasegawa 2017 Table S1 sowing -> transplanting
TMPSB_COVERED = 9.0                                   # ORYZA2000 manual: Japan, plastic cover (Horie)
EVENTS_1999 = {"panicle_initiation": 54.0, "anthesis": 78.0, "maturity": 123.0}   # Sasaki 2005; maturity = sampling
GROUPS = K.GROUPS
SCORE_GROUPS = C.SCORE_GROUPS


def _tr(sbdur, tmpsb, events=False, placement="optionC", source="station", dT=0.0):
    return dict(sbdur=dict(sbdur), tmpsb=tmpsb, events=events, placement=placement, source=source, dT=dT)


# "C_repro" re-runs option C's primary through this file's code (a check of the code, not a result)
DESIGNS = {
    "C_repro": _tr({y: 30 for y in YEARS}, 0.0),
    "N": _tr(SBDUR_EMERGENCE, TMPSB_COVERED),
    "NG": _tr(SBDUR_EMERGENCE, TMPSB_COVERED, placement="sample_dvs"),
    "NE": _tr(SBDUR_EMERGENCE, TMPSB_COVERED, events=True, placement="sample_dvs"),
    "N_s1_tmpsb0": _tr(SBDUR_EMERGENCE, 0.0),
    "N_s2_sowing": _tr(SBDUR_SOWING, TMPSB_COVERED),
    "N_s3_power+1.7C": _tr(SBDUR_EMERGENCE, TMPSB_COVERED, source="power", dT=1.7),
    "NE_s1_tmpsb0": _tr(SBDUR_EMERGENCE, 0.0, events=True, placement="sample_dvs"),
}
SENSITIVITIES = ("N_s1_tmpsb0", "N_s2_sowing", "N_s3_power+1.7C")


def run(arm, year, rates, tr):
    if arm == "A":
        crop = O2K.load_crop()
        crop.update(rates)
        return K.run_port(year, tr["dT"], tr["sbdur"][year], tr["source"], crop=crop, tmpsb=tr["tmpsb"])
    return K.run_shipped(year, tr["dT"], source=tr["source"], params_kw=dict(rates))


def targets(tr, obs):
    """{year: {stage: target day after transplanting}} -- sampling days, or the 1999 events (E1)."""
    t = {y: {s: float(obs[y][s]["dat"]) for s, _ in C.STAGE_DVS} for y in YEARS}
    if tr["events"]:
        t[1999] = dict(EVENTS_1999)
    return t


def calibrate(arm, year, tr, tgt):
    """option C's DRATES-style calibration, to the targets `tgt[year]`."""
    base = C.base_rates(arm)
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


def sample_dvs(own_runs, obs):
    """The DVS each driver's own-year-calibrated run has on each sampling day."""
    return {y: {s: float(np.interp(obs[y][s]["dat"], own_runs[y]["dat"], own_runs[y]["DVS"]))
                for s in K.SCORED} for y in YEARS}


def placed_index(r, sdvs_y):
    """Fractional day index where run r reaches each sample DVS (last day if never)."""
    idx, reached = {}, {}
    for s, lim in sdvs_y.items():
        x = C.crossing(r, lim)
        reached[s] = x is not None
        idx[s] = float(np.interp(x, r["dat"], np.arange(len(r["dat"])))) if x is not None \
            else float(len(r["dat"]) - 1)
    return idx, reached


def evaluate_sample_dvs(runs, obs, root_obs, sdvs):
    """As K.evaluate, with every sample read at its sample DVS (Q2 scored without anthesis)."""
    err = {g: [] for g in GROUPS}
    detail = {}
    for y in YEARS:
        idx, reached = placed_index(runs[y], sdvs[y])
        mv = K.model_values(runs[y], idx)
        root_m = {s: (mv[s]["ROOT"], mv[s]["TOT"]) for s in root_obs} if y == 1999 else None
        qm = K.quantities(mv, root_m)
        qo = K.quantities({s: obs[y][s] for s in K.SCORED}, root_obs if y == 1999 else None)
        for g in GROUPS:
            for s, o in qo[g].items():
                if g == "Q2" and s == "anthesis":
                    continue
                err[g].append(np.log10(max(qm[g][s], 1e-12)) - np.log10(o))
        detail[y] = dict(idx=idx, reached=reached, qm=qm, qo=qo, mv=mv)
    rmse = {("Q2x" if g == "Q2" else g): float(np.sqrt(np.mean(np.square(e)))) for g, e in err.items()}
    bias = {g: float(np.mean(e)) for g, e in err.items()}
    rmse["POOLEDx"] = float(np.mean([rmse[g] for g in SCORE_GROUPS]))
    return rmse, bias, detail


def score(runs, obs, root_obs, tr, own_runs):
    if tr["placement"] == "optionC":
        rmse, bias, detail = K.evaluate(runs, obs, root_obs)
        return rmse, bias, detail, None
    sd = sample_dvs(own_runs, obs)
    rmse, bias, detail = evaluate_sample_dvs(runs, obs, root_obs, sd)
    return rmse, bias, detail, sd


def phenology_errors(runs, tgt):
    return {y: {s: (None if C.crossing(runs[y], lim) is None else C.crossing(runs[y], lim) - tgt[y][s])
                for s, lim in C.STAGE_DVS} for y in YEARS}


def run_design(tr, obs, root_obs, cache):
    key = (tuple(sorted(tr["sbdur"].items())), tr["tmpsb"], tr["events"], tr["source"], tr["dT"])
    if key not in cache:                       # NG shares N's calibration and runs
        tgt = targets(tr, obs)
        cal = {}
        for arm in ARMS:
            own = {y: calibrate(arm, y, tr, tgt) for y in YEARS}
            lr = C.loyo_rates(own, arm)
            cal[arm] = dict(own=own, loyo_rates=lr,
                            runs={y: run(arm, y, lr[y], tr) for y in YEARS},
                            own_runs={y: run(arm, y, own[y][0], tr) for y in YEARS})
        cache[key] = (tgt, cal)
    tgt, cal = cache[key]
    out = dict(targets=tgt)
    for arm in ARMS:
        c = cal[arm]
        rl, bl, dl, sd = score(c["runs"], obs, root_obs, tr, c["own_runs"])
        ro, _, do, _ = score(c["own_runs"], obs, root_obs, tr, c["own_runs"])
        out[arm] = dict(c, rmse=rl, bias=bl, detail=dl, sample_dvs=sd, own_rmse=ro, own_detail=do,
                        phen=phenology_errors(c["runs"], tgt))
    reached = [out["A"]["detail"][y]["reached"]["maturity"] for y in YEARS]
    out["c1"] = C.c1(out["A"]["rmse"], out["B"]["rmse"], reached)
    return out


def c3(phen):
    return all(abs(phen[y]["anthesis"]) <= 5 and phen[y]["maturity"] is not None
               and abs(phen[y]["maturity"]) <= 7 for y in YEARS)


def main(verbose=True):
    obs, root_obs = K.load_obs()
    cache = {}
    res = {name: run_design(tr, obs, root_obs, cache) for name, tr in DESIGNS.items()}
    ref = C.run_treatment(C.TREATMENTS["primary"], obs, root_obs)          # option C as recorded
    unc = K.run_treatment(obs, root_obs, dT=0.0, sbdur=30, source="station")
    q2x_unc = unc["A"]["rmse"]["Q2x"]
    checks = dict(
        repro=all(abs(res["C_repro"][a]["rmse"][g] - ref[a]["rmse"][g]) < 1e-9
                  for a in ARMS for g in SCORE_GROUPS + ("POOLEDx",)),
        shipped_untouched=all(abs(res["N"]["B"]["rmse"][g] - ref["B"]["rmse"][g]) < 1e-9
                              for g in SCORE_GROUPS + ("POOLEDx",)))
    N = res["N"]
    drop = q2x_unc - N["A"]["rmse"]["Q2x"]
    flips = [n for n in SENSITIVITIES if res[n]["c1"]["passed"] != N["c1"]["passed"]]
    out = dict(res=res, ref=ref, checks=checks, q2x_uncal=q2x_unc, d4_drop=drop,
               d1=N["c1"], d2_NE=res["NE"]["c1"], d2_NG=res["NG"]["c1"], flips=flips,
               c3_N=c3(N["A"]["phen"]), c3_NE=c3(res["NE"]["A"]["phen"]))
    if not verbose:
        return out

    print("=" * 78)
    print("OPTION C REDONE -- recorded nursery (N1) and the 1999 event dates (E1); LOYO, station weather")
    print("=" * 78)
    print(f"\n0. Checks: option C primary reproduced by this file's code: {checks['repro']};"
          f" shipped arm untouched by N1: {checks['shipped_untouched']}")

    print("\n1. Calibrated rates per year (own-year), relative to each driver's base values")
    for name in ("C_repro", "N", "NE"):
        for arm in ARMS:
            base = C.base_rates(arm)
            cells = []
            for y in YEARS:
                r, f = res[name][arm]["own"][y]
                cells.append(f"{y} J/I x{f:.3f} P x{r['DVRP'] / base['DVRP']:.3f} R x{r['DVRR'] / base['DVRR']:.3f}")
            print(f"  {name:<8}{arm}  " + " | ".join(cells))

    print("\n2. Phenology prediction, LOYO (model - target, days; NE 1999 = vs the events PI 54 / heading 78)")
    for name in ("C_repro", "N", "NE"):
        for arm in ARMS:
            ph = res[name][arm]["phen"]
            print(f"  {name:<8}{arm}  " + " | ".join(
                f"{y} " + " ".join(f"{'n/r' if v is None else f'{v:+5.1f}'}" for v in ph[y].values())
                for y in YEARS) + "   (PI anth mat)")
    print(f"  C3 (port; |anthesis| <= 5 d, |maturity| <= 7 d every year): N {'adequate' if out['c3_N'] else 'NOT adequate'},"
          f" NE {'adequate' if out['c3_NE'] else 'NOT adequate'}")

    print("\n3. Sample DVS (own-year runs, sample-DVS placement)")
    for name in ("NG", "NE"):
        for arm in ARMS:
            sd = res[name][arm]["sample_dvs"]
            print(f"  {name:<3}{arm}  " + " | ".join(
                f"{y} " + " ".join(f"{sd[y][s]:.3f}" for s in K.SCORED) for y in YEARS)
                + "   (PI anth mid-rip mat)")

    print("\n4. log10 RMSE, LOYO (bias for Q1/Q3/Q4/Q5)")
    print("                     " + "".join(f"{g:>16}" for g in SCORE_GROUPS) + f"{'POOLED':>9}  own-yr")
    for name in DESIGNS:
        for arm in ARMS:
            r, b = res[name][arm]["rmse"], res[name][arm]["bias"]
            cells = []
            for g in SCORE_GROUPS:
                cells.append(f"  {r[g]:6.3f}" + ("         " if g == "Q2x" else f" ({b[g]:+.3f})"))
            print(f"  {name:<16}{arm} " + "".join(cells) + f"  {r['POOLEDx']:6.3f}  {res[name][arm]['own_rmse']['POOLEDx']:6.3f}")

    print("\n5. Verdicts")
    for name in DESIGNS:
        R = res[name]
        d = R["c1"]
        print(f"  {name:<16} pooled A {R['A']['rmse']['POOLEDx']:.3f} vs B {R['B']['rmse']['POOLEDx']:.3f}"
              f"  Q1 {R['A']['rmse']['Q1']:.3f}/{R['B']['rmse']['Q1']:.3f}"
              f"  Q2x {R['A']['rmse']['Q2x']:.3f}/{R['B']['rmse']['Q2x']:.3f}"
              f" | (a) {d['a']} (b) {d['b']} (c) {d['c']} -> {'PASS' if d['passed'] else 'FAIL'}")
    print(f"  D1 (primary N): {'PASSES' if out['d1']['passed'] else 'FAILS'}"
          f" (option C: {'PASS' if ref['c1']['passed'] else 'FAIL'})")
    print(f"  D2: NE {'PASS' if out['d2_NE']['passed'] else 'FAIL'}, NG {'PASS' if out['d2_NG']['passed'] else 'FAIL'}")
    print("  D3: " + (f"sensitive to {', '.join(flips)}" if flips else "no flip under (s1)-(s3)"))
    print(f"  D4: port Q2x {q2x_unc:.3f} (IR72 phenology) -> {N['A']['rmse']['Q2x']:.3f} under N; drop {drop:+.3f} -> "
          + ("largely PHENOLOGICAL" if drop >= 0.10 else "mostly IR72's leaf-death table (DRLVT)"))
    return out


if __name__ == "__main__":
    main()
