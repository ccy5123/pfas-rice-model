#!/usr/bin/env python3
# =============================================================================
# validation/oryza2000_kim_face_phenology.py
# -----------------------------------------------------------------------------
# OPTION C: CALIBRATE BOTH BIOMASS DRIVERS' DEVELOPMENT RATES TO AKITAKOMACHI, THEN
# RE-RUN THE SHIZUKUISHI FACE COMPARISON (follow-up to validation/oryza2000_kim_face.py)
#
# PRE-REGISTRATION -- written and pushed before any driver was run with calibrated rates.
# (Already seen: every result of oryza2000_kim_face.py -- the pre-registered primary verdict
# FAIL, the post-hoc station run, and the stage-matched ratios of both drivers with IR72
# phenology. The expectation below is written with that knowledge and says so.)
#
# Why: in oryza2000_kim_face.py both drivers ran with IR72 phenology. Two things in that
# verdict may be phenology, not biomass structure:
#   * the port's ripening lasted 49-66 d vs observed 35-50 d; ORYZA's leaf death is a rate
#     per DAY, so a longer ripening piles up dead leaves -- part of the port's dead-leaf
#     excess (Q2x 0.328 vs the shipped driver's 0.172) may be this;
#   * the shipped driver never reached maturity, so its "maturity" values were read on
#     Nov 30 at DVS < 2 -- not a stage at all.
# Option C removes the phenology mismatch the way ORYZA itself prescribes for a new
# cultivar (manual Section 8.2, program DRATES: development rates from observed dates of
# sowing, transplanting, panicle initiation, flowering and maturity), then asks whether the
# port beats the shipped driver on organ shapes once both develop like Akitakomachi.
#
# What is calibrated -- ONLY the four development rates of each driver:
#   DVRJ (juvenile, DVS 0-0.4), DVRI (photoperiod-sensitive, 0.4-0.65), DVRP (panicle
#   formation, 0.65-1), DVRR (grain filling, 1-2). Everything else -- partitioning, leaf
#   death, SLA, assimilation, the transplanting shock, PPSE = 0 -- stays IR72.D92 (port) /
#   the shipped OryzaParams (shipped driver). So any difference left is a difference in
#   what the drivers do structurally.
# How (DRATES-style, sequential root-finding on the drivers' own heat units):
#   1. DVRJ and DVRI by ONE common factor on their IR72 values, so DVS 0.65 is reached on
#      the observed PI day. (The juvenile/photoperiod split is not observable from these
#      data: there is no observation at DVS 0.4.)
#   2. DVRP so DVS 1.0 is reached on the observed anthesis day.
#   3. DVRR so DVS 2.0 is reached on the observed maturity day.
#   "Reached" = the fractional day at which DVS crosses the threshold (linear between the
#   daily states), so the root-finding is smooth. Observed days are the sampling DATs at
#   PI / anthesis / maturity (1998 53/81/131, 1999 62/88/123, 2000 56/78/120): sampling
#   at a stage, not the exact event day -- an error of a few days, stated not corrected.
#   The port runs from emergence = transplanting - SBDUR, so its PI target depends on the
#   assumed nursery: SOWING / NURSERY DURATION IS NOT IN THE REPO'S TRANSCRIPTION OF KIM
#   2001/2003, so SBDUR = 30 d is an ASSUMPTION (Japanese mat seedlings), with 20 and 40 d
#   as sensitivity. The seedbed cover (TMPSB) stays 0 as in IR72DS92. The shipped driver
#   starts at transplanting (DVS 0) as it always does.
#
# Out-of-sample by construction -- LEAVE-ONE-YEAR-OUT (primary): for each held-out year the
#   rates are the mean of the two other years' own calibrated rates (common J/I factor
#   averaged; DVRP, DVRR averaged); the held-out year is then run with them and scored.
#   SECONDARY: each year with its OWN calibrated rates (phenology exact by construction --
#   the cleanest test of the organ shapes, but in-sample for phenology).
#
# Data, weather, scoring -- as oryza2000_kim_face.py, with two changes fixed in advance:
#   * weather PRIMARY = the JMA station temperatures + POWER radiation (that file's S4; the
#     observed weather). SENSITIVITY: SBDUR 20 / 40 d (station); POWER cell + 1.7 C (S1).
#   * Q2 is scored WITHOUT the anthesis point (that file's Q2x): the shipped driver has no
#     leaf death before DVS 1, so its anthesis senescent share is exactly 0 and log-scores
#     to infinity; the anthesis point is reported, not scored. POOLED = mean of Q1, Q2x,
#     Q3, Q4, Q5.
#   Scored quantities stay STAGE-MATCHED ratios on each driver's own DVS (Q1 green-leaf
#   course, Q2x senescent share, Q3 leaf share of the plant, Q4 dry-matter timing, Q5 root
#   share). REPORTED, NOT SCORED: phenology prediction error in days (LOYO); the same ratios
#   read on the OBSERVED sampling days (calendar-matched, which adds the phenology error);
#   the calibrated rates themselves and how far they move from IR72.
#
# Decision rules (PRIMARY = LOYO, station weather, SBDUR 30):
#   C1 (the decision) -- with calibrated phenology the independent evidence SUPPORTS the
#      switch to biomass="oryza2000" if ALL hold:
#        (a) port POOLED <= shipped POOLED + 0.05;
#        (b) the port is not worse than the shipped driver by more than 0.05 on Q1 or Q2x;
#        (c) the port reaches DVS 2 by Nov 30 in all three held-out years.
#      Otherwise C1 FAILS and the default stays (option A), with the reason recorded.
#   C2 (what the dead-leaf excess was) -- if the port's Q2x falls by >= 0.10 log from its
#      uncalibrated station value (0.328), the excess is recorded as largely PHENOLOGICAL;
#      otherwise as sitting in IR72's leaf-death table (DRLVT), which this data cannot
#      re-fit (ORYZA's PARAM needs stem and panicle weights on the same dates; Kim has leaf
#      blades, totals and roots only).
#   C3 (phenology prediction, descriptive) -- LOYO |error| <= 5 d at anthesis and <= 7 d at
#      maturity in every year = adequate. Not a gate.
#   C4 Robustness: C1 is re-evaluated under every sensitivity treatment; a flip is reported
#      as "sensitive to <treatment>".
#   C5 This file changes NO default and no parameter file. The calibrated rates are an
#      Akitakomachi-at-Shizukuishi result, not a temperate-japonica parameter set: the PFAS
#      datasets use other cultivars (Koshihikari/Kasalath, Nipponbare, Jinyu, unknown).
#
# EXPECTATION (stated so it can be wrong; from the uncalibrated run and the parameter
#   tables, NOT from any calibrated run): calibration shortens the port's ripening, so its
#   dead-leaf share falls, but its leaf death starts at DVS 0.6 and its DRLVT is steeper than
#   the shipped one, so it stays above the observed floor; the port's Q2x falls by LESS than
#   0.10 (C2: the excess is mostly DRLVT). The shipped driver, once it matures in 35-50 d,
#   lands near the observed senescent share. The port keeps its Q4 lead; the shipped driver
#   keeps its Q3 lead. C1 FAILS on (b). LOYO phenology is adequate (C3). The verdict does not
#   flip under the sensitivities.
#
#   python validation/oryza2000_kim_face_phenology.py
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

import oryza2000 as O2K           # noqa: E402
import oryza_growth as OG         # noqa: E402
import oryza2000_kim_face as K    # noqa: E402

YEARS = K.YEARS
ARMS = ("A", "B")                 # A = ORYZA2000 port, B = shipped oryza_growth
RATE_KEYS = ("DVRJ", "DVRI", "DVRP", "DVRR")
STAGE_DVS = (("panicle_initiation", 0.65), ("anthesis", 1.0), ("maturity", 2.0))
TREATMENTS = {"primary": dict(source="station", dT=0.0, sbdur=30),
              "S_sbdur20": dict(source="station", dT=0.0, sbdur=20),
              "S_sbdur40": dict(source="station", dT=0.0, sbdur=40),
              "S_power+1.7C": dict(source="power", dT=1.7, sbdur=30)}
SCORE_GROUPS = ("Q1", "Q2x", "Q3", "Q4", "Q5")


def base_rates(arm):
    if arm == "A":
        c = O2K.load_crop()
        return {k: float(c[k]) for k in RATE_KEYS}
    p = OG.OryzaParams()
    return {k: float(getattr(p, k)) for k in RATE_KEYS}


def run(arm, year, rates, tr):
    if arm == "A":
        crop = O2K.load_crop()
        crop.update(rates)
        return K.run_port(year, tr["dT"], tr["sbdur"], tr["source"], crop=crop)
    return K.run_shipped(year, tr["dT"], source=tr["source"], params_kw=dict(rates))


def crossing(r, lim):
    """Fractional day after transplanting at which DVS first reaches `lim` (None if never)."""
    d, t = r["DVS"], r["dat"]
    hit = np.where(d >= lim - 1e-12)[0]
    if not len(hit):
        return None
    i = int(hit[0])
    if i == 0:
        return float(t[0])
    f = (lim - d[i - 1]) / (d[i] - d[i - 1])
    return float(t[i - 1] + f * (t[i] - t[i - 1]))


def calibrate(arm, year, tr, obs):
    """DRATES-style: common factor on DVRJ/DVRI -> PI day; DVRP -> anthesis; DVRR -> maturity."""
    base = base_rates(arm)
    rates = dict(base)
    target = {s: float(obs[year][s]["dat"]) for s, _ in STAGE_DVS}

    def miss(lim, tgt):
        def g(r_):
            x = crossing(run(arm, year, r_, tr), lim)
            return (x if x is not None else 1.0e3) - tgt
        return g

    def solve(fun, lo, hi):
        for _ in range(6):
            if fun(lo) * fun(hi) < 0:
                return brentq(fun, lo, hi, xtol=1e-7, rtol=1e-7)
            lo, hi = lo / 2.0, hi * 2.0
        raise RuntimeError(f"no bracket ({arm}, {year})")

    g1 = miss(0.65, target["panicle_initiation"])
    f = solve(lambda x: g1({**rates, "DVRJ": base["DVRJ"] * x, "DVRI": base["DVRI"] * x}), 0.3, 3.0)
    rates.update(DVRJ=base["DVRJ"] * f, DVRI=base["DVRI"] * f)
    g2 = miss(1.0, target["anthesis"])
    rates["DVRP"] = solve(lambda x: g2({**rates, "DVRP": x}), base["DVRP"] / 3, base["DVRP"] * 3)
    g3 = miss(2.0, target["maturity"])
    rates["DVRR"] = solve(lambda x: g3({**rates, "DVRR": x}), base["DVRR"] / 3, base["DVRR"] * 3)
    return rates, f


def loyo_rates(own, arm):
    """Held-out year -> mean of the other two years' calibrated rates (J/I via the factor)."""
    base = base_rates(arm)
    out = {}
    for y in YEARS:
        others = [own[o] for o in YEARS if o != y]
        fbar = float(np.mean([o[1] for o in others]))
        out[y] = {"DVRJ": base["DVRJ"] * fbar, "DVRI": base["DVRI"] * fbar,
                  "DVRP": float(np.mean([o[0]["DVRP"] for o in others])),
                  "DVRR": float(np.mean([o[0]["DVRR"] for o in others]))}
    return out


def calendar_scores(runs, obs, root_obs):
    """The same ratios read on the OBSERVED sampling days (adds the phenology error)."""
    err = {g: [] for g in K.GROUPS}
    for y in YEARS:
        r = runs[y]
        idx = {s: float(np.interp(obs[y][s]["dat"], r["dat"], np.arange(len(r["dat"]))))
               for s in K.SCORED}
        mv = K.model_values(r, idx)
        qm = K.quantities(mv, {s: (mv[s]["ROOT"], mv[s]["TOT"]) for s in root_obs} if y == 1999 else None)
        qo = K.quantities({s: obs[y][s] for s in K.SCORED}, root_obs if y == 1999 else None)
        for g in K.GROUPS:
            for s, o in qo[g].items():
                if g == "Q2" and s == "anthesis":
                    continue
                err[g].append(K.point_error(g, qm[g][s], o))
    rm = {("Q2x" if g == "Q2" else g): float(np.sqrt(np.mean(np.square(e)))) for g, e in err.items()}
    rm["POOLEDx"] = float(np.mean([rm[g] for g in SCORE_GROUPS]))
    return rm


def phenology_errors(runs, obs):
    return {y: {s: (None if crossing(runs[y], lim) is None
                    else crossing(runs[y], lim) - obs[y][s]["dat"]) for s, lim in STAGE_DVS}
            for y in YEARS}


def run_treatment(tr, obs, root_obs):
    out = {}
    for arm in ARMS:
        own = {y: calibrate(arm, y, tr, obs) for y in YEARS}
        lr = loyo_rates(own, arm)
        runs_l = {y: run(arm, y, lr[y], tr) for y in YEARS}
        runs_o = {y: run(arm, y, own[y][0], tr) for y in YEARS}
        rl, bl, dl = K.evaluate(runs_l, obs, root_obs)
        ro, bo, do = K.evaluate(runs_o, obs, root_obs)
        out[arm] = dict(own=own, loyo_rates=lr, runs=runs_l, rmse=rl, bias=bl, detail=dl,
                        own_rmse=ro, own_detail=do, phen=phenology_errors(runs_l, obs),
                        cal=calendar_scores(runs_l, obs, root_obs))
    out["c1"] = K.switch_rule(out["A"]["detail"], out["B"]["detail"])          # AMENDMENT 2 rule
    out["c1_own"] = K.switch_rule(out["A"]["own_detail"], out["B"]["own_detail"])
    return out


def main(verbose=True):
    obs, root_obs = K.load_obs()
    res = {name: run_treatment(tr, obs, root_obs) for name, tr in TREATMENTS.items()}
    P = res["primary"]
    phen = P["A"]["phen"]
    c3 = all(abs(phen[y]["anthesis"]) <= 5 and phen[y]["maturity"] is not None
             and abs(phen[y]["maturity"]) <= 7 for y in YEARS)
    flips = [n for n, R in res.items() if R["c1"]["verdict"] != P["c1"]["verdict"]]
    out = dict(res=res, obs=obs, root_obs=root_obs, c1=P["c1"], c3=c3, flips=flips)
    # POST-HOC: the decision under the other Q2 reading (see K's Q2 READING note)
    out["c1_attached"] = {n: K.switch_rule(R["A"]["detail"], R["B"]["detail"], reading="attached")
                          for n, R in res.items()}
    if not verbose:
        return out

    print("=" * 78)
    print("OPTION C -- both drivers' development rates calibrated to Akitakomachi (DRATES-style),")
    print("Shizukuishi FACE, leave-one-year-out; station weather primary; AMENDMENT 2 scoring")
    print("=" * 78)
    print("\n1. Calibrated rates per year (own-year), relative to each driver's base values")
    for arm in ARMS:
        base = base_rates(arm)
        print(f"  {arm}  base " + ", ".join(f"{k} {base[k]:.6f}" for k in RATE_KEYS))
        for y in YEARS:
            r, f = P[arm]["own"][y]
            print(f"     {y}  J/I x{f:.3f}   DVRP x{r['DVRP'] / base['DVRP']:.3f}   "
                  f"DVRR x{r['DVRR'] / base['DVRR']:.3f}")
    print("\n2. Phenology prediction, leave-one-year-out (model - observed, days)")
    for arm in ARMS:
        for y in YEARS:
            e = P[arm]["phen"][y]
            print(f"    {arm} {y}: " + "  ".join(
                f"{s.split('_')[0]} {'n/r' if v is None else f'{v:+.1f}'}" for s, v in e.items()))
    print(f"    C3 (port; |anthesis| <= 5 d, |maturity| <= 7 d every year): {'adequate' if c3 else 'NOT adequate'}")

    print("\n3. Stage-matched ratios, LOYO, primary (observed | A | B)")
    labels = {"Q1": "GL(s)/GL(anthesis)", "Q2": "SL/(GL+SL) [floor]", "Q3": "GL/(TOT-SL)",
              "Q4": "(TOT-SL)(s)/(TOT-SL)(maturity)", "Q5": "root/(total-SL) (1999)"}
    for g in K.GROUPS:
        print(f"  {g} {labels[g]}" + ("   [anthesis point reported, not scored]" if g == "Q2" else ""))
        for y in YEARS:
            qo, qa, qb = (P["A"]["detail"][y]["qo"][g], P["A"]["detail"][y]["qm"][g],
                          P["B"]["detail"][y]["qm"][g])
            for st in qo:
                print(f"    {y} {st:<19} {qo[st]:7.3f} | {qa[st]:7.3f} | {qb[st]:7.3f}")

    print("\n4. log10 RMSE (bias), primary")
    hdr = "".join(f"{g:>16}" for g in SCORE_GROUPS) + f"{'POOLED':>9}"
    print("              " + hdr)
    for arm in ARMS:
        r, b = P[arm]["rmse"], P[arm]["bias"]
        cells = []
        for g in SCORE_GROUPS:
            bb = b["Q2"] if g == "Q2x" else b[g]
            cells.append(f"  {r[g]:6.3f}" + (f" ({bb:+.3f})" if g != "Q2x" else "         "))
        print(f"    {arm} LOYO   " + "".join(cells) + f"  {r['POOLEDx']:6.3f}")
    for arm in ARMS:
        r = P[arm]["own_rmse"]
        print(f"    {arm} own-yr " + "".join(f"  {r[g]:6.3f}         " for g in SCORE_GROUPS)
              + f"  {r['POOLEDx']:6.3f}")
    for arm in ARMS:
        r = P[arm]["cal"]
        print(f"    {arm} calendar" + "".join(f" {r[g]:6.3f}          " for g in SCORE_GROUPS)
              + f" {r['POOLEDx']:6.3f}")

    print("\n5. Decision (AMENDMENT 2: paired bootstrap of POOLEDx(B) - POOLEDx(A), 95 % interval)")
    for name, R in res.items():
        c = R["c1"]
        print(f"    {name:<13} A {R['A']['rmse']['POOLEDx']:.3f} vs B {R['B']['rmse']['POOLEDx']:.3f}"
              f"  B-A {c['delta']:+.3f} [{c['lo']:+.3f}, {c['hi']:+.3f}] -> {c['verdict']}"
              f"; port matures every year: {c['reached']}")
    c = P["c1_own"]
    print(f"    own-year rates: B-A {c['delta']:+.3f} [{c['lo']:+.3f}, {c['hi']:+.3f}] -> {c['verdict']}")
    print("    " + (f"verdict differs under {', '.join(flips)}" if flips else "same verdict under every sensitivity"))
    print("    (C2 is void under AMENDMENT 2: its premise -- a dead-leaf ERROR -- is removed by the floor)")
    print("   POST-HOC, Q2 read as 'attached' (no blade had detached; Q2 scored symmetrically):")
    for name, c in out["c1_attached"].items():
        print(f"    {name:<13} B-A {c['delta']:+.3f} [{c['lo']:+.3f}, {c['hi']:+.3f}] -> {c['verdict']}")
    return out


if __name__ == "__main__":
    main()
