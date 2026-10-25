#!/usr/bin/env python3
# =============================================================================
# validation/oryza2000_wuxi_partitioning.py
# -----------------------------------------------------------------------------
# A SECOND INDEPENDENT CHECK OF THE BIOMASS DRIVERS -- ORGAN PARTITIONING AT THE WUXI
# RICE FACE (Yang et al. 2006, Field Crops Res. 98:12-19; ambient-CO2 plots, 2001-2003)
#
# PRE-REGISTRATION -- written and pushed before either driver was run on Wuxi weather or
# compared with these data. (Already seen: the transcribed data -- docs/literature_db/raw_si/
# yang2006_wuxi_face_dm.csv -- and Hasegawa et al. 2017's conditions table; every Shizukuishi
# result: 5b, option C, the redo with the recorded nursery and events, and the Gao crop set.
# The weather had not been downloaded when this was written.)
#
# Why: the Shizukuishi data (Kim 2001/2003) have leaf blades, totals and roots only -- no
# stem or panicle weights -- so they cannot test the one thing the PFAS model takes most
# directly from a biomass driver: how dry matter is split between leaf, stem and panicle
# over the season (the grain is the PFAS terminal sink; the straw = stem + leaf). Yang 2006
# reports exactly that split (leaf / stem incl. sheath / spike, summing to 1) plus the
# above-ground dry matter, at five samplings in each of three years, for another japonica
# (Wuxiangjing 14, an early-maturing late japonica) at another site (31.6 N) and lab.
#
# What the data do NOT give (stated, not worked around): heading and maturity DATES (the
# modellers of Hasegawa et al. 2017 received them, but they are not published), roots, and
# whether dead leaves are inside "leaf" or outside the total. Hasegawa et al. 2017 Table S1
# gives Wuxi's mean growth duration, 158 d -- read as sowing to maturity, mean of the years:
# for Shizukuishi the same table's 146 d matches Sasaki's 1999 maturity sampling exactly and
# Kim's maturity samplings within 3-5 d (149-151 d) -- i.e. Wuxi maturity at ~132 days after
# transplanting; and the season mean temperature and radiation per year.
#
# Arms (ambient CO2 per year 354/364/375 umol/mol; 24 hills m-2, 3 plants/hill; sown 18 May,
# transplanted 13 June; lat 31.58 N):
#   A        the ORYZA2000 port, IR72.D92, its own development rates.
#   B        the shipped oryza_growth (anchored -- the anchor scales all organs equally, so
#            every scored ratio is anchor-free), its own development rates.
#   G        the port with the Gao 2021 japonica crop set (as in oryza2000_kim_face_gao_
#            japonica.py, Gao's own development rates).                       [secondary]
#   A+PART   IR72.D92 with Gao's partitioning tables only (FLVTB/FSTTB/FSOTB). [secondary]
#   Port seedbed: emergence 4 d after sowing (assumption) -> SBDUR 22 d; open seedbed
#   (TMPSB 0; late-May Wuxi needs no cover; the paper says only "raised under ambient air").
#
# Weather: NASA POWER daily for the Wuxi cell (31.58 N, 120.50 E; cell elevation ~6 m, site
#   5 m), 2001-03-01..2003-12-31, saved to raw_si/. CHECK fixed in advance: the POWER mean
#   temperature over transplanting -> DAT 132 is compared per year with Hasegawa's season
#   mean (25.2 / 24.7 / 25.4 oC); if any year differs by > 0.5 oC, a sensitivity run shifts
#   that year's Tmin/Tmax by the difference. (The primary run uses POWER as downloaded.)
#
# Placement -- by stage on each driver's OWN DVS (the samplings were "fixed to coincide with"
#   the stages): PI sample at DVS 0.65; heading at DVS 1.0 (heading taken as flowering, as
#   in the Shizukuishi files); mid-ripening at the driver's DVS-1 day + 20 d (2001, 2002) /
#   + 30 d (2003), the paper's "days after heading"; maturity at DVS 2.0. If a driver never
#   reaches DVS 2 before Dec 31, its maturity is read on Dec 31 and flagged. The tillering
#   sample (27-28 DAT) has no stage clock: reported on the calendar day, not scored.
#
# Scored (log10 errors, ambient plots, 3 years pooled; RMSE per group):
#   W1 leaf share  at PI, heading, mid-ripening, maturity
#   W2 spike share at heading, mid-ripening, maturity
#   W3 stem share  at PI, heading, mid-ripening, maturity
#   W4 dry-matter timing AGB(stage)/AGB(maturity) at PI, heading, mid-ripening
#   POOLED = mean of W1-W4.
#   Model organs: leaf = green + dead leaf blades (WLVG + WLVD), stem = WST (sheaths are in
#   ORYZA's stem, as in Yang's), spike = WSO (panicles), total = their sum (PRIMARY reading:
#   dead leaves counted, nothing left out of the total). SENSITIVITY reading: leaf = green
#   only, dead leaves outside the total.
#
# Designs:
#   W   PRIMARY. Own development rates (no phenology data exist to calibrate them).
#   Wp  SECONDARY. Season length pinned with the only phenology information there is: one
#       common factor on DVRJ/DVRI puts DVS 0.65 on the PI sampling day (51/47/45 DAT), and
#       one common factor on DVRP and DVRR together (their ratio kept) puts DVS 2.0 at
#       DAT 132 every year (Hasegawa's mean duration -- the same day for every year is an
#       assumption). Heading then falls where each driver's own DVRP/DVRR ratio puts it.
#   Sensitivities of W: SBDUR 18 / 26 d (port arms); the green-leaf-only reading; the
#   temperature offset (only if the CHECK triggers).
#
# Decision rules:
#   R1 (primary, W) -- the independent partitioning evidence SUPPORTS the port over the
#      shipped driver if ALL hold: (a) A POOLED <= B POOLED + 0.05; (b) A is not worse than
#      B by more than 0.05 on W2 (spike share -- the PFAS grain sink) or W3 (stem share --
#      most of the straw); (c) A reaches DVS 2 before Dec 31 in all three years.
#      Otherwise R1 FAILS, with the failing clause recorded.
#   R2 -- R1 re-evaluated on Wp and under each sensitivity; a flip is reported as
#      "sensitive to <x>".
#   R3 (secondary) -- Gao's japonica partitioning fits better than IR72's if
#      mean(W1, W2, W3) of A+PART <= that of A - 0.05 (design W). Reported for G as well.
#   R4 (descriptive) -- each driver's heading/maturity DAT next to the ~132-DAT maturity;
#      the port's above-ground DM at maturity next to the observed 1674-1829 g m-2 (the
#      shipped driver's absolute level is set by its anchor and is not compared).
#   R5 -- no default and no parameter file changes.
#
# EXPECTATION (stated so it can be wrong; from the IR72DS92 and Shizukuishi results and the
#   parameter tables, NOT from any Wuxi run):
#   * The shipped driver keeps growing its stem after flowering and under-fills panicles (it
#     did on IR72DS92: stem 7095 vs 4243, panicles 6912 vs 9843 kg/ha), so at Wuxi its spike
#     share at mid-ripening and maturity is too low and its stem share too high; the port is
#     closer on W2 and W3. R1 PASSES. The shipped driver does reach DVS 2 at Wuxi (warmer
#     than Shizukuishi), late.
#   * The port with IR72's own rates matures early -- before DAT 115, i.e. >= 15 d before the
#     ~132-DAT maturity (IR72 is a short-duration indica, Wuxiangjing 14 a photoperiod-
#     sensitive late japonica; IR72's PPSE is 0), and its above-ground DM at maturity falls
#     short of the observed.
#   * W1: under the primary (green + dead) reading the port's leaf share at maturity is above
#     the observed (its dead pool is large); W1 is no better for the port than for the
#     shipped driver.
#   * R3: Gao's partitioning is not better than IR72's by 0.05 (both start panicle growth
#     around DVS 0.75 and give a similar spike share by heading).
#   * R2: Wp agrees with W (no flip); the green-only reading does not flip R1.
#
#   python validation/oryza2000_wuxi_partitioning.py
#
# AMENDMENT 2: R1's fixed 0.05 margins are replaced by the paired-bootstrap rule of
# validation/oryza2000_kim_face.py (AMENDMENT 2): statistic POOLED(B) - POOLED(A) over the
# scored points, resampled within W1-W4; "A better" / "B better" / "indistinguishable" by
# the 95 % interval; the port must reach DVS 2. R3 (Gao vs IR72) uses the same rule on W1-W3.
# Expectation: R1 "A better", R3 "indistinguishable".
# =============================================================================
from __future__ import annotations

import csv
import datetime as _dt
import os
import sys

import numpy as np
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, HERE)

import oryza2000 as O2K                         # noqa: E402
import oryza_growth as OG                       # noqa: E402
import oryza2000_kim_face_phenology as C        # noqa: E402
import oryza2000_kim_face as K                  # noqa: E402
import oryza2000_kim_face_gao_japonica as GJ    # noqa: E402

RAW = os.path.join(ROOT, "docs", "literature_db", "raw_si")
YEARS = (2001, 2002, 2003)
LAT, NH, NPLH = 31.58, 24.0, 3.0
CO2 = {2001: 354.0, 2002: 364.0, 2003: 375.0}          # Hasegawa 2017 Table S1, daylight ambient
SOW = {y: _dt.date(y, 5, 18) for y in YEARS}
TRANSPLANT = {y: _dt.date(y, 6, 13) for y in YEARS}
SBDUR = 22                                             # emergence 4 d after sowing (assumption)
MATURITY_DAT = 132.0                                   # 158 d sowing -> maturity (Hasegawa S1) - 26
HASEGAWA_TMEAN = {2001: 25.2, 2002: 24.7, 2003: 25.4}
DAH_MID = {2001: 20, 2002: 20, 2003: 30}               # mid-ripening sampling, days after heading
SCORED = ("panicle_initiation", "heading", "mid_ripening", "maturity")
WGROUPS = ("W1", "W2", "W3", "W4")
ARMS = ("A", "B", "G", "A+PART")
RATE_KEYS = C.RATE_KEYS


def _rows(name):
    with open(os.path.join(RAW, name)) as f:
        return list(csv.DictReader(line for line in f if not line.startswith("#")))


def load_obs():
    """{year: {stage: {dat, AGB, leaf, stem, spike}}} (ambient plots; fractions sum to 1)."""
    obs = {y: {} for y in YEARS}
    for r in _rows("yang2006_wuxi_face_dm.csv"):
        if r["co2"] != "AMB":
            continue
        y = int(r["year"])
        sp = float(r["spike_frac"]) if r["spike_frac"] else 0.0
        d = {"AGB": float(r["aboveground_dm_g_m2"]), "leaf": float(r["leaf_frac"]),
             "stem": float(r["stem_frac"]), "spike": sp}
        if r["sampling"].endswith("DAT"):
            d["dat"] = float(r["sampling"].split()[0])
        obs[y][r["stage"]] = d
    return obs


def _doy(d):
    return d.timetuple().tm_yday


def load_weather(year, dT=0.0):
    rows = sorted((r for r in _rows("nasa_power_wuxi_2001_2003.csv") if r["date"].startswith(str(year))),
                  key=lambda r: int(r["doy"]))
    return dict(doy0=int(rows[0]["doy"]), rad_mj=np.array([float(r["rad_mj"]) for r in rows]),
                tmin=np.array([float(r["tmin"]) for r in rows]) + dT,
                tmax=np.array([float(r["tmax"]) for r in rows]) + dT)


def weather_check():
    """POWER mean temperature, transplanting -> DAT 132, vs Hasegawa's season mean."""
    out = {}
    for y in YEARS:
        w = load_weather(y)
        k = _doy(TRANSPLANT[y]) - w["doy0"]
        sl = slice(k, k + int(MATURITY_DAT) + 1)
        tm = float(np.mean(0.5 * (w["tmin"][sl] + w["tmax"][sl])))
        out[y] = dict(power=tm, hasegawa=HASEGAWA_TMEAN[y], diff=HASEGAWA_TMEAN[y] - tm,
                      rad=float(np.mean(w["rad_mj"][sl])))
    return out


def run_port(year, crop, sbdur=SBDUR, dT=0.0):
    w = load_weather(year, dT)
    tp = _doy(TRANSPLANT[year])
    m = O2K.Management(lat=LAT, emd=tp - sbdur, sbdur=sbdur, nplh=NPLH, nh=NH, co2=CO2[year])
    r = O2K.simulate(w, crop, m, max_days=len(w["rad_mj"]) - (m.emd - w["doy0"]))
    assert r["events"].get("transplanting") == tp
    return {"dat": r["doy"] - tp, "DVS": r["dvs"], "GL": r["wlvg"], "SL": r["wlvd"],
            "ST": r["wst"], "SO": r["wso"]}


def run_shipped(year, params_kw=None, dT=0.0):
    w = load_weather(year, dT)
    tp = _doy(TRANSPLANT[year])
    n = int(w["doy0"] + len(w["rad_mj"]) - tp)                     # transplanting .. Dec 31
    p = OG.OryzaParams(lat_deg=LAT, sow_doy=tp, season=float(n - 1), **(params_kw or {}))
    k = tp - w["doy0"]
    sim = OG.simulate_oryza(p, {key: w[key][k:k + n] for key in ("tmax", "tmin", "rad_mj")})
    dead = np.concatenate([[0.0], np.cumsum(sim["drlv"] * sim["wlv"])[:-1]])
    f = p.shoot_anchor_g_m2 * 10.0 / (sim["wlv"][-1] + sim["wst"][-1] + sim["wso"][-1])
    return {"dat": sim["t"], "DVS": sim["dvs"], "GL": sim["wlv"] * f, "SL": dead * f,
            "ST": sim["wst"] * f, "SO": sim["wso"] * f}


def base_rates(arm):
    return C.base_rates("B") if arm == "B" else GJ.base_rates("A" if arm == "A+PART" else arm)


def run(arm, year, rates=None, sbdur=SBDUR, dT=0.0):
    if arm == "B":
        return run_shipped(year, dict(rates) if rates else None, dT)
    c = GJ.crop(arm)
    if rates:
        c.update(rates)
    return run_port(year, c, sbdur, dT)


def stage_index(r, year):
    """Fractional day index of each sample on the driver's own DVS clock; `reached` flags."""
    ar = np.arange(len(r["dat"]))
    idx, reached = {}, {}
    for st, lim in (("panicle_initiation", 0.65), ("heading", 1.0), ("maturity", 2.0)):
        x = C.crossing(r, lim)
        reached[st] = x is not None
        idx[st] = float(np.interp(x, r["dat"], ar)) if x is not None else float(ar[-1])
    head_day = float(np.interp(idx["heading"], ar, r["dat"]))
    idx["mid_ripening"] = min(float(np.interp(head_day + DAH_MID[year], r["dat"], ar)), idx["maturity"])
    idx["tillering"] = None
    return idx, reached


def organs(r, i, green_only=False):
    g = lambda k: float(np.interp(i, np.arange(len(r[k])), r[k]))
    leaf = g("GL") if green_only else g("GL") + g("SL")
    st, so = g("ST"), g("SO")
    tot = leaf + st + so
    return {"AGB": tot, "leaf": leaf / tot, "stem": st / tot, "spike": so / tot}


def quantities(v):
    q = {g: {} for g in WGROUPS}
    for s in SCORED:
        q["W1"][s] = v[s]["leaf"]
        q["W3"][s] = v[s]["stem"]
    for s in ("heading", "mid_ripening", "maturity"):
        q["W2"][s] = v[s]["spike"]
    for s in ("panicle_initiation", "heading", "mid_ripening"):
        q["W4"][s] = v[s]["AGB"] / v["maturity"]["AGB"]
    return q


def evaluate(runs, obs, green_only=False):
    err = {g: [] for g in WGROUPS}
    detail = {}
    for y in YEARS:
        idx, reached = stage_index(runs[y], y)
        mv = {s: organs(runs[y], idx[s], green_only) for s in SCORED}
        qm, qo = quantities(mv), quantities(obs[y])
        pe = {g: {s: float(np.log10(qm[g][s]) - np.log10(o)) for s, o in qo[g].items()} for g in WGROUPS}
        for g in WGROUPS:
            err[g].extend(pe[g].values())
        dat = {s: float(np.interp(idx[s], np.arange(len(runs[y]["dat"])), runs[y]["dat"]))
               for s in ("heading", "maturity")}
        till = organs(runs[y], float(np.interp(obs[y]["tillering"]["dat"], runs[y]["dat"],
                                               np.arange(len(runs[y]["dat"])))), green_only)
        detail[y] = dict(idx=idx, reached=reached, qm=qm, qo=qo, mv=mv, dat=dat, tillering=till, err=pe)
    rmse = {g: float(np.sqrt(np.mean(np.square(e)))) for g, e in err.items()}
    bias = {g: float(np.mean(e)) for g, e in err.items()}
    rmse["POOLED"] = float(np.mean([rmse[g] for g in WGROUPS]))
    rmse["PART3"] = float(np.mean([rmse[g] for g in ("W1", "W2", "W3")]))
    return rmse, bias, detail


def calibrate_wp(arm, year, obs, sbdur=SBDUR):
    """Wp: common J/I factor -> DVS 0.65 on the PI sampling day; one factor on DVRP and DVRR
    together -> DVS 2.0 at DAT 132."""
    base = base_rates(arm)

    def solve(fun, lo, hi):
        for _ in range(6):
            if fun(lo) * fun(hi) < 0:
                return brentq(fun, lo, hi, xtol=1e-7, rtol=1e-7)
            lo, hi = lo / 2.0, hi * 2.0
        raise RuntimeError(f"no bracket ({arm}, {year})")

    def miss(rates, lim, day):
        x = C.crossing(run(arm, year, rates, sbdur), lim)
        return (x if x is not None else 1.0e3) - day

    f = solve(lambda x: miss({**base, "DVRJ": base["DVRJ"] * x, "DVRI": base["DVRI"] * x}, 0.65,
                             obs[year]["panicle_initiation"]["dat"]), 0.3, 3.0)
    r1 = {**base, "DVRJ": base["DVRJ"] * f, "DVRI": base["DVRI"] * f}
    k = solve(lambda x: miss({**r1, "DVRP": base["DVRP"] * x, "DVRR": base["DVRR"] * x}, 2.0,
                             MATURITY_DAT), 0.3, 3.0)
    return {**r1, "DVRP": base["DVRP"] * k, "DVRR": base["DVRR"] * k}, f, k


def errors(detail, groups=WGROUPS):
    return {g: [(y, st, detail[y]["err"][g][st]) for y in YEARS for st in detail[y]["err"][g]] for g in groups}


def r1(R, a="A", b="B", groups=WGROUPS):
    """AMENDMENT 2: paired bootstrap of POOLED(b) - POOLED(a) over the scored points + a reaches DVS 2."""
    c = K.compare(errors(R[a]["detail"], groups), errors(R[b]["detail"], groups), groups=groups)
    c["reached"] = all(R[a]["detail"][y]["reached"]["maturity"] for y in YEARS)
    return c


def design(obs, arms, rates=None, sbdur=SBDUR, dT=None, green_only=False):
    out = {}
    for arm in arms:
        runs = {y: run(arm, y, (rates or {}).get(arm, {}).get(y), sbdur if arm != "B" else SBDUR,
                       (dT or {}).get(y, 0.0)) for y in YEARS}
        rm, bi, de = evaluate(runs, obs, green_only)
        out[arm] = dict(runs=runs, rmse=rm, bias=bi, detail=de)
    return out


def main(verbose=True):
    obs = load_obs()
    wc = weather_check()
    offset_needed = any(abs(v["diff"]) > 0.5 for v in wc.values())
    res = {"W": design(obs, ARMS)}
    wp_rates = {arm: {} for arm in ("A", "B")}
    wp_fk = {arm: {} for arm in ("A", "B")}
    for arm in ("A", "B"):
        for y in YEARS:
            rr, f, k = calibrate_wp(arm, y, obs)
            wp_rates[arm][y], wp_fk[arm][y] = rr, (f, k)
    res["Wp"] = design(obs, ("A", "B"), rates=wp_rates)
    res["W_sbdur18"] = design(obs, ("A", "B"), sbdur=18)
    res["W_sbdur26"] = design(obs, ("A", "B"), sbdur=26)
    res["W_greenonly"] = design(obs, ("A", "B"), green_only=True)
    if offset_needed:
        res["W_Toffset"] = design(obs, ("A", "B"), dT={y: v["diff"] for y, v in wc.items()})
    verdict = {d: r1(R) for d, R in res.items()}
    W = res["W"]
    r3 = {arm: r1(W, a=arm, b="A", groups=("W1", "W2", "W3")) for arm in ("A+PART", "G")}
    flips = [d for d in verdict if d != "W" and verdict[d]["verdict"] != verdict["W"]["verdict"]]
    out = dict(obs=obs, weather=wc, offset_needed=offset_needed, res=res, verdict=verdict,
               r1=verdict["W"], r3=r3, flips=flips, wp_rates=wp_rates, wp_fk=wp_fk)
    if not verbose:
        return out

    print("=" * 78)
    print("WUXI FACE (Yang 2006) -- leaf / stem / spike partitioning, ambient plots, 2001-2003")
    print("=" * 78)
    print("\n0. Weather check: POWER mean T, transplanting -> DAT 132, vs Hasegawa S1 season mean")
    for y, v in wc.items():
        print(f"    {y}: POWER {v['power']:.2f} C vs {v['hasegawa']:.1f} (diff {v['diff']:+.2f}); POWER rad {v['rad']:.1f} MJ")
    print("    -> " + ("offset sensitivity RUN (a year differs by > 0.5 C)" if offset_needed else
                      "no year differs by > 0.5 C; no offset run"))

    print("\n1. Phenology on each driver's own rates (W) and Wp factors -- DAT of heading / maturity")
    for arm in ARMS:
        d = W[arm]["detail"]
        print(f"    {arm:<7} " + "  ".join(
            f"{y}: {d[y]['dat']['heading']:.0f} / " + (f"{d[y]['dat']['maturity']:.0f}" if d[y]["reached"]["maturity"] else "n/r")
            for y in YEARS) + f"   (maturity ~{MATURITY_DAT:.0f} DAT, Hasegawa mean)")
    for arm in ("A", "B"):
        print(f"    Wp {arm}: " + "  ".join(f"{y} J/I x{wp_fk[arm][y][0]:.3f} P,R x{wp_fk[arm][y][1]:.3f}"
                                         f" heading {res['Wp'][arm]['detail'][y]['dat']['heading']:.0f}" for y in YEARS))
    agb = [W["A"]["detail"][y]["mv"]["maturity"]["AGB"] / 10.0 for y in YEARS]
    print("    port A above-ground DM at maturity (g m-2): " + ", ".join(f"{v:.0f}" for v in agb)
          + "   observed: " + ", ".join(f"{obs[y]['maturity']['AGB']:.0f}" for y in YEARS))

    print("\n2. Shares by stage, design W (observed | A | B | G | A+PART)")
    for key, lab in (("leaf", "leaf"), ("stem", "stem"), ("spike", "spike")):
        print(f"  {lab}")
        for y in YEARS:
            for s in ("tillering",) + SCORED:
                o = obs[y][s][key]
                vals = [W[a]["detail"][y]["tillering"][key] if s == "tillering" else W[a]["detail"][y]["mv"][s][key]
                        for a in ARMS]
                print(f"    {y} {s:<19} {o:6.3f} | " + " | ".join(f"{v:6.3f}" for v in vals)
                      + ("   (calendar, not scored)" if s == "tillering" else ""))

    print("\n3. log10 RMSE (bias)")
    print("                     " + "".join(f"{g:>16}" for g in WGROUPS) + f"{'POOLED':>9}{'W1-3':>8}")
    for d, R in res.items():
        for arm in R:
            rr, b = R[arm]["rmse"], R[arm]["bias"]
            print(f"  {d:<12} {arm:<7}" + "".join(f"  {rr[g]:6.3f} ({b[g]:+.3f})" for g in WGROUPS)
                  + f"  {rr['POOLED']:6.3f}  {rr['PART3']:6.3f}")

    print("\n4. Decisions (AMENDMENT 2: paired bootstrap, 95 % interval)")
    for d, v in verdict.items():
        R = res[d]
        print(f"  R1 {d:<12} POOLED A {R['A']['rmse']['POOLED']:.3f} vs B {R['B']['rmse']['POOLED']:.3f}"
              f"  B-A {v['delta']:+.3f} [{v['lo']:+.3f}, {v['hi']:+.3f}] -> {v['verdict']}"
              f"; port matures every year: {v['reached']}")
    print("  R2: " + (f"verdict differs under {', '.join(flips)}" if flips else "same verdict under Wp and every sensitivity"))
    for arm, v in r3.items():
        print(f"  R3 {arm:<7} vs IR72 (W1-3): IR72-{arm} {v['delta']:+.3f} [{v['lo']:+.3f}, {v['hi']:+.3f}] -> "
              f"{v['verdict'].replace('A better', arm + ' better').replace('B better', 'IR72 better')}")
    return out


if __name__ == "__main__":
    main()
