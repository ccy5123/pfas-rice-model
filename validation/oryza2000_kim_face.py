#!/usr/bin/env python3
# =============================================================================
# validation/oryza2000_kim_face.py
# -----------------------------------------------------------------------------
# AN INDEPENDENT CHECK OF THE ORYZA2000 PORT (and the shipped oryza_growth) ON THE
# SHIZUKUISHI RICE FACE EXPERIMENT (Kim et al. 2001 New Phytol.; Kim et al. 2003 GCB)
#
# PRE-REGISTRATION -- written and pushed before either driver was run on the Shizukuishi
# weather or compared with these measurements. (Already done before this file: the
# measurements were transcribed in an earlier session, so their values have been seen;
# the weather was downloaded and its header read; NO model was run on it.)
#
# Why: validation/oryza2000_ir72ds92.py showed the port reproduces IR72DS92 and beats the
# shipped driver there -- but IR72DS92 is the ORYZA2000 manual's own calibration example,
# so for the port that was a FIDELITY test, not an independent one. Before deciding
# whether biomass="oryza2000" should become the default, this file asks whether the
# port's advantage survives on an experiment neither driver was built or tuned on:
# different cultivar (Akitakomachi, a temperate japonica), site (Shizukuishi, 39.6 N),
# years (1998-2000) and lab. (Bouman & van Laar 2006 Fig. 2b, the other organ-biomass
# series in the repo, is the SAME IR72DS92 experiment and is therefore not used.)
#
# Data (docs/literature_db/raw_si/rice_organ_biomass_literature.csv; ambient-CO2 plots,
# medium N = 8 g N m-2 in 1998, 9 in 1999-2000; 19.05 hills m-2, 3 seedlings per hill;
# transplanted 1998-05-21, 1999-05-20, 2000-05-22):
#   Kim 2003 Table 2   green leaf blade GL, senescent leaf blade SL (g m-2) and green LAI
#                      at 5 samplings per year -- exact, tabulated.
#   Kim 2003 Fig. 2    total biomass TOT (g m-2) at the same samplings -- digitised
#                      (+-40 g m-2); root inclusion implied by the Methods (Kim 2001 Table
#                      1's whole-plant total at PI, 425.9, matches the 1999 point, 435).
#   Kim 2001 Table 1   root and total dry matter at PI and maturity, 1999 -- exact; roots
#                      from a 30 x 17.5 x 20 cm soil block, washed, so root is a FLOOR.
# The five samplings are, in order, tillering, panicle initiation (PI), anthesis,
# mid-ripening and grain maturity (Kim 2003 abstract, read via OpenAlex). ORYZA leaves
# are leaf BLADES (sheaths are in stems), so GL ~ WLVG and SL ~ WLVD; SL is a FLOOR for
# WLVD because detached dead blades are not recovered.
#
# THE CULTIVAR PROBLEM, and how it is handled. Both drivers carry IR72 (tropical indica)
# phenology, so neither is expected to time Akitakomachi's stages; comparing on calendar
# days would score the cultivar, not the driver. So every observation is placed by its
# STAGE, not its date: each driver is read on the day its own DVS first reaches 0.65 (PI),
# 1.0 (anthesis), 2.0 (maturity); mid-ripening is read at the same FRACTION of the
# anthesis->maturity interval as the sampling sat in that year (1998 0.56, 1999 0.49,
# 2000 0.52 -- from the sampling DATs, not from any model). Tillering is NOT scored: it
# has no DVS definition, its totals (50-95 g m-2) carry +-40 g m-2 digitising error, and
# it is dominated by the seedbed assumptions below. Every scored quantity is a RATIO, so
# the shipped driver's anchor scaling and the absolute yield level (cultivar, N, CO2)
# drop out.
#
# Scored quantities (per year; MN, ambient CO2):
#   Q1  green-leaf time course   GL(s)/GL(anthesis), s = PI, mid-ripening, maturity     (9)
#   Q2  senescent share           SL/(GL+SL) at anthesis, mid-ripening, maturity         (9)
#   Q3  leaf share of the plant   GL/TOT at PI, anthesis, mid-ripening, maturity        (12)
#   Q4  dry-matter timing         TOT(s)/TOT(maturity), s = PI, anthesis, mid-ripening   (9)
#   Q5  root share                root/total at PI and maturity, 1999 (Kim 2001)          (2)
#   Model counterparts: GL = WLVG, SL = WLVD (port) / the shipped driver's dead leaves
#   accumulated from its own death rate (drlv * WLV, as in oryza2000_ir72ds92.py), TOT =
#   green + dead leaves + stems (incl. reserves) + panicles + roots.
#   Per group: log10 RMSE over its values. POOLED = the mean of the five group RMSEs (so
#   Q3's 12 values do not outweigh Q5's 2).
# Reported, NOT scored (cultivar-specific or not comparable):
#   phenology in days after transplanting (both drivers are IR72); the anthesis->maturity
#   duration (observed 50 / 35 / 42 d; the shipped DVRR 1.05e-3 vs IR72.D92's 1.784e-3
#   makes this differ between the drivers, but grain-fill thermal time is still a cultivar
#   trait); absolute masses; specific leaf area GLAI/GL; the port's spikelet sterility
#   (IR72's cold-sterility response applies at 39.6 N).
#
# Drivers:
#   A  oryza2000.simulate, IR72.D92 unchanged, Management(lat 39.63, nh 19.05, nplh 3,
#      co2 368 [ambient 1998-2000], emergence = transplanting - SBDUR, SBDUR 30 d, NPLSB
#      1000 and the rest of the seedbed as in IR72DS92). The nursery practice (Japanese
#      mat seedlings, often under cover) is NOT reported in the repo's transcription, so
#      SBDUR is an ASSUMPTION; run until DVS > 2 or the weather ends (Nov 30).
#   B  oryza_growth (shipped OryzaParams; lat_deg 39.63, sow_doy = transplanting, season =
#      days to Nov 30), started at transplanting as it always is, anchored (the PFAS-default
#      configuration; ratios do not depend on the anchor).
#   If a driver has not reached DVS 2 by Nov 30, its maturity (and mid-ripening) values are
#   read on Nov 30 and the failure is reported.
#
# Weather: NASA POWER daily for the Shizukuishi cell (raw_si/nasa_power_shizukuishi_
#   1998_2000.csv), UNCORRECTED -- PRIMARY. The cell sits at 461 m vs a ~200 m paddy and no
#   station days exist here to correct it, so: SENSITIVITY S1 = +1.7 C on Tmin and Tmax
#   (lapse rate 6.5 C/km to ~200 m, itself an assumption); S2 = SBDUR 20 d; S3 = SBDUR
#   40 d (S2/S3 move only A).
#
# Decision rules (PRIMARY treatment):
#   D1 (the decision this file exists for) -- the independent evidence SUPPORTS the switch
#      to biomass="oryza2000" if ALL hold:
#        (a) A's POOLED log10 RMSE <= B's POOLED + 0.05;
#        (b) A is not worse than B by more than 0.05 on Q1 or on Q2 (the green-leaf course
#            and senescence feed the PFAS leaf-loss flux, the term the driver most directly
#            sets);
#        (c) A reaches DVS 2 by Nov 30 in all three years.
#      Otherwise D1 FAILS: the IR72DS92 advantage does not transfer and the default stays,
#      with the reason recorded.
#   D2 Robustness: D1 is re-evaluated under S1, S2, S3; if its verdict flips under any of
#      them it is reported as "sensitive to <treatment>".
#   D3 Descriptive adequacy per group (both drivers): log10 RMSE <= 0.10 good, <= 0.20
#      fair, > 0.20 poor. Not a gate.
#   D4 This file changes NO default. A passing D1 means only that the switch (biomass=
#      "oryza2000" + a season tied to crop maturity + the ORYZA-based re-fits) can be
#      pre-registered as one step, for the user to decide.
#
# EXPECTATION (stated so it can be wrong; from reading the parameter tables, not from runs):
#   both drivers share ORYZA's structure and IR72-like partitioning, so on Q3 and Q5 they
#   land within 0.05 of each other; both over-predict the senescent share (Q2: SL is a
#   floor) and under-predict the green leaf left at maturity (Q1: model ~0.3 of the
#   anthesis value vs observed 0.41-0.59); A is better on Q4 at PI (B grows too slowly
#   early, as on IR72DS92) and on the anthesis->maturity duration (B's ~1.7x longer
#   grain-fill thermal time); B may not reach maturity by late October in at least one
#   year. D1 PASSES and does not flip under S1-S3.
#
#   python validation/oryza2000_kim_face.py
#
# ADDENDUM -- written AFTER the primary and S1-S3 results were seen, BEFORE the station run.
#   What was seen: on the primary (uncorrected) weather NEITHER driver reaches DVS 2 by
#   Nov 30 in any year, so D1(c) fails; under S1 the port matures and D1 passes (D2:
#   "sensitive to S1"). Also seen: the shipped driver has no leaf death before DVS 1, so
#   its senescent share at anthesis is exactly 0 and its Q2 log10 RMSE is infinite -- the
#   pre-registration did not foresee a zero prediction.
#   The header's "no station days exist here" was WRONG (not searched for): the JMA AMeDAS
#   station Shizukuishi (block 0221, 195 m, ~8 km from the POWER point) has daily Tmax/Tmin
#   for 1998-2000 (raw_si/jma_amedas_shizukuishi_1998_2000.csv), and a DEM gives 202 m at
#   the POWER point, confirming the ~200 m assumption. Over May-Oct the station runs
#   +1.1 / +1.45 / +1.2 C (1998/1999/2000) warmer than the cell, i.e. between the primary
#   treatment and S1.
#   POST-HOC treatment S4 = station Tmax/Tmin (JMA ']' insufficient-data days replaced by
#   POWER + that year's mean station-POWER offset) + POWER radiation (AMeDAS measures none);
#   same drivers, rules and scoring. POST-HOC Q2 variant: Q2 without the anthesis point
#   (where the shipped driver predicts exactly 0), so both drivers get a finite score.
#   How they are used: the PRE-REGISTERED verdict is the primary one and stays recorded as
#   such. S4 is reported beside it as the treatment closest to the observed weather; it is
#   post-hoc (its existence was learned after the results), so it cannot overturn the
#   pre-registered verdict by itself -- it informs the user's decision, which D4 leaves open.
# =============================================================================
from __future__ import annotations

import csv
import datetime as _dt
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))

import oryza2000 as O2K      # noqa: E402
import oryza_growth as OG    # noqa: E402

RAW = os.path.join(ROOT, "docs", "literature_db", "raw_si")
YEARS = (1998, 1999, 2000)
TRANSPLANT = {1998: _dt.date(1998, 5, 21), 1999: _dt.date(1999, 5, 20), 2000: _dt.date(2000, 5, 22)}
STAGES = ("tillering", "panicle_initiation", "anthesis", "mid_ripening", "maturity")  # Kim 2003 abstract
SCORED = ("panicle_initiation", "anthesis", "mid_ripening", "maturity")
DVS_AT = {"panicle_initiation": 0.65, "anthesis": 1.0, "maturity": 2.0}
LAT, NH, NPLH, CO2 = 39.63, 19.05, 3.0, 368.0
GROUPS = ("Q1", "Q2", "Q3", "Q4", "Q5")
TREATMENTS = {"primary": dict(dT=0.0, sbdur=30), "S1_lapse+1.7C": dict(dT=1.7, sbdur=30),
              "S2_sbdur20": dict(dT=0.0, sbdur=20), "S3_sbdur40": dict(dT=0.0, sbdur=40)}


def _rows(name):
    with open(os.path.join(RAW, name)) as f:
        return list(csv.DictReader(line for line in f if not line.startswith("#")))


def _dat(s):
    return int(re.search(r"DAT(\d+)", s).group(1))


def load_obs():
    """{year: {stage: {dat, GL, SL, GLAI, TOT}}} for ambient CO2, medium N, plus the Kim 2001
    1999 root/total pairs {stage: (root, total)}."""
    obs = {y: {} for y in YEARS}
    dat_stage = {}
    for r in _rows("rice_organ_biomass_literature.csv"):
        src = r["source"]
        m = re.match(r"Kim2003_GCB_Table2_(\d{4})_AMB", src)
        if m:
            y, dat = int(m.group(1)), _dat(r["DAT_or_date"])
            k = int(re.search(r"sampling (\d)", r["note"]).group(1))
            st = STAGES[k - 1]
            dat_stage[(y, dat)] = st
            key = {"green leaf blade (GLdw)": "GL", "senescent leaf blade (SLdw)": "SL",
                   "green leaf area index (GLAI)": "GLAI"}[r["organ"]]
            obs[y].setdefault(st, {"dat": dat})[key] = float(r["value"])
    for r in _rows("rice_organ_biomass_literature.csv"):
        m = re.match(r"Kim2003_GCB_Fig2_(\d{4})_AMB", r["source"])
        if m and r["N_level"].startswith("MN"):
            y, dat = int(m.group(1)), _dat(r["DAT_or_date"])
            obs[y][dat_stage[(y, dat)]]["TOT"] = float(r["value"])
    root = {}
    for r in _rows("rice_organ_biomass_literature.csv"):
        if r["source"] == "Kim2001_NewPhytol_Table1_AMB" and r["N_level"].startswith("MN"):
            st = "panicle_initiation" if "(PI)" in r["DAT_or_date"] else "maturity"
            root.setdefault(st, {})["root" if r["organ"] == "root" else "total"] = float(r["value"])
    return obs, {st: (v["root"], v["total"]) for st, v in root.items()}


def ripening_fraction(obs_y):
    a, mr, ma = (obs_y[s]["dat"] for s in ("anthesis", "mid_ripening", "maturity"))
    return (mr - a) / (ma - a)


def load_weather(year, dT=0.0):
    rows = [r for r in _rows("nasa_power_shizukuishi_1998_2000.csv") if r["date"].startswith(str(year))]
    rows.sort(key=lambda r: int(r["doy"]))
    return dict(doy0=int(rows[0]["doy"]), rad_mj=np.array([float(r["rad_mj"]) for r in rows]),
                tmin=np.array([float(r["tmin"]) for r in rows]) + dT,
                tmax=np.array([float(r["tmax"]) for r in rows]) + dT)


def _doy(d):
    return d.timetuple().tm_yday


def run_port(year, dT=0.0, sbdur=30):
    w = load_weather(year, dT)
    tp = _doy(TRANSPLANT[year])
    m = O2K.Management(lat=LAT, emd=tp - sbdur, sbdur=sbdur, nplh=NPLH, nh=NH, co2=CO2)
    r = O2K.simulate(w, O2K.load_crop(), m, max_days=len(w["rad_mj"]) - (m.emd - w["doy0"]))
    assert r["events"].get("transplanting") == tp
    return {"dat": r["doy"] - tp, "DVS": r["dvs"], "GL": r["wlvg"], "SL": r["wlvd"],
            "ROOT": r["wrt"], "TOT": r["wlvg"] + r["wlvd"] + r["wst"] + r["wso"] + r["wrt"],
            "LAI": r["lai"], "NSP": r["nsp"], "NGR": r["ngr"], "WSO": r["wso"]}


def run_shipped(year, dT=0.0, anchored=True):
    w = load_weather(year, dT)
    tp = _doy(TRANSPLANT[year])
    n = int(w["doy0"] + len(w["rad_mj"]) - tp)                     # transplanting .. Nov 30
    p = OG.OryzaParams(lat_deg=LAT, sow_doy=tp, season=float(n - 1))
    k = tp - w["doy0"]
    sim = OG.simulate_oryza(p, {key: w[key][k:k + n] for key in ("tmax", "tmin", "rad_mj")})
    dead = np.concatenate([[0.0], np.cumsum(sim["drlv"] * sim["wlv"])[:-1]])
    f = 1.0
    if anchored:
        f = p.shoot_anchor_g_m2 * 10.0 / (sim["wlv"][-1] + sim["wst"][-1] + sim["wso"][-1])
    return {"dat": sim["t"], "DVS": sim["dvs"], "GL": sim["wlv"] * f, "SL": dead * f,
            "ROOT": sim["wrt"] * f, "TOT": (sim["wlv"] + dead + sim["wst"] + sim["wso"] + sim["wrt"]) * f,
            "LAI": sim["lai"], "WSO": sim["wso"] * f}


def stage_index(run, frac_mr):
    """Fractional day index of each stage on the driver's own DVS clock; `reached` flags."""
    idx, reached = {}, {}
    for st, lim in DVS_AT.items():
        hit = np.where(run["DVS"] >= lim - 1e-9)[0]
        reached[st] = bool(len(hit))
        idx[st] = float(hit[0]) if len(hit) else float(len(run["DVS"]) - 1)
    idx["mid_ripening"] = idx["anthesis"] + frac_mr * (idx["maturity"] - idx["anthesis"])
    return idx, reached


def at(run, key, i):
    return float(np.interp(i, np.arange(len(run[key])), run[key]))


def quantities(v, root=None):
    """The scored ratios from per-stage values v[stage] = {GL, SL, TOT} (and Q5 pairs)."""
    q = {g: {} for g in GROUPS}
    for s in ("panicle_initiation", "mid_ripening", "maturity"):
        q["Q1"][s] = v[s]["GL"] / v["anthesis"]["GL"]
    for s in ("anthesis", "mid_ripening", "maturity"):
        q["Q2"][s] = v[s]["SL"] / (v[s]["GL"] + v[s]["SL"])
    for s in SCORED:
        q["Q3"][s] = v[s]["GL"] / v[s]["TOT"]
    for s in ("panicle_initiation", "anthesis", "mid_ripening"):
        q["Q4"][s] = v[s]["TOT"] / v["maturity"]["TOT"]
    if root is not None:
        for s, (r, t) in root.items():
            q["Q5"][s] = r / t
    return q


def model_values(run, idx):
    return {s: {k: at(run, k, idx[s]) for k in ("GL", "SL", "TOT", "ROOT")} for s in SCORED}


def evaluate(runs, obs, root_obs):
    """Per-group log10 errors pooled over years; `runs` = {year: run}."""
    err = {g: [] for g in GROUPS}
    detail = {}
    for y in YEARS:
        frac = ripening_fraction(obs[y])
        idx, reached = stage_index(runs[y], frac)
        mv = model_values(runs[y], idx)
        root_m = None
        if y == 1999:
            root_m = {s: (mv[s]["ROOT"], mv[s]["TOT"]) for s in root_obs}
        qm = quantities(mv, root_m)
        qo = quantities({s: obs[y][s] for s in SCORED}, root_obs if y == 1999 else None)
        for g in GROUPS:
            for s, o in qo[g].items():
                err[g].append(np.log10(qm[g][s]) - np.log10(o))
        detail[y] = dict(idx=idx, reached=reached, qm=qm, qo=qo, mv=mv)
    rmse = {g: float(np.sqrt(np.mean(np.square(e)))) for g, e in err.items()}
    rmse["POOLED"] = float(np.mean([rmse[g] for g in GROUPS]))
    bias = {g: float(np.mean(e)) for g, e in err.items()}
    return rmse, bias, detail


def d1(ra, rb, detail_a):
    a = ra["POOLED"] <= rb["POOLED"] + 0.05
    b = (ra["Q1"] <= rb["Q1"] + 0.05) and (ra["Q2"] <= rb["Q2"] + 0.05)
    c = all(detail_a[y]["reached"]["maturity"] for y in YEARS)
    return dict(a=a, b=b, c=c, passed=a and b and c)


def adequacy(x):
    return "good" if x <= 0.10 else ("fair" if x <= 0.20 else "poor")


def run_treatment(obs, root_obs, dT, sbdur):
    ra_runs = {y: run_port(y, dT, sbdur) for y in YEARS}
    rb_runs = {y: run_shipped(y, dT) for y in YEARS}
    ra, ba, da = evaluate(ra_runs, obs, root_obs)
    rb, bb, db = evaluate(rb_runs, obs, root_obs)
    return dict(A=dict(rmse=ra, bias=ba, detail=da, runs=ra_runs),
                B=dict(rmse=rb, bias=bb, detail=db, runs=rb_runs), d1=d1(ra, rb, da))


def _fmt_stage_row(name, d, run):
    out = []
    for st in ("panicle_initiation", "anthesis", "maturity"):
        i = d["idx"][st]
        dat = at(run, "dat", i)
        out.append(f"{dat:5.0f}{'' if d['reached'][st] else '*'}")
    return f"    {name:<8}" + "  ".join(out)


def main(verbose=True):
    obs, root_obs = load_obs()
    res = {name: run_treatment(obs, root_obs, **kw) for name, kw in TREATMENTS.items()}
    if not verbose:
        return dict(res=res, obs=obs, root_obs=root_obs)
    P = res["primary"]
    print("=" * 78)
    print("ORYZA2000 port (A) vs shipped oryza_growth (B) -- Shizukuishi FACE, Akitakomachi,")
    print("ambient CO2, medium N, 1998-2000; observations placed by STAGE on each driver's DVS")
    print("=" * 78)
    print("\n1. Phenology on each driver's own clock (days after transplanting; * = not reached")
    print("   by Nov 30). REPORTED, NOT SCORED -- both drivers carry IR72 phenology.")
    print("    year      PI  anthesis  maturity | observed sampling DAT (PI / anthesis / maturity)")
    for y in YEARS:
        o = obs[y]
        print(f"  {y}  observed {o['panicle_initiation']['dat']:5d} {o['anthesis']['dat']:5d} "
              f"{o['maturity']['dat']:5d}")
        for arm in ("A", "B"):
            print(_fmt_stage_row(arm, P[arm]["detail"][y], P[arm]["runs"][y]))
    print("\n   anthesis -> maturity duration [d]:")
    for y in YEARS:
        o = obs[y]
        dur_o = o["maturity"]["dat"] - o["anthesis"]["dat"]
        s = []
        for arm in ("A", "B"):
            d = P[arm]["detail"][y]
            s.append(f"{arm} {d['idx']['maturity'] - d['idx']['anthesis']:.0f}"
                     f"{'' if d['reached']['maturity'] else '*'}")
        print(f"    {y}: observed {dur_o}   " + "   ".join(s))

    print("\n2. Scored ratios, primary weather (observed | A | B), per year")
    labels = {"Q1": "GL(s)/GL(anthesis)", "Q2": "SL/(GL+SL)", "Q3": "GL/TOT",
              "Q4": "TOT(s)/TOT(maturity)", "Q5": "root/total (1999)"}
    for g in GROUPS:
        print(f"  {g} {labels[g]}")
        for y in YEARS:
            qo, qa, qb = (P["A"]["detail"][y]["qo"][g], P["A"]["detail"][y]["qm"][g],
                          P["B"]["detail"][y]["qm"][g])
            for s in qo:
                print(f"    {y} {s:<19} {qo[s]:7.3f} | {qa[s]:7.3f} | {qb[s]:7.3f}")

    print("\n3. log10 RMSE per group (bias = mean log10 model/observed), primary weather")
    print("        " + "".join(f"{g:>16}" for g in GROUPS) + f"{'POOLED':>9}")
    for arm in ("A", "B"):
        r, b = P[arm]["rmse"], P[arm]["bias"]
        print(f"    {arm}   " + "".join(f"  {r[g]:.3f} ({b[g]:+.3f})" for g in GROUPS)
              + f"  {r['POOLED']:.3f}")
    for arm in ("A", "B"):
        print(f"    {arm} adequacy (D3): " + ", ".join(f"{g} {adequacy(P[arm]['rmse'][g])}"
                                                 for g in GROUPS))

    print("\n4. D1 / D2 -- does the independent evidence support the switch?")
    for name, R in res.items():
        d = R["d1"]
        print(f"    {name:<15} pooled A {R['A']['rmse']['POOLED']:.3f} vs B {R['B']['rmse']['POOLED']:.3f}"
              f" | (a) {d['a']} (b) {d['b']} (c) {d['c']} -> {'PASS' if d['passed'] else 'FAIL'}")
    v0 = P["d1"]["passed"]
    flips = [n for n, R in res.items() if R["d1"]["passed"] != v0]
    print(f"\n  VERDICT (primary): D1 {'PASSES' if v0 else 'FAILS'}"
          + (f"; sensitive to {', '.join(flips)}" if flips else "; unchanged under S1-S3"))

    print("\n5. Reported, not scored")
    print("   absolute (g m-2): observed TOT(maturity) | A | B(anchored);  GL(anthesis) obs | A | B")
    for y in YEARS:
        o = obs[y]
        va, vb = P["A"]["detail"][y]["mv"], P["B"]["detail"][y]["mv"]
        print(f"    {y}: {o['maturity']['TOT']:6.0f} | {va['maturity']['TOT'] / 10:6.0f} | "
              f"{vb['maturity']['TOT'] / 10:6.0f}     {o['anthesis']['GL']:5.0f} | "
              f"{va['anthesis']['GL'] / 10:5.0f} | {vb['anthesis']['GL'] / 10:5.0f}")
    print("   specific green leaf area GLAI/GL at anthesis (m2 g-1): observed | A | B")
    for y in YEARS:
        o = obs[y]["anthesis"]
        s = []
        for arm in ("A", "B"):
            d, run = P[arm]["detail"][y], P[arm]["runs"][y]
            i = d["idx"]["anthesis"]
            s.append(f"{at(run, 'LAI', i) / (at(run, 'GL', i) / 10):.4f}")
        print(f"    {y}: {o['GLAI'] / o['GL']:.4f} | " + " | ".join(s))
    print("   port spikelet fertility (grains / spikelets at the end of the run):")
    for y in YEARS:
        run = P["A"]["runs"][y]
        print(f"    {y}: {run['NGR'][-1] / max(run['NSP'][-1], 1e-9):.2f}  "
              f"(spikelets {run['NSP'][-1] / 1e4:.0f} m-2)")
    return dict(res=res, obs=obs, root_obs=root_obs)


if __name__ == "__main__":
    main()
