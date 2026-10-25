#!/usr/bin/env python3
# =============================================================================
# validation/oryza2000_ir72ds92.py
# -----------------------------------------------------------------------------
# ORYZA2000 (IR72.D92) PORT vs THE SHIPPED oryza_growth, ON THE IR72DS92 EXPERIMENT
# (handoff step 5, docs/HANDOFF_mass_basis_refit.md)
#
# PRE-REGISTRATION -- written and pushed before either model was run on the 1992 weather
# or compared with the IR72DS92 measurements. (Already done before this file: the port
# was run on SYNTHETIC weather to check it runs and that its carbon balance closes; the
# IR72DS92 measurements were transcribed, so their values have been seen.)
#
# Why: the per-organ biomass M_k(t) drives every PFAS concentration (growth dilution,
# terminal accumulation). The default driver, src/oryza_growth.py, is a REDUCED Level-1
# re-coding whose parameters are "IR72 order of magnitude" -- not IR72.D92 (e.g. stems
# keep 25-30% of post-flowering growth vs 0 in IR72.D92; maturity 116 d after
# transplanting vs 98 d observed). src/oryza2000.py ports ORYZA2000 potential production
# line by line from the manual's Fortran with IR72.D92 verbatim. This file asks how each
# reproduces the experiment the ORYZA2000 manual ships with IR72.D92.
#
# Data: docs/literature_db/raw_si/oryza_ir72ds92_obs.csv (manual Section 7.2, IR72DS2.T92):
# IRRI 1992 dry season, emergence doy 4, transplanted doy 16, 25 hills x 5 plants, 225 kg
# N/ha, irrigated. NOTE the framing: the manual uses this very file as the example input
# of its calibration programs (DRATES/PARAM) and IR72.D92 was "derived from various
# experiments at IRRI". For the PORT this is therefore a FIDELITY test (a faithful port
# with the right weather should reproduce its own calibration experiment), NOT an
# independent validation of ORYZA2000. For the SHIPPED driver it measures how far it is
# from ORYZA2000-IR72 behaviour. The run was N-balance/water-balance with ample inputs,
# i.e. near-potential, not strictly potential.
#
# Weather: the station file PHIL1.992 is not obtainable here (requested; login-gated).
#   PRIMARY until it arrives: NASA POWER daily (raw_si/nasa_power_losbanos_1992.csv)
#   corrected to the station by the 12 days the manual prints (raw_si/
#   irri_phil1_1992_printed.csv): Tmin +1.17 C, Tmax +1.67 C (additive), radiation x1.211
#   (ratio of means). SENSITIVITY: temperature-only correction; no correction. All three
#   are reported. If PHIL1.992 becomes available, the station run SUPERSEDES these (same
#   rules, same scoring) and the POWER runs stay as the record of the fallback.
#
# Arms:
#   A  oryza2000.simulate, IR72.D92, IR72DS92 management, CO2 340, from emergence (doy 4).
#   B  oryza_growth (shipped OryzaParams, lat 14.18), started at transplanting (sow_doy 16,
#      its day 0 = transplanting), the same weather from doy 16, anchored to its 1740 g/m2
#      shoot target (scale_to_anchor=True: the configuration the PFAS default runs on).
#      Its dead leaves are accumulated from its own death rate (drlv * WLV) and scaled by
#      the same anchor factor. SENSITIVITY: B unanchored.
#
# Scoring: main-field dates only (doy 34, 58, 68, 83, 97, 114) -- the doy-16 record is the
# seedbed/transplanting transient (ORYZA writes the pre-dilution state on that day).
# Variables WLVG, WLVD, WST, WSO, WAGT (= green + dead leaves + stem + panicle) and LAI,
# zero observations excluded; per variable log10 RMSE. Phenology: day of year at which DVS
# first reaches 0.65 / 1 / 2 vs observed 58 / 83 / 114.
#
# Decision rules (on the PRIMARY weather):
#   R1 (fidelity of the port) ADEQUATE if |d flowering| <= 5 d, |d maturity| <= 5 d,
#      log10 RMSE(WAGT) <= 0.10 and log10 RMSE <= 0.15 for each of WLVG, WST, WSO.
#      Otherwise NOT ADEQUATE; if it passes under another weather treatment, the verdict
#      is recorded as weather-limited.
#   R2 (A vs B) A fits IR72DS92 better if its |d maturity| is smaller AND its log10 RMSE is
#      smaller on at least 3 of the 4 series WLVG, WST, WSO, WAGT.
#   R3 This file does NOT change any default. Whether biomass="oryza2000" replaces "oryza"
#      as the default is decided with the user after seeing these results (and then needs
#      the ORYZA-based re-fits redone).
#   POST-HOC (labelled, not part of the rules): the size of the PFAS BAF change between the
#      two drivers under the PFAS default climatology -- what a switch would move.
#
# EXPECTATION (stated so it can be wrong): A is ADEQUATE under the primary weather
#   (phenology within ~4 d, WAGT within ~15%); under uncorrected POWER its phenology runs a
#   few days late. B matures ~15-20 d late, keeps growing its stem after flowering (WST
#   over at doy 97/114) and under-fills its panicles before the anchor scale; A beats B by
#   R2.
#
#   python validation/oryza2000_ir72ds92.py
# =============================================================================
from __future__ import annotations

import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))

import oryza2000 as O2K      # noqa: E402
import oryza_growth as OG    # noqa: E402

RAW = os.path.join(ROOT, "docs", "literature_db", "raw_si")
SCORED_DOY = (34, 58, 68, 83, 97, 114)
SERIES = ("WLVG", "WLVD", "WST", "WSO", "WAGT", "LAI")
R2_SERIES = ("WLVG", "WST", "WSO", "WAGT")
OBS_EVENTS = {"panicle_initiation": 58, "flowering": 83, "maturity": 114}
TREATMENTS = ("corrected", "temp_only", "uncorrected")
TRANSPLANT_DOY = 16


def _rows(name):
    with open(os.path.join(RAW, name)) as f:
        return list(csv.DictReader(line for line in f if not line.startswith("#")))


def load_obs():
    obs, events = {}, {}
    for r in _rows("oryza_ir72ds92_obs.csv"):
        if r["variable"].startswith("event_"):
            events[r["variable"][6:]] = int(r["doy"])
        else:
            obs.setdefault(r["variable"], {})[int(r["doy"])] = float(r["value"])
    return obs, events


def station_bias():
    """Station - POWER over the 12 printed station days: (dTmin, dTmax, radiation ratio)."""
    power = {int(r["doy"]): r for r in _rows("nasa_power_losbanos_1992.csv") if r["date"].startswith("1992")}
    st = _rows("irri_phil1_1992_printed.csv")
    dn = np.mean([float(r["tmin"]) - float(power[int(r["doy"])]["tmin"]) for r in st])
    dx = np.mean([float(r["tmax"]) - float(power[int(r["doy"])]["tmax"]) for r in st])
    rr = (np.mean([float(r["rad_kj"]) / 1000.0 for r in st])
          / np.mean([float(power[int(r["doy"])]["rad_mj"]) for r in st]))
    return float(dn), float(dx), float(rr)


def load_weather(treatment="corrected"):
    """1992 daily weather from doy 1 (NASA POWER, corrected to the station as requested)."""
    rows = sorted((r for r in _rows("nasa_power_losbanos_1992.csv")
                   if r["date"] <= "1992-05-31"), key=lambda r: int(r["doy"]))
    rad = np.array([float(r["rad_mj"]) for r in rows])
    tmin = np.array([float(r["tmin"]) for r in rows])
    tmax = np.array([float(r["tmax"]) for r in rows])
    dn, dx, rr = station_bias()
    if treatment in ("corrected", "temp_only"):
        tmin, tmax = tmin + dn, tmax + dx
    if treatment == "corrected":
        rad = rad * rr
    return dict(doy0=int(rows[0]["doy"]), rad_mj=rad, tmin=tmin, tmax=tmax)


def run_port(weather):
    r = O2K.simulate(weather, O2K.load_crop(), O2K.Management())
    doy = r["doy"].astype(int)
    out = {"doy": doy, "DVS": r["dvs"], "WLVG": r["wlvg"], "WLVD": r["wlvd"], "WST": r["wst"],
           "WSO": r["wso"], "WAGT": r["wagt"], "LAI": r["lai"], "cbc_max": r["cbc_max"]}
    return out


def run_shipped(weather, anchored=True):
    p = OG.OryzaParams(lat_deg=14.18, sow_doy=TRANSPLANT_DOY)
    n = int(np.ceil(p.season)) + 1
    k = TRANSPLANT_DOY - weather["doy0"]
    w = {key: weather[key][k:k + n] for key in ("tmax", "tmin", "rad_mj")}
    sim = OG.simulate_oryza(p, w)
    dead = np.concatenate([[0.0], np.cumsum(sim["drlv"] * sim["wlv"])[:-1]])
    f = 1.0
    if anchored:
        shoot_final = sim["wlv"][-1] + sim["wst"][-1] + sim["wso"][-1]
        f = p.shoot_anchor_g_m2 * 10.0 / shoot_final
    doy = (TRANSPLANT_DOY + sim["t"]).astype(int)
    return {"doy": doy, "DVS": sim["dvs"], "WLVG": sim["wlv"] * f, "WLVD": dead * f,
            "WST": sim["wst"] * f, "WSO": sim["wso"] * f,
            "WAGT": (sim["wlv"] + dead + sim["wst"] + sim["wso"]) * f, "LAI": sim["lai"],
            "anchor_factor": f}


def phenology(run):
    out = {}
    for name, lim in (("panicle_initiation", 0.65), ("flowering", 1.0), ("maturity", 2.0)):
        idx = np.where(run["DVS"] >= lim - 1e-9)[0]
        out[name] = int(run["doy"][idx[0]]) if len(idx) else None
    return out


def score(run, obs):
    res = {}
    for v in SERIES:
        err = []
        for d in SCORED_DOY:
            o = obs[v].get(d, 0.0)
            if o <= 0.0:
                continue
            j = np.where(run["doy"] == d)[0]
            m = float(run[v][j[0]]) if len(j) else float(run[v][-1])     # held after maturity
            err.append(np.log10(max(m, 1e-9)) - np.log10(o))
        res[v] = float(np.sqrt(np.mean(np.square(err)))) if err else float("nan")
    return res


def value_at(run, v, d):
    j = np.where(run["doy"] == d)[0]
    return float(run[v][j[0]]) if len(j) else float(run[v][-1])


def r1(ph, sc, obs_events):
    df = abs(ph["flowering"] - obs_events["flowering"]) if ph["flowering"] else 99
    dm = abs(ph["maturity"] - obs_events["maturity"]) if ph["maturity"] else 99
    return (df <= 5 and dm <= 5 and sc["WAGT"] <= 0.10
            and all(sc[v] <= 0.15 for v in ("WLVG", "WST", "WSO")))


def r2(phA, scA, phB, scB, obs_events):
    dA = abs(phA["maturity"] - obs_events["maturity"])
    dB = abs(phB["maturity"] - obs_events["maturity"])
    wins = sum(scA[v] < scB[v] for v in R2_SERIES)
    return dA < dB and wins >= 3, wins


def main(posthoc=True):
    obs, ev = load_obs()
    dn, dx, rr = station_bias()
    print("=" * 92)
    print("ORYZA2000 (IR72.D92) PORT vs SHIPPED oryza_growth on IR72DS92 -- pre-registered in 288de63")
    print("=" * 92)
    print(f"weather: NASA POWER corrected by the 12 printed station days: Tmin {dn:+.2f} C, "
          f"Tmax {dx:+.2f} C, radiation x{rr:.3f}")
    results = {}
    for tr in TREATMENTS:
        w = load_weather(tr)
        A = run_port(w)
        B = run_shipped(w, anchored=True)
        Bu = run_shipped(w, anchored=False)
        res = {}
        for lab, run in (("A", A), ("B", B), ("B_unanchored", Bu)):
            res[lab] = dict(run=run, ph=phenology(run), sc=score(run, obs))
        results[tr] = res

    for tr in TREATMENTS:
        res = results[tr]
        print(f"\n--- weather: {tr} ---")
        print(f"  {'':14}{'PI':>6}{'flower':>8}{'mature':>8}   " + "".join(f"{v:>7}" for v in SERIES)
              + "   (log10 RMSE, scored dates)")
        print(f"  {'observed':14}{ev['panicle_initiation']:>6}{ev['flowering']:>8}{ev['maturity']:>8}")
        for lab in ("A", "B", "B_unanchored"):
            ph, sc = res[lab]["ph"], res[lab]["sc"]
            print(f"  {lab:14}{ph['panicle_initiation'] or -1:>6}{ph['flowering'] or -1:>8}"
                  f"{ph['maturity'] or -1:>8}   " + "".join(f"{sc[v]:>7.3f}" for v in SERIES))

    prim = results["corrected"]
    A, B = prim["A"], prim["B"]
    print("\nPrimary weather -- values on the scored dates (kg/ha; LAI m2/m2):")
    print(f"  {'doy':>4} " + "".join(f"{v:>21}" for v in ("WLVG o/A/B", "WST o/A/B", "WSO o/A/B", "WAGT o/A/B")))
    for d in SCORED_DOY:
        cells = []
        for v in ("WLVG", "WST", "WSO", "WAGT"):
            cells.append(f"{obs[v][d]:>7.0f}{value_at(A['run'], v, d):>7.0f}{value_at(B['run'], v, d):>7.0f}")
        print(f"  {d:>4} " + "".join(cells))
    print(f"  port carbon balance, max relative error: {A['run']['cbc_max']:.1e}; "
          f"shipped anchor factor {B['run']['anchor_factor']:.3f}")

    v_r1 = {tr: r1(results[tr]["A"]["ph"], results[tr]["A"]["sc"], ev) for tr in TREATMENTS}
    v_r2, wins = r2(A["ph"], A["sc"], B["ph"], B["sc"], ev)
    verdict_r1 = ("ADEQUATE" if v_r1["corrected"] else
                  ("NOT ADEQUATE (weather-limited: passes under "
                   + ", ".join(t for t in TREATMENTS if v_r1[t]) + ")" if any(v_r1.values())
                   else "NOT ADEQUATE"))
    print("\n" + "=" * 92)
    print("VERDICT (against the pre-registered rules, primary weather)")
    print("=" * 92)
    print(f"  R1 port fidelity: {verdict_r1}   (by treatment: "
          + ", ".join(f"{t} {'pass' if v_r1[t] else 'fail'}" for t in TREATMENTS) + ")")
    print(f"  R2 port fits IR72DS92 better than the shipped driver: {'YES' if v_r2 else 'NO'} "
          f"(maturity error {abs(A['ph']['maturity'] - ev['maturity'])} vs "
          f"{abs(B['ph']['maturity'] - ev['maturity'])} d; lower RMSE on {wins}/4 series)")
    print("  R3 no default changed")
    post = posthoc_pfas() if posthoc else None
    return dict(results=results, r1=v_r1, r1_verdict=verdict_r1, r2=v_r2, r2_wins=wins,
                bias=(dn, dx, rr), posthoc=post)


POSTHOC_CONGENERS = ("PFBA", "PFOA", "PFOS", "PFDA", "PFDoDA")


def posthoc_pfas():
    """POST-HOC (not part of the rules): what switching the PFAS driver would move.

    Default PFAS scenario (Cwo 1, measured Q_TP, the repo's tropical climatology), harvest
    BAFs for the shipped 'oryza' at 120 d, the port at 120 d (the default season runs past its
    maturity) and the port harvested at its own maturity."""
    import model_api as api
    t = np.linspace(0.0, 120.0, 481)
    t_mat = float(t[np.argmax(O2K.organ_biomass(t, 120.0)["dvs"] >= 2.0)])
    g = dict(scenario="generic")                 # the climatology, not the calibration site-year
    cases = (("oryza, 120 d", dict(g, biomass="oryza", season=120.0)),
             ("oryza2000, 120 d", dict(g, biomass="oryza2000", season=120.0)),
             (f"oryza2000, {t_mat:.0f} d", dict(g, biomass="oryza2000", season=t_mat)))
    print("\n" + "=" * 92)
    print(f"POST-HOC -- PFAS harvest BAF [L/kg fw] under each driver (port matures at t = {t_mat:.0f} d "
          "in the PFAS climatology)")
    print("=" * 92)
    print(f"  {'':8}" + "".join(f"{lab:>27}" for lab, _ in cases) + "      ratio port/oryza (120 d)")
    print(f"  {'':8}" + "".join(f"{'root  straw  grain':>27}" for _ in cases))
    out = {}
    for c in POSTHOC_CONGENERS:
        row = []
        for lab, kw in cases:
            r = api.simulate(c, **kw)
            row.append((r["baf_final"]["root"], r["straw_baf"], r["baf_final"]["grain"]))
        out[c] = row
        ratio = [row[1][i] / row[0][i] for i in range(3)]
        print(f"  {c:8}" + "".join(f"{a:>9.3g}{b:>9.3g}{g:>9.3g}" for a, b, g in row)
              + f"      {ratio[0]:.2f} {ratio[1]:.2f} {ratio[2]:.2f}")
    return dict(t_maturity=t_mat, baf=out)


if __name__ == "__main__":
    main()
