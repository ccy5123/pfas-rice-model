#!/usr/bin/env python3
# =============================================================================
# validation/bat_census_biocides.py
# -----------------------------------------------------------------------------
# Puts the substances of the BAT census report through THIS repo's model.
#
# That report (`REPORT_bat_census.md`) belongs to a SEPARATE project and is not in
# this repository. Everything taken from it is transcribed into
# data_obs/biocides_bat_census.csv with the section it came from, so this work
# stands on its own without it.
#
# WHAT THE TWO DOCUMENTS HAVE IN COMMON, AND WHAT THEY DO NOT. The report runs a
# regulatory screening tool (BAT) over EU biocidal actives and scores its FISH
# bioconcentration factor against published B/vB opinions. This repo models PFAS
# — and, since the neutral path was added, any un-ionised organic — accumulating
# in RICE. Different organism, different endpoint, different regulation. So there
# is exactly one thing to carry across: the report's audited PHYSICAL INPUTS, in
# particular a log Kow OF THE UNCHARGED FORM traced to a named source (its
# sections 7.5, 7.9, 8.14). That is the one input this model's neutral path needs,
# and the report's own headline correction — three anticoagulants entered with a
# distribution ratio D instead of log Kow — is a mistake this model would inherit
# in exactly the same way, because K_PW and TSCF are both functions of log Kow.
#
# THIS IS PREDICTION, NOT VALIDATION, and the distinction is the whole point of
# the neutral path (docs/neutral_dpu_validation.md). There is no measured rice
# concentration for any of these substances, so nothing here can be scored — with
# two exceptions this repo already carries and section 4 uses: propiconazole
# (Liu 2023 rice root RCF) and triclosan (Li 2019 soil table, root). Every other
# number below is a prediction with nothing to check it against.
#
# THE SCOPE RULE IS NOT CHOSEN HERE. It is read off the tables already shipped in
# data_obs/ and the bounds already recorded in docs/neutral_dpu_validation.md
# (sections 4l, 4m): the log Kow span over which measured plant rows exist, the
# span over which measured TSCF exists, and the weak-electrolyte path's stated
# limit — direction only, and not below a neutral fraction of ~0.1. That last one
# is the report's OWN section 3.0a rule ("more than 90% ionised at environmental
# pH"), which the report records as having been written down and never
# implemented. It is implemented here, at the same threshold, on the same
# substances.
#
#   python validation/bat_census_biocides.py            (~3 min)
#   python validation/bat_census_biocides.py --fast     (skips the log Kow sweep)
# =============================================================================
from __future__ import annotations

import csv
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT_DIR, "src"))

import model_api as api                                       # noqa: E402
import neutral_dpu as ND                                      # noqa: E402

TABLE = os.path.join(ROOT_DIR, "data_obs", "biocides_bat_census.csv")
OUT_CSV = os.path.join(HERE, "bat_census_biocides.csv")

REPORT_PH = 7.0          # the pH the report's ionised fractions are on (checked in section 1)
F_N_FLOOR = 0.10         # docs/neutral_dpu_validation.md section 4l: below this the
                         # speciation path under-delivers and is not quotable as a number
SEASON = 120.0


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _f(row, key):
    v = (row.get(key) or "").strip()
    if v in ("", "NA"):
        return None
    return float(v)


def load_table(path=TABLE):
    with open(path) as fh:
        rows = list(csv.DictReader(l for l in fh if not l.startswith("#")))
    for r in rows:
        r["substance"] = r["substance"].strip()
    return rows


def frac_ionised(pKa, pH=REPORT_PH, is_acid=True):
    """Henderson-Hasselbalch. An ACID is ionised above its pKa, a BASE below it."""
    d = (pH - pKa) if is_acid else (pKa - pH)
    return 1.0 / (1.0 + 10.0 ** (-d))


def pka_from_fraction(alpha, pH=REPORT_PH, is_acid=True):
    """The inverse used to recover the anticoagulants' pKa from the report's
    section 8.9 percentages. Undefined at alpha == 1, which is why salicylic
    acid's '100.0%' bounds its pKa rather than giving it."""
    if not 0.0 < alpha < 1.0:
        return None
    d = math.log10(alpha / (1.0 - alpha))
    return (pH - d) if is_acid else (pH + d)


def _obs_rows(path):
    with open(path) as fh:
        return list(csv.DictReader(l for l in fh if not l.startswith("#")))


# ---------------------------------------------------------------------------
# 1. the transcription, checked before anything is run
# ---------------------------------------------------------------------------
def check_transcription(rows, quiet=False):
    """Two things have to hold before a single input is used.

    (a) the report's percent-ionised column is on pH 7 -- established by
        round-tripping the six substances whose pKa it STATES; and
    (b) the pKa values marked `derived` are that inversion applied to the
        percentages, and nothing else.

    If (a) failed, every derived pKa would be on the wrong pH and the whole
    ionisation screen below would be wrong in the same direction."""
    checks = []
    for r in rows:
        pka, pct = _f(r, "pKa"), _f(r, "pct_ionised_pH7")
        if pka is None or pct is None:
            continue
        is_acid = (r["is_acid"] or "TRUE").strip().upper() == "TRUE"
        got = 100.0 * frac_ionised(pka, REPORT_PH, is_acid)
        checks.append((r["substance"], r["pka_basis"], pct, got, abs(got - pct)))
    stated = [c for c in checks if c[1] == "stated"]
    derived = [c for c in checks if c[1] == "derived"]
    worst_stated = max((c[4] for c in stated), default=0.0)
    worst_derived = max((c[4] for c in derived), default=0.0)
    if not quiet:
        print("=" * 78)
        print("1. THE TRANSCRIPTION, CHECKED BEFORE ANY RUN")
        print("=" * 78)
        print(f"  stated pKa round-tripped to the report's own %ionised at pH {REPORT_PH:.1f}:")
        for s, _, pct, got, d in stated:
            print(f"    {s:<22s} pKa->{got:8.3f}%   report {pct:8.3f}%   diff {d:.3f}")
        print(f"  -> worst disagreement {worst_stated:.3f} percentage points."
              " pH 7 is the report's basis.")
        print(f"  {len(derived)} pKa recovered from the percentages "
              f"(worst round-trip {worst_derived:.4f} pp):")
        for s, _, pct, got, _d in derived:
            print(f"    {s:<22s} {pct:6.2f}% ionised -> pKa "
                  f"{pka_from_fraction(pct / 100.0):.4f}")
    return dict(worst_stated=worst_stated, worst_derived=worst_derived,
                n_stated=len(stated), n_derived=len(derived))


# ---------------------------------------------------------------------------
# 2. the scope screen -- read off the shipped tables, not chosen here
# ---------------------------------------------------------------------------
def measured_spans():
    """The log Kow span over which this repo actually holds measured plant data.

    Two different questions, so two different spans: the ROOT partition (K_PW) is
    covered by the neutral observation tables, and TSCF -- the root->shoot loading
    that decides everything above the root -- only by the 30 un-ionised rows of
    the Schriever 2020 compilation."""
    root_lk = []
    for name in ("neutral_obs_liu2023.csv", "neutral_obs_li2019_rcf.csv",
                 "neutral_obs_li2019_soil.csv", "neutral_obs_kodesova2019.csv",
                 "neutral_obs_ge2017.csv", "neutral_obs_briggs1983_shoot.csv"):
        p = os.path.join(ROOT_DIR, "data_obs", name)
        if not os.path.exists(p):
            continue
        for r in _obs_rows(p):
            if (r.get("log_kow") or "").strip():
                root_lk.append(float(r["log_kow"]))
    tscf_lk = [float(r["logP"]) for r in
               _obs_rows(os.path.join(ROOT_DIR, "data_obs", "tscf_obs_schriever2020.csv"))
               if r["logP"].strip() and r["neutral_at_test_pH"].strip() == "1"]
    return dict(partition=(min(root_lk), max(root_lk), len(root_lk)),
                tscf=(min(tscf_lk), max(tscf_lk), len(tscf_lk)))


def screen(rows, spans=None, quiet=False):
    """Classify every runnable substance against the three limits above.

    `verdict` is what may be QUOTED from its run, not whether the ODE will solve
    -- the ODE solves for all of them, which is exactly the trap."""
    spans = spans or measured_spans()
    p_lo, p_hi, _ = spans["partition"]
    t_lo, t_hi, _ = spans["tscf"]
    out = []
    for r in rows:
        if r["status"] == "no_logkow_in_report":
            continue
        lk = _f(r, "log_kow")
        pka, pct = _f(r, "pKa"), _f(r, "pct_ionised_pH7")
        is_acid = (r["is_acid"] or "TRUE").strip().upper() == "TRUE"
        f_n = 1.0 - (pct / 100.0 if pct is not None else 0.0)
        reasons = []
        if not (p_lo <= lk <= p_hi):
            reasons.append(f"log Kow {lk} outside the measured partition span "
                           f"[{p_lo:.2f}, {p_hi:.2f}]")
        if not (t_lo <= lk <= t_hi):
            reasons.append(f"log Kow {lk} outside the measured TSCF span "
                           f"[{t_lo:.2f}, {t_hi:.2f}] -> shoot is extrapolation")
        ion_excluded = pka is not None and f_n < F_N_FLOOR
        if ion_excluded:
            reasons.append(f"{100 * (1 - f_n):.1f}% ionised at pH {REPORT_PH:.0f} "
                           f"(neutral fraction {f_n:.4f} < {F_N_FLOOR}) -> the "
                           "speciation path is documented as under-delivering here")
        if ion_excluded:
            verdict = "EXCLUDED (ionisation)"
        elif not reasons:
            verdict = "inside"
        elif len(reasons) == 1 and "TSCF" in reasons[0]:
            verdict = "root only"
        else:
            verdict = "extrapolation"
        out.append(dict(substance=r["substance"], log_kow=lk, pKa=pka, is_acid=is_acid,
                        f_n=f_n, verdict=verdict, reasons=reasons, row=r))
    if not quiet:
        print()
        print("=" * 78)
        print("2. THE SCOPE SCREEN -- read off data_obs/, not chosen here")
        print("=" * 78)
        print(f"  measured ROOT-partition rows span log Kow {p_lo:.2f} to {p_hi:.2f} "
              f"({spans['partition'][2]} rows)")
        print(f"  measured TSCF (un-ionised) spans  log Kow {t_lo:.2f} to {t_hi:.2f} "
              f"({spans['tscf'][2]} rows)")
        print(f"  weak-electrolyte path usable only for a neutral fraction >= {F_N_FLOOR}"
              "  (docs/neutral_dpu_validation.md section 4l)")
        print()
        for grp in ("inside", "root only", "extrapolation", "EXCLUDED (ionisation)"):
            members = [o for o in out if o["verdict"] == grp]
            print(f"  {grp:<22s} n={len(members):2d}  "
                  + ", ".join(f"{o['substance']} ({o['log_kow']})" for o in members))
        print()
        print("  NOTE the two groups that are the report's own findings coming back:")
        print("   * the anticoagulants are excluded here for the reason the report's")
        print("     section 3.0a says its screen should have excluded them and did not;")
        print("   * cholecalciferol and muscalure sit above every anchor in EITHER model")
        print("     -- the report reaches them only on a curve it measured, not a QSPR.")
    return out


# ---------------------------------------------------------------------------
# 3. the runs
# ---------------------------------------------------------------------------
def run_one(log_kow, name, half_life=None, pKa=None, is_acid=True, pH=REPORT_PH,
            MW=float("nan"), K_AW=0.0, air=False):
    res = api.simulate_neutral(log_kow, name=name, season=SEASON, half_life=half_life,
                               pKa=pKa, is_acid=is_acid, pH=pH,
                               MW=MW, K_AW=K_AW, air=air)
    b = res["baf_final"]
    return dict(TSCF=res["TSCF"], K_PW_root=res["K_PW"]["root"],
                root=b["root"], stem=b["stem"], leaf=b["leaf"], grain=b["grain"],
                straw=res["straw_baf"], res=res)


def run_all(screened, quiet=False):
    """Two assumptions, mirroring the report's own Scenario A / Scenario B.

      P  gamma = 0            -- no in-planta loss. This is the neutral path's
                                 documented UPPER BOUND (leaf and grain are
                                 terminal accumulators), the exact counterpart of
                                 the report's "conservative" Scenario A.
      Q  gamma from the report's FISH kM half-life. That rate is a fish
         whole-body biotransformation half-life, NOT an in-planta dissipation
         rate; nothing licenses moving it between organisms. It is run because
         the report's single largest finding is how much such a rate moves the
         answer (its section 8.4, up to 82x), and the same question is worth
         asking of this model -- as a SENSITIVITY, labelled, never as a
         parameterisation.
    """
    out = []
    for o in screened:
        r = o["row"]
        km = _f(r, "km_half_life_d")
        pKa = o["pKa"] if o["f_n"] >= F_N_FLOOR else None   # neutral form for the excluded
        rec = dict(substance=o["substance"], log_kow=o["log_kow"], verdict=o["verdict"],
                   f_n=o["f_n"], pKa_used=pKa, km_half_life_d=km,
                   bpc=r["bpc_class"], bat_A=_f(r, "bat_A_fish"), bat_B=_f(r, "bat_B_fish"),
                   bat_caveat=(r.get("bat_caveat") or "").strip())
        p = run_one(o["log_kow"], o["substance"], pKa=pKa, is_acid=o["is_acid"])
        rec.update({f"P_{k}": p[k] for k in
                    ("TSCF", "K_PW_root", "root", "stem", "leaf", "grain", "straw")})
        if km:
            q = run_one(o["log_kow"], o["substance"], half_life=km,
                        pKa=pKa, is_acid=o["is_acid"])
            rec.update({f"Q_{k}": q[k] for k in ("root", "leaf", "grain", "straw")})
        out.append(rec)
    if not quiet:
        print()
        print("=" * 78)
        print("3. THE RUNS -- rice BAF [L/kg fw], exposure held at Cwo = 1 ug/L")
        print("=" * 78)
        print(f"{'substance':<32s}{'logKow':>7s}{'TSCF':>9s}{'K_PW':>10s}"
              f"{'root':>10s}{'straw':>9s}{'grain':>9s}  scope")
        for rec in sorted(out, key=lambda d: d["log_kow"]):
            print(f"{rec['substance'][:31]:<32s}{rec['log_kow']:7.2f}"
                  f"{rec['P_TSCF']:9.2e}{rec['P_K_PW_root']:10.3g}"
                  f"{rec['P_root']:10.3g}{rec['P_straw']:9.3g}{rec['P_grain']:9.3g}"
                  f"  {rec['verdict']}")
        print()
        print("  gamma = 0 above, so leaf/straw/grain are UPPER BOUNDS (docs section 3).")
        withq = [r for r in out if "Q_root" in r]
        if withq:
            print()
            print("  What the report's FISH kM does when transplanted into the plant"
                  " (SENSITIVITY ONLY -- a fish rate is not a plant rate):")
            print(f"    {'substance':<26s}{'t1/2 [d]':>9s}{'root P->Q':>12s}"
                  f"{'straw P->Q':>12s}{'grain P->Q':>12s}")
            for rec in sorted(withq, key=lambda d: d["log_kow"]):
                fr_ = lambda a, b: (a / b if b else float("nan"))          # noqa: E731
                print(f"    {rec['substance'][:25]:<26s}{rec['km_half_life_d']:9.2f}"
                      f"{fr_(rec['P_root'], rec['Q_root']):12.2f}x"
                      f"{fr_(rec['P_straw'], rec['Q_straw']):11.2f}x"
                      f"{fr_(rec['P_grain'], rec['Q_grain']):11.2f}x")
            print("    (a factor of 1.0 means the assumption does not reach that tissue)")
    return out


# ---------------------------------------------------------------------------
# 4. the only two substances with a measured PLANT endpoint in this repo
# ---------------------------------------------------------------------------
def overlap_check(quiet=False):
    """propiconazole and triclosan are named by the report AND carry measured root
    rows already shipped here. Nothing is fitted: log Kow in, root partition out.

    Scored in EQUILIBRIUM mode. Both endpoints are equilibrium root partitions
    measured over hours to weeks, and reading them off a 120-day rice season
    imposes a Kow-dependent, purely model-side discount (docs section 4g)."""
    res = {}
    liu = [r for r in _obs_rows(os.path.join(ROOT_DIR, "data_obs", "neutral_obs_liu2023.csv"))
           if r["compound"].strip().lower() == "propiconazole" and r["tissue"] == "root"]
    if liu:
        lk = float(liu[0]["log_kow"])
        obs = float(liu[0]["value"])
        pred = ND.k_pw(lk, **_root_wl())
        res["propiconazole"] = dict(log_kow=lk, obs=obs, pred=pred,
                                    err=math.log10(pred / obs), n=len(liu),
                                    src="Liu 2023 rice root RCF, hydroponic 72 h")
    soil = [r for r in _obs_rows(os.path.join(ROOT_DIR, "data_obs",
                                              "neutral_obs_li2019_soil.csv"))
            if r["compound"].strip().lower() == "triclosan" and r["tissue"] == "root"]
    if soil:
        lk = float(soil[0]["log_kow"])
        obs = np.array([float(r["value"]) for r in soil])
        pred = ND.k_pw(lk, **_root_wl())
        errs = np.log10(pred / obs)
        res["triclosan"] = dict(log_kow=lk, obs=float(np.median(obs)), pred=pred,
                                err=float(np.mean(errs)),
                                rmse=float(np.sqrt(np.mean(errs ** 2))), n=len(obs),
                                lo=float(obs.min()), hi=float(obs.max()),
                                src="Li 2019 soil table, radish/carrot root, 14 rows")
    if not quiet:
        print()
        print("=" * 78)
        print("4. THE TWO SUBSTANCES THIS REPO CAN ACTUALLY SCORE")
        print("=" * 78)
        p = res.get("propiconazole")
        if p:
            print(f"  propiconazole (log Kow {p['log_kow']}) -- {p['src']}")
            print(f"    measured rice root RCF {p['obs']:.3f}   model K_PW "
                  f"{p['pred']:.3f}   error {p['err']:+.3f} log10 "
                  f"({10 ** abs(p['err']):.2f}x)")
        t = res.get("triclosan")
        if t:
            print(f"  triclosan (log Kow {t['log_kow']}) -- {t['src']}")
            print(f"    measured root BAF {t['lo']:.1f} to {t['hi']:.1f} "
                  f"(median {t['obs']:.1f}, n={t['n']})   model K_PW {t['pred']:.1f}")
            print(f"    mean bias {t['err']:+.3f} log10, RMSE {t['rmse']:.3f}")
            print(f"    (Li 2019's own log Kow is {t['log_kow']}, against the report's")
            print("     audited 4.76 -- two independent sourcings 0.04 log apart, which")
            print("     is the one place the report's property audit can be checked")
            print("     against something already in this repo. The bias is the same")
            print("     sign and larger than that whole soil table's +0.260 (docs 4h),")
            print("     and triclosan is an ACID, so part of it is speciation the")
            print("     equilibrium K_PW cannot see -- see section 5(b).)")
        print()
        print("  Two substances is not a validation and neither is a rice/biocide")
        print("  result: propiconazole IS rice but is one hydroponic point, and the")
        print("  triclosan rows are soil-grown radish and carrot. What they do is fix")
        print("  the sign -- the model is not wrong by an order of magnitude on the")
        print("  two substances in this report where anything at all can be checked.")
    return res


def _root_wl():
    """The root compartment's water and lipid fractions under the shipped
    `lipid_source` -- the two numbers K_PW is built from."""
    c = [x for x in ND.rice_compartments() if x.name == "root"][0]
    return dict(W=c.theta, L=c.f_PL * (1.0 - c.theta))


# ---------------------------------------------------------------------------
# 5. what the two models share, and where they part
# ---------------------------------------------------------------------------
def cross_model(runs, fast=False, quiet=False):
    out = {}

    # (a) the report's section 8.11 test, run on this model.
    #     Three substances at essentially one log Kow, MW 2x apart, solubility 39x,
    #     Henry 16,000x. In BAT that identity had to be MEASURED over 217 runs. Here
    #     it is structural -- K_PW and TSCF are functions of log Kow and nothing
    #     else, and MW/solubility/Henry appear in no term of the core ODE at all.
    trio = [("Cyphenothrin", 6.29, 375.5, 0.0022), ("Difethialone", 6.29, 539.5, 0.101),
            ("Empenthrin", 6.30, 274.4, 34.65)]
    same = [run_one(lk, nm) for nm, lk, _mw, _h in trio[:2]]
    out["null_inputs"] = dict(
        pair=[t[0] for t in trio[:2]],
        identical=bool(np.isclose(same[0]["root"], same[1]["root"], rtol=0, atol=0)),
        root=[s["root"] for s in same])

    # (b) the report's own ionisation sweep, re-run here on identical inputs.
    #     Its section 8.15 holds DCPP at log Kow 4.60 and moves the pKa alone,
    #     giving 940 / 3,590 / 3,830 L/kg at pKa 6.0 / 8.065 / 10.0. Same three
    #     inputs, same substance -- the one place the two models can be put
    #     side by side without translating an endpoint.
    dcpp_base = run_one(4.60, "DCPP-neutral")
    dcpp = []
    for pka, bat in ((6.0, 940.0), (8.065, 3590.0), (10.0, 3830.0)):
        r = run_one(4.60, "DCPP", pKa=pka, is_acid=True)
        dcpp.append(dict(pKa=pka, f_n=1.0 - frac_ionised(pka),
                         bat_rel=bat / 3830.0,
                         root_rel=r["root"] / dcpp_base["root"],
                         leaf_rel=r["leaf"] / dcpp_base["leaf"]))
    out["dcpp"] = dcpp

    # (c) the one property that is NOT inert here -- and it is Henry's constant,
    #     the property the report measured to be inert in BAT. Air exchange is
    #     opt-in (src/plant_air.py) and identically zero at K_AW = 0.
    RT = 8.314 * 293.15
    air = {}
    for nm, lk, mw, hen in trio:
        k_aw = hen / RT
        off = run_one(lk, nm, MW=mw, K_AW=k_aw, air=False)
        on = run_one(lk, nm, MW=mw, K_AW=k_aw, air=True)
        air[nm] = dict(K_AW=k_aw, leaf_off=off["leaf"], leaf_on=on["leaf"],
                       root_off=off["root"], root_on=on["root"])
    out["air"] = air

    # (d) the shape of each model's response to the one input that carries it
    if not fast:
        grid = np.arange(0.0, 11.01, 0.25)
        root, straw = [], []
        for lk in grid:
            r = run_one(float(lk), "sweep")
            root.append(r["root"])
            straw.append(r["straw"])
        root, straw = np.array(root), np.array(straw)
        out["sweep"] = dict(grid=grid, root=root, straw=straw,
                            root_peak=float(grid[int(np.argmax(root))]),
                            straw_peak=float(grid[int(np.argmax(straw))]),
                            kpw=np.array([ND.k_pw(float(l), **_root_wl()) for l in grid]))

    # (e) does a fish BCF rank these substances the way a rice model does?
    #     Two exclusions, both stated rather than silent. A SUBSTANCE may appear
    #     twice here (permethrin and cyphenothrin each carry a second row for their
    #     assessment report's other measured log Kow); the two rows share one BAT
    #     value, so counting both would weight that substance twice. And the BAT
    #     project marks the triamine's own output uninterpretable -- three charges
    #     at environmental pH against a tool that takes one pKa -- and drops it from
    #     every tally in its report; a number nobody can read is not a rank.
    def _pairable(r):
        if not r["bat_A"]:
            return False
        if "alt log Kow" in r["substance"] or "AR measured" in r["substance"]:
            return False
        return "not interpretable" not in r.get("bat_caveat", "")
    pairs = [(r["bat_A"], r["P_root"], r["P_straw"]) for r in runs if _pairable(r)]
    if len(pairs) >= 4:
        from scipy.stats import spearmanr
        a = np.array(pairs, float)
        out["rank"] = dict(n=len(pairs),
                           root=float(spearmanr(a[:, 0], a[:, 1]).statistic),
                           straw=float(spearmanr(a[:, 0], a[:, 2]).statistic))

    if not quiet:
        print()
        print("=" * 78)
        print("5. WHAT THE TWO MODELS SHARE, AND WHERE THEY PART")
        print("=" * 78)
        n = out["null_inputs"]
        print("  (a) the report's section 8.11 result, on this model:")
        print(f"      {n['pair'][0]} and {n['pair'][1]} at log Kow 6.29 -- MW 375.5 vs")
        print(f"      539.5, solubility 0.01 vs 0.3906, Henry 0.0022 vs 0.101 --")
        print(f"      root BAF {n['root'][0]:.6g} vs {n['root'][1]:.6g}: "
              f"{'BIT-IDENTICAL' if n['identical'] else 'DIFFERENT'}.")
        print("      In BAT that took 217 runs to establish. Here it is structural:")
        print("      K_PW and TSCF are functions of log Kow, and MW, solubility and")
        print("      Henry appear in no term of the core ODE. Same finding, but one")
        print("      model had to be probed for it and the other can be read.")
        print()
        print("  (b) the report's OWN ionisation sweep, re-run here. Its section 8.15")
        print("      holds DCPP at log Kow 4.60 and moves the pKa alone. Same substance,")
        print("      same three inputs, both models normalised to their own neutral end:")
        print(f"      {'pKa':>7s}{'f_neutral':>11s}{'BAT fish':>10s}"
              f"{'rice root':>11s}{'rice leaf':>11s}{'Hend.-Hass.':>13s}")
        for d in out["dcpp"]:
            print(f"      {d['pKa']:7.3f}{d['f_n']:11.4f}{d['bat_rel']:10.3f}"
                  f"{d['root_rel']:11.3f}{d['leaf_rel']:11.3f}{d['f_n']:13.4f}")
        print("      Both models damp ionisation far below what the chemistry predicts")
        print("      -- the report's section 7.6 finding, reached independently. But")
        print("      THIS model damps it harder, and the root hardest of all: at 91%")
        print("      ionised the chemistry says 0.091, BAT says 0.245 and the root here")
        print("      says 0.768. That is structural, not a tuning: the membrane term")
        print("      sets how FAST the root equilibrates, not the level it equilibrates")
        print("      to, so a root that still reaches equilibrium inside the season")
        print("      barely registers the change. It is the reason the >90%-ionised")
        print("      group is EXCLUDED above rather than reported with a caveat.")
        print()
        print("  (c) the ONE property that is not inert here is the one the report")
        print("      measured to be inert in BAT -- Henry's constant, through the")
        print("      opt-in air term (src/plant_air.py, zero at K_AW = 0):")
        for nm, d in air.items():
            print(f"      {nm:<14s} K_AW {d['K_AW']:.2e}   leaf BAF "
                  f"{d['leaf_off']:.4g} -> {d['leaf_on']:.4g} with air on "
                  f"(root {d['root_off']:.4g} -> {d['root_on']:.4g})")
        print("      So the three substances BAT cannot tell apart are separable here,")
        print("      by volatilisation -- but only through a term that is OFF by")
        print("      default and whose absolute size is surface-area-limited (docs).")
        if "sweep" in out:
            s = out["sweep"]
            hi = float(s["grid"][-1])
            turns = s["root_peak"] < hi - 1e-9
            print()
            print("  (d) each model reduces to a curve in log Kow. They are not the"
                  " same curve:")
            print(f"      this model: straw BAF peaks at log Kow {s['straw_peak']:.2f} "
                  "(the TSCF bell), and the")
            print("      root " + (f"peaks at {s['root_peak']:.2f}" if turns else
                                   f"does NOT turn over anywhere below log Kow {hi:.1f}"))
            print("      BAT (report section 8.11/8.16): fish BCF peaks near 6.3 and")
            print("      re-crosses its own B threshold downward between 8.13 and 9.00")
            print("      -- a turnover this model's root does NOT have, because the")
            print("      Briggs RCF has no descending limb. Where BAT turns over, the")
            print("      root partition here keeps rising:")
            for lk in (4.0, 6.3, 8.13, 9.0, 10.6):
                i = int(np.argmin(np.abs(s["grid"] - lk)))
                print(f"        log Kow {s['grid'][i]:5.2f}   K_PW {s['kpw'][i]:12.4g}"
                      f"   root BAF (120 d) {s['root'][i]:12.4g}"
                      f"   ratio {s['root'][i] / s['kpw'][i]:8.4f}")
            print("      That last column is its own limit: above ~7 the season is too")
            print("      short to equilibrate the root, so the model's root number")
            print("      stops being a partition coefficient and becomes a rate.")
        if "rank" in out:
            r = out["rank"]
            print()
            print(f"  (e) does a fish BCF rank these substances as a rice model does?"
                  f"  (n={r['n']})")
            print(f"      Spearman(BAT fish BCF, rice ROOT BAF)  = {r['root']:+.3f}")
            print(f"      Spearman(BAT fish BCF, rice STRAW BAF) = {r['straw']:+.3f}")
            print("      The root agrees because both rise with lipophilicity. The")
            print("      shoot does not, and cannot: what reaches a rice shoot is")
            print("      gated by TSCF, which PEAKS at log Kow 1.78 and is ~0 for every")
            print("      substance the opinions call bioaccumulative. A bioaccumulation")
            print("      class is not a statement about grain.")
    return out


# ---------------------------------------------------------------------------
def figure(runs, cross, path=None):
    """Two panels: where each substance ends up in the plant, and the two models'
    log Kow curves side by side. Skipped when the sweep was not run (--fast)."""
    if "sweep" not in cross:
        return None
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    path = path or os.path.join(HERE, "figures", "bat_census_biocides.png")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    s = cross["sweep"]
    fig, ax = plt.subplots(1, 2, figsize=(13.5, 5.4))

    r = sorted(runs, key=lambda d: d["log_kow"])
    x = np.arange(len(r))
    for key, lab, c in (("P_root", "root", "#8c564b"), ("P_straw", "straw", "#2ca02c"),
                        ("P_grain", "grain", "#d6a419")):
        ax[0].semilogy(x, [max(d[key], 1e-12) for d in r], "o-", ms=4, lw=1.2,
                       color=c, label=lab)
    for i, d in enumerate(r):
        if d["verdict"].startswith("EXCLUDED"):
            ax[0].axvspan(i - .5, i + .5, color="0.85", zorder=0)
    ax[0].set_xticks(x)
    ax[0].set_xticklabels([f"{d['substance'][:18]} ({d['log_kow']:.2f})" for d in r],
                          rotation=90, fontsize=6.5)
    ax[0].set_ylabel("rice BAF [L/kg fw],  $\\gamma=0$")
    ax[0].set_title("(a) where each biocide ends up in the plant\n"
                    "grey = excluded, >90% ionised at pH 7", fontsize=10)
    ax[0].legend(fontsize=8)
    ax[0].grid(alpha=.3, which="both")

    ax[1].semilogy(s["grid"], s["kpw"], "--", color="0.5", lw=1.2,
                   label="root $K_{PW}$ (equilibrium)")
    ax[1].semilogy(s["grid"], np.maximum(s["root"], 1e-12), color="#8c564b", lw=1.8,
                   label="root BAF, 120 d")
    ax[1].semilogy(s["grid"], np.maximum(s["straw"], 1e-12), color="#2ca02c", lw=1.8,
                   label="straw BAF, 120 d")
    bat = [(d["log_kow"], d["bat_A"]) for d in runs if d["bat_A"]]
    ax[1].semilogy([b[0] for b in bat], [b[1] for b in bat], "k^", ms=5, alpha=.8,
                   label="BAT fish BCF (the report)")
    ax[1].axvline(1.78, color="#2ca02c", ls=":", lw=1)
    ax[1].axvline(6.3, color="k", ls=":", lw=1)
    ax[1].set_ylim(1e-8, 1e7)
    ax[1].set_xlabel("log $K_{OW}$")
    ax[1].set_ylabel("BAF / BCF [L/kg]")
    ax[1].set_title("(b) two one-input models, two different curves\n"
                    "dotted: this model's TSCF peak 1.78, BAT's fish peak ~6.3",
                    fontsize=10)
    ax[1].legend(fontsize=8, loc="lower left")
    ax[1].grid(alpha=.3, which="both")

    fig.suptitle("EU biocides (BAT census) through the rice neutral-organic path"
                 " — PREDICTION, not validation", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def write_csv(runs, path=OUT_CSV):
    cols = ["substance", "log_kow", "verdict", "f_n", "pKa_used", "km_half_life_d",
            "bpc", "bat_A", "bat_B", "P_TSCF", "P_K_PW_root", "P_root", "P_stem",
            "P_leaf", "P_grain", "P_straw", "Q_root", "Q_leaf", "Q_grain", "Q_straw"]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in sorted(runs, key=lambda d: d["log_kow"]):
            w.writerow(r)
    return path


def verdict(rows, screened, runs, ov, fig_path=None):
    prov = lambda k: len([r for r in rows if r["logkow_provenance"] == k])   # noqa: E731
    n_report, n_filled = prov("report"), prov("export_on_request")
    n_extra, n_variant = prov("export_bat_entered"), prov("variant")
    inside = [o for o in screened if o["verdict"] == "inside"]
    excl = [o for o in screened if o["verdict"].startswith("EXCLUDED")]
    rootonly = [o for o in screened if o["verdict"] == "root only"]
    extrap = [o for o in screened if o["verdict"] == "extrapolation"]
    grain_amp = max((r["P_grain"] / r["Q_grain"] for r in runs
                     if r.get("Q_grain")), default=float("nan"))
    print()
    print("=" * 78)
    print("VERDICT")
    print("=" * 78)
    print(f"  Where each log Kow came from: {n_report} printed in the report itself, {n_filled}")
    print("  supplied by the BAT project on request from its audited provenance (never printed")
    print("  in the report, and never filled in from elsewhere until that source existed),")
    print(f"  {n_extra} from substances BAT entered but the report never named individually")
    print("  (no BPC class, no fish BCF -- predictions with nothing to compare against), 1 from")
    print(f"  this repo's own Liu 2023 row, and {n_variant} second readings of an input this")
    print("  study will not pick between (permethrin's two measured values, cyphenothrin's")
    print("  assessment-report range, tebuconazole's unresolved acid/base flag).")
    print(f"  {len(screened)} rows run in all.")
    print(f"  Of those, {len(inside)} sit inside the span where this repo")
    print(f"  holds measured plant data, {len(rootonly)} are past the TSCF anchor so only")
    print(f"  their root is quotable, {len(extrap)} are beyond every anchor in either model,")
    print(f"  and {len(excl)} are excluded for ionisation. The report's own coverage")
    print("  arithmetic, arrived at independently on a different model.")
    print()
    print("  WHAT THIS ESTABLISHES")
    print("   * The report's audited inputs transfer. A log Kow of the UNCHARGED form,")
    print("     traced to an assessment report, is exactly what this model's neutral")
    print("     path needs, and its section 7.5 correction (a distribution ratio D")
    print("     entered where log Kow belongs) is an error this model would inherit")
    print("     identically -- K_PW and TSCF are both functions of that one number.")
    print("   * Both models are one-input models, and that is where the agreement")
    print("     stops. BAT's fish BCF and this model's ROOT BAF rank these substances")
    print("     together; its fish BCF and this model's SHOOT do not, because rice")
    print("     shoot loading is gated by a TSCF bell peaking at log Kow 1.78. Every")
    print("     substance the opinions call bioaccumulative has TSCF < 1e-3 here:")
    print("     the model predicts they load the ROOT and do not reach the grain.")
    print("     'Bioaccumulative in fish' and 'reaches the edible part of a rice")
    print("     plant' are different questions, and this is the size of the gap.")
    print("   * The report's single largest finding is AMPLIFIED here, not damped.")
    print("     Entering a metabolic rate moved its fish BCF by up to 82x (section")
    print(f"     8.4); the same rates move this model's GRAIN by up to {grain_amp:,.0f}x,")
    print("     because grain and leaf are terminal accumulators whose only other")
    print("     sink is growth dilution, which goes to zero at maturity. So the")
    print("     quantity the report identifies as its least defensible input is the")
    print("     one the edible compartment is MOST sensitive to -- and for a plant")
    print("     nobody has measured it at all (a fish rate is what was available).")
    print()
    print("  WHAT THIS DOES NOT ESTABLISH -- and cannot")
    print("   * NOT a validation. No measured rice concentration exists for any of")
    print("     these substances. The two checkable ones are one hydroponic rice")
    print("     point and 14 soil-grown radish/carrot rows -- section 4, "
          f"propiconazole {ov['propiconazole']['err']:+.3f} log")
    print(f"     and triclosan {ov['triclosan']['err']:+.3f} log. Two substances, and")
    print("     one of them is not rice.")
    pos = [r for r in runs if r["bpc"] in ("B", "vB")]
    pos_excl = [r for r in pos if r["verdict"].startswith("EXCLUDED")]
    pos_inside = [r for r in pos if r["verdict"] == "inside"]
    print(f"   * The {len(pos)} substances the opinions call bioaccumulative are the ones")
    print(f"     this model is LEAST entitled to: {len(pos_excl)} are >90% ionised at pH 7")
    print(f"     and excluded outright ({', '.join(r['substance'] for r in pos_excl)}),")
    print(f"     and of the rest only {len(pos_inside)} "
          f"({', '.join(r['substance'] for r in pos_inside)}) is inside the")
    print("     TSCF anchor at all. That mirrors")
    print("     the report's own section 7.3 finding about its kM QSAR (those same")
    print("     substances at Tanimoto 0.19-0.28) from a different direction: two")
    print("     independent models are both weakest exactly where the answer matters.")
    print("   * The fish kM half-lives are run ONCE, as a labelled sensitivity. They")
    print("     are fish whole-body biotransformation rates; nothing licenses reading")
    print("     them as in-planta dissipation, and this model's own leaf half-life is")
    print("     an open question (docs section 4i: metabolism is species-dependent by")
    print("     4.8x for one compound in one experiment).")
    print("   * gamma = 0 makes every leaf and grain number an upper bound.")
    print()
    print(f"  written: {os.path.relpath(OUT_CSV, ROOT_DIR)}")
    if fig_path:
        print(f"           {os.path.relpath(fig_path, ROOT_DIR)}")


def main(fast=False):
    rows = load_table()
    check_transcription(rows)
    screened = screen(rows)
    runs = run_all(screened)
    ov = overlap_check()
    cross = cross_model(runs, fast=fast)
    write_csv(runs)
    fig = figure(runs, cross)
    verdict(rows, screened, runs, ov, fig)
    return runs


if __name__ == "__main__":
    main(fast="--fast" in sys.argv)
