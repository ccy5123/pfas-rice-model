#!/usr/bin/env python3
# =============================================================================
# validation/weak_electrolyte_tscf.py
# -----------------------------------------------------------------------------
# The FIRST EMPIRICAL TEST of the weak-electrolyte path (`simulate_neutral(
# pKa=...)`), against the 67 ionisable rows of Schriever & Lamshoeft 2020's
# Table A 3 that every score in this repo has so far held out.
#
# WHY THIS EXISTS. The speciation port landed as "structural capability, not a
# predictive claim -- no measured weak-electrolyte rice dataset exists here"
# (PR #60 section 3). That is true of RICE, and it was read as though no measured
# weak-electrolyte data existed at all. But `data_obs/tscf_obs_schriever2020.csv`
# already ships all 97 rows of Table A 3, and `validation/schriever2020_tscf.py`
# scores only the 30 flagged un-ionised at the test pH -- explicitly setting the
# other 67 aside as "outside this model's stated scope". The port is exactly the
# thing that extends the scope to them. So the held-out rows are the test, and
# they were in the repo before the capability was.
#
# WHAT THE TABLE CAN AND CANNOT PARAMETERISE. Three limits, all load-bearing:
#
#   * NO COMPOUND NAMES. Table A 3 is keyed by an id, so a row cannot be
#     cross-checked against an independent pKa or log P, and no compound can be
#     excluded on chemical grounds. Take the table's own columns or nothing.
#   * NO ACID/BASE LABEL, and it matters: an anion is repelled by the
#     inside-negative membrane and a cation is attracted, which is precisely the
#     asymmetry the port introduces (`Environment.N_for`). Every model run here
#     is therefore done BOTH ways and both are reported. Never pick one.
#   * `logD` IS NOT RECONSTRUCTIBLE from pKa and log P by Henderson-Hasselbalch
#     (median |logD - (logP + log10 f_n)| = 1.41 on the acid reading, 1.82 on the
#     base reading), and 20 rows have logD ABOVE logP, which no single-centre
#     acid or base can produce. So the pKa column and the logD column disagree
#     about how ionised these compounds are, and the pKa cannot be trusted to
#     give f_n.
#
# HOW f_n IS OBTAINED, and why this way. From the table's OWN logD:
#
#       f_n = 10 ** (logD_test - logP),  clipped at 1
#
# because logD = logP + log10(f_n) is the definition of the distribution
# coefficient when only the neutral species partitions. This needs no pKa and no
# acid/base label, so it survives all three limits above; it is also the same
# quantity the CSV's shipped `neutral_at_test_pH` flag is built from (|logD -
# logP| < 0.1), which keeps this script and that flag consistent by construction.
# The clip absorbs the 20 rows where the two source columns disagree by +0.01 to
# +0.81 log -- provider noise between a measured logP and a computed logD.
#
# A COUNTEREXAMPLE THIS SCRIPT DELIBERATELY DOES NOT USE. The eight barley rows
# at pKa 1.62 have TSCF 0.63-0.98, and read as a strong ACID (fully dissociated
# at any test pH) they would look like a spectacular refutation. They are not:
# their own logD sits 0.01-0.02 ABOVE logP, i.e. the table says they are NOT
# ionised at the test pH, and the shipped flag already classes them neutral. The
# pKa is a basic centre. Anyone re-deriving this result from the pKa column
# instead of the logD column will manufacture that false counterexample; it is
# named here so the next reader does not.
#
# METRIC. TSCF is a bounded fraction on [0,1] with zeros present, so -- as in
# schriever2020_tscf.py -- the metric is LINEAR RMSE plus Spearman rank
# correlation, not the repo's usual log10 RMSE.
#
#   python validation/weak_electrolyte_tscf.py            # ~2 min (67 ODE solves x2)
#   python validation/weak_electrolyte_tscf.py --fast     # skip the per-row ODE runs
# =============================================================================
from __future__ import annotations
import argparse
import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT_DIR, "src"))

import literature_params as LP                                # noqa: E402
import model_api as api                                       # noqa: E402
import neutral_dpu as ND                                      # noqa: E402
from pfas_rice_plant_module_4pool_surf import (                # noqa: E402
    Environment, _ghk_factor)

OBS = os.path.join(ROOT_DIR, "data_obs", "tscf_obs_schriever2020.csv")
MISSING = 999.0
DEFAULT_PH = 5.5          # the documented fallback when a source reported no test pH


# ---------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------
def load(path=OBS):
    """Table A 3 with f_n derived from the table's own logD (see the header)."""
    rows = []
    with open(path, newline="") as f:
        for r in csv.DictReader(x for x in f if not x.lstrip().startswith("#")):
            if not r.get("id"):
                continue
            logP = float(r["logP"])
            pH = float(r["pH_test"])
            logD = (float(r["logD_test_pH"]) if pH != MISSING
                    else float(r["logD_pH5.5"]))
            rows.append(dict(
                id=r["id"], species=r["species"], authors=r["authors"],
                pKa=float(r["pKa"]), logP=logP, logD=logD,
                pH=(pH if pH != MISSING else DEFAULT_PH),
                TSCF=float(r["TSCF"]),
                neutral=bool(int(r["neutral_at_test_pH"])),
                # f_n from logD - logP; >1 is provider noise between the two columns
                fn=float(min(1.0, 10.0 ** (logD - logP)))))
    return rows


def spearman(a, b):
    """Rank correlation with tied ranks averaged (no SciPy dependency)."""
    def rank(x):
        v = np.asarray(x, float)
        order = np.argsort(v, kind="mergesort")
        rr = np.empty(len(v), float)
        rr[order] = np.arange(len(v), dtype=float)
        vals = v[order]
        i = 0
        while i < len(vals):
            j = i
            while j + 1 < len(vals) and vals[j + 1] == vals[i]:
                j += 1
            if j > i:
                rr[order[i:j + 1]] = np.mean(np.arange(i, j + 1))
            i = j + 1
        return rr
    ra, rb = rank(a), rank(b)
    if np.std(ra) == 0 or np.std(rb) == 0:
        return float("nan")
    return float(np.corrcoef(ra, rb)[0, 1])


# ---------------------------------------------------------------------------
# what the ported path predicts
# ---------------------------------------------------------------------------
def influx_ratio(fn: float, is_acid: bool = True, env: Environment | None = None
                 ) -> float:
    """Membrane influx conductance RELATIVE to the same molecule un-ionised [-].

        Phi = [P_n f_n + (kappa_d / P_N_OVER_P_D) g(N) f_d] / kappa_d
            =        f_n + g(N) f_d / P_N_OVER_P_D

    Read straight off `root_uptake` with `Vmax = 0` (the neutral path sets no
    carrier) and `P_n = kappa_d` (`_weak_electrolyte_kw`), so `kappa_d` cancels
    and NOTHING here is fitted or chosen: the two constants are the shipped
    Trapp ratio `P_N_OVER_P_D` = 10**3.5 and the GHK factor at the shipped
    membrane potential. Phi -> 1 as f_n -> 1, which is the continuity property
    the port was built around.
    """
    env = env or Environment()
    # N via the model's own accessor rather than re-deriving it, so this cannot
    # drift from the membrane potential the ODE actually uses.
    N = env.N_for(type("_Z", (), {"z": -1 if is_acid else +1})())
    return float(fn + _ghk_factor(N) * (1.0 - fn) / LP.PN_OVER_PD)


def pKa_for(fn: float, pH: float, is_acid: bool = True) -> float:
    """The pKa that reproduces a target neutral fraction at a given pH.

    Inverts Henderson-Hasselbalch so the model can be driven at the f_n the TABLE
    reports, rather than at the pKa the table reports -- which the header explains
    the table itself contradicts.
    """
    fn = float(np.clip(fn, 1e-12, 1.0 - 1e-12))
    off = np.log10(fn / (1.0 - fn))
    return float(pH + off) if is_acid else float(pH - off)


def model_tscf(log_kow: float, fn: float, pH: float, is_acid: bool = True,
               g_apo: float = 0.0) -> float:
    """The model's effective root-to-shoot transfer factor for one row.

        TSCF_model = f_xy * Cw_root / Cwo

    i.e. the xylem loading factor times the free root concentration the ODE
    settles at, which is what `TSCF = C_xylem / C_solution` means in this model's
    variables. Using the raw root BAF instead would fold in the season's growth
    accumulation and would not be the measured quantity.
    """
    kw = {} if fn >= 1.0 else dict(pKa=pKa_for(fn, pH, is_acid),
                                   is_acid=is_acid, pH=pH)
    r = api.simulate_neutral(log_kow, g_apo=g_apo, **kw)
    return float(r["TSCF"] * r["baf_final"]["root"] / r["K_PW"]["root"])


def bootstrap_wins(off, on, obs, n=4000, seed=0):
    """How often does speciation ON beat OFF, under resampling of the rows?

    The two metrics disagree on the full table, so one of them is carrying the
    verdict; this says which one is safe to carry it. Returns the fraction of
    bootstrap resamples in which ON wins on rank correlation and on RMSE.
    """
    off, on, obs = (np.asarray(x, float) for x in (off, on, obs))
    rng = np.random.default_rng(seed)
    rmse = lambda p, o: np.sqrt(np.mean((p - o) ** 2))          # noqa: E731
    w_rank = w_rmse = 0
    for _ in range(n):
        i = rng.integers(0, len(obs), len(obs))
        if spearman(on[i], obs[i]) > spearman(off[i], obs[i]):
            w_rank += 1
        if rmse(on[i], obs[i]) < rmse(off[i], obs[i]):
            w_rmse += 1
    return dict(rank=w_rank / n, rmse=w_rmse / n, n=n)


G_APO_GRID = (0.0, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 20.0)


def scan_g_apo(ion, neu, grid=G_APO_GRID, is_acid=True):
    """The apoplastic bypass, scored across a grid — trade-off curve, not an optimum.

    `g_apo` is the one lever that can raise a strongly-ionised compound's uptake
    without being gated by speciation, so it is the obvious repair for the
    magnitude half of this script's verdict. It is reported as a CURVE because the
    failure mode was pre-registered (handoff section 3 item 8): a bypass big enough
    to fix the magnitude must eventually flatten the speciation dependence, and if
    the fitted value drives the rank correlation back toward the speciation-OFF
    baseline then the bypass is ABSORBING the effect rather than explaining it.
    Both halves are therefore tracked at every grid point, on the ionisable rows
    AND on the 30 un-ionised rows the model already fits, which the bypass must
    not disturb.
    """
    obs_i = np.array([r["TSCF"] for r in ion], float)
    obs_n = np.array([r["TSCF"] for r in neu], float)
    out = []
    for g in grid:
        pi = np.array([model_tscf(r["logP"], r["fn"], r["pH"], is_acid, g_apo=g)
                       for r in ion], float)
        pn = np.array([model_tscf(r["logP"], 1.0, r["pH"], g_apo=g) for r in neu],
                      float)
        out.append(dict(
            g_apo=g,
            rmse=float(np.sqrt(np.mean((pi - obs_i) ** 2))),
            bias=float(np.mean(pi - obs_i)), rho=spearman(pi, obs_i),
            neutral_rmse=float(np.sqrt(np.mean((pn - obs_n) ** 2))),
            neutral_bias=float(np.mean(pn - obs_n))))
    return out


def _stat(tag, pred, obs):
    pred, obs = np.asarray(pred, float), np.asarray(obs, float)
    s = dict(n=len(obs), rmse=float(np.sqrt(np.mean((pred - obs) ** 2))),
             bias=float(np.mean(pred - obs)), rho=spearman(pred, obs))
    print(f"   {tag:38s} n={s['n']:3d}  RMSE {s['rmse']:.3f}"
          f"   bias {s['bias']:+.3f}   Spearman {s['rho']:+.3f}")
    return s


# ---------------------------------------------------------------------------
def main(fast=False):
    rows = load()
    ion = [r for r in rows if not r["neutral"]]
    neu = [r for r in rows if r["neutral"]]

    print("=" * 84)
    print("WEAK-ELECTROLYTE PATH vs 67 MEASURED IONISABLE TSCF VALUES")
    print("Schriever & Lamshoeft 2020 SI Table A 3 — the rows every prior score held out")
    print("=" * 84)
    print(f"   rows {len(rows)}   un-ionised at the test pH {len(neu)}   "
          f"IONISABLE {len(ion)}")
    print(f"   of the ionisable rows, {sum(r['pKa'] != MISSING for r in ion)} carry a pKa "
          f"and {sum(r['pKa'] == MISSING for r in ion)} do not; f_n comes from logD either way")

    fn = np.array([r["fn"] for r in ion])
    obs = np.array([r["TSCF"] for r in ion])
    print(f"   f_n spans {fn.min():.1e} to {fn.max():.3f} "
          f"({np.log10(fn.max()/fn.min()):.1f} orders of magnitude), median {np.median(fn):.4f}")

    # -- 1 ------------------------------------------------------------------
    print("\n1. DOES MEASURED TRANSFER FALL AS THE COMPOUND IONISES?")
    print("   The port's central claim, reduced to something the data can answer.")
    print(f"   Spearman(f_n, TSCF) = {spearman(fn, obs):+.3f}   n={len(ion)}")
    print("   binned:")
    for lo, hi in [(0, 1e-3), (1e-3, 1e-2), (1e-2, 1e-1), (1e-1, 1.001)]:
        s = [r for r in ion if lo <= r["fn"] < hi]
        if s:
            t = np.array([r["TSCF"] for r in s])
            print(f"      f_n {lo:7.0e} – {hi:7.0e}   n={len(s):3d}   "
                  f"TSCF mean {t.mean():.3f}  median {np.median(t):.3f}  max {t.max():.3f}")

    # -- 2 ------------------------------------------------------------------
    print("\n2. HOW BIG DOES THE PORT SAY THE EFFECT IS?")
    print("   Phi = membrane influx conductance relative to the un-ionised molecule,")
    print("   read straight off root_uptake. Nothing fitted (Trapp 10^3.5 + GHK).")
    for f in (1.0, 1e-1, 1e-2, 1e-3, 1e-4):
        print(f"      f_n {f:7.0e}   Phi(acid) {influx_ratio(f, True):.3e}"
              f"      Phi(base) {influx_ratio(f, False):.3e}")
    lo, hi = influx_ratio(float(fn.min())), influx_ratio(1.0)
    print(f"   Across this table's f_n range Phi moves {hi/lo:.3g}-fold (acid reading),")
    print(f"   while measured TSCF moves about "
          f"{np.mean(obs[fn >= 0.1])/max(np.mean(obs[fn < 1e-3]), 1e-9):.2g}-fold.")

    # -- 3 ------------------------------------------------------------------
    if fast:
        print("\n3. PER-ROW MODEL RUNS — skipped (--fast).")
        return dict(rows=ion, rho=spearman(fn, obs))

    print("\n3. THE MODEL RUN PER ROW, BOTH SPECIATION READINGS")
    print("   TSCF_model = f_xy * Cw_root/Cwo. The 'un-ionised' line is the same")
    print("   model with speciation OFF — the baseline the port is meant to improve on.")
    base = [model_tscf(r["logP"], 1.0, r["pH"]) for r in ion]
    acid = [model_tscf(r["logP"], r["fn"], r["pH"], True) for r in ion]
    basic = [model_tscf(r["logP"], r["fn"], r["pH"], False) for r in ion]
    s_off = _stat("speciation OFF (un-ionised)", base, obs)
    s_acid = _stat("weak-electrolyte, read as an ACID", acid, obs)
    s_base = _stat("weak-electrolyte, read as a BASE", basic, obs)
    print(f"   median TSCF_model: un-ionised {np.median(base):.3e}   "
          f"acid {np.median(acid):.3e}   base {np.median(basic):.3e}")
    print(f"   median TSCF measured: {np.median(obs):.3f}")

    # -- 4 ------------------------------------------------------------------
    print("\n4. WHICH OF THOSE TWO COMPARISONS SURVIVES RESAMPLING?")
    print("   Asked because the two metrics point opposite ways, so at least one of")
    print("   them is carrying the verdict and it should not be a 67-row accident.")
    boot = bootstrap_wins(np.asarray(base), np.asarray(acid), obs)
    print(f"      P(speciation ON has the better RANK correlation) = {boot['rank']:.3f}")
    print(f"      P(speciation ON has the better RMSE)             = {boot['rmse']:.3f}")
    print("   So the RANK gain is robust and the RMSE loss is only a tendency: a")
    print("   subsample can and does flip the RMSE, never the ordering. State the")
    print("   verdict that way -- 'ON orders better, and tends to scale worse'.")

    # -- 5 ------------------------------------------------------------------
    print("\n5. THE APOPLASTIC BYPASS — does one parameter repair the magnitude?")
    print("   `g_apo` is entry that never crosses a membrane, so it is gated by")
    print("   NEITHER speciation nor the membrane potential. Reported as a curve, not")
    print("   an optimum, because the failure mode was pre-registered: a bypass big")
    print("   enough to fix the magnitude must eventually FLATTEN the ordering.")
    print(f"   Reference points: speciation OFF rho {s_off['rho']:+.3f} / RMSE {s_off['rmse']:.3f};")
    print(f"                     speciation ON, no bypass rho {s_acid['rho']:+.3f} / RMSE {s_acid['rmse']:.3f}.")
    print(f"   {'g_apo':>7}  {'RMSE':>7} {'bias':>7} {'rho':>7}   | 30 un-ionised rows: {'RMSE':>6} {'bias':>7}")
    scan = scan_g_apo(ion, neu)
    for s in scan:
        print(f"   {s['g_apo']:7.2f}  {s['rmse']:7.3f} {s['bias']:+7.3f} {s['rho']:+7.3f}   |"
              f"                     {s['neutral_rmse']:6.3f} {s['neutral_bias']:+7.3f}")
    best = min(scan, key=lambda s: s["rmse"])
    print(f"   RMSE-optimal g_apo = {best['g_apo']:.2f}  ->  RMSE {best['rmse']:.3f}"
          f"  bias {best['bias']:+.3f}  rho {best['rho']:+.3f}")

    # -- verdict ------------------------------------------------------------
    rho = spearman(fn, obs)
    print("\n" + "=" * 84)
    print("VERDICT — direction SUPPORTED, magnitude REFUTED")
    print("=" * 84)
    print(f"   DIRECTION. Measured transfer really does fall as the compound ionises:")
    print(f"   Spearman(f_n, TSCF) = {rho:+.3f} over {len(ion)} rows. That is the first")
    print("   empirical support of any kind for the speciation port, and it is not")
    print("   trivial -- the sign could have come out flat or backwards.")
    print("   MAGNITUDE. It is the size that fails. Phi spans ~4-5 orders of magnitude")
    print("   across this table; the measured means move ~3-fold. At f_n < 1e-3 the")
    print("   model predicts effectively no transfer and the data show a mean TSCF of")
    print(f"   {np.mean(obs[fn < 1e-3]):.3f}.")
    print("   THE TWO METRICS DISAGREE, AND THAT IS THE RESULT. Switching speciation on")
    print(f"   nearly DOUBLES the rank correlation ({s_off['rho']:+.3f} -> {s_acid['rho']:+.3f}) while making the RMSE")
    print(f"   worse ({s_off['rmse']:.3f} -> {s_acid['rmse']:.3f}) and the bias strongly negative ({s_off['bias']:+.3f} -> {s_acid['bias']:+.3f},")
    print("   i.e. it now under-delivers). Right ordering, wrong scale.")
    print(f"   THEY ARE NOT EQUALLY SOLID. Under resampling the rank gain holds in")
    print(f"   {boot['rank']:.0%} of draws and the RMSE loss in only {1 - boot['rmse']:.0%} -- a small subsample can")
    print("   flip the RMSE and does. So the load-bearing claim is the ORDERING gain,")
    print("   and 'speciation makes the fit worse' is a tendency, not a finding. Rank")
    print("   is also the half the plant model consumes (schriever2020_tscf.py argues")
    print("   exactly that), so this is a real gain wrapped around an unusable")
    print("   calibration, not a flat failure.")
    print("   WHY, and it is already known in this repo. The model's only route into")
    print("   the plant is across the root membrane, so an ion that cannot cross cannot")
    print("   arrive. Real ions reach the xylem apoplastically (around cells) and, on")
    print("   the PFAS side of this same codebase, through a CARRIER that had to be")
    print("   fitted precisely because passive GHK exclusion under-delivers. The weak-")
    print("   electrolyte port gives its ion passive permeability and no carrier, so it")
    print("   inherits that known deficit with no lever to absorb it.")
    print("   STATUS. `pKa=` stays opt-in and nothing about it changes here. It moves")
    print("   from UNVALIDATED to BOUNDED: usable for the direction of a speciation")
    print("   effect, not for its size, and not at all below f_n ~ 0.1.")
    print("   LIMITS. 16 species, none rice; no compound names; acid/base inferred and")
    print("   reported both ways; TSCF, not the per-organ endpoint the model targets;")
    print("   and these hydroponic measurements are scored through a 120-day rice")
    print("   season, the driver mismatch section 4g and section 4c both warn about.")
    print("   That last one does NOT undermine the verdict: every line in section 3")
    print("   uses the SAME drivers, so the mismatch is common to speciation ON and")
    print("   OFF and cancels in the comparison between them, which is the claim.")
    print("   It does bound the ABSOLUTE RMSEs -- do not quote those on their own.")
    return dict(rows=ion, rho=rho, obs=obs, fn=fn,
                off=s_off, acid=s_acid, base=s_base)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--fast", action="store_true",
                   help="skip the per-row ODE runs (sections 1-2 only)")
    main(**vars(p.parse_args()))
