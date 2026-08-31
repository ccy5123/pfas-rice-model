#!/usr/bin/env python3
# =============================================================================
# validation/entry_vs_sequestration.py
# -----------------------------------------------------------------------------
# DOES THE CHAIN-LENGTH DEPENDENCE BELONG TO ENTRY, OR TO SEQUESTRATION?
# Handoff item B2 (docs/HANDOFF_carrier_vs_bypass.md section 3).
#
# WHY THIS IS THE QUESTION, and why it is worth a run rather than an assertion.
# Two independent fits of DIFFERENT root-entry terms both turned out to need the
# same chain-length correction:
#
#   * the apoplastic BYPASS (validation/carrier_vs_bypass.py section 3): fitted
#     per congener, g_apo trends at corr(n_C, log10 g_apo) = +0.832 over a 25x
#     spread -- which REFUTED theory_anchor.tex's claim that eta is "essentially
#     independent of tail length";
#   * the CARRIER (LC6, runs/pfas-rice-carrier): the per-congener Vmax multiplier
#     needed is ~flat to C10, then 2.0x at C11 and 5.5x at C12.
#
# A requirement common to two different entry terms is unlikely to be a property
# of entry. docs/theory_anchor.tex now records that reading -- that eq. (factor),
# f_xy = eta(E_m, carrier) x phi_free(n), is MIS-PARTITIONED, and what the fits
# absorb into eta belongs in phi_free / B_k (sequestration) -- but records it as a
# reading, with no replacement asserted. This file tests it.
#
# THE TEST. The repo already HAS a resolved sequestration term: the two-pool root
# (model_api.simulate_twopool_seq), whose non-K_PL U-shaped k_seq(n, head_group)
# was fitted to Yamazaki and separates PFOS from PFUnDA at identical K_PL. So:
# fit the ENTRY conductance per congener against a model that has sequestration
# resolved, and see whether the chain-length trend survives.
#
#   If the mis-partition reading is right, the trend should COLLAPSE -- the
#   sequestration term has already taken the chain-length signal, and entry has
#   nothing left to explain. If the trend SURVIVES, the dependence is genuinely
#   at the membrane and the theory doc's factorisation is wrong in a different
#   way than section B3 supposes.
#
# =============================================================================
# PRE-REGISTRATION -- written and committed BEFORE the run.
# =============================================================================
#
# SIX ARMS. Every arm uses the monotone physical f_xy (a-priori, NOT fitted),
# the carrier OFF and the bypass as the entry term, and the SAME demo forcings,
# so only the root model differs:
#
#   S    single pool                     (simulate, lipid loading off)
#   S+L  single pool + lipid loading      (simulate(lipid_loading=True))
#   T    two-pool, k_seq U-shaped + lipid (simulate_twopool_seq, cached fit)
#   C    two-pool with k_seq = 0 + lipid  (kseq_override=0 -- the CONTROL)
#
# each fitted twice: with ONE GLOBAL g_apo, and with g_apo PER CONGENER.
#
# WHY THE CONTROL EXISTS, and it is the thing most likely to be got wrong here:
# T differs from S by TWO changes, not one -- it adds sequestration AND
# lipid-facilitated loading. Attributing a collapse to sequestration without
# separating them would be exactly the error this file is meant to test for. The
# ladder S -> S+L -> T isolates each, and C (sequestration removed, lipid kept)
# is the direct control. **If C flattens the trend as much as T does, the credit
# belongs to LIPID LOADING, not sequestration, and section B3's reading needs
# revising rather than confirming** -- lipid loading is a LOADING term, not a
# binding term, so that outcome would move the chain-length dependence to a third
# place, neither entry nor B_k.
#
# FORCINGS, and why the published +0.832 is NOT the baseline here. The two-pool's
# cached fit is only valid on the DEMO forcings (documented in its docstring),
# while carrier_vs_bypass ran on the measured forcings with the ORYZA biomass. A
# comparison across that difference would be meaningless, so every arm here runs
# on the demo forcings and the single-pool baseline trend is RE-MEASURED under
# matched conditions. If the re-measured S trend does not itself come out clearly
# positive, this test has no baseline to speak from and the run is INCONCLUSIVE
# regardless of what the other arms do -- that is gate G0 below.
#
# DECISION RULES, fixed here:
#   G0. GATE. The S arm must reproduce a clear chain-length trend on these
#       forcings: corr(n_C, log10 g_apo) >= +0.5 AND spread >= 5x. Otherwise
#       INCONCLUSIVE, full stop -- there is nothing for sequestration to remove.
#   R1. SUPPORTED (the mis-partition reading is confirmed) iff BOTH:
#       (a) T's trend is at most HALF of S's on both measures (corr and log
#           spread), AND
#       (b) T's PENALTY for using one global g_apo instead of per-congener is at
#           most half of S's. This is the "removes the NEED" half: a flat fitted
#           parameter means little if a global one still fits badly.
#   R2. REFUTED if T's trend is essentially S's (corr within 20%, spread within
#       2x) -- the chain-length dependence is at entry after all.
#   R3. ATTRIBUTION. If C (no sequestration, lipid kept) flattens as much as T,
#       report the effect as LIPID LOADING rather than sequestration, whatever
#       R1 says. This overrides R1's wording, not its arithmetic.
#   Any close call goes through a bootstrap over congeners, per the standing rule
#   on this thread -- three times now a gap that looked real did not survive one.
#
# WHAT THIS CANNOT DO, stated so it is not claimed later. The two-pool's k_seq
# was FITTED on the same Yamazaki root data that scores it here, so T starts with
# an in-sample advantage that S does not have; the RMSE LEVELS are therefore not
# comparable between arms and are reported only for context. What IS comparable
# is the SHAPE of the fitted entry parameter and each arm's OWN global-vs-
# per-congener penalty, both of which are internal to an arm. The lipid
# conductances also differ between the single-pool fit and the two-pool's cached
# ones, so S+L and T do not share a lipid parameterisation. And nothing here can
# promote k_seq: its gate is the section 5 wet-lab assay, unchanged.
#
#   python validation/entry_vs_sequestration.py          # ~4 min
#   python validation/entry_vs_sequestration.py --fast   # coarser grid, ~2 min
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

import model_api as api                                          # noqa: E402

OBS = os.path.join(ROOT_DIR, "data_obs", "obs_baf_Yamazaki.csv")
KEYS = ("root", "straw", "grain")
SEASON = 120.0
N_T = 241
# demo forcings for every arm -- see the FORCINGS note in the header
FORCING = dict(season=SEASON, measured_forcing=False, biomass="growth_rice")
GATE_CORR = 0.5
GATE_SPREAD = 5.0


def load_obs(path=OBS):
    obs = {}
    with open(path) as f:
        for r in csv.DictReader(f):
            obs.setdefault(r["compound"], {})[r["tissue"]] = float(r["baf"])
    return {k: v for k, v in obs.items() if all(t in v for t in KEYS)}


# ---------------------------------------------------------------------------
# the four root models, behind one signature
# ---------------------------------------------------------------------------
def predict(arm, name, g_apo):
    """(root, straw, grain) BAF for one congener under one arm, at one g_apo.

    Carrier OFF in every arm (vmax_scale=0) so the ENTRY term being fitted is the
    bypass alone and the arms differ only in the ROOT model.
    """
    kw = dict(vmax_scale=0.0, g_apo=g_apo, n_t=N_T)
    if arm == "S":
        r = api.simulate(name, f_xy_source="recommended", **FORCING, **kw)
    elif arm == "S+L":
        r = api.simulate(name, f_xy_source="recommended", lipid_loading=True,
                         **FORCING, **kw)
    elif arm == "T":
        r = api.simulate_twopool_seq(name, **FORCING, **kw)
    elif arm == "C":
        r = api.simulate_twopool_seq(name, kseq_override=0.0, **FORCING, **kw)
    else:
        raise ValueError(arm)
    return (r["baf_final"]["root"], r["straw_baf"], r["baf_final"]["grain"])


def residuals(arm, names, obs, g_apo):
    out = []
    for nm in names:
        gv = g_apo[nm] if isinstance(g_apo, dict) else g_apo
        for t, p in zip(KEYS, predict(arm, nm, gv)):
            out.append(np.log10(max(p, 1e-12)) - np.log10(obs[nm][t]))
    return np.asarray(out, float)


def rmse(arm, names, obs, g_apo):
    return float(np.sqrt(np.mean(residuals(arm, names, obs, g_apo) ** 2)))


def fit_global(arm, names, obs, grid):
    return min(((v, rmse(arm, names, obs, v)) for v in grid), key=lambda vr: vr[1])


def fit_per_congener(arm, names, obs, grid):
    best = {}
    for nm in names:
        best[nm] = min(((v, rmse(arm, [nm], obs, v)) for v in grid),
                       key=lambda vr: vr[1])[0]
    return best


def trend(per, n_C, names):
    """corr(n_C, log10 g_apo) and the max/min spread of the fitted values."""
    v = np.array([per[nm] for nm in names], float)
    n = np.array([n_C[nm] for nm in names], float)
    spread = float(v.max() / max(v.min(), 1e-12))
    corr = float(np.corrcoef(n, np.log10(v))[0, 1]) if v.min() > 0 else float("nan")
    return corr, spread


def _boot_smaller(a, b, n=4000, seed=0):
    """P(RMSE from residual set `a` < that from `b`) resampling congener blocks."""
    rng = np.random.default_rng(seed)
    A, B = a.reshape(-1, 3), b.reshape(-1, 3)          # one row per congener
    rm = lambda r: np.sqrt(np.mean(r ** 2))            # noqa: E731
    w = 0
    for _ in range(n):
        i = rng.integers(0, len(A), len(A))
        w += rm(A[i]) < rm(B[i])
    return w / n


def main(fast=False):
    obs = load_obs()
    names = [n for n in obs if n in api._CONG]
    n_C = {nm: api._CONG[nm]["n_C"] for nm in names}
    names.sort(key=lambda nm: n_C[nm])
    grid = ([0.5, 2.0, 5.0, 20.0, 50.0, 200.0] if fast else
            [0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 500.0])

    print("=" * 84)
    print("IS THE CHAIN-LENGTH DEPENDENCE AT ENTRY, OR IN SEQUESTRATION?")
    print("handoff B2 — fitting the same entry term against four root models")
    print("=" * 84)
    print("   Read the PRE-REGISTRATION in this file's header before the VERDICT.")
    print(f"   {len(names)} congeners x {len(KEYS)} tissues; demo forcings; monotone")
    print("   physical f_xy; carrier OFF in every arm, so only the ROOT model differs.")

    arms = {"S": "single pool", "S+L": "single pool + lipid",
            "T": "two-pool k_seq + lipid", "C": "two-pool k_seq=0 + lipid (CONTROL)"}
    res = {}
    for arm, label in arms.items():
        gv, rg = fit_global(arm, names, obs, grid)
        per = fit_per_congener(arm, names, obs, grid)
        rp = rmse(arm, names, obs, per)
        corr, spread = trend(per, n_C, names)
        res[arm] = dict(g_global=gv, rmse_global=rg, per=per, rmse_per=rp,
                        corr=corr, spread=spread, penalty=rg - rp)
        print(f"\n{arm:4s} — {label}")
        print("      per-congener g_apo: " +
              "  ".join(f"{nm}:{per[nm]:g}" for nm in names))
        print(f"      corr(n_C, log10 g_apo) = {corr:+.3f}   spread = {spread:.1f}x")
        print(f"      one global g_apo = {gv:g} -> RMSE {rg:.3f}; "
              f"per-congener -> {rp:.3f}; penalty {rg - rp:+.3f}")

    # -- gate ----------------------------------------------------------------
    print("\n" + "-" * 84)
    print("G0. GATE — does the single pool show the trend on THESE forcings?")
    s = res["S"]
    ok = (s["corr"] >= GATE_CORR) and (s["spread"] >= GATE_SPREAD)
    print(f"    S: corr {s['corr']:+.3f} (need >= {GATE_CORR:+.2f}), "
          f"spread {s['spread']:.1f}x (need >= {GATE_SPREAD:g}x) -> "
          f"{'PASS' if ok else 'FAIL — INCONCLUSIVE'}")

    # -- rules ---------------------------------------------------------------
    print("\nR1/R2/R3 — the pre-registered comparisons")
    print(f"    {'arm':5s}{'corr':>9s}{'spread':>9s}{'penalty':>10s}"
          f"{'corr vs S':>11s}{'spread vs S':>13s}")
    for arm in arms:
        a = res[arm]
        print(f"    {arm:5s}{a['corr']:+9.3f}{a['spread']:9.1f}{a['penalty']:+10.3f}"
              f"{a['corr'] / s['corr'] if s['corr'] else float('nan'):11.2f}"
              f"{np.log10(a['spread']) / np.log10(s['spread']) if s['spread'] > 1 else float('nan'):13.2f}")
    print("\n    bootstrap, per-congener vs global within each arm "
          "(how much the chain-length freedom is worth):")
    for arm in arms:
        rp = residuals(arm, names, obs, res[arm]["per"])
        rg = residuals(arm, names, obs, res[arm]["g_global"])
        print(f"      {arm:5s} P(per-congener beats global) = "
              f"{_boot_smaller(rp, rg):.3f}")

    # -- verdict -------------------------------------------------------------
    print("\n" + "=" * 84)
    print("VERDICT")
    print("=" * 84)
    print("   (filled in by the results commit — the pre-registration ships first)")
    return res


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--fast", action="store_true", help="coarser g_apo grid")
    main(**vars(p.parse_args()))
