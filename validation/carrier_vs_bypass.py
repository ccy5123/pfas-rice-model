#!/usr/bin/env python3
# =============================================================================
# validation/carrier_vs_bypass.py
# -----------------------------------------------------------------------------
# IS THE MICHAELIS-MENTEN CARRIER NECESSARY? The repo's one structural addition
# to the Trapp ionizable-compound cell model, tested against two levers that are
# already inside that framework.
#
# WHY THIS IS THE QUESTION. docs/theory_anchor.tex is explicit that the
# four-compartment DPU (Rein 2011, Brunetti 2019) is, at the membrane level, the
# Trapp (2000, 2004) ionizable-compound cell model, and its symbol table marks
# only three things "--- (new)": f_xy, eta, and the Michaelis-Menten carrier
# (J_max, K_M). f_xy is a lumped stand-in for Trapp's own detailed root model and
# B_k is his K_RW re-parameterised, so the CARRIER is the single piece of physics
# this repo adds. It was added because Trapp's model is degenerate in the PFAS
# limit -- the ion trap needs a neutral form (f_n ~ 1e-5..1e-7 and pH-invariant),
# and GHK suppresses passive anion influx ~100x at every chain length -- so the
# f_n -> 0 limit predicts almost no uptake while PFAS demonstrably accumulates.
# Measured here: switching the carrier off takes PFOA root BAF 0.479 -> 0.038.
#
# But the carrier has never been compared with an alternative. It was FITTED
# (parameters.json: "fixed during W2 fit"), and docs/neutral_dpu_validation.md
# section 4m has since shown that on the neutral side the same anion-entry gap is
# substantially closable by a PASSIVE apoplastic bypass. So the question is not
# "does the carrier work" but "did the extension have to be made at all".
#
# THREE ARMS, ONE GLOBAL UPTAKE PARAMETER EACH. This is a fair comparison and it
# was worth checking rather than assuming: `Vmax_in` is NOT per-congener. It is a
# single global value (20.0) in parameters.json, so the carrier has exactly the
# same number of free parameters as the two alternatives -- one.
#
#   A. CARRIER          vmax_scale=1, g_apo=0,     E_m=-120 mV   (the incumbent)
#   B. BYPASS           vmax_scale=0, g_apo FIT,   E_m=-120 mV   (new term, but of
#                       the DPU base's own Fickian form: g_apo*(Cwo - Cw))
#   C. DEPOLARISATION   vmax_scale=0, g_apo=0,     E_m FIT       (NO new term at
#                       all -- E_m is already in Trapp's GHK)
#
# Arm C is the strongest form of the question, because it adds nothing whatever:
# parameters.json records E_m_plausible_range_V = [-0.12, -0.09] with the note
# that NH4+ in flooded paddy depolarises the membrane, relaxing e^N from 107 to
# 33. That range is a PRIOR, not a free parameter, and the fit is confined to it.
#
# HELD FIXED ACROSS ARMS so that only the uptake term differs: the monotone
# physical `f_xy_source="recommended"` (a-priori, NOT fitted to Yamazaki), B_k,
# the biomass driver and the measured forcings. The absolute level is therefore
# a-priori-limited for every arm -- the documented a-priori error is ~0.84-0.98
# log10 RMSE -- and NO arm will look good. That is expected and it is not what is
# being measured: the comparison is BETWEEN arms, on identical everything else,
# exactly as section 4l's driver mismatch cancels in an ON-vs-OFF comparison.
#
# PRE-REGISTERED, so the answer cannot be chosen after seeing it:
#
#   1. If arm B needs a chain-length-DEPENDENT g_apo to work, it is a relabelled
#      carrier, not a bypass. theory_anchor.tex states eta (which contains the
#      apoplastic bypass) is "essentially independent of tail length", so a
#      per-congener fit that comes out flat SUPPORTS it and one that has to trend
#      with n REFUTES it. Section 3 fits g_apo per congener to check.
#   2. The carrier arm is the INCUMBENT and its Vmax was fitted on this very
#      dataset. It therefore starts with an in-sample advantage that arms B and C
#      do not have. If B or C merely ties it, that is a win for B/C on parsimony.
#   3. A single global lever is ALREADY KNOWN to fail at long chain: LC6
#      (runs/pfas-rice-carrier, REFUTED) found the per-congener Vmax multiplier
#      needed is ~flat to C10 then 2.0x at C11 and 5.5x at C12. So expect ALL
#      THREE arms to miss C11-C12. If they miss EQUALLY, the long-chain residual
#      is not about the entry mechanism at all and this test cannot speak to it --
#      report that rather than reading a winner out of the long-chain end.
#
#   python validation/carrier_vs_bypass.py           # ~12 min
#   python validation/carrier_vs_bypass.py --fast    # coarse grids, ~4 min
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

import model_api as api                                       # noqa: E402

OBS = os.path.join(ROOT_DIR, "data_obs", "obs_baf_Yamazaki.csv")
KEYS = ("root", "straw", "grain")
SEASON = 120.0
E_M_DEFAULT = -120.0
E_M_PLAUSIBLE = (-120.0, -90.0)      # parameters.json environment.E_m_plausible_range_V


def load_obs(path=OBS):
    obs = {}
    with open(path) as f:
        for r in csv.DictReader(f):
            obs.setdefault(r["compound"], {})[r["tissue"]] = float(r["baf"])
    return {k: v for k, v in obs.items() if all(t in v for t in KEYS)}


def predict(name, *, vmax_scale=1.0, g_apo=0.0, E_m_mV=E_M_DEFAULT):
    """The three tissue BAFs for one congener under one arm's uptake term."""
    r = api.simulate(name, f_xy_source="recommended", season=SEASON,
                     vmax_scale=vmax_scale, g_apo=g_apo, E_m_mV=E_m_mV)
    return (r["baf_final"]["root"], r["straw_baf"], r["baf_final"]["grain"])


def rmse(names, obs, **kw):
    """log10 RMSE over every congener x tissue, the repo's standard PFAS metric."""
    res = []
    for nm in names:
        pred = predict(nm, **kw)
        for t, p in zip(KEYS, pred):
            res.append(np.log10(max(p, 1e-12)) - np.log10(obs[nm][t]))
    return float(np.sqrt(np.mean(np.square(res))))


def per_congener_rmse(names, obs, **kw):
    out = {}
    for nm in names:
        pred = predict(nm, **kw)
        d = [np.log10(max(p, 1e-12)) - np.log10(obs[nm][t])
             for t, p in zip(KEYS, pred)]
        out[nm] = float(np.sqrt(np.mean(np.square(d))))
    return out


def scan(names, obs, key, grid, **fixed):
    """1-D scan of one arm's single parameter. Returns [(value, rmse), ...]."""
    return [(v, rmse(names, obs, **{key: v}, **fixed)) for v in grid]


def _best(curve):
    return min(curve, key=lambda vr: vr[1])


def main(fast=False):
    obs = load_obs()
    names = [n for n in obs if n in api._CONG]
    n_C = {nm: api._CONG[nm]["n_C"] for nm in names}
    names.sort(key=lambda nm: n_C[nm])

    print("=" * 84)
    print("IS THE CARRIER NECESSARY? — the repo's one structural addition to Trapp,")
    print("against two levers already inside his framework")
    print("=" * 84)
    print(f"   {len(names)} congeners x {len(KEYS)} tissues = {len(names)*3} observations")
    print("   Yamazaki 2023; monotone physical f_xy (a-priori, NOT fit); identical")
    print("   forcings and biomass in every arm, so only the uptake term differs.")

    # -- 0 -------------------------------------------------------------------
    print("\n0. THE GAP THE CARRIER FILLS")
    on = rmse(names, obs, vmax_scale=1.0)
    off = rmse(names, obs, vmax_scale=0.0)
    pf_on = predict("PFOA", vmax_scale=1.0)[0]
    pf_off = predict("PFOA", vmax_scale=0.0)[0]
    print(f"   carrier ON  : log10 RMSE {on:.3f}   (PFOA root BAF {pf_on:.4f})")
    print(f"   carrier OFF : log10 RMSE {off:.3f}   (PFOA root BAF {pf_off:.4f})")
    print("   Trapp's model in the PFAS limit, with nothing added, is the OFF line.")

    # -- 1 -------------------------------------------------------------------
    print("\n1. ARM B — APOPLASTIC BYPASS (one global g_apo, carrier off)")
    grid_b = ([0.0, 0.5, 2.0, 5.0, 20.0, 50.0] if fast else
              [0.0, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0])
    curve_b = scan(names, obs, "g_apo", grid_b, vmax_scale=0.0)
    for v, e in curve_b:
        print(f"      g_apo {v:7.2f}   log10 RMSE {e:.3f}")
    b_val, b_rmse = _best(curve_b)
    print(f"   best g_apo = {b_val:g}  ->  {b_rmse:.3f}")

    # -- 2 -------------------------------------------------------------------
    print("\n2. ARM C — MEMBRANE DEPOLARISATION (E_m only; NO new term)")
    print(f"   confined to the recorded plausible range {E_M_PLAUSIBLE} mV.")
    grid_c = ([-120.0, -105.0, -90.0] if fast else
              [-120.0, -114.0, -108.0, -102.0, -96.0, -90.0])
    curve_c = scan(names, obs, "E_m_mV", grid_c, vmax_scale=0.0)
    for v, e in curve_c:
        eN = np.exp(-v / 1000.0 * 96485.0 / (8.314 * 298.15))
        print(f"      E_m {v:7.1f} mV   e^N {eN:6.1f}   log10 RMSE {e:.3f}")
    c_val, c_rmse = _best(curve_c)
    print(f"   best E_m = {c_val:g} mV  ->  {c_rmse:.3f}")

    # -- 3 -------------------------------------------------------------------
    print("\n3. PRE-REGISTERED CHECK 1 — is the fitted bypass chain-length-INDEPENDENT?")
    print("   theory_anchor.tex says eta (which contains the bypass) is 'essentially")
    print("   independent of tail length'. A flat per-congener fit supports that; a")
    print("   trend with n means this is a relabelled carrier, not a bypass.")
    grid_pc = ([0.5, 2.0, 5.0, 20.0, 50.0] if fast else
               [0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0])
    per = {}
    for nm in names:
        cur = [(v, rmse([nm], obs, g_apo=v, vmax_scale=0.0)) for v in grid_pc]
        per[nm] = _best(cur)[0]
        print(f"      {nm:9s} n_C {n_C[nm]:2d}   best g_apo {per[nm]:8.2f}")
    vals = np.array([per[nm] for nm in names], float)
    ncs = np.array([n_C[nm] for nm in names], float)
    spread = float(vals.max() / max(vals.min(), 1e-12))
    r = float(np.corrcoef(ncs, np.log10(vals))[0, 1]) if vals.min() > 0 else float("nan")
    print(f"   spread max/min = {spread:.1f}x ;  corr(n_C, log10 g_apo) = {r:+.3f}")

    # -- 4 -------------------------------------------------------------------
    print("\n4. PRE-REGISTERED CHECK 3 — do the arms differ at the long-chain end?")
    arms = {"A carrier": dict(vmax_scale=1.0),
            "B bypass": dict(vmax_scale=0.0, g_apo=b_val),
            "C depol": dict(vmax_scale=0.0, E_m_mV=c_val)}
    pc = {k: per_congener_rmse(names, obs, **kw) for k, kw in arms.items()}
    print(f"   {'congener':10s}{'n_C':>4s} " + "".join(f"{k:>12s}" for k in arms))
    for nm in names:
        print(f"   {nm:10s}{n_C[nm]:4d} " + "".join(f"{pc[k][nm]:12.3f}" for k in arms))
    longs = [nm for nm in names if n_C[nm] >= 10]
    if longs:
        print("   long chain (n_C >= 10) mean: " +
              "  ".join(f"{k} {np.mean([pc[k][nm] for nm in longs]):.3f}" for k in arms))

    # -- verdict -------------------------------------------------------------
    print("\n" + "=" * 84)
    print("SUMMARY")
    print("=" * 84)
    print(f"   A carrier (incumbent, Vmax fitted on THIS data)  log10 RMSE {on:.3f}")
    print(f"   B bypass  (one global g_apo = {b_val:g})               log10 RMSE {b_rmse:.3f}")
    print(f"   C depol   (E_m = {c_val:g} mV, no new term)          log10 RMSE {c_rmse:.3f}")
    print(f"   nothing   (Trapp's PFAS limit as-is)             log10 RMSE {off:.3f}")
    return dict(names=names, carrier=on, nothing=off, bypass=(b_val, b_rmse),
                depol=(c_val, c_rmse), per_congener_g_apo=per,
                spread=spread, corr_nC=r, arm_rmse=pc)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--fast", action="store_true", help="coarser grids")
    main(**vars(p.parse_args()))
