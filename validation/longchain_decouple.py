#!/usr/bin/env python3
# =============================================================================
# validation/longchain_decouple.py
# -----------------------------------------------------------------------------
# Does a root->shoot DECOUPLING close the long-chain shoot the "complete recipe"
# could not? (Follow-up to validation/longchain_complete.py.)
#
# longchain_complete.py showed the complete recipe (2-pool + LC6 root-matching
# carrier + lipid) closes the long-chain ROOT and GRAIN but OVER-feeds the SHOOT
# (PFDoDA straw 2.3x), because the carrier that fixes the root enlarges the mobile
# pool whose free loading f_xy*Cw_m alone exceeds the observed straw. FINDINGS sec.7
# named the missing piece a root->shoot DECOUPLING: an irreversible root store that
# holds the root burden WITHOUT feeding the xylem.
#
# This tests the simplest such lever on the 2-pool prototype: make the bound store
# IRREVERSIBLE (asymmetric kinetics k_on = ratio*k_off*seq, seq>=1 = stronger
# sequestration into the non-translocating bound pool rb). At a FIXED root-matching
# carrier, scan seq and ask: can it bring the shoot to ~observed while keeping the
# root at ~observed (simultaneous closure within a factor 2)?
#
# The fittable straw is bounded BELOW by the g_xy=0 floor (lipid g_xy>=0 can only
# ADD), so straw_fit_ratio = max(straw_floor/obs, 1) -- no lipid fit is needed for
# the root<->shoot tension; we read the floor directly (fast, deterministic).
#
# RESULT (honest): the lever does NOT decouple. Suppressing the mobile pool to
# protect the shoot LOWERS the internal conc, which RAISES the net uptake gradient,
# so the root burden INFLATES with seq (PFDoDA root 69 -> 481 at seq=10) faster than
# the shoot is relieved (straw 2.3x -> 1.4x). The best balanced point (seq~2) leaves
# BOTH root and shoot ~2.2x off (simultaneity gap ~0.34 log10) -- a partial trade,
# not a clean simultaneous closure. The correct decoupling must instead break the
# UPTAKE<->mobile-conc coupling (e.g. direct deposition into an inert apoplastic
# store that does not raise the influx gradient), a now-sharper open target.
#
#   python validation/longchain_decouple.py
# =============================================================================
import os
import sys
import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "src"))
sys.path.insert(0, HERE)
import twopool_longchain as T  # noqa: E402

K_OFF = 0.02
SERIES = ("PFUnDA", "PFDoDA")           # the two chains the complete recipe over-fed
SEQ_GRID = (1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 10.0)


def simulate3(name, f_xy, g_xy, g_ph, seq=1.0, k_off=K_OFF, vmax_in=None):
    """2-pool root with an IRREVERSIBLE sequestration lever: k_on = ratio*k_off*seq.
    seq=1 reproduces the equilibrium 2-pool (validation/twopool_longchain.simulate2);
    seq>1 traps more burden in the non-translocating bound store rb."""
    c = T.CONG[name]
    cc = T._comps()
    L_Ph = c.get("L_Ph_oryza") or 0.01
    kappa_d = c.get("kappa_d_oryza") or 2.0
    vmax_in = T.carr["Vmax_in"] if vmax_in is None else vmax_in
    cmpd = T.Compound(name=name, K_prot=c["K_prot_Lkg"], K_PL=c["K_PL_Lkg"],
                      K_cw=c["K_cw_wholecw_Lkg"]["root"], kappa_d=kappa_d,
                      Vmax_in=vmax_in, Km_in=T.carr["Km_in"], Vmax_out=T.carr["Vmax_out"],
                      Km_out=T.carr["Km_out"], L_Ph=L_Ph, f_xy=f_xy)
    Bm = T._Bmobile(cc["root"], c)
    Bt = T._Bfull(cc["root"], c)
    ratio = max((Bt - Bm) / Bm, 0.0)
    k_on = ratio * k_off * seq
    Bst = T._Bfull(cc["stem"], c)
    Blf = T._Bfull(cc["leaf"], c)
    S_leaf, S_grain = 20.0, 2.0

    def rhs(tt, C):
        rm, rb, st, lf, gr = C
        Q = float(T.fQ(tt))
        mr, ms, ml, mg = float(T.fMr(tt)), float(T.fMs(tt)), float(T.fMl(tt)), float(T.fMg(tt))
        mur = float(T.dMr(tt)) / mr if mr > 0 else 0.0
        mus = float(T.dMs(tt)) / ms if ms > 0 else 0.0
        mul = float(T.dMl(tt)) / ml if ml > 0 else 0.0
        mug = float(T.dMg(tt)) / mg if mg > 0 else 0.0
        Cw_m, Cw_st, Cw_lf = rm / Bm, st / Bst, lf / Blf
        jR = T.root_uptake(1.0, Cw_m, cmpd, T.env)
        Cw_xyl = f_xy * Cw_m + g_xy * rm
        A3, A4 = S_leaf * ml, S_grain * mg
        sp = A3 / (A3 + A4) if A3 + A4 > 0 else 0.5
        f3, f4 = sp, 1 - sp
        QPh = max(float(T.dMg(tt)) * 10.0 + 0.1 * Q, 0.0)
        CPh = L_Ph * Cw_lf + g_ph * lf
        drm = jR - (Q / mr) * Cw_xyl - k_on * rm + k_off * rb + 0.1 * (QPh / mr) * CPh - mur * rm
        drb = k_on * rm - k_off * rb - mur * rb
        dst = (Q / ms) * (Cw_xyl - Cw_st) - mus * st
        dlf = f3 * (Q / ml) * Cw_st - 1.1 * (QPh / ml) * CPh - mul * lf
        dgr = f4 * (Q / mg) * Cw_st + (QPh / mg) * CPh - mug * gr
        return [drm, drb, dst, dlf, dgr]

    sol = solve_ivp(rhs, (0.0, T.SEASON), np.zeros(5), method="BDF",
                    rtol=1e-5, atol=1e-8, t_eval=[T.SEASON])
    rm, rb, st, lf, gr = sol.y[:, -1]
    mlf, mst = float(T.fMl(T.SEASON)), float(T.fMs(T.SEASON))
    return {"root": rm + rb, "straw": (st * mst + lf * mlf) / (mst + mlf),
            "grain": gr, "rm": rm, "rb": rb}


def _vm_root_match(name):
    """Carrier (Vmax_in) that reproduces the measured root at seq=1 (g_xy=g_ph=0)."""
    c = T.CONG[name]
    o = T.obs[name]
    f_xy = c.get("f_xy_oryza") or c["f_xy_recommended"]
    base = T.carr["Vmax_in"]

    def g(lm):
        return np.log10(max(simulate3(name, f_xy, 0.0, 0.0, 1.0, vmax_in=base * 10 ** lm)["root"], 1e-6)) - np.log10(o["root"])

    a, b = 0.0, 4.0
    if g(a) * g(b) < 0:
        return base * 10 ** brentq(g, a, b, xtol=1e-2, maxiter=40)
    return base * (10 ** a if abs(g(a)) < abs(g(b)) else 10 ** b)


def scan(name):
    """At the fixed root-matching carrier, scan seq; the fittable straw is the
    g_xy=0 floor (lipid can only add), so straw_fit_ratio = max(floor/obs, 1)."""
    c = T.CONG[name]
    o = T.obs[name]
    f_xy = c.get("f_xy_oryza") or c["f_xy_recommended"]
    vm = _vm_root_match(name)
    rows = []
    for seq in SEQ_GRID:
        r = simulate3(name, f_xy, 0.0, 0.0, seq, vmax_in=vm)   # g_xy=0 floor
        root_ratio = r["root"] / o["root"]
        straw_fit_ratio = max(r["straw"] / o["straw"], 1.0)    # lipid g_xy>=0 can only add
        gap = max(abs(math.log10(max(root_ratio, 1e-6))), math.log10(straw_fit_ratio))
        rows.append((seq, root_ratio, straw_fit_ratio, gap, r))
    return vm, rows


def run():
    out = {}
    for name in SERIES:
        vm, rows = scan(name)
        seq1 = rows[0]
        best = min(rows, key=lambda x: x[3])             # min simultaneity gap
        seqmax = rows[-1]
        out[name] = {
            "vm_x_base": vm / T.carr["Vmax_in"],
            "straw_ratio_seq1": seq1[2],
            "root_ratio_seqmax": seqmax[1],
            "root_inflation": seqmax[4]["root"] / seq1[4]["root"],
            "min_gap": best[3],
            "best_seq": best[0],
            "best_root_ratio": best[1],
            "best_straw_ratio": best[2],
            "rows": rows,
        }
    return out


def main():
    res = run()
    for name in SERIES:
        d = res[name]
        o = T.obs[name]
        print(f"\n{name} (obs root {o['root']:.1f} straw {o['straw']:.1f}); "
              f"carrier {d['vm_x_base']:.1f}x base, scan seq (g_xy=0 floor):")
        print(f"  {'seq':>5}{'root x':>9}{'straw x':>9}{'gap(log10)':>12}")
        for seq, rr, sr, gap, _ in d["rows"]:
            print(f"  {seq:>5.1f}{rr:>9.2f}{sr:>9.2f}{gap:>12.3f}")
        print(f"  -> baseline (seq=1) straw {d['straw_ratio_seq1']:.2f}x; "
              f"root inflation seq1->seq{int(SEQ_GRID[-1])} = {d['root_inflation']:.2f}x; "
              f"best simultaneity gap {d['min_gap']:.3f} (factor {10**d['min_gap']:.2f}) "
              f"at seq={d['best_seq']:.1f} (root {d['best_root_ratio']:.2f}x, straw {d['best_straw_ratio']:.2f}x)")
    print("\nVerdict: the irreversible-sequestration lever TRADES root overshoot for shoot relief")
    print("(root inflates as the suppressed mobile pool raises net uptake) but cannot CLOSE both")
    print("within a factor 2 -- best balanced point ~2x on each. Not a clean simultaneous closure.")
    return res


if __name__ == "__main__":
    main()
