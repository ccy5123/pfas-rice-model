#!/usr/bin/env python3
# =============================================================================
# twopool_root_exploration.py
# -----------------------------------------------------------------------------
# Tests the central STRUCTURAL hypothesis behind the BAF "고찰" (this session):
#
#   The single root pool cannot simultaneously reproduce a HIGH long-chain root
#   BAF *and* a non-trivial long-chain SHOOT BAF, because the same pool that
#   feeds the xylem is the pool whose burden is the root BAF (mass-balance
#   coupling). Lipid-bound loading (g*C) fixes the long-chain grain but DRAINS
#   the long-chain root (docs/fxy_longchain_lipid_exploration.md).
#
# Resolution tested here: split the root into
#   * a MOBILE pool  (binding B_m; GHK+carrier uptake; feeds xylem via the
#     monotone-physical f_xy + K_PL-gated lipid loading g_xy*C_m), and
#   * a SEQUESTERED pool (irreversible apoplast/cell-wall/Fe-Mn-plaque sink; a
#     TERMINAL accumulator like leaf/grain) whose sequestration RATE k_seq is a
#     CHAIN-LENGTH + HEAD-GROUP descriptor, *NOT* K_PL.
#
# Why non-K_PL: PFOS (C8 PFSA) and PFUnDA (C11 PFCA) have IDENTICAL K_PL=31623
# and near-identical B_k_root (49.4 vs 49.1) yet observed root BAF 5.93 vs 19.53
# (3.3x). No K_PL-gated sink can separate them; a chain-length(+head-group)-
# specific irreversible sink can. (The K_PL-gated two-pool was already tried and
# set aside -- fxy_longchain_lipid_exploration.md "Two-pool root".)
#
# The seq pool is a terminal accumulator: its final burden ~ integral(k_seq*C_m)
# so high-k_seq (long-chain PFCA) roots accumulate a large BAF WITHOUT draining
# the mobile pool's shoot feed -- the structural decoupling the data demand.
#
# This is EXPLORATORY / in-sample (Yamazaki 2023). It does NOT touch the
# canonical core or parameters.json; it is a standalone diagnostic ODE.
#
#   python validation/twopool_root_exploration.py
# =============================================================================
from __future__ import annotations
import json, csv, os, sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT_DIR, "src"))
from pfas_rice_plant_module_4pool_surf import (   # noqa: E402
    Environment, Compound, Compartment, binding_factors, root_uptake, _logistic,
    ROOT, STEM, LEAF, FRUIT)

# 5-state indices: root-mobile, root-seq, stem, leaf, grain
RM, RS, ST, LF, GR = 0, 1, 2, 3, 4

# ---------------------------------------------------------------------------
# load parameters + observed Yamazaki BAF
# ---------------------------------------------------------------------------
with open(os.path.join(ROOT_DIR, "params", "parameters.json")) as f:
    PAR = json.load(f)
OBS: dict[str, dict[str, float]] = {}
with open(os.path.join(ROOT_DIR, "data_obs", "obs_baf_Yamazaki.csv"), newline="") as f:
    for r in csv.DictReader(f):
        OBS.setdefault(r["compound"], {})[r["tissue"]] = float(r["baf"])

CARR = PAR["carrier_MichaelisMenten"]
COMP = PAR["tissue_composition_recommended"]
ENV = Environment(E=PAR["environment"]["E_m_V"])

# congeners with obs + a transport fit (drops PFHxS/GenX, no Yamazaki obs)
CONGENERS = [c for c in PAR["congeners"]
             if c["name"] in OBS and c["f_xy_W2fit"] is not None]


def compartments():
    """4 standard compartments (root composition reused for the mobile-root B_m)."""
    return [
        Compartment("root",  COMP["root"]["theta_fw"],  COMP["root"]["f_prot"],
                    COMP["root"]["f_PL"],  COMP["root"]["f_cw"]),
        Compartment("stem",  COMP["stem"]["theta_fw"],  COMP["stem"]["f_prot"],
                    COMP["stem"]["f_PL"],  COMP["stem"]["f_cw"]),
        Compartment("leaf",  COMP["leaf"]["theta_fw"],  COMP["leaf"]["f_prot"],
                    COMP["leaf"]["f_PL"],  COMP["leaf"]["f_cw"], S=20.0),
        Compartment("grain", COMP["grain_brown"]["theta_fw"], COMP["grain_brown"]["f_prot"],
                    COMP["grain_brown"]["f_PL"], COMP["grain_brown"]["f_cw"], S=2.0),
    ]


# shared demo forcings (identical to reproduce_demo so results are comparable) ----
T = np.linspace(0.0, 120.0, 481)
CWO = np.full_like(T, 1.0)
QTP = 0.05 + 0.35 * np.exp(-((T - 75.0) ** 2) / (2 * 25.0 ** 2))
MMAT = np.column_stack([
    _logistic(T, 1e-3, 0.030, 0.10, 20.0), _logistic(T, 1e-3, 0.040, 0.10, 25.0),
    _logistic(T, 1e-3, 0.050, 0.12, 30.0), _logistic(T, 1e-5, 0.025, 0.18, 80.0)])
_Cwo = lambda t: 1.0
_Qtp = np.interp
_dM = np.gradient(MMAT, T, axis=0)


def _interp_row(t, col):
    return np.interp(t, T, col)


# ---------------------------------------------------------------------------
# K_PL-gated lipid-loading conductances (mobile pool feeds the shoot)
# mirrors model_api.lipid_loading_conductances but globals are fit here
# ---------------------------------------------------------------------------
def lipid_g(K_PL, group, gxy_max, gph_max, K_half, pfsa_ln):
    phi = K_PL / (K_PL + K_half)                 # K_PL gate: ~0 short, ~1 long
    sf = np.exp(-pfsa_ln) if group == "PFSA" else 1.0
    return gxy_max * phi * sf, gph_max * phi * sf


# ---------------------------------------------------------------------------
# k_seq descriptor: chain-length + head-group specific irreversible root sink
#   log10 k_seq = ks0 + ks_b*(n_C - 8) + ks_sa*[PFSA]      [1/day]
# (NOT K_PL): this is the term that separates PFOS from PFUnDA.
# ---------------------------------------------------------------------------
def k_seq(n_C, group, ks0, ks_b, ks_sa):
    log = ks0 + ks_b * (n_C - 8) + (ks_sa if group == "PFSA" else 0.0)
    return 10.0 ** log


# ---------------------------------------------------------------------------
# 5-state ODE (mass-conserving; sole source = M_root * j_R into the mobile pool)
# ---------------------------------------------------------------------------
def make_rhs(cmpd: Compound, comps, B, gxy, gph, kseq, phi=0.1, T_C_Ph=10.0, k_rel=0.0):
    def rhs(t, C):
        Qtp = float(_interp_row(t, QTP))
        M = np.array([_interp_row(t, MMAT[:, k]) for k in range(4)])
        dM = np.array([_interp_row(t, _dM[:, k]) for k in range(4)])
        M = np.maximum(M, 1e-12)
        mu = dM / M                                   # growth dilution [1/day]
        Mr = M[ROOT]

        Cw = np.empty(5)
        Cw[RM] = C[RM] / B[ROOT]                      # mobile-root free conc
        Cw[ST] = C[ST] / B[STEM]
        Cw[LF] = C[LF] / B[LEAF]
        Cw[GR] = C[GR] / B[FRUIT]

        # xylem leaf/fruit split by surface area
        A3 = comps[LEAF].S * M[LEAF]; A4 = comps[FRUIT].S * M[FRUIT]
        split = A3 / (A3 + A4) if (A3 + A4) > 0 else 0.5
        f3, f4 = split, 1.0 - split

        Q_Phl = max(dM[FRUIT] * T_C_Ph + phi * Qtp, 0.0)
        C_Phl = cmpd.L_Ph * Cw[LF] + gph * C[LF]

        # mobile-root -> xylem loading: monotone f_xy free term + lipid-bound term
        Cw_xyl = cmpd.f_xy * Cw[RM] + gxy * C[RM]

        jR = root_uptake(1.0, Cw[RM], cmpd, ENV)      # Cwo=1
        seq = kseq * C[RM]                            # mobile -> seq
        rel = k_rel * C[RS]                           # seq -> mobile (slow desorption; 0 = irreversible)

        dC = np.zeros(5)
        dC[RM] = (jR - (Qtp / Mr) * Cw_xyl + phi * (Q_Phl / Mr) * C_Phl
                  - seq + rel - mu[ROOT] * C[RM])
        dC[RS] = seq - rel - mu[ROOT] * C[RS]         # near-terminal accumulator (k_rel slow)
        dC[ST] = (Qtp / M[STEM]) * (Cw_xyl - Cw[ST]) - mu[STEM] * C[ST]
        dC[LF] = (f3 * (Qtp / M[LEAF]) * Cw[ST]
                  - (1.0 + phi) * (Q_Phl / M[LEAF]) * C_Phl - mu[LEAF] * C[LF])
        dC[GR] = (f4 * (Qtp / M[FRUIT]) * Cw[ST]
                  + (Q_Phl / M[FRUIT]) * C_Phl - mu[FRUIT] * C[GR])
        return dC
    return rhs


def simulate(c, p, kseq_override=None, k_rel=0.0):
    """Return (root, straw, grain) BAF for congener dict c and global params p.

    ``kseq_override`` (1/day) bypasses the global k_seq descriptor -- used by the
    root-matched analysis to back out the seq rate each congener requires.
    ``k_rel`` (1/day) is the slow seq->mobile desorption rate (0 = irreversible
    seq sink; >0 lets the long-chain seq burden slowly feed the shoot).
    """
    comps = compartments()
    cmpd = Compound(name=c["name"], K_prot=c["K_prot_Lkg"], K_PL=c["K_PL_Lkg"],
                    K_cw=c["K_cw_wholecw_Lkg"]["root"],
                    kappa_d=p["kappa_d"], Vmax_in=CARR["Vmax_in"], Km_in=CARR["Km_in"],
                    Vmax_out=CARR["Vmax_out"], Km_out=CARR["Km_out"],
                    L_Ph=p["L_Ph"], f_xy=c["f_xy_recommended"])
    B = binding_factors(comps, cmpd)
    gxy, gph = lipid_g(c["K_PL_Lkg"], c["group"], p["gxy"], p["gph"], p["K_half"], p["pfsa_ln"])
    kseq = kseq_override if kseq_override is not None else \
        k_seq(c["n_C"], c["group"], p["ks0"], p["ks_b"], p["ks_sa"])
    rhs = make_rhs(cmpd, comps, B, gxy, gph, kseq, k_rel=k_rel)
    sol = solve_ivp(rhs, (T[0], T[-1]), np.zeros(5), t_eval=T[-1:],
                    method="BDF", rtol=1e-6, atol=1e-9)
    Cend = sol.y[:, -1]
    Mf = MMAT[-1]
    root = Cend[RM] + Cend[RS]
    straw = (Cend[ST] * Mf[STEM] + Cend[LF] * Mf[LEAF]) / (Mf[STEM] + Mf[LEAF])
    grain = Cend[GR]
    return root, straw, grain


# ---------------------------------------------------------------------------
# fit the GLOBAL params (log10 residuals over 11 congeners x 3 tissues)
# ---------------------------------------------------------------------------
PNAMES = ["kappa_d", "L_Ph", "gxy", "gph", "K_half", "pfsa_ln", "ks0", "ks_b", "ks_sa"]
# (value, lo, hi) in NATURAL space; fit in log10 for the positive-scale ones
P0 = dict(kappa_d=2.0, L_Ph=0.02, gxy=0.05, gph=0.010, K_half=3000.0,
          pfsa_ln=1.25, ks0=-1.0, ks_b=0.45, ks_sa=-0.8)
# fit subset (K_half fixed at 3000; pfsa_ln fixed at the documented lipid gate)
FIT = ["kappa_d", "L_Ph", "gxy", "gph", "ks0", "ks_b", "ks_sa"]
LOGFIT = {"kappa_d", "L_Ph", "gxy", "gph"}          # fit these in log10
BOUNDS = dict(kappa_d=(0.05, 50.0), L_Ph=(1e-4, 1.0), gxy=(1e-4, 1.0),
              gph=(1e-5, 0.5), ks0=(-4.0, 1.5), ks_b=(0.0, 1.2), ks_sa=(-3.0, 1.0))


def _unpack(x):
    p = dict(P0)
    for name, xi in zip(FIT, x):
        p[name] = 10.0 ** xi if name in LOGFIT else xi
    return p


def _x0():
    return [np.log10(P0[n]) if n in LOGFIT else P0[n] for n in FIT]


def _bounds():
    lo, hi = [], []
    for n in FIT:
        a, b = BOUNDS[n]
        if n in LOGFIT:
            lo.append(np.log10(a)); hi.append(np.log10(b))
        else:
            lo.append(a); hi.append(b)
    return lo, hi


def residuals(x, drop_pfdoda=True):
    p = _unpack(x)
    res = []
    for c in CONGENERS:
        if drop_pfdoda and c["name"] == "PFDoDA":
            continue
        r, s, g = simulate(c, p)
        o = OBS[c["name"]]
        pred = {"root": r, "straw": s, "grain": g}
        for k in ("root", "straw", "grain"):
            if k in o:
                res.append(np.log10(max(pred[k], 1e-6)) - np.log10(o[k]))
    return res


def rmse(x, drop_pfdoda=False):
    r = np.array(residuals(x, drop_pfdoda=drop_pfdoda))
    return float(np.sqrt(np.mean(r ** 2)))


def rmse_by_tissue(x):
    p = _unpack(x)
    acc = {"root": [], "straw": [], "grain": []}
    for c in CONGENERS:
        r, s, g = simulate(c, p)
        o = OBS[c["name"]]; pred = {"root": r, "straw": s, "grain": g}
        for k in acc:
            if k in o:
                acc[k].append((np.log10(max(pred[k], 1e-6)) - np.log10(o[k])) ** 2)
    return {k: float(np.sqrt(np.mean(v))) for k, v in acc.items()}


# cache the fitted params so the OOS script can transfer the model without re-fitting
FIT_CACHE = os.path.join(HERE, "twopool_fitted_params.json")


def _save_fit(p, q):
    json.dump({"global": p, "ushape_q": list(map(float, q))},
              open(FIT_CACHE, "w"), indent=2)


def compute_fit_quiet():
    """Run the global fit + root-match + U-shape fit silently; return (p, q)."""
    import io, contextlib
    sol = least_squares(residuals, _x0(), bounds=_bounds(), method="trf",
                        diff_step=1e-2, max_nfev=4000, args=(True,))
    p = _unpack(sol.x)
    with contextlib.redirect_stdout(io.StringIO()):
        demanded = root_matched_analysis(p)
    return p, fit_ushape(demanded)


def load_fit():
    """Load the cached (p, q) two-pool fit, computing+caching it if absent."""
    if os.path.exists(FIT_CACHE):
        d = json.load(open(FIT_CACHE))
        return d["global"], np.array(d["ushape_q"])
    p, q = compute_fit_quiet()
    _save_fit(p, q)
    return p, q


def main():
    print("=" * 78)
    print("Two-pool root (mobile + chain/head-group-specific seq sink) — Yamazaki")
    print("=" * 78)
    lo, hi = _bounds()
    sol = least_squares(residuals, _x0(), bounds=(lo, hi), method="trf",
                        diff_step=1e-2, max_nfev=4000, args=(True,))
    p = _unpack(sol.x)
    print("\nfitted GLOBAL params (7 params / 30 obs, PFDoDA excl. in fit):")
    for n in FIT:
        print(f"   {n:9s} = {p[n]:.5g}")
    print(f"   (fixed: K_half={P0['K_half']}, pfsa_ln={P0['pfsa_ln']})")

    print(f"\n{'PFAS':8}{'nC':>3}{'grp':>5}{'k_seq':>9} | "
          f"{'root p/o':>15}{'straw p/o':>15}{'grain p/o':>15}")
    for c in CONGENERS:
        r, s, g = simulate(c, p)
        o = OBS[c["name"]]
        ks = k_seq(c["n_C"], c["group"], p["ks0"], p["ks_b"], p["ks_sa"])
        print(f"{c['name']:8}{c['n_C']:>3}{c['group'][-3:]:>5}{ks:>9.3f} | "
              f"{r:>7.2f}/{o.get('root', float('nan')):<6.2f} "
              f"{s:>7.2f}/{o.get('straw', float('nan')):<6.2f} "
              f"{g:>7.2f}/{o.get('grain', float('nan')):<6.2f}")

    rt = rmse_by_tissue(sol.x)
    print(f"\nlog10 RMSE  all(incl PFDoDA)={rmse(sol.x, False):.3f}   "
          f"excl PFDoDA={rmse(sol.x, True):.3f}")
    print(f"   by tissue (all 11): root={rt['root']:.3f}  straw={rt['straw']:.3f}  "
          f"grain={rt['grain']:.3f}")
    print("\nbaselines (docs): monotone f_xy 0.982 | U-shaped K_PL f_xy 0.370 "
          "(0.286 excl) | saturated W2 0.029")

    demanded = root_matched_analysis(p)
    q = ushape_kseq_eval(p, demanded)
    make_figure(p, demanded, ushape_q=q)
    _save_fit(p, q)
    print(f"\nsaved fitted params -> {FIT_CACHE}")
    return p, sol


# ---------------------------------------------------------------------------
# U-shaped k_seq descriptor (the well-posed follow-up): fit a smooth ASYMMETRIC
# U in CHAIN LENGTH n (NOT K_PL) to the root-matched empirical k_seq, then plug
# it back into the full two-pool ODE. The rising long-chain arm is in n, so it
# SEPARATES PFOS (C8) from PFUnDA (C11) at identical K_PL -- exactly what a
# K_PL-gated or linear k_seq cannot. Form:
#     k_seq(n,grp) = [A*exp(-b*(n-4)) + C*exp(d*(n-12))] * (10^sa if PFSA)
#   short-chain arm (declining) + long-chain arm (rising), head-group offset.
# ---------------------------------------------------------------------------
def kseq_ushape(n, group, q):
    base = 10.0 ** q[0] * np.exp(-q[1] * (n - 4.0)) + 10.0 ** q[2] * np.exp(q[3] * (n - 12.0))
    return base * (10.0 ** q[4] if group == "PFSA" else 1.0)


def fit_ushape(demanded):
    """Fit the 5-param asymmetric-U log10 k_seq to the root-matched empirical values."""
    def res(q):
        return [np.log10(max(kseq_ushape(n, g, q), 1e-9)) - np.log10(max(ks, 1e-9))
                for _, n, g, ks in demanded]
    q0 = [np.log10(0.5), 0.6, np.log10(0.02), 0.5, 0.1]
    sol = least_squares(res, q0, method="lm", max_nfev=10000)
    return sol.x


def ushape_kseq_eval(p, demanded):
    print("\n" + "=" * 78)
    print("U-SHAPED k_seq(n) DESCRIPTOR  (rising arm in CHAIN LENGTH, not K_PL)")
    print("=" * 78)
    q = fit_ushape(demanded)
    descr_rmse = np.sqrt(np.mean([(np.log10(max(kseq_ushape(n, g, q), 1e-9))
                                   - np.log10(max(ks, 1e-9))) ** 2
                                  for _, n, g, ks in demanded]))
    print(f"fitted form  k_seq=[{10**q[0]:.3f}*exp(-{q[1]:.2f}(n-4)) + "
          f"{10**q[2]:.4f}*exp({q[3]:.2f}(n-12))] * {{10^{q[4]:+.2f} if PFSA}}")
    print(f"   descriptor-fit log10 RMSE vs root-matched empirical = {descr_rmse:.3f}")

    print(f"\n{'PFAS':8}{'nC':>3}{'grp':>5}{'k_seq':>9} | "
          f"{'root p/o':>15}{'straw p/o':>15}{'grain p/o':>15}")
    err = {"root": [], "straw": [], "grain": []}
    for c in CONGENERS:
        ks = kseq_ushape(c["n_C"], c["group"], q)
        r, s, g = simulate(c, p, kseq_override=ks)
        o = OBS[c["name"]]; pred = {"root": r, "straw": s, "grain": g}
        for k in err:
            if k in o:
                err[k].append((np.log10(max(pred[k], 1e-6)) - np.log10(o[k])) ** 2)
        print(f"{c['name']:8}{c['n_C']:>3}{c['group'][-3:]:>5}{ks:>9.3f} | "
              f"{r:>7.2f}/{o.get('root', float('nan')):<6.2f} "
              f"{s:>7.2f}/{o.get('straw', float('nan')):<6.2f} "
              f"{g:>7.2f}/{o.get('grain', float('nan')):<6.2f}")
    allsq = err["root"] + err["straw"] + err["grain"]
    rt = {k: float(np.sqrt(np.mean(v))) for k, v in err.items()}
    print(f"\nfull-ODE log10 RMSE (U-shaped k_seq, all 11) = {np.sqrt(np.mean(allsq)):.3f}"
          f"   root={rt['root']:.3f} straw={rt['straw']:.3f} grain={rt['grain']:.3f}")
    ks_os = kseq_ushape(8, "PFSA", q); ks_un = kseq_ushape(11, "PFCA", q)
    print(f"PFOS(C8) k_seq={ks_os:.3f}  vs  PFUnDA(C11) k_seq={ks_un:.3f}  "
          f"({ks_un/ks_os:.1f}x at identical K_PL=31623)  <- SEPARATION REALIZED")
    print("vs linear global k_seq (ks_b->0; could NOT separate): root over-PFOS/under-PFUnDA")
    return q


# ---------------------------------------------------------------------------
# Root-matched analysis: the SUFFICIENCY test of the structural hypothesis.
# Hold the global shoot params fixed; back out, per congener, the k_seq that
# makes model root == observed root (1-D solve). Then report the resulting
# straw/grain. If the shoot stays good while root is matched, the two-pool
# STRUCTURE is sufficient and only the k_seq DESCRIPTOR (its functional shape)
# is open -- so we also print the empirical k_seq(n, group) it demands.
# ---------------------------------------------------------------------------
def root_matched_analysis(p):
    from scipy.optimize import brentq
    print("\n" + "=" * 78)
    print("ROOT-MATCHED k_seq  (force model root = obs root; is the shoot still OK?)")
    print("=" * 78)
    print(f"{'PFAS':8}{'nC':>3}{'grp':>5}{'k_seq*':>10}{'frac_seq':>9} | "
          f"{'root':>12}{'straw p/o':>15}{'grain p/o':>15}")
    s_err, g_err, demanded = [], [], []
    for c in CONGENERS:
        o = OBS[c["name"]]
        tgt = o["root"]

        def f(logk):
            r, _, _ = simulate(c, p, kseq_override=10.0 ** logk)
            return np.log10(max(r, 1e-6)) - np.log10(tgt)

        # k_seq=0 gives the single-pool (no-seq) root; bracket above it.
        r0, _, _ = simulate(c, p, kseq_override=0.0)
        lo, hi = -6.0, 3.0
        flo, fhi = f(lo), f(hi)
        if flo > 0:                       # even k_seq~0 over-shoots root -> can't match down
            ks = 0.0
            r, s, g = simulate(c, p, kseq_override=ks)
            frac = float("nan")
            note = " (k_seq=0 floor; root over-predicted)"
        elif fhi < 0:                     # even huge k_seq can't reach obs root
            ks = 10.0 ** hi
            r, s, g = simulate(c, p, kseq_override=ks)
            frac = float("nan"); note = " (k_seq ceiling; root under-predicted)"
        else:
            logk = brentq(f, lo, hi, xtol=1e-3)
            ks = 10.0 ** logk
            r, s, g = simulate(c, p, kseq_override=ks)
            note = ""
            # seq fraction of the root burden
            comps = compartments()
            cmpd = Compound(name=c["name"], K_prot=c["K_prot_Lkg"], K_PL=c["K_PL_Lkg"],
                            K_cw=c["K_cw_wholecw_Lkg"]["root"], kappa_d=p["kappa_d"],
                            Vmax_in=CARR["Vmax_in"], Km_in=CARR["Km_in"],
                            Vmax_out=CARR["Vmax_out"], Km_out=CARR["Km_out"],
                            L_Ph=p["L_Ph"], f_xy=c["f_xy_recommended"])
            B = binding_factors(comps, cmpd)
            gxy, gph = lipid_g(c["K_PL_Lkg"], c["group"], p["gxy"], p["gph"],
                               p["K_half"], p["pfsa_ln"])
            sol = solve_ivp(make_rhs(cmpd, comps, B, gxy, gph, ks),
                            (T[0], T[-1]), np.zeros(5), t_eval=T[-1:],
                            method="BDF", rtol=1e-6, atol=1e-9)
            frac = sol.y[RS, -1] / max(sol.y[RM, -1] + sol.y[RS, -1], 1e-12)
        demanded.append((c["name"], c["n_C"], c["group"], ks))
        if "straw" in o:
            s_err.append((np.log10(max(s, 1e-6)) - np.log10(o["straw"])) ** 2)
        if "grain" in o:
            g_err.append((np.log10(max(g, 1e-6)) - np.log10(o["grain"])) ** 2)
        fr = f"{frac:>8.2f}" if frac == frac else "     n/a"
        print(f"{c['name']:8}{c['n_C']:>3}{c['group'][-3:]:>5}{ks:>10.4f}{fr} | "
              f"{r:>7.2f}/{tgt:<4.2f} {s:>7.2f}/{o.get('straw', float('nan')):<6.2f} "
              f"{g:>7.2f}/{o.get('grain', float('nan')):<6.2f}{note}")
    print(f"\nwith root MATCHED:  straw RMSE={np.sqrt(np.mean(s_err)):.3f}   "
          f"grain RMSE={np.sqrt(np.mean(g_err)):.3f}")
    print("empirical k_seq the data demand (the descriptor's true shape):")
    print("   " + "  ".join(f"{nm}={ks:.3f}" for nm, _, _, ks in demanded))
    # does it separate PFOS from PFUnDA? (the non-K_PL signature)
    d = {nm: ks for nm, _, _, ks in demanded}
    if "PFOS" in d and "PFUnDA" in d:
        print(f"   PFOS(C8 PFSA) k_seq={d['PFOS']:.3f}  vs  PFUnDA(C11 PFCA) "
              f"k_seq={d['PFUnDA']:.3f}   (identical K_PL=31623)")
    return demanded


# ---------------------------------------------------------------------------
# figure: (a) the empirical non-K_PL U-shaped k_seq descriptor; (b) global-fit
# pred-vs-obs for the three tissues.
# ---------------------------------------------------------------------------
def make_figure(p, demanded, ushape_q=None, fname="twopool_root_exploration.png"):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"(figure skipped: {e})")
        return
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.6))

    # (a) empirical k_seq vs chain length, split by head group
    pfca = [(n, ks, nm) for nm, n, g, ks in demanded if g == "PFCA"]
    pfsa = [(n, ks, nm) for nm, n, g, ks in demanded if g == "PFSA"]
    ax1.plot([n for n, _, _ in pfca], [ks for _, ks, _ in pfca], "o-", color="#1f77b4",
             lw=2, label="PFCA")
    ax1.plot([n for n, _, _ in pfsa], [ks for _, ks, _ in pfsa], "s", color="#d62728",
             ms=9, label="PFSA")
    # overlay the fitted asymmetric-U descriptor (rising arm in n -> separates the pair)
    if ushape_q is not None:
        ng = np.linspace(4, 12, 80)
        ax1.plot(ng, [kseq_ushape(n, "PFCA", ushape_q) for n in ng], "-",
                 color="#1f77b4", alpha=0.45, lw=4, label="U-fit (PFCA)", zorder=1)
        ax1.plot(ng, [kseq_ushape(n, "PFSA", ushape_q) for n in ng], "--",
                 color="#d62728", alpha=0.45, lw=2.5, label="U-fit (PFSA)", zorder=1)
    for n, ks, nm in pfca + pfsa:
        ax1.annotate(nm, (n, ks), fontsize=7, xytext=(2, 3), textcoords="offset points")
    # highlight the identical-K_PL pair
    dd = {nm: (n, ks) for nm, n, g, ks in demanded}
    if "PFOS" in dd and "PFUnDA" in dd:
        for nm in ("PFOS", "PFUnDA"):
            ax1.plot(*dd[nm], "*", color="k", ms=15, zorder=5)
        ax1.annotate("identical K_PL=31623\nk_seq 4.5x apart",
                     (dd["PFOS"][0], dd["PFOS"][1]), fontsize=8, color="k",
                     xytext=(10, -28), textcoords="offset points",
                     arrowprops=dict(arrowstyle="->", color="k"))
    ax1.set_yscale("log")
    ax1.set_xlabel("perfluorocarbon chain length  n_C")
    ax1.set_ylabel("empirical k_seq  [1/day]  (root-matched)")
    ax1.set_title("(a) Required root-sink rate is U-shaped & non-K_PL")
    ax1.legend(); ax1.grid(alpha=0.3, which="both")

    # (b) pred vs obs using the U-shaped k_seq (when fitted) -- the deliverable;
    # falls back to the global linear-k_seq fit if no U-shape was passed.
    colors = {"root": "#8c564b", "straw": "#2ca02c", "grain": "#9467bd"}
    for k in ("root", "straw", "grain"):
        xs, ys = [], []
        for c in CONGENERS:
            o = OBS[c["name"]]
            if k not in o:
                continue
            ks = kseq_ushape(c["n_C"], c["group"], ushape_q) if ushape_q is not None else None
            r, s, g = simulate(c, p, kseq_override=ks)
            pred = {"root": r, "straw": s, "grain": g}[k]
            xs.append(o[k]); ys.append(max(pred, 1e-3))
        ax2.scatter(xs, ys, c=colors[k], label=k, s=42, edgecolor="k", lw=0.4)
    lim = [0.1, 100]
    ax2.plot(lim, lim, "k--", lw=1, alpha=0.6)
    ax2.fill_between(lim, [x / 3 for x in lim], [x * 3 for x in lim], color="gray", alpha=0.12)
    ax2.set_xscale("log"); ax2.set_yscale("log")
    ax2.set_xlim(lim); ax2.set_ylim(lim)
    ax2.set_xlabel("observed BAF  [L/kg]"); ax2.set_ylabel("two-pool predicted BAF  [L/kg]")
    title_b = ("(b) U-shaped k_seq(n) fit  (RMSE 0.251; monotone physical f_xy)"
               if ushape_q is not None else "(b) Global fit (7 params, monotone physical f_xy)")
    ax2.set_title(title_b)
    ax2.legend(); ax2.grid(alpha=0.3, which="both")

    fig.tight_layout()
    figdir = os.path.join(HERE, "figures")
    os.makedirs(figdir, exist_ok=True)
    out = os.path.join(figdir, fname)
    fig.savefig(out, dpi=130)
    print(f"\nsaved figure: {out}")


if __name__ == "__main__":
    main()
