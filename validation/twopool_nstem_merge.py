#!/usr/bin/env python3
# =============================================================================
# validation/twopool_nstem_merge.py
# -----------------------------------------------------------------------------
# THE STRUCTURAL MERGE: two-pool sequestration ROOT + redistributed (N-stem+leaf)
# SHOOT -- and the per-organ out-of-sample test it finally makes fair.
#
# Why this exists (docs/HANDOFF_BAF_twopool.md item #2 / Result 7):
#   * `twopool_root_exploration.py` fixes the ROOT (mobile + sequestered pools,
#     non-K_PL U-shaped k_seq) but keeps the basic 4pool PASS-THROUGH stem. Its
#     Tang per-organ OOS therefore scored 1.40 -- WORSE than the single-pool
#     baselines -- because the stalk TF collapses for a reason that has nothing
#     to do with the root mechanism (`twopool_root_oos_tang.py`: two-pool LEAF
#     RMSE was the BEST of all models; only the stalk dragged it up).
#   * `pfas_rice_plant_module_nstem_leaf.py` fixes the SHOOT (transpiration
#     deposition + retention across N stem segments) but keeps ONE root pool, so
#     it cannot hold a high long-chain root BAF while still feeding the shoot.
# Neither half could be tested on the other's endpoint. This merges them
# (`model_api.simulate_twopool_nstem`) and runs the pre-registered sequence:
#
#   Stage 1  TRANSFER the cached two-pool fit onto the new shoot, unchanged
#            -> quantifies how much of that fit was shoot-specific.
#   Stage 2  RE-FIT the merged model on Yamazaki (in-sample), by the SAME staged
#            method as the two-pool arc: 7 globals with a linear k_seq -> per
#            congener root-match -> fit the asymmetric-U k_seq(n) -> plug back.
#            Comparable to the two-pool headline (log10 RMSE 0.251, root 0.156).
#   Stage 3  Tang 2026 per-organ TF (stalk/leaf/endosperm, dw, 0.1 ug/g) OOS,
#            transferring the Stage-2 Yamazaki fit with NO Tang re-fit, against
#            the documented baselines (monotone 1.23 / lipid 0.52 / Tang-refit
#            in-sample 0.52 / two-pool-with-pass-through-stem 1.40).
#
# Honest framing (carry forward): the fit is Yamazaki in-sample; Tang is a single
# independent set; GenX f_xy_recommended (ether offset) is provisional and
# over-predicts; grain/endosperm is structurally under (units doc); PFDoDA is a
# near-MQL outlier and is excluded from the fit (reported in the RMSE).
# EXPLORATORY / opt-in: parameters.json, simulate() and reproduce_demo UNCHANGED.
#
# FORCINGS: the MEASURED family (forcing_rice Q_TP + growth_rice ORYZA-IR72 organ
# biomass), i.e. the same forcings as validation/twopool_root_measured.py (two-pool
# reference in-sample RMSE 0.278) and the same family the Tang run uses. This matters
# structurally: the redistributed shoot is TRANSPIRATION-DEPOSITION fed, so it is far
# more sensitive to Q_TP than the phloem-fed 4pool shoot -- on the demo forcings
# (peak Q_TP 0.40 L/d/hill, ~4x the measured 0.098) the deposition route floods the
# grain and the fit collapses onto its bounds. `--demo` keeps the old comparison.
#
#   python validation/twopool_nstem_merge.py            (re-fit, ~10 min)
#   python validation/twopool_nstem_merge.py --cached   (reuse the cached fit)
#   python validation/twopool_nstem_merge.py --demo     (demo forcings instead)
# =============================================================================
from __future__ import annotations
import json, os, sys
import numpy as np
from scipy.optimize import least_squares, brentq

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT_DIR, "src"))
sys.path.insert(0, HERE)

import forcing_rice as fr                      # measured transpiration Q_TP(t)
from growth_rice import organ_biomass          # ORYZA IR72 organ biomass M(t)
import twopool_root_exploration as TP          # obs, k_seq descriptors, fit machinery
import model_api as api                        # the merged model + Tang endpoints

FIT_CACHE = os.path.join(HERE, "twopool_nstem_fitted_params.json")
OBS, CONGENERS = TP.OBS, TP.CONGENERS
DOSE = "low"                                   # Tang 0.1 ug/g (environmentally closest)
SEASON = 120.0

# shoot levers: the documented nstem_leaf defaults, NOT point-fit to any dataset
SHOOT = dict(N=4, stem_transp_frac=0.45, lam_grain=0.05, retention=0.6)


_TP_DEMO_FORCINGS = (TP.T, TP.CWO, TP.QTP, TP.MMAT, TP._dM)   # saved before any swap


def _restore_demo_forcings():
    TP.T, TP.CWO, TP.QTP, TP.MMAT, TP._dM = _TP_DEMO_FORCINGS


def install_forcings(kind="measured"):
    """Build the driver dict AND install the same forcings into TP, so the merged
    model and the two-pool baseline are compared on identical drivers."""
    if kind == "demo":
        return dict(t=TP.T, Cwo=TP.CWO, Qtp=TP.QTP, M=TP.MMAT)
    t = np.linspace(0.0, SEASON, 241)
    Qtp = fr.Q_TP(t, SEASON)
    b = organ_biomass(t, SEASON)
    M = np.maximum(np.column_stack([b["root"], b["stem"], b["leaf"], b["grain"]]), 1e-4)
    TP.T, TP.CWO, TP.QTP, TP.MMAT = t, np.full_like(t, 1.0), Qtp, M
    TP._dM = np.gradient(M, t, axis=0)
    return dict(t=t, Cwo=np.full_like(t, 1.0), Qtp=Qtp, M=M)


DRIVERS = install_forcings("measured")


# ---------------------------------------------------------------------------
# merged forward model (thin wrapper over the API -- no duplicated ODE)
# ---------------------------------------------------------------------------
def simulate_merged(c, p, kseq, drivers=DRIVERS, **kw):
    """Merged two-pool-root + redistributed-shoot run for congener dict `c`.

    `p` holds the transport globals (kappa_d, L_Ph, gxy/gph lipid gate); `kseq`
    is the root sequestration rate [1/day]. Returns the API result dict.
    """
    gxy, gph = TP.lipid_g(c["K_PL_Lkg"], c["group"], p["gxy"], p["gph"],
                          p["K_half"], p["pfsa_ln"])
    opts = dict(SHOOT, **kw)
    return api.simulate_nstem_leaf(
        c["name"], f_xy_source="recommended", K_cw_organ="root",
        kappa_d_override=p["kappa_d"], L_Ph_override=p["L_Ph"],
        g_xy_override=gxy, g_ph_override=gph,
        k_seq=float(kseq), drivers=drivers, **opts)


def tissues(c, p, kseq, **kw):
    """(root, straw, grain) BAF -- the Yamazaki aggregate endpoints."""
    r = simulate_merged(c, p, kseq, **kw)
    return (r["baf_final"]["root"], r["straw_baf"], r["baf_final"]["grain"])


# ---------------------------------------------------------------------------
# Stage 2 fit machinery (mirrors twopool_root_exploration: same params, bounds,
# staging -- only the shoot underneath differs)
# ---------------------------------------------------------------------------
P0 = dict(TP.P0)
FIT = list(TP.FIT)                     # kappa_d, L_Ph, gxy, gph, ks0, ks_b, ks_sa
LOGFIT, BOUNDS = TP.LOGFIT, TP.BOUNDS


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
        lo.append(np.log10(a) if n in LOGFIT else a)
        hi.append(np.log10(b) if n in LOGFIT else b)
    return lo, hi


def _pred_obs(c, p, kseq):
    r, s, g = tissues(c, p, kseq)
    return {"root": r, "straw": s, "grain": g}, OBS[c["name"]]


def residuals(x, drop_pfdoda=True):
    p = _unpack(x)
    res = []
    for c in CONGENERS:
        if drop_pfdoda and c["name"] == "PFDoDA":
            continue
        ks = TP.k_seq(c["n_C"], c["group"], p["ks0"], p["ks_b"], p["ks_sa"])
        pred, o = _pred_obs(c, p, ks)
        for k in ("root", "straw", "grain"):
            if k in o:
                res.append(np.log10(max(pred[k], 1e-6)) - np.log10(o[k]))
    return res


def rmse_report(p, kseq_of):
    """log10 RMSE overall + by tissue for a (p, k_seq(c)) parameterisation."""
    err = {"root": [], "straw": [], "grain": []}
    rows = []
    for c in CONGENERS:
        ks = kseq_of(c)
        pred, o = _pred_obs(c, p, ks)
        for k in err:
            if k in o:
                err[k].append((np.log10(max(pred[k], 1e-6)) - np.log10(o[k])) ** 2)
        rows.append((c, ks, pred, o))
    allsq = err["root"] + err["straw"] + err["grain"]
    by = {k: float(np.sqrt(np.mean(v))) for k, v in err.items()}
    return float(np.sqrt(np.mean(allsq))), by, rows


def _print_rows(rows):
    print(f"{'PFAS':8}{'nC':>3}{'grp':>5}{'k_seq':>9} | "
          f"{'root p/o':>15}{'straw p/o':>15}{'grain p/o':>15}")
    for c, ks, pred, o in rows:
        nan = float("nan")
        print(f"{c['name']:8}{c['n_C']:>3}{c['group'][-3:]:>5}{ks:>9.3f} | "
              f"{pred['root']:>7.2f}/{o.get('root', nan):<6.2f} "
              f"{pred['straw']:>7.2f}/{o.get('straw', nan):<6.2f} "
              f"{pred['grain']:>7.2f}/{o.get('grain', nan):<6.2f}")


def root_match(p):
    """Back out, per congener, the k_seq that makes merged root == observed root."""
    out = []
    for c in CONGENERS:
        tgt = OBS[c["name"]]["root"]

        def f(logk):
            r, _, _ = tissues(c, p, 10.0 ** logk)
            return np.log10(max(r, 1e-6)) - np.log10(tgt)

        lo, hi = -6.0, 3.0
        flo, fhi = f(lo), f(hi)
        if flo > 0:            # even k_seq~0 over-shoots the root
            ks = 0.0
        elif fhi < 0:          # even a huge sink cannot reach the observed root
            ks = 10.0 ** hi
        else:
            ks = 10.0 ** brentq(f, lo, hi, xtol=1e-3)
        out.append((c["name"], c["n_C"], c["group"], ks))
    return out


# ---------------------------------------------------------------------------
# Stage 3: Tang 2026 per-organ OOS (dry-weight TF, same convention as
# model_api.tang_tf_validation so the baselines are directly comparable)
# ---------------------------------------------------------------------------
def merged_tang_tf(name, p, q, season=150.0):
    """Merged-model per-organ DRY-weight TF on the standard Tang forcings."""
    c = api._CONG[name]
    ks = TP.kseq_ushape(c["n_C"], c["group"], q)
    gxy, gph = TP.lipid_g(c["K_PL_Lkg"], c["group"], p["gxy"], p["gph"],
                          p["K_half"], p["pfsa_ln"])
    r = api.simulate_nstem_leaf(
        name, f_xy_source="recommended", K_cw_organ="root", season=season,
        kappa_d_override=p["kappa_d"], L_Ph_override=p["L_Ph"],
        g_xy_override=gxy, g_ph_override=gph, k_seq=float(ks), **SHOOT)
    froot = 1.0 - api._COMP["root"]["theta_fw"]
    return {org: r["tf_final"][mk] * froot / (1.0 - api._COMP[tk]["theta_fw"])
            for mk, org, tk in api._TANG_ORGANS}


def _rmse(pairs):
    e = [(np.log10(max(pr, 1e-6)) - np.log10(max(ob, 1e-6))) ** 2 for pr, ob in pairs]
    return float(np.sqrt(np.mean(e))) if e else float("nan")


def tang_oos(p, q, dose=DOSE):
    print("\n" + "=" * 86)
    print("STAGE 3 — Tang 2026 per-organ TF, OUT-OF-SAMPLE (no Tang re-fit)")
    print("=" * 86)
    print(f"{'cong':6}{'organ':>11}{'obs':>8}{'merged':>9}{'2pool':>9}"
          f"{'mono':>9}{'lipid':>9}{'refit(IS)':>11}")
    acc = {k: [] for k in ("merged", "2pool", "mono", "lipid", "refit")}
    # the pass-through-stem two-pool baseline is reproduced exactly as
    # twopool_root_oos_tang.py publishes it (its own demo forcings + cached fit)
    saved = (TP.T, TP.CWO, TP.QTP, TP.MMAT, TP._dM)
    _restore_demo_forcings()
    p2, q2 = TP.load_fit()
    import twopool_root_oos_tang as T2
    tf_2pool_all = {nm: T2.twopool_tf(api._CONG[nm], p2, q2) for nm in api.TANG_CONGENERS}
    TP.T, TP.CWO, TP.QTP, TP.MMAT, TP._dM = saved
    for nm in api.TANG_CONGENERS:
        obs = api.tang_observed_tf(nm, dose)
        tf_merged = merged_tang_tf(nm, p, q)
        tf_2pool = tf_2pool_all[nm]
        v_m = api.tang_tf_validation(nm, f_xy_source="recommended", dose=dose)
        v_l = api.tang_tf_validation(nm, f_xy_source="recommended", dose=dose,
                                     lipid_loading=True)
        v_r = api.tang_tf_validation(nm, f_xy_source="recommended", use_refit=True,
                                     dose=dose)
        for _, organ, _ in api._TANG_ORGANS:
            o = obs[organ]
            vals = (tf_merged[organ], tf_2pool[organ], v_m["model_tf"][organ],
                    v_l["model_tf"][organ], v_r["model_tf"][organ])
            for k, v in zip(("merged", "2pool", "mono", "lipid", "refit"), vals):
                acc[k].append((v, o))
            print(f"{nm:6}{organ:>11}{o:>8.3f}" + "".join(f"{v:>9.3f}" for v in vals[:4])
                  + f"{vals[4]:>11.3f}")
    rmses = {k: _rmse(v) for k, v in acc.items()}
    print("\n" + "-" * 86)
    print(f"log10 RMSE vs Tang per-organ TF (dw, dose={dose}):")
    print(f"  MERGED 2pool-root + redistributed shoot (OOS) = {rmses['merged']:.3f}  <- this work")
    print(f"  two-pool root, pass-through stem       (OOS) = {rmses['2pool']:.3f}")
    print(f"  single-pool monotone                   (OOS) = {rmses['mono']:.3f}")
    print(f"  single-pool lipid loading              (OOS) = {rmses['lipid']:.3f}")
    print(f"  single-pool Tang-refit f_xy             (IS) = {rmses['refit']:.3f}")
    organs = [o for _, o, _ in api._TANG_ORGANS]
    print(f"\n  per-organ RMSE      {'stalk':>8}{'leaf':>8}{'endosperm':>11}")
    for k in ("merged", "2pool", "mono", "lipid"):
        vals = [_rmse(acc[k][i::3]) for i in range(len(organs))]
        print(f"    {k:16}{vals[0]:>8.2f}{vals[1]:>8.2f}{vals[2]:>11.2f}")
    return rmses


# ---------------------------------------------------------------------------
def _reference_twopool_fit(forcing):
    """The cached two-pool (pass-through-stem) fit made on the SAME forcings, with
    its published in-sample reference RMSE."""
    if forcing == "demo":
        p, q = TP.load_fit()
        return p, q, "0.251 / root 0.156 / straw 0.260 / grain 0.311 (demo forcings)"
    d = json.load(open(os.path.join(HERE, "twopool_fitted_params_measured.json")))
    return (d["global"], np.array(d["ushape_q"]),
            "0.278 / root 0.154 (measured forcings)")


def main(cached=False, forcing="measured"):
    global DRIVERS
    DRIVERS = install_forcings(forcing)
    print("=" * 86)
    print("STRUCTURAL MERGE — two-pool sequestration root + redistributed (N-stem+leaf) shoot")
    print("=" * 86)
    print(f"forcings = {forcing}  (Q_TP peak {DRIVERS['Qtp'].max():.3f} L/d/hill, "
          f"season {DRIVERS['t'][-1]:.0f} d)")

    # --- Stage 1: transfer the cached two-pool fit onto the new shoot ---------
    p_old, q_old, ref = _reference_twopool_fit(forcing)
    print("\nSTAGE 1 — transfer the cached two-pool fit (fit WITH the pass-through stem)")
    r1, by1, rows1 = rmse_report(
        p_old, lambda c: TP.kseq_ushape(c["n_C"], c["group"], q_old))
    _print_rows(rows1)
    print(f"\n  in-sample log10 RMSE = {r1:.3f}   root={by1['root']:.3f} "
          f"straw={by1['straw']:.3f} grain={by1['grain']:.3f}")
    print(f"  (same fit with its OWN pass-through-stem 4pool shoot: {ref})")
    print("  => the gap is the part of that fit which was SHOOT-specific: the redistributed")
    print("     shoot feeds the grain by panicle deposition + residual xylem, not by the")
    print("     4pool's phloem-only route, so L_Ph/g_ph do not carry over. Hence Stage 2.")

    # --- Stage 2: re-fit the merged model on Yamazaki ------------------------
    if cached and os.path.exists(FIT_CACHE):
        d = json.load(open(FIT_CACHE))
        p, q = d["global"], np.array(d["ushape_q"])
        print(f"\nSTAGE 2 — loaded cached merged fit ({FIT_CACHE})")
    else:
        print("\nSTAGE 2 — re-fitting the merged model on Yamazaki (7 globals, "
              "PFDoDA excluded) ...")
        lo, hi = _bounds()
        sol = least_squares(residuals, _x0(), bounds=(lo, hi), method="trf",
                            diff_step=1e-2, max_nfev=400, args=(True,))
        p = _unpack(sol.x)
        print("  fitted globals:  " + "  ".join(f"{n}={p[n]:.4g}" for n in FIT))
        demanded = root_match(p)
        print("  root-matched empirical k_seq: "
              + "  ".join(f"{nm}={ks:.3f}" for nm, _, _, ks in demanded))
        q = TP.fit_ushape(demanded)
        json.dump({"global": p, "ushape_q": list(map(float, q)),
                   "shoot": SHOOT, "forcing": forcing,
                   "note": "merged two-pool root + nstem_leaf redistributed shoot; "
                           "Yamazaki in-sample; EXPLORATORY (parameters.json UNCHANGED)"},
                  open(FIT_CACHE, "w"), indent=2)
        print(f"  saved -> {FIT_CACHE}")

    print(f"\n  U-shaped k_seq = [{10**q[0]:.3f}*exp(-{q[1]:.2f}(n-4)) + "
          f"{10**q[2]:.4f}*exp({q[3]:.2f}(n-12))] * {{10^{q[4]:+.2f} if PFSA}}")
    r2, by2, rows2 = rmse_report(p, lambda c: TP.kseq_ushape(c["n_C"], c["group"], q))
    _print_rows(rows2)
    print(f"\n  in-sample log10 RMSE (all 11, incl PFDoDA) = {r2:.3f}   "
          f"root={by2['root']:.3f} straw={by2['straw']:.3f} grain={by2['grain']:.3f}")
    ks_os = TP.kseq_ushape(8, "PFSA", q); ks_un = TP.kseq_ushape(11, "PFCA", q)
    print(f"  PFOS(C8) k_seq={ks_os:.3f} vs PFUnDA(C11) k_seq={ks_un:.3f} "
          f"({ks_un/max(ks_os,1e-9):.1f}x at identical K_PL=31623) — non-K_PL separation")

    # --- Stage 3: the fair per-organ OOS ------------------------------------
    rmses = tang_oos(p, q)

    print("\n" + "=" * 86)
    print("VERDICT (data-driven; read with the honest limits below)")
    print("=" * 86)
    best = min(("merged", "2pool", "mono", "lipid"), key=lambda k: rmses[k])
    print(f"  merged {rmses['merged']:.3f} vs pass-through-stem two-pool {rmses['2pool']:.3f} "
          f"(the unfair comparison Result 7 flagged);")
    print(f"  single-pool baselines: monotone {rmses['mono']:.3f}, lipid {rmses['lipid']:.3f}; "
          f"best OOS = {best}.")
    print("  HONEST LIMITS: Yamazaki in-sample fit; Tang is a single independent set of 3")
    print("  congeners (C5-C8) so the long-chain root decoupling is still not exercised;")
    print("  GenX (provisional ether f_xy) over-predicts; endosperm is structurally under.")
    print("  parameters.json / simulate() / reproduce_demo UNCHANGED (opt-in module).")
    return p, q, rmses


if __name__ == "__main__":
    main(cached="--cached" in sys.argv,
         forcing="demo" if "--demo" in sys.argv else "measured")
