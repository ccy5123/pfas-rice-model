#!/usr/bin/env python3
# =============================================================================
# validation/default_switch.py
# -----------------------------------------------------------------------------
# SHOULD THE DEFAULT CONFIGURATION BECOME THE REDISTRIBUTED SHOOT WITH LIPID LOADING?
#
# PRE-REGISTRATION (committed and pushed before the switch is implemented or this
# file's analysis is run; the results are appended below the line at the end).
#
# THE TWO CONFIGURATIONS
#   current default  = simulate(c): the 4-pool core (single straw), free anion,
#                      f_xy_source="recommended" (monotone physical f_xy; L_Ph_rec,
#                      kappa_d_rec from validation/refit_site_year.py Fit B).
#   candidate        = simulate_nstem_leaf(c, lipid_loading=True) EXACTLY AS SCORED by
#                      validation/oos_tang_lipid.py, kim2019_censored.py and
#                      li2025_stage_resolved.py: N=4 stem segments, stem_transp_frac 0.45,
#                      lam_grain 0.05, retention 0.6, K_cw_organ "stem", f_xy_source
#                      "recommended" with the same L_Ph_rec / kappa_d_rec, LIPID_LOADING
#                      (six globals fitted on Yamazaki excl PFDoDA, with the 4-pool core).
#                      NOTHING is re-fitted for the switch: refitting would change the
#                      configuration the transfer scores describe.
#
# WHAT IS ALREADY KNOWN (honest status of this rule)
#   The transfer scores below were computed earlier in this branch, so D1 is a DECLARED
#   decision rule on known numbers, not a blind test. Log10 RMSE, no parameter fitted to
#   these data:
#                                 current default   4-pool lipid   candidate
#     Tang 2026 per-organ TF           1.714            0.992         0.640
#     Kim 2019 grain, point            2.33             0.56          0.47
#     Kim 2019 grain, [0, MDL]         2.05             0.38          0.15
#     Li 2025 shoot TF by stage        1.26             0.90          0.36
#     Li 2025 grain TF by stage        1.68             1.24          0.51
#   NOT yet computed, and therefore genuinely tested here: the candidate on Yamazaki
#   (the two-harvest observable), on all 13 curated congeners, its run time, its own
#   unmeasured-input sensitivity, and the band and inverse it implies.
#
# DECISION RULE
#   D1 (known numbers): on every endpoint above the candidate beats the current default
#      by more than 0.3 log AND is not worse than the 4-pool lipid by more than 0.1.
#   G1 Yamazaki: the candidate's two-harvest dry-basis log10 RMSE over the used cells
#      (root, straw, grain; PFBS grain unused) is not above the current default's (Fit B,
#      1.198). The lipid constants are in-sample here, so this is a guard against a broken
#      configuration, not evidence for it.
#   G2 all 13 curated congeners (12 Yamazaki + GenX) solve, give finite positive BAFs in
#      every tissue, and conserve mass (the nstem_leaf invariant, relative residual < 1e-6).
#   G3 one run takes less than 5 s on this machine (the app re-runs on every change).
#   G4 after implementation, simulate()'s new default reproduces the candidate scores
#      above to 0.005 (it must BE the scored configuration).
#   SWITCH iff D1 and G1-G4 hold. Otherwise the default stays and the reason is recorded.
#
# WHAT A SWITCH CHANGES (decided now, so the run cannot steer it)
#   * simulate() gains shoot="redistributed"|"4pool" (default "redistributed") and
#     lipid_loading defaults to True. shoot="4pool", lipid_loading=False is the old default,
#     kept runnable under that name.
#   * The calibration records keep their configuration explicitly: reproduce_demo.py and
#     refit_site_year.py are the 4-pool fits (Fit A/B and the lipid globals were fitted on
#     the 4-pool core), and the mechanism analyses that were run on the 4-pool core
#     (carrier_vs_bypass, dose_series_carrier, theta_fw_sensitivity, bayesian_inverse_demo,
#     apportionment) are pinned to it rather than silently re-run on another model.
#   * The honesty band (APRIORI_OOS_SCORES) becomes the candidate's transfer errors on the
#     three independent datasets: Tang 0.640, Kim point 0.47, Li shoot 0.36, Li grain 0.51.
#     Yamazaki is left out because the lipid constants were fitted on it. Kim enters with
#     its point score (the conservative one).
#   * The candidate's own unmeasured inputs get the same one-at-a-time sensitivity the
#     current default had (refit_site_year.py UNKNOWNS), on the Yamazaki observable, with
#     the same yardstick (median |log10| difference between the two harvest runs):
#       E_m -90 mV; Km x100; phi 0.05 / 0.2; T_C_Ph 4 / 12; dead blades lost;
#       and the redistributed shoot's three architecture levers, none measured:
#       stem_transp_frac 0.30 / 0.60 (default 0.45), retention 0.40 / 0.80 (default 0.6),
#       lam_grain 0.025 / 0.18 (default 0.05; the grain xylem share the 4-pool sensitivity
#       spans at s_grain 0.2 / 2.0 m2/kg dry). The ranges for stem_transp_frac and
#       retention have no source; they are +/- one third and are labelled as such.
#     s_grain_dry does not apply (the redistributed shoot splits the stream by lam_grain).
#
# EXPECTATIONS (written before the run)
#   E1 D1 holds (it can be read off the table).
#   E2 G1 holds with the candidate's overall Yamazaki RMSE between 0.3 and 0.7; root below
#      0.15 (kappa_d_rec was fitted to root on the 4-pool, and the root barely depends on the
#      shoot); straw below 1.0; PFBA and PFPeA grain over-predicted by more than 3x.
#   E3 G2 and G3 hold (one PFOA run took 2.3 s).
#   E4 the band becomes median(0.640, 0.47, 0.36, 0.51) = 0.49, x3.1.
#   E5 sensitivity: Km consequential in every tissue; stem_transp_frac and retention
#      consequential for straw; lam_grain consequential for grain; E_m consequential.
#   E6 the inverse field check (Kim 2019 PFOA brown rice 0.349 ug/kg against void water
#      0.0787 ug/L; 16x over under the current default) moves to within 2-5x over, because
#      the candidate's PFOA grain BAF is ~5x the current default's.
#
#   python validation/default_switch.py
#
# -----------------------------------------------------------------------------
# OUTCOME (appended after the run; the rule above is unchanged)
#   D1 PASS on all five endpoints (table above; recomputed: 1.714/0.992/0.640,
#      2.329/0.563/0.469, 2.052/0.376/0.154, 1.260/0.896/0.358, 1.682/1.241/0.510).
#   G4 PASS: simulate(shoot="redistributed", lipid_loading=True) is bit-identical to the
#      scored simulate_nstem_leaf runs on all 13 congeners (worst |log10| 0).
#   G1 PASS: Yamazaki 0.380 (root 0.090, straw 0.253, grain 0.619) vs Fit B 1.198.
#   G2 PASS: 13/13 congeners finite and positive; worst mass residual 5e-16.
#   G3 PASS: 1.71 s per run.
#   => SWITCH. DEFAULT_SHOOT = "redistributed", DEFAULT_LIPID_LOADING = True.
#   E1 held. E2 held on its stated terms (0.380; root 0.090; straw 0.253; PFBA 6.6x and PFPeA
#      7.1x over), BUT NOT FORESEEN: the candidate over-predicts Yamazaki grain for nearly every
#      congener (PFHpA 6.9x, PFUnDA 7.3x, PFOA 3.2x, PFNA 2.9x, PFOS 2.3x, PFDA 1.9x; PFDoDA 1.1x),
#      so its in-sample grain error (0.619) is larger than its transfer grain errors (Kim 0.47,
#      Li 0.51). E3 held. E4 held (x3.1). E5 partly wrong: Km (1.27), E_m (0.13), retention
#      (straw 0.11, grain 0.23) and lam_grain (grain 0.063) are consequential, but stem_transp_frac
#      is NOT (straw 0.011 against a 0.036 yardstick); dead blades lost moves straw 0.045; phi and
#      T_C_Ph are not consequential for this shoot. E6 held: the Kim PFOA inverse is 2.7x over
#      (0.214 ug/L, 95% 0.106-0.433) against 16.4x for the old default.
# =============================================================================
"""Pre-registered decision on the default configuration (see the header)."""
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "src"))
sys.path.insert(0, HERE)
import model_api as api            # noqa: E402
import site_scenario as ss         # noqa: E402

CAND = dict(shoot="redistributed", lipid_loading=True)
CURRENT = dict(api.FOUR_POOL_FREE, f_xy_source="recommended")
TISS = ("root", "straw", "grain")
OBS = api.yamazaki_obs()
YCONG = [c for c in api.CONGENERS if c in OBS]
FIT_B = 1.198                       # reproduce_demo.py --rec (current default, used cells)
MARGIN_DEFAULT, MARGIN_LIPID = 0.3, 0.1
SENS = {                            # name -> settings (pre-registered ranges)
    "E_m": [dict(E_m_mV=-90.0)],
    "Km": [dict(km_scale=100.0)],
    "phi": [dict(phi=0.05), dict(phi=0.2)],
    "T_C_Ph": [dict(T_C_Ph=4.0), dict(T_C_Ph=12.0)],
    "dead_blades_lost": [dict(dead_leaf_attached=0.0)],
    "stem_transp_frac": [dict(stem_transp_frac=0.30), dict(stem_transp_frac=0.60)],
    "retention": [dict(retention=0.40), dict(retention=0.80)],
    "lam_grain": [dict(lam_grain=0.025), dict(lam_grain=0.18)],
}


# ---------------------------------------------------------------------------- D1
def transfer_scores():
    import oos_tang_lipid as T
    import kim2019_censored as K
    import li2025_stage_resolved as L
    tang = {n: T._rmse("low", **kw)[0] for n, kw in (
        ("current", dict(shoot="4pool", lipid_loading=False)),
        ("lipid4", dict(shoot="4pool", lipid_loading=True)),
        ("candidate", dict(shoot="redistributed", lipid_loading=True)))}
    kim = K.main(verbose=False) if "verbose" in K.main.__code__.co_varnames else K.main()
    li = L.main(sensitivity=False) if "sensitivity" in L.main.__code__.co_varnames else L.main()
    ksc, lsc = kim["scores"], li["scores"]
    rows = {
        "Tang 2026 per-organ TF": (tang["current"], tang["lipid4"], tang["candidate"]),
        "Kim 2019 grain, point": (ksc["mono"]["point"], ksc["lipid"]["point"], ksc["nstem_lipid"]["point"]),
        "Kim 2019 grain, interval": (ksc["mono"]["interval"], ksc["lipid"]["interval"],
                                     ksc["nstem_lipid"]["interval"]),
        "Li 2025 shoot TF": (lsc["mono"]["rmse"]["shoot"], lsc["lipid"]["rmse"]["shoot"],
                             lsc["nstem_lipid"]["rmse"]["shoot"]),
        "Li 2025 grain TF": (lsc["mono"]["rmse"]["grain"], lsc["lipid"]["rmse"]["grain"],
                             lsc["nstem_lipid"]["rmse"]["grain"]),
    }
    return rows


def d1(rows):
    ok = {k: (cur - cand > MARGIN_DEFAULT) and (cand - lip <= MARGIN_LIPID)
          for k, (cur, lip, cand) in rows.items()}
    return all(ok.values()), ok


# ---------------------------------------------------------------------------- G4
def g4():
    """simulate()'s redistributed path must BE the scored configuration: compare it with
    simulate_nstem_leaf (what the transfer scripts ran) on every curated congener."""
    worst = 0.0
    for c in api.CONGENERS:
        a = api.simulate(c, **CAND)["baf_final_dw"]
        b = api.simulate_nstem_leaf(c, lipid_loading=True)["baf_final_dw"]
        worst = max(worst, max(abs(np.log10(a[k] / b[k])) for k in a))
    return worst


# ---------------------------------------------------------------------------- G1 / sensitivity
def _per_run(args):
    c, kw = args
    cwo = OBS[c]["water_ng_L"] / 1000.0
    out = []
    for h in ss.YAMAZAKI_HARVESTS:
        drv = api.scenario_drivers(Cwo=cwo, harvest=h)
        r = api.simulate(c, drivers=drv, **kw)
        out.append({k: r["baf_final_dw"][k] for k in TISS})
    return c, out


def observable(pool, kw):
    runs = dict(pool.map(_per_run, [(c, kw) for c in YCONG]))
    geo = {c: {k: float(10 ** np.mean([np.log10(r[k]) for r in runs[c]])) for k in TISS}
           for c in YCONG}
    return runs, geo


def residuals(geo):
    res = {}
    for c in YCONG:
        d = OBS[c]
        res[c] = {k: float(np.log10(geo[c][k] / d["baf_dw"][k])) for k in TISS if d["use"][k]}
    return res


def rmse(res, tissues=TISS):
    e = [v[k] for v in res.values() for k in v if k in tissues]
    return float(np.sqrt(np.mean(np.square(e))))


def sensitivity(pool, base_runs, base_geo):
    ref = {k: float(np.median([abs(np.log10(base_runs[c][0][k] / base_runs[c][1][k]))
                               for c in YCONG])) for k in TISS}
    eff = {}
    for name, settings in SENS.items():
        per = []
        for st in settings:
            _, geo = observable(pool, {**CAND, **st})
            per.append({c: {k: float(np.log10(geo[c][k] / base_geo[c][k])) for k in TISS}
                        for c in YCONG})
        med = {k: float(np.median([max(abs(p[c][k]) for p in per) for c in YCONG])) for k in TISS}
        eff[name] = dict(median_abs_delta=med, consequential={k: med[k] > ref[k] for k in TISS})
    return ref, eff


# ---------------------------------------------------------------------------- G2 / G3
def g2():
    bad = []
    worst_mass = 0.0
    from pfas_rice_plant_module_4pool_surf import binding_factors, root_uptake
    drv = api.scenario_drivers(Cwo=1.0)
    for c in api.CONGENERS:
        r = api.simulate(c, **CAND)
        vals = [v for v in r["baf_final_dw"].values()] + [v for v in r["baf_final"].values()]
        if not (r["success"] and np.all(np.isfinite(vals)) and min(vals) > 0):
            bad.append(c)
        cmpd = api._compound_for(c, "recommended", lipid_loading=True, K_cw_organ="stem")
        m = api._run_redistributed(c, cmpd, drv["t"], drv["Cwo"], drv["Qtp"], drv["M"],
                                   drv["leaf_loss"], float(drv["t"][-1]), -120.0, _keep_model=True)
        mod, sol = m["_model"], m["_sol"]
        B = binding_factors(mod.comps, mod.cmpd)
        for ti in (40.0, 80.0, 110.0):
            C = sol.sol(ti); dC = mod.rhs(ti, C)
            M = mod.inputs.M_(ti); dM = mod.inputs.dM_(ti)
            jR = root_uptake(mod.inputs.Cwo_(ti), (C / B)[0], mod.cmpd, mod.env)
            dmass = float(np.sum(dM * C + M * dC))
            src = float(M[0] * jR - mod.inputs.leaf_loss_(ti) * M[mod.LEAF] * C[mod.LEAF])
            worst_mass = max(worst_mass, abs(dmass - src) / max(abs(src), 1e-30))
    return bad, worst_mass


def g3():
    t0 = time.time()
    api.simulate("PFOA", **CAND)
    return time.time() - t0


# ---------------------------------------------------------------------------- E6
def inverse_check():
    truth = 0.0787
    cur = api.estimate_exposure_bayesian("PFOA", {"grain": 0.349}, **CURRENT)
    cand = api.estimate_exposure_bayesian("PFOA", {"grain": 0.349}, **CAND)
    return truth, cur["median"], cand["median"], cand["ci95"]


def main():
    print("=" * 88)
    print("DEFAULT SWITCH -- pre-registered (see header)")
    print("=" * 88)
    rows = transfer_scores()
    ok1, per = d1(rows)
    print("\nD1  transfer scores (log10 RMSE)        current   4-pool lipid   candidate   pass")
    for k, (cur, lip, cand) in rows.items():
        print(f"    {k:34s} {cur:8.3f} {lip:12.3f} {cand:11.3f}   {per[k]}")
    print(f"    D1: {'PASS' if ok1 else 'FAIL'}")

    w4 = g4()
    print(f"\nG4  simulate(shoot='redistributed', lipid_loading=True) vs the scored runs: "
          f"worst |log10| difference {w4:.2e} -> {'PASS' if w4 < 0.005 else 'FAIL'}")

    with ProcessPoolExecutor(max_workers=os.cpu_count() or 2) as pool:
        base_runs, base_geo = observable(pool, CAND)
        res = residuals(base_geo)
        overall = rmse(res)
        per_t = {k: rmse(res, (k,)) for k in TISS}
        ok_g1 = overall <= FIT_B
        print(f"\nG1  Yamazaki two-harvest observable, used cells: candidate {overall:.3f} "
              f"(root {per_t['root']:.3f}, straw {per_t['straw']:.3f}, grain {per_t['grain']:.3f}) "
              f"vs current Fit B {FIT_B:.3f} -> {'PASS' if ok_g1 else 'FAIL'}")
        print(f"    {'congener':8s}" + "".join(f"{k:>18s}" for k in TISS))
        for c in YCONG:
            d = OBS[c]
            print(f"    {c:8s}" + "".join(
                (f"{base_geo[c][k]:9.3g}/{d['baf_dw'][k]:<8.3g}" if d["use"][k] else f"{base_geo[c][k]:9.3g}/{'--':<8s}")
                for k in TISS))
        bad, mass = g2()
        print(f"\nG2  13 congeners: failures {bad or 'none'}; worst mass residual {mass:.1e} "
              f"-> {'PASS' if not bad and mass < 1e-6 else 'FAIL'}")
        tt = g3()
        print(f"G3  one PFOA run {tt:.2f} s -> {'PASS' if tt < 5.0 else 'FAIL'}")
        switch = ok1 and ok_g1 and (not bad and mass < 1e-6) and tt < 5.0 and w4 < 0.005
        print(f"\nVERDICT: {'SWITCH' if switch else 'KEEP THE CURRENT DEFAULT'}")

        band = [rows["Tang 2026 per-organ TF"][2], rows["Kim 2019 grain, point"][2],
                rows["Li 2025 shoot TF"][2], rows["Li 2025 grain TF"][2]]
        med = float(np.median(band))
        print(f"\nE4  band from the candidate's transfer errors {[round(b, 3) for b in band]}: "
              f"median {med:.3f} -> x{10 ** med:.1f}")
        truth, cur_m, cand_m, ci = inverse_check()
        print(f"E6  inverse, Kim PFOA brown rice 0.349 ug/kg (void water {truth} ug/L): current "
              f"{cur_m:.3f} ({cur_m / truth:.1f}x), candidate {cand_m:.3f} ({cand_m / truth:.1f}x; "
              f"95% {ci[0]:.3f}-{ci[1]:.3f})")

        ref, eff = sensitivity(pool, base_runs, base_geo)
        print("\nSENSITIVITY of the candidate (median |log10| change; yardstick = two harvest runs "
              f"root {ref['root']:.3f} / straw {ref['straw']:.3f} / grain {ref['grain']:.3f})")
        for n, e in eff.items():
            m = e["median_abs_delta"]
            print(f"    {n:18s} root {m['root']:.3f}  straw {m['straw']:.3f}  grain {m['grain']:.3f}"
                  f"   consequential: {[k for k, v in e['consequential'].items() if v] or 'none'}")
    return dict(rows=rows, d1=ok1, g4=w4, g1=overall, g1_tissue=per_t, residuals=res,
                g2=(bad, mass), g3=tt, switch=switch, band=band, band_median=med,
                inverse=(truth, cur_m, cand_m, ci), sens=(ref, eff))


if __name__ == "__main__":
    main()
