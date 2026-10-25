#!/usr/bin/env python3
# =============================================================================
# validation/panicle_share_amendment.py
# -----------------------------------------------------------------------------
# SHOULD THE REDISTRIBUTED SHOOT GIVE THE PANICLE A SHARE OF THE XYLEM STREAM THAT
# SCALES WITH ITS SIZE?
#
# PRE-REGISTRATION (committed and pushed before the amendment is implemented or this
# file's analysis is run; the results are appended below the line at the end).
#
# THE DEFECT
#   In the default shoot (NStemLeafModel) the grain receives, once the panicle has formed,
#       residual = (lam_grain + (1 - retention) * (sum(tau) + lam_leaf)) * xylem_flux
#   i.e. the panicle's own transpiration deposit (lam_grain 0.05) PLUS every non-retained
#   share of the stem and leaf deposits -- 0.05 + 0.4 * 0.95 = 0.43 of the root's xylem
#   export, from the day the grain mass leaves its floor (the formation gate gamma ramps
#   0 -> 1 between the floor and 1.5x the floor). Before that the residual goes to the leaf.
#   The panicle is tiny when the gate opens, so its concentration peaks ~2 weeks before
#   heading (PFOA ~3.3x its harvest value, pinned in tests/test_model_api.py::
#   test_grain_formation_gate) and then dilutes. The 4-pool core has no such peak: its grain
#   share follows the grain's dry matter. Scaling lam_grain alone cannot remove the peak
#   (it is 0.05 of the 0.43); the routing of the whole residual has to scale.
#
# THE AMENDMENT (one variant; decided now)
#   grain_share="size": the grain receives w(t) * residual and the leaf (1 - w(t)) * residual,
#   with
#       w(t) = clip( (Mg(t) - Mg_min) / (Mg_max - Mg_min), 0, 1 )
#   where Mg is the grain column of the driver's organ-mass series and Mg_min / Mg_max are
#   that series' minimum and maximum (the same series statistics the formation gate already
#   uses). Grain theta is constant, so the fresh ratio equals the dry-matter ratio. A driver
#   with no floor period (constant mass) gets w = 1, exactly as it gets gamma = 1. The phloem
#   stream to the grain is NOT changed: it is already proportional to the grain's dry growth
#   and stays gated by gamma. No new parameter.
#   Why this form and not the 4-pool's surface split: w = 1 at the crop's mature panicle
#   size, so the harvest-time routing of the default is unchanged and the amendment only
#   moves delivery away from the small panicle -- the defect named above. A surface split
#   (S_grain*Mg / (S_leaf*Ml + S_grain*Mg)) would cut the grain's share at maturity from
#   0.43 to ~0.06 of the xylem export: that is a different model, not this amendment.
#   Mg_max comes from the driver's whole season. Drivers are computed before the PFAS solve,
#   so this is available; it is a normalisation by the crop's own mature panicle mass.
#   grain_share="gate" (the current behaviour) stays runnable.
#
# DECISION RULE (current = the default configuration as shipped: redistributed shoot,
# lipid loading, carrier, grain_share "gate")
#   A1 the defect: with "size", for every one of the 13 curated congeners the default run's
#      grain concentration never exceeds 1.5x its harvest value (the criterion the 4-pool
#      already meets in test_grain_formation_gate).
#   A2 transfer, non-inferiority: on each of the five endpoints the default was chosen on
#      -- Tang 2026 per-organ TF (0.1 ug/g), Kim 2019 grain point and [0, MDL] interval,
#      Li 2025 shoot and grain TF by stage -- "size" is not worse than "gate" by more than
#      0.10 log10 RMSE. These datasets already chose the default, so this is a guard against
#      making the transfer worse, not evidence for the amendment; an improvement is reported
#      but cannot count in its favour.
#   A3 Yamazaki guard: the two-harvest dry-basis observable over the used cells is not worse
#      than "gate" (0.380) by more than 0.10 (in-sample for the lipid constants).
#   A4 all 13 congeners solve with finite positive BAFs and the mass balance closes (the
#      default_switch G2 invariant, relative residual < 1e-6).
#   ADOPT ("size" becomes the default) iff A1-A4 hold. Otherwise the default stays "gate",
#   "size" stays runnable as an option, and the failed item is recorded.
#   A5 (after an adoption) simulate() must reproduce the scored "size" runs to 0.005 log.
#
# WHAT AN ADOPTION CHANGES (decided now, so the run cannot steer it)
#   * model_api.DEFAULT_GRAIN_SHARE = "size"; simulate()/simulate_nstem_leaf()/
#     tang_tf_validation() take grain_share=; NStemLeafModel(grain_share=...).
#   * APRIORI_OOS_SCORES and the honesty band become the "size" scores on the same four
#     entries (Tang, Kim point, Li shoot, Li grain); the selection caveat stays.
#   * Every number the docs, app text and policy briefs quote from the default is
#     regenerated (CLAUDE.md, ASSUMPTIONS, HANDOFF, README, POLICY_BRIEF_KR via
#     policy_brief_runs.py); the young-panicle test changes from "peak present" to "no peak".
#   * Not re-fitted: nothing. The lipid constants and *_rec were fitted on the 4-pool shoot,
#     which this amendment does not touch.
#
# DESCRIPTIVE ONLY (cannot change the decision): the same table for a heading-referenced
#   size, w = min(1, Mg(t) / Mg(heading)) -- the "it transpires little inside the boot"
#   reading -- to show how much the result depends on the reference chosen.
#
# EXPECTATIONS (written before the run)
#   E1 A1 holds for all 13 congeners; no grain maximum before heading.
#   E2 harvest grain BAF falls for every congener, by 1.5-4x (Yamazaki observable PFOA grain
#      1.71 -> 0.4-1.1 L/kg); root changes < 0.01 log; straw rises by less than 0.1 log
#      (the residual that no longer reaches the small panicle goes to the leaf).
#   E3 Yamazaki improves (grain was over-predicted 2-7x): grain RMSE 0.619 -> 0.25-0.45,
#      overall 0.380 -> 0.25-0.35. A3 holds.
#   E4 transfer: Tang worse by 0.05-0.25 (the endosperm is already under-predicted);
#      Kim point within +-0.10; Kim interval worse, 0.15 -> 0.25-0.40 (PFUnDA and PFNA move
#      further below their intervals); Li grain TF changes by up to +-0.15, direction not
#      predicted; Li shoot changes by less than 0.05.
#   E5 net: A1, A3, A4 pass and A2 FAILS (Kim interval and/or Tang) -> the default stays
#      "gate" and "size" is kept as an option.
#   E6 PFOA grain influx through the xylem falls from 80 % to 55-70 % (apportionment).
#   E7 A4 passes (the routing only redistributes the same xylem export between leaf and
#      grain, so the balance still telescopes).
#
#   python validation/panicle_share_amendment.py   (~5 min)
#
# -----------------------------------------------------------------------------
# OUTCOME (appended after the run; the rule above is unchanged)
#   A1 PASS: with "size" no congener's grain exceeds its harvest value (max/harvest 1.00 for
#      all 13; "gate" 3.21-3.82, maximum at day 73-74, heading day 87.0).
#   A2 FAIL on two endpoints (gate -> size, log10 RMSE): Tang 0.640 -> 0.708 (+0.069, ok);
#      Kim point 0.469 -> 0.421 (-0.048, ok); Kim interval 0.154 -> 0.267 (+0.113, FAIL);
#      Li shoot 0.358 -> 0.398 (+0.040, ok); Li grain 0.510 -> 0.940 (+0.430, FAIL).
#   A3 PASS: Yamazaki 0.380 -> 0.254 (root 0.090 -> 0.088, straw 0.253 -> 0.233, grain
#      0.619 -> 0.372).
#   A4 PASS: 13/13 finite and positive; worst mass residual 8e-16.
#   => KEEP grain_share="gate" as the default; "size" stays runnable
#      (simulate(grain_share="size"), simulate_nstem_leaf, tang_tf_validation, apportionment).
#   Expectations: E1 held. E2 mostly held -- grain fell for every congener, by 1.48-2.29x
#      (-0.171..-0.360 log; PFOA 1.48x, just under the stated 1.5x), root moved < 0.01 log
#      (max +0.008), but straw rose by up to 0.104 log (stated < 0.1) and the Yamazaki PFOA
#      grain landed at 1.16 L/kg (stated 0.4-1.1). E3 held (0.254; grain 0.372). E4: Tang
#      +0.069, Kim point -0.048, Kim interval 0.267 and Li shoot +0.040 all inside their stated
#      ranges, but Li grain WRONG -- +0.430 against "up to +-0.15". E5 held in substance: A2 failed
#      (on Kim interval, as expected, and on Li grain, which was not), so the default stays.
#      E6 slightly wrong: PFOA grain xylem share 80.3 % -> 71.2 % (stated 55-70 %). E7 held.
#   NOT FORESEEN, and the reason the default stays: Li 2025 scores the grain from the filling
#      stage (model day ~78), when the "gate" panicle is still near its early maximum; with
#      "size" the young grain holds much less, and the grain/root TF at filling-milk falls
#      well below Li's. Li's field data therefore favour EARLY loading of the panicle over the
#      size-scaled share -- the young-panicle maximum may be partly real, not only an artefact
#      of the routing (the mechanism is not identified; Li's grain is not defined in the paper).
#      Yamazaki, which sees only the harvest, prefers "size". The two datasets pull in opposite
#      directions on the timing of grain loading, and harvest-only data cannot decide it.
#   DESCRIPTIVE heading-referenced size: harvest grain BAF (FRESH, Cwo = 1) between the two
#      (PFOA 1.11 vs gate 1.24 / size 0.83 fresh); max/harvest <= 1.10 for all congeners.
# =============================================================================
"""Pre-registered amendment: panicle-size-scaled grain share (see the header)."""
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "src"))
sys.path.insert(0, HERE)
import model_api as api            # noqa: E402
import default_switch as DS        # noqa: E402  (the Yamazaki observable, used cells)

SHARES = ("gate", "size")
BASE = dict(shoot="redistributed", lipid_loading=True)
MARGIN_TRANSFER, MARGIN_YAMAZAKI, PEAK_MAX = 0.10, 0.10, 1.5
ENDPOINTS = ("Tang 2026 per-organ TF", "Kim 2019 grain, point", "Kim 2019 grain, interval",
             "Li 2025 shoot TF", "Li 2025 grain TF")


def t_heading():
    run = api.site_run()
    return float(np.interp(1.0, run["dvs"], run["t"]))


# ---------------------------------------------------------------------------- A1
def peak_table(share):
    """max grain concentration / harvest value, and the day of the maximum, per congener."""
    out = {}
    for c in api.CONGENERS:
        r = api.simulate(c, grain_share=share, **BASE)
        cg = np.asarray(r["conc"]["grain"], float)
        i = int(np.argmax(cg))
        out[c] = dict(ratio=float(cg[i] / cg[-1]), t_max=float(r["t"][i]),
                      grain=float(r["baf_final_dw"]["grain"]), root=float(r["baf_final_dw"]["root"]),
                      straw=float(r["baf_final_dw"]["straw"]))
    return out


# ---------------------------------------------------------------------------- A2
def transfer_scores(share):
    import oos_tang_lipid as T
    import kim2019_censored as K
    import li2025_stage_resolved as L
    import literature_params as lp
    tang = T._rmse("low", grain_share=share, **BASE)[0]
    b = lp.kim2019_grain_baf_bounds()
    pred = {nm: float(api.simulate_nstem_leaf(nm, Cwo=1.0, lipid_loading=True,
                                              grain_share=share)["baf_final"]["grain"])
            for nm in K.SCORED}
    ep = [K.point_error(pred[nm], b[nm]) for nm in K.SCORED]
    ei = [K.interval_error(pred[nm], b[nm]) for nm in K.SCORED]
    data = L.load_data()
    days = L.stage_days()
    model = {c: L.model_tf(api.simulate_nstem_leaf(c, Cwo=1.0, lipid_loading=True,
                                                   grain_share=share), days)
             for c in L.CONGENERS}
    li = L.score(data, model, np.random.default_rng(0))
    return {ENDPOINTS[0]: tang,
            ENDPOINTS[1]: float(np.sqrt(np.mean(np.square(ep)))),
            ENDPOINTS[2]: float(np.sqrt(np.mean(np.square(ei)))),
            ENDPOINTS[3]: li["rmse"]["shoot"], ENDPOINTS[4]: li["rmse"]["grain"]}, pred


# ---------------------------------------------------------------------------- A4
def mass_check(share):
    """default_switch G2 on the given share: finite positive BAFs on 13 congeners and the
    mass balance (d/dt of the total burden == root uptake - leaf senescence loss)."""
    from pfas_rice_plant_module_4pool_surf import binding_factors, root_uptake
    bad, worst = [], 0.0
    drv = api.scenario_drivers(Cwo=1.0)
    for c in api.CONGENERS:
        r = api.simulate(c, grain_share=share, **BASE)
        vals = list(r["baf_final_dw"].values()) + list(r["baf_final"].values())
        if not (r["success"] and np.all(np.isfinite(vals)) and min(vals) > 0):
            bad.append(c)
        cmpd = api._compound_for(c, "recommended", lipid_loading=True, K_cw_organ="stem")
        m = api._run_redistributed(c, cmpd, drv["t"], drv["Cwo"], drv["Qtp"], drv["M"],
                                   drv["leaf_loss"], float(drv["t"][-1]), -120.0,
                                   grain_share=share, _keep_model=True)
        mod, sol = m["_model"], m["_sol"]
        B = binding_factors(mod.comps, mod.cmpd)
        for ti in (40.0, 60.0, 80.0, 110.0):
            C = sol.sol(ti); dC = mod.rhs(ti, C)
            M = mod.inputs.M_(ti); dM = mod.inputs.dM_(ti)
            jR = root_uptake(mod.inputs.Cwo_(ti), (C / B)[0], mod.cmpd, mod.env)
            dmass = float(np.sum(dM * C + M * dC))
            src = float(M[0] * jR - mod.inputs.leaf_loss_(ti) * M[mod.LEAF] * C[mod.LEAF])
            worst = max(worst, abs(dmass - src) / max(abs(src), 1e-30))
    return bad, worst


# ---------------------------------------------------------------------------- descriptive
def heading_referenced(t_head):
    """DESCRIPTIVE ONLY: w = min(1, Mg(t)/Mg(heading)) -- the 'transpires little inside the
    boot' reading. Built by re-pointing the size weight of a "size" model; harvest grain BAF
    (fresh, at Cwo = 1) and max/harvest per congener."""
    from scipy.interpolate import interp1d
    drv = api.scenario_drivers(Cwo=1.0)
    out = {}
    for c in api.CONGENERS:
        cmpd = api._compound_for(c, "recommended", lipid_loading=True, K_cw_organ="stem")
        m = api._run_redistributed(c, cmpd, drv["t"], drv["Cwo"], drv["Qtp"], drv["M"],
                                   drv["leaf_loss"], float(drv["t"][-1]), -120.0,
                                   grain_share="size", _keep_model=True)
        mod = m["_model"]
        inp = mod.inputs
        gm = inp.M[:, -1]
        glo = float(gm.min())
        gh = float(np.interp(t_head, inp.t, gm))
        w = np.clip((gm - glo) / max(gh - glo, 1e-30), 0.0, 1.0)
        inp._grain_size = interp1d(inp.t, w, kind="linear", bounds_error=False,
                                   fill_value="extrapolate")
        sol = mod.solve(inp.t)
        cg = sol.y[mod.GRAIN]
        out[c] = dict(grain=float(cg[-1]), ratio=float(cg.max() / cg[-1]))
    return out


def _apportion(share):
    a = api.apportionment("PFOA", grain_share=share, **BASE)
    return a["fraction"]["grain"]


def main():
    print("=" * 88)
    print("PANICLE-SIZE-SCALED GRAIN SHARE -- pre-registered (see header)")
    print("=" * 88)
    th = t_heading()
    peaks = {s: peak_table(s) for s in SHARES}
    print(f"\nA1  grain max / harvest value (heading at day {th:.1f}); the rule needs <= {PEAK_MAX}")
    print(f"    {'congener':8s}{'gate':>16s}{'size':>16s}{'grain BAF dw gate -> size':>30s}")
    for c in api.CONGENERS:
        g, z = peaks["gate"][c], peaks["size"][c]
        print(f"    {c:8s}{g['ratio']:8.2f} @{g['t_max']:5.1f}{z['ratio']:8.2f} @{z['t_max']:5.1f}"
              f"{g['grain']:14.3g} -> {z['grain']:<8.3g}({np.log10(z['grain'] / g['grain']):+.3f} log)")
    a1 = all(peaks["size"][c]["ratio"] <= PEAK_MAX for c in api.CONGENERS)
    print(f"    A1: {'PASS' if a1 else 'FAIL'}")
    dr = {k: [np.log10(peaks['size'][c][k] / peaks['gate'][c][k]) for c in api.CONGENERS]
          for k in ("root", "straw", "grain")}
    print("    change, size vs gate (log10, over 13 congeners): "
          + "; ".join(f"{k} {min(v):+.3f}..{max(v):+.3f}" for k, v in dr.items()))

    sc = {s: transfer_scores(s) for s in SHARES}
    print("\nA2  transfer (log10 RMSE)                  gate     size    size-gate  pass")
    a2 = True
    for k in ENDPOINTS:
        g, z = sc["gate"][0][k], sc["size"][0][k]
        ok = (z - g) <= MARGIN_TRANSFER
        a2 &= ok
        print(f"    {k:36s}{g:8.3f}{z:9.3f}{z - g:+10.3f}   {ok}")
    print(f"    A2: {'PASS' if a2 else 'FAIL'}")
    print("    Kim grain BAF predictions (L/kg): "
          + ", ".join(f"{nm} {sc['gate'][1][nm]:.3g}->{sc['size'][1][nm]:.3g}" for nm in sc["gate"][1]))

    with ProcessPoolExecutor(max_workers=os.cpu_count() or 2) as pool:
        yam = {}
        for s in SHARES:
            _, geo = DS.observable(pool, dict(BASE, grain_share=s))
            res = DS.residuals(geo)
            yam[s] = (DS.rmse(res), {k: DS.rmse(res, (k,)) for k in DS.TISS}, geo)
    a3 = yam["size"][0] - yam["gate"][0] <= MARGIN_YAMAZAKI
    print(f"\nA3  Yamazaki observable (used cells): gate {yam['gate'][0]:.3f} "
          f"(root {yam['gate'][1]['root']:.3f}, straw {yam['gate'][1]['straw']:.3f}, grain "
          f"{yam['gate'][1]['grain']:.3f}) -> size {yam['size'][0]:.3f} (root {yam['size'][1]['root']:.3f}, "
          f"straw {yam['size'][1]['straw']:.3f}, grain {yam['size'][1]['grain']:.3f}) -> "
          f"{'PASS' if a3 else 'FAIL'}")
    print(f"    {'congener':8s}{'grain obs':>10s}{'gate':>9s}{'size':>9s}")
    for c in DS.YCONG:
        d = DS.OBS[c]
        print(f"    {c:8s}{d['baf_dw']['grain']:10.3g}{yam['gate'][2][c]['grain']:9.3g}"
              f"{yam['size'][2][c]['grain']:9.3g}" + ("" if d["use"]["grain"] else "   (unused)"))

    bad, worst = mass_check("size")
    a4 = (not bad) and worst < 1e-6
    print(f"\nA4  size: failures {bad or 'none'}; worst mass residual {worst:.1e} -> "
          f"{'PASS' if a4 else 'FAIL'}")

    adopt = a1 and a2 and a3 and a4
    print(f"\nVERDICT: {'ADOPT grain_share=\'size\' as the default' if adopt else 'KEEP grain_share=\'gate\' (size stays an option)'}")

    fg = {s: _apportion(s) for s in SHARES}
    print("\nE6  PFOA grain influx shares (apportionment): "
          + "; ".join(f"{s}: " + ", ".join(f"{k} {v:.1%}" for k, v in fg[s].items()) for s in SHARES))

    hr = heading_referenced(th)
    print("\nDESCRIPTIVE (cannot change the decision): heading-referenced size, w = min(1, Mg/Mg(heading))")
    print(f"    {'congener':8s}{'grain BAF fw, Cwo=1':>22s}{'max/harvest':>13s}")
    for c in api.CONGENERS:
        print(f"    {c:8s}{hr[c]['grain']:22.3g}{hr[c]['ratio']:13.2f}")
    return dict(peaks=peaks, a1=a1, transfer=sc, a2=a2, yamazaki=yam, a3=a3, a4=(bad, worst),
                adopt=adopt, apportion=fg, heading=hr)


if __name__ == "__main__":
    main()
