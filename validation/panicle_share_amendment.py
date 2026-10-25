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
#   python validation/panicle_share_amendment.py
# =============================================================================
"""Pre-registered amendment: panicle-size-scaled grain share (see the header)."""
raise SystemExit("pre-registration only: the amendment is implemented and run after this "
                 "header is committed")
