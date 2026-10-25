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
# =============================================================================
"""Pre-registered decision on the default configuration (see the header)."""

if __name__ == "__main__":
    raise SystemExit("pre-registration only; the analysis is added after this commit is pushed")
