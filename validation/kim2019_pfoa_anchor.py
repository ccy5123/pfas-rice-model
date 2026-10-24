#!/usr/bin/env python3
# =============================================================================
# validation/kim2019_pfoa_anchor.py
# -----------------------------------------------------------------------------
# KIM 2019 PFOA GRAIN BAF: WHICH ESTIMATOR?  (handoff step 4b, docs/HANDOFF_mass_basis_refit.md)
#
# PRE-REGISTRATION -- written and pushed before the censored estimate, the slope, the
# leave-one-out checks and the model-side numbers below were produced.
#
# Why: the repo's Kim PFOA grain BAF (4.43 L/kg) is a ratio of two field AVERAGES, and
# one WWTP-impacted field (A-2, 1760 ng/L) carries 82.8% of the void-water PFOA sum.
# Main-text Table 1 (p.718) gives PFOA per site, so the BAF can be estimated per site.
# Data: docs/literature_db/raw_si/kim2019_table1_sites.csv (transcribed by the user from
# the PDF text layer, checked row by row against the rendered page). Reproduces SI Table
# S6: rice PFOA mean 0.3494 ng/g with 17/30 detects (SI 0.349, 57%), void water 78.73 ng/L
# with 26/27 (SI 78.7, 96%), soil 0.1599 ng/g dw with 13/30 (SI 0.160, 43%).
#
# Already computed (data only, before this file):
#   E1  published ratio of means (rice over 30 sites / water over 27)   4.44 L/kg
#   E2  ratio of means over the 27 sites that have water (ND = 0)      4.18
#   E3  ratio of means without A-2                                     21.5
#   E4  per-site geometric mean, the 14 sites with BOTH detected       27.0 (median 28.3,
#       range 0.98 at A-2 to 597 at D-7)
# Also already seen: kim2019_censored.py (step 4a) PRINTED the four variants' PFOA grain
# BAF as unscored anchor rows, all far below 4.43. Nothing below uses them to choose an
# estimator. And: no shipped parameter is fitted to Kim PFOA (L_Ph_W2fit, L_Ph_oryza and
# LIPID_LOADING are Yamazaki fits; only the literature_params DEMO fits L_Ph to 4.43), so
# the "PFOA* used in the L_Ph fit" exclusion in the OOS scripts has no basis.
#
# Site accounting: 30 sites. 14 have rice and water both detected; 12 have water detected
# and rice < MDL (0.088 ng/g) -> per-site BAF < 0.088/water, a LEFT-CENSORED ratio; C-6
# has both < MDL (no information on the ratio); E-2..E-4 have no water (NS).
#
# Questions and decision rules:
#   Q1 (data) Which estimator? Admissibility, stated before E5 is computed:
#        (a) numerator and denominator come from the same sites;
#        (b) no site is dropped for its detection status or its value;
#        (c) no single site carries more than half the weight.
#      E1 fails (a) and (c); E2 fails (c); E3 fails (b) (A-2 dropped for its value; the
#      WWTP footnote also marks D-1, which E3 keeps); E4 fails (b) (it drops the 12 rice-ND
#      sites, which are low-ratio sites by construction, so it is biased UP).
#      E5 = censored per-site geometric mean: ML fit of a normal to log10(rice/water) over
#      the 26 sites with both a water value and an informative ratio (14 exact, 12
#      left-censored at log10(MDL_rice/water)). E5 is the only candidate that passes all
#      three; the criteria were written knowing that, so the content of the rule is the
#      criteria, stated so they can be disputed. Uncertainty: bootstrap over the 26 sites
#      (2000 resamples) -> 90% interval. Sensitivity (not primary): the Kaplan-Meier median
#      on the flipped scale (the standard estimator for multiple detection limits), and the
#      leave-one-site-out range of E5.
#   Q2 (data) Does the per-site ratio behave like a BAF? A constant BAF means rice is
#      proportional to water: slope 1 of log10 rice on log10 water. Censored (Tobit)
#      regression over the same 26 sites (E5 is its b = 1 special case); bootstrap 90%
#      interval of b. If the interval excludes 1, the harvest void water (a 2-day snapshot)
#      is not proportional to what the plant took up over the season, and E5 is a typical
#      value with that noise in it, not a plant constant.
#   Q3 (model) Each variant's PFOA grain BAF (kim2019_censored.predict: mono, W2, lipid,
#      nstem_lipid; Cwo 1, fresh weight) against E1 and against E5 with its interval.
#   Q3b (model) Is the model linear in exposure over Kim's water range? Slope of log10
#      grain conc on log10 Cwo between 0.0008 and 1.76 ug/L (the lowest and highest PFOA
#      void water), per variant. If ~1, a data slope below 1 cannot come from the carrier.
#   Q4 (data) A-2's share of the void-water SUM-17 (Table 1). Kim's other congeners have
#      no per-site values, so this is the only check on whether their denominators are
#      A-2-weighted too.
#
#   DECISION RULES:
#   D1 kim2019_grain_baf() keeps the published ratio of means for ALL congeners: it is the
#      only estimator the other five have, and a chain-length series on mixed estimators
#      would compare unlike quantities. The per-site PFOA analysis is added beside it, not
#      substituted.
#   D2 The literature_params demo's L_Ph target moves from E1 to E5 if E5's 90% interval
#      excludes E1; otherwise it stays at E1.
#   D3 The "PFOA* used in the fit" labels are removed (factually wrong). PFOA is reported as
#      a separate OOS check (each variant vs E5 and its interval: inside / below / above);
#      it is NOT folded into the pre-registered 5-congener Kim RMSE of step 4a, which stays
#      the record.
#   D4 If Q4's A-2 share exceeds 50%, the absolute LEVEL of every Kim congener is recorded as
#      A-2-weighted (the chain-length SHAPE compares like with like; the level does not).
#
#   EXPECTATION (stated so it can be wrong): E5 lies below E4 (the dropped sites are
#      low-ratio sites), somewhere around 10-20, and its interval excludes E1 (so D2 moves
#      the demo target); b is well below 1 (rice varies ~19x across sites, water ~2000x);
#      the model is close to linear (Q3b slope > 0.9), so the data slope is not the carrier;
#      every variant sits below E5's interval. No expectation for Q4.
#
#   python validation/kim2019_pfoa_anchor.py
# =============================================================================
