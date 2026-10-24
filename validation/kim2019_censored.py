#!/usr/bin/env python3
# =============================================================================
# validation/kim2019_censored.py
# -----------------------------------------------------------------------------
# KIM 2019 GRAIN BAF WITH ITS NON-DETECTS  (handoff step 4, docs/HANDOFF_mass_basis_refit.md)
#
# PRE-REGISTRATION -- written and pushed before any model number below was produced.
#
# Why: every Kim 2019 grain BAF in the repo is a ratio of two field averages in which the
# non-detects were averaged as ZERO (stated in the paper). For the low-detection
# congeners both the numerator (brown rice, n=30) and the denominator (void water, n=27)
# are mostly non-detects -- PFDoDA rice 1/30, water 3/27 -- so the "measured" BAF is one
# corner of a range set by how the non-detects are treated, not a measurement.
#
# Data: SI Tables S5 (MDL) + S6 (min/avg/max/DF) -> docs/literature_db/raw_si/
# kim2019_si_detection.csv (averages identical to kim2019_field_conc.csv). Each non-detect
# lies in [0, MDL], so the true ratio of means lies in
#     [rice_avg / water_hi,  rice_hi / water_avg],   x_hi = x_avg + n_ND * MDL / n.
# Already computed (data only, before this file): PFHpA [0.39, 2.11], PFOA [4.43, 4.92],
# PFNA [2.21, 7.07], PFDA [1.55, 11.3], PFUnDA [19.8, 168], PFDoDA [5.5, 630] L/kg. The
# shipped ND=0 values are the LOWER bound for PFHpA/PFNA/PFDA/PFOA. Also already seen:
# the "reliable DF>=15%" Kim score in oos_crossdataset is ONE congener (PFNA) -- the
# label there says PFHpA+PFNA, but PFHpA's DF is 13%.
#
# Model: grain BAF at harvest, Cwo constant 1, default ORYZA2000 biomass + measured Q_TP,
# fresh organ mass (as oos_crossdataset). Variants: 4-pool monotone f_xy, 4-pool W2,
# 4-pool lipid loading (the three oos_crossdataset scores), plus nstem_leaf + lipid (the
# configuration that won Tang and Li). Scored congeners: PFHpA, PFNA, PFDA, PFUnDA, PFDoDA
# (PFOA stays out, as in oos_crossdataset). Model fresh-weight BAF vs Kim as published
# (brown-rice basis not stated; dw would move Kim by 0.065 log).
#
# Questions and decision rules:
#   Q1 (data) Is the long-chain RISE in Kim's grain BAF robust to the non-detects? A step
#      counts as robust when the lower bound of the longer chain exceeds the upper bound
#      of the shorter one. (Seen already: PFUnDA lo 19.8 > PFDA hi 11.3 and PFNA hi 7.07
#      -> the C11 rise is robust; PFDoDA lo 5.5 < PFDA hi 11.3 -> the C12 level is not.)
#   Q2 (model) log10 RMSE per variant two ways: POINT (vs the ND=0 value, the existing
#      score) and INTERVAL (error 0 inside the bounds, else distance to the nearer bound).
#   VERDICT RULE: a variant's Kim advantage is "robust to the non-detects" only if it beats
#      every other variant by more than 0.10 log10 on BOTH the point and the interval score.
#   EXPECTATION (stated so it can be wrong): lipid's advantage over monotone and W2
#      survives the interval score (monotone sits below even the lower bounds of the long
#      chains); every RMSE shrinks under the interval score; no prediction for nstem_leaf
#      + lipid vs 4-pool lipid.
#
#   python validation/kim2019_censored.py
# =============================================================================
raise SystemExit("pre-registration only; the analysis is added in the next commit")
