#!/usr/bin/env python3
# =============================================================================
# validation/oryza2000_kim_face_phenology.py
# -----------------------------------------------------------------------------
# OPTION C: CALIBRATE BOTH BIOMASS DRIVERS' DEVELOPMENT RATES TO AKITAKOMACHI, THEN
# RE-RUN THE SHIZUKUISHI FACE COMPARISON (follow-up to validation/oryza2000_kim_face.py)
#
# PRE-REGISTRATION -- written and pushed before any driver was run with calibrated rates.
# (Already seen: every result of oryza2000_kim_face.py -- the pre-registered primary verdict
# FAIL, the post-hoc station run, and the stage-matched ratios of both drivers with IR72
# phenology. The expectation below is written with that knowledge and says so.)
#
# Why: in oryza2000_kim_face.py both drivers ran with IR72 phenology. Two things in that
# verdict may be phenology, not biomass structure:
#   * the port's ripening lasted 49-66 d vs observed 35-50 d; ORYZA's leaf death is a rate
#     per DAY, so a longer ripening piles up dead leaves -- part of the port's dead-leaf
#     excess (Q2x 0.328 vs the shipped driver's 0.172) may be this;
#   * the shipped driver never reached maturity, so its "maturity" values were read on
#     Nov 30 at DVS < 2 -- not a stage at all.
# Option C removes the phenology mismatch the way ORYZA itself prescribes for a new
# cultivar (manual Section 8.2, program DRATES: development rates from observed dates of
# sowing, transplanting, panicle initiation, flowering and maturity), then asks whether the
# port beats the shipped driver on organ shapes once both develop like Akitakomachi.
#
# What is calibrated -- ONLY the four development rates of each driver:
#   DVRJ (juvenile, DVS 0-0.4), DVRI (photoperiod-sensitive, 0.4-0.65), DVRP (panicle
#   formation, 0.65-1), DVRR (grain filling, 1-2). Everything else -- partitioning, leaf
#   death, SLA, assimilation, the transplanting shock, PPSE = 0 -- stays IR72.D92 (port) /
#   the shipped OryzaParams (shipped driver). So any difference left is a difference in
#   what the drivers do structurally.
# How (DRATES-style, sequential root-finding on the drivers' own heat units):
#   1. DVRJ and DVRI by ONE common factor on their IR72 values, so DVS 0.65 is reached on
#      the observed PI day. (The juvenile/photoperiod split is not observable from these
#      data: there is no observation at DVS 0.4.)
#   2. DVRP so DVS 1.0 is reached on the observed anthesis day.
#   3. DVRR so DVS 2.0 is reached on the observed maturity day.
#   "Reached" = the fractional day at which DVS crosses the threshold (linear between the
#   daily states), so the root-finding is smooth. Observed days are the sampling DATs at
#   PI / anthesis / maturity (1998 53/81/131, 1999 62/88/123, 2000 56/78/120): sampling
#   at a stage, not the exact event day -- an error of a few days, stated not corrected.
#   The port runs from emergence = transplanting - SBDUR, so its PI target depends on the
#   assumed nursery: SOWING / NURSERY DURATION IS NOT IN THE REPO'S TRANSCRIPTION OF KIM
#   2001/2003, so SBDUR = 30 d is an ASSUMPTION (Japanese mat seedlings), with 20 and 40 d
#   as sensitivity. The seedbed cover (TMPSB) stays 0 as in IR72DS92. The shipped driver
#   starts at transplanting (DVS 0) as it always does.
#
# Out-of-sample by construction -- LEAVE-ONE-YEAR-OUT (primary): for each held-out year the
#   rates are the mean of the two other years' own calibrated rates (common J/I factor
#   averaged; DVRP, DVRR averaged); the held-out year is then run with them and scored.
#   SECONDARY: each year with its OWN calibrated rates (phenology exact by construction --
#   the cleanest test of the organ shapes, but in-sample for phenology).
#
# Data, weather, scoring -- as oryza2000_kim_face.py, with two changes fixed in advance:
#   * weather PRIMARY = the JMA station temperatures + POWER radiation (that file's S4; the
#     observed weather). SENSITIVITY: SBDUR 20 / 40 d (station); POWER cell + 1.7 C (S1).
#   * Q2 is scored WITHOUT the anthesis point (that file's Q2x): the shipped driver has no
#     leaf death before DVS 1, so its anthesis senescent share is exactly 0 and log-scores
#     to infinity; the anthesis point is reported, not scored. POOLED = mean of Q1, Q2x,
#     Q3, Q4, Q5.
#   Scored quantities stay STAGE-MATCHED ratios on each driver's own DVS (Q1 green-leaf
#   course, Q2x senescent share, Q3 leaf share of the plant, Q4 dry-matter timing, Q5 root
#   share). REPORTED, NOT SCORED: phenology prediction error in days (LOYO); the same ratios
#   read on the OBSERVED sampling days (calendar-matched, which adds the phenology error);
#   the calibrated rates themselves and how far they move from IR72.
#
# Decision rules (PRIMARY = LOYO, station weather, SBDUR 30):
#   C1 (the decision) -- with calibrated phenology the independent evidence SUPPORTS the
#      switch to biomass="oryza2000" if ALL hold:
#        (a) port POOLED <= shipped POOLED + 0.05;
#        (b) the port is not worse than the shipped driver by more than 0.05 on Q1 or Q2x;
#        (c) the port reaches DVS 2 by Nov 30 in all three held-out years.
#      Otherwise C1 FAILS and the default stays (option A), with the reason recorded.
#   C2 (what the dead-leaf excess was) -- if the port's Q2x falls by >= 0.10 log from its
#      uncalibrated station value (0.328), the excess is recorded as largely PHENOLOGICAL;
#      otherwise as sitting in IR72's leaf-death table (DRLVT), which this data cannot
#      re-fit (ORYZA's PARAM needs stem and panicle weights on the same dates; Kim has leaf
#      blades, totals and roots only).
#   C3 (phenology prediction, descriptive) -- LOYO |error| <= 5 d at anthesis and <= 7 d at
#      maturity in every year = adequate. Not a gate.
#   C4 Robustness: C1 is re-evaluated under every sensitivity treatment; a flip is reported
#      as "sensitive to <treatment>".
#   C5 This file changes NO default and no parameter file. The calibrated rates are an
#      Akitakomachi-at-Shizukuishi result, not a temperate-japonica parameter set: the PFAS
#      datasets use other cultivars (Koshihikari/Kasalath, Nipponbare, Jinyu, unknown).
#
# EXPECTATION (stated so it can be wrong; from the uncalibrated run and the parameter
#   tables, NOT from any calibrated run): calibration shortens the port's ripening, so its
#   dead-leaf share falls, but its leaf death starts at DVS 0.6 and its DRLVT is steeper than
#   the shipped one, so it stays above the observed floor; the port's Q2x falls by LESS than
#   0.10 (C2: the excess is mostly DRLVT). The shipped driver, once it matures in 35-50 d,
#   lands near the observed senescent share. The port keeps its Q4 lead; the shipped driver
#   keeps its Q3 lead. C1 FAILS on (b). LOYO phenology is adequate (C3). The verdict does not
#   flip under the sensitivities.
#
#   python validation/oryza2000_kim_face_phenology.py
# =============================================================================
from __future__ import annotations

import sys

if __name__ == "__main__":
    sys.exit("pre-registration only -- the calibration is implemented in the next commit")
