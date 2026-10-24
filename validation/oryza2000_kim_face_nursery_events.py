#!/usr/bin/env python3
# =============================================================================
# validation/oryza2000_kim_face_nursery_events.py
# -----------------------------------------------------------------------------
# OPTION C REDONE WITH THE RECORDED NURSERY AND THE 1999 EVENT DATES
# (POST-HOC ADDENDUM to validation/oryza2000_kim_face_phenology.py)
#
# PRE-REGISTRATION -- written and pushed before any driver was run with these settings.
# Already seen, and the expectation below is written knowing it: every result of option C
# (oryza2000_kim_face_phenology.py, re-run after the Kim 2003 Fig. 2 re-read): C1 PASSES,
# pooled 0.222 (port) vs 0.244 (shipped); Q1 0.174/0.287, Q2x 0.256/0.221, Q3 0.148/0.133,
# Q4 0.207/0.306, Q5 0.324/0.274; own-year 0.246/0.265; port J/I factor x1.352/1.092/1.156
# (1998/1999/2000); LOYO port phenology error PI +10.6/-6.3/-2.4, anthesis +10.5/-10.3/+0.7,
# maturity +2.0/-2.1/-3.4 d; C2 drop 0.072; C3 NOT adequate; its SBDUR sensitivity moved
# the port's pooled score 0.233 (20 d) / 0.222 (30 d) / 0.214 (40 d), verdict unchanged.
#
# Why: option C used two inputs that sources transcribed AFTER it ran now contradict.
#   N1. NURSERY. Option C ran the port from emergence = transplanting - 30 d on an OPEN
#       seedbed (SBDUR 30 assumed, TMPSB 0). Kim et al. 2003 state the seedlings were
#       transplanted 14, 23 and 24 d after emergence (1998/1999/2000), raised in trays in
#       PLASTIC CHAMBERS; Sasaki et al. 2005 confirm 23 d and the plastic chambers for 1999
#       (docs/literature_db/raw_si/sasaki2005_shizukuishi_1999.csv). ORYZA's input for a
#       covered seedbed is TMPSB, the rise of the daily maximum temperature under the cover
#       while CROPSTA <= 2; the ORYZA2000 manual gives "In Japan, the maximum temperature is
#       raised by 9 oC (Horie)" (its example file: 0 without cover, 9.5 with). The seedbed
#       matters beyond the emergence date: the port's transplanting shock halts development
#       for SHCKD x TSTR heat units (SHCKD 0.4, TSTR = heat sum at transplanting), and the
#       seedbed heat sets the seedling size at transplanting. The shipped driver starts at
#       transplanting (DVS 0), so N1 changes ONLY the port arm.
#   E1. EVENTS (1999 only). Option C calibrated to the SAMPLING days Kim labels PI /
#       anthesis / maturity (1999: DAT 62/88/123). Sasaki 2005 -- same experiment, same
#       year -- gives the events: PI ~54 DAT, heading 78 DAT (ambient plots), and calls the
#       DAT-62 and DAT-88 samplings "pre-heading" and "early grain filling". So in 1999 the
#       targets were 8 and 10 d late, and the two samples were not at DVS 0.65 / 1.0.
#       No event dates exist for 1998 or 2000 in any source in hand; no maturity date exists
#       for any year (Sasaki's last dry-weight sampling is DAT 119, Kim's DAT 123).
#
# What changes (nothing else does -- rates calibrated the DRATES way, LOYO, station weather,
# Q2 scored without the anthesis point, C1 rule, all exactly as option C):
#   N1  port SBDUR = 14 / 23 / 24 d (Kim 2003), TMPSB = 9 oC (ORYZA manual, Japan).
#   E1  1999 targets PI 54, DVS 1.0 at heading 78 (ORYZA's DVS 1 is flowering; rice flowers
#       ~1-3 d after heading -- taken as equal, stated not corrected), maturity stays the
#       sampling DAT 123. 1998/2000 targets stay the sampling DATs (option C). NOT applied:
#       shifting 1998/2000 by the 1999 lag. The samplings were timed by the observed crop,
#       and 1998 was the coolest year: subtracting 1999's 10 d would put 1998 heading at
#       71 DAT, before the warmer 1999's 78 -- not plausible.
#   PLACEMENT under E1: a 1999 sample is no longer at a target DVS, so every sample is
#       placed at the DVS that the driver's OWN-YEAR-calibrated run has on the sampling day
#       ("sample DVS"), and on any other run (LOYO) it is read where that run reaches the
#       sample DVS. For a sample taken on a target day this is exactly DVS 0.65 / 1.0 / 2.0,
#       so for 1998/2000 PI, anthesis and maturity it is option C's placement. It differs
#       from option C for MID-RIPENING (option C: at the observed day-fraction between
#       anthesis and maturity), so a control isolates the placement from the event dates.
#
# Designs (all LOYO, station weather):
#   N    PRIMARY. Option C with N1 only; option C's placement. The shipped arm is untouched
#        by N1 and must reproduce option C's shipped numbers exactly (a self-check).
#   NG   CONTROL. N with the sample-DVS placement, targets unchanged (= option C's).
#   NE   SECONDARY. N + E1 with the sample-DVS placement.
#   Sensitivities of N: (s1) TMPSB 0 (cover ignored); (s2) SBDUR counted from SOWING,
#        20 / 27 / 29 d (Hasegawa et al. 2017 Table S1 sowing dates; treats sowing as
#        emergence), TMPSB 9; (s3) POWER cell + 1.7 oC (option C's S_power+1.7C), N1.
#   NE is also run under (s1).
#   N is primary because the nursery is recorded for all three years while the events are
#   known for one; NE uses better information for 1999 only and mixes event targets with
#   sampling targets, so it is a check that the 1999 lag does not drive the verdict.
#
# Decision rules:
#   D1 (primary) -- C1 of option C, unchanged, on N: (a) port POOLEDx <= shipped + 0.05;
#      (b) port not worse by > 0.05 on Q1 or Q2x; (c) the port reaches DVS 2 in all three
#      held-out years. If the verdict differs from option C's PASS, option C's verdict is
#      recorded as dependent on the assumed nursery.
#   D2 -- the same rule on NE and on NG. If NE and N disagree, the verdict is recorded as
#      sensitive to the 1999 event dates; NG says whether a disagreement is the placement.
#   D3 -- the same rule under each sensitivity; a flip is reported as "sensitive to <s>".
#   D4 -- option C's C2 on N: does the port's Q2x fall by >= 0.10 from its uncalibrated
#      station value (0.328)?
#   D5 (descriptive) -- LOYO phenology error per year and stage under N and NE; under NE
#      the 1999 errors are against the events (PI 54, heading 78). Option C's C3 threshold
#      (|anthesis| <= 5 d, |maturity| <= 7 d every year) is reported for both.
#   D6 -- no default, no parameter file and no recorded option C number changes; option C
#      stays the record of what it pre-registered, this file sits beside it.
#
# EXPECTATION (stated so it can be wrong; from option C's results and a heat-sum estimate
#   of the seedbeds, NOT from any run with these settings):
#   * N: the calibration absorbs N1 mostly into the port's J/I factor. With TMPSB 9 the
#     1999 and 2000 covered nurseries (23/24 d) accumulate MORE seedbed heat than option C's
#     30 open days, the 14-d 1998 nursery about as much, so the 1999/2000 factors fall and
#     1998 stays the fastest year: the year-to-year spread of the factor does NOT shrink,
#     and the LOYO PI errors do not improve. The port's pooled score moves by < 0.02 (option
#     C's own SBDUR sensitivity moved it 0.019); D1 PASSES; the shipped arm is identical.
#   * D4: the Q2x drop stays < 0.10 (the dead-leaf excess is IR72's leaf-death table).
#   * NE: 1999's targets move 8-10 d earlier, so 1999's J/I factor and DVRP rise toward
#     1998's. LOYO anthesis errors shrink in magnitude: 1999 within 3 d of the true heading,
#     1998 below option C's +10.5, 2000 turns negative. C3 stays NOT adequate (the maturity
#     errors come from DVRR, which E1 does not touch). D2: NE PASSES; NG and N agree within
#     0.01 on the pooled scores.
#   * D3: no flip under (s1)-(s3).
#
#   python validation/oryza2000_kim_face_nursery_events.py
# =============================================================================
from __future__ import annotations

import sys

if __name__ == "__main__":
    sys.exit("pre-registration only -- the analysis is implemented in the next commit")
