#!/usr/bin/env python3
# =============================================================================
# validation/oryza2000_kim_face_gao_japonica.py
# -----------------------------------------------------------------------------
# A JAPONICA CROP SET IN THE ORYZA2000 PORT: Gao et al. 2021 (cold-region japonica,
# Suijing-18, Heilongjiang) transferred to the Shizukuishi FACE data (Akitakomachi)
#
# PRE-REGISTRATION -- written and pushed before the port was run with the Gao crop set on
# any weather. Already seen, and the expectation below is written knowing it: every result
# of oryza2000_kim_face.py (5b), oryza2000_kim_face_phenology.py (option C) and
# oryza2000_kim_face_nursery_events.py (the redo): with the recorded nursery the IR72.D92
# port PASSES option C's rule C1 against the shipped driver (N: pooled 0.220 vs 0.244), but
# with the 1999 event dates it FAILS on (b) (NE: dead-leaf share Q2x 0.281 vs 0.181),
# because it over-predicts the senescent share in every year (+0.10...+0.37 log) -- the
# excess option C's C2 put in IR72's leaf-death table DRLVT, which Kim's data cannot
# re-fit (no stem/panicle weights on the same dates).
#
# Why: the port's remaining weakness on independent data is a CROP-FILE question (IR72.D92
# is a tropical indica), not a code question. Gao et al. 2021 (Agric. Water Manag. 250:106864;
# docs/literature_db/raw_si/gao2021_cold_japonica_oryza_params.csv) publish an ORYZA2000
# crop set calibrated on a temperate (cold-region) japonica: development rates, leaf/stem/
# panicle partitioning, a leaf-death table, an SLA table, RGRLMX, FRPAR, a single KDF and
# COLDMIN. Its calibration data (Qing'an 2018) cannot be re-run here -- the daily weather
# was not published -- so this is a TRANSFER test: does a japonica crop set, calibrated
# on another japonica at another site, fix the port's organ shapes on Kim's data?
#
# Arms (all three with phenology calibrated to Akitakomachi, DRATES-style, exactly as in
# the redo's N design; the calibration overwrites each arm's DVRJ..DVRR, so the arms differ
# ONLY in their non-phenology parameters):
#   A  the port with IR72.D92 (= the redo's arm A).
#   G  the port with Gao's calibrated values for FLVTB, FSTTB, FSOTB, DRLVT, SLATB (switch
#      SWISLA = TABLE), KDF (as a constant KDFTB), RGRLMX, FRPAR, and Gao's DVRJ..DVRR as the
#      base the calibration starts from. Everything Gao does not list stays IR72.D92 (as in
#      Gao's own crop file). NOT applied: COLDMIN (the port has no cold-stop rule; Gao lowers
#      it from 12 to 8 oC, i.e. it binds less, and Shizukuishi's season means stay well above
#      8 oC before maturity). Gao's management values (5 plants/hill, 20 hills/m2) are NOT
#      used -- Kim's management stays (3 plants/hill, 19.05 hills/m2).
#   B  the shipped oryza_growth (= the redo's arm B).
#
# Designs (LOYO, station weather, as the redo):
#   N    PRIMARY. Recorded nursery (SBDUR 14/23/24 d, TMPSB 9 oC), sampling-day targets,
#        option C's stage placement.
#   NE   SECONDARY. N + the 1999 event dates (PI 54, heading 78), sample-DVS placement --
#        the design on which arm A failed.
#   Sensitivities of N (rule C1 for G vs B only): (s1) TMPSB 0; (s2) nursery from sowing.
#   DESCRIPTIVE (reported, no rule):
#     * G with Gao's OWN development rates (no calibration; IR72-like 5b scoring on its own
#       DVS) -- does Gao's phenology transfer? Phenology error in days per year and stage.
#     * One group at a time swapped into IR72.D92 (each with calibrated phenology), design N:
#         A+DRLVT   Gao's leaf-death table only;
#         A+PART    Gao's FLVTB/FSTTB/FSOTB only;
#         A+CANOPY  Gao's SLATB (TABLE), KDF, RGRLMX, FRPAR only.
#       Says which part of the crop set moves which score.
#
# Scoring: option C's stage-matched ratios (Q1 green-leaf course, Q2x senescent share
#   without the anthesis point, Q3 leaf share, Q4 dry-matter timing, Q5 root share),
#   log10 RMSE pooled over years, POOLEDx = mean of the five.
#
# Decision rules:
#   G1 (primary, N) -- the Gao leaf-death table FIXES the dead-leaf excess if
#      Q2x(G) <= Q2x(A) - 0.05; otherwise recorded as NOT fixed by it.
#   G2 -- option C's rule C1 with G in place of A, on N and on NE:
#        (a) G POOLEDx <= B POOLEDx + 0.05; (b) G not worse than B by > 0.05 on Q1 or Q2x;
#        (c) G reaches DVS 2 in all three held-out years.
#      Consequence recorded:
#        passes N AND NE  -> a japonica crop set makes the port's pass robust to the 1999
#                            events: a default switch would need a japonica crop file, not
#                            IR72.D92 (Gao's set is a CANDIDATE, not adopted -- G5);
#        passes N, fails NE -> the same fragility as IR72.D92;
#        fails N           -> the Gao set does not transfer to Akitakomachi.
#   G3 -- no general loss: G POOLEDx <= A POOLEDx + 0.02 on N (reported with G1/G2).
#   G4 -- robustness: G2 on N re-evaluated under (s1), (s2); a flip is reported.
#   G5 -- no default, no parameter file and the port's default crop file (IR72.D92) do not
#      change. The Gao values live in this file's crop override only.
#
# Limits stated in advance: Gao calibrated with ORYZA2000's water- and N-limited modules
#   (PADDY, N balance), so their crop values may compensate for stresses this potential-
#   production port does not simulate; their calibration matched GREEN leaf weight, not dead
#   leaves, so their DRLVT is tuned to green-leaf decline, not to the senescent share scored
#   here; Suijing-18 at 47 N is a different (earlier) cultivar than Akitakomachi at 39.6 N.
#
# EXPECTATION (stated so it can be wrong; from the parameter tables only, NOT from any run
#   with the Gao set):
#   * Gao's leaf death starts at DVS 0.4 (0.010 d-1) and runs ABOVE IR72's up to DVS ~1.3
#     (0.024 vs 0.015 at DVS 1), BELOW it after DVS 1.6 (0.020 vs 0.025-0.050). Dead leaves
#     accumulate in WLVD, so the early extra death offsets the late saving at maturity and
#     raises the share at mid-ripening: G1 NOT fixed (Q2x falls by < 0.05, may rise).
#     A+DRLVT alone behaves like G on Q2x.
#   * Gao's lower leaf fraction at DVS 0.5 (0.50 vs 0.60) and the early leaf death lower the
#     leaf share: the port already under-predicts it (Q3 bias -0.076), so Q3 gets worse.
#   * G3 fails (G pooled > A pooled + 0.02); G2 fails on NE (as A) and its N verdict is not
#     better than A's.
#   * Gao's OWN development rates sit inside the range option C calibrated for Akitakomachi
#     (DVRJ x1.15, DVRP x1.35, DVRR x1.12 of IR72's), so uncalibrated they come within 7 d of
#     the observed anthesis and 10 d of maturity in every year (descriptive).
#
#   python validation/oryza2000_kim_face_gao_japonica.py
# =============================================================================
from __future__ import annotations

import sys

if __name__ == "__main__":
    sys.exit("pre-registration only -- the analysis is implemented in the next commit")
