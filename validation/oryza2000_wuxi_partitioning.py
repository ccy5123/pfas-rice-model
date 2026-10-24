#!/usr/bin/env python3
# =============================================================================
# validation/oryza2000_wuxi_partitioning.py
# -----------------------------------------------------------------------------
# A SECOND INDEPENDENT CHECK OF THE BIOMASS DRIVERS -- ORGAN PARTITIONING AT THE WUXI
# RICE FACE (Yang et al. 2006, Field Crops Res. 98:12-19; ambient-CO2 plots, 2001-2003)
#
# PRE-REGISTRATION -- written and pushed before either driver was run on Wuxi weather or
# compared with these data. (Already seen: the transcribed data -- docs/literature_db/raw_si/
# yang2006_wuxi_face_dm.csv -- and Hasegawa et al. 2017's conditions table; every Shizukuishi
# result: 5b, option C, the redo with the recorded nursery and events, and the Gao crop set.
# The weather had not been downloaded when this was written.)
#
# Why: the Shizukuishi data (Kim 2001/2003) have leaf blades, totals and roots only -- no
# stem or panicle weights -- so they cannot test the one thing the PFAS model takes most
# directly from a biomass driver: how dry matter is split between leaf, stem and panicle
# over the season (the grain is the PFAS terminal sink; the straw = stem + leaf). Yang 2006
# reports exactly that split (leaf / stem incl. sheath / spike, summing to 1) plus the
# above-ground dry matter, at five samplings in each of three years, for another japonica
# (Wuxiangjing 14, an early-maturing late japonica) at another site (31.6 N) and lab.
#
# What the data do NOT give (stated, not worked around): heading and maturity DATES (the
# modellers of Hasegawa et al. 2017 received them, but they are not published), roots, and
# whether dead leaves are inside "leaf" or outside the total. Hasegawa et al. 2017 Table S1
# gives Wuxi's mean growth duration, 158 d -- read as sowing to maturity, mean of the years:
# for Shizukuishi the same table's 146 d matches Sasaki's 1999 maturity sampling exactly and
# Kim's maturity samplings within 3-5 d (149-151 d) -- i.e. Wuxi maturity at ~132 days after
# transplanting; and the season mean temperature and radiation per year.
#
# Arms (ambient CO2 per year 354/364/375 umol/mol; 24 hills m-2, 3 plants/hill; sown 18 May,
# transplanted 13 June; lat 31.58 N):
#   A        the ORYZA2000 port, IR72.D92, its own development rates.
#   B        the shipped oryza_growth (anchored -- the anchor scales all organs equally, so
#            every scored ratio is anchor-free), its own development rates.
#   G        the port with the Gao 2021 japonica crop set (as in oryza2000_kim_face_gao_
#            japonica.py, Gao's own development rates).                       [secondary]
#   A+PART   IR72.D92 with Gao's partitioning tables only (FLVTB/FSTTB/FSOTB). [secondary]
#   Port seedbed: emergence 4 d after sowing (assumption) -> SBDUR 22 d; open seedbed
#   (TMPSB 0; late-May Wuxi needs no cover; the paper says only "raised under ambient air").
#
# Weather: NASA POWER daily for the Wuxi cell (31.58 N, 120.50 E; cell elevation ~6 m, site
#   5 m), 2001-03-01..2003-12-31, saved to raw_si/. CHECK fixed in advance: the POWER mean
#   temperature over transplanting -> DAT 132 is compared per year with Hasegawa's season
#   mean (25.2 / 24.7 / 25.4 oC); if any year differs by > 0.5 oC, a sensitivity run shifts
#   that year's Tmin/Tmax by the difference. (The primary run uses POWER as downloaded.)
#
# Placement -- by stage on each driver's OWN DVS (the samplings were "fixed to coincide with"
#   the stages): PI sample at DVS 0.65; heading at DVS 1.0 (heading taken as flowering, as
#   in the Shizukuishi files); mid-ripening at the driver's DVS-1 day + 20 d (2001, 2002) /
#   + 30 d (2003), the paper's "days after heading"; maturity at DVS 2.0. If a driver never
#   reaches DVS 2 before Dec 31, its maturity is read on Dec 31 and flagged. The tillering
#   sample (27-28 DAT) has no stage clock: reported on the calendar day, not scored.
#
# Scored (log10 errors, ambient plots, 3 years pooled; RMSE per group):
#   W1 leaf share  at PI, heading, mid-ripening, maturity
#   W2 spike share at heading, mid-ripening, maturity
#   W3 stem share  at PI, heading, mid-ripening, maturity
#   W4 dry-matter timing AGB(stage)/AGB(maturity) at PI, heading, mid-ripening
#   POOLED = mean of W1-W4.
#   Model organs: leaf = green + dead leaf blades (WLVG + WLVD), stem = WST (sheaths are in
#   ORYZA's stem, as in Yang's), spike = WSO (panicles), total = their sum (PRIMARY reading:
#   dead leaves counted, nothing left out of the total). SENSITIVITY reading: leaf = green
#   only, dead leaves outside the total.
#
# Designs:
#   W   PRIMARY. Own development rates (no phenology data exist to calibrate them).
#   Wp  SECONDARY. Season length pinned with the only phenology information there is: one
#       common factor on DVRJ/DVRI puts DVS 0.65 on the PI sampling day (51/47/45 DAT), and
#       one common factor on DVRP and DVRR together (their ratio kept) puts DVS 2.0 at
#       DAT 132 every year (Hasegawa's mean duration -- the same day for every year is an
#       assumption). Heading then falls where each driver's own DVRP/DVRR ratio puts it.
#   Sensitivities of W: SBDUR 18 / 26 d (port arms); the green-leaf-only reading; the
#   temperature offset (only if the CHECK triggers).
#
# Decision rules:
#   R1 (primary, W) -- the independent partitioning evidence SUPPORTS the port over the
#      shipped driver if ALL hold: (a) A POOLED <= B POOLED + 0.05; (b) A is not worse than
#      B by more than 0.05 on W2 (spike share -- the PFAS grain sink) or W3 (stem share --
#      most of the straw); (c) A reaches DVS 2 before Dec 31 in all three years.
#      Otherwise R1 FAILS, with the failing clause recorded.
#   R2 -- R1 re-evaluated on Wp and under each sensitivity; a flip is reported as
#      "sensitive to <x>".
#   R3 (secondary) -- Gao's japonica partitioning fits better than IR72's if
#      mean(W1, W2, W3) of A+PART <= that of A - 0.05 (design W). Reported for G as well.
#   R4 (descriptive) -- each driver's heading/maturity DAT next to the ~132-DAT maturity;
#      the port's above-ground DM at maturity next to the observed 1674-1829 g m-2 (the
#      shipped driver's absolute level is set by its anchor and is not compared).
#   R5 -- no default and no parameter file changes.
#
# EXPECTATION (stated so it can be wrong; from the IR72DS92 and Shizukuishi results and the
#   parameter tables, NOT from any Wuxi run):
#   * The shipped driver keeps growing its stem after flowering and under-fills panicles (it
#     did on IR72DS92: stem 7095 vs 4243, panicles 6912 vs 9843 kg/ha), so at Wuxi its spike
#     share at mid-ripening and maturity is too low and its stem share too high; the port is
#     closer on W2 and W3. R1 PASSES. The shipped driver does reach DVS 2 at Wuxi (warmer
#     than Shizukuishi), late.
#   * The port with IR72's own rates matures early -- before DAT 115, i.e. >= 15 d before the
#     ~132-DAT maturity (IR72 is a short-duration indica, Wuxiangjing 14 a photoperiod-
#     sensitive late japonica; IR72's PPSE is 0), and its above-ground DM at maturity falls
#     short of the observed.
#   * W1: under the primary (green + dead) reading the port's leaf share at maturity is above
#     the observed (its dead pool is large); W1 is no better for the port than for the
#     shipped driver.
#   * R3: Gao's partitioning is not better than IR72's by 0.05 (both start panicle growth
#     around DVS 0.75 and give a similar spike share by heading).
#   * R2: Wp agrees with W (no flip); the green-only reading does not flip R1.
#
#   python validation/oryza2000_wuxi_partitioning.py
# =============================================================================
from __future__ import annotations

import sys

if __name__ == "__main__":
    sys.exit("pre-registration only -- the analysis is implemented in the next commit")
