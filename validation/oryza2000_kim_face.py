#!/usr/bin/env python3
# =============================================================================
# validation/oryza2000_kim_face.py
# -----------------------------------------------------------------------------
# AN INDEPENDENT CHECK OF THE ORYZA2000 PORT (and the shipped oryza_growth) ON THE
# SHIZUKUISHI RICE FACE EXPERIMENT (Kim et al. 2001 New Phytol.; Kim et al. 2003 GCB)
#
# PRE-REGISTRATION -- written and pushed before either driver was run on the Shizukuishi
# weather or compared with these measurements. (Already done before this file: the
# measurements were transcribed in an earlier session, so their values have been seen;
# the weather was downloaded and its header read; NO model was run on it.)
#
# Why: validation/oryza2000_ir72ds92.py showed the port reproduces IR72DS92 and beats the
# shipped driver there -- but IR72DS92 is the ORYZA2000 manual's own calibration example,
# so for the port that was a FIDELITY test, not an independent one. Before deciding
# whether biomass="oryza2000" should become the default, this file asks whether the
# port's advantage survives on an experiment neither driver was built or tuned on:
# different cultivar (Akitakomachi, a temperate japonica), site (Shizukuishi, 39.6 N),
# years (1998-2000) and lab. (Bouman & van Laar 2006 Fig. 2b, the other organ-biomass
# series in the repo, is the SAME IR72DS92 experiment and is therefore not used.)
#
# Data (docs/literature_db/raw_si/rice_organ_biomass_literature.csv; ambient-CO2 plots,
# medium N = 8 g N m-2 in 1998, 9 in 1999-2000; 19.05 hills m-2, 3 seedlings per hill;
# transplanted 1998-05-21, 1999-05-20, 2000-05-22):
#   Kim 2003 Table 2   green leaf blade GL, senescent leaf blade SL (g m-2) and green LAI
#                      at 5 samplings per year -- exact, tabulated.
#   Kim 2003 Fig. 2    total biomass TOT (g m-2) at the same samplings -- digitised
#                      (+-40 g m-2); root inclusion implied by the Methods (Kim 2001 Table
#                      1's whole-plant total at PI, 425.9, matches the 1999 point, 435).
#   Kim 2001 Table 1   root and total dry matter at PI and maturity, 1999 -- exact; roots
#                      from a 30 x 17.5 x 20 cm soil block, washed, so root is a FLOOR.
# The five samplings are, in order, tillering, panicle initiation (PI), anthesis,
# mid-ripening and grain maturity (Kim 2003 abstract, read via OpenAlex). ORYZA leaves
# are leaf BLADES (sheaths are in stems), so GL ~ WLVG and SL ~ WLVD; SL is a FLOOR for
# WLVD because detached dead blades are not recovered.
#
# THE CULTIVAR PROBLEM, and how it is handled. Both drivers carry IR72 (tropical indica)
# phenology, so neither is expected to time Akitakomachi's stages; comparing on calendar
# days would score the cultivar, not the driver. So every observation is placed by its
# STAGE, not its date: each driver is read on the day its own DVS first reaches 0.65 (PI),
# 1.0 (anthesis), 2.0 (maturity); mid-ripening is read at the same FRACTION of the
# anthesis->maturity interval as the sampling sat in that year (1998 0.56, 1999 0.49,
# 2000 0.52 -- from the sampling DATs, not from any model). Tillering is NOT scored: it
# has no DVS definition, its totals (50-95 g m-2) carry +-40 g m-2 digitising error, and
# it is dominated by the seedbed assumptions below. Every scored quantity is a RATIO, so
# the shipped driver's anchor scaling and the absolute yield level (cultivar, N, CO2)
# drop out.
#
# Scored quantities (per year; MN, ambient CO2):
#   Q1  green-leaf time course   GL(s)/GL(anthesis), s = PI, mid-ripening, maturity     (9)
#   Q2  senescent share           SL/(GL+SL) at anthesis, mid-ripening, maturity         (9)
#   Q3  leaf share of the plant   GL/TOT at PI, anthesis, mid-ripening, maturity        (12)
#   Q4  dry-matter timing         TOT(s)/TOT(maturity), s = PI, anthesis, mid-ripening   (9)
#   Q5  root share                root/total at PI and maturity, 1999 (Kim 2001)          (2)
#   Model counterparts: GL = WLVG, SL = WLVD (port) / the shipped driver's dead leaves
#   accumulated from its own death rate (drlv * WLV, as in oryza2000_ir72ds92.py), TOT =
#   green + dead leaves + stems (incl. reserves) + panicles + roots.
#   Per group: log10 RMSE over its values. POOLED = the mean of the five group RMSEs (so
#   Q3's 12 values do not outweigh Q5's 2).
# Reported, NOT scored (cultivar-specific or not comparable):
#   phenology in days after transplanting (both drivers are IR72); the anthesis->maturity
#   duration (observed 50 / 35 / 42 d; the shipped DVRR 1.05e-3 vs IR72.D92's 1.784e-3
#   makes this differ between the drivers, but grain-fill thermal time is still a cultivar
#   trait); absolute masses; specific leaf area GLAI/GL; the port's spikelet sterility
#   (IR72's cold-sterility response applies at 39.6 N).
#
# Drivers:
#   A  oryza2000.simulate, IR72.D92 unchanged, Management(lat 39.63, nh 19.05, nplh 3,
#      co2 368 [ambient 1998-2000], emergence = transplanting - SBDUR, SBDUR 30 d, NPLSB
#      1000 and the rest of the seedbed as in IR72DS92). The nursery practice (Japanese
#      mat seedlings, often under cover) is NOT reported in the repo's transcription, so
#      SBDUR is an ASSUMPTION; run until DVS > 2 or the weather ends (Nov 30).
#   B  oryza_growth (shipped OryzaParams; lat_deg 39.63, sow_doy = transplanting, season =
#      days to Nov 30), started at transplanting as it always is, anchored (the PFAS-default
#      configuration; ratios do not depend on the anchor).
#   If a driver has not reached DVS 2 by Nov 30, its maturity (and mid-ripening) values are
#   read on Nov 30 and the failure is reported.
#
# Weather: NASA POWER daily for the Shizukuishi cell (raw_si/nasa_power_shizukuishi_
#   1998_2000.csv), UNCORRECTED -- PRIMARY. The cell sits at 461 m vs a ~200 m paddy and no
#   station days exist here to correct it, so: SENSITIVITY S1 = +1.7 C on Tmin and Tmax
#   (lapse rate 6.5 C/km to ~200 m, itself an assumption); S2 = SBDUR 20 d; S3 = SBDUR
#   40 d (S2/S3 move only A).
#
# Decision rules (PRIMARY treatment):
#   D1 (the decision this file exists for) -- the independent evidence SUPPORTS the switch
#      to biomass="oryza2000" if ALL hold:
#        (a) A's POOLED log10 RMSE <= B's POOLED + 0.05;
#        (b) A is not worse than B by more than 0.05 on Q1 or on Q2 (the green-leaf course
#            and senescence feed the PFAS leaf-loss flux, the term the driver most directly
#            sets);
#        (c) A reaches DVS 2 by Nov 30 in all three years.
#      Otherwise D1 FAILS: the IR72DS92 advantage does not transfer and the default stays,
#      with the reason recorded.
#   D2 Robustness: D1 is re-evaluated under S1, S2, S3; if its verdict flips under any of
#      them it is reported as "sensitive to <treatment>".
#   D3 Descriptive adequacy per group (both drivers): log10 RMSE <= 0.10 good, <= 0.20
#      fair, > 0.20 poor. Not a gate.
#   D4 This file changes NO default. A passing D1 means only that the switch (biomass=
#      "oryza2000" + a season tied to crop maturity + the ORYZA-based re-fits) can be
#      pre-registered as one step, for the user to decide.
#
# EXPECTATION (stated so it can be wrong; from reading the parameter tables, not from runs):
#   both drivers share ORYZA's structure and IR72-like partitioning, so on Q3 and Q5 they
#   land within 0.05 of each other; both over-predict the senescent share (Q2: SL is a
#   floor) and under-predict the green leaf left at maturity (Q1: model ~0.3 of the
#   anthesis value vs observed 0.41-0.59); A is better on Q4 at PI (B grows too slowly
#   early, as on IR72DS92) and on the anthesis->maturity duration (B's ~1.7x longer
#   grain-fill thermal time); B may not reach maturity by late October in at least one
#   year. D1 PASSES and does not flip under S1-S3.
#
#   python validation/oryza2000_kim_face.py
# =============================================================================
from __future__ import annotations

import sys

if __name__ == "__main__":
    sys.exit("pre-registration only -- the comparison is implemented in the next commit")
