#!/usr/bin/env python3
# =============================================================================
# validation/oryza2000_ir72ds92.py
# -----------------------------------------------------------------------------
# ORYZA2000 (IR72.D92) PORT vs THE SHIPPED oryza_growth, ON THE IR72DS92 EXPERIMENT
# (handoff step 5, docs/HANDOFF_mass_basis_refit.md)
#
# PRE-REGISTRATION -- written and pushed before either model was run on the 1992 weather
# or compared with the IR72DS92 measurements. (Already done before this file: the port
# was run on SYNTHETIC weather to check it runs and that its carbon balance closes; the
# IR72DS92 measurements were transcribed, so their values have been seen.)
#
# Why: the per-organ biomass M_k(t) drives every PFAS concentration (growth dilution,
# terminal accumulation). The default driver, src/oryza_growth.py, is a REDUCED Level-1
# re-coding whose parameters are "IR72 order of magnitude" -- not IR72.D92 (e.g. stems
# keep 25-30% of post-flowering growth vs 0 in IR72.D92; maturity 116 d after
# transplanting vs 98 d observed). src/oryza2000.py ports ORYZA2000 potential production
# line by line from the manual's Fortran with IR72.D92 verbatim. This file asks how each
# reproduces the experiment the ORYZA2000 manual ships with IR72.D92.
#
# Data: docs/literature_db/raw_si/oryza_ir72ds92_obs.csv (manual Section 7.2, IR72DS2.T92):
# IRRI 1992 dry season, emergence doy 4, transplanted doy 16, 25 hills x 5 plants, 225 kg
# N/ha, irrigated. NOTE the framing: the manual uses this very file as the example input
# of its calibration programs (DRATES/PARAM) and IR72.D92 was "derived from various
# experiments at IRRI". For the PORT this is therefore a FIDELITY test (a faithful port
# with the right weather should reproduce its own calibration experiment), NOT an
# independent validation of ORYZA2000. For the SHIPPED driver it measures how far it is
# from ORYZA2000-IR72 behaviour. The run was N-balance/water-balance with ample inputs,
# i.e. near-potential, not strictly potential.
#
# Weather: the station file PHIL1.992 is not obtainable here (requested; login-gated).
#   PRIMARY until it arrives: NASA POWER daily (raw_si/nasa_power_losbanos_1992.csv)
#   corrected to the station by the 12 days the manual prints (raw_si/
#   irri_phil1_1992_printed.csv): Tmin +1.17 C, Tmax +1.67 C (additive), radiation x1.211
#   (ratio of means). SENSITIVITY: temperature-only correction; no correction. All three
#   are reported. If PHIL1.992 becomes available, the station run SUPERSEDES these (same
#   rules, same scoring) and the POWER runs stay as the record of the fallback.
#
# Arms:
#   A  oryza2000.simulate, IR72.D92, IR72DS92 management, CO2 340, from emergence (doy 4).
#   B  oryza_growth (shipped OryzaParams, lat 14.18), started at transplanting (sow_doy 16,
#      its day 0 = transplanting), the same weather from doy 16, anchored to its 1740 g/m2
#      shoot target (scale_to_anchor=True: the configuration the PFAS default runs on).
#      Its dead leaves are accumulated from its own death rate (drlv * WLV) and scaled by
#      the same anchor factor. SENSITIVITY: B unanchored.
#
# Scoring: main-field dates only (doy 34, 58, 68, 83, 97, 114) -- the doy-16 record is the
# seedbed/transplanting transient (ORYZA writes the pre-dilution state on that day).
# Variables WLVG, WLVD, WST, WSO, WAGT (= green + dead leaves + stem + panicle) and LAI,
# zero observations excluded; per variable log10 RMSE. Phenology: day of year at which DVS
# first reaches 0.65 / 1 / 2 vs observed 58 / 83 / 114.
#
# Decision rules (on the PRIMARY weather):
#   R1 (fidelity of the port) ADEQUATE if |d flowering| <= 5 d, |d maturity| <= 5 d,
#      log10 RMSE(WAGT) <= 0.10 and log10 RMSE <= 0.15 for each of WLVG, WST, WSO.
#      Otherwise NOT ADEQUATE; if it passes under another weather treatment, the verdict
#      is recorded as weather-limited.
#   R2 (A vs B) A fits IR72DS92 better if its |d maturity| is smaller AND its log10 RMSE is
#      smaller on at least 3 of the 4 series WLVG, WST, WSO, WAGT.
#   R3 This file does NOT change any default. Whether biomass="oryza2000" replaces "oryza"
#      as the default is decided with the user after seeing these results (and then needs
#      the ORYZA-based re-fits redone).
#   POST-HOC (labelled, not part of the rules): the size of the PFAS BAF change between the
#      two drivers under the PFAS default climatology -- what a switch would move.
#
# EXPECTATION (stated so it can be wrong): A is ADEQUATE under the primary weather
#   (phenology within ~4 d, WAGT within ~15%); under uncorrected POWER its phenology runs a
#   few days late. B matures ~15-20 d late, keeps growing its stem after flowering (WST
#   over at doy 97/114) and under-fills its panicles before the anchor scale; A beats B by
#   R2.
#
#   python validation/oryza2000_ir72ds92.py
# =============================================================================
