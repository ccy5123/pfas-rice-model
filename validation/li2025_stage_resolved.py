#!/usr/bin/env python3
# =============================================================================
# validation/li2025_stage_resolved.py
# -----------------------------------------------------------------------------
# LI 2025 STAGE BY STAGE  (handoff step 3, docs/HANDOFF_mass_basis_refit.md)
#
# PRE-REGISTRATION -- written and pushed before any model number below was produced.
#
# Why: the repo's Li comparison (validation/oos_crossdataset.li_tf) sets the model's
# HARVEST-day BAF next to the data's 5-STAGE MEAN concentrations, converted to fresh
# weight with assumed theta. Two mismatches: time (harvest vs season mean) and basis
# (the data are dry weight). The authors' own TF is a third estimator (the arithmetic
# mean of per-stage ratios), 1.0-5.8x above the repo's.
#
# Data: SI Table S6 -> docs/literature_db/raw_si/li2025_stage_conc.csv (ng/g dw; 5
# stages without dates x root/stem/leaf/grain, grain from heading). Transcription
# checked before this was written: the 5-stage means reproduce data_obs/Li2025_BAF_TF.csv
# in all 20 congener x tissue cells, and the mean of per-stage ratios reproduces the
# authors' PFOA stem TF 2.35.
#
# Observable: TF_k = C_k / C_root, dry weight, per stage. It needs no water value (the
# irrigation-water denominator is Li's documented weak point); the root-surface confound
# stays (roots washed with distilled water, no plaque removal).
#
# Model: the dynamic run on the default ORYZA2000 biomass + measured Q_TP, Cwo constant 1;
# model fresh-weight concentration -> dry weight with the binding factor's theta. Stage ->
# model day by ORYZA DVS: tillering 0.4, heading 1.0, filling 1.25, milk 1.5, ripening
# 2.0 (first day it is reached). Variants as in oos_crossdataset (4-pool: monotone f_xy,
# W2 fit, lipid loading) plus nstem_leaf with lipid loading, because Li reports stem and
# leaf separately and the 4-pool stem is a pass-through. Congeners PFBA, PFHxA, PFOA,
# PFBS, PFOS (PFHpA is ND in every tissue). Grain is scored from FILLING on: at heading
# the model's grain-formation gate has barely opened (structural; the value is printed).
#
# Questions and decision rules:
#   Q1 (estimator) How far apart are the three TF estimators on the DATA: ratio of means
#      (repo), arithmetic mean of ratios (authors), geometric mean of ratios? Reported.
#      Scoring uses the geometric mean, the estimator consistent with a log10 metric.
#   Q2 (noise) Across-stage SD of the observed log10 TF vs the SD the reported +/-SDs
#      imply. If the median ratio exceeds 3, the reported SDs do not describe the
#      stage-to-stage variation and per-stage SHAPE results are treated as noise-limited.
#   Q3 (level) log10 RMSE of the stage-geometric-mean TF -- the SAME estimator over the
#      SAME stages on both sides -- per variant, for shoot (stem + leaf) and grain.
#   Q4 (shape) Pooled Spearman between de-meaned model and observed log10 TF trajectories
#      per variant; permutation p (5000 shuffles within each series). A shape signal
#      counts only at p < 0.05.
#   VERDICT RULE: Li moves from "inconclusive" to discriminating only if one variant beats
#      every other by more than 0.10 log10 RMSE on Q3 for BOTH shoot AND grain. Otherwise
#      Li stays inconclusive, and what changes is only how the comparison is built
#      (matched times, dry basis, a stated estimator).
#   EXPECTATION (stated so it can be wrong): Li stays inconclusive; the estimator choice
#      alone moves the observed TF by more than the model variants differ; Q4 finds no
#      shape signal (p > 0.05).
#   SENSITIVITY (fixed in advance): every stage DVS shifted by -0.15 and +0.15 (ripening
#      held at 2.0).
#
#   python validation/li2025_stage_resolved.py
# =============================================================================
raise SystemExit("pre-registration only; the analysis is added in the next commit")
