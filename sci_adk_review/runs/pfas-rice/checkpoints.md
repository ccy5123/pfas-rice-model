# Agent judgment checkpoints

proof/qualitative hypotheses awaiting an in-session agent verdict (no autonomous LLM call). For each, author verdicts/<hyp-id>.json with the chief-over-N trail, then re-enter the loop (sci-adk resolve <run-dir>).

## hyp-anion (qualitative)
- Criterion: GHK electrodiffusion with measured E_m=-120 mV and z=-1 yields an anion-exclusion factor e^N of order 1e2 (the IOC-formulation signature) => support; no exclusion => refute
- Finding: N = z*E_m*F/RT = +4.67 at E_m=-120 mV, z=-1, T=298.15 K => e^N = 106.8 (test_plant_module asserts e^N ~ 107); the GHK term in root_uptake() excludes the anion accordingly.

## hyp-yamazaki (qualitative)
- Criterion: the Yamazaki agreement is OUT-OF-SAMPLE predictive validation (independent data not used to fit the per-congener transport parameters) => support; an in-sample / saturated reproduction => refute
- Finding: reproduce_demo.py: W2 transport fit reproduces Yamazaki root/straw/grain BAFs across 11 congeners at log10 RMSE 0.029 (e.g. PFOA root 0.49/0.49, straw 0.83/0.83 pred/obs). The fit assigns ~3 transport params per congener against 3 tissue observations => SATURATED: reproduction is structurally guaranteed and is NOT out-of-sample prediction (CLAUDE.md §6).
reproduce_demo.py --rec: with the theory/QSPR MONOTONE f_xy (a-priori, NOT fit to the tissue BAFs) the predictive error is log10 RMSE 0.837 vs the saturated W2 fit's 0.029. Straw is off 6-40x (PFBA 45/11, PFBS 33/2.2). The model does NOT a-priori predict the Yamazaki tissue BAFs -- a quantitative confirmation of the REFUTED verdict.
validation/apriori_prediction.py: redistributed-shoot model (nstem_leaf, N stem segments + retention) lowers the a-priori OOS log10 RMSE 0.987 -> 0.951 (same monotone f_xy + drivers) -- MARGINAL. Short-chain straw improves but long-chain straw/grain still collapse (PFDoDA straw 0.35 vs 49.75); needs hysteretic high-B sorption (docs/nstem_gradient_exploration.md). Predictive claim remains refuted.
validation/refit_oryza.py: per-congener (f_xy,L_Ph,kappa_d) re-fit to Yamazaki on the MECHANISTIC ORYZA2000 biomass (the new default) reproduces in-sample at log10 RMSE 0.236 -- SATURATED (DOF 0), i.e. reproduction not prediction. PFDoDA(C12) hits f_xy=1/L_Ph=1/kappa_d=0.01 ceilings yet stays ~4-6x under (root 19 vs 69, grain 7.4 vs 46) -- a structural long-chain floor. Wired as model_api f_xy_source='oryza' so the default-biomass model reproduces.

## hyp-grain (qualitative)
- Criterion: model grain (brown-rice) BAF/TF matches measured within a factor adequate for dietary risk assessment => support; systematic structural under/over-prediction => refute
- Finding: Tang 2026 OOS: PFOA endosperm TF model 0.11 vs measured 0.95 (dw); grain structurally ~3-8x under across congeners; not closable by L_Ph / lipid tuning (docs/tang2026_grain_units_exploration.md).
Kim 2019: grain BAF reproduced only by FORCING L_Ph (0.07 -> 4.43, L_Ph~0.84) as a single-point in-sample anchor; Kim is grain-only, so root/straw TF are unconstrained.
validation/structural_adequacy_fit.py, driven by the MECHANISTIC ORYZA2000 biomass (oryza_growth) + measured Q_TP (forcing_rice) -- NOT the logistic placeholder. CONSTRAINED (DOF>0) fits vs Yamazaki: A f_xy+global L_Ph+global kappa_d (DOF 20): root 0.45 straw 0.18 grain 0.52 overall 0.41; B +per-cong L_Ph (DOF 10): grain 0.36 overall 0.35; C +per-cong kappa_d (DOF 10): root 0.26 straw 0.16 grain 0.51 overall 0.34. => structure reproduces SHOOT translocation (straw ~0.16-0.18, i.e. within ~1.5x) under a constrained fit; root needs per-congener kappa_d (->0.26); grain improves with per-congener L_Ph (->0.36) but keeps a long-chain residual floor. Whole plant within ~factor 2.2 (overall 0.34) at DOF 10 -- vs saturated W2 0.029 (DOF 0) and a-priori 0.84.

## hyp-smiles (qualitative)
- Criterion: the SMILES front-end reproduces the curated measured-parameter model for a KNOWN PFAS via read-across (structure -> same Compound) => support; a mismatch => refute
- Finding: tests/test_pfas_structure.py: 23 passed (RDKit). A canonical SMILES that matches a curated congener rebuilds the SAME Compound from params/parameters.json (measured read-across). CAVEAT: for NOVEL structures f_xy is provisional (QSPR/interpolated), NOT validated — this claim is scoped to known structures only.

## hyp-adequacy (qualitative)
- Criterion: a CONSTRAINED (degrees-of-freedom > 0) calibration of the structure reproduces the measured STRAW (shoot translocation) BAFs across the C4-C12 PFCA/PFSA series => support; only a saturated (DOF 0) per-congener fit reproduces => refute
- Finding: validation/structural_adequacy_fit.py, driven by the MECHANISTIC ORYZA2000 biomass (oryza_growth) + measured Q_TP (forcing_rice) -- NOT the logistic placeholder. CONSTRAINED (DOF>0) fits vs Yamazaki: A f_xy+global L_Ph+global kappa_d (DOF 20): root 0.45 straw 0.18 grain 0.52 overall 0.41; B +per-cong L_Ph (DOF 10): grain 0.36 overall 0.35; C +per-cong kappa_d (DOF 10): root 0.26 straw 0.16 grain 0.51 overall 0.34. => structure reproduces SHOOT translocation (straw ~0.16-0.18, i.e. within ~1.5x) under a constrained fit; root needs per-congener kappa_d (->0.26); grain improves with per-congener L_Ph (->0.36) but keeps a long-chain residual floor. Whole plant within ~factor 2.2 (overall 0.34) at DOF 10 -- vs saturated W2 0.029 (DOF 0) and a-priori 0.84.
