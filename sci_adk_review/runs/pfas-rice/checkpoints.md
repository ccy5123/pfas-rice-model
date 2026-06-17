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

## hyp-grain (qualitative)
- Criterion: model grain (brown-rice) BAF/TF matches measured within a factor adequate for dietary risk assessment => support; systematic structural under/over-prediction => refute
- Finding: Tang 2026 OOS: PFOA endosperm TF model 0.11 vs measured 0.95 (dw); grain structurally ~3-8x under across congeners; not closable by L_Ph / lipid tuning (docs/tang2026_grain_units_exploration.md).
Kim 2019: grain BAF reproduced only by FORCING L_Ph (0.07 -> 4.43, L_Ph~0.84) as a single-point in-sample anchor; Kim is grain-only, so root/straw TF are unconstrained.
validation/structural_adequacy_fit.py: CONSTRAINED fit (per-congener f_xy + global L_Ph=3.2e-3 + global kappa_d=2.05; 13 params / 33 obs, DOF 20) vs Yamazaki. Per-tissue log10 RMSE: STRAW 0.048 (reproduced across the WHOLE series incl. long chains PFUnDA 8.18/8.16, PFDoDA 34/49), root 0.384 (within ~2-3x on one global kappa_d), GRAIN 0.987 (short over ~10x: PFBA 11 vs 1; long under ~75x: PFDoDA 0.6 vs 46). Overall 0.612 (vs saturated W2 0.029 at DOF 0; a-priori 0.84). => structure reproduces SHOOT translocation under a constrained fit; GRAIN does not (needs per-congener phloem loading).

## hyp-smiles (qualitative)
- Criterion: the SMILES front-end reproduces the curated measured-parameter model for a KNOWN PFAS via read-across (structure -> same Compound) => support; a mismatch => refute
- Finding: tests/test_pfas_structure.py: 23 passed (RDKit). A canonical SMILES that matches a curated congener rebuilds the SAME Compound from params/parameters.json (measured read-across). CAVEAT: for NOVEL structures f_xy is provisional (QSPR/interpolated), NOT validated — this claim is scoped to known structures only.

## hyp-adequacy (qualitative)
- Criterion: a CONSTRAINED (degrees-of-freedom > 0) calibration of the structure reproduces the measured STRAW (shoot translocation) BAFs across the C4-C12 PFCA/PFSA series => support; only a saturated (DOF 0) per-congener fit reproduces => refute
- Finding: validation/structural_adequacy_fit.py: CONSTRAINED fit (per-congener f_xy + global L_Ph=3.2e-3 + global kappa_d=2.05; 13 params / 33 obs, DOF 20) vs Yamazaki. Per-tissue log10 RMSE: STRAW 0.048 (reproduced across the WHOLE series incl. long chains PFUnDA 8.18/8.16, PFDoDA 34/49), root 0.384 (within ~2-3x on one global kappa_d), GRAIN 0.987 (short over ~10x: PFBA 11 vs 1; long under ~75x: PFDoDA 0.6 vs 46). Overall 0.612 (vs saturated W2 0.029 at DOF 0; a-priori 0.84). => structure reproduces SHOOT translocation under a constrained fit; GRAIN does not (needs per-congener phloem loading).
