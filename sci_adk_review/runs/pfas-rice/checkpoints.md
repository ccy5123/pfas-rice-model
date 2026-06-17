# Agent judgment checkpoints

proof/qualitative hypotheses awaiting an in-session agent verdict (no autonomous LLM call). For each, author verdicts/<hyp-id>.json with the chief-over-N trail, then re-enter the loop (sci-adk resolve <run-dir>).

## hyp-anion (qualitative)
- Criterion: GHK electrodiffusion with measured E_m=-120 mV and z=-1 yields an anion-exclusion factor e^N of order 1e2 (the IOC-formulation signature) => support; no exclusion => refute
- Finding: N = z*E_m*F/RT = +4.67 at E_m=-120 mV, z=-1, T=298.15 K => e^N = 106.8 (test_plant_module asserts e^N ~ 107); the GHK term in root_uptake() excludes the anion accordingly.

## hyp-yamazaki (qualitative)
- Criterion: the Yamazaki agreement is OUT-OF-SAMPLE predictive validation (independent data not used to fit the per-congener transport parameters) => support; an in-sample / saturated reproduction => refute
- Finding: reproduce_demo.py: W2 transport fit reproduces Yamazaki root/straw/grain BAFs across 11 congeners at log10 RMSE 0.029 (e.g. PFOA root 0.49/0.49, straw 0.83/0.83 pred/obs). The fit assigns ~3 transport params per congener against 3 tissue observations => SATURATED: reproduction is structurally guaranteed and is NOT out-of-sample prediction (CLAUDE.md §6).

## hyp-grain (qualitative)
- Criterion: model grain (brown-rice) BAF/TF matches measured within a factor adequate for dietary risk assessment => support; systematic structural under/over-prediction => refute
- Finding: Tang 2026 OOS: PFOA endosperm TF model 0.11 vs measured 0.95 (dw); grain structurally ~3-8x under across congeners; not closable by L_Ph / lipid tuning (docs/tang2026_grain_units_exploration.md).
Kim 2019: grain BAF reproduced only by FORCING L_Ph (0.07 -> 4.43, L_Ph~0.84) as a single-point in-sample anchor; Kim is grain-only, so root/straw TF are unconstrained.

## hyp-smiles (qualitative)
- Criterion: the SMILES front-end reproduces the curated measured-parameter model for a KNOWN PFAS via read-across (structure -> same Compound) => support; a mismatch => refute
- Finding: tests/test_pfas_structure.py: 23 passed (RDKit). A canonical SMILES that matches a curated congener rebuilds the SAME Compound from params/parameters.json (measured read-across). CAVEAT: for NOVEL structures f_xy is provisional (QSPR/interpolated), NOT validated — this claim is scoped to known structures only.
