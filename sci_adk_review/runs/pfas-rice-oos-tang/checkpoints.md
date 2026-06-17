# Agent judgment checkpoints

proof/qualitative hypotheses awaiting an in-session agent verdict (no autonomous LLM call). For each, author verdicts/<hyp-id>.json with the chief-over-N trail, then re-enter the loop (sci-adk resolve <run-dir>).

## hyp-001 (qualitative)
- Criterion: Expert judgment based on evidence
- Finding: Out-of-sample vs the INDEPENDENT Tang 2026 dataset (theory/QSPR monotone f_xy, NOT fit to Tang): log10 RMSE 1.232 (~17x). Systematic miss -- PFSA way UNDER (PFOS stalk 0.013 vs 0.571, endosperm 0.004 vs 0.944; ~40-200x) and ether OVER (GenX stalk 8.3 vs 0.92, leaf 17.7 vs 1.33; ~10-13x). Only an IN-SAMPLE Tang-refit f_xy reaches 0.519. The model does NOT predict an independent dataset from another dataset's parameters: f_xy is dataset/condition-dependent (PFSA head-group offset, GenX ether offset).
