# Background

The audit concluded the model was not yet usable for dietary risk assessment: the
risk-relevant compartment (brown-rice grain) was under-predicted (H4 REFUTED), and
the longest chains could not be reproduced (the single-pool refit hit ceilings,
~4-6x under). Two developments change this. (1) BREAKTHROUGH
(validation/longchain_closure.py): the long-chain root<->shoot was NOT structurally
unresolvable -- that was an artifact of holding f_xy fixed and only adding the
non-subtractable lipid term. Recognizing two INDEPENDENT physical facts about a
long chain -- a LOW xylem-loading f_xy (strong root retention) and an ENHANCED
active carrier (high uptake) -- the standard 2-pool reproduces C10-C12 root, straw
and grain at log10 RMSE ~0.08 (saturated, DOF 0 = structural adequacy, not a-priori
prediction). (2) The K_PL-gated lipid mechanism predicts the dietary compartment
OUT-OF-SAMPLE: independent Kim 2019 brown-rice grain BAF at log10 RMSE 0.48 (factor
~3), reliable subset 0.20 (factor ~1.6). This run asks whether, together, these make
the model usable as a dietary risk-assessment tool, and at what level.

# Goal

Is the PFAS-rice model usable as a dietary risk-assessment tool, and at what assurance level? A risk-assessment-readiness verdict over grain (the dietary compartment).

Decide, against measured data, whether the model can serve dietary risk assessment
for PFAS in brown rice:

1. Structural coverage: with the 2-pool + free f_xy + active carrier, is there any
   remaining structural blind spot across the C4-C12 congener range (including the
   long chains)?
2. Dietary-compartment prediction: is brown-rice grain predicted OUT-OF-SAMPLE on an
   independent dataset within a screening-adequate factor (~3)?
3. Reliable subset: is the reliable-detection grain prediction within ~factor 1.6?
4. Honest assurance level: is the worst-case grain residual large enough that the
   tool is SCREENING-level (bounded uncertainty), not regulatory-precision?

# Method

Run the live breakthrough experiment (validation/longchain_closure.py) for the
long-chain structural-coverage statistic, and take the grain out-of-sample statistics
from the committed cross-dataset runs (runs/pfas-rice-oos-multidataset for Kim grain;
runs/pfas-rice-oos-lipid for the worst-case PFOS-endosperm residual). Freeze these as
threshold hypotheses; the engine resolves them and renders the paper. All evidence is
measured (vs Yamazaki / Kim / Tang). The verdict is explicitly an assurance-LEVEL
statement, not a binary -- screening vs regulatory -- so the honest worst-case bound
is a first-class hypothesis, not a footnote.

# Expected Output

An engine-rendered run (runs/pfas-rice-risk-readiness) whose claims state whether the
model is usable as a SCREENING-level dietary risk-assessment tool: grain predicted
out-of-sample within ~factor 3 (reliable ~1.6), the long-chain structural blind spot
closed, and the worst-case uncertainty bounded (~5x for PFSA endosperm / ether) so the
tool is screening-grade not regulatory-grade. Reproducible under sci-adk verify, with
the congener-specific uncertainty factors recorded so a risk assessor can use them.
