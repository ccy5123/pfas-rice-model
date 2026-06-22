# Background

The audit established that the free-anion transport model fails out-of-sample
(runs/pfas-rice-oos-tang) while the K_PL-gated lipid-facilitated loading mechanism
generalizes out-of-sample across two clean independent datasets
(runs/pfas-rice-oos-lipid, runs/pfas-rice-oos-multidataset). Separately, on the
in-sample Yamazaki series, the lipid mechanism cuts the whole-series log10 RMSE
from 1.035 (free-anion) to 0.386 (runs/pfas-rice-longchain). The lipid mechanism
is currently OPT-IN (default off). The long-chain follow-ups
(runs/pfas-rice-longchain-complete, -decouple) showed the remaining long-chain
root tradeoff is not yet cleanly resolvable in the single-pool prototype. This run
turns those findings into one engine-adjudicated MODEL-SELECTION verdict: across
all available measured evidence, which transport configuration is best-supported,
and should the lipid mechanism be the recommended configuration?

# Goal

Model selection across all measured evidence: is K_PL-gated lipid-facilitated loading the best-supported transport configuration, and should it be recommended over the free-anion default?

Decide, from the measured datasets (Yamazaki in-sample whole-series; Tang 2026 and
Kim 2019 out-of-sample), whether the lipid mechanism is the consistent best:

1. In-sample: does lipid beat the free-anion default on the whole-series Yamazaki
   fit by a decisive margin?
2. Out-of-sample (Tang): does lipid beat free-anion on the independent Tang dataset?
3. Out-of-sample (Kim): does lipid beat the monotone/free baseline on the
   independent Kim grain dataset?
4. Consistency: is lipid the winner on EVERY measured dataset (not just on average),
   i.e. is the minimum win-margin across all three still decisive?

# Method

Compose a frozen Spec of threshold hypotheses whose statistics are the verified
log10 RMSE values from the committed runs/experiments (in-sample whole-series from
validation/longchain_mechanism.py; OOS from runs/pfas-rice-oos-tang/-lipid and
runs/pfas-rice-oos-multidataset). Classify the evidence as measured (every number
is a comparison against real data). The engine resolves the numeric win-margins
autonomously and renders the paper; an agent-authored LaTeX-safe prose narrative
states the actionable recommendation. No new model is introduced; this is a
selection verdict over results already adjudicated in prior runs.

# Expected Output

An engine-rendered run (runs/pfas-rice-model-selection) whose claims state, factually,
whether the lipid mechanism is the consistent best-supported transport model across
all measured evidence, plus an explicit, honest recommendation: prefer lipid loading
for shoot/grain/out-of-sample prediction; keep the free-anion as the conservative
default in code (lipid stays opt-in) until a reliable 2-pool root removes the
single-pool long-chain root tradeoff. Reproducible under sci-adk verify.
