# Background

This repository's four-compartment dynamic uptake model for PFAS in rice
(Oryza sativa) has been subjected to an adversarial rigor audit with
sci-adk across seven independent runs: the main audit (pfas-rice), a
synthetic-data trap (pfas-rice-trap), a long-chain mechanism
sub-investigation (pfas-rice-longchain), a carrier-QSPR test
(pfas-rice-carrier), and three cross-dataset out-of-sample runs
(pfas-rice-oos-tang, pfas-rice-oos-lipid, pfas-rice-oos-multidataset).
sci-adk itself was created in response to a failure of this project: a
run on an empirical proposal once used synthetic data and the harness
reported "4/4 SUPPORTED" (the rice-failure defect named in
core/validity.py). The empirical investigation is now complete, but its
findings are scattered across seven run records and a Korean narrative
(FINDINGS.md). The per-run papers were rendered by sci-adk but their
non-ASCII bodies did not typeset, and there is no single citable
synthesis. This run consolidates the completed investigation into one
engine-rendered paper whose every headline claim is resolved by the
engine from the verified sub-run statistics, so the consolidation is
itself a reproducible sci-adk artifact rather than a hand-authored
document.

# Goal

Consolidated rigor audit of a four-compartment PFAS rice uptake model: foundations reproduce, naive out-of-sample prediction fails, and a lipid-facilitated loading mechanism generalizes across independent datasets.

The goal is a single, citable, engine-rendered synthesis of the seven
sci-adk runs that states and resolves the cross-run findings:

1. Reproducibility: every committed sub-run re-derives from its frozen
   record under sci-adk verify (exit 0), so the audit record is
   trustworthy.
2. Naive out-of-sample prediction fails: the free-anion model driven by
   theory/QSPR transport parameters not fit to the target predicts an
   independent dataset (Tang 2026) far worse than an in-sample refit.
3. The lipid-facilitated loading mechanism generalizes: a mechanism fit
   on Yamazaki long chains, never fit to the target, restores the
   independent-dataset prediction to the in-sample level.
4. The generalization is robust across a second clean independent
   dataset (Kim 2019 grain), not a single-dataset artifact.
5. Structural adequacy: under a constrained (degrees-of-freedom > 0) fit
   on the mechanistic ORYZA2000 biomass, the structure reproduces shoot
   (straw) translocation across the C4-C12 series.

# Method

Author the consolidation as a frozen Spec of threshold hypotheses whose
statistics are the verified outputs of the seven sub-runs (record
digests and log10 RMSE values), classified honestly: generated for the
reproducibility meta-claim, measured for the out-of-sample and
structural-adequacy claims. The sci-adk engine compiles the Spec,
persists the classified Evidence, resolves the numeric threshold claims
autonomously, and renders the paper draft. An agent-authored,
LaTeX-safe narrative (abstract, introduction, discussion) is supplied as
prose input and injected verbatim into the engine-rendered draft (this
is input, not autonomous generation). sci-adk verify then re-derives
every claim from the record with no LLM call. The consolidation
introduces no new experiments: every number traces to a per-run record,
and the per-run runs remain the authoritative source for their own
claims.

# Expected Output

An engine-rendered consolidated paper (runs/pfas-rice-consolidation/
paper/draft.tex) that ties the seven runs into one narrative arc with a
master ledger summarized in the discussion, centralized honest caveats,
and a verified record-digest table; plus a sci-adk verify reproduction
(exit 0, all claims reproduced) of the consolidation run. This paper
supersedes the hand-authored consolidation and the non-rendering per-run
skeleton drafts as the single readable, citable synthesis.
