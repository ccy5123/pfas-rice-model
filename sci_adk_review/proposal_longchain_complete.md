# Background

The long-chain sub-investigation (runs/pfas-rice-longchain, LC1-LC6) found that
free-anion loading structurally starves the long-chain (C10-C12) shoot (LC1), a
B-independent lipid-facilitated loading term closes most of the gap (LC2), a
single pool does so only by draining the root (LC3), a 2-pool root closes C10/C11
but not PFDoDA (LC4), membrane conductance cannot lift PFDoDA (LC5a) while an
enhanced active carrier can reach its root and grain (LC5b), and the carrier
enhancement is not QSPR-able from chain length (LC6). From these, FINDINGS sec.7
PROPOSED a "complete long-chain resolution = 2-pool (free + lipid-bound) root +
lipid-facilitated loading + enhanced long-chain active carrier." But each lever
was only ever tested in isolation; the proposed combination was never run as one
model and tested for the property that matters -- reproducing the long-chain ROOT
and SHOOT simultaneously. This run tests that proposal.

# Goal

Does the proposed complete long-chain resolution close the long chains? Testing 2-pool + lipid + carrier as one model for simultaneous root and shoot reproduction.

The goal is to combine all three levers into one model and adjudicate, against the
measured Yamazaki BAFs across C8-C12, whether the complete recipe reproduces the
long-chain (nC>=10) ROOT and SHOOT (straw) and grain SIMULTANEOUSLY:

1. Root: with the LC6 root-matching active-carrier multiplier, does the long-chain
   root close?
2. Grain: with the fitted lipid phloem term, does the long-chain grain close?
3. Shoot: does the straw close at the same time, or does the carrier that fixes the
   root over-feed the shoot?
4. Diagnosis: quantify any straw over-feeding for the longest chains (C11-C12).

# Method

Reuse the committed LC 2-pool prototype (validation/twopool_longchain.py) without
changing it. For each congener, find the active-carrier multiplier (Vmax_in) that
reproduces the measured root (the LC6 root-matching multiplier), then, WITH that
carrier, fit the lipid-facilitated loading g_xy (to straw) and g_ph (to grain) on
the 2-pool root. Run validation/longchain_complete.py and report root/straw/grain
pred vs observed (Yamazaki) and the long-chain (nC>=10) log10 RMSE per tissue, plus
the C11-C12 predicted/observed straw ratio. Freeze the resulting statistics as
threshold hypotheses; the sci-adk engine resolves them autonomously and renders the
paper. All evidence is measured (vs Yamazaki). No core model is changed; this is a
prototype synthesis test of an already-recorded mechanism set.

# Expected Output

An engine-rendered run (runs/pfas-rice-longchain-complete) whose claims state,
factually, whether the complete recipe closes the long-chain root, grain and shoot
simultaneously -- with the honest verdict (whatever it is) recorded and reproducible
under sci-adk verify. This refines the FINDINGS sec.7 "complete resolution" claim
from a proposal into a tested result and, if simultaneous closure fails, names the
precise remaining open problem.
