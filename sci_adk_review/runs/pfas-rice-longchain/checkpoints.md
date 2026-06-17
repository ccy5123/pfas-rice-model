# Agent judgment checkpoints

proof/qualitative hypotheses awaiting an in-session agent verdict (no autonomous LLM call). For each, author verdicts/<hyp-id>.json with the chief-over-N trail, then re-enter the loop (sci-adk resolve <run-dir>).

## hyp-lc-freethrottle (qualitative)
- Criterion: free-anion xylem/phloem loading structurally under-predicts long-chain (C10-C12) shoot accumulation -- even at the f_xy=1 ceiling the flux is throttled by Cw=C/B collapse => support; long chains are reachable by free loading => refute
- Finding: validation/longchain_mechanism.py (ORYZA2000 biomass): free-only (monotone f_xy) long-chain (nC>=10) straw+grain log10 RMSE 2.026 (~100x); PFDA straw 0.08 vs 3.46, PFDoDA straw 0.33 vs 49.8 / grain 0.10 vs 45.5. refit_oryza.py: even f_xy=1/L_Ph=1 ceilings leave PFDoDA straw 14.6 vs 49.8 -- the Cw=C/B collapse throttles free loading regardless of f_xy.

## hyp-lc-lipidfix (qualitative)
- Criterion: a B-independent lipid-facilitated bound-loading term (g_xy*C, g_ph*C) closes most of the long-chain shoot gap that free loading cannot => support; no improvement => refute
- Finding: validation/longchain_mechanism.py: the B-independent lipid term (g_xy*C, g_ph*C; K_PL-gated) cuts long-chain straw+grain log10 RMSE 2.026 -> 0.428 (~100x -> ~2.7x; whole series 1.035 -> 0.386). PFDA straw 0.08->5.95 (obs 3.46), PFUnDA straw 0.16->11.15 (obs 8.16), PFOS straw 0.17->5.17 (obs 4.35).
{"acquired": [], "failed": [{"doi": "10.1021/acs.est.4c06734", "error": "no OA PDF"}, {"doi": "10.1021/acs.est.5c11716", "error": "no OA PDF"}, {"doi": "10.1021/acs.est.7b06128", "error": "no OA PDF"}, {"doi": "10.1021/acsestengg.4c00107", "error": "no OA PDF"}, {"doi": "10.48130/newcontam-0025-0007", "error": "no OA PDF"}, {"doi": "10.1007/s40726-020-00168-y", "error": "no OA PDF"}, {"doi": "10.1139/er-2025-0116", "error": "no OA PDF"}], "counts": {"succeeded": 0, "failed": 7}, "normalization": {"pdfs": [], "counts": {"normalized": 0, "already_extractable": 0, "locked": 0, "error": 0}}, "citation_keys": {}, "citation_key_collisions": []}

## hyp-lc-nocost (qualitative)
- Criterion: the lipid mechanism reproduces the long-chain shoot WITHOUT degrading the root or the short/mid-chain fits (a clean single-mechanism fix) => support; it trades off root / leaves a residual => refute
- Finding: validation/longchain_mechanism.py: the single-pool lipid term TRADES OFF root for long chains -- PFUnDA root 20.6->3.9 (obs 19.5), PFDoDA root 159->4.4 (obs 69); and PFDoDA shoot is STILL ~3-4x under (straw 14.7 vs 49.8). A clean cost-free fix is refuted: it needs a 2-pool (free + bound) split and a residual long-chain mechanism.

## hyp-lc-twopool (qualitative)
- Criterion: a 2-pool root (mobile water+protein pool feeding the xylem + a slow-exchanging lipid/cell-wall bound store holding the measured root burden) reproduces the long-chain shoot WHILE keeping the measured root high, closing the LC3 single-pool tradeoff => support; it cannot match root and shoot simultaneously => refute
- Finding: validation/twopool_longchain.py (ORYZA2000): the 2-pool root (mobile + slow lipid/cell-wall bound store; lipid loads from the mobile pool) matches MID-LONG chains root AND shoot SIMULTANEOUSLY -- PFDA(C10) root 3.5 vs 4.2, straw 5.0 vs 3.5, grain 4.1 vs 3.4; PFUnDA(C11) root 9.8 vs 19.5 (~2x), straw 13.4 vs 8.2, grain 3.6 vs 3.1. The single pool could not (LC3 root drain). Shoot RMSE 0.32.
validation/twopool_longchain.py: PFDoDA(C12) FAILS -- root collapses to 1.2 vs 69.3 (mobile pool rm=0.02 starves -> bound store rb only 1.2; straw 8.6 vs 49.7, grain 26 vs 45.5). The mobile-pool soil uptake (jR) is too small to sustain a high bound root AND feed the shoot: a mass-balance/UPTAKE limit, not internal distribution. The 2-pool does not close the longest chain.

## hyp-lc-cond (qualitative)
- Criterion: enhancing the root membrane CONDUCTANCE (kappa_d) closes the PFDoDA long-chain residual => support; conductance has no effect (the GHK anion-exclusion ceiling caps Cw at Cwo/e^N regardless of kappa_d) => refute
- Finding: validation/twopool_longchain.py kappa_d scan (PFDoDA): increasing the membrane conductance 1x->5000x leaves root at 1.0-1.2 vs 69 (unchanged). The GHK term is capped by anion exclusion -- Cw_m -> Cwo/e^N (e^N~107) regardless of kappa_d -- so conductance is NOT the lever.

## hyp-lc-carrier (qualitative)
- Criterion: enhancing the long-chain ACTIVE CARRIER capacity (Vmax) closes the PFDoDA root and grain (the carrier overcomes the anion-exclusion ceiling) => support; the carrier cannot reach it => refute
- Finding: validation/twopool_longchain.py carrier scan (PFDoDA): at Vmax_in ~5x base (20->100) the 2-pool reaches root 62 vs 69 and grain 46 vs 45.5 (straw 102 vs 50, ~2x over). The ACTIVE carrier overcomes the anion-exclusion ceiling the GHK pathway cannot -- the longest-chain residual is a carrier-capacity limit, closable. Consistent with the literature (active carrier-mediated root uptake).
