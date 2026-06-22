# Background

The long-chain "complete resolution" test (runs/pfas-rice-longchain-complete)
showed that combining the three proposed levers (2-pool root + lipid loading +
LC6 root-matching active carrier) into one model closes the long-chain root and
grain but OVER-feeds the shoot: the active carrier that fixes the root enlarges
the mobile root pool, whose free xylem loading f_xy*Cw_m alone exceeds the
observed straw (PFUnDA ~3.3x, PFDoDA ~2.3x over), and the lipid term (g_xy>=0)
cannot subtract. FINDINGS sec.7 named the missing piece a root->shoot DECOUPLING:
an irreversible root store that holds the measured root burden WITHOUT feeding the
xylem in proportion. This run tests the simplest such lever.

# Goal

Does an irreversible root sequestration decouple the long-chain root from the shoot? Testing asymmetric bound-store kinetics for simultaneous root and shoot closure.

The lever: make the 2-pool bound store irreversible by an asymmetric-kinetics
factor seq, k_on = ratio*k_off*seq (seq=1 = the equilibrium 2-pool; seq>1 traps
more burden in the non-translocating bound pool rb). At a fixed root-matching
carrier, scan seq and ask, against the measured Yamazaki BAFs for the long chains
(C11-C12):

1. Does increasing sequestration relieve the shoot over-feed while keeping the root
   at the measured level (a clean simultaneous closure within a factor 2)?
2. Or does it merely trade one tissue for the other?
3. Diagnose the coupling that governs the trade.

# Method

Reuse the committed LC 2-pool prototype (validation/twopool_longchain.py) without
changing it; add only the seq lever in a separate experiment
(validation/longchain_decouple.py). Per congener, fix the active carrier at its
seq=1 root-matching value, then scan seq. The fittable straw is bounded below by
the g_xy=0 floor (lipid g_xy>=0 can only add), so the straw fit ratio is
max(straw_floor/obs, 1) and no lipid fit is needed for the root<->shoot tension.
Report, per seq, the root ratio, the straw fit ratio, and the simultaneity gap
max(|log10 root/obs|, log10 straw_fit_ratio); summarize the baseline (seq=1, which
IS the complete recipe), the root inflation across the scan, and the best (minimum)
simultaneity gap. Freeze these statistics as threshold hypotheses; the engine
resolves them and renders the paper. All evidence is measured (vs Yamazaki). No
core model is changed.

# Expected Output

An engine-rendered run (runs/pfas-rice-longchain-decouple) whose claims state,
factually, whether the irreversible-sequestration lever achieves a clean
simultaneous root+shoot closure for the long chains, or only trades root for shoot
-- with the honest verdict recorded and reproducible under sci-adk verify. If it
fails, the run names the sharper open target (break the uptake<->mobile-conc
coupling, not just the bound-store kinetics).
