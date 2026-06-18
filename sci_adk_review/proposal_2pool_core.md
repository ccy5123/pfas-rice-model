# Background

The breakthrough (validation/longchain_closure.py) showed the long chains ARE
structurally closable: a 2-pool root (mobile + bound) with a LOW per-congener f_xy
(strong root retention) and an ENHANCED active carrier (high uptake) reproduces
C10-C12 root, straw and grain at log10 RMSE ~0.08, where the single-pool core could
not (refit_oryza hit ceilings ~4-6x under). For the model to be usable -- in
particular as the long-chain-capable configuration for dietary risk screening -- the
breakthrough must be a proper, reusable model component, not a validation script. It
has now been PROMOTED to src/pfas_rice_two_pool.py with a clean API and a model_api
hook (simulate_two_pool / close_longchain_2pool), additive to (and leaving unchanged)
the canonical 4pool_surf core. This run verifies that the promoted core component
reproduces the breakthrough.

# Goal

Does the promoted src/ two-pool component reproduce the long-chain breakthrough? Verifying the wired core component, and confirming the two independent physical levers.

Verify, against the measured Yamazaki long chains, that the wired src module:

1. reproduces the long-chain closure (root+straw+grain log10 RMSE at the breakthrough level);
2. requires a LOW f_xy for every long chain (the strong-root-retention lever);
3. requires an ENHANCED active carrier for the longest chain (the high-uptake lever),
   confirming the two levers are independent and both necessary.

# Method

Import src/pfas_rice_two_pool.py and run its close_longchain (the saturated 3-param
structural-adequacy fit: free f_xy -> straw, active carrier -> root, g_ph -> grain)
for PFDA, PFUnDA and PFDoDA against the Yamazaki BAFs. Compute the long-chain
root+straw+grain log10 RMSE and read the fitted f_xy and carrier multipliers. Freeze
these as threshold hypotheses; the sci-adk engine resolves them and renders the paper.
All evidence is measured (vs Yamazaki). The canonical core is unchanged; this verifies
the additive component, and is structural ADEQUACY (reproduction, DOF 0), not a-priori
prediction.

# Expected Output

An engine-rendered run (runs/pfas-rice-2pool-core) confirming the promoted src
component reproduces the breakthrough (long-chain RMSE at the ~0.08 level) and that the
closure needs both a low f_xy (root retention) and an enhanced carrier (uptake) --
making the long-chain-capable model a reusable component for risk screening.
Reproducible under sci-adk verify.
