# Long-chain shoot delivery & the f_xy U-shape — exploration record

> Companion to `nstem_gradient_exploration.md`. Documents the investigation of the
> repo's central f_xy tension so it isn't re-discovered. Diagnostic scripts were
> run against Yamazaki 2023 with the MEASURED forcings (`forcing_rice` + ORYZA
> `growth_rice`); `model_api.simulate` now exposes `f_xy_override`, `L_Ph_override`,
> `kappa_d_override` for this kind of work.

## The tension (what we set out to resolve)
`parameters.json` carries two f_xy: the **monotone** `f_xy_recommended` (TSCF theory,
C4 0.79 → C12 0.003, labelled "physical") and `f_xy_W2fit` (a per-congener fit whose
long-chain values *rise*, labelled "non-physical artifact"). Only the W2 fit reproduces
the data (log10 RMSE 0.029). The repo's framing: the monotone is correct and the W2
rise is a single-compartment artifact.

## Diagnosis (the framing was wrong)
1. **`f_xy_recommended` gives RMSE 0.982** (measured forcing). The error is entirely in
   **straw** (root is fine) and is **two-sided**: short chains *over*-predicted, long
   chains *under*-predicted. Observed straw is **U-shaped** in chain length (high at C4,
   min ~C7–C9, high again at C10–C12); a monotone f_xy can only give a monotone straw.
2. **straw is leaf-dominated** in the model (stem contributes ~0), so this is not the
   within-stem vertical gradient that `nstem` addresses — it is whole-shoot delivery.
3. Inverting for "the f_xy each congener needs to match straw" gives a clear **U-shape**
   (C4 0.24 → C7 0.016 min → C11 0.20). The long-chain rise is **robust at PFDA/PFUnDA**,
   not just the near-MQL PFDoDA outlier. **So the W2 long-chain rise is real physics, and
   `f_xy_recommended`'s monotone decline is the over-extrapolation.**
4. The rise is **not** driven by `K_prot` (flat, 10–17 L/kg) but tracks **`K_PL`**
   (membrane/lipid partition, 42 → 10⁵ L/kg): long-chain PFAS are surfactants that
   partition into membranes and cross the endodermis via the lipid phase even though the
   free-anion electrodiffusion (TSCF) is shut down.

## What works
A **U-shaped f_xy driven by measured `K_PL`** — `f_xy(n) = a·e^(−b(n−4)) + c·K_PL/(K_PL+d)`,
5 GLOBAL params (not the saturated 33 of W2) — plus the existing per-congener `kappa_d`:

| model | RMSE (all 11) | RMSE (excl PFDoDA) | root | straw | grain |
|---|---|---|---|---|---|
| monotone `f_xy_recommended` | 0.982 | — | — | — | — |
| **U-shaped f_xy (5 global, K_PL)** | **0.370** | **0.286** | **0.121** | **0.179** | 0.461 |
| saturated W2 (33 per-congener) | 0.029 | — | — | — | — |

Root and straw are essentially solved with a handful of global, mechanistic, measured-
K_PL-driven params. **This overturns the "use `f_xy_recommended`" guidance** and converts
the "non-physical W2 rise" into a measured mechanism (lipid-facilitated translocation) —
a PFAS-specific term absent from the neutral/PPCP Brunetti DPU framework.

## The residual: grain, and the single root cause
`grain` stays under-predicted (RMSE ~0.46). It traces to the **same root cause** as the
earlier "root vs shoot" tension: the model loads only the **free aqueous fraction** into
xylem (`f_xy·Cw`) and phloem (`L_Ph·Cw`), and `Cw = C/B` collapses (~1/B) for high-binding
long chains — so translocation is bottlenecked exactly where the data says it shouldn't be.

We confirmed this structurally: adding a **B-independent "bound-loading" term** `g·C`
(the lipid-associated fraction rides the membrane phase, bypassing the free-conc
bottleneck) **fixes the long-chain grain** — PFDA grain 0.25 → 3.1 (obs 3.37), PFUnDA → 4.0
(obs 3.1). BUT a *flat* `g·C` over-feeds short/mid chains (their `C` is not small enough),
giving RMSE 0.37–0.49. The lipid loading must be **`K_PL`-gated** (≈0 for short chains, on
for long) — which is exactly why the K_PL-driven U-shaped f_xy (the phenomenological
version of the same physics) scores best (0.286).

## Conclusion / recommendation
- The correct mechanism is **lipid-bound, `K_PL`-gated loading on every membrane step**
  (root influx, xylem, phloem): short chains move as the free anion (low, declining TSCF),
  long chains ride the membrane/lipid fraction. This unifies the root, straw and grain
  long-chain behaviour and is the genuine PFAS-specific extension beyond Brunetti.
- A fully **global** mechanistic fit does **not** reach the saturated W2 (0.029) — expected,
  since W2 is over-parameterised — but cuts the monotone error ~3× (0.98 → ~0.29) with
  4–9 global params, and **solves root + straw outright**.
- **Honest limits:** grain remains the hardest tissue (free-conc bottleneck is sharpest
  there); `kappa_d` could not be globalised without breaking the steep long-chain root
  rise; **PFDoDA is a near-MQL outlier** that caps the achievable RMSE. All of this is
  still **in-sample** (Yamazaki) — it is the discovery of the correct *shape and mechanism*,
  **not validation**.
- **Decisive validation experiment:** measure f_xy directly as the xylem-sap-to-root-water
  ratio per congener (root-pressure exudate) across chain length — the U-shape and its
  `K_PL`-gated lipid arm are a sharp, falsifiable prediction.

## Out-of-sample test (Kim 2019 grain) — the first predictive signal
`validation/oos_crossdataset.py` transfers the Yamazaki-fit models to data they were NOT
fit on. Against **Kim 2019 brown-rice (grain) BAF** (porewater basis; PFOA excluded as it
was used in the L_Ph fit), the lipid mechanism predicts the independent chain-length pattern
— including the long-chain **rise** — far better than the alternatives:

| log10 RMSE vs Kim grain | lipid | monotone | W2 |
|---|---|---|---|
| excl. PFOA (all) | **0.55** | 2.04 | 1.11 |
| reliable only (PFHpA, PFNA; DF≥15%) | **0.23** | 1.91 | 1.43 |

The monotone and W2 models give ~0.05 for the long-chain grain (obs PFUnDA 33, PFDoDA 35);
the lipid model gives O(6) — under the (low-detection-frequency, unreliable) long chains but
the only model that captures the direction and order of magnitude. **This is the project's
first genuine out-of-sample predictive signal, and it specifically supports the
lipid-facilitated grain-loading mechanism.** Li 2025 is inconclusive — its reported BAFs scale
inversely with water quality (group-water denominator), so only the water-independent TF is
usable, and even that is root-surface-confounded; W2 wins straw/root, lipid wins grain/root,
all within factors of 2.5–30. Caveat: still not validation in full — Kim's long chains are
low-DF and lipid under-predicts them; a tissue-/time-resolved dataset (the xylem-/phloem-sap
experiment) remains the decisive test.

## Two-pool root (the root/shoot trade-off) — TESTED, set aside
Lipid loading fixes the long-chain grain but drains the long-chain **root** (mass balance:
the bound pool that feeds the shoot is the same pool whose burden is the root BAF). The
natural fix is a **two-pool root** — a mobile pool (uptake + loading) plus a slow,
sequestered apoplast/cell-wall/Fe-Mn-plaque pool that accumulates the root burden without
feeding the xylem. We built and fit it (5-state ODE: root_mobile, root_seq, stem, leaf,
grain; mass-conserving; `k_seq`/`k_rel` kinetic sequestration):

- **Uniform sequestration** over-inflates the *short*-chain root (PFBS root 2→24, PFOS 6→58):
  it sequesters everything. RMSE 0.67.
- **K_PL-gated sequestration** (only long chains park) still over-predicts mid-chain and PFSA
  root (PFDA 4→45, PFOS 6→108). RMSE 0.51 — both worse than the single-pool lipid model (0.36).

**Why it cannot work (decisive):** the observed root BAF does **not track `K_PL`**. PFOS (C8
PFSA) and PFUnDA (C11 PFCA) have **identical `K_PL` = 31623** yet roots of **5.93 vs 19.53**
(3.3×); and within the PFCAs the root rises ~4×/CF₂ while `K_PL` rises ~2×/CF₂. So no
`K_PL`-gated sequestration can give PFOS a low root and PFUnDA a high one simultaneously.
The very-long-chain root burden is governed by a **separate, chain-length-(not `K_PL`)-specific**
process (irreversible cell-wall/plaque binding, head-group-dependent) — plus PFDoDA is a
near-MQL outlier. Resolving it cleanly would need a per-congener root-sink descriptor, i.e.
over-parameterising the noisiest tail. **Set aside** (cf. the nstem hysteresis avenue). The
honest takeaway: the root sink and the shoot-delivery mechanism are **distinct**, so the
single-pool lipid term keeps the grain fix with a known, documented root trade-off.

## Status — WIRED (opt-in, default off)
The mechanism is now in the canonical core. `Compound` gained `g_xy`/`g_ph` (lipid-bound
loading conductances, default 0 → free-only model recovered exactly; mass conservation and
all tests unchanged). `RiceUptakeModel.rhs` loads `f_xy·Cw + g_xy·C` (xylem) and
`L_Ph·Cw + g_ph·C` (phloem). `model_api.simulate(lipid_loading=True)` switches on the
**`K_PL`-gated** parameterization via `model_api.lipid_loading_conductances(n_C, K_PL, group)`
(`LIPID_LOADING` constants, fit to Yamazaki excl. PFDoDA); `f_xy/L_Ph/kappa_d/g_xy/g_ph`
overrides are exposed for further work. Tests: `tests/test_model_api.py`
(`test_lipid_loading_*`).

**What it buys, honestly.** Turning it on cuts the monotone-`f_xy` error from 0.982 to
~0.36 (excl PFDoDA) and **fixes the long-chain grain** the free-only model could not reach
(PFDA grain 0.04 → ~5, obs 3.37; PFUnDA → ~6, obs 3.1). But it is a **trade-off, not a free
win**: loading the bound pool drains the long-chain **root** (PFUnDA root 2.3 vs obs 19.5),
because a single root pool with one partition cannot hold a high root burden *and* export
heavily to the shoot. The phenomenological `K_PL`-driven U-shaped `f_xy` still scores better
overall (0.286) by keeping the root and tuning `f_xy` up — at the cost of grain. The honest
reading: the bound-loading term is the **mechanistically correct** representation of
long-chain translocation, and it exposes a remaining **structural** limit (the root needs a
sequestered, non-draining sub-pool — e.g. apoplast/plaque `K_surf` or a two-pool root — to
carry the huge long-chain root BAF without starving the shoot). Still in-sample; the
xylem-sap experiment remains the decisive test.
