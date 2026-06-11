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

## Status
Recorded, not yet wired into the canonical core. Wiring the `K_PL`-gated lipid loading into
`pfas_rice_plant_module_4pool_surf` (as an opt-in term, default off) is the natural next
implementation step; `model_api` already carries the override hooks used here.
