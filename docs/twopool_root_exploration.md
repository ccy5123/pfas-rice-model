# Two-pool root — decoupling the root sink from shoot delivery (exploration record)

> Companion to `fxy_longchain_lipid_exploration.md` and `nstem_gradient_exploration.md`.
> Records the BAF-prediction "고찰" session that tested whether splitting the root into a
> mobile (shoot-feeding) pool and a sequestered (burden-holding) pool resolves the central
> mass-balance tension. Script: `validation/twopool_root_exploration.py` (standalone,
> EXPLORATORY, in-sample Yamazaki 2023; the canonical core and `parameters.json` are
> UNCHANGED). Figure: `validation/figures/twopool_root_exploration.png`.

## The tension this addresses

The single root pool cannot simultaneously reproduce a **high long-chain root BAF** and a
**non-trivial long-chain shoot BAF**: the same pool whose burden is the root BAF is the pool
that loads the xylem, so any mechanism that delivers long chains to the shoot
(lipid-bound loading `g·C`) **drains the root** (PFUnDA root 19.5 → ~2.3 when lipid loading
is on; `fxy_longchain_lipid_exploration.md` §"Two-pool root"). The data require both:

| congener | root | straw | grain | straw/root | grain/root |
|---|--:|--:|--:|--:|--:|
| PFOA (C8) | 0.49 | 0.83 | 0.46 | 1.7 | 0.95 |
| PFUnDA (C11) | 19.53 | 8.16 | 3.13 | 0.42 | 0.16 |
| PFDoDA (C12) | 69.28 | 49.75 | 45.51 | 0.72 | 0.66 |
| PFOS (C8 PFSA) | 5.93 | 4.35 | 1.97 | 0.73 | 0.33 |

The very-long chains (PFDoDA) bind hardest (K_PL 66069) yet translocate to shoot with TF ≈
0.7 — a single reversible pool cannot do this.

## Why the sequestration descriptor must be NON-K_PL

The decisive fact: **PFOS (C8 PFSA) and PFUnDA (C11 PFCA) have identical K_PL = 31623 and
near-identical B_k_root (49.4 vs 49.1), yet observed root BAF is 5.93 vs 19.53 (3.3×).**
No K_PL-gated sink can separate them (already confirmed for the K_PL-gated two-pool in the
companion doc). Expressed as the fraction of the binding ceiling B_root each congener reaches:

```
PFBA 0.57  PFHxA 0.18  PFOA 0.12  PFNA 0.06  PFDA 0.17  PFUnDA 0.40  PFDoDA 0.69
PFBS 0.95              PFOS 0.12
```

This is U-shaped in chain length and head-group-split — a real chain-length·head-group
process, not a K_PL one. The hypothesis: an **irreversible apoplast / cell-wall / Fe–Mn-plaque
sink** with a rate `k_seq(n, head_group)` independent of K_PL.

## The model (`twopool_root_exploration.py`)

5-state ODE `[root_mobile, root_seq, stem, leaf, grain]`, mass-conserving, sole source
`M_root·j_R` into the mobile pool, identical demo forcings to `reproduce_demo`:

- **root_mobile**: binding `B_m` (= the measured basis-A B_root); GHK+carrier uptake; loads
  the xylem with the **monotone physical** `f_xy_recommended` free term **+** a K_PL-gated
  lipid-bound term `g_xy·C_m` (mirrors `model_api.lipid_loading_conductances`).
- **root_seq**: receives `k_seq·C_m`, released only by growth dilution → a **terminal
  accumulator** (like leaf/grain). Its final burden ≈ ∫k_seq·C_m so high-`k_seq` chains
  accumulate a large root BAF *without* draining the mobile pool's shoot feed.
- `k_seq(n, group) = 10^(ks0 + ks_b·(n−8) + ks_sa·[PFSA])` — the **non-K_PL** descriptor.

Total root BAF = `C_mobile + C_seq`. 7 global params fit to Yamazaki (11 congeners × 3
tissues, PFDoDA excluded from the fit as a near-MQL outlier).

## Result 1 — global fit (the structure lets us keep the physical f_xy)

| model | params | RMSE (excl PFDoDA) | RMSE (all) | uses physical monotone f_xy? |
|---|--:|--:|--:|:--:|
| monotone `f_xy_recommended` (single pool) | — | — | 0.982 | yes |
| U-shaped K_PL-driven f_xy (single pool) | 5 global | 0.286 | 0.370 | **no** (fitted U-shape) |
| **two-pool root (this work)** | 7 global | **0.257** | 0.305 | **yes** |
| saturated W2 (per-congener f_xy/L_Ph/κ_d) | 33 | — | 0.029 | no |

By tissue (all 11): root 0.366, straw 0.237, grain 0.296.

**The headline result:** the two-pool structure matches the best previous global model
(U-shaped K_PL f_xy, 0.286) **while keeping the monotone, theory-consistent
`f_xy_recommended`** — the straw U-shape now emerges from lipid loading + the root
decoupling, so we no longer need the "non-physical" U-shaped f_xy. This is a genuine
structural improvement: the shoot pattern and the root burden are produced by *separate,
physically-motivated* mechanisms.

**But the smooth global k_seq collapses (`ks_b ≈ 0`).** The optimizer drove the chain-length
slope to zero — k_seq reverted to ~constant — so the chain-length root rise is carried by
`B_m` (≈ K_PL) again, and the model **inherits B's inability to separate PFOS from PFUnDA**:
PFOS root 16.1 (obs 5.9, over) vs PFUnDA 9.5 (obs 19.5, under). Excluding PFDoDA (the
strongest long-chain anchor) from the fit removed the signal that would pull `ks_b` up.

## Result 2 — root-matched k_seq (the sufficiency test + the descriptor's true shape)

Holding the global *shoot* params fixed, back out per congener the `k_seq` that makes model
root **exactly** equal the observed root (1-D solve), then read off the resulting straw/grain.
This separates two questions: *is the structure sufficient?* and *what shape must the
descriptor be?*

| congener | n_C | grp | empirical k_seq | frac in seq pool | straw p/o | grain p/o |
|---|--:|:--:|--:|--:|--:|--:|
| PFBA | 4 | FCA | 0.289 | 0.95 | 6.5 / 11.0 | 0.72 / 1.06 |
| PFPeA | 5 | FCA | 0.090 | 0.85 | 4.0 / 2.74 | 0.52 / 0.26 |
| PFHxA | 6 | FCA | 0.072 | 0.83 | 1.9 / 1.14 | 0.36 / 0.42 |
| PFHpA | 7 | FCA | 0.153 | 0.91 | 0.94 / 1.07 | 0.30 / 0.42 |
| PFOA | 8 | FCA | 0.043 | 0.74 | 0.64 / 0.83 | 0.36 / 0.46 |
| PFNA | 9 | FCA | **0.014** | 0.47 | 1.25 / 0.69 | 1.03 / 0.41 |
| PFDA | 10 | FCA | 0.074 | 0.83 | 3.4 / 3.46 | 3.16 / 3.37 |
| PFUnDA | 11 | FCA | **0.210** | 0.93 | 6.73 / 8.16 | 6.08 / 3.13 |
| PFDoDA | 12 | FCA | **0.490** | 0.97 | 11.0 / 49.75 | 9.07 / 45.51 |
| PFBS | 4 | FSA | 0.466 | 0.97 | 2.30 / 2.18 | 0.24 / 0.25 |
| PFOS | 8 | FSA | **0.047** | 0.75 | 2.42 / 4.35 | 0.81 / 1.97 |

**With root MATCHED:  straw RMSE = 0.255, grain RMSE = 0.307** — essentially unchanged from
the global fit (0.237 / 0.296).

Three conclusions:

1. **The two-pool STRUCTURE is sufficient (the §3 wall is breakable).** Forcing the root to the
   observed value — including PFUnDA 19.5 and PFDoDA 69.3 — does **not** degrade the shoot. The
   seq pool carries 47–97 % of the root burden while the mobile pool keeps feeding the xylem.
   You CAN hold a high long-chain root BAF *and* deliver to the shoot; the mass-balance coupling
   that defeated the single-pool lipid model is removed.

2. **The required root-sink descriptor is real and genuinely non-K_PL.** At identical
   K_PL = 31623, **PFUnDA demands k_seq = 0.210 but PFOS only 0.047 (4.5×)** — the sequestration
   the data require separates the two exactly where K_PL/B cannot. This quantifies the
   chain-length·head-group-specific irreversible root sink hypothesised above.

3. **The descriptor is U-shaped, which is why the monotone global fit collapsed.** The empirical
   k_seq is **not** monotone in chain length — it falls PFBA 0.29 → PFNA 0.014 (min) then rises
   to PFDoDA 0.49 (and PFBS 0.47 ≫ PFOS 0.047 for PFSA), mirroring the U-shaped "fraction of the
   binding ceiling reached". My global form `ks0 + ks_b·(n−8)` was **monotone**, so it structurally
   cannot represent a U and the optimizer correctly drove `ks_b → 0`. **The concrete fix is a
   U-shaped k_seq** (same functional family as the straw `f_xy` U-shape), anchored on these
   empirical values.

> Caveat: PFBA's high k_seq (0.29, 95 % sequestered) is partly an inverse of its high
> `f_xy_recommended = 0.79` draining the mobile root — the *effective* root retention is what is
> U-shaped; k_seq and f_xy are coupled through the mobile-pool balance. PFDoDA shoot stays
> under-predicted (near-MQL outlier; its straw 49.7 is the single largest residual).

## Result 3 — U-shaped `k_seq(n)` realized (the well-posed follow-up, DONE)

The fix flagged by Result 2 is now implemented: an **asymmetric U in chain length** (a declining
short-chain arm + a rising long-chain arm), with the rising arm in **`n`, NOT `K_PL`**, fitted to
the root-matched empirical `k_seq` and plugged back into the full ODE:

```
k_seq(n,grp) = [0.268·e^(−0.52(n−4)) + 0.615·e^(1.35(n−12))] · {10^(+0.18) if PFSA}
```

5 descriptor params (fit to the 11 root-matched values, descriptor-fit log10 RMSE 0.199), with
all global *shoot* params held from Result 1.

| | full-ODE log10 RMSE (all 11, **incl PFDoDA**) | root | straw | grain |
|---|--:|--:|--:|--:|
| linear global `k_seq` (Result 1, `ks_b→0`) | 0.305 | 0.366 | 0.237 | 0.296 |
| **U-shaped `k_seq(n)` (this step)** | **0.251** | **0.156** | 0.260 | 0.311 |

- **The PFOS/PFUnDA separation is realized.** At identical K_PL = 31623 the form gives
  **PFOS (C8) k_seq = 0.054 vs PFUnDA (C11) k_seq = 0.166 (3.1×)** — purely from chain length
  (8 vs 11) plus the head-group offset — so model root is now PFOS 6.61 (obs 5.93) and PFUnDA
  15.86 (obs 19.53), where the linear fit had them backwards (16.1 / 9.5). The non-K_PL signature
  the data demand is reproduced by a smooth descriptor.
- **Root is essentially solved (0.156), now including PFDoDA** (root 82.3 vs obs 69.3 — a mild
  over, vs the linear fit's 17.7 *under*). The rising long-chain arm carries the C10–C12 root rise
  that `B_m`/K_PL could not.
- **Overall 0.251 (all 11) beats the linear two-pool (0.305) and the single-pool U-shaped-K_PL-f_xy
  (0.370 all / 0.286 excl)** — while keeping the **monotone physical `f_xy_recommended`**.
- **The residual is now entirely the very-long-chain SHOOT.** PFDoDA straw 10.45 vs 49.75 and
  grain 8.68 vs 45.51 are the dominant errors; every other tissue is within a factor ~2.5. This is
  the same C12 shoot-delivery floor seen everywhere in the repo (near-MQL PFDoDA + the active-
  carrier capacity limit, cf. PR #21 `runs/pfas-rice-longchain` LC4–LC6) — a *shoot* problem, not a
  root one, so it is outside what `k_seq` (a root term) can fix.

Figure panel (b) plots the U-shaped-`k_seq` predictions (RMSE 0.251); panel (a) overlays the
fitted U-form on the root-matched empirical `k_seq`.

## Honest status

- **What works:** the two-pool *structure* decouples root burden from shoot delivery, lets the
  model keep the monotone physical `f_xy_recommended`, and — with the **U-shaped `k_seq(n)`**
  (Result 3) — reaches log10 RMSE **0.251 across all 11 congeners including PFDoDA** (root 0.156),
  reproducing the non-K_PL **PFOS/PFUnDA separation** (3.1×) that B/K_PL cannot. The root-matched
  test proves the structure is **sufficient**; the U-form realizes it with a parsimonious smooth
  descriptor. Mass-conserving; the canonical core and `parameters.json` are untouched.
- **The remaining residual is the very-long-chain SHOOT, not the root.** PFDoDA straw/grain stay
  ~3–5× under (near-MQL outlier + the active-carrier capacity limit on C12 shoot delivery). This
  is a *shoot*-delivery problem `k_seq` (a root term) cannot address, and is consistent with every
  other long-chain finding in the repo. Promoting the U-shaped `k_seq` into `parameters.json` (it
  is currently exploration-only) would need an out-of-sample check first.
- **Mechanistic reading:** the model now says two *distinct* processes set the root BAF — a
  reversible membrane/lipid partition (`B_m`, ~K_PL, feeds the shoot) and an irreversible
  chain·head-group-specific sink (`k_seq`, non-K_PL, holds the burden). Long-chain PFCAs load the
  sink heavily (PFUnDA/PFDoDA), PFSAs much less at matched chain (PFOS), and short chains barely
  at all in absolute terms.
- **Decisive experiment (unchanged):** per-congener xylem-sap / root-water ratio (direct f_xy)
  **and** a desorption-resistant root-fraction assay (isolating the irreversible `k_seq` pool)
  across chain length and head group. Still in-sample (Yamazaki only); this is mechanism
  discovery and a well-posed next fit, **not** validation.
