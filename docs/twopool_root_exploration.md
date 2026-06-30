# Two-pool root — decoupling the root sink from shoot delivery (exploration record)

> Companion to `fxy_longchain_lipid_exploration.md` and `nstem_gradient_exploration.md`.
> Records the BAF-prediction "고찰" session that tested whether splitting the root into a
> mobile (shoot-feeding) pool and a sequestered (burden-holding) pool resolves the central
> mass-balance tension. Script: `validation/twopool_root_exploration.py` (standalone,
> EXPLORATORY, in-sample Yamazaki 2023; the canonical core and `parameters.json` are
> UNCHANGED). OOS transfer: `validation/twopool_root_oos.py`; long-chain shoot-floor
> diagnostic: `validation/twopool_root_seqrelease.py`; measured-forcing robustness re-fit:
> `validation/twopool_root_measured.py`. Figure: `validation/figures/twopool_root_exploration.png`.

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

## Result 4 — out-of-sample transfer (`validation/twopool_root_oos.py`)

The Yamazaki-fit two-pool U-shaped-`k_seq` model is transferred **without re-fitting** to two
independent datasets. All four models (two-pool, single-pool monotone `f_xy`, single-pool W2,
single-pool lipid) are run on the **same demo forcings** so the comparison is apples-to-apples.
(The fitted params are cached to `validation/twopool_fitted_params.json`.)

**(1) Kim 2019 brown-rice (grain) BAF, porewater basis — the decisive OOS series** (spans
PFHpA→PFDoDA and shows the long-chain grain rise):

| log10 RMSE vs Kim grain (same forcings) | two-pool | monotone | W2 | lipid |
|---|--:|--:|--:|--:|
| all | 0.62 | 1.40 | 0.66 | 1.02 |
| **excl PFOA** | **0.47** | 1.49 | 0.57 | 1.12 |
| reliable (DF ≥ 15 %: PFOA, PFNA) | 0.81 | 1.01 | 0.85 | 0.58 |

- **The two-pool transfers BEST on the clean series (excl PFOA 0.47)** and **captures the
  long-chain grain RISE** that the monotone model collapses on: 2pool PFUnDA 6.1 / PFDoDA 8.7
  vs monotone 0.19 / 0.52 (Kim obs ~33 / 35). This is the project's first OOS signal for the
  *two-pool* model specifically, and on identical forcings it beats the lipid model (which
  over-shoots the long chains here, 118 / 152).
- **Honest limits:** the absolute long-chain grain is still **under** (6–9 vs 33–35; these are
  the lowest-detection-frequency Kim congeners, DF 3–13 %); and Kim PFOA grain (4.43) is ~10×
  Yamazaki's (0.46) — a genuine between-dataset shift the Yamazaki-fit model cannot bridge
  (it under-predicts PFOA grain at 0.36). The headline supports the *mechanism/direction*, not
  the absolute long-chain level.

**(2) Li 2025 grain/root TF (water-independent ratio) — inconclusive, as documented.** Two-pool
grain/root TF RMSE 1.26 vs monotone 0.86: Li's short-chain grain/root TFs are anomalously high
(PFBS 19.3, PFHxA 7.9 — grain ≫ root for short chains), a root-surface / husk-confounding the
clean-root model does not reproduce. Consistent with the companion doc flagging Li as inconclusive.

**OOS conclusion:** the transfer **supports the structure and the long-chain mechanism** (best
clean-dataset transfer; reproduces the independent long-chain rise), but does **not** warrant
promoting the fitted `k_seq` into `parameters.json` — it is a single clean OOS set on demo
forcings, the absolute long-chain grain stays under, and Li is confounded. Keep exploration-only.

## Result 5 — the long-chain shoot floor is a SHOOT-loading ceiling, not a root problem (`twopool_root_seqrelease.py`)

The one residual after Result 3 is the very-long-chain **shoot** (PFDoDA straw 10.5 vs obs 49.8,
grain 8.7 vs 45.5). Since 97 % of the PFDoDA root burden sits in the *irreversible* seq pool, the
natural hypothesis was that the seq pool is a **slow buffer** — a desorption rate `k_rel` would
trickle the huge long-chain burden back to the mobile pool and feed the shoot. Tested by sweeping
`k_rel` (mass-conserving seq→mobile release added to the ODE; `k_rel = 0` recovers Result 3):

| `k_rel` (1/day) | all-11 RMSE | root | straw | grain | PFDoDA root/straw/grain |
|---|--:|--:|--:|--:|--:|
| 0.000 | 0.251 | 0.156 | 0.260 | 0.311 | 82.3 / 10.5 / 8.7 |
| 0.010 | 0.258 | 0.199 | 0.254 | 0.309 | 65.5 / 11.1 / 8.9 |
| 0.050 | 0.336 | 0.434 | 0.242 | 0.303 | 33.2 / 12.4 / 9.5 |
| 0.200 | 0.472 | 0.724 | 0.235 | 0.297 | 11.7 / 13.4 / 10.3 |

**`k_rel` cannot lift the shoot.** As it rises the PFDoDA straw barely moves (10.5 → 13.4) while
the **root collapses** (82 → 12) — the seq pool just equilibrates with the mobile pool and stops
retaining. The shoot is not starved for *mobile-pool burden*; it is limited by the **xylem-loading
capacity** itself. The `g_xy` diagnostic confirms it:

| `g_xy` × | PFDA straw | PFUnDA straw | PFDoDA straw | all-11 RMSE |
|---|--:|--:|--:|--:|
| obs | 3.5 | 8.2 | 49.7 | — |
| ×1 | 3.4 | 6.8 | 10.5 | 0.251 |
| ×4 | 10.1 | 18.6 | 26.5 | 0.469 |
| ×8 | 15.1 | 25.6 | 34.7 | 0.665 |

Reaching PFDoDA straw ≈ 50 needs `g_xy` ≳ ×8 (and even ×8 only gets to 35), which **over-feeds
PFDA/PFUnDA 3–4×** and balloons the RMSE to 0.665. **No smooth, K_PL-gated (QSPR-able) loading term
can selectively lift PFDoDA without wrecking PFDA/PFUnDA** — the long-chain straw needs a
per-congener boost for C12 alone. This independently quantifies PR #21's LC5/LC6 conclusion (an
active-carrier / xylem-loading **capacity limit**, not QSPR-able). Note too that the observed
PFDoDA straw is a **6× jump over PFUnDA for one CF₂** (vs only 3.5× in the root) — physically odd
and consistent with PFDoDA being a near-MQL outlier.

**Conclusion:** the long-chain shoot floor is a **shoot-side xylem-loading ceiling + a near-MQL
outlier**, structurally outside what any *root* term (`k_seq`, `k_rel`) can reach. The two-pool
root model at RMSE 0.251 is therefore at the achievable floor for this structure; the residual is
not a missing root mechanism. `k_rel`/`g_xy` are diagnostics only — the default model keeps
`k_rel = 0` and the Result-3 `g_xy`.

## Result 6 — robustness to MEASURED forcings (`twopool_root_measured.py`)

Everything above used the demo logistic forcings (identical to `reproduce_demo`), whose
transpiration peaks ~5× too high (CLAUDE.md). The decisive robustness check re-fits the whole
model on the **measured forcings** — `forcing_rice.Q_TP` (peak 0.098 L/d/hill, T/ET 0.42) +
`growth_rice` ORYZA-IR72 organ biomass (HI 0.53) — the same forcings the fxy-doc baselines use, so
the Kim-grain OOS becomes directly comparable.

| | in-sample RMSE | root | straw | grain | PFOS/PFUnDA k_seq sep. |
|---|--:|--:|--:|--:|--:|
| two-pool, demo forcings (Result 3) | 0.251 | 0.156 | 0.260 | 0.311 | 3.1× |
| **two-pool, MEASURED forcings** | **0.278** | **0.154** | 0.295 | 0.348 | **4.5×** |
| fxy-doc U-shaped-K_PL-f_xy, measured (single pool) | 0.286 | — | — | — | — |

- **The qualitative result is robust.** On realistic biomass/transpiration the two-pool still
  **solves the root (0.154)**, ties the fxy-doc U-shaped-K_PL-f_xy in-sample (0.278 vs 0.286)
  **while keeping the monotone physical `f_xy_recommended`**, and the **non-K_PL PFOS/PFUnDA
  separation HOLDS and even sharpens to 4.5×** (k_seq 0.031 vs 0.141). The fitted U-shape is the
  same family (rising arm slightly steeper); `g_xy` roughly doubles (0.021 → 0.040) to offset the
  ~4× lower transpiration.
- **Kim 2019 grain OOS, now apples-to-apples** (measured forcings, vs the fxy-doc baselines):

| Kim grain log10 RMSE (measured forcings) | two-pool | lipid | monotone | W2 |
|---|--:|--:|--:|--:|
| excl PFOA | **0.56** | 0.55 | 2.04 | 1.11 |

  The two-pool **ties the best prior model (lipid, 0.55)** on the clean Kim series and crushes
  monotone/W2 — but does so while *also keeping the high long-chain root* that the lipid model
  drains. The long-chain shoot stays under (PFDoDA straw 7.5; sharper than the demo's 10.5 because
  the measured transpiration is lower) — the same Result-5 xylem-loading ceiling.

**Robustness conclusion:** the two-pool structure, the monotone physical `f_xy`, the non-K_PL
U-shaped `k_seq`, the PFOS/PFUnDA separation, and the OOS transfer all survive the switch from demo
to measured forcings. The result is not an artifact of the placeholder transpiration. (Fit cached
to `validation/twopool_fitted_params_measured.json`; still in-sample, `parameters.json` untouched.)

## Honest status

- **What works:** the two-pool *structure* decouples root burden from shoot delivery, lets the
  model keep the monotone physical `f_xy_recommended`, and — with the **U-shaped `k_seq(n)`**
  (Result 3) — reaches log10 RMSE **0.251 across all 11 congeners including PFDoDA** (root 0.156),
  reproducing the non-K_PL **PFOS/PFUnDA separation** (3.1×) that B/K_PL cannot. The root-matched
  test proves the structure is **sufficient**; the U-form realizes it with a parsimonious smooth
  descriptor. Mass-conserving; the canonical core and `parameters.json` are untouched.
- **The remaining residual is the very-long-chain SHOOT, not the root — and it is structural
  (Result 5).** PFDoDA straw/grain stay ~3–5× under. A slow seq-pool release (`k_rel`) cannot lift
  it (the root collapses before the shoot moves), and the `g_xy` diagnostic shows the bottleneck is
  the **xylem-loading capacity**: reaching PFDoDA straw needs `g_xy` ×8, which over-feeds
  PFDA/PFUnDA 3–4× — no smooth/QSPR-able term selectively fixes C12. This is a *shoot*-loading
  ceiling (corroborating PR #21 LC5/LC6) plus a near-MQL outlier, outside what any root term reaches.
- **OOS-checked (Result 4):** the U-shaped `k_seq` transfers best on the clean Kim 2019 grain
  series (excl-PFOA RMSE 0.47) and reproduces the independent long-chain rise, but stays
  exploration-only (single clean OOS set, demo forcings; `parameters.json` untouched).
- **Mechanistic reading:** the model now says two *distinct* processes set the root BAF — a
  reversible membrane/lipid partition (`B_m`, ~K_PL, feeds the shoot) and an irreversible
  chain·head-group-specific sink (`k_seq`, non-K_PL, holds the burden). Long-chain PFCAs load the
  sink heavily (PFUnDA/PFDoDA), PFSAs much less at matched chain (PFOS), and short chains barely
  at all in absolute terms.
- **Decisive experiment (unchanged):** per-congener xylem-sap / root-water ratio (direct f_xy)
  **and** a desorption-resistant root-fraction assay (isolating the irreversible `k_seq` pool)
  across chain length and head group. Still in-sample (Yamazaki only); this is mechanism
  discovery and a well-posed next fit, **not** validation.

## API access (opt-in, handoff item ①)

The two-pool model is now callable through the UI-agnostic API as
`model_api.simulate_twopool_seq(congener, …)`, mirroring `simulate_nstem_leaf` — so the app and
other validation can use it **without changing any default** (the canonical `simulate`,
`reproduce_demo`, and `parameters.json` are untouched). It loads the cached fit from
`validation/twopool_fitted_params.json` (via this module's `load_fit()`/`kseq_ushape()`/`lipid_g()`)
and re-implements the 5-state ODE inside `model_api`, driven by the standard
forcing/`drivers=` machinery, so it returns the **same dict shape** as `simulate()`
(root/stem/leaf/grain conc & BAF series, finals, `straw_baf`, `tf_final`) plus the root
**mobile/seq split** (`conc["root_mobile"|"root_seq"]`, `seq_fraction`) and the two-pool
levers (`k_rel`, `kseq_override`). The reported `root` BAF is the **sum** of the mobile and
sequestered pools.

```python
import model_api as api
r = api.simulate_twopool_seq("PFUnDA")           # defaults reproduce the demo-forcing headline
r["baf_final"]["root"], r["seq_fraction"]    # 15.8, 0.91  (root = mobile + seq)
r["params"]["k_seq"], r["params"]["f_xy"]    # 0.166, monotone physical f_xy_recommended
api.simulate_twopool_seq("PFUnDA", k_rel=0.2)    # Result-5 desorption sweep (root collapses)
```

With the defaults (`measured_forcing=False, season=120`) the wrapper reproduces this record's
headline — overall log10 RMSE **0.251** (root 0.156) and the non-K_PL PFOS/PFUnDA **3.1×**
`k_seq` separation — to within the ~1 % driver-grid difference (cross-impl RMSE 0.014). A drift
guard (`tests/test_model_api.py::test_simulate_twopool_seq_matches_validation_and_rmse`) pins the two
implementations together. Still EXPLORATORY / in-sample: the cached fit is on the demo forcings,
and `twopool_fitted_params_measured.json` (Result 6) is **not** auto-loaded.

> **Not to be confused with** the *carrier* two-pool model (`model_api.simulate_twopool_carrier` /
> `close_longchain_2pool`, `src/pfas_rice_two_pool.py`): that is a separate opt-in two-pool root model
> whose second pool is a reversible **bound store** tuned by carrier/f_xy levers (the saturated
> long-chain closure), whereas THIS model's second pool is an irreversible **sequestration sink**
> (`k_seq`). The `_seq`/`_carrier` suffixes disambiguate the two `model_api` entry points.

---

## Result 7 — Tang 2026 per-organ OOS (`validation/twopool_root_oos_tang.py`): a NEGATIVE/diagnostic result

The handoff (`docs/HANDOFF_BAF_twopool.md` §4.2) flagged Tang 2026 per-organ TF
(stalk/leaf/endosperm, dw) as a natural next OOS for the two-pool. Transferring the
Yamazaki-fit two-pool to Tang **without re-fitting** gives:

| model (OOS, Tang 0.1 µg/g) | overall log10 RMSE | stalk | leaf | endosperm |
|---|---|---|---|---|
| two-pool U-shaped k_seq | **1.40** | 1.89 | **0.38** | 1.46 |
| single-pool monotone f_xy | 1.23 | 1.15 | 0.98 | 1.50 |
| single-pool lipid loading | **0.52** | 0.47 | 0.58 | 0.49 |
| single-pool Tang-refit f_xy (in-sample) | 0.52 | — | — | — |

The two-pool is **WORSE** than even the single-pool monotone, and far worse than the
lipid model (which remains the Tang winner, consistent with `oos_tang_lipid.py`).

**Why — and why it is informative, not a failure of the root mechanism.** Tang
per-organ TF is a **SHOOT-resolution** test (stalk vs leaf vs endosperm), but the
two-pool's entire innovation is in the **ROOT** (mobile/seq split). Its shoot is the
**unmodified basic 4pool** with a **pass-through stem** (PFOA stem conc 0.008 vs leaf
1.14), so the **stalk TF collapses** — the documented over-translocation / empty-stem
defect that `nstem_leaf` (redistributed shoot + transpiration retention) was built to
fix. The per-organ breakdown isolates this cleanly: the two-pool's **leaf RMSE (0.38)
is the BEST of all three models**, and only the stalk drags the overall up. The
single-pool baselines here use `nstem_leaf`, so their stalk is populated — the stalk
comparison is **apples-to-oranges** (different SHOOT model), not a root-mechanism
difference.

Compounding it: Tang's congeners are **C5–C8** (GenX, PFOA, PFOS), so the **long-chain
root decoupling — the two-pool's whole point — is not even exercised**. (GenX further
over-predicts via the provisional ether `f_xy_recommended`, a separate condition-
dependent QSPR issue.)

**Conclusion.** The two-pool ROOT mechanism and the Tang per-organ SHOOT pattern are
largely **orthogonal**; **Tang is not a suitable OOS test of the two-pool root**. The
actionable, structural finding is that a fair per-organ Tang comparison needs the
two-pool **root** merged with the `nstem_leaf` **redistributed shoot** (a future
structural merge). **Kim 2019 grain** (`twopool_root_oos.py`, Result 4 — grain/root,
not shoot-resolved) remains the informative two-pool OOS. This does **not** add support
for promoting the fitted `k_seq` into `parameters.json`.

Guarded by `tests/test_model_api.py::test_twopool_simulate_organs_and_tang_passthrough_diagnosis`
(pins the `simulate`/`simulate_organs` solve-path consistency + the pass-through-stem
diagnosis). `simulate_organs(c, p, kseq_override=, k_rel=)` was added to
`twopool_root_exploration.py` to expose the per-organ (stem/leaf) split on the same
solve path as `simulate`.
