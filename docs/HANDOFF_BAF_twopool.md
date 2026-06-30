# HANDOFF — BAF 예측 고찰: two-pool root + U-shaped k_seq

> Session handoff for the next Claude/dev. Picks up the BAF-prediction "고찰" arc.
> Branch: **`claude/peaceful-babbage-rei1a1`**  ·  PR: **#24** (draft, base `main`).
> Full scientific record: **`docs/twopool_root_exploration.md`** (read this first).
> Everything below is EXPLORATORY / in-sample; the canonical core and
> `params/parameters.json` are **UNCHANGED** (opt-in only).

---

## 0. TL;DR (where we are)

The central BAF wall — *a single root pool cannot hold a high long-chain root BAF
**and** deliver a non-trivial long-chain shoot BAF* — is **broken** by a two-pool
root (mobile shoot-feeding pool + sequestered burden-holding pool). The root sink
descriptor `k_seq` is **non-K_PL and U-shaped in chain length**; it realizes the
PFOS/PFUnDA separation (identical K_PL=31623, root 5.93 vs 19.53) that B/K_PL cannot.
The model reaches log10 RMSE **0.251** (all 11, incl PFDoDA; root 0.156) while keeping
the **monotone physical `f_xy_recommended`**. It transfers out-of-sample to Kim 2019
grain as well as the best prior model, and the whole story is **robust to measured
forcings**. The one remaining residual (very-long-chain shoot) is a diagnosed
**xylem-loading ceiling**, not a root problem.

**Status: a coherent, closed mechanism-discovery arc. NOT validation** (Yamazaki
in-sample fit → OOS transfer; the decisive wet-lab experiment is still missing).

---

## 1. What this session delivered (PR #24, 5 commits on top of `8da1e39`)

| commit | what |
|---|---|
| `10ea7b9` | two-pool 5-state ODE; global fit + **root-matched sufficiency test** |
| `94e18dc` | **U-shaped `k_seq(n)`** realized → PFOS/PFUnDA separation, RMSE 0.251 |
| `b44a8a3` | **OOS transfer** (Kim 2019 grain, Li 2025 TF) |
| `d307726` | long-chain **shoot-floor diagnosis** (k_rel sweep + g_xy diagnostic) |
| `c143adb` | **robustness re-fit on measured forcings** (capstone) |

**New files** (all under `validation/` + one doc):
- `twopool_root_exploration.py` — the core standalone 5-state ODE + fit + U-shape +
  figure. Exposes `simulate(c, p, kseq_override=, k_rel=)`, `kseq_ushape(n, grp, q)`,
  `compute_fit_quiet()`, `load_fit()` (cached), `CONGENERS`, `OBS`.
- `twopool_root_oos.py` — Kim/Li OOS transfer (reuses cached fit).
- `twopool_root_seqrelease.py` — k_rel seq-release sweep + g_xy bottleneck diagnostic.
- `twopool_root_measured.py` — re-fit on `forcing_rice` + `growth_rice` forcings.
- `twopool_fitted_params.json` / `twopool_fitted_params_measured.json` — cached fits.
- `figures/twopool_root_exploration.png` — (a) U-shaped k_seq + fit; (b) pred-vs-obs.
- `docs/twopool_root_exploration.md` — Results 1–6, honest status. **Source of truth.**

CLAUDE.md §6 (two-pool bullet, items 1–7) and §7 (run commands) are updated.

---

## 2. The model (1 paragraph)

5-state ODE `[root_mobile, root_seq, stem, leaf, grain]`, mass-conserving, sole source
`M_root·j_R` into the mobile pool.
- **root_mobile**: binding `B_m` (= measured basis-A B_root); GHK+carrier uptake; loads
  xylem with `f_xy_recommended·Cw_m + g_xy·C_m` (monotone physical f_xy + K_PL-gated
  lipid term).
- **root_seq**: receives `k_seq·C_m` (released only by growth dilution → terminal
  accumulator; optional slow desorption `k_rel`, default 0). Holds the root burden
  without draining the shoot feed.
- `k_seq(n, group)` = **non-K_PL** chain·head-group descriptor. The fitted U-form
  (demo forcings): `k_seq = [0.268·e^(−0.52(n−4)) + 0.615·e^(1.35(n−12))]·{10^+0.18 if PFSA}`.

Root BAF = `C_mobile + C_seq`.

---

## 3. Key results (numbers to trust)

| result | metric |
|---|---|
| U-shaped k_seq, demo forcings | RMSE **0.251** (all 11); root 0.156, straw 0.260, grain 0.311 |
| PFOS/PFUnDA separation | k_seq 0.054 vs 0.166 = **3.1×** (demo), **4.5×** (measured) |
| root-matched empirical k_seq | PFOS 0.047 vs PFUnDA 0.210; **U-shaped in n** (why linear `ks_b→0`) |
| OOS Kim grain (demo forcings) | two-pool excl-PFOA **0.47** = best (mono 1.49, W2 0.57, lipid 1.12) |
| OOS Kim grain (measured forcings) | two-pool excl-PFOA **0.56** = **ties lipid 0.55**; mono 2.04, W2 1.11 |
| measured-forcing in-sample | RMSE **0.278**, root 0.154 (ties fxy-doc U-K_PL-f_xy 0.286) |
| long-chain shoot floor | k_rel can't lift; needs g_xy ×8 → over-feeds PFDA/PFUnDA, RMSE→0.665 |

**Honest caveats** (carry these forward, do not drop):
- In-sample fit is **Yamazaki only**; OOS is a single clean dataset (Kim grain) +
  one confounded one (Li, inconclusive). Demo *and* measured forcings both tested.
- PFDoDA is a **near-MQL outlier** (obs straw is a 6× jump over PFUnDA for one CF₂ vs
  3.5× in root) and caps the achievable RMSE.
- The two-pool ties (not beats) the lipid model on OOS; its advantage is keeping the
  **high long-chain root** that lipid drains, with the **monotone physical f_xy**.

---

## 4. Open questions / candidate next steps (prioritized)

1. **Promotion decision (the live one).** Should the two-pool + U-shaped `k_seq` move
   from exploration-only into the canonical path? Options:
   - (a) Wire it as an **opt-in module** (`model_api.simulate_twopool_seq(...)`, like
     `simulate_nstem_leaf`) — low-risk, makes it usable by the app/other validation
     without changing defaults. **Recommended first concrete step. ✅ DONE** (this
     session): `model_api.simulate_twopool_seq(...)` returns the standard `simulate()` dict
     + the root mobile/seq split; loads the cached fit; reproduces the headline
     (RMSE 0.251, PFOS/PFUnDA 3.1×) with the monotone physical f_xy; drift-guarded by
     `tests/test_model_api.py::test_simulate_twopool_seq_matches_validation_and_rmse`.
     Defaults / `reproduce_demo` / `parameters.json` unchanged. See
     `docs/twopool_root_exploration.md` §"API access".
   - (b) Promote the U-shaped `k_seq` coefficients into `parameters.json` (new fields,
     provenance-tagged). **✅ DECIDED (with the user) — DO NOT promote; keep opt-in.**
     The decision was made after assembling all the evidence item ② produced:
     - the **Tang per-organ OOS adds no support** (#2 below — NEGATIVE/diagnostic);
     - the `k_seq` **mechanism review** (`docs/twopool_kseq_mechanism.md`, #5 below)
       strengthens the *story* (PARTIALLY SUPPORTED) but provides **no direct measured
       root-sorbent dataset** — `k_seq` stays anchored by analogy + superposition;
     - the standing honest caveats remain: **Yamazaki in-sample fit**, a **single clean
       OOS** (Kim grain), **PFDoDA near-MQL**.
     → The two-pool stays an **opt-in** module (`model_api.simulate_twopool_seq`, option
     (a)); defaults / `reproduce_demo` (RMSE 0.029) / `parameters.json` **unchanged**.
     **Promotion gate**: the `twopool_kseq_mechanism.md` §5 rice-root **cell-wall /
     Fe-Mn-plaque batch-sorption + desorption assay** (chain-length × head-group). Until
     that direct measurement exists, promotion is not warranted. (Re-open only with that
     data, or a second independent clean per-organ OOS once the two-pool root is merged
     with the `nstem_leaf` redistributed shoot — see #2.)
2. **More OOS datasets.** Kim + Li done (Result 4). **Tang 2026 per-organ TF** now also
   transferred (`validation/twopool_root_oos_tang.py`, **Result 7**) — but it is a
   **NEGATIVE/diagnostic** result: the two-pool OOS RMSE **1.40** is WORSE than the
   single-pool monotone (1.23) / lipid (0.52), because Tang per-organ is a **SHOOT**
   test and the two-pool's innovation is in the **ROOT** (its shoot is the basic 4pool
   with a **pass-through stem** → the stalk TF collapses; the two-pool's **leaf** RMSE
   0.38 is actually the best of all models). Tang's congeners are C5–C8, so the
   long-chain root decoupling is not even exercised. **Conclusion: Tang is not a fair
   OOS of the two-pool root** — a per-organ Tang test needs the two-pool root merged
   with the `nstem_leaf` redistributed shoot (a future structural merge). Kim grain
   stays the informative two-pool OOS. ✅ DONE (this session).
3. **The decisive experiment (cannot be done in-silico).** Per-congener
   **xylem-sap / root-water ratio** (direct f_xy, root-pressure exudate) **+** a
   **desorption-resistant root-fraction assay** (isolates the irreversible `k_seq` pool)
   across chain length and head group. This is the data gap that would turn discovery
   into validation. Flag it to the user / write it into the data-gap doc.
4. **Long-chain shoot floor (declared structural).** Could try a **carrier-saturation /
   chain-superlinear `g_xy`** term to model the C12 shoot directly — but Result 5 argues
   it is NOT QSPR-able (over-fits one near-MQL outlier). **Low priority / likely a
   negative result**; only pursue if the user wants the shoot side modeled explicitly.
5. **U-shaped `k_seq` provenance.** ✅ **DONE (literature search)** —
   `docs/twopool_kseq_mechanism.md` (deep-research, 17 sources → 25 verified claims).
   **Verdict: PARTIALLY SUPPORTED.** The U-shape is a **superposition of two distinct
   mechanisms** (short-chain electrostatic/anion-exchange arm + long-chain hydrophobic/
   **desorption-resistant** arm; TII→0.98 at C10), which explains the model's two-
   exponential-arm form; the **irreversible sink** and the **PFSA offset** (`10^+0.18`
   fitted ≈ measured **+0.23 log** sulfonate, lignin/soil) are now anchored. **Fe-Mn
   plaque demoted** (outer-sphere/pH-reversible/acidic-only) in favour of cell-wall
   entrapment. **Central data gap**: no rice-root cell-wall / Fe-Mn-plaque PFAS
   coefficient resolved by BOTH chain length AND head group exists — anchored by analogy,
   not a direct dataset. So this strengthens the mechanistic *story* but does NOT by
   itself warrant promoting `k_seq` (the §5 cell-wall/plaque batch-sorption experiment
   is the gate). Remaining open item is now only #1 (promotion decision) + #3 (the
   decisive wet-lab experiment, now scoped in `twopool_kseq_mechanism.md` §5).

---

## 5. How to resume (commands)

```bash
# env
pip install -r requirements.txt           # numpy, scipy (matplotlib for the figure)

# reproduce the whole arc
python validation/twopool_root_exploration.py     # fit + root-match + U-shape + figure (~3 min)
python validation/twopool_root_oos.py             # Kim/Li OOS transfer (reuses cache, ~5 s)
python validation/twopool_root_oos_tang.py        # Tang per-organ OOS (NEGATIVE/diagnostic, ~25 s)
python validation/twopool_root_seqrelease.py      # shoot-floor diagnostic (~20 s)
python validation/twopool_root_measured.py        # measured-forcing robustness re-fit (~3 min)

# the cached fits (delete to force a re-fit):
#   validation/twopool_fitted_params.json           (demo forcings)
#   validation/twopool_fitted_params_measured.json  (measured forcings)
```

To build on the model in code: `import twopool_root_exploration as TP` (from
`validation/`), then `TP.load_fit()` → `(p, q)`, and
`TP.simulate(c, p, kseq_override=TP.kseq_ushape(c["n_C"], c["group"], q), k_rel=0.0)`
for any `c in TP.CONGENERS`.

---

## 6. Status of the arc — all in-silico items CLOSED

The two-pool BAF arc's in-silico work is **complete**. Items ①–③ are resolved:

- **① opt-in module** (`model_api.simulate_twopool_seq`) — ✅ DONE (§4.1a).
- **② OOS extension** — ✅ DONE. Kim/Li (Result 4) + **Tang per-organ**
  (`validation/twopool_root_oos_tang.py`, Result 7; NEGATIVE/diagnostic) + the
  **`k_seq` mechanism literature review** (`docs/twopool_kseq_mechanism.md`;
  PARTIALLY SUPPORTED).
- **③ promotion decision** — ✅ DECIDED with the user (§4.1b): **DO NOT promote;
  keep opt-in.** `parameters.json` / defaults / `reproduce_demo` (RMSE 0.029) UNCHANGED.

**What remains is NOT in-silico** — it is the **promotion gate / decisive experiment**:
the `twopool_kseq_mechanism.md` §5 rice-root **cell-wall / Fe-Mn-plaque batch-sorption +
desorption assay** (chain-length × head-group), and/or a second independent clean
per-organ OOS once the two-pool **root** is merged with the `nstem_leaf` **redistributed
shoot** (the structural merge that would make a per-organ Tang test fair; Result 7).

Carry the standing honest caveats forward on every claim: **Yamazaki in-sample fit,
single clean OOS, PFDoDA near-MQL**. Do not promote `k_seq` into `parameters.json`
without the §5 direct measurement.
