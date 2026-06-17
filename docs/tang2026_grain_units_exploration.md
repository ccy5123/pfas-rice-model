# Tang 2026 grain TF — a fresh/dry units correction, and the structural grain gap

> Exploration note (companion to `VALIDATION_TANG2026_KR.md`,
> `validation/tang2026_validation.py`, `validation/tang2026_fxy_refit.py`).
> Started as an attempt to close the grain/endosperm gap by re-fitting the phloem
> loading `L_Ph`; ended by finding (a) a fresh/dry **units inconsistency** in the
> per-organ TF comparison, and (b) that, once corrected, the grain gap is
> **structural** — `L_Ph` cannot close it.

## TL;DR
1. **Units.** Tang's TF is **dry-weight / dry-weight**; the model's tissue
   concentrations are **fresh-weight** (`C = B_k·Cw`, basis A). TF is independent of
   the *exposure* basis (`Cw` cancels: Tang soil vs model pore water) — but **not** of
   the *tissue-moisture* basis, because the fresh→dry factor `(1−θ_fw)` differs between
   tissues (root θ=0.90 vs grain θ=0.14) and does **not** cancel:
   `TF_dw = TF_fw · (1−θ_root)/(1−θ_tissue)`  →  ×0.59 stem, ×0.45 leaf, **×0.116 grain**.
2. **The old "grain matches" was a fw/dw artifact.** `tang2026_validation.py` previously
   compared the **fresh-weight** model TF to Tang's **dry-weight** TF. That made PFOS
   endosperm read "0.80 ≈ 0.77 (Tang)". Corrected to dw, it is **0.09 vs 0.77 — ~8× under**.
3. **The grain gap is structural.** With the f_xy re-calibrated to Tang (shoot) and the
   shoot distribution restored (nstem_leaf), **no `L_Ph` (even maxed) nor the lipid `g_ph`
   term** lifts grain TF_dw to Tang's endosperm. It is not an `L_Ph` calibration issue.

## 1. The units inconsistency (now fixed)
`tang2026_validation.py.model_tf` returned `baf_final[tissue]/baf_final[root]` — a
fresh-weight ratio — and compared it to Tang's dry-weight S8 TF. The header justified
this by "TF is denominator-free / independent of the exposure basis", which is true for
the **Cw** axis but was wrongly extended to the **fw↔dw** axis. The two other Tang
scripts (`tang2026_fxy_TF_validation.py`, `tang2026_fxy_refit.py`) already dw-converted,
so the repo was internally inconsistent. `model_tf` now applies the dw factor; all three
scripts are on the same (correct, dry-weight) basis.

### Corrected per-organ TF (dw), model vs Tang (across-dose mean)
| congener | tissue | Tang | monotone | W2 | lipid |
|---|---|---|---|---|---|
| PFOA | stalk | 1.45 | 0.02 | 0.01 | 0.05 |
| PFOA | leaf | 1.66 | 2.70 | 1.74 | 5.92 |
| PFOA | **endosperm** | **0.95** | 0.05 | 0.03 | **0.27** |
| PFOS | stalk | 0.58 | 0.00 | 0.05 | 0.23 |
| PFOS | leaf | 0.68 | 0.07 | 0.72 | 2.57 |
| PFOS | **endosperm** | **0.77** | 0.00 | 0.01 | **0.09** |
| GenX | stalk | 1.10 | 0.15 | 0.15 | 0.10 |
| GenX | leaf | 1.38 | 32.2 | 32.2 | 21.5 |
| GenX | **endosperm** | **1.39** | 1.69 | 1.69 | **1.19** |

(log10 RMSE, lipid: all 0.85, leaf+grain 0.73 — the dw fix slightly *helps* overall,
because it pulls the over-predicted leaf down more than it pushes grain/stem down; but
the per-organ story is now honest: **stem under, leaf over, grain badly under** for
PFOA/PFOS. GenX grain is ~ok only because GenX over-translocates everything.)

## 2. The grain ceiling — why `L_Ph` cannot close it
Grain TF_dw at the **Tang-refit f_xy** (PFOA 0.097 / PFOS 0.320 / GenX 0.017), sweeping
the phloem loading and the lipid term (nstem_leaf + ORYZA biomass):

| congener | `L_Ph`=0 | `L_Ph`=0.5 | +lipid `g_ph` | **Tang endosperm (0.1 µg/g)** |
|---|---|---|---|---|
| PFOA | 0.32 | 0.51 | 0.69 | **1.37** |
| PFOS | 0.09 | 0.11 | 0.29 | **0.94** |
| GenX | 0.17 | 0.29 | 0.25 | **1.62** |

And across model structures (refit f_xy, dw):

| congener | 4-pool | 4-pool+lipid | nstem_leaf | nstem+lipid | Tang |
|---|---|---|---|---|---|
| PFOA | 0.18 | 0.61 | 0.33 | 0.69 | 1.37 |
| PFOS | 0.05 | 0.21 | 0.09 | 0.29 | 0.94 |
| GenX | 0.21 | 0.24 | 0.22 | 0.25 | 1.62 |

The best any configuration reaches is ~0.7 (PFOA, nstem+lipid) — still ~2× under; PFOS/GenX
are ~3–6× under. **`L_Ph` saturates well below Tang.**

## 3. Interpretation
- Per **dry** gram, Tang's endosperm holds ≈1.0–1.6× the root. In **fresh** terms that
  is ≈9–14× the root (grain is ~dry, root is ~90% water), which phloem delivery + the
  *measured* grain binding `B_grain` (Chen K_PL / Zhou K_prot) cannot produce.
- **Restoring the shoot distribution starves the grain.** The 4-pool over-delivers to
  grain through its leaf-sink (an artifact); nstem_leaf (correct stem~leaf split) lowers
  grain. You cannot get the shoot split **and** the grain right with the current phloem
  structure — they trade off.
- So grain accumulation needs a **grain-specific mechanism** (e.g. endosperm
  storage-protein/starch sequestration beyond the current `B_grain`, or a phloem-unloading
  enrichment step), **not** an `L_Ph` re-fit. This mirrors the documented root-sequestration
  negative result (`data_inventory_and_gaps.md` §3).

## 4. What changed in the repo
- `validation/tang2026_validation.py` — `model_tf` now dw-converts; figure regenerated.
- `validation/tang2026_nstem_validation.py` — `baseline_tf`/`nstem_tf` now dw-convert.
  This MOVED the headline numbers: shape RMSE 0.84→0.11 ⇒ **0.85→0.39**, overall
  1.28→1.01→0.18 ⇒ **1.53→1.20→0.71** (the nstem_leaf pattern cure is real but more
  modest; the grain under-prediction now shows). `VALIDATION_TANG2026_NSTEM_KR.md` + the
  CLAUDE.md §6 nstem bullet corrected to match.
- `VALIDATION_TANG2026_KR.md` — grain "match" corrected (units note added).
- `validation/tang2026_fxy_refit.py` — note added that `L_Ph` cannot close the grain gap.
- No model or `parameters.json` change; the Tang-refit f_xy stays override-only.

### Exhaustive fw/dw audit (this session)
Every script that compares model output to observed data was checked for the basis on
both sides. Only the **two Tang TF scripts above** had the mismatch (they consume the
**dry-weight** `raw_si/tang2026_doseresponse.csv` and compared a fresh-weight model TF).
All other comparisons are self-consistent:
- `data_obs/obs_baf_{Yamazaki,Li2025}.csv` are **fw-converted** (so `reproduce_demo`,
  `calibration`, `revalidation_crosscheck`, `validation_summary`, `oos_crossdataset`
  compare fw↔fw);
- `validation/S6_Gap4.py` uses the `_fw` (not `_dw`) columns of `Li2025_BAF_TF.csv`;
- `validation/tang2026_fxy_{TF_validation,refit}.py` already dw-converted.

## 5. Open follow-up
A grain-enrichment mechanism is a **data** question first: per-organ subcellular PFAS in
the **filling grain** (endosperm vs aleurone vs pericarp), and phloem-sap [PFAS] at the
peduncle, to decide between (i) higher effective `B_grain` and (ii) active phloem
unloading. Until then, grain TF should be reported as an **order-of-magnitude** quantity
with this known low bias (~3–8× for PFCA/PFSA on a dry-weight basis).
