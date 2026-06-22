# Rice biomass partitioning & the root:shoot anchor

> Status: literature-grounded for the **shoot split**; the **root fraction** is now
> anchored to a site-relevant field value (Japanese flooded paddy, root:shoot ~0.08–0.12,
> §2) but remains condition-dependent (water table, cultivar) so it is given as a range.
> Scope: explains the organ-biomass `M_k(t)` that feeds the burden = `C_k·M_k`
> question and the calibration coupling it creates. No default behaviour changed.

## 1. Why this matters

"Where does PFAS accumulate?" is a **burden** (`C_k·M_k`, µg) question, and burden
depends directly on the per-organ biomass `M_k`. A wrong partitioning silently
biases the burden ranking — and, as shown in §4, the model's transport calibration
turns out to be **entangled** with the biomass assumption.

## 2. What the literature says (maturity, lowland rice)

**Aboveground split — well established, and the model matches it:**

| organ | % of aboveground DM at maturity | source |
|---|---|---|
| panicle / grain | ~45–58% (HI 0.45–0.55, modern) | Ntanos & Koutroubas 2002; Amanullah & Inamullah 2016 |
| stem / culm + sheath | ~24–33% | Amanullah & Inamullah 2016 |
| leaf blade | ~18–24% | Amanullah & Inamullah 2016 |

**Root fraction — the weak link (but now anchored to a site-relevant field value):**

- Root mass fraction declines from ~0.2 (seedling) to ~0.1 (heading) of total plant;
  at **maturity it is ~0.07–0.13** (root:shoot ~0.08–0.15), pushed lower in high-yield
  cultivars by grain-fill dilution (root:shoot negatively correlates with yield).
- A concrete heading anchor: shoot 3 Mg/ha ↔ root ~0.33 Mg/ha → root:shoot ≈ 0.11.
- **Site-relevant anchor (Japanese flooded paddy — same system class as the Yamazaki
  calibration data):** measured root biomass ~70–112 g/m² against an aboveground
  biomass ~910 g/m² at harvest → **root:shoot ≈ 0.08–0.12** (root ~7–11% of total).
  This independently lands on the same ~0.10 and is the closest field analog to the
  model's target system. Sources: Wang et al. 2016 *Sci. Rep.* (10.1038/srep29333,
  paddy root biomass g/m²); Frontiers Plant Sci. 2021 (10.3389/fpls.2021.713814,
  japonica paddy — RDW peaks at filling then *declines* to maturity while SDW keeps
  rising, so root:shoot falls toward harvest).
- Caveat: BGB:AGB is **water-management dependent** (lower under saturated/flooded vs
  drained/peat soils — *Sci. Rep.* 2024, 10.1038/s41598-024-64616-1), and field root
  recovery under-counts deep/broken roots, so measured values are if anything a *floor*.
- **Combine, don't co-locate**: the full per-organ distribution is built from the
  shoot split (above) × an *independent* root:shoot ratio — the two need not come from
  the same paper. Combining gives (total-plant %): root 9–13, grain 46–48, stem 25–26,
  leaf 16–17 for root:shoot 0.10–0.15.

### Citations (DOI)
- Ntanos & Koutroubas 2002, *Field Crops Res.* 74:93–101 — **10.1016/S0378-4290(01)00203-9** (maturity HI 0.47–0.61; aboveground only).
- Amanullah & Inamullah 2016, *Rice Science* 23(2):78–87 — maturity panicle/culm/leaf % (aboveground only).
- Nada & Abogadallah 2018, *Acta Physiol. Plant.* 40:123 — **10.1007/s11738-018-2697-5** (root:shoot is genotype/condition dependent; pot/greenhouse).
- "Response of Grain Yield and Root … to Nitrogen Levels in Paddy Rice", *Front. Plant Sci.* 2021 — **10.3389/fpls.2021.713814** (japonica paddy, tillering→maturity SDW/RDW/root:shoot; RDW peaks at filling then declines to maturity).
- Wang et al. 2016, "Optimizing rice plant photosynthate allocation reduces N₂O emissions from paddy fields", *Sci. Rep.* 6:29333 — **10.1038/srep29333** (paddy root biomass ~70–112 g/m²; HI–N₂O link).
- *Sci. Rep.* 2024 — **10.1038/s41598-024-64616-1** (rice BGB:AGB is water-table/soil dependent).
- Yoshida 1981, *Fundamentals of Rice Crop Science*, IRRI (classic reference for root fraction ~5–10%).
- Model partitioning source: Bouman & van Laar 2006, *Agric. Syst.* 87:249–273 (**10.1016/j.agsy.2004.09.011**); Li et al. 2017, *Agric. For. Meteorol.* 237–238:246–256 (ORYZA v3).

> Honest caveat: none of these is the model's exact target system (e.g. Yamazaki
> Andosol greenhouse). Partitioning is genotype × environment dependent, so the
> literature gives a defensible **range**, not a site truth. A maturity, root-
> inclusive measurement for the target system would close the gap.

## 3. What the model assumes vs literature

| biomass driver | root % of total | root:shoot | HI |
|---|---|---|---|
| `reproduce_demo` placeholder logistic | 23 | **0.30** | **0.07** |
| `growth_rice` (`simulate` default) | 4.7 | **0.049** | 0.51 |
| `oryza_growth` (**app default**, mechanistic) | 6.2 | **0.066** | 0.44 |
| literature (maturity) | 7–13 | 0.08–0.15 | 0.45–0.55 |

The **shoot split matches the literature**; only the **root fraction is a bit low**.
Both `growth_rice` and `oryza_growth` use the **experimental ORYZA2000 IR72 root
partitioning table (FRTTB)** — `FRT = 0.50, 0.25, 0.00` at `DVS = 0, 0.43, 1.0`
(Bouman & van Laar 2006; IRRI field-calibrated). That experimental coefficient already
lifts root:shoot to **0.049** (`growth_rice`) / **0.066** (`oryza_growth`, the app
default), i.e. most of the way to the literature 0.08–0.13 with **no tuning**. The
residual is ORYZA's known under-prediction of post-flowering root (FRTTB → 0 at flowering,
no root turnover modelled) plus, for `growth_rice`, its grain-fill-weighted biomass
logistic. The old `growth_rice` value (0.035) was a *crude* FSH guess (0.45/0.85), now
**aligned to the cited experimental FRTTB**. The `reproduce_demo` placeholder is the
opposite extreme — root *too high* (0.30) and HI *non-physical* (0.07).

> Aligning `growth_rice`'s FSH to the experimental FRTTB shifts the short-chain shoot BAFs
> (~+30 % PFBA straw/grain; long chains <5 %) but leaves `reproduce_demo` untouched (it
> uses the placeholder, RMSE stays 0.029).

### 3.1 The experimental coefficient, then three ways to *force* the literature ratio

**The data-grounded answer is the experimental FRTTB itself** (now used by both drivers →
root:shoot 0.05–0.07, no tuning). To go further and *force* exactly the literature ratio
(~0.10), three options exist — but note a **root-inclusive per-organ biomass *time series*
for the target system is a data gap** (field rice studies omit roots; the closest,
`10.3389/fpls.2021.713814`, gives root + *total* shoot by stage, no organ split), so A is
blocked and B/C are tuned to the *ratio*:

| | method | how | artificiality | status |
|---|---|---|---|---|
| **A** | measured biomass | drive `M(t)` from a measured root-inclusive per-organ series (`measured_biomass.py`) | lowest | **blocked — data gap** |
| **B** | `target_root_shoot=` | scale the root *assimilate-partitioning* fraction `(1−FSH)` (solve for the factor), re-integrate; shoot renormalised, split preserved | medium (a biological allocation parameter) | implemented (opt-in) |
| **C** | `root_shoot=` | post-hoc constant rescale of the root *output* trajectory | highest (multiplies the output) | implemented (opt-in) |

Both **B and C land root:shoot ≈ 0.10**; **B preserves the harvest index better** (HI
0.50 vs C's 0.48) and is the recommended of the two. Neither is the default (`None`) —
forcing the exact ratio is a deliberate, documented choice on top of the experimental
FRTTB. **The principled fix is A**, which needs a root-inclusive maturity biomass series
for a comparable system (still outstanding).

## 4. Calibration coupling — the key finding

`validation/root_shoot_biomass_sensitivity.py` holds the **W2 transport fit**
(`params/parameters.json`) fixed and swaps the biomass driver, scoring the Yamazaki
root/straw/grain BAF (log10 RMSE):

| biomass | root:shoot | HI | W2 RMSE |
|---|---|---|---|
| reproduce placeholder | 0.30 | 0.07 | **0.029** |
| growth_rice default | 0.035 | 0.51 | 0.305 |
| growth_rice R/S=0.10 (lit) | 0.10 | 0.48 | 0.255 |
| growth_rice R/S=0.15 (lit) | 0.15 | 0.46 | 0.251 |
| growth_rice R/S=0.30 | 0.30 | 0.40 | 0.258 |

**The celebrated RMSE 0.029 reproduction is attained only with the non-physical
placeholder biomass** (root:shoot 0.30 *and* HI 0.07). With any realistic HI (~0.5)
the same W2 parameters give RMSE ~0.25–0.31 **regardless of root:shoot** — i.e. it is
the unrealistically low HI (almost no grain), not just the root, that the fit leans
on. The transport fit (`f_xy_W2fit`, `L_Ph_W2fit`, `kappa_d_W2fit`) is therefore
**entangled with a non-physical biomass** and is not valid on a literature-consistent
one *without re-fitting*.

### 4.1 Re-fit on the realistic biomass (`validation/refit_realistic_biomass.py`)

Re-fitting per-congener `f_xy`/`L_Ph`/`kappa_d` to the *same* Yamazaki BAF on the
literature biomass (root:shoot 0.10, HI 0.53; drivers otherwise identical to
reproduce_demo) **restores the reproduction**:

| biomass / parameters | overall log10 RMSE |
|---|---|
| placeholder, W2 fit (reproduce_demo) | 0.029 |
| realistic, W2 fit **un-changed** | ~0.31 |
| realistic, **re-fitted** | **0.103** (≈ **0.017** excl. near-MQL PFDoDA) |

8 of 11 congeners re-fit to RMSE ≈ 0 — realistic biomass is fully fittable; the
"collapse" was only the un-refit parameters. Findings (`params/refit_realistic_biomass.csv`):

- **`f_xy` is systematically LOWER than the W2 fit (~0.4–0.6×)**: PFOA 0.026→0.015,
  PFOS 0.142→0.077, PFDA 0.082→0.061. The root-richer, grain-heavier realistic biomass
  needs less transpiration-stream loading to match straw/grain.
- **The non-monotone (U-shaped) `f_xy` across chain length PERSISTS** (declines to C7,
  rises C8→C12) — so the long-chain shoot rise is **not a biomass artifact**; it
  survives a literature biomass (cf. `docs/fxy_longchain_lipid_exploration.md`).
- **PFDoDA stays unfittable** (RMSE 0.337; near-MQL outlier with f_xy/L_Ph at bounds) —
  the same outlier flagged in reproduce_demo.

OVERRIDE-only: written to `params/refit_realistic_biomass.csv`; `params/parameters.json`
is unchanged pending a decision to promote it.

## 5. Implications for the burden / "leaf dominates" question

- The leaf-dominant burden seen in the default `simulate` (4-pool) is partly the
  documented **leaf-sink runaway** (`nstem_leaf` redistributes it to grain≈leaf≈stem)
  and partly the **root mass being too low** (so root burden is under-counted,
  especially for root-dominated PFOS / long chains). See
  `validation/root_shoot_biomass_sensitivity.py` and the burden sensitivity in chat.
- Correcting root:shoot to ~0.10 roughly triples root burden share (PFOS root
  54%→78%) and ~doubles total uptake (bigger root = more uptake surface).

## 6. Follow-up (open)

1. **Re-fit the transport parameters** (`f_xy`, `L_Ph`, `kappa_d`) on a literature-
   consistent biomass — **DONE** (§4.1, `validation/refit_realistic_biomass.py`):
   reproduction restored (RMSE ~0.017 excl. PFDoDA), `f_xy` ~0.4–0.6× the W2 values,
   U-shape persists. Open decision: whether to **promote** `refit_realistic_biomass.csv`
   to the canonical `parameters.json` (would change all downstream BAFs; needs a
   provenance/versioning step and a re-run of dependent validations).
2. **Close the root data gap**: a maturity, root-inclusive biomass for the target
   system to pin root:shoot (currently a literature range 0.08–0.15).
3. Decide whether to make a literature root:shoot the `growth_rice` default (kept
   `None`/unchanged for now to preserve reproducibility until the re-fit in #1).
