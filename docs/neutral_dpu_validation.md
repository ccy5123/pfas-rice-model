# Neutral-organic (Briggs/Kow) path — implementation and validation status

> `src/neutral_dpu.py` · `validation/neutral_dpu_validation.py` · `tests/test_neutral_dpu.py`
> **Status: IMPLEMENTED and internally consistent with the published QSPRs;
> NOT YET VALIDATED against measured plant data** (blocked on data access — §4).

---

## 1. Why a neutral path exists in a PFAS repo

Every PFAS result here is entangled with PFAS-specific parameters that had to be
**fitted** — `f_xy`, `k_seq`, the lipid conductances. That is exactly why the honest
verdicts across `CLAUDE.md` §6 read *"reproduction, not prediction"*, and why the
a-priori predictive error is log10 RMSE ≈ 0.84–0.95 rather than the 0.03 of a
saturated fit. The framework's skeleton — four compartments, xylem advection,
growth dilution, terminal accumulation — has therefore never been tested apart
from the ionic extension bolted onto it.

A neutral organic is the one case where it can be. Its partitioning (`K_PW`) and
its root→shoot loading (`TSCF`) are **both fixed from outside the model** by
QSPRs on log Kow, with nothing left to tune. A failure on the neutral path would
be a failure of the skeleton; a failure on the PFAS side alone is not.

The neutral base was **derived** in `docs/dpu_model_summary_corrected.tex` but
never implemented: `Compound.fn` was a field pinned at 0.0 that appeared in no
equation, the neutral membrane term existed only as a comment
(`# optional neutral passive term (negligible for PFAS): cmpd.fn * ...`), and
`data_obs/` holds only PFAS.

## 2. The implementation — no new ODE

A neutral compound is the **existing** 4-compartment model with the ionic
machinery removed by physics rather than by special-casing:

| ionic feature | neutral setting | consequence |
|---|---|---|
| valence `z = -1` → `N = zEF/RT = 4.67` | **`z = 0` → `N = 0`** | GHK factor → 1, exclusion factor `e^N` 107 → **1**; the membrane term degenerates *exactly* to passive Fickian diffusion `κ_d·(C_w^o − C_w,root)` |
| saturable carrier (needed to beat exclusion) | **`Vmax = 0`** | neutral molecules cross unaided |
| basis-A binding `θ + (1−θ)(f_prot K_prot + f_PL K_PL + f_cw K_cw)` | `f_prot = f_cw = 0`, `K_PL = a·Kow^b` | the same call returns `W + L·a·Kow^b` = the Trapp/Briggs **K_PW**, term for term |
| fitted `f_xy` | **`f_xy = TSCF(log Kow)`** | the Briggs bell, computed not fitted |
| phloem-fed grain | **phloem OFF** | the base states plainly that "no dissociation, pH-dependent speciation, membrane electrical potential, ion-trap, or **phloem transport** is included" — the phloem is an addition of the *ionisable* extension |

Anchors (already in-repo, `docs/theory_anchor.tex` eqs. briggsT/briggsR; Briggs,
Bromilow & Evans 1982, *Pestic. Sci.* 13:495–504):

```
TSCF = 0.784 · exp[ −(log Kow − 1.78)² / 2.44 ]
RCF  = 0.82 + 10^(0.77·log Kow − 1.52)
```

`RCF` is the same `K_PW` form with `W = 0.82`, `L·a = 10^−1.52`, `b = 0.77`; only
the **product** `L·a` is identifiable, so `a = 1.22` is a convention and `L·a` is
the anchored quantity.

## 3. What has been checked (and what that is worth)

`python validation/neutral_dpu_validation.py`

**§1 Partition adapter — PASS to machine precision (max rel. err 1.3e-16).** The
core's `binding_factors`, given the neutral composition, reproduces Briggs' RCF
exactly across log Kow −1…6. This tests the *adapter*, not Briggs: it establishes
that the neutral path really is the published partition core and not a re-fit of it.

**§2 Kow signature — PASS, with zero fitted parameters.** Scanning log Kow −1…5.5
on the measured forcings:

| log Kow | TSCF | K_PW root | root BAF | straw BAF | straw/root |
|---:|---:|---:|---:|---:|---:|
| −1.00 | 0.033 | 0.90 | 0.87 | 6.8 | 7.8 |
| 0.00 | 0.214 | 0.91 | 0.74 | 34.5 | 46.5 |
| 1.00 | 0.611 | 0.96 | 0.58 | 67.6 | 115.8 |
| **1.78** | **0.784** | 1.13 | 0.62 | 76.4 | **123.2** |
| 2.50 | 0.634 | 1.72 | 1.03 | 68.9 | 66.6 |
| 3.50 | 0.233 | 5.74 | 4.61 | 36.7 | 8.0 |
| 4.50 | 0.038 | 29.4 | 28.2 | 7.5 | 0.27 |
| 5.50 | 0.003 | 168.6 | 168.0 | 0.5 | 0.00 |

The model reproduces the law every uptake study reports — polar compounds
translocate to the shoot, lipophilic ones are retained in the root — and the
**straw/root ratio peaks at exactly log Kow 1.78, the Briggs TSCF peak**, crossing
1 near log Kow 4.5. Nothing was fitted to produce this.

**§3 Scope — metabolism and volatilisation are load-bearing here.** With no
phloem, no air exchange and `γ = 0`, the leaf is the sole xylem terminal and its
only sink is growth dilution, which → 0 at maturity. So a recalcitrant,
non-volatile neutral compound **must** run away — the leaf integrates the whole
transpiration stream. That is arithmetic, not a bug, but it is a real difference
from the PFAS case where `γ ≈ 0` is defensible:

| half-life | γ [1/d] | leaf BAF | root BAF |
|---|---:|---:|---:|
| ∞ (recalcitrant) | 0 | 193.6 | 0.62 |
| 60 d | 0.0116 | 119.5 | 0.62 |
| 21 d | 0.0330 | 59.8 | 0.62 |
| 7 d | 0.0990 | 19.4 | 0.62 |
| 2 d | 0.3466 | 5.0 | 0.61 |

The root is exposure-buffered and barely moves; the leaf spans 40×. **Report
neutral runs with a measured half-life, or state them as an upper bound.**
Volatilisation is not implemented at all (the core ODE has no air terms), so
`neutral_dpu.k_aw_warning()` refuses to let a high-`K_AW` compound run silently.

**§4 Switch — PASS.** `N = 0`, `e^N = 1`, `Vmax = 0`: the ionic machinery is
verifiably off, and `root_uptake` is exactly `κ_d·(C_w^o − C_w,root)`.

**What none of this establishes.** These are checks against published QSPRs and
against the model's own structure. They can falsify the implementation and they
quantify its scope. **They are not validation against measured plant data.**

## 4. The missing piece: measured data — and why it is missing

§5 of the validation script is the comparison harness, and it is **inert until a
measured table is supplied** (`--obs`). The schema is
`data_obs/neutral_obs_template.csv` (which deliberately contains **no data** — its
placeholder rows are refused by the loader so it can never be mistaken for
measurements). The harness converts dry-weight observations onto the model's
fresh-weight basis using the run's own tissue water contents, since `(1−θ)`
differs by tissue and therefore does **not** cancel in a tissue/root ratio — the
same trap that produced the Tang 2026 fresh-vs-dry artifact.

**Why no dataset is included: the session that built this could not obtain one.**
The execution environment's network policy blocks all outbound HTTPS except a
small allowlist. Every publisher, PMC, J-STAGE and Crossref host returns a 403 at
the proxy CONNECT stage. Search worked; fetching did not. **No numeric value in
this repo has been transcribed from a paper that was not already in it**, and
none should be added from search snippets.

### Candidate datasets, ranked (from search metadata only — verify before use)

| # | dataset | why | limitation |
|---|---|---|---|
| 1 | **Liu et al. 2023**, *Sci. Total Environ.* 858:159826, `10.1016/j.scitotenv.2022.159826` | the only rice dataset found with all three of: a compound series spanning ~5 log Kow units, **time-resolved** uptake, and per-organ **root/stem/leaf**. Subcellular fractionation (cell wall / organelle / cell sap) additionally maps onto the binding decomposition | hydroponic (no soil side), **no grain**, seedlings (no growth dilution); the sulfonylureas in it are weak acids and must be dropped |
| 2 | **Inao et al. 2018**, *J. Pestic. Sci.* 43(2):132–141, `10.1584/jpestics.D17-083` (+ companion model paper) | **fully open access** (J-STAGE + PMC6140689); a real paddy with measured **paddy-water and layered-soil time series** — i.e. the `C_w^o(t)` driver the HYDRUS coupling produces, so it can test the coupling, not just the plant | only 2 compounds (cannot test the QSPR); **check first whether "rice plant" is organ-resolved or whole-plant** — that decides whether it can test a 4-compartment model at all |
| 3 | **Ge et al. 2017**, *Environ. Pollut.* 226:479–485, `10.1016/j.envpol.2017.04.043` | small but on-model: soil-treated pots, imidacloprid / thiamethoxam / difenoconazole — a genuine hydrophilic↔lipophilic contrast that should reproduce the §2 ordering | 3 compounds, root+leaf only, no grain, dry-weight basis |
| 4 | **Trapp, McFarlane & Matthies 1994**, *Environ. Toxicol. Chem.* 13:413–422, `10.1002/etc.5620130308` | the canonical 4-compartment DPU validation this framework descends from (bromacil in soybean, reportedly with no parameter adjustment) | soybean not rice; 1994, so values are likely figures-only and need digitising |
| 5 | **Brunetti et al. 2021**, *ES&T* 55:2991–3000, `10.1021/acs.est.0c07420` (PMC8023655) | literally the HYDRUS+DPU coupling this repo mirrors; carbamazepine is genuinely neutral and all four compartments incl. fruit are measured | green pea not rice; one compound, so it tests dynamics not the QSPR |

**Recommendation.** Lead with **Liu 2023** if journal access is available (it is
the only one that can falsify the QSPR and the dynamics simultaneously), and take
**Inao 2018** regardless — it is open access and it is the only candidate that
tests the soil coupling. Reproducing **Brunetti 2021** would be the cleanest
framework sanity check before trusting any rice result.

### How to run it once you have a table

```bash
# transcribe into the schema, then:
python validation/neutral_dpu_validation.py --obs data_obs/<your>.csv
```

The reported log10 RMSE is a **genuine a-priori prediction** — the only one in
this repo, since `K_PW` and `TSCF` come from log Kow alone and nothing is fitted
to the table. Read it against the PFAS side's a-priori error of ≈ 0.84–0.95.

## 5. Honest summary

- The neutral path is **implemented** and is the published Briggs core, verified
  to machine precision, running on the unmodified 4-compartment ODE.
- It reproduces the qualitative Kow law, with the translocation peak landing on
  the Briggs bell maximum, using **zero fitted parameters**.
- It quantifies its own scope: no phloem (faithful to the base), no air exchange
  (so volatile compounds are flagged, not silently run), and an unbounded leaf
  without metabolism (so a half-life is required, not optional).
- **It has not been compared to a single measured plant concentration.** Until it
  is, the correct claim is "the neutral backbone is implemented and internally
  consistent", never "the backbone is validated".
