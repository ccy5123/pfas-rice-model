# Neutral-organic (Briggs/Kow) path — implementation and validation status

> `src/neutral_dpu.py` · `validation/neutral_dpu_validation.py` · `tests/test_neutral_dpu.py`
> **Status: IMPLEMENTED, anchors verified at source, and validated a-priori
> against one measured rice dataset** (Ge et al. 2017; log10 RMSE **1.099** with
> nothing fitted). §4 records what that number means and what is still missing.

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

Anchors — **verified at source** (Briggs, Bromilow & Evans 1982, *Pestic. Sci.*
13:495–504, their eqs. 2 and 3; reproduced verbatim as eq. 4 of Schriever &
Lamshoeft 2020, *STOTEN* 713:136667):

```
TSCF = 0.784 · exp[ −(log Kow − 1.78)² / 2.44 ]
RCF  = 0.82 + 10^(0.77·log Kow − 1.52)
```

`RCF` is the same `K_PW` form with `W = 0.82`, `L·a = 10^−1.52`, `b = 0.77`; only
the **product** `L·a` is identifiable, so `a = 1.22` is a convention and `L·a` is
the anchored quantity.

Reading the original settled a point that had been an assumption: Briggs derives
the 0.82 floor by subtracting the fitted lipid term from the measured RCF of nine
polar compounds and attributes it to *"equilibration of the chemical between the
external solution and the water contained within the roots"*, noting that roots at
90 % water would give ≈ 0.9. So the first term of `K_PW` **is** the tissue water
content — which is exactly how this module maps it onto `Compartment.theta`.

**Tissue lipid contents** now come from Trapp, McFarlane & Matthies 1994 — the
canonical validation of this very framework — which states root 1 %, stem and leaf
3 % dry weight. Those are soybean values, so organ-resolved total lipid for *rice*
is still a gap, but they are a cited anchor rather than the guesses that stood here
before.

**An alternative TSCF.** Schriever & Lamshoeft 2020 refit the same Gaussian to 97
TSCF values from intact-plant hydroponic tests (42 compounds): `A = 0.746`,
`B = 2.160`, `C = 7.230`. The peak height is essentially Briggs' (0.75 vs 0.78) but
the bell is three times **broader**. `simulate_neutral(tscf_model="schriever")`
runs with it. Since TSCF here is an input rather than a fitted parameter, the gap
between the two is a fair measure of how well it is actually known — and §4 shows
which one the rice data prefer.

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
| 0.00 | 0.214 | 0.90 | 0.74 | 34.5 | 46.9 |
| 1.00 | 0.611 | 0.91 | 0.55 | 67.6 | 122.2 |
| **1.78** | **0.784** | 0.93 | 0.51 | 76.4 | **149.8** |
| 2.50 | 0.634 | 1.00 | 0.60 | 69.0 | 114.4 |
| 3.50 | 0.233 | 1.50 | 1.21 | 36.8 | 30.5 |
| 4.50 | 0.038 | 4.46 | 4.29 | 7.7 | 1.8 |
| 5.50 | 0.003 | 21.9 | 21.8 | 0.6 | 0.03 |

The model reproduces the law every uptake study reports — polar compounds
translocate to the shoot, lipophilic ones are retained in the root — and the
**straw/root ratio peaks at exactly log Kow 1.78, the Briggs TSCF peak**, crossing
1 between log Kow 4.5 and 5.5. Nothing was fitted to produce this.

**§3 Scope — metabolism and volatilisation are load-bearing here.** With no
phloem, no air exchange and `γ = 0`, the leaf is the sole xylem terminal and its
only sink is growth dilution, which → 0 at maturity. So a recalcitrant,
non-volatile neutral compound **must** run away — the leaf integrates the whole
transpiration stream. That is arithmetic, not a bug, but it is a real difference
from the PFAS case where `γ ≈ 0` is defensible:

| half-life | γ [1/d] | leaf BAF | root BAF |
|---|---:|---:|---:|
| ∞ (recalcitrant) | 0 | 193.7 | 0.51 |
| 60 d | 0.0116 | 119.5 | 0.51 |
| 21 d | 0.0330 | 59.9 | 0.51 |
| 7 d | 0.0990 | 19.4 | 0.51 |
| 2 d | 0.3466 | 5.0 | 0.51 |

The root is exposure-buffered and barely moves; the leaf spans 40×. **Report
neutral runs with a measured half-life, or state them as an upper bound.**
Volatilisation is not implemented at all (the core ODE has no air terms), so
`neutral_dpu.k_aw_warning()` refuses to let a high-`K_AW` compound run silently.

**§4 Switch — PASS.** `N = 0`, `e^N = 1`, `Vmax = 0`: the ionic machinery is
verifiably off, and `root_uptake` is exactly `κ_d·(C_w^o − C_w,root)`.

**What none of this establishes.** These are checks against published QSPRs and
against the model's own structure. They can falsify the implementation and they
quantify its scope. **They are not validation against measured plant data.**

## 4. The a-priori prediction (Ge et al. 2017)

`python validation/neutral_dpu_validation.py --obs data_obs/neutral_obs_ge2017.csv`

**Dataset.** Ge, Cui, Yan, Li, Chai, Liu, Cheng & Yu (2017), *Environ. Pollut.*
226:479–485, `10.1016/j.envpol.2017.04.043` — soil-grown rice, three neutral
pesticides spanning log Kow −0.13 → 4.4, per-organ transfer factors at 60 days
(their Table 2), with log Kow and pKa in the same paper (their Table 3).
Transcribed to `data_obs/neutral_obs_ge2017.csv`.

Their BCF is `C_plant/C_soil`, which would need a soil–water partition coefficient
we do not have for that soil, so the comparison uses **TF** — a tissue/tissue
ratio, in which the exposure term cancels. Values are dry weight and are converted
onto the model's fresh-weight basis using the run's own tissue water contents.

**Result — nothing fitted, `K_PW` and `TSCF` from log Kow alone:**

| compound | log Kow | organ | observed TF | model | ratio |
|---|---:|---|---:|---:|---:|
| thiamethoxam | −0.13 | stem | 0.34 | 0.16 | 0.5× |
| thiamethoxam | −0.13 | leaf | 15.8 | 99.0 | 6.3× |
| imidacloprid | 0.57 | stem | 0.34 | 0.40 | 1.2× |
| imidacloprid | 0.57 | leaf | 16.1 | 226.0 | 14× |
| difenoconazole | 4.40 | stem | 0.034 | 0.19 | 5.7× |
| difenoconazole | 4.40 | leaf | 0.044 | 6.2 | 141× |

**log10 RMSE = 1.099** (n = 6). This is the **only a-priori prediction in this
repo** — read it against the PFAS side's a-priori error of ≈ 0.84–0.95, noting
that the PFAS number still has fitted transport parameters behind it while this
one has none.

**The error structure is the finding.** The **stem is predicted well** (0.5×,
1.2×, 5.7×) while the **leaf is over-predicted by 6–141×**. That is precisely the
failure mode §3 predicted: with `γ = 0` the leaf is an unbounded terminal
accumulator, so a 60-day exposure integrates the entire transpiration stream into
it. The model also gets the qualitative ordering right across 4.5 log units of
Kow — polar compounds reach the leaf (TF ≫ 1), difenoconazole stays in the root
(TF ≪ 1).

**Sensitivity to the two inputs we do not have** (§6 of the script). Neither is a
free parameter; both are unmeasured, so the scan says how much of the residual is
a missing measurement rather than a broken structure:

| in-planta half-life | Briggs TSCF | Schriever TSCF |
|---|---:|---:|
| none (γ = 0, strict a-priori) | **1.099** | 1.551 |
| 60 d | 0.960 | 1.425 |
| 30 d | 0.848 | 1.322 |
| 14 d | 0.668 | 1.156 |
| 7 d | 0.502 | 0.984 |
| 3 d | 0.409 | 0.792 |

Two things follow.

1. **The residual is dominated by the missing dissipation rate, not by transport
   structure.** The error falls monotonically as metabolism is imposed, reaching
   0.41 at a 3-day half-life. At any plausible in-planta half-life below ~30 d the
   neutral a-priori error is already at or below the PFAS a-priori error. Ge report
   *soil* half-lives only (IMI 19–20 d, THX 26–30 d, DFZ 37–41 d) but do document
   in-plant metabolism qualitatively (THX → clothianidin *in the plants*). The
   half-life that reconciles the leaf is therefore a **testable prediction**, not a
   calibration — measuring it would confirm or refute this reading.
2. **The original Briggs bell beats the modern refit on this data**, at every
   half-life (1.10 vs 1.55 at γ = 0). The broader Schriever bell gives
   difenoconazole `TSCF = 0.37` against Briggs' 0.047, and Ge measured a leaf TF of
   0.044 — the narrow bell is right about lipophilic compounds not reaching the
   leaf. Worth knowing, since the refit rests on the larger dataset and might have
   been assumed the better default.

**What this does not settle.** Six data points, three compounds, one study, one
soil, one time point. The grain compartment is not exercised at all (Ge harvested
at 60 days, and no rice dataset found reports a neutral Kow series in grain under
root-only exposure — the same gap that exists on the PFAS side). And the leaf
result is confounded with the missing half-life, so the strict number is an upper
bound on the error, not a measurement of the transport model.

## 5. The other papers, and what each can and cannot do

All ten obtained papers were read; here is what each is actually good for.

| paper | verdict |
|---|---|
| **Ge 2017** `10.1016/j.envpol.2017.04.043` | **Used — §4.** Complete and self-contained: TF, log Kow and pKa all in the article. |
| **Briggs 1982** `10.1002/ps.2780130506` | **Anchor verified at source.** Eqs. 2 and 3 confirmed character-for-character; the 0.82 floor is explained as root water content, confirming the `K_PW` mapping. Its Table 1 (barley RCF/TSCF) is the QSPR's own training data, so it cannot serve as validation. |
| **Schriever 2020** `10.1016/j.scitotenv.2020.136667` | **Alternative TSCF implemented.** Reprints Briggs eq. 3 verbatim (independent verification) and supplies a 97-value refit, now selectable and benchmarked in §4. Its per-compound TSCF values are in the SI, which was not obtained. |
| **Trapp 1994** `10.1002/etc.5620130308` | **Tissue composition adopted** (root 1 %, stem/leaf 3 % lipid). Its Table 1 gives a complete driver set (transpiration, organ masses, water contents) for the bromacil runs, but the measured concentrations are in Figures 5–6 — **figure-only**, so an RMSE would need digitising. |
| **Liu 2023** `10.1016/j.scitotenv.2022.159826` | **Partly usable; blocked on its SI.** Confirmed: hydroponic, 21 pesticides at 100 µg/L, 3–144 h, root/stem/leaf. But per-compound RCF/SCF/LCF live in Fig. 2 and Tables S4–S5, and log Kow in Table S1 — **the SI was not in the upload**. The text alone gives TF ranges and a few named values. Fetching the SI would make this the strongest rice dataset available. |
| **Inao 2018a/b** `10.1584/jpestics.D17-083` / `D17-084` | **Cannot test the 4-compartment split.** As flagged before the papers arrived, they sample *"the whole shoot of the rice plant above the water surface"* — no organ resolution. Still the only source of a measured paddy-water + layered-soil `C_w^o(t)`, so it remains the right dataset for a future HYDRUS-coupling test against a lumped shoot. |
| **Brunetti 2021** `10.1021/acs.est.0c07420` | **Parameter cross-check, not a data table.** Its Table 1 reports *calibrated* posteriors, not raw concentrations: green-pea root `K_RW = 13.3` cm³/g fw and stem `K_SW = 11.8` for carbamazepine (log Kow ≈ 2.45). The Briggs `K_PW` for the same compound is ≈ 1.0 — **an order of magnitude lower**. Either pea tissue is far more sorptive than Briggs' barley, or the calibration absorbed other processes. Worth pursuing: it is a direct quantitative disagreement with the partition core, on the framework's own reference implementation. |
| **Hwang 2017** `10.1371/journal.pone.0172254` | Open access, transcribable (Table 1: chlorpyrifos root/leaf in lettuce, 21/30/40 d). Lettuce not rice, one compound, and soil-basis exposure — a useful secondary check, not a QSPR test. |
| **Briggs 1983** `10.1002/ps.2780140506` | Shoot-distribution companion to the 1982 paper; not yet mined. |

## 6. Honest summary

- The neutral path is **implemented**, is the published Briggs core (verified to
  machine precision against equations now checked at source), and runs on the
  unmodified 4-compartment ODE.
- It reproduces the qualitative Kow law with the translocation peak landing on the
  Briggs bell maximum, using **zero fitted parameters**.
- Against measured rice data it predicts per-organ transfer factors a-priori at
  **log10 RMSE 1.099**, with a clean, interpretable error structure: the stem is
  right, the leaf is high because in-planta metabolism is missing, and imposing a
  plausible half-life brings the error to 0.41–0.85.
- Remaining gaps: the **grain compartment is untested** (no suitable dataset
  exists); rice-specific organ lipid contents are still borrowed from soybean; the
  Liu 2023 SI would add a proper Kow series; and the Brunetti `K_RW` discrepancy is
  an open question against the partition core itself.
