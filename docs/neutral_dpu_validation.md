# Neutral-organic (Briggs/Kow) path — implementation and validation status

> `src/neutral_dpu.py` · `validation/neutral_dpu_validation.py` · `tests/test_neutral_dpu.py`
> **Status: IMPLEMENTED, anchors verified at source, and validated a-priori
> against two measured rice datasets** — root partition across 14 compounds and
> 5 log-Kow units at log10 RMSE **0.281**, per-organ transfer factors at **0.783**,
> both with nothing fitted. §4–§5 record what those numbers mean and what is left.

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
3 %, on a **fresh-weight** basis (the same basis as `W`, as the `K_PW` sum
requires; see the correction at the end of §4). Those are soybean values, so
organ-resolved total lipid for *rice* is still a gap, but they are a cited anchor
rather than the guesses that stood here before.

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
| 0.00 | 0.214 | 0.91 | 0.75 | 34.5 | 46.3 |
| 1.00 | 0.611 | 0.97 | 0.59 | 67.6 | 114.1 |
| **1.78** | **0.784** | 1.19 | 0.65 | 76.4 | **117.3** |
| 2.50 | 0.634 | 1.93 | 1.16 | 69.0 | 59.5 |
| 3.50 | 0.233 | 6.94 | 5.58 | 36.8 | 6.6 |
| 4.50 | 0.038 | 36.5 | 35.1 | 7.7 | 0.22 |
| 5.50 | 0.003 | 210 | 210 | 0.5 | 0.00 |

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
| ∞ (recalcitrant) | 0 | 193.2 | 0.65 |
| 60 d | 0.0116 | 119.1 | 0.65 |
| 21 d | 0.0330 | 59.5 | 0.65 |
| 7 d | 0.0990 | 19.1 | 0.65 |
| 2 d | 0.3466 | 4.7 | 0.64 |

The root is exposure-buffered and barely moves; the leaf spans 40×. **Report
neutral runs with a measured half-life, or state them as an upper bound.**

**§3b Volatilisation — the second leaf sink, now implemented** (`src/plant_air.py`;
was the largest structural gap in the neutral base). The leaf has exactly two
sinks besides growth dilution, metabolism and volatilisation, and until now only
the first existed — so every volatile compound was an upper bound *by
construction*. The equation set was already written out in
`docs/dpu_model_summary_corrected.tex` `sec:permeability` and is implemented from
it: cuticle (`eq:Pc`), air boundary layer (`eq:Pair`) and aqueous layer
(`eq:Paqua`) in **series** (`eq:Pctot`), in **parallel** with the transpiration-
linked stomatal conductance (`eq:Ps`, `eq:Pp`), driving volatilisation (`eq:Qvol`)
against gaseous uptake (`eq:Qgas`) through the plant–air partition `K_PA`
(`eq:Kpa`). The derivation's own assumptions are kept: **no volatilisation from
roots**, **stem via cuticle only**, **leaf and grain via cuticle + stomata**.

Scanning `K_AW` at fixed log Kow 2.42, MW 131.4, clean air:

| `K_AW` | leaf t½ (d) | leaf BAF | grain BAF | root BAF |
|---|---:|---:|---:|---:|
| **0 (PFAS)** | **∞** | 177 | 8.81 | 1.058 |
| 1e−5 | 6.8e4 | 177 | 8.81 | 1.058 |
| 1e−4 | 787 | 169 | 8.78 | 1.058 |
| 1e−3 | 8.0 | 23.8 | 6.92 | 1.058 |
| 1e−2 | 0.080 | 0.251 | 0.203 | 1.058 |
| 1e−1 | 0.0008 | 0.0025 | 0.0021 | 1.058 |

Three things to read off it. (i) The `K_AW = 0` row is the **PFAS constraint**:
`P_air` and `P_S` are both proportional to `K_AW`, and `P_air` sits in *series*
with the cuticle, so the air pathway is **structurally absent**, not numerically
small — and the core keeps `air=None` as its default, so the terms are never even
evaluated. `reproduce_demo` stays at log10 RMSE **0.029**, verified. (ii) The
**root is invariant** across the whole ladder, as the derivation requires. (iii)
The ladder independently **checks the pre-existing warning threshold**: the
judgement-call `K_AW > 1e−4` in `k_aw_warning` sits just below where the physics
puts the crossover (t½ 787 d at 1e−4, negligible against a season; 8 d at 1e−3,
dominant), so it errs toward flagging early — the safe direction.

The warning therefore changed job rather than disappearing: with air exchange off
it still marks a volatile result as an upper bound, but it now names the remedy
(`simulate_neutral(..., air=True)`) instead of saying the process is not modelled.

**The honest limit is the surface area, not the equations.** The flux scales
linearly with each tissue's specific surface `S` [m²/kg], and `S` entered this repo
as a leaf/grain **ratio** for splitting the xylem stream — only the ratio was ever
load-bearing, so the absolute values (leaf 20, grain 2) have never been calibrated
as areas. Treat an absolute volatilisation magnitude as order-of-magnitude until
they are measured; the `K_AW`/Kow *dependence* is the derivation's.
`AirExchange(S=...)` overrides them. The shipped stem `S` is **0**, so the stem's
cuticular term is inert until a real stem area is supplied — pinned by a test so it
cannot silently become a hidden number. Particle deposition (`eq:Qdep`) remains
unimplemented; `f_particle` only excludes the particle-bound share from the gaseous
uptake, as `eq:Qgas` specifies.

**§4 Switch — PASS.** `N = 0`, `e^N = 1`, `Vmax = 0`: the ionic machinery is
verifiably off, and `root_uptake` is exactly `κ_d·(C_w^o − C_w,root)`.

**What none of this establishes.** These are checks against published QSPRs and
against the model's own structure. They can falsify the implementation and they
quantify its scope. **They are not validation against measured plant data.**

## 4. The a-priori predictions

Two independent rice datasets, both run with **nothing fitted** — `K_PW` and
`TSCF` follow from log Kow alone.

### 4a. Root partition — Liu et al. 2023 (14 compounds, 5 log-Kow units)

`python validation/neutral_dpu_validation.py --obs data_obs/neutral_obs_liu2023.csv`

*Sci. Total Environ.* 858:159826, `10.1016/j.scitotenv.2022.159826`. Hydroponic
rice; the subcellular experiment (1 mg/L, 72 h, 25-day seedlings).

**Reconstruction.** Per-compound tissue concentrations appear only in figures, but
the SI gives subcellular concentrations per tissue (Table S4) and the mass fraction
of each subcellular component (Table S3) — and those fractions sum to exactly 1.00
for every tissue. So the total is their mass-weighted sum, `C = Σ_f w_f·C_f`:
arithmetic on published numbers, not digitising. It is checkable against the
paper's own text — the reconstructed leaf/stem ratios for the neonicotinoids fall
in 1.4–3.2 against the text's "TF_L/S … >2 (2.2–5.2)", and give "rest <1" for the
triazoles as the text also states. The seven **sulfonylureas are excluded**: they
are weak acids, ionised at the solution pH 5.6–5.8, and belong on the ionic path.

**Why root partition is the cleanest possible test:** it needs no assumption about
transpiration, exposure duration, or metabolism, because it is an *equilibrium* —
Briggs' own observation that roots reach "an equilibrium concentration that does
not change further with time". So this isolates the `K_PW` core, on rice, against
a QSPR fitted to barley.

| compound | log Kow | obs RCF | model | ratio |
|---|---:|---:|---:|---:|
| nitenpyram | −0.66 | 0.86 | 0.84 | 0.98 |
| dinotefuran | −0.55 | 1.00 | 0.83 | 0.83 |
| thiamethoxam | −0.13 | 1.13 | 0.77 | 0.68 |
| imidacloprid | 0.57 | 0.73 | 0.64 | 0.88 |
| clothianidin | 0.70 | 1.23 | 0.63 | 0.51 |
| acetamiprid | 0.80 | 0.44 | 0.61 | 1.39 |
| thiacloprid | 1.26 | 3.14 | 0.59 | 0.19 |
| triadimefon | 2.77 | 2.89 | 1.65 | 0.57 |
| myclobutanil | 2.94 | 7.21 | 2.13 | 0.30 |
| epoxiconazole | 3.58 | 8.04 | 6.46 | 0.80 |
| propiconazole | 3.72 | 9.32 | 8.36 | 0.90 |
| flusilazole | 3.87 | 15.4 | 11.0 | 0.72 |
| hexaconazole | 3.90 | 8.32 | 11.7 | 1.40 |
| difenoconazole | 4.40 | 44.9 | 29.2 | 0.65 |

**The reconstructed shoot values, recorded but not used.** The same `S3 × S4`
reconstruction also yields stem and leaf concentrations. They are *not* in the
comparison table — a 72-hour exposure of 25-day seedlings is not comparable to the
season-long growing run the model drives, and root partition is the only endpoint
that transfers because it is an equilibrium. They are kept here so the
reconstruction does not have to be redone if a short-exposure driver is ever built:

| compound | log Kow | root | stem | leaf | stem/root | leaf/root |
|---|---:|---:|---:|---:|---:|---:|
| nitenpyram | −0.66 | 0.86 | 0.50 | 1.15 | 0.58 | 1.34 |
| dinotefuran | −0.55 | 1.00 | 0.33 | 0.88 | 0.33 | 0.88 |
| thiamethoxam | −0.13 | 1.13 | 1.84 | 4.65 | 1.63 | 4.12 |
| imidacloprid | 0.57 | 0.73 | 1.96 | 2.66 | 2.68 | 3.63 |
| clothianidin | 0.70 | 1.23 | 0.41 | 0.89 | 0.33 | 0.72 |
| acetamiprid | 0.80 | 0.44 | 0.94 | 2.99 | 2.12 | 6.77 |
| thiacloprid | 1.26 | 3.14 | 2.75 | 4.62 | 0.88 | 1.47 |
| triadimefon | 2.77 | 2.89 | 2.74 | 1.15 | 0.95 | 0.40 |
| myclobutanil | 2.94 | 7.21 | 13.58 | 14.38 | 1.88 | 1.99 |
| epoxiconazole | 3.58 | 8.04 | 9.04 | 3.96 | 1.13 | 0.49 |
| propiconazole | 3.72 | 9.32 | 10.99 | 4.75 | 1.18 | 0.51 |
| flusilazole | 3.87 | 15.38 | 13.88 | 3.52 | 0.90 | 0.23 |
| hexaconazole | 3.90 | 8.32 | 11.04 | 6.36 | 1.33 | 0.76 |
| difenoconazole | 4.40 | 44.94 | 5.22 | 0.21 | 0.12 | 0.005 |

(Concentrations as tabulated in Table S4, i.e. mg/kg at a 1 mg/L exposure, so the
root column is numerically the RCF used above. Note these shoot ratios are *not*
the paper's own main-text TFs, which come from the separate 100 µg/L kinetics run
and differ by roughly 2× — a different experiment, not a contradiction.)

**log10 RMSE = 0.281** (n = 14). Eleven of fourteen fall within a factor of 1.5,
across a 50-fold range of measured RCF. The worst miss is thiacloprid (5.4× under).
The half-life sensitivity is **flat** (0.281 → 0.294 at a 3-day half-life), which is
exactly right for an equilibrium endpoint and is the counterpart of §4b, where it
is steep.

### 4b. Per-organ transfer — Ge et al. 2017 (3 compounds)

`python validation/neutral_dpu_validation.py --obs data_obs/neutral_obs_ge2017.csv`

*Environ. Pollut.* 226:479–485, `10.1016/j.envpol.2017.04.043`. Soil-grown rice,
per-organ TF at 60 days (their Table 2), log Kow and pKa in the same paper (Table
3). Their BCF is `C_plant/C_soil` and would need a soil–water coefficient we do not
have, so the comparison uses **TF**, in which the exposure term cancels.

| compound | log Kow | organ | observed TF | model | ratio |
|---|---:|---|---:|---:|---:|
| thiamethoxam | −0.13 | stem | 0.34 | 0.17 | 0.5× |
| thiamethoxam | −0.13 | leaf | 15.8 | 98.0 | 6.2× |
| imidacloprid | 0.57 | stem | 0.34 | 0.43 | 1.3× |
| imidacloprid | 0.57 | leaf | 16.1 | 218.6 | 14× |
| difenoconazole | 4.40 | stem | 0.034 | 0.13 | 3.9× |
| difenoconazole | 4.40 | leaf | 0.044 | 0.62 | 14× |

**log10 RMSE = 0.783** (n = 6), against the PFAS side's a-priori error of ≈ 0.84–0.95
— which still has fitted transport parameters behind it, while this has none.

**The error structure is the finding.** The **stem is predicted well** (0.5×, 1.3×,
3.9×) while the **leaf is over-predicted by 6–14×** — precisely the failure mode §3
predicts, since with `γ = 0` the leaf is an unbounded terminal accumulator over a
60-day exposure. The model also gets the ordering right across 4.5 log units: polar
compounds reach the leaf (TF ≫ 1), difenoconazole stays in the root (TF ≪ 1).

**Sensitivity** (§6 of the script) — a scan, not a fit:

| in-planta half-life | Briggs TSCF | Schriever TSCF |
|---|---:|---:|
| none (γ = 0, strict a-priori) | **0.783** | 1.212 |
| 60 d | 0.622 | 1.061 |
| 30 d | 0.488 | 0.933 |
| 14 d | 0.274 | 0.716 |
| **7 d** | **0.210** | 0.486 |
| 3 d | 0.515 | 0.331 |

Three things follow.

1. **The residual is dominated by the missing dissipation rate, not by transport
   structure** — and there is a genuine **minimum at ≈ 7 days**, not a monotone
   slide toward faster loss. So the model does not merely prefer "more metabolism";
   it predicts a specific in-planta half-life, which a measurement could test. Ge
   report *soil* half-lives only (19–41 d) but document in-plant metabolism
   qualitatively (THX → clothianidin *in the plants*).
2. **The original Briggs bell beats the modern refit**, at every half-life. The
   broader Schriever bell gives difenoconazole `TSCF = 0.37` against Briggs' 0.047,
   and Ge measured a leaf TF of 0.044 — the narrow bell is right that lipophilic
   compounds do not reach the leaf. Worth knowing, since the refit rests on the
   larger dataset and might have been assumed the better default.
3. Together with §4a the two datasets **cross-validate the diagnosis**: the
   half-life sensitivity is steep where the endpoint is an accumulator (leaf) and
   flat where it is an equilibrium (root).

### A correction the data forced

The first version of this module converted Trapp's lipid contents as *dry*-weight
fractions. They are **fresh** weight — the same basis as `W`, as `K_PW = W +
L·a·Kow^b` requires, and as Briggs' own barley anchor shows (`L·a = 10^−1.52` ⇒
`L = 2.5 %` alongside `W = 0.82`). The Liu root data caught it immediately:
log10 RMSE **0.605 → 0.198** on the partition term alone, and the Ge prediction
improved from 1.099 → 0.783. Mixing the two bases understates `K_PW` about tenfold
for lipophilic compounds.

### What these numbers do not settle

Two studies, one crop, seventeen compound-tissue pairs. The **grain compartment is
not exercised at all** — no dataset found reports a neutral Kow series in grain
under root-only exposure, the same gap that exists on the PFAS side. Liu's shoot
TFs are not used, because a 72-hour seedling exposure is not comparable to the
season-long growing run the model drives; only its root partition transfers. And
the Ge leaf result stays confounded with the missing half-life, so 0.783 is an
upper bound on the transport error rather than a measurement of it.

## 5. The other papers, and what each can and cannot do

All ten obtained papers were read; here is what each is actually good for.

| paper | verdict |
|---|---|
| **Ge 2017** `10.1016/j.envpol.2017.04.043` | **Used — §4.** Complete and self-contained: TF, log Kow and pKa all in the article. |
| **Briggs 1982** `10.1002/ps.2780130506` | **Anchor verified at source.** Eqs. 2 and 3 confirmed character-for-character; the 0.82 floor is explained as root water content, confirming the `K_PW` mapping. Its Table 1 (barley RCF/TSCF) is the QSPR's own training data, so it cannot serve as validation. |
| **Schriever 2020** `10.1016/j.scitotenv.2020.136667` | **Alternative TSCF implemented.** Reprints Briggs eq. 3 verbatim (independent verification) and supplies a 97-value refit, now selectable and benchmarked in §4. Its per-compound TSCF values are in the SI, which was not obtained. |
| **Trapp 1994** `10.1002/etc.5620130308` | **Tissue composition adopted** (root 1 %, stem/leaf 3 % lipid). Its Table 1 gives a complete driver set (transpiration, organ masses, water contents) for the bromacil runs, but the measured concentrations are in Figures 5–6 — **figure-only**, so an RMSE would need digitising. |
| **Liu 2023** `10.1016/j.scitotenv.2022.159826` | **Used — §4a**, once the SI arrived. Table S1 gives log Kow for all 21 compounds (PubChem-sourced); Tables S3 × S4 reconstruct total tissue concentrations. Its shoot TFs remain unused (72-h seedling exposure, not comparable to a season run). |
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
- Against measured rice data, with nothing fitted, it predicts **root partition**
  across 14 compounds and 5 log-Kow units at **log10 RMSE 0.281**, and **per-organ
  transfer factors** at **0.783** — both better than the PFAS side's a-priori
  0.84–0.95, which still carries fitted transport parameters.
- The error structure is interpretable and cross-validated between the two
  datasets: steep in half-life where the endpoint accumulates (leaf), flat where it
  equilibrates (root). The Ge leaf residual points to a specific in-planta
  half-life of ≈ 7 days — a testable prediction, not a calibration.
- **Air exchange is implemented** (`src/plant_air.py`, §3b), closing the largest
  structural gap in the neutral base: the leaf now has both of its non-growth
  sinks. It is opt-in (`simulate_neutral(air=True)`) because it needs `K_AW` and a
  molar mass, which the strict Kow-only a-priori run does not use — so the 0.281
  and 0.783 above are unaffected, and so is every PFAS number.
- Remaining gaps: the **grain compartment is untested** (no suitable dataset
  exists); rice-specific organ lipid contents are still borrowed from soybean;
  **tissue specific surface areas** are ratios, not measurements, which bounds how
  far an absolute volatilisation flux can be trusted; and the Brunetti `K_RW`
  discrepancy (13.3 vs a Briggs `K_PW` of ~1.0 for carbamazepine) is an open
  question against the partition core itself.
