# Neutral-organic (Briggs/Kow) path — implementation and validation status

> `src/neutral_dpu.py` · `validation/neutral_dpu_validation.py` · `tests/test_neutral_dpu.py`
> **Status: IMPLEMENTED, anchors verified at source, and validated a-priori
> against four measured datasets** — root partition on rice at log10 RMSE **0.281**
> (Liu 2023, n=14), per-organ transfer factors at **0.783** (Ge 2017), a second
> root-partition table over 11 species at **0.598** (Li 2019, n=29, §4d) and
> carbamazepine from three soils at **0.191** (Kodešová 2019, n=21, §4f) — all with
> nothing fitted — plus TSCF tested on its own for the first time (§4e).
> §4–§5 record what those numbers mean and what is left.

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

### 4c. Hwang et al. 2017 — a diagnosis, not a score (lettuce, chlorpyrifos)

`python validation/hwang2017_lettuce.py`

The only **time-resolved per-organ** neutral dataset to hand (3 samplings × 2 soil
levels), **lipophilic** (log Kow 4.01, on the falling limb of the Briggs bell where
Ge's difenoconazole sits), and the only one supplying a **measured `Kd`** — which is
what makes it usable, since it converts a soil residue into the pore-water exposure
the model needs: `Ce(t) = C0·(½)^(t/T)/Kd`.

**It is not on the footing of Liu 2023 / Ge 2017, and its RMSEs must not be quoted
as validation.** Four limits stack, two of them unresolvable from the article:
lettuce not rice and one compound; **Table 1's fresh/dry basis is not stated**
(a ~20× lever at 95 % water); the **growth curve's functional form is not in the
transcription** (`Ig`, `Kg` without the equation — reconstructed as log-log, since
an exponential reading gives 10¹⁹ g); and the roots grew **in soil**, so with
`Koc ≈ 2219` the root endpoint is contact-confounded the way Li 2025 is.

What it nonetheless establishes:

**(i) Table 1 is internally consistent on one basis.** `whole` is the mass-weighted
mean of leaf and root at a root mass fraction of **5.4 ± 0.9 %**, the same at every
sampling (4 usable rows; two are forced to 0 by 1-dp rounding). So the columns can
be read together — and 5.4 % is characteristic of **fresh** weight for lettuce
(dry weight would put ~11 % in the root). Weak evidence, but it is the only
non-circular evidence there is, and it agrees with the convention that produce
residues are reported as-eaten.

**(ii) The basis spans the verdict.** The modelled root cannot exceed its
equilibrium partition — `K_PW = 15.8 L/kg` is a *structural ceiling*:

| reading | measured root BAF | vs ceiling | consequence |
|---|---:|---:|---|
| fresh weight | 44–165 | **2.8–10.4× above** | the model structurally **cannot reach** it |
| dry weight | 4.4–16.5 | 0.28–1.04× | measurement sits **under** it; model over-predicts |

**(iii) The two readings fail on OPPOSITE organs** — the sharpest result here:

| basis | root RMSE | leaf RMSE | |
|---|---:|---:|---|
| fw | 0.711 | **0.489** | leaf predicted well, root unreachable |
| dw | **0.393** | 0.948 | root bracketed, leaf over by 3–24× |

(overall 0.610 fw / 0.726 dw). Neither reading makes the plant coherent: swapping
the basis only moves *which* organ is wrong. So the discrepancy is **not a units
artifact awaiting a footnote** — something about a lipophilic compound in a
soil-grown plant is genuinely missing.

**(iv) It corroborates an already-open problem.** On the fresh reading the measured
root exceeds the Briggs ceiling by 3–10×, and soil contact cannot account for it
(explaining the 10 mg/kg rows outright would need **12–49 %** of the washed root
mass to be soil). That is the **same direction and comparable magnitude** as the
open **Brunetti 2021** disagreement (calibrated pea root `K_RW` = 13.3 vs a Briggs
`K_PW` ≈ 1.0 for carbamazepine, §5). Two independent sightings of one thing: for
**lipophilic compounds in soil-grown plants the Briggs root partition looks too
low**.

**A trap named explicitly.** The half-life scan improves the dry reading a lot
(0.73 → 0.30 at 4 d) and makes the fresh one steadily worse — metabolism can only
lower modelled concentrations, so it helps wherever the model is too high. The fit
therefore "prefers" dry weight. **That is not evidence about the basis**: choosing
the reading that lets the model fit is circular, and it would be choosing against
(i). The basis is an open question about the *paper*, not about the model —
cheapest fix is to ask the authors.

**Not shipped as a `data_obs/` table, deliberately.** The shared `--obs` harness
runs every row on the *rice* drivers (120 d, constant exposure); Hwang is lettuce
over 40 d under a decaying one, so a number from that path would be silently
meaningless. The data lives in the dedicated script's constants instead.

### 4d. Root partition again — Li et al. 2019 (29 out-of-sample rows, 11 species)

`data_obs/neutral_obs_li2019_rcf.csv` · `validation/li2019_rcf_apriori.py`

This is the second independent a-priori test of `K_PW`, and the one that made the
partition term's error legible. It arrived by accident: Li, Chiou, Li & Schnoor
2019 (`10.1016/j.envint.2019.02.020`) was requested from the acquisition queue as
row **A1**, for a **rice root lipid** value. **Its SI does not contain one** — no
rice row appears anywhere in it. What it does contain is Table S1: **48
hydroponic root concentration factors** over 11 species and log Kow −0.57 → 5.41,
the same equilibrium endpoint as Liu 2023 but twice the size, over a wider range,
and including four rice rows.

Eighteen of the 48 are Briggs' own barley — the data the RCF QSPR was fitted to —
so they are shipped marked `subset=calibration` and held out of the score. The
remaining **29 are out-of-sample**; clotrimazole is excluded outright as a weak
base (pKa 6.02), the same speciation rule that removed the sulfonylureas from the
Liu file.

**Result: log10 RMSE 0.598 (n = 29), nothing fitted.** But the number is not the
finding. The finding is that **the error is one offset, not scatter**: every one
of the 11 species is biased the same way, low, from −0.30 (lettuce) to −0.95
(*Plantago*), with rice itself at −0.47. Eleven species do not agree by accident.

**The offset is not kinetic, and that was pre-registered rather than discovered.**
The obs file's header names the one confounder in advance: these are 6 h – 12 d
exposures, and Li et al.'s own model writes RCF = α_pt·`K_PW` with α_pt ≤ 1
falling as Kow rises, because lipophilic compounds equilibrate slowly. A model
with no α_pt should therefore *over*-predict the short, high-Kow rows. It
under-predicts them. Splitting the 29 rows settles which story is right:

| split | n | mean log10 bias |
|---|---|---|
| exposure < 1 d | 1 | −0.124 |
| exposure 1–3 d | 18 | −0.442 |
| exposure > 3 d | 10 | −0.447 |
| log Kow < 2 | 4 | **−0.030** |
| log Kow 2 – 3.5 | 6 | −0.305 |
| log Kow 3.5 – 4.5 | 11 | −0.462 |
| log Kow > 4.5 | 8 | **−0.688** |

**Flat in exposure time** — and inside the high-Kow cell alone, longer exposures
are no better (−0.586 under 3 d vs −0.517 over). So non-equilibrium is ruled out
as the driver, and the deficit sits in the **lipophilic sorption term**: it
vanishes at low Kow, where the water floor `W` dominates and there is no lipid
term to be wrong about.

**Two caveats on the Kow ladder, added after an adversarial re-read (§3c of the
script).** First, the *monotone* reading is not robust: **one study, Namiki 2015,
supplies 10 of the 29 rows and all 10 sit in the top two bins**, as two compounds
measured in five species each. Collapsing every compound × study pair to one row
gives −0.030 / −0.305 / −0.576 / **−0.501** — the top two bins go **flat**, not
still rising. What does survive is that the bias is absent below log Kow 2 and
−0.3 to −0.6 above it, and that **all ten source studies are biased the same
way** (−0.20 to −0.89), so it is not one lab's artifact. Second, and more
importantly for the decision below: an earlier version of this section argued
that a flat increase in `L` was "the wrong instrument" for a Kow-dependent
deficit. **That was wrong.** `K_PW = W + L·a·Kow^b` is dominated by the water
floor at low Kow, so scaling `L` is *inherently* Kow-dependent — restoring the
anchor shifts predictions by +0.045 log at log Kow 1 and +0.391 at log Kow 5,
which is close to the shape of the observed bias and delivers 60–75 % of it.
**So this table is genuine evidence in favour of restoring the anchor.**

**Part of that deficit is internal to the model, not a disagreement with the world.**
`neutral_dpu` anchors on Briggs' RCF, whose lipid term is `L·a = 10^−1.52 =
0.0302`. But `rice_compartments` substitutes a *measured* `L = 0.01` and keeps the
*conventional* `a = 1.22`, which makes the product `0.0122` — **2.5× below the
anchor the module claims**. The module header itself says `L` and `a` are
identifiable only as their product; the rice compartment replaces one factor and
leaves the other. The sharpest way to see the cost is on Briggs' own barley rows:

| root composition | Briggs' barley rows (n=18) | Li 2019 a-priori (n=29) | Liu 2023 rice (n=14) |
|---|---|---|---|
| **shipped** `L·a` = 0.0122, b = 0.77 | 0.266 | 0.541 | **0.206** |
| Briggs anchor `L·a` = 0.0302, b = 0.77 | **0.111** | **0.295** | 0.286 |
| Li/Chiou `f_lip` = 0.010, `K_lip` = 1.27·Kow^1.03 | 0.284 | 0.528 | 0.638 |
| *fitted (2 free params — not a-priori)* | *0.101* | *0.254* | *0.192* |

(equilibrium `K_PW` log10 RMSE; the full-ODE scan in §4 of the script agrees in
direction — Li 0.598 → 0.331 and Ge 0.783 → 0.651 at the anchor, Liu 0.281 → 0.288.)

Sizing the two contributions: restoring the anchor shifts predictions by
**+0.394 log**, against an observed **−0.688** in the worst cell (log Kow > 4.5).
So the anchor accounts for about **57 %** of it, and the remaining ~2× at high
lipophilicity does not come from the lipid fraction at all.

**The default is NOT changed, and the reasoning is the point.** Restoring the
anchor improves two tables and **degrades the only rice one**. An intermediate
`L ≈ 0.015` would fit all three — and would convert the neutral path's single real
claim, that nothing is fitted, into a fitted result. `BRIGGS_ANCHORED_LIPID_FW`
makes the alternative runnable, a test pins the 2.48× discrepancy so it cannot
drift back silently, and the promotion decision is left open.

**What A1 did settle.** The queue's own DEFINITION NOTE said the *definition* of
"lipid" mattered more than any number, and it was right. Verified at source: Li et
al.'s `f_lip` is **fresh weight** (dry-basis reports converted at 90 % root water)
and it enters an RCF expression of exactly this model's form, so their per-crop
values are commensurable with `L`. Their cereals — **barley 1.00 %, wheat
1.10–1.14 %, maize 0.53 %** — bracket the 1 % this repo runs. So the value is
corroborated and its provenance improves from a soybean model run to measured
cereal roots. Rice itself is still unmeasured.

**And what it sharpened.** Li et al. assign Briggs' barley `f_lip ≈ 1.00 %`, while
Briggs' fitted coefficient implies 2.47 %. **Briggs 1982 measures no lipid at
all** — verified at source; he attributes the 0.82 floor to root *water* and never
reports a lipid content — so both figures are inferences from a paper that
measured neither. The physical reading of the 2.5× excess is that Briggs'
coefficient carries **non-lipid sorption**: cell wall and lignin, which this
repo's PFAS side models explicitly as `f_cw·K_cw` and lists as **GAP A**, and
which the neutral composition sets to zero. If that is right, the fix is a
measured neutral-organic cell-wall coefficient, not a larger lipid fraction.

### 4e. TSCF on its own — Schriever & Lamshoeft 2020, all 97 values

`data_obs/tscf_obs_schriever2020.csv` · `validation/schriever2020_tscf.py`

Queue row **A3**, closed. Until now TSCF could only be reached *through* the plant
ODE, in §4b, where its error is inseparable from the unknown in-planta half-life —
which is exactly why the Ge result is recorded as an **upper bound** on transport
error rather than a measurement of it. The SI's Table A 3 removes the model from
the middle. All 97 rows are transcribed; the parse reproduces the SI's own Table
A 1 species tally exactly (barley 51, poplar 15, cattails 8, wheat 6, …).

TSCF is a bounded fraction with zeros present, so the metric is **linear** RMSE
plus rank correlation, not the repo's usual log10.

| QSPR | rows | RMSE | bias | Spearman |
|---|---|---|---|---|
| **`briggs_tscf`** (repo default), 30 un-ionised rows | 30 | 0.310 | **−0.221** | +0.794 |
| `schriever_tscf` — **in-sample**, these are its training data | 30 | 0.234 | −0.094 | +0.744 |
| `briggs_tscf`, barley only | 10 | 0.166 | −0.130 | +0.355 |
| `briggs_tscf`, every other species | 20 | 0.361 | −0.267 | +0.647 |

**The default TSCF under-predicts, by 0.22 on a 0–1 scale** — the *same direction*
as the partition finding in §4d, which is why the two are reported together: both
a-priori inputs of the neutral path are biased low, and this session measured each
one directly. Two qualifications, both load-bearing. The Briggs bell is
out-of-sample here (these are Briggs 1987 / Inoue / Shone & Wood compounds, not
the 1982 set) but **not independent** — same lab, same species, same method
lineage. And `schriever_tscf`'s 0.234 is *reproduction*: it is the floor a fitted
model reaches on its own training set, not a rival result.

One free by-product. On the full 97, Schriever & Lamshoeft's central claim is
confirmed on their own table: switching the descriptor from log P to log D at the
test pH lifts the Briggs bell's rank correlation from **+0.313 to +0.653**. That
does not change this repo — the neutral path is un-ionised by construction — but
it is an independent check that the transcription carries real signal.

### 4f. Carbamazepine from three soils — Kodešová et al. 2019 (21 rows), and the anchor vote

`data_obs/neutral_obs_kodesova2019.csv` · `validation/kodesova2019_carbamazepine.py`

Queue row **A4**, closed by the SI. This table is small but it was fetched to
settle two specific things, and it settles both — in the direction that was not
expected.

**The exposure is the cleanest in the repo.** Three measured quantities, no model
between them: root concentration (SI Table S2, ng/g **dry weight**, basis stated),
soil concentration from the *same pot at the same harvest* (Table S4), and the
paper's own measured Freundlich isotherms per soil (Table 1). Then
`c = (C_soil/K_F)^n` and `RCF = C_root/c`. No mass balance, no pot geometry, no
dissipation model — and carbamazepine is un-ionised across the whole
environmental range (pKa 1.0 / 13.9) with `DT50 > 1000 d` in all three soils, so
the authors' own reason for declining to compute bioaccumulation factors
(repeated dosing plus dissipation) does not bite for this one compound.

The isotherm convention is the single pivotal assumption, so it is argued rather
than assumed. At `c = 1 mg/L` the isotherm gives `S = K_F`, so `K_F/f_oc` is a
`K_oc`: **222 / 189 / 154** across three soils spanning 3.8× in organic carbon —
consistent, and inside carbamazepine's literature band. The alternative unit
reading gives ~1600 everywhere. The derived pore water (0.10–0.70 mg/L) also
lands on the scale of both the applied solution (~0.5–1.0 mg/L) and the sorption
study's calibration range (0.5–10 mg/L).

**Result: log10 RMSE 0.191 (n = 21), nothing fitted.** For scale, on the same
path: Liu 2023 0.281, Li 2019 0.598. **But see §4g** — that 0.191 is flattered by
a scoring artifact; on the appropriate equilibrium basis it is **0.237**, and Liu
2023 at 0.206 is the best a-priori result here.

**It votes against restoring the Briggs anchor**, and it is well-conditioned to
do so: at log Kow 2.25 the two compositions differ by 1.6×.

| root composition | Kodešová n=21 | Liu 2023 rice n=14 | Li 2019 n=29 |
|---|---|---|---|
| **shipped** `L` = 0.010 | **0.191** | **0.281** | 0.598 |
| anchored `L` = 0.0247 | 0.216 | 0.288 | **0.331** |

Driver-free cross-check, in case the rice ODE drivers are doing the work:
measured `RCF_fw` median **1.10** against `K_PW` of **1.56** shipped and **2.53**
anchored — the shipped composition already runs *high* here, and the anchor
doubles the excess. The absolute numbers differ from the table above (the ODE's
root BAF sits below `K_PW` because the xylem drains it, and log Kow 2.25 is near
the TSCF peak where that drain is largest) but **the ordering is identical under
both framings**, which is what the decision rests on. It also survives the
dw→fw lever: across root water contents 0.85–0.95 the anchored composition stays
0.21 log further from the data, because the lever moves both together.

**And it weakens the Brunetti sighting — the one it was fetched for.**

| | carbamazepine root partition |
|---|---|
| Brunetti 2021, *calibrated* green-pea `K_RW` | 13.3 L/kg fw |
| **Kodešová 2019, measured**, 4 plants × 3 soils | **1.10** (0.47–2.49) |
| Briggs `K_PW` as this model runs it | 1.56 |

Brunetti's value is **~12× above a direct measurement of the same compound**,
while Briggs lands within ~1.5×. So that disagreement is far more likely to sit
in Brunetti's calibration — a posterior free to absorb whatever else the model
was missing — than in the partition core. **This table is a counter-example to
the four-sightings group below**, not a fifth member of it.

Honest limits: four leafy/root vegetables and no rice; one compound at one
lipophilicity, so it says nothing about the high-Kow end where §4d's deficit
lives; roots were rinsed, not exhaustively cleaned, though adhering soil would
bias the observation *up* and so works against the conclusion drawn here.

### 4h. The soil half of Li 2019 — 376 rows that flip the sign

`data_obs/neutral_obs_li2019_soil.csv` · `validation/li2019_soil_table.py`

Table S2 of the same paper as §4d: **376 soil-grown root RCFs over 13 crops**,
twelve times the hydroponic table, with each crop's measured root lipid attached.
It was overlooked when the hydroponic table was mined, and it changes the reading.

| Li 2019 table | n | mean log10 bias |
|---|---|---|
| **S1, hydroponic** (§4d) | 29 | **−0.432** — model LOW |
| **S2, soil** | 376 | **+0.260** — model HIGH |

Same authors, same fresh-weight convention, same endpoint definition, opposite
sign. So "the root partition is too low" is a property of **one exposure route**,
not of the partition core — and the larger half argues *against* moving the
default:

| root composition | bias | RMSE |
|---|---|---|
| **shipped** `L` = 0.010 | +0.260 | **0.639** |
| anchored `L` = 0.0247 | +0.645 | 0.873 |

**Two results here do generalise.**

*The composition term works.* Substituting each crop's **own** measured lipid for
the model's rice value — across an 11× spread, radish 0.09 % to wheat 1.14 % —
takes the mean bias from **+0.260 to −0.001** and tightens RMSE 0.639 → 0.461. So
the *form* `K_PW = W + L·a·Kow^b` handles a large range in `L` correctly. (Caveat:
`f_lip` is Li's own model input, so this is consistency inside their framework,
not independent validation.)

*Most of the apparent error is the exposure, not the plant.* `value` here is
**derived** — Li divide the measured soil concentration by `K_om` — and their
Table S3 records which `K_om` was measured and which was estimated:

| `K_om` source | n | bias |
|---|---|---|
| **experimental** | 62 | **+0.033** |
| from a QSPR | 261 | +0.291 |
| unmatched | 53 | +0.371 |

**On the rows with a measured `K_om` the model is essentially unbiased.** The soil
organic-matter gradient says the same (bias +0.499 below 1 % OM, +0.101 above 4 %).
That now holds on three independent tables — these 62 rows, Liu 2023, and
Kodešová 2019, all with a measured or directly known exposure — which relocates
most of the neutral path's apparent partition error onto the **exposure term**,
a soil-side problem rather than a `K_PW` problem.

One framing note. Li's own model writes `RCF_water = α_pt·[f_pw + f_ch K_ch +
f_lip K_lip]` with a median `α_pt` of **0.098** here, i.e. their equilibrium term
is ~10× their measured RCF. This model has no `α_pt` and lands on the measured
values anyway, because Briggs' `b = 0.77` is far flatter than Li's
`K_lip ≈ Kow^1.03` (20× apart at log Kow 5). What Li put in `α_pt`, Briggs put in
the exponent; dividing by `α_pt` here would double-count.

### 4i. Kodešová's leaf half — translocation, and metabolism that was measured

Sections 6 and 7 of `validation/kodesova2019_carbamazepine.py`. The root
comparison in §4f used half the table; the leaf half answers two different
questions, and one of them removes a free parameter.

**Translocation, with the exposure cancelled.** A leaf/root ratio needs no
exposure at all, so it does not inherit the isotherm assumption §4f defends.
Measured (n=21, 20–26 d): median **3.25**, range 0.31–9.05. Model at the matched
harvest: **181** — ~55× over. That is §3's terminal-leaf runaway, measured for
the first time against per-organ data with the exposure divided out.

**But metabolism cannot close it**, and that is the useful part:

| in-planta half-life | ∞ | 30 d | 14 d | 7 d | 3 d | 1.5 d |
|---|---|---|---|---|---|---|
| model leaf/root | 181 | 147 | 119 | 86 | 49 | **28.5** |

Even a 1.5-day half-life leaves the model ~9× above the measurement. So the leaf
excess is **not only the missing γ** — most of it is the driver mismatch, a rice
season's transpiration per unit leaf mass applied to a 340 cm³ pot of lettuce.
**That is a warning about §4b**: the Ge 2017 half-life minimum at ≈7 d may be
absorbing the same mismatch rather than measuring metabolism, so it should not be
quoted as a prediction of an in-planta half-life without that caveat.

**Metabolism, measured rather than fitted.** Kodešová quantified carbamazepine's
four metabolites alongside the parent, so the parent fraction is a direct
observation — the first in this repo:

| | root | leaf |
|---|---|---|
| lamb's lettuce | 0.874 ± 0.030 | **0.169** ± 0.044 |
| radish | 0.963 ± 0.005 | **0.810** ± 0.047 |
| both | 0.919 | 0.489 |

Two consequences. **Transformation is a shoot process**, and it is large — so
`γ = 0` is simply wrong for carbamazepine even though its *soil* DT50 exceeds
1000 d. Soil persistence does not imply plant persistence. And it is strongly
**species-dependent** (0.17 vs 0.81 on the same compound, same soils, same
harvest), so an in-planta half-life is not a compound property — which is exactly
what fitting one to a single dataset assumes. Note also that the metabolites are
not inert: in lamb's lettuce leaves the epoxide *exceeds* the parent (17000 vs
6400 ng/g), so a parent-only model understates the total burden several-fold.

### 4j. Briggs 1983 — the STEM had an anchor too, and it was never read

`validation/briggs1983_stem.py` · `neutral_dpu.briggs_scf`

The companion to the 1982 paper this whole path is built on was sitting unread in
the obtained set. It does for the shoot what 1982 did for the root — measures two
compound series in barley shoots after root uptake and fits a partition of exactly
the `K_PW` form (**verified at source**):

```
log(K_stem/xylem_sap − 0.82) = 0.95·log Kow − 2.05      (eq. 2)
SCF = K_stem/xylem_sap × TSCF                            (eq. 3)
```

The implementation checks itself: computed from the coefficients alone, eq. 3
peaks at **6.4 at log Kow 4.5**, against the paper's own "about 6 … at about
log Kow = 4.5".

**The stem had no anchor in this repo until now** — it runs Trapp 1994's soybean
3 % lipid with the conventional `a = 1.22` and Briggs' *root* exponent 0.77:

| organ | shipped `L·a` | Briggs anchor | ratio | exponent |
|---|---|---|---|---|
| root | 0.0122 | 0.0302 | **0.40** (2.5× below) | 0.77 vs 0.77 — same |
| stem | 0.0366 | 0.0089 | **4.11** (4.1× above) | 0.77 vs **0.95** |

So the two organs were parameterised from unrelated sources and neither was
checked against the anchor that existed for it — the root missing low, the stem
missing high.

**But the consequences differ sharply, and that is the part to quote.** For the
root the coefficient gap *is* the disagreement, because the exponent is shared.
For the stem it largely **cancels** against Briggs' steeper slope:

| log Kow | 0.0 | 1.0 | 1.78 | 2.5 | 3.5 | 4.5 | 5.5 |
|---|---|---|---|---|---|---|---|
| repo SCF ÷ Briggs SCF, log10 | +0.02 | +0.07 | +0.13 | +0.13 | −0.02 | −0.20 | −0.38 |

At most **0.13 log** across log Kow 0–3.5, the range where the TSCF bell delivers
anything at all; the larger deviations sit above 4.5 where TSCF has collapsed and
the stem receives almost nothing — and where Briggs' own text says the predicted
decline "was not tested", so both sides are extrapolating.

**Nothing is changed.** The anchor is one species and is explicitly a shoot *base*
rather than a true stem, and unlike the root case no measured table in this repo
would arbitrate it. It is recorded, exposed as `briggs_scf()`, and pinned by a
test. The stem is a **provenance** problem, not a prediction problem, and it ranks
below both the root question and the exposure-term work.

*(The last sentence was an argument when it was written. §4k now measures it —
and the "no measured table would arbitrate it" clause is what turned out to be
wrong: the arbitrating table was inside the same paper, in a section §4j did not
read.)*

### 4k. Briggs 1983's own table — the stem, scored against measurements

`data_obs/neutral_obs_briggs1983_shoot.csv` · `validation/briggs1983_shoot.py` ·
`tests/test_briggs1983_shoot.py`

§4j compared the repo's stem to Briggs' *fitted equations*. His **Table 1** — the
data behind them — had never been transcribed. It gives, for **16 non-ionised
chemicals** over log Kow −0.57 to 3.7, the distribution of the shoot burden
across leaf blade / central stem / stem base at **24 and 48 h**, after uptake by
11-day-old barley from a **nutrient solution of known concentration**. Section 2.2
gives the section fresh weights (0.54 / 0.40 / 1.2 g, six plants). Those four
published numbers reconstruct the Stem Concentration Factor the article defines:

```
SCF = (total dpm × section %/100) ÷ section g ÷ solution dpm mL⁻¹     [L/kg fw]
```

**The reconstruction is checked before anything is scored with it**, three ways.
Eq. (3) computed from the transcribed coefficients peaks at **SCF 6.39 at log Kow
4.43** against the paper's stated "about 6 … at about 4.5". The reconstructed
central-stem SCF sits on that curve at **bias +0.049** — the article says its own
points "fit quite well" there, so the derived numbers *are* its Fig. 4. And the
stem **base** runs high (+0.188), which the article predicts and explains as
direct contact with the treating solution. Every (compound × harvest) row sums to
100 ± 0.1 %, which caught two digits the PDF's text layer had wrong.

**The result — the repo's first measured stem test, nothing fitted:**

| | n | bias | RMSE |
|---|---|---|---|
| **repo `K_PW` × TSCF**, central stem | 30 | **−0.030** | **0.299** |
| Briggs' own *fitted* eq. (3), same rows | 30 | +0.049 | 0.282 |
| repo, stem base (contact-confounded) | 30 | +0.109 | 0.290 |

The repo's stem predicts measured SCF **as well as the paper's own equation**,
which had this very data to fit — and at the level the root tables reach (Liu
0.281, Li 2019 hydroponic 0.598). **This settles §4j empirically**: the 4.1×
coefficient gap really does cancel in the observable, so it is a provenance
problem and not a prediction problem. The rice culm's composition is still
unmeasured; it is simply not *wrong* anywhere it can be checked.

**What it does NOT do, which is why item 4 was worth running rather than
assuming.** The table was wanted to constrain the **stem/leaf split**
independently of Ge 2017. It cannot. The stem equilibrates and the leaf does
not, so the split is a function of **exposure duration**: the stem's share of the
shoot burden falls by a **median 1.55× (up to 2.11×) between two harvests one day
apart**. Scored as an equilibrium the leaf drifts by **+0.229 at 24 h and +0.586
at 48 h** — the bias *growing* with exposure is the terminal-accumulator
signature §3 describes, and the article states it independently ("leaf amounts
generally increased up to 72 or 96 h"). So the leaf half is a **second,
independent confirmation of the terminal-accumulator structure** — different
species, different exposure from Ge 2017 — and not a split constraint.

Three compounds deviate and the reconstruction reproduces the article's own
diagnoses for each: **aldicarb** 3.6× high (it puts it at "about three times",
from in-planta oxidation to a trapped polar sulphoxide), **aldoxycarb** 3.7× high
(log Kow < 0, where the article says its own TSCF curve runs low), and
**4-(4-bromophenoxy)phenylurea** 6.5× low (never equilibrated — the article
excludes it from its own Fig. 4, and so does the CSV).

**Scope.** Barley, 11-day seedlings, 24–48 h, one laboratory, two chemical
series, and section weights that are the article's *typical* values applied to
every test. It anchors the **partition form** for a stem. It says nothing about
a rice culm. Deliberately **not** wired into the `--obs` harness, for the Hwang
2017 reason (§4c): these are sections of a 48-hour barley seedling and the shared
harness would run them against a 120-day rice season.

### 4l. The weak-electrolyte path, tested — direction supported, magnitude refuted

`validation/weak_electrolyte_tscf.py` · `tests/test_weak_electrolyte_tscf.py`

The speciation port (`simulate_neutral(pKa=…)`) landed labelled *"structural
capability, not a predictive claim — no measured weak-electrolyte **rice**
dataset exists here"*. True of rice, and it was read as though no measured
weak-electrolyte data existed at all. It did: §4e scores only the **30** rows of
Table A 3 flagged un-ionised at the test pH and explicitly sets the other **67**
aside as *"outside this model's stated scope"*. The port is the thing that
extends the scope to them, so **the held-out rows were in the repo before the
capability was**, and they are its first empirical test.

**How `f_n` is obtained, and why not from the pKa.** From the table's own logD:
`f_n = 10^(logD_test − logP)`, clipped at 1 — the definition of the distribution
coefficient when only the neutral species partitions. This needs no pKa and no
acid/base label, which matters because the table has neither reliable. Its pKa
column and its logD column **disagree about how ionised these compounds are**:
Henderson–Hasselbalch on the pKa misses the reported logD by a median **1.41 log**
(acid reading; 1.82 as a base), and 20 rows have logD *above* logP, which no
single-centre acid or base can produce. The logD route is also what the shipped
`neutral_at_test_pH` flag is built from, so script and flag stay consistent.

> **A false counterexample, named so it is not rediscovered.** The eight barley
> rows at **pKa 1.62** have TSCF 0.63–0.98. Read as a strong acid they are fully
> dissociated at any test pH and look like a spectacular refutation. They are
> not: their own logD sits 0.01–0.02 *above* logP, so the table says they are
> un-ionised (the pKa is a basic centre) and the flag already classes them
> neutral. Deriving `f_n` from the pKa column manufactures this artefact.
> Pinned by `test_the_pKa_162_rows_are_NOT_ionisable_the_false_counterexample`.

**What the port predicts.** Read straight off `root_uptake` with `Vmax = 0` and
`P_n = kappa_d`, the membrane influx conductance relative to the same molecule
un-ionised is `Φ = f_n + g(N)·f_d / 10^3.5` — nothing fitted, just the shipped
Trapp ratio and the GHK factor. Over this table's `f_n` range (4.7e−5 → 1, 4.3
orders of magnitude) **Φ moves ~1.6e4-fold**.

**What the data say.** Measured TSCF over the same rows moves about **3-fold**.

| | n | RMSE | bias | Spearman |
|---|---|---|---|---|
| speciation **OFF** (un-ionised) | 67 | **0.272** | +0.023 | +0.284 |
| weak electrolyte, read as an **acid** | 67 | 0.304 | −0.203 | **+0.520** |
| weak electrolyte, read as a **base** | 67 | 0.303 | −0.202 | +0.516 |

Measured transfer really does rise with the neutral fraction —
**Spearman(f_n, TSCF) = +0.480** over 67 rows — and that is the port's first
empirical support of any kind; the sign could have come out flat or backwards.
But the size fails: at `f_n < 1e−3` (13 rows) the model predicts effectively no
transfer and the measured mean is **0.127**.

**The two metrics disagree, and they are not equally solid.** Bootstrap
(n = 4000): the **rank gain holds in 94% of resamples**, the RMSE loss in only
**82%** — a small subsample flips the RMSE, and one did, which is why the guard
test asserts the ordering and the under-delivery but *not* the RMSE. So the
load-bearing claim is the **ordering** gain; *"speciation makes the fit worse"*
is a tendency, not a finding. Rank is also the half the plant model consumes
(§4e's own argument), so this is a real gain wrapped around an unusable
calibration, not a flat failure.

**Why, and it is already known one directory away.** The model's only route into
the plant is across the root membrane, so an ion that cannot cross cannot arrive.
Real ions reach the xylem **apoplastically** — around the cells, not through them
— and on the PFAS side of this same codebase the passive GHK route had to be
supplemented by a **fitted carrier** for exactly this reason (CLAUDE.md §2: net
uptake requires the carrier to overcome electrostatic exclusion). The weak
electrolyte gets passive ionic permeability and no carrier, so it inherits that
known deficit with no lever to absorb it.

**Status.** `pKa=` stays opt-in and nothing about the model changed here. It
moves from **unvalidated** to **bounded**: usable for the *direction* of a
speciation effect, not its size, and not at all below `f_n ≈ 0.1`.

**Limits.** 16 species, none rice; no compound names in Table A 3, so no row can
be cross-checked or excluded on chemical grounds; acid/base inferred and
therefore reported both ways; and TSCF is not the per-organ endpoint the model
targets.

### 4g. A scoring artifact that affects every root number above

`compare_to_obs(..., mode="equilibrium")` · `tests/test_li2019_schriever_tables.py`

All three root-partition tables (§4a, §4d, §4f) measure an **equilibrium** over
24 h – 26 d. Until now they were all scored by running the **120-day rice
season** and reading `baf_final["root"]`. That imposes a discount which is
**Kow-dependent and purely model-side**: the xylem drains the root hardest near
the TSCF peak.

| log Kow | −0.5 | 1.0 | **1.78** | **2.25** | 3.72 | 5.0 |
|---|---|---|---|---|---|---|
| ODE root BAF ÷ `K_PW` | 0.91 | 0.61 | **0.55** | **0.57** | 0.85 | 0.99 |
| in log units | −0.04 | −0.22 | **−0.26** | **−0.24** | −0.07 | −0.01 |

Applying a rice season's drain to a 24-hour barley measurement is an error, not a
modelling choice. Rescored on the appropriate basis:

| table | ODE basis (as published) | **equilibrium basis** |
|---|---|---|
| Liu 2023, rice root, 72 h | 0.281 | **0.206** |
| Li 2019, root, 6 h – 12 d | 0.598 | **0.541** |
| Kodešová 2019, root, 20–26 d | 0.191 | **0.237** |
| Ge 2017, per-organ TF at 60 d | 0.783 | 0.783 — *n/a, no root rows; the season IS the endpoint here* |

**Two things follow.** First, **Kodešová was flattered by the artifact** — it sits
at log Kow 2.25, almost exactly where the discount peaks — so the claim that it
is "the best a-priori result in the repo" is wrong on the appropriate basis;
**Liu 2023 at 0.206 is**. Second, the anchor verdicts of §4f **survive and their
margins widen**:

| root composition | Liu (rice) | Li 2019 | Kodešová |
|---|---|---|---|
| **shipped** `L` = 0.010 | **0.206** | 0.541 | **0.237** |
| anchored `L` = 0.0247 | 0.286 | **0.295** | 0.410 |

So removing the artifact **sharpens the contradiction rather than resolving it**,
which is the honest outcome: the tables disagree more clearly, not less.

**The default stays `mode="ode"`** so no published number moves silently, and the
mode is reported rather than switched. Which basis should headline the neutral
path is a one-line decision left open with the anchor question.

**A consequence worth stating plainly: this document quotes both bases, in
different sections.** §4a, §4d and §4f are ODE-basis; §4h's soil numbers are
equilibrium-basis, because `li2019_soil_table.py` scores `k_pw` directly with no
ODE. Neither is wrong, but an RMSE here means nothing without its mode — the same
table reads **0.549 → 0.670** (ODE) or **0.639 → 0.873** (equilibrium) across the
two lipid readings. `--lipid-source both --mode …` now prints either, so the
labels can be checked rather than trusted.

### What these numbers do not settle

Two studies, one crop, seventeen compound-tissue pairs. The **grain compartment is
not exercised at all** — no dataset found reports a neutral Kow series in grain
under root-only exposure, the same gap that exists on the PFAS side. Liu's shoot
TFs are not used, because a 72-hour seedling exposure is not comparable to the
season-long growing run the model drives; only its root partition transfers. And
the Ge leaf result stays confounded with the missing half-life, so 0.783 is an
upper bound on the transport error rather than a measurement of it.

## 5. The other papers, and what each can and cannot do

Here is what each obtained paper is actually good for. The second block is the
`DPU4OC_add.zip` delivery; `docs/literature_db/Acquisition_Queue.csv` carries the
per-row status and what is still missing.

| paper | verdict |
|---|---|
| **Li 2019** `10.1016/j.envint.2019.02.020` | **Used — §4d**, but not for what it was requested for. Asked for a rice root lipid (queue A1); **contains no rice**. Delivered instead the lipid's operational definition verified at source (fresh weight, 90 % root water) plus Table S1's 48 hydroponic RCF values, now the repo's second a-priori partition test. |
| **Schriever 2020 SI** `10.1016/j.scitotenv.2020.136667` | **Used — §4e.** Queue A3 closed. All 97 TSCF values transcribed, letting the TSCF QSPR be tested with no plant model in between for the first time. |
| **Rosado 2022** `10.3389/fpls.2022.868319` | **Read, deliberately not adopted.** Queue A2. Confirms rice straw lipophilic extractives at **3.4 ± 0.1 % dry basis** (dichloromethane), stated explicitly. Not interchangeable with the model's fresh-weight octanol-equivalent `L`, and converting it needs a straw fresh/dry ratio the paper does not give. |
| **Kodešová 2019** `10.1007/s11356-019-04333-9` | **Used — §4f.** Queue A4, closed: the SI arrived. Carbamazepine root partition, 4 plants × 3 soils, with the exposure derived from the paper's own measured soil concentrations and Freundlich isotherms. **log10 RMSE 0.191 on the ODE basis, 0.237 on the equilibrium basis** — and it does two things to the open questions: it **supersedes the Brunetti sighting** (measures 1.10 where Brunetti calibrated 13.3 for the same compound) and it **votes against restoring the Briggs anchor**. It is *not* the repo's best a-priori result: §4g shows the 0.191 was flattered by the ODE scoring artifact, which peaks almost exactly at this compound's log Kow 2.25. Liu 2023 at **0.206** is. |
| **"McFarlane 1987"** | **Wrong paper.** The supplied file is a constructed-wetland nitrogen isotope study — no bromacil, no plant uptake. Queue A5 remains open. |
| **Deng 2018 / Shi 2025 / Li 2021 / Kondo 2019 / Boulange 2015 / Phong 2008** | **Screened out for queue C1** (a neutral compound in rice **grain** under root-only exposure). Deng samples exactly the right matrices — soil, straw, hull, brown rice — but reports **all final residues below the detection limit**. The rest are subcellular-uptake or water/soil dissipation studies with no grain compartment. The grain remains the one entirely untested compartment. |
| **Gu 2017 / Parida 2021 / Wang 2024 / Li 2019 (PeerJ)** | **Screened out for queue C2** (tissue specific surface area). All are **root** morphology, and the air-exchange term takes no root contribution at all; root area enters this model only through `a_R`, which is lumped into `kappa_d` and is documented as non-separable without inhibitor data. Parida's specific root area is figure-only and dry-basis. |
| **Honda 2023 / Ji 2015** | **Screened out for queue C3** (rice root total lipid). Honda quantifies 120 lipid *species* compositionally (mol %, −P/+P ratios) with no total mass fraction; Ji measures bound suberin diacids, one apoplastic class. Neither gives an `L`. |
| **Ge 2017** `10.1016/j.envpol.2017.04.043` | **Used — §4.** Complete and self-contained: TF, log Kow and pKa all in the article. |
| **Briggs 1982** `10.1002/ps.2780130506` | **Anchor verified at source.** Eqs. 2 and 3 confirmed character-for-character; the 0.82 floor is explained as root water content, confirming the `K_PW` mapping. Its Table 1 (barley RCF/TSCF) is the QSPR's own training data, so it cannot serve as validation. |
| **Schriever 2020** `10.1016/j.scitotenv.2020.136667` | **Alternative TSCF implemented.** Reprints Briggs eq. 3 verbatim (independent verification) and supplies a 97-value refit, now selectable and benchmarked in §4. Its per-compound TSCF values are in the SI, which was not obtained. |
| **Trapp 1994** `10.1002/etc.5620130308` | **Tissue composition adopted** (root 1 %, stem/leaf 3 % lipid). Its Table 1 gives a complete driver set (transpiration, organ masses, water contents) for the bromacil runs, but the measured concentrations are in Figures 5–6 — **figure-only**, so an RMSE would need digitising. |
| **Liu 2023** `10.1016/j.scitotenv.2022.159826` | **Used — §4a**, once the SI arrived. Table S1 gives log Kow for all 21 compounds (PubChem-sourced); Tables S3 × S4 reconstruct total tissue concentrations. Its shoot TFs remain unused (72-h seedling exposure, not comparable to a season run). |
| **Inao 2018a/b** `10.1584/jpestics.D17-083` / `D17-084` | **Cannot test the 4-compartment split.** As flagged before the papers arrived, they sample *"the whole shoot of the rice plant above the water surface"* — no organ resolution. Still the only source of a measured paddy-water + layered-soil `C_w^o(t)`, so it remains the right dataset for a future HYDRUS-coupling test against a lumped shoot. |
| **Brunetti 2021** `10.1021/acs.est.0c07420` | **Largely superseded — see §4f.** Its Table 1 reports *calibrated* posteriors, not raw concentrations: green-pea root `K_RW = 13.3` cm³/g fw for carbamazepine. That looked like the strongest evidence against the partition core until the A4 SI supplied a **direct measurement of the same compound**: 1.10 L/kg fw across 4 plants and 3 soils, against Briggs' 1.56. Brunetti's value is ~12× above the measurement, so the disagreement most likely sits in the calibration rather than in Briggs. |
| **Hwang 2017** `10.1371/journal.pone.0172254` | **Run — §4c.** Its measured `Kd` converts the soil residue to a pore-water exposure, which is what makes it usable at all. Outcome is a *diagnosis, not a score*: the unstated fresh/dry basis spans the verdict, and the two readings fail on **opposite organs**. |
| **Briggs 1983** `10.1002/ps.2780140506` | **Fully mined — §4j (equations) and §4k (Table 1).** Shoot-distribution companion to the 1982 paper. Its Table 1 is now `data_obs/neutral_obs_briggs1983_shoot.csv`: 16 chemicals × 2 harvests × 3 shoot sections, hydroponic with a known solution, reconstructed into the SCF the article itself plots. It supplies the repo's **only measured stem test** (a-priori RMSE **0.299**) and, unexpectedly, a second independent confirmation of the terminal-accumulator leaf. It does **not** supply a stem/leaf split — that quantity moves ~1.6× in 24 h. |

### Is the root partition too low? No — the evidence splits by EXPOSURE ROUTE

This looked, in turn, like one accumulating result, then like a lipophilicity
split, and then like a contradiction between two hydroponic datasets. With Li
2019's soil table (§4h) added it resolves into something simpler and more useful.

| dataset | n | exposure | mean bias |
|---|---|---|---|
| Liu 2023 (rice) | 14 | hydroponic, **known solution** | −0.053 |
| Li 2019 S2, `K_om` **experimental** | 62 | soil, **measured** conversion | **+0.033** |
| Kodešová 2019 | 21 | soil, **measured** isotherm | +0.162 |
| Li 2019 S2, `K_om` from a QSPR | 261 | soil, **estimated** conversion | +0.291 |
| Li 2019 S1 | 29 | hydroponic | **−0.432** |

**Where the exposure is measured or directly known, the model is close to
unbiased** — three independent tables, two exposure routes, bias −0.05 to +0.16.
The large deviations sit where the exposure was *estimated* (+0.29) or, in the S1
case, where it is known but the measurements themselves are the outlier.

That reframes the open question twice over. It is not "is `K_PW` too low" — on the
best-conditioned data it is about right. And most of the neutral path's apparent
partition error is an **exposure-term** problem, which is the soil side, not the
plant side.

**The S1 anomaly stays open.** It is the one hydroponic table with a known
solution that disagrees, by −0.43, and no subgroup explains it: aquatic −0.39 vs
terrestrial −0.46, organochlorines −0.51 vs everything else −0.35, the four rice
rows −0.28, and all ten source studies negative. A plausible mechanism is
root-surface sorption inflating short-exposure hydroponic RCFs, which would bite
hardest for the hydrophobic compounds — consistent with the −0.51 vs −0.35 split
— but that is a hypothesis, not a finding. On propiconazole, the one compound
S1 and Liu 2023 both measured at log Kow 3.72, they report RCF **43.65**
(lettuce) and **9.32** (rice): **4.7× apart**, against an anchor worth 2.4×.

**For the anchor decision this is now three tables against and one for**, and the
one for is the hydroponic half of a paper whose own soil half says the opposite.

A sorbing phase the neutral composition lacks remains a live physical idea, and
the PFAS side carries exactly one (`f_cw·K_cw`, **GAP A**) — a measured
neutral-organic cell-wall coefficient would serve both paths. But it is no longer
needed to explain any of the tables above, so it drops below the exposure-side
work in priority.

### The decision, and why it is a mode rather than a value

**Settled: the default stays `lipid_source="measured"` (`L_root` = 1 % fw).** The
reason is *"the evidence no longer supports moving"*, **not** *"1 % is validated
for rice"* — it is not. No rice organ lipid has been measured; the root value is
corroborated by measured *cereal* roots (Li 2019: barley 1.00, wheat 1.10–1.14,
maize 0.53 %) and stem/leaf are still Trapp's soybean figures.

The alternative is kept **runnable, not deleted**, because the question is open
rather than answered:

```bash
python validation/neutral_dpu_validation.py --lipid-source both   # add --mode equilibrium
```

| table | n | ODE: measured / anchor | better | equilibrium: measured / anchor |
|---|---|---|---|---|
| Liu 2023 (rice) | 14 | **0.281** / 0.288 | measured | **0.206** / 0.286 |
| Ge 2017 (per-organ) | 6 | 0.783 / **0.651** | anchor | 0.783 / **0.651** |
| Li 2019 hydroponic | 29 | 0.598 / **0.331** | anchor | 0.541 / **0.295** |
| Kodešová 2019 | 21 | **0.191** / 0.216 | measured | **0.237** / 0.410 |
| Li 2019 soil | 376 | **0.549** / 0.670 | measured | **0.639** / 0.873 |

Ge has no root rows, so the basis does not apply to it — the two columns agree.

Three points the raw tally hides. The anchor **raises** predicted root uptake, and
the tables it damages are the ones where the model **already runs high** — so the
two directions are not symmetric. The largest table it damages (n=376) is the
**soil half of the same paper** as the one it most improves. And only the product
`L·a` is identifiable, so "raise `L` to 0.0247" and "raise `a` to 3.02" are the
same model — which means *"don't fit the measured `L`"* is **not** an argument
against the anchor. Note `a = 1.22` has no citation anywhere in this repo.

What would settle it is named in §5: a hydroponic **rice** root RCF above
log Kow 3.5.

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
- On a **second, larger** root-partition table (Li 2019, 29 out-of-sample rows,
  11 species, §4d) it reaches **0.598** — worse, and *informatively* so: the error
  is a single offset shared by every species, traced to the model's own
  composition running 2.48× below the Briggs anchor it claims. Read the two
  numbers together: 0.281 is the result on rice over a narrower Kow range, 0.598
  is what happens when the range and the species count go up.
- On the **cleanest exposure in the repo** (Kodešová 2019, §4f — measured soil
  concentration plus the paper's own measured isotherm, same pot, same harvest)
  it reaches **0.191**. That table **reverses two open questions**: it supersedes
  the Brunetti disagreement and it votes against restoring the Briggs anchor,
  opposite to Li 2019. It is **not** the best a-priori result here, though it was
  briefly recorded as such: §4g shows the 0.191 was flattered by the ODE scoring
  artifact, whose Kow-dependent discount peaks almost exactly at carbamazepine's
  log Kow 2.25. On the basis these tables actually measure it is **0.237**, and
  **Liu 2023 at 0.206** is the best.
- **The stem now has a measured test too** (§4k), which it did not when §4j was
  written: against Briggs 1983's own Table 1 — 16 chemicals, hydroponic, exposure
  known — the repo predicts the measured Stem Concentration Factor at **log10
  RMSE 0.299 with nothing fitted**, statistically indistinguishable from the
  paper's own *fitted* equation (0.282). So the 4.1× stem coefficient gap §4j
  found is confirmed to be a **provenance** problem, not a prediction one. The
  same table also confirms the **terminal-accumulator leaf** from a second,
  independent direction (its equilibrium-scored bias grows +0.229 → +0.586
  between the 24 h and 48 h harvests), but it does **not** yield a stem/leaf
  split — that quantity moves a median 1.55× in one day, so it is a property of
  the exposure duration rather than of the compound.
- **Every number here is on `lipid_source="measured"`**, the default — and on
  `mode="ode"` *except* §4h's soil numbers, which are equilibrium-basis. Both are
  named modes rather than buried constants precisely because neither choice is
  settled: `neutral_dpu.LIPID_SOURCES` selects the root-lipid reading (1 % fw
  measured vs Briggs' implied 2.47 %, up to 2.5× apart for lipophilic compounds)
  and `compare_to_obs(mode=…)` selects the scoring basis. **An RMSE in this
  document means nothing without both labels** — the same soil table reads
  0.549 → 0.670 or 0.639 → 0.873 depending only on the basis. Run
  `validation/neutral_dpu_validation.py --lipid-source both [--mode equilibrium]`
  to check any of it rather than taking this document's word for it.
- **The apparent partition error is mostly an EXPOSURE-term problem** (§4h, §5):
  on the three tables where the exposure is measured or directly known the bias is
  −0.05 to +0.16, while the estimated-`K_om` rows sit at +0.29. Li 2019's own soil
  half (n=376) runs the model HIGH, opposite to its hydroponic half (n=29) — so
  the anchor decision is now three tables against and one for, and that one is the
  outlier. Treat the S1 hydroponic deficit as an open anomaly, not as evidence
  about `K_PW`.
- **In-planta metabolism was measured for the first time** (§4i): carbamazepine's
  leaf parent fraction is 0.49 on average and 0.17 in lamb's lettuce, despite a
  soil DT50 > 1000 d. `γ = 0` is not a safe default for neutral organics, and a
  half-life fitted to leaf data absorbs the rice-driver mismatch — which is a
  caveat the Ge 2017 ≈7 d figure must now carry.
- Both a-priori inputs are now measured directly and **both are biased low**:
  the partition by the offset above, and TSCF by **−0.221** on a 0–1 scale
  (§4e, the first test of TSCF in this repo without the plant model in between).
- **The weak-electrolyte path has been tested, and it is no longer merely
  "unvalidated" — it is BOUNDED** (§4l). It shipped labelled *structural
  capability, no measured dataset exists*; that was true of **rice** only. §4e
  scores 30 of Table A 3's 97 rows and holds the other **67 ionisable ones** back
  as out of scope — and the port is precisely what extends the scope to them, so
  the test data were already in the repo. Verdict: **direction supported,
  magnitude refuted.** Measured transfer does rise with the neutral fraction
  (Spearman +0.480, n=67, the port's first empirical support of any kind), and
  switching speciation on nearly doubles the model's rank correlation
  (+0.284 → +0.520) — but its influx conductance `Φ` moves ~1.6e4-fold across
  this table where the measurements move ~3-fold, so it under-delivers badly
  (bias +0.023 → −0.203) and at `f_n < 1e−3` predicts nothing where the measured
  mean is 0.127. **The two metrics are not equally solid**: bootstrapped, the
  rank gain survives 94 % of resamples and the RMSE loss only 82 %, so the
  ordering gain is the claim and *"speciation makes the fit worse"* is a
  tendency. The cause is structural and already known one directory away — the
  model's only entry is transmembrane, while real ions arrive apoplastically and
  the PFAS side needed a **fitted carrier** for the same reason. Usable for the
  direction of a speciation effect, not its size, and not below `f_n ≈ 0.1`.
- The error structure is interpretable and cross-validated between the two
  datasets: steep in half-life where the endpoint accumulates (leaf), flat where it
  equilibrates (root). The Ge leaf residual points to a specific in-planta
  half-life of ≈ 7 days — a testable prediction, not a calibration.
- **Air exchange is implemented** (`src/plant_air.py`, §3b), closing the largest
  structural gap in the neutral base: the leaf now has both of its non-growth
  sinks. It is opt-in (`simulate_neutral(air=True)`) because it needs `K_AW` and a
  molar mass, which the strict Kow-only a-priori run does not use — so the 0.281
  and 0.783 above are unaffected, and so is every PFAS number.
- **Hwang 2017 (§4c) is a diagnosis, not a score.** Its RMSEs (0.610 fw / 0.726 dw)
  are not validation results: the unstated fresh/dry basis spans the verdict and
  the two readings fail on opposite organs. When it was run it looked like an
  independent corroboration of the Brunetti discrepancy — a second sighting of
  "the Briggs root partition is too low for lipophilic compounds in soil-grown
  plants". **§4f then superseded Brunetti**, so that synthesis no longer stands:
  Hwang's fw root exceedance is now a single unresolved sighting whose basis is
  itself unstated, not half of a pattern.
- Remaining gaps: the **grain compartment is untested** (no suitable dataset
  exists); the root lipid is corroborated by measured *cereal* roots but stem and
  leaf are still Trapp's soybean figures, and no rice organ lipid has been
  measured; **tissue specific surface areas** are ratios, not measurements, which
  bounds how far an absolute volatilisation flux can be trusted.
- The **open anomaly is now Li 2019's hydroponic half** (§4d, §5), not the
  partition core. Earlier revisions of this document nominated the Brunetti
  `K_RW` 13.3 vs Briggs ≈1.0 disagreement as the best-evidenced open question
  against `K_PW`; **§4f retired that** by measuring 1.10 for the same compound,
  and §4h then showed the deficit is confined to one *exposure route* rather than
  being a property of the partition. What is left unexplained is a −0.43 offset
  on 29 hydroponic rows that no subgroup accounts for, against a rice table
  (Liu) that is unbiased over the same range. Settling it needs a hydroponic
  **rice** root RCF above log Kow 3.5 — the one cell no table in hand covers.
