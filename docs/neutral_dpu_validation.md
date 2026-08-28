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

**Flat in exposure time, monotone in log Kow** — and inside the high-Kow cell
alone, longer exposures are no better (−0.586 under 3 d vs −0.517 over). So
non-equilibrium is ruled out as the driver, and the deficit sits in the
**lipophilic sorption term**: it vanishes at low Kow, where the water floor `W`
dominates and there is no lipid term to be wrong about.

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

**Result: log10 RMSE 0.191 (n = 21), nothing fitted — the best a-priori result in
this repo.** For scale, on the same path: Liu 2023 0.281, Li 2019 0.598.

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
| **Kodešová 2019** `10.1007/s11356-019-04333-9` | **Used — §4f.** Queue A4, closed: the SI arrived. Carbamazepine root partition, 4 plants × 3 soils, with the exposure derived from the paper's own measured soil concentrations and Freundlich isotherms. **log10 RMSE 0.191 — the best a-priori result in the repo** — and it does two things to the open questions: it **supersedes the Brunetti sighting** (measures 1.10 where Brunetti calibrated 13.3 for the same compound) and it **votes against restoring the Briggs anchor**. |
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
| **Briggs 1983** `10.1002/ps.2780140506` | Shoot-distribution companion to the 1982 paper; not yet mined. |

### Is the root partition too low? The evidence splits by lipophilicity

This looked, until the A4 SI arrived, like one accumulating result. It is not.
The sightings do not agree, and where they disagree is informative.

| sighting | log Kow | direction and size | weight |
|---|---|---|---|
| **Li 2019** (§4d), 29 rows, 11 species | −0.5 → 5.4 | model **low**, and the deficit *grows* with Kow: −0.03 below 2, −0.69 above 4.5 | strongest — n = 29, every species the same way |
| **Hwang 2017** (§4c), fresh-weight reading | 4.01 | measurement 2.8–10.4× **above** the `K_PW` ceiling | one compound, lettuce, unstated basis |
| **the internal anchor** (§4d) | all | shipped `L·a` is 2.48× under Briggs' own | not a sighting — the *cause* of part of the above |
| **Kodešová 2019** (§4f), 21 rows, 4 species | 2.25 | model **high** by ~1.5× — a **counter-example** | cleanest exposure in the repo |
| **Brunetti 2021** calibrated pea `K_RW` = 13.3 | 2.25 | looked ~8× low; **now ~12× above a direct measurement of the same compound** | calibrated posterior, largely **superseded** by §4f |

Two of the five moved when the A4 SI arrived, and both moved the same way:
**Brunetti is no longer evidence about the partition core** — §4f measures 1.10
for the compound Brunetti calibrated 13.3 for — and **Kodešová actively opposes**
the group at moderate lipophilicity.

What survives is narrower, and better posed than "the root partition is too low":

> **The deficit is confined to the lipophilic end.** Below log Kow ≈ 2 there is
> essentially no bias (Li 2019: −0.03; Kodešová at 2.25 is *positive*). Above 4.5
> it reaches −0.69. That is the signature of a missing **sorption** term, not a
> wrong water floor and not a global scale error — and §4d already showed it is
> not kinetic either.

A sorbing phase the neutral composition does not have is the natural reading, and
the PFAS side of this same repo carries exactly one (`f_cw·K_cw`, whole cell wall)
and lists its coefficient as **GAP A**. A measured neutral-organic cell-wall
partition coefficient would serve both paths at once and remains the highest-value
wet-lab item on this side, ahead of the in-planta half-lives.

But note what §4f costs that story: a missing sorption term should not make the
model run *high* anywhere, and at log Kow 2.25 it does. Either the water floor is
slightly generous for these species, or the ODE's xylem drain is too weak near the
TSCF peak. Both are open, and neither is addressed by adding lipid — which is
also why raising `L` to fix the lipophilic end would be the wrong instrument.

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
  it reaches **0.191**, the best a-priori result here. That table also **reverses
  two open questions**: it supersedes the Brunetti disagreement and it votes
  against restoring the Briggs anchor, opposite to Li 2019. The partition
  deficit is therefore **not global** — it is confined to log Kow above ~3.
- Both a-priori inputs are now measured directly and **both are biased low**:
  the partition by the offset above, and TSCF by **−0.221** on a 0–1 scale
  (§4e, the first test of TSCF in this repo without the plant model in between).
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
  the two readings fail on opposite organs. Its value is that it **independently
  corroborates the Brunetti discrepancy below** — for lipophilic compounds in
  soil-grown plants the Briggs root partition looks too low, now seen in two
  unrelated datasets rather than one.
- Remaining gaps: the **grain compartment is untested** (no suitable dataset
  exists); rice-specific organ lipid contents are still borrowed from soybean;
  **tissue specific surface areas** are ratios, not measurements, which bounds how
  far an absolute volatilisation flux can be trusted; and the **root partition for
  lipophilic compounds** — the Brunetti `K_RW` 13.3 vs a Briggs `K_PW` of ~1.0 for
  carbamazepine, now joined by Hwang's 3–10× root exceedance — is an open question
  against the partition core itself, and after §4c the best-evidenced one.
