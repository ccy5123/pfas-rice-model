# EU biocides from the BAT census, put through the rice model

*Record for `validation/bat_census_biocides.py` + `data_obs/biocides_bat_census.csv`.
Run it (`python validation/bat_census_biocides.py`, ~3 min; `--fast` skips the log Kow
sweep) rather than quoting this file — every number below is printed by that script.*

**한 줄 요약.** BAT 보고서에 나오는 물질들을 이 저장소의 **중성 유기물 경로**
(`model_api.simulate_neutral`)에 넣어 벼 뿌리/짚/낟알 축적을 계산했다. 24개 물질이
돌아가지만, 그중 11개만 이 저장소가 측정 데이터를 가진 log Kow 구간 안에 있고,
5개는 이온화 때문에 제외되며, 검증 가능한 물질은 **단 2개**(프로피코나졸·트리클로산)다.
**예측이지 검증이 아니다.**

---

## What was done, and why it is a legitimate thing to do

`REPORT_bat_census.md` — a user-supplied document from a separate project, **not
in this repository** — records a regulatory screening tool (BAT) being run over
every EU biocidal active substance with a published bioaccumulation opinion, and
scored against those opinions. Its endpoint is a **fish** bioconcentration factor.
This repo models accumulation in **rice**. Different organism, different endpoint,
different regulation — so only one thing carries across, and it is the thing the
report spent most of its effort on:

> a **log Kow of the uncharged form**, traced to a named source.

That is precisely the one input this repo's neutral-organic path needs
(`docs/neutral_dpu_validation.md`): `K_PW` and `TSCF` are both functions of it and
nothing is fitted. And the report's headline input correction — three
anticoagulants entered with a **distribution ratio D** at environmental pH where
log Kow belongs (its §7.5) — is an error this model would have inherited in
exactly the same way, and in the same direction.

## The scope rule was not invented here

Three limits, all pre-existing, applied before any result was read:

| limit | value | where it comes from |
|---|---|---|
| measured root-partition span | log Kow **−0.66 … 8.70** (560 rows) | the neutral tables in `data_obs/` |
| measured TSCF span (un-ionised) | log Kow **−1.52 … 5.46** (30 rows) | `data_obs/tscf_obs_schriever2020.csv` |
| weak-electrolyte floor | neutral fraction **≥ 0.10** | `docs/neutral_dpu_validation.md` §4l — the path is bounded to *direction*, not magnitude, and not below f_n ≈ 0.1 |

The third is **the report's own §3.0a rule** — "more than 90% ionised at the pH of
natural water" — which the report records as having been written down before the
runs and then never implemented. It is implemented here, at the same threshold,
and it fires on the same substances.

The script recomputes the first two from the shipped tables at run time, so adding
a measured table widens the scope automatically.

## The census

| | n |
|---|---|
| substances named in the report | 44 |
| — carrying a log Kow the report states | **23** |
| — plus propiconazole, whose log Kow comes from this repo's own Liu 2023 row | 24 |
| inside the measured-data span | **11** |
| past the TSCF anchor → **root only** is quotable | 6 |
| beyond every anchor in either model (cholecalciferol 10.24, muscalure 10.61) | 2 |
| **excluded, >90% ionised at pH 7** | 5 |

The 20 substances with no report-stated log Kow are kept in the table with
`status=no_logkow_in_report` and **not filled in from elsewhere** — inventing an
input is the defect the report's own §8.14 is about.

### The transcription is checked before it is used

The report tabulates a **percent ionised at pH 7** for the anticoagulants (§8.9)
rather than their pKa. pKa is recovered by inverting Henderson–Hasselbalch — and
that inversion is not assumed: it round-trips the six substances whose pKa the
report *states* to within **0.002 percentage points**, which is what pins pH 7 as
its basis. Had that failed, every derived pKa and the whole ionisation screen
would have been wrong in the same direction with nothing to catch it.

## What the model says

Rice BAF at a held exposure, γ = 0, so leaf/straw/grain are **upper bounds**:

| substance | log Kow | TSCF | K_PW root | root | straw | grain | scope |
|---|---:|---:|---:|---:|---:|---:|---|
| Chlorophacinone | 2.46 | 6.5e-01 | 1.86 | 1.19 | 35.6 | 9.76 | inside |
| Propiconazole | 3.72 | 1.7e-01 | 9.83 | 8.69 | 15.4 | 3.37 | inside |
| Tralopyril | 4.08 | 9.0e-02 | 17.8 | 15.7 | 8.87 | 1.66 | inside |
| Permethrin | 4.60 | 3.0e-02 | 43.4 | 42.4 | 4.01 | 0.482 | inside |
| Triclosan | 4.76 | 2.1e-02 | 57.3 | 56.3 | 2.87 | 0.289 | inside |
| Triflumuron | 4.91 | 1.4e-02 | 74.5 | 73.7 | 2.08 | 0.173 | inside |
| Epsilon-metofluthrin | 5.00 | 1.1e-02 | 87.3 | 86.5 | 1.68 | 0.124 | inside |
| Hexaflumuron | 5.68 | 1.5e-03 | 289 | 287 | 0.229 | 0.0058 | root only |
| Difethialone | 6.29 | 1.9e-04 | 851 | 708 | 0.0192 | 1.6e-04 | root only |
| Empenthrin | 6.30 | 1.8e-04 | 867 | 716 | 0.0184 | 1.5e-04 | root only |
| Cholecalciferol | 10.24 | 1.4e-13 | 9.4e+05 | 1550 | 2.3e-14 | 1.7e-19 | extrapolation |

The pattern is one sentence: **everything the opinions call bioaccumulative loads
the root and does not reach the grain**, because rice shoot loading is gated by a
TSCF bell peaking at log Kow 1.78, and every one of those substances sits at
TSCF < 1e-3.

## The two substances that can actually be scored

The only ones in the report for which this repo already holds a measured plant
root row. A-priori — log Kow in, partition out, nothing fitted. Scored as an
**equilibrium** (`docs/neutral_dpu_validation.md` §4g: reading a 24 h–26 d
partition off a 120-day rice season imposes a Kow-dependent model-side discount):

| substance | source | measured | model K_PW | error |
|---|---|---|---|---|
| **Propiconazole** | Liu 2023, **rice** root RCF, hydroponic 72 h | 9.32 | 9.83 | **+0.023 log (1.05×)** |
| **Triclosan** | Li 2019 soil table, radish/carrot root, 14 rows | 5.4–58.9 (median 14.0) | 61.5 | +0.625 log, RMSE 0.689 |

Two substances is not a validation. What they do is fix the sign: on the two
places anything can be checked, the model is not wrong by an order of magnitude.
The triclosan bias is the same sign as, and larger than, that whole soil table's
documented +0.260 (§4h) — and triclosan is an **acid**, so part of the excess is
speciation the equilibrium `K_PW` cannot see (next section).

A free by-product: Li 2019's own log Kow for triclosan is **4.8** against the
report's audited **4.76** — two independent sourcings 0.04 log apart, which is the
one place the report's property audit can be checked against something already
here.

## Where the two models agree, and where they part

**(a) Both are one-input models — but one had to be probed for it.** The report
needed 217 runs to establish that BAT's fish BCF is a function of the entered log
Kow and of nothing else that varies (its §8.11: difethialone and cyphenothrin, MW
539.5 vs 375.5, solubility 39× apart, Henry 16,000× apart, return the same
number). Here the same statement is **structural** — `K_PW` and `TSCF` are
functions of log Kow, and MW, solubility and Henry appear in no term of the core
ODE — and the two substances come back **bit-identical** (root 707.637 both).

**(b) The one property that is not inert here is the one BAT proved inert.**
Henry's constant reaches the answer through the opt-in air term
(`src/plant_air.py`, identically zero at `K_AW = 0`). Empenthrin (K_AW 1.4e-2)
loses about half its leaf burden to volatilisation; cyphenothrin (9e-7) loses
nothing. So the three substances BAT cannot tell apart are separable here — but
only through a term that is off by default and whose absolute magnitude is
surface-area-limited (docs §plant-air).

**(c) The report's own ionisation sweep, re-run on identical inputs.** Its §8.15
holds DCPP at log Kow 4.60 and moves the pKa alone. Both models normalised to
their own neutral end:

| pKa | f_neutral | Henderson–Hasselbalch | BAT fish BCF | rice root | rice leaf |
|---:|---:|---:|---:|---:|---:|
| 6.000 | 0.091 | 0.091 | 0.245 | **0.768** | 0.507 |
| 8.065 | 0.921 | 0.921 | 0.937 | 0.998 | 0.992 |
| 10.000 | 0.999 | 0.999 | 1.000 | 1.000 | 1.000 |

Both models damp ionisation far below what the chemistry predicts — the report's
§7.6 finding, reached here independently. But **this model damps it harder, and
the root hardest of all**, and the reason is structural rather than a tuning: the
membrane term sets how *fast* the root equilibrates, not the level it equilibrates
to, so a root that still reaches equilibrium inside the season barely registers
the change. **That is why the >90%-ionised group is excluded rather than reported
with a caveat.**

**(d) The two log Kow curves are not the same curve.** BAT's fish BCF peaks near
6.3 and re-crosses its own bioaccumulative threshold *downward* between 8.13 and
9.00 (its §8.11/§8.16). This model's root has no such turnover — the Briggs RCF
has no descending limb — and its straw peaks at the TSCF bell, log Kow **1.75**,
about 4.5 log units below BAT's peak. What the root does instead is **saturate
kinetically**: K_PW / root BAF is 0.92 at log Kow 4, 0.85 at 6.25, 0.055 at 8.25
and 0.001 at 10.5. Above ~7 the season is too short to equilibrate the root, so
the model's root number stops being a partition coefficient and becomes a rate.
That is a limit of this model, stated here rather than discovered later.

**(e) A fish BCF does not rank these substances the way a rice model does.**
Spearman against BAT's Scenario A fish BCF over the 20 substances that have one:
**root +0.467, straw −0.466**. The root agrees because both rise with
lipophilicity; the shoot *anti*-correlates, and cannot do otherwise, because what
reaches a rice shoot is gated by a bell that has already collapsed by log Kow 5.
"Bioaccumulative in fish" is not a statement about grain.

**(f) The report's largest finding is amplified here, not damped.** Entering a
metabolic rate moved BAT's fish BCF by up to 82× (its §8.4). Those same fish kM
half-lives, transplanted into the plant as a **labelled sensitivity**, move this
model's **grain** by up to **5,573×** — because grain and leaf are terminal
accumulators whose only other sink is growth dilution, which goes to zero at
maturity. The quantity the report identifies as its least defensible input is the
one the edible compartment is most sensitive to, and for a plant nobody has
measured it at all.

## What this does not establish — and cannot

* **Not a validation.** No measured rice concentration exists for any of these
  substances. The two checkable ones are one hydroponic rice point and 14
  soil-grown radish/carrot rows.
* **The substances that matter are the ones the model is least entitled to.** Of
  the seven the opinions call bioaccumulative, four are >90% ionised at pH 7 and
  excluded outright (brodifacoum, difenacoum, bromadiolone, flocoumafen), and of
  the remaining three only triclosan is inside the TSCF anchor at all. That mirrors
  the report's own §7.3 finding about its kM QSAR — those same substances at
  Tanimoto 0.19–0.28 from anything it was trained on — reached from a different
  direction: **two independent models are both weakest exactly where the answer
  matters.**
* **The fish kM half-lives are not plant parameters.** They are run once, labelled,
  because the question the report makes unavoidable is worth asking of this model
  too. This repo's own in-planta half-life is itself open (§4i: measured parent
  fractions differ 4.8× between species for one compound in one experiment).
* **γ = 0** makes every leaf and grain number an upper bound.
* **The exposure is a held pore-water concentration**, not a soil inventory. Any
  real use of these numbers needs the soil side (`cwo_profile=`, or HYDRUS) and a
  substance-specific Koc — which for these biocides this repo does not have.

## Nothing was changed to do this

`params/parameters.json`, `simulate()`'s defaults and `reproduce_demo.py`
(log10 RMSE 0.029) are untouched; the neutral path's published a-priori numbers
(Liu 0.281 / Ge 0.783) are untouched. This work adds a data table, a validation
script and its guards, and nothing else.

Guards: `tests/test_bat_census_biocides.py` (14) — the pH-7 round-trip, the
derived-pKa inversion, "no input is invented", the ionisation screen firing on the
report's own group, the two a-priori overlap predictions, the structural
inert-input identity, the air-term separation, the damping comparison, and the
caveats themselves.
