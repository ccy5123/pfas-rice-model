# HANDOFF — neutral-organic path: what is done, what is next

> Session handoff for the next Claude/dev.
> **PRs [#56](https://github.com/ccy5123/pfas-rice-model/pull/56) and
> [#57](https://github.com/ccy5123/pfas-rice-model/pull/57) are MERGED** — `main` is at
> `b3b5004`, full suite re-run on the merged tree: **245 pass, 2 skip**.
> Open: **[#58](https://github.com/ccy5123/pfas-rice-model/pull/58)** (draft, branch
> `claude/latest-handoff-review-da78q1`) — the Expert-only neutral tab + this handoff.
> Scientific records: **`docs/neutral_dpu_validation.md`** (neutral path, §3b air / §4c Hwang) and
> **`docs/twopool_root_exploration.md` §Result 8** (structural merge). Read those first.
> `parameters.json`, `simulate()` and `reproduce_demo` (RMSE 0.029) are **UNCHANGED**
> throughout — everything here is additive or opt-in.

---

## 0. TL;DR

**The neutral-organic path is structurally complete and shipped.** §3's three tasks
(A1 → A3 → A2) are done and merged, and the follow-up that could be done in-repo
(§3.1-3, the app tab) is done too:

- **A1 air exchange** — the leaf's second sink now exists, so a volatile compound is
  no longer an upper bound by construction. Structurally zero at `K_AW = 0`.
- **A3 Hwang 2017** — a **diagnosis, not a score**: the article's unstated fresh/dry
  basis spans the verdict, and the two readings fail on *opposite organs*. Its payoff
  is corroborating the open **Brunetti** root-partition disagreement.
- **A2 `model_api.simulate_neutral`** — reachable from the app/validation,
  bit-identical to the standalone module.
- **The app tab** — Expert-only "🧪 Neutral organics" (PR #58).

**What is left is almost entirely BLOCKED ON PAPERS**, and they are now enumerated with
DOIs in **`docs/literature_db/Acquisition_Queue.csv`** (§4). Note the constraint that
shapes everything below: **this environment cannot fetch them** — every academic host is
blocked by the network egress policy. They have to arrive out-of-band.

---

## 1. What this session delivered

| commit | what |
|---|---|
| `d8e7f9a` | **A1** air exchange: `src/plant_air.py`, opt-in `RiceUptakeModel(air=)`, `simulate_neutral(air=True)` |
| `97dbe75` | **A3** Hwang 2017 as a diagnosis: basis spans the verdict; opposite-organ failure; Brunetti corroborated |
| `50b5586` | **A2** `model_api.simulate_neutral(log_kow, …)` + drift guards (purely additive, 83 lines) |
| `bb70b8a` | handoff rewritten for the A1/A3/A2 outcome |
| *(merged)* | `8169257` = PR #56, `b3b5004` = PR #57 |
| `c6e9d8a` | **app**: Expert-only "🧪 Neutral organics" tab (PR #58) |
| `494acc4` | **`docs/literature_db/Acquisition_Queue.csv`** — the wanted papers, with DOIs and `doi_status` |

**New files**: `src/plant_air.py`, `validation/hwang2017_lettuce.py`,
`tests/test_plant_air.py` (18), `tests/test_hwang2017.py` (8),
`docs/literature_db/Acquisition_Queue.csv`.

Full suite on merged `main`: **246 collected, 245 pass, 2 skip** (+29 vs the 217 this
work started from). The 2 skips are optional deps — `emcee` and `sci-adk`.

**⚠️ CI does not test any of this.** `.github/workflows/rigor.yml` runs only
`tests/test_sci_adk_rigor.py`, so a green check says nothing about the model or the app.
Run the full suite locally (~13 min) before claiming green.

**Prior session** (merged as #56): the structural merge and the neutral path's first
implementation — `1a707e6`, `1028e3e`, `3ee864e`, `a8ba9e7`, `f7b870e`, `cc4c54c`.

---

## 2. Numbers to trust

| result | metric |
|---|---|
| **air exchange**, K_AW ladder (log Kow 2.42) | leaf BAF **177 → 0.0025** across K_AW 0 → 0.1; **root invariant**; `K_AW=0` bit-identical to air-off |
| **Hwang 2017**, both basis readings | 0.610 fw / 0.726 dw — **not validation scores**; fw fails the root (0.711) and fits the leaf (0.489), dw the reverse (0.393 / 0.948) |
| Hwang, Table 1 internal consistency | `whole` = mass-weighted mean at root fraction **5.4 ± 0.9 %**, identical at every sampling |
| Hwang, root vs the Briggs ceiling | `K_PW` = 15.8 L/kg; measured root **2.8–10.4× above** (fw) or under it (dw) |
| neutral root partition, Liu 2023, **nothing fitted** | log10 RMSE **0.281** (n=14, log Kow −0.66→4.4); 11/14 within 1.5× |
| neutral per-organ TF, Ge 2017, **nothing fitted** | log10 RMSE **0.783** (n=6); minimum **0.210** at a 7-day half-life |
| partition adapter vs Briggs RCF | machine precision, 1.3e-16 |
| merge, Yamazaki in-sample | 0.301 (root 0.153); un-refit transfer 0.316 vs 0.278 reference |
| merge, Tang per-organ OOS | **0.801** (stalk 0.61, leaf **0.28**) vs 1.398 pass-through stem |
| PFAS a-priori, for comparison | 0.84–0.95 — **and it still has fitted transport behind it** |

**Carry these caveats forward on every claim.** Neutral: two studies, one crop, 20
compound-tissue pairs; the **grain compartment is untested**; rice organ lipid is
borrowed from soybean; the Ge leaf residual is confounded with the missing half-life,
so 0.783 is an *upper bound* on transport error. Merge: Yamazaki in-sample fit; does
not beat single-pool lipid loading on Tang (0.516); the residual is entirely the
endosperm; Tang is three C5–C8 congeners so the long-chain root decoupling is still
unexercised. The **promotion decision did not move** — the gate is still the wet-lab
assay in `docs/twopool_kseq_mechanism.md` §5.

Two more. **Air exchange**: the equations are the derivation's, but the flux scales with
each tissue's specific SURFACE AREA, which this repo has only ever used as a leaf/grain
*ratio* for the xylem split — so absolute volatilisation magnitudes are
order-of-magnitude until the areas are measured, and the shipped stem area is 0.
**Hwang**: lettuce not rice, one compound, an unstated fresh/dry basis and a
reconstructed growth curve — its RMSEs are a diagnosis, never a validation score.

---

## 3. Next tasks (in-silico) — ✅ A1, A3 and A2 are all DONE

> All three tasks below were completed this session; each carries its outcome inline.
> The remaining work is in §3.1 (code-shaped follow-ups these three surfaced) and
> §4 (blocked on data or experiment).

### 3.1 Follow-ups these three surfaced

Ranked. None is large, and none is blocking.

1. **Measure or source tissue specific surface areas** [m²/kg fw] for rice. This is
   the single input that bounds A1: volatilisation flux is linear in it, and the
   values in the repo were never calibrated as areas (`AirExchange(S=...)` overrides
   them; the stem's is 0, so its cuticular term is inert). Cheapest real gain on the
   air side, and it needs a source, not a fit.
2. **Ask Hwang's authors which basis Table 1 uses.** One email closes §4c: it decides
   whether their root measurement sits above or below the model's partition ceiling.
   Note it does NOT dissolve the discrepancy — the two readings fail on opposite
   organs — but it says which organ to chase.
3. ~~A Streamlit tab for the neutral path~~ ✅ **DONE** — the **EXPERT-ONLY**
   "🧪 Neutral organics" tab (`ui/expert.py` tabs[8]). The Simple/Expert question was
   settled the conservative way: Simple is congener-driven and symbol-free while a
   neutral compound is a log Kow, and the neutral **grain is untested**, so putting
   absolute numbers for it on a general-audience screen would be exactly backwards.
   Verified with headless Streamlit + Playwright. Docs: `docs/visualization_tool.md`.
4. **Particle deposition (`eq:Qdep`)** is still unimplemented — deliberately, as a
   separate atmospheric-deposition pathway rather than plant-air equilibrium
   exchange. Only worth doing if a use case needs airborne particulate input.
5. **`NStemLeafModel` has no air hook.** Fine today (it is a PFAS-side model and PFAS
   have no air exchange), but if the neutral path ever needs the redistributed shoot,
   the hook has to be added there too.

### 3.2 The three completed tasks, with outcomes

### A1. Air exchange — volatilisation + gaseous uptake  ✅ **DONE**

> **Implemented**: `src/plant_air.py` (the equation set), an optional `air=` field on
> `RiceUptakeModel` (default `None` → the terms are never evaluated),
> `simulate_neutral(air=True)`, `tests/test_plant_air.py` (18 tests), and §3b of both
> `validation/neutral_dpu_validation.py` and `docs/neutral_dpu_validation.md`.
> `parameters.json`, `simulate()` and `reproduce_demo` (**RMSE 0.029, re-verified**)
> are unchanged. `k_aw_warning` now names the remedy instead of refusing.
>
> Results: the leaf, an unbounded terminal accumulator without metabolism, is now
> bounded by volatilisation — leaf BAF 177 → 0.0025 across `K_AW` 0 → 0.1, while the
> **root is invariant** (no volatilisation below ground) and `K_AW = 0` is bit-identical
> to air-off. The ladder also *checks* the old judgement-call warning threshold
> (`K_AW > 1e-4`): the physics puts the crossover just above it (leaf t½ 787 d at 1e-4
> vs 8 d at 1e-3), so it errs toward flagging early.
>
> **The residual limit is the surface area, not the equations** — the flux scales with
> each tissue's specific surface `S`, which this repo has only ever used as a leaf/grain
> *ratio* for the xylem split, so absolute magnitudes are order-of-magnitude until `S`
> is measured. `AirExchange(S=...)` overrides it; the shipped stem `S = 0` leaves the
> stem term inert (pinned by a test). Particle deposition (`eq:Qdep`) is still not
> implemented, deliberately — it is a separate deposition pathway, not plant–air
> equilibrium exchange.

**Original brief, for the record.** The largest structural gap in the neutral base,
and fully specified in-repo. The core ODE had **no air terms at all**. Today
`neutral_dpu.k_aw_warning()` merely *refuses* a high-`K_AW` compound. That is
defensible for PFAS (`K_AW ≈ 0` — CLAUDE.md §2 turns air exchange off deliberately)
but it truncates the neutral path, which is meant to cover volatile organics.

Everything needed is written out in `docs/dpu_model_summary_corrected.tex`
§`sec:permeability` (line ~263) onward — note these equations assume **SI base units**
(m, s, g/mol), so unit conversion into the repo's day/L/kg convention is part of the job:

| eq. label | line | what |
|---|---|---|
| `eq:Pc` | 268 | cuticle permeability `10^(0.704 log Kow − 11.2)` |
| `eq:Pair` | 273 | air boundary-layer conductance |
| `eq:Paqua` | 278 | aqueous/apoplast permeability |
| `eq:Pctot` | 284 | the three in **series** |
| `eq:Ps` | 289 | stomatal conductance, tied to the transpiration stream |
| `eq:Csat` | 293 | saturation water-vapour concentration |
| `eq:Pp` | 299 | cuticle + stomata in **parallel** |
| `eq:Kpa` | 314 | plant–air partition `K_PA = K_PW·ρ/K_AW` |
| `eq:Qvol` | 322 | volatilisation flux |
| `eq:Qgas` | 330 | gaseous uptake |

Modelling assumptions the tex states explicitly: **no volatilisation from roots**;
**stem via cuticle only**; **leaf and fruit via cuticle + stomata in parallel**.

Watch the singularity: `eq:Ps` carries `1/(1−φ)` which blows up as relative humidity
→ 1. The tex flags it — the product stays finite because `Q_TP → 0` in that limit, but
guard it numerically (cap `φ < 1`).

**Design constraint:** keep it off the PFAS path. The cleanest route is an optional
term in the neutral module (or an opt-in flag on the core) that is identically zero
when `K_AW = 0`, so every existing PFAS number is untouched — the same discipline as
`k_seq=0` in the merge. Then replace `k_aw_warning`'s refusal with an actual run, and
add a test that a volatile compound now loses mass through the leaf while a
non-volatile one is bit-identical to today.

**Why it matters beyond completeness:** §3 of the validation shows the leaf is an
unbounded terminal accumulator without metabolism *or* volatilisation. Right now only
one of those two sinks exists, so every volatile compound is an upper bound by
construction.

### A3. Hwang 2017 — a second per-organ test, with a growth + exposure model  ✅ **DONE**

> **Run**: `validation/hwang2017_lettuce.py` + `tests/test_hwang2017.py` (8) + §4c of
> `docs/neutral_dpu_validation.md`. Outcome is a **diagnosis, not a score** — and the
> pre-registered caveat below is exactly what decided it.
>
> **The basis was NOT resolvable.** This session had no way to re-read the article
> (the environment's network policy blocks every academic host), but that changes
> nothing: the brief already recorded that the basis *is not stated in the text*, so
> the previous session had the PDF and found it absent. It is a question for the
> authors, not for a re-read.
>
> **What the run establishes.** (i) Table 1 is internally consistent on ONE basis —
> `whole` is the mass-weighted mean of leaf and root at a root mass fraction of
> **5.4 ± 0.9 %**, identical at every sampling — and 5.4 % is characteristic of
> **fresh** weight for lettuce (dry would be ~11 %). (ii) The modelled root cannot
> exceed its equilibrium partition, `K_PW = 15.8 L/kg`, so the basis flips the
> verdict: read fresh the measurement is **2.8–10.4× above that ceiling**
> (unreachable), read dry it sits **under** it. (iii) The sharpest result — **the two
> readings fail on OPPOSITE organs** (fw: leaf 0.489 / root 0.711; dw: root 0.393 /
> leaf 0.948), so the basis decides *where* the model is wrong, not *whether*; the
> discrepancy is not a units artifact awaiting a footnote. (iv) Soil contact cannot
> explain the root exceedance (it would need 12–49 % of washed root mass to be soil).
>
> **Why it was worth doing anyway**: the fresh-reading root exceedance is the **same
> direction and magnitude as the open Brunetti 2021 `K_RW` disagreement**, so the
> "Briggs root partition is too low for lipophilic compounds in soil-grown plants"
> problem now has **two independent sightings** instead of one. That promotes it to
> the best-evidenced open question against the partition core.
>
> **A trap named in the code and docs**: the half-life scan improves the dry reading
> (0.73 → 0.30) and worsens the fresh one, so the fit "prefers" dry weight. Using
> that to choose the basis would be **circular** and would override the only
> non-circular evidence (i). Do not let a later session quietly do it.
>
> **Deliberately NOT shipped as a `data_obs/` CSV**: the shared `--obs` harness would
> run it on the *rice* drivers (120 d, constant exposure) and return a silently
> meaningless number. The data stays in the script's constants.
>
> `Tp` was scanned, never adopted — the authors' caveat (§6) holds.

**Original brief, for the record.**

`10.1371/journal.pone.0172254` (PLOS ONE, open access). Chlorpyrifos in lettuce from
treated soil. **All values needed are recorded in §6 below**, because the uploaded
PDFs live in an ephemeral session scratchpad and will be gone.

Why it is worth doing: it is the only **time-resolved per-organ** neutral dataset to
hand (3 sampling times × 2 soil levels), it is **lipophilic** (log Kow 4.01, where Ge's
difenoconazole sits), and unusually it supplies a measured `Kd` — so soil concentration
converts to a **pore-water exposure**, which is what the model actually wants. It also
reports a plant half-life, letting the Ge diagnosis be checked where metabolism is not
a free unknown.

Honest limits to state up front: **lettuce, not rice**; one compound; the fresh/dry
basis of Table 1 is **not stated in the text** — resolve it before trusting the
absolute numbers (this is exactly the trap that produced the Tang fresh/dry artifact
and the lipid-basis error in `a8ba9e7`); and the plant half-life is an *assumed* model
input, not their measurement (§6).

### A2. Expose the neutral path through `model_api`  ✅ **DONE**

> **`model_api.simulate_neutral(log_kow, ...)`** — mirrors the `simulate_nstem_leaf` /
> `simulate_twopool_seq` pattern: the same driver machinery (`drivers=`, `biomass=`,
> `measured_forcing=`) and the same result-dict contract as `simulate()`
> (`conc`/`baf`/`baf_final`/`straw`/`straw_baf`/`tf_final`/`cwo_ref`/`season`/`M`/
> `params`), plus the neutral-only diagnostics (`K_PW`, `TSCF`, `rcf_briggs`,
> `air_summary`). Air exchange and phloem stay opt-in; no default changes anywhere.
>
> The first argument is a **log Kow**, not a congener name — a neutral compound has
> no congener analogue, and lipophilicity is the one input the QSPRs need.
>
> **Drift guards** (`tests/test_model_api.py`): the wrapper is bit-identical to
> `neutral_dpu.simulate_neutral` on the same forcings, through BOTH the built-in
> path and `drivers=`, so the published a-priori numbers cannot silently stop
> describing what the API returns. Also pinned: the `simulate()`-shaped contract,
> `N = 0` / `e^N = 1` (ionic machinery off), and that air exchange is opt-in and
> bit-identical at `K_AW = 0`.
>
> **Streamlit tab: still not built, deliberately** — as the brief says, that needs a
> decision about where a neutral compound belongs in the Simple/Expert split, and
> the app is PFAS-only today. The API is what unblocks it.

**Original brief, for the record.**

`neutral_dpu` is currently standalone — every other capability in the repo is reachable
through `model_api`, so the app and other validation cannot use it. Small, low-risk,
and it makes the work usable. Mirror the `simulate_nstem_leaf` / `simulate_twopool_seq`
pattern: same result-dict shape, `drivers=` support, no change to any default. Whether
to add a Streamlit tab is a separate call — the app is currently PFAS-only and the
Simple/Expert split would need a decision about where a neutral compound belongs.

---

## 4. Blocked on data or experiment (not code)

**Fetchable literature — the full list is now a file**:
**`docs/literature_db/Acquisition_Queue.csv`** (priority, DOI, `doi_status`, what to
extract from each, and why it matters). Papers already in hand are deliberately not in
it — those are in `docs/neutral_dpu_validation.md` §5.

Headlines, so this section still reads on its own:

| # | what it unlocks | DOI | status |
|---|---|---|---|
| **A1** | **rice ROOT lipid** — `K_PW` is linear in it, so it conditions the 0.281 directly. Get **the SI** | `10.1016/j.envint.2019.02.020` | search-only |
| A2 | rice STRAW lipid, basis stated (3.4 % dw lipophilic extractives) — **open access** | `10.3389/fpls.2022.868319` | search-only |
| A3 | the 97 per-compound TSCF values — **SI only**, the article is in hand | `10.1016/j.scitotenv.2020.136667` | verified |
| A4 | the measured concentrations behind Brunetti's `K_RW` = 13.3 | `10.1007/s11356-019-04333-9` | unverified |
| A5 | the bromacil data behind Trapp 1994 (figure-only there) | *(none — McFarlane 1987, JEQ 16(4):372–376)* | none |
| B1/B2 | the 2 long-chain PFAS papers never obtained (5 of 7 were) | `10.1021/acs.est.5c11716`, `10.1139/er-2025-0116` | verified |
| C1–C3 | rice **grain** neutral data, rice **specific surface areas**, a fallback root-lipid source | *(no paper identified — search specs)* | none |

**⚠ Two things to carry with this list.**

*It cannot be fetched from inside a session.* Every academic host is blocked by this
environment's network egress policy — `journals.plos.org`, PMC, Europe PMC, `doi.org`
and Crossref were all checked and refused. `WebSearch` works (that is how A1/A2 were
found) but `WebFetch` does not, so **nothing here was verified at source**: `doi_status`
says exactly how far each one got. These have to arrive out-of-band, the way
`DPU4OC.zip` did.

*The definition of "lipid" matters more than the number.* The `L` in
`K_PW = W + L·a·Kow^b` is an **operationally defined octanol-like phase** calibrated so
the expression reproduces measured RCF — Briggs' own barley anchor `L·a = 10^−1.52`
implies `L ≈ 2.5 %` of **fresh** weight. Measured rice straw spans **0.14–1.0 % dw** as
crude fat (ether extract) but **3.4 % dw** as dichloromethane lipophilic extractives
(A2): a 3–24× spread that is a difference of *definition*, not of measurement quality.
Substituting a proximate crude-fat value therefore changes what the parameter **means**,
not just its value, and could be worse than the current cited-soybean state. That is why
**A1 outranks A2** — its lipid is defined in an RCF-modelling context, the same
operational meaning as ours. Same class of trap as the three basis errors this repo has
already hit.

**Wet lab** (cannot be fetched): in-planta half-lives for the Ge compounds (the model
predicts ≈7 d — a falsifiable prediction); rice-root cell-wall / Fe-Mn-plaque
batch-sorption + desorption assay across chain length × head group (**the `k_seq`
promotion gate**); per-congener xylem-sap / root-water ratio (direct `f_xy`); a direct
`K_cw` measurement (GAP A).

---

## 5. Housekeeping

- **#56 and #57 are merged**; `main` = `b3b5004`. **#58 is an open draft** (the Expert
  neutral tab + this handoff). Marking it ready or merging is the user's call — do not
  flip it unasked.
- **CI only runs `tests/test_sci_adk_rigor.py`** (`.github/workflows/rigor.yml`). A green
  check does NOT mean a change is tested — none of this session's work would have been
  caught by it. Run `pytest -q` locally (~13 min) before claiming green.
- **Test counts**: CLAUDE.md still says 217; the tree is at **246 collected, 245 pass, 2
  skip**. Worth refreshing next time CLAUDE.md is edited anyway.
- **The network egress policy blocks every academic host** — PLOS, PMC, Europe PMC,
  `doi.org`, Crossref all refused. `WebSearch` works; `WebFetch` does not. So no paper
  can be verified at source from inside a session, and §4's `doi_status` column exists
  for exactly that reason. Do not silently upgrade a `search-only` DOI to verified.
- `docs/HANDOFF_BAF_twopool.md` §6 item ④ records the merge outcome; keep it in sync if
  the two-pool story moves again.

---

## 6. Recorded data for A3 (the PDFs are ephemeral)

From Hwang, Lee & Kim 2017, PLOS ONE 12(2):e0172254 — **transcribed from the article,
not re-derived**. Chlorpyrifos (CP) in lettuce grown on treated soil.

**Table 1 — measured residues (mg/kg), mean ± SD, n = 3.** Fresh/dry basis **not
stated in the text — resolve before use.**

| soil level (mg/kg) | day | leaf | root | whole |
|---|---|---|---|---|
| 10 | 21 | 0.5 ± 0.05 | 7.5 ± 0.66 | 0.9 ± 0.07 |
| 10 | 30 | 0.4 ± 0.03 | 5.6 ± 0.38 | 0.6 ± 0.01 |
| 10 | 40 | 0.1 ± 0.00 | 1.8 ± 0.12 | 0.2 ± 0.01 |
| 20 | 21 | 0.8 ± 0.05 | 2.1 ± 0.14 | 0.8 ± 0.04 |
| 20 | 30 | 0.3 ± 0.02 | 3.6 ± 0.05 | 0.5 ± 0.02 |
| 20 | 40 | 0.1 ± 0.00 | 0.4 ± 0.00 | 0.1 ± 0.00 |

**Table 2 — model parameters**, by nominal treatment (10 / 20 mg/kg):

| parameter | unit | 10 | 20 | provenance |
|---|---|---|---|---|
| initial soil concentration `C0` | mg/kg | 15.2 | 24.9 | measured |
| half-life in soil `T` | d | 17.2 | 7.9 | see caveat |
| **half-life in plant `Tp`** | d | **8.7** | **8.0** | **see caveat** |
| soil–water distribution `Kd` | mL/g | 82.1 | 82.1 | measured (adsorption isotherm) |
| `Koc` | – | 2218.9 | 2218.9 | derived |
| `Kow` | – | 1.02×10⁴ (**log Kow 4.01**) | same | literature |
| transpiration stream `Qw` | mL/d | 46.8 | 46.8 | measured (potometer) |
| log initial plant weight `Ig` | g | 0.3062 | 0.3092 | fitted growth curve (r > 0.99) |
| plant growth constant `Kg` | – | 1.1020 | 1.2031 | fitted growth curve |

**Caveat, stated by the authors:** *"Most of the parameters were obtained from the
laboratory experiments, except for half-life values of CP in soil and lettuce."* So
`Tp` is **not** their measurement — do not present it as an independent confirmation of
the ≈7 d the Ge sensitivity scan predicts. It is a plausible literature-scale value
sitting near that prediction, and no more than that.

Their exposure model is `Ce(t) = C0·(1/2)^(t/T) / Kd`, i.e. soil concentration decays
first-order and divides by `Kd` to give the soil-solution concentration — that is the
`C_w^o(t)` driver to feed `simulate_neutral(drivers=...)`.

---

## 7. How to resume

```bash
pip install -r requirements.txt

# the neutral path
python validation/neutral_dpu_validation.py                                   # structural checks + the air ladder (3b)
python validation/neutral_dpu_validation.py --obs data_obs/neutral_obs_liu2023.csv   # 0.281
python validation/neutral_dpu_validation.py --obs data_obs/neutral_obs_ge2017.csv    # 0.783
python src/neutral_dpu.py                                                     # per-compound demo
python src/plant_air.py                                                       # permeabilities + volatilisation half-lives
python validation/hwang2017_lettuce.py                                        # A3 — read the VERDICT block

# the merge
python validation/twopool_nstem_merge.py --cached      # reuse the fit (~1 min)

# in code / in the app
python -c "import sys; sys.path.insert(0,'src'); import model_api as a; \
           print(a.simulate_neutral(2.45, name='carbamazepine', half_life=7.0)['baf_final'])"
streamlit run app.py            # Expert mode -> the "🧪 Neutral organics" tab

pytest -q                                              # the full suite, ~13 min
```

**Resume prompt — pick by what arrived.**

*If papers arrived* (the likely case — see §4 and
`docs/literature_db/Acquisition_Queue.csv`):

> The paper(s) in `docs/literature_db/Acquisition_Queue.csv` row <A1/A2/…> have arrived
> at <path>. Extract the value, record the source and its **stated weight basis** in the
> literature DB, and only then decide whether it should replace the current
> `neutral_dpu` default. Read the DEFINITION NOTE at the bottom of that CSV first: a
> proximate crude-fat number is **not** interchangeable with the `L` in
> `K_PW = W + L·a·Kow^b`, which is an operationally calibrated octanol-like phase. If the
> value changes `K_PW`, re-run the Liu 2023 and Ge 2017 a-priori comparisons and report
> the movement honestly — those two numbers are the repo's only genuine a-priori results.

*If nothing arrived*, the honest answer is that the neutral arc is at a natural stopping
point and the remaining in-repo items (§3.1-4, §3.1-5) are low-value. Better uses of a
session: ask the user which of §4 they can supply, or move to a different track
(`docs/HANDOFF_BAF_twopool.md`, the PFAS side).
