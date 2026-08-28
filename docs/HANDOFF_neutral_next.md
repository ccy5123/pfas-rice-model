# HANDOFF — neutral-organic path: what is done, what is next

> Session handoff for the next Claude/dev.
> Branch: **`claude/latest-handoff-docs-aq5b80`** · PR: **[#56](https://github.com/ccy5123/pfas-rice-model/pull/56)** (draft, base `main`).
> Scientific records: **`docs/neutral_dpu_validation.md`** (neutral path) and
> **`docs/twopool_root_exploration.md` §Result 8** (structural merge). Read those first.
> `parameters.json`, `simulate()` and `reproduce_demo` (RMSE 0.029) are **UNCHANGED**
> throughout — everything here is additive or opt-in.

---

## 0. TL;DR

Two things landed. (a) The **structural merge** — the two-pool sequestration root
grafted onto the redistributed shoot — which confirmed the standing diagnosis that
the two-pool's Tang failure was a *shoot* artifact (per-organ OOS **1.398 → 0.801**).
(b) The **neutral-organic path**, the framework's Briggs/Kow base, implemented for the
first time and validated **a-priori** against two measured rice datasets (**0.281**
root partition, **0.783** per-organ TF, nothing fitted).

**What is next is code, not data**: three self-contained tasks (§3), of which the
first — air exchange — is the largest remaining structural gap in the neutral base,
and its equations are already written out in the repo.

---

## 1. What this session delivered

| commit | what |
|---|---|
| `1a707e6` | structural merge: `NStemLeafModel(k_seq=, k_rel=)`, `simulate_twopool_nstem`, RHS perf (2× faster) |
| `1028e3e` | Result 8 — merge confirms the Result 7 diagnosis; docs + guards |
| `3ee864e` | neutral path implemented; Briggs anchors verified at source; Trapp lipids; Schriever TSCF |
| `a8ba9e7` | Liu 2023 SI → root partition a-priori 0.281; **fresh-vs-dry-weight correction** |
| `f7b870e` | doc fix: Trapp lipids are fresh weight |
| `cc4c54c` | loose-end audit: `conc`-endpoint bug + four stale/incorrect prose items |

**New files**: `src/neutral_dpu.py`, `validation/neutral_dpu_validation.py`,
`tests/test_neutral_dpu.py`, `tests/test_twopool_nstem_merge.py`,
`validation/twopool_nstem_merge.py`, `data_obs/neutral_obs_{liu2023,ge2017,template}.csv`,
`docs/neutral_dpu_validation.md`.

Full suite: **217 collected, 215 pass, 2 skip** (HYDRUS-engine skips).

---

## 2. Numbers to trust

| result | metric |
|---|---|
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

---

## 3. Next tasks (in-silico, do these in order)

### A1. Air exchange — volatilisation + gaseous uptake  ← start here

**The largest structural gap in the neutral base, and fully specified in-repo.** The
core ODE has **no air terms at all** (verify: `grep -niE "volatil|K_AW|P_air|Q_VOL"
src/pfas_rice_plant_module_4pool_surf.py` returns nothing). Today
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

### A3. Hwang 2017 — a second per-organ test, with a growth + exposure model

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

### A2. Expose the neutral path through `model_api`

`neutral_dpu` is currently standalone — every other capability in the repo is reachable
through `model_api`, so the app and other validation cannot use it. Small, low-risk,
and it makes the work usable. Mirror the `simulate_nstem_leaf` / `simulate_twopool_seq`
pattern: same result-dict shape, `drivers=` support, no change to any default. Whether
to add a Streamlit tab is a separate call — the app is currently PFAS-only and the
Simple/Expert split would need a decision about where a neutral compound belongs.

---

## 4. Blocked on data or experiment (not code)

**Fetchable literature** (ranked; full rationale in the session transcript and
`docs/neutral_dpu_validation.md` §5):

1. **Rice organ total lipid, fresh weight** — currently Trapp's *soybean* values.
   Scales `K_PW` linearly, so it directly conditions the 0.281. Cheapest real gain.
   Must come from a source that states its basis unambiguously.
2. **A neutral compound measured in rice GRAIN under root-only exposure** — the one
   completely untested compartment, and the food-safety endpoint. Look for submerged /
   nursery-box application trials reporting brown rice + straw; **foliar application is
   useless here** (tests deposition, not root uptake).
3. **Schriever 2020 SI** — the 97 per-compound TSCF values, to test the TSCF QSPR
   directly rather than through the plant model.
4. **Kodešová et al. 2019 ESPR** `10.1007/s11356-019-04333-9` *(DOI unverified)* — the
   measured concentrations behind Brunetti 2021, to settle the open `K_RW` = 13.3 vs
   Briggs `K_PW` ≈ 1.0 disagreement.
5. **McFarlane, Pfleeger & Fletcher 1987**, *J. Environ. Qual.* 16(4):372–376
   *(unverified)* — the bromacil measurements behind Trapp 1994, which are figure-only
   in Trapp.

**Wet lab** (cannot be fetched): in-planta half-lives for the Ge compounds (the model
predicts ≈7 d — a falsifiable prediction); rice-root cell-wall / Fe-Mn-plaque
batch-sorption + desorption assay across chain length × head group (**the `k_seq`
promotion gate**); per-congener xylem-sap / root-water ratio (direct `f_xy`); a direct
`K_cw` measurement (GAP A).

---

## 5. Housekeeping

- **PR #56 is still a draft.** CI green, mergeable, no review threads. Marking it ready
  for review is the user's call — do not flip it unasked.
- Test counts in CLAUDE.md were refreshed to 217 in `cc4c54c`; re-check if you add tests.
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
python validation/neutral_dpu_validation.py                                   # structural checks
python validation/neutral_dpu_validation.py --obs data_obs/neutral_obs_liu2023.csv   # 0.281
python validation/neutral_dpu_validation.py --obs data_obs/neutral_obs_ge2017.csv    # 0.783
python src/neutral_dpu.py                                                     # per-compound demo

# the merge
python validation/twopool_nstem_merge.py --cached      # reuse the fit (~1 min)
python validation/twopool_nstem_merge.py               # full re-fit (~40 min)

pytest tests/test_neutral_dpu.py tests/test_twopool_nstem_merge.py -q
```

**Resume prompt (example):**

> Implement task A1 from `docs/HANDOFF_neutral_next.md`: air exchange (volatilisation
> + gaseous uptake) for the neutral path, following the equations in
> `docs/dpu_model_summary_corrected.tex` §`sec:permeability` (eqs. `Pc`…`Qgas`). Keep
> it identically zero at `K_AW = 0` so every PFAS number and `reproduce_demo` (RMSE
> 0.029) is untouched; mind the SI-unit convention and the `1/(1−φ)` singularity.
> Replace the `k_aw_warning` refusal with a real run, and add tests that a volatile
> compound loses leaf mass while a non-volatile one is bit-identical to today.
