# HANDOFF — neutral-organic path: what is done, what is next

> Session handoff for the next Claude/dev.
> PRs [#56](https://github.com/ccy5123/pfas-rice-model/pull/56) and
> [#57](https://github.com/ccy5123/pfas-rice-model/pull/57) are **merged**; `main` is at
> `b3b5004`. Open: **[#58](https://github.com/ccy5123/pfas-rice-model/pull/58)** (draft) —
> the Expert neutral tab, the acquisition queue, and now the work below, which is stacked
> on the same branch.
> Scientific records: **`docs/neutral_dpu_validation.md`** (neutral path — §4d Li 2019
> partition, §4e TSCF, §4f Kodešová, §3b air, §4c Hwang) and **`docs/twopool_root_exploration.md`
> §Result 8** (structural merge). Read those first.
> `parameters.json`, `simulate()` and `reproduce_demo` (RMSE 0.029) are **UNCHANGED**
> throughout — everything here is additive or opt-in.

---

## 0. TL;DR

**The papers arrived** (`DPU4OC_add.zip`), and the previous handoff's blocking item is
cleared. `docs/literature_db/Acquisition_Queue.csv` is now a **record** of what each row
actually contained, with a per-row status, rather than a want-list.

The headline is not the two new datasets. It is what they exposed:

- **Both a-priori inputs were measured directly for the first time**, and both look
  biased low on the tables that spread widest in lipophilicity: the partition by a
  species-independent offset (§4d), TSCF by −0.221 on a 0–1 scale (§4e). But see the
  A4 result below — the partition half of that does **not** hold at moderate Kow.
- **Part of the partition offset is the model's own doing.** `neutral_dpu` anchors on
  Briggs' RCF (`L·a = 0.0302`) but the rice compartment runs `0.0122` — **2.48× under** —
  because it substituted a measured lipid and kept the conventional octanol correction,
  when only their *product* is identifiable. On Briggs' own barley rows that costs log10
  RMSE **0.266 against the anchor's 0.111**.
- **The default was not changed, and that is the decision, not an omission.** Restoring
  the anchor helps Li 2019 and hurts both soil-grown tables (Liu, Kodešová), and the two
  sides disagree by more than the change is worth. Fitting an intermediate value would
  make the path's single real claim — that nothing is fitted — false.

**Three queue rows did not deliver what they asked for**, and that is recorded rather
than papered over: A1 contains no rice, the supplied "McFarlane 1987" is the wrong paper
entirely, and every C1/C2/C3 candidate screens out for a stated reason.

**The A4 SI then arrived mid-session, and an adversarial re-read followed.**
Kodešová 2019 gives the cleanest exposure in the repo and **votes the other way**
(0.191 shipped vs 0.216 anchored); it also **supersedes the Brunetti sighting**,
measuring 1.10 for the compound Brunetti calibrated 13.3 for. The re-read
(`666459f`) then **retracted two of this session's own claims** — see §3 item 1 —
and what is left is that **the datasets contradict each other** at log Kow
3.5–4.5 by more than the anchor is worth. The anchor question is therefore **not
decidable from the data in hand**; §4f and §5 of the validation doc carry it.

---

## 1. What this session delivered

| commit | what |
|---|---|
| `1019461` | `data_obs/neutral_obs_li2019_rcf.csv` (48 rows) + `data_obs/tscf_obs_schriever2020.csv` (97 rows), `validation/li2019_rcf_apriori.py`, `validation/schriever2020_tscf.py`, the `subset` filter, `BRIGGS_ANCHORED_LIPID_FW`, `tests/test_li2019_schriever_tables.py` (11) |
| `d4dd996` | the queue rewritten as a delivery record; `docs/neutral_dpu_validation.md` §4d, §4e and the four-sightings synthesis |
| `2de5619` | section 3b — the kinetic explanation ruled out; stale test counts refreshed |
| *(this commit)* | queue A4 closed: `data_obs/neutral_obs_kodesova2019.csv`, `validation/kodesova2019_carbamazepine.py`, §4f, the sightings synthesis rewritten |

**Prior sessions on this branch**: `c6e9d8a` (Expert neutral tab), `494acc4` (the queue),
`9216dd1` (the previous handoff). Merged earlier: `d8e7f9a` air exchange, `97dbe75`
Hwang, `50b5586` `model_api.simulate_neutral`.

**⚠️ CI does not test any of this.** `.github/workflows/rigor.yml` runs only
`tests/test_sci_adk_rigor.py`, so a green check says nothing about the model. Run the
full suite locally before claiming green: **262 collected, 260 pass, 2 skip** (~10 min)
on this branch; the 2 skips are the optional `emcee` and `sci-adk` deps.

---

## 2. Numbers to trust

| result | metric |
|---|---|
| **Li 2019 root partition**, 29 out-of-sample rows, 11 species, nothing fitted | log10 RMSE **0.598**; every species biased low, −0.30 → −0.95 |
| the same, held-out Briggs barley rows (in-sample, for contrast) | 0.266 shipped / **0.111** at the anchor |
| **TSCF direct**, 30 un-ionised rows, Briggs bell | RMSE **0.310**, bias **−0.221**, Spearman +0.794 |
| the same, `schriever_tscf` — **in-sample, its own training set** | 0.234, bias −0.094 |
| logD vs logP on the full 97 (Schriever's own claim, checked) | rank corr +0.313 → **+0.653** |
| root-lipid scan, full ODE — Li 2019 / Liu 2023 / Ge 2017 | L=0.010 **0.598 / 0.281 / 0.783**; L=0.0247 (anchor) 0.331 / **0.288** / 0.651 |
| Liu 2023 root partition (unchanged) | 0.281 |
| Ge 2017 per-organ TF (unchanged) | 0.783; minimum 0.210 at a 7-day half-life |
| Hwang 2017, both basis readings | 0.610 fw / 0.726 dw — **a diagnosis, not scores** |
| **Kodešová 2019 carbamazepine**, n=21, 4 plants × 3 soils, nothing fitted | log10 RMSE 0.191 on the ODE basis, **0.237** on the appropriate equilibrium basis (§4g) |
| the repo's best a-priori root result, appropriate basis | **Liu 2023, 0.206** (not Kodešová) |
| the same, driver-free (equilibrium `K_PW`) | measured `RCF_fw` median **1.10** vs 1.56 shipped / 2.53 anchored |
| Brunetti's pea `K_RW` vs that measurement | 13.3 vs **1.10** — ~12× above, so **superseded** |
| merge, Tang per-organ OOS | 0.801 (stalk 0.61, leaf 0.28) |

**Caveats that must travel with every one of these.** Li 2019 is out-of-sample but the
model's root is parameterised for *rice* and only 4 of the 29 rows are rice; the rest is
a deliberate test of Briggs' own claim that root lipophilic character does not vary much
across species, and Li's own crop table (0.10 %–1.14 %) says it does. Exposure times run
6 h to 12 d — that confounder was pre-registered in the file and has now been **checked
and ruled out** (§3 item 2), so it is no longer a live caveat. The TSCF table is 16 species, none
of them rice, and barley is half of it; the Briggs bell is out-of-sample there but **not
independent** (same lab, same species, same method lineage). The neutral **grain
compartment remains entirely untested** — that has not changed. Kodešová is four
leafy/root vegetables at a single log Kow of 2.25, so it says nothing about the
high-Kow end where §4d's deficit actually lives; its one pivotal assumption (the
Freundlich unit reading) is defended on the implied `Koc` and flagged in the data file.

---

## 3. Next tasks (in-silico)

Ranked. The first is the only one that is both significant and unblocked.

1. **The anchor question is NOT currently decidable, and that is the finding.**
   Evidence: `docs/neutral_dpu_validation.md` §4d–§4f; reproduce with
   `python validation/li2019_rcf_apriori.py` (§3c and §3d especially) and
   `python validation/kodesova2019_carbamazepine.py`.

   Two props this handoff previously rested on were **retracted** in `666459f`,
   and a later session must not reinstate them:
   - the Li 2019 bias is *not* reliably monotone in log Kow — Namiki 2015 alone
     supplies 10 of the 29 rows, all in the top two bins;
   - raising `L` is **not** "the wrong instrument" for a Kow-dependent deficit.
     `K_PW = W + L·a·Kow^b` is water-floor-dominated at low Kow, so scaling `L`
     is inherently Kow-dependent. Li 2019 is genuine evidence **for** the anchor.

   What blocks the decision is that **the datasets contradict each other**. At
   log Kow 3.5–4.5, Li 2019 says the model is ~3× low (−0.462, n=11) and Liu 2023
   — rice, the same hydroponic endpoint — says it is exact (−0.008, n=5). On
   propiconazole, the one compound both measured, RCF is 43.65 (lettuce) vs 9.32
   (rice): **4.7× apart, against an anchor worth 2.4×**. So the question is not
   "is the partition too low" but **"which hydroponic dataset describes a rice
   root above log Kow 3.5"**, and only new data settles it (§4, item 1).

   Note also that `a = 1.22` carries **no citation anywhere in the repo**. Only
   the product `L·a` is identifiable, so "raise `L` to 0.0247" and "raise `a` to
   3.02" are the *same model* — which means "don't fit `L`, it is measured" is
   not by itself an argument against the anchor.

   Three outcomes remain defensible and the choice is the user's:
   - *keep `L = 0.01`* (current) — defensible as "do not move without evidence",
     supported mainly by Liu, which is rice and is exact where the anchor would
     do most damage. Cost: a documented 2.48× inconsistency with the anchor.
   - *restore the anchor* via `ND.BRIGGS_ANCHORED_LIPID_FW` — internally
     consistent, and Li 2019 supports it; worse on Liu and Kodešová.
   - *add the missing term* — a measured neutral-organic cell-wall coefficient
     (§4). It would fix the level independently of either dataset, but it is no
     longer needed to *explain* the Li 2019 shape.
   Do **not** quietly pick an intermediate fitted `L`: that destroys the path's
   one real claim, that nothing is fitted.
2. ~~Check the equilibration confounder on the Li table~~ ✅ **DONE** — section 3b of
   `validation/li2019_rcf_apriori.py`. The bias is **flat in exposure time** (−0.442 at
   1–3 d vs −0.447 at > 3 d; inside the high-Kow cell alone, −0.586 vs −0.517) and
   **monotone in log Kow** (−0.03 → −0.31 → −0.46 → −0.69). Non-equilibrium is ruled
   out as the driver, and the deficit is in the lipophilic sorption term — it vanishes
   at low Kow where the water floor dominates. Restoring the anchor is worth +0.394 log
   against an observed −0.688 in the worst cell, so it accounts for ~57 % and the rest
   is not lipid. Pinned by a test.
3. ~~The ODE-vs-equilibrium scoring artifact~~ ✅ **DONE** —
   `compare_to_obs(mode="equilibrium")`, §4g of the validation doc. Root tables
   measure an equilibrium but were scored through the 120-day rice season, which
   discounts the root by a Kow-dependent 0.01–0.26 log. Rescored: Liu
   0.281→**0.206**, Li 2019 0.598→**0.541**, Kodešová 0.191→**0.237**. Kodešová
   was flattered and is *not* the repo's best a-priori result; Liu is. The anchor
   verdicts survive with **wider** margins. Default stays `"ode"`; which basis
   should headline is a one-line decision, open with the anchor question.
4. **Mine Briggs 1983** (`10.1002/ps.2780140506`), the shoot-distribution companion,
   already in hand and never opened. It is the only in-hand source that could constrain
   the stem/leaf split independently of Ge 2017.
5. **Tissue specific surface areas** for rice — still the single input bounding the air
   term, still unsourced. Note the C2 papers supplied this round do **not** help (§4).
6. `NStemLeafModel` has no air hook; particle deposition (`eq:Qdep`) is still
   deliberately unimplemented. Both unchanged from the last handoff, both low value.

---

## 4. Blocked on data or experiment (not code)

**Read `docs/literature_db/Acquisition_Queue.csv` first** — it now carries a `status`
column per row and, for the rows that arrived, what they actually contained.

**A4 is closed** — the SI arrived and is worked up in §4f. What remains is
lower-value, and the two search specs are now the binding gaps.

*Previous request, for the record.* *The SI of Kodešová et al. 2019,*
`10.1007/s11356-019-04333-9` *(Tables S2 and S5).* The article arrived and gives the
chemical side — carbamazepine log Kow 2.25, pKa 1.0/13.9 so un-ionised throughout, and
Freundlich `KF` for three soils, which is the soil→pore-water conversion the exposure
needs; the measured root and leaf concentrations were SI-only. **It delivered exactly
what was hoped for, and then contradicted the expectation**: a-priori log10 RMSE 0.191,
the best in the repo, but voting *against* the anchor and superseding rather than
confirming the Brunetti sighting.

Still outstanding, lower value: the real McFarlane, Pfleeger & Fletcher 1987 (J. Environ.
Qual. 16(4):372–376 — the file supplied under that name is a different paper); the two
2025 long-chain PFAS papers (B1/B2); and the C1–C3 search specs, whose candidates this
round all screened out for reasons now recorded in the queue.

**Wet lab, ranked by what it would settle.**

1. **A neutral-organic cell-wall / non-lipid sorption coefficient for rice root.** New
   top item. Four independent sightings now say the root partition is too low (§5 of the
   validation doc), the anchor accounts for ~2.5× of it, and the residue is consistent
   with a sorbing phase the neutral composition sets to zero. The PFAS side needs the
   same measurement — it is GAP A there — so **one experiment serves both paths**.
2. Rice **root total lipid**, basis stated, measured the way Li/Chiou's crop values were
   (queue C3). This is what A1 was supposed to supply.
3. In-planta half-lives for the Ge compounds — the model predicts ≈ 7 d, a falsifiable
   prediction.
4. Unchanged from before: the `k_seq` promotion gate (rice-root cell-wall / Fe–Mn-plaque
   batch sorption + desorption across chain length × head group — note this overlaps
   item 1 and could be one campaign); per-congener xylem-sap / root-water ratio.

---

## 5. Housekeeping

- **#58 is an open draft** carrying the Expert neutral tab, the queue and this session's
  work. Marking it ready or merging is the user's call — do not flip it unasked.
- **CI only runs `tests/test_sci_adk_rigor.py`.** A green check does not mean a change is
  tested. Run `pytest -q` locally (~15 min).
- **`subset` is a load-bearing convention now.** `data_obs/*.csv` may carry a `subset`
  column; `compare_to_obs` scores `subset="apriori"` by default and files without the
  column are untouched. That inertness is what keeps Liu 0.281 / Ge 0.783 bit-identical,
  and it is pinned by a test. Do not "simplify" it away.
- **Do not silently upgrade a `doi_status`.** Several rows moved `search-only → verified`
  this session because the article itself arrived and was read. That is the only
  legitimate route.
- The papers live in an ephemeral session scratchpad and are **not committed** (copyright).
  Everything extracted from them is in `data_obs/` or the queue, with provenance.

---

## 6. How to resume

```bash
pip install -r requirements.txt

# the two new tests, and the diagnosis
python validation/li2019_rcf_apriori.py --fast    # a-priori 0.598 + the anchor table
python validation/schriever2020_tscf.py           # TSCF alone: 0.310, bias -0.221
python validation/neutral_dpu_validation.py --obs data_obs/neutral_obs_li2019_rcf.csv

# the unchanged baselines — these two must not move
python validation/neutral_dpu_validation.py --obs data_obs/neutral_obs_liu2023.csv   # 0.281
python validation/neutral_dpu_validation.py --obs data_obs/neutral_obs_ge2017.csv    # 0.783
python reproduce_demo.py                                                             # 0.029

pytest -q
```

**Resume prompt — pick by what arrived.**

*If the Kodešová SI arrived*:

> The SI of `10.1007/s11356-019-04333-9` is at `<path>`. Build a per-organ neutral obs
> table from Tables S2/S5 (carbamazepine only — atenolol and sulfamethoxazole are
> ionisable and belong on the ionic path), converting soil concentration to pore water
> with the paper's own Freundlich `KF` per soil, and record the weight basis the SI
> states. Then run it as an a-priori test and report it against Liu 0.281 / Li 0.598.
> It bears on the four-sightings question in `docs/neutral_dpu_validation.md` §5, so say
> explicitly whether it is a fifth sighting or a counter-example.

*If nothing new arrived*, take §3 item 1 to the user as a decision, then §3 items 2–3,
which are both small and both unblocked. Do not manufacture work on §3 items 4–5.
