# HANDOFF — neutral-organic path: what is done, what is next

> Session handoff for the next Claude/dev.
> **PRs [#56](https://github.com/ccy5123/pfas-rice-model/pull/56),
> [#57](https://github.com/ccy5123/pfas-rice-model/pull/57),
> [#58](https://github.com/ccy5123/pfas-rice-model/pull/58) and
> [#59](https://github.com/ccy5123/pfas-rice-model/pull/59) are MERGED**;
> [#60](https://github.com/ccy5123/pfas-rice-model/pull/60) (the `lipid_source` mode +
> the §5/§6 consistency fix, closing §3 items 1–2) is the only thing still open.
> Scientific record: **`docs/neutral_dpu_validation.md`** — §4a Liu, §4b Ge, §4c Hwang,
> §4d Li 2019 hydroponic, §4e TSCF, §4f Kodešová, §4g the scoring artifact, §4h Li 2019
> soil, §4i Kodešová's leaf, §4j Briggs 1983 stem, §5 the synthesis. **Read §5 first.**
> `parameters.json`, `simulate()` and `reproduce_demo` (RMSE 0.029) are **UNCHANGED**
> throughout — everything on this arc is additive or opt-in, and no PFAS number moved.
> Full suite on the merged tree: **278 collected, 276 pass, 2 skip** (~12 min); the two
> skips are the optional `emcee` and `sci-adk` deps.

---

## 0. TL;DR

The acquisition queue's papers arrived and have been **mined out**. Four measured
tables now ship, and the question the last handoff was built around has **dissolved
rather than been answered**.

**It is not "the root partition is too low". It is "the exposure term is the weak
part".** Where the exposure is measured or directly known, the model is close to
unbiased on three independent tables; the large deviations sit where the exposure was
*estimated*:

| dataset | n | exposure | bias |
|---|---|---|---|
| Liu 2023 (rice) | 14 | hydroponic, **known** | −0.053 |
| Li 2019 soil, `K_om` **experimental** | 62 | soil, **measured** | **+0.033** |
| Kodešová 2019 | 21 | soil, **measured isotherm** | +0.162 |
| Li 2019 soil, `K_om` from a QSPR | 261 | soil, **estimated** | +0.291 |
| Li 2019 hydroponic | 29 | hydroponic, known | **−0.432** ← the open anomaly |

**The anchor decision is now CLOSED — and it closed as a mode, not a value.** Li 2019's
hydroponic half was the only thing arguing for restoring the Briggs anchor; its own
**soil half, twelve times larger, argues the opposite**. Three tables against, one for.
The default stays **`lipid_source="measured"` (`L` = 0.01)**, recorded with the reason
*"no evidence to move"*, **not** *"0.01 is validated"* — and the alternative is kept
runnable as a named mode rather than deleted, because the question is open, not
answered. Reproduce the whole trade-off with one command:

```bash
python validation/neutral_dpu_validation.py --lipid-source both   # add --mode equilibrium
```

---

## 1. What this arc delivered

| commit | what |
|---|---|
| `1019461` | `data_obs/neutral_obs_li2019_rcf.csv` (48) + `tscf_obs_schriever2020.csv` (97), `li2019_rcf_apriori.py`, `schriever2020_tscf.py`, the `subset` filter, `BRIGGS_ANCHORED_LIPID_FW` |
| `d4dd996` | the acquisition queue rewritten as a delivery record; docs §4d, §4e |
| `2de5619` | §3b — the kinetic explanation for the Li offset ruled out |
| `78d5794` | queue A4 closed: `neutral_obs_kodesova2019.csv` (21), `kodesova2019_carbamazepine.py`, §4f |
| `666459f` | **two of this arc's own claims retracted** — see §2 |
| `a3b07fa` | `compare_to_obs(mode="equilibrium")` — the ODE scoring artifact, §4g |
| `784ef4d` | `neutral_obs_li2019_soil.csv` (376) + `li2019_soil_table.py` + Kodešová's leaf half, §4h/§4i |
| `5426237` | Briggs 1983 read: `briggs_scf()`, `briggs1983_stem.py`, §4j |

Earlier on the same arc: `d8e7f9a` air exchange, `97dbe75` Hwang, `50b5586`
`model_api.simulate_neutral`, `c6e9d8a` the Expert neutral tab, `494acc4` the queue.

**Then, closing the two open decisions**, the root-lipid anchor became a named mode
`lipid_source="measured"|"briggs_anchor"` in the repo's own idiom (`f_xy_source`,
`cwo_profile`, `biomass`, `tscf_model`), threaded `rice_compartments` →
`neutral_dpu.simulate_neutral` → `model_api.simulate_neutral` → `compare_to_obs`, with
`compare_lipid_sources()` / `--lipid-source both` scoring all five shipped tables under
both readings. Nothing moved (Liu 0.281 / Ge 0.783 / `reproduce_demo` 0.029 re-verified).
The two prior anchor tests were switched off their global `TRAPP1994_LIPID_FW` mutation
onto the mode.

**And a defect found while reading in:** §5's table and the whole of §6 still carried two
of the three claims §2 lists as retracted — "Kodešová 0.191 is the best a-priori result"
and the Brunetti-based "best-evidenced open question" synthesis — even though §4f and §4g
retract both in the body. They survived in the **Honest summary**, the section most likely
to be read on its own. Now corrected, along with the mixed-basis problem in §3 item 2.
Treat that as the standing risk in this document: the body gets updated, the summary
does not.

**⚠️ CI does not test any of this.** `.github/workflows/rigor.yml` runs only
`tests/test_sci_adk_rigor.py`. Run `pytest -q` locally before claiming green.

---

## 2. Things this arc got wrong and then corrected — do not reinstate them

Three claims were made, written into the docs, and then retracted on re-reading. Each
is pinned by a test now. A later session that "rediscovers" any of them is going
backwards.

1. **"The Li 2019 bias is monotone in log Kow."** Not robust. **Namiki 2015 alone
   supplies 10 of the 29 rows, all in the top two Kow bins**, as two compounds × five
   species. Collapsing compound × study flattens the top two bins (−0.576, −0.501).
   What survives: no bias below log Kow 2, −0.3 to −0.6 above, and all ten source
   studies biased the same way.
2. **"A flat lipid rise is the wrong instrument for a Kow-dependent deficit."**
   Wrong. `K_PW = W + L·a·Kow^b` is water-floor-dominated at low Kow, so scaling `L`
   is *inherently* Kow-dependent (+0.045 log at log Kow 1, +0.391 at 5). Li 2019
   hydroponic is genuine evidence **for** the anchor, not against it.
3. **"Kodešová 0.191 is the best a-priori result in the repo."** It was flattered by
   a scoring artifact (§4g): root tables measure an equilibrium but were scored
   through the 120-day rice season, which discounts the root by a Kow-dependent
   factor peaking at 0.55 near log Kow 1.78 — right where Kodešová sits. On the
   appropriate basis it is **0.237**, and **Liu 2023 at 0.206** is the best.

Also worth carrying: **`a = 1.22` has no citation anywhere in the repo.** Only the
product `L·a` is identifiable, so "raise `L` to 0.0247" and "raise `a` to 3.02" are
the same model — which means "don't fit the measured `L`" is *not* an argument
against the anchor.

---

## 3. Next tasks

Ranked. **Items 1 and 2 of the previous handoff are DONE** — see §1.

1. ~~Put the anchor decision to rest~~ **DONE.** Default stays `lipid_source="measured"`;
   the alternative is a named mode, and the 3-vs-1 evidence is a command rather than a
   doc claim. Recorded in `docs/neutral_dpu_validation.md` §5.
2. ~~Decide which basis headlines the neutral path~~ **DONE, as "quote both".**
   `mode="ode"` stays the default so no published number moved, and every quoted number
   now names its basis. This also caught a real inconsistency: **the document was already
   mixing bases across sections** — §4a/§4d/§4f are ODE, §4h's soil numbers are
   equilibrium — so the same table read 0.549→0.670 or 0.639→0.873 depending on a mode
   nobody was stating. Both are now printed by the same command.
3. **Chase the Li 2019 hydroponic anomaly** (§5). It is the one open disagreement and
   no subgroup explains it: aquatic −0.39 vs terrestrial −0.46, organochlorines −0.51
   vs everything else −0.35, the four rice rows −0.28, all ten studies negative. On
   propiconazole — the one compound it and Liu both measured at log Kow 3.72 — they
   report RCF **43.65** (lettuce) vs **9.32** (rice). Hypothesis on the table:
   root-surface sorption inflating short hydroponic RCFs, which would bite hardest
   for hydrophobic compounds and is consistent with the −0.51/−0.35 split. Testing it
   needs the source studies' washing protocols, i.e. new papers (§4).
4. **Mine what is left of Briggs 1983** — its per-section shoot concentrations (stem
   base / central / leaf blade) are the only in-hand data that could constrain the
   stem/leaf split independently of Ge 2017. §4j used only its fitted equations.
5. **Tissue specific surface areas** for rice — still the one input bounding the air
   term, still unsourced. The C2 papers supplied this round do **not** help (they are
   root morphology; the air term takes no root contribution).
6. `NStemLeafModel` has no air hook; particle deposition (`eq:Qdep`) is deliberately
   unimplemented. Both low value.

---

## 4. Blocked on data or experiment

`docs/literature_db/Acquisition_Queue.csv` carries per-row status. **A2, A3, A4 are
closed; A1 is PARTIAL** (it had no rice); **A5, B1, B2, C1, C2, C3 remain open**.

**The one request that would move something.** *A hydroponic rice root RCF above
log Kow 3.5.* That is the exact cell no table has — Liu tops out at 4.4 with only five
rows there, Li 2019's rice rows are all below 3.1, Kodešová is a single point at 2.25 —
and it is what would settle item 3. Note this is a **different** request from the
queue's C3 ("rice root total lipid"), which is now the lower priority of the two.

Also still wanted, unchanged: the real McFarlane, Pfleeger & Fletcher 1987 (the file
supplied under that name is a different paper); the two 2025 long-chain PFAS papers
(B1/B2); a neutral compound in rice **grain** under root-only exposure (C1 — the one
entirely untested compartment; every candidate this round had residues below LOD).

**Wet lab**, re-ranked by what the new results imply:

1. **A measured pore-water concentration alongside tissue** for any neutral compound
   in rice. §5 now locates most of the apparent partition error in the exposure term,
   so this beats any further partition work.
2. In-planta half-lives — but note §4i: they are **species-dependent** (carbamazepine
   leaf parent fraction 0.17 in lamb's lettuce, 0.81 in radish, same compound, same
   soils, same harvest), so a single fitted value is not a compound property.
3. A neutral-organic cell-wall coefficient (PFAS **GAP A**, serves both paths). It has
   dropped in priority: it is no longer needed to explain any table.
4. Unchanged: the `k_seq` promotion gate; per-congener xylem-sap / root-water ratio.

---

## 5. Housekeeping

- **CI only runs `tests/test_sci_adk_rigor.py`.** A green check means nothing here.
- **`subset` and `mode` are load-bearing conventions.** `compare_to_obs` scores
  `subset="apriori"` and `mode="ode"` by default; files without a `subset` column are
  untouched, which is what keeps Liu 0.281 / Ge 0.783 bit-identical. Both are pinned
  by tests. Do not "simplify" either away.
- **Do not silently upgrade a `doi_status`.** `verified` means the article itself was
  opened and read.
- The papers live in an ephemeral scratchpad and are **not committed** (copyright).
  Everything extracted is in `data_obs/` or the queue, with provenance.
- Two assumptions are flagged in the data files rather than buried, and a later
  session should attack these first if it doubts a result: Kodešová's **Freundlich
  unit reading** (defended on the implied `K_oc` 222/189/154 across three soils; the
  alternative gives ~1600 and would reverse that table's vote), and that Kodešová's
  **mass balance does not close** (measured soil holds 25–33 % of the applied load).

---

## 6. How to resume

```bash
pip install -r requirements.txt

# the four new tables and what they show
python validation/li2019_soil_table.py            # 376 soil rows: the sign flip, ~5 s
python validation/li2019_rcf_apriori.py --fast    # 29 hydroponic rows + the anchor
python validation/kodesova2019_carbamazepine.py   # exposure, anchor vote, leaf, metabolism
python validation/schriever2020_tscf.py           # TSCF alone, 97 values
python validation/briggs1983_stem.py              # the stem anchor, ~1 s

# the closed anchor decision, as a command rather than a claim (~6 min)
python validation/neutral_dpu_validation.py --lipid-source both
python validation/neutral_dpu_validation.py --lipid-source both --mode equilibrium

# the baselines that must not move
python validation/neutral_dpu_validation.py --obs data_obs/neutral_obs_liu2023.csv  # 0.281
python validation/neutral_dpu_validation.py --obs data_obs/neutral_obs_ge2017.csv   # 0.783
python reproduce_demo.py                                                            # 0.029

pytest -q                                          # 278 collected, 276 pass, 2 skip
```

**Resume prompt.**

> The neutral-organic arc is merged, its papers are mined out, and **both open
> decisions are now closed** (§3 items 1–2): the root lipid stays at `L = 0.01` as the
> default of a named `lipid_source` mode, and the docs quote both scoring bases with
> every number labelled. Read `docs/neutral_dpu_validation.md` §5 and §2 of
> `docs/HANDOFF_neutral_next.md` before touching anything — §2 lists three claims this
> arc made and then retracted, and re-deriving any of them is going backwards. Note the
> failure mode that produced two of them: the body of the validation doc gets corrected
> and the **Honest summary does not**, so check §6 against §4 whenever you change a
> result. Remaining work, honestly ranked: §3 item 3 (the Li 2019 hydroponic anomaly —
> blocked, needs a hydroponic rice root RCF above log Kow 3.5) then item 4 (Briggs
> 1983's per-section shoot data, already in hand and the only unmined thing here). The
> PFAS side is a separate arc — `docs/HANDOFF_BAF_twopool.md`.
