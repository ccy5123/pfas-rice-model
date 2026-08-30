# HANDOFF — the carrier-vs-bypass question: what was settled, and the two jobs left

> Session handoff. **[PR #62](https://github.com/ccy5123/pfas-rice-model/pull/62)**
> carries the pre-registration (`fb0abed`) and the results (`9217314`) as separate
> commits, in that order, on purpose.
> Record: **`validation/carrier_vs_bypass.py`** — its header holds the
> pre-registration and its VERDICT block the outcome. Read the header first; the
> three pre-registered items are what the result is judged against.
> `params/parameters.json`, `simulate()` defaults and `reproduce_demo` (0.029) are
> **UNCHANGED**.
>
> **The next session is to do A1 then B** (see §3), decided with the user.

---

## 0. TL;DR

**The question.** `docs/theory_anchor.tex` states that the four-compartment DPU
(Rein 2011, Brunetti 2019) is, *at the membrane level*, the **Trapp (2000, 2004)
ionizable-compound cell model**, and marks only three things `— (new)`: `f_xy`,
`η`, and the **Michaelis–Menten carrier**. `f_xy` is a lumped stand-in for
Trapp's own detailed root model and `B_k` is his `K_RW` re-parameterised — so the
**carrier is the one piece of physics this repo adds**, it was fitted
(*"fixed during W2 fit"*), and it had never been compared with an alternative.

**Three arms, one global parameter each**, on Yamazaki root/straw/grain with the
a-priori monotone `f_xy` held fixed so only the uptake term differs:

| | log10 RMSE |
|---|---|
| nothing (Trapp's PFAS limit as-is) | 2.640 |
| **A carrier** (incumbent, `Vmax` fitted on this data) | **1.035** |
| **B bypass** (one global `g_apo` = 20) | **0.996** |
| **C depolarisation** (`E_m` = −90 mV, **no new term**) | **2.289** |

**Three findings, in order of how much they should change what you do.**

1. **An addition IS necessary — arm C is REFUTED.** Depolarising the membrane to
   the far end of its own recorded plausible range (`e^N` 107 → 33) buys almost
   nothing and stays far worse than either addition. The lever already inside
   Trapp's framework cannot carry PFAS, so **the repo's extension of it is
   justified in existence.** That was the strongest form of the question.
2. **WHICH addition is undecided.** A vs B is 0.04 log units over 33
   observations; bootstrapped, **P(bypass beats carrier) = 0.749** — not a
   difference. And B cannot claim parsimony: `Vmax_in` is a **single global
   20.0**, not per-congener (checked, not assumed — an earlier framing of this
   test got that wrong).
3. **The bypass's own distinguishing claim FAILS — pre-registered item 1 is
   REFUTED.** Fitted per congener, `g_apo` trends strongly with chain length
   (**corr +0.832, spread 25×**), where `theory_anchor.tex` says η — which
   contains the apoplastic bypass ε — is *"essentially independent of tail
   length"*. On the criterion written before the run, that is a relabelled
   carrier, not a water-flow bypass. **This is the first time the theory doc's
   η claim has been contradicted by data.**

**And the convergence is the real finding.** LC6 found the per-congener `Vmax`
multiplier needed is ~flat to C10 then 2.0× at C11 and 5.5× at C12; per-congener
`g_apo` here runs 2–5 for C4–C8 then 20–50 for C9–C12. **Two different entry
terms need the same chain-length correction, of similar size** — and a
requirement common to both is not a property of entry. It belongs to
`B_k`/`φ_free` (sequestration), exactly where the two-pool work put it.
Pre-registered item 3 says it from the other side and was **CONFIRMED**: A and B
miss the long chain about equally (1.640 vs 1.439), so no winner was read out of
that end.

---

## 1. What this delivered

| commit | what |
|---|---|
| `fb0abed` | **pre-registration**, pushed before the run existed: three arms, three failure modes, `model_api.simulate(vmax_scale=, g_apo=)` in the existing override idiom (defaults bit-identical) |
| `9217314` | **results**, judged against those three items |

Upstream on the same thread: `000a4e0` added `g_apo` itself (the mechanism
§4l pointed at), and `dfe8106` is the weak-electrolyte test that motivated it.

---

## 2. Things to not re-derive

1. **"The bypass has fewer free parameters than the carrier."** False. `Vmax_in`
   is global (20.0), not per-congener. Both arms are one parameter. This was
   asserted before being checked and is wrong.
2. **"η is essentially independent of tail length" (`theory_anchor.tex`).** Not
   supported: the per-congener `g_apo` fit correlates with `n_C` at **+0.832**.
   The doc is a derivation and has not been edited — see §3 B3.
3. **"B beats A."** 0.996 vs 1.035 is not a result (bootstrap 0.749). Quote them
   as indistinguishable or not at all.
4. **"The dose series refutes the carrier."** Over-read. The pre-registered rule
   was **not met** (1 of 3). PFOS refutes; PFOA is too noisy; GenX cannot speak —
   it is saturated at every dose, so its "not below" is not evidence *for* the
   carrier either. The defensible word is **disfavoured**.
5. **The gate in `dose_series_carrier.py` was written across congeners.** It
   should have been per congener; running it exposed that. GenX passes the global
   gate on PFOS's back while being individually non-informative.
6. The absolute levels here are **a-priori-limited** (monotone `f_xy`, not fit)
   and are NOT comparable to `reproduce_demo`'s fitted 0.029, nor to the
   documented a-priori ~0.84 which runs on different drivers.

---

## 3. The two jobs left, as decided with the user

### A1 — keep the carrier as default, add the bypass as a named mode — **DONE**

Shipped as `model_api.UPTAKE_MODES` / `simulate(uptake="carrier"|"bypass")`,
default `"carrier"` and bit-identical to the shipped solve (pinned by
`tests/test_carrier_vs_bypass.py::test_uptake_mode_default_is_the_shipped_carrier`);
explicit `vmax_scale=`/`g_apo=` still win over the mode, so the scans in
`validation/carrier_vs_bypass.py` keep working unchanged. The reason the carrier
stays default — **by default, not by evidence** — is recorded in the code beside
the modes, not only here. Original brief below.



The repo's own idiom for an unsettled question (`lipid_source`, `f_xy_source`,
`cwo_profile`, `biomass`, `tscf_model`, `mode`): the default does not move, and
the alternative is **kept runnable rather than buried in a validation script**.

- `model_api.simulate(uptake="carrier"|"bypass")`, default `"carrier"`.
- Record the reason honestly, as `lipid_source` did: the carrier keeps its place
  **by default, not by evidence** — A and B are indistinguishable on the only
  dataset that has been asked.
- `parameters.json`, defaults and `reproduce_demo` must not move.

### B — the follow-up, ranked

**B1 — DONE.** `validation/dose_series_carrier.py` (pre-registration `f705de3`
committed before the results, as asked). Outcome, in the order it should be
quoted:

- **The pre-registered rule is NOT met** — 1 of 3 congeners, where 2 were
  required. On the rule as written the carrier is **not refuted**.
- **The tally is misleading and the run says why.** Only ONE congener is
  well-conditioned for this test. **GenX** sits at `Cwo/Km` = 465 at the *lowest*
  dose, so the carrier is saturated throughout and the two arms differ by 7% —
  failure mode G1 for GenX alone. That is **a flaw in the pre-registration
  itself: the gate should have been per congener**, not across them.
- **PFOA is inconclusive**: observed 1.33× *is* below the 1.48× midpoint, but the
  bootstrap over Tang's SDs is 0.671, and it is Kd-limited at the low end.
- **PFOS is the well-conditioned case and it refutes**: the only compound whose
  exposure actually crosses `Km` inside the series. Carrier predicts a **6.28×**
  BCF decline; measured is **1.17×**; bootstrap 1.000; holds in 8 of 9 Kd
  combinations.
- **The endpoint test agrees and does not depend on `Kd` at all**: the model
  confirms entry magnitude divides out of TF (dose-invariant to 3 dp in both
  arms), so a carrier *cannot* produce a TF trend — yet PFOA's TF falls 2.1–2.3×
  where its BCF falls 1.33×. The dose response is in translocation, and is
  bigger than the entire uptake signal.
- **The durable result is a BOUND, and it is POST-HOC (labelled as such in the
  file, never quote it as a passed pre-registered test).** To be as flat as
  measured the carrier needs `Km` ≥ **500 µg/L, 100× the fitted 5** — and PFOA and
  PFOS land on the *same* bound despite ~6× different pore water. A carrier that
  linear across its whole exposure range **is** the bypass term, with `Vmax/Km`
  as its conductance. So the series bounds the carrier into the bypass's
  functional form rather than choosing between them.

Honest summary: **disfavoured, not refuted.** Nothing adopted — `Km` is not
re-fitted from three congeners in one soil.

**B3 — DONE.** The η contradiction is recorded in `docs/theory_anchor.tex` as a
paragraph after eq. (factor). The equation is **left as written on purpose**: a
chain-length requirement common to two *different* entry terms is unlikely to be
a property of entry, so the reading is that the factorisation is
**mis-partitioned** — what the fits absorb into η belongs in `φ_free`/`B_k`. No
replacement is asserted; that is B2, gated on the wet-lab assay.

**Original B1 brief, kept for the record:**

**B1 (recommended). Separate the two with a DOSE SERIES.** This is the point:
§0 finding 2 is a limitation of the *observable*, not of the question. A
saturable carrier and a linear bypass differ exactly in how they respond to
concentration, and the data are already in the repo —
`docs/literature_db/raw_si/tang2026_doseresponse.csv`: PFOA/PFOS/GenX × **5 soil
doses spanning 1000×** (0.1 → 100 µg/g), BCF and TF with SDs.

First look, **not a result**: PFOA BCF runs 0.240 / 0.225 / 0.191 / 0.209 / 0.181
across that 1000× range — about 1.3-fold. A carrier at `Km_in = 5 µg/L` must pass
from linear into saturation somewhere in that span, which should cost far more
than 1.3×. That is *suggestive* against the carrier and it is exactly why the run
is worth doing rather than assuming.

**Pre-register before running it** (this is the lesson of `g_apo`, twice over):
- The confounder is stated in CLAUDE.md already — **Tang's TF declines with dose
  because of toxicity**, not only because of saturation. A carrier and a toxic
  response both predict a decline, so the discriminator is the **magnitude and
  shape**, not the sign. Decide in advance what decline is "carrier-sized"
  (compute it from `Km_in` and the pore water implied by each dose via `Kd`)
  versus what is toxicity-sized, and say how they will be told apart.
- Say in advance what the three congeners must do **together**: a saturation
  effect should track each compound's own pore-water concentration, a toxic one
  should track dose more uniformly.
- Bootstrap any close comparison. Twice this session a gap that looked real did
  not survive resampling.

**B2. Chase the convergence.** Both entry terms need the same chain-length
correction ⇒ it lives in `B_k`/`φ_free`, not at the membrane. Test whether fixing
the sequestration term removes the need for a chain-length-dependent entry term
at all. More fundamental than B1 and longer; it entangles with the `k_seq`
promotion decision, whose gate is the §5 wet-lab assay, so it cannot be closed
in-repo.

**B3. Record the η contradiction in `theory_anchor.tex`.** Small. It is a
derivation document and was not edited here on purpose. Do it alongside B1 or B2
rather than on its own — it is bookkeeping, not new knowledge.

---

## 4. How to resume

```bash
pip install -r requirements.txt

python validation/carrier_vs_bypass.py --fast    # ~8 min; drop --fast for finer grids
python validation/dose_series_carrier.py         # ~6 min; --fast skips the Kd sweep
python reproduce_demo.py                         # 0.029, must not move
pytest -q                                        # CI runs this too now
```

**What is left after this session: B2 only** (plus the wet-lab gate it shares
with `k_seq`). A1, B1 and B3 are done. B2 is the one the two independent results
now point at from opposite sides — the bypass fit and the carrier fit need the
*same* chain-length correction, and the dose series says the entry term is not
where the concentration response lives either. Test whether fixing the
sequestration term (`B_k`/`φ_free`) removes the need for a chain-length-dependent
entry term at all. It cannot be closed in-repo: its gate is the §5 rice-root
cell-wall / Fe-Mn-plaque batch-sorption and desorption assay.


**Resume prompt.**

> The carrier-vs-bypass test is done and its record is
> `validation/carrier_vs_bypass.py` — read the header (the pre-registration)
> before the VERDICT, and §2 of `docs/HANDOFF_carrier_vs_bypass.md` for three
> things already checked that must not be re-derived. Outcome in one line: an
> addition to Trapp's ionizable cell model **is** necessary (depolarisation alone
> is refuted), but **which** addition is undecided (carrier vs bypass, bootstrap
> 0.749), and the bypass's own distinguishing claim — chain-length independence,
> asserted by `theory_anchor.tex` — was refuted at corr +0.832.
>
> **A1, B1 and B3 are DONE — do not start them again.** A1 is
> `simulate(uptake="carrier"|"bypass")`; B1 is
> `validation/dose_series_carrier.py` (read its header before its VERDICT); B3 is
> the paragraph after eq. (factor) in `theory_anchor.tex`. B1's outcome in one
> line: the pre-registered rule was **not met** (1 of 3 congeners, needed 2)
> because only ONE congener is well-conditioned — but that one, PFOS, refutes the
> carrier decisively (6.28× predicted vs 1.17× measured, bootstrap 1.000), so the
> honest verdict is **disfavoured, not refuted**, and the durable result is the
> POST-HOC bound `Km` ≥ 500 µg/L (100× the fitted value), at which the carrier is
> linear across the whole span and is therefore the bypass term in disguise.
>
> **What is left is B2**, and both of this arc's results now point at it from
> opposite sides: the two entry terms need the *same* chain-length correction
> (so it is not a property of entry), and the dose series says the concentration
> response is not in the entry term either. Test whether fixing the sequestration
> term (`B_k`/`φ_free`) removes the need for a chain-length-dependent entry term.
> It cannot be closed in-repo — its gate is the §5 wet-lab assay it shares with
> the `k_seq` promotion decision.
>
> §2 of this file lists things already checked that must not be re-derived; add
> to them that **the pre-registered gate in `dose_series_carrier.py` was written
> across congeners when it should have been per congener** — GenX is saturated at
> every dose and its "not below" is not evidence for the carrier.
>
> The neutral-organic arc is a separate handoff (`HANDOFF_neutral_next.md`); the
> PFAS two-pool arc is `HANDOFF_BAF_twopool.md`.
