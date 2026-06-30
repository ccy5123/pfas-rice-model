# Mechanistic anchor for the U-shaped, head-group-dependent root sequestration rate `k_seq`

> Literature investigation (handoff item ②, `docs/HANDOFF_BAF_twopool.md` §4.5):
> does the phenomenological two-pool `k_seq` — an **irreversible** root-sequestration
> rate that is **U-shaped in chain length** with a **PFSA head-group offset**, fit to
> data and currently with NO mechanistic justification — have a physical basis?
>
> Method: a fan-out, adversarially-verified literature synthesis (deep-research
> harness; 17 sources → 66 candidate claims → 25 verified by 3-vote panels). This
> doc distills the verified findings and maps them onto the model's `k_seq`. It is
> a **literature/theory** deliverable; `parameters.json` and the model are UNCHANGED.

---

## 0. Verdict — PARTIALLY SUPPORTED (and it explains the model's two-arm form)

The phenomenological `k_seq` is **partially anchored**, and the way it is anchored is
itself informative. The two pieces that the model needed but could not justify —
the **irreversibility** of the sink and the **PFSA offset** — are now both supported
by independent literature. The **U-shape in chain length is NOT a single mechanism**;
it is a **superposition of two distinct sorption mechanisms** acting at opposite ends
of the chain-length range. That superposition is exactly why the model's fitted
`k_seq` form is a **sum of two exponential arms**:

```
k_seq(n, group) = [ 0.268·e^(−0.52·(n−4))   +   0.615·e^(1.35·(n−12)) ] · {10^(+0.18) if PFSA}
                    \__ short-chain arm __/     \__ long-chain arm __/      \__ head-group __/
                    decays from C4               rises toward C12            sulfonate offset
```

| model element | mechanistic anchor (verified) | confidence |
|---|---|---|
| **irreversible sink** (terminal accumulator, not a reversible K) | desorption-resistant / hysteretic PFAS retention even with near-linear isotherms; TII rises with chain length, ≈0.98 (near-complete) for C10 PFDA | HIGH |
| **long-chain rising arm** `e^(1.35(n−12))` | the same hysteresis/irreversibility index increases **monotonically** with chain length (hydrophobic entrapment) | HIGH |
| **short-chain arm** `e^(−0.52(n−4))` | **biphasic** soil Koc: short chains sorb MORE than a long-chain QSPR predicts, via **electrostatic / anion-exchange** to inorganic surfaces (a DIFFERENT mechanism from the long-chain arm) | HIGH |
| **PFSA offset** `10^(+0.18)` | sulfonate head group sorbs more than the chain-matched carboxylate to lignin / cell-wall / soil: **+0.23 log** Kd (lignin) / Koc (soil), non-lipid, electrostatic | HIGH |
| (NOT the mechanism) Fe-Mn root plaque as the irreversible sink | ferrihydrite binding is **outer-sphere / electrostatic, pH-reversible**, monotonic, C6+ only, acidic-pH only → weakest candidate | HIGH |

**Quantitative corroboration of the PFSA offset.** The model's fitted PFSA offset on
`k_seq` is `10^(+0.18)` = **+0.18 log**. The independently **measured** sulfonate-vs-
carboxylate offset is **+0.23 log** (Higgins & Luthy soil Koc; lignin Kd). A
fitted +0.18 against a measured +0.23 is a genuine, not-circular agreement — the
model recovered a head-group offset of the right sign and nearly the right magnitude
from BAF data alone.

**The central caveat (also the decisive experiment).** No measured **rice-root
cell-wall** or **rice Fe-Mn-plaque** PFAS sorption coefficient resolved by **both**
chain length **and** head group exists in the surveyed literature. Every anchor above
is an **analog matrix** (lignin, ferrihydrite, sediment, bulk-soil OC). So `k_seq`
is anchored by **analogy + superposition**, not by a direct root-sorbent dataset.

---

## 1. The two arms are two mechanisms (why a single linear `k_seq` failed)

The model's history (`docs/twopool_root_exploration.md` Result 2–3) found that a
**linear** global `k_seq(n)` collapsed (`ks_b → 0`) and only a **U-shape** — a
decaying short-chain arm plus a rising long-chain arm — fit the root-matched values.
At the time this was purely empirical. The literature now gives each arm a separate
physical origin:

- **Long-chain arm (rising, C9–C12).** PFAS retention is genuinely **non-equilibrium
  and desorption-resistant** even when the sorption isotherm is near-linear: a fraction
  is entrapped and resists desorption (strong hysteresis), and the **thermodynamic
  irreversibility index (TII) increases with chain length**, reaching **≈0.98 for the
  C10 carboxylate PFDA** (urban-reservoir sediments; Chemosphere sciencedirect.com/.../pii/S0045653515302502).
  This both (i) **justifies an irreversible sink term** (a reversible partition K
  cannot represent it) and (ii) supplies the **rising long-chain arm**.
  → Hydrophobic/size-driven entrapment; head-group-agnostic.

- **Short-chain arm (the high C4–C5 limb).** Anionic-PFAS organic-carbon-normalized
  **log Koc is biphasic** in molecular size: short chains sorb **more** than the long-
  chain hydrophobic QSPR predicts, because of additional **electrostatic / anion-
  exchange** interactions with inorganic surfaces that hydrophobic partitioning does
  not capture (multi-study soil synthesis; Environments `10.3390/environments10100175`).
  → A **different** (polar/electrostatic) mechanism from the long-chain arm.

So the model's `k_seq` U-shape is a **superposition of a short-chain electrostatic
arm and a long-chain hydrophobic/hysteretic arm**, not one process with a minimum.
The adversarial panel **refuted** a specific proposed molecular split (Al-oxide for
short / OM+air-water for long; vote 0-3), so the exact partition of the two arms is
unsettled — but the **two-mechanism superposition itself is supported**.

---

## 2. The PFSA head-group offset is anchored (non-lipid, electrostatic)

The model separates PFOS (C8 sulfonate) from PFUnDA (C11 carboxylate) — identical
measured `K_PL` (31623) yet root BAF 5.93 vs 19.53 — with a **non-K_PL** `k_seq` plus
a PFSA multiplier `10^(+0.18)`. The literature supports a head-group offset on a
**non-lipid** sink, independent of pure hydrophobicity:

- **Lignin** (a rice cell-wall polymer): sulfonate PFAS (PFOS, PFBS) show higher Kd
  than carboxylate analogs; the **sulfonate moiety adds +0.23 log Kd**
  (J. Hazard. Mater. sciencedirect.com/.../pii/S0304389424025950).
- **Soil**: Higgins & Luthy quantify **+0.23 log Koc** for the sulfonate moiety beyond
  the carboxylate at equal chain length (the same offset already used elsewhere in
  this model's `f_xy` head-group term; Soil Research `10.1071/SR22183`).
- **Kaolinite (MD)**: deprotonated **sulfonate head groups coordinate more strongly**
  to the surface hydroxyl layer than carboxylate at equal chain length
  (ES&T `10.1021/acs.est.5c01046`).

This anchors the **sign and approximate magnitude** of the model's PFSA `k_seq`
offset. Caveat: the measured offset is **monotonic and small** (~0.23 log) — it
supplies the **PFSA-vs-PFCA separation**, not the U-shape; and the matrices are
lignin/clay/soil, not directly rice-root cell wall.

---

## 3. Fe-Mn root plaque is the WEAKEST candidate for the irreversible sink

It is tempting to attribute an irreversible paddy-root sink to the **Fe/Mn plaque**
(oxyhydroxides deposited on the root surface under flooded redox). The evidence argues
**against** plaque as the desorption-resistant `k_seq` sink:

- PFAS binding to **ferrihydrite** (the plaque analog) is **outer-sphere /
  electrostatic**: S K-edge XANES shows the sulfonate sulfur keeps its +V state with
  **no inner-sphere surface complexes**; sorption tracks surface charge / zeta-
  potential, is **inversely related to pH**, and is suppressed by competing phosphate
  (ES&T `10.1021/acs.est.0c01646`; PMC7745537).
- It is **monotonic and C6+-only**: at pH 4, sorption rises from 31–43% (short chains)
  to 60% (PFOA) to 100% (PFDA) — short chains sorb **less**, the opposite of the
  short-chain arm. And it needs **acidic pH** (zeta-potential > +20–25 mV); **flooded
  paddies are often circumneutral**, where ferrihydrite PFAS sorption is low.
- Outer-sphere binding is **reversible / pH-switchable**, not a desorption-resistant
  sink.

**Nuance (a plausible PFCA-only plaque component).** ATR-FTIR work indicates PFCA
**carboxylates can form inner-sphere Fe-carboxylate complexes** while PFOS sulfonate
stays outer-sphere — so a head-group-specific plaque contribution is plausible **for
PFCAs**, not for the sulfonate. And rice Fe-plaque IS a proven **anion-selective,
partially irreversible "molecular sieve"** in another system: it sequesters arsenate
desorption-resistantly (inner-sphere bidentate-binuclear) while remobilizing phosphate
(labile/outer-sphere) (ES&T `10.1021/acs.est.5c16125`). This establishes that plaque
**can** host a selective, desorption-resistant sink **in principle** — but arsenate
binds by inner-sphere ligand exchange (different chemistry from PFAS), so it is a
**precedent/analogy**, not direct PFAS evidence.

→ The irreversible sink is better attributed to **cell-wall / hydrophobic-entrapment**
(§1–§2) than to Fe-Mn plaque, though a minor PFCA-specific inner-sphere plaque term
cannot be excluded.

---

## 4. Apoplast / Casparian strip — a translocation barrier, not a quantified sink

The **apoplastic / Casparian-strip** barrier is a physically real, species-variable
control on root→shoot delivery (radish, which lacks root Casparian strips, transfers
PFAS to leaves more readily than pak choi; Sci. Total Environ.
sciencedirect.com/.../pii/S0048969724053555), and long-chain PFAS (C7–C14) are preferentially
**retained below ground** while short chains (C4–C6) move up (Environ. Int.
sciencedirect.com/.../pii/S0160412021002671). But the evidence supports it as a **translocation
barrier** (consistent with the model's existing `f_xy` and the equilibrium hydrophobic/
size framing), **not** as a quantified irreversible root **sink** with head-group or
chain-length-resolved coefficients. The panel **refuted** the stronger framing that the
Casparian strip is *the* mechanism for chain-length-dependent root retention distinct
from lipid partitioning (vote 1-2), and refuted a clean U-shaped C4–C8 PFCA root trend
in one source (0-3). So this is consistent with the model's split — **`f_xy` =
translocation barrier, `k_seq` = irreversible sink** — but does not anchor `k_seq`.

---

## 5. The central data gap = the decisive experiment

**No measured rice-root cell-wall PFAS sorption coefficient, and no rice Fe-Mn-plaque
PFAS partition coefficient, resolved by BOTH chain length AND head group, was found.**
All chain-length/head-group-resolved coefficients are for **analog** matrices (lignin
Kd, ferrihydrite sorption %, bulk-soil Koc). The model's `k_seq` is therefore anchored
by **analogy + superposition**, not by a direct root-sorbent dataset.

The decisive experiment (extends the data-gap list in `docs/HANDOFF_BAF_twopool.md`
§4.3 and `docs/twopool_root_exploration.md`):

1. **Isolated rice-root cell-wall fraction** (pectin / hemicellulose / lignin)
   **batch sorption** of a PFCA + PFSA chain-length series (C4–C12) → cell-wall Kd(n,
   head group), with a **desorption (hysteresis / TII) step** to measure
   irreversibility per congener. This directly tests the §1 short-arm-vs-long-arm
   superposition and the §2 head-group offset in the actual root matrix.
2. **Rice Fe-Mn plaque** (DCB-extractable plaque, or plaque-coated roots) PFAS
   sorption + desorption across the same series at **paddy-relevant pH** (circumneutral,
   not just acidic) → tests the §3 plaque verdict and any PFCA-specific inner-sphere term.

These convert `k_seq` from a **descriptor with a mechanistic analogy** into a
**directly anchored** rate — and would settle the §3 (this PR §4.2 / §Result 7)
promotion question.

---

## 6. Implications for the model

- **`k_seq` is now mechanistically defensible as a SUPERPOSITION**: a short-chain
  electrostatic/anion-exchange arm + a long-chain hydrophobic/desorption-resistant arm,
  with a sulfonate head-group offset — each anchored (HIGH confidence) in analog
  matrices. The two-exponential-arm functional form is no longer purely ad-hoc.
- **The irreversibility of the sink** (terminal accumulator, not a reversible K) is
  justified by measured desorption hysteresis (TII→0.98 at C10).
- **The PFSA offset** `10^(+0.18)` ≈ the measured **+0.23 log** sulfonate offset —
  a quantitative, non-circular corroboration.
- **Fe-Mn plaque** is demoted as the irreversible-sink candidate (cell-wall entrapment
  is the better attribution); a minor PFCA-specific plaque term stays open.
- **Promotion decision (§3 / PR #39 §4.2) unchanged**: this strengthens the *mechanistic
  story* but does NOT add a direct measured root-sorbent dataset, so it does not by
  itself warrant promoting the fitted `k_seq` into `parameters.json`. The §5 cell-wall /
  plaque batch-sorption experiment is the gate.
- Everything here is **literature/theory**; `parameters.json`, the model math, and
  `reproduce_demo` (RMSE 0.029) are **UNCHANGED**.

---

## Sources (verified, 3-vote adversarial panels)

| topic | source | DOI / URL |
|---|---|---|
| sulfonate +0.23 log Kd on lignin | J. Hazard. Mater. 2024 | sciencedirect.com/.../pii/S0304389424025950 |
| kaolinite MD, sulfonate vs carboxylate | ES&T 2025 | `10.1021/acs.est.5c01046` |
| Higgins & Luthy soil Koc, +0.23 log sulfonate | Soil Research | `10.1071/SR22183` |
| ferrihydrite, S K-edge XANES outer-sphere | ES&T 2020 | `10.1021/acs.est.0c01646` |
| ferrihydrite sorption % vs chain length / pH | PMC | https://pmc.ncbi.nlm.nih.gov/articles/PMC7745537/ |
| PFOS outer-sphere surface complexation | Geochem. Trans. 2025 | `10.1186/s12932-025-00105-2` |
| desorption hysteresis, TII=0.98 (PFDA) | Chemosphere 2015 | sciencedirect.com/.../pii/S0045653515302502 |
| biphasic soil Koc (short-chain electrostatic) | Environments 2023 | `10.3390/environments10100175` |
| rice plaque arsenate/phosphate molecular sieve | ES&T 2025 | `10.1021/acs.est.5c16125` |
| Casparian strip / radish translocation barrier | Sci. Total Environ. 2024 | sciencedirect.com/.../pii/S0048969724053555 |
| below-ground long-chain retention | Environ. Int. 2021 | sciencedirect.com/.../pii/S0160412021002671 |

*Generated by the deep-research harness (fan-out search → fetch → 3-vote adversarial
verify → synthesize). Confidence labels and vote tallies are from the verification
panels; claims that failed verification (e.g., a clean U-shaped C4–C8 root trend, a
specific two-sorbent split, electrostatics-alone for ferrihydrite) were dropped.*
