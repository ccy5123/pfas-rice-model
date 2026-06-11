# GAP B — Root→shoot xylem-loading factor `f_xy(n)` (DELIVERABLE)

**Status: CLOSED.** Recommended values in `params/f_xy_recommended.csv` and
`params/parameters.json` (`f_xy_recommended`). `f_xy ∈ (0,1]` is the TSCF analog: the fraction of
root-internal free aqueous PFAS loaded into the ascending xylem. Theory in `docs/theory_anchor.tex`.

## Recommended values (MONOTONE)

`f_xy(n)` is a **monotonic decline** with perfluorocarbon chain length — *not* a Briggs bell.

| congener | n_C | f_xy_rec | congener | n_C | f_xy_rec |
|---|--:|--:|---|--:|--:|
| PFBA | 4 | 0.790 | PFDA | 10 | 0.0051 |
| PFPeA | 5 | 0.470 | PFUnDA | 11 | 0.0039 |
| PFHxA | 6 | 0.216 | PFDoDA | 12 | 0.0030 |
| PFHpA | 7 | 0.098 | PFBS (SA) | 4 | 0.176 |
| PFOA | 8 | 0.040 | PFHxS (SA) | 6 | 0.048 |
| PFNA | 9 | 0.014 | PFOS (SA) | 8 | 0.0089 |

**Closed form** (bounded logistic): `f_xy(n) = η / (1 + 10^(s·(n − n₀)))`, η = 1, s ≈ 0.372,
n₀ ≈ 4.74; equivalently `logit f_xy = 4.061 − 0.857·n` (OLS, R² = 0.97). Slope bracket
β ∈ [0.6, 1.1] C⁻¹ (LFER prediction), pooled empirical 0.67. Polar ceiling f_xy(PFBA) = 0.79 ≈
Briggs neutral maximum 0.784. Head group: **PFSA = PFCA · e^(−1.5)** (≈ +1.8–2.2 "effective
carbons"; sign uncertain, sensitivity term). The exponential long-chain tail and the short-chain
plateau are the two asymptotes of this single logistic.

## Theory (why monotone — Task A)

The DPU (Rein 2011 / Brunetti 2019) is the Trapp (2000/2004) ionizable-compound cell model at the
membrane. Factor `f_xy(n) = η(E_m, carrier) · φ_free(n)` with `φ_free(n) = θ_R / B_k(n)`. Since
every partition coefficient grows log-linearly with n (Collander/Briggs LFER) and one sorptive term
dominates at long chains, `φ_free` — and hence `f_xy` — is logistic in n. The **Briggs bell's rising
limb is removed** for a permanent anion by three independent facts: (i) f_d ≈ 1 → no ion trap;
(ii) inside-negative E_m → electrostatic anion exclusion (e^N ≈ 107) at all n; (iii) rising chain →
rising B_k (~3 orders C4→C12) → sequestration. Only the sorption-driven declining limb survives.

## Decisions

- **F1 — bell/U-shaped RCF belongs to root *sorption*, not f_xy. VINDICATED.** RCF and f_xy are
  different quantities (Trapp's two-process decomposition: RCF ~ sorption eq.5; TSCF ~ central-
  cylinder loading eqs.29–32). A bell/U RCF is a root-compartment property; f_xy stays monotone.
- **F4 — the IOC pH term enters mainly through E_m, not speciation (Task B).** With pKa ≤ 1,
  f_n ≈ 10⁻⁵–10⁻⁷ across pH 5.5–7.0 and is pH-invariant in absolute terms → ion-trap lever off. The
  live lever is **E_m**: NH₄⁺ depolarisation (−120 → −90 mV in flooded anaerobic paddy) relaxes
  e^N 107 → 33, raising passive anion influx ~2.5×. The pH-responsive congener **PFDA** is explained
  by E_m modulation (and possibly pH-dependent cell-wall Donnan charge in K_cw) in the steep low-f_xy
  regime — *not* by its own protonation. **Key data gap: in-situ paddy E_m** (−90 to −120 mV ⇒
  ~2.5–3× passive-influx spread; the dominant quantitative uncertainty in the pH term).

## ★ Reconciliation with the S6 transport fit (resolves H7 §7.2)

A separate S6 fit (`W2_transport_fit.csv`, `f_xy_W2fit`) re-fit f_xy to the Yamazaki root/straw/grain
BAF. It **reproduces** Yamazaki (full-ODE log10 RMSE 0.029) but is **non-monotone**: it rises for
C10+ (PFDoDA f_xy_W2fit = 0.67 vs f_xy_rec = 0.003, ~200×). This contradicts theory.

**Adjudication — the recommended monotone values are physical; the W2 long-chain rise is an
artifact:**
1. Theory **requires** monotone decline (above).
2. The fit is **saturated** (3 param / 3 obs per congener) — cost ≈ 0 is guaranteed, not validation;
   any B-structure mismatch is absorbed into the free parameters.
3. The absorbed mismatch is the **stem accumulation gradient**: Yamazaki upper-stem short-chain
   concentrations span 16→793 pg/g, which the single mass-weighted "straw" compartment cannot host;
   the fit inflates long-chain f_xy to compensate.
4. **Independent confirmation:** TF = tissue/tissue is **water-denominator-independent**, so it is
   immune to the Li2025 pore-water artifact that wrecks BAF. The cross-field TF_straw is **monotone-
   decreasing in BOTH** Yamazaki (16.2→0.42) and Li2025 (2.73→0.69) — exactly the theory's prediction
   (`validation/S6_Gap4*`).

**⇒ Use `f_xy_recommended`.** To use it *directly in prediction*, the stem compartment needs a
multi-height refinement (open modeling item); until then `f_xy_W2fit` is a structure-compensating
fit for reproducing Yamazaki only.

## Empirical anchoring (the values' provenance)

Absolute level = Felizeter solution-normalised TSCF (TSCF 0.05–0.25; PFBA ceiling ≈ 0.79); chain-
length **shape** = rice (JHM) / maize (Krippner) patterns; theory consistency = Trapp. Root uptake
is carrier-mediated/partly active (Michaelis–Menten, inhibitor-sensitive: Wen 2013/2014, Liu 2023
slow anion channel) — hence the GHK + carrier root influx (not pure passive). **Open:** absolute
scale is provisional pending measured Q_TP(t), M(t); a rice-direct solution-normalised TSCF and
xylem-sap / split-root data are absent.

## Sources

| key | reference | role | DOI |
|---|---|---|---|
| Trapp 2000 | *Pest Manag. Sci.* 56:767 | ionizable root-uptake/TSCF model (origin) | 10.1002/1526-4998(200009)56:9<767::AID-PS198>3.0.CO;2-Q |
| Trapp 2004 | *ESPR* 11:33 | neutral/ionic plant uptake | 10.1065/espr2003.08.169 |
| Briggs 1982 | *Pestic. Sci.* 13:495 | TSCF–logKow bell; RCF slope | 10.1002/ps.2780130506 |
| Rein 2011 | *SAR QSAR Environ. Res.* 22:191 | dynamic plant uptake (DPU concept) | 10.1080/1062936X.2010.548829 |
| Brunetti 2019 | *Water Resour. Res.* 55:8967 | DPU module for HYDRUS (framework) | 10.1029/2019WR025432 |
| Felizeter 2014/2012 | *JAFC* 62:3334 / *EST* 46:11735 | solution-normalised TSCF (absolute anchor) | 10.1021/jf500674j / 10.1021/es302398u |
| Krippner 2014/2015 | *Chemosphere* 94:85 / *JAFC* 63:3646 | maize chain-length / pH shape | 10.1016/j.chemosphere.2013.09.018 / 10.1021/acs.jafc.5b00012 |
| Liu 2023 | *ES&T* 57:8739 | slow anion-channel (carrier evidence) | 10.1021/acs.est.3c00504 |
| Yamazaki 2023 | *EST* 57 (es2c08767) | rice tissue distribution (validation) | 10.1021/acs.est.2c08767 |
| Li 2025 | *JHM* 492 | paddy field (cross-field TF; BAF water-caveat) | 10.1016/j.jhazmat.2025.138256 |
