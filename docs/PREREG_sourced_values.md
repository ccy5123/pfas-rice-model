# Pre-registration: no value without a source

Committed before any replacement run. The user's rule (2026-10-25): the tool may contain no
value chosen without a source. The inventory is `params/value_provenance.csv`; this file fixes
how each row is resolved and how the results will be judged, so the runs cannot steer either.

## 1. How a value becomes admissible (user's choice: staged)

1. **Literature.** A value read in a source (DOI or URL, page/table, verbatim quote), rice first,
   then other cereals with the species stated. A reported range is kept as a range.
2. **Fit**, only when the data identify it: the 95 % bootstrap interval of the fitted value
   (refit over resampled cells) is narrower than the physically possible range, and the fitted
   configuration is then scored on data not used to fit it.
3. **Simplify.** If neither, the mechanism (or the bound) that needs the value is removed, or
   reduced to the form the data do identify (e.g. a saturable carrier with unknown Km, Vmax
   becomes a linear entry when the dose response is linear over the data).
4. **Range.** If a value is still needed, the model runs at both ends of its physically
   possible range and reports the outputs as an interval.

Numerical values (floors, tolerances, grids, segment counts) are admissible only if the results
**converge**: moving the value toward its limit changes no reported number by 0.0005 or more in
log10 (half the last digit of a 3-decimal log10 score). A value whose results keep moving as it
approaches its limit is structural and must be fixed in the model (`validation/numerics_convergence.py`).

## 2. How a comparison is decided (user's choice: statistical)

Fixed margins (0.10, 0.3, 0.05/0.10, 1.5x, 90 %) are retired. A comparison of configuration B
against A on one endpoint uses `validation/decision_stats.decide`: paired bootstrap over the
scored items, 95 % percentile interval of RMSE(B) - RMSE(A); "better" if the interval is below 0,
"worse" if above, "tie" otherwise; the resample count doubles until the verdict and the interval
ends (3 decimals) are stable.

- A **default change** needs "better" on the endpoint it targets and "not worse" on every other
  scored endpoint (Yamazaki observable and the transfer endpoints).
- Datasets that chose the default (Tang, Kim, Li) keep their selection caveat; a dataset that
  was used to fit a value is never used to score it.

## 3. Order of work

1. Literature values in (composition, phloem, shoot shares, carrier, site management, soil).
2. Structural corrections that need no new number: each organ's own K_cw in the redistributed
   shoot; the inverse's measurement sigma = the default's transfer error
   (`APRIORI_LOG10_RMSE`) instead of 0.15; per-hill tables need the planting density as input.
3. Convergence fixes (the grain gate; see §4).
4. Re-fit the transport and lipid parameters on the new inputs, without the L_Ph upper bound
   (it binds for seven congeners and has no source).
5. Re-run every scored endpoint and the past decisions under §2; record the outcome below.

## 4. Expectations (written before the runs; misses will be recorded)

- E1 The gate's young-panicle peak does not converge as the organ-mass floor goes to 0
  (it grew 1.41x -> 9.65x over floors 1e-3 -> 1e-5 kg), while harvest values move less.
- E2 The BDF tolerances, the harvest sub-daily grid and the inner 1e-9 floor converge.
- E3 No rice Michaelis-Menten Km/Vmax for PFAS exists in the literature; the carrier becomes a
  linear entry.
- E4 Without the L_Ph upper bound, the long-chain L_Ph fits rise above 1.
- E5 Re-sourced organ composition moves the Yamazaki scores by less than the transfer-band
  width (log10 0.49), because the lipid and transport fits re-absorb it.
- E6 At least one shoot lever (retention, stem_transp_frac, lam_grain) has neither a source nor
  an identifiable fit and ends as a range output.

## 5. Amendment 1 -- the literature results and the structural changes (committed before the re-fit)

Four literature searches (root entry, organ composition, shoot water routing, site/soil) and the
Yamazaki 2023 SI (obtained from the ACS figshare mirror) returned values; each adopted value is in
`docs/literature_db/sourced_inputs.csv` with its quote, and the quotes marked `main` were re-read
against the source text. What they change:

### 5.1 Adopted values (literature; a range runs at its midpoint, rule 5.4)
| input | was | now | source |
|---|---|---|---|
| heading | harvest - 30.5 d | Kasalath 2020-08-03, Koshihikari 2020-08-11 (mean run: their mean) | Yamazaki SI p.S3 |
| nursery after emergence (SBDUR) | 23 d | 25-26 d (sown 04-13, transplanted 05-11, germination 2-3 d) | Yamazaki SI p.S2-S3; Ibaraki calendar (germination 28-30 C, 2-3 d) |
| seedbed warming (TMPSB) | 0 | 0-5.5 C (managed 15-25 C against an ambient seedbed mean of 14.6 C) | Ibaraki calendar; JMA Tsukuba 2020 |
| planting density | 19.05 hills m-2, 3 plants | 18.5 hills m-2, 4-5 plants (prefecture standard; the field's own is unreported) | Ibaraki calendar |
| irrigation-water exposure | back-calculated | Table S16 values | Yamazaki SI p.S31 |
| organ composition | untraceable table | root theta 0.713; f_prot 0.04-0.08; f_PL 0.019-0.022; f_cw 0.11-0.94 (no source: physical bound) / stem theta 0.645-0.774 (straw/shoot, no stem-only source); f_prot 0.034-0.11; f_PL 0.0019-0.0050 (straw ether extract x leaf PL share); f_cw 0.62-0.81 (straw NDF) / leaf theta 0.659-0.731; f_prot 0.088-0.319 (IR72DS92 leaf N x 6.25, maturity-tillering); f_PL 0.0085-0.0097; f_cw 0.62-0.81 (straw NDF; no blade source) / brown rice theta 0.20-0.22 (on the plant at harvest); f_prot 0.083-0.097; f_PL 0.0007-0.0017 (diacyl-total); f_cw 0.034-0.045 (NDF) | sourced_inputs.csv (composition rows) |
| phi | 0.1 | 0.02-0.10 (daytime phloem/xylem volume flow, non-cereal MRI; no cereal value exists) | Windt et al. 2006 |
| T_C,Ph | 10 L/kg | 4.4-7.8 L/kg (rice sieve-tube sucrose 177-225 g/L and 573.8 mM, with and without ORYZA's CRGSO) | Chino et al. 1988; Hayashi & Chino 1990 |
| panicle share of transpiration lam_grain | 0.05 constant | 0 before heading; 0.09-0.33 at heading (DVS 1); 0-0.14 at maturity (DVS 2); linear in DVS between (the form is an assumption, stated) | Sugimoto 1973 Table 2, 1975; Roberts & Batchelor 1983; Takahashi & Kurata 2007 |
| sheath + culm share stem_transp_frac | 0.45 | 0.22-0.36 (post-heading data, extended to the whole season: no pre-heading value exists) | Sugimoto 1973 Table 2, 1975 |

### 5.2 Structural changes (no new number)
- **Root entry**: the GHK + Michaelis-Menten pair (Km, Vmax with no source; uptake is linear at
  environmental levels in every series found) is replaced by the linear equilibrium entry
  `lam_eq * (r_eq * Cwo - Cw_root)`: `r_eq` fitted per congener (one root cell identifies it),
  `lam_eq` numerical (converged; an exploration put the change between 1e3 and 1e4 below 2e-4 log).
  E_m, kappa_d and the four carrier constants leave the default model.
- **Retention = 1**: solute delivered to an evaporating organ stays there (the module's own
  physics: a dead end); export from a leaf is the phloem term, modelled separately. The residual
  stream to the panicle becomes its own transpiration deposit, so `grain_share` (gate / size)
  no longer exists as a choice.
- **The 4-pool's leaf/grain xylem split** uses the same sourced lam_grain(DVS) share instead of the
  unsourced surfaces (20 and 0.5116 m2/kg).
- **L_Ph**: no upper bound (it bound seven congeners and had no source); the lower floor 1e-5 must be
  shown equal to 0 at reporting precision.
- **Removed from the API**: the carrier and bypass uptake modes (their constants have no source);
  the analyses that used them move to `validation/archive/`.
- Already done before this amendment: one converged mass floor and a step formation gate; each
  organ's own K_cw; the inverse's sigma = model transfer error (+ the user's measurement error).

### 5.3 Re-fit and decisions
- Re-fit on the redistributed shoot (the default) with the new inputs at their midpoints:
  per congener (r_eq, L_Ph); the six lipid-loading globals across congeners (PFDoDA included unless
  its fit runs to a bound, which is then recorded). The 4-pool and the free-anion variants are
  re-fitted the same way so the comparisons in 5.3 are like with like.
- Default configuration: the current default (redistributed + lipid) stays unless another fitted
  variant is "better" (section 2) on at least one transfer endpoint and "not worse" on all;
  otherwise the result is recorded as it is.
- The honesty band becomes the median transfer error of whatever default results.

### 5.4 Range inputs
A value known only as a sourced range runs at its midpoint. The output carries the interval from
the ends whenever an end moves a reported number by 0.0005 log10 or more: one input at a time, fits
held fixed, the log deviations summed per direction (a conservative envelope, checked against the
corner runs on the default scenario). `model_api.INPUT_RANGE_BAND` holds it for the default
scenario; the Simple view shows it next to the transfer band.

### 5.5 Expectations (written before the re-fit; misses will be recorded)
- E7 Every r_eq fit reproduces its root cell (saturated), and the lam_eq scan is converged by 1e4.
- E8 With retention 1 and lam_grain from the literature, the grain over-prediction of the default
  (2-7x) shrinks for most congeners before any re-fit.
- E9 After the re-fit, at least one long-chain L_Ph exceeds 1.
- E10 The Yamazaki observable is no worse than the current 0.380 (the fit has the same number of
  per-congener parameters).
- E11 The composition ranges dominate the input-range band for short chains (low K_PL), the shoot
  shares for the grain.
- E12 The transfer scores move by less than the current band (0.49 log) but the ranking of the
  variants does not change.
- E13 The input-range band on the grain is wider than log10 0.1 for most congeners.

## 6. Outcome

(appended after the runs)
