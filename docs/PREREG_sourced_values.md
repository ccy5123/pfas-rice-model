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

### 5.6 Addendum (committed before the re-fit; no fit had been run)
Corrections found while implementing 5.1:
- **TMPSB 0-5.0 C**, not 0-5.5. ORYZA's TMPSB raises the daily *maximum*. The ambient maximum over
  the seedbed days averages 20.0 C (JMA Tsukuba 2020, emergence to transplanting) and the Ibaraki
  calendar keeps the seedbed at 15-25 C (greening) and 10-25 C (hardening). So the cover adds
  between 0 and 25 - 20.0 = 5.0 C. The printed 0-5.5 compared a band midpoint with a mean
  temperature. Midpoint 2.5.
- **SBDUR midpoint 25.5** runs as the drivers of the 25-d and 26-d runs, linearly interpolated.
  ORYZA counts whole days, and the two ends move straw and grain by up to 0.006 log10, so either
  end alone is not admissible.
- **Observations built from the SI** (`docs/literature_db/tools/yamazaki_obs_from_si.py`):
  - water = Table S16 (it was back-calculated);
  - straw = the band concentrations weighted by each band's share of straw dry mass. The shares
    come from the stratified cutting of Tanaka et al. 2022 (Fig. 3, 10 and 49 d after heading).
    These two profiles bracket Yamazaki's harvest and are the range ends; the observation is
    their mean. The straw was an unweighted band mean, which has no source;
  - brown-rice milling share and < MQL cells at their midpoints, with both ends kept as `lo`/`hi`.

Fit protocol (`validation/refit_sourced.py`):
- **Variants**:
  - `redistributed_lipid` (the default);
  - `redistributed_free` (monotone `f_xy_recommended`);
  - `redistributed_fit` (f_xy fitted to straw; a reproduction);
  - `4pool_lipid`;
  - `4pool_free`.

  Each variant's fitted values are stored under its own key (parameters.json `r_eq__<v>`,
  `L_Ph__<v>`, `f_xy__<v>`, `lipid_loading[shoot]`). No fitted number is written in the code.
- **r_eq** in closed form: the model is linear in the state with zero initial state, so every
  tissue scales with r_eq. The script checks this numerically.
- **L_Ph >= 0** solves the grain/root ratio (log10 tolerance 1e-4) with no upper bound. If the
  model's ratio at L_Ph = 0 already reaches the observed one, then L_Ph = 0 (the physical bound;
  the 1e-5 floor is gone). A ratio that no L_Ph reaches is recorded as such.
- **`*_fit`**: f_xy in (0, 1] solves the straw/root ratio, alternating with L_Ph, to the same
  tolerance.
- **PFBS** has no grain value. Its L_Ph is the range [0, L_MQL], where L_MQL makes the model
  grain equal the censoring bound (the obs grain `hi`). It runs at the midpoint, and the range
  goes into the input-range band.
- **Per-congener interval**: each cultivar fitted on its own observation with its own harvest
  run. This is the bootstrap over the two cultivars.
- **Lipid constants per shoot**: profiled least squares. The per-congener r_eq and L_Ph are
  re-solved at every evaluation, on the residuals they leave: every congener's straw, plus the
  grain of any congener at L_Ph = 0.
  - Start from the current constants.
  - Physical bounds only: 0 < a <= 1, g >= 0, K_half > 0; b and pfsa_ln are free.
  - PFDoDA is included.
- **Identification by backward elimination** (rule 1.3). Each term in turn goes to its bound and
  the rest is refitted. The candidates are:
  - g_ph = 0;
  - g_xy = 0;
  - gate "none" (K_half -> 0) and gate "linear" (K_half -> infinity, g per L/kg);
  - pfsa_ln = 0;
  - b = 0.

  A term is removed when `decide(full, reduced)` on the profiled per-cell errors is not "worse"
  for the reduced model. Among several removable terms, the one whose removal raises the RMSE
  least goes first. Repeat until every remaining removal makes the fit "worse". Scoring on data
  not used to fit (rule 1.2) is the transfer endpoints, as in 5.3.
- **Read-across** for a record without a fit (GenX, a SMILES structure): the median of the fitted
  congeners' values. A fitted L_Ph can be exactly 0, so the old geometric mean no longer applies.

Expectations added before the fit:
- E14 g_ph is removed in both lipid variants: the per-congener L_Ph carries the grain.
- E15 At least one of PFBA/PFPeA sits at L_Ph = 0 in every variant.
- E16 The per-cultivar r_eq interval exceeds a factor of 2 for PFBA, PFPeA, PFHpA and PFBS, whose
  root values differ 5.8x, 3.6x, 3.1x and 2.4x between the cultivars (Tables S18/S19), and stays
  within a factor of 2 for the others.

### 5.7 Addendum (committed before the transfer runs; the re-fit was running, no transfer score computed)
- **Kim 2019 grain on the dry basis.** The Kim endpoints compared the model's FRESH grain BAF (at
  harvest moisture) with Kim's ng/g. Kim freeze-dried the grain at -30 C before dehusking and
  weighing 2 g for extraction (main p.715-716; `Sampling_Protocols.csv`), so the reported mass is
  dry mass. Every Kim comparison (`transfer_sourced.py`, `kim2019_censored.py`,
  `kim2019_pfoa_anchor.py`) now uses `baf_final_dw`. At the sourced grain water (0.21) this
  raises the model by log10(1/0.79) = 0.10 against the old comparison. The inverse
  (`estimate_exposure_bayesian`) already read Kim's value as dry.
- No expectation is attached. The shift is the same 0.10 log10 for every congener and variant, and
  the sign of the re-fitted variants' Kim errors is not known yet.

### 5.8 Addendum (committed before the neutral-path re-runs)
- **One water content per organ on the neutral path too.** Amendment 1 made the organ water
  `theta_fw` the sourced one (`params/rice_composition_sourced.csv`: root 0.713, stem 0.7095,
  leaf 0.695, grain 0.21). The dry -> fresh mass conversion (`organ_mass`) reads it, but the
  neutral path's own copy `neutral_dpu.RICE_WATER` (0.90/0.83/0.78/0.14, source unresolved:
  `value_provenance.csv` neu_water) was left. `tests/test_organ_mass.py` requires the copies to
  agree. `RICE_WATER` now reads `organ_mass.THETA_FW`. Found by a diagnostic run of the Liu 2023
  table with only the mass conversion changed: 0.213 -> 0.226.
- The neutral validation tables are re-run on this and their scores replace the old ones,
  whatever they are. The neutral lipid fractions are fresh-weight values reported next to root
  waters of 0.90-0.94. Restating them at the sourced, drier root needs a dry-basis lipid source
  and goes to amendment 2. No expectation is attached to the scores.

### 5.9 Addendum (committed before the transfer runs; no transfer score of amendment 1 computed)
- **Li 2025 stage -> model day, from sources instead of chosen DVS.** The Li endpoints enter the
  honesty band, and `li2025_stage_resolved.STAGE_DVS` (tillering 0.4, heading 1.0, filling 1.25,
  milk 1.5, ripening 2.0) had no source (`value_provenance.csv` tr_li_dvs). Li et al. 2025 sampled
  five stages without dates, in this order: tillering, heading, filling, milk-ripening, ripening
  (main text p.2). The mapping now uses only what sources fix:
  - heading = DVS 1. ORYZA's DVS 1 is flowering, and 50 % heading precedes anthesis by about 2 d
    (Counce et al. 2000 R3 -> R4, as reproduced in the Arkansas Rice Production Handbook MP192
    ch. 2);
  - ripening = DVS 2, physiological maturity, the end of the season;
  - tillering lies between transplanting and panicle initiation (DVS 0.65). Tillering ends at
    jointing, which coincides with panicle initiation (CMA agrometeorological observation
    criteria, as quoted by Li R. et al. 2018; MP192 ch. 2). It runs at the midpoint of
    [DVS at transplanting on the run, 0.65];
  - filling and milk-ripening lie between heading and ripening, in that order; no source found
    fixes them closer for Li's site. They run at the centre of that ordered range,
    1 + 1/3 and 1 + 2/3: the centroid of the region 1 < filling < milk < 2. As a check, the one
    dated japonica record found (Nanjing; Li R. et al. 2018, Table 2) puts filling 14-20 d after
    heading start; the script prints that span's DVS on the run, but does not use it.
- The fixed +-0.15 DVS sensitivity is replaced by the range ends: tillering at both ends;
  (filling, milk) at the corners of the ordered region, (1, 1), (1, 2) and (2, 2). The score
  reported is the midpoint run; the ends are printed next to it.
- No expectation is attached.

## 6. Outcome

(appended after the runs)
