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

### 5.10 Addendum (committed before the neutral-path re-runs on it; amendment 1's outcome is recorded in section 6)
Addendum 5.8 put the sourced rice organ water into the neutral path but left the lipid fractions on the
basis they were reported on (fresh weight at ~90 % water), and converted every table's dry-weight
observations with the rice organ water whatever the plant (section 6, 5.8 diagnosis). Choices made by
the user on 2026-10-26: where no species or organ value exists, use Li et al. 2019's stated convention
(0.90 water for fresh roots; their 13-crop root lipid list as a range); Liu 2023's unstated basis is
scored both ways.

- **Lipid is carried on a dry basis** and restated at the organ's water: L_fw = L_dw·(1 − W). A source's
  fresh lipid is turned to dry with that source's own water (Li 2019: 0.90, so L_dw = 10·L_fw).
- **Rice compartments (the tool's neutral path)**: W = the sourced rice θ (unchanged); L_dw root =
  the cereal range of Li 2019 (maize 0.53, barley 1.00, wheat 1.10-1.14 % fw at 0.90 -> 5.3-11.4 % dw);
  stem and leaf = rice straw ether extract 0.9-2.1 % DM (Feedipedia, n = 161); grain = brown rice crude
  fat 1.9-3.3 % DW (FAO 1993 Table 14). Trapp 1994's soybean 1 % / 3 % leave the default; the
  `briggs_anchor` mode stays as the named alternative.
- **Each validation table runs on its own plant's composition**, and its dry-weight observations are
  converted with that same water (the ODE's dry -> fresh organ mass too, so K_PW and M share one water):

  | table | basis | water | dry lipid |
  |---|---|---|---|
  | Liu 2023 rice root, hydroponic 25 d | fw (inferred) and dw, both reported | root 0.824-0.888 (hydroponic rice seedlings: Choe 2026, Ya 2026) | cereal range |
  | Ge 2017 rice TF, soil pots 60 d | dw (stated) | root 0.845-0.851 (Saleh 2025, soil pots 56 d); leaf 0.731 (Hoang 2019, tillering) | rice compartments |
  | Kodesova 2019, 4 vegetables | dw (stated) | radish 0.928-0.945 (Ievinsh 2023 Table 3); lamb's lettuce, spinach, arugula 0.90 (Li 2019 convention) | spinach 3.4 %, radish 1.0 % (Li 2019); lamb's lettuce, arugula: Li's 13-crop range 1.0-11.4 % |
  | Li 2019 hydroponic | fw (Li, at 0.90) | 0.90 (Li) | barley 10.0 %, pumpkin (Cucurbita pepo) 7.0 %; rice cereal range; other species Li's 13-crop range |
  | Li 2019 soil | fw (Li, at 0.90) | 0.90 (Li) | each row's own f_lip (Li Table S2) / 0.10 |
  | Briggs 1983 barley shoot and stem, 11 d | fw | barley stem + leaf 0.861-0.875 (Wang 2017 Table 1; 34-d hydroponic plants) | barley straw ether extract 0.6-3.8 % DM (Feedipedia, n = 16) |

  Organs a table gives no value for take the rice compartments' values. Ge's stem rows are not scored:
  no rice stem water at the vegetative stage was found (the ODE's stem runs on the rice θ).
  Li's 0.90 is stated for roots, so it is not used for Briggs' shoots, which have a barley shoot value.
- **Ranges**: every score is reported at the midpoints, with its minimum and maximum over the corners
  (each range at its low and high end) next to it.
- The pre-5.8 code paths (rice water 0.90/0.83/0.78/0.14, unresolved) do not come back.
- **Expectations** (misses will be recorded): N1 the Li 2019 hydroponic and soil scores move to within
  0.05 of their pre-5.8 values (0.546, 0.550), since both return to Li's own basis; N2 Liu (fw reading)
  and Kodesova score better than their 5.8 values (0.256, 0.469); N3 the default rice root K_PW rises
  for lipophilic compounds (L_fw 1 % -> 2.4 %), close to the Briggs anchor (2.47 %); N4 no expectation for
  Ge (3 leaf rows) or Briggs.

## 6. Outcome

Runs of 2026-10-26 on amendment 1 (records: `validation/refit_sourced.json`,
`validation/transfer_sourced.json`, `params/input_range_band.json`, `params/neutral_scores.json`;
the checks of E7-E16 are read off them by `validation/prereg_a1_outcome.py`).

Numerics changed during the runs, both under the rule of section 1 (no result moved at reporting
precision except where stated):
- The plant ODE is affine in its state, so it is now integrated by exponential steps
  (`src/linear_expm.py`; 16 steps per driver interval, A frozen at the step midpoint) instead of BDF.
  Against BDF at rtol 1e-8 every reported number is unchanged (`numerics_convergence.py`); each
  model's vectorised (A, b) is held equal to its right-hand side (`tests/test_linear_expm.py`).
  Intervals where a rate changes fast get more steps (`linear_expm.GRADE` 0.2; 0.05 moves
  nothing). A run is ~10x faster; the full re-fit takes ~30 min instead of ~8 h.
- In the `*_fit` variant the alternation of the straw (f_xy) and grain (L_Ph) solves did not settle
  for PFDA (grain residual -0.367 recorded as solved); a joint solve now follows a non-settling
  alternation (PFDA f_xy 0.602, L_Ph 2.81, both cells exact).
- `LAM_EQ` 1e4 was not converged: 1e4 -> 1e5 moved the BAFs by 0.0008 log10. Raised to 1e6
  (1e5 -> 1e6: 0.00008; 1e6 -> 1e7: 0.00001) and the re-fit repeated.
- `hydrus_drivers` dropped the drivers' DVS, so `simulate(drivers=hydrus_drivers(...))` failed on
  the default shoot (the panicle share follows DVS): fixed, with a test.

Section 4 (expectations before amendment 1):
- E1 HELD. The old ramp's young-panicle peak grew as the floor fell (1.41x -> 9.65x); the gate is now
  a step at the floor and converges (floors 1e-8 to 1e-12: no number moves). Under the sourced
  stage share (0 before heading) the PFOA panicle has no pre-heading peak at all (max at harvest).
- E2 HELD. BDF tolerances, the harvest sub-daily grid and the inner floor converge (and the
  exponential integrator's step count).
- E3 HELD. No rice Km/Vmax for PFAS exists; the carrier became the linear entry (5.2).
- E4 HELD. Unbounded, PFDA L_Ph = 6.44 and PFDoDA 13.2 (default); PFUnDA sits at 0.
- E5 HELD, confounded. The Yamazaki observable moved 0.380 -> 0.159, less than 0.49, but with every
  amendment-1 change at once, not composition alone.
- E6 HELD. retention became derived (1); the stem and panicle shares are sourced ranges carried in the
  band; no lever is a free fit.

Section 5.5 / 5.6:
- E7 PARTLY WRONG. Every r_eq reproduces its root cell (residual 0); the lam_eq scan was NOT converged
  at 1e4 (above).
- E8 HELD. Before the re-fit, retention 1 and the sourced panicle share cut the grain over-prediction
  for 7 of 9 congeners (PFBA 6.6x -> 2.24x, PFPeA 7.1x -> 2.61x, PFHpA 6.9x -> 2.36x, PFOA 3.2x -> 2.60x,
  PFNA 2.9x -> 2.67x, PFDA 1.9x -> 1.87x, PFOS 2.3x -> 2.00x); PFUnDA (7.3x -> 8.05x) and PFDoDA (1.1x ->
  1.76x) rose.
- E9 HELD (PFDA 6.44, PFDoDA 13.2). E10 HELD (0.159; root 0.000, straw 0.253, grain 0.103).
- E11 HELD. Composition leads the root and straw band for PFBA, PFPeA, PFHxA, PFBS; a shoot share
  (mostly the panicle share at heading) leads the grain band for 11/11.
- E12 WRONG. Tang 0.729 (+0.089), Kim point 0.489 (+0.020), Li shoot 0.386 (+0.028) moved less than 0.49,
  but Li grain moved +0.509 (0.510 -> 1.019): the grain is now under-predicted on Li (mean -0.98 log).
  The default stays (5.3): no variant is better on a transfer endpoint and not worse on all.
- E13 HELD (grain band > 0.1 log10 for 11/11; PFOA grain -0.38/+0.36). E14 HELD (g_ph removed in both
  lipid variants; b, pfsa_ln and the K_PL gate were removed too, so both lipid shoots ship one f_xy
  = a and one g_xy for every chain). E15 HELD (PFBA and PFPeA at L_Ph = 0 in every variant). E16 HELD
  (spans PFBA 5.76x, PFPeA 3.65x, PFHpA 3.12x, PFBS 2.40x; all others < 2x).

Scores (log10 RMSE; default = redistributed_lipid): Yamazaki 0.159 (in-sample); Tang 0.729; Kim point
0.489, interval 0.371; Li shoot 0.386, grain 1.019; honesty band = their median 0.489 (x3.1). Variants:
redistributed_free 1.266 / 1.757 / 1.967 / 1.710 / 1.343 / 1.405; redistributed_fit 0.152 / 0.695 / 0.554
/ 0.444 / 0.510 / 1.183; 4pool_lipid 0.157 / 1.030 / 0.504 / 0.385 / 1.014 / 1.027; 4pool_free 1.359 /
3.693 / 1.934 / 1.669 / 3.299 / 1.378 (Yamazaki / Tang / Kim point / interval / Li shoot / grain).
In the free variants, L_Ph of an unreachable long chain is the ratio's maximum; where the ratio rises
monotonically to an asymptote below the observation the stored value is merely large (up to ~1e9) and
any larger value gives the same grain.

Addenda:
- 5.7 (Kim on the dry basis): applied; no expectation.
- 5.8 (one water per organ on the neutral path): applied. Scores Liu 0.256 (equilibrium 0.229), Ge
  0.702, Kodesova 0.469 (0.392), Li hydroponic 0.568, Li soil 0.549 (0.639), Briggs stem 0.313 -- all
  worse than before (0.213, 0.311, 0.216, 0.546, 0.550, 0.299). Diagnosis: a basis mismatch the addendum
  missed -- the neutral lipid fractions are fresh-weight fractions reported at ~90 % water, and the
  comparison converts other species' dry-weight observations with the rice root's water. Restating
  the neutral composition on one basis is amendment 2.
- 5.9 (Li stages from sources): applied; the range ends are printed by
  `validation/li2025_stage_resolved.py`.
- Reopened (`params/value_provenance.csv` tis_theta): organ water falls with stage (root 0.76-0.89
  vegetative, 0.51-0.77 at maturity); a constant theta is amendment 2.

### Addendum 5.10 outcome (2026-10-26; record: `params/neutral_scores.json`)
Scores at the midpoints, log10 RMSE, [min, max] over the range corners; equilibrium basis in brackets:
Liu 2023 0.212 [0.208, 0.278] (0.204) · Ge 2017 0.631 [0.568, 0.724] (leaf rows only) · Kodesova 2019
0.157 [0.160, 0.193] (0.172) · Li 2019 hydroponic 0.677 [0.525, 1.115] · Li 2019 soil 0.418 (each row's own
lipid; no range) (0.461) · Briggs 1983 barley stem 0.357 [0.317, 0.417]. Liu read as dry weight: 0.821.
- N1 WRONG. Li hydroponic 0.677 (+0.131 against 0.546) and Li soil 0.418 (-0.132 against 0.550) both
  moved by more than 0.05, in opposite directions. Not foreseen: the hydroponic species without a value in
  Li's list run at the 13-crop range midpoint, and the hydroponic half is under-predicted (mean bias
  -0.56); Briggs' fitted root term scores it 0.297.
- N2 HELD. Liu 0.212 (5.8: 0.256), Kodesova 0.157 (5.8: 0.469).
- N3 HELD. The rice root runs L_fw 2.40 % (Briggs anchor 2.47 %).
- N4 (no expectation): Ge 0.631 (leaf 2.5-5.8x over without metabolism), Briggs stem 0.357 (the
  model 1.8x under on average).
- Liu's fresh-weight reading (0.212) fits far better than the dry-weight one (0.821).

## 7. Amendment 2a -- organ water, the grain basis and the crop parameters (committed before its runs)

Amendment 1 left three default-scope rows open (`tests/test_value_provenance.py` OPEN: `tis_theta`,
`crop_cultivar`, `obs_exposure_time`), and an external review of the model (P1-P9, read against the
branch on 2026-10-26; `docs/amendment2/draft.md`) found further defects in the default model and the
measured-table path. Choices made by the user on 2026-10-26:
- amendment 2 is registered in two parts. **2a** (this section): the default model and its re-fit.
  **2b** (registered later, before its own runs): the neutral path, the site soil and the
  irrigation-source exposure, the app values and the SMILES `K_cw` (section 7.9);
- organ water follows the development stage, linear between sourced stage anchors;
- the rachis share is searched for before this commit, and run as a range if no source exists;
- the exposure stays constant over the season in 2a; its seasonal shape is built in 2b with the
  irrigation-water scenario.

Two further choices were proposed in the session and not answered by the user before this commit; they
stand as written unless changed by an addendum before any 2a run: Li 2025's grain is scored on both the
brown-rice and the rough-grain basis (7.2); Choi 2018 is a check, not a scored endpoint (7.2, E21). The
destination of the panicle's transpired solute was asked too and is settled by the sources read (7.2).

The rules of sections 1-2 and the fit protocol of 5.6 are unchanged. Nothing in this section has
been run on the amended model; the prototypes in `docs/amendment2/prototypes/` were built on the
pre-amendment-1 tree, were not fitted, and two of them (0001, 0003) no longer apply.

### 7.1 Organ water by development stage (`tis_theta`; review P2)
- theta_fw(DVS) per organ, linear in DVS between sourced anchors and flat outside them
  (`params/rice_theta_stage_sourced.csv`; the linear form is a stated assumption, as for the
  panicle's transpiration share in 5.1). Anchors at the range midpoints:

  | organ | vegetative (DVS <= 0.65) | anchor 2 | sources (`sourced_inputs.csv`) |
  |---|---|---|---|
  | root | 0.8245 (0.761-0.888; 4 studies, hydroponic and flooded pots) | DVS 2: 0.6365 (0.506-0.767; 6 studies, pots at maturity) | tis_theta_root |
  | stem | 0.818 (0.726-0.910; whole shoot, no stem-only value) | DVS 2: 0.7095 (0.645-0.774; whole shoot 10-45 DAF and straw) | tis_theta_stem, tis_theta_root (shoot columns) |
  | leaf | 0.731 (tillering, 183 accessions) | DVS 2: 0.659 (flag leaf at maturity) | tis_theta_leaf |
  | grain | see 7.2 (the anchors must describe the pool the grain becomes) | | tis_theta_grain |

- **One theta per organ and time for every use.** The fresh mass (`organ_mass`), the binding
  factor B_k(t) and the harvest operator read the same theta(t). With a moving theta,
  `dM_fresh * (1 - theta)` is no longer the dry growth, so the grain's phloem stream takes the dry
  growth from the dry driver; the harvest operator divides burden by dry matter at the harvest time.
- The burden P = C * M_fresh is the conserved quantity; the free concentration is
  `Cw = P / (M_dry * (theta/(1-theta) + X))` with X = sum f*K on the dry basis. Mass conservation
  of both shoots, the vectorised (A, b) against the right-hand side, the dry growth and the harvest
  dry basis are tested (prototype `tests/test_theta_stage.py`, carried over).
- **Band**: each anchor's ends, one at a time (rule 5.4). The amendment-1 constant theta rows leave
  `INPUT_RANGE_BAND`.

### 7.2 The grain pool is brown rice; the panicle's transpired solute stays in the husk (review P3)
The grain pool is WSO, ORYZA's whole panicle ("Storage organs (panicles)", Bouman et al. 2001 p.169;
"Measured panicle dry weight (storage organ)(WSO)" p.161), while its composition, its water and every
grain observation (Yamazaki, Kim, Tang) are brown rice. Read for this section (main session, full text,
2026-10-26; `sourced_inputs.csv` grain_*):
- **Where xylem water goes in the panicle.** Zee 1972 (Aust. J. Biol. Sci. 25:411, p.411): dye fed to the
  culm xylem reached "the vascular bundles of the rachilla, glumes, lemmas, and paleas" within minutes,
  "but no dye was detected in the tracheary elements of the pericarp bundles unless the caryopses were
  exposed by carefully separating the lemma and palea ... and the plant allowed to transpire in the dye
  solution for 3-5 hr" (the xylem is continuous, unlike wheat, but carries no flow to the intact
  caryopsis). Oparka & Gates 1981 (Planta 152:388, p.393-394): caryopses took up tracer only under high
  humidity, "under low humidity solutions moved instead into the palea and lemma"; "most of the water
  which enters the rice caryopsis during grain filling arrives by mass flow in the phloem", and water
  leaves the caryopsis during filling, possibly to "the xylem of the palea and lemma prior to
  transpiration".
- **What the husk holds.** Choi et al. 2018 (J. Appl. Biol. Chem. 61:191, Table 1; cv. Sindongjin,
  greenhouse pots of spiked soil, harvest 138 d after transplanting, husk and brown rice separated by a
  dehuller): husk / brown rice = PFOA 8.8 (1 mg/kg soil) and 4.2 (5 mg/kg), PFOS 37 and 47; PFOS exceeds
  PFOA in straw and husk but not in brown rice. Contrast for a phloem-mobile element: Cd in rice, phloem
  ranked "rachis > brown rice > rice husk" (Luo et al. 2022, Ecotoxicol. Environ. Saf. 240:113683).

Changes (no new number):
- **Two panicle pools.** Brown rice (the grain pool: composition `grain_brown`, the observations) and the
  panicle remainder (husk + rachis). Brown rice dry mass = (WSO − R)·(1 − h(t)), the remainder takes the
  rest of WSO.
- **The panicle's transpiration deposit goes to the remainder**, which keeps it (retention 1, as for
  every evaporating organ in 5.2). **Brown rice receives the phloem stream only**, drawn by its own dry
  growth (`T_C_Ph` per kg brown-rice DM). The xylem route into the caryopsis is removed (rule 1.3): the
  sources above give no flow into the intact caryopsis and no number for one. The remainder loads no
  phloem (in the model only the leaf does); it is a terminal pool.
- **Husk share of rough grain h(t)** (dry). h = 1 at DVS 1 (the caryopsis starts at fertilisation), then
  the per-grain hull and brown-rice weights of Ansari et al. 2017 (Bangladesh Agron. J. 19(2):125,
  10.3329/baj.v19i2.31861, Fig. 1, Nipponbare; heading, 10, 20, 30 d after heading and maturity), digitised
  and overlay-checked against the figure before the runs (`raw_si/ansari2017_fig1_digitised.csv`), mapped
  to DVS on the run's own heading date; linear in DVS between them (the stated form). Maturity range:
  Yoshida 1981 Table 1.1 (rough rice, 11 varieties 0.186-0.254; japonica 0.206, 0.209) with Ansari's
  Nipponbare value; the run uses the midpoint and the band carries the ends.
- **Rachis R: no source** (searched English, Japanese, Chinese, Korean; Matsui & Kagata 2002 weighs
  spikelets only). Range (rule 1.4): ORYZA puts all panicle growth after DVS 0.95 into grains (Bouman et
  al. 2001 p.57), so the rachis is complete by then and cannot exceed the panicle at DVS 0.95 (on the
  site-year mean run 11.6 % of WSO at harvest; that panicle also holds the hulls formed before
  flowering, so the upper end over-states R). R = [0, WSO(DVS 0.95)], run at the midpoint, ends in the
  band. A rachis weight from a source replaces the range (lead: Kobayashi & Horie 1994, Jpn. J. Crop Sci.
  63:193).
- **Grain water before harvest enters no dry-basis result**: the brown-rice equation has no term in its
  own free concentration (deposit, phloem from the leaf, growth dilution only), and the phloem stream reads
  the dry growth. So brown rice runs on its sourced harvest water (0.20-0.22) for the fresh display and
  the harvest operator; no filling-stage water is needed. A test holds that changing the grain's
  pre-harvest theta moves no dry-basis number. (The 0.58 "early filling" anchor of the draft is whole
  rough grain, the wrong basis for the brown-rice pool; it is not used.)
- **Observation bases**: Yamazaki and Kim = brown rice; Tang endosperm = brown rice ×
  `endosperm_white_factor` (unchanged). Li 2025 "paddy rice / grains" is undefined (it may include the
  husk): scored on both bases, brown rice and rough grain (brown + remainder, the rachis included because
  the pools are not separable), each reported, the range as in 5.10. Choi 2018 is a check, not a scored
  endpoint (section 2 is unchanged): it is not used to fit anything.
- **Band**: h at maturity, R, and the phloem draw of the remainder (0 = shipped; the other end lets the
  husk + rachis dry growth draw phloem sap too, as for brown rice).

### 7.3 The measured-table path (review P4, P5)
- **P4 leaf loss.** `drivers_from_tables` passes no leaf loss, so a falling measured green-leaf mass
  concentrates the leaf (mu < 0), the artefact the ORYZA path removed (CLAUDE.md, leaf senescence).
  A table that gives green leaf mass only gets `leaf_loss = max(0, -dM_leaf/dt) / M_leaf` (the same
  identity as ORYZA's D/M_leaf = DRLV: a falling green mass is read as death, not as export). A table
  with a dead-leaf column uses it. No new number.
- **P5 interpolation.** `np.interp` of measured masses makes mu = (dM/dt)/M a step function that
  jumps at every sampling date. The form between samples is not measured, so the path runs two
  forms and reports both: piecewise-linear in M (as now) and monotone cubic (PCHIP, continuous mu).
  The midpoint convention does not apply to a choice of form; the app shows the linear run and the
  spread. No scored endpoint uses this path, so no score moves.

### 7.4 Planting density in the generic scenario (review P7)
- 25 hills/m2 (0.04 m2/hill) has no source and is written five times: `forcing_rice`,
  `growth_rice` and `oryza_growth` (`AREA_PER_HILL_M2`), `oryza2000` (`HA_PER_HILL` and
  `Management.nh`), read again by `soil_hydrus`. Checked on 2026-10-26 (generic scenario, PFOA,
  120 d): with the transpiration and the ORYZA2000 mass both at 1/18.5 m2/hill every BAF is unchanged
  (log10 difference 0.000000), because Q_TP and M per hill scale together. With the transpiration
  copy changed alone, the stem BAF moves +0.131 and the grain +0.122 log10. So the value only sets
  per-hill absolute amounts, but the copies can drift apart and then it sets the BAFs.
- One constant in one module, at the sourced 18.5 hills/m2 (Ibaraki standard, 5.1), read by every
  driver; a test that the generic BAFs are invariant to it and that no second copy exists; every
  per-hill display states the density (or shows per m2).

### 7.5 Japonica crop sets as a crop-parameter range (`crop_cultivar`; review P6)
- No ORYZA set exists for Koshihikari or Kasalath (searched; `sourced_inputs.csv`
  crop_cultivar_search). Two published japonica sets run beside IR72, each with phenology
  recalibrated to Yamazaki's measured heading and harvest (as IR72 is now):
  - Gao et al. 2021 (Suijing-18, Heilongjiang; `raw_si/gao2021_cold_japonica_oryza_params.csv`);
  - Wu et al. 2022 (Liaoning; ORYZA v3 partitioning FLV/FST/FSO to DVS 1.2 and SLA; the rest IR72, as
    the authors state; `raw_si/wu2022_liaoning_oryza_params.csv`, to be verified against the PDF
    before the runs). The Kaiyuan partitioning duplicates Donggang (a possible copy error in the
    paper): Kaiyuan is not used.
- A categorical input has no midpoint: the default stays IR72 (the shipped run and the re-fit) and the
  two sets enter the band as ends (their log deviations from the IR72 run, per direction).

### 7.6 Documentation only (review P8, P9)
- P8: the model root is a root without plaque (`K_surf` off, not identifiable); every root comparison
  states whether the observed root kept its plaque (`Sampling_Protocols.csv`).
- P9: the stem's dry mass falls after heading (port: -13 % from its maximum to harvest; remobilisation),
  so mu < 0 concentrates the stem. This is the crop model's carbon balance, kept; its size on the
  stem BAF is quantified and recorded (theta(t) changes it).

### 7.7 Re-fit and decisions
- Re-fit every variant by protocol 5.6 on 7.1-7.2 (the corrections of 7.3-7.4 do not touch the
  default scenario). Fitted values replace amendment 1's under the same keys.
- 7.1-7.4 are corrections of a defect (a value or form without a source, or a wrong basis), adopted
  whatever the scores do, as in amendment 1. The default variant is then re-decided by 5.3, and the
  honesty band is the median transfer error of the resulting default.
- Every scored endpoint and the input-range band are re-run and recorded in section 8.

### 7.8 Expectations (written before the runs; misses will be recorded)
- E17 The root fit absorbs the harvest root water in closed form. With entry at equilibrium,
  `conc_dw(root) = r_eq * Cwo * (theta/(1-theta) + X)`, so the new r_eq / amendment-1 r_eq is
  `(2.484 + X)/(1.751 + X)` (theta 0.713 -> 0.6365): PFBA 1.13, PFPeA 1.10, PFHxA 1.06, PFHpA 1.03,
  PFBS 1.03, PFOA 1.02, others < 1.01 (X from the shipped binding). Each holds within +-0.02.
  (The draft's E-a, ">20 % for short chains", was a guess this calculation contradicts; it is
  replaced.)
- E18 After the re-fit the Yamazaki observable moves by less than 0.05 log10 (0.159).
- E19 The root's maturity water anchor is among the three largest contributors to the root band
  for PFBA, PFPeA, PFHxA and PFBS.
- E20 With the deposit moved to the husk, every congener with a grain cell fits L_Ph > 0 (amendment 1
  left PFBA, PFPeA, PFHpA and PFUnDA at 0 with the grain still over), and the straw and root fits do not
  move (the husk is not in Yamazaki's straw).
- E21 Choi 2018 check (the model on the site-year drivers, Choi's sites not modelled): the model's husk +
  rachis / brown-rice dry concentration ratio exceeds 1 for PFOA and PFOS, and is larger for PFOS than for
  PFOA (Choi: 4.2-8.8 vs 37-47). The log10 distance to Choi's range is recorded.
- E22 Li 2025 grain scores better on the rough-grain basis than on the brown-rice basis (the default
  under-predicts Li's grain by -0.98 log on average; the husk raises the rough-grain value).
- E23 The crop-set ends move straw and grain by >= 0.0005 log10 for every congener (they enter the
  band), and the panicle's transpiration share at heading stays the largest grain-band contributor
  (as E11).
- E24 Every BAF of the generic scenario is invariant to the planting density (a test).

### 7.9 Registered later (2b)
Neutral path (organ water theta(t) or constant, one value for K_PW and the organ mass; dry-basis
lipid; `neu_lipid_a` 1.22, `neu_leafexp` 0.95, `neu_surface` 20/2, `neu_kaw_cut` 1e-4,
`mem_ph` 6.5 -> paddy pore-water pH 6.5-7.0; volatile-compound data); site soil and the
irrigation-source scenario (`soil_*` rows, `obs_exposure_time`); app (`app_density`, `smiles_kcw`,
KNHANES intake and body weight, the white/brown basis of the intake fraction); analysis rows
(`tr_drivers`, `tr_li_grain_from`, `calib_sigma`; `fxy_recommended`).

### 7.10 Addendum (committed before any amendment-2a run; found while digitising Ansari Fig. 1)
- `raw_si/ansari2017_fig1_digitised.csv`, checked by `docs/literature_db/tools/ansari2017_fig1_check.py`
  (axis ticks re-fitted; every marker standing alone re-detected; merged markers checked on the overlay).
  Nipponbare hull / brown rice per grain [mg]: heading 3.02 / -; 10 DAH 3.80 / 7.5 (merged, 7.0-8.0);
  20 DAH 3.84 / 16.67; 30 DAH 4.05 / 20.50; maturity 3.80 / 19.8. Husk share of rough grain h:
  10 DAH 0.336, 20 DAH 0.187, 30 DAH 0.165, maturity 0.161. Panel B's y label is misprinted ("A hull
  weight"); the caption (brown rice) is followed.
- **The form of h(t), made explicit.** 7.2 put the maturity anchor at the midpoint of a range that also
  contains indica values (0.161-0.254, midpoint 0.2075), above Nipponbare's own 30-DAH value (0.165): read
  literally, h would rise again at the end of filling, which no series shows. Instead the Nipponbare course
  is scaled so that its maturity value equals the maturity anchor, h(t) = h_N(t) * h_mat / h_N(maturity),
  for t >= 10 DAH; h = 1 at DVS 1, linear in DVS to the 10-DAH value. h_mat keeps the range of 7.2
  (midpoint 0.2075; ends 0.161 and 0.254 in the band).
- DAH -> DVS: on each cultivar's own run, the DVS reached the given number of days after its measured
  heading date; maturity = DVS 2.

### 7.11 Addendum (committed before any amendment-2a run; found while implementing 7.1-7.2)
- **Stem vegetative anchor re-read at source**: 0.726-0.906, not 0.726-0.910 (midpoint 0.816, not
  0.818). The prototype's 0.910 is Ya et al. 2026 CYQH FD4 control, 1 - 0.16/1.71 = 0.906, rounded.
  The anchor's four sources are now in `sourced_inputs.csv` (`tis_theta_shoot_veg`: Saleh 2025
  SG 65 0.726 and Giza 177 0.737; Ya 2026 controls 0.813-0.906; Choe 2026 0.873; Haruma 2023 rice
  0.748 -- the same paper's Miscanthus rows, ~0.9, are not rice), read in full from Europe PMC.
- **The grain is not staged** (7.2): `params/rice_theta_stage_sourced.csv` holds root, stem and
  leaf only; brown rice keeps its harvest water (0.21). The neutral path keeps the constant water
  and the whole panicle until 2b (`simulate_neutral` passes `organ_water="constant"`,
  `grain_basis="panicle"`).
- **The husk + rachis pool is carried per kg DRY** with its dry mass as the driver (it is
  terminal: nothing leaves it, so it needs no binding factor and no water). Its mass from the
  stated forms (7.2, 7.10) dips by up to ~5 % between Ansari anchors; a shrinking terminal pool
  concentrates by the same amount. Kept (no sourced form to correct it with); it reaches only the
  husk and rough-grain values (Li rough basis, the Choi check).
- **Implementation choices, no new number**: `model_api.simulate(organ_water="stage"|"constant",
  grain_basis="brown"|"panicle", husk_phloem=False)`; defaults are the amendment (7.7); the
  amendment-1 forms stay runnable. The stage anchors' band ends pass as
  `scenario_kw={"theta_anchors": ...}`, the panicle's as `scenario_kw={"panicle_kw": {...}}`.

### 7.12 Addendum (committed before any amendment-2a run; found while implementing 7.3-7.5)
- **7.4, the density value: the generic crop run's own, not 18.5.** 7.4 put the one constant at the
  sourced 18.5 hills/m2. But the generic scenario's density is not only a per-hill conversion: the
  ORYZA2000 port runs the IR72DS92 experiment file verbatim, whose NH = 25 hills/m2 is a crop input
  (establishment density), and a per-hill conversion at another density than the run's would be
  inconsistent. So the one definition is `oryza2000.GENERIC_HILLS_PER_M2 = Management().nh`
  (inherited with the crop file), from which `AREA_PER_HILL_M2` / `HA_PER_HILL` derive and which
  `forcing_rice`, `growth_rice` and `oryza_growth` import; the site-year keeps its sourced 18.5
  (`site_scenario.NH`). What 7.4 asked for holds: one definition, no second copy (a test greps for
  one), and the generic BAFs are invariant to the value (a test moves it to 1/18.5; relative change
  < 1e-6, the remainder being the fixed 1e-9 kg organ-mass floor).
- **7.5, which Wu sets**: Huanren, Donggang, Liaozhong and Dawa (Kaiyuan excluded, as registered).
  With Gao 2021 that is five band "ends"; the band takes the most negative and the most positive
  deviation from the IR72 run over them (`crop_sets.BAND_SETS`).
- **7.3, the measured-table grain**: a growth table gives no husk mass, so its grain column is taken
  as given and receives the panicle's transpiration deposit (the amendment-1 grain); stated in
  `drivers_from_tables`. Its interpolation form: `interp="linear"` (shown) and `"pchip"`; the
  Expert view prints the monotone-cubic run's deviation.
- **The app's driver cache** now carries `theta` and `husk` (`ui/common._drivers_tuple`): without
  them a staged fresh mass ran with constant binding and the panicle deposit went back to brown rice.
- Before the re-fit, the tests that compare a fitted score against its record
  (`test_kim2019_censored`, `test_plots::test_fig_congener_compare_builds`) fail: with the deposit in
  the husk, the congeners amendment 1 left at L_Ph = 0 have no grain at all until re-fitted (E20).

### 7.13 Addendum (committed while the amendment-2a re-fit runs; no transfer score of 2a computed)
- **Which Li 2025 grain basis decides.** 7.2 scores Li's undefined "grains" on two bases. The
  decision rule (section 2) and the honesty band use the BROWN-rice basis -- the model's grain pool
  and the basis of every other grain endpoint (Yamazaki, Kim, Tang) -- under the same endpoint name
  `li_grain`. The rough-grain basis is reported beside it (`li_grain_rough` in
  `validation/transfer_sourced.json`) and is the one E22 is about; it enters no decision.
- **Choi 2018** runs as `validation/choi2018_husk_check.py` (E21): the model's husk + rachis pool
  against Choi's husk, on the site-year drivers; nothing is fitted to it.
- The amendment-1 records are kept in `validation/archive/amendment1/` (refit and transfer JSON, the
  band, the transfer scores, the fitted parameter files) before the re-fit overwrites them.

## 8. Outcome (amendment 2a)

Runs of 2026-10-26 on amendment 2a (records: `validation/refit_sourced.json`,
`validation/transfer_sourced.json`, `validation/input_range_band__redistributed_lipid.json` (the band
E19 and E23 were registered on), `params/input_range_band.json` (the band of the default chosen
below); the checks of E17-E23 are read off them by `validation/prereg_a2_outcome.py`; the amendment-1
records are in `validation/archive/amendment1/`).

Re-fit (log10 RMSE over 32 cells: all / straw / grain; root 0.000 in every variant):
redistributed_lipid 0.154 / 0.262 / 0.011 (amendment 1: 0.159 / 0.253 / 0.103) · redistributed_fit
0.116 / 0.156 / 0.127 · 4pool_lipid 0.151 / 0.258 / 0.000 · redistributed_free 1.261 / 1.639 / 1.461 ·
4pool_free 1.309 / 1.794 / 1.394. Backward elimination left `f_xy = a` 0.0659 and `g_xy` 0.00712 on
redistributed_lipid (0.0635 / 0.00730 on 4pool_lipid); b, g_ph, pfsa_ln and the K_PL gate were
removed again. In redistributed_fit PFDoDA alone is out of reach (f_xy at 1; straw -0.52, grain -0.40).

Transfer (log10 RMSE; Tang / Kim point / Kim interval / Li shoot / Li grain on brown rice; Li grain on
rough grain in brackets, no decision):

| variant | Tang | Kim point | Kim interval | Li shoot | Li grain [rough] |
|---|---|---|---|---|---|
| redistributed_lipid (incumbent) | 0.687 | 0.554 | 0.444 | 0.364 | 1.361 [0.864] |
| **redistributed_fit** | 0.686 | 0.562 | 0.443 | 0.467 | 1.334 [1.042] |
| 4pool_lipid | 0.916 | 0.559 | 0.445 | 0.937 | 1.373 [0.885] |
| redistributed_free | 1.738 | 2.118 | 1.866 | 1.279 | 1.863 [1.463] |
| 4pool_free | 1.876 | 2.037 | 1.777 | 1.494 | 1.815 [1.437] |

**Decision (5.3).** redistributed_fit is "better" than the incumbent on the Li grain (delta -0.027,
95 % [-0.045, -0.003]) and "not worse" on every other endpoint (Tang -0.001 [-0.084, +0.086], Kim point
+0.008 [-0.002, +0.038], Kim interval -0.001 [-0.002, 0.000], Li shoot +0.103 [-0.035, +0.191]; all
ties); no other variant is better anywhere. (The first run, before the numerics correction below, gave
the same decision: -0.029 [-0.046, -0.006].) The rule switches the default to **redistributed_fit**,
and the user applied it (2026-10-26). Recorded with it: the margin is small, the Li shoot point score
is worse (0.364 -> 0.467, inside the tie), and the chosen variant fits f_xy per congener to Yamazaki's
straw, so on Yamazaki it is a reproduction, not a prediction. The honesty band is its transfer median,
0.562 log10 (x3.65; amendment 1: 0.489, x3.1). `model_api` now resolves an unspecified loading to the
default shoot's (`resolve_loading`); the single-straw core has no `_fit` variant and keeps the lipid
loading (`SHOOT_DEFAULT_LOADING`).

Against amendment 1 the grain transfer is worse: Kim point 0.489 -> 0.554 (incumbent) and Li grain
1.019 -> 1.361. The brown-rice pool takes no xylem: what the panicle transpires stays in the husk, so
the brown-rice value is lower than amendment 1's whole panicle; on the rough basis Li grain is 0.864.

Section 7.8 (written for the incumbent default; read on its records):
- E17 HELD. Every r_eq moved by the closed-form factor within 0.001 (PFBA 1.129, PFPeA 1.097, PFHxA
  1.064, PFHpA 1.034, PFBS 1.032, PFOA 1.016, the rest <= 1.006).
- E18 HELD. 0.159 -> 0.154.
- E19 HELD. root.theta@DVS2 is third, after root.f_cw and the lignin share, for PFBA, PFPeA, PFHxA, PFBS.
- E20 WRONG. Every congener with a grain cell now fits L_Ph > 0 in every variant (amendment 1 left five at
  0), and the root does not move; but the straw does, by up to 0.086 log10. Not foreseen: the
  expectation held the lipid constants fixed, but they are shared by all chains and re-fit with the
  grain (a 0.0538 -> 0.0659, g_xy 0.00494 -> 0.00712), which moves every straw cell.
- E21 HELD in direction, short in size. husk + rachis / brown rice, dry: PFOA 2.71 (Choi 4.2-8.8; -0.187
  log10), PFOS 4.74 (37-47; -0.892). On redistributed_fit: PFOA 2.91 (-0.156), PFOS 3.31 (-1.048).
- E22 HELD. Li grain 1.361 brown, 0.864 rough.
- E23 WRONG. The crop sets enter the straw and grain band for 11/11, but they, not the panicle share at
  heading, lead the grain band (crop_set 6/11, phi 4/11, PFBS's L_Ph range 1/11; the panicle share 0/11).
  Not foreseen: the japonica sets move the heading-to-maturity partitioning that sets the brown rice's
  dry growth, hence its phloem draw, and with the deposit in the husk the panicle share no longer
  reaches the grain at all.
- E24 HELD (`tests/test_amendment2a.py`: relative change < 1e-6).

Band (`input_range_band`, median over 11 congeners, down / up): on the incumbent, root ÷1.29 / ×1.29,
straw ÷1.87 / ×2.44, grain ÷3.75 / ×6.28 (amendment 1: ÷1.27 / ×1.25, ÷1.35 / ×1.32, ÷2.54 / ×2.07); the
crop sets lead the straw band for 11/11. On the chosen default (redistributed_fit; the band the app
and `ASSUMPTIONS.md` carry): root ÷1.29 / ×1.29, straw ÷1.40 / ×1.89, grain ÷2.97 / ×3.65; the crop sets
lead the straw band 11/11 and the grain band 7/11 (phi 3/11, PFBS's L_Ph range 1/11), root composition
(f_cw 6/11, f_PL 5/11) the root band. The envelope covers the all-corner run in 38/66 cases (incumbent
39/66; amendment 1 53/66): with the crop sets in, the inputs interact more than the one-at-a-time
envelope assumes.

Numerics and recorded deviations (section 1: a number that moves as a numerical setting approaches
its limit is fixed in the model):
- **The grain was not converged in the organ-mass floor.** `numerics_convergence.py` on the records
  above: floor 1e-12 kg moved the grain by 0.0011 log10 (1e-15: up to 0.0023), growing with log(floor).
  Cause: brown rice is 0 at DVS 1 by construction (h = 1), and in the 2020-09-01 run the calibrated DVS
  at the heading node is 1 + 2.3e-9, so brown rice there is 2.2e-11 kg -- below the 1e-9 floor, above
  1e-12. The grain formation gate (a step at the floor) therefore opened a day earlier at a lower floor,
  and the first day of filling took a full day of inflow instead of half. Two changes, no new number:
  (i) brown rice's phloem draw is not gated (it is proportional to its own dry growth, so an organ not
  growing draws nothing; the gate stays on the panicle basis, whose xylem share needs it); (ii) the
  terminal pools (grain, husk) are integrated as burden y = M*C (`linear_expm.burden_affine`;
  rhs and affine stay in concentration form), so their growth dilution cancels exactly -- without the
  gate, the concentration form has a relaxation time M_floor/dM ~ 1e-15 d at heading that BDF cannot
  step. After both, floor 1e-15 moves no number by more than 1e-5, and brown rice's pre-harvest water
  moves its dry value by 0 (it was 0.0012 log10; `tests/test_amendment2a.py` now holds it to 1e-9).
  (iii) Without the gate, the central-difference growth assigned half of the first day's growth to the
  day before, when brown rice still sits on the floor: its burden there was right but its
  concentration (burden / floor) was ~1e5 x the harvest value at that one node
  (`test_grain_formation_gate`). The growth that draws a terminal pool's inflow (the grain's dM, its
  dry growth, the husk's growth) is now the slope of the linearly interpolated mass on each interval
  (`linear_expm.interval_slopes`); this moved every Yamazaki cell by <= 0.0002 log10, below reporting
  precision, so the re-fit was not repeated again. The transfer records and the brief were: the Li
  grain, scored on the filling-stage concentrations, moved by -0.003 in every variant (the table
  above is the repeat), and the decision did not change.
  At the fitted values the Yamazaki observable moved 0.1160 -> 0.1158 (redistributed_fit); the re-fit,
  the transfer runs, both bands and the brief were repeated on the corrected numerics, and the numbers
  in this section are those of the repeat (the first run's are in git history).
- `numerics_convergence.py` on the final state: the floor (1e-8 to 1e-12), the substeps (32, 64), the
  step grading (0.05), BDF at rtol 1e-8, the harvest grid (96/day), N (8, 16) and lam_eq (1e5, 1e7) move
  every reported number by < 0.0005 (largest 0.0002, BDF).
- The husk + rachis mass dips by up to ~5 % between Ansari anchors (7.11); kept.

## 9. Amendment 2b -- the remaining values without a source (committed before its runs)

The user's principle for this amendment (2026-10-27): what the tool runs on must be scientifically
correct and supported by a source; plans and procedure come second. Decision rules (§1, §2, 5.3, 5.4)
and the default variant (`redistributed_fit`) are unchanged. Every item below is a correction of a
value or form without a source and is adopted whatever the scores do (as 7.7). Sources are rows of
`docs/literature_db/sourced_inputs.csv` (ids given); papers not readable here are listed in
`fulltext_needed.md`. The literature for this section was read before this text was written; no
model run of 2b had been made.

### 9.1 Neutral path
- (a) **One plant for both paths.** The rice tool runs the neutral path on the forms of 7.1-7.2: organ
  water by development stage, brown rice plus the husk/rachis pool, the crop's DVS. `neutral_forms`
  stays for back-comparison only. The validation tables keep their own plants' composition (5.10); their
  drivers carry the crop's DVS (`growth_rice` provides it).
- (b) **Leaf/panicle split of the stem's xylem stream** = the sourced panicle share lam_grain(DVS) (5.2).
  The neutral path used the transpiring-surface split with `S_XYLEM_DRY` (leaf 20, grain 0.5116), a
  literal the code calls "NOT measured" and the provenance table did not list (new row `neu_xylem_split`).
- (c) **The brown rice is phloem-fed** (7.2: the panicle's transpired solute stays in the husk). With
  the phloem off (the neutral base) the neutral brown rice receives nothing; the neutral phloem-sap
  concentration (passive equilibration with the leaf's free water, as the DPU fruit models assume) is
  not yet read at source (Trapp 2007 SAR QSAR 18:367; Trapp 2011 SAR QSAR 22:191 -- requested). Until it
  is, the neutral BROWN RICE is reported as not predicted (NaN); rough grain and husk are reported.
  A phloem run (opt-in) uses the sourced phi and T_C_Ph midpoints (5.1); the unsourced defaults
  0.1 / 10 L/kg are removed (new row `neu_phloem_defaults`).
- (d) **Lipid-sorption exponent per organ.** Trapp, Matthies & McFarlane 1994 p.417: "0.77 for roots and
  stem [3,4], and 0.95 for leaves [5]" (Briggs 1982 macerated roots; Briggs 1983 macerated shoots). The
  ODE used 0.77 for every organ (the leaf's 0.95 reached only the equilibrium diagnostics). Root and stem
  0.77, leaf 0.95; the grain has no measured value and runs at the midpoint of the two measured cereal-
  tissue values, 0.77-0.95 (5.4). Rows `neu_leafexp`, `neu_leafexp2`.
- (e) **a = 1.22** is the octanol density correction 1/0.822 (Trapp & Matthies 1996, ES&T 30:360); the
  plant-lipid = octanol factor (x1) is the stated assumption of that source ("no values are known").
  Value unchanged; status derived (`neu_lipid_a`).
- (f) **Air-exchange surface.** Both faces of the green leaf (cuticle and boundary layer exchange on every
  face the air reaches; the permeabilities are per area of that face): 2 x the crop model's leaf area
  (ORYZA LAI x ground area per hill where the run records it; otherwise SLA(DVS) x green-leaf dry mass
  with the IR72.D92 SLA function), divided by the fresh leaf mass at each step. Stem: no sourced cuticle
  area, no exchange. Panicle: its exposed surface is the hull, i.e. the husk/rachis pool, which has no
  partition coefficient in the model, so no exchange (the sourced panicle areas -- Chang 2020, Meera
  2019, `neu_surface_grain*` -- are recorded for when a husk partition is sourced). The constants
  20 / 2 m2/kg are removed. Air exchange now runs with the staged water (B(t)).
- (g) **The K_AW threshold 1e-4 is removed.** The warning is raised for any K_AW > 0 with air off; the
  size of the sink is what `air=True` computes.
- (h) **Pore-water pH** 6.5-7.0 (Yoshida 1981 p.173; Moraes 1973), run at 6.75 (`pw_ph`, `soil_ph_*`).
- (i) No measured volatile-compound plant dataset was found; air exchange stays an opt-in sensitivity.

### 9.2 Soil: the HYDRUS paddy scenario (Tsukuba Andosol)
- Profile: plough layer 0-15 cm ("耕深15cmの確保", Ibaraki calendar), subsoil below. Bulk density: plough
  layer 0.61 g/cm3 (Koyama & Hayashi 2017, Tsukuba Andosol paddy topsoil); subsoil 0.70 / 0.57 (Tournebize
  2006 Table 3, the two Tsukuba Andisol layers) as range ends.
- Retention (van Genuchten, single porosity): the Tsukuba Andisol undisturbed cores of Tournebize 2006
  Table 3 -- (theta_s, alpha, n) = (0.709, 0.0050, 1.068) and (0.717, 0.0098, 1.156) as range ends, for
  both layers (transfer from the subsoil to the puddled plough layer stated); theta_r not printed there:
  0-0.030 (Rudiyanto 2013 fixed 0; Tokumoto 2010 fitted 0.030); Ks = the measured 11.1 / 93 cm/d of the
  same layers; l = 0.5 (Mualem's general value; no paddy value exists, Rudiyanto's fits are unstable).
- Water management: top boundary = the water input; under ponding the profile percolation is the sourced
  0.41-0.67 cm/d (Kondo 2019, Ibaraki Andosol paddy), imposed as the input excess over ET. No plough pan
  is modelled: no pan Ks is measured anywhere (Tournebize's is calibrated, in alluvial soil) and with the
  sourced percolation imposed it is not needed. Flooded from transplanting; midsummer drainage from
  35-40 DAT for 15-20 d (range ends; midpoints 37.5 / 17.5); flooded again to 30 d after heading
  ("出穂後30日頃まで水を入れて"), which at the Yamazaki site is harvest.
- Dispersivity 0.2-1 cm (Tokumoto 2005 columns; Tournebize 2006 field tracer 1 cm), midpoint 0.6.
- f_oc 0.037-0.038 (Yamazaki SI Table S5); Kd = Koc(QSPR) x f_oc as before.
- The root zone's exposure = the root-weighted mean pore water over the profile, replacing the single node
  at 10 cm (no source for a depth). The root profile (linear to 30 cm) has no source in the repo
  (new row `soil_root_beta`); it is checked against Beyrouty et al. 1997 if obtainable, else run as is
  and flagged.
- Profile depth, node spacing and the k_leach cap: convergence (§1).
- The analytic `flooded` shape's per-congener k_leach table (`params/cwo_kleach.csv`) is recalibrated on
  this scenario.

### 9.3 The irrigation-source scenario (`obs_exposure_time`)
The Yamazaki field's PFAS arrive with irrigation water and rain onto a soil that already holds them.
Scenario: initial pore water from the measured pre-season soil (Yamazaki SI Table S16, water extract +
methanol extract after water, pg/g dw, 0-15 cm) at equilibrium with the scenario's Kd; the top input =
irrigation water (S16) and rain (S16 R1-R5 mean, wet-deposition volumes from the site's JMA weather),
flux-weighted. Output: root-zone pore water / irrigation concentration over the season, per congener.
- **Out-of-sample check** (not used to set anything): the predicted post-season 0-15 cm total soil
  concentration against the measured post-season soil (S16 "After cultivation"), per congener, log10
  ratio recorded. Single composite samples, so no statistical decision; recorded.
- The default observable stays on the irrigation grab held constant (the BAF's own denominator). The
  plant model is also run on the scenario's Cw(t) (normalised to the grab, the `cwo_profile` convention)
  and the change of each Yamazaki cell is recorded; under 5.4 it enters the input-range band when it
  moves a number by 0.0005 log10 or more.

### 9.4 Analysis rows
- `fxy_recommended`: **untraceable.** Its coefficients (logit 4.061 - 0.857 n) cannot be recovered from any
  recorded point set nor from the shipped table (4.377 / -0.913), and the cited absolute anchor contradicts
  it: Felizeter 2012/2014 solution-normalised TSCF is roughly flat at 0.05-0.23 from C4 to C11 in cabbage,
  zucchini and tomato (lettuce U-shaped), against the curve's 15-35x lower PFDA. f_xy is not a TSCF, so
  Felizeter cannot replace it either. The `*_free` configurations that read it are withdrawn from scoring
  and decisions (kept runnable, labelled). The Felizeter points are recorded (`raw_si/felizeter201*_tscf.csv`).
- `tr_drivers`: no dataset reports its own site calendar and weather in the main text (Tang: no dates;
  Kim: anonymous regions, harvest window only; Li: coordinates and year, no stage dates). The SI holding
  them is requested; until read the transfers stay on the Yamazaki drivers, documented as a limitation.
- `tr_li_grain_from`: grain scored from heading vs from filling, both recorded; reported if they differ.
- `calib_sigma`: the default 0.3 is removed from the synthetic-recovery machinery; sigma is required.

### 9.5 App
- `app_density` (per-volume reporting): removed; it entered no result and its values had no source.
- `smiles_kcw` and the K_cw of the curated congeners (`bind_kcw`: the polysaccharide scale has no
  traceable source, the lignin value is a Kraft lignin): the K_cw side project (`docs/kcw_side_project/`)
  decides them; they are recorded as open here and will be a later addendum.

### 9.6 Runs
- PFAS default: no input of `simulate()` changes in 9.1-9.5 (pH reaches f_d only; f_d at pKa 0.5 moves
  in the 7th decimal). No re-fit. Every scored Yamazaki/transfer endpoint is re-run to show it.
- Re-run: `validation/neutral_scores.py` (old vs new), the HYDRUS scenario per congener, the k_leach
  recalibration, the irrigation-source scenario and its post-season soil check, numerics convergence of
  the soil inputs, and the input-range band.

### 9.7 Expectations (written before the runs; misses will be recorded)
- E25 Every PFAS `simulate()` output of the default variant is unchanged to 4 decimals in log10, and the
  Yamazaki observable is unchanged.
- E26 The equilibrium-basis neutral root scores do not change (3 decimals); the season-ODE root scores move
  by less than 0.02.
- E27 The Ge 2017 per-organ score changes (leaf exponent, staged water); the direction is not predicted.
- E28 The neutral leaf air-exchange surface at heading is 14-17 m2/kg fresh (2 x SLA 24.5 x (1 - theta_leaf)).
- E29 Retardation falls (rho_b/theta_s 0.61/0.71 against 1.3/0.46) and percolation rises (0.54 against
  0.30), so every congener's pore water leaches faster than in the old scenario: the recalibrated k_leach
  rises for every congener with a non-zero value.
- E30 In the irrigation-source scenario the root-zone pore water of C4-C6 PFCAs tracks the irrigation
  concentration within the first flooding period, while C >= 10 stays set by the soil.
- E31 The post-season soil check: the model's short-chain post/pre ratio is below 1 as measured (PFBA
  water extract 82 -> 40 pg/g), and its PFOA-PFDA ratios are within 0.3 log10 of 1 as measured.
  (Selection note: the measured post-season values were seen when Table S16 was read for the initial
  condition, before this was written; E31 is therefore not a blind expectation.)

### 9.8 Addendum (committed before any soil run; found while building 9.2)
9.2 dropped the plough pan. That was wrong: a flooded paddy's plough layer is saturated under standing
water, and without a restricting layer the imposed percolation drains the Andosol at unit gradient, where
the measured retention curves leave it at theta/theta_s = 0.998 (Tournebize layer 1) to 0.89 (layer 2) --
up to 11 % below saturation, which enters retardation and pore velocity. The pan returns, with nothing
chosen that the data do not identify:
- below the 15-cm plough layer; thickness 2 cm (Tournebize 2006's; only Ks/thickness is identifiable from
  percolation, so the thickness itself does not change the flux);
- its Ks is solved so that steady percolation under the ponded depth equals the sourced percolation
  (0.41-0.67 cm/d, midpoint 0.54), with the plough layer saturated and the subsoil at unit gradient;
  ponded depth 2-3 cm ("浅水(2～3cm)", Ibaraki calendar), midpoint 2.5;
- under ponding the top input is then the water the field takes (ET + percolation) and the pan, not an
  imposed flux, sets the percolation; during midsummer drainage the input is rain only.

### 9.9 Addendum (committed before any soil run; found while implementing 9.2-9.3)
- **The subsoil's initial PFAS is not measured** (Yamazaki sampled 0-15 cm). Its pore water runs from 0 to
  the plough layer's value as range ends, at the midpoint (half); both ends are reported.
- **Flooded water balance**: the field takes ET + percolation each day, rain first and irrigation for the
  rest; rain beyond that leaves over the bund (not stored as ponding). While drained, all rain enters.
- **Rain concentration** = the mean of Yamazaki's five rain samples (R1-R5); a value below its MQL is read as
  half the MQL and flagged (none of the 12 model compounds is affected in the irrigation or soil columns).
- The scenario needs the crop's transpiration and evaporation; there is no default (the old gaussian
  transpiration and the flat 0.10 cm/d evaporation had no source). Generic forcing: `forcing_rice`
  transpiration and its evaporation complement (1 - f_T) Kc ET0 (Nay Htoon 2018 split); the site: the
  ORYZA2000 Tp and Ep of the Yamazaki run, JMA Tsukuba 2020 daily rain (`raw_si/jma_tsukuba_2020_precip.csv`).
- `top_bc=3` of the old scenario was surface RUN-OFF, so the old "flooded" paddy never ponded; the
  scenario now uses the surface layer (`top_bc=2`).

### 9.10 Addendum (committed before any scored soil run; found in the first test runs)
- **Root water uptake.** The old scenario kept HYDRUS's generic Feddes reduction (p0 = -10 cm), which stops
  uptake in saturated soil -- an upland crop's anaerobiosis. It never acted there because the old top boundary
  ran water off and the soil never saturated; with ponding it halves the transpiration and the ponded depth
  climbs to 17 cm. The crop side of the model is ORYZA2000 in potential production (no water stress), and the
  transpiration driving the soil run is that crop's (site) or one measured in flooded paddies (generic), so a
  reduction in HYDRUS would contradict its own driver. No reduction above the wilting point: uptake at the
  potential rate from ponding (p0 = popt above any ponded head) down to -15000 cm (-1.5 MPa), where it stops.
  No rice-specific Feddes set was found (searched: HYDRUS database, SWAP crop files).
- **Grid.** Under the steady state the pan's base is unsaturated over 0.04-0.37 cm (`pan_ks`), so the node
  spacing there is dx/100 (graded grid, `_graded_profile`); convergence (§1) covers it.
- **Pan solve.** The pan's Ks is solved with the pan's own K(h) (its base desaturates), not as a saturated
  layer: the saturated form under-estimated the resistance and gave 0.21 instead of 0.54 cm/d. Result:
  Ks_pan 0.046-0.056 cm/d over 2 cm (the order of Tournebize 2006's independently calibrated 0.08).

### 9.11 Addendum (after the registered k_leach recalibration; before anything is built on it)
The registered recalibration ran (E29 holds: k_leach rose for every congener). It also showed that the
analytic `flooded` shape cannot represent the new scenario: 9 of 13 congeners end at the scan ceiling 0.15
(`soil_kleach_max`), i.e. the result keeps moving with a numerical bound, which §1 calls structural. Cause,
checked in the HYDRUS profiles: clean water drives a leaching FRONT down from the surface (PFNA, Kd 12.8: top
5 cm flushed, 10 cm at 0.64, >= 15 cm untouched at harvest; solute balance error 0.015 %), and the root zone is
weighted to the top; a well-mixed reservoir with first-order loss of its dissolved fraction cannot take that
shape for any k_leach. Correction, no new number: the `flooded` profile and the soil-inventory mode read the
scenario's own HYDRUS shapes, stored per congener as c(t)/c(0) on the site's season (`params/cwo_shapes_hydrus.csv`,
built by `validation/cwo_shape_table.py`); a structure without a stored shape (SMILES, neutral) interpolates in
log10 Kd between the stored ones. `k_leach`, its per-congener table and its cap are removed. Outside the site's
season the shape is held at its end values (stated).

### 9.12 Addendum (soil numerics; before the shape table and the irrigation check are rebuilt on it)
The §1 convergence check of the soil scenario failed as first built: the root-zone shape moved with the
node spacing (mean c/c0 0.456-0.480 for PFOA over dx 0.5-0.1, non-monotone) and so did the water flow
(season percolation 0.44-0.53 cm/d, ponded depth at harvest 1.6-6.3 cm). Three causes, each found in the
HYDRUS output and fixed without a new number:
- **The grid refined the saturated plough layer.** Its elements are highly conductive (Ks 11-93 cm/d), and
  HYDRUS's head iteration (TolH 1 cm) stops with a large flux residual there: water balance errors grew to
  93 % on a 0.01-cm grid. The grid is now refined only inside the pan and the top of the subsoil (where the
  steady head climbs ~250 cm within ~0.2 cm, `_pan_profile`); the plough layer keeps `dx`.
- **The pan's thickness depended on the grid.** HYDRUS averages K arithmetically between nodes, so the
  element joining the plough layer to the pan carries no pan resistance; with the pan's material starting
  one node below 15 cm, the pan was thinner by that element (an O(dx) error, 2.4 % of the percolation at
  dx/5 = 0.05 cm). The pan's material now spans the nodes at 15 and 17 cm exactly, as `pan_ks` assumes.
- **The hydraulic function near saturation.** The Tsukuba Andisol retention has n = 1.068-1.156; for
  1.0 < n < 1.3 the van Genuchten-Mualem conductivity falls steeply just below saturation (K(-2 cm) is 14 %
  of Ks at the midpoint set), which the HYDRUS-1D manual (4.17, pp.23-24) says can change the
  conductivity function "considerably" and for which it recommends its model with an air-entry value of
  -2 cm (Vogel & Cislerova 1988 form, theta_m extrapolated; row `soil_hentry`). The scenario uses that model,
  and `pan_ks` solves the pan with the same function (Ks_pan 0.0256 cm/d at the midpoints, h below the pan
  -256 cm). The primary source (Vogel et al. 2001) is closed access and was not read.
Result (PFBA, PFOA, PFDoDA): at dx 0.25 cm (pan base 0.005 cm, dx_deep 1 cm, 100 cm, TolH 1 cm, Crank-
Nicolson) every check -- dx 0.5 / 0.125, pan 0.01 / 0.002, the near-pan spacing, dx_deep 0.5, 150 cm, TolH
0.1 -- moves the season-mean shape by <= 0.0004 log10 (§1: 0.0005); water balance errors 0.003-0.07 %. The
fully implicit time weighting moves PFBA by 0.0023 log10; Crank-Nicolson (second order) is kept. The shape
table (9.11) and the irrigation-source check (9.3) are rebuilt on these settings; their earlier numbers
(dx 0.5, upstream weighting, fully implicit) are superseded.

### 9.13 Addendum (water after midsummer drainage; found re-reading the calendar for 9.12)
9.2 read the Ibaraki calendar as "flooded again to 30 d after heading". The calendar says the field is kept in
INTERMITTENT irrigation after midsummer drainage ("中干し後は、根の活力を維持するため、引き続き間断かんがいを
行う"; "出穂後30日までの間断かんがい"; "出穂後30日頃まで水を入れて"), and Tournebize 2006 (the NIAES Tsukuba paddy)
the same ("intermittent irrigation was performed until 1 week before the harvest"; 0-5 cm ponded water depth).
Neither gives the wetting/drying cycle, and none was found. What is measured is the water passing the profile
in that stage: Tournebize's stage-3 percolation 0.8 (1998) and 0.3 (1999) cm/d, which brackets the flooded
0.41-0.67 cm/d used. The scenario therefore keeps the flooded water balance after midsummer drainage (the field
takes ET + percolation; the solute leaves with the percolation) as its base, and runs as the other end a field
that receives no irrigation after midsummer drainage (rain only), reported beside the base in the irrigation-
source check. No cycle is chosen. (Calendar re-read 2026-10-27: sha256 63e299e9...)

### 9.14 Outcome of amendment 2b (records named; nothing fitted)
**Expectations.**
- E25 **held**: every PFAS `simulate()` output of the default and the Yamazaki observable are bit-identical
  (137 cells, max |Δlog10| 0), re-checked after the last code change of 2b.
- E26 **held**: the neutral root scores did not move (equilibrium and season ODE, 3 decimals).
- E27 **held** (direction not predicted): Ge 2017 per-organ 0.631 → 0.606 (`params/neutral_scores.json`).
- E28 **held**: the leaf air-exchange surface at heading is 14.1 m2/kg fresh (season 6.1-17.2).
- E29 **held** for the registered recalibration (every k_leach rose), but that surrogate was then removed (9.11).
- E30 **missed**: at the end of the first flooding period the C4-C6 root-zone pore water is still 1.7-2.0 x the
  irrigation concentration (it starts at the soil's own pore water, 2.8-4.8 x at day 15), not tracking it.
- E31 **partly missed**: the short chains fall over the season as measured (predicted post/pre -0.59 to -0.91
  log10 against -0.03 to -0.30 measured, so too far), PFNA -0.25 and PFDA -0.07 are within 0.3 of 1, PFOA
  (-0.60) is not. Post-season soil, predicted/measured: C4-C8 -0.46 to -0.62 log10 (over-leached), PFNA-PFDoDA
  and PFOS -0.15 to +0.01; the subsoil end (0 / half / same) moves none of them by more than 0.01; the
  rain-only end after midsummer drainage (9.13) moves C5-C7 and PFBS further down (to -1.39 for PFBS) and
  the rest by <= 0.1 (`validation/irrigation_source_record.json`).

**Deviations from the registration, with the reason.**
- The irrigation-source scenario is **not merged into the input-range band** (9.3 said it would enter it
  when it moves a cell by >= 0.0005). It fails its own out-of-sample check for C4-C8: it leaches those
  chains several-fold faster than the field shows, so the Koc-QSPR Kd is too low for this Andosol (or a
  sorption process is missing) and its pore-water course for them is not supported. Its plant deltas are
  recorded (root/straw/grain +0.00 to +0.27 log10 for C4-C8, +0.47 to +1.18 for C9-C12 and PFOS, which
  pass the soil check) and `obs_exposure_time` stays open.
- The k_leach surrogate, its table and its cap were removed (9.11), not recalibrated.
- The neutral **rough grain** is not predicted with the phloem off either (9.1c reported it): it carries the
  brown rice's burden, which is unknown, so with the brown rice at 0 it was only a lower bound. The husk +
  rachis pool is reported.
- Downstream products rebuilt on 2b: the biocide census and its briefing now report the husk + rachis pool
  (dry basis) where they reported a neutral "grain", and say that the brown rice is not predicted; the
  half-equilibrium point moved 6.26 -> 6.40 (root water at harvest 0.6365); the rank correlations are
  unchanged (+0.904 / -0.777). Two numbers that had already disagreed with the census record before 2b
  (propiconazole 1.03x, triclosan +0.624) were corrected to the record (1.24x, -0.197). In the PFAS brief the
  flooded-schedule grain is now per unit season-mean pore water (the harvest-day pore water of a short chain
  leaches to ~0, so a BAF on it is meaningless).
- `tr_li_grain_from` is not invariant (brown rice from heading 3.33 vs from filling 1.33 log10 on the default;
  the variant ranking is the same): filling stays, both are recorded (`validation/li_grain_from_record.json`).

**What 2b leaves open**: `obs_exposure_time` (default scope), `soil_root_beta`, `tr_drivers`, `smiles_kcw` and
the curated K_cw (side project), the neutral phloem-sap concentration (Trapp 2007, 2011 requested).
The input-range band is unchanged (no input of `simulate()` changed, E25).

### 9.15 Addendum (the neutral brown rice; committed before its runs)
9.1c left the neutral brown rice unpredicted because the neutral phloem-sap concentration had not been read at
source. Trapp 2007 and Rein, Legind & Trapp 2011 are now read (rows `neu_phloem_xylem_eq`, `neu_phloem_fruit_flow`,
`neu_phloem_downward`):
- Trapp 2007 p.374: "no difference is made between chemical concentration in phloem and in xylem ... For the less
  polar neutral compounds, xylem and phloem concentrations will equilibrate within a short distance".
- Rein 2011 p.194, eqs 4-5: the phloem flux to the fruit, T_C,Ph x the fruit's dry growth, is taken "under the
  assumption of phase equilibrium to xylem flux": subtracted from the flow to the leaves and added to the flow
  to the fruit.
- Rein 2011 p.208: downward phloem transport is irrelevant for the less polar neutrals (it matters for weak
  acids, by the ion trap, and for very polar neutrals, which leave the sieve tubes slowly).
Change (neutral path only): the brown rice is fed by a phloem stream at the stem's xylem-sap concentration, of
volume T_C,Ph (the rice value, 6.1 L/kg, 4.4-7.8, PREREG 5.1) x its dry growth, taken out of the leaf's xylem
inflow (never more than that inflow); no downward stream (phi = 0, as before); the husk keeps the panicle's
transpiration deposit. The brown rice is then predicted by default, and the rough grain with it. `phloem=True`
(the leaf-source stream with the weak-acid trap) stays as it was. No log Kow below which a neutral counts as
"very polar" is given in either source, so none is set: the equality is reported as the sources' assumption and
its polar-end limit is stated in the outputs. The PFAS path is untouched (the core's default is unchanged).
Expectations:
- E32 Every PFAS `simulate()` output of the default and the Yamazaki observable are bit-identical.
- E33 No neutral validation table has a grain row, so no neutral score changes by more than 0.01 log10; the root
  scores do not change at 3 decimals.
- E34 The neutral brown rice (dry) is below the husk + rachis pool (dry) for every log Kow in the biocide census:
  per kg dry, the brown rice draws 6.1 L of sap, the hull the panicle's transpiration share.
- E35 The census rank correlations (root +0.904, straw -0.777) do not change at 3 decimals.

### 9.16 Outcome of 9.15
- E32 **held**: the PFAS default and the Yamazaki observable are bit-identical (137 cells, max |Δlog10| 0).
- E33 **missed**: Ge 2017 per-organ moved 0.606 → 0.593 (0.013 > 0.01). The grain's stream is taken out of the
  leaf's xylem inflow, so the leaf's TF falls. The root scores did not move (no downward stream).
- E34 **held**: brown rice (dry) < husk + rachis (dry) on all 117 census rows (triclosan 0.119 against 1.64).
- E35 **held**: census rank correlations +0.904 / −0.777.
- The two solvers agree on the new stream (exponential integrator against BDF, ≤ 5e-6 log10; tested).
- 9.14's neutral outcome is superseded where it says the brown and rough grain are not predicted: they are, on
  the sources' assumption (less polar neutrals; no polar-end log Kow is given). The biocide census and briefing
  report the brown rice again, beside the husk + rachis pool (dry); the fish-kM sensitivity moves the brown
  rice by up to 8,344×. `grain_feed=None` keeps the unpredicted reading.

### 9.17 Record (the requested SIs; no model change)
The Tang 2026, Kim 2019 and Felizeter 2012/2014 SIs were read (uploaded by the user). `tr_drivers` cannot be resolved from
the papers: neither Tang's SI nor Kim's gives a sowing, transplanting, heading or harvest date, and Kim's regions are
anonymised without coordinates. The transfers stay on the calibration site-year's drivers, now as a limitation of the data,
not a pending request (row `tr_drivers` closed). Felizeter's SIs give the experiment dates and tissue tables but no TSCF
table; TSCF (their eq. 4) uses a modelled solution concentration, so the figure-read TSCFs stand. The repo's earlier
transcriptions were checked against the SIs: Tang S7/S8, 75 cells, one SD corrected (PFOS 0.1 ug/g BCF 0.025 -> 0.0254; no
code reads it); Kim S5/S6, 51 rows, all equal. New transcriptions: `raw_si/tang2026_si_table_s2_soil.csv`,
`raw_si/kim2019_si_table_s2_sites.csv`, `raw_si/felizeter2014_si_table_s3_dates.csv`.
