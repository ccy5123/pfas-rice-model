# Archived analyses (not maintained, not current results)

These scripts, caches and fitted tables were computed on drivers, an observation operator or
a mass basis that were later found to be wrong:

- **drivers**: a placeholder logistic or a 120-day climatology with transpiration borrowed
  from other sites (`forcing_rice`, about half the Yamazaki site-year's seasonal total) and
  the reduced ORYZA re-coding (`oryza_growth`), instead of the calibration site-year
  (`src/site_scenario.py`: JMA Tsukuba 2020, the ORYZA2000 port on the paper's dates,
  ORYZA2000 transpiration);
- **observation operator**: the model's GREEN-leaf straw compared with a harvested straw that
  holds dead blades, and fresh-weight observations built from ASSUMED water contents
  (now `src/harvest.py` and the dry-basis `data_obs/obs_baf_Yamazaki.csv`);
- **mass basis**: several were fitted while crop-model dry matter was used as the ODE's fresh
  organ mass (fixed by `src/organ_mass.py`).

They are kept as the record of how the model developed (CLAUDE.md section 6 cites them), not
as results. They are not run by the tests and may not run against the current API
(`mass_basis`, the W2 / ORYZA-reduced transport fields and the two-pool wrappers were removed).
The current transport fit is `validation/refit_sourced.py`.

Notable groups:

| files | what it was |
|---|---|
| `twopool_*`, `longchain_*`, `pfas_rice_two_pool.py`, `twopool_*fitted_params*.json` | the two-pool root explorations (sequestration and carrier). On fresh organ mass the single-pool lipid model already beat both, and the two-pool's distinguishing evidence (the non-K_PL PFOS/PFUnDA separation) did not survive; they were not re-fitted on the corrected scenario. |
| `S6_*`, `W2_transport_fit.csv`, `apriori_prediction.py`, `absolute_fxy_fit*` | the original W2 transport fit on the placeholder logistic and its a-priori comparisons. |
| `refit_oryza.py`, `refit_oryza*.csv`, `structural_adequacy_fit.py`, `refit_realistic_biomass*`, `root_shoot_biomass_sensitivity.py`, `biomass_shape_compare.py`, `oryza_growth_validation.py`, `mass_drivers_plot.py` | fits and driver comparisons on the reduced ORYZA / growth_rice drivers. |
| `fit_lipid_loading.py`, `lipid_loading_fit.json` | the lipid-loading fit on growth_rice drivers (later part of `refit_site_year.py`, now of `refit_sourced.py`). |
| `oos_tang.py`, `tang2026_validation.py`, `tang2026_nstem_validation.py`, `tang2026_fxy_TF_validation.py`, `validation_summary.py`, `revalidation_crosscheck.py` | out-of-sample reports and summary figures built on the superseded fits. |
| `gate_numerics_sensitivity.py` | the first look at the grain gate's 1e-4 kg floor and 1.5x ramp (PFOA young-panicle peak 1.41-9.65x). Superseded by `validation/numerics_convergence.py`, which showed neither value was converged; the model now uses one 1e-9 kg floor and a step gate (2026-10-25). |

## Retired by the no-unsourced-value amendment (2026-10-25)

`docs/PREREG_sourced_values.md` (amendment 1) replaced every model value that had no source.
These analyses ran on the replaced structure and are kept as the record of the decisions they
made; their numbers are not current results.

| files | what it was | why retired |
|---|---|---|
| `refit_site_year.py`, `refit_site_year.json`, `refit_site_year.csv` | the pre-registered transport fit on the site-year (Fit A/B, lipid globals; amendments 1-2) | fitted kappa_d with the GHK + Michaelis-Menten root (Km, Vmax without a source), L_Ph bounded at 1 (no source), the 4-pool's unsourced transpiring surfaces and shoot shares; replaced by `validation/refit_sourced.py` |
| `carrier_vs_bypass.py`, `dose_series_carrier.py` | carrier vs apoplastic bypass; the Tang dose series against the fitted carrier | the carrier and bypass modes left the model (their constants have no source; the root entry is now linear) |
| `default_switch.py` | the pre-registered switch to the redistributed shoot + lipid loading | its candidates ran with retention 0.6, lam_grain 0.05 and stem share 0.45 (no source); the default is re-decided under PREREG 5.3 |
| `panicle_share_amendment.py` | "gate" vs "size" routing of the residual xylem stream | with retention 1 there is no residual stream to route (PREREG 5.2) |
| `tang2026_fxy_refit.py` | an in-sample f_xy refit to Tang | on the replaced structure; the API override it fed (`TANG_REFIT_FXY`) was removed |
| `nstem_gradient_check.py` | the N-segment stem gradient against Yamazaki's height bands | ran with the carrier constants |
| `default_sensitivity.py` | the default's unmeasured inputs, one at a time against the two-harvest yardstick | its inputs now have sourced ranges; `validation/input_range_band.py` carries them into the output (PREREG 5.4) |
| `oos_crossdataset.py`, `oos_multidataset.py` | the Kim / Tang / Li transfer summary of the 4-pool variants | compared Kim on the fresh basis (Kim freeze-dried its grain: PREREG 5.7) and ran a 4-pool `fit` variant that no longer exists; `validation/transfer_sourced.py` scores every re-fitted variant on every transfer endpoint |
| `pfas_rice_plant_module_nstem.py`, `test_nstem.py` | N serial stem segments with the grain at the end of the stem stream (a multi-height mixer) and its kinetic variant | used only by the archived `nstem_gradient_check.py`; its phloem defaults (phi 0.1, T_C_Ph 10) had no source. The redistributed shoot (`pfas_rice_plant_module_nstem_leaf.py`) is the N-segment model the tool runs |
| `theta_fw_sensitivity.py` | organ water theta_fw one organ at a time; fixed 0.05/0.10 decision margins | its ranges came from the retired `rice_tissue_params.csv` (no source) and its margins were fixed; the sourced theta ranges are in `validation/input_range_band.py`. The model properties it found (shoot theta cancels on the dry basis, root theta is the sensitive one, the panicle share) are pinned in `tests/test_theta_fw_sensitivity.py` |
| `li2019_rcf_apriori.py` | the Li 2019 hydroponic table and the Briggs-anchor diagnosis (the rice root's 1 % fw lipid 2.5x below Briggs' fitted term) | the 1 % was a fresh fraction reported at 0.90 water run at another water; PREREG 5.10 carries lipid on a dry basis (the rice root now runs 2.40 % fw, near the anchor) and runs Li's table on Li's own water and species lipid. Its score comes from `neutral_dpu_validation.compare_to_obs` (`validation/neutral_scores.py`) |
