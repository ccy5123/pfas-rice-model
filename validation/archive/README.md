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
The current transport fit is `validation/refit_site_year.py`.

Notable groups:

| files | what it was |
|---|---|
| `twopool_*`, `longchain_*`, `pfas_rice_two_pool.py`, `twopool_*fitted_params*.json` | the two-pool root explorations (sequestration and carrier). On fresh organ mass the single-pool lipid model already beat both, and the two-pool's distinguishing evidence (the non-K_PL PFOS/PFUnDA separation) did not survive; they were not re-fitted on the corrected scenario. |
| `S6_*`, `W2_transport_fit.csv`, `apriori_prediction.py`, `absolute_fxy_fit*` | the original W2 transport fit on the placeholder logistic and its a-priori comparisons. |
| `refit_oryza.py`, `refit_oryza*.csv`, `structural_adequacy_fit.py`, `refit_realistic_biomass*`, `root_shoot_biomass_sensitivity.py`, `biomass_shape_compare.py`, `oryza_growth_validation.py`, `mass_drivers_plot.py` | fits and driver comparisons on the reduced ORYZA / growth_rice drivers. |
| `fit_lipid_loading.py`, `lipid_loading_fit.json` | the lipid-loading fit on growth_rice drivers (now part of `refit_site_year.py`). |
| `oos_tang.py`, `tang2026_validation.py`, `tang2026_nstem_validation.py`, `tang2026_fxy_TF_validation.py`, `validation_summary.py`, `revalidation_crosscheck.py` | out-of-sample reports and summary figures built on the superseded fits. |
