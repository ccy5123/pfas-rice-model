# Assumption register

Every input the shipped model runs on, with its source. The rule (2026-10-25): no value may be
chosen without a source (`docs/PREREG_sourced_values.md`). The machine-readable inventory is
`params/value_provenance.csv` (guarded by `tests/test_value_provenance.py`); the quotes behind
each sourced value are in `docs/literature_db/sourced_inputs.csv`.

How a value is admitted (PREREG section 1): a value read in a source; else a fit the data identify
(scored on data not used to fit it); else the mechanism that needs it is removed or reduced;
else the model runs at the midpoint of the sourced range and reports the range. Numerical values
(floors, tolerances, grids) are admitted only where the results converge at reporting precision
(0.0005 log10; `validation/numerics_convergence.py`).

**Effects of the ranges** are not typed here: they are in `params/input_range_band.json`
(`validation/input_range_band.py`, PREREG 5.4) -- per congener and tissue the log10 envelope
when each range input moves to its ends with the fits held, and the inputs that contribute most.

## 1. Calibration site-year (Yamazaki et al. 2023, NARO Tsukuba 2020) -- `src/site_scenario.py`

| input | value | status | source |
|---|---|---|---|
| weather | JMA Tsukuba 2020 daily + hourly | measured | `raw_si/jma_tsukuba_2020.csv` |
| transplanting | 2020-05-11 | measured | Yamazaki SI p.S2-S3 |
| harvests | Kasalath 2020-09-01, Koshihikari 2020-09-10 | measured | Yamazaki SI p.S3 |
| heading | Kasalath 2020-08-03, Koshihikari 2020-08-11 (mean run: their mean) | measured | Yamazaki SI p.S3 |
| harvest = physiological maturity (DVS 2) | -- | convention | the phenology is calibrated so that DVS 2 falls on the harvest day |
| phenology calibration | one factor on DVRJ/DVRI/DVRP -> heading; DVRR -> harvest | fitted | ORYZA2000 manual section 8.2 (DRATES) |
| crop parameters | IR72.D92 verbatim | **open** | no ORYZA crop set for Koshihikari or Kasalath has been published (`sourced_inputs` crop_cultivar_search); the japonica sets of Gao et al. 2021 and Wu et al. 2022 are the next amendment's range |
| nursery after emergence (SBDUR) | 25.5 d (range 25-26) | derived | sown 04-13 (Yamazaki), emergence 2-3 d after sowing (Ibaraki calendar); ORYZA counts whole days, so the midpoint runs the 25- and 26-d drivers interpolated |
| seedbed cover (TMPSB, on the daily maximum) | 2.5 C (range 0-5.0) | derived | Ibaraki calendar seedbed 15-25 / 10-25 C against an ambient maximum of 20.0 C (JMA Tsukuba 2020) |
| planting density | 18.5 hills m-2 | measured | Ibaraki prefecture standard (the field's own is unreported) |
| plants per hill | 4.5 (range 4-5) | measured | Ibaraki prefecture standard |
| CO2 | 414 ppm | measured | Mauna Loa 2020 annual mean |
| reference ET terms | FAO-56 Penman-Monteith radiation / aerodynamic terms | inherited | the definition the ORYZA2000 manual cites |
| transpiration | ORYZA eq. 4.4, actual = potential | inherited | continuous flooding except a midsummer drainage (Yamazaki SI p.S3) |

## 2. Observation (Yamazaki 2023) -- `data_obs/obs_baf_Yamazaki.csv`

Built by `docs/literature_db/tools/yamazaki_obs_from_si.py` from the SI transcription
`raw_si/yamazaki2023_si_s16_s18_s19.csv`; where the SI leaves a value open the table carries the
midpoint and both ends (`lo`, `hi`).

| item | value | status | source |
|---|---|---|---|
| tissue basis | dry weight (pg/g dw) | measured | Tables S18, S19 |
| exposure | irrigation water, Table S16 (one grab, 2020-06-08), held over the season | measured value; **time course open** | the seasonal course of the water is not reported |
| two cultivars | geometric mean (observation and model alike) | convention | each cultivar runs on its own heading and harvest |
| straw | height bands weighted by their share of straw dry mass; mean of the 10 and 49 d-after-heading profiles (range = each profile) | measured | Tanaka et al. 2022, Rice 15:10, Fig. 3 (digitised; `raw_si/tanaka2022_fig3_layers.csv`) |
| brown rice | white share 0.88 (range 0.86-0.90) of white + bran | measured | SI p.S3 |
| < MQL cells | midpoint of [0, MQL) with both ends | derived | Table S6 MQLs; among the used cells only PFOS indica bran |
| PFBS grain | not scored (white and bran < MQL in both cultivars); its censoring bound sets PFBS's L_Ph range | -- | |
| dead blades in the straw | all attached | measured (qualitative) | rice and most grasses do not shed dead blades (McGlone et al. 2014; Chang & Bardenas 1965); no fraction is published -- all-lost is a reported sensitivity |

## 3. Plant model

| constant | value | status | source / note |
|---|---|---|---|
| organ composition (theta_fw, f_prot, f_PL, f_cw, lignin share) | midpoints of `params/rice_composition_sourced.csv` | measured; physical bounds where no rice source exists (root and leaf f_cw upper ends; lignin shares of root, leaf and grain) | `sourced_inputs` tis_* rows |
| root theta_fw over the season | 0.713 (one maturity study) | **open** | later search: maturity pot roots 0.51-0.77 (five studies), seedling to vegetative 0.76-0.89 (hydroponic and flooded soil). The model's theta is constant over the season; how to represent the seasonal fall (a range, or theta(t)) is amendment 2 |
| K_PL, K_prot | per congener | measured | Chen 2025 (SSLM), Zhou 2025 |
| K_cw | per organ: (1 - s) K_poly + s K_lignin | derived | GAP A deliverable (Guo 2025, Mel 2024); s = the organ's lignin share |
| root entry | linear equilibrium `lam_eq (r_eq Cwo - Cw_root)` | structural | the GHK + Michaelis-Menten pair had no sourced constants; uptake is linear in exposure in every series found (Blaine 2014 and others) |
| r_eq | per congener and variant | fitted | Yamazaki root (closed form); per-cultivar interval recorded |
| lam_eq | 1e6 per day per kg | numerical | converged (`numerics_convergence.py`; 1e4 was not: 0.0008 log10) |
| L_Ph | per congener and variant, >= 0, no upper bound | fitted | Yamazaki grain; 0 where the xylem alone already over-feeds the grain; PFBS: midpoint of [0, L at the censoring bound] |
| lipid loading constants | per shoot model; terms not identified are removed | fitted | Yamazaki straw (and grain at L_Ph = 0), backward elimination by the bootstrap rule (`validation/refit_sourced.py`) |
| free-anion f_xy (`*_free` variants only) | monotone logistic | **untraceable**: its coefficients cannot be recovered from any recorded point set, and Felizeter 2012/2014's TSCFs (flat 0.05-0.23 from C4 to C11) contradict it; the `*_free` variants are withdrawn from scoring and decisions (PREREG 9.4) | `raw_si/felizeter201*_tscf.csv` |
| phi (phloem stream to the root) | 0.06 (range 0.02-0.10) | measured (non-cereal) | Windt et al. 2006 daytime phloem/xylem volume flow; no cereal value exists |
| T_C,Ph (phloem sap per kg grain DM) | 6.1 L/kg (range 4.4-7.8) | derived | rice sieve-tube sucrose (Chino 1988; Hayashi & Chino 1990) |
| panicle share of transpiration | 0 before heading; 0.21 (0.09-0.33) at DVS 1; 0.07 (0-0.14) at DVS 2; linear in DVS between | measured anchors; the linear form between them is stated as an assumption | Sugimoto 1973, 1975; Roberts & Batchelor 1983; Takahashi & Kurata 2007 |
| sheath + culm share | 0.29 (range 0.22-0.36) | measured (post-heading, extended to the season) | Sugimoto 1973, 1975 |
| retention | 1 | derived | solute delivered to an evaporating organ stays (a dead end); export is the phloem term |
| stem segments N | 4 | numerical | converged |
| organ-mass floor | 1e-9 kg; the grain formation gate is a step at the floor | numerical | converged (`numerics_convergence.py`) |
| metabolism gamma | 0 | structural | PFAAs are terminal, persistent products (Gredelj 2020; Just 2022; Arp 2024) |
| air exchange | 0 | structural | the anions do not volatilise (Vierke 2013; Johansson 2017) |
| root plaque pool K_surf | 0 (off) | structural | not identifiable: no dataset removed plaque |

## 4. Fits -- `validation/refit_sourced.py` (PREREG 5.3, 5.6)

Variants, each with its own fitted values: `redistributed_fit` (the default since amendment 2a:
f_xy fitted per congener to the straw as well, so on Yamazaki a reproduction; chosen by the
transfer rule, PREREG section 8), `redistributed_lipid` (amendment 1's default), `4pool_lipid`
(the single-straw core has no `_fit` variant; unspecified, it runs `4pool_lipid`), and
`redistributed_free` / `4pool_free`, which are fitted and runnable but withdrawn from scoring and
decisions (`model_api.WITHDRAWN_VARIANTS`; PREREG 9.4: their f_xy is untraceable). Results: `validation/refit_sourced.json`; parameters:
`params/refit_sourced.csv`, `params/lipid_loading_sourced.json` -> `parameters.json`.

## 5. Transfer datasets

| dataset | drivers used | why |
|---|---|---|
| Tang 2026 (screenhouse, Nipponbare, Hangzhou soil) | calibration site-year | the paper gives timing relative to exposure only (heading ~month 4, harvest day 150), no year or dates |
| Kim 2019 (30 Korean fields) | calibration site-year | harvest window Oct 6-18 2016 only; fields anonymised |
| Li 2025 (Tianjin, Jinyu Japonica No.22) | calibration site-year, stages by ORYZA DVS | samples July-October 2023; no stage dates. Stage -> DVS from sources (PREREG 5.9): heading 1 and ripening 2 by definition, tillering the midpoint of [transplanting, panicle initiation 0.65], filling and milk-ripening the centre (4/3, 5/3) of their ordered range; the range ends are printed |

Decision rule for every comparison: `validation/decision_stats.decide` (paired bootstrap, 95 %
interval; PREREG section 2). Scores: `validation/transfer_sourced.json` and
`params/transfer_scores.json` (the honesty band).

## 6. Soil side (HYDRUS-1D, flooded profile) -- `src/soil_hydrus.py` (not the default exposure)

| input | value | status |
|---|---|---|
| f_oc | 0.02 | **open**: Yamazaki SI Table S5 gives 3.7-3.8 % TOC (next amendment) |
| bulk density | 1.3 g/cm3 | **open**: 0.61 measured at a Tsukuba Andosol paddy (Koyama & Hayashi 2017) |
| percolation | 0.30 cm/d | **open**: 0.41-0.67 cm/d at an Ibaraki Andosol paddy (Kondo et al. 2019) |
| flooding | continuous to day 90 | **open**: midsummer drainage 35-40 d after transplanting for 15-20 d (Ibaraki calendar) |
| profile, dispersivity, node spacing | 100 cm, 5 cm, 2 cm | **open**: literature range or convergence |

## 7. App

| input | value | status | source |
|---|---|---|---|
| rice intake | 113.9 g/day white rice (range 101.68-126.13) | measured | KNHANES 2023 adult age-group means |
| body weight | 64.5 kg | measured | KNHANES 2013-2016 adults (Yoon et al. 2022) |
| rice eaten | white rice at 14 % moisture, from the dry brown-rice value and the white-rice composition | derived | `model_api.white_rice_conc` |
| intake signal | green below / red at or above the TWI | structural | EFSA 2020 TWI (the only sourced threshold) |
| tissue density (per-volume view) | removed (PREREG 9.5) | -- | it entered no result and had no source |
| SMILES K_cw | 20 L/kg placeholder | **open** | no source; decided by the K_cw side project (`docs/kcw_side_project/`, PREREG 9.5) |
