# Assumption register

Every input the shipped model runs on, with its source, or — where there is none — marked
**unmeasured** with the range it is run over and what that range does. Kept in step with the
code; the numbers in the *effect* column come from the scripts named there.

Conventions: *measured* = a value read from a cited table; *derived* = computed from measured
values by a stated rule; *convention* = a modelling choice with a stated reason;
*unmeasured* = no value exists for this crop/site; the default is stated, and the range is run.

**Effect** columns (sections 1–3): the median over the 11 Yamazaki congeners of
|Δ log10 BAF| (root / straw / grain, dry basis) when the input is moved to the end of its
range with the fitted transport parameters held fixed (`validation/refit_site_year.py`,
"Unknown inputs"). The yardstick is the model difference between the two observed harvest
runs (two cultivars, nine days apart): **0.003 / 0.033 / 0.055**. An input that moves a
tissue by more than that is marked **consequential** (bold) and was re-fitted at both ends;
the re-fit column gives the median (range) change of the fitted parameters.

## 1. Calibration site-year (Yamazaki et al. 2023, NARO Tsukuba 2020) — `src/site_scenario.py`

| input | value | status | source / range |
|---|---|---|---|
| weather | JMA Tsukuba 2020 daily + hourly | measured | `raw_si/jma_tsukuba_2020.csv` (radiation = sum of hourly, vapour pressure = mean of hourly) |
| transplanting | 2020-05-11 | measured | Yamazaki SI p.S2 |
| harvests | 2020-09-01 and 2020-09-10 (two cultivars) | measured | Yamazaki SI p.S3 |
| heading | harvest − 30.5 d | derived | the paper states 30–31 d; 30 / 31 d: effect 0.000 / 0.004 / 0.008 |
| harvest = physiological maturity (DVS 2) | — | convention | IR72's own grain-fill rate kept instead: 0.000 / 0.001 / 0.002 |
| phenology calibration | one factor on DVRJ/DVRI/DVRP → heading; DVRR → harvest | convention | ORYZA2000 manual §8.2 (DRATES); the juvenile/photoperiod split is not observable from two dates |
| crop parameters | IR72.D92 verbatim | convention | the cultivars (Kasalath, Koshihikari) have no ORYZA crop file; IR72.D92 fits IR72DS92 and is closer than the shipped driver on the Wuxi partitioning check |
| two cultivars → two runs, geometric mean | — | convention | the observations are the cultivar geometric mean; the single mean-date run (the app default): 0.000 / 0.003 / 0.010 |
| nursery (SBDUR) | 23 d after emergence, open seedbed | unmeasured | Kim 2003 nurseries 14/23/24 d; 14 / 30 d: 0.001 / 0.010 / 0.027; covered seedbed (TMPSB 9 °C): 0.002 / 0.017 / 0.051 |
| planting density | 19.05 hills m⁻², 3 plants hill⁻¹ | unmeasured | Japanese machine spacing 30 × 15–22 cm; 15 / 22 hills m⁻²: 0.001 / 0.008 / 0.026; 5 plants hill⁻¹: 0.002 / 0.016 / 0.048 |
| CO₂ | 414 ppm | measured (global) | Mauna Loa 2020 annual mean; 400 / 430: 0.000 / 0.007 / 0.012 |
| reference ET terms | FAO-56 Penman–Monteith radiation / aerodynamic terms | convention | the ORYZA2000 manual does not print its Penman routine; FAO-56 is the definition it cites |
| field reflection | ORYZA eq. 4.2 with RFS 0.05 under standing water | measured (manual) | manual §4.1 |
| transpiration | ORYZA eq. 4.4, actual = potential | convention | continuous flooding except a midsummer drainage; no water stress modelled |
| wind height | 10.6 m → 2 m (FAO-56 eq. 47) | measured | JMA station metadata |

## 2. Observation (Yamazaki 2023) — `data_obs/obs_baf_Yamazaki.csv`, `src/harvest.py`

| item | value | status | note |
|---|---|---|---|
| tissue basis | dry weight (pg/g dw) | measured | the fresh table it replaces was these × assumed water contents; inverted exactly |
| exposure | one irrigation-water grab, 2020-06-08, held constant | measured value, convention for time | the seasonal course of the water is unknown |
| straw | unweighted mean of four height bands | measured as reported | band masses not reported; the band range bounds any mass-weighted value (`straw_lo/hi`) |
| brown rice | 0.88 white + 0.12 bran | derived | white:brown 86–90 % by mass (SI p.S3) |
| PFBS grain | not used | — | wholly imputed (< MQL) |
| other < MQL cells | MQL/2 upstream | unknown which | needs SI S18/S19 |
| σ | 0.4 log10 | unmeasured | n = 1 composite per tissue |
| dead blades in straw / leaf | all kept (attached = 1) | unmeasured | ORYZA's WAGT convention and a whole-plant harvest; all lost: **straw 0.174** (root, grain 0); re-fit f_xy ×1.41 (1.00–1.53), L_Ph ×0.73 (0.45–1.00) |
| dead tissue exchange | none after death | convention | no wash-off, no remobilisation |

## 3. Plant model constants

| constant | value | status | range run / note |
|---|---|---|---|
| θ_fw (root/stem/leaf/grain) | 0.90 / 0.83 / 0.78 / 0.14 | estimates | `params/rice_tissue_params.csv` ranges; `validation/theta_fw_sensitivity.py`: root INTERMEDIATE, stem/leaf/grain SMALL (pre-registered rule; dry-basis metrics) |
| composition f_prot, f_PL, f_cw | `rice_tissue_params.csv` | measured / estimate per row | |
| K_PL, K_prot | Chen 2025, Zhou 2025 | measured | |
| K_cw | GAP A deliverable | derived | no measured coefficient exists |
| E_m | −120 mV | measured (Wang 1994) | −90 mV (NH₄⁺ depolarisation): **0.128 / 0.126 / 0.127**; re-fit κ_d ×1.46 (1.33–2.29), f_xy and L_Ph unchanged |
| carrier Vmax_in / Km / Vmax_out | 20 / 5 µg L⁻¹ / 8 | **unmeasured, fixed** | Km × 100 (the dose-series lower bound ≥ 500 µg/L): **1.26 / 1.26 / 1.26**; re-fit: κ_d → ~0 (the now-linear carrier takes over), f_xy ×3.9 (1.0–30), L_Ph ×1.8 (1.0–3300). The level of every fitted transport parameter depends on this unmeasured constant. |
| φ (phloem recirculation) | 0.1 | unmeasured | 0.05 / 0.2: 0.000 / 0.028 / **0.129**; re-fit L_Ph ×1.5 / ×0.6 |
| T_C,Ph (phloem flow per kg grain DM) | 10 L/kg | literature order of magnitude | 4 / 12 L/kg (sap 100–300 g/L sugar + grain water): 0.000 / 0.006 / 0.053 |
| grain xylem surface | 0.5116 m² per kg DM (~6 % of the stream during filling on the site-year) | unmeasured | 0.2 (~2.5 %) / 2.0 (~18 %): 0.000 / **0.041 / 0.235**; re-fit L_Ph ×1.37 (up to ×104) / ×0.16 (down to the floor) |
| leaf xylem surface | 20 m² per kg DM | convention (SLA reading) | |
| grain formation gate | mass leaves its floor | convention | |
| root surface / plaque pool K_surf | 0 | convention | no dataset removed plaque; absorbed by κ_d in the fit |
| f_xy bounds | ≤ 1 | physical | |
| L_Ph, κ_d bounds | [1e-5, 1], [1e-2, 50] | inherited | as every earlier fit |

## 4. Transport fit — `validation/refit_site_year.py` (pre-registered, amendment 1)

| item | value | status |
|---|---|---|
| `*_fit` | (f_xy, L_Ph, κ_d) per congener, saturated | fitted (reproduction) |
| `*_rec` | (L_Ph, κ_d) with the monotone f_xy fixed | fitted; straw a prediction |
| PFBS L_Ph | geometric mean of the others | convention (no grain value) |
| GenX / SMILES records | geometric mean of the fitted congeners' L_Ph, κ_d | convention |
| lipid globals | six, on the `*_rec` L_Ph/κ_d, PFDoDA excluded | fitted |

Fit A reproduces the 32 used cells at log10 RMSE 0.121 (root 0.000, straw 0.034, grain
0.213; saturated -- a reproduction). PFBA and PFPeA sit at the L_Ph floor with their grain
over; PFDoDA at f_xy = 1 with its straw and grain under. Fit B's straw (a prediction) is at
log10 RMSE 1.372. The lipid globals reach 0.214 over 29 cells (g_ph at its lower bound).

**Reading the consequential inputs together**: the fitted (f_xy, L_Ph, κ_d) are conditional
on the carrier constants (Km above all), E_m, φ and the grain xylem surface, none of which is
measured for rice. A fit that reproduces Yamazaki does not pin them; a transfer test does
not either, because each transfer dataset also runs on borrowed drivers (section 5).

## 5. Transfer datasets (out-of-sample)

| dataset | drivers used | why |
|---|---|---|
| Tang 2026 (pot, Nipponbare) | calibration site-year | Tang's site and dates are not in the repo |
| Kim 2019 (30 Korean fields) | calibration site-year | no dated crop record per field |
| Li 2025 (one field, Jinyu Japonica No.22) | calibration site-year, stages by ORYZA DVS | stage dates not reported |
| neutral tables (Liu, Ge, Kodešová, Li 2019, …) | generic (`forcing_rice` + `growth_rice`, `neutral_dpu_validation.drivers`) | other species/sites; the neutral validation is about the partition core |
| Schriever TSCF weak-electrolyte rows | calibration site-year (`model_api.simulate_neutral` default) | the comparison is speciation ON vs OFF on the same drivers |

## 6. Soil side (HYDRUS-1D, flooded profile) — `src/soil_hydrus.py`

| input | value | status |
|---|---|---|
| Kd | Koc(chain length) × f_oc, f_oc 0.02 | QSPR (measured anchors), f_oc unmeasured for the site |
| flooding | continuous to day 90, then drainage | convention (Yamazaki: continuous except midsummer drainage; dates not reported) |
| percolation | 0.30 cm/d | convention |
| profile | 100 cm, ρ_b 1.3, dispersivity 5 cm | convention |

## 7. External data that would replace an assumption

- Yamazaki 2023 SI (S16 water, S18/S19 raw values and MQL flags, planting density, nursery).
- A panicle / caryopsis transpiration share (grain xylem surface).
- Phloem sap volume per grain dry matter for rice (T_C,Ph).
- Organ water content over a season (θ_fw).
- A rice-root cell-wall / plaque PFAS sorption coefficient (K_cw, K_surf).
