# Assumption register

Every input the shipped model runs on, with its source, or — where there is none — marked
**unmeasured** with the range it is run over and what that range does. Kept in step with the
code; the numbers in the *effect* column come from the scripts named there.

Conventions: *measured* = a value read from a cited table; *derived* = computed from measured
values by a stated rule; *convention* = a modelling choice with a stated reason;
*unmeasured* = no value exists for this crop/site; the default is stated, and the range is run.

## 1. Calibration site-year (Yamazaki et al. 2023, NARO Tsukuba 2020) — `src/site_scenario.py`

| input | value | status | source / range |
|---|---|---|---|
| weather | JMA Tsukuba 2020 daily + hourly | measured | `raw_si/jma_tsukuba_2020.csv` (radiation = sum of hourly, vapour pressure = mean of hourly) |
| transplanting | 2020-05-11 | measured | Yamazaki SI p.S2 |
| harvests | 2020-09-01 and 2020-09-10 (two cultivars) | measured | Yamazaki SI p.S3 |
| heading | harvest − 30.5 d | derived | the paper states 30–31 d; sensitivity 30 / 31 d |
| harvest = physiological maturity (DVS 2) | — | convention | sensitivity: IR72's own grain-fill rate kept (harvest not at DVS 2) |
| phenology calibration | one factor on DVRJ/DVRI/DVRP → heading; DVRR → harvest | convention | ORYZA2000 manual §8.2 (DRATES); the juvenile/photoperiod split is not observable from two dates |
| crop parameters | IR72.D92 verbatim | convention | the cultivars (Kasalath, Koshihikari) have no ORYZA crop file; IR72.D92 fits IR72DS92 and is closer than the shipped driver on the Wuxi partitioning check |
| two cultivars → two runs, geometric mean | — | convention | the observations are the cultivar geometric mean; the single mean-date run is the app default and is run as a sensitivity |
| nursery (SBDUR) | 23 d after emergence, open seedbed | unmeasured | Kim 2003 nurseries 14/23/24 d; sensitivity 14 / 30 d and a covered seedbed (TMPSB 9 °C) |
| planting density | 19.05 hills m⁻², 3 plants hill⁻¹ | unmeasured | Japanese machine spacing 30 × 15–22 cm; sensitivity 15 / 22 hills m⁻², 5 plants hill⁻¹ |
| CO₂ | 414 ppm | measured (global) | Mauna Loa 2020 annual mean; sensitivity 400 / 430 |
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
| dead blades in straw / leaf | all kept (attached = 1) | unmeasured | ORYZA's WAGT convention and a whole-plant harvest; sensitivity 0 (all lost) |
| dead tissue exchange | none after death | convention | no wash-off, no remobilisation |

## 3. Plant model constants

| constant | value | status | range run / note |
|---|---|---|---|
| θ_fw (root/stem/leaf/grain) | 0.90 / 0.83 / 0.78 / 0.14 | estimates | `params/rice_tissue_params.csv` ranges; `validation/theta_fw_sensitivity.py` |
| composition f_prot, f_PL, f_cw | `rice_tissue_params.csv` | measured / estimate per row | |
| K_PL, K_prot | Chen 2025, Zhou 2025 | measured | |
| K_cw | GAP A deliverable | derived | no measured coefficient exists |
| E_m | −120 mV | measured (Wang 1994) | −90 mV (NH₄⁺ depolarisation) |
| carrier Vmax_in / Km / Vmax_out | 20 / 5 µg L⁻¹ / 8 | **unmeasured, fixed** | Km × 100 (the dose-series lower bound ≥ 500 µg/L) |
| φ (phloem recirculation) | 0.1 | unmeasured | 0.05 / 0.2 |
| T_C,Ph (phloem flow per kg grain DM) | 10 L/kg | literature order of magnitude | 4 / 12 L/kg (sap 100–300 g/L sugar + grain water) |
| grain xylem surface | 0.5116 m² per kg DM (~4 % of the stream at mid-filling) | unmeasured | 0.2 (~2 %) / 2.0 (~14 %) |
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

The effect of each unmeasured input on the fitted model is in the refit output (section
"Unknown inputs"); see CLAUDE.md §6 for the numbers.

## 5. Transfer datasets (out-of-sample)

| dataset | drivers used | why |
|---|---|---|
| Tang 2026 (pot, Nipponbare) | calibration site-year | Tang's site and dates are not in the repo |
| Kim 2019 (30 Korean fields) | calibration site-year | no dated crop record per field |
| Li 2025 (one field, Jinyu Japonica No.22) | calibration site-year, stages by ORYZA DVS | stage dates not reported |
| neutral tables (Liu, Ge, Kodešová, …) | generic (`forcing_rice` + `growth_rice`) | other species/sites; the neutral validation is about the partition core |

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
