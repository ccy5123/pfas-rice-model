# Assumption register

Every input the shipped model runs on, with its source, or — where there is none — marked
**unmeasured** with the range it is run over and what that range does. Kept in step with the
code; the numbers in the *effect* column come from the scripts named there.

Conventions: *measured* = a value read from a cited table; *derived* = computed from measured
values by a stated rule; *convention* = a modelling choice with a stated reason;
*unmeasured* = no value exists for this crop/site; the default is stated, and the range is run.

**Effect** (sections 1–3): the median over the 11 Yamazaki congeners of |Δ log10 BAF|
(root / straw / grain, dry basis) when the input is moved to the end of its range with the
fitted parameters held fixed. Two configurations are reported, because the transport
parameters were fitted on one and the model now runs on the other:
- **default** — the redistributed shoot + lipid loading (`validation/default_sensitivity.py`;
  yardstick = the two observed harvest runs' model difference, **0.004 / 0.036 / 0.037**);
- **4-pool** — the single-straw free-anion fits (Fit A held; `validation/refit_site_year.py`,
  "Unknown inputs"; yardstick **0.003 / 0.027 / 0.067**). Inputs consequential there were
  re-fitted at both ends; "re-fit" gives the median (range) change of the fitted parameters.
An effect above its yardstick is **consequential** (bold).

## 1. Calibration site-year (Yamazaki et al. 2023, NARO Tsukuba 2020) — `src/site_scenario.py`

| input | value | status | source / range |
|---|---|---|---|
| weather | JMA Tsukuba 2020 daily + hourly | measured | `raw_si/jma_tsukuba_2020.csv` (radiation = sum of hourly, vapour pressure = mean of hourly) |
| transplanting | 2020-05-11 | measured | Yamazaki SI p.S2 |
| harvests | 2020-09-01 and 2020-09-10 (two cultivars) | measured | Yamazaki SI p.S3 |
| heading | harvest − 30.5 d | derived | the paper states 30–31 d; 30 / 31 d: default 0.000 / 0.004 / 0.003; 4-pool 0.000 / 0.003 / 0.009 |
| harvest = physiological maturity (DVS 2) | — | convention | IR72's own grain-fill rate kept instead: default 0.000 / 0.001 / 0.001; 4-pool 0.000 / 0.001 / 0.002 |
| phenology calibration | one factor on DVRJ/DVRI/DVRP → heading; DVRR → harvest | convention | ORYZA2000 manual §8.2 (DRATES); the juvenile/photoperiod split is not observable from two dates |
| crop parameters | IR72.D92 verbatim | convention | the cultivars (Kasalath, Koshihikari) have no ORYZA crop file; IR72.D92 fits IR72DS92 and is closer than the shipped driver on the Wuxi partitioning check |
| two cultivars → two runs, geometric mean | — | convention | the observations are the cultivar geometric mean; the single mean-date run (the app default): default 0.000 / 0.004 / 0.005; 4-pool 0.000 / 0.003 / 0.008 |
| nursery (SBDUR) | 23 d after emergence, open seedbed | unmeasured | Kim 2003 nurseries 14/23/24 d; 14 / 30 d: default 0.001 / 0.008 / 0.004, 4-pool 0.001 / 0.008 / 0.025; covered seedbed (TMPSB 9 °C): default 0.002 / 0.014 / 0.005, 4-pool 0.002 / 0.014 / 0.048 |
| planting density | 19.05 hills m⁻², 3 plants hill⁻¹ | unmeasured | Japanese machine spacing 30 × 15–22 cm; 15 / 22 hills m⁻²: default 0.001 / 0.006 / 0.005, 4-pool 0.001 / 0.006 / 0.024; 5 plants hill⁻¹: default 0.002 / 0.013 / 0.004, 4-pool 0.001 / 0.014 / 0.046 |
| CO₂ | 414 ppm | measured (global) | Mauna Loa 2020 annual mean; 400 / 430: default 0.000 / 0.007 / 0.008; 4-pool 0.000 / 0.006 / 0.008 |
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
| dead blades in straw / leaf | all kept (attached = 1) | unmeasured | ORYZA's WAGT convention and a whole-plant harvest; all lost: default **straw 0.045** (root, grain 0); 4-pool **straw 0.174**, re-fit f_xy ×1.31 (1.00–1.53), L_Ph ×1.00 (0.44–1.00) |
| dead tissue exchange | none after death | convention | no wash-off, no remobilisation |

## 3. Plant model constants

| constant | value | status | range run / note |
|---|---|---|---|
| θ_fw (root/stem/leaf/grain) | 0.90 / 0.83 / 0.78 / 0.14 | estimates, **source not traceable** | `params/rice_tissue_params.csv` ranges (its `source_key` column -- `S_physiol`, `S_grainMC`, ... -- resolves to no reference anywhere in the repo, so the rows labelled "measured" there are not checkable either); `validation/theta_fw_sensitivity.py`: root INTERMEDIATE, stem/leaf/grain SMALL (pre-registered rule; dry-basis metrics) |
| composition f_prot, f_PL, f_cw | `rice_tissue_params.csv` | **source not traceable** | the file labels each row measured or estimate, but none of its source keys resolves to a reference in the repo; no sensitivity has been run on f_prot, f_PL or f_cw |
| K_PL, K_prot | Chen 2025, Zhou 2025 | measured | |
| K_cw | GAP A deliverable | derived | no measured coefficient exists |
| E_m | −120 mV | measured (Wang 1994) | −90 mV (NH₄⁺ depolarisation): default **0.128 / 0.126 / 0.127**; 4-pool **0.128 / 0.126 / 0.127**, re-fit κ_d ×1.46 (1.33–2.29), f_xy and L_Ph unchanged |
| carrier Vmax_in / Km / Vmax_out | 20 / 5 µg L⁻¹ / 8 | **unmeasured, fixed, no source** (the values are chosen, not taken from any study) | Km × 100 (the dose-series lower bound ≥ 500 µg/L): default **1.27 / 1.27 / 1.27**; 4-pool **1.26 / 1.26 / 1.26**, re-fit: κ_d → ~0 (the now-linear carrier takes over), f_xy ×3.9 (1.0–28), L_Ph ×1.1 (1.0–9000). The level of every fitted transport parameter depends on this unmeasured constant. |
| φ (phloem recirculation to the root) | 0.1 | unmeasured, no source | 0.05 / 0.2: default 0.000 / 0.026 / 0.002 (not consequential); 4-pool 0.001 / **0.048** / 0.026 (straw, through the root's recirculated burden), re-fit f_xy ×0.94 / ×1.11 (up to ×5.7) |
| T_C,Ph (phloem flow per kg grain DM) | 10 L/kg | literature order of magnitude | 4 / 12 L/kg (sap 100–300 g/L sugar + grain water): default 0.000 / 0.009 / 0.017 (not consequential: PFOA's default grain is ~80 % xylem-fed, `apportionment`); 4-pool 0.000 / 0.012 / **0.160** — it sets the whole grain phloem volume there; re-fit L_Ph to its bound or ×0.85 |
| grain xylem surface (4-pool only) | 0.5116 m² per kg DM (~6 % of the leaf + grain stream from panicle formation to harvest on the site-year, transpiration-weighted; 8 % from heading) | unmeasured | 0.2 (~2.5 %) / 2.0 (~18 %): 4-pool 0.000 / **0.044 / 0.280**; re-fit L_Ph ×1.31 (up to ×300) / ×0.15 (down to the floor). The default shoot splits the stream by `lam_grain` instead (below) |
| redistributed shoot: `lam_grain` (grain share of the stream) | 0.05 | unmeasured | 0.025 / 0.18 (the share the 4-pool range spans): default 0.000 / 0.034 / **0.063** |
| redistributed shoot: grain routing of the residual xylem stream (`grain_share`) | "gate": all of it once the panicle forms | structural choice, tested | "size" (share = panicle mass / mature mass): Yamazaki 0.380 → 0.254, Kim interval 0.154 → 0.267, Li grain TF 0.510 → 0.940; kept "gate" by the pre-registered rule (`validation/panicle_share_amendment.py`) |
| redistributed shoot: `retention` (share of the deposited solute an organ keeps) | 0.6 | unmeasured, no source | ± one third, 0.40 / 0.80: default 0.000 / **0.110 / 0.225** |
| redistributed shoot: `stem_transp_frac` (stem share of the shoot stream) | 0.45 | unmeasured, no source | ± one third, 0.30 / 0.60: default 0.000 / 0.011 / 0.008 |
| redistributed shoot: stem segments N | 4 | convention | the configuration scored in every transfer test |
| leaf xylem surface (4-pool only) | 20 m² per kg DM | convention, **value not sourced** | the original S read as an SLA; not a rice SLA measurement -- IR72.D92's own SLA is 24–45 m² per kg (ASLA 0.0024 + BSLA 0.0025·e^(−4.5(DVS−0.14)) ha/kg, capped at SLAMAX 0.0045); not run over a range |
| grain formation gate | opens as the grain mass leaves the 1e-4 kg organ-mass floor, fully open at 1.5 × the floor (one day later on the site-year) | **numerical choice, no source** | `validation/gate_numerics_sensitivity.py` (ramp 1.1 / 3 / 10 ×, floor 1e-5 / 1e-3 kg): the PFOA young-panicle peak moves 1.41–9.65× its harvest value (shipped 3.29×) -- its size is set by these two numbers, not by data; the scores barely move: Yamazaki 0.367–0.384 (shipped 0.380), grain BAFs median \|Δ log10\| ≤ 0.024, transfer endpoints ≤ 0.038 (Li grain 0.498–0.548) |
| root surface / plaque pool K_surf | 0 | convention | no dataset removed plaque; absorbed by κ_d in the fit |
| f_xy bounds | ≤ 1 | physical | |
| f_xy head-group offset (free anion, `f_xy_recommended`) | PFSA = PFCA × e^−1.1, ether × e^−0.7 | derived from TF ratios | PFOS/PFOA TF of Tang 2026 (0.26) and Yamazaki 2023 (0.43) -- so Tang's PFOS is not independent of the free-anion f_xy; the default (lipid loading) uses the fitted `pfsa_ln` instead. GenX from Tang alone (provisional) |
| organ-mass floor | 1e-4 kg (`organ_mass.organ_matrix`) | numerical choice | keeps the ODE off zero mass; also sets when the grain gate opens (row above) |
| L_Ph, κ_d bounds | [1e-5, 1], [1e-2, 50] | convention | L_Ph = 1 means sieve-tube sap at the leaf's free concentration (no active concentrating); the upper bound BINDS for PFDA, PFUnDA, PFOS and PFDoDA, whose grain stays under -- whether phloem loading can concentrate PFAS is unmeasured |

## 4. Transport fit — `validation/refit_site_year.py` (pre-registered, amendment 1)

| item | value | status |
|---|---|---|
| `*_fit` | (f_xy, L_Ph, κ_d) per congener, saturated | fitted (reproduction) |
| `*_rec` | (L_Ph, κ_d) with the monotone f_xy fixed | fitted; straw a prediction |
| PFBS L_Ph | geometric mean of the others | convention (no grain value) |
| GenX / SMILES records | geometric mean of the fitted congeners' L_Ph, κ_d | convention |
| lipid globals | six, on the `*_rec` L_Ph/κ_d, PFDoDA excluded | fitted |

Fit A reproduces the 32 used cells at log10 RMSE 0.162 (root 0.000, straw 0.035, grain
0.287; saturated -- a reproduction). PFBA and PFPeA sit at the L_Ph floor with their grain
over (xylem-fed); PFDA, PFUnDA and PFOS reach the L_Ph upper bound (1) with their grain
under, and PFDoDA sits at f_xy = L_Ph = 1 with straw and grain under: with the grain drawing
only its growth sap, the free anion cannot carry the mid/long-chain grain. Fit B's straw (a
prediction) is at log10 RMSE 1.401. The lipid globals reach 0.207 over 29 cells.

**The default configuration** (redistributed shoot + lipid loading; `validation/default_switch.py`)
adds no fitted parameter: it runs the `*_rec` L_Ph/κ_d and the lipid globals, all fitted on the
4-pool. On the Yamazaki observable it scores 0.380 (root 0.090, straw 0.253, grain 0.619) and
over-predicts the grain of most congeners (PFPeA 7.1×, PFUnDA 7.3×, PFOA 3.2×).

**Reading the consequential inputs together**: the fitted (f_xy, L_Ph, κ_d) are conditional
on the carrier constants (Km above all), E_m, φ, T_C,Ph and the grain xylem surface, none of which is
measured for rice; the default's results are, in addition, conditional on the shoot levers
`retention` and `lam_grain`, which are not measured either. A fit that reproduces Yamazaki does not pin them; a transfer test does
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
| Kd | Koc(chain length) × f_oc, f_oc 0.02 | QSPR (measured anchors); f_oc unmeasured for the site, no range run |
| flooding | continuous to day 90, then drainage | convention (Yamazaki: continuous except midsummer drainage; dates not reported) |
| percolation | 0.30 cm/d | chosen, no source, no range run |
| profile | 100 cm, ρ_b 1.3, dispersivity 5 cm | chosen, no source, no range run |

## 7. External data that would replace an assumption

- Yamazaki 2023 SI (S16 water, S18/S19 raw values and MQL flags, planting density, nursery).
- A panicle / caryopsis transpiration share (grain xylem surface).
- Phloem sap volume per grain dry matter for rice (T_C,Ph).
- Organ water content over a season (θ_fw).
- A rice-root cell-wall / plaque PFAS sorption coefficient (K_cw, K_surf).

## 8. Decision thresholds

The pre-registered analyses decide by margins. The margins were fixed before each run (so a
result could not tune them), but no study supplies their values: they are choices. For each
decision, the margin at which the verdict would change, from the recorded scores:

| decision (script) | margin used | recorded result | verdict holds for |
|---|---|---|---|
| default switch D1, beat the old default on every transfer endpoint (`default_switch.py`) | > 0.3 log | smallest gain 0.90 (Li shoot 1.26 → 0.36) | any margin < 0.90 |
| default switch D1, within 0.1 of the 4-pool lipid variant | 0.1 log | better on all five endpoints | any margin ≥ 0 |
| panicle amendment A2, transfer non-inferiority of `"size"` (`panicle_share_amendment.py`) | 0.10 log | Kim interval +0.113, Li grain +0.430 | FAIL for any margin < 0.43 (from 0.12 only Li grain fails) |
| Kim robust winner between the two lipid variants (`kim2019_censored.py`) | 0.10 log on point and interval | point 0.56 vs 0.47 (0.09), interval 0.38 vs 0.15 | **close**: a margin below ~0.09 names the default |
| θ_fw: constant θ kept if every metric < 0.05, θ(t) built if any > 0.10 (`theta_fw_sensitivity.py`) | 0.05 / 0.10 log | root 0.053 (Ge) → INTERMEDIATE | **close**: a lower cut of 0.054 would make root SMALL |
| panicle amendment A1, `"size"` has no young-panicle peak | max ≤ 1.5 × harvest value | `"size"` max 1.00× (`"gate"` 3.2–3.8×) | any cut ≥ 1.00; the `"gate"` peak size itself depends on the gate numerics (§3, grain formation gate) |

