## 7. Amendment 2 -- what amendment 1 left open (DRAFT; commit before its runs)

Default scope (named in tests/test_value_provenance.py OPEN):
- crop_cultivar: no ORYZA set for Koshihikari/Kasalath exists. Run the published japonica sets
  (Gao et al. 2021 Heilongjiang, already transcribed; Wu et al. 2022 Liaoning ORYZA v3
  partitioning + SLA) beside IR72 as a crop-parameter range on the site-year (phenology
  recalibrated to the measured dates in each); carry the spread in the input-range band.
- obs_exposure_time: the one irrigation-water grab (2020-06-08) held over the season. Options:
  the measured rainwater series (Table S16 R1-R5, Table S17 fluxes) as a seasonal shape? No --
  rain is not the irrigation water. The honest form is a range over the shapes the data allow:
  constant (default) vs the HYDRUS paddy shape driven by the measured soil before/after
  (Table S16 soil columns). To decide.

Other scopes:
- soil (HYDRUS / flooded profile): f_oc from Yamazaki SI Table S5 (3.7-3.8 % TOC, 0-15 cm);
  bulk density 0.61 g/cm3 (Koyama & Hayashi 2017, Tsukuba Andosol); percolation 0.41-0.67 cm/d
  (Kondo et al. 2019, Ibaraki Andosol paddy); midsummer drainage 35-40 DAT for 15-20 d (Ibaraki
  calendar); dispersivity: literature range; profile depth and node spacing: convergence;
  _KLEACH_MAX and the 0.02 fallback: remove or derive.
- neutral path: RICE_WATER -> the sourced composition theta; RICE_SURFACE -> lam_grain(DVS);
  leaf sorption exponent 0.95 and the K_AW volatility cut 1e-4: source or remove; pore-water pH
  6.5 (weak electrolytes): source a paddy pore-water pH range.
- transfer analyses: Tang/Kim/Li on their own drivers where the papers give site and dates; the
  Li stage -> DVS mapping through a published staging system (Counce et al. 2000) or a range.
- app: exposure defaults (rice intake, body weight) from the Korean national survey (KNHANES)
  with the year stated; the per-volume tissue density: source or remove the per-volume view;
  SMILES K_cw placeholder: derive by chain length (GAP A method) or refuse.
- free variants: record the points of the f_xy_recommended OLS fit or refit from the cited TSCF
  tables.
- calibration.py synthetic sigma 0.3; bayesian_inverse_demo SIGMA 0.05: state as synthetic-noise
  demo settings or derive.

Transfer calendars (searched 2026-10-25; full texts + SI read):
- Li 2025: Tianjin (38 34 N, 116 43 E), samples July-October 2023, five stages, no dates, no
  staging system cited ("filling" listed before "milk-ripening").
- Kim 2019: harvest period Oct 6-18 2016, 30 anonymised Korean fields, japonica, no cultivar.
- Tang 2026 (JHM 502:141017): screenhouse, "time of rice cultivation was consistent with natural
  sowing", Hangzhou soil; months 2-3 tillering/jointing, months 4-5 heading to maturity, harvest
  at exposure day 150; year not given.
=> none gives sowing/heading dates; own-driver runs would need assumed calendars. Options: keep
   the borrowed site-year drivers (stated) or run climatological calendars as a range.

Exposure (app) -- found 2026-10-25:
- rice intake: KNHANES 2023 (KHIDI 다소비식품 섭취량 통계, data.go.kr 15160399) 멥쌀, adults
  19-29 / 30-49 / 50-64 y: 101.68 / 106.17 / 126.13 g/day (both sexes); all ages 117.30.
  (Statistics Korea 양곡소비량조사 2025: 147.7 g/day household raw grain, excl. eating out.)
- body weight: KNHANES 2013-2016 adults 64.5 kg (Yoon et al. 2022 JESEE 33:1013).
- BASIS BUG: intake_fraction is fed the model's BROWN-rice concentration at HARVEST moisture
  (theta 0.21), while the intake is milled white rice (멥쌀). Fix: white-rice concentration from
  the dry brown value x endosperm_white_factor (sourced white composition) at the white rice's own
  moisture (storage standard 0.14, sourced tis_theta_grain).
Paddy pore-water pH: 6.5-6.7 within 3 weeks of submergence (Yoshida 1981 p.173 citing
  Ponnamperuma 1965); 6.5-7.0 (Moraes 1973) -> replace "6.5 (range 5-7)".

Neutral path consequence of amendment 1 (found 2026-10-25): organ_mass now converts dry -> fresh
with the sourced theta (root 0.713, one pot study at maturity; was 0.90 unsourced), while
neutral_dpu.RICE_WATER (root 0.90) still sets the water term of K_PW -> the two disagree inside the
neutral path. Liu 2023 a-priori score moved 0.213 -> 0.226 on that change alone. Amendment 2:
one theta per organ for both uses; search more rice ROOT water contents (only one source found).

Agent report 2026-10-25 (staging + crop sets):
- Li stage mapping: DONE in amendment 1 addendum 5.9 (sourced definitions + ordered-range centroid).
- Wu et al. 2022 AEER 20(1):633 (OA, aloki.hu/pdf/2001_633651.pdf): prints DVRJ/DVRI/DVRP/DVRR (5 sites),
  FST/FLV/FSO to DVS 1.2 and SLA (5 sites); NOT FSHTB, leaf death, RGRL, SPGF, WGRMX ("bulk of our model
  parameters were taken from IR72"). Kaiyuan partitioning = Donggang (possible copy error).
- No ORYZA set for Koshihikari/other Japanese japonica published (Akitakomachi AgMIP set on request only).
=> crop_cultivar range: IR72 (shipped) vs Gao 2021 vs Wu 2022 partitioning+SLA (IR72 for the rest), phenology
   recalibrated to the measured dates in each.

Agent report 2026-10-25 (water/lipid/density) -- recorded in sourced_inputs.csv:
- ROOT theta: maturity pots 0.506/0.527/0.590/0.622/0.767 (+0.713), vegetative 0.761/0.824/0.845-0.851/0.867-0.888.
  -> theta falls over the season (all organs: grain 0.58 early filling -> 0.20 harvest; leaf 0.731 tillering ->
  0.659 flag leaf maturity). With fast linear root entry, Cw_root = r_eq*Cwo whatever theta; root theta enters via
  the fitted r_eq (root dw obs at HARVEST -> maturity theta) and via the lipid term g*C_root (season theta).
  DECIDE: theta(t) per organ from stage anchors (structural; state = burden per dry mass is cleanest) vs constant
  at the all-stage midpoint with the range in the band. Needs a refit.
- Neutral lipid dry basis: brown rice 1.9-4.0 % DW; straw 0.9-2.1 % DW (mean 1.4); leaf ~1-3 % DW (forage EE);
  ROOT: no total-lipid measurement (phospholipid floor 1.9-2.2 % DW only). -> restate L per dry mass; root range
  [phospholipid floor, ?] needs an upper source.
- Density: grain particle 1.13-1.31 (rough..brown, FAO), root ~0.65-1.00 (porosity-derived), culm/leaf none ->
  remove the per-volume view for stem/leaf or the whole view.
- Neutral validation: compare_to_obs converts dw observations of OTHER species / seedling experiments (Kodesova,
  Li 2019, Briggs barley) to fresh with the RICE model's water (pre-existing flaw, exposed by 5.8). Fix: compare on
  the basis each table reports (dw model output via the harvest operator) or use each paper's own water content.

## theta(t) design (prototype in scratchpad/a2; NOT run or fitted before the amendment is committed)
Physics: with burden P = C*M_fresh, Cw = P / (M_dry*(theta/(1-theta) + X)). The fresh formulation stays valid
with theta(t) iff M_fresh(t) and B(t) use the same theta(t); the grain phloem stream needs DRY growth taken
from the dry driver (not dM_fresh*(1-theta), wrong when theta moves). Harvest operator: dry = fresh*(1-theta(t_h)).
Anchors (DVS; linear in DVS between anchors -- the form is the stated assumption, as for the panicle share):
- root: vegetative (DVS <= 0.65): mid of [0.761, 0.888]; maturity (DVS 2): mid of [0.506, 0.767] (6 studies incl. 0.713)
- grain: DVS 1: 0.58 (Yoshida: initial phase of filling); DVS 2: mid of [0.20, 0.22]
- leaf: DVS <= 0.65: 0.731 (tillering, 183 accessions); DVS 2: 0.659 (flag leaf at maturity)
- stem: whole-shoot proxy -- vegetative mid of [0.726, 0.91] (soil pots 56 DAS; hydroponic seedlings); DVS 2: 0.729
  (straw at maturity); the 10-45 DAF whole-shoot 0.645-0.774 between (check, or a range input).
Band: anchors' range ends one at a time (input_range_band). Refit required (r_eq absorbs harvest root theta).

## DRAFT prereg text, theta(t) (commit before the amendment-2 refit)
- Organ water follows the development stage: theta(DVS) linear between sourced anchors, flat outside
  (params/rice_theta_stage_sourced.csv). Root 0.8245 (vegetative, range 0.761-0.888) -> 0.6365 at maturity
  (0.506-0.767); stem 0.818 (0.726-0.91) -> 0.7095 (0.645-0.774); leaf 0.731 -> 0.659; grain 0.58 at DVS 1 ->
  0.21 (0.20-0.22). Fresh mass and the binding factor read the same theta(t); the grain phloem stream uses the dry
  driver's growth; the harvest operator divides burden by dry matter at harvest.
- Implemented and tested in scratchpad/a2 (tests/test_theta_stage.py: mass conservation of both shoots, flux
  mirror, dry growth, harvest dry basis).
- Re-fit: all variants (protocol 5.6 unchanged). Band: each anchor's ends.
- Expectations: E-a: short-chain r_eq (PFBA, PFPeA, PFHxA, PFBS) rise > 20 % vs amendment 1, long chains
  (PFDA, PFUnDA, PFDoDA, PFOS) < 10 % (the harvest root water term 1.75 vs 2.48 matters only where binding is small).
  E-b: the root's maturity water anchor is among the three largest contributors to the root band for the short chains.
  E-c: the Yamazaki observable of the default moves by < 0.05 log10 (the fits absorb the root water).
- Neutral path: NOT in this part (no rice root lipid measurement; restating L per dry mass needs one) -> 2b.
- SMILES (app): compound_from_smiles defaults L_Ph 0.005, kappa_d 0.5, Vmax/Km 20/5/8/5 have no source -- check
  whether simulate_from_smiles uses them at all (transport for a novel structure is the fitted congeners' median);
  remove if dead. K_cw 20 placeholder: derive per organ from the GAP A ladder (params/Kcw_Klignin_params_v2.csv:
  per_CF2_slope 0.1 log, headgroup_offset 0.4; docs/DELIVERABLE_GAP_A_Kcw.md) and the organ lignin share.
- Prototypes DONE in scratchpad/a2 (local commits): theta(t) (691621c, b67ad93), crop sets (ca4751c:
  src/crop_sets.py, raw_si/wu2022_liaoning_oryza_params.csv verified against the PDF by the main session,
  site_scenario crop_set=, band categorical input).
- SMILES prototype DONE (copy): literature_params.kcw_pools/kcw_wholecw (GAP A ladder, +0.10/CF2 beyond the
  table); compound_from_smiles carries no transport constants; simulate_from_smiles gives each organ its K_cw.
- Prototype side effect to resolve in 2b: the neutral path's fresh mass now follows theta(t) while its K_PW water
  is constant -> tests/test_weak_electrolyte_tscf.py::test_a_small_bypass_helps_and_a_large_one_absorbs flips its
  last (ordering) assertion in the copy. The neutral composition must move to theta(t) + dry-basis lipid together.
- literature_compound/_environment/_paddy_soil: unused by the tool (tests only), carry unsourced transport defaults
  -> remove in amendment 2 and give the tests explicit synthetic compounds.
- obs_exposure_time: soil_hydrus.run_paddy_hydrus models a CONTAMINATED SOIL under CLEAN irrigation
  (top_conc 0) -- the opposite of Yamazaki (irrigation water is the source). An irrigation-source scenario
  needs: top_conc = irrigation water, initial pore water from the pre-season soil (Table S16 soil), and sourced
  hydraulics (van Genuchten thr/ths/alpha/n/Ks for the Tsukuba Andosol paddy -- now 0.08/0.46/0.016/1.37/15
  without source), root distribution (beta to 30 cm, no source), esoil 0.10 cm/d (no source), plus the found
  f_oc 3.7-3.8 %, bulk density 0.61, percolation 0.41-0.67 cm/d, drainage 35-40 DAT for 15-20 d.
  -> needs a soil-hydraulics literature search (Andosol paddy) before it can be built.

Agent report 2026-10-26 (soil side) -- for the irrigation-source exposure scenario:
- Site soil confirmed: Yamazaki main text "Andosol (volcanic ash soil) paddy field", 36 1'28.6"N 140 6'25.0"E, TC 3.8 %
  (the agent's Gray Lowland caveat concerns the NARC Kannondai paddy, a different field).
- Hydraulics: no open measured van Genuchten set for a Kanto Andosol PADDY; upland Andisol theta_s 0.68-0.75, two-step
  retention, air entry ~-15 cm (Tokumoto 2005/2007); in-situ puddled Andosol paddy Ks 9-16 mm/d by layer and percolation
  36-53 mm/d with no pan (Yoshida & Nishida 2016, doi:10.11408/jjsidre.84.10_847). Site percolation 4.1-6.7 mm/d implies a
  tighter pan whose Ks is unpublished (Tournebize 2006 NIAES Tsukuba, paywalled) -> pan Ks FITTED to the percolation range.
- Roots: ~90 % of root length within 20 cm (Beyrouty 1997), few below 35 cm (Kato & Okami 2011), Andosol japonica
  0-5/5-10/10-15 cm = 47/31/21 % (Shinoto & Otani 2020).
- ET 3.5-4.4 mm/d (Mase, Tsukuba, Ikawa 2017); T/ET ~0.85-0.9 at full canopy (isotope, secondary); pond 3-10 cm;
  Mase drainage DAT 55-73 vs Ibaraki calendar 35-40 (check which applies -> range).
- Dispersivity: 0.2-0.3 cm lab (Tokumoto 2005); cm-scale undisturbed -> range 0.2-5 cm.
- No PFAS Kd for any Japanese Andosol (Eun 2024 says none known) -> Koc QSPR x f_oc 3.7-3.8 %, run over the QSPR range.

## External review (MODEL_REVIEW_FOR_AGENT.md, read 2026-10-26; written against main 0bc6282)
Verified against the branch working tree:
- P1 (dry M in fresh equations): fixed on the branch (organ_mass dry->fresh; phloem on dry growth; xylem split on dry). Remainder = P2.
- P2 theta constant: valid -> theta(t) prototype (scratchpad/a2). Do NOT use the review's placeholder series (unsourced); use the sourced stage table.
- P3 grain mass = WSO (whole panicle incl. hulls/rachis) but composition/theta = brown rice; obs compared with brown rice: VALID, new. Plan: grain pool = brown rice mass (WSO x sourced brown-rice share), hull+rachis to a separate pool or stem/leaf with a stated basis; pre-register; report Tang/Kim/Yamazaki grain change.
- P4 leaf_loss: site-year/ORYZA path has it from DVS 0.6 (IR72.D92 DRLVT) + dead blades kept in straw (harvest operator). App measured-table path (drivers_from_tables): supplied M gets NO leaf_loss: VALID for that path.
- P5 np.interp of measured M -> stepwise mu: VALID for the measured-table path.
- P6 IR72: valid, already open (crop sets in amendment 2).
- P7 AREA_PER_HILL_M2 = 0.04 in forcing_rice/growth_rice/oryza_growth (generic scenario; also the neutral validation drivers): VALID (unsourced). Site-year uses 18.5; measured path has no default.
- P8 plaque: K_surf off, documented as not identifiable; add the explicit comparison basis (model root vs plaque-included root).
- P9 stem mass falls 13 % after heading (port; day 91 max -> harvest) -> mu<0 concentrates the stem; not documented: VALID (document; quantify; theta(t) changes it).
- P7 checked 2026-10-26: generic simulate("PFOA", scenario="generic", season=120) with AREA_PER_HILL 0.04 vs 1/18.5 ->
  every BAF identical (log10 diff 0.000000): Q_TP/hill and M/hill both scale with area/hill and cancel. Only per-hill
  absolute amounts (burden ug/hill, mass kg/hill displays) depend on it. A2: one constant, a test of the invariance,
  per-hill displays labelled with the density or shown per m2.
- Numerics (A1 runs): exponential integrator (linear_expm, 16 substeps), LAM_EQ 1e6 (1e4 not converged).

## P3 literature (agent report 2026-10-26; to verify + transcribe into sourced_inputs.csv)
- Hull share of rough rice (dry): 18.6-25.4 %, mean 20.8 % (11 varieties; japonica Hiderishirazu 20.9, Rikuto Norin 21
  20.6) -- Yoshida 1981 Fundamentals, Table 1.1 pp.6-7 (full). FAO Juliano 1993 Rice in Human Nutrition ch.3: ~20 %
  (16-28 %) of rough rice. Japonica per-grain: Nipponbare hull ~3.85 mg vs brown ~19.7 mg at maturity (~16 %; Ansari
  2016 Bangladesh Agron J 19(2):125, Fig.1 read); Hokkaido japonica hull 4.00+-0.23 mg/grain vs brown 22.5 g/1000 (~15 %;
  Mizuno 1987 JJSSPN 58(5):586, Figs 4-5).
- Hull is NOT set at heading: hull dry weight rises ~30-75 % after heading and plateaus by 10-20 DAH (Mizuno 1987: 2.3
  mg at heading -> 4.0 at 20 DAH -> 4.1 maturity, "籾殻重は出穂後20日目でほとんど最終重量に達し"; Ansari 2016: 3.0 ->
  3.8 (10 DAH) -> 3.85 -> 4.05 -> 3.85). Hull share of grain DW: ~34 % (10 DAH), ~19 % (20 DAH), ~16 % (maturity).
- Rachis share of panicle DW: NOT FOUND. Unfilled spikelets by weight: NOT FOUND (by number: 85 % filled expected,
  Yoshida Table 1.15; Korean japonica 95.1 % (90.6-97.9) ripened grains, Mo et al. 2012 AJCS 6(3):436 Table 2).
- Water: whole rough grain ~58 % early filling -> <=20 % at maturity (Yoshida p.59 Fig.1.49); brown rice 71-81 % at 4
  DAF -> 23-24 % at 20 DAF (Siddiquee 2003, indica). Hull-only water: not found.
- Hull composition: crude protein 2.0-2.8, fat 0.3-0.8, NDF 66-74, ash 13.2-21.0 (FAO Juliano 1993 Table 14, at 14 %
  moisture); Si 25.3 -> 62.7 g/kg DM heading -> maturity (Ansari Table 2); hull lignin only from Feedipedia (unattributed,
  "do not quote") -> lignin share of hull: physical bound or another source.
- Consequence for the design: grain pool = brown rice mass = WSO x (1 - hull share(t) - rachis share); rachis has no
  source -> either a range (physical) or keep it in the panicle pool with the hull. Hull pool composition sourced except
  lignin. Hull share time course from Mizuno/Ansari (japonica) as anchors (heading, 10, 20 DAH, maturity).
