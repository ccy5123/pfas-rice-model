# Full texts needed (not obtainable from this environment)

Rule (user, CLAUDE.md §0): a paper that could not be read in full is reported to the user as a DOI
list; the user obtains it out of band. Abstract- or snippet-only values are marked as such.
Open-access copies were tried first (OpenAlex `oa_url`, J-STAGE, publisher PDF). Remove a row once
the text has been read, and record what was taken from it in `sourced_inputs.csv`.

## Open

| DOI | paper | needed for | why it failed |
|---|---|---|---|
| -- | Iwatsuki 1932, 農業及園芸 7:64 | seasonal rice root profile to 89 cm (10.2) | no DOI; cited only |

## Not obtainable (user, 2026-11-04) -- the consequence stays as recorded
| item | needed for | what stays |
|---|---|---|
| Sörengård 2020 Fig. 1 raw values of the hidden markers (authors) | PFBS (3 of 4 sorbents) and PFBA (2) K_cw lower end | lowest READ value, flagged `min_may_be_lower` (PREREG 11.7) |
| Eun et al. 2022 Fig. 3 per-congener drainage water (authors) | in-season exposure of C4-C8 | `obs_exposure_time` open (11.5) |
| IRRI ORYZA download folder / Jing & Bouman files (Nipponbare, Takanari) | Japanese-cultivar crop set | IR72 runs, published japonica sets as band ends (11.6, 11.8) |
| AgMIP-Rice Akitakomachi parameters (Hasegawa) | same | same; DSSAT `Shizukuishi.CRP` stays out (undocumented) |

Leads, not yet requested (ask before relying on them):
- (Kobayashi & Horie 1994 read 2026-11-04: counts only, no rachis weight.)

## Obtained 2026-11-04 (user, Google Drive `papers (5).zip`, sha256 466246d847f6cbfafc328e2c7763aeaf669d2888ae081c3fbec4b2bc93840441)
Read in full (Sörengård by the main session; the rest by literature subagents, PREREG 11.7-11.8). PDFs not committed.
- 10.1016/j.jece.2020.103744 Sörengård 2020 (sorengard2020.pdf 4e614b14b55df97f36f591cb58dcaf53b2823d693222172cacaee71abfea2928; SI open,
  Elsevier mmc1.docx 0ff4b962118d0cdbf070db0c6164eb44be0da252cdf93d28a7211dd89c75e168) -> `raw_si/sorengard2020_fig1_all_compounds.csv`
- 10.1016/j.agwat.2025.110041 Zhang 2026 (zhang2026.pdf 4dc3dcebeec18377be886461617b2287d2b3e78dc6b169a18e4911c32a978c96) -> `raw_si/zhang2026_*.csv`
- 10.1016/j.scitotenv.2024.177377 Wu 2024 (5255230b23ba26c4d89b58a3b307d8cf7649924811cf9ce2c95aab270e002f18) -> `raw_si/wu2024_text_values.csv` (SI not obtained)
- 10.1016/j.envpol.2026.127828 Wu 2026 (b25e36762e69b09579b6ef4283ba9bafb5556b7bddb803fc4451f2c507ba9364) -> `raw_si/wu2026_text_values.csv` (SI not obtained)
- 10.1021/acs.est.4c03525 Chen 2024 (22e15eeac6d709ba8d1523993a797333005e572803b148495b6a9fc73a7af2c7) -> `raw_si/chen2024_fig2_stage_labels.csv` (SI not obtained)
- 10.3791/64400 Liu, Zhou & Zhu 2022 (b1d18cd92581b100d394b5ea33228908585d910f45b8b04f1749f70648e82905) -> `raw_si/liu2022_jove_splitroot_text.csv`
- 10.1016/j.fcr.2004.10.003 Bannayan 2005 (0fe9a7895047615b8a8bdc53b43436ebf66ef864b8343d7298a809d06e8ac18b): no parameter values printed
- 10.1016/j.agsy.2008.06.005 Jing 2008 (619e6cff94f8b749d28ab9fd203a0d8ca30c365424f6839cff885a62686e8811): no values ("available from the authors")
- 10.1111/gcb.12758 Li 2015 (li2015.pdf 3696889ad78c839e5eddec64d892f06455a01c8912439fa92391a7e74c6c7708; SI docx
  0658b14f366397edb711c1120a3a69138226fee120b33bf2cf07579d19c3f442) -> `raw_si/li2015_agmip_s3_oryza_params.csv`, `shizukuishi_crp_vs_li2015_s3.csv`
- 10.1002/agj2.20258 Lu 2020 (f7f1b27753f77a635d27d84e9b87c0e20d1cc82dc957a87c0476cdc8b64b3314) -> `raw_si/lu2020_jiangsu_oryza_v3_params.csv`
SIs obtained 2026-11-04 (user, `papers_SI.zip` sha256 967b5d2e72e13176063f74c018dcb93345e2600d8e553d894e7cdf66ad3147f9; read by a subagent):
Wu2024_SI1.docx b70489c12cedf4c5ed25b45c0aade663ff33853651e2e85f33c42a333766799e; Wu2026_SI1.doc 55cfbe7fd555b30f35b72e0b7017d75d174de2f87fc887455d6d38c127a0a54f;
Chen2024_SI1.pdf d5cec0ae37d47b892af906731f3bd6495707057bb7af1320851bbe5371e1f7b4; Zhang2026_SI1.docx 6124f8df9fc52001503dd38f00804a67d524fd6aa5cc278f61e3b99327a20302.
Also read: 10.1626/jcs.63.193 Kobayashi & Horie 1994 (open access J-STAGE scan, 2e18480c9866ab434b082013d5bd633443c7744f8822b251fc2cec4d4422fbee; no SI exists):
no rachis dry weight. Nothing in these changes an open question (PREREG 11.11).

## Obtained (2026-10-26, uploaded by the user or open access) — read in the main session
- 10.1626/pps.5.198 Matsui & Kagata 2002 (no rachis weight: spikelets only)
- 10.3839/jabc.2018.028 Choi et al. 2018 (husk vs brown rice PFOA/PFOS; `raw_si/choi2018_husk_brown.csv`)
- 10.1007/BF00387436, 10.1007/BF00385354 Oparka & Gates 1981a,b (caryopsis water and assimilate routes)
- 10.1071/BI9720411 Zee 1972 (spikelet xylem)
- 10.1007/BF02705120 Krishnan & Dayanandan 2003; 10.3724/SP.J.1006.2010.01198 Yu et al. 2010
- 10.1016/j.ecoenv.2022.113683 Luo et al. 2022 (Cd by panicle organ)
- 10.1007/s10333-005-0028-0 Tournebize et al. 2006 (PCPF-SWMS; for amendment 2b soil — not yet read)
- 10.1021/acs.est.5c11716 Muensterman et al. 2026 ES&T 60:1207 + correction 10.1021/acs.est.6c07668 (eq. 4:
  log kIAM = 0.046 CHI_IAM + 0.42) — biomimetic K_PL / K_prot for PFAS; read 2026-10-27 with its SI (below)
- 10.3329/baj.v19i2.31861 Ansari et al. 2017 (open access; hull and brown-rice weight in Fig. 1, to digitise)
- 10.1626/pps.15.258 Fukuoka et al. 2012; 10.1007/s11270-007-9425-x Rahman et al. 2007 (open access)

## Obtained (2026-10-27, uploaded by the user or open access) -- read in the main session or by subagents
- Trapp, Matthies & McFarlane 1994 (10.1002/etc.5620130308); Trapp & Matthies 1995 (10.1021/es00009a027) and its
  correction 1996 (10.1021/es951008p); Trapp 2004 (10.1065/espr2003.08.169) -- neutral constants (PREREG 9.1)
- Felizeter 2012 / 2014 main texts; Tang 2026, Kim 2019, Li 2025 main texts (`tr_drivers`)
- Tournebize 2006 (read this time); Matsumoto, Yoshida & Nishida 2016 (10.11408/jjsidre.84.10_847)
- Tokumoto 2005 (10.3178/jjshwr.18.401, open access); Tokumoto 2010 (10.1016/j.coldregions.2010.03.005);
  Rudiyanto 2013 thesis (hdl.handle.net/10076/13878, open access); Moraes 1973 (Embrapa, open access)
- Yamazaki 2023 SI (figshare, open access): Tables S5, S16 newly transcribed
- Chang 2020 (10.1186/s13007-020-00633-1), Meera 2019 (10.1007/s13197-019-03631-x), Setter 1995, Dharni 2022 -- panicle area
- Malaviya 2025 (10.1139/er-2025-0116): the accepted manuscript (R1) is open at hdl.handle.net/1807/151243
- HYDRUS-1D manuals 4.17 (sha256 286f4aedaf30064fbc723aba3f0e7f9e9b45c1ed933b04b184888a8bd9462755) and 4.08
  (sha256 20c29f1bb6d0567cb5f63a725a877c9d5aa277ed284a6255c718730f0fc03847), pc-progress.com, open: the air-entry
  model 3 (`soil_hentry`).
- Uploaded by the user 2026-10-27 (papers_6.zip), read in the main session: Trapp 2007 (10.1080/10629360701303693;
  sha256 a34d2a7f...), Rein, Legind & Trapp 2011 (10.1080/1062936X.2010.548829; sha256 943a0511...) -- the neutral
  phloem sap (PREREG 9.15); Vogel, van Genuchten & Cislerova 2001 (10.1016/S0309-1708(00)00037-3; sha256 b007f8f9...)
  -- `soil_hentry`. The same upload's Tang 2026, Kim 2019 and Felizeter 2012/2014 files are the main texts
  (sha256 fb2dbd93..., c2faf194..., c70492c6..., d79aa6f8...); their SIs are still needed (table above).
- SIs uploaded by the user 2026-10-27 (papers_SI.zip), read in the main session: Tang 2026 SI (docx, sha256 deb3d35b...),
  Kim 2019 SI (docx, sha256 352711e3...), Felizeter 2012 SI (sha256 9270ca0e...), Felizeter 2014 SI (sha256 fe7bdafb...).
  None of them gives a crop calendar for Tang or Kim (`tr_drivers` cannot be resolved from the papers); Felizeter's SIs give
  dates and tissue tables but no TSCF table. Transcriptions: `raw_si/tang2026_si_table_s2_soil.csv`,
  `raw_si/kim2019_si_table_s2_sites.csv`, `raw_si/felizeter2014_si_table_s3_dates.csv`; checks: Tang S7/S8 (75 cells; one SD
  0.025 -> 0.0254 corrected) and Kim S5/S6 (51 rows, all equal) against the repo's earlier transcriptions.
## Read 2026-11-03 (open access; amendment 3, PREREG 10)
- Yamazaki 2023 SI re-read for the extraction method (p.S6; sha256 e2e53c9a..., same file as before).
- Morita, Suga & Yamazaki 1988 (10.1626/jcs.57.438; sha256 8c2e5493...; main session, Table 2), Kang 1994
  (10.1626/jcs.63.118; 400a1206...), Tanaka 1985 (10.1626/jcs.54.379; 18aab139...), Kawata 1980 (10.1626/jcs.49.311;
  7930828e...), Arai-Sanoh 2014 (10.1038/srep05563; 514d66e2...), Kawata 1963 (10.1626/jcs.32.163; 70414c82...) -- roots.
- Tyree 1979 (10.1104/pp.63.2.367; d37f7afc...), Kleier 1988 (10.1104/pp.86.3.803; 1353450f...), Hsu, Kleier &
  Melander 1988 (10.1104/pp.86.3.811; 3e07b090...) -- phloem mobility.
- Campos-Pereira 2020 (PMC7745537), Oliver 2020 (accepted MS), Fabregat-Palau 2025 (PMC12020356), Evich 2025
  (PMC11962083), Eun 2022 (PMC9408263), Eun 2026 (10.5917/jagh.68.59), Yamazaki 2024 (10.2116/bunsekikagaku.73.31),
  Endo 2026 (10.5985/jec.36.s9) -- PFAS sorption (subagent; no Andosol Kd in any).

SHA-256 of every file: in the session record (PREREG 9 commit message).

Requested 2026-11-03 for items 4-5 (PREREG 10.5):
| DOI | paper | needed for | why it failed |
|---|---|---|---|
| -- | ORYZA v3 Nipponbare / Takanari crop files (Espe 2016 p.126) | a Japanese-cultivar crop set | IRRI download folder unreadable |

New open (2026-11-03, after the user's papers (4).zip):
| DOI | paper | needed for | why |
|---|---|---|---|
| -- | Iwatsuki 1932, 農業及園芸 7:64 (岩槻信次「水稲の根に関する研究」) | root profile (pot study, per Kawata 1963) | no DOI; low priority |

## Obtained 2026-11-03 (user, Google Drive `papers (4).zip`, sha256 ddb7bff0c31867f341194589631155fa2a3f9911ab5f72004a34a506e52c1291)
Read in the main session (values set) or by subagents (records): bromilow1990 2f9c0d53..., hsu1990 c856b6b3..., hsu1996
c0a6bd63..., kleier1994 5eb01238..., kleier1996 f18e0cca..., grayson1990 fe44eb98..., brudenell1995 1ab4b058..., grimm1995
0daeda7a..., zhivotovskycastello2014 f785d89c..., hayashi1990 87bb6899..., fukumorita1982 c55a8e0f..., tcherkez2023 e718b07d...,
jiang2003 ab219215..., pate1983 b0c4d755..., yamaguchi1990 7b24de2c..., nguyen2020 89f92567..., campospereira2023 0e7344b0...,
brusseau2023 67782d65..., fabregatpalau2021 1682c066..., kookana2022 03e3480b..., li2018 f62352fb..., yamazaki2019 3f39d762...,
wen2014 3ee40c99..., attia2026 d51ed3c8..., krippner2015 a19ccff6..., stahl2009 3539c5f3..., Liu 2019 SI (docx) b637df14....
- SIs uploaded by the user 2026-11-03 (papers_SI.zip, sha256 5f98567b...): Nguyen 2020 SI 9ea6c31c..., Fabregat-Palau 2021 SI
  f23d9409..., Campos-Pereira 2023 Appendix A a21915e5..., Yamazaki 2019 SI (docx) 6d545c00... -- read by a subagent (PREREG 10.7d).
## Obtained (2026-10-27, uploaded by the user) — read in the K_cw side-project session
Main texts and SI read; extracted values in `raw_si/` (file named per entry); reading log in
`docs/kcw_side_project/READING_LOG_KR.md`.
- 10.1016/j.jhazmat.2024.136016 Mel et al. 2024 (lignin Kd, 4 PFAS; SI Table S1 is an image) → `raw_si/mel2024_lignin_kd.csv`
- 10.1021/acs.est.5c11716 Muensterman et al. 2026 + 10.1021/acs.est.6c07668 (SI Tables S12, S13) → `raw_si/muensterman2026_tableS12_units.csv`, `raw_si/muensterman2026_tableS13_kiam.csv`
- 10.1016/j.scitotenv.2020.139383 Wang et al. 2020 (*Alisma*; SI Tables S10, S11, S13, S15) → `raw_si/wang2020_table*.csv`
- 10.1021/acs.est.3c00504 Liu et al. 2023 (wheat; root fractions only in main-text Fig. 3a, to digitise; SI has no table of them)
- 10.1021/acs.est.1c01175 Yu et al. 2021 (lettuce PFOS; fractions in Fig. 1 and SI Fig. S6, to digitise)
- 10.1016/j.envint.2021.106642 Mei et al. 2021 (review; no cell-wall component Kd — not the K_cw_poly source)
- 10.1038/s41467-025-65191-3 Guo et al. 2025 (open access; Source Data Fig. 3b) → `raw_si/guo2025_fig3b_subcellular.csv`

## Not found
- The "Prompt2 handoff" that came with `params/Kcw_Klignin_params_v2.csv` (source of the K_cw_poly absolute scale, the
  0.40 sulfonate offset and the 0.10/CF2 slope): not in the repo, its history, the importing session's transcript,
  Google Drive or Notion; the user could not obtain it (2026-10-27). K_cw_poly is untraceable until a source is found.

## Obtained (2026-10-27, second batch, uploaded by the user) — read in the K_cw side-project session
- 10.1016/j.chemosphere.2022.134167 Campos-Pereira et al. 2022 (main text; SI not obtained) → `raw_si/campos2022_table1_nmr.csv`
- 10.1016/j.jece.2020.103744 Sörengård et al. 2020 (main text; SI not obtained) → `raw_si/sorengard2020_fig1_plant_sorbents.csv`
- 10.1016/j.emcon.2022.01.001 Mokwena & Okonkwo 2022 (main text): raw maize tassel at pH 2 and mg/L — not usable for K_poly; nothing extracted
- Makselon 2019 SLU MSc thesis (open access, stud.epsilon.slu.se/15180) → `raw_si/makselon2019_peatOi_*.csv`
- 10.1139/er-2025-0116 Malaviya et al. 2025, Environ. Rev. — open accepted manuscript found 2026-10-27 (University of Toronto
  repository, utoronto.scholaris.ca bitstream b9636f7f-77fd-4aef-b9f9-eb82fef09469; SHA-256 65ae5c0a59fb26ec7fa3d878baddbe30e930811c587919c488d4b00145e28091).
  Checked only for cell-wall K data (none). NOT yet read for its purpose here (the long-chain lipid mechanism); not the version of record.
- hdl 2376/111903 Tian 2019, WSU PhD dissertation "Novel adsorbents for per- and polyfluoroalkyl substances" (uploaded by the user to Drive,
  same file as the repository copy, SHA-256 1f57b8a82043ffb03623f544eb64856a8fc815469b237685e5828797073e08f3; read in full by a subagent, key
  sentences re-checked by the main session). No cellulose at all; unmodified Kraft lignin (Indulin AT) showed "no adsorption" of PFOA (§4.3.2,
  p. 116) and "nearly no adsorption" of PFBA/PFBS (§5.3.2.1, p. 144) at 50 mg/L sorbent — a non-detect that bounds Kd only to < ~1000 L/kg.
  Not the source of Mel 2024's "quaternized cellulose" values; nothing extracted to raw_si (no numbers for a native material).
