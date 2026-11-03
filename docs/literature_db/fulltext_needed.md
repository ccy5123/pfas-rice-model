# Full texts needed (not obtainable from this environment)

Rule (user, CLAUDE.md §0): a paper that could not be read in full is reported to the user as a DOI
list; the user obtains it out of band. Abstract- or snippet-only values are marked as such.
Open-access copies were tried first (OpenAlex `oa_url`, J-STAGE, publisher PDF). Remove a row once
the text has been read, and record what was taken from it in `sourced_inputs.csv`.

## Open

| DOI | paper | needed for | why it failed |
|---|---|---|---|
| 10.1021/acs.est.0c05705 | Nguyen et al. 2020 ES&T | per-soil, per-congener log Kd, AEC, oxalate Fe/Al, pH series (PREREG 10.1) | ACS 403 (Cloudflare); no PMC copy |
| 10.1016/j.chemosphere.2023.138133 | Campos-Pereira et al. 2023 Chemosphere | Kd vs pH, oxalate Fe/Al (10.1) | Elsevier 403; SLU/DiVA repository failed |
| 10.3390/environments10100175 | Brusseau 2023 Environments | short- vs long-chain Kd vs soil properties (10.1) | MDPI 403; repository 403 |
| 10.1016/j.scitotenv.2021.149343 | Fabregat-Palau et al. 2021 STOTEN | Koc and Kmin coefficients by chain length (10.1) | Elsevier 403 |
| 10.1071/SR22183 | Kookana et al. 2022 Soil Research | soil properties governing PFAS sorption (AEC, oxides) (10.1) | CSIRO blocked |
| 10.1016/j.scitotenv.2018.01.167 | Li, Oliver & Kookana 2018 STOTEN | compiled Kd dataset (10.1) | Elsevier 403 |
| 10.1017/S0043174500056575 | Bromilow, Chamberlain & Evans 1990 Weed Sci | log Kow-pKa map of phloem mobility; polar end for neutrals (10.3) | Cambridge paywall |
| 10.1017/S0043174500056587 | Hsu & Kleier 1990 Weed Sci | optimum log Kow vs membrane permeability (10.3) | Cambridge paywall |
| 10.1093/jxb/47.Special_Issue.1265 | Hsu & Kleier 1996 J Exp Bot | review; any stated cut (10.3) | OUP paywall |
| 10.1002/ps.2780300108; 10.1002/ps.2780420102 | Grayson & Kleier 1990; Kleier 1994 Pestic Sci | phloem mobility IV, V (10.3) | Wiley paywall |
| 10.1017/s0043174500094637 | Kleier 1996 Weed Sci | phloem mobility VII (10.3) | Cambridge paywall |
| 10.1007/bf00024777 | Brudenell 1995 Plant Growth Regul | tabulated properties vs mobility (10.3) | Springer paywall |
| 10.1007/BF00239934 | Grimm et al. 1995 Planta | polar compounds retained/leaking in Cyclamen petiole (10.3) | Springer bot challenge |
| 10.1071/fp13267 | Zhivotovsky Castello 2014 Funct Plant Biol | xylem discontinuity at the cereal grain (grain route, item 4) | no open copy |
| -- | Iwatsuki 1932, 農業及園芸 7:64 | seasonal rice root profile to 89 cm (10.2) | no DOI; cited only |

Leads, not yet requested (ask before relying on them):
- Kobayashi & Horie 1994, Jpn. J. Crop Sci. 63:193-199 (spikelet and rachis-branch differentiation) —
  possible rachis weights for PREREG 7.2 (cited by Matsui & Kagata 2002).

## Obtained (2026-10-26, uploaded by the user or open access) — read in the main session
- 10.1626/pps.5.198 Matsui & Kagata 2002 (no rachis weight: spikelets only)
- 10.3839/jabc.2018.028 Choi et al. 2018 (husk vs brown rice PFOA/PFOS; `raw_si/choi2018_husk_brown.csv`)
- 10.1007/BF00387436, 10.1007/BF00385354 Oparka & Gates 1981a,b (caryopsis water and assimilate routes)
- 10.1071/BI9720411 Zee 1972 (spikelet xylem)
- 10.1007/BF02705120 Krishnan & Dayanandan 2003; 10.3724/SP.J.1006.2010.01198 Yu et al. 2010
- 10.1016/j.ecoenv.2022.113683 Luo et al. 2022 (Cd by panicle organ)
- 10.1007/s10333-005-0028-0 Tournebize et al. 2006 (PCPF-SWMS; for amendment 2b soil — not yet read)
- 10.1021/acs.est.5c11716 Muensterman et al. 2026 ES&T 60:1207 + correction 10.1021/acs.est.6c07668 (eq. 4:
  log kIAM = 0.046 CHI_IAM + 0.42) — biomimetic K_PL / K_prot for PFAS; not yet read (binding work, not 2a)
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
