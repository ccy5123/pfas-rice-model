# Full texts needed (not obtainable from this environment)

Rule (user, CLAUDE.md §0): a paper that could not be read in full is reported to the user as a DOI
list; the user obtains it out of band. Abstract- or snippet-only values are marked as such.
Open-access copies were tried first (OpenAlex `oa_url`, J-STAGE, publisher PDF). Remove a row once
the text has been read, and record what was taken from it in `sourced_inputs.csv`.

## Open

| DOI | paper | needed for | why it failed |
|---|---|---|---|
| 10.1016/j.jhazmat.2025.141017 (SI) | Tang et al. 2026 J Hazard Mater, Supplementary Tables S1-S8 | soil properties (Table S2), the experiment's dates (`tr_drivers`, PREREG 9.4) | main text read; SI not supplied |
| 10.1016/j.scitotenv.2019.03.240 (SI) | Kim et al. 2019 Sci Total Environ, SI Table S2 | site regions, distances, soil and water properties (`tr_drivers`) | main text read; SI not supplied |
| 10.1021/es302398u (SI), 10.1021/jf500674j (SI) | Felizeter et al. 2012 / 2014, SI | concentration, growth and transpiration tables; experiment dates | main texts read; SI not supplied |
| 10.1080/10629360701303693 | Trapp 2007, "Fruit Tree model for uptake of organic compounds from soil and air", SAR QSAR Environ Res 18:367 | the neutral phloem-sap concentration (PREREG 9.1c) | closed access |
| 10.1080/1062936X.2010.548829 | Trapp 2011 (Rein, Legind & Trapp), "New concepts for dynamic plant uptake models", SAR QSAR Environ Res 22:191 | same | closed access |

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
SHA-256 of every file: in the session record (PREREG 9 commit message).
