# Task 3: ORYZA crop sets for Japanese cultivars (search 2026-11-04, subagent)

## Obtained file
- `dl_shizukuishi/Shizukuishi_4.6.CRP` (identical copy `Shizukuishi_4.7.5.11_prev.CRP`)
  - SHA-256 baf76586a2e13ecee8a7a4548044bc2e63cd8057b15b28cbbb4f984303a1e502
  - URL: https://raw.githubusercontent.com/allfed/mink/main/basics_15jun22/sge_Mink3daily/actual_program_4.6_whpt_08jul16/not_potatoes_stuff/Shizukuishi.CRP
    (same file under actual_program_4.7.5.11/previous_CULetc/)
  - Context: DSSAT 4.6/4.7 ORYZA cultivar file RIORZ046/047.CUL lists `IB6002 IR72-Shizukuishi . Shizukuishi.CRP`
    (also in agmip/gbuilder2d Data/DSSAT47/Genotype/RIORZ047.CUL). EXPNO blank; no reference given.
  - Header still the IR72.D92 header (copied). CO2 = 368 ppm (Shizukuishi FACE ambient, cf. Hasegawa 2017's
    367). Likely the AgMIP-Rice Shizukuishi (Akitakomachi) calibration, but that is an inference: no document names it.
  - vs IR72.D92 (repo raw_si/oryza_ir72_d92_params.csv): UNCHANGED TBD/TOD/TMD, DVRJ/I/P/R (0.000773/0.000758/
    0.000784/0.001784 = IR72), SHCKD, SHCKL, SSGATB, SPGF 64900, WGRMX 2.49e-5.
    CHANGED: RGRLMX 0.00955 (IR72 0.0085), RGRLMN 0.00462 (0.004); SLATB table 0:0.00470;0.42:0.00222;0.65:0.00219;
    0.80:0.00219;1.0:0.00219;1.2:0.00195;2.0:0.00142;2.5:0.00053; KDFTB 0:0.462;0.65:0.541;1.0:0.632;2.5:0.671;
    EFFTB 0:0.523;10:0.633;40:0.435;60:0.185; FSTR 0.0917 (0.20); FSHTB 0:0.582;0.42:0.954;1:1;2.5:1;
    FLVTB 0:0.614;0.42:0.351;0.65:0.448;0.80:0.235;1.0:0;...; FSTTB 0.002:0.386;0.42:0.649;0.65:0.552;0.80:0.632;
    1.0:0.299;1.2:0.045;2.0:0; FSOTB 0:0;0.42:0;0.65:0;0.80:0.133;1.0:0.701;1.2:0.955;2.0:1;
    DRLVT 0:0;0.6:0;1.0:0.0030;1.6:0.0056;2.0:0.0163;2.5:0.0224; plus COLDMIN 5, COLDEAD 8, SHADET 0.918,
    ZRTMCD 0.56, FSWTD 0.346, N params (FNTRT, NFLVI, FNLVI, NMAXSO, NMAXLT). FLV+FST+FSO = 1 at every DVS.
  - DSSAT's own IR72.crp (same repo, sha256 c0e7ce82...) is itself an ORYZA v3 re-tune, not IR72.D92.

## Read in full, no usable values
- Jing Q 2007 WUR PhD thesis (10.18174/31932; edepot.wur.nl/31932; sha256 8f8b3280...): ch.3 calibrates ORYZA2000
  for Nipponbare, Takanari, IR72 (Iwate, Nagano, Shimane, Kyoto, Nanjing, Chiangmai, Ubon, Taoyuan; 2002 data):
  DVR, SLA, FLV/FST/FSO calibrated, rest IR72; "A complete list of all model input parameter values is available
  from the authors" (p.32). Same for Wuxiangjing9 (ch.4, p.53). = the IRRI "calibrated varieties" page entry, and
  almost surely the Nipponbare/Takanari files Espe 2016 used. Published as Jing et al. 2008 Agric Syst 98:177
  (10.1016/j.agsy.2008.06.005).
- Espe 2016 FCR 193:123 (sha256 69ee4867...): p.126 only says the Nipponbare/Takanari files were starting points.
- Hasegawa 2017 Sci Rep 7:14858 + SI (sha256 4729a289..., 0772803c...): no parameter values; parameters
  "reported to the AgMIP-Rice team leader".
- Li T 2015 GCB 21:1328 (igsnrr OA copy): Shizukuishi = Akitakomachi; values in Table S3 (calibration method) at
  Wiley (403).
- Li T 2017 AFM (ORYZA v3; Europe PMC XML): no Japanese cultivar.
- Hyun & Kim 2019 KJAFM 21:42 (Korean japonica Shindongjin; sha256 ea90bdb2...): DVR quantiles (Table 3-4) and
  growth-parameter quantiles (Table 5) from yield-only data; authors say the yield fit is unreliable (R2 < 0.1).
- APSIM Classic Oryza.xml: Australian japonica (amaroo, langi, quest, riverina): DVR/MOPP/PPSE only.
- IRRI calibrated-varieties page: lists IR72/Nipponbare/Takanari (China, Japan, Thailand) with no values.

## Not obtained
- IRRI ORYZA Google Drive folder 0BwRcFC9_t3sQVnNRaXZMM3BYelk (401; Drive MCP: not found); v2.13 / v3 installers
  behind Google Forms.
- Bannayan et al. 2005 FCR 93:237 (10.1016/j.fcr.2004.10.003): ORYZA2000 Akitakomachi Shizukuishi FACE; SLA and
  root:shoot corrected (conference-abstract snippet only); closed.
- Jing et al. 2008 Agric Syst 98:177 (10.1016/j.agsy.2008.06.005) (thesis chapter; likely the same "from the authors").
- Li T 2015 GCB SI Table S3 (Wiley 403).
- Lu B et al. 2020 Agron J (10.1002/agj2.20258): Jiangsu japonica, ORYZA v3 potential production (Wiley 403).
- Paleari et al. 2022 GCB (10.1111/gcb.16087): AgMIP-Rice trait ensemble, unimi copy restricted.

## Null searches
- GitHub code: 'Koshihikari ORYZA', 'Koshihikari DVRJ', 'Nipponbare CROP DVRJ', 'Takanari DVRJ', 'Kasalath DVRJ',
  'japonica DVRJ RGRLMX', 'Nipponbare/Takanari/Koshihikari/Akitakomachi WGRMX', 'Nipponbare FSOTB',
  'Koshihikari SLATB', filename:nipponbare.crp / takanari.crp / NIPPONBARE.DAT / takanari.dat, 'Akitakomachi DVRJ',
  filename:Nanjing.CRP: 0 relevant hits. All 45+71 hits for 'DVRJ extension:crp' / 'DVRJ RGRLMX FLVTB' listed;
  the only Japanese-site file is Shizukuishi.CRP.
- J-STAGE API (ORYZA2000, ORYZA, ORYZA v3): no Japanese ORYZA calibration paper (one CSSJ abstract on Colombian
  varieties). CiNii OpenSearch 502. Web (EN/JA): no Koshihikari/Kasalath ORYZA set.
