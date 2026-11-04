# Task 2 notes — literature search 2026-11-04 (subagent; repo read-only)

## (A) Seasonal exposure of C4-C8 at the Yamazaki 2023 site

### Read in full
1. Eun, Yamazaki, Pan, Taniyasu, Noborio, Yamashita 2022, IJERPH 19:10379, 10.3390/ijerph191610379 (PMC9408263).
   Europe PMC full-text XML sha256 9f1f2021383f321b9d514796a850c3ffccc8935c0d7c5e16141718738daea3a9; figures from mdpi-res.com.
   Listed in the repo as "Eun 2022 (PMC9408263)" under the 2026-11-03 sorption subagent, but NOTHING was extracted from it.
   - Site: Meiji Univ. paddy lysimeter (Kawasaki), 2 x 2 m, 1.8 m deep; low-humic Andosol 0-35 cm over Kanto loam; clear roof
     (no rain); 16.7 L spike of simulated contaminated water after puddling (21 Jun 2015), then non-drinking tap water.
     Leachate from drains at 5/40/75/110/145/180/240 cm, May-Nov 2015 (Fig. 3 dates 18 Jun - 8 Oct). It is the same
     lysimeter programme as Yamazaki 2019, NOT the 2023 Tsukuba field.
   - Table 1 (read): tap irrigation water pg/L PFBS 352, PFBA 858, PFPeA 391, PFHxA 647, PFHpA 528, PFOA 1166,
     PFNA 498, PFHxS 1721, PFOS 1240; sum PFAS 7634 pg/L. Spike ng/L PFBA 14,754, PFPeA 15,100, PFHxA 7873, PFHpA 3198,
     PFOA 19,156, PFBS 1484, PFHxS 702,885, PFOS 395,909.
   - No percolation rate, no water balance, no Kd, no TC by layer beyond "~4 % upper, 2-3 % lower", no tables of soil or
     leachate values: soil by depth only in Fig. 2 (pre/post, 5-145 cm), leachate only in Fig. 3 (stacked bars, legend
     colours repeat, so per-congener values cannot be read) and Fig. 4 (PFHxS/PFOS).
   - By-eye figure readings (NOT table values, +-20 %): Fig. 3, 5 cm leachate before the spike (18-19 Jun) sum PFAS
     ~100 and ~60 ng/L, i.e. ~10x the tap irrigation water (7.6 ng/L); 240 cm on 18 Jun ~40 ng/L. Fig. 2 pre-experiment
     PFBA/PFPeA/PFHxA/PFHpA/PFOA present from 5 to ~85 cm at 100s of pg/g dw (resident pool through the profile).
     Post-experiment short chains at 5 cm are higher than pre (the spike), so the profile cannot show seasonal depletion.
   - Text p. Results: "shorter chains (C <= 4), such as PFBA and PFPeA, showed insignificant fluctuations".

2. MAFF (農林水産省) 2025, 「農業環境中のPFOS、PFOAのコメへの移行、蓄積性について」 (FY2024 results).
   https://www.maff.go.jp/j/syouan/seisaku/PFAS/attach/pdf/index-9.pdf sha256 a8bb800887c5fbc49df11564fe68430ce675991c1ab1f96c7c58d276329fad05
   15 paddy plots (region withheld) irrigated from rivers over the provisional guideline; soil 0-15 cm before tillage and at
   harvest; brown rice at harvest. PFOS soil 640-6600 ng/kg dw, brown rice < 9 ng/kg (n 10); PFOA soil 510-18,000,
   brown rice < 2-20 ng/kg (median 5, n 63); brown/soil max PFOS 0.005, PFOA 0.004 (p.2-3). Linear isomers only.

3. MAFF regulatory-science project 「農産物中PFASの分析法の確立、農地土壌、水等からのPFAS移行特性の解明」 (FY2023-2027),
   cumulative progress report https://www.maff.go.jp/j/syouan/seisaku/regulatory_science/attach/pdf/kadai_jisshi-53.pdf
   sha256 61ce7b51e985497a5c0be211524775d176363c47e308666b024b14f1c8142593 (overview kadai_jisshi-9.pdf sha256 724347a8...).
   FY2024 section (report p.16-24; pdftotext lines 1702-1980): 30 PFAS, paddy A (13 plots, PFOA-dominated), paddy B (2 plots,
   PFOS-dominated); soil 0-15 cm 5 points/plot spring (before) and autumn (after); paddy SURFACE water (田面水) at the outlet,
   one grab in early season and one in summer; brown rice at harvest. Soil type and region not given.
   - Soil stable over the season (text p.17, Fig. 6, Fig. 13 caption: "no large change in species or concentration between
     spring and autumn"): c-1 PFOA 14,076 -> 15,704 ng/kg dw, PFOS 462 -> 510, PFNA 1,084 -> 1,230, PFHxS 73 -> 57;
     B-a PFOS 4,514 -> 5,070, PFOA 292 -> 359, PFNA 951 -> 1,006, PFHxS 1,276 -> 946. Short-chain soil values only in stacked
     bars (too thin to read).
   - Surface water c-1: spring PFOA 2,812, summer 391 ng/L (channel works upstream); "soil showed no clear effect; PFAS
     strongly adsorbed and retained in the soil" (p.18).
   - B-a surface water spring -> summer: PFHxS 43 -> 44, PFOS 70 -> 65, PFOA 9 -> 13, PFNA 28 -> 25 ng/L; short-chain PFCAs
     ~2x, PFBA 10 -> 44 ng/L (p.22); B-b PFBA 9 -> 29 (Fig. 13, bar label). "reason unknown, few samples".
   - Brown rice: PFBA dominant (A: median 86, max 325 ng/kg; B: median 86), PFHxA (A median 20), HFPO-DA (A median 30).
   - p.21: the authors state the uptake should be assessed on the dissolved PFAS in soil water ("今後明らかにすべき課題") —
     no pore water measured.
   FY2025 section: upland crops only (no paddy).

### Abstract-only (full text not obtained)
- Zhang J., Hou R., Su A., Wang M. et al. 2025 Agric Water Manag, 10.1016/j.agwat.2025.110041 (OA per OpenAlex):
  MIP-DGT + HR-Peeper pore water, 2D-DIFS fit of PFBA and PFOA (R, Tc, k-1, Kd), resupply R declines with irrigation
  cycles; surface C_DGT PFBA 6.76 +- 1.44, PFOA 7.21 +- 2.14 ng/L (ABSTRACT). ScienceDirect 403; SSRN preprints
  10.2139/ssrn.5157625 / 5125235 "Content Blocked". Soil not Andosol (Chinese paddy; type unknown).
- Wu J.-y. et al. 2024 STOTEN 10.1016/j.scitotenv.2024.177377 (irrigation-induced N-Fe redox; C4-C9 released into overlying
  water) — closed; snippet only.
- Wu J. et al. 2026 Environ Pollut 10.1016/j.envpol.2026.127828 (rice residue decomposition; pore and overlying water as
  sinks, drained vs waterlogged) — closed; snippet only.
- Wu J.-y. et al. 2025 Environ Pollut 10.1016/j.envpol.2025.127193 (quinolones; SWI redistribution) — closed.
- Wu J.-y. et al. 2024 Environ Pollut 10.1016/j.envpol.2024.125580 (microbial contribution to SWI migration, waterlogged
  paddy) — closed.
- Li F. et al. 2025 Acta Agric Univ Jiangxiensis 10.3724/aauj.2025147 (short-chain immobilisation by amendments, lab
  incubation; OA PDF on sciengine) — abstract read only; not a field series, low priority.
- Yamashita, Yamazaki, Taniyasu, Hanari, Yeung 2024 Chemosphere 10.1016/j.chemosphere.2024.143073 — activated carbon /
  biochar; abstract cites "an earlier lysimeter study" (= Eun 2022). Not fetched (ScienceDirect blocked); not needed.

### Null searches (A)
- OpenAlex author works since 2022-06 for Eriko Yamazaki (A5040460113), Sachi Taniyasu (A5018485133), Heesoo Eun
  (A5050083483), Nobuyoshi Yamashita (A5023459603): no follow-up of the 2023 Tsukuba field (no in-season water, percolation,
  soil by depth, TOP). Only Eun 2022 (lysimeter), Yamazaki 2024 Bunseki Kagaku, Eun 2026 JAGH, Yamashita 2024 (all known).
- Web/J-STAGE (Japanese): 有機フッ素化合物 水田 黒ボク土 浸透 灌漑水 イネ; PFAS 水田 玄米 もみ殻 土壌溶液 — nothing beyond
  MAFF and known items; JAGH 68:35 (Nojima, analytics review) and JEC 36:s34 (Nomiyama, monitoring review) not relevant.
- TOP / precursor assay in paddy or flooded rice soil: none found (agricultural-soil dTOP work exists for upland soils,
  Rastatt-type; not paddy).
- PFAS sorption/desorption kinetics on allophanic/Andisol soils: none found (only herbicide kinetics on volcanic soils).
- Paddy pore-water PFAS time series through a season: none with per-congener values in a readable full text
  (Eun 2022 figure-only; Zhang 2025 / Wu 2024-2026 closed).

## (B) Grain loading route and timing

### Checked
- Tang et al. 2026 JHM (10.1016/j.jhazmat.2025.141017; PFOA/PFOS/GenX "across the full growth cycle", Zhejiang) is the
  paper behind the Chinese news items; already in the repo (Tang 2026). News page gives no stage-resolved panicle organs.
- Li 2025 JHM (10.1016/j.jhazmat.2025.138256): already in the repo; "grains" by stage from heading, husk undefined.
- MAFF FY2024 (above): brown rice only, at harvest; PFOS < LOQ in brown rice where soil and surface water PFOS were highest
  (B-a soil 4,514-5,070 ng/kg, water 65-70 ng/L) — consistent with PFOS being excluded from brown rice (supports a high
  husk/brown or a low brown-rice PFOS, but husk not measured).

### Abstract-only
- Chen X. et al. 2024 ES&T 10.1021/acs.est.4c03525 (wheat, full life cycle near a fluorochemical park): "PFAS concentrations
  in wheat decreased in the stem and ear but increased in the leaves as wheat matured"; PFOA foliar uptake ~25 % (ABSTRACT).
  The only PFAS cereal-ear time course found; wheat, not rice. Closed (ACS).
- Liu S., Zhou J., Zhu L. 2022 JoVE 10.3791/64400 (wheat split-root, hydroponic 7 d): PFOA and PFOS move shoot -> root via
  phloem and are released to the untreated root's solution (ABSTRACT). Closed.

### Mechanism context only (not PFAS values; not read this session)
- Carey et al. 2011 New Phytol 10.1111/j.1469-8137.2011.03789.x: arsenic species phloem transport from flag leaf to grain
  during filling (flag-leaf feeding) — title/abstract level.
- Repo already holds Luo 2022 (Cd by panicle organ), Zee 1972, Oparka & Gates 1981.

### Null searches (B)
- OpenAlex title/abstract: (PFAS|perfluoro*) AND rice AND (husk|hull|panicle|glume|bran); ... AND (grain filling|growth
  stage|milk); perfluoroalkyl phloem; PFAS phloem; perfluoroalkyl split-root; perfluoroalkyl wheat ear; perfluorinated rice
  grain translocation (0 hits).
- Web: PFAS rice husk brown rice grain filling time course; perfluoroalkyl rice foliar phloem grain/panicle/husk; PFOS PFOA
  rice growth stages hydroponic heading filling husk; Chinese 全氟化合物 水稻 稻壳 糙米 灌浆期. No rice study samples husk /
  brown / rachis at several times, none uses tracer, girdling, flag-leaf removal or detached panicle.
- PFOS husk/brown discrepancy (Yamazaki 2019 >= 38-800 vs 2023 field 1.4; Choi 2018 37-47): no new rice husk/brown pair found.
