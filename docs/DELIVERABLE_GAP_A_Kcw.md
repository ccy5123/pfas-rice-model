# GAP A — Cell-wall–water partition coefficient `K_cw` (DELIVERABLE)

**Status: CLOSED.** Recommended values in `params/Kcw_Klignin_params_v2.csv` (14 congeners ×
intrinsic poly/lignin + whole-cw resolved to 5 organs). Consumed by `B_k` as the `f_cw·K_cw` term.

## Recommended values

`K_cw` enters the 4-pool `B_k` as **whole-cell-wall per organ**, `K_cw_wholecw_<organ>` =
`(f_poly·K_poly + f_lignin·K_lignin)/(f_poly+f_lignin)` [L/kg pool-dw]. Intrinsic anchors
(rec, PFOA): `K_cw_poly = 4.0`, `K_lignin = 17`. Whole-cw (PFOA): root 6.6, stem 6.53, leaf 5.86,
husk 7.69, grain_brown 7.71. Chain trend: log-linear, `+0.10/CF₂`; sulfonate head-group offset
`+0.40 log` over the carboxylate of equal chain. Full table per congener in the CSV; consolidated
into `params/parameters.json` (`K_cw_wholecw_Lkg`, `K_cw_poly_rec_Lkg`, `K_lignin_rec_Lkg`).

## Method / verdict

No rice-direct batch-sorption `K_cw` exists, so values are **anchored, not back-calculated from
in-planta data**:
- **Lignin** — direct anchor: Mel et al. 2024 (JHM) measured PFAS–lignin Kd at pH 6.5, I = 20 mM
  (`K_lignin` rec; Mel-interp C4–C8, Mel-extrap 0.10/CF₂ for >C8).
- **Polysaccharide** — component + DFT: Guo et al. 2025 (Nat. Commun., Fig. 3f) M06-2X/def2-TZVPP
  DFT association free energies give a **relative** site ladder (Glu > GalA > Xyl > Gal; mechanism =
  COO⁻···sugar-OH H-bond, perfluoro tail non-binding). Absolute scale from measured component Kd,
  **not** literal ΔG→Kd (rejected: implicit-solvent association ΔG under-predicts macroscopic
  partitioning by 2–4 orders).
- **Rice organ weighting** — rice_tissue package poly/lignin dw fractions (root 0.40/0.10,
  stem 0.58/0.14, leaf 0.48/0.08, husk 0.48/0.19, grain 0.025/0.010).

**Rejections (kept out of `K_cw`):** literal ΔG→Kd conversion; soil-OM per-CF₂ slopes (0.60/0.83 —
incompatible with the non-binding-tail mechanism, gives `K_lignin`≈1500 vs measured 35 for PFOS);
soil-OC Koc extrapolation. Surface sorption / Fe-Mn plaque / suberin / silica are a **separate
sorption term outside `K_cw`** (double-counting otherwise) — in the model as `Compound.K_surf`
(default 0; see GAP B / surface note).

## Experimental design (to replace the anchor)

Batch-sorption isotherm on **isolated rice-root cell-wall fractions** (cellulose+hemicellulose+
pectin pool, silica-removed) at paddy porewater conditions (pH 5.5–7.0, I ≈ 1–20 mM), per congener
C4–C14, to yield a directly-measured `K_cw_poly`. Secondary: rice cw monosaccharide composition
(Guo Table 11 rice panel) to re-anchor the grass-specific polysaccharide intrinsic.

## Limitations

- **Direct `K_cw_poly` measurement absent — the long-term weakest point.** Current value is
  DFT-mechanism + dicot-component anchored; grass (rice) composition (pectin-poor, arabinoxylan-rich)
  could shift it slightly but Glu-dominance offsets — qualitative caveat only, kept at the dicot anchor.
- Organ whole-cw spread is ±15–20%; the 4-pool module uses a single `K_cw` per compound (root
  whole-cw) — acceptable since cell wall is a minor pool for the membrane-dominated long chains.
- PFSA only to C8 (PFOS); C9–C14 PFSA via slope 0.10 + offset 0.40 extrapolation (flag if used).

## Sources

| key | reference | role | DOI |
|---|---|---|---|
| Mel 2024 | Mel, Lau, Hockaday 2024, *JHM* 480:136016 | K_lignin direct anchor (pH 6.5) | 10.1016/j.jhazmat.2024.136016 |
| Guo 2025 | Guo et al. 2025, *Nat. Commun.* 16:10283 | DFT ΔG ladder + mechanism + monosaccharide | 10.1038/s41467-025-65191-3 |
| rice_tissue | user package | organ poly/lignin dw fractions | — |
| Liu 2023 | Liu et al. 2023, *ES&T* 57:8739 | context (surface-sorption dominance → K excluded) | 10.1021/acs.est.3c00504 |
