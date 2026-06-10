# Literature parameter database

Curated empirical parameters and calibration/validation targets for the
PFAS–rice compartmental uptake model, organised by the six data categories of
the model (C1–C6) plus an overall source shortlist and a gap analysis.

## Files
- `PFAS_rice_parameter_database.xlsx` — source workbook (one sheet per tab).
- `*.csv` — per-sheet exports (diff-able, machine-readable; regenerate from the
  xlsx with `openpyxl`).

## Tabs → model terms
| tab | tier | model term it feeds |
|-----|------|---------------------|
| `C1_Rice_Tissue_BAF`      | Tier-1 | calibration targets: `B_k`, `g_in/g_out`, `f_xy`, `Π` (tissue BAFs) |
| `C2_Growth_Water_Forcing` | Tier-0 | forcing functions `M_k(t)`, `Q_TP(t)` |
| `C3_Soil_Sorption` / `C3_Chainlength_QSPR` | — | soil sub-model `K_F`, `n`, `Koc` (`soil_paddy.py`) |
| `C4_Tissue_Binding_QSPR`  | Tier-3 | binding factor `B_k` (`K_prot`, `K_PL`, `K_cw`) |
| `C5_Membrane_Translocation` | Tier-2 | GHK `E_m`, carrier `V_max/K_m`, TSCF `f_xy`, phloem `L_Ph` |
| `C6_Physchem`             | Tier-0 | `f_d` (pKa), `P_d`, `D_aq` |

## How it is consumed
`src/literature_params.py` transcribes the **verified** quantities (QSPR slopes,
the measured Koc anchors, `f_d`, `E_m`) into builders that plug into
`pfas_rice_plant_module` and `soil_paddy`. Run `python src/literature_params.py`
for a demo of the QSPRs and an end-to-end literature-parametrised run.

## ⚠️ DOI status (no-fabrication rule)
Each row carries a `DOI_status`:
- `verified` — DOI string was observed verbatim during the literature search.
- `UNVERIFIED` — a lead that was **not** confirmed this session; **confirm the
  DOI before citing** (some C2 crop-model/FAO refs, the Liu 2017 per-CF2 slope,
  and the Wang 1994 rice `E_m` DOI are in this state).

## Known gaps (see `Gap_Analysis`)
Compartment-resolved TF under flooding, anoxic/flooded sorption isotherms,
rice-specific carrier `V_max/K_m`, a PFAS phloem→grain transporter, and the
absolute per-congener `K_prot`/`K_PL`/`K_cw` intercepts (slopes are known;
intercepts must be read from the cited SI or fitted).
