# Literature parameter database

Curated empirical parameters and calibration/validation targets for the
PFAS–rice compartmental uptake model, organised by the six data categories of
the model (C1–C6) plus an overall source shortlist and a gap analysis.

## Files
- `PFAS_rice_parameter_database.xlsx` — source workbook (one sheet per tab).
- `*.csv` — per-sheet exports (diff-able, machine-readable; regenerate from the
  xlsx with `openpyxl`).
- `raw_si/` — per-congener numeric tables extracted from the cited papers'
  **Supporting Information** (Chen 2025 K_MW/HSA K_D, Droge 2019 K_MW, Zhou 2025
  BSA K_A, Kim 2019 field concentrations + grain BAF). These populate the
  *measured* values in `src/literature_params.py`. See `raw_si/README.md`.

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
`src/literature_params.py` transcribes the **verified** quantities — QSPR slopes,
the measured Koc anchors, `f_d`, `E_m`, the **measured per-congener `K_PL`/`K_prot`**
(from `raw_si/`), and the **Kim 2019 grain BAF** calibration data — into builders
that plug into `pfas_rice_plant_module` and `soil_paddy`. Run
`python src/literature_params.py` for a demo of the QSPRs, an end-to-end
literature-parametrised run, and a Tier-1 fit of `L_Ph` to the Kim 2019 PFOA grain BAF.

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
