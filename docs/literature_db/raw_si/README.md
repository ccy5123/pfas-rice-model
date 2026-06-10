# Raw SI-derived data (provenance)

Numeric tables extracted from the **Supporting Information** of the cited papers,
used to populate `src/literature_params.py`. Only small derived parameter tables
are kept here — the publisher PDFs/DOCX themselves are **not** committed
(copyright); regenerate these CSVs from the SI files if needed.

| file | source | DOI | what it holds |
|------|--------|-----|---------------|
| `chen2025_kmw_hsa.csv` | Chen et al. 2025, ES&T (Table S5) | `10.1021/acs.est.4c06734` | per-PFAS **log K_MW** (membrane–water, L/kg lipid, pH 7.0) and **K_D** (HSA, µmol/L, pH 7.4), 60 PFAS |
| `droge2019_kmw.csv` | Droge 2019, ES&T (SI Tables) | `10.1021/acs.est.8b05052` | **log K_MW,SSLM** (measured) and COSMOmic log K_DMPC–W — cross-check of Chen's K_MW |
| `zhou2025_kprow.csv` | Zhou et al. 2025, Ecotox. Environ. Saf. (**main text Table 1**) | `10.1016/j.ecoenv.2025.117902` | dialysis **log K_prow** [L/kg] for 4 proteins (BSA, FAF-BSA, C-phycocyanin, **soy protein isolate**), 7 PFAS — the protein–water partition the model uses (`K_prot`) |
| `zhou2025_bsa_ka.csv` | Zhou et al. 2025 (SI Table S4) | `10.1016/j.ecoenv.2025.117902` | **BSA** binding constant K_A [L/mol] and n (300 K), 6 PFAS — fluorescence quenching (a binding constant, not a partition coefficient) |
| `kim2019_field_conc.csv` | Kim et al. 2019, STOTEN (Table 4) | `10.1016/j.scitotenv.2019.03.240` | per-congener field averages: pore (void) water [ng/L], paddy soil [ng/g dw], brown rice [ng/g], rice DF [%] |
| `kim2019_grain_baf.csv` | derived from Kim 2019 | — | brown-rice (grain) **BAF** [L/kg], pore-water basis — `calibration.load_baf_csv`-ready |

## Notes / caveats
- **K_MW units**: Chen and Droge both use the SSLM/TRANSIL method; their PFCA
  values agree within ~0.2 log → confirmed **L/kg lipid**. `K_PL = 10**logK_MW`.
- **K_prot**: the model uses the **dialysis-measured `K_prow`** (Zhou 2025
  **Table 1**, `zhou2025_kprow.csv`) — `protein="soy"` (soy protein isolate, the
  grain storage-protein analog) for plant tissues, `"bsa"` for an animal
  reference. The Chen HSA `K_D` and Zhou BSA `K_A` are **binding constants**, not
  partition coefficients; the single-site `K_D` → partition route OVERESTIMATES
  by ~50× (BSA PFOA: dialysis 110 vs `K_D`-derived ~5850), so they are kept for
  reference only (`k_prot_albumin()`).
- **Kim 2019 BAFs** are **field ensemble averages** (pore water n=27; soil & rice
  n=30), **not** same-site paired, and brown-rice averages for low detection-
  frequency congeners are uncertain — treat as approximate calibration anchors.
  Only congeners with detectable brown rice are included.
