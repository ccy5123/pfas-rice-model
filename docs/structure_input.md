# Structure (SMILES) input — parameterising any PFAS from its chemistry

`src/pfas_structure.py` lets a PFAS **chemical structure** (a SMILES string) be the
model input, mapping it onto a `Compound` so the uptake model runs for *any* PFAS,
not only the 12 curated congeners + GenX. This is the "option 3" front end: a
structure → parameters adapter built on **RDKit**.

## What it is (and isn't)

It is **mechanistic read-across + a fragment QSPR**, *not* a black-box ML model
(there is no training set to fit one honestly). Pipeline:

```
SMILES --(RDKit)--> structural descriptors --map--> Compound
        parse + SMARTS    {n_perfluoroC, head_group,      (1) MEASURED read-across if the
                           n_ether_O, n_CF3, branched,         structure matches a curated
                           MW, formula, is_linear}             congener  -> parameters.json
                                                          (2) else QSPR (per-CF2 slope +
                                                              head-group offset; ether/
                                                              sulfonamide PROVISIONAL)
```

| Parameter | How the structure sets it |
|---|---|
| `K_PL`, `K_prot`, `K_cw` (binding) | curated value (known congener) **or** chain-length QSPR / measured anchor |
| `f_d` (speciation) | head-group pKa (carboxylate/sulfonate → anion, `f_d≈1`) |
| `f_xy` (translocation) | **not structure-derivable** — curated monotone `f_xy_recommended` (known) or the PFCA monotone series interpolated on `n_perfluoroC` × the head-group offset (novel, PROVISIONAL) |
| `L_Ph`, `kappa_d`, carrier | fitted defaults (not from structure) |

The structure fixes binding + speciation + the head-group **ordering** of `f_xy`;
the **absolute** translocation scale stays fitted (same caveat as the whole model).

## Usage

```python
import model_api as api
from pfas_structure import compound_from_smiles, descriptors

# 1. structure -> parameters (+ a descriptor record with provenance notes)
cmpd, d = compound_from_smiles("OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F")
print(d.n_perfluoroC, d.head_group, d.matched_name)     # 7 carboxylate PFOA
print(cmpd.K_PL, cmpd.f_xy)                              # 1905.0 0.04   (curated)

# 2. structure -> full simulation (4-compartment ODE)
res = api.simulate_from_smiles("OC(=O)" + "C(F)(F)"*11 + "C(F)(F)F", season=150.0)  # PFTrDA C13
print(res["provisional"], res["baf_final"])             # True {...}  (novel, runs from SMILES)
```

`compound_from_smiles` returns `(Compound, Descriptors)`; `Descriptors.notes`
records what was measured vs predicted vs assumed. `simulate_from_smiles` returns
the usual `simulate(...)` dict plus `descriptors` and `provisional`.

## Read-across vs novel (the honest boundary)

- **Known congener** (canonical-SMILES or, for linear chains, `(n_perfluoroC,
  head_group)` match): uses the **curated `params/parameters.json`** values — a
  SMILES-built PFOA reproduces the named PFOA exactly (`tests/test_pfas_structure.py`).
- **Novel structure**: binding from the QSPR (measured `K_PL` anchor where it exists);
  `f_xy` from the monotone series × head-group offset; **flagged `provisional`**.
  Non-linear (ether/branched) structures are explicitly outside the carboxylate/
  sulfonate QSPR calibration domain — the GenX BCF over-prediction is the cautionary
  example (an ether-specific Koc/QSPR is the documented follow-up).
- **Sulfonamides / neutral species**: detected and flagged — the model assumes a
  *permanent anion* (`f_d≈1`), which they violate, so their speciation is approximate.

## Scope of the parser

The descriptors cover the PFAS subset: perfluoroalkyl/-ether carboxylic & sulfonic
acids and sulfonamides, linear and branched. Head groups are matched by SMARTS;
`n_perfluoroC` = carbons bearing ≥1 F; `n_ether_O` = backbone `-O-` bridging two
carbons. Validated against the 12 congeners + GenX (descriptor recovery + read-across
consistency). RDKit is an **optional** dependency (`pip install -r
requirements-structure.txt`); the tests skip when it is absent.

## Next steps

1. **Ether/sulfonamide fragment terms** — extend the QSPR (`koc`/`k_pl`) with
   per-ether-O and head-group contributions so novel PFECAs/PFSAs are predicted, not
   approximated as carboxylate (closes the GenX Koc gap).
2. **Measured anchors** — wire more measured `K_PL`/`K_prot`/`Koc` so read-across
   covers more of the structure space (currently the curated 13).
3. **Optional true QSAR/ML** — a trained SMILES→property model would replace the
   QSPR fallback, but needs a measured PFAS training set (a data-gap task).
