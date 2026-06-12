"""
Structure (SMILES) -> model parameters adapter for PFAS  [RDKit]
================================================================

The **option-3** front end: let a PFAS *chemical structure* (a SMILES string)
be the input and map it onto a :class:`pfas_rice_plant_module.Compound`, so the
uptake model can be run for *any* PFAS, not only the hand-curated congeners.

Honesty / scope
---------------
This is **mechanistic read-across + a fragment QSPR**, NOT a black-box ML model
(there is no training set to fit one honestly).  The pipeline is::

    SMILES --(RDKit)--> structural descriptors {n_perfluoroC, head_group,
                         n_ether_O, n_CF3, branched, MW, ...}
           --map-->     Compound parameters, by:
                          (1) MEASURED read-across -- if the (canonical) structure
                              matches a congener with a measured K_PL/K_prot/Koc
                              (Chen2025, Zhou2025, Milinovic2015), use the lab value;
                          (2) else the literature_params QSPR -- per-CF2 slope +
                              head-group offset (carboxylate/sulfonate calibrated;
                              ether/sulfonamide PROVISIONAL).

The translocation parameters (``f_xy``, ``L_Ph``, ``kappa_d``, carrier
``Vmax/Km``) are NOT structure-derivable: the structure fixes only the head-group
*ordering*/offset; absolute values stay fitted.  For a known congener ``f_xy``
is taken from the curated ``params/parameters.json`` (monotone ``f_xy_recommended``);
for a novel structure it is the carboxylate monotone series (interpolated on
n_perfluoroC) times the head-group offset, flagged PROVISIONAL.  Every Compound
carries descriptor ``notes`` recording what was measured vs predicted vs assumed.

Requires RDKit (``pip install -r requirements-structure.txt``).  Validated against
the 12 calibrated congeners + GenX in ``tests/test_pfas_structure.py``.
"""
from __future__ import annotations

import json
import math
import os
from dataclasses import dataclass, field

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors as _rdDesc
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")          # silence parse warnings (we handle None)
except ImportError as exc:                  # pragma: no cover
    raise ImportError(
        "pfas_structure requires RDKit.  Install it with:\n"
        "    pip install -r requirements-structure.txt   (or: pip install rdkit)"
    ) from exc

import literature_params as L
from pfas_rice_plant_module import Compound

_PARAMS_JSON = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                            "params", "parameters.json")

# ---------------------------------------------------------------------------
# 0. canonical SMILES for the known congeners (for read-across + tests)
# ---------------------------------------------------------------------------
def _pfca(n_perfluoro: int) -> str:
    return "OC(=O)" + "C(F)(F)" * (n_perfluoro - 1) + "C(F)(F)F"


def _pfsa(n_perfluoro: int) -> str:
    return "OS(=O)(=O)" + "C(F)(F)" * (n_perfluoro - 1) + "C(F)(F)F"


KNOWN_SMILES: dict[str, str] = {
    "PFBA": _pfca(3), "PFPeA": _pfca(4), "PFHxA": _pfca(5), "PFHpA": _pfca(6),
    "PFOA": _pfca(7), "PFNA": _pfca(8), "PFDA": _pfca(9), "PFUnDA": _pfca(10),
    "PFDoDA": _pfca(11),
    "PFBS": _pfsa(4), "PFHxS": _pfsa(6), "PFOS": _pfsa(8),
    # GenX / HFPO-DA: C3F7-O-CF(CF3)-COOH  (an ether-PFCA, 5 perfluoro-C)
    "GenX": "OC(=O)C(F)(OC(F)(F)C(F)(F)C(F)(F)F)C(F)(F)F",
}

# GenX measured K_PL lives in build_parameters.py (not the literature_params KMW
# dict); mirror it so read-across covers the 13th congener.
_KMW_LOG_EXTRA = {"GenX": math.log10(117.5)}     # Chen 2025 HFPO-DA, log K_MW = 2.07

# SMARTS for the (anionic) PFAS head groups, tried in priority order.
_HEAD_SMARTS = [
    ("carboxylate", "[CX3](=[OX1])[OX2H1,OX1-]"),
    ("sulfonate",   "[SX4](=[OX1])(=[OX1])[OX2H1,OX1-]"),
    ("sulfonamide", "[SX4](=[OX1])(=[OX1])[NX3]"),
    ("phosphonate", "[PX4](=[OX1])([OX2H1,OX1-])[OX2H1,OX1-]"),
]
_HEAD_PATTERNS = [(nm, Chem.MolFromSmarts(sm)) for nm, sm in _HEAD_SMARTS]


def _canon(smiles: str) -> str | None:
    m = Chem.MolFromSmiles(smiles)
    return None if m is None else Chem.MolToSmiles(m)


_CANON_KNOWN = {c: nm for nm, smi in KNOWN_SMILES.items() if (c := _canon(smi))}


# ---------------------------------------------------------------------------
# 1. f_xy monotone series from params/parameters.json (grounding, not a fit here)
# ---------------------------------------------------------------------------
def _load_congener_records():
    """Curated congener records (by name) + the PFCA f_xy(nPFC) monotone series,
    both from params/parameters.json -- the read-across / grounding source."""
    with open(_PARAMS_JSON) as f:
        cong = json.load(f)["congeners"]
    by_name, pfca_pts = {}, []
    for c in cong:
        by_name[c["name"]] = c
        if c["group"] == "PFCA":
            pfca_pts.append((c["n_C"] - 1, float(c["f_xy_recommended"])))   # PFCA: nPFC = C-1
    pfca_pts.sort()
    return by_name, pfca_pts


_CONG_REC, _FXY_PFCA = _load_congener_records()


def _fxy_carboxylate(n_perfluoroC: float) -> float:
    """Monotone carboxylate f_xy(nPFC), log-interpolated from parameters.json PFCA."""
    xs = [p[0] for p in _FXY_PFCA]
    ys = [math.log10(p[1]) for p in _FXY_PFCA]
    if n_perfluoroC <= xs[0]:
        return float(10.0 ** ys[0])
    if n_perfluoroC >= xs[-1]:
        return float(10.0 ** ys[-1])
    for k in range(1, len(xs)):
        if n_perfluoroC <= xs[k]:
            t = (n_perfluoroC - xs[k - 1]) / (xs[k] - xs[k - 1])
            return float(10.0 ** (ys[k - 1] + t * (ys[k] - ys[k - 1])))
    return float(10.0 ** ys[-1])


# ---------------------------------------------------------------------------
# 2. structural descriptors (RDKit)
# ---------------------------------------------------------------------------
@dataclass
class Descriptors:
    """PFAS-relevant structural descriptors extracted from a SMILES via RDKit."""
    smiles: str
    canonical_smiles: str
    formula: str
    mol_weight: float
    n_C: int
    n_F: int
    n_perfluoroC: int            # carbons bearing >= 1 F (the QSPR chain length)
    n_CF2: int
    n_CF3: int
    n_ether_O: int               # backbone -O- bridging two carbons
    head_group: str              # carboxylate / sulfonate / sulfonamide / phosphonate / unknown
    branched: bool               # a carbon with >= 3 carbon neighbours
    is_linear: bool              # straight perfluoroalkyl acid (QSPR-calibrated domain)
    matched_name: str | None = None
    notes: list[str] = field(default_factory=list)

    @property
    def transport_class(self) -> str:
        """Head-group class used for the f_xy offset.  Ether backbone -> 'ether'."""
        return "ether" if self.n_ether_O > 0 else self.head_group


def descriptors(smiles: str) -> Descriptors:
    """Parse a SMILES into :class:`Descriptors` using RDKit."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles!r}")
    canon = Chem.MolToSmiles(mol)

    cf = {}                                   # C atom idx -> bonded F count
    for a in mol.GetAtoms():
        if a.GetSymbol() == "C":
            cf[a.GetIdx()] = sum(1 for nb in a.GetNeighbors() if nb.GetSymbol() == "F")
    n_C = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "C")
    n_F = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "F")
    n_perfluoroC = sum(1 for v in cf.values() if v >= 1)
    n_CF2 = sum(1 for v in cf.values() if v == 2)
    n_CF3 = sum(1 for v in cf.values() if v == 3)

    # ether O: oxygen with exactly two carbon neighbours, both single bonds
    n_ether_O = 0
    for a in mol.GetAtoms():
        if a.GetSymbol() == "O" and a.GetDegree() == 2:
            nbrs = list(a.GetNeighbors())
            if all(nb.GetSymbol() == "C" for nb in nbrs) and all(
                    mol.GetBondBetweenAtoms(a.GetIdx(), nb.GetIdx()).GetBondTypeAsDouble() == 1.0
                    for nb in nbrs):
                n_ether_O += 1

    head = "unknown"
    for nm, patt in _HEAD_PATTERNS:
        if patt is not None and mol.HasSubstructMatch(patt):
            head = nm
            break

    branched = any(
        sum(1 for nb in a.GetNeighbors() if nb.GetSymbol() == "C") >= 3
        for a in mol.GetAtoms() if a.GetSymbol() == "C")
    is_linear = (head in ("carboxylate", "sulfonate") and n_ether_O == 0
                 and n_CF3 <= 1 and not branched)

    d = Descriptors(
        smiles=smiles, canonical_smiles=canon,
        formula=_rdDesc.CalcMolFormula(mol) if hasattr(_rdDesc, "CalcMolFormula")
        else Chem.rdMolDescriptors.CalcMolFormula(mol),
        mol_weight=float(_rdDesc.MolWt(mol)),
        n_C=n_C, n_F=n_F, n_perfluoroC=n_perfluoroC, n_CF2=n_CF2, n_CF3=n_CF3,
        n_ether_O=n_ether_O, head_group=head, branched=branched, is_linear=is_linear,
    )
    d.matched_name = _match_known(d)
    if head not in ("carboxylate", "sulfonate"):
        d.notes.append(f"head group '{head}': the model assumes a PERMANENT ANION (f_d~1); "
                       "sulfonamides/neutral species violate this -> speciation is APPROXIMATE")
    if not d.is_linear:
        d.notes.append("non-linear/ether/branched: outside the carboxylate/sulfonate QSPR "
                       "calibration domain -> binding/Koc are PROVISIONAL")
    return d


def _match_known(d: Descriptors) -> str | None:
    """Read-across.  Linear PFCA/PFSA: match by (n_perfluoroC, head) OR canonical
    SMILES.  Non-linear (ether/branched): EXACT canonical-SMILES match only."""
    if d.canonical_smiles in _CANON_KNOWN:
        return _CANON_KNOWN[d.canonical_smiles]
    if d.is_linear:
        for name, (_tot, npf, hg) in L.SPECIES.items():
            if npf == d.n_perfluoroC and hg == d.head_group:
                return name
    return None


# ---------------------------------------------------------------------------
# 3. structure -> Compound
# ---------------------------------------------------------------------------
def compound_from_smiles(smiles: str, *, name: str | None = None,
                         pH: float = L.PADDY_PH, plant_protein: bool = True,
                         f_xy: float | None = None, L_Ph: float = 0.005,
                         kappa_d: float = 0.5, Vmax_in: float = 20.0, Km_in: float = 5.0,
                         Vmax_out: float = 8.0, Km_out: float = 5.0,
                         kcw_anchor: float = L.KCW_ANCHOR_LKG):
    """Build a :class:`Compound` from a SMILES string.

    Returns ``(compound, descriptors)``.  Binding (K_PL, K_prot) and speciation
    (f_d) come from MEASURED read-across when the structure matches a known
    congener, else from the literature_params QSPR.  ``f_xy`` is taken from the
    curated monotone series (exact for a known congener; interpolated x head-group
    offset for a novel one -- PROVISIONAL) unless supplied explicitly.
    """
    d = descriptors(smiles)
    known = name or d.matched_name
    npf = d.n_perfluoroC
    acid_hg = d.head_group if d.head_group in L.PKA else "carboxylate"
    fd = float(L.f_d(L.PKA[acid_hg], pH))                 # head-group pKa; ether acid = carboxylate
    rec = _CONG_REC.get(known)

    if rec is not None:
        # --- known congener: read across the CURATED params (exact consistency) ---
        K_PL = float(rec["K_PL_Lkg"]); kpl_src = f"curated parameters.json ({known})"
        K_prot = float(rec["K_prot_Lkg"]); kprot_src = kpl_src
        kcw_anchor = float(rec["K_cw_wholecw_Lkg"]["root"])
        f_xy_rec = float(rec["f_xy_recommended"]); fxy_src = f"curated f_xy_recommended ({known})"
    else:
        # --- novel structure: QSPR (measured anchor where available) ---
        if known in L.KMW_CHEN2025_LOG:
            K_PL = float(10.0 ** L.KMW_CHEN2025_LOG[known]); kpl_src = f"measured (Chen2025, {known})"
        elif known in _KMW_LOG_EXTRA:
            K_PL = float(10.0 ** _KMW_LOG_EXTRA[known]); kpl_src = f"measured (Chen2025, {known})"
        else:
            binding_hg = d.head_group if d.head_group in L.KPL_PER_CF2 else "carboxylate"
            K_PL = L.k_pl(npf, binding_hg, n_ether_O=d.n_ether_O)
            kpl_src = (f"QSPR (CF2 slope, {binding_hg}"
                       + (f", {d.n_ether_O}x ether {L.KPL_ETHER_LOG_OFFSET:+.2f}log)" if d.n_ether_O else ")"))
        K_prot = L.k_prot(npf, plant=plant_protein); kprot_src = "QSPR (chain factor)"
        off_class = d.transport_class if d.transport_class in L.FXY_HEADGROUP_LN_OFFSET else "carboxylate"
        f_xy_rec = L.f_xy_headgroup(_fxy_carboxylate(npf), off_class)
        fxy_src = (f"monotone PFCA(nPFC={npf}) x exp({L.FXY_HEADGROUP_LN_OFFSET[off_class]:+.1f})"
                   f" [{off_class}] (PROVISIONAL)")

    if f_xy is None:
        f_xy = f_xy_rec
    else:
        fxy_src = "user-supplied"

    label = known or f"PF{npf}{'S' if d.head_group == 'sulfonate' else 'A'}*"
    provisional = (not d.is_linear) or (d.matched_name is None and name is None)
    cmpd = Compound(
        name=label, K_prot=K_prot, K_PL=K_PL, K_cw=kcw_anchor,
        kappa_d=kappa_d, Vmax_in=Vmax_in, Km_in=Km_in,
        Vmax_out=Vmax_out, Km_out=Km_out, L_Ph=L_Ph, f_xy=f_xy,
        fd=fd, fn=0.0,
    )
    d.notes.append(f"K_PL: {kpl_src};  K_prot: {kprot_src};  f_d({acid_hg});  f_xy: {fxy_src}")
    if provisional:
        d.notes.append("OVERALL: PROVISIONAL (novel/non-calibrated structure or no measured "
                       "anchor); binding predicted by QSPR, translocation is a head-group estimate")
    return cmpd, d


# ---------------------------------------------------------------------------
# 4. demo
# ---------------------------------------------------------------------------
def _demo():
    print("SMILES -> RDKit descriptors -> Compound  (read-across measured; else QSPR)\n")
    tests = [
        ("PFOA", None), ("PFOS", None), ("GenX", None),
        # novel ether-PFCA (ADONA-like) -- NOT GenX, must stay novel/QSPR
        (None, "OC(=O)C(F)(F)OC(F)(F)C(F)(F)OC(F)(F)C(F)(F)F"),
        # perfluorooctane sulfonamide (FOSA-like) -- speciation warning expected
        (None, "NS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"),
    ]
    for nm, smi in tests:
        smi = smi or KNOWN_SMILES[nm]
        c, d = compound_from_smiles(smi, name=nm)
        print(f"# {nm or '(novel)'}\n  {d.canonical_smiles}  [{d.formula}, MW {d.mol_weight:.1f}]")
        print(f"  descriptors: nPFC={d.n_perfluoroC} head={d.head_group} ether_O={d.n_ether_O} "
              f"CF3={d.n_CF3} branched={d.branched} linear={d.is_linear} match={d.matched_name}")
        print(f"  Compound: K_PL={c.K_PL:.0f}  K_prot={c.K_prot:.0f}  f_d={c.fd:.3f}  f_xy={c.f_xy:.4f}")
        for ln in d.notes:
            print(f"    - {ln}")
        print()


if __name__ == "__main__":
    _demo()
