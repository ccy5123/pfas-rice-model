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
    from rdkit.Chem import Crippen as _rdCrippen
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")          # silence parse warnings (we handle None)
except ImportError as exc:                  # pragma: no cover
    raise ImportError(
        "pfas_structure requires RDKit.  Install it with:\n"
        "    pip install -r requirements-structure.txt   (or: pip install rdkit)"
    ) from exc

import literature_params as L
from pfas_rice_plant_module import Compound

# Minimum perfluorinated carbons for the PFAS (permanent-anion) branch.  Below this a
# carboxylate/sulfonate is an ordinary ionizable organic, not a PFAS -- see
# Descriptors.compound_class.  3 keeps PFBA/PFPrA (the shortest curated congeners) on
# the PFAS branch while sending 2,4-D, benzoic acid etc. to the neutral branch.
PFAS_MIN_PERFLUORO_C = 3

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
    logKow_crippen: float | None = None      # RDKit Crippen MolLogP (neutral-branch QSPR input)
    notes: list[str] = field(default_factory=list)

    @property
    def transport_class(self) -> str:
        """Head-group class used for the f_xy offset.  Ether backbone -> 'ether'."""
        return "ether" if self.n_ether_O > 0 else self.head_group

    @property
    def compound_class(self) -> str:
        """``'PFAS'`` or ``'organic'`` -- which model branch the structure belongs to.

        PFAS means a **permanently dissociated** anion: a strong acid head group on a
        perfluorinated backbone (the CF2 chain is what drives the pKa to ~0).  A
        carboxylic acid with NO perfluorination is an ordinary weak acid (2,4-D has
        head_group 'carboxylate' and n_perfluoroC 0), so both conditions are required.
        Anything else routes to the neutral / weak-electrolyte DPU branch, where the
        caller may still supply a ``pKa``.
        """
        return ("PFAS" if (self.n_perfluoroC >= PFAS_MIN_PERFLUORO_C
                           and self.head_group in ("carboxylate", "sulfonate"))
                else "organic")


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
        logKow_crippen=float(_rdCrippen.MolLogP(mol)),
    )
    d.matched_name = _match_known(d)
    if d.compound_class == "organic":
        # NOT a permanent anion -> the neutral / weak-electrolyte DPU branch handles it
        # (phase 3).  This used to be a bare "assumption violated" flag with nowhere to go.
        d.notes.append(
            f"compound class: ORGANIC (head '{head}', {d.n_perfluoroC} perfluorinated C) -- "
            "not a permanent anion, so the PFAS branch does not apply.  Use "
            "neutral_compound_from_smiles() / simulate_from_smiles(): binding comes from the "
            f"Briggs lipid term at Crippen logKow={d.logKow_crippen:.2f} (an ESTIMATE -- supply "
            "a measured logKow when you have one), and a pKa makes it a weak electrolyte.")
    elif head not in ("carboxylate", "sulfonate"):
        d.notes.append(f"head group '{head}': the model assumes a PERMANENT ANION (f_d~1); "
                       "sulfonamides violate this -> speciation is APPROXIMATE")
    if not d.is_linear and d.compound_class == "PFAS":
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
# 3b. structure -> NEUTRAL / WEAK-ELECTROLYTE Compound   (phase 3)
# ---------------------------------------------------------------------------
def neutral_compound_from_smiles(smiles: str, *, name: str | None = None,
                                 logKow: float | None = None,
                                 pKa: float | None = None, is_acid: bool = True,
                                 pH: float = L.PADDY_PH, P_n: float | None = None,
                                 L_Ph: float = 1.0, f_xy: float | None = None,
                                 K_AW: float = 0.0):
    """Build a NEUTRAL / WEAK-ELECTROLYTE :class:`Compound` from a SMILES string.

    Returns ``(compound, descriptors)`` -- the same contract as
    :func:`compound_from_smiles`, for the other half of the compound spectrum.

    This branch is in one way *easier* than the PFAS one: the whole neutral DPU
    parameterisation keys off a single descriptor, ``log K_ow``, and RDKit's Crippen
    ``MolLogP`` estimates it directly from the structure.  There is no PFAS analogue
    of that -- a PFAS ``K_ow`` is not even well defined (the molecule is both
    hydrophobic and hydrophilic and is permanently ionised), which is why the PFAS
    branch had to build read-across plus a fragment QSPR instead.

    ``logKow`` : pass a MEASURED value to override Crippen.  Recommended when one
        exists: Crippen is an atom-contribution estimate (typically within ~0.5-1 log
        unit, worse for polar heterocycles), and every downstream neutral parameter --
        ``K_lip`` (Briggs), TSCF, ``K_oc`` -- is a function of it, so its error
        propagates everywhere.  The descriptor notes record which source was used.
    ``pKa`` : ``None`` (default) = strictly neutral.  Supplying it makes the compound a
        weak electrolyte, with ``(f_n, f_d)`` from :func:`literature_params.speciation`
        at the SOIL ``pH`` and both membrane pathways running in parallel.  pKa is NOT
        predicted from structure here -- RDKit has no reliable pKa model, so an
        unsupplied pKa means "treat as neutral", which is stated rather than guessed.
    """
    d = descriptors(smiles)
    if logKow is None:
        logKow = float(d.logKow_crippen)
        kow_src = f"RDKit Crippen MolLogP = {logKow:.2f} (ESTIMATE)"
    else:
        logKow = float(logKow)
        kow_src = f"user-supplied measured logKow = {logKow:.2f}"

    cmpd = L.neutral_compound(
        logKow, name=name or d.formula, pKa=pKa, is_acid=is_acid, pH=pH,
        P_n=L.PN_DEFAULT if P_n is None else float(P_n),
        L_Ph=L_Ph, f_xy=f_xy, K_AW=K_AW)
    cmpd.logKow = logKow
    cmpd.mol_weight = float(d.mol_weight)

    spec = ("strictly NEUTRAL (pKa not supplied -> f_n=1, f_d=0)" if pKa is None else
            f"WEAK {'ACID' if is_acid else 'BASE'} (pKa={pKa:g} at pH={pH:g} -> "
            f"f_n={cmpd.fn:.3f}, f_d={cmpd.fd:.3f})")
    d.notes.append(f"NEUTRAL branch: logKow from {kow_src};  {spec};  "
                   f"K_lip=Briggs a*Kow^b;  f_xy=Briggs TSCF bell"
                   + ("" if f_xy is None else " (user-supplied)"))
    if d.compound_class == "PFAS":
        d.notes.append("WARNING: this structure looks like a PFAS (permanent anion) but was "
                       "run on the NEUTRAL branch -- the Briggs Kow correlations are not "
                       "calibrated for perfluorinated anions.  Use compound_from_smiles().")
    return cmpd, d


def compound_from_smiles_auto(smiles: str, **kw):
    """Dispatch a SMILES to the right branch by :attr:`Descriptors.compound_class`.

    Returns ``(compound, descriptors, compound_class)``.  ``compound_class`` is
    ``'PFAS'`` or ``'organic'``; PFAS goes to :func:`compound_from_smiles`, anything
    else to :func:`neutral_compound_from_smiles`.  Keyword arguments are forwarded to
    whichever builder is chosen, so pass only that branch's arguments (e.g. ``pKa=``
    only makes sense on the neutral branch).

    Supplying ``pKa`` forces the neutral/weak-electrolyte branch even for a
    PFAS-looking structure -- that is the honest reading of an explicit pKa, and it is
    how a user models e.g. a fluorotelomer or sulfonamide that is NOT fully dissociated.
    """
    d = descriptors(smiles)
    cls = "organic" if kw.get("pKa") is not None else d.compound_class
    if cls == "PFAS":
        c, d = compound_from_smiles(smiles, **kw)
    else:
        c, d = neutral_compound_from_smiles(smiles, **kw)
    return c, d, cls


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
