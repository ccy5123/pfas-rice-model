"""Tests for the SMILES -> parameters structure adapter (src/pfas_structure.py).

Skips cleanly when RDKit is not installed (it is an optional dependency:
`pip install -r requirements-structure.txt`)."""
import os
import sys

import pytest

pytest.importorskip("rdkit", reason="pfas_structure requires RDKit")

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

import model_api as api                         # noqa: E402
import pfas_structure as ps                     # noqa: E402

# expected (n_perfluoroC, head_group) for the calibrated congeners + GenX
EXPECT = {
    "PFBA": (3, "carboxylate"), "PFPeA": (4, "carboxylate"), "PFHxA": (5, "carboxylate"),
    "PFHpA": (6, "carboxylate"), "PFOA": (7, "carboxylate"), "PFNA": (8, "carboxylate"),
    "PFDA": (9, "carboxylate"), "PFUnDA": (10, "carboxylate"), "PFDoDA": (11, "carboxylate"),
    "PFBS": (4, "sulfonate"), "PFHxS": (6, "sulfonate"), "PFOS": (8, "sulfonate"),
    "GenX": (5, "carboxylate"),
}


@pytest.mark.parametrize("name", list(EXPECT))
def test_descriptor_recovery(name):
    """Parser recovers n_perfluoroC + head group + read-across name for every congener."""
    d = ps.descriptors(ps.KNOWN_SMILES[name])
    npf, hg = EXPECT[name]
    assert d.n_perfluoroC == npf, f"{name}: nPFC {d.n_perfluoroC} != {npf}"
    assert d.head_group == hg
    assert d.matched_name == name


def test_linear_vs_nonlinear_flags():
    assert ps.descriptors(ps.KNOWN_SMILES["PFOA"]).is_linear is True
    assert ps.descriptors(ps.KNOWN_SMILES["PFOS"]).is_linear is True
    gen = ps.descriptors(ps.KNOWN_SMILES["GenX"])
    assert gen.is_linear is False and gen.n_ether_O == 1 and gen.transport_class == "ether"


def test_known_compound_uses_curated_params():
    """A known structure reads across the CURATED parameters.json values exactly."""
    for name in ("PFOA", "PFOS", "GenX"):
        rec = api._CONG[name]
        cmpd, d = ps.compound_from_smiles(ps.KNOWN_SMILES[name])
        assert cmpd.K_PL == pytest.approx(rec["K_PL_Lkg"], rel=1e-6)
        assert cmpd.K_prot == pytest.approx(rec["K_prot_Lkg"], rel=1e-6)
        assert cmpd.f_xy == pytest.approx(rec["f_xy_recommended"], rel=1e-6)
        assert d.matched_name == name


def test_alternate_smiles_writing_still_matches():
    """A differently-written but equivalent SMILES canonicalises to the same congener."""
    alt = "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O"   # PFOA reversed
    d = ps.descriptors(alt)
    assert d.matched_name == "PFOA" and d.n_perfluoroC == 7


def test_novel_ether_not_misread_as_genx():
    """A different ether-PFCA (two backbone O) must NOT read across to GenX."""
    d = ps.descriptors("OC(=O)C(F)(F)OC(F)(F)C(F)(F)OC(F)(F)C(F)(F)F")
    assert d.matched_name is None
    assert d.n_ether_O == 2 and d.is_linear is False


def test_sulfonamide_speciation_warning():
    """Sulfonamide (not a permanent anion) is detected and flagged."""
    d = ps.descriptors("NS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F")
    assert d.head_group == "sulfonamide"
    assert any("PERMANENT ANION" in n for n in d.notes)


def test_invalid_smiles_raises():
    with pytest.raises(ValueError):
        ps.descriptors("not_a_smiles((")


def test_fxy_grounded_in_parameters_json():
    """Novel f_xy is interpolated from the curated PFCA monotone series + head-group offset."""
    # a novel C7 (nPFC=6) carboxylate -> between PFHxA(0.216) and PFHpA(0.098) anchors
    cmpd, d = ps.compound_from_smiles("OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F")
    assert d.matched_name == "PFHpA"          # actually a known one -> curated
    # genuinely novel: C13 PFCA (nPFC=12) below the PFDoDA anchor 0.003
    cmpd2, d2 = ps.compound_from_smiles("OC(=O)" + "C(F)(F)" * 11 + "C(F)(F)F")
    assert d2.matched_name is None and d2.n_perfluoroC == 12
    assert cmpd2.f_xy <= 0.0031               # monotone decline, clamped at long-chain anchor


def test_simulate_from_smiles_matches_canonical():
    """simulate_from_smiles on a known structure == simulate on the named congener."""
    for name in ("PFOA", "PFOS"):
        a = api.simulate(name, season=150.0)
        b = api.simulate_from_smiles(ps.KNOWN_SMILES[name], season=150.0)
        assert b["provisional"] is False
        for k in a["baf_final"]:
            assert b["baf_final"][k] == pytest.approx(a["baf_final"][k], rel=1e-6)


def test_simulate_from_smiles_runs_novel():
    """A novel PFAS (not in the curated set) runs end-to-end from SMILES."""
    r = api.simulate_from_smiles("OC(=O)" + "C(F)(F)" * 11 + "C(F)(F)F", season=150.0)  # C13
    assert r["success"] and r["provisional"] is True
    assert r["baf_final"]["root"] > 0 and r["descriptors"].n_perfluoroC == 12
