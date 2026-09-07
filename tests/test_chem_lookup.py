"""CompTox (CTX) lookup — offline tests.

Every test here injects a fake transport (`_get_fn`), so nothing touches the
network and nothing needs an API key: CI has neither. What is pinned is the part
that would silently corrupt a run if it broke — the identifier classification, the
EXPERIMENTAL-over-predicted preference, the Henry -> K_AW conversion, the
mapping onto `simulate_neutral` kwargs, and the no-key / no-network degradation.

The response SHAPES below are the ones the parser is written against; if EPA's
change, `python src/chem_lookup.py <compound> --raw` on a machine with the key
shows what actually arrives.
"""
import sys
import os

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import chem_lookup as cl  # noqa: E402


# --- fake transport --------------------------------------------------------
SEARCH = [{"dtxsid": "DTXSID3022125", "preferredName": "Carbamazepine",
           "casrn": "298-46-4", "smiles": "NC(=O)N1c2ccccc2C=Cc2ccccc21",
           "inchikey": "FFGPTBGBLSHEPO-UHFFFAOYSA-N"}]
DETAIL = {"dtxsid": "DTXSID3022125", "averageMass": 236.269,
          "smiles": "NC(=O)N1c2ccccc2C=Cc2ccccc21"}
PROPS = [
    # the same property twice, experimental and predicted, deliberately with the
    # PREDICTED row first so the ordering cannot pass by accident
    {"propertyId": "logP", "propType": "predicted", "value": 3.39, "modelName": "OPERA"},
    {"propertyId": "LogKow: Octanol-Water", "propType": "experimental", "value": 2.45,
     "source": "PHYSPROP", "unit": ""},
    {"propertyId": "Henry's Law constant", "propType": "predicted", "value": 1.08e-10,
     "unit": "atm-m3/mol", "modelName": "OPERA"},
    {"propertyId": "logKoc", "propType": "predicted", "value": 2.4, "modelName": "OPERA"},
    {"propertyId": "Melting Point", "propType": "experimental", "value": 190.0, "unit": "°C"},
]


def fake_get(payloads=None, fail=None):
    """Build a `_get_fn` returning canned payloads keyed by a path fragment."""
    payloads = payloads if payloads is not None else {
        "/chemical/search/": SEARCH,
        "/chemical/detail/": [DETAIL],
        "/chemical/property/search": PROPS,
    }

    def _get(path, key, base_url=cl.DEFAULT_BASE_URL, params=None):
        if fail:
            return None, fail
        for frag, payload in payloads.items():
            if frag in path:
                return payload, None
        return None, "not found (HTTP 404)"
    return _get


# --- identifiers -----------------------------------------------------------
def test_identifier_kind_distinguishes_name_cas_and_structure():
    """A NAME must never be read as a structure: 'carbamazepine' is a valid SMILES
    string to a permissive parser, and mis-classifying it would send the wrong
    query to the API."""
    assert cl.identifier_kind("carbamazepine") == "name"
    assert cl.identifier_kind("Perfluorooctanoic acid") == "name"
    assert cl.identifier_kind("298-46-4") == "casrn"
    assert cl.identifier_kind("DTXSID3022125") == "dtxsid"
    assert cl.identifier_kind("FFGPTBGBLSHEPO-UHFFFAOYSA-N") == "inchikey"
    rdkit = pytest.importorskip("rdkit", reason="SMILES detection needs RDKit")
    assert rdkit is not None
    assert cl.identifier_kind("NC(=O)N1c2ccccc2C=Cc2ccccc21") == "smiles"
    assert cl.identifier_kind("OC(=O)C(F)(F)C(F)(F)F") == "smiles"


def test_smiles_is_resolved_through_its_inchikey():
    """CTX has no SMILES search, so a structure query goes via the InChIKey."""
    pytest.importorskip("rdkit")
    assert cl.smiles_to_inchikey("NC(=O)N1c2ccccc2C=Cc2ccccc21").startswith("FFGPTBGBLSHEPO")
    assert cl.smiles_to_inchikey("not a molecule") == ""


# --- provenance ------------------------------------------------------------
def test_experimental_value_beats_the_predicted_one():
    """The neutral path's a-priori numbers are on MEASURED log Kow, so when both
    exist the experimental value must win — and say so."""
    r = cl.lookup("carbamazepine", key="k", _get_fn=fake_get())
    assert r.ok and r.dtxsid == "DTXSID3022125"
    p = r.props["log_kow"]
    assert p.value == pytest.approx(2.45)          # NOT the OPERA 3.39
    assert p.is_experimental and "predicted" not in p.badge()
    # a property with only a predicted row is still returned, but labelled
    assert r.props["log_koc"].source == "predicted"
    assert "OPERA" in r.props["log_koc"].badge()


def test_henry_becomes_a_dimensionless_kaw_or_nothing():
    """K_AW switches the leaf's volatilisation sink on; a wrong unit conversion is
    worse than a missing value, so an unrecognised unit returns None."""
    assert cl.henry_to_kaw(1.0, "atm-m3/mol") == pytest.approx(1.0 / 0.024465, rel=1e-3)
    assert cl.henry_to_kaw(101325.0, "Pa-m3/mol") == pytest.approx(1.0 / 0.024465, rel=1e-3)
    assert cl.henry_to_kaw(0.5, "dimensionless") == pytest.approx(0.5)
    assert cl.henry_to_kaw(1.0, "mmHg") is None
    r = cl.lookup("carbamazepine", key="k", _get_fn=fake_get())
    assert r.props["K_AW"].value == pytest.approx(1.08e-10 / 0.024465, rel=1e-3)


def test_neutral_kwargs_maps_onto_simulate_neutral_and_omits_half_life():
    """What the lookup hands the model — and what it must NOT: the in-planta
    half-life is not a dashboard property (and per Kodesova not even a compound
    constant), so it is never auto-filled."""
    r = cl.lookup("carbamazepine", key="k", _get_fn=fake_get())
    nk = r.neutral_kwargs()
    assert nk["log_kow"] == pytest.approx(2.45)
    assert nk["MW"] == pytest.approx(236.269, rel=1e-4)
    assert nk["Koc"] == pytest.approx(10 ** 2.4, rel=1e-6)
    assert "half_life" not in nk
    assert "pKa" not in nk                        # none in the fixture -> strictly neutral
    # and the kwargs are accepted by the real entry point
    api = pytest.importorskip("model_api")
    res = api.simulate_neutral(nk["log_kow"], MW=nk["MW"], K_AW=nk["K_AW"],
                               air=False, half_life=7.0, season=60.0, n_t=61)
    assert res["baf_final"]["root"] > 0


def test_pka_is_only_filled_when_one_was_found():
    """`pKa=None` is the strictly-neutral, bit-identical path — a missing pKa must
    stay missing rather than defaulting to some number."""
    props = list(PROPS) + [{"propertyId": "pKa_a", "propType": "predicted",
                            "value": 13.9, "modelName": "OPERA"}]
    r = cl.lookup("carbamazepine", key="k",
                  _get_fn=fake_get({"/chemical/search/": SEARCH,
                                    "/chemical/detail/": [DETAIL],
                                    "/chemical/property/search": props}))
    nk = r.neutral_kwargs()
    assert nk["pKa"] == pytest.approx(13.9) and nk["is_acid"] is True


# --- degradation -----------------------------------------------------------
def test_no_key_and_no_network_degrade_to_a_note_not_an_exception():
    """The app must stay usable with manual entry when the lookup cannot run."""
    r = cl.resolve("carbamazepine", key="", _get_fn=fake_get())
    assert not r.ok and "key" in r.note.lower()

    r = cl.lookup("carbamazepine", key="k", _get_fn=fake_get(fail="ConnectionError: blocked"))
    assert not r.ok and "blocked" in r.note

    r = cl.lookup("not-a-real-compound-xyz", key="k",
                  _get_fn=fake_get({"/chemical/detail/": [DETAIL]}))
    assert not r.ok and r.note


def test_response_envelope_and_vocabulary_changes_degrade_safely():
    """A wrapped payload is still read; a property this module does not know about
    yields a MISSING property, never a wrong one attached to the right slot.

    Matching is by name SUBSTRING on purpose — the CTX property vocabulary is not
    stable enough to hard-code ids — so the guarantee is one-directional: an
    unrelated name is ignored, but a renamed logP that still contains "logp" WILL
    be picked up (which is the intent, not a defect)."""
    wrapped = {"data": SEARCH}
    renamed = [{"propertyId": "Bioconcentration factor", "propType": "experimental", "value": 9.9}]
    r = cl.lookup("carbamazepine", key="k",
                  _get_fn=fake_get({"/chemical/search/": wrapped,
                                    "/chemical/detail/": {"data": [DETAIL]},
                                    "/chemical/property/search": renamed}))
    assert r.ok and r.mol_weight == pytest.approx(236.269, rel=1e-4)
    assert "log_kow" not in r.props and r.neutral_kwargs().get("log_kow") is None


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
