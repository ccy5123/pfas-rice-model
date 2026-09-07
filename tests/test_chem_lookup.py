"""CompTox (CTX) lookup — offline tests, written against the REAL API's shape.

Every test injects a fake transport (`_get_fn`), so nothing touches the network and
nothing needs an API key: CI has neither. The fixtures reproduce the response shapes
and the five silent-failure traps recorded in the working notes (see the header of
`src/chem_lookup.py`), because those are exactly the failures that produce empty or
WRONG columns rather than an exception:

  2. string `"NaN"` in an experimental row — must not displace a good prediction;
  3. camelCase on the property endpoints vs snake_case inside the fate endpoint;
  4. provenance comes from the ENDPOINT, never from a payload field;
  5. Henry's-law units — the wrong reading is a clean 5.006 log-unit offset.

Trap 1 (batch bodies are newline-separated text) is not covered because this module
makes single lookups only.
"""
import sys
import os

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import chem_lookup as cl  # noqa: E402


SEARCH = [{"dtxsid": "DTXSID3022125", "preferredName": "Carbamazepine",
           "casrn": "298-46-4", "smiles": "NC(=O)N1c2ccccc2C=Cc2ccccc21",
           "inchikey": "FFGPTBGBLSHEPO-UHFFFAOYSA-N"}]
DETAIL = {"dtxsid": "DTXSID3022125", "averageMass": 236.269,
          "smiles": "NC(=O)N1c2ccccc2C=Cc2ccccc21"}
# camelCase, and NOTHING in the rows says "experimental" — the URL does (trap 3, 4)
EXPERIMENTAL = [
    {"propName": "LogKow: Octanol-Water", "propValue": 2.45, "propUnit": "",
     "sourceName": "PHYSPROP"},
    {"propName": "Melting Point", "propValue": 190.0, "propUnit": "°C", "sourceName": "PHYSPROP"},
]
PREDICTED = [
    {"propName": "LogP: Octanol-Water", "propValue": 3.39, "propUnit": "",
     "sourceName": "OPERA", "adConclusionGlobal": "Inside", "adValueGlobal": 0.83},
    {"propName": "Henry's Law Constant", "propValue": 1.08e-10, "propUnit": "atm-m3/mol",
     "sourceName": "OPERA", "adConclusionGlobal": "Inside"},
    {"propName": "pKa_a", "propValue": 13.9, "propUnit": "", "sourceName": "OPERA",
     "adConclusionGlobal": "Inside"},
]
# the fate endpoint NESTS its records and uses snake_case, mixing both types (trap 3).
# The VALUE here is the real one the live API returns for carbamazepine: 549.541
# **L/kg**, i.e. LINEAR — not the log10 the property is usually written in.
FATE = {"dtxsid": "DTXSID3022125",
        "predictedFateData": [
            {"prop_name": "Soil Adsorption Coefficient (Koc)", "prop_value": 549.541,
             "prop_unit": "L/kg", "prop_type": "predicted", "source_name": "OPERA",
             "ad_conclusion_global": "Inside"}]}


def fake_get(payloads=None, fail=None):
    """`_get_fn` returning canned payloads keyed by a path fragment."""
    payloads = payloads if payloads is not None else {
        "/chemical/search/": SEARCH,
        "/chemical/detail/": [DETAIL],
        "/chemical/property/experimental/": EXPERIMENTAL,
        "/chemical/property/predicted/": PREDICTED,
        "/chemical/fate/": [FATE],
    }

    def _get(path, key, base_url=cl.DEFAULT_BASE_URL, params=None):
        if fail:
            return None, fail
        for frag, payload in payloads.items():
            if frag in path:
                return payload, None
        return None, "not found (HTTP 404)"
    return _get


# --- base URL --------------------------------------------------------------
def test_base_url_is_the_live_host_not_the_retired_one():
    """`api-ccte.epa.gov` no longer resolves; a default pointing there fails with a
    DNS error that no offline test can catch, so the constant itself is pinned."""
    assert cl.DEFAULT_BASE_URL == "https://comptox.epa.gov/ctx-api"
    assert "api-ccte" not in cl.DEFAULT_BASE_URL
    assert cl.base_url() == cl.DEFAULT_BASE_URL
    assert cl.base_url("https://example.test/x/") == "https://example.test/x"


# --- identifiers -----------------------------------------------------------
def test_identifier_kind_distinguishes_name_cas_and_structure():
    """A NAME must never be read as a structure: 'carbamazepine' parses as a valid
    SMILES to a permissive parser, and mis-classifying it sends the wrong query."""
    assert cl.identifier_kind("carbamazepine") == "name"
    assert cl.identifier_kind("Perfluorooctanoic acid") == "name"
    assert cl.identifier_kind("298-46-4") == "casrn"
    assert cl.identifier_kind("DTXSID3022125") == "dtxsid"
    assert cl.identifier_kind("FFGPTBGBLSHEPO-UHFFFAOYSA-N") == "inchikey"
    pytest.importorskip("rdkit", reason="SMILES detection needs RDKit")
    assert cl.identifier_kind("NC(=O)N1c2ccccc2C=Cc2ccccc21") == "smiles"
    assert cl.identifier_kind("OC(=O)C(F)(F)C(F)(F)F") == "smiles"


def test_smiles_is_resolved_through_its_inchikey():
    """CTX has no SMILES search, so a structure query goes via the InChIKey."""
    pytest.importorskip("rdkit")
    assert cl.smiles_to_inchikey("NC(=O)N1c2ccccc2C=Cc2ccccc21").startswith("FFGPTBGBLSHEPO")
    assert cl.smiles_to_inchikey("not a molecule") == ""


# --- trap 4: provenance is the endpoint, not the payload -------------------
def test_provenance_comes_from_the_endpoint_and_experimental_wins():
    """Neither fixture row says "experimental" anywhere — the only thing that knows
    is which URL returned it. The measured log Kow (2.45) must win over the OPERA
    prediction (3.39), and both must be labelled correctly."""
    r = cl.lookup("carbamazepine", key="k", _get_fn=fake_get())
    assert r.ok and r.dtxsid == "DTXSID3022125"
    p = r.props["log_kow"]
    assert p.value == pytest.approx(2.45) and p.is_experimental
    assert "PHYSPROP" in p.badge() and "predicted" not in p.badge()
    # a property only the predicted endpoint has is still returned, and labelled
    assert r.props["henry"].source == "predicted"
    assert "OPERA" in r.props["henry"].badge()


def test_fate_endpoint_is_read_with_snake_case_and_its_own_type_field():
    """The fate payload NESTS its records and names fields snake_case — and it is
    the one endpoint that carries the type itself. Koc lives there."""
    r = cl.lookup("carbamazepine", key="k", _get_fn=fake_get())
    koc = r.props["koc"]
    assert koc.value == pytest.approx(549.541) and koc.source == "predicted"
    assert r.neutral_kwargs()["Koc"] == pytest.approx(549.541)


# --- trap 5 again: the LOG-vs-LINEAR scale ----------------------------------
def test_koc_is_taken_as_linear_unless_the_row_says_log():
    """The sharpest form of the unit trap, and the one the live API actually sprang:
    CTX sends Koc LINEAR (549.541 L/kg) while the property is normally written
    "log Koc", so assuming log10 is a 10**549 overflow. The scale is read off the
    row (unit first, then the name); a genuinely log-valued row is converted."""
    assert not cl.is_log_scale("L/kg", "Soil Adsorption Coefficient (Koc)")
    assert cl.is_log_scale("Log10 unitless", "LogKow")
    assert cl.to_linear(549.541, "L/kg", "koc") == pytest.approx(549.541)
    assert cl.to_linear(2.4, "log10 L/kg", "koc") == pytest.approx(10 ** 2.4)
    assert cl.to_log10(2.45, "Log10 unitless", "LogKow") == pytest.approx(2.45)
    assert cl.to_log10(282.0, "", "Octanol-Water Partition Coefficient") == pytest.approx(2.45,
                                                                                         abs=1e-2)
    # ...and the REVERSE misread is guarded: a row that merely omits the word "log"
    # must not have its already-log 2.45 turned into 0.389.
    m = cl.LOGKOW_PLAUSIBLE_MAX
    assert cl.to_log10(2.45, "", "Octanol-Water Partition Coefficient", m) == pytest.approx(2.45)
    assert cl.to_log10(282.0, "", "Octanol-Water Partition Coefficient", m) == pytest.approx(
        2.45, abs=1e-2)
    unlogged = [{"propName": "Octanol-Water Partition Coefficient", "propValue": 2.45,
                 "sourceName": "PHYSPROP"}]
    r0 = cl.lookup("carbamazepine", key="k",
                   _get_fn=fake_get({"/chemical/search/": SEARCH, "/chemical/detail/": [DETAIL],
                                     "/chemical/property/experimental/": unlogged}))
    assert r0.props["log_kow"].value == pytest.approx(2.45)
    # end to end: a log-valued fate row lands as the same LINEAR L/kg the model wants
    logged = {"dtxsid": "DTXSID3022125", "predictedFateData": [
        {"prop_name": "log Koc", "prop_value": 2.74, "prop_unit": "log10 L/kg",
         "prop_type": "predicted"}]}
    r = cl.lookup("carbamazepine", key="k",
                  _get_fn=fake_get({"/chemical/search/": SEARCH, "/chemical/detail/": [DETAIL],
                                    "/chemical/fate/": [logged]}))
    assert r.neutral_kwargs()["Koc"] == pytest.approx(10 ** 2.74, rel=1e-6)


def test_an_implausible_koc_is_refused_rather_than_passed_on():
    """A Koc that could only come from a scale misread must not reach the soil model
    (the crash this replaced produced Koc = 4.2e+31 for benzoic acid)."""
    huge = {"dtxsid": "DTXSID3022125", "predictedFateData": [
        {"prop_name": "Soil Adsorption Coefficient (Koc)", "prop_value": 1e30,
         "prop_unit": "L/kg", "prop_type": "predicted"}]}
    r = cl.lookup("carbamazepine", key="k",
                  _get_fn=fake_get({"/chemical/search/": SEARCH, "/chemical/detail/": [DETAIL],
                                    "/chemical/fate/": [huge]}))
    assert "koc" in r.props and "Koc" not in r.neutral_kwargs()


def test_a_vapour_pressure_row_is_not_served_as_the_log_kow():
    """'log p' is a substring of 'log pvap', so a permissive substring match would
    hand the soil/plant model a vapour pressure in the log Kow slot."""
    exp = [{"propName": "LogPvap: Vapor Pressure", "propValue": -7.6, "sourceName": "PHYSPROP"}]
    r = cl.lookup("carbamazepine", key="k",
                  _get_fn=fake_get({"/chemical/search/": SEARCH, "/chemical/detail/": [DETAIL],
                                    "/chemical/property/experimental/": exp,
                                    "/chemical/property/predicted/": PREDICTED}))
    assert r.props["log_kow"].value == pytest.approx(3.39)      # the real prediction, not -7.6


def test_a_pka_row_with_no_stated_centre_is_still_used_but_flagged():
    """OPERA labels its two centres (pKa_a / pKa_b), but a row named just "pKa" is
    a real value whose centre the payload does not state — use it and leave the
    acid/base choice to the caller rather than dropping the compound's ionisation."""
    pred = [{"propName": "pKa", "propValue": 4.2, "sourceName": "OPERA"}]
    r = cl.lookup("benzoic acid", key="k",
                  _get_fn=fake_get({"/chemical/search/": SEARCH, "/chemical/detail/": [DETAIL],
                                    "/chemical/property/predicted/": pred}))
    assert "pka_unlabelled" in r.props
    nk = r.neutral_kwargs()
    assert nk["pKa"] == pytest.approx(4.2) and nk["is_acid"] is True
    # a LABELLED centre wins and the unlabelled fallback stays out of the way
    r2 = cl.lookup("carbamazepine", key="k", _get_fn=fake_get())
    assert "pka_unlabelled" not in r2.props and r2.props["pka_acidic"].value == pytest.approx(13.9)


# --- trap 2: string nulls ---------------------------------------------------
def test_a_string_nan_measurement_never_displaces_a_usable_prediction():
    """float("NaN") SUCCEEDS, so an experimental row of "NaN" would otherwise rank
    first (experimental beats predicted) and blank out the property."""
    for null in ("NaN", "null", "N/A", "-Infinity", ""):
        exp = [{"propName": "LogKow: Octanol-Water", "propValue": null, "sourceName": "PHYSPROP"}]
        r = cl.lookup("carbamazepine", key="k",
                      _get_fn=fake_get({"/chemical/search/": SEARCH,
                                        "/chemical/detail/": [DETAIL],
                                        "/chemical/property/experimental/": exp,
                                        "/chemical/property/predicted/": PREDICTED}))
        p = r.props["log_kow"]
        assert p.value == pytest.approx(3.39), f"{null!r} displaced the prediction"
        assert p.source == "predicted"
    assert cl._num("NaN") is None and cl._num("inf") is None and cl._num(float("nan")) is None
    assert cl._num("2.45") == pytest.approx(2.45)


# --- trap 5: units ----------------------------------------------------------
def test_henry_becomes_a_dimensionless_kaw_or_nothing():
    """Reading Pa-m3/mol as atm-m3/mol is a clean 5.006 log-unit error that still
    looks plausible, so each unit is converted explicitly and an unrecognised one
    returns None — a wrong K_AW silently switches the leaf's volatilisation sink."""
    import math
    assert cl.henry_to_kaw(1.0, "atm-m3/mol") == pytest.approx(1.0 / 0.024465, rel=1e-3)
    pa = cl.henry_to_kaw(1.0, "Pa-m3/mol")
    assert math.log10(cl.henry_to_kaw(1.0, "atm-m3/mol") / pa) == pytest.approx(5.006, abs=1e-3)
    assert cl.henry_to_kaw(0.5, "dimensionless") == pytest.approx(0.5)
    assert cl.henry_to_kaw(1.0, "mmHg") is None
    r = cl.lookup("carbamazepine", key="k", _get_fn=fake_get())
    assert r.props["K_AW"].value == pytest.approx(1.08e-10 / 0.024465, rel=1e-3)


# --- the OPERA applicability domain ----------------------------------------
def test_outside_the_applicability_domain_reaches_the_badge_and_loses_the_tie():
    """A prediction its own model calls out-of-domain is not merely uncertain, so it
    is labelled — and an in-domain candidate for the same property outranks it."""
    pred = [{"propName": "LogP", "propValue": 9.9, "sourceName": "OPERA",
             "adConclusionGlobal": "Outside the global applicability domain"},
            {"propName": "LogP: Octanol-Water", "propValue": 3.39, "sourceName": "OPERA",
             "adConclusionGlobal": "Inside"}]
    r = cl.lookup("carbamazepine", key="k",
                  _get_fn=fake_get({"/chemical/search/": SEARCH, "/chemical/detail/": [DETAIL],
                                    "/chemical/property/predicted/": pred}))
    p = r.props["log_kow"]
    assert p.value == pytest.approx(3.39) and not p.outside_ad
    bad = cl.Prop(value=9.9, source="predicted", origin="OPERA", ad="Outside the domain")
    assert bad.outside_ad and "OUTSIDE" in bad.badge()


# --- the model-facing contract ---------------------------------------------
def test_neutral_kwargs_maps_onto_simulate_neutral_and_omits_half_life():
    """What the lookup hands the model — and what it must NOT: the in-planta
    half-life is not a dashboard property (and per Kodesova not even a compound
    constant), so it is never auto-filled."""
    r = cl.lookup("carbamazepine", key="k", _get_fn=fake_get())
    nk = r.neutral_kwargs()
    assert nk["log_kow"] == pytest.approx(2.45)
    assert nk["MW"] == pytest.approx(236.269, rel=1e-4)
    assert nk["pKa"] == pytest.approx(13.9) and nk["is_acid"] is True
    assert "half_life" not in nk
    api = pytest.importorskip("model_api")
    res = api.simulate_neutral(nk["log_kow"], MW=nk["MW"], K_AW=nk["K_AW"],
                               air=False, half_life=7.0, season=60.0, n_t=61)
    assert res["baf_final"]["root"] > 0


def test_pka_is_absent_when_none_was_found():
    """`pKa=None` is the strictly-neutral, bit-identical path — a missing pKa must
    stay missing rather than defaulting to a number."""
    r = cl.lookup("carbamazepine", key="k",
                  _get_fn=fake_get({"/chemical/search/": SEARCH, "/chemical/detail/": [DETAIL],
                                    "/chemical/property/experimental/": EXPERIMENTAL}))
    assert "pKa" not in r.neutral_kwargs()


# --- the pKa is the one field that changes WHICH MODEL RUNS ----------------
def test_speciation_note_agrees_with_the_model_and_separates_its_two_warnings():
    """A filled pKa reroutes the run onto the weak-electrolyte path, so it is
    reported apart from the generic provenance badge — and the two things that can
    be wrong are INDEPENDENT:

    * `f_n` below the tested floor is a MODEL-scope statement, true however good the
      pKa is — benzoic acid's own MEASURED 4.18 lands at f_n 0.005;
    * a PREDICTED pKa is a provenance problem — carbamazepine is the in-repo
      counterexample, where OPERA's 5.07 and the measured 13.9 that this repo's
      Kodešová 2019 table is built on are nine log units apart.
    """
    LP = pytest.importorskip("literature_params")
    for pka, acid in ((5.07, True), (13.9, True), (4.18, True), (9.0, False)):
        f_n, _ = cl.speciation_note(pka, acid)
        assert f_n == pytest.approx(LP.speciation(pka, cl.PH_ROOT_ZONE, acid)[0], abs=1e-12)

    car_opera, txt_pred = cl.speciation_note(5.07, True, predicted=True)
    car_meas, _ = cl.speciation_note(13.9, True, predicted=True)
    assert car_opera == pytest.approx(0.0358, abs=1e-3) and car_meas == pytest.approx(1.0)
    assert "PREDICTED" in txt_pred and "BELOW" in txt_pred        # both fire here

    # a MEASURED pKa still trips the floor (it is about the model, not the source)
    benzoic, txt_meas = cl.speciation_note(4.18, True, predicted=False)
    assert benzoic < cl.F_N_TESTED_FLOOR
    assert "BELOW" in txt_meas and "PREDICTED" not in txt_meas

    # ...and a well-inside pKa from a measurement gets neither warning
    _, quiet = cl.speciation_note(13.9, True, predicted=False)
    assert "WARNING" not in quiet


def test_the_cli_note_reads_the_pka_provenance_from_the_lookup():
    """The predicted-pKa warning must key on where the pKa came from, not on
    whether anything else in the record was predicted."""
    r = cl.lookup("carbamazepine", key="k", _get_fn=fake_get())   # pKa_a is OPERA
    assert r.props["pka_acidic"].source == "predicted"
    assert "PREDICTED" in cl._pka_note(r, r.neutral_kwargs())

    exp_pka = [{"propName": "pKa_a", "propValue": 4.18, "sourceName": "OPERA"}]
    r2 = cl.lookup("benzoic acid", key="k",
                   _get_fn=fake_get({"/chemical/search/": SEARCH, "/chemical/detail/": [DETAIL],
                                     "/chemical/property/experimental/": exp_pka}))
    assert r2.props["pka_acidic"].source == "experimental"
    note = cl._pka_note(r2, r2.neutral_kwargs())
    assert "PREDICTED" not in note and "BELOW" in note


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
    yields a MISSING property, never a wrong one in the right slot. Matching is by
    name substring on purpose — the vocabulary is not stable enough to hard-code —
    so the guarantee is one-directional: an unrelated name is ignored."""
    r = cl.lookup("carbamazepine", key="k",
                  _get_fn=fake_get({"/chemical/search/": {"data": SEARCH},
                                    "/chemical/detail/": {"data": [DETAIL]},
                                    "/chemical/property/experimental/":
                                        [{"propName": "Bioconcentration factor",
                                          "propValue": 9.9}]}))
    assert r.ok and r.mol_weight == pytest.approx(236.269, rel=1e-4)
    assert "log_kow" not in r.props and r.neutral_kwargs().get("log_kow") is None


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
