"""The FastAPI web backend (server/) -- the React front end's API.

What is pinned here is the CONTRACT that lets the web UI replace the Streamlit one
without changing a single number: the default scenario returns exactly what
`model_api.simulate()` returns, every data source the Streamlit sidebar offers is
reachable, figures/downloads come out of the same `plots`/`model_api` builders,
and the one deliberate difference (driver dicts keep the ORYZA `leaf_loss`, which
the Streamlit cache tuple drops) is real and in the documented direction.
"""
import json
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))     # the `server` package

fastapi = pytest.importorskip("fastapi")
from fastapi.testclient import TestClient  # noqa: E402

import model_api as api  # noqa: E402


@pytest.fixture(scope="module")
def client():
    import os
    os.environ["PFAS_WEB_WARMUP"] = "0"
    from server.app import app
    with TestClient(app) as c:
        yield c


def test_meta_lists_the_model_vocabulary(client):
    m = client.get("/api/meta").json()
    assert m["congeners"] == api.CONGENERS
    assert m["default_uptake"] == api.DEFAULT_UPTAKE
    assert set(m["labels"]) == {"en", "ko"}
    assert len(m["scenarios_ko"]) == 4 and len(m["presets"]["ko"]) == 3
    assert m["strings"]["disclaimer"]["ko"]
    assert m["default_k_leach"]["PFBA"] == pytest.approx(api.default_k_leach("PFBA"))


def test_default_run_is_bit_identical_to_model_api(client):
    d = client.post("/api/run", json={}).json()
    ref = api.simulate("PFOA")
    assert d["congener"] == "PFOA"
    for k in api.TISSUES:
        assert d["baf_final"][k] == ref["baf_final"][k]
        assert np.array_equal(np.asarray(d["conc"][k], float), ref["conc"][k])
    assert d["straw_baf"] == ref["straw_baf"]
    assert d["finals"]["grain"] == float(ref["conc"]["grain"][-1])
    assert d["intake"] == pytest.approx(api.intake_fraction(float(ref["conc"]["grain"][-1]),
                                                            congener="PFOA"), nan_ok=True)
    ms = api.metric_series(ref, "conc")
    assert d["schematic"]["conc"]["cmax"] == ms["cmax"]


def test_mechanism_switches_reach_the_model(client):
    d = client.post("/api/run", json={"compound": {"congener": "PFOS"},
                                      "mechanism": {"uptake": "bypass", "lipid_loading": True,
                                                    "E_m": -130}}).json()
    ref = api.simulate("PFOS", uptake="bypass", lipid_loading=True, E_m_mV=-130.0)
    assert d["baf_final"]["root"] == ref["baf_final"]["root"]
    assert d["params"]["uptake"] == "bypass"


@pytest.mark.parametrize("mode", ["tables", "csv", "inventory", "biomon"])
def test_every_data_source_runs(client, mode):
    r = client.post("/api/run", json={"exposure": {"mode": mode}})
    assert r.status_code == 200, r.text
    d = r.json()
    assert np.isfinite(d["finals"]["grain"]) and d["effective"]["drivers"] == (mode != "biomon")
    if mode == "biomon":
        assert d["bio_baf"] and d["measured_bio"]["Cwo"] == 1.0
    if mode == "inventory":
        assert d["has_isotherm"]


def test_bad_tables_fall_back_with_a_message(client):
    d = client.post("/api/run", json={"lang": "ko", "exposure": {
        "mode": "tables", "tables": {"cwo": {"day": [0], "Cwo": [1]}}}}).json()
    assert d["effective"]["drivers"] is False
    assert d["messages"] and "표를 읽지 못해" in d["messages"][0]["text"]


def test_driver_modes_keep_the_oryza_leaf_loss(client):
    """The deliberate difference from the Streamlit cache tuple: the soil-inventory
    drivers carry ORYZA's leaf senescence rate and it must reach the model, so the
    leaf is not spuriously concentrated by a shrinking mass."""
    from server import service as S
    from server.schemas import ExposureSpec, Scenario
    b = S.build(Scenario(exposure=ExposureSpec(mode="inventory")))
    assert "leaf_loss" in b.drivers
    d = client.post("/api/run", json={"exposure": {"mode": "inventory"}}).json()
    stripped = api.simulate("PFOA", drivers={k: b.drivers[k] for k in ("t", "Cwo", "Qtp", "M")})
    assert d["baf_final"]["leaf"] < stripped["baf_final"]["leaf"]


def test_neutral_compound_class(client):
    d = client.post("/api/run", json={"compound": {"kind": "neutral",
                                                   "neutral": {"log_kow": 2.45}}}).json()
    ref = api.simulate_neutral(2.45, name="carbamazepine", half_life=7.0, g_apo=0.0)
    assert d["kind"] == "neutral" and d["obs"] == {}
    assert d["TSCF"] == ref["TSCF"] and d["baf_final"]["root"] == ref["baf_final"]["root"]
    assert d["koc"] == pytest.approx(float(__import__("literature_params").koc_neutral(2.45)))


def test_smiles_structure(client):
    if not api.rdkit_available():
        pytest.skip("RDKit not installed")
    pfoa = "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"
    d = client.post("/api/run", json={"compound": {"kind": "smiles", "smiles": pfoa}}).json()
    assert d["congener"] == "PFOA" and d["descriptors"]["matched_name"] == "PFOA"
    assert d["baf_final"]["root"] == api.simulate("PFOA")["baf_final"]["root"]
    bad = client.post("/api/run", json={"compound": {"kind": "smiles", "smiles": "xx"}})
    assert bad.status_code == 400


@pytest.mark.parametrize("name,opts", [
    ("tissue", {}), ("burden", {}), ("drivers", {}), ("soil_profile", {}),
    ("buildup_plain", {}), ("where_plain", {}), ("baf", {}), ("chain", {"key": "B_root"}),
    ("intake_gauge", {"rice_g": 200, "bw": 50}), ("compare", {"congeners": ["PFBA", "PFOA"]}),
])
def test_figures_build(client, name, opts):
    r = client.post("/api/figure", json={"name": name, "opts": opts, "lang": "ko"})
    assert r.status_code == 200, r.text
    out = r.json()
    assert out["figure"]["data"], name
    if name == "intake_gauge":
        assert out["data"]["info"]["rice_intake_g_day"] == 200


def test_unknown_figure_is_a_404(client):
    assert client.post("/api/figure", json={"name": "nope"}).status_code == 404


def test_inverse_streams_progress_then_a_result(client):
    truth = api.simulate("PFOA", Cwo=2.0)
    meas = {"root": float(truth["conc"]["root"][-1]), "grain": float(truth["conc"]["grain"][-1])}
    lines = []
    with client.stream("POST", "/api/inverse", json={"measured": meas}) as r:
        for line in r.iter_lines():
            if line:
                lines.append(json.loads(line))
    assert lines[0]["type"] == "progress" and lines[-1]["type"] == "result"
    est = lines[-1]["estimate"]
    assert est["median"] == pytest.approx(2.0, rel=0.05)
    assert est["ci95"][0] < 2.0 < est["ci95"][1]
    assert lines[-1]["figure"]["data"]


def test_inverse_refuses_a_smiles_compound(client):
    with client.stream("POST", "/api/inverse", json={
            "scenario": {"compound": {"kind": "smiles", "smiles": "OC(=O)C(F)(F)F"}},
            "measured": {"root": 1.0}}) as r:
        msg = json.loads(next(r.iter_lines()))
    assert msg["type"] == "error"


@pytest.mark.parametrize("kind", ["summary_csv", "timeseries_csv", "summary_md", "map_html"])
def test_downloads(client, kind):
    r = client.post(f"/api/download/{kind}", json={"lang": "ko"})
    assert r.status_code == 200 and "attachment" in r.headers["content-disposition"]
    if kind == "summary_csv":
        assert r.text == api.summary_csv(api.simulate("PFOA"), api.observed_baf("PFOA"))


def test_small_endpoints(client):
    assert client.get("/api/k-leach", params={"congener": "PFBA"}).json()["k_leach"] == \
        pytest.approx(api.default_k_leach("PFBA"))
    t = client.get("/api/tables/defaults", params={"season": 120, "Cwo": 2}).json()
    assert t["cwo"]["Cwo"] == [2.0] * 4 and len(t["growth"]["day"]) == 7
    assert client.get("/api/health").json()["ok"]


def test_scenario_rejects_unknown_fields(client):
    assert client.post("/api/run", json={"exposure": {"bogus": 1}}).status_code == 422
