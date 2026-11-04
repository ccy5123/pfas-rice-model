"""The FastAPI web backend (server/) -- the React front end's API.

What is pinned here is the CONTRACT that lets the web UI replace the Streamlit one
without changing a single number: the default scenario returns exactly what
`model_api.simulate()` returns, every data source the Streamlit sidebar offers is
reachable, the transport variant and its fall-backs are the sidebar's, figures and
downloads come out of the same `plots` / `model_api` / `report_html` builders, and
the intake reads white rice as the Streamlit Simple view does.
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
    os.environ["PFAS_WEB_WARMUP"] = "0"
    from server.app import app
    with TestClient(app) as c:
        yield c


def test_meta_lists_the_model_vocabulary(client):
    m = client.get("/api/meta").json()
    assert m["congeners"] == api.CONGENERS
    assert m["default_shoot"] == api.DEFAULT_SHOOT
    assert m["default_f_xy_source"] == api.DEFAULT_F_XY_SOURCE
    assert m["default_lipid_loading"] == api.DEFAULT_LIPID_LOADING
    assert m["default_variant"] == api.variant_key()
    assert m["withdrawn_variants"] == api.WITHDRAWN_VARIANTS
    assert [k for k, _ in m["crops"]][0] == m["default_crop"] == "site2020"
    assert m["site_season"] == pytest.approx(api.site_season())
    assert set(m["labels"]) == {"en", "ko"}
    assert len(m["scenarios_ko"]) == 4 and len(m["presets"]["ko"]) == 3
    assert m["strings"]["disclaimer"]["ko"]
    assert m["intake"]["rice_g"] == api.DEFAULT_RICE_INTAKE_G_PER_DAY
    assert m["apriori_log10_rmse"] == api.APRIORI_LOG10_RMSE
    # nothing of the pre-amendment-1 API is advertised any more
    for gone in ("uptake_modes", "default_uptake", "default_k_leach", "tissue_density",
                 "organ_labels", "tang_refit_fxy"):
        assert gone not in m


def test_default_run_is_bit_identical_to_model_api(client):
    d = client.post("/api/run", json={}).json()
    ref = api.simulate("PFOA")
    assert d["congener"] == "PFOA"
    for k in api.TISSUES:
        assert d["baf_final"][k] == ref["baf_final"][k]
        assert np.array_equal(np.asarray(d["conc"][k], float), ref["conc"][k])
    assert d["straw_baf"] == ref["straw_baf"]
    assert d["baf_final_dw"] == ref["baf_final_dw"]
    assert d["finals"]["grain"] == float(ref["conc"]["grain"][-1])
    # the intake is eaten white rice, as the Streamlit Simple view computes it
    white = api.white_rice_conc(ref)
    assert d["finals"]["white"] == white
    assert d["intake"] == pytest.approx(api.intake_fraction(white, congener="PFOA"), nan_ok=True)
    assert d["variant"] == {"key": api.variant_key(), "withdrawn": None}
    ms = api.metric_series(ref, "conc")
    assert d["schematic"]["conc"]["cmax"] == ms["cmax"]


def test_variant_switches_reach_the_model(client):
    d = client.post("/api/run", json={"compound": {"congener": "PFOS"},
                                      "mechanism": {"lipid_loading": True,
                                                    "fxy_source": "recommended"}}).json()
    ref = api.simulate("PFOS", lipid_loading=True, f_xy_source="recommended")
    assert d["baf_final"]["root"] == ref["baf_final"]["root"]
    assert d["params"]["lipid_loading"] is True and d["variant"]["key"] == "redistributed_lipid"
    assert not d["messages"]


@pytest.mark.parametrize("mech,variant", [({"shoot": "4pool"}, "4pool_free"),
                                          ({"lipid_loading": True}, "redistributed_lipid")])
def test_fit_falls_back_as_the_sidebar_does(client, mech, variant):
    """f_xy 'fit' exists for the free anion on the redistributed shoot only: the sidebar
    runs 'recommended' with a warning, and a withdrawn variant is flagged."""
    d = client.post("/api/run", json={"mechanism": mech}).json()
    assert d["variant"]["key"] == variant
    assert d["messages"] and d["messages"][0]["level"] == "warning"
    assert d["variant"]["withdrawn"] == api.WITHDRAWN_VARIANTS.get(variant)
    assert d["params"]["withdrawn"] == api.WITHDRAWN_VARIANTS.get(variant)


def test_removed_parameters_are_refused(client):
    for body in ({"mechanism": {"E_m": -120}}, {"mechanism": {"uptake": "carrier"}},
                 {"mechanism": {"km_scale": 1}}, {"exposure": {"k_leach": 0.02}},
                 {"exposure": {"measured_forcing": False}},
                 {"exposure": {"tables": {"density": {"root": 1.0}}}}):
        assert client.post("/api/run", json=body).status_code == 422, body


def test_generic_crop_and_flooded_shape(client):
    d = client.post("/api/run", json={"mechanism": {"crop": "oryza2000"},
                                      "exposure": {"season": 130, "cwo_profile": "flooded"}}).json()
    ref = api.simulate("PFOA", scenario="generic", biomass="oryza2000", season=130.0,
                       cwo_profile="flooded")
    assert d["effective"]["season"] == 130.0
    assert np.array_equal(np.asarray(d["Cwo"], float), ref["Cwo"])
    assert d["baf_final_dw"] == ref["baf_final_dw"]


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
    if mode == "tables":                  # both interpolation forms ran (PREREG 7.3)
        assert set(d["interp_spread"]) == {"root", "straw", "grain"}


def test_inventory_seeds_are_the_sidebar_seeds(client):
    import soil_hydrus as sh
    seed = client.get("/api/inventory/defaults", params={"congener": "PFOA"}).json()
    assert seed["C_total"] == pytest.approx(sh.post_season_soil("PFOA")[0] / 1000.0)
    assert seed["K_F"] == pytest.approx(api.site_kd("PFOA"))
    d = client.post("/api/run", json={"exposure": {"mode": "inventory"}}).json()
    t = np.linspace(0.0, api.site_season(), 241)
    Cwo, _ = api.pore_water_from_inventory(t, seed["C_total"], K_F=seed["K_F"], n=1.0,
                                           congener="PFOA")
    ref = api.simulate("PFOA", drivers=api.drivers_from_arrays(
        t, Cwo, season=api.site_season(), scenario="yamazaki2020", biomass="oryza2000"))
    assert d["baf_final_dw"] == ref["baf_final_dw"]


def test_bad_tables_fall_back_with_a_message(client):
    d = client.post("/api/run", json={"lang": "ko", "exposure": {
        "mode": "tables", "tables": {"cwo": {"day": [0], "Cwo": [1]}}}}).json()
    assert d["effective"]["drivers"] is False
    assert d["messages"] and "표를 읽지 못해" in d["messages"][0]["text"]


def test_neutral_compound_class(client):
    d = client.post("/api/run", json={"compound": {"kind": "neutral",
                                                   "neutral": {"log_kow": 2.45}}}).json()
    ref = api.simulate_neutral(2.45, name="carbamazepine", half_life=7.0, g_apo=0.0,
                               scenario="yamazaki2020", biomass="oryza2000",
                               Koc=float(__import__("literature_params").koc_neutral(2.45)))
    assert d["kind"] == "neutral" and d["obs"] == {} and d["variant"] is None
    assert d["TSCF"] == ref["TSCF"] and d["baf_final"]["root"] == ref["baf_final"]["root"]
    assert d["koc"] == pytest.approx(float(__import__("literature_params").koc_neutral(2.45)))
    # no white-rice value on the neutral path: the intake falls back to the brown rice
    assert d["finals"]["white"] == d["finals"]["grain"]
    off = client.post("/api/run", json={"compound": {"kind": "neutral",
                                                     "neutral": {"grain_feed": None}}}).json()
    assert off["grain_predicted"] is False


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
    ("congener_compare", {}), ("tang_tf", {"dose": "low"}),
])
def test_figures_build(client, name, opts):
    r = client.post("/api/figure", json={"name": name, "opts": opts, "lang": "ko"})
    assert r.status_code == 200, r.text
    out = r.json()
    assert out["figure"]["data"], name
    if name == "intake_gauge":
        assert out["data"]["info"]["rice_intake_g_day"] == 200
        assert out["data"]["white"] == api.white_rice_conc(api.simulate("PFOA"))


def test_yamazaki_observable_follows_the_view(client):
    """Expert: the observable follows the f_xy source and the mechanism switches; Simple:
    the f_xy source only (ui/expert, ui/simple)."""
    body = {"name": "baf_observable", "opts": {"view": "expert"},
            "scenario": {"mechanism": {"lipid_loading": True, "fxy_source": "recommended"}}}
    out = client.post("/api/figure", json=body).json()
    assert out["data"]["mech"] == {"f_xy_source": "recommended", "lipid_loading": True,
                                   "shoot": "redistributed"}
    import plots
    ref = {"congener": "PFOA", "baf_final_dw": api.yamazaki_harvest_baf("PFOA", **out["data"]["mech"])}
    assert out["figure"] == json.loads(plots.fig_baf(ref, api.observed_baf("PFOA")).to_json())
    body["opts"]["view"] = "simple"
    assert client.post("/api/figure", json=body).json()["data"]["mech"] == {
        "f_xy_source": "recommended"}


def test_unknown_figure_is_a_404(client):
    assert client.post("/api/figure", json={"name": "nope"}).status_code == 404


def test_cwo_preview_is_the_stored_flooded_shape(client):
    import plots
    out = client.post("/api/cwo-preview", json={"compound": {"congener": "PFBA"}, "level": 2.0}).json()
    ref = plots.fig_cwo_profile("PFBA", level=2.0, profile="flooded", season=api.site_season())
    assert out["figure"] == json.loads(ref.to_json())


def test_inverse_streams_progress_then_a_result(client):
    truth = api.simulate("PFOA", Cwo=2.0)
    meas = {"root": float(truth["conc_dw"]["root"][-1]),
            "grain": float(truth["conc_dw"]["grain"][-1])}
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


def test_inverse_rsd_is_the_lognormal_sigma():
    from server import service as S
    from ui.common import meas_sigma_log10
    for rsd in (0.0, 20.0, 150.0):
        assert S.meas_sigma_log10(rsd) == meas_sigma_log10(rsd)


@pytest.mark.parametrize("kind", ["summary_csv", "timeseries_csv", "summary_md", "summary_html",
                                  "map_html"])
def test_downloads(client, kind):
    r = client.post(f"/api/download/{kind}", json={"lang": "ko"})
    assert r.status_code == 200 and "attachment" in r.headers["content-disposition"]
    if kind == "summary_csv":
        assert r.text == api.summary_csv(api.simulate("PFOA"), api.observed_baf("PFOA"))
    if kind == "summary_html":
        assert "낟알(현미, 수확 시)" in r.text


def test_small_endpoints(client):
    t = client.get("/api/tables/defaults", params={"Cwo": 2}).json()
    assert t["cwo"]["Cwo"] == [2.0] * 4 and len(t["growth"]["day"]) == 7
    assert t["growth"]["day"][-1] == round(api.site_season())
    g = client.get("/api/tables/defaults", params={"crop": "growth_rice", "season": 100}).json()
    assert g["growth"]["day"][-1] == 100.0
    assert client.get("/api/inventory/defaults").json()["K_F"] == 1.0          # a SMILES compound
    assert client.get("/api/health").json()["ok"]


def test_scenario_rejects_unknown_fields(client):
    assert client.post("/api/run", json={"exposure": {"bogus": 1}}).status_code == 422


def test_streamlit_summary_uses_the_shared_builder():
    """ui/common._summary_html delegates to report_html.summary_html (one builder for both
    front ends), with main's white-rice intake and brown-rice label."""
    pytest.importorskip("streamlit")
    from types import SimpleNamespace
    import report_html
    from ui import common
    res = api.simulate("PFOA")
    cfg = SimpleNamespace(res=res, congener="PFOA", preset_label="중간")
    body, why = common._summary_html(cfg, lang="ko")
    assert why is None and "낟알(현미, 수확 시)" in body.decode("utf-8")
    ref = report_html.summary_html(res, congener="PFOA", disclaimer=common._DISCLAIMER_KO,
                                   scenario="중간", lang="ko")
    strip = lambda b: __import__("re").sub(rb"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}",
                                           b"", b)                          # noqa: E731
    assert strip(body) == strip(ref)
    white = api.white_rice_conc(res)
    assert f"{api.intake_fraction(white, congener='PFOA')['percent']:.0f}%" in body.decode("utf-8")
