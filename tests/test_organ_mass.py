"""Organ mass basis: crop-model DRY matter -> the FRESH mass the plant ODE runs on.

The biomass drivers (ORYZA2000, growth_rice, measured tables) all produce dry matter;
PlantInputs.M and B_k are fresh weight. Until `organ_mass` existed the dry matter went
into M unconverted, inflating every shoot concentration 4-6x relative to the root.
These tests lock in the conversion, its single source of theta_fw (the binding factor's),
and that mass_basis="dry" still reproduces the numbers published before the fix.
"""
import json
import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

import organ_mass as om  # noqa: E402
import model_api as api  # noqa: E402
import oryza_growth as og  # noqa: E402


def test_theta_is_the_binding_factors_theta():
    comp = json.load(open(os.path.join(ROOT, "params", "parameters.json")))["tissue_composition_recommended"]
    for organ, key in (("root", "root"), ("stem", "stem"), ("leaf", "leaf"), ("grain", "grain_brown")):
        assert om.THETA_FW[organ] == comp[key]["theta_fw"]
        assert om.fresh_factor(organ) == pytest.approx(1.0 / (1.0 - comp[key]["theta_fw"]))
        assert om.fresh_factor(organ, "dry") == 1.0
    with pytest.raises(ValueError):
        om.fresh_factor("leaf", "wet")


def test_to_model_mass_keeps_relative_rates():
    b = {"root": np.ones(3), "stem": np.ones(3), "leaf": np.ones(3), "grain": np.ones(3),
         "leaf_death_rate": np.full(3, 0.02)}
    m = om.to_model_mass(b)
    assert np.allclose(m["leaf"], 1.0 / (1.0 - om.THETA_FW["leaf"]))
    assert m["leaf_death_rate"] is b["leaf_death_rate"]          # a 1/day rate: untouched
    assert np.allclose(om.to_model_mass(b, "dry")["root"], 1.0)


def test_default_simulate_runs_on_fresh_organ_mass():
    t = np.linspace(0.0, 120.0, 241)
    dry = og.organ_biomass_oryza(t, p=og.OryzaParams(season=120.0))
    r = api.simulate("PFOA")
    assert r["params"]["mass_basis"] == "fresh"
    for i, organ in enumerate(api.TISSUES):
        expect = np.maximum(dry[organ] * om.fresh_factor(organ), 1e-4)
        assert np.allclose(r["M"][:, i], expect)


def test_dry_basis_reproduces_the_pre_fix_numbers():
    r = api.simulate("PFOA", mass_basis="dry")
    f = r["baf_final"]
    assert f["root"] == pytest.approx(0.479, rel=2e-3)
    assert f["leaf"] == pytest.approx(2.5225, rel=2e-3)
    assert f["grain"] == pytest.approx(0.1476, rel=2e-3)


def test_oryza_refit_matches_its_own_mass_basis():
    """f_xy_source='oryza' has one fit per basis: *_oryza (fresh, the default) and
    *_oryza_dry (the pre-fix fit). Each reproduces Yamazaki PFOA on its own basis."""
    o = api.observed_baf("PFOA")
    for mb, key in (("fresh", "f_xy_oryza"), ("dry", "f_xy_oryza_dry")):
        r = api.simulate("PFOA", f_xy_source="oryza", mass_basis=mb)
        assert r["params"]["f_xy"] == pytest.approx(api._CONG["PFOA"][key])
        assert abs(np.log10(r["straw_baf"]) - np.log10(o["straw"])) < 0.05
        assert abs(np.log10(r["baf_final"]["grain"]) - np.log10(o["grain"])) < 0.05
    assert api._CONG["PFOA"]["f_xy_oryza"] > 3 * api._CONG["PFOA"]["f_xy_oryza_dry"]


def test_fix_moves_the_shoot_not_the_root():
    d = api.simulate("PFOA", mass_basis="dry")["baf_final"]
    w = api.simulate("PFOA")["baf_final"]
    assert w["root"] == pytest.approx(d["root"], rel=0.05)        # uptake per unit mass: ~unchanged
    assert w["leaf"] < d["leaf"] / 3.0                            # terminal accumulator: ~1/(1-0.78)
    assert w["grain"] < d["grain"] / 3.0                          # phloem-fed from the diluted leaf


def test_placeholder_forcing_is_untouched():
    r = api.simulate("PFOA", measured_forcing=False)
    assert r["params"]["mass_basis"] is None                      # reproduce_demo's drivers have no basis
    assert r["baf_final"]["grain"] == pytest.approx(0.7532, rel=2e-3)


def test_driver_builders_convert():
    t = np.linspace(0.0, 120.0, 61)
    _, M = api.measured_forcing(t)
    _, Md = api.measured_forcing(t, mass_basis="dry")
    assert np.allclose(M[-1] / Md[-1], [om.fresh_factor(o) for o in om.ORGANS])
    drv = og.oryza_drivers(n_t=61)
    drv_d = og.oryza_drivers(n_t=61, mass_basis="dry")
    assert drv["M"][-1, 2] / drv_d["M"][-1, 2] == pytest.approx(om.fresh_factor("leaf"))


def test_measured_biomass_table_is_dry_in_fresh_out():
    import io
    import measured_biomass as mb
    table = io.StringIO("day,root,stem,leaf,grain\n0,1,1,1,0\n60,10,20,15,0\n120,12,40,10,50\n")
    drv = mb.biomass_drivers(mb.load_biomass_table(table), "g/hill", n_t=13)
    assert drv["M"][-1, 1] == pytest.approx(0.040 * om.fresh_factor("stem"), rel=1e-6)
    assert drv["M"][-1, 3] == pytest.approx(0.050 * om.fresh_factor("grain"), rel=1e-6)


def test_every_copy_of_theta_agrees():
    """theta_fw is read in three places -- the binding factor (model_api._COMP), the
    dry->fresh mass conversion (organ_mass.THETA_FW) and the neutral path's own
    hard-coded neutral_dpu.RICE_WATER. They must not drift apart."""
    import neutral_dpu as ND
    for organ, key in (("root", "root"), ("stem", "stem"), ("leaf", "leaf"),
                       ("grain", "grain_brown")):
        assert ND.RICE_WATER[organ] == om.THETA_FW[organ] == api._COMP[key]["theta_fw"]


def test_the_xylem_split_follows_the_mass_basis_of_M():
    """The leaf/grain xylem split is the one place that needs DRY matter, so M carries
    its basis: builders tag their drivers, an untagged M is fresh (§8), the built-in
    measured forcing is on `mass_basis`, and the placeholder logistic is basis-free."""
    from pfas_rice_plant_module_4pool_surf import PlantInputs
    t = np.linspace(0.0, 120.0, 241)
    base = api.simulate("PFOA")["baf_final"]
    d = api.drivers_from_arrays(t, np.ones_like(t))          # crop M -> tagged on its basis
    assert d["mass_basis"] == "fresh"
    assert api.simulate("PFOA", drivers=d)["baf_final"]["grain"] == pytest.approx(base["grain"])
    untagged = {k: v for k, v in d.items() if k != "mass_basis"}
    assert api.simulate("PFOA", drivers=untagged)["baf_final"] == api.simulate(
        "PFOA", drivers=d)["baf_final"]
    # at the shipped theta S_XYLEM_DRY reproduces the old fresh-mass split by construction
    as_given = api.simulate("PFOA", drivers={**d, "mass_basis": None})["baf_final"]
    assert as_given["grain"] == pytest.approx(base["grain"], rel=1e-4)
    # where the tag matters: dry-matter M (the pre-fix records) must be split as given
    dd = api.drivers_from_arrays(t, np.ones_like(t), mass_basis="dry")
    assert dd["mass_basis"] == "dry"
    dry = api.simulate("PFOA", drivers=dd)["baf_final"]["grain"]
    assert dry == pytest.approx(api.simulate("PFOA", mass_basis="dry")["baf_final"]["grain"])
    # what the tag buys: with M fresh, the split sees only the DRY matter, so changing an
    # organ's water (theta) does not move it; split as given, it would
    from dataclasses import replace
    from pfas_rice_plant_module_4pool_surf import xylem_leaf_share, FRUIT
    Md = om.organ_matrix(api._biomass_fn("oryza")(t, 120.0), "dry")[-1]
    share = {}
    for th in (0.14, 0.60):
        comps = [replace(c, theta=th) if i == FRUIT else c
                 for i, c in enumerate(api._compartments())]
        Mf = Md / (1.0 - np.array([c.theta for c in comps]))
        share[th] = (xylem_leaf_share(comps, Mf, "fresh"), xylem_leaf_share(comps, Mf, None))
    assert share[0.60][0] == pytest.approx(share[0.14][0], rel=1e-12)
    assert 1.0 - share[0.60][1] > 1.5 * (1.0 - share[0.14][1])   # grain share 0.118 -> 0.223
    with pytest.raises(ValueError):
        PlantInputs(t=t, Cwo=np.ones_like(t), Qtp=np.ones_like(t), M=d["M"], M_basis="wet")
