"""Organ mass basis: crop-model DRY matter -> the FRESH mass the plant ODE runs on.

The biomass drivers (the ORYZA2000 port, growth_rice, measured tables) produce dry matter;
PlantInputs.M and B_k are fresh weight. These tests lock in the conversion, its single
source of theta_fw (the binding factor's), and the one place that needs dry matter back --
the leaf/grain xylem split.
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


def test_to_model_mass_keeps_relative_rates():
    b = {"root": np.ones(3), "stem": np.ones(3), "leaf": np.ones(3), "grain": np.ones(3),
         "leaf_death_rate": np.full(3, 0.02)}
    m = om.to_model_mass(b)
    assert np.allclose(m["leaf"], 1.0 / (1.0 - om.THETA_FW["leaf"]))
    assert m["leaf_death_rate"] is b["leaf_death_rate"]          # a 1/day rate: untouched


def test_default_simulate_runs_on_fresh_site_year_mass():
    import site_scenario as ss
    run = ss.yamazaki2020_mean_run()
    r = api.simulate("PFOA")
    assert r["params"]["mass_basis"] == "fresh" and r["params"]["scenario"] == "yamazaki2020"
    for i, organ in enumerate(api.TISSUES):
        expect = np.maximum(run[organ] * om.fresh_factor(organ), om.MASS_FLOOR)
        assert np.allclose(r["M"][:, i], expect)


def test_placeholder_forcing_is_refused():
    with pytest.raises(ValueError):                 # the placeholder logistics were removed
        api.simulate("PFOA", measured_forcing=False)


def test_driver_builders_convert():
    t = np.linspace(0.0, 120.0, 61)
    _, M = api.measured_forcing(t, scenario="generic", biomass="oryza")
    dry = og.organ_biomass_oryza(t, p=og.OryzaParams(season=120.0))
    assert np.allclose(M[-1], [max(dry[o][-1] * om.fresh_factor(o), om.MASS_FLOOR) for o in om.ORGANS])
    drv = og.oryza_drivers(n_t=61)
    assert drv["mass_basis"] == "fresh"
    assert drv["M"][-1, 2] == pytest.approx(dry["leaf"][-1] * om.fresh_factor("leaf"), rel=1e-6)


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
    its basis: builders tag their drivers (fresh), an untagged M is fresh (§8), and the
    placeholder logistic is basis-free (None, split on M as given)."""
    from dataclasses import replace
    from pfas_rice_plant_module_4pool_surf import PlantInputs, xylem_leaf_share, FRUIT
    d = api.scenario_drivers()
    assert d["mass_basis"] == "fresh"
    untagged = {k: v for k, v in d.items() if k != "mass_basis"}
    assert api.simulate("PFOA", drivers=untagged)["baf_final"] == api.simulate(
        "PFOA", drivers=d)["baf_final"]
    # what the tag buys: with M fresh the split sees only the DRY matter, so changing an
    # organ's water (theta) does not move it; split as given, it would
    t = d["t"]
    Md = np.array([d["M"][-1, i] / om.fresh_factor(o) for i, o in enumerate(om.ORGANS)])
    share = {}
    for th in (0.14, 0.60):
        comps = [replace(c, theta=th) if i == FRUIT else c
                 for i, c in enumerate(api._compartments())]
        Mf = Md / (1.0 - np.array([c.theta for c in comps]))
        share[th] = (xylem_leaf_share(comps, Mf, "fresh"), xylem_leaf_share(comps, Mf, None))
    assert share[0.60][0] == pytest.approx(share[0.14][0], rel=1e-12)
    assert 1.0 - share[0.60][1] > 1.5 * (1.0 - share[0.14][1])
    with pytest.raises(ValueError):
        PlantInputs(t=t, Cwo=np.ones_like(t), Qtp=np.ones_like(t), M=d["M"], M_basis="dry")
