"""Tests for measured_biomass: unit conversion, table load, M(t), drivers, simulate."""
import io
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))

import measured_biomass as mb          # noqa: E402
import model_api as api                # noqa: E402


def test_unit_conversions():
    v = np.array([1000.0])
    assert mb.to_kg_per_hill(v, "g/hill")[0] == pytest.approx(1.0)
    assert mb.to_kg_per_hill(v, "mg/hill")[0] == pytest.approx(1e-3)
    assert mb.to_kg_per_hill(v, "kg/hill")[0] == pytest.approx(1000.0)
    # g/plant with 3 plants/hill -> kg/hill
    assert mb.to_kg_per_hill(np.array([10.0]), "g/plant", plants_per_hill=3)[0] == pytest.approx(0.03)
    # area-based: 1 t/ha at 25 hills/m2 = 1000 kg/ha /1e4 /25 = 0.004 kg/hill
    assert mb.to_kg_per_hill(np.array([1.0]), "t/ha", hills_per_m2=25)[0] == pytest.approx(0.004)
    assert mb.to_kg_per_hill(np.array([100.0]), "g/m2", hills_per_m2=25)[0] == pytest.approx(0.004)


def test_g_plant_requires_density():
    with pytest.raises(ValueError):
        mb.to_kg_per_hill(np.array([1.0]), "g/plant")


def test_unknown_units():
    with pytest.raises(ValueError):
        mb.to_kg_per_hill(np.array([1.0]), "stones/acre")


def _demo_csv():
    return io.StringIO(
        "day,root,stem,leaf,grain\n"
        "0,0.05,0.02,0.03,0\n"
        "40,0.6,1.2,1.4,0\n"
        "80,1.0,5.5,4.0,3.0\n"
        "150,1.1,7.0,4.5,12.0\n")


def test_load_table_and_matrix():
    day, organs = mb.load_biomass_table(_demo_csv())
    assert day[0] == 0.0 and day[-1] == 150.0
    assert set(k for k, v in organs.items() if v is not None) == {"root", "stem", "leaf", "grain"}
    t = np.linspace(0, 150, 61)
    M = mb.biomass_matrix(t, day, {k: v for k, v in organs.items()})
    assert M.shape == (61, 4)
    # monotone-ish grain rises; final grain interpolated to the last sample
    assert M[-1, 3] == pytest.approx(12.0, rel=1e-6)
    assert np.all(M > 0)


def test_missing_root_reconstructed_from_ratio():
    csv = io.StringIO("day,stem,leaf,grain\n0,0.02,0.03,0\n150,7.0,4.5,12.0\n")
    day, organs = mb.load_biomass_table(csv)
    assert organs["root"] is None
    t = np.linspace(0, 150, 11)
    M = mb.biomass_matrix(t, day, {k: v for k, v in organs.items()}, root_shoot_ratio=0.1)
    shoot_final = M[-1, 1] + M[-1, 2] + M[-1, 3]
    assert M[-1, 0] == pytest.approx(0.1 * shoot_final, rel=1e-6)


def test_biomass_drivers_and_simulate():
    drv = mb.biomass_drivers(_demo_csv(), units="g/plant", plants_per_hill=3,
                             season=150.0, n_t=121)
    assert set(drv) == {"t", "Cwo", "Qtp", "M"}
    assert drv["M"].shape == (121, 4)
    assert drv["t"][-1] == pytest.approx(150.0)
    # the drivers must run through the canonical ODE and conserve the structure
    res = api.simulate("PFOA", drivers=drv)
    assert res["success"]
    assert all(np.isfinite(v) for v in res["baf_final"].values())
    assert res["t"][-1] == pytest.approx(150.0)
