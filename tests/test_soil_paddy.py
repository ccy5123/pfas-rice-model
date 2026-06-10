"""Tests for the paddy soil Freundlich sub-model and input adapters."""
import numpy as np
import pytest

from soil_paddy import (
    FreundlichSoil, PaddyRedox, inputs_from_soil, load_inputs_csv,
    example_paddy_redox,
)
from pfas_rice_plant_module import PlantInputs


def test_freundlich_roundtrip_nonlinear():
    """pore_water(total(Cw)) must recover Cw for a non-linear isotherm."""
    soil = FreundlichSoil(K_F=3.0, n=0.8, theta_g=0.4)
    for Cw in (0.0, 0.01, 1.0, 7.5, 100.0):
        assert soil.pore_water(float(soil.total(Cw))) == pytest.approx(Cw, rel=1e-6, abs=1e-9)


def test_freundlich_linear_closed_form():
    """n=1 reduces to a linear K_d: Cw = C_T/(K_F+theta_g)."""
    soil = FreundlichSoil(K_F=2.0, n=1.0, theta_g=0.5)
    C_T = 10.0
    assert soil.pore_water(C_T) == pytest.approx(C_T / (2.0 + 0.5))


def test_pore_water_monotonic_and_zero():
    soil = FreundlichSoil(K_F=1.5, n=0.9, theta_g=0.3)
    assert soil.pore_water(0.0) == 0.0
    cw = [soil.pore_water(ct) for ct in (1.0, 5.0, 20.0, 100.0)]
    assert all(b > a for a, b in zip(cw, cw[1:]))            # strictly increasing


def test_nonlinear_sorption_differs_from_linear():
    """A genuine Freundlich (n<1) must not behave like a linear K_d."""
    nonlin = FreundlichSoil(K_F=3.0, n=0.7, theta_g=0.4)
    # secant Kd decreases with concentration when n<1
    assert nonlin.Kd_eff(0.1) > nonlin.Kd_eff(10.0)


def test_redox_flooding_increases_bioavailability():
    """Default paddy: flooding weakens sorption -> higher pore-water conc."""
    redox = example_paddy_redox(K_F_drained=2.0, K_F_flooded=1.0, n=0.85)
    C_T = 5.0
    cw_flooded = redox.flooded.pore_water(C_T)
    cw_drained = redox.drained.pore_water(C_T)
    assert cw_flooded > cw_drained


def test_redox_series_selects_per_timestep():
    redox = example_paddy_redox()
    C_total = np.array([5.0, 5.0, 5.0])
    flooded = np.array([True, False, True])
    cw = redox.pore_water_series(C_total, flooded)
    assert cw[0] == pytest.approx(redox.flooded.pore_water(5.0))
    assert cw[1] == pytest.approx(redox.drained.pore_water(5.0))
    with pytest.raises(ValueError):
        redox.pore_water_series(np.array([1.0, 2.0]), np.array([True]))


def test_inputs_from_soil_builds_plantinputs():
    t = np.linspace(0.0, 10.0, 11)
    C_total = np.full_like(t, 5.0)
    Qtp = np.full_like(t, 0.1)
    M = np.column_stack([np.full_like(t, 0.01)] * 4)
    redox = example_paddy_redox()
    flooded = t < 5.0
    inp = inputs_from_soil(t, C_total, Qtp, M, redox, flooded=flooded)
    assert isinstance(inp, PlantInputs)
    # flooded part should be more bioavailable than drained part
    assert inp.Cwo_(2.0) > inp.Cwo_(8.0)


def test_inputs_from_soil_requires_flood_schedule_for_redox():
    t = np.linspace(0.0, 1.0, 3)
    with pytest.raises(ValueError):
        inputs_from_soil(t, np.ones_like(t), np.ones_like(t),
                         np.ones((3, 4)), example_paddy_redox())


def test_load_inputs_csv_roundtrip(tmp_path):
    p = tmp_path / "soil.csv"
    p.write_text(
        "t,Cwo,Qtp,M_root,M_stem,M_leaf,M_grain\n"
        "0,1.0,0.05,0.001,0.001,0.001,0.0\n"
        "60,1.2,0.30,0.02,0.03,0.04,0.01\n"
        "120,0.8,0.10,0.03,0.04,0.05,0.025\n"
    )
    inp = load_inputs_csv(str(p))
    assert isinstance(inp, PlantInputs)
    assert inp.Cwo_(0.0) == pytest.approx(1.0)
    assert inp.M_(120.0)[3] == pytest.approx(0.025, rel=1e-6)


def test_load_inputs_csv_missing_column(tmp_path):
    p = tmp_path / "bad.csv"
    p.write_text("t,Cwo,Qtp\n0,1,0.1\n")
    with pytest.raises(ValueError):
        load_inputs_csv(str(p))
