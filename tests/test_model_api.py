"""Tests for the high-level model API (src/model_api.py) used by the Streamlit app.

Head-less: exercises simulate() across congeners and scenario toggles, so the app
is covered without launching Streamlit."""
import os

import numpy as np
import pytest

import model_api as api

_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def test_congener_list_and_chain_table():
    assert len(api.CONGENERS) == 13          # 12 calibrated + GenX (ether-PFAS, provisional)
    assert "GenX" in api.CONGENERS
    rows = api.chain_table()
    assert len(rows) == 13
    assert all({"name", "n_C", "group", "K_PL", "f_xy_recommended"} <= set(r) for r in rows)


def test_simulate_all_congeners_run_and_finite():
    for nm in api.CONGENERS:
        r = api.simulate(nm)
        assert r["success"]
        assert set(r["baf_final"]) == set(api.TISSUES)
        assert all(np.isfinite(v) and v >= 0 for v in r["baf_final"].values())
        assert r["straw_baf"] >= 0 and np.isfinite(r["straw_baf"])
        assert all(b > 0 for b in r["B_k"].values())


def test_simulate_structure_and_units():
    r = api.simulate("PFOA", Cwo=1.0)
    # Cwo=1 -> tissue conc equals BAF at final time
    assert r["conc"]["root"][-1] == pytest.approx(r["baf_final"]["root"], rel=1e-6)
    assert r["t"][0] == 0.0 and r["t"][-1] == pytest.approx(120.0)


def test_membrane_potential_lever():
    # more negative E_m -> stronger anion exclusion (larger e^N)
    r_neg = api.simulate("PFOA", E_m_mV=-140.0)
    r_pos = api.simulate("PFOA", E_m_mV=-90.0)
    assert r_neg["eN"] > r_pos["eN"]


def test_forcing_and_fxy_toggles_run():
    a = api.simulate("PFOA", measured_forcing=True, f_xy_source="recommended")
    b = api.simulate("PFOA", measured_forcing=False, f_xy_source="W2fit")
    assert a["success"] and b["success"]
    # the two f_xy sources differ for a long chain (recommended monotone vs W2 fit)
    rec = api.simulate("PFDA", f_xy_source="recommended")["params"]["f_xy"]
    w2 = api.simulate("PFDA", f_xy_source="W2fit")["params"]["f_xy"]
    assert rec != w2


def test_observed_baf_lookup():
    o = api.observed_baf("PFOA")
    assert {"root", "straw", "grain"} <= set(o)
    assert api.observed_baf("not-a-pfas") == {}


def test_unknown_congener_raises():
    with pytest.raises(KeyError):
        api.simulate("not-a-pfas")


def test_result_carries_drivers_and_baf_series():
    r = api.simulate("PFOA")
    assert r["Cwo"].shape == r["t"].shape == r["Qtp"].shape
    assert r["M"].shape == (len(r["t"]), 4)
    assert set(r["baf"]) == set(api.TISSUES)
    assert r["season"] == pytest.approx(120.0)
    # baf series matches the final BAF at t[-1]
    assert r["baf"]["root"][-1] == pytest.approx(r["baf_final"]["root"], rel=1e-9)


def test_metric_series_and_schematic_values():
    r = api.simulate("PFOA")
    ms = api.metric_series(r, "conc")
    assert set(ms["data"]) == set(api.TISSUES) and ms["cmax"] > ms["cmin"]
    sv = api.schematic_values(r, "conc", -1)
    assert {"root", "stem", "leaf", "grain", "straw"} <= set(sv["values"])
    assert sv["Cwo"] > 0 and np.isfinite(sv["values"]["leaf"])
    # baf metric also builds
    assert np.isfinite(api.schematic_values(r, "baf", 0)["values"]["root"])


def test_drivers_override_runs_and_sets_season():
    t = np.linspace(0.0, 90.0, 91)
    Cwo = np.full_like(t, 2.0)
    drv = api.drivers_from_arrays(t, Cwo)               # Qtp/M default to measured
    assert drv["M"].shape == (91, 4)
    r = api.simulate("PFOA", drivers=drv)
    assert r["success"] and r["season"] == pytest.approx(90.0)
    assert r["Cwo"][0] == pytest.approx(2.0)


def test_pore_water_from_inventory_drained_and_flooded():
    t = np.linspace(0.0, 120.0, 121)
    Cwo_d, soil = api.pore_water_from_inventory(t, 5.0, K_F=2.0, n=0.85, theta_g=0.35)
    assert np.all(Cwo_d > 0) and hasattr(soil, "K_F")
    # flooded (higher water content) dilutes -> same inventory gives LOWER pore water
    flooded = np.ones_like(t, bool)
    Cwo_f, redox = api.pore_water_from_inventory(
        t, 5.0, K_F=2.0, n=0.85, theta_g=0.35, theta_g_flooded=0.60, flooded=flooded)
    assert Cwo_f[10] < Cwo_d[10]


def test_load_driver_csv_example():
    path = os.path.join(_ROOT, "examples", "hydrus_drivers_example.csv")
    drv = api.load_driver_csv(path)
    assert drv["t"][0] == 0.0 and drv["M"].shape[1] == 4
    r = api.simulate("PFOA", drivers=drv)
    assert r["success"]


def test_biomonitoring_baf_and_csv():
    baf = api.baf_from_measurement({"root": 0.49, "straw": 0.83, "grain": 0.46}, 1.0)
    assert baf["straw"] == pytest.approx(0.83)
    assert api.baf_from_measurement({"root": 1.0}, 0.0) == {}     # no pore water -> {}
    bio = api.load_biomonitoring_csv(os.path.join(_ROOT, "examples", "biomonitoring_example.csv"))
    assert {"root", "straw", "grain"} <= set(bio["conc"]) and bio["Cwo"] == pytest.approx(1.0)


def test_lipid_loading_off_matches_baseline():
    """lipid_loading=False must recover the free-only model exactly (g=0)."""
    for nm in ("PFBA", "PFOA", "PFDA"):
        off = api.simulate(nm, lipid_loading=False, measured_forcing=True)
        assert off["params"]["g_xy"] == 0.0 and off["params"]["g_ph"] == 0.0
        rec = api.simulate(nm, f_xy_source="recommended", measured_forcing=True)
        assert off["baf_final"]["grain"] == pytest.approx(rec["baf_final"]["grain"], rel=1e-9)


def test_lipid_loading_rescues_longchain_grain():
    """The opt-in lipid term must lift the long-chain shoot the free-only model
    structurally starves (PFDA grain ~0.04 -> O(1)+), staying finite/mass-sane."""
    off = api.simulate("PFDA", lipid_loading=False, measured_forcing=True)
    on = api.simulate("PFDA", lipid_loading=True, measured_forcing=True)
    assert on["success"]
    assert on["baf_final"]["grain"] > 20 * off["baf_final"]["grain"]   # 0.04 -> few
    assert all(np.isfinite(v) and v >= 0 for v in on["baf_final"].values())


def test_lipid_conductances_are_KPL_gated():
    """g_xy/g_ph ~0 for short chains, rise with K_PL for long chains; f_xy declines."""
    short = api.lipid_loading_conductances(4, 42.0, "PFCA")      # PFBA
    long = api.lipid_loading_conductances(11, 5.6e4, "PFCA")     # PFUnDA-like
    assert short[1] < 0.01 and short[2] < 0.01                  # bound loading off (short)
    assert long[1] > short[1] and long[2] > short[2]            # rises with K_PL
    assert long[0] < short[0]                                   # free TSCF declines
    # PFSA head-group factor lowers all three vs the CF2-matched PFCA
    assert api.lipid_loading_conductances(8, 1e4, "PFSA")[1] < api.lipid_loading_conductances(8, 1e4, "PFCA")[1]


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))
