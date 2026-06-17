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


def test_nstem_leaf_biomass_fn_override():
    """simulate_nstem_leaf accepts a custom biomass driver; the two drivers
    (growth_rice vs ORYZA) give different leaf transfer factors. Driver-explicit on
    both sides so it is robust to the default (now ORYZA2000)."""
    import oryza_growth as og
    import growth_rice as gr2
    ob = lambda t, s: og.organ_biomass_oryza(t, p=og.OryzaParams(season=s))
    growth = api.simulate_nstem_leaf("PFOA", Cwo=1.0, biomass_fn=gr2.organ_biomass)
    oryza = api.simulate_nstem_leaf("PFOA", Cwo=1.0, biomass_fn=ob)
    assert oryza["success"] and growth["success"]
    assert set(oryza["tf_final"]) == set(growth["tf_final"])
    # different biomass trajectory -> different leaf transfer factor
    assert abs(oryza["tf_final"]["leaf"] - growth["tf_final"]["leaf"]) > 1e-3
    assert all(np.isfinite(v) for v in oryza["tf_final"].values())
    # the no-arg default now resolves to ORYZA2000
    assert api.simulate_nstem_leaf("PFOA", Cwo=1.0)["tf_final"]["leaf"] == pytest.approx(
        oryza["tf_final"]["leaf"], rel=1e-9)


def test_tang_tf_validation_and_observed():
    """Tang per-organ TF (dry weight) for the 3 Tang congeners; None for others."""
    assert api.tang_tf_validation("PFNA") is None
    assert api.tang_observed_tf("PFNA") == {}
    for c in api.TANG_CONGENERS:
        v = api.tang_tf_validation(c)
        assert set(v["organs"]) == {"stalk", "leaf", "endosperm"}
        assert set(v["model_tf"]) == set(v["organs"])
        assert v["tang_tf"] and all(t > 0 for t in v["tang_tf"].values())
        assert all(np.isfinite(m) and m >= 0 for m in v["model_tf"].values())
        assert api.tang_tf_validation(c, use_refit=True)["f_xy"] == api.TANG_REFIT_FXY[c]
    # GenX refit LOWERS f_xy (0.233 -> 0.017, the documented ~12x over-prediction)
    assert (api.tang_tf_validation("GenX", use_refit=True)["f_xy"]
            < api.tang_tf_validation("GenX")["f_xy"])
    # Tang TF declines with dose -> 0.1 ug/g (low) value exceeds the across-dose mean
    low, mean = api.tang_observed_tf("PFOA", "low"), api.tang_observed_tf("PFOA", "mean")
    assert low["stalk"] > mean["stalk"]


def test_simulate_biomass_driver():
    """Biomass driver is selectable; DEFAULT == oryza (mechanistic ORYZA2000);
    growth_rice differs and both stay finite."""
    base = api.simulate("PFOA")
    gr_ = api.simulate("PFOA", biomass="growth_rice")
    oz = api.simulate("PFOA", biomass="oryza")
    assert base["baf_final"] == oz["baf_final"]                  # default is ORYZA2000
    assert all(np.isfinite(v) and v >= 0 for v in gr_["baf_final"].values())
    assert not np.allclose(np.asarray(gr_["M"][-1]), np.asarray(oz["M"][-1]))   # different M(t)
    assert oz["baf_final"]["grain"] != gr_["baf_final"]["grain"]


def test_oryza_leaf_senescence_loss():
    """ORYZA exposes a leaf death rate; the leaf senescence-loss flux (-leaf_loss*C) removes
    the spurious senescence concentration. growth_rice (no senescence) is unaffected."""
    import oryza_growth as og
    drv = og.oryza_drivers("PFOA", Cwo=1.0, season=120.0, p=og.OryzaParams(season=120.0))
    assert "leaf_loss" in drv and float(np.max(drv["leaf_loss"])) > 0.0   # rate exposed
    leaf_fix = api.simulate("PFOA", drivers=drv)["baf_final"]["leaf"]
    leaf_nofix = api.simulate(
        "PFOA", drivers={k: v for k, v in drv.items() if k != "leaf_loss"})["baf_final"]["leaf"]
    assert leaf_fix < leaf_nofix                              # senescence loss lowers the leaf
    # the default path is now ORYZA2000 (carries the senescence loss); growth_rice differs
    assert api.simulate("PFOA")["baf_final"] == api.simulate("PFOA", biomass="oryza")["baf_final"]


def test_grain_formation_gate():
    """The grain takes NO PFAS before the panicle forms (no pre-flowering spike), then
    rises to its harvest value; a constant-mass driver still loads the grain (gate=1)."""
    r = api.simulate("PFOA", biomass="oryza", season=120.0)
    t, cg = r["t"], np.asarray(r["conc"]["grain"])
    assert np.all(cg[t < 45] < 1e-3)                 # ~0 well before flowering (~d66)
    assert cg[-1] > 1e-3                              # accumulates by harvest (terminal sink)
    assert cg.max() <= cg[-1] * 1.5                   # NO pre-formation spike above harvest
    # constant-mass driver (e.g. HYDRUS CSV): grain must still load
    t2 = np.linspace(0, 120, 121)
    Mc = np.tile([0.003, 0.02, 0.013, 0.03], (121, 1))
    rc = api.simulate("PFOA", drivers=api.drivers_from_arrays(t2, np.ones(121), M=Mc))
    assert rc["baf_final"]["grain"] > 0.0
