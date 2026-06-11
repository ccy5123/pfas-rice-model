"""Tests for the high-level model API (src/model_api.py) used by the Streamlit app.

Head-less: exercises simulate() across congeners and scenario toggles, so the app
is covered without launching Streamlit."""
import numpy as np
import pytest

import model_api as api


def test_congener_list_and_chain_table():
    assert len(api.CONGENERS) == 12
    rows = api.chain_table()
    assert len(rows) == 12
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


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))


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
