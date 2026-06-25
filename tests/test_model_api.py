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


def test_oryza_refit_reproduces():
    """f_xy_source='oryza' (per-congener transport re-fit on the ORYZA2000 biomass,
    validation/refit_oryza.py) reproduces the Yamazaki tissue BAFs under the new
    default biomass -- the W2 fit (placeholder driver) does not."""
    o = api.observed_baf("PFOA")
    r = api.simulate("PFOA", f_xy_source="oryza", biomass="oryza")
    # PFOA reproduces to well within a factor of ~1.4 in log space on the ORYZA refit
    assert abs(np.log10(r["baf_final"]["root"]) - np.log10(o["root"])) < 0.15
    assert abs(np.log10(r["straw_baf"]) - np.log10(o["straw"])) < 0.15
    assert api._CONG["PFOA"]["f_xy_oryza"] > 0          # refit fields are present
    # the three sources resolve to different f_xy (recommended / W2 / oryza)
    fxy = {s: api.simulate("PFOA", f_xy_source=s)["params"]["f_xy"]
           for s in ("recommended", "W2fit", "oryza")}
    assert len({round(v, 4) for v in fxy.values()}) == 3


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


# --- sequestration two-pool root model (EXPLORATORY, opt-in) --------------------
def test_simulate_twopool_seq_structure_and_keys():
    """Opt-in sequestration two-pool run returns the SAME tissue keys as simulate(), plus the
    root mobile/seq split; root BAF == mobile + seq; finite and non-negative."""
    r = api.simulate_twopool_seq("PFOA")
    assert r["success"]
    assert set(api.TISSUES) <= set(r["baf_final"])
    assert all(np.isfinite(v) and v >= 0 for v in r["baf_final"].values())
    # reported root = mobile + sequestered pool
    assert r["baf_final"]["root"] == pytest.approx(
        r["baf_final"]["root_mobile"] + r["baf_final"]["root_seq"], rel=1e-6)
    assert 0.0 <= r["seq_fraction"] <= 1.0
    assert r["params"]["k_seq"] > 0 and r["params"]["k_rel"] == 0.0
    with pytest.raises(KeyError):
        api.simulate_twopool_seq("not-a-pfas")


def test_simulate_twopool_seq_matches_validation_and_rmse():
    """Drift guard: the model_api wrapper reproduces the standalone validation
    endpoints (within the ~1% driver-grid difference) AND the documented in-sample
    headline (overall log10 RMSE ~0.251, root ~0.156, PFOS/PFUnDA k_seq ~3x at
    identical K_PL) -- all with the MONOTONE physical f_xy_recommended."""
    (p, q), TP = api._twopool_seq()
    sq, err = [], {"root": [], "straw": [], "grain": []}
    for c in TP.CONGENERS:
        ks = TP.kseq_ushape(c["n_C"], c["group"], q)
        vr, vs, vg = TP.simulate(c, p, kseq_override=ks)         # standalone (endpoint)
        r = api.simulate_twopool_seq(c["name"])                      # wrapper (full series)
        ap = {"root": r["baf_final"]["root"], "straw": r["straw_baf"],
              "grain": r["baf_final"]["grain"]}
        for a, b in ((vr, ap["root"]), (vs, ap["straw"]), (vg, ap["grain"])):
            sq.append((np.log10(max(a, 1e-6)) - np.log10(max(b, 1e-6))) ** 2)
        o = api.observed_baf(c["name"])
        for k in err:
            if k in o:
                err[k].append((np.log10(max(ap[k], 1e-6)) - np.log10(o[k])) ** 2)
    assert np.sqrt(np.mean(sq)) < 0.05                            # cross-impl consistency
    overall = np.sqrt(np.mean(err["root"] + err["straw"] + err["grain"]))
    assert overall == pytest.approx(0.251, abs=0.02)              # documented headline
    assert np.sqrt(np.mean(err["root"])) == pytest.approx(0.156, abs=0.02)
    # monotone PHYSICAL f_xy (not a fitted U-shape) + non-K_PL PFOS/PFUnDA separation
    pfos, pfunda = api.simulate_twopool_seq("PFOS"), api.simulate_twopool_seq("PFUnDA")
    assert pfos["params"]["f_xy"] == api._CONG["PFOS"]["f_xy_recommended"]
    assert pfunda["params"]["k_seq"] / pfos["params"]["k_seq"] == pytest.approx(3.1, abs=0.4)


def test_simulate_twopool_seq_krel_drains_root_to_shoot():
    """A slow seq->mobile desorption (k_rel>0) releases the long-chain root burden:
    the root BAF falls (the irreversible sink stops retaining). Documented Result 5."""
    base = api.simulate_twopool_seq("PFUnDA", k_rel=0.0)
    rel = api.simulate_twopool_seq("PFUnDA", k_rel=0.2)
    assert rel["baf_final"]["root"] < base["baf_final"]["root"]
    assert rel["seq_fraction"] < base["seq_fraction"]


# ---------------------------------------------------------------------------
# time-varying pore-water exposure shape (cwo_profile)
# ---------------------------------------------------------------------------
def test_cwo_profile_constant_is_default():
    """The default (cwo_profile='constant') holds Cwo flat == the BAF-reproduction
    convention; passing it explicitly must reproduce the default run exactly."""
    a = api.simulate("PFOA", Cwo=1.0)
    b = api.simulate("PFOA", Cwo=1.0, cwo_profile="constant")
    assert np.allclose(a["Cwo"], 1.0)                         # flat
    assert a["Cwo"].std() == 0.0
    assert b["baf_final"]["grain"] == pytest.approx(a["baf_final"]["grain"], rel=1e-9)


def test_cwo_profile_flooded_shape_no_engine():
    """The analytic 'flooded' shape (Freundlich dilution+leaching) needs NO HYDRUS
    engine and is congener-resolved: a SHORT chain (low K_F -> large dissolved
    fraction) leaches to a steep decline, a LONG chain stays buffered (~flat). The
    season-mean is preserved (== level), so Cwo stays the average exposure."""
    short = api.simulate("PFBA", Cwo=2.0, cwo_profile="flooded")["Cwo"]
    long = api.simulate("PFDoDA", Cwo=2.0, cwo_profile="flooded")["Cwo"]
    # mean preserved == level for both
    assert np.mean(short) == pytest.approx(2.0, rel=1e-6)
    assert np.mean(long) == pytest.approx(2.0, rel=1e-6)
    # short chain declines (leaches), long chain is essentially flat (buffered)
    assert short[-1] < 0.5 * short[0]                         # steep decline
    assert long[-1] > 0.95 * long[0]                          # buffered
    assert short.std() > long.std()                          # short more time-variable


def test_cwo_profile_flooded_leach_rate_steepens_decline():
    """A larger leaching rate makes the short-chain pore water decline faster
    (end/start ratio smaller); a knob exposed via cwo_kw."""
    gentle = api.simulate("PFBA", Cwo=1.0, cwo_profile="flooded",
                          cwo_kw=dict(k_leach=0.01))["Cwo"]
    harsh = api.simulate("PFBA", Cwo=1.0, cwo_profile="flooded",
                         cwo_kw=dict(k_leach=0.08))["Cwo"]
    assert harsh[-1] / harsh[0] < gentle[-1] / gentle[0]


@pytest.mark.skipif(not api.hydrus_available(),
                    reason="HYDRUS-1D engine/phydrus not available")
def test_cwo_profile_hydrus_shape():
    """cwo_profile='hydrus' uses a real HYDRUS-1D run's Cwᵒ(t) shape (engine needed),
    normalised to season-mean == level; short chains leach (time-variable)."""
    r = api.simulate("PFBA", Cwo=1.0, cwo_profile="hydrus")
    assert np.mean(r["Cwo"]) == pytest.approx(1.0, rel=0.05)
    assert r["Cwo"].std() > 0.1                               # genuinely time-varying
    assert np.all(np.isfinite(r["conc"]["grain"]))


@pytest.mark.skipif(not api.hydrus_available(),
                    reason="HYDRUS-1D engine/phydrus not available")
def test_cwo_profile_flooded_matches_hydrus_direction():
    """The engine-free analytic 'flooded' shape reproduces the HYDRUS-1D DIRECTION:
    a short chain leaches (declines) under BOTH, a long chain stays buffered (~flat)
    under BOTH. (Quantitative k_leach calibration: validation/cwo_profile_check.py.)"""
    t = np.linspace(0.0, 120.0, 121)
    for cong, declines in (("PFBA", True), ("PFDoDA", False)):
        c = api._CONG[cong]
        hyd = api.cwo_profile_series(t, 1.0, "hydrus", n_C=c["n_C"], group=c["group"],
                                     congener=cong)
        flo = api.cwo_profile_series(t, 1.0, "flooded", n_C=c["n_C"], group=c["group"],
                                     congener=cong)
        if declines:                                         # short chain: both fall
            assert hyd[-1] < 0.5 * hyd[0] and flo[-1] < 0.7 * flo[0]
        else:                                                # long chain: both ~flat
            assert hyd[-1] > 0.9 * hyd[0] and flo[-1] > 0.95 * flo[0]
