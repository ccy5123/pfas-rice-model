"""Tests for the high-level model API (src/model_api.py) used by the Streamlit app.

Head-less: exercises simulate() across congeners and scenario toggles, so the app
is covered without launching Streamlit."""
import os

import numpy as np
import pytest

from conftest import SOLVE_INVARIANCE_REL

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
    assert r["t"][0] == 0.0 and r["t"][-1] == pytest.approx(api.site_season())


def test_membrane_potential_lever():
    # more negative E_m -> stronger anion exclusion (larger e^N)
    r_neg = api.simulate("PFOA", E_m_mV=-140.0)
    r_pos = api.simulate("PFOA", E_m_mV=-90.0)
    assert r_neg["eN"] > r_pos["eN"]


def test_forcing_and_fxy_toggles_run():
    a = api.simulate("PFOA", measured_forcing=True, f_xy_source="recommended")
    b = api.simulate("PFOA", measured_forcing=False, f_xy_source="fit")
    assert a["success"] and b["success"]
    # the two f_xy sources differ for a long chain (recommended monotone vs the fit)
    rec = api.simulate("PFDA", f_xy_source="recommended")["params"]["f_xy"]
    fit = api.simulate("PFDA", f_xy_source="fit")["params"]["f_xy"]
    assert rec != fit


def test_observed_baf_lookup():
    o = api.observed_baf("PFOA")
    assert {"root", "straw", "grain"} <= set(o)
    assert api.observed_baf("not-a-pfas") == {}
    assert set(api.observed_baf("PFBS")) == {"root", "straw"}        # imputed grain left out
    d = api.yamazaki_obs("PFOA")
    assert d["straw_lo"] < o["straw"] < d["straw_hi"] and d["water_ng_L"] > 0


def test_unknown_congener_raises():
    with pytest.raises(KeyError):
        api.simulate("not-a-pfas")


def test_result_carries_drivers_and_baf_series():
    r = api.simulate("PFOA")
    assert r["Cwo"].shape == r["t"].shape == r["Qtp"].shape
    assert r["M"].shape == (len(r["t"]), 4)
    assert set(r["baf"]) == set(api.TISSUES)
    assert r["season"] == pytest.approx(api.scenario_drivers()["t"][-1])
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


def test_drivers_from_arrays_respects_the_crop_selection():
    """drivers_from_arrays fills Qtp/M from the selected crop drivers: the site-year by
    default (its ORYZA leaf senesces and carries a death rate), the generic drivers with a
    biomass driver otherwise (growth_rice's leaf is monotonic and carries none)."""
    t = np.linspace(0.0, 110.0, 111)
    li = api.TISSUES.index("leaf")
    ds = api.drivers_from_arrays(t, np.ones_like(t))
    dg = api.drivers_from_arrays(t, np.ones_like(t), season=110.0, scenario="generic",
                                 biomass="growth_rice")
    Ms, Mg = ds["M"][:, li], dg["M"][:, li]
    site = api.scenario_drivers()
    assert np.allclose(ds["Qtp"], np.interp(t, site["t"], site["Qtp"]))
    assert "leaf_loss" in ds and Ms[-1] < 0.7 * Ms.max()    # the ORYZA leaf senesces
    assert "leaf_loss" not in dg and Mg[-1] >= 0.95 * Mg.max()  # growth_rice monotonic


def test_load_driver_csv_without_M_uses_the_crop_selection():
    """A driver CSV that omits the M_* columns falls back to the selected crop drivers --
    the CSV mode honours the sidebar crop choice like the soil-inventory mode."""
    import io
    rows = "\n".join(f"{d},1.0,0.05" for d in range(0, 111, 5))
    csv = "t,Cwo,Qtp\n" + rows
    li = api.TISSUES.index("leaf")
    ds = api.load_driver_csv(io.StringIO(csv))
    dg = api.load_driver_csv(io.StringIO(csv), scenario="generic", biomass="growth_rice")
    assert "leaf_loss" in ds and ds["M"][:, li][-1] < 0.7 * ds["M"][:, li].max()
    assert "leaf_loss" not in dg and dg["M"][:, li][-1] >= 0.95 * dg["M"][:, li].max()
    assert np.allclose(ds["Qtp"], 0.05)                       # a supplied Qtp is kept


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


def test_export_csv_helpers():
    """The UI download helpers produce well-formed CSV strings (pure, head-less)."""
    res = api.simulate("PFOA")
    obs = api.observed_baf("PFOA")

    summ = api.summary_csv(res, obs)
    lines = summ.strip().splitlines()
    assert lines[0] == "tissue,model_BAF_L_per_kg,final_conc_ug_per_kg,observed_BAF,measured_BAF"
    assert [ln.split(",")[0] for ln in lines[1:]] == ["root", "straw", "grain"]
    # each data row has exactly 5 fields and the model BAF parses back to the result
    root_row = lines[1].split(",")
    assert len(root_row) == 5
    assert float(root_row[1]) == pytest.approx(res["baf_final"]["root"], rel=1e-5)  # written as %.6g
    # measured column is blank without a biomonitoring dict; populated with one
    assert root_row[4] == ""
    summ_m = api.summary_csv(res, obs, {"root": 0.47}).splitlines()
    assert summ_m[1].split(",")[4] == "0.47"

    ts = api.timeseries_csv(res)
    head = ts.splitlines()[0].split(",")
    assert head[:3] == ["t", "Cwo", "Qtp"]
    assert "conc_grain" in head and "M_root" in head
    assert len(ts.strip().splitlines()) == len(res["t"]) + 1     # header + one row per time


def test_drivers_from_tables_growth_and_cwo():
    """User tables -> drivers: growth (fresh g/hill) sets M, Cwᵒ table is absolute and
    interpolated, and either table alone falls back to the built-in value."""
    growth = dict(day=[0, 40, 80, 120], root=[1, 80, 150, 160], stem=[1, 120, 300, 320],
                  leaf=[2, 140, 400, 300], grain=[0, 0, 200, 900])      # fresh g/hill
    cwo = dict(day=[0, 60, 120], Cwo=[2.0, 1.5, 0.5])
    drv = api.drivers_from_tables(growth, cwo, growth_units="g/hill")
    assert set(("t", "Cwo", "Qtp", "M")) <= set(drv)
    assert drv["t"][0] == 0.0 and drv["t"][-1] == pytest.approx(120.0)
    assert drv["Cwo"][0] == pytest.approx(2.0) and drv["Cwo"][-1] == pytest.approx(0.5)
    # g/hill -> kg/hill: final grain 900 g -> 0.9 kg
    assert drv["M"][-1, 3] == pytest.approx(0.9, rel=1e-6)
    # unsorted rows are handled (sorted by day internally)
    shuf = dict(day=[120, 0, 60], Cwo=[0.5, 2.0, 1.5])
    d2 = api.drivers_from_tables(growth, shuf, growth_units="g/hill")
    assert d2["Cwo"][0] == pytest.approx(2.0) and d2["Cwo"][-1] == pytest.approx(0.5)
    # Cwᵒ table only -> M falls back to the biomass driver (4 organ columns, leaf_loss attached)
    only_cwo = api.drivers_from_tables(None, cwo, Cwo_const=1.0)
    assert only_cwo["M"].shape[1] == 4 and "leaf_loss" in only_cwo
    # growth table only -> flat Cwo_const series
    only_g = api.drivers_from_tables(growth, None, growth_units="g/hill", Cwo_const=3.0)
    assert np.allclose(only_g["Cwo"], 3.0)
    # the drivers run through the ODE
    r = api.simulate("PFOA", drivers=drv)
    assert r["success"] and np.isfinite(r["straw_baf"])
    with pytest.raises(ValueError):
        api.drivers_from_tables(None, None)


def test_estimate_exposure_bayesian_recovers_and_brackets():
    """Synthetic recovery: generate tissue conc at a known Cwᵒ, then the Bayesian
    inverse should recover that level (within a few %) with the truth inside the
    95% credible interval, and reproduce the inputs at the MAP."""
    truth = 2.5
    r = api.simulate("PFOA", Cwo=truth)
    dw = r["conc_dw"]                                   # the inverse takes DRY-weight data
    meas = {"root": dw["root"][-1], "straw": dw["straw"][-1], "grain": dw["grain"][-1]}
    est = api.estimate_exposure_bayesian("PFOA", meas)
    assert est["median"] == pytest.approx(truth, rel=0.10)        # MAP recovers the level
    lo, hi = est["ci95"]
    assert lo < truth < hi                                        # truth inside 95% CI
    assert lo < est["median"] < hi
    assert set(est["used_tissues"]) == {"root", "straw", "grain"}
    # the model reproduces the entered measurements at the best estimate
    for t_ in est["used_tissues"]:
        assert est["model_fit"][t_] == pytest.approx(meas[t_], rel=0.20)
    # grid is monotone-increasing Cwᵒ with a normalised density peaking at 1
    g = np.asarray(est["grid"]["Cwo"])
    assert np.all(np.diff(g) > 0) and max(est["grid"]["density"]) == pytest.approx(1.0)
    # a single tissue still yields an estimate; empty input is rejected
    one = api.estimate_exposure_bayesian("PFOA", {"grain": meas["grain"]})
    assert one["median"] > 0 and one["used_tissues"] == ["grain"]
    with pytest.raises(ValueError):
        api.estimate_exposure_bayesian("PFOA", {"root": 0.0})


def test_estimate_exposure_bayesian_runs_under_the_selected_mechanism():
    """The inverse must use the SAME model as the forward run, not silently fall back
    to the shipped one: the app passes the sidebar's mechanism switches straight
    through. The bypass takes up more anion per unit exposure, so explaining the same
    measured tissue needs a LOWER Cwᵒ than the carrier does."""
    meas = {"grain": api.simulate("PFOA", Cwo=2.0)["conc_dw"]["grain"][-1]}
    carrier = api.estimate_exposure_bayesian("PFOA", meas, uptake="carrier")
    bypass = api.estimate_exposure_bayesian("PFOA", meas, uptake="bypass")
    assert carrier["median"] == pytest.approx(2.0, rel=0.10)       # default is unchanged
    assert bypass["median"] < carrier["median"]
    # lipid loading is likewise honoured (it moves the grain hard, so the level moves)
    lipid = api.estimate_exposure_bayesian("PFOA", meas, lipid_loading=True)
    assert lipid["median"] != pytest.approx(carrier["median"], rel=0.05)


def test_estimate_exposure_progress_callback():
    """The optional progress callback (for the app's live bar) fires once per ODE
    solve with monotone done<=total, and a throwing callback never breaks the fit."""
    calls = []
    r = api.simulate("PFOA", Cwo=1.0)
    meas = {"grain": r["conc"]["grain"][-1]}
    est = api.estimate_exposure_bayesian("PFOA", meas, progress=lambda d, t: calls.append((d, t)))
    assert calls, "progress callback was never called"
    dones = [d for d, _ in calls]
    totals = {t for _, t in calls}
    assert totals == {api._ESTIMATE_MAX_SOLVES}            # total is the documented bound
    assert dones == sorted(dones)                          # monotone non-decreasing
    assert max(dones) <= api._ESTIMATE_MAX_SOLVES          # never exceeds the bound
    assert len(calls) == est["n_eval"]                     # one tick per ODE solve
    # a callback that raises must be swallowed (the estimate still succeeds)
    def _boom(d, t):
        raise RuntimeError("ui blew up")
    est2 = api.estimate_exposure_bayesian("PFOA", meas, progress=_boom)
    assert est2["median"] == pytest.approx(est["median"], rel=1e-9)


def test_predictive_band_and_uncertainty_factor():
    """The general-audience honesty band: a ×/÷ factor from the documented a-priori
    predictive log10 RMSE (~7×), symmetric in log space, with graceful NaN edges."""
    f = api.uncertainty_factor()
    assert f == pytest.approx(10.0 ** api.APRIORI_LOG10_RMSE)
    assert api.APRIORI_LOG10_RMSE == pytest.approx(
        float(np.median(list(api.APRIORI_OOS_SCORES.values()))))
    assert f > 10.0                     # the default model is more than 10x off out of sample
    b = api.predictive_band(2.0)
    assert b["factor"] == pytest.approx(f)
    assert b["lo"] == pytest.approx(2.0 / f) and b["hi"] == pytest.approx(2.0 * f)
    assert b["lo"] < b["value"] < b["hi"]
    # geometric symmetry: value is the geomean of lo/hi
    assert np.sqrt(b["lo"] * b["hi"]) == pytest.approx(b["value"])
    # a custom (smaller) RMSE narrows the band; zero -> no band (factor 1)
    assert api.uncertainty_factor(0.30) < f
    assert api.uncertainty_factor(0.0) == pytest.approx(1.0)
    # non-finite / negative -> NaN lo/hi but a finite factor (no crash on edge values)
    for bad in (float("nan"), -1.0):
        e = api.predictive_band(bad)
        assert not np.isfinite(e["lo"]) and not np.isfinite(e["hi"])
        assert np.isfinite(e["factor"])
    # value 0 is a valid (degenerate) band, not an error
    z = api.predictive_band(0.0)
    assert z["lo"] == 0.0 and z["hi"] == 0.0


def test_intake_fraction_efsa_twi():
    """The EFSA-TWI intake reference (policy hook): a pure conversion of a grain
    concentration to a % of the group tolerable weekly intake, with the group
    membership flag and traffic-light signal, and graceful NaN handling."""
    # unit check: 0.4 µg/kg == 0.4 ng/g; 150 g/day, 60 kg -> 7 ng/kg bw/week
    # = 7/4.4 = ~159 % of the 4.4 ng/kg bw/week TWI.
    info = api.intake_fraction(0.4, congener="PFOA",
                               rice_intake_g_day=150.0, body_weight_kg=60.0)
    assert info["weekly_intake_ng_per_kg_bw"] == pytest.approx(0.4 * 150 * 7 / 60)
    assert info["twi"] == pytest.approx(api.EFSA_TWI_NG_PER_KG_BW_WEEK)
    assert info["percent"] == pytest.approx(100.0 * (0.4 * 150 * 7 / 60) / 4.4, rel=1e-9)
    assert info["in_group"] is True and info["signal"] == "red"
    # linearity in grain conc and rice intake; inverse in body weight
    a = api.intake_fraction(1.0, congener="PFOA")
    b = api.intake_fraction(2.0, congener="PFOA")
    assert b["percent"] == pytest.approx(2 * a["percent"])
    heavier = api.intake_fraction(1.0, congener="PFOA", body_weight_kg=120.0)
    assert heavier["percent"] == pytest.approx(a["percent"] * 60.0 / 120.0)
    # traffic-light thresholds (10 % / 100 %)
    assert api.intake_fraction(1e-4, congener="PFOA")["signal"] == "green"
    # the four group members flagged in, everything else out; None if unspecified
    for g in api.EFSA_TWI_CONGENERS:
        assert api.intake_fraction(1.0, congener=g)["in_group"] is True
    assert api.intake_fraction(1.0, congener="PFBA")["in_group"] is False
    assert api.intake_fraction(1.0)["in_group"] is None
    # NaN / negative / zero-bw -> NaN percent, signal None, no crash
    for bad in (float("nan"), -1.0):
        e = api.intake_fraction(bad, congener="PFOA")
        assert not np.isfinite(e["percent"]) and e["signal"] is None
    assert api.intake_fraction(1.0, congener="PFOA", body_weight_kg=0.0)["signal"] is None


def test_summary_report_md():
    """The one-page handout summary is a pure Markdown string carrying the tissue
    numbers (+band), the EFSA-TWI line, the caveats and the sources, in both langs."""
    res = api.simulate("PFOA", Cwo=1.0)
    info = api.intake_fraction(float(res["conc"]["grain"][-1]), congener="PFOA")
    md = api.summary_report_md(res, congener="PFOA", obs=api.observed_baf("PFOA"),
                               intake=info, scenario="중간", lang="ko")
    assert isinstance(md, str) and md.startswith("# ")
    assert "PFOA" in md and "EFSA" in md and "출처" in md
    assert "MRL" in md                                        # no-legal-limit caveat present
    assert f"{info['percent']:.0f}%" in md                    # the EFSA % line
    # English variant renders too, and flags an out-of-group congener as reference
    info_pfba = api.intake_fraction(1.0, congener="PFBA")
    en = api.summary_report_md(res, congener="PFBA", intake=info_pfba, lang="en")
    assert "Sources:" in en and "reference" in en             # PFBA not in the EFSA four
    # tolerates a missing intake (no EFSA line, still builds)
    assert api.summary_report_md(res, congener="PFOA").startswith("# ")


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
    """simulate_nstem_leaf runs on the scenario drivers by default; a custom biomass driver
    switches it to the generic drivers and gives a different leaf transfer factor."""
    import growth_rice as gr2
    default = api.simulate_nstem_leaf("PFOA", Cwo=1.0)
    growth = api.simulate_nstem_leaf("PFOA", Cwo=1.0, biomass_fn=gr2.organ_biomass)
    assert default["success"] and growth["success"]
    assert default["season"] == pytest.approx(api.scenario_drivers()["t"][-1])
    assert abs(default["tf_final"]["leaf"] - growth["tf_final"]["leaf"]) > 1e-3
    assert all(np.isfinite(v) for v in default["tf_final"].values())


def test_tang_tf_validation_and_observed():
    """Tang per-organ TF (dry weight, harvest operator) for the 3 Tang congeners."""
    assert api.tang_tf_validation("PFNA") is None
    assert api.tang_observed_tf("PFNA") == {}
    for c in api.TANG_CONGENERS:
        v = api.tang_tf_validation(c)
        assert set(v["organs"]) == {"stalk", "leaf", "endosperm"}
        assert set(v["model_tf"]) == set(v["organs"])
        assert v["tang_tf"] and all(t > 0 for t in v["tang_tf"].values())
        assert all(np.isfinite(m) and m >= 0 for m in v["model_tf"].values())
        assert api.tang_tf_validation(c, use_refit=True)["f_xy"] == api.TANG_REFIT_FXY[c]
    # Tang TF declines with dose -> 0.1 ug/g (low) value exceeds the across-dose mean
    low, mean = api.tang_observed_tf("PFOA", "low"), api.tang_observed_tf("PFOA", "mean")
    assert low["stalk"] > mean["stalk"]


def test_scenarios():
    """The default drivers are the calibration site-year; the generic drivers take a
    biomass driver and a season; an unknown scenario or biomass driver raises."""
    base = api.simulate("PFOA")
    assert base["params"]["scenario"] == "yamazaki2020"
    assert base["season"] == pytest.approx(api.scenario_drivers()["t"][-1])
    gr_ = api.simulate("PFOA", scenario="generic", biomass="growth_rice")
    oz = api.simulate("PFOA", scenario="generic", biomass="oryza2000")
    assert gr_["season"] == 120.0 and not np.allclose(gr_["M"][-1], oz["M"][-1])
    assert all(np.isfinite(v) and v >= 0 for v in gr_["baf_final"].values())
    with pytest.raises(ValueError):
        api.simulate("PFOA", scenario="tsukuba1999")
    with pytest.raises(ValueError):
        api.simulate("PFOA", scenario="generic", biomass="sucros")


def test_fit_source_reproduces_yamazaki():
    """f_xy_source='fit' (validation/refit_site_year.py) reproduces the Yamazaki dry-basis
    observation for a congener that sits off the bounds; 'recommended' keeps the monotone
    f_xy and reproduces root and grain only -- and only where L_Ph is off its bounds:
    from PFOA (C8) up the monotone f_xy leaves the grain short even at L_Ph = 1, so the
    'recommended' check uses PFHpA (off every bound in both fits)."""
    o = api.observed_baf("PFOA")
    fit = api.yamazaki_harvest_baf("PFOA", f_xy_source="fit")
    for k in ("root", "straw", "grain"):
        assert abs(np.log10(fit[k] / o[k])) < 0.01
    o7 = api.observed_baf("PFHpA")
    rec = api.yamazaki_harvest_baf("PFHpA", f_xy_source="recommended")
    assert abs(np.log10(rec["root"] / o7["root"])) < 0.01
    assert abs(np.log10(rec["grain"] / o7["grain"])) < 0.01
    short = api.yamazaki_harvest_baf("PFOA", f_xy_source="recommended")
    assert short["grain"] < o["grain"]          # at the L_Ph bound, grain under
    p = {s: api._transport_defaults(api._CONG["PFOA"], s, False) for s in api.F_XY_SOURCES}
    assert p["recommended"][0] == api._CONG["PFOA"]["f_xy_recommended"]
    assert p["fit"][0] == api._CONG["PFOA"]["f_xy_fit"]
    with pytest.raises(ValueError):
        api._transport_defaults(api._CONG["PFOA"], "W2fit", False)


def test_leaf_senescence_loss():
    """The ORYZA leaf death rate rides on the drivers; the leaf ODE loses -leaf_loss*C (the
    dead blade takes its PFAS) and the harvest operator keeps those blades in leaf/straw."""
    drv = api.scenario_drivers()
    assert float(np.max(drv["leaf_loss"])) > 0.0
    fix = api.simulate("PFOA", drivers=drv)
    nofix = api.simulate("PFOA", drivers={k: v for k, v in drv.items() if k != "leaf_loss"})
    assert fix["baf_final"]["leaf"] < nofix["baf_final"]["leaf"]
    assert fix["dead_leaf"]["burden"][-1] > 0 and nofix["dead_leaf"]["burden"][-1] == 0


def test_grain_formation_gate():
    """The grain takes NO PFAS before the panicle forms (no pre-flowering spike), then
    rises to its harvest value; a constant-mass driver still loads the grain (gate=1)."""
    r = api.simulate("PFOA")
    run = api.site_run()
    t_head = float(np.interp(1.0, run["dvs"], run["t"]))
    t, cg = r["t"], np.asarray(r["conc"]["grain"])
    assert np.all(cg[t < t_head - 25] < 1e-3)         # ~0 well before heading
    assert cg[-1] > 1e-3                              # accumulates by harvest (terminal sink)
    assert cg.max() <= cg[-1] * 1.5                   # NO pre-formation spike above harvest
    t2 = np.linspace(0, 120, 121)
    Mc = np.tile([0.003, 0.02, 0.013, 0.03], (121, 1))
    rc = api.simulate("PFOA", drivers=api.drivers_from_arrays(t2, np.ones(121), M=Mc))
    assert rc["baf_final"]["grain"] > 0.0


def test_simulate_nstem_leaf_k_seq_off_by_default_and_drivers_path():
    """The optional sequestered root pool is inert unless asked for, and the `drivers=`
    entry point reproduces the built-in scenario it mirrors."""
    base = api.simulate_nstem_leaf("PFOA")
    assert base["params"]["k_seq"] == 0.0 and base["seq_fraction"] == 0.0
    assert "root_seq" not in base["conc"]
    drv = api.scenario_drivers()
    same = api.simulate_nstem_leaf("PFOA", drivers=drv)
    for k in ("root", "stem", "leaf", "grain"):
        assert same["baf_final"][k] == pytest.approx(base["baf_final"][k], rel=1e-6)


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


def test_default_k_leach_is_hydrus_calibrated():
    """flooded defaults k_leach to a per-congener value CALIBRATED to HYDRUS-1D:
    short chains leach (larger k_leach), long chains are buffered (~0). A novel/SMILES
    compound (no table entry) falls back to k_leach(log10 Koc), clipped to [0, max]."""
    short, long = api.default_k_leach("PFBA"), api.default_k_leach("PFDoDA")
    assert short > 0.02 and long == pytest.approx(0.0, abs=1e-9) and short > long
    # PFBA's flooded decline with the auto (calibrated) default leaches harder than the old flat 0.02
    auto = api.simulate("PFBA", Cwo=1.0, cwo_profile="flooded")["Cwo"]
    flat = api.simulate("PFBA", Cwo=1.0, cwo_profile="flooded", cwo_kw=dict(k_leach=0.02))["Cwo"]
    assert auto[-1] / auto[0] < flat[-1] / flat[0]
    nov = api.default_k_leach(None, n_C=8, group="PFCA")          # novel -> Koc fallback
    assert 0.0 <= nov <= api._KLEACH_MAX


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
    t = api.scenario_drivers()["t"]                       # the site-year grid (its HYDRUS run)
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


# ---------------------------------------------------------------------------
# NEUTRAL-organic path exposed through the API (handoff A2)
# ---------------------------------------------------------------------------
def test_simulate_neutral_matches_the_standalone_module():
    """The wrapper must be a pass-through, not a second implementation.

    `model_api.simulate_neutral` and `neutral_dpu.simulate_neutral` have to agree
    BIT FOR BIT on the same forcings -- otherwise the published a-priori numbers
    (Liu 0.281, Ge 0.783) would silently stop describing what the API returns.
    Both entry points are checked: the built-in forcings and an explicit
    `drivers=` dict.
    """
    import neutral_dpu as ND

    lk, season, n_t = 2.45, 120.0, 241
    t = np.linspace(0.0, season, n_t)
    Qtp, M = api.measured_forcing(t, season, biomass="growth_rice", scenario="generic")
    drv = dict(t=t, Cwo=np.ones(n_t), Qtp=Qtp, M=M)
    ref = ND.simulate_neutral(ND.NeutralCompound("probe", lk), drv)

    built_in = api.simulate_neutral(lk, name="probe", biomass="growth_rice",
                                    season=season, n_t=n_t, scenario="generic")
    explicit = api.simulate_neutral(lk, name="probe", drivers=drv)
    for k in api.TISSUES:
        assert np.array_equal(built_in["conc"][k], ref["conc"][k])
        assert np.array_equal(explicit["conc"][k], ref["conc"][k])


def test_simulate_neutral_returns_the_simulate_contract():
    """Same result-dict shape as `simulate()`, so plots/validation can treat the
    two interchangeably -- plus the neutral diagnostics that have no PFAS analogue."""
    import neutral_dpu as ND

    r = api.simulate_neutral(2.45, name="carbamazepine", half_life=7.0)
    pfas = api.simulate("PFOA")
    for k in ("t", "conc", "baf", "baf_final", "straw", "straw_baf", "cwo_ref",
              "season", "M", "Cwo", "Qtp", "params", "success"):
        assert k in r and k in pfas
    assert set(r["conc"]) >= set(api.TISSUES)
    assert set(r["tf_final"]) >= set(api.TISSUES)      # neutral/nstem_leaf endpoint
    # nothing fitted: partitioning and loading both follow from log Kow
    assert r["TSCF"] == pytest.approx(ND.briggs_tscf(2.45))
    assert r["K_PW"]["root"] > 0 and r["rcf_briggs"] > 0
    assert r["params"]["gamma"] == pytest.approx(np.log(2.0) / 7.0)
    # z = 0: the ionic machinery is off (PFAS: N = 4.67, e^N ~ 107)
    assert r["N"] == pytest.approx(0.0, abs=1e-12)
    assert r["eN"] == pytest.approx(1.0, abs=1e-12)
    assert pfas["N"] > 4.0


def test_simulate_neutral_air_is_opt_in_and_zero_for_pfas_like_compounds():
    """Air exchange must stay opt-in through the API too, and must be identically
    zero at K_AW = 0 -- the property that keeps every PFAS number untouched."""
    kw = dict(name="probe", MW=131.4, biomass="growth_rice", season=120.0, n_t=241,
              scenario="generic")
    off = api.simulate_neutral(2.42, K_AW=0.0, **kw)
    on_zero = api.simulate_neutral(2.42, K_AW=0.0, air=True, **kw)
    volatile = api.simulate_neutral(2.42, K_AW=0.1, air=True, **kw)
    assert off["air"] is False and on_zero["air"] is True
    for k in api.TISSUES:
        assert np.array_equal(off["conc"][k], on_zero["conc"][k])
    # a genuinely volatile compound loses the leaf; the root never volatilises
    assert volatile["baf_final"]["leaf"] < 1e-3 * off["baf_final"]["leaf"]
    # Two DISTINCT solves agreeing by physics, not by construction -- case 3 in
    # tests/conftest.py. The exact np.array_equal above is the strong guard here;
    # this one can only hold to the solver's accuracy.
    assert volatile["baf_final"]["root"] == pytest.approx(
        off["baf_final"]["root"], rel=SOLVE_INVARIANCE_REL)
