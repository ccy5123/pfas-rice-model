"""Tests for the real HYDRUS-1D soil coupling (src/soil_hydrus.py).

The HYDRUS-invoking tests skip automatically when the compiled executable or
phydrus is unavailable (e.g. CI without gfortran), so the suite stays green
everywhere; build the engine to exercise them (see src/soil_hydrus.py docstring).
"""
import numpy as np
import pytest

import soil_hydrus as sh

pytest.importorskip("phydrus")
hydrus = pytest.mark.skipif(not sh.hydrus_available(),
                            reason="HYDRUS-1D executable not built")


# ---- pure parameter logic (no engine needed) -------------------------------
def test_paddy_kd_increases_with_chain_length():
    kd = [sh.paddy_kd(n, "PFCA") for n in (4, 6, 8, 10, 12)]
    assert all(b > a for a, b in zip(kd, kd[1:]))           # strictly increasing
    assert kd[0] < 0.1 < kd[-1]                             # spans >3 orders


def test_paddy_kd_sulfonate_above_matched_carboxylate():
    # PFSA sorbs more than the CF2-matched PFCA (same perfluoro-C count)
    assert sh.paddy_kd(8, "PFSA") > sh.paddy_kd(9, "PFCA")  # PFOS(8 CF) vs PFNA(8 CF)


def test_soil_inputs_are_the_sourced_midpoints():
    """PREREG 9.2: every scenario input is a sourced value or the midpoint of a sourced range."""
    soil = sh.andosol_soil()
    for k, (lo, hi) in sh.SOIL_RANGES.items():
        assert soil[k] == pytest.approx(0.5 * (lo + hi))
    a, b = sh.ANDOSOL_SETS["ash1"], sh.ANDOSOL_SETS["ash2"]
    assert soil["ths"] == pytest.approx(0.5 * (a["ths"] + b["ths"]))
    assert sh.FOC == pytest.approx(0.0375) and sh.RHO_PLOUGH == 0.61 and sh.PLOUGH_DEPTH == 15.0


def test_pan_ks_reproduces_the_sourced_percolation():
    """PREREG 9.8: the pan's Ks makes the steady flux equal the sourced percolation; integrating
    the steady head through the pan with that Ks returns its thickness and the plough layer's head."""
    for ret in ("ash1", "mid", "ash2"):
        soil = sh.andosol_soil(ret)
        ks, h_s = sh.pan_ks(soil, return_head=True)
        assert 0.0 < ks < soil["percolation"]
        top = sh._pan_profile(soil, ks, h_s, np.array([sh.PAN_THICKNESS]))[0]
        assert top == pytest.approx(soil["pond"] + sh.PLOUGH_DEPTH
                                    - soil["percolation"] * sh.PLOUGH_DEPTH / soil["Ks"], rel=1e-3)


def test_irrigation_source_initial_pore_water():
    """PREREG 9.3: the plough layer's initial pore water holds the measured pre-season soil at
    equilibrium, C_T = c (Kd + theta_s / rho); the subsoil runs at half (9.9)."""
    soil = sh.andosol_soil()
    kw, _ = sh.yamazaki_source_kw("PFOA", 3.6, soil)
    pre, _post = sh.post_season_soil("PFOA")
    assert kw["c_init"] * (3.6 + soil["ths"] / sh.RHO_PLOUGH) == pytest.approx(pre)
    assert kw["c_sub"] == pytest.approx(0.5 * kw["c_init"])
    assert kw["c_irrigation"] == pytest.approx(8.61) and kw["source"] == "irrigation"


def test_air_entry_conductivity():
    """PREREG 9.12: the hydraulic function has HYDRUS model 3's air entry (-2 cm): saturated above
    it, continuous at it, falling below; h_entry = 0 is the plain van Genuchten-Mualem."""
    soil = sh.andosol_soil()
    vg = (soil["ths"], soil["thr"], soil["alpha"], soil["n"], soil["Ks"])
    assert sh._vg_K(-1.0, *vg) == pytest.approx(soil["Ks"])
    assert sh._vg_K(-2.0 - 1e-9, *vg) == pytest.approx(soil["Ks"], rel=1e-6)
    assert sh._vg_K(-50.0, *vg) < sh._vg_K(-10.0, *vg) < soil["Ks"]
    m = 1 - 1 / soil["n"]
    se = (1 + (soil["alpha"] * 2.0) ** soil["n"]) ** -m
    plain = soil["Ks"] * se ** 0.5 * (1 - (1 - se ** (1 / m)) ** m) ** 2
    assert sh._vg_K(-2.0, *vg, h_entry=0.0) == pytest.approx(plain)
    assert plain < 0.2 * soil["Ks"]                    # the near-saturation drop model 3 removes


def test_after_drainage_options():
    with pytest.raises(ValueError, match="after_drainage"):
        sh.run_paddy_hydrus(1.0, 10.0, after_drainage="wet", tpot_fn=np.ones_like, evap_fn=np.ones_like)


# ---- engine runs -----------------------------------------------------------
import forcing_rice as _fr                                    # noqa: E402
_GEN = dict(tpot_fn=lambda d: _fr.transpiration_mm_d(d, 60.0),
            evap_fn=lambda d: _fr.evaporation_mm_d(d, 60.0))
@hydrus
def test_run_paddy_short_chain_leaches():
    """A weakly-sorbed short chain must leach: pore water falls during flooding."""
    res = sh.run_paddy_hydrus(sh.paddy_kd(4, "PFCA"), season=60.0, **_GEN)
    assert np.all(np.isfinite(res.Cw)) and np.all(np.isfinite(res.vroot))
    assert res.Cw[0] == pytest.approx(1.0, abs=0.05)        # initial pore water
    assert res.Cw.min() < 0.5 * res.Cw[0]                   # leached away
    assert res.vroot.max() > 0                              # roots take up water


@hydrus
def test_run_paddy_long_chain_buffered():
    """A strongly-sorbed long chain stays buffered: pore water ~ flat."""
    res = sh.run_paddy_hydrus(sh.paddy_kd(12, "PFCA"), season=60.0, **_GEN)
    assert res.Cw.min() > 0.9 * res.Cw[0]                   # barely moves


@hydrus
def test_inputs_from_hydrus_normalised_and_shaped():
    inp, res = sh.inputs_from_hydrus(8, "PFCA", season=60.0, Cwo_ref=2.0, n_t=121,
                                     scenario="generic")
    assert inp.Cwo.shape == inp.t.shape == (121,)
    assert np.mean(inp.Cwo) == pytest.approx(2.0, rel=1e-6)  # normalised to ref
    assert inp.M.shape == (121, 4)
    assert np.all(inp.M > 0) and np.all(np.isfinite(inp.Qtp))


@hydrus
def test_inputs_from_hydrus_honours_the_crop_selection():
    """The HYDRUS coupling builds the crop side from the selected scenario: the site-year
    by default (the ORYZA leaf senesces and carries its death rate as leaf_loss), the
    generic drivers with a biomass driver otherwise (growth_rice: monotonic leaf, no rate)."""
    import site_scenario as ss
    li = 2                                                   # root,stem,leaf,grain
    io_, _ = sh.inputs_from_hydrus(8, "PFCA")
    assert io_.t[-1] == pytest.approx(ss.yamazaki2020_mean_run()["t"][-1])
    ig, _ = sh.inputs_from_hydrus(8, "PFCA", season=120.0, n_t=121, biomass="growth_rice",
                                  scenario="generic")
    Lo, Lg = io_.M[:, li], ig.M[:, li]
    assert io_.leaf_loss is not None and Lo[-1] < 0.7 * Lo.max()      # ORYZA senesces
    assert ig.leaf_loss is None and Lg[-1] >= 0.95 * Lg.max()         # growth_rice monotonic


@hydrus
def test_hydrus_qtp_matches_the_crop_transpiration_when_unstressed():
    """qtp_from_hydrus (default): HYDRUS's actual root water uptake reproduces the crop's
    potential transpiration (the site-year's ORYZA2000 Tp) all season -- the soil run imposes
    no stress the stress-free crop driver does not have (PREREG 9.10); the ponded paddy holds its
    water and percolates at the sourced rate."""
    import site_scenario as ss
    run = ss.yamazaki2020_mean_run()
    inp, res = sh.inputs_from_hydrus(8, "PFCA")
    on = inp.t >= 1.0                                     # HYDRUS starts on day 1
    assert np.allclose(inp.Qtp[on], run["Qtp"][on], rtol=0.02, atol=1e-3)


@hydrus
def test_hydrus_drivers_run_the_default_configuration():
    """The default shoot needs the drivers' DVS (the panicle's transpiration share follows it);
    hydrus_drivers once dropped it, and simulate(drivers=hydrus_drivers(...)) raised."""
    import model_api as api
    drv, _ = api.hydrus_drivers("PFOA", Cwo_ref=1.0)
    assert "dvs" in drv and len(drv["dvs"]) == len(drv["t"])
    r = api.simulate("PFOA", drivers=drv)
    assert r["success"]


@hydrus
def test_water_flow_is_converged_and_balanced(tmp_path):
    """PREREG 9.12: the pan's material spans the 15 and 17 cm nodes, so HYDRUS reproduces the pan
    solve's percolation at the start (within 3 %), and the water balance closes (< 0.5 %)."""
    import re
    soil = sh.andosol_soil()
    sh.run_paddy_hydrus(500.0, 20.0, workspace=str(tmp_path), keep=True, **_GEN)
    rows = [l.split() for l in open(tmp_path / "T_LEVEL.OUT").read().splitlines()
            if l.strip()[:1].isdigit()]
    v_bot = -float(rows[0][5])
    assert v_bot == pytest.approx(soil["percolation"], rel=0.03)
    bal = re.findall(r"WatBalR\s*\[%\]\s*([-\d.E+]+)", open(tmp_path / "BALANCE.OUT").read())
    assert abs(float(bal[-1])) < 0.5
