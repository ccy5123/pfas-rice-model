"""Tests for the ORYZA2000 Level-1 (Python) growth core, oryza_growth.py."""
import numpy as np
import pytest

import oryza_growth as oz


def test_phenology_and_canopy_anchors():
    """IR72 potential run hits flowering/maturity/LAI/HI field anchors."""
    sim = oz.simulate_oryza()
    t = sim["dvs"]
    # DVS is non-decreasing and spans emergence -> maturity
    assert np.all(np.diff(sim["dvs"]) >= -1e-9)
    assert sim["dvs"][-1] == pytest.approx(2.0, abs=0.05)
    i_fl = int(np.argmax(sim["dvs"] >= 1.0))
    i_mat = int(np.argmax(sim["dvs"] >= 2.0))
    assert 55 <= sim["t"][i_fl] <= 75           # flowering ~ day 65
    assert 105 <= sim["t"][i_mat] <= 120         # maturity ~ day 115
    assert 4.5 <= sim["lai"].max() <= 7.5        # realistic peak LAI
    shoot = sim["wlv"][-1] + sim["wst"][-1] + sim["wso"][-1]
    assert 0.40 <= sim["wso"][-1] / shoot <= 0.52    # harvest index


def test_organ_biomass_scaled_to_anchor():
    """Anchored output reproduces the IR72 shoot anchor (~1740 g/m^2 -> kg/hill)."""
    t = np.linspace(0, 120.0, 241)
    b = oz.organ_biomass_oryza(t, scale_to_anchor=True)
    assert {"root", "stem", "leaf", "grain"} <= set(b)
    for k in ("root", "stem", "leaf", "grain"):
        assert b[k].shape == t.shape and np.all(b[k] > 0)
    # leaf senescence loss RATE [1/d] is also exposed (>=0, nonzero after flowering)
    assert b["leaf_death_rate"].shape == t.shape and np.all(b["leaf_death_rate"] >= 0)
    assert np.max(b["leaf_death_rate"]) > 0
    shoot = b["stem"][-1] + b["leaf"][-1] + b["grain"][-1]
    target = oz.OryzaParams().shoot_anchor_g_m2 * 10.0 * oz.HA_PER_HILL   # kg/hill
    assert shoot == pytest.approx(target, rel=1e-2)


def test_leaf_senescence_and_grain_fill():
    """Mechanistic signatures absent from the logistic reconstruction:
    leaves senesce (decline from peak) and grain fills only after flowering."""
    t = np.linspace(0, 120.0, 481)
    b = oz.organ_biomass_oryza(t)
    # leaf declines from its mid-season peak by harvest (senescence)
    assert b["leaf"][-1] < 0.85 * b["leaf"].max()
    # grain ~0 well before flowering, large at harvest, and monotone while filling
    assert b["grain"][np.argmin(abs(t - 30))] < 0.02 * b["grain"][-1]
    assert b["grain"][-1] > 0.3 * (b["stem"][-1] + b["leaf"][-1] + b["grain"][-1])
    i0 = int(np.argmin(abs(t - 80)))
    assert np.all(np.diff(b["grain"][i0:]) >= -1e-12)


def test_M_matrix_and_drivers_shapes():
    t = np.linspace(0, 120.0, 200)
    M = oz.M_matrix_oryza(t)
    assert M.shape == (200, 4) and np.all(M > 0)
    dr = oz.oryza_drivers("PFOA", Cwo=1.0, n_t=200)
    assert set(("t", "Cwo", "Qtp", "M")).issubset(dr)
    assert dr["M"].shape == (200, 4) and np.all(dr["M"] > 0)
    assert np.allclose(dr["Cwo"], 1.0)


def test_end_to_end_through_pfas_ode():
    """ORYZA biomass drivers run through the full 4-pool PFAS ODE."""
    api = pytest.importorskip("model_api")
    dr = oz.oryza_drivers("PFOA", Cwo=1.0)
    r = api.simulate("PFOA", drivers=dr)
    for k in ("root", "grain"):
        assert np.isfinite(r["baf_final"][k]) and r["baf_final"][k] > 0
    assert np.isfinite(r["straw_baf"]) and r["straw_baf"] > 0


def test_weather_override_increases_potential_biomass():
    """More radiation -> larger UNSCALED potential shoot (model responds to weather)."""
    t = np.linspace(0, 120.0, 241)
    n = int(np.ceil(oz.OryzaParams().season)) + 1
    lo = dict(tmax=np.full(n, 30.0), tmin=np.full(n, 22.0), rad_mj=np.full(n, 14.0))
    hi = dict(tmax=np.full(n, 30.0), tmin=np.full(n, 22.0), rad_mj=np.full(n, 24.0))
    s_lo = oz.organ_biomass_oryza(t, weather=lo, scale_to_anchor=False)
    s_hi = oz.organ_biomass_oryza(t, weather=hi, scale_to_anchor=False)
    shoot_lo = s_lo["stem"][-1] + s_lo["leaf"][-1] + s_lo["grain"][-1]
    shoot_hi = s_hi["stem"][-1] + s_hi["leaf"][-1] + s_hi["grain"][-1]
    assert shoot_hi > 1.15 * shoot_lo


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))
