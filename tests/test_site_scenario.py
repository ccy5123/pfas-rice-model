"""The Yamazaki 2020 site-year scenario (src/site_scenario.py) and its FAO-56 / ORYZA2000 water
terms (src/fao56.py).

Pinned: (1) FAO-56 reproduces the manual's worked Example 18 (Brussels, 6 July: ET0 3.9 mm/d) and
its wind conversion (Example 14: 3.2 m/s at 10 m -> 2.4 m/s); (2) ORYZA2000 eqs. 4.2-4.4 at their
limits; (3) the JMA Tsukuba series loads contiguous over the season; (4) each plant run puts DVS 1
on the heading day and DVS 2 on the harvest day the paper states, closes its carbon balance, and
gives transpiration in L/d/hill with a closed-canopy T/ET of 0.55-0.70.
"""
import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

import fao56 as F               # noqa: E402
import site_scenario as S       # noqa: E402


def test_fao56_example_18():
    lat, doy, tmax, tmin = 50 + 48 / 60, 187, 21.5, 12.3
    ea = (F.es(tmin) * 84 + F.es(tmax) * 63) / 200
    assert ea == pytest.approx(1.409, abs=0.001)
    ra = F.extraterrestrial_radiation(doy, lat)
    assert ra == pytest.approx(41.09, abs=0.01)
    n_max = 24 / np.pi * np.arccos(-np.tan(np.deg2rad(lat))
                                   * np.tan(0.409 * np.sin(2 * np.pi * doy / 365 - 1.39)))
    rs = (0.25 + 0.5 * 9.25 / n_max) * ra
    assert rs == pytest.approx(22.07, abs=0.01)
    u2 = F.wind_2m(10 / 3.6, 10)
    assert u2 == pytest.approx(2.078, abs=0.001)
    assert F.net_radiation(rs, tmax, tmin, ea, doy, lat, 100) == pytest.approx(13.28, abs=0.01)
    et0, rad, aero = F.et0_terms(tmax, tmin, rs, ea, u2, doy, lat, 100)
    assert et0 == pytest.approx(3.9, abs=0.05) and rad + aero == pytest.approx(et0)
    assert F.wind_2m(3.2, 10) == pytest.approx(2.4, abs=0.01)


def test_oryza_water_equations_at_their_limits():
    assert F.field_reflection(0.0) == pytest.approx(0.05)
    assert F.field_reflection(0.0, flooded=False) == pytest.approx(0.25)
    assert F.field_reflection(50.0) == pytest.approx(0.25)
    assert F.oryza_transpiration(4.0, 1.0, 0.0) == 0.0
    assert F.oryza_evaporation(4.0, 1.0, 0.0) == pytest.approx(5.0)
    assert F.oryza_transpiration(4.0, 1.0, 50.0) == pytest.approx(6.0)     # ETAE x min(2, LAI)


def test_tsukuba_weather_loads():
    w = S.load_tsukuba_2020()
    assert w["doy0"] == 75 and len(w["rad_mj"]) == 231
    assert np.all(w["tmax"] >= w["tmin"]) and np.all(w["rad_mj"] > 0) and np.all(w["u2"] > 0)


@pytest.fixture(scope="module")
def runs():
    return S.yamazaki2020_runs()


def test_plant_runs_hit_the_reported_dates(runs):
    tp = S._doy(S.YAMAZAKI_TRANSPLANT)
    for r, h, hd in zip(runs, S.YAMAZAKI_HARVESTS, S.YAMAZAKI_HEADINGS):
        hv = S._doy(h)
        assert r["harvest_doy"] == hv and r["heading_doy"] == S._doy(hd)     # SI p.S3
        t_head = np.interp(1.0, r["dvs"], r["t"])
        assert tp + t_head == pytest.approx(r["heading_doy"], abs=0.05)
        assert r["dvs"][-1] == pytest.approx(2.0, abs=1e-5) and r["t"][-1] == hv - tp
        assert abs(r["cbc_max"]) < 1e-10


def test_fractional_nursery_interpolates_the_whole_day_runs():
    """ORYZA counts the seedbed in whole days; the sourced midpoint (25.5 d) is the two
    neighbouring runs' drivers interpolated, and ORYZA itself refuses a fraction."""
    w = S.load_tsukuba_2020()
    h = S.YAMAZAKI_HARVESTS[0]
    a, b = S.plant_run(h, w=w, sbdur=25), S.plant_run(h, w=w, sbdur=26)
    m = S.plant_run(h, w=w, sbdur=25.5)
    assert S.SBDUR == 25.5 and m["sbdur_runs"] == (25, 26, 0.5)
    for k in ("root", "stem", "leaf", "grain", "Qtp", "dvs", "leaf_death_rate"):
        assert np.allclose(m[k], 0.5 * (a[k] + b[k]), rtol=0, atol=1e-15)
    assert np.array_equal(m["t"], a["t"])
    with pytest.raises(ValueError):
        S._mgmt(25.5)


def test_water_terms(runs):
    for r in runs:
        tp, ep = r["Tp"].sum(), r["Ep"].sum()
        assert 0.55 < tp / (tp + ep) < 0.70
        assert np.allclose(r["Qtp"], r["Tp"] / S.NH)              # mm/d x m2/hill = L/d/hill
        assert 6.0 < r["Tp"].max() < 8.0
        assert np.all(r["Ep"] >= 0) and np.all(np.diff(r["dvs"]) >= -1e-12)


def test_an_explicit_season_is_refused_not_ignored():
    """The site scenario sets its own season (transplanting to harvest). A different
    explicit season used to be ignored silently -- the caller believed they ran 60 days
    and got the site-year. It must raise; the site's own value and None are accepted."""
    import model_api as api
    site = api.site_season()
    with pytest.raises(ValueError, match="does not apply"):
        api.simulate("PFOA", season=60.0)
    with pytest.raises(ValueError, match="does not apply"):
        api.simulate_neutral(2.45, season=120.0)
    with pytest.raises(ValueError, match="does not apply"):
        api.simulate_nstem_leaf("PFOA", season=150.0)
    assert api.simulate("PFOA", season=site)["season"] == pytest.approx(site)
    g = api.simulate("PFOA", season=90.0, scenario="generic")       # generic: season applies
    assert g["season"] == pytest.approx(90.0)
