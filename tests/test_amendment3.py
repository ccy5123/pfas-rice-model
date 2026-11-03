"""Amendment 3 (PREREG 10): site desorption Kd, Morita root profile, Artacho crop set, exposure shape hook."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))
import crop_sets as CS          # noqa: E402
import soil_hydrus as sh        # noqa: E402


def test_site_kd_reproduces_the_extraction_ratio():
    r = sh._s16("PFOA")
    R = r["pre_wtmeoh"][0] / r["pre_water"][0]
    m = sh.EXTRACT_WET_G * (1 - sh.EXTRACT_MOISTURE["pre"])
    rr, w0 = sh.EXTRACT_WATER_ML / m, sh.EXTRACT_WET_G * sh.EXTRACT_MOISTURE["pre"] / m
    for vr in (0.0, 0.4):
        K = sh.site_desorption_kd("PFOA", vr)
        F = (vr + K) / (w0 + rr + K) * ((vr + K) / (vr + rr + K)) ** 2
        assert abs(F / (1 - F) - R) < 1e-8


def test_site_kd_bounds():
    vmax = sh.extraction_vr_max()
    assert abs(vmax - 0.9631) < 1e-3
    assert sh.site_desorption_kd("PFBA", vmax) == 0.0
    assert np.isnan(sh.site_desorption_kd("PFBA", vmax + 0.1))
    ks = [sh.site_desorption_kd("PFOA", v) for v in (0.0, 0.5 * vmax, vmax)]
    assert ks[0] > ks[1] > ks[2] > 0                     # V_r = 0 gives the largest Kd


def test_morita_root_profile():
    z = -np.array([1.0, 7.0, 15.0, 25.0, 31.0])
    assert np.allclose(sh.root_density(z, "morita_C"), [7.6, 8.0, 3.1, 1.4, 0.0])
    assert np.allclose(sh.root_density(z, "linear30"), [1 - 1 / 30, 1 - 7 / 30, 0.5, 1 / 6, 0.0])
    assert sh.ROOT_PROFILE_DEFAULT == "morita_C"


def test_artacho_crop_set():
    o = CS.overrides("artacho2011")
    assert "artacho2011" in CS.BAND_SETS
    assert o["DVRJ"] == 0.002126 and o["TOD"] == 28.0 and o["WGRMX"] == 3.043e-5
    assert np.isclose(np.interp(1.0, *o["FSOTB"]), 0.70)


def test_exposure_shape_none_is_the_observable():
    import model_api as api
    b = api.yamazaki_harvest_baf("PFDA", runs="mean")
    c = api.yamazaki_harvest_baf("PFDA", runs="mean", cwo_shape=(np.array([0.0, 200.0]), np.ones(2)))
    for k in b:
        assert abs(np.log10(b[k] / c[k])) < 1e-9
    assert api.site_kd_exposure_shape("PFBA") is None    # not admitted (PREREG 10.6)
    assert api.site_kd_exposure_shape("PFDA") is not None
