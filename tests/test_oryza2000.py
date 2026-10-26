"""ORYZA2000 potential production, ported from the manual (src/oryza2000.py), and its
pre-registered IR72DS92 comparison (validation/oryza2000_ir72ds92.py, 288de63).

Pinned: (1) routines against values the equations fix exactly (SUBDD, SASTRO, the GPPARGET
branch continuity); (2) the carbon balance closes, as ORYZA's own SUBCBC demands; (3) the
crop-stage logic -- emergence -> 12-d seedbed -> transplanting on doy 16, the dilution by
NPLH*NH/NPLSB, no development during the transplanting shock; (4) the pre-registered verdict
(R1 ADEQUATE under every weather treatment, R2 the port beats the shipped driver); (5) the
model_api mode is opt-in and the default is untouched.
"""
import contextlib
import io
import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "validation"))

import oryza2000 as O            # noqa: E402
import oryza2000_ir72ds92 as V   # noqa: E402


def _synthetic(n=220):
    d = np.arange(n)
    return dict(doy0=1, rad_mj=19 + 3 * np.sin(2 * np.pi * (d - 20) / 120.0),
                tmin=23 + np.sin(2 * np.pi * (d - 30) / 120.0),
                tmax=31 + 1.5 * np.sin(2 * np.pi * (d - 30) / 120.0))


def test_routines_against_exact_values():
    assert O.subdd(30.0, 30.0, 8.0, 30.0, 42.0) == pytest.approx(22.0)       # constant at TOD
    assert O.subdd(5.0, 5.0, 8.0, 30.0, 42.0) == 0.0                         # below TBD
    assert O.sastro(80, 0.0)[2] == pytest.approx(12.0, abs=0.01)              # equinox, equator
    gp = (340.0, 0.4, 0.5, 1.0)
    hi = O.gppar_get(0.0, 0.0, 0.45, gp)[0]
    lo = O.gppar_get(0.0, 0.0, 0.45, (340.0, 0.4, 0.5 - 1e-9, 1.0))[0]
    assert hi == pytest.approx(lo, abs=2e-3)               # branches meet (to the printed 68.33's rounding)
    assert hi == pytest.approx(20.5, abs=0.01)                                # 9.5 + 22*0.5 at 340 ppm
    crop = O.load_crop()
    assert (crop["DVRJ"], crop["DVRR"], crop["SPGF"]) == (0.000773, 0.001784, 64900.0)
    assert list(crop["FSTTB"][1]) == [0.40, 0.40, 0.70, 0.40, 0.0, 0.0]       # no stem growth after DVS 1.2


def test_carbon_balance_and_crop_stages():
    r = O.simulate(_synthetic(), O.load_crop(), O.Management())
    assert r["cbc_max"] < 1e-9                                                # SUBCBC
    assert r["events"]["transplanting"] == 16                                 # emergence 4 + 12-d seedbed
    i3 = int(np.where(r["cropsta"] == 3)[0][0])
    pltr = 5 * 25 / 1000
    # transplanting day: the organs are diluted by PLTR, and -- as printed -- that day's growth is
    # still computed on the seedbed canopy, so one seedbed day of growth is added on top
    added = r["wlvg"][i3 + 1] - pltr * r["wlvg"][i3]
    seedbed_day = r["wlvg"][i3] - r["wlvg"][i3 - 1]
    assert 0.5 < added / seedbed_day < 2.0
    # transplanting shock: development stops for SHCKD * TS(transplant) degree-days
    dvs_after = r["dvs"][i3 + 1:i3 + 4]
    assert np.ptp(dvs_after) == pytest.approx(0.0, abs=1e-12)
    assert r["dvs"][-1] > 2.0


def test_organ_biomass_bridge():
    t = np.linspace(0.0, 120.0, 481)
    b = O.organ_biomass(t, 120.0)
    t_mat = t[np.argmax(b["dvs"] >= 2.0)]
    assert t_mat == pytest.approx(93.0, abs=1.0)
    assert np.all(b["leaf_death_rate"][t >= t_mat] == 0.0)                   # frozen after maturity
    assert b["grain"][-1] > b["leaf"][-1] and b["grain"][0] < 1e-6


@pytest.fixture(scope="module")
def verdict():
    with contextlib.redirect_stdout(io.StringIO()):
        return V.main(posthoc=False)


def test_station_bias_from_the_printed_days(verdict):
    dn, dx, rr = verdict["bias"]
    assert (dn, dx, rr) == pytest.approx((1.17, 1.67, 1.211), abs=0.006)


def test_the_pre_registered_verdict(verdict):
    assert verdict["r1_verdict"] == "ADEQUATE"
    assert all(verdict["r1"].values())                                       # every weather treatment
    assert verdict["r2"] and verdict["r2_wins"] == 4
    prim = verdict["results"]["corrected"]
    assert prim["A"]["ph"] == {"panicle_initiation": 56, "flowering": 82, "maturity": 115}
    assert prim["B"]["ph"]["maturity"] == 126
    assert prim["A"]["sc"]["WAGT"] == pytest.approx(0.024, abs=0.005)
    assert prim["B"]["sc"]["WAGT"] == pytest.approx(0.236, abs=0.01)


def test_model_api_generic_biomass_modes():
    """On the generic drivers the port and the reduced driver are both selectable and give
    different crops; an unknown name raises instead of falling back silently. (The default
    scenario runs the port on the Yamazaki site-year -- tests/test_site_scenario.py.)"""
    import model_api as api
    port = api.simulate("PFOA", scenario="generic", season=120.0, biomass="oryza2000")
    reduced = api.simulate("PFOA", scenario="generic", season=120.0, biomass="oryza")
    assert not np.allclose(port["M"][-1], reduced["M"][-1])
    with pytest.raises(ValueError):
        api.simulate("PFOA", scenario="generic", season=120.0, biomass="ORYZA2000")


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
