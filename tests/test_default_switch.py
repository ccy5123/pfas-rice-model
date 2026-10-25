"""The default configuration after the pre-registered switch (validation/default_switch.py).

Pinned: (1) `simulate()` runs the redistributed shoot with lipid loading and reports it;
(2) that path IS the configuration the transfer tests scored (`simulate_nstem_leaf` with
lipid loading) -- G4 of the switch, here on two congeners; (3) the previous default stays
reachable as `FOUR_POOL_FREE`; (4) the source apportionment of the default shoot (PFOA: the
grain is ~80 % xylem-fed, the root ~3/4 soil uptake and ~1/4 phloem recirculation).
"""
import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

import model_api as api  # noqa: E402


def test_the_default_is_the_redistributed_shoot_with_lipid_loading():
    assert (api.DEFAULT_SHOOT, api.DEFAULT_LIPID_LOADING) == ("redistributed", True)
    p = api.simulate("PFOA")["params"]
    assert p["shoot"] == "redistributed" and p["lipid_loading"] is True
    q = api.simulate("PFOA", **api.FOUR_POOL_FREE)["params"]
    assert q["shoot"] == "4pool" and q["lipid_loading"] is False
    with pytest.raises(ValueError):
        api.simulate("PFOA", shoot="two-pool")


@pytest.mark.parametrize("congener", ["PFOA", "PFOS"])
def test_the_default_is_the_scored_configuration(congener):
    a = api.simulate(congener)["baf_final_dw"]
    b = api.simulate_nstem_leaf(congener, lipid_loading=True)["baf_final_dw"]
    for k in a:
        assert a[k] == b[k], k


def test_apportionment_of_the_default_shoot():
    ap = api.apportionment("PFOA")
    g = ap["fraction"]["grain"]
    assert g["xylem_from_stem"] + g["phloem_from_leaf"] == pytest.approx(1.0)
    assert ap["grain_phloem_fraction"] == pytest.approx(0.197, abs=0.02)
    r = ap["fraction"]["root"]
    assert r["soil_uptake"] == pytest.approx(0.735, abs=0.02)
    assert np.isfinite(ap["delivered_to_shoot"]) and ap["success"]
