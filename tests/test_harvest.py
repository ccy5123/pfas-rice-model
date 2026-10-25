"""The harvest observation operator (src/harvest.py).

Pinned: (1) with no leaf death the operator is the plain dry-basis conversion and a mass-weighted
straw at each time (the old straw series used the FINAL masses at every time); (2) the dead pool
on a case with a closed form; (3) mass balance on a real solve -- the burden in the four living
organs plus the dead blades equals the cumulative root uptake (the model's sole source), so the
dead pool holds exactly what the leaf ODE loses; (4) attached=0 is the green-only straw.
"""
import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

import harvest as hv            # noqa: E402
import model_api as api         # noqa: E402
from pfas_rice_plant_module_4pool_surf import root_uptake, binding_factors, Environment  # noqa: E402

THETA = dict(root=0.9, stem=0.8, leaf=0.75, grain=0.14)


def test_no_leaf_death_is_the_dry_conversion():
    t = np.linspace(0, 10, 11)
    M = np.column_stack([np.full(11, 1.0), np.linspace(1, 3, 11), np.linspace(2, 1, 11), np.full(11, 0.5)])
    conc = dict(root=np.full(11, 2.0), stem=np.full(11, 1.0), leaf=np.full(11, 4.0), grain=np.full(11, 3.0))
    o = hv.harvest_observables(t, conc, M, hv.dry_fractions(THETA, "fresh"))
    assert np.allclose(o["conc_dw"]["root"], 2.0 / 0.1) and np.allclose(o["conc_dw"]["grain"], 3.0 / 0.86)
    s = (1.0 * M[:, 1] + 4.0 * M[:, 2]) / (0.2 * M[:, 1] + 0.25 * M[:, 2])
    assert np.allclose(o["conc_dw"]["straw"], s)                 # masses at each time
    assert np.allclose(o["dead_burden"], 0.0) and np.allclose(o["conc_dw"]["leaf"], 16.0)


def test_dead_pool_closed_form():
    t = np.linspace(0, 20, 201)
    C, Ml, k, d = 5.0, 2.0, 0.03, 0.25
    b, m = hv.dead_leaf_pool(t, np.full_like(t, C), np.full_like(t, Ml), np.full_like(t, k), d)
    assert np.allclose(b, k * C * Ml * t) and np.allclose(m, k * Ml * d * t)


@pytest.mark.parametrize("config, kcw", [(dict(api.FOUR_POOL_FREE), "root"),
                                         (dict(shoot="redistributed", lipid_loading=True), "stem")])
def test_mass_balance_with_dead_blades(config, kcw):
    # constant organ masses (no growth-dilution term) and a constant leaf death rate: the burden
    # in the living organs plus the dead blades must equal the cumulative root uptake -- for the
    # 4-pool core and for the default (redistributed shoot + lipid loading, whose compound takes
    # the stem cell-wall K everywhere, as in every scored redistributed run)
    t = np.linspace(0.0, 60.0, 1201)
    M = np.tile([0.05, 0.10, 0.03, 0.04], (len(t), 1))
    drv = dict(t=t, Cwo=np.ones_like(t), Qtp=np.full_like(t, 0.2), M=M,
               leaf_loss=np.full_like(t, 0.02), mass_basis="fresh")
    r = api.simulate("PFOA", drivers=drv, **config)
    assert r["dead_leaf"]["burden"][-1] > 0
    living = sum(r["conc"][k] * M[:, i] for i, k in enumerate(api.TISSUES))
    total = living + r["dead_leaf"]["burden"]
    cmpd = api._compound_for("PFOA", "recommended", K_cw_organ=kcw,
                             lipid_loading=config["lipid_loading"])
    B = binding_factors(api._compartments(), cmpd)
    env = Environment(E=-0.120)
    j = np.array([root_uptake(1.0, cr / B[0], cmpd, env) for cr in r["conc"]["root"]]) * M[:, 0]
    uptake = np.concatenate([[0.0], np.cumsum(0.5 * (j[1:] + j[:-1]) * np.diff(t))])
    # compared as increments after day 1: the root fills within hours, faster than the grid
    k = int(np.searchsorted(t, 1.0))
    assert total[-1] - total[k] == pytest.approx(uptake[-1] - uptake[k], rel=1e-4)
    # the dead blades hold what the leaf lost at the death rate: integral of k * C_leaf * M_leaf
    lost = np.trapezoid(0.02 * r["conc"]["leaf"] * M[:, 2], t)
    assert r["dead_leaf"]["burden"][-1] == pytest.approx(lost, rel=1e-3)
    # the dead blades' dry mass grows at the death rate times the leaf's dry matter
    assert r["dead_leaf"]["dry"][-1] == pytest.approx(0.02 * 0.03 * (1 - api._COMP["leaf"]["theta_fw"]) * 60, rel=1e-6)


def test_attached_zero_is_green_only():
    a1 = api.simulate("PFOA")
    a0 = api.simulate("PFOA", dead_leaf_attached=0.0)
    assert a0["dead_leaf"]["burden"][-1] == 0.0
    assert a0["baf_final_dw"]["leaf"] == pytest.approx(a0["baf_final_dw"]["leaf_green"])
    assert a1["baf_final_dw"]["leaf_green"] == pytest.approx(a0["baf_final_dw"]["leaf_green"])
    assert a1["baf_final_dw"]["straw"] != pytest.approx(a0["baf_final_dw"]["straw"])
    with pytest.raises(ValueError):
        api.simulate("PFOA", dead_leaf_attached=1.5)
