"""theta_fw sensitivity (validation/theta_fw_sensitivity.py, pre-registered).

Pinned: (1) a shoot organ's water fraction cancels on the DRY basis -- it moves the
fresh-weight BAF but not the dry one; (2) the ROOT's water fraction is the sensitive
one, because its water term theta/(1-theta) ~ 9 dominates and the root is every TF's
denominator; (3) the decision the run forced: the 4-pool core split the xylem between
leaf and grain by S*M with M FRESH, so the grain's water set its share of the
transpiration stream. It now splits on DRY matter with S_XYLEM_DRY, which is theta-free
and, at the shipped theta, the same split as before.
"""
import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "validation"))

import theta_fw_sensitivity as T  # noqa: E402


def _bafs(organ, value, name="PFOA"):
    with T.theta(organ, value):
        return T.pfas_bafs(name)


def test_the_context_manager_restores_every_copy():
    import neutral_dpu as ND
    before = (T.api._COMP["stem"]["theta_fw"], T.om.THETA_FW["stem"], ND.RICE_WATER["stem"])
    with T.theta("stem", 0.5):
        assert T.api._COMP["stem"]["theta_fw"] == T.om.THETA_FW["stem"] == ND.RICE_WATER["stem"] == 0.5
    assert (T.api._COMP["stem"]["theta_fw"], T.om.THETA_FW["stem"], ND.RICE_WATER["stem"]) == before


def test_shoot_theta_cancels_on_the_dry_basis():
    lo_fw, lo_dw = _bafs("stem", 0.80)
    hi_fw, hi_dw = _bafs("stem", 0.88)
    d = lambda a, b: abs(np.log10(a) - np.log10(b))
    assert d(lo_dw["straw"], hi_dw["straw"]) < 0.02
    assert d(lo_fw["straw"], hi_fw["straw"]) > 2 * d(lo_dw["straw"], hi_dw["straw"])


def test_root_theta_is_the_sensitive_one():
    _, lo = _bafs("root", 0.88)
    _, hi = _bafs("root", 0.92)
    assert abs(np.log10(hi["root"]) - np.log10(lo["root"])) > 0.02     # 4.68 -> 5.17 dw


def test_the_xylem_split_runs_on_dry_matter_and_is_theta_free():
    from pfas_rice_plant_module_4pool_surf import xylem_leaf_share, S_XYLEM_DRY, LEAF, FRUIT
    shares, M = {}, {}
    for th in (0.14, 0.60):
        with T.theta("grain", th):
            M[th] = T.api.simulate("PFOA")["M"][-1]
            shares[th] = 1.0 - xylem_leaf_share(T.api._compartments(), M[th], "fresh")
    assert shares[0.60] == pytest.approx(shares[0.14], rel=1e-9)          # grain water: no effect
    # at the shipped theta it is the split the fresh-mass S*M (20:2) used to give
    ml, mg = M[0.14][2], M[0.14][3]
    assert shares[0.14] == pytest.approx(2.0 * mg / (20.0 * ml + 2.0 * mg), rel=1e-4)
    assert S_XYLEM_DRY[LEAF] == 20.0 and S_XYLEM_DRY[FRUIT] == 0.5116


def test_the_grain_share_of_the_transpiration_stream():
    """The grain xylem surface S_XYLEM_DRY[grain] = 0.5116 m2/kg dry is unmeasured: on the
    site-year it sends ~6% of the transpiration stream to the grain during filling; the
    original S = 2 read per kg dry would send ~18%. (validation/refit_site_year.py runs
    0.2 / 2.0 as a sensitivity.)"""
    import site_scenario as ss
    from pfas_rice_plant_module_4pool_surf import PlantInputs
    run = ss.yamazaki2020_mean_run()
    t, Q = run["t"], run["Qtp"]
    Md = np.maximum(np.column_stack([run[o] for o in ("root", "stem", "leaf", "grain")]), 1e-4)
    gate = PlantInputs(t=t, Cwo=np.ones_like(t), Qtp=Q, M=Md)._grain_gate(t)
    fill = gate > 0.5

    def filling_share(S_grain):
        g = S_grain * Md[:, 3] / (20.0 * Md[:, 2] + S_grain * Md[:, 3]) * gate
        return np.trapezoid((g * Q)[fill], t[fill]) / np.trapezoid(Q[fill], t[fill])
    now, s2 = filling_share(0.5116), filling_share(2.0)
    assert 0.05 < now < 0.07 and 0.16 < s2 < 0.20, (now, s2)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
