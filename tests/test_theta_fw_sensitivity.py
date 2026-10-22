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


def test_the_rejected_reading_sends_three_times_more_of_the_stream_to_the_grain():
    """Reading the original grain S=2 as per kg DRY was option A as first proposed: the
    grain then takes ~14% of the transpiration stream during filling against ~4% now
    (nstem_leaf assumes 5%), and 8 of 11 Yamazaki grains fall out of reach of the
    oryza re-fit (RMSE 0.283 -> 0.394). Pinned so the reason for 0.5116 stays checkable."""
    import forcing_rice as fr
    from pfas_rice_plant_module_4pool_surf import PlantInputs
    t = np.linspace(0.0, 120.0, 241)
    Q = fr.Q_TP(t, 120.0)
    Md = T.om.organ_matrix(T.api._biomass_fn("oryza")(t, 120.0), "dry")
    gate = PlantInputs(t=t, Cwo=np.ones_like(t), Qtp=Q, M=Md)._grain_gate(t)
    fill = gate > 0.5

    def filling_share(S_grain):
        g = S_grain * Md[:, 3] / (20.0 * Md[:, 2] + S_grain * Md[:, 3]) * gate
        return np.trapezoid((g * Q)[fill], t[fill]) / np.trapezoid(Q[fill], t[fill])
    now, rejected = filling_share(0.5116), filling_share(2.0)
    assert 0.03 < now < 0.06 and 0.12 < rejected < 0.16, (now, rejected)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
