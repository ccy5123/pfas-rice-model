"""theta_fw sensitivity (validation/theta_fw_sensitivity.py, pre-registered).

Pinned: (1) a shoot organ's water fraction cancels on the DRY basis -- it moves the
fresh-weight BAF but not the dry one; (2) the ROOT's water fraction is the sensitive
one, because its water term theta/(1-theta) ~ 9 dominates and the root is every TF's
denominator; (3) the open issue the run exposed: the 4-pool core splits the xylem
between leaf and grain by S*M with M FRESH, so the grain's water fraction sets its share
of the transpiration stream. On dry mass that share is theta-free. When the basis of S
is decided, (3) is the test that has to change.
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


def test_open_issue_the_xylem_split_runs_on_fresh_grain_mass():
    S_leaf, S_grain = 20.0, 2.0
    shares = {}
    for th in (0.14, 0.60):
        with T.theta("grain", th):
            r = T.api.simulate("PFOA")
            ml, mg = r["M"][-1, 2], r["M"][-1, 3]
            fresh = S_grain * mg / (S_leaf * ml + S_grain * mg)
            dl, dg = ml * (1 - T.om.THETA_FW["leaf"]), mg * (1 - th)
            shares[th] = (fresh, S_grain * dg / (S_leaf * dl + S_grain * dg))
    assert shares[0.60][0] > 1.5 * shares[0.14][0]                       # 0.118 -> 0.223
    assert shares[0.60][1] == pytest.approx(shares[0.14][1], rel=1e-9)    # dry: 0.343 both


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
