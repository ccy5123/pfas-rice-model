"""Organ water theta_fw: properties of the model, pinned.

(1) A shoot organ's water fraction cancels on the DRY basis: it moves the fresh-weight BAF
but not the dry one. (2) The ROOT's water fraction is the sensitive one: its water term
theta/(1-theta) is large and the root is every TF's denominator. (3) The 4-pool's leaf/grain
xylem split is the panicle's sourced share of transpiration, so no organ's water sets it.
theta is changed everywhere it is read at once (model_api._COMP, organ_mass.THETA_FW,
neutral_dpu.RICE_WATER). The analysis that first found these (validation/archive/
theta_fw_sensitivity.py) ran on ranges without a source; the sourced ranges are carried by
validation/input_range_band.py.
"""
import contextlib
import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

import model_api as api       # noqa: E402
import neutral_dpu as ND      # noqa: E402
import organ_mass as om       # noqa: E402

KEY = {"root": "root", "stem": "stem", "leaf": "leaf", "grain": "grain_brown"}


@contextlib.contextmanager
def theta(organ, value):
    """Set one organ's theta_fw everywhere it is read, then restore it."""
    key = KEY[organ]
    saved = (api._COMP[key]["theta_fw"], om.THETA_FW[organ], ND.RICE_WATER[organ])
    api._COMP[key]["theta_fw"] = om.THETA_FW[organ] = ND.RICE_WATER[organ] = float(value)
    try:
        yield
    finally:
        api._COMP[key]["theta_fw"], om.THETA_FW[organ], ND.RICE_WATER[organ] = saved


def _bafs(organ, value, name="PFOA"):
    """Root/straw/grain BAF on the fresh and on the dry basis (harvest operator), with the
    organ's water held at `value` through the season: its stage anchors (PREREG 7.1) both at
    `value`, and the constant copies too."""
    anchors = ({organ: [float(value)] * len(om.THETA_STAGE[organ]["dvs"])}
               if organ in om.THETA_STAGE else None)            # the grain is not staged
    with theta(organ, value):
        r = api.simulate(name, **api.FOUR_POOL_FREE, scenario_kw={"theta_anchors": anchors})
    C = {k: r["conc"][k][-1] / r["cwo_ref"] for k in api.TISSUES}
    fw = dict(root=C["root"], straw=r["straw_baf"], grain=C["grain"])
    dw = {k: r["baf_final_dw"][k] for k in ("root", "straw", "grain")}
    return fw, dw


def test_the_context_manager_restores_every_copy():
    before = (api._COMP["stem"]["theta_fw"], om.THETA_FW["stem"], ND.RICE_WATER["stem"])
    with theta("stem", 0.5):
        assert api._COMP["stem"]["theta_fw"] == om.THETA_FW["stem"] == ND.RICE_WATER["stem"] == 0.5
    assert (api._COMP["stem"]["theta_fw"], om.THETA_FW["stem"], ND.RICE_WATER["stem"]) == before


def test_shoot_theta_cancels_on_the_dry_basis():
    lo_fw, lo_dw = _bafs("stem", 0.645)              # the sourced stem range
    hi_fw, hi_dw = _bafs("stem", 0.774)
    d = lambda a, b: abs(np.log10(a) - np.log10(b))
    assert d(lo_dw["straw"], hi_dw["straw"]) < 0.02
    assert d(lo_fw["straw"], hi_fw["straw"]) > 2 * d(lo_dw["straw"], hi_dw["straw"])


def test_root_theta_is_the_sensitive_one():
    _, lo = _bafs("root", 0.713)                     # the sourced value
    _, hi = _bafs("root", 0.90)                      # the unsourced value it replaced
    assert abs(np.log10(hi["root"]) - np.log10(lo["root"])) > 0.02


def test_the_grain_share_of_the_transpiration_stream():
    """Both shoots give the panicle its sourced share of transpiration (PREREG 5.1): 0 before
    heading, the heading and maturity anchors at DVS 1 and 2, linear in DVS between -- not
    the unmeasured transpiring surfaces the 4-pool used before (leaf 20, grain 0.5116 m2/kg)."""
    from pfas_rice_plant_module_4pool_surf import panicle_share
    run = api.site_run()
    h, m = api.LAM_GRAIN_DVS
    sh = np.array([panicle_share(d, h, m) for d in run["dvs"]])
    assert np.all(sh[run["dvs"] < 1.0] == 0.0)
    assert sh[int(np.argmax(run["dvs"] >= 1.0))] == pytest.approx(h, abs=0.01)   # first day at DVS >= 1
    assert sh[-1] == pytest.approx(m, abs=1e-6)
    assert (h, m) == (api.range_mid("lam_grain_heading"), api.range_mid("lam_grain_maturity"))


def test_grain_water_does_not_move_the_4pool_grain_on_the_dry_basis():
    """The panicle share replaced the surface split whose grain term read the grain's
    water; the dry grain BAF no longer depends on grain theta through the split."""
    _, lo = _bafs("grain", 0.20)                     # the sourced grain range
    _, hi = _bafs("grain", 0.22)
    assert abs(np.log10(hi["grain"]) - np.log10(lo["grain"])) < 0.01


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
