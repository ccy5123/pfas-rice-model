"""What Tang 2026's "endosperm" is compared against (model_api.TANG_ENDOSPERM).

The model's grain pool is the whole caryopsis (brown rice); Tang's endosperm is most
likely milled white rice, but the milling method is not reported. So "brown" stays the
default and "white" is the named alternative, which re-partitions the grain's free
concentration into a white-rice composition with membrane phospholipid counted the way
the grain_brown row counts it.

Pinned: the default is unchanged; the white composition is derived as documented; and
the finding -- the white reading LOWERS the model endosperm 1.5-4x, so the documented
grain under-prediction gets WORSE, not better (lipid OOS 0.640 -> 0.831).
"""
import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

import model_api as api  # noqa: E402


def _rmse(**kw):
    e = []
    for c in api.TANG_CONGENERS:
        v = api.tang_tf_validation(c, dose="low", **kw)
        e += [np.log10(v["model_tf"][o]) - np.log10(v["tang_tf"][o]) for o in v["organs"]]
    return float(np.sqrt(np.mean(np.square(e))))


def test_default_is_brown_and_unchanged():
    a = api.tang_tf_validation("PFOA", dose="low")
    b = api.tang_tf_validation("PFOA", dose="low", endosperm="brown")
    assert a["endosperm"] == "brown"
    assert a["model_tf"] == b["model_tf"]
    with pytest.raises(ValueError):
        api.tang_tf_validation("PFOA", endosperm="husked")


def test_white_composition_counts_membrane_pl_like_the_brown_row():
    w = api._grain_white_composition()
    brown_pl = api._COMP["grain_brown"]["f_PL"]
    assert w["f_PL"] == pytest.approx(0.2 * brown_pl / 0.88)
    assert w["f_PL"] < brown_pl                       # bran/embryo PL milled away
    assert (w["theta_fw"], w["f_prot"], w["f_cw"]) == (0.11, 0.075, 0.006)


def test_white_lowers_the_endosperm_most_for_the_high_kpl_sulfonate():
    f = {c: api.endosperm_white_factor(c) for c in api.TANG_CONGENERS}
    assert all(v < 0.7 for v in f.values()), f
    assert f["PFOS"] < f["PFOA"] < f["GenX"], f       # K_PL 31623 > 1905 > 118
    v = api.tang_tf_validation("PFOS", dose="low", endosperm="white")
    vb = api.tang_tf_validation("PFOS", dose="low")
    assert v["model_tf"]["endosperm"] == pytest.approx(vb["model_tf"]["endosperm"] * f["PFOS"])
    assert v["model_tf"]["stalk"] == vb["model_tf"]["stalk"]   # only the endosperm moves


def test_the_white_reading_makes_the_grain_gap_worse():
    assert _rmse(lipid_loading=True) == pytest.approx(0.640, abs=0.02)
    assert _rmse(lipid_loading=True, endosperm="white") == pytest.approx(0.831, abs=0.03)


def test_the_honesty_band_uses_the_default_shoot_on_tang():
    """The band's Tang entry must be the DEFAULT configuration's score (the 4-pool core,
    free anion). It once held the redistributed shoot's 1.496, which the Tang check runs
    unless told otherwise -- the two shoot models score differently on the same data."""
    key = [k for k in api.APRIORI_OOS_SCORES if k.startswith("Tang")][0]
    four = _rmse(shoot="4pool")
    assert api.APRIORI_OOS_SCORES[key] == pytest.approx(four, abs=0.005)
    assert _rmse() == pytest.approx(1.496, abs=0.01)           # redistributed shoot, free anion
    assert four > _rmse() + 0.1
    with pytest.raises(ValueError):
        api.tang_tf_validation("PFOA", shoot="two-pool")


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
