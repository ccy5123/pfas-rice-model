"""What Tang 2026's "endosperm" is compared against (model_api.TANG_ENDOSPERM).

The model's grain pool is the whole caryopsis (brown rice); whether Tang's endosperm was
dehulled only or milled is not reported. So "brown" stays the
default and "white" is the named alternative, which re-partitions the grain's free
concentration into a white-rice composition with membrane phospholipid counted the way
the grain_brown row counts it.

Pinned: the default is unchanged; the white composition is the sourced grain_white row
(params/rice_composition_sourced.csv, midpoints); the white reading lowers the model
endosperm (most for the high-K_PL sulfonate); and the band's Tang entry is the default
configuration's score, the same as validation/transfer_sourced.json records.
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


def test_white_composition_is_the_sourced_row():
    w = api._grain_white_composition()
    assert w == api._COMP["grain_white"]
    rng = api.PARAMS["tissue_composition_ranges"]["grain_white"]
    for q in ("theta_fw", "f_prot", "f_PL", "f_cw"):
        assert w[q] == pytest.approx(0.5 * (rng[q][0] + rng[q][1]), abs=1e-6)
    assert w["f_PL"] < api._COMP["grain_brown"]["f_PL"]          # bran/embryo PL milled away


def test_white_lowers_the_endosperm():
    f = {c: api.endosperm_white_factor(c) for c in api.TANG_CONGENERS}
    assert all(v < 1.0 for v in f.values()), f
    # with the sourced composition the milled-away phospholipid no longer orders them: PFOS and
    # PFOA (K_PL 31623 vs 1905) within 0.5 %, GenX (118) the least lowered
    assert abs(f["PFOS"] - f["PFOA"]) < 0.005 and f["GenX"] > max(f["PFOS"], f["PFOA"]), f
    v = api.tang_tf_validation("PFOS", dose="low", endosperm="white")
    vb = api.tang_tf_validation("PFOS", dose="low")
    assert v["model_tf"]["endosperm"] == pytest.approx(vb["model_tf"]["endosperm"] * f["PFOS"])
    assert v["model_tf"]["stalk"] == vb["model_tf"]["stalk"]   # only the endosperm moves


def test_the_honesty_band_uses_the_default_configuration_on_tang():
    """The band's Tang entry must be the DEFAULT configuration's score, recomputed here, and
    every variant's Tang score must match the transfer record."""
    import json
    key = [k for k in api.APRIORI_OOS_SCORES if k.startswith("Tang")][0]
    default = _rmse()
    assert api.APRIORI_OOS_SCORES[key] == pytest.approx(default, abs=1e-4)
    with open(os.path.join(ROOT, "validation", "transfer_sourced.json")) as f:
        rec = json.load(f)["errors"]
    for v, kw in (("redistributed_fit", {}), ("redistributed_lipid", dict(lipid_loading=True)),
                  ("redistributed_free", dict(lipid_loading=False)),
                  ("4pool_free", dict(shoot="4pool", lipid_loading=False))):
        assert _rmse(**kw) == pytest.approx(float(np.sqrt(np.mean(np.square(rec[v]["tang"])))),
                                            abs=1e-6), v
    with pytest.raises(ValueError):
        api.tang_tf_validation("PFOA", shoot="two-pool")


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
