"""The Yamazaki observation table is built from the SI transcription, not edited by hand."""
import csv
import math
import os
import sys

import pytest

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "src"))
sys.path.insert(0, os.path.join(_ROOT, "docs", "literature_db", "tools"))

import model_api as api                 # noqa: E402
import yamazaki_obs_from_si as Y        # noqa: E402


def test_table_is_the_builders_output():
    rows = Y.build()
    with open(Y.OBS) as f:
        tab = list(csv.DictReader(line for line in f if not line.startswith("#")))
    assert len(tab) == len(rows) == 3 * len(Y.CONGENERS)
    for r, (c, t, mid, lo, hi, w, u) in zip(tab, rows):
        assert (r["compound"], r["tissue"], int(r["use"])) == (c, t, u)
        for key, v in (("baf_dw", mid), ("lo", lo), ("hi", hi), ("water_ng_L", w)):
            assert float(r[key]) == pytest.approx(v, rel=1e-5, abs=1e-12)


def test_exposure_is_table_s16_and_root_the_cultivar_geomean():
    si = Y.read_si()
    assert si["PFOA"]["water"] == 8.61 and si["PFBA"]["water"] == 7.41      # S16, p.S31
    d = api.yamazaki_obs("PFOA")
    assert d["water_ng_L"] == 8.61
    assert d["baf_dw"]["root"] == pytest.approx(math.sqrt(39 * 45) / 8.61, rel=1e-5)


def test_open_values_carry_their_ends():
    d = api.yamazaki_obs()
    for c, r in d.items():
        for t in ("root", "straw", "grain"):
            assert r["lo"][t] <= r["baf_dw"][t] <= r["hi"][t]
        assert r["lo"]["root"] == r["hi"]["root"]          # every root value is measured
    assert not d["PFBS"]["use"]["grain"] and d["PFBS"]["lo"]["grain"] == 0.0


def test_straw_band_shares_bracket_the_harvest():
    sh = Y.straw_shares()
    assert sorted(sh) == [10, 49]                           # days after heading
    for w in sh.values():
        assert sum(w) == pytest.approx(1.0)
        assert all(0.15 < x < 0.35 for x in w)
