"""Transcriptions of the second wave of rice-growth sources (docs/literature_db/raw_si/):
Yang et al. 2006 (Wuxi FACE), Gao et al. 2021 (cold-region japonica ORYZA2000 set + Fig. 3),
Sasaki et al. 2005 (the 1999 Shizukuishi FACE crop) and Hasegawa et al. 2017 Table S1.

Pinned: each table's own internal consistency -- the checks that caught transcription errors
elsewhere in this repo -- and the one cross-source fact the Sasaki paper adds: Kim 2003's 1999
"PI" and "anthesis" samplings were taken ~8 and ~10 days after the events themselves.
"""
import csv
import os
from collections import defaultdict

import pytest

RAW = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "docs", "literature_db", "raw_si")


def _rows(name):
    with open(os.path.join(RAW, name)) as f:
        return list(csv.DictReader(line for line in f if not line.startswith("#")))


def test_yang2006_fractions_sum_to_one_and_dm_grows():
    rows = _rows("yang2006_wuxi_face_dm.csv")
    assert len(rows) == 30                                         # 3 years x 2 CO2 x 5 stages
    by = defaultdict(list)
    for r in rows:
        tot = float(r["leaf_frac"]) + float(r["stem_frac"]) + float(r["spike_frac"] or 0)
        assert tot == pytest.approx(1.0, abs=0.0015)
        by[(r["year"], r["co2"])].append(float(r["aboveground_dm_g_m2"]))
        assert (r["spike_frac"] == "") == (int(r["stage_no"]) <= 2)   # no panicle before heading
    for dm in by.values():
        assert dm == sorted(dm)                                    # monotone accumulation
    face = {k: v for k, v in by.items() if k[1] == "FACE"}
    amb = {k: v for k, v in by.items() if k[1] == "AMB"}
    for (y, _), v in face.items():                                 # FACE > AMB at every stage (Table 1)
        assert all(f > a for f, a in zip(v, amb[(y, "AMB")]))


def test_gao2021_parameters_and_phenology():
    p = {r["parameter"]: r for r in _rows("gao2021_cold_japonica_oryza_params.csv")}
    assert (p["DVRJ"]["calibrated"], p["DVRR"]["calibrated"]) == ("0.000890", "0.002000")
    assert p["DVRJ"]["default_as_printed"] == "0.000773"           # the IR72.D92 value
    assert p["DRLVT"]["calibrated"].startswith("0.00:0.000;0.40:0.010")
    # partitioning tables: leaf + stem + panicle = 1 at every DVS node
    tabs = [dict(tuple(map(float, kv.split(":"))) for kv in p[k]["calibrated"].split(";"))
            for k in ("FLVTB", "FSTTB", "FSOTB")]
    for dvs in tabs[0]:
        assert sum(t[dvs] for t in tabs) == pytest.approx(1.0)
    ev = {k: int(p[k]["calibrated"]) for k in ("EMD", "IDOYTR", "IDOYPI", "IDOYFL", "IDOYM")}
    assert ev == {"EMD": 111, "IDOYTR": 136, "IDOYPI": 191, "IDOYFL": 214, "IDOYM": 261}
    ph = _rows("gao2021_qingan_phenology.csv")
    assert sum(int(r["days"]) for r in ph if r["year"] == "2018") == 125   # Table 2 total
    assert sum(int(r["days"]) for r in ph if r["year"] == "2019") == 128


def test_gao2021_fig3_organs_sum_to_green_biomass():
    rows = _rows("gao2021_fig3_digitised.csv")
    d = defaultdict(dict)
    for r in rows:
        if r["value"]:
            d[(r["year"], r["date"])][r["series"]] = float(r["value"])
    checked = 0
    for k, v in d.items():
        if all(s in v for s in ("D-AGB", "DWST", "DWGL", "DWSO")) and v["D-AGB"] > 1000:
            rel = (v["DWST"] + v["DWGL"] + v["DWSO"]) / v["D-AGB"] - 1
            if k == ("2018", "2018-07-26"):
                assert 0.10 < rel < 0.18                            # the one overlapped date, flagged
            else:
                assert -0.06 < rel < 0.0
            checked += 1
    assert checked == 12
    assert len({r["date"] for r in rows if r["year"] == "2018" and r["series"] == "LAI"}) == 9
    assert len({r["date"] for r in rows if r["year"] == "2019" and r["series"] == "LAI"}) == 10


def test_sasaki2005_dates_and_the_kim_sampling_lag():
    s = {r["item"]: r for r in _rows("sasaki2005_shizukuishi_1999.csv") if r["item"] != "lai_sampling"}
    assert int(s["panicle_initiation"]["dat"]) == 54 and int(s["heading_ambient"]["dat"]) == 78
    assert int(s["transplanting"]["doy_1999"]) - int(s["emergence"]["doy_1999"]) == 23
    lai = [int(r["dat"]) for r in _rows("sasaki2005_shizukuishi_1999.csv") if r["item"] == "lai_sampling"]
    kim = [int(r["DAT_or_date"].split("DAT")[1].split()[0]) for r in _rows("rice_organ_biomass_literature.csv")
           if r["source"] == "Kim2003_GCB_Table2_1999_AMB" and r["organ"].startswith("green leaf blade")]
    assert lai == kim == [32, 62, 88, 105, 123]                    # the same five samplings
    assert kim[1] - 54 == 8 and kim[2] - 78 == 10                  # "PI" / "anthesis" lag the events


def test_hasegawa2017_nursery_durations():
    import datetime as dt
    h = defaultdict(dict)
    for r in _rows("hasegawa2017_face_conditions.csv"):
        h[(r["site"], r["year"])][r["item"]] = r["value"]
    days = []
    for y in ("1998", "1999", "2000"):
        so = dt.date.fromisoformat(f"{y}-{h[('Shizukuishi', y)]['sowing']}")
        tp = dt.date.fromisoformat(f"{y}-{h[('Shizukuishi', y)]['transplanting']}")
        days.append((tp - so).days)
    assert days == [20, 27, 29]
    assert h[("Wuxi", "2001")]["Tmean_season"] == "25.2"
