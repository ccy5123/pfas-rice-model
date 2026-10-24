"""Three earlier figure digitisations re-checked against the article bitmaps
(docs/literature_db/tools/fig_overlay_checks.py; the PDFs are not in the repo, so these tests pin the
data the check left behind, not the image processing).

Kim 2003 Fig. 2 (ambient-CO2 total biomass): re-read values, the missing 1999 HN point, and the hidden
markers read from their crescents. Li 2025 Fig. 2a (irrigation-water composition): the shares close
against SI Table S5 class sums on the four clean stages, reproduce the printed PFOA/PFBA ranges, and
the verified water sits BESIDE the old values in data_obs. Tang 2026 Fig. 1: the control bars.
"""
import csv
import os
import re
from collections import defaultdict

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "docs", "literature_db", "raw_si")


def _rows(path):
    with open(path) as f:
        return list(csv.DictReader(line for line in f if not line.startswith("#")))


def test_kim2003_fig2_values_after_the_overlay_check():
    tot = {}
    for r in _rows(os.path.join(RAW, "rice_organ_biomass_literature.csv")):
        m = re.match(r"Kim2003_GCB_Fig2_(\d{4})_AMB", r["source"])
        if m:
            tot[(int(m.group(1)), r["N_level"][:2], int(re.search(r"DAT(\d+)", r["DAT_or_date"]).group(1)))] = float(r["value"])
            assert "overlay-verified" in r["basis"]
    assert len(tot) == 28 and (1999, "HN", 89) in tot                      # the point the first pass missed
    assert tot[(1998, "MN", 25)] == 23 and tot[(1999, "MN", 32)] == 47 and tot[(2000, "MN", 29)] == 75
    assert tot[(2000, "MN", 100)] == 1440 and tot[(1999, "LN", 61)] == 394
    for (y, nl), dats in {(1998, "MN"): (25, 53, 81, 109, 131), (1999, "MN"): (32, 62, 88, 105, 123),
                          (2000, "MN"): (29, 56, 78, 100, 120)}.items():
        v = [tot[(y, nl, d)] for d in dats]
        assert v == sorted(v)                                              # monotone accumulation
    # Kim 2001 Table 1 (1999, tabulated): the re-read figure runs 1-6 % above it, never below
    for (nl, d, table) in (("LN", 61, 372.07), ("MN", 62, 425.94), ("HN", 63, 452.53),
                           ("LN", 124, 1375.89), ("MN", 123, 1495.64), ("HN", 124, 1543.54)):
        assert 1.0 < tot[(1999, nl, d)] / table < 1.065


def test_li2025_fig2a_shares_close_against_the_tabulated_class_sums():
    sh = defaultdict(dict)
    for r in _rows(os.path.join(RAW, "li2025_fig2a_water_shares.csv")):
        sh[r["stage"]][r["congener"]] = (float(r["share_pct"]), float(r["conc_ngL"]), float(r["stage_sum_pfas_ngL"]))
    cls = defaultdict(dict)
    for r in _rows(os.path.join(RAW, "li2025_stage_water_soil.csv")):
        if r["medium"].startswith("irrigation"):
            cls[r["stage"]][r["class"]] = float(r["mean"])
    for st in ("tillering", "heading", "filling", "milk", "ripening"):
        t = sum(cls[st].values())
        pfca = sum(sh[st][c][0] for c in ("PFBA", "PFHxA", "PFHpA", "PFOA"))
        pfsa = sum(sh[st][c][0] for c in ("PFEtS", "PFBS", "PFOS"))
        if st == "ripening":                                               # the figure's own anomaly (6:2 FTS 4.45 %)
            assert 100 * cls[st]["PFCAs"] / t - pfca == pytest.approx(4.04, abs=0.05)
        else:
            assert abs(pfca - 100 * cls[st]["PFCAs"] / t) < 0.4
            assert abs(pfsa - 100 * cls[st]["PFSAs"] / t) < 0.25
        assert sh[st]["PFOA"][2] == pytest.approx(t, abs=0.01)
    pfoa = [sh[s]["PFOA"][0] for s in sh]; pfba = [sh[s]["PFBA"][0] for s in sh]
    assert min(pfoa) == pytest.approx(59.3, abs=0.15) and max(pfoa) == pytest.approx(70.2, abs=0.1)   # text p.3
    assert min(pfba) == pytest.approx(18.6, abs=0.1) and max(pfba) == pytest.approx(31.1, abs=0.1)
    mean = {c: sum(sh[s][c][1] for s in sh) / 5 for c in ("PFBA", "PFHxA", "PFOA", "PFBS", "PFOS")}
    assert mean["PFHxA"] == pytest.approx(1.21, abs=0.01) and mean["PFOS"] == pytest.approx(0.71, abs=0.01)
    # data_obs keeps the old record and carries the verified values beside it
    old = {r["congener"]: r for r in _rows(os.path.join(ROOT, "data_obs", "Li2025_BAF_TF.csv"))}
    for c, r in old.items():
        assert float(r["water_ngL_fig2a_verified"]) == pytest.approx(mean[c], abs=0.006)
        ratio = float(r["water_ngL"]) / float(r["water_ngL_fig2a_verified"])
        assert float(r["BAF_root_fw_verified"]) == pytest.approx(float(r["BAF_root_fw"]) * ratio, rel=1e-2)   # water rounded to 0.01
    assert old["PFHxA"]["water_ngL"] == "2.5" and old["PFOS"]["water_ngL"] == "0.4"   # record untouched
    for r in _rows(os.path.join(ROOT, "data_obs", "obs_baf_Li2025.csv")):
        ratio = float(old[r["compound"]]["water_ngL"]) / float(r["water_ngL_fig2a_verified"])
        assert float(r["baf_verified"]) == pytest.approx(float(r["baf"]) * ratio, rel=1e-2)


def test_tang2026_fig1_control_bars():
    v = {r["quantity"]: float(r["value_gfw_per_plant"]) for r in _rows(os.path.join(RAW, "tang2026_harvest_biomass.csv"))}
    assert v["whole_plant_biomass_control"] == 33.47 and v["ear_panicle_control"] == 6.86
    assert v["harvest_index_fw"] == pytest.approx(6.86 / 33.47, abs=0.001)
