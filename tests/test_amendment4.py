"""PREREG §11 (amendment 4): sourced whole-wall K_cw in parameters.json; Lu 2020 crop sets in the band."""
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import crop_sets as CS          # noqa: E402
import kcw_sourced as K         # noqa: E402

P = json.load(open(ROOT / "params" / "parameters.json"))
CONG = {c["name"]: c for c in P["congeners"]}


def _n_head(c):
    if c["group"] == "PFCA":
        return c["n_C"] - 1, "carboxylate", 0
    if c["group"] == "PFSA":
        return c["n_C"], "sulfonate", 0
    return 5, "carboxylate", 1                       # GenX: 5 perfluorinated C, 1 ether O


@pytest.mark.parametrize("name", sorted(CONG))
def test_parameters_carry_the_sourced_kcw(name):
    c = CONG[name]
    k = K.kcw_sourced(*_n_head(c)[:2], n_ether_O=_n_head(c)[2])
    vals = set(c["K_cw_wholecw_Lkg"].values())
    assert len(vals) == 1                            # one value, every organ
    assert vals.pop() == pytest.approx(k.rec, abs=1e-4)
    assert c["K_cw_case"] == k.case
    assert "K_cw_poly_rec_Lkg" not in c and "K_lignin_rec_Lkg" not in c


def test_cases_of_the_decision_tree():
    case = {n: CONG[n]["K_cw_case"] for n in CONG}
    assert case["PFHpA"] == case["GenX"] == "removed"
    assert {case[n] for n in ("PFBA", "PFPeA", "PFHxA")} == {"lower_only"}
    assert {case[n] for n in ("PFOA", "PFNA", "PFDA", "PFUnDA", "PFDoDA", "PFBS", "PFHxS", "PFOS")} == {"domain"}


def test_lu2020_sets_join_the_band():
    assert {f"lu2020_{s}" for s in ("es1", "es2", "es4", "es9")} <= set(CS.BAND_SETS)
    o = CS.overrides("lu2020_es9")
    assert o["SLATB"][1][-1] == pytest.approx(0.0025)     # ES9's own DVS-2 row (Table 2 layout, PREREG 11.8)
    assert o["DVRJ"] == pytest.approx(0.000545)
