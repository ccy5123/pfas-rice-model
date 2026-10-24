"""Re-fits on the FRESH organ-mass basis, each pinned next to its DRY-mass record.

`src/organ_mass.py` fixed a dry-matter-as-fresh-mass bug in every biomass driver. The
fits made before it were kept (the `*_dry` tables / caches, selected by
mass_basis='dry') and re-made on fresh mass. Each test below pins BOTH so neither the
record nor the re-fit can drift silently:

  * TANG_REFIT_FXY   (validation/tang2026_fxy_refit.py)   in-sample Tang 0.519 -> 0.327;
                      PFOS at the f_xy = 1 ceiling on fresh mass.
  * LIPID_LOADING    (validation/fit_lipid_loading.py)    Yamazaki excl PFDoDA 0.372 -> 0.177;
                      Tang OOS 0.516 -> 0.520 (the generalisation survives).
  * two-pool seq, measured forcings (validation/twopool_root_measured.py): the non-K_PL
                      PFOS/PFUnDA sink separation (4.5x on dry mass) does NOT carry over --
                      on fresh mass both sinks are ~0.
"""
import json
import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "validation"))

import model_api as api  # noqa: E402


def _tang_rmse(mass_basis, **kw):
    e = []
    for c in api.TANG_CONGENERS:
        v = api.tang_tf_validation(c, f_xy_source="recommended", dose="low",
                                   mass_basis=mass_basis, **kw)
        e += [np.log10(v["model_tf"][o]) - np.log10(v["tang_tf"][o]) for o in v["organs"]]
    return float(np.sqrt(np.mean(np.square(e))))


# --- TANG_REFIT_FXY --------------------------------------------------------------
def test_tang_refit_table_follows_the_mass_basis():
    assert api.tang_refit_fxy("dry") is api.TANG_REFIT_FXY_DRY
    assert api.tang_refit_fxy() is api.TANG_REFIT_FXY
    assert api.TANG_REFIT_FXY_DRY == {"PFOA": 0.097, "PFOS": 0.320, "GenX": 0.017}
    for mb, table in (("dry", api.TANG_REFIT_FXY_DRY), ("fresh", api.TANG_REFIT_FXY)):
        v = api.tang_tf_validation("PFOA", use_refit=True, mass_basis=mb)
        assert v["f_xy"] == table["PFOA"] and v["refit_fxy"] == table["PFOA"]
    # on fresh mass PFOS needs the free-anion f_xy at its ceiling (unrestricted DPU)
    assert api.TANG_REFIT_FXY["PFOS"] >= 0.99
    with pytest.raises(ValueError):
        api.tang_refit_fxy("wet")


def test_tang_refit_in_sample_rmse_on_both_bases():
    assert _tang_rmse("dry", use_refit=True) == pytest.approx(0.519, abs=0.01)
    assert _tang_rmse("fresh", use_refit=True) == pytest.approx(0.327, abs=0.02)


# --- LIPID_LOADING ---------------------------------------------------------------
def test_lipid_constants_follow_the_mass_basis():
    assert api.lipid_loading_params("dry") is api.LIPID_LOADING_DRY
    assert api.lipid_loading_params() is api.LIPID_LOADING
    assert api.LIPID_LOADING_DRY == dict(a=0.25, b=0.50, g_xy=0.05, g_ph=0.010,
                                         K_half=3000.0, pfsa_ln=1.25)
    d = api.simulate("PFDA", lipid_loading=True, mass_basis="dry")["params"]
    f = api.simulate("PFDA", lipid_loading=True)["params"]
    exp_d = api.lipid_loading_conductances(10, api._CONG["PFDA"]["K_PL_Lkg"], "PFCA", "dry")
    exp_f = api.lipid_loading_conductances(10, api._CONG["PFDA"]["K_PL_Lkg"], "PFCA")
    assert (d["f_xy"], d["g_xy"], d["g_ph"]) == pytest.approx(exp_d)
    assert (f["f_xy"], f["g_xy"], f["g_ph"]) == pytest.approx(exp_f)
    assert f["g_xy"] != pytest.approx(d["g_xy"])


def test_lipid_fit_yamazaki_on_both_bases():
    """The fit script's protocol reproduces the documented ~0.36 on dry mass with the
    stored dry constants, and the fresh re-fit reaches 0.177 -- root 0.488 -> ~0.09:
    the single-pool root/shoot trade-off was largely the dry-mass M."""
    import fit_lipid_loading as FL
    s_dry = FL.score(FL.predict_all(api.LIPID_LOADING_DRY, mass_basis="dry"))
    s_fresh = FL.score(FL.predict_all(api.LIPID_LOADING, mass_basis="fresh"))
    assert s_dry["excl"] == pytest.approx(0.372, abs=0.01)
    assert s_dry["root"] == pytest.approx(0.488, abs=0.02)
    assert s_fresh["excl"] == pytest.approx(0.177, abs=0.01)
    assert s_fresh["root"] < 0.12


def test_lipid_tang_oos_survives_the_refit():
    """Tang OOS (constants never fit on Tang): dry record 0.516, fresh re-fit 0.520 --
    both far below the free anion (1.232 dry / 1.483 fresh)."""
    assert _tang_rmse("dry", lipid_loading=True) == pytest.approx(0.516, abs=0.01)
    lip = _tang_rmse("fresh", lipid_loading=True)
    assert lip == pytest.approx(0.520, abs=0.02)
    assert lip < _tang_rmse("fresh") - 0.5


# --- two-pool seq on measured forcings -------------------------------------------
def test_twopool_measured_caches_per_basis():
    import twopool_root_exploration as TP
    here = os.path.join(ROOT, "validation")
    dry = json.load(open(os.path.join(here, "twopool_fitted_params_measured_dry.json")))
    fresh = json.load(open(os.path.join(here, "twopool_fitted_params_measured.json")))
    assert dry["mass_basis"] == "dry" and fresh["mass_basis"] == "fresh"
    qd, qf = np.array(dry["ushape_q"]), np.array(fresh["ushape_q"])
    # dry record: the non-K_PL separation at identical K_PL (PFOS C8 PFSA vs PFUnDA C11)
    ratio = TP.kseq_ushape(11, "PFCA", qd) / TP.kseq_ushape(8, "PFSA", qd)
    assert ratio == pytest.approx(4.5, abs=0.5)
    # fresh: neither sink is in use, so there is no separation to speak of
    assert TP.kseq_ushape(8, "PFSA", qf) < 1e-3
    assert TP.kseq_ushape(11, "PFCA", qf) < 1e-3


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
