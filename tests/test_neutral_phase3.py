"""Neutral / weak-electrolyte extension -- phase 3 (peripheral integration).

Phase 3 does not add model physics; it routes the EXISTING neutral parameterisation
through the three places that were still PFAS-only:

  1. soil sorption  -- Koc for the flooded/HYDRUS exposure shape came from the PFAS
     chain-length QSPR, so a neutral had no soil branch at all.
  2. the SMILES adapter -- neutral structures were flagged as "assumption violated"
     with nowhere to go.
  3. (app) a compound-class selector, covered by the UI smoke test.

So these tests check DISPATCH and DIRECTION, not new equations, plus the standing
invariant that the PFAS branch is untouched.

See docs/NEUTRAL_DPU_EXTENSION_DESIGN_KR.md section 6 and
docs/HANDOFF_neutral_extension.md section 5.
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

import literature_params as lp
import model_api as api
import soil_hydrus as sh


# --------------------------------------------------------------------------
# 1. soil Koc -- the neutral branch
# --------------------------------------------------------------------------
def test_paddy_kd_neutral_branch_is_karickhoff():
    """Kd = 0.41*Kow*f_oc, and n_C is ignored once logKow is given."""
    f_oc = 0.02
    assert sh.paddy_kd(logKow=3.0, f_oc=f_oc) == pytest.approx(0.41 * 1e3 * f_oc)
    # the neutral branch must not consult the chain-length QSPR at all
    assert sh.paddy_kd(4, "PFCA", f_oc, logKow=3.0) == sh.paddy_kd(12, "PFSA", f_oc, logKow=3.0)


def test_paddy_kd_still_needs_one_of_the_two_descriptors():
    with pytest.raises(TypeError):
        sh.paddy_kd()


def test_paddy_kd_pfas_branch_unchanged():
    """The PFAS default path must be bit-identical to the chain-length QSPR."""
    assert sh.paddy_kd(8, "PFCA") == lp.koc_to_KF(lp.koc(7, "carboxylate"), 0.02)


def test_neutral_k_leach_falls_with_hydrophobicity():
    """A polar neutral leaches; a hydrophobic one stays buffered (k_leach -> 0).

    This is the same physics that makes short-chain PFAS leach and long-chain PFAS
    buffer -- only the Koc source differs.
    """
    ks = [api.default_k_leach(logKow=x) for x in (1.0, 2.0, 3.0, 4.0, 5.0, 6.0)]
    assert all(a > b or (a == b == 0.0) for a, b in zip(ks, ks[1:])), ks
    assert ks[0] > 0.02 and ks[-1] == pytest.approx(0.0, abs=1e-9)
    assert all(0.0 <= k <= 0.15 for k in ks)          # clipped to the calibrated range


def test_default_k_leach_pfas_table_still_wins():
    """A curated congener keeps its HYDRUS-calibrated table value."""
    assert api.default_k_leach("PFOA") == api.default_k_leach("PFOA", logKow=3.0)


# --------------------------------------------------------------------------
# 2. the exposure shape for a neutral compound
# --------------------------------------------------------------------------
def test_flooded_profile_preserves_mean_exposure_for_a_neutral():
    t = np.linspace(0.0, 120.0, 241)
    for lk in (1.0, 3.0, 5.0):
        c = api.cwo_profile_series(t, level=2.5, profile="flooded", logKow=lk)
        assert np.mean(c) == pytest.approx(2.5, rel=1e-9)


def test_flooded_shape_direction_tracks_logKow():
    """Polar neutral -> steep decline; hydrophobic neutral -> ~flat (buffered)."""
    t = np.linspace(0.0, 120.0, 241)
    ratio = [api.cwo_profile_series(t, profile="flooded", logKow=lk)[-1] /
             api.cwo_profile_series(t, profile="flooded", logKow=lk)[0]
             for lk in (1.0, 3.0, 5.0)]
    assert ratio[0] < 0.3, ratio          # polar leaches away
    assert ratio[-1] > 0.95, ratio        # hydrophobic stays put
    assert ratio[0] < ratio[1] < ratio[2]


def test_simulate_neutral_accepts_the_exposure_shape():
    """Phase 3 acceptance: simulate_neutral(cwo_profile='flooded') runs, and
    'constant' is exactly the pre-phase-3 default."""
    base = api.simulate_neutral(3.0, Cwo=1.0)
    same = api.simulate_neutral(3.0, Cwo=1.0, cwo_profile="constant")
    assert base["straw_baf"] == same["straw_baf"]

    fl = api.simulate_neutral(3.0, Cwo=1.0, cwo_profile="flooded")
    assert fl["success"]
    assert np.mean(fl["Cwo"]) == pytest.approx(1.0, rel=1e-9)
    assert not np.allclose(fl["Cwo"], base["Cwo"])     # the shape really is different


def test_flooded_exposure_matters_more_for_a_polar_neutral():
    """A buffered (hydrophobic) compound barely notices the shape; a leaching one does.

    Guards the wiring: if logKow were NOT reaching the Koc branch both would shift by
    the same amount.
    """
    def shift(lk):
        a = api.simulate_neutral(lk, Cwo=1.0)["straw_baf"]
        b = api.simulate_neutral(lk, Cwo=1.0, cwo_profile="flooded")["straw_baf"]
        return abs(b - a) / a
    assert shift(1.0) > 10.0 * shift(5.0)


# --------------------------------------------------------------------------
# 3. the SMILES adapter -- compound-class dispatch
# --------------------------------------------------------------------------
PFOA_SMILES = "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"
ATRAZINE = "CCNc1nc(Cl)nc(NC(C)C)n1"
D24 = "OC(=O)COc1ccc(Cl)cc1Cl"          # 2,4-D: carboxylate head, NOT perfluorinated


def _ps():
    return pytest.importorskip("pfas_structure", reason="RDKit not installed")


def test_compound_class_separates_pfas_from_ordinary_acids():
    ps = _ps()
    assert ps.descriptors(PFOA_SMILES).compound_class == "PFAS"
    # 2,4-D is the discriminating case: 'carboxylate' head but zero perfluorinated C,
    # so a head-group test alone would wrongly call it a permanent anion.
    d = ps.descriptors(D24)
    assert d.head_group == "carboxylate" and d.n_perfluoroC == 0
    assert d.compound_class == "organic"
    assert ps.descriptors(ATRAZINE).compound_class == "organic"


def test_crippen_logkow_is_populated_and_ordered():
    ps = _ps()
    assert ps.descriptors(ATRAZINE).logKow_crippen < ps.descriptors(PFOA_SMILES).logKow_crippen


def test_neutral_smiles_builds_a_briggs_compound():
    ps = _ps()
    c, d = ps.neutral_compound_from_smiles(ATRAZINE)
    assert (c.fn, c.fd) == (1.0, 0.0)                  # strictly neutral without a pKa
    assert c.K_PL == 0.0 and c.K_prot == 0.0           # PFAS pools off
    assert c.K_lip > 0.0                                # Briggs lipid term carries binding
    assert c.logKow == pytest.approx(d.logKow_crippen)


def test_measured_logkow_overrides_crippen():
    ps = _ps()
    c, d = ps.neutral_compound_from_smiles(ATRAZINE, logKow=2.61)
    assert c.logKow == 2.61
    assert c.logKow != d.logKow_crippen          # really overrode, not coincidence
    assert any("user-supplied measured" in n for n in d.notes)


def test_pKa_makes_it_a_weak_electrolyte():
    ps = _ps()
    c, _ = ps.neutral_compound_from_smiles(D24, pKa=2.73)
    assert c.fd > 0.99 and c.fn < 0.01                 # fully dissociated at paddy pH
    assert c.kappa_d == pytest.approx(lp.PN_DEFAULT / lp.PN_OVER_PD)


def test_auto_dispatch_routes_each_class():
    ps = _ps()
    assert ps.compound_from_smiles_auto(PFOA_SMILES)[2] == "PFAS"
    assert ps.compound_from_smiles_auto(ATRAZINE)[2] == "organic"
    # an explicit pKa forces the weak-electrolyte branch even for a PFAS structure
    assert ps.compound_from_smiles_auto(PFOA_SMILES, pKa=0.5)[2] == "organic"


def test_simulate_from_smiles_dispatches_and_completes():
    """Phase 3 acceptance: SMILES -> parameters -> ODE finishes for both classes."""
    _ps()
    pfas = api.simulate_from_smiles(PFOA_SMILES)
    assert pfas["compound_class"] == "PFAS" and pfas["success"]

    neu = api.simulate_from_smiles(ATRAZINE)
    assert neu["compound_class"] == "organic" and neu["success"]
    assert neu["baf_final"]["root"] > 0.0


def test_provisional_flag_tracks_the_logkow_source():
    """A Crippen estimate is provisional; a measured logKow is not. logKow drives
    every neutral parameter, so its source is the thing that decides."""
    _ps()
    assert api.simulate_from_smiles(ATRAZINE)["provisional"] is True
    assert api.simulate_from_smiles(ATRAZINE, logKow=2.61)["provisional"] is False


def test_known_pfas_from_smiles_is_unaffected_by_the_dispatcher():
    """The PFAS branch must return exactly what it returned before phase 3."""
    _ps()
    assert (api.simulate_from_smiles(PFOA_SMILES)["baf_final"]["root"]
            == api.simulate("PFOA")["baf_final"]["root"])


# --------------------------------------------------------------------------
# 4. the app wiring (phase 3 item 3)
# --------------------------------------------------------------------------
# The UI itself is verified with headless Streamlit + Playwright (repo convention);
# what is worth locking in pytest is the CONTRACT between the sidebar and the runner,
# because a renamed cfg field fails silently at render time rather than at import.
def test_sidebar_exports_the_compound_class_contract():
    import inspect
    st = pytest.importorskip("streamlit")          # noqa: F841
    import ui.sidebar
    src = inspect.getsource(ui.sidebar.build)
    for fld in ("compound_class", "neutral_smiles", "logKow", "pKa",
                "is_acid", "K_AW", "compound_name"):
        assert f"cfg.{fld}" in src, f"sidebar no longer exports cfg.{fld}"


def test_simple_mode_stays_pfas_only():
    """The general-audience view must not reach the unvalidated neutral branch."""
    import inspect
    pytest.importorskip("streamlit")
    import ui.sidebar
    src = inspect.getsource(ui.sidebar.build)
    # the shared defaults block (used as-is by Simple mode) pins PFAS
    assert 'compound_class = "PFAS"' in src


def test_run_model_dispatches_on_compound_class():
    import inspect
    pytest.importorskip("streamlit")
    import ui.common
    src = inspect.getsource(ui.common.run_model)
    assert 'compound_class' in src and '_simulate_neutral' in src


def test_ui_modules_import():
    """Catches syntax/name errors in the UI package without a Streamlit runtime."""
    pytest.importorskip("streamlit")
    import ui.common, ui.expert, ui.i18n, ui.sidebar, ui.simple    # noqa: F401


def test_sulfonamide_now_has_somewhere_to_go():
    """BEHAVIOUR CHANGE: a perfluoroalkyl sulfonamide (FOSA) is not a permanent anion
    (pKa ~6), and before phase 3 it was flagged as violating the PFAS assumption while
    still being run on the PFAS branch. It now classifies as 'organic' and routes to
    the neutral / weak-electrolyte model, which is where that chemistry belongs.
    The flag itself is unchanged -- the note still says NOT a PERMANENT ANION."""
    ps = _ps()
    fosa = "NS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"
    d = ps.descriptors(fosa)
    assert d.head_group == "sulfonamide"
    assert d.compound_class == "organic"
    assert any("PERMANENT ANION" in n for n in d.notes)
    assert api.simulate_from_smiles(fosa)["compound_class"] == "organic"
    # the PFAS branch is still reachable explicitly, for comparison against the old path
    assert api.simulate_from_smiles(fosa, compound_class="PFAS")["compound_class"] == "PFAS"
