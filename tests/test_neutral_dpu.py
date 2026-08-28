"""
Tests for the NEUTRAL-organic (Briggs/Kow) path -- src/neutral_dpu.py.

The neutral path's whole value is that nothing in it is fitted: partitioning and
root->shoot loading both follow from log Kow via published QSPRs. These tests lock
that in -- the adapter really evaluates Briggs, the ionic machinery really is off,
and the emergent Kow ordering is the one the literature reports -- plus the scope
limits that must not be forgotten (no phloem, no air exchange, unbounded leaf
without metabolism).
"""
import os
import sys

import numpy as np
import pytest

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "src"))
sys.path.insert(0, os.path.join(_ROOT, "validation"))

from pfas_rice_plant_module_4pool_surf import binding_factors  # noqa: E402
import neutral_dpu as ND  # noqa: E402


@pytest.fixture(scope="module")
def drv():
    import neutral_dpu_validation as V
    return V.drivers()


# --------------------------------------------------------------------------
# the published QSPRs
# --------------------------------------------------------------------------
def test_briggs_tscf_is_the_published_bell():
    """Peak 0.784 at log Kow 1.78, falling away symmetrically in log Kow."""
    assert ND.briggs_tscf(1.78) == pytest.approx(0.784, rel=1e-9)
    grid = np.linspace(-2, 7, 400)
    vals = [ND.briggs_tscf(x) for x in grid]
    assert grid[int(np.argmax(vals))] == pytest.approx(1.78, abs=0.03)
    assert ND.briggs_tscf(1.78 - 2.0) == pytest.approx(ND.briggs_tscf(1.78 + 2.0), rel=1e-9)
    assert ND.briggs_tscf(7.0) < 0.01 and ND.briggs_tscf(-2.0) < 0.01


def test_partition_adapter_reproduces_briggs_rcf_exactly():
    """The core's basis-A binding, given the neutral composition, IS Briggs' RCF.

    This is what makes the neutral path a use of the published partition core
    rather than a re-fit of it -- so it must hold to machine precision, not
    approximately.
    """
    root = ND.briggs_root_compartment()
    for lk in (-1.0, 0.0, 1.78, 3.0, 6.0):
        c = ND.neutral_compound(ND.NeutralCompound("probe", lk))
        assert float(binding_factors([root], c)[0]) == pytest.approx(
            ND.briggs_rcf(lk), rel=1e-12)


def test_k_pw_is_water_plus_lipid_term():
    assert ND.k_pw(2.0, W=0.8, L=0.0) == pytest.approx(0.8)
    assert ND.k_pw(3.0, W=0.8, L=0.02, a=1.0, b=1.0) == pytest.approx(0.8 + 0.02 * 1000)


# --------------------------------------------------------------------------
# the ionic machinery is off -- by physics, not by special-casing
# --------------------------------------------------------------------------
def test_zero_valence_removes_exclusion_and_carrier():
    env = ND.neutral_environment()
    assert env.N == pytest.approx(0.0, abs=1e-15)
    assert np.exp(env.N) == pytest.approx(1.0, abs=1e-15)   # PFAS: e^N ~ 107
    c = ND.neutral_compound(ND.NeutralCompound("probe", 2.0))
    assert c.Vmax_in == 0.0 and c.Vmax_out == 0.0
    assert c.K_prot == 0.0 and c.K_cw == 0.0                # Briggs uses lipid only
    assert c.f_xy == pytest.approx(ND.briggs_tscf(2.0))     # TSCF, not a fitted f_xy


def test_membrane_term_degenerates_to_passive_diffusion():
    """With z=0 the GHK term must be exactly kappa_d*(Cwo - Cw_root)."""
    from pfas_rice_plant_module_4pool_surf import root_uptake
    c = ND.neutral_compound(ND.NeutralCompound("probe", 2.0, kappa_d=3.0))
    env = ND.neutral_environment()
    for cwo, cw in ((1.0, 0.0), (1.0, 0.4), (0.2, 0.9)):
        assert root_uptake(cwo, cw, c, env) == pytest.approx(3.0 * (cwo - cw), rel=1e-12)


# --------------------------------------------------------------------------
# emergent behaviour
# --------------------------------------------------------------------------
def test_kow_ordering_signature(drv):
    """Zero fitted parameters, yet the run must reproduce the law every uptake
    study reports: polar -> shoot-dominated, lipophilic -> root-retained, with the
    turnover set by the Briggs bell."""
    tf, root, lks = {}, {}, (0.0, 1.78, 3.5, 5.5)
    for lk in lks:
        r = ND.simulate_neutral(ND.NeutralCompound(f"p{lk}", lk), drv)
        tf[lk] = r["straw_baf"] / r["baf_final"]["root"]
        root[lk] = r["baf_final"]["root"]
    assert tf[1.78] > tf[0.0]                    # rising limb of the bell
    assert tf[1.78] > tf[3.5] > tf[5.5]          # falling limb
    assert tf[5.5] < 1.0                         # lipophilic: root dominates
    assert tf[0.0] > 1.0                         # polar: shoot dominates
    assert root[5.5] > root[3.5] > root[0.0]     # root loading rises with K_PW


def test_leaf_is_unbounded_without_metabolism(drv):
    """Scope limit that must not be forgotten: with no phloem, no air exchange and
    gamma=0, the leaf is a terminal accumulator fed by the whole transpiration
    stream. A realistic half-life must bound it."""
    fast = ND.rice_compartments(gammas=dict.fromkeys(
        ("root", "stem", "leaf", "grain"), np.log(2) / 7.0))
    free = ND.simulate_neutral(ND.NeutralCompound("p", 1.78), drv)
    met = ND.simulate_neutral(ND.NeutralCompound("p", 1.78), drv, comps=fast)
    # a 7-day half-life cuts the leaf by ~10x (194 -> 19); the root, which is
    # exposure-buffered rather than accumulating, barely moves
    assert free["baf_final"]["leaf"] > 5 * met["baf_final"]["leaf"]
    assert met["baf_final"]["root"] == pytest.approx(free["baf_final"]["root"], rel=0.1)


def test_phloem_is_off_by_default(drv):
    """The neutral base explicitly excludes phloem transport; turning it on is an
    opt-in departure and must visibly change the grain."""
    base = ND.simulate_neutral(ND.NeutralCompound("p", 2.0), drv)
    ph = ND.simulate_neutral(ND.NeutralCompound("p", 2.0), drv, phloem=True)
    assert base["phloem"] is False and ph["phloem"] is True
    assert ph["baf_final"]["grain"] != pytest.approx(base["baf_final"]["grain"], rel=1e-3)


def test_volatile_compound_is_flagged():
    """Air exchange is not implemented -- a volatile compound must not run silently."""
    assert ND.k_aw_warning(ND.NeutralCompound("nonvolatile", 2.0, K_AW=1e-7)) is None
    w = ND.k_aw_warning(ND.NeutralCompound("volatile", 2.0, K_AW=0.1))
    assert w and "NOT modelled" in w


def test_bad_water_content_rejected():
    with pytest.raises(ValueError):
        ND.neutral_compartment("root", W=1.0, L=0.01)


def test_schriever_refit_is_the_same_form_and_reproduces_briggs_paper():
    """Schriever & Lamshoeft 2020 refit the SAME Gaussian to 97 modern TSCF values
    and reprint Briggs' equation verbatim as their eq. 4 -- so the two must agree
    on peak height while the refit is much broader."""
    assert ND.schriever_tscf(ND.SCHRIEVER_TSCF["B"]) == pytest.approx(0.746, rel=1e-9)
    # near-identical maximum, far broader bell
    assert abs(ND.SCHRIEVER_TSCF["A"] - ND.TSCF_MAX) < 0.05
    assert ND.SCHRIEVER_TSCF["C"] > 2 * ND.TSCF_WIDTH
    # so a lipophilic compound is predicted to translocate far more under the refit
    assert ND.tscf(4.4, "schriever") > 5 * ND.tscf(4.4, "briggs")
    with pytest.raises(ValueError):
        ND.tscf(2.0, "nope")


def test_lipid_defaults_are_the_cited_trapp1994_values():
    """The lipid contents must stay traceable to Trapp 1994 (root 1%, stem/leaf 3%
    dry weight), not drift back to unsourced guesses."""
    assert ND.TRAPP1994_LIPID_DW["root"] == 0.01
    assert ND.TRAPP1994_LIPID_DW["stem"] == ND.TRAPP1994_LIPID_DW["leaf"] == 0.03
    comps = {c.name: c for c in ND.rice_compartments()}
    # stored as fresh-weight: f_PL*(1-theta) == dw fraction * (1-theta)
    fw_root = comps["root"].f_PL * (1.0 - comps["root"].theta)
    assert fw_root == pytest.approx(0.01 * (1.0 - ND.RICE_WATER["root"]), rel=1e-9)


# --------------------------------------------------------------------------
# the measured-data harness
# --------------------------------------------------------------------------
def test_ge2017_apriori_prediction():
    """The repo's first genuine A-PRIORI prediction: Ge et al. 2017 per-organ
    transfer factors for three neutral pesticides spanning log Kow -0.13 to 4.4,
    with NOTHING fitted (K_PW and TSCF both follow from log Kow alone).

    Pins the headline (log10 RMSE ~1.10) and the error STRUCTURE that interprets
    it: the stem is predicted well, the leaf is over-predicted, because the strict
    a-priori run has no in-planta metabolism and the leaf is therefore an
    unbounded terminal accumulator (validation section 3)."""
    import neutral_dpu_validation as V
    path = os.path.join(_ROOT, "data_obs", "neutral_obs_ge2017.csv")
    drv = V.drivers()
    rmse = V.compare_to_obs(path, drv, quiet=True)
    assert rmse == pytest.approx(1.099, abs=0.05)
    # the error is dominated by the MISSING half-life, not by transport structure:
    # imposing a realistic in-planta dissipation must reduce it monotonically
    errs = [V.compare_to_obs(path, drv, half_life=h, quiet=True)
            for h in (None, 30.0, 14.0, 7.0, 3.0)]
    assert errs == sorted(errs, reverse=True)
    assert errs[-1] < 0.55
    # and the ORIGINAL Briggs bell beats the broader modern refit on this data
    assert V.compare_to_obs(path, drv, tscf_model="schriever", quiet=True) > rmse


def test_obs_template_carries_no_data():
    """The shipped template must be schema-only: its placeholder rows are refused
    so the template can never be mistaken for measurements."""
    import neutral_dpu_validation as V
    path = os.path.join(_ROOT, "data_obs", "neutral_obs_template.csv")
    assert V.load_neutral_obs(path) == []


def test_obs_loader_reads_and_converts(tmp_path):
    import neutral_dpu_validation as V
    p = tmp_path / "obs.csv"
    p.write_text(
        "compound,log_kow,tissue,value,basis,endpoint\n"
        "probe,2.0,root,1.5,fw,baf\n"
        "probe,2.0,leaf,3.0,dw,baf\n")
    rows = V.load_neutral_obs(str(p))
    assert len(rows) == 2 and rows[0]["log_kow"] == 2.0 and rows[1]["basis"] == "dw"
    with pytest.raises(ValueError):
        bad = tmp_path / "bad.csv"
        bad.write_text("compound,log_kow,tissue\nprobe,2.0,root\n")
        V.load_neutral_obs(str(bad))


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
