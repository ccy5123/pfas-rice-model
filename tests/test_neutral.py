"""Neutral / weak-electrolyte (DPU base) extension -- phase 1.

Two jobs:
  A. GUARD the PFAS path.  Every field the extension adds defaults to the value at
     which its term vanishes identically, so the anion results must not move at all.
  B. LOCK the structural claims of the neutral extension (speciation limits, the
     Briggs TSCF bell and where the tissue peak actually lands, root equilibrium,
     and the ion trap's kinetic -- not thermodynamic -- switch-off).

See docs/NEUTRAL_DPU_EXTENSION_DESIGN_KR.md.
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

import literature_params as lp
import model_api as api
from pfas_rice_plant_module_4pool_surf import (
    Compartment, Compound, Environment, binding_factors, root_uptake)

# Golden PFAS BAFs (root, stem, leaf, grain), captured from the code BEFORE the
# neutral extension landed and verified bit-identical after it.  These move only if
# the PFAS model itself changes.
GOLDEN_PFAS = {
    ("PFOA", "recommended"): (0.4757274987026183, 0.015321669867345053,
                              2.3412053079362147, 0.1438949489120824),
    ("PFOS", "W2fit"): (5.496261343207211, 0.47116670808650096,
                        7.0621989630642865, 0.4634163827272172),
    ("PFBA", "oryza"): (0.37987808562007025, 0.15231292788197617,
                        48.030924423109205, 2.5034485779765236),
}

_UPTAKE_KW = dict(K_prot=0.0, K_PL=0.0, K_cw=0.0, Km_in=1.0, Km_out=1.0,
                  L_Ph=1.0, f_xy=1.0)


# --------------------------------------------------------------------------
# A. PFAS regression guards
# --------------------------------------------------------------------------
@pytest.mark.parametrize("key", sorted(GOLDEN_PFAS))
def test_pfas_results_are_unchanged_to_the_last_bit(key):
    """The extension must not perturb the anion path at all."""
    congener, source = key
    r = api.simulate(congener, f_xy_source=source, biomass="growth_rice", n_t=121)
    got = tuple(r["baf_final"][k] for k in ("root", "stem", "leaf", "grain"))
    assert got == GOLDEN_PFAS[key], f"{congener}/{source} moved: {got}"


def test_pfas_ignores_the_neutral_permeability():
    """f_n = 0 must kill the neutral branch, whatever P_n is."""
    js = [root_uptake(2.0, 0.5,
                      Compound("pfas", kappa_d=0.5, Vmax_in=20.0, Vmax_out=8.0,
                               fn=0.0, fd=1.0, P_n=P_n, **_UPTAKE_KW),
                      Environment())
          for P_n in (0.0, 1.0, 1e6)]
    assert js[0] == js[1] == js[2]


def test_binding_is_unchanged_when_the_lipid_term_is_off():
    """f_lip and K_lip both default to 0, so B_k must be identical."""
    comps = [Compartment("root", 0.90, 0.07, 0.015, 0.50),
             Compartment("stem", 0.83, 0.05, 0.005, 0.72)]
    cmpd = Compound("x", K_prot=50.0, K_PL=100.0, K_cw=7.0, kappa_d=0.5,
                    Vmax_in=20.0, Km_in=5.0, Vmax_out=8.0, Km_out=5.0, L_Ph=0.005)
    expected = np.array([c.theta + (1 - c.theta) * (c.f_prot * 50.0 + c.f_PL * 100.0
                                                    + c.f_cw * 7.0) for c in comps])
    assert (binding_factors(comps, cmpd) == expected).all()


def test_environment_N_for_falls_back_to_environment_valence():
    env = Environment()
    assert env.N_for(None) == env.N
    cmpd = Compound("x", kappa_d=0.5, Vmax_in=0.0, Vmax_out=0.0, **_UPTAKE_KW)
    assert env.N_for(cmpd) == env.N                      # Compound.z defaults to None


# --------------------------------------------------------------------------
# B. Structural claims of the neutral extension
# --------------------------------------------------------------------------
def test_neutral_does_not_feel_the_membrane_potential():
    """f_d = 0 kills the GHK branch, so the valence becomes irrelevant (D7).

    This is why the speciation switch is (f_n, f_d) and not z: a weak acid is a
    neutral molecule AND an anion at once, which one global valence cannot express.
    """
    env = Environment()
    js = [root_uptake(2.0, 0.5,
                      Compound("neutral", kappa_d=0.5, Vmax_in=0.0, Vmax_out=0.0,
                               fn=1.0, fd=0.0, P_n=10.0, z=z, **_UPTAKE_KW), env)
          for z in (-2, -1, 0, 1, 2)]
    assert all(j == js[0] for j in js)
    assert js[0] == pytest.approx(10.0 * (2.0 - 0.5))    # plain Fickian


def test_speciation_covers_the_whole_spectrum():
    fn, fd = lp.speciation(-3.0, 6.5)                     # PFSA-like
    assert fd > 0.999999 and fn < 1e-9
    fn, fd = lp.speciation(10.0, 6.5)                     # essentially neutral
    assert fn > 0.999
    fn, fd = lp.speciation(6.5, 6.5)                       # pKa == pH
    assert fn == pytest.approx(0.5) and fd == pytest.approx(0.5)
    # a weak BASE is protonated (cation) at low pH -- the opposite direction
    fn_a, _ = lp.speciation(8.0, 6.5, is_acid=True)
    fn_b, _ = lp.speciation(8.0, 6.5, is_acid=False)
    assert fn_a > 0.9 and fn_b < 0.1


def test_ion_trap_strong_acid_limit_is_the_pH_ratio_not_unity():
    """The trap switches off KINETICALLY for a permanent anion, not thermodynamically.

    The equilibrium factor tends to 10**(dpH), NOT 1 -- it is derived assuming the
    neutral species carries the transport.  What actually vanishes for PFAS is the
    permeability-weighted neutral fraction that would establish it.
    """
    assert lp.ion_trap_factor(-3.0, 7.2, 8.0) == pytest.approx(10.0 ** 0.8, rel=1e-3)
    assert lp.ion_trap_factor(4.0, 7.2, 8.0) == pytest.approx(6.3, abs=0.1)
    assert lp.neutral_pathway_ratio(4.0, 7.2) > 1.0        # herbicide: neutral carries it
    assert lp.neutral_pathway_ratio(-3.0, 7.2) < 1e-5      # PFAS: pathway is dead
    assert (lp.neutral_pathway_ratio(4.0, 7.2)
            / lp.neutral_pathway_ratio(-3.0, 7.2)) > 1e6


def test_briggs_relationships():
    assert lp.briggs_tscf(lp.TSCF_LOGKOW_OPT) == pytest.approx(lp.TSCF_MAX)
    assert lp.briggs_tscf(1.0) < lp.TSCF_MAX and lp.briggs_tscf(5.0) < 0.02
    assert lp.briggs_klip(0.0) == pytest.approx(lp.BRIGGS_A)
    # the fresh-weight "octanol equivalent" is NOT a measured dry-weight fraction
    assert lp.f_lip_from_fresh_weight(0.025, 0.90) == pytest.approx(0.25)


def test_neutral_root_reaches_the_briggs_equilibrium():
    """With the default fast-exchange P_n the root approaches C = B * C_w^o."""
    for logKow in (1.0, 3.0, 5.0):
        r = api.simulate_neutral(logKow, n_t=121)
        assert r["success"]
        assert r["baf_final"]["root"] / r["B_k"]["root"] > 0.95


def test_neutral_kinetic_limit_is_below_equilibrium():
    """P_n is a RATE: lowering it delays the root without changing its target."""
    fast = api.simulate_neutral(2.0, n_t=121)
    slow = api.simulate_neutral(2.0, P_n=1.0, n_t=121)
    assert slow["B_k"]["root"] == fast["B_k"]["root"]          # same equilibrium
    assert slow["baf_final"]["root"] < 0.5 * fast["baf_final"]["root"]


def test_tissue_concentration_peak_is_right_of_the_tscf_bell():
    """conc ~ TSCF * B and B rises with K_ow, so the tissue peak shifts RIGHT.

    Guards against wiring up TSCF but forgetting the lipid binding term, which
    would put the straw peak back at the TSCF optimum of 1.78.
    """
    grid = np.arange(1.0, 5.01, 0.5)
    straw = [api.simulate_neutral(float(x), n_t=121)["straw_baf"] for x in grid]
    peak = float(grid[int(np.argmax(straw))])
    assert 3.0 <= peak <= 4.5, f"straw peak at {peak}, expected right of the bell"
    assert peak > lp.TSCF_LOGKOW_OPT + 1.0


def test_weak_acid_runs_both_membrane_pathways():
    r = api.simulate_neutral(2.0, pKa=4.0, n_t=121)
    assert r["success"]
    assert 0.0 < r["params"]["fn"] < 1.0
    assert r["params"]["fn"] + r["params"]["fd"] == pytest.approx(1.0)
    assert r["params"]["z"] == -1                          # acid -> anion
    assert api.simulate_neutral(2.0, pKa=8.0, is_acid=False, n_t=121)["params"]["z"] == +1


def test_simulate_neutral_matches_the_simulate_contract():
    """Same dict shape as simulate(), so plots/exports work unchanged."""
    neutral = api.simulate_neutral(2.0, n_t=121)
    pfas = api.simulate("PFOA", n_t=121)
    shared = set(pfas) - {"params"}
    assert shared <= set(neutral)
    for key in ("conc", "baf", "baf_final", "B_k"):
        assert set(neutral[key]) == set(pfas[key]) == set(api.TISSUES)
    assert neutral["briggs"]["tscf"] == pytest.approx(lp.briggs_tscf(2.0))


def test_metabolism_hook_reduces_tissue_concentration():
    """gamma was always in the model but pinned at 0 for recalcitrant PFAS."""
    base = api.simulate_neutral(2.0, n_t=121)
    meta = api.simulate_neutral(2.0, gamma=0.05, n_t=121)
    assert meta["baf_final"]["root"] < base["baf_final"]["root"]
    assert meta["straw_baf"] < base["straw_baf"]
