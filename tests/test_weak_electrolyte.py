"""Guards for the weak-electrolyte speciation ported from PR #54.

WHAT WAS PORTED AND WHY. Two independent neutral-organic implementations existed:
this repo's `src/neutral_dpu.py` (a separate module reaching a neutral compound by
`z=0`, carrying all five measured validation tables) and PR #54 (extending the PFAS
core in place with an `(fn, fd)` speciation pair). They overlapped on air exchange,
`simulate_neutral` and the Briggs lipid partition -- but #54 could do one thing the
merged path could not: a WEAK ELECTROLYTE, which is a neutral molecule and an ion at
the same time and so cannot be described by a single valence. Only that capability
was ported; the duplicated parts were dropped in favour of what is already here.

The tests below are ordered by what would hurt most if it broke:

  1. the PFAS path is unmoved -- the whole change is worthless if not;
  2. the NEUTRAL path is unmoved, so every published a-priori number
     (Liu 0.281, Ge 0.783, the stem 0.299) still describes what the code does;
  3. the two modes are CONTINUOUS at their boundary, which is the property that
     makes `pKa=` safe to add rather than a second, subtly different model;
  4. the physics the port exists for: acid/base asymmetry, and the ion trap's
     kinetic (not thermodynamic) switch-off for a permanent anion.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

import literature_params as LP  # noqa: E402
import model_api as api  # noqa: E402
import neutral_dpu as ND  # noqa: E402
import pfas_rice_plant_module_4pool_surf as P4  # noqa: E402
from model_api import PHI, T_C_PH   # the sourced phloem constants  # noqa: E402


# --------------------------------------------------------------------------
# 1-2. nothing that already worked may move
# --------------------------------------------------------------------------
# Recorded from the pre-port tree. The speciation terms are written so they vanish
# IDENTICALLY on the PFAS path (P_n defaults to 0 and `0.0 + x == x`), so a leaked
# term is what this guards against.
#
# WHY THIS IS NOT EXACT `==`, WHICH IS HOW IT WAS FIRST WRITTEN. These are outputs
# of an ADAPTIVE STIFF ODE SOLVE, and those are not reproducible across machines:
# BLAS/LAPACK build, CPU instruction set (FMA contraction) and SciPy version all
# move the result. Putting the suite in CI proved it on the first run -- the same
# commit was built on two GitHub runners and passed on one, failed on the other.
#
# WHERE THE TOLERANCE COMES FROM, and it is NOT from observed samples. The first
# attempt set rel=1e-6 by measuring one pair of runners (~1e-8 apart) -- and CI
# then produced a 4.6e-6 disagreement on the very next commit, because two samples
# do not bound a distribution. This is the same mistake as fitting a claim to a
# subsample; see docs/HANDOFF_neutral_next.md section 2.
#
# The principled bound is the SOLVER's own tolerance. `solve_ivp(..., method="BDF",
# rtol=1e-6, atol=1e-9)` defines the trajectory only to ~1e-6 relative per step, so
# two runs whose arithmetic differs may take different step sequences and land
# anywhere within the accumulated error of a 120-day integration -- necessarily
# LOOSER than 1e-6, not tighter. rel=1e-3 sits three orders above the solver's rtol
# and still an order BELOW the smallest leak that could matter: a speciation term
# that failed to vanish moves these numbers by percent, and even the smallest
# bypass conductance on the g_apo grid (0.02) moves them by more than 1e-3.
#
# THE REAL BIT-EXACT GUARD IS ELSEWHERE, and that is the point. Golden constants
# compared across machines can never be exact; an algebraic identity checked in one
# process can. `test_root_uptake_neutral_term_is_identically_absent_by_default`
# reconstructs j_R without the added terms and demands exact `==` -- that is what
# actually proves nothing leaked. These constants are a coarse backstop, and the
# tolerance is set so they cannot become a source of false alarms.
#
# NOT ALL rel=1e-9 IN THIS SUITE IS SAFE, and the first version of this comment
# said it was ("they compare two results in the same process -- deterministic").
# Too broad: same process is not the criterion, same SOLVE is. Two runs that build
# the same Compound by different routes are exact; two runs that differ in one
# state and are compared because the physics says another state cannot move are
# not, because the coupled system takes different adaptive steps. CI failed on
# exactly those two assertions the run after this file was fixed. The three cases
# are set out in tests/conftest.py, which now owns that rule.
GOLDEN_REL = 1e-3

# Goldens on two fitted variants (calibration site-year, fresh organ mass, harvest at the mean
# date) and on the DEFAULT configuration (redistributed_fit), regenerated for amendment 2a
# (docs/PREREG_sourced_values.md sections 7-8: staged organ water, brown-rice grain, the
# terminal pools as burden, re-fit). Regenerate them only when the model itself is deliberately changed.
PFAS_GOLDEN = {
    ("PFOA", "4pool_free"): {"root": 1.7687909049187536, "grain": 0.20021364082263676},
    ("PFOS", "redistributed_fit"): {"root": 21.569014721901613, "grain": 1.8083563015803177},
    ("PFBA", "4pool_free"): {"root": 2.4728430001475936, "grain": 0.9640453467296496},
}
_VARIANT_KW = {"4pool_free": dict(api.FOUR_POOL_FREE),
               "redistributed_fit": dict(f_xy_source="fit", lipid_loading=False)}


@pytest.mark.parametrize("key,expected", list(PFAS_GOLDEN.items()))
def test_pfas_path_is_unmoved(key, expected):
    """The speciation terms must not move the PFAS path (a permanent anion, pKa None)."""
    cong, v = key
    r = api.simulate(cong, **_VARIANT_KW[v])
    for tissue, want in expected.items():
        assert r["baf_final"][tissue] == pytest.approx(want, rel=GOLDEN_REL), tissue


PFAS_GOLDEN_DEFAULT = {
    "PFOA": {"root": 1.7689644924094703, "grain": 0.424960408461994},
    "PFBA": {"root": 2.4730429003189482, "grain": 0.9653576963292074},
}


@pytest.mark.parametrize("cong,expected", list(PFAS_GOLDEN_DEFAULT.items()))
def test_pfas_path_is_unmoved_in_the_default(cong, expected):
    """The same backstop on the default configuration (redistributed shoot, per-congener f_xy)."""
    r = api.simulate(cong)
    for tissue, want in expected.items():
        assert r["baf_final"][tissue] == pytest.approx(want, rel=GOLDEN_REL), tissue


def test_neutral_path_is_unmoved():
    """`pKa=None` must reach the code that produced the published neutral RMSEs.

    Note the two assertions carry DIFFERENT tolerances, on purpose: `baf_final` is
    an ODE output and so is machine-dependent in its last digits (see the comment
    on GOLDEN_REL), while `K_PW` is closed-form arithmetic — `W + L*a*Kow**b` —
    with no solver in it, so it is exactly reproducible and stays an exact `==`.
    Do not "tidy" the second one onto approx: it is the control showing the
    tolerance above is about the solver and not about the partition core.
    """
    r = api.simulate_neutral(2.45, name="carbamazepine", half_life=7.0)
    # PREREG 9.1: the root at its harvest water (0.6365, staged), lipid restated there
    assert r["baf_final"]["root"] == pytest.approx(2.8015531388751085, rel=GOLDEN_REL)
    assert r["K_PW"]["root"] == 3.4878512580341674
    # the compound the neutral path builds must carry NO speciation fields
    core = ND.neutral_compound(ND.NeutralCompound("x", 2.45))
    assert core.pKa is None and core.P_n == 0.0 and core.z is None
    assert (core.fn, core.fd) == (1.0, 1.0)


def test_root_uptake_neutral_term_is_identically_absent_by_default():
    """The added `j_n` term must not perturb an anion by so much as a float ulp."""
    env = LP.literature_environment()
    cmpd = LP.literature_compound("PFOA")
    assert cmpd.P_n == 0.0
    j = P4.root_uptake(1.0, 0.01, cmpd, env)
    # reconstruct without the neutral term and demand exact equality
    N = env.N_for(cmpd)
    import numpy as np
    g = P4._ghk_factor(N)
    j_ed = cmpd.kappa_d * g * (cmpd.fd * 1.0 - cmpd.fd * np.exp(N) * 0.01)
    j_carr = (cmpd.Vmax_in * 1.0 / (cmpd.Km_in + 1.0)
              - cmpd.Vmax_out * 0.01 / (cmpd.Km_out + 0.01))
    assert j == j_ed + j_carr


def test_N_for_falls_back_to_environment_z():
    env = LP.literature_environment()
    assert env.N_for(None) == env.N
    assert env.N_for(LP.literature_compound("PFOA")) == env.N   # z is None -> fallback


# --------------------------------------------------------------------------
# 3. the two modes meet continuously
# --------------------------------------------------------------------------
def test_weak_acid_far_above_its_pKa_reproduces_the_neutral_path():
    """THE safety property. A weak acid with pKa far above the pH is, physically,
    the same molecule as the neutral one -- so the two code paths must agree there.
    If they diverge, `pKa=` is a second model rather than an extension of this one,
    and `P_n = kappa_d` in `_weak_electrolyte_kw` is the choice that guarantees it."""
    kw = dict(name="x")
    neutral = api.simulate_neutral(2.45, **kw)
    almost = api.simulate_neutral(2.45, pKa=12.0, **kw)
    for tissue in ("root", "stem", "leaf"):
        assert almost["baf_final"][tissue] == pytest.approx(
            neutral["baf_final"][tissue], rel=1e-3), tissue
    # the brown rice is not predicted on either (phloem off, PREREG 9.1c); the husk is
    assert np.isnan(almost["baf_final"]["grain"]) and np.isnan(neutral["baf_final"]["grain"])
    assert almost["baf_final_dw"]["husk"] == pytest.approx(
        neutral["baf_final_dw"]["husk"], rel=1e-3)


def test_uptake_falls_monotonically_as_the_acid_gets_stronger():
    """As pKa drops the compound loses the fast neutral route AND gains anion
    exclusion, so root uptake must fall monotonically towards the PFAS-like limit."""
    kw = dict(name="x")
    roots = [api.simulate_neutral(2.45, pKa=p, **kw)["baf_final"]["root"]
             for p in (12.0, 8.0, 6.5, 4.5, -3.0)]
    assert all(a > b for a, b in zip(roots, roots[1:])), roots
    assert roots[0] / roots[-1] > 1000.0        # a large, physical collapse


# --------------------------------------------------------------------------
# 4. the physics the port exists for
# --------------------------------------------------------------------------
def test_a_weak_base_is_a_cation_and_is_not_excluded():
    """The reason valence had to move onto the COMPOUND. At the same pKa an acid
    (anion, z=-1) is excluded by the inside-negative membrane while a base (cation,
    z=+1) is attracted -- one global `Environment.z` cannot say both.

    The SIGN is the claim: the acid sits far below the neutral form, the base does not. "Not
    excluded" is checked against the neutral form of the same compound (the season ODE keeps
    every compound somewhat below its partition ceiling K_PW, so K_PW is not the yardstick)."""
    kw = dict(name="x", pKa=4.5)
    acid = api.simulate_neutral(2.45, is_acid=True, **kw)["baf_final"]["root"]
    base = api.simulate_neutral(2.45, is_acid=False, **kw)["baf_final"]["root"]
    neutral = api.simulate_neutral(2.45, name="x")["baf_final"]["root"]
    assert base > 3.0 * acid
    assert base == pytest.approx(neutral, rel=0.02)          # the cation is not excluded


def test_speciation_spans_the_whole_range():
    fn, fd = LP.speciation(4.0, 7.2)
    assert fd == pytest.approx(0.99937, abs=1e-4)
    assert LP.speciation(-3.3, 7.2)[1] == pytest.approx(1.0)      # PFAS: all anion
    assert LP.speciation(12.0, 7.2)[0] == pytest.approx(1.0, abs=1e-4)   # ~all neutral
    # a weak BASE dissociates the other way round
    assert LP.speciation(9.0, 7.2, is_acid=False)[1] > 0.9


def test_the_ion_trap_does_not_switch_off_thermodynamically():
    """The correction PR #54 exists for, pinned. Lambda is an EQUILIBRIUM ratio
    derived assuming the neutral species carries transport, so as pKa falls it tends
    to 10**(delta pH) = 6.31, NOT to 1. Multiplying L_Ph by it would hand a permanent
    anion a spurious ~6.3x phloem enrichment. What actually switches the trap off is
    kinetic -- f_n -> 0 -- which is `neutral_pathway_ratio`."""
    lam_acid = LP.ion_trap_factor(4.0, 7.2, 8.0)
    lam_pfas = LP.ion_trap_factor(-3.3, 7.2, 8.0)
    assert lam_acid == pytest.approx(6.31, abs=0.05)
    assert lam_pfas == pytest.approx(10.0 ** (8.0 - 7.2), abs=0.01)   # -> 6.31, not 1
    # the KINETIC statement is the one that separates them: ~7 orders of magnitude
    assert LP.neutral_pathway_ratio(4.0, 7.2) > 1.0
    assert LP.neutral_pathway_ratio(-3.3, 7.2) < 1e-6


def test_phloem_trap_is_off_for_pfas_by_construction():
    """Not a numerical coincidence: the PFAS path sets neither a compound pKa nor a
    leaf pH, and the factor returns `L_Ph` from the branch before any arithmetic."""
    res = api.simulate("PFOA")                          # smoke: PFAS still runs
    assert res["baf_final"]["root"] > 0

    comps = ND.rice_compartments()
    assert all(c.pH is None for c in comps)
    cmpd = LP.literature_compound("PFOA")
    assert cmpd.pKa is None
    env = LP.literature_environment()
    m = P4.RiceUptakeModel(env=env, cmpd=cmpd, comps=comps,
                           inputs=_dummy_inputs(), phi=PHI, T_C_Ph=T_C_PH)
    assert m.phloem_loading_factor() == cmpd.L_Ph        # exact, not approx


def test_phloem_trap_enriches_a_weak_acid_when_it_is_configured():
    comps = list(ND.rice_compartments())
    from dataclasses import replace
    comps[P4.LEAF] = replace(comps[P4.LEAF], pH=P4.LEAF_CYTOSOL_PH)
    cmpd = ND.neutral_compound(ND.NeutralCompound("x", 2.45, pKa=4.5))
    m = P4.RiceUptakeModel(env=LP.literature_environment(), cmpd=cmpd,
                           comps=comps, inputs=_dummy_inputs(), phi=PHI, T_C_Ph=T_C_PH)
    # L_Ph is 0 on the neutral base, so the trap route is the whole factor
    assert cmpd.L_Ph == 0.0
    assert m.phloem_loading_factor() > 1.0


def _dummy_inputs():
    import numpy as np
    t = np.linspace(0.0, 10.0, 11)
    M = np.full((11, 4), 0.01)
    return P4.PlantInputs(t=t, Cwo=np.ones_like(t), Qtp=np.full_like(t, 0.1), M=M)
