"""Make the plant module importable from tests without installing the package.

Also fixes the tolerance for one recurring kind of assertion -- see below.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))


# ---------------------------------------------------------------------------
# How tight may an assertion about ODE output be?
# ---------------------------------------------------------------------------
# Three kinds of comparison appear in this suite and they have DIFFERENT floors.
# Getting this wrong produced a run of CI failures that all passed locally, so
# the rule is written down once here instead of guessed per test:
#
#   1. TWO READS OF ONE SOLVE -- e.g. `baf["root"][-1]` vs `baf_final["root"]`,
#      or two `simulate()` calls whose arguments build the SAME Compound by
#      different routes (`lipid_loading=False` vs `f_xy_source="recommended"`).
#      Same RHS -> same adaptive steps -> the same floats. Exact `==`, or an
#      essentially decorative rel=1e-9, is correct and should stay.
#
#   2. A TERM THAT VANISHES STRUCTURALLY -- e.g. air exchange at K_AW = 0, or
#      `g_apo = 0`. The added expression is identically zero, `x + 0.0 == x`, so
#      the RHS is unchanged bit-for-bit. These SHOULD be exact (`np.array_equal`)
#      and that is the strongest guard in the suite. Do not loosen them.
#
#   3. TWO DISTINCT SOLVES COMPARED BECAUSE THE PHYSICS SAYS THEY AGREE -- e.g.
#      "air exchange must not move the root" at K_AW = 0.1 vs air off. Here the
#      equations agree exactly but the NUMERICS need not: the leaf state differs,
#      so `solve_ivp` selects a different step sequence for the coupled system,
#      and the root is re-integrated along it. The result can only agree to the
#      solver's own accuracy.
#
# `SOLVE_INVARIANCE_REL` is for case 3, and it is derived rather than observed --
# picking it from a couple of observed differences is what went wrong twice.
# `solve_ivp(method="BDF", rtol=1e-6, atol=1e-9)` defines the trajectory to ~1e-6
# relative per step, so no claim about two different solves can be tighter than
# that; one order of headroom for accumulation over a season gives 1e-5.
#
# It stays far more sensitive than it needs to be: a leak of the kind these tests
# guard against is enormous, not marginal. Air exchange at K_AW = 0.1 takes the
# leaf BAF from 177 to 0.0025 -- if any of that reached the root it would show up
# thousands of times above this bound.
SOLVE_INVARIANCE_REL = 1e-5
