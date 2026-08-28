"""
Tests for the redistributed-shoot module (N stem segments + explicit leaf).

Locks in:
  * exact mass conservation (sole source = root membrane flux M_root*j_R for
    gamma=0), independent of the retention efficiency;
  * the transpiration split must close to 1 (structural invariant);
  * the structural CURE of the Tang over-translocation: vs the single-straw core,
    the redistributed model fills the stem (stalk TF up) and de-monopolizes the
    leaf (leaf burden fraction down) -- and more retention sends more to the stalk.
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

from pfas_rice_plant_module_nstem_leaf import (
    NStemLeafModel, PlantInputsNL, make_stem_leaf_compartments, split_from_stem_frac)
from pfas_rice_plant_module_4pool_surf import (
    Environment, Compound, binding_factors, root_uptake, _logistic)

N = 4


def _setup(retention=0.6, k_seq=0.0, k_rel=0.0):
    t = np.linspace(0.0, 150.0, 481)
    Qtp = 0.02 + 0.10 * np.exp(-((t - 95.0) ** 2) / (2 * 30.0 ** 2))
    stem = _logistic(t, 1e-3, 0.060, 0.11, 30.0)
    leaf = _logistic(t, 1e-3, 0.020, 0.11, 30.0)
    M = np.column_stack([_logistic(t, 1e-3, 0.030, 0.10, 20.0)]
                        + [stem / N] * N
                        + [leaf, _logistic(t, 1e-5, 0.035, 0.18, 90.0)])
    inputs = PlantInputsNL(t=t, Cwo=np.full_like(t, 1.0), Qtp=Qtp, M=M)
    comps = make_stem_leaf_compartments(
        N, dict(theta=0.83, f_prot=0.05, f_PL=0.005, f_cw=0.72),
        dict(theta=0.90, f_prot=0.07, f_PL=0.015, f_cw=0.50),
        dict(theta=0.78, f_prot=0.10, f_PL=0.010, f_cw=0.56),
        dict(theta=0.14, f_prot=0.09, f_PL=0.003, f_cw=0.035))
    cmpd = Compound(name="PFOA", K_prot=12.3, K_PL=1905.0, K_cw=6.5, kappa_d=2.0,
                    Vmax_in=20.0, Km_in=5.0, Vmax_out=8.0, Km_out=5.0, L_Ph=0.01, f_xy=0.04)
    tau, lam_leaf, lam_grain = split_from_stem_frac(N, 0.45, lam_grain=0.05)
    m = NStemLeafModel(env=Environment(), cmpd=cmpd, comps=comps, inputs=inputs,
                       tau=tau, lam_leaf=lam_leaf, lam_grain=lam_grain, retention=retention,
                       k_seq=k_seq, k_rel=k_rel)
    return t, m


def _mass_residual(m, t, ti=70.0):
    sol = m.solve(t)
    B = binding_factors(m.comps, m.cmpd)
    C = sol.sol(ti); dC = m.rhs(ti, C)
    M = m.inputs.M_(ti); dM = m.inputs.dM_(ti)
    jR = root_uptake(m.inputs.Cwo_(ti), (C[:m.n_comp] / B)[0], m.cmpd, m.env)
    n = m.n_comp
    dmass = float(np.sum(dM * C[:n] + M * dC[:n]))
    if m.two_pool:            # the sequestered pool shares the root mass
        dmass += float(dM[0] * C[m.SEQ] + M[0] * dC[m.SEQ])
    return dmass, float(M[0] * jR)


@pytest.mark.parametrize("retention", [0.0, 0.6, 1.0])
def test_mass_conservation(retention):
    t, m = _setup(retention)
    dmass, src = _mass_residual(m, t)
    assert dmass == pytest.approx(src, rel=1e-6, abs=1e-9)


def test_transpiration_split_must_close_to_one():
    t, m = _setup()
    tau = m.tau.copy()
    with pytest.raises(AssertionError):
        NStemLeafModel(env=m.env, cmpd=m.cmpd, comps=m.comps, inputs=m.inputs,
                       tau=tau, lam_leaf=0.9, lam_grain=0.9)   # sums > 1


def test_runs_finite_and_nonnegative():
    t, m = _setup()
    sol = m.solve(t)
    assert sol.success and np.all(np.isfinite(sol.y)) and np.all(sol.y >= -1e-9)


def test_more_retention_sends_more_to_stalk():
    t, m_lo = _setup(retention=0.2)
    _, m_hi = _setup(retention=1.0)
    C_lo = m_lo.solve(t).y[:, -1]; C_hi = m_hi.solve(t).y[:, -1]
    stalk_lo = m_lo.stem_aggregate(C_lo, t[-1]) / C_lo[0]
    stalk_hi = m_hi.stem_aggregate(C_hi, t[-1]) / C_hi[0]
    assert stalk_hi > stalk_lo            # retention accumulates the stalk


def test_cures_over_translocation_vs_single_straw():
    """The redistribution must fill the stem and de-monopolize the leaf."""
    import model_api as api
    base = api.simulate("PFOA", f_xy_source="recommended", measured_forcing=True, season=150.0, n_t=361)
    nl = api.simulate_nstem_leaf("PFOA", retention=0.6, stem_transp_frac=0.45, season=150.0)
    b_root = base["baf_final"]["root"]
    base_stalk_tf = base["baf_final"]["stem"] / b_root
    base_leaf_tf = base["baf_final"]["leaf"] / b_root
    nl_stalk_tf = nl["tf_final"]["stem"]; nl_leaf_tf = nl["tf_final"]["leaf"]
    # stem filled (was a near-empty pass-through), leaf no longer runs away
    assert nl_stalk_tf > 10 * base_stalk_tf
    assert nl_leaf_tf < 0.5 * base_leaf_tf
    # PFOA lands at O(1) across stalk/leaf/grain (Tang ~ 1.45/1.66/0.95)
    for k in ("stem", "leaf", "grain"):
        assert 0.3 < nl["tf_final"][k] < 3.5


# --------------------------------------------------------------------------
# STRUCTURAL MERGE: optional two-pool (mobile + sequestered) root
# --------------------------------------------------------------------------
def test_k_seq_zero_is_exactly_the_single_pool_model():
    """The merge must be inert when off: no extra state, identical trajectory."""
    t, m0 = _setup()
    assert m0.two_pool is False and m0.n_states == m0.n_comp and m0.SEQ is None
    y = m0.solve(t).y
    assert y.shape[0] == m0.n_comp
    np.testing.assert_allclose(m0.root_total(y), y[0])


@pytest.mark.parametrize("k_seq,k_rel", [(0.05, 0.0), (0.30, 0.0), (0.10, 0.02)])
def test_mass_conservation_two_pool(k_seq, k_rel):
    """Sequestration only MOVES solute inside the root -- the sole source is still
    the root membrane flux, so the total burden balance must still close."""
    t, m = _setup(k_seq=k_seq, k_rel=k_rel)
    assert m.two_pool and m.n_states == m.n_comp + 1
    dmass, src = _mass_residual(m, t)
    assert dmass == pytest.approx(src, rel=1e-6, abs=1e-9)


def test_k_seq_holds_root_burden_and_starves_the_shoot():
    """The point of the two-pool root: sequestration raises the root and lowers
    the shoot, and the sequestered pool holds the bulk of the root burden."""
    t, m0 = _setup()
    _, m1 = _setup(k_seq=0.30)
    C0 = m0.solve(t).y[:, -1]
    C1 = m1.solve(t).y[:, -1]
    assert m1.root_total(C1) > m0.root_total(C0)          # root holds more
    assert C1[m1.LEAF] < C0[m0.LEAF]                      # less delivered upward
    assert C1[m1.SEQ] > C1[m1.ROOT]                       # burden sits in the sink
    # ... and desorption (k_rel) gives some of it back to the shoot
    _, m2 = _setup(k_seq=0.30, k_rel=0.05)
    C2 = m2.solve(t).y[:, -1]
    assert C2[m2.LEAF] > C1[m1.LEAF]


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
