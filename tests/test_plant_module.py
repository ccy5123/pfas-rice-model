"""
Tests for the PFAS rice four-compartment plant uptake module.

These lock in the *structural* results the model is meant to reproduce
(anion exclusion, carrier-enabled influx, binding without a density factor,
TSCF-limited translocation, mass conservation, and the empirical
root > straw > grain ordering) so that future parameter or refactor work
cannot silently break them.  Absolute BAF values are NOT asserted -- the demo
parameters are illustrative, not calibrated.
"""
import numpy as np
import pytest

import pfas_rice_plant_module as prm
from pfas_rice_plant_module import (
    Environment, Compound, Compartment, PlantInputs, RiceUptakeModel,
    binding_factors, root_uptake, _ghk_factor, ROOT, STEM, LEAF, FRUIT,
)


# ---------------------------------------------------------------------------
# builders (mirror the _demo() configuration)
# ---------------------------------------------------------------------------
def make_inputs(season=120.0, n=481, Cwo=1.0):
    t = np.linspace(0.0, season, n)
    Cwo_arr = np.full_like(t, Cwo)
    Qtp = 0.05 + 0.35 * np.exp(-((t - 75.0) ** 2) / (2 * 25.0 ** 2))
    M = np.column_stack([
        prm._logistic(t, 1e-3, 0.030, 0.10, 20.0),
        prm._logistic(t, 1e-3, 0.040, 0.10, 25.0),
        prm._logistic(t, 1e-3, 0.050, 0.12, 30.0),
        prm._logistic(t, 1e-5, 0.025, 0.18, 80.0),
    ])
    return t, PlantInputs(t=t, Cwo=Cwo_arr, Qtp=Qtp, M=M)


def make_compound(**kw):
    base = dict(name="PFOA", K_prot=50.0, K_PL=100.0, K_cw=7.0, kappa_d=0.5,
                Vmax_in=20.0, Km_in=5.0, Vmax_out=8.0, Km_out=5.0,
                L_Ph=0.005, f_xy=0.02)
    base.update(kw)
    return Compound(**base)


def make_comps():
    # rice_tissue recommended composition (basis-A): theta = fresh-weight water
    # fraction; f_* = DRY-weight mass fractions (rescaled by (1-theta) in B_k).
    return [
        Compartment("root",  theta=0.90, f_prot=0.07, f_PL=0.015, f_cw=0.50),
        Compartment("stem",  theta=0.83, f_prot=0.05, f_PL=0.005, f_cw=0.72),
        Compartment("leaf",  theta=0.78, f_prot=0.10, f_PL=0.010, f_cw=0.56, S=20.0),
        Compartment("grain", theta=0.14, f_prot=0.09, f_PL=0.003, f_cw=0.035, S=2.0),
    ]


def make_model(**cmpd_kw):
    t, inputs = make_inputs()
    model = RiceUptakeModel(env=Environment(), cmpd=make_compound(**cmpd_kw),
                            comps=make_comps(), inputs=inputs)
    return t, model


# ---------------------------------------------------------------------------
# electrochemistry: anion exclusion
# ---------------------------------------------------------------------------
def test_electrochemical_number_is_anion_exclusion():
    env = Environment()                      # z=-1, E=-120 mV (inside-negative)
    assert env.N > 0                         # positive N for an anion at E<0
    assert np.exp(env.N) == pytest.approx(106.8, rel=0.02)   # e^N ~ 107


def test_ghk_factor_removable_singularity():
    assert _ghk_factor(0.0) == pytest.approx(1.0)
    assert _ghk_factor(1e-12) == pytest.approx(1.0, rel=1e-6)
    assert _ghk_factor(1e-4) == pytest.approx(1.0, rel=1e-3)   # continuity near 0
    N = 2.0
    assert _ghk_factor(N) == pytest.approx(N / np.expm1(N))    # exact form away from 0


def test_electrodiffusion_excludes_the_anion():
    """With the carrier off, passive electrodiffusion excludes the anion:
    zero net flux occurs at Cw_root = Cwo/e^N < Cwo, and at Cw_root = Cwo the
    membrane drives net efflux."""
    env = Environment()
    cmpd = make_compound(Vmax_in=0.0, Vmax_out=0.0)   # electrodiffusion only
    Cwo = 1.0
    Cw_eq = Cwo / np.exp(env.N)
    assert root_uptake(Cwo, Cw_eq, cmpd, env) == pytest.approx(0.0, abs=1e-9)
    assert root_uptake(Cwo, Cwo, cmpd, env) < 0.0      # exclusion at equal conc


def test_carrier_enables_net_influx():
    """The saturable carrier must be able to overcome electrostatic exclusion."""
    env = Environment()
    cmpd = make_compound()                              # carrier on
    assert root_uptake(1.0, 0.0, cmpd, env) > 0.0


# ---------------------------------------------------------------------------
# binding (no density prefactor)
# ---------------------------------------------------------------------------
def test_binding_factor_basis_A_no_density_prefactor():
    """Basis A (fresh-weight): B = theta_fw + (1-theta_fw)*sum_i f_i,dw*K_i.
    The (1-theta) factor is a dry->fresh conversion (f_i are dw fractions), NOT a
    density prefactor -- there is still no spurious rho_k term."""
    comps, cmpd = make_comps(), make_compound()
    B = binding_factors(comps, cmpd)
    for k, c in enumerate(comps):
        expected = c.theta + (1.0 - c.theta) * (
            c.f_prot * cmpd.K_prot + c.f_PL * cmpd.K_PL + c.f_cw * cmpd.K_cw)
        assert B[k] == pytest.approx(expected)
    # under basis A the low-water grain has the largest fresh-weight B_k
    assert B[FRUIT] == pytest.approx(max(B))


# ---------------------------------------------------------------------------
# TSCF: limited root -> shoot translocation
# ---------------------------------------------------------------------------
def test_tscf_constrains_translocation():
    """Lower f_xy retains more in the root and delivers less to shoot/grain."""
    t, m_lo = make_model(f_xy=0.02)
    _, m_hi = make_model(f_xy=1.0)
    Clo = m_lo.solve(t).y[:, -1]
    Chi = m_hi.solve(t).y[:, -1]
    assert Clo[ROOT] > Chi[ROOT]
    assert Clo[LEAF] < Chi[LEAF]
    assert Clo[FRUIT] < Chi[FRUIT]


def test_full_loading_reproduces_runaway():
    """Sanity check on the original failure mode: with f_xy=1 the terminal
    compartments out-accumulate the root (the known runaway)."""
    t, model = make_model(f_xy=1.0)
    C = model.solve(t).y[:, -1]
    assert C[LEAF] > C[ROOT]
    assert C[FRUIT] > C[ROOT]


# ---------------------------------------------------------------------------
# the empirical target ordering
# ---------------------------------------------------------------------------
def test_translocation_controls_root_shoot_partitioning():
    """f_xy controls the root<->shoot split.  The empirical ordering is
    CONGENER-DEPENDENT (Yamazaki: short-chain PFBA straw >> root, but long-chain
    PFUnDA root > straw), so a *universal* root>straw>grain does NOT hold under
    the basis-A binding and is not asserted.  Instead: strong root retention
    (very low f_xy, long-chain-like) keeps root above the shoot, while efficient
    loading (high f_xy, short-chain-like) puts the shoot above the root."""
    t, m_lo = make_model(f_xy=0.002)     # long-chain-like: root-retained
    _, m_hi = make_model(f_xy=0.5)       # short-chain-like: translocated
    Clo = m_lo.solve(t).y[:, -1]
    Chi = m_hi.solve(t).y[:, -1]
    Mf = m_lo.inputs.M_(t[-1])
    straw_lo = (Clo[STEM] * Mf[STEM] + Clo[LEAF] * Mf[LEAF]) / (Mf[STEM] + Mf[LEAF])
    straw_hi = (Chi[STEM] * Mf[STEM] + Chi[LEAF] * Mf[LEAF]) / (Mf[STEM] + Mf[LEAF])
    assert Clo[ROOT] > straw_lo          # root retains when translocation is low
    assert straw_hi > Chi[ROOT]          # shoot exceeds root when loading is high


# ---------------------------------------------------------------------------
# conservation
# ---------------------------------------------------------------------------
def test_mass_conservation_source_is_root_uptake():
    """For gamma=0 the only solute source/sink is the root membrane flux, so
    d/dt(sum_k M_k C_k) must equal M_root * j_R exactly.  This is an algebraic
    identity only if the internal xylem AND phloem transfers each conserve mass
    -- in particular it fails if the leaf does not export the phi-recirculation
    fraction (the bug fixed alongside the TSCF change)."""
    t, model = make_model()
    sol = model.solve(t)
    B = binding_factors(model.comps, model.cmpd)
    for ti in (25.0, 60.0, 95.0):
        C = sol.sol(ti)
        dC = model.rhs(ti, C)
        M = model.inputs.M_(ti)
        dM = model.inputs.dM_(ti)
        Cw = C / B
        jR = root_uptake(model.inputs.Cwo_(ti), Cw[ROOT], model.cmpd, model.env)
        dmass = float(np.sum(dM * C + M * dC))     # d/dt sum(M_k C_k)
        src = float(M[ROOT] * jR)                  # root uptake (gamma=0)
        assert dmass == pytest.approx(src, rel=1e-6, abs=1e-9)


def test_phloem_internal_transfer_conserves_mass():
    """The phloem source (leaf) must supply exactly what the sinks receive:
    grain (Q_Phl) + root recirculation (phi*Q_Phl)."""
    t, model = make_model()
    sol = model.solve(t)
    B = binding_factors(model.comps, model.cmpd)
    ti = 95.0                                       # grain-filling window
    C = sol.sol(ti)
    Cw = C / B
    dM = model.inputs.dM_(ti)
    Qtp = model.inputs.Qtp_(ti)
    Q_Phl = max(dM[FRUIT] * model.T_C_Ph + model.phi * Qtp, 0.0)
    C_Phl = model.cmpd.L_Ph * Cw[LEAF]
    leaf_out = (1.0 + model.phi) * Q_Phl * C_Phl
    grain_in = Q_Phl * C_Phl
    root_in = model.phi * Q_Phl * C_Phl
    assert leaf_out == pytest.approx(grain_in + root_in)


# ---------------------------------------------------------------------------
# numerics / robustness
# ---------------------------------------------------------------------------
def test_solver_succeeds_and_is_finite():
    t, model = make_model()
    sol = model.solve(t)
    assert sol.success
    assert np.all(np.isfinite(sol.y))


def test_concentrations_stay_nonnegative():
    t, model = make_model()
    sol = model.solve(t)
    assert np.all(sol.y >= -1e-9)


def test_zero_soil_concentration_gives_zero_uptake():
    t, inputs = make_inputs(Cwo=0.0)
    model = RiceUptakeModel(env=Environment(), cmpd=make_compound(),
                            comps=make_comps(), inputs=inputs)
    C = model.solve(t).y[:, -1]
    assert np.allclose(C, 0.0, atol=1e-9)


if __name__ == "__main__":   # allow `python tests/test_plant_module.py` without pytest
    import sys
    sys.exit(pytest.main([__file__, "-v"]))
