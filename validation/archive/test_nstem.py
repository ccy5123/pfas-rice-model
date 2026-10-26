"""
Tests for the N-segment ("multi-height stem") modules.

Locks in the structural invariants: exact mass conservation (sole source is the
root membrane flux M_root*j_R for gamma=0) for both the equilibrium and the
kinetic stem, and that the kinetic model recovers the equilibrium one as the
radial conductance k_rad -> infinity.
"""
import numpy as np
import pytest

from pfas_rice_plant_module_nstem import (
    NStemModel, NStemKineticModel, PlantInputsN, make_stem_compartments)
from pfas_rice_plant_module_4pool_surf import (
    Environment, Compound, binding_factors, root_uptake, _logistic)

N = 4


def _setup():
    t = np.linspace(0.0, 120.0, 481)
    Qtp = 0.05 + 0.35 * np.exp(-((t - 75.0) ** 2) / (2 * 25.0 ** 2))
    straw = _logistic(t, 1e-3, 0.09, 0.11, 27.0)
    M = np.column_stack([_logistic(t, 1e-3, 0.030, 0.10, 20.0)] + [straw / N] * N
                        + [_logistic(t, 1e-5, 0.025, 0.18, 80.0)])
    inputs = PlantInputsN(t=t, Cwo=np.full_like(t, 1.0), Qtp=Qtp, M=M)
    comps = make_stem_compartments(
        N, dict(theta=0.83, f_prot=0.05, f_PL=0.005, f_cw=0.72),
        dict(theta=0.90, f_prot=0.07, f_PL=0.015, f_cw=0.50),
        dict(theta=0.14, f_prot=0.09, f_PL=0.003, f_cw=0.035))
    cmpd = Compound(name="PFOA", K_prot=12.3, K_PL=1905.0, K_cw=6.6, kappa_d=2.0,
                    Vmax_in=20.0, Km_in=5.0, Vmax_out=8.0, Km_out=5.0, L_Ph=0.01, f_xy=0.04)
    tau = np.array([0.30, 0.28, 0.24, 0.18]); tau = tau / tau.sum() * 0.85
    return t, inputs, comps, cmpd, tau


def _mass_source_residual(model, t, ti=70.0):
    sol = model.solve(t)
    B = binding_factors(model.comps, model.cmpd)
    C = sol.sol(ti); dC = model.rhs(ti, C)
    M = model.inputs.M_(ti); dM = model.inputs.dM_(ti)
    jR = root_uptake(model.inputs.Cwo_(ti), (C / B)[0], model.cmpd, model.env)
    dmass = float(np.sum(dM * C + M * dC))          # d/dt sum_k M_k C_k
    return dmass, float(M[0] * jR)                  # vs sole source M_root*j_R


def test_equilibrium_nstem_conserves_mass():
    t, inputs, comps, cmpd, tau = _setup()
    m = NStemModel(env=Environment(), cmpd=cmpd, comps=comps, inputs=inputs, tau=tau)
    dmass, src = _mass_source_residual(m, t)
    assert dmass == pytest.approx(src, rel=1e-6, abs=1e-9)


def test_kinetic_nstem_conserves_mass():
    t, inputs, comps, cmpd, tau = _setup()
    m = NStemKineticModel(env=Environment(), cmpd=cmpd, comps=comps,
                          inputs=inputs, tau=tau, k_rad=0.5)
    dmass, src = _mass_source_residual(m, t)
    assert dmass == pytest.approx(src, rel=1e-6, abs=1e-9)


def test_kinetic_recovers_equilibrium_at_large_k_rad():
    t, inputs, comps, cmpd, tau = _setup()
    eq = NStemModel(env=Environment(), cmpd=cmpd, comps=comps, inputs=inputs, tau=tau)
    kin = NStemKineticModel(env=Environment(), cmpd=cmpd, comps=comps,
                            inputs=inputs, tau=tau, k_rad=1e6)
    Ceq = eq.solve(t).y[:, -1]
    Ckin = kin.solve(t).y[:, -1]
    assert np.allclose(Ceq, Ckin, rtol=2e-2, atol=1e-6)


def test_nstem_runs_and_is_finite():
    t, inputs, comps, cmpd, tau = _setup()
    for m in (NStemModel(env=Environment(), cmpd=cmpd, comps=comps, inputs=inputs, tau=tau),
              NStemKineticModel(env=Environment(), cmpd=cmpd, comps=comps, inputs=inputs,
                                tau=tau, k_rad=0.5)):
        sol = m.solve(t)
        assert sol.success and np.all(np.isfinite(sol.y)) and np.all(sol.y >= -1e-9)


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))
