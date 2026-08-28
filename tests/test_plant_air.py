"""
Tests for plant-air exchange -- src/plant_air.py (volatilisation + gaseous uptake).

Two things have to hold at once, and they pull in opposite directions:

  * the term must be REAL -- a volatile compound now loses leaf mass, an ambient
    air concentration now loads the shoot from clean soil, and both follow the
    equation set in docs/dpu_model_summary_corrected.tex sec:permeability;
  * the term must be ABSENT for PFAS -- `K_AW = 0` has to give not "small" but
    exactly zero, and the core's default (`air=None`) must not evaluate it at all,
    so every existing number and `reproduce_demo`'s RMSE 0.029 are untouched.

The unit bridge is the other thing worth pinning: the published correlations are
SI (m/s, g/mol) while the model runs in day/L/kg/ug, so the tests below check the
conversions against the raw formulas rather than against the implementation.
"""
import os
import sys

import numpy as np
import pytest

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "src"))
sys.path.insert(0, os.path.join(_ROOT, "validation"))

import plant_air as PA  # noqa: E402
import neutral_dpu as ND  # noqa: E402
from pfas_rice_plant_module_4pool_surf import (  # noqa: E402
    ROOT, STEM, LEAF, FRUIT, RiceUptakeModel, binding_factors,
)


@pytest.fixture(scope="module")
def drv():
    import neutral_dpu_validation as V
    return V.drivers()


def _setup(log_kow=2.0, MW=150.0, **kw):
    """An AirExchange plus the compartments / K_PW it acts on."""
    comps = ND.rice_compartments()
    core = ND.neutral_compound(ND.NeutralCompound("probe", log_kow))
    K_PW = binding_factors(comps, core)
    return PA.AirExchange(MW=MW, log_kow=log_kow, **kw), comps, K_PW


# --------------------------------------------------------------------------
# the published equations, and the SI -> day/L/kg bridge
# --------------------------------------------------------------------------
def test_cuticle_permeability_is_the_published_correlation():
    """eq:Pc, `10^(0.704 logKow - 11.2)` m/s, expressed in m/day."""
    for lk in (0.0, 2.0, 4.5):
        assert PA.p_cuticle(lk) == pytest.approx(
            10.0 ** (0.704 * lk - 11.2) * 86400.0, rel=1e-12)


def test_air_boundary_layer_is_proportional_to_kaw():
    """eq:Pair. Linearity in K_AW is what makes the whole pathway vanish for a
    perfluorinated anion -- it is not an approximation, it is the equation."""
    assert PA.p_air_boundary(1e-3, 150.0) == pytest.approx(
        np.sqrt(300.0) * 1e-3 / (200.0 * np.sqrt(150.0)) * 86400.0, rel=1e-12)
    assert PA.p_air_boundary(0.0, 150.0) == 0.0
    assert PA.p_air_boundary(2e-3, 150.0) == pytest.approx(
        2 * PA.p_air_boundary(1e-3, 150.0), rel=1e-12)


def test_aqueous_permeability_scales_oxygen_by_molar_mass():
    """eq:Paqua: oxygen's diffusivity over the path length, scaled by (MW/32)^-0.5."""
    assert PA.p_aqueous(32.0, z_path=1e-3) == pytest.approx(
        PA.D_O2_WATER / 1e-3 * 86400.0, rel=1e-12)
    assert PA.p_aqueous(128.0) == pytest.approx(0.5 * PA.p_aqueous(32.0), rel=1e-12)
    with pytest.raises(ValueError):
        PA.p_aqueous(150.0, z_path=0.0)


def test_saturation_vapour_concentration():
    """eq:Csat: ~0.023 kg/m^3 at 25 C, the textbook value."""
    assert PA.c_h2o_sat(298.15) == pytest.approx(0.0230, abs=0.0005)
    assert PA.c_h2o_sat(288.15) < PA.c_h2o_sat(298.15)   # colder air holds less


def test_series_resistances_are_limited_by_the_slowest():
    """eq:Pctot. In series the total is below every member and dominated by the
    smallest -- and a zero member (K_AW=0) makes the total exactly zero rather
    than a division by zero."""
    tot = PA.p_cuticle_total(log_kow=2.0, K_AW=1e-3, MW=150.0)
    members = (PA.p_cuticle(2.0), PA.p_air_boundary(1e-3, 150.0), PA.p_aqueous(150.0))
    assert tot < min(members)
    assert tot == pytest.approx(1.0 / sum(1.0 / m for m in members), rel=1e-12)
    assert PA.p_cuticle_total(log_kow=2.0, K_AW=0.0, MW=150.0) == 0.0


def test_stomata_and_cuticle_act_in_parallel():
    """eq:Pp: parallel conductances ADD (unlike the series inside eq:Pctot)."""
    air, _, _ = _setup(K_AW=1e-2)
    q = 1e-4
    assert air.p_plant(q, stomata=True) == pytest.approx(
        air.p_cuticular() + PA.p_stomata(1e-2, 150.0, q, air.rh, air.T), rel=1e-12)
    assert air.p_plant(q, stomata=False) == pytest.approx(air.p_cuticular(), rel=1e-12)


def test_humidity_pole_is_guarded():
    """eq:Ps carries 1/(1-phi), singular at saturation. The tex notes the product
    stays finite because Q_TP -> 0 there; the factor itself is capped so the
    arithmetic cannot overflow."""
    assert np.isfinite(PA.p_stomata(1e-2, 150.0, 1e-4, rh=1.0))
    assert PA.p_stomata(1e-2, 150.0, 1e-4, rh=1.0) == pytest.approx(
        PA.p_stomata(1e-2, 150.0, 1e-4, rh=PA.RH_MAX), rel=1e-12)
    assert PA.p_stomata(1e-2, 150.0, 1e-4, rh=0.9) > PA.p_stomata(1e-2, 150.0, 1e-4, rh=0.5)
    with pytest.raises(ValueError):
        PA.p_stomata(1e-2, 150.0, 1e-4, rh=-0.1)


def test_molar_mass_is_required():
    """Every air-side equation needs MW, so a compound without one must be refused
    rather than run with a silently wrong permeability."""
    with pytest.raises(ValueError, match="molar mass"):
        PA.AirExchange(K_AW=1e-3, MW=float("nan"), log_kow=2.0)
    with pytest.raises(ValueError):
        PA.AirExchange(K_AW=-1.0, MW=150.0, log_kow=2.0)


# --------------------------------------------------------------------------
# the derivation's modelling assumptions
# --------------------------------------------------------------------------
def test_no_exchange_from_roots():
    """'No volatilisation is assumed from the roots' -- the root is below ground."""
    air, comps, K_PW = _setup(K_AW=0.5, C_air=1e3)
    f = air.flux(C=np.ones(4), M=np.full(4, 0.05), K_PW=K_PW, comps=comps,
                 Qtp=0.1, xyl_share=(0.9, 0.1))
    assert f[ROOT] == 0.0


def test_stem_is_cuticle_only_and_inert_without_a_surface_area():
    """'Only the cuticular pathway is considered for the stem.' And the shipped
    rice compartments give the stem S = 0, so the term stays inert until a real
    stem surface area is supplied -- a documented scope limit, pinned here so it
    cannot silently become a hidden number."""
    air, comps, K_PW = _setup(K_AW=0.5)
    assert comps[STEM].S == 0.0
    rows = air.summary(comps, K_PW, M=np.full(4, 0.05), Qtp=0.1)
    assert rows["stem"]["P_stomatal"] == 0.0        # no stomata on the stem
    assert rows["leaf"]["P_stomatal"] > 0.0
    assert rows["stem"]["rate"] == 0.0              # S = 0 -> inert
    with_area = PA.AirExchange(K_AW=0.5, MW=150.0, log_kow=2.0,
                               S=dict(stem=5.0, leaf=20.0, grain=2.0))
    assert with_area.summary(comps, K_PW, M=np.full(4, 0.05), Qtp=0.1)["stem"]["rate"] > 0.0


# --------------------------------------------------------------------------
# the K_AW = 0 constraint: PFAS must be untouched
# --------------------------------------------------------------------------
def test_core_defaults_to_no_air_exchange():
    """The hook must be opt-in: with `air=None` the term is not merely zero, it is
    never evaluated -- which is what keeps every PFAS result bit-identical."""
    assert "air" in RiceUptakeModel.__dataclass_fields__
    assert RiceUptakeModel.__dataclass_fields__["air"].default is None


def test_flux_is_identically_zero_at_kaw_zero():
    air, comps, K_PW = _setup(K_AW=0.0, C_air=1e6)
    f = air.flux(C=np.full(4, 100.0), M=np.full(4, 0.05), K_PW=K_PW, comps=comps,
                 Qtp=0.1, xyl_share=(0.9, 0.1))
    assert np.array_equal(f, np.zeros(4))           # exactly zero, not approximately
    assert air.p_plant(1e-4, stomata=True) == 0.0   # both pathways vanish
    assert np.isinf(air.k_pa(K_PW[LEAF]))           # partition undefined, flux still fine


def test_enabling_air_cannot_perturb_a_nonvolatile_run(drv):
    """The PFAS-safety property end to end: switching air exchange ON for a
    compound with K_AW = 0 must reproduce the air-off trajectory bit for bit."""
    c = ND.NeutralCompound("nonvolatile", 2.0, MW=200.0, K_AW=0.0)
    off = ND.simulate_neutral(c, drv)
    on = ND.simulate_neutral(c, drv, air=True)
    assert on["air"] is True and off["air"] is False
    for k in ("root", "stem", "leaf", "grain"):
        assert np.array_equal(off["conc"][k], on["conc"][k])


# --------------------------------------------------------------------------
# the term is real
# --------------------------------------------------------------------------
def test_volatilisation_bounds_the_leaf(drv):
    """The point of the exercise. Without metabolism the leaf is an unbounded
    terminal accumulator (neutral_dpu_validation section 3) -- volatilisation is
    the second sink, so the leaf must fall monotonically as K_AW rises, while the
    root (no air exchange) does not move."""
    def run(kaw):
        return ND.simulate_neutral(
            ND.NeutralCompound("p", 2.42, MW=131.4, K_AW=kaw), drv, air=True)

    leaf = [run(k)["baf_final"]["leaf"] for k in (0.0, 1e-5, 1e-3, 1e-1)]
    assert leaf[0] > leaf[1] > leaf[2] > leaf[3]
    assert leaf[3] < 1e-3 * leaf[0]             # a truly volatile compound is stripped
    roots = [run(k)["baf_final"]["root"] for k in (0.0, 1e-1)]
    assert roots[0] == pytest.approx(roots[1], rel=1e-9)


def test_gaseous_uptake_loads_the_shoot_from_clean_soil(drv):
    """eq:Qgas is a SOURCE: with no soil exposure at all, ambient air alone must
    load the leaf -- and only the above-ground organs."""
    c = ND.NeutralCompound("p", 2.0, MW=150.0, K_AW=1e-3)
    clean = dict(drv, Cwo=np.zeros(len(drv["t"])))
    r = ND.simulate_neutral(c, clean, air=True, air_kw=dict(C_air=10.0))
    assert r["conc"]["leaf"][-1] > 0.0
    assert r["conc"]["root"][-1] == pytest.approx(0.0, abs=1e-12)
    # and the particle-bound fraction is excluded from the gaseous uptake
    half = ND.simulate_neutral(c, clean, air=True,
                               air_kw=dict(C_air=10.0, f_particle=0.5))
    assert half["conc"]["leaf"][-1] == pytest.approx(0.5 * r["conc"]["leaf"][-1], rel=1e-6)


def test_equilibrium_with_the_air_gives_zero_net_flux():
    """Volatilisation and gaseous uptake must balance exactly at Henry's law on the
    tissue's FREE aqueous concentration, `C_air = 1000 * K_AW * C / K_PW`.

    This is the test that pins the unit bridge: the 1000 is the m^3 -> L factor
    between a per-cubic-metre air concentration and the model's per-litre
    aqueous one, and getting it wrong is a thousandfold error that nothing else
    here would catch.
    """
    comps = ND.rice_compartments()
    core = ND.neutral_compound(ND.NeutralCompound("probe", 2.0))
    K_PW = binding_factors(comps, core)
    C = np.array([1.0, 2.0, 3.0, 4.0])
    K_AW = 1e-3
    air = PA.AirExchange(K_AW=K_AW, MW=150.0, log_kow=2.0,
                         C_air=PA.M3_TO_L * K_AW * C[LEAF] / K_PW[LEAF])
    f = air.flux(C=C, M=np.array([0.03, 0.04, 0.05, 0.02]), K_PW=K_PW, comps=comps,
                 Qtp=0.1, xyl_share=(0.9, 0.1))
    assert f[LEAF] == pytest.approx(0.0, abs=1e-12)
    # the grain sits at a different concentration, so it is NOT at equilibrium
    assert f[FRUIT] < 0.0


def test_volatilisation_half_life_is_reported_and_ordered():
    """`summary` turns the term into the same currency as the metabolic half-life
    the neutral path already scans, so the two sinks can be compared."""
    comps = ND.rice_compartments()
    core = ND.neutral_compound(ND.NeutralCompound("probe", 2.0))
    K_PW = binding_factors(comps, core)
    M, Q = np.array([0.03, 0.04, 0.05, 0.02]), 0.1
    hl = [PA.AirExchange(K_AW=k, MW=150.0, log_kow=2.0)
          .summary(comps, K_PW, M=M, Qtp=Q)["leaf"]["half_life"]
          for k in (0.0, 1e-4, 1e-2)]
    assert np.isinf(hl[0])                       # PFAS: no air pathway at all
    assert hl[1] > hl[2] > 0.0                   # more volatile -> faster loss


def test_warning_is_replaced_by_a_run(drv):
    """The pre-A1 behaviour was to REFUSE a volatile compound with a warning. With
    air exchange available the warning belongs only to the air-off run, and the
    remedy it names must actually work."""
    c = ND.NeutralCompound("volatile", 2.0, MW=150.0, K_AW=0.1)
    off = ND.simulate_neutral(c, drv)
    on = ND.simulate_neutral(c, drv, air=True)
    assert off["warning"] and "air=True" in off["warning"]
    assert on["warning"] is None
    assert on["air_summary"]["leaf"]["half_life"] > 0.0
    assert off["air_summary"] is None
    # the warning's premise: without air terms the shoot really is an upper bound
    assert off["baf_final"]["leaf"] > on["baf_final"]["leaf"]


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
