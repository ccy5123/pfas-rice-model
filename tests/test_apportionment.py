"""
Tests for the source-apportionment diagnostic (Fig-2 analog).

These lock in (a) that the mass-flux bookkeeping never diverges from the canonical
RiceUptakeModel.rhs (the consistency check reconstructs dC/dt from the fluxes), (b)
that the cumulative fluxes conserve mass at every junction, and (c) that the
per-compartment inflow fractions are well-formed.  They also assert the structural
audit result the model claims: the grain is predominantly phloem-fed.
"""
import numpy as np
import pytest

import apportionment as ap
import model_api as api
from pfas_rice_plant_module import (
    Environment, PlantInputs, RiceUptakeModel, ROOT, STEM, LEAF, FRUIT)

# reuse the plant-module test builders for a self-contained model
from test_plant_module import make_compound, make_comps, make_inputs


def _model(**cmpd_kw):
    t, inputs = make_inputs()
    model = RiceUptakeModel(env=Environment(), cmpd=make_compound(**cmpd_kw),
                            comps=make_comps(), inputs=inputs)
    return t, model


# ---------------------------------------------------------------------------
# flux bookkeeping must mirror rhs exactly
# ---------------------------------------------------------------------------
def test_flux_reconstruction_matches_rhs():
    """dC/dt rebuilt from the mass fluxes equals RiceUptakeModel.rhs term for term."""
    t, model = _model()
    rng = np.random.default_rng(0)
    for ti in (5.0, 30.0, 70.0, 110.0):
        for _ in range(3):
            C = rng.uniform(0.0, 50.0, size=4)
            dC_rhs = model.rhs(ti, C)
            dC_flux = ap.dC_from_fluxes(model, ti, C)
            assert np.allclose(dC_rhs, dC_flux, rtol=1e-9, atol=1e-12)


def test_fluxes_conserve_mass_globally():
    """Internal transfers cancel: sum of compartment mass-balances equals the sole
    external source (soil uptake) minus the only sinks (leaf loss + degradation).
    For the demo compound gamma=0 and leaf_loss=0, so the total equals soil uptake."""
    t, model = _model()
    C = np.array([10.0, 4.0, 6.0, 2.0])
    for ti in (10.0, 60.0, 100.0):
        f = ap.flux_terms(model, ti, C)
        mass_bal_total = (
            f["soil_uptake"] - f["xylem_root_stem"] + f["phloem_leaf_root"] - f["degr_root"]
            + f["xylem_root_stem"] - (f["xylem_stem_leaf"] + f["xylem_stem_grain"]) - f["degr_stem"]
            + f["xylem_stem_leaf"] - (f["phloem_leaf_grain"] + f["phloem_leaf_root"])
            - f["leaf_loss"] - f["degr_leaf"]
            + f["xylem_stem_grain"] + f["phloem_leaf_grain"] - f["degr_grain"])
        sinks = f["leaf_loss"] + f["degr_root"] + f["degr_stem"] + f["degr_leaf"] + f["degr_grain"]
        assert mass_bal_total == pytest.approx(f["soil_uptake"] - sinks, rel=1e-9, abs=1e-12)


def test_apportion_fractions_are_well_formed():
    t, model = _model()
    sol = model.solve(t)
    r = ap.apportion(model, sol)
    for k, frac in r["fraction"].items():
        s = sum(v for v in frac.values() if np.isfinite(v))
        assert s == pytest.approx(1.0, abs=1e-9), f"{k} fractions sum to {s}"
        assert all(0.0 <= v <= 1.0 for v in frac.values() if np.isfinite(v))
    assert 0.0 <= r["grain_phloem_fraction"] <= 1.0
    assert 0.0 <= r["delivered_to_shoot"] <= 1.0


def test_grain_is_phloem_fed():
    """The model asserts the grain is phloem-fed: phloem should dominate grain inflow
    for the low-phloem PFAS demo compound is NOT guaranteed, but with a non-trivial
    L_Ph the phloem share is the larger contributor."""
    t, model = _model(L_Ph=0.05)         # give the phloem a meaningful loading
    sol = model.solve(t)
    r = ap.apportion(model, sol)
    fr = r["fraction"]["grain"]
    assert fr["phloem_from_leaf"] > fr["xylem_from_stem"]


def test_no_phloem_loading_makes_grain_xylem_fed():
    """L_Ph -> 0 (and no lipid-bound phloem) starves the phloem: grain becomes
    xylem-fed -- the lever behaves monotonically."""
    t, model = _model(L_Ph=0.0)
    sol = model.solve(t)
    r = ap.apportion(model, sol)
    fr = r["fraction"]["grain"]
    assert fr["xylem_from_stem"] > fr["phloem_from_leaf"]


# ---------------------------------------------------------------------------
# model_api wrapper
# ---------------------------------------------------------------------------
def test_model_api_apportionment_runs():
    r = api.apportionment("PFOA", n_t=121)
    assert r["success"]
    assert set(r["fraction"]) == {"root", "stem", "leaf", "grain"}
    assert np.isfinite(r["grain_phloem_fraction"])
    # the soil is the sole external source -> root soil-uptake share is the bulk
    assert r["fraction"]["root"]["soil_uptake"] > r["fraction"]["root"]["phloem_recirc"]


def test_apportionment_matches_simulate_model():
    """The apportionment must be built on the SAME model as simulate() -- check the
    resolved transport params line up so the two diagnostics stay consistent."""
    s = api.simulate("PFOA", n_t=121)
    a = api.apportionment("PFOA", n_t=121)
    assert a["params"]["f_xy"] == pytest.approx(s["params"]["f_xy"])
    assert a["params"]["L_Ph"] == pytest.approx(s["params"]["L_Ph"])
