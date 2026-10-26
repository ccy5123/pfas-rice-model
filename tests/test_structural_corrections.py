"""Corrections made under the no-unsourced-value rule (docs/PREREG_sourced_values.md section 3)
that need no new number: each organ's own cell-wall K, the inverse's error model, the even
(and inert) stem draw-off split, and the unused observation sigma."""
import os
import sys

import numpy as np
import pytest

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "src"))

import model_api as api  # noqa: E402


def test_each_organ_binds_with_its_own_cell_wall_k():
    c = api._CONG["PFOA"]
    own = api.kcw_by_organ(c)
    assert own == {o: float(c["K_cw_wholecw_Lkg"][o]) for o in api.TISSUES}
    assert len(set(own.values())) > 1                       # the organs really differ
    assert set(api.kcw_by_organ(c, "stem").values()) == {float(c["K_cw_wholecw_Lkg"]["stem"])}
    for shoot in ("redistributed", "4pool"):
        r = api.simulate("PFOA", shoot=shoot)
        assert r["params"]["K_cw"] == own


def test_binding_factor_uses_the_compartment_k():
    from pfas_rice_plant_module_4pool_surf import Compartment, Compound, binding_factors
    cmpd = Compound(name="x", K_prot=10.0, K_PL=100.0, K_cw=5.0, kappa_d=1.0, Vmax_in=0.0,
                    Km_in=1.0, Vmax_out=0.0, Km_out=1.0, L_Ph=0.1, f_xy=0.1)
    a = Compartment("a", 0.8, 0.1, 0.01, 0.5)
    b = Compartment("b", 0.8, 0.1, 0.01, 0.5, K_cw=20.0)
    Ba, Bb = binding_factors([a, b], cmpd)
    assert Bb - Ba == pytest.approx(0.2 * 0.5 * (20.0 - 5.0))


def test_stem_draw_off_split_is_even_and_inert_for_every_reported_number(monkeypatch):
    """No source gives how the stalk's transpiration divides with height, so default_tau
    splits it evenly; with retention 1 nothing reported depends on the split (segments export
    nothing): a height gradient moves no reported number at reporting precision (0.0005 log10)."""
    import pfas_rice_plant_module_nstem_leaf as NL
    assert np.allclose(NL.default_tau(4, 0.29), 0.29 / 4)
    base = api.simulate("PFOA")["baf_final_dw"]
    monkeypatch.setattr(NL, "default_tau",
                        lambda N, f: np.linspace(1.3, 0.7, N) / np.linspace(1.3, 0.7, N).sum() * f)
    grad = api.simulate("PFOA")["baf_final_dw"]
    for k in ("root", "stem", "leaf", "straw", "grain"):     # below reporting precision
        assert abs(np.log10(grad[k] / base[k])) < 0.0005


def test_inverse_error_is_model_error_plus_the_users_measurement_error():
    r = api.estimate_exposure_bayesian("PFOA", {"grain": 0.349})
    assert r["sigma_log10"] == pytest.approx(api.APRIORI_LOG10_RMSE)
    r2 = api.estimate_exposure_bayesian("PFOA", {"grain": 0.349}, meas_sigma_log10=0.2)
    assert r2["sigma_log10"] == pytest.approx(np.hypot(0.2, api.APRIORI_LOG10_RMSE))
    assert api.NEUTRAL_APRIORI_LOG10_RMSE == pytest.approx(
        np.median(list(api.NEUTRAL_APRIORI_OOS_SCORES.values())))


def test_observations_carry_no_assumed_sigma():
    d = api.yamazaki_obs("PFOA")
    assert "sigma" not in d
    with open(os.path.join(_ROOT, "data_obs", "obs_baf_Yamazaki.csv")) as f:
        header = next(line for line in f if not line.startswith("#"))
    assert "sigma" not in header


def test_per_area_biomass_needs_the_fields_density():
    import measured_biomass as mb
    with pytest.raises(ValueError):
        mb.to_kg_per_hill(np.array([100.0]), "g/m2")
    assert mb.to_kg_per_hill(np.array([100.0]), "g/m2", hills_per_m2=20.0)[0] == pytest.approx(0.005)
