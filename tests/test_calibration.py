"""Tests for the Tier-1 calibration machinery."""
import numpy as np
import pytest

from pfas_rice_plant_module import (
    Environment, Compound, Compartment, PlantInputs, RiceUptakeModel,
    _logistic, ROOT,
)
import calibration as cal
from calibration import (
    Param, ObservedBAF, predict_bafs, set_param, get_param, calibrate,
    load_baf_csv,
)


def build_model(**cmpd_kw):
    t = np.linspace(0.0, 120.0, 481)
    Cwo = np.full_like(t, 1.0)
    Qtp = 0.05 + 0.35 * np.exp(-((t - 75.0) ** 2) / (2 * 25.0 ** 2))
    M = np.column_stack([
        _logistic(t, 1e-3, 0.030, 0.10, 20.0), _logistic(t, 1e-3, 0.040, 0.10, 25.0),
        _logistic(t, 1e-3, 0.050, 0.12, 30.0), _logistic(t, 1e-5, 0.025, 0.18, 80.0)])
    base = dict(name="PFOA", K_prot=50.0, K_PL=100.0, K_cw=20.0, kappa_d=0.5,
                Vmax_in=20.0, Km_in=5.0, Vmax_out=8.0, Km_out=5.0, L_Ph=0.005, f_xy=0.02)
    base.update(cmpd_kw)
    comps = [Compartment("root", 0.70, 0.05, 0.010, 0.30),
             Compartment("stem", 0.80, 0.01, 0.005, 0.08),
             Compartment("leaf", 0.80, 0.03, 0.020, 0.04, S=20.0),
             Compartment("grain", 0.15, 0.08, 0.010, 0.10, S=2.0)]
    return RiceUptakeModel(env=Environment(), cmpd=Compound(**base), comps=comps,
                           inputs=PlantInputs(t=t, Cwo=Cwo, Qtp=Qtp, M=M)), t


# ---------------------------------------------------------------------------
def test_param_log_and_linear_mapping():
    p = Param("f_xy", 1e-3, 1.0, log=True)
    assert p.to_value(p.to_x(0.02)) == pytest.approx(0.02)
    assert p.bounds_x == (pytest.approx(-3.0), pytest.approx(0.0))
    q = Param("phi", 0.0, 0.5, log=False)
    assert q.to_value(q.to_x(0.1)) == pytest.approx(0.1)


def test_set_get_param_roundtrip():
    model, _ = build_model()
    set_param(model, "f_xy", 0.07)
    assert get_param(model, "f_xy") == pytest.approx(0.07)
    set_param(model, "phi", 0.2)
    assert model.phi == pytest.approx(0.2)
    with pytest.raises(KeyError):
        set_param(model, "not_a_param", 1.0)


def test_observed_baf_validates_tissue():
    with pytest.raises(ValueError):
        ObservedBAF("xylem", 1.0)


def test_predict_bafs_keys_and_ordering():
    model, t = build_model()
    b = predict_bafs(model, t)
    assert set(b) == {"root", "stem", "leaf", "straw", "grain"}
    assert b["root"] > b["straw"] > b["grain"]      # demo ordering


def test_calibrate_fits_and_does_not_mutate_input():
    model, t = build_model()
    # observations generated from a *different* parameter set
    truth = {"f_xy": 0.03, "L_Ph": 0.01, "kappa_d": 0.6}
    gen, _ = build_model(f_xy=truth["f_xy"], L_Ph=truth["L_Ph"], kappa_d=truth["kappa_d"])
    obs = [ObservedBAF(tis, v) for tis, v in predict_bafs(gen, t).items()
           if tis in ("root", "straw", "grain")]
    params = [Param("f_xy", 1e-3, 1.0), Param("L_Ph", 1e-4, 0.5), Param("kappa_d", 1e-2, 5.0)]
    before = get_param(model, "f_xy")
    res = calibrate(model, params, obs, x0=[0.02, 0.005, 0.5])   # local, fast
    # input model is untouched (calibration works on a copy)
    assert get_param(model, "f_xy") == before
    # the fit reproduces the observations closely
    for tis in ("root", "straw", "grain"):
        assert res.predicted[tis] == pytest.approx(res.observed[tis], rel=0.02)


def test_synthetic_recovery_local():
    """Well-posed 3-param/3-obs recovery converges to the truth (local solver)."""
    truth, res = cal.synthetic_recovery(noise=0.0, global_search=False)
    assert res.success
    assert res.values["f_xy"] == pytest.approx(truth["f_xy"], rel=0.10)
    assert res.values["L_Ph"] == pytest.approx(truth["L_Ph"], rel=0.15)
    assert res.values["kappa_d"] == pytest.approx(truth["kappa_d"], rel=0.15)


def test_synthetic_recovery_with_noise():
    """Under 10% observation noise the recovered Tier-1 params stay in the right
    ballpark (within a factor of two of the truth)."""
    truth, res = cal.synthetic_recovery(noise=0.10, seed=3, global_search=False)
    assert res.success
    for k in ("f_xy", "L_Ph", "kappa_d"):
        assert 0.5 * truth[k] < res.values[k] < 2.0 * truth[k]


def test_load_baf_csv(tmp_path):
    p = tmp_path / "baf.csv"
    p.write_text("compound,tissue,baf,sigma\n"
                 "PFOA,root,9.0,0.2\nPFOA,straw,4.3,0.3\nPFOA,grain,2.3,0.3\n")
    obs = load_baf_csv(str(p))
    assert len(obs) == 3
    assert obs[0].tissue == "root" and obs[0].value == pytest.approx(9.0)
    assert obs[0].sigma == pytest.approx(0.2)
    assert obs[2].tissue == "grain"
