"""The exponential integrator (src/linear_expm.py): the plant ODE is affine in its state, so
each model's vectorised `affine(ts)` must equal the (A, b) read off its own `rhs`, and the
stepped solution must agree with BDF at reporting precision."""
import numpy as np
import pytest
from scipy.linalg import expm

import linear_expm as LE
import model_api as api
import pfas_rice_plant_module_4pool_surf as P4
import pfas_rice_plant_module_nstem_leaf as NL

KW = dict(f_xy_override=0.054, g_xy_override=0.0049, g_ph_override=0.001, L_Ph_override=0.47,
          r_eq_override=0.09, lipid_loading=False)


def _captured(shoot, fn=None, **extra):
    """The model object simulate() (or `fn`) builds, captured at solve."""
    cap = []
    cls = NL.NStemLeafModel if shoot == "redistributed" else P4.RiceUptakeModel
    orig = cls.solve

    def grab(self, t, C0=None):
        cap.append((self, t))
        return orig(self, t, C0)
    cls.solve = grab
    try:
        if fn is None:
            api.simulate("PFOS", shoot=shoot, **{**KW, **extra})
        else:
            fn("PFOS", **extra)
    finally:
        cls.solve = orig
    return cap[0]


@pytest.mark.parametrize("shoot,two_pool", [("redistributed", False), ("4pool", False),
                                            ("redistributed", True)])
def test_affine_equals_the_rhs(shoot, two_pool):
    if two_pool:                                   # the sequestered root pool (nstem_leaf)
        m, t = _captured(shoot, fn=api.simulate_nstem_leaf, k_seq=0.02, k_rel=0.003, **KW)
        assert m.two_pool
    else:
        m, t = _captured(shoot)
    ts = np.linspace(t[0], t[-1], 301)
    A, b = m.affine(ts)
    for j in range(0, len(ts), 5):
        A2, b2 = LE.affine_parts(m.rhs, ts[j], A.shape[1])
        assert np.allclose(A[j], A2, rtol=1e-12, atol=1e-12 * max(1.0, np.abs(A2).max()))
        assert np.allclose(b[j], b2, rtol=1e-12, atol=1e-12)


def test_expm_batch_matches_scipy():
    rng = np.random.default_rng(0)
    A = rng.normal(size=(50, 9, 9)) * rng.uniform(0, 30, size=(50, 1, 1))
    A[:, np.arange(9), np.arange(9)] -= rng.uniform(0, 500, size=(50, 9))   # stiff decay
    E = LE.expm_batch(A)
    for a, e in zip(A, E):
        ref = expm(a)
        assert np.allclose(e, ref, rtol=1e-9, atol=1e-12 * max(1.0, np.abs(ref).max()))


@pytest.mark.parametrize("shoot", ["redistributed", "4pool"])
def test_expm_agrees_with_bdf(shoot):
    out = {}
    for solver in ("expm", "bdf"):
        P4.SOLVER = solver
        try:
            r = api.simulate("PFOA", shoot=shoot, **KW)
        finally:
            P4.SOLVER = "expm"
        out[solver] = r
    for k in ("root", "straw", "grain"):
        d = abs(np.log10(out["expm"]["baf_final_dw"][k] / out["bdf"]["baf_final_dw"][k]))
        assert d < 5e-4, (k, d)                    # the reporting precision


def test_a_carrier_falls_back_to_bdf():
    m, t = _captured("4pool")
    m.cmpd.Vmax_in = 1.0
    assert m.affine(t[:3]) is None
    assert LE.root_coefficients(P4.root_uptake, m.cmpd, m.env) is None
