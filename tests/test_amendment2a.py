"""Amendment 2a (docs/PREREG_sourced_values.md section 7): organ water by development stage and the
brown-rice grain with a husk + rachis pool. These pin the structure, not a fitted number:
the stage anchors are the sourced midpoints; fresh mass and binding read one theta(t); the
vectorised system equals the right-hand side; mass is conserved with the extra pool; brown rice
takes no xylem; its pre-harvest water moves no dry-basis number; the panicle split sums to WSO."""
import numpy as np
import pytest

import linear_expm as LE
import model_api as api
import organ_mass as om
import panicle as pn
import pfas_rice_plant_module_4pool_surf as P4
import pfas_rice_plant_module_nstem_leaf as NL
import site_scenario as ss


def _captured(shoot, fn=None, **kw):
    cap = []
    cls = NL.NStemLeafModel if shoot == "redistributed" else P4.RiceUptakeModel
    orig = cls.solve

    def grab(self, t, C0=None):
        cap.append(self)
        return orig(self, t, C0)
    cls.solve = grab
    try:
        (fn or api.simulate)("PFOS", **({} if fn else {"shoot": shoot}), **kw)
    finally:
        cls.solve = orig
    return cap[0]


def test_stage_anchors_are_the_sourced_midpoints():
    for organ, a in om.THETA_STAGE.items():
        mid = 0.5 * (a["low"] + a["high"])
        assert np.allclose(om.theta_at(organ, a["dvs"]), mid)
        assert om.theta_at(organ, 0.0) == pytest.approx(mid[0])       # flat outside
        assert om.theta_at(organ, 3.0) == pytest.approx(mid[-1])
    assert "grain" not in om.THETA_STAGE                              # brown rice: not staged
    assert om.theta_at("root", 2.0) == pytest.approx(0.6365)


def test_fresh_mass_and_binding_read_one_theta():
    d = api.scenario_drivers()
    run = api.site_run()
    th = d["theta"]
    brown, husk = pn.split(run["t"], run["grain"], run["dvs"])
    assert np.allclose(d["husk"], husk) and np.allclose(brown + husk, run["grain"])
    for i, o in enumerate(api.TISSUES):
        dry = brown if o == "grain" else run[o]
        assert np.allclose(d["M"][:, i], np.maximum(dry / (1.0 - th[:, i]), om.MASS_FLOOR))
    m = _captured("4pool")
    for j in (5, 60, 110):
        B = m._binding(m.inputs.extras_at(float(d["t"][j]))[0])
        assert np.allclose(B, th[j] + (1.0 - th[j]) * P4.dry_binding(m.comps, m.cmpd))


@pytest.mark.parametrize("shoot,kw", [("redistributed", {}), ("4pool", {}),
                                      ("redistributed", {"husk_phloem": True}),
                                      ("4pool", {"husk_phloem": True})])
def test_affine_equals_the_rhs_with_the_husk_pool(shoot, kw):
    m = _captured(shoot, **kw)
    assert m.husk and m._staged
    t = m.inputs.t
    ts = np.linspace(t[0], t[-1], 241)
    A, b = m.affine(ts)
    for j in range(0, len(ts), 6):
        A2, b2 = LE.affine_parts(m.rhs, ts[j], A.shape[1])
        assert np.allclose(A[j], A2, rtol=1e-12, atol=1e-12 * max(1.0, np.abs(A2).max()))
        assert np.allclose(b[j], b2, rtol=1e-12, atol=1e-12)


def test_affine_equals_the_rhs_two_pool_root_with_husk():
    m = _captured("redistributed", fn=api.simulate_nstem_leaf, k_seq=0.02, k_rel=0.003)
    assert m.two_pool and m.husk and m.HUSK == m.SEQ + 1
    ts = np.linspace(m.inputs.t[0], m.inputs.t[-1], 121)
    A, b = m.affine(ts)
    for j in range(0, len(ts), 10):
        A2, b2 = LE.affine_parts(m.rhs, ts[j], A.shape[1])
        assert np.allclose(A[j], A2, rtol=1e-12, atol=1e-12 * max(1.0, np.abs(A2).max()))


@pytest.mark.parametrize("shoot", ["redistributed", "4pool"])
def test_mass_is_conserved_with_the_husk_pool(shoot):
    """Shoot burden at harvest (the husk per kg DRY included) + senescence loss == xylem export
    from the root - phloem return to it. (The root side is not checked this way: with the stiff
    linear entry the net soil uptake is a small difference of large terms, and its trapezoid along
    the solution is noise; the shoot side carries the new pool.) Tolerance: the daily trapezoid."""
    r = api.apportionment("PFOA", shoot=shoot)
    s = api.simulate("PFOA", shoot=shoot)
    husk = api.scenario_drivers()["husk"]
    shoot_burden = sum(s["conc"][o][-1] * s["M"][-1, i]
                       for i, o in enumerate(api.TISSUES) if o != "root")
    shoot_burden += s["conc_dw"]["husk"][-1] * husk[-1]
    c = r["cum"]
    out = c.get("xylem_root_out", c.get("xylem_root_stem")) - c["phloem_leaf_root"]
    assert shoot_burden + c["leaf_loss"] == pytest.approx(out, rel=3e-3)


@pytest.mark.parametrize("shoot", ["redistributed", "4pool"])
def test_brown_rice_takes_no_xylem(shoot):
    fr = api.apportionment("PFOA", shoot=shoot)["fraction"]
    assert fr["grain"]["xylem_from_stem"] == 0.0 and fr["grain"]["phloem_from_leaf"] == 1.0
    assert fr["husk"]["phloem_from_leaf"] == 0.0            # only the leaf loads, only brown draws


def test_pre_harvest_grain_water_moves_no_dry_basis_number():
    """PREREG 7.2: the brown-rice equation has no term in its own free concentration, so a
    different grain water before harvest (with the fresh mass it implies) changes nothing dry."""
    d = api.scenario_drivers()
    th = d["theta"].copy()
    dry_g = d["M"][:, 3] * (1.0 - th[:, 3])
    th[:, 3] = np.interp(d["dvs"], [1.0, 2.0], [0.6, 0.21])    # wet early in filling, harvest as is
    M = d["M"].copy()
    M[:, 3] = np.maximum(dry_g / (1.0 - th[:, 3]), om.MASS_FLOOR)
    alt = dict(d, theta=th, M=M)
    for shoot in ("redistributed", "4pool"):
        a = api.simulate("PFOA", shoot=shoot, drivers=d)["baf_final_dw"]
        b = api.simulate("PFOA", shoot=shoot, drivers=alt)["baf_final_dw"]
        # analytically equal; the daily driver representation (M linear between nodes, dM by
        # central differences) leaves ~0.0012 log10 (recorded in PREREG section 8)
        assert abs(np.log10(b["grain"] / a["grain"])) < 0.002
        assert b["root"] == pytest.approx(a["root"], rel=1e-9)


def test_husk_share_follows_ansari_scaled_to_the_maturity_anchor():
    run = api.site_run()
    h = pn.husk_share(run["t"], run["dvs"])
    assert np.all(h[run["dvs"] <= 1.0] == 1.0)
    assert h[-1] == pytest.approx(pn.H_MAT)                    # harvest at DVS 2
    i10 = int(np.argmin(abs(run["t"] - (np.interp(1.0, run["dvs"], run["t"]) + 10))))
    scale = pn.H_MAT / pn.H_NIPPONBARE["maturity"]
    assert h[i10] == pytest.approx(pn.H_NIPPONBARE["10DAH"] * scale, rel=0.02)
    assert pn.H_NIPPONBARE["10DAH"] == pytest.approx(3.80 / (3.80 + 7.5))


def test_amendment_one_forms_stay_runnable():
    a = api.simulate("PFOA", organ_water="constant", grain_basis="panicle")
    assert "husk" not in a["conc_dw"] and a["params"]["grain_basis"] == "panicle"
    with pytest.raises(ValueError):
        api.simulate("PFOA", organ_water="wet")
