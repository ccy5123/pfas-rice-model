"""
Source apportionment for the 4-compartment PFAS rice model
==========================================================

The PFAS analog of dynamiCROP's **Figure 2** (Pang et al. 2020, Environ. Pollut.
256:113285): for each plant compartment, decompose the cumulative PFAS mass
DELIVERED by each transport pathway over the season.

Why it differs from dynamiCROP's Fig. 2
---------------------------------------
dynamiCROP splits a compartment's residue by INITIAL SOURCE pool (air/soil/leaf-
surface/fruit-surface) using linear superposition -- valid because that model is a
linear, constant-coefficient ODE driven only by an initial condition (a single
application pulse). Our model is driven by a CONTINUOUS soil source C_w^o(t) through
the root entry, with time-varying coefficients, so the initial-source superposition does
not apply. Everything here originates in the
pore water -> root uptake j_R; from there PFAS moves

    root --(xylem)--> stem --(xylem)--> leaf, grain
    leaf --(phloem)--> grain,  leaf --(phloem)--> root (recirculation)

so the meaningful decomposition is FLUX-BASED and by transport PATHWAY. The
headline audit number is the GRAIN phloem fraction: the model asserts the grain is
phloem-fed (loading L_Ph, no weak-acid ion-trap), and this quantifies how much of
the grain burden actually arrives via phloem vs the residual xylem.

API
---
flux_terms(model, t, C)
    Instantaneous inter-compartment MASS fluxes [ug/day] at a state, mirroring
    RiceUptakeModel.rhs term-by-term. Reported at the MASS level (conc-rate * M_k),
    where the growth-dilution term -mu*C cancels (it redistributes concentration,
    not mass), so the returned fluxes conserve mass at every junction
    (root xylem-out == stem xylem-in, etc.). The tests reconstruct dC/dt from these
    and check it against rhs() to guarantee the two never diverge.
apportion(model, sol)
    Integrate the fluxes along a solved trajectory -> cumulative pathway masses and
    per-compartment inflow fractions.
"""
from __future__ import annotations
import numpy as np

from pfas_rice_plant_module_4pool_surf import (
    RiceUptakeModel, binding_factors, root_uptake, xylem_leaf_share, panicle_share,
    phloem_flows, MASS_FLOOR, ROOT, STEM, LEAF, FRUIT)

_trapz = getattr(np, "trapezoid", None) or np.trapz   # numpy>=2 renames trapz

# ordered flux names (mass fluxes [ug/day])
FLUX_NAMES = (
    "soil_uptake",        # soil pore water -> root   (the sole external source)
    "xylem_root_stem",    # root -> stem  (transpiration stream, f_xy-loaded)
    "xylem_stem_leaf",    # stem -> leaf
    "xylem_stem_grain",   # stem -> grain
    "xylem_stem_husk",    # stem -> husk + rachis (PREREG 7.2: the panicle's xylem; 0 without it)
    "phloem_leaf_grain",  # leaf -> grain (phloem)
    "phloem_leaf_husk",   # leaf -> husk (a band end of 7.2; 0 by default)
    "phloem_leaf_root",   # leaf -> root  (phloem recirculation)
    "leaf_loss",          # leaf -> out   (senescence / shedding, ORYZA driver only)
    "degr_root", "degr_stem", "degr_leaf", "degr_grain",   # metabolism (PFAS ~ 0)
)


def flux_terms(model: RiceUptakeModel, t: float, C: np.ndarray) -> dict:
    """Instantaneous inter-compartment MASS fluxes [ug/day] at state (t, C).

    Mirrors RiceUptakeModel.rhs term-by-term but reports each transfer as a mass
    flux (the conc-rate term times M_k). At the mass level the growth-dilution
    term -mu*C cancels exactly, so these fluxes conserve mass at every junction.
    """
    inp = model.inputs
    Cwo, Qtp, M, dM, leaf_loss, gam, dvs = inp.drivers_at(t)     # the drivers rhs() reads
    th, gdry, Mh, dMh, gh = inp.extras_at(t)
    M = np.maximum(M, MASS_FLOOR)
    cmpd = model.cmpd

    B = model._binding(th)
    Cw = C[:4] / B

    if model.lam_grain_dvs is not None:                          # the split rhs() uses
        f4 = panicle_share(dvs, *model.lam_grain_dvs)
        f3 = 1.0 - f4
    else:
        split = xylem_leaf_share(model.comps, M, inp.M_basis, model.s_grain_dry)
        f3, f4 = split, 1.0 - split

    Q_g, Q_r = phloem_flows(model.comps, dM, Qtp, model.T_C_Ph, model.phi, inp.M_basis)
    if model._staged:
        Q_g = max(gdry * model.T_C_Ph, 0.0)
    Q_h = max(dMh * model.T_C_Ph, 0.0) * gh if model.husk_phloem else 0.0
    C_Phl = model.phloem_loading_factor() * Cw[LEAF] + cmpd.g_ph * C[LEAF]
    Cw_xyl = cmpd.f_xy * Cw[ROOT] + cmpd.g_xy * C[ROOT]
    jR = root_uptake(Cwo, Cw[ROOT], cmpd, model.env)
    g = [c.gamma for c in model.comps]
    wx = gh if model.husk else gam

    return dict(
        soil_uptake=jR * M[ROOT],
        xylem_root_stem=Qtp * Cw_xyl,
        xylem_stem_leaf=(f3 + (1.0 - wx) * f4) * Qtp * Cw[STEM],
        xylem_stem_grain=0.0 if model.husk else gam * f4 * Qtp * Cw[STEM],
        xylem_stem_husk=wx * f4 * Qtp * Cw[STEM] if model.husk else 0.0,
        phloem_leaf_husk=Q_h * C_Phl,
        phloem_leaf_grain=gam * Q_g * C_Phl,
        phloem_leaf_root=Q_r * C_Phl,
        leaf_loss=leaf_loss * C[LEAF] * M[LEAF],
        degr_root=g[ROOT] * C[ROOT] * M[ROOT],
        degr_stem=g[STEM] * C[STEM] * M[STEM],
        degr_leaf=g[LEAF] * C[LEAF] * M[LEAF],
        degr_grain=g[FRUIT] * C[FRUIT] * M[FRUIT],
    )


def dC_from_fluxes(model: RiceUptakeModel, t: float, C: np.ndarray) -> np.ndarray:
    """Reconstruct dC/dt from the mass fluxes (consistency mirror of rhs).

    d(m_k)/dt = [mass inflow - outflow]_k, and dC_k/dt = d(m_k)/dt / M_k - mu_k*C_k
    (the -mu*C growth-dilution term is re-added at the concentration level). Used by
    the tests to assert flux_terms() never diverges from RiceUptakeModel.rhs().
    """
    f = flux_terms(model, t, C)
    _, _, M, dM, _, _, _ = model.inputs.drivers_at(t)
    _, _, Mh, dMh, _ = model.inputs.extras_at(t)
    M = np.maximum(M, MASS_FLOOR)
    mu = dM / M
    mass_bal = np.array([
        f["soil_uptake"] - f["xylem_root_stem"] + f["phloem_leaf_root"] - f["degr_root"],
        f["xylem_root_stem"] - (f["xylem_stem_leaf"] + f["xylem_stem_grain"]
                                + f["xylem_stem_husk"]) - f["degr_stem"],
        f["xylem_stem_leaf"] - (f["phloem_leaf_grain"] + f["phloem_leaf_root"]
                                + f["phloem_leaf_husk"]) - f["leaf_loss"] - f["degr_leaf"],
        f["xylem_stem_grain"] + f["phloem_leaf_grain"] - f["degr_grain"],
    ])
    dC = mass_bal / M - mu * C[:4]
    if model.husk:
        dC = np.append(dC, (f["xylem_stem_husk"] + f["phloem_leaf_husk"]) / Mh
                       - dMh / Mh * C[model.HUSK])
    return dC


def apportion(model: RiceUptakeModel, sol) -> dict:
    """Integrate inter-compartment mass fluxes along a solved trajectory.

    Returns
    -------
    dict with
      cum      : {flux_name: cumulative mass over the season [ug]}
      inflow   : per-compartment cumulative inflow by pathway [ug]
      fraction : per-compartment inflow fractions (sum to 1 per compartment)
      delivered_to_shoot : fraction of cumulative soil uptake that left the root
                           (xylem root->stem), i.e. the realised root->shoot transfer
      grain_phloem_fraction : the headline audit number (phloem share of grain inflow)
    """
    t = sol.t
    series = {name: np.empty(len(t)) for name in FLUX_NAMES}
    for i, ti in enumerate(t):
        f = flux_terms(model, float(ti), sol.y[:, i])
        for name in FLUX_NAMES:
            series[name][i] = f[name]
    cum = {name: float(_trapz(series[name], t)) for name in FLUX_NAMES}

    inflow = {
        "root":  {"soil_uptake": cum["soil_uptake"], "phloem_recirc": cum["phloem_leaf_root"]},
        "stem":  {"xylem_from_root": cum["xylem_root_stem"]},
        "leaf":  {"xylem_from_stem": cum["xylem_stem_leaf"]},
        "grain": {"xylem_from_stem": cum["xylem_stem_grain"],
                  "phloem_from_leaf": cum["phloem_leaf_grain"]},
    }
    if model.husk:
        inflow["husk"] = {"xylem_from_stem": cum["xylem_stem_husk"],
                          "phloem_from_leaf": cum["phloem_leaf_husk"]}
    fraction = {}
    for k, d in inflow.items():
        tot = sum(max(v, 0.0) for v in d.values())
        fraction[k] = {p: (max(v, 0.0) / tot if tot > 0 else float("nan"))
                       for p, v in d.items()}

    gtot = inflow["grain"]["xylem_from_stem"] + inflow["grain"]["phloem_from_leaf"]
    grain_phloem = inflow["grain"]["phloem_from_leaf"] / gtot if gtot > 0 else float("nan")
    delivered = cum["xylem_root_stem"] / cum["soil_uptake"] if cum["soil_uptake"] > 0 else float("nan")
    return dict(t_span=(float(t[0]), float(t[-1])), cum=cum, inflow=inflow,
                fraction=fraction, grain_phloem_fraction=grain_phloem,
                delivered_to_shoot=delivered)


# ---------------------------------------------------------------------------
# demo: run the apportionment audit across a few congeners
# ---------------------------------------------------------------------------
def _demo():
    import model_api as api
    print("Source apportionment audit (Fig-2 analog) -- cumulative pathway delivery\n")
    print(f"{'congener':8s} {'n_C':>3s} {'grain:xylem':>11s} {'grain:phloem':>12s} "
          f"{'root->shoot':>11s}")
    for cong in ("PFBA", "PFOA", "PFOS", "PFDoDA"):
        if cong not in api._CONG:
            continue
        r = api.apportionment(cong)
        fr_g = r["fraction"]["grain"]
        print(f"{cong:8s} {api._CONG[cong]['n_C']:>3d} "
              f"{fr_g['xylem_from_stem']:>11.3f} {fr_g['phloem_from_leaf']:>12.3f} "
              f"{r['delivered_to_shoot']:>11.4f}")
    print("\n(grain phloem fraction near 1 => phloem-fed, as the model asserts;"
          "\n root->shoot = fraction of soil uptake that leaves the root via xylem)")


if __name__ == "__main__":
    _demo()


# ---------------------------------------------------------------------------
# The redistributed shoot (pfas_rice_plant_module_nstem_leaf.NStemLeafModel)
# ---------------------------------------------------------------------------
def flux_terms_redistributed(model, t, C):
    """Instantaneous MASS fluxes [ug/day] of the redistributed-shoot model at (t, C),
    term by term as in NStemLeafModel.rhs. The xylem is not passed stem -> leaf -> grain
    here: each organ draws its own transpiration share, so the grain's xylem inflow is the
    panicle draw plus the non-retained residual (formation-gated), and before the panicle
    forms that residual goes to the leaf."""
    from pfas_rice_plant_module_4pool_surf import root_uptake as _ru
    inp = model.inputs
    Cwo, Q, M, dM, leaf_loss, gam, dvs = inp.drivers_at(t)       # the drivers rhs() reads
    th, gdry, Mh, dMh, gh = inp.extras_at(t)
    M = np.maximum(M, MASS_FLOOR)
    n = model.n_comp
    B = th + (1.0 - th) * model._X if model._staged else model._B
    Cw = C[:n] / B
    leaf, grain, r = model.LEAF, model.GRAIN, model.retention
    if not model._staged:
        gdry = dM[grain] * (1.0 - model.comps[grain].theta)
    Q_g = max(gdry * model.T_C_Ph, 0.0)
    Q_r = max(model.phi * Q, 0.0)
    Q_h = max(dMh * model.T_C_Ph, 0.0) * gh if model.husk_phloem else 0.0
    C_Phl = model.cmpd.L_Ph * Cw[leaf] + model.cmpd.g_ph * C[leaf]
    jR = _ru(Cwo, Cw[model.ROOT], model.cmpd, model.env)
    xyl = Q * (model.cmpd.f_xy * Cw[model.ROOT] + model.cmpd.g_xy * C[model.ROOT])
    if model.lam_grain_dvs is not None:                          # the split rhs() uses
        lg = panicle_share(dvs, *model.lam_grain_dvs)
        ll = 1.0 - model.tau.sum() - lg
    else:
        lg, ll = model.lam_grain, model.lam_leaf
    residual = (lg + (1.0 - r) * (model.tau.sum() + ll)) * xyl
    wx = gh if model.husk else gam
    return dict(
        soil_uptake=float(M[model.ROOT] * jR),
        xylem_root_out=float(xyl),
        xylem_to_stem=float(r * model.tau.sum() * xyl),
        xylem_to_leaf=float(r * ll * xyl + (1.0 - wx) * residual),
        xylem_to_grain=0.0 if model.husk else float(gam * residual),
        xylem_to_husk=float(wx * residual) if model.husk else 0.0,
        phloem_leaf_husk=float(Q_h * C_Phl),
        phloem_leaf_grain=float(gam * Q_g * C_Phl),
        phloem_leaf_root=float(Q_r * C_Phl),
        leaf_loss=float(leaf_loss * M[leaf] * C[leaf]),
    )


def apportion_redistributed(model, sol) -> dict:
    """`apportion` for the redistributed shoot: cumulative pathway masses and per-compartment
    inflow fractions, with the same keys for the root and the grain as the 4-pool version."""
    t = sol.t
    names = ("soil_uptake", "xylem_root_out", "xylem_to_stem", "xylem_to_leaf",
             "xylem_to_grain", "xylem_to_husk", "phloem_leaf_grain", "phloem_leaf_husk",
             "phloem_leaf_root", "leaf_loss")
    series = {k: np.empty(len(t)) for k in names}
    for i, ti in enumerate(t):
        f = flux_terms_redistributed(model, float(ti), sol.y[:, i])
        for k in names:
            series[k][i] = f[k]
    cum = {k: float(_trapz(series[k], t)) for k in names}
    inflow = {
        "root":  {"soil_uptake": cum["soil_uptake"], "phloem_recirc": cum["phloem_leaf_root"]},
        "stem":  {"xylem_from_root": cum["xylem_to_stem"]},
        "leaf":  {"xylem_from_root": cum["xylem_to_leaf"]},
        "grain": {"xylem_from_stem": cum["xylem_to_grain"],
                  "phloem_from_leaf": cum["phloem_leaf_grain"]},
    }
    if model.husk:
        inflow["husk"] = {"xylem_from_stem": cum["xylem_to_husk"],
                          "phloem_from_leaf": cum["phloem_leaf_husk"]}
    fraction = {}
    for k, d in inflow.items():
        tot = sum(max(v, 0.0) for v in d.values())
        fraction[k] = {p: (max(v, 0.0) / tot if tot > 0 else float("nan")) for p, v in d.items()}
    gtot = cum["xylem_to_grain"] + cum["phloem_leaf_grain"]
    return dict(t_span=(float(t[0]), float(t[-1])), cum=cum, inflow=inflow, fraction=fraction,
                grain_phloem_fraction=(cum["phloem_leaf_grain"] / gtot if gtot > 0 else float("nan")),
                delivered_to_shoot=(cum["xylem_root_out"] / cum["soil_uptake"]
                                    if cum["soil_uptake"] > 0 else float("nan")))

