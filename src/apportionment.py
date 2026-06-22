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
application pulse). Our model is (a) NONLINEAR (GHK electrodiffusion + Michaelis-
Menten carrier) and (b) driven by a CONTINUOUS soil source C_w^o(t), so the
initial-source superposition does not apply. Everything here originates in the
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
    RiceUptakeModel, binding_factors, root_uptake, ROOT, STEM, LEAF, FRUIT)

_trapz = getattr(np, "trapezoid", None) or np.trapz   # numpy>=2 renames trapz

# ordered flux names (mass fluxes [ug/day])
FLUX_NAMES = (
    "soil_uptake",        # soil pore water -> root   (the sole external source)
    "xylem_root_stem",    # root -> stem  (transpiration stream, f_xy-loaded)
    "xylem_stem_leaf",    # stem -> leaf
    "xylem_stem_grain",   # stem -> grain
    "phloem_leaf_grain",  # leaf -> grain (phloem)
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
    Cwo = inp.Cwo_(t)
    Qtp = inp.Qtp_(t)
    M = np.maximum(inp.M_(t), 1e-12)
    cmpd = model.cmpd

    B = binding_factors(model.comps, cmpd)
    Cw = C / B

    A3 = model.comps[LEAF].S * M[LEAF]
    A4 = model.comps[FRUIT].S * M[FRUIT]
    split = A3 / (A3 + A4) if (A3 + A4) > 0 else 0.5
    f3, f4 = split, 1.0 - split

    Q_Phl = max(inp.dM_(t)[FRUIT] * model.T_C_Ph + model.phi * Qtp, 0.0)
    C_Phl = cmpd.L_Ph * Cw[LEAF] + cmpd.g_ph * C[LEAF]
    Cw_xyl = cmpd.f_xy * Cw[ROOT] + cmpd.g_xy * C[ROOT]
    jR = root_uptake(Cwo, Cw[ROOT], cmpd, model.env)
    gam = inp.grain_gate_(t)
    g = [c.gamma for c in model.comps]

    return dict(
        soil_uptake=jR * M[ROOT],
        xylem_root_stem=Qtp * Cw_xyl,
        xylem_stem_leaf=(f3 + (1.0 - gam) * f4) * Qtp * Cw[STEM],
        xylem_stem_grain=gam * f4 * Qtp * Cw[STEM],
        phloem_leaf_grain=gam * Q_Phl * C_Phl,
        phloem_leaf_root=model.phi * Q_Phl * C_Phl,
        leaf_loss=inp.leaf_loss_(t) * C[LEAF] * M[LEAF],
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
    M = np.maximum(model.inputs.M_(t), 1e-12)
    mu = model.inputs.dM_(t) / M
    mass_bal = np.array([
        f["soil_uptake"] - f["xylem_root_stem"] + f["phloem_leaf_root"] - f["degr_root"],
        f["xylem_root_stem"] - (f["xylem_stem_leaf"] + f["xylem_stem_grain"]) - f["degr_stem"],
        f["xylem_stem_leaf"] - (f["phloem_leaf_grain"] + f["phloem_leaf_root"])
        - f["leaf_loss"] - f["degr_leaf"],
        f["xylem_stem_grain"] + f["phloem_leaf_grain"] - f["degr_grain"],
    ])
    return mass_bal / M - mu * C


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
