"""
What a harvest measures: the observation operator from the plant ODE to a tissue sample.

The ODE carries GREEN leaves. When a crop driver kills leaves (ORYZA's leaf death rate, passed
as `leaf_loss`), the leaf ODE loses `leaf_loss * C_leaf` -- the dead blade takes its PFAS with
it. That is right for the green-leaf pool, but a harvested shoot still holds those blades: rice
leaves dry on the culm rather than abscise, ORYZA counts them in the above-ground biomass
(WAGT = WLVG + WLVD + WST + WSO), and every dataset the model is scored on sampled "stems
including leaves" (Yamazaki 2023), "leaves" (Tang 2026, Li 2025) or straw, without separating
green from dead. Dropping the dead pool made the model's straw and leaf a GREEN-tissue
concentration compared against a green + dead measurement.

Here the dead pool is reconstructed exactly from the solve: its burden accumulates at the rate
the leaf ODE loses it, `leaf_loss * C_leaf * M_leaf`, and its dry mass at `leaf_loss * M_leaf`
times the leaf's dry fraction. Dead tissue takes part in no further exchange (no uptake, no
export, no wash-off): it keeps the concentration the blade had when it died.

`attached` is the fraction of dead blades still on the plant at harvest. It is NOT measured for
any dataset here. 1.0 (default) is ORYZA's own convention and the reading of a whole-plant
harvest; 0.0 reproduces the old convention (dead blades gone). Kim 2003's senescent share at
maturity (0.36-0.47) sits below the port's (0.70-0.79), which is consistent either with blades
lost before sampling or with the port killing leaves too fast -- so both ends are run as a
sensitivity wherever a harvest is scored.

Concentrations against dry-weight data (all PFAS tissue tables in use report ng/g dw) are formed
on the dry basis: burden / dry mass. That needs the dry fraction of each LIVING organ (the
binding factor's theta_fw) but no water content of dead tissue. The fresh-basis straw, kept for
display, counts dead blades at their dry mass (a dead blade's water is not modelled).
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import cumulative_trapezoid

ORGANS = ("root", "stem", "leaf", "grain")
REFINE_PER_DAY = 24          # sub-daily grid for the dead-pool integral when a dense solution exists


def dead_leaf_pool(t, C_leaf, M_leaf, leaf_loss, dry_fraction_leaf, C_leaf_fn=None,
                   M_leaf_fn=None, leaf_loss_fn=None):
    """Cumulative dead-leaf burden [µg/hill] and dry mass [kg/hill] on the grid `t`.

    `C_leaf`, `M_leaf` and `leaf_loss` are arrays on `t`; the optional `*_fn` callables (dense
    solution, driver interpolants) let the integral run on a sub-daily grid."""
    t = np.asarray(t, float)
    if leaf_loss is None:
        z = np.zeros_like(t)
        return z, z.copy()
    if C_leaf_fn is not None and len(t) > 1:
        n = max(int(np.ceil((t[-1] - t[0]) * REFINE_PER_DAY)), len(t))
        tf = np.union1d(np.linspace(t[0], t[-1], n + 1), t)
        cl = np.asarray(C_leaf_fn(tf), float)
        ml = np.asarray(M_leaf_fn(tf), float) if M_leaf_fn is not None else np.interp(tf, t, M_leaf)
        ll = (np.asarray(leaf_loss_fn(tf), float) if leaf_loss_fn is not None
              else np.interp(tf, t, np.asarray(leaf_loss, float)))
    else:
        tf, cl, ml, ll = t, np.asarray(C_leaf, float), np.asarray(M_leaf, float), np.asarray(leaf_loss, float)
    ll = np.maximum(ll, 0.0)
    burden = cumulative_trapezoid(ll * cl * ml, tf, initial=0.0)
    dry = cumulative_trapezoid(ll * ml * dry_fraction_leaf, tf, initial=0.0)
    return np.interp(t, tf, burden), np.interp(t, tf, dry)


def harvest_observables(t, conc, M, dry_fraction, leaf_loss=None, attached=1.0,
                        C_leaf_fn=None, M_leaf_fn=None, leaf_loss_fn=None):
    """Dry-basis tissue concentrations as a harvest sample measures them.

    conc          dict root/stem/leaf/grain -> concentration series [µg/kg on M's basis]
    M             (n, 4) organ mass [root, stem, leaf, grain] on the ODE's basis
    dry_fraction  dict organ -> dry matter per unit of M (1 - theta_fw)

    Returns dict with conc_dw {root, stem, leaf_green, leaf, grain, straw}, straw_fw (dead blades
    at their dry mass), dead_burden, dead_dry, attached."""
    M = np.asarray(M, float)
    a = float(attached)
    if not 0.0 <= a <= 1.0:
        raise ValueError("attached must lie in [0, 1]")
    burden = {o: np.asarray(conc[o], float) * M[:, i] for i, o in enumerate(ORGANS)}
    dry = {o: M[:, i] * dry_fraction[o] for i, o in enumerate(ORGANS)}
    db, dd = dead_leaf_pool(t, conc["leaf"], M[:, 2], leaf_loss, dry_fraction["leaf"],
                            C_leaf_fn, M_leaf_fn, leaf_loss_fn)
    db, dd = a * db, a * dd
    conc_dw = {o: burden[o] / dry[o] for o in ("root", "stem", "grain")}
    conc_dw["leaf_green"] = burden["leaf"] / dry["leaf"]
    conc_dw["leaf"] = (burden["leaf"] + db) / (dry["leaf"] + dd)
    conc_dw["straw"] = (burden["stem"] + burden["leaf"] + db) / (dry["stem"] + dry["leaf"] + dd)
    straw_fw = (burden["stem"] + burden["leaf"] + db) / (M[:, 1] + M[:, 2] + dd)
    return dict(conc_dw=conc_dw, straw_fw=straw_fw, dead_burden=db, dead_dry=dd, attached=a)


def dry_fractions(theta, M_basis="fresh"):
    """Dry matter per unit of M for each organ: 1 - theta_fw. (The basis-free placeholder,
    M_basis None, is read as fresh too.)"""
    return {o: 1.0 - float(theta[o]) for o in ORGANS}
