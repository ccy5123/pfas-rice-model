"""
organ_mass.py — crop-model DRY matter -> the FRESH organ mass the plant ODE runs on
===================================================================================

Every biomass driver in this repo produces DRY matter [kg DM/hill]: the ORYZA2000
carbon balance (`oryza_growth`), the IR72 partitioning reconstruction (`growth_rice`)
and the measured-table ingester (`measured_biomass`, which takes organ dry weights).
That is the right unit for a crop model -- assimilation, HI and the 1740 g/m^2 anchor
are all dry weights.

The plant ODE is FRESH-weight throughout: `PlantInputs.M` is "tissue fresh mass",
the binding factor B_k = theta + (1-theta)*(f_prot*K_prot + f_PL*K_PL + f_cw*K_cw) is
per kg FRESH, and every observation is converted to fresh weight before comparison.
Until this module existed the dry matter went into M unconverted, so a terminal
organ's concentration (burden / M) was a burden per kg DRY matter reported as per kg
fresh: leaf x1/(1-0.78)=4.5, stem x5.9 and grain x1.16, against x10 for the root,
whose concentration is set by uptake per unit mass and therefore barely moves. Net
effect on the default run: shoot BAFs 4-6x high relative to the root
(PFOA leaf 2.52 -> 0.61, grain 0.148 -> 0.036; root 0.479 -> 0.488).

    fresh = dry / (1 - theta_fw)

with theta_fw the SAME per-organ water fractions the binding factor uses
(`params/parameters.json` -> `tissue_composition_recommended`, grain = `grain_brown`),
so the mass basis and the concentration basis cannot drift apart.

Driver dicts carry the tag `mass_basis="fresh"` (MASS_BASIS) so the one place that needs
dry matter -- the leaf/grain xylem split -- can take it back out (`xylem_leaf_share`).

LIMITATIONS (docs/sampling_protocols_audit_KR.md): theta_fw is held constant over
the season (grain 0.14 is harvest moisture; filling grain is ~0.5-0.6 water, mature
straw dries to 0.1-0.5), and the root/stem/leaf values are estimates -- no rice
dataset in the repo, nor any found in the literature audit, reports organ water
content over a season.
"""
from __future__ import annotations

import json
import os

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ORGANS = ("root", "stem", "leaf", "grain")
MASS_BASIS = "fresh"          # the tag a driver dict carries for an M built here
# model organ -> tissue_composition_recommended key (the grain pool is brown rice)
_COMP_KEY = {"root": "root", "stem": "stem", "leaf": "leaf", "grain": "grain_brown"}


def _load_theta():
    with open(os.path.join(_ROOT, "params", "parameters.json")) as f:
        comp = json.load(f)["tissue_composition_recommended"]
    return {o: float(comp[k]["theta_fw"]) for o, k in _COMP_KEY.items()}


THETA_FW = _load_theta()


def fresh_factor(organ):
    """Multiplier dry -> fresh organ mass: 1/(1-theta_fw)."""
    return 1.0 / (1.0 - THETA_FW[organ])


def to_model_mass(b, theta=None):
    """Convert a crop-model biomass dict {root,stem,leaf,grain: DRY kg/hill, ...} to FRESH
    organ mass. Non-organ keys (e.g. ORYZA's `leaf_death_rate`, a RELATIVE rate that a
    constant factor leaves unchanged) are passed through untouched. `theta` overrides the
    organ water {organ: theta} (the neutral validation tables run on their own plant's
    water, PREREG 5.10)."""
    th = dict(THETA_FW, **(theta or {}))
    out = dict(b)
    for o in ORGANS:
        if o in b and b[o] is not None:
            out[o] = np.asarray(b[o], float) * (1.0 / (1.0 - th[o]))
    return out


# The ODE divides by organ mass, so an organ the crop model has not formed yet (the grain
# before panicle initiation) sits on a floor. The floor is numerical, not a datum: it is
# admissible because the results do not depend on it -- harvest values are identical for
# every floor <= 1e-6 kg and the young-panicle time course for every floor <= 1e-9 kg
# (validation/numerics_convergence.py). 1e-4 kg, used before, was not converged.
MASS_FLOOR = 1e-9            # kg


def organ_matrix(b, floor=MASS_FLOOR, theta=None):
    """(n, 4) FRESH mass matrix [root, stem, leaf, grain], floored at `MASS_FLOOR`."""
    m = to_model_mass(b, theta)
    return np.maximum(np.column_stack([m[o] for o in ORGANS]), floor)
