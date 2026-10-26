# =============================================================================
# neutral_composition.py -- the neutral path's composition, per plant, from sources
# -----------------------------------------------------------------------------
# PREREG 5.10. K_PW = W + L*a*Kow^b needs the organ's water W and its FRESH lipid L on the same
# basis. Lipid is carried here on a DRY basis (a source's fresh lipid turned dry with that
# source's own water) and restated at the organ's water: L_fw = L_dw * (1 - W). Each
# validation table runs on its own plant's water and lipid
# (params/neutral_composition_sourced.csv); an organ a table gives no value for takes the rice
# compartments' (water = organ_mass.THETA_FW, lipid = the "rice" rows).
#
# A range runs at its midpoint; `corner` moves every range of one quantity to an end
# ("low" / "high"), so a score can be reported with its envelope (PREREG 5.4).
# =============================================================================
from __future__ import annotations

import csv
import itertools
import os

import organ_mass as _om

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PATH = os.path.join(_ROOT, "params", "neutral_composition_sourced.csv")
ORGANS = ("root", "stem", "leaf", "grain")
MID = ("mid", "mid")                       # (water corner, lipid corner)
CORNERS = tuple(itertools.product(("low", "high"), repeat=2))


def _load():
    out = {}
    with open(PATH) as f:
        for r in csv.DictReader(l for l in f if not l.startswith("#")):
            out[(r["plant"], r["organ"], r["quantity"])] = (float(r["low"]), float(r["high"]),
                                                             r["source"])
    return out


TABLE = _load()


def value(plant, organ, quantity, corner="mid"):
    lo, hi, _ = TABLE[(plant, organ, quantity)]
    return {"low": lo, "high": hi, "mid": 0.5 * (lo + hi)}[corner]


def rice_lipid_dw(corner="mid"):
    return {o: value("rice", o, "lipid_dw", corner) for o in ORGANS}


RICE_LIPID_DW = rice_lipid_dw()


def fresh_lipid(lipid_dw, waters):
    """L_fw = L_dw * (1 - W), organ by organ."""
    return {o: lipid_dw[o] * (1.0 - waters[o]) for o in lipid_dw}


# Plant labels in the tables -> composition rows (PREREG 5.10). Li 2019's hydroponic species
# without a value in Li's list take the cereal range (rice) or the 13-crop range (the rest).
LI_SPECIES = {"barley": "barley", "hordeum vulgare": "barley", "cucurbita pepo": "pumpkin",
              "rice": "cereal", "rice seedling": "cereal"}
TABLE_PLANT = {"liu2023": {"root": "rice_hydroponic_seedling"},
               "ge2017": {"root": "rice_soil_vegetative", "leaf": "rice_soil_vegetative"},
               "briggs1983": {"stem": "barley_seedling", "leaf": "barley_seedling"}}


def spec(table, species=None, corner=MID, f_lip=None):
    """(waters, fresh lipids) {organ: value} for one table row.

    table   : "rice" (the tool), "liu2023", "ge2017", "kodesova2019", "li2019_rcf",
              "li2019_soil", "briggs1983".
    species : the row's plant, for the per-species tables.
    corner  : (water corner, lipid corner), each "mid" | "low" | "high".
    f_lip   : Li 2019 soil's per-row fresh lipid (at Li's 0.90 water).
    """
    wc, lc = corner
    W = dict(_om.THETA_FW)
    Ldw = rice_lipid_dw(lc)
    if table in TABLE_PLANT:
        for organ, plant in TABLE_PLANT[table].items():
            if (plant, organ, "water") in TABLE:
                W[organ] = value(plant, organ, "water", wc)
            if (plant, organ, "lipid_dw") in TABLE:
                Ldw[organ] = value(plant, organ, "lipid_dw", lc)
    elif table == "kodesova2019":
        W["root"] = value(species, "root", "water", wc)
        Ldw["root"] = value(species, "root", "lipid_dw", lc)
    elif table in ("li2019_rcf", "li2019_soil"):
        W["root"] = value("li2019", "root", "water", wc)
        if table == "li2019_soil":
            Ldw["root"] = float(f_lip) / (1.0 - value("li2019", "root", "water"))
        else:
            plant = LI_SPECIES.get(str(species).strip().lower(), "crop")
            Ldw["root"] = value(plant, "root", "lipid_dw", lc)
    elif table != "rice":
        raise ValueError(f"unknown neutral table {table!r}")
    return W, fresh_lipid(Ldw, W)
