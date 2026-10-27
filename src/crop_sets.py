"""
crop_sets.py -- published ORYZA crop parameter sets for temperate japonica, as overrides on IR72.D92
===================================================================================================

No ORYZA crop set has been published for Koshihikari or Kasalath (the Yamazaki cultivars;
docs/literature_db/sourced_inputs.csv crop_cultivar_search). The site-year therefore runs IR72.D92
verbatim (the only complete published set) and carries the published japonica sets as the range of
the crop parameters:

  gao2021        Gao et al. 2021 Agric. Water Manag. 250:106864, Table 3 (cold-region japonica,
                 Heilongjiang): FLVTB, FSTTB, FSOTB, DRLVT, SLATB (table), KDF, RGRLMX, FRPAR, DVR*;
                 COLDMIN is not applied (the port has no cold-stop rule).
  wu2022_<site>  Wu et al. 2022 Appl. Ecol. Environ. Res. 20:633, Tables 3-4 (five Liaoning sites,
                 ORYZA v3): DVR*, FSTTB, FLVTB, FSOTB, SLATB (table).

Everything a paper does not list stays IR72's, as in the papers' own crop files. The development
rates of every set are re-calibrated to the measured heading and harvest dates on the site-year
(`site_scenario.calibrate_rates`), so a set enters through its partitioning, leaf area and leaf death.
"""
from __future__ import annotations

import csv
import os

import numpy as np

_RAW = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                    "docs", "literature_db", "raw_si")
RATE_KEYS = ("DVRJ", "DVRI", "DVRP", "DVRR")
WU_SITES = ("huanren", "donggang", "liaozhong", "kaiyuan", "dawa")
# Kaiyuan's printed partitioning equals Donggang's (a possible copy error in the paper), so it is
# not used (PREREG 7.5); the band runs IR72 against the other five sets
WU_USED = tuple(s for s in WU_SITES if s != "kaiyuan")
CROP_SETS = ("IR72", "gao2021") + tuple(f"wu2022_{s}" for s in WU_USED)
BAND_SETS = CROP_SETS[1:]


def _table(v):
    pairs = [p.split(":") for p in v.split(";")]
    return np.array([float(a) for a, _ in pairs]), np.array([float(b) for _, b in pairs])


def _rows(name):
    with open(os.path.join(_RAW, name)) as f:
        return list(csv.DictReader(l for l in f if not l.startswith("#")))


def _gao():
    rows = {r["parameter"]: r["calibrated"] for r in _rows("gao2021_cold_japonica_oryza_params.csv")}
    out = {k: _table(rows[k]) for k in ("FLVTB", "FSTTB", "FSOTB", "DRLVT", "SLATB")}
    out["SWISLA"] = "TABLE"
    out["KDFTB"] = (np.array([0.0, 2.5]), np.array([float(rows["KDF"])] * 2))   # a single KDF
    for k in ("RGRLMX", "FRPAR") + RATE_KEYS:
        out[k] = float(rows[k])
    return out


def _wu(site):
    r = next(r for r in _rows("wu2022_liaoning_oryza_params.csv") if r["site"].lower() == site)
    out = {k: _table(r[k]) for k in ("FSTTB", "FLVTB", "FSOTB", "SLATB")}
    out["SWISLA"] = "TABLE"
    out.update({k: float(r[k]) for k in RATE_KEYS})
    return out


def overrides(name="IR72"):
    """Crop-dict overrides on IR72.D92 for a published set (empty for IR72 itself)."""
    if name == "IR72":
        return {}
    if name == "gao2021":
        return _gao()
    if name.startswith("wu2022_") and name[7:] in WU_USED:
        return _wu(name[7:])
    raise ValueError(f"unknown crop set {name!r}; expected one of {CROP_SETS}")
