"""
panicle.py -- the crop model's panicle (ORYZA WSO) split into brown rice and the rest
=====================================================================================

ORYZA's storage organ WSO is the whole panicle ("Storage organs (panicles)", Bouman et al. 2001
p.169), while the grain the model is scored on is brown rice (PREREG 7.2). The panicle is split:

    brown rice  = (WSO - R) * (1 - h(t))
    remainder   = WSO - brown rice            (husk + rachis + the panicle built before flowering)

h(t), the husk share of rough grain (dry), is 1 at DVS 1 (the caryopsis starts at fertilisation)
and then follows the Nipponbare per-grain hull and brown-rice weights of Ansari et al. 2017 Fig. 1
(`docs/literature_db/raw_si/ansari2017_fig1_digitised.csv`, checked against the figure), scaled so
that its maturity value equals the maturity anchor `h_mat` (PREREG 7.10), linear in DVS between the
anchors. Days after heading map to DVS on the run's own grid (heading = DVS 1: the site-year runs
are calibrated to reach DVS 1 on the measured heading date); an anchor that falls at or beyond
DVS 2 (a harvest within 30 d of heading) is dropped in favour of the maturity value.

R, the rachis, has no source (PREREG 7.2): ORYZA sends all panicle growth after DVS 0.95 into
grains, so the rachis is complete by then and cannot exceed the panicle at DVS 0.95. R is run
over [0, WSO(DVS 0.95)] as `rachis_frac` in [0, 1] of it (midpoint 0.5); before DVS 0.95 it is
that fraction of the growing panicle.
"""
from __future__ import annotations

import csv
import os

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_ANSARI = os.path.join(_ROOT, "docs", "literature_db", "raw_si", "ansari2017_fig1_digitised.csv")

# husk share of rough grain at maturity: Yoshida 1981 Table 1.1 (11 varieties, 0.186-0.254) and
# Ansari 2017 Nipponbare (0.161); the run uses the midpoint, the band the ends (PREREG 7.2)
H_MAT_RANGE = (0.161, 0.254)
RACHIS_FRAC_RANGE = (0.0, 1.0)          # R / WSO(DVS 0.95): no source, the physical range
DVS_RACHIS = 0.95                       # ORYZA: GGR = GSO after DVS 0.95 (Bouman et al. 2001 p.57)


def _mid(r):
    return 0.5 * (r[0] + r[1])


H_MAT = _mid(H_MAT_RANGE)
RACHIS_FRAC = _mid(RACHIS_FRAC_RANGE)


def _load_ansari():
    with open(_ANSARI) as f:
        rows = [r for r in csv.DictReader(line for line in f if not line.startswith("#"))]
    v = {(r["panel"], r["time"]): float(r["value_mg"]) for r in rows}
    out = {}
    for time in ("10DAH", "20DAH", "30DAH", "maturity"):
        hull, brown = v[("A", time)], v[("B", time)]
        out[time] = hull / (hull + brown)
    return out


H_NIPPONBARE = _load_ansari()           # {10DAH, 20DAH, 30DAH, maturity: husk share}
DAH = {"10DAH": 10.0, "20DAH": 20.0, "30DAH": 30.0}


def husk_share(t, dvs, h_mat=H_MAT):
    """Husk share of rough grain on the grid `t` (days) with development stage `dvs`."""
    t = np.asarray(t, float)
    dvs = np.asarray(dvs, float)
    if dvs.max() < 1.0:                              # no heading in the run: no caryopsis
        return np.ones_like(t)
    t_head = float(np.interp(1.0, dvs, t))
    scale = h_mat / H_NIPPONBARE["maturity"]
    pts = [(1.0, 1.0)]
    for key, days in DAH.items():
        d = float(np.interp(t_head + days, t, dvs))
        if t_head + days <= t[-1] and d < 2.0:
            pts.append((d, H_NIPPONBARE[key] * scale))
    pts.append((2.0, float(h_mat)))
    x, y = zip(*sorted(pts))
    return np.interp(dvs, x, y)


def split(t, wso, dvs, h_mat=H_MAT, rachis_frac=RACHIS_FRAC):
    """(brown rice, remainder) DRY mass on the grid from the panicle dry mass `wso`."""
    t, wso, dvs = (np.asarray(a, float) for a in (t, wso, dvs))
    w95 = float(np.interp(DVS_RACHIS, dvs, wso)) if dvs.max() >= DVS_RACHIS else None
    R = rachis_frac * (np.minimum(wso, w95) if w95 is not None else wso)
    R = np.where(dvs >= DVS_RACHIS, rachis_frac * w95 if w95 is not None else R, R)
    brown = np.maximum(wso - R, 0.0) * (1.0 - husk_share(t, dvs, h_mat))
    return brown, wso - brown


def split_biomass(b, t, h_mat=H_MAT, rachis_frac=RACHIS_FRAC):
    """A dry biomass dict whose "grain" is brown rice, with the panicle in "panicle" and the
    remainder in "husk" (all kg DM/hill). `b` must carry "dvs"; without it the dict is returned
    unchanged (no stage, no split)."""
    if b.get("dvs") is None:
        return dict(b)
    out = dict(b)
    brown, rest = split(t, b["grain"], b["dvs"], h_mat, rachis_frac)
    out.update(panicle=np.asarray(b["grain"], float), grain=brown, husk=rest)
    return out
