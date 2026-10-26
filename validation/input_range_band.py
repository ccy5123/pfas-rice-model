#!/usr/bin/env python3
# =============================================================================
# validation/input_range_band.py
# -----------------------------------------------------------------------------
# THE INPUT-RANGE BAND of the default configuration (docs/PREREG_sourced_values.md 5.4).
#
# Every input known only as a sourced range runs at its midpoint. Here each one is moved to
# each end of its range, ONE AT A TIME, with the fitted parameters held (r_eq, L_Ph, lipid
# constants of the default variant), on the Yamazaki two-harvest observable (dry basis,
# root / straw / grain, 11 congeners). An end that moves a reported number by >= 0.0005 log10
# counts; per congener and tissue the log deviations are summed per direction (down: the
# negative deviations, up: the positive ones) -- a conservative envelope. It is checked
# against the corner runs: per congener, tissue and direction, every input at the end that
# moved that number that way, all at once.
#
# Inputs: the shoot and phloem shares (model_api.RANGE_INPUTS: stem_transp_frac,
# lam_grain at heading and at maturity, phi, T_C_Ph), the site inputs (sbdur, tmpsb, nplh;
# the crop is re-run), every organ-composition cell with a range (f_prot, f_PL, f_cw, lignin
# share -> K_cw, the grain's harvest water; parameters.json tissue_composition_ranges), the
# organ water's stage anchors one anchor at a time (PREREG 7.1; root, stem and leaf water no
# longer enter as constants), the husk share at maturity, the rachis and the husk's phloem draw
# (7.2), the published japonica crop sets (7.5; a categorical input: every set is an "end", the
# band takes the extreme deviations), and PFBS's L_Ph range [0, L_MQL] (its grain is censored).
#
#   python validation/input_range_band.py     (writes params/input_range_band.json)
# =============================================================================
from __future__ import annotations

import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
import model_api as api                 # noqa: E402
import organ_mass as om                 # noqa: E402
import panicle as pn                    # noqa: E402
import crop_sets as CS                  # noqa: E402

PRECISION = 0.0005
TISS = ("root", "straw", "grain")
OBS = api.yamazaki_obs()
CONG = [c for c in api.CONGENERS if c in OBS]
ORGAN_OF = {"root": "root", "stem": "stem", "leaf": "leaf", "grain_brown": "grain"}
OUT = os.path.join(ROOT, "params", "input_range_band.json")


def inputs():
    """[(name, (low setting, high setting))]; a setting is a dict of overrides."""
    R = api.RANGE_INPUTS
    out = []
    for k in ("stem_transp_frac", "phi", "T_C_Ph"):
        out.append((k, tuple({"sim": {k: v}} for v in R[k])))
    h, m = api.LAM_GRAIN_DVS
    out.append(("lam_grain_heading", tuple({"sim": {"lam_grain": (v, m)}}
                                           for v in R["lam_grain_heading"])))
    out.append(("lam_grain_maturity", tuple({"sim": {"lam_grain": (h, v)}}
                                            for v in R["lam_grain_maturity"])))
    for k in ("sbdur", "tmpsb", "nplh"):
        out.append((k, tuple({"site": {k: v}} for v in R[k])))
    rng = api.PARAMS["tissue_composition_ranges"]
    for organ in ORGAN_OF:
        for q, (lo, hi) in rng[organ].items():
            if q == "theta_fw" and ORGAN_OF[organ] in om.THETA_STAGE:
                continue                        # staged: its anchors below (PREREG 7.1)
            if hi > lo:
                out.append((f"{organ}.{q}", tuple({"comp": (organ, q, v)} for v in (lo, hi))))
    for organ, a in om.THETA_STAGE.items():
        for i, d in enumerate(a["dvs"]):
            if a["high"][i] > a["low"][i]:
                out.append((f"{organ}.theta@DVS{d:g}",
                            tuple({"anchor": (organ, i, v)} for v in (a["low"][i], a["high"][i]))))
    out.append(("husk_share_maturity", tuple({"panicle": {"h_mat": v}} for v in pn.H_MAT_RANGE)))
    out.append(("rachis", tuple({"panicle": {"rachis_frac": v}} for v in pn.RACHIS_FRAC_RANGE)))
    out.append(("husk_phloem", ({"sim": {"husk_phloem": False}}, {"sim": {"husk_phloem": True}})))
    out.append(("crop_set", tuple({"site": {"crop_set": cs}} for cs in CS.BAND_SETS)))
    out.append(("PFBS.L_Ph", ({"pfbs_L": 0.0}, {"pfbs_L": "mql"})))
    return out


def _apply_comp(organ, q, v):
    """Patch one composition cell (and, for the lignin share, every congener's organ K_cw; for
    theta_fw also the dry -> fresh conversion of the organ mass, which reads the same theta)."""
    saved = (dict(api._COMP[organ]), {c: dict(api._CONG[c]["K_cw_wholecw_Lkg"]) for c in api._CONG},
             dict(om.THETA_FW))
    if q == "theta_fw":
        om.THETA_FW[ORGAN_OF[organ]] = v
    if q == "lignin_share":
        o = ORGAN_OF[organ]
        for c, rec in api._CONG.items():
            kp, kl = rec.get("K_cw_poly_rec_Lkg"), rec.get("K_lignin_rec_Lkg")
            if kp is not None and kl is not None:
                rec["K_cw_wholecw_Lkg"][o] = (1.0 - v) * kp + v * kl
    else:
        api._COMP[organ][q] = v
    return saved


def _restore(organ, saved):
    api._COMP[organ].clear()
    api._COMP[organ].update(saved[0])
    om.THETA_FW.clear()
    om.THETA_FW.update(saved[2])
    for c, kc in saved[1].items():
        api._CONG[c]["K_cw_wholecw_Lkg"].clear()
        api._CONG[c]["K_cw_wholecw_Lkg"].update(kc)


def observable(job):
    """log10 {root, straw, grain} of the two-run observable for one congener and settings."""
    c, settings = job
    sim, site, patched, anchors, pan = {}, {}, [], {}, {}
    for st in settings:
        sim.update(st.get("sim", {}))
        site.update(st.get("site", {}))
        pan.update(st.get("panicle", {}))
        if "anchor" in st:                      # merged per organ, so a corner moves several
            organ, i, v = st["anchor"]
            a = om.THETA_STAGE[organ]
            vals = anchors.setdefault(organ, list(0.5 * (a["low"] + a["high"])))
            vals[i] = v
        if "comp" in st:
            organ, q, v = st["comp"]
            patched.append((organ, _apply_comp(organ, q, v)))
        if "pfbs_L" in st and c == "PFBS":
            v = st["pfbs_L"]
            sim["L_Ph_override"] = (api._CONG["PFBS"]["L_Ph_MQL__" + api.variant_key()]
                                    if v == "mql" else v)
    if anchors:
        site["theta_anchors"] = anchors
    if pan:
        site["panicle_kw"] = pan
    try:
        b = api.yamazaki_harvest_baf(c, scenario_kw=site or None, **sim)
    finally:
        for organ, saved in reversed(patched):
            _restore(organ, saved)
    return c, {k: float(np.log10(b[k])) for k in TISS}


def main():
    ins = inputs()
    with ProcessPoolExecutor(max_workers=os.cpu_count() or 4) as pool:
        base = dict(pool.map(observable, [(c, ()) for c in CONG]))
        jobs = [(c, (st,)) for _, ends in ins for st in ends for c in CONG]
        res = list(pool.map(observable, jobs))
        it = iter(res)
        delta = {}                                      # input -> [lo end, hi end] -> c -> k
        for name, ends in ins:
            delta[name] = []
            for _ in ends:
                d = {}
                for _c in CONG:
                    c, v = next(it)
                    d[c] = {k: v[k] - base[c][k] for k in TISS}
                delta[name].append(d)
        band, attrib = {}, {}
        for c in CONG:
            band[c], attrib[c] = {}, {}
            for k in TISS:
                dn = up = 0.0
                contrib = {}
                for name in delta:
                    ds = [d_[c][k] for d_ in delta[name]]
                    if max(abs(x) for x in ds) < PRECISION:
                        continue
                    dn += min(0.0, *ds)
                    up += max(0.0, *ds)
                    contrib[name] = max(abs(x) for x in ds)
                band[c][k] = [dn, up]
                attrib[c][k] = dict(sorted(contrib.items(), key=lambda z: -z[1])[:5])
        # corner check: every counted input at the end that moved this number this way
        corner_jobs, keys = [], []
        for c in CONG:
            for k in TISS:
                for sgn in (-1, 1):
                    sts = []
                    for name, ends in ins:
                        ds = [d_[c][k] for d_ in delta[name]]
                        if max(abs(x) for x in ds) < PRECISION:
                            continue
                        i = int(np.argmin(ds)) if sgn < 0 else int(np.argmax(ds))
                        if sgn * ds[i] > 0:
                            sts.append(ends[i])
                    corner_jobs.append((c, tuple(sts)))
                    keys.append((c, k, sgn))
        corner = {}
        for (c, k, sgn), (_, v) in zip(keys, pool.map(observable, corner_jobs)):
            corner.setdefault(c, {}).setdefault(k, [0.0, 0.0])[0 if sgn < 0 else 1] = \
                v[k] - base[c][k]
    print(f"{'congener':8} " + "  ".join(f"{k:>22}" for k in TISS))
    for c in CONG:
        print(f"{c:8} " + "  ".join(
            f"[{band[c][k][0]:+.3f}, {band[c][k][1]:+.3f}] ({corner[c][k][0]:+.3f},{corner[c][k][1]:+.3f})"
            for k in TISS))
    cover = [(abs(corner[c][k][i]) <= abs(band[c][k][i]) + PRECISION)
             for c in CONG for k in TISS for i in (0, 1)]
    print(f"\nenvelope covers the corner run in {sum(cover)}/{len(cover)} cases")
    for c in CONG:
        print(f"  {c:7} largest contributors: " + "; ".join(
            f"{k}: " + ", ".join(f"{n} {v:.3f}" for n, v in list(attrib[c][k].items())[:3])
            for k in TISS))
    with open(OUT, "w") as f:
        json.dump(dict(source="validation/input_range_band.py (PREREG 5.4)",
                       variant=api.variant_key(), units="log10 [down, up] around the default",
                       band=band, corner=corner, attribution=attrib,
                       inputs=[n for n, _ in ins], base=base), f, indent=1)
    print(f"[written] {os.path.relpath(OUT, ROOT)}")


if __name__ == "__main__":
    main()
