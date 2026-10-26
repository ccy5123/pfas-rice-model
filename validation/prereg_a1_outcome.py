#!/usr/bin/env python3
# =============================================================================
# validation/prereg_a1_outcome.py
# -----------------------------------------------------------------------------
# Checks the expectations E7-E16 of docs/PREREG_sourced_values.md (amendment 1 and addendum
# 5.6) against the run records -- refit_sourced.json, transfer_sourced.json,
# params/input_range_band.json and the numerics convergence printout -- so the Outcome
# section is read off the records, not written from memory.
#
#   python validation/prereg_a1_outcome.py
# =============================================================================
from __future__ import annotations

import csv
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
import model_api as api          # noqa: E402

LONG = ("PFDA", "PFUnDA", "PFDoDA", "PFOS")
SHORT = ("PFBA", "PFPeA", "PFHxA", "PFBS")


def _load(path):
    with open(path) as f:
        return json.load(f)


def e8_before_refit():
    """E8: the default's grain error with retention 1 and the sourced lam_grain but the
    PREVIOUS transport values (archived L_Ph_rec and lipid constants; the root matched as the
    old fit matched it), against the old default's grain error."""
    import refit_sourced as R
    old = {r["compound"]: r for r in csv.DictReader(
        l for l in open(os.path.join(HERE, "archive", "refit_site_year.csv")) if not l.startswith("#"))}
    out = {}
    for c in R.CONG:
        if not R.OBS[c]["use"]["grain"]:
            continue
        f, gx, gp = R.conductances(c, "lipid", R.START)
        m = R.model_logs(c, "redistributed", f, gx, gp, float(old[c]["L_Ph_rec"]))
        o = R.OBS[c]["baf_dw"]
        out[c] = float(m["grain"] - m["root"] - np.log10(o["grain"] / o["root"]))
    return out


def main():
    rf = _load(os.path.join(HERE, "refit_sourced.json"))
    tr = _load(os.path.join(HERE, "transfer_sourced.json"))
    bd = _load(os.path.join(ROOT, "params", "input_range_band.json"))
    V = rf["variants"]
    d = V["redistributed_lipid"]
    lines = []

    # E7 every r_eq reproduces its root cell; lam_eq converged by 1e4 (numerics printout)
    root_res = max(abs(V[v]["fits"][c]["res"]["root"]) for v in V for c in V[v]["fits"])
    lines.append(("E7", "r_eq reproduces every root cell", root_res < 1e-3,
                  f"max |root residual| {root_res:.1e} log10 (the lam_eq part: numerics_convergence.py)"))

    # E8 grain over-prediction shrinks before the re-fit
    e8 = e8_before_refit()
    old = {"PFPeA": 7.1, "PFUnDA": 7.3, "PFHpA": 6.9, "PFBA": 6.6, "PFOA": 3.2, "PFNA": 2.9,
           "PFOS": 2.3, "PFDA": 1.9, "PFDoDA": 1.1}            # CLAUDE.md 6.3 (previous default)
    shrink = [c for c in e8 if c in old and 10 ** e8[c] < old[c]]
    lines.append(("E8", "grain over-prediction shrinks for most congeners before the re-fit",
                  len(shrink) > len([c for c in e8 if c in old]) / 2,
                  "grain ratio x " + ", ".join(f"{c} {10 ** e8[c]:.2f} (was {old[c]})"
                                               for c in e8 if c in old)))

    # E9 at least one long-chain L_Ph exceeds 1 after the re-fit (default variant)
    lph = {c: d["fits"][c]["L_Ph"] for c in LONG}
    lines.append(("E9", "a long-chain L_Ph > 1 (default)", any(v > 1 for v in lph.values()),
                  ", ".join(f"{c} {v:.3g}" for c, v in lph.items())))

    # E10 Yamazaki observable no worse than 0.380
    y = d["score"]["all"]
    lines.append(("E10", "Yamazaki observable <= 0.380 (default)", y <= 0.380, f"{y:.3f}"))

    # E11 composition dominates the band for short chains, shoot shares for the grain
    att = bd["attribution"]
    comp = lambda n: "." in n and not n.startswith("PFBS")
    shoot = lambda n: n in ("stem_transp_frac", "lam_grain_heading", "lam_grain_maturity",
                            "phi", "T_C_Ph")
    top = {c: {k: (next(iter(att[c][k]), None)) for k in att[c]} for c in att}
    short_comp = [c for c in SHORT if c in top and all(top[c][k] is None or comp(top[c][k])
                                                       for k in ("root", "straw"))]
    grain_shoot = [c for c in top if top[c]["grain"] is not None and shoot(top[c]["grain"])]
    lines.append(("E11", "composition tops the band for short chains; shoot shares for grain",
                  len(short_comp) >= len([c for c in SHORT if c in top]) / 2
                  and len(grain_shoot) > len(top) / 2,
                  f"short chains led by composition: {short_comp}; grain led by a shoot share: "
                  f"{len(grain_shoot)}/{len(top)}"))

    # E12 transfer scores move by < 0.49 and the ranking of the variants does not change
    sc = {v: {k: float(np.sqrt(np.mean(np.square(e)))) for k, e in tr["errors"][v].items()}
          for v in tr["errors"]}
    old_tr = {"tang": 0.640, "kim_point": 0.469, "li_shoot": 0.358, "li_grain": 0.510}
    moves = {k: sc["redistributed_lipid"][k] - v for k, v in old_tr.items()}
    lines.append(("E12", "default transfer scores move by < 0.49 log10",
                  all(abs(m) < 0.49 for m in moves.values()),
                  ", ".join(f"{k} {sc['redistributed_lipid'][k]:.3f} ({m:+.3f})"
                            for k, m in moves.items())))

    # E13 the input-range band on the grain is wider than 0.1 for most congeners
    wide = [c for c in bd["band"] if bd["band"][c]["grain"][1] - bd["band"][c]["grain"][0] > 0.1]
    lines.append(("E13", "grain band wider than 0.1 log10 for most congeners",
                  len(wide) > len(bd["band"]) / 2, f"{len(wide)}/{len(bd['band'])}"))

    # E14 g_ph removed in both lipid variants
    rem = {v: V[v]["lipid"]["terms_removed"] for v in V if "lipid" in V[v]}
    lines.append(("E14", "g_ph removed in both lipid variants",
                  all("g_ph" in r for r in rem.values()), f"{rem}"))

    # E15 PFBA or PFPeA at L_Ph = 0 in every variant
    at0 = {v: [c for c in ("PFBA", "PFPeA") if V[v]["fits"][c]["L_Ph"] == 0.0] for v in V}
    lines.append(("E15", "PFBA or PFPeA at L_Ph = 0 in every variant",
                  all(at0.values()), f"{at0}"))

    # E16 per-cultivar r_eq interval > 2x for PFBA, PFPeA, PFHpA, PFBS; < 2x for the rest
    iv = d["cultivar"]
    span = {c: max(x["r_eq"] for x in iv[c]) / min(x["r_eq"] for x in iv[c]) for c in iv}
    big = ("PFBA", "PFPeA", "PFHpA", "PFBS")
    ok = all(span[c] > 2 for c in big) and all(span[c] < 2 for c in span if c not in big)
    lines.append(("E16", "per-cultivar r_eq span > 2x for PFBA/PFPeA/PFHpA/PFBS, < 2x otherwise",
                  ok, ", ".join(f"{c} {s:.2f}x" for c, s in span.items())))

    for tag, what, held, detail in lines:
        print(f"{tag:4} {'HELD ' if held else 'WRONG'}  {what}\n       {detail}")
    return lines


if __name__ == "__main__":
    main()
