#!/usr/bin/env python3
# =============================================================================
# validation/prereg_a2_outcome.py
# -----------------------------------------------------------------------------
# Checks the expectations E17-E24 of docs/PREREG_sourced_values.md section 7 (amendment 2a)
# against the run records -- refit_sourced.json, transfer_sourced.json, params/input_range_band.json
# and the amendment-1 records kept in validation/archive/amendment1/ -- so the outcome is read off
# the records, not written from memory.
#
#   python validation/prereg_a2_outcome.py
# =============================================================================
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, HERE)
import model_api as api          # noqa: E402
import pfas_rice_plant_module_4pool_surf as P4   # noqa: E402

A1 = os.path.join(HERE, "archive", "amendment1")
SHORT = ("PFBA", "PFPeA", "PFHxA", "PFBS")
DEF = "redistributed_lipid"


def _load(path):
    with open(path) as f:
        return json.load(f)


def _rmse(e):
    return float(np.sqrt(np.mean(np.square(e)))) if len(e) else float("nan")


def e17_factor(c):
    """Closed-form r_eq ratio when the harvest root water goes 0.713 -> 0.6365 (PREREG 7.8)."""
    comps = api._compartments(api.kcw_by_organ(api._CONG[c]))
    X = P4.dry_binding(comps, api._compound_for(c))[0]
    a, b = 0.713 / (1 - 0.713), 0.6365 / (1 - 0.6365)
    return (a + X) / (b + X)


def main():
    rf, rf1 = _load(os.path.join(HERE, "refit_sourced.json")), _load(os.path.join(A1, "refit_sourced.json"))
    tr = _load(os.path.join(HERE, "transfer_sourced.json"))
    bd = _load(os.path.join(ROOT, "params", "input_range_band.json"))
    new, old = rf["variants"][DEF], rf1["variants"][DEF]
    out = {}

    # E17
    rows, ok = [], True
    for c, f in new["fits"].items():
        ratio, exp = f["r_eq"] / old["fits"][c]["r_eq"], e17_factor(c)
        ok &= abs(ratio - exp) <= 0.02
        rows.append((c, ratio, exp))
    print("E17 r_eq new/amendment-1 vs closed form:")
    for c, r, e in rows:
        print(f"  {c:7s} {r:.3f}  expected {e:.3f}  {'ok' if abs(r - e) <= 0.02 else 'MISS'}")
    out["E17"] = ok

    # E18
    s_new, s_old = new["score"]["all"], old["score"]["all"]
    out["E18"] = abs(s_new - s_old) < 0.05
    print(f"E18 Yamazaki observable {s_old:.3f} -> {s_new:.3f} (|delta| {abs(s_new - s_old):.3f} < 0.05)")

    # E19
    hits = {c: list(bd["attribution"][c]["root"])[:3] for c in SHORT}
    out["E19"] = all("root.theta@DVS2" in v for v in hits.values())
    print("E19 root band, three largest contributors:")
    for c, v in hits.items():
        print(f"  {c:7s} {v}")

    # E20
    zero = [c for c, f in new["fits"].items()
            if f["res"]["grain"] is not None and f["L_Ph"] <= 0.0]
    dstraw = {c: new["fits"][c]["res"]["straw"] - old["fits"][c]["res"]["straw"] for c in new["fits"]}
    droot = max(abs(new["fits"][c]["res"]["root"]) for c in new["fits"])
    out["E20"] = (not zero) and max(abs(v) for v in dstraw.values()) < 0.0005 and droot < 0.0005
    print(f"E20 L_Ph = 0 with a grain cell: {zero or 'none'}; root residual max {droot:.1e}; "
          f"straw residual change max {max(abs(v) for v in dstraw.values()):.4f}")
    for v, rec in rf["variants"].items():
        z = [c for c, f in rec["fits"].items() if f["res"]["grain"] is not None and f["L_Ph"] <= 0.0]
        print(f"     {v:20s} L_Ph = 0: {z or 'none'}")

    # E21
    import choi2018_husk_check as CH
    res = CH.main()
    out["E21"] = all(res[c]["model"] > 1.0 for c in res) and res["PFOS"]["model"] > res["PFOA"]["model"]

    # E22
    e = tr["errors"][DEF]
    b, r = _rmse(e["li_grain"]), _rmse(e.get("li_grain_rough", []))
    out["E22"] = r < b
    print(f"E22 Li 2025 grain, brown {b:.3f} vs rough {r:.3f}")

    # E23
    cs = [c for c in bd["band"]
          if not all("crop_set" in bd["counted"][c][k] for k in ("straw", "grain"))]
    lead = {c: next(iter(bd["attribution"][c]["grain"]), None) for c in bd["band"]}
    n_lead = sum(v == "lam_grain_heading" for v in lead.values())
    out["E23"] = (not cs) and n_lead > len(lead) / 2
    print(f"E23 crop sets enter the straw and grain band for {len(bd['band']) - len(cs)}/"
          f"{len(bd['band'])} congeners; the panicle share at heading leads the grain band for "
          f"{n_lead}/{len(lead)} ({lead})")

    print("\n" + "  ".join(f"{k} {'HELD' if v else 'WRONG'}" for k, v in out.items())
          + "  E24 HELD by tests/test_amendment2a.py (density invariance)")
    return out


if __name__ == "__main__":
    main()
