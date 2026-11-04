#!/usr/bin/env python3
"""Snapshot of the PFAS default outputs, to check that a change moved them only where registered.

Cells: for every congener in parameters.json, `simulate(c)["baf_final_dw"]` (each organ / pool) and the Yamazaki
two-run observable `yamazaki_harvest_baf(c)` (root, straw, grain). Writes JSON; `--compare old.json new.json`
prints the cell count and the largest |log10 change| per congener.
Usage: python validation/default_snapshot.py out.json  |  python validation/default_snapshot.py --compare a.json b.json
"""
import json
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))


def snapshot():
    import model_api as api
    out = {}
    for c in api.CONGENERS:
        r = api.simulate(c)
        for k, v in r["baf_final_dw"].items():
            out[f"{c}|simulate|{k}"] = float(v)
        if c in api.yamazaki_obs():
            for k, v in api.yamazaki_harvest_baf(c).items():
                out[f"{c}|yamazaki|{k}"] = float(v)
    return out


def compare(a, b):
    A, B = (json.load(open(p)) for p in (a, b))
    keys = sorted(set(A) & set(B))
    d = {}
    for k in keys:
        x, y = A[k], B[k]
        if x > 0 and y > 0 and math.isfinite(x) and math.isfinite(y):
            d[k] = math.log10(y / x)
    print(f"{len(keys)} cells; max |dlog10| {max(abs(v) for v in d.values()):.4f}")
    per = {}
    for k, v in d.items():
        c = k.split("|")[0]
        per.setdefault(c, []).append((abs(v), k, v))
    for c, vs in per.items():
        m = max(vs)
        print(f"  {c:7s} max {m[2]:+.4f} ({m[1].split('|', 1)[1]})")


if __name__ == "__main__":
    if sys.argv[1] == "--compare":
        compare(sys.argv[2], sys.argv[3])
    else:
        json.dump(snapshot(), open(sys.argv[1], "w"), indent=1)
