#!/usr/bin/env python3
# =============================================================================
# validation/theta_fw_sensitivity.py
# -----------------------------------------------------------------------------
# HOW MUCH DOES THE ORGAN WATER FRACTION theta_fw MATTER?  (handoff step 2)
#
# PRE-REGISTRATION -- written before any number below was produced.
#
# Why: theta_fw is held constant over the season and its values (root 0.90, stem 0.83,
# leaf 0.78, grain 0.14) cite source keys that resolve nowhere in the repo. Real rice
# organs dry during the season (filling grain ~0.5-0.6 -> 0.14 at harvest; mature straw
# 0.1-0.5). Making theta a function of time would make B_k time-dependent -- an ODE
# change. This script measures whether that change is worth making.
#
# Where theta enters, and how it is varied here: (i) B_k = theta + (1-theta)*X, (ii) the
# crop model's dry matter -> fresh M = M_dry/(1-theta) (src/organ_mass.py), (iii) every
# dw<->fw comparison. All three read the SAME value, so one organ's theta is changed
# everywhere at once: model_api._COMP, organ_mass.THETA_FW and neutral_dpu.RICE_WATER.
# Ranges: params/rice_tissue_params.csv low/high (root .88-.92, stem .80-.88,
# leaf .72-.82, grain .14-.60), one organ at a time.
#
# Metrics are read where the DATA are dry weight, so an observation cannot move with the
# assumed theta: Tang 2026 per-organ TF (lipid loading, 0.1 ug/g) and Ge 2017 neutral
# per-organ TF. Model BAFs are also printed on both bases for three congeners.
#
# DECISION RULE: if every dry-side metric moves < 0.05 log10 RMSE across the ranges,
# constant theta stays and theta(t) is not built. If any moves > 0.10, theta(t) is
# worth building for that organ. In between: record it and decide with the user.
#
# EXPECTATION (analytic, stated so it can be wrong): an organ's burden is
# Cw * M_dry * (theta/(1-theta) + X). Where X (binding per kg dry) is large, theta
# nearly cancels on a dry basis. The ROOT's water term theta/(1-theta) ~ 9 dominates for
# short chains, and the root is every TF's denominator, so root theta should matter most.
#
#   python validation/theta_fw_sensitivity.py
# =============================================================================
from __future__ import annotations

import contextlib
import io
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, HERE)

import model_api as api            # noqa: E402
import organ_mass as om            # noqa: E402
import neutral_dpu as ND           # noqa: E402

# organ -> (model_api composition key, low, high) from params/rice_tissue_params.csv
RANGES = {"root": ("root", 0.88, 0.92), "stem": ("stem", 0.80, 0.88),
          "leaf": ("leaf", 0.72, 0.82), "grain": ("grain_brown", 0.14, 0.60)}
PFAS = ("PFOA", "PFOS", "PFDoDA")
GE2017 = os.path.join(ROOT, "data_obs", "neutral_obs_ge2017.csv")
MOVE_SMALL, MOVE_LARGE = 0.05, 0.10


@contextlib.contextmanager
def theta(organ, value):
    """Set one organ's theta_fw everywhere it is read, then restore it."""
    key = RANGES[organ][0]
    saved = (api._COMP[key]["theta_fw"], om.THETA_FW[organ], ND.RICE_WATER[organ])
    api._COMP[key]["theta_fw"] = om.THETA_FW[organ] = ND.RICE_WATER[organ] = float(value)
    try:
        yield
    finally:
        api._COMP[key]["theta_fw"], om.THETA_FW[organ], ND.RICE_WATER[organ] = saved


def pfas_bafs(name):
    """Root/straw/grain BAF on the fresh basis and on the dry basis (harvest operator)."""
    r = api.simulate(name)
    C = {k: r["conc"][k][-1] / r["cwo_ref"] for k in api.TISSUES}
    fw = dict(root=C["root"], straw=r["straw_baf"], grain=C["grain"])
    dw = {k: r["baf_final_dw"][k] for k in ("root", "straw", "grain")}
    return fw, dw


def tang_lipid_rmse():
    e = []
    for c in api.TANG_CONGENERS:
        v = api.tang_tf_validation(c, dose="low", lipid_loading=True)
        e += [np.log10(v["model_tf"][o]) - np.log10(v["tang_tf"][o]) for o in v["organs"]]
    return float(np.sqrt(np.mean(np.square(e))))


def ge2017_rmse():
    import neutral_dpu_validation as NV
    with contextlib.redirect_stdout(io.StringIO()):
        return float(NV.compare_to_obs(GE2017, quiet=True))


def evaluate():
    out = dict(tang=tang_lipid_rmse(), ge=ge2017_rmse())
    for nm in PFAS:
        out[nm] = pfas_bafs(nm)
    return out


def main():
    print("=" * 84)
    print("theta_fw SENSITIVITY — one organ at a time over params/rice_tissue_params.csv ranges")
    print("=" * 84)
    base = evaluate()
    rows = [("baseline", None, None, base)]
    for organ, (_, lo, hi) in RANGES.items():
        for v in (lo, hi):
            with theta(organ, v):
                rows.append((organ, v, None, evaluate()))

    print(f"\n{'scenario':16}{'Tang lipid':>11}{'Ge 2017':>9}   "
          + "   ".join(f"{nm} root/straw/grain (dw BAF)" for nm in PFAS[:1]))
    for organ, v, _, m in rows:
        tag = "baseline" if organ is None or organ == "baseline" else f"{organ} {v:.2f}"
        d = m["PFOA"][1]
        print(f"{tag:16}{m['tang']:>11.3f}{m['ge']:>9.3f}   "
              f"{d['root']:.3f} / {d['straw']:.3f} / {d['grain']:.3f}")

    print("\nmove in each dry-side metric vs baseline (max |Δ| over an organ's low/high):")
    verdict = {}
    for organ in RANGES:
        ms = [m for o, _, _, m in rows if o == organ]
        dt = max(abs(m["tang"] - base["tang"]) for m in ms)
        dg = max(abs(m["ge"] - base["ge"]) for m in ms)
        dbaf = max(abs(np.log10(m[nm][1][t]) - np.log10(base[nm][1][t]))
                   for m in ms for nm in PFAS for t in ("root", "straw", "grain"))
        fbaf = max(abs(np.log10(m[nm][0][t]) - np.log10(base[nm][0][t]))
                   for m in ms for nm in PFAS for t in ("root", "straw", "grain"))
        worst = max(dt, dg)
        verdict[organ] = ("SMALL" if worst < MOVE_SMALL else
                          "LARGE" if worst > MOVE_LARGE else "INTERMEDIATE")
        print(f"  {organ:6} Tang {dt:.3f}  Ge {dg:.3f}  |  PFAS BAF log10: dw {dbaf:.3f}, "
              f"fw {fbaf:.3f}   -> {verdict[organ]}")

    print("\nper-congener BAF (fresh | dry) at baseline:")
    for nm in PFAS:
        fw, dw = base[nm]
        print(f"  {nm:7} fw root/straw/grain {fw['root']:.3f}/{fw['straw']:.3f}/{fw['grain']:.3f}"
              f"   dw {dw['root']:.3f}/{dw['straw']:.3f}/{dw['grain']:.3f}")

    print("\n" + "=" * 84)
    print("VERDICT (against the pre-registered rule)")
    print("=" * 84)
    for organ, v in verdict.items():
        print(f"  {organ:6} {v}")
    print("  SMALL = constant theta stays for that organ; LARGE = theta(t) worth building;")
    print("  INTERMEDIATE = record and decide with the user. Model math UNCHANGED either way.")
    xylem_split_note()
    return rows, verdict


def xylem_split_note():
    """POST-HOC (not pre-registered): the xylem split and grain theta.

    The first run of this script found that grain theta moved the 4-pool grain BAF
    through the leaf/grain xylem split, which then ran on FRESH mass (S*M), so the
    organs' water set their transpiration shares. Decided with the user: the split runs
    on DRY matter with S_XYLEM_DRY (leaf 20, grain 0.5116 per kg dry; see the core), so
    grain theta no longer moves it. This prints the check."""
    from pfas_rice_plant_module_4pool_surf import xylem_leaf_share
    print("\nPOST-HOC — the xylem split runs on DRY matter (S_XYLEM_DRY), so grain theta no "
          "longer moves it")
    for th in (RANGES["grain"][1], RANGES["grain"][2]):
        with theta("grain", th):
            r = api.simulate("PFOA")
            share = 1.0 - xylem_leaf_share(api._compartments(), r["M"][-1], "fresh")
            print(f"  grain theta {th:.2f}: grain share of the xylem at harvest {share:.3f}; "
                  f"PFOA grain BAF (dw) {r['conc']['grain'][-1] / (1 - th) / r['cwo_ref']:.3f}")
    print("  (Reading the original grain S=2 as per kg DRY instead would give 0.343 at harvest,")
    print("  ~14% during filling: docs/HANDOFF_mass_basis_refit.md section 2.)")


if __name__ == "__main__":
    main()
