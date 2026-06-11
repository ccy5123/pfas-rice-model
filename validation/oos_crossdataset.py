#!/usr/bin/env python3
"""
Out-of-sample cross-dataset test: does the K_PL-gated lipid mechanism (fit to
Yamazaki 2023) PREDICT independent data without refitting?
=========================================================================

The transport parameters are fit to Yamazaki 2023 only. This script transfers
three model variants -- monotone `f_xy_recommended`, the saturated per-congener
W2 fit, and this session's K_PL-gated lipid loading (`lipid_loading=True`) -- to
two datasets they were NOT fit on and compares the predictions.

Targets (and why each is/ isn't usable):
  (1) Kim 2019 brown-rice (grain) BAF, porewater basis -- the cleaner test.
      PFOA was used to fit L_Ph, so it is flagged and excluded from the score;
      the long chains have low rice detection frequency (DF) and are unreliable.
  (2) Li 2025 (Tianjin field) -- only the WATER-INDEPENDENT tissue ratios (TF)
      are usable: the reported BAFs scale inversely with water quality (PFOS
      "poor 0.3-1.6%" -> BAF 250; PFOA "good" -> BAF ~2), i.e. the group-water
      denominator is unreliable. Even the TF carries a root surface-sorption
      confound, so Li is expected to be inconclusive.

Headline result: the lipid mechanism predicts the Kim grain chain-length pattern
(including the long-chain RISE the monotone/W2 models structurally miss) far
better out-of-sample -- the project's first genuine predictive signal. Li does
not discriminate. Run: python validation/oos_crossdataset.py
"""
import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
import model_api as api          # noqa: E402
import literature_params as lp   # noqa: E402

MODES = ("lipid", "mono", "W2")


def predict(nm, mode):
    if mode == "lipid":
        r = api.simulate(nm, lipid_loading=True, measured_forcing=True)
    elif mode == "mono":
        r = api.simulate(nm, f_xy_source="recommended", measured_forcing=True)
    else:
        r = api.simulate(nm, f_xy_source="W2fit", measured_forcing=True)
    return {"root": r["baf_final"]["root"], "straw": r["straw_baf"], "grain": r["baf_final"]["grain"]}


def kim_grain():
    """(1) Kim 2019 grain BAF — out-of-sample (PFOA* used in the L_Ph fit)."""
    kim = lp.kim2019_grain_baf("porewater")
    DF = {"PFHpA": 13, "PFOA": 57, "PFNA": 20, "PFDA": 6.7, "PFUnDA": 13, "PFDoDA": 3.3}
    print("=== (1) Kim 2019 grain BAF [porewater]  (OUT-OF-SAMPLE; * = used in fit) ===")
    print(f"{'PFAS':7}{'DF%':>5}{'obs':>9}{'lipid':>9}{'mono':>9}{'W2':>9}")
    rows = []
    for nm, o in kim.items():
        if nm not in api.CONGENERS:
            continue
        p = {m: predict(nm, m)["grain"] for m in MODES}
        print(f"{nm:6}{'*' if nm=='PFOA' else ' '}{DF.get(nm,0):>5.0f}{o:>9.2f}"
              f"{p['lipid']:>9.2f}{p['mono']:>9.2f}{p['W2']:>9.2f}")
        rows.append((nm, o, p, DF.get(nm, 0)))

    def rmse(mode, hiDF=False):
        e = [(np.log10(max(p[mode], 1e-6)) - np.log10(o)) ** 2
             for nm, o, p, df in rows if nm != "PFOA" and (not hiDF or df >= 15)]
        return float(np.sqrt(np.mean(e))) if e else float("nan")
    print("log10 RMSE (excl PFOA)      : " + "  ".join(f"{m}={rmse(m):.2f}" for m in MODES))
    print("log10 RMSE (DF>=15%, no PFOA): " + "  ".join(f"{m}={rmse(m, True):.2f}" for m in MODES)
          + "   [reliable: PFHpA, PFNA]")


def li_tf():
    """(2) Li 2025 water-independent tissue ratios (TF) — the only usable Li signal."""
    obs = {}
    with open(os.path.join(ROOT, "data_obs", "obs_baf_Li2025.csv")) as f:
        for r in csv.DictReader(f):
            obs.setdefault(r["compound"], {})[r["tissue"]] = float(r["baf"])
    print("\n=== (2) Li 2025 TF = tissue/tissue (WATER-INDEPENDENT) ===")
    print(f"{'PFAS':7}{'straw/root o':>13}{'lipid':>7}{'mono':>7}{'W2':>7}"
          f"   {'grain/root o':>13}{'lipid':>7}{'mono':>7}{'W2':>7}")
    es = {m: [] for m in MODES}
    eg = {m: [] for m in MODES}
    for nm in ("PFBA", "PFHxA", "PFOA", "PFBS", "PFOS"):
        o = obs[nm]
        P = {m: predict(nm, m) for m in MODES}
        ts_o = o["straw"] / o["root"]
        tg_o = o.get("grain", np.nan) / o["root"] if "grain" in o else np.nan
        ts = {m: P[m]["straw"] / P[m]["root"] for m in MODES}
        tg = {m: P[m]["grain"] / P[m]["root"] for m in MODES}
        for m in MODES:
            es[m].append((np.log10(ts[m]) - np.log10(ts_o)) ** 2)
            if not np.isnan(tg_o):
                eg[m].append((np.log10(max(tg[m], 1e-6)) - np.log10(tg_o)) ** 2)
        print(f"{nm:7}{ts_o:>13.2f}{ts['lipid']:>7.2f}{ts['mono']:>7.2f}{ts['W2']:>7.2f}"
              f"   {tg_o:>13.2f}{tg['lipid']:>7.2f}{tg['mono']:>7.2f}{tg['W2']:>7.2f}")
    print("log10 RMSE TF straw/root: " + "  ".join(f"{m}={np.sqrt(np.mean(es[m])):.2f}" for m in MODES))
    print("log10 RMSE TF grain/root: " + "  ".join(f"{m}={np.sqrt(np.mean(eg[m])):.2f}" for m in MODES))
    print("(Li is field/group-water/surface-confounded -> expected to be inconclusive.)")


if __name__ == "__main__":
    kim_grain()
    li_tf()
