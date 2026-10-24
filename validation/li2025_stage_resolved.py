#!/usr/bin/env python3
# =============================================================================
# validation/li2025_stage_resolved.py
# -----------------------------------------------------------------------------
# LI 2025 STAGE BY STAGE  (handoff step 3, docs/HANDOFF_mass_basis_refit.md)
#
# PRE-REGISTRATION -- written and pushed before any model number below was produced.
#
# Why: the repo's Li comparison (validation/oos_crossdataset.li_tf) sets the model's
# HARVEST-day BAF next to the data's 5-STAGE MEAN concentrations, converted to fresh
# weight with assumed theta. Two mismatches: time (harvest vs season mean) and basis
# (the data are dry weight). The authors' own TF is a third estimator (the arithmetic
# mean of per-stage ratios), 1.0-5.8x above the repo's.
#
# Data: SI Table S6 -> docs/literature_db/raw_si/li2025_stage_conc.csv (ng/g dw; 5
# stages without dates x root/stem/leaf/grain, grain from heading). Transcription
# checked before this was written: the 5-stage means reproduce data_obs/Li2025_BAF_TF.csv
# in all 20 congener x tissue cells, and the mean of per-stage ratios reproduces the
# authors' PFOA stem TF 2.35.
#
# Observable: TF_k = C_k / C_root, dry weight, per stage. It needs no water value (the
# irrigation-water denominator is Li's documented weak point); the root-surface confound
# stays (roots washed with distilled water, no plaque removal).
#
# Model: the dynamic run on the default ORYZA2000 biomass + measured Q_TP, Cwo constant 1;
# model fresh-weight concentration -> dry weight with the binding factor's theta. Stage ->
# model day by ORYZA DVS: tillering 0.4, heading 1.0, filling 1.25, milk 1.5, ripening
# 2.0 (first day it is reached). Variants as in oos_crossdataset (4-pool: monotone f_xy,
# W2 fit, lipid loading) plus nstem_leaf with lipid loading, because Li reports stem and
# leaf separately and the 4-pool stem is a pass-through. Congeners PFBA, PFHxA, PFOA,
# PFBS, PFOS (PFHpA is ND in every tissue). Grain is scored from FILLING on: at heading
# the model's grain-formation gate has barely opened (structural; the value is printed).
#
# Questions and decision rules:
#   Q1 (estimator) How far apart are the three TF estimators on the DATA: ratio of means
#      (repo), arithmetic mean of ratios (authors), geometric mean of ratios? Reported.
#      Scoring uses the geometric mean, the estimator consistent with a log10 metric.
#   Q2 (noise) Across-stage SD of the observed log10 TF vs the SD the reported +/-SDs
#      imply. If the median ratio exceeds 3, the reported SDs do not describe the
#      stage-to-stage variation and per-stage SHAPE results are treated as noise-limited.
#   Q3 (level) log10 RMSE of the stage-geometric-mean TF -- the SAME estimator over the
#      SAME stages on both sides -- per variant, for shoot (stem + leaf) and grain.
#   Q4 (shape) Pooled Spearman between de-meaned model and observed log10 TF trajectories
#      per variant; permutation p (5000 shuffles within each series). A shape signal
#      counts only at p < 0.05.
#   VERDICT RULE: Li moves from "inconclusive" to discriminating only if one variant beats
#      every other by more than 0.10 log10 RMSE on Q3 for BOTH shoot AND grain. Otherwise
#      Li stays inconclusive, and what changes is only how the comparison is built
#      (matched times, dry basis, a stated estimator).
#   EXPECTATION (stated so it can be wrong): Li stays inconclusive; the estimator choice
#      alone moves the observed TF by more than the model variants differ; Q4 finds no
#      shape signal (p > 0.05).
#   SENSITIVITY (fixed in advance): every stage DVS shifted by -0.15 and +0.15 (ripening
#      held at 2.0).
#
#   python validation/li2025_stage_resolved.py
# =============================================================================
from __future__ import annotations

import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))

import model_api as api            # noqa: E402
import oryza_growth as og          # noqa: E402

DATA = os.path.join(ROOT, "docs", "literature_db", "raw_si", "li2025_stage_conc.csv")
STAGES = ("tillering", "heading", "filling", "milk", "ripening")
STAGE_DVS = (0.4, 1.0, 1.25, 1.5, 2.0)
CONGENERS = ("PFBA", "PFHxA", "PFOA", "PFBS", "PFOS")
SHOOT, GRAIN = ("stem", "leaf"), "grain"
GRAIN_FROM = STAGES.index("filling")          # grain scored from filling on (pre-registered)
VARIANTS = ("mono", "W2", "lipid", "nstem_lipid")
SEASON, N_T = 120.0, 241
THETA = {"root": api._COMP["root"]["theta_fw"], "stem": api._COMP["stem"]["theta_fw"],
         "leaf": api._COMP["leaf"]["theta_fw"], "grain": api._COMP["grain_brown"]["theta_fw"]}
NOISE_LIMIT, MARGIN, ALPHA, N_PERM = 3.0, 0.10, 0.05, 5000


# ---------------------------------------------------------------------------- data
def load_data(path=DATA):
    """{(compound, tissue): (mean[5], sd[5])} in stage order; NaN where ND / not sampled."""
    out = {}
    with open(path) as f:
        for r in csv.DictReader(l for l in f if not l.startswith("#")):
            key = (r["compound"], r["tissue"])
            m, s = out.setdefault(key, (np.full(5, np.nan), np.full(5, np.nan)))
            i = int(r["stage_order"]) - 1
            if r["nd"] == "0":
                m[i], s[i] = float(r["mean_ng_g_dw"]), float(r["sd_ng_g_dw"])
    return out


def scored_stages(tissue):
    return range(GRAIN_FROM, 5) if tissue == GRAIN else range(5)


def observed_tf(data, c, k):
    """Per-stage TF_k = C_k / C_root (dw) and its log10 SD implied by the reported SDs."""
    mk, sk = data[(c, k)]
    mr, sr = data[(c, "root")]
    tf = mk / mr
    sig = 0.4343 * np.sqrt((sk / mk) ** 2 + (sr / mr) ** 2)
    return tf, sig


def estimators(data, c, k):
    """Q1: the three TF estimators on the data (repo, authors, geometric)."""
    mk, _ = data[(c, k)]
    mr, _ = data[(c, "root")]
    tf, _ = observed_tf(data, c, k)
    ok = ~np.isnan(tf)
    if not ok.any():
        return None
    return dict(ratio_of_means=float(np.nanmean(mk) / np.nanmean(mr)),
                mean_of_ratios=float(np.mean(tf[ok])),
                geo_mean_ratio=float(np.exp(np.mean(np.log(tf[ok])))))


# ---------------------------------------------------------------------------- model
def stage_days(dvs_targets=STAGE_DVS, season=SEASON):
    sim = og.simulate_oryza(og.OryzaParams(season=season))
    t, dvs = np.asarray(sim["t"], float), np.asarray(sim["dvs"], float)
    return np.array([t[np.argmax(dvs >= d)] if d >= dvs.max() else float(np.interp(d, dvs, t))
                     for d in dvs_targets])


def run_variant(c, v):
    if v == "nstem_lipid":
        return api.simulate_nstem_leaf(c, Cwo=1.0, season=SEASON, n_t=N_T, lipid_loading=True,
                                       biomass_fn=api._biomass_fn("oryza"))
    kw = dict(mono=dict(f_xy_source="recommended"), W2=dict(f_xy_source="W2fit"),
              lipid=dict(lipid_loading=True))[v]
    return api.simulate(c, Cwo=1.0, season=SEASON, n_t=N_T, **kw)


def model_tf(res, days):
    """Per-stage model TF (dry weight) at the given days, for stem/leaf/grain."""
    t = np.asarray(res["t"], float)
    dw = {k: np.interp(days, t, res["conc"][k]) / (1.0 - THETA[k]) for k in THETA}
    return {k: dw[k] / dw["root"] for k in ("stem", "leaf", "grain")}


# ---------------------------------------------------------------------------- scoring
def _ranks(x):
    r = np.empty(len(x)); r[np.argsort(x)] = np.arange(len(x)); return r


def _spearman(a, b):
    return float(np.corrcoef(_ranks(a), _ranks(b))[0, 1])


def score(data, model, rng):
    """Q3 level errors (shoot, grain) and Q4 pooled shape correlation for one variant."""
    lev = {"shoot": [], "grain": []}
    series = []
    for c in CONGENERS:
        for k in SHOOT + (GRAIN,):
            if (c, k) not in data:
                continue
            tf_o, _ = observed_tf(data, c, k)
            idx = [i for i in scored_stages(k) if not np.isnan(tf_o[i])]
            if len(idx) < 2:
                continue
            lo, lm = np.log10(tf_o[idx]), np.log10(np.maximum(model[c][k][idx], 1e-12))
            lev["grain" if k == GRAIN else "shoot"].append(lm.mean() - lo.mean())
            series.append((lo - lo.mean(), lm - lm.mean()))
    rmse = {g: float(np.sqrt(np.mean(np.square(e)))) for g, e in lev.items()}
    o = np.concatenate([s[0] for s in series]); m = np.concatenate([s[1] for s in series])
    rho = _spearman(m, o)
    null = np.empty(N_PERM)
    for j in range(N_PERM):
        op = np.concatenate([rng.permutation(s[0]) for s in series])
        null[j] = _spearman(m, op)
    p = float((np.sum(null >= rho) + 1) / (N_PERM + 1))
    return dict(rmse=rmse, rho=rho, p=p, n_series=len(series), n_points=len(o))


def verdict(scores):
    """Pre-registered rule: one variant beats every other by > MARGIN on BOTH groups."""
    for v in VARIANTS:
        others = [u for u in VARIANTS if u != v]
        if all(scores[u]["rmse"][g] - scores[v]["rmse"][g] > MARGIN
               for u in others for g in ("shoot", "grain")):
            return v
    return None


# ---------------------------------------------------------------------------- report
def main(sensitivity=True):
    rng = np.random.default_rng(0)
    data = load_data()
    print("=" * 88)
    print("LI 2025 STAGE BY STAGE — pre-registered (see header); TF = tissue/root, dry weight")
    print("=" * 88)

    print("\nQ1 — three TF estimators on the DATA (stem/leaf over 5 stages, grain over its 4):")
    print(f"  {'series':14}{'ratio of means':>16}{'mean of ratios':>16}{'geo. mean':>11}{'spread':>8}")
    spreads = []
    for c in CONGENERS:
        for k in SHOOT + (GRAIN,):
            e = estimators(data, c, k) if (c, k) in data else None
            if e is None:
                continue
            vals = list(e.values()); spreads.append(max(vals) / min(vals))
            print(f"  {c + ' ' + k:14}{e['ratio_of_means']:>16.3f}{e['mean_of_ratios']:>16.3f}"
                  f"{e['geo_mean_ratio']:>11.3f}{spreads[-1]:>7.1f}x")
    print(f"  median spread {np.median(spreads):.1f}x, max {max(spreads):.1f}x")

    print("\nQ2 — is the stage-to-stage scatter explained by the reported SDs?")
    ratios = []
    for c in CONGENERS:
        for k in SHOOT + (GRAIN,):
            if (c, k) not in data:
                continue
            tf, sig = observed_tf(data, c, k)
            idx = [i for i in scored_stages(k) if not np.isnan(tf[i])]
            if len(idx) < 3:
                continue
            across = float(np.std(np.log10(tf[idx]), ddof=1))
            within = float(np.sqrt(np.mean(sig[idx] ** 2)))
            ratios.append(across / within)
            print(f"  {c + ' ' + k:14} across-stage SD {across:.2f} log  vs SD-implied {within:.2f}"
                  f"  -> {across / within:5.1f}x")
    noise_limited = bool(np.median(ratios) > NOISE_LIMIT)
    print(f"  median {np.median(ratios):.1f}x -> "
          f"{'NOISE-LIMITED: per-stage shape results are not used for a verdict' if noise_limited else 'the SDs describe the scatter'}")

    days = stage_days()
    print("\nstage -> model day (ORYZA DVS): "
          + ", ".join(f"{s} {d:.0f}" for s, d in zip(STAGES, days)))
    runs = {v: {c: run_variant(c, v) for c in CONGENERS} for v in VARIANTS}

    def evaluate(stage_day):
        model = {v: {c: model_tf(runs[v][c], stage_day) for c in CONGENERS} for v in VARIANTS}
        return model, {v: score(data, model[v], rng) for v in VARIANTS}

    model, scores = evaluate(days)
    print("\nQ3/Q4 — same estimator (stage geometric mean) over the same stages on both sides:")
    print(f"  {'variant':12}{'shoot RMSE':>11}{'grain RMSE':>11}{'shape rho':>11}{'perm p':>8}")
    for v in VARIANTS:
        s = scores[v]
        print(f"  {v:12}{s['rmse']['shoot']:>11.2f}{s['rmse']['grain']:>11.2f}"
              f"{s['rho']:>11.2f}{s['p']:>8.3f}")
    print("\n  per-series stage-geometric-mean TF (obs | " + " | ".join(VARIANTS) + "):")
    for c in CONGENERS:
        for k in SHOOT + (GRAIN,):
            if (c, k) not in data:
                continue
            tf_o, _ = observed_tf(data, c, k)
            idx = [i for i in scored_stages(k) if not np.isnan(tf_o[i])]
            if len(idx) < 2:
                continue
            gm = lambda x: float(np.exp(np.mean(np.log(np.maximum(x[idx], 1e-12)))))
            print(f"  {c + ' ' + k:14}{gm(tf_o):>8.3f} | "
                  + " | ".join(f"{gm(model[v][c][k]):.3f}" for v in VARIANTS))
    print("\n  grain at HEADING (not scored; the model's grain-formation gate has barely opened):")
    for c in CONGENERS:
        tf_o, _ = observed_tf(data, c, GRAIN) if (c, GRAIN) in data else (np.full(5, np.nan), None)
        print(f"  {c:6} obs {tf_o[1]:.3f}  model " + "  ".join(
            f"{v} {model[v][c]['grain'][1]:.2g}" for v in VARIANTS))

    win = verdict(scores)
    sens = {}
    if sensitivity:
        print("\nSENSITIVITY — stage DVS shifted (ripening held at 2.0):")
        for shift in (-0.15, +0.15):
            dv = tuple(d + shift if d < 2.0 else d for d in STAGE_DVS)
            _, sc = evaluate(stage_days(dv))
            sens[shift] = sc
            print(f"  shift {shift:+.2f}: " + "  ".join(
                f"{v} {sc[v]['rmse']['shoot']:.2f}/{sc[v]['rmse']['grain']:.2f} "
                f"(rho {sc[v]['rho']:+.2f}, p {sc[v]['p']:.2f})" for v in VARIANTS)
                + f"   -> verdict {verdict(sc) or 'none'}")

    post = posthoc(data, runs, days)

    print("\n" + "=" * 88)
    print("VERDICT (against the pre-registered rule)")
    print("=" * 88)
    print(f"  Q2 noise: {'noise-limited' if noise_limited else 'not noise-limited'} "
          f"(median {np.median(ratios):.1f}x)")
    print(f"  Q3 discriminating variant (> {MARGIN} on shoot AND grain): {win or 'none'}"
          f" -> Li {'discriminates' if win else 'stays INCONCLUSIVE'}")
    sig = [v for v in VARIANTS if scores[v]["p"] < ALPHA]
    print(f"  Q4 shape signal at p < {ALPHA}: {', '.join(sig) if sig else 'none'}")
    return dict(spread_median=float(np.median(spreads)), spread_max=float(max(spreads)),
                noise_ratio=float(np.median(ratios)), noise_limited=noise_limited,
                days=days, scores=scores, winner=win, sensitivity=sens, posthoc=post)


def level_rmse(data, model, estimator):
    """Q3 with another TF estimator applied identically to model and data (post-hoc)."""
    lev = {"shoot": [], "grain": []}
    for c in CONGENERS:
        for k in SHOOT + (GRAIN,):
            if (c, k) not in data:
                continue
            tf_o, _ = observed_tf(data, c, k)
            idx = [i for i in scored_stages(k) if not np.isnan(tf_o[i])]
            if len(idx) < 2:
                continue
            if estimator == "mean_of_ratios":
                o, m = np.mean(tf_o[idx]), np.mean(model[c][k][idx])
            else:                                   # ratio of means over the same stages
                mk, _ = data[(c, k)]
                mr, _ = data[(c, "root")]
                o = np.mean(mk[idx]) / np.mean(mr[idx])
                m = np.mean(model[c]["_dw"][k][idx]) / np.mean(model[c]["_dw"]["root"][idx])
            lev["grain" if k == GRAIN else "shoot"].append(np.log10(max(m, 1e-12)) - np.log10(o))
    return {g: float(np.sqrt(np.mean(np.square(e)))) for g, e in lev.items()}


def posthoc(data, runs, days):
    """POST-HOC (not pre-registered), reported beside the verdict, never used for it:
    (a) does the Q3 ranking survive the other two estimators, applied to both sides?
    (b) is the winner's advantage the redistributed SHOOT or the lipid term? -- the
        same nstem_leaf shoot with the free anion only (no lipid loading)."""
    runs = dict(runs)
    runs["nstem_mono"] = {c: api.simulate_nstem_leaf(c, Cwo=1.0, season=SEASON, n_t=N_T,
                                                     biomass_fn=api._biomass_fn("oryza"))
                          for c in CONGENERS}
    out = {}
    print("\nPOST-HOC (not pre-registered; not used for the verdict)")
    print("  (a) Q3 with the other estimators, applied identically to model and data:")
    for est in ("mean_of_ratios", "ratio_of_means"):
        row = {}
        for v in VARIANTS:
            model = {}
            for c in CONGENERS:
                t = np.asarray(runs[v][c]["t"], float)
                dw = {k: np.interp(days, t, runs[v][c]["conc"][k]) / (1.0 - THETA[k]) for k in THETA}
                model[c] = {**model_tf(runs[v][c], days), "_dw": dw}
            row[v] = level_rmse(data, model, est)
        out[est] = row
        print(f"    {est:15}" + "  ".join(f"{v} {row[v]['shoot']:.2f}/{row[v]['grain']:.2f}"
                                          for v in VARIANTS)
              + f"   -> {verdict({v: {'rmse': row[v]} for v in VARIANTS}) or 'none'}")
    rng = np.random.default_rng(0)
    sc = {v: score(data, {c: model_tf(runs[v][c], days) for c in CONGENERS}, rng)["rmse"]
          for v in ("nstem_mono", "nstem_lipid", "lipid")}
    out["decomposition"] = sc
    floor = {"shoot": [], "grain": []}
    for c in CONGENERS:
        for k in SHOOT + (GRAIN,):
            if (c, k) not in data:
                continue
            tf_o, _ = observed_tf(data, c, k)
            idx = [i for i in scored_stages(k) if not np.isnan(tf_o[i])]
            if len(idx) >= 2:
                se = np.std(np.log10(tf_o[idx]), ddof=1) / np.sqrt(len(idx))
                floor["grain" if k == GRAIN else "shoot"].append(se ** 2)
    out["noise_floor"] = {g: float(np.sqrt(np.mean(v))) for g, v in floor.items()}
    print("  noise floor of Q3 -- RMSE a perfect model would score, from the across-stage scatter"
          f" alone: shoot {out['noise_floor']['shoot']:.2f} / grain {out['noise_floor']['grain']:.2f}")
    print("  (b) shoot structure vs lipid term (shoot/grain RMSE, geometric mean):")
    for v, lab in (("lipid", "4-pool + lipid"), ("nstem_mono", "nstem_leaf, free anion only"),
                   ("nstem_lipid", "nstem_leaf + lipid")):
        print(f"    {lab:30} {sc[v]['shoot']:.2f} / {sc[v]['grain']:.2f}")
    return out


if __name__ == "__main__":
    main()
