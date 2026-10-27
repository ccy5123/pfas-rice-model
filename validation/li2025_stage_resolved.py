#!/usr/bin/env python3
# =============================================================================
# validation/li2025_stage_resolved.py
# -----------------------------------------------------------------------------
# LI 2025 STAGE BY STAGE  (handoff step 3, docs/HANDOFF_mass_basis_refit.md)
#
# PRE-REGISTRATION -- written and pushed before any model number below was produced.
#
# Why: the repo's Li comparison (validation/archive/oos_crossdataset.li_tf) sets the model's
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
# Model: the dynamic run on the default drivers (the calibration site-year; Li's own site
# drivers are not in the repo), Cwo constant 1; dry-weight concentrations from the harvest
# operator (src/harvest.py: leaf = green + attached dead blades at that time). Stage ->
# model day by the ORYZA2000 DVS of that run, from sources (docs/PREREG_sourced_values.md
# 5.9; until then 0.4 / 1.0 / 1.25 / 1.5 / 2.0 were chosen): heading DVS 1, ripening DVS 2,
# tillering the midpoint of [DVS at transplanting, 0.65], filling and milk-ripening the centre
# of their ordered range between heading and ripening (4/3, 5/3). Variants as in the archived
# oos_crossdataset (4-pool: monotone f_xy,
# saturated "fit", lipid loading) plus nstem_leaf with lipid loading, because Li reports stem and
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
#   SENSITIVITY (PREREG 5.9): the range ends -- tillering at both ends of its range;
#      (filling, milk) at the corners (1, 1), (1, 2), (2, 2) of their ordered range.
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

DATA = os.path.join(ROOT, "docs", "literature_db", "raw_si", "li2025_stage_conc.csv")
STAGES = ("tillering", "heading", "filling", "milk", "ripening")
PI_DVS = 0.65                                  # panicle initiation (ORYZA2000 definition)


def tillering_range():
    """[DVS at transplanting on the default run, panicle initiation]: tillering lies between
    them (PREREG 5.9)."""
    return float(np.asarray(api.site_run()["dvs"], float)[0]), PI_DVS


def stage_dvs(tillering=None, filling=1.0 + 1.0 / 3.0, milk=1.0 + 2.0 / 3.0):
    """The stage DVS (PREREG 5.9): heading 1 and ripening 2 by definition; tillering the
    midpoint of `tillering_range()`; filling and milk-ripening the centroid of
    1 < filling < milk < 2."""
    if tillering is None:
        tillering = 0.5 * sum(tillering_range())
    return (float(tillering), 1.0, float(filling), float(milk), 2.0)
CONGENERS = ("PFBA", "PFHxA", "PFOA", "PFBS", "PFOS")
SHOOT, GRAIN = ("stem", "leaf"), "grain"
GRAIN_FROM = STAGES.index("filling")          # grain scored from filling on (pre-registered)
VARIANTS = ("mono", "fit", "lipid", "nstem_lipid")
# The original rule's fixed numbers (noise ratio 3, margin 0.10) are retired
# (docs/PREREG_sourced_values.md section 2): Q2 is a chi-square test of the stage scatter
# against the reported SDs, the Q3 verdict a paired bootstrap (decision_stats.decide); both at
# the conventional 0.05 level. N_PERM is numerical (the permutation p's resolution).
ALPHA, N_PERM = 0.05, 5000


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
def stage_days(dvs_targets=None):
    """Model day of each stage DVS on the default drivers' ORYZA2000 run (`stage_dvs()` by
    default)."""
    dvs_targets = stage_dvs() if dvs_targets is None else dvs_targets
    run = api.site_run()
    t, dvs = np.asarray(run["t"], float), np.asarray(run["dvs"], float)
    return np.array([t[np.argmax(dvs >= d - 1e-9)] if d >= dvs.max() - 1e-9
                     else float(np.interp(d, dvs, t)) for d in dvs_targets])


VARIANT_KW = dict(mono=dict(f_xy_source="recommended", **api.FOUR_POOL_FREE),
                  fit=dict(f_xy_source="fit", lipid_loading=False),
                  lipid=dict(lipid_loading=True, shoot="4pool"),
                  nstem_lipid=dict(lipid_loading=True, shoot="redistributed"))


def run_variant(c, v):
    """One variant with its own fitted values (validation/refit_sourced.py)."""
    return api.simulate(c, Cwo=1.0, **VARIANT_KW[v])


def _dw_at(res, days):
    t = np.asarray(res["t"], float)
    return {k: np.interp(days, t, res["conc_dw"][k]) for k in ("root", "stem", "leaf", "grain")}


def model_tf(res, days, grain="grain"):
    """Per-stage model TF (dry weight, harvest operator) at the given days. `grain` = "grain"
    (brown rice) or "rough" (brown + husk + rachis): Li's "paddy rice / grains" is undefined
    (PREREG 7.2), so it is scored on both."""
    dw = _dw_at(res, days)
    if grain == "rough":
        dw["grain"] = np.interp(days, res["t"], res["conc_dw"]["rough"])
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
    err = {g: list(map(float, e)) for g, e in lev.items()}
    o = np.concatenate([s[0] for s in series]); m = np.concatenate([s[1] for s in series])
    rho = _spearman(m, o)
    null = np.empty(N_PERM)
    for j in range(N_PERM):
        op = np.concatenate([rng.permutation(s[0]) for s in series])
        null[j] = _spearman(m, op)
    p = float((np.sum(null >= rho) + 1) / (N_PERM + 1))
    return dict(rmse=rmse, err=err, rho=rho, p=p, n_series=len(series), n_points=len(o))


def verdict(scores):
    """One variant that decide() calls "better" than every other on BOTH groups, or None."""
    import decision_stats as D
    for v in VARIANTS:
        others = [u for u in VARIANTS if u != v]
        if all(D.decide(np.array(scores[u]["err"][g]), np.array(scores[v]["err"][g]))["verdict"]
               == "better" for u in others for g in ("shoot", "grain")):
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
    ratios, chi2, dof = [], 0.0, 0
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
            lt = np.log10(tf[idx])
            chi2 += float(np.sum(((lt - lt.mean()) / sig[idx]) ** 2))
            dof += len(idx) - 1
            print(f"  {c + ' ' + k:14} across-stage SD {across:.2f} log  vs SD-implied {within:.2f}"
                  f"  -> {across / within:5.1f}x")
    from scipy.stats import chi2 as _chi2
    p_noise = float(_chi2.sf(chi2, dof))
    noise_limited = bool(p_noise < ALPHA)
    print(f"  median {np.median(ratios):.1f}x; chi-square {chi2:.0f} on {dof} df, p {p_noise:.2g} -> "
          f"{'NOISE-LIMITED: per-stage shape results are not used for a verdict' if noise_limited else 'the SDs describe the scatter'}")

    days = stage_days()
    print("\nstage -> model day (ORYZA DVS, PREREG 5.9): "
          + ", ".join(f"{s} {d:.0f} (DVS {v:.2f})" for s, d, v in zip(STAGES, days, stage_dvs())))
    run = api.site_run()
    t_r, dvs_r = np.asarray(run["t"], float), np.asarray(run["dvs"], float)
    t_head = float(np.interp(1.0, dvs_r, t_r))
    chk = [float(np.interp(t_head + dd, t_r, dvs_r)) for dd in (14.0, 20.0)]
    print(f"  check (not used): filling 14-20 d after heading start in Nanjing (Li R. et al. 2018)"
          f" is DVS {chk[0]:.2f}-{chk[1]:.2f} on this run; the centroid uses {stage_dvs()[2]:.2f}")
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
        lo, hi = tillering_range()
        ends = {f"tillering {lo:.2f}": stage_dvs(tillering=lo),
                f"tillering {hi:.2f}": stage_dvs(tillering=hi)}
        ends.update({f"filling/milk {f:g}/{m:g}": stage_dvs(filling=f, milk=m)
                     for f, m in ((1.0, 1.0), (1.0, 2.0), (2.0, 2.0))})
        print("\nSENSITIVITY — the stage ranges' ends (PREREG 5.9):")
        for label, dv in ends.items():
            _, sc = evaluate(stage_days(dv))
            sens[label] = sc
            print(f"  {label:22}: " + "  ".join(
                f"{v} {sc[v]['rmse']['shoot']:.2f}/{sc[v]['rmse']['grain']:.2f} "
                f"(rho {sc[v]['rho']:+.2f}, p {sc[v]['p']:.2f})" for v in VARIANTS)
                + f"   -> verdict {verdict(sc) or 'none'}")

    post = posthoc(data, runs, days)

    print("\n" + "=" * 88)
    print("VERDICT (against the pre-registered rule)")
    print("=" * 88)
    print(f"  Q2 noise: {'noise-limited' if noise_limited else 'not noise-limited'} "
          f"(median {np.median(ratios):.1f}x)")
    print(f"  Q3 discriminating variant (bootstrap 'better' on shoot AND grain): {win or 'none'}"
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
    runs["nstem_mono"] = {c: api.simulate(c, Cwo=1.0, shoot="redistributed", lipid_loading=False)
                          for c in CONGENERS}
    out = {}
    print("\nPOST-HOC (not pre-registered; not used for the verdict)")
    print("  (a) Q3 with the other estimators, applied identically to model and data:")
    for est in ("mean_of_ratios", "ratio_of_means"):
        row = {}
        for v in VARIANTS:
            model = {}
            for c in CONGENERS:
                model[c] = {**model_tf(runs[v][c], days), "_dw": _dw_at(runs[v][c], days)}
            row[v] = level_rmse(data, model, est)
        out[est] = row
        print(f"    {est:15}" + "  ".join(f"{v} {row[v]['shoot']:.2f}/{row[v]['grain']:.2f}"
                                          for v in VARIANTS))
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
