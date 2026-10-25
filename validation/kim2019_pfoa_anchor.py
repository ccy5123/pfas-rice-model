#!/usr/bin/env python3
# =============================================================================
# validation/kim2019_pfoa_anchor.py
# -----------------------------------------------------------------------------
# KIM 2019 PFOA GRAIN BAF: WHICH ESTIMATOR?  (handoff step 4b, docs/HANDOFF_mass_basis_refit.md)
#
# PRE-REGISTRATION -- written and pushed before the censored estimate, the slope, the
# leave-one-out checks and the model-side numbers below were produced.
#
# Why: the repo's Kim PFOA grain BAF (4.43 L/kg) is a ratio of two field AVERAGES, and
# one WWTP-impacted field (A-2, 1760 ng/L) carries 82.8% of the void-water PFOA sum.
# Main-text Table 1 (p.718) gives PFOA per site, so the BAF can be estimated per site.
# Data: docs/literature_db/raw_si/kim2019_table1_sites.csv (transcribed by the user from
# the PDF text layer, checked row by row against the rendered page). Reproduces SI Table
# S6: rice PFOA mean 0.3494 ng/g with 17/30 detects (SI 0.349, 57%), void water 78.73 ng/L
# with 26/27 (SI 78.7, 96%), soil 0.1599 ng/g dw with 13/30 (SI 0.160, 43%).
#
# Already computed (data only, before this file):
#   E1  published ratio of means (rice over 30 sites / water over 27)   4.44 L/kg
#   E2  ratio of means over the 27 sites that have water (ND = 0)      4.18
#   E3  ratio of means without A-2                                     21.5
#   E4  per-site geometric mean, the 14 sites with BOTH detected       27.0 (median 28.3,
#       range 0.98 at A-2 to 597 at D-7)
# Also already seen: kim2019_censored.py (step 4a) PRINTED the four variants' PFOA grain
# BAF as unscored anchor rows, all far below 4.43. Nothing below uses them to choose an
# estimator. And: no shipped parameter is fitted to Kim PFOA (L_Ph_fit, L_Ph_rec and
# LIPID_LOADING are Yamazaki fits; only the literature_params DEMO fits L_Ph to 4.43), so
# the "PFOA* used in the L_Ph fit" exclusion in the OOS scripts has no basis.
#
# Site accounting: 30 sites. 14 have rice and water both detected; 12 have water detected
# and rice < MDL (0.088 ng/g) -> per-site BAF < 0.088/water, a LEFT-CENSORED ratio; C-6
# has both < MDL (no information on the ratio); E-2..E-4 have no water (NS).
#
# Questions and decision rules:
#   Q1 (data) Which estimator? Admissibility, stated before E5 is computed:
#        (a) numerator and denominator come from the same sites;
#        (b) no site is dropped for its detection status or its value;
#        (c) no single site carries more than half the weight.
#      E1 fails (a) and (c); E2 fails (c); E3 fails (b) (A-2 dropped for its value; the
#      WWTP footnote also marks D-1, which E3 keeps); E4 fails (b) (it drops the 12 rice-ND
#      sites, which are low-ratio sites by construction, so it is biased UP).
#      E5 = censored per-site geometric mean: ML fit of a normal to log10(rice/water) over
#      the 26 sites with both a water value and an informative ratio (14 exact, 12
#      left-censored at log10(MDL_rice/water)). E5 is the only candidate that passes all
#      three; the criteria were written knowing that, so the content of the rule is the
#      criteria, stated so they can be disputed. Uncertainty: bootstrap over the 26 sites
#      (2000 resamples) -> 90% interval. Sensitivity (not primary): the Kaplan-Meier median
#      on the flipped scale (the standard estimator for multiple detection limits), and the
#      leave-one-site-out range of E5.
#   Q2 (data) Does the per-site ratio behave like a BAF? A constant BAF means rice is
#      proportional to water: slope 1 of log10 rice on log10 water. Censored (Tobit)
#      regression over the same 26 sites (E5 is its b = 1 special case); bootstrap 90%
#      interval of b. If the interval excludes 1, the harvest void water (a 2-day snapshot)
#      is not proportional to what the plant took up over the season, and E5 is a typical
#      value with that noise in it, not a plant constant.
#   Q3 (model) Each variant's PFOA grain BAF (kim2019_censored.predict: mono, W2, lipid,
#      nstem_lipid; Cwo 1, fresh weight) against E1 and against E5 with its interval.
#   Q3b (model) Is the model linear in exposure over Kim's water range? Slope of log10
#      grain conc on log10 Cwo between 0.0008 and 1.76 ug/L (the lowest and highest PFOA
#      void water), per variant. If ~1, a data slope below 1 cannot come from the carrier.
#   Q4 (data) A-2's share of the void-water SUM-17 (Table 1). Kim's other congeners have
#      no per-site values, so this is the only check on whether their denominators are
#      A-2-weighted too.
#
#   DECISION RULES:
#   D1 kim2019_grain_baf() keeps the published ratio of means for ALL congeners: it is the
#      only estimator the other five have, and a chain-length series on mixed estimators
#      would compare unlike quantities. The per-site PFOA analysis is added beside it, not
#      substituted.
#   D2 The literature_params demo's L_Ph target moves from E1 to E5 if E5's 90% interval
#      excludes E1; otherwise it stays at E1.
#   D3 The "PFOA* used in the fit" labels are removed (factually wrong). PFOA is reported as
#      a separate OOS check (each variant vs E5 and its interval: inside / below / above);
#      it is NOT folded into the pre-registered 5-congener Kim RMSE of step 4a, which stays
#      the record.
#   D4 If Q4's A-2 share exceeds 50%, the absolute LEVEL of every Kim congener is recorded as
#      A-2-weighted (the chain-length SHAPE compares like with like; the level does not).
#
#   EXPECTATION (stated so it can be wrong): E5 lies below E4 (the dropped sites are
#      low-ratio sites), somewhere around 10-20, and its interval excludes E1 (so D2 moves
#      the demo target); b is well below 1 (rice varies ~19x across sites, water ~2000x);
#      the model is close to linear (Q3b slope > 0.9), so the data slope is not the carrier;
#      every variant sits below E5's interval. No expectation for Q4.
#
#   python validation/kim2019_pfoa_anchor.py
# =============================================================================
from __future__ import annotations

import contextlib
import csv
import io
import os
import sys

import numpy as np
from scipy import optimize, stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, HERE)

import literature_params as lp     # noqa: E402
import model_api as api            # noqa: E402
import kim2019_censored as K       # noqa: E402

TABLE = os.path.join(ROOT, "docs", "literature_db", "raw_si", "kim2019_table1_sites.csv")
MDL_RICE = lp.KIM2019_DETECTION["PFOA"][0]        # 0.088 ng/g
N_BOOT, SEED, LEVEL = 2000, 20261024, 0.90
CWO_LO, CWO_HI = 0.0008, 1.76                     # ug/L: lowest / highest PFOA void water


def load_sites(path=TABLE):
    with open(path) as f:
        return list(csv.DictReader(line for line in f if not line.startswith("#")))


def _num(x):
    return None if x in ("ND", "NS") else float(x)


def ratio_data(rows):
    """Per-site log10 water / log10 rice / censoring, over the sites that inform the ratio.

    Kept: water detected (26 sites). Rice detected -> exact; rice ND -> left-censored at the
    MDL. Dropped: C-6 (water ND as well: no information on the ratio) and E-2..E-4 (NS)."""
    out = []
    for r in rows:
        w, ri = r["water_PFOA_ngL"], r["rice_PFOA_ngg"]
        if w in ("ND", "NS"):
            continue
        out.append(dict(site=r["site"], wwtp=r["wwtp_near"] == "1",
                        x=np.log10(float(w) / 1000.0),                    # log10 ug/L
                        y=np.log10(MDL_RICE if ri == "ND" else float(ri)),  # log10 ug/kg
                        cens=ri == "ND"))
    return out


def estimators(rows):
    """E1-E4: the candidates computed before the pre-registration (data only)."""
    water = np.array([np.nan if r["water_PFOA_ngL"] == "NS" else (_num(r["water_PFOA_ngL"]) or 0.0)
                      for r in rows])
    rice = np.array([_num(r["rice_PFOA_ngg"]) or 0.0 for r in rows])
    has_w = ~np.isnan(water)
    not_a2 = np.array([r["site"] != "A-2" for r in rows])
    both = [(r["site"], _num(r["rice_PFOA_ngg"]), _num(r["water_PFOA_ngL"])) for r in rows
            if _num(r["rice_PFOA_ngg"]) and _num(r["water_PFOA_ngL"])]
    b = np.array([1000.0 * ri / w for _, ri, w in both])
    return dict(
        E1=1000.0 * rice.mean() / np.nanmean(water),
        E2=1000.0 * rice[has_w].mean() / water[has_w].mean(),
        E3=1000.0 * rice[not_a2].mean() / np.nanmean(water[not_a2]),
        E4=float(np.exp(np.mean(np.log(b)))), E4_median=float(np.median(b)),
        E4_range=(float(b.min()), float(b.max())), n_both=len(b),
        a2_water_share=float(1760.0 / np.nansum(water)))


def _tobit_nll(theta, x, y, cens, slope):
    if slope is None:
        a, b, ls = theta
    else:
        (a, ls), b = theta, slope
    s = np.exp(ls)
    mu = a + b * x
    z = (y - mu) / s
    return -(np.sum(stats.norm.logpdf(z[~cens]) - ls) + np.sum(stats.norm.logcdf(z[cens])))


def tobit(x, y, cens, slope=None):
    """Left-censored normal ML of y = a + b x. slope=1 gives E5 (a = mean log10 ratio)."""
    a0 = np.mean(y - (x if slope is not None else 0.0))
    th0 = [a0, np.log(np.std(y) + 0.1)] if slope is not None else [np.mean(y), 0.0,
                                                                     np.log(np.std(y) + 0.1)]
    res = optimize.minimize(_tobit_nll, th0, args=(x, y, cens, slope), method="Nelder-Mead",
                            options=dict(xatol=1e-8, fatol=1e-10, maxiter=20000))
    if slope is None:
        a, b, ls = res.x
    else:
        (a, ls), b = res.x, slope
    return dict(a=float(a), b=float(b), sigma=float(np.exp(ls)))


def km_median(v, cens):
    """Kaplan-Meier median of left-censored values (flip to right-censored, Helsel)."""
    z = -np.asarray(v, float)
    ev = ~np.asarray(cens)
    S = 1.0
    for t in np.unique(z[ev]):
        n = np.sum(z >= t)
        d = np.sum((z == t) & ev)
        S *= 1.0 - d / n
        if S <= 0.5:
            return float(-t)
    return float("nan")


def anchor(d):
    x = np.array([r["x"] for r in d])
    y = np.array([r["y"] for r in d])
    c = np.array([r["cens"] for r in d])
    lr = y - x                                             # log10 BAF [L/kg]
    e5 = tobit(x, y, c, slope=1.0)
    reg = tobit(x, y, c)
    rng = np.random.default_rng(SEED)
    boot_e5, boot_b = [], []
    for _ in range(N_BOOT):
        i = rng.integers(0, len(d), len(d))
        if c[i].all():
            continue
        boot_e5.append(tobit(x[i], y[i], c[i], slope=1.0)["a"])
        boot_b.append(tobit(x[i], y[i], c[i])["b"])
    q = [(1 - LEVEL) / 2, (1 + LEVEL) / 2]
    loo = [tobit(np.delete(x, k), np.delete(y, k), np.delete(c, k), slope=1.0)["a"]
           for k in range(len(d))]
    k_lo, k_hi = int(np.argmin(loo)), int(np.argmax(loo))
    return dict(E5=10 ** e5["a"], E5_sigma_log10=e5["sigma"],
                E5_ci=tuple(10 ** np.quantile(boot_e5, q)),
                km_median=10 ** km_median(lr, c),
                loo=((d[k_lo]["site"], 10 ** loo[k_lo]), (d[k_hi]["site"], 10 ** loo[k_hi])),
                slope=reg["b"], slope_ci=tuple(np.quantile(boot_b, q)),
                n=len(d), n_cens=int(c.sum()))


def profile_interval(d, level=LEVEL):
    """POST-HOC: profile-likelihood interval of E5 (not the pre-registered bootstrap)."""
    x = np.array([r["x"] for r in d])
    y = np.array([r["y"] for r in d])
    c = np.array([r["cens"] for r in d])

    def prof(a):
        return optimize.minimize_scalar(lambda ls: _tobit_nll([a, ls], x, y, c, 1.0),
                                        bounds=(-3.0, 2.0), method="bounded").fun
    a_hat = tobit(x, y, c, slope=1.0)["a"]
    best, crit = prof(a_hat), stats.chi2.ppf(level, 1) / 2.0
    lo = optimize.brentq(lambda a: prof(a) - best - crit, a_hat - 2.0, a_hat)
    hi = optimize.brentq(lambda a: prof(a) - best - crit, a_hat, a_hat + 2.0)
    return 10 ** lo, 10 ** hi


def model_baf(v, cwo=1.0):
    """Harvest PFOA grain BAF [L/kg fw] at constant Cwo (as kim2019_censored.predict)."""
    if v == "nstem_lipid":
        r = api.simulate_nstem_leaf("PFOA", Cwo=cwo, lipid_loading=True)
    else:
        kw = dict(mono=dict(f_xy_source="recommended"), fit=dict(f_xy_source="fit"),
                  lipid=dict(lipid_loading=True))[v]
        r = api.simulate("PFOA", Cwo=cwo, **kw)
    return float(r["baf_final"]["grain"])


def where(p, lo, hi):
    return "inside" if lo <= p <= hi else ("below" if p < lo else "above")


def main():
    rows = load_sites()
    est = estimators(rows)
    d = ratio_data(rows)
    an = anchor(d)
    w17 = [float(r["water_sum17_ngL"]) for r in rows if r["water_sum17_ngL"] != "NS"]
    a2_sum17 = float([r for r in rows if r["site"] == "A-2"][0]["water_sum17_ngL"]) / sum(w17)

    print("=" * 84)
    print("KIM 2019 PFOA GRAIN BAF: WHICH ESTIMATOR? — pre-registered in bedadbf (see header)")
    print("=" * 84)
    print("\nQ1 — candidate estimators [L/kg, pore-water basis]")
    print(f"  E1 published ratio of means (rice/30, water/27)      {est['E1']:7.2f}   fails a, c")
    print(f"  E2 ratio of means, the 27 sites with water           {est['E2']:7.2f}   fails c "
          f"(A-2 = {100 * est['a2_water_share']:.1f}% of water)")
    print(f"  E3 ratio of means without A-2                        {est['E3']:7.2f}   fails b")
    print(f"  E4 per-site geomean, both detected (n={est['n_both']})          {est['E4']:7.2f}   "
          f"fails b (median {est['E4_median']:.1f}, range {est['E4_range'][0]:.2f}-"
          f"{est['E4_range'][1]:.0f})")
    lo, hi = an["E5_ci"]
    print(f"  E5 censored per-site geomean (n={an['n']}, {an['n_cens']} censored)   {an['E5']:7.2f}   "
          f"admissible; {int(100 * LEVEL)}% bootstrap [{lo:.2f}, {hi:.2f}], "
          f"sd {an['E5_sigma_log10']:.2f} log10")
    print(f"     sensitivity: Kaplan-Meier median {an['km_median']:.2f}; leave-one-site-out "
          f"{an['loo'][0][1]:.2f} (without {an['loo'][0][0]}) to {an['loo'][1][1]:.2f} "
          f"(without {an['loo'][1][0]})")

    slo, shi = an["slope_ci"]
    print(f"\nQ2 — censored slope of log10 rice on log10 water: b = {an['slope']:.2f} "
          f"[{slo:.2f}, {shi:.2f}]  (a constant BAF means b = 1)")
    print(f"     -> the interval {'EXCLUDES' if shi < 1 or slo > 1 else 'includes'} 1")

    print("\nQ3 — model PFOA grain BAF (Cwo 1, fresh weight) vs E1 and E5")
    with contextlib.redirect_stdout(io.StringIO()):
        pred = {v: model_baf(v) for v in K.VARIANTS}
    for v in K.VARIANTS:
        p = pred[v]
        print(f"  {v:12}{p:8.3f}   vs E1 {np.log10(p / est['E1']):+.2f} log   vs E5 "
              f"{np.log10(p / an['E5']):+.2f} log   ({where(p, lo, hi)} E5's interval)")

    print(f"\nQ3b — model slope of log10 grain conc on log10 Cwo, {CWO_LO}-{CWO_HI} ug/L")
    with contextlib.redirect_stdout(io.StringIO()):
        slope_m = {v: 1.0 + np.log10(model_baf(v, CWO_HI) / model_baf(v, CWO_LO))
                   / np.log10(CWO_HI / CWO_LO) for v in K.VARIANTS}
    print("  " + "   ".join(f"{v} {s:.3f}" for v, s in slope_m.items()))

    print(f"\nQ4 — A-2 share of the void-water SUM-17: {100 * a2_sum17:.1f}% "
          f"(PFOA {100 * est['a2_water_share']:.1f}%)")

    d2 = not (lo <= est["E1"] <= hi)
    d4 = a2_sum17 > 0.5
    print("\n" + "=" * 84)
    print("DECISIONS (against the pre-registered rules)")
    print("=" * 84)
    print("  D1 kim2019_grain_baf() keeps the published ratio of means for all congeners;"
          " the per-site PFOA estimate is added beside it")
    print(f"  D2 demo L_Ph target: {'E1 -> E5' if d2 else 'stays E1'} "
          f"(E1 {est['E1']:.2f} is {'outside' if d2 else 'inside'} [{lo:.2f}, {hi:.2f}])")
    print("  D3 'PFOA* used in the fit' labels removed; PFOA reported as a separate OOS check: "
          + ", ".join(f"{v} {where(pred[v], lo, hi)}" for v in K.VARIANTS))
    print(f"  D4 Kim series LEVEL A-2-weighted: {'YES' if d4 else 'no'} "
          f"(A-2 SUM-17 share {100 * a2_sum17:.1f}%)")
    plo, phi = profile_interval(d)
    print("\nPOST-HOC (not part of the rule)")
    print(f"  D2 is decided by a 3% margin (E1 {est['E1']:.2f} vs bootstrap low {lo:.2f}). The bootstrap"
          " low end stays 4.54-4.79 over 10 seeds (E1 outside in all), but the profile-likelihood")
    print(f"  {int(100 * LEVEL)}% interval is [{plo:.2f}, {phi:.2f}], which CONTAINS E1: robust to the seed,"
          " not to the interval method. What does not depend on it: E5 is ~2x E1 and ~2.7x below E4.")
    print("  Why the means are low: A-2, the highest-water site, has the lowest ratio (0.98) -- the same"
          f" b < 1 pattern; a ratio of means is weighted toward high-water sites.")
    return dict(est=est, anchor=an, pred=pred, slope_model=slope_m, a2_sum17=a2_sum17,
                d2=d2, d4=d4, profile=(plo, phi))


if __name__ == "__main__":
    main()
