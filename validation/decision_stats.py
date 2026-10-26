"""Decision rule without a chosen margin: a paired bootstrap of the score difference.

Every pre-registered comparison in this repo used to decide by a fixed margin (0.10 log10,
0.3 log10, ...). A margin is a number with no source. This module replaces it with the
sampling uncertainty of the comparison itself: the scored items (congeners, or
congener x tissue cells) are resampled with replacement, both configurations' log10 RMSE is
recomputed on the SAME resample, and the 95 % percentile interval of the difference decides
(Efron & Tibshirani 1993, An Introduction to the Bootstrap, ch. 13; the conventional
two-sided 0.05 level).

    better      the interval of RMSE(B) - RMSE(A) lies entirely below 0
    worse       ... entirely above 0
    tie         the interval contains 0 (the data cannot tell them apart)

Non-inferiority is "not worse"; a change of default needs "better" on the endpoint it
targets and "not worse" on every other scored endpoint. The number of resamples is not a
tuning knob: `decide` doubles it until the verdict and the interval ends (to 3 decimals)
stop changing.

Interval-censored items (an observation known only as [lo, hi]) score 0 inside the
interval and the log distance to the nearer end outside it, as in kim2019_censored.
"""
import numpy as np

LEVEL = 0.95
B0 = 2000          # first resample count; doubled until the decision is stable


def rmse(err):
    err = np.asarray(err, float)
    return float(np.sqrt(np.mean(np.square(err))))


def _boot_diff(err_a, err_b, n_boot, rng):
    a2, b2 = np.square(err_a), np.square(err_b)
    n = len(a2)
    idx = rng.integers(0, n, size=(n_boot, n))
    return np.sqrt(b2[idx].mean(axis=1)) - np.sqrt(a2[idx].mean(axis=1))


def paired_interval(err_a, err_b, n_boot=B0, seed=0):
    """95 % percentile interval of RMSE(b) - RMSE(a) over item resamples."""
    err_a = np.asarray(err_a, float)
    err_b = np.asarray(err_b, float)
    if err_a.shape != err_b.shape or err_a.ndim != 1 or len(err_a) < 2:
        raise ValueError("need two equal-length 1-D error vectors with >= 2 items")
    d = _boot_diff(err_a, err_b, n_boot, np.random.default_rng(seed))
    lo, hi = np.quantile(d, [(1 - LEVEL) / 2, 1 - (1 - LEVEL) / 2])
    return float(lo), float(hi)


def verdict(lo, hi):
    return "better" if hi < 0 else ("worse" if lo > 0 else "tie")


def decide(err_a, err_b, seed=0, max_boot=256000):
    """Verdict of B against A, with the resample count doubled until stable.

    Returns dict(delta=RMSE(b)-RMSE(a), lo, hi, verdict, n_boot)."""
    n = B0
    lo, hi = paired_interval(err_a, err_b, n, seed)
    while True:
        n2 = 2 * n
        lo2, hi2 = paired_interval(err_a, err_b, n2, seed + 1)
        stable = (verdict(lo, hi) == verdict(lo2, hi2)
                  and round(lo, 3) == round(lo2, 3) and round(hi, 3) == round(hi2, 3))
        n, lo, hi = n2, lo2, hi2
        if stable or n >= max_boot:
            break
    return dict(delta=rmse(err_b) - rmse(err_a), lo=lo, hi=hi, verdict=verdict(lo, hi),
                n_boot=n)
