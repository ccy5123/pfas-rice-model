#!/usr/bin/env python3
"""
Bayesian inverse demo — recover the EXPOSURE from tissue concentrations
=======================================================================

Question (from the modelling discussion): given a time series of tissue
concentrations C_k(t) and the growth curve M_k(t), can we run a Bayesian
parameter estimation to infer the transpiration stream Q_TP(t) and the
pore-water exposure C_w^o(t)?

This demo answers it concretely on synthetic data, and — more importantly —
EXPOSES THE STRUCTURAL IDENTIFIABILITY of the problem, which is the real
subtlety:

  * `qtp_scale` (a multiplicative transpiration scale on the measured FAO-56
    shape) and `f_xy` (the root->xylem loading factor) enter the SHOOT only
    through their PRODUCT `qtp_scale * f_xy`. So from tissue concentrations they
    are NOT separately identifiable — the posterior collapses onto a ridge
    (correlation ~ -1).  To get Q_TP absolutely you must fix f_xy independently
    (xylem-sap), or vice versa.
  * `cwo_level` and the root-uptake conductance `kappa_d` trade off the same way
    (the linear-regime influx is ~ g_in * C_w^o), another ridge.
  * With the TRANSPORT parameters held fixed, `(qtp_scale, cwo_level)` IS
    identifiable — `cwo_level` sets the overall level (root), `qtp_scale` the
    shoot/root ratio — so the two act on different directions of the data.

Method: a Laplace (Gauss-Newton) posterior around the MAP fit in log space.
Observations are the final root/straw/grain BAFs (positive, span orders → log
space, as in src/calibration.py). The posterior covariance ~ (JᵀJ)⁻¹·σ² makes
the ridge quantitative: a near-±1 parameter correlation + an inflated marginal
std == non-identifiable. Optionally cross-checked with an emcee MCMC if present.

Run:  python validation/bayesian_inverse_demo.py
"""
import os
import sys

import numpy as np
from scipy.optimize import least_squares

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import model_api as api          # noqa: E402
import forcing_rice as fr        # noqa: E402
import growth_rice as gr         # noqa: E402

CONGENER = "PFOA"
SEASON = 120.0
N_T = 161
SIGMA = 0.05                      # log10 observation noise / weight
_RNG = np.random.default_rng(0)

# transport baseline (recommended monotone f_xy + its kappa_d) for CONGENER
_c = api._CONG[CONGENER]
_BASE_FXY, _, _BASE_KD, _, _ = api._transport_defaults(_c, "recommended", False)

_T = np.linspace(0.0, SEASON, N_T)
_BASE_Q = fr.Q_TP(_T, SEASON)
_b = gr.organ_biomass(_T, SEASON)
_M = np.maximum(np.column_stack([_b["root"], _b["stem"], _b["leaf"], _b["grain"]]), 1e-4)

DEFAULTS = dict(qtp_scale=1.0, cwo_level=1.0, f_xy_mult=1.0, kappa_d_mult=1.0)


def predict(knobs):
    """Final root/straw/grain BAF for a set of exposure/transport knobs."""
    k = {**DEFAULTS, **knobs}
    drivers = dict(t=_T, Cwo=np.full_like(_T, k["cwo_level"]),
                   Qtp=k["qtp_scale"] * _BASE_Q, M=_M)
    r = api.simulate(CONGENER, drivers=drivers,
                     f_xy_override=_BASE_FXY * k["f_xy_mult"],
                     kappa_d_override=_BASE_KD * k["kappa_d_mult"])
    return np.array([r["baf_final"]["root"], r["straw_baf"], r["baf_final"]["grain"]])


def _jacobian_at(free, point):
    """d(log10 pred_i)/d(log10 θ_j) at `point` (full knobs dict), central differences.
    This is the FISHER-information Jacobian -- a property of the MODEL at the true
    operating point, so the identifiability diagnostic does not depend on the noise
    realisation (unlike a Jacobian read off a noisy MAP)."""
    base = np.array([np.log10(point[k]) for k in free])

    def logpred(logp):
        knobs = {**point}
        for k, lp in zip(free, logp):
            knobs[k] = float(10.0 ** lp)
        return np.log10(np.maximum(predict(knobs), 1e-12))

    h, J = 1e-3, np.zeros((3, len(free)))
    for j in range(len(free)):
        hi, lo = base.copy(), base.copy()
        hi[j] += h
        lo[j] -= h
        J[:, j] = (logpred(hi) - logpred(lo)) / (2 * h)
    return J / SIGMA


def run_scenario(name, free, truth=None, noise=True):
    """Recover `free` from synthetic root/straw/grain obs (MAP, log space) AND report
    the identifiability ridge from the Fisher Jacobian AT TRUTH (noise-independent)."""
    truth = {**DEFAULTS, **(truth or {})}
    obs = predict({k: truth[k] for k in DEFAULTS})
    if noise:
        obs = obs * np.exp(_RNG.normal(0, SIGMA, size=obs.shape))

    def resid(logp):
        knobs = {k: truth[k] for k in DEFAULTS}
        for k, lp in zip(free, logp):
            knobs[k] = float(10.0 ** lp)
        pred = np.maximum(predict(knobs), 1e-12)
        return (np.log10(pred) - np.log10(obs)) / SIGMA

    x0 = np.array([np.log10(0.5)] * len(free))            # start AWAY from truth (=0)
    res = least_squares(resid, x0, method="lm", diff_step=1e-3, xtol=1e-12, ftol=1e-12)
    J = _jacobian_at(free, {k: truth[k] for k in DEFAULTS})   # identifiability AT TRUTH
    JTJ = J.T @ J
    try:
        cov = np.linalg.inv(JTJ) * SIGMA ** 2             # Laplace posterior covariance
    except np.linalg.LinAlgError:
        cov = np.linalg.pinv(JTJ) * SIGMA ** 2
    sd = np.sqrt(np.clip(np.diag(cov), 0, None))
    with np.errstate(invalid="ignore", divide="ignore"):
        corr = cov / np.outer(sd, sd)
    recovered = {k: float(10.0 ** x) for k, x in zip(free, res.x)}
    cond = float(np.linalg.cond(JTJ))
    off = float(corr[0, 1]) if len(free) == 2 else float("nan")
    prod = float(np.prod([recovered[k] for k in free]))   # the product of the free knobs

    print(f"\n[{name}]  free = {free}")
    for k in free:
        print(f"    {k:12s} truth {truth[k]:6.3f}  ->  recovered {recovered[k]:7.3f}"
              f"   (±{sd[free.index(k)]:.2f} log10, local)")
    print(f"    posterior corr(θ1,θ2) = {off:+.3f}   cond(JᵀJ) = {cond:.1f}"
          f"   -> {'RIDGE: non-identifiable' if abs(off) > 0.95 else 'identifiable'}")
    if abs(off) > 0.95:
        if abs(np.log10(max(prod, 1e-9))) < 0.1:          # a clean product invariant (linear shoot advection)
            print(f"    BUT the PRODUCT θ1·θ2 = {prod:.3f} (truth 1.000) is recovered "
                  "-> only the product is constrained")
        else:                                             # nonlinear uptake -> ridge without an exact product
            print("    ridge in the data, but root uptake is NONLINEAR (GHK + saturable carrier) so "
                  "there is no exact product invariant -> the individual values just slide along the "
                  "ridge (here biased off truth) and are unconstrained")
    return dict(name=name, free=free, truth=truth, recovered=recovered,
                corr=off, cond=cond, sd=sd.tolist(), product=prod)


def emcee_posterior(free, truth=None, nstep=600, nwalk=16, discard=200):
    """Full MCMC posterior for the `free` params (affine-invariant ensemble sampler).
    Returns a dict with the flat chain, per-param medians/16-84% bands, the product
    posterior and the sample correlation -- a sampling cross-check of the Laplace
    result. None if emcee is not installed."""
    try:
        import emcee
    except Exception:
        print("\n(emcee not installed — Laplace posterior only; "
              "`pip install emcee` for the full-MCMC cross-check)")
        return None
    truth = {**DEFAULTS, **(truth or {})}
    obs = predict({k: truth[k] for k in DEFAULTS}) * np.exp(_RNG.normal(0, SIGMA, 3))

    def logpost(logp):
        if np.any(np.abs(logp) > 1.5):                    # broad flat prior, ±1.5 dex
            return -np.inf
        knobs = {k: truth[k] for k in DEFAULTS}
        for k, lp in zip(free, logp):
            knobs[k] = float(10.0 ** lp)
        pred = np.maximum(predict(knobs), 1e-12)
        return -0.5 * float(np.sum(((np.log10(pred) - np.log10(obs)) / SIGMA) ** 2))

    nd = len(free)
    p0 = _RNG.normal(0, 0.05, size=(nwalk, nd))
    s = emcee.EnsembleSampler(nwalk, nd, logpost)
    s.run_mcmc(p0, nstep, progress=False)
    chain = s.get_chain(discard=discard, flat=True)        # log10 space
    lin = 10.0 ** chain
    med = {k: float(np.median(lin[:, i])) for i, k in enumerate(free)}
    lo = {k: float(np.percentile(lin[:, i], 16)) for i, k in enumerate(free)}
    hi = {k: float(np.percentile(lin[:, i], 84)) for i, k in enumerate(free)}
    corr = float(np.corrcoef(chain.T)[0, 1]) if nd == 2 else float("nan")
    prod = float(np.median(np.prod(lin, axis=1)))
    out = dict(free=free, n=int(chain.shape[0]), median=med, lo=lo, hi=hi,
               corr=corr, product=prod, truth=truth)

    print(f"\n[emcee MCMC]  free={free}  ({out['n']} samples)")
    for k in free:
        print(f"    {k:12s} median {med[k]:6.3f}  [{lo[k]:.3f}, {hi[k]:.3f}] (16-84%)"
              f"   truth {truth[k]:.3f}")
    print(f"    sample corr = {corr:+.3f}   product (med) = {prod:.3f}"
          f"   -> {'RIDGE (only the product is pinned)' if abs(corr) > 0.95 else 'identifiable'}")
    return out


def maybe_emcee(nstep=120, nwalk=8):
    """Run the MCMC cross-check on BOTH the well-posed case (A) and the product
    ridge (B), confirming the Laplace verdicts with a real sampled posterior.
    OPT-IN (the forward ODE is ~0.7 s/eval, so a chain is minutes): only invoked by
    the `--emcee` flag."""
    a = emcee_posterior(["qtp_scale", "cwo_level"], nstep=nstep, nwalk=nwalk)  # recovers
    if a is None:
        return None
    emcee_posterior(["qtp_scale", "f_xy_mult"], nstep=nstep, nwalk=nwalk)      # ridge
    return a


def run(do_emcee=False):
    print("Bayesian inverse demo — infer exposure from tissue C(t)\n"
          f"  congener={CONGENER}, obs = final root/straw/grain BAF, "
          f"σ={SIGMA} log10")
    results = [
        run_scenario("A  WELL-POSED: exposure, transport fixed",
                     ["qtp_scale", "cwo_level"]),
        run_scenario("B  RIDGE: transpiration vs root->xylem loading",
                     ["qtp_scale", "f_xy_mult"]),
        run_scenario("C  RIDGE: pore water vs root-uptake conductance",
                     ["cwo_level", "kappa_d_mult"]),
    ]
    if do_emcee:
        maybe_emcee()
    else:
        print("\n(Laplace posterior above; pass --emcee for the full-MCMC cross-check "
              "— a few minutes, the ODE is ~0.7 s/sample.)")
    print("\nTakeaway: (A) with transport fixed, Q_TP-scale & Cwo-level are the BEST-"
          "conditioned (cond ~90) and recover from tissue C(t) -- they act on different "
          "directions (Cwo sets the root level, Q_TP the shoot/root ratio). (B) Q_TP & "
          "f_xy collapse onto a product ridge (corr ~-1, cond ~500): only the PRODUCT "
          "Q_TP·f_xy is constrained. (C) Cwo & root-uptake conductance are even more "
          "degenerate (cond ~1e5) but, uptake being nonlinear, with no clean product "
          "invariant. So pinning Q_TP or Cwo absolutely needs an independent measurement "
          "(xylem sap / a pore-water probe), exactly as the model notes. Multi-compartment "
          "data (root vs straw vs grain) only PARTIALLY break B, because Q_TP also sets the "
          "intra-shoot advection rates that f_xy does not.")
    return results


if __name__ == "__main__":
    run(do_emcee="--emcee" in sys.argv)
