"""
Tier-1 calibration machinery for the PFAS rice uptake model
===========================================================

Fits the BAF-identifiable (Tier-1) parameters of the plant module to observed
tissue bioaccumulation factors (root / straw / grain, optionally stem & leaf
separately).  The objective is a weighted least-squares fit in *log space*
(BAFs span orders of magnitude and are positive), optimised with
``scipy.optimize`` under box constraints -- no extra dependencies.

What is identifiable
--------------------
Per the model report, BAF data constrain the lumped groups ``B_k``,
``g_in/g_out``, ``f_xy`` (TSCF), ``Pi = Q_Phl*L_Ph/Q_TP`` and ``phi``.  In the
low-concentration (linear) limit the channel (``kappa_d``) and carrier
(``Vmax_in``) enter the root BAF only through the lumped influx conductance
``g_in``, so their split is poorly constrained by BAF data and is best resolved
with inhibitor / concentration-series experiments.  (Note the passive channel
alone *excludes* the anion -- root BAF < 1 -- so observed accumulation requires
the carrier.)

Usage
-----
    params = [Param("f_xy", 1e-3, 1.0), Param("L_Ph", 1e-4, 0.5),
              Param("kappa_d", 1e-2, 5.0), Param("K_prot", 1.0, 500.0)]
    obs = [ObservedBAF("root", 9.0), ObservedBAF("straw", 4.3),
           ObservedBAF("grain", 2.3)]
    result = calibrate(model, params, obs)

``model`` is a configured :class:`RiceUptakeModel`; calibration mutates a copy.
"""
from __future__ import annotations

import copy
from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares, differential_evolution

from pfas_rice_plant_module import (
    RiceUptakeModel, ROOT, STEM, LEAF, FRUIT,
)

# parameter name -> object it lives on ("cmpd" attribute or "model" attribute)
_CMPD_PARAMS = {"f_xy", "L_Ph", "kappa_d", "Vmax_in", "Vmax_out",
                "Km_in", "Km_out", "K_prot", "K_PL", "K_cw"}
_MODEL_PARAMS = {"phi", "T_C_Ph"}

TISSUES = ("root", "stem", "leaf", "straw", "grain")


@dataclass
class Param:
    """A free parameter with box bounds, optimised in log space by default."""
    name: str
    low: float
    high: float
    log: bool = True

    def to_x(self, value: float) -> float:
        return np.log10(value) if self.log else value

    def to_value(self, x: float) -> float:
        return float(10.0 ** x) if self.log else float(x)

    @property
    def bounds_x(self) -> tuple[float, float]:
        return (self.to_x(self.low), self.to_x(self.high))


@dataclass
class ObservedBAF:
    """An observed tissue BAF (C_tissue/C_w^o) with a log-space 1-sigma."""
    tissue: str
    value: float
    sigma: float = 0.3            # ~ relative error in log10 space
    compound: str = ""            # optional label (chain-length series)

    def __post_init__(self):
        if self.tissue not in TISSUES:
            raise ValueError(f"tissue must be one of {TISSUES}, got {self.tissue!r}")


# ---------------------------------------------------------------------------
# parameter <-> model plumbing
# ---------------------------------------------------------------------------
def set_param(model: RiceUptakeModel, name: str, value: float) -> None:
    if name in _CMPD_PARAMS:
        setattr(model.cmpd, name, value)
    elif name in _MODEL_PARAMS:
        setattr(model, name, value)
    else:
        raise KeyError(f"unknown calibratable parameter {name!r}; "
                       f"known: {sorted(_CMPD_PARAMS | _MODEL_PARAMS)}")


def get_param(model: RiceUptakeModel, name: str) -> float:
    if name in _CMPD_PARAMS:
        return float(getattr(model.cmpd, name))
    if name in _MODEL_PARAMS:
        return float(getattr(model, name))
    raise KeyError(f"unknown calibratable parameter {name!r}")


def apply_params(model: RiceUptakeModel, params: list[Param], theta_x) -> None:
    for p, x in zip(params, theta_x):
        set_param(model, p.name, p.to_value(x))


# ---------------------------------------------------------------------------
# prediction
# ---------------------------------------------------------------------------
def predict_bafs(model: RiceUptakeModel, t: np.ndarray | None = None) -> dict:
    """Final-time tissue BAFs, including the mass-weighted ``straw``."""
    if t is None:
        t = model.inputs.t
    sol = model.solve(t)
    C = sol.y[:, -1]
    Cwo = model.inputs.Cwo_(t[-1])
    Mf = model.inputs.M_(t[-1])
    straw = (C[STEM] * Mf[STEM] + C[LEAF] * Mf[LEAF]) / (Mf[STEM] + Mf[LEAF])
    return {
        "root": C[ROOT] / Cwo, "stem": C[STEM] / Cwo, "leaf": C[LEAF] / Cwo,
        "straw": straw / Cwo, "grain": C[FRUIT] / Cwo,
    }


def _residuals(theta_x, model, params, obs, t):
    apply_params(model, params, theta_x)
    try:
        pred = predict_bafs(model, t)
    except Exception:
        return np.full(len(obs), 1e3)        # penalise failed integrations
    r = []
    for o in obs:
        p = max(pred[o.tissue], 1e-12)
        r.append((np.log10(p) - np.log10(o.value)) / o.sigma)
    return np.asarray(r)


# ---------------------------------------------------------------------------
# calibration
# ---------------------------------------------------------------------------
@dataclass
class CalibrationResult:
    params: list[Param]
    x: np.ndarray                      # best-fit in optimiser (log) space
    values: dict                       # name -> physical value
    cost: float                        # 0.5 * sum(residual^2)
    predicted: dict                    # tissue -> BAF at the optimum
    observed: dict                     # tissue -> observed BAF
    success: bool
    message: str = ""

    def report(self) -> str:
        lines = [f"calibration {'OK' if self.success else 'FAILED'} "
                 f"(cost={self.cost:.4g}): {self.message}", "  parameters:"]
        for p in self.params:
            lines.append(f"    {p.name:10s} = {self.values[p.name]:.4g}  "
                         f"[{p.low:g}, {p.high:g}]")
        lines.append("  fit (tissue: pred vs obs):")
        for tis, ob in self.observed.items():
            lines.append(f"    {tis:6s} {self.predicted[tis]:8.3f}  vs {ob:8.3f}")
        return "\n".join(lines)


def calibrate(model: RiceUptakeModel, params: list[Param], obs: list[ObservedBAF],
              t: np.ndarray | None = None, x0: np.ndarray | None = None,
              global_search: bool = False, seed: int | None = None,
              de_maxiter: int = 60) -> CalibrationResult:
    """Fit ``params`` to observed tissue BAFs ``obs``.

    Works on a deep copy of ``model`` (the input is left untouched).  With
    ``global_search=True`` a ``differential_evolution`` sweep seeds a final
    local ``least_squares`` polish -- recommended when ``x0`` is unknown.
    """
    work = copy.deepcopy(model)
    if t is None:
        t = work.inputs.t
    lo = np.array([p.bounds_x[0] for p in params])
    hi = np.array([p.bounds_x[1] for p in params])

    if global_search:
        de = differential_evolution(
            lambda xx: 0.5 * float(np.sum(_residuals(xx, work, params, obs, t) ** 2)),
            bounds=list(zip(lo, hi)), seed=seed, tol=1e-7, polish=False,
            maxiter=de_maxiter, init="sobol")
        x_start = de.x
    elif x0 is not None:
        x_start = np.array([p.to_x(v) for p, v in zip(params, x0)])
    else:
        x_start = 0.5 * (lo + hi)

    # diff_step must be large enough that the BAF change from perturbing a
    # parameter exceeds the ODE solver's tolerance floor (~1e-6); the default
    # sqrt(eps) step is swamped by integration noise and stalls the optimiser.
    res = least_squares(_residuals, x_start, bounds=(lo, hi),
                        args=(work, params, obs, t), method="trf",
                        diff_step=1e-2, xtol=1e-10, ftol=1e-10)
    apply_params(work, params, res.x)
    pred = predict_bafs(work, t)
    return CalibrationResult(
        params=params, x=res.x,
        values={p.name: p.to_value(x) for p, x in zip(params, res.x)},
        cost=float(res.cost), predicted={o.tissue: pred[o.tissue] for o in obs},
        observed={o.tissue: o.value for o in obs}, success=bool(res.success),
        message=str(res.message))


# ---------------------------------------------------------------------------
# data loading
# ---------------------------------------------------------------------------
def load_baf_csv(path: str) -> list[ObservedBAF]:
    """Load observed BAFs from a CSV with header ``compound,tissue,baf[,sigma]``."""
    out: list[ObservedBAF] = []
    data = np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding="utf-8")
    data = np.atleast_1d(data)
    has_sigma = "sigma" in (data.dtype.names or ())
    for row in data:
        out.append(ObservedBAF(
            tissue=str(row["tissue"]).strip(), value=float(row["baf"]),
            sigma=float(row["sigma"]) if has_sigma else 0.3,
            compound=str(row["compound"]).strip() if "compound" in data.dtype.names else ""))
    return out


# ---------------------------------------------------------------------------
# demos: synthetic recovery + identifiability
# ---------------------------------------------------------------------------
def _demo_model() -> tuple[RiceUptakeModel, np.ndarray]:
    from pfas_rice_plant_module import (
        Environment, Compound, Compartment, PlantInputs, _logistic)
    t = np.linspace(0.0, 120.0, 481)
    Cwo = np.full_like(t, 1.0)
    Qtp = 0.05 + 0.35 * np.exp(-((t - 75.0) ** 2) / (2 * 25.0 ** 2))
    M = np.column_stack([
        _logistic(t, 1e-3, 0.030, 0.10, 20.0), _logistic(t, 1e-3, 0.040, 0.10, 25.0),
        _logistic(t, 1e-3, 0.050, 0.12, 30.0), _logistic(t, 1e-5, 0.025, 0.18, 80.0)])
    inputs = PlantInputs(t=t, Cwo=Cwo, Qtp=Qtp, M=M)
    cmpd = Compound(name="PFOA", K_prot=50.0, K_PL=100.0, K_cw=20.0, kappa_d=0.5,
                    Vmax_in=20.0, Km_in=5.0, Vmax_out=8.0, Km_out=5.0, L_Ph=0.005, f_xy=0.02)
    comps = [Compartment("root", 0.70, 0.05, 0.010, 0.30),
             Compartment("stem", 0.80, 0.01, 0.005, 0.08),
             Compartment("leaf", 0.80, 0.03, 0.020, 0.04, S=20.0),
             Compartment("grain", 0.15, 0.08, 0.010, 0.10, S=2.0)]
    return RiceUptakeModel(env=Environment(), cmpd=cmpd, comps=comps, inputs=inputs), t


def synthetic_recovery(noise=0.0, seed=0, global_search=True):
    """Generate pseudo-observations from known params and recover them.

    A well-posed (determined) problem: three Tier-1 parameters that act on
    distinct ratios -- ``kappa_d`` (overall uptake -> root level), ``f_xy``
    (translocation -> straw/root) and ``L_Ph`` (grain feeding -> grain/straw)
    -- fit to the three observed BAFs root/straw/grain."""
    model, t = _demo_model()
    truth = {"f_xy": 0.02, "L_Ph": 0.005, "kappa_d": 0.7}
    for k, v in truth.items():
        set_param(model, k, v)
    pred = predict_bafs(model, t)
    rng = np.random.default_rng(seed)
    obs = []
    for tis in ("root", "straw", "grain"):
        val = pred[tis] * (np.exp(rng.normal(0, noise)) if noise else 1.0)
        obs.append(ObservedBAF(tis, val, sigma=max(noise, 0.1)))
    params = [Param("f_xy", 1e-3, 1.0), Param("L_Ph", 1e-4, 0.5),
              Param("kappa_d", 1e-2, 5.0)]
    result = calibrate(model, params, obs, global_search=global_search, seed=seed)
    return truth, result


if __name__ == "__main__":
    print("=== synthetic parameter recovery (noise-free, local solver) ===")
    truth, result = synthetic_recovery(noise=0.0, global_search=False)
    print(result.report())
    print("  truth:", {k: round(v, 4) for k, v in truth.items()})
    print("  (pass global_search=True for a differential-evolution global sweep)")

    print("\n=== with 10% log-normal observation noise ===")
    truth, result = synthetic_recovery(noise=0.10, seed=3, global_search=False)
    print("  recovered: " + ", ".join(f"{k}={result.values[k]:.4g}" for k in truth)
          + "  (truth: " + ", ".join(f"{k}={v:.4g}" for k, v in truth.items()) + ")")
