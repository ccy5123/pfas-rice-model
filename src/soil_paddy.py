"""
Paddy soil sub-model and input adapters for the PFAS rice uptake model
======================================================================

This module supplies the *soil side* of the Method A (loose, one-way) coupling:
it turns a soil PFAS inventory into the pore-water free concentration
``C_w^o(t)`` that drives the plant module, and provides loaders that map
external soil-model output (HYDRUS-1D / Phydrus, or a generic CSV) onto the
three arrays the plant module consumes (``Cwo``, ``Qtp``, ``M``).

Why a Freundlich isotherm
-------------------------
In flooded paddy soil PFAS partitions between pore water (the bioavailable,
free anion ``C_w``) and the solid phase.  Sorption of PFAS is routinely
non-linear (Freundlich exponent ``n`` typically 0.7--1.0), and it is *redox
dependent*: flooding drives the soil anaerobic, which changes organic-matter
conformation, iron-oxide surfaces and pH, and hence the sorption strength
``K_F``.  A linear ``K_d`` (the ``n=1`` special case) cannot capture either
effect.  The model therefore stores the solid-phase load as

    S(C_w) = K_F * C_w**n            [ug / kg dry soil]

and the total PFAS per kg dry soil (solid + the water held in that soil) as

    C_T(C_w) = K_F * C_w**n + theta_g * C_w

with ``theta_g`` the gravimetric water content [L water / kg dry soil].
Given a total inventory ``C_T`` the bioavailable pore-water concentration is
recovered by inverting this strictly-increasing relation (Section ``pore_water``).

Unit system (matches ``pfas_rice_plant_module``)
------------------------------------------------
    aqueous conc    ug/L        (C_w, C_w^o)
    solid-phase     ug/kg dry   (S, C_T)
    water content   L/kg dry    (theta_g, gravimetric)
    K_F             ug^(1-n) L^n / kg     (so K_F*C_w**n is ug/kg)
    time            day
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from scipy.optimize import brentq

from pfas_rice_plant_module import PlantInputs


# ---------------------------------------------------------------------------
# Freundlich sorption
# ---------------------------------------------------------------------------
@dataclass
class FreundlichSoil:
    """Freundlich pore-water <-> solid partitioning for one redox state.

    Parameters
    ----------
    K_F : Freundlich capacity [ug^(1-n) L^n / kg dry].
    n   : Freundlich exponent [-] (1.0 = linear K_d).
    theta_g : gravimetric water content [L water / kg dry soil].
    name : label (e.g. "drained/aerobic", "flooded/anaerobic").
    """
    K_F: float
    n: float = 1.0
    theta_g: float = 0.30
    name: str = "soil"

    def sorbed(self, Cw: float | np.ndarray) -> float | np.ndarray:
        """Solid-phase load S = K_F * C_w**n  [ug/kg dry]."""
        Cw = np.asarray(Cw, dtype=float)
        return self.K_F * np.power(np.clip(Cw, 0.0, None), self.n)

    def total(self, Cw: float | np.ndarray) -> float | np.ndarray:
        """Total PFAS per kg dry soil = sorbed + water-held  [ug/kg dry]."""
        Cw = np.asarray(Cw, dtype=float)
        return self.sorbed(Cw) + self.theta_g * Cw

    def Kd_eff(self, Cw: float) -> float:
        """Effective (secant) partition S/C_w = K_F*C_w**(n-1) [L/kg].
        Concentration-dependent unless n == 1."""
        if Cw <= 0.0:
            return float("inf") if self.n < 1.0 else self.K_F
        return self.K_F * Cw ** (self.n - 1.0)

    def pore_water(self, C_total: float) -> float:
        """Invert C_T(C_w) = K_F C_w**n + theta_g C_w for the free conc C_w.

        C_T is strictly increasing in C_w (>=0), so the root is unique."""
        if C_total <= 0.0:
            return 0.0
        if abs(self.n - 1.0) < 1e-12:                 # linear limit, closed form
            return C_total / (self.K_F + self.theta_g)
        # bracket [0, hi], expanding hi until total(hi) >= C_total
        hi = max(C_total / max(self.theta_g, 1e-12), 1.0)
        for _ in range(60):
            if self.total(hi) >= C_total:
                break
            hi *= 2.0
        return float(brentq(lambda c: float(self.total(c)) - C_total, 0.0, hi,
                            xtol=1e-12, rtol=1e-10))

    def pore_water_series(self, C_total: Sequence[float]) -> np.ndarray:
        """Vectorised :meth:`pore_water` over a total-inventory time series."""
        return np.array([self.pore_water(float(ct)) for ct in np.asarray(C_total)])


@dataclass
class PaddyRedox:
    """Two-state (drained/flooded) redox switch for paddy sorption.

    A boolean flooding schedule selects the anaerobic (flooded) or aerobic
    (drained) :class:`FreundlichSoil` at each time step.  Anaerobic conditions
    are taken to *weaken* sorption by default (lower ``K_F`` -> more
    bioavailable PFAS), an illustrative choice the user can override.
    """
    drained: FreundlichSoil
    flooded: FreundlichSoil

    def soil_at(self, is_flooded: bool) -> FreundlichSoil:
        return self.flooded if is_flooded else self.drained

    def pore_water_series(self, C_total: Sequence[float],
                          flooded: Sequence[bool]) -> np.ndarray:
        C_total = np.asarray(C_total, dtype=float)
        flooded = np.asarray(flooded, dtype=bool)
        if C_total.shape != flooded.shape:
            raise ValueError("C_total and flooded must have the same shape")
        return np.array([self.soil_at(f).pore_water(float(ct))
                         for ct, f in zip(C_total, flooded)])


# ---------------------------------------------------------------------------
# soil -> plant input builders
# ---------------------------------------------------------------------------
def inputs_from_soil(t: np.ndarray, C_total, Qtp, M,
                     soil: FreundlichSoil | PaddyRedox,
                     flooded=None) -> PlantInputs:
    """Build :class:`PlantInputs` from a soil total-inventory scenario.

    ``C_total`` [ug/kg dry] is converted to the pore-water free concentration
    ``C_w^o(t)`` via Freundlich inversion; ``Qtp`` (transpiration) and ``M``
    (tissue masses) pass through unchanged.
    """
    t = np.asarray(t, dtype=float)
    if isinstance(soil, PaddyRedox):
        if flooded is None:
            raise ValueError("PaddyRedox requires a `flooded` boolean schedule")
        Cwo = soil.pore_water_series(C_total, flooded)
    else:
        Cwo = soil.pore_water_series(C_total)
    return PlantInputs(t=t, Cwo=Cwo, Qtp=np.asarray(Qtp, float), M=np.asarray(M, float))


def load_inputs_csv(path: str) -> PlantInputs:
    """Load external soil/growth time series into :class:`PlantInputs`.

    Expected columns (header row, comma-separated)::

        t, Cwo, Qtp, M_root, M_stem, M_leaf, M_grain

    Mapping from a HYDRUS-1D / Phydrus run (Method A):
        * ``t``   -> output times [day]
        * ``Cwo`` -> dissolved solute concentration at the root-zone node
                     (HYDRUS ``Conc`` in ``Obs_Node.out`` / ``solute1.out``)
        * ``Qtp`` -> actual root water uptake / transpiration
                     (HYDRUS ``vRoot`` or ``T_pot``/``T_act`` in ``T_Level.out``)
        * ``M_*`` -> tissue fresh masses from the growth sub-model [kg]
    """
    data = np.genfromtxt(path, delimiter=",", names=True)
    cols = data.dtype.names
    required = ("t", "Cwo", "Qtp", "M_root", "M_stem", "M_leaf", "M_grain")
    missing = [c for c in required if c not in cols]
    if missing:
        raise ValueError(f"CSV missing required columns: {missing}; found {cols}")
    t = np.atleast_1d(data["t"]).astype(float)
    M = np.column_stack([data["M_root"], data["M_stem"],
                         data["M_leaf"], data["M_grain"]]).astype(float)
    return PlantInputs(t=t, Cwo=np.atleast_1d(data["Cwo"]).astype(float),
                       Qtp=np.atleast_1d(data["Qtp"]).astype(float), M=M)


# ---------------------------------------------------------------------------
# library values (illustrative, NOT calibrated)
# ---------------------------------------------------------------------------
def example_paddy_redox(K_F_drained=2.0, K_F_flooded=1.0, n=0.85,
                        theta_g=0.45) -> PaddyRedox:
    """A plausible-but-illustrative paddy soil: non-linear sorption that
    weakens on flooding (anaerobic -> more bioavailable PFAS)."""
    return PaddyRedox(
        drained=FreundlichSoil(K_F=K_F_drained, n=n, theta_g=theta_g, name="drained/aerobic"),
        flooded=FreundlichSoil(K_F=K_F_flooded, n=n, theta_g=theta_g, name="flooded/anaerobic"),
    )


def _demo():
    """Soil -> plant demo: a flooding schedule modulates bioavailability."""
    from pfas_rice_plant_module import (
        Environment, Compound, Compartment, RiceUptakeModel, binding_factors,
        _logistic, ROOT, STEM, LEAF, FRUIT,
    )
    season = 120.0
    t = np.linspace(0.0, season, 481)
    # constant total soil inventory; paddy flooded early, drained for harvest
    C_total = np.full_like(t, 5.0)                 # ug/kg dry soil
    flooded = t < 90.0                             # mid-season drainage at day 90
    Qtp = 0.05 + 0.35 * np.exp(-((t - 75.0) ** 2) / (2 * 25.0 ** 2))
    M = np.column_stack([
        _logistic(t, 1e-3, 0.030, 0.10, 20.0),
        _logistic(t, 1e-3, 0.040, 0.10, 25.0),
        _logistic(t, 1e-3, 0.050, 0.12, 30.0),
        _logistic(t, 1e-5, 0.025, 0.18, 80.0),
    ])
    redox = example_paddy_redox()
    inputs = inputs_from_soil(t, C_total, Qtp, M, redox, flooded=flooded)

    cmpd = Compound(name="PFOA", K_prot=50.0, K_PL=100.0, K_cw=20.0, kappa_d=0.5,
                    Vmax_in=20.0, Km_in=5.0, Vmax_out=8.0, Km_out=5.0,
                    L_Ph=0.005, f_xy=0.02)
    comps = [
        Compartment("root",  theta=0.70, f_prot=0.05, f_PL=0.010, f_cw=0.30),
        Compartment("stem",  theta=0.80, f_prot=0.01, f_PL=0.005, f_cw=0.08),
        Compartment("leaf",  theta=0.80, f_prot=0.03, f_PL=0.020, f_cw=0.04, S=20.0),
        Compartment("grain", theta=0.15, f_prot=0.08, f_PL=0.010, f_cw=0.10, S=2.0),
    ]
    model = RiceUptakeModel(env=Environment(), cmpd=cmpd, comps=comps, inputs=inputs)
    sol = model.solve(t)

    print(f"paddy soil (Freundlich n={redox.drained.n}): "
          f"K_F drained={redox.drained.K_F}, flooded={redox.flooded.K_F}")
    print(f"pore-water C_w^o: flooded={inputs.Cwo_(45):.4f}  drained={inputs.Cwo_(110):.4f} ug/L "
          f"(total soil = {C_total[0]} ug/kg dry)")
    Cend = sol.y[:, -1]
    Mf = inputs.M_(t[-1])
    straw = (Cend[STEM] * Mf[STEM] + Cend[LEAF] * Mf[LEAF]) / (Mf[STEM] + Mf[LEAF])
    print("final tissue conc [ug/kg]: "
          f"root={Cend[ROOT]:.3f}  straw={straw:.3f}  grain={Cend[FRUIT]:.3f}")


if __name__ == "__main__":
    _demo()
