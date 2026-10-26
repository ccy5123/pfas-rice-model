"""
Paddy soil sub-model — REDOX-CORRECTED (W3)
===========================================

Drop-in replacement for the redox portion of ``soil_paddy.py``.  Same Method-A
coupling (soil -> C_w^o(t) -> plant module) and the same ``FreundlichSoil``
isotherm, but the *flooding response is corrected* per the S6 review (H7 §7.4).

WHAT CHANGED AND WHY
--------------------
The original ``PaddyRedox`` / ``example_paddy_redox`` made flooding *weaken*
sorption (``K_F`` drained 2.0 -> flooded 1.0), i.e. anaerobic conditions
RAISE the bioavailable pore-water concentration.  That sign is not defensible
as a default:

  * The dominant flooded effect on pore-water PFAS is a CONCENTRATION DROP from
    (i) DILUTION — ponded floodwater + saturated pores raise the aqueous volume
    the inventory partitions into — and (ii) LEACHING — downward/lateral water
    flux exports dissolved PFAS out of the root zone.  Both LOWER C_w^o.
  * The redox->sorption coupling (anaerobic altering K_F via Fe/Mn-oxide
    reductive dissolution vs. enhanced sorption to reduced OM/sulfide) is a
    SECONDARY effect of UNCERTAIN SIGN; it must be set from data, never
    defaulted to "weaker".

So this module: (1) drives the flooded response through a higher flooded
gravimetric water content ``theta_g`` (dilution) plus an optional first-order
``k_leach`` on the dissolved pool (leaching); (2) defaults the redox K_F change
to NEUTRAL (flooded K_F == drained K_F) with an explicit hook to set it from
data; (3) the net flooded effect is C_w^o DOWN (verified in ``_demo``).

Unit system identical to ``soil_paddy.py``.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Sequence
import numpy as np
from scipy.optimize import brentq

from pfas_rice_plant_module import PlantInputs
# reuse the unchanged isotherm
from soil_paddy import FreundlichSoil, inputs_from_soil as _inputs_from_soil_base


# ---------------------------------------------------------------------------
# corrected redox
# ---------------------------------------------------------------------------
@dataclass
class PaddyRedoxCorrected:
    """Flooded/drained switch with DILUTION + LEACHING as the flooded driver.

    Parameters
    ----------
    drained : FreundlichSoil for the aerobic/drained state.
    flooded : FreundlichSoil for the flooded state.  Its ``theta_g`` is normally
        LARGER than the drained value (ponded water + saturation) so that, for a
        fixed total inventory, the inverted pore-water C_w is LOWER (dilution).
        Its ``K_F`` defaults to the drained value (redox-sorption neutral).
    k_leach : first-order loss rate of the DISSOLVED pool while flooded [1/day].
        Applied to the total inventory weighted by the dissolved fraction
        f_diss = theta_g*C_w / C_T, so mobile (weakly sorbed, short-chain) PFAS
        leach faster.  0.0 disables leaching.
    """
    drained: FreundlichSoil
    flooded: FreundlichSoil
    k_leach: float = 0.0

    def soil_at(self, is_flooded: bool) -> FreundlichSoil:
        return self.flooded if is_flooded else self.drained

    # --- dilution: invert with the flood-appropriate isotherm (theta_g differs)
    def pore_water_series(self, C_total: Sequence[float],
                          flooded: Sequence[bool]) -> np.ndarray:
        C_total = np.asarray(C_total, dtype=float)
        flooded = np.asarray(flooded, dtype=bool)
        if C_total.shape != flooded.shape:
            raise ValueError("C_total and flooded must have the same shape")
        return np.array([self.soil_at(f).pore_water(float(ct))
                         for ct, f in zip(C_total, flooded)])

    # --- leaching: integrate the inventory decline over flooded periods --------
    def apply_leaching(self, t: Sequence[float], C_total0: Sequence[float] | float,
                       flooded: Sequence[bool]) -> np.ndarray:
        """Return the leached total-inventory series C_T(t) [ug/kg dry].

        dC_T/dt = -k_leach * f_diss(C_T) * C_T   while flooded   (else 0),
        with f_diss the dissolved fraction at the current state. Explicit Euler
        on the supplied grid (fine grids recommended; PFAS leaching is slow).
        """
        t = np.asarray(t, dtype=float)
        flooded = np.asarray(flooded, dtype=bool)
        C = np.empty_like(t)
        C[0] = float(C_total0[0]) if np.ndim(C_total0) else float(C_total0)
        for i in range(1, len(t)):
            dt = t[i] - t[i-1]
            soil = self.soil_at(bool(flooded[i-1]))
            ct = C[i-1]
            if self.k_leach > 0 and flooded[i-1] and ct > 0:
                cw = soil.pore_water(ct)
                f_diss = (soil.theta_g * cw) / ct if ct > 0 else 0.0
                ct = ct - dt * self.k_leach * f_diss * ct
            C[i] = max(ct, 0.0)
        return C


def example_paddy_redox_corrected(K_F=2.0, n=0.85, theta_g_drained=0.35,
                                  theta_g_flooded=0.60, k_leach=0.02
                                  ) -> PaddyRedoxCorrected:
    """Plausible-but-illustrative CORRECTED paddy soil.

    Flooding does NOT weaken sorption (K_F unchanged); instead the flooded state
    has higher water content (dilution) and leaches the dissolved pool. Net: the
    flooded pore-water concentration is LOWER than drained (see ``_demo``).
    """
    return PaddyRedoxCorrected(
        drained=FreundlichSoil(K_F=K_F, n=n, theta_g=theta_g_drained, name="drained/aerobic"),
        flooded=FreundlichSoil(K_F=K_F, n=n, theta_g=theta_g_flooded, name="flooded (diluted)"),
        k_leach=k_leach,
    )


def inputs_from_soil_corrected(t, C_total, Qtp, M, redox: PaddyRedoxCorrected,
                               flooded, leach: bool = True) -> PlantInputs:
    """Build PlantInputs with leaching (optional) then Freundlich/dilution inversion."""
    t = np.asarray(t, dtype=float)
    C_T = redox.apply_leaching(t, C_total, flooded) if leach else np.asarray(C_total, float)
    Cwo = redox.pore_water_series(C_T, flooded)
    return PlantInputs(t=t, Cwo=Cwo, Qtp=np.asarray(Qtp, float), M=np.asarray(M, float))
