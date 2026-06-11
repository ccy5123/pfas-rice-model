"""
Rice transpiration forcing Q_TP(t)  (task 2, part 1 — the xylem-flow driver)
============================================================================

Builds the plant transpiration stream Q_TP(t) [L/day per hill] for the
N-segment uptake model from measured paddy crop-physiology, NOT from PFAS data.

Method (FAO-56 dual crop coefficient — transpiration only)
----------------------------------------------------------
In a flooded paddy a large share of evapotranspiration is EVAPORATION from the
ponded water, which does NOT carry the (non-volatile) PFAS anion up the xylem.
Only the TRANSPIRATION component drives xylem advection, so we use the basal
(transpiration) crop coefficient Kcb:

    T(t) [mm/day] = Kcb(t) * ET0          (FAO-56 dual: ET = (Kcb+Ke)*ET0)
    Q_TP(t) [L/day/hill] = T(t) * area_per_hill[m^2]      (1 mm over 1 m^2 = 1 L)

Calibration anchors (measured, cited):
  * ET0 ~ 3.39 mm/day (mean; stagewise 3.42/3.56/.../2.97) — Kumari et al. 2022,
    transplanted puddled-rice lysimeter (Agronomy 12:2850, doi 10.3390/agronomy12112850).
  * actual stagewise Kc 1.13/1.27/1.23/0.93 (initial/dev/mid/late) — Kumari et al. 2022.
  * transpiration/evaporation split: seasonal soil evaporation = 57.9% of water
    loss (=> T/ET ~ 42%), evaporation fraction 96.6% (bare, initial) -> 43.3%
    (closed canopy) — Nay Htoon et al. 2018, PLoS ONE (FAO-56 dual-Kc partition).
  * Kcb (rice, basal/transpiration coefficient): ini 0.15, mid 1.10, end 0.55
    — FAO-56 (Allen et al. 1998), Table 17; consistent with the Kumari Kc once
    the evaporation coefficient Ke is added (Kc_ini>>Kcb_ini because early ET is
    evaporation-dominated, exactly the Nay Htoon result).

So T(t) rises from a low base (transplant/bare water, transpiration ~ a few % of
ET) to a mid-season peak at canopy closure, then declines at maturity/drainage.
Stage boundaries are fractions of the season (FAO-56 transplanted-rice shape).

This module fixes the *shape and magnitude* of Q_TP(t); the absolute per-hill
value still scales with planting density (area_per_hill) and is paired with the
per-segment biomass M_s(t) (task 2 part 2) for the nstem absolute f_xy fit.
"""
from __future__ import annotations
import numpy as np

# Kumari 2022 lysimeter ACTUAL stagewise crop coefficient Kc (initial/dev/mid/late)
KC_STAGE = (1.13, 1.27, 1.23, 0.93)
# stage boundaries as fraction of season (transplant->harvest): end of
# initial / development / mid-season  (late = remainder).  FAO-56 transplanted rice.
STAGE_FRAC = (0.17, 0.40, 0.75)
ET0_MEAN_MM_D = 3.39                    # Kumari 2022 lysimeter mean ET0-PM [mm/day]
# transpiration fraction f_T = T/ET (Nay Htoon 2018: evaporation 96.6%->43.3% from
# bare to closed canopy => f_T 0.034 (initial) -> 0.567 (full canopy); seasonal ~0.42).
# Flooded-paddy open-water evaporation is large, so f_T must be anchored to this
# MEASURED partition (the FAO-56 default Kcb over-states paddy transpiration ~1.5x).
F_T_INI, F_T_FULL = 0.034, 0.567
# planting density: ~25 hills/m^2 (transplanted paddy) -> 0.04 m^2 per hill
AREA_PER_HILL_M2 = 0.04


def kc_curve(t: np.ndarray, season: float) -> np.ndarray:
    """Stagewise actual Kc(t) (Kumari 2022), piecewise-linear across FAO-56 stages."""
    f = np.clip(np.asarray(t, float) / season, 0.0, 1.0)
    f_ini, f_dev, f_mid = STAGE_FRAC
    xp = [0.0, f_ini, (f_ini + f_dev) / 2, f_dev, (f_dev + f_mid) / 2, f_mid, 1.0]
    yp = [KC_STAGE[0], KC_STAGE[0], (KC_STAGE[0] + KC_STAGE[1]) / 2, KC_STAGE[1],
          KC_STAGE[2], KC_STAGE[2], KC_STAGE[3]]
    return np.interp(f, xp, yp)


def ft_curve(t: np.ndarray, season: float) -> np.ndarray:
    """Transpiration fraction f_T(t)=T/ET (Nay Htoon 2018): rises with canopy cover
    from F_T_INI (bare/flooded) to F_T_FULL (closed canopy), then flat."""
    f = np.clip(np.asarray(t, float) / season, 0.0, 1.0)
    f_ini, f_dev, _ = STAGE_FRAC
    # canopy closes by the end of the development stage (transplanted rice) -> f_T
    # reaches its full-canopy value at f_dev, then stays flat (seasonal T/ET ~ 0.42)
    return np.interp(f, [0.0, f_ini, f_dev, 1.0], [F_T_INI, F_T_INI, F_T_FULL, F_T_FULL])


def transpiration_mm_d(t: np.ndarray, season: float = 120.0,
                       et0: float = ET0_MEAN_MM_D) -> np.ndarray:
    """Canopy transpiration T(t) [mm/day] = f_T(t) * Kc(t) * ET0 (evaporation excluded)."""
    return ft_curve(t, season) * kc_curve(t, season) * et0


def Q_TP(t: np.ndarray, season: float = 120.0, et0: float = ET0_MEAN_MM_D,
         area_per_hill: float = AREA_PER_HILL_M2) -> np.ndarray:
    """Transpiration stream Q_TP(t) [L/day per hill] for the plant ODE."""
    return transpiration_mm_d(t, season, et0) * area_per_hill


def seasonal_T_over_ET(season: float = 120.0, et0: float = ET0_MEAN_MM_D) -> float:
    """Seasonal transpiration fraction T/ET (sanity check vs Nay Htoon ~0.42).
    Uses FAO-56 Kc (=Kcb+Ke) approximated by the Kumari stagewise actual Kc."""
    t = np.linspace(0, season, 1000)
    T = transpiration_mm_d(t, season, et0)
    ET = kc_curve(t, season) * et0
    _trapz = getattr(np, "trapezoid", None) or np.trapz       # numpy 2.x renamed trapz
    return float(_trapz(T, t) / _trapz(ET, t))


if __name__ == "__main__":
    season = 120.0
    t = np.linspace(0, season, 481)
    Q = Q_TP(t, season)
    T = transpiration_mm_d(t, season)
    print(f"FAO-56 dual-Kc rice transpiration forcing (ET0={ET0_MEAN_MM_D} mm/d, "
          f"{int(1/AREA_PER_HILL_M2)} hills/m^2):")
    print(f"  peak transpiration T = {T.max():.2f} mm/day "
          f"(at t={t[T.argmax()]:.0f} d, mid-season canopy closure)")
    print(f"  peak Q_TP = {Q.max():.3f} L/day/hill   (placeholder demo used ~0.40 -> ~{0.40/Q.max():.1f}x high)")
    print(f"  seasonal mean Q_TP = {Q.mean():.3f} L/day/hill")
    print(f"  seasonal T/ET = {seasonal_T_over_ET(season):.2f}  (Nay Htoon 2018: ~0.42)")
    for d in (10, 30, 60, 90, 115):
        a = np.array([d])
        print(f"    day {d:3d}: Kc={float(kc_curve(a,season)[0]):.2f}  f_T={float(ft_curve(a,season)[0]):.2f}  "
              f"T={float(transpiration_mm_d(a,season)[0]):.2f} mm/d  "
              f"Q_TP={float(Q_TP(a,season)[0]):.3f} L/d/hill")
