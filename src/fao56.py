"""
FAO-56 Penman-Monteith reference evapotranspiration, split into its radiation and
aerodynamic terms, and ORYZA2000's potential transpiration / soil-water evaporation.

Allen RG, Pereira LS, Raes D, Smith M (1998) Crop evapotranspiration. FAO Irrigation and
Drainage Paper 56 (equations cited as "FAO-56 eq. N"). The ORYZA2000 manual (Bouman et al.
2001, Section 4.1) defines its reference rate ETD by the FAO definition and computes

    Tp   = ETRD * (1 - exp(-0.5 LAI)) + ETAE * min(2, LAI)          (manual eq. 4.4)
    Ep   = (ETRD + ETAE) * exp(-0.5 LAI)                              (manual eq. 4.3)
    RF   = RFS * exp(-0.5 LAI) + 0.25 * (1 - exp(-0.5 LAI))           (manual eq. 4.2)
    RFS  = 0.05 under standing water > 5 mm, else 0.25*(1 - 0.5 theta/theta_s)

with ETRD / ETAE the radiation and drying-power terms of the Penman rate. The manual does not
print its Penman routine (SETPMD; "given by van Kraalingen and Stol (1997) and not repeated
here"), so the two terms here are the radiation and aerodynamic terms of the FAO-56
Penman-Monteith equation (FAO-56 eq. 6), with the field reflection RF of manual eq. 4.2 in
the net shortwave radiation (RF is the input the manual computes for its Penman call).
Units: temperatures C, radiation MJ m-2 d-1, vapour pressure kPa, wind m s-1, ET mm d-1.
"""
from __future__ import annotations

import numpy as np

SIGMA = 4.903e-9            # Stefan-Boltzmann [MJ K-4 m-2 d-1]
GSC = 0.0820                # solar constant [MJ m-2 min-1]


def psychrometric_constant(elev_m):
    """FAO-56 eq. 7-8 [kPa C-1]."""
    p = 101.3 * ((293.0 - 0.0065 * elev_m) / 293.0) ** 5.26
    return 0.665e-3 * p


def es(t):
    """Saturation vapour pressure, FAO-56 eq. 11 [kPa]."""
    return 0.6108 * np.exp(17.27 * t / (t + 237.3))


def slope_svp(t):
    """Slope of the saturation vapour pressure curve, FAO-56 eq. 13 [kPa C-1]."""
    return 4098.0 * es(t) / (t + 237.3) ** 2


def extraterrestrial_radiation(doy, lat_deg):
    """Ra, FAO-56 eq. 21-25 [MJ m-2 d-1]."""
    phi = np.deg2rad(lat_deg)
    j = np.asarray(doy, float)
    dr = 1.0 + 0.033 * np.cos(2 * np.pi * j / 365.0)
    delta = 0.409 * np.sin(2 * np.pi * j / 365.0 - 1.39)
    ws = np.arccos(np.clip(-np.tan(phi) * np.tan(delta), -1.0, 1.0))
    return 24 * 60 / np.pi * GSC * dr * (ws * np.sin(phi) * np.sin(delta)
                                         + np.cos(phi) * np.cos(delta) * np.sin(ws))


def wind_2m(uz, z):
    """Wind speed at 2 m from a measurement at height z [m], FAO-56 eq. 47."""
    return uz * 4.87 / np.log(67.8 * z - 5.42)


def net_radiation(rs, tmax, tmin, ea, doy, lat_deg, elev_m, albedo=0.23):
    """Rn = Rns - Rnl, FAO-56 eq. 37-40 [MJ m-2 d-1]; `albedo` may be an array."""
    ra = extraterrestrial_radiation(doy, lat_deg)
    rso = (0.75 + 2e-5 * elev_m) * ra
    rns = (1.0 - albedo) * rs
    rel = np.clip(rs / rso, 0.0, 1.0)
    rnl = SIGMA * ((tmax + 273.16) ** 4 + (tmin + 273.16) ** 4) / 2.0 \
        * (0.34 - 0.14 * np.sqrt(ea)) * (1.35 * rel - 0.35)
    return rns - rnl


def et0_terms(tmax, tmin, rs, ea, u2, doy, lat_deg, elev_m, albedo=0.23):
    """FAO-56 eq. 6 split into (ET0, radiation term, aerodynamic term) [mm d-1]; G = 0 (daily)."""
    tmean = (tmax + tmin) / 2.0
    d = slope_svp(tmean)
    g = psychrometric_constant(elev_m)
    rn = net_radiation(rs, tmax, tmin, ea, doy, lat_deg, elev_m, albedo)
    den = d + g * (1.0 + 0.34 * u2)
    rad = 0.408 * d * rn / den
    aero = g * 900.0 / (tmean + 273.0) * u2 * ((es(tmax) + es(tmin)) / 2.0 - ea) / den
    return rad + aero, rad, aero


def field_reflection(lai, flooded=True, rfs_soil=0.25):
    """ORYZA2000 manual eq. 4.2: RF of the crop over standing water (RFS 0.05) or soil."""
    rfs = 0.05 if flooded else rfs_soil
    k = np.exp(-0.5 * np.asarray(lai, float))
    return rfs * k + 0.25 * (1.0 - k)


def oryza_transpiration(etrd, etae, lai):
    """Potential transpiration Tp, ORYZA2000 manual eq. 4.4 [mm d-1]."""
    lai = np.asarray(lai, float)
    return etrd * (1.0 - np.exp(-0.5 * lai)) + etae * np.minimum(2.0, lai)


def oryza_evaporation(etrd, etae, lai):
    """Potential soil/water evaporation Ep, ORYZA2000 manual eq. 4.3 [mm d-1]."""
    return (etrd + etae) * np.exp(-0.5 * np.asarray(lai, float))
