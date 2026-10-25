"""
Site-year scenario for the PFAS calibration dataset: Yamazaki et al. 2023 (NARO Tsukuba paddy,
2020), driven by the measured weather and the dates the paper reports.

Why: the PFAS transport parameters are fitted to Yamazaki's tissue concentrations, so the model
must be run under Yamazaki's conditions -- not a generic climatology with a fixed 120-d season, a
constant ET0 from another site and a transpiration shape stretched over that season. Here every
driver comes from one source chain:

  weather     JMA Tsukuba 2020 (raw_si/jma_tsukuba_2020.csv): Tmax/Tmin, global radiation, vapour
              pressure, wind (10.6 m -> 2 m by FAO-56 eq. 47).
  biomass     the ORYZA2000 port (src/oryza2000.py, IR72.D92) on that weather, transplanted on the
              paper's date (2020-05-11), development rates calibrated to the paper's dates: one
              common factor on DVRJ/DVRI/DVRP puts DVS 1 (flowering) on the heading day (harvest
              minus 30-31 d, as the paper states), DVRR puts DVS 2 on the harvest day.
  water flow  ORYZA2000's own potential transpiration (manual eq. 4.4) from the FAO-56
              Penman-Monteith radiation/aerodynamic terms of the same weather and the port's own
              LAI (src/fao56.py); a flooded paddy is not water-limited, so actual = potential.
  mass basis  dry matter -> fresh organ mass through organ_mass (the ODE runs on fresh mass).

Yamazaki harvested two cultivars, on 2020-09-01 and 2020-09-10, each 30-31 d after heading, and
reports the geometric mean of the two; the model is run once per harvest date (same crop file,
different phenology) and the two runs are combined the same way.

Inputs the paper does not report, stated (see docs/ASSUMPTIONS.md): nursery length and cover
(SBDUR 23 d after emergence, open seedbed -- the median of the Shizukuishi FACE nurseries, Kim
2003; sensitivity in validation/yamazaki_site_year.py), planting density (19.05 hills m-2,
3 plants/hill -- Japanese machine spacing, as in Kim 2003), ambient CO2 (414 ppm, the 2020 Mauna
Loa annual mean), heading = harvest - 30.5 d (midpoint of the stated 30-31 d), harvest =
physiological maturity (DVS 2).
"""
from __future__ import annotations

import csv
import datetime as _dt
import os

import numpy as np
from scipy.optimize import brentq

import fao56 as F
import oryza2000 as O2K
import organ_mass as om

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "docs", "literature_db", "raw_si")

# JMA Tsukuba (Aerological Observatory) -- ame_master_20260324.csv
TSUKUBA_LAT = 36.0 + 3.4 / 60.0
TSUKUBA_ELEV_M = 25.0
TSUKUBA_ANEMOMETER_M = 10.6

# Yamazaki et al. 2023 SI p.S2-S3
YAMAZAKI_TRANSPLANT = _dt.date(2020, 5, 11)
YAMAZAKI_HARVESTS = (_dt.date(2020, 9, 1), _dt.date(2020, 9, 10))
DAYS_HEADING_TO_HARVEST = 30.5          # "30-31 d after heading": midpoint of the stated range

# not reported by Yamazaki -- stated defaults (docs/ASSUMPTIONS.md)
SBDUR = 23                              # seedbed days after emergence (Kim 2003 median: 14/23/24)
TMPSB = 0.0                             # open seedbed
NH, NPLH = 19.05, 3.0                   # hills m-2, plants hill-1 (30 x 17.5 cm, as Kim 2003)
CO2_2020 = 414.0                        # ppm, Mauna Loa 2020 annual mean
RATE_KEYS = ("DVRJ", "DVRI", "DVRP", "DVRR")


def _doy(d):
    """Day of year of a date; a number is taken as a (possibly fractional) day of year."""
    if isinstance(d, (int, float, np.floating, np.integer)):
        return float(d)
    return d.timetuple().tm_yday


# the default single run: harvest at the midpoint of the two harvest dates (fractional day).
# The calibration itself runs both harvests and combines them as the data are combined.
YAMAZAKI_HARVEST_MEAN_DOY = sum(_doy(h) for h in YAMAZAKI_HARVESTS) / 2.0


def load_tsukuba_2020():
    """Daily weather (doy0 + arrays) in the port's format, plus ea [kPa] and u2 [m/s]."""
    with open(os.path.join(RAW, "jma_tsukuba_2020.csv")) as f:
        rows = list(csv.DictReader(line for line in f if not line.startswith("#")))
    rows.sort(key=lambda r: int(r["doy"]))
    doys = np.array([int(r["doy"]) for r in rows])
    assert np.all(np.diff(doys) == 1), "weather series must be contiguous"
    g = lambda k: np.array([float(r[k]) for r in rows])
    return dict(doy0=int(doys[0]), rad_mj=g("rad_mj"), tmin=g("tmin"), tmax=g("tmax"),
                ea_kpa=g("ea_kpa"), u2=F.wind_2m(g("wind_10p6m"), TSUKUBA_ANEMOMETER_M))


def _mgmt(sbdur=SBDUR, tmpsb=TMPSB, nh=NH, nplh=NPLH, co2=CO2_2020):
    tp = _doy(YAMAZAKI_TRANSPLANT)
    return O2K.Management(lat=TSUKUBA_LAT, emd=tp - sbdur, sbdur=sbdur, nplh=nplh, nh=nh, co2=co2,
                          tmpsb=tmpsb)


def _run(w, rates, mgmt):
    crop = O2K.load_crop()
    crop.update(rates)
    return O2K.simulate(w, crop, mgmt, max_days=len(w["rad_mj"]) - (mgmt.emd - w["doy0"]))


def _crossing(r, lim):
    """Fractional day of year at which DVS first reaches `lim` (None if never)."""
    d, t = r["dvs"], r["doy"].astype(float)
    hit = np.where(d >= lim - 1e-12)[0]
    if not len(hit):
        return None
    i = int(hit[0])
    if i == 0:
        return float(t[0])
    return float(t[i - 1] + (lim - d[i - 1]) / (d[i] - d[i - 1]) * (t[i] - t[i - 1]))


def calibrate_rates(w, heading_doy, harvest_doy, mgmt):
    """Common factor on DVRJ/DVRI/DVRP -> DVS 1 on heading_doy; DVRR -> DVS 2 on harvest_doy."""
    base = {k: float(O2K.load_crop()[k]) for k in RATE_KEYS}

    def solve(fun, lo, hi):
        for _ in range(6):
            if fun(lo) * fun(hi) < 0:
                return brentq(fun, lo, hi, xtol=1e-8, rtol=1e-8)
            lo, hi = lo / 2.0, hi * 2.0
        raise RuntimeError("no bracket")

    def miss(rates, lim, day):
        x = _crossing(_run(w, rates, mgmt), lim)
        return (x if x is not None else 1.0e4) - day

    f = solve(lambda x: miss({**base, **{k: base[k] * x for k in ("DVRJ", "DVRI", "DVRP")}}, 1.0,
                             heading_doy), 0.3, 3.0)
    rates = {**base, **{k: base[k] * f for k in ("DVRJ", "DVRI", "DVRP")}}
    rates["DVRR"] = solve(lambda x: miss({**rates, "DVRR": x}, 2.0, harvest_doy),
                          base["DVRR"] / 3, base["DVRR"] * 3)
    return rates, f


def plant_run(harvest, w=None, sbdur=SBDUR, tmpsb=TMPSB, nh=NH, nplh=NPLH, co2=CO2_2020,
              days_heading_to_harvest=DAYS_HEADING_TO_HARVEST, harvest_is_maturity=True):
    """One cultivar's season: daily main-field states from transplanting to harvest.

    Returns t (days after transplanting), dry organ mass [kg/hill], LAI, leaf death rate, and the
    water terms (ET0, ETRD, ETAE, Tp, Ep [mm/d]; Qtp [L/d/hill])."""
    w = w or load_tsukuba_2020()
    mgmt = _mgmt(sbdur, tmpsb, nh, nplh, co2)
    tp, hv = _doy(YAMAZAKI_TRANSPLANT), _doy(harvest)
    heading = hv - days_heading_to_harvest
    if harvest_is_maturity:
        rates, f = calibrate_rates(w, heading, hv, mgmt)
    else:                                   # sensitivity: IR72 grain-filling rate kept
        base = {k: float(O2K.load_crop()[k]) for k in RATE_KEYS}
        rates, f = calibrate_rates(w, heading, hv, mgmt)
        rates["DVRR"] = base["DVRR"]
    r = _run(w, rates, mgmt)
    sta = r["cropsta"]
    i3 = int(np.where(sta == 3)[0][0])
    main = np.where(sta >= 4)[0]
    pltr = mgmt.nplh * mgmt.nh / mgmt.nplsb
    td = np.concatenate([[0.0], r["doy"][main] - tp]).astype(float)
    t = np.arange(0.0, np.floor(hv - tp) + 1.0)
    if t[-1] < hv - tp - 1e-9:                         # fractional harvest day: end exactly there
        t = np.append(t, hv - tp)
    ha_per_hill = 1.0 / (nh * 1.0e4)
    out = dict(t=t, rates=rates, factor=f, heading_doy=heading, harvest_doy=hv)
    for name, key in (("root", "wrt"), ("stem", "wst"), ("leaf", "wlvg"), ("grain", "wso")):
        v = np.concatenate([[r[key][i3] * pltr], r[key][main]])
        out[name] = np.maximum(np.interp(t, td, v * ha_per_hill), 1e-9)      # held after maturity
    lai = np.interp(t, td, np.concatenate([[r["lai"][i3] * pltr], r["lai"][main]]))
    drlv = np.concatenate([[r["drlv"][i3]], r["drlv"][main]])
    out["leaf_death_rate"] = np.where(t > td[-1], 0.0, np.maximum(np.interp(t, td, drlv), 0.0))
    out["dvs"] = np.interp(t, td, np.concatenate([[r["dvs"][i3]], r["dvs"][main]]))
    out["lai"] = lai
    k = int(tp - w["doy0"]) + np.floor(t).astype(int)
    doy = tp + t
    alb = F.field_reflection(lai, flooded=True)
    et0, etrd, etae = F.et0_terms(w["tmax"][k], w["tmin"][k], w["rad_mj"][k], w["ea_kpa"][k],
                                  w["u2"][k], doy, TSUKUBA_LAT, TSUKUBA_ELEV_M, albedo=alb)
    out["ET0_ref"] = F.et0_terms(w["tmax"][k], w["tmin"][k], w["rad_mj"][k], w["ea_kpa"][k],
                                 w["u2"][k], doy, TSUKUBA_LAT, TSUKUBA_ELEV_M)[0]
    out["Tp"] = F.oryza_transpiration(etrd, etae, lai)
    out["Ep"] = F.oryza_evaporation(etrd, etae, lai)
    out["Qtp"] = out["Tp"] * (1.0 / nh)                  # mm d-1 x m2 per hill = L d-1 hill-1
    out["cbc_max"] = r["cbc_max"]
    return out


def drivers_from_run(run, Cwo=1.0):
    """A `model_api.simulate(drivers=...)` dict from one plant run (fresh organ mass)."""
    b = {k: run[k] for k in ("root", "stem", "leaf", "grain")}
    return dict(t=run["t"], Cwo=np.full_like(run["t"], float(Cwo)), Qtp=run["Qtp"],
                M=om.organ_matrix(b), leaf_loss=run["leaf_death_rate"],
                mass_basis=om.MASS_BASIS)


def yamazaki2020_runs(**kw):
    """The two plant runs (harvest 2020-09-01 and 2020-09-10)."""
    w = load_tsukuba_2020()
    return [plant_run(h, w=w, **kw) for h in YAMAZAKI_HARVESTS]


def yamazaki2020_mean_run(**kw):
    """The default single run: harvest at the midpoint of the two harvest dates."""
    return plant_run(YAMAZAKI_HARVEST_MEAN_DOY, **kw)
