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
              common factor on DVRJ/DVRI/DVRP puts DVS 1 (flowering) on the measured heading
              day, DVRR puts DVS 2 on the harvest day.
  water flow  ORYZA2000's own potential transpiration (manual eq. 4.4) from the FAO-56
              Penman-Monteith radiation/aerodynamic terms of the same weather and the port's own
              LAI (src/fao56.py); a flooded paddy is not water-limited, so actual = potential.
  mass basis  dry matter -> fresh organ mass through organ_mass (the ODE runs on fresh mass).

Yamazaki harvested two cultivars, Kasalath (heading 2020-08-03, harvest 2020-09-01) and
Koshihikari (heading 2020-08-11, harvest 2020-09-10), and reports the geometric mean of the two;
the model is run once per cultivar (same crop file, each with its own measured heading and
harvest) and the two runs are combined the same way.

Inputs the paper does not report take sourced values (docs/literature_db/sourced_inputs.csv):
nursery 25-26 d after emergence (sown 04-13, transplanted 05-11, emergence 2-3 d after sowing
by the Ibaraki calendar; midpoint 25.5: the drivers of the 25- and 26-d runs interpolated, as
ORYZA counts whole days), seedbed cover 0-5.0 C on the daily maximum
(midpoint 2.5), planting density 18.5 hills m-2 with 4-5 plants per hill (midpoint 4.5; the
Ibaraki prefecture standard), ambient CO2 414 ppm (2020 Mauna Loa annual mean); harvest =
physiological maturity (DVS 2). The ranges are carried into the output (model_api.RANGE_INPUTS).
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
# heading of the cultivar harvested on the same index (Kasalath, Koshihikari; SI p.S3)
YAMAZAKI_HEADINGS = (_dt.date(2020, 8, 3), _dt.date(2020, 8, 11))
YAMAZAKI_SOWING = _dt.date(2020, 4, 13)  # nursery boxes (SI p.S2)

# sourced values and ranges (docs/literature_db/sourced_inputs.csv; PREREG amendment 1)
# seedbed days after emergence: sown 04-13, transplanted 05-11 (28 d), emergence 2-3 d after
# sowing (Ibaraki calendar) -> 25-26 d; midpoint. ORYZA counts whole days, so a fractional
# SBDUR runs the two whole-day neighbours and interpolates their drivers (`plant_run`).
SBDUR_RANGE = (25, 26)
SBDUR = 0.5 * (SBDUR_RANGE[0] + SBDUR_RANGE[1])
# seedbed cover: ORYZA's TMPSB raises the daily MAXIMUM. The Ibaraki calendar keeps the
# seedbed at 15-25 C (greening) and 10-25 C (hardening); the ambient maximum over the seedbed
# days averages 20.0 C (JMA Tsukuba 2020), already inside the band, so the cover can add
# 0 (no warming needed) up to 25 - 20.0 = 5.0 C (the band's top). Midpoint.
TMPSB_RANGE = (0.0, 5.0)
TMPSB = 0.5 * (TMPSB_RANGE[0] + TMPSB_RANGE[1])
# planting: Ibaraki prefecture standard for Koshihikari (18 cm hill spacing, 4-5 seedlings);
# the NARO field's own spacing is unreported
NH = 18.5
NPLH_RANGE = (4.0, 5.0)
NPLH = 0.5 * (NPLH_RANGE[0] + NPLH_RANGE[1])
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
YAMAZAKI_HEADING_MEAN_DOY = sum(_doy(h) for h in YAMAZAKI_HEADINGS) / 2.0


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


def _mgmt(sbdur, tmpsb=TMPSB, nh=NH, nplh=NPLH, co2=CO2_2020):
    if float(sbdur) != int(sbdur):
        raise ValueError("ORYZA's seedbed duration is whole days (plant_run interpolates)")
    sbdur = int(sbdur)
    tp = _doy(YAMAZAKI_TRANSPLANT)
    return O2K.Management(lat=TSUKUBA_LAT, emd=tp - sbdur, sbdur=sbdur, nplh=nplh, nh=nh, co2=co2,
                          tmpsb=tmpsb)


def _run(w, rates, mgmt, crop_set="IR72"):
    import crop_sets as CS
    crop = O2K.load_crop()
    crop.update(CS.overrides(crop_set))      # a published japonica set (src/crop_sets.py)
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


def calibrate_rates(w, heading_doy, harvest_doy, mgmt, crop_set="IR72"):
    """Common factor on DVRJ/DVRI/DVRP -> DVS 1 on heading_doy; DVRR -> DVS 2 on harvest_doy
    (starting from the crop set's own rates)."""
    import crop_sets as CS
    base = {k: float(O2K.load_crop()[k]) for k in RATE_KEYS}
    base.update({k: v for k, v in CS.overrides(crop_set).items() if k in RATE_KEYS})

    def solve(fun, lo, hi):
        for _ in range(6):
            if fun(lo) * fun(hi) < 0:
                return brentq(fun, lo, hi, xtol=1e-8, rtol=1e-8)
            lo, hi = lo / 2.0, hi * 2.0
        raise RuntimeError("no bracket")

    def miss(rates, lim, day):
        x = _crossing(_run(w, rates, mgmt, crop_set), lim)
        return (x if x is not None else 1.0e4) - day

    f = solve(lambda x: miss({**base, **{k: base[k] * x for k in ("DVRJ", "DVRI", "DVRP")}}, 1.0,
                             heading_doy), 0.3, 3.0)
    rates = {**base, **{k: base[k] * f for k in ("DVRJ", "DVRI", "DVRP")}}
    rates["DVRR"] = solve(lambda x: miss({**rates, "DVRR": x}, 2.0, harvest_doy),
                          base["DVRR"] / 3, base["DVRR"] * 3)
    return rates, f


def heading_for(harvest):
    """The measured heading day of year for a harvest: the cultivar's own (Yamazaki SI p.S3);
    the mean-date run takes the mean of the two."""
    hv = _doy(harvest) if not isinstance(harvest, (int, float)) else float(harvest)
    for h, hd in zip(YAMAZAKI_HARVESTS, YAMAZAKI_HEADINGS):
        if abs(_doy(h) - hv) < 1e-9:
            return float(_doy(hd))
    if abs(hv - YAMAZAKI_HARVEST_MEAN_DOY) < 1e-9:
        return YAMAZAKI_HEADING_MEAN_DOY
    raise ValueError(f"no measured heading for harvest day {hv}; pass heading_doy")


def plant_run(harvest, w=None, sbdur=SBDUR, tmpsb=TMPSB, nh=NH, nplh=NPLH, co2=CO2_2020,
              heading_doy=None, harvest_is_maturity=True, crop_set="IR72"):
    """One cultivar's season: daily main-field states from transplanting to harvest.

    Returns t (days after transplanting), dry organ mass [kg/hill], LAI, leaf death rate, and the
    water terms (ET0, ETRD, ETAE, Tp, Ep [mm/d]; Qtp [L/d/hill]).

    ORYZA counts the seedbed in whole days (`Management.sbdur` is an integer). A fractional
    `sbdur` runs the two whole-day neighbours and returns their series linearly interpolated
    (same time grid: transplanting to harvest), the driver at the midpoint of a sourced range."""
    w = w or load_tsukuba_2020()
    lo = int(np.floor(sbdur))
    if abs(sbdur - lo) > 1e-12:
        frac = float(sbdur - lo)
        kw = dict(w=w, tmpsb=tmpsb, nh=nh, nplh=nplh, co2=co2, heading_doy=heading_doy,
                  harvest_is_maturity=harvest_is_maturity, crop_set=crop_set)
        a, b = plant_run(harvest, sbdur=lo, **kw), plant_run(harvest, sbdur=lo + 1, **kw)
        out = {k: ((1.0 - frac) * a[k] + frac * b[k]
                   if isinstance(a[k], np.ndarray) and k != "t" else a[k]) for k in a}
        out["cbc_max"] = max(abs(a["cbc_max"]), abs(b["cbc_max"]))
        out["sbdur_runs"] = (lo, lo + 1, frac)
        return out
    sbdur = lo
    mgmt = _mgmt(sbdur, tmpsb, nh, nplh, co2)
    tp, hv = _doy(YAMAZAKI_TRANSPLANT), _doy(harvest)
    heading = heading_for(harvest) if heading_doy is None else float(heading_doy)
    if harvest_is_maturity:
        rates, f = calibrate_rates(w, heading, hv, mgmt, crop_set)
    else:                                   # sensitivity: IR72 grain-filling rate kept
        base = {k: float(O2K.load_crop()[k]) for k in RATE_KEYS}
        rates, f = calibrate_rates(w, heading, hv, mgmt, crop_set)
        rates["DVRR"] = base["DVRR"]
    r = _run(w, rates, mgmt, crop_set)
    sta = r["cropsta"]
    i3 = int(np.where(sta == 3)[0][0])
    main = np.where(sta >= 4)[0]
    pltr = mgmt.nplh * mgmt.nh / mgmt.nplsb
    td = np.concatenate([[0.0], r["doy"][main] - tp]).astype(float)
    t = np.arange(0.0, np.floor(hv - tp) + 1.0)
    if t[-1] < hv - tp - 1e-9:                         # fractional harvest day: end exactly there
        t = np.append(t, hv - tp)
    ha_per_hill = 1.0 / (nh * 1.0e4)
    out = dict(t=t, rates=rates, factor=f, heading_doy=heading, harvest_doy=hv, crop_set=crop_set)
    for name, key in (("root", "wrt"), ("stem", "wst"), ("leaf", "wlvg"), ("grain", "wso")):
        v = np.concatenate([[r[key][i3] * pltr], r[key][main]])
        out[name] = np.interp(t, td, v * ha_per_hill)     # crop model's own value; held after maturity
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


def drivers_from_run(run, Cwo=1.0, organ_water="constant", grain_basis="panicle",
                     theta_anchors=None, panicle_kw=None):
    """A `model_api.simulate(drivers=...)` dict from one plant run (fresh organ mass).

    organ_water "stage": root, stem and leaf water follow the run's DVS (PREREG 7.1) and the dict
    carries the theta(t) matrix that made M fresh ("theta"); "constant": the composition table's.
    grain_basis "brown": the grain column is brown rice and the panicle remainder's DRY mass is in
    "husk" (PREREG 7.2, `panicle.split_biomass`); "panicle": the grain is ORYZA's WSO."""
    import panicle as pn
    b = {k: run[k] for k in ("root", "stem", "leaf", "grain")}
    b["dvs"] = np.asarray(run["dvs"], float)
    if grain_basis == "brown":
        b = pn.split_biomass(b, run["t"], **(panicle_kw or {}))
    elif grain_basis != "panicle":
        raise ValueError(f"grain_basis must be 'brown' or 'panicle', got {grain_basis!r}")
    if organ_water not in ("stage", "constant"):
        raise ValueError(f"organ_water must be 'stage' or 'constant', got {organ_water!r}")
    stage = organ_water == "stage"
    d = dict(t=run["t"], Cwo=np.full_like(run["t"], float(Cwo)), Qtp=run["Qtp"],
             M=om.organ_matrix(b, stage=stage, anchors=theta_anchors),
             leaf_loss=run["leaf_death_rate"], mass_basis=om.MASS_BASIS, dvs=b["dvs"])
    if stage:
        d["theta"] = om.theta_matrix(b, len(run["t"]), stage=True, anchors=theta_anchors)
    if grain_basis == "brown":
        d["husk"] = np.asarray(b["husk"], float)
        d["panicle"] = np.asarray(b["panicle"], float)
    return d


def yamazaki2020_runs(**kw):
    """The two plant runs (harvest 2020-09-01 and 2020-09-10)."""
    w = load_tsukuba_2020()
    return [plant_run(h, w=w, **kw) for h in YAMAZAKI_HARVESTS]


def yamazaki2020_mean_run(**kw):
    """The default single run: harvest at the midpoint of the two harvest dates."""
    return plant_run(YAMAZAKI_HARVEST_MEAN_DOY, **kw)


def check_season(season, t, where="simulate"):
    """The site scenario sets its own season (transplanting to harvest, ``t[-1]``). An
    explicit ``season`` that differs is refused rather than ignored: ignoring it would run
    the site-year while the caller believes it ran a different season length."""
    site = float(t[-1])
    if season is not None and abs(float(season) - site) > 1e-6:
        raise ValueError(f"{where}: scenario 'yamazaki2020' runs transplanting to harvest "
                         f"({site:g} d), so season={float(season):g} does not apply; pass "
                         "scenario='generic' to choose a season length")
