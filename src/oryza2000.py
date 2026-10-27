"""
oryza2000.py — ORYZA2000 potential production, ported from the published Fortran
=================================================================================

A line-by-line Python port of the ORYZA2000 potential-production crop model as printed in

    Bouman BAM, Kropff MJ, Tuong TP, Wopereis MCS, ten Berge HFM, van Laar HH (2001).
    ORYZA2000: modeling lowland rice. IRRI, Los Banos. Chapter 3 (ORYZA1 and its
    subroutines), Sections 2.2 (MODELS: crop-stage logic and call order), 4.3 (WNOSTRESS)
    and 5.3 (NNOSTRESS).  http://books.irri.org/9712201716_content.pdf

with the cultivar file IR72.D92 transcribed from Section 7.3
(docs/literature_db/raw_si/oryza_ir72_d92_params.csv).  Routine and variable names follow
the Fortran so each block can be checked against the listing it came from.

Ported: PHENOL / SUBDD / SUBCD (phenology with the transplanting shock), SGPCDT / SASTRO /
SSKYC / SGPC1 / SGPL / SRDPRF / GPPARGET (canopy photosynthesis with sunlit/shaded leaves,
the canopy leaf-N profile and the CO2 effects), maintenance and growth respiration, stem
reserves (FSTR / LSTR / LRSTR), partitioning with the grain sink limit, SUBGRN (spikelets
and their cold/heat sterility), SUBLAI2 (exponential -> SLA leaf area growth, seedbed,
dilution at transplanting, shock), leaf death, the transplanting dilution of every organ,
and the carbon-balance check (SUBCBC).

Not ported, on purpose: everything the potential-production mode switches off.  WNOSTRESS
sets LRSTRS = LDSTRS = LESTRS = PCEW = 1 and NNOSTRESS sets NSLLV = RNSTRS = 1, so the
drought (ET, WSTRESS, PADDY, IRRIG) and nitrogen (NCROP, NSOIL) modules never enter a
potential run; their effect factors are kept in the code at their potential value so the
equations read as printed.  SGPC2 (IACC = 2) is not ported: the standard version uses
SGPC1.  Root depth (ZRT) is not tracked: it feeds only the soil-water balance.

One reading had to be chosen.  GPPARGET prints  ``Amax = 9.5+(22.*SLNI)*cREDFT*AmaxCO2``
for SLNI >= 0.5 and ``68.33*(SLNI-0.2)*cREDFT*AmaxCO2`` below it.  The two branches meet at
SLNI = 0.5 only if the temperature and CO2 factors multiply the whole first branch, i.e.
``(9.5+22*SLNI)*REDFT*AmaxCO2``; that continuous form is the default
(``amax_form="continuous"``) and ``"as_printed"`` keeps the literal text.  They differ only
when REDFT*AmaxCO2 != 1, i.e. below 20 oC daytime temperature or away from 340 ppm.

States are in kg dry matter per hectare.  ``simulate`` returns the daily record ORYZA writes
in its rate section (state at the START of each day).
"""
from __future__ import annotations

import csv
import math
import os
from dataclasses import dataclass

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
PARAMS_CSV = os.path.join(ROOT, "docs", "literature_db", "raw_si", "oryza_ir72_d92_params.csv")

# 3-point Gauss (SGPCDT, SGPC1, SGPL)
GSX = (0.112702, 0.500000, 0.887298)
GSW = (0.277778, 0.444444, 0.277778)
PI = 3.1415927
DEGTRAD = 0.017453292


# ---------------------------------------------------------------------------
# Crop file and management
# ---------------------------------------------------------------------------
def load_crop(path: str = PARAMS_CSV) -> dict:
    """IR72.D92 as a dict: scalars -> float, tables -> (xs, ys) arrays, switches -> str."""
    out = {}
    with open(path) as f:
        for r in csv.DictReader(line for line in f if not line.startswith("#")):
            v = r["value"]
            if ":" in v:
                pairs = [p.split(":") for p in v.split(";")]
                out[r["parameter"]] = (np.array([float(a) for a, _ in pairs]),
                                       np.array([float(b) for _, b in pairs]))
            else:
                try:
                    out[r["parameter"]] = float(v)
                except ValueError:
                    out[r["parameter"]] = v
    return out


@dataclass
class Management:
    """Experimental-data-file settings (Section 7.2). Defaults = IR72DS92 (IR72DS2.T92)."""
    lat: float = 14.18          # LAT
    emd: int = 4                # EMD, day of emergence (in the seedbed for transplanted rice)
    sbdur: int = 12             # SBDUR, seedbed duration [d]
    nplh: float = 5.0           # NPLH, plants per hill
    nh: float = 25.0            # NH, hills per m2
    nplsb: float = 1000.0       # NPLSB, plants per m2 in the seedbed
    nplds: float = 200.0        # NPLDS, plants per m2 direct-seeded
    lape: float = 0.0001        # LAPE, initial leaf area per plant [m2]
    estab: str = "TRANSPLANT"   # ESTAB
    co2: float = 340.0          # CO2 [ppm]
    tmpsb: float = 0.0          # TMPSB, seedbed cover temperature increase [oC]
    dvsi: float = 0.0           # DVSI
    wlvgi: float = 0.0          # WLVGI
    wsti: float = 0.0           # WSTI
    wrti: float = 0.0           # WRTI
    wsoi: float = 0.0           # WSOI


def lint2(tab, x: float) -> float:
    """TTUTIL LINT2: linear interpolation, constant beyond the end points."""
    return float(np.interp(x, tab[0], tab[1]))


def notnul(x: float) -> float:
    return x if x != 0.0 else 1.0


# ---------------------------------------------------------------------------
# Phenology (3.2.1)
# ---------------------------------------------------------------------------
def subdd(tmax: float, tmin: float, tbd: float, tod: float, tmd: float) -> float:
    """SUBDD: daily heat units from 24 hourly temperatures (Eqns 3.2-3.4)."""
    tm = (tmax + tmin) / 2.0
    tt = 0.0
    for i in range(1, 25):
        td = tm + 0.5 * abs(tmax - tmin) * math.cos(0.2618 * float(i - 14))
        if tbd < td < tmd:
            if td > tod:
                td = tod - (td - tod) * (tod - tbd) / (tmd - tod)
            tt += (td - tbd) / 24.0
    return tt


# ---------------------------------------------------------------------------
# Canopy photosynthesis (3.2.2)
# ---------------------------------------------------------------------------
def sastro(idoy: int, lat: float):
    """SASTRO: solar constant, daylength and the integrals of sin(solar height)."""
    doy = float(idoy)
    dec = -math.asin(math.sin(23.45 * DEGTRAD) * math.cos(2.0 * PI * (doy + 10.0) / 365.0))
    sinld = math.sin(DEGTRAD * lat) * math.sin(dec)
    cosld = math.cos(DEGTRAD * lat) * math.cos(dec)
    aob = sinld / cosld
    if aob < -1.0:
        dayl, daylp, zzcos, zzsin = 0.0, 0.0, 0.0, 1.0
    elif aob > 1.0:
        dayl, daylp, zzcos, zzsin = 24.0, 24.0, 0.0, -1.0
    else:
        dayl = 12.0 * (1.0 + 2.0 * math.asin(aob) / PI)
        daylp = 12.0 * (1.0 + 2.0 * math.asin((-math.sin(-4.0 * DEGTRAD) + sinld) / cosld) / PI)
        zza = PI * (12.0 + dayl) / 24.0
        zzcos, zzsin = math.cos(zza), math.sin(zza)
    dsinb = 2.0 * 3600.0 * (dayl * 0.5 * sinld - 12.0 * cosld * zzcos / PI)
    dsinbe = 2.0 * 3600.0 * (dayl * (0.5 * sinld + 0.2 * sinld ** 2 + 0.1 * cosld ** 2)
                             - (12.0 * cosld * zzcos + 9.6 * sinld * cosld * zzcos
                                + 2.4 * cosld ** 2 * zzcos * zzsin) / PI)
    solcon = 1370.0 * (1.0 + 0.033 * math.cos(2.0 * PI * doy / 365.0))
    angot = solcon * dsinb
    return solcon, angot, dayl, daylp, dsinb, dsinbe, sinld, cosld


def sskyc(hour, solcon, frpar, dsinbe, sinld, cosld, rdd):
    """SSKYC: sine of solar height and direct/diffuse PAR [W m-2] at `hour`."""
    sinb = sinld + cosld * math.cos((hour - 12.0) * 0.2617993)
    rdpdf = rdpdr = 0.0
    if sinb > 0.0:
        tmpr1 = rdd * sinb * (1.0 + 0.4 * sinb) / dsinbe
        atmtr = tmpr1 / (solcon * sinb)
        if atmtr <= 0.22:
            frdif = 1.0
        elif atmtr <= 0.35:
            frdif = 1.0 - 6.4 * (atmtr - 0.22) ** 2
        else:
            frdif = 1.47 - 1.66 * atmtr
        frdif = max(frdif, 0.15 + 0.85 * (1.0 - math.exp(-0.1 / sinb)))
        rdpdf = tmpr1 * frpar * frdif
        rdpdr = tmpr1 * frpar * (1.0 - frdif)
    return sinb, rdpdr, rdpdf


def gppar_get(gai, gaid, eff_in, gp, amax_form="continuous"):
    """GPPARGET: leaf AMAX from the canopy leaf-N profile, CO2 and temperature."""
    co2, knf, nflv, redft = gp
    amaxco2 = max(0.0, 49.57 / 34.26 * (1.0 - math.exp(-0.208 * (co2 - 60.0) / 49.57)))
    if gai > 0.01 and knf > 0.0:
        slni = nflv * gai * knf * math.exp(-knf * gaid) / (1.0 - math.exp(-knf * gai))
    else:
        slni = nflv
    if slni >= 0.5:
        if amax_form == "as_printed":
            amax = 9.5 + (22.0 * slni) * redft * amaxco2
        else:
            amax = (9.5 + 22.0 * slni) * redft * amaxco2
    else:
        amax = max(0.0, 68.33 * (slni - 0.2) * redft * amaxco2)
    return amax, eff_in


def srdprf(gaid, cslv, sinb, ecpdf, rdpdr, rdpdf):
    """SRDPRF: absorbed PAR of shaded leaves, of leaves perpendicular to the beam, sunlit fraction."""
    tmpr1 = math.sqrt(1.0 - cslv)
    rflh = (1.0 - tmpr1) / (1.0 + tmpr1)
    rfls = rflh * 2.0 / (1.0 + 2.0 * sinb)
    clustf = ecpdf / (0.8 * tmpr1)
    ecpbl = (0.5 / sinb) * clustf
    ecptd = ecpbl * tmpr1
    rapdfl = (1.0 - rflh) * rdpdf * ecpdf * math.exp(-ecpdf * gaid)
    raptdl = (1.0 - rfls) * rdpdr * ecptd * math.exp(-ecptd * gaid)
    rapddl = (1.0 - cslv) * rdpdr * ecpbl * math.exp(-ecpbl * gaid)
    rapshl = rapdfl + raptdl - rapddl
    rapppl = (1.0 - cslv) * rdpdr / sinb
    fslla = clustf * math.exp(-ecpbl * gaid)
    return rapshl, rapppl, fslla


def sgpl(cslv, eff1, ecpdf, gai, gaid, sinb, rdpdr, rdpdf, gp, amax_form):
    """SGPL: gross assimilation and absorbed PAR per unit leaf area at depth GAID."""
    rapshl, rapppl, fslla = srdprf(gaid, cslv, sinb, ecpdf, rdpdr, rdpdf)
    amax2, eff2 = gppar_get(gai, gaid, eff1, gp, amax_form)
    gpshl = amax2 * (1.0 - math.exp(-rapshl * eff2 / amax2)) if amax2 > 0.0 else 0.0
    gpsll = rapsll = 0.0
    for x, w in zip(GSX, GSW):
        tmpr1 = rapshl + rapppl * x
        gpsll = gpsll + amax2 * (1.0 - math.exp(-tmpr1 * eff2 / amax2)) * w if amax2 > 0.0 else 0.0
        rapsll += tmpr1 * w
    gpl = fslla * gpsll + (1.0 - fslla) * gpshl
    rapl = fslla * rapsll + (1.0 - fslla) * rapshl
    return gpl, rapl


def sgpc1(cslv, eff, ecpdf, gai, sinb, rdpdr, rdpdf, gp, amax_form):
    """SGPC1: instantaneous canopy assimilation [kg CO2 ha-1 h-1] and absorbed PAR [W m-2]."""
    gpc = rapc = 0.0
    for x, w in zip(GSX, GSW):
        gpl, rapl = sgpl(cslv, eff, ecpdf, gai, gai * x, sinb, rdpdr, rdpdf, gp, amax_form)
        gpc += gpl * w
        rapc += rapl * w
    return gpc * gai, rapc * gai


def sgpcdt(idoy, lat, rdd, frpar, cslv, eff, ecpdf, gai, gp, amax_form="continuous"):
    """SGPCDT: daily gross canopy assimilation DTGA [kg CO2 ha-1 d-1] (IACC = 1).

    Returns (DTGA, RAPCDT [J m-2 d-1], DAYL [h]).  ``rdd`` is daily global radiation in
    J m-2 d-1 (the WEATHER system converts the file's kJ to J)."""
    solcon, _angot, dayl, _daylp, _dsinb, dsinbe, sinld, cosld = sastro(idoy, lat)
    gpcdt = rapcdt = 0.0
    if gai > 0.0:
        for x, w in zip(GSX, GSW):
            hour = 12.0 + dayl * 0.5 * x
            sinb, rdpdr, rdpdf = sskyc(hour, solcon, frpar, dsinbe, sinld, cosld, rdd)
            if sinb <= 0.0:
                continue
            gpc, rapc = sgpc1(cslv, eff, ecpdf, gai, sinb, rdpdr, rdpdf, gp, amax_form)
            gpcdt += gpc * w
            rapcdt += rapc * w
    return gpcdt * dayl, rapcdt * dayl * 3600.0, dayl


# ---------------------------------------------------------------------------
# The ORYZA1 daily loop (Sections 2.2 call order, 3.1-3.3)
# ---------------------------------------------------------------------------
_RECORD = ("doy", "cropsta", "dvs", "lai", "wlvg", "wlvd", "wsts", "wstr", "wst", "wso",
           "wrt", "wagt", "wrr", "ngr", "nsp", "dtga", "gcr", "llv", "drlv", "ts", "tslv",
           "cbc_rel")


def simulate(weather: dict, crop: dict | None = None, mgmt: Management | None = None,
             max_days: int = 300, amax_form: str = "continuous") -> dict:
    """Run ORYZA2000 potential production from emergence until DVS > 2 (or `max_days`).

    weather: dict with ``doy0`` (day of year of the first entry) and daily arrays ``rad_mj``
    (global radiation, MJ m-2 d-1), ``tmin`` and ``tmax`` (oC) that cover emergence onward.
    Returns daily arrays (state at the start of each simulated day, as ORYZA1 writes them)
    plus ``events`` (day of year of transplanting, panicle initiation, flowering, maturity)
    and ``cbc_max`` (largest relative carbon-balance error, the SUBCBC check)."""
    c = crop if crop is not None else load_crop()
    m = mgmt or Management()
    doy0 = int(weather["doy0"])
    rad = np.asarray(weather["rad_mj"], float) * 1.0e6          # MJ -> J m-2 d-1
    tmn_w = np.asarray(weather["tmin"], float)
    tmx_w = np.asarray(weather["tmax"], float)
    k0 = m.emd - doy0
    if k0 < 0:
        raise ValueError("weather starts after the emergence day")

    # --- NNOSTRESS / WNOSTRESS (potential production) ---
    nflv, nsllv, rnstrs = c["NFLVI"], 1.0, 1.0
    lrstrs = ldstrs = lestrs = pcew = 1.0

    # --- reinitialization at emergence (ORYZA1, CROPSTA = 1) ---
    cropsta = 1
    dvs, wlvg, wlvd = m.dvsi, m.wlvgi, 0.0
    wsts, wstr, wso, wrt = m.wsti, 0.0, m.wsoi, m.wrti
    lai = m.lape * (m.nplsb if m.estab == "TRANSPLANT" else m.nplds)
    ts = tslv = wrr = ngr = nsp = dae = tnass = pwrr = 0.0
    gcr_prev = ggr_prev = 0.0
    dayl_prev = 12.0
    # SUBLAI2 / PHENOL / SUBGRN saved variables
    x_lai, testl, testset = 1.0, False, 0.00001
    wlvgexs = laiexs = wlvgexp = laiexp = 0.0
    tslvtr = tshckl = tstr = 0.0
    coldtt = tfert = ntfert = 0.0
    grains = False
    events = {}
    rec = {k: [] for k in _RECORD}
    cbc_max = 0.0

    for step in range(max_days):
        k = k0 + step
        if k >= len(rad):
            raise ValueError("weather series too short for the crop cycle")
        idoy = (doy0 - 1 + k) % 365 + 1
        wst = wsts + wstr

        # ======== rate calculations (ITASK = 2) ========
        tmpcov = m.tmpsb if cropsta <= 2 else 0.0
        tmax = tmx_w[k] + tmpcov
        tmin = tmn_w[k]
        tav = (tmin + tmax) / 2.0
        tavd = (tmax + tav) / 2.0
        dtr = rad[k]

        # phenology: SUBDD, PHENOL (uses yesterday's DAYL, as in ORYZA1), drought factor
        hu = subdd(tmax, tmin, c["TBD"], c["TOD"], c["TMD"])
        if dvs < 0.40:
            dvr = c["DVRJ"] * hu
        elif dvs < 0.65:
            dl = dayl_prev + 0.9
            ppfac = 1.0 if dl < c["MOPP"] else 1.0 - (dl - c["MOPP"]) * c["PPSE"]
            dvr = c["DVRI"] * hu * min(1.0, max(0.0, ppfac))
        elif dvs < 1.00:
            dvr = c["DVRP"] * hu
        else:
            dvr = c["DVRR"] * hu
        if cropsta == 3:
            tstr = ts
        tshckd = c["SHCKD"] * tstr
        if cropsta > 3 and ts < (tstr + tshckd):
            dvr = 0.0
        dvew = lestrs + dvs * (1.0 - lestrs) if dvs < 1.0 else 1.0
        dvr *= dvew

        # CO2 assimilation
        co2eff = ((1.0 - math.exp(-0.00305 * m.co2 - 0.222))
                  / (1.0 - math.exp(-0.00305 * c["CO2REF"] - 0.222)))
        eff = lint2(c["EFFTB"], tavd) * co2eff
        lairol = lai * (0.5 * lrstrs + 0.5)
        sai = lint2(c["SSGATB"], dvs) * wst
        alai = lairol + 0.5 * sai
        kdf = lint2(c["KDFTB"], dvs)
        redft = lint2(c["REDFTT"], tavd)
        knf = lint2(c["KNFTB"], dvs)
        dtga, _rapcdt, dayl = sgpcdt(idoy, m.lat, dtr, c["FRPAR"], c["SCP"], eff, kdf, alai,
                                     (m.co2, knf, nflv, redft), amax_form)
        dayl_prev = dayl
        dtga *= pcew

        # partitioning (+ drought effect, + grain sink limit from yesterday's rates)
        fsh = lint2(c["FSHTB"], dvs)
        if dvs < 1.0:
            fsh *= lestrs
        frt = 1.0 - fsh
        flv = lint2(c["FLVTB"], dvs)
        fst = lint2(c["FSTTB"], dvs)
        fso = lint2(c["FSOTB"], dvs)
        if grains and ggr_prev >= (pwrr - wrr):
            fso = max(0.0, (pwrr - wrr) / notnul(gcr_prev * fsh))
            fst = 1.0 - fso - flv

        # loss of green leaves and stem reserves
        drlv = lint2(c["DRLVT"], dvs)
        llv = nsllv * wlvg * drlv
        lstr = 0.0 if dvs - 1.0 < 0.0 else wstr / c["TCLSTR"]          # INSW(DVS-1., 0., WSTR/TCLSTR)

        # respiration and crop growth rate
        teff = c["Q10"] ** ((tav - c["TREF"]) / 10.0)
        mndvs = wlvg / notnul(wlvg + wlvd)
        rmcr = (wlvg * c["MAINLV"] + wst * c["MAINST"] + wso * c["MAINSO"]
                + wrt * c["MAINRT"]) * teff * mndvs
        crgcr = (fsh * (c["CRGLV"] * flv + c["CRGST"] * fst * (1.0 - c["FSTR"])
                        + c["CRGSTR"] * c["FSTR"] * fst + c["CRGSO"] * fso)
                 + c["CRGRT"] * frt)
        gcr = ((dtga * 30.0 / 44.0) - rmcr + (lstr * c["LRSTR"] * c["FCSTR"] * 30.0 / 12.0)) / crgcr

        # transplanting dilution and organ growth rates (3.2.7)
        pltr = m.nplh * m.nh / m.nplsb if cropsta == 3 else 1.0
        rwlvg1 = wlvg * (1.0 - pltr)
        gst1 = wsts * (1.0 - pltr)
        rwstr1 = wstr * (1.0 - pltr)
        grt1 = wrt * (1.0 - pltr)
        grt = gcr * frt - grt1
        glv = gcr * fsh * flv - rwlvg1
        rwlvg = glv - llv
        gst = gcr * fsh * fst * (1.0 - c["FSTR"]) - gst1
        gstr = gcr * fsh * fst * c["FSTR"] - rwstr1
        rwstr = gstr - lstr
        gso = gcr * fsh * fso
        ggr = gso if dvs > 0.95 else 0.0

        # SUBGRN: spikelets and grains
        tincr = 5.0 * (1.0 - lrstrs) * 1.6
        gnsp = gcr * c["SPGF"] if 0.65 <= dvs <= 1.0 else 0.0
        if 0.75 <= dvs <= 1.2:
            coldtt += max(0.0, 22.0 - (tav - tincr))
        if 0.96 <= dvs <= 1.2:
            tfert += tmax + tincr
            ntfert += 1.0
        if dvs >= 1.2 and not grains:
            grains = True
            sf1 = min(1.0, max(0.0, 1.0 - (4.6 + 0.054 * coldtt ** 1.56) / 100.0))
            tf = tfert / ntfert if ntfert > 0 else tmax
            sf2 = min(1.0, max(0.0, 1.0 / (1.0 + math.exp(0.853 * (tf - 36.6)))))
            gngr = nsp * min(sf1, sf2)
        else:
            gngr = 0.0

        # leaf area growth: SUBDD for leaves, SLA, SUBLAI2 (transplanted / direct-seeded)
        hulv = subdd(tmax, tmin, c["TBLV"], 30.0, 42.0)
        if c.get("SWISLA", "FUNCTION") == "TABLE" and "SLATB" in c:
            sla = lint2(c["SLATB"], dvs)
        else:
            sla = min(c["SLAMAX"], c["ASLA"] + c["BSLA"] * math.exp(c["CSLA"] * (dvs - c["DSLA"])))
        dldr = 0.0                                               # no drought leaf death
        rgrl = c["RGRLMX"] - (1.0 - rnstrs) * (c["RGRLMX"] - c["RGRLMN"])
        if m.estab == "TRANSPLANT":
            if cropsta < 3:                                      # 1. seedbed
                if lai < 1.0:
                    glai = lai * rgrl * hulv
                    wlvgexs, laiexs = wlvg, lai
                else:
                    test = abs((lai / notnul(wlvg)) - sla) / sla
                    if not testl and test < testset:
                        testl = True
                    if testl:
                        glai = ((wlvg + rwlvg) * sla) - lai
                    else:
                        glai1 = ((wlvg + rwlvg - wlvgexs) * sla + laiexs) - lai
                        glai2 = ((wlvg + rwlvg) * sla) - lai
                        glai = (glai1 + x_lai * glai2) / (x_lai + 1.0)
                        x_lai += 1.0
            elif cropsta == 3:                                   # 2. transplanting
                tslvtr = tslv
                tshckl = c["SHCKL"] * tslvtr
                glai = (lai * m.nh * m.nplh / m.nplsb) - lai
            else:                                                # 3. main field
                if tslv < (tslvtr + tshckl):
                    glai = 0.0
                elif lai < 1.0 and dvs < 1.0:
                    glai = lestrs * lai * rgrl * hulv
                    wlvgexp, laiexp = wlvg, lai
                else:
                    test = abs((lai / notnul(wlvg)) - sla) / sla
                    if not testl and test < testset:
                        testl = True
                    if testl:
                        glai = ((wlvg + rwlvg - dldr) * sla) - lai
                    else:
                        glai1 = ((wlvg + rwlvg - dldr - wlvgexp) * sla + laiexp) - lai
                        glai2 = ((wlvg + rwlvg - dldr) * sla) - lai
                        glai = (glai1 + x_lai * glai2) / (x_lai + 1.0)
                        x_lai += 1.0
        else:                                                    # direct-seeded
            if lai < 1.0 and dvs < 1.0:
                glai = lai * rgrl * hulv * lestrs
                wlvgexp, laiexp = wlvg, lai
            else:
                test = abs((lai / notnul(wlvg)) - sla) / sla
                if not testl and test < testset:
                    testl = True
                if testl:
                    glai = ((wlvg + rwlvg - dldr) * sla) - lai
                else:
                    glai1 = ((wlvg + rwlvg - dldr - wlvgexp) * sla + laiexp) - lai
                    glai2 = ((wlvg + rwlvg - dldr) * sla) - lai
                    glai = (glai1 + x_lai * glai2) / (x_lai + 1.0)
                    x_lai += 1.0

        # growth respiration and the net CO2 assimilation (3.2.11)
        co2rt = 44.0 / 12.0 * (c["CRGRT"] * 12.0 / 30.0 - c["FCRT"])
        co2lv = 44.0 / 12.0 * (c["CRGLV"] * 12.0 / 30.0 - c["FCLV"])
        co2st = 44.0 / 12.0 * (c["CRGST"] * 12.0 / 30.0 - c["FCST"])
        co2str = 44.0 / 12.0 * (c["CRGSTR"] * 12.0 / 30.0 - c["FCSTR"])
        co2so = 44.0 / 12.0 * (c["CRGSO"] * 12.0 / 30.0 - c["FCSO"])
        rgcr = ((grt + grt1) * co2rt + (glv + rwlvg1) * co2lv + (gst + gst1) * co2st
                + gso * co2so + (gstr + rwstr1) * co2str
                + (1.0 - c["LRSTR"]) * lstr * c["FCSTR"] * 44.0 / 12.0)
        ctrans = rwlvg1 * c["FCLV"] + gst1 * c["FCST"] + rwstr1 * c["FCSTR"] + grt1 * c["FCRT"]
        rtnass = ((dtga * 30.0 / 44.0 - rmcr) * 44.0 / 30.0) - rgcr - (ctrans * 44.0 / 12.0)

        # carbon balance check on the states at the start of this day (SUBCBC)
        ckcin = ((wlvg + wlvd - m.wlvgi) * c["FCLV"] + (wsts - m.wsti) * c["FCST"]
                 + wstr * c["FCSTR"] + (wrt - m.wrti) * c["FCRT"] + wso * c["FCSO"])
        ckcfl = tnass * (12.0 / 44.0)
        cbc_rel = abs(ckcin - ckcfl) / max(abs(ckcin), abs(ckcfl), 1e-9)
        if ckcin > 1.0:
            cbc_max = max(cbc_max, cbc_rel)

        # record the start-of-day state (ORYZA1 writes its output in the rate section)
        for key, val in (("doy", idoy), ("cropsta", cropsta), ("dvs", dvs), ("lai", lai),
                         ("wlvg", wlvg), ("wlvd", wlvd), ("wsts", wsts), ("wstr", wstr),
                         ("wst", wst), ("wso", wso), ("wrt", wrt),
                         ("wagt", wlvg + wlvd + wst + wso), ("wrr", wrr), ("ngr", ngr),
                         ("nsp", nsp), ("dtga", dtga), ("gcr", gcr), ("llv", llv),
                         ("drlv", drlv * nsllv), ("ts", ts), ("tslv", tslv),
                         ("cbc_rel", cbc_rel)):
            rec[key].append(val)
        if cropsta == 3:
            events.setdefault("transplanting", idoy)
        for name, lim in (("panicle_initiation", 0.65), ("flowering", 1.0), ("maturity", 2.0)):
            if name not in events and dvs + dvr >= lim > dvs:
                events[name] = idoy + 1                          # state crosses during this day

        # NNOSTRESS rate (called after ORYZA1: NFLV takes effect the next day)
        nflv = lint2(c["NFLVTB"], dvs)

        # ======== integration (ITASK = 3) ========
        ts += hu
        tslv += hulv
        dvs += dvr
        wlvg += rwlvg - dldr
        wlvd += llv + dldr
        wsts += gst
        wstr += rwstr
        wso += gso
        wrt += grt
        wrr += ggr
        ngr += gngr
        nsp += gnsp
        dae += 1.0
        tnass += rtnass
        lai += glai
        pwrr = ngr * c["WGRMX"]
        gcr_prev, ggr_prev = gcr, ggr

        # MODELS: crop stage update after the integration
        if cropsta == 3:
            cropsta = 4
        if cropsta == 2 and dae == float(m.sbdur):
            cropsta = 3
        if cropsta == 1:
            cropsta = 2 if m.estab == "TRANSPLANT" else 4

        if dvs > 2.0:                                            # TERMNL
            idoy_end = (doy0 - 1 + k + 1) % 365 + 1
            for key, val in (("doy", idoy_end), ("cropsta", cropsta), ("dvs", dvs),
                             ("lai", lai), ("wlvg", wlvg), ("wlvd", wlvd), ("wsts", wsts),
                             ("wstr", wstr), ("wst", wsts + wstr), ("wso", wso), ("wrt", wrt),
                             ("wagt", wlvg + wlvd + wsts + wstr + wso), ("wrr", wrr),
                             ("ngr", ngr), ("nsp", nsp), ("dtga", 0.0), ("gcr", 0.0),
                             ("llv", 0.0), ("drlv", lint2(c["DRLVT"], dvs) * nsllv),
                             ("ts", ts), ("tslv", tslv), ("cbc_rel", cbc_rel)):
                rec[key].append(val)
            break

    out = {k: np.asarray(v, float) for k, v in rec.items()}
    out["events"] = events
    out["cbc_max"] = cbc_max
    return out


# ---------------------------------------------------------------------------
# Bridge to the PFAS model: per-organ dry biomass [kg/hill] on days after transplanting
# ---------------------------------------------------------------------------
# The generic scenario's planting density: the crop run's own NH (IR72DS92, inherited verbatim with
# the crop file), so the per-hill conversion and the run agree. Defined here ONCE; forcing_rice,
# growth_rice and oryza_growth read it (PREREG 7.4, 7.12). Every BAF is invariant to it (Q_TP and M
# per hill scale together); only per-hill absolute amounts depend on it. The site-year scenario
# uses its own sourced density (site_scenario.NH, 18.5 hills/m2).
GENERIC_HILLS_PER_M2 = Management().nh
AREA_PER_HILL_M2 = 1.0 / GENERIC_HILLS_PER_M2
HA_PER_HILL = AREA_PER_HILL_M2 / 1.0e4


def default_weather(transplant_doy: int = 180, sbdur: int = 12, n_days: int = 200) -> dict:
    """The repo's representative tropical wet-season climatology (oryza_growth._weather),
    indexed from the transplanting day and extended back over the seedbed days."""
    d = np.arange(-sbdur, n_days - sbdur, dtype=float)
    return dict(doy0=transplant_doy - sbdur,
                rad_mj=19.0 + 3.0 * np.sin(2 * np.pi * (d - 20) / 120.0),
                tmax=31.0 + 1.5 * np.sin(2 * np.pi * (d - 30) / 120.0),
                tmin=23.0 + 1.0 * np.sin(2 * np.pi * (d - 30) / 120.0))


def organ_biomass(t, season: float = 120.0, weather: dict | None = None,
                  mgmt: Management | None = None, transplant_doy: int = 180,
                  crop: dict | None = None, amax_form: str = "continuous") -> dict:
    """Per-organ DRY biomass [kg/hill] on `t` = days after transplanting (root/stem/leaf/grain).

    root = WRT, stem = WST (structural + reserves), leaf = WLVG (green leaves; dead leaves
    have left the plant), grain = WSO (panicles).  t = 0 is the main-field state right after
    transplanting (the seedbed state diluted by NPLH*NH/NPLSB).  After maturity the states
    are held and the leaf death rate is 0.  ``leaf_death_rate`` [1/d] = DRLV(DVS) = LLV/WLVG,
    the rate the PFAS leaf ODE removes with the dead leaf (as for oryza_growth).

    No post-hoc scaling: ORYZA2000 with IR72.D92 is calibrated to give attainable IR72 yields
    at IRRI (validation/oryza2000_ir72ds92.py).  ``weather`` (optional) is in `simulate`'s
    format (doy0 + daily arrays covering emergence onward); default = the repo's tropical
    climatology (`default_weather`)."""
    t = np.asarray(t, float)
    m = mgmt or Management(emd=transplant_doy - 12)
    trans = m.emd + m.sbdur
    w = weather or default_weather(trans, m.sbdur, int(np.ceil(max(season, t.max()))) + m.sbdur + 60)
    r = simulate(w, crop, m, amax_form=amax_form)
    sta = r["cropsta"]
    i3 = int(np.where(sta == 3)[0][0])
    pltr = m.nplh * m.nh / m.nplsb
    main = np.where(sta >= 4)[0]
    td = np.concatenate([[0.0], r["doy"][main] - trans])
    td = np.where(td < 0, td + 365, td)                    # year wrap
    organs = {"root": "wrt", "stem": "wst", "leaf": "wlvg", "grain": "wso"}
    out = {}
    for name, key in organs.items():
        v = np.concatenate([[r[key][i3] * pltr], r[key][main]])
        out[name] = np.interp(t, td, v * HA_PER_HILL)     # crop model's own value (0 before formation)
    drlv = np.concatenate([[r["drlv"][i3]], r["drlv"][main]])
    if r["dvs"][-1] > 2.0:
        drlv[-1] = 0.0
    t_mat = td[-1]
    out["leaf_death_rate"] = np.where(t >= t_mat, 0.0, np.maximum(np.interp(t, td, drlv), 0.0))
    out["dvs"] = np.interp(t, td, np.concatenate([[r["dvs"][i3]], r["dvs"][main]]))
    return out
