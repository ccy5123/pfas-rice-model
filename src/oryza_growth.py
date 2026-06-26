"""
oryza_growth.py — ORYZA2000 potential-production (Level 1) growth core, in Python
=================================================================================

A faithful re-implementation of the ORYZA2000 / ORYZA(v3) **potential production**
crop-growth mechanism (carbon-driven), as opposed to the lightweight DVS-partition-
on-a-logistic reconstruction in `growth_rice.py`.  It actually runs the ORYZA carbon
balance on a daily step:

    radiation + temperature
        -> SUCROS astronomy (daylength, solar height)            [ASTRO]
        -> Gaussian day x canopy gross CO2 assimilation          [TOTASS/ASSIM]
        -> maintenance + growth respiration -> net assimilate
        -> DVS-driven partitioning (root/leaf/stem/panicle)      [IR72 tables]
        -> SLA-based leaf-area growth (+ juvenile RGRL, senescence)
        -> grain filling, maturity at DVS=2

so the per-organ biomass M_s(t) RESPONDS to weather (latitude, sowing date, daily
radiation, Tmax/Tmin) instead of being imposed.  Output is per-organ fresh-free dry
biomass [kg/hill] on a requested time grid, plus LAI(t) and DVS(t).

This is NOT the IRRI ORYZA executable (a Windows binary needing a full weather/crop/
management input deck); it is the published ORYZA2000 Level-1 equation set
re-coded in Python.  Structure & defaults follow:
  * Bouman & van Laar 2006, Agric. Syst. 87:249-273 (ORYZA2000 description/eval)   10.1016/j.agsy.2004.07.011
  * Bouman et al. 2001, "ORYZA2000: modeling lowland rice" (IRRI/WUR manual)
  * Goudriaan & van Laar 1994 (SUCROS canopy assimilation; ASTRO/TOTASS)
  * Li et al. 2017 (ORYZA v3)                                                       10.1016/j.agrformet.2017.02.025
IR72 crop parameters are the ORYZA standard-set order of magnitude, anchored so the
potential run reproduces the IR72 field anchors (shoot ~1740 g/m^2, HI ~0.5,
flowering ~DVS1, LAImax ~5-7).  Parameter provenance is flagged inline.

Use:
    t = np.linspace(0,120,481)
    b = organ_biomass_oryza(t)             # dict root/stem/leaf/grain [kg/hill]
    drivers = oryza_drivers("PFOA")        # -> model_api.simulate(drivers=...)
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np

# planting density convention shared with forcing_rice / growth_rice
AREA_PER_HILL_M2 = 0.04                       # 25 hills/m^2 (transplanted paddy)
HA_PER_HILL = AREA_PER_HILL_M2 / 1.0e4        # kg/ha -> kg/hill

# 3-point Gaussian (Goudriaan) for day & canopy integration
_XG = np.array([0.1127, 0.5, 0.8873])
_WG = np.array([0.2778, 0.4444, 0.2778])


# ---------------------------------------------------------------------------
# Crop / site parameters (ORYZA2000 IR72 standard-set order of magnitude)
# ---------------------------------------------------------------------------
@dataclass
class OryzaParams:
    # --- site / sowing ---
    lat_deg: float = 14.18                    # IRRI Los Banos (IR72 anchor context)
    sow_doy: int = 180                        # day-of-year of transplanting
    season: float = 120.0                     # days transplant->maturity (cap)

    # --- phenology: development rates DVR [1/(deg C d)], Tbase=8C (ORYZA IR72) ---
    Tbase: float = 8.0
    DVRJ: float = 8.00e-4                     # juvenile  (DVS 0   -> 0.4)
    DVRI: float = 7.90e-4                     # photoper. (DVS 0.4 -> 0.65)
    DVRP: float = 8.00e-4                     # panicle   (DVS 0.65-> 1.0)
    DVRR: float = 1.05e-3                     # repro/fill(DVS 1.0 -> 2.0); grain fill ~50 d

    # --- assimilation (SUCROS/ORYZA) ---
    AMAXTB_dvs: tuple = (0.0, 1.0, 1.6, 2.0)  # max leaf assim vs DVS
    AMAXTB_val: tuple = (50.0, 50.0, 38.0, 14.0)   # kg CO2 ha-1 leaf h-1 (IR72-ish)
    EFF: float = 0.45                         # initial light-use eff [kg CO2 ha-1 h-1 /(J m-2 s-1)]
    KDF: float = 0.50                         # diffuse extinction coeff [-]
    SCV: float = 0.20                         # leaf scattering coeff (PAR) [-]
    # temperature reduction factor on AMAX (REDFTT vs Tavg) -- ORYZA table
    REDFTT_T: tuple = (8.0, 15.0, 20.0, 30.0, 37.0, 43.0)
    REDFTT_F: tuple = (0.0, 0.55, 0.85, 1.0, 1.0, 0.0)

    # --- respiration ---
    MAINLV: float = 0.020                     # maintenance [kg CH2O kg-1 d-1]
    MAINST: float = 0.015
    MAINSO: float = 0.003
    MAINRT: float = 0.010
    Q10: float = 2.0
    Tref: float = 25.0
    # growth conversion efficiency [kg DM / kg CH2O] per organ
    CVL: float = 0.72
    CVST: float = 0.69
    CVSO: float = 0.71
    CVR: float = 0.72

    # --- leaf area ---
    RGRL: float = 0.0075                      # juvenile rel. leaf-area growth [ (ha/ha)/(deg C d) ]
    SLATB_dvs: tuple = (0.0, 0.16, 0.33, 0.65, 1.0, 2.5)
    SLATB_val: tuple = (0.0037, 0.0037, 0.0028, 0.0024, 0.0020, 0.0016)  # SLA [ha kg-1 leaf]
    LAII: float = 0.0085                      # initial LAI at transplant (ha/ha)
    # leaf death (senescence) rate after flowering [d-1] vs DVS
    DRLVT_dvs: tuple = (1.0, 1.6, 2.0)
    DRLVT_val: tuple = (0.0, 0.015, 0.045)

    # --- partitioning vs DVS (ORYZA FSH then FLV/FST/FSO; IR72) ---
    FSH_dvs: tuple = (0.0, 0.43, 1.0, 2.5)    # shoot fraction of total growth
    FSH_val: tuple = (0.50, 0.75, 1.00, 1.00)
    FLV_dvs: tuple = (0.0, 0.5, 0.75, 1.0, 1.2, 2.5)
    FLV_val: tuple = (0.60, 0.55, 0.50, 0.05, 0.00, 0.00)
    FST_val: tuple = (0.40, 0.45, 0.50, 0.50, 0.30, 0.25)
    FSO_val: tuple = (0.00, 0.00, 0.00, 0.45, 0.70, 0.75)

    # --- initial total biomass (transplanted seedling) [kg DM/ha] ---
    WLVI: float = 54.0
    WSTI: float = 22.0
    WRTI: float = 20.0

    # --- IR72 field anchor for a light post-hoc scale (potential -> attainable) ---
    shoot_anchor_g_m2: float = 1740.0         # Bouman 2006 IR72 aboveground


def _interp(x, xp, fp):
    return float(np.interp(x, xp, fp))


# ---------------------------------------------------------------------------
# SUCROS astronomy + daily canopy gross assimilation
# ---------------------------------------------------------------------------
def astro(doy: int, lat_deg: float):
    """Daylength and solar-height integrals (Goudriaan & van Laar 1994, ASTRO)."""
    rad = np.pi / 180.0
    dec = -23.45 * np.cos(2.0 * np.pi * (doy + 10) / 365.0) * rad
    sinld = np.sin(lat_deg * rad) * np.sin(dec)
    cosld = np.cos(lat_deg * rad) * np.cos(dec)
    aob = np.clip(sinld / cosld, -1.0, 1.0)
    dayl = 12.0 * (1.0 + 2.0 * np.arcsin(aob) / np.pi)
    dsinbe = 3600.0 * (dayl * (sinld + 0.4 * (sinld ** 2 + 0.5 * cosld ** 2))
                       + 12.0 * cosld * (2.0 + 3.0 * 0.4 * sinld)
                       * np.sqrt(1.0 - aob ** 2) / np.pi)
    return dayl, sinld, cosld, dsinbe


def daily_assim(doy: int, dtr_jm2d: float, lai: float, amax: float, p: OryzaParams) -> float:
    """Daily gross canopy CO2 assimilation DTGA [kg CO2 ha-1 d-1] (SUCROS TOTASS).

    Gaussian integration over the day (3 pts) and the canopy depth (3 pts), with a
    single-stream absorbed-PAR profile (1-SCV) * KDF * exp(-KDF*L).  This is the
    ORYZA assimilation structure with a reduced (diffuse-equivalent) radiation
    submodel -- it preserves the LAI saturation + daylength/radiation response.
    """
    if lai <= 0.0 or amax <= 0.0:
        return 0.0
    dayl, sinld, cosld, dsinbe = astro(doy, p.lat_deg)
    if dayl <= 0.0:
        return 0.0
    kbl = (1.0 - p.SCV)
    dtga = 0.0
    for i in range(3):
        hour = 12.0 + dayl * 0.5 * _XG[i]
        sinb = max(0.0, sinld + cosld * np.cos(2.0 * np.pi * (hour - 12.0) / 24.0))
        # instantaneous PAR flux at the top of the canopy [J m-2 s-1] (~W/m2)
        par = 0.5 * dtr_jm2d * sinb * (1.0 + 0.4 * sinb) / dsinbe
        fgros = 0.0
        for j in range(3):
            laic = lai * _XG[j]
            parl = kbl * p.KDF * par * np.exp(-p.KDF * laic)   # absorbed PAR by layer
            fgl = amax * (1.0 - np.exp(-p.EFF * parl / max(amax, 1e-9)))
            fgros += fgl * _WG[j]
        fgros *= lai                                            # integrate over LAI
        dtga += fgros * _WG[i]
    dtga *= dayl                                                # h -> per day
    return dtga


# ---------------------------------------------------------------------------
# Daily integration of the ORYZA Level-1 carbon balance
# ---------------------------------------------------------------------------
def simulate_oryza(p: OryzaParams | None = None, weather: dict | None = None):
    """Run the daily ORYZA Level-1 loop. Returns a dict of daily arrays:
    t, dvs, lai, wrt, wlv, wst, wso, wagt(shoot), and the partition fractions.
    `weather` (optional) provides arrays `tmax`,`tmin`,`rad_mj` of length ceil(season)+1;
    otherwise a representative tropical wet-season climatology is generated."""
    p = p or OryzaParams()
    n = int(np.ceil(p.season)) + 1
    t = np.arange(n, dtype=float)
    tmax, tmin, rad = _weather(p, weather, n)

    # state [kg DM/ha]
    wrt, wlv, wst, wso = p.WRTI, p.WLVI, p.WSTI, 0.0
    lai = p.LAII
    dvs = 0.0
    out = {k: np.zeros(n) for k in
           ("dvs", "lai", "wrt", "wlv", "wst", "wso", "wagt",
            "fsh", "flv", "fst", "fso", "cgr", "drlv")}

    for d in range(n):
        tavg = 0.5 * (tmax[d] + tmin[d])
        doy = p.sow_doy + d
        # --- assimilation ---
        amax = _interp(dvs, p.AMAXTB_dvs, p.AMAXTB_val) \
            * _interp(tavg, p.REDFTT_T, p.REDFTT_F)
        dtga = daily_assim(doy, rad[d] * 1.0e6, lai, amax, p)     # MJ->J
        # --- respiration (maintenance, temperature-corrected) ---
        teff = p.Q10 ** ((tavg - p.Tref) / 10.0)
        rm = (p.MAINLV * wlv + p.MAINST * wst + p.MAINSO * wso
              + p.MAINRT * wrt) * teff
        asrc = max(0.0, dtga * 30.0 / 44.0 - rm)                  # net CH2O [kg/ha/d]
        # --- partitioning ---
        fsh = _interp(dvs, p.FSH_dvs, p.FSH_val)
        flv = _interp(dvs, p.FLV_dvs, p.FLV_val)
        fst = _interp(dvs, p.FLV_dvs, p.FST_val)
        fso = _interp(dvs, p.FLV_dvs, p.FSO_val)
        s = flv + fst + fso
        flv, fst, fso = (flv / s, fst / s, fso / s) if s > 0 else (0.0, 0.0, 0.0)
        # organ-weighted conversion efficiency (CH2O -> DM)
        cvf = fsh * (flv * p.CVL + fst * p.CVST + fso * p.CVSO) + (1 - fsh) * p.CVR
        gcrop = asrc * cvf                                        # total DM growth [kg/ha/d]
        grt = (1 - fsh) * gcrop
        gshoot = fsh * gcrop
        glv, gst, gso = flv * gshoot, fst * gshoot, fso * gshoot

        # --- leaf area: juvenile exponential, else SLA-limited; minus senescence ---
        sla = _interp(dvs, p.SLATB_dvs, p.SLATB_val)
        hu = max(0.0, tavg - p.Tbase)
        if dvs < 1.0 and lai < 1.0:
            glai = lai * p.RGRL * hu                              # source-independent juvenile
        else:
            glai = sla * glv                                     # carbon-limited
        drlv = _interp(dvs, p.DRLVT_dvs, p.DRLVT_val)            # death rate [d-1]
        dlai = drlv * lai
        dlv = drlv * wlv

        # record (start-of-day state)
        out["dvs"][d], out["lai"][d] = dvs, lai
        out["wrt"][d], out["wlv"][d] = wrt, wlv
        out["wst"][d], out["wso"][d] = wst, wso
        out["wagt"][d] = wlv + wst + wso
        out["fsh"][d], out["flv"][d] = fsh, flv
        out["fst"][d], out["fso"][d] = fst, fso
        out["cgr"][d] = gcrop
        out["drlv"][d] = drlv                                    # leaf death rate [1/d]

        # --- integrate (Euler, dt=1 d) ---
        wrt += grt
        wlv += glv - dlv
        wst += gst
        wso += gso
        lai = max(1e-6, lai + glai - dlai)
        # --- phenology ---
        if dvs < 0.4:
            dvr = p.DVRJ
        elif dvs < 0.65:
            dvr = p.DVRI
        elif dvs < 1.0:
            dvr = p.DVRP
        else:
            dvr = p.DVRR
        dvs += dvr * hu
        if dvs >= 2.0:                                            # maturity: stop filling
            for k in range(d + 1, n):
                for key, val in (("dvs", 2.0), ("lai", lai), ("wrt", wrt),
                                 ("wlv", wlv), ("wst", wst), ("wso", wso),
                                 ("wagt", wlv + wst + wso)):
                    out[key][k] = val
            break

    out["t"] = t
    return out


def _weather(p: OryzaParams, weather: dict | None, n: int):
    """Daily Tmax/Tmin [C] and global radiation [MJ m-2 d-1].
    Representative tropical wet-season climatology (overridable via `weather`)."""
    if weather is not None:
        tmax = np.asarray(weather["tmax"], float)
        tmin = np.asarray(weather["tmin"], float)
        rad = np.asarray(weather["rad_mj"], float)
        assert len(tmax) >= n and len(tmin) >= n and len(rad) >= n, "weather too short"
        return tmax[:n], tmin[:n], rad[:n]
    d = np.arange(n)
    # mild seasonal drift; tropical lowland rice (IRRI-like), no strong seasonality
    rad = 19.0 + 3.0 * np.sin(2 * np.pi * (d - 20) / 120.0)       # ~16-22 MJ/m2/d
    tmax = 31.0 + 1.5 * np.sin(2 * np.pi * (d - 30) / 120.0)
    tmin = 23.0 + 1.0 * np.sin(2 * np.pi * (d - 30) / 120.0)
    return tmax, tmin, rad


# ---------------------------------------------------------------------------
# Public API: per-organ biomass on a time grid + driver builder
# ---------------------------------------------------------------------------
def organ_biomass_oryza(t: np.ndarray, p: OryzaParams | None = None,
                        weather: dict | None = None, scale_to_anchor: bool = True):
    """Per-organ dry biomass [kg/hill] on grid `t` (root/stem/leaf/grain).

    Runs the ORYZA Level-1 loop and interpolates onto `t`.  If `scale_to_anchor`,
    multiply all organs by a single factor so final shoot matches the IR72 field
    anchor (potential -> attainable; preserves the mechanistic *shape* and HI)."""
    p = p or OryzaParams()
    sim = simulate_oryza(p, weather)
    td = sim["t"]
    kg_hill = HA_PER_HILL
    organs = {"root": sim["wrt"], "leaf": sim["wlv"], "stem": sim["wst"], "grain": sim["wso"]}
    if scale_to_anchor:
        shoot_final = sim["wlv"][-1] + sim["wst"][-1] + sim["wso"][-1]    # kg/ha
        target = p.shoot_anchor_g_m2 * 10.0                              # g/m2 -> kg/ha
        f = target / shoot_final if shoot_final > 0 else 1.0
        organs = {k: v * f for k, v in organs.items()}
    out = {k: np.maximum(np.interp(t, td, v * kg_hill), 1e-9) for k, v in organs.items()}
    # leaf senescence loss RATE [1/day] (NOT a mass -> no kg_hill / anchor scaling). The
    # PFAS leaf ODE adds -leaf_death_rate*C so the dead leaf carries its PFAS away
    # (D/M_leaf = drlv exactly), cancelling the spurious senescence concentration.
    out["leaf_death_rate"] = np.maximum(np.interp(t, td, sim["drlv"]), 0.0)
    return out


def M_matrix_oryza(t: np.ndarray, p: OryzaParams | None = None, weather: dict | None = None):
    """Mass matrix [len(t), 4] = columns [root, stem, leaf, grain] [kg/hill]."""
    b = organ_biomass_oryza(t, p, weather)
    return np.column_stack([b["root"], b["stem"], b["leaf"], b["grain"]])


def oryza_drivers(congener: str = "PFOA", Cwo: float = 1.0, season: float = 120.0,
                  n_t: int = 241, p: OryzaParams | None = None, weather: dict | None = None,
                  Qtp=None):
    """Build a `drivers` dict for model_api.simulate(drivers=...) with ORYZA biomass.

    Cwo is held constant (so conc==BAF at Cwo=1); Q_TP defaults to the measured
    forcing_rice transpiration; M is the mechanistic ORYZA biomass."""
    import forcing_rice as fr
    t = np.linspace(0.0, season, n_t)
    b = organ_biomass_oryza(t, p, weather)
    M = np.column_stack([b["root"], b["stem"], b["leaf"], b["grain"]])
    Qtp = fr.Q_TP(t, season) if Qtp is None else np.asarray(Qtp, float)
    return dict(t=t, Cwo=np.full_like(t, float(Cwo)), Qtp=Qtp, M=M,
                leaf_loss=b["leaf_death_rate"])


if __name__ == "__main__":
    p = OryzaParams()
    sim = simulate_oryza(p)
    i_flower = int(np.argmax(sim["dvs"] >= 1.0))
    i_mat = int(np.argmax(sim["dvs"] >= 2.0)) or (len(sim["t"]) - 1)
    shoot = sim["wlv"] + sim["wst"] + sim["wso"]
    print("ORYZA2000 Level-1 (Python) — IR72 potential run")
    print(f"  flowering (DVS=1) ~ day {sim['t'][i_flower]:.0f};  maturity (DVS=2) ~ day {sim['t'][i_mat]:.0f}")
    print(f"  LAI max           = {sim['lai'].max():.2f}  @ day {sim['t'][int(np.argmax(sim['lai']))]:.0f}")
    print(f"  peak CGR          = {sim['cgr'].max():.0f} kg DM/ha/d")
    print(f"  final shoot       = {shoot[-1]:.0f} kg/ha ({shoot[-1]/10:.0f} g/m^2; anchor 1740)")
    print(f"  final grain (WSO) = {sim['wso'][-1]:.0f} kg/ha;  HI = {sim['wso'][-1]/shoot[-1]:.2f}")
    print(f"  final root (WRT)  = {sim['wrt'][-1]:.0f} kg/ha;  root:shoot = {sim['wrt'][-1]/shoot[-1]:.2f}")
    # anchored per-hill biomass at a few days
    t = np.linspace(0, p.season, 481)
    b = organ_biomass_oryza(t, p)
    print("\nanchored biomass [kg/hill]:")
    for dday in (20, 45, 70, 95, 120):
        i = int(np.argmin(abs(t - dday)))
        print(f"   day {dday:3d}: root={b['root'][i]:.4f} stem={b['stem'][i]:.4f} "
              f"leaf={b['leaf'][i]:.4f} grain={b['grain'][i]:.4f}")
