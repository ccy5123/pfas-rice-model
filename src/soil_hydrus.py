"""
Soil side via REAL HYDRUS-1D (Method A, loose one-way coupling)
==============================================================

This wires the genuine HYDRUS-1D solver (built from ``external/hydrus_source``
and driven through ``phydrus``) into the plant module's ``PlantInputs``.  It is
the soil half of Method A:

    HYDRUS-1D (Richards + advection-dispersion + Kd sorption + root water uptake)
        --->  C_w^o(t)  [pore-water PFAS at the root zone]
        --->  Q_TP(t)   [actual root water uptake / transpiration]
        --->  PlantInputs  --->  the 4-compartment plant ODE

Why a real soil model (vs. the constant ``Cwo`` placeholder or the analytic
``soil_paddy`` Freundlich inversion): the pore-water trajectory is strongly
**congener-dependent**.  Weakly-sorbed short chains (small Kd) leach quickly
under the paddy water regime, so C_w^o(t) falls during flooding and rebounds as
the soil drains at the end of the season; strongly-sorbed long chains (large Kd)
are buffered and stay nearly flat.  A constant Cwo cannot represent either.

Sorption is a LINEAR Kd isotherm (HYDRUS ``ks`` with ``beta=1``).  A Freundlich
exponent <1 makes ds/dc -> inf as c -> 0 at the clean-irrigation boundary and the
HYDRUS solute solver fails to converge there; linear Kd is the robust, standard
choice and still gives the full congener-resolved retardation R = 1 + rho*Kd/theta.
Kd is taken from the C3 Koc QSPR (``literature_params.koc``) x f_oc.

Requires the compiled ``hydrus`` executable:
    git submodule update --init external/hydrus_source
    cp external/hydrus_source/makefile external/hydrus_source/source/
    (cd external/hydrus_source/source && make)     # needs gfortran

Units inside HYDRUS: length cm, time day, mass mg (so C_w is mg/cm^3).  The
plant ODE is linear in Cwo, so only the SHAPE and congener-to-congener contrast
of C_w^o(t) matter; ``inputs_from_hydrus`` normalises the series to a chosen
season-mean reference (default 1.0 ug/L) so it is directly comparable to the
constant-Cwo baseline (same mean exposure, realistic temporal structure).
"""
from __future__ import annotations

import os
import re
import shutil
import tempfile
from dataclasses import dataclass

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
HYDRUS_EXE = os.path.join(_ROOT, "external", "hydrus_source", "source", "hydrus")

import literature_params as lp  # noqa: E402  (C3 Koc QSPR)


# ---------------------------------------------------------------------------
# availability + per-congener sorption
# ---------------------------------------------------------------------------
def hydrus_available(exe: str = HYDRUS_EXE) -> bool:
    """True if the compiled HYDRUS-1D executable is present and runnable."""
    return os.path.isfile(exe) and os.access(exe, os.X_OK)


def _has_fortran_source(source_dir: str) -> bool:
    return os.path.isdir(source_dir) and any(
        f.lower().endswith(".for") for f in os.listdir(source_dir))


def build_hydrus_engine(exe: str = HYDRUS_EXE):
    """Best-effort, in-place build of the HYDRUS-1D engine (for Streamlit Cloud / a
    fresh clone). Fetches the source (git submodule, else clone), copies the makefile
    into ``source/`` and runs ``make`` with gfortran. Returns ``(ok, log_lines)``.

    Needs ``gfortran`` + ``make`` on PATH — on Streamlit Cloud add a ``packages.txt``
    with ``gfortran`` and ``make``; ``phydrus`` must be installed (requirements.txt)."""
    import subprocess
    log: list[str] = []

    def _say(m):
        log.append(m)

    if hydrus_available(exe):
        _say("HYDRUS-1D engine already built.")
        return True, log
    if shutil.which("gfortran") is None or shutil.which("make") is None:
        _say("gfortran/make not found on PATH. On Streamlit Cloud add a `packages.txt` "
             "containing `gfortran` and `make`, then redeploy.")
        return False, log

    src_dir = os.path.join(_ROOT, "external", "hydrus_source")
    source = os.path.join(src_dir, "source")
    if not _has_fortran_source(source):
        _say("fetching HYDRUS-1D source…")
        subprocess.run(["git", "submodule", "update", "--init", "external/hydrus_source"],
                       cwd=_ROOT, capture_output=True, text=True)
        if not _has_fortran_source(source):                 # submodule not available -> clone
            try:
                if os.path.isdir(src_dir) and not os.listdir(src_dir):
                    os.rmdir(src_dir)
                subprocess.run(["git", "clone", "--depth", "1",
                                "https://github.com/phydrus/source_code", src_dir],
                               capture_output=True, text=True, check=False)
            except Exception as e:                          # noqa: BLE001
                _say(f"clone failed: {e}")
    if not _has_fortran_source(source):
        _say("could not obtain the HYDRUS-1D FORTRAN source (no network / submodule).")
        return False, log

    mk = os.path.join(src_dir, "makefile")                   # top-level makefile -> source/
    if os.path.isfile(mk):
        shutil.copy(mk, os.path.join(source, "makefile"))
    _say("compiling with gfortran (this takes ~1 minute)…")
    r = subprocess.run(["make"], cwd=source, capture_output=True, text=True)
    if hydrus_available(exe):
        _say("✓ HYDRUS-1D engine built successfully.")
        return True, log
    _say("build did not produce the executable. make output (tail):")
    _say((r.stderr or r.stdout or "")[-1500:])
    return False, log


def paddy_kd(n_C: int, group: str = "PFCA", f_oc: float | None = None) -> float:
    """Linear soil distribution coefficient Kd [L/kg] for a congener.

    Kd = Koc(chain length, head group) * f_oc, with Koc from the C3 QSPR
    (Higgins & Luthy slope, Milinovic PFOA anchor) in ``literature_params``.
    PFCA C_n has n-1 perfluorinated carbons; PFSA C_n has n.
    """
    head = "sulfonate" if group.upper() == "PFSA" else "carboxylate"
    n_pfc = n_C - 1 if group.upper() == "PFCA" else n_C
    return lp.koc_to_KF(lp.koc(n_pfc, head), FOC if f_oc is None else f_oc)


# ---------------------------------------------------------------------------
# robust parsers for the HYDRUS ASCII output (phydrus 0.2 reader is brittle)
# ---------------------------------------------------------------------------
def _parse_obs_node(path: str):
    """Return (t, Conc, theta) at the FIRST observation node from OBS_NODE.OUT."""
    lines = open(path).read().splitlines()
    hdr = next(i for i, l in enumerate(lines) if l.strip().startswith("time"))
    rows = []
    for l in lines[hdr + 1:]:
        s = l.split()
        if len(s) >= 5 and re.match(r"^-?\d", s[0]):
            rows.append([float(x) for x in s])
    a = np.array(rows)
    # columns: time, then (h, theta, Temp, Conc) per node
    return a[:, 0], a[:, 4], a[:, 2]


def _parse_tlevel(path: str):
    """Return (t, vRoot) — actual root water uptake [cm/day] — from T_LEVEL.OUT."""
    lines = open(path).read().splitlines()
    hdr = next(i for i, l in enumerate(lines) if l.strip().startswith("Time"))
    names = lines[hdr].split()
    j_t, j_v = names.index("Time"), names.index("vRoot")
    rows = []
    for l in lines[hdr + 2:]:
        s = l.split()
        if len(s) >= len(names) and re.match(r"^-?\d", s[0]):
            rows.append([float(x) for x in s[:len(names)]])
    a = np.array(rows)
    return a[:, j_t], a[:, j_v]


# ---------------------------------------------------------------------------
# the HYDRUS paddy run
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# the paddy soil: the Tsukuba Andosol of the Yamazaki field (PREREG 9.2, 9.8)
# ---------------------------------------------------------------------------
# Every number is a row of docs/literature_db/sourced_inputs.csv; a sourced range runs at its
# midpoint (PREREG 5.4) and its ends are the band. Retention: the two original Tsukuba Andisol
# layers under the NIAES paddy, undisturbed cores (Tournebize et al. 2006, Table 3; theta_r not
# printed there) -- the end sets, and their scalar midpoint for the run.
ANDOSOL_SETS = {
    "ash1": dict(ths=0.709, alpha=0.0050, n=1.068, Ks=11.1, rho_sub=0.70),
    "ash2": dict(ths=0.717, alpha=0.0098, n=1.156, Ks=93.0, rho_sub=0.57),
}
SOIL_RANGES = {
    "thr": (0.0, 0.030),             # Rudiyanto 2013 fixed 0; Tokumoto 2010 fitted 0.030
    "percolation": (0.41, 0.67),     # cm/d, Kondo et al. 2019 (Ibaraki Andosol paddy)
    "pond": (2.0, 3.0),              # cm, shallow water after establishment (Ibaraki calendar)
    "disper": (0.2, 1.0),            # cm, Tokumoto 2005 columns; Tournebize 2006 field tracer
    "drain_start": (35.0, 40.0),     # DAT, midsummer drainage (Ibaraki calendar)
    "drain_days": (15.0, 20.0),      # d
    "f_oc": (0.037, 0.038),          # Yamazaki SI Table S5 (total C 3.7-3.8 %, 0-15 cm)
}
RHO_PLOUGH = 0.61                    # g/cm3, Tsukuba Andosol paddy topsoil (Koyama & Hayashi 2017)
PLOUGH_DEPTH = 15.0                  # cm, "耕深15cmの確保" (Ibaraki calendar)
PAN_THICKNESS = 2.0                  # cm, Tournebize 2006; only Ks/thickness is identified (9.8)
L_MUALEM = 0.5                       # Mualem's general value; no paddy value exists
REFLOOD_TO_HEADING_DAYS = 30.0       # water kept in to ~30 d after heading (Ibaraki calendar)
FOC = 0.5 * sum(SOIL_RANGES["f_oc"])


def andosol_soil(retention="mid", **ends):
    """The soil dict of the scenario: `retention` "ash1" / "ash2" (an end) or "mid" (the scalar
    midpoint of the two), every SOIL_RANGES input at its midpoint unless given in `ends`."""
    if retention == "mid":
        a, b = ANDOSOL_SETS["ash1"], ANDOSOL_SETS["ash2"]
        ret = {k: 0.5 * (a[k] + b[k]) for k in a}
    else:
        ret = dict(ANDOSOL_SETS[retention])
    soil = dict(ret, **{k: 0.5 * (lo + hi) for k, (lo, hi) in SOIL_RANGES.items()})
    soil.update(ends)
    return soil


def _vg_K(h, ths, thr, alpha, n, Ks, l=L_MUALEM):
    """Mualem-van Genuchten conductivity [cm/d] at pressure head h [cm] (h <= 0)."""
    m = 1.0 - 1.0 / n
    Se = (1.0 + np.abs(alpha * np.minimum(h, 0.0)) ** n) ** (-m)
    return Ks * Se ** l * (1.0 - (1.0 - Se ** (1.0 / m)) ** m) ** 2


def pan_ks(soil: dict, return_head: bool = False):
    """The plough pan's Ks [cm/d] that makes steady percolation under the ponded depth equal the
    sourced percolation q (PREREG 9.8). Steady Darcy-Buckingham flow, z up, q down:
    q = K(h) (dh/dz + 1). Below the pan the subsoil is at unit gradient, K_sub(h_s) = q; the plough
    layer above is saturated, so the pan's top sits at h_top = pond + L - q L / Ks. Across the pan,
    with its own K(h) (its lower part desaturates under the subsoil's suction), the thickness is
        d = integral from h_s to h_top of dh / (q / K_pan(h) - 1),
    and Ks_pan is the root of that condition. (The retention is the soil's; only Ks differs.)"""
    from scipy.optimize import brentq
    q = float(soil["percolation"])
    L, d = PLOUGH_DEPTH, PAN_THICKNESS
    vg = (soil["ths"], soil["thr"], soil["alpha"], soil["n"])
    if q >= soil["Ks"]:
        raise ValueError("the percolation exceeds the subsoil's Ks: no unsaturated steady state")
    h_s = brentq(lambda h: _vg_K(h, *vg, soil["Ks"]) - q, -1e6, -1e-9)
    h_top = soil["pond"] + L - q * L / soil["Ks"]
    hs = np.linspace(h_s, 0.0, 4001)

    def thickness(ks):
        if ks >= q:
            return np.inf                              # a pan at least as conductive as q: no resistance
        g = 1.0 / (q / _vg_K(hs, *vg, ks) - 1.0)       # dz/dh, unsaturated part (K < q throughout)
        unsat = float(np.sum(0.5 * (g[1:] + g[:-1]) * np.diff(hs)))
        return unsat + h_top / (q / ks - 1.0)          # + the saturated part, h 0 -> h_top

    ks = brentq(lambda k: thickness(k) - d, q * 1e-6, q * (1.0 - 1e-9), xtol=1e-14, rtol=1e-12)
    return (ks, h_s) if return_head else ks


def _pan_profile(soil, ks, h_s, z_rel):
    """The steady head across the pan at heights z_rel [cm above its bottom] (`pan_ks`)."""
    q = float(soil["percolation"])
    vg = (soil["ths"], soil["thr"], soil["alpha"], soil["n"])
    hs = np.linspace(h_s, 0.0, 4001)
    g = 1.0 / (q / _vg_K(hs, *vg, ks) - 1.0)
    zu = np.concatenate([[0.0], np.cumsum(0.5 * (g[1:] + g[:-1]) * np.diff(hs))])
    zr = np.asarray(z_rel, float)
    # saturated above the head-zero height: h rises linearly at dh/dz = q/ks - 1
    return np.where(zr <= zu[-1], np.interp(zr, zu, hs), (zr - zu[-1]) * (q / ks - 1.0))


@dataclass
class PaddyResult:
    t: np.ndarray        # output times [day]
    Cw: np.ndarray       # root-zone pore water: the root-weighted mean over the profile [input units]
    theta: np.ndarray    # root-weighted mean water content [-]
    vroot: np.ndarray    # actual transpiration / root water uptake [cm/day]
    Kd: float            # linear sorption coeff used [L/kg]
    soil: dict | None = None             # the soil inputs of the run
    plough_total: np.ndarray | None = None  # 0-15 cm total soil conc [c-units x L/kg]: Kd c + theta c / rho
    ks_pan: float | None = None


def _parse_nod_inf(path: str):
    """{time: (z, h, theta, conc)} from NOD_INF.OUT (the profile at each print time)."""
    out, t, rows = {}, None, []
    for line in open(path).read().splitlines():
        s = line.split()
        if line.strip().startswith("Time:"):
            if t is not None and rows:
                out[t] = np.array(rows)
            t, rows = float(s[1]), []
        elif t is not None and len(s) >= 12 and re.match(r"^\d+$", s[0]):
            # Node Depth Head Moisture K C Flux Sink Kappa v/KsTop Temp Conc Sorb
            rows.append([float(s[1]), float(s[2]), float(s[3]), float(s[11])])
    if t is not None and rows:
        out[t] = np.array(rows)
    return {k: (v[:, 0], v[:, 1], v[:, 2], v[:, 3]) for k, v in out.items()}


def _graded_profile(ps, depth, dx):
    """The node grid: `dx` [cm] away from the pan, dx/5 within 2.5 cm of it, and dx/100 across the
    pan's base and the top of the subsoil, where the steady head climbs from the subsoil's suction to
    saturation within 0.04-0.37 cm (`pan_ks`). Numerical (converged, PREREG 9.2)."""
    zb = PLOUGH_DEPTH + PAN_THICKNESS
    edges = [(0.0, PLOUGH_DEPTH - 1.0, dx), (PLOUGH_DEPTH - 1.0, zb - 0.5, dx / 5),
             (zb - 0.5, zb + 0.5, dx / 100), (zb + 0.5, zb + 3.0, dx / 5), (zb + 3.0, abs(depth), dx)]
    x = np.unique(np.round(np.concatenate([np.arange(a, b, h) for a, b, h in edges]
                                          + [[abs(depth)]]), 6))
    if len(x) > 1001:
        raise ValueError(f"{len(x)} nodes: HYDRUS-1D allows 1001 (NumNPD); coarsen dx")
    base = ps.create_profile(top=0, bot=-1.0, dx=1.0, h=0.0, mat=3, lay=1, conc=0.0)
    row = base.iloc[[0]]
    prof = row.loc[row.index.repeat(len(x))].reset_index(drop=True)
    prof.index = prof.index + 1
    prof["x"] = -x
    return prof


def run_paddy_hydrus(Kd: float, season: float = 120.0, *, soil: dict | None = None,
                     depth: float = 100.0, dx: float = 0.5,
                     tpot_fn=None, evap_fn=None, rain_fn=None, heading: float | None = None,
                     source: str = "soil", c_init: float = 1.0, c_sub: float | None = None,
                     c_irrigation: float = 0.0, c_rain: float = 0.0,
                     exe: str = HYDRUS_EXE, workspace: str | None = None,
                     keep: bool = False) -> PaddyResult:
    """One season of the Tsukuba Andosol paddy in HYDRUS-1D for a linear Kd [L/kg] (PREREG 9.2).

    Profile: the 15-cm plough layer (bulk density 0.61), a 2-cm plough pan whose Ks makes steady
    percolation under the ponded depth equal the sourced percolation (`pan_ks`), the subsoil to
    `depth` with free drainage; retention from the Tsukuba Andisol cores (`andosol_soil`). The top
    is an atmospheric boundary with a surface layer, so water ponds. Water calendar: flooded from
    transplanting, midsummer drainage (`drain_start`, `drain_days`), flooded again to 30 d after
    `heading` (or to the end of the season if no heading is given). While flooded the field takes
    ET + percolation, rain first and irrigation for the rest (rain beyond it leaves over the bund);
    while drained the input is the rain. `tpot_fn`, `evap_fn`, `rain_fn`: day -> mm/day.

    source "soil": the soil holds the PFAS (pore water `c_init` everywhere) and the water is clean
    (the generic contaminated-field case). source "irrigation": pore water `c_init` in the plough
    layer and `c_sub` below (default `c_init`), and the inflow carries `c_irrigation` and `c_rain`
    (flux-weighted each day) -- the Yamazaki field (PREREG 9.3). The exposure is the root-weighted
    mean pore water over the profile; the root profile (linear to 0 at 30 cm) has no source in the
    repo (`soil_root_beta`, flagged).
    """
    import phydrus as ps
    try:
        ps.set_log_level("ERROR")
    except Exception:
        pass
    if not hydrus_available(exe):
        raise FileNotFoundError(
            f"HYDRUS executable not found/built at {exe}. Build it with gfortran "
            "(see module docstring).")
    if source not in ("soil", "irrigation"):
        raise ValueError(f"source must be 'soil' or 'irrigation', got {source!r}")
    soil = dict(andosol_soil() if soil is None else soil)
    kpan, h_sub = pan_ks(soil, return_head=True)
    soil["Ks_pan"] = kpan

    if tpot_fn is None or evap_fn is None:
        raise ValueError("pass the crop's transpiration and evaporation (tpot_fn, evap_fn; "
                         "mm/day): the scenario has no default for either (PREREG 9.2)")
    days = np.arange(1, int(np.ceil(season)) + 1, dtype=float)
    mm = lambda fn: np.maximum(np.asarray(fn(days), float), 0.0) / 10.0   # mm/d -> cm/d
    tpot, evap = mm(tpot_fn), mm(evap_fn)
    rain = mm(rain_fn) if rain_fn is not None else np.zeros_like(days)    # generic: no weather
    d0, d1 = soil["drain_start"], soil["drain_start"] + soil["drain_days"]
    end_flood = season if heading is None else min(season, heading + REFLOOD_TO_HEADING_DAYS)
    flooded = ((days <= d0) | (days > d1)) & (days <= end_flood)
    demand = evap + tpot + soil["percolation"]
    rain_in = np.where(flooded, np.minimum(rain, demand), rain)
    irr = np.where(flooded, demand - rain_in, 0.0)
    prec = rain_in + irr
    with np.errstate(invalid="ignore", divide="ignore"):
        ctop = np.where(prec > 0, (rain_in * c_rain + irr * c_irrigation) / prec, 0.0)
    if source == "soil":
        ctop = np.zeros_like(ctop)

    ws = workspace or tempfile.mkdtemp(prefix="hydrus_paddy_")
    try:
        ml = ps.Model(exe_name=exe, ws_name=ws, name="paddy", mass_units="mg",
                      length_unit="cm", time_unit="days", print_screen=False)
        prints = days[days <= season]
        if len(prints) % 6 == 0:          # phydrus writes an empty line after a full last row of 6,
            prints = np.concatenate([[0.5], prints])   # which HYDRUS misreads: one more print time
        ml.add_time_info(tinit=0, tmax=season, dt=1e-4, dtmin=1e-7, dtmax=0.25,
                         print_times=True, print_array=prints)
        ml.add_waterflow(model=0, top_bc=2, bot_bc=4, maxit=20, tolth=1e-3, tolh=1)
        ml.add_solute_transport(model=0, top_bc=-1, bot_bc=0, epsi=1.0, lupw=True, pecr=2)

        m = ml.get_empty_material_df(n=3)
        for mat, ks, rho in ((1, soil["Ks"], RHO_PLOUGH), (2, kpan, soil["rho_sub"]),
                             (3, soil["Ks"], soil["rho_sub"])):
            m.loc[mat, ("water", "thr")] = soil["thr"]
            m.loc[mat, ("water", "ths")] = soil["ths"]
            m.loc[mat, ("water", "Alfa")] = soil["alpha"]
            m.loc[mat, ("water", "n")] = soil["n"]
            m.loc[mat, ("water", "Ks")] = ks
            m.loc[mat, ("water", "l")] = L_MUALEM
            m.loc[mat, ("solute", "bulk.d")] = rho
            m.loc[mat, ("solute", "DisperL")] = soil["disper"]
            m.loc[mat, ("solute", "frac")] = 1.0
            m.loc[mat, ("solute", "mobile_wc")] = 0.0
        ml.add_material(m)

        zb = PLOUGH_DEPTH + PAN_THICKNESS
        prof = _graded_profile(ps, depth, dx)
        z = prof["x"].values                                   # 0 .. -depth
        prof.loc[z >= -PLOUGH_DEPTH, "Mat"] = 1
        prof.loc[(z < -PLOUGH_DEPTH) & (z >= -zb), "Mat"] = 2
        # puddled at transplanting, at the steady state `pan_ks` solves: the plough layer saturated
        # and hydrostatic under the ponded depth, the head falling across the pan, the subsoil at
        # its unit-gradient head
        prof["h"] = np.where(z >= -PLOUGH_DEPTH, soil["pond"] - z,
                             np.where(z >= -zb, _pan_profile(soil, kpan, h_sub, zb + z), h_sub))
        c_below = c_init if c_sub is None else c_sub
        prof["Conc"] = np.where(z >= -PLOUGH_DEPTH, c_init, c_below)
        beta = np.clip(1.0 + z / 30.0, 0, None)                 # roots: surface -> 0 at -30 cm
        beta[z < -30] = 0.0
        prof["Beta"] = beta
        ml.add_profile(prof)
        # no reduction above the wilting point (PREREG 9.10): uptake at the potential rate from
        # any ponded head (p0 = popt above it) down to -15000 cm, where it stops
        wilt = -15000.0
        ml.add_root_uptake(model=0, p0=1.0e4, poptm=[1.0e4 - 1.0] * 3, p2h=wilt + 2.0,
                           p2l=wilt + 1.0, p3=wilt, r2h=0.5, r2l=0.1)

        import pandas as pd
        atm = pd.DataFrame({"tAtm": days, "Prec": prec, "rSoil": evap, "rRoot": tpot,
                            "hCritA": 1e5, "rB": 0.0, "hB": 0.0, "ht": 0.0, "tTop": 0.0,
                            "tBot": 0.0, "Ampl": 0.0, "cTop": ctop, "cBot": 0.0})
        ml.add_atmospheric_bc(atm, hcrits=1e4, tatm=0.0, prec=0.0, rsoil=0.0, rroot=0.0,
                              hcrita=1e5, rb=0.0, hb=0.0, ht=0.0, ttop=0.0, tbot=0.0, ampl=0.0)

        sdf = ml.get_empty_solute_df()
        sdf.loc[:, "ks"] = float(Kd)     # linear isotherm: s = Kd * c
        sdf.loc[:, "beta"] = 1.0
        sdf.loc[:, "nu"] = 0.0
        ml.add_solute(sdf, difw=0.43, difg=0.0, top_conc=0.0, bot_conc=0.0)

        ml.write_input()
        ml.simulate()

        err = os.path.join(ws, "Error.msg")
        if os.path.isfile(err) and os.path.getsize(err) > 0:
            raise RuntimeError("HYDRUS did not converge: " + open(err).read().strip())

        prof_t = _parse_nod_inf(os.path.join(ws, "NOD_INF.OUT"))
        ts = np.array(sorted(prof_t))
        w = None
        Cw, th, plough = [], [], []
        for tt in ts:
            zz, hh, tq, cc = prof_t[tt]
            if w is None:                 # root density x each node's control volume
                edges = np.concatenate([[zz[0]], 0.5 * (zz[1:] + zz[:-1]), [zz[-1]]])
                w = np.clip(1.0 + zz / 30.0, 0, None) * np.abs(np.diff(edges))
                w = w / w.sum()
                top = zz >= -PLOUGH_DEPTH
                vt = np.abs(np.diff(edges))[top] / np.abs(np.diff(edges))[top].sum()
            Cw.append(float(np.sum(w * cc)))
            th.append(float(np.sum(w * tq)))
            plough.append(float(np.sum(vt * (Kd * cc[top] + tq[top] * cc[top] / RHO_PLOUGH))))
        tv, vroot = _parse_tlevel(os.path.join(ws, "T_LEVEL.OUT"))
        vroot = np.interp(ts, tv, vroot)
        return PaddyResult(t=ts, Cw=np.array(Cw), theta=np.array(th), vroot=vroot, Kd=float(Kd),
                           soil=soil, plough_total=np.array(plough), ks_pan=float(kpan))
    finally:
        if not keep:
            shutil.rmtree(ws, ignore_errors=True)


def tsukuba_rain_mm(doy):
    """JMA Tsukuba 2020 daily precipitation [mm] on day(s) of year `doy`
    (raw_si/jma_tsukuba_2020_precip.csv; the Yamazaki field's rain, PREREG 9.3)."""
    import csv as _csv
    path = os.path.join(_ROOT, "docs", "literature_db", "raw_si", "jma_tsukuba_2020_precip.csv")
    rows = list(_csv.DictReader(l for l in open(path) if not l.startswith("#")))
    d = np.array([int(r["doy"]) for r in rows], float)
    p = np.array([float(r["prec_mm"]) for r in rows])
    return np.interp(np.floor(np.asarray(doy, float)), d, p, left=0.0, right=0.0)


def _s16(compound):
    """Yamazaki SI Table S16 row for `compound` (raw_si/yamazaki2023_si_s16_soil_rain.csv);
    "<x" (below the MQL x) is read as x/2 and flagged."""
    import csv as _csv
    path = os.path.join(_ROOT, "docs", "literature_db", "raw_si", "yamazaki2023_si_s16_soil_rain.csv")
    for r in _csv.DictReader(l for l in open(path) if not l.startswith("#")):
        if r["compound"] == compound:
            val = lambda v: (0.5 * float(v[1:]), True) if v.startswith("<") else (float(v), False)
            return {k: val(v) for k, v in r.items() if k != "compound"}
    raise KeyError(f"{compound} not in Yamazaki Table S16")


def yamazaki_source_kw(compound, Kd, soil=None, sub="mid"):
    """`run_paddy_hydrus` keywords for the irrigation-source scenario of the Yamazaki field
    (PREREG 9.3), in pore-water units ng/L: irrigation water and mean rain (Table S16), and the
    initial plough-layer pore water from the pre-season soil (water extract + methanol extract
    after it, pg/g dw) at equilibrium, C_T = c (Kd + theta_s / rho). The subsoil below 15 cm was
    not sampled: `sub` "same" (= the plough layer), "zero", or "mid" (half, the run)."""
    soil = andosol_soil() if soil is None else soil
    r = _s16(compound)
    rain = [r[f"rain_R{i}"][0] for i in range(1, 6)]
    C_T = r["pre_water"][0] + r["pre_wtmeoh"][0]                 # pg/g dw = ng/kg
    c0 = C_T / (float(Kd) + soil["ths"] / RHO_PLOUGH)            # ng/L
    c_sub = {"same": c0, "zero": 0.0, "mid": 0.5 * c0}[sub]
    flags = sorted(k for k, (_, f) in r.items() if f)
    return dict(source="irrigation", c_init=c0, c_sub=c_sub,
                c_irrigation=r["irrigation_ng_L"][0], c_rain=float(np.mean(rain))), flags


def post_season_soil(compound):
    """Measured 0-15 cm soil before and after the season (water + methanol-after-water, pg/g dw)."""
    r = _s16(compound)
    return (r["pre_water"][0] + r["pre_wtmeoh"][0], r["post_water"][0] + r["post_wtmeoh"][0])


def inputs_from_hydrus(congener_n_C: int = None, group: str = "PFCA", *, season: float | None = None,
                       Cwo_ref: float = 1.0, f_oc: float | None = None,
                       qtp_from_hydrus: bool = True, area_per_hill_m2: float | None = None,
                       n_t: int = 241, biomass: str = "oryza", Kd: float | None = None,
                       scenario: str = "yamazaki2020", **run_kw):
    """Build a :class:`PlantInputs` whose ``Cwo`` AND ``Q_TP`` come from a real
    HYDRUS-1D paddy run for the congener.

    scenario "yamazaki2020" (default): the crop side is the calibration site-year
        (`site_scenario`): the soil run is driven by its ORYZA2000 potential transpiration
        over its season (transplanting to harvest), ``M`` is its fresh organ mass and the
        area per hill is 1/density; ``n_t``, ``biomass`` and ``area_per_hill_m2`` do not
        apply, and an explicit ``season`` other than the site's is an error.
    scenario "generic": ``forcing_rice`` transpiration over ``season`` days (default 120)
        and the selected ``biomass`` driver (alternative drivers).

    The HYDRUS pore-water series is normalised to season-mean ``Cwo_ref`` so the
    average exposure matches a constant-Cwo run -- the difference is purely the
    realistic temporal structure (leaching during flooding, rebound on drainage).

    qtp_from_hydrus : if True (default), the soil run is driven by the MEASURED
        transpiration (``forcing_rice.transpiration_mm_d``) and Q_TP(t) is the
        HYDRUS actual root water uptake ``vRoot`` [cm/day] converted to L/day/hill
        via ``area_per_hill_m2`` -- so it equals the measured forcing when the
        paddy is unstressed and only diverges under soil-water limitation. If
        False, the measured ``forcing_rice.Q_TP`` is used directly (soil run uses
        the default gaussian potential transpiration).
    Kd : pass a linear sorption coefficient [L/kg] directly instead of deriving it
        from the PFAS chain-length Koc QSPR. This is what a NEUTRAL organic needs --
        `paddy_kd` is a per-CF2 group contribution and has no meaning for one -- and
        it is also the way in for a measured Kd. When given, ``congener_n_C``/``group``
        are unused (and ``congener_n_C`` may be None).
    Returns (PlantInputs, PaddyResult).
    """
    from pfas_rice_plant_module_4pool_surf import PlantInputs
    import forcing_rice as fr
    import organ_mass as om

    site = None
    if scenario == "yamazaki2020":
        import site_scenario as ss
        site = ss.yamazaki2020_mean_run()
        ss.check_season(season, site["t"], "inputs_from_hydrus")
        season = float(site["t"][-1])
        area_per_hill_m2 = 1.0 / ss.NH
    elif scenario != "generic":
        raise ValueError(f"unknown scenario {scenario!r}")
    else:
        season = 120.0 if season is None else float(season)
    if area_per_hill_m2 is None:
        area_per_hill_m2 = fr.AREA_PER_HILL_M2

    if Kd is None:
        if congener_n_C is None:
            raise ValueError("give either congener_n_C (PFAS) or an explicit Kd [L/kg]")
        Kd = paddy_kd(congener_n_C, group, f_oc)
    Kd = float(Kd)
    # drive the soil run with the crop's own water terms (mm/day): its potential transpiration and
    # evaporation, and at the site its rain and heading (PREREG 9.2)
    if site is not None:
        tp = ss._doy(ss.YAMAZAKI_TRANSPLANT)
        run_kw.setdefault("tpot_fn", lambda d: np.interp(d, site["t"], site["Tp"]))
        run_kw.setdefault("evap_fn", lambda d: np.interp(d, site["t"], site["Ep"]))
        run_kw.setdefault("rain_fn", lambda d: tsukuba_rain_mm(tp + d - 1.0))
        run_kw.setdefault("heading", float(ss.YAMAZAKI_HEADING_MEAN_DOY - tp))
    else:
        run_kw.setdefault("tpot_fn", lambda d: fr.transpiration_mm_d(d, season))
        run_kw.setdefault("evap_fn", lambda d: fr.evaporation_mm_d(d, season))
    res = run_paddy_hydrus(Kd, season=season, **run_kw)

    if site is not None:
        t = np.asarray(site["t"], float)
        b = {k: site[k] for k in ("root", "stem", "leaf", "grain")}
        leaf_loss = site["leaf_death_rate"]
        q_direct = site["Qtp"]
        dvs = np.asarray(site["dvs"], float)
    else:
        t = np.linspace(0.0, season, n_t)
        if biomass == "oryza":
            import oryza_growth as og
            b = og.organ_biomass_oryza(t, p=og.OryzaParams(season=season))
        else:
            import growth_rice as gr
            b = gr.organ_biomass(t, season)
        leaf_loss = b.get("leaf_death_rate")          # only the senescing ORYZA driver supplies it
        q_direct = fr.Q_TP(t, season)
        dvs = np.asarray(b["dvs"], float)
    Cw_h = np.interp(t, res.t, res.Cw)
    Cw_h = Cw_h * (Cwo_ref / np.mean(Cw_h))          # normalise to mean = Cwo_ref
    M = om.organ_matrix(b)                            # crop DRY matter -> ODE FRESH mass

    if qtp_from_hydrus:
        # vRoot [cm/day] over the surface -> volumetric uptake per hill [L/day]:
        # vRoot * area[cm^2] / 1000[cm^3/L] = vRoot * (area_m2 * 1e4) / 1e3
        Qtp = np.maximum(np.interp(t, res.t, res.vroot), 0.0) * area_per_hill_m2 * 10.0
    else:
        Qtp = q_direct

    return PlantInputs(t=t, Cwo=Cw_h, Qtp=Qtp, M=M, leaf_loss=leaf_loss, M_basis=om.MASS_BASIS,
                       dvs=dvs), res


# ---------------------------------------------------------------------------
def _demo():
    if not hydrus_available():
        print("HYDRUS executable not built — see module docstring. Skipping demo.")
        return
    print(f"HYDRUS-1D paddy soil model  (exe: {HYDRUS_EXE})")
    print(f"{'cong':8}{'nC':>3}{'Kd':>10}{'Cw0':>8}{'Cw_min':>8}{'Cw_end':>8}"
          f"{'buffered?':>11}")
    for name, n, grp in [("PFBA", 4, "PFCA"), ("PFHxA", 6, "PFCA"),
                         ("PFOA", 8, "PFCA"), ("PFNA", 9, "PFCA"),
                         ("PFOS", 8, "PFSA"), ("PFDoDA", 12, "PFCA")]:
        import forcing_rice as fr
        Kd = paddy_kd(n, grp)
        r = run_paddy_hydrus(Kd, tpot_fn=lambda d: fr.transpiration_mm_d(d, 120.0),
                             evap_fn=lambda d: fr.evaporation_mm_d(d, 120.0))
        cw = r.Cw / r.Cw[0]
        print(f"{name:8}{n:>3}{Kd:>10.3f}{cw[0]:>8.2f}{cw.min():>8.2f}{cw[-1]:>8.2f}"
              f"{('flat' if cw.min() > 0.9 else 'leached'):>11}")


if __name__ == "__main__":
    _demo()
