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


def paddy_kd(n_C: int, group: str = "PFCA", f_oc: float = 0.02) -> float:
    """Linear soil distribution coefficient Kd [L/kg] for a congener.

    Kd = Koc(chain length, head group) * f_oc, with Koc from the C3 QSPR
    (Higgins & Luthy slope, Milinovic PFOA anchor) in ``literature_params``.
    PFCA C_n has n-1 perfluorinated carbons; PFSA C_n has n.
    """
    head = "sulfonate" if group.upper() == "PFSA" else "carboxylate"
    n_pfc = n_C - 1 if group.upper() == "PFCA" else n_C
    return lp.koc_to_KF(lp.koc(n_pfc, head), f_oc)


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
@dataclass
class PaddyResult:
    t: np.ndarray        # output times [day]
    Cw: np.ndarray       # pore-water conc at root zone [mg/cm^3, raw HYDRUS]
    theta: np.ndarray    # water content at root zone [-]
    vroot: np.ndarray    # actual transpiration / root water uptake [cm/day]
    Kd: float            # linear sorption coeff used [L/kg]


def run_paddy_hydrus(Kd: float, season: float = 120.0, *, depth: float = 100.0,
                     dx: float = 2.0, bulk_density: float = 1.3, disper: float = 5.0,
                     flood_until: float = 90.0, percolation: float = 0.30,
                     obs_depth: float = 10.0, tpot_fn=None, exe: str = HYDRUS_EXE,
                     workspace: str | None = None, keep: bool = False) -> PaddyResult:
    """Build and run a one-season paddy HYDRUS-1D model for a given Kd.

    Scenario: a contaminated paddy soil (uniform initial pore-water conc = 1)
    under continuous flooding to ``flood_until`` days (irrigation = ET + a
    ``percolation`` excess, clean water -> leaching) then drainage to harvest.
    Linear Kd sorption, no decay (PFAS recalcitrant).  Returns daily series.

    tpot_fn : optional callable mapping a day array -> potential transpiration
        [cm/day] for the atmospheric ``rRoot``.  Pass the MEASURED transpiration
        (e.g. ``forcing_rice.transpiration_mm_d/10``) so HYDRUS's actual root
        water uptake ``vRoot`` carries the measured crop-physiology shape (plus
        any soil-water-stress feedback).  Default is an illustrative gaussian.
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

    ws = workspace or tempfile.mkdtemp(prefix="hydrus_paddy_")
    try:
        ml = ps.Model(exe_name=exe, ws_name=ws, name="paddy", mass_units="mg",
                      length_unit="cm", time_unit="days", print_screen=False)
        ml.add_time_info(tinit=0, tmax=season, dt=1e-4, dtmin=1e-7, dtmax=0.5,
                         print_times=True)
        ml.add_waterflow(model=0, top_bc=3, bot_bc=4, maxit=20, tolth=1e-3, tolh=1)
        # fully-implicit + upstream weighting -> robust solute solve at the
        # clean-water (c -> 0) top boundary
        ml.add_solute_transport(model=0, top_bc=-1, bot_bc=0, epsi=1.0, lupw=True,
                                pecr=2)

        m = ml.get_empty_material_df(n=1)
        m.loc[1, ("water", "thr")] = 0.08
        m.loc[1, ("water", "ths")] = 0.46
        m.loc[1, ("water", "Alfa")] = 0.016
        m.loc[1, ("water", "n")] = 1.37
        m.loc[1, ("water", "Ks")] = 15.0
        m.loc[1, ("water", "l")] = 0.5
        m.loc[1, ("solute", "bulk.d")] = bulk_density
        m.loc[1, ("solute", "DisperL")] = disper
        m.loc[1, ("solute", "frac")] = 1.0
        m.loc[1, ("solute", "mobile_wc")] = 0.0
        ml.add_material(m)

        prof = ps.create_profile(top=0, bot=-abs(depth), dx=dx, h=-25.0, mat=1,
                                 lay=1, conc=1.0)
        z = prof["x"].values
        beta = np.clip(1.0 + z / 30.0, 0, None)   # roots: surface -> 0 at -30 cm
        beta[z < -30] = 0.0
        prof["Beta"] = beta
        ml.add_profile(prof)
        ml.add_obs_nodes([-abs(obs_depth)])
        ml.add_root_uptake(model=0, poptm=[-25.0])

        days = np.arange(1, int(season) + 1)
        if tpot_fn is not None:
            tpot = np.maximum(np.asarray(tpot_fn(days), float), 0.0)   # [cm/day]
        else:
            tpot = 0.05 + 0.55 * np.exp(-((days - 0.62 * season) ** 2) / (2 * (season / 5.5) ** 2))
        esoil = 0.10 * np.ones_like(days, float)
        prec = np.where(days <= flood_until, esoil + tpot + percolation, 0.0)
        import pandas as pd
        atm = pd.DataFrame({"tAtm": days.astype(float), "Prec": prec, "rSoil": esoil,
                            "rRoot": tpot, "hCritA": 1e5, "rB": 0.0, "hB": 0.0,
                            "ht": 0.0, "tTop": 0.0, "tBot": 0.0, "Ampl": 0.0})
        ml.add_atmospheric_bc(atm, hcrits=1.0, tatm=0.0, prec=0.0, rsoil=0.0,
                              rroot=0.0, hcrita=1e5, rb=0.0, hb=0.0, ht=0.0,
                              ttop=0.0, tbot=0.0, ampl=0.0)

        sdf = ml.get_empty_solute_df()
        sdf.loc[1, "ks"] = float(Kd)     # linear isotherm: s = Kd * c
        sdf.loc[1, "beta"] = 1.0
        sdf.loc[1, "nu"] = 0.0
        ml.add_solute(sdf, difw=0.43, difg=0.0, top_conc=0.0, bot_conc=0.0)

        ml.write_input()
        ml.simulate()

        err = os.path.join(ws, "Error.msg")
        if os.path.isfile(err) and os.path.getsize(err) > 0:
            raise RuntimeError("HYDRUS did not converge: "
                               + open(err).read().strip())

        t, Cw, theta = _parse_obs_node(os.path.join(ws, "OBS_NODE.OUT"))
        tv, vroot = _parse_tlevel(os.path.join(ws, "T_LEVEL.OUT"))
        vroot = np.interp(t, tv, vroot)
        return PaddyResult(t=t, Cw=Cw, theta=theta, vroot=vroot, Kd=float(Kd))
    finally:
        if not keep:
            shutil.rmtree(ws, ignore_errors=True)


# ---------------------------------------------------------------------------
# coupling: HYDRUS soil -> PlantInputs
# ---------------------------------------------------------------------------
def inputs_from_hydrus(congener_n_C: int, group: str = "PFCA", *, season: float = 120.0,
                       Cwo_ref: float = 1.0, f_oc: float = 0.02,
                       qtp_from_hydrus: bool = True, area_per_hill_m2: float | None = None,
                       n_t: int = 241, **run_kw):
    """Build a :class:`PlantInputs` whose ``Cwo`` AND ``Q_TP`` come from a real
    HYDRUS-1D paddy run for the congener, with growth ``M`` from ORYZA.

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
    Returns (PlantInputs, PaddyResult).
    """
    from pfas_rice_plant_module_4pool_surf import PlantInputs
    import growth_rice as gr
    import forcing_rice as fr

    if area_per_hill_m2 is None:
        area_per_hill_m2 = fr.AREA_PER_HILL_M2

    Kd = paddy_kd(congener_n_C, group, f_oc)
    # drive the soil run with the measured transpiration shape (mm/day -> cm/day)
    tpot_fn = (lambda d: fr.transpiration_mm_d(d, season) / 10.0) if qtp_from_hydrus else None
    res = run_paddy_hydrus(Kd, season=season, tpot_fn=tpot_fn, **run_kw)

    t = np.linspace(0.0, season, n_t)
    Cw_h = np.interp(t, res.t, res.Cw)
    Cw_h = Cw_h * (Cwo_ref / np.mean(Cw_h))          # normalise to mean = Cwo_ref

    b = gr.organ_biomass(t, season)
    M = np.maximum(np.column_stack([b["root"], b["stem"], b["leaf"], b["grain"]]), 1e-4)

    if qtp_from_hydrus:
        # vRoot [cm/day] over the surface -> volumetric uptake per hill [L/day]:
        # vRoot * area[cm^2] / 1000[cm^3/L] = vRoot * (area_m2 * 1e4) / 1e3
        Qtp = np.maximum(np.interp(t, res.t, res.vroot), 0.0) * area_per_hill_m2 * 10.0
    else:
        Qtp = fr.Q_TP(t, season)

    return PlantInputs(t=t, Cwo=Cw_h, Qtp=Qtp, M=M), res


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
        Kd = paddy_kd(n, grp)
        r = run_paddy_hydrus(Kd)
        cw = r.Cw / r.Cw[0]
        print(f"{name:8}{n:>3}{Kd:>10.3f}{cw[0]:>8.2f}{cw.min():>8.2f}{cw[-1]:>8.2f}"
              f"{('flat' if cw.min() > 0.9 else 'leached'):>11}")


if __name__ == "__main__":
    _demo()
