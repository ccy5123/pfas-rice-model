"""
measured_biomass.py — ingest a MEASURED per-organ rice biomass table -> M(t) driver
===================================================================================

Turns an observed per-organ biomass time series (e.g. Tang 2026 JHM SI, Nipponbare
paddy, full growth cycle) into the model's mass driver M(t) [kg/hill, (n,4)] and a
ready-to-run `drivers` dict for `model_api.simulate(drivers=...)`.  This is the
data-grounded alternative to the modelled biomass (`growth_rice`, `oryza_growth`):
it pins the absolute compartment-mass scale (and hence the f_xy absolute scale,
task #7) to a real experiment instead of an IR72 reconstruction.

Input table: sampling `day` (days after transplanting) + organ dry weights for
root, stem(culm), leaf, grain(panicle).  Several reporting units are accepted and
converted to **kg fresh-weight-free dry mass per hill** (the model's M unit):

    units            extra args needed
    -----            -----------------
    kg/hill          (none)
    g/hill           (none)
    mg/hill          (none)
    g/plant          plants_per_hill
    g/m2             hills_per_m2
    kg/ha | t/ha     hills_per_m2

The model's M is per-hill because Q_TP is per-hill (the dilution term is Q_TP/M and
the growth-dilution sink is (dM/dt)/M), so the per-hill conversion must be consistent
with the transpiration driver.  If the source reports water use / transpiration, pass
it as `Qtp`; otherwise the modelled `forcing_rice.Q_TP` is used (absolute-scale caveat:
mixing a measured biomass with a modelled transpiration is only consistent up to the
T/ET and density assumptions -- see CLAUDE.md task #7).
"""
from __future__ import annotations
import numpy as np

ORGANS = ("root", "stem", "leaf", "grain")


def to_kg_per_hill(values: np.ndarray, units: str, *, plants_per_hill: float | None = None,
                   hills_per_m2: float = 25.0) -> np.ndarray:
    """Convert an organ biomass column to kg(dry)/hill from the given reporting units."""
    v = np.asarray(values, float)
    u = units.strip().lower()
    if u == "kg/hill":
        return v
    if u == "g/hill":
        return v / 1.0e3
    if u == "mg/hill":
        return v / 1.0e6
    if u == "g/plant":
        if plants_per_hill is None:
            raise ValueError("units='g/plant' needs plants_per_hill")
        return v * plants_per_hill / 1.0e3
    if u == "g/m2":
        return v / 1.0e3 / hills_per_m2
    if u in ("kg/ha",):
        return v / 1.0e4 / hills_per_m2
    if u in ("t/ha", "ton/ha", "mg/ha"):           # 1 t/ha = 1000 kg/ha
        return v * 1.0e3 / 1.0e4 / hills_per_m2
    raise ValueError(f"unknown units {units!r}; supported: kg/hill,g/hill,mg/hill,"
                     f"g/plant,g/m2,kg/ha,t/ha")


def load_biomass_table(path_or_buffer):
    """Read a CSV with columns day,root,stem,leaf,grain (root optional).
    Returns (day, {organ: array}); missing organ columns -> None."""
    data = np.genfromtxt(path_or_buffer, delimiter=",", names=True)
    cols = data.dtype.names or ()
    daycol = "day" if "day" in cols else ("t" if "t" in cols else None)
    if daycol is None:
        raise ValueError(f"biomass table needs a 'day' (or 't') column; found {cols}")
    day = np.atleast_1d(data[daycol]).astype(float)
    organs = {o: (np.atleast_1d(data[o]).astype(float) if o in cols else None) for o in ORGANS}
    if all(organs[o] is None for o in ORGANS):
        raise ValueError(f"no organ columns found; need any of {ORGANS}; found {cols}")
    return day, organs


def biomass_matrix(t, day, organs_kg_hill, *, root_shoot_ratio: float | None = None):
    """Interpolate measured per-organ biomass [kg/hill] onto grid `t` -> M (n,4)
    in column order [root, stem, leaf, grain].

    A measured organ series is interpolated (and held flat beyond its last sample).
    If `root` is absent, it is reconstructed from `root_shoot_ratio` * (stem+leaf+grain)
    when given, else floored to a small value (root then plays no transport role)."""
    t = np.asarray(t, float)
    cols = {}
    for o in ORGANS:
        v = organs_kg_hill.get(o)
        if v is not None:
            cols[o] = np.interp(t, day, v)            # flat extrapolation at ends
    if "root" not in cols:
        if root_shoot_ratio is not None:
            shoot = sum(cols[o] for o in ("stem", "leaf", "grain") if o in cols)
            cols["root"] = root_shoot_ratio * shoot
        else:
            cols["root"] = np.full_like(t, 1e-6)
    M = np.column_stack([np.maximum(cols.get(o, np.full_like(t, 1e-9)), 1e-9) for o in ORGANS])
    return M


def biomass_drivers(path_or_table, units, *, season=None, n_t=241, Cwo=1.0,
                    plants_per_hill=None, hills_per_m2=25.0, root_shoot_ratio=None,
                    Qtp=None):
    """One-call: measured biomass table -> `drivers` dict for simulate(drivers=...).

    `path_or_table` is a CSV path/buffer, or a (day, organs) tuple from
    `load_biomass_table`.  `units` is the reporting unit of the table (see
    `to_kg_per_hill`).  `Qtp` (L/d/hill on the output grid) overrides the modelled
    transpiration; otherwise `forcing_rice.Q_TP` is used.  Returns the same dict
    shape as `model_api.drivers_from_arrays`.
    """
    if isinstance(path_or_table, tuple):
        day, organs_raw = path_or_table
    else:
        day, organs_raw = load_biomass_table(path_or_table)
    organs_kg = {o: (to_kg_per_hill(v, units, plants_per_hill=plants_per_hill,
                                    hills_per_m2=hills_per_m2) if v is not None else None)
                 for o, v in organs_raw.items()}
    season = float(np.max(day)) if season is None else float(season)
    t = np.linspace(0.0, season, n_t)
    M = biomass_matrix(t, day, organs_kg, root_shoot_ratio=root_shoot_ratio)
    Cwo_arr = np.full_like(t, float(Cwo))
    if Qtp is None:
        import forcing_rice as fr
        Qtp = fr.Q_TP(t, season)
    else:
        Qtp = np.asarray(Qtp, float)
    return dict(t=t, Cwo=Cwo_arr, Qtp=Qtp, M=M)


if __name__ == "__main__":
    import io
    # tiny synthetic demo (replace with the Tang 2026 SI table)
    demo = io.StringIO(
        "day,root,stem,leaf,grain\n"
        "0,0.05,0.02,0.03,0\n"
        "40,0.6,1.2,1.4,0\n"
        "80,1.0,5.5,4.0,3.0\n"
        "150,1.1,7.0,4.5,12.0\n")
    drv = biomass_drivers(demo, units="g/plant", plants_per_hill=3, season=150.0)
    Mf = drv["M"][-1]
    print("measured-biomass drivers built:")
    print(f"  grid: {drv['t'][0]:.0f}..{drv['t'][-1]:.0f} d, n={len(drv['t'])}")
    print(f"  final M [kg/hill]: root={Mf[0]:.4f} stem={Mf[1]:.4f} "
          f"leaf={Mf[2]:.4f} grain={Mf[3]:.4f}")
    print("  -> pass to model_api.simulate(drivers=drv)")
