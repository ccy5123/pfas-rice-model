"""
Rice biomass forcing M_s(t)  (task 2, part 2 — the compartment-mass driver)
===========================================================================

Generates per-organ dry-biomass trajectories (root, stem, leaf, panicle/grain)
over the season from the ORYZA2000 / ORYZA(v3) DVS-based assimilate-partitioning
scheme, then maps them onto the N-segment uptake model's compartments.  This is
the second half of the task-2 forcing (the first is `forcing_rice.Q_TP`); together
they fix the gradient-flip crossover `B* ~ Q_s/(M_s*mu_s)` and the absolute f_xy.

Source: Bouman & van Laar 2006 (Agric. Syst. 87:249-273) and Li et al. 2017
(ORYZA v3, Agric. For. Meteorol. 237-238:246-256).  IR72 partitioning tables
(read from Li 2017 SI Table S1; cross-checked FLV+FST+FSO = 1 at each DVS):

    DVS:          0      0.5    0.75   1.0    1.2    2.5
    FLV (leaf):   0.52   0.45   0.58   0.03   0.02   0.01
    FST (stem):   0.48   0.55   0.34   0.43   0.07   0.00
    FSO (panicle):0      0      0      0.54   0.91   0.99   (storage organ; 0 pre-flowering)
    FSH (shoot frac of total): 0.45 (DVS0) -> 0.85 (DVS0.5) -> ~1.0 (DVS>=1); FRT = 1-FSH

Convention: daily assimilate is split root vs shoot by FSH, then the shoot share
is split leaf/stem/panicle by FLV/FST/FSO.  Biomass = time-integral of the splits.

The N-segment model folds leaves into each stem segment (Yamazaki lumps
"stem incl. leaves"), so a stem(+leaf) segment mass = (M_stem+M_leaf)/N.
"""
from __future__ import annotations
import numpy as np

# ORYZA IR72 partitioning tables (Li 2017 SI S1)
_DVS = np.array([0.0, 0.5, 0.75, 1.0, 1.2, 2.5])
_FLV = np.array([0.52, 0.45, 0.58, 0.03, 0.02, 0.01])
_FST = np.array([0.48, 0.55, 0.34, 0.43, 0.07, 0.00])
_FSO = np.array([0.00, 0.00, 0.00, 0.54, 0.91, 0.99])
_FSH_DVS = np.array([0.0, 0.5, 1.0, 2.5])
_FSH = np.array([0.45, 0.85, 1.00, 1.00])

# IR72 final aboveground biomass ~17.4 t/ha = 1740 g/m^2 (Bouman 2006); add roots.
WSHOOT_MAX_G_M2 = 1740.0
AREA_PER_HILL_M2 = 0.04                  # 25 hills/m^2 (matches forcing_rice)
G_M2_TO_KG_HILL = AREA_PER_HILL_M2 / 1000.0


def dvs_of_t(t: np.ndarray, season: float, t_flower_frac: float = 0.54) -> np.ndarray:
    """Development stage DVS(t): 0 at transplant, 1 at flowering, 2 at maturity."""
    f = np.asarray(t, float) / season
    tf = t_flower_frac
    return np.where(f < tf, f / tf, 1.0 + (f - tf) / (1.0 - tf))


def _logistic(t, M0, Mmax, k, t0):
    return Mmax / (1.0 + (Mmax / M0 - 1.0) * np.exp(-k * (t - t0)))


def organ_biomass(t: np.ndarray, season: float = 120.0,
                  wshoot_max_g_m2: float = WSHOOT_MAX_G_M2,
                  root_shoot: float | None = None):
    """Per-organ dry biomass [kg/hill] over t: returns dict root/stem/leaf/grain.

    Total assimilate follows a logistic; daily increments are partitioned by the
    ORYZA IR72 DVS tables and integrated.  Scaled so final shoot ~ wshoot_max.

    root_shoot : if given, rescale the WHOLE root trajectory so the final
        root:shoot mass ratio equals this value (the root time-shape is preserved).
        The DVS scheme drives root partitioning to ~0 after flowering, leaving an
        unrealistically low final root:shoot (~0.035); literature lowland-rice
        maturity root:shoot is ~0.08-0.13 (root ~7-12% of total plant; declines
        from ~0.2 at seedling, lower for high-yield cultivars due to grain-fill
        dilution -- Frontiers Plant Sci. 2021, 10.3389/fpls.2021.713814; Yoshida
        1981 IRRI). Default None preserves the original behaviour (reproducibility).
    """
    t = np.asarray(t, float)
    dvs = dvs_of_t(t, season)
    # total (root+shoot) biomass: standard logistic with midpoint at flowering
    # (~50% of biomass by DVS=1, ~99% by maturity) -> realistic harvest index.
    tmid = 0.54 * season
    r = 12.0 / season                    # ~0.10/day for a 120-day season
    Wtot = (wshoot_max_g_m2 / 0.9) / (1.0 + np.exp(-r * (t - tmid)))
    dW = np.clip(np.gradient(Wtot, t), 0.0, None)
    fsh = np.interp(dvs, _FSH_DVS, _FSH)
    flv = np.interp(dvs, _DVS, _FLV)
    fst = np.interp(dvs, _DVS, _FST)
    fso = np.interp(dvs, _DVS, _FSO)
    # normalize the shoot split to sum 1 (guards the small DVS=0.75 rounding)
    s = flv + fst + fso; s[s == 0] = 1.0
    flv, fst, fso = flv / s, fst / s, fso / s
    dshoot = fsh * dW
    dt = t[1] - t[0]
    out = {"root": np.cumsum((1 - fsh) * dW) * dt,
           "leaf": np.cumsum(flv * dshoot) * dt,
           "stem": np.cumsum(fst * dshoot) * dt,
           "grain": np.cumsum(fso * dshoot) * dt}
    # renormalize to the target final shoot biomass, then -> kg/hill
    shoot_final = out["leaf"][-1] + out["stem"][-1] + out["grain"][-1]
    scale = (wshoot_max_g_m2 / shoot_final) * G_M2_TO_KG_HILL
    res = {k: v * scale for k, v in out.items()}
    if root_shoot is not None:
        sh_f = res["stem"][-1] + res["leaf"][-1] + res["grain"][-1]
        rf = res["root"][-1]
        if rf > 0:
            res["root"] = res["root"] * (float(root_shoot) * sh_f / rf)
    return res


def M_for_nstem(t: np.ndarray, N: int = 4, season: float = 120.0):
    """Mass matrix for NStemModel: columns [root, stem_1..stem_N, grain] [kg/hill].
    Stem segments fold in leaves (Yamazaki convention) and share (M_stem+M_leaf)/N.
    A floor avoids division-by-zero in the ODE before emergence."""
    b = organ_biomass(t, season)
    straw = b["stem"] + b["leaf"]
    seg = np.maximum(straw / N, 1e-9)
    cols = [np.maximum(b["root"], 1e-9)] + [seg] * N + [np.maximum(b["grain"], 1e-9)]
    return np.column_stack(cols)


if __name__ == "__main__":
    season = 120.0
    t = np.linspace(0, season, 481)
    b = organ_biomass(t, season)
    Mf = {k: v[-1] for k, v in b.items()}
    shoot = Mf["stem"] + Mf["leaf"] + Mf["grain"]
    print("ORYZA IR72 biomass (final, kg/hill):",
          {k: round(v, 4) for k, v in Mf.items()})
    print(f"  final shoot = {shoot:.4f} kg/hill ({shoot/G_M2_TO_KG_HILL:.0f} g/m^2; target ~1740)")
    print(f"  root:shoot = {Mf['root']/shoot:.2f}  (typical flooded rice ~0.1-0.2)")
    print(f"  grain:shoot (harvest index) = {Mf['grain']/shoot:.2f}  (modern rice ~0.45-0.5)")
    # grain onset
    g = b["grain"]; onset = t[np.argmax(g > 0.01 * g[-1])]
    print(f"  grain fill onset ~ day {onset:.0f} (DVS=1 at day {0.54*season:.0f})")
    for d in (20, 45, 70, 95, 120):
        i = np.argmin(abs(t - d))
        print(f"   day {d:3d} (DVS {dvs_of_t(np.array([d]),season)[0]:.2f}): "
              f"root={b['root'][i]:.4f} stem={b['stem'][i]:.4f} "
              f"leaf={b['leaf'][i]:.4f} grain={b['grain'][i]:.4f} kg/hill")
