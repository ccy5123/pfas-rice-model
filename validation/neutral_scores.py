#!/usr/bin/env python3
# =============================================================================
# validation/neutral_scores.py
# -----------------------------------------------------------------------------
# The NEUTRAL path's a-priori scores, recomputed and written to params/neutral_scores.json
# so that model_api.NEUTRAL_APRIORI_OOS_SCORES (the neutral inverse's model error) and the
# app read them from a record instead of numbers typed into the code.
#
# Nothing is fitted. Each score is the one its own script reports (log10 RMSE, the default
# lipid reading, season-ODE basis):
#   Liu 2023 rice root ............ neutral_dpu_validation.compare_to_obs (liu2023)
#   Ge 2017 per-organ TF .......... neutral_dpu_validation.compare_to_obs (ge2017)
#   Kodesova 2019 carbamazepine ... compare_to_obs (kodesova2019)  (kodesova2019_carbamazepine.py)
#   Li 2019 hydroponic ............ compare_to_obs (li2019_rcf)
#   Li 2019 soil, through the ODE . compare_to_obs (li2019_soil)
#   Briggs 1983 barley stem ....... briggs1983_shoot.section2: repo K_PW x TSCF, central stem
# plus the equilibrium-basis scores of the three root-partition tables (compare_to_obs
# mode="equilibrium"), which the docs quote next to the ODE ones.
# Composition (PREREG 5.10): each table runs on its own plant's sourced water and lipid; each
# score is reported at the midpoints with its [min, max] over the range corners, and Liu 2023
# (basis unstated) is also scored read as dry weight.
#
#   python validation/neutral_scores.py
# =============================================================================
from __future__ import annotations

import contextlib
import io
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, HERE)
import neutral_dpu_validation as V      # noqa: E402
import briggs1983_shoot as BS           # noqa: E402
import neutral_composition as NC        # noqa: E402

OUT = os.path.join(ROOT, "params", "neutral_scores.json")
TABLES = {
    "Liu 2023 rice root": "data_obs/neutral_obs_liu2023.csv",
    "Ge 2017 per-organ TF": "data_obs/neutral_obs_ge2017.csv",
    "Kodesova 2019 carbamazepine root": "data_obs/neutral_obs_kodesova2019.csv",
    "Li 2019 hydroponic": "data_obs/neutral_obs_li2019_rcf.csv",
    "Li 2019 soil, through the ODE": "data_obs/neutral_obs_li2019_soil.csv",
}
EQUILIBRIUM = ("Liu 2023 rice root", "Kodesova 2019 carbamazepine root",
               "Li 2019 soil, through the ODE")


def _briggs_stem(corner=NC.MID):
    rows = BS.load()
    return BS._score(rows, "stem_central", lambda x: BS.repo_scf(x, corner))


def main():
    drv = V.drivers()
    scores, equilibrium, corners = {}, {}, {}
    with contextlib.redirect_stdout(io.StringIO()):
        for label, rel in TABLES.items():
            path = os.path.join(ROOT, rel)
            scores[label] = float(V.compare_to_obs(path, drv, quiet=True))
            cs = [float(V.compare_to_obs(path, drv, quiet=True, corner=c)) for c in NC.CORNERS]
            corners[label] = [min(cs), max(cs)]
        for label in EQUILIBRIUM:
            equilibrium[label] = float(V.compare_to_obs(os.path.join(ROOT, TABLES[label]), drv,
                                                        quiet=True, mode="equilibrium"))
        liu_dw = float(V.compare_to_obs(os.path.join(ROOT, TABLES["Liu 2023 rice root"]), drv,
                                        quiet=True, read_as="dw"))
        n, bias, e = _briggs_stem()
        cs = [_briggs_stem(c)[2] for c in NC.CORNERS]
    scores["Briggs 1983 barley stem"] = float(e)
    corners["Briggs 1983 barley stem"] = [float(min(cs)), float(max(cs))]
    for k, v in scores.items():
        print(f"  {k:36s} {v:.3f}  [{corners[k][0]:.3f}, {corners[k][1]:.3f}]"
              + (f"   (equilibrium basis {equilibrium[k]:.3f})" if k in equilibrium else ""))
    print(f"  Liu 2023 read as dry weight          {liu_dw:.3f}")
    with open(OUT, "w") as f:
        json.dump(dict(source="validation/neutral_scores.py (a priori; nothing fitted)",
                       units="log10 RMSE", scores={k: round(v, 4) for k, v in scores.items()},
                       equilibrium_basis={k: round(v, 4) for k, v in equilibrium.items()},
                       range_corners={k: [round(a, 4), round(b, 4)] for k, (a, b) in corners.items()},
                       liu2023_read_as_dw=round(liu_dw, 4),
                       composition="PREREG 5.10: params/neutral_composition_sourced.csv",
                       briggs_stem_n=int(n), briggs_stem_bias=round(float(bias), 4)),
                  f, indent=1)
    print(f"[written] {os.path.relpath(OUT, ROOT)}")
    return scores


if __name__ == "__main__":
    main()
