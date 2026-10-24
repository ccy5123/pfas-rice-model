"""Two-pool root model (the long-chain breakthrough), promoted to src/.

Locks in: (a) the module runs and returns finite positive BAFs; (b) the saturated
structural-adequacy fit (close_longchain) reproduces the PFDoDA (C12) long chain that
the single-pool core could NOT (refit_oryza hit ceilings ~4-6x under) -- root and grain
within a factor; (c) the canonical core path is unchanged. Scoped to PFDoDA to bound
runtime (the closure re-matches the carrier per f_xy step).
"""
import os
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))


def test_simulate_runs():
    import pfas_rice_two_pool as tp
    r = tp.simulate("PFDoDA", f_xy=0.4, vmax_in=100.0)
    for k in ("root", "straw", "grain", "rm", "rb"):
        assert k in r and r[k] == r[k] and r[k] >= 0.0      # finite, non-negative
    assert r["root"] > 0.0
    # the bound store holds most of the long-chain root burden (the 2-pool point)
    assert r["rb"] > r["rm"]


def test_breakthrough_closes_pfdoda():
    """The breakthrough: with a LOW free f_xy + an ENHANCED carrier, the 2-pool reproduces
    PFDoDA root AND grain simultaneously -- which the single-pool core cannot.

    Pinned to the pre-fix DRY organ mass. On the corrected FRESH mass (src/organ_mass.py)
    the reading does not hold -- see test_model_api_hook: the free f_xy runs to its bound
    and the carrier enhancement shrinks to ~1.5x, so this is a dry-basis result."""
    import pfas_rice_two_pool as tp
    import csv
    obs = {}
    with open(os.path.join(ROOT, "data_obs", "obs_baf_Yamazaki.csv")) as f:
        for row in csv.DictReader(f):
            obs.setdefault(row["compound"], {})[row["tissue"]] = float(row["baf"])
    d = tp.close_longchain("PFDoDA", obs["PFDoDA"], mass_basis="dry")
    sim, o = d["sim"], obs["PFDoDA"]
    assert sim["root"] == pytest.approx(o["root"], rel=0.10)     # root closes (~70 vs 69)
    assert sim["grain"] == pytest.approx(o["grain"], rel=0.20)   # grain closes (~46 vs 45.5)
    assert d["f_xy"] < 0.6                                       # LOW f_xy (strong root retention)
    assert d["carrier_x"] > 1.5                                  # ENHANCED carrier


def test_model_api_hook():
    import model_api as api
    d = api.close_longchain_2pool("PFDoDA")
    assert d["sim"]["root"] > 0 and 0 < d["f_xy"] <= 1.0
    # FRESH organ mass (the default), re-examined: the free f_xy sits at its bound (1)
    # and the carrier needs only ~1.5x (dry: f_xy 0.43, 3.2x), so "low f_xy + enhanced
    # carrier" was a dry-mass reading. Root and grain still close; the straw does NOT
    # (21 vs 49.7 even unrestricted) -- the C12 shoot-loading ceiling remains.
    o = api.observed_baf("PFDoDA")
    assert d["f_xy"] == pytest.approx(1.0, abs=1e-3)
    assert d["carrier_x"] < 2.0
    assert d["sim"]["root"] == pytest.approx(o["root"], rel=0.10)
    assert d["sim"]["grain"] == pytest.approx(o["grain"], rel=0.20)
    assert d["sim"]["straw"] < 0.5 * o["straw"]
    # and the canonical core path is still intact (default 4pool_surf simulate)
    base = api.simulate("PFOA")
    assert base["baf_final"]["root"] > 0
