"""Tests for the real HYDRUS-1D soil coupling (src/soil_hydrus.py).

The HYDRUS-invoking tests skip automatically when the compiled executable or
phydrus is unavailable (e.g. CI without gfortran), so the suite stays green
everywhere; build the engine to exercise them (see src/soil_hydrus.py docstring).
"""
import numpy as np
import pytest

import soil_hydrus as sh

pytest.importorskip("phydrus")
hydrus = pytest.mark.skipif(not sh.hydrus_available(),
                            reason="HYDRUS-1D executable not built")


# ---- pure parameter logic (no engine needed) -------------------------------
def test_paddy_kd_increases_with_chain_length():
    kd = [sh.paddy_kd(n, "PFCA") for n in (4, 6, 8, 10, 12)]
    assert all(b > a for a, b in zip(kd, kd[1:]))           # strictly increasing
    assert kd[0] < 0.1 < kd[-1]                             # spans >3 orders


def test_paddy_kd_sulfonate_above_matched_carboxylate():
    # PFSA sorbs more than the CF2-matched PFCA (same perfluoro-C count)
    assert sh.paddy_kd(8, "PFSA") > sh.paddy_kd(9, "PFCA")  # PFOS(8 CF) vs PFNA(8 CF)


# ---- engine runs -----------------------------------------------------------
@hydrus
def test_run_paddy_short_chain_leaches():
    """A weakly-sorbed short chain must leach: pore water falls during flooding."""
    res = sh.run_paddy_hydrus(sh.paddy_kd(4, "PFCA"), season=60.0)
    assert np.all(np.isfinite(res.Cw)) and np.all(np.isfinite(res.vroot))
    assert res.Cw[0] == pytest.approx(1.0, abs=0.05)        # initial pore water
    assert res.Cw.min() < 0.5 * res.Cw[0]                   # leached away
    assert res.vroot.max() > 0                              # roots take up water


@hydrus
def test_run_paddy_long_chain_buffered():
    """A strongly-sorbed long chain stays buffered: pore water ~ flat."""
    res = sh.run_paddy_hydrus(sh.paddy_kd(12, "PFCA"), season=60.0)
    assert res.Cw.min() > 0.9 * res.Cw[0]                   # barely moves


@hydrus
def test_inputs_from_hydrus_normalised_and_shaped():
    inp, res = sh.inputs_from_hydrus(8, "PFCA", season=60.0, Cwo_ref=2.0, n_t=121)
    assert inp.Cwo.shape == inp.t.shape == (121,)
    assert np.mean(inp.Cwo) == pytest.approx(2.0, rel=1e-6)  # normalised to ref
    assert inp.M.shape == (121, 4)
    assert np.all(inp.M > 0) and np.all(np.isfinite(inp.Qtp))
