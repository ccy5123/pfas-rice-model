"""
Reproduction guards for the STRUCTURAL MERGE (two-pool sequestration root +
redistributed N-stem+leaf shoot) -- validation/twopool_nstem_merge.py.

Pins the three published numbers so the merged model cannot drift silently:
  * Yamazaki in-sample log10 RMSE 0.301 (root 0.153) on the measured forcings,
    i.e. the merge costs ~0.02 against the pass-through-stem two-pool's own 0.278;
  * the non-K_PL PFOS/PFUnDA k_seq separation survives the merge (~4x at
    identical K_PL = 31623);
  * Tang 2026 per-organ OOS 0.801 -- vs 1.398 for the same root mechanism behind
    the pass-through stem. That drop is the whole point of the merge: it confirms
    the Result-7 diagnosis (the two-pool's Tang failure was a SHOOT artifact, not
    a root-mechanism failure) and is driven by the stalk (1.89 -> 0.61).

These use the CACHED merged fit (validation/twopool_nstem_fitted_params.json);
re-fitting is the ~40 min job the validation script runs.
"""
import os
import sys

import numpy as np
import pytest

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "src"))
sys.path.insert(0, os.path.join(_ROOT, "validation"))

pytest.importorskip("scipy")
import twopool_nstem_merge as MG  # noqa: E402


@pytest.fixture(scope="module")
def fit():
    MG.DRIVERS = MG.install_forcings("measured")
    import json
    d = json.load(open(MG.FIT_CACHE))
    return d["global"], np.array(d["ushape_q"])


def test_merged_in_sample_yamazaki_rmse(fit):
    """The merged model keeps the two-pool's in-sample quality (0.301 vs 0.278),
    and the root -- the compartment the two-pool mechanism exists for -- is
    essentially untouched by the shoot swap (0.153 vs 0.154)."""
    p, q = fit
    overall, by, _ = MG.rmse_report(
        p, lambda c: MG.TP.kseq_ushape(c["n_C"], c["group"], q))
    assert overall == pytest.approx(0.301, abs=0.03)
    assert by["root"] == pytest.approx(0.153, abs=0.02)
    assert by["straw"] == pytest.approx(0.319, abs=0.05)
    assert by["grain"] == pytest.approx(0.382, abs=0.05)


def test_non_kpl_separation_survives_the_merge(fit):
    """PFOS (C8 PFSA) and PFUnDA (C11 PFCA) share K_PL = 31623 exactly, so no
    K_PL-gated term can separate their root BAFs (5.93 vs 19.53). The chain +
    head-group k_seq still does, after the shoot was replaced."""
    _, q = fit
    ks_pfos = MG.TP.kseq_ushape(8, "PFSA", q)
    ks_pfunda = MG.TP.kseq_ushape(11, "PFCA", q)
    assert ks_pfunda / ks_pfos == pytest.approx(4.0, abs=0.6)


def test_merged_tang_per_organ_oos(fit):
    """The fair per-organ OOS the merge unlocks: transferring the Yamazaki fit to
    Tang (no Tang re-fit) drops the error from 1.398 (same root, pass-through
    stem) to ~0.80, and the stalk -- the organ the pass-through stem collapsed --
    carries the improvement. Still short of single-pool lipid loading (0.516),
    with the residual in the endosperm (a documented structural under-prediction).
    """
    p, q = fit
    pairs, by_organ = [], {}
    for nm in MG.api.TANG_CONGENERS:
        obs = MG.api.tang_observed_tf(nm, MG.DOSE)
        tf = MG.merged_tang_tf(nm, p, q)
        for _, organ, _ in MG.api._TANG_ORGANS:
            pairs.append((tf[organ], obs[organ]))
            by_organ.setdefault(organ, []).append((tf[organ], obs[organ]))
    assert MG._rmse(pairs) == pytest.approx(0.801, abs=0.06)
    assert MG._rmse(by_organ["stalk"]) == pytest.approx(0.61, abs=0.10)   # was 1.89
    assert MG._rmse(by_organ["leaf"]) == pytest.approx(0.28, abs=0.10)    # was 0.38
    # the residual is the grain, not the shoot
    assert MG._rmse(by_organ["endosperm"]) > MG._rmse(by_organ["stalk"])


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
