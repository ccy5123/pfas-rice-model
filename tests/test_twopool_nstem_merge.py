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

These use the CACHED merged fits -- one per organ-mass basis (src/organ_mass.py):
validation/twopool_nstem_fitted_params_dry.json is the record above (the biomass
driver's DRY matter used as M); validation/twopool_nstem_fitted_params.json is the
re-fit on FRESH mass, where the staged sequestration fit degrades in-sample
(0.301 -> 0.515: the root-match floors k_seq at 0 for most congeners, so the
U-shape is pulled by a few points) while the Tang OOS improves (0.801 -> 0.602),
still behind single-pool lipid loading (0.520). Re-fitting is the ~10 min job the
validation script runs.
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


def _use_basis(basis):
    """Point the module at one organ-mass basis (its cache, forcings and Tang
    transfer all follow MG.MASS_BASIS); returns the cached (p, q)."""
    import json
    MG.MASS_BASIS = basis
    MG.DRIVERS = MG.install_forcings("measured")   # explicit; import must not do this
    d = json.load(open(MG.fit_cache(basis)))
    return d["global"], np.array(d["ushape_q"])


@pytest.fixture
def fit():
    """The DRY-mass record (in-sample 0.301, Tang OOS 0.801)."""
    saved = MG.MASS_BASIS, MG.DRIVERS
    yield _use_basis("dry")
    MG.MASS_BASIS, MG.DRIVERS = saved


@pytest.fixture
def fit_fresh():
    saved = MG.MASS_BASIS, MG.DRIVERS
    yield _use_basis("fresh")
    MG.MASS_BASIS, MG.DRIVERS = saved


def test_import_does_not_mutate_shared_forcings():
    """Regression: this module used to swap `twopool_root_exploration`'s forcing
    globals at IMPORT time. pytest imports every test module during collection, so
    that silently corrupted an unrelated test (test_model_api's two-pool drift
    guard reads those globals) -- and only when the full suite ran, never in
    isolation. Importing must leave them exactly as the two-pool arc set them."""
    import twopool_root_exploration as TP
    assert TP.T[-1] == pytest.approx(120.0)
    assert TP.QTP.max() == pytest.approx(0.40, abs=0.01)     # demo, not measured 0.098
    assert MG.install_forcings("measured")["Qtp"].max() < 0.15   # ... and it stays demo
    assert TP.QTP.max() == pytest.approx(0.40, abs=0.01)


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


def test_caches_are_per_basis_and_default_fresh():
    assert MG.fit_cache("dry").endswith("twopool_nstem_fitted_params_dry.json")
    assert MG.fit_cache("fresh").endswith("twopool_nstem_fitted_params.json")
    import json
    assert json.load(open(MG.fit_cache("fresh")))["mass_basis"] == "fresh"


def test_fresh_refit_in_sample_and_tang(fit_fresh):
    """On FRESH organ mass the merged fit degrades in-sample (the staged U-shape is
    fitted to root-matched k_seq that floor at 0 for most congeners) but transfers to
    Tang better than the dry record -- and still not better than lipid loading."""
    p, q = fit_fresh
    overall, by, _ = MG.rmse_report(
        p, lambda c: MG.TP.kseq_ushape(c["n_C"], c["group"], q))
    assert overall == pytest.approx(0.515, abs=0.04)
    assert overall > 0.301 + 0.1                     # worse than the dry record
    pairs = []
    for nm in MG.api.TANG_CONGENERS:
        obs = MG.api.tang_observed_tf(nm, MG.DOSE)
        tf = MG.merged_tang_tf(nm, p, q)
        pairs += [(tf[o], obs[o]) for _, o, _ in MG.api._TANG_ORGANS]
    assert MG._rmse(pairs) == pytest.approx(0.602, abs=0.06)
    # the non-K_PL PFOS/PFUnDA separation does not carry over: PFOS's sink is ~0
    assert MG.TP.kseq_ushape(8, "PFSA", q) < 1e-2


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
