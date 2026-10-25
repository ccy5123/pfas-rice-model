"""Li 2025 stage by stage (validation/li2025_stage_resolved.py, pre-registered in 4df3047).

Pinned: (1) the SI Table S6 transcription reproduces the repo's 5-stage means and the
authors' own TF; (2) the TF estimator alone moves the observed TF up to 5.8x; (3) the
stage-to-stage scatter is ~4.5x what the reported SDs imply, so per-stage shape is
noise-limited; (4) the pre-registered verdict -- the redistributed shoot with lipid
loading discriminates (the expectation was that Li would stay inconclusive), and the
post-hoc checks that it survives the other estimators and needs BOTH the shoot
structure and the lipid term.
"""
import contextlib
import csv
import io
import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "validation"))

import li2025_stage_resolved as L  # noqa: E402


@pytest.fixture(scope="module")
def result():
    with contextlib.redirect_stdout(io.StringIO()):
        return L.main(sensitivity=False)


def test_transcription_reproduces_the_repo_means_and_the_authors_tf():
    data = L.load_data()
    with open(os.path.join(ROOT, "data_obs", "Li2025_BAF_TF.csv")) as f:
        old = {r["congener"]: r for r in csv.DictReader(f)}
    for c in L.CONGENERS:
        for k in ("root", "stem", "leaf", "grain"):
            m = data[(c, k)][0]
            ref = old[c][f"C_{k}_dw"]
            if ref:
                assert np.nanmean(m) == pytest.approx(float(ref), abs=1.5e-3), (c, k)
            else:
                assert np.all(np.isnan(m)), (c, k)
    tf, _ = L.observed_tf(data, "PFOA", "stem")
    assert np.mean(tf) == pytest.approx(2.35, abs=0.005)          # the authors' published TF


def test_the_estimator_alone_moves_the_observed_tf(result):
    assert result["spread_median"] == pytest.approx(3.0, abs=0.2)
    assert result["spread_max"] == pytest.approx(5.8, abs=0.2)


def test_stage_scatter_exceeds_the_reported_sds(result):
    assert result["noise_limited"]
    assert result["noise_ratio"] == pytest.approx(4.5, abs=0.3)
    assert all(s["p"] > L.ALPHA for s in result["scores"].values())   # no shape signal


def test_the_pre_registered_verdict(result):
    sc = result["scores"]
    assert result["winner"] == "nstem_lipid"
    assert sc["nstem_lipid"]["rmse"]["shoot"] == pytest.approx(0.36, abs=0.02)
    assert sc["nstem_lipid"]["rmse"]["grain"] == pytest.approx(0.51, abs=0.02)
    assert sc["lipid"]["rmse"]["shoot"] == pytest.approx(0.90, abs=0.02)
    assert sc["lipid"]["rmse"]["grain"] == pytest.approx(1.24, abs=0.03)


def test_posthoc_estimators_and_decomposition(result):
    post = result["posthoc"]
    for est in ("mean_of_ratios", "ratio_of_means"):
        row = post[est]
        assert L.verdict({v: {"rmse": row[v]} for v in L.VARIANTS}) == "nstem_lipid"
    dec = post["decomposition"]
    for g in ("shoot", "grain"):
        assert dec["nstem_lipid"][g] < dec["nstem_mono"][g] - 0.3      # the lipid term matters
        assert dec["nstem_lipid"][g] < dec["lipid"][g] - 0.3           # so does the shoot
        assert dec["nstem_lipid"][g] > post["noise_floor"][g]          # but not at the floor


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
