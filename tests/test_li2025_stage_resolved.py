"""Li 2025 stage by stage (validation/li2025_stage_resolved.py, pre-registered in 4df3047).

Pinned: (1) the SI Table S6 transcription reproduces the repo's 5-stage means and the
authors' own TF; (2) the TF estimator alone moves the observed TF up to 5.8x; (3) the
stage-to-stage scatter exceeds what the reported SDs imply (chi-square), so per-stage shape
is noise-limited; (4) the stage -> model-day mapping is the sourced one (PREREG 5.9); (5) the
level scores equal the transfer record of the same fitted variants
(validation/transfer_sourced.json).
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
    assert result["noise_ratio"] == pytest.approx(4.5, abs=0.3)      # data only


def test_the_stage_mapping_is_the_sourced_one():
    """PREREG 5.9: heading and ripening by the DVS definitions, tillering at the midpoint of
    [transplanting, panicle initiation], filling and milk-ripening at the centroid of their
    ordered range -- no chosen DVS."""
    lo, hi = L.tillering_range()
    til, head, fill, milk, ripe = L.stage_dvs()
    assert (head, ripe) == (1.0, 2.0) and hi == 0.65 and lo < hi
    assert til == pytest.approx(0.5 * (lo + hi))
    assert (fill, milk) == pytest.approx((4.0 / 3.0, 5.0 / 3.0))
    assert np.all(np.diff(L.stage_days()) > 0)


def test_level_scores_match_the_transfer_record(result):
    import json
    with open(os.path.join(ROOT, "validation", "transfer_sourced.json")) as f:
        rec = json.load(f)["errors"]
    same = {"nstem_lipid": "redistributed_lipid", "lipid": "4pool_lipid", "mono": "4pool_free",
            "fit": "redistributed_fit"}
    sc = result["scores"]
    for v, key in same.items():
        for g in ("shoot", "grain"):
            r = float(np.sqrt(np.mean(np.square(rec[key][f"li_{g}"]))))
            assert sc[v]["rmse"][g] == pytest.approx(r, abs=1e-6), (v, g)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
