"""validation/decision_stats.py: the margin-free decision rule."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                "validation"))
import decision_stats as D   # noqa: E402


def test_identical_errors_tie():
    e = np.array([0.1, -0.2, 0.3, 0.05, -0.1])
    r = D.decide(e, e)
    assert r["verdict"] == "tie" and r["lo"] == 0.0 and r["hi"] == 0.0


def test_uniformly_smaller_errors_are_better():
    a = np.array([0.8, -0.9, 1.0, 0.7, -0.85, 0.95])
    r = D.decide(a, a / 4)
    assert r["verdict"] == "better" and r["hi"] < 0 and r["delta"] < 0
    assert D.decide(a / 4, a)["verdict"] == "worse"


def test_one_lucky_item_is_a_tie():
    a = np.array([0.30, 0.30, 0.30, 0.30, 0.30])
    b = np.array([0.30, 0.30, 0.30, 0.30, 0.00])
    assert D.decide(a, b)["verdict"] == "tie"


def test_rejects_mismatched_input():
    with pytest.raises(ValueError):
        D.paired_interval([0.1, 0.2], [0.1])
