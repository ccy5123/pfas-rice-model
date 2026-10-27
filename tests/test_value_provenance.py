"""The value inventory (params/value_provenance.csv) guards the rule that no value the default
model runs on is chosen without a source (docs/PREREG_sourced_values.md)."""
import csv
import os

import pytest

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PATH = os.path.join(_ROOT, "params", "value_provenance.csv")
STATUSES = {"measured", "derived", "fitted", "inherited", "numerical", "structural",
            "convention", "chosen", "untraceable", "decision"}
ACTIONS = {"keep", "literature", "fit", "simplify", "range", "converge", "bootstrap"}
UNRESOLVED = {"chosen", "untraceable", "decision"}


def _rows():
    with open(PATH, encoding="utf-8") as f:
        return list(csv.DictReader(line for line in f if not line.startswith("#")))


def test_inventory_is_well_formed():
    rows = _rows()
    ids = [r["id"] for r in rows]
    assert len(ids) == len(set(ids)), "duplicate ids"
    for r in rows:
        assert r["status"] in STATUSES, r["id"]
        assert r["action"] in ACTIONS, r["id"]
        assert r["done"] in ("yes", "no"), r["id"]


def test_no_unsourced_value_in_the_default_model():
    open_rows = [r["id"] for r in _rows()
                 if r["scope"] == "default" and r["status"] in UNRESOLVED and r["done"] != "yes"]
    assert not open_rows, f"default-scope values without a source: {open_rows}"


# Default-scope rows whose action is not finished yet, named one by one: each is an open item
# of the next amendment (docs/PREREG_sourced_values.md). A row joins only when new evidence
# reopens it, with the reason here.
OPEN = {
    "obs_exposure_time",    # exposure held at the one irrigation-water grab over the season; its
                            # seasonal shape is built with the irrigation-source soil scenario (2b)
}


def test_every_open_default_row_is_named():
    pending = {r["id"] for r in _rows() if r["scope"] == "default" and r["done"] != "yes"}
    assert pending <= OPEN, f"unfinished default-scope rows not in OPEN: {sorted(pending - OPEN)}"
    assert OPEN <= {r["id"] for r in _rows()}, "OPEN names a row that does not exist"
