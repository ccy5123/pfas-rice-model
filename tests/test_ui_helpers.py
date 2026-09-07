"""Sidebar widget-range guards — the repo's first tests that touch `ui/`.

WHY THIS MODULE EXISTS. Streamlit raises `StreamlitValueBelowMinError` /
`StreamlitValueAboveMaxError` when a widget's `value` falls outside its
`min_value`/`max_value`, and `app.py` builds the sidebar at the TOP of the script — so
ONE out-of-range default aborts the entire page before anything renders, and the user
cannot even edit the offending field to recover. A live report of exactly that: looking
up water (log Kow −1.38) put Karickhoff's Koc at 0.019 L/kg into a field whose minimum
was 0.1, and the whole app died with a redacted error.

The bug was not the one input. The soil-Koc field's range has to cover what the panel
ITSELF can produce from its own log Kow field, and it did not — at either end.
"""
import sys
import os

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

pytest.importorskip("streamlit", reason="the sidebar helpers need Streamlit")

import literature_params as LP           # noqa: E402
from ui.common import clamp_seed         # noqa: E402
import ui.sidebar as SB                  # noqa: E402


# The log Kow field's own bounds, and the soil-Koc field's, as the sidebar declares
# them. Kept here as the contract these tests enforce; if the widget changes, change
# these and the assertions below will tell you whether the ranges still line up.
LOGKOW_MIN, LOGKOW_MAX = -2.0, 8.0
KOC_MIN, KOC_MAX = 0.0, SB._KOC_MAX


def test_the_koc_field_covers_every_koc_its_own_log_kow_field_can_produce():
    """The regression. `koc_neutral` is applied to whatever the log Kow field holds, so
    the Koc field must accept the image of that whole interval — otherwise a legitimate
    compound (water at the bottom, anything very lipophilic at the top) crashes the
    page. The old [0.1, 1e6] failed BOTH ends."""
    lo = LP.koc_neutral(LOGKOW_MIN)
    hi = LP.koc_neutral(LOGKOW_MAX)
    assert lo == pytest.approx(0.00474, rel=1e-2)
    assert hi == pytest.approx(3.68e7, rel=1e-2)
    assert KOC_MIN <= lo and hi <= KOC_MAX, (
        f"soil Koc field [{KOC_MIN:g}, {KOC_MAX:g}] does not cover Karickhoff over "
        f"log Kow [{LOGKOW_MIN}, {LOGKOW_MAX}] = [{lo:g}, {hi:g}] — the crash is back")
    # the specific value that was reported: water
    water = LP.koc_neutral(-1.38)
    assert water == pytest.approx(0.0195, rel=1e-2)
    assert KOC_MIN <= water <= KOC_MAX
    assert not (0.1 <= water), "water's Koc really is below the old 0.1 minimum"


def test_the_koc_ceiling_matches_what_the_lookup_will_hand_over():
    """A Koc the lookup is willing to return must be displayable, or the fix just moves
    the crash from Karickhoff to CompTox."""
    cl = pytest.importorskip("chem_lookup")
    assert KOC_MAX >= cl.KOC_MAX_LKG


def test_clamp_seed_never_lets_a_looked_up_value_leave_the_range():
    """Out-of-range seeds are clamped and REPORTED, not silently accepted (which would
    change a scientific input without saying so) and not raised (which would abort the
    page the user needs in order to fix it)."""
    v, note = clamp_seed(0.0195, 0.1, 1e6)
    assert v == 0.1 and note and "below" in note and "0.0195" in note
    v, note = clamp_seed(3.7e7, 0.1, 1e6)
    assert v == 1e6 and note and "above" in note
    v, note = clamp_seed(2.45, -2.0, 8.0)
    assert v == 2.45 and note is None                 # in range -> untouched, no noise
    v, note = clamp_seed(float("nan"), 0.0, 1.0)      # NaN has no sensible clamp
    assert v == 0.0 and note and "NaN" in note
    # the bounds themselves are inclusive
    assert clamp_seed(-2.0, -2.0, 8.0) == (-2.0, None)
    assert clamp_seed(8.0, -2.0, 8.0) == (8.0, None)


def test_every_seeded_sidebar_field_is_clamped():
    """A guard on the PATTERN, not one field: each looked-up property the neutral panel
    seeds a bounded widget with must pass through `clamp_seed` first. Catches a future
    field added without the guard."""
    import inspect
    lines = inspect.getsource(SB._neutral_panel).splitlines()
    for prop in ("log_kow", "MW", "K_AW", "pKa", "Koc"):
        seed = f'found.get("{prop}"'
        hits = [i for i, l in enumerate(lines) if seed in l]
        assert hits, f"{prop} is no longer seeded from the lookup — update this test"
        # the clamp may sit a line or two below the seed (Koc derives a default first),
        # so look at a small window rather than the one line
        window = "\n".join(lines[hits[0]:hits[0] + 4])
        assert "clamp_seed" in window, (
            f"the {prop} seed reaches a widget without clamp_seed — an out-of-range "
            f"lookup would abort the whole app")


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
