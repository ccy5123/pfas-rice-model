"""The Streamlit live-HYDRUS path passes only inputs the soil model takes (PREREG 9.2-9.13 fixed the water
calendar and the percolation; the old flood/percolation sliders raised a TypeError), and the sidebar no longer
offers the removed placeholder crop logistics (`simulate(measured_forcing=False)` refuses them)."""
import inspect
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT)]
import model_api as api        # noqa: E402
import soil_hydrus as sh       # noqa: E402

pytest.importorskip("streamlit")
from ui.common import _hydrus_drivers_cached, crop_kw   # noqa: E402


def test_sidebar_offers_no_removed_input():
    src = (ROOT / "ui" / "sidebar.py").read_text()
    for gone in ("flood_until", "Percolation excess", "Crop drivers on (off"):
        assert gone not in src


def test_cached_helper_takes_only_soil_model_inputs():
    params = inspect.signature(_hydrus_drivers_cached).parameters
    assert "flood_until" not in params and "percolation" not in params
    assert "after_drainage" in inspect.signature(sh.run_paddy_hydrus).parameters


@pytest.mark.skipif(not api.hydrus_available(), reason="HYDRUS-1D engine not built")
@pytest.mark.parametrize("after", ["supplied", "rain"])
def test_the_sidebar_call_runs(after):
    d, _ = api.hydrus_drivers("PFOA", season=None, f_oc=sh.FOC, Kd=None, after_drainage=after,
                              **crop_kw("site2020"))
    assert len(d["t"]) == len(d["Cwo"]) > 10


def test_placeholder_logistics_stay_refused():
    with pytest.raises(ValueError):
        api.simulate("PFOA", measured_forcing=False)
