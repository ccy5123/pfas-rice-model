"""Head-less tests for the Plotly figure builders (src/plots.py)."""
import pytest

pytest.importorskip("plotly")          # app-only dependency

import plotly.graph_objects as go      # noqa: E402
import model_api as api                # noqa: E402
import plots                           # noqa: E402


def test_all_figures_build():
    res = api.simulate("PFOA")
    obs = api.observed_baf("PFOA")
    rows = api.chain_table()
    figs = [
        plots.fig_tissue(res),
        plots.fig_baf(res, obs),
        plots.fig_baf(res, {}),                       # no-obs branch
        plots.fig_chain(rows, "PFOA", "K_PL"),
        plots.fig_chain(rows, "PFOS", "f_xy_recommended"),
        plots.fig_forcings(res["t"], 120.0),
        plots.fig_compare({n: api.simulate(n) for n in ("PFBA", "PFOA", "PFOS")}, "straw"),
    ]
    for f in figs:
        assert isinstance(f, go.Figure)
        assert len(f.data) >= 1
        assert f.layout.title.text


def test_chain_log_axis_for_partition_keys():
    rows = api.chain_table()
    assert plots.fig_chain(rows, "PFOA", "K_PL").layout.yaxis.type == "log"
    assert plots.fig_chain(rows, "PFOA", "K_cw_root").layout.yaxis.type == "linear"


def test_plant_schematic_and_soil_figures_build():
    import numpy as np
    res = api.simulate("PFOA")
    figs = [
        plots.fig_schematic_from_res(res, "conc", -1),
        plots.fig_schematic_from_res(res, "baf", 0),
        plots.fig_soil_profile(res),
        plots.fig_drivers(res),
    ]
    for f in figs:
        assert isinstance(f, go.Figure) and len(f.data) >= 1 and f.layout.title.text
    # the schematic draws the rice silhouette as layout shapes
    sch = plots.fig_schematic_from_res(res, "conc", -1)
    assert len(sch.layout.shapes) > 20            # soil + roots + leaves + culms + grain beads


def test_schematic_animated_has_frames():
    res = api.simulate("PFOA")
    fig = plots.fig_schematic_animated(res, "conc", n_frames=8)
    assert len(fig.frames) >= 2
    assert fig.layout.updatemenus and fig.layout.sliders


def test_schematic_handles_straw_only_values():
    # biomonitoring case: only root/straw/grain measured (no stem/leaf)
    vals = {"root": 0.49, "straw": 0.83, "grain": 0.46}
    fig = plots.fig_plant_schematic(vals, cmin=0.0, cmax=0.83, label="conc", Cwo=1.0)
    assert isinstance(fig, go.Figure) and len(fig.layout.shapes) > 20


def test_soil_profile_heatmap_and_isotherm():
    import numpy as np
    res = api.simulate("PFOA")
    prof = dict(time=res["t"], depth=np.linspace(0, 40, 9),
                conc=np.outer(np.linspace(1, 0.5, 9), np.ones_like(res["t"])))
    fig = plots.fig_soil_profile(res, profile=prof)
    assert any(d.type == "heatmap" for d in fig.data)
    _, soil = api.pore_water_from_inventory(res["t"], 5.0)
    iso = plots.fig_isotherm(soil, Cwo_now=1.0)
    assert isinstance(iso, go.Figure) and len(iso.data) >= 1


def test_biomon_compare_builds():
    mb = {"root": 0.49, "straw": 0.83, "grain": 0.46}
    fig = plots.fig_biomon_compare(mb, {"root": 0.47, "straw": 0.9, "grain": 0.15})
    assert isinstance(fig, go.Figure) and len(fig.data) == 2


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))


def test_fig_tang_tf_builds():
    import model_api as api, plots
    val = api.tang_tf_validation("PFOA")
    valr = api.tang_tf_validation("PFOA", use_refit=True)
    fig = plots.fig_tang_tf(val, valr)
    assert len(fig.data) == 3                       # Tang + model + refit
    assert fig.layout.yaxis.type == "log"
    # single-arg form (no refit bar) also builds
    assert len(plots.fig_tang_tf(val).data) == 2


def test_fig_mass_builds():
    import model_api as api, plots
    fig = plots.fig_mass(api.simulate("PFOA", biomass="oryza"))
    assert len(fig.data) == 5                                    # 4 organs + whole plant
    assert "kg" in fig.layout.yaxis.title.text
