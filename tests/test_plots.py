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


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))
