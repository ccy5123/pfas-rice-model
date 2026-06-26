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


def test_fig_baf_extra_overlay():
    """The optional `extra` overlay (e.g. the two-pool seq model) adds one bar per
    series alongside the core model and observed; absent -> backward-compatible."""
    res = api.simulate("PFUnDA")
    obs = api.observed_baf("PFUnDA")
    base = plots.fig_baf(res, obs)
    over = plots.fig_baf(res, obs, extra={"two-pool (seq)": {"root": 15.8, "straw": 6.8, "grain": 6.5}})
    names = [tr.name for tr in over.data]
    assert "two-pool (seq)" in names and "model (4-pool core)" in names and "Yamazaki 2023" in names
    assert len(over.data) == len(base.data) + 1               # exactly one extra series


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


def test_fig_burden_builds():
    import model_api as api, plots
    fig = plots.fig_burden(api.simulate("PFOA", biomass="oryza"))
    assert len(fig.data) == 5                                    # 4 organs + whole plant
    assert "µg" in fig.layout.yaxis.title.text                  # PFAS mass (burden), not biomass


def test_fig_where_plain_band_error_bars():
    """band=True overlays the a-priori predictive uncertainty as asymmetric error
    bars (×/÷ ~7); default off keeps the bar bare (backward-compatible)."""
    res = api.simulate("PFOA")
    bare = plots.fig_where_plain(res)
    assert bare.data[0].error_y.array is None       # no band by default
    banded = plots.fig_where_plain(res, lang="ko", band=True)
    ey = banded.data[0].error_y
    assert ey.array is not None and ey.symmetric is False
    # the band matches model_api.predictive_band for each bar (root/straw/grain)
    vals = [res["conc"]["root"][-1], res["straw"][-1], res["conc"]["grain"][-1]]
    for i, v in enumerate(vals):
        b = api.predictive_band(float(v))
        assert ey.array[i] == pytest.approx(b["hi"] - float(v))
        assert ey.arrayminus[i] == pytest.approx(float(v) - b["lo"])


def test_fig_baf_korean_variant():
    """lang='ko' localises the BAF axis/title/legend; English stays the default."""
    res = api.simulate("PFOA")
    obs = api.observed_baf("PFOA")
    ko = plots.fig_baf(res, obs, lang="ko")
    assert "축적 배수" in ko.layout.yaxis.title.text and "BAF" not in ko.layout.yaxis.title.text
    assert list(ko.data[0].x) == ["뿌리", "짚", "낟알"]
    names = {tr.name for tr in ko.data}
    assert "모델" in names and "Yamazaki 2023 (실측)" in names
    # English default unchanged
    en = plots.fig_baf(res, obs)
    assert en.layout.yaxis.title.text == "BAF [L/kg]"
    assert "model (4-pool core)" in {tr.name for tr in en.data}


def test_plain_language_figures_build():
    """Simple-mode plain-language builders use friendly tissue names + jargon-free
    titles (no 'BAF'/'Cwᵒ'/'f_xy' symbols leaking into the general-audience view)."""
    import model_api as api, plots
    res = api.simulate("PFOA")
    bld = plots.fig_buildup_plain(res)
    whr = plots.fig_where_plain(res)
    for f in (bld, whr):
        assert isinstance(f, go.Figure) and len(f.data) >= 1 and f.layout.title.text
        for sym in ("BAF", "Cwᵒ", "f_xy", "eᴺ", "µg/kg"):     # no expert jargon in titles
            assert sym not in f.layout.title.text
    # build-up uses friendly tissue names, not the raw symbols
    names = {tr.name for tr in bld.data}
    assert {"Roots", "Stems", "Leaves", "Grain"} <= names
    # the "where it ends up" bar has the three plant parts (friendly labels)
    assert list(whr.data[0].x) == ["Roots", "Straw", "Grain"]


def test_fig_exposure_posterior_builds():
    """The Bayesian inverse posterior plot builds with a log-x density + 95% band."""
    import model_api as api, plots
    r = api.simulate("PFOA", Cwo=1.5)
    est = api.estimate_exposure_bayesian("PFOA", {"root": r["conc"]["root"][-1],
                                                  "grain": r["conc"]["grain"][-1]})
    fig = plots.fig_exposure_posterior(est)
    assert isinstance(fig, go.Figure) and len(fig.data) >= 1
    assert fig.layout.xaxis.type == "log"
    assert "soil water" in fig.layout.xaxis.title.text
    assert fig.layout.title.text


def test_plain_figures_korean_variant():
    """lang='ko' renders Korean labels in the Simple-mode builders; English stays default."""
    import model_api as api, plots
    res = api.simulate("PFOA")
    bld = plots.fig_buildup_plain(res, lang="ko")
    assert {tr.name for tr in bld.data} >= {"뿌리", "줄기", "잎", "낟알"}
    assert "축적" in bld.layout.title.text
    whr = plots.fig_where_plain(res, lang="ko")
    assert list(whr.data[0].x) == ["뿌리", "짚", "낟알"]
    mp = plots.fig_schematic_from_res(res, "conc", -1, lang="ko")
    assert "지도" in mp.layout.title.text
    est = api.estimate_exposure_bayesian("PFOA", {"grain": res["conc"]["grain"][-1]})
    assert "토양수" in plots.fig_exposure_posterior(est, lang="ko").layout.xaxis.title.text
    # English default unchanged
    assert "Roots" in {tr.name for tr in plots.fig_buildup_plain(res).data}


def test_fig_cwo_profile_builds():
    import numpy as np, plots
    fig = plots.fig_cwo_profile("PFBA", level=1.0, profile="flooded")
    assert [d.name for d in fig.data] == ["constant", "flooded"]
    flat, shaped = np.array(fig.data[0].y), np.array(fig.data[1].y)
    assert np.allclose(flat, 1.0)                                # constant baseline is flat
    assert shaped[-1] < shaped[0]                                # short chain leaches (declines)
    # long chain stays buffered (~flat); SMILES fallback (congener=None) still builds
    long = np.array(plots.fig_cwo_profile("PFDoDA", profile="flooded").data[1].y)
    assert long[-1] > 0.95 * long[0]
    assert len(plots.fig_cwo_profile(None, profile="flooded").data) == 2
