"""One-page HTML run summary (slide / handout download) -- UI-agnostic.

Moved out of ``ui/common.py`` so the Streamlit app and the FastAPI web backend
(``server/``) build the SAME file: the plant accumulation map, the EFSA-TWI
intake gauge, the key tissue numbers (with the coarse a-priori band) and the
caveats/sources, as a single offline HTML file (Plotly embedded, so it opens
without a network). Pure: no Streamlit import.
"""
from __future__ import annotations

import numpy as np

import model_api as api
import plots


def summary_html(res, *, congener, disclaimer, scenario=None, lang="ko"):
    """Return the one-page summary as UTF-8 ``bytes``.

    ``disclaimer`` is the (markdown) disclaimer text for ``lang``; ``**`` emphasis
    is stripped. ``scenario`` is an optional short scenario label. Raises on a
    malformed ``res`` -- callers decide how to degrade.
    """
    ko = lang != "en"
    cong = congener
    grain_c = float(res["conc"]["grain"][-1])
    info = api.intake_fraction(grain_c, congener=cong)
    map_fig = plots.fig_schematic_from_res(res, "conc", -1, lang=lang)
    gauge = plots.fig_intake_gauge(info, lang=lang)
    # embed plotly.js once (offline-safe), then the gauge without re-embedding
    map_div = map_fig.to_html(full_html=False, include_plotlyjs=True,
                              default_height="480px")
    gauge_div = gauge.to_html(full_html=False, include_plotlyjs=False,
                              default_height="300px")
    rows = [("뿌리" if ko else "Roots", float(res["conc"]["root"][-1])),
            ("짚(줄기+잎)" if ko else "Straw", float(res["straw"][-1])),
            ("낟알(먹는 쌀)" if ko else "Grain", grain_c)]
    trs = ""
    for name, v in rows:
        b = api.predictive_band(v)
        trs += (f"<tr><td>{name}</td><td>{v:.2g}</td>"
                f"<td>{b['lo']:.2g}–{b['hi']:.2g}</td></tr>")
    efsa = ""
    if np.isfinite(info.get("percent", float("nan"))):
        grp = ("" if info.get("in_group")
               else ("(참고 — EFSA 4종 합 미포함)" if ko else "(reference — not in the EFSA four)"))
        efsa = (f"<p><b>EFSA 안전기준(TWI) 대비 약 {info['percent']:.0f}%</b> {grp} "
                f"— 하루 {info['rice_intake_g_day']:.0f} g 쌀 섭취 가정.</p>" if ko else
                f"<p><b>~{info['percent']:.0f}% of the EFSA TWI</b> {grp} "
                f"— assuming {info['rice_intake_g_day']:.0f} g rice/day.</p>")
    disc = disclaimer.replace("**", "")
    title = (f"PFAS–벼 흡수 요약 — {cong}" if ko else f"PFAS–Rice uptake summary — {cong}")
    scen_line = (f"<p style='color:#555'>시나리오: {scenario}</p>" if (ko and scenario)
                 else (f"<p style='color:#555'>Scenario: {scenario}</p>" if scenario else ""))
    th = ("부위,예측 농도 (µg/kg),대략 범위" if ko else "Tissue,Predicted (µg/kg),Rough range").split(",")
    src = ("출처: EFSA 2020 그룹 TWI 4.4 ng/kg 체중/주 (doi:10.2903/j.efsa.2020.6223) · "
           "쌀 섭취량 KOSIS/국민건강영양조사 · 실측 Yamazaki et al. 2023." if ko else
           "Sources: EFSA 2020 TWI 4.4 ng/kg bw/week (doi:10.2903/j.efsa.2020.6223); "
           "rice intake KOSIS/KNHANES; observed Yamazaki et al. 2023.")
    html = (
        f"<!doctype html><html lang='{'ko' if ko else 'en'}'><head><meta charset='utf-8'>"
        f"<meta name='viewport' content='width=device-width,initial-scale=1'>"
        f"<title>{title}</title><style>"
        "body{font-family:'Malgun Gothic',-apple-system,Segoe UI,Arial,sans-serif;"
        "max-width:960px;margin:24px auto;padding:0 16px;color:#222}"
        "h1{font-size:1.5rem;margin:.2rem 0}"
        "table{border-collapse:collapse;margin:8px 0}"
        "td,th{border:1px solid #ccc;padding:5px 12px;text-align:left}"
        "th{background:#f5f5f5}"
        ".row{display:flex;gap:24px;flex-wrap:wrap;align-items:flex-start}"
        ".caveat{background:#fff8e1;border-left:4px solid #f0ad4e;padding:10px 14px;"
        "border-radius:4px;margin:12px 0;font-size:.92rem}"
        ".src{color:#666;font-size:.82rem;margin-top:10px}</style></head><body>"
        f"<h1>🌾 {title}</h1>{scen_line}"
        f"<div class='row'><div style='flex:2;min-width:340px'>{map_div}</div>"
        f"<div style='flex:1;min-width:280px'>{gauge_div}"
        f"<table><tr><th>{th[0]}</th><th>{th[1]}</th><th>{th[2]}</th></tr>{trs}</table>"
        f"{efsa}</div></div>"
        f"<div class='caveat'>⚠ {disc}</div>"
        f"<p class='src'>{src}</p></body></html>")
    return html.encode("utf-8")
