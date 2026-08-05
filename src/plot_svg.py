# -*- coding: utf-8 -*-
"""Inline-SVG rice-plant accumulation map (no Plotly).

A clean, bold rice-plant glyph rendered as a self-contained SVG string, with each
organ (roots+seed / leaf blades / central culm / grain ear) shaded by a heat
colormap of its PFAS value. SVG gives true round stroke caps and pixel-crisp
vectors, so it matches the icon design far better than Plotly path shapes.

Pure functions (no Streamlit import): `plant_svg(...)` takes the per-organ values,
`plant_svg_from_res(res, t_index)` builds them from a `model_api` result — the SVG
analogue of `plots.fig_schematic_from_res`. The app renders the returned string
with `st.components.v1.html`.
"""
from __future__ import annotations

import numpy as np

import model_api as api

# warm agricultural heat ramp (matches plots._HEAT): cream→wheat→gold→ochre→terracotta
_STOPS = [(0.0, (244, 247, 238)), (0.22, (234, 223, 155)), (0.5, (231, 178, 76)),
          (0.75, (217, 138, 68)), (1.0, (180, 86, 46))]
# a more saturated sky background for the plant card
_SKY = "#8FCFEE"
_SOIL = "#7A5A34"
_FONT = "'Malgun Gothic','맑은 고딕',sans-serif"

_LAB = {"ko": {"root": "뿌리", "stem": "줄기", "leaf": "잎", "grain": "낟알",
               "conc": "농도", "high": "높음", "low": "낮음", "pw": "토양수"},
        "en": {"root": "Roots", "stem": "Stem", "leaf": "Leaf", "grain": "Grain",
               "conc": "level", "high": "high", "low": "low", "pw": "pore water"}}


def _heat(f):
    f = max(0.0, min(1.0, f))
    for i in range(len(_STOPS) - 1):
        a, ca = _STOPS[i]
        b, cb = _STOPS[i + 1]
        if f <= b:
            t = 0.0 if b == a else (f - a) / (b - a)
            r = [round(ca[k] + (cb[k] - ca[k]) * t) for k in range(3)]
            return f"rgb({r[0]},{r[1]},{r[2]})"
    return "rgb(180,86,46)"


def _frac(v, cmin, cmax):
    if v is None or not np.isfinite(v):
        return None
    return 0.5 if cmax <= cmin else float(np.clip((v - cmin) / (cmax - cmin), 0.0, 1.0))


def _col(v, cmin, cmax):
    f = _frac(v, cmin, cmax)
    return "#dcdcdc" if f is None else _heat(f)


def _fmt(v):
    return "—" if (v is None or not np.isfinite(v)) else f"{v:.3g}"


def plant_svg(values, *, cmin, cmax, cwo=None, lang="ko", labels=True,
              bg=_SKY, W=520, H=560):
    """Return a self-contained SVG string of the rice-plant map (no outlines)."""
    L = _LAB.get(lang, _LAB["ko"])
    root = values.get("root")
    stem = values.get("stem", values.get("straw"))
    leaf = values.get("leaf", values.get("straw"))
    grain = values.get("grain")
    cr, cs, cl, cg = (_col(root, cmin, cmax), _col(stem, cmin, cmax),
                      _col(leaf, cmin, cmax), _col(grain, cmin, cmax))

    def lab(x, y, key, val, anchor="start"):                 # annotation box (no outline)
        if not labels:
            return ""
        bx = x - 4 if anchor == "start" else x - 88
        tx = x if anchor == "start" else x - 84
        return (f'<g><rect x="{bx}" y="{y-15}" width="92" height="30" rx="8" '
                f'fill="rgba(255,255,255,0.86)"/>'
                f'<text x="{tx}" y="{y-2}" font-size="13" font-weight="700" fill="#2f3b33">{L[key]}</text>'
                f'<text x="{tx}" y="{y+12}" font-size="12" fill="#586257">{_fmt(val)}</text></g>')

    lg = "".join(f'<stop offset="{p*100:.0f}%" stop-color="rgb{c}"/>' for p, c in _STOPS)
    soil = f'<path d="M0,470 Q130,462 200,468 Q290,474 520,468 L520,{H} L0,{H} Z" fill="{_SOIL}"/>'
    cwo_txt = ""
    if cwo is not None and np.isfinite(cwo):
        cwo_txt = (f'<g><rect x="332" y="500" width="152" height="34" rx="9" fill="rgba(70,48,26,0.92)"/>'
                   f'<text x="408" y="514" text-anchor="middle" font-size="12" fill="#f3ead4">{L["pw"]}</text>'
                   f'<text x="408" y="528" text-anchor="middle" font-size="12" fill="#f3ead4">PFAS={cwo:.3g} µg/L</text></g>')

    # The icon is translated DOWN so the roots sit inside the soil; the junction
    # circle is background-coloured (in front of roots / behind stem) to mask the
    # overlap. No element carries an outline.
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="100%" '
        f'height="{H}" font-family="{_FONT}" style="max-width:{W}px;display:block;margin:auto">'
        f'<rect width="{W}" height="{H}" rx="20" fill="{bg}"/>'
        f'{soil}'
        f'<g transform="translate(4,70)">'
        f'<g stroke="{cr}" stroke-width="18" stroke-linecap="round" fill="none">'
        f'<path d="M246,394 A42,42 0 0 1 198.9,435.7"/>'
        f'<path d="M266,394 A42,42 0 0 0 313.1,435.7"/>'
        f'<path d="M256,395 L256,452"/></g>'
        f'<circle cx="256" cy="378" r="24" fill="{bg}"/>'
        f'<g stroke="{cs}" stroke-width="30" stroke-linecap="round" fill="none">'
        f'<path d="M256,378 L256,170 A80,80 0 0 1 296,100.7"/></g>'
        f'<g stroke="{cl}" stroke-width="30" stroke-linecap="round" fill="none">'
        f'<path d="M178.1,172.6 L216.9,211.4"/><path d="M333.3,215.5 L294.5,254.3"/>'
        f'<path d="M178.3,258.2 L217.1,297.0"/><path d="M333.3,296.7 L294.5,335.5"/></g>'
        f'<g fill="{cg}">'
        f'<ellipse cx="349.3" cy="63.5" rx="28.8" ry="18.7" transform="rotate(-22 349.3 63.5)"/>'
        f'<ellipse cx="398.7" cy="101.3" rx="28.8" ry="18.7" transform="rotate(27 398.7 101.3)"/>'
        f'<ellipse cx="352.3" cy="123.1" rx="28.8" ry="18.7" transform="rotate(52 352.3 123.1)"/>'
        f'<ellipse cx="409.9" cy="157.3" rx="28.8" ry="18.7" transform="rotate(51 409.9 157.3)"/></g>'
        f'</g>'
        f'{cwo_txt}'
        f'{lab(150, 150, "grain", grain, "end")}'
        f'{lab(60, 300, "leaf", leaf)}'
        f'{lab(60, 392, "stem", stem)}'
        f'{lab(150, 520, "root", root, "end")}'
        f'<defs><linearGradient id="lg" x1="0" y1="1" x2="0" y2="0">{lg}</linearGradient></defs>'
        f'<rect x="490" y="170" width="16" height="230" rx="6" fill="url(#lg)"/>'
        f'<text x="498" y="160" font-size="12" fill="#33414b" text-anchor="middle">{L["conc"]}</text>'
        f'<text x="486" y="176" font-size="11" fill="#33414b" text-anchor="end">{L["high"]}</text>'
        f'<text x="486" y="400" font-size="11" fill="#33414b" text-anchor="end">{L["low"]}</text>'
        f'</svg>')


def plant_svg_from_res(res, t_index=-1, *, lang="ko", labels=True, bg=_SKY):
    """Build the SVG plant map from a `model_api` result at one time index."""
    sv = api.schematic_values(res, "conc", t_index)
    return plant_svg(sv["values"], cmin=sv["cmin"], cmax=sv["cmax"],
                     cwo=sv.get("Cwo"), lang=lang, labels=labels, bg=bg)


if __name__ == "__main__":
    r = api.simulate("PFOA", Cwo=1.0, measured_forcing=True, biomass="oryza")
    open("plant_svg_out.svg", "w").write(plant_svg_from_res(r, -1, lang="ko"))
    print("wrote plant_svg_out.svg")
