# -*- coding: utf-8 -*-
"""Inline-SVG rice-plant accumulation map (no Plotly).

A clean rice-plant glyph rendered as a self-contained SVG string, styled to the
shared design: a calm warm-paper canvas + a shallow soil band (no vivid sky/heavy
soil), each organ (roots+seed / leaf blades / central culm / grain ear) shaded by
a heat colormap of its PFAS value, white label cards with a colour swatch + value
+ a "최고" chip on the top-accumulating organ joined by leader lines, a rounded
concentration legend and a dark pore-water tag.

The junction circle between stem and root is filled with the BACKGROUND colour
(no outline), sitting in front of the roots and behind the stem, so it just masks
the overlap. No compartment carries an outline.

Pure functions (no Streamlit). `plant_svg_from_res(res, t_index)` builds the values
from a `model_api` result — the SVG analogue of `plots.fig_schematic_from_res`.
"""
from __future__ import annotations

import numpy as np

import model_api as api

# heat ramp (cream → wheat → gold → ochre → terracotta), matches plots._HEAT / design
_STOPS = [(0.0, (244, 247, 238)), (0.22, (234, 223, 155)), (0.5, (231, 178, 76)),
          (0.75, (217, 138, 68)), (1.0, (180, 86, 46))]

# design palette
_SKY = "#8FCFEE"         # canvas background (sky)
_PAPER = "#F7F3EC"       # (kept for reference / paper variant)
_SOIL = "#B98A4E"        # soil band (deeper earth brown)
_SOIL_DK = "#95672F"
_SPECK = "#6F4A24"
_INK = "#211E18"
_MUTED = "#5C554A"
_BORDER = "#E4DCCE"
_LEADER = "#C9B48C"
_PWBG = "#4A3520"        # pore-water chip
_PWTX = "#F3EAD4"
_CHIP_TOP_BG, _CHIP_TOP_TX = "#FBE4DE", "#B23A2E"     # "최고"
_CHIP_SAFE_BG, _CHIP_SAFE_TX = "#DCEFE6", "#0E6B4F"   # "안전"
_CHIP_WARN_BG, _CHIP_WARN_TX = "#FAEBD1", "#9A5A00"
_CHIP_DANG_BG, _CHIP_DANG_TX = "#FADEDA", "#B23A2E"
_FONT = "'Malgun Gothic','맑은 고딕',sans-serif"
_MONO = "'DM Mono','SFMono-Regular',Consolas,'Malgun Gothic',monospace"

_LAB = {"ko": {"root": "뿌리", "stem": "줄기", "leaf": "잎", "grain": "낟알",
               "conc": "농도", "high": "높음", "low": "낮음", "pw": "공극수", "top": "최고",
               "safe": "안전", "warn": "주의", "dang": "초과"},
        "en": {"root": "Roots", "stem": "Stem", "leaf": "Leaf", "grain": "Grain",
               "conc": "level", "high": "high", "low": "low", "pw": "pore water", "top": "top",
               "safe": "safe", "warn": "caution", "dang": "over"}}
_SIGCHIP = {"green": ("safe", _CHIP_SAFE_BG, _CHIP_SAFE_TX),
            "amber": ("warn", _CHIP_WARN_BG, _CHIP_WARN_TX),
            "red": ("dang", _CHIP_DANG_BG, _CHIP_DANG_TX)}


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


def _leader(x1, y1, x2, y2):
    return (f'<path d="M{x1},{y1} L{x2},{y2}" fill="none" stroke="{_LEADER}" '
            f'stroke-width="1.4"/>')


def _card(x, y, name, val, swatch, chip=None):
    """White label card (colour swatch + name + mono value + optional chip). x,y = top-left."""
    w, h = 168, 48
    out = [f'<g><rect x="{x}" y="{y}" width="{w}" height="{h}" rx="12" fill="#FFFFFF" '
           f'stroke="{_BORDER}" stroke-width="1"/>',
           f'<rect x="{x+10}" y="{y+11}" width="7" height="{h-22}" rx="3.5" fill="{swatch}"/>',
           f'<text x="{x+26}" y="{y+20}" font-size="13" font-weight="700" fill="{_INK}" '
           f'font-family="{_FONT}">{name}</text>',
           f'<text x="{x+26}" y="{y+38}" font-size="14" fill="{_MUTED}" '
           f'font-family="{_MONO}">{val}</text>']
    if chip:
        label, cbg, ctx = chip
        cw = 34 + 8 * max(0, len(label) - 2)
        out.append(f'<rect x="{x+w-cw-10}" y="{y+9}" width="{cw}" height="18" rx="9" fill="{cbg}"/>')
        out.append(f'<text x="{x+w-10-cw/2:.0f}" y="{y+22}" font-size="11" font-weight="700" '
                   f'fill="{ctx}" text-anchor="middle" font-family="{_FONT}">{label}</text>')
    out.append('</g>')
    return "".join(out)


def plant_svg(values, *, cmin, cmax, cwo=None, lang="ko", labels=True,
              grain_signal=None, bg=_SKY, W=760, H=560):
    """Return a self-contained SVG string of the rice-plant map (design style)."""
    L = _LAB.get(lang, _LAB["ko"])
    root = values.get("root")
    stem = values.get("stem", values.get("straw"))
    leaf = values.get("leaf", values.get("straw"))
    grain = values.get("grain")
    cr, cs, cl, cg = (_col(root, cmin, cmax), _col(stem, cmin, cmax),
                      _col(leaf, cmin, cmax), _col(grain, cmin, cmax))
    # which organ accumulates most (gets the "최고" chip)
    named = {"root": root, "stem": stem, "leaf": leaf, "grain": grain}
    finite = {k: float(v) for k, v in named.items() if v is not None and np.isfinite(v)}
    top_organ = max(finite, key=finite.get) if finite else None

    tx, ty = 120, 70                       # plant translate (native 512-box → frame)
    soil_top = 476
    P = []
    P.append(f'<rect x="0" y="0" width="{W}" height="{H}" rx="20" fill="{bg}"/>')
    # shallow soil band (no outline) + a few faint speckles; roots dip in
    P.append(f'<path d="M0,{soil_top} Q190,{soil_top-8} 380,{soil_top-2} '
             f'Q560,{soil_top+4} {W},{soil_top-6} L{W},{H} L0,{H} Z" fill="{_SOIL}"/>')
    P.append(f'<path d="M0,{soil_top+34} Q220,{soil_top+26} 430,{soil_top+32} '
             f'Q600,{soil_top+38} {W},{soil_top+30} L{W},{H} L0,{H} Z" fill="{_SOIL_DK}" opacity="0.6"/>')
    for sx, sy, r in [(150, 520, 3), (240, 536, 2.4), (470, 520, 3), (610, 540, 2.4), (330, 545, 2.4)]:
        P.append(f'<ellipse cx="{sx}" cy="{sy}" rx="{r}" ry="{r*0.7:.1f}" fill="{_SPECK}" opacity="0.5"/>')

    # ---- the plant (native 512 coords; NO outlines; junction circle = bg) ----
    P.append(f'<g transform="translate({tx},{ty})">')
    P.append(f'<g stroke="{cr}" stroke-width="18" stroke-linecap="round" fill="none">'
             f'<path d="M246,394 A42,42 0 0 1 198.9,435.7"/>'
             f'<path d="M266,394 A42,42 0 0 0 313.1,435.7"/>'
             f'<path d="M256,395 L256,452"/></g>')
    # junction circle: BACKGROUND colour, no outline, in front of roots / behind stem
    P.append(f'<circle cx="256" cy="378" r="24" fill="{bg}"/>')
    P.append(f'<g stroke="{cs}" stroke-width="30" stroke-linecap="round" fill="none">'
             f'<path d="M256,378 L256,170 A80,80 0 0 1 296,100.7"/></g>')
    P.append(f'<g stroke="{cl}" stroke-width="30" stroke-linecap="round" fill="none">'
             f'<path d="M178.1,172.6 L216.9,211.4"/><path d="M333.3,215.5 L294.5,254.3"/>'
             f'<path d="M178.3,258.2 L217.1,297.0"/><path d="M333.3,296.7 L294.5,335.5"/></g>')
    P.append(f'<g fill="{cg}">'
             f'<ellipse cx="349.3" cy="63.5" rx="28.8" ry="18.7" transform="rotate(-22 349.3 63.5)"/>'
             f'<ellipse cx="398.7" cy="101.3" rx="28.8" ry="18.7" transform="rotate(27 398.7 101.3)"/>'
             f'<ellipse cx="352.3" cy="123.1" rx="28.8" ry="18.7" transform="rotate(52 352.3 123.1)"/>'
             f'<ellipse cx="409.9" cy="157.3" rx="28.8" ry="18.7" transform="rotate(51 409.9 157.3)"/></g>')
    P.append('</g>')

    # ---- leader lines ON TOP of the plant, ending at each organ's near edge with
    #      a small connector dot (so they read as clean connectors, not crossing) ----
    if labels:
        ends = {"grain": (472, 190), "leaf": (298, 328), "stem": (361, 372), "root": (368, 520)}
        starts = {"grain": (198, 174), "leaf": (198, 306), "stem": (198, 388), "root": (198, 500)}
        for k in ("grain", "leaf", "stem", "root"):
            (sx, sy), (ex, ey) = starts[k], ends[k]
            P.append(_leader(sx, sy, ex, ey))
            P.append(f'<circle cx="{ex}" cy="{ey}" r="3" fill="{_LEADER}"/>')

    # ---- label cards (on top) ----
    if labels:
        def chip_for(key, color_v):
            if key == top_organ:
                return (L["top"], _CHIP_TOP_BG, _CHIP_TOP_TX)
            if key == "grain" and grain_signal in _SIGCHIP:
                nm, cbg, ctx = _SIGCHIP[grain_signal]
                return (L[nm], cbg, ctx)
            return None
        P.append(_card(30, 150, L["grain"], _fmt(grain), cg, chip_for("grain", cg)))
        P.append(_card(30, 282, L["leaf"], _fmt(leaf), cl, chip_for("leaf", cl)))
        P.append(_card(30, 364, L["stem"], _fmt(stem), cs, chip_for("stem", cs)))
        P.append(_card(30, 476, L["root"], _fmt(root), cr, chip_for("root", cr)))

    # ---- pore-water tag (dark chip) ----
    if cwo is not None and np.isfinite(cwo):
        P.append(f'<g><rect x="{W-232}" y="{H-70}" width="196" height="42" rx="14" fill="{_PWBG}"/>'
                 f'<text x="{W-134}" y="{H-52}" text-anchor="middle" font-size="12.5" '
                 f'font-weight="700" fill="{_PWTX}" font-family="{_FONT}">{L["pw"]}</text>'
                 f'<text x="{W-134}" y="{H-36}" text-anchor="middle" font-size="12" '
                 f'fill="{_PWTX}" font-family="{_MONO}">PFAS = {cwo:.3g} µg/L</text></g>')

    # ---- concentration legend inside a rounded white card (right) ----
    lg = "".join(f'<stop offset="{p*100:.0f}%" stop-color="rgb{c}"/>' for p, c in _STOPS)
    gx = W - 48                                        # gradient x inside the card
    P.append(f'<defs><linearGradient id="pfaslg" x1="0" y1="1" x2="0" y2="0">{lg}</linearGradient></defs>')
    P.append(f'<rect x="{W-98}" y="160" width="82" height="282" rx="16" fill="#FFFFFF" '
             f'stroke="{_BORDER}" stroke-width="1"/>')
    P.append(f'<text x="{gx+7}" y="186" font-size="12" fill="{_MUTED}" text-anchor="middle" '
             f'font-family="{_FONT}">{L["conc"]}</text>')
    P.append(f'<rect x="{gx}" y="206" width="14" height="200" rx="7" fill="url(#pfaslg)"/>')
    P.append(f'<text x="{gx-9}" y="212" font-size="11" fill="{_MUTED}" text-anchor="end" '
             f'font-family="{_FONT}">{L["high"]}</text>')
    P.append(f'<text x="{gx-9}" y="406" font-size="11" fill="{_MUTED}" text-anchor="end" '
             f'font-family="{_FONT}">{L["low"]}</text>')

    # rounded on ALL corners (clip the square-cornered soil band too), and fill the
    # container width (no centred max-width cap).
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" '
            f'preserveAspectRatio="xMidYMid meet" font-family="{_FONT}" '
            f'style="width:100%;height:auto;display:block">'
            f'<defs><clipPath id="pfascard"><rect x="0" y="0" width="{W}" height="{H}" rx="24"/></clipPath></defs>'
            f'<g clip-path="url(#pfascard)">' + "".join(P) + '</g></svg>')


def plant_svg_from_res(res, t_index=-1, *, lang="ko", labels=True, bg=_SKY, grain_signal=None):
    """Build the SVG plant map from a `model_api` result at one time index."""
    sv = api.schematic_values(res, "conc", t_index)
    return plant_svg(sv["values"], cmin=sv["cmin"], cmax=sv["cmax"], cwo=sv.get("Cwo"),
                     lang=lang, labels=labels, bg=bg, grain_signal=grain_signal)


if __name__ == "__main__":
    r = api.simulate("PFOA", Cwo=1.0, measured_forcing=True, biomass="oryza")
    open("plant_svg_out.svg", "w").write(plant_svg_from_res(r, -1, lang="ko"))
    print("wrote plant_svg_out.svg")
