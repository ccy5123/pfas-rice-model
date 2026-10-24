"""Overlay check of raw_si/gao2021_fig3_digitised.csv against the original Gao et al. 2021 Fig. 3.

The figure is the only bitmap on page 8 of the article PDF (doi 10.1016/j.agwat.2021.106864; the PDF is
not in the repo - copyright).

    python docs/literature_db/tools/gao2021_fig3_check.py <gao2021.pdf | fig3.png> [--out DIR]

(1) re-fits the axis calibration on every tick mark and compares it with the constants the CSV uses;
(2) checks marker centres against the error bars: an SD bar is symmetric, so the midpoint of its two caps is
    the plotted value whatever the marker's shape (only bars whose caps are unambiguous are used: a cap must
    END its bar, and no other marker of the same date may sit between the marker and the cap);
(3) writes overlay PNGs - original on the left, the same crop with every CSV value's full NOMINAL marker
    outline and reading range drawn on the right. A correct value makes the outline coincide with the
    visible marker, including markers that are partly hidden.
Needs Pillow, numpy, scipy (+ pymupdf when given the PDF).
"""
import argparse
import csv
import datetime
import os

import numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "..", "raw_si", "gao2021_fig3_digitised.csv")

# calibration used by the CSV: x ticks on the 20th of Apr..Sep; y(0) and y(top) per panel row
XT = {2018: (171.5, 324.0, 483.0, 634.0, 792.5, 950.5), 2019: (1040.0, 1190.5, 1349.5, 1502.5, 1660.0, 1819.0)}
YB = {"bio": (659.5, 105.0, 16000.0), "LAI": (1249.0, 694.0, 6.0)}
PANELS = {"a": ((150, 200), (90, 680), "bio"), "c": ((1015, 1065), (90, 680), "bio"),
          "b": ((150, 200), (680, 1270), "LAI"), "d": ((1015, 1065), (680, 1270), "LAI")}
SHAPE = {"D-AGB": "square", "LAI": "square", "DWSO": "circle", "DWST": "tri_up", "DWGL": "tri_down"}
INK = {"D-AGB": (0, 90, 0), "LAI": (0, 90, 0), "DWSO": (150, 0, 0), "DWST": (170, 80, 0), "DWGL": (0, 0, 140)}
FLAGCOL = {"clear": (255, 0, 255), "partial": (0, 170, 255), "near_axis": (0, 150, 0), "on_axis_zero": (110, 110, 110)}


def load_image(path):
    if path.lower().endswith(".pdf"):
        import pymupdf
        doc = pymupdf.open(path)
        pix = pymupdf.Pixmap(doc, doc[7].get_images(full=True)[0][0])
        return Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
    return Image.open(path).convert("RGB")


def rows():
    with open(CSV) as f:
        return list(csv.DictReader(line for line in f if not line.startswith("#")))


def to_px(year, date, series, value):
    doys = [datetime.date(year, m, 20).timetuple().tm_yday for m in (4, 5, 6, 7, 8, 9)]
    x = float(np.interp(datetime.date.fromisoformat(date).timetuple().tm_yday, doys, XT[year]))
    y0, yt, vt = YB["LAI" if series == "LAI" else "bio"]
    return x, y0 - float(value) / vt * (y0 - yt), (y0 - yt) / vt


def _runs(v):
    out, i = [], 0
    while i < len(v):
        if v[i]:
            j = i
            while j < len(v) and v[j]:
                j += 1
            out.append((i, j - 1)); i = j
        else:
            i += 1
    return out


def tick_fit(im):
    dark = im.max(axis=2) < 120
    print("(1) axis calibration re-fitted on the tick marks")
    for name, (xr, yr, kind) in PANELS.items():
        ax = max(range(*xr), key=lambda x: max((b - a for a, b in _runs(dark[yr[0]:yr[1], x])), default=0))
        col = dark[yr[0]:yr[1], ax - 6:ax - 1].sum(axis=1)
        ticks = np.array([(a + b) / 2 + yr[0] for a, b in _runs(col >= 3)])
        step = 2000 if kind == "bio" else 1
        vals = np.arange(len(ticks))[::-1] * step
        p = np.polyfit(ticks, vals, 1)
        y0, yt, vt = YB[kind]
        print(f"   panel {name}: {len(ticks)} ticks  y(0) {-p[1] / p[0]:.2f} (used {y0})  "
              f"y({vt:g}) {(vt - p[1]) / p[0]:.2f} (used {yt})  max tick residual {abs(vals - np.polyval(p, ticks)).max() / abs(p[0]):.2f} px")


def cap_check(im, data):
    R, G, B = im[..., 0], im[..., 1], im[..., 2]
    mx = im.max(axis=2)
    dark = (mx < 140) & (np.ptp(im, axis=2) < 40)
    masks = {"D-AGB": (R < 60) & (G > 90) & (G < 170) & (B < 60), "DWSO": (R > 200) & (G < 90) & (B < 90),
             "DWST": (R > 200) & (G > 100) & (G < 190) & (B < 90), "DWGL": (R < 90) & (G < 90) & (B > 200)}
    masks["LAI"] = masks["D-AGB"]
    anycol = masks["D-AGB"] | masks["DWSO"] | masks["DWST"] | masks["DWGL"]

    def blob(s, px, py):
        x0, y0 = int(px) - 25, int(py) - 25
        lab, n = ndimage.label(masks[s][y0:y0 + 51, x0:x0 + 51]); best = None
        for i in range(1, n + 1):
            ys, xs = np.where(lab == i)
            if len(xs) < 8:
                continue
            cx, cy = (xs.min() + xs.max()) / 2 + x0, (ys.min() + ys.max()) / 2 + y0
            d = np.hypot(cx - px, cy - py)
            if best is None or d < best[0]:
                best = (d, cx, cy, ys.min() + y0, ys.max() + y0)
        return best

    def caps(bx, ylo, yhi):
        hit = []
        for y in range(ylo, yhi):
            for b in (bx - 1, bx, bx + 1):
                if not dark[y, b]:
                    continue
                l, r = b, b
                while dark[y, l - 1] and b - l < 20:
                    l -= 1
                while dark[y, r + 1] and r - b < 20:
                    r += 1
                if 6 <= r - l + 1 <= 16 and abs((l + r) / 2 - bx) <= 3:
                    hit.append(y); break
        out = []
        for a, b in _runs(np.isin(np.arange(ylo, yhi), hit)):
            grp = np.arange(ylo + a, ylo + b + 1)
            w = np.array([255 - mx[y, bx - 4:bx + 5].mean() for y in grp])
            out.append(float(np.dot(grp, w) / w.sum()))
        return out

    def ends_bar(c, d, bx):           # d=-1: cap above the marker
        yc = int(round(c))
        return (any(dark[yc - d * k, bx - 1:bx + 2].any() for k in (2, 3, 4))
                and not any(dark[yc + d * k, bx - 1:bx + 2].any() for k in (3, 4, 5)))

    def continuous(bx, ya, yb):
        gap = 0
        for y in range(int(min(ya, yb)) + 1, int(max(ya, yb))):
            gap = 0 if (dark[y, bx - 1:bx + 2].any() or anycol[y, bx - 2:bx + 3].any()) else gap + 1
            if gap > 3:
                return False
        return True

    pts = []
    for r in data:
        if r["flag"] != "clear":
            continue
        s, yr = r["series"], int(r["year"])
        px, py, sc = to_px(yr, r["date"], s, r["value"])
        b = blob(s, px, py)
        if b:
            pts.append(dict(r=r, s=s, yr=yr, sc=sc, cx=b[1], cy=b[2], top=b[3], bot=b[4]))
    off = {}
    for p in pts:
        same = [q for q in pts if q is not p and q["yr"] == p["yr"] and q["r"]["date"] == p["r"]["date"]
                and (q["s"] == "LAI") == (p["s"] == "LAI")]
        cand = list(range(int(round(p["cx"])) - 2, int(round(p["cx"])) + 3))
        bx = max(cand, key=lambda c: dark[max(0, p["top"] - 60):p["top"] - 1, c].sum() + dark[p["bot"] + 2:p["bot"] + 60, c].sum())
        cs = caps(bx, max(1, p["top"] - 300), min(im.shape[0] - 1, p["bot"] + 300))
        up = [c for c in cs if c < p["top"] - 0.5 and ends_bar(c, -1, bx)]
        dn = [c for c in cs if c > p["bot"] + 0.5 and ends_bar(c, +1, bx)]
        if not up or not dn:
            continue
        up, dn = max(up), min(dn)
        if any(min(c, p["cy"]) < q["cy"] < max(c, p["cy"]) for q in same for c in (up, dn)):
            continue
        if continuous(bx, up, p["top"]) and continuous(bx, p["bot"], dn):
            off.setdefault(SHAPE[p["s"]], []).append((up + dn) / 2 - p["cy"])
    allv = np.concatenate([np.array(v) for v in off.values()])
    print("(2) error-bar cap midpoint minus marker centre (px; + = midpoint lower, value reads high)")
    for k, v in off.items():
        print(f"   {k:9s} n={len(v):2d}  mean {np.mean(v):+.2f}  max|.| {np.max(np.abs(v)):.1f}")
    print(f"   all       n={len(allv):2d}  mean {allv.mean():+.2f} +- {allv.std(ddof=1):.2f}  max|.| {np.abs(allv).max():.1f}")


def outline(d, s, X, Y, S, x0, y0):
    T = lambda x, y: ((x - x0 + 0.5) * S, (y - y0 + 0.5) * S)
    if s in ("D-AGB", "LAI"):
        h = 6.5 if s == "D-AGB" else 7.0
        d.polygon([T(X - h, Y - h), T(X + h, Y - h), T(X + h, Y + h), T(X - h, Y + h)], outline=INK[s])
    elif s == "DWST":
        d.polygon([T(X, Y - 8), T(X + 9.5, Y + 8), T(X - 9.5, Y + 8)], outline=INK[s])
    elif s == "DWGL":
        d.polygon([T(X - 9, Y - 7.5), T(X + 9, Y - 7.5), T(X, Y + 7.5)], outline=INK[s])
    else:
        d.ellipse([T(X - 8, Y - 8), T(X + 8, Y + 8)], outline=INK[s])


def overlays(img, data, out):
    crops = {"panel_a_2018": ((165, 100, 1010, 665), 1), "panel_c_2019": ((1035, 100, 1878, 665), 1),
             "panel_b_2018_LAI": ((165, 690, 1010, 1255), 1), "panel_d_2019_LAI": ((1035, 690, 1878, 1255), 1),
             "zoom_2018_early": ((360, 626, 500, 664), 6), "zoom_2019_early": ((1215, 610, 1400, 664), 5)}
    os.makedirs(out, exist_ok=True)
    for name, ((x0, y0, x1, y1), S) in crops.items():
        c = img.crop((x0, y0, x1, y1)).resize(((x1 - x0) * S, (y1 - y0) * S), Image.NEAREST)
        orig = c.copy(); d = ImageDraw.Draw(c)
        for r in data:
            s, yr = r["series"], int(r["year"])
            X, Y, _ = to_px(yr, r["date"], s, r["value"])
            if not (x0 - 10 <= X < x1 + 10 and y0 - 10 <= Y < y1 + 10):
                continue
            outline(d, s, X, Y, S, x0, y0)
            _, ylo, _ = to_px(yr, r["date"], s, r["value_lo"]); _, yhi, _ = to_px(yr, r["date"], s, r["value_hi"])
            xs = (X - x0 + 0.5) * S
            d.line([xs, (yhi - y0 + 0.5) * S, xs, (ylo - y0 + 0.5) * S], fill=FLAGCOL[r["flag"]], width=1)
        w = Image.new("RGB", (c.width * 2 + 10, c.height), (255, 255, 255))
        w.paste(orig, (0, 0)); w.paste(c, (c.width + 10, 0))
        w.save(os.path.join(out, f"gao2021_fig3_overlay_{name}.png"))
    print(f"(3) overlays written to {out}/gao2021_fig3_overlay_*.png "
          "(left: original; right: nominal outline + reading range of every CSV value)")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("source", help="the article PDF, or the Fig. 3 bitmap as PNG")
    ap.add_argument("--out", default=".", help="directory for the overlay PNGs")
    a = ap.parse_args()
    img = load_image(a.source)
    assert img.size == (1900, 1367), f"unexpected bitmap size {img.size}"
    im = np.asarray(img).astype(int)
    data = rows()
    tick_fit(im)
    cap_check(im, data)
    overlays(img, data, a.out)


if __name__ == "__main__":
    main()
