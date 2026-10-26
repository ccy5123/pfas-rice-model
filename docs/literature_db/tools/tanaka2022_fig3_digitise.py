"""Digitise Tanaka et al. 2022, Rice 15:10 (10.1186/s12284-022-00557-5; CC BY 4.0), Fig. 3 A and C:
the wild-type (cv. Notohikari, plant height 115.8 cm; Kawatabi, Miyagi, 2021) canopy cut in 20-cm
layers from the soil surface at 10 and 49 days after heading -- dry weight [g m-2] of leaf blades
(green, left axis), sheaths + stems (dark grey), panicles (orange) and dead organs (white) per layer.

Usage: extract the figure page image from the open-access PDF (pdfimages -j, page 5), then
    python docs/literature_db/tools/tanaka2022_fig3_digitise.py p5-000.jpg
Calibration: the axis tick marks of the 1595 x 1670 px page image (left axis 80 g m-2 at x 157.0,
0 at 497.75; right axis 600 g m-2 at x 834.0, 0 at 498.0; height ticks per panel below). Each
layer is read on 18 rows away from its error bars; the median edge is kept. Output:
docs/literature_db/raw_si/tanaka2022_fig3_layers.csv (an independent reading by a literature
subagent agreed within 0.01 on every straw share).
"""
import sys
import numpy as np
from PIL import Image
X0 = 497.75                                           # the zero line (left ticks 497.5, right 498)
PX_LEFT = (497.75 - 157.0) / 80.0                     # px per g m-2, leaf-blade axis
PX_RIGHT = (834.0 - 498.0) / 600.0                    # px per g m-2, other organs
YT = {"A": [659.5, 554.0, 449.0, 344.5, 240.5, 135.5, 30.0],      # 0,20,...,120 cm
      "C": [1516.7, 1411.3, 1306.0, 1201.5, 1097.5, 992.5, 887.0]}

def cls(p):
    r, g, b = p
    if r > 200 and g > 200 and b > 200: return "w"
    if r < 90 and g < 90 and b < 90: return "k"
    if r > 200 and 80 < g < 170 and b < 80: return "o"
    if g > 170 and g - r > 20 and g - b > 20: return "g"
    return "?"

def runs(row, x0):
    out, cur, s = [], row[0], 0
    for i, c in enumerate(row[1:], 1):
        if c != cur:
            out.append([cur, s + x0, i - 1 + x0]); cur, s = c, i
    out.append([cur, s + x0, len(row) - 1 + x0])
    return out

def right_side(y):
    row = [cls(im[y, x]) for x in range(499, 852)]
    rr = runs(row, 499)
    # the background starts at the first white run longer than 20 px (text further right,
    # such as the panel labels, is not part of the bar)
    ib = next(i for i, r in enumerate(rr) if r[0] == "w" and r[2] - r[1] > 20)
    bg = rr[ib]
    rr = rr[:ib + 1]
    if ib == 0:
        return X0, None, X0
    total = (rr[-2][1] + rr[-2][2]) / 2 if rr[-2][0] == "?" else bg[1] - 0.5
    kend = oend = None
    for i, (c, a, b) in enumerate(rr[:-1]):
        nxt = rr[i + 1]
        edge = (nxt[1] + nxt[2]) / 2 if nxt[0] == "?" else b + 0.5
        if c == "k" and b - a > 1: kend = edge
        if c == "o" and b - a > 1: oend = edge
    kend = X0 if kend is None else kend
    return kend, oend, total

def left_side(y):
    row = [cls(im[y, x]) for x in range(140, 497)]
    rr = runs(row, 140)
    g = [r for r in rr if r[0] == "g" and r[2] - r[1] > 1]
    if not g:
        return X0
    i = rr.index(g[0])
    prev = rr[i - 1]
    return (prev[1] + prev[2]) / 2 if prev[0] == "?" else g[0][1] - 0.5

def main(path):
    global im
    im = np.asarray(Image.open(path).convert("RGB")).astype(int)
    out = {}
    for panel in ("A", "C"):
        yt = YT[panel]
        for band in range(6):
            ylo, yhi = yt[band], yt[band + 1]          # ylo is the lower edge (larger pixel y)
            h = ylo - yhi
            rows = [int(round(yhi + f * h)) for f in np.r_[np.linspace(0.18, 0.38, 9), np.linspace(0.62, 0.82, 9)]]
            L, K, O, T = [], [], [], []
            for y in rows:
                L.append(left_side(y))
                k, o, t = right_side(y)
                K.append(k); O.append(o if o is not None else k); T.append(t)
            md = lambda v: float(np.median(v))
            blade = (X0 - md(L)) / PX_LEFT
            stem = (md(K) - X0) / PX_RIGHT
            pan = (md(O) - md(K)) / PX_RIGHT
            dead = (md(T) - md(O)) / PX_RIGHT
            spread = max(np.ptp(L) / PX_LEFT, np.ptp(T) / PX_RIGHT)
            out[(panel, band)] = (blade, stem, pan, dead)
            print(f"{panel} {20*band:>3}-{20*band+20:<3} blade {blade:6.1f}  stem {stem:6.1f}  panicle {pan:6.1f}"
                  f"  dead {dead:5.1f}   (row spread {spread:.1f})")
    for panel in ("A", "C"):
        straw = [sum(out[(panel, b)][i] for i in (0, 1, 3)) for b in range(6)]
        bands = [straw[0], straw[1], straw[2], sum(straw[3:])]
        tot = sum(bands)
        print(panel, "straw per band (0-20, 20-40, 40-60, >60):", [round(v, 1) for v in bands],
              "shares", [round(v / tot, 4) for v in bands], "total", round(tot, 1))
    return out


def write_csv(out, path):
    import csv
    with open(path, "w", newline="") as f:
        f.write("# Tanaka et al. 2022, Rice 15:10 (10.1186/s12284-022-00557-5; CC BY 4.0), Fig. 3 A (10 DAH) and C (49 DAH),\n"
                "# wild type (cv. Notohikari), 15 g N m-2 plots, mean of three plots of four hills. Layers cut every 20 cm\n"
                "# from the soil surface (\"clipped every 20 cm vertically\", Methods). Dry weight [g m-2], digitised by\n"
                "# docs/literature_db/tools/tanaka2022_fig3_digitise.py (reading about +-0.5 g m-2 blades, +-2 g m-2 others).\n")
        w = csv.writer(f)
        w.writerow(["dah", "layer_cm", "leaf_blade", "sheath_stem", "panicle", "dead"])
        for panel, dah in (("A", 10), ("C", 49)):
            for b in range(6):
                v = out[(panel, b)]
                w.writerow([dah, f"{20 * b}-{20 * b + 20}"] + [f"{max(x, 0.0):.1f}" for x in v])


if __name__ == "__main__":
    o = main(sys.argv[1])
    if len(sys.argv) > 2:
        write_csv(o, sys.argv[2])
