"""Overlay check of raw_si/ansari2017_fig1_digitised.csv against the original Ansari et al. 2017 Fig. 1.

Fig. 1 (p.135 of Bangladesh Agron. J. 19(2):125-140, doi 10.3329/baj.v19i2.31861, open access) is two
bitmaps on page 11 of the article PDF: (A) hull and (B) brown-rice dry weight per grain [mg] for four
cultivars at heading, 10, 20, 30 d after heading and maturity (B has no heading point). Panel B's y label
is printed "A hull weight (mg)"; the caption says brown rice, and its values (13-35 mg) are brown-rice
weights, so the caption is followed.

    python docs/literature_db/tools/ansari2017_fig1_check.py <ansari2017.pdf> [--out DIR]

(1) re-fits the axis calibration on the tick marks and compares it with the constants the CSV uses;
(2) re-detects every Nipponbare (open triangle) marker standing alone and checks the CSV value against it;
    a marker merged with another cultivar's is checked only against the CSV reading range;
(3) writes an overlay PNG per panel: the original with the CSV value's nominal triangle outline and its
    reading range drawn on it. Needs Pillow, numpy, scipy and pymupdf (or poppler's pdfimages).
"""
import argparse
import csv
import os

import numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "..", "raw_si", "ansari2017_fig1_digitised.csv")

# calibration used by the CSV: pixel row of y = 0 and of the top tick, the top tick's value, and the
# category centres (midpoints between the x-axis boundary ticks), on the PDF's own bitmaps
CAL = {"A": dict(image=0, y0=402, ytop=32, vtop=8.0, x={"heading": 154.0, "10DAH": 255.5, "20DAH": 357.0,
                                                          "30DAH": 458.5, "maturity": 559.0}),
       "B": dict(image=2, y0=404, ytop=32, vtop=50.0, x={"10DAH": 178.5, "20DAH": 302.0, "30DAH": 426.0,
                                                          "maturity": 550.0})}
HALF = 5          # half height/width of the plotted triangle [px]


def load_panels(pdf):
    try:
        import fitz
    except ImportError:                           # poppler's pdfimages writes the same bitmaps
        import subprocess
        import tempfile
        tmp = tempfile.mkdtemp()
        subprocess.run(["pdfimages", "-f", "11", "-l", "11", "-png", pdf, os.path.join(tmp, "f")], check=True)
        return {k: Image.open(os.path.join(tmp, f"f-{c['image']:03d}.png")).convert("L") for k, c in CAL.items()}
    doc = fitz.open(pdf)
    imgs = doc[10].get_images(full=True)          # page 11: four entries (two images + their masks)
    out = {}
    for key, c in CAL.items():
        pix = fitz.Pixmap(doc, imgs[c["image"]][0])
        if pix.n - pix.alpha > 3 or pix.colorspace is None or pix.colorspace.n != 3:
            pix = fitz.Pixmap(fitz.csRGB, pix)
        out[key] = Image.frombytes("RGB", (pix.width, pix.height), pix.samples).convert("L")
    return out


def ticks(a, x_axis_col, y_axis_row):
    """Tick rows on the y axis (short marks to its right) and boundary columns on the x axis."""
    col = np.where(a[:, x_axis_col])[0]
    rows = [y for y in range(col.min(), col.max() + 1) if a[y, x_axis_col + 1:x_axis_col + 5].sum() >= 3]
    xs = [x for x in range(a.shape[1]) if a[y_axis_row - 5:y_axis_row, x].sum() >= 3]
    return rows, xs


def is_triangle(rows):
    """An open triangle: narrow apex rows, a full-width base (the last two rows)."""
    return len(rows) >= 8 and rows[0] <= 3 and rows[-1] >= 9 and rows[-2] >= 9 and rows[1] <= 4


def alone_markers(a, c, xc):
    x = int(round(xc))
    win = a[c["ytop"] - 5:c["y0"] - 2, x - 8:x + 9]
    lab, n = ndimage.label(win)
    found = []
    for k in range(1, n + 1):
        ys, _ = np.where(lab == k)
        if len(ys) < 12 or ys.max() - ys.min() + 1 > 11:        # merged blobs are taller
            continue
        rows = [int((lab[y] == k).sum()) for y in range(ys.min(), ys.max() + 1)]
        if is_triangle(rows):
            yc = 0.5 * (ys.min() + ys.max()) + c["ytop"] - 5
            found.append((c["y0"] - yc) / (c["y0"] - c["ytop"]) * c["vtop"])
    return found


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pdf")
    ap.add_argument("--out", default=".")
    args = ap.parse_args()
    panels = load_panels(args.pdf)
    with open(CSV) as f:
        rows = [r for r in csv.DictReader(line for line in f if not line.startswith("#"))]
    ok = True
    for key, c in CAL.items():
        a = np.array(panels[key]) < 160
        colsum = a.sum(0)
        xaxis = int(np.argmax(colsum >= 0.9 * (c["y0"] - c["ytop"])))   # leftmost: the y axis, not the frame
        tr, _ = ticks(a, xaxis, c["y0"])
        print(f"panel {key}: y ticks at rows {tr} (CSV: 0 at {c['y0']}, {c['vtop']:g} at {c['ytop']})")
        ok &= tr[0] == c["ytop"] and tr[-1] == c["y0"]
        im = panels[key].convert("RGB")
        d = ImageDraw.Draw(im)
        for r in (r for r in rows if r["panel"] == key):
            v, lo, hi = float(r["value_mg"]), float(r["value_lo"]), float(r["value_hi"])
            xc = c["x"][r["time"]]
            ypx = lambda val: c["y0"] - val / c["vtop"] * (c["y0"] - c["ytop"])
            y = ypx(v)
            d.polygon([(xc, y - HALF), (xc - HALF, y + HALF), (xc + HALF, y + HALF)], outline=(255, 0, 255))
            d.line([(xc + 9, ypx(lo)), (xc + 9, ypx(hi))], fill=(0, 150, 255), width=2)
            if r["flag"] == "clear":
                det = alone_markers(a, c, xc)
                hit = [m for m in det if abs(m - v) <= max(0.5 * (hi - lo), 0.02 * c["vtop"])]
                print(f"  {r['time']:9s} CSV {v:6.2f} detected {['%.2f' % m for m in det]} "
                      f"{'OK' if hit else 'MISMATCH'}")
                ok &= bool(hit)
            else:
                print(f"  {r['time']:9s} CSV {v:6.2f} [{lo:.2f}, {hi:.2f}] ({r['flag']}: check the overlay)")
        path = os.path.join(args.out, f"ansari2017_fig1{key}_overlay.png")
        im.save(path)
        print(f"  overlay -> {path}")
    print("ALL CHECKS PASS" if ok else "CHECK FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
