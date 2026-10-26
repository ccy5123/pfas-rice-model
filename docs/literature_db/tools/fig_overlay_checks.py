"""Overlay checks of three earlier figure digitisations against the article bitmaps.

    python docs/literature_db/tools/fig_overlay_checks.py kim2003  <Kim2003 GCB PDF>  [--out DIR]
    python docs/literature_db/tools/fig_overlay_checks.py li2025   <Li2025 JHM PDF>   [--out DIR]
    python docs/literature_db/tools/fig_overlay_checks.py tang2026 <Tang2026 JHM PDF> [--out DIR]

The PDFs are not in the repo (copyright). Each check re-fits the axes on the figure's own tick marks,
re-reads the plotted values, prints them next to what the repo stores, and writes an overlay PNG (the
reading drawn on the original) so the placement can be inspected by eye.

kim2003  Fig. 2 (p.7 bitmap): total biomass, open circles = ambient CO2 (stored in
         raw_si/rice_organ_biomass_literature.csv), filled = FACE. The bitmap is an inverted mask.
         Rings are found by an annulus template; a ring hidden under its FACE disk is read from the
         crescent it leaves beyond the disk (ring outer edge - 23 px radius).
li2025   Fig. 2a (p.4 bitmap): irrigation-water composition, stacked bars split at their separator
         lines; the shares are checked against SI Table S5 class sums (raw_si/li2025_stage_water_soil.csv).
tang2026 Fig. 1a/1b (p.4 bitmap): harvest biomass and ear, control bars read at their top border.
Needs pymupdf, Pillow, numpy, scipy.
"""
import argparse
import csv
import os
import re

import numpy as np
from PIL import Image, ImageDraw, ImageOps
from scipy import ndimage
from scipy.signal import fftconvolve

RAW = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "raw_si")


def page_image(pdf, page_index):
    import pymupdf
    doc = pymupdf.open(pdf)
    pix = pymupdf.Pixmap(doc, doc[page_index].get_images(full=True)[0][0])
    if pix.n == 1:
        return Image.frombytes("L", (pix.width, pix.height), pix.samples)
    if pix.n > 3:
        pix = pymupdf.Pixmap(pymupdf.csRGB, pix)
    return Image.frombytes("RGB", (pix.width, pix.height), pix.samples)


def runs(v):
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


def rows_of(name):
    with open(os.path.join(RAW, name)) as f:
        return list(csv.DictReader(line for line in f if not line.startswith("#")))


# ------------------------------------------------------------------------------------------- Kim 2003
KIM_PANEL = {(0, 0): (1998, "HN"), (0, 1): (1999, "HN"), (0, 2): (2000, "HN"), (1, 0): (1998, "MN"),
             (1, 1): (1999, "MN"), (1, 2): (2000, "MN"), (2, 1): (1999, "LN"), (2, 2): (2000, "LN")}


def kim2003(pdf, out):
    img = page_image(pdf, 6).convert("L")
    assert img.size == (3567, 2406), img.size
    ink = np.asarray(img).astype(int) > 128
    cols, rows_ = [(0, 1112), (1224, 2338), (2448, 3567)], [(0, 720), (820, 1560), (1690, 2406)]
    frames, cal = {}, {}
    for ri, (y0, y1) in enumerate(rows_):
        for ci, (x0, x1) in enumerate(cols):
            sub = ink[y0:y1, x0:x1]
            vc = [x for x in range(sub.shape[1]) if max((b - a for a, b in runs(sub[:, x])), default=0) > 0.68 * (y1 - y0)]
            hr = [y for y in range(sub.shape[0]) if max((b - a for a, b in runs(sub[y, :])), default=0) > 0.68 * (x1 - x0)]
            if not vc or not hr:
                continue
            L = np.mean([x for x in vc if x < sub.shape[1] / 2]) + x0
            R = np.mean([x for x in vc if x > sub.shape[1] / 2]) + x0
            T = np.mean([y for y in hr if y < sub.shape[0] / 2]) + y0
            B = np.mean([y for y in hr if y > sub.shape[0] / 2]) + y0
            frames[(ri, ci)] = (L, R, T, B)
    print("(1) per-panel axis calibration: frame = 0..2500 g m-2 and 0..150 d; minor ticks every 100 g m-2")
    for k, (L, R, T, B) in frames.items():
        ys, vs = [T, B], [2500.0, 0.0]
        xin = int(L) + 3
        tk = [y for y in range(int(T) + 3, int(B) - 2) if (lambda rr: rr and rr[0][0] == 0 and 13 <= rr[0][1] + 1 <= 38)(runs(ink[y, xin:xin + 60]))]
        for a, b in runs(np.isin(np.arange(ink.shape[0]), tk)):
            y = (a + b) / 2; v = round(2500 * (B - y) / (B - T) / 100) * 100
            ys.append(y); vs.append(v)
        p = np.polyfit(ys, vs, 1)
        q = np.polyfit([L, R], [0.0, 150.0], 1)
        cal[k] = (p, q)
        print(f"   {KIM_PANEL[k]}: 1 px = {abs(p[0]):.2f} g m-2, {len(ys)} ticks, max tick residual "
              f"{abs(np.array(vs) - np.polyval(p, ys)).max():.1f} g m-2")
    yy, xx = np.mgrid[-24:25, -24:25]; r = np.hypot(yy, xx)
    ann = ((r <= 23.5) & (r >= 17)).astype(float); ann /= ann.sum()
    core = (r <= 12).astype(float); core /= core.sum()
    A = fftconvolve(ink.astype(float), ann, mode="same"); C = fftconvolve(ink.astype(float), core, mode="same")
    lab, n = ndimage.label(A > 0.80)
    found = []
    for i in range(1, n + 1):
        ys_, xs_ = np.where(lab == i); j = np.argmax(A[ys_, xs_]); y, x = ys_[j], xs_[j]
        for k, (L, R, T, B) in frames.items():
            if L < x < R and T < y < B:
                p, q = cal[k]
                found.append((KIM_PANEL[k], "FACE" if C[y, x] > 0.75 else "AMB", float(np.polyval(q, x)), float(np.polyval(p, y)), x, y))
    stored = {}
    for rr in rows_of("rice_organ_biomass_literature.csv"):
        m = re.match(r"Kim2003_GCB_Fig2_(\d{4})_AMB", rr["source"])
        if m:
            stored[(int(m.group(1)), rr["N_level"][:2], int(re.search(r"DAT(\d+)", rr["DAT_or_date"]).group(1)))] = float(rr["value"])
    print("(2) AMB rings found vs stored (g m-2); hidden rings are read from their crescent, see the CSV notes")
    for (yr, nl), kind, dat, val, x, y in sorted(found):
        if kind != "AMB":
            continue
        key = min((k for k in stored if k[0] == yr and k[1] == nl), key=lambda k: abs(k[2] - dat))
        print(f"   {yr} {nl} DAT {dat:5.1f}: ring {val:6.0f}   stored {stored[key]:6.0f} (DAT{key[2]})")
    rgb = ImageOps.invert(img).convert("RGB"); d = ImageDraw.Draw(rgb)
    for (yr, nl), kind, dat, val, x, y in found:
        d.ellipse([x - 23.5, y - 23.5, x + 23.5, y + 23.5], outline=(0, 160, 0) if kind == "AMB" else (220, 0, 0), width=3)
    key_of = {v: k for k, v in KIM_PANEL.items()}
    for (yr, nl, dat), v in stored.items():
        p, q = cal[key_of[(yr, nl)]]
        X, Y = (dat - q[1]) / q[0], (v - p[1]) / p[0]
        d.line([X - 30, Y, X + 30, Y], fill=(0, 0, 255), width=3); d.line([X, Y - 30, X, Y + 30], fill=(0, 0, 255), width=3)
    path = os.path.join(out, "kim2003_fig2_overlay.png"); rgb.save(path)
    print(f"(3) overlay: {path} (green ring = AMB found, red = FACE, blue cross = stored AMB value)")


# ------------------------------------------------------------------------------------------- Li 2025
LI_ORDER = ["L-PFPrS", "PFEtS", "PFBA", "PFHxA", "PFHpA", "PFOA", "PFBS", "PFOS", "6:2 FTS", "PFECHS",
            "8:2 diPAP", "6:2 diPAP", "OBS", "6:2 Cl-PFESA", "8:2 Cl-PFESA"]
LI_STAGES = ["tillering", "heading", "filling", "milk", "ripening"]


def li2025(pdf, out):
    img = page_image(pdf, 3).convert("RGB")
    assert img.size == (1000, 631), img.size
    im = np.asarray(img).astype(int); dark = im.max(axis=2) < 110
    legend_rows = [r for r in runs((im[:, 850].min(axis=1) < 235) & (im[:, 850].max(axis=1) > 60)) if r[1] - r[0] >= 8][:15]
    leg = {nm: im[(a + b) // 2, 845:860].mean(axis=0) for nm, (a, b) in zip(LI_ORDER[::-1], legend_rows)}
    ay = max(range(230, 300), key=lambda y: max((b - a for a, b in runs(dark[y, :820])), default=0))
    tx = np.array([(a + b) / 2 for a, b in runs(dark[ay + 2:ay + 7, :].sum(axis=0) >= 3) if 140 < (a + b) / 2 < 810])
    cal = np.polyfit(tx, np.arange(len(tx)) * 20.0, 1)
    print(f"(1) x axis: {len(tx)} ticks, {1 / cal[0]:.2f} px per %, residual "
          f"{abs(np.arange(len(tx)) * 20.0 - np.polyval(cal, tx)).max():.3f} %")
    bars = [r for r in runs(np.array([np.abs(im[y, 500] - leg["PFOA"]).sum() < 40 for y in range(ay)])) if r[1] - r[0] > 10]
    tot, cls = {}, {}
    for rr in rows_of("li2025_stage_water_soil.csv"):
        if rr["medium"].startswith("irrigation"):
            tot[rr["stage"]] = tot.get(rr["stage"], 0) + float(rr["mean"]); cls.setdefault(rr["stage"], {})[rr["class"]] = float(rr["mean"])
    print("(2) shares (%) per stage, and class sums vs SI Table S5")
    draw = ImageDraw.Draw(img)
    for st, (y0, y1) in zip(LI_STAGES, bars):
        med = np.median(im[y0 + 3:y1 - 2, 140:805], axis=0)
        sep = [i + 140 for i, c in enumerate(med) if c.max() < 170 and c.max() - c.min() < 120]
        seps = [(a + b) / 2 + 140 for a, b in runs(np.isin(np.arange(665), np.array(sep) - 140))]
        shares = {}
        for left, right in zip(seps[:-1], seps[1:]):
            xs = np.arange(int(np.floor(left)) + 1, int(np.ceil(right)))
            if len(xs) == 0:
                continue
            band = im[y0 + 3:y1 - 2, xs]; band = band[:, 1:-1] if band.shape[1] > 3 else band
            c = np.median(band.reshape(-1, 3), axis=0)
            nm = min(leg, key=lambda k: np.abs(leg[k] - c).sum())
            if np.abs(leg[nm] - c).sum() < 40:
                shares[nm] = shares.get(nm, 0) + (right - left) * cal[0]
            for xb in (left, right):
                draw.line([xb, y0 - 3, xb, y1 + 3], fill=(0, 0, 255))
        pfca = sum(shares.get(c, 0) for c in ("PFBA", "PFHxA", "PFHpA", "PFOA"))
        pfsa = sum(shares.get(c, 0) for c in ("PFEtS", "PFBS", "PFOS"))
        print(f"   {st:9s} " + " ".join(f"{k} {v:.2f}" for k, v in shares.items()))
        print(f"   {'':9s} PFCAs {pfca:.2f} vs S5 {100 * cls[st]['PFCAs'] / tot[st]:.2f}   PFSAs {pfsa:.2f} vs S5 {100 * cls[st]['PFSAs'] / tot[st]:.2f}")
    path = os.path.join(out, "li2025_fig2a_overlay.png")
    img.resize((img.width * 3, img.height * 3), Image.NEAREST).save(path)
    print(f"(3) overlay: {path} (blue lines = segment boundaries used)")
    print("   stored per-stage values and 5-stage means: raw_si/li2025_fig2a_water_shares.csv")


# ------------------------------------------------------------------------------------------- Tang 2026
def tang2026(pdf, out):
    img = page_image(pdf, 3).convert("RGB")
    assert img.size == (2098, 1651), img.size
    im = np.asarray(img).astype(int); dark = im.max(axis=2) < 120
    stored = {r["quantity"]: float(r["value_gfw_per_plant"]) for r in rows_of("tang2026_harvest_biomass.csv")}
    d = ImageDraw.Draw(img)
    for name, axx, step, x_end, key in (("1a biomass", 81, 10.0, 690, "whole_plant_biomass_control"),
                                         ("1b ear", 801, 2.0, 1400, "ear_panicle_control")):
        tk = [y for y in range(480) if (lambda rr: rr and rr[0][0] == 0 and 6 <= rr[0][1] + 1 <= 12)(runs(dark[y, axx:axx + 15]))]
        ty = np.array([(a + b) / 2 for a, b in runs(np.isin(np.arange(480), tk))])
        bottom = max(r for r in range(400, 480) if dark[r, axx + 50:axx + 300].mean() > 0.9) - 0.5
        yy = np.r_[ty, bottom]; vv = np.r_[np.arange(len(ty), 0, -1) * step, 0.0]
        cal = np.polyfit(yy, vv, 1)
        ctrl = im[300, axx + 41]                        # inside the first (control) bar
        isb = np.abs(im - ctrl).sum(axis=2) < 60
        bars = [r for r in runs(np.isin(np.arange(im.shape[1]), [x for x in range(axx + 5, x_end) if isb[100:470, x].sum() > 60])) if r[1] - r[0] >= 10]
        vals = []
        for a, b in bars:
            x = a + 3
            fill_top = np.where(isb[60:470, x])[0].min() + 60
            border = [y for y in range(fill_top - 6, fill_top) if dark[y, x]]
            top = np.mean(border) if border else fill_top
            vals.append(float(np.polyval(cal, top)))
            d.line([a - 8, top, b + 8, top], fill=(255, 0, 255), width=2)
        print(f"   Fig {name}: {len(ty)} ticks, residual {abs(vv - np.polyval(cal, yy)).max():.3f}; control bars "
              + ", ".join(f"{v:.2f}" for v in vals) + f"  | stored {stored[key]}")
    path = os.path.join(out, "tang2026_fig1_overlay.png"); img.crop((0, 0, 1400, 520)).save(path)
    print(f"   overlay: {path} (magenta = the top border read)")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("figure", choices=("kim2003", "li2025", "tang2026"))
    ap.add_argument("pdf")
    ap.add_argument("--out", default=".")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    {"kim2003": kim2003, "li2025": li2025, "tang2026": tang2026}[a.figure](a.pdf, a.out)


if __name__ == "__main__":
    main()
