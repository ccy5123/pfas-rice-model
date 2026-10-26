"""Build data_obs/obs_baf_Yamazaki.csv and data_obs/yamazaki_stem_height.csv from the SI
transcription (docs/literature_db/raw_si/yamazaki2023_si_s16_s18_s19.csv).

Every observed number is a tissue concentration (Tables S18/S19, pg/g dw) over the
irrigation-water concentration (Table S16, ng/L = pg/mL) -> BAF [L/kg dw]. The two cultivars
are combined by geometric mean, as the model's two harvest runs are
(model_api.yamazaki_harvest_baf). Where the SI leaves a value open, the observation is the
midpoint of what the SI allows and `lo`/`hi` hold the ends:

  grain  brown rice = white share x white + (1 - white share) x bran, white share 0.86-0.90
         ("weight ratio of white rice to brown rice ranging from 86-90%", SI p.S3; midpoint 0.88);
         a value below the MQL is anywhere in [0, MQL) (midpoint MQL/2).
  straw  "stems, including leaves" in four height bands, whose dry masses Yamazaki does not report.
         The straw is the band concentrations weighted by each band's share of straw dry mass,
         from the stratified cutting of a japonica canopy in 20-cm layers (Tanaka et al. 2022,
         Fig. 3; raw_si/tanaka2022_fig3_layers.csv) at 10 and 49 days after heading -- the two
         profiles bracket Yamazaki's harvest (29-30 d after heading) and are the range's ends;
         the observation is their mean.
  root   one measured value per cultivar.

`use` = 0 for a cell with no measured value at all (PFBS grain: white and bran both < MQL).

Run: python docs/literature_db/tools/yamazaki_obs_from_si.py  (rewrites the two tables).
"""
from __future__ import annotations

import csv
import math
import os

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
SRC = os.path.join(ROOT, "docs", "literature_db", "raw_si", "yamazaki2023_si_s16_s18_s19.csv")
LAYERS = os.path.join(ROOT, "docs", "literature_db", "raw_si", "tanaka2022_fig3_layers.csv")
OBS = os.path.join(ROOT, "data_obs", "obs_baf_Yamazaki.csv")
STEM = os.path.join(ROOT, "data_obs", "yamazaki_stem_height.csv")

CONGENERS = ("PFBA", "PFPeA", "PFHxA", "PFHpA", "PFOA", "PFNA", "PFDA", "PFUnDA", "PFDoDA",
             "PFBS", "PFOS")
BANDS = ("stem_0_20", "stem_20_40", "stem_40_60", "stem_gt60")
WHITE_SHARE = (0.86, 0.90)          # SI p.S3


def _cell(s):
    """(low, mid, high) of one table cell: a number, or '<MQL' -> (0, MQL/2, MQL)."""
    s = s.strip()
    if s.startswith("<"):
        q = float(s[1:])
        return 0.0, 0.5 * q, q
    v = float(s)
    return v, v, v


def _gm(a, b):
    return math.sqrt(a * b)


def read_si(path=SRC):
    out = {}
    with open(path, encoding="utf-8") as f:
        for r in csv.DictReader(line for line in f if not line.startswith("#")):
            d = out.setdefault(r["compound"], {"water": float(r["irrigation_ng_L"]), "cv": {}})
            d["cv"][r["cultivar"]] = {k: _cell(v) for k, v in r.items()
                                      if k not in ("compound", "irrigation_ng_L", "cultivar")}
    return out


def brown(cv):
    """(low, mid, high) brown-rice concentration of one cultivar."""
    w = cv["white_rice"]
    b = cv["bran"]
    mid_s = 0.5 * (WHITE_SHARE[0] + WHITE_SHARE[1])
    ends = [s * w[i] + (1.0 - s) * b[j] for s in WHITE_SHARE for i in (0, 2) for j in (0, 2)]
    return min(ends), mid_s * w[1] + (1.0 - mid_s) * b[1], max(ends)


def straw_shares(path=LAYERS):
    """{dah: (w_0-20, w_20-40, w_40-60, w_>60)}: each band's share of straw dry mass (leaf
    blades + sheaths/stems + dead organs; panicles excluded)."""
    by = {}
    with open(path, encoding="utf-8") as f:
        for r in csv.DictReader(line for line in f if not line.startswith("#")):
            lo = int(r["layer_cm"].split("-")[0])
            band = min(lo // 20, 3)
            straw = float(r["leaf_blade"]) + float(r["sheath_stem"]) + float(r["dead"])
            w = by.setdefault(int(r["dah"]), [0.0] * 4)
            w[band] += straw
    return {d: tuple(x / sum(w) for x in w) for d, w in by.items()}


def straw(cv, shares):
    """(low, mid, high) straw concentration of one cultivar: the two profiles are the ends,
    their mean the observation."""
    v = [sum(wi * cv[k][1] for wi, k in zip(w, BANDS)) for w in shares.values()]
    wm = [sum(w[i] for w in shares.values()) / len(shares) for i in range(4)]
    mid = sum(wi * cv[k][1] for wi, k in zip(wm, BANDS))
    return min(v), mid, max(v)


def build(si=None):
    si = si or read_si()
    shares = straw_shares()
    rows = []
    for c in CONGENERS:
        d = si[c]
        w = d["water"]
        cvs = list(d["cv"].values())
        assert len(cvs) == 2
        root = [cv["root"] for cv in cvs]
        assert all(r[0] == r[2] for r in root), f"{c} root below MQL"
        r_mid = _gm(root[0][1], root[1][1]) / w
        rows.append((c, "root", r_mid, r_mid, r_mid, w, 1))
        sb = [straw(cv, shares) for cv in cvs]
        s_lo, s_mid, s_hi = (_gm(sb[0][i], sb[1][i]) / w for i in range(3))
        rows.append((c, "straw", s_mid, s_lo, s_hi, w, 1))
        g = [brown(cv) for cv in cvs]
        g_lo, g_mid, g_hi = (_gm(g[0][i], g[1][i]) / w for i in range(3))
        measured = any(cv["white_rice"][0] == cv["white_rice"][2]
                       or cv["bran"][0] == cv["bran"][2] for cv in cvs)
        rows.append((c, "grain", g_mid, g_lo, g_hi, w, 1 if measured else 0))
    return rows


def stem_bands(si=None):
    si = si or read_si()
    out = []
    for c in CONGENERS:
        cvs = list(si[c]["cv"].values())
        out.append([c] + [_gm(cvs[0][k][1], cvs[1][k][1]) for k in ("root",) + BANDS])
    return out


HEADER = """# Yamazaki et al. 2023 (10.1021/acs.est.2c08767), NARO Tsukuba Andosol paddy FIELD 2020, Kasalath + Koshihikari.
# BUILT by docs/literature_db/tools/yamazaki_obs_from_si.py from raw_si/yamazaki2023_si_s16_s18_s19.csv (do not edit).
# baf_dw = tissue [pg/g DRY weight] (SI Tables S18 indica, S19 japonica) / irrigation water [pg/mL] (Table S16, one grab
#   2020-06-08) [L/kg dw]; the two cultivars combined by geometric mean (as the model's two harvest runs).
# root  = whole root system (tap + Milli-Q rinse; no plaque removal).
# straw = the four height bands weighted by their share of straw dry mass (Tanaka et al. 2022 Fig. 3, 10 and 49 d after
#         heading; raw_si/tanaka2022_fig3_layers.csv): baf_dw at the mean of the two profiles, lo/hi at each profile.
# grain = brown rice = s x white + (1 - s) x bran, s = 0.86-0.90 (SI p.S3); baf_dw at s = 0.88 and <MQL at MQL/2 (the
#         midpoints), lo/hi over both ends of s and of [0, MQL].
# use=0: no measured value at all (PFBS grain: white and bran < MQL in both cultivars).
"""


def write(rows=None, path=OBS):
    rows = rows or build()
    with open(path, "w", newline="", encoding="utf-8") as f:
        f.write(HEADER)
        wr = csv.writer(f, lineterminator="\n")
        wr.writerow(["compound", "tissue", "baf_dw", "lo", "hi", "water_ng_L", "use"])
        for c, t, mid, lo, hi, w, u in rows:
            wr.writerow([c, t, f"{mid:.6g}", f"{lo:.6g}", f"{hi:.6g}", f"{w:g}", u])
    with open(STEM, "w", newline="", encoding="utf-8") as f:
        f.write("# Yamazaki 2023 (10.1021/acs.est.2c08767) SI Table S18 (Indica) + S19 (Japonica): geometric mean of the two\n"
                "# cultivars, PFAS tissue concentration [pg/g dw], stem(+leaf) bands by height. BUILT by\n"
                "# docs/literature_db/tools/yamazaki_obs_from_si.py (do not edit).\n")
        wr = csv.writer(f, lineterminator="\n")
        wr.writerow(["pfas", "root", "stem_0_20", "stem_20_40", "stem_40_60", "stem_gt60"])
        for r in stem_bands():
            wr.writerow([r[0]] + [f"{v:.1f}" for v in r[1:]])


if __name__ == "__main__":
    write()
    print(f"[written] {os.path.relpath(OBS, ROOT)}, {os.path.relpath(STEM, ROOT)}")
