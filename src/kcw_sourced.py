"""Whole-cell-wall K_cw from measured plant-derived sorbents (PREREG §11, amendment 4).

One K_cw [L/kg cell-wall dw] per compound, the same in every organ:

- upper end: Makselon 2019 Sphagnum peat Oi, mean Kd of the runs at pH >= 5.8 (>= 4 runs),
  computed here from ``raw_si/makselon2019_peatOi_tableA9_logkoc.csv`` (Kd = Koc x 0.449);
  log-linear in the perfluorinated-carbon count n between measured chains of one head group;
- lower end: Sörengård 2020 Fig. 1 (pH 7.5), the lowest log Kd of the four plant-derived
  sorbents (peat, granulated bark, composted bark, composted wood fibre);
- representative: the geometric mean of the ends (both ends are input-range band ends).

Outside the domain (PREREG 11.2 decision tree): (i) a Sörengård value but no peat upper end ->
K_cw = that value (``kcw_lower_only``, no band); (ii) neither -> the cell-wall term is removed
(K_cw = 0, ``kcw_removed``), with the nearest measured chain's representative reported.
Ethers: only inside the domain, PFCA/PFSA value of the same n, log offset per ether O 0
(band -0.49 ... +0.49). Any head group other than carboxylate / sulfonate: case (ii).
"""
from __future__ import annotations

import csv
import io
import math
from dataclasses import dataclass, field
from pathlib import Path

RAW = Path(__file__).resolve().parents[1] / "docs" / "literature_db" / "raw_si"
MAKSELON_CSV = RAW / "makselon2019_peatOi_tableA9_logkoc.csv"
SORENGARD_WIDE_CSV = RAW / "sorengard2020_fig1_plant_sorbents.csv"     # PFOA, PFOS (2026-10-27)
SORENGARD_LONG_CSV = RAW / "sorengard2020_fig1_all_compounds.csv"      # every compound (PREREG 11.7 reading)

MAKSELON_OC_FRACTION = 0.449          # Table 1: organic C 44.90 % dw
MAKSELON_MIN_PH = 5.8
MAKSELON_MIN_RUNS = 4
PLANT_SORBENTS = ("Peat", "Granulated bark", "Composted bark", "Composted wood fiber")
ETHER_LOG_OFFSET_RANGE = (-0.49, 0.49)   # per ether O: Chen 2025 K_MW / Muensterman 2026 CHI log D
ETHER_LOG_OFFSET_REC = 0.0

# n = perfluorinated carbons (pfas_structure.descriptors.n_perfluoroC)
PFCA = {3: "PFBA", 4: "PFPeA", 5: "PFHxA", 6: "PFHpA", 7: "PFOA", 8: "PFNA", 9: "PFDA",
        10: "PFUnDA", 11: "PFDoDA", 12: "PFTrDA", 13: "PFTeDA"}
PFSA = {4: "PFBS", 5: "PFPeS", 6: "PFHxS", 7: "PFHpS", 8: "PFOS", 9: "PFNS", 10: "PFDS"}
HEADS = {"carboxylate": PFCA, "sulfonate": PFSA}
# the compounds Sörengård 2020 Fig. 1 plots (side-project PREREG 2.2): asking for one of them
# before its reading is in the repo is an error, not a silent "outside the domain"
SORENGARD_REGISTERED = ("PFBA", "PFPeA", "PFHxA", "PFHpA", "PFOA", "PFNA", "PFDA", "PFUnDA",
                        "PFDoDA", "PFTeDA", "PFBS", "PFHxS", "PFOS")


class PendingReading(RuntimeError):
    """A registered Sörengård compound whose Fig. 1 value has not been read yet."""


def _rows(path):
    return list(csv.DictReader(io.StringIO("".join(l for l in open(path, encoding="utf-8")
                                                   if not l.startswith("#")))))


def makselon_upper() -> dict[str, float]:
    """{compound: mean Kd [L/kg dw] of the peat Oi runs at pH >= 5.8} (>= 4 runs only)."""
    rows = _rows(MAKSELON_CSV)
    out = {}
    for k in rows[0]:
        if not k.startswith("logKoc_"):
            continue
        v = [10.0 ** float(r[k]) * MAKSELON_OC_FRACTION for r in rows
             if r[k] and float(r["pH"]) >= MAKSELON_MIN_PH]
        if len(v) >= MAKSELON_MIN_RUNS:
            out[k[len("logKoc_"):]] = sum(v) / len(v)
    return out


def sorengard_lower() -> dict[str, dict]:
    """{compound: {"Kd", "sorbent", "low_conf", "min_may_be_lower"}} -- the lowest READ log Kd of the four
    plant-derived sorbents (PREREG 11.7); None where the compound is read but has no value (not plotted)."""
    vals: dict[str, list] = {}
    hid: dict[str, list] = {}
    if SORENGARD_LONG_CSV.exists():
        for r in _rows(SORENGARD_LONG_CSV):
            if r["sorbent"] not in PLANT_SORBENTS:
                continue
            vals.setdefault(r["pfas"], [])
            try:
                v = float(r["log_Kd"])
            except ValueError:
                if r["log_Kd"] == "hidden":
                    hid.setdefault(r["pfas"], []).append(float(r["hidden_low"]))
                continue
            vals[r["pfas"]].append((v, r["sorbent"], "low_conf" in (r.get("flag") or "")))
    else:
        for r in _rows(SORENGARD_WIDE_CSV):
            if r["sorbent"] not in PLANT_SORBENTS:
                continue
            for c in ("PFOA", "PFOS"):
                vals.setdefault(c, []).append(
                    (float(r[f"log_Kd_{c}"]), r["sorbent"], f"{c} low_conf" in (r.get("flag") or "")))
    out = {}
    for c, v in vals.items():
        if not v:
            out[c] = None                      # read, but no value (not plotted / not quantifiable)
            continue
        lo = min(v)
        out[c] = {"Kd": 10.0 ** lo[0], "sorbent": lo[1], "low_conf": lo[2],
                  "min_may_be_lower": any(h < lo[0] for h in hid.get(c, []))}
    return out


@dataclass
class KcwSourced:
    rec: float
    low: float
    high: float
    case: str                 # "domain" | "lower_only" | "removed"
    flags: list = field(default_factory=list)
    provenance: str = ""
    nearest_rec: float | None = None     # case "removed": the nearest measured chain's representative


def _upper_at(n, head, upper):
    names = HEADS[head]
    pts = sorted((k, upper[v]) for k, v in names.items() if v in upper)
    if not pts or n < pts[0][0] or n > pts[-1][0]:
        return None
    for (n0, k0), (n1, k1) in zip(pts, pts[1:] + [pts[-1]]):
        if n0 <= n <= n1:
            if n1 == n0:
                return k0
            w = (n - n0) / (n1 - n0)
            return 10.0 ** ((1 - w) * math.log10(k0) + w * math.log10(k1))
    return None


def _lower_at(n, head, lower):
    name = HEADS[head].get(n)
    if name is None:
        return None, None
    if name in SORENGARD_REGISTERED and name not in lower:
        raise PendingReading(f"{name}: Sörengård 2020 Fig. 1 value not read yet "
                             "(docs/literature_db/fulltext_needed.md)")
    return lower.get(name), name


def kcw_sourced(n: int, head: str, n_ether_O: int = 0, branched: bool = False,
                upper=None, lower=None, _nearest=True) -> KcwSourced:
    """K_cw of a PFAS with n perfluorinated carbons and head group `head` (PREREG 11.2)."""
    upper = makselon_upper() if upper is None else upper
    lower = sorengard_lower() if lower is None else lower
    flags = ["branched: n counted, no correction"] if branched else []
    if head not in HEADS:
        return KcwSourced(0.0, 0.0, 0.0, "removed", flags + ["kcw_removed: head group outside"],
                          "PREREG 11.2 (ii)")
    hi = _upper_at(n, head, upper)
    lo_rec, name = _lower_at(n, head, lower)
    lo = None if lo_rec is None else lo_rec["Kd"]
    if hi is not None and lo is not None:
        rec = math.sqrt(hi * lo)
        a, b = sorted((lo, hi))
        if lo_rec["low_conf"]:
            flags.append(f"lower end set by a low_conf reading ({lo_rec['sorbent']})")
        if lo_rec.get("min_may_be_lower"):
            flags.append("min_may_be_lower (a hidden marker could lie below the lowest read value)")
        if n_ether_O:
            e0, e1 = ETHER_LOG_OFFSET_RANGE
            f = 10.0 ** (ETHER_LOG_OFFSET_REC * n_ether_O)
            a, b, rec = a * 10.0 ** (e0 * n_ether_O), b * 10.0 ** (e1 * n_ether_O), rec * f
            flags.append(f"{n_ether_O} ether O: offset 0 (band {e0:+.2f}..{e1:+.2f}/O), linear in O (assumed)")
        interp = HEADS[head].get(n) not in upper
        return KcwSourced(rec, a, b, "domain",
                          flags + (["upper end interpolated (log-linear in n)"] if interp else []),
                          f"PREREG 11.2: geomean(Makselon 2019 peat Oi pH>=5.8 {hi:.3g}, "
                          f"Sörengård 2020 {lo_rec['sorbent']} {lo:.3g})")
    if lo is not None and not n_ether_O:
        if lo_rec["low_conf"]:
            flags.append(f"set by a low_conf reading ({lo_rec['sorbent']})")
        if lo_rec.get("min_may_be_lower"):
            flags.append("min_may_be_lower (a hidden marker could lie below the lowest read value)")
        return KcwSourced(lo, lo, lo, "lower_only", flags + ["kcw_lower_only"],
                          f"PREREG 11.2 (i): Sörengård 2020 {lo_rec['sorbent']} {lo:.3g} (one-sided)")
    near = _nearest_rec(n, head, upper, lower) if _nearest else None
    return KcwSourced(0.0, 0.0, 0.0, "removed", flags + ["kcw_removed"], "PREREG 11.2 (ii)", near)


def _nearest_rec(n, head, upper, lower):
    for m in sorted(HEADS[head], key=lambda m: (abs(m - n), m)):
        try:
            r = kcw_sourced(m, head, upper=upper, lower=lower, _nearest=False)
        except PendingReading:
            continue
        if r.case == "domain":
            return r.rec
    return None
