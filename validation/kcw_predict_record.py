#!/usr/bin/env python3
"""K_cw side project, step 4 (docs/kcw_side_project/PREREG_KR.md §3, §5; main PREREG 11.2): a RECORD, no decision.

Predicted root cell-wall share of the tissue burden, (1-theta) f_cw K_cw / B_root (basis A; rice root composition at
its midpoints, the range over every corner of the K_cw ends, the ether offset ends and the composition ranges), against
the subcellular fractions of Guo 2025 (O. rosea root), Liu 2023 (wheat root, no humic acid), Yu 2021 (lettuce root
PFOS) and Wang 2020 (Alisma root, 16 d). K_PL: Chen 2025 measured where listed, else the literature_params QSPR (ether
offset included); K_prot: the existing path (curated for a known congener, else the QSPR); K_cw: src/kcw_sourced.py.
Records the direction pairs X1-X5 and Spearman rho per dataset. §3.3: the SMILES model impact (root / straw / brown-rice
dry BAF) for five new structures and PFOA/PFOS controls at K_cw = 20 (old placeholder), the rule (rec, low, high) and,
where Mel 2024 measured the compound, the Kraft-lignin value.
Writes validation/kcw_predict_record.json.
"""
import csv
import io
import itertools
import json
import math
import os
import sys

import numpy as np
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
import kcw_sourced as KS                  # noqa: E402
import literature_params as L             # noqa: E402
import model_api as api                   # noqa: E402
from pfas_structure import compound_from_smiles, descriptors   # noqa: E402

RAW = os.path.join(ROOT, "docs", "literature_db", "raw_si")
OUT = os.path.join(HERE, "kcw_predict_record.json")

CF2 = "C(F)(F)"
SMILES = {
    "PFBA": "OC(=O)" + CF2 * 2 + "C(F)(F)F",
    "PFOA": "OC(=O)" + CF2 * 6 + "C(F)(F)F",
    "PFNA": "OC(=O)" + CF2 * 7 + "C(F)(F)F",
    "PFHxS": "OS(=O)(=O)" + CF2 * 5 + "C(F)(F)F",
    "PFOS": "OS(=O)(=O)" + CF2 * 7 + "C(F)(F)F",
    "HFPO-DA": "OC(=O)C(F)(OC(F)(F)C(F)(F)C(F)(F)F)C(F)(F)F",
    "HFPO-TrA": "OC(=O)C(F)(C(F)(F)F)OC(F)(F)C(F)(C(F)(F)F)OC(F)(F)C(F)(F)C(F)(F)F",
    "HFPO-TeA": ("OC(=O)C(F)(C(F)(F)F)OC(F)(F)C(F)(C(F)(F)F)OC(F)(F)C(F)(C(F)(F)F)"
                 "OC(F)(F)C(F)(F)C(F)(F)F"),
    "6:2 Cl-PFESA": "OS(=O)(=O)" + CF2 + CF2 + "O" + CF2 * 5 + "C(F)(F)Cl",
    "wH-PFBA": "OC(=O)" + CF2 * 2 + "C(F)F",
    "wH-PFNA": "OC(=O)" + CF2 * 7 + "C(F)F",
}
CHEN_ALIAS = {"6:2 Cl-PFESA": "9Cl-PF3ONS"}          # Chen 2025's name for 6:2 Cl-PFESA
NEW_STRUCTURES = {                                   # §3.3
    "ADONA": "OC(=O)C(F)OC(F)(F)C(F)(F)C(F)(F)OC(F)(F)C(F)(F)C(F)(F)F",   # 4,8-dioxa-3H-PFNA
    "PFNS (C9 PFSA)": "OS(=O)(=O)" + CF2 * 8 + "C(F)(F)F",
    "br-PFOA (5-CF3)": "OC(=O)" + CF2 * 4 + "C(F)(C(F)(F)F)C(F)(F)F",
    "HFPO-TrA": SMILES["HFPO-TrA"],
    "6:2 Cl-PFESA": SMILES["6:2 Cl-PFESA"],
    "PFOA (control)": SMILES["PFOA"],
    "PFOS (control)": SMILES["PFOS"],
}
MEL_KRAFT = {"PFBA": 10 ** 0.88, "PFOA": 10 ** 1.22, "PFBS": 10 ** 1.43, "PFOS": 10 ** 1.54}


def _rows(name):
    with open(os.path.join(RAW, name), encoding="utf-8") as f:
        return list(csv.DictReader(io.StringIO("".join(l for l in f if not l.startswith("#")))))


def chen_logkmw():
    return {r["pfas"]: float(r["logKMW_pH7.0"]) for r in _rows("chen2025_kmw_hsa.csv")}


def binding_inputs(name):
    """{n, head, n_ether_O, branched, K_PL, K_PL_src, K_prot, kcw (KcwSourced)}."""
    d = descriptors(SMILES[name])
    cmpd, _ = compound_from_smiles(SMILES[name])
    chen = chen_logkmw()
    key = CHEN_ALIAS.get(name, name)
    if key in chen:
        K_PL, src = 10.0 ** chen[key], f"Chen 2025 ({key})"
    elif name.startswith("HFPO") or name.startswith("wH") or d.matched_name is None:
        hg = d.head_group if d.head_group in L.KPL_PER_CF2 else "carboxylate"
        K_PL, src = L.k_pl(d.n_perfluoroC, hg, n_ether_O=d.n_ether_O), "QSPR (literature_params.k_pl)"
    else:
        K_PL, src = cmpd.K_PL, "curated"
    k = KS.kcw_sourced(d.n_perfluoroC, d.head_group, n_ether_O=d.n_ether_O, branched=not d.is_linear)
    return dict(n=d.n_perfluoroC, head=d.head_group, n_ether_O=d.n_ether_O, branched=not d.is_linear,
                K_PL=K_PL, K_PL_src=src, K_prot=cmpd.K_prot, kcw=k)


def root_composition():
    mid = dict(api.PARAMS["tissue_composition_recommended"]["root"])
    rng = {q: tuple(v) for q, v in api.PARAMS["tissue_composition_ranges"]["root"].items()
           if q in ("theta_fw", "f_prot", "f_PL", "f_cw")}
    return mid, rng


def cw_share(comp, K_prot, K_PL, K_cw):
    th = comp["theta_fw"]
    cw = (1 - th) * comp["f_cw"] * K_cw
    B = th + (1 - th) * (comp["f_prot"] * K_prot + comp["f_PL"] * K_PL + comp["f_cw"] * K_cw)
    return cw / B


def predict(name):
    b = binding_inputs(name)
    mid, rng = root_composition()
    k = b["kcw"]
    rec = cw_share(mid, b["K_prot"], b["K_PL"], k.rec)
    vals = []
    for kcw in (k.low, k.rec, k.high):
        for corner in itertools.product(*[rng[q] for q in rng]):
            comp = dict(mid, **dict(zip(rng, corner)))
            vals.append(cw_share(comp, b["K_prot"], b["K_PL"], kcw))
    return dict(rec=rec, lo=min(vals), hi=max(vals), K_cw=k.rec, K_cw_range=[k.low, k.high],
                K_cw_case=k.case, K_cw_flags=k.flags, K_PL=b["K_PL"], K_PL_src=b["K_PL_src"],
                K_prot=b["K_prot"], n=b["n"], head=b["head"], n_ether_O=b["n_ether_O"],
                branched=b["branched"])


def observations():
    guo = {r["pfas"]: float(r["f_cell_wall"]) for r in _rows("guo2025_fig3b_subcellular.csv")}
    liu = {r["pfas"]: float(r["f_F1_cell_wall"]) for r in _rows("liu2023_fig3a_digitised.csv")
           if r["treatment"] == "PFAS"}
    yu = [float(r["f_F1_cell_wall"]) for r in _rows("yu2021_fig1b_digitised.csv")]
    # Wang 2020 Table S11, root at 16 d (cell walls %): PFOA 29.1, PFOS 47.7
    with open(os.path.join(RAW, "wang2020_tables11_proportion.csv"), encoding="utf-8") as f:
        lines = [l for l in csv.reader(f)]
    i = next(j for j, l in enumerate(lines) if l and l[0] == "Root")
    row16 = next(l for l in lines[i + 1:] if l and l[0] == "16")
    pct = lambda s: float(s.split("±")[0]) / 100.0          # noqa: E731
    wang = {"PFOA": pct(row16[3]), "PFOS": pct(row16[7])}
    return {"Guo2025": guo, "Liu2023": liu, "Yu2021_PFOS": {"PFOS": float(np.mean(yu))},
            "Wang2020_16d": wang}


def pair(pred, obs, a, b, expected_model):
    """Record the model and observed order of a vs b (a < b ?)."""
    m = pred[a]["rec"] < pred[b]["rec"]
    o = obs[a] < obs[b]
    return dict(pair=f"{a} vs {b}", model=f"{a} {'<' if m else '>='} {b}",
                observed=f"{a} {'<' if o else '>='} {b}", agree=bool(m == o),
                model_as_expected=bool(m == expected_model))


def record_pairs(pred, obs):
    g, l, w = obs["Guo2025"], obs["Liu2023"], obs["Wang2020_16d"]
    out = {}
    out["X1"] = pair(pred, g, "PFOA", "PFBA", True)                     # expected model PFOA < PFBA
    out["X2"] = dict(Guo=pair(pred, g, "PFOS", "PFOA", True), Liu=pair(pred, l, "PFOS", "PFOA", True),
                     Wang=pair(pred, w, "PFOS", "PFOA", True))
    order_m = sorted(l, key=lambda c: pred[c]["rec"])
    order_o = sorted(l, key=lambda c: l[c])
    rho = spearmanr([pred[c]["rec"] for c in l], [l[c] for c in l]).statistic
    out["X3"] = dict(model_order=order_m, observed_order=order_o, spearman_rho=float(rho),
                     expected="model low Cl-PFESA and PFOS, high PFHxS; rho < 0",
                     model_as_expected=bool(set(order_m[:2]) == {"6:2 Cl-PFESA", "PFOS"}
                                            and order_m[-1] == "PFHxS"),
                     rho_negative=bool(rho < 0))
    m = [pred[c]["rec"] for c in ("HFPO-DA", "HFPO-TrA", "HFPO-TeA")]
    o = [g[c] for c in ("HFPO-DA", "HFPO-TrA", "HFPO-TeA")]
    out["X4"] = dict(model=dict(zip(("HFPO-DA", "HFPO-TrA", "HFPO-TeA"), m)),
                     observed=dict(zip(("HFPO-DA", "HFPO-TrA", "HFPO-TeA"), o)),
                     model_DA_gt_TrA_gt_TeA=bool(m[0] > m[1] > m[2]),
                     observed_DA_gt_TrA_gt_TeA=bool(o[0] > o[1] > o[2]))
    lo, hi = sorted((pred["PFBA"]["rec"], pred["PFOA"]["rec"]))
    olo, ohi = sorted((g["PFBA"], g["PFOA"]))
    out["X5"] = dict(model_between=bool(lo <= pred["HFPO-DA"]["rec"] <= hi),
                     observed_between=bool(olo <= g["HFPO-DA"] <= ohi))
    rhos = {}
    for ds in ("Guo2025", "Liu2023"):
        cs = [c for c in obs[ds] if c in pred]
        rhos[ds] = dict(n=len(cs), rho=float(spearmanr([pred[c]["rec"] for c in cs],
                                                         [obs[ds][c] for c in cs]).statistic))
    out["spearman"] = rhos
    return out


def model_impact():
    """§3.3: dry-basis root / straw / brown-rice BAF of the SMILES structures (site-year default)."""
    out = {}
    for label, smi in NEW_STRUCTURES.items():
        d = descriptors(smi)
        k = KS.kcw_sourced(d.n_perfluoroC, d.head_group, n_ether_O=d.n_ether_O, branched=not d.is_linear)
        runs = {"old_placeholder_20": 20.0, "rule_rec": k.rec, "rule_low": k.low, "rule_high": k.high}
        known = d.matched_name
        if known in MEL_KRAFT:
            runs["kraft_mel2024"] = MEL_KRAFT[known]
        res = {}
        for run, kcw in runs.items():
            cmpd, _ = compound_from_smiles(smi, kcw_anchor=kcw)
            group = ("PFSA" if d.head_group == "sulfonate" else "ether" if d.n_ether_O else "PFCA")
            rec = {"name": f"{label}|{run}", "n_C": d.n_C, "group": group,
                   "K_prot_Lkg": cmpd.K_prot, "K_PL_Lkg": cmpd.K_PL,
                   "K_cw_wholecw_Lkg": {"root": kcw}, "f_xy_recommended": cmpd.f_xy}
            r = api.simulate(rec["name"], record=rec)
            b = r["baf_final_dw"]
            res[run] = dict(K_cw=kcw, root=b["root"], straw=b["straw"],
                            brown_rice=b["grain"])
        out[label] = dict(n=d.n_perfluoroC, head=d.head_group, n_ether_O=d.n_ether_O,
                          matched=known, K_cw_case=k.case, runs=res)
    return out


def main():
    pred = {n: predict(n) for n in SMILES}
    obs = observations()
    rec = dict(protocol="docs/kcw_side_project/PREREG_KR.md §3 + §5 (correction 1); main PREREG 11.2/11.7",
               note="a record: no pass mark; root composition = rice (the other plants' composition is unknown)",
               predictions=pred, observations=obs, pairs=record_pairs(pred, obs),
               model_impact=model_impact())
    with open(OUT, "w") as f:
        json.dump(rec, f, indent=1, default=float)
    for n, p in pred.items():
        print(f"{n:13s} K_cw {p['K_cw']:7.3g} ({p['K_cw_case']:10s}) K_PL {p['K_PL']:9.4g} "
              f"pred cw {p['rec']:.3f} [{p['lo']:.3f}, {p['hi']:.3f}]  obs Guo "
              f"{obs['Guo2025'].get(n, float('nan')):.3f}  Liu {obs['Liu2023'].get(n, float('nan')):.3f}")
    print(json.dumps(rec["pairs"], indent=1))
    for k, v in rec["model_impact"].items():
        print(k, v["K_cw_case"], {r: {t: round(x, 3) for t, x in d.items()} for r, d in v["runs"].items()})
    print(f"[written] {os.path.relpath(OUT, ROOT)}")


if __name__ == "__main__":
    main()
