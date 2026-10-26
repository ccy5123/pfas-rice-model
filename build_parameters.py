#!/usr/bin/env python3
# =============================================================================
# build_parameters.py — assemble the CANONICAL consolidated parameter set
# -----------------------------------------------------------------------------
# Merges the closed deliverables into one machine-readable file for downstream
# (Claude Code) use:
#   K_prot(n), K_PL(n)         <- params/Bk_table_S5.csv     (Zhou2025, Chen2025)
#   K_cw whole-cw per organ    <- params/Kcw_Klignin_params_v2.csv  (GAP A, Prompt2)
#   f_xy_recommended(n)        <- params/Bk_table_S5.csv f_xy (= S4/theory_anchor, MONOTONE)
#   r_eq__<v>, L_Ph__<v>, f_xy__<v> <- params/refit_sourced.csv (validation/refit_sourced.py; Yamazaki 2020
#                                   site-year; one set per variant v = <shoot>_<lipid|free|fit>)
#   lipid_loading[shoot]       <- params/lipid_loading_sourced.json (same script)
#   tissue composition         <- params/rice_composition_sourced.csv (midpoints of sourced ranges)
#   K_cw whole-cw per organ    <- (1 - s) K_poly + s K_lignin, s = the organ's lignin share (same table)
#   B_k(n) per organ (basis-A) <- recomputed here (self-consistent)
# Emits: params/parameters.json, params/f_xy_recommended.csv
# Run:   python build_parameters.py   (from package root)
# =============================================================================
import csv, json, os

HERE = os.path.dirname(os.path.abspath(__file__))
P = os.path.join(HERE, "params")
def rd(name):
    with open(os.path.join(P, name), newline="") as f:
        return list(csv.DictReader(f))

S5  = rd("Bk_table_S5.csv")
KCW = {r["pfas"]: r for r in rd("Kcw_Klignin_params_v2.csv")}
def rd_commented(name):
    with open(os.path.join(P, name), newline="") as f:
        return list(csv.DictReader(line for line in f if not line.startswith("#")))

RC = rd_commented("rice_composition_sourced.csv")

ORGANS = ["root", "stem", "leaf", "grain_brown"]   # grain_brown == model "grain"
QUANT = ("theta_fw", "f_prot", "f_PL", "f_cw", "lignin_share")
ranges = {}
for r in RC:
    ranges.setdefault(r["organ"], {})[r["quantity"]] = (float(r["low"]), float(r["high"]))

def mid(o, q):
    lo, hi = ranges[o][q]
    return round(0.5 * (lo + hi), 6)

# tissue composition = midpoints of the sourced ranges (whole cell wall f_cw)
composition = {o: {q: mid(o, q) for q in ("theta_fw", "f_prot", "f_PL", "f_cw")}
               for o in ORGANS + ["grain_white"]}
lignin_share = {o: mid(o, "lignin_share") for o in ORGANS + ["grain_white"]}

def kcw_whole(k_poly, k_lig, o):
    """Whole-cell-wall K_cw of organ o: lignin-share-weighted mix of the two pools."""
    s = lignin_share[o]
    return round((1.0 - s) * k_poly + s * k_lig, 4)

def Bk_basisA(comp, K_prot, K_PL, K_cw):
    th = comp["theta_fw"]
    return round(th + (1 - th) * (comp["f_prot"] * K_prot
                                  + comp["f_PL"] * K_PL
                                  + comp["f_cw"] * K_cw), 3)

congeners = []
for r in S5:
    name = r["PFAS"]; k = KCW[name]
    K_prot = float(r["K_prot"]); K_PL = float(r["K_PL"])
    Kcw_org = {o: kcw_whole(float(k["K_cw_poly_rec"]), float(k["K_lignin_rec"]), o) for o in ORGANS}
    Bk = {("grain" if o == "grain_brown" else o):
          Bk_basisA(composition[o], K_prot, K_PL, Kcw_org[o]) for o in ORGANS}
    congeners.append({
        "name": name, "n_C": int(r["n_C"]),
        "group": "PFCA" if r["group"] == "CA" else "PFSA",
        "K_prot_Lkg": K_prot, "K_PL_Lkg": K_PL,
        "K_cw_wholecw_Lkg": {("grain" if o == "grain_brown" else o): Kcw_org[o] for o in ORGANS},
        "K_cw_poly_rec_Lkg": float(k["K_cw_poly_rec"]),
        "K_lignin_rec_Lkg": float(k["K_lignin_rec"]),
        "f_xy_recommended": float(r["f_xy"]),                  # MONOTONE (S4/theory) — USE THIS
        "B_k_basisA_Lkg_fw": Bk,
    })

# --- GenX (HFPO-DA), an ETHER-PFAS — PROVISIONAL addition (not in the calibrated 12) ---
# Binding ~ PFPeA: GenX Chen2025 log K_MW = 2.07 ~ PFPeA 2.02, so the membrane K_PL is
# the MEASURED Chen2025 HFPO-DA value and K_prot / K_cw are borrowed from PFPeA (Zhou2025
# has no GenX). Translocation f_xy = carboxylate f_xy at nPFC=5 (PFPeA) × the ether
# head-group offset exp(-0.7) (Tang 2026): GenX is short-chain so it stays MORE shoot-
# mobile than PFOA despite the ether penalty (f_xy 0.23 ≫ PFOA 0.04). Soil Kd uses the
# carboxylate Koc at nPFC=5 (weakly sorbed). Added for ether-PFAS support + the Tang OOS test.
import math
_pe = next(x for x in S5 if x["PFAS"] == "PFPeA")
_pe_kcw = KCW["PFPeA"]
GENX_ETHER_LN_OFFSET = -0.7                                   # literature_params.FXY_HEADGROUP_LN_OFFSET["ether"]
_gx_kprot, _gx_kpl = float(_pe["K_prot"]), 117.5             # K_prot≈PFPeA; K_PL=Chen2025 HFPO-DA
_gx_kcw = {o: kcw_whole(float(_pe_kcw["K_cw_poly_rec"]), float(_pe_kcw["K_lignin_rec"]), o)
           for o in ORGANS}
_gx_fxy = round(float(_pe["f_xy"]) * math.exp(GENX_ETHER_LN_OFFSET), 4)
congeners.append({
    "name": "GenX", "n_C": 5, "group": "ether",
    "K_prot_Lkg": _gx_kprot, "K_PL_Lkg": _gx_kpl,
    "K_cw_wholecw_Lkg": {("grain" if o == "grain_brown" else o): _gx_kcw[o] for o in ORGANS},
    "K_cw_poly_rec_Lkg": float(_pe_kcw["K_cw_poly_rec"]),
    "K_lignin_rec_Lkg": float(_pe_kcw["K_lignin_rec"]),
    "f_xy_recommended": _gx_fxy,
    "B_k_basisA_Lkg_fw": {("grain" if o == "grain_brown" else o):
                          Bk_basisA(composition[o], _gx_kprot, _gx_kpl, _gx_kcw[o]) for o in ORGANS},
    "provisional": "ether-PFAS (HFPO-DA/GenX): K_PL measured (Chen2025), K_prot/K_cw≈PFPeA (matched K_MW), "
                   "f_xy=PFCA(nPFC5)·exp(-0.7) (Tang2026 ether offset); NOT calibrated to BAF data.",
})

# --- merge the transport fit (validation/refit_sourced.py) if present ---
# One set per variant v = <shoot>_<lipid|free|fit>: r_eq__v, L_Ph__v (and f_xy__v for "fit"),
# fitted to Yamazaki 2023 root and grain (fit: also straw) on the 2020 Tsukuba site-year with
# that shoot and loading form; the lipid constants per shoot in lipid_loading_sourced.json.
_refit_path = os.path.join(P, "refit_sourced.csv")
_lipid_path = os.path.join(P, "lipid_loading_sourced.json")
lipid_loading = {}
if os.path.exists(_refit_path):
    _refit = {r["compound"]: r for r in rd_commented("refit_sourced.csv")}
    for cg in congeners:
        rr = _refit.get(cg["name"])
        if rr:
            for k, v in rr.items():
                if k != "compound" and v not in ("", None):
                    cg[k] = round(float(v), 8)
if os.path.exists(_lipid_path):
    with open(_lipid_path) as f:
        lipid_loading = json.load(f)["lipid_loading"]

doc = {
    "_meta": {
        "title": "PFAS–rice 4-compartment uptake model — consolidated parameters",
        "model": "IOC extension of DPU/Trapp; 12 calibrated congeners (PFCA C4–C12 + PFSA C4/C6/C8) + GenX (ether-PFAS, provisional)",
        "binding_basis": "A (fresh-weight): B_k = theta_fw + (1-theta_fw)*sum_i f_i,dw * K_i  [L/kg fw]",
        "units": {"K_*": "L/kg pool-dw", "B_k": "L/kg fw", "theta": "L/kg fw", "f_*": "kg/kg dw"},
        "provenance": {
            "K_PL": "Chen2025 SSLM", "K_prot": "Zhou2025",
            "K_cw": "GAP A / Prompt2 v2 (Guo2025 Fig.3f DFT + Mel2024 lignin anchor; rice whole-cw per organ)",
            "f_xy_recommended": "GAP B / S4 + docs/theory_anchor.tex (Trapp+Briggs LFER: MONOTONE logistic)",
            "transport": "r_eq__<v>, L_Ph__<v> (and f_xy__<v> for v = *_fit): fitted per congener to Yamazaki 2023 root and grain (fit: also straw) on the 2020 Tsukuba site-year with that shoot and loading form; validation/refit_sourced.py, params/refit_sourced.csv (pre-registered: docs/PREREG_sourced_values.md)",
            "lipid_loading": "per shoot model, fitted with the lipid variant's per-congener values (params/lipid_loading_sourced.json); terms the data do not identify are removed",
            "composition": "params/rice_composition_sourced.csv: midpoints of sourced ranges (docs/literature_db/sourced_inputs.csv); cells without a rice source carry physical bounds",
        },
        "WARNING_f_xy": "f_xy_recommended is the monotone a-priori value; f_xy__*_fit is a saturated reproduction fit (one value per congener per dataset) and is not a transferable property.",
    },
    "environment": {
        "z": -1, "f_d": 1.0, "f_n": 0.0, "T_K": 298.15,
        "note": "fully dissociated anion. The root entry is the linear equilibrium term (PREREG amendment 1); the membrane potential and the Michaelis-Menten carrier left the model (no source for their constants).",
    },
    "f_xy_functional_form": {
        "form": "logistic", "expression": "f_xy(n) = eta / (1 + 10^(s*(n - n0)))",
        "logit_natural_log_fit": "logit f_xy = 4.061 - 0.857*n   (OLS, R^2=0.97)",
        "eta_ceiling": 1.0, "beta_per_C_pooled": 0.67, "beta_per_C_OLS": 0.86,
        "s_base10": 0.372, "n0": 4.74,
        "headgroup": "PFSA = PFCA * exp(-1.1)  (Tang2026 + Yamazaki2023 TF: PFOS/PFOA ~0.26-0.43; sign CONFIRMED PFSA<PFCA; was placeholder exp(-1.5))",
        "ceiling_check": "f_xy(PFBA)=0.79 ≈ Briggs neutral max 0.784",
        "source": "docs/theory_anchor.tex (Task A)",
    },
    "tissue_composition_recommended": composition,
    "tissue_composition_ranges": {o: {q: list(v) for q, v in d.items()} for o, d in ranges.items()},
    "congeners": congeners,
    "lipid_loading": lipid_loading,
}

with open(os.path.join(P, "parameters.json"), "w") as f:
    json.dump(doc, f, indent=2, ensure_ascii=False)

# flat f_xy table for quick reference
with open(os.path.join(P, "f_xy_recommended.csv"), "w", newline="") as f:
    wr = csv.writer(f)
    wr.writerow(["pfas", "n_C", "group", "f_xy_recommended_monotone", "f_xy_fit_site_year"])
    for c in congeners:
        wr.writerow([c["name"], c["n_C"], c["group"], c["f_xy_recommended"],
                     c.get("f_xy__redistributed_fit")])

print("[written] params/parameters.json  (%d congeners)" % len(congeners))
print("[written] params/f_xy_recommended.csv")
print("\nf_xy_recommended (monotone)  vs  f_xy__redistributed_fit (site-year reproduction fit):")
for c in congeners:
    ft = (f"{c['f_xy__redistributed_fit']:.4f}" if c.get('f_xy__redistributed_fit') is not None
          else "  ND ")
    print(f"  {c['name']:7}(C{c['n_C']:>2}) {c['group']}  rec={c['f_xy_recommended']:.4f}"
          f"   fit={ft}   B_root={c['B_k_basisA_Lkg_fw']['root']:.2f}")
