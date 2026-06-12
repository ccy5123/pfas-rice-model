#!/usr/bin/env python3
# =============================================================================
# build_parameters.py — assemble the CANONICAL consolidated parameter set
# -----------------------------------------------------------------------------
# Merges the closed deliverables into one machine-readable file for downstream
# (Claude Code) use:
#   K_prot(n), K_PL(n)         <- params/Bk_table_S5.csv     (Zhou2025, Chen2025)
#   K_cw whole-cw per organ    <- params/Kcw_Klignin_params_v2.csv  (GAP A, Prompt2)
#   f_xy_recommended(n)        <- params/Bk_table_S5.csv f_xy (= S4/theory_anchor, MONOTONE)
#   f_xy_W2fit(n), L_Ph, kappa_d <- params/W2_transport_fit.csv (S6 transport fit; long-chain NON-physical)
#   tissue composition         <- params/rice_tissue_params.csv (recommended)
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
W2  = {r["congener"]: r for r in rd("W2_transport_fit.csv")}
RT  = rd("rice_tissue_params.csv")

ORGANS = ["root", "stem", "leaf", "grain_brown"]   # grain_brown == model "grain"
def rt(o, p):
    for r in RT:
        if r["organ"] == o and r["parameter"] == p:
            return float(r["value_recommended"])
    return None

# tissue composition (recommended): 4-pool f_cw = polysaccharide + lignin
composition = {}
for o in ORGANS:
    composition[o] = {
        "theta_fw":  rt(o, "theta_fw"),
        "f_prot":    rt(o, "f_prot"),
        "f_PL":      rt(o, "f_PL_membrane"),
        "f_cw":      round(rt(o, "f_cw_polysaccharide") + rt(o, "lignin"), 4),
    }

def Bk_basisA(comp, K_prot, K_PL, K_cw):
    th = comp["theta_fw"]
    return round(th + (1 - th) * (comp["f_prot"] * K_prot
                                  + comp["f_PL"] * K_PL
                                  + comp["f_cw"] * K_cw), 3)

congeners = []
for r in S5:
    name = r["PFAS"]; k = KCW[name]
    K_prot = float(r["K_prot"]); K_PL = float(r["K_PL"])
    Kcw_org = {o: float(k[f"K_cw_wholecw_{o}"]) for o in ORGANS}
    Bk = {("grain" if o == "grain_brown" else o):
          Bk_basisA(composition[o], K_prot, K_PL, Kcw_org[o]) for o in ORGANS}
    w = W2.get(name, {})
    congeners.append({
        "name": name, "n_C": int(r["n_C"]),
        "group": "PFCA" if r["group"] == "CA" else "PFSA",
        "K_prot_Lkg": K_prot, "K_PL_Lkg": K_PL,
        "K_cw_wholecw_Lkg": {("grain" if o == "grain_brown" else o): Kcw_org[o] for o in ORGANS},
        "K_cw_poly_rec_Lkg": float(k["K_cw_poly_rec"]),
        "K_lignin_rec_Lkg": float(k["K_lignin_rec"]),
        "f_xy_recommended": float(r["f_xy"]),                  # MONOTONE (S4/theory) — USE THIS
        "f_xy_W2fit": float(w["f_xy_fit"]) if w else None,     # transport-fit (long-chain artifact)
        "L_Ph_W2fit": float(w["L_Ph_fit"]) if w else None,
        "kappa_d_W2fit": float(w["kappa_d_fit"]) if w else None,
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
_gx_kcw = {o: float(_pe_kcw[f"K_cw_wholecw_{o}"]) for o in ORGANS}
_gx_fxy = round(float(_pe["f_xy"]) * math.exp(GENX_ETHER_LN_OFFSET), 4)
congeners.append({
    "name": "GenX", "n_C": 5, "group": "ether",
    "K_prot_Lkg": _gx_kprot, "K_PL_Lkg": _gx_kpl,
    "K_cw_wholecw_Lkg": {("grain" if o == "grain_brown" else o): _gx_kcw[o] for o in ORGANS},
    "K_cw_poly_rec_Lkg": float(_pe_kcw["K_cw_poly_rec"]),
    "K_lignin_rec_Lkg": float(_pe_kcw["K_lignin_rec"]),
    "f_xy_recommended": _gx_fxy,
    "f_xy_W2fit": None, "L_Ph_W2fit": None, "kappa_d_W2fit": None,
    "B_k_basisA_Lkg_fw": {("grain" if o == "grain_brown" else o):
                          Bk_basisA(composition[o], _gx_kprot, _gx_kpl, _gx_kcw[o]) for o in ORGANS},
    "provisional": "ether-PFAS (HFPO-DA/GenX): K_PL measured (Chen2025), K_prot/K_cw≈PFPeA (matched K_MW), "
                   "f_xy=PFCA(nPFC5)·exp(-0.7) (Tang2026 ether offset); NOT calibrated to BAF data.",
})

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
            "f_xy_W2fit": "S6 W2 transport fit to Yamazaki (SATURATED 3param/3obs; long-chain values NON-physical, see docs/DELIVERABLE_GAP_B_fxy.md)",
            "composition": "rice_tissue package (recommended anchors)",
        },
        "WARNING_f_xy": "Use f_xy_recommended (monotone). f_xy_W2fit rises spuriously for C10+ (single-straw-compartment entanglement, H7 §7.2); theory + cross-field TF both require monotone decline.",
    },
    "environment": {
        "E_m_V": -0.120, "E_m_plausible_range_V": [-0.120, -0.090],
        "z": -1, "f_d": 1.0, "f_n": 0.0, "T_K": 298.15,
        "note": "fully-dissociated anion; E_m is the live pH/redox lever (Task B); in-situ paddy E_m unmeasured (key gap).",
    },
    "carrier_MichaelisMenten": {
        "Vmax_in": 20.0, "Km_in": 5.0, "Vmax_out": 8.0, "Km_out": 5.0,
        "note": "fixed during W2 fit; carrier overcomes GHK anion exclusion (e^N≈107).",
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
    "congeners": congeners,
    "status": {
        "GAP_A_Kcw": "CLOSED (docs/DELIVERABLE_GAP_A_Kcw.md)",
        "GAP_B_fxy": "CLOSED (docs/DELIVERABLE_GAP_B_fxy.md; H7 §7.2 entanglement resolved by theory)",
        "validation": "Yamazaki full-ODE reproduced (log10 RMSE 0.029); cross-field TF (water-independent) confirms monotone f_xy; see docs/H8_handoff_S6_final.md",
        "open_data_limited": [
            "rice (not wheat) per-congener root subcellular -> alpha/QC1 point estimate",
            "reliable per-congener pore-water OR hydroponic RCF -> surface test + f_xy absolute",
            "measured Q_TP(t), M(t) -> f_xy absolute scale",
            "direct K_cw_poly + rice cw monosaccharide composition (long-term weakest point)",
            "in-situ paddy E_m (-90..-120 mV -> ~2.5-3x passive-influx spread)",
        ],
        "config_decisions": {
            "root_theta": "0.90 recommended (rice_tissue/Liu measured 0.90-0.92); 0.70 was an early choice",
            "root_f_PL": "0.015 recommended (estimate, range 0.01-0.02; dominates long-chain B_root)",
            "grain_theta": "0.14 (harvest/dry) vs 0.30 (filling) — stage-dependent",
        },
    },
}

with open(os.path.join(P, "parameters.json"), "w") as f:
    json.dump(doc, f, indent=2, ensure_ascii=False)

# flat f_xy table for quick reference
with open(os.path.join(P, "f_xy_recommended.csv"), "w", newline="") as f:
    wr = csv.writer(f)
    wr.writerow(["pfas", "n_C", "group", "f_xy_recommended_monotone", "f_xy_W2fit_nonphysical_longchain"])
    for c in congeners:
        wr.writerow([c["name"], c["n_C"], c["group"], c["f_xy_recommended"], c["f_xy_W2fit"]])

print("[written] params/parameters.json  (%d congeners)" % len(congeners))
print("[written] params/f_xy_recommended.csv")
print("\nf_xy_recommended (monotone)  vs  f_xy_W2fit (long-chain artifact):")
for c in congeners:
    w2 = f"{c['f_xy_W2fit']:.4f}" if c['f_xy_W2fit'] is not None else "  ND "
    print(f"  {c['name']:7}(C{c['n_C']:>2}) {c['group']}  rec={c['f_xy_recommended']:.4f}"
          f"   W2={w2}   B_root={c['B_k_basisA_Lkg_fw']['root']:.2f}")
