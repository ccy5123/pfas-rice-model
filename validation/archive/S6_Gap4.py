#!/usr/bin/env python3
# =============================================================================
# S6 — Gap4: B_k / TF 독립 검증  (Yamazaki full-chain + Li2025 cross-field)
# -----------------------------------------------------------------------------
# 두 갈래:
#  (1) FULL-ODE 재현: first-principles basis-A B_k(n)[S5 binding + rice_tissue 조성]
#      + W2 transport(f_xy,L_Ph,kappa_d) 로 4-compartment ODE를 돌려 Yamazaki
#      root/straw/grain BAF를 재현 (foundation coherence 검증).
#  (2) CROSS-FIELD TF: TF=tissue/tissue 라 water 분모 상쇄 → Li2025 water artifact
#      면역. Yamazaki vs Li2025 의 TF_straw(n) 단조성(=TSCF 신호) 대조.
#
# 입력: Bk_table_S5.csv, W2_transport_fit.csv, Kcw_Klignin_params_v2.csv,
#       rice_tissue_params.csv, obs_baf_{Yamazaki,Li2025}.csv, Li2025_BAF_TF.csv
# 산출: S6_Gap4_ode_repro.csv, S6_Gap4_TF_crossfield.csv, S6_Gap4.png
# =============================================================================
import csv, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pfas_rice_plant_module_4pool_surf import (
    Environment, Compound, Compartment, RiceUptakeModel, PlantInputs,
    _logistic, ROOT, STEM, LEAF, FRUIT)

def rd(p):
    with open(p, newline="") as f: return list(csv.DictReader(f))

S5  = rd("prompt3_pkg/Prompt3_Kpl_downweight/inputs/Bk_table_S5.csv")
W2  = {r["congener"]: r for r in rd("W2_transport_fit.csv")}
KCW = {r["pfas"]: r for r in rd("Kcw_Klignin_params_v2.csv")}
RTd = rd("rice_tissue_pkg/rice_tissue/rice_tissue_params.csv")
def rt(o, p):
    for r in RTd:
        if r["organ"] == o and r["parameter"] == p: return float(r["value_recommended"])
def obs_baf(path):
    d = {}
    for r in rd(path): d.setdefault(r["compound"], {})[r["tissue"]] = float(r["baf"])
    return d
YA = obs_baf("obs_baf_Yamazaki.csv"); LI = obs_baf("obs_baf_Li2025.csv")
LITF = {r["congener"]: r for r in rd("Li2025_BAF_TF.csv")}

# rice_tissue rec 조성 (= module demo 와 동일)
comps = [
    Compartment("root",  rt("root","theta_fw"),  rt("root","f_prot"),  rt("root","f_PL_membrane"),  rt("root","f_cw_polysaccharide")+rt("root","lignin")),
    Compartment("stem",  rt("stem","theta_fw"),  rt("stem","f_prot"),  rt("stem","f_PL_membrane"),  rt("stem","f_cw_polysaccharide")+rt("stem","lignin")),
    Compartment("leaf",  rt("leaf","theta_fw"),  rt("leaf","f_prot"),  rt("leaf","f_PL_membrane"),  rt("leaf","f_cw_polysaccharide")+rt("leaf","lignin"), S=20.0),
    Compartment("grain", rt("grain_brown","theta_fw"), rt("grain_brown","f_prot"), rt("grain_brown","f_PL_membrane"), rt("grain_brown","f_cw_polysaccharide")+rt("grain_brown","lignin"), S=2.0),
]

# 공통 드라이버 (W2/W1 데모 프로파일; Cwo=1 → C=BAF)
season = 120.0; t = np.linspace(0.0, season, 481)
Cwo = np.full_like(t, 1.0)
Qtp = 0.05 + 0.35*np.exp(-((t-75.0)**2)/(2*25.0**2))
M = np.column_stack([_logistic(t,1e-3,0.030,0.10,20.0), _logistic(t,1e-3,0.040,0.10,25.0),
                     _logistic(t,1e-3,0.050,0.12,30.0), _logistic(t,1e-5,0.025,0.18,80.0)])
inputs = PlantInputs(t=t, Cwo=Cwo, Qtp=Qtp, M=M)

# ---- (1) FULL-ODE 재현 ----------------------------------------------------
print("="*92)
print("(1) FULL-ODE 재현 — basis-A B_k + W2 transport → Yamazaki BAF (Cwo=1)")
print("="*92)
print(f"{'PFAS':8}{'f_xy':>7}{'kappa_d':>9} | {'root pred/obs':>16}{'straw pred/obs':>17}{'grain pred/obs':>17}")
rows1 = []
for r5 in S5:
    p = r5["PFAS"]
    if p not in W2: continue
    w = W2[p]
    cmpd = Compound(name=p, K_prot=float(r5["K_prot"]), K_PL=float(r5["K_PL"]),
                    K_cw=float(KCW[p]["K_cw_wholecw_root"]),
                    kappa_d=float(w["kappa_d_fit"]), Vmax_in=20.0, Km_in=5.0,
                    Vmax_out=8.0, Km_out=5.0, L_Ph=float(w["L_Ph_fit"]), f_xy=float(w["f_xy_fit"]))
    model = RiceUptakeModel(env=Environment(), cmpd=cmpd, comps=comps, inputs=inputs)
    sol = model.solve(t); C = sol.y[:, -1]; Mf = inputs.M_(t[-1])
    pr_root = C[ROOT]; pr_straw = (C[STEM]*Mf[STEM]+C[LEAF]*Mf[LEAF])/(Mf[STEM]+Mf[LEAF]); pr_grain = C[FRUIT]
    o = YA.get(p, {})
    o_root, o_straw, o_grain = o.get("root", np.nan), o.get("straw", np.nan), o.get("grain", np.nan)
    print(f"{p:8}{float(w['f_xy_fit']):>7.3f}{float(w['kappa_d_fit']):>9.3f} | "
          f"{pr_root:>7.2f}/{o_root:<7.2f}{pr_straw:>8.2f}/{o_straw:<7.2f}{pr_grain:>8.2f}/{o_grain:<7.2f}")
    rows1.append([p, r5["n_C"], round(pr_root,3), round(o_root,3), round(pr_straw,3),
                  round(o_straw,3), round(pr_grain,3), round(o_grain,3)])
# residual
arr = np.array([[r[2],r[3],r[4],r[5],r[6],r[7]] for r in rows1], float)
lr = np.log10(np.clip(arr,1e-6,None))
rmse = np.sqrt(np.nanmean([(lr[:,0]-lr[:,1])**2,(lr[:,2]-lr[:,3])**2,(lr[:,4]-lr[:,5])**2]))
print(f"\n  log10 RMSE(pred vs obs, root+straw+grain) = {rmse:.3f}  "
      f"(W2 saturated fit → ≈0 기대; reconstruction coherence 확인)")
with open("S6_Gap4_ode_repro.csv","w",newline="") as f:
    wr=csv.writer(f); wr.writerow(["pfas","n_C","root_pred","root_obs","straw_pred","straw_obs","grain_pred","grain_obs"])
    wr.writerows(rows1)
print("[written] S6_Gap4_ode_repro.csv")

# ---- (2) CROSS-FIELD TF (water-independent) -------------------------------
print("\n" + "="*92)
print("(2) CROSS-FIELD TF_straw = C_straw/C_root  (water 분모 상쇄 → Li2025 artifact 면역)")
print("="*92)
nC = {r["PFAS"]: int(r["n_C"]) for r in S5}
def ya_tf(p, tis):
    o = YA.get(p, {})
    return o.get(tis, np.nan)/o["root"] if "root" in o and tis in o else np.nan
overlap = ["PFBA","PFHxA","PFOA","PFBS","PFOS"]
print(f"{'PFAS':8}{'nC':>3}{'Yamazaki TF_straw':>19}{'Li2025 TF_stem':>16}{'Li2025 TF_grain':>17}")
rows2 = []
allp = sorted(set(list(YA)+list(LI)), key=lambda x: (nC.get(x,99), x))
for p in allp:
    ya_s = ya_tf(p, "straw")
    li_s = float(LITF[p]["TF_stem_fw"]) if p in LITF and LITF[p]["TF_stem_fw"] else np.nan
    li_g = float(LITF[p]["TF_grain_fw"]) if p in LITF and LITF[p]["TF_grain_fw"] else np.nan
    if np.isnan(ya_s) and np.isnan(li_s): continue
    f=lambda v: f"{v:.2f}" if not np.isnan(v) else "  -"
    print(f"{p:8}{nC.get(p,0):>3}{f(ya_s):>19}{f(li_s):>16}{f(li_g):>17}")
    rows2.append([p, nC.get(p,0), round(ya_s,3) if not np.isnan(ya_s) else "",
                  round(li_s,3) if not np.isnan(li_s) else "", round(li_g,3) if not np.isnan(li_g) else ""])
with open("S6_Gap4_TF_crossfield.csv","w",newline="") as f:
    wr=csv.writer(f); wr.writerow(["pfas","n_C","Yamazaki_TF_straw","Li2025_TF_stem","Li2025_TF_grain"])
    wr.writerows(rows2)
print("[written] S6_Gap4_TF_crossfield.csv")
# 단조성 (PFCA)
ca_ya = [(nC[p], ya_tf(p,"straw")) for p in YA if p in nC and S5[[r['PFAS'] for r in S5].index(p)]['group']=='CA' and not np.isnan(ya_tf(p,'straw'))]
ca_ya.sort()
print(f"\n  Yamazaki TF_straw PFCA C4→C12: " + " → ".join(f"{v:.2f}" for _,v in ca_ya) + "  (단조 감소 = TSCF)")
print("  Li2025 TF_stem  C4/C6/C8(PFBA/PFHxA/PFOA): "
      f"{LITF['PFBA']['TF_stem_fw']} → {LITF['PFHxA']['TF_stem_fw']} → {LITF['PFOA']['TF_stem_fw']}  (동일 감소 경향)")

# ---- figures --------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(11.4, 4.4))
# (A) ode pred vs obs
lab=["root","straw","grain"]; col=["#2980b9","#27ae60","#e67e22"]
for j,(lo,hi) in enumerate([(2,3),(4,5),(6,7)]):
    pr=[r[lo] for r in rows1]; ob=[r[hi] for r in rows1]
    ax[0].scatter(ob, pr, s=42, color=col[j], edgecolor="k", lw=.4, label=lab[j], zorder=3)
lims=[5e-2, 2e2]; ax[0].plot(lims,lims,"k--",lw=1,alpha=.6)
ax[0].set_xscale("log"); ax[0].set_yscale("log"); ax[0].set_xlim(lims); ax[0].set_ylim(lims)
ax[0].set_xlabel("observed BAF (Yamazaki)"); ax[0].set_ylabel("ODE predicted BAF")
ax[0].set_title("(1) full-ODE reproduction"); ax[0].legend(fontsize=8); ax[0].grid(alpha=.3,which="both")
# (B) cross-field TF_straw vs chain
yx=[(nC[p], ya_tf(p,"straw")) for p in YA if not np.isnan(ya_tf(p,"straw")) and p in nC]; yx.sort()
ax[1].plot([a for a,_ in yx],[b for _,b in yx],"o-",color="#34495e",lw=2,label="Yamazaki TF_straw")
lx=[(nC[p], float(LITF[p]["TF_stem_fw"])) for p in LITF if LITF[p]["TF_stem_fw"]]; lx.sort()
ax[1].plot([a for a,_ in lx],[b for _,b in lx],"s--",color="#c0392b",lw=2,label="Li2025 TF_stem")
ax[1].axhline(1,ls=":",color="gray"); ax[1].set_yscale("log")
ax[1].set_xlabel("chain length n$_C$"); ax[1].set_ylabel("TF (straw/root)")
ax[1].set_title("(2) TF↓ with chain — both fields (water-independent)")
ax[1].legend(fontsize=8); ax[1].grid(alpha=.3,which="both")
fig.tight_layout(); fig.savefig("S6_Gap4.png", dpi=140); print("[written] S6_Gap4.png")
