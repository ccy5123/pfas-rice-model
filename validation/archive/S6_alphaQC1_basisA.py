#!/usr/bin/env python3
# =============================================================================
# S6 — alpha / QC1 RE-ANALYSIS on basis-A  (정정 basis에서 막 down-weight 재검)
# -----------------------------------------------------------------------------
# H7가 무효화한 Prompt3 v3(naive basis: B = theta + sum f_i K_i)를
# basis-A 로 옮긴다:
#       B_k = theta_fw + (1 - theta_fw) * sum_i f_{i,dw} * K_i      [L/kg fw]
# (module binding_factors / rice_tissue NOTES path-A 와 동일)
#
# 입력(검증 파일에서 직접 로드):
#   K_PL(n), K_prot(n)        <- inputs/Bk_table_S5.csv
#   K_cw whole-cw per organ   <- Kcw_Klignin_params_v2.csv
#   기관 조성(theta,f_prot,f_PL,f_cw)  <- rice_tissue_params.csv (recommended)
#   관측 root BAF (fw)        <- obs_baf_Yamazaki.csv
#
# 산출:
#   (1) S6_Bk_basisA_allorgan.csv     basis-A B_k(n) root/stem/leaf/grain (alpha=1)
#   (2) S6_alphaQC1_root_basisA.csv   root pool-share + alpha 식별성(min/F2)
#   (3) S6_alphaQC1_root.png          막분율 & alpha 경계 figure
# =============================================================================
import csv, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- Liu2023 wheat 근부 subcellular target bands (Prompt3 v3 동일) ----------
F1 = (47.7, 59.0)   # cell wall
F2 = (18.3, 43.0)   # membrane + organelle
F3 = ( 6.0, 34.0)   # soluble / cytosol

# ---- 1. 입력 로드 ----------------------------------------------------------
def read_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))

S5   = read_csv("prompt3_pkg/Prompt3_Kpl_downweight/inputs/Bk_table_S5.csv")
KCW  = {r["pfas"]: r for r in read_csv("Kcw_Klignin_params_v2.csv")}
RT   = read_csv("rice_tissue_pkg/rice_tissue/rice_tissue_params.csv")
OBSr = read_csv("obs_baf_Yamazaki.csv")

# congener 순서 = S5 (12종, PFHxS 포함)
pfas = [r["PFAS"] for r in S5]
nC   = np.array([int(r["n_C"]) for r in S5])
K_PL = np.array([float(r["K_PL"])   for r in S5])
K_pr = np.array([float(r["K_prot"]) for r in S5])
K_cw = np.array([float(KCW[p]["K_cw_wholecw_root"]) for p in pfas])  # rice whole-cw root

# 관측 root BAF (fw) — 11종(PFHxS ND)
obs_root = {}
for r in OBSr:
    if r["tissue"] == "root":
        obs_root[r["compound"]] = float(r["baf"])
obs_arr = np.array([obs_root.get(p, np.nan) for p in pfas])

# 기관 조성(recommended): 4-pool f_cw = poly + lignin
def rt(organ, param):
    for r in RT:
        if r["organ"] == organ and r["parameter"] == param:
            return float(r["value_recommended"])
    return np.nan

ORG = {}
for organ, key in [("root","root"),("stem","stem"),("leaf","leaf"),("grain_brown","grain")]:
    th  = rt(organ, "theta_fw")
    fpr = rt(organ, "f_prot")
    fPL = rt(organ, "f_PL_membrane")
    fcw = rt(organ, "f_cw_polysaccharide") + rt(organ, "lignin")
    ORG[key] = dict(theta=th, f_prot=fpr, f_PL=fPL, f_cw=fcw)
print("기관 조성 (rice_tissue recommended, 4-pool f_cw = poly+lignin):")
for k, v in ORG.items():
    print(f"  {k:6}  theta={v['theta']:.2f}  f_prot={v['f_prot']:.3f} "
          f"f_PL={v['f_PL']:.4f}  f_cw={v['f_cw']:.3f}")

# ---- 2. basis-A binding & 분해 --------------------------------------------
def comp_basisA(theta, fpr, fPL, fcw, Kcw, aPL=1.0):
    """basis-A pool 분해. 반환 (pools[N,4]=[water,prot,mem,cw], B[N])."""
    w  = np.full_like(K_PL, theta)
    pr = (1 - theta) * fpr  * K_pr
    me = (1 - theta) * aPL * fPL * K_PL
    c  = (1 - theta) * fcw  * Kcw
    B  = w + pr + me + c
    return np.vstack([w, pr, me, c]).T, B

def pct(P, B): return P / B[:, None] * 100.0
def band(v, lo, hi): return "OK" if lo <= v <= hi else ("LOW" if v < lo else "HIGH")
def res(v, lo, hi):  return 0.0 if lo <= v <= hi else min(abs(v-lo), abs(v-hi))

# whole-cw per organ (장쇄까지 organ별)
KCW_ORG = {k: np.array([float(KCW[p][f"K_cw_wholecw_{k if k!='grain' else 'grain_brown'}"])
                        for p in pfas]) for k in ORG}

# ---- (1) basis-A B_k(n) all-organ (alpha=1, W1 재생성) ---------------------
Ball = {}
for k, v in ORG.items():
    _, B = comp_basisA(v["theta"], v["f_prot"], v["f_PL"], v["f_cw"], KCW_ORG[k], 1.0)
    Ball[k] = B
with open("S6_Bk_basisA_allorgan.csv", "w", newline="") as f:
    wr = csv.writer(f)
    wr.writerow(["pfas","n_C","group","B_root","B_stem","B_leaf","B_grain",
                 "theta_root","theta_stem","theta_leaf","theta_grain"])
    for i, p in enumerate(pfas):
        wr.writerow([p, nC[i], S5[i]["group"],
                     round(Ball["root"][i],3), round(Ball["stem"][i],3),
                     round(Ball["leaf"][i],3), round(Ball["grain"][i],3),
                     ORG["root"]["theta"], ORG["stem"]["theta"],
                     ORG["leaf"]["theta"], ORG["grain"]["theta"]])
print("\n[written] S6_Bk_basisA_allorgan.csv")
print("  basis-A B_root (theta=0.90):  " +
      "  ".join(f"{p}={Ball['root'][i]:.2f}" for i, p in enumerate(pfas) if p in ("PFOA","PFOS","PFDoDA")))

# ---- (A) 막분율(QC1) basis-A 에서 지속하는가 (root, alpha=1) ----------------
r = ORG["root"]
P0, B0 = comp_basisA(r["theta"], r["f_prot"], r["f_PL"], r["f_cw"], K_cw, 1.0)
p0 = pct(P0, B0); s_mem0 = P0[:, 2] / B0
print("\n" + "="*78)
print("(A) basis-A root 막분율 (alpha=1, theta=0.90) — QC1 지속 여부")
print("="*78)
print(f"{'PFAS':8}{'B_root':>9}{'water%':>8}{'prot%':>7}{'MEM%':>7}{'cw%':>7}")
for i, p in enumerate(pfas):
    print(f"{p:8}{B0[i]:>9.2f}{p0[i,0]:>8.1f}{p0[i,1]:>7.1f}{p0[i,2]:>7.1f}{p0[i,3]:>7.1f}")
long = [i for i,p in enumerate(pfas) if nC[i] >= 8 and S5[i]['group']=='CA']
print(f"  → 장쇄(C8+ PFCA) 막분율 범위: {p0[long,2].min():.0f}–{p0[long,2].max():.0f}%  "
      f"(naive basis 75–98%와 동일 정성 — 막 지배 지속)")

# ---- (B) Prompt3 v3 '권장구성' 잔차를 basis-A 로 (naive 2.3 와 대조) --------
def alpha_to_membrane(theta, fpr, fPL, fcw, Kcw, s_t):
    P, B = comp_basisA(theta, fpr, fPL, fcw, Kcw, 1.0)
    s = P[:, 2] / B
    return np.minimum((s_t/(1-s_t)) * ((1-s)/s), 1.0)

def total_resid(theta, s_t=0.30):
    a = alpha_to_membrane(theta, r["f_prot"], r["f_PL"], r["f_cw"], K_cw, s_t)
    P, B = comp_basisA(theta, r["f_prot"], r["f_PL"], r["f_cw"], K_cw, a)
    p = pct(P, B); cw, mem, f3 = p[:,3], p[:,2], p[:,0]+p[:,1]
    f1r = np.mean([res(v,*F1) for v in cw])
    f2r = np.mean([res(v,*F2) for v in mem])
    f3r = np.mean([res(v,*F3) for v in f3])
    return f1r+f2r+f3r, cw.min(), cw.max(), a
print("\n" + "="*78)
print("(B) v3 권장구성(f_cw=0.50 + rice K_cw + alpha→mem0.30)을 basis-A 로 재계산")
print("    — naive basis 총잔차 2.3 (v3 결론)과 대조")
print("="*78)
print(f"{'theta':>6}{'tot_resid(|F1|+|F2|+|F3|)':>26}{'cw_min%':>9}{'cw_max%':>9}")
for theta in [0.70, 0.90, 0.92]:
    tr, cmn, cmx, _ = total_resid(theta)
    print(f"{theta:>6.2f}{tr:>26.1f}{cmn:>9.0f}{cmx:>9.0f}   (F1 하한 {F1[0]})")
print(f"  → basis-A 에서 cw%가 F1 밴드(47.7–59.0)에 도달 못함 = v3의 'f_cw가 QC1 해결' 불성립.")
print(f"    (1-theta) 인자가 cw/prot/mem 결합 pool 전체를 축소하므로 막→0.30 보정 후에도 cw share가 부족.")

# ---- (C) ceiling 기반 alpha 식별성 (Yamazaki sub-equilibrium) ---------------
#   B_root(alpha) >= obs_root_BAF  (관측은 ceiling 이하 = 흡수 한계)
#   alpha_min : B_root(alpha)=obs 를 푸는 값. <=0 이면 비막 floor 가 이미 obs 이상 (하한 없음).
floor_nonmem = B0 - P0[:, 2]           # alpha=0 일 때 B_root (막 제외)
denom = (1 - r["theta"]) * r["f_PL"] * K_PL
alpha_min = (obs_arr - floor_nonmem) / denom
alpha_min_clip = np.clip(alpha_min, 0.0, None)
# alpha_F2 : Liu 막→0.30 anchor (basis-A)
alpha_F2 = alpha_to_membrane(r["theta"], r["f_prot"], r["f_PL"], r["f_cw"], K_cw, 0.30)

print("\n" + "="*78)
print("(C) alpha 식별성: ceiling 하한(alpha_min) vs Liu-F2 anchor(alpha_F2), root")
print("="*78)
print(f"{'PFAS':8}{'obs_BAF':>9}{'B_root(a=1)':>12}{'obs/ceil':>10}"
      f"{'alpha_min':>11}{'alpha_F2':>10}{'note':>8}")
rows = []
for i, p in enumerate(pfas):
    if np.isnan(obs_arr[i]):
        amn, ratio, note = np.nan, np.nan, "ND"
    else:
        ratio = obs_arr[i] / B0[i]
        if alpha_min[i] <= 0:
            amn, note = 0.0, "floor≥obs"          # 막 하한 제약 없음
        elif alpha_min[i] > 1:
            amn, note = np.nan, "INFEAS"           # ceiling<obs (surface 신호)
        else:
            amn, note = alpha_min[i], "bound"
    obs_s = f"{obs_arr[i]:.2f}" if not np.isnan(obs_arr[i]) else "  ND"
    amn_s = f"{amn:.3f}" if not np.isnan(amn) else "   -"
    print(f"{p:8}{obs_s:>9}{B0[i]:>12.2f}{(ratio if not np.isnan(ratio) else float('nan')):>10.3f}"
          f"{amn_s:>11}{alpha_F2[i]:>10.3f}{note:>8}")
    rows.append([p, nC[i], S5[i]["group"], obs_arr[i], round(B0[i],3),
                 round(ratio,3) if not np.isnan(ratio) else "",
                 round(amn,4) if not np.isnan(amn) else "", round(alpha_F2[i],4), note,
                 round(s_mem0[i]*100,1), round(p0[i,0],1), round(p0[i,1],1),
                 round(p0[i,2],1), round(p0[i,3],1)])

with open("S6_alphaQC1_root_basisA.csv", "w", newline="") as f:
    wr = csv.writer(f)
    wr.writerow(["pfas","n_C","group","obs_root_BAF_fw","B_root_basisA_a1","obs_over_ceiling",
                 "alpha_min_ceiling","alpha_F2_LiuAnchor","feasibility_note",
                 "mem_share_a1_%","water_%","prot_%","mem_%","cw_%"])
    wr.writerows(rows)
print("\n[written] S6_alphaQC1_root_basisA.csv")

# ---- (3) figure ------------------------------------------------------------
ca = [i for i in range(len(pfas)) if S5[i]["group"] == "CA"]
fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
ax[0].plot(nC[ca], s_mem0[ca]*100, "o-", color="#c0392b", lw=2, label="membrane share (α=1)")
ax[0].axhspan(F2[0], F2[1], color="#2980b9", alpha=0.15, label="Liu F2 band (18–43%)")
ax[0].set_xlabel("perfluorocarbon chain length n$_C$"); ax[0].set_ylabel("share of B$_{root}$ [%]")
ax[0].set_title("(A) basis-A membrane domination (PFCA, root)"); ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)
amn_plot = np.array([alpha_min_clip[i] if not np.isnan(alpha_min[i]) and alpha_min[i]<=1 else np.nan for i in ca])
ax[1].plot(nC[ca], alpha_F2[ca], "s-", color="#2980b9", lw=2, label=r"$\alpha_{F2}$ (Liu mem→0.30)")
ax[1].plot(nC[ca], amn_plot, "^--", color="#27ae60", lw=2, label=r"$\alpha_{min}$ (ceiling≥obs)")
ax[1].fill_between(nC[ca], amn_plot, 1.0, color="#27ae60", alpha=0.10)
ax[1].set_yscale("log"); ax[1].set_ylim(1e-3, 2)
ax[1].set_xlabel("perfluorocarbon chain length n$_C$"); ax[1].set_ylabel(r"$\alpha$ (K$_{PL}$ down-weight)")
ax[1].set_title("(B) α bounds: transport vs subcellular anchor"); ax[1].legend(fontsize=8); ax[1].grid(alpha=.3, which="both")
fig.tight_layout(); fig.savefig("S6_alphaQC1_root.png", dpi=140)
print("[written] S6_alphaQC1_root.png")
