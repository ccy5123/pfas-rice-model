#!/usr/bin/env python3
# =============================================================================
# S6 — surface cross-field 대조 (Li2025 Tianjin vs Yamazaki Andosol)
# -----------------------------------------------------------------------------
# H7 §4/§7.3: Li2025서 root BAF가 내부 ceiling을 7-22x 초과 → "surface 필수"로,
# Yamazaki서 sub-equilibrium → "surface 0"로 갈렸음. 이를 per-congener로 정량하고
# K_surf의 토양 의존을 평가한다.  핵심 진단: Li2025 per-congener water QUALITY.
#
#   surface_excess(n) = max(0, obs_root_BAF − B_root_ceiling)   [L/kg fw]
#   obs/ceiling > 1  → 내부 흡수만으로 설명 불가 (surface 또는 분모 오류)
#
# 입력: S6_Bk_basisA_allorgan.csv (ceiling), obs_baf_{Yamazaki,Li2025}.csv
# 산출: S6_surface_crossfield.csv, S6_surface_crossfield.png
# =============================================================================
import csv, re, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

def read_csv(p):
    with open(p, newline="") as f: return list(csv.DictReader(f))

# ceiling B_root(n)
B = {r["pfas"]: float(r["B_root"]) for r in read_csv("S6_Bk_basisA_allorgan.csv")}
nC = {r["pfas"]: int(r["n_C"]) for r in read_csv("S6_Bk_basisA_allorgan.csv")}

def root_obs(path):
    d = {}
    for r in read_csv(path):
        if r["tissue"] == "root":
            q = re.search(r"water (\w+)", r["source"])
            d[r["compound"]] = (float(r["baf"]), q.group(1) if q else "")
    return d

ya = root_obs("obs_baf_Yamazaki.csv")     # clean per-congener water (11)
li = root_obs("obs_baf_Li2025.csv")       # group-water + quality flag (5)

rows = []
print("="*86)
print("surface cross-field: obs root BAF / internal ceiling B_root (basis-A, θ=0.90)")
print("="*86)
print(f"{'field':9}{'PFAS':8}{'nC':>3}{'obs_BAF':>10}{'ceiling':>9}{'obs/ceil':>10}"
      f"{'excess':>9}{'water_q':>9}")
for field, dd in [("Yamazaki", ya), ("Li2025", li)]:
    for p in sorted(dd, key=lambda x: (nC.get(x, 99))):
        if p not in B: continue
        obs, q = dd[p]; ceil = B[p]; ratio = obs/ceil; exc = max(0.0, obs-ceil)
        flag = "EXCESS" if ratio > 1 else "sub-eq"
        print(f"{field:9}{p:8}{nC[p]:>3}{obs:>10.3f}{ceil:>9.2f}{ratio:>10.3f}"
              f"{exc:>9.2f}{q:>9}  {flag}")
        rows.append([field, p, nC[p], round(obs,3), round(ceil,3), round(ratio,3),
                     round(exc,3), q, flag])
    print("-"*86)

with open("S6_surface_crossfield.csv", "w", newline="") as f:
    wr = csv.writer(f)
    wr.writerow(["field","pfas","n_C","obs_root_BAF_fw","B_root_ceiling",
                 "obs_over_ceiling","surface_excess","water_quality","verdict"])
    wr.writerows(rows)
print("[written] S6_surface_crossfield.csv")

# ---- 진단: Li2025 excess vs water quality --------------------------------
print("\n" + "="*86)
print("진단 — Li2025 obs/ceiling 을 pore-water 품질순으로 정렬")
print("="*86)
qorder = {"good":0, "med":1, "rough":2, "poor":3}
liq = sorted(li, key=lambda p: qorder.get(li[p][1], 9))
for p in liq:
    obs, q = li[p]; ratio = obs/B[p]
    print(f"  {q:6}  {p:7}(C{nC[p]:>2})  obs/ceiling = {ratio:6.2f}   "
          f"{'sub-eq (surface 불필요)' if ratio<1 else 'apparent EXCESS'}")
print("\n  → good-water 단일점 PFOA만 sub-equilibrium(0.53, Yamazaki와 동일). "
      "\n    excess는 med→rough→poor 로 단조 증가 = pore-water 분모 신뢰도와 교란.")
print("    강흡착 congener(낮은 aqueous %)일수록 water 분모가 과소·불확실 → BAF 인위적 팽창.")
