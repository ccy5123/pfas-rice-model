#!/usr/bin/env python3
# =============================================================================
# validation/policy_brief_runs.py
# -----------------------------------------------------------------------------
# Produces every number in `docs/POLICY_BRIEF_KR.md` -- the briefing written so a
# slide-building agent can turn it into a deck for a policy council without
# re-deriving anything or over-claiming.
#
# WHY A SCRIPT AND NOT JUST A DOCUMENT. A briefing whose numbers cannot be
# regenerated is a briefing nobody can check a year later. Every figure quoted in
# the report is printed here, written to `docs/policy_brief_results.csv`, and
# pinned by `tests/test_policy_brief.py`.
#
# WHICH MODEL. The DEFAULT configuration -- the redistributed shoot with K_PL-gated
# lipid loading (`model_api.DEFAULT_SHOOT` / `DEFAULT_LIPID_LOADING`, adopted by the
# pre-registered rule in validation/default_switch.py) -- on the calibration
# site-year. The previous default (4-pool, free anion: `model_api.FOUR_POOL_FREE`) is
# computed beside it wherever the briefing contrasts the two. All tissue values are
# DRY basis (as the measurements are reported) except the grain "as eaten" used for
# the intake comparison.
#
# WHAT IT RUNS
#   1. the 13 curated PFAS: where each ends up, and what a flooded schedule does
#   2. one substance end to end (PFOA at Kim 2019's field pore water, and at 50 ng/L)
#   3. the Korean FIELD dataset with paired pore water and brown rice (Kim 2019),
#      scored with its non-detects (validation/kim2019_censored.py)
#   4. the inverse question (grain measurement -> pore water), against a known
#      synthetic truth and against that field measurement
#   5. structure (SMILES) input, including substances outside the curated list
#   6. the EU biocide census, read back from validation/bat_census_biocides.csv
#
#   python validation/policy_brief_runs.py            (~6 min, writes CSV + 2 figures)
#   python validation/policy_brief_runs.py --fast     (skips the figures)
# =============================================================================
from __future__ import annotations

import contextlib
import csv
import io
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT_DIR, "src"))
sys.path.insert(0, HERE)

import model_api as api                     # noqa: E402
import literature_params as lp              # noqa: E402

OUT_CSV = os.path.join(ROOT_DIR, "docs", "policy_brief_results.csv")
FIG_DIR = os.path.join(HERE, "figures")
PREV = dict(api.FOUR_POOL_FREE)             # the previous default, for contrast

# A realistic Korean paddy pore-water level for the worked example: the PFOA field
# average of Kim et al. 2019 (78.7 ng/L = 0.0787 ug/L).
KIM_PFOA_POREWATER = 0.0787

_ROWS: list[tuple] = []


def rec(section, key, value, unit="", note=""):
    _ROWS.append((section, key, value, unit, note))
    return value


# ---------------------------------------------------------------------------
# 1. the 13 curated PFAS -- where each one ends up
# ---------------------------------------------------------------------------
def section1(quiet=False):
    """BAF at a constant 1 ug/L pore water (the generic table), the flooded-schedule
    grain, and -- for the Yamazaki congeners -- the field observable (two harvest
    runs at the irrigation-water level) against the measurement, for the default and
    the previous default."""
    out = []
    obs_all = api.yamazaki_obs()
    for c in api.chain_table():
        n = c["name"]
        base = api.simulate(n)
        flood = api.simulate(n, cwo_profile="flooded")
        prev = api.simulate(n, **PREV)
        d = dict(name=n, n_C=c["n_C"], group=c["group"],
                 root=base["baf_final_dw"]["root"], straw=base["baf_final_dw"]["straw"],
                 grain=base["baf_final_dw"]["grain"],
                 grain_prev=prev["baf_final_dw"]["grain"],
                 grain_flooded=flood["baf_final_dw"]["grain"],
                 obs_grain=None, yz_grain=None, yz_grain_prev=None)
        if n in obs_all:
            o = api.observed_baf(n)
            d["obs_grain"] = o.get("grain")
            d["yz_grain"] = api.yamazaki_harvest_baf(n)["grain"]
            d["yz_grain_prev"] = api.yamazaki_harvest_baf(n, **PREV)["grain"]
        d["flood_ratio"] = d["grain_flooded"] / d["grain"]
        out.append(d)
        for k in ("root", "straw", "grain", "grain_prev", "grain_flooded"):
            rec("1_congeners", f"{n}.{k}", d[k], "L/kg dw")
        for k in ("obs_grain", "yz_grain", "yz_grain_prev"):
            if d[k] is not None:
                rec("1_congeners", f"{n}.{k}", d[k], "L/kg dw",
                    "Yamazaki 2023 field (measured)" if k == "obs_grain" else
                    "model, Yamazaki observable")
    short = [d for d in out if d["n_C"] <= 6]
    long_ = [d for d in out if d["n_C"] >= 9]
    for tag, key in (("base", "grain"), ("prev", "grain_prev")):
        rs = float(np.mean([d[key] for d in short]))
        rl = float(np.mean([d[key] for d in long_]))
        rec("1_congeners", f"mean_grain_shortchain_C4_C6_{tag}", rs, "L/kg dw")
        rec("1_congeners", f"mean_grain_longchain_C9_C12_{tag}", rl, "L/kg dw")
    yz = [d for d in out if d["obs_grain"]]
    over = {d["name"]: d["yz_grain"] / d["obs_grain"] for d in yz}
    under_prev = {d["name"]: d["obs_grain"] / d["yz_grain_prev"] for d in yz}
    for nm, f in over.items():
        rec("1_congeners", f"{nm}.yamazaki_grain_model_over_obs", f, "x")
    for nm, f in under_prev.items():
        rec("1_congeners", f"{nm}.yamazaki_grain_obs_over_prev", f, "x")
    if not quiet:
        print("=" * 78)
        print("1. 13종 PFAS — 물질마다 도달 부위가 다르다  [BAF, L/kg 건물중, 공극수 1 µg/L]")
        print("=" * 78)
        print(f"{'물질':9s}{'C':>3s}{'계열':>8s}{'뿌리':>9s}{'짚':>9s}{'낟알':>9s}"
              f"{'낟알(이전)':>11s}{'담수낟알':>10s}")
        for d in sorted(out, key=lambda x: (x["group"], x["n_C"])):
            print(f"{d['name']:9s}{d['n_C']:3d}{d['group']:>8s}{d['root']:9.3f}"
                  f"{d['straw']:9.3f}{d['grain']:9.3f}{d['grain_prev']:11.3f}"
                  f"{d['grain_flooded']:10.3f}")
        print("\n  Yamazaki 현장 관측(관개수 농도, 두 수확) — 현미 BAF: 실측 | 기본 | 이전 기본")
        for d in yz:
            print(f"    {d['name']:8s} {d['obs_grain']:8.3g} | {d['yz_grain']:8.3g} | "
                  f"{d['yz_grain_prev']:8.3g}   (기본/실측 {over[d['name']]:.1f}배)")
    return out


# ---------------------------------------------------------------------------
# 2. one substance, end to end
# ---------------------------------------------------------------------------
def section2(quiet=False, cwo=KIM_PFOA_POREWATER, tag="2_pfoa"):
    r = api.simulate("PFOA", Cwo=cwo)
    grain = float(r["conc"]["grain"][-1])                 # as eaten (grain water 14 %)
    dw = r["conc_dw"]
    ap = api.apportionment("PFOA", Cwo=cwo)["fraction"]
    d = dict(cwo=cwo, root=float(dw["root"][-1]), straw=float(dw["straw"][-1]), grain=grain,
             grain_dw=float(dw["grain"][-1]), band=api.predictive_band(grain),
             intake=api.intake_fraction(grain, congener="PFOA"), paths=ap)
    rec(tag, "cwo", cwo, "ug/L")
    for k, u in (("root", "ug/kg dw"), ("straw", "ug/kg dw"), ("grain", "ug/kg as eaten"),
                 ("grain_dw", "ug/kg dw")):
        rec(tag, f"conc.{k}", d[k], u)
    rec(tag, "grain.band_lo", d["band"]["lo"], "ug/kg as eaten")
    rec(tag, "grain.band_hi", d["band"]["hi"], "ug/kg as eaten")
    rec(tag, "grain.band_factor", d["band"]["factor"], "x")
    rec(tag, "efsa_twi_percent", d["intake"]["percent"], "%",
        "EFSA 2020 4종 합 TWI 4.4 ng/kg bw/주, 쌀 150 g/일, 체중 60 kg")
    rec(tag, "grain_via_xylem", ap["grain"]["xylem_from_stem"], "-")
    rec(tag, "grain_via_phloem", ap["grain"]["phloem_from_leaf"], "-")
    rec(tag, "root_from_soil", ap["root"]["soil_uptake"], "-")
    if not quiet:
        print()
        print("=" * 78)
        print(f"2. 한 물질을 끝까지 — PFOA, 공극수 {cwo * 1000:.1f} ng/L")
        print("=" * 78)
        print(f"  뿌리 {d['root']:.4f} · 짚 {d['straw']:.4f} µg/kg 건물중 · "
              f"낟알 {d['grain']:.4f} µg/kg (먹는 상태)")
        print(f"  낟알 불확실 범위 {d['band']['lo']:.4f} – {d['band']['hi']:.4f}"
              f"  (×/÷{d['band']['factor']:.1f}, 전이 오차 중앙값 log10 {d['band']['log10_rmse']})")
        print(f"  EFSA 주간섭취허용량 대비 {d['intake']['percent']:.1f}%")
        print(f"  낟알 유입: 물관 {ap['grain']['xylem_from_stem'] * 100:.0f}%"
              f" · 체관 {ap['grain']['phloem_from_leaf'] * 100:.0f}%;"
              f" 뿌리 유입: 토양 {ap['root']['soil_uptake'] * 100:.0f}%")
    return d


# ---------------------------------------------------------------------------
# 3. the Korean field dataset, with its non-detects
# ---------------------------------------------------------------------------
def section3(quiet=False):
    import kim2019_censored as K
    with contextlib.redirect_stdout(io.StringIO()):
        res = K.main()
    b, preds, scores = res["bounds"], res["preds"], res["scores"]
    rows = []
    for nm in K.CHAIN:
        pw, soil, rice, df = lp.KIM2019_FIELD[nm]
        lo, pub, hi = b[nm]
        rows.append(dict(name=nm, porewater_ngL=pw, brownrice_ugkg=rice, DF_pct=df,
                         obs_baf=pub, lo=lo, hi=hi,
                         base=preds["nstem_lipid"][nm], prev=preds["mono"][nm]))
        rec("3_kim2019", f"{nm}.porewater", pw, "ng/L", "Kim 2019 field mean (ND as 0)")
        rec("3_kim2019", f"{nm}.brownrice", rice, "ug/kg", "Kim 2019 field mean (ND as 0); basis not stated")
        rec("3_kim2019", f"{nm}.detection_freq", df, "%", "brown rice, 30 fields")
        rec("3_kim2019", f"{nm}.obs_grain_baf", pub, "L/kg")
        rec("3_kim2019", f"{nm}.obs_lo", lo, "L/kg", "every non-detect at 0")
        rec("3_kim2019", f"{nm}.obs_hi", hi, "L/kg", "every non-detect at the MDL")
        rec("3_kim2019", f"{nm}.model_grain_baf", preds["nstem_lipid"][nm], "L/kg",
            "pore water 1 ug/L, harvest moisture")
        rec("3_kim2019", f"{nm}.model_grain_baf_prev", preds["mono"][nm], "L/kg")
    s, p = scores["nstem_lipid"], scores["mono"]
    rec("3_kim2019", "rmse_point_default", s["point"], "log10", "PFOA excluded")
    rec("3_kim2019", "rmse_interval_default", s["interval"], "log10", "PFOA excluded")
    rec("3_kim2019", "rmse_point_prev", p["point"], "log10", "PFOA excluded")
    rec("3_kim2019", "rmse_interval_prev", p["interval"], "log10", "PFOA excluded")
    g, glo, ghi = lp.KIM2019_PFOA_SITE_BAF
    rec("3_kim2019", "pfoa_site_geomean", g, "L/kg", "censored per-site geomean (kim2019_pfoa_anchor.py)")
    rec("3_kim2019", "pfoa_site_geomean_lo90", glo, "L/kg")
    rec("3_kim2019", "pfoa_site_geomean_hi90", ghi, "L/kg")
    if not quiet:
        print()
        print("=" * 78)
        print("3. 한국 논 현장 자료 — Kim et al. 2019 (공극수·현미 쌍, 불검출 포함)")
        print("=" * 78)
        print(f"{'물질':9s}{'검출률%':>8s}{'실측 BAF':>10s}{'[0, MDL] 구간':>18s}{'기본':>9s}{'이전':>9s}")
        for r in rows:
            print(f"{r['name']:9s}{r['DF_pct']:8.1f}{r['obs_baf']:10.2f}"
                  f"   [{r['lo']:6.2f}, {r['hi']:7.2f}]{r['base']:9.3g}{r['prev']:9.3g}")
        print(f"\n  log10 RMSE (PFOA 제외, 점 / 구간): 기본 {s['point']:.2f} / {s['interval']:.2f}"
              f"   ·   이전 기본 {p['point']:.2f} / {p['interval']:.2f}")
    return rows, dict(point=s["point"], interval=s["interval"],
                      prev_point=p["point"], prev_interval=p["interval"])


# ---------------------------------------------------------------------------
# 4. the inverse question
# ---------------------------------------------------------------------------
def section4(quiet=False):
    truth = KIM_PFOA_POREWATER
    r = api.simulate("PFOA", Cwo=truth)
    dw = r["conc_dw"]
    meas = {"root": float(dw["root"][-1]), "straw": float(dw["straw"][-1]),
            "grain": float(dw["grain"][-1])}
    syn = api.estimate_exposure_bayesian("PFOA", meas)
    field = api.estimate_exposure_bayesian("PFOA", {"grain": 0.349})
    prev = api.estimate_exposure_bayesian("PFOA", {"grain": 0.349}, **PREV)
    deck = api.estimate_exposure_bayesian("PFOA", {"grain": 0.5})
    rec("4_inverse", "deck_example.grain", 0.5, "ug/kg dw")
    rec("4_inverse", "deck_example.median", deck["median"], "ug/L")
    rec("4_inverse", "deck_example.ci95_lo", deck["ci95"][0], "ug/L")
    rec("4_inverse", "deck_example.ci95_hi", deck["ci95"][1], "ug/L")
    rec("4_inverse", "synthetic.truth", truth, "ug/L")
    rec("4_inverse", "synthetic.median", syn["median"], "ug/L")
    rec("4_inverse", "synthetic.ci95_lo", syn["ci95"][0], "ug/L")
    rec("4_inverse", "synthetic.ci95_hi", syn["ci95"][1], "ug/L")
    rec("4_inverse", "field.measured_brownrice", 0.349, "ug/kg", "Kim 2019 PFOA 현미 평균")
    rec("4_inverse", "field.median", field["median"], "ug/L")
    rec("4_inverse", "field.ci95_lo", field["ci95"][0], "ug/L")
    rec("4_inverse", "field.ci95_hi", field["ci95"][1], "ug/L")
    rec("4_inverse", "field.overestimate_factor", field["median"] / truth, "x")
    rec("4_inverse", "field.overestimate_factor_prev", prev["median"] / truth, "x")
    if not quiet:
        print()
        print("=" * 78)
        print("4. 거꾸로 풀기 — 쌀 측정값에서 논 공극수 추정")
        print("=" * 78)
        print(f"  (가) 정답을 아는 경우: 참값 {truth:.4f} → 추정 {syn['median']:.4f} µg/L"
              f"  (95% {syn['ci95'][0]:.4f}–{syn['ci95'][1]:.4f})")
        print(f"  (나) 현장값: 현미 0.349 µg/kg → {field['median']:.3f} µg/L"
              f"  (95% {field['ci95'][0]:.3f}–{field['ci95'][1]:.3f}); 실측 {truth} → "
              f"{field['median'] / truth:.1f}배 (이전 기본 {prev['median'] / truth:.1f}배)")
    return syn, field


# ---------------------------------------------------------------------------
# 5. structure input
# ---------------------------------------------------------------------------
def section5(quiet=False):
    cases = [("PFOA", "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
              "목록에 있는 물질 — 구조로 넣어도 실측 파라미터를 그대로 불러온다"),
             ("PFPeS", "OS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
              "C5 술폰산 — 목록에 없음, QSPR 추정(잠정)"),
             ("ADONA형 에터", "OC(=O)C(F)(F)OC(F)(F)C(F)(F)OC(F)(F)C(F)(F)F",
              "에터 결합 신규 구조 — QSPR 추정(잠정)")]
    out = []
    if not api.rdkit_available():
        if not quiet:
            print("\n5. 구조(SMILES) 입력 — RDKit 미설치, 건너뜀")
        return out
    for label, smi, note in cases:
        r = api.simulate_from_smiles(smi)
        b = r["baf_final_dw"]
        out.append(dict(label=label, matched=r.get("congener"), provisional=bool(r.get("provisional")),
                        root=b["root"], straw=b["straw"], grain=b["grain"], note=note))
        for k in ("root", "straw", "grain"):
            rec("5_smiles", f"{label}.{k}", b[k], "L/kg dw", note)
        rec("5_smiles", f"{label}.provisional", float(bool(r.get("provisional"))), "-")
    if not quiet:
        print()
        print("=" * 78)
        print("5. 구조(SMILES)로 넣기 — 목록에 없는 물질도 계산된다")
        print("=" * 78)
        for d in out:
            flag = "잠정(QSPR)" if d["provisional"] else "실측 파라미터"
            print(f"  {d['label']:16s} 인식 {d['matched']:8s} [{flag}]  뿌리 {d['root']:.3f}"
                  f" · 짚 {d['straw']:.3f} · 낟알 {d['grain']:.3f}")
    return out


# ---------------------------------------------------------------------------
# 6. the EU biocide census (read back, not re-run)
# ---------------------------------------------------------------------------
def section6(quiet=False):
    path = os.path.join(HERE, "bat_census_biocides.csv")
    if not os.path.exists(path):
        if not quiet:
            print("\n6. EU 살생물질 — validation/bat_census_biocides.py 를 먼저 실행하세요")
        return None
    rows = list(csv.DictReader(open(path)))
    excl = [r for r in rows if r["verdict"].startswith("EXCLUDED")]
    inside = [r for r in rows if r["verdict"] == "inside"]
    rec("6_biocides", "n_runnable", len(rows), "행")
    rec("6_biocides", "n_inside_scope", len(inside), "행")
    rec("6_biocides", "n_excluded_ionisation", len(excl), "행")
    # the rank correlation and the bioaccumulative group's TSCF come from the biocide figure
    # script's own record (validation/policy_biocides_figures.py), read back, not recomputed
    fig_csv = os.path.join(ROOT_DIR, "docs", "policy_biocides_results.csv")
    if os.path.exists(fig_csv):
        got = {(r["section"], r["key"]): r for r in csv.DictReader(open(fig_csv, encoding="utf-8"))}
        for (sec, key), name in [(("3_rank", "rho_root"), "rho_fishBCF_root"),
                                 (("3_rank", "rho_straw"), "rho_fishBCF_straw"),
                                 (("3_rank", "n"), "rho_n"),
                                 (("4_bioacc", "max_tscf"), "max_tscf_bioaccumulative")]:
            if (sec, key) in got:
                rec("6_biocides", name, float(got[(sec, key)]["value"]), got[(sec, key)]["unit"],
                    "from docs/policy_biocides_results.csv")
    if not quiet:
        print()
        print("=" * 78)
        print("6. EU 살생물질 목록 — 규제 분류와 쌀 노출은 다른 질문")
        print("=" * 78)
        print(f"  계산 {len(rows)}행 · 적용범위 내 {len(inside)}행 · 이온화로 제외 {len(excl)}행")
    return dict(n=len(rows), inside=len(inside), excluded=len(excl))


# ---------------------------------------------------------------------------
def figures(cong, kim):
    """Two figures, English labels only (no Korean font ships with matplotlib here).
    The briefing carries the Korean captions."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(FIG_DIR, exist_ok=True)
    BASE, PREVC, OBS = "#00808F", "#9A8F7A", "#2B2B2B"

    # (a) Yamazaki field grain BAF vs chain length: measured, default, previous default
    pfca = sorted([d for d in cong if d["group"] == "PFCA" and d["obs_grain"]],
                  key=lambda x: x["n_C"])
    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    x = [d["n_C"] for d in pfca]
    ax.plot(x, [d["obs_grain"] for d in pfca], "o", color=OBS, ms=8, mec="white", mew=1.2,
            label="measured, Yamazaki 2023 (paddy field, Tsukuba)", zorder=3)
    ax.plot(x, [d["yz_grain"] for d in pfca], "s-", color=BASE, lw=1.8, ms=6,
            label="model, default (redistributed shoot + lipid loading)")
    ax.plot(x, [d["yz_grain_prev"] for d in pfca], "o--", color=PREVC, lw=1.4, ms=5, mfc="white",
            label="model, previous default (4-pool, free anion)")
    ax.set_yscale("log")
    ax.set_xlabel("perfluorocarbon chain length  (n C)")
    ax.set_ylabel("brown-rice BAF  [L/kg dry]  (log scale)")
    ax.set_title("Grain BAF along the PFCA series: the default follows the long-chain rise\n"
                 "but sits about 2-7x above most measured values (in-sample for its lipid constants)",
                 fontsize=10.5)
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    ax.grid(alpha=.25, which="both")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.tight_layout()
    p1 = os.path.join(FIG_DIR, "policy_grain_by_chain.png")
    fig.savefig(p1, dpi=160)
    plt.close(fig)

    # (b) Kim 2019: measured range [0, MDL] vs model
    fig, ax = plt.subplots(figsize=(6.4, 5.0))
    names = [r["name"] for r in kim]
    xi = np.arange(len(names))
    lo = np.array([r["lo"] for r in kim]); hi = np.array([r["hi"] for r in kim])
    ax.vlines(xi, lo, hi, color=OBS, lw=6, alpha=.25, label="measured, non-detects at 0 to MDL")
    ax.plot(xi, [r["obs_baf"] for r in kim], "o", color=OBS, ms=7, label="measured, published (ND = 0)")
    ax.plot(xi, [r["base"] for r in kim], "s", color=BASE, ms=7, label="model, default")
    ax.plot(xi, [r["prev"] for r in kim], "o", color=PREVC, ms=6, mfc="white", label="model, previous default")
    ax.set_xticks(xi); ax.set_xticklabels(names)
    ax.set_yscale("log")
    ax.set_ylabel("brown-rice BAF  [L/kg]")
    ax.set_title("Korean paddy fields (Kim 2019): grain BAF with its non-detects", fontsize=10.5)
    ax.legend(frameon=False, fontsize=8.5, loc="upper left")
    ax.grid(alpha=.25, which="both", axis="y")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.tight_layout()
    p2 = os.path.join(FIG_DIR, "policy_kim2019_grain.png")
    fig.savefig(p2, dpi=160)
    plt.close(fig)
    return p1, p2


def write_csv(path=OUT_CSV):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["section", "key", "value", "unit", "note"])
        for r in _ROWS:
            w.writerow(r)
    return path


def main(fast=False):
    rec("0_config", "shoot", api.DEFAULT_SHOOT, "-")
    rec("0_config", "lipid_loading", float(api.DEFAULT_LIPID_LOADING), "-")
    rec("0_config", "band_log10", api.APRIORI_LOG10_RMSE, "log10")
    cong = section1()
    section2()
    section2(cwo=0.05, tag="2_pfoa_50ngL")              # the deck's worked example
    kim, _ = section3()
    section4()
    section5()
    section6()
    p = write_csv()
    print()
    print("=" * 78)
    print(f"  기록: {os.path.relpath(p, ROOT_DIR)}  ({len(_ROWS)} 값)")
    if not fast:
        for f in figures(cong, kim):
            print(f"        {os.path.relpath(f, ROOT_DIR)}")
    print("  보고서: docs/POLICY_BRIEF_KR.md")
    return _ROWS


if __name__ == "__main__":
    main(fast="--fast" in sys.argv)
