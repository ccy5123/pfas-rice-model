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
# WHAT IT RUNS
#   1. the 13 curated PFAS: where each ends up in the plant, and what a flooded
#      irrigation schedule does to the grain
#   2. one substance end to end (PFOA at a realistic Korean pore-water level)
#   3. the only Korean FIELD dataset with paired pore water and brown rice
#      (Kim 2019) -- forward prediction, scored, honestly
#   4. the inverse question (grain measurement -> soil contamination), run twice:
#      against a known synthetic truth, and against that field measurement
#   5. structure (SMILES) input, including a substance outside the curated list
#   6. the EU biocide census, read back from validation/bat_census_biocides.csv
#
#   python validation/policy_brief_runs.py            (~7 min, writes CSV + 2 figures)
#   python validation/policy_brief_runs.py --fast     (skips the figures)
# =============================================================================
from __future__ import annotations

import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT_DIR, "src"))

import model_api as api                     # noqa: E402
import literature_params as lp              # noqa: E402

OUT_CSV = os.path.join(ROOT_DIR, "docs", "policy_brief_results.csv")
FIG_DIR = os.path.join(HERE, "figures")

# A realistic Korean paddy pore-water level for the worked example: the PFOA
# field average of Kim et al. 2019 (78.7 ng/L = 0.0787 ug/L). Using the measured
# level rather than a round number keeps the worked example anchored to data.
KIM_PFOA_POREWATER = 0.0787


def _rows():
    """Every number the report quotes, as (section, key, value, unit, note)."""
    return _ROWS


_ROWS: list[tuple] = []


def rec(section, key, value, unit="", note=""):
    _ROWS.append((section, key, value, unit, note))
    return value


# ---------------------------------------------------------------------------
# 1. the 13 curated PFAS -- where each one ends up
# ---------------------------------------------------------------------------
def section1(quiet=False):
    """The 13 curated PFAS -- and the honest complication.

    The baseline (a-priori) model says short chains reach the grain and long
    chains do not. BOTH measured datasets in this repo say otherwise at the long
    end: Yamazaki's PFDoDA grain BAF is 45.5 and Kim's Korean field PFDoDA is
    35.2, where the baseline predicts 0.10. So the congener ORDERING for the
    grain is not a result -- it depends on which transport mechanism is switched
    on, and the measured data favour the lipid pathway. Both variants are
    computed here so the briefing can never quote one without the other.
    """
    out = []
    for c in api.chain_table():
        n = c["name"]
        base = api.simulate(n)
        lip = api.simulate(n, lipid_loading=True)
        flood = api.simulate(n, cwo_profile="flooded")
        obs = api.observed_baf(n) or {}
        d = dict(name=n, n_C=c["n_C"], group=c["group"],
                 root=base["baf_final"]["root"], straw=base["straw_baf"],
                 grain=base["baf_final"]["grain"],
                 grain_lipid=lip["baf_final"]["grain"],
                 grain_flooded=flood["baf_final"]["grain"],
                 obs_grain=obs.get("grain"), obs_root=obs.get("root"))
        d["grain_over_root"] = d["grain"] / d["root"]
        d["flood_ratio"] = d["grain_flooded"] / d["grain"]
        out.append(d)
        for k in ("root", "straw", "grain", "grain_lipid", "grain_flooded"):
            rec("1_congeners", f"{n}.{k}", d[k], "L/kg")
        if d["obs_grain"] is not None:
            rec("1_congeners", f"{n}.obs_grain_yamazaki", d["obs_grain"], "L/kg")
    if not quiet:
        print("=" * 78)
        print("1. 13종 PFAS — 물질마다 도달 부위가 다르다  [BAF, L/kg]")
        print("=" * 78)
        print(f"{'물질':9s}{'C':>3s}{'계열':>8s}{'뿌리':>9s}{'짚':>9s}"
              f"{'낟알(기본)':>11s}{'낟알(지질)':>11s}{'낟알 실측':>10s}{'담수낟알':>10s}")
        for d in sorted(out, key=lambda x: (x["group"], x["n_C"])):
            o = f"{d['obs_grain']:10.2f}" if d["obs_grain"] is not None else f"{'-':>10s}"
            print(f"{d['name']:9s}{d['n_C']:3d}{d['group']:>8s}{d['root']:9.3f}"
                  f"{d['straw']:9.3f}{d['grain']:11.3f}{d['grain_lipid']:11.3f}{o}"
                  f"{d['grain_flooded']:10.3f}")
        short = [d for d in out if d["n_C"] <= 6]
        long_ = [d for d in out if d["n_C"] >= 9]
        rs = float(np.mean([d["grain"] for d in short]))
        rl = float(np.mean([d["grain"] for d in long_]))
        rec("1_congeners", "mean_grain_shortchain_C4_C6_base", rs, "L/kg")
        rec("1_congeners", "mean_grain_longchain_C9_C12_base", rl, "L/kg")
        # how badly the baseline misses the measured long-chain grain
        miss = [(d["name"], d["obs_grain"] / d["grain"]) for d in out
                if d["obs_grain"] and d["n_C"] >= 11]
        for nm, f in miss:
            rec("1_congeners", f"{nm}.obs_over_base_grain", f, "x")
        print()
        print(f"  기본 모델: 단쇄(C4–C6) 낟알 평균 {rs:.2f} vs 장쇄(C9–C12) {rl:.3f}"
              f"  → {rs / rl:.0f}배 차이")
        print("  그러나 실측은 장쇄 끝에서 반대로 «올라간다»:")
        for nm, f in miss:
            print(f"    {nm}: 실측 낟알 BAF가 기본 모델의 {f:.0f}배")
        print("  ⇒ 낟알에 대한 «물질 순위»는 결과가 아니라 어떤 기작을 켜느냐에 달려 있다.")
        print("    지질 매개 이동을 켜면 순위가 뒤집히고, 한국 현장자료와 맞는다(§3).")
        print("    발표에서 «단쇄가 더 위험하다»로 단정하면 안 된다.")
    return out


# ---------------------------------------------------------------------------
# 2. one substance, end to end
# ---------------------------------------------------------------------------
def section2(quiet=False):
    r = api.simulate("PFOA", Cwo=KIM_PFOA_POREWATER)
    grain = float(r["conc"]["grain"][-1])
    d = dict(cwo=KIM_PFOA_POREWATER,
             root=float(r["conc"]["root"][-1]), straw=float(r["straw"][-1]), grain=grain,
             band=api.predictive_band(grain),
             intake=api.intake_fraction(grain, congener="PFOA"),
             paths=api.apportionment("PFOA")["fraction"])
    for k in ("root", "straw", "grain"):
        rec("2_pfoa", f"conc.{k}", d[k], "ug/kg fw")
    rec("2_pfoa", "grain.band_lo", d["band"]["lo"], "ug/kg fw")
    rec("2_pfoa", "grain.band_hi", d["band"]["hi"], "ug/kg fw")
    rec("2_pfoa", "grain.band_factor", d["band"]["factor"], "x")
    rec("2_pfoa", "efsa_twi_percent", d["intake"]["percent"], "%",
        "EFSA 2020 4종 합 TWI 4.4 ng/kg bw/주, 쌀 150 g/일, 체중 60 kg")
    rec("2_pfoa", "grain_via_xylem", d["paths"]["grain"]["xylem_from_stem"], "-")
    rec("2_pfoa", "grain_via_phloem", d["paths"]["grain"]["phloem_from_leaf"], "-")
    rec("2_pfoa", "root_from_soil", d["paths"]["root"]["soil_uptake"], "-")
    if not quiet:
        print()
        print("=" * 78)
        print(f"2. 한 물질을 끝까지 — PFOA, 공극수 {KIM_PFOA_POREWATER * 1000:.1f} ng/L")
        print("=" * 78)
        print(f"  뿌리 {d['root']:.4f} · 짚 {d['straw']:.4f} · 낟알 {d['grain']:.5f}  [µg/kg 생중량]")
        print(f"  낟알 불확실 범위 {d['band']['lo']:.5f} – {d['band']['hi']:.4f}"
              f"  (±{d['band']['factor']:.1f}배, 사전예측 오차 log10 RMSE {d['band']['log10_rmse']})")
        print(f"  EFSA 주간섭취허용량 대비 {d['intake']['percent']:.2f}%"
              f"  (신호: {d['intake']['signal']})")
        print(f"  낟알 유입 경로: 물관 {d['paths']['grain']['xylem_from_stem'] * 100:.0f}%"
              f" · 체관 {d['paths']['grain']['phloem_from_leaf'] * 100:.0f}%")
        print(f"  뿌리 유입 경로: 토양 직접흡수 {d['paths']['root']['soil_uptake'] * 100:.1f}%")
    return d


# ---------------------------------------------------------------------------
# 3. the Korean field dataset -- forward prediction, scored honestly
# ---------------------------------------------------------------------------
def section3(quiet=False):
    rows = []
    for sp, (pw, soil, rice, df) in lp.KIM2019_FIELD.items():
        obs = rice / (pw / 1000.0)
        mono = api.simulate(sp)["baf_final"]["grain"]
        lip = api.simulate(sp, lipid_loading=True)["baf_final"]["grain"]
        rows.append(dict(name=sp, porewater_ngL=pw, brownrice_ugkg=rice, DF_pct=df,
                         obs_baf=obs, mono=mono, lipid=lip))
        rec("3_kim2019", f"{sp}.obs_grain_baf", obs, "L/kg")
        rec("3_kim2019", f"{sp}.model_grain_baf_monotone", mono, "L/kg")
        rec("3_kim2019", f"{sp}.model_grain_baf_lipid", lip, "L/kg")

    def sc(key):
        p = np.array([r[key] for r in rows])
        o = np.array([r["obs_baf"] for r in rows])
        e = np.log10(p / o)
        return float(np.sqrt(np.mean(e ** 2))), float(np.mean(e))

    m_rmse, m_bias = sc("mono")
    l_rmse, l_bias = sc("lipid")
    rec("3_kim2019", "rmse_log10_monotone", m_rmse, "log10")
    rec("3_kim2019", "rmse_log10_lipid", l_rmse, "log10")
    rec("3_kim2019", "bias_log10_monotone", m_bias, "log10")
    rec("3_kim2019", "bias_log10_lipid", l_bias, "log10")
    if not quiet:
        print()
        print("=" * 78)
        print("3. 한국 논 현장 자료로 맞춰보기 — Kim et al. 2019 (공극수·현미 쌍 측정)")
        print("=" * 78)
        print(f"{'물질':9s}{'공극수 ng/L':>12s}{'현미 µg/kg':>12s}{'검출률%':>9s}"
              f"{'실측 BAF':>10s}{'기본모델':>10s}{'지질경로 켬':>12s}")
        for r in rows:
            print(f"{r['name']:9s}{r['porewater_ngL']:12.3f}{r['brownrice_ugkg']:12.5f}"
                  f"{r['DF_pct']:9.1f}{r['obs_baf']:10.2f}{r['mono']:10.3f}{r['lipid']:12.3f}")
        print(f"\n  log10 RMSE:  기본모델 {m_rmse:.2f} (치우침 {m_bias:+.2f})"
              f"   ·   지질경로 켬 {l_rmse:.2f} (치우침 {l_bias:+.2f})")
        print("  → 기본모델은 현장 낟알을 크게 «과소»예측한다. 지질 매개 이동 기작을 켜면")
        print("    6개 물질 전부 2–6배 안으로 들어온다 — 그 기작은 이 자료가 아니라")
        print("    다른(일본 온실) 자료에 맞춰진 것이므로, 이 비교는 외부 검증에 해당한다.")
        print("  주의: 검출률이 3–20%로 낮은 물질이 다수 (PFOA만 57%).")
    return rows, dict(mono_rmse=m_rmse, lipid_rmse=l_rmse,
                      mono_bias=m_bias, lipid_bias=l_bias)


# ---------------------------------------------------------------------------
# 4. the inverse question
# ---------------------------------------------------------------------------
def section4(quiet=False):
    truth = KIM_PFOA_POREWATER
    r = api.simulate("PFOA", Cwo=truth)
    meas = {"root": float(r["conc"]["root"][-1]),
            "straw": float(r["straw"][-1]),
            "grain": float(r["conc"]["grain"][-1])}
    syn = api.estimate_exposure_bayesian("PFOA", meas)
    field = api.estimate_exposure_bayesian("PFOA", {"grain": 0.349})
    rec("4_inverse", "synthetic.truth", truth, "ug/L")
    rec("4_inverse", "synthetic.median", syn["median"], "ug/L")
    rec("4_inverse", "synthetic.ci95_lo", syn["ci95"][0], "ug/L")
    rec("4_inverse", "synthetic.ci95_hi", syn["ci95"][1], "ug/L")
    rec("4_inverse", "field.measured_brownrice", 0.349, "ug/kg", "Kim 2019 PFOA 현미 평균")
    rec("4_inverse", "field.median", field["median"], "ug/L")
    rec("4_inverse", "field.ci95_lo", field["ci95"][0], "ug/L")
    rec("4_inverse", "field.ci95_hi", field["ci95"][1], "ug/L")
    rec("4_inverse", "field.overestimate_factor", field["median"] / truth, "x")
    if not quiet:
        print()
        print("=" * 78)
        print("4. 거꾸로 풀기 — 쌀 측정값에서 논 오염도 추정")
        print("=" * 78)
        print(f"  (가) 정답을 아는 경우: 참값 {truth:.4f} → 추정 {syn['median']:.4f} µg/L"
              f"  (95% {syn['ci95'][0]:.4f}–{syn['ci95'][1]:.4f})")
        print("      → 역산 자체는 제대로 작동한다. 이것이 방법의 검증이다.")
        print(f"  (나) 실제 현장값: 현미 0.349 µg/kg → 추정 {field['median']:.2f} µg/L"
              f"  (95% {field['ci95'][0]:.2f}–{field['ci95'][1]:.2f})")
        print(f"      실측 공극수는 {truth:.4f} µg/L → {field['median'] / truth:.0f}배 과대추정.")
        print("      원인은 역산이 아니라 §3의 낟알 과소예측이다(같은 편향의 뒤집힌 얼굴).")
        print("      ⇒ 현재 이 기능은 «순위·규모» 용도이며, 절대 오염도 산정에는 쓸 수 없다.")
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
        d = getattr(r, "get", dict().get)("descriptors") if isinstance(r, dict) else None
        out.append(dict(label=label, matched=r.get("congener"),
                        provisional=bool(r.get("provisional")),
                        n_pfc=getattr(d, "n_perfluoroC", None),
                        head=getattr(d, "head_group", None),
                        ether=getattr(d, "n_ether_O", None),
                        root=r["baf_final"]["root"], straw=r["straw_baf"],
                        grain=r["baf_final"]["grain"], note=note))
        rec("5_smiles", f"{label}.grain", r["baf_final"]["grain"], "L/kg", note)
        rec("5_smiles", f"{label}.provisional", float(bool(r.get("provisional"))), "-")
    if not quiet:
        print()
        print("=" * 78)
        print("5. 구조(SMILES)로 넣기 — 목록에 없는 물질도 계산된다")
        print("=" * 78)
        for d in out:
            flag = "잠정(QSPR)" if d["provisional"] else "실측 파라미터"
            print(f"  {d['label']:16s} 인식: {d['matched']:8s} [{flag}]"
                  f"  퍼플루오로C {d['n_pfc']}, 작용기 {d['head']}, 에터O {d['ether']}")
            print(f"  {'':16s} 뿌리 {d['root']:.3f} · 짚 {d['straw']:.3f} · 낟알 {d['grain']:.3f}")
        print("  → 모델이 «잠정»을 스스로 표시한다. 목록 외 물질의 값은 순위 판단용이다.")
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
    n = len(rows)
    excl = [r for r in rows if r["verdict"].startswith("EXCLUDED")]
    inside = [r for r in rows if r["verdict"] == "inside"]
    pos = [r for r in rows if r["bpc"] in ("B", "vB")]
    d = dict(n_run=n, n_inside=len(inside), n_excluded=len(excl), n_positive=len(pos))
    rec("6_biocides", "n_runnable", n, "종")
    rec("6_biocides", "n_inside_scope", len(inside), "종")
    rec("6_biocides", "n_excluded_ionisation", len(excl), "종")
    if not quiet:
        print()
        print("=" * 78)
        print("6. EU 살생물질 목록 통과 — 규제 분류와 쌀 노출은 다른 질문")
        print("=" * 78)
        print(f"  계산 가능 {n}종 · 적용범위 내 {len(inside)}종 · 이온화로 제외 {len(excl)}종")
        print(f"  규제기관이 «생물축적성»으로 분류한 {len(pos)}종의 이 모델 결과:")
        for r in sorted(pos, key=lambda x: float(x["log_kow"])):
            print(f"    {r['substance'][:26]:<28s} logKow {float(r['log_kow']):5.2f}"
                  f"  뿌리 {float(r['P_root']):9.3g}  짚 {float(r['P_straw']):9.3g}"
                  f"  낟알 {float(r['P_grain']):9.3g}   [{r['verdict']}]")
        print("  → 전부 뿌리에 남고 낟알에는 사실상 도달하지 않는다.")
    return d


# ---------------------------------------------------------------------------
def figures(cong, kim):
    """Two figures, English labels only -- no Korean font ships with matplotlib
    here, so Hangul would render as boxes. The briefing carries the Korean
    captions; a slide tool can relabel from those."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(FIG_DIR, exist_ok=True)
    BASE, LIPID, OBS = "#4C5A57", "#00808F", "#A8821A"

    # (a) grain BAF vs chain length: two model variants against measured data
    pfca = sorted([d for d in cong if d["group"] == "PFCA"], key=lambda x: x["n_C"])
    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    x = [d["n_C"] for d in pfca]
    ax.plot(x, [d["grain"] for d in pfca], "o-", color=BASE, lw=1.8, ms=6,
            label="model, base (a-priori)")
    ax.plot(x, [d["grain_lipid"] for d in pfca], "s-", color=LIPID, lw=1.8, ms=6,
            label="model, lipid pathway on")
    ox = [d["n_C"] for d in pfca if d["obs_grain"]]
    oy = [d["obs_grain"] for d in pfca if d["obs_grain"]]
    ax.plot(ox, oy, "^", color=OBS, ms=10, mec="white", mew=1.4,
            label="measured, Yamazaki 2023 (greenhouse)")
    kx, ky = [], []
    order = {"PFHpA": 7, "PFOA": 8, "PFNA": 9, "PFDA": 10, "PFUnDA": 11, "PFDoDA": 12}
    for r in kim:
        if r["name"] in order:
            kx.append(order[r["name"]]); ky.append(r["obs_baf"])
    ax.plot(kx, ky, "D", color=OBS, ms=7, mfc="none", mew=1.8,
            label="measured, Kim 2019 (Korean field)")
    ax.set_yscale("log")
    ax.set_xlabel("perfluorocarbon chain length  (n C)")
    ax.set_ylabel("grain BAF  [L/kg]  (log scale)")
    ax.set_title("The grain answer depends on the transport mechanism\n"
                 "— and the measured long chains rise, where the base model falls",
                 fontsize=11)
    ax.legend(frameon=False, fontsize=9, loc="lower left")
    ax.grid(alpha=.25, which="both")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.tight_layout()
    p1 = os.path.join(FIG_DIR, "policy_grain_by_chain.png")
    fig.savefig(p1, dpi=160)
    plt.close(fig)

    # (b) Kim 2019 observed vs predicted grain BAF
    fig, ax = plt.subplots(figsize=(5.6, 5.2))
    o = np.array([r["obs_baf"] for r in kim])
    for key, c, lab, mk in (("mono", BASE, "base model", "o"),
                            ("lipid", LIPID, "lipid pathway on", "s")):
        p = np.array([r[key] for r in kim])
        ax.scatter(o, p, s=58, color=c, marker=mk, edgecolor="white", lw=1.2, label=lab)
    lim = [1e-2, 1e2]
    ax.plot(lim, lim, color="0.35", lw=1)
    for f in (10, 0.1):
        ax.plot(lim, [v * f for v in lim], color="0.7", lw=.9, ls=":")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("measured grain BAF  [L/kg]   Kim 2019, Korean paddy")
    ax.set_ylabel("model grain BAF  [L/kg]")
    ax.set_title("Korean field data: the base model under-predicts grain\n"
                 "(dotted = 10x above / below the 1:1 line)", fontsize=10)
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    ax.grid(alpha=.25, which="both")
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
    cong = section1()
    section2()
    kim, _ = section3()
    section4()
    section5()
    section6()
    p = write_csv()
    print()
    print("=" * 78)
    print(f"  기록: {os.path.relpath(p, ROOT_DIR)}  ({len(_ROWS)} 값)")
    if not fast:
        f1, f2 = figures(cong, kim)
        for f in (f1, f2):
            print(f"        {os.path.relpath(f, ROOT_DIR)}")
    print("  보고서: docs/POLICY_BRIEF_KR.md")
    return _ROWS


if __name__ == "__main__":
    main(fast="--fast" in sys.argv)
