#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""정책활용협의회 브리핑 — EU 살생물질 census + 적용범위(AD) 그림/수치.

`docs/POLICY_BRIEF_BIOCIDES_KR.md` 가 인용하는 모든 수치와 그림을 만든다.
원본 계산은 `validation/bat_census_biocides.py` 이고, 이 파일은 그 결과를
**정책 발표용으로 다시 읽어 그리는** 층이다 — 모델을 다시 정의하지 않는다.

    python validation/policy_biocides_figures.py           (~4분, 그림 4장 + CSV)
    python validation/policy_biocides_figures.py --fast    (log Kow 스윕 생략, 그림 3장)

주의 — 이 파일이 새로 계산하는 것은 **log Kow 스윕 하나뿐**이다(그림 4).
나머지는 전부 `validation/bat_census_biocides.csv` 를 읽는다. 그 CSV 가 없으면
먼저 `python validation/bat_census_biocides.py` 를 돌려야 한다. 그렇게 나눈 이유는
발표용 그림을 고치려고 117행 census 를 다시 돌리는 일이 없게 하기 위해서다.

그림 라벨은 **영어**다. 기존 정책 그림(`policy_brief_runs.figures`)과 같은 규칙이고,
한글 폰트가 없는 환경에서도 두부(□□□) 없이 같은 PNG 가 나오게 하기 위해서다.
슬라이드를 만들 때 한글로 다시 그려도 되고, 그 경우 `_labels()` 의 문자열만 바꾸면 된다.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(HERE)
FIG_DIR = os.path.join(HERE, "figures")
CENSUS_CSV = os.path.join(HERE, "bat_census_biocides.csv")
OUT_CSV = os.path.join(ROOT_DIR, "docs", "policy_biocides_results.csv")
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT_DIR, "src"))

# 판정 4등급 — bat_census_biocides.screen() 이 붙이는 문자열 그대로
V_IN, V_ROOT = "inside", "root only"
V_EXT, V_EXCL = "extrapolation", "EXCLUDED (ionisation)"
ORDER = [V_IN, V_ROOT, V_EXT, V_EXCL]
COLOR = {V_IN: "#2b7a3d", V_ROOT: "#c98200", V_EXT: "#8a8a8a", V_EXCL: "#b02020"}


def _plain_log(ax, axis="both"):
    """로그 축 눈금을 mathtext 없이 ASCII 로 찍는다.

    matplotlib 의 기본 로그 포매터는 `$10^{-9}$` 를 mathtext 로 그리고 mathtext 의
    마이너스는 U+2212 인데, 한글을 가진 CJK 폰트에는 그 글리프가 없어 두부가 된다.
    폰트를 바꾸는 대신 포매터를 바꾸는 쪽이 어느 폰트에서도 안전하다."""
    from matplotlib.ticker import FuncFormatter, LogLocator

    def fmt(v, _pos):
        if v <= 0:
            return ""
        e = int(round(np.log10(v)))
        if abs(v - 10.0 ** e) > 1e-9 * max(1.0, v):
            return ""
        if -2 <= e <= 3:
            return f"{10.0 ** e:g}"
        return f"1e{e}"

    for a in ((ax.xaxis, ax.yaxis) if axis == "both" else
              ((ax.xaxis,) if axis == "x" else (ax.yaxis,))):
        a.set_major_locator(LogLocator(base=10.0))
        a.set_major_formatter(FuncFormatter(fmt))
        a.set_minor_formatter(FuncFormatter(lambda *_: ""))


def _labels():
    """그림 문자열 — 영어. 한글로 다시 그리려면 여기만 바꾸면 된다."""
    return dict(
        logkow="log Kow  (more lipophilic →)",
        ad_title="The applicability domain is two nested bands, not one",
        ad_sub="117 EU biocide rows placed on this repo's measured-data spans",
        band_root="root partition — 560 measured rows",
        band_shoot="shoot (TSCF) — 30 measured rows",
        v={V_IN: "inside", V_ROOT: "root only", V_EXT: "beyond both anchors",
           V_EXCL: "refused: too ionised"},
        n_sub="substances",
        tscf_title="Why the shoot domain is the narrower one",
        tscf_sub="root->shoot loading (TSCF) peaks at log Kow 1.78 and collapses",
        tscf_y="TSCF  (root -> xylem loading)",
        tscf_peak="peak log Kow 1.78",
        tscf_edge="shoot data edge 5.46",
        bioacc_ok="bioaccumulative in fish, and this model answers",
        bioacc_no="bioaccumulative in fish, but refused here (too ionised)",
        rank_title="A fish BCF says nothing about the grain",
        rank_sub="same 60 substances, same log Kow -- with the root, against the straw",
        bat_x="BAT fish BCF  (Scenario A, L/kg)",
        root_y="this model, root BAF  (L/kg)",
        straw_y="this model, straw BAF  (L/kg)",
        rho="Spearman rho",
        ratio_title="Above log Kow ~6.8 the root number changes meaning",
        ratio_sub="one season (120 d) cannot equilibrate the root -- past the red line "
                  "the root BAF is a rate, not a partition coefficient",
        ratio_y="root BAF / K_PW  (1 = equilibrated)",
        ratio_eq="equilibrium (safe to read as a partition coefficient)",
        nominal="nominal AD edge 8.70  (measured data reaches here)",
        effective="effective edge {edge:.2f}  (root BAF = half of K_PW)",
    )


# ---------------------------------------------------------------------------
def load_census():
    """bat_census_biocides.csv 를 읽는다. 이 파일은 census 를 다시 돌리지 않는다."""
    if not os.path.exists(CENSUS_CSV):
        raise SystemExit(
            "validation/bat_census_biocides.csv 가 없습니다.\n"
            "먼저: python validation/bat_census_biocides.py")
    rows = []
    with open(CENSUS_CSV, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            def f(k):
                v = (r.get(k) or "").strip()
                try:
                    return float(v)
                except ValueError:
                    return None
            rows.append(dict(substance=r["substance"], verdict=r["verdict"],
                             log_kow=f("log_kow"), bpc=(r.get("bpc") or "").strip(),
                             bat_A=f("bat_A"), tscf=f("P_TSCF"), kpw=f("P_K_PW_root"),
                             root=f("P_root"), straw=f("P_straw"), grain=f("P_grain")))
    return rows


def spans():
    import bat_census_biocides as BC
    return BC.measured_spans()


def _rho(x, y):
    """스피어만 — census 스크립트와 같은 정의를 쓰되 여기서 다시 계산하지 않고
    순위 상관을 직접 구한다(의존성 최소화)."""
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    rx, ry = rank(x), rank(y)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5
    return num / den if den else float("nan")


def rank_pairs(rows):
    """순위상관에 들어가는 짝. census 스크립트와 **같은 배제 규칙**을 쓴다:
    어류 BCF 가 없는 행(수집시트 52종 포함), 해석 불가 판정을 받은 행, 그리고
    한 물질을 두 번 세게 만드는 «두 번째 log Kow» 행."""
    out = []
    for r in rows:
        if r["bat_A"] is None:
            continue
        if "alt log Kow" in r["substance"] or "AR measured" in r["substance"]:
            continue
        if "diamine" in r["substance"].lower():   # BAT 출력이 해석 불가로 표시된 행
            continue
        out.append(r)
    return out


# ---------------------------------------------------------------------------
def fig_domain(rows, sp, L, path):
    """그림 1 — AD 지도. 두 겹 밴드 위에 117행을 4등급 색으로 얹는다."""
    import matplotlib.pyplot as plt
    r_lo, r_hi, r_n = sp["partition"]
    t_lo, t_hi, t_n = sp["tscf"]
    fig, (ax, axh) = plt.subplots(2, 1, figsize=(11.0, 7.2), sharex=True,
                                  gridspec_kw=dict(height_ratios=[1.25, 1.0],
                                                   hspace=0.16))
    # 위 패널 — 두 겹 밴드 + 등급별 rug (각 등급이 자기 줄을 갖는다)
    ax.barh(4.0, r_hi - r_lo, left=r_lo, height=0.62, color="#cfe3d4",
            edgecolor="#2b7a3d", lw=1.4)
    ax.barh(3.1, t_hi - t_lo, left=t_lo, height=0.62, color="#dfe6f5",
            edgecolor="#33529e", lw=1.4)
    ax.text(r_lo + 0.15, 4.0, L["band_root"], va="center", fontsize=9.5,
            color="#1d5c2c")
    ax.text(t_lo + 0.15, 3.1, L["band_shoot"], va="center", fontsize=9.5,
            color="#2a3f78")
    ax.text(r_hi + 0.15, 4.0, f"{r_hi:.2f}", va="center", fontsize=9.5,
            color="#2b7a3d", fontweight="bold")
    ax.text(t_hi + 0.15, 3.1, f"{t_hi:.2f}", va="center", fontsize=9.5,
            color="#33529e", fontweight="bold")
    rug_y = {V_IN: 2.1, V_ROOT: 1.5, V_EXT: 0.9, V_EXCL: 0.3}
    for v in ORDER:
        xs = [r["log_kow"] for r in rows if r["verdict"] == v]
        ax.plot(xs, [rug_y[v]] * len(xs), "|", ms=13, mew=1.6, color=COLOR[v])
    ax.axvline(t_hi, color="#33529e", ls="--", lw=1.0, alpha=0.6, zorder=0)
    ax.axvline(r_hi, color="#2b7a3d", ls="--", lw=1.0, alpha=0.6, zorder=0)
    ax.set_ylim(-0.15, 4.6)
    # 등급 이름은 축 밖(눈금 라벨)에 둔다 -- 축 안에 쓰면 왼쪽 끝 물질과 겹친다
    ax.set_yticks([rug_y[v] for v in ORDER])
    ax.set_yticklabels([f"{L['v'][v]}  n="
                        f"{len([r for r in rows if r['verdict'] == v])}"
                        for v in ORDER], fontsize=9.5)
    for lab, v in zip(ax.get_yticklabels(), ORDER):
        lab.set_color(COLOR[v])
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)

    # 아래 패널 — 등급별 누적 히스토그램
    bins = np.arange(-3.5, 11.5, 0.5)
    bottom = np.zeros(len(bins) - 1)
    for v in ORDER:
        xs = [r["log_kow"] for r in rows if r["verdict"] == v]
        h, _ = np.histogram(xs, bins=bins)
        axh.bar(bins[:-1], h, width=0.48, bottom=bottom, align="edge",
                color=COLOR[v])
        bottom += h
    axh.axvline(t_hi, color="#33529e", ls="--", lw=1.2)
    axh.axvline(r_hi, color="#2b7a3d", ls="--", lw=1.2)
    axh.set_xlabel(L["logkow"], fontsize=10.5)
    axh.set_ylabel(L["n_sub"], fontsize=10.5)
    axh.set_ylim(0, bottom.max() * 1.45)
    # 범례 없음 -- 위 패널의 y 눈금이 이미 같은 4등급을 같은 색으로 이름 붙였고,
    # 두 패널 사이에 범례를 끼우면 위 패널의 rug 와 겹친다.
    axh.set_ylim(0, bottom.max() * 1.08)
    axh.set_xlim(-3.5, 11.3)
    for a in (ax, axh):
        a.spines[["top", "right"]].set_visible(False)
    fig.suptitle(L["ad_title"], fontsize=14, x=0.055, ha="left", y=0.985)
    fig.text(0.055, 0.943, L["ad_sub"], fontsize=10, color="#555")
    fig.subplots_adjust(top=0.90)
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    return path


def fig_tscf(rows, sp, L, path):
    """그림 2 — TSCF 종형곡선. 지상부 AD 가 왜 좁은지를 한 장으로."""
    import matplotlib.pyplot as plt
    import neutral_dpu as ND
    t_hi = sp["tscf"][1]
    x_hi = 8.8
    g = np.linspace(-3.0, x_hi, 400)
    tscf = np.array([ND.tscf(float(x)) for x in g])   # 기본 = Briggs 종형
    fig, ax = plt.subplots(figsize=(10.0, 5.4))
    ax.plot(g, tscf, lw=2.2, color="#33529e", zorder=3)
    ax.axvline(1.78, color="#c98200", ls="--", lw=1.3)
    ax.text(1.72, 3e-9, L["tscf_peak"], rotation=90, ha="right", va="bottom",
            fontsize=9, color="#c98200")
    ax.axvline(t_hi, color="#888", ls=":", lw=1.4)
    ax.text(t_hi - 0.08, 3e-9, L["tscf_edge"], rotation=90, ha="right",
            va="bottom", fontsize=9, color="#666")
    # 7종을 «모델이 답한 것»과 «이온화로 거부한 것»으로 나눠 찍는다. 후자의 TSCF 는
    # 곡선 위치를 보이려고 그린 것이지 인용해도 되는 값이 아니다.
    bio = [r for r in rows if r["bpc"] in ("B", "vB") and r["tscf"]]
    ok = [r for r in bio if not r["verdict"].startswith("EXCLUDED")]
    no = [r for r in bio if r["verdict"].startswith("EXCLUDED")]
    if ok:
        ax.plot([r["log_kow"] for r in ok], [r["tscf"] for r in ok], "o", ms=9,
                color="#b02020", zorder=6, label=f"{L['bioacc_ok']}  (n={len(ok)})")
    if no:
        ax.plot([r["log_kow"] for r in no], [r["tscf"] for r in no], "o", ms=9,
                mfc="none", mec="#b02020", mew=1.8, zorder=6,
                label=f"{L['bioacc_no']}  (n={len(no)})")
    # 축 안에는 곡선과 두 개의 세로 주석이 이미 있어 범례가 들어갈 자리가 없다.
    # 축 바깥(제목 아래)에 둔다 -- bbox_inches="tight" 가 여백을 늘려 준다.
    ax.legend(frameon=False, fontsize=9.5, loc="lower left",
              bbox_to_anchor=(0.0, 1.01), ncol=1, handletextpad=0.6)
    ax.set_yscale("log")
    _plain_log(ax, "y")
    ax.set_ylim(1e-9, 3.0)
    ax.set_xlim(-3.0, x_hi)
    ax.set_xlabel(L["logkow"], fontsize=10.5)
    ax.set_ylabel(L["tscf_y"], fontsize=10.5)
    ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle(L["tscf_title"], fontsize=14, x=0.055, ha="left", y=1.10)
    fig.text(0.055, 1.045, L["tscf_sub"], fontsize=10, color="#555")
    fig.subplots_adjust(top=0.88)
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    return path


def fig_rank(pairs, L, path):
    """그림 3 — 어류 BCF vs 벼 뿌리 / 벼 짚, 같은 60종."""
    import matplotlib.pyplot as plt
    x = [p["bat_A"] for p in pairs]
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.9))
    for ax, key, ylab, col in ((axes[0], "root", L["root_y"], "#2b7a3d"),
                               (axes[1], "straw", L["straw_y"], "#b02020")):
        y = [p[key] for p in pairs]
        rho = _rho(x, y)
        ax.plot(x, y, "o", ms=6, alpha=0.75, color=col)
        ax.set_xscale("log")
        ax.set_yscale("log")
        _plain_log(ax, "both")
        ax.set_xlabel(L["bat_x"], fontsize=10)
        ax.set_ylabel(ylab, fontsize=10)
        ax.set_title(f"{L['rho']} = {rho:+.3f}   (n={len(pairs)})", fontsize=11.5)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle(L["rank_title"], fontsize=13, x=0.02, ha="left")
    fig.text(0.02, 0.925, L["rank_sub"], fontsize=9.5, color="#555")
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    return path


def fig_ratio(sweep, sp, L, path):
    """그림 4 — 뿌리 BAF / K_PW. 명목 AD 와 실질 AD 가 갈리는 지점."""
    import matplotlib.pyplot as plt
    r_hi = sp["partition"][1]
    g, ratio = sweep["grid"], sweep["ratio"]
    edge = half_meaning_edge(sweep)
    fig, ax = plt.subplots(figsize=(10.0, 5.4))
    ax.axvspan(edge, g[-1], color="#f6e2e2", alpha=0.8, zorder=0)
    ax.plot(g, ratio, lw=2.4, color="#2b7a3d", zorder=3)
    ax.axhline(1.0, color="#999", ls=":", lw=1.2)
    ax.text(g[0] + 0.1, 1.015, L["ratio_eq"], fontsize=9, color="#666", va="bottom")
    # 두 경계는 세로 라벨로 -- 곡선/그늘과 겹치지 않는 유일한 배치
    ax.axvline(edge, color="#b02020", lw=1.8, zorder=4)
    ax.text(edge - 0.10, 0.03, L["effective"].format(edge=edge), rotation=90,
            ha="right", va="bottom", fontsize=9.5, color="#b02020")
    ax.axvline(r_hi, color="#2b7a3d", ls="--", lw=1.4, zorder=4)
    ax.text(r_hi - 0.10, 0.03, L["nominal"], rotation=90, ha="right",
            va="bottom", fontsize=9.5, color="#2b7a3d")
    ax.set_ylim(-0.02, 1.10)
    ax.set_xlim(g[0], g[-1])
    ax.set_xlabel(L["logkow"], fontsize=10.5)
    ax.set_ylabel(L["ratio_y"], fontsize=10.5)
    ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle(L["ratio_title"], fontsize=14, x=0.055, ha="left", y=0.995)
    fig.text(0.055, 0.925, L["ratio_sub"], fontsize=10, color="#555")
    fig.subplots_adjust(top=0.85)
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
def run_sweep(lo=0.0, hi=10.5, n=36, quiet=False):
    """그림 4 가 필요로 하는 유일한 새 계산: log Kow 격자 위의 뿌리 BAF / K_PW.

    census 의 §5(d) 와 같은 양을 같은 함수(`run_one`)로 구한다 — 값이 갈라질 수
    없도록 재구현하지 않는다."""
    import bat_census_biocides as BC
    grid = np.linspace(lo, hi, n)
    root, kpw = [], []
    for i, lk in enumerate(grid):
        r = BC.run_one(float(lk), f"sweep{i}")
        root.append(r["root"])
        kpw.append(r["K_PW_root"])
        if not quiet:
            print(f"    log Kow {lk:5.2f}  K_PW {r['K_PW_root']:11.4g}  "
                  f"root {r['root']:11.4g}  "
                  f"ratio {r['root'] / r['K_PW_root']:7.4f}")
    root, kpw = np.array(root), np.array(kpw)
    return dict(grid=grid, root=root, kpw=kpw, ratio=root / kpw)


def half_meaning_edge(sweep):
    """뿌리 BAF 가 K_PW 의 절반으로 떨어지는 log Kow — «분배계수로 읽어도 되는»
    구간의 실질 상한. 보고서가 인용하는 값이라 여기서 한 번만 정의한다."""
    g, ratio = sweep["grid"], sweep["ratio"]
    for i in range(1, len(g)):
        if ratio[i] < 0.5 <= ratio[i - 1]:
            f = (ratio[i - 1] - 0.5) / (ratio[i - 1] - ratio[i])
            return float(g[i - 1] + f * (g[i] - g[i - 1]))
    return float("nan")


# ---------------------------------------------------------------------------
def _carry_over(path, section):
    """기존 CSV 에서 한 section 을 그대로 읽어 온다 (--fast 가 지우지 않도록)."""
    if not os.path.exists(path):
        return []
    out = []
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r.get("section") == section:
                note = r.get("note", "")
                if "이전 실행에서 옮김" not in note:
                    note = (note + " [--fast: 이전 실행에서 옮김]").strip()
                out.append(dict(section=r["section"], key=r["key"],
                                value=r["value"], unit=r.get("unit", ""),
                                note=note))
    return out


def write_csv(rows, sp, pairs, sweep, path=OUT_CSV):
    recs = []

    def rec(section, key, value, unit="", note=""):
        recs.append(dict(section=section, key=key, value=value, unit=unit, note=note))

    r_lo, r_hi, r_n = sp["partition"]
    t_lo, t_hi, t_n = sp["tscf"]
    rec("1_ad", "root_span_lo", r_lo, "log Kow", "측정 뿌리분배 자료의 하한")
    rec("1_ad", "root_span_hi", r_hi, "log Kow", "측정 뿌리분배 자료의 상한")
    rec("1_ad", "root_span_n", r_n, "행", "data_obs/ 중성 테이블 6종 합계")
    rec("1_ad", "shoot_span_lo", t_lo, "log Kow", "비이온화 TSCF 측정 하한")
    rec("1_ad", "shoot_span_hi", t_hi, "log Kow", "비이온화 TSCF 측정 상한 — 실질 한계")
    rec("1_ad", "shoot_span_n", t_n, "행", "Schriever 2020 비이온화 행")
    rec("1_ad", "ionisation_gate_fn", 0.10, "-", "중성분율 하한 (보고서 §3.0a 기준)")
    for v in ORDER:
        sel = [r for r in rows if r["verdict"] == v]
        rec("1_ad", f"n_{v.split()[0].lower().strip('(')}", len(sel), "행", v)
    rec("1_ad", "n_total", len(rows), "행", "돌아간 전체 행")

    ext = [r for r in rows if r["verdict"] == V_EXT]
    rec("2_edges", "n_below_low_edge", len([r for r in ext if r["log_kow"] < r_lo]),
        "종", "저 log Kow 이탈 — 이번에 새로 생긴 이탈 방향")
    rec("2_edges", "n_above_high_edge", len([r for r in ext if r["log_kow"] > r_hi]),
        "종", "고 log Kow 이탈")

    x = [p["bat_A"] for p in pairs]
    rec("3_rank", "n", len(pairs), "종", "어류 BCF 가 존재하는 물질만")
    rec("3_rank", "rho_root", round(_rho(x, [p["root"] for p in pairs]), 3), "-", "")
    rec("3_rank", "rho_straw", round(_rho(x, [p["straw"] for p in pairs]), 3), "-", "")

    bio = [r for r in rows if r["bpc"] in ("B", "vB")]
    rec("4_bioacc", "n_classified", len(bio), "종", "규제기관 생물축적성 분류")
    rec("4_bioacc", "max_tscf", max((r["tscf"] for r in bio if r["tscf"]), default=0.0),
        "-", "이들 중 최대 TSCF — 낟알 도달 여부의 척도")

    if sweep is not None:
        rec("5_meaning", "half_equilibrium_logkow", round(half_meaning_edge(sweep), 2),
            "log Kow", "뿌리 BAF 가 K_PW 의 절반으로 떨어지는 지점")
        for lk in (4.0, 6.25, 8.25, 10.5):
            i = int(np.argmin(np.abs(sweep["grid"] - lk)))
            rec("5_meaning", f"ratio_at_logkow_{lk:g}",
                round(float(sweep["ratio"][i]), 4), "-",
                "뿌리 BAF / K_PW (1 = 평형)")
    else:
        # --fast 는 스윕을 돌리지 않는다. 그렇다고 이미 계산돼 있던 5_meaning 행을
        # 지우면, 브리핑이 인용하는 6.83 이 «빠른 재실행» 한 번에 조용히 사라진다.
        # 그래서 기존 값을 그대로 옮겨 적고, 새로 계산한 것이 아님을 표시한다.
        recs.extend(_carry_over(path, "5_meaning"))

    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["section", "key", "value", "unit", "note"])
        w.writeheader()
        for r in recs:
            w.writerow(r)
    return path


# ---------------------------------------------------------------------------
def main(fast=False):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    L = _labels()

    rows = load_census()
    sp = spans()
    pairs = rank_pairs(rows)
    os.makedirs(FIG_DIR, exist_ok=True)

    print(f"\ncensus {len(rows)}행 읽음. 판정:",
          {v: len([r for r in rows if r['verdict'] == v]) for v in ORDER})
    print(f"순위상관에 들어가는 물질 {len(pairs)}종 "
          f"(어류 BCF 없는 {len(rows) - len(pairs)}행 제외)")

    made = []
    made.append(fig_domain(rows, sp, L, os.path.join(FIG_DIR, "policy_biocide_domain.png")))
    made.append(fig_tscf(rows, sp, L, os.path.join(FIG_DIR, "policy_biocide_tscf.png")))
    made.append(fig_rank(pairs, L, os.path.join(FIG_DIR, "policy_biocide_rank.png")))

    sweep = None
    if not fast:
        print("\nlog Kow 스윕 (그림 4 — 유일하게 새로 계산하는 부분):")
        sweep = run_sweep()
        made.append(fig_ratio(sweep, sp, L,
                              os.path.join(FIG_DIR, "policy_biocide_meaning.png")))
        print(f"  → 뿌리 BAF 가 K_PW 의 절반이 되는 지점: "
              f"log Kow {half_meaning_edge(sweep):.2f}")

    p = write_csv(rows, sp, pairs, sweep)
    print("\n작성됨:")
    for m in made:
        print("  " + os.path.relpath(m, ROOT_DIR))
    print("  " + os.path.relpath(p, ROOT_DIR))
    if fast:
        kept = len(_carry_over(p, "5_meaning"))
        print(f"\n(--fast: 그림 4 는 생략. 5_meaning 수치는 "
              + (f"이전 실행 값 {kept}행을 그대로 옮겼습니다)" if kept
                 else "아직 없습니다 -- --fast 없이 한 번 돌리세요)"))
    return made, p


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fast", action="store_true",
                    help="log Kow 스윕(그림 4)을 생략한다")
    main(fast=ap.parse_args().fast)
