"""General-audience (Korean) Simple-mode view. Split out of app.py (HANDOFF P3-1)."""
import streamlit as st
import streamlit.components.v1 as components

import model_api as api
import plots
import plot_svg

from ui.common import (_DISCLAIMER_KO, _cong_label_ko, _nearest_index, _simulate,
                       _render_inverse_estimator, _glossary_md, _png_bytes, _html_bytes,
                       _summary_html, signal_badge_html)


def render(cfg):
    """Render the Simple (general-audience, Korean) view from a populated cfg."""
    congener = cfg.congener
    res = cfg.res
    obs = cfg.obs
    bio_baf = cfg.bio_baf
    preset_word = cfg.preset_word
    E_m = cfg.E_m
    fxy_source = cfg.fxy_source
    biomass = cfg.biomass
    # ---- hands-on guidance (for the try-it-yourself audience) ---------------
    with st.expander("💡 처음이신가요? 3단계로 해보세요", expanded=False):
        st.markdown(
            "1. 왼쪽에서 **빠른 시나리오**를 하나 고르세요 (또는 '✏️ 직접 설정'에서 물질·오염도 선택).\n"
            "2. 아래 **🗺️ 어디로 가나** 지도에서 색이 진한 부위 = PFAS가 많이 쌓이는 곳을 보세요.\n"
            "3. **🍚 밥으로 먹으면** 탭에서 안전기준 대비 수준을, **⚖️ 물질 비교**에서 물질별 차이를 확인하세요.\n\n"
            "모든 수치는 **대략적 예측**입니다 — 부위·물질 간 **상대 비교**로 읽어 주세요.")
    # ---- plain-language headline -------------------------------------------
    grain_c = float(res["conc"]["grain"][-1])
    root_c = float(res["conc"]["root"][-1])
    straw_c = float(res["straw"][-1])
    tops = {"roots": root_c, "straw (stems + leaves)": straw_c, "grain": grain_c}
    where_most = max(tops, key=tops.get)

    # Honest a-priori predictive band (×/÷ ~7): the model's out-of-sample error is
    # large, so every absolute number is shown with a range + a several-fold caveat
    # rather than a precise-looking single figure (HANDOFF P1).
    bands = {k: api.predictive_band(v) for k, v in
             (("root", root_c), ("straw", straw_c), ("grain", grain_c))}
    _fold = api.uncertainty_factor()

    def _rng_ko(b):
        return f"대략 {b['lo']:.2g}–{b['hi']:.2g} µg/kg"

    st.subheader(f"{_cong_label_ko(congener)}")

    # Headline RESULT STRIP first (big signal + one-line takeaway), so the answer
    # lands before the numbers. Colour-blind-safe badge (colour + shape + label).
    _where_ko = {"roots": "뿌리", "straw (stems + leaves)": "짚(줄기+잎)", "grain": "낟알"}[where_most]
    _intake = api.intake_fraction(grain_c, congener=congener)
    _pct_txt = f"기준의 {_intake['percent']:.0f}%" if _intake["in_group"] else "참고"
    _badge = signal_badge_html(_intake["signal"], ko=True, extra=_pct_txt)
    _cls = {"green": "safe", "amber": "warn", "red": "dang"}.get(_intake["signal"], "")
    st.markdown(
        f"<div class='pfas-result {_cls}'>{_badge}"
        f"<div class='pfas-result-text'>대부분의 PFAS는 <b>{_where_ko}</b>에 남고, 먹는 "
        f"<b>낟알</b>엔 약 <b>{grain_c:.2g} µg/kg</b> 들어 있을 것으로 추정합니다.</div></div>",
        unsafe_allow_html=True)

    # then the three tissue numbers as cards
    m1, m2, m3 = st.columns(3)
    m1.metric("🥕 뿌리 속", f"{root_c:.2g} µg/kg", _rng_ko(bands["root"]), delta_color="off")
    m2.metric("🌿 짚(줄기+잎) 속", f"{straw_c:.2g} µg/kg", _rng_ko(bands["straw"]), delta_color="off")
    m3.metric("🌾 낟알(먹는 쌀) 속", f"{grain_c:.2g} µg/kg", _rng_ko(bands["grain"]), delta_color="off")

    st.caption(f"신호등은 EFSA 건강기반 안전기준 대비 수준입니다(색·모양·라벨로 표시). 숫자는 "
               f"**대략적 예측**(실측과 ~{_fold:.0f}배 차이 가능)이며 법적 식품기준이 아닙니다. "
               f"자세한 비교는 🍚·⚖️·🔬 탭에서.")

    s_tabs = st.tabs(["🗺️ 어디로 가나", "📈 시간에 따른 축적", "📊 얼마나 쌓이나",
                      "⚖️ 물질 비교", "🍚 밥으로 먹으면", "🔎 거꾸로 추정", "🔬 신뢰도 & 안내"])

    # ---- Simple tab 1: the plant + soil map --------------------------------
    with s_tabs[0]:
        day = st.slider("이앙 후 일수", float(res["t"][0]), float(res["t"][-1]),
                        float(res["t"][-1]), 1.0,
                        help="슬라이더를 옮기면 그 날짜의 축적 상태로 색이 바뀝니다.")
        ti = _nearest_index(res["t"], day)
        # grain safety chip (안전/주의/초과) — only for the EFSA-4 congeners, so a
        # borrowed reference signal is never shown as a food-safety verdict.
        _grain_sig = _intake["signal"] if _intake.get("in_group") else None
        components.html(plot_svg.plant_svg_from_res(res, ti, lang="ko",
                                                    grain_signal=_grain_sig), height=560)
        st.markdown(
            f"<span class='pfas-caveat'>ⓘ 대략적 예측 · 실측과 ~{api.uncertainty_factor():.0f}배 "
            f"차이 가능</span>", unsafe_allow_html=True)
        st.caption("**색이 진할수록 그 부위에 PFAS가 많습니다.** 날짜 슬라이더로 한 철 축적을 "
                   "보세요. (짚 = 줄기+잎 평균)")

    # ---- Simple tab 2: build-up over time ----------------------------------
    with s_tabs[1]:
        st.plotly_chart(plots.fig_buildup_plain(res, lang="ko"), width="stretch")
        st.caption("한 철 동안 각 식물 부위의 PFAS 농도가 어떻게 변하는지. **낟알**은 형성된 뒤(개화 무렵)부터 "
                   "PFAS를 흡수하기 시작해 수확까지 계속 쌓입니다. 곡선은 **대략적 모델 예측**으로, "
                   "실측과 수배 차이날 수 있습니다.")

    # ---- Simple tab 3: how much builds up ----------------------------------
    with s_tabs[2]:
        st.plotly_chart(plots.fig_where_plain(res, lang="ko", band=True), width="stretch")
        st.caption(f"수확 시 각 부위의 PFAS 농도. 보통 뿌리에 가장 많이 남고, 먹는 낟알까지 "
                   f"얼마나 도달하는지는 화학물질에 따라 다릅니다. 막대의 **회색 오차선**은 모델의 "
                   f"대략적 예측 불확실성(실측과 약 {_fold:.0f}배까지 차이날 수 있음)입니다.")
        if obs:
            with st.expander("🔬 실제 측정값과 비교 (Yamazaki 2023)"):
                st.plotly_chart(plots.fig_baf(res, obs, lang="ko"), width="stretch")
                st.caption("막대는 모델의 축적 배수를 출판된 온실 벼 연구(Yamazaki et al. 2023)의 측정값과 "
                           "비교한 것입니다. 막대가 비슷할수록 이 화학물질에 대해 모델이 실제 데이터와 잘 맞습니다.")

    # ---- Simple tab 4: cross-chemical comparison (policy hook C) -------------
    with s_tabs[3]:
        st.markdown("#### 물질에 따라 쌓이는 곳이 다릅니다")
        _rep = [c for c in ["PFBA", "PFHxA", "PFOA", "PFDA", "PFDoDA", "PFOS"]
                if c in api.CONGENERS]
        try:
            _cmp = {c: _simulate(c, Cwo=cfg.Cwo_const, season=cfg.season,
                                 measured_forcing=cfg.measured, E_m_mV=E_m,
                                 f_xy_source=fxy_source, biomass=biomass) for c in _rep}
            st.plotly_chart(plots.fig_congener_compare(_cmp, lang="ko", order=_rep),
                            width="stretch")
            st.caption(
                f"같은 오염 수준(토양수 {cfg.Cwo_const:g} µg/L)에서 물질별 수확기 농도. "
                "**짧은사슬(PFBA·PFHxA)**은 짚·낟알까지 잘 이동해 먹는 부위에 도달하기 쉽고, "
                "**긴사슬(PFDA·PFDoDA)**은 뿌리에 강하게 잔류해 낟알로는 덜 갑니다. "
                "PFOS는 대표적 술폰산 계열입니다. 값 차이가 매우 커서 세로축은 **로그(배수) 눈금**이며, "
                "막대 높이는 **대략적 모델 예측**입니다.")
        except Exception as e:                                  # noqa: BLE001
            st.warning(f"물질 비교를 계산할 수 없습니다: {e}")

    # ---- Simple tab 5: on the dinner table (EFSA TWI intake reference) -------
    with s_tabs[4]:
        st.markdown("#### 이 쌀을 먹으면 안전기준 대비 어느 정도일까요? (참고)")
        cA, cB = st.columns(2)
        rice_g = cA.slider("하루 쌀 섭취량 [g/일]", 30, 400,
                           int(api.DEFAULT_RICE_INTAKE_G_PER_DAY), 10,
                           help="국민건강영양조사·KOSIS 기준 국내 1인 1일 쌀 섭취량은 약 150 g입니다.")
        bw = cB.slider("체중 [kg]", 20, 100, int(api.DEFAULT_BODY_WEIGHT_KG), 5,
                       help="1인 기준 체중. 어린이 등 체중이 작을수록 체중당 섭취 비율이 높아집니다.")
        info = api.intake_fraction(grain_c, congener=congener,
                                   rice_intake_g_day=float(rice_g), body_weight_kg=float(bw))
        _gb = signal_badge_html(info["signal"], ko=True,
                                extra=(f"기준의 {info['percent']:.0f}%" if info["in_group"] else "참고"))
        st.markdown(_gb, unsafe_allow_html=True)
        st.plotly_chart(plots.fig_intake_gauge(info, lang="ko"), width="stretch", theme=None)
        _band = api.predictive_band(grain_c)
        _lo = api.intake_fraction(_band["lo"], congener=congener,
                                  rice_intake_g_day=float(rice_g), body_weight_kg=float(bw))
        _hi = api.intake_fraction(_band["hi"], congener=congener,
                                  rice_intake_g_day=float(rice_g), body_weight_kg=float(bw))
        mA, mB = st.columns([1, 2])
        mA.metric("주간 섭취량 / 안전기준", f"{info['percent']:.0f}%",
                  f"대략 {_lo['percent']:.0f}–{_hi['percent']:.0f}%", delta_color="off")
        mB.markdown(
            f"낟알 **{grain_c:.2g} µg/kg**인 쌀을 하루 **{rice_g} g** 드시면 PFAS 주간 섭취량은 "
            f"EFSA 안전기준(**{info['twi']} ng/kg 체중/주**, 4종 합)의 약 **{info['percent']:.0f}%** "
            f"수준입니다.")
        if info["in_group"] is False:
            st.warning(f"⚠ **{congener}**는 EFSA 4종(PFOA·PFOS·PFNA·PFHxS)에 **미포함** — 위 %는 "
                       f"같은 기준을 빌린 **참고치**입니다.")
        with st.expander("ℹ️ 이 비교를 읽는 법 (중요)"):
            st.markdown(
                "- 쌀에 대한 **법적 PFAS 기준(MRL)은 EU·한국 모두 없습니다.** 이 화면은 식품기준 초과가 "
                "아니라 **건강기반 섭취기준(EFSA TWI) 환산 참고치**입니다.\n"
                "- EFSA 기준은 **4종 합**인데 여기선 한 물질만 비교 → **‘이 물질 하나만으로도’ 보수적 가정**.\n"
                "- **쌀만 섭취원으로 가정**한 단순 계산 (실제론 물·타 식품에서도 섭취).\n"
                "- 낟알 농도 예측 자체가 실측과 **수배 차이날 수 있습니다.**")
        st.caption("출처: EFSA 2020 그룹 TWI 4.4 ng/kg 체중/주 (PFOA·PFOS·PFNA·PFHxS, "
                   "doi:10.2903/j.efsa.2020.6223) · 쌀 섭취량 KOSIS/국민건강영양조사. "
                   "🔴빨강 = 기준의 100% 이상, 🟡노랑 = 10–100%, 🟢초록 = 10% 미만.")

    # ---- Simple tab 6: work backwards (Bayesian inverse estimate) -----------
    with s_tabs[5]:
        _render_inverse_estimator(congener, E_m_mV=E_m, f_xy_source=fxy_source,
                                  biomass=biomass, key="inv_simple", simple=True)

    # ---- Simple tab 7: model trust (validation status) + about & glossary ---
    with s_tabs[6]:
        st.markdown("### 🔬 이 모델을 얼마나 믿을 수 있나요?")
        st.markdown(
            "이 모델은 **정성적 경향(어디에 더 쌓이는가)** 은 잘 맞지만 **절대 수치는 대략치**입니다. "
            "발표·정책 판단 시 아래 검증 상태를 함께 봐 주세요.")
        st.markdown(
            "| 항목 | 상태 | 근거 / 설명 |\n"
            "|---|---|---|\n"
            "| 음이온 배제(뿌리 흡수 장벽) · 질량보존 | 🟢 **구조적으로 검증** | "
            "막전위 기반 물리(eᴺ≈107); 모델이 반드시 만족 |\n"
            "| 짧은사슬 PFAS 이행 **방향** | 🟢 **실측과 일치** | Yamazaki 2023 온실 벼 |\n"
            "| 토양→공극수 결합 | 🟢 **실제 HYDRUS-1D 엔진 연동** | Method A (일방향) |\n"
            f"| **절대 농도**(예측력) | 🟡 **대략치 (실측과 약 ×{api.uncertainty_factor():.0f})** | "
            "미보정 예측 오차 log10 RMSE ≈0.85 |\n"
            "| 긴사슬(C10~C12) 짚·낟알 축적 | 🔴 **탐색적·미해결** | 문헌·메커니즘 공백 |\n"
            "| 물질별 이행계수 f_xy 절대값 | 🟡 **데이터셋 의존** | 논·품종·조건마다 다름 |\n")
        st.caption("🟢 검증됨 · 🟡 대략치/조건의존 · 🔴 미해결(연구 진행 중). 절대 수치보다 **부위 간 "
                   "상대 비교**와 **물질 간 경향**을 신뢰하세요.")
        if obs:
            st.markdown("**실제 측정값과의 비교 (Yamazaki 2023 온실 벼)** — 막대가 비슷할수록 이 물질에 "
                        "대해 모델이 데이터와 잘 맞습니다.")
            st.plotly_chart(plots.fig_baf(res, obs, lang="ko"), width="stretch", key="baf_trust")
        else:
            st.caption(f"{congener}에 대한 공개 실측 BAF가 없어 이 물질은 모델 예측만 표시됩니다.")
        st.info(
            "**정직한 한계** — 이 모델의 검증은 대부분 **한 실험(Yamazaki, in-sample)** 에 맞춘 것이고, "
            "독립 데이터에 대한 예측 검증은 제한적입니다(단일 clean 데이터셋). 결정적 확증에는 아직 "
            "**현장 실험 데이터**(부위별 물관액/토양수 비율 등)가 필요합니다. 따라서 **정책 스크리닝·"
            "우선순위·교육용**으로는 유용하지만, **개별 지점의 규제 판단 근거로 단독 사용은 부적절**합니다.")
        st.divider()
        st.markdown(
            "### 이 도구가 하는 일\n"
            "벼가 흙에서 물과 녹아 있는 화학물질을 빨아들여 줄기 위로 올리고 낟알에 저장하는 과정을 "
            "**메커니즘 모델**로 계산합니다. PFAS 화학물질과 논의 오염 정도를 주면, 한 철 동안 뿌리·짚·"
            "먹는 낟알에 쌓이는 양을 추정합니다.\n\n"
            "### 보는 법\n"
            "- **🗺️ 어디로 가나** — 식물 그림; 색이 뜨거울수록 PFAS가 많음.\n"
            "- **📈 시간에 따른 축적** — 이앙부터 수확까지 농도 변화.\n"
            "- **📊 얼마나 쌓이나** — 최종 농도와 실제 측정값과의 비교.\n"
            "- **⚖️ 물질 비교** — 여러 PFAS 물질이 뿌리·짚·낟알에 쌓이는 정도를 한눈에.\n"
            "- **🍚 밥으로 먹으면** — 예측된 낟알 농도를 EFSA 건강기반 섭취기준에 환산한 **참고** 비교.\n"
            "- **🔎 거꾸로 추정** — 실험실 측정값이 있으면 토양수 오염도를 (불확실성 범위와 함께 — "
            "베이지안 추정) 역추정.\n"
            "- **🔬 신뢰도 & 안내** — 지금 이 화면. 모델 검증 상태와 용어.\n\n"
            "### 쉬운 용어 사전")
        st.markdown(_glossary_md(ko=True))
        st.warning(_DISCLAIMER_KO)
        st.caption("수식·파라미터·토양 결합(HYDRUS-1D)·구조(SMILES) 입력이 필요하면 사이드바의 "
                   "**전문가/고급 모드**를 켜세요.")

    # ---- downloads (Simple) ------------------------------------------------
    with st.expander("⬇️ 결과 내려받기 · 발표용 요약", expanded=False):
        # Presentation one-pager: plant map + EFSA gauge + numbers + caveats,
        # a single offline HTML file (open in a browser, print to PDF for a slide).
        one_pager, why = _summary_html(cfg, lang="ko")
        if one_pager is not None:
            st.download_button("🖼️ 발표용 요약 1장 (HTML)", one_pager,
                               file_name=f"{congener}_요약.html", mime="text/html",
                               help="식물 지도 + 안전기준 게이지 + 핵심 수치 + 한계를 한 파일로. "
                                    "브라우저에서 열어 PDF로 인쇄하면 슬라이드 한 장이 됩니다.")
        _intake_md = api.intake_fraction(grain_c, congener=congener)
        st.download_button("📄 요약 텍스트 (Markdown)",
                           api.summary_report_md(res, congener=congener, obs=obs,
                                                 intake=_intake_md, scenario=preset_word, lang="ko"),
                           file_name=f"{congener}_요약.md", mime="text/markdown")
        cda, cdb = st.columns(2)
        cda.download_button("요약 표 (CSV)", api.summary_csv(res, obs, bio_baf),
                            file_name=f"{congener}_summary.csv", mime="text/csv")
        cdb.download_button("전체 시계열 (CSV)", api.timeseries_csv(res),
                            file_name=f"{congener}_timeseries.csv", mime="text/csv")
        _map_fig = plots.fig_schematic_from_res(res, "conc", -1, lang="ko")
        png, why = _png_bytes(_map_fig)
        if png is not None:
            st.download_button("식물 지도 (PNG)", png, file_name=f"{congener}_map.png", mime="image/png")
        else:
            # No kaleido/Chrome -> offer an interactive HTML instead (always works).
            html, _ = _html_bytes(_map_fig)
            if html is not None:
                st.download_button("식물 지도 (대화형 HTML)", html,
                                   file_name=f"{congener}_map.html", mime="text/html")
            st.caption("정적 **PNG** 내보내기는 선택 패키지 `kaleido`가 필요합니다 — "
                       "`pip install kaleido && plotly_get_chrome` 후 다시 실행하세요. "
                       "그동안 위 **대화형 HTML**(브라우저에서 확대·툴팁 가능)과 CSV는 그대로 받을 수 있습니다.")
