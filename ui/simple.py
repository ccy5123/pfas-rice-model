"""General-audience (Korean) Simple-mode view. Split out of app.py (HANDOFF P3-1)."""
import streamlit as st

import model_api as api
import plots

from ui.common import (_DISCLAIMER_KO, _cong_label_ko, _nearest_index,
                       _render_inverse_estimator, _glossary_md, _png_bytes, _html_bytes)


def render(cfg):
    """Render the Simple (general-audience, Korean) view from a populated cfg."""
    congener = cfg.congener
    res = cfg.res
    obs = cfg.obs
    bio_baf = cfg.bio_baf
    preset_word = cfg.preset_word
    use_custom_tables = cfg.use_custom_tables
    E_m = cfg.E_m
    fxy_source = cfg.fxy_source
    biomass = cfg.biomass
    # ---- plain-language headline -------------------------------------------
    grain_c = float(res["conc"]["grain"][-1])
    root_c = float(res["conc"]["root"][-1])
    straw_c = float(res["straw"][-1])
    grain_baf = res["baf_final"]["grain"]
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
    m1, m2, m3 = st.columns(3)
    m1.metric("뿌리 속", f"{root_c:.2g} µg/kg", _rng_ko(bands["root"]), delta_color="off")
    m2.metric("짚(줄기+잎) 속", f"{straw_c:.2g} µg/kg", _rng_ko(bands["straw"]), delta_color="off")
    m3.metric("낟알(먹는 쌀) 속", f"{grain_c:.2g} µg/kg", _rng_ko(bands["grain"]), delta_color="off")
    st.caption(f"각 수치 아래 범위는 모델의 **대략적 예측 불확실성**입니다 — 실제 측정값과 약 "
               f"**{_fold:.0f}배**까지 차이날 수 있습니다 (예측값이며 실측이 아닙니다).")

    _where_ko = {"roots": "뿌리", "straw (stems + leaves)": "짚(줄기+잎)", "grain": "낟알"}[where_most]
    _lead = ("입력하신 **성장 + 오염 표**를 바탕으로, " if use_custom_tables
             else f"선택하신 **{preset_word}** 오염 수준에서, ")
    gb = bands["grain"]
    # BAF<1 reads awkwardly as "약 0.2배"; phrase it as "lower than the soil water" (P2-4).
    if grain_baf >= 1.0:
        _baf_phrase = f"토양수 농도의 약 **{grain_baf:.1f}배**"
    elif grain_baf > 0:
        _baf_phrase = f"토양수 농도보다 **낮음**(약 1/{1.0 / grain_baf:.0f} 수준)"
    else:
        _baf_phrase = "토양수 농도보다 **매우 낮음**"
    st.info(
        _lead +
        f"이 모델은 벼 **낟알**에 {congener}가 약 **{grain_c:.2g} µg/kg** "
        f"(대략 {gb['lo']:.2g}–{gb['hi']:.2g}) 들어 있을 것으로 추정합니다 "
        f"({_baf_phrase}). 대부분의 화학물질은 **{_where_ko}**에 남습니다. "
        f"이 값은 **대략적 모델 예측**이라 실측과 수배 차이날 수 있습니다.")
    st.caption("예시용 모델 추정치이며 — 식품안전·건강 판단이 아닙니다.")

    s_tabs = st.tabs(["🗺️ 어디로 가나", "📈 시간에 따른 축적",
                      "📊 얼마나 쌓이나", "🔎 거꾸로 추정", "ℹ️ 안내 & 용어"])

    # ---- Simple tab 1: the plant + soil map --------------------------------
    with s_tabs[0]:
        cc1, cc2 = st.columns([1, 1])
        animate = cc1.checkbox("▶ 한 철 재생", value=False,
                               help="벼가 자라는 동안 PFAS가 하루하루 쌓이는 모습을 봅니다.")
        day = cc2.slider("이앙 후 일수", float(res["t"][0]), float(res["t"][-1]),
                         float(res["t"][-1]), 1.0, disabled=animate)
        if animate:
            st.plotly_chart(plots.fig_schematic_animated(res, "conc", lang="ko"),
                            width="stretch", theme=None)
        else:
            ti = _nearest_index(res["t"], day)
            st.plotly_chart(plots.fig_schematic_from_res(res, "conc", ti, obs=None, lang="ko"),
                            width="stretch", theme=None)
        st.caption("벼와 논 흙을 실제 비율로 그렸습니다. **색이 진할수록(뜨거울수록) 그 부위에 PFAS 농도가 높습니다.** "
                   "날짜 슬라이더를 끌거나 ▶를 눌러, 이앙부터 수확까지 화학물질이 어디에 쌓이는지 보세요. "
                   "색은 **부위별 농도**라서 잎 농도가 높으면 잎이 가장 뜨겁게 보일 수 있습니다 — "
                   "위 요약의 '대부분 ○○'은 뿌리·짚·낟알 기준이고, **짚 = 줄기와 잎의 평균**입니다.")

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

    # ---- Simple tab 4: work backwards (Bayesian inverse estimate) -----------
    with s_tabs[3]:
        _render_inverse_estimator(congener, E_m_mV=E_m, f_xy_source=fxy_source,
                                  biomass=biomass, key="inv_simple", simple=True)

    # ---- Simple tab 5: about & glossary ------------------------------------
    with s_tabs[4]:
        st.markdown(
            "### 이 도구가 하는 일\n"
            "벼가 흙에서 물과 녹아 있는 화학물질을 빨아들여 줄기 위로 올리고 낟알에 저장하는 과정을 "
            "**메커니즘 모델**로 계산합니다. PFAS 화학물질과 논의 오염 정도를 주면, 한 철 동안 뿌리·짚·"
            "먹는 낟알에 쌓이는 양을 추정합니다.\n\n"
            "### 보는 법\n"
            "- **🗺️ 어디로 가나** — 식물 그림; 색이 뜨거울수록 PFAS가 많음.\n"
            "- **📈 시간에 따른 축적** — 이앙부터 수확까지 농도 변화.\n"
            "- **📊 얼마나 쌓이나** — 최종 농도와 실제 측정값과의 비교.\n"
            "- **🔎 거꾸로 추정** — 실험실 측정값이 있으면 토양수 오염도를 (불확실성 범위와 함께 — "
            "베이지안 추정) 역추정.\n\n"
            "### 쉬운 용어 사전")
        st.markdown(_glossary_md(ko=True))
        st.warning(_DISCLAIMER_KO)
        st.caption("수식·파라미터·토양 결합(HYDRUS-1D)·구조(SMILES) 입력이 필요하면 사이드바의 "
                   "**전문가/고급 모드**를 켜세요.")

    # ---- downloads (Simple) ------------------------------------------------
    with st.expander("⬇️ 결과 내려받기"):
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
