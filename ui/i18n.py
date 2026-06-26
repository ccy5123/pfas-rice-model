"""Lightweight i18n string table for the UI (HANDOFF P3-2).

Centralizes the **bilingual** UI copy that used to be scattered through app.py as
inline ``"…한국어…" if ko else "…English…"`` ternaries, so adding/adjusting a
language is a table edit instead of a code hunt. Use ``t(key, lang, **fmt)``.

Scope: the *shared* widgets that render in BOTH modes (the inverse-estimator and
custom-tables panels) plus the paired data constants (disclaimer, congener labels,
contamination presets, glossary). The mode-specific page copy lives, by design, in
its own single-language module — ``ui/simple.py`` (Korean) and ``ui/expert.py``
(English) — and ``plots.py`` keeps its own ``lang=`` argument for figure labels.
"""

LANGS = ("en", "ko")


def t(key, lang="en", **fmt):
    """Look up ``key`` for ``lang`` (falling back to English), and ``str.format`` it
    with ``fmt`` when keyword arguments are given."""
    entry = STRINGS[key]
    s = entry.get(lang) or entry["en"]
    return s.format(**fmt) if fmt else s


# --- paired data constants (keyed by congener / preset label) ----------------
CONGENER_LABELS = {
    "en": {
        "PFBA":   "PFBA — short-chain acid (C4)",
        "PFPeA":  "PFPeA — short-chain acid (C5)",
        "PFHxA":  "PFHxA — short-chain acid (C6)",
        "PFHpA":  "PFHpA — medium-chain acid (C7)",
        "PFOA":   "PFOA — common 'forever chemical' acid (C8)",
        "PFNA":   "PFNA — long-chain acid (C9)",
        "PFDA":   "PFDA — long-chain acid (C10)",
        "PFUnDA": "PFUnDA — long-chain acid (C11)",
        "PFDoDA": "PFDoDA — long-chain acid (C12)",
        "PFBS":   "PFBS — short-chain sulfonate (C4)",
        "PFHxS":  "PFHxS — sulfonate (C6)",
        "PFOS":   "PFOS — common 'forever chemical' sulfonate (C8)",
        "GenX":   "GenX — newer PFOA replacement (ether)",
    },
    "ko": {
        "PFBA":   "PFBA — 단쇄 카복실산 (C4)",
        "PFPeA":  "PFPeA — 단쇄 카복실산 (C5)",
        "PFHxA":  "PFHxA — 단쇄 카복실산 (C6)",
        "PFHpA":  "PFHpA — 중쇄 카복실산 (C7)",
        "PFOA":   "PFOA — 대표적 '영원한 화학물질' 카복실산 (C8)",
        "PFNA":   "PFNA — 장쇄 카복실산 (C9)",
        "PFDA":   "PFDA — 장쇄 카복실산 (C10)",
        "PFUnDA": "PFUnDA — 장쇄 카복실산 (C11)",
        "PFDoDA": "PFDoDA — 장쇄 카복실산 (C12)",
        "PFBS":   "PFBS — 단쇄 술폰산 (C4)",
        "PFHxS":  "PFHxS — 술폰산 (C6)",
        "PFOS":   "PFOS — 대표적 '영원한 화학물질' 술폰산 (C8)",
        "GenX":   "GenX — PFOA 대체물질 (에터)",
    },
}

# Contamination presets. en: {label: µg/L}. ko: {label: (µg/L, short word reused in
# the headline sentence)}. Kept in this shape because common.py exposes them as the
# legacy _PRESETS / _PRESETS_KO that app.py already consumes.
PRESETS = {
    "en": {
        "Low — lightly contaminated (0.1 µg/L)": 0.1,
        "Medium — moderately contaminated (1 µg/L)": 1.0,
        "High — heavily contaminated (10 µg/L)": 10.0,
    },
    "ko": {
        "낮음 — 약하게 오염 (0.1 µg/L)": (0.1, "낮은"),
        "중간 — 보통 오염 (1 µg/L)": (1.0, "중간"),
        "높음 — 심하게 오염 (10 µg/L)": (10.0, "높은"),
    },
}


# --- string templates (use t(key, lang, **fmt)) ------------------------------
STRINGS = {
    # ---- disclaimer (top banner + footer) ----
    "disclaimer": {
        "en": ("**Research & educational model — illustrative estimates only.** "
               "This is **not** a regulatory, food-safety, or health determination. "
               "Do **not** use it for real exposure or safety decisions."),
        "ko": ("**연구·교육용 모델 — 예시 추정치일 뿐입니다.** "
               "규제·식품안전·건강 판단이 **아니며**, 실제 노출·안전 결정에 **사용하지 마세요**."),
    },
    # ---- plain-language glossary (About tab + Simple expander) ----
    "glossary": {
        "en": (
            "| Term you'll see | What it means in plain words |\n"
            "|---|---|\n"
            "| **PFAS** | A family of long-lasting synthetic 'forever chemicals'. |\n"
            "| **Pore-water level** *(Cwᵒ)* | How much PFAS is dissolved in the soil water around the roots. |\n"
            "| **Build-up factor** *(BAF)* | How many times more concentrated the PFAS is in the plant tissue than in the soil water. A factor of 2 means the tissue holds twice the water's level. |\n"
            "| **Roots / Straw / Grain** | The plant parts. *Straw* = the stems + leaves together. *Grain* = the edible brown rice. |\n"
            "| **Concentration** *(µg/kg)* | Micrograms of PFAS per kilogram of plant tissue. |\n"
            "| **Congener** | One specific PFAS chemical (e.g. PFOA, PFOS). Longer carbon chains generally stick to the plant more. |\n"
            "| **Uptake / translocation** | How the chemical gets into the roots and then moves up into the shoot and grain. |\n"
            "| **Bayesian estimate** | Working backwards from a measurement to the most likely cause, **with an uncertainty range** instead of a single number. |\n"),
        "ko": (
            "| 보게 될 용어 | 쉬운 설명 |\n"
            "|---|---|\n"
            "| **PFAS** | 잘 분해되지 않는 인공 '영원한 화학물질' 무리. |\n"
            "| **토양수 농도** *(Cwᵒ)* | 뿌리 주변 토양수에 녹아 있는 PFAS의 양. |\n"
            "| **축적 배수** *(BAF)* | 식물 조직이 토양수보다 PFAS를 몇 배나 진하게 모았는지. 2배면 토양수의 두 배. |\n"
            "| **뿌리 / 짚 / 낟알** | 식물 부위. *짚* = 줄기 + 잎, *낟알* = 먹는 현미. |\n"
            "| **농도** *(µg/kg)* | 식물 조직 1 kg당 PFAS 마이크로그램. |\n"
            "| **화학종(congener)** | 특정 PFAS 하나(예: PFOA, PFOS). 탄소 사슬이 길수록 대체로 식물에 더 잘 달라붙음. |\n"
            "| **흡수 / 이동** | 화학물질이 뿌리로 들어가 줄기·낟알로 올라가는 과정. |\n"
            "| **베이지안 추정** | 측정값에서 거꾸로 가장 가능성 높은 원인을 찾되, 하나의 숫자가 아니라 **불확실성 범위**까지 제시. |\n"),
    },

    # ---- inverse estimator ("work backwards") panel ----
    "inv.not_curated": {
        "en": ("This works with the curated chemicals — pick one of them in the sidebar "
               "(it needs the calibrated model, which a custom SMILES structure doesn't have)."),
        "ko": ("선별된 13종 화학물질에서만 작동합니다 — 사이드바에서 하나를 고르세요 "
               "(보정된 모델이 필요하며, 직접 입력한 SMILES 구조에는 없습니다)."),
    },
    "inv.intro": {
        "en": ("Already have a **lab result** for rice grown on contaminated land? Enter what was "
               "measured in the plant and this estimates **how contaminated the soil water likely "
               "was** — working the model backwards, with an **uncertainty range** (a Bayesian "
               "estimate), not just a single number."),
        "ko": ("오염된 땅에서 자란 벼의 **실험실 측정값**이 있으신가요? 식물에서 측정된 값을 입력하면 "
               "모델을 거꾸로 돌려 **토양수가 얼마나 오염됐을지**를 — 하나의 숫자가 아니라 "
               "**불확실성 범위**까지(베이지안 추정) — 추정합니다."),
    },
    "inv.in_root": {"en": "Measured in roots [µg/kg]", "ko": "뿌리 측정값 [µg/kg]"},
    "inv.in_straw": {"en": "Measured in straw (stems+leaves) [µg/kg]", "ko": "짚(줄기+잎) 측정값 [µg/kg]"},
    "inv.in_grain": {"en": "Measured in grain [µg/kg]", "ko": "낟알 측정값 [µg/kg]"},
    "inv.precision_label": {"en": "How precise are the measurements?", "ko": "측정값이 얼마나 정밀한가요?"},
    "inv.precision_help": {
        "en": "Sets the measurement+model uncertainty used in the estimate.",
        "ko": "추정에 쓰이는 측정+모델 불확실성을 설정합니다.",
    },
    "inv.estimate_btn": {"en": "📐 Estimate the contamination level", "ko": "📐 오염 수준 추정하기"},
    "inv.enter_first": {
        "en": "Enter at least one measured tissue concentration above, then press the button.",
        "ko": "위에 측정값을 하나 이상 입력한 뒤 버튼을 누르세요.",
    },
    "inv.press_estimate": {
        "en": "Press **Estimate** to run (or re-run after changing a value).",
        "ko": "**추정** 버튼을 누르세요 (값을 바꾼 뒤에는 다시 누르세요).",
    },
    "inv.preparing": {"en": "Preparing the estimate…", "ko": "오염 수준 추정 준비 중…"},
    "inv.running": {
        "en": "Running the model backwards… (step {done}/{total})",
        "ko": "모델을 거꾸로 계산하는 중… ({done}/{total} 단계)",
    },
    "inv.error": {"en": "Could not run the estimate: {e}", "ko": "추정을 실행할 수 없습니다: {e}"},
    "inv.metric_label": {"en": "Most likely soil-water level", "ko": "가장 가능성 높은 토양수 수준"},
    "inv.range95": {"en": "95% range {lo:.2g}–{hi:.2g}", "ko": "95% 범위 {lo:.2g}–{hi:.2g}"},
    "inv.range_unconstrained": {"en": "range unconstrained", "ko": "범위 미확정"},
    "inv.summary_lead": {
        "en": ("Given your measurements of **{congener}**, the model estimates the PFAS dissolved "
               "in the soil water was most likely **{med:.3g} µg/L**"),
        "ko": ("입력하신 **{congener}** 측정값으로 보면, 토양수에 녹아 있던 PFAS는 가장 가능성 높게 "
               "**{med:.3g} µg/L**"),
    },
    "inv.summary_ci": {
        "en": ", and we're 95% confident it was between **{lo:.2g}** and **{hi:.2g} µg/L**.",
        "ko": "이며, **{lo:.2g}~{hi:.2g} µg/L** 사이일 확률이 95%입니다.",
    },
    "inv.summary_no_ci": {
        "en": " (the measurements didn't pin down a clear range).",
        "ko": "입니다 (측정값으로는 뚜렷한 범위가 좁혀지지 않았습니다).",
    },
    "inv.summary_tail": {
        "en": " The spread of the curve below is the uncertainty.",
        "ko": " 아래 곡선의 퍼짐이 불확실성입니다.",
    },
    "inv.fit_row": {
        "en": "{name}: you {meas:.3g} vs model {model:.3g}",
        "ko": "{name}: 입력 {meas:.3g} vs 모델 {model:.3g}",
    },
    "inv.fit_caption": {
        "en": ("At the best estimate the model reproduces your inputs — {rows} (µg/kg). "
               "This estimates only the **contamination level**; it assumes the model's "
               "plant-uptake is right and can't separately tell apart water uptake vs how "
               "the chemical moves inside the plant (that needs a sap/soil-water measurement)."),
        "ko": ("최적 추정에서 모델이 입력값을 재현합니다 — {rows} (µg/kg). "
               "이는 **오염 수준**만 추정하며, 모델의 흡수 거동이 옳다고 가정합니다 — 특정 논의 실측이 아닌 예시입니다."),
    },

    # ---- custom tables (growth + Cwᵒ) panel ----
    "ct.intro": {
        "en": ("Enter your **own season** as two tables — **growth = organ FRESH weight** and "
               "**pore water = the absolute PFAS dissolved in the soil water (µg/L)** over time. "
               "Edit cells, add/remove rows, or paste from a spreadsheet; the days are interpolated "
               "onto the model timeline. Leave either table at its default to use the built-in value."),
        "ko": ("**나만의 한 철**을 두 표로 입력하세요 — **성장 = 기관별 신선중**, **공극수 = 토양수에 "
               "녹아 있는 PFAS 절대 농도(µg/L)**의 시간 변화. 셀을 편집·행 추가/삭제하거나 스프레드시트에서 "
               "붙여넣으세요; 날짜는 모델 타임라인으로 보간됩니다. 어느 한 표를 기본값으로 두면 내장값을 씁니다."),
    },
    "ct.growth_caption": {
        "en": "🌱 Growth — organ FRESH weight per hill", "ko": "🌱 성장 — 포기당 기관별 신선중",
    },
    "ct.growth_upload": {
        "en": "…or upload a growth CSV (day,root,stem,leaf,grain)",
        "ko": "…또는 성장 CSV 업로드 (day,root,stem,leaf,grain)",
    },
    "ct.growth_units": {"en": "Growth weight units", "ko": "성장 무게 단위"},
    "ct.cwo_caption": {
        "en": "💧 Pore-water contamination Cwᵒ(t) — absolute µg/L",
        "ko": "💧 공극수 오염 Cwᵒ(t) — 절대 µg/L",
    },
    "ct.cwo_upload": {"en": "…or upload a Cwᵒ CSV (day,Cwo)", "ko": "…또는 Cwᵒ CSV 업로드 (day,Cwo)"},
    "ct.density_md": {
        "en": ("**Compartment density** ρ [kg/L, fresh] — links the entered weight to tissue "
               "volume (rice leaf/culm hold air spaces ⇒ < 1; grain is denser ⇒ > 1)."),
        "ko": ("**구획 밀도** ρ [kg/L, 신선] — 입력한 무게를 조직 부피와 연결합니다 "
               "(벼 잎/줄기는 통기조직으로 < 1, 낟알은 더 조밀해 > 1)."),
    },
    "ct.read_error": {
        "en": "Using the default scenario — couldn't read the tables: {e}",
        "ko": "표를 읽지 못해 기본 시나리오를 사용합니다: {e}",
    },
    "ct.implied_volume": {
        "en": ("Implied end-of-season organ **volume** (fresh mass ÷ density): {rows}. "
               "The transport model integrates on fresh mass; density sets the mass↔volume scale."),
        "ko": ("추정된 수확기 기관 **부피** (신선중 ÷ 밀도): {rows}. "
               "수송 모델은 신선중(질량)으로 적분하며, 밀도는 질량↔부피 환산용입니다."),
    },
    "ct.panel_title": {
        "en": "📋 Your data tables — growth curve + pore-water contamination",
        "ko": "📋 내 데이터 표 — 성장 곡선 + 토양수 오염",
    },

    # ---- header (title + intro) ----
    "header.intro1": {
        "en": ("**PFAS are long-lasting 'forever chemicals'.** This tool estimates how much PFAS "
               "a rice plant takes up from a contaminated paddy's water/soil and **where it ends "
               "up** — roots, straw (stems+leaves), and the edible **grain**."),
        "ko": ("**PFAS는 잘 분해되지 않는 '영원한 화학물질'입니다.** 이 도구는 벼가 오염된 논의 물·흙에서 "
               "PFAS를 얼마나 흡수하고 **어디에 쌓이는지** — 뿌리, 짚(줄기+잎), 먹는 **낟알** — 를 추정합니다."),
    },
    "header.intro2": {
        "en": ("👉 **Start here:** pick a **chemical** and a **contamination level** on the left, then "
               "look at the **🗺️ Where it goes** map below. No chemistry background needed."),
        "ko": ("👉 **여기서 시작:** 왼쪽에서 **화학물질**과 **오염 수준**을 고른 뒤, 아래 "
               "**🗺️ 어디로 가나** 지도를 보세요. 화학 배경지식은 필요 없습니다."),
    },
    "header.expert_caption": {
        "en": ("Mechanistic 4-compartment dynamic model for permanently-anionic PFAS in paddy rice "
               "(IOC extension of the Trapp/Brunetti DPU). Parameters are measured/cited "
               "(docs/literature_db). Charts are interactive — hover, zoom, toggle. Outputs illustrative."),
        "ko": ("Mechanistic 4-compartment dynamic model for permanently-anionic PFAS in paddy rice "
               "(IOC extension of the Trapp/Brunetti DPU). Parameters are measured/cited "
               "(docs/literature_db). Charts are interactive — hover, zoom, toggle. Outputs illustrative."),
    },

    # ---- footer ----
    "footer.links": {
        "en": "[Source code]({repo}) · [Documentation]({docs})",
        "ko": "[소스 코드]({repo}) · [문서]({docs})",
    },
    "footer.cite": {
        "en": ("How to cite: PFAS–Rice Compartmental Uptake Model (IOC extension of the "
               "Trapp/Brunetti DPU framework), 2026."),
        "ko": ("인용: PFAS–Rice Compartmental Uptake Model (Trapp/Brunetti DPU의 "
               "이온성유기화합물 확장), 2026."),
    },
}

# organ labels for the custom-tables density inputs / volume readout (not a t() string)
ORGAN_LABELS = {
    "en": {"root": "root", "stem": "stem", "leaf": "leaf", "grain": "grain"},
    "ko": {"root": "뿌리", "stem": "줄기", "leaf": "잎", "grain": "낟알"},
}
