# HANDOFF — App interface rework for a GENERAL audience

> Next-session work order. Read this first, then execute. The model/science is done;
> this is **UI/UX only** (`app.py`, `src/plots.py`, captions). No `parameters.json`,
> no `src/model_api.py` math changes unless a label needs a passthrough.

## Decision (from the user, this session)
- **Target audience = GENERAL / broad**: policy makers, undergraduates, the public —
  NOT domain experts. The current app is built for environmental-science researchers
  (dense jargon, 8 tabs, 5 modes, expert sliders). It must be **re-pitched for
  non-experts** while keeping the expert depth available behind a toggle.
- The user wants this done **next session** (this doc is the plan, no code yet).

## Guiding principle: progressive disclosure
Default to a **Simple mode** (plain language, one clear question → one clear answer,
minimal controls). Let experts opt into the current full UI via an **"Expert /
advanced"** toggle. Nothing is deleted — it's layered.

---

## Interface review findings (this session) — re-scoped for the general audience

### CORE (do first — these block a general-audience release)

1. **Onboarding / "what is this".** A first-time non-expert lands on the plant–soil
   heatmap (`tabs[0]`) with no framing. Add a short **plain-language intro** at the
   top (1–3 sentences) + a "Start here" hint:
   > *"PFAS are long-lasting 'forever chemicals'. This tool estimates how much PFAS a
   > rice plant takes up from contaminated paddy water/soil, and where it ends up
   > (roots, straw, grain). Pick a chemical and a contamination level on the left."*
   - Location: after `st.title`/`st.caption` (app.py:41–44), before the sidebar logic
     renders. Consider a dismissible `st.info` or an always-on intro card.

2. **Prominent research/educational DISCLAIMER (policy/public critical).** Today it's
   one word ("Outputs illustrative", app.py:44). For policy/public, make it a visible,
   unmissable banner:
   > *"Research & educational model — illustrative estimates only. NOT a regulatory,
   > food-safety, or health determination. Do not use for real exposure decisions."*
   - Location: top of page (a styled `st.warning`/`st.error` band) AND repeated in the
     footer. Keep it on every mode.

3. **Plain-language relabeling / glossary.** The UI is full of expert symbols. For the
   default (Simple) view, relabel; keep the symbol in parentheses or under Expert:
   | current (expert) | plain-language label (Simple mode) |
   |---|---|
   | `BAF [L/kg]` | "Concentration factor — how much it builds up vs water" |
   | `Cwᵒ` / "Pore-water Cwᵒ" | "PFAS in the soil water" |
   | `Q_TP` | "Plant water uptake" |
   | `f_xy` / "Root→shoot loading" | "How easily it moves up into the shoot" |
   | `E_m` (membrane potential) | hide in Simple (Expert only) |
   | "anion exclusion eᴺ" metric (app.py:407) | hide in Simple (Expert only) |
   | `B_k`, `κ_d`, `L_Ph`, `K_PL/K_prot/K_cw` | Expert only |
   | "formation-gated / DPU-consistent / monotone TSCF / W2 fit / log10 RMSE / two-pool seq" | strip from Simple captions; keep in Expert/About |
   - Add a small **glossary** (expander or the About tab, reworded in plain language).
   - Rewrite the dense captions (app.py:466–468, 476–482, 491–492, 525–531, 603–615)
     into one plain sentence each for Simple mode.

4. **Simplify the default control set (Expert toggle).** A non-expert does not need:
   SMILES input, `E_m`, `f_xy source`, `biomass driver`, `cwo_profile`/`k_leach`,
   "Run HYDRUS-1D (live)", "Soil inventory", the two-pool overlay, Chain-length and
   Tang-TF tabs. **Simple mode** should expose only:
   - compound: the curated congener dropdown (maybe a friendlier subset / names),
   - one "how contaminated?" input (the pore-water level, plain-labelled, maybe a
     low/medium/high preset instead of a µg/L number),
   - the **map**, **tissue dynamics**, and a plain **BAF-vs-observed** read-out.
   - Everything else → behind `st.toggle("Expert / advanced controls")` in the sidebar
     and extra tabs shown only when Expert is on.

5. **Biomonitoring-mode confusion (carry-over bug).** In Biomonitoring mode the map
   uses the measured data, but **Tissue dynamics (`tabs[1]`) and Soil & drivers
   (`tabs[2]`) still show an unrelated default model run** (drivers=None → Cwo=1).
   Either hide those tabs in Biomonitoring mode or clearly mark them "model reference,
   not your measured data". (For Simple mode, Biomonitoring may be Expert-only anyway.)

### SHOULD (polish for release)

6. **Export / download.** Policy users want artifacts. Add: a **"Download results
   (CSV)"** button (the BAF table + the driver/tissue series) and, where kaleido is
   present, a figure PNG export. Graceful when kaleido absent.

7. **Consistency.**
   - `season` slider only exists in parametric mode; other modes hardcode 120
     (app.py:310,327). Either expose it everywhere or drop it from the Simple view.
   - **Two different `k_leach` sliders**: parametric flooded (0–0.15, per-congener
     calibrated default, app.py:231) vs soil-inventory (0–0.1, flat 0.02, app.py:330).
     Unify the range/default (the soil-inventory one should also use
     `api.default_k_leach`). Expert-only anyway.

8. **Footer / metadata.** Add a footer with: app version, repo + docs links, a
   "How to cite" line, and the disclaimer (again). None exist today.

### NICE (later)

9. Tab count (8) is high and wraps on narrow screens — in Simple mode show ~3 tabs;
   group the niche ones (Chain-length, Tang-TF, Compare) under Expert.
10. Colorblind-safe colormap option (currently YlOrRd only, `plots._HEAT`).
11. Mobile/narrow reflow check (4-col metrics app.py:400, multi-col controls).
12. Friendlier congener names (e.g. "PFOA (C8 carboxylic acid)") in the dropdown.

---

## Suggested execution order (next session)
1. Add the **Expert toggle** scaffolding in the sidebar + gate the advanced controls/
   tabs on it (biggest structural change; do first).
2. Add the **intro card + disclaimer banner + footer** (cheap, high impact).
3. **Plain-language relabel** the Simple-mode controls, metrics, and captions; add a
   glossary expander.
4. Fix the **Biomonitoring-mode** non-map tabs (hide or mark).
5. Add **CSV/PNG export**.
6. Consistency cleanups (`season`, the two `k_leach` sliders).
7. Re-screenshot (headless Streamlit + Playwright; Chromium is pre-installed at
   `/opt/pw-browsers/chromium`) to verify the Simple-mode landing reads cleanly to a
   non-expert, and that Expert mode still has everything.

## Acceptance criteria
- A non-expert can land, understand what the tool does in <15 s, pick a chemical +
  contamination level, and read a plain-language result — **without meeting a single
  unexplained symbol** (BAF/Cwᵒ/f_xy/eᴺ/B_k…) in the default view.
- The research/educational disclaimer is visible on every screen.
- Expert mode restores 100% of today's controls/tabs (nothing lost).
- `tests/test_model_api.py` / `test_plots.py` still pass (UI-agnostic core untouched);
  add `test_plots` coverage for any new plain-language plot builder.
- `parameters.json`, `reproduce_demo` (RMSE 0.029), and the model math are UNCHANGED.

## Resume prompt for the next session
> "Rework `app.py` for a GENERAL audience per `docs/HANDOFF_app_general_audience.md`:
> add a Simple/Expert toggle (Simple = default), an intro card + a prominent
> research/educational disclaimer + a footer, plain-language relabels and a glossary,
> fix the Biomonitoring-mode non-map tabs, and add CSV/PNG export. Keep all expert
> controls behind the toggle; do not touch the model math or parameters.json. Verify
> with a headless Streamlit + Playwright screenshot of the Simple landing."

## Notes / pointers
- All compute is UI-agnostic in `src/model_api.py`; plot builders in `src/plots.py`
  (head-less tested). Simple mode should reuse them — relabeling is at the `app.py` /
  caption layer, plus maybe a couple of plain-language `plots.fig_*` wrappers.
- Deploy: Streamlit Cloud from `main` (see `docs/deploy.md`); reboot to pick up changes.
- This handoff complements `docs/HANDOFF_BAF_twopool.md` (science handoff); this one is
  purely the interface/audience rework.
