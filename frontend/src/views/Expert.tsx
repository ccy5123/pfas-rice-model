// Full research (English) view -- the port of ui/expert.py.

import { useState } from "react";
import { download, svgToPng } from "../api";
import { CustomTables } from "../components/CustomTables";
import { InversePanel } from "../components/InversePanel";
import { MeasuredPlantMap, PlantMap } from "../components/PlantMap";
import { ServerFigure } from "../components/Plot";
import { Alert, Button, Caption, Checkbox, Columns, Expander, Md, Metric, Radio, Select, Spinner, Tabs } from "../components/ui";
import { fmtF, fmtG } from "../format";
import { t } from "../i18n";
import { plantSvgFromRun } from "../plantSvg";
import type { Meta, RunResult, Scenario } from "../types";
import type { ExpertState } from "./ExpertSidebar";
import { modeLabel } from "./ExpertSidebar";

const ALL_TABS = ["🗺️ Plant & soil map", "📈 Tissue dynamics", "🟫 Soil & drivers", "📊 BAF vs observed",
                  "🔗 Chain-length trends", "⚖️ Compare congeners", "✅ Tang TF (OOS)", "🔎 Inverse (Bayesian)",
                  "ℹ️ About / coupling"];
const PFAS_ONLY = new Set(["📊 BAF vs observed", "🔗 Chain-length trends", "⚖️ Compare congeners", "✅ Tang TF (OOS)"]);

export function ExpertHeader({ meta, neutral }: { meta: Meta; neutral: boolean }) {
  return (
    <>
      <h1>🌾 PFAS in Rice — Uptake Explorer</h1>
      <Alert kind="warning">{t(meta, "disclaimer", "en")}</Alert>
      <Caption>{neutral
        ? "**Neutral-organic mode** — the same 4-compartment dynamic plant-uptake model run as the Briggs/Kow base (`z = 0`: no anion exclusion, no carrier, nothing fitted; `K_PW` and `TSCF` follow from log Kow). Switch the compound class in the sidebar to return to PFAS."
        : t(meta, "header.expert_caption", "en")}</Caption>
    </>
  );
}

const BIOMON_NOTE = "⚠ This chart shows a **model reference run** (Cwᵒ=1), **not** your measured biomonitoring data. Your measured values are on the **map** and **BAF** tabs.";

function TangTab({ meta, scn, run }: { meta: Meta; scn: Scenario; run: RunResult }) {
  const [dose, setDose] = useState<"across-dose mean" | "0.1 µg/g (env-closest)">("across-dose mean");
  const [oryza, setOryza] = useState(false);
  const [lipid, setLipid] = useState(false);
  const [fx, setFx] = useState<{ f_xy: number; f_xy_refit: number } | null>(null);
  const c = run.congener;
  if (scn.compound.kind !== "congener" || !meta.tang_congeners.includes(c)) {
    return <Alert kind="info">{`Tang 2026 covers **${meta.tang_congeners.join(", ")}** only — pick one of these curated congeners (sidebar) to see the comparison.`}</Alert>;
  }
  const verdict: Record<string, string> = {
    GenX: "GenX over-predicted ~12× by the provisional 0.233; the refit ≈ the documented 0.013.",
    PFOS: "PFOS is **dataset-dependent** — Yamazaki W2 0.142 vs Tang ~0.32; report the **range with conditions**, not one value.",
    PFOA: "PFOA f_xy is dose-sensitive (0.064 across-dose mean → 0.097 at 0.1 µg/g).",
  };
  return (
    <>
      <div className="cols" style={{ gridTemplateColumns: "2fr 1fr 1fr" }}>
        <Radio label="Tang dose" horizontal options={["across-dose mean", "0.1 µg/g (env-closest)"] as const} value={dose} onChange={setDose} />
        <Checkbox label="ORYZA biomass driver" checked={oryza} onChange={setOryza}
                  help="Drive the shoot model with the mechanistic ORYZA2000 biomass instead of the logistic (slower)." />
        <Checkbox label="lipid loading (OOS)" checked={lipid} onChange={setLipid}
                  help="Turn on the K_PL-gated lipid-facilitated loading. Its constants were fit on YAMAZAKI and nothing here is re-fit to Tang, so this is a genuine out-of-sample prediction — and it is the one that works: Tang per-organ log10 RMSE 1.232 (free anion) → 0.516, matching the in-sample Tang refit (0.519). See validation/oos_tang_lipid.py." />
      </div>
      <ServerFigure scenario={scn} name="tang_tf" lang="en"
                    opts={{ dose: dose.startsWith("0.1") ? "low" : "mean", biomass: oryza ? "oryza" : "growth_rice", lipid }}
                    onData={(d) => setFx(d as any)} />
      {lipid && (
        <Alert kind="success">{"**The purple bars are a prediction, not a fit.** The green refit bars were calibrated ON these measurements; the lipid mechanism was calibrated on Yamazaki and transferred here untouched, yet lands in the same place (0.516 vs 0.519 across Tang's three congeners). The dominant free-anion failure — PFOS, the high-K_PL sulfonate, ~40–200× under — is fixed at the mechanism level (stalk 0.013 → 0.620 vs Tang 0.571). Residuals it does NOT fix: GenX (provisional ether f_xy) and the PFOS endosperm (~5× under)."}</Alert>
      )}
      {fx && (
        <Md>{`- **f_xy**: recommended **${fmtG(fx.f_xy, 3)}** → Tang-refit **${fmtG(fx.f_xy_refit, 3)}**. ${verdict[c]}\n` +
             "- **Dry-weight basis** — `TF_dw = TF_fw·(1−θ_root)/(1−θ_tissue)`; comparing the model's fresh-weight TF to Tang's dry-weight TF (without this factor) is a units error that flatters the grain ~8×.\n" +
             "- **Grain/endosperm is structurally under-predicted ~3–8×** and is **not** closable by `L_Ph`/lipid — a phloem-delivery + dry-weight limit, not a calibration gap."}</Md>
      )}
    </>
  );
}

function BafTab({ scn, run }: { scn: Scenario; run: RunResult }) {
  const [tp, setTp] = useState(true);
  const [mg, setMg] = useState(false);
  const curated = scn.compound.kind === "congener";
  if (run.bio_baf) {
    return (
      <>
        <ServerFigure scenario={scn} name="biomon_compare" lang="en" />
        <Caption>Measured biomonitoring BAF (tissue conc ÷ pore water) vs the model prediction.</Caption>
      </>
    );
  }
  const overlays = curated && (tp || mg);
  return (
    <>
      {curated && (
        <Columns n={2}>
          <Checkbox label="Overlay the two-pool (seq) model — EXPLORATORY" checked={tp} onChange={setTp}
                    help="The sequestration two-pool root model (model_api.simulate_twopool_seq; docs/twopool_root_exploration.md): a mobile pool + an irreversible non-K_PL k_seq sink. Shown at its calibrated operating point (Cwo=1, season≈120, demo forcings) so it is comparable to the fixed Yamazaki bars. In-sample / opt-in; parameters.json and the canonical core are unchanged." />
          <Checkbox label="…and the MERGED root+shoot model" checked={mg} onChange={setMg}
                    help="The structural merge (model_api.simulate_twopool_nstem): the same two-pool sequestration ROOT driving the redistributed N-stem+leaf SHOOT, instead of the 4-pool pass-through stem. Swapping the whole shoot costs only ~0.04 log on Yamazaki (root RMSE unchanged) — the root mechanism is separable — and it fixes the two-pool's Tang stalk collapse (per-organ OOS 1.398→0.801). EXPLORATORY, same calibrated operating point." />
        </Columns>
      )}
      {scn.mechanism.lipid_loading && (
        <Alert kind="info">{"**Lipid-facilitated loading is ON**, so the blue 'core' bar is the K_PL-gated mechanism, not the shipped free-anion model. It is the repo's strongest out-of-sample result elsewhere (Tang, Kim), but against these Yamazaki bars it is *in-sample* — its constants were fit on them."}</Alert>
      )}
      <ServerFigure scenario={scn} name="baf" lang="en" opts={{ twopool_seq: curated && tp, twopool_nstem: curated && mg }} />
      {overlays && (
        <Caption>{"🧪 **two-pool (seq)** is EXPLORATORY / in-sample (Yamazaki fit): it captures the long-chain **root** BAF and the PFOS/PFUnDA split the single-pool core misses, while keeping the monotone physical f_xy — overall log10 RMSE 0.251. Run at its calibrated point, so it does **not** track the sidebar f_xy/Cwᵒ/biomass (unlike the 4-pool core bar). The **merged** overlay is that same root driving the redistributed shoot instead of the pass-through stem (Yamazaki 0.301 vs 0.278 — the shoot swap costs ~0.04 log, so the root mechanism is separable). The *carrier* two-pool (`close_longchain_2pool`, saturated DOF-0 closure) is API-only — it reproduces the observed bars by construction and is too slow (~1 min) to render live."}</Caption>
      )}
      {Object.keys(run.obs).length === 0 ? (
        <Alert kind="info">No Yamazaki BAF for this congener (model prediction only).</Alert>
      ) : (
        <Expander title="ℹ️ What are the Yamazaki 2023 conditions? (and how to match them)">
          <Md>{"**Yamazaki et al. 2023**, *Environ. Sci. Technol.* **57** ([doi:10.1021/acs.est.2c08767](https://doi.org/10.1021/acs.est.2c08767)).\n\n" +
               "- **Design**: greenhouse **pot** study, Japanese **Andosol** soil; each congener spiked **individually** with clean irrigation water; grown a **full cycle** to maturity.\n" +
               "- **Cultivars**: the plotted observed BAFs are the **geomean of Indica + Japonica** (SI tables S16/S18/S19), fresh-weight.\n" +
               "- **BAF definition**: **tissue conc ÷ pore-water (soil-solution) conc** [L/kg], for **root / straw (stem+leaf) / grain (brown rice)** — the same definition the model reports.\n\n" +
               "**These observed points are fixed measurements** — they do **not** move when you change the sidebar. The model was calibrated to reproduce them at one operating point, so the overlay is a like-for-like comparison **only at those settings**:"}</Md>
          <Md>{"| sidebar control | match-Yamazaki value |\n|---|---|\n| Data source | **Model (parametric)**, *Measured forcings* on |\n| Pore-water Cwᵒ | **1.0 µg/L** (so tissue conc = BAF) |\n| Root→shoot f_xy | **W2 fit (reproduces Yamazaki)** |\n| E_m | **−120 mV** (default) · Season ~**120 d** |"}</Md>
          <Caption>Changing E_m / f_xy / Cwᵒ / season, or a non-parametric source (HYDRUS, soil inventory…), moves the model AWAY from the Yamazaki experiment — then the overlay is a reference trend, not a calibrated match. (Yamazaki = clean per-congener Andosol pot study; your scenario may differ in soil, exposure and crop.)</Caption>
        </Expander>
      )}
    </>
  );
}

const ABOUT_MD = `
### Six ways to drive the plant model

| Mode | PFAS driver Cwᵒ(t) comes from | Q_TP, M(t) | When to use |
|---|---|---|---|
| **Model (parametric)** | a constant you set (flat, or the analytic flooded shape) | measured FAO-56 / ORYZA | quick what-ifs, teaching |
| **Custom tables** | a Cwᵒ(t) table you type or paste | your own growth table, or measured | you have field numbers but no soil model |
| **HYDRUS / CSV drivers** | a HYDRUS-1D / Phydrus run (CSV) | from the same CSV, or measured | you have a calibrated soil-water-solute model |
| **Run HYDRUS-1D (live)** | a real HYDRUS-1D run executed here | HYDRUS root uptake + ORYZA | you want the engine to run in-app (needs it built) |
| **Soil inventory** | inverting a soil load with a Freundlich isotherm | measured | you know total soil PFAS but not pore water |
| **Biomonitoring** | a measured pore-water / soil-solution value | — (not needed) | you have field tissue + water concentrations |

**Live HYDRUS-1D** runs the genuine engine (built from \`external/hydrus_source\`) via \`phydrus\`:
a one-season paddy model (Richards flow + advection–dispersion + linear Kd sorption + root
water uptake) gives the **congener-dependent** pore water — weakly-sorbed short chains leach
under flooding (Cwᵒ falls), strongly-sorbed long chains stay buffered (flat). Kd comes from the
Koc(chain-length) QSPR. The FORTRAN source is **vendored in the repo** (no submodule to
initialise): copy \`external/hydrus_source/makefile\` into \`source/\`, run \`make\` there (needs
gfortran), and \`pip install phydrus\`; the mode auto-detects the engine.

### HYDRUS-1D coupling (Method A, one-way) — inputs & outputs

**HYDRUS *inputs* you set up (soil side):** soil-hydraulic parameters (van Genuchten),
the atmospheric boundary condition (precip / irrigation / evaporation), the root
water-uptake distribution, solute-transport parameters (a linear K_d **or** the
Freundlich K_F, n; dispersivity), and the initial / boundary PFAS concentration.

**HYDRUS *outputs* this tool consumes** (map them into the driver CSV columns):

| CSV column | HYDRUS-1D source | meaning |
|---|---|---|
| \`t\` | output times | day after transplant |
| \`Cwo\` | \`Conc\` at the root-zone node (\`Obs_Node.out\` / \`solute1.out\`) | pore-water free anion [µg/L] |
| \`Qtp\` | \`vRoot\` (root water uptake) / \`T_act\` (\`T_Level.out\`) | transpiration stream [L/day] |
| \`M_root,M_stem,M_leaf,M_grain\` | a plant **growth** sub-model (not HYDRUS) | organ fresh mass [kg] |

The plant ODE is solved in Python; HYDRUS is not modified (the *tight* Method B —
embedding the root-uptake term \`j_R\` in the HYDRUS FORTRAN — is future work).

### Biomonitoring — when HYDRUS is unnecessary
If you already have **measured tissue concentrations** and a **measured pore-water
(or soil-solution) concentration**, the bioaccumulation factor is just
\`BAF = C_tissue / Cwᵒ\` — no transport simulation required. Use this mode to read
BAFs straight off field data and to overlay them on the model for a sanity check.

### Plain-language glossary
`;

export function ExpertView({ meta, st, set, scn, run }: {
  meta: Meta; st: ExpertState; set: (s: ExpertState) => void; scn: Scenario; run: RunResult;
}) {
  // `scn` is the DEBOUNCED scenario the run was computed for (figures/downloads use it);
  // `st` is the live sidebar state (edits made here go through it).
  const neutral = scn.compound.kind === "neutral";
  const names = neutral ? ALL_TABS.filter((n) => !PFAS_ONLY.has(n)) : ALL_TABS;
  const [tabName, setTabName] = useState(names[0]);
  const active = names.includes(tabName) ? tabName : names[0];
  const [chainKey, setChainKey] = useState<string>("K_PL");
  const [cmpTissue, setCmpTissue] = useState<"root" | "straw" | "grain">("straw");
  const [busy, setBusy] = useState(false);
  const p = run.params;
  const isBiomon = scn.exposure.mode === "biomon" && !!run.measured_bio;
  const mech = [...(scn.mechanism.uptake !== meta.default_uptake ? [`uptake: **${scn.mechanism.uptake}**`] : []),
                ...(scn.mechanism.lipid_loading ? ["**lipid loading ON**"] : [])];
  const desc = run.descriptors;
  const dl = async (kind: string) => { setBusy(true); try { await download(kind, scn, "en"); } finally { setBusy(false); } };
  const tissueMetric = (tis: "root" | "straw" | "grain") => {
    const pred = tis === "straw" ? run.straw_baf : run.baf_final[tis];
    let sub = tis in run.obs ? `obs ${fmtF(run.obs[tis], 2)}` : "no obs";
    if (run.bio_baf && tis in run.bio_baf) sub = `measured ${fmtF(run.bio_baf[tis], 2)}`;
    return <Metric key={tis} label={`${tis} BAF [L/kg]`} value={fmtF(pred, 2)} sub={sub} />;
  };

  return (
    <div>
      {scn.exposure.mode === "tables" && (
        <Expander open title={t(meta, "ct.panel_title", "en")}>
          <CustomTables meta={meta} lang="en" tables={st.scn.exposure.tables} run={run}
                        onChange={(tb) => set({ ...st, scn: { ...st.scn, exposure: { ...st.scn.exposure, tables: tb } } })}
                        season={scn.exposure.season} biomass={scn.mechanism.biomass} Cwo={scn.exposure.Cwo} />
        </Expander>
      )}
      <h2>{neutral ? `${run.congener}  (neutral organic · log Kow ${fmtF(p.log_kow, 2)})  ·  source: ${modeLabel(scn.exposure.mode)}`
                   : `${run.congener}  (C${p.n_C} ${p.group})  ·  source: ${modeLabel(scn.exposure.mode)}`}</h2>
      {mech.length > 0 && !neutral && (
        <Alert kind="warning">{`Non-default mechanism — ${mech.join(" · ")}. These are EXPLORATORY switches (sidebar → ⚙️ Mechanism); every number on this screen moves with them, and neither is the shipped calibration.`}</Alert>
      )}
      <Columns n={4}>
        {tissueMetric("root")}{tissueMetric("straw")}{tissueMetric("grain")}
        {neutral ? <Metric label="TSCF" value={fmtF(run.TSCF!, 3)} sub={`K_PW root ${fmtF(run.K_PW!.root, 2)}`} />
                 : <Metric label="anion exclusion eᴺ" value={fmtF(run.eN!, 0)} sub={`N=${fmtF(run.N!, 2)}`} />}
      </Columns>

      {neutral && (
        <>
          <Expander open title="🧪 Neutral organic — what this run is, and its scope">
            <Md>{"The framework's **Briggs/Kow base**, run on the *same* 4-compartment ODE with `z = 0`: the GHK factor → 1, anion exclusion `eᴺ` falls 107 → **1**, the membrane term degenerates exactly to passive diffusion, and the carrier is off. **Nothing here is fitted** — `K_PW` and `TSCF` both follow from log Kow via published QSPRs, so this is the one setting where the DPU *backbone* is exposed without the fitted PFAS transport behind it."}</Md>
            <Md>{`- a-priori vs measured rice: **log10 RMSE 0.281** root partition (Liu 2023) and **0.783** per-organ TF (Ge 2017) — on the season-ODE basis this run uses; on the equilibrium basis proper to a root-partition measurement Liu is **0.206**. Expect roughly a **2–6× error**, with the ordering far more reliable than the level.\n` +
                 `- \`TSCF\` = **${fmtF(run.TSCF!, 3)}** (${p.tscf_model} QSPR) · \`K_PW\` root **${fmtF(run.K_PW!.root, 2)}** / stem ${fmtF(run.K_PW!.stem, 2)} / leaf ${fmtF(run.K_PW!.leaf, 2)} L/kg · Briggs barley RCF at this log Kow **${fmtF(run.rcf_briggs!, 2)}** · lipid reading \`${p.lipid_source}\`` +
                 (run.half_life ? ` · half-life **${fmtG(run.half_life, 3)} d**` : " · **no metabolism** (γ=0)")}</Md>
            {!run.half_life && <Alert kind="warning">γ = 0: the leaf is an **unbounded terminal accumulator**, so the leaf and straw numbers are an UPPER BOUND, not a prediction. Set an in-planta half-life in the sidebar.</Alert>}
            <Caption>Scope: the **grain compartment is UNTESTED** for neutrals (no dataset reports a neutral Kow series in grain under root-only exposure) — read it as a structural output. Stem/leaf lipid are still Trapp 1994's soybean values (the root is corroborated by measured cereal roots). This is a RICE model: 120-day season, measured rice transpiration, ORYZA biomass — another crop needs its own composition and drivers, not just a log Kow.</Caption>
          </Expander>
          {run.warning && <Alert kind="warning">{run.warning}</Alert>}
          {run.air_summary && (() => {
            const hl = run.air_summary.leaf.half_life;
            return <Caption>{`Air exchange ON — leaf volatilisation half-life **${hl === null || hl === undefined ? "∞" : fmtG(hl, 3) + " d"}** (snapshot at final mass/transpiration), leaf P_plant ${fmtG(run.air_summary.leaf.P_plant, 3)} m/d. The flux scales with each tissue's specific SURFACE AREA, which this repo has only ever used as a leaf/grain *ratio* — treat absolute magnitudes as order-of-magnitude.`}</Caption>;
          })()}
        </>
      )}

      {desc && (
        <Expander open title={`🧬 Structure → parameters  ·  ${desc.matched_name ? "read-across (matches a curated congener)" : "QSPR (novel structure)"}`}>
          {run.provisional
            ? <Alert kind="warning">⚠ **PROVISIONAL** — novel / non-calibrated structure. Binding (K_PL/K_prot) is from the QSPR and translocation `f_xy` is a head-group estimate, so treat the magnitudes as indicative (ordering is more reliable than absolute level).</Alert>
            : <Alert kind="success">{`Matches **${desc.matched_name}** — uses the curated, measured-anchored parameters.`}</Alert>}
          <Columns n={2}>
            <Md>{`**Structure** \`${desc.formula}\`, MW ${fmtF(desc.mol_weight, 1)}\n\n- perfluoro-C: **${desc.n_perfluoroC}** · head group: **${desc.head_group}**\n- ether-O: ${desc.n_ether_O} · CF₃: ${desc.n_CF3} · branched: ${desc.branched ? "True" : "False"} · linear: ${desc.is_linear ? "True" : "False"}\n- read-across match: **${desc.matched_name ?? "none (novel)"}**`}</Md>
            <Md>{`**Model parameters**\n\n- K_PL = ${fmtF(p.K_PL, 0)} · K_prot = ${fmtF(p.K_prot, 0)} · K_cw = ${fmtF(p.K_cw, 0)}  [L/kg]\n- f_xy = ${fmtG(p.f_xy, 4)} · L_Ph = ${fmtG(p.L_Ph, 3)} · κ_d = ${fmtG(p.kappa_d, 3)}\n- B_k: ` +
                 Object.entries(run.B_k ?? {}).map(([k, v]) => `${k} ${fmtF(v, 1)}`).join(", ")}</Md>
          </Columns>
          {desc.notes.map((n) => <Caption key={n}>{"• " + n}</Caption>)}
        </Expander>
      )}

      <Tabs names={names} active={names.indexOf(active)} onChange={(i) => setTabName(names[i])} />
      <div className="tab-body">
        {active === "🗺️ Plant & soil map" && (isBiomon ? (
          <>
            <MeasuredPlantMap conc={run.measured_bio!.conc} cwo={run.measured_bio!.Cwo} lang="en" />
            <Caption>Compartments coloured by the MEASURED concentration [µg/kg]; the tag shows the pore-water Cwᵒ [µg/L]. Stem/leaf share the straw colour when only straw is reported.</Caption>
          </>
        ) : (
          <>
            <PlantMap run={run} lang="en" showMetric />
            <Caption>Each compartment is filled by its accumulation on a shared colour scale — drag the day slider to watch where PFAS builds up (numbers in µg/kg, or L/kg for BAF; the tag is pore-water µg/L). Leaf is xylem-terminal, grain is phloem-fed; the root retains the anion (low f_xy).</Caption>
          </>
        ))}
        {active === "📈 Tissue dynamics" && (
          <>
            {isBiomon && <Alert kind="warning">{BIOMON_NOTE}</Alert>}
            <ServerFigure scenario={scn} name="tissue" lang="en" />
            <ServerFigure scenario={scn} name="burden" lang="en" />
            {neutral ? (
              <Md>{"**K_PW [L/kg fw]** — " + Object.entries(run.K_PW ?? {}).map(([k, v]) => `${k}: ${fmtF(v, 2)}`).join(", ") +
                   `  ·  TSCF=${fmtG(run.TSCF!, 4)} (computed, not fitted), κ_d=${fmtG(p.kappa_d, 3)}, γ=${fmtG(p.gamma, 4)} 1/d`}</Md>
            ) : (
              <Md>{"**B_k [L/kg fw]** — " + Object.entries(run.B_k ?? {}).map(([k, v]) => `${k}: ${fmtF(v, 2)}`).join(", ") +
                   `  ·  f_xy=${fmtG(p.f_xy, 4)}, L_Ph=${fmtG(p.L_Ph, 4)}, κ_d=${fmtG(p.kappa_d, 3)}`}</Md>
            )}
            <Caption>Top: tissue **concentration** C_k(t) [µg/kg] (intensive). Bottom: **PFAS mass** (burden) = C_k(t)·M_k(t) [µg/hill] (extensive) — where the chemical actually ends up. A tissue can be high-concentration yet low-mass (small organ), so the two views differ; the terminal leaf/grain keep gaining mass as the organ grows. (Organ *biomass* M_k(t) is in the **🟫 Soil & drivers** tab.) The **grain takes up no PFAS until it forms (~flowering)** and then accumulates — its loading is *formation-gated* (the panicle is absent before then, so no solute enters it; DPU-consistent), and the empty pre-set period is not plotted.</Caption>
          </>
        )}
        {active === "🟫 Soil & drivers" && (
          <>
            {isBiomon && <Alert kind="warning">{BIOMON_NOTE}</Alert>}
            <ServerFigure scenario={scn} name="drivers" lang="en" />
            {run.has_isotherm && <ServerFigure scenario={scn} name="isotherm" lang="en" />}
            <ServerFigure scenario={scn} name="soil_profile" lang="en" />
            <Caption>Cwᵒ(t) is the only PFAS-specific driver; Q_TP(t) and M(t) are crop physiology (measured FAO-56 transpiration + the selected biomass driver) reused across modes.</Caption>
          </>
        )}
        {active === "📊 BAF vs observed" && <BafTab scn={scn} run={run} />}
        {active === "🔗 Chain-length trends" && (() => {
          const ring = meta.congeners.includes(run.congener) ? run.congener
            : desc?.matched_name && meta.congeners.includes(desc.matched_name) ? desc.matched_name : "PFOA";
          return (
            <>
              <Select label="Parameter" options={["K_PL", "K_prot", "K_cw_root", "f_xy_recommended", "B_root", "B_grain"]}
                      value={chainKey} onChange={setChainKey} />
              <ServerFigure scenario={scn} name="chain" lang="en" opts={{ key: chainKey }} />
              {!meta.congeners.includes(run.congener) && (
                <Caption>{`'${run.congener}' is a novel structure (not in the curated 13); the ring marks **${ring}** for reference. Its own parameters are in the 🧬 panel above.`}</Caption>
              )}
            </>
          );
        })()}
        {active === "⚖️ Compare congeners" && (
          <>
            <Radio label="Tissue" horizontal options={["root", "straw", "grain"] as const} value={cmpTissue} onChange={setCmpTissue} />
            {st.compare.length ? <ServerFigure scenario={scn} name="compare" lang="en" opts={{ congeners: st.compare, tissue: cmpTissue }} />
                               : <Alert kind="info">Select congeners in the sidebar to compare.</Alert>}
          </>
        )}
        {active === "✅ Tang TF (OOS)" && (
          <>
            <Md>**Tang 2026** (flooded paddy, Nipponbare, 150 d; PFOA/PFOS/GenX) — per-organ transfer factor **TF = C_organ/C_root**, an *out-of-sample* check of the root→shoot loading `f_xy` (only Tang's head-group *sign* went into the model build).</Md>
            <TangTab meta={meta} scn={scn} run={run} />
            <Caption>Source: docs/literature_db/raw_si/tang2026_doseresponse.csv (SI S7/S8). Model: redistributed-shoot `simulate_nstem_leaf`; the f_xy refit is OVERRIDE-only (parameters.json unchanged). Details: docs/VALIDATION_TANG2026_NSTEM_KR.md · docs/tang2026_grain_units_exploration.md.</Caption>
          </>
        )}
        {active === "🔎 Inverse (Bayesian)" && (
          <>
            <Md>**Bayesian inverse** — estimate the pore-water exposure Cwᵒ from measured tissue concentrations (Laplace posterior in log Cwᵒ; the well-posed direction of `validation/bayesian_inverse_demo.py`). Transport is held at the sidebar defaults.</Md>
            <InversePanel meta={meta} scenario={scn} lang="en"
                          enabled={neutral || (scn.compound.kind === "congener" && meta.congeners.includes(scn.compound.congener))} />
            <Caption>Identifiability: only the EXPOSURE level is estimated here. From tissue data alone Q_TP·f_xy is a product ridge and Cwᵒ vs root-uptake conductance is degenerate, so pinning transport absolutely needs an independent measurement (xylem sap / pore-water probe). See docs + validation/bayesian_inverse_demo.py.</Caption>
          </>
        )}
        {active === "ℹ️ About / coupling" && (
          <>
            <Md>{ABOUT_MD}</Md>
            <Md>{t(meta, "glossary", "en")}</Md>
            <Caption>Open item (recorded, not modelled): long-chain PFCA shoot accumulation (hysteretic sorption) — see docs/nstem_gradient_exploration.md.</Caption>
          </>
        )}
      </div>

      <Expander title="⬇️ Download results">
        <div className="downloads">
          <Button onClick={() => dl("summary_csv")} disabled={busy}>BAF summary (CSV)</Button>
          <Button onClick={() => dl("timeseries_csv")} disabled={busy}>Driver + tissue time series (CSV)</Button>
          <Button disabled={busy} onClick={() => svgToPng(plantSvgFromRun(run, -1, { lang: "en" }), `${run.congener}_map.png`)}>Plant map (PNG)</Button>
          <Button onClick={() => dl("map_html")} disabled={busy}>Plant map (interactive HTML)</Button>
          <Button onClick={() => dl("summary_html")} disabled={busy}>One-page summary (HTML)</Button>
        </div>
        {busy && <Spinner label="Preparing the file…" />}
      </Expander>
    </div>
  );
}
