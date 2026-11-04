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
// dropped (not blanked) for a neutral compound: statements about the PFAS series
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
  const [lipid, setLipid] = useState(false);
  const [fx, setFx] = useState<{ f_xy: number; rmse: number } | null>(null);
  const c = run.congener;
  if (scn.compound.kind !== "congener" || !meta.tang_congeners.includes(c)) {
    return <Alert kind="info">{`Tang 2026 covers **${meta.tang_congeners.join(", ")}** only — pick one of these curated congeners (sidebar) to see the comparison.`}</Alert>;
  }
  return (
    <>
      <div className="cols" style={{ gridTemplateColumns: "2fr 1fr" }}>
        <Radio label="Tang dose" horizontal options={["across-dose mean", "0.1 µg/g (env-closest)"] as const} value={dose} onChange={setDose} />
        <Checkbox label="lipid loading (amendment-1 default)" checked={lipid} onChange={setLipid}
                  help="Add the same shoot with the lipid loading instead of the per-congener f_xy. The grey bars are the default model, whose fitted constants come from YAMAZAKI with nothing fitted to Tang (a transfer, not a fit)." />
      </div>
      <ServerFigure scenario={scn} name="tang_tf" lang="en"
                    opts={{ dose: dose.startsWith("0.1") ? "low" : "mean", lipid }}
                    onData={(d) => setFx(d as any)} />
      {fx && (
        <Md>{`- **${c}, this dose**: log10 RMSE over the three organs **${fmtF(fx.rmse, 3)}** (computed now; nothing was fitted to Tang).\n` +
             "- **Dry-weight basis** — Tang's TF is dry weight; the model's is formed from the harvest operator's dry-basis concentrations (`conc_dw`)."}</Md>
      )}
    </>
  );
}

function BafTab({ meta, scn, run }: { meta: Meta; scn: Scenario; run: RunResult }) {
  const curated = scn.compound.kind === "congener";
  const hasObs = Object.keys(run.obs).length > 0;
  const kw = run.effective.sim_kw ?? {};
  if (run.bio_baf) {
    return (
      <>
        <ServerFigure scenario={scn} name="biomon_compare" lang="en" />
        <Caption>Measured biomonitoring BAF (tissue conc ÷ pore water) vs the model prediction.</Caption>
      </>
    );
  }
  void meta;
  return (
    <>
      {kw.lipid_loading ? (
        <Alert kind="info">{"**Lipid-facilitated loading is ON**, so against these Yamazaki bars the model is *in-sample* — its two lipid constants were fit on them. Its out-of-sample record is on Tang, Kim and Li."}</Alert>
      ) : kw.f_xy_source === "fit" ? (
        <Alert kind="info">{"**f_xy 'fit'** (the default): f_xy, L_Ph and r_eq were fitted per congener to these Yamazaki bars, so here the model is a *reproduction*. Its out-of-sample record is on Tang, Kim and Li."}</Alert>
      ) : null}
      {hasObs && curated ? (
        <>
          <ServerFigure scenario={scn} name="baf_observable" lang="en" opts={{ view: "expert" }} />
          <Caption>Model bar = the Yamazaki 2020 site-year run the way the data were taken (dry weight; straw with attached dead blades; two harvests combined as the two cultivars are; exposure = the congener's irrigation-water level). It follows the f_xy source and mechanism settings, not the sidebar exposure scenario.</Caption>
        </>
      ) : (
        <ServerFigure scenario={scn} name="baf" lang="en" />
      )}
      {!hasObs ? (
        <Alert kind="info">No Yamazaki BAF for this congener (model prediction only).</Alert>
      ) : (
        <Expander title="ℹ️ What are the Yamazaki 2023 conditions? (and how to match them)">
          <Md>{"**Yamazaki et al. 2023**, *Environ. Sci. Technol.* **57** ([doi:10.1021/acs.est.2c08767](https://doi.org/10.1021/acs.est.2c08767)).\n\n" +
               "- **Design**: an open **paddy FIELD** (NARO, Tsukuba; 10 × 50 m, Andosol), continuously flooded except the midsummer drainage, **ambient** exposure (irrigation water, rain, airborne particles) — **nothing was spiked**; grown a **full cycle**, harvested ~30 d after heading.\n" +
               "- **Cultivars**: the plotted observed BAFs are the **geomean of Kasalath (indica) + Koshihikari (japonica)** (SI tables S16/S18/S19), one composite sample per tissue.\n" +
               "- **BAF definition**: the paper computes **no BAF**. This repo builds it as tissue conc (**dry weight**, as reported) ÷ the **single irrigation-water grab** (no pore water was measured), for **root / straw / grain**. Straw = unweighted mean of the four stem+leaf height bands; grain = brown rice **reconstructed** as 0.88·white + 0.12·bran (brown rice itself was not analysed); PFBS grain is wholly imputed and not shown. See `docs/sampling_protocols_audit_KR.md`.\n\n" +
               "**The model bar is that same observation, computed**: the 2020 Tsukuba site-year (JMA weather, the ORYZA2000 port on the paper's transplanting / heading / harvest dates), two harvest runs combined as the two cultivars are, exposure = the congener's irrigation-water level, dry weight with the dead blades still on the straw. It follows the f_xy source and mechanism settings; the sidebar's exposure scenario does not enter it. With f_xy source **fit** it reproduces the bars by construction (3 fitted parameters per congener); otherwise the straw is a prediction."}</Md>
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
  // the effective transport (after the sidebar's fall-backs) against the default model
  const kw = run.effective.sim_kw ?? {};
  const mech = neutral ? [] : [
    ...(kw.lipid_loading === meta.default_lipid_loading ? [] : [`**lipid loading ${kw.lipid_loading ? "ON" : "OFF"}**`]),
    ...(kw.lipid_loading || kw.f_xy_source === meta.default_f_xy_source ? [] : [`f_xy: **${kw.f_xy_source}**`]),
    ...(kw.shoot === meta.default_shoot ? [] : [`shoot: **${kw.shoot}**`]),
  ];
  const desc = run.descriptors;
  const dl = async (kind: string) => { setBusy(true); try { await download(kind, scn, "en"); } finally { setBusy(false); } };
  // DRY basis on both sides: the Yamazaki observations are dry weight
  const dw = run.baf_final_dw;
  const tissueMetric = (tis: "root" | "straw" | "grain") => {
    const pred = dw ? dw[tis] : (tis === "straw" ? run.straw_baf : run.baf_final[tis]);
    let sub = tis in run.obs && dw ? `obs ${fmtF(run.obs[tis], 2)}` : "no obs";
    if (run.bio_baf && tis in run.bio_baf) sub = `measured ${fmtF(run.bio_baf[tis], 2)}`;
    return <Metric key={tis} label={`${tis} BAF [L/kg ${dw ? "dw" : "fw"}]`} value={fmtF(pred, 2)} sub={sub} />;
  };
  const kcwRoot = p.K_cw && typeof p.K_cw === "object" ? p.K_cw.root : p.K_cw;
  const ns = meta.neutral_scores, neq = meta.neutral_equilibrium_scores;
  const nsVals = Object.values(ns);
  const nsLo = nsVals.length ? Math.min(...nsVals) : NaN, nsHi = nsVals.length ? Math.max(...nsVals) : NaN;

  return (
    <div>
      {scn.exposure.mode === "tables" && (
        <Expander open title={t(meta, "ct.panel_title", "en")}>
          <CustomTables meta={meta} lang="en" tables={st.scn.exposure.tables} run={run}
                        onChange={(tb) => set({ ...st, scn: { ...st.scn, exposure: { ...st.scn.exposure, tables: tb } } })}
                        season={scn.exposure.season} crop={scn.mechanism.crop} Cwo={scn.exposure.Cwo} />
        </Expander>
      )}
      <h2>{neutral ? `${run.congener}  (neutral organic · log Kow ${fmtF(p.log_kow, 2)})  ·  source: ${modeLabel(scn.exposure.mode)}`
                   : `${run.congener}  (C${p.n_C} ${p.group})  ·  source: ${modeLabel(scn.exposure.mode)}`}</h2>
      {mech.length > 0 && (
        <Alert kind="warning">{`Non-default mechanism — ${mech.join(" · ")}. An alternative fitted variant (sidebar → ⚙️ Mechanism); every number on this screen moves with it, and it is not the default model.`}</Alert>
      )}
      {run.variant?.withdrawn && (
        <Alert kind="warning">{`Variant **${run.variant.key}** is withdrawn from scoring and decisions — ${run.variant.withdrawn}. It stays runnable and labelled (model_api.WITHDRAWN_VARIANTS).`}</Alert>
      )}
      <Columns n={4}>
        {tissueMetric("root")}{tissueMetric("straw")}{tissueMetric("grain")}
        {neutral ? <Metric label="TSCF" value={fmtF(run.TSCF!, 3)} sub={`K_PW root ${fmtF(run.K_PW!.root, 2)}`} />
                 : <Metric label="root entry r_eq" value={fmtG(p.r_eq, 3)} sub="fitted to Yamazaki's root" />}
      </Columns>
      {run.interp_spread && (
        <Caption>{"Growth table interpolated piecewise-linearly (shown). The monotone-cubic form moves the dry BAF by " +
          (["root", "straw", "grain"] as const).map((k) => `${k} ${run.interp_spread![k] >= 0 ? "+" : ""}${fmtF(run.interp_spread![k], 3)}`).join(", ") +
          " log10 — the spread from the unmeasured form between samples."}</Caption>
      )}

      {neutral && (
        <>
          <Expander open title="🧪 Neutral organic — what this run is, and its scope">
            <Md>{"The framework's **Briggs/Kow base**, run on the *same* 4-compartment ODE with `z = 0`: no charge, so the membrane term is exactly passive diffusion and there is no carrier. **Nothing here is fitted** — `K_PW` and `TSCF` both follow from log Kow via published QSPRs, so this is the one setting where the DPU *backbone* is exposed without the fitted PFAS transport behind it."}</Md>
            <Md>{(ns["Liu 2023 rice root"] !== undefined && ns["Ge 2017 per-organ TF"] !== undefined
                   ? `- a-priori vs measured rice: **log10 RMSE ${fmtF(ns["Liu 2023 rice root"], 3)}** root partition (Liu 2023) and **${fmtF(ns["Ge 2017 per-organ TF"], 3)}** per-organ TF (Ge 2017) — on the season-ODE basis this run uses; on the equilibrium basis proper to a root-partition measurement Liu is **${fmtF(neq["Liu 2023 rice root"], 3)}**. Across the ${nsVals.length} tables it has been scored on the error is log10 ${fmtF(nsLo, 2)}–${fmtF(nsHi, 2)} (**×${fmtF(10 ** nsLo, 1)}–${fmtF(10 ** nsHi, 1)}**), with the ordering more reliable than the level (validation/neutral_scores.py).\n`
                   : "- a-priori scores: run validation/neutral_scores.py.\n") +
                 `- \`TSCF\` = **${fmtF(run.TSCF!, 3)}** (${p.tscf_model} QSPR) · \`K_PW\` root **${fmtF(run.K_PW!.root, 2)}** / stem ${fmtF(run.K_PW!.stem, 2)} / leaf ${fmtF(run.K_PW!.leaf, 2)} L/kg · Briggs barley RCF at this log Kow **${fmtF(run.rcf_briggs!, 2)}** · lipid reading \`${p.lipid_source}\`` +
                 (run.half_life ? ` · half-life **${fmtG(run.half_life, 3)} d**` : " · **no metabolism** (γ=0)")}</Md>
            {!run.half_life && <Alert kind="warning">γ = 0: the leaf is an **unbounded terminal accumulator**, so the leaf and straw numbers are an UPPER BOUND, not a prediction. Set an in-planta half-life in the sidebar.</Alert>}
            {run.grain_predicted === false && <Alert kind="info">The neutral brown rice is **not predicted** in this run (grain feed off, or log Kow below 0 where neutrals are retained in the phloem; PREREG 9.15, 10.7c).</Alert>}
            <Caption>Scope: the **grain compartment is UNTESTED** for neutrals (no dataset reports a neutral Kow series in grain under root-only exposure) — read it as a structural output. Organ lipid is sourced on a dry basis (root: measured cereal roots; stem/leaf: rice straw; grain: brown rice). This is a RICE model: rice composition and rice crop drivers (by default the Tsukuba 2020 site-year) — another crop needs its own composition and drivers, not just a log Kow.</Caption>
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
            <Md>{`**Model parameters**\n\n- K_PL = ${fmtF(p.K_PL, 0)} · K_prot = ${fmtF(p.K_prot, 0)} · K_cw (root) = ${fmtF(kcwRoot, 0)}  [L/kg]\n- f_xy = ${fmtG(p.f_xy, 4)} · L_Ph = ${fmtG(p.L_Ph, 3)} · r_eq = ${fmtG(p.r_eq, 3)}\n- B_k: ` +
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
                   `  ·  f_xy=${fmtG(p.f_xy, 4)}, L_Ph=${fmtG(p.L_Ph, 4)}, root entry r_eq=${fmtG(p.r_eq, 3)}`}</Md>
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
            <Caption>Cwᵒ(t) is the only PFAS-specific driver; Q_TP(t) and M(t) are crop physiology from the selected crop drivers (default: the Tsukuba 2020 site-year — the ORYZA2000 port and its transpiration on JMA weather), reused across modes.</Caption>
          </>
        )}
        {active === "📊 BAF vs observed" && <BafTab meta={meta} scn={scn} run={run} />}
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
            <Caption>Source: docs/literature_db/raw_si/tang2026_doseresponse.csv (SI S7/S8). Model: the default configuration on the calibration site-year drivers (Tang's site and dates are not in the repo). Scoring script: validation/oos_tang_lipid.py.</Caption>
          </>
        )}
        {active === "🔎 Inverse (Bayesian)" && (
          <>
            <Md>**Bayesian inverse** — estimate the pore-water exposure Cwᵒ from measured tissue concentrations (Laplace posterior in log Cwᵒ; the well-posed direction of `validation/bayesian_inverse_demo.py`). Transport is held at the sidebar defaults.</Md>
            <InversePanel meta={meta} scenario={scn} lang="en" view="expert"
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
