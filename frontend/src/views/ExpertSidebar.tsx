// Expert sidebar -- the port of the Expert branch of ui/sidebar.py (+ its neutral,
// CompTox and mechanism panels). It edits one `Scenario` object; the App debounces it
// and asks the server for a run only when it actually changes.

import { useEffect, useState } from "react";
import { buildHydrus, ctxLookup, getCwoPreview, getInventoryDefaults, getMolSvg } from "../api";
import { Plot } from "../components/Plot";
import {
  Alert, Button, Caption, Checkbox, Columns, Expander, FileText, MultiSelect, NumberField, Radio, Select,
  Slider, Spinner, TextField,
} from "../components/ui";
import { fmtF, fmtG } from "../format";
import { useAsync } from "../hooks";
import { congenerLabel } from "../i18n";
import type { Crop, CtxResult, Meta, Mode, NeutralSpec, RunResult, Scenario } from "../types";

export const MODES: [Mode, string][] = [
  ["parametric", "Model (parametric)"], ["tables", "Custom tables (Cwᵒ + growth)"],
  ["csv", "HYDRUS / CSV drivers"], ["hydrus", "Run HYDRUS-1D (live)"],
  ["inventory", "Soil inventory → pore water"], ["biomon", "Biomonitoring (measured tissue)"],
];
export const modeLabel = (m: Mode) => MODES.find(([k]) => k === m)?.[1] ?? m;
const CLASSES = ["Curated congener", "SMILES (structure)", "Neutral organic (log Kow)"] as const;
const classOf = (k: string) => (k === "smiles" ? CLASSES[1] : k === "neutral" ? CLASSES[2] : CLASSES[0]);
const FXY = ["fit (reproduces Yamazaki)", "recommended (monotone, a priori)"] as const;
const SHOOT = ["redistributed shoot (default)", "single-straw 4-pool"] as const;
const SHAPES = ["constant (flat)", "flooded (the paddy's HYDRUS shape)"] as const;

export interface ExpertState {
  scn: Scenario;
  compare: string[];
  smilesExample: string;
  biomonInput: "Manual" | "Upload CSV (tissue,conc[,Cwo])";
  csvName: string | null;
  biomonCsvName: string | null;
  lookup: { query: string; result: CtxResult | null; notes: string[]; koc: number | null;
            pka: number | null; pkaMeasured: boolean } | null;
}

/** Henderson–Hasselbalch (literature_params.speciation): [f_n, f_d]. */
export function speciation(pKa: number, pH: number, isAcid: boolean): [number, number] {
  const fd = isAcid ? 1 / (1 + Math.pow(10, pKa - pH)) : 1 / (1 + Math.pow(10, pH - pKa));
  return [1 - fd, fd];
}
export const karickhoff = (meta: Meta, logKow: number) =>
  Math.pow(10, meta.karickhoff.intercept + meta.karickhoff.slope * logKow);

/** Clamp a looked-up DEFAULT into a field's range and SAY so (ui/common.clamp_seed). */
function clampSeed(v: number, lo: number, hi: number): [number, string | null] {
  if (!Number.isFinite(v)) return [lo, "came back as NaN"];
  if (v < lo) return [lo, `came back as ${fmtG(v)}, below this field's minimum ${fmtG(lo)}`];
  if (v > hi) return [hi, `came back as ${fmtG(v)}, above this field's maximum ${fmtG(hi)}`];
  return [v, null];
}

/** The quoted a-priori neutral scores (Liu, Ge), when the record carries both. */
function quotedNeutral(meta: Meta): string {
  const s = meta.neutral_scores;
  const liu = s["Liu 2023 rice root"], ge = s["Ge 2017 per-organ TF"];
  return liu !== undefined && ge !== undefined ? ` (Liu ${fmtF(liu, 3)}, Ge ${fmtF(ge, 3)})` : "";
}

// ---------------------------------------------------------------- CompTox lookup
function LookupPanel({ meta, st, set }: { meta: Meta; st: ExpertState; set: (s: ExpertState) => void }) {
  const [q, setQ] = useState(st.lookup?.query ?? "");
  const [busy, setBusy] = useState(false);
  const r = st.lookup?.result;
  const LABEL: Record<string, string> = { log_kow: "log Kow", K_AW: "K_AW", koc: "soil Koc",
    pka_acidic: "pKa (acidic)", pka_basic: "pKa (basic)",
    pka_unlabelled: "pKa (centre not stated — check acid/base)" };
  const go = async () => {
    const query = q.trim();
    if (!query) return;
    setBusy(true);
    try {
      const res = await ctxLookup(query);
      const n: NeutralSpec = { ...st.scn.compound.neutral };
      const notes: string[] = [];
      let koc: number | null = null, pka: number | null = null;
      if (res.ok) {
        const kw = res.kwargs;
        const seed = (field: string, v: number, lo: number, hi: number) => {
          const [c, note] = clampSeed(v, lo, hi);
          if (note) notes.push(`**${field}** ${note} — clamped to ${fmtG(c)}.`);
          return c;
        };
        if (kw.log_kow !== undefined) n.log_kow = seed("log Kow", kw.log_kow, -2, 8);
        n.name = res.name || kw.name || query;
        if (kw.MW !== undefined) n.MW = seed("MW", kw.MW, 1, 1e5);
        if (kw.K_AW !== undefined) n.K_AW = seed("K_AW", kw.K_AW, 0, 10);
        if (kw.Koc !== undefined) { koc = seed("soil Koc", kw.Koc, 0, meta.koc_max); n.Koc = koc; }
        else n.Koc = null;
        if (kw.pKa !== undefined) {
          pka = seed("pKa", kw.pKa, -5, 14);
          n.pKa = pka;
          n.is_acid = kw.is_acid !== false;
          // only a MEASURED pKa switches the model path by itself
          n.weak_electrolyte = res.pka_measured;
        }
      }
      set({ ...st, scn: { ...st.scn, compound: { ...st.scn.compound, neutral: n } },
            lookup: { query, result: res, notes, koc, pka, pkaMeasured: res.pka_measured } });
    } finally { setBusy(false); }
  };
  const predicted: string[] = [], outAd: string[] = [];
  const rows: string[] = [];
  if (r?.ok) {
    for (const [pk, badge] of Object.entries(r.badges)) {
      if (["henry", "water_solubility", "melting_point"].includes(pk)) continue;
      rows.push(`- **${LABEL[pk] ?? pk}** — ${badge}`);
      if (r.sources[pk] !== "experimental") predicted.push(LABEL[pk] ?? pk);
      if (r.outside_ad[pk]) outAd.push(LABEL[pk] ?? pk);
    }
    if (r.kwargs.MW) rows.push(`- **MW** ${fmtF(r.kwargs.MW, 2)} g/mol`);
  }
  const pkaKw = r?.ok ? r.kwargs.pKa : undefined;
  return (
    <Expander title="🔎 Look up the compound (EPA CompTox)">
      {!meta.available.ctx && (
        <Caption>No EPA CompTox key on the server (`CTX_API_KEY`) — the lookup will report that, and every field below stays manual.</Caption>
      )}
      <div className="row">
        <input type="text" value={q} placeholder="carbamazepine  ·  298-46-4  ·  NC(=O)N1c2ccccc2…"
               onChange={(e) => setQ(e.target.value)} onKeyDown={(e) => { if (e.key === "Enter") go(); }} />
        <Button onClick={go} disabled={busy || !q.trim()}>Look up</Button>
      </div>
      {busy && <Spinner label="Looking the compound up on CompTox…" />}
      {r && !r.ok && <Alert kind="warning">{`No fill: ${r.note}`}</Alert>}
      {r?.ok && (
        <>
          <Alert kind="success">{`**${r.name || st.lookup?.query}** · ${r.dtxsid}${r.casrn ? ` · CAS ${r.casrn}` : ""}`}</Alert>
          {st.lookup?.notes.map((n) => <Alert key={n} kind="warning">{n}</Alert>)}
          <Caption>{rows.length ? rows.join("\n") : "_no physicochemical rows returned_"}</Caption>
          {outAd.length > 0 && (
            <Alert kind="error">{`**Outside its own applicability domain: ${outAd.join(", ")}.** OPERA reports these as beyond the chemical space the model was trained on — that is not 'uncertain', it is out of scope. Replace them by hand or leave the corresponding term off.`}</Alert>
          )}
          {predicted.length > 0 && (
            <Alert kind="warning">{`**Predicted (not measured): ${predicted.join(", ")}.** Filled in anyway and marked — but this path's published a-priori errors${quotedNeutral(meta)} are on MEASURED log Kow, so replace it with a measured value when you have one before quoting them for this run.`}</Alert>
          )}
          {pkaKw !== undefined && (() => {
            const [fn] = speciation(pkaKw, 6.5, r.kwargs.is_acid !== false);
            const head = `A pKa (${fmtG(pkaKw)}, ${r.kwargs.is_acid !== false ? "acid" : "base"}) came back — at the default root-zone pH 6.5 that is **f_n = ${fmtG(fn, 3)}**.`;
            return r.pka_measured
              ? <Alert kind="info">{head + " It is MEASURED, so the run is now on the **weak-electrolyte** path rather than the strictly neutral one. Untick the ⚗️ box to force the strictly neutral path."}</Alert>
              : <Alert kind="warning">{head + " It is **PREDICTED**, and a pKa is the one field that changes WHICH MODEL RUNS — so it is filled in but **left switched off**: the run stays on the strictly neutral path until you tick the ⚗️ box yourself. **Carbamazepine is why**: CompTox's OPERA gives an acidic pKa of **5.07** (→ f_n 0.036, 96% ionised) while the measured value this repo's own best-conditioned table (Kodešová 2019, §4f) is built on is **13.9** — un-ionised everywhere (f_n = 1.00). Nine log units apart, and only the measured one reproduces the published a-priori result. Check it against a source, then switch it on."}</Alert>;
          })()}
          <Caption>The **in-planta half-life is never filled** — it is not a dashboard property, and Kodesova 2019 measured the surviving parent fraction varying 4.8× BETWEEN SPECIES for one compound, so it is not a compound constant. Set it yourself.</Caption>
        </>
      )}
    </Expander>
  );
}

// ---------------------------------------------------------------- neutral organic
function NeutralPanel({ meta, st, set }: { meta: Meta; st: ExpertState; set: (s: ExpertState) => void }) {
  const n = st.scn.compound.neutral;
  const setN = (patch: Partial<NeutralSpec>) =>
    set({ ...st, scn: { ...st.scn, compound: { ...st.scn.compound, neutral: { ...n, ...patch } } } });
  const kar = karickhoff(meta, n.log_kow);
  const kocShown = n.Koc ?? kar;
  const [fn, fd] = speciation(n.pKa, n.pH, n.is_acid);
  const found = st.lookup?.result?.ok ? st.lookup : null;
  const lip = meta.neutral_lipid_fw;
  return (
    <>
      <LookupPanel meta={meta} st={st} set={set} />
      <NumberField label="log Kow" value={n.log_kow} min={-2} max={8} step={0.05}
                   onChange={(v) => setN({ log_kow: v, Koc: found?.koc ?? null })}
                   help="The one required input. Drives K_PW and TSCF. The measured tables behind this path span about -0.7 to 5.4; outside that you are extrapolating." />
      <TextField label="compound name" value={n.name} commitOnBlur onChange={(v) => setN({ name: v || "neutral" })} />
      <NumberField label="in-planta half-life [d]  (0 = none)" value={n.half_life} min={0} max={365} step={1}
                   onChange={(v) => setN({ half_life: v })}
                   help="STRONGLY recommended. With no metabolism the leaf is an unbounded terminal accumulator, so a run without it is an UPPER BOUND. Do not treat one number as a compound constant: Kodesova measured the parent fraction varying 4.8x BETWEEN SPECIES for the same compound, and soil persistence is not plant persistence." />
      <Radio label="TSCF QSPR" horizontal options={["briggs", "schriever"] as const} value={n.tscf_model}
             onChange={(v) => setN({ tscf_model: v })}
             help="Briggs 1982 (narrow bell, peak 0.784 @ logKow 1.78) vs the Schriever 2020 refit (97 values, ~3x broader). TSCF is an INPUT here, not a fitted parameter, so the gap between them measures how well it is known." />
      <Expander title="🧪 Composition, phloem, air">
        <Radio label="root lipid reading" horizontal options={["measured", "briggs_anchor"] as const}
               value={n.lipid_source} onChange={(v) => setN({ lipid_source: v })}
               help={`'measured' (default) = measured cereal-root lipid (Li 2019), carried on a dry basis and restated at the rice root's water: ${fmtF(lip.measured * 100, 2)}% fresh weight; 'briggs_anchor' = ${fmtF(lip.briggs_anchor * 100, 2)}%, what Briggs' 1982 regression implies — he measured no lipid at all.`} />
        <Checkbox label="phloem ON (departs from the base)" checked={n.phloem} onChange={(v) => setN({ phloem: v })}
                  help="The neutral base excludes phloem transport — it is an addition of the ionisable extension. Turning it on drives the small terminal grain hard: a statement about assuming unrestricted loading, not a prediction." />
        <Checkbox label="brown rice fed by the xylem sap (model default)" checked={n.grain_feed === "xylem"}
                  onChange={(v) => setN({ grain_feed: v ? "xylem" : null })}
                  help="simulate_neutral's grain_feed: the neutral BROWN RICE takes a sap stream at the xylem-sap concentration, T_C,Ph × its dry growth, out of the leaf's xylem inflow (Trapp 2007; Rein 2011; PREREG 9.15). Off = the brown rice is not predicted. Below log Kow 0 it is not predicted either way (neutrals retained in the phloem, PREREG 10.7c)." />
        <Checkbox label="plant–air exchange (volatilisation + gaseous uptake)" checked={n.air}
                  onChange={(v) => setN({ air: v })}
                  help="Off by default: needs K_AW and a molar mass, which the strict Kow-only a-priori run does not use. Identically zero at K_AW = 0. Leaving it off for a volatile compound makes the leaf an upper bound BY CONSTRUCTION." />
        {n.air && (
          <>
            <Columns n={2}>
              <NumberField label="MW [g/mol]" value={n.MW} min={1} max={1e5} step={1} onChange={(v) => setN({ MW: v })} />
              <NumberField label="K_AW [-]" value={n.K_AW} min={0} max={10} step={1e-4}
                           format={(v) => v.toExponential(2)} onChange={(v) => setN({ K_AW: v })}
                           help="Dimensionless Henry's-law constant (converted from the Henry constant when looked up). 0 = the air pathway is structurally absent, not just small." />
            </Columns>
            <NumberField label="ambient C_air [µg/m³]" value={n.C_air} min={0} step={0.1}
                         onChange={(v) => setN({ C_air: v })} help="0 = clean air, i.e. volatilisation only." />
          </>
        )}
      </Expander>
      <Expander title="⚗️ Weak electrolyte (pKa) + apoplastic bypass">
        <Checkbox label="this compound is an acid / base (has a pKa)" checked={n.weak_electrolyte}
                  onChange={(v) => setN({ weak_electrolyte: v })}
                  help="OFF = the strictly neutral path. ON = the compound is a neutral molecule AND an ion at once, weighted by Henderson-Hasselbalch; the ion feels the GHK membrane term, the neutral species does not. A lookup fills the pKa either way, but only a MEASURED one turns this on for you — a predicted pKa is left for you to switch on deliberately, because it decides which model runs." />
        {found?.pka !== null && found?.pka !== undefined && !found.pkaMeasured && !n.weak_electrolyte && (
          <Caption>{`↑ a **predicted** pKa ${fmtG(found.pka)} is filled in below and ready — tick the box if you want the weak-electrolyte path.`}</Caption>
        )}
        {n.weak_electrolyte && (
          <>
            <Columns n={2}>
              <NumberField label="pKa" value={n.pKa} min={-5} max={14} step={0.1} onChange={(v) => setN({ pKa: v })} />
              <NumberField label="root-zone pH" value={n.pH} min={3} max={10} step={0.1} onChange={(v) => setN({ pH: v })} />
            </Columns>
            <Radio label="acid or base" horizontal options={["acid", "base"] as const}
                   value={n.is_acid ? "acid" : "base"} onChange={(v) => setN({ is_acid: v === "acid" })} />
            <Caption>{`→ neutral fraction **f_n = ${fmtG(fn, 3)}**, ionic **f_d = ${fmtG(fd, 3)}**. ` +
              (n.is_acid ? "The ion is an ANION: excluded by the inside-negative membrane (the PFAS case at f_n→0)."
                : "The ion is a CATION: ATTRACTED by the inside-negative membrane, not excluded — at the same pKa the base reaches the neutral form's root BAF, the acid stays far below it (tests/test_weak_electrolyte.py).")}</Caption>
            {fn < meta.f_n_tested_floor && (
              <Alert kind="error">{`**f_n = ${fmtG(fn, 3)} is BELOW the ≈${fmtG(meta.f_n_tested_floor)} floor where this path was tested** (§4l: direction supported, magnitude REFUTED — under it the model predicts almost nothing where a measured transfer remains, because its only entry is transmembrane while a real ion also arrives apoplastically). Read this run as a lower bound on uptake, not a prediction.`}</Alert>
            )}
          </>
        )}
        <NumberField label="apoplastic bypass g_apo  [L/kg/d]" value={n.g_apo} min={0} max={50} step={0.5}
                     onChange={(v) => setN({ g_apo: v })}
                     help="A route AROUND the membrane, so it feels neither speciation nor GHK. 0 = structurally absent (the default — NOTHING is adopted)." />
        <Caption>Speciation is TESTED and BOUNDED on Schriever's ionisable rows: the direction is supported (speciation raises the rank correlation), the magnitude is not (a negative bias; ~nothing predicted at very low f_n). A bypass g_apo is not adopted: no value of it is supported by a source. validation/weak_electrolyte_tscf.py; §4l–4m of docs/neutral_dpu_validation.md.</Caption>
      </Expander>
      <NumberField label="soil Koc  [L/kg]" value={kocShown} min={0} max={meta.koc_max}
                   format={(v) => fmtG(v, 4)} onChange={(v) => setN({ Koc: v })}
                   help={`Used by the flooded Cwᵒ(t) shape and the live HYDRUS run (Kd = Koc·f_oc). Default ${fmtF(found?.koc ?? kar, 1)} from ${found?.koc != null ? "CompTox" : "Karickhoff"} (Karickhoff would give ${fmtF(kar, 1)}) — PROVISIONAL: no table in this repo scores a predicted Koc, and Li 2019's soil half shows the sorption term moves the result (its bias differs between rows with an estimated and a measured K_om; validation/li2019_soil_table.py). Type a MEASURED Kd/f_oc here whenever you have one.`} />
    </>
  );
}

// ---------------------------------------------------------------- PFAS mechanism
function MechanismPanel({ st, set }: { st: ExpertState; set: (s: ExpertState) => void }) {
  const m = st.scn.mechanism;
  const setM = (patch: Partial<Scenario["mechanism"]>) =>
    set({ ...st, scn: { ...st.scn, mechanism: { ...m, ...patch } } });
  return (
    <Expander title="⚙️ Mechanism (advanced)">
      <Radio label="Shoot model" options={SHOOT} value={m.shoot === "redistributed" ? SHOOT[0] : SHOOT[1]}
             onChange={(v) => setM({ shoot: v.startsWith("redistributed") ? "redistributed" : "4pool" })}
             help="Redistributed: N stem segments + leaf, each organ keeping the solute its own transpiration deposits; the sheath+culm and panicle shares of transpiration are sourced ranges (Sugimoto 1973/1975 and others), run at their midpoints. 4-pool: one straw pool. Each was fitted on its own." />
      <Checkbox label="Lipid-facilitated loading" checked={m.lipid_loading}
                onChange={(v) => setM({ lipid_loading: v })}
                help="Adds the B-independent bound-loading term g·C: the free anion Cw = C/B starves high-binding long chains, but the membrane-lipid-bound pool can still ride the transpiration stream. Its constants are fitted on Yamazaki only (terms the data do not identify are removed); on other datasets it is a transfer. It supplies its own f_xy, so the 'Root→shoot f_xy' choice applies only with it off." />
      <Caption>Defaults = the default model. Changing a switch moves every tab (map, dynamics, BAF), never `parameters.json`.</Caption>
    </Expander>
  );
}

function SmilesPanel({ meta, st, set }: { meta: Meta; st: ExpertState; set: (s: ExpertState) => void }) {
  const smiles = st.scn.compound.smiles ?? "";
  const { data: mol } = useAsync(smiles && meta.available.rdkit ? () => getMolSvg(smiles) : null, smiles);
  return (
    <>
      <Select label="Example structure" options={meta.example_smiles.map(([l]) => l)} value={st.smilesExample}
              help="Known → reads the curated parameters; NOVEL → QSPR (provisional)."
              onChange={(l) => set({ ...st, smilesExample: l, scn: { ...st.scn, compound: { ...st.scn.compound,
                smiles: meta.example_smiles.find(([x]) => x === l)![1] } } })} />
      <TextField label="SMILES" multiline commitOnBlur value={smiles}
                 onChange={(v) => set({ ...st, scn: { ...st.scn, compound: { ...st.scn.compound, smiles: v.trim() } } })} />
      {!meta.available.rdkit ? (
        <Alert kind="warning">{"RDKit is not installed — the SMILES mode needs it.\n\n`pip install rdkit` (or `-r requirements-structure.txt`). Meanwhile use **Curated congener**."}</Alert>
      ) : mol?.svg ? (
        <><div className="mol" dangerouslySetInnerHTML={{ __html: mol.svg }} /><Caption>structure (RDKit)</Caption></>
      ) : mol?.error ? <Caption>{`⚠ ${mol.error}`}</Caption> : null}
    </>
  );
}

// ---------------------------------------------------------------- exposure (section 3)
function ScenarioSection({ meta, st, set, onHydrusBuilt }: {
  meta: Meta; st: ExpertState; set: (s: ExpertState) => void; onHydrusBuilt: () => void;
}) {
  const e = st.scn.exposure;
  const c = st.scn.compound;
  const crop = st.scn.mechanism.crop;
  const setE = (patch: Partial<Scenario["exposure"]>, extra: Partial<ExpertState> = {}) =>
    set({ ...st, ...extra, scn: { ...st.scn, exposure: { ...e, ...patch } } });
  const neutral = c.kind === "neutral";
  const [building, setBuilding] = useState<{ busy: boolean; log: string[] | null; ok?: boolean }>({ busy: false, log: null });
  const previewBody = { compound: c, crop, level: e.Cwo, season: e.season };
  const preview = useAsync(e.mode === "parametric" && e.cwo_profile === "flooded"
                           ? () => getCwoPreview(previewBody) : null, JSON.stringify(previewBody) + e.mode + e.cwo_profile);
  const koc = neutral ? (c.neutral.Koc ?? karickhoff(meta, c.neutral.log_kow)) : undefined;
  const invKey = JSON.stringify([c.kind, c.congener, koc, e.mode]);
  const inv = useAsync(e.mode === "inventory"
    ? () => getInventoryDefaults(neutral ? { log_kow: c.neutral.log_kow, koc }
                                 : c.kind === "congener" ? { congener: c.congener } : {})
    : null, invKey);
  switch (e.mode) {
    case "parametric":
      return (
        <>
          <NumberField label="Pore-water Cwᵒ  [µg/L]" value={e.Cwo} min={0} step={0.1} onChange={(v) => setE({ Cwo: v })}
                       help="Free anion conc. driving root uptake. 1.0 → tissue conc equals BAF." />
          {crop !== "site2020" && (
            <Slider label="Season length  [days]" min={90} max={160} step={5} value={e.season}
                    onChange={(v) => setE({ season: v })} />
          )}
          <Radio label="Pore-water Cwᵒ(t) shape" options={SHAPES}
                 value={e.cwo_profile === "constant" ? SHAPES[0] : SHAPES[1]}
                 onChange={(v) => setE({ cwo_profile: v.startsWith("constant") ? "constant" : "flooded" })}
                 help="constant → Cwᵒ held flat (tissue conc == BAF). flooded → the Tsukuba Andosol paddy's root-zone pore water from stored HYDRUS-1D runs (contaminated soil, clean water): short chains LEACH as a front passes the root zone, long chains stay buffered; the season-MEAN is held at Cwᵒ, so only the time shape changes. No engine needed (for a live run use the 'Run HYDRUS-1D (live)' data source)." />
          {e.cwo_profile === "flooded" && preview.data && <Plot figure={preview.data.figure} height={230} />}
        </>
      );
    case "tables":
      return (
        <>
          <Caption>Enter the **growth** and **pore-water Cwᵒ(t)** tables in the main panel → (editable grids). Q_TP defaults to the measured transpiration; omit either table to fall back to the built-in value.</Caption>
          <NumberField label="Fallback constant Cwᵒ [µg/L]" value={e.Cwo} min={0} step={0.1} onChange={(v) => setE({ Cwo: v })}
                       help="Used where the Cwᵒ table is missing, and to seed it." />
        </>
      );
    case "csv":
      return (
        <>
          <Caption>CSV columns: `t, Cwo, Qtp, M_root, M_stem, M_leaf, M_grain` (Qtp/M optional → measured forcings). See the **About** tab for the HYDRUS map.</Caption>
          <FileText label="Driver CSV (HYDRUS-1D / Phydrus output)"
                    onText={(txt, name) => setE({ csv_text: txt }, { csvName: name })} />
          {st.csvName && <Caption>{`📎 ${st.csvName}`}</Caption>}
          <Checkbox label="Use bundled example" checked={e.csv_text === null}
                    onChange={(v) => setE({ csv_text: v ? null : "" }, { csvName: null })} />
        </>
      );
    case "hydrus":
      if (!meta.available.hydrus) {
        return (
          <>
            <Alert kind="warning">HYDRUS-1D engine not built in this environment yet — build it once (compiles the FORTRAN solver with gfortran; ~1 min, cached).</Alert>
            <Button disabled={building.busy} onClick={async () => {
              setBuilding({ busy: true, log: null });
              try { const r = await buildHydrus(); setBuilding({ busy: false, log: r.log, ok: r.ok }); if (r.ok) onHydrusBuilt(); }
              catch (ex) { setBuilding({ busy: false, log: [String(ex)], ok: false }); }
            }}>⚙ Build the HYDRUS-1D engine now</Button>
            {building.busy && <Spinner label="Compiling HYDRUS-1D with gfortran…" />}
            {building.log && (building.ok ? <Alert kind="success">✓ Engine built — loading the live mode…</Alert>
              : <><Alert kind="error">Build failed — details below.</Alert><pre className="code">{building.log.join("\n")}</pre></>)}
            <Caption>The FORTRAN source is VENDORED in the repo (no submodule to init): copy `external/hydrus_source/makefile` into `source/`, run `make` there (gfortran), `pip install phydrus`. Until built, the tool falls back to the parametric model.</Caption>
          </>
        );
      }
      return (
        <>
          <Caption>Runs a **real HYDRUS-1D** paddy model (Richards + advection–dispersion + linear Kd + root uptake) → Cwᵒ(t), Q_TP(t). Cached per setting.</Caption>
          <Slider label="Soil organic carbon f_oc" min={0.005} max={0.05} step={0.0025} value={e.f_oc ?? meta.hydrus_foc}
                  fmt={(v) => fmtF(v, 4)} onChange={(v) => setE({ f_oc: v })}
                  help={(neutral ? "Kd = Koc·f_oc → retardation R = 1+ρKd/θ. Koc is the value in the compound panel above."
                    : "Kd = Koc(chain length)·f_oc → per-congener retardation R = 1+ρKd/θ.") +
                    ` Default ${fmtG(meta.hydrus_foc)}: the site's (Yamazaki SI Table S5).`} />
          <Caption>The water calendar (flooding, midsummer drainage, intermittent irrigation) and the percolation are the site's sourced ones (soil_hydrus, PREREG 9.2–9.14); they are not free settings.</Caption>
        </>
      );
    case "inventory": {
      const seed = inv.data;
      return (
        <>
          <Caption>A total soil load at transplanting is split between sorbed and dissolved at the flooded plough layer's water (Tsukuba Andosol), and the pore water then follows the paddy's own HYDRUS-1D shape for this sorption.</Caption>
          <NumberField label="Total soil inventory  [µg/kg dry]" value={e.C_total ?? seed?.C_total ?? 0} min={0} max={1e4} step={0.01}
                       format={(v) => fmtG(v, 4)} onChange={(v) => setE({ C_total: v })}
                       help={seed?.C_total_sourced ? "Default: the Yamazaki field's measured pre-season soil (SI Table S16, water + methanol extracts)."
                         : "Enter the measured inventory (no default for this compound)."} />
          <NumberField label="Isotherm K_F  [L/kg]" value={e.K_F ?? seed?.K_F ?? 1} min={1e-4} max={1e5}
                       format={(v) => fmtG(v, 4)} onChange={(v) => setE({ K_F: v })}
                       help="Default: the compound's linear Kd on the site soil (Koc × f_oc 3.75 %)." />
          <Slider label="Isotherm exponent n" min={0.6} max={1} step={0.01} value={e.n_F} fmt={(v) => fmtF(v, 2)}
                  onChange={(v) => setE({ n_F: v })} help="1 = linear, the form the soil scenario runs on." />
        </>
      );
    }
    case "biomon": {
      const b = e.biomon ?? { conc: { root: 0.49, straw: 0.83, grain: 0.46 }, Cwo: 1.0 };
      const setB = (patch: Partial<typeof b>) => setE({ biomon: { ...b, ...patch }, biomon_csv_text: null });
      return (
        <>
          <Caption>Enter MEASURED tissue concentrations + the pore-water/soil-solution Cwᵒ. BAF is read straight off the data — no HYDRUS run needed.</Caption>
          <Radio label="Input" horizontal options={["Manual", "Upload CSV (tissue,conc[,Cwo])"] as const} value={st.biomonInput}
                 onChange={(v) => set({ ...st, biomonInput: v, scn: { ...st.scn, exposure: { ...e,
                   biomon: v === "Manual" ? b : null, biomon_csv_text: null } }, biomonCsvName: null })} />
          {st.biomonInput === "Manual" ? (
            <>
              <NumberField label="Pore-water Cwᵒ  [µg/L]" value={b.Cwo ?? 1} min={1e-6} max={1e4} step={0.1}
                           format={(v) => fmtF(v, 4)} onChange={(v) => setB({ Cwo: v })} />
              {(["root", "straw", "grain"] as const).map((k) => (
                <NumberField key={k} label={`${k} conc  [µg/kg]`} value={b.conc[k] ?? 0} min={0} max={1e6} step={0.1}
                             onChange={(v) => setB({ conc: { ...b.conc, [k]: v } })} />
              ))}
            </>
          ) : (
            <>
              <FileText label="Biomonitoring CSV" onText={(txt, name) => setE({ biomon_csv_text: txt, biomon: null }, { biomonCsvName: name })} />
              {st.biomonCsvName ? <Caption>{`📎 ${st.biomonCsvName}`}</Caption>
                : <Caption>Using examples/biomonitoring_example.csv (Yamazaki PFOA) until a file is chosen.</Caption>}
            </>
          )}
        </>
      );
    }
  }
  return null;
}

export function ExpertSidebar({ meta, st, set, run, onHydrusBuilt }: {
  meta: Meta; st: ExpertState; set: (s: ExpertState) => void; run: RunResult | null; onHydrusBuilt: () => void;
}) {
  const scn = st.scn;
  const c = scn.compound;
  const neutral = c.kind === "neutral";
  const setScn = (patch: Partial<Scenario>) => set({ ...st, scn: { ...scn, ...patch } });
  useEffect(() => {        // a new compound reseeds the soil-inventory defaults (Streamlit widget identity)
    if (scn.exposure.C_total !== null || scn.exposure.K_F !== null)
      set({ ...st, scn: { ...scn, exposure: { ...scn.exposure, C_total: null, K_F: null } } });
  }, [c.kind, c.congener, c.smiles, c.neutral.log_kow, c.neutral.Koc]);   // eslint-disable-line react-hooks/exhaustive-deps
  const cropLabel = (k: Crop) => meta.crops.find(([x]) => x === k)?.[1] ?? k;
  return (
    <>
      <h3>1 · Data source</h3>
      <Radio label="How is the soil exposure supplied?" options={MODES.map(([, l]) => l)} value={modeLabel(scn.exposure.mode)}
             onChange={(l) => setScn({ exposure: { ...scn.exposure, mode: MODES.find(([, x]) => x === l)![0] } })}
             help="Ways to feed the plant model. 'Custom tables' lets you type/paste your own growth curve + time-varying Cwᵒ; 'Run HYDRUS-1D (live)' executes the real engine (if built); biomonitoring needs no soil model." />
      <h3>2 · Compound</h3>
      <Radio label="Compound class" options={CLASSES} value={classOf(c.kind)}
             onChange={(v) => {
               const kind = v === CLASSES[1] ? "smiles" : v === CLASSES[2] ? "neutral" : "congener";
               const smiles = kind === "smiles" && !c.smiles
                 ? meta.example_smiles.find(([l]) => l === st.smilesExample)![1] : c.smiles;
               setScn({ compound: { ...c, kind, smiles } });
             }}
             help="WHICH COMPOUND CLASS the run is about, not just which molecule — the first two are PFAS (permanently dissociated anion; root entry fitted per congener), the third is the Briggs/Kow neutral base (z=0, nothing fitted). Everything downstream — map, dynamics, drivers, exposure modes, inverse, downloads — follows this choice." />
      {c.kind === "congener" && (
        <Select label="PFAS congener" options={meta.congeners} value={c.congener}
                format={(x) => congenerLabel(meta, x, "en")} onChange={(v) => setScn({ compound: { ...c, congener: v } })} />
      )}
      {c.kind === "smiles" && <SmilesPanel meta={meta} st={st} set={set} />}
      {neutral && <NeutralPanel meta={meta} st={st} set={set} />}
      {!neutral && (
        <Radio label="Root→shoot loading f_xy" options={FXY}
               value={scn.mechanism.fxy_source === "fit" ? FXY[0] : FXY[1]}
               onChange={(v) => setScn({ mechanism: { ...scn.mechanism, fxy_source: v.startsWith("recommended") ? "recommended" : "fit" } })}
               help="Applies only with lipid loading OFF (⚙️ Mechanism below; off by default); lipid loading supplies its own f_xy. fit (DEFAULT): f_xy fitted per congener to Yamazaki's straw with L_Ph and r_eq to its grain and root (a saturated reproduction; chosen by the transfer tests, PREREG 5.3). recommended: the monotone physical f_xy, with L_Ph and r_eq fitted given it — straw is a prediction. validation/refit_sourced.py." />
      )}
      <Radio label="Crop drivers (transpiration + organ mass)" options={meta.crops.map(([, l]) => l)}
             value={cropLabel(scn.mechanism.crop)}
             onChange={(l) => setScn({ mechanism: { ...scn.mechanism, crop: meta.crops.find(([, x]) => x === l)![0] } })}
             help="DEFAULT: the calibration site-year — JMA Tsukuba 2020 weather, the ORYZA2000 port (IR72.D92) transplanted on Yamazaki's date with phenology on its heading/harvest dates, ORYZA2000 transpiration (FAO-56 terms); the season is transplanting to harvest. The transport parameters are fitted on it. 'generic climatology' = forcing_rice transpiration + a biomass driver over the season slider (alternatives for sensitivity; the reduced ORYZA matures 12 d late on its own calibration experiment, growth_rice has no leaf death). Ignored when a driver CSV supplies M." />
      {!neutral && <MechanismPanel st={st} set={set} />}
      <h3>3 · Scenario</h3>
      <ScenarioSection meta={meta} st={st} set={set} onHydrusBuilt={onHydrusBuilt} />
      {run && run.messages.length > 0 && run.mode === scn.exposure.mode &&
        run.messages.map((m, i) => <Alert key={i} kind={m.level}>{m.text}</Alert>)}
      <hr />
      {!neutral && (
        <MultiSelect label="Compare congeners (overlay)" options={meta.congeners} value={st.compare}
                     onChange={(v) => set({ ...st, compare: v })} help="Shown in the 'Compare' tab." />
      )}
    </>
  );
}
