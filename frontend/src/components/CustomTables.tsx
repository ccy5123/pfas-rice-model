// Editable growth + pore-water tables (+ per-compartment density) -- the
// st.data_editor panel of ui/common._render_custom_tables. Until a table is edited it
// is sent as `null`, which the server seeds with exactly the default the grid shows.

import { useEffect, useState } from "react";
import { getTableDefaults } from "../api";
import { fmtF } from "../format";
import { t } from "../i18n";
import type { Lang, Meta, RunResult, Table, TablesSpec, Tissue } from "../types";
import { TISSUES } from "../types";
import { Alert, Button, Caption, Columns, FileText, Md, NumberField, Select } from "./ui";

function parseCsv(text: string): Table {
  const lines = text.split(/\r?\n/).filter((l) => l.trim());
  if (!lines.length) return {};
  const sep = lines[0].includes("\t") ? "\t" : ",";
  const head = lines[0].split(sep).map((h) => h.trim());
  const out: Table = Object.fromEntries(head.map((h) => [h, []]));
  for (const ln of lines.slice(1)) {
    const cells = ln.split(sep);
    head.forEach((h, i) => {
      const v = parseFloat(cells[i]);
      out[h].push(Number.isFinite(v) ? v : null);
    });
  }
  return out;
}

function Grid({ table, onChange }: { table: Table; onChange: (t: Table) => void }) {
  const cols = Object.keys(table);
  const n = Math.max(0, ...cols.map((c) => table[c].length));
  const set = (c: string, i: number, raw: string) => {
    const v = raw.trim() === "" ? null : parseFloat(raw);
    const next: Table = { ...table, [c]: [...table[c]] };
    next[c][i] = v !== null && Number.isFinite(v) ? v : null;
    onChange(next);
  };
  const paste = (c0: number, r0: number, text: string): boolean => {
    const rows = text.split(/\r?\n/).filter((l) => l.length);
    if (rows.length < 2 && !rows[0]?.includes("\t")) return false;
    const next: Table = Object.fromEntries(cols.map((c) => [c, [...table[c]]]));
    rows.forEach((row, dr) => row.split(/\t|,/).forEach((cell, dc) => {
      const c = cols[c0 + dc];
      if (!c) return;
      const v = parseFloat(cell);
      while (next[c].length <= r0 + dr) next[c].push(null);
      next[c][r0 + dr] = Number.isFinite(v) ? v : null;
    }));
    const len = Math.max(...cols.map((c) => next[c].length));
    for (const c of cols) while (next[c].length < len) next[c].push(null);
    onChange(next);
    return true;
  };
  return (
    <div className="grid-wrap">
      <table className="grid">
        <thead><tr>{cols.map((c) => <th key={c}>{c}</th>)}<th /></tr></thead>
        <tbody>
          {Array.from({ length: n }, (_, i) => (
            <tr key={i}>
              {cols.map((c, j) => (
                <td key={c}>
                  <input value={table[c][i] ?? ""} inputMode="decimal"
                         onChange={(e) => set(c, i, e.target.value)}
                         onPaste={(e) => { if (paste(j, i, e.clipboardData.getData("text"))) e.preventDefault(); }} />
                </td>
              ))}
              <td><button className="row-del" title="remove row"
                          onClick={() => onChange(Object.fromEntries(cols.map((c) => [c, table[c].filter((_, k) => k !== i)])))}>×</button></td>
            </tr>
          ))}
        </tbody>
      </table>
      <button className="row-add" onClick={() => onChange(Object.fromEntries(cols.map((c) => [c, [...table[c], null]])))}>+ row</button>
    </div>
  );
}

export function CustomTables({ meta, lang, tables, onChange, season, biomass, Cwo, run }: {
  meta: Meta; lang: Lang; tables: TablesSpec; onChange: (t: TablesSpec) => void;
  season: number; biomass: string; Cwo: number; run: RunResult | null;
}) {
  const [defaults, setDefaults] = useState<{ growth: Table; cwo: Table } | null>(null);
  useEffect(() => {
    let alive = true;
    getTableDefaults(season, biomass, Cwo).then((d) => { if (alive) setDefaults(d); });
    return () => { alive = false; };
  }, [season, biomass, Cwo]);
  const organ = meta.organ_labels[lang];
  const growth = tables.growth ?? defaults?.growth ?? null;
  const cwo = tables.cwo ?? defaults?.cwo ?? null;
  const readErr = run?.messages.find((m) => m.level === "warning" && run.mode === "tables");
  return (
    <div>
      <Md>{t(meta, "ct.intro", lang)}</Md>
      <Columns n={2}>
        <div>
          <Caption>{t(meta, "ct.growth_caption", lang)}</Caption>
          <FileText label={t(meta, "ct.growth_upload", lang)}
                    onText={(txt) => onChange({ ...tables, growth: parseCsv(txt) })} />
          {growth && <Grid table={growth} onChange={(g) => onChange({ ...tables, growth: g })} />}
          <Select label={t(meta, "ct.growth_units", lang)} value={tables.growth_units}
                  options={["g/hill", "kg/hill", "g/m2", "kg/ha", "t/ha"] as const}
                  onChange={(u) => onChange({ ...tables, growth_units: u })} />
        </div>
        <div>
          <Caption>{t(meta, "ct.cwo_caption", lang)}</Caption>
          <FileText label={t(meta, "ct.cwo_upload", lang)}
                    onText={(txt) => onChange({ ...tables, cwo: parseCsv(txt) })} />
          {cwo && <Grid table={cwo} onChange={(c) => onChange({ ...tables, cwo: c })} />}
        </div>
      </Columns>
      {(tables.growth || tables.cwo) && (
        <Button onClick={() => onChange({ ...tables, growth: null, cwo: null })}>
          {lang === "ko" ? "↺ 기본 표로 되돌리기" : "↺ Reset both tables to the defaults"}
        </Button>
      )}
      <Md>{t(meta, "ct.density_md", lang)}</Md>
      <Columns n={4}>
        {TISSUES.map((o: Tissue) => (
          <NumberField key={o} label={`ρ ${organ[o]}`} min={0.05} max={2} step={0.05}
                       value={tables.density[o]}
                       onChange={(v) => onChange({ ...tables, density: { ...tables.density, [o]: v } })} />
        ))}
      </Columns>
      {readErr && <Alert kind="warning">{readErr.text}</Alert>}
      {run && run.effective.drivers && run.mode === "tables" && (() => {
        // implied end-of-season organ volume = fresh mass / density (client-side: instant)
        const last = run.t.length - 1;
        const rows = TISSUES.map((o) =>
          `${organ[o]} ${fmtF((run.M[o][last] / Math.max(tables.density[o], 1e-6)) * 1e3, 0)} mL`).join(" · ");
        return <Caption>{t(meta, "ct.implied_volume", lang, { rows })}</Caption>;
      })()}
    </div>
  );
}
