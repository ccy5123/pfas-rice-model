// Minimal form + layout primitives (the Streamlit widget set, as plain React).

import { useEffect, useId, useState, type ReactNode } from "react";
import ReactMarkdown from "react-markdown";
import remarkGfm from "remark-gfm";

export function Md({ children, className }: { children: string; className?: string }) {
  return (
    <div className={"md " + (className ?? "")}>
      <ReactMarkdown remarkPlugins={[remarkGfm]}
        components={{ a: ({ node: _n, ...p }) => <a {...p} target="_blank" rel="noreferrer" /> }}>
        {children}
      </ReactMarkdown>
    </div>
  );
}

export function Caption({ children }: { children: ReactNode }) {
  return <div className="caption">{typeof children === "string" ? <Md>{children}</Md> : children}</div>;
}

export function Alert({ kind, children }: { kind: "info" | "success" | "warning" | "error"; children: ReactNode }) {
  return (
    <div className={`alert ${kind}`} role={kind === "error" ? "alert" : "status"}>
      {typeof children === "string" ? <Md>{children}</Md> : children}
    </div>
  );
}

export function Expander({ title, children, open = false }: { title: ReactNode; children: ReactNode; open?: boolean }) {
  const [isOpen, setOpen] = useState(open);
  return (
    <details className="expander" open={isOpen} onToggle={(e) => setOpen((e.target as HTMLDetailsElement).open)}>
      <summary>{title}</summary>
      {isOpen && <div className="expander-body">{children}</div>}
    </details>
  );
}

export function Tabs({ names, active, onChange }: { names: string[]; active: number; onChange: (i: number) => void }) {
  return (
    <div className="tabs" role="tablist">
      {names.map((n, i) => (
        <button key={n} role="tab" aria-selected={i === active} className={i === active ? "tab active" : "tab"}
                onClick={() => onChange(i)}>{n}</button>
      ))}
    </div>
  );
}

export function Metric({ label, value, sub }: { label: string; value: string; sub?: string }) {
  return (
    <div className="metric">
      <div className="metric-label">{label}</div>
      <div className="metric-value">{value}</div>
      {sub && <div className="metric-sub">{sub}</div>}
    </div>
  );
}

export function Field({ label, help, children }: { label: ReactNode; help?: string; children: ReactNode }) {
  return (
    <div className="field">
      <div className="field-label">{label}{help && <span className="help" title={help}>?</span>}</div>
      {children}
    </div>
  );
}

export function Radio<T extends string>({ label, options, value, onChange, help, horizontal, cards, format }: {
  label?: ReactNode; options: readonly T[]; value: T; onChange: (v: T) => void; help?: string;
  horizontal?: boolean; cards?: boolean; format?: (v: T) => ReactNode;
}) {
  const name = useId();
  const body = (
    <div className={`radio ${horizontal ? "horizontal" : ""} ${cards ? "cards" : ""}`} role="radiogroup">
      {options.map((o) => (
        <label key={o} className={o === value ? "checked" : ""}>
          <input type="radio" name={name} checked={o === value} onChange={() => onChange(o)} />
          <span>{format ? format(o) : o}</span>
        </label>
      ))}
    </div>
  );
  return label ? <Field label={label} help={help}>{body}</Field> : body;
}

export function Select<T extends string>({ label, options, value, onChange, help, format }: {
  label: ReactNode; options: readonly T[]; value: T; onChange: (v: T) => void; help?: string;
  format?: (v: T) => string;
}) {
  return (
    <Field label={label} help={help}>
      <select value={value} onChange={(e) => onChange(e.target.value as T)}>
        {options.map((o) => <option key={o} value={o}>{format ? format(o) : o}</option>)}
      </select>
    </Field>
  );
}

export function Slider({ label, min, max, step, value, onChange, help, fmt }: {
  label: ReactNode; min: number; max: number; step: number; value: number; onChange: (v: number) => void;
  help?: string; fmt?: (v: number) => string;
}) {
  return (
    <Field label={<>{label} <b className="slider-val">{fmt ? fmt(value) : value}</b></>} help={help}>
      <input type="range" min={min} max={max} step={step} value={value}
             onChange={(e) => onChange(parseFloat(e.target.value))} />
    </Field>
  );
}

/** Number input that commits on blur/Enter (no solve per keystroke), clamped to [min,max]. */
export function NumberField({ label, value, onChange, min, max, step, help, format }: {
  label: ReactNode; value: number; onChange: (v: number) => void; min?: number; max?: number;
  step?: number; help?: string; format?: (v: number) => string;
}) {
  const show = (v: number) => (format ? format(v) : String(v));
  const [text, setText] = useState(show(value));
  useEffect(() => setText(show(value)), [value]);   // eslint-disable-line react-hooks/exhaustive-deps
  const commit = () => {
    const v = parseFloat(text);
    if (!Number.isFinite(v)) { setText(show(value)); return; }
    let c = v;
    if (min !== undefined) c = Math.max(min, c);
    if (max !== undefined) c = Math.min(max, c);
    if (c !== value) onChange(c);
    setText(show(c));
  };
  return (
    <Field label={label} help={help}>
      <input type="number" value={text} step={step ?? "any"} min={min} max={max}
             onChange={(e) => setText(e.target.value)} onBlur={commit}
             onKeyDown={(e) => { if (e.key === "Enter") commit(); }} />
    </Field>
  );
}

export function TextField({ label, value, onChange, help, placeholder, multiline, commitOnBlur }: {
  label: ReactNode; value: string; onChange: (v: string) => void; help?: string; placeholder?: string;
  multiline?: boolean; commitOnBlur?: boolean;
}) {
  const [text, setText] = useState(value);
  useEffect(() => setText(value), [value]);
  const props = {
    value: text, placeholder,
    onChange: (e: { target: { value: string } }) => { setText(e.target.value); if (!commitOnBlur) onChange(e.target.value); },
    onBlur: () => { if (commitOnBlur && text !== value) onChange(text); },
  };
  return (
    <Field label={label} help={help}>
      {multiline ? <textarea rows={3} {...props} /> :
        <input type="text" {...props}
               onKeyDown={(e) => { if (e.key === "Enter" && commitOnBlur) onChange(text); }} />}
    </Field>
  );
}

export function Checkbox({ label, checked, onChange, help }: {
  label: ReactNode; checked: boolean; onChange: (v: boolean) => void; help?: string;
}) {
  return (
    <label className="checkbox">
      <input type="checkbox" checked={checked} onChange={(e) => onChange(e.target.checked)} />
      <span>{label}</span>
      {help && <span className="help" title={help}>?</span>}
    </label>
  );
}

export function Toggle({ label, checked, onChange, help }: {
  label: ReactNode; checked: boolean; onChange: (v: boolean) => void; help?: string;
}) {
  return (
    <label className="toggle" title={help}>
      <input type="checkbox" role="switch" checked={checked} onChange={(e) => onChange(e.target.checked)} />
      <span className="track"><span className="thumb" /></span>
      <span>{label}</span>
    </label>
  );
}

export function Button({ children, onClick, primary, disabled }: {
  children: ReactNode; onClick: () => void; primary?: boolean; disabled?: boolean;
}) {
  return <button className={primary ? "btn primary" : "btn"} onClick={onClick} disabled={disabled}>{children}</button>;
}

export function MultiSelect({ label, options, value, onChange, help }: {
  label: ReactNode; options: string[]; value: string[]; onChange: (v: string[]) => void; help?: string;
}) {
  return (
    <Field label={label} help={help}>
      <div className="chips">
        {options.map((o) => {
          const on = value.includes(o);
          return (
            <button key={o} className={on ? "chip on" : "chip"} aria-pressed={on}
                    onClick={() => onChange(on ? value.filter((x) => x !== o) : [...value, o])}>{o}</button>
          );
        })}
      </div>
    </Field>
  );
}

/** File picker that hands back the file's TEXT (CSV uploads go to the server as JSON). */
export function FileText({ label, onText, accept = ".csv" }: {
  label: ReactNode; onText: (text: string, name: string) => void; accept?: string;
}) {
  return (
    <Field label={label}>
      <input type="file" accept={accept} onChange={async (e) => {
        const f = e.target.files?.[0];
        if (f) onText(await f.text(), f.name);
      }} />
    </Field>
  );
}

export function Spinner({ label }: { label?: string }) {
  return <div className="spinner"><span className="dot" />{label}</div>;
}

export function Progress({ frac, label }: { frac: number; label: string }) {
  return (
    <div className="progress">
      <div className="progress-bar" style={{ width: `${Math.round(Math.min(1, frac) * 100)}%` }} />
      <span>{label}</span>
    </div>
  );
}

export function Columns({ children, n }: { children: ReactNode; n?: number }) {
  return <div className="cols" style={n ? { gridTemplateColumns: `repeat(${n}, minmax(0, 1fr))` } : undefined}>{children}</div>;
}
