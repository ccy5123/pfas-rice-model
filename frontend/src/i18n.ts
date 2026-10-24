// The bilingual copy shared with the Streamlit app comes from ui/i18n.py via /api/meta,
// so both front ends say exactly the same thing. Page-specific copy lives in the views.
import { pyFormat } from "./format";
import type { Lang, Meta } from "./types";

export function t(meta: Meta, key: string, lang: Lang, args?: Record<string, unknown>): string {
  const entry = meta.strings[key];
  if (!entry) return key;
  const s = entry[lang] || entry.en;
  return args ? pyFormat(s, args) : s;
}

export const congenerLabel = (meta: Meta, name: string, lang: Lang) => meta.labels[lang]?.[name] ?? name;

export const PLAIN_KO: Record<string, string> = { root: "뿌리", stem: "줄기", leaf: "잎", grain: "낟알", straw: "짚" };

// safety signal: colour + SHAPE + label (colour-blind-safe), as ui/common.signal_badge_html
const SIGNAL: Record<string, [string, string, string, string]> = {
  green: ["safe", "●", "안전", "Safe"],
  amber: ["warn", "▲", "주의", "Caution"],
  red: ["dang", "⬢", "초과", "Exceeds"],
};
export function signalBadge(signal: string | null | undefined, ko = true, extra = "") {
  if (!signal || !(signal in SIGNAL)) return null;
  const [cls, shape, k, e] = SIGNAL[signal];
  return { cls, text: `${shape} ${ko ? `${k} · ${e}` : e}${extra ? " " + extra : ""}` };
}
