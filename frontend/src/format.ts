// Python-compatible number formatting, so every number the web UI prints reads
// exactly like the Streamlit app / the Python-built strings it is compared with
// ("{:.3g}" -> 0.479, 1.23e+03; "{:.0f}"). Python rounds the EXACT binary value
// half-to-even; JS's toFixed/toExponential round exact ties up (1234.5 -> 1235
// where Python says 1234), so the rounding is done here on an exact expansion.

/** Round a digit string to `keep` digits, half-to-even on exact ties.
 *  Returns [digits, carried] where carried=true means it overflowed (999 -> 1000). */
function roundDigits(head: string, tail: string): [string, boolean] {
  const first = tail.charCodeAt(0) - 48;
  const restNonZero = /[1-9]/.test(tail.slice(1));
  let up = false;
  if (first > 5 || (first === 5 && restNonZero)) up = true;
  else if (first === 5) up = (head.charCodeAt(head.length - 1) - 48) % 2 === 1;   // exact tie
  if (!up) return [head, false];
  const d = head.split("").map((c) => c.charCodeAt(0) - 48);
  let i = d.length - 1;
  while (i >= 0) {
    if (d[i] === 9) { d[i] = 0; i--; } else { d[i]++; break; }
  }
  if (i < 0) return ["1" + d.join(""), true];
  return [d.join(""), false];
}

/** Significant digits of |x| (x > 0) rounded to p digits: [digits, decimal exponent]. */
function sigDigits(x: number, p: number): [string, number] {
  const [m, e] = x.toExponential(Math.min(100, p + 24)).split("e");
  const all = m.replace(".", "");
  let exp = parseInt(e, 10);
  let [head, carried] = roundDigits(all.slice(0, p), all.slice(p));
  if (carried) { head = head.slice(0, p); exp += 1; }
  return [head, exp];
}

const stripZeros = (s: string) => (s.includes(".") ? s.replace(/\.?0+$/, "") : s);

/** Python `format(v, ".{p}g")`. */
export function fmtG(v: number | null | undefined, p = 6): string {
  if (v === null || v === undefined) return "None";
  if (Number.isNaN(v)) return "nan";
  if (!Number.isFinite(v)) return v > 0 ? "inf" : "-inf";
  if (v === 0) return Object.is(v, -0) ? "-0" : "0";
  const prec = p === 0 ? 1 : p;
  const sign = v < 0 ? "-" : "";
  const [d, e] = sigDigits(Math.abs(v), prec);
  if (e < -4 || e >= prec) {
    const mant = stripZeros(d[0] + (d.length > 1 ? "." + d.slice(1) : ""));
    const ae = Math.abs(e);
    return `${sign}${mant}e${e < 0 ? "-" : "+"}${ae < 10 ? "0" + ae : ae}`;
  }
  let s: string;
  if (e >= 0) s = d.slice(0, e + 1) + (d.length > e + 1 ? "." + d.slice(e + 1) : "");
  else s = "0." + "0".repeat(-e - 1) + d;
  return sign + stripZeros(s);
}

/** Python `format(v, ".{d}f")`. */
export function fmtF(v: number | null | undefined, d = 0): string {
  if (v === null || v === undefined) return "None";
  if (Number.isNaN(v)) return "nan";
  if (!Number.isFinite(v)) return v > 0 ? "inf" : "-inf";
  const a = Math.abs(v);
  if (a >= 1e21) return v.toFixed(d);
  const [ip, fp = ""] = a.toFixed(Math.min(100, d + 25)).split(".");
  const [r, carried] = roundDigits(ip + fp.slice(0, d), fp.slice(d) || "0");
  const intLen = ip.length + (carried ? 1 : 0);
  const s = d > 0 ? r.slice(0, intLen) + "." + r.slice(intLen) : r;
  return (v < 0 || Object.is(v, -0) ? "-" : "") + s;
}

/** A tiny `str.format` for the shared i18n templates: {name}, {name:.3g}, {name:.0f}. */
export function pyFormat(template: string, args: Record<string, unknown>): string {
  return template.replace(/\{(\w+)(?::([^}]*))?\}/g, (whole, name: string, spec?: string) => {
    if (!(name in args)) return whole;
    const v = args[name];
    if (spec && typeof v === "number") {
      const g = /^\.(\d+)g$/.exec(spec);
      if (g) return fmtG(v, parseInt(g[1], 10));
      const f = /^\.(\d+)f$/.exec(spec);
      if (f) return fmtF(v, parseInt(f[1], 10));
    }
    return String(v);
  });
}

/** Python's round(): half-to-even on exact ties (used by the SVG colour ramp). */
export function roundHalfEven(x: number): number {
  const fl = Math.floor(x);
  const diff = x - fl;
  if (diff > 0.5) return fl + 1;
  if (diff < 0.5) return fl;
  return fl % 2 === 0 ? fl : fl + 1;
}

export const isNum = (v: unknown): v is number => typeof v === "number" && Number.isFinite(v);
