// Typed client for the FastAPI backend. Every model computation happens server-side;
// identical requests are de-duplicated here (in-flight + a small LRU), so revisiting a
// tab or re-selecting a value is instant.

import type { CtxResult, Estimate, FigurePayload, Lang, Meta, RunResult, Scenario, Table } from "./types";

export class ApiError extends Error {
  status: number;
  constructor(status: number, message: string) { super(message); this.status = status; }
}

async function request<T>(path: string, init?: RequestInit): Promise<T> {
  const r = await fetch(path, init);
  if (!r.ok) {
    let detail = r.statusText;
    try {
      const j = await r.json();
      detail = typeof j.detail === "string" ? j.detail : JSON.stringify(j.detail);
    } catch { /* not JSON */ }
    throw new ApiError(r.status, detail);
  }
  return (await r.json()) as T;
}

const post = (body: unknown, signal?: AbortSignal): RequestInit => ({
  method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body), signal,
});

// ---- small promise LRU keyed by the request body ----
class PromiseCache<T> {
  private m = new Map<string, Promise<T>>();
  constructor(private max = 120) {}
  get(key: string, make: () => Promise<T>): Promise<T> {
    const hit = this.m.get(key);
    if (hit) { this.m.delete(key); this.m.set(key, hit); return hit; }
    const p = make();
    this.m.set(key, p);
    p.catch(() => this.m.delete(key));       // never cache a failure
    while (this.m.size > this.max) this.m.delete(this.m.keys().next().value as string);
    return p;
  }
}

const runCache = new PromiseCache<RunResult>(40);
const figCache = new PromiseCache<FigurePayload>(160);

export const scenarioKey = (s: Scenario) => JSON.stringify(s);

export function getMeta(): Promise<Meta> {
  return request<Meta>("/api/meta");
}

export function runScenario(s: Scenario): Promise<RunResult> {
  const key = scenarioKey(s);
  return runCache.get(key, () => request<RunResult>("/api/run", post(s)));
}

export function getFigure(s: Scenario, name: string, lang: Lang,
                          opts: Record<string, unknown> = {}): Promise<FigurePayload> {
  const body = { scenario: s, name, lang, opts };
  return figCache.get(JSON.stringify(body), () => request<FigurePayload>("/api/figure", post(body)));
}

export function getCwoPreview(body: unknown): Promise<FigurePayload> {
  return figCache.get("cwo:" + JSON.stringify(body),
                      () => request<FigurePayload>("/api/cwo-preview", post(body)));
}

export function getTableDefaults(season: number, biomass: string, Cwo: number):
    Promise<{ growth: Table; cwo: Table; density: Record<string, number> }> {
  const q = new URLSearchParams({ season: String(season), biomass, Cwo: String(Cwo) });
  return request(`/api/tables/defaults?${q}`);
}

export function getKLeach(params: { congener?: string; log_kow?: number; koc?: number }):
    Promise<{ k_leach: number }> {
  const q = new URLSearchParams();
  for (const [k, v] of Object.entries(params)) if (v !== undefined && v !== null) q.set(k, String(v));
  return request(`/api/k-leach?${q}`);
}

export function getMolSvg(smiles: string): Promise<{ svg?: string; error?: string }> {
  return request(`/api/smiles/svg?${new URLSearchParams({ smiles })}`);
}

export function ctxLookup(q: string): Promise<CtxResult> {
  return request(`/api/ctx?${new URLSearchParams({ q })}`);
}

export function buildHydrus(): Promise<{ ok: boolean; log: string[] }> {
  return request("/api/hydrus/build", { method: "POST" });
}

/** Stream the Bayesian inverse: onProgress(done,total) per ODE solve, resolves with the result. */
export async function runInverse(
  body: { scenario: Scenario; measured: Record<string, number>; sigma_log10: number; lang: Lang },
  onProgress: (done: number, total: number) => void,
  signal?: AbortSignal,
): Promise<{ estimate: Estimate; figure: FigurePayload["figure"] }> {
  const r = await fetch("/api/inverse", post(body, signal));
  if (!r.ok || !r.body) throw new ApiError(r.status, r.statusText);
  const reader = r.body.getReader();
  const dec = new TextDecoder();
  let buf = "";
  for (;;) {
    const { value, done } = await reader.read();
    if (done) break;
    buf += dec.decode(value, { stream: true });
    let nl: number;
    while ((nl = buf.indexOf("\n")) >= 0) {
      const line = buf.slice(0, nl).trim();
      buf = buf.slice(nl + 1);
      if (!line) continue;
      const msg = JSON.parse(line);
      if (msg.type === "progress") onProgress(msg.done, msg.total);
      else if (msg.type === "error") throw new ApiError(400, msg.detail);
      else if (msg.type === "result") return { estimate: msg.estimate, figure: msg.figure };
    }
  }
  throw new ApiError(500, "the estimate stream ended without a result");
}

function filenameFrom(disposition: string | null, fallback: string): string {
  if (!disposition) return fallback;
  const star = /filename\*=UTF-8''([^;]+)/i.exec(disposition);
  if (star) return decodeURIComponent(star[1]);
  const plain = /filename="([^"]+)"/i.exec(disposition);
  return plain ? plain[1] : fallback;
}

export function saveBlob(blob: Blob, name: string) {
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = name;
  document.body.appendChild(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 2000);
}

export async function download(kind: string, s: Scenario, lang: Lang, scenarioLabel?: string | null) {
  const r = await fetch(`/api/download/${kind}`,
                        post({ scenario: s, lang, scenario_label: scenarioLabel ?? null }));
  if (!r.ok) throw new ApiError(r.status, r.statusText);
  saveBlob(await r.blob(), filenameFrom(r.headers.get("Content-Disposition"), kind));
}

/** Render an SVG string to a PNG (client-side; no kaleido/Chrome needed on the server). */
export function svgToPng(svg: string, name: string, width = 1520, height = 1120): Promise<void> {
  return new Promise((resolve, reject) => {
    const fixed = svg.replace("style=\"width:100%;height:auto;display:block\"",
                              `width="${width}" height="${height}"`);
    const img = new Image();
    img.onload = () => {
      const c = document.createElement("canvas");
      c.width = width;
      c.height = height;
      const ctx = c.getContext("2d");
      if (!ctx) { reject(new Error("no canvas")); return; }
      ctx.drawImage(img, 0, 0, width, height);
      c.toBlob((b) => { if (b) { saveBlob(b, name); resolve(); } else reject(new Error("toBlob")); },
               "image/png");
    };
    img.onerror = () => reject(new Error("could not rasterise the SVG"));
    img.src = "data:image/svg+xml;charset=utf-8," + encodeURIComponent(fixed);
  });
}
