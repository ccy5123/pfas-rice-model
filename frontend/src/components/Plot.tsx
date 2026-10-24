// Plotly figure renderer. plotly.js is loaded lazily (one ~1.4 MB gzip chunk), so the
// first paint -- headline numbers + the SVG plant map -- does not wait for it.
// Figures arrive as the JSON the Python `plots.py` builders produce (same template,
// colours, hover text as the Streamlit app); only the font colour follows the theme.

import { useEffect, useRef, useState } from "react";
import type { FigurePayload, Lang, Scenario } from "../types";
import { getFigure } from "../api";
import { useAsync, useDarkMode } from "../hooks";
import { Alert, Spinner } from "./ui";

let plotlyPromise: Promise<any> | null = null;
export const loadPlotly = () => (plotlyPromise ??= import("plotly.js-dist-min").then((m) => m.default));

export function Plot({ figure, height }: { figure: FigurePayload["figure"]; height?: number }) {
  const ref = useRef<HTMLDivElement>(null);
  const dark = useDarkMode();
  const [err, setErr] = useState<string | null>(null);
  useEffect(() => {
    let alive = true;
    const el = ref.current;
    loadPlotly().then((Plotly) => {
      if (!alive || !el) return;
      const ink = dark ? "#ECE6D9" : "#2a2a2a";
      const layout = {
        ...figure.layout,
        font: { ...(figure.layout?.font ?? {}), color: ink },
        autosize: true,
        ...(height ? { height } : {}),
      };
      Plotly.react(el, figure.data, layout, { responsive: true, displaylogo: false });
    }).catch((e) => setErr(String(e)));
    return () => { alive = false; };
  }, [figure, dark, height]);
  useEffect(() => {
    const el = ref.current;
    if (!el) return;
    const ro = new ResizeObserver(() => {
      if (plotlyPromise) plotlyPromise.then((P) => { if ((el as any)._fullLayout) P.Plots.resize(el); });
    });
    ro.observe(el);
    return () => {
      ro.disconnect();
      if (plotlyPromise) plotlyPromise.then((P) => P.purge(el));
    };
  }, []);
  if (err) return <Alert kind="error">{`Could not load the chart library: ${err}`}</Alert>;
  return <div ref={ref} className="plot" />;
}

/** Fetch a named server-built figure for the scenario and render it. */
export function ServerFigure({ scenario, name, lang, opts, height, onData }: {
  scenario: Scenario; name: string; lang: Lang; opts?: Record<string, unknown>; height?: number;
  onData?: (d: Record<string, any>) => void;
}) {
  const key = JSON.stringify([scenario, name, lang, opts ?? {}]);
  const { data, error, loading } = useAsync(() => getFigure(scenario, name, lang, opts ?? {}), key);
  useEffect(() => { if (data && onData) onData(data.data); }, [data]);   // eslint-disable-line react-hooks/exhaustive-deps
  if (error && !data) return <Alert kind="warning">{error}</Alert>;
  return (
    <div className={loading ? "fig loading" : "fig"}>
      {loading && <Spinner label={lang === "ko" ? "그래프 계산 중…" : "Computing…"} />}
      {data && <Plot figure={data.figure} height={height} />}
      {!data && !loading && null}
      {!data && loading && <div className="fig-skeleton" />}
    </div>
  );
}
