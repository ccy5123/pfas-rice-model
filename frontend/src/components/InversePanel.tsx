// "Work backwards": Bayesian estimate of the pore-water Cwᵒ from measured tissue
// concentrations (model_api.estimate_exposure_bayesian), streamed with live progress.
// Same panel for Simple (Korean) and Expert (English), like ui/common._render_inverse_estimator.

import { useMemo, useRef, useState } from "react";
import { runInverse } from "../api";
import { fmtG } from "../format";
import { PLAIN_KO, t } from "../i18n";
import type { Estimate, FigurePayload, Lang, Meta, Scenario } from "../types";
import { Plot } from "./Plot";
import { Alert, Button, Caption, Columns, Md, Metric, NumberField, Progress, Radio } from "./ui";

export function InversePanel({ meta, scenario, lang, enabled = true }: {
  meta: Meta; scenario: Scenario; lang: Lang; enabled?: boolean;
}) {
  const [vals, setVals] = useState({ root: 0, straw: 0, grain: 0 });
  const prec = meta.inverse_precision[lang];
  const [precLabel, setPrecLabel] = useState(prec[0][0]);
  const sigma = prec.find(([l]) => l === precLabel)?.[1] ?? 0.15;
  const [state, setState] = useState<{
    sig: string; prog: [number, number] | null; est: Estimate | null; fig: FigurePayload["figure"] | null;
    error: string | null; running: boolean;
  } | null>(null);
  const abort = useRef<AbortController | null>(null);
  const have = vals.root > 0 || vals.straw > 0 || vals.grain > 0;
  const sig = useMemo(() => JSON.stringify([scenario, vals, sigma]), [scenario, vals, sigma]);

  if (!enabled) return <Alert kind="info">{t(meta, "inv.not_curated", lang)}</Alert>;

  const estimate = async () => {
    abort.current?.abort();
    const ac = new AbortController();
    abort.current = ac;
    setState({ sig, prog: [0, 1], est: null, fig: null, error: null, running: true });
    const measured: Record<string, number> = {};
    for (const [k, v] of Object.entries(vals)) if (v > 0) measured[k] = v;
    try {
      const r = await runInverse({ scenario, measured, sigma_log10: sigma, lang },
                                 (d, n) => setState((s) => (s && s.sig === sig ? { ...s, prog: [d, n] } : s)),
                                 ac.signal);
      setState({ sig, prog: null, est: r.estimate, fig: r.figure, error: null, running: false });
    } catch (e) {
      if (ac.signal.aborted) return;
      setState({ sig, prog: null, est: null, fig: null, error: e instanceof Error ? e.message : String(e),
                 running: false });
    }
  };

  const current = state && state.sig === sig ? state : null;
  const est = current?.est;
  const lo = est?.ci95?.[0], hi = est?.ci95?.[1];
  const hasCi = lo !== null && lo !== undefined && hi !== null && hi !== undefined;
  const names = lang === "ko" ? PLAIN_KO : {};
  return (
    <div>
      <Md>{t(meta, "inv.intro", lang)}</Md>
      <Columns n={3}>
        <NumberField label={t(meta, "inv.in_root", lang)} value={vals.root} min={0} max={1e6} step={0.1}
                     onChange={(v) => setVals({ ...vals, root: v })} />
        <NumberField label={t(meta, "inv.in_straw", lang)} value={vals.straw} min={0} max={1e6} step={0.1}
                     onChange={(v) => setVals({ ...vals, straw: v })} />
        <NumberField label={t(meta, "inv.in_grain", lang)} value={vals.grain} min={0} max={1e6} step={0.1}
                     onChange={(v) => setVals({ ...vals, grain: v })} />
      </Columns>
      <Radio label={t(meta, "inv.precision_label", lang)} help={t(meta, "inv.precision_help", lang)}
             horizontal options={prec.map(([l]) => l)} value={precLabel} onChange={setPrecLabel} />
      <Button primary disabled={!have || !!current?.running} onClick={estimate}>
        {t(meta, "inv.estimate_btn", lang)}
      </Button>
      {!have && <Caption>{t(meta, "inv.enter_first", lang)}</Caption>}
      {have && !current && <Caption>{t(meta, "inv.press_estimate", lang)}</Caption>}
      {current?.running && current.prog && (
        <Progress frac={current.prog[0] / Math.max(1, current.prog[1])}
                  label={current.prog[0] === 0 ? t(meta, "inv.preparing", lang)
                         : t(meta, "inv.running", lang, { done: current.prog[0], total: current.prog[1] })} />
      )}
      {current?.error && <Alert kind="error">{current.error}</Alert>}
      {est && (
        <div className="inverse-result">
          <div className="cols" style={{ gridTemplateColumns: "1fr 2fr" }}>
            <Metric label={t(meta, "inv.metric_label", lang)} value={`${fmtG(est.median, 3)} µg/L`}
                    sub={hasCi ? t(meta, "inv.range95", lang, { lo, hi }) : t(meta, "inv.range_unconstrained", lang)} />
            <Md>{t(meta, "inv.summary_lead", lang, { congener: est.congener, med: est.median }) +
                 (hasCi ? t(meta, "inv.summary_ci", lang, { lo, hi }) : t(meta, "inv.summary_no_ci", lang)) +
                 t(meta, "inv.summary_tail", lang)}</Md>
          </div>
          {current?.fig && <Plot figure={current.fig} />}
          <Caption>{t(meta, "inv.fit_caption", lang, {
            rows: est.used_tissues.map((k) => t(meta, "inv.fit_row", lang, {
              name: (names as Record<string, string>)[k] ?? k, meas: est.measured[k], model: est.model_fit[k] })).join(" · "),
          })}</Caption>
        </div>
      )}
    </div>
  );
}
