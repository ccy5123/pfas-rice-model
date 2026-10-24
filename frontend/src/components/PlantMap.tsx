import { useEffect, useMemo, useState } from "react";
import { nearestIndex, plantSvg, plantSvgFromRun } from "../plantSvg";
import type { Lang, RunResult } from "../types";
import { Radio, Slider } from "./ui";

/** The season-scrubbable plant map. Everything here is client-side: moving the day
 *  slider or switching conc/BAF re-draws the SVG without a server round trip. */
export function PlantMap({ run, lang, grainSignal, showMetric, onSvg }: {
  run: RunResult; lang: Lang; grainSignal?: string | null; showMetric?: boolean;
  onSvg?: (svg: string) => void;
}) {
  const t0 = run.t[0], t1 = run.t[run.t.length - 1];
  const [day, setDay] = useState(t1);
  const [metric, setMetric] = useState<"conc" | "baf">("conc");
  useEffect(() => { setDay(t1); }, [t1]);
  const ti = nearestIndex(run.t, day);
  const svg = useMemo(() => plantSvgFromRun(run, ti, { lang, metric, grainSignal }),
                      [run, ti, lang, metric, grainSignal]);
  useEffect(() => { onSvg?.(svg); }, [svg]);   // eslint-disable-line react-hooks/exhaustive-deps
  return (
    <div>
      {showMetric && (
        <Radio label="Colour by" horizontal options={["concentration", "BAF"] as const}
               value={metric === "baf" ? "BAF" : "concentration"}
               onChange={(v) => setMetric(v === "BAF" ? "baf" : "conc")} />
      )}
      <Slider label={lang === "ko" ? "이앙 후 일수" : "Day after transplant"} min={t0} max={t1} step={1}
              value={day} onChange={setDay} fmt={(v) => v.toFixed(0)}
              help={lang === "ko" ? "슬라이더를 옮기면 그 날짜의 축적 상태로 색이 바뀝니다." : undefined} />
      <div className="plant-svg" dangerouslySetInnerHTML={{ __html: svg }} />
    </div>
  );
}

/** Biomonitoring: compartments coloured by the MEASURED concentrations. */
export function MeasuredPlantMap({ conc, cwo, lang }: { conc: Record<string, number>; cwo: number | null; lang: Lang }) {
  const finite = Object.values(conc).filter((v) => Number.isFinite(v));
  const cmax = finite.length ? Math.max(...finite) : 1;
  const svg = plantSvg(conc, { cmin: 0, cmax, cwo, lang });
  return <div className="plant-svg" dangerouslySetInnerHTML={{ __html: svg }} />;
}
