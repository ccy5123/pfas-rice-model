// App shell: sidebar (Simple ↔ Expert) + main view. All interaction state lives here
// in the browser; the server is asked for a run only when the debounced scenario
// actually changes, and each tab fetches its own figures only when it is opened.

import { useEffect, useMemo, useState } from "react";
import { getMeta, runScenario, scenarioKey } from "./api";
import { ErrorBoundary } from "./components/ErrorBoundary";
import { Alert, Caption, Spinner, Toggle } from "./components/ui";
import { useAsync, useDebounced } from "./hooks";
import { t } from "./i18n";
import type { Meta, Scenario } from "./types";
import { defaultScenario } from "./types";
import { ExpertHeader, ExpertView } from "./views/Expert";
import type { ExpertState } from "./views/ExpertSidebar";
import { ExpertSidebar } from "./views/ExpertSidebar";
import { initialSimple, SimpleHeader, SimpleSidebar, SimpleView, simpleScenario } from "./views/Simple";
import type { SimpleState } from "./views/Simple";

function initialExpert(meta: Meta): ExpertState {
  return { scn: defaultScenario("en"), compare: ["PFBA", "PFOA", "PFDA", "PFOS"].filter((c) => meta.congeners.includes(c)),
           smilesExample: meta.example_smiles[0][0], biomonInput: "Manual", csvName: null, biomonCsvName: null,
           advMech: false, lookup: null };
}

function Footer({ meta, expert }: { meta: Meta; expert: boolean }) {
  const lang = expert ? "en" : "ko";
  return (
    <footer>
      <hr />
      <Caption>{t(meta, "disclaimer", lang)}</Caption>
      <div className="cols" style={{ gridTemplateColumns: "repeat(3, minmax(0,1fr))" }}>
        <Caption>{`**PFAS–Rice Uptake Model** · v${meta.app_version}`}</Caption>
        <Caption>{t(meta, "footer.links", lang, { repo: meta.repo_url, docs: meta.docs_url })}</Caption>
        <Caption>{t(meta, "footer.cite", lang)}</Caption>
      </div>
    </footer>
  );
}

export default function App() {
  const [metaVersion, setMetaVersion] = useState(0);
  const metaState = useAsync(getMeta, String(metaVersion));
  const meta = metaState.data;
  if (!meta) {
    return (
      <div className="boot">
        {metaState.error ? <Alert kind="error">{`Cannot reach the model server: ${metaState.error}`}</Alert>
                         : <Spinner label="🌾 Loading…" />}
      </div>
    );
  }
  return <Loaded meta={meta} onMetaRefresh={() => setMetaVersion((v) => v + 1)} />;
}

function Loaded({ meta, onMetaRefresh }: { meta: Meta; onMetaRefresh: () => void }) {
  const [expert, setExpert] = useState(() => new URLSearchParams(location.search).get("mode") === "expert");
  const [simple, setSimple] = useState<SimpleState>(() => initialSimple(meta));
  const [ex, setEx] = useState<ExpertState>(() => initialExpert(meta));
  const [navOpen, setNavOpen] = useState(false);

  useEffect(() => {
    document.documentElement.lang = expert ? "en" : "ko";
    const u = new URL(location.href);
    if (expert) u.searchParams.set("mode", "expert"); else u.searchParams.delete("mode");
    history.replaceState(null, "", u);
  }, [expert]);

  const simpleDerived = useMemo(() => simpleScenario(meta, simple), [meta, simple]);
  const scenario = expert ? ex.scn : simpleDerived.scn;
  const debounced = useDebounced(scenario, 300);
  const key = scenarioKey(debounced);
  const run = useAsync(() => runScenario(debounced), key);
  const fresh = run.data && !run.loading;
  // Views render with the scenario that PRODUCED the shown result, never the one still
  // loading -- otherwise a class switch (PFAS -> neutral) reads neutral-only fields off
  // the old PFAS payload for the length of a solve.
  const shownScn: Scenario = run.dataKey ? JSON.parse(run.dataKey) : debounced;

  return (
    <div className={`app ${expert ? "expert" : "simple"}`}>
      <aside className={navOpen ? "sidebar open" : "sidebar"}>
        <Toggle label="🔬 전문가/고급 모드 (Expert / advanced)" checked={expert} onChange={setExpert}
                help="끄면(기본): 쉬운 한국어 화면 (화학물질 + 오염 수준만 선택). 켜면: 전체 연구용 인터페이스(영어) — 6가지 노출 모드, SMILES 구조 입력, 메커니즘 스위치, 모든 모델 파라미터." />
        {expert ? <ExpertSidebar meta={meta} st={ex} set={setEx} run={run.data} onHydrusBuilt={onMetaRefresh} />
                : <SimpleSidebar meta={meta} st={simple} set={setSimple} />}
      </aside>
      <button className="nav-toggle" onClick={() => setNavOpen((o) => !o)} aria-expanded={navOpen}>
        {navOpen ? "✕" : expert ? "⚙ Settings" : "⚙ 설정"}
      </button>
      <main>
        {expert ? <ExpertHeader meta={meta} neutral={ex.scn.compound.kind === "neutral"} /> : <SimpleHeader meta={meta} />}
        {run.error && <Alert kind="error">{run.error}</Alert>}
        {!run.data && run.loading && (
          <Spinner label={expert ? "Running the model…" : "🌾 벼 한 철 축적을 계산하는 중…"} />
        )}
        {run.data && (
          <div className={fresh ? "results" : "results stale"}>
            {run.loading && <div className="stale-badge"><Spinner label={expert ? "Updating…" : "다시 계산 중…"} /></div>}
            <ErrorBoundary resetKey={run.dataKey ?? ""}>
            {expert
              ? <ExpertView meta={meta} st={ex} set={setEx} scn={shownScn} run={run.data} />
              : <SimpleView meta={meta} scn={shownScn} run={run.data} word={simpleDerived.word}
                            label={simpleDerived.label} st={simple} set={setSimple} />}
            </ErrorBoundary>
          </div>
        )}
        <Footer meta={meta} expert={expert} />
      </main>
    </div>
  );
}
