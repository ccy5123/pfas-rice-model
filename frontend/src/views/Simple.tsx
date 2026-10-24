// General-audience (Korean) view -- the port of ui/simple.py + the Simple sidebar.

import { useState } from "react";
import { download, getFigure, svgToPng } from "../api";
import { CustomTables } from "../components/CustomTables";
import { InversePanel } from "../components/InversePanel";
import { PlantMap } from "../components/PlantMap";
import { Plot, ServerFigure } from "../components/Plot";
import { Alert, Button, Caption, Checkbox, Columns, Expander, Md, Metric, Radio, Select, Slider, Spinner, Tabs } from "../components/ui";
import { fmtF, fmtG } from "../format";
import { useAsync } from "../hooks";
import { congenerLabel, signalBadge, t } from "../i18n";
import { plantSvgFromRun } from "../plantSvg";
import type { Meta, RunResult, Scenario, TablesSpec } from "../types";
import { defaultScenario } from "../types";

export const MANUAL = "✏️ 직접 설정";

export interface SimpleState {
  scen: string;          // a one-click scenario label, or MANUAL
  congener: string;
  preset: string;        // a contamination preset label
  useTables: boolean;
  tables: TablesSpec;
}

export function initialSimple(meta: Meta): SimpleState {
  return { scen: meta.scenarios_ko[0][0], congener: "PFOA", preset: meta.presets.ko[1][0],
           useTables: false, tables: defaultScenario("ko").exposure.tables };
}

/** Simple sidebar state -> (scenario, short level word, label for the handout). */
export function simpleScenario(meta: Meta, st: SimpleState): { scn: Scenario; word: string; label: string } {
  const s = defaultScenario("ko");
  let cong: string, Cwo: number, word: string;
  if (st.scen === MANUAL) {
    const p = meta.presets.ko.find(([l]) => l === st.preset) ?? meta.presets.ko[1];
    cong = st.congener; Cwo = p[1]; word = p[2];
  } else {
    const sc = meta.scenarios_ko.find(([l]) => l === st.scen) ?? meta.scenarios_ko[0];
    cong = sc[1]; Cwo = sc[2]; word = sc[3];
  }
  s.compound.congener = cong;
  s.exposure.Cwo = Cwo;
  if (st.useTables) { s.exposure.mode = "tables"; s.exposure.tables = st.tables; }
  return { scn: s, word, label: st.scen === MANUAL ? st.preset : st.scen };
}

export function SimpleSidebar({ meta, st, set }: { meta: Meta; st: SimpleState; set: (s: SimpleState) => void }) {
  const scenLabels = [MANUAL, ...meta.scenarios_ko.map(([l]) => l)];
  const sc = meta.scenarios_ko.find(([l]) => l === st.scen);
  return (
    <>
      <h3>① 시나리오</h3>
      <Radio label="빠른 시나리오" cards options={scenLabels} value={st.scen} onChange={(v) => set({ ...st, scen: v })}
             help="정책 상황별로 화학물질과 오염 수준을 한 번에 설정합니다. '직접 설정'을 고르면 물질·오염도를 각각 지정할 수 있습니다." />
      {st.scen === MANUAL ? (
        <>
          <h3>② 화학물질 선택</h3>
          <Select label="PFAS 화학물질" options={meta.congeners} value={st.congener}
                  format={(c) => congenerLabel(meta, c, "ko")} onChange={(v) => set({ ...st, congener: v })}
                  help="특정 '영원한 화학물질' 하나. PFOA·PFOS가 가장 잘 알려져 있고, 사슬이 길수록 대체로 식물에 더 잘 달라붙습니다." />
          <h3>③ 오염 정도</h3>
          <Radio label="논의 오염 수준" options={meta.presets.ko.map(([l]) => l)} value={st.preset}
                 onChange={(v) => set({ ...st, preset: v })}
                 help="토양수에 녹아 있는 PFAS의 양. 높을수록 식물로 더 많이 들어갑니다." />
        </>
      ) : sc && (
        <>
          <Alert kind="success">{`**${sc[0]}**\n\n${sc[4]}`}</Alert>
          <Caption>{`→ 물질 **${congenerLabel(meta, sc[1], "ko")}** · 오염 **${sc[3]}** (${fmtG(sc[2])} µg/L). '✏️ 직접 설정'에서 바꿀 수 있습니다.`}</Caption>
        </>
      )}
      <Checkbox label="📋 내 데이터 표 사용" checked={st.useTables} onChange={(v) => set({ ...st, useTables: v })}
                help="성장 곡선과 시간에 따른 토양수 오염 수준을 편집 가능한 표로 직접 입력합니다 (본문에 나타납니다)." />
      <Caption>전체 연구용 인터페이스가 필요하면 위의 **전문가/고급 모드**를 켜세요.</Caption>
    </>
  );
}

const TABS = ["🗺️ 어디로 가나", "📈 시간에 따른 축적", "📊 얼마나 쌓이나", "⚖️ 물질 비교",
              "🍚 밥으로 먹으면", "🔎 거꾸로 추정", "🔬 신뢰도 & 안내"];

export function SimpleHeader({ meta }: { meta: Meta }) {
  return (
    <>
      <div className="hero">
        <div className="hero-kicker">논 PFAS 벼 축적 시뮬레이터</div>
        <h1 className="hero-title">🌾 벼의 어디에, 얼마나 쌓일까?</h1>
        <div className="hero-sub">논의 물·흙에 녹은 <b>'영원한 화학물질' PFAS</b>가 벼의 <b>뿌리·짚·먹는 낟알</b>에 얼마나
          쌓이는지 한눈에 보여줍니다.</div>
      </div>
      <div className="disc">⚠️ 연구·교육용 예시 추정치입니다. 규제·식품안전·건강 판단이 <b>아니며</b>, 실제 노출·안전 결정에 사용하지 마세요.</div>
      <Caption>{t(meta, "header.intro2", "ko")}</Caption>
    </>
  );
}

function IntakeTab({ meta, scn, run }: { meta: Meta; scn: Scenario; run: RunResult }) {
  const [rice, setRice] = useState(Math.round(meta.intake.rice_g));
  const [bw, setBw] = useState(Math.round(meta.intake.bw));
  const opts = { rice_g: rice, bw };
  const key = JSON.stringify([scn, rice, bw]);
  const { data } = useAsync(() => getFigure(scn, "intake_gauge", "ko", opts), key);
  const info = data?.data.info;
  const grain = run.finals.grain;
  const badge = info ? signalBadge(info.signal, true,
    info.in_group ? `기준의 ${fmtF(info.percent, 0)}%` : "참고") : null;
  return (
    <div>
      <h4>이 쌀을 먹으면 안전기준 대비 어느 정도일까요? (참고)</h4>
      <Columns n={2}>
        <Slider label="하루 쌀 섭취량 [g/일]" min={30} max={400} step={10} value={rice} onChange={setRice}
                help="국민건강영양조사·KOSIS 기준 국내 1인 1일 쌀 섭취량은 약 150 g입니다." />
        <Slider label="체중 [kg]" min={20} max={100} step={5} value={bw} onChange={setBw}
                help="1인 기준 체중. 어린이 등 체중이 작을수록 체중당 섭취 비율이 높아집니다." />
      </Columns>
      {badge && <span className={`badge ${badge.cls}`}>{badge.text}</span>}
      {data ? <Plot figure={data.figure} /> : <Spinner label="계산 중…" />}
      {info && (
        <>
          <div className="cols" style={{ gridTemplateColumns: "1fr 2fr" }}>
            <Metric label="주간 섭취량 / 안전기준" value={`${fmtF(info.percent, 0)}%`}
                    sub={`대략 ${fmtF(data!.data.lo_percent, 0)}–${fmtF(data!.data.hi_percent, 0)}%`} />
            <Md>{`낟알 **${fmtG(grain, 2)} µg/kg**인 쌀을 하루 **${rice} g** 드시면 PFAS 주간 섭취량은 EFSA 안전기준(**${fmtG(info.twi)} ng/kg 체중/주**, 4종 합)의 약 **${fmtF(info.percent, 0)}%** 수준입니다.`}</Md>
          </div>
          {info.in_group === false && (
            <Alert kind="warning">{`⚠ **${run.congener}**는 EFSA 4종(PFOA·PFOS·PFNA·PFHxS)에 **미포함** — 위 %는 같은 기준을 빌린 **참고치**입니다.`}</Alert>
          )}
        </>
      )}
      <Expander title="ℹ️ 이 비교를 읽는 법 (중요)">
        <Md>{"- 쌀에 대한 **법적 PFAS 기준(MRL)은 EU·한국 모두 없습니다.** 이 화면은 식품기준 초과가 아니라 **건강기반 섭취기준(EFSA TWI) 환산 참고치**입니다.\n" +
             "- EFSA 기준은 **4종 합**인데 여기선 한 물질만 비교 → **‘이 물질 하나만으로도’ 보수적 가정**.\n" +
             "- **쌀만 섭취원으로 가정**한 단순 계산 (실제론 물·타 식품에서도 섭취).\n" +
             "- 낟알 농도 예측 자체가 실측과 **수배 차이날 수 있습니다.**"}</Md>
      </Expander>
      <Caption>{"출처: EFSA 2020 그룹 TWI 4.4 ng/kg 체중/주 (PFOA·PFOS·PFNA·PFHxS, doi:10.2903/j.efsa.2020.6223) · 쌀 섭취량 KOSIS/국민건강영양조사. 🔴빨강 = 기준의 100% 이상, 🟡노랑 = 10–100%, 🟢초록 = 10% 미만."}</Caption>
    </div>
  );
}

export function SimpleView({ meta, scn, run, word, label, st, set }: {
  meta: Meta; scn: Scenario; run: RunResult; word: string; label: string;
  st: SimpleState; set: (s: SimpleState) => void;
}) {
  const [tab, setTab] = useState(0);
  const [busy, setBusy] = useState<string | null>(null);
  const fold = run.uncertainty_factor;
  const { root, straw, grain } = run.finals;
  const tops: [string, number][] = [["뿌리", root], ["짚(줄기+잎)", straw], ["낟알", grain]];
  const whereKo = tops.reduce((a, b) => (b[1] > a[1] ? b : a))[0];
  const intake = run.intake;
  const badge = signalBadge(intake.signal, true, intake.in_group ? `기준의 ${fmtF(intake.percent, 0)}%` : "참고");
  const cls = badge?.cls ?? "";
  const rng = (b: { lo: number | null; hi: number | null }) => `대략 ${fmtG(b.lo, 2)}–${fmtG(b.hi, 2)} µg/kg`;
  const grainSignal = intake.in_group ? intake.signal : null;
  const dl = async (kind: string, lbl?: string | null) => {
    setBusy(kind);
    try { await download(kind, scn, "ko", lbl); } finally { setBusy(null); }
  };

  return (
    <div>
      {st.useTables && (
        <Expander open title={t(meta, "ct.panel_title", "ko")}>
          <CustomTables meta={meta} lang="ko" tables={st.tables} onChange={(tb) => set({ ...st, tables: tb })}
                        season={scn.exposure.season} biomass={scn.mechanism.biomass} Cwo={scn.exposure.Cwo} run={run} />
        </Expander>
      )}
      <Expander title="💡 처음이신가요? 3단계로 해보세요">
        <Md>{"1. 왼쪽에서 **빠른 시나리오**를 하나 고르세요 (또는 '✏️ 직접 설정'에서 물질·오염도 선택).\n" +
             "2. 아래 **🗺️ 어디로 가나** 지도에서 색이 진한 부위 = PFAS가 많이 쌓이는 곳을 보세요.\n" +
             "3. **🍚 밥으로 먹으면** 탭에서 안전기준 대비 수준을, **⚖️ 물질 비교**에서 물질별 차이를 확인하세요.\n\n" +
             "모든 수치는 **대략적 예측**입니다 — 부위·물질 간 **상대 비교**로 읽어 주세요."}</Md>
      </Expander>
      <h2>{congenerLabel(meta, run.congener, "ko")}</h2>
      <div className={`result ${cls}`}>
        {badge && <span className={`badge ${badge.cls}`}>{badge.text}</span>}
        <div className="result-text">대부분의 PFAS는 <b>{whereKo}</b>에 남고, 먹는 <b>낟알</b>엔 약{" "}
          <b>{fmtG(grain, 2)} µg/kg</b> 들어 있을 것으로 추정합니다.</div>
      </div>
      <Columns n={3}>
        <Metric label="🥕 뿌리 속" value={`${fmtG(root, 2)} µg/kg`} sub={rng(run.bands.root)} />
        <Metric label="🌿 짚(줄기+잎) 속" value={`${fmtG(straw, 2)} µg/kg`} sub={rng(run.bands.straw)} />
        <Metric label="🌾 낟알(먹는 쌀) 속" value={`${fmtG(grain, 2)} µg/kg`} sub={rng(run.bands.grain)} />
      </Columns>
      <Caption>{`신호등은 EFSA 건강기반 안전기준 대비 수준입니다(색·모양·라벨로 표시). 숫자는 **대략적 예측**(실측과 ~${fmtF(fold, 0)}배 차이 가능)이며 법적 식품기준이 아닙니다. 자세한 비교는 🍚·⚖️·🔬 탭에서.`}</Caption>

      <Tabs names={TABS} active={tab} onChange={setTab} />
      <div className="tab-body">
        {tab === 0 && (
          <div>
            <PlantMap run={run} lang="ko" grainSignal={grainSignal} />
            <span className="caveat">ⓘ 대략적 예측 · 실측과 ~{fmtF(fold, 0)}배 차이 가능</span>
            <Caption>{"**색이 진할수록 그 부위에 PFAS가 많습니다.** 카드 숫자는 그 부위 **1 kg당 PFAS 농도**(µg/kg)이고, 맨 아래 **공극수**는 논물·흙물 1 L당 농도(µg/L)입니다. 날짜 슬라이더로 한 철 축적을 보세요. (짚 = 줄기+잎 평균)"}</Caption>
          </div>
        )}
        {tab === 1 && (
          <div>
            <ServerFigure scenario={scn} name="buildup_plain" lang="ko" />
            <Caption>{"한 철 동안 각 식물 부위의 PFAS 농도가 어떻게 변하는지. **낟알**은 형성된 뒤(개화 무렵)부터 PFAS를 흡수하기 시작해 수확까지 계속 쌓입니다. 곡선은 **대략적 모델 예측**으로, 실측과 수배 차이날 수 있습니다."}</Caption>
          </div>
        )}
        {tab === 2 && (
          <div>
            <ServerFigure scenario={scn} name="where_plain" lang="ko" />
            <Caption>{`수확 시 각 부위의 PFAS 농도. 보통 뿌리에 가장 많이 남고, 먹는 낟알까지 얼마나 도달하는지는 화학물질에 따라 다릅니다. 막대의 **회색 오차선**은 모델의 대략적 예측 불확실성(실측과 약 ${fmtF(fold, 0)}배까지 차이날 수 있음)입니다.`}</Caption>
            {Object.keys(run.obs).length > 0 && (
              <Expander title="🔬 실제 측정값과 비교 (Yamazaki 2023)">
                <ServerFigure scenario={scn} name="baf" lang="ko" />
                <Caption>막대는 모델의 축적 배수를 출판된 온실 벼 연구(Yamazaki et al. 2023)의 측정값과 비교한 것입니다. 막대가 비슷할수록 이 화학물질에 대해 모델이 실제 데이터와 잘 맞습니다.</Caption>
              </Expander>
            )}
          </div>
        )}
        {tab === 3 && (
          <div>
            <h4>물질에 따라 쌓이는 곳이 다릅니다</h4>
            <ServerFigure scenario={scn} name="congener_compare" lang="ko" />
            <Caption>{`같은 오염 수준(토양수 ${fmtG(scn.exposure.Cwo)} µg/L)에서 물질별 수확기 농도. **짧은사슬**(PFBA·PFHxA)은 짚·낟알까지 잘 이동해 먹는 부위에 도달하기 쉽고, **긴사슬**(PFDA·PFDoDA)은 뿌리에 강하게 잔류해 낟알로는 덜 갑니다. PFOS는 대표적 술폰산 계열입니다. 값 차이가 매우 커서 세로축은 **로그(배수) 눈금**이며, 막대 높이는 **대략적 모델 예측**입니다.`}</Caption>
          </div>
        )}
        {tab === 4 && <IntakeTab meta={meta} scn={scn} run={run} />}
        {tab === 5 && <InversePanel meta={meta} scenario={scn} lang="ko" />}
        {tab === 6 && (
          <div>
            <h3>🔬 이 모델을 얼마나 믿을 수 있나요?</h3>
            <Md>이 모델은 **정성적 경향(어디에 더 쌓이는가)** 은 잘 맞지만 **절대 수치는 대략치**입니다. 발표·정책 판단 시 아래 검증 상태를 함께 봐 주세요.</Md>
            <Md>{"| 항목 | 상태 | 근거 / 설명 |\n|---|---|---|\n" +
                 "| 음이온 배제(뿌리 흡수 장벽) · 질량보존 | 🟢 **구조적으로 검증** | 막전위 기반 물리(eᴺ≈107); 모델이 반드시 만족 |\n" +
                 "| 짧은사슬 PFAS 이행 **방향** | 🟢 **실측과 일치** | Yamazaki 2023 온실 벼 |\n" +
                 "| 토양→공극수 결합 | 🟢 **실제 HYDRUS-1D 엔진 연동** | Method A (일방향) |\n" +
                 `| **절대 농도**(예측력) | 🟡 **대략치 (실측과 약 ×${fmtF(fold, 0)})** | 미보정 예측 오차 log10 RMSE ≈0.85 |\n` +
                 "| 긴사슬(C10~C12) 짚·낟알 축적 | 🔴 **탐색적·미해결** | 문헌·메커니즘 공백 |\n" +
                 "| 물질별 이행계수 f_xy 절대값 | 🟡 **데이터셋 의존** | 논·품종·조건마다 다름 |\n"}</Md>
            <Caption>🟢 검증됨 · 🟡 대략치/조건의존 · 🔴 미해결(연구 진행 중). 절대 수치보다 **부위 간 상대 비교**와 **물질 간 경향**을 신뢰하세요.</Caption>
            {Object.keys(run.obs).length > 0 ? (
              <>
                <Md>**실제 측정값과의 비교 (Yamazaki 2023 온실 벼)** — 막대가 비슷할수록 이 물질에 대해 모델이 데이터와 잘 맞습니다.</Md>
                <ServerFigure scenario={scn} name="baf" lang="ko" />
              </>
            ) : <Caption>{`${run.congener}에 대한 공개 실측 BAF가 없어 이 물질은 모델 예측만 표시됩니다.`}</Caption>}
            <Alert kind="info">{"**정직한 한계** — 이 모델의 검증은 대부분 **한 실험(Yamazaki, in-sample)** 에 맞춘 것이고, 독립 데이터에 대한 예측 검증은 제한적입니다(단일 clean 데이터셋). 결정적 확증에는 아직 **현장 실험 데이터**(부위별 물관액/토양수 비율 등)가 필요합니다. 따라서 **정책 스크리닝·우선순위·교육용**으로는 유용하지만, **개별 지점의 규제 판단 근거로 단독 사용은 부적절**합니다."}</Alert>
            <hr />
            <Md>{"### 이 도구가 하는 일\n벼가 흙에서 물과 녹아 있는 화학물질을 빨아들여 줄기 위로 올리고 낟알에 저장하는 과정을 **메커니즘 모델**로 계산합니다. PFAS 화학물질과 논의 오염 정도를 주면, 한 철 동안 뿌리·짚·먹는 낟알에 쌓이는 양을 추정합니다.\n\n" +
                 "### 보는 법\n- **🗺️ 어디로 가나** — 식물 그림; 색이 뜨거울수록 PFAS가 많음.\n- **📈 시간에 따른 축적** — 이앙부터 수확까지 농도 변화.\n- **📊 얼마나 쌓이나** — 최종 농도와 실제 측정값과의 비교.\n- **⚖️ 물질 비교** — 여러 PFAS 물질이 뿌리·짚·낟알에 쌓이는 정도를 한눈에.\n- **🍚 밥으로 먹으면** — 예측된 낟알 농도를 EFSA 건강기반 섭취기준에 환산한 **참고** 비교.\n- **🔎 거꾸로 추정** — 실험실 측정값이 있으면 토양수 오염도를 (불확실성 범위와 함께 — 베이지안 추정) 역추정.\n- **🔬 신뢰도 & 안내** — 지금 이 화면. 모델 검증 상태와 용어.\n\n### 쉬운 용어 사전"}</Md>
            <Md>{t(meta, "glossary", "ko")}</Md>
            <Alert kind="warning">{t(meta, "disclaimer", "ko")}</Alert>
            <Caption>수식·파라미터·토양 결합(HYDRUS-1D)·구조(SMILES) 입력이 필요하면 사이드바의 **전문가/고급 모드**를 켜세요.</Caption>
          </div>
        )}
      </div>

      <Expander title="⬇️ 결과 내려받기 · 발표용 요약">
        <div className="downloads">
          <Button onClick={() => dl("summary_html", label)} disabled={!!busy}>🖼️ 발표용 요약 1장 (HTML)</Button>
          <Button onClick={() => dl("summary_md", word)} disabled={!!busy}>📄 요약 텍스트 (Markdown)</Button>
          <Button onClick={() => dl("summary_csv")} disabled={!!busy}>요약 표 (CSV)</Button>
          <Button onClick={() => dl("timeseries_csv")} disabled={!!busy}>전체 시계열 (CSV)</Button>
          <Button disabled={!!busy} onClick={() => svgToPng(plantSvgFromRun(run, -1, { lang: "ko", grainSignal }),
                                                           `${run.congener}_map.png`)}>식물 지도 (PNG)</Button>
          <Button onClick={() => dl("map_html")} disabled={!!busy}>식물 지도 (대화형 HTML)</Button>
        </div>
        {busy && <Spinner label="파일 만드는 중…" />}
        <Caption>발표용 요약은 식물 지도 + 안전기준 게이지 + 핵심 수치 + 한계를 한 파일로 담습니다. 브라우저에서 열어 PDF로 인쇄하면 슬라이드 한 장이 됩니다.</Caption>
      </Expander>
    </div>
  );
}
