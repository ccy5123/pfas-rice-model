"""Long-chain root->shoot DECOUPLING run -- engine-rendered, built by sci-adk.

Tests whether an irreversible root sequestration (asymmetric bound-store kinetics,
k_on = ratio*k_off*seq) decouples the long-chain root from the shoot -- the missing
piece FINDINGS sec.7 named after the complete-resolution run over-fed the shoot.

Runs the live experiment (validation/longchain_decouple.py), freezes its verified
statistics into a threshold-hypothesis Spec, records measured Evidence, lets the
engine resolve the numeric claims and render the paper (narrative = LaTeX-safe prose
INPUT, not hand-authored, not LLM-generated), and verifies.

Run:    python sci_adk_review/build_longchain_decouple.py
Then:   sci-adk verify sci_adk_review/runs/pfas-rice-longchain-decouple
"""
from __future__ import annotations

import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from sci_adk.core.spec import (
    Spec, RawProposal, Hypothesis, DecisionRule, MethodPlan, TargetClaim,
    HypothesisMode, DecisionRuleKind,
)
from sci_adk.core.evidence import (
    EvidenceItem, Provenance, Result, Bearing, BearingDirection, EvidenceKind,
)
from sci_adk.core.parser import ProposalParser
from sci_adk.loop.compiler import ResearchCompiler
from sci_adk.loop.checkpoint_loop import run_checkpoint_loop
from sci_adk.loop.verify import verify_run
from sci_adk.render.prose import PaperProse

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PROPOSAL = HERE / "proposal_longchain_decouple.md"
NOW = datetime.now(timezone.utc)
SPEC_ID = "pfas-rice-longchain-decouple"

sys.path.insert(0, str(ROOT / "validation"))
sys.path.insert(0, str(ROOT / "src"))
import longchain_decouple as LC  # noqa: E402


def _code_ref() -> str:
    try:
        sha = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=HERE).decode().strip()
        return f"pfas-rice-model@{sha}"
    except Exception:
        return "pfas-rice-model@HEAD"


CODE_REF = _code_ref()

# Live experiment ONCE; thresholds are frozen, the points are the actual outputs.
_RES = LC.run()
_D = _RES["PFDoDA"]                       # the decisive C12 case
_U = _RES["PFUnDA"]                       # corroboration
_SEQ1_GAP = _D["rows"][0][3]             # seq=1 gap = the complete-recipe baseline
_IMPROVE = round(_SEQ1_GAP - _D["min_gap"], 3)

R_INFLATE = ("at a fixed root-matching carrier, increasing the irreversible sequestration seq "
             "INFLATES the long-chain (PFDoDA) root by > 2x (the lever pushes the root the wrong "
             "way) => support; root does not inflate => refute")
R_NOTCLEAN = ("the lever does NOT achieve a clean simultaneous root+shoot closure: the best "
              "(minimum over seq) simultaneity gap stays above 0.30 log10 (factor 2) => support; "
              "gap <= 0.30 (both within a factor 2) => refute")
R_MARGINAL = ("the best simultaneity gap improves on the complete-recipe baseline (seq=1) by less "
              "than 0.05 log10 -- a negligible gain, so the sequestration lever does not materially "
              "help => support; improvement >= 0.05 => refute")


def _raw_proposal() -> RawProposal:
    s = ProposalParser()._extract_sections(PROPOSAL.read_text(encoding="utf-8"))
    return RawProposal(background=s["background"], goal=s["goal"],
                       method=s["method"], expected_output=s["expected_output"])


def build_spec() -> Spec:
    hyps = [
        Hypothesis(
            id="hyp-decouple-inflates-root",
            statement="At a fixed root-matching carrier, increasing the irreversible bound-store "
                      "sequestration (seq) INFLATES the long-chain (PFDoDA) root rather than "
                      "relieving the shoot: root inflates by more than a factor 2 across the scan.",
            mode=HypothesisMode.CONFIRMATORY,
            decision_rule=DecisionRule(kind=DecisionRuleKind.THRESHOLD, expression=R_INFLATE,
                params={"statistic": "point", "op": ">", "value": 2.0, "combine": "latest"}),
            referent="empirical",
        ),
        Hypothesis(
            id="hyp-decouple-not-clean",
            statement="The irreversible-sequestration lever does NOT achieve a clean simultaneous "
                      "root+shoot closure for the long chains: the best (minimum over seq) "
                      "simultaneity gap stays above 0.30 log10 (i.e. root and shoot cannot both be "
                      "brought within a factor 2 at once).",
            mode=HypothesisMode.CONFIRMATORY,
            decision_rule=DecisionRule(kind=DecisionRuleKind.THRESHOLD, expression=R_NOTCLEAN,
                params={"statistic": "point", "op": ">", "value": 0.30, "combine": "latest"}),
            referent="empirical",
        ),
        Hypothesis(
            id="hyp-decouple-marginal",
            statement="The lever helps only marginally: the best simultaneity gap improves on the "
                      "complete-recipe baseline (seq=1) by less than 0.05 log10, so asymmetric "
                      "bound-store kinetics do not materially resolve the long-chain root<->shoot "
                      "tension.",
            mode=HypothesisMode.CONFIRMATORY,
            decision_rule=DecisionRule(kind=DecisionRuleKind.THRESHOLD, expression=R_MARGINAL,
                params={"statistic": "point", "op": "<", "value": 0.05, "combine": "latest"}),
            referent="empirical",
        ),
    ]
    tcs = [
        TargetClaim(id="tc-inflate", statement="Sequestration inflates the long-chain root.",
                    answers="hyp-decouple-inflates-root"),
        TargetClaim(id="tc-notclean", statement="No clean simultaneous closure.",
                    answers="hyp-decouple-not-clean"),
        TargetClaim(id="tc-marginal", statement="The lever helps only marginally.",
                    answers="hyp-decouple-marginal"),
    ]
    method = MethodPlan(approaches=[
        "add a seq lever (k_on = ratio*k_off*seq) to the 2-pool prototype (separate file)",
        "fix the active carrier at its seq=1 root-matching value; scan seq",
        "read the g_xy=0 straw floor (lipid can only add) -> straw fit ratio, simultaneity gap",
        "freeze baseline/inflation/best-gap as thresholds; engine resolves and renders; verify",
    ])
    return Spec(id=SPEC_ID, created_at=NOW, version=1, raw_proposal=_raw_proposal(),
                hypotheses=hyps, method=method, target_claims=tcs)


def _ev(id_, kind, ds, result, bears_on, env):
    return EvidenceItem(id=id_, created_at=NOW, spec_id=SPEC_ID, kind=kind,
                        provenance=Provenance(code_ref=CODE_REF, data_source=ds, environment=env),
                        result=result, bears_on=bears_on)


def evidence(spec, workspace):
    items = []
    items.append(_ev(
        "evi-decouple-literature", EvidenceKind.LITERATURE, None,
        Result(type="qualitative", finding=(
            "cited prior work (LC records): long-chain PFAS are strongly retained in roots vs "
            "mobile short chains (newcontam 2025; Adu 2024 ML 10.1021/acsestengg.4c00107); Chen "
            "2025 membrane K_MW (10.1021/acs.est.4c06734). Follow-up to runs/pfas-rice-longchain-"
            "complete; not a novelty claim.")),
        [], "follow-up to runs/pfas-rice-longchain-complete"))

    items.append(_ev(
        "evi-decouple-inflate", EvidenceKind.EXPERIMENT_RUN, "measured",
        Result(type="quantitative", point=round(_D["root_inflation"], 2), finding=(
            f"validation/longchain_decouple.py vs Yamazaki: at the fixed root-matching carrier "
            f"({_D['vm_x_base']:.1f}x base), increasing seq 1->10 INFLATES the PFDoDA root "
            f"{_D['root_inflation']:.2f}x (root ratio {_D['rows'][0][1]:.2f} -> {_D['root_ratio_seqmax']:.2f}) "
            f"while the shoot over-feed only relaxes (straw {_D['straw_ratio_seq1']:.2f}x -> "
            f"{_D['rows'][-1][2]:.2f}x). Mechanism: suppressing the mobile pool lowers the internal "
            f"conc, which RAISES the net uptake gradient, so the bound store grows -- the lever "
            f"pushes the root the wrong way. PFUnDA corroborates (root inflation {_U['root_inflation']:.2f}x).")),
        [Bearing(target_id="hyp-decouple-inflates-root", direction=BearingDirection.SUPPORTS)],
        "validation/longchain_decouple.py"))

    items.append(_ev(
        "evi-decouple-notclean", EvidenceKind.EXPERIMENT_RUN, "measured",
        Result(type="quantitative", point=round(_D["min_gap"], 3), finding=(
            f"validation/longchain_decouple.py: the best (minimum over seq) simultaneity gap for "
            f"PFDoDA is {_D['min_gap']:.3f} log10 (factor {10**_D['min_gap']:.2f}) at seq={_D['best_seq']:.1f} "
            f"(root {_D['best_root_ratio']:.2f}x, straw {_D['best_straw_ratio']:.2f}x) -- above 0.30, so root "
            f"and shoot CANNOT both be brought within a factor 2 at once. PFUnDA best gap "
            f"{_U['min_gap']:.3f} (factor {10**_U['min_gap']:.2f}).")),
        [Bearing(target_id="hyp-decouple-not-clean", direction=BearingDirection.SUPPORTS)],
        "validation/longchain_decouple.py"))

    items.append(_ev(
        "evi-decouple-marginal", EvidenceKind.EXPERIMENT_RUN, "measured",
        Result(type="quantitative", point=_IMPROVE, finding=(
            f"validation/longchain_decouple.py: seq=1 IS the complete recipe (gap {_SEQ1_GAP:.3f}); the "
            f"best seq lowers the PFDoDA gap only to {_D['min_gap']:.3f}, an improvement of {_IMPROVE:.3f} "
            f"log10 (< 0.05) -- negligible. So asymmetric bound-store kinetics do NOT materially "
            f"resolve the long-chain root<->shoot tension. The correct decoupling must break the "
            f"uptake<->mobile-conc coupling (e.g. deposit uptake into an inert apoplastic store that "
            f"does not raise the influx gradient), a now-sharper open target.")),
        [Bearing(target_id="hyp-decouple-marginal", direction=BearingDirection.SUPPORTS)],
        "validation/longchain_decouple.py"))

    return items


PROSE = PaperProse(
    abstract=(
        "The complete-resolution run (runs/pfas-rice-longchain-complete) closed the long-chain "
        "root and grain but over-fed the shoot, because the active carrier that fixes the root "
        "enlarges the mobile root pool whose free xylem loading exceeds the observed straw. "
        "FINDINGS section 7 named the missing piece a root-to-shoot decoupling: an irreversible "
        "root store that holds the root burden without feeding the xylem. This run tests the "
        "simplest such lever on the 2-pool prototype -- an asymmetric bound-store kinetics factor "
        "seq (k_on = ratio times k_off times seq) that traps more burden in the non-translocating "
        "bound pool. At a fixed root-matching carrier we scan seq against the measured Yamazaki "
        "long chains. Result, honest: the lever does NOT decouple. Suppressing the mobile pool to "
        "protect the shoot lowers the internal concentration, which RAISES the net uptake gradient, "
        "so the root burden inflates with seq (PFDoDA root 1x to 6.9x at seq=10) faster than the "
        "shoot is relieved (straw 2.3x to 1.6x). The best balanced point (seq about 2) leaves both "
        "root and shoot about 2.2x off -- a simultaneity gap of about 0.34 log10, only about 0.02 "
        "better than the complete recipe. So asymmetric bound-store kinetics trade root for shoot "
        "but do not achieve a clean simultaneous closure. The correct decoupling must break the "
        "uptake-to-mobile-concentration coupling itself."
    ),
    introduction=(
        "This is a direct follow-up to the long-chain complete-resolution test. That run reduced "
        "FINDINGS section 7's proposed three-lever resolution to a precise open problem -- the "
        "carrier that fixes the root structurally over-feeds the shoot -- and named the missing "
        "mechanism a root-to-shoot decoupling. sci-adk's discipline is to test a named mechanism, "
        "not assume it. Here the mechanism is the simplest concrete form of decoupling available in "
        "the 2-pool prototype: make the bound store irreversible so it holds root burden without "
        "translocating. The experiment reuses the committed prototype unchanged and adds only the "
        "seq lever in a separate file; the engine resolves the numeric claims from the live run. "
        "The narrative is agent-authored prose input; the engine renders. No core model is changed."
    ),
    discussion=(
        "The verdict is an honest near-negative that sharpens the open problem. The asymmetric "
        "bound-store lever does provide a knob -- increasing seq traps more burden in the "
        "non-translocating pool and does relieve the shoot somewhat -- but it cannot deliver a "
        "clean simultaneous closure, because the model couples uptake to the mobile concentration: "
        "the Goldman-Hodgkin-Katz and carrier influx grow as the internal (mobile) concentration "
        "falls, so any attempt to suppress the mobile pool (to protect the shoot) increases net "
        "uptake and therefore inflates the root. The two objectives -- low mobile concentration "
        "(for the shoot) and a finite measured root burden -- are in direct tension under this "
        "coupling. The best the lever achieves is a balanced compromise around a factor 2.2 on each "
        "tissue, a negligible 0.02-log improvement over the complete recipe. The implication is "
        "specific and constructive: the correct root-to-shoot decoupling must break the "
        "uptake-to-mobile-concentration coupling itself -- for example, route long-chain uptake "
        "directly into an inert apoplastic or cell-wall store that holds the measured root burden "
        "WITHOUT contributing to the mobile aqueous concentration that both feeds the xylem and "
        "sets the influx gradient. That is a different object from asymmetric kinetics on an "
        "equilibrating bound pool, and it is the now-sharper target for closing the longest chains. "
        "This is a prototype, in-sample synthesis test; the core model is unchanged, and the three "
        "claims reproduce under sci-adk verify. It consolidates into FINDINGS.md section 7."
    ),
)


def main():
    spec = build_spec()
    run_dir = HERE / "runs" / spec.id
    run_dir.mkdir(parents=True, exist_ok=True)
    run_checkpoint_loop(run_dir=run_dir, spec=spec, experiment=evidence, workspace_dir=HERE)
    result = ResearchCompiler(workspace_dir=HERE).compile(
        PROPOSAL.read_text(encoding="utf-8"), spec=spec, experiment=evidence, prose=PROSE)
    print(f"=== LONG-CHAIN DECOUPLE RUN  '{result.spec.id}' ===")
    print(f"  PFDoDA: root inflation {_D['root_inflation']:.2f}x | best gap {_D['min_gap']:.3f} "
          f"(factor {10**_D['min_gap']:.2f}) at seq {_D['best_seq']:.1f} | improvement over complete {_IMPROVE:.3f}")
    print(f"  evidence: {len(result.evidence)} | claims: {len(result.claims)}")
    for c in result.claims:
        print(f"  {c.answers:28s} -> {c.status.value.upper():9s} | {c.confidence.basis[:62]}")
    print(f"  paper: {result.paper_path}")
    report = verify_run(run_dir)
    print(f"\n=== VERIFY (digest sha256:{report.digest[:16]}...) -> all reproduced: {report.all_reproduced} ===")
    for o in report.outcomes:
        rd = o.rederived_status.value if o.rederived_status else "n/a"
        print(f"  {o.hypothesis_id:28s} -> {o.result:11s} (recorded={o.recorded_status.value}, re-derived={rd})")


if __name__ == "__main__":
    main()
