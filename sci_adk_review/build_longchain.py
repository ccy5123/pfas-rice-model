"""
sci-adk long-chain mechanism sub-investigation for the PFAS-rice model.
========================================================================

Frames the long-chain (C10-C12) shortfall as a frozen pre-registration Spec
(proposal_longchain.md) and adjudicates three hypotheses with sci-adk:

  LC1 free-anion loading structurally under-predicts long-chain shoot   -> SUPPORTED
  LC2 a B-independent lipid bound-loading term closes most of the gap    -> SUPPORTED
  LC3 it does so WITHOUT degrading root / the rest (a clean fix)         -> REFUTED

Evidence are the ACTUAL model outputs (validation/longchain_mechanism.py on the
ORYZA2000 biomass + validation/refit_oryza.py), classified measured (vs Yamazaki).
Run:  python sci_adk_review/build_longchain.py
Then: sci-adk verify sci_adk_review/runs/pfas-rice-longchain
"""
from __future__ import annotations

import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from sci_adk.core.spec import (
    Spec, RawProposal, Hypothesis, DecisionRule, MethodPlan, TargetClaim,
    HypothesisMode, DecisionRuleKind,
)
from sci_adk.core.evidence import (
    EvidenceItem, Provenance, Result, Bearing, BearingDirection, EvidenceKind,
)
from sci_adk.core.claim import ConfidenceLevel
from sci_adk.core.parser import ProposalParser
from sci_adk.loop.verdict import VerdictTrail, PanelVerdict, ChiefVerdict, VerdictProvenance
from sci_adk.loop.checkpoint_loop import run_checkpoint_loop
from sci_adk.loop.verify import verify_run

HERE = Path(__file__).resolve().parent
PROPOSAL = HERE / "proposal_longchain.md"
NOW = datetime(2026, 6, 17, tzinfo=timezone.utc)


def _git_commit():
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=HERE.parent, text=True).strip()
    except Exception:
        return "unknown"

CODE_REF = _git_commit()

R_LC1 = ("free-anion xylem/phloem loading structurally under-predicts long-chain "
         "(C10-C12) shoot accumulation -- even at the f_xy=1 ceiling the flux is "
         "throttled by Cw=C/B collapse => support; long chains are reachable by free "
         "loading => refute")
R_LC2 = ("a B-independent lipid-facilitated bound-loading term (g_xy*C, g_ph*C) closes "
         "most of the long-chain shoot gap that free loading cannot => support; no "
         "improvement => refute")
R_LC3 = ("the lipid mechanism reproduces the long-chain shoot WITHOUT degrading the "
         "root or the short/mid-chain fits (a clean single-mechanism fix) => support; "
         "it trades off root / leaves a residual => refute")


def build_spec():
    sec = ProposalParser()._extract_sections(PROPOSAL.read_text(encoding="utf-8"))
    raw = RawProposal(background=sec["background"], goal=sec["goal"],
                      method=sec["method"], expected_output=sec["expected_output"])
    H = lambda i, s, e: Hypothesis(id=i, statement=s, mode=HypothesisMode.CONFIRMATORY,
                                   decision_rule=DecisionRule(kind=DecisionRuleKind.QUALITATIVE,
                                                              expression=e), referent="empirical")
    hyps = [
        H("hyp-lc-freethrottle",
          "Free-anion loading structurally under-predicts long-chain (C10-C12) shoot "
          "accumulation; even at the f_xy=1 ceiling the loaded flux is throttled by the "
          "Cw=C/B collapse, so straw/grain cannot reach the observed values.", R_LC1),
        H("hyp-lc-lipidfix",
          "A B-independent lipid-facilitated bound-loading term (g_xy*C_root, g_ph*C_leaf) "
          "closes most of the long-chain shoot gap that free loading cannot.", R_LC2),
        H("hyp-lc-nocost",
          "The lipid mechanism reproduces the long-chain shoot without degrading the root "
          "or the short/mid-chain fits (a clean single-mechanism fix).", R_LC3),
    ]
    tcs = [TargetClaim(id="tc-lc1", statement="Free loading is the long-chain shortfall cause.",
                       answers="hyp-lc-freethrottle"),
           TargetClaim(id="tc-lc2", statement="Lipid bound-loading closes the long-chain shoot.",
                       answers="hyp-lc-lipidfix"),
           TargetClaim(id="tc-lc3", statement="The lipid fix is cost-free.",
                       answers="hyp-lc-nocost")]
    method = MethodPlan(approaches=["free-only vs lipid on ORYZA2000 biomass vs Yamazaki",
                                    "Cw=C/B collapse diagnosis", "refit f_xy=1 ceiling check"])
    return Spec(id="pfas-rice-longchain", created_at=NOW, version=1, raw_proposal=raw,
                hypotheses=hyps, method=method, target_claims=tcs)


def _ev(id_, ds, point, finding, bears, env):
    return EvidenceItem(id=id_, created_at=NOW, spec_id="pfas-rice-longchain",
                        kind=EvidenceKind.EXPERIMENT_RUN,
                        provenance=Provenance(code_ref=CODE_REF, data_source=ds, environment=env),
                        result=Result(type="quantitative", point=point, finding=finding),
                        bears_on=bears)


def _lit():
    """Prior-work LITERATURE record (sci-adk: discovery via the agent's web_search is
    upstream; this records the result). paperforge OA-PDF acquisition is unavailable
    here (private [tools]), so the DOIs are recorded without local PDFs."""
    dois = [
        {"doi": "10.1021/acs.est.4c06734",
         "note": "Chen2025 ES&T: membrane-water partition +0.36/CF2 (C4-C16) RISES for long "
                 "chains while protein (HSA) PLATEAUS ~C6-C8 -> lipid pool dominates long-chain "
                 "partitioning. Direct mechanistic basis for LC2 (B-independent lipid loading)."},
        {"doi": "10.1021/acs.est.5c11716",
         "note": "Biomimetic chromatography membrane-water + protein-water partition for PFAS (LC2)."},
        {"doi": "10.1021/acs.est.7b06128",
         "note": "Chain-length-dependent tissue distribution (membrane vs protein) in crucian carp."},
        {"doi": "10.1021/acsestengg.4c00107",
         "note": "ML plant uptake/translocation: MW dominates RCF/SCF/TF; PFCA log BCF concave, "
                 "PFSA rises with chain length (LC1 chain-length dependence)."},
        {"doi": "10.48130/newcontam-0025-0007",
         "note": "Soil-plant systems review: Casparian strip restricts long-chain (C>=7 PFCA, "
                 ">=6 PFSA) translocation; long chains root-retained (LC1)."},
        {"doi": "10.1007/s40726-020-00168-y",
         "note": "Current Pollution Reports: PFAS plant uptake by chain length / functional group."},
        {"doi": "10.1139/er-2025-0116",
         "note": "Critical review: PFAS uptake, translocation, toxicity in plants (context)."},
    ]
    finding = json.dumps({
        "acquired": dois, "failed": [],
        "corroboration": "LC1 corroborated (Casparian long-chain translocation barrier; long "
                         "chains root-retained). LC2 mechanism corroborated (membrane/lipid "
                         "partition keeps rising with chain length while protein plateaus, so a "
                         "lipid-facilitated bound pool -- not protein -- carries long chains; "
                         "phospholipids facilitate anion transfer to the lipid phase). NOT a "
                         "novelty claim: these mechanisms are literature-established."})
    return EvidenceItem(id="evi-lc-literature", created_at=NOW, spec_id="pfas-rice-longchain",
                        kind=EvidenceKind.LITERATURE,
                        provenance=Provenance(code_ref=CODE_REF, data_source=None,
                                              environment="agent web_search (paperforge OA acquisition unavailable)"),
                        result=Result(type="qualitative", finding=finding), bears_on=[])


def evidence(spec, ws):
    B = Bearing
    return [
        _lit(),
        _ev("evi-lc-free", "measured", 2.026,
            "validation/longchain_mechanism.py (ORYZA2000 biomass): free-only (monotone "
            "f_xy) long-chain (nC>=10) straw+grain log10 RMSE 2.026 (~100x); PFDA straw "
            "0.08 vs 3.46, PFDoDA straw 0.33 vs 49.8 / grain 0.10 vs 45.5. refit_oryza.py: "
            "even f_xy=1/L_Ph=1 ceilings leave PFDoDA straw 14.6 vs 49.8 -- the Cw=C/B "
            "collapse throttles free loading regardless of f_xy.",
            [B(target_id="hyp-lc-freethrottle", direction=BearingDirection.SUPPORTS)],
            "validation/longchain_mechanism.py + validation/refit_oryza.py (real runs)"),
        _ev("evi-lc-lipid", "measured", 0.428,
            "validation/longchain_mechanism.py: the B-independent lipid term (g_xy*C, "
            "g_ph*C; K_PL-gated) cuts long-chain straw+grain log10 RMSE 2.026 -> 0.428 "
            "(~100x -> ~2.7x; whole series 1.035 -> 0.386). PFDA straw 0.08->5.95 (obs "
            "3.46), PFUnDA straw 0.16->11.15 (obs 8.16), PFOS straw 0.17->5.17 (obs 4.35).",
            [B(target_id="hyp-lc-lipidfix", direction=BearingDirection.SUPPORTS)],
            "validation/longchain_mechanism.py (real run; ORYZA2000)"),
        _ev("evi-lc-cost", "measured", None,
            "validation/longchain_mechanism.py: the single-pool lipid term TRADES OFF "
            "root for long chains -- PFUnDA root 20.6->3.9 (obs 19.5), PFDoDA root "
            "159->4.4 (obs 69); and PFDoDA shoot is STILL ~3-4x under (straw 14.7 vs 49.8). "
            "A clean cost-free fix is refuted: it needs a 2-pool (free + bound) split and a "
            "residual long-chain mechanism.",
            [B(target_id="hyp-lc-nocost", direction=BearingDirection.REFUTES)],
            "validation/longchain_mechanism.py (real run; ORYZA2000)"),
    ]


def _trail(hyp, rub, d, basis, panel):
    return VerdictTrail(hypothesis_id=hyp, rule_kind="qualitative", rubric_expression=rub,
                        panel=panel, chief=ChiefVerdict(direction=d, level=ConfidenceLevel.STRONG, basis=basis),
                        provenance=VerdictProvenance(spec_version=1, timestamp=NOW.isoformat(),
                                                     agent_ids=["claude-in-session"]))


def verdicts():
    SUP, REF, ST, MO = (BearingDirection.SUPPORTS, BearingDirection.REFUTES,
                        ConfidenceLevel.STRONG, ConfidenceLevel.MODERATE)
    return {
        "hyp-lc-freethrottle": _trail(
            "hyp-lc-freethrottle", R_LC1, SUP,
            "Decisive: free-only long-chain straw+grain RMSE 2.026 (~100x under), and the "
            "ORYZA re-fit hits f_xy=1/L_Ph=1 ceilings yet PFDoDA straw stays 14.6 vs 49.8 "
            "-- the loaded flux scales with Cw=C/B, which collapses as B grows, so the "
            "free-anion shoot cannot reach the long chains. Cause confirmed. Literature "
            "corroborates: the Casparian strip restricts long-chain (C>=7 PFCA) translocation "
            "and long chains are root-retained (soil-plant reviews, evi-lc-literature).",
            [PanelVerdict(direction=SUP, level=ST, basis="f_xy=1 ceiling still under -> not an f_xy value issue but the Cw=C/B throttle."),
             PanelVerdict(direction=SUP, level=ST, basis="free-only long-chain RMSE 2.03 (~100x) vs measured Yamazaki; corroborated by the long-chain root-retention literature.")]),
        "hyp-lc-lipidfix": _trail(
            "hyp-lc-lipidfix", R_LC2, SUP,
            "Decisive: the B-independent lipid term cuts the long-chain straw+grain RMSE "
            "2.026 -> 0.428 (~100x -> ~2.7x) and the whole series 1.035 -> 0.386, recovering "
            "PFDA/PFUnDA/PFOS shoot the free model starved. The mechanism is the right "
            "direction. SCOPE: in-sample (K_PL-gated fit), and PFDoDA shoot is still ~3-4x "
            "under -- a residual floor remains. Literature corroborates the MECHANISM "
            "(evi-lc-literature): Chen2025 (ES&T 10.1021/acs.est.4c06734) shows the "
            "membrane-water partition keeps rising +0.36/CF2 for long chains while protein "
            "(HSA) plateaus ~C6-C8, so the lipid (membrane) pool -- not protein -- is what "
            "carries long chains; phospholipids facilitate anion transfer to the lipid phase.",
            [PanelVerdict(direction=SUP, level=ST, basis="long-chain straw+grain RMSE 0.43 vs free 2.03 -- ~5x closer in log space."),
             PanelVerdict(direction=SUP, level=MO, basis="mechanism literature-grounded: membrane partition rises while protein plateaus for long chains (Chen2025).")]),
        "hyp-lc-nocost": _trail(
            "hyp-lc-nocost", R_LC3, REF,
            "Decisive: the single-pool lipid term reproduces shoot at the COST of root -- "
            "PFUnDA root 20.6->3.9, PFDoDA root 159->4.4 (now under vs obs 69/19.5) -- and "
            "PFDoDA shoot is still ~3-4x under. Not a clean cost-free fix; it needs a 2-pool "
            "(free + lipid-bound) split so the bound pool feeds the shoot without draining "
            "the root, plus a residual long-chain mechanism.",
            [PanelVerdict(direction=REF, level=ST, basis="long-chain root degrades (single-pool drains root into the bound shoot flux)."),
             PanelVerdict(direction=REF, level=ST, basis="PFDoDA shoot residual ~3-4x under even with lipid on.")]),
    }


def main():
    spec = build_spec()
    run_dir = HERE / "runs" / spec.id
    run_dir.mkdir(parents=True, exist_ok=True)
    vdir = run_dir / "verdicts"; vdir.mkdir(parents=True, exist_ok=True)
    for h, tr in verdicts().items():
        (vdir / f"{h}.json").write_text(json.dumps(tr.model_dump(mode="json"), indent=2, ensure_ascii=False),
                                        encoding="utf-8")
    res = run_checkpoint_loop(run_dir=run_dir, spec=spec, experiment=evidence, workspace_dir=HERE)
    print(f"=== LONG-CHAIN RUN '{spec.id}' ({res.iterations} iter) ===")
    for c in res.claims:
        print(f"  {c.answers:22s} -> {c.status.value.upper():9s} | {c.confidence.basis[:80]}")
    rep = verify_run(run_dir)
    print(f"\nverify: all_reproduced={rep.all_reproduced}  digest sha256:{rep.digest[:16]}…")
    for o in rep.outcomes:
        print(f"  {o.hypothesis_id:22s} {o.result}")


if __name__ == "__main__":
    main()
