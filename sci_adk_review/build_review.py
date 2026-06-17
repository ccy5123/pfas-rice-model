"""
sci-adk rigor review of the PFAS–rice compartmental uptake model.
=================================================================

This script runs the PFAS–rice research proposal (``proposal.md``) through the
**sci-adk** rigor/verification ADK (https://github.com/ccy5123/sci-adk) as a
*referee*: it freezes the four-pane proposal into a pre-registration Spec, records
the model's results as honestly-classified Evidence, has the in-session agent author
chief-over-N Verdicts, and lets sci-adk's DecisionEngine + evidence-validity gate
decide which Claims survive.

WHY this exists
---------------
sci-adk was built to fix what its own source calls the "rice-failure defect"
(``src/sci_adk/core/validity.py``): *"a run on an EMPIRICAL proposal used SYNTHETIC
data and the harness reported 4/4 SUPPORTED."* Its Provenance docstring names
"the rice numbers" as the canonical ``synthetic_proxy`` example. The four-pane
proposal (Background / Goal / Method / Expected Output) is exactly sci-adk's Spec
input — so this review applies the tool to the very project that motivated it.

The honest finding (see FINDINGS.md): the model's STRUCTURE is sound (mass
conservation + anion exclusion → SUPPORTED, formal), but its EMPIRICAL predictive
claims do NOT survive the rigor gate (Yamazaki agreement is a saturated in-sample
fit; grain is structurally under-predicted → REFUTED). The bundled demo BAFs are
``synthetic_proxy`` and the validity gate HALTS rather than certify them (the trap
run reproduces — and now refuses — the original rice-failure).

Run:  python sci_adk_review/build_review.py
Then: sci-adk verify sci_adk_review/runs/pfas-rice
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
from sci_adk.core.validity import ValidityHalt
from sci_adk.loop.verdict import VerdictTrail, PanelVerdict, ChiefVerdict, VerdictProvenance
from sci_adk.loop.checkpoint_loop import run_checkpoint_loop
from sci_adk.loop.verify import verify_run

HERE = Path(__file__).resolve().parent
PROPOSAL = HERE / "proposal.md"

# Provenance anchor: the pfas-rice-model commit these results were produced from.
def _git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=HERE.parent, text=True
        ).strip()
    except Exception:
        return "unknown"

CODE_REF = _git_commit()
NOW = datetime(2026, 6, 17, tzinfo=timezone.utc)   # fixed for a reproducible record digest


def _raw_proposal() -> RawProposal:
    """Use sci-adk's own parser to lift the four panes out of proposal.md verbatim."""
    sections = ProposalParser()._extract_sections(PROPOSAL.read_text(encoding="utf-8"))
    return RawProposal(
        background=sections["background"],
        goal=sections["goal"],
        method=sections["method"],
        expected_output=sections["expected_output"],
    )


# --- frozen DecisionRule expressions (rubric R; copied into each verdict trail) -----
R_MASS = ("max abs relative mass-balance residual over the season is < 1e-5 => support "
          "(mass-conserving to solver tolerance); >= 1e-5 => refute")
R_ANION = ("GHK electrodiffusion with measured E_m=-120 mV and z=-1 yields an "
           "anion-exclusion factor e^N of order 1e2 (the IOC-formulation signature) "
           "=> support; no exclusion => refute")
R_YAMA = ("the Yamazaki agreement is OUT-OF-SAMPLE predictive validation (independent "
          "data not used to fit the per-congener transport parameters) => support; an "
          "in-sample / saturated reproduction => refute")
R_GRAIN = ("model grain (brown-rice) BAF/TF matches measured within a factor adequate "
           "for dietary risk assessment => support; systematic structural "
           "under/over-prediction => refute")
R_SOIL = ("per-congener soil sorption (Koc chain-length QSPR) spans > 2 log10 units "
          "across the congener set, so a single constant pore-water Cwo cannot "
          "represent all congeners (the soil-coupling rationale) => support; < 2 log10 "
          "spread => refute")
R_SMILES = ("the SMILES front-end reproduces the curated measured-parameter model for a "
            "KNOWN PFAS via read-across (structure -> same Compound) => support; a "
            "mismatch => refute")
R_ADEQ = ("a CONSTRAINED (degrees-of-freedom > 0) calibration of the structure "
          "reproduces the measured STRAW (shoot translocation) BAFs across the C4-C12 "
          "PFCA/PFSA series => support; only a saturated (DOF 0) per-congener fit "
          "reproduces => refute")
R_DEMO = ("the bundled demonstration BAFs establish that the model predicts field PFAS "
          "bioaccumulation in rice => support")


# ============================================================================
# 1) The honest Spec  (the heuristic parser cannot infer referent class / rule
#    kind / non-circularity, so — exactly as the README's "capability supplies a
#    pre-built Spec" path intends — we author them deliberately from proposal.md.)
# ============================================================================
def build_spec() -> Spec:
    raw = _raw_proposal()
    hypotheses = [
        Hypothesis(
            id="hyp-mass",
            statement="The four-compartment uptake ODE conserves PFAS mass: total "
                      "compartment burden equals cumulative net root influx minus "
                      "xylem/phloem exports and growth dilution, to solver tolerance.",
            mode=HypothesisMode.CONFIRMATORY,
            decision_rule=DecisionRule(
                kind=DecisionRuleKind.THRESHOLD, expression=R_MASS,
                params={"statistic": "point", "op": "<", "value": 1e-5, "combine": "latest"},
            ),
            referent="formal",
            non_circularity=(
                "the residual is an independent post-hoc accounting (sum of C_k*M_k "
                "vs the time-integrated root influx minus exports and growth dilution), "
                "computed OUTSIDE the BDF integrator; a conserved total is therefore not "
                "presupposed by the solver. Verified by "
                "tests/test_plant_module.py::test_mass_conservation_source_is_root_uptake."
            ),
        ),
        Hypothesis(
            id="hyp-anion",
            statement="An inside-negative root plasmalemma electrostatically excludes "
                      "the PFAS anion (electrochemical factor e^N ~ 1e2), so passive "
                      "diffusion alone cannot drive net uptake.",
            mode=HypothesisMode.CONFIRMATORY,
            decision_rule=DecisionRule(kind=DecisionRuleKind.QUALITATIVE, expression=R_ANION),
            referent="formal",
            non_circularity=(
                "e^N = exp(z*E_m*F/RT) is a physical consequence of the membrane "
                "potential and valence, fixed by Tier-0 inputs (E_m, z), not fitted to "
                "any uptake observation; root_uptake()'s GHK term uses it directly."
            ),
        ),
        Hypothesis(
            id="hyp-yamazaki",
            statement="Calibration to the Yamazaki (2023) dataset constitutes "
                      "out-of-sample predictive validation of the model's tissue BAFs.",
            mode=HypothesisMode.CONFIRMATORY,
            decision_rule=DecisionRule(kind=DecisionRuleKind.QUALITATIVE, expression=R_YAMA),
            referent="empirical",
        ),
        Hypothesis(
            id="hyp-grain",
            statement="The model predicts brown-rice grain PFAS accumulation accurately "
                      "enough to support dietary risk assessment.",
            mode=HypothesisMode.CONFIRMATORY,
            decision_rule=DecisionRule(kind=DecisionRuleKind.QUALITATIVE, expression=R_GRAIN),
            referent="empirical",
        ),
        Hypothesis(
            id="hyp-soil",
            statement="Per-congener soil sorption (Koc chain-length QSPR) is strongly "
                      "congener-dependent, so the one-way soil coupling cannot be "
                      "replaced by a single constant pore-water concentration.",
            mode=HypothesisMode.CONFIRMATORY,
            decision_rule=DecisionRule(
                kind=DecisionRuleKind.THRESHOLD, expression=R_SOIL,
                params={"statistic": "point", "op": ">", "value": 2.0, "combine": "latest"},
            ),
            referent="formal",
            non_circularity=(
                "the Koc spread is computed from the independently-anchored Higgins-Luthy "
                "+0.55/CF2 QSPR (literature_params.koc), NOT from any plant BAF; the "
                "congener-dependence of soil retardation is a property of the soil "
                "sub-model, not assumed by the plant ODE."
            ),
        ),
        Hypothesis(
            id="hyp-smiles",
            statement="The SMILES structure front-end parameterizes a known PFAS by "
                      "measured read-across, reproducing the curated congener record.",
            mode=HypothesisMode.CONFIRMATORY,
            decision_rule=DecisionRule(kind=DecisionRuleKind.QUALITATIVE, expression=R_SMILES),
            referent="formal",
            non_circularity=(
                "read-across reproducibility is a software property verified by an "
                "independent test (tests/test_pfas_structure.py) that builds the Compound "
                "from a canonical SMILES and checks it equals the curated "
                "params/parameters.json record; it asserts the structure->parameter "
                "mapping is FAITHFUL, not that the parameters are physically correct."
            ),
        ),
        Hypothesis(
            id="hyp-adequacy",
            statement="Driven by the mechanistic ORYZA2000 biomass (oryza_growth) + measured "
                      "transpiration, the model's translocation structure reproduces the "
                      "measured straw (shoot) BAFs across the C4-C12 PFCA/PFSA series under a "
                      "CONSTRAINED (DOF>0) calibration -- per-congener f_xy + a single shared "
                      "L_Ph and kappa_d (13 params / 33 obs, DOF 20), NOT the saturated W2 fit.",
            mode=HypothesisMode.CONFIRMATORY,
            decision_rule=DecisionRule(kind=DecisionRuleKind.QUALITATIVE, expression=R_ADEQ),
            referent="empirical",
        ),
    ]
    target_claims = [
        TargetClaim(id="tc-mass", statement="The transport ODE is mass-conserving.",
                    answers="hyp-mass"),
        TargetClaim(id="tc-anion", statement="The IOC anion-exclusion mechanism holds.",
                    answers="hyp-anion"),
        TargetClaim(id="tc-yamazaki",
                    statement="The model is predictively validated against Yamazaki.",
                    answers="hyp-yamazaki"),
        TargetClaim(id="tc-grain",
                    statement="The model predicts grain accumulation for risk assessment.",
                    answers="hyp-grain"),
        TargetClaim(id="tc-soil",
                    statement="The soil coupling is congener-resolved (not a constant Cwo).",
                    answers="hyp-soil"),
        TargetClaim(id="tc-smiles",
                    statement="SMILES read-across reproduces the curated model.",
                    answers="hyp-smiles"),
        TargetClaim(id="tc-adequacy",
                    statement="The structure reproduces shoot translocation via a constrained fit.",
                    answers="hyp-adequacy"),
    ]
    method = MethodPlan(approaches=[
        "four-compartment dynamic ODE (BDF)",
        "GHK + Michaelis-Menten hybrid root uptake",
        "basis-A binding factor from measured K_PL / K_prot",
        "f_xy xylem loading + L_Ph phloem grain feed",
        "calibration/validation vs Yamazaki 2023, Tang 2026, Kim 2019",
    ])
    return Spec(
        id="pfas-rice", created_at=NOW, version=1, raw_proposal=raw,
        hypotheses=hypotheses, method=method, target_claims=target_claims,
    )


# ============================================================================
# 2) Evidence  (honestly classified: measured | generated | synthetic_proxy).
#    Numbers are the model's ACTUAL outputs (reproduce_demo.py, the test suite,
#    and the documented Tang/Kim validations), keyed to the repo commit.
# ============================================================================
def _ev(id_, kind, data_source, result, bears_on, env):
    return EvidenceItem(
        id=id_, created_at=NOW, spec_id="pfas-rice", kind=kind,
        provenance=Provenance(code_ref=CODE_REF, data_source=data_source, environment=env),
        result=result, bears_on=bears_on,
    )


def evidence(spec, workspace):
    items = []

    # Prior-work record: the literature the model is actually built on (closes the
    # spec-time prior-work reminder honestly; bears on no hypothesis).
    dois = ["10.1016/j.jhazmat.2025.141017",  # Tang 2026
            "10.1016/j.scitotenv.2019.03.240", # Kim 2019
            "10.1029/2019WR025432",            # Brunetti 2019 (IOC DPU)
            "10.1021/acs.est.0c07420"]         # Brunetti 2021
    items.append(_ev(
        "evi-literature", EvidenceKind.LITERATURE, None,
        Result(type="qualitative",
               finding=json.dumps({"acquired": [{"doi": d} for d in dois], "failed": []})),
        [], "literature review (Yamazaki2023, Tang2026, Kim2019, Brunetti2019/2021, "
            "Higgins-Luthy Koc, Chen2025 K_PL, Zhou2025 K_prot)"))

    # H1 mass conservation — GENERATED (a formal property of the ODE), numeric.
    items.append(_ev(
        "evi-mass", EvidenceKind.EXPERIMENT_RUN, "generated",
        Result(type="quantitative", point=1e-6,
               finding="mass-balance residual bounded by rel 1e-6 / abs 1e-9; "
                       "tests/test_plant_module.py::test_mass_conservation_source_is_root_uptake "
                       "PASSES (2 passed)"),
        [Bearing(target_id="hyp-mass", direction=BearingDirection.SUPPORTS)],
        "pytest tests/test_plant_module.py -k mass (numpy/scipy)"))

    # H2 anion exclusion — GENERATED (physical consequence of E_m, z; not fitted).
    items.append(_ev(
        "evi-anion", EvidenceKind.EXPERIMENT_RUN, "generated",
        Result(type="qualitative",
               finding="N = z*E_m*F/RT = +4.67 at E_m=-120 mV, z=-1, T=298.15 K "
                       "=> e^N = 106.8 (test_plant_module asserts e^N ~ 107); the GHK "
                       "term in root_uptake() excludes the anion accordingly."),
        [Bearing(target_id="hyp-anion", direction=BearingDirection.SUPPORTS)],
        "pfas_rice_plant_module_4pool_surf + test_plant_module"))

    # H3 Yamazaki — MEASURED comparison, but the fit is SATURATED (in-sample) ->
    # it REFUTES the *predictive-validation* claim.
    items.append(_ev(
        "evi-yamazaki", EvidenceKind.EXPERIMENT_RUN, "measured",
        Result(type="quantitative", point=0.029,
               finding="reproduce_demo.py: W2 transport fit reproduces Yamazaki "
                       "root/straw/grain BAFs across 11 congeners at log10 RMSE 0.029 "
                       "(e.g. PFOA root 0.49/0.49, straw 0.83/0.83 pred/obs). The fit "
                       "assigns ~3 transport params per congener against 3 tissue "
                       "observations => SATURATED: reproduction is structurally "
                       "guaranteed and is NOT out-of-sample prediction (CLAUDE.md §6)."),
        [Bearing(target_id="hyp-yamazaki", direction=BearingDirection.REFUTES)],
        "reproduce_demo.py (real run; numpy/scipy); data_obs/Yamazaki"))

    # H3 (loop iteration 2): the ACTUAL a-priori prediction. Using the theory/QSPR
    # monotone f_xy (NOT fit to the tissue BAFs; reproduce_demo.py --rec), the genuine
    # out-of-sample predictive error is log10 RMSE 0.837 -- ~29x worse than the saturated
    # 0.029; straw is off 6-40x. The engine REFUTED "0.029 = validation"; the agent then
    # ran the real prediction and recorded its honest error (the loop driving forward).
    items.append(_ev(
        "evi-yamazaki-apriori", EvidenceKind.EXPERIMENT_RUN, "measured",
        Result(type="quantitative", point=0.837,
               finding="reproduce_demo.py --rec: with the theory/QSPR MONOTONE f_xy "
                       "(a-priori, NOT fit to the tissue BAFs) the predictive error is "
                       "log10 RMSE 0.837 vs the saturated W2 fit's 0.029. Straw is off "
                       "6-40x (PFBA 45/11, PFBS 33/2.2). The model does NOT a-priori "
                       "predict the Yamazaki tissue BAFs -- a quantitative confirmation "
                       "of the REFUTED verdict."),
        [Bearing(target_id="hyp-yamazaki", direction=BearingDirection.REFUTES)],
        "reproduce_demo.py --rec (real run; numpy/scipy)"))

    # H3 (loop iteration 3): a bounded MODEL-IMPROVEMENT attempt. The a-priori error
    # is straw-dominated, so the redistributed-shoot model (nstem_leaf) is tried with
    # the SAME monotone f_xy + drivers. Result: a-priori OOS RMSE 0.987 -> 0.951 (a
    # real but MARGINAL gain); short-chain straw improves (PFBA 14.8->10.7 vs 11.0)
    # but the long-chain straw/grain collapse remains (PFDoDA straw 0.35 vs 49.75,
    # ~140x) -- the documented hysteretic-sorption gap. An honest negative: the
    # improvement does not change the REFUTED verdict.
    items.append(_ev(
        "evi-yamazaki-improve", EvidenceKind.EXPERIMENT_RUN, "measured",
        Result(type="quantitative", point=0.951,
               finding="validation/apriori_prediction.py: redistributed-shoot model "
                       "(nstem_leaf, N stem segments + retention) lowers the a-priori "
                       "OOS log10 RMSE 0.987 -> 0.951 (same monotone f_xy + drivers) -- "
                       "MARGINAL. Short-chain straw improves but long-chain straw/grain "
                       "still collapse (PFDoDA straw 0.35 vs 49.75); needs hysteretic "
                       "high-B sorption (docs/nstem_gradient_exploration.md). Predictive "
                       "claim remains refuted."),
        [Bearing(target_id="hyp-yamazaki", direction=BearingDirection.REFUTES)],
        "validation/apriori_prediction.py (real run; model_api simulate vs nstem_leaf)"))

    # H4 grain — MEASURED comparisons; the model is structurally UNDER -> REFUTES.
    items.append(_ev(
        "evi-grain-tang", EvidenceKind.EXPERIMENT_RUN, "measured",
        Result(type="quantitative", point=0.11,
               finding="Tang 2026 OOS: PFOA endosperm TF model 0.11 vs measured 0.95 "
                       "(dw); grain structurally ~3-8x under across congeners; not "
                       "closable by L_Ph / lipid tuning "
                       "(docs/tang2026_grain_units_exploration.md)."),
        [Bearing(target_id="hyp-grain", direction=BearingDirection.REFUTES)],
        "validation/tang2026_nstem_validation.py; raw_si/tang2026_doseresponse.csv"))
    items.append(_ev(
        "evi-grain-kim", EvidenceKind.EXPERIMENT_RUN, "measured",
        Result(type="qualitative",
               finding="Kim 2019: grain BAF reproduced only by FORCING L_Ph "
                       "(0.07 -> 4.43, L_Ph~0.84) as a single-point in-sample anchor; "
                       "Kim is grain-only, so root/straw TF are unconstrained."),
        [Bearing(target_id="hyp-grain", direction=BearingDirection.REFUTES)],
        "literature_params.kim2019_grain_baf(); raw_si/kim2019_*"))

    # H5 soil coupling — GENERATED (a computed property of the Koc QSPR), numeric.
    # Compute the congener Koc spread LIVE from the model's own QSPR.
    import math, sys as _sys
    _sys.path.insert(0, str(HERE.parent / "src"))
    from literature_params import koc
    koc_short = koc(n_perfluoroC=4, head_group="carboxylate")    # PFBA
    koc_long = koc(n_perfluoroC=12, head_group="carboxylate")    # PFDoDA
    spread = math.log10(koc_long / koc_short)
    items.append(_ev(
        "evi-soil", EvidenceKind.EXPERIMENT_RUN, "generated",
        Result(type="quantitative", point=round(spread, 3),
               finding=f"literature_params.koc (Higgins-Luthy +0.55/CF2 QSPR): "
                       f"Koc(PFBA C4)={koc_short:.2f} vs Koc(PFDoDA C12)={koc_long:.0f} "
                       f"L/kg => {spread:.2f} log10 spread (~{koc_long/koc_short:.0f}x). "
                       f"Soil retardation R=1+rho*Kd/theta is therefore strongly "
                       f"congener-dependent; one constant Cwo cannot represent all "
                       f"congeners (CLAUDE.md: short chains leach, long chains buffer)."),
        [Bearing(target_id="hyp-soil", direction=BearingDirection.SUPPORTS)],
        "src/literature_params.koc (live)"))

    # H-adequacy + H4(grain): the CONSTRAINED (DOF 20) structural-adequacy fit.
    # One measured experiment bearing on TWO hypotheses: it SUPPORTS shoot
    # translocation adequacy (straw RMSE 0.048) and REFUTES grain (RMSE 0.987).
    items.append(_ev(
        "evi-adequacy", EvidenceKind.EXPERIMENT_RUN, "measured",
        Result(type="quantitative", point=0.184,
               finding="validation/structural_adequacy_fit.py, driven by the MECHANISTIC "
                       "ORYZA2000 biomass (oryza_growth) + measured Q_TP (forcing_rice) -- "
                       "NOT the logistic placeholder. CONSTRAINED (DOF>0) fits vs Yamazaki:"
                       " A f_xy+global L_Ph+global kappa_d (DOF 20): root 0.45 straw 0.18 "
                       "grain 0.52 overall 0.41; B +per-cong L_Ph (DOF 10): grain 0.36 "
                       "overall 0.35; C +per-cong kappa_d (DOF 10): root 0.26 straw 0.16 "
                       "grain 0.51 overall 0.34. => structure reproduces SHOOT translocation "
                       "(straw ~0.16-0.18, i.e. within ~1.5x) under a constrained fit; root "
                       "needs per-congener kappa_d (->0.26); grain improves with per-congener "
                       "L_Ph (->0.36) but keeps a long-chain residual floor. Whole plant "
                       "within ~factor 2.2 (overall 0.34) at DOF 10 -- vs saturated W2 0.029 "
                       "(DOF 0) and a-priori 0.84."),
        [Bearing(target_id="hyp-adequacy", direction=BearingDirection.SUPPORTS),
         Bearing(target_id="hyp-grain", direction=BearingDirection.REFUTES)],
        "validation/structural_adequacy_fit.py (ORYZA2000 biomass; real run; numpy/scipy)"))

    # H6 SMILES read-across — GENERATED (a software/formal property), qualitative.
    items.append(_ev(
        "evi-smiles", EvidenceKind.EXPERIMENT_RUN, "generated",
        Result(type="qualitative",
               finding="tests/test_pfas_structure.py: 23 passed (RDKit). A canonical "
                       "SMILES that matches a curated congener rebuilds the SAME "
                       "Compound from params/parameters.json (measured read-across). "
                       "CAVEAT: for NOVEL structures f_xy is provisional (QSPR/"
                       "interpolated), NOT validated — this claim is scoped to known "
                       "structures only."),
        [Bearing(target_id="hyp-smiles", direction=BearingDirection.SUPPORTS)],
        "pytest tests/test_pfas_structure.py (RDKit); src/pfas_structure.py"))

    return items


# ============================================================================
# 3) Verdict trails  (the in-session agent's chief-over-N judgments; the engine
#    refuses a binding qualitative verdict without a well-formed trail — F2 gate).
# ============================================================================
def _trail(hyp_id, rubric, direction, level, chief_basis, panel):
    return VerdictTrail(
        hypothesis_id=hyp_id, rule_kind="qualitative", rubric_expression=rubric,
        panel=panel,
        chief=ChiefVerdict(direction=direction, level=level, basis=chief_basis),
        provenance=VerdictProvenance(spec_version=1, timestamp=NOW.isoformat(),
                                     agent_ids=["claude-in-session"]),
    )


def verdicts():
    SUP, REF = BearingDirection.SUPPORTS, BearingDirection.REFUTES
    STRONG, MOD = ConfidenceLevel.STRONG, ConfidenceLevel.MODERATE
    return {
        "hyp-anion": _trail(
            "hyp-anion", R_ANION, SUP, STRONG,
            "Both panelists agree the e^N~1e2 exclusion is a closed-form consequence of "
            "Tier-0 inputs (E_m, z), independently asserted by test_plant_module; it is a "
            "structural/formal property, not a fitted one, so the IOC signature holds.",
            [PanelVerdict(direction=SUP, level=STRONG,
                          basis="e^N=106.8 follows analytically from N=zE_mF/RT; matches "
                                "the GHK term in root_uptake()."),
             PanelVerdict(direction=SUP, level=STRONG,
                          basis="the anti-exclusion carrier (Vmax/Km) is required to "
                                "overcome e^N, confirming passive diffusion alone is "
                                "insufficient — exactly the rule's criterion.")]),
        "hyp-yamazaki": _trail(
            "hyp-yamazaki", R_YAMA, REF, STRONG,
            "Decisive panel reasoning: the W2 transport parameters are FIT to the same "
            "Yamazaki tissues they then 'predict' (~3 params / 3 obs per congener), so "
            "RMSE 0.029 is saturated in-sample reproduction. The a-priori test "
            "(monotone f_xy, --rec) gives the REAL predictive error log10 RMSE 0.837 "
            "(straw off 6-40x) -- the model does not predict out of sample. The "
            "predictive-validation hypothesis is refuted quantitatively.",
            [PanelVerdict(direction=REF, level=STRONG,
                          basis="pred ~ obs to 2 decimals (0.49/0.49) is the signature of "
                                "a saturated fit, not prediction."),
             PanelVerdict(direction=REF, level=MOD,
                          basis="CLAUDE.md itself flags RMSE 0.029 as 'a saturated W2 fit "
                                "... reproduction is guaranteed, NOT predictive "
                                "validation' and notes the ordering is congener-dependent.")]),
        "hyp-grain": _trail(
            "hyp-grain", R_GRAIN, REF, STRONG,
            "Decisive panel reasoning: against MEASURED data the model under-predicts "
            "grain by ~3-8x (Tang PFOA endosperm 0.11 vs 0.95) and only matches Kim by "
            "forcing a single-point L_Ph anchor. A structural under-prediction of the "
            "risk-relevant compartment cannot support dietary risk assessment => refute.",
            [PanelVerdict(direction=REF, level=STRONG,
                          basis="Tang OOS shows a structural grain floor not closable by "
                                "L_Ph/lipid (docs/tang2026_grain_units_exploration.md)."),
             PanelVerdict(direction=REF, level=STRONG,
                          basis="the only 'match' (Kim) is an in-sample forced anchor on "
                                "grain alone; no independent grain prediction exists.")]),
        "hyp-smiles": _trail(
            "hyp-smiles", R_SMILES, SUP, STRONG,
            "Both panelists agree this is a faithful structure->parameter mapping, an "
            "independently tested software property (23/23). It is SUPPORTED only for "
            "KNOWN structures (read-across); novel-structure prediction is out of scope "
            "and remains provisional.",
            [PanelVerdict(direction=SUP, level=STRONG,
                          basis="tests/test_pfas_structure.py passes 23/23: a SMILES-built "
                                "known PFAS equals the curated record."),
             PanelVerdict(direction=SUP, level=MOD,
                          basis="the rule is scoped to read-across reproducibility, which "
                                "is deterministic and verified; it does not over-claim "
                                "novel-structure prediction.")]),
        "hyp-adequacy": _trail(
            "hyp-adequacy", R_ADEQ, SUP, STRONG,
            "Decisive panel reasoning: driven by the user's MECHANISTIC ORYZA2000 biomass "
            "(not the placeholder), the CONSTRAINED DOF-20 fit reaches straw RMSE ~0.18 "
            "(0.16 at DOF 10) across the full C4-C12 PFCA/PFSA range with per-congener f_xy "
            "+ shared L_Ph/kappa_d. So the translocation structure (GHK exclusion + f_xy "
            "TSCF + binding) reproduces shoot accumulation under a non-saturated calibration. "
            "The whole plant reaches overall ~0.34 (within ~factor 2.2) at DOF 10; grain "
            "keeps a residual long-chain floor (separately refuted).",
            [PanelVerdict(direction=SUP, level=STRONG,
                          basis="straw log10 RMSE ~0.16-0.18 at DOF 20/10 on realistic "
                                "ORYZA2000 biomass is a genuine (non-saturated) goodness-of-"
                                "fit -- structural adequacy, not expressiveness."),
             PanelVerdict(direction=SUP, level=MOD,
                          basis="root improves to 0.26 with per-congener kappa_d and the "
                                "whole plant to overall 0.34 at DOF 10; the shoot "
                                "translocation claim holds, scoped away from the grain floor.")]),
    }


# ============================================================================
# Drivers
# ============================================================================
def _write_verdicts(run_dir: Path, trails: dict) -> None:
    vdir = run_dir / "verdicts"
    vdir.mkdir(parents=True, exist_ok=True)
    for hyp_id, trail in trails.items():
        (vdir / f"{hyp_id}.json").write_text(
            json.dumps(trail.model_dump(mode="json"), indent=2, ensure_ascii=False),
            encoding="utf-8")


def run_main():
    spec = build_spec()
    run_dir = HERE / "runs" / spec.id
    run_dir.mkdir(parents=True, exist_ok=True)
    _write_verdicts(run_dir, verdicts())        # author verdicts BEFORE the loop
    result = run_checkpoint_loop(run_dir=run_dir, spec=spec, experiment=evidence,
                                 workspace_dir=HERE)
    print(f"\n=== MAIN RUN  '{spec.id}'  ({result.iterations} iteration(s)) ===")
    for c in result.claims:
        print(f"  {c.answers:14s} -> {c.status.value.upper():9s} | {c.confidence.basis[:90]}")
    if result.unresolved:
        print(f"  unresolved: {result.unresolved}")

    report = verify_run(run_dir)
    print(f"\n=== VERIFY  (headless re-derivation; record digest sha256:{report.digest[:16]}…) ===")
    for o in report.outcomes:
        rd = o.rederived_status.value if o.rederived_status else "n/a"
        print(f"  {o.hypothesis_id:14s} -> {o.result:11s} (recorded={o.recorded_status.value}, "
              f"re-derived={rd})")
    print(f"  all reproduced: {report.all_reproduced}")


def run_trap():
    """The 'rice-failure' reproduction: feed the bundled demo BAFs (synthetic_proxy)
    to an empirical hypothesis and watch sci-adk's validity gate REFUSE to certify it."""
    raw = _raw_proposal()
    hyp = Hypothesis(
        id="hyp-demo",
        statement="The model's bundled demonstration BAFs establish that it predicts "
                  "field PFAS bioaccumulation in rice.",
        mode=HypothesisMode.EXPLORATORY,
        decision_rule=DecisionRule(kind=DecisionRuleKind.QUALITATIVE, expression=R_DEMO),
        referent="empirical",
    )
    spec = Spec(id="pfas-rice-trap", created_at=NOW, version=1, raw_proposal=raw,
                hypotheses=[hyp], method=MethodPlan(approaches=["bundled _demo()"]),
                target_claims=[TargetClaim(id="tc-demo",
                                           statement="The demo proves field prediction.",
                                           answers="hyp-demo")])
    run_dir = HERE / "runs" / spec.id
    run_dir.mkdir(parents=True, exist_ok=True)
    _write_verdicts(run_dir, {
        "hyp-demo": _trail(
            "hyp-demo", R_DEMO, BearingDirection.SUPPORTS, ConfidenceLevel.STRONG,
            "Naive reading: the demo prints root>straw>grain BAFs, so the model 'predicts' "
            "bioaccumulation.",
            [PanelVerdict(direction=BearingDirection.SUPPORTS, level=ConfidenceLevel.STRONG,
                          basis="the demo runs and prints plausible BAFs.")]),
    })

    def demo_evidence(_spec, _ws):
        return [_ev(
            "evi-demo", EvidenceKind.EXPERIMENT_RUN, "synthetic_proxy",
            Result(type="qualitative",
                   finding="bundled _demo() BAFs (root 24.x / straw / grain) from "
                           "pfas_rice_plant_module_4pool_surf with ILLUSTRATIVE, NOT "
                           "calibrated parameters (Cwo=1 ug/L constant; PFOA-like "
                           "placeholder K_*/f_xy). Comment in source: 'values are "
                           "illustrative, NOT calibrated.'"),
            [Bearing(target_id="hyp-demo", direction=BearingDirection.SUPPORTS)],
            "pfas_rice_plant_module_4pool_surf._demo()")]

    print(f"\n=== TRAP RUN  '{spec.id}'  (the rice-failure reproduction) ===")
    try:
        run_checkpoint_loop(run_dir=run_dir, spec=spec, experiment=demo_evidence,
                            workspace_dir=HERE)
        print("  !! NO HALT — the gate failed to catch the synthetic_proxy evidence")
    except ValidityHalt as e:
        print(f"  HALT on '{e.hypothesis_id}':")
        print(f"  {e.reason}")
        (run_dir / "VALIDITY_HALT.txt").write_text(
            f"evidence-validity halt on '{e.hypothesis_id}':\n{e.reason}\n", encoding="utf-8")
        print("  -> recorded to VALIDITY_HALT.txt (no Claim was written; the ungrounded "
              "empirical result is REFUSED, not self-certified)")


if __name__ == "__main__":
    run_main()
    run_trap()
