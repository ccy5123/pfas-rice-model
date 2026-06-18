"""
Always-on over-claim guard: re-derive the committed sci-adk rigor review and
assert the model's claims stay honestly classified.

This institutionalizes the sci-adk loop (https://github.com/ccy5123/sci-adk) as a
repo gate. It re-runs sci-adk's headless `verify` over the COMMITTED run
(`sci_adk_review/runs/pfas-rice`) -- no LLM, no model re-run -- and fails if:
  * any recorded claim no longer reproduces from the record, OR
  * an EMPIRICAL predictive claim (Yamazaki predictive-validation / grain risk
    assessment) is ever marked SUPPORTED -- the exact over-claim the review found
    to be unjustified (saturated in-sample fit; a-priori RMSE ~0.84-0.95), and the
    "rice-failure defect" sci-adk exists to prevent, OR
  * the synthetic_proxy validity-HALT artifact (the demo-BAF refusal) disappears.

Skips cleanly when sci-adk is not installed (like the RDKit / HYDRUS tests), so it
never blocks the suite -- but when sci-adk IS present it is a hard gate.

Regenerate the run after a model change with:  python sci_adk_review/build_review.py
"""
import json
from pathlib import Path

import pytest

sci_adk = pytest.importorskip("sci_adk", reason="sci-adk not installed (pip install -e <sci-adk>)")

from sci_adk.loop.verify import verify_run  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
RUN = ROOT / "sci_adk_review" / "runs" / "pfas-rice"
RUN_LC = ROOT / "sci_adk_review" / "runs" / "pfas-rice-longchain"
RUN_CARRIER = ROOT / "sci_adk_review" / "runs" / "pfas-rice-carrier"
RUN_OOS = ROOT / "sci_adk_review" / "runs" / "pfas-rice-oos-tang"
RUN_OOS_LIPID = ROOT / "sci_adk_review" / "runs" / "pfas-rice-oos-lipid"
RUN_OOS_MULTI = ROOT / "sci_adk_review" / "runs" / "pfas-rice-oos-multidataset"
RUN_LC_COMPLETE = ROOT / "sci_adk_review" / "runs" / "pfas-rice-longchain-complete"
RUN_LC_DECOUPLE = ROOT / "sci_adk_review" / "runs" / "pfas-rice-longchain-decouple"
RUN_CONSOLIDATION = ROOT / "sci_adk_review" / "runs" / "pfas-rice-consolidation"
RUN_SELECTION = ROOT / "sci_adk_review" / "runs" / "pfas-rice-model-selection"
TRAP = ROOT / "sci_adk_review" / "runs" / "pfas-rice-trap"

# Empirical predictive claims that MUST NOT be SUPPORTED (the over-claim guard).
_PREDICTIVE = {"hyp-yamazaki", "hyp-grain"}
# Formal/computational claims that SHOULD hold (structure/mechanism/tooling).
_FORMAL_SUPPORTED = {"hyp-mass", "hyp-anion", "hyp-soil", "hyp-smiles"}


def _claim_status():
    out = {}
    for p in (RUN / "claims").glob("claim-*.json"):
        d = json.loads(p.read_text(encoding="utf-8"))
        out[d["answers"]] = d["status"]
    return out


@pytest.mark.skipif(not RUN.exists(), reason="committed sci-adk run absent")
def test_record_reproduces():
    """sci-adk verify: every recorded claim re-derives from the record (exit-0 path)."""
    report = verify_run(RUN)
    assert report.all_reproduced, (
        "sci-adk verify: a recorded claim DIVERGED/UNRESOLVED -- regenerate with "
        "`python sci_adk_review/build_review.py` and review:\n"
        + "\n".join(f"  {o.hypothesis_id}: {o.result}" for o in report.outcomes))


@pytest.mark.skipif(not RUN.exists(), reason="committed sci-adk run absent")
def test_empirical_predictive_claims_not_supported():
    """The over-claim guard: predictive claims must stay REFUTED/PROPOSED, never SUPPORTED."""
    status = _claim_status()
    for hyp in _PREDICTIVE:
        assert status.get(hyp) != "supported", (
            f"{hyp} is SUPPORTED -- an empirical predictive claim was certified. The "
            f"Yamazaki fit is saturated in-sample (a-priori RMSE ~0.84-0.95) and grain is "
            f"structurally under-predicted; certifying these is the rice-failure defect. "
            f"Fix the model/evidence, do not relax this gate.")


@pytest.mark.skipif(not RUN.exists(), reason="committed sci-adk run absent")
def test_formal_claims_supported():
    """Sanity: the structural/mechanistic/tooling claims remain SUPPORTED."""
    status = _claim_status()
    for hyp in _FORMAL_SUPPORTED:
        assert status.get(hyp) == "supported", f"{hyp} regressed to {status.get(hyp)!r}"


@pytest.mark.skipif(not RUN_LC.exists(), reason="long-chain sub-investigation absent")
def test_longchain_run_reproduces():
    """The long-chain mechanism sub-investigation re-derives from its record (LC1/LC2
    SUPPORTED, LC3 REFUTED)."""
    report = verify_run(RUN_LC)
    assert report.all_reproduced, "\n".join(f"  {o.hypothesis_id}: {o.result}" for o in report.outcomes)


@pytest.mark.skipif(not RUN_CARRIER.exists(), reason="carrier-QSPR sub-investigation absent")
def test_carrier_run_reproduces():
    """The carrier-QSPR sub-investigation (compiled via the `sci-adk run` CLI) re-derives
    from its record (hyp-001 REFUTED -- the long-chain carrier enhancement is not QSPR-able)."""
    report = verify_run(RUN_CARRIER)
    assert report.all_reproduced, "\n".join(f"  {o.hypothesis_id}: {o.result}" for o in report.outcomes)


@pytest.mark.skipif(not RUN_OOS.exists(), reason="OOS Tang cross-dataset run absent")
def test_oos_tang_run_reproduces():
    """The out-of-sample cross-dataset test (compiled via the `sci-adk run` CLI) re-derives
    from its record (hyp-001 REFUTED -- theory params do NOT predict the independent Tang
    dataset out-of-sample: OOS RMSE 1.23 vs in-sample 0.52). This is the project's central
    predictive-validation result on data NOT used to fit, and must stay REFUTED."""
    report = verify_run(RUN_OOS)
    assert report.all_reproduced, "\n".join(f"  {o.hypothesis_id}: {o.result}" for o in report.outcomes)


@pytest.mark.skipif(not RUN_OOS_LIPID.exists(), reason="OOS lipid-generalization run absent")
def test_oos_lipid_run_reproduces():
    """The lipid-mechanism out-of-sample GENERALIZATION test (compiled via the `sci-adk run`
    CLI) re-derives from its record (hyp-001 SUPPORTED -- the K_PL-gated lipid loading fit on
    YAMAZAKI, NOT Tang, drops the independent-dataset Tang OOS RMSE 1.23 -> 0.52, matching the
    in-sample refit; the project's first strong cross-dataset OOS predictive success). This is
    the mechanism generalizing, NOT in-sample reproduction -- it is a genuine predictive
    success and may legitimately be SUPPORTED (cf. the hyp-yamazaki/grain over-claim guard)."""
    report = verify_run(RUN_OOS_LIPID)
    assert report.all_reproduced, "\n".join(f"  {o.hypothesis_id}: {o.result}" for o in report.outcomes)


@pytest.mark.skipif(not RUN_OOS_MULTI.exists(), reason="multi-dataset OOS robustness run absent")
def test_oos_multidataset_run_reproduces():
    """The multi-dataset OOS robustness test (compiled via the `sci-adk run` CLI) re-derives
    from its record (hyp-001 SUPPORTED -- the lipid mechanism's OOS generalization holds across
    BOTH clean independent datasets: Tang per-organ TF 0.52<1.23 and Kim grain 0.48<2.05, not a
    Tang artifact; Li 2025 confounded/inconclusive as pre-registered). Strengthens the n=3 Tang
    result to a robust multi-dataset cross-validation."""
    report = verify_run(RUN_OOS_MULTI)
    assert report.all_reproduced, "\n".join(f"  {o.hypothesis_id}: {o.result}" for o in report.outcomes)


@pytest.mark.skipif(not RUN_LC_COMPLETE.exists(), reason="long-chain complete run absent")
def test_longchain_complete_run_reproduces():
    """The long-chain COMPLETE-resolution test (build_longchain_complete.py: combines the 3
    proposed levers -- 2-pool + lipid + LC6 root-matching carrier -- into ONE model) re-derives
    from its record: 4/4 SUPPORTED -- root closes (RMSE 0.002), grain closes (0.23), but the
    SHOOT does NOT close (straw RMSE 0.39) because the root-fixing carrier over-feeds the shoot
    (PFDoDA straw 2.27x). FINDINGS sec.7's 'complete resolution' is NOT a simultaneous closure;
    the long chains need a root->shoot decoupling. This is a true factual finding, not a
    predictive over-claim."""
    report = verify_run(RUN_LC_COMPLETE)
    assert report.all_reproduced, "\n".join(f"  {o.hypothesis_id}: {o.result}" for o in report.outcomes)


@pytest.mark.skipif(not RUN_LC_DECOUPLE.exists(), reason="long-chain decouple run absent")
def test_longchain_decouple_run_reproduces():
    """The root->shoot DECOUPLING test (build_longchain_decouple.py: an irreversible bound-store
    sequestration lever seq) re-derives from its record: 3/3 SUPPORTED -- the lever INFLATES the
    root (PFDoDA 6.94x) instead of relieving the shoot, does NOT achieve clean within-factor-2
    simultaneous closure (best gap 0.336), and improves on the complete recipe only marginally
    (0.020 log10). Honest near-negative: asymmetric bound-store kinetics are the wrong lever; the
    fix must break the uptake<->mobile-conc coupling. A true factual finding, not an over-claim."""
    report = verify_run(RUN_LC_DECOUPLE)
    assert report.all_reproduced, "\n".join(f"  {o.hypothesis_id}: {o.result}" for o in report.outcomes)


@pytest.mark.skipif(not RUN_CONSOLIDATION.exists(), reason="consolidation run absent")
def test_consolidation_run_reproduces():
    """The engine-rendered consolidation (built by `build_consolidation.py`: a frozen
    threshold-hypothesis Spec whose statistics ARE the verified sub-run outputs; the paper is
    rendered by sci-adk, the narrative supplied as LaTeX-safe prose INPUT, not hand-authored)
    re-derives from its record: 5/5 SUPPORTED -- reproducibility, naive-OOS-fails,
    lipid-generalizes, lipid-robust, structural-adequacy. The synthesis claims are restatements
    of the sub-run records, so they may legitimately be SUPPORTED."""
    report = verify_run(RUN_CONSOLIDATION)
    assert report.all_reproduced, "\n".join(f"  {o.hypothesis_id}: {o.result}" for o in report.outcomes)


@pytest.mark.skipif(not RUN_SELECTION.exists(), reason="model-selection run absent")
def test_model_selection_run_reproduces():
    """The transport model-selection verdict (build_model_selection.py) re-derives from its
    record: 4/4 SUPPORTED -- the K_PL-gated lipid mechanism is the CONSISTENT best across every
    measured dataset (in-sample Yamazaki margin 0.65, Tang OOS 0.72, Kim OOS 1.57; min 0.65),
    so it is the recommended transport configuration (kept opt-in pending a reliable 2-pool
    root). A selection over already-adjudicated results; the wins are genuine, so SUPPORTED."""
    report = verify_run(RUN_SELECTION)
    assert report.all_reproduced, "\n".join(f"  {o.hypothesis_id}: {o.result}" for o in report.outcomes)


def test_synthetic_proxy_halt_recorded():
    """The demo-BAF (synthetic_proxy) refusal artifact must persist (rice-failure block)."""
    halt = TRAP / "VALIDITY_HALT.txt"
    assert halt.exists(), "the synthetic_proxy validity-HALT artifact is missing"
    assert "synthetic_proxy" in halt.read_text(encoding="utf-8")
