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


def test_synthetic_proxy_halt_recorded():
    """The demo-BAF (synthetic_proxy) refusal artifact must persist (rice-failure block)."""
    halt = TRAP / "VALIDITY_HALT.txt"
    assert halt.exists(), "the synthetic_proxy validity-HALT artifact is missing"
    assert "synthetic_proxy" in halt.read_text(encoding="utf-8")
