#!/usr/bin/env bash
# Re-run the sci-adk rigor review end-to-end and re-verify (the "always-on" loop,
# run locally). Regenerates runs/ from proposal.md + the live model outputs, then
# headlessly re-derives belief and runs the over-claim guard.
#
# Requires: sci-adk installed (pip install -e <sci-adk>), plus numpy scipy pytest
# (and rdkit for the SMILES evidence). Run from the repo root:
#   bash sci_adk_review/run_rigor.sh
set -euo pipefail
cd "$(dirname "$0")/.."

echo "== 1. regenerate the rigor review (compile + evidence + verdicts + loop) =="
python sci_adk_review/build_review.py

echo "== 2. headless re-verification (no LLM; exit 0 iff all claims reproduce) =="
sci-adk verify sci_adk_review/runs/pfas-rice

echo "== 3. over-claim guard (predictive claims must not be SUPPORTED) =="
pytest tests/test_sci_adk_rigor.py -q

echo "== rigor loop OK =="
