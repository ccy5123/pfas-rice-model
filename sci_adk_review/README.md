# sci-adk rigor audit — index

Adversarial rigor audit of this repository's four-compartment PFAS rice
uptake model, run with **sci-adk** (https://github.com/ccy5123/sci-adk):
*agents propose; the engine judges by frozen criteria; no
self-certification.*

## Read this first

- **`docs/sci_adk_rigor_review.tex`** — single consolidated, citable
  manuscript (EN). All seven runs, one master ledger, the narrative
  arc, centralized caveats, verified digests. **Provenance:** this is an
  *agent-authored cross-run synthesis* (hand-written from the records),
  **not** a sci-adk render — sci-adk is the referee/renderer, not the
  author. The engine's own paper output is the per-run, deterministically
  rendered `runs/*/paper/draft.tex` (frozen-record artifacts; the
  consolidation supersedes them only as readable prose, since their
  non-ASCII bodies did not render). Every figure is traceable to a
  verified record.
- **`FINDINGS.md`** — the authoritative Korean narrative (also
  agent-authored from the records).

## The seven runs (all `sci-adk verify` → exit 0)

| run | scope | headline verdict |
|---|---|---|
| `runs/pfas-rice` | main audit: H1–H7 | formal SUPPORTED; predictive (H3/H4) REFUTED; structural-adequacy (H7) SUPPORTED |
| `runs/pfas-rice-trap` | synthetic-data trap | **HALT** — `synthetic_proxy` on an empirical hypothesis, no claim written |
| `runs/pfas-rice-longchain` | LC1–LC5b mechanism | lipid-facilitated loading is the right direction; PFDoDA (C12) needs an active carrier |
| `runs/pfas-rice-carrier` | LC6 carrier QSPR | REFUTED — carrier enhancement not QSPR-able from chain length |
| `runs/pfas-rice-oos-tang` | cross-dataset OOS (free-anion) | REFUTED — OOS 1.23 vs in-sample 0.52 |
| `runs/pfas-rice-oos-lipid` | OOS with lipid mechanism | SUPPORTED — OOS 1.23 → 0.52 (mechanism generalizes) |
| `runs/pfas-rice-oos-multidataset` | robustness (Tang + Kim + Li) | SUPPORTED — robust across two clean independent datasets |

## Specs / drivers

- `proposal.md` — frozen four-pane pre-registration for the main run.
- `proposal_*.md` — pre-registrations for the sub-investigation runs.
- `build_review.py`, `build_longchain.py` — reproducible run drivers.
- `run_rigor.sh` — local full regenerate + verify.

## Reproduce

```bash
pip install -e /path/to/sci-adk          # or PYTHONPATH=sci-adk/src
pip install numpy scipy pytest rdkit     # for the model-output evidence
python sci_adk_review/build_review.py
for r in pfas-rice pfas-rice-longchain pfas-rice-carrier \
         pfas-rice-oos-tang pfas-rice-oos-lipid pfas-rice-oos-multidataset; do
  sci-adk verify sci_adk_review/runs/$r   # exit 0, all claims REPRODUCED
done
```

Standing guards in `tests/test_sci_adk_rigor.py` re-verify the committed
runs on every push and fail if an empirical predictive claim is ever
promoted to SUPPORTED without basis.
