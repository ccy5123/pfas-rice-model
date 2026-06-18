# sci-adk rigor audit — index

Adversarial rigor audit of this repository's four-compartment PFAS rice
uptake model, run with **sci-adk** (https://github.com/ccy5123/sci-adk):
*agents propose; the engine judges by frozen criteria; no
self-certification.*

## Read this first

- **`runs/pfas-rice-consolidation/paper/draft.tex`** — single
  consolidated, citable synthesis paper (EN), **rendered by sci-adk**.
  Built by `build_consolidation.py`, which freezes a threshold-hypothesis
  Spec whose statistics ARE the verified sub-run outputs; the engine
  resolves the numeric claims and renders the paper, with the narrative
  supplied as LaTeX-safe **prose input** (not hand-authored, not
  LLM-generated). All runs, a master ledger (in the discussion),
  centralized caveats, verified digests. Provenance is honest by
  construction: `\author{sci-adk (deterministic render)}`; every figure
  traces to a verified record; 5/5 claims reproduce under `sci-adk
  verify`.
- **`FINDINGS.md`** — the authoritative Korean narrative.

## The runs (all `sci-adk verify` → exit 0)

| run | scope | headline verdict |
|---|---|---|
| `runs/pfas-rice` | main audit: H1–H7 | formal SUPPORTED; predictive (H3/H4) REFUTED; structural-adequacy (H7) SUPPORTED |
| `runs/pfas-rice-trap` | synthetic-data trap | **HALT** — `synthetic_proxy` on an empirical hypothesis, no claim written |
| `runs/pfas-rice-longchain` | LC1–LC5b mechanism | lipid-facilitated loading is the right direction; PFDoDA (C12) needs an active carrier |
| `runs/pfas-rice-carrier` | LC6 carrier QSPR | REFUTED — carrier enhancement not QSPR-able from chain length |
| `runs/pfas-rice-oos-tang` | cross-dataset OOS (free-anion) | REFUTED — OOS 1.23 vs in-sample 0.52 |
| `runs/pfas-rice-oos-lipid` | OOS with lipid mechanism | SUPPORTED — OOS 1.23 → 0.52 (mechanism generalizes) |
| `runs/pfas-rice-oos-multidataset` | robustness (Tang + Kim + Li) | SUPPORTED — robust across two clean independent datasets |
| `runs/pfas-rice-longchain-complete` | the 3-lever "complete resolution" as one model | 4/4 SUPPORTED — root+grain close but shoot does NOT (carrier over-feeds): not a simultaneous closure |
| `runs/pfas-rice-consolidation` | engine-rendered synthesis paper | 5/5 SUPPORTED — reproducibility, naive-OOS-fails, lipid-generalizes, lipid-robust, structural-adequacy |

## Specs / drivers

- `proposal.md` — frozen four-pane pre-registration for the main run.
- `proposal_*.md` — pre-registrations (incl. `proposal_consolidation.md`,
  `proposal_longchain_complete.md`).
- `build_review.py`, `build_longchain.py`, `build_longchain_complete.py`,
  `build_consolidation.py` — reproducible run drivers (the engine
  compiles/judges/renders); `validation/longchain_complete.py` is the
  live experiment behind the complete-resolution run.
- `run_rigor.sh` — local full regenerate + verify.

## Reproduce

```bash
pip install -e /path/to/sci-adk          # or PYTHONPATH=sci-adk/src
pip install numpy scipy pytest rdkit     # for the model-output evidence
python sci_adk_review/build_review.py             # main run
python sci_adk_review/build_longchain_complete.py # long-chain "complete resolution" test
python sci_adk_review/build_consolidation.py      # engine-rendered synthesis paper
for r in pfas-rice pfas-rice-longchain pfas-rice-carrier \
         pfas-rice-oos-tang pfas-rice-oos-lipid pfas-rice-oos-multidataset \
         pfas-rice-longchain-complete pfas-rice-consolidation; do
  sci-adk verify sci_adk_review/runs/$r   # exit 0, all claims REPRODUCED
done
```

Standing guards in `tests/test_sci_adk_rigor.py` re-verify the committed
runs on every push and fail if an empirical predictive claim is ever
promoted to SUPPORTED without basis.
