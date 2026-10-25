# Phase 106 run log

- 2026-10-23 00:53  DESIGN.md frozen, SHA-256 f074a8e3b0b587b1f63b19d4af5b6a28afe1a5af134603752480a62730042cc7. No experiment code exists yet.
- 2026-10-23 00:55  p106.py written; `validate` run: V1 (2,500 graphs, 74,400 regions, 0 mismatches against p101_core), V2-V5 all pass.
- 2026-10-23 00:58  `enumerate` run: 163 s; distinct two-logical profiles n=2..6: 2, 7, 29, 127, 735; classical n=2..6: 10, 26, 63, 131, 258.
- 2026-10-23 01:06  `analyse` run. Primary cell (n=6, p=0.2): gain 0 at every w; P1 nominal (tie at w=0.5 only), P2 fail, P3 fail (as predicted), P4 holds (constant 0), P5 fail.
- 2026-10-23 01:06  **Deviation 1 (bug fix, classical twin).** ALL_cl contained codes that never store bit B (e.g. three copies of A), which SEP_cl and every quantum code exclude (quantum validity S(R)=2 means all qubits recover both logicals). Fixed in `load_profiles`: a classical code must determine both bits from all coordinates. Unfixed numbers kept in `results_classical_unrestricted.json` / `analyse_output_classical_unrestricted.txt` (D3: both reported). No primary verdict changes; the fix only removes spurious classical gains at lopsided w.
- 2026-10-23 01:06  Post-hoc check `p106_random_check.py`: 300,000 random labelled 8-vertex graphs, every (p, w): none beats the exhaustive n=6 optimum (largest excess 0).
- 2026-10-23 01:06  Note: the cost-per-qubit table in `results.json` labels families 'shared'/'separate'; an early printout abbreviated both to 's' (display only, no effect on results).
- 2026-10-23 01:06  D4 exploratory n=7 started (`p106_explore7.py`, search for ALL(7), exact SEP(7)).
- 2026-10-23 01:09  D4 exploratory n=7 done (3 min): gain p=0.1: +0.031 (w=0.5) to 0 (w>=0.9); p=0.2: +0.0006 at w=0.5, 0 after; p=0.3: 0. Written to explore7.json.
- 2026-10-23 01:09  RESULTS.md written. No file outside phase106/ edited; nothing committed. Python caches in phase106/__pycache__ only.
- 2026-10-23 17:53 Folder moved to qec/notebooks/goal_codes/phase106 (from llm-project/Colab). No content changed.
- 2026-10-23 19:59 Parent folder renamed: qec/notebooks/goal_codes -> qec/notebooks/what_to_keep. Earlier lines here keep the old name.
- 2026-10-23 20:02 Rename reverted: the folder is qec/notebooks/goal_codes again.
