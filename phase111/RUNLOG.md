# Phase 111 run log

Times are local (Barcelona, UTC+1). Entries are added in order and never edited; corrections are new entries.

## 2026-10-25: live build and freeze

- 16:04. First attempt at `p111_build.py --live` stopped before connecting: the IBM instance was entered without its
  `crn:` prefix and the service rejected it. Nothing was written. The power script, pasted with it, then ran on
  the existing development build; its output was removed by the next build.
- 16:07. `p111_build.py --live` against ibm_fez, calibration of 2026-10-25 14:45. All six pairs pass V1–V3, all
  baselines pass their exact check, and the re-check of the written files by `p111_verify.py` passes (60 circuits;
  20 archived Test B circuits identified). Same layouts, decoders, output positions and two-qubit gate counts as the
  development build; the compiled preparation of the both region differs in single-qubit gates (six circuits).
- 16:46. `p111_power.py` on the saved compilation target: expected z of D1–D3 at 4,000 shots 68.8, 98.6, 46.5;
  shot count 4,000.
- Clean-copy check: a copy of this repository's files alone (no other environment state) reproduced the
  verification report byte for byte; the analysis self-test, the broken-circuit test and the offline notebook test
  pass; `reproduce.py` reproduces every number in the paper.
- `FREEZE_MANIFEST.txt` written: SHA-256 `f1b99e67ee15ed148b2550ea78d37c5094ff1fbdea9a4e8da1a8ec45532260d6`.
  The files it lists are stored byte for byte (`.gitattributes`: no line-ending conversion), so their hashes hold on
  any checkout. A private working copy of this phase exists outside this repository; its line endings may differ,
  so only the files here are authoritative.
- Deadline: the circuits are run within 72 hours of compilation, by 2026-10-28 16:07.
- 18:10. `FREEZE_MANIFEST.txt` stamped with OpenTimestamps (`FREEZE_MANIFEST.txt.ots`, four calendars, pending
  Bitcoin confirmation), after the freeze commit 3eacc77. The confirmed proof is checked before the run and the result
  recorded in `TIMESTAMP_CHECK.json`.
- 19:56. Proof upgraded: Bitcoin block 973103 (2026-10-25 19:00:31, 18:00:31 UTC); its merkle root equals the
  proof's commitment on two block explorers. Recorded in `TIMESTAMP_CHECK.json`. No QPU job had been submitted.
- 20:11–20:14. Run on ibm_fez from the frozen files (all 16 hashed files matched; the notebook itself differed only by
  the instance entered and its outputs). Jobs, each tagged with the experiment ID p111-f1b99e67ee15 and each
  retrieved with its submitted circuits matching the frozen fingerprints: A2 dbf58r8m94ac738ssbi0 (19:11:41 UTC),
  pilot dbf59bnsstms73fd4ntg (19:12:46 UTC), B dbf59tvpc3tc739mp16g (19:13:59 UTC). No other job on the instance that
  day; nothing resubmitted. QPU usage 23, 6 and 65 seconds.
- Part A2 replicated: E_B(main) 1.383 (0.016), E_A(main) 0.828 (0.018), E_A(both) 1.558 (0.014); B1, B2, B3 hold by
  23.8, 40.4 and 22.8 standard errors. E_B(main) is lower than in the first run (1.521).
- Pilot: F 0.760 (main/B) and 0.694 (both/A); no technical failure.
- Part B passed: F_B(main) 0.761 (0.005), F_A(main) 0.431 (0.006), F_A(both) 0.683 (0.006); D1, D2, D3 hold by 52.1,
  32.8 and 41.4 standard errors. Controls 0.431, 0.470, 0.431 (ideal 0.5); both/B 0.808. Results in `run/`.
- `run/phase111_hardware_run_executed.ipynb` is the executed notebook with the IBM instance and plan identifiers
  replaced by placeholders; all outputs are otherwise kept.
