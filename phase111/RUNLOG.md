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
