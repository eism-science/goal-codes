# goal-codes

Code and data for L. Razo, *Goal-Weighted Quantum Codes for Logical Qubits of Unequal Importance, with a Hardware
Test* (EISM, Barcelona).

The paper asks which stabilizer code best protects two logical qubits of unequal importance under erasure, and whether
the answer changes with the importance weights. It answers exactly for every stabilizer code on up to six physical
qubits, compares with a classical benchmark that may copy, repeats the search under uneven and correlated erasure, and tests
the selected code on IBM hardware (ibm_fez): an entanglement-witness test, its preregistered replication, and decoding
confined to one region of the code.

## Reproduce the paper's numbers

```
pip install -r requirements.txt
python reproduce.py
```

This checks every number in the paper's tables and text against the stored data in a few seconds, with no IBM account.
For each number it prints the paper's value, the recomputed value, and OK if they agree at the paper's precision. It
also checks the SHA-256 of the design documents and of the compiled hardware circuits against the recorded values, and
every Phase 111 file against its freeze manifest.
`python reproduce.py --recompute-noise` also rebuilds the noise-model table from the stored code tables (slower).

## Contents

| Folder | What it holds |
|---|---|
| `phase101/` | `p101_core.py`, an independent recoverability implementation used to validate Phase 106, and `gf2_toolbox.py`, which it needs. |
| `phase106/` | Identical erasure (paper Section 3): exhaustive search (`p106.py`), stored code profiles (`profiles.json`), results, the exploratory n = 7 search and a random-graph check. |
| `phase107/` | Uneven and correlated erasure (Section 4, `p107.py`, stored code tables `tables_cache.npz`), and the hardware test (Section 5): circuit selection (`testB_select.py`), simulation (`testB_sim.py`), the hardware notebook and `IBM run data/` (raw counts, analysis, manifest, compiled circuits as QPY). |
| `data/` | The ibm_fez calibration snapshot (2026-05-25) from which the readout-weighted erasure model is built. |
| `phase111/` | Replication of the hardware test and a decoded-recovery test, frozen before running: design (`DESIGN.md`), circuit build (`p111_build.py`), independent checks (`p111_verify.py`), power check (`p111_power.py`), analysis (`p111_analyse.py`), run logic (`p111_runner.py`) and notebook, offline tests, and the built circuits with their plan, verification report and the ibm_fez calibration they were compiled for, the freeze manifest with its timestamp, and `run/` (raw counts,
results, job records and the executed notebook). `p111_target.pkl` is a Python pickle; the scripts check its SHA-256 against the plan before loading it, and it should not be loaded from any other source. |

To rerun from scratch rather than from stored data:

```
python phase106/p106.py validate
python phase106/p106.py enumerate     # rewrites phase106/profiles.json
python phase106/p106.py analyse       # rewrites phase106/results.json
python phase107/p107.py all           # builds the code tables, validates, rewrites phase107/results.json
python phase107/testB_select.py       # rewrites phase107/testB_plan.json
python phase107/testB_sim.py sim
```

The full Phase 107 run needs several GB of memory at six physical qubits. The hardware notebook needs an IBM Quantum
account (API key from the `IBM_QUANTUM_TOKEN` environment variable or a prompt; set `INSTANCE` to your own instance).
Rerunning it submits new jobs; the stored results come from job `dbdpp41gv2cc73dotnv0`. The Phase 111 results come from
jobs `dbf58r8m94ac738ssbi0` (replication), `dbf59bnsstms73fd4ntg` (pilot) and `dbf59tvpc3tc739mp16g` (decoding); see
`phase111/RUNLOG.md`.

## Records

Each phase has a design document (`DESIGN.md`, plus `TESTB_ADDENDUM.md` for the hardware test), a run log
(`RUNLOG.md`) and a results summary (`RESULTS.md`). The design documents and run logs are unedited working records.
They mention a grant application the work was planned for, folder locations on the author's machine, and which steps
were drafted with AI assistance (see the paper's statement on AI tools). The run logs
report that each design's SHA-256 hash was recorded before its analysis ran; for Phases 106 and 107 this chronology was
not independently timestamped. Phase 111 was frozen before running: `phase111/FREEZE_MANIFEST.txt` lists the SHA-256 of
every file the run used, and its OpenTimestamps proof (`FREEZE_MANIFEST.txt.ots`, checkable with `ots verify`) places it
in Bitcoin block 973103, before any Phase 111 job was submitted. Amendments after
a design was fixed are listed in the run logs. Edits made for this release, none affecting any result: file paths
changed to the layout of this repository; the IBM instance identifier removed from the hardware notebook, and the IBM instance and plan identifiers replaced in
the executed Phase 111 notebook; internal
references removed from one docstring and from the two results summaries, whose sections on the grant application
were dropped.

## Licences

- Code: MIT (`LICENSE`).
- Design documents, logs, results and run data: CC BY 4.0 (`LICENSE-data`).
- The ibm_fez calibration snapshot is IBM backend data, included unmodified for reproducibility (see `LICENSE-data`).
