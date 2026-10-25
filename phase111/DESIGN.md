# Phase 111 design: replication and decoded recovery on hardware

**Author:** Luis Razo (EISM, Barcelona). Drafted with AI assistance (Claude) under the author's direction.
**Status:** FROZEN, 2026-10-25, with the live build described under "Frozen build" (the timestamp of
`FREEZE_MANIFEST.txt` is recorded in `RUNLOG.md`). Draft 2 adopted all six fixes of a first review of the design (every review in this phase was by ChatGPT, directed and assessed by the author); draft 3 aligned the text with
the built and tested procedure (placement rule, shared layout per region, routed baseline, software versions, power
result) and added stated limitations; draft 4 adopted the blocking fixes of a second review, of the code
(enforced checks, a separate verifier, power check on the live calibration, jobs bound to the frozen manifest), and
the fixes of its follow-up (checks and power against the saved compilation target, run logic in a hashed module,
timestamp record, submitted-circuit check, Test B circuit identity).
Freezing happens after the circuits are built and verified (Part C): the design and every artifact it fixes are
listed with their SHA-256 in `FREEZE_MANIFEST.txt`, committed to the code repository, and the manifest is
timestamped with OpenTimestamps before any QPU job. Changes after freezing are logged in `RUNLOG.md` as amendments,
with their reasons, and re-stamped before any further QPU job.

## Purpose

Phases 106 and 107 support the paper *Goal-Weighted Quantum Codes for Logical Qubits of Unequal Importance, with a
Hardware Test*. A review of the paper accepted the numerical results and the entanglement witnesses, and raised three
points this phase answers:

1. **One run.** Part A repeats the hardware test.
2. **Chronology.** The earlier designs' hashes were self-recorded. This design and its artifacts are timestamped
   independently before anything runs. A timestamp shows the commitment existed by that time; it does not by itself
   show that no earlier run took place, so every IBM job ID from the account during this phase, including failed or
   cancelled jobs, is listed in `RUNLOG.md`.
3. **Entanglement, not recovery.** Part B decodes each logical qubit from a region and measures how well it returns.

The claim Part B can support: *a region-local Clifford decoder recovers the target logical qubit, and we measured its
entanglement fidelity on superconducting hardware.* It does not test robustness to an actual qubit-loss event.

## Part A: replication

### A1. Search (laptop, exact)

Rerun from scratch with the released code: `p106.py validate`, `p106.py enumerate`, `p106.py analyse`, `p107.py all`
(software versions as in `requirements.txt`; recorded in `RUNLOG.md`).

- **A1a, scientific invariants (exact):** the sets of optimal profiles, the recoverability classification of the
  Test B code, the optimal J values, the switch point at n = 6, p = 0.1, and every quantity checked by `reproduce.py`.
- **A1b, numerical agreement:** every numerical field of the regenerated `profiles.json` and both `results.json`
  agrees with the stored files within 1e-12, after canonical ordering (records sorted by their profile tuples; tied
  representatives compared as sets).
- **Pass:** A1a and A1b hold. Byte equality of files is not required. Any difference is reported with its cause.

### A2. Hardware witness test (ibm_fez)

- **A2a, circuit replay (primary):** resubmit the archived compiled circuits (`testB_transpiled.qpy`, SHA-256 in the
  first run's manifest), with the same physical layout, measurement mapping and 4,000 shots each.
- **A2b, independent compilation (only if A2a cannot run):** if the backend rejects the archived circuits because its
  target has changed, regenerate them from `testB_plan.json` with the first run's compilation procedure and seeds,
  verify them noiselessly, and run those. This is reported as a new compilation, not a circuit replay.
- **Criteria (unchanged from Test B):** B1 E_B(main) > 1, B2 E_A(both) > 1, B3 E_B(main) − E_A(main) > 0, each by at
  least three standard errors on raw counts. Replication means all three hold.
- **Reported, not tested:** all six witnesses with 95% intervals, beside the first run's values and their differences.

## Part B: decoded recovery

### What is measured

For each of the six (region, target) pairs of Test B, a decoder acting only on the four physical qubits of the region
moves the target's information onto one qubit s of the region. The reference R and s are measured in three bases, and
the entanglement fidelity with |Φ+⟩ is

F = (1 + c_XX − c_YY + c_ZZ) / 4,

where c_PP is the measured expectation of P_R P_s (decoders include the Pauli correction that makes |Φ+⟩ the ideal
state). The other logical qubit's reference is not acted on, so F is the entanglement fidelity of the target's
effective channel with the other logical input maximally mixed. The average fidelity (2F + 1)/3 is reported for the
same channel: averaged over the target's pure inputs, with the other logical input maximally mixed; it is not an
average over arbitrary joint inputs.

### Ideal values (from `testB_plan.json`, checked in the first review; ideal F recomputed on the compiled circuits by `p111_verify.py`)

| Region | Target | I(R : S) | Partners inside the region | Ideal F |
|---|---|---|---|---|
| main [2,3,6,7] | B | 2 | X: Y2Y3Z6X7, Z: X2Z6 | 1 |
| main [2,3,6,7] | A | 1 | Z: X3Z7 | 1/2 |
| both [2,3,4,5] | A | 2 | X: X2Z3Z4, Z: X4Z5 | 1 |
| both [2,3,4,5] | B | 2 | X: Z2X3Z5, Z: X3Y4Y5 | 1 |
| neither [3,4,5,7] | A | 1 | Z: X3Z7 | 1/2 |
| neither [3,4,5,7] | B | 1 | Z: Z4X5Z7 | 1/2 |

For I = 1 the reference–region state is separable, so no decoder on the region can exceed F = 1/2; the control
decoders below reach it, so the controls are decoded as well as the ideal task allows.

### Decoder synthesis (deterministic)

For each pair and each choice of output qubit s in S (four choices):

1. **Tableau.** Recoverable pairs: a Clifford on S mapping P_S → X_s and Q_S → Z_s. Controls: mapping Q_S → Z_s.
   The tableau is completed to a full 4-qubit Clifford by symplectic Gram–Schmidt over the Pauli basis in a fixed
   order (X then Z on qubits of S in increasing index), implemented in `p111_build.py`.
2. **Signs.** A final Pauli on s and, for the measurement, the ideal signs computed by statevector make |Φ+⟩ (or, for
   controls, ⟨Z_R Z_s⟩ = +1) the ideal output.
3. **Synthesis.** `qiskit.quantum_info.Clifford(...).to_circuit()` (qiskit 2.0.0).
4. **Compilation.** See the physical procedure below, with transpiler seeds 0–9.

Candidates: 4 output qubits × 10 seeds for each pair and each kept layout. Cost: total two-qubit gates (preparation
plus decoder) of the final compiled ZZ circuit; ties broken by depth, then lowest s, then lowest seed. The selected
decoder is the minimum **among these enumerated candidates**, not a proven minimum over all Clifford decoders. Selection uses only these
deterministic criteria; no noisy simulation is used to choose decoders. Every candidate and its cost is archived in
`p111_candidates.json`.

### Physical procedure and locality

- **Layout search (per region).** The compiler's own layouts leave the regions disconnected on the heavy-hex
  device, so the region is placed explicitly. For each 4-qubit path of the device (306 on ibm_fez), the region's
  qubits are placed on it in the 2 orderings that put the most graph edges inside the region onto device edges (ties
  by order); each other code qubit is placed, in order of its number of graph neighbours in the region, on the free
  device qubit with the smallest summed distance to its already placed neighbours (ties by index). The preparation
  is compiled with that initial layout (optimization level 3, seeds 0–2), and kept if the region's four qubits are
  still connected after routing. The 5 kept layouts with the fewest preparation two-qubit gates (ties by path,
  ordering, seed) are the candidates.
- **One layout per region.** Both targets of a region use the same layout and so the same preparation circuit,
  which makes D3 a comparison between targets on identical circuits up to the decoder. The layout chosen is the
  candidate minimising the sum of the two targets' best costs, ties by rank. The two targets still have different
  decoders, so D3 compares the two specified end-to-end procedures on the same preparation, not an intrinsic code
  property isolated from hardware; for each pair the preparation and total two-qubit gate counts, depth and physical
  output qubit are reported.
- **Decoder compilation.** The decoder is compiled separately on the coupling map restricted to the physical
  positions of S at the end of the preparation (taken from the preparation's final layout), so routing can use only
  those qubits. It is appended after a barrier.
- **Measurement.** Only R and s are measured, at their actual final physical positions; nothing else is measured.
- **Runtime options.** Dynamical decoupling and twirling are disabled; no other error mitigation is applied by the
  service. Options are recorded with the job.
- **Job identity.** All run logic is in `p111_runner.py`, frozen and hashed; the notebook checks that file's hash and
  only calls it. The experiment ID is `p111-` plus the first 12 hex digits of the timestamped manifest's SHA-256.
  Each job carries the tag `<experiment ID>-<step>` (A2, pilot, B). Immediately before each submission the runner
  re-checks every frozen file's hash, that today's target supports every circuit of the step, the step's frozen shot
  count, and that the previous step was retrieved. Before waiting it saves a job record (job ID, experiment ID, tag,
  ordered circuit names and fingerprints, shots, versions, timestamp record). If IBM holds more than one job with a
  step's tag, the run stops for manual resolution. Retrieval rejects a job whose tag, backend, record or per-circuit
  shot count does not match, and compares the circuits the service recorded as the job's input with the frozen
  fingerprints. If the service does not return them, a job submitted by this procedure is analysed with that fact
  recorded; a job found only through its tag is not analysed. The service's own record is the only independent
  evidence of the submitted circuits; where it is not available this is stated.

### Validity checks before any QPU job (all must pass)

- **V1, locality, on the final submitted ISA circuits.** Tracking the logical-to-physical map through every routing
  step: after the preparation boundary (a single barrier), every gate is a native gate acting only on the physical
  positions of S, except on R, whose gates together must equal the specified basis rotation (checked as an
  operator); no SWAP, reset or gate on more than two qubits; no qubit outside S ∪ {R} has any operation after the
  boundary; there are exactly two measurements, the last operations on their qubits, and none before the boundary.
- **V2, noiseless correctness.** The final ISA circuits, simulated exactly (statevector, after removing measurements),
  give the ideal F of the table above for every pair to within 1e-9.
- **V3, mapping.** The state at the boundary is the code's graph state on the recorded physical layout (all eight
  stabilizers equal +1 there, and any qubit used for routing is back in |0⟩), which checks the layout independently of
  the compiler's own layout record; clbit 0 measures R and clbit 1 the recorded output position in S.
- **Where the checks run.** The builder runs V1–V3 on each candidate it keeps, and then `p111_verify.py` reruns them
  on the QPY file as written (plus the calibration and baseline circuits, and the archived Test B circuits: each is,
  in order, the setting its manifest name says, checked on its exact pre-measurement state, with clbit v measuring
  logical qubit v) and writes `p111_verify_report.json`, which also records a fingerprint of every circuit.
  Instruction support is checked against the backend target saved at compilation (`p111_target.pkl`: gates, coupling
  map, error rates, durations, T1/T2; hash recorded in the plan), so the frozen build is checked against the live
  target it was compiled for. The mathematical checks do not depend on the target. Any failure of a primary pair, a baseline or the re-check stops the
  build with a nonzero exit and writes only `p111_plan_FAILED.json`, which cannot be frozen.
- **Failure consequences, fixed now.** If V1–V3 cannot be met for a **primary pair** (main/B, main/A, both/A) under
  the procedure above, Part B is not run; the failure is reported and any redesign is a new, separately frozen design.
  If they cannot be met for a **secondary pair** (both/B, neither/A, neither/B), that pair is recorded in the plan as
  not run, its circuits are left out, and the primary criteria are unaffected. No pair is dropped silently.

### Estimators and criteria (raw counts primary)

- **Correlator:** with n shots and k shots of even parity, ĉ = (2k − n)/n (unsmoothed). Its standard error uses
  p̃ = (k + 1)/(n + 2): SE(ĉ) = 2√(p̃(1 − p̃)/n). Smoothing affects only the standard error.
- **Fidelity:** F̂ = (1 + s_XX ĉ_XX + s_YY ĉ_YY + s_ZZ ĉ_ZZ)/4 with the fixed signs s (from the ideal state; for
  |Φ+⟩, s_XX = s_ZZ = +1, s_YY = −1), and SE(F̂) = ¼ √(SE_XX² + SE_YY² + SE_ZZ²). The three correlators are measured
  in separate circuits, so their errors are independent.
- **D1:** F_B(main) − 1/2 > 3 SE. **D2:** F_A(both) − 1/2 > 3 SE. **D3:** F_B(main) − F_A(main) > 3 SE, with the two
  standard errors combined in quadrature.
- **Joint rule:** Part B passes only if D1, D2 and D3 all hold. No multiplicity correction is applied, since all three
  must hold; z-values use a normal approximation and are reported.
- **Reported, not tested:** all six F with standard errors and 95% intervals; (2F + 1)/3; readout-corrected values. A
  control at or below 1/2 is consistent with the prediction; it is not proof that the target cannot be recovered,
  which is the code's exact theoretical property.

### Readout correction (secondary)

- For every distinct physical measured pair (R, s) used, four calibration circuits (prepare 00, 01, 10, 11 and
  measure the pair with the same measurement settings), in the same job, 4,000 shots each.
- Correction by inverting that pair's 4 × 4 assignment matrix. Its condition number is reported; if it exceeds 10,
  corrected values for that pair are not reported. Quasi-probabilities outside [0, 1] are reported, not clipped.
- Uncertainty by bootstrap over data and calibration counts (10,000 resamples). A resample whose matrix exceeds the
  condition cap or is singular is not used; the number used is reported, and below 99% no corrected SE is reported.
  The bootstrap covers counting noise only, not drift or a wrong assignment model.
- A corrected value never replaces a failed raw criterion.

### Baseline

For every distinct physical measured pair (R, s), a Bell pair prepared on R and its neighbour on the shortest
device path to s, then moved to s by SWAPs through the path's qubits (which start and end in |0⟩; no SWAP when R and
s are adjacent), and measured on R and s with the same basis rotations, 4,000 shots per correlator. It gives the two-qubit preparation-and-measurement
reference for that pair. It does not decompose the shortfall of F below 1 into decoding and other losses, because
the encoded circuits have a different gate history; it is reported as a reference only (path, SWAP count, physical
qubits, two-qubit gates), and no decoded fidelity is normalised by it. Its qubits may differ from those of the decoded
circuits.

### Pilot

A fixed pilot with no outcome-dependent abort: the six correlator circuits of D1 and D2 (XX, YY, ZZ for main/B and
both/A) plus their calibration circuits, 1,000 shots each. Its purpose is to measure QPU time and check the pipeline.
The full job is not submitted only for a technical failure, fixed now: the job errors, returns the wrong number of
results, or V3's bit mapping does not match the returned data. Pilot data are reported and not pooled with the full
run.

## Part C: before freezing

1. Build decoders and circuits (`p111_build.py`); write `p111_plan.json` (supports, decoders, layouts, output
   positions) and `p111_candidates.json`.
2. Run V1–V3 on the written files (`p111_verify.py`, called by the builder; it can also be run alone).
3. Power check (`p111_power.py`): density-matrix simulation of the final circuits, 200,000 shots per circuit
   (expectation error below 0.002), with the noise model built from the saved compilation target, so the circuits and
   their power check use the same calibration: for the frozen build, the live ibm_fez calibration at compilation;
   FakeFez only for the development build. `p111_power.json` records the target's hash. Shots are the
   only quantity this check may change: the smallest of 4,000, 8,000 and 16,000 at which the expected z of D1, D2 and
   D3 are all at least 5. If none reaches 5, the script stops and Part B is not run without a design amendment. The
   notebook takes the shot count from `p111_power.json`. It does not change decoders.
4. Write the offline analysis (`p111_analyse.py`), tested on simulated counts.
5. Review of this design, the circuits and the notebook.
6. Freeze: `FREEZE_MANIFEST.txt` lists the SHA-256 of DESIGN.md, `p111_build.py`, `p111_verify.py`,
   `p111_power.py`, `p111_analyse.py`, the run notebook, `p111_plan.json` (which records software versions),
   `p111_runner.py`, `p111_candidates.json`, the final ISA circuits (`p111_circuits.qpy`), the compilation target
   (`p111_target.pkl`), `p111_verify_report.json`,
   `p111_power.json`, the ibm_fez calibration snapshot used for compilation (`p111_backend_snapshot.json`), and the
   archived Test B files (`testB_transpiled.qpy`, `testB_manifest.json`, `testB_plan.json`). Runtime options are fixed
   in the notebook. Commit; stamp the manifest with OpenTimestamps; commit the proof. The whole procedure is rerun
   once from a clean checkout of the repository before freezing.
7. Before the run: the OpenTimestamps proof is verified outside the notebook, once the stamp is confirmed in a
   Bitcoin block, and the result is written to `TIMESTAMP_CHECK.json` (manifest hash, attestation time, block,
   method), kept with the proof and committed. The runner refuses to start without it, or if its manifest hash
   differs, and every job record includes it, so the attestation time can be compared with the job creation times.
   The runner checks every file listed in the manifest (the notebook itself, which Colab rewrites, is checked in the
   final audit against its frozen copy instead; it contains no run logic), asks for the manifest's SHA-256 from the
   timestamped commit, and checks that the verification report passes for these exact files and target, that the
   circuit names and fingerprints match it in order, that the power check is for these circuits and this
   calibration, that the software versions match the build, and that today's target supports every instruction. Run within 72 hours of compilation; if the backend rejects the frozen circuits,
   the A2b rule applies to Part A, and for Part B the circuits are recompiled by the same procedure, re-verified,
   logged as an amendment and re-stamped before submission.

## Software, development build and frozen build

- Python 3.12.3; qiskit 2.6.1, qiskit-aer 0.18.0, qiskit-ibm-runtime 0.51.0 (the versions that ran Test B; QPY
  files need a Qiskit at least as new as their writer), numpy 2.5.4, networkx 3.6.1. Recorded in `p111_plan.json`.
- The **development build** targets the FakeFez model of ibm_fez. It fixed and tested the procedure: all six pairs
  pass V1–V3 exactly; the power check gives expected z of about 87, 88 and 54 for D1–D3 at 4,000 shots (the model is
  optimistic; Test B's hardware witnesses fell below its simulation), so 4,000 shots are kept; and the notebook ran
  end to end against a simulated service (three submissions, none repeated after reruns or lost local records).
- The **frozen build** runs the same procedure against the live ibm_fez target (`p111_build.py --live`) and saves the
  backend's calibration (`p111_backend_snapshot.json`). Its circuits may differ from the development build's; they
  must pass V1–V3 and the power rule (on the live calibration) before freezing, and only they are run; the notebook
  refuses a development build.
- **Frozen build, 2026-10-25.** Compiled at 16:07 (UTC+1) against the live ibm_fez target, calibration of
  2026-10-25 14:45 (UTC+1). It selected the same layouts, decoders, output positions and two-qubit gate counts as the
  development build; the compiled preparation of the both region differs in single-qubit gates (six circuits). All
  six pairs pass V1–V3, all baselines pass their exact check, and `p111_verify.py` passes on the written files (60 circuits; the
  20 archived Test B circuits identified). Power check on that calibration: expected F 0.812 (main/B), 0.466
  (main/A), 0.874 (both/A); expected z of D1–D3 at 4,000 shots 68.8, 98.6 and 46.5, so 4,000 shots. To be run
  within 72 hours of compilation.
- SamplerV2 is deprecated in qiskit-ibm-runtime 0.51.0 but available; its availability is checked before the run.

## Hardware budget (provisional)

Part A: 20 circuits × 4,000 shots (about 25 QPU seconds in the first run). Part B: 18 decoding circuits, one baseline
of three circuits and one calibration set of four circuits per distinct measured pair, × 4,000 shots (more if the
power check raises shots). The estimate is computed from the final compiled circuits' durations in Part C, with a
margin for the pilot and calibration changes; until then, 30–60 QPU seconds for Part B is a plausible range, not a
verified figure. Part A runs first, then the Part B pilot, then Part B, as separate jobs, each submitted once with
its job record saved before waiting.

## What a result will and will not show

- A Part B pass shows that a decoder confined to the region recovers the target with entanglement fidelity above the
  classical bound, on the same region that does not yield the other target.
- It does not show recovery from a physical loss event, protection over time, or worst-case recovery of joint inputs.
- The qubits outside the region are not acted on after preparation, but they remain physically present; always-on
  coupling (crosstalk) between them and the region cannot be ruled out and is a stated limitation.
- Readout calibration circuits are run with the other qubits idle, unlike the data circuits.
- A Part B failure with a Part A pass means the entanglement is present but this device does not decode it well
  enough with these circuits; both outcomes are reported.
