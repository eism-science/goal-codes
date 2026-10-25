# Phase 107 Test A results

Run 2026-10-23, locally, exact (every stabilizer code with n ≤ 6 physical qubits, at every placement, under each noise
model). Design `DESIGN.md` (hash in `RUNLOG.md`). Code `p107.py`; output `results.json`, `run_output.txt`.

## In one paragraph

The bid rule passes, but only just and only through the stylised correlated-noise model. With correlated erasure on
a ring (M2) at the primary cell (6 qubits, mean erasure 0.2), knowing the goal is worth 0.0148 on average over lopsided
goals, against exactly 0 under identical noise. The mechanism is the Phase 106 low-erasure switch, moved to the
primary noise level: at a balanced goal a shared code protecting both logicals equally (0.8708 each) is best; from
w ≈ 0.7 the best design becomes the 5-qubit code on A plus a bare qubit on B. With noise copied from the device
(M1), the goal is worth 0.0049 (median of 20 qubit paths on ibm_fez, range 0.0009–0.0166, 2 of 20 above 0.01), so P1
fails there. The effect does not survive at the other noise levels: M2 gives 0.0024 at erasure 0.1 and 0 at 0.3.

## Predictions (primary cell n = 6, p̄ = 0.2)

| Prediction | M1 device-uneven | M2 correlated | M1s stress |
|---|---|---|---|
| P1 V̄ ≥ 2·V̄(M0) and ≥ 0.01 | **fails** (0.0049 vs 0) | **holds** (0.0148 vs 0) | fails (0) |
| P2 more design switches than M0 | holds (2 vs 0) | holds (1 vs 0) | fails (0) |
| P3 quantum V̄ > classical V̄ | holds (0.0049 vs 0.0002) | holds (0.0148 vs 0.0002) | fails (0 vs 0.0001) |

**Bid rule (P1 for M1 or M2): PASS, via M2 only.**

My recorded expectation was the reverse: P1 a coin flip for M1 and less likely for M2. M2 is the model that passed.

## Value of knowing the goal, 6 qubits (quantum | classical), descriptive

| p̄ | M0 identical | M1 device (median) | M1s stress | M2 correlated |
|---|---|---|---|---|
| 0.1 | 0.0031 \| 0 | 0.0027 \| 0 | 0.0005 \| 0 | 0.0024 \| 0 |
| 0.2 | 0 \| 0.0004 | 0.0049 \| 0.0002 | 0 \| 0.0001 | **0.0148** \| 0.0002 |
| 0.3 | 0 \| 0.0027 | 0.0030 \| 0.0010 | 0 \| 0.0005 | 0 \| 0.0012 |

For n ≤ 5 the quantum value is 0 in every model and cell except one (M1s, n = 5, p̄ = 0.1: 0.0088).

Under M2 at the primary cell, V(w) grows steadily from 0 at w ≤ 0.65 to 0.0365 at w = 0.99. Under the device paths it
appears only for very lopsided goals (w ≥ 0.9, up to 0.017 at 0.99).

## Reading

- **Where the goal matters in quantum codes, it does so through sharing.** Every non-zero value comes from a shared
  code being best at a balanced goal and a separate design being best at a lopsided one. Noise that makes sharing pay
  (low erasure in Phase 106, correlated erasure here) is what makes the goal matter.
- **Quantum beats classical where the effect exists.** In the cells where the quantum value is non-zero, the classical
  value is near zero, and vice versa at high erasure, where classical parity-style sharing helps and quantum sharing
  does not.
- **Uneven noise alone does not do it.** The strong-spread stress model gives 0: when some qubits are much better than
  others, the best design puts the valued logical on the good qubits at every goal.
- **Fragility.** The pass rests on one stylised model at one noise level. A reviewer would reasonably call this a
  narrow window, not a general effect.

## Validity

- **V1:** reproduces Phase 106 (worst difference 5.6e-16).
- **V2:** every distribution and marginal is exact.
- **V3:** 300 random labelled codes, none above the optimum, all matched to the enumerated set, with probabilities
  agreeing with an independent event-by-event computation.
- **V4:** no negative sharing gain.
- The first full run was stopped by the system for low memory. It was restarted with a smaller chunk size, and the
  results are unaffected (RUNLOG).

## Limits

- Six qubits, erasure noise only, and stabilizer codes only.
- M2 is a stylised ring model with uniform marginals, not a fitted device model. M1 uses readout errors as a stand-in
  for per-qubit loss.
- "Goal-blind" means the best balanced-goal design, given the best member of its tie class.

## Interpretation

Test A passes the pre-set rule. But the evidence that realistic noise
makes the goal matter is thin:
- one stylised model;
- one noise level;
- device-derived noise fails the bar.

The honest claim is narrower than "realistic noise makes goals matter": **correlated noise can make shared codes win
at balanced goals, and then the goal decides between shared and separate designs.** Whether real device correlations
do this is an open question that hardware characterisation could answer.

---

# Test B: hardware result (2026-10-23)

Run on **ibm_fez** (IBM Open plan) with the notebook
`phase107_testB_hardware_run.ipynb`, under the addendum (`TESTB_ADDENDUM.md`) and the
procedure changes logged in `RUNLOG.md`. Pilot job `dbdpohk24mrc739iq8sg` (4 QPU s); full job 20 circuits ×
4,000 shots (23 QPU s). Raw counts, analysis, manifest and compiled circuits are in `IBM run data/` (full job `dbdpp41gv2cc73dotnv0`;
the compiled circuits' SHA-256 matches the manifest; an independent re-analysis of the raw counts reproduces every
witness exactly).

| Region / target | E raw (SE) | Ideal | E corrected |
|---|---|---|---|
| main [2,3,6,7] / B | 1.521 (0.014) | 2 | 1.618 |
| main [2,3,6,7] / A | 0.812 (0.019) | 1 | 0.956 |
| both [2,3,4,5] / A | 1.550 (0.014) | 2 | 1.859 |
| both [2,3,4,5] / B | 1.548 (0.014) | 2 | 1.641 |
| neither [3,4,5,7] / A | 0.811 (0.018) | 1 | 0.955 |
| neither [3,4,5,7] / B | 0.835 (0.018) | 1 | 0.890 |

**Criteria fixed before the run (raw, primary): all pass.** B1: E_B(main) = 1.521, 35.9 SE above 1. B2: E_A(both) =
1.550, 39.2 SE above 1. B3: E_B(main) − E_A(main) = 0.709, 30.0 SE above 0. Readout-corrected values also pass.

**What it shows.** On real hardware, the same 4-qubit piece of the code retains entanglement with B's reference
(E > 1, certified by the anticommuting-pair witness) and not detectably with A's (E = 0.81, as predicted for a
classical trace), while A's entanglement is retained by another piece. Every witness lands in the predicted class:
the three "recoverable" cases at 1.52–1.55, the three "classical trace" cases at 0.81–0.84, with uniform loss to
hardware noise (ideal 2 and 1). This certifies the predicted entanglement asymmetry, not full recoverability or
non-recoverability; the latter is the exact theoretical prediction for this code.

**Limits.** One code, one device, one run; 8 qubits; normal-approximation z-scores (adequate here, not exact
p-values). The circuit is Phase 106's [[6,2]] code, not a holographic code.
