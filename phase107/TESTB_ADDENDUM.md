# Phase 107 Test B addendum (frozen before any noisy simulation or hardware job)

Written 2026-10-23 after `testB_select.py` chose the wiring and regions (noise-free, deterministic), before
`testB_sim.py` existed. Hash in `RUNLOG.md`. This fixes what DESIGN.md left to an addendum.

## Code, wiring and regions (from `testB_plan.json`)

- **Code:** the Phase 106 shared [[6,2]] optimum at n = 6, p = 0.1. Eight qubits: references A (qubit 0), B (qubit 1),
  and physical qubits 2–7.
- **Wiring:** prepared as a graph state, using the member of its local-complementation orbit (3,088 graphs, the 60 with
  fewest edges compiled with 5 seeds each) that compiles to the fewest two-qubit gates on ibm_fez's layout. That is
  11 edges, 20 two-qubit gates after routing, depth 21. Recoverability is unchanged by local unitaries.
- **Change from DESIGN.md:** in this code every 4-qubit region that recovers A also recovers B, so the planned
  "A but not B" flip control does not exist. The regions are therefore, with I the mutual information in bits:
  - **main [2, 3, 4, 7]:** B fully recoverable (I = 2); A only classically correlated (I = 1). Same region, two targets.
  - **both [2, 3, 4, 5]:** A and B both fully recoverable. This control shows A is recoverable on this hardware, so a
    failure on the main region is not just "A is noisier".
  - **neither [2, 3, 5, 6]:** A and B both only classically correlated (I = 1).
  - The first region of each kind in lexicographic order was taken.

## What is measured

For each (region, target) pair:
- **Correlators:** three of the form σ_ref ⊗ P, one for each σ ∈ {X, Y, Z}. P is the lowest-weight stabilizer partner,
  truncated to the region when it extends outside it.
- **Ideal values:** ±1 when the partner lies inside the region, 0 otherwise (signs from the statevector).
- **Circuits:** 18 settings × 4,000 shots, plus 2 readout-calibration circuits (all 0, all 1) × 4,000 shots.

Two statistics, with c_σ the measured correlator times its ideal sign:

- **E = c_X + c_Z** (entanglement witness): P_X and P_Z anticommute on the region, so E > 1 certifies entanglement
  between the reference and the region, i.e. quantum recovery. Ideal values:
  - 2 when the target is fully recoverable;
  - 1 when it is classically correlated only (here A on main and on neither, where only the Z partner lies inside);
  - for B on neither, only Z lies inside, so 1 as well.
- **W = (1 + c_X + c_Y + c_Z)/4** (recovery score): ideal 1 for full, 0.5 for classical only. Reported, not tested.

Standard errors: Var⟨P⟩ = (1 − ⟨P⟩²)/shots, so Var(E) is the sum of the two variances.

## Pass criteria (hardware run; raw counts are primary, readout-corrected values are secondary)

- **B1:** E_B(main) > 1 by at least 3 standard errors. The region holds B's quantum information.
- **B2:** E_A(both) > 1 by at least 3 standard errors. A's quantum information is recoverable elsewhere on the same
  run.
- **B3:** E_B(main) − E_A(main) > 0 by at least 3 standard errors. On the same region, B beats A.
- **Pass:** B1, B2 and B3 all hold. E_A(main) and both "neither" values are expected at or below 1. They are reported,
  but absence of entanglement cannot be certified by a witness, so they are not tested.

## Simulation gate before hardware

- **What runs:** the same 20 circuits on Qiskit Aer with the FakeFez noise model (qiskit-ibm-runtime 0.37.0), at the
  same shot counts, with transpile seed and layout as chosen above, plus a noise-free run.
- **Go:** the noisy simulation passes B1–B3 on raw counts.
- **No-go:** it fails. Then Test B moves to an all-to-all device (IonQ) and is not run on IBM.

## Hardware run (only after the simulation gate passes and the user approves)

- **Pilot:** one small job (2 circuits × 1,000 shots) to measure actual QPU usage.
- **Full run:** within the user's 10-minute allowance, on the backend the user chooses, transpiled the same way for that
  backend.
- **No changes after seeing data;** any deviation is logged.
