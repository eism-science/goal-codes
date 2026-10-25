# Phase 106 results

Run 2026-10-23, overnight and unattended, locally. Design `DESIGN.md` (frozen before any code existed; hash in
`RUNLOG.md`). Code `p106.py`; outputs `profiles.json`, `results.json`, `analyse_output.txt`. Every number below is
exact (exhaustive over all stabilizer codes with n ≤ 6 physical qubits, up to local Cliffords and qubit relabelling),
except the exploratory n = 7 section.

## In one paragraph

In the primary cell (6 qubits, erasure 0.2), sharing never pays: for every goal, from evenly split to 0.99 on one
qubit, the best design is the 5-qubit code for one logical and a bare qubit for the other. So the main question (does
sharing survive a lopsided goal?) is moot there, because there is nothing to survive. Where sharing does pay (lower
erasure, or 4 to 5 qubits), it shrinks as the goal becomes lopsided and, at 6 qubits and erasure 0.1, disappears
entirely beyond w ≈ 0.82, where the best design switches from a shared 6-qubit code to the 5-qubit code plus a bare
qubit. That switch is a genuine goal-driven change of architecture, found exactly, but in a secondary cell. The one
exception is 4 qubits, where no separate design can protect anything, so the shared [[4,2,2]] code wins at every goal.

## Frozen predictions (primary cell n = 6, p = 0.2)

| Prediction | Verdict | Numbers |
|---|---|---|
| P1 optimal design changes across w | **holds nominally only** | The only "change" is after w = 0.50, where the tie class contains both orientations (5-qubit code on A or on B); from 0.55 on, one design. No substantive change. |
| P2 sharing gain at w = 0.5 > 0.01 | **fails** | gain 0 (2e-16) |
| P3 sharing gain at w = 0.9 > 0.01 (predicted to fail) | **fails, as predicted**, but vacuously | gain 0 at every w, so it did not shrink: it was never there |
| P4 gain non-increasing in w | holds (trivially) | constant 0 |
| P5 quantum gain > classical gain at w = 0.5 | **fails** | quantum 0, classical 0.0049 |

The best design in the primary cell at every w: the [[5,1,3]] code on A (P_A = 0.9421) and a bare qubit on B
(P_B = 0.8). J at w = 0.5 is 0.8710 for both the best of all codes and the best separate design.

## Sharing gain across cells (quantum | classical), descriptive

| n, p | w = 0.5 | 0.7 | 0.9 | 0.99 | where the best design changes |
|---|---|---|---|---|---|
| 4, 0.1 | 0.048 \| 0.004 | 0.048 \| 0.005 | 0.048 \| 0.007 | 0.048 \| 0 | never ([[4,2,2]] always) |
| 5, 0.1 | 0.024 \| 0.004 | 0.014 \| 0.003 | 0.005 \| 0.001 | 0.0005 \| 0.001 | never |
| 6, 0.1 | **0.025** \| 0.001 | 0.010 \| 0.001 | **0** \| 0.001 | 0 \| 0.0001 | **at w = 0.823: shared [[6,2]] → [[5,1,3]] + bare** |
| 4, 0.2 | 0.019 \| 0.010 | 0.019 \| 0.016 | 0.019 \| 0.010 | 0.019 \| 0 | never |
| 5, 0.2 | 0.010 \| 0.015 | 0.006 \| 0.008 | 0.002 \| 0.005 | 0.0002 \| 0.0003 | never |
| 6, 0.2 | **0** \| 0.005 | 0 \| 0.005 | 0 \| 0.003 | 0 \| 0.001 | never (tie only at 0.5) |
| any n, 0.3 | 0 \| up to 0.04 | 0 | 0 | 0 | never |
| n ≤ 3 | quantum 0 everywhere | | | | |

Reading:

- **Sharing is quantum-specific only in a narrow window.** At erasure 0.1 the quantum gain is 5 to 30 times the
  classical one (n = 4 to 6). At erasure 0.2 they are comparable, and at 0.3 no quantum code beats a bare qubit
  except the 5-qubit code, while classical sharing (parity bits) still helps.
- **Where sharing pays, a lopsided goal removes it,** unless no single-logical code is available at that size (n = 4).
  This is the direction P3 predicted, seen in the secondary cells; in the primary cell it is vacuous.
- **The goal does change the architecture,** exactly and not by tie-break, at n = 6 and p = 0.1. The shared code there
  is a [[6,2]] code with unequal protection: neither logical recovers from any 3 qubits; A recovers from 14 of the 15
  four-qubit sets and B from 12 (profiles in `results.json`). Below w = 0.823 it beats the 5-qubit code plus a bare
  qubit; above, it loses.

## Cost per qubit (best size from 2 to 6, descriptive)

- p = 0.1, λ = 0.02: the shared [[4,2,2]] code is best for every w up to 0.95; at w = 0.99 the best is 6 qubits,
  separate ([[5,1,3]] + bare). The goal changes both the size and the architecture.
- p = 0.2 and 0.3, λ = 0.02: two bare qubits for w ≤ 0.55; the 5-qubit code plus a bare qubit for w ≥ 0.6. The goal
  decides whether protection is worth four extra qubits.
- λ = 0.05: two bare qubits at every goal. λ ≤ 0.01: always 6 qubits.

## Validity

- V1 to V5 passed before the grid was computed: graph-state recoverability matched `p101_core.Code.recoverable` on
  2,500 random codes and 74,400 regions with no mismatch; [[5,1,3]], [[4,2,2]] and bare-qubit profiles reproduced;
  133,632 graphs enumerated for n = 6; no single-logical profile exceeds the no-cloning bound.
- Post hoc: 300,000 random labelled 8-vertex graphs never beat the exhaustive n = 6 optimum at any (p, w).
- Deviation 1 (bug fix, logged): the classical "all codes" family first included codes that never store bit B. Fixed so
  both bits must be stored, as in every quantum code. No primary verdict changes. The unfixed run is kept
  (`results_classical_unrestricted.json`).

## Limits

- Stabilizer codes only, independent erasure only, n ≤ 6 exhaustive. Real hardware noise is not erasure.
- The primary cell was chosen before running; that it turned out to be the one cell where a perfect single-logical
  code fits exactly (5 + 1 = 6) is part of the result, not a reason to replace it.
- The thresholds and the "w ≈ 0.82" crossover are properties of these small codes; no scaling claim is made.

## Exploratory: n = 7 (search, not exhaustive)

`p106_explore7.py`, `explore7.json`; 3 minutes. The separate side is exact (single-logical codes up to 6 qubits are
exhaustive); the all-codes side is the best found by 100,000 random 9-vertex graphs, the extensions of every n = 6
optimum, and greedy edge-flip ascent. A positive gain is real; a zero gain means only that none was found.

| p | w = 0.5 | 0.6 | 0.7 | 0.8 | 0.9 | 0.99 |
|---|---|---|---|---|---|---|
| 0.1 | +0.031 | +0.022 | +0.014 | +0.006 | 0 | 0 |
| 0.2 | +0.0006 | 0 | 0 | 0 | 0 | 0 |
| 0.3 | 0 | 0 | 0 | 0 | 0 | 0 |

The same pattern holds at 7 qubits. At erasure 0.2, sharing appears at a balanced goal but is negligible (0.0006,
far below 0.01) and is gone by w = 0.6. At 0.1 it is larger at a balanced goal and gone by 0.9. The primary null is not
an artefact of 6 = 5 + 1. The best separate design at n = 7 is no better than at n = 6 at any p, because no
6-qubit single-logical code beats the 5-qubit code with a spare qubit, and no 2-qubit code beats a bare qubit.
