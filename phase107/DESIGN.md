# Phase 107 (design, frozen before any run): does realistic noise make the goal matter more?

Status: written 2026-10-23 by Claude, at the user's request, before any Phase 107 code existed. Frozen by the SHA-256
recorded in `RUNLOG.md`. Anything changed after freezing is listed in `RUNLOG.md` as a deviation, with its reason.
Two tests are designed here. Test A runs now (laptop, exact). Test B (hardware) is designed but runs only if Test A
passes its bid rule and only with the user's go-ahead; its hardware parameters are fixed in a later addendum before
any hardware job is submitted.

## Why

Phase 106 (`../phase106/RESULTS.md`) found, under independent and identical erasure, that the goal changes the best
code design only in a narrow window (6 qubits, erasure 0.1, switch at w = 0.823), and that in the primary cell
(6 qubits, erasure 0.2) knowing the goal is worth nothing: the best design for a balanced goal is also best for every
lopsided goal. Identical noise on every qubit is the case where code structure matters least. Real devices have
uneven and correlated noise. If the goal matters only under identical noise in narrow windows, "learning what to keep
when nothing can be copied" is a niche and should not be the centre of the QML bid. If uneven or correlated noise makes
the goal matter clearly more, it can be.

The quantum-specific reason the goal matters here (and not merely for estimating averages, where classical shadows
let targets be chosen after measuring): for an unknown quantum state no two disjoint sets of qubits can both recover
the same logical qubit, while a classical scheme can let both recover the same bit. Which qubits keep which target is
a choice that cannot be deferred.

## Setting (as Phase 106, except the noise)

- Two logical qubits, A (valued) and B. Goal weight w on A, 1 − w on B, w ≥ 0.5.
- n physical qubits. A noise model gives a probability for every surviving set S of physical qubits (the others are
  erased). P_X = Σ_S rec_X(S) · Pr(S survives exactly). J(w) = w·P_A + (1 − w)·P_B.
- rec_X(S) as in Phase 106: purified stabilizer codes as graph states, I(R_X : S) = 2 via cut-ranks over GF(2).
- Under uneven or correlated noise the placement of a code on the physical qubits matters, so every code is evaluated
  at every placement (all n! relabellings of the physical qubits). Everything is exact: no sampling, no search.

## Families (all exhaustive, n ≤ 6)

- **ALL(n):** every [[n,2]] stabilizer code (graphs on n + 2 vertices as in Phase 106, every ordered reference pair,
  valid when S(R_A R_B) = 2), every placement.
- **SEP(n):** a [[m,1]] code for A on m qubits and a [[n − m,1]] code for B on the rest, every single-logical code,
  every split and placement.
- **Classical twin (ALL_cl, SEP_cl):** as Phase 106 after its Deviation 1: every multiset of the 7 coordinate types
  that stores both bits, every placement; SEP_cl uses only copies of x_A and x_B.

## Noise models (all with the same mean erasure probability p̄ per qubit as the identical baseline)

- **M0, identical:** every qubit erased independently with probability p̄ (Phase 106).
- **M1, device-uneven:** independent erasure with p_i = p̄ · e_i / mean(e), where e_i are the readout errors of a
  path of n coupled qubits on ibm_fez, from the calibration snapshot
  `qec/data/iteration_depth_30/ibm_fez/calibration_snapshot.json` (2026-05-25). Twenty paths are drawn with seed
  107: start at a uniformly random qubit, extend by a uniformly random unvisited cz-neighbour of the current end,
  restart if stuck, until n = 6 qubits. A path is rejected and redrawn if max(e)/mean(e) ≥ 3.1 (so that every p_i stays
  below 0.95 at p̄ = 0.3); rejections are logged. For n < 6 the first n qubits of each path are used. The primary
  statistic for M1 is the median over the twenty paths.
- **M1s, stylised strong spread:** independent, p_i = p̄ · f_i with f evenly spaced from 0.1 to 1.9 (stress case).
- **M2, correlated:** qubits on a ring (n edges; one edge for n = 2). Each qubit is erased independently with
  probability q, and each ring edge independently erases both its qubits with probability r, with
  1 − q = sqrt(1 − p̄) and (1 − r)^deg = sqrt(1 − p̄), deg the ring degree (2 for n ≥ 3, 1 for n = 2). Every qubit
  then has marginal erasure p̄, half of its survival hazard coming from shared events.

## Quantities

- J*(w): the best J over ALL(n). The optimal set at w is the tie class of (P_A, P_B) pairs within 1e-12 of J*(w).
- **Value of knowing the goal**, V(w) = J*(w) − max over the w = 0.5 tie class of J(w). The goal-blind design is the
  best design for a balanced goal, and it is allowed the best member of its tie class (the generous choice, which
  makes V harder to be large). V̄ = mean of V(w) over w ∈ {0.55, 0.60, …, 0.95, 0.99} (10 points).
- **Design switches:** the number of adjacent pairs in the w grid from 0.55 to 0.99 where the optimal (P_A, P_B) tie
  class changes. Starting at 0.55 excludes the nominal w = 0.5 orientation tie that made Phase 106's P1 hold vacuously.
- **Sharing gain** G(w) = max over ALL(n) − max over SEP(n) (descriptive, as Phase 106).
- The same three quantities for the classical twin.

Grid: w ∈ {0.50, 0.55, …, 0.95, 0.99}; p̄ ∈ {0.1, 0.2, 0.3}; n ∈ {2, …, 6}. **Primary cell: n = 6, p̄ = 0.2.**

## Predictions (primary cell; threshold 0.01 in J, as Phases 101 and 106)

For each realistic model M ∈ {M1, M2} separately (M1 by its path median):

- **P1, the goal matters more.** V̄_M ≥ 2 · V̄_M0 and V̄_M ≥ 0.01.
- **P2, the goal changes the design more often.** switches_M > switches_M0.
- **P3, more for quantum than classical.** V̄_M (quantum) > V̄_M (classical).

M1s is reported with the same three statistics as a stress case and does not enter the verdicts.

My own expectation, recorded before running: P1 roughly a coin flip for M1 and less likely for M2 (correlated
erasure on a ring with uniform marginals keeps much of the symmetry); P2 likely for M1; P3 uncertain.

## Bid decision rule (fixed now)

- **Pass:** P1 holds for M1 or M2. Then Test B (hardware target swap) is prepared and the "what to keep" vision stays
  a candidate centrepiece.
- **Fail:** P1 fails for both. The goal's value stays narrow under realistic noise at this size; the vision is not the
  centrepiece, and Test B is not run for the bid (it may still be run for the Top-Down paper).
- P2 and P3 are reported but do not change the decision.

## Validation before the grid is computed

- **V1.** With M0, J*(w), max SEP and the classical equivalents reproduce `../phase106/results.json` for every n,
  p̄, w within 1e-9.
- **V2.** Every noise distribution sums to 1 within 1e-12, and every per-qubit marginal equals its specified value
  (p_i for M0, M1, M1s; p̄ for M2) within 1e-12.
- **V3.** Independent path: for 300 random labelled codes at n = 6 (random graph, random ordered reference pair, valid
  only), P_A and P_B under M1 (path 0) and M2 at p̄ = 0.2 are recomputed region by region from fresh cut-ranks and
  the explicit noise formula; none may lie above the computed optimum for any w, and for codes that are in the
  enumerated set the values must match exactly.
- **V4.** SEP(n) ⊂ ALL(n): the sharing gain is never below −1e-12.

## Test B (hardware), outline only; parameters fixed in an addendum before any job

- The Phase 106 [[6,2]] code with unequal protection, each logical maximally entangled with a reference qubit
  (8 qubits). Keep a 4-qubit region that recovers A and not B; a control region that recovers B and not A; a region
  that recovers neither. Measure the reference–logical correlations needed to certify recovery (X, Y, Z logical
  correlations) for each region and target.
- Prediction: for the same region, the recovery fidelity of A clearly exceeds that of B, and the difference flips for
  the control region. First in simulation with the device noise model, then on IBM hardware within the user's 10-minute
  allowance after a small pilot job measures actual usage.

## Decision rules for running unattended

- D1. If V1–V4 fail: find and fix the bug, rerun all validations, log it. Predictions do not change.
- D2. If the n = 6 computation is projected to exceed 3 hours: reduce M1 to ten paths (logged as a deviation) before
  changing anything else.
- D3. Anything unforeseen: take the option harder for the goal to matter, log it as a deviation, run both where
  practical.
- D4. No existing file outside `phase107/` is edited by the run, and nothing is committed.

## Disclosure

Written by the same assistant that runs it. Seen before writing: Phase 106's full results (above), Phase 101's results,
and three external reviews of other bid directions. Nothing in this setting has been computed.
