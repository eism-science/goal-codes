# Phase 106 (design, frozen before any run): does sharing still pay when the goal is lopsided?

Status: written 2026-10-23 by Claude during an unattended overnight session at the user's request, before any
experiment code existed. Frozen by the SHA-256 recorded in `RUNLOG.md`. Anything changed after freezing is listed in
`RUNLOG.md` as a deviation, with its reason.

## Why

Phase 101 (run 2026-10-12, `../phase101/RESULTS.md`) found that the goal did not change the best design. A review
on 2026-10-23 found the null was largely built in: every design had to use all 20 qubits for 6 logicals with no
cost for size, the three goals tested were all near-uniform (another design only wins when one logical carries about
0.89 of the weight), Q was decided by a tie-break between labels, and the quantum separate-code baseline was
annealed while the classical one was exact, which favours a P4 pass.

Once qubits can be traded between targets, it is close to certain that the best design changes with the goal: you
spend more on what you value more. That is a manipulation check, not a finding. The question that can genuinely fail
is the one Phase 101's P4 raised: **encoding two logical qubits jointly beats protecting each separately at a balanced
goal; does that advantage survive when the goal is lopsided?** And is the answer different from the classical case,
where copying is allowed?

## Setting

- Two logical qubits, A (the valued one) and B. Goal weight w on A and 1 − w on B, with w ≥ 0.5 (by symmetry).
- n physical qubits, each erased independently with probability p. No blocks.
- P_X = probability that logical X is recoverable from the qubits that survive. Logical X is recoverable from a
  region S when I(R_X : S) = 2 in the purified code.
- Objective J(w) = w·P_A + (1 − w)·P_B. With a cost per qubit, J_λ = J − λ·n.
- Every quantity is exact: profiles are integer counts c_X[s] of surviving sets of size s that recover X, and
  P_X(p) = Σ_s c_X[s]·(1 − p)^s·p^(n − s). No sampling and no search in the confirmatory part.

## Families (all exhaustive)

- **Quantum, all codes, ALL(n):** every [[n,2]] stabilizer code. Every stabilizer state is local-Clifford equivalent
  to a graph state, and recoverability (a mutual information) does not change under local unitaries, so the purified
  codes are, without loss, graph states on n + 2 vertices with two labelled reference vertices (A and B) and
  S(R) = 2. Erasure is the same on every physical qubit, so graphs are needed only up to isomorphism. For n + 2 ≤ 7,
  the networkx graph atlas lists every graph up to isomorphism. For n = 6 (8 vertices), every 8-vertex graph is a
  7-vertex atlas graph plus one vertex with one of the 128 neighbourhoods (delete any vertex and relabel), so the 1,044
  × 128 extensions cover every case. Every ordered pair of vertices is tried as (A, B). Entropies are cut-ranks of the
  adjacency matrix over GF(2).
- **Quantum, separate, SEP(n):** a [[m_A,1]] code for A and a [[m_B,1]] code for B, m_A + m_B = n, m_A, m_B ≥ 1. Each
  single-logical code is exhaustive in the same way (graphs on m + 1 ≤ 7 vertices). SEP(n) ⊂ ALL(n), so the gain
  below is never negative.
- **Classical twin:** two bits. Each of the n coordinates stores any non-constant Boolean function of (x_A, x_B) (7
  types up to negation, linear and non-linear). ALL_cl(n): every multiset of types. SEP_cl(n): copies of x_A and copies
  of x_B only. A bit is recoverable from the surviving coordinates when they determine it. Copying is allowed here and
  forbidden in the quantum case; that is the contrast.

Sharing gain: **G_n(w, p) = max_ALL(n) J − max_SEP(n) J**, at matched n, so λ plays no part in it. The classical gain
G_cl is defined the same way.

## Grid

- w ∈ {0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.99}.
- p ∈ {0.1, 0.2, 0.3}; primary p = 0.2.
- n ∈ {2, …, 6}; primary n = 6 (the largest size that is exhaustive).
- λ ∈ {0, 0.01, 0.02, 0.05} for the descriptive cost-per-qubit analysis (best n from 2 to 6 per w).

## Predictions (primary cell n = 6, p = 0.2; threshold 0.01 in J, as in Phase 101)

- **P1, manipulation check.** The set of optimal ALL(6) profiles changes at least once across the w grid. Expected:
  holds. If it fails, this setting cannot show goal-dependence, and that is reported first.
- **P2, sharing at a balanced goal.** G_6(0.5) > 0.01. Expected: holds (Phase 101's P4 direction).
- **P3, the main question, sharing under a lopsided goal.** G_6(0.9) > 0.01. **Prediction: fails.** The gain shrinks
  below 0.01 by w = 0.9, because a lopsided goal is best served by putting nearly everything into one single-logical
  code. This is a genuine guess, stated so that either outcome is informative. If P3 holds, sharing pays even when one
  target dominates, which is the stronger result for the bid.
- **P4, shape.** G_6(w) is non-increasing in w across the grid. Expected: holds.
- **P5, quantum against classical.** G_6(0.5) > G_cl,6(0.5). Expected: holds (Phase 101 found a classical shared gain
  of 0.0048 or 0 where the quantum gain was about 0.02).

Reported for every n and p, descriptively: G(w); the optimal profiles and where they change (breakpoints); the regret
of keeping the w = 0.5 optimum at each w; under J_λ, the best n and whether the best design is shared, per w.

Comparison rule: designs are compared by their profile pair (c_A, c_B), never by labels. "The optimum" means the
whole tie class of optimal profiles. Equality tolerance in J is 1e-12.

## Validation before the grid is computed

- **V1.** Graph-state recoverability equals `phase101/p101_core.Code.recoverable` (three external reviews; matched
  I(R:A) = 2 on 17,184 cases on 2026-10-23) on at least 2,000 random graphs with n from 2 to 6, every region.
- **V2.** Known profiles: a bare qubit recovers only from itself; [[5,1,3]] recovers from every 3 of 5; [[4,2,2]]
  recovers each logical from every 3 of 4 and from no 2; a Bell pair between two physical qubits plus nothing else is
  not a valid code (S(R) check).
- **V3.** Enumeration size: 1,044 seven-vertex atlas graphs, 133,632 extensions for n = 6.
- **V4.** Classical: copies of a bit give P = 1 − p^m; the parity coordinate alone recovers neither bit.
- **V5.** The best single-logical profiles never exceed the no-cloning bound (for each complementary pair of regions,
  at most one recovers the logical).

## Decision rules for running unattended (fixed now)

- D1. If V1–V5 fail: find and fix the bug, rerun all validations, log it. Predictions do not change.
- D2. If the n = 6 enumeration is projected to take over 3 hours: the primary cell becomes n = 5 (exhaustive), n = 6
  becomes exploratory, and this is logged as a deviation.
- D3. Anything unforeseen: take the option that is harder for sharing to pass, log it as a deviation, and run both
  options where practical.
- D4. Optional, only if time remains and labelled exploratory: n = 7 by local search over graphs, seeded with the
  n = 6 optima, with the single-logical [[7,1]] side searched the same way.
- D5. No existing file outside `phase106/` is edited, and nothing is committed.

## Disclosure

The predictions were written by the same assistant that will run the experiment. Facts seen before writing them:
Phase 101's results and its review (above), and Phases 98, 100 and 105 (classical). Nothing from this setting has been
computed.
