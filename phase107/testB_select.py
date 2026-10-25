"""Phase 107 Test B, step 1: choose the code wiring and the regions (no noise, no hardware).

The code is Phase 106's shared [[6,2]] optimum at n = 6, p = 0.1 (8-vertex graph state, references A = vertex 0,
B = vertex 1). Recoverability is unchanged by local unitaries, so any graph in its local-complementation orbit is the
same code up to local Cliffords; we pick the orbit member that compiles to the fewest two-qubit gates on ibm_fez's
heavy-hex layout, then choose the three regions and the stabilizer correlators to measure.
"""
import itertools
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'phase106'))
import p106  # noqa: E402

V, REF_A, REF_B = 8, 0, 1
PHYS = list(range(2, 8))


def local_complement(adj, v):
    adj = list(adj)
    nb = [u for u in range(V) if (adj[v] >> u) & 1]
    for i, u in enumerate(nb):
        for w in nb[i + 1:]:
            adj[u] ^= 1 << w
            adj[w] ^= 1 << u
    return tuple(adj)


def orbit(adj):
    seen = {tuple(adj)}
    frontier = [tuple(adj)]
    while frontier:
        nxt = []
        for g in frontier:
            for v in range(V):
                h = local_complement(g, v)
                if h not in seen:
                    seen.add(h)
                    nxt.append(h)
        frontier = nxt
    return list(seen)


def edges(adj):
    return [(u, v) for u in range(V) for v in range(u + 1, V) if (adj[u] >> v) & 1]


def mi(cr, X, Y):
    return int(cr[X]) + int(cr[Y]) - int(cr[X | Y])


def stabilizers(adj):
    """All 256 stabilizer elements as (sign, pauli string over V qubits) with letters I/X/Y/Z."""
    out = []
    for sub in range(1 << V):
        x = z = 0
        for v in range(V):
            if (sub >> v) & 1:
                x |= 1 << v
                z ^= adj[v]
        s = ''.join('Y' if (x >> v) & 1 and (z >> v) & 1 else 'X' if (x >> v) & 1 else 'Z' if (z >> v) & 1 else 'I'
                    for v in range(V))
        out.append((graph_expectation(adj, s), s))
    return out


def graph_state(adj):
    from qiskit.quantum_info import Statevector
    from qiskit import QuantumCircuit
    qc = QuantumCircuit(V)
    qc.h(range(V))
    for u, v in edges(adj):
        qc.cz(u, v)
    return Statevector(qc)


_SV = {}


def graph_expectation(adj, s):
    """<psi|P|psi> for Pauli string s (s[v] acts on qubit v), computed on the statevector."""
    from qiskit.quantum_info import Pauli
    key = tuple(adj)
    if key not in _SV:
        _SV[key] = graph_state(adj)
    return float(np.real(_SV[key].expectation_value(Pauli(s[::-1]))))


def correlators(stabs, ref, region):
    """For each Pauli letter on the reference: the lowest-weight stabilizer element whose reference letter is that
    letter. Returns letter -> (sign, string, contained_in_region)."""
    allowed = set(region) | {ref}
    other = REF_B if ref == REF_A else REF_A
    best = {}
    for sign, s in stabs:
        L = s[ref]
        if L == 'I' or s[other] != 'I':
            continue
        supp = [v for v in range(V) if s[v] != 'I']
        inside = set(supp) <= allowed
        key = (not inside, len(supp))
        if L not in best or key < best[L][0]:
            best[L] = (key, sign, s, inside)
    out = {}
    for L, b in best.items():
        s = b[2]
        trunc = ''.join(s[v] if v in allowed else 'I' for v in range(V))
        out[L] = (b[1], trunc, b[3])
    return out


def main():
    with open(os.path.join(HERE, '..', 'phase106', 'results.json')) as f:
        ex = json.load(f)['cells']['n=6,p=0.1']['rows'][0]['opt_example']
    adj0 = ex['adj']
    assert ex['A'] == REF_A and ex['B'] == REF_B
    orb = orbit(adj0)
    print('orbit size', len(orb))

    from qiskit import QuantumCircuit, transpile
    from qiskit_ibm_runtime.fake_provider import FakeFez
    backend = FakeFez()

    def circuit(adj):
        qc = QuantumCircuit(V)
        qc.h(range(V))
        for u, v in edges(adj):
            qc.cz(u, v)
        return qc

    cands = sorted(orb, key=lambda g: len(edges(g)))[:60]
    scored = []
    for g in cands:
        best = None
        for seed in range(5):
            t = transpile(circuit(g), backend, optimization_level=3, seed_transpiler=seed)
            n2 = sum(1 for ins in t.data if ins.operation.num_qubits == 2)
            if best is None or n2 < best[0]:
                best = (n2, t.depth(), seed)
        scored.append((best[0], best[1], len(edges(g)), g, best[2]))
    scored.sort(key=lambda x: (x[0], x[1]))
    n2, depth, ne, adj, seed = scored[0]
    print(f'chosen wiring: {ne} edges, {n2} two-qubit gates after routing, depth {depth} (original graph: '
          f'{len(edges(adj0))} edges)')

    cr = p106.cutranks(np.array(adj, dtype=np.int64), V)
    stabs = stabilizers(adj)
    rows = []
    for R in itertools.combinations(PHYS, 4):
        m = sum(1 << v for v in R)
        rows.append({'region': list(R), 'I_A': mi(cr, 1 << REF_A, m), 'I_B': mi(cr, 1 << REF_B, m)})
    # In this code B's recoverable 4-regions contain A's, so there is no "A but not B" region. Regions used:
    # main = B recovers fully (I = 2), A only classically (I = 1); both = both recover; neither = both I = 1.
    kinds = {'main_B_not_A': (1, 2), 'both': (2, 2), 'neither': (1, 1)}
    print('4-qubit regions by (I_A, I_B):', {k: sum((r['I_A'], r['I_B']) == v for r in rows) for k, v in kinds.items()})
    chosen = {}
    for name, v in kinds.items():
        match = [r['region'] for r in rows if (r['I_A'], r['I_B']) == v]
        if match:
            chosen[name] = match[0]
    plan = {'adj': list(adj), 'edges': edges(adj), 'transpile_seed': seed, 'two_qubit_gates': n2, 'depth': depth,
            'regions': chosen, 'region_table': rows, 'correlators': {}}
    for name, R in chosen.items():
        for ref, lab in ((REF_A, 'A'), (REF_B, 'B')):
            c = correlators(stabs, ref, R)
            plan['correlators'][f'{name}/{lab}'] = {L: {'sign': s, 'pauli': p, 'inside_region': ins}
                                                     for L, (s, p, ins) in c.items()}
            print(name, lab, {L: (s, p, ins) for L, (s, p, ins) in c.items()})
    with open(os.path.join(HERE, 'testB_plan.json'), 'w') as f:
        json.dump(plan, f, indent=1)


if __name__ == '__main__':
    main()
