"""Phase 111 Part B: build the region-local decoders and the hardware circuits, and verify them (V1-V3).

Usage:
    python p111_build.py            build against FakeFez (development; the frozen build uses the live backend)
    python p111_build.py --live     build against the live ibm_fez target (the frozen build)

Writes p111_plan.json, p111_candidates.json, p111_circuits.qpy and p111_target.pkl (the compilation target: gates,
coupling map, error rates, T1/T2; used by p111_verify.py and p111_power.py), then re-checks them with p111_verify.py and writes
p111_verify_report.json. Any failure of a primary pair, a baseline or the re-check stops the build with a nonzero exit
and writes only p111_plan_FAILED.json. Procedure: DESIGN.md, "Decoder synthesis" and
"Physical procedure and locality".
"""
import hashlib
import io
import json
import os
import sys

import numpy as np
from qiskit import ClassicalRegister, QuantumCircuit, qpy, transpile
from qiskit.converters import circuit_to_dag, dag_to_circuit
from qiskit.quantum_info import Clifford, Operator, Pauli, Statevector
from qiskit.transpiler import CouplingMap, PassManager
from qiskit.transpiler.passes import BasisTranslator, Optimize1qGatesDecomposition
from qiskit.circuit.equivalence_library import SessionEquivalenceLibrary

import p111_verify as PV

HERE = os.path.dirname(os.path.abspath(__file__))
PLAN_B = os.path.join(HERE, '..', 'phase107', 'testB_plan.json')
V = 8
REF = {'A': 0, 'B': 1}
SEEDS = range(10)
BASES = ('XX', 'YY', 'ZZ')
PAIRS = [('main_B_not_A', 'B'), ('main_B_not_A', 'A'), ('both', 'A'),
         ('both', 'B'), ('neither', 'A'), ('neither', 'B')]
PRIMARY = {('main_B_not_A', 'B'), ('main_B_not_A', 'A'), ('both', 'A')}
LAYOUTS = {}


# ---------------------------------------------------------------- symplectic helpers (local Paulis as letter lists)

def to_xz(letters):
    x = np.array([L in 'XY' for L in letters], dtype=np.uint8)
    z = np.array([L in 'ZY' for L in letters], dtype=np.uint8)
    return np.concatenate([x, z])


def omega(a, b):
    n = len(a) // 2
    return int((a[:n] @ b[n:] + a[n:] @ b[:n]) % 2)


def complete(fixed_pairs, single, n):
    """Symplectic basis of n qubits. fixed_pairs: [(x_image, z_image)] given; single: a z_image needing a partner
    (or None). Candidates in fixed order X_0, Z_0, X_1, Z_1, ... Returns the list of pairs, given ones first."""
    pairs = list(fixed_pairs)

    def reduce(v):
        v = v.copy()
        for a, b in pairs:
            v = (v + omega(v, b) * a + omega(v, a) * b) % 2
        return v

    basis = []
    for q in range(n):
        for k in (q, n + q):
            e = np.zeros(2 * n, dtype=np.uint8)
            e[k] = 1
            basis.append(e)
    if single is not None:
        z = reduce(single)
        x = next(reduce(e) for e in basis if omega(reduce(e), z) == 1)
        pairs.append((x, z))
    while len(pairs) < n:
        cands = [reduce(e) for e in basis]
        a = next(c for c in cands if c.any())
        b = next(c for c in cands if omega(a, c) == 1)
        pairs.append((a, b))
    return pairs


def decoder_clifford(P, Q, s, n=4):
    """Clifford D on n qubits with D P D^dag = +-X_s (if P given) and D Q D^dag = +-Z_s."""
    assert len(Q) == n and to_xz(Q).any(), ('Q must be a nonidentity Pauli on the region', Q)
    if P is not None:
        assert len(P) == n and to_xz(P).any() and omega(to_xz(P), to_xz(Q)) == 1, ('P, Q must anticommute', P, Q)
    fixed = [(to_xz(P), to_xz(Q))] if P is not None else []
    pairs = complete(fixed, to_xz(Q) if P is None else None, n)
    flat = [v for pair in pairs for v in pair]                 # must be a symplectic basis: omega = [[0,1],[1,0]]^n
    gram = np.array([[omega(a, b) for b in flat] for a in flat])
    assert np.array_equal(gram, np.kron(np.eye(n, dtype=int), np.array([[0, 1], [1, 0]]))), 'not a symplectic basis'
    order = [s] + [q for q in range(n) if q != s]
    tab = np.zeros((2 * n, 2 * n + 1), dtype=bool)
    for (xa, za), q in zip(pairs, order):
        tab[q, :2 * n] = xa                       # image of X_q
        tab[n + q, :2 * n] = za                   # image of Z_q
    C = Clifford(tab)                             # C X_q C^dag = pair images
    D = C.adjoint()
    for img, want in ((P, 'X'), (Q, 'Z')):
        if img is None:
            continue
        label = ''.join(reversed(img))
        out = Pauli(label).evolve(D, frame='s')
        target = ['I'] * n
        target[s] = want
        assert out.to_label().lstrip('-') == ''.join(reversed(target)), (img, out, s)
    return D


# ---------------------------------------------------------------- layouts: put the region on a connected path

N_ORDER = 2          # orderings of the region on each path kept (best edge overlap first)
LAYOUT_SEEDS = range(3)
KEEP = 5             # connected layouts kept per region (fewest preparation two-qubit gates)


def layout_candidates(plan, region, backend):
    """Initial layouts placing the region on each 4-qubit path of the device (orderings ranked by how many graph
    edges inside the region land on device edges), the other four qubits greedily next to their graph neighbours.
    Deterministic. Returns [(prep_ISA, initial_layout, seed)] whose region is connected after routing, best first."""
    import itertools
    import networkx as nx
    G = nx.Graph([tuple(e) for e in backend.coupling_map.get_edges()])
    dist = dict(nx.all_pairs_shortest_path_length(G, cutoff=6))
    gedges = [tuple(e) for e in plan['edges']]
    nbrs = {q: {v for e in gedges for v in e if q in e and v != q} for q in range(V)}
    paths = set()
    for u in G:
        for v in G[u]:
            for w in G[v]:
                if w == u:
                    continue
                for x in G[w]:
                    if x not in (u, v):
                        paths.add(min((u, v, w, x), (x, w, v, u)))
    others = [q for q in range(V) if q not in region]
    found = []
    for T in sorted(paths):
        scored = []
        for perm in itertools.permutations(region):
            pos = dict(zip(perm, T))
            score = sum(1 for a, b in gedges if a in pos and b in pos and G.has_edge(pos[a], pos[b]))
            scored.append((-score, perm))
        for _, perm in sorted(scored)[:N_ORDER]:
            pos = dict(zip(perm, T))
            for q in sorted(others, key=lambda q: (-len(nbrs[q] & set(region)), q)):
                placed = [pos[v] for v in nbrs[q] if v in pos] or list(T)
                free = [n for n in G if n not in pos.values()]
                pos[q] = min(free, key=lambda n: (sum(dist[n].get(m, 99) for m in placed), n))
            init = [pos[q] for q in range(V)]
            for seed in LAYOUT_SEEDS:
                tp = transpile(prep(plan), backend, optimization_level=3, seed_transpiler=seed,
                               initial_layout=init)
                fin = tp.layout.final_index_layout(filter_ancillas=True)
                phys = [fin[q] for q in region]
                if nx.is_connected(G.subgraph(phys)):
                    n2 = sum(1 for ins in tp.data if ins.operation.num_qubits == 2)
                    found.append((n2, T, perm, seed, tp, init))
    found.sort(key=lambda f: (f[0], f[1], f[2], f[3]))
    return [(f[4], f[5], f[3]) for f in found[:KEEP]], len(found)


# ---------------------------------------------------------------- circuits

def prep(plan):
    qc = QuantumCircuit(V)
    qc.h(range(V))
    for u, v in plan['edges']:
        qc.cz(u, v)
    return qc


def rotate(qc, q, basis_letter):
    if basis_letter == 'X':
        qc.h(q)
    elif basis_letter == 'Y':
        qc.sdg(q)
        qc.h(q)


def ideal_corr(plan, region, ref, D, s, corr=''):
    """<sigma_R sigma_s> for sigma in X, Y, Z after prep + decoder (+ Pauli correction on s), on the 8 logical qubits."""
    qc = prep(plan)
    qc.append(D.to_circuit(), region)
    for g in corr:
        getattr(qc, g.lower())(region[s])
    sv = Statevector(qc)
    out = {}
    for L in 'XYZ':
        lab = ['I'] * V
        lab[ref] = L
        lab[region[s]] = L
        out[L] = float(np.real(sv.expectation_value(Pauli(''.join(reversed(lab))))))
    return out


def correction(c, control):
    """Pauli on s making the ideal state |Phi+> (XX = ZZ = +1), or ZZ = +1 for a control."""
    g = ''
    if c['Z'] < 0:
        g += 'X'
    if not control and c['X'] < 0:
        g += 'Z'
    return g


ONEQ = None


def oneq_pm(backend):
    global ONEQ
    if ONEQ is None:
        basis = [n for n in backend.operation_names if n in ('rz', 'sx', 'x', 'cz')]
        ONEQ = PassManager([BasisTranslator(SessionEquivalenceLibrary, basis),
                            Optimize1qGatesDecomposition(basis=basis)])
    return ONEQ


def build_pair(plan, backend, region_key, target):
    region = plan['regions'][region_key]
    ref = REF[target]
    cor = plan['correlators'][f'{region_key}/{target}']
    local = lambda pauli: [pauli[q] for q in region]
    control = not cor['X']['inside_region']
    P = None if control else local(cor['X']['pauli'])
    Q = local(cor['Z']['pauli'])
    cmap = backend.coupling_map
    edges = set(map(tuple, cmap.get_edges()))
    cands = []
    if region_key not in LAYOUTS:
        LAYOUTS[region_key] = layout_candidates(plan, region, backend)
    layouts, n_found = LAYOUTS[region_key]
    for pseed, (tp, init, lseed) in enumerate(layouts):
        pos = tp.layout.final_index_layout(filter_ancillas=True)          # virtual -> physical after routing
        phys = [pos[q] for q in region]
        sub = [(i, j) for i in range(4) for j in range(4) if i != j and (phys[i], phys[j]) in edges]
        connected = CouplingMap(sub).is_connected() if sub and len({a for e in sub for a in e}) == 4 else False
        n2_prep = sum(1 for ins in tp.data if ins.operation.num_qubits == 2)
        if not connected:
            cands.append({'layout_rank': pseed, 'connected': False, 'phys_region': phys})
            continue
        for s in range(4):
            D = decoder_clifford(P, Q, s)
            c0 = ideal_corr(plan, region, ref, D, s)
            corr = correction(c0, control)
            for dseed in SEEDS:
                circs, ok = {}, True
                for basis in BASES:
                    dec = QuantumCircuit(4)
                    dec.append(D.to_circuit(), range(4))
                    for g in corr:
                        getattr(dec, g.lower())(s)
                    rotate(dec, s, basis[1])
                    td = transpile(dec, coupling_map=CouplingMap(sub), basis_gates=['cz', 'rz', 'sx', 'x'],
                                   optimization_level=3, seed_transpiler=dseed, initial_layout=[0, 1, 2, 3])
                    s_end = td.layout.final_index_layout()[s] if td.layout else s
                    full = tp.copy()
                    full.barrier()
                    full.compose(td, qubits=phys, inplace=True)
                    rotate(full, pos[ref], basis[0])
                    full = oneq_pm(backend).run(full)
                    creg = ClassicalRegister(2, 'meas')
                    full.add_register(creg)
                    full.measure(pos[ref], creg[0])
                    full.measure(phys[s_end], creg[1])
                    circs[basis] = full
                n2 = sum(1 for ins in circs['ZZ'].data if ins.operation.num_qubits == 2)
                depth = circs['ZZ'].depth()
                cands.append({'layout_rank': pseed, 'initial_layout': init, 'layout_seed': lseed, 'connected': True, 'phys_region': phys, 's': s,
                              'decoder_seed': dseed, 'correction': corr, 'ideal_logical': c0,
                              'two_qubit_gates': n2, 'two_qubit_prep': n2_prep, 'depth': depth,
                              'layout_final': [pos[q] for q in range(V)], 'ref_phys': pos[ref],
                              'out_phys': phys[s_end], '_circs': circs})
    return cands


def best_in(cands, rank):
    """Best candidate for one pair under layout `rank`: fewest total two-qubit gates (preparation + decoder), then
    depth, then output qubit, then decoder seed."""
    good = [c for c in cands if c.get('connected') and c['layout_rank'] == rank]
    return min(good, key=lambda c: (c['two_qubit_gates'], c['depth'], c['s'], c['decoder_seed'])) if good else None


def choose_region(cands_by_target):
    """One layout per region, shared by both targets (so a comparison between targets on a region uses the same
    preparation): the layout minimising the sum of the two targets' best total two-qubit gates, ties by rank."""
    ranks = sorted({c['layout_rank'] for cs in cands_by_target.values() for c in cs if c.get('connected')})
    scored = []
    for r in ranks:
        bests = {t: best_in(cs, r) for t, cs in cands_by_target.items()}
        if all(bests.values()):
            scored.append((sum(b['two_qubit_gates'] for b in bests.values()), r, bests))
    return min(scored, key=lambda x: (x[0], x[1]))[2] if scored else {t: None for t in cands_by_target}


# ---------------------------------------------------------------- calibration and baseline circuits

def aux_circuits(plan_out, backend):
    """For each distinct measured physical pair (R, s): four calibration circuits and a direct Bell-pair baseline in
    three bases. The baseline prepares the Bell pair on R and its neighbour on the shortest device path to s, then
    moves the second half to s with SWAPs through the path's qubits, which start and end in |0>; when R and s are
    adjacent no SWAP is needed. Returns (names, circuits, record)."""
    import networkx as nx
    G = nx.Graph([tuple(e) for e in backend.coupling_map.get_edges()])
    names, circs, rec = [], [], {}
    pairs = []
    for name, r in plan_out['pairs'].items():
        if r.get('status') != 'built':
            continue
        R, s_ = r['ref_phys'], r['out_phys']
        r['physical_pair'] = f'{R}-{s_}'
        if (R, s_) not in pairs:
            pairs.append((R, s_))
    n = backend.num_qubits
    for R, s_ in pairs:
        key = f'{R}-{s_}'
        for ab in ('00', '01', '10', '11'):
            qc = QuantumCircuit(n)
            if ab[0] == '1':
                qc.x(R)
            if ab[1] == '1':
                qc.x(s_)
            creg = ClassicalRegister(2, 'meas')
            qc.add_register(creg)
            qc.measure(R, creg[0])
            qc.measure(s_, creg[1])
            names.append(f'calib/{key}/{ab}')
            circs.append(qc)
        path = nx.shortest_path(G, R, s_)
        checks = {}
        for basis in BASES:
            qc = QuantumCircuit(n)
            qc.h(R)
            qc.cx(R, path[1])
            for a, b in zip(path[1:-1], path[2:]):
                qc.swap(a, b)
            qc.barrier()
            rotate(qc, R, basis[0])
            rotate(qc, s_, basis[1])
            qc = oneq_pm(backend).run(qc)
            creg = ClassicalRegister(2, 'meas')
            qc.add_register(creg)
            qc.measure(R, creg[0])
            qc.measure(s_, creg[1])
            val = PV.exact_parity(qc)
            want = {'XX': 1.0, 'YY': -1.0, 'ZZ': 1.0}[basis]
            bad = [ins.operation.name for ins in qc.data if ins.operation.name not in ('barrier', 'measure')
                   and not backend.target.instruction_supported(
                       ins.operation.name, tuple(qc.find_bit(q).index for q in ins.qubits))]
            checks[basis] = {'exact': val, 'ok': abs(val - want) < 1e-9 and not bad}
            names.append(f'baseline/{key}/{basis}')
            circs.append(qc)
        rec[key] = {'path': path, 'swaps': len(path) - 2,
                    'two_qubit_gates': sum(1 for ins in circs[-1].data if ins.operation.num_qubits == 2),
                    'checks': checks}
        print(f"pair {key}: path {path}, baseline 2q gates {rec[key]['two_qubit_gates']}, "
              f"exact {'PASS' if all(c['ok'] for c in checks.values()) else 'FAIL'}", flush=True)
    return names, circs, rec


def write_outputs(plan_out, circuits, names, all_cands):
    for c, n in zip(circuits, names):
        c.name = n                                # the QPY carries each circuit's name; the notebook checks the order
    buf = io.BytesIO()
    qpy.dump(circuits, buf)
    plan_out['circuit_names'] = names
    plan_out['qpy_sha256'] = hashlib.sha256(buf.getvalue()).hexdigest()
    with open(os.path.join(HERE, 'p111_circuits.qpy'), 'wb') as f:
        f.write(buf.getvalue())
    with open(os.path.join(HERE, 'p111_plan.json'), 'w') as f:
        json.dump(plan_out, f, indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
    with open(os.path.join(HERE, 'p111_candidates.json'), 'w') as f:
        json.dump(all_cands, f, default=lambda o: o.item() if hasattr(o, 'item') else str(o))


def get_backend(live):
    """FakeFez for development; with live=True the current ibm_fez target (token from IBM_QUANTUM_TOKEN, instance
    from IBM_QUANTUM_INSTANCE), whose calibration is saved as p111_backend_snapshot.json."""
    if not live:
        from qiskit_ibm_runtime.fake_provider import FakeFez
        return FakeFez(), None
    from qiskit_ibm_runtime import QiskitRuntimeService
    service = QiskitRuntimeService(channel='ibm_cloud', token=os.environ['IBM_QUANTUM_TOKEN'],
                                   instance=os.environ['IBM_QUANTUM_INSTANCE'])
    backend = service.backend('ibm_fez')
    props = backend.properties()
    if props is None:
        sys.exit('STOP: the live backend returned no calibration data; no build written')
    data = json.dumps(props.to_dict(), default=str, sort_keys=True).encode()
    with open(os.path.join(HERE, 'p111_backend_snapshot.json'), 'wb') as f:
        f.write(data)
    return backend, {'file': 'p111_backend_snapshot.json', 'sha256': hashlib.sha256(data).hexdigest(),
                     'last_update_date': str(props.last_update_date)}


def stop(msg, plan_out):
    """A failed build is written only as p111_plan_FAILED.json, never as a freezeable p111_plan.json."""
    with open(os.path.join(HERE, 'p111_plan_FAILED.json'), 'w') as f:
        json.dump(plan_out, f, indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
    sys.exit(f'STOP: {msg}. Part B is not run under this design (DESIGN.md, "Failure consequences"); '
             'details in p111_plan_FAILED.json')


def main():
    import qiskit, qiskit_aer, qiskit_ibm_runtime
    live = '--live' in sys.argv
    backend, snapshot = get_backend(live)
    with open(PLAN_B) as f:
        plan = json.load(f)
    plan_out = {'backend': backend.name, 'live_target': live, 'calibration_snapshot': snapshot,
                'target': PV.save_target(backend.target), 'pairs': {},
                'versions': {'python': sys.version.split()[0], 'qiskit': qiskit.__version__,
                             'qiskit_aer': qiskit_aer.__version__, 'qiskit_ibm_runtime': qiskit_ibm_runtime.__version__,
                             'numpy': np.__version__}}
    for f in ('p111_plan_FAILED.json', 'p111_verify_report.json', 'p111_plan.json', 'p111_circuits.qpy',
              'p111_candidates.json', 'p111_power.json'):          # no output of an earlier build survives
        if os.path.exists(os.path.join(HERE, f)):
            os.remove(os.path.join(HERE, f))
    all_cands, circuits, names = {}, [], []
    chosen = {}
    for key in dict.fromkeys(k for k, _ in PAIRS):
        targets = [t for k, t in PAIRS if k == key]
        cands = {t: build_pair(plan, backend, key, t) for t in targets}
        for t, cs in cands.items():
            all_cands[f'{key}/{t}'] = [{k: v for k, v in c.items() if k != '_circs'} for c in cs]
        for t, b in choose_region(cands).items():
            chosen[(key, t)] = b
    for key, target in PAIRS:
        name = f'{key}/{target}'
        primary = (key, target) in PRIMARY
        best = chosen[(key, target)]
        if best is None:
            plan_out['pairs'][name] = {'status': 'not run: no connected layout', 'primary': primary}
            print(f'{name}: NO CONNECTED LAYOUT', flush=True)
            continue
        rec = {k: v for k, v in best.items() if k != '_circs'}
        rec.update({'primary': primary, 'control': not plan['correlators'][name]['X']['inside_region']})
        checks, vals = {}, {}
        for basis, circ in best['_circs'].items():
            bad, vals[basis] = PV.check_pair(circ, rec, basis, plan, key, target, backend.target)
            checks[basis] = {'failures': bad, 'V2_exact': vals[basis]}
        F = (1 + vals['XX'] - vals['YY'] + vals['ZZ']) / 4
        allok = not any(c['failures'] for c in checks.values())
        rec.update({'status': 'built' if allok else 'not run: failed V1-V3', 'checks': checks, 'ideal_F_compiled': F})
        plan_out['pairs'][name] = rec
        print(f"{name}: layout rank {best['layout_rank']}, s={best['s']}, decoder seed {best['decoder_seed']}, "
              f"2q gates {best['two_qubit_gates']} (prep {best['two_qubit_prep']}), depth {best['depth']}, "
              f"F_ideal {F:.6f}, V1-V3 {'PASS' if allok else 'FAIL'}", flush=True)
        if allok:
            for basis, circ in best['_circs'].items():
                circuits.append(circ)
                names.append(f'{name}/{basis}')
    bad = [n for n, r in plan_out['pairs'].items() if r['primary'] and r['status'] != 'built']
    if bad:
        stop(f'primary pair(s) {bad} failed', plan_out)
    an, ac, rec = aux_circuits(plan_out, backend)
    plan_out['baselines'] = rec
    if not all(c['ok'] for r in rec.values() for c in r['checks'].values()):
        stop('a baseline circuit failed its exact check', plan_out)
    write_outputs(plan_out, circuits + ac, names + an, all_cands)
    report = PV.verify_all()                      # independent re-check of the files as written
    with open(os.path.join(HERE, 'p111_verify_report.json'), 'w') as f:
        json.dump(report, f, indent=1)
    if not report['pass']:
        os.replace(os.path.join(HERE, 'p111_plan.json'), os.path.join(HERE, 'p111_plan_FAILED.json'))
        sys.exit(f"STOP: verification of the written files failed: {report['failures']}")
    print('build written and verified:', len(circuits + ac), 'circuits')


if __name__ == '__main__':
    main()
