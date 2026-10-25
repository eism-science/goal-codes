"""Phase 107 Test B: build the 20 circuits and run them in simulation (noise-free and FakeFez noise), or on hardware.

Criteria are fixed in TESTB_ADDENDUM.md. Usage:
    python testB_sim.py sim            noise-free and FakeFez-noise runs, writes testB_sim_results.json
"""
import json
import math
import os
import sys

import numpy as np
from qiskit import QuantumCircuit, transpile

HERE = os.path.dirname(os.path.abspath(__file__))
V = 8
SHOTS = 4000


def load_plan():
    with open(os.path.join(HERE, 'testB_plan.json')) as f:
        return json.load(f)


def prep(plan):
    qc = QuantumCircuit(V)
    qc.h(range(V))
    for u, v in plan['edges']:
        qc.cz(u, v)
    return qc


def settings(plan):
    """List of (name, pauli string over virtual qubits, ideal sign, inside_region)."""
    out = []
    for key, cor in plan['correlators'].items():
        for L in ('X', 'Y', 'Z'):
            c = cor[L]
            out.append((f'{key}/{L}', c['pauli'], 1 if c['sign'] > 0 else -1, c['inside_region']))
    return out


def measure_circuit(plan, pauli):
    qc = prep(plan)
    for v, L in enumerate(pauli):
        if L == 'X':
            qc.h(v)
        elif L == 'Y':
            qc.sdg(v)
            qc.h(v)
    qc.measure_all()
    return qc


def calib_circuits():
    c0 = QuantumCircuit(V)
    c0.measure_all()
    c1 = QuantumCircuit(V)
    c1.x(range(V))
    c1.measure_all()
    return [c0, c1]


def parity_expectation(counts, support):
    tot = sum(counts.values())
    s = 0
    for bits, n in counts.items():
        b = bits.replace(' ', '')[::-1]                   # b[v] = bit of virtual qubit v
        par = sum(int(b[v]) for v in support) % 2
        s += n * (1 - 2 * par)
    return s / tot


def analyse(counts_list, plan, sets):
    c0, c1 = counts_list[-2], counts_list[-1]
    e0 = [1 - parity_expectation(c0, [v]) for v in range(V)]   # 2*P(1|0)
    e1 = [1 + parity_expectation(c1, [v]) for v in range(V)]   # 2*P(0|1)
    shrink = [1 - e0[v] / 2 - e1[v] / 2 for v in range(V)]
    vals = {}
    for (name, pauli, sign, inside), counts in zip(sets, counts_list[:-2]):
        support = [v for v, L in enumerate(pauli) if L != 'I']
        raw = sign * parity_expectation(counts, support)
        corr = raw / np.prod([shrink[v] for v in support])
        vals[name] = {'raw': raw, 'corrected': float(corr), 'ideal': float(sign * sign if inside else 0),
                      'se': max(math.sqrt(max(1 - raw ** 2, 0) / SHOTS), 1e-9)}
    out = {'correlators': vals, 'readout_shrink': shrink, 'stats': {}}
    for key in plan['correlators']:
        for kind in ('raw', 'corrected'):
            cx, cy, cz = (vals[f'{key}/{L}'][kind] for L in 'XYZ')
            se = math.sqrt(vals[f'{key}/X']['se'] ** 2 + vals[f'{key}/Z']['se'] ** 2)
            out['stats'].setdefault(key, {})[kind] = {'E': cx + cz, 'E_se': se, 'W': (1 + cx + cy + cz) / 4}
    s = out['stats']
    tests = {}
    for kind in ('raw', 'corrected'):
        eb, ea_both, ea = s['main_B_not_A/B'][kind], s['both/A'][kind], s['main_B_not_A/A'][kind]
        d_se = math.sqrt(eb['E_se'] ** 2 + ea['E_se'] ** 2)
        tests[kind] = {
            'B1': {'value': eb['E'], 'z': (eb['E'] - 1) / eb['E_se'], 'holds': eb['E'] - 1 > 3 * eb['E_se']},
            'B2': {'value': ea_both['E'], 'z': (ea_both['E'] - 1) / ea_both['E_se'],
                   'holds': ea_both['E'] - 1 > 3 * ea_both['E_se']},
            'B3': {'value': eb['E'] - ea['E'], 'z': (eb['E'] - ea['E']) / d_se, 'holds': eb['E'] - ea['E'] > 3 * d_se},
        }
        tests[kind]['pass'] = all(tests[kind][k]['holds'] for k in ('B1', 'B2', 'B3'))
    out['tests'] = tests
    return out


def build(plan, backend):
    sets = settings(plan)
    circs = [measure_circuit(plan, p) for _, p, _, _ in sets] + calib_circuits()
    t0 = transpile(prep(plan), backend, optimization_level=3, seed_transpiler=plan['transpile_seed'])
    layout = [t0.layout.initial_layout[q] for q in prep(plan).qubits]
    tc = transpile(circs, backend, optimization_level=3, initial_layout=layout,
                   seed_transpiler=plan['transpile_seed'])
    n2 = [sum(1 for ins in c.data if ins.operation.num_qubits == 2) for c in tc]
    return sets, tc, layout, n2


def sim():
    from qiskit_aer import AerSimulator
    from qiskit_ibm_runtime.fake_provider import FakeFez
    plan = load_plan()
    backend = FakeFez()
    sets, tc, layout, n2 = build(plan, backend)
    print('layout', layout, 'two-qubit gates per circuit', sorted(set(n2)))
    res = {'layout': [int(x) for x in layout], 'two_qubit_gates': n2}
    for label, simulator in (('ideal', AerSimulator()), ('fakefez', AerSimulator.from_backend(backend))):
        job = simulator.run(tc, shots=SHOTS, seed_simulator=1072)
        counts = [job.result().get_counts(i) for i in range(len(tc))]
        out = analyse(counts, plan, sets)
        res[label] = out
        print(f'--- {label}')
        for key, st in out['stats'].items():
            print(f"  {key:22s} E raw {st['raw']['E']:+.3f} (se {st['raw']['E_se']:.3f})  corr {st['corrected']['E']:+.3f}"
                  f"   W raw {st['raw']['W']:.3f} corr {st['corrected']['W']:.3f}")
        for kind, t in out['tests'].items():
            print(f"  {kind}: " + '  '.join(f"{k} {t[k]['value']:+.3f} z={t[k]['z']:.1f} {'ok' if t[k]['holds'] else 'FAIL'}"
                                            for k in ('B1', 'B2', 'B3')) + f"  -> {'PASS' if t['pass'] else 'FAIL'}")
    with open(os.path.join(HERE, 'testB_sim_results.json'), 'w') as f:
        json.dump(res, f, indent=1, default=float)


if __name__ == '__main__':
    {'sim': sim}[sys.argv[1]]()
