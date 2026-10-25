"""Phase 111: independent checks V1-V3 of the saved circuits (DESIGN.md, "Validity checks before any QPU job").

Usage:
    python p111_verify.py           check p111_circuits.qpy against p111_plan.json; writes p111_verify_report.json

Works on the QPY as saved, not on objects held by the builder. Checks:
- package: QPY hash and circuit names (each circuit's own name) match the plan, in order;
- V1 locality: after the single barrier, two-qubit gates only inside the region's physical positions, no SWAP, reset or
  gate on more than two qubits; on the reference only the specified basis rotation (checked as an operator); no other
  qubit touched; exactly two measurements, last on their qubits;
- V2 exact: noiseless <sigma sigma> of the measured pair equals the ideal value to 1e-9;
- V3 mapping: the state before the barrier is the code's graph state on the recorded 8-qubit physical layout (checked
  through its stabilizers, independently of the compiler's layout record), and clbit 0 / clbit 1 measure the recorded
  reference and output positions;
- calibration and baseline circuits: gates, measurement mapping, and the baseline's exact Bell correlators;
- Part A2: each archived Test B circuit is, in order, the setting its manifest name says: the state before measurement
  equals the graph state with that setting's basis rotations (or the calibration state) on its final layout, and
  clbit i measures logical qubit i's final position.
Instruction support is checked against the backend target saved at compilation (`p111_target.pkl`, whose SHA-256 the
plan records), so a live build is checked against the live target it was compiled for; the run notebook checks again
against the target on the day. Exits with status 1 on any failure.
"""
import ast
import hashlib
import json
import os
import pickle
import sys

import numpy as np
from qiskit import QuantumCircuit, qpy
from qiskit.converters import circuit_to_dag, dag_to_circuit
from qiskit.quantum_info import Operator, Pauli, Statevector

HERE = os.path.dirname(os.path.abspath(__file__))
PLAN_B = os.path.join(HERE, '..', 'phase107', 'testB_plan.json')
TESTB_QPY = os.path.join(HERE, '..', 'phase107', 'IBM run data', 'testB_transpiled.qpy')
TESTB_MANIFEST = os.path.join(HERE, '..', 'phase107', 'IBM run data', 'testB_manifest.json')
V = 8
REF = {'A': 0, 'B': 1}
BASES = ('XX', 'YY', 'ZZ')
GATES = {'cz', 'rz', 'sx', 'x'}
PRIMARY = {'main_B_not_A/B', 'main_B_not_A/A', 'both/A'}


def idx(circ, q):
    return circ.find_bit(q).index


# ---------------------------------------------------------------- the compilation target, saved with the build

TARGET_FILE = 'p111_target.pkl'


def save_target(target, folder=HERE):
    """Pickle the backend target used for compilation (gates, coupling map, error rates, durations, T1/T2)."""
    data = pickle.dumps(target)
    with open(os.path.join(folder, TARGET_FILE), 'wb') as f:
        f.write(data)
    return {'file': TARGET_FILE, 'sha256': hashlib.sha256(data).hexdigest()}


def load_target(plan, folder=HERE):
    """The saved compilation target, after checking its hash against the plan. Our own file: only load it from the
    frozen package."""
    data = open(os.path.join(folder, plan['target']['file']), 'rb').read()
    if hashlib.sha256(data).hexdigest() != plan['target']['sha256']:
        raise ValueError(f"{plan['target']['file']} differs from the plan's record")
    return pickle.loads(data)


def snapshot_backend(target, name):
    """A backend object around a saved target, for Aer's noise model; it cannot run jobs."""
    from qiskit.providers import BackendV2, Options

    class SnapshotBackend(BackendV2):
        def __init__(self):
            super().__init__(name=name)

        @property
        def target(self):
            return target

        @property
        def max_circuits(self):
            return None

        @classmethod
        def _default_options(cls):
            return Options()

        def run(self, *args, **kwargs):
            raise RuntimeError('snapshot backend: cannot run jobs')

    return SnapshotBackend()


def fingerprint(circ):
    """SHA-256 of a circuit's instruction list (gate, parameters to 1e-10, qubit and clbit indices), for comparing a
    submitted circuit with the frozen one independently of serialisation details."""
    rows = [(ins.operation.name, [round(float(x), 10) for x in ins.operation.params],
             [idx(circ, q) for q in ins.qubits], [idx(circ, c) for c in ins.clbits]) for ins in circ.data]
    return hashlib.sha256(json.dumps([circ.num_qubits, circ.num_clbits, rows]).encode()).hexdigest()


def rotation_op(letter):
    qc = QuantumCircuit(1)
    if letter == 'X':
        qc.h(0)
    elif letter == 'Y':
        qc.sdg(0)
        qc.h(0)
    return Operator(qc)


def measurements(circ):
    """clbit -> qubit for every measurement; fails on more than one per clbit."""
    meas = {}
    for ins in circ.data:
        if ins.operation.name == 'measure':
            c = idx(circ, ins.clbits[0])
            if c in meas:
                raise ValueError(f'clbit {c} measured twice')
            meas[c] = idx(circ, ins.qubits[0])
    return meas


def measures_last(circ):
    """True if nothing acts on a measured qubit after its measurement."""
    done = set()
    for ins in circ.data:
        qs = {idx(circ, q) for q in ins.qubits}
        if ins.operation.name == 'measure':
            done |= qs
        elif ins.operation.name != 'barrier' and qs & done:
            return False
    return True


def unsupported(circ, target):
    return [(ins.operation.name, tuple(idx(circ, q) for q in ins.qubits)) for ins in circ.data
            if ins.operation.name not in ('barrier', 'measure')
            and not target.instruction_supported(ins.operation.name, tuple(idx(circ, q) for q in ins.qubits))]


def check_locality(circ, region_phys, ref_phys, out_phys, ref_letter):
    """V1 and the measurement part of V3. Returns a list of failures (empty if it passes)."""
    bad = []
    bars = [i for i, ins in enumerate(circ.data) if ins.operation.name == 'barrier']
    if len(bars) != 1:
        return [f'{len(bars)} barriers, expected 1']
    region = set(region_phys)
    if ref_phys in region:
        bad.append('reference inside region')
    if out_phys not in region:
        bad.append(f'output {out_phys} not in region {sorted(region)}')
    for ins in circ.data[:bars[0]]:
        if ins.operation.name in ('measure', 'reset'):
            bad.append(f'{ins.operation.name} before the boundary')
    ref_gates = QuantumCircuit(1)
    for ins in circ.data[bars[0] + 1:]:
        name, qs = ins.operation.name, [idx(circ, q) for q in ins.qubits]
        if name == 'measure':
            continue
        if name not in GATES or len(qs) > 2:
            bad.append(f'{name} on {qs} after the boundary')
        elif qs == [ref_phys]:
            ref_gates.append(ins.operation, [0])
        elif not set(qs) <= region:
            bad.append(f'{name} on {qs} outside the region after the boundary')
    if not Operator(ref_gates).equiv(rotation_op(ref_letter)):
        bad.append(f'reference gates are not the {ref_letter} basis rotation')
    try:
        meas = measurements(circ)
    except ValueError as e:
        return bad + [str(e)]
    if meas != {0: ref_phys, 1: out_phys} or circ.num_clbits != 2:
        bad.append(f'measurement map {meas}, expected {{0: {ref_phys}, 1: {out_phys}}}')
    if not measures_last(circ):
        bad.append('operation after a measurement')
    return bad


def active(circ, upto=None):
    """Circuit without measurements and barriers (optionally only the first `upto` instructions) on its non-idle
    qubits, and the map physical index -> position in it."""
    body = circ.copy_empty_like()
    for ins in circ.data[:upto]:
        if ins.operation.name not in ('measure', 'barrier'):
            body.append(ins)
    body.remove_final_measurements()
    dag = circuit_to_dag(body)
    idle = [w for w in dag.idle_wires() if w in dag.qubits]
    dag.remove_qubits(*idle)
    keep = [q for q in body.qubits if q not in idle]
    return dag_to_circuit(dag), {idx(body, q): i for i, q in enumerate(keep)}


def expect(sv, where, paulis):
    """<prod of physical Paulis>; paulis: physical index -> letter. Qubits not in the circuit are |0>."""
    lab = ['I'] * sv.num_qubits
    for p, L in paulis.items():
        if p not in where:
            if L in 'XY':
                return 0.0
            continue
        lab[where[p]] = L
    return float(np.real(sv.expectation_value(Pauli(''.join(reversed(lab))))))


def exact_parity(circ):
    """V2: exact <Z Z> on the two measured qubits (basis rotations are in the circuit), on its active qubits."""
    meas = measurements(circ)
    small, where = active(circ)
    return expect(Statevector(small), where, {meas[0]: 'Z', meas[1]: 'Z'})


def check_layout(circ, layout, edges):
    """V3 layout: the state before the barrier has stabilizers X_v Z_N(v) = +1 on the physical positions layout[v],
    and any other qubit used for routing is back in |0>."""
    bar = next(i for i, ins in enumerate(circ.data) if ins.operation.name == 'barrier')
    small, where = active(circ, bar)
    if not set(layout) <= set(where):
        return [f'layout {sorted(layout)} not all active before the boundary {sorted(where)}']
    sv = Statevector(small)
    bad = [f'routing qubit {p} not back in |0> at the boundary' for p in where if p not in layout
           and abs(expect(sv, where, {p: 'Z'}) - 1) > 1e-9]
    for v in range(V):
        stab = {layout[v]: 'X'}
        for a, b in edges:
            if v in (a, b):
                stab[layout[b if a == v else a]] = 'Z'
        val = expect(sv, where, stab)
        if abs(val - 1) > 1e-9:
            bad.append(f'stabilizer of vertex {v} on the layout = {val:.6f}')
    return bad


def ideal(control, basis):
    return ({'XX': 0.0, 'YY': 0.0, 'ZZ': 1.0} if control else {'XX': 1.0, 'YY': -1.0, 'ZZ': 1.0})[basis]


def check_pair(circ, rec, basis, plan_b, region_key, target, tgt):
    layout = rec['layout_final']
    region = plan_b['regions'][region_key]
    bad = []
    if [layout[q] for q in region] != rec['phys_region']:
        bad.append('recorded region positions differ from the layout')
    if layout[REF[target]] != rec['ref_phys']:
        bad.append('recorded reference position differs from the layout')
    val = float('nan')
    try:
        bad += check_locality(circ, rec['phys_region'], rec['ref_phys'], rec['out_phys'], basis[0])
        bad += check_layout(circ, layout, plan_b['edges'])
        bad += [f'unsupported {u}' for u in unsupported(circ, tgt)]
        val = exact_parity(circ)
        if abs(val - ideal(rec['control'], basis)) > 1e-9:
            bad.append(f'V2 exact {val:.6f}, ideal {ideal(rec["control"], basis)}')
    except Exception as e:                        # a malformed circuit fails; it never stops the report
        bad.append(f'check could not run: {e!r}')
    return bad, val


def check_calib(circ, R, s, ab, tgt):
    bad = [f'unsupported {u}' for u in unsupported(circ, tgt)]
    xs = sorted(idx(circ, ins.qubits[0]) for ins in circ.data if ins.operation.name == 'x')
    other = [ins.operation.name for ins in circ.data if ins.operation.name not in ('x', 'measure', 'barrier')]
    want = sorted([R] * (ab[0] == '1') + [s] * (ab[1] == '1'))
    if xs != want or other:
        bad.append(f'preparation x on {xs}, other gates {other}; expected x on {want}')
    if measurements(circ) != {0: R, 1: s} or circ.num_clbits != 2:
        bad.append(f'measurement map {measurements(circ)}')
    return bad


def check_baseline(circ, R, s, basis, path, tgt):
    bad = [f'unsupported {u}' for u in unsupported(circ, tgt)]
    used = {idx(circ, q) for ins in circ.data if ins.operation.name not in ('barrier',) for q in ins.qubits}
    if not used <= set(path):
        bad.append(f'qubits {sorted(used - set(path))} outside the recorded path')
    if measurements(circ) != {0: R, 1: s} or circ.num_clbits != 2 or not measures_last(circ):
        bad.append(f'measurement map {measurements(circ)}')
    val = exact_parity(circ)
    want = {'XX': 1.0, 'YY': -1.0, 'ZZ': 1.0}[basis]
    if abs(val - want) > 1e-9:
        bad.append(f'exact {val:.6f}, ideal {want}')
    return bad


def restricted(circ, extra=()):
    """The circuit's gates (no measurements or barriers) on the qubits they use plus `extra`; and the map physical
    index -> position."""
    used = sorted({idx(circ, q) for ins in circ.data if ins.operation.name not in ('measure', 'barrier')
                   for q in ins.qubits} | set(extra))
    where = {p: i for i, p in enumerate(used)}
    qc = QuantumCircuit(len(used))
    for ins in circ.data:
        if ins.operation.name not in ('measure', 'barrier'):
            qc.append(ins.operation, [where[idx(circ, q)] for q in ins.qubits])
    return qc, where


def testB_names(man):
    names = man['circuit_names']
    return ast.literal_eval(names) if isinstance(names, str) else names


def testB_expected(plan_b, name, layout, where, n):
    """The ideal pre-measurement state of Test B circuit `name` (as built in phase107/testB_sim.py) on the active
    qubits of the compiled circuit: logical v at position where[layout[v]], other qubits in |0>."""
    qc = QuantumCircuit(n)
    pos = [where[layout[v]] for v in range(V)]
    if name == 'calib/1':
        qc.x(pos)
    elif name != 'calib/0':
        key, L = name.rsplit('/', 1)
        qc.h(pos)
        for u, v in plan_b['edges']:
            qc.cz(pos[u], pos[v])
        for v, P in enumerate(plan_b['correlators'][key][L]['pauli']):
            if P == 'X':
                qc.h(pos[v])
            elif P == 'Y':
                qc.sdg(pos[v])
                qc.h(pos[v])
    return Statevector(qc)


def check_A2(qpy_path, manifest_path, plan_b):
    """The archived Test B circuits: hash as in their manifest; circuit i is the setting named i-th in the manifest
    (its pre-measurement state checked exactly); clbit v measures logical qubit v's final position."""
    man = json.load(open(manifest_path))
    names = testB_names(man)
    data = open(qpy_path, 'rb').read()
    bad = [] if hashlib.sha256(data).hexdigest() == man['qpy_sha256'] else ['testB QPY hash differs from its manifest']
    circs = qpy.load(open(qpy_path, 'rb'))
    if len(circs) != len(names):
        return bad + [f'{len(circs)} circuits, {len(names)} names'], len(circs), []
    mapping = []
    for i, (c, name) in enumerate(zip(circs, names)):
        fin = c.layout.final_index_layout(filter_ancillas=True)
        if measurements(c) != {v: fin[v] for v in range(V)}:
            bad.append(f'Test B circuit {i} ({name}): measurement map differs from its final layout')
            continue
        small, where = restricted(c, [fin[v] for v in range(V)])
        if not Statevector(small).equiv(testB_expected(plan_b, name, fin, where, small.num_qubits)):
            bad.append(f'Test B circuit {i}: not the setting {name}')
        mapping.append({'index': i, 'name': name, 'physical': [fin[v] for v in range(V)]})
    return bad, len(circs), mapping


def guarded(check, *args):
    try:
        return check(*args)
    except Exception as e:
        return [f'check could not run: {e!r}']


def verify_all(plan_path=None, qpy_path=None, plan_b_path=PLAN_B):
    plan_path = plan_path or os.path.join(HERE, 'p111_plan.json')
    qpy_path = qpy_path or os.path.join(HERE, 'p111_circuits.qpy')
    plan = json.load(open(plan_path))
    plan_b = json.load(open(plan_b_path))
    raw = open(qpy_path, 'rb').read()
    circs = qpy.load(open(qpy_path, 'rb'))
    names = plan['circuit_names']
    report = {'qpy_sha256': hashlib.sha256(raw).hexdigest(),
              'plan_sha256': hashlib.sha256(open(plan_path, 'rb').read()).hexdigest(),
              'target_sha256': plan['target']['sha256'], 'live_target': plan['live_target'], 'failures': {},
              'pairs': {}}
    fail = report['failures']
    tgt = load_target(plan, os.path.dirname(os.path.abspath(plan_path)))
    report['fingerprints'] = {n: fingerprint(c) for n, c in zip(names, circs)}
    if report['qpy_sha256'] != plan['qpy_sha256']:
        fail['package'] = ['QPY hash differs from the plan']
    if [c.name for c in circs] != names:
        fail.setdefault('package', []).append('circuit names in the QPY differ from the plan, or are out of order')
    by_name = dict(zip(names, circs))
    for name, rec in plan['pairs'].items():
        if rec.get('status') != 'built':
            if name in PRIMARY:
                fail[name] = [f'primary pair not built: {rec.get("status")}']
            continue
        region_key, target = name.split('/')
        vals = {}
        for basis in BASES:
            c = by_name.get(f'{name}/{basis}')
            if c is None:
                fail.setdefault(name, []).append(f'{basis} circuit missing')
                continue
            bad, vals[basis] = check_pair(c, rec, basis, plan_b, region_key, target, tgt)
            if bad:
                fail[f'{name}/{basis}'] = bad
        if len(vals) == 3:
            report['pairs'][name] = {'exact': vals, 'F': (1 + vals['XX'] - vals['YY'] + vals['ZZ']) / 4}
    for key, b in plan.get('baselines', {}).items():
        R, s = map(int, key.split('-'))
        for ab in ('00', '01', '10', '11'):
            c = by_name.get(f'calib/{key}/{ab}')
            bad = ['missing'] if c is None else guarded(check_calib, c, R, s, ab, tgt)
            if bad:
                fail[f'calib/{key}/{ab}'] = bad
        for basis in BASES:
            c = by_name.get(f'baseline/{key}/{basis}')
            bad = ['missing'] if c is None else guarded(check_baseline, c, R, s, basis, b['path'], tgt)
            if bad:
                fail[f'baseline/{key}/{basis}'] = bad
    used = {r['physical_pair'] for r in plan['pairs'].values() if r.get('status') == 'built'}
    if used != set(plan.get('baselines', {})):
        fail['aux'] = [f'calibrated pairs {sorted(plan.get("baselines", {}))} differ from measured pairs {sorted(used)}']
    if os.path.exists(TESTB_QPY):
        bad, n, mapping = check_A2(TESTB_QPY, TESTB_MANIFEST, plan_b)
        report['A2_circuits'] = n
        report['A2_mapping'] = mapping
        report['A2_fingerprints'] = [fingerprint(c) for c in qpy.load(open(TESTB_QPY, 'rb'))]
        if bad:
            fail['A2'] = bad
    else:
        fail['A2'] = [f'{TESTB_QPY} not found']
    report['n_circuits'] = len(circs)
    report['pass'] = not fail
    return report


def main():
    report = verify_all()
    with open(os.path.join(HERE, 'p111_verify_report.json'), 'w') as f:
        json.dump(report, f, indent=1)
    for k, v in report['failures'].items():
        print(f'FAIL {k}: {v}')
    for k, r in report['pairs'].items():
        print(f"{k:16s} exact F {r['F']:.6f}")
    print(f"{report['n_circuits']} circuits, Test B archive {report.get('A2_circuits')} circuits:",
          'ALL CHECKS PASS' if report['pass'] else 'CHECKS FAILED')
    sys.exit(0 if report['pass'] else 1)


if __name__ == '__main__':
    main()
