"""Negative tests of p111_verify: each deliberately broken circuit must fail."""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p111_verify as PV
from qiskit import qpy
from qiskit.circuit import CircuitInstruction
from qiskit.circuit.library import CZGate, XGate, RZGate, Measure, Reset
from qiskit_ibm_runtime.fake_provider import FakeFez

H = os.path.dirname(os.path.abspath(__file__))
plan = json.load(open(os.path.join(H, 'p111_plan.json')))
pb = json.load(open(PV.PLAN_B))
circs = dict(zip(plan['circuit_names'], qpy.load(open(os.path.join(H, 'p111_circuits.qpy'), 'rb'))))
tgt = FakeFez().target
name, basis = 'main_B_not_A/B', 'XX'
rec = plan['pairs'][name]
base = circs[f'{name}/{basis}']
region, ref, out = rec['phys_region'], rec['ref_phys'], rec['out_phys']
outside = next(q for q in rec['layout_final'] if q not in region and q != ref)
other_s = next(q for q in region if q != out)
bad_ok, ok = PV.check_pair(base, rec, basis, pb, 'main_B_not_A', 'B', tgt)
assert not bad_ok, bad_ok
print('unbroken circuit passes')


def mutate(fn):
    c = base.copy()
    fn(c)
    return PV.check_pair(c, rec, basis, pb, 'main_B_not_A', 'B', tgt)[0]


def last_meas(c, clbit):
    return next(i for i, ins in enumerate(c.data) if ins.operation.name == 'measure'
                and c.find_bit(ins.clbits[0]).index == clbit)


def redirect(c, clbit, q):
    i = last_meas(c, clbit)
    c.data[i] = CircuitInstruction(Measure(), [c.qubits[q]], c.data[i].clbits)


bar = lambda c: next(i for i, ins in enumerate(c.data) if ins.operation.name == 'barrier')
cases = {
    'output measured on another region qubit': lambda c: redirect(c, 1, other_s),
    'reference measured on an outside qubit': lambda c: redirect(c, 0, outside),
    'gate on an outside qubit after boundary': lambda c: c.data.insert(bar(c) + 1, CircuitInstruction(XGate(), [c.qubits[outside]])),
    'cz from region to outside after boundary': lambda c: c.data.insert(bar(c) + 1, CircuitInstruction(CZGate(), [c.qubits[region[0]], c.qubits[outside]])),
    'extra rz on the reference': lambda c: c.data.insert(bar(c) + 1, CircuitInstruction(RZGate(0.3), [c.qubits[ref]])),
    'reset in the region after boundary': lambda c: c.data.insert(bar(c) + 1, CircuitInstruction(Reset(), [c.qubits[region[0]]])),
    'extra measurement of an outside qubit': lambda c: c.data.insert(bar(c) + 1, CircuitInstruction(Measure(), [c.qubits[outside]], [c.clbits[1]])),
    'gate on output after its measurement': lambda c: c.append(XGate(), [out]),
    'preparation altered (x on a code qubit)': lambda c: c.data.insert(bar(c), CircuitInstruction(XGate(), [c.qubits[rec['layout_final'][5]]])),
}
fails = 0
for k, fn in cases.items():
    bad = mutate(fn)
    print(('caught ' if bad else 'MISSED ') + k + ':', bad[:2])
    fails += not bad
print('ALL BROKEN CIRCUITS CAUGHT' if not fails else f'{fails} MISSED')
