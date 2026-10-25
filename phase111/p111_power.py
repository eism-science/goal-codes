"""Phase 111 Part C step 3: power check. Noisy density-matrix simulation of the built Part B circuits; expected
fidelities and expected z of D1-D3 (DESIGN.md: raise shots, up to 16,000, until every expected z is at least 5). Does
not change decoders. Its shot count is the one the run notebook uses.

The noise model is built from the backend target saved with the build (p111_target.pkl, hash checked against the
plan): for the frozen build, the live ibm_fez calibration at compilation time, so the circuits and their power check
use the same calibration; for the development build, FakeFez. Exits with status 1, writing no shot count, if no
allowed shot count reaches z = 5: then Part B is not run without a design amendment."""
import json
import math
import os
import sys

from qiskit import qpy
from qiskit_aer import AerSimulator

import p111_verify as PV

HERE = os.path.dirname(os.path.abspath(__file__))
SIM_SHOTS = 200_000
SIGNS = {'XX': 1, 'YY': -1, 'ZZ': 1}


def main():
    plan = json.load(open(os.path.join(HERE, 'p111_plan.json')))
    circs = qpy.load(open(os.path.join(HERE, 'p111_circuits.qpy'), 'rb'))
    names = plan['circuit_names']
    target = PV.load_target(plan)
    backend = PV.snapshot_backend(target, plan['backend'])
    snap = plan.get('calibration_snapshot') or {}
    source = {'model': ('live ibm_fez calibration at compilation' if plan['live_target']
                        else 'FakeFez snapshot (development reference, not current device conditions)'),
              'target_sha256': plan['target']['sha256'], 'calibration_snapshot_sha256': snap.get('sha256'),
              'last_update_date': snap.get('last_update_date')}
    sim = AerSimulator.from_backend(backend, method='density_matrix')
    res = sim.run(circs, shots=SIM_SHOTS, seed_simulator=111).result()
    c = {}
    for i, n in enumerate(names):
        if n.count('/') != 2 or n.startswith(('calib', 'baseline')):
            continue
        cnt = res.get_counts(i)
        tot = sum(cnt.values())
        even = sum(v for k, v in cnt.items() if (k.replace(' ', '').count('1') % 2) == 0)
        c[n] = (2 * even - tot) / tot
    F = {}
    for pair in plan['pairs']:
        if all(f'{pair}/{b}' in c for b in SIGNS):
            F[pair] = (1 + sum(SIGNS[b] * c[f'{pair}/{b}'] for b in SIGNS)) / 4

    def se(pair, shots):
        return 0.25 * math.sqrt(sum(1 - c[f'{pair}/{b}'] ** 2 for b in SIGNS) / shots)

    print('expected fidelities under the', source['model'])
    for p, f in F.items():
        print(f'  {p:16s} F {f:.3f}')
    out = {'noise_source': source, 'qpy_sha256': plan['qpy_sha256'], 'sim_shots': SIM_SHOTS, 'expected_F': F, 'z': {}}
    for shots in (4000, 8000, 16000):
        z1 = (F['main_B_not_A/B'] - 0.5) / se('main_B_not_A/B', shots)
        z2 = (F['both/A'] - 0.5) / se('both/A', shots)
        z3 = (F['main_B_not_A/B'] - F['main_B_not_A/A']) / math.hypot(se('main_B_not_A/B', shots),
                                                                        se('main_B_not_A/A', shots))
        out['z'][shots] = {'D1': z1, 'D2': z2, 'D3': z3}
        print(f'  shots {shots:6d}: expected z  D1 {z1:.1f}  D2 {z2:.1f}  D3 {z3:.1f}')
    chosen = next((s for s in (4000, 8000, 16000) if min(out['z'][s].values()) >= 5), None)
    out['shots_rule'] = chosen
    json.dump(out, open(os.path.join(HERE, 'p111_power.json'), 'w'), indent=1)
    if chosen is None:
        sys.exit('STOP: no shot count up to 16,000 reaches z = 5; Part B is not run without a design amendment')
    print('shots by the design rule:', chosen)


if __name__ == '__main__':
    main()
