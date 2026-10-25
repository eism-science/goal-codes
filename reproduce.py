"""Reproduce every number in the paper from the data stored in this repository.

    python reproduce.py                    offline check of all reported numbers (about a minute; no IBM account)
    python reproduce.py --recompute-noise  also rebuild the noise-model table from the stored code tables (slower)

Each line prints the value in the paper, the value recomputed here, and OK if they agree at the paper's precision.
Nothing in the repository is overwritten; the slow option writes to a temporary folder.
"""
import hashlib
import json
import os
import sys
import tempfile

import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(ROOT, 'phase106'))
sys.path.insert(0, os.path.join(ROOT, 'phase107'))
sys.path.insert(0, os.path.join(ROOT, 'phase111'))

FAILS = []


def check(label, paper, value, digits):
    """Agreement at the paper's precision: |value - paper| <= half a unit in the last printed digit."""
    ok = abs(value - paper) <= 0.5 * 10 ** -digits + 1e-12
    if not ok:
        FAILS.append(label)
    print(f'  {label:52s} paper {paper:<9.{digits}f} recomputed {value:<12.6f} {"OK" if ok else "MISMATCH"}')


# ---------------------------------------------------------------- provenance

def sha256(path):
    with open(path, 'rb') as f:
        return hashlib.sha256(f.read()).hexdigest()


def provenance():
    print('Design documents and compiled circuits (SHA-256 against the recorded values)')
    recorded = {
        'phase106/DESIGN.md': 'f074a8e3b0b587b1f63b19d4af5b6a28afe1a5af134603752480a62730042cc7',
        'phase107/DESIGN.md': '831b2a4dbd69026d44112064b9c757d60d0af8f540282e9acfaf3aac4a3c5670',
        'phase107/TESTB_ADDENDUM.md': '17190b858817ece3a2d0a590686b6243761cd45c2e6af025308a096cbe57a933',
    }
    with open(os.path.join(ROOT, 'phase107', 'IBM run data', 'testB_manifest.json')) as f:
        recorded['phase107/IBM run data/testB_transpiled.qpy'] = json.load(f)['qpy_sha256']
    for rel, h in recorded.items():
        ok = sha256(os.path.join(ROOT, rel)) == h
        if not ok:
            FAILS.append(rel)
        print(f'  {rel:52s} {"matches" if ok else "DOES NOT MATCH"}')
    freeze()


def freeze():
    """Phase 111: every file in the freeze manifest matches it, and the manifest is the one whose Bitcoin timestamp is
    recorded in TIMESTAMP_CHECK.json (the proof itself, FREEZE_MANIFEST.txt.ots, is checked with OpenTimestamps)."""
    d = os.path.join(ROOT, 'phase111')
    lines = open(os.path.join(d, 'FREEZE_MANIFEST.txt'), encoding='utf8').read().splitlines()
    paths = {l[4:].split(': ', 1)[0]: l[4:].split(': ', 1)[1] for l in lines if l.startswith('#   ')}
    hashes = {l.split(None, 1)[1]: l.split(None, 1)[0] for l in lines if l.strip() and not l.startswith('#')}
    bad = [n for n, h in hashes.items() if sha256(os.path.join(ROOT, paths[n])) != h]
    with open(os.path.join(d, 'TIMESTAMP_CHECK.json')) as f:
        ts = json.load(f)
    stamped = ts['manifest_sha256'] == sha256(os.path.join(d, 'FREEZE_MANIFEST.txt'))
    if bad or not stamped:
        FAILS.append('phase111 freeze')
    print(f'  {"phase111: " + str(len(hashes)) + " frozen files against FREEZE_MANIFEST.txt":52s} '
          f'{"match" if not bad else "DO NOT MATCH: " + str(bad)}')
    print(f'  {"phase111: manifest timestamped (TIMESTAMP_CHECK.json)":52s} '
          f'{"Bitcoin block " + str(ts["bitcoin_block"]) + ", " + ts["attested_utc"] if stamped else "DOES NOT MATCH"}')


# ---------------------------------------------------------------- Section 3: identical erasure (Table 1)

def identical():
    import p106
    q2, q1, c2 = p106.load_profiles()

    def gains(n, w, p):
        jq = p106.best_two(q2[n][0], n, w, p)[0]
        js = p106.sep_quantum(q1, n, w, p)[0]
        jc = p106.best_two(c2[n][0], n, w, p)[0]
        jcs = p106.sep_classical(n, w, p)[0]
        return float(jq - js), float(jc - jcs)

    print('\nTable 1: sharing gain under identical erasure (quantum | classical benchmark)')
    table = {  # (n, p): {w: (quantum, classical, digits)}
        (4, 0.1): {0.5: (0.048, 0.004, 3), 0.7: (0.048, 0.005, 3), 0.9: (0.048, 0.007, 3), 0.99: (0.048, 0, 3)},
        (5, 0.1): {0.5: (0.024, 0.004, 3), 0.7: (0.014, 0.003, 3), 0.9: (0.005, 0.001, 3), 0.99: (0.0005, 0.001, (4, 3))},
        (6, 0.1): {0.5: (0.025, 0.001, 3), 0.7: (0.010, 0.001, 3), 0.9: (0, 0.001, 3), 0.99: (0, 0.0001, 4)},
        (6, 0.2): {0.5: (0, 0.005, 3), 0.7: (0, 0.005, 3), 0.9: (0, 0.003, 3), 0.99: (0, 0.001, 3)},
    }
    for (n, p), row in table.items():
        for w, (gq, gc, d) in row.items():
            dq, dc = d if isinstance(d, tuple) else (d, d)
            q, c = gains(n, w, p)
            check(f'n={n} p={p} w={w} quantum', gq, q, dq)
            check(f'n={n} p={p} w={w} classical', gc, c, dc)
    for n in (4, 5, 6):
        for w in (0.5, 0.7, 0.9, 0.99):
            check(f'n={n} p=0.3 w={w} quantum', 0, gains(n, w, 0.3)[0], 3)
    worst_cl = max(gains(n, w, 0.3)[1] for n in (4, 5, 6) for w in (0.5, 0.7, 0.9, 0.99))
    ok = worst_cl <= 0.04 + 1e-12
    FAILS.extend([] if ok else ['p=0.3 classical <= 0.04'])
    print(f'  {"p=0.3 classical gain at most 0.04":52s} recomputed max {worst_cl:.4f} {"OK" if ok else "MISMATCH"}')

    # exact switch point, n = 6, p = 0.1: the shared optimum minus the separate optimum is piecewise linear in w
    lo, hi = 0.5, 0.99
    for _ in range(60):
        mid = (lo + hi) / 2
        (lo, hi) = (mid, hi) if gains(6, mid, 0.1)[0] > 1e-12 else (lo, mid)
    check('n=6 p=0.1 switch point w', 0.823, hi, 3)

    print('  n=6 p=0.2: optimum is the five-qubit code on A beside a bare qubit on B')
    check('  P_A of the five-qubit code at p=0.2', 0.942, p106.best_one(q1, 5, 0.2), 3)
    check('  P_B of a bare qubit at p=0.2', 0.8, p106.best_one(q1, 1, 0.2), 3)


# ---------------------------------------------------------------- Section 4: uneven and correlated erasure (Table 2)

def noise(results):
    c = results['cells']
    print('\nTable 2: mean value of knowing the goal, six physical qubits (quantum | classical benchmark)')
    table = {
        0.1: {'M0': (0.0031, 0), 'M1': (0.0027, 0), 'M1s': (0.0005, 0), 'M2': (0.0024, 0)},
        0.2: {'M0': (0, 0.0004), 'M1': (0.0049, 0.0002), 'M1s': (0, 0.0001), 'M2': (0.0148, 0.0002)},
        0.3: {'M0': (0, 0.0027), 'M1': (0.0030, 0.0010), 'M1s': (0, 0.0005), 'M2': (0, 0.0012)},
    }
    names = {'M0': 'identical', 'M1': 'readout-weighted (median)', 'M1s': 'stylized uneven', 'M2': 'correlated'}
    for p, row in table.items():
        for M, (vq, vc) in row.items():
            cell = c[f'n=6,p={p}'][M]
            check(f'p={p} {names[M]} quantum', vq, cell['Vbar_q'], 4)
            check(f'p={p} {names[M]} classical', vc, cell['Vbar_c'], 4)
    prim = c['n=6,p=0.2']
    m2 = prim['M2']['detail'][0]['q_all']['rows']
    print('  Section 4 text')
    check('  correlated, p=0.2: P_A = P_B at the balanced goal', 0.871, m2[0]['ties'][0][0], 3)
    check('  correlated, p=0.2: V(w) at w=0.99', 0.036, m2[-1]['V'], 3)
    paths = prim['M1']['Vbar_q_all_paths']
    check('  readout-weighted, p=0.2: lowest of twenty paths', 0.001, min(paths), 3)
    check('  readout-weighted, p=0.2: highest of twenty paths', 0.017, max(paths), 3)


def recompute_noise():
    import p107
    d = np.load(os.path.join(ROOT, 'phase107', 'tables_cache.npz'))
    T = {n: {f: d[f'{n}_{f}'] for f in ('q_all', 'q_sep', 'c_all', 'c_sep')} for n in p107.N_GRID}
    out = tempfile.mkdtemp(prefix='goal_codes_')
    p107.HERE = out                                  # results.json goes to the temporary folder
    p107.analyse(T)
    with open(os.path.join(out, 'results.json')) as f:
        return json.load(f)


# ---------------------------------------------------------------- Section 5: hardware (Table 3)

def bits_of(key, V):
    b = key.replace(' ', '')[::-1]                   # b[v] = bit of virtual qubit v
    return [int(b[v]) for v in range(V)]


def hardware():
    run = os.path.join(ROOT, 'phase107', 'IBM run data')
    with open(os.path.join(run, 'testB_raw_counts.json')) as f:
        raw = json.load(f)
    with open(os.path.join(run, 'testB_manifest.json')) as f:
        plan = json.load(f)['plan']                  # the plan as run; equals testB_plan.json up to sign rounding
    V = len(plan['adj'])
    counts = dict(zip(raw['circuit_names'], raw['counts']))

    # readout model from the two calibration circuits, for the secondary corrected column
    c0, c1 = counts['calib/0'], counts['calib/1']
    n0, n1 = sum(c0.values()), sum(c1.values())
    p01 = [sum(n for k, n in c0.items() if bits_of(k, V)[v] == 1) / n0 for v in range(V)]
    p10 = [sum(n for k, n in c1.items() if bits_of(k, V)[v] == 0) / n1 for v in range(V)]
    a = [1 - p01[v] - p10[v] for v in range(V)]
    b = [p10[v] - p01[v] for v in range(V)]

    def corr(key, L):
        c = plan['correlators'][key][L]
        sign = 1 if c['sign'] > 0 else -1
        support = [v for v, x in enumerate(c['pauli']) if x != 'I']
        cts = counts[f'{key}/{L}']
        n = sum(cts.values())
        k = sum(m for s, m in cts.items() if sum(bits_of(s, V)[v] for v in support) % 2 == 0)
        mean = (2 * k - n) / n
        pt = (k + 1) / (n + 2)
        se = 2 * np.sqrt(pt * (1 - pt) / n)
        cor = sum(m * np.prod([((1 - 2 * bits_of(s, V)[v]) - b[v]) / a[v] for v in support])
                  for s, m in cts.items()) / n
        return sign * mean, se, sign * cor

    print('\nTable 3: witness values on ibm_fez, recomputed from the raw counts')
    table = {  # (region key, target): (E raw, SE, E corrected)
        ('main_B_not_A', 'B'): (1.521, 0.014, 1.618), ('main_B_not_A', 'A'): (0.812, 0.019, 0.956),
        ('both', 'A'): (1.550, 0.014, 1.859), ('both', 'B'): (1.548, 0.014, 1.641),
        ('neither', 'A'): (0.811, 0.018, 0.955), ('neither', 'B'): (0.835, 0.018, 0.890),
    }
    E = {}
    for (region, target), (er, se, ec) in table.items():
        x, z = corr(f'{region}/{target}', 'X'), corr(f'{region}/{target}', 'Z')
        E[(region, target)] = (x[0] + z[0], float(np.hypot(x[1], z[1])))
        label = f'{region.split("_")[0]} / {target}'
        check(f'{label} E raw', er, E[(region, target)][0], 3)
        check(f'{label} SE', se, E[(region, target)][1], 3)
        check(f'{label} E readout-corrected', ec, x[2] + z[2], 3)
    print('  Criteria (raw counts, standard errors above the bound)')
    eb, ea, eab = E[('main_B_not_A', 'B')], E[('main_B_not_A', 'A')], E[('both', 'A')]
    check('  B1: E_B(main) above 1', 36, (eb[0] - 1) / eb[1], 0)
    check('  B2: E_A(both) above 1', 39, (eab[0] - 1) / eab[1], 0)
    check('  B3: E_B(main) above E_A(main)', 30, (eb[0] - ea[0]) / np.hypot(eb[1], ea[1]), 0)


def decoding():
    """Section 5, second run: decoded recovery and the replicated witness test, from the stored raw counts."""
    import p111_analyse as AN
    AN.N_BOOT = 200                                  # readout-corrected values are not in the paper
    d = os.path.join(ROOT, 'phase111')
    with open(os.path.join(d, 'p111_plan.json')) as f:
        plan = json.load(f)
    with open(os.path.join(d, 'run', 'B_raw_counts.json')) as f:
        b = AN.analyse_B(json.load(f)['counts'], plan)
    print('\nSection 5, second run: entanglement fidelity after region-local decoding (raw counts)')
    for name, (F, se) in {'main_B_not_A/B': (0.761, 0.005), 'both/A': (0.683, 0.006),
                          'main_B_not_A/A': (0.431, 0.006)}.items():
        check(f'{name} F', F, b['pairs'][name]['F'], 3)
        check(f'{name} SE', se, b['pairs'][name]['F_se'], 3)
    check('  abstract: recovered by the main region', 0.76, b['pairs']['main_B_not_A/B']['F'], 2)
    check('  abstract: not recovered by the main region', 0.43, b['pairs']['main_B_not_A/A']['F'], 2)
    check('  abstract: recovered by the both region', 0.68, b['pairs']['both/A']['F'], 2)
    for k, z in {'D1': 52, 'D2': 33, 'D3': 41}.items():
        check(f'  {k} standard errors', z, b['tests'][k]['z'], 0)
    ideal_half = [b['pairs'][n]['F'] for n in ('main_B_not_A/A', 'neither/A', 'neither/B')]
    check('  ideal-1/2 cases: lowest', 0.43, min(ideal_half), 2)
    check('  ideal-1/2 cases: highest', 0.47, max(ideal_half), 2)
    check('  B from the both region', 0.81, b['pairs']['both/B']['F'], 2)
    with open(os.path.join(d, 'run', 'A2_raw_counts.json')) as f:
        a2 = AN.analyse_A2(json.load(f)['counts'], json.load(open(os.path.join(ROOT, 'phase107', 'testB_plan.json'))))
    print('  Replication of the witness test (archived circuits)')
    check('  E_B(main)', 1.38, a2['witnesses']['main_B_not_A/B']['E'], 2)
    for k, z in {'B1': 24, 'B2': 40, 'B3': 23}.items():
        check(f'  {k} standard errors', z, a2['tests'][k]['z'], 0)


if __name__ == '__main__':
    provenance()
    identical()
    if '--recompute-noise' in sys.argv:
        noise(recompute_noise())
    else:
        with open(os.path.join(ROOT, 'phase107', 'results.json')) as f:
            noise(json.load(f))
    hardware()
    decoding()
    print('\nALL NUMBERS REPRODUCE' if not FAILS else f'\n{len(FAILS)} MISMATCHES: {FAILS}')
    sys.exit(1 if FAILS else 0)
