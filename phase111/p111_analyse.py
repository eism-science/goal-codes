"""Phase 111 analysis: decoded-recovery fidelities (Part B) and the replicated witness test (Part A2), from raw counts.

Usage:
    python p111_analyse.py B counts.json          Part B (or the pilot): fidelities, D1-D3, readout-corrected values
    python p111_analyse.py A2 counts.json         Part A2: witnesses and B1-B3, as in Test B
    python p111_analyse.py selftest               checks the estimators on synthetic counts with known answers

counts.json maps circuit names to Qiskit count dictionaries. Circuit names (fixed by p111_plan.json):
    Part B   '<region>/<target>/<XX|YY|ZZ>'            two measured bits: bit 0 = reference R, bit 1 = output s
             'calib/<R phys>-<s phys>/<ab>'            prepared a on R, b on s (ab in 00, 01, 10, 11)
             'baseline/<R phys>-<s phys>/<XX|YY|ZZ>'   direct Bell pair on that physical pair
    Part A2  the 20 Test B names ('main_B_not_A/B/X', ..., 'calib/0', 'calib/1')
Estimators and criteria: DESIGN.md, "Estimators and criteria" and "Readout correction". Raw counts decide every
criterion; corrected values are secondary.
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SIGNS = {'XX': 1, 'YY': -1, 'ZZ': 1}        # |Phi+>: F = (1 + c_XX - c_YY + c_ZZ) / 4
Z95 = 1.959964
N_BOOT = 10_000
COND_MAX = 10.0


# ---------------------------------------------------------------- raw estimators

def bits(key, n):
    """Qiskit key -> list of n bits, index i = classical bit i."""
    b = key.replace(' ', '')[::-1]
    return [int(b[i]) for i in range(n)]


def correlator(counts, support, n=2):
    """Raw parity expectation c = (2k - n_shots)/n_shots, k = even-parity shots; SE = 2 sqrt(p~(1-p~)/n_shots),
    p~ = (k + 1)/(n_shots + 2). Smoothing enters only the standard error."""
    N = sum(counts.values())
    k = sum(c for key, c in counts.items() if sum(bits(key, n)[i] for i in support) % 2 == 0)
    pt = (k + 1) / (N + 2)
    return (2 * k - N) / N, 2 * math.sqrt(pt * (1 - pt) / N), N


def fidelity(cs):
    """cs: basis -> (c, se). Returns F, SE(F)."""
    F = (1 + sum(SIGNS[b] * cs[b][0] for b in SIGNS)) / 4
    return F, 0.25 * math.sqrt(sum(cs[b][1] ** 2 for b in SIGNS))


# ---------------------------------------------------------------- readout correction (secondary)

def probs(counts):
    """Two-bit counts -> probability vector over outcomes (r, s) indexed 2*r + s."""
    v = np.zeros(4)
    for key, c in counts.items():
        r, s = bits(key, 2)
        v[2 * r + s] += c
    return v / v.sum()


def assignment(calib):
    """calib: 'ab' -> counts. A[measured, prepared], prepared index 2*a + b."""
    A = np.zeros((4, 4))
    for ab, counts in calib.items():
        A[:, 2 * int(ab[0]) + int(ab[1])] = probs(counts)
    return A


PARITY = np.array([1, -1, -1, 1])           # (-1)^(r + s) for outcomes 00, 01, 10, 11


def corrected_corr(counts, A):
    return float(PARITY @ np.linalg.solve(A, probs(counts)))


def corrected_fidelity(data, calib, rng):
    """data: basis -> counts; calib: 'ab' -> counts. Returns (F, SE by bootstrap, condition number, resamples used).
    A resample whose assignment matrix exceeds the condition cap or is singular is not used; if fewer than 99% are
    used, the SE is not reported."""
    A = assignment(calib)
    cond = float(np.linalg.cond(A))
    if cond > COND_MAX:
        return None, None, cond, 0
    F = (1 + sum(SIGNS[b] * corrected_corr(data[b], A) for b in SIGNS)) / 4
    boot = []
    for _ in range(N_BOOT):
        Ab = np.zeros((4, 4))
        for ab, counts in calib.items():
            N = sum(counts.values())
            Ab[:, 2 * int(ab[0]) + int(ab[1])] = rng.multinomial(N, probs(counts)) / N
        if np.linalg.cond(Ab) > COND_MAX:
            continue
        try:
            fb = 1.0
            for b in SIGNS:
                N = sum(data[b].values())
                pb = rng.multinomial(N, probs(data[b])) / N
                fb += SIGNS[b] * float(PARITY @ np.linalg.solve(Ab, pb))
            boot.append(fb / 4)
        except np.linalg.LinAlgError:
            continue
    se = float(np.std(boot, ddof=1)) if len(boot) >= 0.99 * N_BOOT else None
    return F, se, cond, len(boot)


# ---------------------------------------------------------------- Part B

def analyse_B(counts, plan, rng=None):
    rng = rng or np.random.default_rng(111)
    out = {'pairs': {}, 'baselines': {}}
    for name, rec in plan['pairs'].items():
        if rec.get('status') != 'built' or not all(f'{name}/{b}' in counts for b in SIGNS):
            continue
        cs = {b: correlator(counts[f'{name}/{b}'], [0, 1])[:2] for b in SIGNS}
        F, se = fidelity(cs)
        pair = rec['physical_pair']
        calib = {ab: counts[f'calib/{pair}/{ab}'] for ab in ('00', '01', '10', '11') if f'calib/{pair}/{ab}' in counts}
        Fc = sec = cond = nboot = None
        if len(calib) == 4:
            Fc, sec, cond, nboot = corrected_fidelity({b: counts[f'{name}/{b}'] for b in SIGNS}, calib, rng)
        out['pairs'][name] = {
            'correlators': {b: {'value': cs[b][0], 'se': cs[b][1]} for b in SIGNS},
            'F': F, 'F_se': se, 'F_95': [F - Z95 * se, F + Z95 * se], 'F_avg': (2 * F + 1) / 3,
            'physical_pair': pair, 'F_corrected': Fc, 'F_corrected_se': sec, 'assignment_condition': cond,
            'bootstrap_used': nboot,
            'ideal_F': 0.5 if rec.get('control') else 1.0, 'primary': rec.get('primary'),
        }
        base = {b: counts.get(f'baseline/{pair}/{b}') for b in SIGNS}
        if all(base.values()) and pair not in out['baselines']:
            bF, bse = fidelity({b: correlator(base[b], [0, 1])[:2] for b in SIGNS})
            out['baselines'][pair] = {'F': bF, 'F_se': bse}
    p = out['pairs']
    tests = {}
    if 'main_B_not_A/B' in p:
        tests['D1'] = (p['main_B_not_A/B']['F'] - 0.5, p['main_B_not_A/B']['F_se'])
    if 'both/A' in p:
        tests['D2'] = (p['both/A']['F'] - 0.5, p['both/A']['F_se'])
    if 'main_B_not_A/B' in p and 'main_B_not_A/A' in p:
        tests['D3'] = (p['main_B_not_A/B']['F'] - p['main_B_not_A/A']['F'],
                       math.hypot(p['main_B_not_A/B']['F_se'], p['main_B_not_A/A']['F_se']))
    out['tests'] = {k: {'value': v, 'se': s, 'z': v / s, 'holds': v / s > 3} for k, (v, s) in tests.items()}
    out['pass'] = len(tests) == 3 and all(t['holds'] for t in out['tests'].values())
    return out


# ---------------------------------------------------------------- Part A2 (Test B witnesses, raw)

def analyse_A2(counts, testb_plan):
    V = len(testb_plan['adj'])
    E = {}
    for key, cor in testb_plan['correlators'].items():
        cs = {}
        for L in ('X', 'Z'):
            c = cor[L]
            support = [v for v, x in enumerate(c['pauli']) if x != 'I']
            val, se, _ = correlator(counts[f'{key}/{L}'], support, n=V)
            cs[L] = ((1 if c['sign'] > 0 else -1) * val, se)
        E[key] = {'E': cs['X'][0] + cs['Z'][0], 'E_se': math.hypot(cs['X'][1], cs['Z'][1])}
    eb, ea, eab = E['main_B_not_A/B'], E['main_B_not_A/A'], E['both/A']
    tests = {'B1': (eb['E'] - 1, eb['E_se']), 'B2': (eab['E'] - 1, eab['E_se']),
             'B3': (eb['E'] - ea['E'], math.hypot(eb['E_se'], ea['E_se']))}
    out = {'witnesses': E, 'tests': {k: {'value': v, 'se': s, 'z': v / s, 'holds': v / s > 3}
                                     for k, (v, s) in tests.items()}}
    out['pass'] = all(t['holds'] for t in out['tests'].values())
    return out


# ---------------------------------------------------------------- self-test on synthetic counts

def synth_counts(c, shots, rng, flip=(0.0, 0.0)):
    """Two-qubit counts with <ZZ> = c (marginals uniform), then independent readout flips (p01 = p10 = flip[i])."""
    p_even = (1 + c) / 2
    out = {}
    for _ in range(1):
        n_even = rng.binomial(shots, p_even)
        r = np.concatenate([rng.integers(0, 2, n_even), rng.integers(0, 2, shots - n_even)])
        s = np.concatenate([r[:n_even], 1 - r[n_even:]])
        r = np.where(rng.random(shots) < flip[0], 1 - r, r)
        s = np.where(rng.random(shots) < flip[1], 1 - s, s)
        for a, b in zip(r, s):
            key = f'{b}{a}'
            out[key] = out.get(key, 0) + 1
    return out


def synth_calib(ab, shots, rng, flip):
    a, b = int(ab[0]), int(ab[1])
    r = np.where(rng.random(shots) < flip[0], 1 - a, a)
    s = np.where(rng.random(shots) < flip[1], 1 - b, b)
    out = {}
    for x, y in zip(r, s):
        out[f'{y}{x}'] = out.get(f'{y}{x}', 0) + 1
    return out


def selftest():
    rng = np.random.default_rng(1)
    shots, flip = 4000, (0.02, 0.03)
    truth = {'main_B_not_A/B': 0.80, 'main_B_not_A/A': 0.45, 'both/A': 0.78}
    plan = {'pairs': {}}
    counts = {}
    for i, (name, F) in enumerate(truth.items()):
        control = name == 'main_B_not_A/A'
        # correlators consistent with F: recoverable c_XX = c_ZZ = -c_YY = (4F - 1)/3; control ZZ only
        cs = {'XX': 0.0, 'YY': 0.0, 'ZZ': 4 * F - 1} if control else \
             {'XX': (4 * F - 1) / 3, 'YY': -(4 * F - 1) / 3, 'ZZ': (4 * F - 1) / 3}
        pair = f'{10 + i}-{20 + i}'
        for b, c in cs.items():
            counts[f'{name}/{b}'] = synth_counts(c, shots, rng, flip)
        for ab in ('00', '01', '10', '11'):
            counts[f'calib/{pair}/{ab}'] = synth_calib(ab, shots, rng, flip)
        plan['pairs'][name] = {'status': 'built', 'control': control, 'primary': True,
                               'physical_pair': pair}
    global N_BOOT
    N_BOOT = 500
    out = analyse_B(counts, plan, np.random.default_rng(2))
    shrink = (1 - 2 * flip[0]) * (1 - 2 * flip[1])
    ok = True
    for name, F in truth.items():
        r = out['pairs'][name]
        raw_expected = 0.25 + (F - 0.25) * shrink
        print(f"  {name:16s} true {F:.3f}  raw {r['F']:.4f} (expected {raw_expected:.4f}, se {r['F_se']:.4f})  "
              f"corrected {r['F_corrected']:.4f} (se {r['F_corrected_se']:.4f}, cond {r['assignment_condition']:.2f})")
        ok &= abs(r['F'] - raw_expected) < 4 * r['F_se'] and abs(r['F_corrected'] - F) < 4 * r['F_corrected_se']
    for k, t in out['tests'].items():
        print(f"  {k}: value {t['value']:+.4f}, z {t['z']:.1f}, {'holds' if t['holds'] else 'fails'}")
    # adversarial cases: bit order, a singular calibration, a missing correlator, a missing calibration
    ok &= bits('01', 2) == [1, 0] and bits('0000 0001', 8) == [1, 0, 0, 0, 0, 0, 0, 0]
    ok &= correlator({'01': 7}, [0])[0] == -1 and correlator({'01': 7}, [1])[0] == 1
    flat = {ab: {'00': 100} for ab in ('00', '01', '10', '11')}
    ok &= corrected_fidelity({b: {'00': 10} for b in SIGNS}, flat, rng)[0] is None
    for drop in ('main_B_not_A/A/YY', 'calib/10-20/00'):
        o = analyse_B({k: v for k, v in counts.items() if k != drop}, plan, np.random.default_rng(2))
        ok &= (not o['pass']) if drop.endswith('YY') else o['pairs']['main_B_not_A/B']['F_corrected'] is None
    print('SELFTEST PASS' if ok and out['pass'] else 'SELFTEST FAIL')
    return ok and out['pass']


def main():
    cmd = sys.argv[1]
    if cmd == 'selftest':
        sys.exit(0 if selftest() else 1)
    with open(sys.argv[2]) as f:
        counts = json.load(f)
    if cmd == 'B':
        with open(os.path.join(HERE, 'p111_plan.json')) as f:
            out = analyse_B(counts, json.load(f))
    elif cmd == 'A2':
        with open(os.path.join(HERE, '..', 'phase107', 'testB_plan.json')) as f:
            out = analyse_A2(counts, json.load(f))
    else:
        raise SystemExit(__doc__)
    print(json.dumps(out, indent=1, default=float))


if __name__ == '__main__':
    main()
