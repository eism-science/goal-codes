"""Phase 107 Test A: does realistic noise make the goal matter more? (design: DESIGN.md, frozen; log: RUNLOG.md)

Recoverability as Phase 106 (graph states, cut-ranks). New here: the noise is not the same on every qubit, so a code's
placement on the physical qubits matters. A code is stored as two 2^n-bit tables rec_A, rec_B over surviving regions S
(bit S set when the logical is recoverable from S); a noise model is a vector d[S] = Pr(exactly S survives); placement
pi gives P_X = sum_S rec_X(S) d[pi(S)]. Every code is evaluated at all n! placements.

Usage:
    python p107.py validate
    python p107.py analyse
"""
import itertools
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'phase106'))
import p106  # noqa: E402

W_GRID = p106.W_GRID
P_GRID = p106.P_GRID
N_GRID = p106.N_GRID
W_LOPSIDED = [w for w in W_GRID if w > 0.5]
THRESH = 0.01
TOL = 1e-12
CALIB = os.path.join(HERE, '..', 'data', 'iteration_depth_30', 'ibm_fez',
                     'calibration_snapshot.json')
N_PATHS = int(os.environ.get('P107_N_PATHS', '20'))
MAX_RATIO = 3.1


# ---------------------------------------------------------------- noise models

def dist_independent(p):
    """d[S] for independent erasure probabilities p (bit i of S set = qubit i survives)."""
    n = len(p)
    d = np.ones(1 << n)
    for S in range(1 << n):
        for i in range(n):
            d[S] *= (1 - p[i]) if (S >> i) & 1 else p[i]
    return d


def ring_edges(n):
    if n == 2:
        return [(0, 1)]
    return [(i, (i + 1) % n) for i in range(n)]


def m2_params(n, pbar):
    deg = 1 if n == 2 else 2
    s = math.sqrt(1 - pbar)
    return 1 - s, 1 - s ** (1 / deg)


def dist_m2(n, pbar):
    """Ring-correlated erasure: qubits erased independently w.p. q, each ring edge erases both ends w.p. r."""
    q, r = m2_params(n, pbar)
    E = ring_edges(n)
    d = np.zeros(1 << n)
    for F in range(1 << len(E)):
        pf = 1.0
        cov = 0
        for k, (u, v) in enumerate(E):
            if (F >> k) & 1:
                pf *= r
                cov |= (1 << u) | (1 << v)
            else:
                pf *= 1 - r
        for S in range(1 << n):
            if S & cov:
                continue
            ps = 1.0
            for i in range(n):
                if (S >> i) & 1:
                    ps *= 1 - q
                elif not (cov >> i) & 1:
                    ps *= q
            d[S] += pf * ps
    return d


def device_paths():
    """Twenty 6-qubit paths on ibm_fez (seed 107) and their readout errors, per DESIGN.md."""
    with open(CALIB) as f:
        cal = json.load(f)
    ro = [[p['value'] for p in q if p['name'] == 'readout_error'][0] for q in cal['qubits']]
    nb = {i: set() for i in range(len(ro))}
    for g in cal['gates']:
        if g['gate'] == 'cz':
            u, v = g['qubits']
            nb[u].add(v)
            nb[v].add(u)
    rng = np.random.default_rng(107)
    paths, rejected = [], 0
    while len(paths) < N_PATHS:
        path = [int(rng.integers(len(ro)))]
        while len(path) < 6:
            options = sorted(nb[path[-1]] - set(path))
            if not options:
                path = [int(rng.integers(len(ro)))]
                continue
            path.append(int(options[int(rng.integers(len(options)))]))
        e = np.array([ro[i] for i in path])
        if e.max() / e.mean() >= MAX_RATIO:
            rejected += 1
            continue
        paths.append((path, e))
    return paths, rejected


def noise_models(n, pbar, paths):
    """Name -> list of d vectors (M1 has one per path)."""
    out = {'M0': [dist_independent([pbar] * n)]}
    m1 = []
    for _, e in paths:
        ee = e[:n]
        m1.append(dist_independent(pbar * ee / ee.mean()))
    out['M1'] = m1
    out['M1s'] = [dist_independent(pbar * np.linspace(0.1, 1.9, n))]
    out['M2'] = [dist_m2(n, pbar)]
    return out


# ---------------------------------------------------------------- code tables

def pack(rec):
    """rec: (G, 2^n) bool -> (G,) uint64 tables."""
    w = np.left_shift(np.uint64(1), np.arange(rec.shape[1], dtype=np.uint64))
    return (rec.astype(np.uint64) * w).sum(axis=1, dtype=np.uint64)


def quantum_two_tables(n):
    V = n + 2
    g = p106.graphs_on(V)
    cr = p106.cutranks_many(g, V)
    tabs = []
    for a, b in itertools.permutations(range(V), 2):
        S, _ = p106.pair_indices(V, (a, b))
        base = cr[:, S].astype(np.int64)
        recA = cr[:, S | (1 << a)].astype(np.int64) == base - 1
        recB = cr[:, S | (1 << b)].astype(np.int64) == base - 1
        valid = cr[:, (1 << a) | (1 << b)] == 2
        tabs.append(np.stack([pack(recA[valid]), pack(recB[valid])], axis=1))
    return np.unique(np.concatenate(tabs), axis=0)


def quantum_one_tables(m):
    if m == 1:
        return np.array([0b10], dtype=np.uint64)            # bare qubit: recoverable only from {qubit 0}
    V = m + 1
    g = p106.graphs_on(V)
    cr = p106.cutranks_many(g, V)
    tabs = []
    for a in range(V):
        S, _ = p106.pair_indices(V, (a,))
        rec = cr[:, S | (1 << a)].astype(np.int64) == cr[:, S].astype(np.int64) - 1
        valid = cr[:, 1 << a] == 1
        tabs.append(pack(rec[valid]))
    return np.unique(np.concatenate(tabs))


def combine(tA, mA, tB, mB):
    """Separate design: A's table on qubits 0..mA-1, B's on mA..mA+mB-1."""
    n = mA + mB
    lowmask = (1 << mA) - 1
    out = np.zeros((len(tA) * len(tB), 2), dtype=np.uint64)
    k = 0
    for x in tA:
        rA = 0
        for S in range(1 << n):
            if (int(x) >> (S & lowmask)) & 1:
                rA |= 1 << S
        for y in tB:
            rB = 0
            for S in range(1 << n):
                if (int(y) >> (S >> mA)) & 1:
                    rB |= 1 << S
            out[k] = (rA, rB)
            k += 1
    return np.unique(out, axis=0)


def quantum_sep_tables(n, ones):
    return np.unique(np.concatenate([combine(ones[mA], mA, ones[n - mA], n - mA) for mA in range(1, n)]), axis=0)


def classical_two_tables(n, sep=False):
    types = [p106.T_A, p106.T_B] if sep else p106.classical_types()
    tabs = []
    for combo in itertools.combinations_with_replacement(types, n):
        if sep and (p106.T_A not in combo or p106.T_B not in combo):
            continue
        rA = rB = 0
        for S in range(1 << n):
            sub = [combo[i] for i in range(n) if (S >> i) & 1]
            rA |= p106._determines(sub, p106.T_A) << S
            rB |= p106._determines(sub, p106.T_B) << S
        full = (1 << n) - 1
        if not ((rA >> full) & 1 and (rB >> full) & 1):  # Phase 106 Deviation 1: both bits stored
            continue
        tabs.append((rA, rB))
    return np.unique(np.array(tabs, dtype=np.uint64), axis=0)


# ---------------------------------------------------------------- evaluation

def perm_maps(n):
    perms = list(itertools.permutations(range(n)))
    M = np.zeros((1 << n, len(perms)), dtype=np.int64)
    for k, pi in enumerate(perms):
        for S in range(1 << n):
            T = 0
            for i in range(n):
                if (S >> i) & 1:
                    T |= 1 << pi[i]
            M[S, k] = T
    return M


def unpack(t, n):
    return ((t[:, None] >> np.arange(1 << n, dtype=np.uint64)) & np.uint64(1)).astype(np.float64)


def candidates(tables, n, d, pm, chunk=2000):
    """All (P_A, P_B) pairs that are optimal (within 1e-9) for some w in the grid, over codes x placements."""
    D = d[pm]                                              # (2^n, n!)
    cand = []
    for s in range(0, len(tables), chunk):
        t = tables[s:s + chunk]
        PA = unpack(t[:, 0], n) @ D
        PB = unpack(t[:, 1], n) @ D
        for w in W_GRID:
            J = w * PA + (1 - w) * PB
            sel = J >= J.max() - 1e-9
            cand.append(np.stack([PA[sel], PB[sel]], axis=1))
    c = np.concatenate(cand)
    return np.unique(np.round(c, 13), axis=0)


def summarise(cand):
    """J*(w), optimal tie classes, value of knowing the goal V(w) and its mean, switches."""
    rows = []
    half = None
    for w in W_GRID:
        J = w * cand[:, 0] + (1 - w) * cand[:, 1]
        jmax = J.max()
        ties = cand[J >= jmax - TOL]
        if w == 0.5:
            half = ties
        rows.append({'w': w, 'J': float(jmax), 'ties': [tuple(float(x) for x in r) for r in ties]})
    for r in rows:
        w = r['w']
        r['V'] = r['J'] - float((w * half[:, 0] + (1 - w) * half[:, 1]).max())
    lop = [r for r in rows if r['w'] > 0.5]
    switches = sum(1 for i in range(1, len(lop)) if lop[i]['ties'] != lop[i - 1]['ties'])
    vbar = float(np.mean([r['V'] for r in lop]))
    return {'rows': rows, 'Vbar': vbar, 'switches': switches}


def best_J(cand):
    return {w: float((w * cand[:, 0] + (1 - w) * cand[:, 1]).max()) for w in W_GRID}


# ---------------------------------------------------------------- validation

def build_tables():
    t0 = time.time()
    ones = {m: quantum_one_tables(m) for m in range(1, max(N_GRID))}
    T = {}
    for n in N_GRID:
        T[n] = {'q_all': quantum_two_tables(n), 'q_sep': quantum_sep_tables(n, ones),
                'c_all': classical_two_tables(n), 'c_sep': classical_two_tables(n, sep=True)}
        print(f'n={n}: ' + ', '.join(f'{k} {len(v)}' for k, v in T[n].items()) + f'  ({time.time() - t0:.0f}s)',
              flush=True)
    return T, ones


def validate(T=None):
    ok = True
    if T is None:
        T, _ = build_tables()
    paths, rejected = device_paths()
    print(f'device paths: {len(paths)} accepted, {rejected} rejected')

    # V2: distributions
    for n in N_GRID:
        for pbar in P_GRID:
            for name, ds in noise_models(n, pbar, paths).items():
                for k, d in enumerate(ds):
                    ok &= abs(d.sum() - 1) < 1e-12
                    marg = np.array([sum(d[S] for S in range(1 << n) if not (S >> i) & 1) for i in range(n)])
                    if name == 'M0':
                        want = np.full(n, pbar)
                    elif name == 'M1':
                        e = paths[k][1][:n]
                        want = pbar * e / e.mean()
                    elif name == 'M1s':
                        want = pbar * np.linspace(0.1, 1.9, n)
                    else:
                        want = np.full(n, pbar)
                    ok &= bool(np.all(np.abs(marg - want) < 1e-12))
    print('V2 distributions and marginals:', ok)

    # V1: identical noise reproduces Phase 106
    with open(os.path.join(HERE, '..', 'phase106', 'results.json')) as f:
        r106 = json.load(f)
    worst = 0.0
    for n in N_GRID:
        pm = perm_maps(n)
        for pbar in P_GRID:
            d = dist_independent([pbar] * n)
            got = {k: best_J(candidates(T[n][k], n, d, pm)) for k in ('q_all', 'q_sep', 'c_all', 'c_sep')}
            for row in r106['cells'][f'n={n},p={pbar}']['rows']:
                w = row['w']
                for k, key in (('q_all', 'J_all'), ('q_sep', 'J_sep'), ('c_all', 'J_all_cl'), ('c_sep', 'J_sep_cl')):
                    worst = max(worst, abs(got[k][w] - row[key]))
    print(f'V1 Phase 106 reproduced, worst difference {worst:.2e}')
    ok &= worst < 1e-9

    # V3: independent recomputation for random labelled codes at n = 6
    n, V, pbar = 6, 8, 0.2
    rng = np.random.default_rng(1071)
    pm = perm_maps(n)
    e = paths[0][1]
    p_m1 = pbar * e / e.mean()
    q, r = m2_params(n, pbar)
    E = ring_edges(n)
    models = {'M1': noise_models(n, pbar, paths)['M1'][0], 'M2': dist_m2(n, pbar)}
    opt = {k: best_J(candidates(T[n]['q_all'], n, d, pm)) for k, d in models.items()}
    allset = set(map(tuple, T[n]['q_all'].tolist()))
    # explicit distributions, computed independently of dist_independent / dist_m2
    full = (1 << n) - 1
    ex1 = np.array([np.prod([(1 - p_m1[i]) if (S >> i) & 1 else p_m1[i] for i in range(n)]) for S in range(1 << n)])
    ex2 = np.zeros(1 << n)
    for ev in range(1 << (n + len(E))):                    # every qubit-erasure and edge-erasure event
        er, pe = 0, 1.0
        for i in range(n):
            if (ev >> i) & 1:
                er |= 1 << i
                pe *= q
            else:
                pe *= 1 - q
        for k2, (u, v) in enumerate(E):
            if (ev >> (n + k2)) & 1:
                er |= (1 << u) | (1 << v)
                pe *= r
            else:
                pe *= 1 - r
        ex2[full ^ er] += pe
    explicit = {'M1': ex1, 'M2': ex2}
    tested = matched = above = mismatch = 0
    while tested < 300:
        adj = [0] * V
        for u in range(V):
            for v in range(u + 1, V):
                if rng.random() < rng.uniform(0.2, 0.8):
                    adj[u] |= 1 << v
                    adj[v] |= 1 << u
        a, b = (int(x) for x in rng.choice(V, 2, replace=False))
        cr = p106.cutranks(np.array(adj, dtype=np.int64), V)
        if cr[(1 << a) | (1 << b)] != 2:
            continue
        tested += 1
        phys = [v for v in range(V) if v not in (a, b)]
        recA, recB = [], []
        for S in range(1 << n):
            m = sum(1 << phys[i] for i in range(n) if (S >> i) & 1)
            recA.append(int(cr[m | (1 << a)]) == int(cr[m]) - 1)
            recB.append(int(cr[m | (1 << b)]) == int(cr[m]) - 1)
        # explicit noise, region by region
        PA = {k: float(np.dot(recA, explicit[k])) for k in explicit}
        PB = {k: float(np.dot(recB, explicit[k])) for k in explicit}
        for key in explicit:
            for w in W_GRID:
                if w * PA[key] + (1 - w) * PB[key] > opt[key][w] + 1e-12:
                    above += 1
            # the noise-model code agrees with the explicit computation for this code's own placement
            mismatch += abs(float(np.dot(recA, models[key])) - PA[key]) > 1e-12
        # is some relabelling of this code in the enumerated set?
        sA = np.nonzero(recA)[0]
        sB = np.nonzero(recB)[0]
        for k in range(pm.shape[1]):
            pA = int(sum(1 << int(x) for x in pm[sA, k]))
            pB = int(sum(1 << int(x) for x in pm[sB, k]))
            if (pA, pB) in allset:
                matched += 1
                break
    print(f'V3 {tested} random codes: {above} above optimum, {matched} matched to enumerated set, '
          f'{mismatch} probability mismatches')
    ok &= above == 0 and matched == tested and mismatch == 0

    # V4: SEP never beats ALL
    neg = 0
    for n in N_GRID:
        pm = perm_maps(n)
        for pbar in P_GRID:
            for name, ds in noise_models(n, pbar, paths).items():
                for d in ds[:3]:
                    a_ = best_J(candidates(T[n]['q_all'], n, d, pm))
                    s_ = best_J(candidates(T[n]['q_sep'], n, d, pm))
                    neg += sum(a_[w] - s_[w] < -1e-12 for w in W_GRID)
    print('V4 negative sharing gains:', neg)
    ok &= neg == 0
    print('ALL VALIDATIONS PASS' if ok else 'VALIDATION FAILURE')
    return ok


# ---------------------------------------------------------------- analysis

def analyse(T=None):
    t0 = time.time()
    if T is None:
        T, _ = build_tables()
    paths, rejected = device_paths()
    out = {'paths': [{'qubits': p, 'readout_error': e.tolist()} for p, e in paths], 'paths_rejected': rejected,
           'cells': {}}
    for n in N_GRID:
        pm = perm_maps(n)
        for pbar in P_GRID:
            cell = {}
            for name, ds in noise_models(n, pbar, paths).items():
                per = []
                for d in ds:
                    res = {}
                    for fam in ('q_all', 'q_sep', 'c_all', 'c_sep'):
                        c = candidates(T[n][fam], n, d, pm)
                        res[fam] = summarise(c) if fam.endswith('all') else {'J': best_J(c)}
                    res['gain_q'] = {w: res['q_all']['rows'][i]['J'] - res['q_sep']['J'][w]
                                     for i, w in enumerate(W_GRID)}
                    res['gain_c'] = {w: res['c_all']['rows'][i]['J'] - res['c_sep']['J'][w]
                                     for i, w in enumerate(W_GRID)}
                    per.append(res)
                vq = [x['q_all']['Vbar'] for x in per]
                vc = [x['c_all']['Vbar'] for x in per]
                sw = [x['q_all']['switches'] for x in per]
                cell[name] = {
                    'Vbar_q': float(np.median(vq)), 'Vbar_c': float(np.median(vc)),
                    'switches_q': float(np.median(sw)), 'switches_c': float(np.median([x['c_all']['switches']
                                                                                         for x in per])),
                    'Vbar_q_all_paths': vq, 'switches_q_all_paths': sw,
                    'gain_q_w05': float(np.median([x['gain_q'][0.5] for x in per])),
                    'gain_c_w05': float(np.median([x['gain_c'][0.5] for x in per])),
                    'detail': per if name != 'M1' else per[:1],
                }
            out['cells'][f'n={n},p={pbar}'] = cell
            print(f'n={n} p={pbar} ' + '  '.join(
                f"{k}: V={v['Vbar_q']:.4f}/cl {v['Vbar_c']:.4f} sw={v['switches_q']:.0f} G.5={v['gain_q_w05']:+.4f}"
                for k, v in cell.items()) + f'  ({time.time() - t0:.0f}s)', flush=True)
    prim = out['cells']['n=6,p=0.2']
    v0, s0 = prim['M0']['Vbar_q'], prim['M0']['switches_q']
    preds = {}
    for M in ('M1', 'M2', 'M1s'):
        c = prim[M]
        preds[M] = {
            'P1_goal_matters_more': {'value': [c['Vbar_q'], v0],
                                     'holds': c['Vbar_q'] >= 2 * v0 and c['Vbar_q'] >= THRESH},
            'P2_more_switches': {'value': [c['switches_q'], s0], 'holds': c['switches_q'] > s0},
            'P3_quantum_gt_classical': {'value': [c['Vbar_q'], c['Vbar_c']], 'holds': c['Vbar_q'] > c['Vbar_c']},
        }
    out['predictions'] = preds
    out['bid_rule_pass'] = preds['M1']['P1_goal_matters_more']['holds'] or preds['M2']['P1_goal_matters_more']['holds']
    with open(os.path.join(HERE, 'results.json'), 'w') as f:
        json.dump(out, f, indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
    for M, pr in preds.items():
        print(M, pr)
    print('BID RULE PASS' if out['bid_rule_pass'] else 'BID RULE FAIL')


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'all':
        T, _ = build_tables()
        if validate(T):
            analyse(T)
    else:
        {'validate': validate, 'analyse': analyse}[cmd]()
