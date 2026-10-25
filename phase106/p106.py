"""Phase 106: does sharing still pay when the goal is lopsided? (design: DESIGN.md, frozen; log: RUNLOG.md)

Codes are purified stabilizer states, taken without loss as graph states (local-Clifford equivalence leaves every
mutual information unchanged). A graph on V vertices is a list of adjacency bitmasks; the entropy of a vertex subset
X is the GF(2) rank of the adjacency submatrix Gamma[X, complement of X] (the cut-rank). Logical X with reference
vertex r is recoverable from a physical region S when I(r : S) = 2, i.e. cr(S + r) = cr(S) - 1 (given cr(r) = 1).

Usage:
    python p106.py validate      V1-V5
    python p106.py enumerate     exhaustive profiles, writes profiles_q.npz and profiles_cl.json
    python p106.py analyse       grid, predictions, writes results.json
"""
import itertools
import json
import math
import os
import sys
import time

import numpy as np
from numba import njit

HERE = os.path.dirname(os.path.abspath(__file__))
W_GRID = [0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.99]
P_GRID = [0.1, 0.2, 0.3]
N_GRID = [2, 3, 4, 5, 6]
LAMBDAS = [0.0, 0.01, 0.02, 0.05]
THRESH = 0.01
TOL = 1e-12


# ---------------------------------------------------------------- cut-rank

@njit(cache=True)
def _rank_rows(rows, nrows):
    """GF(2) rank of nrows uint16 bitmask rows (destroys the array)."""
    r = 0
    for i in range(nrows):
        v = rows[i]
        if v == 0:
            continue
        p = v & (-v)
        r += 1
        for j in range(i + 1, nrows):
            if rows[j] & p:
                rows[j] ^= v
    return r


@njit(cache=True)
def cutranks(adj, V):
    """cr[X] for every subset X of V vertices (adj: int64 adjacency bitmasks)."""
    full = (1 << V) - 1
    out = np.zeros(1 << V, dtype=np.int8)
    rows = np.zeros(V, dtype=np.int64)
    for X in range(1 << V):
        comp = full ^ X
        k = 0
        for v in range(V):
            if (X >> v) & 1:
                rows[k] = adj[v] & comp
                k += 1
        out[X] = _rank_rows(rows, k)
    return out


@njit(cache=True)
def cutranks_many(adjs, V):
    G = adjs.shape[0]
    out = np.zeros((G, 1 << V), dtype=np.int8)
    for g in range(G):
        out[g] = cutranks(adjs[g], V)
    return out


def popcount_table(V):
    return np.array([bin(x).count('1') for x in range(1 << V)], dtype=np.int64)


# ---------------------------------------------------------------- profiles

def pair_indices(V, refs):
    """For physical = all vertices minus refs: the subset masks S (as V-bit masks) and their sizes."""
    phys = [v for v in range(V) if v not in refs]
    masks, sizes = [], []
    for bits in range(1 << len(phys)):
        m = 0
        for i, v in enumerate(phys):
            if (bits >> i) & 1:
                m |= 1 << v
        masks.append(m)
        sizes.append(bin(bits).count('1'))
    return np.array(masks, dtype=np.int64), np.array(sizes, dtype=np.int64)


def profiles_two(cr, V):
    """cr: (G, 2^V). For every ordered reference pair (a, b): validity and the count profiles c_A, c_B (length n+1).
    Returns dict (a, b) -> (valid bool (G,), cA (G, n+1), cB (G, n+1))."""
    n = V - 2
    out = {}
    for a, b in itertools.permutations(range(V), 2):
        S, sz = pair_indices(V, (a, b))
        onehot = np.zeros((len(S), n + 1), dtype=np.int64)
        onehot[np.arange(len(S)), sz] = 1
        base = cr[:, S].astype(np.int64)
        recA = cr[:, S | (1 << a)].astype(np.int64) == base - 1
        recB = cr[:, S | (1 << b)].astype(np.int64) == base - 1
        valid = (cr[:, (1 << a) | (1 << b)] == 2)
        out[(a, b)] = (valid, recA.astype(np.int64) @ onehot, recB.astype(np.int64) @ onehot)
    return out


def profiles_one(cr, V):
    """Single-logical codes: reference a, physical = the rest. Returns list of (valid (G,), c (G, m+1))."""
    m = V - 1
    res = []
    for a in range(V):
        S, sz = pair_indices(V, (a,))
        onehot = np.zeros((len(S), m + 1), dtype=np.int64)
        onehot[np.arange(len(S)), sz] = 1
        rec = cr[:, S | (1 << a)].astype(np.int64) == cr[:, S].astype(np.int64) - 1
        valid = cr[:, 1 << a] == 1
        res.append((valid, rec.astype(np.int64) @ onehot))
    return res


def prob(c, p):
    """P(p) for count profile(s) c with last axis s = 0..n."""
    c = np.asarray(c, dtype=np.float64)
    n = c.shape[-1] - 1
    s = np.arange(n + 1)
    return c @ ((1 - p) ** s * p ** (n - s))


# ---------------------------------------------------------------- graphs

def atlas_graphs(V):
    from networkx.generators.atlas import graph_atlas_g
    out = []
    for g in graph_atlas_g():
        if g.number_of_nodes() != V:
            continue
        adj = [0] * V
        for u, v in g.edges():
            adj[u] |= 1 << v
            adj[v] |= 1 << u
        out.append(adj)
    return np.array(out, dtype=np.int64)


def graphs_on(V):
    """Every graph on V vertices up to isomorphism (atlas for V <= 7; atlas-7 plus one vertex for V = 8)."""
    if V <= 7:
        return atlas_graphs(V)
    assert V == 8
    base = atlas_graphs(7)
    out = np.zeros((len(base) * 128, 8), dtype=np.int64)
    k = 0
    for adj in base:
        for nb in range(128):
            a = list(adj) + [nb]
            for v in range(7):
                if (nb >> v) & 1:
                    a[v] |= 1 << 7
            out[k] = a
            k += 1
    return out


# ---------------------------------------------------------------- classical twin

INPUTS = [(0, 0), (0, 1), (1, 0), (1, 1)]          # (x_A, x_B), index 2*x_A + x_B


def classical_types():
    seen, out = set(), []
    for t in range(1, 15):
        c = min(t, 15 ^ t)
        if c not in seen:
            seen.add(c)
            out.append(c)
    return sorted(out)


T_A, T_B = 0b1100, 0b1010                          # truth tables of x_A and x_B (bit i = value on input i)


def _determines(tables, target):
    for u, v in itertools.combinations(range(4), 2):
        if ((target >> u) & 1) != ((target >> v) & 1):
            if all(((t >> u) & 1) == ((t >> v) & 1) for t in tables):
                return False
    return True


def classical_profile(types):
    n = len(types)
    cA, cB = [0] * (n + 1), [0] * (n + 1)
    for bits in range(1 << n):
        S = [types[i] for i in range(n) if (bits >> i) & 1]
        s = len(S)
        cA[s] += _determines(S, T_A)
        cB[s] += _determines(S, T_B)
    return cA, cB


# ---------------------------------------------------------------- validation

def validate():
    sys.path.insert(0, os.path.join(HERE, '..', 'phase101'))
    import p101_core as core
    rng = np.random.default_rng(106)
    ok = True

    # V1: graph-state recoverability against p101_core.Code.recoverable
    checked, mism = 0, 0
    graphs_tested = 0
    while graphs_tested < 2500:
        n = int(rng.integers(2, 7))
        V = n + 2
        adj = [0] * V
        for u in range(V):
            for v in range(u + 1, V):
                if rng.random() < 0.5:
                    adj[u] |= 1 << v
                    adj[v] |= 1 << u
        cr = cutranks(np.array(adj, dtype=np.int64), V)
        if cr[(1 << n) | (1 << (n + 1))] != 2:
            continue
        graphs_tested += 1
        rows = []
        for v in range(V):
            r = 1 << (2 * v)
            for u in range(V):
                if (adj[v] >> u) & 1:
                    r |= 1 << (2 * u + 1)
            rows.append(r)
        code = core.Code(rows, n, 2, [list(range(n))])
        code.check()
        for bits in range(1 << n):
            region = [q for q in range(n) if (bits >> q) & 1]
            ref = code.recoverable(core.qmask(region))
            base = int(cr[bits])
            mine = [int(cr[bits | (1 << (n + i))]) == base - 1 for i in range(2)]
            checked += 1
            mism += mine != list(ref)
    print(f'V1: {graphs_tested} graphs, {checked} regions, {mism} mismatches')
    ok &= mism == 0

    # V2: known profiles
    def prof1(adj, V, a):
        cr = cutranks(np.array(adj, dtype=np.int64), V)
        return [list(x[0]) for x in [profiles_one(cr[None, :], V)[a][1]]][0], cr
    # bare qubit: Bell pair between reference 0 and physical 1
    c, _ = prof1([0b10, 0b01], 2, 0)
    ok &= c == [0, 1]
    print('V2 bare qubit', c)
    # [[5,1,3]]: the 5-qubit code's purification is LC-equivalent to the 6-vertex graph (pentagon + centre)
    adj = [0] * 6
    for i in range(5):
        for j in ((i + 1) % 5, (i - 1) % 5, 5):
            adj[i] |= 1 << j
            adj[j] |= 1 << i
    c, _ = prof1(adj, 6, 5)
    want = [0, 0, 0, 10, 5, 1]
    print('V2 [[5,1,3]]', c, 'want', want)
    ok &= c == want
    # [[4,2,2]]: search all 6-vertex graphs for a code whose two logicals each recover from every 3 of 4 and no 2
    g6 = graphs_on(6)
    cr6 = cutranks_many(g6, 6)
    pr = profiles_two(cr6, 6)
    found = any(bool(np.any(v & np.all(A == [0, 0, 0, 4, 1], axis=1) & np.all(B == [0, 0, 0, 4, 1], axis=1)))
                for (v, A, B) in pr.values())
    print('V2 [[4,2,2]] profile found among 6-vertex graphs:', found)
    ok &= found
    # Bell pair between two physical qubits, references isolated: invalid
    cr = cutranks(np.array([0b10, 0b01, 0, 0], dtype=np.int64), 4)
    print('V2 invalid code rejected:', cr[0b1100] != 2)
    ok &= cr[0b1100] != 2

    # V3: enumeration size
    n7 = len(atlas_graphs(7))
    n8 = len(graphs_on(8))
    print('V3 atlas-7', n7, 'extensions', n8)
    ok &= n7 == 1044 and n8 == 133632

    # V4: classical
    for m in range(1, 5):
        cA, _ = classical_profile([T_A] * m)
        ok &= abs(prob(cA, 0.2) - (1 - 0.2 ** m)) < 1e-12
    cA, cB = classical_profile([0b0110])
    print('V4 parity alone', cA, cB)
    ok &= cA == [0, 0] and cB == [0, 0]
    print('V4 classical types', [bin(t) for t in classical_types()])
    ok &= len(classical_types()) == 7

    # V5: single-logical profiles under the no-cloning bound
    for m in range(1, 7):
        g = graphs_on(m + 1)
        cr = cutranks_many(g, m + 1)
        for valid, c in profiles_one(cr, m + 1):
            c = c[valid]
            for s in range(m + 1):
                if 2 * s < m:
                    bad = np.any(c[:, s] + c[:, m - s] > math.comb(m, s))
                    ok &= not bad
                elif 2 * s == m:
                    ok &= not np.any(2 * c[:, s] > math.comb(m, s))
    print('V5 no-cloning bound respected')
    print('ALL VALIDATIONS PASS' if ok else 'VALIDATION FAILURE')
    return ok


# ---------------------------------------------------------------- enumeration

def enumerate_all():
    t0 = time.time()
    res = {'two': {}, 'one': {}}
    reps = {}
    for n in N_GRID:
        V = n + 2
        g = graphs_on(V)
        t = time.time()
        cr = cutranks_many(g, V)
        pr = profiles_two(cr, V)
        keys = {}
        for (a, b), (valid, A, B) in pr.items():
            idx = np.nonzero(valid)[0]
            AB = np.concatenate([A[idx], B[idx]], axis=1)
            uniq, first = np.unique(AB, axis=0, return_index=True)
            for row, f in zip(uniq, first):
                k = tuple(int(x) for x in row)
                if k not in keys:
                    gi = int(idx[f])
                    keys[k] = {'adj': [int(x) for x in g[gi]], 'A': a, 'B': b}
        res['two'][n] = keys
        print(f'n={n}: {len(g)} graphs, {len(keys)} distinct (cA,cB) profiles, {time.time() - t:.1f}s', flush=True)
    for m in range(1, max(N_GRID)):
        V = m + 1
        if V == 1:
            res['one'][m] = {(0, 1): None}
            continue
        g = graphs_on(V)
        cr = cutranks_many(g, V)
        keys = set()
        for valid, c in profiles_one(cr, V):
            for row in np.unique(c[valid], axis=0):
                keys.add(tuple(int(x) for x in row))
        res['one'][m] = {k: None for k in keys}
        print(f'[[{m},1]]: {len(keys)} distinct profiles')
    # m = 1 single qubit: bare
    res['one'][1] = {(0, 1): None}

    cl = {}
    types = classical_types()
    for n in N_GRID:
        prof = {}
        for combo in itertools.combinations_with_replacement(types, n):
            cA, cB = classical_profile(list(combo))
            k = tuple(cA + cB)
            if k not in prof:
                prof[k] = list(combo)
        cl[n] = prof
        print(f'classical n={n}: {len(prof)} distinct profiles')
    out = {
        'quantum_two': {str(n): [[list(k), v] for k, v in d.items()] for n, d in res['two'].items()},
        'quantum_one': {str(m): [list(k) for k in d] for m, d in res['one'].items()},
        'classical_two': {str(n): [[list(k), v] for k, v in d.items()] for n, d in cl.items()},
    }
    with open(os.path.join(HERE, 'profiles.json'), 'w') as f:
        json.dump(out, f)
    print(f'enumeration done in {time.time() - t0:.0f}s')


# ---------------------------------------------------------------- analysis

def load_profiles():
    with open(os.path.join(HERE, 'profiles.json')) as f:
        d = json.load(f)
    q2 = {int(n): (np.array([k for k, _ in v]), [r for _, r in v]) for n, v in d['quantum_two'].items()}
    q1 = {int(m): np.array(v) for m, v in d['quantum_one'].items()}
    c2 = {int(n): (np.array([k for k, _ in v]), [r for _, r in v]) for n, v in d['classical_two'].items()}
    if not CLASSICAL_UNRESTRICTED:
        # Deviation 1 (RUNLOG): a classical code must store both bits (all n coordinates determine x_A and x_B),
        # as every quantum code does (S(R) = 2); otherwise ALL_cl contains codes that drop B, which SEP_cl and the
        # quantum families do not.
        for n in list(c2):
            AB, reps = c2[n]
            keep = (AB[:, n] == 1) & (AB[:, 2 * n + 1] == 1)
            c2[n] = (AB[keep], [r for r, k in zip(reps, keep) if k])
    return q2, q1, c2


CLASSICAL_UNRESTRICTED = os.environ.get('P106_CLASSICAL_UNRESTRICTED') == '1'


def best_two(AB, n, w, p):
    PA, PB = prob(AB[:, :n + 1], p), prob(AB[:, n + 1:], p)
    J = w * PA + (1 - w) * PB
    jmax = J.max()
    ties = np.nonzero(J >= jmax - TOL)[0]
    return jmax, ties, PA, PB


def best_one(q1, m, p):
    return float(prob(q1[m], p).max())


def sep_quantum(q1, n, w, p):
    best, arg = -1, None
    for mA in range(1, n):
        mB = n - mA
        J = w * best_one(q1, mA, p) + (1 - w) * best_one(q1, mB, p)
        if J > best + TOL:
            best, arg = J, (mA, mB)
    return best, arg


def sep_classical(n, w, p):
    best, arg = -1, None
    for mA in range(1, n):
        mB = n - mA
        J = w * (1 - p ** mA) + (1 - w) * (1 - p ** mB)
        if J > best + TOL:
            best, arg = J, (mA, mB)
    return best, arg


def canon(AB, ties, n):
    return sorted(tuple(int(x) for x in AB[t]) for t in ties)


def analyse():
    q2, q1, c2 = load_profiles()
    out = {'cells': {}, 'cost': {}}
    for p in P_GRID:
        for n in N_GRID:
            AB, reps = q2[n]
            ABc, repc = c2[n]
            rows = []
            j_half_ties = None
            for w in W_GRID:
                jq, ties, PA, PB = best_two(AB, n, w, p)
                js, sarg = sep_quantum(q1, n, w, p)
                jc, tiesc, _, _ = best_two(ABc, n, w, p)
                jcs, csarg = sep_classical(n, w, p)
                if w == 0.5:
                    j_half_ties = ties
                regret = jq - float((w * PA + (1 - w) * PB)[j_half_ties].max())
                rows.append({
                    'w': w, 'J_all': jq, 'J_sep': js, 'gain': jq - js, 'sep_alloc': sarg,
                    'opt_profiles': canon(AB, ties, n), 'opt_example': reps[int(ties[0])],
                    'PA': float(PA[ties[0]]), 'PB': float(PB[ties[0]]),
                    'regret_of_w05_optimum': regret,
                    'J_all_cl': jc, 'J_sep_cl': jcs, 'gain_cl': jc - jcs, 'sep_alloc_cl': csarg,
                    'opt_profiles_cl': canon(ABc, tiesc, n), 'opt_example_cl': repc[int(tiesc[0])],
                })
            changes = [W_GRID[i] for i in range(1, len(W_GRID))
                       if rows[i]['opt_profiles'] != rows[i - 1]['opt_profiles']]
            out['cells'][f'n={n},p={p}'] = {'rows': rows, 'breakpoints_after': changes}
        # cost per qubit: best n per w
        for lam in LAMBDAS:
            res = []
            for w in W_GRID:
                best = None
                for n in N_GRID:
                    jq, ties, _, _ = best_two(q2[n][0], n, w, p)
                    js, _ = sep_quantum(q1, n, w, p)
                    for fam, j in (('shared', jq), ('separate', js)):
                        v = j - lam * n
                        if best is None or v > best[0] + TOL:
                            best = (v, n, 'separate' if js >= jq - TOL else fam)
                res.append({'w': w, 'J_lambda': best[0], 'n': best[1], 'family': best[2]})
            out['cost'][f'p={p},lambda={lam}'] = res

    prim = out['cells']['n=6,p=0.2']
    r = {x['w']: x for x in prim['rows']}
    g = [x['gain'] for x in prim['rows']]
    preds = {
        'P1_profile_changes': {'value': prim['breakpoints_after'], 'holds': len(prim['breakpoints_after']) > 0},
        'P2_gain_w0.5': {'value': r[0.5]['gain'], 'holds': r[0.5]['gain'] > THRESH},
        'P3_gain_w0.9': {'value': r[0.9]['gain'], 'holds': r[0.9]['gain'] > THRESH, 'predicted': 'fails'},
        'P4_nonincreasing': {'value': g, 'holds': all(g[i + 1] <= g[i] + TOL for i in range(len(g) - 1))},
        'P5_quantum_gt_classical_w0.5': {'value': [r[0.5]['gain'], r[0.5]['gain_cl']],
                                         'holds': r[0.5]['gain'] > r[0.5]['gain_cl']},
    }
    out['predictions'] = preds
    name = 'results_classical_unrestricted.json' if CLASSICAL_UNRESTRICTED else 'results.json'
    with open(os.path.join(HERE, name), 'w') as f:
        json.dump(out, f, indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
    for k, v in preds.items():
        print(k, v)
    print()
    for key, cell in out['cells'].items():
        print(key, 'breakpoints after', cell['breakpoints_after'])
        for x in cell['rows']:
            print(f"  w={x['w']:.2f}  all={x['J_all']:.5f} sep={x['J_sep']:.5f} gain={x['gain']:+.5f}"
                  f"  | cl all={x['J_all_cl']:.5f} sep={x['J_sep_cl']:.5f} gain={x['gain_cl']:+.5f}"
                  f"  PA={x['PA']:.4f} PB={x['PB']:.4f} regret={x['regret_of_w05_optimum']:.5f}")


if __name__ == '__main__':
    {'validate': validate, 'enumerate': enumerate_all, 'analyse': analyse}[sys.argv[1]]()
