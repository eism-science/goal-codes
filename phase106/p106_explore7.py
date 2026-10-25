"""Phase 106, exploratory (design rule D4): n = 7. Not confirmatory.

SEP(7) is exact: single-logical [[m,1]] codes for m <= 6 are exhaustive (graphs on m + 1 <= 7 vertices, atlas).
ALL(7) is searched, not enumerated: random 9-vertex graphs, then greedy single-edge-flip ascent from the best starts
and from every n = 6 optimum extended by one vertex (all 256 neighbourhoods), per (p, w). A positive gain found is
real (it is a valid code); a zero gain is only "none found".
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p106

W = [0.5, 0.6, 0.7, 0.8, 0.9, 0.99]
EDGES = [(u, v) for u in range(9) for v in range(u + 1, 9)]


def single_best(m, p):
    if m == 1:
        return 1 - p
    g = p106.graphs_on(m + 1)
    cr = p106.cutranks_many(g, m + 1)
    best = 0.0
    for valid, c in p106.profiles_one(cr, m + 1):
        best = max(best, float(p106.prob(c[valid], p).max()))
    return best


def profiles(adjs):
    return p106.profiles_two(p106.cutranks_many(adjs, 9), 9)


def evaluate(adjs, w, p, pr=None):
    """Best J over ordered reference pairs, per graph (invalid codes score -1)."""
    if pr is None:
        pr = profiles(adjs)
    best = np.full(len(adjs), -1.0)
    for (a, b), (v, A, B) in pr.items():
        J = np.where(v, w * p106.prob(A, p) + (1 - w) * p106.prob(B, p), -1.0)
        best = np.maximum(best, J)
    return best, None


def flip(adj, u, v):
    a = adj.copy()
    a[u] ^= 1 << v
    a[v] ^= 1 << u
    return a


def ascend(adj, w, p):
    J, _ = evaluate(adj[None, :], w, p)
    J = J[0]
    while True:
        nb = np.array([flip(adj, u, v) for u, v in EDGES])
        Jn, _ = evaluate(nb, w, p)
        i = int(np.argmax(Jn))
        if Jn[i] <= J + 1e-12:
            return adj, J
        adj, J = nb[i], Jn[i]


def main():
    rng = np.random.default_rng(1067)
    q2, _, _ = p106.load_profiles()
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'profiles.json')) as f:
        reps6 = [r for _, r in json.load(f)['quantum_two']['6']]
    G = 100_000
    E = rng.random((G, len(EDGES))) < rng.uniform(0.2, 0.8, size=(G, 1))
    rand = np.zeros((G, 9), dtype=np.int64)
    for e, (u, v) in enumerate(EDGES):
        rand[E[:, e], u] |= 1 << v
        rand[E[:, e], v] |= 1 << u
    pr_rand = profiles(rand)
    print('random profiles computed', flush=True)
    out = []
    for p in p106.P_GRID:
        s = {m: single_best(m, p) for m in range(1, 7)}
        for w in W:
            sep = max(w * s[m] + (1 - w) * s[7 - m] for m in range(1, 7))
            Jr, _ = evaluate(rand, w, p, pr_rand)
            starts = [rand[i] for i in np.argsort(-Jr)[:20]]
            # extend the n = 6 optima at this (p, w) by one physical vertex, every neighbourhood
            jq, ties, _, _ = p106.best_two(q2[6][0], 6, w, p)
            for t in ties[:4]:
                r = reps6[int(t)]
                base = list(r['adj'])
                for nb in range(0, 256, 1):
                    a = base + [nb]
                    for v in range(8):
                        if (nb >> v) & 1:
                            a[v] |= 1 << 8
                    starts.append(np.array(a, dtype=np.int64))
            ext = np.array(starts[20:])
            Je, _ = evaluate(ext, w, p) if len(ext) else (np.array([]), None)
            best_start = float(max(Jr.max(), Je.max() if len(Je) else -1))
            best = best_start
            for st in starts[:20] + [ext[i] for i in np.argsort(-Je)[:10]]:
                _, J = ascend(st, w, p)
                best = max(best, float(J))
            row = {'p': p, 'w': w, 'J_sep_exact': sep, 'J_all_found': best, 'gain_found': best - sep}
            out.append(row)
            print(f"p={p} w={w:.2f}  sep(exact)={sep:.5f}  all(found)={best:.5f}  gain={best - sep:+.5f}", flush=True)
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'explore7.json'), 'w') as f:
        json.dump(out, f, indent=1)


if __name__ == '__main__':
    main()
