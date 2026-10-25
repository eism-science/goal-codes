"""Post-hoc check (logged in RUNLOG): random labelled 8-vertex graphs never beat the exhaustive n = 6 optimum."""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p106

q2, q1, _ = p106.load_profiles()
rng = np.random.default_rng(1061)
G = 300_000
iu = [(u, v) for u in range(8) for v in range(u + 1, 8)]
E = rng.random((G, len(iu))) < rng.uniform(0.2, 0.8, size=(G, 1))
adjs = np.zeros((G, 8), dtype=np.int64)
for e, (u, v) in enumerate(iu):
    adjs[E[:, e], u] |= 1 << v
    adjs[E[:, e], v] |= 1 << u
cr = p106.cutranks_many(adjs, 8)
pr = p106.profiles_two(cr, 8)
worst = 0.0
for p in p106.P_GRID:
    for w in p106.W_GRID:
        jex = p106.best_two(q2[6][0], 6, w, p)[0]
        jr = max(float((w * p106.prob(A[v], p) + (1 - w) * p106.prob(B[v], p)).max()) for (v, A, B) in pr.values())
        worst = max(worst, jr - jex)
print(f'{G} random graphs; largest excess of random over exhaustive optimum: {worst:.3e}')
