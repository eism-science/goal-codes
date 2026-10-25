import numpy as np
from itertools import combinations
from collections import Counter

# ── GF(2) core ──

def pauli_to_zx(p):
    n = len(p)
    z = np.zeros(n, dtype=np.uint8)
    x = np.zeros(n, dtype=np.uint8)
    for i, c in enumerate(p):
        if c == 'X': x[i] = 1
        elif c == 'Z': z[i] = 1
        elif c == 'Y': z[i] = x[i] = 1
    return z, x

def zx_to_row(z, x):
    return np.concatenate([z, x]).astype(np.uint8) % 2

def row_to_pauli(v, n):
    z, x = v[:n] % 2, v[n:] % 2
    m = {(0,0): 'I', (0,1): 'X', (1,0): 'Z', (1,1): 'Y'}
    return ''.join(m[(z[i], x[i])] for i in range(n))

def gf2_rref(A):
    A = (A.copy() % 2).astype(np.uint8)
    r, c = A.shape
    piv, row = [], 0
    for col in range(c):
        if row >= r: break
        pivot = None
        for i in range(row, r):
            if A[i, col]: pivot = i; break
        if pivot is None: continue
        if pivot != row: A[[row, pivot]] = A[[pivot, row]]
        for i in range(r):
            if i != row and A[i, col]: A[i] ^= A[row]
        piv.append(col); row += 1
    return A, piv

def gf2_rank(A):
    return len(gf2_rref(A)[1])

def gf2_nullspace(A):
    A_rref, pivs = gf2_rref(A)
    _, c = A_rref.shape
    pivset = set(pivs)
    free = [j for j in range(c) if j not in pivset]
    basis = []
    for f in free:
        x = np.zeros(c, dtype=np.uint8); x[f] = 1
        for i, pcol in enumerate(pivs):
            s = 0
            for j in free:
                if A_rref[i, j] and x[j]: s ^= 1
            x[pcol] = s
        basis.append(x)
    return basis

def gf2_solve(A, b):
    A = (A.copy() % 2).astype(np.uint8)
    b = (b.copy() % 2).astype(np.uint8)
    m, k = A.shape
    Ab = np.hstack([A, b.reshape(-1, 1)])
    row = 0; piv = [-1] * m
    for col in range(k):
        if row >= m: break
        pivot = None
        for i in range(row, m):
            if Ab[i, col]: pivot = i; break
        if pivot is None: continue
        if pivot != row: Ab[[row, pivot]] = Ab[[pivot, row]]
        for i in range(m):
            if i != row and Ab[i, col]: Ab[i] ^= Ab[row]
        piv[row] = col; row += 1
    for i in range(m):
        if not Ab[i, :k].any() and Ab[i, k]: return False, None
    x = np.zeros(k, dtype=np.uint8)
    for i in range(m):
        if piv[i] != -1: x[piv[i]] = Ab[i, k]
    return True, x

def symp_inner(v, w, n):
    return (np.dot(v[n:] % 2, w[:n] % 2) + np.dot(v[:n] % 2, w[n:] % 2)) % 2

def select_independent_rows(rows):
    if rows.size == 0: return rows
    R, piv = gf2_rref(rows)
    return R[:len(piv)] % 2

def restrict_operator_to_region(op_row, region, S_rows, n):
    region = set(region)
    outside = [q for q in range(n) if q not in region]
    S_out = np.hstack([S_rows[:, outside], S_rows[:, [q+n for q in outside]]]) % 2
    op_out = np.concatenate([op_row[outside], op_row[[q+n for q in outside]]]) % 2
    ok, a = gf2_solve(S_out.T, op_out)
    if not ok: return False, None
    op2 = op_row.copy() % 2
    for i, ai in enumerate(a):
        if ai: op2 ^= S_rows[i]
    if int(((op2[outside] | op2[[q+n for q in outside]]).sum())) != 0:
        return False, None
    return True, op2

def induced_stabilizers_on_subset(S_full, keep, n_full):
    keep = sorted(keep)
    T = [q for q in range(n_full) if q not in set(keep)]
    S_T = np.hstack([S_full[:, T], S_full[:, [q+n_full for q in T]]]) % 2
    null = gf2_nullspace(S_T.T % 2)
    if len(null) == 0:
        return np.zeros((0, 2*len(keep)), dtype=np.uint8)
    S_K = np.hstack([S_full[:, keep], S_full[:, [q+n_full for q in keep]]]) % 2
    rows = []
    for y in null:
        y = y.astype(np.uint8) % 2
        v = (y @ S_K) % 2
        if v.any(): rows.append(v.astype(np.uint8))
    if not rows:
        return np.zeros((0, 2*len(keep)), dtype=np.uint8)
    return np.vstack(rows) % 2

# ── Stabilizer projection (handles anti-commuting generators) ──

def stabilizer_project(S, g, n):
    S = S.copy() % 2; g = g.copy() % 2
    anti = [i for i in range(len(S)) if symp_inner(S[i], g, n)]
    if not anti:
        combined = np.vstack([S, g.reshape(1, -1)]) % 2
        if gf2_rank(combined) > gf2_rank(S):
            return np.vstack([S, g.reshape(1, -1)]) % 2
        return S
    j = anti[0]
    for i in anti[1:]:
        S[i] = (S[i] ^ S[j]) % 2
    S[j] = g
    return S % 2

# ── Pauli embedding ──

def embed_pauli(p, off, N):
    return 'I'*off + p + 'I'*(N - off - len(p))

def bell_pair_gens(q1, q2, N):
    xs = ['I']*N; zs = ['I']*N
    xs[q1] = xs[q2] = 'X'; zs[q1] = zs[q2] = 'Z'
    return [''.join(xs), ''.join(zs)]

def single_Z(q, N):
    p = ['I']*N; p[q] = 'Z'; return ''.join(p)

print('GF(2) utilities loaded.')