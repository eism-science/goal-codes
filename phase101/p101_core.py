"""Phase 101 core: purified stabilizer codes, exact erasure recoverability, access profiles, joint-encoding partition.

Representation. A code with n physical qubits and k logical qubits is held as its purified state: n + k qubits
(physical 0..n-1, references n..n+k-1), described by N = n + k independent commuting Pauli generators. Each generator
is a Python int: bit 2q is the X part and bit 2q+1 the Z part of qubit q. Signs are not tracked; nothing below depends
on them (all quantities are entropies of the stabilizer state).

Recoverability. Logical i is recoverable from a physical region A when I(R_i : A) = 2, equivalently (for a pure
stabilizer state) when no element of the stabilizer group supported on A-complement plus the references acts
non-trivially on R_i once restricted to A-union-references. Concretely: eliminate the generators restricted to
A-union-R on the A columns; the rows whose A part vanishes span the group elements with no A support; logical i is
recoverable iff none of them touches reference i. One elimination per region gives the answer for every logical.
This is the criterion that an element with X, respectively Z, on reference i, identity on the other references
and support otherwise in A exists.

Joint encoding. For a pure state the subsets with zero entropy form a Boolean algebra, so there is a finest product
decomposition (its atoms). For a stabilizer state, every row of the reduced row echelon form of the generator matrix
lies inside one atom (the row space is the direct sum of the atoms' row spaces on disjoint columns, pivot sets add,
and an RREF vector vanishes on every other pivot). The connected components of the RREF row supports are therefore
exactly the atoms. RREF is unique, so the result does not depend on the generator choice; atoms are a property of the
state, so local Cliffords do not change them. Logical qubits whose references share an atom are jointly encoded.
"""
import itertools
import math
import os
import random

import numpy as np


# ---------------------------------------------------------------- bit helpers

def qmask(qubits):
    m = 0
    for q in qubits:
        m |= 3 << (2 * q)
    return m


def support(v, nq):
    return [q for q in range(nq) if (v >> (2 * q)) & 3]


def gf2_rank(rows, mask=None):
    basis = {}
    r = 0
    for v in rows:
        if mask is not None:
            v &= mask
        while v:
            p = v & -v
            if p in basis:
                v ^= basis[p]
            else:
                basis[p] = v
                r += 1
                break
    return r


def entropy(rows, qubits):
    """Von Neumann entropy (bits) of a pure stabilizer state on a subset of qubits."""
    qubits = list(qubits)
    return gf2_rank(rows, qmask(qubits)) - len(qubits)


EVEN = int('01' * 128, 2)                       # bit 2q set for q < 128: the X parts


def commute(a, b):
    """Symplectic product of two Paulis as ints: 0 if they commute."""
    xa, za = a & EVEN, (a >> 1) & EVEN
    xb, zb = b & EVEN, (b >> 1) & EVEN
    return (bin(xa & zb).count('1') + bin(za & xb).count('1')) % 2


def pauli_to_int(p, offset=0):
    v = 0
    for i, c in enumerate(p):
        q = offset + i
        if c in 'XY':
            v |= 1 << (2 * q)
        if c in 'ZY':
            v |= 1 << (2 * q + 1)
    return v


# ---------------------------------------------------------------- codes

class Code:
    """Purified stabilizer code: rows (N ints), n physical qubits, k references, block layout of physical qubits."""

    def __init__(self, rows, n, k, blocks):
        self.rows, self.n, self.k, self.blocks = list(rows), n, k, [list(b) for b in blocks]
        self.N = n + k
        self.ref = [n + i for i in range(k)]
        self.ref_masks = [qmask([r]) for r in self.ref]
        self.R = qmask(self.ref)

    def copy(self):
        return Code(self.rows, self.n, self.k, self.blocks)

    def check(self):
        assert len(self.rows) == self.N, 'wrong number of generators'
        assert gf2_rank(self.rows) == self.N, 'generators not independent'
        for a, b in itertools.combinations(self.rows, 2):
            assert commute(a, b) == 0, 'generators do not commute'
        assert entropy(self.rows, self.ref) == self.k, 'references not maximally mixed'
        return True

    def recoverable(self, region_mask):
        """List of booleans: is logical i recoverable from the physical qubits in region_mask (a qmask)."""
        keep = region_mask | self.R
        basis = {}
        nulls = 0
        for r in self.rows:
            v = r & keep
            a = v & region_mask
            while a:
                p = a & -a
                if p in basis:
                    v ^= basis[p]
                    a = v & region_mask
                else:
                    basis[p] = v
                    break
            else:
                nulls |= v
        return [not (nulls & m) for m in self.ref_masks]


def initial_code(n, k, blocks):
    """Logical i in a Bell pair with physical qubit i; other physical qubits in |0>."""
    rows = []
    for i in range(k):
        rows.append((1 << (2 * i)) | (1 << (2 * (n + i))))                                         # X_i X_ri
        rows.append((1 << (2 * i + 1)) | (1 << (2 * (n + i) + 1)))                                     # Z_i Z_ri
    for j in range(k, n):
        rows.append(1 << (2 * j + 1))                                                                  # Z_j
    return Code(rows, n, k, blocks)


# ---------------------------------------------------------------- Clifford gates on a code (in place)

def _bits(v, q):
    return (v >> (2 * q)) & 1, (v >> (2 * q + 1)) & 1


def _set(v, q, x, z):
    v &= ~(3 << (2 * q))
    return v | (x << (2 * q)) | (z << (2 * q + 1))


def gate_h(code, q):
    for i, v in enumerate(code.rows):
        x, z = _bits(v, q)
        code.rows[i] = _set(v, q, z, x)


def gate_s(code, q):
    for i, v in enumerate(code.rows):
        x, z = _bits(v, q)
        code.rows[i] = _set(v, q, x, z ^ x)


def gate_cnot(code, c, t):
    for i, v in enumerate(code.rows):
        xc, zc = _bits(v, c)
        xt, zt = _bits(v, t)
        v = _set(v, t, xt ^ xc, zt)
        code.rows[i] = _set(v, c, xc, zc ^ zt)


def gate_swap(code, a, b):
    for i, v in enumerate(code.rows):
        xa, za = _bits(v, a)
        xb, zb = _bits(v, b)
        v = _set(v, a, xb, zb)
        code.rows[i] = _set(v, b, xa, za)


def random_local_clifford(code, q, rng):
    for _ in range(rng.randint(1, 6)):
        (gate_h if rng.random() < 0.5 else gate_s)(code, q)


def random_two_qubit_clifford(code, a, b, rng):
    """Random product of local Cliffords and CNOTs on (a, b); not exactly uniform, but reaches the whole group."""
    for _ in range(rng.randint(1, 3)):
        random_local_clifford(code, a, rng)
        random_local_clifford(code, b, rng)
        if rng.random() < 0.5:
            gate_cnot(code, a, b)
        else:
            gate_cnot(code, b, a)
    random_local_clifford(code, a, rng)
    random_local_clifford(code, b, rng)


def random_code(n, k, blocks, depth, rng):
    code = initial_code(n, k, blocks)
    for _ in range(depth):
        a, b = rng.sample(range(n), 2)
        random_two_qubit_clifford(code, a, b, rng)
    return code


def scramble_generators(code, rng, steps=200):
    """Random invertible row operations: same state, different generators."""
    rows = code.rows
    for _ in range(steps):
        i, j = rng.sample(range(len(rows)), 2)
        rows[i] ^= rows[j]
    rng.shuffle(rows)


# ---------------------------------------------------------------- erasures

def block_region(code, block_ids):
    return qmask([q for b in block_ids for q in code.blocks[b]])


def access_profile(code):
    """Per logical, for m = 1..B blocks: fraction of m-block subsets whose union recovers it (within-block loss off)."""
    B = len(code.blocks)
    prof = [[0.0] * B for _ in range(code.k)]
    for m in range(1, B + 1):
        subs = list(itertools.combinations(range(B), m))
        for sub in subs:
            rec = code.recoverable(block_region(code, sub))
            for i in range(code.k):
                prof[i][m - 1] += rec[i] / len(subs)
    return prof


THRESHOLD = (0.0, 0.0, 1.0, 1.0, 1.0)
LOCAL = (0.2, 0.4, 0.6, 0.8, 1.0)


def depth_index(profile):
    dl = math.dist(profile, LOCAL)
    dt = math.dist(profile, THRESHOLD)
    return 0.0 if dl + dt == 0 else (dl - dt) / (dl + dt)


def exact_patterns(n, blocks, p_b, p_q):
    """All survival patterns (as qmasks of surviving physical qubits) with their probabilities."""
    per_block = []
    for blk in blocks:
        opts = {}
        for keep in itertools.product([0, 1], repeat=len(blk)):
            s = sum(keep)
            pr = (1 - p_b) * (1 - p_q) ** s * p_q ** (len(blk) - s)
            if s == 0:
                pr += p_b
            opts[qmask([q for q, kp in zip(blk, keep) if kp])] = pr
        per_block.append(list(opts.items()))
    out = []
    for combo in itertools.product(*per_block):
        m, pr = 0, 1.0
        for mm, pp in combo:
            m |= mm
            pr *= pp
        if pr > 0:
            out.append((m, pr))
    return out


def sample_patterns(blocks, p_b, p_q, size, rng):
    out = []
    for _ in range(size):
        m = 0
        for blk in blocks:
            if rng.random() >= p_b:
                for q in blk:
                    if rng.random() >= p_q:
                        m |= 3 << (2 * q)
        out.append(m)
    return out


def recovery_probs(code, patterns, weights=None):
    tot = [0.0] * code.k
    for j, m in enumerate(patterns):
        w = 1.0 if weights is None else weights[j]
        for i, r in enumerate(code.recoverable(m)):
            if r:
                tot[i] += w
    norm = len(patterns) if weights is None else 1.0
    return [t / norm for t in tot]


def objective(code, goal, patterns, weights=None):
    return sum(w * p for w, p in zip(goal, recovery_probs(code, patterns, weights)))


# ---------------------------------------------------------------- joint encoding

def rref(rows):
    """Reduced row echelon form over GF(2), pivot = lowest set bit; unique for a given row space."""
    out = []
    for v in rows:
        for b in out:
            if v & (b & -b):
                v ^= b
        if v:
            p = v & -v
            out = [b ^ v if b & p else b for b in out]
            out.append(v)
    return sorted(out)


def atoms(code):
    """Finest product decomposition of the purified state, as a sorted list of sorted qubit lists."""
    parent = list(range(code.N))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for v in rref(code.rows):
        s = support(v, code.N)
        for q in s[1:]:
            ra, rb = find(s[0]), find(q)
            if ra != rb:
                parent[rb] = ra
    groups = {}
    for q in range(code.N):
        groups.setdefault(find(q), []).append(q)
    return sorted(sorted(g) for g in groups.values())


def logical_groups(code):
    """Partition of the logical qubits into jointly encoded groups (references sharing an atom)."""
    out = []
    for a in atoms(code):
        g = [r - code.n for r in a if r >= code.n]
        if g:
            out.append(g)
    return sorted(out)


# ---------------------------------------------------------------- known constructions (built with the qec toolbox)

TOOLBOX = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'gf2_toolbox.py')


def _load_toolbox():
    import contextlib, io
    ns = {}
    with contextlib.redirect_stdout(io.StringIO()):
        exec(open(TOOLBOX, encoding='utf-8').read(), ns)
    return ns


def _from_toolbox_rows(T, S_rows, nq_total, keep_order):
    """Convert toolbox (z|x) rows on nq_total qubits, restricted to keep_order, into int rows on len(keep_order)."""
    pos = {q: i for i, q in enumerate(keep_order)}
    ind = T['select_independent_rows'](T['induced_stabilizers_on_subset'](S_rows, sorted(keep_order), nq_total))
    nk = len(keep_order)
    sorted_keep = sorted(keep_order)
    rows = []
    for r in ind:
        v = 0
        for j, q in enumerate(sorted_keep):
            z, x = int(r[j]), int(r[j + nk])
            v |= x << (2 * pos[q])
            v |= z << (2 * pos[q] + 1)
        rows.append(v)
    return rows


def _contract(T, tensor_gens, bell_pairs, N):
    raw = []
    for gens, off in tensor_gens:
        for g in gens:
            raw.append(T['embed_pauli'](g, off, N))
    first = len(raw)
    for a, b in bell_pairs:
        raw += T['bell_pair_gens'](a, b, N)
    rows = np.array([T['zx_to_row'](*T['pauli_to_zx'](g)) for g in raw], dtype=np.uint8)
    S = rows[:first].copy()
    for g in rows[first:]:
        S = T['stabilizer_project'](S, g, N)
    return T['select_independent_rows'](S)


FIVE = ['XZZXI', 'IXZZX', 'XIXZZ', 'ZXIXZ']
FIVE_TENSOR = [g + 'I' for g in FIVE] + ['XXXXX' + 'X', 'ZZZZZ' + 'Z']       # legs 0-4, leg 5 = logical


def happy_patch_20_6():
    """The depth-2 open-bulk patch of target_swap_step1.py: blocks of 4 boundary qubits; logical 0 = centre."""
    T = _load_toolbox()
    gid = lambda t, leg: 6 * t + leg
    N = 36
    S = _contract(T, [(FIVE_TENSOR, 6 * t) for t in range(6)], [(gid(0, i), gid(1 + i, 0)) for i in range(5)], N)
    phys = [gid(t, leg) for t in range(1, 6) for leg in (1, 2, 3, 4)]
    refs = [gid(0, 5)] + [gid(t, 5) for t in range(1, 6)]
    rows = _from_toolbox_rows(T, S, N, phys + refs)
    return Code(rows, 20, 6, [list(range(4 * b, 4 * b + 4)) for b in range(5)])


C422_TENSOR = ['XXXXII', 'ZZZZII', 'XXIIXI', 'ZIZIZI', 'XIXIIX', 'ZZIIIZ']   # 4 physical, leg 4 = inner, leg 5 = outer


def concatenated_20_6():
    """Five-qubit code for the centre, each of its qubits fed with an outer logical into a [[4,2,2]] block."""
    T = _load_toolbox()
    N = 36                                       # 0-5: five-qubit tensor (5 = centre); 6 + 6b ..: block b
    tensors = [(FIVE_TENSOR, 0)] + [(C422_TENSOR, 6 + 6 * b) for b in range(5)]
    S = _contract(T, tensors, [(b, 6 + 6 * b + 4) for b in range(5)], N)
    phys = [6 + 6 * b + j for b in range(5) for j in range(4)]
    refs = [5] + [6 + 6 * b + 5 for b in range(5)]
    rows = _from_toolbox_rows(T, S, N, phys + refs)
    return Code(rows, 20, 6, [list(range(4 * b, 4 * b + 4)) for b in range(5)])


def five_qubit_code():
    rows = [pauli_to_int(g) for g in FIVE_TENSOR]
    return Code(rows, 5, 1, [[q] for q in range(5)])


# ---------------------------------------------------------------- stop rule for exact scoring (review F1, 2026-10-09)

def stop_z(n_remaining, alpha=0.05):
    """Standard-error multiplier for stopping exact scoring with n_remaining designs left unscored. Each remaining
    design's exact J exceeds its validation J plus z standard errors with probability at most 1 - Phi(z) (normal
    approximation, 10,000 validation patterns). With z = Phi^-1(1 - alpha / n_remaining), a union bound gives
    probability at most alpha that any remaining design does (Bonferroni, valid under any dependence between designs)."""
    from statistics import NormalDist
    return NormalDist().inv_cdf(1 - alpha / max(1, n_remaining))
