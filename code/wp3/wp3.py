#!/usr/bin/env python
"""wp3.py -- WP3 certificate build for the K3 Kummer approximant: residual and Ricci bounds OUTSIDE the cores,
uniform in the neck size a in (0, 0.03] (DEFINITIONS-AND-LEMMAS-v2.md, cut-off window [0.25 a, a] of WP2a, so the
region 'distance >= a from every bubble' is exactly where P_a = psi).

Certified outputs (a in (0, A0], three a-pieces, one touching 0):
  C_out = sup |F| / a^8          (F = log det g, q = eps^2 Q(z; eps), L4 (d))
  R_out = sup |Ric|_g / a^6      (Ric = -i d dbar F)
  e_out = sup |E|_{g0}           (by-product, positivity norm of 1.6)
  a certified lower bound for sup |F8bulk| (point evaluation, Ewald in arb)
Method (details in BUILD-WP3.md):
  * exact homogeneous Taylor polynomials of psi_1, psi_2 at the bubble (lattice moment sums, Gegenbauer form),
    near shells (0 < |l| <= 1, 168 points) in closed form, far part (|l| > 1) as a degree-22 polynomial;
  * explicit tails: A2 (Lie-norm Gegenbauer bounds) + Cauchy on polydiscs for derivatives (homogeneity trick);
  * NECK (|z| <= R_NR): scale-covariant variables (omega = z/|z|, u = eps/|z|^4 in [0,1], r = |z|); leading
    Ricci term in closed form, branch and bound with centred forms; corrections by interval jets;
  * BULK (|z| >= R_NF): Taylor models of order 4 in 4 real variables around box centres (all cancellations of
    the bubble at 0 happen inside exact coefficient sums), Ricci by interval 2-jets;
  * quotient by G_J (order 192) with WP2a's Dirichlet cone; V0 by its 24 facets.
Modes
  --controls OUT.json   J1-J4 controls, mutation tests, NaN/inf guards, timing of a REDUCED cover
  --main OUT.json       certified run; refuses unless <lab>/wp3/GO-WP3 exists; chunked (exit 3 = resume)
Run:  python <lab-tools>/memcap.py 8 python <lab>/wp3/wp3.py --controls <lab>/wp3/controls_wp3.json
Rigor: arb (python-flint, 128 bits) for constants; numpy intervals with one-ulp outward rounding for +,-,*,/,sqrt;
matrix products in midpoint-radius form with the gamma_n error bound (any summation order).
"""
import os, sys, json, time, math, argparse, itertools, pickle
import numpy as np
from flint import arb, ctx, fmpq

ctx.prec = 128
HERE = os.path.dirname(os.path.abspath(__file__))
GO_FILE = os.path.join(HERE, 'GO-WP3')
STATE_DIR = os.path.join(HERE, 'state_wp3')
A0 = 0.03
A_PIECES = [(0.0, 0.01), (0.01, 0.02), (0.02, 0.03)]
R_NF = 0.03          # F:     neck |z| <= R_NF ((15/16) c6 + corrections), bulk |z| >= R_NF (Taylor models)
R_NR = 0.08          # Ricci: neck |z| <= R_NR (order-r^2 closed form B&B + O(r^4) rest), bulk |z| >= R_NR (2-jets)
RM = 8.0             # radius of the explicit lattice sums for the moment tables (tails by Lemma A4)
DF1, DF3 = 22, 20    # degrees of the far polynomials (nu = 1, 3)
DN1, DN3 = 18, 16    # degrees of the neck polynomials (nu = 1, 3)
KP = 8               # order of the source Taylor expansions at box centres
KT = 4               # order of the Taylor models
COV = math.sqrt(2) / 4   # covering radius of L (1/(2 sqrt 2))
MUT = {'drop_l0_cancel': False, 'drop_N1': False}   # mutation switches (controls only)

# =====================================================================================================
# 1. interval arithmetic (I from WP2a) + midpoint-radius matrix products
# =====================================================================================================
_NINF, _PINF = -np.inf, np.inf
UR = 2.0 ** -53
_EPS52 = 2.0 ** -52
def _dn(a):
    """a value <= nextafter(a, -inf): a - (|a| 2^-52 + 1e-300) (the subtrahend is >= one ulp of a, rounding is
       monotone, so the result is <= the next float below a).  Cheaper than np.nextafter."""
    return a - (np.abs(a) * _EPS52 + 1e-300)
def _up(a):
    return a + (np.abs(a) * _EPS52 + 1e-300)

class I:
    """interval [lo, hi] (numpy arrays); every op rounds outward by one ulp."""
    __slots__ = ('lo', 'hi')
    def __init__(self, lo, hi=None):
        self.lo = np.asarray(lo, dtype=float)
        self.hi = self.lo if hi is None else np.asarray(hi, dtype=float)
    @staticmethod
    def c(v):
        """enclosure of a python float/int/Fraction/arb constant"""
        if isinstance(v, arb):
            return I(float(_dn(float(v.lower()))), float(_up(float(v.upper()))))
        if isinstance(v, (int, np.integer)) and abs(int(v)) < 2 ** 53:
            return I(float(v))
        f = float(v)
        if isinstance(v, float) or f == v:
            return I(f)
        return I(float(_dn(f)), float(_up(f)))
    def __add__(s, o):
        if isinstance(o, I): return I(_dn(s.lo + o.lo), _up(s.hi + o.hi))
        if not isinstance(o, (int, float, np.floating, np.integer, np.ndarray)): return NotImplemented
        o = np.asarray(o, float); return I(_dn(s.lo + o), _up(s.hi + o))
    __radd__ = __add__
    def __neg__(s): return I(-s.hi, -s.lo)
    def __sub__(s, o):
        if isinstance(o, I): return I(_dn(s.lo - o.hi), _up(s.hi - o.lo))
        o = float(o); return I(_dn(s.lo - o), _up(s.hi - o))
    def __rsub__(s, o): return (-s) + o
    def __mul__(s, o):
        if isinstance(o, I):
            a = s.lo * o.lo; b = s.lo * o.hi; c = s.hi * o.lo; d = s.hi * o.hi
            return I(_dn(np.minimum(np.minimum(a, b), np.minimum(c, d))), _up(np.maximum(np.maximum(a, b), np.maximum(c, d))))
        if isinstance(o, np.ndarray):
            return s * I(o, o)
        if not isinstance(o, (int, float, np.floating, np.integer)): return NotImplemented
        o = float(o)
        if o >= 0: return I(_dn(s.lo * o), _up(s.hi * o))
        return I(_dn(s.hi * o), _up(s.lo * o))
    __rmul__ = __mul__
    def __truediv__(s, o):
        if isinstance(o, I): return s * irecip(o)
        return s * irecip(I(float(o)))
    def __getitem__(s, idx): return I(s.lo[idx], s.hi[idx])
    def mag(s): return np.maximum(np.abs(s.lo), np.abs(s.hi))
    def mig(s): return np.where((s.lo <= 0) & (s.hi >= 0), 0.0, np.minimum(np.abs(s.lo), np.abs(s.hi)))
    def width(s): return s.hi - s.lo
    @property
    def shape(s): return s.lo.shape
    def __repr__(s): return 'I(%s, %s)' % (s.lo, s.hi)

def isq(x):
    if not isinstance(x, I): return x * x
    a = x.lo * x.lo; b = x.hi * x.hi
    z = (x.lo <= 0) & (x.hi >= 0)
    return I(np.where(z, 0.0, _dn(np.minimum(a, b))), _up(np.maximum(a, b)))
def isqrt(x):
    return I(np.maximum(_dn(np.sqrt(np.maximum(x.lo, 0.0))), 0.0), _up(np.sqrt(np.maximum(x.hi, 0.0))))
def irecip(x):
    if np.any((x.lo <= 0) & (x.hi >= 0)):
        raise FloatingPointError('irecip: interval contains 0')
    return I(_dn(1.0 / x.hi), _up(1.0 / x.lo))
def iabs(x):
    z = (x.lo <= 0) & (x.hi >= 0)
    return I(np.where(z, 0.0, np.minimum(np.abs(x.lo), np.abs(x.hi))), np.maximum(np.abs(x.lo), np.abs(x.hi)))
def imax0(x): return I(np.maximum(x.lo, 0.0), np.maximum(x.hi, 0.0))
def isect(a, b): return I(np.maximum(a.lo, b.lo), np.minimum(a.hi, b.hi))
def iwhere(m, a, b): return I(np.where(m, a.lo, b.lo), np.where(m, a.hi, b.hi))
def ihull(a, b): return I(np.minimum(a.lo, b.lo), np.maximum(a.hi, b.hi))
def ipow(x, k):
    out = I(np.ones_like(x.lo))
    for _ in range(k): out = out * x
    return out
def sym(m): return I(-np.asarray(m, float), np.asarray(m, float))

def check_finite(name, *ivs):
    for k, x in enumerate(ivs):
        lo, hi = np.asarray(x.lo), np.asarray(x.hi)
        if not (np.all(np.isfinite(lo)) and np.all(np.isfinite(hi))):
            raise FloatingPointError('%s[%d]: non-finite interval endpoint' % (name, k))
        if np.any(lo > hi):
            raise FloatingPointError('%s[%d]: empty interval (lo > hi)' % (name, k))
def check_finite_arr(name, *arrs):
    for k, a in enumerate(arrs):
        if not np.all(np.isfinite(a)):
            raise FloatingPointError('%s[%d]: non-finite value' % (name, k))

def a2f(x):
    if not x.is_finite():
        raise FloatingPointError('non-finite arb ball %s' % x)
    return float(_dn(float(x.lower()))), float(_up(float(x.upper())))
def aup(x): return a2f(x)[1]
def ai(x):
    lo, hi = a2f(x); return I(lo, hi)

def gam(n):
    """gamma_{n} bound for a floating dot product of length n (plus one product rounding), with margin"""
    k = n + 2
    return 1.01 * k * UR / (1.0 - k * UR)
def to_mr(x):
    m = 0.5 * x.lo + 0.5 * x.hi
    r = np.maximum(_up(x.hi - m), _up(m - x.lo))
    return m, np.maximum(r, 0.0)
def from_mr(m, r):
    return I(_dn(m - r), _up(m + r))
def imatmul(A, B):
    """A: I (n, k) ; B: I (k, p) or exact ndarray -> I (n, p).  Rump: |fl(AmBm) - AmBm| <= gamma_k |Am||Bm|."""
    Am, Ar = to_mr(A)
    if isinstance(B, I):
        Bm, Br = to_mr(B)
    else:
        Bm, Br = np.asarray(B, float), None
    k = Am.shape[-1]; g = gam(k)
    aAm = np.abs(Am); aBm = np.abs(Bm)
    m = Am @ Bm
    r = g * (aAm @ aBm) + Ar @ aBm
    if Br is not None:
        r = r + (aAm + Ar) @ Br
    r = r * (1.0 + 4.0 * g) + 1e-300
    return from_mr(m, r)
def ieinsum(spec, A, B, n):
    """einsum of two I arrays; n = number of terms in each output sum."""
    Am, Ar = to_mr(A); Bm, Br = to_mr(B); g = gam(n)
    m = np.einsum(spec, Am, Bm)
    aAm = np.abs(Am); aBm = np.abs(Bm)
    r = g * np.einsum(spec, aAm, aBm) + np.einsum(spec, Ar, aBm) + np.einsum(spec, aAm + Ar, Br)
    r = r * (1.0 + 4.0 * g) + 1e-300
    return from_mr(m, r)
def isum(x, axis):
    m, r = to_mr(x); n = m.shape[axis]; g = gam(n)
    S = m.sum(axis); R = r.sum(axis) + g * np.abs(m).sum(axis)
    return from_mr(S, R * (1.0 + 4.0 * g) + 1e-300)
def upsum(a, axis=None):
    """upper bound of a sum of nonnegative floats"""
    a = np.asarray(a, float); n = a.size if axis is None else a.shape[axis]
    return a.sum(axis) * (1.0 + 4.0 * gam(n)) + 1e-300

# =====================================================================================================
# 2. lattice L = Hurwitz/2 (points l = m/4, m in Z^4 all even or all odd), lattice sums, Lemma A4
# =====================================================================================================
def L_points_m(Nmax, Nmin=1):
    M = int(math.isqrt(Nmax)) + 1
    out = []
    for par in (0, 1):
        vals = np.array([v for v in range(-M, M + 1) if v % 2 == par], dtype=np.int64)
        a, b, c = np.meshgrid(vals, vals, vals, indexing='ij'); a = a.ravel(); b = b.ravel(); c = c.ravel()
        q3 = a * a + b * b + c * c
        for m0 in vals:
            N = m0 * m0 + q3
            sel = (N >= Nmin) & (N <= Nmax)
            if sel.any():
                out.append(np.stack([np.full(int(sel.sum()), m0), a[sel], b[sel], c[sel]], 1))
    m = np.concatenate(out)
    return m, (m * m).sum(1)

NEAR_M, NEAR_N = L_points_m(16, 1)          # 0 < |l| <= 1 : 168 points
NEAR = NEAR_M.astype(float) / 4.0
MIN_VEC = NEAR[NEAR_N == 4]                 # the 24 minimal vectors (|l| = 1/2): V0 facets x.v <= 1/8

def a4_arb(T, p):
    """Lemma A4: sum_{|x-l| >= T} |x-l|^-p <= 1024 sum_j C(3,j) 4^(j-3) (T-1/2)^(j+1-p)/(p-j-1), T > 1/2, p > 4."""
    T = arb(T)
    assert float(T.lower()) > 0.5 and p > 4
    s = arb(0)
    for j in range(4):
        s += math.comb(3, j) * arb(4) ** (j - 3) * (T - arb(0.5)) ** (j + 1 - p) / (p - j - 1)
    return 1024 * s

class LatticeSums:
    """Z_full[p] = sum_{l != 0} |l|^-p, Z_far[p] = sum_{|l| > 1} |l|^-p, Z_tail[p] = sum_{|l| > RM} |l|^-p (arb),
       W_far(p, rho0) >= sum_{|l| > 1} (|l| - rho0)^-p  (upper floats, rho0 on a grid)."""
    def __init__(self, pmax=80):
        Nmax = int(16 * RM * RM)
        m, N = L_points_m(Nmax, 1)
        self.Ns, self.cnt = np.unique(N, return_counts=True)
        self.pmax = pmax
        self.Z_full = {}; self.Z_far = {}; self.Z_tail = {}
        for p in range(6, pmax + 1):
            tail = a4_arb(RM, p)
            fin_full = arb(0); fin_far = arb(0)
            for Nv, c in zip(self.Ns.tolist(), self.cnt.tolist()):
                t = c * (arb(16) / Nv) ** (arb(p) / 2)
                fin_full += t
                if Nv > 16: fin_far += t
            self.Z_full[p] = fin_full + arb(0).union(tail)
            self.Z_far[p] = fin_far + arb(0).union(tail)
            self.Z_tail[p] = arb(0).union(tail)
        self.rho_grid = np.round(np.arange(0, 0.46, 0.01), 2)
        self.Wfar = {}
        sq = np.sqrt(self.Ns.astype(float)) / 4.0
        far = self.Ns > 16
        for p in range(6, pmax + 1):
            row = []
            for rho in self.rho_grid:
                fin = float(np.sum(self.cnt[far] * (_dn(sq[far] - rho)) ** (-float(p))))
                fin = fin * (1 + 1e-12)
                tl = aup(a4_arb(RM - rho, p))
                row.append(_up(fin + tl))
            self.Wfar[p] = np.array(row)
    def wfar(self, p, rho0):
        """vectorized: upper bound of sum_{|l|>1} |c - l|^-p for |c| <= rho0 (array)"""
        k = np.ceil(np.asarray(rho0) / 0.01 - 1e-9).astype(int)
        if np.any(k >= len(self.rho_grid)): raise ValueError('rho0 outside table')
        return self.Wfar[p][np.maximum(k, 0)]

# =====================================================================================================
# 3. symmetry group G_J and Dirichlet cone (WP2a, reused)
# =====================================================================================================
def enumerate_W():
    mins = [tuple(2 * s * (np.arange(4) == i)) for i in range(4) for s in (1, -1)]
    mins += [tuple(v) for v in itertools.product((1, -1), repeat=4)]
    M = np.array(mins, dtype=np.int64)
    B = np.array([[2, 0, 0, 0], [0, 2, 0, 0], [0, 0, 2, 0], [1, 1, 1, 1]], dtype=np.int64)
    G = B @ B.T
    cand = [[]]
    for i in range(4):
        nxt = []
        for c in cand:
            for j in range(24):
                if M[j] @ M[j] != G[i, i]: continue
                if all(M[c[k]] @ M[j] == G[k, i] for k in range(i)):
                    nxt.append(c + [j])
        cand = nxt
    mats = []
    for c in cand:
        Vv = M[c]
        g = Vv.T.astype(float) @ np.linalg.inv(B.T.astype(float))
        g = np.round(g * 2) / 2
        assert np.allclose(g @ B.T, Vv.T) and np.allclose(g @ g.T, np.eye(4))
        mats.append(g)
    return mats
JMAT = np.array([[0, -1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, 1, 0]], float)
P_GEN = np.array([10.0, 4.0, 2.0, 1.0])
def group_GJ():
    W = enumerate_W()
    GJ = [g for g in W if np.allclose(g @ JMAT, JMAT @ g) or np.allclose(g @ JMAT, -JMAT @ g)]
    hol = sum(1 for g in GJ if np.allclose(g @ JMAT, JMAT @ g))
    return W, GJ, hol
def cone_normals(GJ):
    Ls = []
    for g in GJ:
        v = P_GEN - g @ P_GEN
        if np.any(v != 0): Ls.append(v)
    return np.array(Ls)
FREE = np.array([[j for j in range(4) if j != k] for k in range(4)])
def full_y(face, ylo, yhi):
    n = len(face); k = face // 2; sg = np.where(face % 2 == 0, 1.0, -1.0)
    Ylo = np.empty((n, 4)); Yhi = np.empty((n, 4)); r = np.arange(n)
    Ylo[r, k] = sg; Yhi[r, k] = sg
    fr = FREE[k]
    for a in range(3):
        Ylo[r, fr[:, a]] = ylo[:, a]; Yhi[r, fr[:, a]] = yhi[:, a]
    return Ylo, Yhi
def cone_keep_dirs(Lm, face, ylo, yhi):
    n = len(face); keep = np.ones(n, bool)
    Ylo, Yhi = full_y(face, ylo, yhi)
    for s in range(0, n, 20000):
        a = Ylo[s:s + 20000][:, None, :] * Lm[None]; b = Yhi[s:s + 20000][:, None, :] * Lm[None]
        mx = _up(np.maximum(a, b).sum(-1))
        keep[s:s + 20000] = np.all(mx >= 0, axis=1)
    return keep
def x_tight(Ylo, Yhi):
    """exact-range enclosure of x_i = y_i/|y| over the box (WP2a)."""
    mn = np.where((Ylo <= 0) & (Yhi >= 0), 0.0, np.minimum(Ylo * Ylo, Yhi * Yhi))
    mx = np.maximum(Ylo * Ylo, Yhi * Yhi)
    xs = []
    for i in range(4):
        o = [j for j in range(4) if j != i]
        Rmin = I(_dn(mn[:, o].sum(1)), _up(mn[:, o].sum(1)))
        Rmax = I(_dn(mx[:, o].sum(1)), _up(mx[:, o].sum(1)))
        yl, yh = Ylo[:, i], Yhi[:, i]
        Rl = iwhere(yl < 0, Rmin, Rmax); Rh = iwhere(yh >= 0, Rmin, Rmax)
        fl = I(yl) * irecip(isqrt(I(yl) * I(yl) + Rl))
        fh = I(yh) * irecip(isqrt(I(yh) * I(yh) + Rh))
        xs.append(I(fl.lo, fh.hi))
    return xs
def initial_dir_boxes(Lm, ny=8):
    g = np.linspace(-1, 1, ny + 1)
    cells = np.array(list(itertools.product(range(ny), repeat=3)))
    ylo = g[cells]; yhi = g[cells + 1]
    faces = np.repeat(np.arange(8), len(cells)); ylo = np.tile(ylo, (8, 1)); yhi = np.tile(yhi, (8, 1))
    keep = cone_keep_dirs(Lm, faces, ylo, yhi)
    return faces[keep], ylo[keep], yhi[keep]

# =====================================================================================================
# 4. c6 (WP2a code, reused): Ewald at eta and direct sum; enclosures intersected after an overlap check
# =====================================================================================================
def L_shell_groups(R):
    N_max = 16 * R * R
    out = {}
    M = int(4 * R) + 2
    for par in (0, 1):
        vals = np.array([v for v in range(-M, M + 1) if v % 2 == par])
        a, b, c = np.meshgrid(vals, vals, vals, indexing='ij')
        rest = (a * a + b * b + c * c).ravel()
        cnt = np.bincount(rest)
        for m0 in vals:
            q0 = int(m0 * m0)
            lim = N_max - q0
            for Nr in np.nonzero(cnt)[0]:
                if Nr >= lim: break
                N = q0 + int(Nr)
                if N == 0: continue
                out[(q0, N)] = out.get((q0, N), 0) + int(cnt[Nr])
    return out
def Jb(k, R, beta):
    den = 2 * beta - max(k - 1, 0) / (R * R)
    assert den > 0
    return R ** (k - 1) * (-beta * R * R).exp() / den
def c6_ewald(eta, T=3.5, K=10.0):
    eta_a = arb(eta); pi = arb.pi(); V = arb(1) / 32
    total = eta_a ** 4 / 24
    grp = L_shell_groups(T)
    fact = [arb(math.factorial(i)) for i in range(8)]
    for (q0, N), cnt in grp.items():
        Rl = arb(fmpq(N, 16)); l0sq = arb(fmpq(q0, 16))
        e = (-eta_a * Rl).exp()
        def Im(m):
            return fact[m] * e * sum((eta_a * Rl) ** i / fact[i] for i in range(m + 1)) / Rl ** (m + 1)
        cl = (64 * l0sq ** 3 * Im(6) - 480 * l0sq ** 2 * Im(5) + 720 * l0sq * Im(4) - 120 * Im(3)) / 720
        total += cnt * cl
    M = int(K / 2) + 1
    rng = np.arange(-M, M + 1)
    n = np.array(np.meshgrid(rng, rng, rng, rng, indexing='ij')).reshape(4, -1).T
    n = n[(n.sum(1) % 2 == 0)]
    k2 = 4 * (n * n).sum(1); k0sq = 4 * n[:, 0] ** 2
    sel = (k2 > 0) & (k2 < K * K)
    keys, cnts = np.unique(np.stack([k0sq[sel], k2[sel]], 1), axis=0, return_counts=True)
    rec = arb(0)
    for (a, b), c in zip(keys.tolist(), cnts.tolist()):
        rec += c * (-pi * pi * b / eta_a).exp() / b * (2 * pi) ** 6 * arb(a) ** 3
    total += -rec / (V * 720)
    hab = {0: 120, 2: 720, 4: 480, 6: 64}
    cq = {}
    for j, a in hab.items():
        m = (6 + j) // 2
        for i in range(m + 1):
            q = 2 * i - 8
            cq[q] = cq.get(q, arb(0)) + a * fact[m] / fact[i] * eta_a ** i
    Rr = arb(T) - arb(0.5)
    assert float(Rr.lower()) ** 2 > 3.0 / eta
    real_tail = 1024 * sum(c * sum(math.comb(3, jj) * arb(4) ** (jj - 3) * Jb(q + jj, Rr, eta_a) for jj in range(4))
                           for q, c in cq.items()) / 720
    Rk = arb(K) - 2 * arb(2).sqrt(); bet = pi * pi / eta_a
    assert float((Rk * Rk).lower()) >= 4 * eta / (2 * math.pi ** 2) and float((2 * bet * Rk * Rk).lower()) > 6
    rec_tail = (2 * pi) ** 6 / (V * 720) * sum(math.comb(3, jj) * arb(2).sqrt() ** (3 - jj) * Jb(4 + jj, Rk, bet) for jj in range(4))
    P6 = total + arb(0, 1) * (real_tail + rec_tail)
    return P6 * 8 / 3
def c6_certified():
    e6 = c6_ewald(6.0); e8 = c6_ewald(8.0)
    if not e6.overlaps(e8): raise FloatingPointError('c6: Ewald enclosures at eta 6 and 8 disagree')
    return e6.intersection(e8), e6, e8

# =====================================================================================================
# 5. multivariate polynomials in 4 real variables
# =====================================================================================================
def monomial_list(D):
    out = []
    for d in range(D + 1):
        for a in range(d, -1, -1):
            for b in range(d - a, -1, -1):
                for c in range(d - a - b, -1, -1):
                    out.append((a, b, c, d - a - b - c))
    return out
MON = monomial_list(DF1)                      # global ordering by degree; first nmon(K) entries = degree <= K
MON_E = np.array(MON, dtype=np.int64)
MON_DEG = MON_E.sum(1)
MON_IDX = {m: i for i, m in enumerate(MON)}
def nmon(D): return math.comb(D + 4, 4)
IDX4 = -np.ones((DF1 + 1,) * 4, dtype=np.int64)
IDX4[MON_E[:, 0], MON_E[:, 1], MON_E[:, 2], MON_E[:, 3]] = np.arange(len(MON))
FACT = [math.factorial(k) for k in range(40)]
def mfact(e): return FACT[e[0]] * FACT[e[1]] * FACT[e[2]] * FACT[e[3]]

def monomial_values(x, D, mr=False):
    """x: float array (n,4) (exact) -> x^gamma for |gamma| <= D: I (n, nmon(D)), or (mid, rad) if mr
       (relative error <= (D+3) u per entry: D-1 + 3 roundings)."""
    n = len(x); E = MON_E[:nmon(D)]
    P = [np.ones((n, D + 1)) for _ in range(4)]
    for i in range(4):
        for k in range(1, D + 1):
            P[i][:, k] = P[i][:, k - 1] * x[:, i]
    v = P[0][:, E[:, 0]] * P[1][:, E[:, 1]] * P[2][:, E[:, 2]] * P[3][:, E[:, 3]]
    r = np.abs(v) * ((D + 4) * UR * 1.01) + 1e-300
    if mr: return v, r
    return I(_dn(v - r), _up(v + r))
def monomial_ranges(xs, D):
    """xs: list of 4 I (n,) coordinate ranges -> I (n, nmon(D)) of ranges of x^gamma (natural, exact per monomial)."""
    n = len(xs[0].lo); E = MON_E[:nmon(D)]
    PW = []
    for i in range(4):
        cols_lo = [np.ones(n)]; cols_hi = [np.ones(n)]
        for k in range(1, D + 1):
            lo, hi = xs[i].lo, xs[i].hi
            a = lo ** k; b = hi ** k
            if k % 2 == 1:
                plo, phi = np.minimum(a, b), np.maximum(a, b)
            else:
                z = (lo <= 0) & (hi >= 0)
                plo = np.where(z, 0.0, np.minimum(a, b)); phi = np.maximum(a, b)
            r = (k + 1) * UR * 1.01
            cols_lo.append(_dn(plo - np.abs(plo) * r)); cols_hi.append(_up(phi + np.abs(phi) * r))
        PW.append(I(np.stack(cols_lo, 1), np.stack(cols_hi, 1)))
    out = PW[0][:, E[:, 0]] * PW[1][:, E[:, 1]]
    out = out * PW[2][:, E[:, 2]]
    out = out * PW[3][:, E[:, 3]]
    return out

def poly_eval(coef, x):
    """coef: I (nmon(D),) ; x: (n,4) floats -> I (n,)"""
    D = int(MON_DEG[len(coef.lo) - 1])
    return imatmul(monomial_values(x, D), I(coef.lo[:, None], coef.hi[:, None]))[:, 0]

_DMC = {}
def deriv_maps(Kfrom, beta_max):
    key = (Kfrom, beta_max)
    if key not in _DMC: _DMC[key] = _deriv_maps(Kfrom, beta_max)
    return _DMC[key]
def _deriv_maps(Kfrom, beta_max):
    """for every multi-index beta with |beta| <= beta_max: arrays (src, dst, factor) mapping Taylor coefficients of
       degree <= Kfrom to those of d^beta (factor = gamma!/(gamma-beta)!)."""
    betas = [b for b in MON if sum(b) <= beta_max]
    maps = {}
    nK = nmon(Kfrom)
    for b in betas:
        src = []; dst = []; fac = []
        for i in range(nK):
            g = MON[i]
            if all(g[k] >= b[k] for k in range(4)):
                gm = tuple(g[k] - b[k] for k in range(4))
                src.append(i); dst.append(MON_IDX[gm])
                f = 1
                for k in range(4):
                    f *= FACT[g[k]] // FACT[gm[k]]
                fac.append(float(f))
        maps[b] = (np.array(src, np.int64), np.array(dst, np.int64), np.array(fac))
    return betas, maps

# ---- Gegenbauer coefficients and lattice moment sums --------------------------------------------------
def gegen_c(nu, d, j):
    """c^(nu)_{d,j}: G_d^(nu)(x;e) = sum_j c_{d,j} (x.e)^(d-2j) |x|^(2j)  (exact rational as Fraction)."""
    from fractions import Fraction as Fr
    rising = 1
    for t in range(d - j): rising *= (nu + t)
    return Fr((-1) ** j * 2 ** (d - 2 * j) * rising, FACT[j] * FACT[d - 2 * j])

def even_alphas(D):
    return [a for a in MON if sum(a) <= D and all(x % 2 == 0 for x in a)]

class Moments:
    """mu(alpha, p) = sum_{l in S} l^alpha |l|^-p for a set S of shells (S = near: 4 <= |m|^2 <= 16, or
       far: 16 < |m|^2 <= 16 RM^2); even alpha only (odd moments vanish by the sign symmetries of L)."""
    def __init__(self, which, D):
        if which == 'near':
            m, N = NEAR_M, NEAR_N
        else:
            m, N = L_points_m(int(16 * RM * RM), 17)
        self.Ns, inv = np.unique(N, return_inverse=True)
        sorted_alphas = sorted({tuple(sorted(a, reverse=True)) for a in even_alphas(D)})
        mf = m.astype(float)
        pw = [[mf[:, i] ** k for k in range(D + 1)] for i in range(4)]
        self.PS = {}
        for a in sorted_alphas:
            v = pw[0][a[0]] * pw[1][a[1]] * pw[2][a[2]] * pw[3][a[3]]
            self.PS[a] = np.bincount(inv, weights=v, minlength=len(self.Ns))
        self.cache = {}
    def mu(self, alpha, p):
        """I enclosure of sum_l l^alpha |l|^-p (l = m/4)."""
        a = tuple(sorted(alpha, reverse=True)); key = (a, p)
        if key in self.cache: return self.cache[key]
        ps = self.PS[a]
        k = sum(a)
        # weight per shell: 4^-|alpha| (16/N)^(p/2), p even; computed by repeated multiplication
        w = (16.0 / self.Ns.astype(float))
        wp = np.ones_like(w)
        for _ in range(p // 2): wp = wp * w
        if p % 2: raise ValueError('odd p')
        rel = (p + k + 30) * UR * 1.01 + 1e-13          # PS nonneg sums: rel. error <= ~1e-13 (n <= 1e6 terms)
        val = float(np.sum(ps * wp)) * 0.25 ** k
        out = I(_dn(val * (1 - rel)), _up(val * (1 + rel)))
        self.cache[key] = out
        return out

def gegen_poly(nu, d, mom):
    """coefficients (I over all monomials of degree <= d, nonzero only at degree d) of
       sum_{l in S} G_d^(nu)(x; lhat) |l|^(-2nu-d)."""
    coef_lo = np.zeros(nmon(d)); coef_hi = np.zeros(nmon(d))
    gams = [g for g in MON if sum(g) == d and all(x % 2 == 0 for x in g)]
    rows = []; vals = []
    for g in gams:
        terms = []
        for j in range(d // 2 + 1):
            c = gegen_c(nu, d, j); m = d - 2 * j; p = 2 * nu + 2 * d - 2 * j
            for b in itertools.product(*[range(gk // 2 + 1) for gk in g]):
                if sum(b) != j: continue
                al = tuple(g[k] - 2 * b[k] for k in range(4))
                from fractions import Fraction as Fr
                w = Fr(c) * Fr(FACT[m], mfact(al)) * Fr(FACT[j], mfact(b))
                terms.append((w, al, p))
        acc = I(0.0)
        for (w, al, p) in terms:
            acc = acc + I.c(w) * mom.mu(al, p)
        i = MON_IDX[g]
        coef_lo[i] = acc.lo; coef_hi[i] = acc.hi
    return I(coef_lo, coef_hi)

def H6_poly():
    """H6 = 6 m3 - (15/2) p2 m2 + (15/8) p2^3 as exact coefficients (degree-6 monomials)."""
    from collections import defaultdict
    c = defaultdict(lambda: 0)
    # m3 = sum x_i^6 ; p2 m2 = sum_i sum_j x_i^2 x_j^4 ; p2^3 = (sum x_i^2)^3
    for i in range(4):
        e = [0] * 4; e[i] = 6; c[tuple(e)] += 6
    for i in range(4):
        for j in range(4):
            e = [0] * 4; e[i] += 2; e[j] += 4; c[tuple(e)] += -7.5
    for i in range(4):
        for j in range(4):
            for k in range(4):
                e = [0] * 4; e[i] += 2; e[j] += 2; e[k] += 2; c[tuple(e)] += 1.875
    v = np.zeros(nmon(6))
    for e, val in c.items(): v[MON_IDX[e]] = val
    return v
def s_poly(D):
    v = np.zeros(nmon(D))
    for i in range(4):
        e = [0] * 4; e[i] = 2; v[MON_IDX[tuple(e)]] = 1.0
    return v

def pad(p, D):
    """I polynomial padded to nmon(D)"""
    n = nmon(D); lo = np.zeros(n); hi = np.zeros(n)
    k = min(len(p.lo), n)
    lo[:k] = p.lo[:k]; hi[:k] = p.hi[:k]
    return I(lo, hi)

# =====================================================================================================
# 6. explicit tail constants (A2 + Cauchy / homogeneity)
# =====================================================================================================
_BDM = {}
def Bdm(d, m):
    """B(d,m) >= max_{|beta|=m} |d^beta P(x)| / (C |x|^(d-m)) for P homogeneous of degree d with
       |P(y)| <= C Lie(y)^d on C^4 (x real): m! (2 sqrt2 (d-m)/m)^m (d/(d-m))^d ; m = 0: 1 ; m = d: d! (2 sqrt 2)^d.
       (Cauchy on the polydisc of radius rho = m|x|/(2 sqrt2 (d-m)) around x; Lie(x + w) <= |x| + 2 sqrt2 rho.)
       Evaluated in arb, returned rounded up."""
    key = (d, m)
    if key in _BDM: return _BDM[key]
    if m == 0: v = 1.0
    elif m > d: v = 0.0
    elif m == d: v = aup(arb(FACT[d]) * (2 * arb(2).sqrt()) ** d)
    else: v = aup(arb(FACT[m]) * (2 * arb(2).sqrt() * (d - m) / m) ** m * (arb(d) / (d - m)) ** d)
    _BDM[key] = v
    return v
def binom_nu(d, nu):
    """|G_d^(nu)| <= binom(d+2nu-1, d) Lie^d"""
    return math.comb(d + 2 * nu - 1, d)

def gegen_tail_coeffs(nu, m, d0, nterms=60):
    """coefficients for sum_{d >= d0} B(d,m) binom(d+2nu-1,d) t^(d-m) (explicit part) and the geometric tail
       parameters: returns (coef array over d = d0..d0+nterms-1, crude tail factor function)."""
    cs = np.array([Bdm(d, m) * binom_nu(d, nu) for d in range(d0, d0 + nterms)])
    return cs
def gegen_tail(t, nu, m, d0, nterms=60):
    """upper bound (float array) of sum_{d >= d0} B(d,m) binom(d+2nu-1,d) t^(d-m), for 0 <= t <= 1/2.
       Beyond d0+nterms: B(d,m) <= m!(2 sqrt2 e^2/m)^m d^m (d >= 2m) and the term ratio
       ((d+1)/d)^m (d+2nu)/(d+1) t is decreasing in d, so a geometric bound applies."""
    t = np.asarray(t, float)
    if np.any(t > 0.5): raise ValueError('gegen_tail: t > 1/2')
    cs = gegen_tail_coeffs(nu, m, d0, nterms)
    s = np.zeros_like(t)
    for k in range(nterms - 1, -1, -1):
        s = s * t + cs[k]
    s = s * upow(t, d0 - m)
    dl = d0 + nterms
    if m == 0:
        crude = binom_nu(dl, nu)
    else:
        crude = aup(arb(FACT[m]) * (2 * arb(2).sqrt() * arb(1).exp() ** 2 / m) ** m * arb(dl) ** m * binom_nu(dl, nu))
    q = ((dl + 1) / dl) ** m * (dl + 2 * nu) / (dl + 1) * t
    rest = crude * upow(t, dl - m) / (1 - q)
    return (s + rest) * (1 + 1e-10)

# =====================================================================================================
# 7. the polynomial data (built once, cached): neck polynomials, far polynomials, far T-matrices
# =====================================================================================================
class Data:
    def __init__(self, log=print, c6=None):
        t0 = time.time()
        self.LS = LatticeSums(pmax=2 * 3 + DF1 + 2 * 3 + 70)
        log('  lattice sums (Z, W tables) %.1f s' % (time.time() - t0))
        if c6 is None:
            c6, _, _ = c6_certified()
        self.c6_arb = c6; self.c6 = ai(c6)
        log('  c6 = %s (%.1f s)' % (c6.str(20), time.time() - t0))
        self.mom_near = Moments('near', DF1)
        self.mom_far = Moments('far', DF1)
        log('  moment tables %.1f s' % (time.time() - t0))
        # near and far homogeneous parts
        self.Pnear = {1: {}, 3: {}}; self.Pfar = {1: {}, 3: {}}
        for d in range(0, DF1 + 1, 2):
            self.Pnear[1][d] = gegen_poly(1, d, self.mom_near)
            if d >= 8: self.Pfar[1][d] = gegen_poly(1, d, self.mom_far)
        for d in range(0, DF3 + 1, 2):
            self.Pnear[3][d] = gegen_poly(3, d, self.mom_near)
            if d >= 4: self.Pfar[3][d] = gegen_poly(3, d, self.mom_far)
        # cross-check of the moment machinery on d = 6: far moments (with RM tail) vs c6 H6 - near6
        self.pi2 = ai(arb.pi() ** 2)
        H6 = H6_poly()
        P6far_c6 = pad(I(H6) * self.c6, 6) - self.Pnear[1][6]
        P6far_mom = gegen_poly(1, 6, self.mom_far)
        tail6 = binom_nu(6, 1) * aup(self.LS.Z_tail[8])          # coefficientwise bound is not needed; see check
        self.check_P6 = {'max_abs_diff_coef': float(np.max(np.abs(0.5 * (P6far_c6.lo + P6far_c6.hi) - 0.5 * (P6far_mom.lo + P6far_mom.hi)))),
                         'RM_tail_sup_bound_on_sphere': tail6}
        self.Pfar[1][6] = P6far_c6
        self.Pfar[1][4] = I(np.zeros(nmon(4)))
        P2 = np.zeros(nmon(2)); P2s = s_poly(2)
        self.Pfar[1][2] = I(P2s) * (16.0 * self.pi2)                # P_2[S_far] = 16 pi^2 |x|^2 (L3 (c))
        self.Pfar[3][2] = -self.Pnear[3][2]                          # so that P_2[psi_2] = 0 exactly (L3 (c))
        # full (neck) parts:  psi1 = sum_d P1[d], psi2 = sum_d P2[d]
        self.P1 = {6: pad(I(H6) * (self.c6 * -0.5), 6)}
        for d in range(8, DN1 + 1, 2):
            self.P1[d] = (self.Pnear[1][d] + self.Pfar[1][d]) * -0.5
        self.Pnear[5] = {}; self.Pfar[5] = {}; self.P3 = {}
        for d in range(4, DN3 + 1, 2):
            self.Pnear[5][d] = gegen_poly(5, d, self.mom_near); self.Pfar[5][d] = gegen_poly(5, d, self.mom_far)
            self.P3[d] = (self.Pnear[5][d] + self.Pfar[5][d]) * (-1.0 / 80.0)     # psi_3 = (1/16)[-(1/5) sum |z-l|^-10 + 2 Z_3 |z|^2]
        self.P2 = {}
        for d in range(4, DN3 + 1, 2):
            self.P2[d] = (self.Pnear[3][d] + self.Pfar[3][d]) * (1.0 / 24.0)
        log('  polynomials %.1f s' % (time.time() - t0))
        # far polynomials (sum over degrees) for the bulk sources: F1 = S_far (nu = 1), F3 (nu = 3, P_2 adjusted)
        self.F1 = pad(self.Pfar[1][2], DF1)
        for d in range(4, DF1 + 1, 2): self.F1 = self.F1 + pad(self.Pfar[1][d], DF1)
        self.F3 = pad(self.Pfar[3][2], DF3)
        for d in range(4, DF3 + 1, 2): self.F3 = self.F3 + pad(self.Pfar[3][d], DF3)
        self.T1 = self.shift_matrix(self.F1, DF1)
        self.T3 = self.shift_matrix(self.F3, DF3)
        log('  shift matrices %.1f s' % (time.time() - t0))
        # tail constants: C_d bounds |P_d(y)| <= C_d Lie(y)^d
        Z = lambda tab, p: aup(tab[p]) if p >= 6 else float('nan')
        self.Cfar_D = {nu: {d: binom_nu(d, nu) * Z(self.LS.Z_far, 2 * nu + d) for d in range(4, 99)} for nu in (1, 3)}
        self.Cfar_R = {1: {d: (binom_nu(d, 1) * Z(self.LS.Z_tail, 2 + d) if d >= 8 else 0.0) for d in range(0, 99)},
                       3: {d: (binom_nu(d, 3) * Z(self.LS.Z_tail, 6 + d) if d >= 4 else 0.0) for d in range(0, 99)}}
        self.Cfull = {nu: {d: binom_nu(d, nu) * Z(self.LS.Z_full, 2 * nu + d) for d in range(4, 99)} for nu in (1, 3)}
        self.Cfull[5] = {d: binom_nu(d, 5) * Z(self.LS.Z_full, 10 + d) / 80.0 for d in range(4, 95)}
        self.Cfar_R[5] = {d: binom_nu(d, 5) * Z(self.LS.Z_tail, 10 + d) / 80.0 for d in range(4, 95)}
        # R-tail for the neck polynomials (moment sums over |l| <= RM), degrees >= 8 (nu=1) and >= 4 (nu=3)
        self.Z3 = ai(self.LS.Z_full[12]); self.Z2 = ai(self.LS.Z_full[8])
        self.seconds = time.time() - t0

    def shift_matrix(self, F, D):
        """T[mu, beta] = coef_{mu+beta} binom(mu+beta, beta): Taylor coefficients at c = C(c) @ T, C = (c^mu)."""
        nb = nmon(KP); nm = nmon(D)
        lo = np.zeros((nm, nb)); hi = np.zeros((nm, nb))
        E = MON_E[:nm]
        for jb in range(nb):
            b = MON[jb]
            G = E + np.array(b)
            ok = G.sum(1) <= D
            gi = np.where(ok, IDX4[np.minimum(G[:, 0], DF1), np.minimum(G[:, 1], DF1), np.minimum(G[:, 2], DF1), np.minimum(G[:, 3], DF1)], -1)
            ok = ok & (gi >= 0)
            binv = np.ones(nm)
            for k in range(4):
                for t in range(b[k]):
                    binv *= (G[:, k] - t)
                binv /= FACT[b[k]]
            binv = np.where(ok, binv, 0.0)          # exact integers (< 2^53)
            ci_lo = np.where(ok, F.lo[np.maximum(gi, 0)], 0.0); ci_hi = np.where(ok, F.hi[np.maximum(gi, 0)], 0.0)
            prod = I(ci_lo, ci_hi) * binv
            lo[:, jb] = prod.lo; hi[:, jb] = prod.hi
        return I(lo, hi)


# =====================================================================================================
# 8. bulk sources: Taylor coefficients at box centres (order KP) + remainder bounds for derivatives 0..4
# =====================================================================================================
_PJ = {}
def _pj_tables():
    """per degree p: alpha indices (exact degree p), parent index (degree p-1) and variable; matrices M_{p,j}."""
    if _PJ: return _PJ
    for p in range(KP + 1):
        al = [i for i in range(nmon(KP)) if MON_DEG[i] == p]
        par = []; var = []
        if p > 0:
            prev = [i for i in range(nmon(KP)) if MON_DEG[i] == p - 1]
            ppos = {MON[i]: k for k, i in enumerate(prev)}
            for i in al:
                e = MON[i]; v = next(k for k in range(4) if e[k] > 0)
                pe = list(e); pe[v] -= 1
                par.append(ppos[tuple(pe)]); var.append(v)
        _PJ[('al', p)] = al; _PJ[('par', p)] = np.array(par, np.int64); _PJ[('var', p)] = np.array(var, np.int64)
        for j in range((KP - p) // 2 + 1):
            M = np.zeros((len(al), nmon(KP)))
            for k, i in enumerate(al):
                a = MON[i]
                for b in MON:
                    if sum(b) != j: continue
                    g = tuple(a[t] + 2 * b[t] for t in range(4))
                    M[k, MON_IDX[g]] += float(FACT[p] // mfact(a) * (FACT[j] // mfact(b)))
                    if sum(b) > j: break
            _PJ[('M', p, j)] = M
    return _PJ

def kappa_pj(nu, p, j):
    from fractions import Fraction as Fr
    rising = 1
    for t in range(p + j): rising *= (nu + t)
    return Fr((-1) ** (p + j) * 2 ** p * rising, FACT[p] * FACT[j])

def upow(x, k):
    """upper bound of x^k for x >= 0 (float arrays), by k-1 multiplications rounded up"""
    out = np.ones_like(np.asarray(x, float))
    for _ in range(k): out = _up(out * x)
    return out

def mr_matmul(Am, Ar, B):
    """(Am +- Ar) @ B with B an I matrix or exact ndarray -> I"""
    if isinstance(B, I):
        Bm, Br = to_mr(B)
    else:
        Bm, Br = np.asarray(B, float), None
    k = Am.shape[-1]; g = gam(k)
    aAm = np.abs(Am); aBm = np.abs(Bm)
    m = Am @ Bm
    r = g * (aAm @ aBm) + Ar @ aBm
    if Br is not None: r = r + (aAm + Ar) @ Br
    return from_mr(m, r * (1.0 + 4.0 * g) + 1e-300)

def near_taylor(c, pts, nu, Kp, dmax):
    """Taylor coefficients (degree <= Kp) at c of sum_{l in pts} |c + delta - l|^(-2 nu), and remainder bounds
       rem[:, m] >= sup_{|delta| <= dmax} max_{|beta| = m} |d^beta (f - T_Kp)|, m = 0..4 (Gegenbauer tails, A2 +
       homogeneity-Cauchy).  Requires dmax <= |c - l|/2 for all l.
       Arithmetic: plain floats with a running relative error bound (first order, x 1.01), sums in Rump form:
       y = fl(c - l) (rel. err. u), y^alpha (2p u), A = |y|^2 (7u), A^-k (9k u), kappa (u), products (u)."""
    T = _pj_tables()
    n = len(c); m = len(pts)
    y = c[:, None, :] - pts[None, :, :]
    A = (y * y).sum(2)
    if np.any(A <= 0): raise FloatingPointError('near_taylor: centre on a lattice point')
    Ainv = 1.0 / A
    pw = {0: np.ones((n, m))}
    for k in range(1, nu + Kp + 1): pw[k] = pw[k - 1] * Ainv
    cm = np.zeros((n, nmon(Kp))); cr = np.zeros((n, nmon(Kp)))
    Y = np.ones((n, m, 1))
    g = gam(m)
    for p in range(Kp + 1):
        if p > 0:
            Y = Y[:, :, T[('par', p)]] * y[:, :, T[('var', p)]]
        aY = np.abs(Y)
        for j in range((Kp - p) // 2 + 1):
            k = nu + p + j
            kap = float(kappa_pj(nu, p, j))
            w = pw[k] * kap
            rel = (9 * k + 2 * p + 4) * UR * 1.01
            Wm = np.einsum('nl,nla->na', w, Y)
            Wr = np.einsum('nl,nla->na', np.abs(w), aY) * (rel + g) * 1.01 + 1e-300
            Mpj = T[('M', p, j)][:, :nmon(Kp)]
            G = gam(Mpj.shape[0])
            cm = cm + Wm @ Mpj
            cr = cr + (Wr + G * np.abs(Wm)) @ np.abs(Mpj) * (1 + 4 * G) + UR * np.abs(cm) * 2
    coef = from_mr(cm, cr * (1 + 1e-12) + 1e-300)
    ylo = np.sqrt(A) * (1 - 16 * UR)
    t = _up(dmax[:, None] / ylo)
    if np.any(t > 0.5): raise ValueError('near_taylor: box too large for the Gegenbauer remainder')
    rem = np.zeros((n, 5))
    yinv = _up(1.0 / ylo)
    for mm in range(5):
        tl = gegen_tail(t, nu, mm, Kp + 1)
        rem[:, mm] = upsum(tl * upow(yinv, 2 * nu + mm), axis=1)
    return coef, rem

def far_coef_tail(Dd, nu, rho, m):
    """bound for |Taylor coefficient of order m| at a real point |c| <= rho of the far function minus its
       computed polynomial: degrees > DF (Z_far) and the lattice tail |l| > RM of the computed degrees."""
    DF = DF1 if nu == 1 else DF3
    tot = np.zeros_like(rho)
    dlast = 98 - 2 * nu
    for d in range(4, dlast + 1, 2):
        Cd = (Dd.Cfar_D[nu][d] if d > DF else Dd.Cfar_R[nu][d])
        if Cd == 0.0: continue
        if m > d: continue
        tot = tot + (Bdm(d, m) / FACT[m]) * Cd * upow(rho, d - m)
    # beyond dlast: ratio of successive terms <= ((d+2)/d)^(m+2nu-1) (rho/1.118)^2 (two degrees) * e^m-type growth;
    # crude geometric bound with ratio q <= 0.5 for rho <= 0.45
    q = ((dlast + 2) / dlast) ** (m + 2 * nu + 2) * (rho / 1.118) ** 2 * 1.2
    last = (Bdm(dlast, m) / FACT[m]) * Dd.Cfar_D[nu][dlast] * upow(rho, dlast - m)
    return _up((tot + last * q / (1 - q)) * (1 + 1e-10))

def far_taylor(Dd, c, nu, dmax, Cmon):
    """Taylor coefficients at c of the far function (nu = 1: S_far, nu = 3: F3) to order KP, coefficient tails
       added as intervals, remainder bounds rem[:, m] (m = 0..4) from the per-lattice-point Gegenbauer tails."""
    n = len(c)
    if nu == 1:
        coef = mr_matmul(Cmon[0], Cmon[1], Dd.T1)
    else:
        coef = mr_matmul(Cmon[0][:, :nmon(DF3)], Cmon[1][:, :nmon(DF3)], Dd.T3)
    rho = np.sqrt((c * c).sum(1)) * (1 + 1e-12)
    for mdeg in range(KP + 1):
        cols = np.where(MON_DEG[:nmon(KP)] == mdeg)[0]
        tb = far_coef_tail(Dd, nu, rho, mdeg)
        coef.lo[:, cols] = _dn(coef.lo[:, cols] - tb[:, None]); coef.hi[:, cols] = _up(coef.hi[:, cols] + tb[:, None])
    rem = np.zeros((n, 5))
    rho0 = rho + 1e-9
    for mm in range(5):
        s = np.zeros(n)
        for d in range(KP + 1, KP + 41):
            s = s + Bdm(d, mm) * binom_nu(d, nu) * upow(dmax, d - mm) * Dd.LS.wfar(2 * nu + d, rho0)
        d = KP + 41
        crude = Bdm(d, mm) * binom_nu(d, nu) * upow(dmax, d - mm) * Dd.LS.wfar(2 * nu + d, rho0)
        q = ((d + 1 - mm) / (d - mm)) ** mm * ((d + 2 * nu) / (d + 1)) * dmax / (1.118 - rho0) * (1 + 1e-9)
        if np.any(q >= 0.5): raise ValueError('far_taylor: box too large')
        rem[:, mm] = _up((s + crude / (1 - q)) * (1 + 1e-10))
    return coef, rem

def eps_tail_bounds(Dd, c, dmax, eps_hi, kmin=3):
    """B[:, m] >= sup over the box of max_{|beta| = m} |d^beta (psi~_2 - psi_2)|, m = 0..4, where
       psi~_2 - psi_2 = sum_{k>=3} eps^(k-2) psi_k (L3 (a), (h)); A6 (f) for the derivatives of |x|^(2-4k)."""
    n = len(c); B = np.zeros((n, 5))
    if np.all(eps_hi == 0): return B
    yv = c[:, None, :] - NEAR[None]
    ylo = np.sqrt((yv * yv).sum(2)) - dmax[:, None]
    if np.any(ylo <= 0.1): raise ValueError('eps_tail: box too close to a lattice point')
    rho = np.sqrt((c * c).sum(1)) + dmax
    for mm in range(5):
        tot = np.zeros(n)
        for k in range(kmin, 9):
            # b_k = binom(1/2, k)
            from fractions import Fraction as Fr
            bkf = Fr(1)
            for t in range(k): bkf *= (Fr(1, 2) - t)
            bkf /= FACT[k]
            bk = abs(float(bkf)) * (1 + 1e-15)
            nu2 = 2 * k - 1; p = 2 * nu2 + mm
            Msum = (ylo ** (-float(p))).sum(1) + Dd.LS.wfar(p, rho + 1e-9) if p <= Dd.LS.pmax else None
            der = FACT[mm] * math.comb(mm + 2 * nu2 - 1, mm) * (2 * math.sqrt(2)) ** mm
            term = eps_hi ** (k - 2) * bk * (der * Msum / (2 * k - 1))
            Zk = aup(Dd.LS.Z_full[4 * k]) if 4 * k <= Dd.LS.pmax else None
            quad = [rho ** 2, 2 * rho, 2.0 * np.ones(n), np.zeros(n), np.zeros(n)][mm]
            term = term + eps_hi ** (k - 2) * bk * (k - 1) * Zk * quad
            tot = tot + term
        # k > 8: eps |y|^-4 <= 0.03^4/0.15^4 < 1e-3 per step, binomial growth < 4: geometric with ratio < 0.01
        B[:, mm] = _up(tot * (1 + 0.02) + 1e-300)
    return B

def bulk_sources(Dd, c, h, eps_hi):
    """psi_1 and psi~_2 Taylor data at the box centres c (cubes of half-width h)."""
    dmax = 2.0 * h * (1 + 1e-12)
    n1, r1n = near_taylor(c, NEAR, 1, KP, dmax)
    n3, r3n = near_taylor(c, NEAR, 3, KP, dmax)
    Cmon = monomial_values(c, DF1, mr=True)
    f1, r1f = far_taylor(Dd, c, 1, dmax, Cmon)
    f3, r3f = far_taylor(Dd, c, 3, dmax, Cmon)
    S1 = n1 + f1; S3 = n3 + f3
    # s = |c|^2 + 2 c.delta + |delta|^2
    s_t = I(np.zeros((len(c), nmon(KP))))
    cc = I(c[:, 0]) * I(c[:, 0]) + I(c[:, 1]) * I(c[:, 1]) + I(c[:, 2]) * I(c[:, 2]) + I(c[:, 3]) * I(c[:, 3])
    s_t.lo[:, 0] = cc.lo; s_t.hi[:, 0] = cc.hi
    for i in range(4):
        e = [0] * 4; e[i] = 1; k = MON_IDX[tuple(e)]
        s_t.lo[:, k] = 2 * c[:, i]; s_t.hi[:, k] = 2 * c[:, i]
        e = [0] * 4; e[i] = 2; k = MON_IDX[tuple(e)]
        s_t.lo[:, k] = 1.0; s_t.hi[:, k] = 1.0
    a1 = s_t * (8.0 * Dd.pi2) - S1 * 0.5
    rem1 = _up((r1n + r1f) * 0.5 * (1 + 1e-12))
    a2 = S3 * (1.0 / 24.0)
    rem2 = _up((r3n + r3f) / 24.0 * (1 + 1e-12)) + eps_tail_bounds(Dd, c, dmax, eps_hi)
    check_finite('bulk_sources', a1, a2); check_finite_arr('bulk_sources_rem', rem1, rem2)
    return a1, rem1, a2, rem2

# =====================================================================================================
# 9. Taylor models (order KT = 4 in 4 real variables, cube boxes |delta_i| <= h)
# =====================================================================================================
NK = nmon(KT)
_TMT = {}
def _tm_tables():
    if _TMT: return _TMT
    deg = MON_DEG[:NK]
    PI_, PJ_, PK_ = [], [], []
    H = np.zeros((NK, NK))
    for i in range(NK):
        for j in range(NK):
            g = tuple(MON[i][k] + MON[j][k] for k in range(4))
            if sum(g) <= KT:
                PI_.append(i); PJ_.append(j); PK_.append(MON_IDX[g])
            else:
                H[i, j] = 1.0
    S = np.zeros((len(PI_), NK)); S[np.arange(len(PI_)), PK_] = 1.0
    even = np.array([all(x % 2 == 0 for x in MON[i]) for i in range(nmon(KP))])
    _TMT.update(dict(PI=np.array(PI_), PJ=np.array(PJ_), S=S, H=H, deg=deg, even=even))
    return _TMT
def hpow(h, D):
    """h^k for k = 0..D, (n, D+1), rounded up"""
    out = np.ones((len(h), D + 1))
    for k in range(1, D + 1): out[:, k] = _up(out[:, k - 1] * h)
    return out

class TM:
    __slots__ = ('coef', 'rem', 'hp')
    def __init__(s, coef, rem, hp):
        s.coef = coef; s.rem = np.asarray(rem, float); s.hp = hp
    def mags(s):
        return s.coef.mag() * s.hp[:, MON_DEG[:NK]]
    def bound(s):
        return upsum(s.mags(), 1) + s.rem
    def __add__(s, o):
        if isinstance(o, TM): return TM(s.coef + o.coef, _up(s.rem + o.rem), s.hp)
        c = I(s.coef.lo.copy(), s.coef.hi.copy())
        if isinstance(o, np.ndarray): o = I(o)
        if isinstance(o, I):
            c.lo[:, 0] = _dn(c.lo[:, 0] + o.lo); c.hi[:, 0] = _up(c.hi[:, 0] + o.hi)
        else:
            o = float(o); c.lo[:, 0] = _dn(c.lo[:, 0] + o); c.hi[:, 0] = _up(c.hi[:, 0] + o)
        return TM(c, s.rem, s.hp)
    __radd__ = __add__
    def __neg__(s): return TM(-s.coef, s.rem, s.hp)
    def __sub__(s, o): return s + (-o)
    def __rsub__(s, o): return (-s) + o
    def __mul__(s, o):
        if isinstance(o, TM):
            T = _tm_tables()
            P = s.coef[:, T['PI']] * o.coef[:, T['PJ']]
            coef = imatmul(P, T['S'])
            am = s.mags(); bm = o.mags()
            high = upsum((am @ T['H']) * bm, 1) * (1 + 4 * gam(NK))
            Ab = upsum(am, 1); Bb = upsum(bm, 1)
            rem = _up((high + Ab * o.rem + Bb * s.rem + s.rem * o.rem) * (1 + 1e-12))
            return TM(coef, rem, s.hp)
        if isinstance(o, I):
            return TM(s.coef * I(o.lo[:, None], o.hi[:, None]), _up(s.rem * o.mag()), s.hp)
        o = float(o)
        return TM(s.coef * o, _up(s.rem * abs(o)), s.hp)
    __rmul__ = __mul__
    def range(s):
        T = _tm_tables(); d = MON_DEG[:NK]
        hd = s.hp[:, d]
        a = s.coef
        ev = T['even'][:NK]
        lo = np.where(ev[None, :], np.minimum(0.0, a.lo * hd), -a.mag() * hd)
        hi = np.where(ev[None, :], np.maximum(0.0, a.hi * hd), a.mag() * hd)
        lo[:, 0] = a.lo[:, 0]; hi[:, 0] = a.hi[:, 0]
        S = isum(I(_dn(lo), _up(hi)), 1)
        return I(_dn(S.lo - s.rem), _up(S.hi + s.rem))
    def at_centre(s):
        return I(_dn(s.coef.lo[:, 0] - s.rem), _up(s.coef.hi[:, 0] + s.rem))

def poly_tail_bound(coef, hp, dfrom):
    """sum over monomials of degree > dfrom of |a| h^deg (upper)"""
    D = int(MON_DEG[coef.shape[1] - 1])
    sel = MON_DEG[:coef.shape[1]] > dfrom
    if not np.any(sel): return np.zeros(coef.shape[0])
    return upsum(coef.mag()[:, sel] * hp[:, MON_DEG[:coef.shape[1]][sel]], 1)
def tm_from_poly(coef, rem, hp):
    """TM of order KT from a polynomial (I coefficients, any degree) plus a remainder bound"""
    extra = poly_tail_bound(coef, hp, KT)
    return TM(coef[:, :NK], _up(rem + extra), hp)
def poly_range_sym(coef, hp):
    """enclosure of a polynomial over the box: a_0 +- sum_{gamma != 0} |a| h^deg"""
    t = poly_tail_bound(coef, hp, 0)
    return I(_dn(coef.lo[:, 0] - t), _up(coef.hi[:, 0] + t))

_DM = {}
def dmaps():
    if not _DM:
        betas, maps = deriv_maps(KP, 4)
        _DM['betas'] = betas; _DM['maps'] = maps
    return _DM
def deriv_poly(a, beta):
    src, dst, fac = dmaps()['maps'][beta]
    m = sum(beta)
    n = a.shape[0]
    lo = np.zeros((n, nmon(KP - m))); hi = np.zeros((n, nmon(KP - m)))
    v = a[:, src] * fac
    lo[:, dst] = v.lo; hi[:, dst] = v.hi
    return I(lo, hi)
def E2(a, b):
    e = [0, 0, 0, 0]; e[a] += 1; e[b] += 1; return tuple(e)
PAIRS = [(a, b) for a in range(4) for b in range(a, 4)]

# =====================================================================================================
# 10. the Q formula (generic: TM or Jet arithmetic)
# =====================================================================================================
def N_of(D, Aq):
    """N[psi] = (1/4) sum_ab (x_a x_b + (Jx)_a (Jx)_b) d_a d_b psi ; D[(a,b)] second derivatives, Aq[(a,b)] the quadratics"""
    tot = None
    for (a, b) in PAIRS:
        w = 0.25 if a == b else 0.5
        t = (Aq[(a, b)] * D[(a, b)]) * w
        tot = t if tot is None else tot + t
    return tot
def L_of(D):
    return (D[(0, 0)] + D[(1, 1)] + D[(2, 2)] + D[(3, 3)]) * 0.25
def K_of(D):
    """complex Hessian entries: k11, k22 (real), k12 = kre + i kim"""
    k11 = (D[(0, 0)] + D[(1, 1)]) * 0.25
    k22 = (D[(2, 2)] + D[(3, 3)]) * 0.25
    kre = (D[(0, 2)] + D[(1, 3)]) * 0.25
    kim = (D[(0, 3)] - D[(1, 2)]) * 0.25
    return k11, k22, kre, kim
def det_of(K):
    k11, k22, kre, kim = K
    return k11 * k22 - kre * kre - kim * kim
def quad_A(X):
    JX = [-X[1], X[0], -X[3], X[2]]
    return {(a, b): X[a] * X[b] + JX[a] * JX[b] for (a, b) in PAIRS}
def Q_formula(D1, D2, Aq, sinv3, theta, eps, det_scale=None, sinv4=None):
    """Q = theta ((N[psi1] + eps N[psi2]) s^-3 + Lap psi2) + det(K1 + eps K2)   (L4 (d), rewritten with
       (1/theta - theta) = u theta).  Mutations (controls only):
         drop_N1        : omit the N[psi1]/s^3 term of the cancelled form
         drop_l0_cancel : keep the l = 0 term (1/4) s^-4 of tr d dbar Psi_2 without its partner det H_0 (L4 (e))"""
    N1 = N_of(D1, Aq); N2 = N_of(D2, Aq); L2 = L_of(D2)
    Mt = {k: D1[k] + D2[k] * eps for k in PAIRS}
    det = det_of(K_of(Mt))
    if det_scale is not None: det = det * det_scale
    inner = N2 * eps
    if not MUT['drop_N1']: inner = N1 + inner
    if MUT['drop_l0_cancel'] and sinv4 is not None: L2 = L2 + sinv4 * 0.25
    Q = theta * (inner * sinv3 + L2) + det
    return Q, Mt

# =====================================================================================================
# 11. interval 2-jets (value, gradient, Hessian in the 4 real coordinates)
# =====================================================================================================
class Jet:
    __slots__ = ('v', 'g', 'H')
    def __init__(s, v, g, H): s.v = v; s.g = g; s.H = H
    @staticmethod
    def const(v, n):
        z = I(np.zeros(n)); return Jet(v, [z] * 4, {k: z for k in PAIRS})
    def __add__(s, o):
        if isinstance(o, Jet): return Jet(s.v + o.v, [a + b for a, b in zip(s.g, o.g)], {k: s.H[k] + o.H[k] for k in PAIRS})
        return Jet(s.v + o, s.g, s.H)
    __radd__ = __add__
    def __neg__(s): return Jet(-s.v, [-a for a in s.g], {k: -s.H[k] for k in PAIRS})
    def __sub__(s, o): return s + (-o)
    def __rsub__(s, o): return (-s) + o
    def __mul__(s, o):
        if isinstance(o, Jet):
            g = [s.g[a] * o.v + s.v * o.g[a] for a in range(4)]
            H = {(a, b): s.H[(a, b)] * o.v + s.v * o.H[(a, b)] + s.g[a] * o.g[b] + s.g[b] * o.g[a] for (a, b) in PAIRS}
            return Jet(s.v * o.v, g, H)
        return Jet(s.v * o, [a * o for a in s.g], {k: s.H[k] * o for k in PAIRS})
    __rmul__ = __mul__
def jet_x(xs):
    n = len(xs[0].lo); one = I(np.ones(n)); z = I(np.zeros(n))
    return [Jet(xs[a], [one if b == a else z for b in range(4)], {k: z for k in PAIRS}) for a in range(4)]
def jet_radial(phi0, phi1, phi2, xs):
    """jet of phi(s), s = |x|^2: grad 2 phi' x_a, Hess 4 phi'' x_a x_b + 2 phi' delta_ab"""
    g = [phi1 * xs[a] * 2.0 for a in range(4)]
    H = {}
    for (a, b) in PAIRS:
        t = phi2 * (xs[a] * xs[b]) * 4.0
        if a == b: t = t + phi1 * 2.0
        H[(a, b)] = t
    return Jet(phi0, g, H)
def theta_derivs(s, eps):
    """theta(s) = (1 + eps s^-2)^(-1/2) and its first two s-derivatives (I)"""
    si = irecip(s); w = eps * si * si
    b = irecip(isqrt(w + 1.0))          # (1+w)^(-1/2)
    b3 = b * b * b; b5 = b3 * b * b
    t0 = b
    t1 = eps * si * si * si * b3
    t2 = (eps * si * si * si * si * b3) * (-3.0) + (eps * eps * si * si * si * si * si * si * b5) * 3.0
    return t0, t1, t2
def jet_from_derivs(dr, key, n):
    """dr[beta] = I enclosure of d^beta psi over the box (|beta| <= 4); jet of d_a d_b psi"""
    a, b = key
    base = E2(a, b)
    def add(e, c):
        e = list(e); e[c] += 1; return tuple(e)
    v = dr[base]
    g = [dr[add(base, c)] for c in range(4)]
    H = {(c, d): dr[add(add(base, c), d)] for (c, d) in PAIRS}
    return Jet(v, g, H)
def hessian_complex(J):
    """complex Hessian (2x2 Hermitian: m11, m22 real, m12 = re + i im) and complex gradient of a real jet"""
    H = J.H
    m11 = (H[(0, 0)] + H[(1, 1)]) * 0.25; m22 = (H[(2, 2)] + H[(3, 3)]) * 0.25
    mre = (H[(0, 2)] + H[(1, 3)]) * 0.25; mim = (H[(0, 3)] - H[(1, 2)]) * 0.25
    G1 = (J.g[0] * 0.5, -J.g[1] * 0.5); G2 = (J.g[2] * 0.5, -J.g[3] * 0.5)
    return (m11, m22, mre, mim), (G1, G2)
def g0_norm(X, xs, theta):
    """upper bound of the g0-norm (1.6) of a Hermitian form X = (x11, x22, xre, xim) at points x (I coords) with
       radial g0-eigenvalue theta and tangential 1/theta: spectral radius of [[X_rr/theta, X_rt],[., theta X_tt]]."""
    x11, x22, xre, xim = X
    z1 = (xs[0], xs[1]); z2 = (xs[2], xs[3])
    a1 = isq(xs[0]) + isq(xs[1]); a2 = isq(xs[2]) + isq(xs[3]); s = a1 + a2
    # P = z1 conj(z2)
    pre = xs[0] * xs[2] + xs[1] * xs[3]; pim = xs[1] * xs[2] - xs[0] * xs[3]
    rP = xre * pre - xim * pim                      # Re(X12 z1 conj z2)
    Xrr = (x11 * a1 + x22 * a2 + rP * 2.0) * irecip(s)
    Xtt = (x11 * a2 + x22 * a1 - rP * 2.0) * irecip(s)
    # X_rt = (x22 - x11) z1 z2 + X12 z1^2 - conj(X12) z2^2
    q_re = xs[0] * xs[2] - xs[1] * xs[3]; q_im = xs[0] * xs[3] + xs[1] * xs[2]          # z1 z2
    s1_re = isq(xs[0]) - isq(xs[1]); s1_im = xs[0] * xs[1] * 2.0                          # z1^2
    s2_re = isq(xs[2]) - isq(xs[3]); s2_im = xs[2] * xs[3] * 2.0                          # z2^2
    d = x22 - x11
    rt_re = d * q_re + (xre * s1_re - xim * s1_im) - (xre * s2_re + xim * s2_im)
    rt_im = d * q_im + (xre * s1_im + xim * s1_re) - (xre * s2_im - xim * s2_re)
    si = irecip(s)
    rt_re = rt_re * si; rt_im = rt_im * si
    ti = irecip(theta)
    A = Xrr * ti; Bt = Xtt * theta
    tr = (A + Bt) * 0.5; df = (A - Bt) * 0.5
    rad = isqrt(isq(df) + isq(rt_re) + isq(rt_im))
    return _up(iabs(tr).hi + rad.hi)
def frob_K(K):
    k11, k22, kre, kim = K
    return isqrt(isq(k11) + isq(k22) + (isq(kre) + isq(kim)) * 2.0).hi

def ricci_from_jet(Qj, Qv_mag, eps, xs, theta):
    """g0-norm of d dbar Ftilde, Ftilde = log(1 + eps^2 Q)/eps^2:  M/(1+eps^2 Q) - eps^2 G Gbar/(1+eps^2 Q)^2"""
    (m11, m22, mre, mim), (G1, G2) = hessian_complex(Qj)
    e2 = eps * eps
    den = 1.0 - _up(e2.hi * Qv_mag)
    if np.any(den <= 0.5): raise FloatingPointError('eps^2 Q not small')
    nM = g0_norm((m11, m22, mre, mim), xs, theta) / den
    gn = (isq(G1[0]) + isq(G1[1]) + isq(G2[0]) + isq(G2[1])).hi           # |G|^2 (Euclidean); g0^{-1} <= 1/theta
    nG = _up(e2.hi * gn / (den * den) / iabs(theta).lo)
    return _up(nM + nG)


# =====================================================================================================
# 12. bulk evaluator: Q by Taylor models, Ricci and e by interval 2-jets
# =====================================================================================================
BINOM_MH = {1: -0.5, 2: 0.375, 3: -0.3125, 4: 0.2734375, 5: -0.24609375}      # binom(-1/2, k), exact dyadics
def theta_tm(sinv2, eps):
    w = sinv2 * eps
    wr = imax0(w.range())
    w0 = 0.5 * (wr.lo + wr.hi)
    v = w - w0
    vb = v.bound()
    b = irecip(isqrt(I(w0) + 1.0))
    b2 = b * b
    out = TM(I(np.zeros_like(w.coef.lo), np.zeros_like(w.coef.hi)), np.zeros(len(w0)), w.hp) + b
    vk = None; bk = b
    for k in range(1, KT + 1):
        vk = v if vk is None else vk * v
        bk = bk * b2
        out = out + vk * (bk * BINOM_MH[k])
    tail = abs(BINOM_MH[KT + 1]) * upow(vb, KT + 1)
    return TM(out.coef, _up(out.rem + tail * (1 + 1e-12)), out.hp)

def deriv_ranges(a, rem, hp, betas):
    out = {}
    for b in betas:
        r = poly_range_sym(deriv_poly(a, b), hp)
        m = rem[:, sum(b)]
        out[b] = I(_dn(r.lo - m), _up(r.hi + m))
    return out

def bulk_eval(Dd, c, h, eps_hi, want_R=True):
    """c (n,4) centres, h (n,) half-widths (cubes), eps_hi (n,) upper end of eps = a^4 on the box.
       Returns dict: Q range (I), |Q| lower bound at the centre, R upper bound (|Ric|_g/a^6), e upper bound."""
    n = len(c)
    dmax = 2.0 * h * (1 + 1e-12)
    a1, rem1, a2, rem2 = bulk_sources(Dd, c, h, eps_hi)
    hp = hpow(h, KP)
    X = []
    for a in range(4):
        co = np.zeros((n, NK)); co[:, 0] = c[:, a]
        e = [0] * 4; e[a] = 1; co[:, MON_IDX[tuple(e)]] = 1.0
        X.append(TM(I(co), np.zeros(n), hp))
    Aq = quad_A(X)
    D1 = {k: tm_from_poly(deriv_poly(a1, E2(*k)), rem1[:, 2], hp) for k in PAIRS}
    D2 = {k: tm_from_poly(deriv_poly(a2, E2(*k)), rem2[:, 2], hp) for k in PAIRS}
    zero = np.zeros((1, 4))
    c3, r3 = near_taylor(c, zero, 3, KT, dmax)
    c2, r2 = near_taylor(c, zero, 2, KT, dmax)
    sinv3 = TM(c3, r3[:, 0], hp); sinv2 = TM(c2, r2[:, 0], hp)
    eps = I(np.zeros(n), eps_hi)
    theta = theta_tm(sinv2, eps)
    sinv4 = sinv2 * sinv2 if MUT['drop_l0_cancel'] else None
    Q, Mt = Q_formula(D1, D2, Aq, sinv3, theta, eps, sinv4=sinv4)
    Qr = Q.range(); Qc = Q.at_centre()
    check_finite('bulk_eval Q', Qr, Qc)
    out = {'Q': Qr, 'Qmag': Qr.mag(), 'Qc_mig': Qc.mig(), 'Qc': Qc}
    # positivity norm e <= eps |K1 + eps K2~|_F / theta_min
    K = K_of(Mt)
    rc = np.maximum(np.sqrt((c * c).sum(1)) * (1 - 1e-12) - dmax, 1e-300)
    th_min = _dn(1.0 / np.sqrt(1.0 + eps_hi / rc ** 4) * (1 - 1e-12))          # theta >= 1/sqrt(1 + eps/s^2)
    fr = frob_K(tuple(k.range() for k in K))
    out['e'] = _up(eps_hi * fr / th_min)
    if want_R:
        betas = dmaps()['betas']
        dr1 = deriv_ranges(a1, rem1, hp, betas); dr2 = deriv_ranges(a2, rem2, hp, betas)
        xs = [I(_dn(c[:, a] - h), _up(c[:, a] + h)) for a in range(4)]
        Xj = jet_x(xs); Aqj = quad_A(Xj)
        D1j = {k: jet_from_derivs(dr1, k, n) for k in PAIRS}
        D2j = {k: jet_from_derivs(dr2, k, n) for k in PAIRS}
        s = isq(xs[0]) + isq(xs[1]) + isq(xs[2]) + isq(xs[3])
        si = irecip(s)
        s3j = jet_radial(si * si * si, si * si * si * si * (-3.0), si * si * si * si * si * 12.0, xs)
        t0, t1, t2 = theta_derivs(s, eps)
        thj = jet_radial(t0, t1, t2, xs)
        Qj, _ = Q_formula(D1j, D2j, Aqj, s3j, thj, eps)
        e_ok = np.minimum(out['e'], 0.5)
        Rn = ricci_from_jet(Qj, out['Qmag'], eps, xs, t0)
        out['R'] = _up(np.sqrt(eps_hi) * Rn / (1.0 - e_ok))
        out['R'][out['e'] >= 0.5] = np.inf
    check_finite_arr('bulk_eval e', out['e'])
    return out

# =====================================================================================================
# 13. double-precision reference (Ewald, eta = 6; used in controls only)
# =====================================================================================================
class DoubleRef:
    def __init__(self, Z2, eta=6.0, rcut=3.6, kcut=5.2):
        B = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [.5, .5, .5, .5]]) / 2.0
        def pts(Bm, rad):
            nn = int(math.ceil(rad / 0.2)) + 2; r = np.arange(-nn, nn + 1)
            cc = np.array(list(itertools.product(r, repeat=4)), float) @ Bm
            return cc[(cc ** 2).sum(1) <= rad * rad + 1e-12]
        self.LS = pts(B, rcut); self.j0 = int(np.argmin((self.LS ** 2).sum(1)))
        KS = pts(np.linalg.inv(B).T, kcut); KS = KS[(KS ** 2).sum(1) > 1e-9]
        self.KS = KS; K2 = (KS ** 2).sum(1); self.KC = np.exp(-math.pi ** 2 * K2 / eta) / K2 * 32.0
        self.eta = eta; self.Z2 = Z2; self.B = B
        L1 = self.LS[(self.LS ** 2).sum(1) > 1e-12]
        self.L2far = pts(B, 6.0); self.L2far = self.L2far[(self.L2far ** 2).sum(1) > 1e-12]
        self.Z2t = float(((self.L2far ** 2).sum(1) ** -4.0).sum())     # truncated consistently (isotropic tails cancel)
    def pieces(self, z):
        """real Hessians of psi_1 and psi_2 at z (relative to the bubble at 0)"""
        eta = self.eta; out1 = np.empty((len(z), 4, 4)); out2 = np.empty((len(z), 4, 4))
        for i0 in range(0, len(z), 10):
            x = z[i0:i0 + 10]; m = len(x)
            d = x[:, None, :] - self.LS[None]; s = (d ** 2).sum(2)
            s_ = np.where(s < 1e-20, 1.0, s); e = np.exp(-eta * s_)
            vp = -e * (eta / s_ + 1 / s_ ** 2); vpp = e * (eta ** 2 / s_ + 2 * eta / s_ ** 2 + 2 / s_ ** 3)
            s0 = s[:, self.j0]
            vp0 = np.empty(m); vpp0 = np.empty(m)
            for k in range(m):
                sv = s0[k]
                if eta * sv < 0.1:
                    a1 = 0.0; a2 = 0.0
                    for nn in range(2, 25):
                        cf = (-eta) ** nn / math.factorial(nn)
                        a1 += cf * (nn - 1) * sv ** (nn - 2)
                        if nn >= 3: a2 += cf * (nn - 1) * (nn - 2) * sv ** (nn - 3)
                    vp0[k] = a1; vpp0[k] = a2
                else:
                    ee = math.exp(-eta * sv); em1 = math.expm1(-eta * sv)
                    vp0[k] = -eta * ee / sv - em1 / sv ** 2
                    vpp0[k] = eta ** 2 * ee / sv + 2 * eta * ee / sv ** 2 + 2 * em1 / sv ** 3
            vp[:, self.j0] = vp0; vpp[:, self.j0] = vpp0
            H = 4 * np.einsum('nm,nmi,nmj->nij', vpp, d, d)
            tr = 2 * vp.sum(1)
            for i in range(4): H[:, i, i] += tr
            ph = np.cos(2 * math.pi * x @ self.KS.T)
            H += -4 * math.pi ** 2 * np.einsum('nk,k,ki,kj->nij', ph, self.KC, self.KS, self.KS)
            H1 = -0.5 * H
            for i in range(4): H1[:, i, i] += 16 * math.pi ** 2
            out1[i0:i0 + m] = H1
            # psi_2 = (1/24) sum_{l != 0} |z-l|^-6 - (Z2/8)|z|^2: Hessian of |y|^-6 = -6 |y|^-8 (2 delta) + 48 |y|^-10 y y^T
            dd = x[:, None, :] - self.L2far[None]; q = (dd ** 2).sum(2)
            H2 = 48 * np.einsum('nm,nmi,nmj->nij', q ** -5, dd, dd)
            tr2 = -6 * (q ** -4).sum(1)
            for i in range(4): H2[:, i, i] += tr2
            H2 = H2 / 24.0
            for i in range(4): H2[:, i, i] -= self.Z2t / 4.0
            out2[i0:i0 + m] = H2
        return out1, out2
    def Q(self, z, eps=0.0):
        H1, H2 = self.pieces(z)
        x = z
        Jx = np.stack([-x[:, 1], x[:, 0], -x[:, 3], x[:, 2]], 1)
        A = np.einsum('ni,nj->nij', x, x) + np.einsum('ni,nj->nij', Jx, Jx)
        N1 = 0.25 * np.einsum('nij,nij->n', A, H1); N2 = 0.25 * np.einsum('nij,nij->n', A, H2)
        L2 = 0.25 * np.einsum('nii->n', H2)
        M = H1 + eps * H2
        k11 = 0.25 * (M[:, 0, 0] + M[:, 1, 1]); k22 = 0.25 * (M[:, 2, 2] + M[:, 3, 3])
        kre = 0.25 * (M[:, 0, 2] + M[:, 1, 3]); kim = 0.25 * (M[:, 0, 3] - M[:, 1, 2])
        s = (x * x).sum(1); th = s / np.sqrt(s * s + eps)
        return th * ((N1 + eps * N2) / s ** 3 + L2) + (k11 * k22 - kre ** 2 - kim ** 2)


# =====================================================================================================
# 14. NECK, leading Ricci term in closed form.  At |z| = r, eps = u r^4 (u = eps/s^2 in [0,1]) and r -> 0:
#   r^2 d dbar_z Q -> M_lead = d dbar_x [ theta_u(s) A(x) s^-3 ] at the unit vector x = omega,
#   A = N[P_6[psi_1]] = -(5 c6/2) H6, theta_u(s) = (1 + u/s^2)^(-1/2).  In the frame e_r = z, e_t = (-zb2, zb1):
#   M_rr = (phi'' + 7 phi' + 5 phi) A,  M_tt = (phi' - 5 phi) A,  M_rt = phi' conj(t.D) + phi K_rt,
#   phi = theta_u s^-3 at s = 1 (derivatives in s), D_i = d_i A, K_rt = sum K_ij z_i conj(t_j)  (N[H6] = 5 H6,
#   Lap H6 = 0, Euler).  g0 = theta (radial), 1/theta (tangential) => |M|_g0 = |a A|/2 + sqrt(b^2 A^2/4 + |M_rt|^2).
#   L(omega, u) = sqrt(u) |M_lead|_g0  is the a -> 0 limit of |Ric|_g / a^6 at |z| = a u^(-1/4).
# =====================================================================================================
class Du:
    """first-order forward AD over intervals (4 variables: 3 face coordinates + u)."""
    __slots__ = ('v', 'g')
    def __init__(s, v, g): s.v = v; s.g = g
    def __add__(s, o):
        if isinstance(o, Du): return Du(s.v + o.v, [a + b for a, b in zip(s.g, o.g)])
        return Du(s.v + o, s.g)
    __radd__ = __add__
    def __neg__(s): return Du(-s.v, [-a for a in s.g])
    def __sub__(s, o):
        if isinstance(o, Du): return Du(s.v - o.v, [a - b for a, b in zip(s.g, o.g)])
        return Du(s.v - o, s.g)
    def __rsub__(s, o): return (-s) + o
    def __mul__(s, o):
        if isinstance(o, Du): return Du(s.v * o.v, [a * o.v + s.v * b for a, b in zip(s.g, o.g)])
        return Du(s.v * o, [a * o for a in s.g])
    __rmul__ = __mul__
def du_recip(x):
    r = irecip(x.v); r2 = r * r
    return Du(r, [-(a * r2) for a in x.g])
def du_sqrt(x):
    """requires x.v > 0 for the gradient"""
    v = isqrt(x.v)
    ok = x.v.lo > 0
    inv2 = irecip(iwhere(ok, v, I(np.ones_like(v.lo)))) * 0.5
    return Du(v, [a * inv2 for a in x.g]), ok
def gsq(x): return x * x if not isinstance(x, I) else isq(x)

def h6_parts(x):
    """H6 data at a UNIT vector x (generic arithmetic): z1, z2 (complex pairs), H, D1, D2 (complex), K11 (real),
       K12 (complex);  H = H6 = (3/8) Re E, D_i = d_i H, K_ij = d_i dbar_j H (WP2a angular_invariants)."""
    x0, x1, x2, x3 = x
    def cm(a, b): return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])
    def ca(a, b): return (a[0] + b[0], a[1] + b[1])
    def cs(a, b): return (a[0] - b[0], a[1] - b[1])
    def ck(a, k): return (a[0] * k, a[1] * k)
    def cc(a): return (a[0], -a[1])
    z1 = (x0, x1); z2 = (x2, x3)
    z1_2 = cm(z1, z1); z2_2 = cm(z2, z2); z1_4 = cm(z1_2, z1_2); z2_4 = cm(z2_2, z2_2)
    z1_3 = cm(z1_2, z1); z2_3 = cm(z2_2, z2); z1_5 = cm(z1_4, z1); z2_5 = cm(z2_4, z2)
    A1 = cs(z1_5, ck(cm(z1, z2_4), 5.0)); A2 = cs(z2_5, ck(cm(z1_4, z2), 5.0))
    d1A1 = ck(cs(z1_4, z2_4), 5.0)
    d2A1 = ck(cm(z1, z2_3), -20.0); d1A2 = ck(cm(z1_3, z2), -20.0)
    d2A2 = (-d1A1[0], -d1A1[1])
    zb1 = cc(z1); zb2 = cc(z2)
    E = ca(cm(zb1, A1), cm(zb2, A2))
    H = E[0] * 0.375
    D1 = ck(ca(ca(cm(zb1, d1A1), cm(zb2, d1A2)), cc(A1)), 0.1875)
    D2 = ck(ca(ca(cm(zb1, d2A1), cm(zb2, d2A2)), cc(A2)), 0.1875)
    K11 = d1A1[0] * 0.375
    K12 = ck(ca(d1A2, cc(d2A1)), 0.1875)
    return z1, z2, H, D1, D2, K11, K12

def lead_parts(x, u, ops):
    """returns (sqrt(u)), H, a, b, |V|^2 pieces for L = sqrt(u) (5/2) c6 [|a H|/2 + sqrt(b^2 H^2/4 + |V|^2)]
       (H-normalised: A = -(5 c6/2) H).  ops: dict with 'recip', 'sqrt' for the arithmetic type."""
    def cm(a, b): return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])
    def cc(a): return (a[0], -a[1])
    z1, z2, H, D1, D2, K11, K12 = h6_parts(x)
    zb1 = cc(z1); zb2 = cc(z2)
    tD = (cm(zb1, D2)[0] - cm(zb2, D1)[0], cm(zb1, D2)[1] - cm(zb2, D1)[1])         # t.D, t = (-zb2, zb1)
    z1z2 = cm(z1, z2); z1s = cm(z1, z1); z2s = cm(z2, z2)
    Krt = (K11 * z1z2[0] * (-2.0) + cm(K12, z1s)[0] - cm(cc(K12), z2s)[0],
           K11 * z1z2[1] * (-2.0) + cm(K12, z1s)[1] - cm(cc(K12), z2s)[1])
    t0 = ops['recip'](ops['sqrt'](u + 1.0))
    t02 = t0 * t0; t03 = t02 * t0; t05 = t03 * t02
    th1 = u * t03
    th2 = u * t03 * (-3.0) + u * u * t05 * 3.0
    ph = t0; ph1 = th1 - t0 * 3.0; ph2 = th2 - th1 * 6.0 + t0 * 12.0
    P = ph2 + ph1 * 7.0 + ph * 5.0; Qq = ph1 - ph * 5.0
    ti = ops['recip'](t0)
    a = P * ti + t0 * Qq; b = P * ti - t0 * Qq
    ctD = (tD[0], -tD[1])
    V = (ph1 * ctD[0] + ph * Krt[0], ph1 * ctD[1] + ph * Krt[1])
    return H, a, b, V

def L_point(x, u, c6):
    """float evaluation (controls)"""
    ops = {'recip': lambda v: 1.0 / v, 'sqrt': np.sqrt}
    H, a, b, V = lead_parts(x, u, ops)
    return np.sqrt(u) * 2.5 * c6 * (np.abs(a * H) / 2 + np.sqrt(b * b * H * H / 4 + V[0] ** 2 + V[1] ** 2))

def lead_eval(face, ylo, yhi, ulo, uhi, c6I):
    """enclosure of L over boxes (face, y-box, u-box): natural extension intersected with the centred form;
       rigorous point value at the box centre (lower bound); split dimension."""
    n = len(face)
    Ylo, Yhi = full_y(face, ylo, yhi)
    fr = FREE[face // 2]; r = np.arange(n)
    xt = x_tight(Ylo, Yhi)
    Yi = [I(Ylo[:, i], Yhi[:, i]) for i in range(4)]
    nn = isq(Yi[0]) + isq(Yi[1]) + isq(Yi[2]) + isq(Yi[3])
    rr = irecip(isqrt(nn)); r3 = rr * rr * rr
    z = I(np.zeros(n))
    xd = []
    for i in range(4):
        g = []
        for a in range(3):
            ja = fr[:, a]
            Yj = I(Ylo[r, ja], Yhi[r, ja])
            gi = iwhere(ja == i, rr, z) - (Yi[i] * Yj) * r3
            g.append(gi)
        g.append(z)
        xd.append(Du(xt[i], g))
    U = Du(I(ulo, uhi), [z, z, z, I(np.ones(n))])
    def dsqrt(v):
        out, ok = du_sqrt(v); return out
    H, a, b, V = lead_parts(xd, U, {'recip': du_recip, 'sqrt': dsqrt})
    k = c6I * 2.5
    aH = a * H; bH = b * H
    disc = bH * bH * 0.25 + V[0] * V[0] + V[1] * V[1]
    sq_u, ok_u = du_sqrt(U)
    sq_d, ok_d = du_sqrt(disc)
    # natural enclosure
    nat = sq_u.v * (iabs(aH.v) * 0.5 + isqrt(imax0(disc.v))) * k
    # point value at the centre
    yc = 0.5 * (ylo + yhi); uc = 0.5 * (ulo + uhi)
    Yc = full_y(face, yc, yc)[0]
    nr = I(Yc[:, 0]) * I(Yc[:, 0]) + I(Yc[:, 1]) * I(Yc[:, 1]) + I(Yc[:, 2]) * I(Yc[:, 2]) + I(Yc[:, 3]) * I(Yc[:, 3])
    rc = irecip(isqrt(nr))
    xc = [I(Yc[:, i]) * rc for i in range(4)]
    Hc, ac, bc, Vc = lead_parts(xc, I(uc), {'recip': irecip, 'sqrt': isqrt})
    dc = isq(bc * Hc) * 0.25 + isq(Vc[0]) + isq(Vc[1])
    Lc = isqrt(I(uc)) * (iabs(ac * Hc) * 0.5 + isqrt(imax0(dc))) * k
    # centred form where differentiable
    sgn_ok = (aH.v.lo > 0) | (aH.v.hi < 0)
    ok = sgn_ok & ok_d & ok_u
    sg = np.where(aH.v.lo > 0, 1.0, -1.0)
    hs = [0.5 * (yhi[:, a] - ylo[:, a]) for a in range(3)] + [0.5 * (uhi - ulo)]
    Lc_ = Lc
    out = nat
    smear = np.zeros((n, 4))
    grad = []
    for d in range(4):
        gd = (sq_u.g[d] * (iabs(aH.v) * 0.5 + sq_d.v) + sq_u.v * (aH.g[d] * sg * 0.5 + sq_d.g[d])) * k
        grad.append(gd)
        smear[:, d] = gd.mag() * hs[d]
    cf = Lc_
    for d in range(4):
        cf = cf + grad[d] * I(-hs[d], hs[d])
    out = iwhere(ok, isect(nat, cf), nat)
    if np.any(out.lo > out.hi): raise FloatingPointError('lead_eval: empty intersection (bug)')
    widths = np.stack(hs, 1)
    smear = np.where(ok[:, None], smear, widths * np.array([1.0, 1.0, 1.0, 0.5]))
    check_finite('lead_eval', out, Lc)
    return out.hi, Lc.lo, np.argmax(smear, 1)

def lead_new_state(Lm, nu_cells=16):
    face, ylo, yhi = initial_dir_boxes(Lm)
    k = len(face)
    U = np.linspace(0, 1, nu_cells + 1)
    rep = np.repeat(np.arange(nu_cells), k)
    B = {'face': np.tile(face, nu_cells), 'ylo': np.tile(ylo, (nu_cells, 1)), 'yhi': np.tile(yhi, (nu_cells, 1)),
         'ulo': U[rep], 'uhi': U[rep + 1]}
    B['ub'] = np.full(len(B['face']), np.inf)
    return {'boxes': B, 'L': 0.0, 'Lpt': None, 'leaves': {'ub': np.zeros(0), 'ulo': np.zeros(0), 'uhi': np.zeros(0)},
            'it': 0, 'n_eval': 0, 'seconds': 0.0, 'done': False, 'U_current': np.inf}

def lead_run(Lm, c6I, st, tol_rel, budget_s, chunk=20000, log=print, max_active=2_000_000):
    t0 = time.time()
    while True:
        B = st['boxes']; nact = len(B['face'])
        if nact == 0: st['done'] = True; break
        if time.time() - t0 > budget_s: break
        order = np.argsort(-B['ub'])
        take = order[:chunk]; rest = order[chunk:]
        sub = {kk: v[take] for kk, v in B.items()}
        ub, lb, dim = lead_eval(sub['face'], sub['ylo'], sub['yhi'], sub['ulo'], sub['uhi'], c6I)
        st['n_eval'] += len(take)
        j = int(np.argmax(lb))
        if lb[j] > st['L']:
            st['L'] = float(lb[j])
            st['Lpt'] = [int(sub['face'][j]), (0.5 * (sub['ylo'][j] + sub['yhi'][j])).tolist(), float(0.5 * (sub['ulo'][j] + sub['uhi'][j]))]
        L = st['L']; thr = L * (1 + tol_rel)
        fin = ub <= thr
        lv = st['leaves']
        st['leaves'] = {'ub': np.concatenate([lv['ub'], ub[fin]]), 'ulo': np.concatenate([lv['ulo'], sub['ulo'][fin]]),
                        'uhi': np.concatenate([lv['uhi'], sub['uhi'][fin]])}
        act = ~fin
        S = {kk: v[act] for kk, v in sub.items()}; d = dim[act]; ubA = ub[act]
        lo_k = {kk: v.copy() for kk, v in S.items()}; hi_k = {kk: v.copy() for kk, v in S.items()}
        for a in range(3):
            sel = d == a
            mid = 0.5 * (S['ylo'][:, a] + S['yhi'][:, a])
            lo_k['yhi'][sel, a] = mid[sel]; hi_k['ylo'][sel, a] = mid[sel]
        sel = d == 3
        mid = 0.5 * (S['ulo'] + S['uhi'])
        lo_k['uhi'][sel] = mid[sel]; hi_k['ulo'][sel] = mid[sel]
        kids = {kk: np.concatenate([lo_k[kk], hi_k[kk]]) for kk in S}
        kids['ub'] = np.concatenate([ubA, ubA])
        keep = cone_keep_dirs(Lm, kids['face'], kids['ylo'], kids['yhi'])
        kids = {kk: v[keep] for kk, v in kids.items()}
        if len(rest):
            R = {kk: v[rest] for kk, v in B.items()}
            st['boxes'] = {kk: np.concatenate([R[kk], kids[kk]]) for kk in kids}
        else:
            st['boxes'] = kids
        st['it'] += 1
        nb = len(st['boxes']['face'])
        st['U_current'] = float(max(st['boxes']['ub'].max() if nb else 0.0, st['leaves']['ub'].max() if len(st['leaves']['ub']) else 0.0, L))
        if st['it'] % 10 == 0 or nb == 0:
            log('  lead it %d: active %d, L %.10g, U %.10g (gap %.2e), %.0f s' % (st['it'], nb, L, st['U_current'], st['U_current'] / max(L, 1e-300) - 1, time.time() - t0))
        if nb > max_active: log('  lead: active set too large, stopping chunk'); break
    st['seconds'] += time.time() - t0
    return st

def lead_sup_on(st, ulo, uhi):
    """upper bound of sup L over u in [ulo, uhi] (all directions) from the finished leaves (and active boxes)."""
    lv = st['leaves']; B = st['boxes']
    ubs = []
    sel = (lv['uhi'] >= ulo) & (lv['ulo'] <= uhi)
    if np.any(sel): ubs.append(lv['ub'][sel].max())
    if len(B['face']):
        sel = (B['uhi'] >= ulo) & (B['ulo'] <= uhi)
        if np.any(sel): ubs.append(B['ub'][sel].max())
    return max(ubs) if ubs else 0.0

# =====================================================================================================
# 15. NECK corrections: everything except the leading term, by interval 2-jets at unit directions
# =====================================================================================================
def neck_tables(Dd):
    """derivative-coefficient matrices: for every degree d and multi-index beta (|beta| <= 4) the coefficients
       of d^beta P_d as columns over the monomials of degree <= DN1."""
    betas, maps = deriv_maps(DN1, 4)
    nm = nmon(DN1)
    def mats(Pdict):
        out = {}
        for d, P in Pdict.items():
            P = pad(P, DN1)
            lo = np.zeros((nm, len(betas))); hi = np.zeros((nm, len(betas)))
            for jb, b in enumerate(betas):
                src, dst, fac = maps[b]
                v = P[src] * fac
                lo[dst, jb] = v.lo; hi[dst, jb] = v.hi
            out[d] = I(lo, hi)
        return out
    return {'betas': betas, 'P1': mats(Dd.P1), 'P2': mats(Dd.P2)}

def neck_dir_ranges(NT, xs):
    """per degree: I (n, nbeta) ranges of d^beta P_d over the direction boxes (natural, monomial ranges)"""
    MRg = monomial_ranges(xs, DN1)
    out1 = {d: imatmul(MRg, M) for d, M in NT['P1'].items()}
    out2 = {d: imatmul(MRg, M) for d, M in NT['P2'].items()}
    return out1, out2

def rpow_I(rlo, rhi, k):
    """[rlo^k, rhi^k] for 0 <= rlo <= rhi (k >= 0)"""
    if k == 0: return I(1.0)
    lo = 1.0; hi = 1.0
    for _ in range(k): lo = float(_dn(lo * rlo)); hi = float(_up(hi * rhi))
    return I(max(lo, 0.0), hi)

def neck_tail_bound(Dd, nu, m, r_hi, shift):
    """bound for |d^beta (psi_poly_tail)| at |x| = 1, |beta| = m, where the neck polynomial of family nu
       (1: psi_1, 3: psi_2, 5: psi_3 with its prefactor -1/80) misses (i) all degrees > DN (A2 with the full lattice) and
       (ii) the lattice tail |l| > RM of the computed degrees (moment sums are over |l| <= RM); weight r^(d - shift)."""
    DN = DN1 if nu == 1 else DN3
    dmin = 8 if nu == 1 else 4
    tot = 0.0
    for d in range(DN + 2, (93 if nu == 5 else 97), 2):
        tot += _up(Bdm(d, m) * Dd.Cfull[nu][d] * float(upow(np.array(r_hi), d - shift)))
    for d in range(dmin, DN + 1, 2):
        tot += _up(Bdm(d, m) * Dd.Cfar_R[nu][d] * float(upow(np.array(r_hi), d - shift)))
    # degrees beyond 96: ratio per two degrees <= (r_hi/0.5)^2 * 1.5 < 0.05 for r_hi <= 0.08
    return _up(tot * 1.2)

def make_rest_lookup(rest_cells, scale=1.0):
    keys = np.array(list(rest_cells.keys()), float)
    vals = np.array([rest_cells[tuple(k)]['Xr'] for k in rest_cells], float)
    def rest(rlo, rhi, ulo, uhi):
        out = np.zeros(len(rlo))
        covered = np.zeros(len(rlo), bool)
        for (a, b, c, d), v in zip(keys, vals):
            sel = (rhi >= a) & (rlo <= b) & (uhi >= c) & (ulo <= d)
            out = np.where(sel, np.maximum(out, v), out)
            covered |= sel & (rlo >= a) & (rhi <= b) & (ulo >= c) & (uhi <= d)
        # a box must be covered by the union of the cells it meets; cells tile [0,R_NR] x [0,umax(r)]
        return out
    return rest



# =====================================================================================================
# 15b. NECK Ricci to order r^2 in closed form (centred-form branch and bound in (face y, u, r))
#   q~(x) = theta_u(s) [N[psi1^(r)] s^-3 + u N[psi2^(r)] s^-3 + Lap psi2^(r)] + r^8 det(...)
#   psi1^(r) = P6 + r^2 P8 + O(r^4),  psi2^(r) = r^2 c42 |x|^4 + O(r^4)  (L3 (e), (f)).
#   Frame components at the unit vector x (phi = theta_u s^-3, primes = d/ds at s = 1):
#     P6 (A = N[P6] = -(5c6/2) H6, degree 6): rr (phi''+7phi'+5phi) A, tt (phi'-5phi) A, rt phi' conj(t.dA) + phi K_rt(A)
#     P8 (A8 = N[P8], degree 8):              rr (phi''+9phi') A8 + phi N[A8], tt phi' A8 + phi (Lap A8 - N[A8]),
#                                             rt phi' conj(t.dA8) + phi K_rt(A8)
#     psi2: radial g(s) = c42 theta_u(s)(4u/s + 6s): rr g' + g'', tt g', rt 0.
# =====================================================================================================
def _polymul_lin(p, coef_by_var):
    """p (I over monomials of degree <= D) times a linear form sum_a coef_a x_a -> degree D+1"""
    D = int(MON_DEG[p.shape[0] - 1]); n2 = nmon(D + 1)
    lo = np.zeros(n2); hi = np.zeros(n2)
    out = I(lo, hi)
    for a, ca in coef_by_var.items():
        if ca == 0: continue
        idx = []
        for i in range(p.shape[0]):
            e = list(MON[i]); e[a] += 1; idx.append(MON_IDX[tuple(e)])
        v = p * float(ca)
        tmp_lo = np.zeros(n2); tmp_hi = np.zeros(n2)
        tmp_lo[idx] = v.lo; tmp_hi[idx] = v.hi
        out = out + I(tmp_lo, tmp_hi)
    return out
def _polymul_quad(p, terms):
    """p times a quadratic sum c x_a x_b (terms: list of (c, a, b))"""
    D = int(MON_DEG[p.shape[0] - 1]); n2 = nmon(D + 2)
    out = I(np.zeros(n2))
    for (c, a, b) in terms:
        idx = []
        for i in range(p.shape[0]):
            e = list(MON[i]); e[a] += 1; e[b] += 1; idx.append(MON_IDX[tuple(e)])
        v = p * float(c)
        lo = np.zeros(n2); hi = np.zeros(n2); lo[idx] = v.lo; hi[idx] = v.hi
        out = out + I(lo, hi)
    return out
def _pderiv(p, beta):
    D = int(MON_DEG[p.shape[0] - 1])
    betas, maps = deriv_maps(D, sum(beta))
    src, dst, fac = maps[beta]
    m = sum(beta)
    lo = np.zeros(nmon(max(D - m, 0))); hi = np.zeros(nmon(max(D - m, 0)))
    v = p[src] * fac
    lo[dst] = v.lo; hi[dst] = v.hi
    return I(lo, hi)
def _padto(p, D): return pad(p, D)
def poly_N(p):
    D = int(MON_DEG[p.shape[0] - 1])
    Jterms = {0: [(1, 1)], 1: [(-1, 0)], 2: [(1, 3)], 3: [(-1, 2)]}      # (Jx)_a = sign * x_b
    out = I(np.zeros(nmon(D)))
    for (a, b) in PAIRS:
        w = 0.25 if a == b else 0.5
        terms = [(w, a, b)]
        sa, ia = Jterms[a][0]; sb, ib = Jterms[b][0]
        terms.append((w * sa * sb, ia, ib))
        out = out + _padto(_polymul_quad(_pderiv(p, E2(a, b)), terms), D)
    return out
def poly_Lap(p):
    D = int(MON_DEG[p.shape[0] - 1])
    out = I(np.zeros(nmon(D - 2)))
    for a in range(4):
        out = out + _pderiv(p, E2(a, a)) * 0.25
    return out
def frame_polys(A):
    """for a real homogeneous polynomial A of degree d: A, N[A], Lap A, Re/Im (t.dA), Re/Im K_rt(A), all padded to
       degree d (as functions on the unit sphere; Lap A has degree d-2)."""
    D = int(MON_DEG[A.shape[0] - 1])
    d = [_pderiv(A, tuple(int(k == a) for k in range(4))) for a in range(4)]
    tre = (_padto(_polymul_lin(d[0], {2: -1.0}), D) + _padto(_polymul_lin(d[1], {3: 1.0}), D)
           + _padto(_polymul_lin(d[2], {0: 1.0}), D) + _padto(_polymul_lin(d[3], {1: -1.0}), D)) * 0.5
    tim = (_padto(_polymul_lin(d[1], {2: 1.0}), D) + _padto(_polymul_lin(d[0], {3: 1.0}), D)
           + _padto(_polymul_lin(d[3], {0: -1.0}), D) + _padto(_polymul_lin(d[2], {1: -1.0}), D)) * 0.5
    dd = {k: _pderiv(A, E2(*k)) for k in PAIRS}
    K11 = (dd[(0, 0)] + dd[(1, 1)]) * 0.25; K22 = (dd[(2, 2)] + dd[(3, 3)]) * 0.25
    Kre = (dd[(0, 2)] + dd[(1, 3)]) * 0.25; Kim = (dd[(0, 3)] - dd[(1, 2)]) * 0.25
    dK = K22 - K11
    q = lambda P, terms: _padto(_polymul_quad(P, terms), D)
    z1z2_re = [(1, 0, 2), (-1, 1, 3)]; z1z2_im = [(1, 0, 3), (1, 1, 2)]
    z1s_re = [(1, 0, 0), (-1, 1, 1)]; z1s_im = [(2, 0, 1)]
    z2s_re = [(1, 2, 2), (-1, 3, 3)]; z2s_im = [(2, 2, 3)]
    kre = (q(dK, z1z2_re) + q(Kre, z1s_re) - q(Kim, z1s_im) - q(Kre, z2s_re) - q(Kim, z2s_im))
    kim = (q(dK, z1z2_im) + q(Kre, z1s_im) + q(Kim, z1s_re) - q(Kre, z2s_im) + q(Kim, z2s_re))
    NA = poly_N(A); LA = _padto(poly_Lap(A), D)
    return [A, NA, LA, tre, tim, kre, kim]

class NeckPolys:
    """precomputed P8-part polynomials and their first partials, as columns over monomials of degree <= 8"""
    def __init__(self, Dd):
        A8 = poly_N(pad(Dd.P1[8], 8))
        fp = frame_polys(A8)
        self.names = ['A', 'NA', 'LA', 'TDre', 'TDim', 'KRre', 'KRim']
        cols_lo = []; cols_hi = []
        nm = nmon(8)
        for P in fp:
            P = pad(P, 8)
            cols_lo.append(P.lo); cols_hi.append(P.hi)
            for a in range(4):
                g = pad(_pderiv(P, tuple(int(k == a) for k in range(4))), 8)
                cols_lo.append(g.lo); cols_hi.append(g.hi)
        self.M = I(np.stack(cols_lo, 1), np.stack(cols_hi, 1))          # (nm, 35)
        self.c42 = ai(Dd.LS.Z_full[10] / 4)                               # c42 = Z_{5/2}/4 (L3 (f))
        # cross-check: P_4[psi_2] = c42 |x|^4
        p4 = Dd.P2[4]; s2 = np.zeros(nmon(4))
        for i in range(4):
            for j in range(4):
                e = [0] * 4; e[i] += 2; e[j] += 2; s2[MON_IDX[tuple(e)]] += 1.0
        diff = p4 - I(s2) * self.c42
        self.check_c42 = float(np.max(np.abs(diff.lo) + np.abs(diff.hi)))

def frame_terms(Dd):
    """list of (kind, d, rexp, uexp, k, polys): every polynomial part of q~ beyond order r^2, as a homogeneous A of
       degree deg with a radial factor phi = theta_u s^-k:
         psi1 degree d >= 10 : A = N[P1_d],   k = 3, factor r^(d-6)
         psi2 degree d >= 6  : A = N[P2_d],   k = 3, factor u r^(d-2);   A = Lap P2_d, k = 0, factor r^(d-2)."""
    out = []
    for d in range(10, DN1 + 1, 2):
        A = poly_N(pad(Dd.P1[d], d)); out.append(('psi1N', d, d - 6, 0, 3, d, frame_polys(A)))
    for d in range(6, DN3 + 1, 2):
        A = poly_N(pad(Dd.P2[d], d)); out.append(('psi2N', d, d - 2, 1, 3, d, frame_polys(A)))
        A = poly_Lap(pad(Dd.P2[d], d)); out.append(('psi2L', d, d - 2, 0, 0, d - 2, frame_polys(A)))
    # eps psi_3 inside psi~_2 (eps = u r^4): scaled psi~_2 contains u sum_d r^(d+2) P3_d
    for d in range(4, DN3 + 1, 2):
        A = poly_N(pad(Dd.P3[d], d)); out.append(('psi3N', d, d + 2, 2, 3, d, frame_polys(A)))
        A = poly_Lap(pad(Dd.P3[d], d)); out.append(('psi3L', d, d + 2, 1, 0, d - 2, frame_polys(A)))
    return out

def frame_ranges(terms, xs):
    """ranges of the 7 frame polynomials of every term over direction boxes (natural, monomial ranges)"""
    cols_lo = []; cols_hi = []
    for t in terms:
        for P in t[6]:
            P = pad(P, DN1); cols_lo.append(P.lo); cols_hi.append(P.hi)
    M = I(np.stack(cols_lo, 1), np.stack(cols_hi, 1))
    out = []
    n = len(xs[0].lo)
    for s0 in range(0, n, 2048):
        sub = [I(x.lo[s0:s0 + 2048], x.hi[s0:s0 + 2048]) for x in xs]
        out.append(imatmul(monomial_ranges(sub, DN1), M))
    return I(np.concatenate([o.lo for o in out]), np.concatenate([o.hi for o in out]))

def neck_rest(Dd, NT, Lm, cells, nsub=2, nfine=4, log=print):
    """For every cell (rlo, rhi, ulo, uhi), sup over all directions of:
         Cf = |q~ - theta n6|                         (F, used for |z| <= R_NF; interval 2-jets, all terms)
         Xr = |r^2 d dbar Ftilde - M_(order r^2)|_g0 = (exact frame polynomials of the psi1 degrees 10..DN1 and psi2
              degrees 6..DN3, on fine direction boxes) + (interval 2-jets of the tails beyond DN, the eps tail of
              psi~_2, the det term and the nonlinearity of log(1 + eps^2 Q))
         e  = |E|_g0."""
    t0 = time.time()
    # ---- polynomial parts by exact frame polynomials on fine direction boxes
    terms = frame_terms(Dd)
    ff, fyl, fyh = initial_dir_boxes(Lm, ny=8 * nfine)
    FYlo, FYhi = full_y(ff, fyl, fyh)
    fxs = x_tight(FYlo, FYhi)
    FR = frame_ranges(terms, fxs)
    log('  neck rest: %d frame terms on %d fine direction boxes, %.1f s' % (len(terms), len(ff), time.time() - t0))
    # ---- jets for the small remainders
    face, ylo, yhi = initial_dir_boxes(Lm, ny=8 * nsub)
    n = len(face)
    Ylo, Yhi = full_y(face, ylo, yhi)
    xs = x_tight(Ylo, Yhi)
    R1, R2 = neck_dir_ranges(NT, xs)
    betas = NT['betas']; bidx = {b: k for k, b in enumerate(betas)}
    Xj = jet_x(xs); Aqj = quad_A(Xj)
    one = I(np.ones(n))
    s3j = jet_radial(one, one * (-3.0), one * 12.0, xs)
    out = {}
    for (rlo, rhi, ulo, uhi) in cells:
        U = I(ulo, uhi)
        # (a) polynomial parts, exact frame components
        t0u = irecip(isqrt(U + 1.0)); t03 = t0u * t0u * t0u; t05 = t03 * t0u * t0u
        th1 = U * t03; th2 = U * t03 * (-3.0) + U * U * t05 * 3.0
        acc = [I(np.zeros(len(ff))) for _ in range(4)]        # rr, tt, re, im
        col = 0
        for (kind, d, rexp, uexp, k, deg, _) in terms:
            ph = t0u; ph1 = th1 - t0u * float(k); ph2 = th2 - th1 * (2.0 * k) + t0u * float(k * (k + 1))
            fac = rpow_I(rlo, rhi, rexp)
            for _ in range(uexp): fac = fac * U
            A_, NA, LA, TDre, TDim, KRre, KRim = [FR[:, col + j] for j in range(7)]
            col += 7
            rr = ((ph2 + ph1 * float(deg + 1)) * A_ + NA * ph) * fac
            tt = (ph1 * A_ + (LA - NA) * ph) * fac
            re = (ph1 * TDre + ph * KRre) * fac
            im = (ph1 * (-TDim) + ph * KRim) * fac
            acc = [acc[0] + rr, acc[1] + tt, acc[2] + re, acc[3] + im]
        ti = irecip(t0u)
        tr = acc[0] * ti + t0u * acc[1]; df = acc[0] * ti - t0u * acc[1]
        poly_part = float(_up(iabs(tr).hi * 0.5 + isqrt(isq(df) * 0.25 + isq(acc[2]) + isq(acc[3])).hi).max())
        # (b) jets: tails beyond DN, eps tail (psi~_2), det, nonlinear
        def ranges(d1_from, d2_from, kmin):
            dr1c = {}; dr2 = {}
            for b in betas:
                m = sum(b); kb = bidx[b]
                acc1 = I(np.zeros(n))
                for d in range(d1_from, DN1 + 1, 2):
                    acc1 = acc1 + R1[d][:, kb] * rpow_I(rlo, rhi, d - 6)
                tb = neck_tail_bound(Dd, 1, m, rhi, 6)
                dr1c[b] = I(_dn(acc1.lo - tb), _up(acc1.hi + tb))
                acc2 = I(np.zeros(n))
                for d in range(d2_from, DN3 + 1, 2):
                    acc2 = acc2 + R2[d][:, kb] * rpow_I(rlo, rhi, d - 2)
                tb = neck_tail_bound(Dd, 3, m, rhi, 2)
                # eps tail of psi~_2 in scaled units: d_x^beta [D(r x)/r^2] = r^(m-2) (d^beta D)(r x); every term of
                # D = sum_{k>=3} eps^(k-2) psi_k is at least linear in eps = u r^4, hence
                # sup_{r <= rhi} r^(m-2) B_m(eps = u r^4) <= rhi^(m-2) B_m(eps = uhi rhi^4)   (also for m < 2)
                if rhi > 0 and uhi > 0:
                    eh = float(_up(uhi * float(upow(np.array(rhi), 4))))
                    epsB = eps_tail_bounds(Dd, np.zeros((1, 4)), np.array([rhi]), np.array([eh]), kmin=kmin)[0, m]
                    tb = tb + float(_up(epsB * rhi ** (m - 2) * (1 + 1e-12)))
                    if kmin == 4:      # psi_3 handled by frame polynomials up to DN3; its tails here (weight u r^(d+2))
                        tb = tb + float(_up(uhi * neck_tail_bound(Dd, 5, m, rhi, -2)))
                dr2[b] = I(_dn(acc2.lo - tb), _up(acc2.hi + tb))
            return ({k: jet_from_derivs(dr1c, k, n) for k in PAIRS}, {k: jet_from_derivs(dr2, k, n) for k in PAIRS})
        D1a, D2a = ranges(8, 4, 3)                    # all terms beyond P6 (for F and for the full q~)
        D1t, D2t = ranges(DN1 + 2, DN3 + 2, 4)        # tails beyond DN, RM tails, eps^2 tail, psi_3 tails only
        dr6 = {b: R1[6][:, bidx[b]] for b in betas}
        D6 = {k: jet_from_derivs(dr6, k, n) for k in PAIRS}
        t0_, t1_, t2_ = theta_derivs(one, U * one)
        thj = jet_radial(t0_, t1_, t2_, xs)
        r8 = rpow_I(rlo, rhi, 8)
        Qlin_a, _ = Q_formula(D1a, D2a, Aqj, s3j, thj, U * one, det_scale=I(0.0))
        Mt = {k: D6[k] + D1a[k] + D2a[k] * (U * one) for k in PAIRS}
        Kf = K_of(Mt)
        detj = det_of(Kf) * r8
        Qa = Qlin_a + detj
        Qlin_t, _ = Q_formula(D1t, D2t, Aqj, s3j, thj, U * one, det_scale=I(0.0))
        Qt = Qlin_t + detj
        lead = thj * (N_of(D6, Aqj) * s3j)
        Qfull = lead + Qa
        eps = I(0.0, float(_up(uhi * float(rpow_I(rlo, rhi, 4).hi))))          # eps = u r^4
        qmag = Qfull.v.mag()
        e2 = eps.hi ** 2
        den = 1.0 - _up(e2 * qmag)
        if np.any(den < 0.5): raise FloatingPointError('neck_rest: eps^2 q not small')
        Mtl, _ = hessian_complex(Qt)
        Mf, (G1, G2) = hessian_complex(Qfull)
        nt = g0_norm(Mtl, xs, t0_)
        nf = g0_norm(Mf, xs, t0_)
        gn = (isq(G1[0]) + isq(G1[1]) + isq(G2[0]) + isq(G2[1])).hi
        jet_part = _up(nt + nf * _up(e2 * qmag) / den + e2 * gn / (den * den) / iabs(t0_).lo)
        Xr = float(_up(poly_part + float(jet_part.max())))
        Kv = tuple(k.v for k in Kf)
        e = float(_up(uhi * r8.hi * frob_K(Kv) / iabs(t0_).lo).max())
        out[(rlo, rhi, ulo, uhi)] = {'Xr': Xr, 'Xr_poly': poly_part, 'Xr_jets': float(jet_part.max()),
                                     'Cf': float(Qa.v.mag().max()), 'e': e}
    log('  neck rest: %d cells, %.1f s' % (len(cells), time.time() - t0))
    return out

def neck_parts(xd_or_xc, u, r, Pv, c6I, c42, ops):
    """generic (I or Du) assembly of sqrt(u) |M|_g0, M = M_P6 + r^2 (M_P8 + M_psi2).  Pv: dict of P8-part values."""
    def cm(a, b): return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])
    def cc(a): return (a[0], -a[1])
    z1, z2, H, D1, D2, K11, K12 = h6_parts(xd_or_xc)
    zb1 = cc(z1); zb2 = cc(z2)
    tD = (cm(zb1, D2)[0] - cm(zb2, D1)[0], cm(zb1, D2)[1] - cm(zb2, D1)[1])
    z1z2 = cm(z1, z2); z1s = cm(z1, z1); z2s = cm(z2, z2)
    Krt = (K11 * z1z2[0] * (-2.0) + cm(K12, z1s)[0] - cm(cc(K12), z2s)[0],
           K11 * z1z2[1] * (-2.0) + cm(K12, z1s)[1] - cm(cc(K12), z2s)[1])
    t0 = ops['recip'](ops['sqrt'](u + 1.0))
    t02 = t0 * t0; t03 = t02 * t0; t05 = t03 * t02
    th1 = u * t03
    th2 = u * t03 * (-3.0) + u * u * t05 * 3.0
    ph = t0; ph1 = th1 - t0 * 3.0; ph2 = th2 - th1 * 6.0 + t0 * 12.0
    kap = c6I * (-2.5)
    r2 = r * r
    # P6
    A = H * kap
    Mrr = (ph2 + ph1 * 7.0 + ph * 5.0) * A
    Mtt = (ph1 - ph * 5.0) * A
    Mre = (ph1 * tD[0] + ph * Krt[0]) * kap
    Mim = (ph1 * (-tD[1]) + ph * Krt[1]) * kap
    # P8 part
    Mrr = Mrr + r2 * ((ph2 + ph1 * 9.0) * Pv['A'] + ph * Pv['NA'])
    Mtt = Mtt + r2 * (ph1 * Pv['A'] + ph * (Pv['LA'] - Pv['NA']))
    Mre = Mre + r2 * (ph1 * Pv['TDre'] + ph * Pv['KRre'])
    Mim = Mim + r2 * (ph1 * (-Pv['TDim']) + ph * Pv['KRim'])
    # psi2 part (radial)
    g1 = th1 * (u * 4.0 + 6.0) + t0 * (6.0 - u * 4.0)
    g2 = th2 * (u * 4.0 + 6.0) + th1 * (6.0 - u * 4.0) * 2.0 + t0 * u * 8.0
    Mrr = Mrr + r2 * (g1 + g2) * c42
    Mtt = Mtt + r2 * g1 * c42
    ti = ops['recip'](t0)
    tr = Mrr * ti + t0 * Mtt
    df = Mrr * ti - t0 * Mtt
    disc = df * df * 0.25 + Mre * Mre + Mim * Mim
    return tr, disc

def neck2_parts(x, u, a, Pv, c6I, c42, ops, lead_w=1.0):
    """generic (I or Du): tr, disc of the g0-weighted 2x2 matrix  W = sqrt(u) M_lead + a^2 M_2  (sqrt(u) r^2 = a^2),
       so that |W|_g0 = |tr|/2 + sqrt(disc)."""
    def cm(p, q): return (p[0] * q[0] - p[1] * q[1], p[0] * q[1] + p[1] * q[0])
    def cc(p): return (p[0], -p[1])
    z1, z2, H, D1, D2, K11, K12 = h6_parts(x)
    zb1 = cc(z1); zb2 = cc(z2)
    tD = (cm(zb1, D2)[0] - cm(zb2, D1)[0], cm(zb1, D2)[1] - cm(zb2, D1)[1])
    z1z2 = cm(z1, z2); z1s = cm(z1, z1); z2s = cm(z2, z2)
    Krt = (K11 * z1z2[0] * (-2.0) + cm(K12, z1s)[0] - cm(cc(K12), z2s)[0],
           K11 * z1z2[1] * (-2.0) + cm(K12, z1s)[1] - cm(cc(K12), z2s)[1])
    t0 = ops['recip'](ops['sqrt'](u + 1.0))
    t02 = t0 * t0; t03 = t02 * t0; t05 = t03 * t02
    th1 = u * t03
    th2 = u * t03 * (-3.0) + u * u * t05 * 3.0
    ph = t0; ph1 = th1 - t0 * 3.0; ph2 = th2 - th1 * 6.0 + t0 * 12.0
    kap = c6I * (-2.5)
    su = ops['sqrt'](u) * lead_w
    a2 = a * a
    A = H * kap
    rrL = (ph2 + ph1 * 7.0 + ph * 5.0) * A
    ttL = (ph1 - ph * 5.0) * A
    reL = (ph1 * tD[0] + ph * Krt[0]) * kap
    imL = (ph1 * (-tD[1]) + ph * Krt[1]) * kap
    g1 = th1 * (u * 4.0 + 6.0) + t0 * (6.0 - u * 4.0)
    g2 = th2 * (u * 4.0 + 6.0) + th1 * (6.0 - u * 4.0) * 2.0 + t0 * u * 8.0
    rr2 = (ph2 + ph1 * 9.0) * Pv['A'] + ph * Pv['NA'] + (g1 + g2) * c42
    tt2 = ph1 * Pv['A'] + ph * (Pv['LA'] - Pv['NA']) + g1 * c42
    re2 = ph1 * Pv['TDre'] + ph * Pv['KRre']
    im2 = ph1 * (-Pv['TDim']) + ph * Pv['KRim']
    rr = su * rrL + a2 * rr2; tt = su * ttL + a2 * tt2
    re = su * reL + a2 * re2; im = su * imL + a2 * im2
    ti = ops['recip'](t0)
    tr = rr * ti + t0 * tt
    df = rr * ti - t0 * tt
    disc = df * df * 0.25 + re * re + im * im
    return tr, disc

def neck2_eval(NP, boxes, c6I, lead_w=1.0):
    """enclosure of |sqrt(u) M_lead + a^2 M_2|_g0 over boxes (face, y, u, a): natural extension intersected with the
       centred form in (y1, y2, y3, u, a); rigorous value at the centre; split dimension."""
    face, ylo, yhi = boxes['face'], boxes['ylo'], boxes['yhi']
    ulo, uhi, alo, ahi = boxes['ulo'], boxes['uhi'], boxes['alo'], boxes['ahi']
    n = len(face)
    Ylo, Yhi = full_y(face, ylo, yhi)
    fr = FREE[face // 2]; ar = np.arange(n)
    xt = x_tight(Ylo, Yhi)
    Yi = [I(Ylo[:, i], Yhi[:, i]) for i in range(4)]
    nn = isq(Yi[0]) + isq(Yi[1]) + isq(Yi[2]) + isq(Yi[3])
    rr_ = irecip(isqrt(nn)); r3 = rr_ * rr_ * rr_
    z = I(np.zeros(n)); one = I(np.ones(n))
    dx = []
    for i in range(4):
        g = []
        for k in range(3):
            ja = fr[:, k]
            Yj = I(Ylo[ar, ja], Yhi[ar, ja])
            g.append(iwhere(ja == i, rr_, z) - (Yi[i] * Yj) * r3)
        dx.append(g)
    xd = [Du(xt[i], dx[i] + [z, z]) for i in range(4)]
    U = Du(I(ulo, uhi), [z, z, z, one, z])
    Aa = Du(I(alo, ahi), [z, z, z, z, one])
    V = imatmul(monomial_ranges(xt, 8), NP.M)
    # centre point (needed for the centred forms of the P8 values and of the total)
    yc = 0.5 * (ylo + yhi); uc = 0.5 * (ulo + uhi); ac = 0.5 * (alo + ahi)
    Yc = full_y(face, yc, yc)[0]
    nr = I(Yc[:, 0]) * I(Yc[:, 0]) + I(Yc[:, 1]) * I(Yc[:, 1]) + I(Yc[:, 2]) * I(Yc[:, 2]) + I(Yc[:, 3]) * I(Yc[:, 3])
    rcp = irecip(isqrt(nr))
    xc = [I(Yc[:, i]) * rcp for i in range(4)]
    Vc = imatmul(monomial_ranges(xc, 8), NP.M)
    hy = [0.5 * (yhi[:, k] - ylo[:, k]) for k in range(3)]
    Pv = {}
    for k, nm_ in enumerate(NP.names):
        grads = [V[:, 5 * k + 1] * dx[0][j] + V[:, 5 * k + 2] * dx[1][j] + V[:, 5 * k + 3] * dx[2][j] + V[:, 5 * k + 4] * dx[3][j] for j in range(3)]
        val = Vc[:, 5 * k]
        for j in range(3): val = val + grads[j] * I(-hy[j], hy[j])         # mean value form in y
        val = isect(val, V[:, 5 * k])
        Pv[nm_] = Du(val, grads + [z, z])
    def dsqrt(v): return du_sqrt(v)[0]
    upos = ulo > 0
    Us = Du(iwhere(upos, U.v, one), U.g)            # gradient of sqrt(u) only meaningful for u_lo > 0
    tr, disc = neck2_parts(xd, Us, Aa, Pv, c6I, NP.c42, {'recip': du_recip, 'sqrt': dsqrt}, lead_w)
    # natural enclosure with the true u-range (sqrt of [ulo, uhi] including 0)
    trn, dn_ = neck2_parts(xt, U.v, Aa.v, {k: v.v for k, v in Pv.items()}, c6I, NP.c42, {'recip': irecip, 'sqrt': isqrt}, lead_w)
    nat = iabs(trn) * 0.5 + isqrt(imax0(dn_))
    sq_d, ok_d = du_sqrt(disc)
    Pc = {nm_: Vc[:, 5 * kk] for kk, nm_ in enumerate(NP.names)}
    trc, dcc = neck2_parts(xc, I(uc), I(ac), Pc, c6I, NP.c42, {'recip': irecip, 'sqrt': isqrt}, lead_w)
    Lc = iabs(trc) * 0.5 + isqrt(imax0(dcc))
    sgn_ok = (tr.v.lo > 0) | (tr.v.hi < 0)
    ok = sgn_ok & ok_d & upos
    sg = np.where(tr.v.lo > 0, 1.0, -1.0)
    hs = [0.5 * (yhi[:, k] - ylo[:, k]) for k in range(3)] + [0.5 * (uhi - ulo), 0.5 * (ahi - alo)]
    cf = Lc; smear = np.zeros((n, 5))
    for d in range(5):
        gd = tr.g[d] * sg * 0.5 + sq_d.g[d]
        cf = cf + iwhere(ok, gd, z) * I(-hs[d], hs[d])
        smear[:, d] = np.where(ok, gd.mag() * hs[d], 0.0)
    out = iwhere(ok, isect(nat, cf), nat)
    if np.any(out.lo > out.hi): raise FloatingPointError('neck2_eval: empty intersection (bug)')
    widths = np.stack(hs, 1) * np.array([1.0, 1.0, 1.0, 0.5, 0.5 / A0])
    smear = np.where(ok[:, None], smear, widths)
    check_finite('neck2_eval', out, Lc)
    return out.hi, Lc.lo, np.argmax(smear, 1)

def neck2_new_state(Lm, nu_cells=16, pieces=None, ny=16):
    face, ylo, yhi = initial_dir_boxes(Lm, ny=ny)
    k = len(face)
    rows = []
    for (alo, ahi) in (A_PIECES if pieces is None else pieces):
        U = np.linspace(0, 1, nu_cells + 1)
        for j in range(nu_cells): rows.append((U[j], U[j + 1], alo, ahi))
    rows = np.array(rows); m = len(rows)
    B = {'face': np.tile(face, m), 'ylo': np.tile(ylo, (m, 1)), 'yhi': np.tile(yhi, (m, 1)),
         'ulo': np.repeat(rows[:, 0], k), 'uhi': np.repeat(rows[:, 1], k), 'alo': np.repeat(rows[:, 2], k), 'ahi': np.repeat(rows[:, 3], k)}
    B['ub'] = np.full(len(B['face']), np.inf)
    lv = {kk: np.zeros(0) for kk in ('ub', 'ulo', 'uhi', 'alo', 'ahi')}
    return {'boxes': B, 'L': 0.0, 'Lpt': None, 'leaves': lv, 'it': 0, 'n_eval': 0, 'seconds': 0.0, 'done': False, 'U_current': np.inf}

def neck2_run(NP, Lm, c6I, st, tol_rel, budget_s, chunk=20000, log=print, max_active=3_000_000, lead_w=1.0):
    t0 = time.time()
    while True:
        B = st['boxes']; nact = len(B['face'])
        if nact == 0: st['done'] = True; break
        if time.time() - t0 > budget_s: break
        order = np.argsort(-B['ub'])
        take = order[:chunk]; rest_ = order[chunk:]
        sub = {kk: v[take] for kk, v in B.items()}
        ub, lb, dim = neck2_eval(NP, sub, c6I, lead_w)
        st['n_eval'] += len(ub)
        j = int(np.argmax(lb))
        if lb[j] > st['L']:
            st['L'] = float(lb[j])
            st['Lpt'] = [int(sub['face'][j]), (0.5 * (sub['ylo'][j] + sub['yhi'][j])).tolist(), float(0.5 * (sub['ulo'][j] + sub['uhi'][j])), float(0.5 * (sub['alo'][j] + sub['ahi'][j]))]
        L = st['L']; thr = L * (1 + tol_rel)
        fin = ub <= thr
        lv = st['leaves']
        st['leaves'] = {kk: np.concatenate([lv[kk], (ub if kk == 'ub' else sub[kk])[fin]]) for kk in lv}
        act = ~fin
        S = {kk: v[act] for kk, v in sub.items()}; d = dim[act]; ubA = ub[act]
        lo_k = {kk: v.copy() for kk, v in S.items()}; hi_k = {kk: v.copy() for kk, v in S.items()}
        for k in range(3):
            sel = d == k
            mid = 0.5 * (S['ylo'][:, k] + S['yhi'][:, k])
            lo_k['yhi'][sel, k] = mid[sel]; hi_k['ylo'][sel, k] = mid[sel]
        for (dd_, l_, h_) in ((3, 'ulo', 'uhi'), (4, 'alo', 'ahi')):
            sel = d == dd_
            mid = 0.5 * (S[l_] + S[h_])
            lo_k[h_][sel] = mid[sel]; hi_k[l_][sel] = mid[sel]
        kids = {kk: np.concatenate([lo_k[kk], hi_k[kk]]) for kk in S}
        kids['ub'] = np.concatenate([ubA, ubA])
        keep = cone_keep_dirs(Lm, kids['face'], kids['ylo'], kids['yhi'])
        kids = {kk: v[keep] for kk, v in kids.items()}
        if len(rest_):
            R = {kk: v[rest_] for kk, v in B.items()}
            st['boxes'] = {kk: np.concatenate([R[kk], kids[kk]]) for kk in kids}
        else:
            st['boxes'] = kids
        st['it'] += 1
        nb = len(st['boxes']['face'])
        st['U_current'] = float(max(st['boxes']['ub'].max() if nb else 0.0, st['leaves']['ub'].max() if len(st['leaves']['ub']) else 0.0, L))
        if st['it'] % 10 == 0 or nb == 0:
            log('  neck2 it %d: active %d, L %.10g, U %.10g (gap %.2e), %.0f s' % (st['it'], nb, L, st['U_current'], st['U_current'] / max(L, 1e-300) - 1, time.time() - t0))
        if nb > max_active: log('  neck2: active set too large, stopping chunk'); break
    st['seconds'] += time.time() - t0
    return st

def neck2_piece_sup(st, piece):
    """upper bound of sup |sqrt(u) M_lead + a^2 M_2|_g0 over a in the piece (all leaves and pending boxes whose a-range
       meets the piece)."""
    alo, ahi = piece
    out = 0.0
    for src in (st['leaves'], st['boxes']):
        if len(src['ub']) == 0: continue
        sel = (src['ahi'] > alo) & (src['alo'] <= ahi)
        if np.any(sel): out = max(out, float(src['ub'][sel].max()))
    return out

def neck_M2_sup(NP, Lm, nfine=4, nu_cells=16, log=print):
    """S2 >= sup over directions and u in [0,1] of |M_2|_g0, where M = M_lead + r^2 M_2 + (rest) is the order-r^2
       expansion of r^2 d dbar_z Q (section 15b): M_2 = P8 part + psi2 (c42 |x|^4) part.  Since sqrt(u) r^2 = a^2,
       sqrt(u) |M|_g0 <= L(omega, u) + a^2 S2 + sqrt(u) |rest|_g0 (the g0-norm is a norm).  Natural interval
       evaluation of the exact frame polynomials on fine direction boxes."""
    t0 = time.time()
    face, ylo, yhi = initial_dir_boxes(Lm, ny=8 * nfine)
    Ylo, Yhi = full_y(face, ylo, yhi)
    xs = x_tight(Ylo, Yhi)
    cols = [5 * k for k in range(len(NP.names))]
    Mv = I(NP.M.lo[:, cols], NP.M.hi[:, cols])
    V = imatmul(monomial_ranges(xs, 8), Mv)
    P = {nm_: V[:, k] for k, nm_ in enumerate(NP.names)}
    S2 = 0.0; arg = None
    for j in range(nu_cells):
        U = I(j / nu_cells, (j + 1) / nu_cells)
        t0u = irecip(isqrt(U + 1.0)); t03 = t0u * t0u * t0u; t05 = t03 * t0u * t0u
        th1 = U * t03; th2 = U * t03 * (-3.0) + U * U * t05 * 3.0
        ph = t0u; ph1 = th1 - t0u * 3.0; ph2 = th2 - th1 * 6.0 + t0u * 12.0
        g1 = th1 * (U * 4.0 + 6.0) + t0u * (6.0 - U * 4.0)
        g2 = th2 * (U * 4.0 + 6.0) + th1 * (6.0 - U * 4.0) * 2.0 + t0u * U * 8.0
        rr = (ph2 + ph1 * 9.0) * P['A'] + ph * P['NA'] + (g1 + g2) * NP.c42
        tt = ph1 * P['A'] + ph * (P['LA'] - P['NA']) + g1 * NP.c42
        re = ph1 * P['TDre'] + ph * P['KRre']
        im = ph1 * (-P['TDim']) + ph * P['KRim']
        ti = irecip(t0u)
        tr = rr * ti + t0u * tt; df = rr * ti - t0u * tt
        rho = _up(iabs(tr).hi * 0.5 + isqrt(isq(df) * 0.25 + isq(re) + isq(im)).hi)
        k = int(np.argmax(rho))
        if rho[k] > S2: S2 = float(rho[k]); arg = (j, int(face[k]), (0.5 * (ylo[k] + yhi[k])).tolist())
    log('  neck: S2 = %.6g on %d fine direction boxes x %d u-cells (%.1f s)' % (S2, len(face), nu_cells, time.time() - t0))
    return S2, arg


# =====================================================================================================
# 16. bulk Ricci by 2-jets of Taylor models (centred; used where the natural 2-jets are too wide)
# =====================================================================================================
def binom_gen(p, k):
    """binom(-p, k) for p = j/2 (exact Fraction)"""
    from fractions import Fraction as Fr
    v = Fr(1)
    for t in range(k): v *= (Fr(-p) - t)
    return v / FACT[k]
def pow_tm(w, p):
    """TM of (1 + w)^(-p), w >= 0 (p = 1/2, 3/2, 5/2): Taylor series around w0 + Lagrange remainder."""
    from fractions import Fraction as Fr
    wr = imax0(w.range())
    w0 = 0.5 * (wr.lo + wr.hi)
    v = w - w0
    vb = v.bound()
    b = irecip(isqrt(I(w0) + 1.0))                  # (1+w0)^(-1/2)
    base = I(np.ones_like(w0))
    for _ in range(int(round(2 * p))): base = base * b      # (1+w0)^(-p)
    inv = b * b                                      # (1+w0)^(-1)
    out = TM(I(np.zeros_like(w.coef.lo), np.zeros_like(w.coef.hi)), np.zeros(len(w0)), w.hp) + base
    vk = None; bk = base
    for k in range(1, KT + 1):
        vk = v if vk is None else vk * v
        bk = bk * inv
        out = out + vk * (bk * I.c(binom_gen(Fr(p), k)))
    tail = abs(float(binom_gen(Fr(p), KT + 1))) * (1 + 1e-12) * upow(vb, KT + 1)
    return TM(out.coef, _up(out.rem + tail), out.hp)

def bulk_ricci_tm(Dd, c, h, eps_hi, a1, rem1, a2, rem2, Qmag, e):
    n = len(c); hp = hpow(h, KP); dmax = 2.0 * h * (1 + 1e-12)
    X = []
    for a in range(4):
        co = np.zeros((n, NK)); co[:, 0] = c[:, a]
        ee = [0] * 4; ee[a] = 1; co[:, MON_IDX[tuple(ee)]] = 1.0
        X.append(TM(I(co), np.zeros(n), hp))
    def tmd(a, rem, b):
        return tm_from_poly(deriv_poly(a, b), rem[:, sum(b)], hp)
    def add1(e, k):
        e = list(e); e[k] += 1; return tuple(e)
    def jetD(a, rem, key):
        base = E2(*key)
        v = tmd(a, rem, base)
        g = [tmd(a, rem, add1(base, k)) for k in range(4)]
        H = {(k, l): tmd(a, rem, add1(add1(base, k), l)) for (k, l) in PAIRS}
        return Jet(v, g, H)
    zero = TM(I(np.zeros((n, NK))), np.zeros(n), hp)
    one = zero + 1.0
    Xj = [Jet(X[a], [one if b == a else zero for b in range(4)], {k: zero for k in PAIRS}) for a in range(4)]
    Aqj = quad_A(Xj)
    D1j = {k: jetD(a1, rem1, k) for k in PAIRS}
    D2j = {k: jetD(a2, rem2, k) for k in PAIRS}
    sv = {}
    for nu in (2, 3, 4, 5, 6):
        cf, rr = near_taylor(c, np.zeros((1, 4)), nu, KT, dmax)
        sv[nu] = TM(cf, rr[:, 0], hp)
    eps = I(np.zeros(n), eps_hi)
    w = sv[2] * eps
    b1 = pow_tm(w, 0.5); b3 = pow_tm(w, 1.5); b5 = pow_tm(w, 2.5)
    th0 = b1
    th1 = sv[3] * b3 * eps
    th2 = (sv[4] * b3 * eps) * (-3.0) + (sv[6] * b5 * (eps * eps)) * 3.0
    thj = jet_radial(th0, th1, th2, X)
    s3j = jet_radial(sv[3], sv[4] * (-3.0), sv[5] * 12.0, X)
    Qj, _ = Q_formula(D1j, D2j, Aqj, s3j, thj, eps)
    (m11, m22, mre, mim), (G1, G2) = hessian_complex(Qj)
    Mr = tuple(t.range() for t in (m11, m22, mre, mim))
    xs = [I(_dn(c[:, a] - h), _up(c[:, a] + h)) for a in range(4)]
    thr = th0.range()
    e2 = eps * eps
    den = 1.0 - _up(e2.hi * Qmag)
    nM = g0_norm(Mr, xs, thr) / den
    Gr = [G1[0].range(), G1[1].range(), G2[0].range(), G2[1].range()]
    gn = (isq(Gr[0]) + isq(Gr[1]) + isq(Gr[2]) + isq(Gr[3])).hi
    nG = _up(e2.hi * gn / (den * den) / iabs(thr).lo)
    R = _up(np.sqrt(eps_hi) * _up(nM + nG) / (1.0 - np.minimum(e, 0.5)))
    R[e >= 0.5] = np.inf
    return R

# =====================================================================================================
# 17. bulk branch and bound (F everywhere in the bulk, Ricci where |z| >= R_NR)
# =====================================================================================================
def bulk_keep(c, h, Lm, rmin):
    """box may meet {x in V0, x in the G_J cone, |x| >= rmin}"""
    hv = h[:, None]
    v = c @ MIN_VEC.T - hv * np.abs(MIN_VEC).sum(1)[None, :]
    k1 = np.all(_dn(v) <= 0.125, axis=1)
    w = c @ Lm.T + hv * np.abs(Lm).sum(1)[None, :]
    k2 = np.all(_up(w) >= 0, axis=1)
    k3 = _up(np.sqrt((c * c).sum(1)) + 2.0 * h) >= rmin
    return k1 & k2 & k3
def bulk_valid(c, h):
    """size rule for the Taylor models: 2h <= |c|/4 and 2h <= |c - l|/4 for the near lattice points"""
    r = np.sqrt((c * c).sum(1))
    dl = np.sqrt(((c[:, None, :] - NEAR[None]) ** 2).sum(2)).min(1)
    return (2.0 * h <= 0.25 * r) & (2.0 * h <= 0.25 * dl)
def split16(c, h):
    off = np.array(list(itertools.product((-0.5, 0.5), repeat=4)))
    cc = (c[:, None, :] + off[None, :, :] * h[:, None, None]).reshape(-1, 4)
    hh = np.repeat(h * 0.5, 16)
    return cc, hh

def bulk_new_state(Lm):
    c = np.zeros((1, 4)); h = np.array([0.375])          # cube [-0.375, 0.375]^4 contains V0 (radius 0.3536)
    return {'c': c, 'h': h, 'ub': np.array([np.inf]), 'LF': 0.0, 'LFpt': None,
            'leafF': 0.0, 'leafR': 0.0, 'leafE': 0.0, 'leafF_pt': None, 'leafR_pt': None,
            'n_eval': 0, 'n_tmjet': 0, 'n_leaves': 0, 'it': 0, 'seconds': 0.0, 'done': False}

def bulk_run(Dd, Lm, st, tolF, Rcap, budget_s, chunk=256, log=print):
    t0 = time.time()
    if st['LF'] == 0.0:
        # seed the lower bound with a certified point value next to the facet centre (1/4, 0, 0, 0), where the sampled
        # extreme of |F8bulk| sits (only affects efficiency; the reported bound is the max over final boxes)
        cs = np.array([[0.2490234375, 0.0, 0.0, 0.0]])
        o = bulk_eval(Dd, cs, np.array([2.0 ** -14]), np.array([A0 ** 4]), want_R=False)
        st['LF'] = float(o['Qc_mig'][0]); st['LFpt'] = cs[0].tolist()
    while True:
        nact = len(st['h'])
        if nact == 0: st['done'] = True; break
        if time.time() - t0 > budget_s: break
        order = np.argsort(-st['ub'])
        take = order[:chunk]; rest = order[chunk:]
        c = st['c'][take]; h = st['h'][take]
        keep = bulk_keep(c, h, Lm, R_NF)
        c = c[keep]; h = h[keep]
        valid = bulk_valid(c, h)
        newc = []; newh = []; newub = []
        if np.any(~valid):
            cc, hh = split16(c[~valid], h[~valid])
            newc.append(cc); newh.append(hh); newub.append(np.full(len(hh), np.inf))
        c = c[valid]; h = h[valid]
        if len(h):
            r = np.sqrt((c * c).sum(1))
            eps_hi = _up(upow(np.minimum(A0, r + 2.0 * h), 4))
            needR = (r + 2.0 * h) >= R_NR
            o = bulk_eval(Dd, c, h, eps_hi, want_R=True)
            st['n_eval'] += len(h)
            Qm = o['Qmag']
            Fub = _up(Qm / (1.0 - eps_hi * eps_hi * Qm))
            Flb = o['Qc_mig']
            j = int(np.argmax(Flb))
            if Flb[j] > st['LF']: st['LF'] = float(Flb[j]); st['LFpt'] = c[j].tolist()
            R = np.where(needR, o['R'], 0.0)
            # where the natural 2-jet bound exceeds the cap, use the Taylor-model 2-jets
            sel = needR & (R > Rcap)
            if np.any(sel):
                a1, rem1, a2, rem2 = bulk_sources(Dd, c[sel], h[sel], eps_hi[sel])
                Rt = bulk_ricci_tm(Dd, c[sel], h[sel], eps_hi[sel], a1, rem1, a2, rem2, Qm[sel], o['e'][sel])
                R[sel] = np.minimum(R[sel], Rt); st['n_tmjet'] += int(sel.sum())
            LF = st['LF']
            doneF = Fub <= LF * (1 + tolF)
            doneR = R <= Rcap
            fin = doneF & doneR
            if np.any(fin):
                k = np.argmax(np.where(fin, Fub, -1));
                if Fub[k] > st['leafF']: st['leafF'] = float(Fub[k]); st['leafF_pt'] = c[k].tolist()
                k = np.argmax(np.where(fin, R, -1))
                if R[k] > st['leafR']: st['leafR'] = float(R[k]); st['leafR_pt'] = c[k].tolist()
                st['leafE'] = max(st['leafE'], float(o['e'][fin].max()))
                st['n_leaves'] += int(fin.sum())
            sp = ~fin
            if np.any(sp):
                if np.any(h[sp] < 1e-7): raise FloatingPointError('bulk: box size underflow (no convergence)')
                cc, hh = split16(c[sp], h[sp])
                newc.append(cc); newh.append(hh); newub.append(np.repeat(np.where(doneF[sp], 0.0, Fub[sp]), 16))
        # pruning of waiting boxes is not done by F alone (Ricci may still need them); keep all
        if len(rest) or newc:
            st['c'] = np.concatenate([st['c'][rest]] + newc) if newc else st['c'][rest]
            st['h'] = np.concatenate([st['h'][rest]] + newh) if newh else st['h'][rest]
            st['ub'] = np.concatenate([st['ub'][rest]] + newub) if newub else st['ub'][rest]
        else:
            st['c'] = np.zeros((0, 4)); st['h'] = np.zeros(0); st['ub'] = np.zeros(0)
        st['it'] += 1
        if st['it'] % 20 == 0 or len(st['h']) == 0:
            log('  bulk it %d: active %d, evaluated %d (tm-jets %d), LF %.8g, leaf F %.8g, leaf R %.6g, %.0f s' % (
                st['it'], len(st['h']), st['n_eval'], st['n_tmjet'], st['LF'], st['leafF'], st['leafR'], time.time() - t0))
    st['seconds'] += time.time() - t0
    return st


# =====================================================================================================
# 18. Ewald point evaluator in arb (independent route for psi_1; Lemma A6 tails, order 2)
# =====================================================================================================
def ewald_Sreg_hessian(x, eta, T, K):
    """arb 4x4 real Hessian of S_reg = S - 1/|x|^2 at the point x (exact floats), Ewald split eta.
       Real-space terms |x - l| < T (l = 0 replaced by (e^{-eta s} - 1)/s), reciprocal 0 < |k| < K, tails A6."""
    from flint import arb_series
    eta_a = arb(eta); pi = arb.pi(); V = arb(1) / 32
    xa = [arb(float(v)) for v in x]
    H = [[arb(0) for _ in range(4)] for _ in range(4)]
    rx = math.sqrt(sum(v * v for v in x))
    Nmax = int(16 * (T + rx + 0.01) ** 2) + 1
    m, N = L_points_m(Nmax, 1)
    lp = m.astype(float) / 4.0
    dist = np.sqrt(((np.array(x)[None, :] - lp) ** 2).sum(1))
    sel = dist < T + 1e-9
    pts = [None] + [tuple(int(v) for v in row) for row in m[sel]]
    for p in pts:
        if p is None:
            d = xa
        else:
            d = [xa[i] - arb(fmpq(p[i], 4)) for i in range(4)]
        s = d[0] * d[0] + d[1] * d[1] + d[2] * d[2] + d[3] * d[3]
        if p is None:
            Ss = arb_series([s, 1], prec=3)
            Vs = ((-eta_a * Ss).exp() - 1) / Ss
            co = Vs.coeffs() + [arb(0)] * 3
            f1 = co[1]; f2 = 2 * co[2]
        else:
            e = (-eta_a * s).exp()
            f1 = -e * (eta_a / s + 1 / (s * s))
            f2 = e * (eta_a * eta_a / s + 2 * eta_a / (s * s) + 2 / (s * s * s))
        for i in range(4):
            for j in range(i, 4):
                t = 4 * f2 * d[i] * d[j]
                if i == j: t += 2 * f1
                H[i][j] += t
    # reciprocal: L* = {k in (2Z)^4 : (k0+k1+k2+k3)/2 even}
    Mk = int(K / 2) + 1
    rng = np.arange(-Mk, Mk + 1)
    nvec = np.array(np.meshgrid(rng, rng, rng, rng, indexing='ij')).reshape(4, -1).T
    nvec = nvec[(nvec.sum(1) % 2 == 0)]
    k2 = 4 * (nvec * nvec).sum(1)
    nvec = nvec[(k2 > 0) & (k2 < K * K)]
    for nv in nvec:
        kv = [2 * int(v) for v in nv]
        kk = sum(v * v for v in kv)
        ph = 2 * pi * sum(kv[i] * xa[i] for i in range(4))
        w = -4 * pi * pi / V * (-pi * pi * kk / eta_a).exp() / kk * ph.cos()
        for i in range(4):
            for j in range(i, 4):
                H[i][j] += w * kv[i] * kv[j]
    # tails (Lemma A6, n = 2): G_2(r) = e^{-eta r^2} (10 r^-4 + 10 eta r^-2 + 4 eta^2)
    Rr = arb(T) - arb(0.5)
    assert float(Rr.lower()) ** 2 > 2.0 / eta
    cq = {-4: arb(10), -2: 10 * eta_a, 0: 4 * eta_a * eta_a}
    real_tail = 1024 * sum(c * sum(math.comb(3, jj) * arb(4) ** (jj - 3) * Jb(q + jj, Rr, eta_a) for jj in range(4)) for q, c in cq.items())
    Rk = arb(K) - 2 * arb(2).sqrt(); bet = pi * pi / eta_a
    assert float((2 * bet * Rk * Rk).lower()) > 2
    rec_tail = (2 * pi) ** 2 / V * sum(math.comb(3, jj) * arb(2).sqrt() ** (3 - jj) * Jb(jj, Rk, bet) for jj in range(4))
    tail = real_tail + rec_tail
    out = [[None] * 4 for _ in range(4)]
    for i in range(4):
        for j in range(i, 4):
            out[i][j] = out[j][i] = H[i][j] + arb(0, 1) * tail
    return out, tail

def point_F8bulk(Dd, x, eta, T, K):
    """certified F8bulk(x) = Q(x; 0): psi_1 Hessian by Ewald (arb), psi_2 Hessian by the polynomial route."""
    Hs, tail = ewald_Sreg_hessian(x, eta, T, K)
    pi2 = arb.pi() ** 2
    H1 = {}
    for (a, b) in PAIRS:
        v = -Hs[a][b] / 2
        if a == b: v += 16 * pi2
        H1[(a, b)] = ai(v)
    c = np.array([x], float)
    a1, rem1, a2, rem2 = bulk_sources(Dd, c, np.array([1e-9]), np.array([0.0]))
    H2 = {}
    for (a, b) in PAIRS:
        g = E2(a, b); k = MON_IDX[g]
        f = 2.0 if a == b else 1.0
        v = a2[:, k] * f
        H2[(a, b)] = I(_dn(v.lo - rem2[:, 2]), _up(v.hi + rem2[:, 2]))[0]
    xs = [I(float(v)) for v in x]
    Aq = quad_A(xs)
    N1 = N_of(H1, Aq); L2 = L_of(H2)
    det = det_of(K_of(H1))
    s = isq(xs[0]) + isq(xs[1]) + isq(xs[2]) + isq(xs[3])
    si = irecip(s)
    return N1 * si * si * si + L2 + det, float(aup(tail))

def lower_bound_point(ref):
    """candidate for the extreme of |F8bulk|: double-precision scan of the coordinate axis up to the Voronoi facet
       x0 = 1/4 (the sampled extreme sits there; the axis is symmetric about 1/4); returned as a dyadic point."""
    ts = np.linspace(0.2, 0.25, 501)
    z = np.zeros((len(ts), 4)); z[:, 0] = ts
    v = np.abs(ref.Q(z))
    t = float(ts[int(np.argmax(v))])
    t = float(np.round(t * 2 ** 30) / 2 ** 30)
    return [t, 0.0, 0.0, 0.0]

# =====================================================================================================
# 19. assembly
# =====================================================================================================
R_CELLS = [(0.0, 0.01), (0.01, 0.02), (0.02, 0.03), (0.03, 0.04), (0.04, 0.05), (0.05, 0.06), (0.06, 0.07), (0.07, R_NR)]
NU_CELLS = 8
def neck_cells():
    """(rlo, rhi, ulo, uhi) cells tiling {0 <= r <= R_NR, 0 <= u <= min(1, (A0/r_lo)^4)} (same tiling as the B&B)"""
    cells = []
    for (a, b) in R_CELLS:
        umax = 1.0 if a <= A0 else min(1.0, (A0 / a) ** 4 * (1 + 1e-9))
        U = np.linspace(0, umax, NU_CELLS + 1)
        for j in range(NU_CELLS): cells.append((a, b, float(U[j]), float(U[j + 1])))
    return cells
def admissible(cell, piece):
    """is there (r, u) in the cell with a = u^(1/4) r in the piece (a_lo, a_hi]?"""
    rlo, rhi, ulo, uhi = cell; alo, ahi = piece
    amin = (ulo ** 0.25) * rlo * (1 - 1e-12); amax = (uhi ** 0.25) * rhi * (1 + 1e-12)
    return amax > alo and amin <= ahi

def assemble(Dd, lead_st, s2_st, rest, bulk_st, lower):
    """per a-piece: C_out = max(neck, bulk), R_out = max(neck, bulk); neck Ricci (sqrt(u) r^2 = a^2, |.|_g0 a norm):
       R <= [ L + a_hi^2 S2 + sqrt(u_hi) Xr(cell) ] / (1 - e)   (sections 14, 15, 15b)."""
    c6hi = float(Dd.c6.hi)
    Lhi = max(lead_st['U_current'], lead_st['L']); S2 = max(s2_st['U_current'], s2_st['L'])
    res = {'pieces': [], 'lead_L_interval': [lead_st['L'], Lhi], 'lead_argmax(face,y,u)': lead_st['Lpt'],
           'S2_interval': [s2_st['L'], S2], 'S2_argmax(face,y,u,a)': s2_st['Lpt']}
    bulkF = bulk_st['leafF']; bulkR = bulk_st['leafR']; bulkE = bulk_st['leafE']
    for piece in A_PIECES:
        eps_p = piece[1] ** 4
        Cf = 0.0; e_neck = 0.0; Xmax = 0.0
        for cell, v in rest.items():
            if not admissible(cell, piece): continue
            e_neck = max(e_neck, v['e'])
            Xmax = max(Xmax, math.sqrt(cell[3]) * v['Xr'])
            if cell[1] <= R_NF + 1e-12: Cf = max(Cf, v['Cf'])
        QnF = _up(15.0 / 16.0 * c6hi + Cf)                 # |theta n6| <= (5 c6/2)(3/8) = (15/16) c6 (L3 (e))
        Cn = _up(QnF / (1 - eps_p * eps_p * QnF))           # |log(1+x)| <= |x|/(1-|x|)
        a2S2 = _up(piece[1] ** 2 * S2 * (1 + 1e-12))
        Rn = _up((Lhi + a2S2 + Xmax) / (1 - min(e_neck, 0.5)))
        res['pieces'].append({'a_interval': list(piece), 'C_out_upper': max(Cn, bulkF), 'C_neck': Cn, 'C_neck_corr_max': Cf,
                              'C_bulk': bulkF, 'R_out_upper': max(Rn, bulkR), 'R_neck': Rn,
                              'R_neck_terms': {'L_lead': Lhi, 'a2_S2': a2S2, 'sqrt_u_rest_max': Xmax, 'e_neck': e_neck},
                              'R_bulk': bulkR, 'e_upper': max(e_neck, bulkE)})
    res['C_out_upper'] = max(p['C_out_upper'] for p in res['pieces'])
    res['R_out_upper'] = max(p['R_out_upper'] for p in res['pieces'])
    res['e_out_upper'] = max(p['e_upper'] for p in res['pieces'])
    res['C_out_lower'] = bulk_st['LF']; res['C_out_lower_point'] = bulk_st['LFpt']
    res['R_out_lower_a_to_0'] = lead_st['L']           # a -> 0 limit of R along r = a u^(-1/4): lower bound (piece touching 0)
    res['F8bulk_lower'] = lower
    return res

# =====================================================================================================
# 20. state handling, main (certified run) and controls
# =====================================================================================================
def save_npz(path, st, arrays):
    meta = {k: v for k, v in st.items() if k not in arrays}
    payload = {}
    for k in arrays:
        v = st[k]
        if isinstance(v, dict):
            for kk, vv in v.items(): payload[k + '__' + kk] = vv
        else:
            payload[k] = v
    np.savez(path, meta=json.dumps(meta, default=float), **payload)
def load_npz(path, arrays):
    z = np.load(path, allow_pickle=False)
    st = json.loads(str(z['meta']))
    for k in arrays:
        if k in z.files: st[k] = z[k]
        else: st[k] = {kk.split('__', 1)[1]: z[kk] for kk in z.files if kk.startswith(k + '__')}
    return st

LEAD_TOL = 1e-4; S2_TOL = 1e-2; BULK_TOLF = 1e-3; RCAP_FRAC = 0.8

def main(out_path, budget_s):
    if not os.path.exists(GO_FILE):
        print('REFUSED: %s does not exist (the lead creates it after reading the controls).' % GO_FILE)
        sys.exit(2)
    os.makedirs(STATE_DIR, exist_ok=True)
    t0 = time.time()
    def log(s): print(s, flush=True)
    def left(): return budget_s - (time.time() - t0)
    rpath = os.path.join(STATE_DIR, 'main_results.json')
    R = json.load(open(rpath)) if os.path.exists(rpath) else {}
    def save(): json.dump(R, open(rpath, 'w'), indent=1, default=float)
    Dd = Data(log=log)
    R['c6'] = Dd.c6_arb.str(25); R['data_check_P6'] = Dd.check_P6
    W, GJ, hol = group_GJ(); Lm = cone_normals(GJ)
    if len(GJ) != 192: raise RuntimeError('G_J order')
    NP = NeckPolys(Dd); R['check_c42'] = NP.check_c42
    # stage 1: O(r^4) rest of the neck, per cell
    if 'rest' not in R:
        NT = neck_tables(Dd)
        rc = neck_rest(Dd, NT, Lm, neck_cells(), nsub=2, log=log)
        R['rest'] = [[list(k), v] for k, v in rc.items()]; save()
        log('main: neck rest done (%.0f s)' % (time.time() - t0))
    rest_cells = {tuple(k): v for k, v in R['rest']}
    rest = make_rest_lookup(rest_cells)
    # stage 2: neck Ricci.  R <= [L(omega,u) + a^2 S2 + sqrt(u) Xr]/(1-e), L = sqrt(u)|M_lead|_g0, S2 >= sup |M_2|_g0
    LKEYS = ['boxes', 'leaves']
    lp = os.path.join(STATE_DIR, 'lead.npz')
    if not R.get('lead_done'):
        st = load_npz(lp, LKEYS) if os.path.exists(lp) else lead_new_state(Lm)
        while not st['done'] and left() > 40:
            st = lead_run(Lm, Dd.c6, st, LEAD_TOL, min(60.0, left() - 30), log=log)
            save_npz(lp, st, LKEYS)
        R['lead_done'] = bool(st['done']); save()
        if not st['done']:
            log('main: INCOMPLETE (lead), rerun the same command to resume'); sys.exit(3)
    lead_st = load_npz(lp, LKEYS)
    sp = os.path.join(STATE_DIR, 's2.npz')
    if not R.get('s2_done'):
        st = load_npz(sp, LKEYS) if os.path.exists(sp) else neck2_new_state(Lm, 16, pieces=[(1.0, 1.0)], ny=8)
        while not st['done'] and left() > 40:
            st = neck2_run(NP, Lm, Dd.c6, st, S2_TOL, min(60.0, left() - 30), log=log, lead_w=0.0)
            save_npz(sp, st, LKEYS)
        R['s2_done'] = bool(st['done']); save()
        if not st['done']:
            log('main: INCOMPLETE (S2), rerun the same command to resume'); sys.exit(3)
    s2_st = load_npz(sp, LKEYS)
    nk_st = lead_st
    # stage 3: bulk
    bp = os.path.join(STATE_DIR, 'bulk.npz')
    if not R.get('bulk_done'):
        st = load_npz(bp, ['c', 'h', 'ub']) if os.path.exists(bp) else bulk_new_state(Lm)
        Rcap = RCAP_FRAC * nk_st['L']
        while not st['done'] and left() > 40:
            st = bulk_run(Dd, Lm, st, BULK_TOLF, Rcap, min(60.0, left() - 30), log=log)
            save_npz(bp, st, ['c', 'h', 'ub'])
        R['bulk_done'] = bool(st['done']); save()
        if not st['done']:
            log('main: INCOMPLETE (bulk), rerun the same command to resume'); sys.exit(3)
    bulk_st = load_npz(bp, ['c', 'h', 'ub'])
    # stage 4: certified point lower bound for sup |F8bulk|
    if 'lower' not in R:
        if left() < 120: log('main: INCOMPLETE (lower bound pending), rerun'); sys.exit(3)
        ref = DoubleRef(0.5 * (float(Dd.Z2.lo) + float(Dd.Z2.hi)))
        x = lower_bound_point(ref)
        v6, t6 = point_F8bulk(Dd, x, 6.0, 2.45, 7.0)
        v8, t8 = point_F8bulk(Dd, x, 8.0, 2.3, 8.0)
        if not (v6.hi >= v8.lo and v8.hi >= v6.lo): raise FloatingPointError('lower bound: Ewald enclosures disagree')
        v = isect(v6, v8)
        R['lower'] = {'point': x, 'F8bulk_eta6': [float(v6.lo), float(v6.hi)], 'F8bulk_eta8': [float(v8.lo), float(v8.hi)],
                      'abs_lower': float(v.mig())}
        save()
    res = assemble(Dd, lead_st, s2_st, rest_cells, bulk_st, R['lower'])
    out = {'status': 'COMPLETE', 'A0': A0, 'pieces': A_PIECES, 'R_NF': R_NF, 'R_NR': R_NR, 'results': res,
           'lead': {k: lead_st[k] for k in ('L', 'U_current', 'Lpt', 'n_eval', 'seconds')},
           's2': {k: s2_st[k] for k in ('L', 'U_current', 'Lpt', 'n_eval', 'seconds')},
           'bulk': {k: bulk_st[k] for k in ('LF', 'LFpt', 'leafF', 'leafF_pt', 'leafR', 'leafR_pt', 'leafE', 'n_eval', 'n_tmjet', 'n_leaves', 'seconds')},
           'c6': R['c6'], 'check_c42': R['check_c42'], 'check_P6': R['data_check_P6'],
           'rest_cells': R['rest']}
    json.dump(out, open(out_path, 'w'), indent=1, default=float)
    log('main: COMPLETE -> %s (%.0f s this chunk)' % (out_path, time.time() - t0))
    sys.exit(0)

def controls(out_path):
    t_start = time.time()
    res = {'mode': 'controls', 'python_flint': __import__('flint').__version__, 'numpy': np.__version__}
    def log(s): print(s, flush=True)
    Dd = Data(log=log)
    res['data_seconds'] = Dd.seconds; res['c6'] = Dd.c6_arb.str(25); res['check_P6_far_moments_vs_c6'] = Dd.check_P6
    W, GJ, hol = group_GJ(); Lm = cone_normals(GJ)
    res['GJ'] = {'W': len(W), 'GJ': len(GJ), 'holomorphic': hol}
    # ---------------- J1: double-precision F8bulk at the 396 check points (DEFINITIONS L4 (e) CHECK (1))
    REF = 14569.136077584635          # <lab>/wp1/check_num.json 'bulk_cell'[0] (full precision)
    Z2m = 0.5 * (float(Dd.Z2.lo) + float(Dd.Z2.hi))
    ref = DoubleRef(Z2m)
    B = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [.5, .5, .5, .5]]) / 2.0
    rg = np.random.default_rng(7)
    cell = rg.random((400, 4)) @ B
    nb = np.vstack([np.zeros((1, 4)), NEAR, NEAR + NEAR[:1] * 0])
    near = np.array(sorted({tuple(np.round(v, 6)) for v in (np.array(list(itertools.product(range(-3, 4), repeat=4)), float) @ B) if np.linalg.norm(v) < 1.6}))
    dd = cell[:, None, :] - near[None]; j = (dd ** 2).sum(2).argmin(1); zc = cell - near[j]
    keep = np.linalg.norm(zc, axis=1) >= 0.1
    Fd = ref.Q(zc[keep])
    sup = float(np.abs(Fd).max())
    res['J1'] = {'n_points': int(keep.sum()), 'sup_double': sup, 'reference': REF, 'rel_diff': abs(sup / REF - 1),
                 'PASS_1e-6': bool(abs(sup / REF - 1) <= 1e-6)}
    log('J1 double F8bulk: n = %d, sup %.10g vs %.10g, rel %.2e -> %s' % (keep.sum(), sup, REF, abs(sup / REF - 1), res['J1']['PASS_1e-6']))
    # J1b: certified engine (point boxes) contains the double values at those points
    zs = zc[keep]
    ob = bulk_eval(Dd, zs, np.full(len(zs), 1e-7), np.zeros(len(zs)), want_R=False)
    tol = 1e-9 * np.abs(Fd)
    inside = (ob['Q'].lo <= Fd + tol) & (Fd - tol <= ob['Q'].hi)
    res['J1b_engine_contains_double'] = {'n': int(len(zs)), 'violations': int((~inside).sum()),
                                          'max_width': float((ob['Q'].hi - ob['Q'].lo).max()),
                                          'max_rel_mid_diff': float(np.max(np.abs(0.5 * (ob['Q'].lo + ob['Q'].hi) - Fd) / np.abs(Fd)))}
    log('J1b engine point enclosures vs double: %s' % json.dumps(res['J1b_engine_contains_double']))
    # ---------------- J2: Ewald eta = 6 vs eta = 8 (arb, independent psi_1 route) at 100 points; and vs the engine
    t0 = time.time()
    pts = zs[:100]
    n_ov = 0; n_eng = 0; widths = []
    for p in pts:
        x = [float(v) for v in p]
        v6, _ = point_F8bulk(Dd, x, 6.0, 2.45, 7.0)
        v8, _ = point_F8bulk(Dd, x, 8.0, 2.3, 8.0)
        e1 = bulk_eval(Dd, np.array([x]), np.array([1e-9]), np.array([0.0]), want_R=False)['Q'][0]
        ov = (v6.hi >= v8.lo) and (v8.hi >= v6.lo)
        n_ov += int(ov)
        n_eng += int((e1.hi >= v6.lo) and (v6.hi >= e1.lo) and (e1.hi >= v8.lo) and (v8.hi >= e1.lo))
        widths.append([float(v6.hi - v6.lo), float(v8.hi - v8.lo), float(e1.hi - e1.lo)])
    widths = np.array(widths)
    res['J2'] = {'n_points': len(pts), 'eta6_eta8_overlap': n_ov, 'engine_overlaps_both': n_eng,
                 'median_width_eta6_eta8_engine': np.median(widths, 0).tolist(), 'max_width': widths.max(0).tolist(),
                 'PASS': bool(n_ov == len(pts) and n_eng == len(pts)), 'seconds': time.time() - t0}
    log('J2 Ewald eta 6 vs 8 (arb) and engine: %s' % json.dumps(res['J2']))
    # ---------------- J3: mutations
    t0 = time.time()
    rngm = np.random.default_rng(13)
    # boxes on the shell |z| = r_c = 0.1 (bulk side) and near the sup point; reference values inside
    dirs = rngm.normal(size=(64, 4)); dirs /= np.linalg.norm(dirs, axis=1)[:, None]
    cs = np.round(dirs * 0.1 * 2 ** 20) / 2 ** 20; hs = np.full(len(cs), 2 ** -10)
    pin = cs + rngm.uniform(-1, 1, size=cs.shape) * hs[:, None]
    vref = ref.Q(pin)
    mut = {}
    for name in ('none', 'drop_l0_cancel', 'drop_N1'):
        MUT['drop_l0_cancel'] = name == 'drop_l0_cancel'; MUT['drop_N1'] = name == 'drop_N1'
        try:
            o = bulk_eval(Dd, cs, hs, np.zeros(len(cs)), want_R=False)
            viol = int(np.sum((vref < o['Q'].lo) | (vref > o['Q'].hi)))
            mut[name] = {'containment_violations': viol, 'sup_bound': float(o['Qmag'].max()),
                         'flag_blowup_over_1e5': bool(o['Qmag'].max() > 1e5)}
        finally:
            MUT['drop_l0_cancel'] = False; MUT['drop_N1'] = False
    mut['PASS'] = bool(mut['none']['containment_violations'] == 0 and not mut['none']['flag_blowup_over_1e5']
                       and (mut['drop_l0_cancel']['flag_blowup_over_1e5'] or mut['drop_l0_cancel']['containment_violations'] > 0)
                       and (mut['drop_N1']['flag_blowup_over_1e5'] or mut['drop_N1']['containment_violations'] > 0))
    # lead closed form vs finite differences of the leading function (and a mutated formula)
    c6f = 0.5 * (float(Dd.c6.lo) + float(Dd.c6.hi))
    def H6f(x):
        p2 = (x ** 2).sum(-1); m2 = (x ** 4).sum(-1); m3 = (x ** 6).sum(-1)
        return 6 * m3 - 7.5 * p2 * m2 + 1.875 * p2 ** 3
    def Gl(x, u):
        s = (x ** 2).sum(-1); th = 1 / np.sqrt(1 + u / s ** 2)
        return th * (-2.5 * c6f) * H6f(x) / s ** 3
    def Lfd(x, u, h=1e-4):
        E = np.eye(4); Hm = np.zeros((4, 4))
        for i in range(4):
            for jj in range(4):
                Hm[i, jj] = (Gl(x + h * E[i] + h * E[jj], u) - Gl(x + h * E[i] - h * E[jj], u) - Gl(x - h * E[i] + h * E[jj], u) + Gl(x - h * E[i] - h * E[jj], u)) / (4 * h * h)
        m11 = (Hm[0, 0] + Hm[1, 1]) / 4; m22 = (Hm[2, 2] + Hm[3, 3]) / 4
        m12 = ((Hm[0, 2] + Hm[1, 3]) + 1j * (Hm[0, 3] - Hm[1, 2])) / 4
        M = np.array([[m11, m12], [np.conj(m12), m22]])
        z = np.array([x[0] + 1j * x[1], x[2] + 1j * x[3]]); s = (abs(z) ** 2).sum(); th = 1 / np.sqrt(1 + u / s ** 2)
        g0 = (1 / th) * np.eye(2) + (th - 1 / th) * np.outer(np.conj(z), z) / s
        return np.sqrt(u) * np.abs(np.linalg.eigvals(np.linalg.solve(g0, M))).max()
    worst = 0.0
    for _ in range(40):
        x = rngm.normal(size=4); x /= np.linalg.norm(x); u = rngm.uniform(0.05, 1)
        worst = max(worst, abs(L_point(list(x), u, c6f) / Lfd(x, u) - 1))
    mut['lead_closed_form_vs_fd_max_rel'] = worst
    mut['PASS'] = bool(mut['PASS'] and worst < 1e-5)
    mut['seconds'] = time.time() - t0
    res['J3_mutations'] = mut
    log('J3 mutations: %s' % json.dumps(mut))
    # ---------------- NaN / inf guards
    g = {}
    def raises(fn):
        try:
            fn(); return False
        except (FloatingPointError, ValueError, AssertionError):
            return True
    g['check_finite_nan'] = raises(lambda: check_finite('t', I(np.array([np.nan]), np.array([1.0]))))
    g['check_finite_inf'] = raises(lambda: check_finite('t', I(np.array([0.0]), np.array([np.inf]))))
    g['a2f_nonfinite_arb'] = raises(lambda: a2f(arb('nan')))
    g['bulk_eval_nan_centre'] = raises(lambda: bulk_eval(Dd, np.array([[np.nan, 0.1, 0.1, 0.1]]), np.array([1e-4]), np.array([0.0])))
    g['near_taylor_on_lattice_point'] = raises(lambda: near_taylor(np.array([[0.5, 0, 0, 0]]), NEAR, 1, KP, np.array([1e-4])))
    g['near_taylor_box_too_large'] = raises(lambda: near_taylor(np.array([[0.1, 0, 0, 0]]), np.zeros((1, 4)), 3, KT, np.array([0.06])))
    g['irecip_zero'] = raises(lambda: irecip(I(np.array([-1.0]), np.array([1.0]))))
    g['gegen_tail_t_too_large'] = raises(lambda: gegen_tail(np.array([0.6]), 1, 2, 9))
    g['all_PASS'] = all(g.values())
    res['guards'] = g
    log('guards: %s' % json.dumps(g))
    # ---------------- J3b: the order-r^2 neck model against finite differences of the double reference Q(z; eps)
    t0 = time.time()
    NP = NeckPolys(Dd)
    res['check_c42'] = NP.check_c42
    def model(x, u, r, drop_P8=False):
        xc = [I(np.array([float(v)])) for v in x]
        Vc = imatmul(monomial_ranges(xc, 8), NP.M)
        Pc = {nm_: Vc[:, 5 * k] * (0.0 if drop_P8 else 1.0) for k, nm_ in enumerate(NP.names)}
        tr, disc = neck_parts(xc, I(np.array([u])), I(np.array([r])), Pc, Dd.c6, NP.c42, {'recip': irecip, 'sqrt': isqrt})
        return float((np.sqrt(u) * (abs(0.5 * (tr.lo + tr.hi)) * 0.5 + np.sqrt(0.5 * (disc.lo + disc.hi))))[0])
    def truth(x, u, r, h):
        eps = u * r ** 4; z0 = r * np.array(x); E = np.eye(4); Hm = np.zeros((4, 4))
        pts = [z0 + si * h * E[i] + sj * h * E[j] for i in range(4) for j in range(4) for si, sj in ((1, 1), (1, -1), (-1, 1), (-1, -1))]
        vals = ref.Q(np.array(pts), eps).reshape(4, 4, 4)
        for i in range(4):
            for j in range(4):
                Hm[i, j] = (vals[i, j, 0] - vals[i, j, 1] - vals[i, j, 2] + vals[i, j, 3]) / (4 * h * h)
        m11 = (Hm[0, 0] + Hm[1, 1]) / 4; m22 = (Hm[2, 2] + Hm[3, 3]) / 4
        m12 = ((Hm[0, 2] + Hm[1, 3]) + 1j * (Hm[0, 3] - Hm[1, 2])) / 4
        M = np.array([[m11, m12], [np.conj(m12), m22]])
        zz = np.array([z0[0] + 1j * z0[1], z0[2] + 1j * z0[3]]); s = (abs(zz) ** 2).sum(); th = s / np.sqrt(s * s + eps)
        g0 = (1 / th) * np.eye(2) + (th - 1 / th) * np.outer(np.conj(zz), zz) / s
        return np.sqrt(eps) * np.abs(np.linalg.eigvals(np.linalg.solve(g0, M))).max()
    j3b = {}
    rngb = np.random.default_rng(2)
    for r in (0.03, 0.06, 0.1):
        w2 = 0.0; w0 = 0.0; wm = 0.0
        for _ in range(6):
            x = rngb.normal(size=4); x /= np.linalg.norm(x); u = rngb.uniform(0.2, 1.0) * min(1, (A0 / r) ** 4)
            t = truth(x, u, r, r * 2e-3)
            w2 = max(w2, abs(model(x, u, r) / t - 1)); w0 = max(w0, abs(model(x, u, 0.0) / t - 1))
            wm = max(wm, abs(model(x, u, r, drop_P8=True) / t - 1))
        j3b['r=%g' % r] = {'order_r2_model': w2, 'leading_only': w0, 'mutated_drop_P8': wm}
    j3b['PASS'] = bool(j3b['r=0.03']['order_r2_model'] < 1e-3 and j3b['r=0.06']['order_r2_model'] < 1e-2
                       and j3b['r=0.06']['mutated_drop_P8'] > 3 * j3b['r=0.06']['order_r2_model'])
    j3b['seconds'] = time.time() - t0
    res['J3b_neck_model_vs_fd'] = j3b
    log('J3b neck model vs finite differences: %s' % json.dumps(j3b))
    # ---------------- J4 (diagnostic): not a certified integral in this build (see BUILD-WP3.md)
    # ---------------- timing of REDUCED covers (controls only)
    tim = {}
    t0 = time.time(); NT = neck_tables(Dd)
    cells = neck_cells()
    rc = neck_rest(Dd, NT, Lm, cells, nsub=1, log=log)
    tim['neck_rest_nsub1'] = {'seconds': time.time() - t0, 'max_Xr_by_rcell': {}, 'max_Cf_rle_RNF': max(v['Cf'] for k, v in rc.items() if k[1] <= R_NF + 1e-12)}
    for (a, b) in R_CELLS:
        tim['neck_rest_nsub1']['max_Xr_by_rcell']['%g-%g' % (a, b)] = max(v['Xr'] for k, v in rc.items() if k[0] == a)
    log('timing neck rest (nsub 1): %s' % json.dumps(tim['neck_rest_nsub1']))
    rest = make_rest_lookup(rc)
    t0 = time.time(); stn = lead_new_state(Lm)
    stn = lead_run(Lm, Dd.c6, stn, 1e-3, 100.0, log=log)
    t1 = time.time(); s2 = neck2_new_state(Lm, 16, pieces=[(1.0, 1.0)], ny=8)
    s2 = neck2_run(NP, Lm, Dd.c6, s2, 3e-2, 60.0, log=log, lead_w=0.0)
    Xmax = {str(p): max(math.sqrt(k[3]) * v['Xr'] for k, v in rc.items() if admissible(k, p)) for p in A_PIECES}
    Lh = max(stn['U_current'], stn['L']); S2h = max(s2['U_current'], s2['L'])
    tim['lead_tol1e-3_100s_and_S2_tol3e-2_60s'] = {'lead_done': stn['done'], 'L': stn['L'], 'U': stn['U_current'], 'lead_evals': stn['n_eval'],
        'lead_seconds': t1 - t0, 'lead_argmax(face,y,u)': stn['Lpt'], 'S2_done': s2['done'], 'S2_L': s2['L'], 'S2_U': s2['U_current'],
        'S2_evals': s2['n_eval'], 'S2_seconds': time.time() - t1, 'sqrt_u_rest_max_by_piece': Xmax,
        'R_neck_preview_by_piece': {str(p): Lh + p[1] ** 2 * S2h + Xmax[str(p)] for p in A_PIECES}}
    log('timing lead + S2: %s' % json.dumps(tim['lead_tol1e-3_100s_and_S2_tol3e-2_60s']))
    t0 = time.time(); stb = bulk_new_state(Lm)
    stb = bulk_run(Dd, Lm, stb, 3e-2, RCAP_FRAC * stn['L'], 80.0, log=log)
    tim['bulk_reduced_tolF_3e-2_100s'] = {k: stb[k] for k in ('done', 'LF', 'LFpt', 'leafF', 'leafR', 'leafE', 'n_eval', 'n_tmjet', 'n_leaves', 'it')}
    tim['bulk_reduced_tolF_3e-2_100s']['active_left'] = int(len(stb['h']))
    tim['bulk_reduced_tolF_3e-2_100s']['active_hist_h'] = {('%.2e' % k): int(v) for k, v in zip(*np.unique(stb['h'], return_counts=True))}
    tim['bulk_reduced_tolF_3e-2_100s']['seconds'] = time.time() - t0
    log('timing bulk (reduced): %s' % json.dumps(tim['bulk_reduced_tolF_3e-2_100s']))
    res['timing_reduced'] = tim
    try:
        import psutil; res['peak_rss_GB_process'] = psutil.Process().memory_info().peak_wset / 1e9
    except Exception:
        res['peak_rss_GB_process'] = None
    res['seconds_total'] = time.time() - t_start
    res['ALL_CONTROLS_PASS'] = bool(res['J1']['PASS_1e-6'] and res['J1b_engine_contains_double']['violations'] == 0
                                    and res['J2']['PASS'] and res['J3_mutations']['PASS'] and res['J3b_neck_model_vs_fd']['PASS'] and g['all_PASS'])
    json.dump(res, open(out_path, 'w'), indent=1, default=float)
    log('controls written to %s (%.0f s); ALL_CONTROLS_PASS = %s' % (out_path, res['seconds_total'], res['ALL_CONTROLS_PASS']))

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--controls', metavar='OUT_JSON'); g.add_argument('--main', metavar='OUT_JSON')
    ap.add_argument('--budget', type=float, default=530.0, help='seconds per --main invocation (capped at 530)')
    a = ap.parse_args()
    if a.controls:
        controls(a.controls)
    else:
        main(a.main, min(a.budget, 530.0))
