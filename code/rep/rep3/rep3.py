"""rep3.py -- blind replication seat R3 (K3 Kummer approximant, one bubble core).

Object (DEFINITIONS-AND-LEMMAS-v2 secs 1-2, window changed to [a/4, a]):
  Phi_a = phi_EH,a(|z|^2) + psi - chi_a psi_oo,  chi_a = 1 - st((4|z|/a - 1)/3)
  F = log det g  (g = dd^bar Phi_a in z; det g_EH = 1),  Ric = -i dd^bar F,
  |Ric|_g = spectral radius (operator norm w.r.t. g) of the endomorphism g^{-1} Ric.
Scaled chart zeta = z/a: chart potential a^-2 Phi_a(a zeta) = phi_1(sigma) + a^8 Ph(zeta),
  Ph = a^-10 [P_a(a zeta) - P_a(0)],  P_a = psi - chi_a psi_oo.
  q = det(g)/det(g0) - 1 ; F/a^8 = log1p(q)/a^8 ; |Ric|_g/a^6 = chart spectral radius / a^8.

Method (seat R3): exact Taylor polynomial of psi in complex monomials z^alpha zbar^beta
(lattice sums, closed-form Wirtinger pairing formula); odd-odd part = monomials with |alpha|,|beta|
odd; 4-jets of the FULL potential by multivariate Taylor arithmetic in (z1, z2, w1, w2), w = zbar;
log det(g0 + a^8 E) = log det g0 + log(1 + a^8 tr(adj(g0)E)/det g0 + a^16 det E/det g0) (exact 2x2
identity; det g0 computed); Ric = -dd^bar of that jet. Near E: resolved chart u = zeta1^2, w = zeta2/zeta1.
Rigorous part: same pipeline with first-order affine (Taylor-model) complex numbers on boxes.

Usage:  python rep3.py            all stages (resumes from cache/), writes rep3_result.json
        python rep3.py <stage>    stage in: coeffs controls sample interval report
"""
import sys, os, json, time, math, itertools
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "cache")
os.makedirs(CACHE, exist_ok=True)
V = 1.0 / 32.0

def log(*a):
    print(*a, flush=True)

# ==================================================================================
# 1. lattice L = Hurwitz/2, constants
# ==================================================================================
def lattice_points(R, exclude_zero=True):
    M = int(math.ceil(2 * R)) + 1
    r = np.arange(-M, M + 1, dtype=np.float64)
    out = []
    for shift in (0.0, 0.5):
        g = np.stack(np.meshgrid(r, r, r, r, indexing="ij"), -1).reshape(-1, 4)
        p = (g + shift) / 2.0
        n2 = (p * p).sum(1)
        keep = n2 <= R * R + 1e-12
        if exclude_zero:
            keep &= n2 > 1e-12
        out.append(p[keep])
    return np.concatenate(out)

def Zkappa(kappa):
    import mpmath as mp
    k = mp.mpf(kappa)
    return 24 * mp.mpf(4) ** (2 * k) * (1 - mp.mpf(2) ** (1 - 2 * k)) * mp.zeta(2 * k) * mp.zeta(2 * k - 1)

def binom_half(k):
    b = 1.0
    for j in range(k):
        b *= (0.5 - j) / (j + 1)
    return b

def A4_tail(T, p):
    """sum_{l in L, |l| >= T} |l|^-p <= 1024 sum_j C(3,j) 4^(j-3) (T-1/2)^(j+1-p)/(p-j-1)  (spec A4, x=0)."""
    return 1024.0 * sum(math.comb(3, j) * 4.0 ** (j - 3) * (T - 0.5) ** (j + 1 - p) / (p - j - 1) for j in range(4))

# ==================================================================================
# 2. Taylor coefficients of |z - l|^(-2 nu) at 0 (REP3-NOTES sec 2):
#    coeff of z^alpha w^beta = sum_{gamma<=min(alpha,beta)} f^(n)(|l|^2) (-lbar)^(alpha-gamma) (-l)^(beta-gamma)
#                              / (gamma! (alpha-gamma)! (beta-gamma)!),  n = |alpha|+|beta|-|gamma|, f = s^-nu
# ==================================================================================
def monomials(d):
    out = []
    for a1 in range(d + 1):
        for a2 in range(d + 1 - a1):
            for b1 in range(d + 1 - a1 - a2):
                out.append((a1, a2, b1, d - a1 - a2 - b1))
    return out

def lattice_coeffs(nu, d, pts, chunk=200000):
    from math import factorial as fa
    monos = monomials(d)
    acc = {m: 0j for m in monos}
    for c0 in range(0, len(pts), chunk):
        p = pts[c0:c0 + chunk]
        l1 = p[:, 0] + 1j * p[:, 1]; l2 = p[:, 2] + 1j * p[:, 3]
        s0 = (p * p).sum(1)
        A1 = [np.ones_like(l1)]; A2 = [np.ones_like(l1)]; B1 = [np.ones_like(l1)]; B2 = [np.ones_like(l1)]
        for _ in range(d):
            A1.append(A1[-1] * (-np.conj(l1))); A2.append(A2[-1] * (-np.conj(l2)))
            B1.append(B1[-1] * (-l1)); B2.append(B2[-1] * (-l2))
        F = []
        for n in range(d + 1):
            rising = 1.0
            for j in range(n):
                rising *= (nu + j)
            F.append(((-1) ** n) * rising * s0 ** (-nu - n))
        for m in monos:
            a1, a2, b1, b2 = m
            tot = 0j
            for g1 in range(min(a1, b1) + 1):
                for g2 in range(min(a2, b2) + 1):
                    n = d - g1 - g2
                    cf = 1.0 / (fa(g1) * fa(g2) * fa(a1 - g1) * fa(a2 - g2) * fa(b1 - g1) * fa(b2 - g2))
                    tot += cf * np.sum(F[n] * A1[a1 - g1] * A2[a2 - g2] * B1[b1 - g1] * B2[b2 - g2])
            acc[m] += tot
    return acc

def radius_for(nu, d):
    return 8.0 if (nu == 1 and d <= 8) else 5.0

DMAX = 18   # explicit Taylor degree, sampling
DRIG = 12   # explicit Taylor degree, rigorous stage (higher degrees -> rigorous tail)
KMAX = 3    # eps^k pieces explicit

def stage_coeffs():
    fn = os.path.join(CACHE, "coeffs.json")
    if os.path.exists(fn):
        return json.load(open(fn))
    t0 = time.time()
    pts_big = lattice_points(8.0)
    pts_small = pts_big[(pts_big ** 2).sum(1) <= 25.0 + 1e-12]
    log("lattice points R=8:", len(pts_big), " R=5:", len(pts_small))
    table = {}
    for k in range(1, KMAX + 1):
        nu = 1 if k == 1 else 2 * k - 1
        pref = -0.5 if k == 1 else binom_half(k) / (1 - 2 * k)
        for d in range(2 if k >= 2 else 4, DMAX + 1, 2):
            R = radius_for(nu, d)
            c = lattice_coeffs(nu, d, pts_big if R == 8.0 else pts_small)
            table["%d_%d" % (k, d)] = {",".join(map(str, m)): [pref * v.real, pref * v.imag] for m, v in c.items()}
            log("  k=%d d=%2d R=%g done  (%.1fs)" % (k, d, R, time.time() - t0))
    json.dump(table, open(fn, "w"))
    return table

def load_poly(table, k, d):
    tab = table["%d_%d" % (k, d)]
    exps = np.array([list(map(int, s.split(","))) for s in tab.keys()], dtype=np.int64)
    cf = np.array([v[0] + 1j * v[1] for v in tab.values()])
    return exps, cf

def H6_poly():
    """H6 = (3/8) Re[zb1 z1^5 + zb2 z2^5 - 5 zb1 z1 z2^4 - 5 zb2 z1^4 z2] (spec L3(e)), exps (z1,z2,w1,w2)."""
    t = {(5, 0, 1, 0): 3 / 16, (1, 0, 5, 0): 3 / 16, (0, 5, 0, 1): 3 / 16, (0, 1, 0, 5): 3 / 16,
         (1, 4, 1, 0): -15 / 16, (1, 0, 1, 4): -15 / 16, (4, 1, 0, 1): -15 / 16, (0, 1, 4, 1): -15 / 16}
    return np.array(list(t.keys()), dtype=np.int64), np.array(list(t.values()), dtype=complex)

def absz4_poly():
    t = {(2, 0, 2, 0): 1.0, (0, 2, 0, 2): 1.0, (1, 1, 1, 1): 2.0}
    return np.array(list(t.keys()), dtype=np.int64), np.array(list(t.values()), dtype=complex)

def c6_enclosure(R=16.0):
    """c6 = (8/3) sum_{l!=0} U6(l0/|l|)/|l|^8 ; finite sum + rigorous tail (8/3)*7*A4(R,8)."""
    fn = os.path.join(CACHE, "c6.json")
    if os.path.exists(fn):
        return json.load(open(fn))
    M = int(2 * R) + 1
    tot = 0.0; absum = 0.0
    r = np.arange(-M, M + 1, dtype=np.float64)
    for shift in (0.0, 0.5):
        for i0 in r:
            g = np.stack(np.meshgrid(r, r, r, indexing="ij"), -1).reshape(-1, 3)
            p = np.concatenate([np.full((len(g), 1), i0), g], 1)
            p = (p + shift) / 2.0
            n2 = (p * p).sum(1)
            keep = (n2 <= R * R) & (n2 > 1e-12)
            p = p[keep]; n2 = n2[keep]
            t2 = p[:, 0] ** 2 / n2
            U6 = 64 * t2 ** 3 - 80 * t2 ** 2 + 24 * t2 - 1
            terms = U6 / n2 ** 4
            tot += terms.sum(); absum += np.abs(terms).sum()
    c6 = 8.0 / 3.0 * tot
    rad = 8.0 / 3.0 * (7.0 * A4_tail(R, 8) + 1e-12 * absum)
    out = {"c6_mid": c6, "c6_rad": rad, "R": R}
    json.dump(out, open(fn, "w"))
    return out

# ==================================================================================
# 3. Interval helpers (outward rounding by nextafter; exp padded by 2 ulp)
# ==================================================================================
def _dn(x): return np.nextafter(x, -np.inf)
def _up(x): return np.nextafter(x, np.inf)
def iv(lo, hi=None):
    lo = np.asarray(lo, dtype=np.float64)
    return (lo, lo if hi is None else np.asarray(hi, dtype=np.float64))
def iadd(a, b): return (_dn(a[0] + b[0]), _up(a[1] + b[1]))
def isub(a, b): return (_dn(a[0] - b[1]), _up(a[1] - b[0]))
def imul(a, b):
    p = np.stack(np.broadcast_arrays(a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1]))
    return (_dn(p.min(0)), _up(p.max(0)))
def iscal(c, a):
    c = float(c)
    return (_dn(c * a[0]), _up(c * a[1])) if c >= 0 else (_dn(c * a[1]), _up(c * a[0]))
def irec(a):
    if np.any((a[0] <= 0) & (a[1] >= 0)):
        raise ValueError("irec: interval contains 0")
    return (_dn(1.0 / a[1]), _up(1.0 / a[0]))
def isqrt(a): return (_dn(np.sqrt(np.maximum(a[0], 0.0))), _up(np.sqrt(a[1])))
def iexp(a): return (_dn(_dn(np.exp(a[0]))), _up(_up(np.exp(a[1]))))
def ihull(a, b): return (np.minimum(a[0], b[0]), np.maximum(a[1], b[1]))

def ut_var(lo, hi, K):
    z = np.zeros_like(lo)
    return [iv(lo, hi), iv(np.ones_like(lo))] + [iv(z) for _ in range(K - 1)]
def ut_const(c, like, K):
    z = np.zeros_like(like)
    return [iv(np.full_like(like, c))] + [iv(z) for _ in range(K)]
def ut_add(a, b): return [iadd(x, y) for x, y in zip(a, b)]
def ut_sub(a, b): return [isub(x, y) for x, y in zip(a, b)]
def ut_scal(c, a): return [iscal(c, x) for x in a]
def ut_mul(a, b):
    K = len(a) - 1; out = []
    for k in range(K + 1):
        s = imul(a[0], b[k])
        for i in range(1, k + 1):
            s = iadd(s, imul(a[i], b[k - i]))
        out.append(s)
    return out
def ut_rec(a):
    K = len(a) - 1; y0 = irec(a[0]); y = [y0]
    for k in range(1, K + 1):
        s = imul(a[1], y[k - 1])
        for i in range(2, k + 1):
            s = iadd(s, imul(a[i], y[k - i]))
        y.append(iscal(-1.0, imul(s, y0)))
    return y
def ut_sqrt(a):
    K = len(a) - 1; y0 = isqrt(a[0]); y = [y0]; inv2 = irec(iscal(2.0, y0))
    for k in range(1, K + 1):
        s = a[k]
        for i in range(1, k):
            s = isub(s, imul(y[i], y[k - i]))
        y.append(imul(s, inv2))
    return y
def ut_exp(a):
    K = len(a) - 1; y = [iexp(a[0])]
    for k in range(1, K + 1):
        s = imul(a[1], y[k - 1])
        for i in range(2, k + 1):
            s = iadd(s, iscal(float(i), imul(a[i], y[k - i])))
        y.append(imul(s, iv(_dn(np.full_like(a[0][0], 1.0 / k)), _up(np.full_like(a[0][0], 1.0 / k)))))
    return y
def ut_compose(Tj, fk):
    K = len(Tj) - 1
    z = np.zeros_like(Tj[0][0])
    H = [iv(z)] + Tj[1:]
    res = [fk[0]] + [iv(z) for _ in range(K)]
    P = H
    for j in range(1, K + 1):
        res = ut_add(res, [imul(fk[j], p) for p in P])
        P = ut_mul(P, H)
    return res

# ==================================================================================
# 4. cut-off profile c(sigma) = 1 - chi_a = st(T), T = (sqrt(sigma) - 1/4)/(3/4)
# ==================================================================================
WIN_LO, WIN_HI = 0.25, 1.0
TE = 0.01
KC = 5

def _S_polys(mmax):
    P = np.polynomial.polynomial
    S = [None, np.array([1.0])]
    for m in range(1, mmax):
        S.append(P.polyadd(P.polymul([1.0, -2.0], S[m]), P.polymul([0.0, 1.0, -1.0], P.polyder(S[m]))))
    return S

def S_sup_bounds(mmax=KC, ncell=4000):
    S = _S_polys(mmax); out = [None]
    e = np.linspace(0.0, 1.0, ncell + 1); U = (e[:-1], e[1:])
    for m in range(1, mmax + 1):
        c = S[m]
        acc = iv(np.full(ncell, c[-1]))
        for co in c[-2::-1]:
            acc = iadd(imul(acc, U), iv(np.full(ncell, co)))
        out.append(float(np.maximum(np.abs(acc[0]), np.abs(acc[1])).max()) * (1 + 1e-12))
    return out

def _lah(k, m):
    return math.comb(k - 1, m - 1) * math.factorial(k) // math.factorial(m)

_ENDB = None
def end_bounds():
    """b_k >= sup_{0<t<=TE} |st^(k)(t)|, k=0..KC."""
    global _ENDB
    if _ENDB is None:
        M = S_sup_bounds(KC)
        b = [math.e ** 2 * math.exp(-1.0 / TE)]
        for k in range(1, KC + 1):
            s = 0.0
            for m in range(1, k + 1):
                assert TE <= 1.0 / (k + m)
                s += M[m] * 2 ** m * _lah(k, m) * TE ** (-(k + m))
            b.append(math.e ** 2 * math.exp(-1.0 / TE) * s * (1 + 1e-12))
        _ENDB = (b, M)
    return _ENDB[0]

def _sig_of_T(T):
    return (WIN_LO + (WIN_HI - WIN_LO) * T) ** 2

def c_jet_interval(slo, shi, K=KC):
    """enclosures of c_k = c^(k)(sigma)/k!, k=0..K, valid for all sigma in [slo, shi]."""
    n = slo.shape
    def Tjet(lo, hi):
        sq = ut_sqrt(ut_var(lo, hi, K))
        return ut_scal(1.0 / (WIN_HI - WIN_LO), ut_sub(sq, ut_const(WIN_LO, lo, K)))
    def merge(cur, new, mask):
        if not np.any(mask):
            return cur
        newm = [(np.where(mask, x[0], np.inf), np.where(mask, x[1], -np.inf)) for x in new]
        return newm if cur is None else [ihull(a, b) for a, b in zip(cur, newm)]
    s_te, s_1te = _sig_of_T(TE), _sig_of_T(1 - TE)
    s0, s1 = WIN_LO ** 2, WIN_HI ** 2
    zero = np.zeros(n); out = None
    out = merge(out, [iv(zero)] * (K + 1), slo <= s0)
    out = merge(out, [iv(np.ones(n))] + [iv(zero)] * K, shi >= s1)
    b = end_bounds()
    sym = [iv(np.full(n, -b[j] / math.factorial(j)), np.full(n, b[j] / math.factorial(j))) for j in range(1, K + 1)]
    lo = np.maximum(slo, s0); hi = np.minimum(shi, s_te * (1 + 1e-12)); m = lo <= hi
    if np.any(m):
        lo2 = np.where(m, lo, s_te); hi2 = np.where(m, hi, s_te)
        out = merge(out, ut_compose(Tjet(lo2, hi2), [iv(zero, np.full(n, b[0]))] + sym), m)
    lo = np.maximum(slo, s_1te * (1 - 1e-12)); hi = np.minimum(shi, s1); m = lo <= hi
    if np.any(m):
        lo2 = np.where(m, lo, s_1te); hi2 = np.where(m, hi, s_1te)
        out = merge(out, ut_compose(Tjet(lo2, hi2), [iv(np.full(n, 1.0 - b[0]), np.ones(n))] + sym), m)
    lo = np.maximum(slo, s_te * (1 - 1e-12)); hi = np.minimum(shi, s_1te * (1 + 1e-12)); m = lo <= hi
    if np.any(m):
        mid = 0.5 * (s_te + s_1te)
        lo2 = np.where(m, lo, mid); hi2 = np.where(m, hi, mid)
        T = Tjet(lo2, hi2)
        g = ut_sub(ut_rec(T), ut_rec(ut_sub(ut_const(1.0, lo2, K), T)))
        st = ut_rec(ut_add(ut_const(1.0, lo2, K), ut_exp(g)))
        out = merge(out, st, m)
    return out

def p1_jet_interval(slo, shi, K=KC):
    """Taylor coefficients of p(sigma) = sqrt(sigma^2+1)/sigma = phi_1'(sigma)."""
    x = ut_var(slo, shi, K)
    fh = ut_sqrt(ut_add(ut_mul(x, x), ut_const(1.0, slo, K)))
    return ut_mul(fh, ut_rec(x))

class RadialTable:
    """rigorous enclosures of c_k and p_k (k = 0..KC) on fine sigma cells; range queries by sparse tables
    (naive interval Taylor arithmetic over a wide sigma range overestimates badly; on cells of width 1e-5 it is tight)."""
    def __init__(self, s0=0.02, s1=2.0, ds=1e-5):
        self.s0, self.ds = s0, ds
        n = int(math.ceil((s1 - s0) / ds))
        e = s0 + ds * np.arange(n + 1)
        lo, hi = e[:-1], e[1:]
        cj = c_jet_interval(lo, hi); pj = p1_jet_interval(lo, hi)
        self.n = n
        self.tabs = []
        for J in (cj, pj):
            levels = []
            for (l, h) in J:
                L = [(l, h)]
                w = 1
                while 2 * w <= n:
                    pl, ph = L[-1]
                    L.append((np.minimum(pl[:-w], pl[w:]), np.maximum(ph[:-w], ph[w:])))
                    w *= 2
                levels.append(L)
            self.tabs.append(levels)
    def _q(self, levels, i0, i1):
        out = []
        ln = i1 - i0 + 1
        j = np.floor(np.log2(ln)).astype(int)
        for L in levels:
            lo = np.empty(len(i0)); hi = np.empty(len(i0))
            for jj in np.unique(j):
                m = j == jj
                a = i0[m]; b = i1[m] - (1 << jj) + 1
                lo[m] = np.minimum(L[jj][0][a], L[jj][0][b]); hi[m] = np.maximum(L[jj][1][a], L[jj][1][b])
            out.append((lo, hi))
        return out
    def ranges(self, slo, shi):
        i0 = np.floor((slo - self.s0) / self.ds).astype(int) - 1
        i1 = np.floor((shi - self.s0) / self.ds).astype(int) + 1
        if np.any(i0 < 0) or np.any(i1 >= self.n):
            raise ValueError("sigma outside radial table")
        return self._q(self.tabs[0], i0, i1), self._q(self.tabs[1], i0, i1)

_RT = None
def radial_table():
    global _RT
    if _RT is None:
        _RT = RadialTable()
    return _RT

# ==================================================================================
# 5. affine complex numbers: (c, r): c[..., 0, :] centre, c[..., 1:, :] noise coefficients
#    (box coordinates e_k in [-1,1]), r >= 0 disc radius.  Point mode: one slot, r = 0.
# ==================================================================================
U_RND = 1.2e-16

def af_mag(x):
    c, r = x
    return np.abs(c[..., 0, :]) + np.abs(c[..., 1:, :]).sum(-2) + r

def af_mul(x, y):
    c1, r1 = x; c2, r2 = y
    a0 = c1[..., 0:1, :]; b0 = c2[..., 0:1, :]
    cen = a0 * b0
    lin = a0 * c2[..., 1:, :] + b0 * c1[..., 1:, :]
    if lin.shape[-2] == 0:
        c = cen
    else:
        shp = np.broadcast_shapes(cen.shape[:-2], lin.shape[:-2])
        c = np.concatenate([np.broadcast_to(cen, shp + cen.shape[-2:]), np.broadcast_to(lin, shp + lin.shape[-2:])], axis=-2)
    s1 = np.abs(c1[..., 1:, :]).sum(-2); s2 = np.abs(c2[..., 1:, :]).sum(-2)
    R = np.abs(a0[..., 0, :]) * r2 + np.abs(b0[..., 0, :]) * r1 + (s1 + r1) * (s2 + r2)
    mag = np.abs(cen[..., 0, :]) + np.abs(lin).sum(-2)
    return c, R * (1 + 8 * U_RND) + 8 * U_RND * mag

def af_add(x, y):
    c = x[0] + y[0]
    return c, x[1] + y[1] + 2 * U_RND * af_mag((c, 0.0))
def af_sub(x, y):
    c = x[0] - y[0]
    return c, x[1] + y[1] + 2 * U_RND * af_mag((c, 0.0))
def af_scale(x, s):
    c = x[0] * s
    return c, x[1] * abs(s) + 2 * U_RND * af_mag((c, 0.0))
def af_iscale(x, siv):
    """x times a real interval constant given as (mid, rad)."""
    m, rd = siv
    c = x[0] * m
    return c, x[1] * abs(m) + rd * af_mag(x) + 2 * U_RND * af_mag((c, 0.0))
def af_conj(x): return np.conj(x[0]), x[1]
def af_const(val, nn, n):
    c = np.zeros((nn, n), complex); c[0] = val
    return c, np.zeros(n)
def af_real_range(x):
    c, r = x
    m = c[..., 0, :].real
    w = np.abs(c[..., 1:, :].real).sum(-2) + r
    return m - w, m + w
def af_from_interval_and_slope(cen_iv, slope_iv, D):
    """value = [cen_iv] + slope*(sigma - sigma_c), slope in slope_iv, D = sigma - sigma_c (affine, centre 0)."""
    cm = 0.5 * (cen_iv[0] + cen_iv[1]); cr = 0.5 * (cen_iv[1] - cen_iv[0])
    sm = 0.5 * (slope_iv[0] + slope_iv[1]); sr = 0.5 * (slope_iv[1] - slope_iv[0])
    c = D[0] * sm
    c[0] = cm
    Dmag = af_mag(D)
    r = np.abs(sm) * D[1] + sr * Dmag + cr * (1 + 4 * U_RND) + 4 * U_RND * (np.abs(cm) + np.abs(sm) * Dmag)
    return c, r
TIGHT_FAC = 0.3
def tighter(x, rng):
    """per box, the tighter of two valid enclosures of the same real quantity: the mean-value affine form x
    or the plain range interval rng (from the radial table)."""
    c, r = x
    wA = np.abs(c[1:].real).sum(0) + r
    m = 0.5 * (rng[0] + rng[1]); rd = 0.5 * (rng[1] - rng[0])
    rd = rd * (1 + 4 * U_RND) + 4 * U_RND * np.abs(m)
    use = rd < TIGHT_FAC * wA
    c2 = c.copy(); c2[:, use] = 0; c2[0, use] = m[use]
    r2 = np.where(use, rd, r)
    return c2, r2

def af_func_real(x, f, dfrange):
    """f applied to a real-valued affine x by the mean-value form; dfrange(lo, hi) encloses f' on [lo,hi]."""
    c, r = x
    c0 = c[0].real
    lo, hi = af_real_range(x)
    dlo, dhi = dfrange(lo, hi)
    D = (np.concatenate([np.zeros_like(c[:1]), c[1:].real.astype(complex)], 0), r)
    fc = f(c0)
    return af_from_interval_and_slope((_dn(_dn(fc)), _up(_up(fc))), (dlo, dhi), D)

# ==================================================================================
# 6. jets in nv variables, order <= N, affine coefficients: (c (nc, nn, n), r (nc, n))
# ==================================================================================
def mindices(nv, N):
    out = []
    for o in range(N + 1):
        for t in itertools.product(range(o + 1), repeat=nv):
            if sum(t) == o:
                out.append(t)
    return out

class JetSpace:
    def __init__(self, N, nv=4):
        self.N = N; self.nv = nv
        self.idx = mindices(nv, N)
        self.pos = {t: i for i, t in enumerate(self.idx)}
        self.nc = len(self.idx)
        I, J, K = [], [], []
        for i, a in enumerate(self.idx):
            for j, b in enumerate(self.idx):
                s = tuple(x + y for x, y in zip(a, b))
                if sum(s) <= N:
                    I.append(i); J.append(j); K.append(self.pos[s])
        o = np.argsort(K, kind="stable")
        self.I = np.array(I)[o]; self.J = np.array(J)[o]; self.K = np.array(K)[o]
        self.starts = np.searchsorted(self.K, np.arange(self.nc))
        self.cnt = np.bincount(self.K, minlength=self.nc)
    def unit(self, v):
        t = [0] * self.nv; t[v] = 1
        return self.pos[tuple(t)]
    def mul(self, A, B):
        c, R = af_mul((A[0][self.I], A[1][self.I]), (B[0][self.J], B[1][self.J]))
        cs = np.add.reduceat(c, self.starts, axis=0)
        Rs = np.add.reduceat(R, self.starts, axis=0)
        mag = np.add.reduceat(af_mag((c, R)), self.starts, axis=0)
        return cs, Rs + mag * U_RND * (self.cnt[:, None] + 2)
    def scal(self, s, A):
        return af_mul((s[0][None], s[1][None]), A)
    def zeros(self, nn, n):
        return np.zeros((self.nc, nn, n), complex), np.zeros((self.nc, n))
    def const_part(self, A):
        return A[0][0], A[1][0]
    def nonconst(self, A):
        c = A[0].copy(); r = A[1].copy(); c[0] = 0; r[0] = 0
        return c, r
    def series(self, fk, H):
        """sum_k fk[k] H^k, H with zero constant term, fk affine scalars (fk[0] may be None)."""
        nn, n = H[0].shape[1], H[0].shape[2]
        acc = self.zeros(nn, n)
        if fk[0] is not None:
            acc[0][0] += fk[0][0]; acc[1][0] += fk[0][1]
        P = H
        for k in range(1, len(fk)):
            if k > 1:
                P = self.mul(P, H)
            acc = af_add(acc, self.scal(fk[k], P))
        return acc
    def series_pw(self, fk, Ppow):
        nn, n = Ppow[1][0].shape[1], Ppow[1][0].shape[2]
        acc = self.zeros(nn, n)
        if fk[0] is not None:
            acc[0][0] += fk[0][0]; acc[1][0] += fk[0][1]
        for k in range(1, len(fk)):
            acc = af_add(acc, self.scal(fk[k], Ppow[k]))
        return acc

def mixed_maps(JS4, JS2, off=2):
    """(i,j) -> (source indices in JS4, factors) giving the order-2 jet of d_{v_i} d_{v_{off+j}} f."""
    maps = {}
    for i in range(2):
        for j in range(2):
            src, fac = [], []
            for g in JS2.idx:
                t = list(g); t[i] += 1; t[off + j] += 1
                src.append(JS4.pos[tuple(t)]); fac.append((g[i] + 1) * (g[off + j] + 1))
            maps[(i, j)] = (np.array(src), np.array(fac, dtype=float))
    return maps

def extract(J4, mp):
    src, fac = mp
    return J4[0][src] * fac[:, None, None], J4[1][src] * fac[:, None]

# ==================================================================================
# 7. polynomials (scaled units) and their jets via a Taylor-shift matrix
# ==================================================================================
def build_polys(table, a_lo, a_hi, Dexp, rigorous):
    """Ph = c(sigma) A + B ; A = odd-odd monomials, B = even-even; coefficients (mid, rad) for a in [a_lo,a_hi]."""
    am = 0.5 * (a_lo + a_hi)
    acc = {}
    def put(ex, cm, cr):
        for e, m_, r_ in zip(map(tuple, ex), cm, cr):
            if e in acc:
                acc[e][0] += m_; acc[e][1] += r_
            else:
                acc[e] = [m_, r_]
    info = []
    c6 = c6_enclosure() if rigorous else None
    for k in range(1, KMAX + 1):
        for d in range(4, Dexp + 1, 2):
            if k == 1 and d == 4:
                continue                                   # P_4[psi_1] = 0 (spec L3(d))
            e = 4 * k + d - 10
            ex, cf = load_poly(table, k, d)
            isoo = ((ex[:, 0] + ex[:, 1]) % 2) == 1
            keep = isoo if d % 4 == 2 else ~isoo
            ex, cf = ex[keep], cf[keep]
            if rigorous and (k, d) == (1, 6):
                ex, h6 = H6_poly()
                put(ex, -0.5 * c6["c6_mid"] * h6, 0.5 * c6["c6_rad"] * np.abs(h6) + 1e-15 * np.abs(h6) * c6["c6_mid"])
                continue
            if rigorous and (k, d) == (2, 4):
                ex, z4 = absz4_poly()
                c42 = float(Zkappa(2.5)) / 4.0
                put(ex, c42 * z4 * am ** e, c42 * np.abs(z4) * (max(a_hi ** e - am ** e, am ** e - a_lo ** e) + 1e-14 * am ** e))
                continue
            span = max(abs(a_hi ** e - am ** e), abs(am ** e - a_lo ** e))
            put(ex, cf * am ** e, np.abs(cf) * (span + 1e-9 * am ** e) if rigorous else np.zeros(len(cf)))
    keys = list(acc.keys())
    ex = np.array(keys, dtype=np.int64)
    cm = np.array([acc[k][0] for k in keys]); cr = np.array([acc[k][1] for k in keys])
    oo = ((ex[:, 0] + ex[:, 1]) % 2) == 1
    return {"A": (ex[oo], cm[oo], cr[oo]), "B": (ex[~oo], cm[~oo], cr[~oo])}

def shift_matrix(JS, ex, cm, cr):
    """jet coefficient gamma = sum_m coef_m prod binom(m_i, gamma_i) x^(m - gamma)."""
    cols = {}; rows = []; vals_m = []; vals_r = []; cidx = []
    for m, a, b in zip(map(tuple, ex), cm, cr):
        for gi, g in enumerate(JS.idx):
            if all(gg <= mm for gg, mm in zip(g, m)):
                mp_ = tuple(mm - gg for gg, mm in zip(g, m))
                f = 1.0
                for gg, mm in zip(g, m):
                    f *= math.comb(mm, gg)
                if mp_ not in cols:
                    cols[mp_] = len(cols)
                rows.append(gi); cidx.append(cols[mp_]); vals_m.append(a * f); vals_r.append(b * f)
    Tm = np.zeros((JS.nc, len(cols)), complex); Tr = np.zeros((JS.nc, len(cols)))
    np.add.at(Tm, (rows, cidx), vals_m); np.add.at(Tr, (rows, cidx), vals_r)
    mex = np.array(list(cols.keys()), dtype=np.int64).reshape(-1, JS.nv)
    return Tm, Tr, mex

def monomial_values(vars4, mex):
    """affine values of prod v_i^{mex_i}, vars4 = list of 4 affine scalars (v3, v4 conj of v1, v2 or not)."""
    D = int(mex.max()) if len(mex) else 0
    pw = []
    for v in vars4:
        nn, n = v[0].shape
        p = [af_const(1.0, nn, n), v]
        for _ in range(2, D + 1):
            p.append(af_mul(p[-1], v))
        pw.append(p)
    def stackp(p, idx):
        return np.stack([p[i][0] for i in idx]), np.stack([p[i][1] for i in idx])
    X = af_mul(stackp(pw[0], mex[:, 0]), stackp(pw[1], mex[:, 1]))
    Y = af_mul(stackp(pw[2], mex[:, 2]), stackp(pw[3], mex[:, 3]))
    return af_mul(X, Y)

def poly_jet(T, mono):
    Tm, Tr, _ = T
    c, r = mono
    nm, nn, n = c.shape
    cj = (Tm @ c.reshape(nm, nn * n)).reshape(-1, nn, n)
    if nn == 1 and not np.any(r) and not np.any(Tr):
        return cj, np.zeros((Tm.shape[0], n))
    mag = af_mag(mono)
    aT = np.abs(Tm)
    rj = aT @ r + Tr @ mag + (aT @ mag) * (nm + 2) * U_RND
    return cj, rj

# ==================================================================================
# 8. pipelines
# ==================================================================================
JS4 = JetSpace(4); JS2 = JetSpace(2)
MAP_Z = mixed_maps(JS4, JS2)

def a8_of(a_lo, a_hi):
    lo, hi = a_lo ** 8, a_hi ** 8
    return 0.5 * (lo + hi), 0.5 * (hi - lo) * (1 + 1e-14) + 1e-30

def finish(G, E, a8, rigorous, detG_is_one):
    """G, E: dicts (i,j) -> order-2 jets of g0 and E entries.
    Q = [tr(adj(g0) E) + a8 det E]/det g0  (exact 2x2 identity: det(g0 + a8 E) = det g0 (1 + a8 Q)).
    detG_is_one: use the proved identity det g0 = 1 (zeta coordinates, DEFS 1.3) instead of computing it.
    Returns F/a^8 (affine), Ric/a^8 matrix (affine), g matrix, det g (affine), det g0 check."""
    M = JS2.mul
    trA = af_add(af_sub(M(G[1, 1], E[0, 0]), M(G[0, 1], E[1, 0])), af_sub(M(G[0, 0], E[1, 1]), M(G[1, 0], E[0, 1])))
    detE = af_sub(M(E[0, 0], E[1, 1]), M(E[0, 1], E[1, 0]))
    num = af_add(trA, af_iscale(detE, a8) if rigorous else af_scale(detE, a8[0]))
    detG = af_sub(M(G[0, 0], G[1, 1]), M(G[0, 1], G[1, 0])) if not rigorous else None
    if detG_is_one:
        Q = num
    else:
        x0 = JS2.const_part(detG); H = JS2.nonconst(detG)
        f0 = (1.0 / x0[0], x0[1] * 0); f1 = (-1.0 / x0[0] ** 2, x0[1] * 0); f2 = (1.0 / x0[0] ** 3, x0[1] * 0)
        Q = M(num, JS2.series([f0, f1, f2], H))
    Q0 = JS2.const_part(Q)
    HQ = JS2.nonconst(Q)
    if rigorous:
        am, ar = a8; ahi = am + ar
        qlo, qhi = af_real_range(Q0)
        qmax = np.maximum(np.abs(qlo), np.abs(qhi))
        x = ahi * qmax
        bad = ~(x < 0.1)
        x = np.where(bad, 0.0, x)
        Q0r = (Q0[0].real.astype(complex), Q0[1])
        sq = af_mul(Q0r, Q0r)
        l0 = af_sub(Q0r, af_iscale(sq, (0.5 * am, 0.5 * ar)))
        l0 = (l0[0], l0[1] + ahi ** 2 * qmax ** 3 / (3 * (1 - x) ** 3))
        l1 = af_sub(af_const(1.0, *Q0[0].shape), af_iscale(Q0r, a8))
        l1 = (l1[0], l1[1] + x ** 2 / (1 - x))
        l2 = af_const(-0.5 * am, *Q0[0].shape)
        l2 = (l2[0], l2[1] + 0.5 * ar + 0.5 * ahi * (1 / (1 - x) ** 2 - 1))
        dg = af_add(af_const(1.0, *Q0[0].shape), af_iscale(Q0r, a8))      # det g = det g0 (1 + a8 Q), det g0 = 1
        l0 = (l0[0], np.where(bad, np.inf, l0[1])); l1 = (l1[0], np.where(bad, np.inf, l1[1]))
        dg = (dg[0], np.where(bad, np.inf, dg[1]))
    else:
        a = a8[0]; q0 = Q0[0].real
        l0 = (np.log1p(a * q0) / a + 0j, Q0[1] * 0)
        l1 = (1.0 / (1 + a * q0) + 0j, Q0[1] * 0)
        l2 = (-0.5 * a / (1 + a * q0) ** 2 + 0j, Q0[1] * 0)
        dg = None
    L = JS2.series([l0, l1, l2], HQ)
    Ric = {}
    for i in range(2):
        for j in range(2):
            k = JS2.pos[tuple(1 if (v == i or v == 2 + j) else 0 for v in range(4))]
            Ric[i, j] = af_scale((L[0][k], L[1][k]), -1.0)
    g = {}
    for i in range(2):
        for j in range(2):
            e0 = JS2.const_part(E[i, j]); g0 = JS2.const_part(G[i, j])
            g[i, j] = af_add(g0, af_iscale(e0, a8) if rigorous else af_scale(e0, a8[0]))
    if dg is None:
        dg = af_sub(af_mul(g[0, 0], g[1, 1]), af_mul(g[0, 1], g[1, 0]))
    detG0 = None if detG is None else JS2.const_part(detG)[0][0].real
    return (l0[0], l0[1]), Ric, g, dg, detG0

def spectral(g, R, rigorous, dg):
    """max |eigenvalue| of g^-1 R (2x2): roots of dg x^2 - tau x + dR, tau = tr(adj(g) R)."""
    tau = af_sub(af_add(af_mul(g[0, 0], R[1, 1]), af_mul(g[1, 1], R[0, 0])), af_add(af_mul(g[0, 1], R[1, 0]), af_mul(g[1, 0], R[0, 1])))
    dR = af_sub(af_mul(R[0, 0], R[1, 1]), af_mul(R[0, 1], R[1, 0]))
    disc = af_sub(af_mul(tau, tau), af_scale(af_mul(dg, dR), 4.0))
    if not rigorous:
        t = tau[0][0].real; dd = disc[0][0].real; g_ = dg[0][0].real
        return (np.abs(t) + np.sqrt(np.maximum(dd, 0))) / (2 * g_)
    tlo, thi = af_real_range(tau); dlo, dhi = af_real_range(disc); glo, ghi = af_real_range(dg)
    num = _up(np.maximum(np.abs(tlo), np.abs(thi)) + _up(np.sqrt(np.maximum(dhi, 0.0))))
    with np.errstate(all="ignore"):
        out = _up(_up(num / _dn(2 * glo)) * (1 + 1e-14))
    return np.where((glo > 0) & np.isfinite(out), out, np.inf)

def zeta_pipeline(z1, z2, polys_T, a_lo, a_hi, rigorous, tails=None, want_detg0=False, sig_override=None, want_hint=False):
    """z1, z2 affine scalars (c (nn,n), r (n,)) in scaled zeta coords; returns (F/a^8, |Ric|/a^6) values or bounds.
    sig_override: affine sigma = |zeta|^2 supplied directly (polar boxes: sigma = rho^2 exactly)."""
    nn, n = z1[0].shape
    w1, w2 = af_conj(z1), af_conj(z2)
    sig = af_add(af_mul(z1, w1), af_mul(z2, w2)) if sig_override is None else sig_override
    sc = sig[0][0].real
    if rigorous:
        slo, shi = af_real_range(sig)
        slo = np.maximum(slo, 1e-6)
        D = (np.concatenate([np.zeros((1, n), complex), sig[0][1:].real.astype(complex)], 0), sig[1])
        cR, pR = radial_table().ranges(slo, shi)
        cC = c_jet_interval(sc, sc); pC = p1_jet_interval(sc, sc)
        ck = [tighter(af_from_interval_and_slope(cC[k], iscal(k + 1, cR[k + 1]), D), cR[k]) for k in range(5)]
        ph = [None] + [tighter(af_from_interval_and_slope(iscal(1.0 / k, pC[k - 1]), pR[k], D), iscal(1.0 / k, pR[k - 1])) for k in range(1, 5)]
    else:
        cC = c_jet_interval(sc, sc); pC = p1_jet_interval(sc, sc)
        mid = lambda x: (0.5 * (x[0] + x[1]))
        ck = [(np.array([mid(cC[k]) + 0j]), np.zeros(n)) for k in range(5)]
        ph = [None] + [(np.array([mid(pC[k - 1]) / k + 0j]), np.zeros(n)) for k in range(1, 5)]
    # Delta = sigma(zeta0 + dz) - sigma0 as a jet
    Dl = JS4.zeros(nn, n)
    for i, (zi, wi) in enumerate(((z1, w1), (z2, w2))):
        Dl[0][JS4.unit(i)] = wi[0]; Dl[1][JS4.unit(i)] = wi[1]
        Dl[0][JS4.unit(2 + i)] = zi[0]; Dl[1][JS4.unit(2 + i)] = zi[1]
        t = [0] * 4; t[i] = 1; t[2 + i] = 1
        Dl[0][JS4.pos[tuple(t)]][0] = 1.0
    Pw = [None, Dl]
    for k in range(2, 5):
        Pw.append(JS4.mul(Pw[-1], Dl))
    PhiEH = JS4.series_pw(ph, Pw)
    Cj = JS4.series_pw(ck, Pw)
    mono = monomial_values([z1, z2, w1, w2], polys_T["mex"])
    Aj = poly_jet(polys_T["A"], mono)
    Bj = poly_jet(polys_T["B"], mono)
    Ph = af_add(JS4.mul(Cj, Aj), Bj)
    if tails is not None:
        tau_oo, tau_ee = tails
        magC = af_mag(Cj)
        add = np.add.reduceat(magC[JS4.I] * tau_oo[JS4.J][:, None], JS4.starts, axis=0)
        Ph = (Ph[0], Ph[1] + add * (1 + 1e-12) + tau_ee[:, None])
    G = {k: extract(PhiEH, MAP_Z[k]) for k in MAP_Z}
    E = {k: extract(Ph, MAP_Z[k]) for k in MAP_Z}
    a8 = a8_of(a_lo, a_hi)
    Fa8, Ric, g, dg, detG0 = finish(G, E, a8, rigorous, detG_is_one=rigorous)
    rho = spectral(g, Ric, rigorous, dg)
    if rigorous:
        flo, fhi = af_real_range(Fa8)
        UF = _up(np.maximum(np.abs(flo), np.abs(fhi)))
        if want_hint:
            sens = sum(np.abs(Ric[k][0][1:]) for k in Ric).real + 0 * rho
            sensF = np.abs(Fa8[0][1:])
            return UF, rho, sens, sensF
        return UF, rho
    if want_detg0:
        return Fa8[0][0].real, rho, detG0
    return Fa8[0][0].real, rho

def prepare_T(polys, mode_rig):
    """shift matrices for A and B over a common monomial basis."""
    TA = shift_matrix(JS4, *polys["A"]); TB = shift_matrix(JS4, *polys["B"])
    allm = {}
    for mex in (TA[2], TB[2]):
        for m in map(tuple, mex):
            if m not in allm:
                allm[m] = len(allm)
    mex = np.array(list(allm.keys()), dtype=np.int64)
    def remap(T):
        Tm, Tr, mx = T
        Tm2 = np.zeros((Tm.shape[0], len(allm)), complex); Tr2 = np.zeros((Tm.shape[0], len(allm)))
        cols = [allm[tuple(m)] for m in mx]
        Tm2[:, cols] = Tm; Tr2[:, cols] = Tr
        return Tm2, Tr2, mex
    return {"A": remap(TA), "B": remap(TB), "mex": mex}

def point_af(x):
    """x: (n,4) real points -> affine scalars zeta1, zeta2 in point mode."""
    z1 = (x[:, 0] + 1j * x[:, 1])[None, :]; z2 = (x[:, 2] + 1j * x[:, 3])[None, :]
    n = x.shape[0]
    return (z1, np.zeros(n)), (z2, np.zeros(n))

def box_af(m, h):
    """box centre m (n,4), half widths h (n,4) -> affine zeta1, zeta2 with 4 noise symbols."""
    n = m.shape[0]
    z1 = np.zeros((5, n), complex); z2 = np.zeros((5, n), complex)
    z1[0] = m[:, 0] + 1j * m[:, 1]; z1[1] = h[:, 0]; z1[2] = 1j * h[:, 1]
    z2[0] = m[:, 2] + 1j * m[:, 3]; z2[3] = h[:, 2]; z2[4] = 1j * h[:, 3]
    return (z1, np.zeros(n)), (z2, np.zeros(n))

# ---------------- resolved chart (u, w), point mode only -----------------------
JS4c = JS4; MAP_C = MAP_Z

def jet_func(JS, J, kind):
    c0 = J[0][0][0]
    H = JS.nonconst(J); n = c0.shape[0]
    fk = []
    for k in range(5):
        if kind == "sqrt":
            v = binom_half(k) * c0 ** (0.5 - k)
        elif kind == "log":
            v = np.log(c0) if k == 0 else ((-1) ** (k + 1)) / (k * c0 ** k)
        fk.append((v[None, :], np.zeros(n)))
    return JS.series(fk, H)

def chart_pipeline(u, w, polyB_chart_T, a, with_E=True):
    """point mode in chart variables (u, w, ubar, wbar). Returns (F/a^8, |Ric|/a^6, F_EH-check, RicEH-check)."""
    n = u.shape[0]
    def var(i, val):
        J = JS4.zeros(1, n); J[0][0][0] = val; J[0][JS4.unit(i)][0] = 1.0
        return J
    U = var(0, u); W = var(1, w); Ub = var(2, np.conj(u)); Wb = var(3, np.conj(w))
    one = JS4.zeros(1, n); one[0][0][0] = 1.0
    uu = JS4.mul(U, Ub)
    ww1 = af_add(one, JS4.mul(W, Wb))
    sc = JS4.mul(uu, JS4.mul(ww1, ww1))
    fh = jet_func(JS4, af_add(one, sc), "sqrt")
    Phi = af_sub(af_add(fh, jet_func(JS4, ww1, "log")), jet_func(JS4, af_add(one, fh), "log"))
    if with_E:
        mono = monomial_values([(u[None, :] + 0j, np.zeros(n)), (w[None, :] + 0j, np.zeros(n)),
                                (np.conj(u)[None, :], np.zeros(n)), (np.conj(w)[None, :], np.zeros(n))], polyB_chart_T[2])
        Ph = poly_jet(polyB_chart_T, mono)
    else:
        Ph = JS4.zeros(1, n)
    G = {k: extract(Phi, MAP_C[k]) for k in MAP_C}
    E = {k: extract(Ph, MAP_C[k]) for k in MAP_C}
    a8 = (a ** 8, 0.0)
    Fa8, Ric, g, dg, detG0 = finish(G, E, a8, False, detG_is_one=False)
    rho = spectral(g, Ric, False, dg)
    # Q4-type control: log det g0 must be pluriharmonic in the chart
    detG = af_sub(JS2.mul(G[0, 0], G[1, 1]), JS2.mul(G[0, 1], G[1, 0]))
    ld = jet_func2(detG)
    ric_eh = max(np.abs(ld[0][JS2.pos[(1, 0, 1, 0)]][0]).max(), np.abs(ld[0][JS2.pos[(0, 1, 0, 1)]][0]).max(), np.abs(ld[0][JS2.pos[(1, 0, 0, 1)]][0]).max())
    return Fa8[0][0].real, rho, ric_eh

def jet_func2(J):
    c0 = J[0][0][0]; H = JS2.nonconst(J); n = c0.shape[0]
    fk = [(np.log(c0)[None, :], np.zeros(n)), ((1 / c0)[None, :], np.zeros(n)), ((-0.5 / c0 ** 2)[None, :], np.zeros(n))]
    return JS2.series(fk, H)

def chart_poly_T(polys):
    """B (even-even) monomials z^a w^b -> u^{(a1+a2)/2} w^{a2} ubar^{(b1+b2)/2} wbar^{b2}."""
    ex, cm, cr = polys["B"]
    assert np.all((ex[:, 0] + ex[:, 1]) % 2 == 0) and np.all((ex[:, 2] + ex[:, 3]) % 2 == 0)
    ex2 = np.stack([(ex[:, 0] + ex[:, 1]) // 2, ex[:, 1], (ex[:, 2] + ex[:, 3]) // 2, ex[:, 3]], 1)
    return shift_matrix(JS4, ex2, cm, cr)

# ==================================================================================
# 9. independent checks in mpmath (different formulation: direct log det of a finite-
#    difference Hessian of the full potential; Ricci by finite differences of F)
# ==================================================================================
def mp_direct(xpt, a, polys, win_lo=None, dps=100, h_in="1e-30", h_out="1e-8"):
    import mpmath as mp
    mp.mp.dps = dps
    wl = mp.mpf(WIN_LO if win_lo is None else win_lo); wh = mp.mpf(WIN_HI)
    a = mp.mpf(a); a8 = a ** 8
    terms = []
    for key in ("A", "B"):
        ex, cm, _ = polys[key]
        for e, c in zip(ex, cm):
            terms.append((key, tuple(int(t) for t in e), mp.mpc(c.real, c.imag)))
    D = max([max(t[1]) for t in terms] + [1])
    def h(t):
        return mp.e ** (-1 / t) if t > 0 else mp.mpf(0)
    def st(t):
        return h(t) / (h(t) + h(1 - t))
    def Phi(x):
        z1 = mp.mpc(x[0], x[1]); z2 = mp.mpc(x[2], x[3])
        s = x[0] ** 2 + x[1] ** 2 + x[2] ** 2 + x[3] ** 2
        fh = mp.sqrt(s * s + 1)
        phi1 = fh + mp.log(s) - mp.log(1 + fh)
        c = st((mp.sqrt(s) - wl) / (wh - wl))
        pw = []
        for v in (z1, z2, mp.conj(z1), mp.conj(z2)):
            p = [mp.mpf(1)]
            for _ in range(D):
                p.append(p[-1] * v)
            pw.append(p)
        A = mp.mpc(0); B = mp.mpc(0)
        for key, e, cf in terms:
            m = cf * pw[0][e[0]] * pw[1][e[1]] * pw[2][e[2]] * pw[3][e[3]]
            if key == "A":
                A += m
            else:
                B += m
        return phi1 + a8 * mp.re(c * A + B)
    def hess(f, x, hh):
        H = [[None] * 4 for _ in range(4)]
        f0 = f(x)
        for i in range(4):
            for j in range(i, 4):
                def sh(di, dj):
                    y = list(x); y[i] += di * hh; y[j] += dj * hh
                    return f(y)
                if i == j:
                    H[i][i] = (sh(1, 1) - 2 * f0 + sh(-1, -1)) / (4 * hh * hh)
                else:
                    H[i][j] = H[j][i] = (sh(1, 1) - sh(1, -1) - sh(-1, 1) + sh(-1, -1)) / (4 * hh * hh)
        return H
    def wirt(H):
        g = [[None] * 2 for _ in range(2)]
        for i in range(2):
            for j in range(2):
                xi, yi, xj, yj = 2 * i, 2 * i + 1, 2 * j, 2 * j + 1
                g[i][j] = (H[xi][xj] + H[yi][yj] + 1j * (H[xi][yj] - H[yi][xj])) / 4
        return g
    hin = mp.mpf(h_in); hout = mp.mpf(h_out)
    def Ffun(x):
        g = wirt(hess(Phi, x, hin))
        return mp.log(mp.re(g[0][0] * g[1][1] - g[0][1] * g[1][0]))
    x = [mp.mpf(float(t)) for t in xpt]
    g = wirt(hess(Phi, x, hin))
    F = mp.log(mp.re(g[0][0] * g[1][1] - g[0][1] * g[1][0]))
    Rm = wirt(hess(Ffun, x, hout))
    Rm = [[-Rm[i][j] for j in range(2)] for i in range(2)]
    dg = mp.re(g[0][0] * g[1][1] - g[0][1] * g[1][0])
    tau = mp.re(g[0][0] * Rm[1][1] + g[1][1] * Rm[0][0] - g[0][1] * Rm[1][0] - g[1][0] * Rm[0][1])
    dR = mp.re(Rm[0][0] * Rm[1][1] - Rm[0][1] * Rm[1][0])
    rho = (abs(tau) + mp.sqrt(max(tau * tau - 4 * dg * dR, 0))) / (2 * dg)
    return float(F / a8), float(rho / a8)

def mp_psi_direct(x, a, dps=50, Rreal=3.6, Kmax=9.0, Rra=5.0):
    """psi(x) - psi(0) from the definitions: (mu-1)|x|^2 + sum_{l!=0} R_a(|x-l|^2) - (eps/2) S_reg(x),
    S by Ewald (eta = 6), R_a(s) = phi_a(s) - s + eps/(2 s) (closed form), mu - 1 = 8 pi^2 eps - delta."""
    import mpmath as mp
    mp.mp.dps = dps
    a = mp.mpf(a); eps = a ** 4; eta = mp.mpf(6)
    x = [mp.mpf(float(t)) for t in x]
    Vm = mp.mpf(1) / 32
    def Ra(s):
        f = mp.sqrt(s * s + eps)
        return f + a * a * mp.log(s) - a * a * mp.log(a * a + f) - s + eps / (2 * s)
    def Sreg(y):
        tot = mp.mpf(0)
        for l in Lreal:
            d2 = sum((y[i] - l[i]) ** 2 for i in range(4))
            tot += mp.e ** (-eta * d2) / d2
        s0 = sum(t * t for t in y)
        tot += (mp.e ** (-eta * s0) - 1) / s0 if s0 > 0 else -eta
        for k in Krec:
            k2 = sum(t * t for t in k)
            tot += mp.e ** (-mp.pi ** 2 * k2 / eta) / k2 * mp.cos(2 * mp.pi * sum(k[i] * y[i] for i in range(4))) / Vm
        return tot - mp.pi ** 2 / (Vm * eta)
    def Rsum(y):
        tot = mp.mpf(0)
        for l in Lra:
            tot += Ra(sum((y[i] - l[i]) ** 2 for i in range(4)))
        return tot
    Lreal = [[mp.mpf(float(t)) for t in l] for l in lattice_points(Rreal)]
    Lra = [[mp.mpf(float(t)) for t in l] for l in lattice_points(Rra)]
    M = int(Kmax / 2) + 1
    Krec = []
    for k in itertools.product(range(-M, M + 1), repeat=4):
        kk = [2 * t for t in k]
        if any(kk) and (sum(kk) // 2) % 2 == 0 and sum(t * t for t in kk) <= Kmax ** 2:
            Krec.append([mp.mpf(t) for t in kk])
    delta = mp.mpf(0)
    for k in range(2, 12):
        delta += (1 - k) * mp.binomial(mp.mpf(1) / 2, k) * eps ** k * Zkappa(k)
    mu1 = 8 * mp.pi ** 2 * eps - delta
    zero = [mp.mpf(0)] * 4
    s = sum(t * t for t in x)
    val = mu1 * s + (Rsum(x) - Rsum(zero)) - eps / 2 * (Sreg(x) - Sreg(zero))
    return val

def poly_psi(table, x, a, Dmax=DMAX):
    """sum_{k<=KMAX, 4<=d<=Dmax} eps^k P_{k,d}(x)  (z-units, real x)."""
    eps = a ** 4
    z1 = x[0] + 1j * x[1]; z2 = x[2] + 1j * x[3]
    v = np.array([z1, z2, np.conj(z1), np.conj(z2)])
    tot = 0j
    for k in range(1, KMAX + 1):
        for d in range(4, Dmax + 1, 2):
            if k == 1 and d == 4:
                continue
            ex, cf = load_poly(table, k, d)
            tot += eps ** k * np.sum(cf * np.prod(v[None, :] ** ex, axis=1))
    return tot.real

# ==================================================================================
# 10. stages
# ==================================================================================
def make_sampler(table, a, Dexp=16):
    polys = build_polys(table, a, a, Dexp, False)
    T = prepare_T(polys, False)
    TB = chart_poly_T(polys)
    def ev(x, batch=2000):
        x = np.atleast_2d(np.asarray(x, dtype=float))
        n = len(x)
        F = np.empty(n); Rr = np.empty(n)
        r = np.linalg.norm(x, axis=1)
        mz = r >= 0.2
        for idx in (np.nonzero(mz)[0],):
            for b in range(0, len(idx), batch):
                ii = idx[b:b + batch]
                F[ii], Rr[ii] = zeta_pipeline(*point_af(x[ii]), T, a, a, False)
        idx = np.nonzero(~mz)[0]
        if len(idx):
            xc = x[idx].copy()
            z1 = xc[:, 0] + 1j * xc[:, 1]; z2 = xc[:, 2] + 1j * xc[:, 3]
            sw = np.abs(z2) > np.abs(z1)
            z1, z2 = np.where(sw, z2, z1), np.where(sw, z1, z2)
            u = z1 ** 2
            w = np.where(np.abs(z1) > 0, z2 / np.where(np.abs(z1) > 0, z1, 1), 0)
            for b in range(0, len(idx), batch):
                sl = slice(b, b + batch)
                f_, r_, _ = chart_pipeline(u[sl], w[sl], TB, a)
                F[idx[sl]] = f_; Rr[idx[sl]] = r_
        return F, Rr
    def ev_chart(u, w):
        f_, r_, ric_eh = chart_pipeline(u, w, TB, a)
        return f_, r_
    return ev, ev_chart, polys, T, TB

def sample_points(seed=20261113):
    rng = np.random.default_rng(seed)
    d = rng.normal(size=(10000, 4)); d /= np.linalg.norm(d, axis=1)[:, None]
    S1 = d * rng.uniform(0, 1, size=(10000, 1))
    d = rng.normal(size=(8000, 4)); d /= np.linalg.norm(d, axis=1)[:, None]
    S2 = d * rng.uniform(0, 1, size=(8000, 1)) ** 0.25
    dirs = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [1, 1, 0, 0], [1, 0, 1, 0], [1, 0, 0, 1],
                     [1, 1, 1, 1], [1, 1, 1, -1], [1, 1, 1, 0], [2, 1, 0, 0], [1, 0, 1, 1]], dtype=float)
    dirs /= np.linalg.norm(dirs, axis=1)[:, None]
    rr = np.linspace(0.004, 1.0, 250)
    S3 = (dirs[:, None, :] * rr[None, :, None]).reshape(-1, 4)
    X = np.concatenate([S1, S2, S3])
    # exceptional curve and its neighbourhood, chart points (u, w), |w| <= 1
    m = 600
    w = np.sqrt(rng.uniform(0, 1, m)) * np.exp(2j * np.pi * rng.uniform(0, 1, m))
    u = np.where(np.arange(m) < 300, 0.0, 1e-4 * rng.uniform(0, 1, m) * np.exp(2j * np.pi * rng.uniform(0, 1, m)))
    return X, u.astype(complex), w

def pattern_refine(f, seeds, s0=0.02, smin=2e-6, rng=None, logp=None):
    """synchronous pattern search maximising f over |x| <= 1. Returns history of (step, max, argmax)."""
    rng = rng or np.random.default_rng(7)
    X = seeds.copy(); v = f(X)
    hist = []
    s = s0
    while s >= smin:
        for _ in range(200):
            dirs = np.concatenate([np.eye(4), -np.eye(4)])
            rd = rng.normal(size=(8, 4)); rd /= np.linalg.norm(rd, axis=1)[:, None]
            dirs = np.concatenate([dirs, rd, -rd])
            C = X[:, None, :] + s * dirs[None, :, :]
            nr = np.linalg.norm(C, axis=2)
            C = np.where(nr[..., None] > 1, C / nr[..., None], C)
            vc = f(C.reshape(-1, 4)).reshape(len(X), -1)
            j = np.argmax(vc, axis=1)
            best = vc[np.arange(len(X)), j]
            imp = best > v * (1 + 1e-15)
            if not np.any(imp):
                break
            X[imp] = C[np.arange(len(X)), j][imp]; v[imp] = best[imp]
        k = int(np.argmax(v))
        hist.append((s, float(v[k]), X[k].tolist()))
        if logp:
            logp("    step %.1e  max %.10e  at |x|=%.6f" % (s, v[k], np.linalg.norm(X[k])))
        s *= 0.5
    return hist, X, v

def pick_seeds(X, vals, k=20, sep=0.03):
    order = np.argsort(-vals)
    chosen = []
    for i in order:
        if all(np.linalg.norm(X[i] - X[j]) > sep for j in chosen):
            chosen.append(i)
        if len(chosen) == k:
            break
    return np.array(chosen)

def converged_value(hist, tol=1e-4):
    """first step level at which successive maxima agree to tol (relative)."""
    for (s1, m1, _), (s2, m2, _) in zip(hist, hist[1:]):
        if abs(m2 - m1) <= tol * abs(m2):
            return s2, m2
    return None, None

def stage_sample(alist=(0.03, 0.01)):
    table = stage_coeffs()
    X, uE, wE = sample_points()
    out = {"n_points_zeta": int(len(X)), "n_points_chart_E": int(len(uE))}
    for a in alist:
        fa = os.path.join(CACHE, "sample_a%g.json" % a)
        if os.path.exists(fa):
            out["a=%g" % a] = json.load(open(fa))
            continue
        t0 = time.time()
        ev, evc, polys, T, TB = make_sampler(table, a)
        F, Rr = ev(X)
        Fe, Re = evc(uE, wE)
        log("a=%g: sampled %d + %d points in %.1fs" % (a, len(X), len(uE), time.time() - t0))
        res = {}
        r = np.linalg.norm(X, axis=1)
        win = (r >= WIN_LO) & (r <= WIN_HI)
        res["sampled_max_absF_over_a8"] = float(max(np.abs(F).max(), np.abs(Fe).max()))
        res["sampled_max_Ric_over_a6"] = float(max(Rr.max(), Re.max()))
        res["sampled_max_absF_window"] = float(np.abs(F[win]).max())
        res["sampled_max_Ric_window"] = float(Rr[win].max())
        res["inner_|zeta|<1/4_max_absF"] = float(max(np.abs(F[r < WIN_LO]).max(), np.abs(Fe).max()))
        res["inner_|zeta|<1/4_max_Ric"] = float(max(Rr[r < WIN_LO].max(), Re.max()))
        res["on_E_F_over_a8_range"] = [float(Fe[:300].min()), float(Fe[:300].max())]
        res["on_E_Ric_over_a6_max"] = float(Re[:300].max())
        for name, vals, sgn in (("absF", np.abs(F), 1), ("Ric", Rr, 1)):
            seeds = X[pick_seeds(X, vals)]
            log("  refine %s from 20 seeds (a=%g), seed max %.8e" % (name, a, vals.max()))
            f = (lambda Y: np.abs(ev(Y)[0])) if name == "absF" else (lambda Y: ev(Y)[1])
            hist, Xf, vf = pattern_refine(f, seeds, logp=log)
            s_conv, m_conv = converged_value(hist)
            k = int(np.argmax(vf)); xb = Xf[k]
            Fb, Rb = ev(xb[None])
            res["refined_" + name] = {"max": float(vf[k]), "argmax_x": xb.tolist(), "argmax_|zeta|": float(np.linalg.norm(xb)),
                                      "F_over_a8_at_argmax": float(Fb[0]), "Ric_over_a6_at_argmax": float(Rb[0]),
                                      "history": [(h[0], h[1]) for h in hist], "first_step_with_1e-4_agreement": s_conv,
                                      "seed_values_top20": sorted(vals[pick_seeds(X, vals)].tolist(), reverse=True)}
            # truncation-degree convergence at the argmax
            dv = {}
            for Dx in (12, 14, 16, 18):
                evx = make_sampler(table, a, Dx)[0]
                dv[str(Dx)] = [float(t) for t in (lambda r_: (r_[0][0], r_[1][0]))(evx(xb[None]))]
            res["refined_" + name]["degree_convergence_F_Ric"] = dv
            log("   -> %s max %.10e at x=%s (|zeta|=%.6f), D-conv %s" % (name, vf[k], np.round(xb, 6), np.linalg.norm(xb), dv))
        out["a=%g" % a] = res
        json.dump(res, open(fa, "w"), indent=1)
        log("  a=%g done %.1fs" % (a, time.time() - t0))
    return out


# ==================================================================================
# 11. rigorous stage: window zone 1/4 <= |zeta| <= 1, a in [0.029, 0.03]
# ==================================================================================
A_LO, A_HI = 0.029, 0.03
TOL_R, TOL_F = 0.03, 0.05

def tail_bounds(a_hi, rho_m, Dexp=None, dmax=80, kmax=12):
    """bounds tau[gamma] on jet coefficients (scaled zeta units) of everything not in the explicit polynomial:
    lattice truncation (A2 + A4), degrees > Dexp (A2, L3(i)), eps^k with k > KMAX (L3(i)); Cauchy estimate A3."""
    Dexp = DRIG if Dexp is None else Dexp
    orders = np.array([sum(g) for g in JS4.idx])
    def fac(n, d):
        if n == 0:
            return rho_m ** d
        if d > n:
            dl = n * rho_m / (2.0 * (d - n))
            return dl ** (-n) * (rho_m + 2 * dl) ** d
        return 2.0 ** d if d == n else 0.0
    def Cgen(k, d):
        if k == 1:
            return 0.5 * (d + 1) * float(Zkappa((d + 2) / 4.0))
        return abs(binom_half(k)) / (2 * k - 1) * math.comb(d + 4 * k - 3, d) * float(Zkappa((4 * k - 2 + d) / 4.0))
    pieces = []
    for k in range(1, KMAX + 1):
        for d in range(4, Dexp + 1, 2):
            if (k, d) in ((1, 4), (1, 6), (2, 4)):
                continue
            nu = 1 if k == 1 else 2 * k - 1
            pref = 0.5 if k == 1 else abs(binom_half(k) / (1 - 2 * k))
            pieces.append((k, d, pref * math.comb(d + 2 * nu - 1, d) * A4_tail(radius_for(nu, d), 2 * nu + d)))
    for k in range(1, kmax + 1):
        for d in range(4, dmax + 1, 2):
            if d > Dexp or k > KMAX:
                pieces.append((k, d, Cgen(k, d)))
    tau_oo = np.zeros(JS4.nc); tau_ee = np.zeros(JS4.nc)
    last = 0.0
    for k, d, C in pieces:
        e = 4 * k + d - 10
        for n in range(5):
            v = C * a_hi ** e * fac(n, d)
            if d % 4 == 2:
                tau_oo[orders == n] += v
            else:
                tau_ee[orders == n] += v
            if d == dmax or k == kmax:
                last = max(last, v)
    assert last < 1e-25, last
    pad = 1 + 1e-6      # rounding in this bookkeeping and the omitted d > dmax, k > kmax (terms decay by > 1e3 per step)
    return tau_oo * pad + 1e-25, tau_ee * pad + 1e-25

def region_keep(m, h):
    lo = m - h; hi = m + h
    mn = np.maximum(0, np.maximum(lo, -hi)); mx = np.maximum(np.abs(lo), np.abs(hi))
    ok = np.all(hi >= 0, axis=1)
    ok &= hi[:, 0] >= lo[:, 1]
    ok &= (mx[:, 0] ** 2 + mx[:, 1] ** 2) >= (mn[:, 2] ** 2 + mn[:, 3] ** 2)
    ok &= (mn ** 2).sum(1) <= WIN_HI ** 2
    ok &= (mx ** 2).sum(1) >= WIN_LO ** 2
    return ok

def rig_eval(m, h, TR, a_lo=A_LO, a_hi=A_HI):
    """Cartesian boxes (used only in probes; the production cover uses polar boxes)."""
    rho_m = float((np.linalg.norm(np.abs(m) + h, axis=1)).max())
    tails = tail_bounds(a_hi, rho_m)
    z1, z2 = box_af(m, h)
    return zeta_pipeline(z1, z2, TR, a_lo, a_hi, True, tails=tails)

def polar_af(face, rc, hr, uc, hu):
    """polar box: rho in rc +- hr, direction n = y/|y| with y_face = 1 and the other three y_j in uc_j +- hu_j (j in order).
    Noise symbols: e1 = rho, e2..e4 = u1..u3. Returns zeta1, zeta2 affine and sigma = rho^2 affine."""
    n = len(rc)
    others = [[j for j in range(4) if j != f] for f in range(4)]
    rho = np.zeros((5, n), complex); rho[0] = rc; rho[1] = hr
    rho = (rho, np.zeros(n))
    yaf = []
    for comp in range(4):
        c = np.zeros((5, n), complex); c[0] = np.where(face == comp, 1.0, 0.0)
        for f in range(4):
            sel = face == f
            if comp != f:
                jpos = others[f].index(comp)
                c[0, sel] = uc[sel, jpos]
                c[2 + jpos, sel] = hu[sel, jpos]
        yaf.append((c, np.zeros(n)))
    ny2 = af_mul(yaf[0], yaf[0])
    for comp in range(1, 4):
        ny2 = af_add(ny2, af_mul(yaf[comp], yaf[comp]))
    ny2 = (ny2[0].real.astype(complex), ny2[1])
    inv = af_func_real(ny2, lambda t: 1.0 / np.sqrt(t), lambda lo, hi: (_dn(-0.5 * np.maximum(lo, 1e-300) ** -1.5), _up(-0.5 * hi ** -1.5)))
    x = [af_mul(rho, af_mul(yaf[comp], inv)) for comp in range(4)]
    z1 = (x[0][0] + 1j * x[1][0], x[0][1] + x[1][1])
    z2 = (x[2][0] + 1j * x[3][0], x[2][1] + x[3][1])
    sig = af_mul(rho, rho)
    sig = (sig[0].real.astype(complex), sig[1])
    return z1, z2, sig

def polar_dirs(face, uc):
    n = len(face); y = np.zeros((n, 4))
    others = [[j for j in range(4) if j != f] for f in range(4)]
    for f in range(4):
        sel = face == f
        y[sel, f] = 1.0
        for jpos, j in enumerate(others[f]):
            y[sel, j] = uc[sel, jpos]
    return y / np.linalg.norm(y, axis=1)[:, None]

def polar_keep(face, rc, hr, uc, hu):
    """conservative test: box meets {x0 >= x1, x0^2 + x1^2 >= x2^2 + x3^2} (with the positive orthant of directions
    this is a fundamental domain of the order-64 subgroup of G_J) and 1/4 <= rho <= 1."""
    n = len(face)
    others = [[j for j in range(4) if j != f] for f in range(4)]
    ylo = np.zeros((n, 4)); yhi = np.zeros((n, 4))
    for f in range(4):
        sel = face == f
        ylo[sel, f] = 1.0; yhi[sel, f] = 1.0
        for jpos, j in enumerate(others[f]):
            ylo[sel, j] = np.maximum(uc[sel, jpos] - hu[sel, jpos], 0.0); yhi[sel, j] = uc[sel, jpos] + hu[sel, jpos]
    nmin = ylo / np.linalg.norm(yhi, axis=1)[:, None]; nmax = yhi / np.linalg.norm(ylo, axis=1)[:, None]
    ok = nmax[:, 0] >= nmin[:, 1]
    ok &= (nmax[:, 0] ** 2 + nmax[:, 1] ** 2) >= (nmin[:, 2] ** 2 + nmin[:, 3] ** 2)
    ok &= (rc - hr <= WIN_HI) & (rc + hr >= WIN_LO)
    return ok

def rig_eval_polar(face, rc, hr, uc, hu, TR, a_lo=A_LO, a_hi=A_HI):
    rho_m = float((rc + hr).max())
    tails = tail_bounds(a_hi, rho_m)
    z1, z2, sig = polar_af(face, rc, hr, uc, hu)
    return zeta_pipeline(z1, z2, TR, a_lo, a_hi, True, tails=tails, sig_override=sig, want_hint=True)


def stage_interval(time_limit=530.0, NB=1500):
    """adaptive cover of the window zone by polar boxes (rho x direction cell on a face of the cube), modulo the
    order-64 subgroup of G_J {diag(i^a, i^b), swap, conj}; rigorous affine evaluation; resumable (cache/)."""
    fs = os.path.join(CACHE, "interval_state.npz")
    fdone = os.path.join(CACHE, "interval_result.json")
    if os.path.exists(fdone):
        return json.load(open(fdone))
    t0 = time.time()
    table = stage_coeffs()
    polys = build_polys(table, A_LO, A_HI, DRIG, True)
    TR = prepare_T(polys, True)
    ev03 = make_sampler(table, A_HI)[0]
    if os.path.exists(fs):
        S = dict(np.load(fs, allow_pickle=True))
        Q = {k: S[k] for k in ("face", "rc", "hr", "uc", "hu")}
        st = json.loads(str(S["stats"]))
    else:
        s03 = json.load(open(os.path.join(CACHE, "sample_a0.03.json")))
        ev29 = make_sampler(table, A_LO)[0]
        orb = json.load(open(os.path.join(CACHE, "orbits.json")))["a=0.03"]
        xR = np.array(orb["Ric"]["argmax_x"]); xF = np.array(orb["absF"]["argmax_x"])
        MR = max(ev03(xR[None])[1][0], ev29(xR[None])[1][0], s03["sampled_max_Ric_window"], orb["Ric"]["max"])
        MF = max(abs(ev03(xF[None])[0][0]), abs(ev29(xF[None])[0][0]), s03["sampled_max_absF_window"], orb["absF"]["max"])
        nr, nu = 12, 8
        rcs = WIN_LO + (WIN_HI - WIN_LO) * (np.arange(nr) + 0.5) / nr
        ucs = (np.arange(nu) + 0.5) / nu
        F_, R_, U1, U2, U3 = np.meshgrid(np.arange(4), rcs, ucs, ucs, ucs, indexing="ij")
        Q = {"face": F_.ravel().astype(int), "rc": R_.ravel(), "hr": np.full(F_.size, 0.5 * (WIN_HI - WIN_LO) / nr),
             "uc": np.stack([U1.ravel(), U2.ravel(), U3.ravel()], 1), "hu": np.full((F_.size, 3), 0.5 / nu)}
        k = polar_keep(Q["face"], Q["rc"], Q["hr"], Q["uc"], Q["hu"])
        Q = {kk: v[k] for kk, v in Q.items()}
        st = {"thrR": (1 + TOL_R) * MR, "thrF": (1 + TOL_F) * MF, "MR_seed": MR, "MF_seed": MF,
              "n_eval": 0, "n_acc": 0, "maxUB_R": 0.0, "maxUB_F": 0.0, "boxUB_R": None, "boxUB_F": None,
              "maxLB_R": 0.0, "maxLB_F": 0.0, "n_unresolved": 0, "pvol_acc": 0.0, "pvol_init": float((2 * Q["hr"] * np.prod(2 * Q["hu"], 1)).sum()),
              "runs": 0, "elapsed": 0.0, "n_init": int(len(Q["rc"]))}
        log("thresholds: R %.6e (seed %.6e), F %.6e (seed %.6e); initial polar boxes %d" % (st["thrR"], MR, st["thrF"], MF, len(Q["rc"])))
    st["runs"] += 1
    it = 0
    def take(Q, n):
        a = {k: v[-n:] for k, v in Q.items()}; b = {k: v[:-n] for k, v in Q.items()}
        return a, b
    def cat(Q, B):
        return {k: np.concatenate([Q[k], B[k]]) for k in Q}
    while len(Q["rc"]) and time.time() - t0 < time_limit:
        B, Q = take(Q, min(NB, len(Q["rc"])))
        UF, UR, sens, sensF = rig_eval_polar(B["face"], B["rc"], B["hr"], B["uc"], B["hu"], TR)
        xc = polar_dirs(B["face"], B["uc"]) * B["rc"][:, None]
        Fc, Rc = ev03(xc)
        st["n_eval"] += len(UR)
        st["maxLB_R"] = max(st["maxLB_R"], float(Rc.max())); st["maxLB_F"] = max(st["maxLB_F"], float(np.abs(Fc).max()))
        ok = (UR <= st["thrR"]) & (UF <= st["thrF"])
        tiny = (B["hr"] < 2e-5) & (B["hu"].max(1) < 1e-4)
        acc = ok | tiny
        st["n_unresolved"] += int((tiny & ~ok).sum())
        if np.any(acc):
            st["n_acc"] += int(acc.sum())
            st["pvol_acc"] += float((2 * B["hr"][acc] * np.prod(2 * B["hu"][acc], 1)).sum())
            iR = int(np.argmax(np.where(acc, UR, -1))); iF = int(np.argmax(np.where(acc, UF, -1)))
            if UR[iR] > st["maxUB_R"]:
                st["maxUB_R"] = float(UR[iR]); st["boxUB_R"] = [int(B["face"][iR]), float(B["rc"][iR]), float(B["hr"][iR]), B["uc"][iR].tolist(), B["hu"][iR].tolist()]
            if UF[iF] > st["maxUB_F"]:
                st["maxUB_F"] = float(UF[iF]); st["boxUB_F"] = [int(B["face"][iF]), float(B["rc"][iF]), float(B["hr"][iF]), B["uc"][iF].tolist(), B["hu"][iF].tolist()]
        rej = ~acc
        if np.any(rej):
            w = np.where((UR > st["thrR"])[None, :], sens / st["thrR"], 0) + np.where((UF > st["thrF"])[None, :], sensF.real / st["thrF"], 0)
            w = w / np.maximum(w.sum(0), 1e-300)
            yn = np.sqrt(1 + (B["uc"] ** 2).sum(1))
            wid = np.concatenate([2 * B["hr"][None, :], (2 * B["hu"] * (B["rc"] / yn)[:, None]).T], 0)
            wid = wid / wid.sum(0)
            # split rule (fixed aspect ratio, calibrated in a probe): radial if 8 hr >= max hu, else the widest angular coordinate
            kdim = np.where(8 * B["hr"] >= B["hu"].max(1), 0, 1 + np.argmax(B["hu"], axis=1))[rej]
            C = {k: v[rej].copy() for k, v in B.items()}
            C2 = {k: v.copy() for k, v in C.items()}
            idx = np.arange(len(kdim))
            rad = kdim == 0
            C["hr"][rad] *= 0.5; C2["hr"][rad] *= 0.5
            C["rc"][rad] -= C["hr"][rad]; C2["rc"][rad] += C2["hr"][rad]
            ang = ~rad
            jj = kdim[ang] - 1; ia = idx[ang]
            C["hu"][ia, jj] *= 0.5; C2["hu"][ia, jj] *= 0.5
            C["uc"][ia, jj] -= C["hu"][ia, jj]; C2["uc"][ia, jj] += C2["hu"][ia, jj]
            Cn = cat(C, C2)
            keep = polar_keep(Cn["face"], Cn["rc"], Cn["hr"], Cn["uc"], Cn["hu"])
            Q = cat(Q, {k: v[keep] for k, v in Cn.items()})
        it += 1
        if it % 5 == 0 and len(Q["hr"]):
            kk = int(np.argmin(Q["hr"]))
            log("      smallest-hr box: face %d rc %.6f hr %.2e uc %s hu %s" % (Q["face"][kk], Q["rc"][kk], Q["hr"][kk], np.round(Q["uc"][kk], 5), Q["hu"][kk]))
        if it % 5 == 0:
            log("  [%6.1fs] evaluated %d accepted %d queue %d  maxUB_R %.6e maxUB_F %.6e  LB_R %.6e LB_F %.6e  pvol_acc %.4f of %.4f  min hr %.1e hu %.1e" % (
                time.time() - t0, st["n_eval"], st["n_acc"], len(Q["rc"]), st["maxUB_R"], st["maxUB_F"], st["maxLB_R"], st["maxLB_F"], st["pvol_acc"], st["pvol_init"],
                Q["hr"].min() if len(Q["hr"]) else 0, Q["hu"].min() if len(Q["hu"]) else 0))
    st["elapsed"] += time.time() - t0
    if len(Q["rc"]) == 0:
        st["done"] = True
        json.dump(st, open(fdone, "w"), indent=1)
        log("interval stage DONE: rigorous sup |Ric|/a^6 <= %.6e, sup |F|/a^8 <= %.6e (evaluated %d boxes)" % (st["maxUB_R"], st["maxUB_F"], st["n_eval"]))
        if os.path.exists(fs):
            os.remove(fs)
        return st
    np.savez(fs, stats=json.dumps(st), **Q)
    log("checkpoint: queue %d, evaluated %d, pvol_acc %.4f of %.4f; rerun 'python rep3.py interval' to resume" % (len(Q["rc"]), st["n_eval"], st["pvol_acc"], st["pvol_init"]))
    return None

# ==================================================================================
# 12. W-orbit refinement, controls, report
# ==================================================================================
def W_orbit_24(v):
    """24 W(F4)-images of a special direction: axis type (+-e_i, (+-1,+-1,+-1,+-1)/2) or (e_i +- e_j)/sqrt2 type."""
    r = np.linalg.norm(v)
    u = v / r
    out = []
    if np.sum(np.abs(u) > 1e-3) in (1, 4):
        for i in range(4):
            for sg in (1, -1):
                e = np.zeros(4); e[i] = sg; out.append(e)
        for sg in itertools.product((1, -1), repeat=4):
            out.append(np.array(sg) / 2.0)
    else:
        for i in range(4):
            for j in range(i + 1, 4):
                for s1 in (1, -1):
                    for s2 in (1, -1):
                        e = np.zeros(4); e[i] = s1; e[j] = s2; out.append(e / math.sqrt(2))
    return np.array(out) * r

def stage_orbits():
    """refine from all 24 W-images of the refined argmax (W-images are G_J-inequivalent in general)."""
    fn = os.path.join(CACHE, "orbits.json")
    if os.path.exists(fn):
        return json.load(open(fn))
    table = stage_coeffs()
    out = {}
    for a in (0.03, 0.01):
        sa = json.load(open(os.path.join(CACHE, "sample_a%g.json" % a)))
        ev = make_sampler(table, a)[0]
        res = {}
        for name in ("absF", "Ric"):
            x0 = np.array(sa["refined_" + name]["argmax_x"])
            u = np.round(x0 / np.abs(x0).max())
            seeds = W_orbit_24(u / np.linalg.norm(u) * np.linalg.norm(x0))
            f = (lambda Y: np.abs(ev(Y)[0])) if name == "absF" else (lambda Y: ev(Y)[1])
            v0 = f(seeds)
            hist, Xf, vf = pattern_refine(f, seeds, s0=0.004, smin=2e-6)
            j = int(np.argmax(vf))
            classes = sorted(set(np.round(vf, 3).tolist()))
            res[name] = {"orbit_values_unrefined": v0.tolist(), "refined_values": vf.tolist(),
                         "distinct_refined_values": classes, "max": float(vf[j]), "argmax_x": Xf[j].tolist(),
                         "argmax_|zeta|": float(np.linalg.norm(Xf[j])), "history": [(h_[0], h_[1]) for h_ in hist],
                         "F_Ric_at_argmax": [float(t[0]) for t in ev(Xf[j][None])]}
            log("a=%g %s: orbit refined max %.10e at %s |zeta|=%.6f ; distinct values %s" % (a, name, vf[j], np.round(Xf[j], 6), np.linalg.norm(Xf[j]), classes))
        out["a=%g" % a] = res
    json.dump(out, open(fn, "w"), indent=1)
    return out

def eh_zeta_check(x):
    """Q4 in the jet engine: det and dd^bar log det of the computed Eguchi-Hanson Hessian jets (zeta coords)."""
    z1, z2 = point_af(x)
    nn, n = z1[0].shape
    w1, w2 = af_conj(z1), af_conj(z2)
    sc = (np.abs(z1[0][0]) ** 2 + np.abs(z2[0][0]) ** 2)
    pC = p1_jet_interval(sc, sc)
    ph = [None] + [(np.array([0.5 * (pC[k - 1][0] + pC[k - 1][1]) / k + 0j]), np.zeros(n)) for k in range(1, 5)]
    Dl = JS4.zeros(nn, n)
    for i, (zi, wi) in enumerate(((z1, w1), (z2, w2))):
        Dl[0][JS4.unit(i)] = wi[0]; Dl[0][JS4.unit(2 + i)] = zi[0]
        t = [0] * 4; t[i] = 1; t[2 + i] = 1
        Dl[0][JS4.pos[tuple(t)]][0] = 1.0
    Pw = [None, Dl]
    for k in range(2, 5):
        Pw.append(JS4.mul(Pw[-1], Dl))
    PhiEH = JS4.series_pw(ph, Pw)
    G = {k: extract(PhiEH, MAP_Z[k]) for k in MAP_Z}
    detG = af_sub(JS2.mul(G[0, 0], G[1, 1]), JS2.mul(G[0, 1], G[1, 0]))
    ld = jet_func2(detG)
    d0 = np.abs(detG[0][0][0] - 1).max()
    ric = max(np.abs(ld[0][JS2.pos[t]][0]).max() for t in ((1, 0, 1, 0), (0, 1, 0, 1), (1, 0, 0, 1)))
    return float(d0), float(ric)

def stage_controls(args):
    fn = os.path.join(CACHE, "controls.json")
    C = json.load(open(fn)) if os.path.exists(fn) else {}
    table = stage_coeffs()
    def save():
        json.dump(C, open(fn, "w"), indent=1)
    if "coeffs" not in C:
        c = {}
        ex, cf = load_poly(table, 1, 4); c["max|P4[psi1]| (L3d: 0)"] = float(np.abs(cf).max())
        ex, cf = load_poly(table, 1, 6); hx, hc = H6_poly(); d = dict(zip(map(tuple, ex), cf))
        rat = np.array([d[tuple(e)] / v for e, v in zip(hx, hc)])
        c["P6[psi1]/H6 ratio spread"] = float(np.ptp(rat.real)); c["c6 from R=8 coefficients (= -2 ratio)"] = float(-2 * rat.real.mean())
        c["max non-H6 coefficient in P6[psi1]"] = float(max(abs(v) for k_, v in d.items() if k_ not in set(map(tuple, hx))))
        c6 = c6_enclosure(); c["c6 rigorous enclosure (R=16 + A4 tail)"] = [c6["c6_mid"] - c6["c6_rad"], c6["c6_mid"] + c6["c6_rad"]]
        ex, cf = load_poly(table, 1, 10); c["max|P10[psi1]| (spec L3(g): 0)"] = float(np.abs(cf).max())
        ex, cf = load_poly(table, 2, 2); d = dict(zip(map(tuple, ex), cf))
        c["P2[psi2] |z|^2 coefficient incl. Z2 term, relative (L3c: 0; R=5 truncation)"] = float(abs(d[(1, 0, 1, 0)] + binom_half(2) * float(Zkappa(2))) / abs(binom_half(2) * float(Zkappa(2))))
        ex, cf = load_poly(table, 2, 4); d = dict(zip(map(tuple, ex), cf))
        c["P4[psi2] |z|^4 coefficient (R=5) vs Z_{5/2}/4"] = [float(d[(2, 0, 2, 0)].real), float(Zkappa(2.5)) / 4]
        worst = 0.0
        for k in range(1, KMAX + 1):
            for dd in range(4, DMAX + 1, 2):
                if (k, dd) in ((1, 4), (1, 10)):
                    continue          # blocks that vanish identically (L3(d), L3(g)); reported above
                ex, cf = load_poly(table, k, dd)
                oo = ((ex[:, 0] + ex[:, 1]) % 2) == 1
                wrong = oo if dd % 4 == 0 else ~oo
                if np.any(wrong):
                    worst = max(worst, float(np.abs(cf[wrong]).max() / np.abs(cf[~wrong]).max()))
        c["max wrong-class / right-class coefficient ratio (spec 1.5: d=2 mod 4 odd-odd, d=0 mod 4 even-even)"] = worst
        C["coeffs"] = c; save(); log("coeff controls", c)
    if "psi_direct" not in C:
        res = []
        for x in ([0.1, 0.0, 0.0, 0.0], [0.05, 0.0, 0.0, 0.0], [0.05, 0.04, -0.06, 0.03]):
            dv = float(mp_psi_direct(x, 0.03)); pv = poly_psi(table, np.array(x), 0.03)
            res.append({"x": x, "direct": dv, "poly": pv, "rel": abs(dv - pv) / abs(pv)})
            log("psi direct vs poly", res[-1])
        C["psi_direct"] = res; save()
    if "Q4_EH" not in C:
        rng = np.random.default_rng(3)
        x = rng.normal(size=(200, 4)); x /= np.linalg.norm(x, axis=1)[:, None]; x *= rng.uniform(0.2, 1.0, size=(200, 1))
        d0, ric = eh_zeta_check(x)
        polys = build_polys(table, 0.03, 0.03, 14, False)
        TB = chart_poly_T(polys)
        u = rng.normal(size=50) * 0.01 + 1j * rng.normal(size=50) * 0.01; u[:10] = 0
        w = rng.uniform(-1, 1, 50) + 1j * rng.uniform(-1, 1, 50)
        fz, rz, ric_c = chart_pipeline(u, w, TB, 0.03, with_E=False)
        polys0 = {"A": (polys["A"][0][:0], polys["A"][1][:0], polys["A"][2][:0]), "B": (polys["B"][0][:0], polys["B"][1][:0], polys["B"][2][:0])}
        mpF, mpR = mp_direct([0.7, 0.2, -0.3, 0.1], 0.03, polys0)
        C["Q4_EH"] = {"zeta: max|det g0 - 1|": d0, "zeta: max|dd^bar log det g0|": ric,
                      "chart: max|F|/a^8 with E off": float(np.abs(fz).max()), "chart: max |Ric|/a^6 with E off": float(np.abs(rz).max()),
                      "chart: max|dd^bar log det g0^ch| (pluriharmonic)": float(ric_c),
                      "mp direct (FD Hessian + log det), lattice terms off: raw F, raw chart spectral radius": [mpF * 0.03 ** 8, mpR * 0.03 ** 8]}
        save(); log("Q4", C["Q4_EH"])
    if "mp_direct" not in C:
        res = []
        orb = stage_orbits()
        for a in (0.03, 0.01):
            polys = build_polys(table, a, a, 14, False)
            T = prepare_T(polys, False)
            pts = [orb["a=%g" % a]["Ric"]["argmax_x"], orb["a=%g" % a]["absF"]["argmax_x"], [0.5, 0.3, -0.4, 0.2], [0.27, 0.05, 0.03, 0.01]]
            for x in pts:
                F1, r1 = zeta_pipeline(*point_af(np.array([x])), T, a, a, False)
                F2, r2 = mp_direct(x, a, polys)
                res.append({"a": a, "x": x, "jet_F/a8": float(F1[0]), "mp_F/a8": F2, "jet_Ric/a6": float(r1[0]), "mp_Ric/a6": r2,
                            "relF": abs(F1[0] - F2) / abs(F2), "relRic": abs(r1[0] - r2) / abs(r2)})
                log("jet vs mp", res[-1])
        C["mp_direct"] = res; save()
    if "chart_overlap" not in C:
        res = {}
        for a in (0.03, 0.01):
            ev, evc, polys, T, TB = make_sampler(table, a)
            rng = np.random.default_rng(5)
            x = rng.normal(size=(300, 4)); x /= np.linalg.norm(x, axis=1)[:, None]; x *= rng.uniform(0.12, 0.25, size=(300, 1))
            F1, r1 = zeta_pipeline(*point_af(x), T, a, a, False)
            z1 = x[:, 0] + 1j * x[:, 1]; z2 = x[:, 2] + 1j * x[:, 3]
            sw = np.abs(z2) > np.abs(z1); z1, z2 = np.where(sw, z2, z1), np.where(sw, z1, z2)
            F2, r2, _ = chart_pipeline(z1 ** 2, z2 / z1, TB, a)
            Fe, Re, _ = chart_pipeline(np.zeros(3, complex), np.array([0, 0.5, 1j]), TB, a)
            res["a=%g" % a] = {"max rel diff F (zeta vs chart, 0.12<|zeta|<0.25)": float(np.max(np.abs(F1 - F2) / np.abs(F2))),
                               "max rel diff Ric": float(np.max(np.abs(r1 - r2) / np.abs(r2))),
                               "F/a^10 on E (a->0 limit Z_{5/2} = %.6f, spec L4(c))" % float(Zkappa(2.5)): float(Fe[0] / a ** 2)}
        C["chart_overlap"] = res; save(); log("chart", res)
    if "symmetry" not in C:
        ev = make_sampler(table, 0.03)[0]
        rng = np.random.default_rng(9)
        x = rng.normal(size=(100, 4)); x /= np.linalg.norm(x, axis=1)[:, None]; x *= rng.uniform(0.3, 1.0, size=(100, 1))
        F0, R0 = ev(x)
        q = np.array([0.5, 0.5, 0.5, 0.5])
        def rmul(y):
            a0, a1, a2, a3 = y.T; b0, b1, b2, b3 = q
            return np.stack([a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3, a0 * b1 + a1 * b0 + a2 * b3 - a3 * b2,
                             a0 * b2 - a1 * b3 + a2 * b0 + a3 * b1, a0 * b3 + a1 * b2 - a2 * b1 + a3 * b0], 1)
        def lmul(y):
            b0, b1, b2, b3 = y.T; a0, a1, a2, a3 = q
            return np.stack([a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3, a0 * b1 + a1 * b0 + a2 * b3 - a3 * b2,
                             a0 * b2 - a1 * b3 + a2 * b0 + a3 * b1, a0 * b3 + a1 * b2 - a2 * b1 + a3 * b0], 1)
        maps = {"z1->i z1 (in G_J)": lambda y: np.stack([-y[:, 1], y[:, 0], y[:, 2], y[:, 3]], 1),
                "swap z1,z2 (in G_J)": lambda y: y[:, [2, 3, 0, 1]],
                "conj both (in G_J)": lambda y: np.stack([y[:, 0], -y[:, 1], y[:, 2], -y[:, 3]], 1),
                "x0<->x1 only (W, not G_J)": lambda y: y[:, [1, 0, 2, 3]],
                "right mult by (1+i+j+k)/2 (W; commutes with J, so in G_J)": rmul,
                "left mult by (1+i+j+k)/2 (W, not G_J)": lmul}
        res = {}
        for k_, f in maps.items():
            F1, R1 = ev(f(x))
            res[k_] = {"max rel |dF|": float(np.max(np.abs(F1 - F0) / np.abs(F0))), "max rel |dRic|": float(np.max(np.abs(R1 - R0) / np.abs(R0)))}
        C["symmetry"] = res; save(); log("symmetry", res)
    if "E_refined" not in C:
        res = {}
        for a in (0.03, 0.01):
            ev, evc, polys, T, TB = make_sampler(table, a)
            g = np.linspace(-1, 1, 61)
            W = (g[:, None] + 1j * g[None, :]).ravel(); W = W[np.abs(W) <= 1.0]
            out = [evc(np.zeros(len(W[i:i + 2000]), complex), W[i:i + 2000]) for i in range(0, len(W), 2000)]
            Fv = np.concatenate([o[0] for o in out]); Rv = np.concatenate([o[1] for o in out])
            k = int(np.argmax(Rv)); w0 = W[k]; st_ = 0.01; best = Rv[k]
            while st_ > 1e-7:
                cand = w0 + st_ * np.exp(2j * np.pi * np.arange(8) / 8)
                f_, r_ = evc(np.zeros(8, complex), cand)
                j = int(np.argmax(r_))
                if r_[j] > best:
                    best, w0 = r_[j], cand[j]
                else:
                    st_ /= 2
            res["a=%g" % a] = {"max_Ric_over_a6_on_E": float(best), "argmax_w": [float(w0.real), float(w0.imag)],
                               "max_Ric_over_a8_on_E": float(best / a ** 2), "F_over_a10_on_E": float(Fv[k] / a ** 2)}
        C["E_refined"] = res; save(); log("E refined", res)
    if "old_window" not in C:
        global WIN_LO
        WIN_LO = 0.5
        try:
            ev = make_sampler(table, 0.01)[0]
            X, _, _ = sample_points()
            F, Rr = ev(X)
            out = {}
            for name, vals in (("absF", np.abs(F)), ("Ric", Rr)):
                seeds = X[pick_seeds(X, vals)]
                f = (lambda Y: np.abs(ev(Y)[0])) if name == "absF" else (lambda Y: ev(Y)[1])
                hist, Xf, vf = pattern_refine(f, seeds, smin=1e-5)
                j = int(np.argmax(vf))
                out[name] = {"max": float(vf[j]), "argmax_x": Xf[j].tolist(), "|zeta|": float(np.linalg.norm(Xf[j]))}
            C["old_window"] = {"window": "[a/2, a] (spec v2 design), a = 0.01, positive control", **out}
        finally:
            WIN_LO = 0.25
        save(); log("old window", C["old_window"])
    return C

def stage_report():
    S = stage_sample(); O = stage_orbits(); C = stage_controls([])
    fdone = os.path.join(CACHE, "interval_result.json")
    I = json.load(open(fdone)) if os.path.exists(fdone) else None
    def fmt(v, n=9):
        return ("%." + str(n) + "e") % v
    def fmt_up(v, n=6):
        """decimal string >= v (upper bounds are rounded up, never to nearest)."""
        e = math.floor(math.log10(abs(v)))
        q = 10.0 ** (e - n)
        return ("%." + str(n) + "e") % (math.ceil(v / q * (1 + 1e-15)) * q)
    res = {"seat": "R3",
           "object": "K3 Kummer approximant near one bubble, cut-off window [a/4, a] (C-infinity step), core |zeta| <= 1",
           "norm_convention": "F = log det(g_{i jbar}) with g_{i jbar} = d_i dbar_j Phi_a (det g_EH = 1, flat reference dz1^dz2); |Ric|_g = max |eigenvalue| of the endomorphism g^{-1} Ric (operator norm w.r.t. g), Ric_{i jbar} = -d_i dbar_j F"}
    for a in (0.03, 0.01):
        sa = S["a=%g" % a]; oa = O["a=%g" % a]
        best = {}
        for name in ("absF", "Ric"):
            r1 = sa["refined_" + name]; r2 = oa[name]
            best[name] = r2 if r2["max"] >= r1["max"] else r1
        res["a=%g" % a] = {
            "n_sample_points": S["n_points_zeta"] + S["n_points_chart_E"],
            "sampled_max_absF_over_a8": fmt(sa["sampled_max_absF_over_a8"]),
            "sampled_max_Ric_over_a6": fmt(sa["sampled_max_Ric_over_a6"]),
            "refined_max_absF_over_a8": fmt(best["absF"]["max"]),
            "refined_argmax_absF_zeta_x0x1x2x3": [round(t, 6) for t in best["absF"]["argmax_x"]],
            "refined_argmax_absF_|zeta|": "%.6f" % best["absF"]["argmax_|zeta|"],
            "refined_max_Ric_over_a6": fmt(best["Ric"]["max"]),
            "refined_argmax_Ric_zeta_x0x1x2x3": [round(t, 6) for t in best["Ric"]["argmax_x"]],
            "refined_argmax_Ric_|zeta|": "%.6f" % best["Ric"]["argmax_|zeta|"],
            "inner_|zeta|<1/4_max_absF_over_a8": fmt(sa["inner_|zeta|<1/4_max_absF"], 6),
            "inner_|zeta|<1/4_max_Ric_over_a6": fmt(sa["inner_|zeta|<1/4_max_Ric"], 6),
            "on_E_F_over_a8": fmt(sa["on_E_F_over_a8_range"][1]),
            "on_E_max_Ric_over_a6_sampled": fmt(sa["on_E_Ric_over_a6_max"], 6),
            "on_E_max_Ric_over_a6_refined": fmt(C["E_refined"]["a=%g" % a]["max_Ric_over_a6_on_E"], 9)}
    if I is not None:
        res["rigorous_window"] = {"zone": "1/4 <= |zeta| <= 1, every a in [0.029, 0.03]",
                                  "Ric_over_a6_upper_bound": fmt_up(I["maxUB_R"]), "absF_over_a8_upper_bound": fmt_up(I["maxUB_F"]), "machine_values_R_F": [I["maxUB_R"], I["maxUB_F"]],
                                  "boxes_evaluated": I["n_eval"], "boxes_accepted": I["n_acc"], "unresolved_boxes": I["n_unresolved"],
                                  "max_box_centre_value_Ric_a0.03": "%.9e" % I["maxLB_R"], "max_box_centre_value_absF_a0.03": "%.9e" % I["maxLB_F"],
                                  "acceptance_thresholds_R_F": [I["thrR"], I["thrF"]], "accepted_param_volume_vs_initial": [I["pvol_acc"], I["pvol_init"]], "elapsed_s": I["elapsed"], "chunks": I["runs"]}
        sa = S["a=0.03"]; oa = O["a=0.03"]
        mR = max(sa["sampled_max_Ric_window"], oa["Ric"]["max"], sa["refined_Ric"]["max"], I["maxLB_R"])
        mF = max(sa["sampled_max_absF_window"], oa["absF"]["max"], sa["refined_absF"]["max"], I["maxLB_F"])
        res["Q5_samples_below_bound"] = {"max window sample Ric/a6 (a=0.03)": "%.9e" % mR, "bound_R": fmt_up(I["maxUB_R"]), "ok_R": bool(mR <= I["maxUB_R"]),
                                          "max window sample |F|/a8 (a=0.03)": "%.9e" % mF, "bound_F": fmt_up(I["maxUB_F"]), "ok_F": bool(mF <= I["maxUB_F"])}
    else:
        res["rigorous_window"] = "NOT FINISHED (resume with: python rep3.py)"
    res["controls"] = C
    json.dump(res, open(os.path.join(HERE, "rep3_result.json"), "w"), indent=1)
    log(json.dumps({k: v for k, v in res.items() if k != "controls"}, indent=1))
    return res

if __name__ == "__main__":
    st = sys.argv[1] if len(sys.argv) > 1 else "all"
    if st == "coeffs":
        stage_coeffs()
    elif st == "sample":
        stage_sample(tuple(float(t) for t in sys.argv[2:]) or (0.03, 0.01))
    elif st == "interval":
        stage_interval()
    elif st == "orbits":
        stage_orbits()
    elif st == "controls":
        stage_controls(sys.argv[2:])
    elif st == "report":
        stage_report()
    elif st == "all":
        stage_coeffs(); stage_sample(); stage_orbits(); stage_controls([])
        r = stage_interval()
        if r is None:
            log("interval stage not finished in this run (checkpointed); rerun 'python rep3.py' to resume")
        stage_report()
