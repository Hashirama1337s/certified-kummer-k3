#!/usr/bin/env python
"""wp2b.py -- WP2b certificate build for the K3 Kummer approximant (DEFINITIONS-AND-LEMMAS-v2.md, window
|zeta| in [0.25, 1] with the C-infinity step, c(sigma) = st((sqrt(sigma) - 1/4)/(3/4))).

Certifies, for a in sub-intervals covering (0, a0], a0 = 0.03, upper bounds over the WHOLE manifold of
  e(a)          = sup |E|_{g0}                                    (g0-operator norm, Fix 4)
  a^2 ||dR||(a) = a^2 sup over g0-unit X, Y of |dR(X,Xbar,Y,Ybar)|  (L4 (g))
  a^2 B+(a)     = a^2 sup Bisec(omega_a)                           (L4 (g))
  a^2 C0(a)     = -a^2 inf over g-orthogonal pairs of Bisec         (mirror lemma, BUILD-WP2b.md)
and a certified LOWER bound for a^2 B+ (and a^2 C0) at a = 0.03 from the point (u, w) = (0, 0) on E.

Zones (nearest bubble at 0): Z1 |zeta| <= 1/4 in the (u, w) chart (contains E); Z2 1/4 <= |zeta| <= 1 in zeta;
NB (neck + bulk) |z| >= a inside the Voronoi cell V0, in z.  All core/neck bounds are a-scaled (a enters only as
an interval through explicit series), the bulk cover is a-independent (eps = a^4 enters at assembly).

Modes
  --controls OUT.json   K1 (lattice off), K2 (2000 random points per a in {0.01, 0.03} vs a double-precision
                        reference built from probe1/k1c.py + k1d.py with the C-infinity window), K3 (lattice x10),
                        NaN/inf guards, plus CHECKs (chamber, Ewald vs k1c, Z_kappa vs closed form, Lemma C).
  --main OUT.json       refuses unless <lab>/wp2/GO-WP2b exists; resumable (state_wp2b/), exit 3 = resume.
Run:  python <lab-tools>/memcap.py 8 python <lab>/wp2/wp2b.py --controls <lab>/wp2/controls_wp2b.json
"""
import os, sys, json, time, math, argparse, itertools, pickle
sys.dont_write_bytecode = True                      # importing the sealed wp2a.py must not write files
import numpy as np
from flint import arb, arb_series, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import wp2a as WA                                   # sealed: read-only import (interval class, step, Lemma C, W)
ctx.prec = 128

GO_FILE = os.path.join(HERE, 'GO-WP2b')
STATE_DIR = os.path.join(HERE, 'state_wp2b')
I, isq, isqrt, irecip, iabs = WA.I, WA.isq, WA.isqrt, WA.irecip, WA.iabs
# Faster outward rounding than np.nextafter (same role, larger margin): every IEEE +, -, *, /, sqrt result q of an
# exact value p satisfies |q - p| <= 2^-53 |p| (normal range) or <= 2^-1075 (subnormal); moving q by
# 2^-51 |q| + 2^-1060 (one more correctly rounded subtraction) lands strictly on the safe side of p (factor-4 margin
# over the needed 2^-53 relative, and the absolute term covers zero and subnormals).  Installed in-memory in the
# imported wp2a module too (its interval class calls these); wp2a.py itself is not modified.
_E51 = 2.0 ** -51; _TINY = 2.0 ** -1060
def _dn(x): return x - (np.abs(x) * _E51 + _TINY)
def _up(x): return x + (np.abs(x) * _E51 + _TINY)
WA._dn, WA._up = _dn, _up
R1W, WW = 0.25, 0.75                                # window: T = (|zeta| - R1W)/WW
ETA = 16.0                                          # Ewald split used here (A6 holds for any fixed eta > 0)
VINV = 32.0                                         # 1/V, V = covolume of L = 1/32
U53 = 2.0 ** -53

# =====================================================================================================
# 0. small interval utilities (all outward rounded)
# =====================================================================================================
def cI(x):
    """interval constant from an arb ball or a float (one ulp outward)."""
    if isinstance(x, arb):
        lo, hi = WA.a2f(x); return I(lo, hi)
    return I(_dn(float(x)), _up(float(x)))

def isum(x, axis):
    """sum of an interval array along an axis; float summation error bounded by gamma_n * sum |a_i|
    (Higham, any summation order), with a 1% margin and an absolute 1e-300 underflow allowance."""
    n = x.lo.shape[axis]
    g = 1.01 * n * U53
    slo = np.sum(x.lo, axis=axis); shi = np.sum(x.hi, axis=axis)
    elo = g * np.sum(np.abs(x.lo), axis=axis) + 1e-300
    ehi = g * np.sum(np.abs(x.hi), axis=axis) + 1e-300
    return I(_dn(slo - elo), _up(shi + ehi))

def imulp(g, m, exact=True):
    """interval g times a float array m.  exact=True: m is the exact value; exact=False: m = fl(exact product) with
    relative error <= 2^-53, covered by an extra 2^-51 relative widening."""
    a = g.lo * m; b = g.hi * m
    lo = np.minimum(a, b); hi = np.maximum(a, b)
    if not exact:
        lo = lo - np.abs(lo) * _E51; hi = hi + np.abs(hi) * _E51
    return I(_dn(lo), _up(hi))

def idot(g, M, inexact=None):
    """sum over the middle axis of interval g (nc, n) times exact float M (nc, n, q) or (n, q) -> I (nc, q).
    Positive/negative parts of M give lo = g.lo M+ - g.hi M-, hi = g.hi M+ - g.lo M-; each matrix product is a
    dot product with |fl(x.y) - x.y| <= gamma_n |x|.|y| (Higham Thm 3.1, any order/FMA), bounded by
    1.02 (n+3) u sum |g| |M| (+ u sum |g||M| for columns flagged inexact: M itself carries one rounding)."""
    Mp = np.maximum(M, 0.0); Mm = np.maximum(-M, 0.0); Ma = np.abs(M)
    gm = np.maximum(np.abs(g.lo), np.abs(g.hi))
    if M.ndim == 2:
        A1 = g.lo @ Mp; B1 = g.hi @ Mm; A2 = g.hi @ Mp; B2 = g.lo @ Mm; E = gm @ Ma
    else:
        mm = lambda a, b: np.matmul(a[:, None, :], b)[:, 0, :]
        A1 = mm(g.lo, Mp); B1 = mm(g.hi, Mm); A2 = mm(g.hi, Mp); B2 = mm(g.lo, Mm); E = mm(gm, Ma)
    n = g.lo.shape[1]
    err = 1.02 * (n + 3) * U53 * E + 1e-300 * (n + 1)
    if inexact is not None:
        err = err + 1.01 * U53 * E * inexact[None, :]
    err = _up(err)
    return I(_dn(_dn(A1 - B1) - err), _up(_up(A2 - B2) + err))

def imag_(x):
    return np.maximum(np.abs(x.lo), np.abs(x.hi))

def iscale(x, c):
    return x * float(c)

def ipow(x, n):
    out = None
    for _ in range(n):
        out = x if out is None else out * x
    return out

# rigorous exp on nonpositive arguments: exp(y)^(2^M), y = x 2^-M, Taylor degree 8 + Lagrange remainder
_EXPM = 16
_INVF = [cI(arb(1) / math.factorial(k)) for k in range(9)]
def _exp_point(x):
    """enclosure of exp(x) for a float array x <= 0 (|x| <= 700)."""
    y = I(np.ldexp(x, -_EXPM))                       # exact power-of-two scaling
    p = _INVF[8]
    for k in range(7, -1, -1):
        p = p * y + _INVF[k]
    ay = np.abs(np.ldexp(x, -_EXPM))
    rem = _up(_up(ay ** 9) / 362880.0 * 1.01)        # |y|^9/9! e^|y|, |y| <= 0.011
    p = I(_dn(p.lo - rem), _up(p.hi + rem))
    p = I(np.maximum(p.lo, 0.0), p.hi)
    for _ in range(_EXPM):
        p = isq(p)
    return p

def iexp_neg(x):
    """exp of an interval x with x.hi <= 0 (monotone: endpoints)."""
    if np.any(x.hi > 0):
        raise FloatingPointError('iexp_neg: positive argument')
    if np.any(x.lo < -700):
        raise FloatingPointError('iexp_neg: argument below -700')
    return I(_exp_point(x.lo).lo, _exp_point(x.hi).hi)

_TWOPI = cI(2 * arb.pi())
_SC_INV = [cI(arb(1) / math.factorial(k)) for k in range(18)]
def cos_sin_2pi(t):
    """enclosures of cos(2 pi t), sin(2 pi t) for exact floats t in [-1/2, 1/2]: Taylor at phi = 2 pi t / 8
    (|phi| <= pi/8, degree 16/17 with remainder |phi|^18/18!), then three double-angle steps."""
    ph = _TWOPI * I(np.ldexp(t, -3))
    p2 = isq(ph)
    c = _SC_INV[16]; s = _SC_INV[17]
    for k in range(14, -1, -2):
        c = c * p2 * (-1.0) + _SC_INV[k]
    for k in range(15, 0, -2):
        s = s * p2 * (-1.0) + _SC_INV[k]
    s = s * ph
    rem = 1e-21                                       # (pi/8)^18/18! = 3.0e-24, margin
    c = I(_dn(c.lo - rem), _up(c.hi + rem)); s = I(_dn(s.lo - rem), _up(s.hi + rem))
    for _ in range(3):
        s, c = (s * c) * 2.0, 1.0 - isq(s) * 2.0
    return c, s

def check_finite(name, *xs):
    for k, x in enumerate(xs):
        if isinstance(x, I):
            WA.check_finite(name, x)
        else:
            a = np.asarray(x, float)
            if not np.all(np.isfinite(a)):
                raise FloatingPointError('%s[%d]: non-finite value' % (name, k))

def aup(x):
    if not x.is_finite():
        raise FloatingPointError('non-finite arb %s' % x)
    return float(_up(float(x.upper())))

def adn(x):
    if not x.is_finite():
        raise FloatingPointError('non-finite arb %s' % x)
    return float(_dn(float(x.lower())))

def aball(lo, hi):
    return arb(lo).union(arb(hi))

# =====================================================================================================
# 1. lattice L = Lambda/2 (Hurwitz/2), Z_kappa, A4 tails, Cauchy-series constants of L3 (i)
# =====================================================================================================
def lattice_points(R):
    """all l in L with |l| <= R (as float array, exact quarters); l = m/4, m in Z^4 all-even or all-odd."""
    M = int(math.floor(4 * R)) + 1
    rng = np.arange(-M, M + 1)
    out = []
    for m0 in rng:
        g = np.stack(np.meshgrid([m0], rng, rng, rng, indexing='ij'), -1).reshape(-1, 4)
        par = g % 2
        ok = np.all(par == par[:, :1], axis=1)
        g = g[ok]
        g = g[(g * g).sum(1) <= 16 * R * R]
        out.append(g)
    m = np.concatenate(out)
    return m / 4.0

def A4_tail(T, p):
    """Lemma A4: sum_{|x-l| >= T} |x-l|^-p <= 1024 sum_j C(3,j) 4^(j-3) (T-1/2)^(j+1-p)/(p-j-1)  (arb)."""
    T = arb(T) if not isinstance(T, arb) else T
    assert T > arb(0.5) and p > 4
    return 1024 * sum(math.comb(3, j) * arb(4) ** (j - 3) * (T - arb(0.5)) ** (j + 1 - p) / (p - j - 1) for j in range(4))

_ZCACHE = {}
def Z_enclosures(Nmax=1600):
    """certified Z_kappa = sum_{l != 0} |l|^(-4 kappa) for kappa in {2, 2.5, 3}: exact shell counts of |m|^2 = N
    (|l|^2 = N/16) for N <= Nmax, plus the A4 tail beyond |l| = sqrt((Nmax+1)/16)."""
    if Nmax in _ZCACHE:
        return _ZCACHE[Nmax]
    R = math.sqrt(Nmax / 16.0)
    m = np.rint(lattice_points(R) * 4).astype(np.int64)
    N = (m * m).sum(1)
    cnt = np.bincount(N, minlength=Nmax + 1)
    out = {}
    T = (arb(Nmax + 1) / 16).sqrt()
    for kap in (2.0, 2.5, 3.0):
        s = arb(0)
        for n in range(1, Nmax + 1):
            if cnt[n]:
                s += int(cnt[n]) * (arb(n) / 16) ** (-2 * arb(kap))
        tail = A4_tail(T, 4 * kap)
        out[kap] = s + arb(0).union(tail)
    out['shell_counts_to_100'] = {int(n): int(cnt[n]) for n in range(1, 101) if cnt[n]}
    _ZCACHE[Nmax] = out
    return out

def Z_closed(kap):
    """CHECK only (shell-count formula is STANDARD in general): 24 4^(2k) (1 - 2^(1-2k)) zeta(2k) zeta(2k-1)."""
    k = arb(kap)
    return 24 * arb(4) ** (2 * k) * (1 - arb(2) ** (1 - 2 * k)) * (2 * k).zeta() * (2 * k - 1).zeta()

def binom_half(k):
    """|binom(1/2, k)| exactly as arb."""
    v = arb(1)
    for i in range(k):
        v = v * (arb(1) / 2 - i) / (i + 1)
    return abs(v)

class Series:
    """Upper bounds for the holomorphic extensions (Lie(y) <= rho, zeta units) of the scaled potential pieces
       A_oo = sum_{d = 2 mod 4, d >= 6} [a^(d-6) P_d[psi_1] + a^(d-2) P_d[psi~_2]]
       A_ee = sum_{d = 0 mod 4}         [a^(d-6) P_d[psi_1] (d >= 8) + a^(d-2) P_d[psi~_2] (d >= 4)]
       using L3 (i), Z_{kappa + D} <= 16^D Z_kappa (|l| >= 1/2) and the binomial-series lemma of BUILD-WP2b.md.
       a is an arb ball (may contain 0)."""
    def __init__(self, Z):
        self.Z = Z
    def MA(self, rho, a, dmin):
        """sum_{d >= dmin} a^(d-6) C1(d) rho^d,  C1(d) = (d+1)/2 Z_{(d+2)/4}."""
        Z0 = self.Z[(dmin + 2) / 4.0].upper()
        Z0 = arb(Z0)
        y = 2 * a * rho
        if not (y.upper() < 1): raise FloatingPointError('MA: 2 a rho >= 1')
        D = dmin
        return Z0 / 2 * a ** (D - 6) * rho ** D * ((D + 1) / (1 - y) + y / (1 - y) ** 2)
    def MB(self, rho, a, dmin):
        """sum_{d >= dmin} a^(d-2) C2(d) rho^d, C2(d) = sum_{k>=2} |b_k| eps^(k-2) C(d+4k-3,d) Z_{(4k-2+d)/4}/(2k-1)."""
        Z0 = arb(self.Z[(6 + dmin) / 4.0].upper())
        y = 2 * a * rho
        if not (y.upper() < 1): raise FloatingPointError('MB: 2 a rho >= 1')
        eps = a ** 4
        D = dmin
        def Tk(k):
            return (binom_half(k) / (2 * k - 1) * eps ** (k - 2) * arb(2) ** (4 * k - 8)
                    * math.comb(D + 4 * k - 3, D) * (1 - y) ** (-(4 * k - 2 + D)))
        q = eps * 16 * arb(math.comb(D + 13, D)) / math.comb(D + 9, D) * (1 - y) ** (-4)
        if not (q.upper() < 1): raise FloatingPointError('MB: q >= 1')
        return Z0 * a ** (D - 2) * rho ** D * (Tk(2) + Tk(3) / (1 - q))
    def Moo(self, rho, a):
        return self.MA(rho, a, 6) + self.MB(rho, a, 6)
    def Mee(self, rho, a):
        return self.MA(rho, a, 8) + self.MB(rho, a, 4)

def cauchy_deriv(Mfun, r, m, a, deltas=(0.05, 0.1, 0.15, 0.2, 0.3, 0.45, 0.7, 1.0, 1.5, 2.2)):
    """||D^m Q(x)|| (real symmetric m-linear norm, |x| <= r) <= delta^-m M(r + m delta), minimised over delta
    (multi-direction Cauchy estimate; Lie(x + sum t_i v_i) <= |x| + m delta)."""
    if m == 0:
        return Mfun(arb(r), a)
    best = None
    for d in deltas:
        rho = arb(r) + m * arb(d)
        if (2 * a * rho).upper() >= 0.999:
            continue
        v = Mfun(rho, a) / arb(d) ** m
        if best is None or v.upper() < best.upper():
            best = v
    if best is None:
        raise FloatingPointError('cauchy_deriv: no admissible delta')
    return best

# =====================================================================================================
# 2. the step: c^(k)(sigma), k = 0..4, on a T-piece (arb).  c(sigma) = st(T), T = (sqrt(sigma) - 1/4)/(3/4)
# =====================================================================================================
def c_derivs(T0, T1, step='cinf'):
    """arb balls containing c^(k)(sigma), k = 0..4, for all sigma with T in [T0, T1] (0 <= T0 < T1 <= 1).
    Interior pieces: natural arb series enclosure.  Pieces touching T = 0 or 1 (width <= 1/8): Lemma C (C1),
    (C5), (C6) for st^(k), then Faa di Bruno with the T(sigma) derivatives."""
    s0 = arb(R1W) + arb(WW) * aball(T0, T1)
    sig = s0 * s0
    X = arb_series([sig, 1], prec=5)
    Ts = (X.sqrt() - arb(R1W)) / arb(WW)
    endp = (step == 'cinf') and (T0 <= 0 or T1 >= 1)
    if not endp:
        C = WA.st_of_series(Ts, step)
        cf = C.coeffs() + [arb(0)] * 5
        return [cf[k] * math.factorial(k) for k in range(5)]
    assert T1 - T0 <= 0.125
    if T0 <= 0:
        Bk = WA.lemmaC_bounds(T1)
        S = [arb(0).union(Bk[0])] + [Bk[k].union(-Bk[k]) for k in range(1, 5)]
    else:
        Bk = WA.lemmaC_bounds(1 - T0)
        S = [(1 - Bk[0]).union(arb(1))] + [Bk[k].union(-Bk[k]) for k in range(1, 5)]
    tc = Ts.coeffs() + [arb(0)] * 6
    T1d, T2d, T3d, T4d = tc[1], 2 * tc[2], 6 * tc[3], 24 * tc[4]
    c0 = S[0]; c1 = S[1] * T1d; c2 = S[2] * T1d ** 2 + S[1] * T2d
    c3 = S[3] * T1d ** 3 + 3 * S[2] * T1d * T2d + S[1] * T3d
    c4 = S[4] * T1d ** 4 + 6 * S[3] * T1d ** 2 * T2d + S[2] * (3 * T2d ** 2 + 4 * T1d * T3d) + S[1] * T4d
    return [c0, c1, c2, c3, c4]

def radial_norms(cd, r):
    """real tensor norms of D^j c(|x|^2), j = 0..4, at |x| <= r, from |c^(k)| (Faa di Bruno, D^3 sigma = 0)."""
    m = [abs(v) for v in cd]
    r = arb(r)
    return [m[0], 2 * r * m[1], 4 * r * r * m[2] + 2 * m[1], 8 * r ** 3 * m[3] + 12 * r * m[2],
            16 * r ** 4 * m[4] + 48 * r * r * m[3] + 12 * m[2]]

# =====================================================================================================
# 3. pointwise curvature lemma (L4 (g) and its mirror): combination of the jet norms
# =====================================================================================================
def dR_bound(N2, N3, N4, Gam):
    """||dR||_{g0} <= N4 + N3 (Gam + N3)/(1-e) + Gam N3/(1-e) + Gam^2 e/(1-e),  e = N2 < 1  (arb, BUILD lemma D)."""
    if not (N2.upper() < 1):
        raise FloatingPointError('e >= 1: positivity not certified')
    q = 1 / (1 - N2)
    return N4 + N3 * (Gam + N3) * q + Gam * N3 * q + Gam * Gam * N2 * q

def bis_bound(a2K, a2D, e):
    """a^2 Bisec <= (a^2 K + a^2 ||dR||)/(1 - e)^2 and a^2 Bisec >= -(same)  (L1, L4 (g), mirror lemma)."""
    return (a2K + a2D) / (1 - e) ** 2

def C0_lower(a2D, e):
    """certified lower bound for a^2 C0 at a point of E where K = 2/a^2 (BUILD lemma M2)."""
    c = e / (1 - e)
    num = 2 - a2D - (2 * c + c * c) * (2 + a2D)
    return num / ((1 + e) ** 2 * (1 + c) ** 2)

# =====================================================================================================
# 4. zone Z2: 1/4 <= |zeta| <= 1, zeta coordinates, 1-D cover in T
# =====================================================================================================
def zone_Z2(SER, a, nT=128, scale=1.0, log=None):
    """returns per-piece arb upper bounds (lists) and the maxima.  a: arb ball (a-interval)."""
    out = dict(e=arb(0), a2D=arb(0), a2bis=arb(0), pieces=[])
    edges = [k / nT for k in range(nT + 1)]
    a8 = a ** 8
    sc = arb(scale)
    for k in range(nT):
        T0, T1 = edges[k], edges[k + 1]
        r0 = R1W + WW * T0; r1 = R1W + WW * T1
        cd = c_derivs(T0, T1)
        Dc = radial_norms(cd, r1)
        # Cauchy bounds of D^m A_oo, D^m A_ee at |zeta| <= r1
        Doo = [cauchy_deriv(SER.Moo, r1, m, a) for m in range(5)]
        Dee = [cauchy_deriv(SER.Mee, r1, m, a) for m in range(5)]
        n = {}
        for kk in (2, 3, 4):
            n[kk] = sum(math.comb(kk, j) * Dc[j] * Doo[kk - j] for j in range(kk + 1)) + Dee[kk]
        sig = aball(r0 * r0, r1 * r1)
        fh = (sig * sig + 1).sqrt()
        th = arb(adn(sig / fh))                       # smallest gh0 eigenvalue sigma/fh (lower end)
        Gam = arb(aup(arb(3).sqrt() / (sig * fh ** arb(1.5))))   # ||d gh0||_{g0} (exact EH frame components)
        a2K = arb(aup(2 / (sig * sig + 1) ** arb(1.5)))
        e = sc * a8 * n[2] / th
        N3 = sc * a8 * n[3] / th ** arb(1.5)
        N4 = sc * a8 * n[4] / th ** 2
        e_u = arb(aup(e)); D = dR_bound(e_u, arb(aup(N3)), arb(aup(N4)), Gam)
        bb = bis_bound(a2K, D, e_u)
        rec = dict(T=(T0, T1), e=aup(e), a2D=aup(D), a2bis=aup(bb), a2K=aup(a2K), n2=aup(n[2]), n3=aup(n[3]),
                   n4=aup(n[4]), theta_lo=adn(th))
        out['pieces'].append(rec)
        for key in ('e', 'a2D', 'a2bis'):
            if rec[key] > float(out[key].upper()):
                out[key] = arb(rec[key])
    return out

# =====================================================================================================
# 5. zone Z1: |zeta| <= 1/4, (u, w) chart, 2-D cover in (|u|, |w|) (EH chart quantities by order-3 jets)
# =====================================================================================================
EXP3 = [e for d in range(4) for e in itertools.product(range(4), repeat=4) if sum(e) == d]
EIDX3 = {e: i for i, e in enumerate(EXP3)}

class TPS:
    """truncated (order 3) jets in 4 real variables with interval coefficients over a batch of boxes."""
    def __init__(self, c):
        self.c = c                                   # dict exponent -> I
    @staticmethod
    def var(i, base):
        e0 = (0, 0, 0, 0); e1 = tuple(1 if j == i else 0 for j in range(4))
        n = len(base.lo)
        return TPS({e0: base, e1: I(np.ones(n))})
    @staticmethod
    def const(v, n):
        return TPS({(0, 0, 0, 0): v if isinstance(v, I) else I(np.full(n, float(v)))})
    def __add__(s, o):
        c = dict(s.c)
        for k, v in o.c.items():
            c[k] = c[k] + v if k in c else v
        return TPS(c)
    def __mul__(s, o):
        if not isinstance(o, TPS):
            return TPS({k: v * o for k, v in s.c.items()})
        c = {}
        for k1, v1 in s.c.items():
            for k2, v2 in o.c.items():
                k = tuple(x + y for x, y in zip(k1, k2))
                if sum(k) > 3: continue
                p = v1 * v2
                c[k] = c[k] + p if k in c else p
        return TPS(c)
    def compose(s, fd):
        """f(s) given fd = [f(s0), f'(s0), f''(s0), f'''(s0)] as I over the base range of s."""
        e0 = (0, 0, 0, 0)
        d = TPS({k: v for k, v in s.c.items() if k != e0})
        out = TPS({e0: fd[0]}); pw = None
        for k in range(1, 4):
            pw = d if pw is None else pw * d
            out = out + pw * (fd[k] * (1.0 / math.factorial(k)) if k < 3 else fd[k] * cI(arb(1) / 6))
        return out
    def deriv(s, alpha):
        """partial derivative of multi-index alpha at the base point (interval)."""
        f = 1
        for a_ in alpha: f *= math.factorial(a_)
        v = s.c.get(tuple(alpha))
        n = len(next(iter(s.c.values())).lo)
        return (v * float(f)) if v is not None else I(np.zeros(n))

def wirt(Phi, ops):
    """complex Wirtinger derivative (re, im intervals) of the real jet Phi.  ops: list of (p, s): p = 0 (u = x0 + i x1)
    or 1 (w = x2 + i x3); s = -1 holomorphic d_p, +1 antiholomorphic dbar_p;  d = (d_a + s i d_b)/2."""
    n = len(ops)
    re = None; im = None
    for choice in itertools.product((0, 1), repeat=n):
        alpha = [0, 0, 0, 0]; coef = 1.0; m = 0
        for (p, s), ch in zip(ops, choice):
            alpha[2 * p + ch] += 1
            if ch: coef *= s; m += 1
        coef *= 0.5 ** n
        # i^m
        ph = m % 4
        v = Phi.deriv(alpha) * coef
        if ph == 2 or ph == 3: v = -v
        if ph % 2 == 0:
            re = v if re is None else re + v
        else:
            im = v if im is None else im + v
    z = I(np.zeros(len(re.lo if re is not None else im.lo)))
    return (re if re is not None else z, im if im is not None else z)

def cmag(c):
    return isqrt(isq(c[0]) + isq(c[1]))

def eh_chart(Ulo, Uhi, Wlo, Whi):
    """EH (a = 1) chart potential Phi = Psi(t) + log(1 + |w|^2), t = |u|^2 (1 + |w|^2)^2, Psi(t) = sqrt(1+t)
    - log(1 + sqrt(1+t)) (pluriharmonic log|u| dropped).  Returns ginv diagonal (g^{u ubar}, g^{w wbar}), the
    g0-norm Gamma of d g0 (chart coordinates) and det g."""
    n = len(Ulo)
    X = [TPS.var(0, I(Ulo, Uhi)), TPS.var(1, I(np.zeros(n))), TPS.var(2, I(Wlo, Whi)), TPS.var(3, I(np.zeros(n)))]
    uu = X[0] * X[0] + X[1] * X[1]
    q = X[2] * X[2] + X[3] * X[3]
    Wt = q + TPS.const(1.0, n)
    t = uu * Wt * Wt
    t0 = t.c[(0, 0, 0, 0)]; t0 = I(np.maximum(t0.lo, 0.0), t0.hi)
    S = isqrt(t0 + 1.0)
    oS = S + 1.0
    P1 = irecip(oS * 2.0)
    P2 = -irecip(S * isq(oS) * 4.0)
    P3 = (irecip(isq(S) * isq(oS)) + irecip(S * oS * isq(oS)) * 2.0) * irecip(S * 8.0)
    Psi = t.compose([I(np.zeros(n)), P1, P2, P3])
    W0 = Wt.c[(0, 0, 0, 0)]
    L = Wt.compose([I(np.zeros(n)), irecip(W0), -irecip(isq(W0)), irecip(isq(W0) * W0) * 2.0])
    Phi = Psi + L
    g = {}
    for i in (0, 1):
        for j in (0, 1):
            g[(i, j)] = wirt(Phi, [(i, -1), (j, +1)])
    det = g[(0, 0)][0] * g[(1, 1)][0] - isq(g[(0, 1)][0]) - isq(g[(0, 1)][1])
    if np.any(det.lo <= 0):
        raise FloatingPointError('EH chart metric: det not certified positive')
    gi = [g[(1, 1)][0] * irecip(det), g[(0, 0)][0] * irecip(det)]
    for v in gi:
        if np.any(v.lo <= 0): raise FloatingPointError('EH chart: g^ii not positive')
    sg = [isqrt(v) for v in gi]
    Gam = None
    for k in (0, 1):
        for i in (0, 1):
            for j in (0, 1):
                d = wirt(Phi, [(k, -1), (i, -1), (j, +1)])
                term = cmag(d) * sg[k] * sg[i] * sg[j]
                Gam = term if Gam is None else Gam + term
    return gi, Gam, det

def chart_entry_bound(H, A, ru, rw):
    """A7: |d^H dbar^A (G o pi)| <= a1! a2! b1! b2! ru^-(a1+b1) rw^-(a2+b2) sup|G|; returns the factor (float)."""
    a1 = H.count(0); a2 = H.count(1); b1 = A.count(0); b2 = A.count(1)
    return (math.factorial(a1) * math.factorial(a2) * math.factorial(b1) * math.factorial(b2)
            * ru ** (-(a1 + b1)) * rw ** (-(a2 + b2)))

RADII = (0.01, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0)
def zone_Z1(SER, a, nU=16, nW=16, scale=1.0, point=None):
    """a: arb ball.  point = (U, W) evaluates the degenerate box at that point (used for the lower bound)."""
    if point is None:
        Ue = np.linspace(0, 1.0 / 16, nU + 1); We = np.linspace(0, 1.0, nW + 1)
        Ulo, Wlo = np.meshgrid(Ue[:-1], We[:-1], indexing='ij'); Uhi, Whi = np.meshgrid(Ue[1:], We[1:], indexing='ij')
        Ulo, Uhi, Wlo, Whi = Ulo.ravel(), Uhi.ravel(), Wlo.ravel(), Whi.ravel()
        keep = Ulo * (1 + Wlo * Wlo) <= 1.0 / 16
        Ulo, Uhi, Wlo, Whi = Ulo[keep], Uhi[keep], Wlo[keep], Whi[keep]
    else:
        Ulo = Uhi = np.array([float(point[0])]); Wlo = Whi = np.array([float(point[1])])
    gi, Gam, det = eh_chart(Ulo, Uhi, Wlo, Whi)
    check_finite('Z1 EH', gi[0], gi[1], Gam, det)
    sg = [np.sqrt(gi[0].hi) * (1 + 1e-15), np.sqrt(gi[1].hi) * (1 + 1e-15)]
    sc = arb(scale)
    a8 = a ** 8
    ahi = a.upper()
    nb = len(Ulo)
    # sup over D_rhoP of |A_ee| for each radius pair (box upper corner)
    res = dict(e=np.zeros(nb), N3=np.zeros(nb), N4=np.zeros(nb))
    combos = {'e': [([i], [j]) for i in (0, 1) for j in (0, 1)],
              'N3': [([k, i], [j]) for k in (0, 1) for i in (0, 1) for j in (0, 1)],
              'N4': [([k, i], [l, j]) for k in (0, 1) for i in (0, 1) for l in (0, 1) for j in (0, 1)]}
    best = {key: np.full(nb, np.inf) for key in combos}
    for ru in RADII:
        for rw in RADII:
            rho2 = 2 * (Uhi + ru) * (1 + (Whi + rw) ** 2)
            rho = np.sqrt(rho2) * (1 + 1e-14)
            # admissibility: a rho < 1/2 - a  (and 2 a rho < 1)
            if np.any(float(ahi) * rho >= 0.5 - float(ahi)) or np.any(2 * float(ahi) * rho >= 0.999):
                continue
            for b in range(nb):
                M = aup(SER.Mee(arb(float(rho[b])), a))
                for key, lst in combos.items():
                    tot = 0.0
                    for (H, A) in lst:
                        fac = chart_entry_bound(H, A, ru, rw)
                        w = 1.0
                        for idx in H + A: w *= sg[idx][b]
                        tot += fac * w
                    v = tot * M * (1 + 1e-12)
                    if v < best[key][b]: best[key][b] = v
    out = dict(e=arb(0), a2D=arb(0), a2bis=arb(0), boxes=nb)
    Gu = Gam.hi
    worst = None
    for b in range(nb):
        e = sc * a8 * arb(best['e'][b]); N3 = sc * a8 * arb(best['N3'][b]); N4 = sc * a8 * arb(best['N4'][b])
        e_u = arb(aup(e))
        D = dR_bound(e_u, arb(aup(N3)), arb(aup(N4)), arb(float(_up(Gu[b]))))
        bb = bis_bound(arb(2), D, e_u)
        for key, v in (('e', e_u), ('a2D', D), ('a2bis', bb)):
            if v.upper() > out[key].upper(): out[key] = arb(aup(v))
        if point is not None:
            out['point'] = dict(e=aup(e_u), a2D=aup(D))
            out['lowerB'] = (2 - D) / (1 + e_u) ** 2
            out['lowerC0'] = C0_lower(D, e_u)
    out['Gamma_max'] = float(np.max(Gu)); out['ginv_max'] = [float(np.max(gi[0].hi)), float(np.max(gi[1].hi))]
    return out

# =====================================================================================================
# 6. neck + bulk: W-chamber cover of V0, Ewald at box centres, order-3 Taylor model, adaptive
# =====================================================================================================
ALPHA = np.array([[0, 1, -1, 0], [0, 0, 1, -1], [0, 0, 0, 1], [0.5, -0.5, -0.5, -0.5]])   # simple roots of W(F4)
MINV = lattice_points(0.5); MINV = MINV[(MINV * MINV).sum(1) > 0]          # 24 minimal vectors, |l| = 1/2

def complex_structures():
    W = WA.enumerate_W(); J = WA.JMAT
    Js = []
    for g in W:
        Jp = g @ J @ g.T
        v = Jp.flatten(); nz = v[np.abs(v) > 1e-12][0]
        if nz < 0: Jp = -Jp
        if not any(np.allclose(Jp, q) for q in Js): Js.append(np.round(Jp * 2) / 2)
    assert len(Js) == 6
    return Js

C2 = list(itertools.combinations_with_replacement(range(4), 2))
C3 = list(itertools.combinations_with_replacement(range(4), 3))
C4 = list(itertools.combinations_with_replacement(range(4), 4))
IDX = {c: i for L_ in (C2, C3, C4) for i, c in enumerate(L_)}
def mult(c):
    from collections import Counter
    m = math.factorial(len(c))
    for v in Counter(c).values(): m //= math.factorial(v)
    return m
def key(*ix): return tuple(sorted(ix))

def kappa_norm(nu, n):
    """elementary bound ||D^n |x|^(-2 nu)|| <= kappa(nu, n) |x|^(-2 nu - n)  (Faa di Bruno + Cauchy-Schwarz)."""
    s = 0
    for k in range(n // 2 + 1):
        poch = 1
        for i in range(n - k): poch *= (nu + i)
        s += math.factorial(n) * poch * 2 ** (n - 2 * k) // (math.factorial(k) * math.factorial(n - 2 * k))
    return s

def ewald_tails(T, K, eta=ETA):
    """per-entry tail bounds (arb) for orders n = 2, 3, 4: real space |x-l| >= T, reciprocal |k| >= K (Lemma A6)."""
    eta_a = arb(eta)
    R = arb(T) - arb(0.5)
    def J(k, Rr, beta):
        if not ((2 * beta * Rr * Rr).lower() > k - 1): raise ValueError('A6 (e) condition fails')
        return Rr ** (k - 1) * (-beta * Rr * Rr).exp() / (2 * beta - max(k - 1, 0) / (Rr * Rr))
    hcoef = {2: {2: 4, 0: 2}, 3: {3: 8, 1: 12}, 4: {4: 16, 2: 48, 0: 12}}
    if not (R > (2 / eta_a).sqrt()): raise ValueError('A6: T - 1/2 <= sqrt(2/eta)')
    real = {}
    for n in (2, 3, 4):
        tot = arb(0)
        for j, anj in hcoef[n].items():
            mj = (n + j) // 2
            for i in range(mj + 1):
                cq = anj * arb(math.factorial(mj)) / math.factorial(i) * eta_a ** i
                q = 2 * i - n - 2
                tot += cq * sum(math.comb(3, jp) * arb(4) ** (jp - 3) * J(q + jp, R, eta_a) for jp in range(4))
        real[n] = 1024 * tot
    Rk = arb(K) - 2 * arb(2).sqrt()
    beta = arb.pi() ** 2 / eta_a
    rec = {}
    for n in (2, 3, 4):
        if not ((Rk * Rk) >= max(n - 2, 0) * eta_a / (2 * arb.pi() ** 2)): raise ValueError('A6 recip monotone')
        if not ((2 * beta * Rk * Rk) > n): raise ValueError('A6 recip condition')
        rec[n] = (2 * arb.pi()) ** n * VINV * sum(math.comb(3, jp) * arb(2).sqrt() ** (3 - jp) * J(n - 2 + jp, Rk, beta)
                                                  for jp in range(4))
    return real, rec

class Ewald:
    """certified D^2, D^3, D^4 of psi_1 = 8 pi^2 |x|^2 - S_reg/2 at exact (dyadic) points |x| <= 0.4."""
    def __init__(self, T=None, K=None, tol=1e-7, eta=ETA):
        self.eta = eta
        if T is None:
            T = 1.0
            while True:
                try:
                    re, _ = ewald_tails(T, 30.0, eta)
                    if max(aup(v) for v in re.values()) <= tol: break
                except ValueError:
                    pass
                T += 0.05
        if K is None:
            K = 5.0
            while True:
                try:
                    _, rc = ewald_tails(T, K, eta)
                    if max(aup(v) for v in rc.values()) <= tol: break
                except ValueError:
                    pass
                K += 0.25
        self.T, self.K = T, K
        re, rc = ewald_tails(T, K, eta)
        self.tail = {n: 0.5 * (aup(re[n]) + aup(rc[n])) for n in (2, 3, 4)}
        Ls = lattice_points(T + 0.45)
        self.L = Ls[(Ls * Ls).sum(1) > 0]
        M = int(K / 2) + 1                                              # L* = 2 D4 = {k in (2Z)^4 : sum k/2 even}
        ev = 2 * np.arange(-M, M + 1)
        ks = np.stack(np.meshgrid(ev, ev, ev, ev, indexing='ij'), -1).reshape(-1, 4).astype(float)
        k2_ = (ks * ks).sum(1)
        ks = ks[(k2_ > 0) & (k2_ < K * K) & (np.mod(ks.sum(1) / 2, 2) == 0)]
        self.k = ks
        k2 = (ks * ks).sum(1).astype(np.int64)
        cache = {}
        lo = np.empty(len(ks)); hi = np.empty(len(ks))
        for i, q in enumerate(k2):
            if q not in cache:
                cache[q] = WA.a2f((-arb.pi() ** 2 * q / arb(eta)).exp() / q)
            lo[i], hi[i] = cache[q]
        self.ck = I(lo, hi)
        self.c16pi2 = cI(16 * arb.pi() ** 2)
        self.tp = {n: cI((2 * arb.pi()) ** n) for n in (2, 3, 4)}
    def Fderivs(self, s):
        """F(s) = e^(-eta s)/s: F^(j)(s) = (-1)^j e^(-eta s) sum_i C(j,i) eta^(j-i) i! s^(-i-1), j = 1..4."""
        ETA = self.eta
        E = iexp_neg(s * (-ETA))
        inv = irecip(s)
        pw = [inv]
        for _ in range(4): pw.append(pw[-1] * inv)
        out = []
        for j in range(1, 5):
            acc = None
            for i in range(j + 1):
                t = pw[i] * float(math.comb(j, i) * ETA ** (j - i) * math.factorial(i))
                acc = t if acc is None else acc + t
            v = E * acc
            out.append(v if j % 2 == 0 else -v)
        return out
    def vderivs_arb(self, s):
        """v(s) = (e^(-eta s) - 1)/s, derivatives j = 1..4 at an exact point s >= 0 (arb)."""
        s = arb(s); eta = arb(self.eta)
        out = []
        for j in range(1, 5):
            if (eta * s) < 1:
                acc = arb(0); M = 60
                for mm in range(j + 1, M + 1):
                    acc += (-eta) ** mm / math.factorial(mm) * math.factorial(mm - 1) / math.factorial(mm - 1 - j) * s ** (mm - 1 - j)
                mm = M + 1
                t = abs(eta ** mm / math.factorial(mm) * math.factorial(mm - 1) / math.factorial(mm - 1 - j) * s ** (mm - 1 - j))
                out.append(acc + arb(0, 2 * t.upper()))
            else:
                F = sum(math.comb(j, i) * eta ** (j - i) * math.factorial(i) * s ** (-i - 1) for i in range(j + 1))
                F = (-eta * s).exp() * F * (-1) ** j
                out.append(F - (-1) ** j * math.factorial(j) * s ** (-j - 1))
        return out
    @staticmethod
    def radial_sums(g, y):
        """g: [g1..g4] I arrays (nc, m); y: float array (nc, m, 4) of exact multiples of 2^-13 with |y| < 4, so that
        monomials of degree <= 3 are exact doubles (|y 2^13| < 2^15, products < 2^45); degree-4 monomials carry one
        rounding (handled in imulp).  Returns D2, D3, D4 (unique components, I arrays (nc,)) of sum_m g(|y_m|^2)."""
        if not (np.all(np.ldexp(y, 13) == np.rint(np.ldexp(y, 13))) and np.all(np.abs(y) < 4)):
            raise ValueError('radial_sums: coordinates not on the 2^-13 grid')
        Y = [y[..., i] for i in range(4)]
        YY = {c: Y[c[0]] * Y[c[1]] for c in C2}
        YYY = {c: YY[c[:2]] * Y[c[2]] for c in C3}
        ones = np.ones(y.shape[:2])
        M2 = np.stack([ones] + Y + [YY[c] for c in C2], -1)                 # exact
        M3 = np.stack([YY[c] for c in C2] + [YYY[c] for c in C3], -1)       # exact
        M4 = np.stack([YYY[c[:3]] * Y[c[3]] for c in C4], -1)               # one rounding each
        S1 = idot(g[0], ones[..., None]); R2_ = idot(g[1], M2); R3_ = idot(g[2], M3)
        R4_ = idot(g[3], M4, inexact=np.ones(len(C4)))
        tk = lambda R, k: I(R.lo[:, k], R.hi[:, k])
        S1 = tk(S1, 0); S2 = tk(R2_, 0)
        S2y = [tk(R2_, 1 + i) for i in range(4)]
        S2yy = {c: tk(R2_, 5 + q) for q, c in enumerate(C2)}
        S3yy = {c: tk(R3_, q) for q, c in enumerate(C2)}
        S3yyy = {c: tk(R3_, 10 + q) for q, c in enumerate(C3)}
        S4y4 = {c: tk(R4_, q) for q, c in enumerate(C4)}
        D2 = {}
        for (i, j) in C2:
            v = S2yy[(i, j)] * 4.0
            if i == j: v = v + S1 * 2.0
            D2[(i, j)] = v
        D3 = {}
        for (i, j, k) in C3:
            v = S3yyy[(i, j, k)] * 8.0
            extra = None
            for (p, q, r) in ((i, j, k), (i, k, j), (j, k, i)):
                if p == q:
                    extra = S2y[r] if extra is None else extra + S2y[r]
            if extra is not None: v = v + extra * 4.0
            D3[(i, j, k)] = v
        D4 = {}
        for c in C4:
            i, j, k, l = c
            v = S4y4[c] * 16.0
            ex = None
            for (p, q, r, s_) in ((i, j, k, l), (i, k, j, l), (i, l, j, k), (j, k, i, l), (j, l, i, k), (k, l, i, j)):
                if p == q:
                    t = S3yy[key(r, s_)]
                    ex = t if ex is None else ex + t
            if ex is not None: v = v + ex * 8.0
            cnt = sum(1 for (p, q, r, s_) in ((i, j, k, l), (i, k, j, l), (i, l, j, k)) if p == q and r == s_)
            if cnt: v = v + S2 * float(4 * cnt)
            D4[c] = v
        return D2, D3, D4
    def centres(self, X, chunk=64):
        """X: (n,4) exact dyadic centres.  Returns D2, D3, D4 of psi_1 (dicts of I arrays (n,))."""
        n = len(X)
        out = [{c: [] for c in C2}, {c: [] for c in C3}, {c: [] for c in C4}]
        T2 = self.T * self.T * (1 + 1e-9)
        for c0 in range(0, n, chunk):
            x = X[c0:c0 + chunk]; nc = len(x)
            # real space l != 0: only l within T of the chunk's bounding box can pass the mask
            blo = x.min(0); bhi = x.max(0)
            gap = np.maximum(0.0, np.maximum(blo[None] - self.L, self.L - bhi[None]))
            Lc = self.L[(gap * gap).sum(1) < T2 * (1 + 1e-9)]
            y = x[:, None, :] - Lc[None]
            s_f = (y * y).sum(2)
            mask = s_f < T2
            s = I(_dn(_dn(s_f)), _up(_up(s_f)))                 # |y|^2: 4 exact squares, three rounded additions
            s = I(np.where(mask, s.lo, 1.0), np.where(mask, s.hi, 1.0))
            g = self.Fderivs(s)
            g = [I(np.where(mask, v.lo, 0.0), np.where(mask, v.hi, 0.0)) for v in g]
            R2, R3, R4 = self.radial_sums(g, y)
            # l = 0 term v(|x|^2)
            gv = [[], [], [], []]
            for xi in x:
                vv = self.vderivs_arb(arb(xi[0]) ** 2 + arb(xi[1]) ** 2 + arb(xi[2]) ** 2 + arb(xi[3]) ** 2)
                for j in range(4): gv[j].append(WA.a2f(vv[j]))
            gv = [I(np.array([p[0] for p in col])[:, None], np.array([p[1] for p in col])[:, None]) for col in gv]
            V2, V3, V4 = self.radial_sums(gv, x[:, None, :])
            # reciprocal
            if not np.all(np.ldexp(x, 13) == np.rint(np.ldexp(x, 13))):
                raise ValueError('centres must lie on the 2^-13 grid (exact k.x and monomials)')
            kx = x @ self.k.T                                   # exact (integers times dyadics)
            t = kx - np.rint(kx)
            cs, sn = cos_sin_2pi(t)
            ckb = I(np.broadcast_to(self.ck.lo, cs.lo.shape), np.broadcast_to(self.ck.hi, cs.lo.shape))
            Cw = cs * ckb; Sw = sn * ckb
            kk = self.k
            K2m = np.stack([kk[:, c[0]] * kk[:, c[1]] for c in C2], -1)             # integer products: exact
            K3m = np.stack([kk[:, c[0]] * kk[:, c[1]] * kk[:, c[2]] for c in C3], -1)
            K4m = np.stack([kk[:, c[0]] * kk[:, c[1]] * kk[:, c[2]] * kk[:, c[3]] for c in C4], -1)
            Q2 = idot(Cw, K2m); Q3 = idot(Sw, K3m); Q4 = idot(Cw, K4m)
            Rc2 = {c: I(Q2.lo[:, q], Q2.hi[:, q]) for q, c in enumerate(C2)}
            Rc3 = {c: I(Q3.lo[:, q], Q3.hi[:, q]) for q, c in enumerate(C3)}
            Rc4 = {c: I(Q4.lo[:, q], Q4.hi[:, q]) for q, c in enumerate(C4)}
            for c in C2:
                v = (R2[c] + V2[c] + Rc2[c] * self.tp[2] * (-VINV)) * (-0.5)
                if c[0] == c[1]: v = v + I(np.full(nc, self.c16pi2.lo), np.full(nc, self.c16pi2.hi))
                t_ = self.tail[2]; out[0][c].append(I(_dn(v.lo - t_), _up(v.hi + t_)))
            for c in C3:
                v = (R3[c] + V3[c] + Rc3[c] * self.tp[3] * VINV) * (-0.5)
                t_ = self.tail[3]; out[1][c].append(I(_dn(v.lo - t_), _up(v.hi + t_)))
            for c in C4:
                v = (R4[c] + V4[c] + Rc4[c] * self.tp[4] * VINV) * (-0.5)
                t_ = self.tail[4]; out[2][c].append(I(_dn(v.lo - t_), _up(v.hi + t_)))
        res = []
        for d in out:
            res.append({c: I(np.concatenate([q.lo for q in v]), np.concatenate([q.hi for q in v])) for c, v in d.items()})
        return res

NEAR_R = 2.0
class NearSums:
    """S_p(box) >= sum_{l != 0} |xi - l|^-p for all xi in the box (near l: distance to the box; far: A4)."""
    def __init__(self, rmax_global=0.42):
        Ls = lattice_points(NEAR_R)
        self.L = Ls[(Ls * Ls).sum(1) > 0]
        self.Tt = NEAR_R - rmax_global
        self.tail = {}
    def tail_p(self, p):
        if p not in self.tail:
            self.tail[p] = aup(A4_tail(arb(self.Tt), p))
        return self.tail[p]
    def __call__(self, lo, hi, ps):
        """lo, hi: (n,4) box corners. returns dict p -> upper bounds (n,), and dmin (n,) lower bounds."""
        n = len(lo)
        res = {p: np.zeros(n) for p in ps}; dmin = np.full(n, np.inf)
        for c0 in range(0, n, 256):
            l_ = lo[c0:c0 + 256][:, None, :]; h_ = hi[c0:c0 + 256][:, None, :]
            L = self.L[None]
            gap = np.maximum(0.0, np.maximum(l_ - L, L - h_))
            gap = _dn(gap)
            d2 = _dn(_dn(gap[..., 0] ** 2 + gap[..., 1] ** 2) + _dn(gap[..., 2] ** 2 + gap[..., 3] ** 2))
            d2 = _dn(d2 * (1 - 4 * U53))
            if np.any(d2 <= 0): raise FloatingPointError('box touches a lattice point')
            dmin[c0:c0 + 256] = np.sqrt(d2.min(1)) * (1 - 4 * U53)
            inv = _up(1.0 / np.sqrt(d2) * (1 + 4 * U53))      # >= 1/dist
            pw = {1: inv}
            for k in range(2, max(ps) + 1):
                pw[k] = _up(pw[k - 1] * inv)                      # repeated rounded-up products
            for p in ps:
                res[p][c0:c0 + 256] = _up(pw[p].sum(1) * (1 + 1.01 * self.L.shape[0] * U53))
        for p in ps:
            res[p] = _up(res[p] + self.tail_p(p))
        return res, dmin

def box_filter(lo, hi):
    """keep boxes that may meet chamber cap V0 (conservative: each half-space test separately)."""
    keep = np.ones(len(lo), bool)
    for a_ in ALPHA:
        mx = np.where(a_ >= 0, hi, lo) @ a_
        keep &= mx >= -1e-12
    for l in MINV:
        mn = np.where(l >= 0, lo, hi) @ l
        keep &= mn <= 0.125 + 1e-12
    return keep

def initial_boxes(h0):
    n = int(round(0.375 / (2 * h0)))
    e = np.arange(n + 1) * 2 * h0
    g = np.stack(np.meshgrid(e[:-1], e[:-1], e[:-1], e[:-1], indexing='ij'), -1).reshape(-1, 4)
    lo = g; hi = g + 2 * h0
    k = box_filter(lo, hi)
    return lo[k], hi[k]

class BulkEval:
    def __init__(self, ew, ns, Js):
        self.ew, self.ns, self.Js = ew, ns, Js
    def evaluate(self, lo, hi):
        """per box: U (upper |K1|_op), Lc (rigorous lower at centre), n3, n4 (upper real norms of D^3, D^4 psi_1),
        S_p for psi~_2, dmin, rmin, rmax."""
        n = len(lo)
        c = 0.5 * (lo + hi); h = 0.5 * (hi[:, 0] - lo[:, 0])
        D2, D3, D4 = self.ew.centres(c)
        S, dmin = self.ns(lo, hi, (7, 8, 9, 10, 12, 13, 14))
        D5 = 0.5 * kappa_norm(1, 5) * S[7] * (1 + 1e-12)        # sup ||D^5 psi_1|| over the box
        dl = 2.0 * h * (1 + 1e-15)                            # |delta| <= 2h (4-D half diagonal)
        rho = _up(D5 * dl ** 3 / 6.0 * (1 + 1e-12))
        H = I(-h, h)
        H2 = I(np.zeros(n), _up(h * h)); Hx = I(-_up(h * h), _up(h * h))
        Tm = {}
        for (i, j) in C2:
            v = D2[(i, j)]
            for k in range(4):
                v = v + D3[key(i, j, k)] * H
            for k in range(4):
                v = v + D4[key(i, j, k, k)] * H2 * 0.5
                for l in range(k + 1, 4):
                    v = v + D4[key(i, j, k, l)] * Hx
            Tm[(i, j)] = v
        def full(Mu, i, j): return Mu[key(i, j)]
        U = np.zeros(n); Lc = np.zeros(n)
        for J in self.Js:
            for which, MM in (('box', Tm), ('ctr', D2)):
                M = {}
                for (i, j) in C2:
                    v = full(MM, i, j)
                    acc = None
                    for p in range(4):
                        if J[p, i] == 0: continue
                        for q in range(4):
                            if J[q, j] == 0: continue
                            t = full(MM, p, q) * float(J[p, i] * J[q, j])
                            acc = t if acc is None else acc + t
                    M[(i, j)] = (v + acc) * 0.5
                tr = M[(0, 0)] + M[(1, 1)] + M[(2, 2)] + M[(3, 3)]
                f2 = None
                for (i, j) in C2:
                    t = isq(M[(i, j)]) * (1.0 if i == j else 2.0)
                    f2 = t if f2 is None else f2 + t
                arg = f2 - isq(tr) * 0.25
                op = iabs(tr) * 0.25 + isqrt(arg) * 0.5         # ||M||_op for M commuting with J
                if which == 'box':
                    U = np.maximum(U, _up(op.hi * 0.5 + rho * 0.5))
                else:
                    Lc = np.maximum(Lc, _dn(op.lo * 0.5))
        # the lower bound only counts centres inside V0 (x.l <= 1/8 for the 24 minimal l; exact for dyadic data)
        inV0 = np.all(c @ MINV.T <= 0.125, axis=1)
        Lc = np.where(inV0, Lc, 0.0)
        # D^3 over box (HS) and D^4
        f3 = None
        for c in C3:
            v = D3[c]
            for k in range(4):
                v = v + D4[key(*c, k)] * H
            t = isq(v) * float(mult(c))
            f3 = t if f3 is None else f3 + t
        n3 = _up(np.sqrt(f3.hi) * (1 + 1e-15) + 0.5 * D5 * dl ** 2 * (1 + 1e-12))
        f4 = None
        for c in C4:
            t = isq(D4[c]) * float(mult(c))
            f4 = t if f4 is None else f4 + t
        n4 = _up(np.sqrt(f4.hi) * (1 + 1e-15) + D5 * dl * (1 + 1e-12))
        # radii
        rmin2 = _dn((np.maximum(0, np.maximum(lo, -hi)) ** 2).sum(1))
        rmax2 = _up((np.maximum(np.abs(lo), np.abs(hi)) ** 2).sum(1) * (1 + 1e-15))
        check_finite('bulk', U, Lc, n3, n4, *[S[p] for p in S], dmin)
        return dict(U=U, Lc=Lc, n3=n3, n4=n4, S=S, dmin=dmin, rmin=np.sqrt(np.maximum(rmin2, 0.0)) * (1 - 1e-15),
                    rmax=_up(np.sqrt(rmax2) * (1 + 1e-15)), h=h)

def bulk_cover(state, tol, h0, hmin, budget_s, log=print, max_boxes=400000):
    """adaptive cover: split boxes whose U exceeds (1+tol) * best lower bound, until none (or h <= hmin).
    state (dict) is updated in place; returns True when finished."""
    t0 = time.time()
    if 'active' not in state:
        lo, hi = initial_boxes(h0)
        state.update(active=(lo, hi), final={k: [] for k in ('lo', 'hi', 'U', 'n3', 'n4', 'S8', 'S9', 'S10', 'S12', 'S13', 'S14', 'dmin', 'rmin', 'rmax')},
                     L=0.0, it=0, evaluated=0)
    BE = state['_BE']
    CH = 3000                                         # boxes per work unit (~30 s): keeps every invocation bounded
    while len(state['active'][0]):
        alo_, ahi_ = state['active']
        if len(alo_) > CH:                            # evaluate a chunk; the rest stays active (order irrelevant)
            lo, hi = alo_[:CH], ahi_[:CH]; rest = (alo_[CH:], ahi_[CH:])
        else:
            lo, hi = alo_, ahi_; rest = (np.zeros((0, 4)), np.zeros((0, 4)))
        r = BE.evaluate(lo, hi)
        state['evaluated'] += len(lo)
        state['L'] = max(state['L'], float(r['Lc'].max()))
        thr = state['L'] * (1 + tol)
        split = (r['U'] > thr) & (r['h'] > hmin)
        fin = ~split
        F = state['final']
        F['lo'].append(lo[fin]); F['hi'].append(hi[fin]); F['U'].append(r['U'][fin]); F['n3'].append(r['n3'][fin])
        F['n4'].append(r['n4'][fin]); F['dmin'].append(r['dmin'][fin]); F['rmin'].append(r['rmin'][fin]); F['rmax'].append(r['rmax'][fin])
        for p in (8, 9, 10, 12, 13, 14): F['S%d' % p].append(r['S'][p][fin])
        lo_s, hi_s = lo[split], hi[split]
        if len(lo_s):
            m = 0.5 * (lo_s + hi_s)
            nlo = []; nhi = []
            for bits in itertools.product((0, 1), repeat=4):
                b = np.array(bits, bool)
                nlo.append(np.where(b, m, lo_s)); nhi.append(np.where(b, hi_s, m))
            nlo = np.concatenate(nlo); nhi = np.concatenate(nhi)
            k = box_filter(nlo, nhi)
            nlo, nhi = nlo[k], nhi[k]
        else:
            nlo = np.zeros((0, 4)); nhi = np.zeros((0, 4))
        state['active'] = (np.concatenate([rest[0], nlo]), np.concatenate([rest[1], nhi])); state['it'] += 1
        log('  bulk unit %d: evaluated %d, split %d -> %d new, active %d, L = %.6f, max U(final so far) = %.6f, %.1f s'
            % (state['it'], len(lo), int(split.sum()), len(nlo), len(state['active'][0]), state['L'],
               max([float(np.max(u)) for u in F['U'] if len(u)] + [0.0]), time.time() - t0))
        if len(state['active'][0]) > max_boxes:
            raise RuntimeError('bulk cover: too many boxes')
        if time.time() - t0 > budget_s and len(state['active'][0]):
            return False
    return True

def bulk_final(state):
    F = state['final']
    return {k: np.concatenate(v) if len(v) else np.zeros(0) for k, v in F.items()}

def zone_NB(fin, Z, a_lo, a_hi, scale=1.0):
    """assembly of the neck/bulk bounds for a in [a_lo, a_hi] (per box, I arithmetic)."""
    n = len(fin['U'])
    alo = float(a_lo); ahi = float(a_hi)
    act = fin['rmax'] >= alo
    if not np.any(act):
        return dict(e=0.0, a2D=0.0, a2bis=0.0, boxes=0)
    f = {k: v[act] for k, v in fin.items()}
    m = len(f['U'])
    rstar = np.maximum(f['rmin'], alo)
    sig = I(_dn(np.maximum(1.0, _dn(_dn(rstar * rstar) / _up(ahi * ahi)))))
    fh = isqrt(isq(sig) + 1.0)
    th = I(_dn((sig * irecip(fh)).lo))
    gam = I(_up((irecip(sig * fh * isqrt(fh)) * cI(arb(3).sqrt())).hi))
    a2K = I(_up((irecip(fh * fh * fh) * 2.0).hi))
    eps = I(np.full(m, _up(ahi ** 4 * (1 + 1e-15))))
    a = I(np.full(m, ahi)); a2 = I(np.full(m, _up(ahi * ahi)))
    Z2 = float(Z[2.0].upper()); Z3 = 16 * Z2 * (1 + 1e-15)
    def psi2(nn, S6n, S10n, dmin):
        """upper bound of ||D^n psi~_2|| (sum over k >= 2, tail geometric)."""
        t2 = kappa_norm(3, nn) / 24.0 * S6n + (Z2 / 8.0 * 2.0 if nn == 2 else 0.0)
        t3 = (kappa_norm(5, nn) / 5.0 * S10n + (2 * Z3 * 2.0 if nn == 2 else 0.0)) / 16.0
        q = eps.hi * np.maximum(3.0 / dmin ** 4, 24.0)
        if np.any(q >= 0.5): raise FloatingPointError('psi~_2 series ratio')
        return _up((t2 + eps.hi * t3 / (1 - q)) * (1 + 1e-12))
    k2 = 0.5 * psi2(2, f['S8'], f['S12'], f['dmin']); m3 = psi2(3, f['S9'], f['S13'], f['dmin']); m4 = psi2(4, f['S10'], f['S14'], f['dmin'])
    sc = float(scale)
    e = (eps * I(f['U']) + isq(eps) * I(k2)) * sc * irecip(th)
    N3 = (eps * I(f['n3']) + isq(eps) * I(m3)) * sc * irecip(th * isqrt(th))
    N4 = (eps * I(f['n4']) + isq(eps) * I(m4)) * sc * irecip(isq(th))
    if np.any(e.hi >= 1): raise FloatingPointError('NB: e >= 1')
    q = irecip(1.0 - I(e.hi))
    a2D = a2 * I(N4.hi) + (a * gam * I(N3.hi) * 2.0 + a2 * isq(I(N3.hi))) * q + isq(gam) * I(e.hi) * q
    bis = (a2K + a2D) * irecip(isq(1.0 - I(e.hi)))
    check_finite('NB', e, a2D, bis)
    return dict(e=float(e.hi.max()), a2D=float(a2D.hi.max()), a2bis=float(bis.hi.max()), boxes=int(m),
                K1max=float(np.max(f['U'])), n3max=float(np.max(f['n3'])), n4max=float(np.max(f['n4'])),
                e_over_eps=float((e.hi / eps.hi).max()))

# =====================================================================================================
# 7. assembly
# =====================================================================================================
def a_upper(alo, ahi):
    """every a-dependent factor of the core/neck bounds is a nonnegative product of a^k and (1 - 2 a rho)^(-m)
    (0 <= a <= a_hi, 2 a_hi rho < 1); for such expressions the upper endpoint of the interval evaluation over
    [a_lo, a_hi] is the evaluation at a_hi (interval product/sum/reciprocal rules for nonnegative intervals), which
    is computed here in arb.  No monotonicity of the bounded quantity itself is assumed."""
    return arb(ahi)

def certify(Z, fin, a_intervals, nT, nU, nW, scale=1.0, log=print):
    SER = Series(Z)
    rows = []
    for (alo, ahi) in a_intervals:
        a = a_upper(alo, ahi)
        t0 = time.time()
        z2 = zone_Z2(SER, a, nT=nT, scale=scale)
        z1 = zone_Z1(SER, a, nU=nU, nW=nW, scale=scale)
        nb = zone_NB(fin, Z, alo, ahi, scale=scale)
        e = max(aup(z1['e']), aup(z2['e']), nb['e'])
        a2D = max(aup(z1['a2D']), aup(z2['a2D']), nb['a2D'])
        bis = max(aup(z1['a2bis']), aup(z2['a2bis']), nb['a2bis'])
        row = dict(a=[alo, ahi], e=e, a2_dR=a2D, a2_Bplus_upper=bis, a2_C0_upper=bis,
                   zones=dict(Z1=dict(e=aup(z1['e']), a2D=aup(z1['a2D']), a2bis=aup(z1['a2bis']), boxes=z1['boxes'],
                                      Gamma_max=z1['Gamma_max'], ginv_max=z1['ginv_max']),
                              Z2=dict(e=aup(z2['e']), a2D=aup(z2['a2D']), a2bis=aup(z2['a2bis']),
                                      worst_piece=max(z2['pieces'], key=lambda r: r['a2bis'])),
                              NB=nb),
                   seconds=time.time() - t0)
        check_finite('assembly', e, a2D, bis)
        rows.append(row)
        log('  a in [%g, %g]: e <= %.4e, a^2||dR|| <= %.4e, a^2 B+ <= %.12f (Z1 %.3e/%.3e, Z2 %.3e/%.3e, NB %.3e/%.3e) %.1fs'
            % (alo, ahi, e, a2D, bis, aup(z1['e']), aup(z1['a2D']), aup(z2['e']), aup(z2['a2D']), nb['e'], nb['a2D'], row['seconds']))
    return rows

def lower_bound_E(Z, a_val=0.03, scale=1.0):
    SER = Series(Z)
    z = zone_Z1(SER, arb(a_val), point=(0.0, 0.0), scale=scale)
    return dict(a=a_val, a2_Bplus_lower=adn(z['lowerB']), a2_C0_lower=adn(z['lowerC0']), point=z['point'])

# =====================================================================================================
# 8. double-precision reference (probe1/k1c.py + k1d.py, read-only) with the C-infinity window
# =====================================================================================================
def st_derivs_double(t):
    """st, st', ..., st'''' at 0 < t < 1 (Lemma C (C2), logistic form)."""
    t = np.asarray(t, float)
    g = 1 / t - 1 / (1 - t)
    with np.errstate(over="ignore", invalid="ignore"):
        u = np.where(g > 0, np.exp(-np.minimum(g, 700)) / (1 + np.exp(-np.minimum(g, 700))), 1 / (1 + np.exp(np.maximum(g, -700))))
    Q1 = u * (1 - u); Q2 = Q1 * (1 - 2 * u); Q3 = Q1 * (1 - 6 * u + 6 * u * u); Q4 = Q1 * (1 - 14 * u + 36 * u * u - 24 * u ** 3)
    ug, ugg, uggg, ugggg = -Q1, Q2, -Q3, Q4
    g1 = -1 / t ** 2 - 1 / (1 - t) ** 2; g2 = 2 / t ** 3 - 2 / (1 - t) ** 3
    g3 = -6 / t ** 4 - 6 / (1 - t) ** 4; g4 = 24 / t ** 5 - 24 / (1 - t) ** 5
    s1 = ug * g1
    s2 = ugg * g1 ** 2 + ug * g2
    s3 = uggg * g1 ** 3 + 3 * ugg * g1 * g2 + ug * g3
    s4 = ugggg * g1 ** 4 + 6 * uggg * g1 ** 2 * g2 + ugg * (3 * g2 ** 2 + 4 * g1 * g3) + ug * g4
    return [u, s1, s2, s3, s4]

def reference_setup():
    sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'probe1'))
    import k1c, k1d
    def eta_cinf(z, a, lam, mode='lam'):
        n = len(z); out = np.zeros((n, k1c.NM))
        if mode == 'zero': return out
        out[:, 0] = 1.0
        if mode == 'one': return out
        r = np.linalg.norm(z, axis=1); r1 = R1W * a; r2 = a
        out[r >= r2, 0] = 0.0
        mid = (r > r1) & (r < r2)
        if mid.any():
            zm = z[mid]
            X = k1c.tvar(zm); S = sum(k1c.tmul(Xi, Xi) for Xi in X); s0 = S[:, 0]
            rd = [s0 ** 0.5, 0.5 * s0 ** -0.5, -0.25 * s0 ** -1.5, 0.375 * s0 ** -2.5, -0.9375 * s0 ** -3.5]
            Rt = k1c.tcompose(rd, S)
            T = (Rt - np.eye(1, k1c.NM)[0] * r1) / (r2 - r1)
            fd = st_derivs_double(T[:, 0])
            out[mid] = -k1c.tcompose(fd, T)
            out[mid, 0] += 1.0
        return out
    k1d.eta_tps = eta_cinf                                   # in-memory patch only; k1d.py is not modified
    return k1c, k1d

def reference_points(a, n_each, rng, k1c):
    pts = {}
    def dirs(m):
        d = rng.normal(size=(m, 4)); return d / np.linalg.norm(d, axis=1)[:, None]
    pts['Z1'] = a * (0.001 + 0.249 * rng.random(n_each))[:, None] * dirs(n_each)
    pts['Z2'] = a * (0.25 + 0.75 * rng.random(n_each))[:, None] * dirs(n_each)
    pts['neck'] = (a + (0.1 - a) * rng.random(n_each))[:, None] * dirs(n_each)
    bl = []
    while sum(len(b) for b in bl) < n_each:
        x = rng.uniform(-0.6, 0.6, size=(4 * n_each, 4))
        z = k1c.to_nearest(x, 'D4')
        z = z[np.linalg.norm(z, axis=1) >= 0.1]
        bl.append(z)
    pts['bulk'] = np.concatenate(bl)[:n_each]
    return pts

def reference_eval(k1c, k1d, z, a, rng):
    model = k1d.get_model(a)
    rest, _ = k1d.rest_jets(z, model, 0.5, a)            # region A (model) for |z| < a, B (lattice sums) beyond
    o = k1c.evaluate_eh(z, a, rest, bis=True, rng=rng, detail=True)
    Rf, Gf = o['Rf'], o['Gf']
    n = len(z)
    Ysp = np.zeros((n, 2, 2), complex); Ysp[:, 0, 0] = 1; Ysp[:, 1, 1] = 1
    mn = -k1c.bis_max(-Rf, Gf, Ysp, rng)[0]
    ev = np.linalg.eigvalsh(Gf - np.eye(2)[None])
    e = np.abs(ev).max(1)
    s = (z * z).sum(1); f = np.sqrt(s * s + a ** 4); A = 2 * a ** 4 / f ** 3
    dR = Rf - A[:, None, None, None, None] * k1c.EH_PATTERN[None]
    W = rng.normal(size=(64, 2, 2)) + 1j * rng.normal(size=(64, 2, 2)); W /= np.linalg.norm(W, axis=2)[..., None]
    W = np.concatenate([W, np.array([[[1, 0], [1, 0]], [[0, 1], [0, 1]], [[1, 0], [0, 1]]], complex)])
    vals = np.abs(np.einsum('nijkl,mi,mj,mk,ml->nm', dR, W[:, 0], np.conj(W[:, 0]), W[:, 1], np.conj(W[:, 1])))
    return dict(max=o['BIS'], min=mn, e=e, dR=vals.max(1), a2K=a * a * A)

# =====================================================================================================
# 9. controls and main
# =====================================================================================================
def logger(path):
    fh = open(path, 'a', encoding='utf-8')
    def log(msg):
        line = time.strftime('%H:%M:%S ') + msg
        print(line, flush=True); fh.write(line + '\n'); fh.flush()
    return log

CTRL = dict(h0=1 / 128, hmin=1 / 1024, tol=0.5, nT=64, nU=16, nW=16, a_int=[(0.0, 0.01), (0.01, 0.02), (0.02, 0.03)])
MAIN = dict(h0=1 / 128, hmin=1 / 4096, tol=0.1, nT=256, nU=32, nW=32,
            a_int=[(0.0, 0.005), (0.005, 0.01), (0.01, 0.015), (0.015, 0.02), (0.02, 0.025), (0.025, 0.03)])

def build_bulk(cfg, log, state=None, budget_s=1e9):
    if state is None: state = {}
    ew = Ewald(tol=1e-5); ns = NearSums(); Js = complex_structures()
    state['_BE'] = BulkEval(ew, ns, Js)
    log('  Ewald: eta = %g, T = %.2f (%d real terms), K = %.2f (%d reciprocal terms), tails %s'
        % (ETA, ew.T, len(ew.L), ew.K, len(ew.k), {n: '%.1e' % v for n, v in ew.tail.items()}))
    done = bulk_cover(state, cfg['tol'], cfg['h0'], cfg['hmin'], budget_s, log=log)
    return state, done, ew

def controls(out_path):
    T0 = time.time()
    log = logger(os.path.join(HERE, 'controls_wp2b.log'))
    log('=== wp2b controls ===')
    res = dict(config=CTRL, eta=ETA)
    # --- CHECK: Z_kappa
    Z = Z_enclosures()
    res['Z'] = {str(k): Z[k].str(12) for k in (2.0, 2.5, 3.0)}
    res['Z_closed_form_check'] = {str(k): bool(Z[k].overlaps(Z_closed(k))) for k in (2.0, 2.5, 3.0)}
    hc = Z['shell_counts_to_100']
    def sod(N): return sum(d for d in range(1, N + 1) if N % d == 0 and d % 2 == 1)
    res['shell_count_check_N_le_100'] = all(hc.get(16 * N // 4 if False else 4 * N, 0) == 24 * sod(N) for N in range(1, 26))
    log('  Z: %s  closed-form overlap %s' % (res['Z'], res['Z_closed_form_check']))
    # --- CHECK: chamber and complex structures
    W = np.array(WA.enumerate_W()); rng = np.random.default_rng(5)
    X = rng.normal(size=(3000, 4)); cnt = set()
    for x in X:
        cnt.add(int(np.sum(np.all((W @ x) @ ALPHA.T > 1e-12, axis=1))))
    res['chamber_images_per_point'] = sorted(cnt)
    res['n_complex_structures'] = len(complex_structures())
    # --- CHECK: EH Gamma formula and chart quantities
    sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'probe1')); import k1c
    s = np.array([0.0625, 0.25, 1.0, 4.0])
    J0, lam, _ = k1c.eh_w_jets(s, 1.0)
    dg = J0[1]
    hs = np.sqrt(sum(abs(dg[:, i, j, k] * lam[:, k] * lam[:, i] * lam[:, j]) ** 2 for i in range(2) for j in range(2) for k in range(2)))
    res['EH_Gamma_HS_vs_formula'] = [float(x) for x in hs / (math.sqrt(3) / (s * (s * s + 1) ** 0.75))]
    # --- CHECK: st derivative double formula vs arb
    tt = [0.1, 0.5, 0.83, 0.95]
    dd = st_derivs_double(np.array(tt))
    worst = 0.0
    for i, t in enumerate(tt):
        S = WA.st_of_series(arb_series([arb(t), 1], prec=5), 'cinf').coeffs()
        for k in range(5):
            ref = float((S[k] * math.factorial(k)).mid())
            worst = max(worst, abs(dd[k][i] - ref) / max(1.0, abs(ref)))
    res['st_double_vs_arb_maxrel'] = worst
    # --- bulk cover (a-independent)
    log('  bulk cover (control resolution: tol %.2f, h0 %.4g, hmin %.4g)' % (CTRL['tol'], CTRL['h0'], CTRL['hmin']))
    state, done, ew = build_bulk(CTRL, log)
    fin = bulk_final(state)
    res['bulk'] = dict(boxes_final=int(len(fin['U'])), evaluated=state['evaluated'], K1_upper=float(fin['U'].max()),
                       K1_lower_centre=state['L'], n3_max=float(fin['n3'].max()), n4_max=float(fin['n4'].max()),
                       ewald=dict(T=ew.T, K=ew.K, n_real=len(ew.L), n_recip=len(ew.k), tails=ew.tail))
    # --- CHECK: Ewald centre enclosures vs k1c double-precision lattice jets (a tiny: psi ~ eps psi_1)
    pts = np.array([[0.125, 0.125, 0.125, 0.0625], [0.25, 0.0, 0.0, 0.0], [0.1875, 0.0625, 0.0625, 0.0], [0.0625, 0.03125, 0.0, 0.0]])
    D2, D3, D4 = ew.centres(pts)
    aa = 1e-3
    j2, j3, j4 = k1c.jets_B(pts, aa, 'D4', k1c.mu_delta('D4', aa)[0])
    bad = 0; worst = 0.0
    for b in range(len(pts)):
        for c in C2:
            ref = j2[b][c] / aa ** 4; v = D2[c]
            worst = max(worst, abs(ref - 0.5 * (v.lo[b] + v.hi[b])) / max(1.0, abs(ref)))
            bad += not (v.lo[b] - 1e-6 * max(1, abs(ref)) <= ref <= v.hi[b] + 1e-6 * max(1, abs(ref)))
        for c in C4:
            ref = j4[b][c] / aa ** 4; v = D4[c]
            worst = max(worst, abs(ref - 0.5 * (v.lo[b] + v.hi[b])) / max(1.0, abs(ref)))
            bad += not (v.lo[b] - 1e-6 * max(1, abs(ref)) <= ref <= v.hi[b] + 1e-6 * max(1, abs(ref)))
    res['ewald_vs_k1c'] = dict(maxrel=worst, outside_1e6=bad)
    log('  Ewald vs k1c: max rel %.2e, outside %d' % (worst, bad))
    # --- certified bounds (control resolution) for scale 1, 0, 10
    cert = {}
    for name, scl in (('base', 1.0), ('K1_lattice_off', 0.0), ('K3_lattice_x10', 10.0)):
        log('  certify scale %g' % scl)
        cert[name] = certify(Z, fin, CTRL['a_int'], CTRL['nT'], CTRL['nU'], CTRL['nW'], scale=scl, log=log)
    res['certified'] = cert
    lb = {name: lower_bound_E(Z, 0.03, scale=scl) for name, scl in (('base', 1.0), ('K1_lattice_off', 0.0))}
    res['lower_bound_a003'] = lb
    # K1
    k1 = cert['K1_lattice_off']
    res['K1'] = dict(e_all_zero=all(r['e'] <= 1e-300 for r in k1),
                     Bplus_interval_a003=[lb['K1_lattice_off']['a2_Bplus_lower'], k1[-1]['a2_Bplus_upper']],
                     C0_interval_a003=[lb['K1_lattice_off']['a2_C0_lower'], k1[-1]['a2_C0_upper']])
    res['K1']['PASS'] = bool(res['K1']['e_all_zero'] and res['K1']['Bplus_interval_a003'][0] <= 2 <= res['K1']['Bplus_interval_a003'][1]
                             and res['K1']['C0_interval_a003'][0] <= 2 <= res['K1']['C0_interval_a003'][1])
    log('  K1: %s' % res['K1'])
    # K3
    base, x10 = cert['base'], cert['K3_lattice_x10']
    ratios = [dict(a=b['a'], e=x['e'] / b['e'], a2dR=x['a2_dR'] / b['a2_dR']) for b, x in zip(base, x10)]
    res['K3'] = dict(ratios=ratios, PASS=all(r['e'] >= 5 and r['a2dR'] >= 5 for r in ratios))
    log('  K3: %s' % res['K3'])
    # K2
    k1c_, k1d_ = reference_setup()
    rngr = np.random.default_rng(2026)
    k2 = {}
    allpass = True
    for aval, idx in ((0.01, 0), (0.03, 2)):
        row = base[idx]
        assert row['a'][0] <= aval <= row['a'][1]
        pts = reference_points(aval, 500, rngr, k1c_)
        zt = row['zones']
        zone_of = {'Z1': zt['Z1'], 'Z2': zt['Z2'], 'neck': zt['NB'], 'bulk': zt['NB']}
        info = {}
        for reg, z in pts.items():
            r = reference_eval(k1c_, k1d_, z, aval, rngr)
            zb = zone_of[reg]
            bound = (r['a2K'] + zb['a2D']) / (1 - zb['e']) ** 2
            okmax = np.all(aval ** 2 * r['max'] <= bound)
            okmin = np.all(aval ** 2 * r['min'] >= -bound)
            oke = np.all(r['e'] <= zb['e']); okd = np.all(aval ** 2 * r['dR'] <= zb['a2D'])
            nan = not (np.all(np.isfinite(r['max'])) and np.all(np.isfinite(r['min'])))
            info[reg] = dict(n=len(z), max_a2bis=float(aval ** 2 * r['max'].max()), min_a2bis=float(aval ** 2 * r['min'].min()),
                             max_slack_upper=float((bound - aval ** 2 * r['max']).min()),
                             max_e_sample=float(r['e'].max()), e_cert=zb['e'], max_a2dR_sample=float(aval ** 2 * r['dR'].max()),
                             a2dR_cert=zb['a2D'], ok=bool(okmax and okmin and oke and okd and not nan))
            allpass &= info[reg]['ok']
            log('  K2 a=%g %s: %s' % (aval, reg, info[reg]))
        k2[str(aval)] = info
    res['K2'] = dict(regions=k2, PASS=bool(allpass))
    # NaN / inf guards
    g = {}
    def expect_raise(fn):
        try: fn(); return False
        except (FloatingPointError, ValueError): return True
    g['nan_interval'] = expect_raise(lambda: check_finite('t', I(np.array([np.nan]), np.array([1.0]))))
    g['inf_interval'] = expect_raise(lambda: check_finite('t', I(np.array([0.0]), np.array([np.inf]))))
    g['nonfinite_arb'] = expect_raise(lambda: aup(arb('nan')))
    g['exp_positive_arg'] = expect_raise(lambda: iexp_neg(I(np.array([0.1]), np.array([0.2]))))
    g['e_ge_1'] = expect_raise(lambda: dR_bound(arb(1.5), arb(1), arb(1), arb(1)))
    g['series_divergent'] = expect_raise(lambda: Series(Z).MA(arb(20.0), arb(0.03), 6))
    g['nan_bulk_record'] = expect_raise(lambda: zone_NB({k: (v.copy() if k != 'U' else np.full_like(v, np.nan)) for k, v in fin.items()}, Z, 0.0, 0.01))
    g['nan_box'] = expect_raise(lambda: state['_BE'].evaluate(np.array([[np.nan, 0, 0, 0.0]]), np.array([[0.01, 0.01, 0.01, 0.01]])))
    res['guards'] = dict(g, PASS=all(g.values()))
    # interval exp/cos sanity vs arb
    xs = np.array([-0.0, -1e-3, -0.7, -5.3, -41.0, -66.6])
    E = iexp_neg(I(xs)); okx = all(E.lo[i] <= float(arb(xs[i]).exp().mid()) <= E.hi[i] for i in range(len(xs)))
    ts = np.array([-0.5, -0.3125, 0.0, 0.123046875, 0.375, 0.5])
    cs, sn = cos_sin_2pi(ts)
    okc = all(cs.lo[i] <= float((2 * arb.pi() * ts[i]).cos().mid()) <= cs.hi[i] and sn.lo[i] <= float((2 * arb.pi() * ts[i]).sin().mid()) <= sn.hi[i] for i in range(len(ts)))
    res['interval_exp_cos_check'] = dict(exp=bool(okx), cos_sin=bool(okc), exp_width_rel=float(np.max((E.hi - E.lo) / np.maximum(E.hi, 1e-300))))
    res['ALL_CONTROLS_PASS'] = bool(res['K1']['PASS'] and res['K2']['PASS'] and res['K3']['PASS'] and res['guards']['PASS']
                                    and all(res['Z_closed_form_check'].values()) and res['chamber_images_per_point'] == [1]
                                    and res['n_complex_structures'] == 6 and res['ewald_vs_k1c']['outside_1e6'] == 0
                                    and okx and okc)
    res['seconds'] = time.time() - T0
    json.dump(res, open(out_path, 'w'), indent=1, default=float)
    log('controls done in %.1f s, ALL_CONTROLS_PASS = %s' % (res['seconds'], res['ALL_CONTROLS_PASS']))
    return res

def main(out_path, budget_s=520):
    if not os.path.exists(GO_FILE):
        print('refusing: %s does not exist (the lead creates it)' % GO_FILE); sys.exit(2)
    T0 = time.time()
    os.makedirs(STATE_DIR, exist_ok=True)
    log = logger(os.path.join(HERE, 'run_wp2b.log'))
    sp = os.path.join(STATE_DIR, 'bulk_state.pkl')
    state = {}
    if os.path.exists(sp):
        with open(sp, 'rb') as fh: state = pickle.load(fh)
        log('resuming bulk cover at iteration %d' % state.get('it', 0))
    if not state.get('done'):
        state, done, ew = build_bulk(MAIN, log, state, budget_s=budget_s)
        state.pop('_BE', None)
        state['done'] = done
        with open(sp + '.tmp', 'wb') as fh: pickle.dump(state, fh)
        os.replace(sp + '.tmp', sp)
        if not done:
            log('bulk cover incomplete after %.0f s: exit 3 (rerun the same command)' % (time.time() - T0)); sys.exit(3)
        state['ewald'] = dict(T=ew.T, K=ew.K, tails=ew.tail)
        with open(sp, 'wb') as fh: pickle.dump(state, fh)
        if time.time() - T0 > budget_s * 0.6:
            log('bulk complete; assembly deferred to the next invocation: exit 3'); sys.exit(3)
    Z = Z_enclosures()
    fin = bulk_final(state)
    rows = certify(Z, fin, MAIN['a_int'], MAIN['nT'], MAIN['nU'], MAIN['nW'], log=log)
    lb = lower_bound_E(Z, 0.03)
    out = dict(status='COMPLETE', config=MAIN, eta=ETA, Z={str(k): Z[k].str(15) for k in (2.0, 2.5, 3.0)},
               bulk=dict(boxes_final=int(len(fin['U'])), evaluated=state['evaluated'], K1_upper=float(fin['U'].max()),
                         K1_lower_centre=state['L'], ewald=state.get('ewald')),
               rows=rows, lower_bound_a003=lb, seconds_last_invocation=time.time() - T0)
    json.dump(out, open(out_path, 'w'), indent=1, default=float)
    log('main COMPLETE: %s' % json.dumps(dict(rows=[dict(a=r['a'], e=r['e'], a2dR=r['a2_dR'], a2B=r['a2_Bplus_upper']) for r in rows], lb=lb)))

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--controls'); ap.add_argument('--main'); ap.add_argument('--budget', type=float, default=520)
    args = ap.parse_args()
    if args.controls: controls(args.controls)
    elif args.main: main(args.main, args.budget)
    else: ap.print_help()
