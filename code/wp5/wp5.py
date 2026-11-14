#!/usr/bin/env python
"""wp5.py -- WP5 certified assembly of the Route 2 analysis chain for the K3 Kummer approximant omega_a.

Sources: ANALYSIS-CHAIN.md (Theorem 2, Prop C0, Lemma G, continuity closure), DEFINITIONS-AND-LEMMAS-v2.md (L3, L4,
L5, 4.1-4.3) with the WP2a window c(sigma) = st((sqrt(sigma) - 1/4)/(3/4)), certified inputs of WP2a, WP2b, WP3.

Stage A (new certificate): for a-intervals covering (0, 0.03] (a enters only as an arb ball, one interval touches 0),
  C_core_F >= sup_{|zeta| <= 1} |F| / a^8   and   C_core_R >= sup_{|zeta| <= 1} |Ric(omega_a)|_g / a^6,
  with every higher-order remainder bounded explicitly (BUILD-WP5.md section 2).
Stage B (assembly): Theorem 2 of ANALYSIS-CHAIN in scaled arb arithmetic, smallness of the continuity method, the
  symmetric bound d2(a) = M2 e^eta - 1, the leading coefficient lim d2/a^3 and a0.

Modes
  --controls OUT.json   X1-X4, NaN/inf guards, reduced Stage A/B; records the sha256 of this file.
  --main OUT.json       refuses unless <lab>/wp5/GO-WP5 exists; resumable (state_wp5/), exit 3 = rerun.
Run:  python <lab-tools>/memcap.py 8 py -3.11 <lab>/wp5/wp5.py --controls <lab>/wp5/controls_wp5.json
"""
import os, sys, json, time, math, argparse, hashlib, re
sys.dont_write_bytecode = True                      # read-only imports of sealed files must not write bytecode
import numpy as np
from flint import arb, ctx

K3 = os.environ.get('K3_ROOT') or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(K3, 'wp2'))
import wp2b as WB                                   # sealed, read-only: Z_kappa, lattice points, A4, step derivatives
ctx.prec = 160                                      # >= 128 bits everywhere
GO_FILE = os.path.join(HERE, 'GO-WP5')
STATE_DIR = os.path.join(HERE, 'state_wp5')
A0_STR = '0.03'                                     # L5 holds for a <= 0.03

class Refuse(Exception):
    """raised whenever a hypothesis of a lemma cannot be certified (the code refuses to give a number)."""

# =====================================================================================================
# 0. arb helpers
# =====================================================================================================
def fin(x, name='value'):
    if isinstance(x, arb):
        if not x.is_finite():
            raise FloatingPointError('%s: non-finite arb %s' % (name, x))
    else:
        if not np.all(np.isfinite(np.asarray(x, float))):
            raise FloatingPointError('%s: non-finite %s' % (name, x))
    return x

def up(x):
    """float >= every point of the arb ball x."""
    fin(x)
    v = float(x.upper())
    return float(np.nextafter(v, np.inf)) if v != 0 else 0.0

def dn(x):
    fin(x)
    v = float(x.lower())
    return float(np.nextafter(v, -np.inf)) if v != 0 else 0.0

def U(x):
    """the arb point upper(x) (an exact number >= the ball)."""
    fin(x)
    return arb(x.upper())

def hull(*xs):
    out = None
    for v in xs:
        v = v if isinstance(v, arb) else arb(v)
        out = v if out is None else out.union(v)
    return out

def nonneg(x):
    """x is known (mathematically) to be >= 0: drop the negative part of the ball."""
    fin(x)
    if x.lower() < 0:
        return arb(0).union(arb(x.upper()))
    return x

def ssqrt(x):
    """sqrt of a quantity known to be >= 0 (ball lower ends below 0 come only from midpoint-radius rounding)."""
    fin(x)
    if x.lower() >= 0:
        return x.sqrt()
    return arb(0).union(arb(x.upper()).sqrt())

class AInt:
    """a in [lo, hi], 0 <= lo < hi, endpoints given as decimal strings (arb balls containing the decimals).
    Powers are hulls of the endpoint powers (x -> x^k is increasing on x >= 0)."""
    def __init__(self, lo, hi):
        self.slo, self.shi = str(lo), str(hi)
        self.lo = arb(self.slo); self.hi = arb(self.shi)
        if not (self.lo.lower() >= 0 or self.slo in ('0', '0.0')):
            raise Refuse('a interval must have lo >= 0')
        if self.slo in ('0', '0.0'):
            self.lo = arb(0)
        self._p = {}
    def p(self, k):
        if k not in self._p:
            self._p[k] = hull(self.lo ** k, self.hi ** k) if k > 0 else arb(1)
        return self._p[k]
    @property
    def ball(self):
        return self.p(1)
    def key(self):
        return '[%s, %s]' % (self.slo, self.shi)

def frac(s):
    from fractions import Fraction
    return Fraction(s)

# monotone special functions on nonnegative balls (exact functions, monotonicity is a property of the function)
def _mono(f, z, at0):
    """f monotone on [0, inf): enclosure from the two endpoint values (no dependency blow-up on wide balls)."""
    z = nonneg(z)
    lo = arb(z.lower()) if z.lower() > 0 else arb(0)
    hi = arb(z.upper())
    vlo = f(lo) if lo > 0 else arb(at0)
    vhi = f(hi) if hi > 0 else arb(at0)
    return hull(vlo, vhi)

def phi1(z):   # (e^z - 1)/z, increasing, = 1 at 0
    return _mono(lambda t: t.expm1() / t, z, 1)

def g1(z):     # (1 - e^-z)/z, decreasing, = 1 at 0
    return _mono(lambda t: -(-t).expm1() / t, z, 1)

def ell(q):    # -log(1-q)/q, increasing, = 1 at 0 (needs q < 1)
    if not (q.upper() < 1):
        raise Refuse('C_q a^8 >= 1')
    return _mono(lambda t: -(-t).log1p() / t, q, 1)

# =====================================================================================================
# 1. inputs from the certified result files (echoed with source and seal ids)
# =====================================================================================================
SRC = {
    'wp2a': (K3 + '/wp2/result_wp2a.json', K3 + '/wp2/VERDICT-WP2a.md'),
    'wp2b': (K3 + '/wp2/result_wp2b.json', K3 + '/wp2/VERDICT-WP2b.md'),
    'wp3': (K3 + '/wp3/result_wp3.json', K3 + '/wp3/VERDICT-WP3.md'),
}

def seal_ids(path):
    t = open(path, encoding='utf-8').read()
    first = t.splitlines()[0]
    return dict(verdict_file=path, claim=re.findall(r'claim ([0-9a-f]{10})', first),
                sealed=re.findall(r'sealed ([0-9a-f]{10})', first), header=first.strip())

def sha_file(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()

REL_MARGIN = arb('1e-12')                           # applied to values read back from decimal strings / floats

def load_inputs():
    out = {}
    wa = json.load(open(SRC['wp2a'][0]))
    c6 = hull(arb(wa['results']['c6']['lower']), arb(wa['results']['c6']['upper']))
    P_up = arb(wa['results']['P_cinf']['upper_float']) * (1 + REL_MARGIN)
    R1_up = arb(wa['results']['R_cinf']['upper_c6=1']) * (1 + REL_MARGIN)
    if wa['windows']['cinf'] != [0.25, 1.0] or wa['status'] != 'COMPLETE':
        raise Refuse('WP2a result is not the C-infinity [0.25, 1] window certificate')
    R6 = U(R1_up * c6)                                  # sup over the core of |ddbar F8|_{gh0}
    F8 = U(3 * c6 / 16 * P_up)                          # sup over the core of |F8| = (c6/2)(3/8) P*
    out['wp2a'] = dict(source=SRC['wp2a'][0], seal=seal_ids(SRC['wp2a'][1]), sha256=sha_file(SRC['wp2a'][0]),
                       c6=c6, P_up=P_up, R1_up=R1_up, R6star=R6, supF8=F8,
                       echo=dict(c6=[wa['results']['c6']['lower'], wa['results']['c6']['upper']],
                                 P_cinf_upper=wa['results']['P_cinf']['upper_float'],
                                 R_cinf_upper_c6_1=wa['results']['R_cinf']['upper_c6=1'],
                                 R6star_used=up(R6), supF8_used=up(F8)))
    wb = json.load(open(SRC['wp2b'][0]))
    if wb['status'] != 'COMPLETE':
        raise Refuse('WP2b incomplete')
    rows = [dict(a=(str(r['a'][0]), str(r['a'][1])), e=r['e'], c0s=r['a2_C0_upper'], bs=r['a2_Bplus_upper'])
            for r in wb['rows']]
    out['wp2b'] = dict(source=SRC['wp2b'][0], seal=seal_ids(SRC['wp2b'][1]), sha256=sha_file(SRC['wp2b'][0]),
                       rows=rows)
    w3 = json.load(open(SRC['wp3'][0]))
    if w3['status'] != 'COMPLETE':
        raise Refuse('WP3 incomplete')
    pcs = [dict(a=(str(p['a_interval'][0]), str(p['a_interval'][1])), C_out=p['C_out_upper'], R_out=p['R_out_upper'],
                e_out=p['e_upper']) for p in w3['results']['pieces']]
    out['wp3'] = dict(source=SRC['wp3'][0], seal=seal_ids(SRC['wp3'][1]), sha256=sha_file(SRC['wp3'][0]),
                      pieces=pcs, e_out=w3['results']['e_out_upper'])
    return out

def containing(rows, ai):
    """the input row whose a-interval contains [ai.lo, ai.hi] (decimal comparison)."""
    lo, hi = frac(ai.slo), frac(ai.shi)
    for r in rows:
        if frac(r['a'][0]) <= lo and hi <= frac(r['a'][1]):
            return r
    raise Refuse('no certified input row contains a in %s' % ai.key())

# =====================================================================================================
# 2. exact polynomials (computer algebra, rational arithmetic) and Bombieri norms
# =====================================================================================================
def exact_polys():
    """H6 (doc L3 (e)); the W(F4)-invariant harmonic octic H8n with H8n(e0) = 1 and the PROOF that such octics form a
    line (so P_8[psi_1] = c8 H8n, c8 = P_8[psi_1](e0)); N[H8n]; the proof that W(F4)-invariant harmonic polynomials
    of degree 10 vanish (so P_10[psi_1] = 0); Bombieri norms.  Invariance is imposed for the signed permutations
    (spanning set: polynomials in p2, m2, m3, m4, m5) and for right multiplication by w = (1+i+j+k)/2; both are in
    W = Aut(L)."""
    import sympy as sp
    x = sp.symbols('x0:4', real=True)
    p2 = sum(v ** 2 for v in x)
    m = lambda k: sum(v ** (2 * k) for v in x)
    h = sp.Rational(1, 2)
    def qmul(a, b):
        a0, a1, a2, a3 = a; b0, b1, b2, b3 = b
        return (a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3, a0 * b1 + a1 * b0 + a2 * b3 - a3 * b2,
                a0 * b2 - a1 * b3 + a2 * b0 + a3 * b1, a0 * b3 + a1 * b2 - a2 * b1 + a3 * b0)
    xw = qmul(x, (h, h, h, h))
    # w is orthogonal: |x w| = |x|
    assert sp.expand(sum(v ** 2 for v in xw) - p2) == 0
    def inv_harm_space(basis):
        cs = sp.symbols('c0:%d' % len(basis))
        P = sum(c * b for c, b in zip(cs, basis))
        lap = sp.expand(sum(sp.diff(P, v, 2) for v in x))
        Pw = sp.expand(P.subs(dict(zip(x, xw)), simultaneous=True) - P)
        eqs = []
        for e in (lap, Pw):
            if e != 0:
                eqs += sp.Poly(e, *x).coeffs()
        sol = list(sp.linsolve(eqs, cs))[0]
        free = sorted({s for e in sol for s in e.free_symbols}, key=str)
        polys = []
        for f in free:
            sub = {g: (1 if g == f else 0) for g in free}
            polys.append(sp.expand(sum(sp.sympify(c).subs(sub) * b for c, b in zip(sol, basis))))
        polys = [q for q in polys if q != 0]
        return polys
    B8 = [p2 ** 4, p2 ** 2 * m(2), m(2) ** 2, p2 * m(3), m(4)]
    sp8 = inv_harm_space(B8)
    # the space is spanned by sp8; check its dimension via the rank of the coefficient vectors
    M8 = sp.Matrix([[sp.Poly(q, *x).coeff_monomial(mon) for mon in sp.Poly(sum(sp8), *x).monoms()] for q in sp8])
    dim8 = M8.rank()
    if dim8 != 1:
        raise Refuse('W(F4)-invariant harmonic octics: dimension %d != 1' % dim8)
    H = sp8[0]
    e0 = {x[0]: 1, x[1]: 0, x[2]: 0, x[3]: 0}
    H8n = sp.expand(H / H.subs(e0))
    B10 = [p2 ** 5, p2 ** 3 * m(2), p2 * m(2) ** 2, p2 ** 2 * m(3), m(2) * m(3), p2 * m(4), m(5)]
    sp10 = inv_harm_space(B10)
    p10_zero = all(sp.expand(q) == 0 for q in sp10)
    if not p10_zero:
        raise Refuse('W(F4)-invariant harmonic decics do not vanish')
    H6 = sp.expand(6 * m(3) - sp.Rational(15, 2) * p2 * m(2) + sp.Rational(15, 8) * p2 ** 3)
    # W-invariant sextics: B4-invariant (p2^3, p2 m2, m3) and w-invariant -> a plane containing p2^3 and H6
    cs6 = sp.symbols('d0:3')
    P6g = cs6[0] * p2 ** 3 + cs6[1] * p2 * m(2) + cs6[2] * m(3)
    e6 = sp.expand(P6g.subs(dict(zip(x, xw)), simultaneous=True) - P6g)
    sol6 = list(sp.linsolve(sp.Poly(e6, *x).coeffs(), cs6))[0]
    free6 = sorted({s for e in sol6 for s in e.free_symbols}, key=str)
    sp6 = [sp.expand(sum(sp.sympify(c).subs({g_: (1 if g_ == f else 0) for g_ in free6}) * b
                         for c, b in zip(sol6, [p2 ** 3, p2 * m(2), m(3)]))) for f in free6]
    monos6 = sp.Poly(sum(sp6) + p2 ** 3 + H6, *x).monoms()
    rk = lambda L_: sp.Matrix([[sp.Poly(q, *x).coeff_monomial(mo) for mo in monos6] for q in L_]).rank()
    if not (rk(sp6) == 2 and rk(sp6 + [p2 ** 3, H6]) == 2 and rk([p2 ** 3, H6]) == 2):
        raise Refuse('W-invariant sextics are not span(p2^3, H6)')
    H6_u = H6.subs({x[0]: 1, x[1]: 1, x[2]: 0, x[3]: 0}) / 8                 # (1,1,0,0)/sqrt2
    H6_v = H6.subs({x[0]: 2, x[1]: 1, x[2]: 0, x[3]: 0}) / 125               # (2,1,0,0)/sqrt5
    assert H6.subs({x[0]: 1, x[1]: 0, x[2]: 0, x[3]: 0}) == sp.Rational(3, 8) and H6_u == -sp.Rational(3, 8)
    J = [-x[1], x[0], -x[3], x[2]]
    grad = lambda f: [sp.diff(f, v) for v in x]
    Eu = lambda f: sp.expand(sum(a * b for a, b in zip(x, grad(f))))
    DJ = lambda f: sp.expand(sum(a * b for a, b in zip(J, grad(f))))
    Nop = lambda f: sp.expand((Eu(Eu(f)) + DJ(DJ(f))) / 4)
    NH8 = Nop(H8n)
    assert sp.expand(sum(sp.diff(H8n, v, 2) for v in x)) == 0           # harmonic: Lap H8n = 0
    assert sp.expand(sum(sp.diff(H6, v, 2) for v in x)) == 0
    assert sp.simplify(Nop(p2 ** 2) - 4 * p2 ** 2) == 0                 # N on bidegree (2,2) is 4
    def bomb2(P):
        Pp = sp.Poly(P, *x); n = Pp.total_degree(); s = sp.Integer(0)
        for mon, c in Pp.terms():
            s += c ** 2 * sp.prod([sp.factorial(k) for k in mon]) / sp.factorial(n)
        return sp.nsimplify(s)
    def to_arb_sqrt(q):
        q = sp.Rational(q)
        return (arb(int(q.p)) / arb(int(q.q))).sqrt()
    out = dict(H8n=str(H8n), dim8=int(dim8), p10_zero=bool(p10_zero),
               BH8=to_arb_sqrt(bomb2(H8n)), BN8=to_arb_sqrt(bomb2(NH8)), BH6=to_arb_sqrt(bomb2(H6)),
               H8n_e0=str(H8n.subs(e0)), H6_e0=str(H6.subs(e0)))
    out['_lamb'] = dict(H8n=sp.lambdify(x, H8n, 'numpy'), NH8=sp.lambdify(x, NH8, 'numpy'),
                        H6=sp.lambdify(x, H6, 'numpy'))
    out['_sym'] = dict(x=x, H8n=H8n, NH8=NH8, H6=H6)
    out['H6_v'] = arb(int(sp.Rational(H6_v).p)) / int(sp.Rational(H6_v).q)
    t0 = time.time()
    out['SN8'] = [sphere_sup_deriv(NH8, x, mm) for mm in range(3)]       # sup_S ||D^m N[H8n]||_F
    out['SH6'] = [sphere_sup_deriv(H6, x, mm) for mm in range(5)]        # sup_S ||D^m H6||_F
    out['sphere_seconds'] = time.time() - t0
    return out

# ---- rigorous sup over the unit sphere of ||D^m P||_F for a homogeneous polynomial P (box cover of the cube surface)
def _ipow(lo, hi, a):
    """interval power of exact floats (outward rounded)."""
    if a == 0:
        return np.ones_like(lo), np.ones_like(hi)
    pl, ph = lo ** a, hi ** a
    if a % 2 == 1:
        L, H = pl, ph
    else:
        L = np.where(lo >= 0, pl, np.where(hi <= 0, ph, 0.0))
        H = np.maximum(pl, ph)
    return np.nextafter(L, -np.inf), np.nextafter(H, np.inf)

def _imul(al, ah, bl, bh):
    c = np.stack([al * bl, al * bh, ah * bl, ah * bh])
    return np.nextafter(c.min(0), -np.inf), np.nextafter(c.max(0), np.inf)

def sphere_sup_deriv(P, x, m, g=32):
    """certified upper bound of sup_{|w| = 1} ||D^m P(w)||_F (Frobenius with multiplicities m!/alpha!, >= the
    operator norm) for homogeneous P of degree n > m with rational coefficients.  The unit sphere is the radial
    image of the cube surface |y|_inf = 1 (8 faces, each cut into g^3 boxes); each entry d^alpha P (homogeneous of
    degree n - m) is enclosed on each box by interval arithmetic and divided by (min |y| on the box)^(n - m)."""
    import sympy as sp
    from fractions import Fraction
    n = sp.Poly(P, *x).total_degree()
    entries = []
    for alpha in [al for al in np.ndindex(*(m + 1,) * 4) if sum(al) == m]:
        d = P
        for i, k in enumerate(alpha):
            if k: d = sp.diff(d, x[i], k)
        d = sp.expand(d)
        if d == 0:
            continue
        w = math.factorial(m)
        for k in alpha: w //= math.factorial(k)
        terms = []
        for mon, c in sp.Poly(d, *x).terms():
            fr = Fraction(int(sp.Rational(c).p), int(sp.Rational(c).q))
            cl = float(fr); ch = cl
            if Fraction(cl) != fr:
                cl = float(np.nextafter(cl, -np.inf)); ch = float(np.nextafter(ch, np.inf))
            terms.append((mon, cl, ch))
        entries.append((w, terms))
    deg = n - m
    edges = np.linspace(-1.0, 1.0, g + 1)                               # exact dyadics for g a power of 2
    J = np.array(list(np.ndindex(g, g, g)))
    best = 0.0
    for face in range(4):
        for sgn in (-1.0, 1.0):
            others = [i for i in range(4) if i != face]
            lo = np.zeros((len(J), 4)); hi = np.zeros((len(J), 4))
            lo[:, face] = hi[:, face] = sgn
            for c, i in enumerate(others):
                lo[:, i] = edges[J[:, c]]; hi[:, i] = edges[J[:, c] + 1]
            pw = {(i, a): _ipow(lo[:, i], hi[:, i], a) for i in range(4) for a in range(n + 1)}
            tot = np.zeros(len(J))
            for (w, terms) in entries:
                sl = np.zeros(len(J)); sh = np.zeros(len(J))
                for mon, cl, ch in terms:
                    ml, mh = pw[(0, mon[0])]
                    for i in (1, 2, 3):
                        ml, mh = _imul(ml, mh, *pw[(i, mon[i])])
                    tl, th = _imul(ml, mh, np.full(len(J), cl), np.full(len(J), ch))
                    sl = np.nextafter(sl + tl, -np.inf); sh = np.nextafter(sh + th, np.inf)
                mag = np.maximum(np.abs(sl), np.abs(sh))
                tot = np.nextafter(tot + np.nextafter(w * mag * mag, np.inf), np.inf)
            ymin2 = np.ones(len(J))
            for i in others:
                ymin2 = ymin2 + np.where((lo[:, i] <= 0) & (hi[:, i] >= 0), 0.0, np.minimum(lo[:, i] ** 2, hi[:, i] ** 2))
            s = np.nextafter(np.sqrt(np.nextafter(ymin2, -np.inf)), -np.inf)       # IEEE sqrt, rounded down
            pden = np.ones(len(J))
            for _ in range(deg):
                pden = np.nextafter(pden * s, -np.inf)
            val = np.nextafter(np.sqrt(tot), np.inf) / pden
            best = max(best, float(np.nextafter(val.max(), np.inf)))
    return arb(best) * (1 + arb('1e-12'))

def ff(n, m):
    """falling factorial n!/(n-m)!"""
    v = 1
    for i in range(m):
        v *= (n - i)
    return v

# =====================================================================================================
# 3. lattice constants: Z_kappa (wp2b, certified), c8 = P_8[psi_1](e0), c6 (direct-sum control), c42~
# =====================================================================================================
def U8(t2):
    return 256 * t2 ** 4 - 448 * t2 ** 3 + 240 * t2 ** 2 - 40 * t2 + 1

def U6(t2):
    return 64 * t2 ** 3 - 80 * t2 ** 2 + 24 * t2 - 1

def lattice_constants(R=10.0, log=print):
    """c8 := P_8[psi_1](e0) = -(1/2) sum_{l != 0} U_8(l_0/|l|) |l|^-10 (doc A2 (a), L4 (c)) and the control
    c6 = (8/3) sum U_6(l_0/|l|) |l|^-8 (doc L3 (e)); exact multiplicities of (|m|^2, m_0^2), l = m/4, |l| <= R, plus
    Lemma A4 tails with |U_8| <= 9, |U_6| <= 7 on [-1, 1]."""
    t0 = time.time()
    Z = WB.Z_enclosures(1600)
    m = np.rint(WB.lattice_points(R) * 4).astype(np.int64)
    N = (m * m).sum(1)
    m = m[N > 0]; N = N[N > 0]
    t = m[:, 0] ** 2
    key = N * (1 << 20) + t
    uq, cnt = np.unique(key, return_counts=True)
    S8 = arb(0); S6 = arb(0)
    for k, c in zip(uq.tolist(), cnt.tolist()):
        n = k >> 20; tt = k & ((1 << 20) - 1)
        t2 = arb(tt) / n
        S8 += c * U8(t2) * (arb(16) / n) ** 5
        S6 += c * U6(t2) * (arb(16) / n) ** 4
    Nmax = int(round(16 * R * R))
    T = (arb(Nmax + 1) / 16).sqrt()
    tail8 = 9 * WB.A4_tail(T, 10); tail6 = 7 * WB.A4_tail(T, 8)
    c8 = -(S8 + hull(-tail8, tail8)) / 2
    c6d = arb(8) / 3 * (S6 + hull(-tail6, tail6))
    # P_6[psi_2](y) = (1/24) sum_{l != 0} |y|^6 C_6^(3)(y.lhat/|y|) |l|^-12 (doc L3 (a), A2 (b)); |C_6^(3)| <= 462
    tail12 = 462 * WB.A4_tail(T, 12) / 24
    def p6psi2(num):                                     # y = num-weighted direction: t^2 = (w.m)^2/(|w|^2 N)
        w = np.array(num, dtype=np.int64); ww = int((w * w).sum())
        tt = (m @ w) ** 2
        kk = N * (1 << 24) + tt
        u2, c2 = np.unique(kk, return_counts=True)
        s = arb(0)
        for k, c in zip(u2.tolist(), c2.tolist()):
            n = k >> 24; t_ = k & ((1 << 24) - 1)
            s += c * C63(arb(t_) / (ww * n)) * (arb(16) / n) ** 6
        return s / 24 + hull(-tail12, tail12)
    v1 = p6psi2((1, 0, 0, 0)); v2 = p6psi2((1, 1, 0, 0)); v3 = p6psi2((2, 1, 0, 0))
    alpha2 = (v1 + v2) / 2; beta2 = (v1 - v2) * 4 / 3
    log('  lattice sums |l| <= %g (%d points, %d classes): c8 = %s, c6(direct) = %s, P6[psi_2] = alpha2 p2^3 + beta2 H6 with '
        'alpha2 = %s, beta2 = %s  [%.1fs]' % (R, len(m), len(uq), c8.str(12), c6d.str(12), alpha2.str(10), beta2.str(10),
                                              time.time() - t0))
    return dict(Z=Z, c8=c8, c6_direct=c6d, npts=int(len(m)), alpha2=alpha2, beta2=beta2, v3=v3)

def C63(t2):
    """Gegenbauer C_6^(3)(t) = 1792 t^6 - 1680 t^4 + 360 t^2 - 10 as a polynomial in t^2 (coefficients from sympy)."""
    return _C63[0] + _C63[1] * t2 + _C63[2] * t2 ** 2 + _C63[3] * t2 ** 3

def _c63_coeffs():
    import sympy as sp
    t = sp.symbols('t')
    P = sp.Poly(sp.expand(sp.gegenbauer(6, 3, t)), t)
    cs = [P.coeff_monomial(t ** (2 * k)) for k in range(4)]
    assert all(P.coeff_monomial(t ** (2 * k + 1)) == 0 for k in range(3))
    assert sp.expand(sp.gegenbauer(6, 3, 1)) == 462
    return [arb(int(c)) for c in cs]
_C63 = _c63_coeffs()

# =====================================================================================================
# 4. homogeneous Lie bounds and their derivative norms (Stage A remainders)
# =====================================================================================================
DMAX = 40

def Bm(d, m):
    """multi-direction Cauchy for a homogeneous P of degree d with |P(y)| <= C Lie(y)^d:
    ||D^m P(x)|| <= C Bm(d, m) |x|^(d-m), Bm = (d-m)^m (d/(d-m))^d (m < d), d^d (m = d), 0 (m > d)."""
    if m > d:
        return arb(0)
    if m == d:
        return arb(d) ** d
    return arb(d - m) ** m * (arb(d) / (d - m)) ** d

class Lie:
    def __init__(self, Z):
        self.Zb_ = {2.0: arb(Z[2.0].upper()), 2.5: arb(Z[2.5].upper()), 3.0: arb(Z[3.0].upper())}
    def Zb(self, kap):
        """Z_kappa <= 16^(kappa - b) Z_b for every base b <= kappa (|l| >= 1/2)."""
        best = None
        for b, zb in self.Zb_.items():
            if b <= kap + 1e-12:
                v = arb(16) ** (arb(kap) - arb(b)) * zb if kap != b else zb
                if best is None or v.upper() < best.upper():
                    best = v
        if best is None:
            raise Refuse('Zb: kappa below the certified bases')
        return best
    def C1(self, d):
        """|P_d[psi_1](y)| <= C1 Lie(y)^d (doc L3 (i))."""
        return arb(d + 1) / 2 * self.Zb((d + 2) / 4.0)
    def C2(self, d, eps, parts=False):
        """|P_d[psi~_2](y)| <= C2 Lie(y)^d (doc L3 (i), second line), k >= 3 summed geometrically."""
        def Tk(k):
            return (WB.binom_half(k) / (2 * k - 1) * eps ** (k - 2) * math.comb(d + 4 * k - 3, d)
                    * self.Zb((4 * k - 2 + d) / 4.0))
        q = 16 * eps * (arb(d + 10) / 10) ** 4
        if not (q.upper() < 1):
            raise Refuse('C2: k-series ratio >= 1')
        return (Tk(2), Tk(3) / (1 - q)) if parts else Tk(2) + Tk(3) / (1 - q)

def homog_sum(L, ai, kind, par, dmin, m, r, base=None):
    """sum over degrees d >= dmin, d = par mod 4, of coef(d) Bm(d, m) r^(d - m)   (base None), or
    r^(d - base) (scaled: the result times r^(base - m) bounds ||D^m||).  coef = a^(d-6) C1(d) (psi_1) or
    a^(d-2) C2(d) (psi~_2).  Terms d <= DMAX explicit; the tail is geometric (ratio bound in BUILD-WP5.md)."""
    eps = ai.p(4)
    r = arb(r)
    tot = arb(0)
    d = dmin
    ex = (lambda d: d - m) if base is None else (lambda d: d - base)
    while d <= DMAX:
        if d >= m:
            coef = ai.p(d - 6) * L.C1(d) if kind == 1 else ai.p(d - 2) * L.C2(d, eps)
            e = ex(d)
            if e < 0:
                raise Refuse('homog_sum: negative r exponent')
            tot += coef * Bm(d, m) * r ** e
        d += 4
    d0 = d
    r4 = r ** 4
    if kind == 1:
        coef0 = ai.p(d0 - 6) * arb(d0 + 1) / 2 * L.Zb((d0 + 2) / 4.0)
        q = ai.p(4) * r4 * 16 * arb(d0 + 5) / (d0 + 1) * 6
    else:
        coef0 = ai.p(d0 - 2) * arb(2) ** (d0 + 2) * L.Zb((6 + d0) / 4.0) / (1 - 256 * eps)
        q = ai.p(4) * r4 * 256 * 6
    if not (q.upper() < 1):
        raise Refuse('homog_sum: tail ratio >= 1')
    tail = coef0 * Bm(d0, m) * r ** ex(d0) / (1 - q)
    return tot + tail

def Aee(L, ai, m, r, base=None):
    """A_ee' = sum_{d = 0 mod 4} [a^(d-6) P_d[psi_1] (d >= 12) + a^(d-2) P_d[psi~_2] (d >= 8)]."""
    return homog_sum(L, ai, 1, 0, 12, m, r, base) + homog_sum(L, ai, 2, 0, 8, m, r, base)

def Aoo(C, ai, m, r):
    """A_oo' = sum_{d = 2 mod 4} [a^(d-6) P_d[psi_1] (d >= 14; P_10[psi_1] = 0) + a^(d-2) P_d[psi~_2] (d >= 6)].
    The d = 6 term a^4 P_6[psi~_2] = a^4 [alpha2 p2^3 + beta2 H6] + a^4 sum_{k >= 3} eps^(k-2) P_6[psi_k] is split:
    the psi_2 part exactly (radial p2^3, sphere-sup norms of H6), the k >= 3 part by its Lie bound."""
    L = C['L']
    r = arb(r); sig = r * r
    rad = WB.radial_norms([sig ** 3, 3 * sig ** 2, 6 * sig, arb(6), arb(0)], r)          # ||D^m |x|^6||, |x| <= r
    exact6 = C['alpha2abs'] * rad[m] + C['beta2abs'] * C['SH6'][m] * r ** (6 - m)
    _, T3 = L.C2(6, ai.p(4), parts=True)
    k3 = T3 * Bm(6, m) * r ** (6 - m)
    return (ai.p(4) * (exact6 + k3) + homog_sum(L, ai, 2, 2, 10, m, r) + homog_sum(L, ai, 1, 2, 14, m, r))

def c42_tilde(L, ai, c42):
    """|c42~(eps)| <= c42 + sum_{k >= 3} eps^(k-2) |b_k| C(4k+1, 4) Z_{k+1/2}/(2k-1)  (A2 (b), L3 (f))."""
    eps = ai.p(4)
    T3 = eps * WB.binom_half(3) * math.comb(13, 4) * L.Zb(3.5) / 5
    q = 16 * eps * arb('3.85')
    if not (q.upper() < 1):
        raise Refuse('c42~ tail')
    return U(c42 + T3 / (1 - q))

# =====================================================================================================
# 5. Stage A: the core |zeta| <= 1 (Z1: |zeta| <= 1/4 with homogeneous scaling; Z2: 1/4 <= |zeta| <= 1)
# =====================================================================================================
R1W, WW = 0.25, 0.75

def radB(sig):
    """radial factors (scaled): B1 = 1/(sigma^2 fh): |B1| = r^-4 b0, ||D B1|| <= r^-5 b1, ||D^2 B1|| <= r^-6 b2;
    B2 = sigma/fh: |B2| = r^2/fh, ||D B2|| = r g1_, ||D^2 B2|| <= g2_."""
    fh = (1 + sig * sig).sqrt()
    b0 = 1 / fh
    b1 = 2 * (2 + 3 * sig * sig) / fh ** 3
    b2 = 4 * (6 + 15 * sig ** 2 + 12 * sig ** 4) / fh ** 5 + 2 * (2 + 3 * sig ** 2) / fh ** 3
    g1_ = 2 / fh ** 3
    g2_ = 12 * sig ** 2 / fh ** 5 + 2 / fh ** 3
    return fh, b0, b1, b2, g1_, g2_

def sigma2_part(sig, r):
    """P = |zeta|^4 exactly: Q1 = (4 + 6 sigma^2)/fh, ||D Q1|| = 2 r sigma (8 + 6 sigma^2)/fh^3,
    |ddbar Q1|_{gh0} = max((16 + 16 s^2 + 6 s^4)/fh^4, s^2 (8 + 6 s^2)/fh^4) <= 16."""
    fh = (1 + sig * sig).sqrt()
    q = (4 + 6 * sig * sig) / fh
    dq = 2 * r * sig * (8 + 6 * sig * sig) / fh ** 3
    rad = (16 + 16 * sig ** 2 + 6 * sig ** 4) / fh ** 4
    tan = sig ** 2 * (8 + 6 * sig ** 2) / fh ** 4
    ric = hull(U(rad), U(tan))
    return q, dq, U(ric)

def zone_Z1(C, ai, npc):
    """|zeta| <= 1/4: c = 0, Ph = a^2 (c8 H8n + c42~ |zeta|^4) + A_ee'.  Returns per-piece arb upper bounds."""
    out = []
    a2, a8 = ai.p(2), ai.p(8)
    pm = [Aee(C['L'], ai, m, R1W, base=8) for m in range(5)]          # ||D^m A_ee'|| <= pm r^(8-m), r <= 1/4
    n = C['nN8']                                                       # ||D^m N[H8n]|| <= n_m r^(8-m)
    c8, c42t = C['c8abs'], C['c42t']
    rad4 = [arb(1), arb(4), arb(12), arb(24), arb(24)]
    for k in range(npc):
        r0, r1 = R1W * k / npc, R1W * (k + 1) / npc
        r = hull(arb(r0), arb(r1)); sig = hull(arb(r0) ** 2, arb(r1) ** 2)
        R1_4 = arb(r1) ** 4
        fh, b0, b1, b2, g1_, g2_ = radB(sig)
        p0, p1, p2_, p3, p4 = pm
        # generic degree-8-scaled piece (A_ee'): |Q1| <= r^4 [..], ||DQ1|| <= r^3 [..], ric = (fh/2)[..]
        qR = R1_4 * (p2_ * b0 / 2 + R1_4 * p2_ / fh)
        dqR = p3 / 2 * b0 + p2_ * b0 + p2_ / 2 * b1 + R1_4 * (p3 / fh + p2_ * g1_)          # times r^3
        ricR = fh / 2 * ((p4 + 4 * p3 + 2 * p2_) / 2 * b0 + 2 * (p3 / 2 + p2_) * b1 + p2_ / 2 * b2
                         + R1_4 * (p4 / fh + 2 * p3 * g1_ + p2_ * g2_))
        # H8 part through N[H8n] (Lap H8n = 0)
        qH = R1_4 * n[0] * b0
        dqH = n[1] * b0 + n[0] * b1                                                        # times r^3
        ricH = fh / 2 * (n[2] * b0 + 2 * n[1] * b1 + n[0] * b2)
        qS, dqS, ricS = sigma2_part(sig, r)
        dqS3 = 2 * (8 + 6 * sig * sig) / fh ** 3                                           # dqS / r^3
        # full potential derivative norms scaled by r^(4-m) (det and positivity)
        km = [a2 * (c42t * rad4[mm] + c8 * C['BH8'] * ff(8, mm) * R1_4) + pm[mm] * R1_4 for mm in range(5)]
        detQ = R1_4 * km[2] ** 2 / 4
        ddet3 = km[2] * km[3] / 2                                                           # ||D det|| / r^3
        ricdet = fh / 4 * (km[3] ** 2 + km[2] * km[4])
        e = a8 * km[2] * fh / 2
        Q = a2 * (c42t * qS + c8 * qH) + qR + a8 * detQ
        DQ3 = a2 * (c42t * dqS3 + c8 * dqH) + dqR + a8 * ddet3                             # ||DQ|| <= r^3 DQ3
        ricLin = a2 * (ricS * c42t + c8 * ricH) + ricR + a8 * ricdet
        out.append(dict(zone='Z1', r=(r0, r1), Q=U(Q), ricLin=U(ricLin), grad2=U(R1_4 * fh * DQ3 ** 2 / 2), e=U(e),
                        parts=dict(ric_c42=up(a2 * ricS * c42t), ric_H8=up(a2 * c8 * ricH), ric_rest=up(ricR),
                                   ric_det=up(a8 * ricdet), Q_L2=up(a2 * (c42t * qS + c8 * qH)), Q_rest=up(qR))))
    return out

def zone_Z2(C, ai, nT):
    """1/4 <= |zeta| <= 1 in zeta coordinates, nT pieces in T = (|zeta| - 1/4)/(3/4)."""
    out = []
    a2, a8 = ai.p(2), ai.p(8)
    c8, c42t, c6 = C['c8abs'], C['c42t'], C['c6']
    L = C['L']
    for k in range(nT):
        T0, T1 = k / nT, (k + 1) / nT
        r0, r1 = R1W + WW * T0, R1W + WW * T1
        sig = hull(arb(r0) ** 2, arb(r1) ** 2); R = arb(r1)
        fh = (1 + sig * sig).sqrt()
        th = arb(r0) ** 2 / (1 + arb(r0) ** 4).sqrt()                  # smallest gh0 eigenvalue sigma/fh (increasing)
        th = arb(th.lower())
        cd = WB.c_derivs(T0, T1)
        Dc = WB.radial_norms(cd, r1)                                    # ||D^j c(|x|^2)||, |x| <= r1
        Ao = [Aoo(C, ai, m, r1) for m in range(5)]
        Ae = [Aee(L, ai, m, r1) for m in range(5)]
        pm = [sum(math.comb(m, j) * Dc[j] * Ao[m - j] for j in range(m + 1)) + Ae[m] for m in range(5)]
        H6m = [C['nH6'][j] * R ** (6 - j) for j in range(5)]
        l0m = [c6 / 2 * sum(math.comb(m, j) * Dc[j] * H6m[m - j] for j in range(m + 1)) for m in range(5)]
        # unscaled radial factors
        B1 = 1 / (sig ** 2 * fh)
        DB1 = 2 * R * (2 + 3 * sig ** 2) / (sig ** 3 * fh ** 3)
        D2B1 = 4 * R ** 2 * (6 + 15 * sig ** 2 + 12 * sig ** 4) / (sig ** 4 * fh ** 5) + 2 * (2 + 3 * sig ** 2) / (sig ** 3 * fh ** 3)
        B2 = sig / fh; DB2 = 2 * R / fh ** 3; D2B2 = 12 * sig ** 2 / fh ** 5 + 2 / fh ** 3
        def Q1(p):
            Nn = R ** 2 * p[2] / 2; N1 = R ** 2 / 2 * p[3] + R * p[2]; N2 = (R ** 2 * p[4] + 4 * R * p[3] + 2 * p[2]) / 2
            q = Nn * B1 + p[2] * B2
            dq = N1 * B1 + Nn * DB1 + p[3] * B2 + p[2] * DB2
            d2q = N2 * B1 + 2 * N1 * DB1 + Nn * D2B1 + p[4] * B2 + 2 * p[3] * DB2 + p[2] * D2B2
            return q, dq, d2q
        qR, dqR, d2qR = Q1(pm)
        _, dF8, _ = Q1(l0m)                                             # ||D F8|| (F8 = Q1[L0])
        nH = [C['nN8'][mm] * R ** (8 - mm) for mm in range(3)]
        qH = nH[0] * B1; dqH = nH[1] * B1 + nH[0] * DB1; d2qH = nH[2] * B1 + 2 * nH[1] * DB1 + nH[0] * D2B1
        qS, dqS, ricS = sigma2_part(sig, R)
        rad4 = [R ** 4, 4 * R ** 3, 12 * R ** 2, 24 * R, arb(24)]
        km = [l0m[mm] + a2 * (c42t * rad4[mm] + c8 * C['BH8'] * ff(8, mm) * R ** (8 - mm)) + pm[mm] for mm in range(5)]
        E0, E1, E2 = km[2] / 2, km[3] / 2, km[4] / 2
        detQ = E0 ** 2; ddet = 2 * E0 * E1; d2det = 2 * E1 ** 2 + 2 * E0 * E2
        e = a8 * E0 / th
        Q = C['supF8'] + a2 * (c42t * qS + c8 * qH) + qR + a8 * detQ
        DQ = dF8 + a2 * (c42t * dqS + c8 * dqH) + dqR + a8 * ddet
        ricLin = C["R6star"] + a2 * (ricS * c42t + c8 * d2qH / (2 * th)) + d2qR / (2 * th) + a8 * d2det / (2 * th)
        out.append(dict(zone='Z2', r=(r0, r1), Q=U(Q), ricLin=U(ricLin), grad2=U(DQ ** 2 / (2 * th)), e=U(e),
                        parts=dict(ric_c42=up(a2 * ricS * c42t), ric_H8=up(a2 * c8 * d2qH / (2 * th)),
                                   ric_rest=up(d2qR / (2 * th)), ric_det=up(a8 * d2det / (2 * th)),
                                   Q_L2=up(a2 * (c42t * qS + c8 * qH)), Q_rest=up(qR), gradF8=up(dF8))))
    return out

def stageA(C, ai, npc, nT):
    """C_core_F >= sup |F|/a^8 and C_core_R >= sup |Ric|_g/a^6 over |zeta| <= 1, for every a in ai."""
    t0 = time.time()
    pcs = zone_Z1(C, ai, npc) + zone_Z2(C, ai, nT)
    a8 = ai.p(8)
    Qmax = max((p['Q'] for p in pcs), key=lambda v: float(v.upper()))
    den = 1 - a8 * Qmax
    if not (den.lower() > 0):
        raise Refuse('a^8 sup|Q| >= 1')
    CF = U(Qmax / den)
    worst = None; CR = arb(0)
    for p in pcs:
        if not (p['e'].upper() < 1):
            raise Refuse('positivity: e >= 1 on a core piece')
        R = (p['ricLin'] / den + a8 * p['grad2'] / den ** 2) / (1 - p['e'])
        p['ric'] = up(R)
        if R.upper() > CR.upper():
            CR = U(R); worst = p
    for p in pcs:
        for kk in ('Q', 'ricLin', 'grad2', 'e'):
            p[kk] = up(p[kk])
    rem_R = up(CR - C['R6star'])
    return dict(a=[ai.slo, ai.shi], C_core_F=up(CF), C_core_R=up(CR), Qmax=up(Qmax),
                remainder_R_over_R6star=float(rem_R / up(C['R6star'])),
                remainder_F_over_supF8=float((up(CF) - up(C['supF8'])) / up(C['supF8'])),
                worst_R_piece=dict(zone=worst['zone'], r=worst['r'], parts=worst['parts']),
                e_core=max(p['e'] for p in pcs), pieces=pcs, seconds=time.time() - t0)

def core_constants(inp, lat, polys, ai_for_c42=None):
    L = Lie(lat['Z'])
    c42 = arb(lat['Z'][2.5].upper()) / 4                              # c42 = Z_{5/2}/4 (doc L3 (f))
    def mn(a_, b_):
        return a_ if a_.upper() <= b_.upper() else b_
    # derivative-norm constants at |x| = 1: the smaller of the Bombieri bound (n)_m [P]_B and the sphere-sup bound
    nN8 = [U(mn(ff(8, mm) * polys['BN8'], polys['SN8'][mm])) for mm in range(3)]
    nH6 = [U(mn(ff(6, mm) * polys['BH6'], polys['SH6'][mm])) for mm in range(5)]
    return dict(L=L, c42=c42, c8abs=U(abs(lat['c8'])), c6=U(inp['wp2a']['c6']), BN8=U(polys['BN8']),
                BH8=U(polys['BH8']), BH6=U(polys['BH6']), R6star=inp['wp2a']['R6star'], supF8=inp['wp2a']['supF8'],
                nN8=nN8, nH6=nH6, SH6=nH6, alpha2abs=U(abs(lat['alpha2'])), beta2abs=U(abs(lat['beta2'])))

def stageA_interval(Cbase, ai, npc, nT):
    C = dict(Cbase)
    C['c42t'] = c42_tilde(C['L'], ai, C['c42'])
    return stageA(C, ai, npc, nT)

# =====================================================================================================
# 6. Lemma G constant: Gammahat(k) by the c0/Doeblin formula (ANALYSIS-CHAIN 3.2-3.4) in arb; Gamma_expl
# =====================================================================================================
_S_GRID = np.linspace(1e-4, 1 - 1e-4, 2001)
_GH = {}

def gh_cert(k, N=4000, Tcap=1.6):
    """certified upper bound U(k) >= Gammahat(k') for every 0 <= k' <= k (D = 1):
    U = sum_{i<j} h (1 - c0lo(t_i)) + T (1 - c0lo(T))/c0lo(T), t_i = i h, T = j h, with
    c0lo(t) = exp(-(sqrt P + sqrt Q)^2) at one admissible s (chosen in double precision), P = 1/(4t(1-s)) - 2 log s,
    Q = 2 sqrt2 k t (1-s).  U is nondecreasing in k (fixed s, t), hence the bound for all k' <= k."""
    kf = float(k)
    if kf < 0 or not math.isfinite(kf):
        raise Refuse('gh_cert: bad k')
    ck = (kf, N, Tcap)
    if ck in _GH:
        return _GH[ck]
    h = Tcap / N
    tt = np.arange(1, N + 1) * h
    s = _S_GRID[None, :]
    P = 1 / (4 * tt[:, None] * (1 - s)) - 2 * np.log(s)
    Q = 2 * math.sqrt(2) * kf * tt[:, None] * (1 - s)
    val = (np.sqrt(P) + np.sqrt(Q)) ** 2
    jj = np.argmin(val, axis=1)
    c0f = np.exp(-val[np.arange(N), jj])
    c0all = np.concatenate([[0.0], c0f])
    left = np.concatenate([[0.0], np.cumsum(h * (1 - c0all[:-1]))])
    tot = left[1:] + tt * (1 - c0f) / c0f
    j = int(np.argmin(tot)) + 1                                         # T = j h
    ha = arb(h)                                                          # h is an exact double
    ka = arb(kf)
    r2 = arb(2).sqrt()
    def c0lo(i):
        t = arb(i) * ha
        sv = arb(float(_S_GRID[jj[i - 1]]))
        Pv = 1 / (4 * t * (1 - sv)) - 2 * sv.log()
        Qv = 2 * r2 * ka * t * (1 - sv)
        v = (Pv.sqrt() + nonneg(Qv).sqrt()) ** 2
        return arb((-v).exp().lower())
    S = ha                                                               # i = 0: c0 = 0
    for i in range(1, j):
        S += ha * (1 - c0lo(i))
    cT = c0lo(j)
    if not (cT > 0):
        raise Refuse('gh_cert: c0(T) not positive')
    Tt = arb(j) * ha
    Ub = U(S + Tt * (1 - cT) / cT)
    res = dict(U=Ub, T=float(j * h), c0T=float(cT.mid()), k=kf, float_value=float(tot[j - 1]))
    _GH[ck] = res
    return res

def gamma_expl(K, D):
    """Gamma_expl(K, D) = (3 D^2/4) exp((sqrt(1 + 2 log(3/2)) + sqrt(K D^2/sqrt2))^2) (ANALYSIS-CHAIN 3.2)."""
    K = nonneg(K)
    return 3 * D * D / 4 * (((1 + 2 * (arb(3) / 2).log()).sqrt() + ssqrt(K * D * D / arb(2).sqrt())) ** 2).exp()

def make_gamma(mode, N=4000):
    def G(K, D):
        if mode == 'expl':
            return U(gamma_expl(K, D)), dict(mode='expl')
        k = up(nonneg(K) * D * D)
        g = gh_cert(k, N=N)
        v1 = U(D * D * g['U'])
        v2 = U(gamma_expl(K, D))
        best = v1 if v1.upper() <= v2.upper() else v2
        return best, dict(mode='cert', k_up=k, Gh_up=up(g['U']), D2Gh=up(v1), Gexpl=up(v2))
    return G

# =====================================================================================================
# 7. Stage B: Theorem 2 in scaled arb arithmetic (valid on an a-interval touching 0)
# =====================================================================================================
def D_L5(ai, e):
    """L5: diam(X, g_R) <= sqrt(2 (1+e)) [4a + (17/16)^(1/4) (1/sqrt2 + 12 (pi - 2) a)], 0 < a <= 0.03, e < 1."""
    e = arb(e) if not isinstance(e, arb) else e
    if not (e.upper() < 1):
        raise Refuse('L5 needs e < 1 (positivity of omega_a); e = %s' % e)
    if frac(ai.shi) > frac(A0_STR):
        raise Refuse('L5 holds only for a <= 0.03')
    a = ai.ball
    return (2 * (1 + e)).sqrt() * (4 * a + (arb(17) / 16) ** (arb(1) / 4) * (1 / arb(2).sqrt() + 12 * (arb.pi() - 2) * a))

def theorem2(ai, Cq, CR, c0s, D, gamma, Lam0=2, iters=6, xi_override=None):
    """Theorem 2 (ANALYSIS-CHAIN 9) with Prop C0 (4), Lemma G (3) and the continuity closure (7), every quantity
    scaled by its power of a (BUILD-WP5.md 3.2).  Inputs are upper bounds valid for every a in ai:
    sup|F| <= Cq a^8 (so eta = -2 log(1 - Cq a^8) >= osc F), |Ric| <= CR a^6 = kappa, a^2 C0 <= c0s, diam <= D."""
    a = ai.ball
    a2, a3, a5, a6, a7, a8 = ai.p(2), ai.p(3), ai.p(5), ai.p(6), ai.p(7), ai.p(8)
    Cq, CR, c0s = arb(Cq), arb(CR), arb(c0s)
    q = Cq * a8
    Eh = 2 * Cq * ell(q)                       # eta / a^8
    eta = a8 * Eh
    if not (eta.exp().upper() < 3):
        raise Refuse('Prop C0 needs e^eta < 3')
    den = 1 - eta.expm1() / 2
    if not (den.lower() > 0):
        raise Refuse('Prop C0 denominator')
    e1h = Eh * g1(eta / 2)                     # eta1/a^8, eta1 = 2(1 - e^(-eta/2))
    e2h = Eh * phi1(eta / 2)                   # eta2/a^8, eta2 = 2(e^(eta/2) - 1)
    if not ((4 * c0s - 2 * CR * a8).lower() > 0):
        raise Refuse('Theorem 2 needs 4 C0 > sigma = 2 kappa')
    dh = 4 * Eh * g1(eta) + 2 * CR / c0s * (-eta).exp()           # delta / a^8
    kappa = CR * a6
    def step(Lam, Kfac):
        Lam = arb(Lam)
        G1, i1 = gamma(kappa, D)                                   # Gamma(kappa, D) on g_a
        G2, i2 = gamma(Kfac * Lam * kappa, ssqrt(Lam) * D)         # Gamma(K_t, sqrt(Lam) D) on g_t
        oh = (2 * e1h * G1 + 2 * e2h * G2) / den                   # o / a^8
        if xi_override is not None:
            xi = float(xi_override)
        else:
            xi = float(ssqrt(dh / (4 * c0s * oh)).mid())         # explicit choice x = 1 + sqrt(delta/(4 C0 o))
        if not (xi > 0 and math.isfinite(xi)):
            raise Refuse('x choice')
        xa = arb(xi)
        x = 1 + a * xa
        rho = dh / (2 * (ssqrt(xa * xa + x * a6 * dh) + xa))     # (R - 1)/a^7, exact
        Lh = a2 * Eh + a * rho + c0s * (2 + a * xa) * oh           # log(1 + eps2)/a^6 (log R <= R - 1)
        epsh = Lh * phi1(a6 * Lh)                                  # eps2 / a^6
        m1h = a3 * epsh + ssqrt(2 * epsh + a6 * epsh ** 2 + a2 * Eh * g1(eta))   # (M2 - 1)/a^3
        d3 = a5 * Eh * phi1(eta) + eta.exp() * m1h                 # d2 / a^3, d2 = M2 e^eta - 1
        M = 1 + a3 * m1h
        Uv = M * eta.exp()                                         # max(upper, 1/lower) = M e^eta
        fin(Uv, 'theorem2'); fin(d3, 'theorem2')
        return dict(Lam=up(Lam), oh=up(oh), o=up(a8 * oh), xi=xi, eps2_over_a6=up(epsh), M=up(M), U=up(Uv),
                    d2_over_a3=up(d3), d2=up(a3 * d3), G1=up(G1), G2=up(G2), g1info=i1, g2info=i2)
    cont = step(Lam0, Lam0)                                        # K_t = kappa' = Lam0 kappa, diameter sqrt(Lam0) D
    closes = cont['U'] < Lam0
    res = dict(a=[ai.slo, ai.shi], inputs=dict(Cq=up(Cq), CR=up(CR), c0s=up(c0s), D=up(D)),
               eta=up(eta), eta_over_a8=up(Eh), delta_over_a8=up(dh), continuity=cont, closes=bool(closes))
    if not closes:
        return res
    Lam = arb(Lam0); seq = []
    fsteps = []
    for _ in range(iters):
        st = step(Lam, 0)                                          # t = 1: Ric(omega_CY) = 0, K_t = 0
        fsteps.append(st)
        Lam = arb(st['U']); seq.append(st['U'])
    res.update(final=fsteps[-1], Lambda_seq=seq, d2=fsteps[-1]['d2'], d2_over_a3=fsteps[-1]['d2_over_a3'],
               smallness=cont['U'])
    return res

def stageB_interval(inp, SA, ai, gamma, extra_e=None):
    """inputs for ai (echoed), then Theorem 2."""
    rb = containing(inp['wp2b']['rows'], ai)
    p3 = containing(inp['wp3']['pieces'], ai)
    Cq = max(SA['C_core_F'], p3['C_out'])
    CR = max(SA['C_core_R'], p3['R_out'])
    es = [rb['e'], p3['e_out'], inp['wp3']['e_out'], SA['e_core']] + ([extra_e] if extra_e is not None else [])
    e = max(es)
    D = D_L5(ai, arb(e))
    r = theorem2(ai, Cq, CR, rb['c0s'], D, gamma)
    r['input_sources'] = dict(C_core_F=SA['C_core_F'], C_out=p3['C_out'], C_core_R=SA['C_core_R'], R_out=p3['R_out'],
                              c0s=rb['c0s'], wp2b_row=list(rb['a']), wp3_piece=list(p3['a']),
                              e_wp2b=rb['e'], e_wp3=p3['e_out'], e_core_stageA=SA['e_core'], e_used=e,
                              D=up(D), Cq_used=Cq, CR_used=CR)
    return r

# =====================================================================================================
# 8. double-precision float path (X1) and reference sampling (X4)
# =====================================================================================================
TABLE_R2_L5_EXPL = {'1e-3': 2.931e-6, '3e-3': 8.303e-5, '1e-2': 3.594e-3}   # ANALYSIS-CHAIN sec. 0, row "R2, D_L5, expl"

def x1_float_path(log):
    """brief inputs: C_q = 2.7e4, kappa = 2.7e6 a^6, B+ = C0 = 2/a^2, D = sqrt(1 + 140 a^4)(1.0153 + 25.33 a),
    Gamma = Gamma_expl (wp4 explicit.py); run through the SAME theorem2 code at point values of a."""
    sys.path.insert(0, os.path.join(K3, 'wp4'))
    rows = {}
    try:
        import explicit as EX                                      # read-only double reference
    except Exception as ex_:
        EX = None
        log('  (wp4 explicit.py not importable: %s)' % ex_)
    G = make_gamma('expl')
    ok = True
    for s, ref in TABLE_R2_L5_EXPL.items():
        ai = AInt(s, s)
        a = arb(s)
        D = (1 + 140 * a ** 4).sqrt() * (arb('1.0153') + arb('25.33') * a)
        r = theorem2(ai, arb('2.7e4'), arb('2.7e6'), arb(2), D, G)
        d = r['d2']
        refw = EX.run(float(s), EX.r2, 'L5') if EX is not None else None
        rel = abs(d - ref) / ref
        relw = abs(d - refw) / refw if refw else None
        rows[s] = dict(d2_float_path=d, table=ref, rel_table=rel, wp4_explicit=refw, rel_wp4=relw)
        ok = ok and rel <= 1e-3
        log('  X1 a=%s: d2 = %.6e  table %.4g (rel %.2e)  wp4 explicit.py %s (rel %s)'
            % (s, d, ref, rel, '%.6e' % refw if refw else '-', '%.2e' % relw if relw else '-'))
    return dict(rows=rows, pass_=bool(ok))

def x4_reference(SAmap, ai_list, log, n_each=300, seed=5):
    """double-precision full F and full Ricci (probe1 k1c/k1d with the C-infinity window, as patched by wp2b)
    at core points |zeta| <= 1, against the Stage A bounds of the a-interval containing a."""
    k1c, k1d = WB.reference_setup()
    rng = np.random.default_rng(seed)
    out = {}
    ok = True
    for a in (0.01, 0.03):
        ai = None
        for x in ai_list:
            if frac(x.slo) < frac(repr(a)) <= frac(x.shi):
                ai = x
        SA = SAmap[ai.key()]
        d = rng.normal(size=(2 * n_each, 4)); d /= np.linalg.norm(d, axis=1)[:, None]
        rz = np.concatenate([0.001 + 0.249 * rng.random(n_each), 0.25 + 0.75 * rng.random(n_each)])
        z = a * rz[:, None] * d
        ax = np.array([[1.0, 0, 0, 0], [0, 1.0, 0, 0], [1, 1, 0, 0] / np.sqrt(2), [0.5, 0.5, 0.5, 0.5]])
        extra = np.concatenate([a * rr * ax for rr in (0.7432, 0.8830, 0.95, 0.1)])
        z = np.concatenate([z, extra])
        model = k1d.get_model(a)
        rest, _ = k1d.rest_jets(z, model, 0.5, a)
        o = k1c.evaluate_eh(z, a, rest, bis=False)
        Fv = np.abs(o['F']) / a ** 8; Rv = o['RIC'] / a ** 6
        if not (np.all(np.isfinite(Fv)) and np.all(np.isfinite(Rv))):
            raise FloatingPointError('X4: non-finite reference values')
        fmax, rmax = float(Fv.max()), float(Rv.max())
        good = fmax <= SA['C_core_F'] and rmax <= SA['C_core_R']
        ok = ok and good
        out[repr(a)] = dict(n_points=int(len(z)), interval=[ai.slo, ai.shi], sampled_max_F_over_a8=fmax,
                            C_core_F=SA['C_core_F'], sampled_max_Ric_over_a6=rmax, C_core_R=SA['C_core_R'],
                            argmax_F_abs_zeta=float(np.linalg.norm(z[np.argmax(Fv)]) / a),
                            argmax_R_abs_zeta=float(np.linalg.norm(z[np.argmax(Rv)]) / a),
                            ratio_F=fmax / SA['C_core_F'], ratio_R=rmax / SA['C_core_R'], dominates=bool(good))
        log('  X4 a=%g (%d pts, interval %s): max|F|/a^8 = %.2f <= %.2f ; max|Ric|/a^6 = %.6e <= %.6e : %s'
            % (a, len(z), ai.key(), fmax, SA['C_core_F'], rmax, SA['C_core_R'], good))
    return dict(rows=out, pass_=bool(ok))

# =====================================================================================================
# 9. drivers
# =====================================================================================================
def logger(path):
    fh = open(path, 'a', encoding='utf-8')
    def log(msg):
        line = time.strftime('%H:%M:%S ') + msg
        print(line, flush=True); fh.write(line + '\n'); fh.flush()
    return log

def a_intervals(width_str, n):
    from fractions import Fraction
    w = Fraction(width_str)
    out = []
    for i in range(n):
        lo, hi = w * i, w * (i + 1)
        out.append(AInt(_dec(lo), _dec(hi)))
    if frac(out[-1].shi) != frac(A0_STR):
        raise Refuse('a-intervals must end at 0.03')
    return out

def _dec(fr):
    from decimal import Decimal
    v = Decimal(fr.numerator) / Decimal(fr.denominator)
    s = format(v.normalize(), 'f')
    return '0' if s in ('0', '0.0') else s

def setup(log):
    t0 = time.time()
    inp = load_inputs()
    for k in ('wp2a', 'wp2b', 'wp3'):
        log('  input %s: %s (sha256 %s...), verdict seal %s claim %s'
            % (k, inp[k]['source'], inp[k]['sha256'][:12], inp[k]['seal']['sealed'], inp[k]['seal']['claim']))
    log('  WP2a: c6 in %s ; P* <= %.15g ; R*(c6=1) <= %.15g ; R6* <= %.6f ; sup|F8| <= %.6f'
        % (inp['wp2a']['c6'].str(20), float(inp['wp2a']['P_up'].upper()), float(inp['wp2a']['R1_up'].upper()),
           up(inp['wp2a']['R6star']), up(inp['wp2a']['supF8'])))
    for r in inp['wp2b']['rows']:
        log('  WP2b row a in [%s, %s]: e <= %.6e, a^2 C0 <= %.13f' % (r['a'][0], r['a'][1], r['e'], r['c0s']))
    for p in inp['wp3']['pieces']:
        log('  WP3 piece a in [%s, %s]: C_out <= %.4f, R_out <= %.4f, e_out <= %.6e' % (p['a'][0], p['a'][1], p['C_out'], p['R_out'], p['e_out']))
    polys = exact_polys()
    log('  exact polynomials: dim(W(F4)-inv harmonic octics) = %d, decics vanish = %s, [H8n]_B = %s, [N H8n]_B = %s, [H6]_B = %s'
        % (polys['dim8'], polys['p10_zero'], polys['BH8'].str(10), polys['BN8'].str(10), polys['BH6'].str(10)))
    lat = lattice_constants(10.0, log)
    ov = lat['c6_direct'].overlaps(inp['wp2a']['c6'])
    log('  control: c6 direct sum %s overlaps WP2a Ewald c6: %s' % (lat['c6_direct'].str(12), ov))
    if not ov:
        raise Refuse('c6 direct sum does not overlap the WP2a certified c6')
    pred = lat['alpha2'] + lat['beta2'] * polys['H6_v']
    ov3 = pred.overlaps(lat['v3'])
    log('  control: P6[psi_2] at (2,1,0,0)/sqrt5: direct %s vs alpha2 + beta2 H6 = %s, overlap %s'
        % (lat['v3'].str(10), pred.str(10), ov3))
    if not ov3:
        raise Refuse('P6[psi_2] is not alpha2 p2^3 + beta2 H6 on the test point')
    log('  sphere sups (box cover, g = 32): ||D^m N[H8n]||_F <= %s (Bombieri %s) ; ||D^m H6||_F <= %s (Bombieri %s) [%.1fs]'
        % ([round(up(v), 4) for v in polys['SN8']], [round(up(ff(8, k) * polys['BN8']), 1) for k in range(3)],
           [round(up(v), 4) for v in polys['SH6']], [round(up(ff(6, k) * polys['BH6']), 1) for k in range(5)],
           polys['sphere_seconds']))
    C = core_constants(inp, lat, polys)
    log('  c42 = Z_{5/2}/4 <= %.6f ; |c8| <= %.6f  [setup %.1fs]' % (up(C['c42']), up(C['c8abs']), time.time() - t0))
    return inp, lat, polys, C

def run_stageA(C, ais, npc, nT, log, cache=None):
    res = {}
    for ai in ais:
        if cache is not None and ai.key() in cache:
            res[ai.key()] = cache[ai.key()]; continue
        SA = stageA_interval(C, ai, npc, nT)
        res[ai.key()] = SA
        w = SA['worst_R_piece']
        log('  A a in %s: C_core_F <= %.4f (sup|F8| part %.4f, rel. rest %.2e) ; C_core_R <= %.2f (R6* + %.3e rel.; worst %s r=%s: c42 %.1f, H8 %.1f, rest %.1f, det %.1e) ; e_core <= %.2e [%.2fs]'
            % (ai.key(), SA['C_core_F'], up(C['supF8']), SA['remainder_F_over_supF8'], SA['C_core_R'],
               SA['remainder_R_over_R6star'], w['zone'], tuple(round(v, 4) for v in w['r']), w['parts']['ric_c42'],
               w['parts']['ric_H8'], w['parts']['ric_rest'], w['parts']['ric_det'], SA['e_core'], SA['seconds']))
        if cache is not None:
            cache[ai.key()] = SA
    return res

def run_stageB(inp, SAmap, ais, gamma, log, cache=None):
    res = {}
    for ai in ais:
        if cache is not None and ai.key() in cache:
            res[ai.key()] = cache[ai.key()]; continue
        r = stageB_interval(inp, SAmap[ai.key()], ai, gamma)
        res[ai.key()] = r
        if r['closes']:
            log('  B a in %s: Cq %.2f CR %.4e c0s %.10f D %.5f e %.2e | smallness M2(2,2k)e^eta <= %.6f < 2 | Lambda* %.8f | d2 <= %.6e, d2/a^3 <= %.3f'
                % (ai.key(), r['input_sources']['Cq_used'], r['input_sources']['CR_used'], r['input_sources']['c0s'],
                   r['input_sources']['D'], r['input_sources']['e_used'], r['smallness'], r['Lambda_seq'][-2],
                   r['d2'], r['d2_over_a3']))
        else:
            log('  B a in %s: continuity step does NOT close (bound %.6f)' % (ai.key(), r['continuity']['U']))
        if cache is not None:
            cache[ai.key()] = r
    return res

def leading(inp, SAmap, ais, gamma, log):
    """lim_{a->0} d2(a)/a^3: (i) the scaled bound evaluated at a = 0 exactly (the limit; formulas continuous at 0);
    (ii) the certified sup of d2/a^3 over (0, 1e-4] (no continuity argument)."""
    first = ais[0]
    SA = SAmap[first.key()]
    out = {}
    for tag, ai in (('a=0', AInt('0', '0')), ('(0,1e-4]', AInt('0', '0.0001'))):
        rb = containing(inp['wp2b']['rows'], first); p3 = containing(inp['wp3']['pieces'], first)
        Cq = max(SA['C_core_F'], p3['C_out']); CR = max(SA['C_core_R'], p3['R_out'])
        e = max(rb['e'], p3['e_out'], inp['wp3']['e_out'], SA['e_core'])
        D = D_L5(ai, arb(e))
        r = theorem2(ai, Cq, CR, rb['c0s'], D, gamma)
        out[tag] = dict(d2_over_a3=r['d2_over_a3'], Lambda_seq=r['Lambda_seq'], Cq=Cq, D=up(D), c0s=rb['c0s'],
                        Gamma1=r['final']['G1'], Gamma2=r['final']['G2'])
        log('  leading coefficient %s: d2/a^3 <= %.4f  (Cq %.2f, D %.6f, a^2 C0 %.13f, Gamma1 %.5f, Gamma2 %.5f)'
            % (tag, r['d2_over_a3'], Cq, up(D), rb['c0s'], r['final']['G1'], r['final']['G2']))
    # closed-form check 8 D0 sqrt(c0s/2 ... ): 2 sqrt(c0s * 8 Cq D0^2 gamma)
    return out

def a0_from(B, ais):
    a0 = '0'
    for ai in ais:
        r = B[ai.key()]
        if r['closes'] and r['smallness'] < 2:
            a0 = ai.shi
        else:
            break
    return a0

def jsonable(o):
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items() if not str(k).startswith('_')}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, arb):
        return [dn(o), up(o)]
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.bool_):
        return bool(o)
    return o

def check_all_finite(o, path='out'):
    if isinstance(o, dict):
        for k, v in o.items(): check_all_finite(v, path + '.' + str(k))
    elif isinstance(o, (list, tuple)):
        for i, v in enumerate(o): check_all_finite(v, '%s[%d]' % (path, i))
    elif isinstance(o, float):
        if not math.isfinite(o):
            raise FloatingPointError('non-finite number at %s' % path)

CTRL = dict(width='0.005', n=6, npc=8, nT=64, ghN=2000)
MAIN = dict(width='0.001', n=30, npc=16, nT=256, ghN=4000)

def controls(out_path):
    T0 = time.time()
    log = logger(os.path.join(HERE, 'controls_wp5.log'))
    log('=== wp5 controls === sha256(wp5.py) = %s' % sha_file(os.path.abspath(__file__)))
    res = dict(sha256_wp5=sha_file(os.path.abspath(__file__)), config=CTRL)
    inp, lat, polys, C = setup(log)
    res['inputs'] = {k: dict(source=inp[k]['source'], seal=inp[k]['seal'], sha256=inp[k]['sha256']) for k in inp}
    res['inputs']['wp2a']['echo'] = inp['wp2a']['echo']
    res['inputs']['wp2b']['rows'] = inp['wp2b']['rows']; res['inputs']['wp3']['pieces'] = inp['wp3']['pieces']
    res['inputs']['wp3']['e_out'] = inp['wp3']['e_out']
    res['constants'] = dict(c8=jsonable(lat['c8']), c6_direct=jsonable(lat['c6_direct']), c42=jsonable(C['c42']),
                            BH8=jsonable(polys['BH8']), BN8=jsonable(polys['BN8']), BH6=jsonable(polys['BH6']),
                            H8n=polys['H8n'], dim8=polys['dim8'], p10_zero=polys['p10_zero'])
    ais = a_intervals(CTRL['width'], CTRL['n'])
    log('-- Stage A (control resolution: Z1 %d pieces, Z2 nT %d)' % (CTRL['npc'], CTRL['nT']))
    t = time.time()
    SAmap = run_stageA(C, ais, CTRL['npc'], CTRL['nT'], log)
    res['stageA'] = SAmap; res['stageA_seconds'] = time.time() - t
    # ---- X2: Gammahat(0)
    t = time.time()
    g0 = gh_cert(0.0, N=CTRL['ghN']); g0m = gh_cert(0.0, N=MAIN['ghN'])
    x2 = dict(Gh0_upper_ctrl=up(g0['U']), Gh0_upper_main=up(g0m['U']), T=g0m['T'], c0T=g0m['c0T'],
              dist_from_4p5308=up(g0m['U']) - 4.5308, Gexpl_0_1=up(gamma_expl(arb(0), arb(1))),
              lower_ref=0.0702, upper_ref=4.5871)
    x2['pass_'] = bool(0.0702 <= x2['Gh0_upper_main'] <= 4.5871 and x2['Gexpl_0_1'] <= 4.58711)
    log('  X2: Gammahat(0) <= %.6f (N=%d), %.6f (N=%d) at T = %.4f, c0(T) = %.5f ; distance from 4.5308: %+.2e ; Gamma_expl(0,1) = %.6f : %s [%.1fs]'
        % (x2['Gh0_upper_ctrl'], CTRL['ghN'], x2['Gh0_upper_main'], MAIN['ghN'], x2['T'], x2['c0T'],
           x2['dist_from_4p5308'], x2['Gexpl_0_1'], x2['pass_'], time.time() - t))
    res['X2'] = x2
    # ---- X1
    log('-- X1 float path (brief inputs, same theorem2 code)')
    res['X1'] = x1_float_path(log)
    # ---- Stage B (control)
    log('-- Stage B (control resolution)')
    G = make_gamma('cert', N=CTRL['ghN'])
    t = time.time()
    B = run_stageB(inp, SAmap, ais, G, log)
    res['stageB'] = B; res['stageB_seconds'] = time.time() - t
    res['leading'] = leading(inp, SAmap, ais, G, log)
    res['a0'] = a0_from(B, ais)
    log('  a0 (largest a <= 0.03 with smallness certified on every interval below it) = %s' % res['a0'])
    # ---- X3 mutations
    log('-- X3 mutations')
    first = ais[0]; SA = SAmap[first.key()]
    rb = containing(inp['wp2b']['rows'], first); p3 = containing(inp['wp3']['pieces'], first)
    Cq = max(SA['C_core_F'], p3['C_out']); CR = max(SA['C_core_R'], p3['R_out'])
    e = max(rb['e'], p3['e_out'], inp['wp3']['e_out'], SA['e_core'])
    a0i = AInt('0', '0')
    D = D_L5(a0i, arb(e))
    base = theorem2(a0i, Cq, CR, rb['c0s'], D, G)['d2_over_a3']
    dbl = theorem2(a0i, 2 * Cq, CR, rb['c0s'], D, G)['d2_over_a3']
    ratio = dbl / base
    refused = False; msg = ''
    try:
        stageB_interval(inp, SA, first, G, extra_e=1.0)
    except Refuse as ex_:
        refused = True; msg = str(ex_)
    x3 = dict(leading=base, leading_doubled_Cq=dbl, ratio=ratio, sqrt2=math.sqrt(2), e_ge_1_refused=refused, refusal=msg)
    x3['pass_'] = bool(abs(ratio - math.sqrt(2)) < 0.01 and refused)
    log('  X3: doubling C_q: leading %.4f -> %.4f, ratio %.6f (sqrt2 = %.6f); forcing e = 1: refused = %s (%s) : %s'
        % (base, dbl, ratio, math.sqrt(2), refused, msg, x3['pass_']))
    res['X3'] = x3
    # ---- X4 reference sampling
    log('-- X4 double-precision samples of full F and full Ricci in the core')
    t = time.time()
    res['X4'] = x4_reference(SAmap, ais, log)
    res['X4']['seconds'] = time.time() - t
    # ---- NaN/inf guards
    log('-- NaN/inf guards')
    guards = {}
    def expect(name, fn, exc=(FloatingPointError, Refuse, ValueError, AssertionError)):
        try:
            fn(); guards[name] = False
        except exc:
            guards[name] = True
    expect('nan_arb', lambda: fin(arb('nan')))
    expect('inf_float', lambda: fin(float('inf')))
    expect('Cq_a8_ge_1', lambda: theorem2(AInt('0.02', '0.03'), arb('1e20'), arb(1e6), arb(2), arb(1.8), G))
    expect('D_beyond_0.03', lambda: D_L5(AInt('0.03', '0.031'), arb(0)))
    expect('e_eq_1', lambda: D_L5(AInt('0.01', '0.02'), arb(1)))
    expect('negative_k', lambda: gh_cert(-1.0))
    expect('nan_input_theorem2', lambda: theorem2(AInt('0.01', '0.02'), arb(float('nan')), arb(1e6), arb(2), arb(1.8), G))
    expect('non_finite_output_check', lambda: check_all_finite({'x': float('nan')}))
    res['guards'] = guards
    log('  guards: %s' % guards)
    check_all_finite(jsonable(res))
    res['ALL_CONTROLS_PASS'] = bool(res['X1']['pass_'] and x2['pass_'] and x3['pass_'] and res['X4']['pass_']
                                    and all(guards.values()))
    res['seconds'] = time.time() - T0
    log('ALL_CONTROLS_PASS = %s  (%.1f s)' % (res['ALL_CONTROLS_PASS'], res['seconds']))
    json.dump(jsonable(res), open(out_path, 'w'), indent=1)
    return res

def main(out_path, budget_s=520):
    if not os.path.exists(GO_FILE):
        print('refusing: %s does not exist (only the lead creates it)' % GO_FILE)
        sys.exit(2)
    main_body(out_path, budget_s, MAIN, STATE_DIR, os.path.join(HERE, 'run_wp5.log'))

def main_body(out_path, budget_s, cfg, state_dir, log_path):
    """the production run (called by main after the GO check); cfg = MAIN for the sealed run."""
    T0 = time.time()
    os.makedirs(state_dir, exist_ok=True)
    log = logger(log_path)
    st_path = os.path.join(state_dir, 'main_state.json')
    state = json.load(open(st_path)) if os.path.exists(st_path) else dict(A={}, B={}, cfg=cfg)
    if state.get('cfg') != cfg:
        raise Refuse('state in %s was written with a different configuration' % state_dir)
    log('=== wp5 main (resume: %d A, %d B intervals done) sha256 %s' % (len(state['A']), len(state['B']),
                                                                       sha_file(os.path.abspath(__file__))))
    inp, lat, polys, C = setup(log)
    ais = a_intervals(cfg['width'], cfg['n'])
    for ai in ais:
        if ai.key() in state['A']:
            continue
        SA = stageA_interval(C, ai, cfg['npc'], cfg['nT'])
        state['A'][ai.key()] = jsonable(SA)
        log('  A a in %s: C_core_F <= %.4f ; C_core_R <= %.2f ; e_core <= %.2e' % (ai.key(), SA['C_core_F'], SA['C_core_R'], SA['e_core']))
        json.dump(state, open(st_path, 'w'))
        if time.time() - T0 > budget_s:
            log('budget reached; rerun the same command (exit 3)'); sys.exit(3)
    G = make_gamma('cert', N=cfg['ghN'])
    for ai in ais:
        if ai.key() in state['B']:
            continue
        r = stageB_interval(inp, state['A'][ai.key()], ai, G)
        state['B'][ai.key()] = jsonable(r)
        log('  B a in %s: closes %s smallness %.6f d2 <= %.6e d2/a^3 <= %.3f' % (ai.key(), r['closes'], r.get('smallness', float('nan')),
                                                                              r.get('d2', float('nan')), r.get('d2_over_a3', float('nan'))))
        json.dump(state, open(st_path, 'w'))
        if time.time() - T0 > budget_s:
            log('budget reached; rerun the same command (exit 3)'); sys.exit(3)
    lead = leading(inp, state['A'], ais, G, log)
    B = state['B']
    a0 = a0_from(B, ais)
    g0 = gh_cert(0.0, N=cfg['ghN'])
    out = dict(status='COMPLETE', sha256_wp5=sha_file(os.path.abspath(__file__)), config=cfg,
               inputs={k: dict(source=inp[k]['source'], seal=inp[k]['seal'], sha256=inp[k]['sha256']) for k in inp},
               constants=dict(c8=jsonable(lat['c8']), c6_direct=jsonable(lat['c6_direct']), c42=jsonable(C['c42']),
                              BH8=jsonable(polys['BH8']), BN8=jsonable(polys['BN8']), Gh0_upper=up(g0['U'])),
               stageA={k: {kk: v for kk, v in s.items()} for k, s in state['A'].items()}, stageB=B,
               leading=jsonable(lead), a0=a0,
               summary=dict(C_core_F_max=max(s['C_core_F'] for s in state['A'].values()),
                            C_core_R_max=max(s['C_core_R'] for s in state['A'].values()),
                            d2_at_0p03=B[ais[-1].key()].get('d2'), smallness_max=max(b.get('smallness', 9) for b in B.values()),
                            leading_coefficient_upper=lead['a=0']['d2_over_a3'],
                            leading_coefficient_sup_0_to_1em4=lead['(0,1e-4]']['d2_over_a3'], a0=a0))
    check_all_finite(jsonable(out))
    json.dump(jsonable(out), open(out_path, 'w'), indent=1)
    log('main COMPLETE: a0 = %s, leading <= %.4f, d2(0.03) <= %.6e  (%.1f s this invocation)'
        % (a0, lead['a=0']['d2_over_a3'], out['summary']['d2_at_0p03'], time.time() - T0))

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--controls'); ap.add_argument('--main'); ap.add_argument('--budget', type=float, default=520)
    args = ap.parse_args()
    if args.controls:
        controls(args.controls)
    elif args.main:
        main(args.main, args.budget)
    else:
        ap.print_help()
