#!/usr/bin/env python
"""wp2a.py -- WP2a certificate build for the K3 Kummer core (DEFINITIONS-AND-LEMMAS-v2.md).

Rigorous enclosures (ball arithmetic: python-flint arb at 128 bits; outward-rounded IEEE interval
arithmetic in numpy for the angular cover) of
  P*  = sup |p(sigma)| over the cut-off window sigma in [r1^2, r2^2]            (L4 (b))
  R*  = sup R6(zeta) over the core shell r1 <= |zeta| <= r2                      (L4 (f), Fix 1 units |Ric|_g/a^6)
  c6  = (8/3) sum_{l != 0} U_6(l_0/|l|)/|l|^8                                   (L3 (e))
for two windows:
  'cinf'    : C-infinity step st = h(t)/(h(t)+h(1-t)), h = exp(-1/t), window |zeta| in [0.25, 1]  (lead decision)
  'quintic' : quintic step st5, window |zeta| in [0.5, 1]  (v1/K1d setting, known-answer control)

Modes
  --controls OUT.json   double-precision reproduction, mutation tests, NaN/inf checks, timing of a reduced
                        cover on the CONTROL window only
  --main OUT.json       the certified run; refuses unless <lab>/wp2/GO-WP2a exists; chunked and
                        resumable (state in <lab>/wp2/state/); rerun until status == COMPLETE
Run:  python <lab-tools>/memcap.py 8 python <lab>/wp2/wp2a.py --controls <lab>/wp2/controls.json

Rigor notes (details in BUILD-WP2a.md)
 * Every arb ball is converted to floats with an outward nextafter step (float() of an exact arb endpoint
   rounds to nearest, so one ulp outward is a valid bound).
 * numpy interval ops use only +, -, *, /, sqrt (IEEE-754 correctly rounded, round-to-nearest) followed by one
   ulp outward (np.nextafter).  No transcendental function is evaluated in floats.
 * exp() only enters through the step and is evaluated in arb (series arithmetic); the end pieces of the
   C-infinity step (T within END_B of 0 or 1) use Lemma C (C1), (C5), (C6).
 * R6 is reduced analytically to |H| * Bh(sigma) + sqrt(H^2 Gam + U A11 + k2 A00 + w A10) with four angular
   invariants (H = H6, U = |dH|^2 - (Im z.dH)^2, k2 = |ddbar H|_F^2, w = Re z^T (ddbar H) dH) on the unit
   sphere and five radial functions; directions are covered modulo G_J (order 192) by a Dirichlet cone.
"""
import os, sys, json, time, math, argparse, itertools
import numpy as np
from flint import arb, arb_series, ctx, fmpq

ctx.prec = 128
HERE = os.path.dirname(os.path.abspath(__file__))
GO_FILE = os.path.join(HERE, 'GO-WP2a')
STATE_DIR = os.path.join(HERE, 'state')
WINDOWS = {'cinf': (0.25, 1.0), 'quintic': (0.5, 1.0)}      # step -> (r1, r2) in |zeta|
END_B = 1.0 / 64            # initial T pieces; pieces touching T = 0 or 1 use Lemma C (needs width <= 1/8)
MUT = {'st': 1.0, 'lemmaC': 1.0}   # mutation factors (controls only): scale derivative enclosures of st

# =====================================================================================================
# 1. outward-rounded interval arithmetic on numpy arrays
# =====================================================================================================
_NINF, _PINF = -np.inf, np.inf
def _dn(a): return np.nextafter(a, _NINF)
def _up(a): return np.nextafter(a, _PINF)

class I:
    """interval [lo, hi] (numpy arrays); every op rounds outward by one ulp."""
    __slots__ = ('lo', 'hi')
    def __init__(self, lo, hi=None):
        self.lo = np.asarray(lo, dtype=float)
        self.hi = self.lo if hi is None else np.asarray(hi, dtype=float)
    def __add__(s, o):
        if isinstance(o, I): return I(_dn(s.lo + o.lo), _up(s.hi + o.hi))
        o = float(o); return I(_dn(s.lo + o), _up(s.hi + o))
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
        o = float(o)          # constants used in this file are dyadic, hence exact doubles
        if o >= 0: return I(_dn(s.lo * o), _up(s.hi * o))
        return I(_dn(s.hi * o), _up(s.lo * o))
    __rmul__ = __mul__
    def take(s, idx): return I(s.lo[idx], s.hi[idx])
    def mag(s): return np.maximum(np.abs(s.lo), np.abs(s.hi))

def isq(x):
    if not isinstance(x, I): return x * x
    a = x.lo * x.lo; b = x.hi * x.hi
    z = (x.lo <= 0) & (x.hi >= 0)
    return I(np.where(z, 0.0, _dn(np.minimum(a, b))), _up(np.maximum(a, b)))
def isqrt(x):
    """sqrt of the nonnegative part (callers guarantee the true value is >= 0)."""
    return I(np.maximum(_dn(np.sqrt(np.maximum(x.lo, 0.0))), 0.0), _up(np.sqrt(np.maximum(x.hi, 0.0))))
def irecip(x):
    """1/x for intervals not containing 0 (caller guarantees)."""
    return I(_dn(1.0 / x.hi), _up(1.0 / x.lo))
def iabs(x):
    z = (x.lo <= 0) & (x.hi >= 0)
    return I(np.where(z, 0.0, np.minimum(np.abs(x.lo), np.abs(x.hi))), np.maximum(np.abs(x.lo), np.abs(x.hi)))
def imax0(x): return I(np.maximum(x.lo, 0.0), np.maximum(x.hi, 0.0))
def isect(a, b): return I(np.maximum(a.lo, b.lo), np.minimum(a.hi, b.hi))
def iwhere(m, a, b): return I(np.where(m, a.lo, b.lo), np.where(m, a.hi, b.hi))

class Du:
    """first-order forward AD over intervals: value v (I) and gradient g (list of I)."""
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
def dsq(x):
    if isinstance(x, Du): return Du(isq(x.v), [(x.v * a) * 2.0 for a in x.g])
    return isq(x)

def check_finite(name, *ivs):
    """NaN / inf / inverted-interval guard; raises."""
    for k, x in enumerate(ivs):
        lo, hi = np.asarray(x.lo), np.asarray(x.hi)
        if not (np.all(np.isfinite(lo)) and np.all(np.isfinite(hi))):
            raise FloatingPointError('%s[%d]: non-finite interval endpoint' % (name, k))
        if np.any(lo > hi):
            raise FloatingPointError('%s[%d]: empty interval (lo > hi): enclosure inconsistency' % (name, k))

def a2f(x):
    """arb ball -> outward float pair; raises on non-finite."""
    if not x.is_finite():
        raise FloatingPointError('non-finite arb ball %s' % x)
    lo = float(x.lower()); hi = float(x.upper())
    return _dn(lo), _up(hi)

# =====================================================================================================
# 2. angular invariants of H6 on the unit sphere (generic: floats, arb, I, Du)
# =====================================================================================================
def angular_invariants(x, sq):
    """x = (x0, x1, x2, x3) a UNIT vector (z1 = x0 + i x1, z2 = x2 + i x3).  Returns H, U, k2, w with
       H  = H6(x) = (3/8) Re E,  E = zb1 A1 + zb2 A2,  A1 = z1^5 - 5 z1 z2^4,  A2 = z2^5 - 5 z1^4 z2   (L3 (e))
       D_i = d_{z_i} H,  K_ij = d_i dbar_j H,  rho = z.D
       U  = |D|^2 - (Im rho)^2,  k2 = sum |K_ij|^2,  w = Re sum_ij z_i K_ij D_j."""
    x0, x1, x2, x3 = x
    def cm(a, b): return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])
    def ca(a, b): return (a[0] + b[0], a[1] + b[1])
    def cs(a, b): return (a[0] - b[0], a[1] - b[1])
    def ck(a, k): return (a[0] * k, a[1] * k)
    def cc(a): return (a[0], -a[1])
    def csq(a): return (sq(a[0]) - sq(a[1]), (a[0] * a[1]) * 2.0)
    z1 = (x0, x1); z2 = (x2, x3)
    z1_2 = csq(z1); z2_2 = csq(z2); z1_4 = csq(z1_2); z2_4 = csq(z2_2)
    z1_3 = cm(z1_2, z1); z2_3 = cm(z2_2, z2); z1_5 = cm(z1_4, z1); z2_5 = cm(z2_4, z2)
    A1 = cs(z1_5, ck(cm(z1, z2_4), 5.0)); A2 = cs(z2_5, ck(cm(z1_4, z2), 5.0))
    d1A1 = ck(cs(z1_4, z2_4), 5.0)                    # d_1 A_1 ; d_2 A_2 = -d_1 A_1
    d2A1 = ck(cm(z1, z2_3), -20.0); d1A2 = ck(cm(z1_3, z2), -20.0)
    d2A2 = (-d1A1[0], -d1A1[1])
    zb1 = cc(z1); zb2 = cc(z2)
    E = ca(cm(zb1, A1), cm(zb2, A2))
    H = E[0] * 0.375
    D1 = ck(ca(ca(cm(zb1, d1A1), cm(zb2, d1A2)), cc(A1)), 0.1875)
    D2 = ck(ca(ca(cm(zb1, d2A1), cm(zb2, d2A2)), cc(A2)), 0.1875)
    K11 = d1A1[0] * 0.375                             # (3/16)(d1A1 + conj d1A1); K22 = -K11 (H harmonic)
    K12 = ck(ca(d1A2, cc(d2A1)), 0.1875)              # (3/16)(d_1 A_2 + conj(d_2 A_1)); K21 = conj K12
    K21 = cc(K12)
    rho = ca(cm(z1, D1), cm(z2, D2))
    U = sq(D1[0]) + sq(D1[1]) + sq(D2[0]) + sq(D2[1]) - sq(rho[1])
    k2 = sq(K11) * 2.0 + (sq(K12[0]) + sq(K12[1])) * 2.0
    KD1 = ca(ck(D1, K11), cm(K12, D2))
    KD2 = ca(cm(K21, D1), ck(D2, -K11))
    w = cm(z1, KD1)[0] + cm(z2, KD2)[0]
    return H, U, k2, w

def fsq(v): return v * v

# =====================================================================================================
# 3. radial functions (arb).  c(sigma) = st(T), T = (sqrt(sigma) - r1)/w, w = r2 - r1   (L4 (b), 1.5)
# =====================================================================================================
def tball(lo, hi):
    if lo == hi: return arb(lo)
    return arb(lo).union(arb(hi))

def st_of_series(Ts, step):
    """st(T) for an arb_series T (constant term strictly inside (0,1) for cinf)."""
    if step == 'quintic':
        return Ts ** 5 * (126 - 420 * Ts + 540 * Ts ** 2 - 315 * Ts ** 3 + 70 * Ts ** 4)
    g = 1 / Ts - 1 / (1 - Ts)                            # Lemma C (C2): st = 1/(1 + e^g)
    if Ts.coeffs()[0].mid() <= 0.5:                      # g >= 0 side: st = e^-g/(1 + e^-g), denominator in [1, 2]
        E = (-g).exp()
        return E / (1 + E)
    return 1 / (1 + g.exp())                             # g <= 0 side: denominator in [1, 2]

def _mutate(S, f):
    if f == 1.0: return S
    c = S.coeffs(); n = S.prec
    return arb_series([c[0]] + [f * v for v in c[1:]], prec=n)

def lah(k, m): return math.comb(k - 1, m - 1) * math.factorial(k) // math.factorial(m)

def lemmaC_bounds(b, kmax=4):
    """Lemma C (C5)+(C6): for 0 < t <= b <= 1/8 and 1 <= k <= 4:
       |st^(k)(t)| <= e^2 e^(-1/b) sum_m 2^m L(k,m) b^-(k+m);  0 <= st(t) <= e^2 e^(-1/b).  Returns arb list."""
    assert 0 < b <= 0.125
    B = arb(b); E = arb(1).exp() ** 2 * (-1 / B).exp()
    out = [E]
    for k in range(1, kmax + 1):
        out.append(MUT['lemmaC'] * E * sum(2 ** m * lah(k, m) * B ** (-(k + m)) for m in range(1, k + 1)))
    return out

def radial(step, r1, r2, Tlo, Thi, L=6):
    """Returns dict of arb balls valid for every T in [Tlo, Thi] (Tlo == Thi: a point):
       sigma-series coefficients of p and of the R6 radial functions (c6 = 1).
       keys: p0, p1 (p and dp/dsigma), and for R6: Bh, Gam, A11, A00, A10 (value) and their d/dT (if L >= 6
       and not an end piece), plus 'grad' flag."""
    w = r2 - r1
    Tb = tball(Tlo, Thi)
    s0 = arb(r1) + arb(w) * Tb
    sig0 = s0 * s0
    X = arb_series([sig0, 1], prec=L)
    Ts = (X.sqrt() - arb(r1)) / arb(w)
    endpiece = (step == 'cinf') and (Tlo <= 0 or Thi >= 1)
    if not endpiece:
        C = _mutate(st_of_series(Ts, step), MUT['st'])
        grad = L >= 6
    else:
        # Lemma C: st^(k) balls on the end piece, then Faa di Bruno in sigma (natural enclosure, no gradient)
        if Tlo <= 0:
            b = Thi; Bk = lemmaC_bounds(b)
            S = [arb(0).union(Bk[0])] + [arb(0, 0).union(Bk[k]).union(-Bk[k]) for k in range(1, 5)]
        else:
            b = 1 - Tlo; Bk = lemmaC_bounds(b)
            S = [(1 - Bk[0]).union(arb(1))] + [Bk[k].union(-Bk[k]) for k in range(1, 5)]
        tc = Ts.coeffs() + [arb(0)] * 6
        T1, T2, T3, T4 = tc[1], 2 * tc[2], 6 * tc[3], 24 * tc[4]
        c0 = S[0]; c1 = S[1] * T1; c2 = S[2] * T1 ** 2 + S[1] * T2
        c3 = S[3] * T1 ** 3 + 3 * S[2] * T1 * T2 + S[1] * T3
        c4 = S[4] * T1 ** 4 + 6 * S[3] * T1 ** 2 * T2 + S[2] * (3 * T2 ** 2 + 4 * T1 * T3) + S[1] * T4
        L = 5
        X = arb_series([sig0, 1], prec=L)
        C = arb_series([c0, c1, c2 / 2, c3 / 6, c4 / 24], prec=L)
        grad = False
    F = (X * X + 1).sqrt()
    C1 = C.derivative(); C2 = C1.derivative()
    Bq = 5 * C + (7 * X + 8 * X * X * X) * C1 + X * X * (1 + X * X) * C2
    p = X * Bq / F                                        # L4 (b)
    if L < 5:
        pc = p.coeffs() + [arb(0)] * 2
        return {'p0': pc[0], 'p1': pc[1] if p.prec >= 2 else None, 'dsig_dT': 2 * arb(w) * s0}
    Pi = -Bq / (2 * X * X * F)                         # Pi = -(c6/2) p / sigma^3 with c6 = 1
    Pi1 = Pi.derivative(); Pi2 = Pi1.derivative()
    a0 = X * X * Pi; a1 = X * X * X * Pi1; a2 = X * X * X * X * Pi2
    beta = X * (8 * a1 + a2) / F + (7 * a1 + a2 + 5 * a0) / (X * F)
    Gam = beta * beta / 4 - 16 * a1 * a1 - a1 * a2 + 5 * a0 * a2
    A11 = a1 * a1; A00 = a0 * a0 / 2; A10 = 2 * a1 * a0
    def co(S, k):
        c = S.coeffs(); return c[k] if k < len(c) else arb(0)
    out = {'p0': co(p, 0), 'p1': co(p, 1) if p.prec >= 2 else None}
    b0 = co(beta, 0)
    out['beta'] = b0
    out['Bh'] = abs(b0) / 2
    out['Gam'] = co(Gam, 0); out['A11'] = co(A11, 0); out['A00'] = co(A00, 0); out['A10'] = co(A10, 0)
    out['sig'] = sig0
    dsig_dT = 2 * arb(w) * s0                             # d sigma / dT
    out['dsig_dT'] = dsig_dT
    if grad and beta.prec >= 2:
        sgn_ok = (b0 > 0) or (b0 < 0)
        out['dGam'] = co(Gam, 1) * dsig_dT; out['dA11'] = co(A11, 1) * dsig_dT
        out['dA00'] = co(A00, 1) * dsig_dT; out['dA10'] = co(A10, 1) * dsig_dT
        if sgn_ok:
            out['dBh'] = (co(beta, 1) * dsig_dT / 2) * (1 if b0 > 0 else -1)
        out['grad'] = bool(sgn_ok)
    else:
        out['grad'] = False
    return out

def radial_T(step, r1, r2, Tb, N):
    """All radial quantities as Taylor series in tau, T = Tb + tau (Tb an arb ball or point inside (0,1)
       for cinf).  sigma-derivatives by the chain rule inside the series algebra: d/dsigma = (d/dtau)/(dsigma/dtau).
       Coefficient k of the result encloses f^(k)(T)/k! for every T in the ball Tb.
       Lengths: p -> N-2, R6 radial functions -> N-4."""
    w = r2 - r1
    Ts = arb_series([Tb, 1], prec=N)
    S = _mutate(st_of_series(Ts, step), MUT['st'])        # c(sigma(T)) = st(T)
    X = (arb(r1) + arb(w) * Ts) ** 2                     # sigma(T)
    Xd = X.derivative()
    Cs = S.derivative() / Xd                             # c'(sigma)
    Css = Cs.derivative() / Xd                           # c''(sigma)
    F = (X * X + 1).sqrt()
    Bq = 5 * S + (7 * X + 8 * X * X * X) * Cs + X * X * (1 + X * X) * Css
    out = {'p': X * Bq / F}                              # L4 (b)
    if N - 4 >= 1:
        Pi = -Bq / (2 * X * X * F)                       # Pi = -(c6/2) p / sigma^3, c6 = 1
        Pis = Pi.derivative() / Xd; Piss = Pis.derivative() / Xd
        a0 = X * X * Pi; a1 = X * X * X * Pis; a2 = X * X * X * X * Piss
        beta = X * (8 * a1 + a2) / F + (7 * a1 + a2 + 5 * a0) / (X * F)
        out.update({'beta': beta, 'Gam': beta * beta / 4 - 16 * a1 * a1 - a1 * a2 + 5 * a0 * a2,
                    'A11': a1 * a1, 'A00': a0 * a0 / 2, 'A10': 2 * a1 * a0})
    return out

def _co(S, k):
    c = S.coeffs()
    if k >= S.prec: raise ValueError('series too short')
    return c[k] if k < len(c) else arb(0)

def taylor_model(fc, fb, h):
    """fc: series at the centre (point), fb: series over the ball, h = half-width (arb).
       value(T)  in f0 + f1 t + f2 t^2 + f3(xi) t^3,       |t| <= h   (Lagrange remainder, K = 3)
       deriv(T)  in f1 + 2 f2 t + 3 f3(xi) t^2"""
    t = arb(0, 1) * h; t2 = arb(0).union(h * h)
    f0, f1, f2 = _co(fc, 0), _co(fc, 1), _co(fc, 2)
    R3 = _co(fb, 3)
    val = f0 + f1 * t + f2 * t2 + R3 * (arb(0, 1) * h ** 3)
    der = f1 + 2 * f2 * t + 3 * R3 * t2
    # natural enclosures from the ball series as a cross-check / tightening
    v_nat, d_nat = _co(fb, 0), _co(fb, 1)
    if not (val.overlaps(v_nat) and der.overlaps(d_nat)):
        raise FloatingPointError('taylor model and natural enclosure disjoint (bug)')
    return val.intersection(v_nat), der.intersection(d_nat)

RKEYS = ['Bh', 'Gam', 'A11', 'A00', 'A10']
def radial_record(step, r1, r2, Tlo, Thi):
    """float record for an R6 T-interval: value over interval (10) | d/dT over interval (10) | centre (10) | flag."""
    rec = np.full(31, np.nan)
    if step == 'cinf' and (Tlo <= 0 or Thi >= 1):          # Lemma C end piece: natural values only
        rb = radial(step, r1, r2, Tlo, Thi, L=6)
        rc = rb
        for i, k in enumerate(RKEYS):
            lo, hi = a2f(rb[k])
            if k in ('Bh', 'A11', 'A00'): lo = max(lo, 0.0)
            rec[2 * i], rec[2 * i + 1] = lo, hi
        Tc = 0.5 * (Tlo + Thi)
        if 0 < Tc < 1:
            rc = radial_T(step, r1, r2, arb(Tc), 5)
            vals = {'Bh': abs(_co(rc['beta'], 0)) / 2}
            for k in RKEYS[1:]: vals[k] = _co(rc[k], 0)
            for i, k in enumerate(RKEYS):
                clo, chi = a2f(vals[k])
                if k in ('Bh', 'A11', 'A00'): clo = max(clo, 0.0)
                rec[20 + 2 * i], rec[21 + 2 * i] = clo, chi
        else:
            rec[20:30] = rec[0:10]
        rec[30] = 0.0
        return rec
    Tc = 0.5 * (Tlo + Thi)                                   # dyadic grid -> exact
    h = (arb(Thi) - arb(Tlo)) / 2
    fc = radial_T(step, r1, r2, arb(Tc), 8)
    fb = radial_T(step, r1, r2, tball(Tlo, Thi), 8)
    bval, bder = taylor_model(fc['beta'], fb['beta'], h)
    sgn = 1 if bval > 0 else (-1 if bval < 0 else 0)
    vals = {'Bh': abs(bval) / 2}; ders = {'Bh': bder * sgn / 2}
    cen = {'Bh': abs(_co(fc['beta'], 0)) / 2}
    for k in RKEYS[1:]:
        vals[k], ders[k] = taylor_model(fc[k], fb[k], h)
        cen[k] = _co(fc[k], 0)
    for i, k in enumerate(RKEYS):
        lo, hi = a2f(vals[k]); clo, chi = a2f(cen[k])
        if k in ('Bh', 'A11', 'A00'): lo = max(lo, 0.0); clo = max(clo, 0.0)
        rec[2 * i], rec[2 * i + 1] = lo, hi
        rec[20 + 2 * i], rec[21 + 2 * i] = clo, chi
        if sgn != 0 or k != 'Bh':
            rec[10 + 2 * i], rec[11 + 2 * i] = a2f(ders[k])
    rec[30] = 1.0 if sgn != 0 else 0.0
    return rec

# =====================================================================================================
# 4. P* : one-variable branch and bound in arb
# =====================================================================================================
def p_enclosure(step, r1, r2, Tlo, Thi):
    if step == 'cinf' and (Tlo <= 0 or Thi >= 1):           # Lemma C end piece
        return radial(step, r1, r2, Tlo, Thi, L=4)['p0']
    Tc = 0.5 * (Tlo + Thi); h = (arb(Thi) - arb(Tlo)) / 2
    fc = radial_T(step, r1, r2, arb(Tc), 6)['p']
    fb = radial_T(step, r1, r2, tball(Tlo, Thi), 6)['p']
    return taylor_model(fc, fb, h)[0]

def p_point(step, r1, r2, T):
    if step == 'cinf' and T <= 0:
        return arb(0)                                        # c = c' = c'' = 0 at T = 0
    if step == 'cinf' and T >= 1:
        sg = arb(r2) ** 2; return 5 * sg / (sg * sg + 1).sqrt()   # c = 1, c' = c'' = 0 at T = 1
    return _co(radial_T(step, r1, r2, arb(T), 3)['p'], 0)

def pstar(step, tol_rel=1e-11, max_iter=400, log=None):
    """sup over sigma in [r1^2, r2^2] of |p|: lower = max of point evaluations, upper = max over the cover."""
    r1, r2 = WINDOWS[step]
    n0 = int(round(1 / END_B))
    ivs = [(k * END_B, (k + 1) * END_B) for k in range(n0)]
    L = arb(0); Lpt = None; Uk = arb(0); n_eval = 0; t0 = time.time()
    for it in range(max_iter):
        encs = []
        for (a, b) in ivs:
            u = abs(p_enclosure(step, r1, r2, a, b)).upper(); n_eval += 1
            m = 0.5 * (a + b)
            pv = abs(p_point(step, r1, r2, m)).lower(); n_eval += 1
            if pv > L: L, Lpt = pv, m
            encs.append((u, a, b))
        thr = (L * (1 + arb(tol_rel))).upper()
        new = []
        for (u, a, b) in encs:
            if u < L: continue
            if u <= thr:
                if u > Uk: Uk = u
                continue
            m = 0.5 * (a + b); new += [(a, m), (m, b)]
        ivs = new
        if log and (it % 10 == 0 or not ivs):
            log('  P* %s it %d: %d intervals active, L = %s' % (step, it, len(ivs), L.str(16)))
        if not ivs: break
    U = Uk if Uk > L else L
    return {'lower': L.str(20, radius=False), 'upper': U.str(20, radius=False),
            'lower_float': float(_dn(float(L))), 'upper_float': float(_up(float(U))),
            'argmax_T': Lpt, 'argmax_abs_zeta': (r1 + (r2 - r1) * Lpt) if Lpt is not None else None,
            'complete': not ivs, 'evaluations': n_eval, 'iterations': it + 1, 'seconds': time.time() - t0}

# =====================================================================================================
# 5. symmetry group G_J and the Dirichlet cone
# =====================================================================================================
def enumerate_W():
    """Aut(L), L = Lambda/2 (Hurwitz/2), as exact 4x4 matrices with entries in (1/2)Z."""
    from fractions import Fraction as Fr
    mins = [tuple(2 * s * (np.arange(4) == i)) for i in range(4) for s in (1, -1)]
    mins += [tuple(v) for v in itertools.product((1, -1), repeat=4)]
    M = np.array(mins, dtype=np.int64)                       # 4 x minimal vectors of L (|.|^2 = 4)
    B = np.array([[2, 0, 0, 0], [0, 2, 0, 0], [0, 0, 2, 0], [1, 1, 1, 1]], dtype=np.int64)   # 4 x basis
    G = B @ B.T
    idx = [np.arange(24)] * 4
    cand = [[]]
    for i in range(4):          # depth-first with Gram constraints
        nxt = []
        for c in cand:
            for j in range(24):
                if M[j] @ M[j] != G[i, i]: continue
                if all(M[c[k]] @ M[j] == G[k, i] for k in range(i)):
                    nxt.append(c + [j])
        cand = nxt
    Binv = np.linalg.inv(B.astype(float))
    mats = []
    for c in cand:
        V = M[c]                                  # rows = images of basis rows
        g = V.T.astype(float) @ np.linalg.inv(B.T.astype(float))
        g = np.round(g * 2) / 2
        assert np.allclose(g @ B.T, V.T) and np.allclose(g @ g.T, np.eye(4))
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

def cone_keep(Lm, face, ylo, yhi):
    """True where the box may meet the cone {x : <x, l> >= 0 for all l in Lm} (exact dyadic data)."""
    n = len(face); keep = np.ones(n, bool)
    Ylo, Yhi = full_y(face, ylo, yhi)
    for s in range(0, n, 20000):
        a = Ylo[s:s + 20000][:, None, :] * Lm[None]; b = Yhi[s:s + 20000][:, None, :] * Lm[None]
        mx = _up(np.maximum(a, b).sum(-1))
        keep[s:s + 20000] = np.all(mx >= 0, axis=1)
    return keep

FREE = np.array([[j for j in range(4) if j != k] for k in range(4)])
def full_y(face, ylo, yhi):
    n = len(face); k = face // 2; sg = np.where(face % 2 == 0, 1.0, -1.0)
    Ylo = np.empty((n, 4)); Yhi = np.empty((n, 4)); r = np.arange(n)
    Ylo[r, k] = sg; Yhi[r, k] = sg
    fr = FREE[k]
    for a in range(3):
        Ylo[r, fr[:, a]] = ylo[:, a]; Yhi[r, fr[:, a]] = yhi[:, a]
    return Ylo, Yhi

# =====================================================================================================
# 6. R6 box evaluator (centered form + natural extension)
# =====================================================================================================
def x_tight(Ylo, Yhi):
    """exact-range enclosure of x_i = y_i/|y| over the box (monotone in y_i, and in |y_j|)."""
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

def x_point(Y):
    n = I(Y[:, 0]) * I(Y[:, 0]) + I(Y[:, 1]) * I(Y[:, 1]) + I(Y[:, 2]) * I(Y[:, 2]) + I(Y[:, 3]) * I(Y[:, 3])
    r = irecip(isqrt(n))
    return [I(Y[:, i]) * r for i in range(4)]

def eval_boxes(rec, face, ylo, yhi, hT):
    """rec: (n,31) radial records; returns upper bounds, rigorous lower bounds at the box centres, the
       chosen split dimension (0 = T, 1..3 = free y), and the number of boxes using the centred form."""
    n = len(face)
    Ylo, Yhi = full_y(face, ylo, yhi)
    fr = FREE[face // 2]; r = np.arange(n)
    # --- Du inputs x_i over the box: value = tight range, gradient d x_i / d y_free(a)
    xt = x_tight(Ylo, Yhi)
    Yi = [I(Ylo[:, i], Yhi[:, i]) for i in range(4)]
    nn = isq(Yi[0]) + isq(Yi[1]) + isq(Yi[2]) + isq(Yi[3])
    rr = irecip(isqrt(nn)); r3 = rr * rr * rr
    xd = []
    for i in range(4):
        g = []
        for a in range(3):
            ja = fr[:, a]
            Yj = I(Ylo[r, ja], Yhi[r, ja])
            dij = (ja == i)
            gi = iwhere(dij, rr, I(np.zeros(n))) - (Yi[i] * Yj) * r3
            g.append(gi)
        xd.append(Du(xt[i], g))
    H, U, k2, w = angular_invariants(xd, dsq)
    Uv = imax0(U.v); k2v = imax0(k2.v)
    nat = lambda j: I(rec[:, 2 * j], rec[:, 2 * j + 1])
    dT = lambda j: I(np.nan_to_num(rec[:, 10 + 2 * j]), np.nan_to_num(rec[:, 11 + 2 * j]))
    cen = lambda j: I(rec[:, 20 + 2 * j], rec[:, 21 + 2 * j])
    Bh, Gm, A11, A00, A10 = [nat(j) for j in range(5)]
    H2 = isq(H.v)
    disc = H2 * Gm + Uv * A11 + k2v * A00 + w.v * A10
    if np.any(disc.hi < 0): raise FloatingPointError('disc enclosure entirely negative (impossible)')
    t1n = iabs(H.v) * Bh; t2n = isqrt(disc)
    Rn = t1n + t2n
    # --- centre values (rigorous point evaluation)
    Yc = 0.5 * (Ylo + Yhi)
    Hc, Uc, k2c, wc = angular_invariants(x_point(Yc), isq)
    Bhc, Gmc, A11c, A00c, A10c = [cen(j) for j in range(5)]
    discc = isq(Hc) * Gmc + imax0(Uc) * A11c + imax0(k2c) * A00c + wc * A10c
    t1c = iabs(Hc) * Bhc; t2c = isqrt(discc)
    Rc = t1c + t2c
    # --- gradients over the box
    hy = 0.5 * (yhi - ylo)
    hTs = I(-hT, hT); hys = [I(-hy[:, a], hy[:, a]) for a in range(3)]
    gradT_ok = rec[:, 30] > 0
    Hsign_ok = (H.v.lo > 0) | (H.v.hi < 0)
    disc_ok = disc.lo > 0
    sH = np.where(H.v.lo > 0, 1.0, -1.0)
    # term1 = |H| Bh
    ok1 = Hsign_ok & gradT_ok
    aH = H.v * sH
    g1T = aH * iwhere(gradT_ok, dT(0), I(np.zeros(n)))
    g1y = [(H.g[a] * sH) * Bh for a in range(3)]
    gy1 = g1y[0] * hys[0] + g1y[1] * hys[1] + g1y[2] * hys[2]
    t1cf = t1c + g1T * hTs + gy1
    # mixed form: natural in T (radial values over the T-interval), centred in the direction variables
    t1yc = iabs(Hc) * Bh
    t1m = t1yc + gy1
    t1 = iwhere(ok1, isect(t1n, t1cf), t1n)
    t1 = iwhere(Hsign_ok, isect(t1, t1m), t1)
    # term2 = sqrt(disc)
    ok2 = disc_ok & gradT_ok
    dsafe = iwhere(disc_ok, disc, I(np.ones(n)))
    inv2s = irecip(isqrt(dsafe) * 2.0)
    g2T = (H2 * dT(1) + Uv * dT(2) + k2v * dT(3) + w.v * dT(4)) * inv2s
    g2y = [((H.v * H.g[a]) * 2.0 * Gm + U.g[a] * A11 + k2.g[a] * A00 + w.g[a] * A10) * inv2s for a in range(3)]
    gy2 = g2y[0] * hys[0] + g2y[1] * hys[1] + g2y[2] * hys[2]
    t2cf = t2c + g2T * hTs + gy2
    t2yc = isqrt(isq(Hc) * Gm + imax0(Uc) * A11 + imax0(k2c) * A00 + wc * A10)
    t2m = t2yc + gy2
    t2 = iwhere(ok2, isect(t2n, t2cf), t2n)
    t2 = iwhere(disc_ok, isect(t2, t2m), t2)
    if np.any(t1.lo > t1.hi) or np.any(t2.lo > t2.hi):
        raise FloatingPointError('empty intersection of natural and centred enclosures (bug)')
    Rbox = isect(t1 + t2, Rn)
    check_finite('eval_boxes', Rbox, Rc)
    # --- split heuristic: smear (|gradient| x half-width) where available; for T without a gradient, the
    #     half-width of the T-only enclosure at the centre direction
    Twid = 0.5 * ((t1yc.hi - t1yc.lo) + (t2yc.hi - t2yc.lo))
    smT = np.where(ok1, g1T.mag() * hT, np.where(gradT_ok, 0.0, 0.5 * Twid)) + np.where(ok2, g2T.mag() * hT, np.where(gradT_ok, 0.0, 0.5 * Twid))
    smy = [np.where(Hsign_ok, g1y[a].mag(), 0) * hy[:, a] + np.where(disc_ok, g2y[a].mag(), 0) * hy[:, a] for a in range(3)]
    sm = np.stack([smT] + smy, axis=1)
    noinfo = ~(Hsign_ok | disc_ok)
    widths = np.stack([hT * 2.0] + [hy[:, a] for a in range(3)], axis=1)
    sm = np.where(noinfo[:, None], widths, sm)
    sm = np.where(np.isfinite(sm), sm, widths)
    dim = np.argmax(sm, axis=1)
    return Rbox.hi, Rc.lo, dim, int(np.sum(ok1 | ok2)), Rbox

# =====================================================================================================
# 7. R* branch and bound (resumable)
# =====================================================================================================
class RadCache:
    def __init__(self, step):
        self.step = step; self.r1, self.r2 = WINDOWS[step]; self.d = {}; self.n_arb = 0
    def get(self, Tlo, Thi):
        out = np.empty((len(Tlo), 31))
        keys = list(zip(Tlo.tolist(), Thi.tolist()))
        for i, k in enumerate(keys):
            v = self.d.get(k)
            if v is None:
                v = radial_record(self.step, self.r1, self.r2, k[0], k[1]); self.d[k] = v; self.n_arb += 1
            out[i] = v
        return out

def initial_boxes(Lm, nT=64, ny=8):
    Tl = np.arange(nT) / nT; Th = (np.arange(nT) + 1) / nT
    g = np.linspace(-1, 1, ny + 1)
    cells = np.array(list(itertools.product(range(ny), repeat=3)))
    ylo = g[cells]; yhi = g[cells + 1]
    faces = np.repeat(np.arange(8), len(cells)); ylo = np.tile(ylo, (8, 1)); yhi = np.tile(yhi, (8, 1))
    keep = cone_keep(Lm, faces, ylo, yhi)
    faces, ylo, yhi = faces[keep], ylo[keep], yhi[keep]
    nA = len(faces)
    TT = np.repeat(np.arange(nT), nA)
    return {'Tlo': Tl[TT], 'Thi': Th[TT], 'face': np.tile(faces, nT), 'ylo': np.tile(ylo, (nT, 1)), 'yhi': np.tile(yhi, (nT, 1))}

def rstar_run(step, Lm, state, tol_rel, budget_s, chunk=40000, max_active=3_000_000, log=print, tag=''):
    """state: dict with box arrays + L, Lpt, Ukept, it, nboxes_evaluated. Runs until done or budget."""
    t0 = time.time(); rc = RadCache(step)
    while True:
        B = state['boxes']
        nact = len(B['face'])
        if nact == 0:
            state['done'] = True; break
        if time.time() - t0 > budget_s: break
        # best-first: process the boxes with the largest previous upper bounds first
        order = np.argsort(-B['ub']) if 'ub' in B else np.arange(nact)
        take = order[:chunk]; rest = order[chunk:]
        sub = {k: v[take] for k, v in B.items()}
        hT = 0.5 * (sub['Thi'] - sub['Tlo'])
        rec = rc.get(sub['Tlo'], sub['Thi'])
        ub, lbc, dim, ncf, _ = eval_boxes(rec, sub['face'], sub['ylo'], sub['yhi'], hT)
        state['n_eval'] += len(take)
        j = int(np.argmax(lbc))
        if lbc[j] > state['L']:
            state['L'] = float(lbc[j])
            state['Lpt'] = [float(0.5 * (sub['Tlo'][j] + sub['Thi'][j])), int(sub['face'][j]),
                            (0.5 * (sub['ylo'][j] + sub['yhi'][j])).tolist()]
        L = state['L']; thr = L * (1 + tol_rel)
        alive = ub >= L
        conv = alive & (ub <= thr)
        if np.any(conv): state['Ukept'] = max(state['Ukept'], float(ub[conv].max()))
        act = alive & ~conv
        # split
        S = {k: v[act] for k, v in sub.items()}; d = dim[act]; ubA = ub[act]
        lo_kids = {k: v.copy() for k, v in S.items()}; hi_kids = {k: v.copy() for k, v in S.items()}
        isT = d == 0
        mT = 0.5 * (S['Tlo'] + S['Thi'])
        lo_kids['Thi'] = np.where(isT, mT, S['Thi']); hi_kids['Tlo'] = np.where(isT, mT, S['Tlo'])
        for a in range(3):
            sel = d == a + 1
            m = 0.5 * (S['ylo'][:, a] + S['yhi'][:, a])
            lo_kids['yhi'][sel, a] = m[sel]; hi_kids['ylo'][sel, a] = m[sel]
        kids = {k: np.concatenate([lo_kids[k], hi_kids[k]]) for k in S}
        kids['ub'] = np.concatenate([ubA, ubA])
        keep = cone_keep(Lm, kids['face'], kids['ylo'], kids['yhi'])
        kids = {k: v[keep] for k, v in kids.items()}
        R = {k: v[rest] for k, v in B.items()} if len(rest) else None
        if R is not None:
            R = {k: v[R['ub'] >= L] for k, v in R.items()}   # prune waiting boxes too
            state['boxes'] = {k: np.concatenate([R[k], kids[k]]) for k in kids}
        else:
            state['boxes'] = kids
        state['it'] += 1
        nb = len(state['boxes']['face'])
        maxub = float(max(state['boxes']['ub'].max() if nb else 0.0, state['Ukept']))
        state['U_current'] = max(maxub, L)
        if state['it'] % 5 == 0 or nb == 0:
            log('  R* %s%s it %d: active %d, L %.10g, U %.10g (gap %.2e), centred %.0f%%, arb radial %d, %.0f s' % (
                step, tag, state['it'], nb, L, state['U_current'], state['U_current'] / L - 1, 100.0 * ncf / max(1, len(take)), rc.n_arb, time.time() - t0))
        if nb > max_active:
            log('  R* %s: active set above %d, stopping this chunk' % (step, max_active)); break
    state['seconds'] += time.time() - t0
    return state

def point_value(rc, T, face, y):
    """rigorous enclosure of R6 (c6 = 1) at one point (T, direction y on the given cube face)."""
    Ta = np.array([T]); ya = np.array([y], float)
    rec = rc.get(Ta, Ta)
    _, lb, _, _, Rb = eval_boxes(rec, np.array([face]), ya, ya, np.zeros(1))
    return float(lb[0]), float(Rb.hi[0])

SEED_DIRS = [(0, (0.0, 0.0, 0.0)), (0, (1.0, 0.0, 0.0)), (0, (1.0, 1.0, 1.0))]   # axis, (1,1,0,0), (1,1,1,1)
def seed_lower(step, rc, nT=256):
    """lower bound seed: rigorous point values on a T grid along three special directions, then a
       ternary refinement in T (only a lower bound is claimed, so no optimality is needed)."""
    best = (0.0, None)
    for face, y in SEED_DIRS:
        for T in (np.arange(1, nT) / nT):
            v = point_value(rc, float(T), face, y)[0]
            if v > best[0]: best = (v, (float(T), face, list(y)))
    T0, face, y = best[1]; a, b = T0 - 1.0 / nT, T0 + 1.0 / nT
    for _ in range(60):
        m1, m2 = a + (b - a) / 3, b - (b - a) / 3
        v1 = point_value(rc, m1, face, y)[0]; v2 = point_value(rc, m2, face, y)[0]
        if v1 > best[0]: best = (v1, (m1, face, y))
        if v2 > best[0]: best = (v2, (m2, face, y))
        if v1 < v2: a = m1
        else: b = m2
    return best

def rstar_new_state(step, Lm):
    B = initial_boxes(Lm)
    B['ub'] = np.full(len(B['face']), np.inf)
    L, pt = seed_lower(step, RadCache(step))
    return {'boxes': B, 'L': L, 'Lpt': [pt[0], pt[1], list(pt[2])], 'L_seed': L, 'Ukept': 0.0, 'it': 0, 'n_eval': 0,
            'seconds': 0.0, 'done': False, 'U_current': np.inf, 'step': step}

def describe_point(step, Lpt):
    r1, r2 = WINDOWS[step]; T, face, y = Lpt
    Y = np.empty(4); k = face // 2; Y[k] = 1.0 if face % 2 == 0 else -1.0; Y[FREE[k]] = y
    return {'T': T, 'abs_zeta': r1 + (r2 - r1) * T, 'direction': (Y / np.linalg.norm(Y)).round(8).tolist()}

def save_state(path, st):
    B = st['boxes']; meta = {k: v for k, v in st.items() if k != 'boxes'}
    np.savez(path, meta=json.dumps(meta), **B)
def load_state(path):
    z = np.load(path, allow_pickle=False)
    meta = json.loads(str(z['meta']))
    meta['boxes'] = {k: z[k] for k in ('Tlo', 'Thi', 'face', 'ylo', 'yhi', 'ub')}
    return meta

# =====================================================================================================
# 8. c6 enclosures
# =====================================================================================================
def L_shell_groups(R):
    """counts of l in L, 0 < |l| < R, grouped by (m0^2, N) with l = m/4, N = |m|^2 (m all even or all odd)."""
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

def c6_direct(R=20.0):
    """(8/3) sum U6(l0/|l|)|l|^-8 over 0 < |l| < R exactly (rational), tail by Lemma A4 (p = 8) x max|U6| = 7."""
    grp = L_shell_groups(R)
    tot = fmpq(0)
    for (q0, N), cnt in grp.items():
        t2 = fmpq(q0, N)
        U6 = 64 * t2 ** 3 - 80 * t2 ** 2 + 24 * t2 - 1
        tot += cnt * U6 * fmpq(16, N) ** 4
    p = 8; T = arb(R)
    a4 = 1024 * sum(math.comb(3, j) * arb(4) ** (j - 3) * (T - arb(0.5)) ** (j + 1 - p) / (p - j - 1) for j in range(4))
    tail = 7 * a4
    P6 = arb(tot) + arb(0, 1) * tail
    return {'c6': (P6 * 8 / 3), 'tail_on_c6': (tail * 8 / 3), 'n_points': sum(grp.values()), 'R': R}

def Jb(k, R, beta):
    """Lemma A6 (e): int_R^inf r^k exp(-beta r^2) dr <= R^(k-1) e^(-beta R^2)/(2 beta - max(k-1,0)/R^2)."""
    den = 2 * beta - max(k - 1, 0) / (R * R)
    assert den > 0
    return R ** (k - 1) * (-beta * R * R).exp() / den

def c6_ewald(eta, T=3.5, K=10.0):
    """P_6[S_reg](e0) = (3/8) c6 from the Ewald form at split parameter eta:
         l = 0 real-space term: eta^4/24;  l != 0: (1/720) int_eta^inf s^3 H_6(sqrt(s) l0) e^{-s|l|^2} ds;
         reciprocal: -(1/V) sum_{k in L*\\0} e^{-pi^2|k|^2/eta} |k|^-2 (2 pi k0)^6/720.
       Tails: Lemma A6 (a),(c),(e) for the pure 6th derivative in direction e0 (one Hermite factor, h_6)."""
    eta_a = arb(eta); pi = arb.pi(); V = arb(1) / 32
    total = eta_a ** 4 / 24
    # real-space finite part, |l| < T
    grp = L_shell_groups(T)
    fact = [arb(math.factorial(i)) for i in range(8)]
    for (q0, N), cnt in grp.items():
        Rl = arb(fmpq(N, 16)); l0sq = arb(fmpq(q0, 16))
        e = (-eta_a * Rl).exp()
        def Im(m):
            return fact[m] * e * sum((eta_a * Rl) ** i / fact[i] for i in range(m + 1)) / Rl ** (m + 1)
        cl = (64 * l0sq ** 3 * Im(6) - 480 * l0sq ** 2 * Im(5) + 720 * l0sq * Im(4) - 120 * Im(3)) / 720
        total += cnt * cl
    # reciprocal finite part: k = 2n, n in Z^4, sum n even, 0 < |k| < K
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
    # tails
    hab = {0: 120, 2: 720, 4: 480, 6: 64}                 # |coefficients| of H_6
    cq = {}
    for j, a in hab.items():
        m = (6 + j) // 2
        for i in range(m + 1):
            q = 2 * i - 8
            cq[q] = cq.get(q, arb(0)) + a * fact[m] / fact[i] * eta_a ** i
    Rr = arb(T) - arb(0.5)
    assert float(Rr.lower()) ** 2 > 3.0 / eta and float(Rr.lower()) ** 2 >= 2.0 / eta   # A6 (d),(e) conditions
    real_tail = 1024 * sum(c * sum(math.comb(3, jj) * arb(4) ** (jj - 3) * Jb(q + jj, Rr, eta_a) for jj in range(4))
                           for q, c in cq.items()) / 720
    Rk = arb(K) - 2 * arb(2).sqrt(); bet = pi * pi / eta_a
    assert float((Rk * Rk).lower()) >= 4 * eta / (2 * math.pi ** 2) and float((2 * bet * Rk * Rk).lower()) > 6
    rec_tail = (2 * pi) ** 6 / (V * 720) * sum(math.comb(3, jj) * arb(2).sqrt() ** (3 - jj) * Jb(4 + jj, Rk, bet) for jj in range(4))
    tail = real_tail + rec_tail
    P6 = total + arb(0, 1) * tail
    return {'c6': P6 * 8 / 3, 'real_tail': real_tail, 'rec_tail': rec_tail, 'eta': eta, 'T': T, 'K': K,
            'n_real_groups': len(grp), 'n_recip_groups': len(keys)}

# =====================================================================================================
# 9. double-precision reference (re-typed from wp1/v2/v2_profile.py, used only in --controls)
# =====================================================================================================
def ref_module():
    import sympy as sp
    ts = sp.symbols('t')
    q5 = ts ** 5 * (126 - 420 * ts + 540 * ts ** 2 - 315 * ts ** 3 + 70 * ts ** 4)
    q5d = [sp.lambdify(ts, sp.diff(q5, ts, k), 'numpy') for k in range(5)]
    def quintic_all(t):
        t = np.asarray(t, float); inside = (t > 0) & (t < 1); res = []
        for k in range(5):
            v = np.asarray(q5d[k](np.clip(t, 0, 1)) * np.ones_like(t), float)
            v = np.where(t >= 1, 1.0, np.where(t <= 0, 0.0, v)) if k == 0 else np.where(inside, v, 0.0)
            res.append(v)
        return res
    U = sp.symbols('u'); Q = [None, U - U ** 2]
    for m in range(1, 4): Q.append(sp.expand(sp.diff(Q[m], U) * U * (1 - U)))
    Qf = [None] + [sp.lambdify(U, Q[m], 'numpy') for m in range(1, 5)]
    def cinf_all(t):
        t = np.asarray(t, float); inside = (t > 0) & (t < 1); tt = np.where(inside, t, 0.5)
        g = 1 / tt - 1 / (1 - tt)
        with np.errstate(over='ignore'): u = 1.0 / (1.0 + np.exp(g))
        gd = [None] + [((-1) ** j * math.factorial(j) * tt ** (-j - 1) - math.factorial(j) * (1 - tt) ** (-j - 1)) for j in range(1, 5)]
        D = [None] + [(-1) ** m * Qf[m](u) for m in range(1, 5)]
        u1 = D[1] * gd[1]; u2 = D[2] * gd[1] ** 2 + D[1] * gd[2]
        u3 = D[3] * gd[1] ** 3 + 3 * D[2] * gd[1] * gd[2] + D[1] * gd[3]
        u4 = D[4] * gd[1] ** 4 + 6 * D[3] * gd[1] ** 2 * gd[2] + D[2] * (3 * gd[2] ** 2 + 4 * gd[1] * gd[3]) + D[1] * gd[4]
        u0 = np.where(t >= 1, 1.0, np.where(t <= 0, 0.0, u))
        return [u0] + [np.where(inside, np.nan_to_num(v), 0.0) for v in (u1, u2, u3, u4)]
    sg = sp.symbols('sigma', positive=True); cs = sp.symbols('c0:5'); Cf = sp.Function('C')(sg)
    fh = sp.sqrt(sg ** 2 + 1)
    pexpr = (sg / fh) * (5 * Cf + (7 * sg + 8 * sg ** 3) * sp.diff(Cf, sg) + sg ** 2 * (1 + sg ** 2) * sp.diff(Cf, sg, 2))
    def subsC(e):
        for k in range(4, -1, -1): e = e.subs(sp.diff(Cf, sg, k) if k else Cf, cs[k])
        return e
    p_l = sp.lambdify((sg,) + cs, subsC(pexpr), 'numpy')
    Piexpr = -sp.Rational(1, 2) * pexpr / sg ** 3                     # c6 = 1
    Pi_l = [sp.lambdify((sg,) + cs, subsC(sp.diff(Piexpr, sg, k)), 'numpy') for k in range(3)]
    def c_sigma(sig, sf, r1, r2):
        rho = np.sqrt(sig); w = r2 - r1; T = (rho - r1) / w; s = sf(T)
        T1 = 0.5 * sig ** -0.5 / w; T2 = -0.25 * sig ** -1.5 / w; T3 = 0.375 * sig ** -2.5 / w; T4 = -0.9375 * sig ** -3.5 / w
        return (s[0], s[1] * T1, s[2] * T1 ** 2 + s[1] * T2, s[3] * T1 ** 3 + 3 * s[2] * T1 * T2 + s[1] * T3,
                s[4] * T1 ** 4 + 6 * s[3] * T1 ** 2 * T2 + s[2] * (3 * T2 ** 2 + 4 * T1 * T3) + s[1] * T4)
    z1, z2, w1, w2 = sp.symbols('z1 z2 w1 w2')
    Eh = w1 * z1 ** 5 + w2 * z2 ** 5 - 5 * w1 * z1 * z2 ** 4 - 5 * w2 * z1 ** 4 * z2
    Eb = z1 * w1 ** 5 + z2 * w2 ** 5 - 5 * z1 * w1 * w2 ** 4 - 5 * z2 * w1 ** 4 * w2
    H6c = sp.Rational(3, 16) * (Eh + Eb); Z = (z1, z2); Wv = (w1, w2); args = (z1, z2, w1, w2)
    H6_l = sp.lambdify(args, H6c, 'numpy')
    dH = [sp.lambdify(args, sp.diff(H6c, Z[i]), 'numpy') for i in range(2)]
    dbH = [sp.lambdify(args, sp.diff(H6c, Wv[j]), 'numpy') for j in range(2)]
    ddH = [[sp.lambdify(args, sp.diff(H6c, Z[i], Wv[j]), 'numpy') for j in range(2)] for i in range(2)]
    steps = {'quintic': quintic_all, 'cinf': cinf_all}
    def p_of(rho, step):
        r1, r2 = WINDOWS[step]; sig = np.asarray(rho, float) ** 2
        return p_l(sig, *c_sigma(sig, steps[step], r1, r2))
    def R6(x, step):
        r1, r2 = WINDOWS[step]; x = np.atleast_2d(x)
        a1 = x[:, 0] + 1j * x[:, 1]; a2 = x[:, 2] + 1j * x[:, 3]; b1 = np.conj(a1); b2 = np.conj(a2)
        zz = (a1, a2); zb = (b1, b2); sig = np.abs(a1) ** 2 + np.abs(a2) ** 2
        cc = c_sigma(sig, steps[step], r1, r2)
        P0, P1, P2 = [Pi_l[k](sig, *cc) for k in range(3)]
        h = H6_l(a1, a2, b1, b2)
        M = np.zeros((len(x), 2, 2), complex); G = np.zeros((len(x), 2, 2), complex)
        fhv = np.sqrt(sig ** 2 + 1); ph2 = (sig / fhv - fhv / sig) / sig
        for i in range(2):
            for j in range(2):
                M[:, i, j] = (P1 * (i == j) + P2 * zb[i] * zz[j]) * h + P1 * (zb[i] * dbH[j](a1, a2, b1, b2) + zz[j] * dH[i](a1, a2, b1, b2)) + P0 * ddH[i][j](a1, a2, b1, b2)
                G[:, i, j] = (fhv / sig) * (i == j) + ph2 * zb[i] * zz[j]
        return np.abs(np.linalg.eigvals(np.linalg.solve(G, M))).max(1)
    return {'p_of': p_of, 'R6': R6, 'steps': steps, 'c_sigma': c_sigma, 'Pi_l': Pi_l}

# =====================================================================================================
# 10. controls
# =====================================================================================================
REF_JSON = os.path.join(os.path.dirname(HERE), 'wp1', 'v2', 'v2_profile.json')

def st_all_arb(t):
    """accurate st..st'''' (arb series, 128 bits) on an array of t in (0,1); used to give the re-typed
       double reference an accurate step (its double logistic formulas lose ~1e-6 near T -> 1)."""
    t = np.atleast_1d(np.asarray(t, float)); out = np.zeros((5, len(t)))
    save = MUT['st']; MUT['st'] = 1.0
    try:
        for i, tv in enumerate(t):
            if tv <= 0: out[0, i] = 0.0; continue
            if tv >= 1: out[0, i] = 1.0; continue
            c = st_of_series(arb_series([arb(float(tv)), 1], prec=5), 'cinf').coeffs() + [arb(0)] * 5
            out[:, i] = [float(c[k].mid()) * math.factorial(k) for k in range(5)]
    finally:
        MUT['st'] = save
    return list(out)

def random_boxes(rng, n, wy_exp=(3, 9), wT_exp=(7, 14)):
    face = rng.integers(0, 8, n)
    wy = 2.0 ** -rng.integers(*wy_exp, n)
    c = np.floor(rng.uniform(-1, 1 - wy[:, None], (n, 3)) / wy[:, None]) * wy[:, None]
    wT = 2.0 ** -rng.integers(*wT_exp, n)
    Tlo = np.floor(rng.uniform(0, 1 - wT) / wT) * wT
    return face, c, c + wy[:, None], Tlo, Tlo + wT

def containment_test(step, ref, rng, n_boxes=300, n_pts=15, rel_tol=1e-9):
    """random boxes; reference double R6 (accurate step) at random points inside must lie in the enclosure."""
    r1, r2 = WINDOWS[step]
    face, ylo, yhi, Tlo, Thi = random_boxes(rng, n_boxes)
    rc = RadCache(step); rec = rc.get(Tlo, Thi)
    _, _, _, _, Rb = eval_boxes(rec, face, ylo, yhi, 0.5 * (Thi - Tlo))
    Ylo, Yhi = full_y(face, ylo, yhi); bad = 0; worst = -np.inf
    for _ in range(n_pts):
        Y = Ylo + rng.uniform(size=(n_boxes, 4)) * (Yhi - Ylo); x = Y / np.linalg.norm(Y, axis=1)[:, None]
        T = Tlo + rng.uniform(size=n_boxes) * (Thi - Tlo); rho = r1 + (r2 - r1) * T
        v = ref['R6'](rho[:, None] * x, step)
        tol = rel_tol * np.maximum(np.abs(v), 1.0)
        out = np.maximum(v - Rb.hi, Rb.lo - v) / np.maximum(np.abs(v), 1e-300)
        bad += int(np.sum((v > Rb.hi + tol) | (v < Rb.lo - tol)))
        worst = max(worst, float(out.max()))
    badp = 0
    for k in range(60):
        a, b = Tlo[k], Thi[k]
        e = p_enclosure(step, r1, r2, a, b)
        for T in rng.uniform(a, b, 5):
            pv = ref['p_of'](np.array([r1 + (r2 - r1) * T]), step)[0]
            if not (float(e.lower()) - 1e-9 * max(1, abs(pv)) <= pv <= float(e.upper()) + 1e-9 * max(1, abs(pv))):
                badp += 1
    return {'R6_points_tested': n_boxes * n_pts, 'R6_violations': bad, 'R6_max_excess_rel(negative=inside)': worst,
            'p_points_tested': 300, 'p_violations': badp}

def controls(out_path):
    import mpmath as mp
    from scipy.optimize import minimize_scalar
    t_start = time.time()
    res = {'mode': 'controls', 'python_flint': __import__('flint').__version__, 'arb_prec_bits': ctx.prec}
    def log(s): print(s, flush=True)
    ref = ref_module()
    ref_acc = ref_module(); ref_acc['steps']['cinf'] = st_all_arb    # same reference, accurate step
    refj = json.load(open(REF_JSON))
    # ---------------- C1: double-precision reproduction of sup|p| (quintic window [0.5, 1])
    REF = refj['sup_p_quintic'][0]; REFarg = refj['sup_p_quintic'][1]
    r1, r2 = WINDOWS['quintic']
    rr = np.linspace(r1, r2, 400001); pv = np.abs(ref['p_of'](rr, 'quintic')); i = int(pv.argmax())
    o = minimize_scalar(lambda r: -abs(float(ref['p_of'](np.array([r]), 'quintic')[0])), bounds=(rr[i] - 1e-3, rr[i] + 1e-3),
                        method='bounded', options={'xatol': 1e-12})
    retyped = -o.fun
    f = lambda T: -abs(float(p_point('quintic', r1, r2, float(T)).mid()))
    Tg = np.linspace(0, 1, 2001)[1:-1]; j = int(np.argmin([f(T) for T in Tg]))
    o2 = minimize_scalar(f, bounds=(Tg[j] - 1e-3, Tg[j] + 1e-3), method='bounded', options={'xatol': 1e-13})
    prod = -o2.fun
    res['C1_sup_p_quintic'] = {'reference_json': REF, 'reference_argmax_zeta': REFarg,
                               'retyped_reference_double': retyped, 'production_formula_double': prod,
                               'production_argmax_zeta': r1 + (r2 - r1) * o2.x,
                               'diff_retyped': abs(retyped - REF), 'diff_production': abs(prod - REF),
                               'PASS_1e-9': bool(abs(prod - REF) <= 1e-9 and abs(retyped - REF) <= 1e-9)}
    log('C1 sup|p| quintic: reference %.15f  production(double) %.15f  diff %.2e  -> %s' % (
        REF, prod, abs(prod - REF), res['C1_sup_p_quintic']['PASS_1e-9']))
    # ---------------- C2: closed-form R6 (invariants) vs eigen-solver reference
    rng = np.random.default_rng(11); c2 = {}
    for step in ('quintic', 'cinf'):
        r1, r2 = WINDOWS[step]; worst = {'retyped': 0.0, 'retyped_accurate_step': 0.0}
        for _ in range(300):
            T = float(rng.uniform(0.01, 0.99)); x = rng.normal(size=4); x /= np.linalg.norm(x)
            rd = radial_T(step, r1, r2, arb(T), 5)
            Bh = abs(float(_co(rd['beta'], 0).mid())) / 2
            Gm, A11, A00, A10 = [float(_co(rd[k], 0).mid()) for k in RKEYS[1:]]
            H, U, k2, w = angular_invariants(list(x), fsq)
            mine = abs(H) * Bh + math.sqrt(max(0.0, H * H * Gm + U * A11 + k2 * A00 + w * A10))
            rho = r1 + (r2 - r1) * T
            for nm, rf in (('retyped', ref), ('retyped_accurate_step', ref_acc)):
                worst[nm] = max(worst[nm], abs(mine / rf['R6'](rho * x[None], step)[0] - 1))
        c2[step] = worst
    ax = lambda rho: 15476.9366 * ref['R6'](np.array([[rho, 0, 0, 0]]), 'quintic')[0]
    c2['quintic_axis_0.75362_x_c6doc'] = ax(0.75362); c2['quintic_axis_0.755_x_c6doc'] = ax(0.755)
    c2['doc_values'] = {'0.75362': 2.63343e6, '0.755': 2.63261e6}
    res['C2_R6_formula'] = c2
    log('C2 R6 closed form vs eigen-solver: %s' % json.dumps({k: v for k, v in c2.items() if k in ('quintic', 'cinf')}))
    log('   quintic axis x c6(doc): %.6g at 0.75362 (doc 2.63343e6), %.6g at 0.755 (doc 2.63261e6)' % (
        c2['quintic_axis_0.75362_x_c6doc'], c2['quintic_axis_0.755_x_c6doc']))
    # ---------------- C3: symmetry group, cone, invariance
    Wg, GJ, hol = group_GJ(); Lm = cone_normals(GJ)
    orbit = {tuple(g @ P_GEN) for g in GJ}
    inv_GJ = 0.0; inv_other = 0.0
    others = [g for g in Wg if not any(np.array_equal(g, h) for h in GJ)]
    for _ in range(40):
        T = float(rng.uniform(0.2, 0.9)); x = rng.normal(size=4); x /= np.linalg.norm(x)
        rho = 0.5 + 0.5 * T; base = ref['R6'](rho * x[None], 'quintic')[0]
        vals = ref['R6'](rho * np.array([g @ x for g in GJ]), 'quintic'); inv_GJ = max(inv_GJ, float(np.abs(vals / base - 1).max()))
        vo = ref['R6'](rho * np.array([g @ x for g in others[:96]]), 'quintic'); inv_other = max(inv_other, float(np.abs(vo / base - 1).max()))
    xs = rng.normal(size=(3000, 4)); xs /= np.linalg.norm(xs, axis=1)[:, None]
    covered = 0
    for x in xs:
        imgs = np.array([g @ x for g in GJ])
        if np.any(np.all(imgs @ Lm.T >= -1e-12, axis=1)): covered += 1
    B0 = initial_boxes(Lm)
    res['C3_symmetry'] = {'W_order': len(Wg), 'GJ_order': len(GJ), 'GJ_holomorphic': hol, 'orbit_of_p': len(orbit),
                          'cone_normals': int(len(Lm)), 'R6_max_rel_change_under_GJ': inv_GJ,
                          'R6_max_rel_change_under_W_minus_GJ': inv_other, 'random_dirs_with_image_in_cone': '%d/3000' % covered,
                          'initial_direction_boxes_kept': int(len(B0['face']) // 64), 'of': 4096}
    log('C3 |W| = %d, |G_J| = %d (%d hol), orbit %d; R6 change under G_J %.1e, under W minus G_J %.2e; cone coverage %d/3000; %d of 4096 direction boxes kept' % (
        len(Wg), len(GJ), hol, len(orbit), inv_GJ, inv_other, covered, len(B0['face']) // 64))
    # ---------------- C4: mutation tests
    mp.mp.dps = 60
    hst = lambda x: mp.exp(-1 / x) / (mp.exp(-1 / x) + mp.exp(-1 / (1 - x)))
    tpts = [0.02, 0.05, 0.1, 0.2, 0.3, 0.45, 0.5, 0.55, 0.7, 0.835, 0.9, 0.95, 0.98]
    mpv = {t: [mp.diff(hst, mp.mpf(t), k) for k in range(5)] for t in tpts}
    c4 = {'logistic_path': {}, 'lemmaC_endpiece': {}}
    for fac in (1.0, 0.5):
        viol = []
        MUT['st'] = fac
        try:
            for t in tpts:
                c = _mutate(st_of_series(arb_series([arb(t).union(arb(t + 1e-12)), 1], prec=5), 'cinf'), MUT['st']).coeffs()
                for k in range(5):
                    enc = c[k] * math.factorial(k)
                    tv = arb(mp.nstr(mpv[t][k], 45))
                    if not enc.overlaps(tv): viol.append((t, k))
        finally:
            MUT['st'] = 1.0
        c4['logistic_path']['factor_%g' % fac] = {'violations': len(viol), 'tested': 5 * len(tpts), 'first': viol[:6]}
    ratios = []
    for b in (1 / 64, 1 / 32, 1 / 16, 1 / 8):
        Bk = lemmaC_bounds(b)
        for t in (b, 0.9 * b, 0.75 * b, 0.5 * b, 0.25 * b):
            for k in range(1, 5):
                tv = abs(mp.diff(hst, mp.mpf(t), k))
                if tv > 0: ratios.append(float(tv / mp.mpf(Bk[k].upper().str(40, radius=False))))
    worst_ratio = max(ratios)
    for fac in (1.0, 0.5, 0.25, 0.125):
        c4['lemmaC_endpiece']['factor_%g' % fac] = {'violations': int(sum(r > fac for r in ratios)), 'tested': len(ratios)}
    c4['lemmaC_endpiece']['max_true_over_bound'] = worst_ratio
    c4['lemmaC_endpiece']['note'] = ('Lemma C (C5) has slack >= %.2f on end pieces b <= 1/8, so halving the bound cannot be '
                                     'detected by any point comparison; the detector fires for factors below %.3f' % (1 / worst_ratio, worst_ratio))
    rngc = np.random.default_rng(5)
    c4['pipeline_unmutated'] = {s: containment_test(s, ref_acc, rngc) for s in ('quintic', 'cinf')}
    MUT['st'] = 0.5
    try:
        c4['pipeline_mutated_st_derivs_x0.5'] = {s: containment_test(s, ref_acc, np.random.default_rng(5)) for s in ('quintic', 'cinf')}
    finally:
        MUT['st'] = 1.0
    res['C4_mutation'] = c4
    log('C4 logistic path: %s' % json.dumps(c4['logistic_path']))
    log('   Lemma C end piece: %s' % json.dumps(c4['lemmaC_endpiece']))
    log('   pipeline unmutated: %s' % json.dumps(c4['pipeline_unmutated']))
    log('   pipeline mutated (st derivatives x 0.5): %s' % json.dumps(c4['pipeline_mutated_st_derivs_x0.5']))
    # ---------------- C5: NaN / inf guards
    c5 = {}
    def raises(fn):
        try:
            fn(); return False
        except (FloatingPointError, ValueError, AssertionError):
            return True
    c5['check_finite_on_nan'] = raises(lambda: check_finite('t', I(np.array([np.nan]), np.array([1.0]))))
    c5['check_finite_on_inf'] = raises(lambda: check_finite('t', I(np.array([0.0]), np.array([np.inf]))))
    c5['a2f_on_nonfinite_arb'] = raises(lambda: a2f(arb('nan')))
    rc = RadCache('quintic'); recq = rc.get(np.array([0.5]), np.array([0.5 + 2 ** -10]))
    c5['eval_boxes_nan_box'] = raises(lambda: eval_boxes(recq, np.array([0]), np.array([[np.nan, 0, 0]]), np.array([[0.1, 0.1, 0.1]]), np.array([2 ** -11])))
    rbad = recq.copy(); rbad[0, 2] = np.nan
    c5['eval_boxes_nan_radial'] = raises(lambda: eval_boxes(rbad, np.array([0]), np.array([[0, 0, 0]]), np.array([[0.1, 0.1, 0.1]]), np.array([2 ** -11])))
    c5['lemmaC_endpiece_too_wide_refused'] = raises(lambda: lemmaC_bounds(0.2))
    c5['all_PASS'] = all(c5.values())
    res['C5_nan_inf'] = c5
    log('C5 NaN/inf guards: %s' % json.dumps(c5))
    # ---------------- C6: timing on the CONTROL window (quintic, [0.5, 1])
    c6t = {}
    t0 = time.time(); pq = pstar('quintic', tol_rel=1e-11); c6t['P_quintic_certified_seconds'] = time.time() - t0
    c6t['P_quintic_certified'] = {k: pq[k] for k in ('lower', 'upper', 'argmax_abs_zeta', 'complete', 'evaluations')}
    log('C6 P* quintic (control window) certified in %.1f s: [%s, %s] at |zeta| = %.6f' % (
        c6t['P_quintic_certified_seconds'], pq['lower'], pq['upper'], pq['argmax_abs_zeta']))
    t0 = time.time(); st = rstar_new_state('quintic', Lm); c6t['R_seed_seconds'] = time.time() - t0
    st = rstar_run('quintic', Lm, st, tol_rel=1e-4, budget_s=300, log=log, tag=' (reduced cover, tol 1e-4)')
    c6t['R_quintic_reduced'] = {'tol_rel': 1e-4, 'done': st['done'], 'lower_c6=1': st['L'], 'upper_c6=1': st['U_current'],
                                'seed_lower': st['L_seed'], 'iterations': st['it'], 'boxes_evaluated': st['n_eval'],
                                'seconds': st['seconds'] + c6t['R_seed_seconds'], 'argmax': describe_point('quintic', st['Lpt']),
                                'times_c6_doc_15476.9366': [st['L'] * 15476.9366, st['U_current'] * 15476.9366]}
    log('C6 R* quintic reduced cover: done %s, [%.10g, %.10g] (c6 = 1) x c6doc = [%.7g, %.7g], %d boxes, %.0f s' % (
        st['done'], st['L'], st['U_current'], st['L'] * 15476.9366, st['U_current'] * 15476.9366, st['n_eval'],
        c6t['R_quintic_reduced']['seconds']))
    t0 = time.time(); c6_ewald(6.0, T=2.5, K=7.0); c6t['c6_ewald_reduced_seconds'] = time.time() - t0
    t0 = time.time(); c6_direct(R=8.0); c6t['c6_direct_R8_seconds'] = time.time() - t0
    c6t['c6_direct_R20_estimate_seconds'] = c6t['c6_direct_R8_seconds'] * (20 / 8) ** 4
    res['C6_timing'] = c6t
    try:
        import psutil; res['peak_rss_GB_process'] = psutil.Process().memory_info().peak_wset / 1e9
    except Exception:
        res['peak_rss_GB_process'] = None
    res['seconds_total'] = time.time() - t_start
    res['ALL_CONTROLS_PASS'] = bool(res['C1_sup_p_quintic']['PASS_1e-9'] and c5['all_PASS']
                                    and c4['logistic_path']['factor_1']['violations'] == 0
                                    and c4['logistic_path']['factor_0.5']['violations'] > 0
                                    and all(v['R6_violations'] == 0 and v['p_violations'] == 0 for v in c4['pipeline_unmutated'].values())
                                    and all(v['R6_violations'] > 0 and v['p_violations'] > 0 for v in c4['pipeline_mutated_st_derivs_x0.5'].values())
                                    and c4['lemmaC_endpiece']['factor_1']['violations'] == 0
                                    and res['C3_symmetry']['GJ_order'] == 192 and res['C3_symmetry']['orbit_of_p'] == 192)
    json.dump(res, open(out_path, 'w'), indent=1, default=str)
    log('controls written to %s (%.0f s); ALL_CONTROLS_PASS = %s' % (out_path, res['seconds_total'], res['ALL_CONTROLS_PASS']))

# =====================================================================================================
# 11. main (certified run; needs GO-WP2a; chunked and resumable)
# =====================================================================================================
R_TOL = 1e-8; P_TOL = 1e-12

def main(out_path, budget_s):
    if not os.path.exists(GO_FILE):
        print('REFUSED: %s does not exist (the lead creates it after reading the controls).' % GO_FILE)
        sys.exit(2)
    os.makedirs(STATE_DIR, exist_ok=True)
    t0 = time.time()
    def log(s): print(s, flush=True)
    rpath = os.path.join(STATE_DIR, 'main_results.json')
    R = json.load(open(rpath)) if os.path.exists(rpath) else {}
    def save(): json.dump(R, open(rpath, 'w'), indent=1, default=str)
    if 'c6' not in R:
        e6 = c6_ewald(6.0); e8 = c6_ewald(8.0); d = c6_direct(20.0)
        ok = e6['c6'].overlaps(e8['c6']) and e6['c6'].overlaps(d['c6']) and e8['c6'].overlaps(d['c6'])
        if not ok: raise FloatingPointError('c6 enclosures disagree')
        c6 = e6['c6'].intersection(e8['c6'])
        R['c6'] = {'eta6': e6['c6'].str(30), 'eta8': e8['c6'].str(30), 'direct_R20': d['c6'].str(15),
                   'eta6_tails': [e6['real_tail'].str(5), e6['rec_tail'].str(5)], 'eta8_tails': [e8['real_tail'].str(5), e8['rec_tail'].str(5)],
                   'certified': c6.str(30), 'lower': c6.lower().str(30, radius=False), 'upper': c6.upper().str(30, radius=False),
                   'all_three_overlap': True}
        save(); log('c6 = %s' % R['c6']['certified'])
    for step in ('quintic', 'cinf'):
        k = 'P_' + step
        if k not in R:
            R[k] = pstar(step, tol_rel=P_TOL, log=log); save(); log('%s: [%s, %s]' % (k, R[k]['lower'], R[k]['upper']))
    Wg, GJ, hol = group_GJ(); Lm = cone_normals(GJ)
    for step in ('quintic', 'cinf'):
        k = 'R_' + step
        if k in R and R[k].get('done'): continue
        sp = os.path.join(STATE_DIR, 'R_%s.npz' % step)
        st = load_state(sp) if os.path.exists(sp) else rstar_new_state(step, Lm)
        while not st['done'] and time.time() - t0 < budget_s:
            st = rstar_run(step, Lm, st, tol_rel=R_TOL, budget_s=min(60.0, budget_s - (time.time() - t0)), log=log)
            save_state(sp, st)
        R[k] = {'done': st['done'], 'lower_c6=1': st['L'],
                'upper_c6=1': max(st['U_current'], st['L']) if st['done'] else st['U_current'],
                'tol_rel': R_TOL, 'iterations': st['it'], 'boxes_evaluated': st['n_eval'], 'seconds': st['seconds'],
                'argmax': describe_point(step, st['Lpt'])}
        save()
        if not st['done']: break
    complete = all(('R_' + s) in R and R['R_' + s]['done'] for s in ('quintic', 'cinf'))
    out = {'status': 'COMPLETE' if complete else 'INCOMPLETE (rerun the same command to resume)', 'windows': WINDOWS,
           'arb_prec_bits': ctx.prec, 'results': R}
    if complete:
        c6lo = arb(R['c6']['lower']); c6hi = arb(R['c6']['upper'])
        for step in ('quintic', 'cinf'):
            P = R['P_' + step]; Rr = R['R_' + step]
            Plo, Phi = arb(P['lower']), arb(P['upper'])
            Rlo = (arb(Rr['lower_c6=1']) * c6lo).lower(); Rhi = (arb(Rr['upper_c6=1']) * c6hi).upper()
            m5lo = Plo if Plo > 5 else arb(5); m5hi = Phi if Phi > 5 else arb(5)
            F8lo = (3 * c6lo / 16 * m5lo).lower(); F8hi = (3 * c6hi / 16 * m5hi).upper()
            out[step] = {'P_star_window': [P['lower'], P['upper']], 'P_star_argmax_abs_zeta': P['argmax_abs_zeta'],
                         'sup_all_sigma_of_abs_p': 'max(5, P*) (5 is the limit sigma -> inf, not attained)',
                         'R_star': [Rlo.str(12, radius=False), Rhi.str(12, radius=False)],
                         'R_star_argmax': Rr['argmax'],
                         'sup_abs_F8': [F8lo.str(12, radius=False), F8hi.str(12, radius=False)]}
    json.dump(out, open(out_path, 'w'), indent=1, default=str)
    log('main: %s -> %s (%.0f s this chunk)' % (out['status'], out_path, time.time() - t0))
    sys.exit(0 if complete else 3)

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--controls', metavar='OUT_JSON'); g.add_argument('--main', metavar='OUT_JSON')
    ap.add_argument('--budget', type=float, default=540.0, help='seconds per --main invocation (capped at 570)')
    a = ap.parse_args()
    if a.controls:
        controls(a.controls)
    else:
        main(a.main, min(a.budget, 570.0))
