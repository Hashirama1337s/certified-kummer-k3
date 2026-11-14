#!/usr/bin/env python
"""s2.py -- Stage 2 re-assembly of the certified Kummer K3 bound on (0, 0.03] with two improved ingredients.

Sources: stage2/STAGE2-DESIGN.md (Lemma L5', Lemma Z', Lemma G'), ANALYSIS-CHAIN-v2.2.md (Theorem 2, Prop C0,
Lemma G 3.2-3.4, continuity 7), wp5/wp5.py + BUILD-WP5.md (Stage B engine).  wp5.py is imported READ-ONLY (bytecode
writing disabled); only the diameter function (L5 -> L5') and the Lemma G constant (Doeblin Gamma -> Gamma') are
replaced.  theorem2(), load_inputs(), containing(), AInt and the arb helpers are wp5's own, unchanged.

Inputs: ONLY the existing certified files wp2/result_wp2a.json (via wp5.load_inputs), wp2/result_wp2b.json,
wp3/result_wp3.json and wp5/result_wp5.json (Stage A values C_core_F, C_core_R, e_core per a-interval).

(1) L5':  D(a, e; rho) = sqrt(2(1+e)) [2 rho + (1 + a^4/rho^4)^(1/4) (1/sqrt2 + 2 (pi - 2) rho)],  0 < rho <= 1/sqrt32,
    any fixed rho is an upper bound for diam(X, g_R); rho per interval from a double optimiser, then evaluated in arb.
(2) Lemma G' with Lemma Z':  Gamma'(K, D) = D^2 Gh'(k = K D^2),
    Gh'(k) = int_0^inf phi,  phi(t) = min(1 - max_{s<=t} c0(s), inf_{t0<=t} B(t0) e^{-lam (t - t0)}, 1),
    lam = pi^2 / cosh^3(sqrt(k/12))  (Lemma Z' in units D = 1),
    B(t0) = (1 + tau/t0)^(2 alpha) exp(2 sqrt2 alpha k tau/(alpha - 1) + alpha/(4 tau)) - 1 at one (alpha, tau) per t0.
    Certified as a left Riemann sum of a pointwise upper bound of the nonincreasing function 1 - m(t), plus the
    spectral tail S(t_N)/lam beyond the last grid point, S(t_N) = B(t0*) e^{-lam (t_N - t0*)} the spectral term itself
    (never phi(t_N) when that is set by the Doeblin term or by 1; REFEREE-S2 item 1).
    Monotone in k: all choices (s, tau, alpha, t0, switch) are made at k_hi = upper(K D^2) and evaluated at k_hi with
    lam = Lemma Z'(k_hi); at fixed choices the value is nondecreasing in k, so it bounds Gh'(k) on [0, k_hi]
    (REFEREE-S2 item 2: never a larger-K B with a smaller-K lam).
(3) Theorem 2 on the 30 WP5 a-intervals with L5' and Gamma'.

Modes
  --controls OUT.json   Z1-Z4, guards, control-resolution Stage B on 6 intervals with float re-evaluation.
  --main OUT.json       refuses unless stage2/GO-S2 exists; resumable (stage2/state_s2/), exit 3 = rerun.
Run:  python <lab-tools>/memcap.py 8 py -3.11 <lab>/stage2/s2.py --controls <lab>/stage2/controls_s2.json
"""
import os, sys, json, time, math, argparse, hashlib
sys.dont_write_bytecode = True
import numpy as np
from flint import arb, ctx

K3 = os.environ.get('K3_ROOT') or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(K3, 'wp5'))
import wp5 as W                                      # sealed, read-only (imports wp2b read-only as well)
ctx.prec = 160
from wp5 import Refuse, AInt, fin, up, U, nonneg, ssqrt, frac

GO_FILE = os.path.join(HERE, 'GO-S2')
STATE_DIR = os.path.join(HERE, 'state_s2')
WP5_RESULT = K3 + '/wp5/result_wp5.json'
A_MAX = '0.05'                                       # DEFINITIONS global hypothesis a <= 0.05 (L5' has no range of its own)
RHO_MAX_SQ_INV = 32                                  # rho <= 1/sqrt(32)  <=>  32 rho^2 <= 1
PI = math.pi

def sha_file(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()

# =====================================================================================================
# 1. Lemma L5' (STAGE2-DESIGN a.2; PROVED there, NOT REFEREED)
# =====================================================================================================
def D_float(a, e, rho):
    return math.sqrt(2 * (1 + e)) * (2 * rho + (1 + a ** 4 / rho ** 4) ** 0.25 * (1 / math.sqrt(2) + 2 * (PI - 2) * rho))

def rho_choice(a_hi, e):
    """double-precision optimiser of D(a_hi, e; rho) over 0 < rho <= 1/sqrt32 (any admissible rho is valid)."""
    if a_hi <= 0:
        return 1e-9                                  # a = 0: D -> sqrt(1+e) as rho -> 0; any small rho is admissible
    from scipy.optimize import minimize_scalar
    rmax = 1 / math.sqrt(32) * (1 - 1e-12)
    r = minimize_scalar(lambda x: D_float(a_hi, e, x), bounds=(a_hi * 1e-3, rmax), method='bounded',
                        options=dict(xatol=1e-12 * max(a_hi, 1e-6)))
    rho = float(r.x)
    # guard the optimiser with a coarse log grid (choose the better of the two)
    g = np.geomspace(a_hi * 1e-2, rmax, 4001)
    vg = np.array([D_float(a_hi, e, x) for x in g])
    if vg.min() < D_float(a_hi, e, rho):
        rho = float(g[int(np.argmin(vg))])
    return rho

def D_L5p(ai, e, rho):
    """L5': diam(X, g_R) <= sqrt(2(1+e)) [2 rho + (1 + a^4/rho^4)^(1/4) (1/sqrt2 + 2(pi-2) rho)] for every a in ai,
    with rho a fixed exact number, 0 < rho <= 1/sqrt32, e < 1, a <= 0.05.  For fixed rho every term is nonnegative
    and increasing in a, and the evaluation is the arb inclusion over the a-ball (no monotonicity of the inf used)."""
    e = arb(e) if not isinstance(e, arb) else e
    fin(e, 'e')
    if not (e.upper() < 1):
        raise Refuse("L5' needs e < 1 (positivity of omega_a); e = %s" % e)
    if not (e.lower() >= 0):
        raise Refuse("L5' needs e >= 0")
    if not (isinstance(rho, float) and math.isfinite(rho) and rho > 0):
        raise Refuse("L5' needs a finite rho > 0")
    ra = arb(rho)                                    # rho is an exact double
    if not ((RHO_MAX_SQ_INV * ra * ra).upper() <= 1):
        raise Refuse("L5' needs rho <= 1/sqrt(32)")
    if frac(ai.shi) > frac(A_MAX):
        raise Refuse('a <= 0.05 required (DEFINITIONS global hypothesis)')
    a4 = ai.p(4)
    fac = (1 + a4 / ra ** 4).sqrt().sqrt()
    return (2 * (1 + e)).sqrt() * (2 * ra + fac * (1 / arb(2).sqrt() + 2 * (arb.pi() - 2) * ra))

# =====================================================================================================
# 2. Lemma Z' and Lemma G' (STAGE2-DESIGN b.2-b.3; PROVED there, NOT REFEREED)
# =====================================================================================================
S_GRID = W._S_GRID                                   # the WP5 s-grid (2001 points) for c0
TAU = np.geomspace(1e-4, 20.0, 6001)                 # tau grid of the double optimiser of B(t0)
R2F = math.sqrt(2.0)

def lam_hat_arb(k):
    """Lemma Z' in units D = 1: lambda1 D^2 >= pi^2 / cosh^3(sqrt(k/12)) (arb ball; use its lower end)."""
    ka = nonneg(arb(k))
    return arb.pi() ** 2 / ((ka / 12).sqrt().cosh() ** 3)

def _c0_float(tt, k):
    s = S_GRID[None, :]
    P = 1 / (4 * tt[:, None] * (1 - s)) - 2 * np.log(s)
    Q = 2 * R2F * k * tt[:, None] * (1 - s)
    val = (np.sqrt(P) + np.sqrt(Q)) ** 2
    jj = np.argmin(val, axis=1)
    return np.exp(-val[np.arange(len(tt)), jj]), jj

def _B_float(t0s, k):
    """per t0: (B value, tau, alpha) of the double optimiser; alpha = 1 at k = 0 (Li-Yau Cor 2.1 limit)."""
    out = []
    for t0 in t0s:
        g = 2 * np.log1p(TAU / t0) + 1 / (4 * TAU)
        if k == 0:
            j = int(np.argmin(g)); out.append((math.expm1(g[j]), float(TAU[j]), 1.0)); continue
        c = 2 * R2F * k * TAU
        v = (np.sqrt(g) + np.sqrt(c)) ** 2
        j = int(np.argmin(v))
        al = 1 + math.sqrt(c[j] / g[j])
        if al - 1 < 1e-12:
            al = 1 + 1e-12
        vv = al * g[j] + c[j] * al / (al - 1)
        out.append((math.expm1(vv), float(TAU[j]), float(al)))
    return out

def _plan(k, lam, h, n, sub, spectral):
    """double-precision plan: c0 running max, B(t0) on the sub-grid, choice of bound per grid point."""
    tt = np.arange(0, n + 1) * h
    c0f = np.zeros(n + 1); jj = np.zeros(n + 1, int)
    c0f[1:], jj[1:] = _c0_float(tt[1:], k)
    runidx = np.zeros(n + 1, int)
    best = 0
    for i in range(1, n + 1):
        if c0f[i] >= c0f[best]:
            best = i
        runidx[i] = best
    plan = dict(tt=tt, c0f=c0f, jj=jj, runidx=runidx)
    if not spectral:
        return plan
    t0i = np.arange(sub, n + 1, sub)
    Bl = _B_float(tt[t0i], k)
    Bf = np.array([b[0] for b in Bl])
    logB = np.log(Bf)
    # second bound at t_i: min over t0_j <= t_i of log B_j - lam (t_i - t0_j)
    M = logB[None, :] + lam * tt[t0i][None, :] - lam * tt[:, None]
    M[tt[:, None] < tt[t0i][None, :] - 1e-15] = np.inf
    jb = np.argmin(M, axis=1)
    b2 = np.exp(M[np.arange(n + 1), jb])
    b1 = 1 - c0f[runidx]
    b1[0] = 1.0
    phi = np.minimum(np.minimum(b1, b2), 1.0)
    tail = b2[n] / lam
    plan.update(t0i=t0i, Bl=Bl, jb=jb, b1=b1, b2=b2, phi=phi, float_total=float(h * phi[:n].sum() + tail))
    return plan

_GHP = {}

def ghp_cert(k, h=1e-3, Tend=2.0, sub=5, lam_override=None, spectral=True, Bfun=None, want_float_only=False):
    """certified upper bound for Gh'(k') for every 0 <= k' <= k (D = 1), Lemma G' + Lemma Z'.
    spectral=False: Doeblin tail instead (reproduces wp5 gh_cert; control Z3a).
    lam_override: an exact lower bound for lambda1 D^2 (controls Z2/Z4); Bfun: mutation hook (control Z2c).
    Returns dict(U=arb upper point, float_value, ...)."""
    kf = float(k)
    if not (math.isfinite(kf) and kf >= 0):
        raise Refuse('ghp_cert: bad k')
    n = int(round(Tend / h))
    key = (kf, h, n, sub, None if lam_override is None else float(lam_override), spectral, Bfun is None, want_float_only)
    if Bfun is None and key in _GHP:
        return _GHP[key]
    if lam_override is not None:
        lam_a = arb(float(lam_override))
        if not (math.isfinite(float(lam_override)) and float(lam_override) > 0):
            raise Refuse('ghp_cert: lambda1 must be > 0')
    else:
        lam_a = arb(lam_hat_arb(kf).lower())         # exact number <= the Lemma Z' bound
    lam_f = float(lam_a)
    if not (lam_f > 0 and math.isfinite(lam_f)):
        raise Refuse('ghp_cert: lambda1 must be > 0')
    if not spectral:
        # Doeblin form of Lemma G (wp5 gh_cert, own implementation): sum_{i<j} h (1 - c0(t_i)) + T (1 - c0(T))/c0(T)
        P = _plan(kf, lam_f, h, n, sub, False)
        c0f = P['c0f']; tt = P['tt']
        left = np.concatenate([[0.0], np.cumsum(h * (1 - c0f[:-1]))])
        with np.errstate(divide='ignore'):
            tot = left[1:] + tt[1:] * (1 - c0f[1:]) / c0f[1:]
        j = int(np.argmin(tot)) + 1
        ha = arb(h); S = ha
        for i in range(1, j):
            S += ha * (1 - _c0lo(i, ha, kf, P['jj']))
        cT = _c0lo(j, ha, kf, P['jj'])
        if not (cT > 0):
            raise Refuse('c0(T) not positive')
        Ub = U(S + arb(j) * ha * (1 - cT) / cT)
        res = dict(U=Ub, float_value=float(tot[j - 1]), T=float(j * h), k=kf, mode='doeblin')
        _GHP[key] = res
        return res
    P = _plan(kf, lam_f, h, n, sub, True)
    if Bfun is not None:                              # mutation: replace the on-diagonal bound (control only)
        tt = P['tt']; t0i = P['t0i']
        Bm = np.array([max(Bfun(t), 1e-300) for t in tt[t0i]])
        M = np.log(Bm)[None, :] + lam_f * tt[t0i][None, :] - lam_f * tt[:, None]
        M[tt[:, None] < tt[t0i][None, :] - 1e-15] = np.inf
        jb = np.argmin(M, axis=1); b2 = np.exp(M[np.arange(n + 1), jb])
        b2[np.array([Bfun(tt[t0i][j]) for j in jb]) <= 0] = 0.0
        phi = np.minimum(np.minimum(P['b1'], b2), 1.0)
        return dict(float_value=float(h * phi[:n].sum() + b2[n] / lam_f), mode='mutation')
    if want_float_only:
        res = dict(float_value=P['float_total'], k=kf, lam=lam_f)
        _GHP[key] = res
        return res
    ha = arb(h)
    ka = arb(kf)
    tt = P['tt']; t0i = P['t0i']; Bl = P['Bl']
    Bup = {}
    def Bcert(j):
        if j not in Bup:
            t0 = arb(int(t0i[j])) * ha
            _, tau, al = Bl[j]
            ta = arb(tau); aa = arb(al)
            if kf == 0 and al == 1.0:
                v = (1 + ta / t0) ** 2 * (1 / (4 * ta)).exp() - 1
            else:
                if not (aa > 1):
                    raise Refuse('alpha must be > 1 for k > 0')
                v = (2 * aa * (1 + ta / t0).log() + 2 * arb(2).sqrt() * aa * ka * ta / (aa - 1) + aa / (4 * ta)).exp() - 1
            Bup[j] = U(nonneg(v))
        return Bup[j]
    c0c = {}
    def c0c_get(i):
        if i not in c0c:
            c0c[i] = _c0lo(i, ha, kf, P['jj'])
        return c0c[i]
    S = arb(0)
    n1 = n2 = n3 = 0
    for i in range(n):
        if i == 0:
            v = arb(1); n3 += 1
        elif P['b1'][i] <= P['b2'][i]:
            if P['b1'][i] >= 1:
                v = arb(1); n3 += 1
            else:
                v = U(1 - c0c_get(int(P['runidx'][i]))); n1 += 1
        else:
            j = int(P['jb'][i])
            if not (int(t0i[j]) <= i):
                raise Refuse('spectral bound needs t0 <= t')
            v = U(Bcert(j) * (-(lam_a * (arb(i) - arb(int(t0i[j]))) * ha)).exp()); n2 += 1
            if v.upper() > 1:
                v = arb(1)
        S += v
    S = S * ha
    jN = int(P['jb'][n])
    tail = U(Bcert(jN) * (-(lam_a * (arb(n) - arb(int(t0i[jN]))) * ha)).exp() / lam_a)
    Ub = U(S + tail)
    fin(Ub, 'Gh_prime')
    res = dict(U=Ub, float_value=P['float_total'], k=kf, lam_lower=lam_f, h=h, Tend=n * h, sub=sub,
               n_c0=n1, n_spec=n2, n_one=n3, tail=up(tail), t0_last=float(tt[t0i[jN]]),
               switch_t=float(tt[int(np.argmax(P['b2'] < P['b1']))]))
    _GHP[key] = res
    return res

def _c0lo(i, ha, kf, jj):
    """c0 lower bound at t_i = i h (one s per t, alpha optimised exactly; ANALYSIS-CHAIN 3.2-3.3 (ii))."""
    t = arb(i) * ha
    sv = arb(float(S_GRID[jj[i]]))
    Pv = 1 / (4 * t * (1 - sv)) - 2 * sv.log()
    Qv = 2 * arb(2).sqrt() * arb(kf) * t * (1 - sv)
    v = (Pv.sqrt() + nonneg(Qv).sqrt()) ** 2
    return arb((-v).exp().lower())

def make_gamma_prime(h=1e-3, Tend=2.0, sub=5):
    """drop-in replacement of wp5.make_gamma: Gamma'(K, D) = D^2 Gh'_cert(upper(K D^2)) (Lemma G' + Lemma Z')."""
    def G(K, D):
        fin(K, 'K'); fin(D, 'D')
        if not (D.lower() > 0):
            raise Refuse("Gamma': D must be > 0")
        k = up(nonneg(K) * D * D)
        g = ghp_cert(k, h=h, Tend=Tend, sub=sub)
        v = U(D * D * g['U'])
        return v, dict(mode='prime', k_up=k, Ghp_up=up(g['U']), lam_hat_lower=g['lam_lower'], D2Ghp=up(v),
                       Gexpl=up(W.gamma_expl(K, D)))
    return G

def gamma_prime_float(K, D, h=1e-3, Tend=2.0, sub=5):
    k = K * D * D
    return D * D * ghp_cert(k, h=h, Tend=Tend, sub=sub, want_float_only=True)['float_value']

# =====================================================================================================
# 3. Stage B with L5' and Gamma'; independent double-precision Theorem 2 (unscaled formulas of ANALYSIS-CHAIN)
# =====================================================================================================
def load_stageA():
    d = json.load(open(WP5_RESULT))
    if d.get('status') != 'COMPLETE':
        raise Refuse('wp5 result incomplete')
    return d['stageA'], d

def interval_inputs(inp, SA, ai):
    rb = W.containing(inp['wp2b']['rows'], ai)
    p3 = W.containing(inp['wp3']['pieces'], ai)
    Cq = max(SA['C_core_F'], p3['C_out'])
    CR = max(SA['C_core_R'], p3['R_out'])
    e = max(rb['e'], p3['e_out'], inp['wp3']['e_out'], SA['e_core'])
    return dict(Cq=Cq, CR=CR, c0s=rb['c0s'], e=e, wp2b_row=list(rb['a']), wp3_piece=list(p3['a']),
                C_core_F=SA['C_core_F'], C_out=p3['C_out'], C_core_R=SA['C_core_R'], R_out=p3['R_out'])

def theorem2_float(a, Cq, CR, c0s, D, gfun, Lam0=2.0, iters=6, Kfac=2.0):
    """double precision, unscaled, written directly from ANALYSIS-CHAIN 2, 4, 5.2, 6, 7, 9 at one point a > 0.
    Same conventions as the certificate: continuity step with K_t = Kfac Lam0 kappa, then iters t = 1 steps."""
    eta = -2 * math.log1p(-Cq * a ** 8)
    eta1 = 2 * (-math.expm1(-eta / 2)); eta2 = 2 * math.expm1(eta / 2)
    kappa = CR * a ** 6; sigma = 2 * kappa; C0 = c0s / a ** 2
    delta = 4 * (-math.expm1(-eta)) + sigma / C0 * math.exp(-eta)
    den = 1 - math.expm1(eta) / 2
    def step(Lam, K2):
        G1 = gfun(kappa, D); G2 = gfun(K2, math.sqrt(Lam) * D)
        o = (2 * eta1 * G1 + 2 * eta2 * G2) / den
        x = 1 + math.sqrt(delta / (4 * C0 * o))
        L = eta + math.log((1 + x + math.sqrt((x - 1) ** 2 + x * delta)) / (2 * x)) + C0 * (1 + x) * o
        eps2 = math.expm1(L)
        m1 = eps2 + math.sqrt(eps2 * (2 + eps2) + (-math.expm1(-eta)))    # M2 - 1
        Uv = (1 + m1) * math.exp(eta)
        d2 = math.expm1(eta) + m1 * math.exp(eta)
        return Uv, d2, G1, G2
    Uc, _, _, _ = step(Lam0, Kfac * Lam0 * kappa)
    if not Uc < Lam0:
        return dict(U=Uc, closes=False)
    Lam = Lam0
    for _ in range(iters):
        Uv, d2, G1, G2 = step(Lam, 0.0)
        Lam = Uv
    return dict(U=Uc, closes=True, d2=d2, d2_over_a3=d2 / a ** 3, G1=G1, G2=G2)

def stageB_s2(inp, SA, ai, gamma, gopts):
    I = interval_inputs(inp, SA, ai)
    a_hi = float(ai.shi)
    rho = rho_choice(a_hi, I['e'])
    Dball = D_L5p(ai, arb(I['e']), rho)
    D = U(Dball)
    t = time.time()
    r = W.theorem2(ai, I['Cq'], I['CR'], I['c0s'], D, gamma)
    r['seconds'] = time.time() - t
    I.update(rho=rho, D=up(D), D_float_at_hi=D_float(a_hi, I['e'], rho),
             D_old_L5=(up(W.D_L5(ai, arb(I['e']))) if frac(ai.shi) <= frac('0.03') else None))
    r['input_sources'] = I
    # float re-evaluation of Theorem 2 at a = a_hi with the same inputs (certified must be >= float)
    if a_hi > 0 and r['closes']:
        gf = lambda K, Dd: gamma_prime_float(K, Dd, **gopts)
        fr = theorem2_float(a_hi, I['Cq'], I['CR'], I['c0s'], D_float(a_hi, I['e'], rho), gf)
        r['float_recheck'] = dict(a=a_hi, U_float=fr['U'], d2_float=fr.get('d2'), d2_over_a3_float=fr.get('d2_over_a3'),
                                  cert_U_ge_float=bool(r['smallness'] >= fr['U'] * (1 - 1e-12)),
                                  cert_d2_ge_float=bool(fr['closes'] and r['d2'] >= fr['d2'] * (1 - 1e-12)),
                                  ratio_d2=(r['d2'] / fr['d2']) if fr.get('d2') else None,
                                  ratio_U=r['smallness'] / fr['U'])
    return r

def leading_s2(inp, SAmap, gamma):
    """limsup d2/a^3 at a = 0 (the a = 0 value of the scaled bound, rho = 1e-9) and the certified sup over (0, 1e-4]."""
    first = AInt('0', '0.001')
    I = interval_inputs(inp, SAmap[first.key()], first)
    out = {}
    for tag, ai, rho in (('a=0', AInt('0', '0'), 1e-9), ('(0,1e-4]', AInt('0', '0.0001'), rho_choice(1e-4, I['e']))):
        D = U(D_L5p(ai, arb(I['e']), rho))
        r = W.theorem2(ai, I['Cq'], I['CR'], I['c0s'], D, gamma)
        out[tag] = dict(d2_over_a3=r['d2_over_a3'], Lambda_seq=r['Lambda_seq'], rho=rho, D=up(D), Cq=I['Cq'],
                        c0s=I['c0s'], Gamma1=r['final']['G1'], Gamma2=r['final']['G2'], smallness=r['smallness'],
                        g1info=r['final']['g1info'])
    return out

def a_intervals30():
    return [AInt(W._dec(frac('0.001') * i), W._dec(frac('0.001') * (i + 1))) for i in range(30)]

# =====================================================================================================
# 4. Control Z1: crossed-ball brute force for L5' step (3)
# =====================================================================================================
def lattice_L(R=1.5):
    """L = Lambda/2, Lambda = Z^4 u (Z^4 + (1/2,1/2,1/2,1/2)) (DEFINITIONS 1.1), points within R of (1/4,1/4,1/4,1/4)."""
    r = np.arange(-6, 8) / 2.0
    g = np.array(np.meshgrid(r, r, r, r, indexing='ij')).reshape(4, -1).T
    pts = np.concatenate([g, g + 0.25])
    c = np.full(4, 0.25)
    return pts[np.linalg.norm(pts - c, axis=1) <= R]

def crossing_counts(P, Dd, ell, Lp, rho):
    """for segments p + t d, t in [0, ell]: number of l with dist(l, segment) < rho (geometric), and number with
    perpendicular distance < rho and axial coordinate in (0, ell) (the design's criterion); endpoint-in-ball flags."""
    w = Lp[None, :, :] - P[:, None, :]
    tax = np.einsum('nij,nj->ni', w, Dd)
    w2 = np.einsum('nij,nij->ni', w, w)
    perp2 = w2 - tax ** 2
    tc = np.clip(tax, 0, ell[:, None])
    dist2 = w2 - 2 * tc * tax + tc ** 2
    geo = (dist2 < rho ** 2).sum(1)
    axial = ((perp2 < rho ** 2) & (tax > 0) & (tax < ell[:, None])).sum(1)
    Q = P + ell[:, None] * Dd
    dq = np.linalg.norm(Lp[None, :, :] - Q[:, None, :], axis=2).min(1)
    dp = np.sqrt(w2.min(1))
    return geo, axial, (dp >= rho) & (dq >= rho)

def control_Z1(log, n_target=100000, seed=11):
    rng = np.random.default_rng(seed)
    Lp = lattice_L()
    minv = Lp[np.abs(np.linalg.norm(Lp - 0, axis=1) - 0.5) < 1e-12]
    out = {}
    ok = True
    for rho in (0.04, 0.06, 0.10, 1 / math.sqrt(32)):
        mx_geo = mx_ax = 0; acc = 0; hist = np.zeros(8, int); att = 0
        while acc < n_target:
            m = 2000
            P = rng.random((m, 4)) * 0.5
            Dd = rng.normal(size=(m, 4)); Dd /= np.linalg.norm(Dd, axis=1)[:, None]
            ell = np.where(rng.random(m) < 0.5, 1 / math.sqrt(2), rng.random(m) / math.sqrt(2))
            geo, ax, good = crossing_counts(P, Dd, ell, Lp, rho)
            att += m
            geo, ax = geo[good], ax[good]
            acc += int(good.sum())
            mx_geo = max(mx_geo, int(geo.max())); mx_ax = max(mx_ax, int(ax.max()))
            hist += np.bincount(np.minimum(geo, 7), minlength=8)
        # adversarial: through 0 and a minimal vector v (|v| = 1/2), start just outside B(0, rho), length 1/sqrt2
        m = 20000
        v = minv[rng.integers(0, len(minv), m)]
        dv = v / np.linalg.norm(v, axis=1)[:, None]
        pert = rng.normal(size=(m, 4)) * 0.02
        Dd = dv + pert; Dd /= np.linalg.norm(Dd, axis=1)[:, None]
        off = rng.normal(size=(m, 4)); off -= np.einsum('ni,ni->n', off, Dd)[:, None] * Dd
        off /= np.linalg.norm(off, axis=1)[:, None]
        off *= (rng.random(m) * 0.9 * rho)[:, None]
        back = rho * (1.0001 + 0.5 * rng.random(m))
        P = off - back[:, None] * Dd
        ell = np.full(m, 1 / math.sqrt(2))
        # translate P into the box frame: shift by a lattice vector of (1/2)Z^4 is a symmetry of L; Lp covers |l - c| <= 1.5
        geo, ax, good = crossing_counts(P, Dd, ell, Lp, rho)
        adv_max = int(geo[good].max()) if good.any() else -1
        adv_two = int((geo[good] == 2).sum())
        passed = (mx_geo <= 2 and mx_ax <= 2 and adv_max <= 2)
        attained2 = bool(hist[2] > 0 or adv_two > 0)
        ok = ok and passed and attained2
        out['%.6f' % rho] = dict(rho=rho, random_accepted=acc, random_attempted=att, max_geometric=mx_geo, max_axial=mx_ax,
                                 hist_geometric=hist.tolist(), adversarial_n=int(good.sum()), adversarial_max=adv_max,
                                 adversarial_count2=adv_two, N_le_2=bool(passed), count2_attained=attained2)
        log('  Z1 rho=%.5f: %d random segments (endpoints outside balls): max crossed %d (geometric), %d (axial criterion); hist %s ; adversarial %d: max %d, count 2 seen %d times'
            % (rho, acc, mx_geo, mx_ax, hist.tolist()[:4], int(good.sum()), adv_max, adv_two))
    return dict(rows=out, n_lattice_points=int(len(Lp)), n_minimal_vectors=int(len(minv)), pass_=bool(ok))

def diagnostic_rho022(log, n_target=100000, seed=13):
    """instrument check (NOT a sealed control): at rho = 0.22 > 1/sqrt32 the count must reach 3, so the counter can fail."""
    rng = np.random.default_rng(seed)
    Lp = lattice_L()
    out = {}
    for rho in (0.22, 0.235):
        mx_u = mx_c = 0; n_u = n_c = 0; hu = np.zeros(8, int); hc = np.zeros(8, int)
        while n_u < n_target:
            m = 2000
            P = rng.random((m, 4)) * 0.5
            Dd = rng.normal(size=(m, 4)); Dd /= np.linalg.norm(Dd, axis=1)[:, None]
            ell = np.full(m, 1 / math.sqrt(2))
            geo, ax, good = crossing_counts(P, Dd, ell, Lp, rho)
            n_u += m; mx_u = max(mx_u, int(geo.max())); hu += np.bincount(np.minimum(geo, 7), minlength=8)
            n_c += int(good.sum())
            if good.any():
                mx_c = max(mx_c, int(geo[good].max())); hc += np.bincount(np.minimum(geo[good], 7), minlength=8)
        out['%.3f' % rho] = dict(rho=rho, n_all=n_u, max_all_segments=mx_u, hist_all=hu.tolist(), n_endpoints_outside=n_c,
                                 max_endpoints_outside=mx_c, hist_endpoints_outside=hc.tolist())
        log("  DIAG rho=%.3f (outside L5' range), length 1/sqrt2: all %d segments: max %d hist %s ; endpoints outside balls (%d): max %d hist %s"
            % (rho, n_u, mx_u, hu.tolist()[:5], n_c, mx_c, hc.tolist()[:5]))
    out['count3_reached_at_0.22_all_segments'] = bool(out['0.220']['max_all_segments'] >= 3)
    out['count3_reached_at_0.22_endpoints_outside'] = bool(out['0.220']['max_endpoints_outside'] >= 3)
    out['count3_reached_at_0.235_endpoints_outside'] = bool(out['0.235']['max_endpoints_outside'] >= 3)
    return out

# =====================================================================================================
# 5. controls
# =====================================================================================================
def logger(path):
    fh = open(path, 'a', encoding='utf-8')
    def log(msg):
        line = time.strftime('%H:%M:%S ') + msg
        print(line, flush=True); fh.write(line + '\n'); fh.flush()
    return log

def theta_torus_minus1(t):
    """V H(z,z,t) - 1 on R^4/Z^4 (V = 1): theta(t)^4 - 1, theta(t) = sum_k exp(-4 pi^2 k^2 t)."""
    k = np.arange(1, 60)
    th = 1 + 2 * np.exp(-4 * PI ** 2 * np.outer(np.atleast_1d(t), k ** 2)).sum(1)
    return th ** 4 - 1

CTRL_G = dict(h=2e-3, Tend=2.0, sub=5)
MAIN_G = dict(h=1e-3, Tend=2.0, sub=5)
CTRL_INTERVALS = ['[0, 0.001]', '[0.002, 0.003]', '[0.009, 0.01]', '[0.019, 0.02]', '[0.028, 0.029]', '[0.029, 0.03]']

def controls(out_path):
    T0 = time.time()
    log = logger(os.path.splitext(out_path)[0] + '.log')
    me = sha_file(os.path.abspath(__file__))
    log('=== s2 controls === sha256(s2.py) = %s' % me)
    res = dict(sha256_s2=me, sha256_wp5=sha_file(K3 + '/wp5/wp5.py'), ctrl_grid=CTRL_G, main_grid=MAIN_G)
    inp = W.load_inputs()
    SAmap, wp5res = load_stageA()
    res['inputs'] = {k: dict(source=inp[k]['source'], seal=inp[k]['seal'], sha256=inp[k]['sha256']) for k in inp}
    res['inputs']['wp5_result'] = dict(source=WP5_RESULT, sha256=sha_file(WP5_RESULT), sha256_wp5_in_result=wp5res['sha256_wp5'])
    log('  inputs: %s' % {k: v['sha256'][:12] for k, v in res['inputs'].items()})
    # ---- Z1
    log('-- Z1 crossed-ball brute force (L5\' step 3)')
    t = time.time(); res['Z1'] = control_Z1(log); res['Z1']['seconds'] = time.time() - t
    res['DIAG_rho022'] = diagnostic_rho022(log)          # reported separately, not part of ALL_CONTROLS_PASS
    # ---- Z3a: Doeblin mode reproduces Gammahat(0) <= 4.530841 (wp5 grid h = 1.6/4000)
    log('-- Z3a spectral tail off (Doeblin tail, wp5 grid) reproduces Gammahat(0)')
    t = time.time()
    z3a = ghp_cert(0.0, h=1.6 / 4000, Tend=1.6, spectral=False)
    ref = 4.530840314434004
    rel = abs(up(z3a['U']) - ref) / ref
    res['Z3a'] = dict(Gh0_doeblin_own=up(z3a['U']), T=z3a['T'], wp5_Gh0_upper=ref, rel=rel, pass_=bool(rel <= 1e-5 and up(z3a['U']) <= 4.530841))
    log('  Z3a: own Doeblin Gammahat(0) <= %.9f (T = %.4f) vs wp5 %.9f: rel %.2e : %s [%.1fs]'
        % (up(z3a['U']), z3a['T'], ref, rel, res['Z3a']['pass_'], time.time() - t))
    # ---- Z3b: old L5 + old Gamma reproduce wp5 d2 on [0.029, 0.03]
    log('-- Z3b old L5 + old Gamma (wp5 functions) on [0.029, 0.03] from result_wp5.json Stage A')
    t = time.time()
    ai = AInt('0.029', '0.03')
    Gold = W.make_gamma('cert', N=4000)
    rold = W.stageB_interval(inp, SAmap[ai.key()], ai, Gold)
    refB = wp5res['stageB'][ai.key()]
    rel_d2 = abs(rold['d2'] - refB['d2']) / refB['d2']; rel_U = abs(rold['smallness'] - refB['smallness']) / refB['smallness']
    res['Z3b'] = dict(d2=rold['d2'], d2_wp5=refB['d2'], rel_d2=rel_d2, U=rold['smallness'], U_wp5=refB['smallness'], rel_U=rel_U,
                      pass_=bool(rel_d2 <= 1e-6 and rel_U <= 1e-6))
    log('  Z3b: d2 = %.12g vs wp5 %.12g (rel %.1e); U = %.10f vs %.10f (rel %.1e) : %s [%.1fs]'
        % (rold['d2'], refB['d2'], rel_d2, rold['smallness'], refB['smallness'], rel_U, res['Z3b']['pass_'], time.time() - t))
    # ---- Z2: flat torus R^4/Z^4 (V = 1, D = 1, K = 0): bound >= 0.0702 (sharp)
    log('-- Z2 flat torus R^4/Z^4 (V = 1, D = 1, K = 0)')
    t = time.time()
    g0c = ghp_cert(0.0, **CTRL_G); g0m = ghp_cert(0.0, **MAIN_G)
    # (b) torus-exact inputs (theta^4 - 1 on the diagonal, lambda1 = 4 pi^2), float, must still be >= 0.0702
    tt = np.arange(0, 20001) * 1e-4
    c0t, _ = _c0_float(tt[1:], 0.0)
    b1 = np.concatenate([[1.0], 1 - np.maximum.accumulate(c0t)])
    b2 = theta_torus_minus1(tt)
    exact_in = float(1e-4 * np.minimum(np.minimum(b1, b2), 1)[:-1].sum())
    # (c) mutation: Euclidean (4 pi t)^-2 - 1 (a LOWER bound of V H - 1) in place of B -> must fall below 0.0702
    mut = ghp_cert(0.0, lam_override=4 * PI ** 2, Bfun=lambda s: (4 * PI * s) ** -2 - 1, **MAIN_G)['float_value']
    z2 = dict(Ghp0_ctrl=up(g0c['U']), Ghp0_main=up(g0m['U']), float_main=g0m['float_value'], switch_t=g0m['switch_t'],
              torus_exact_inputs_float=exact_in, mutation_euclidean_float=mut, sharp=0.0702)
    z2['pass_'] = bool(z2['Ghp0_main'] >= 0.0702 and z2['Ghp0_ctrl'] >= 0.0702 and exact_in >= 0.0702 and mut < 0.0702)
    res['Z2'] = z2
    log("  Z2: Gh'(0) <= %.6f (ctrl grid), %.6f (main grid; float %.6f; switch to spectral bound at t = %.3f) >= 0.0702 ; torus-exact inputs (float) %.5f >= 0.0702 ; mutation Euclidean on-diagonal (float) %.5f < 0.0702 detected : %s [%.1fs]"
        % (z2['Ghp0_ctrl'], z2['Ghp0_main'], z2['float_main'], z2['switch_t'], exact_in, mut, z2['pass_'], time.time() - t))
    # ---- Z4: halving lambda1 raises Gh'(0); also Gh' grid table vs Doeblin
    log("-- Z4 halving lambda1 raises Gh'(0); Gh' vs Doeblin on a k-table")
    t = time.time()
    lam_lo = float(arb(lam_hat_arb(0.0).lower()))
    half = ghp_cert(0.0, lam_override=lam_lo / 2, **MAIN_G)
    quarter = ghp_cert(0.0, lam_override=lam_lo / 4, **MAIN_G)
    z4 = dict(Ghp0=up(g0m['U']), Ghp0_half_lambda=up(half['U']), Ghp0_quarter_lambda=up(quarter['U']))
    table = {}
    for kk in (1e-4, 1e-3, 1e-2, 3e-2, 0.1, 0.3):
        gp = ghp_cert(kk, **MAIN_G)
        gd = W.gh_cert(kk, N=2000)
        table[repr(kk)] = dict(Ghp=up(gp['U']), Ghp_float=gp['float_value'], lam_hat=gp['lam_lower'], Gh_doeblin_wp5=up(gd['U']),
                               prime_le_doeblin=bool(up(gp['U']) <= up(gd['U'])), cert_ge_float=bool(up(gp['U']) >= gp['float_value'] * (1 - 1e-12)), rel_cert_minus_float=(up(gp['U']) - gp['float_value']) / gp['float_value'])
        log("  k=%g: Gh'_cert <= %.6f (float %.6f, lam_hat >= %.6f) ; Doeblin (wp5 gh_cert N=2000) %.6f" % (kk, up(gp['U']), gp['float_value'], gp['lam_lower'], up(gd['U'])))
    mono = [table[repr(kk)]['Ghp'] for kk in (1e-4, 1e-3, 1e-2, 3e-2, 0.1, 0.3)]
    z4['k_table'] = table
    z4['monotone_in_k'] = bool(all(x <= y for x, y in zip([z4['Ghp0']] + mono[:-1], mono)))
    z4['pass_'] = bool(z4['Ghp0_half_lambda'] > z4['Ghp0'] and z4['Ghp0_quarter_lambda'] > z4['Ghp0_half_lambda']
                       and z4['monotone_in_k'] and all(v['prime_le_doeblin'] and v['cert_ge_float'] for v in table.values()))
    res['Z4'] = z4
    log("  Z4: Gh'(0) %.6f -> lambda/2: %.6f -> lambda/4: %.6f ; monotone in k %s : %s [%.1fs]"
        % (z4['Ghp0'], z4['Ghp0_half_lambda'], z4['Ghp0_quarter_lambda'], z4['monotone_in_k'], z4['pass_'], time.time() - t))
    # ---- grid sensitivity (tolerance justification)
    t = time.time()
    fine = ghp_cert(0.0, h=5e-4, Tend=2.0, sub=5)
    long_ = ghp_cert(0.0, h=1e-3, Tend=3.0, sub=5)
    sub2 = ghp_cert(0.0, h=1e-3, Tend=2.0, sub=2)
    res['grid_sensitivity'] = dict(h2e3=up(g0c['U']), h1e3=up(g0m['U']), h5e4=up(fine['U']), Tend3=up(long_['U']), sub2=up(sub2['U']))
    log("  grid: Gh'(0) h=2e-3 %.6f, h=1e-3 %.6f, h=5e-4 %.6f ; Tend=3 %.6f ; t0 sub-grid 2 %.6f [%.1fs]"
        % (up(g0c['U']), up(g0m['U']), up(fine['U']), up(long_['U']), up(sub2['U']), time.time() - t))
    # ---- D(a) table (L5')
    dt = {}
    for a in (1e-4, 1e-3, 0.01, 0.02, 0.03):
        e = 1.3682e-4
        rho = rho_choice(a, e)
        Dv = up(D_L5p(AInt(repr(a), repr(a)), arb(e), rho))
        dt[repr(a)] = dict(rho=rho, D=Dv, D_old=up(W.D_L5(AInt(repr(a), repr(a)), arb(e))))
        log("  L5' a=%g: rho = %.5f, D <= %.6f (old L5 %.6f)" % (a, rho, Dv, dt[repr(a)]['D_old']))
    res['D_table'] = dt
    # ---- control-resolution Stage B on 6 intervals + float re-evaluation
    log('-- Stage B (control grid h=%g) on 6 of the 30 intervals, float re-evaluation at a_hi' % CTRL_G['h'])
    t = time.time()
    G = make_gamma_prime(**CTRL_G)
    B = {}
    fl_ok = True
    for key in CTRL_INTERVALS:
        lo, hi = key.strip('[]').split(', ')
        ai = AInt(lo, hi)
        r = stageB_s2(inp, SAmap[key], ai, G, CTRL_G)
        B[key] = r
        fr = r.get('float_recheck', {})
        fl_ok = fl_ok and r['closes'] and fr.get('cert_U_ge_float', False) and fr.get('cert_d2_ge_float', False)
        log("  B %s: rho %.5f D %.6f (old %.6f) | U %.6f | d2 <= %.6e, d2/a^3 <= %.3f | float at a_hi: U %.6f d2 %.6e (cert/float %.6f) | G1 %.4f G2 %.4f [%.1fs]"
            % (key, r['input_sources']['rho'], r['input_sources']['D'], r['input_sources']['D_old_L5'], r['smallness'], r['d2'],
               r['d2_over_a3'], fr.get('U_float', float('nan')), fr.get('d2_float', float('nan')), fr.get('ratio_d2', float('nan')),
               r['final']['G1'], r['final']['G2'], r['seconds']))
    lead = leading_s2(inp, SAmap, G)
    for tag, v in lead.items():
        log("  leading %s: d2/a^3 <= %.4f (rho %.3g, D %.7f, Gamma1 %.6f, Gamma2 %.6f)" % (tag, v['d2_over_a3'], v['rho'], v['D'], v['Gamma1'], v['Gamma2']))
    res['stageB_ctrl'] = B; res['leading_ctrl'] = lead
    res['float_recheck_pass'] = bool(fl_ok)
    res['stageB_ctrl_seconds'] = time.time() - t
    # closed-form check of the leading coefficient: 8 D0 sqrt(Cq Gh'(0)) (c0s -> 2)
    I0 = interval_inputs(inp, SAmap['[0, 0.001]'], AInt('0', '0.001'))
    cf = 8 * lead['a=0']['D'] * math.sqrt(I0['Cq'] * up(g0c['U']) * I0['c0s'] / 2)
    res['leading_closed_form_ctrl'] = cf
    log('  closed form 8 D0 sqrt(Cq Gh\'(0) c0s/2) = %.4f vs certified a=0 value %.4f' % (cf, lead['a=0']['d2_over_a3']))
    # ---- guards
    log('-- NaN/inf guards')
    guards = {}
    def expect(name, fn, exc=(FloatingPointError, Refuse, ValueError, AssertionError, ZeroDivisionError)):
        try:
            fn(); guards[name] = False
        except exc:
            guards[name] = True
    expect('nan_arb', lambda: fin(arb('nan')))
    expect('inf_float', lambda: fin(float('inf')))
    expect('rho_gt_rho_max', lambda: D_L5p(AInt('0.01', '0.02'), arb(0), 1 / math.sqrt(32) * (1 + 1e-9)))
    expect('rho_zero', lambda: D_L5p(AInt('0.01', '0.02'), arb(0), 0.0))
    expect('rho_nan', lambda: D_L5p(AInt('0.01', '0.02'), arb(0), float('nan')))
    expect('a_gt_0.05', lambda: D_L5p(AInt('0.05', '0.051'), arb(0), 0.05))
    expect('e_eq_1', lambda: D_L5p(AInt('0.01', '0.02'), arb(1), 0.03))
    expect('e_nan', lambda: D_L5p(AInt('0.01', '0.02'), arb('nan'), 0.03))
    expect('negative_k', lambda: ghp_cert(-1.0))
    expect('nan_k', lambda: ghp_cert(float('nan')))
    expect('lambda_nonpositive', lambda: ghp_cert(0.0, lam_override=0.0))
    expect('nan_K_in_Gamma_prime', lambda: make_gamma_prime(**CTRL_G)(arb('nan'), arb(1)))
    expect('inf_D_in_Gamma_prime', lambda: make_gamma_prime(**CTRL_G)(arb(0), arb('inf')))
    expect('nan_input_theorem2', lambda: W.theorem2(AInt('0.01', '0.02'), arb(float('nan')), arb(1e6), arb(2), arb(1.2), G))
    expect('non_finite_output_check', lambda: W.check_all_finite({'x': float('nan')}))
    res['guards'] = guards
    log('  guards: %s' % guards)
    W.check_all_finite(W.jsonable(res))
    res['ALL_CONTROLS_PASS'] = bool(res['Z1']['pass_'] and res['Z2']['pass_'] and res['Z3a']['pass_'] and res['Z3b']['pass_']
                                    and res['Z4']['pass_'] and res['float_recheck_pass'] and all(guards.values()))
    res['seconds'] = time.time() - T0
    log('ALL_CONTROLS_PASS = %s  (%.1f s)' % (res['ALL_CONTROLS_PASS'], res['seconds']))
    json.dump(W.jsonable(res), open(out_path, 'w'), indent=1)
    return res

# =====================================================================================================
# 6. main (refuses without GO-S2; resumable)
# =====================================================================================================
def main(out_path, budget_s=520):
    if not os.path.exists(GO_FILE):
        print('refusing: %s does not exist (only the lead creates it)' % GO_FILE)
        sys.exit(2)
    budget_s = min(float(budget_s), 520.0)
    T0 = time.time()
    os.makedirs(STATE_DIR, exist_ok=True)
    log = logger(os.path.join(HERE, 'run_s2.log'))
    st_path = os.path.join(STATE_DIR, 'main_state.json')
    cfg = dict(grid=MAIN_G, n=30, width='0.001')
    state = json.load(open(st_path)) if os.path.exists(st_path) else dict(B={}, cfg=cfg)
    if state.get('cfg') != cfg:
        raise Refuse('state written with a different configuration')
    me = sha_file(os.path.abspath(__file__))
    log('=== s2 main (resume: %d B intervals done) sha256 %s' % (len(state['B']), me))
    inp = W.load_inputs()
    SAmap, wp5res = load_stageA()
    G = make_gamma_prime(**MAIN_G)
    ais = a_intervals30()
    for ai in ais:
        if ai.key() in state['B']:
            continue
        r = stageB_s2(inp, SAmap[ai.key()], ai, G, MAIN_G)
        state['B'][ai.key()] = W.jsonable(r)
        fr = r.get('float_recheck', {})
        log('  B a in %s: rho %.5f D %.6f | closes %s U %.6f | d2 <= %.6e d2/a^3 <= %.3f | float d2 %.6e cert>=float %s'
            % (ai.key(), r['input_sources']['rho'], r['input_sources']['D'], r['closes'], r.get('smallness', float('nan')),
               r.get('d2', float('nan')), r.get('d2_over_a3', float('nan')), fr.get('d2_float', float('nan')), fr.get('cert_d2_ge_float')))
        json.dump(state, open(st_path, 'w'))
        if time.time() - T0 > budget_s:
            log('budget reached; rerun the same command (exit 3)'); sys.exit(3)
    lead = leading_s2(inp, SAmap, G)
    B = state['B']
    for tag, v in lead.items():
        log('  leading %s: d2/a^3 <= %.4f' % (tag, v['d2_over_a3']))
    g0 = ghp_cert(0.0, **MAIN_G)
    ok_float = all(b.get('float_recheck', {}).get('cert_d2_ge_float') and b.get('float_recheck', {}).get('cert_U_ge_float') for b in B.values())
    out = dict(status='COMPLETE', sha256_s2=me, sha256_wp5=sha_file(K3 + '/wp5/wp5.py'), config=cfg,
               inputs={k: dict(source=inp[k]['source'], seal=inp[k]['seal'], sha256=inp[k]['sha256']) for k in inp},
               wp5_result=dict(source=WP5_RESULT, sha256=sha_file(WP5_RESULT)),
               Ghp0_upper=up(g0['U']), stageB=B, leading=W.jsonable(lead),
               summary=dict(all_close=all(b['closes'] for b in B.values()),
                            smallness_max=max(b.get('smallness', 9) for b in B.values()),
                            U_0p03=B['[0.029, 0.03]'].get('smallness'), d2_0p03=B['[0.029, 0.03]'].get('d2'),
                            d2_0p01=B['[0.009, 0.01]'].get('d2'), d2_over_a3_0p03=B['[0.029, 0.03]'].get('d2_over_a3'),
                            limsup_d2_over_a3=lead['a=0']['d2_over_a3'], sup_d2_over_a3_0_1em4=lead['(0,1e-4]']['d2_over_a3'],
                            float_recheck_all=bool(ok_float)))
    W.check_all_finite(W.jsonable(out))
    json.dump(W.jsonable(out), open(out_path, 'w'), indent=1)
    log('main COMPLETE: max U %.6f, d2(0.03) <= %.6e, limsup <= %.4f, sup (0,1e-4] <= %.4f (%.1f s)'
        % (out['summary']['smallness_max'], out['summary']['d2_0p03'], out['summary']['limsup_d2_over_a3'],
           out['summary']['sup_d2_over_a3_0_1em4'], time.time() - T0))

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
