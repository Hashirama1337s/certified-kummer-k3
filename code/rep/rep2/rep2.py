"""rep2.py - blind replication seat R2: leading bulk residual F8bulk of the Kummer K3 approximant.

Object (spec DEFINITIONS-AND-LEMMAS-v2.md L4(e), second expression, re-derived in REP2-NOTES.md):
    F8bulk(z) = det( ddbar Psi1 ) + tr( ddbar Psi2 ),
    Psi1 = 8 pi^2 |z|^2 - S/2,   tr ddbar Psi2 = (1/4) [ sum_{l in L} |z-l|^-8 - Z2 ],
with L = Hurwitz/2, V = 1/32, ddbar = complex Hessian (d_i dbar_j), z1 = x0 + i x1, z2 = x2 + i x3.

Formulation used here (different arrangement from L4(e)'s N[psi1]/s^3 + Lap psi2 + det form):
  * every complex Hessian is built from radial one-variable derivatives:
        d_i dbar_j g(|w|^2) = g'(s) delta_ij + g''(s) conj(w_i) w_j ,
    and for plane waves d_i dbar_j cos(2 pi k.x) = -pi^2 cos(2 pi k.x) conj(kc_i) kc_j;
  * the |z-l|^-8 lattice sum is itself Ewald-split (generalised Ewald / Epstein, incomplete gamma),
    not summed directly;
  * rigorous part: arb balls, exact rational lattice data grouped by shells, explicit tails derived in
    REP2-NOTES.md; two Ewald parameters eta = 4 and eta = 9 (the spec's eta = 6 deliberately not used);
  * double part: cancelled arrangement F = -tr(H0 K1) + det K1 + (1/4)(sum_{l!=0}|z-l|^-8 - Z2),
    H0 = ddbar(-1/(2s)), K1 = ddbar(8 pi^2 s - S_reg/2), evaluated with eta = 10 (check eta = 6.5).

Run:  python <lab-tools>/memcap.py 8 py -3.11 rep2.py
"""
import json, math, sys, time, os
from fractions import Fraction
from math import gcd
import numpy as np
from flint import arb, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
ctx.prec = 160
PI = arb.pi()
V = Fraction(1, 32)
T0 = time.time()
LOG = []


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.append(s)


# ----------------------------------------------------------------------------------------------
# A. lattice data
# L = (1/2)(Z^4 u Z^4 + h), h = (1/2,1/2,1/2,1/2).  Dual L* = {k : k.l in Z for all l in L}
#   l = e_i/2  -> k_i even;  l = (1,1,1,1)/4 -> sum k_i in 4Z.  So L* = {2m : m in Z^4, sum m even}.
# ----------------------------------------------------------------------------------------------
BASIS = np.array([[0.5, 0, 0, 0], [0, 0.5, 0, 0], [0, 0, 0.5, 0], [0.25, 0.25, 0.25, 0.25]])
HALF = np.array([0.5, 0.5, 0.5, 0.5])


def L_points_int(Dscale, R):
    """integer coordinates (times Dscale, Dscale divisible by 4) of all l in L with |l| <= R."""
    nmax = int(math.ceil(2 * R)) + 1
    r = np.arange(-nmax, nmax + 1, dtype=np.int64)
    g = np.stack(np.meshgrid(r, r, r, r, indexing="ij"), -1).reshape(-1, 4)
    out = []
    c0 = g * (Dscale // 2)                       # l = n/2
    c1 = (2 * g + 1) * (Dscale // 4)             # l = (n + h)/2 = (2n+1)/4
    for c in (c0, c1):
        s = (c * c).sum(1)
        out.append(c[s <= (R * Dscale) ** 2])
    return np.concatenate(out)


def Ls_points_int(K):
    """integer vectors k of L* with 0 < |k| <= K."""
    mmax = int(math.ceil(K / 2)) + 1
    r = np.arange(-mmax, mmax + 1, dtype=np.int64)
    g = np.stack(np.meshgrid(r, r, r, r, indexing="ij"), -1).reshape(-1, 4)
    g = g[(g.sum(1) % 2) == 0]
    k = 2 * g
    s = (k * k).sum(1)
    return k[(s > 0) & (s <= K * K)]


def reduce_to_V0(x):
    """map points x (N,4) to their representative in the Voronoi cell V0 of L (nearest lattice point)."""
    a = np.round(2 * x) / 2
    b = (np.round(2 * x - 0.5) + 0.5) / 2
    da = ((x - a) ** 2).sum(-1)
    db = ((x - b) ** 2).sum(-1)
    l = np.where((da <= db)[..., None], a, b)
    return x - l


# ----------------------------------------------------------------------------------------------
# B. rigorous arb evaluation at a rational point
# ----------------------------------------------------------------------------------------------
def J3_bound(R, beta):
    """int_R^inf r^3 exp(-beta r^2) dr <= R^2 exp(-beta R^2)/(2 beta - 2/R^2)  (by parts; needs beta R^2 > 1)."""
    R = arb(R); beta = arb(beta)
    assert float((beta * R * R).mid()) > 1
    return R * R * (-beta * R * R).exp() / (2 * beta - 2 / (R * R))


def tails(eta, T, K):
    """upper bounds (arb) for:
       tS  : any complex-Hessian entry of sum_{|x-l|>=T} f_eta(x-l), f_eta = exp(-eta s)/s
       t8r : sum_{|x-l|>=T} Gamma(4, eta s)/(6 s^4)
       tSk : any complex-Hessian entry of the reciprocal S-sum over |k| >= K
       t8k : reciprocal |x|^-8 sum over |k| >= K
    Derivations in REP2-NOTES.md (ball packing: L radius 1/4, L* radius sqrt 2)."""
    eta = arb(eta)
    R = arb(T) - arb(1) / 2
    assert float(R.mid()) >= 1
    c125 = arb(125) / 64                       # (r + 1/4)^3 <= (5/4)^3 r^3 for r >= 1
    phiS = eta * eta + 3 * eta + 3             # e^{eta r^2} (|g'| + s |g''|) <= this for r >= 1
    phi8 = 1 + eta + eta * eta / 2 + eta ** 3 / 6   # e^{eta r^2} Gamma(4,eta r^2)/(6 r^8) <= this for r >= 1
    tS = 1024 * c125 * phiS * J3_bound(R, eta)
    t8r = 1024 * c125 * phi8 * J3_bound(R, eta)
    Rk = arb(K) - 2 * arb(2).sqrt()
    assert float(Rk.mid()) >= 2
    beta = PI * PI / eta
    ck = (1 + arb(2).sqrt() / 2) ** 3         # (r + sqrt2)^3 <= (1 + sqrt2/2)^3 r^3 for r >= 2
    Vb = arb(1) / 32
    tSk = ck * PI * PI / Vb * J3_bound(Rk, beta)  # entry <= pi^2 exp(-pi^2 r^2/eta)/V
    t8k = ck * eta ** 3 / (6 * Vb * Rk * Rk) * J3_bound(Rk, beta)   # <= eta^3/(6 V r^2) exp(-pi^2 r^2/eta)
    return tS, t8r, tSk, t8k


def choose_radii(eta, target):
    T = 2.0
    while True:
        tS, t8r, _, _ = tails(eta, T, 10.0)
        if float(tS.upper()) < target and float(t8r.upper()) < target:
            break
        T += 0.1
    K = 5.0
    while True:
        _, _, tSk, t8k = tails(eta, T, K)
        if float(tSk.upper()) < target and float(t8k.upper()) < target:
            break
        K += 0.1
    return round(T, 2), round(K, 2)


def ball(r):
    return arb(0, r)


def gamma_m2(y):
    return y.gamma_upper(-2)


def Z2_ewald(eta, target=1e-30):
    """Z2 = sum_{l != 0} |l|^-8 by the generalised Ewald split at x = 0 (rigorous)."""
    T, K = choose_radii(eta, target)
    D = 4
    lp = L_points_int(D, T)
    s = (lp * lp).sum(1)
    s = s[s > 0]
    vals, cnt = np.unique(s, return_counts=True)
    etab = arb(eta)
    tot = arb(0)
    for sv, c in zip(vals.tolist(), cnt.tolist()):
        sa = arb(sv) / (D * D)
        y = etab * sa
        tot += c * (-y).exp() * (1 + y + y * y / 2 + y ** 3 / 6) / sa ** 4
    tot += -etab ** 4 / 24                                    # l = 0: lim (Gamma(4,eta s)/6 - 1)/s^4
    tot += PI * PI * etab * etab * 32 / 12                    # k = 0: pi^2 eta^2/(12 V)
    kp = Ls_points_int(K)
    ks = (kp * kp).sum(1)
    kv, kc = np.unique(ks, return_counts=True)
    for kk, c in zip(kv.tolist(), kc.tolist()):
        cc = PI * PI * kk
        tot += c * PI * PI * 32 / 6 * cc * cc * gamma_m2(cc / etab)
    _, t8r, _, t8k = tails(eta, T, K)
    tot += ball(t8r.upper()) + ball(t8k.upper())
    return tot, T, K


def F8bulk_arb(zf, eta, Z2, target=1e-26):
    """rigorous enclosure of F8bulk at rational z = zf (4 Fractions), uncancelled arrangement
       det(8 pi^2 I - H_S/2) + (1/4)(sum_l |z-l|^-8 - Z2), H_S = complex Hessian of S."""
    T, K = choose_radii(eta, target)
    den = 1
    for q in zf:
        den = den * q.denominator // gcd(den, q.denominator)
    D = 4 * den
    zi = np.array([int(q * D) for q in zf], dtype=np.int64)
    zn = math.sqrt(sum(float(q) ** 2 for q in zf))
    lp = L_points_int(D, T + zn + 0.01)
    w = zi[None, :] - lp
    s = (w * w).sum(1)
    keep = s.astype(np.float64) <= (T * D) ** 2 * (1 + 1e-12)
    # keep all |w| <= T (extra points just beyond T only make the truncation smaller; tails are for |w| >= T)
    w = w[keep]; s = s[keep]
    assert (s > 0).all()
    a0, b0, a1, b1 = w[:, 0], w[:, 1], w[:, 2], w[:, 3]
    P00 = a0 * a0 + b0 * b0
    P11 = a1 * a1 + b1 * b1
    P01r = a0 * a1 + b0 * b1
    P01i = a0 * b1 - b0 * a1
    vals, inv, cnt = np.unique(s, return_inverse=True, return_counts=True)
    G = len(vals)
    S00 = np.zeros(G, np.int64); np.add.at(S00, inv, P00)
    S11 = np.zeros(G, np.int64); np.add.at(S11, inv, P11)
    S01r = np.zeros(G, np.int64); np.add.at(S01r, inv, P01r)
    S01i = np.zeros(G, np.int64); np.add.at(S01i, inv, P01i)
    etab = arb(eta)
    D2 = arb(D * D)
    h00 = arb(0); h11 = arb(0); h01r = arb(0); h01i = arb(0); s8 = arb(0)
    for gi in range(G):
        sa = arb(int(vals[gi])) / D2
        e = (-etab * sa).exp()
        gp = -e * (etab / sa + 1 / (sa * sa))
        gpp = e * (etab * etab / sa + 2 * etab / (sa * sa) + 2 / sa ** 3)
        c = int(cnt[gi])
        h00 += c * gp + gpp * int(S00[gi]) / D2
        h11 += c * gp + gpp * int(S11[gi]) / D2
        h01r += gpp * int(S01r[gi]) / D2
        h01i += gpp * int(S01i[gi]) / D2
        y = etab * sa
        s8 += c * e * (1 + y + y * y / 2 + y ** 3 / 6) / sa ** 4
    # reciprocal sums
    kp = Ls_points_int(K)
    kz = (kp * zi[None, :]).sum(1)              # k.z * D  (integer)
    ks = (kp * kp).sum(1)
    for idx in range(len(kp)):
        k = kp[idx]
        kk = int(ks[idx])
        cth = (arb(2 * int(kz[idx])) / D).cos_pi()      # cos(2 pi k.z)
        coefS = 32 * (-PI * PI * kk / etab).exp() / kk   # (1/V) exp(-pi^2|k|^2/eta)/|k|^2
        k0, k1, k2, k3 = (int(t) for t in k)
        q00 = k0 * k0 + k1 * k1
        q11 = k2 * k2 + k3 * k3
        q01r = k0 * k2 + k1 * k3                         # conj(kc_0) kc_1 = (k0 - i k1)(k2 + i k3)
        q01i = k0 * k3 - k1 * k2
        f = -coefS * cth * PI * PI                     # theta = 2 pi k.x: d_i theta = pi conj(kc_i)
        h00 += f * q00; h11 += f * q11; h01r += f * q01r; h01i += f * q01i
        cc = PI * PI * kk
        s8 += PI * PI * 32 / 6 * cc * cc * gamma_m2(cc / etab) * cth
    s8 += PI * PI * etab * etab * 32 / 12
    tS, t8r, tSk, t8k = tails(eta, T, K)
    eS = ball((tS + tSk).upper())
    h00 += eS; h11 += eS; h01r += eS; h01i += eS
    s8 += ball((t8r + t8k).upper())
    trH = h00 + h11
    # K = 8 pi^2 I - H_S/2 ;  det K = K00 K11 - |K01|^2
    K00 = 8 * PI * PI - h00 / 2
    K11 = 8 * PI * PI - h11 / 2
    K01r = -h01r / 2; K01i = -h01i / 2
    detK = K00 * K11 - (K01r * K01r + K01i * K01i)
    tr2 = (s8 - Z2) / 4
    F = detK + tr2
    return dict(F=F, detK=detK, tr2=tr2, trHS=trH, T=T, K=K, nL=len(s), nK=len(kp), shells=G)


# ----------------------------------------------------------------------------------------------
# C. double-precision evaluator (cancelled arrangement), vectorised
# ----------------------------------------------------------------------------------------------
from scipy.special import exp1, gammainc


class Fast:
    def __init__(self, eta, Z2f, Treal=None, Krec=None):
        self.eta = eta
        if Treal is None:
            Treal = math.sqrt(46.0 / eta) + 0.36      # exp(-eta r^2) eta^2 < ~1e-18 beyond, plus cell radius
        if Krec is None:
            Krec = math.sqrt(46.0 * eta) / math.pi
        lp = L_points_int(4, Treal).astype(np.float64) / 4
        lp = lp[(lp * lp).sum(1) > 0]
        self.L = lp
        kp = Ls_points_int(Krec).astype(np.float64)
        self.k = kp
        kk = (kp * kp).sum(1)
        c = math.pi ** 2 * kk
        y = c / eta
        g2 = 0.5 * (np.exp(-y) * (y ** -2 - y ** -1) + exp1(y))         # Gamma(-2, y)
        self.coefS = 32.0 * np.exp(-math.pi ** 2 * kk / eta) / kk
        self.coef8 = math.pi ** 2 * 32 / 6 * c * c * g2
        self.kc0 = kp[:, 0] + 1j * kp[:, 1]
        self.kc1 = kp[:, 2] + 1j * kp[:, 3]
        self.Z2 = Z2f
        self.const8 = math.pi ** 2 * eta ** 2 * 32 / 12

    def v_derivs(self, s):
        eta = self.eta
        y = eta * s
        vp = np.empty_like(s); vpp = np.empty_like(s)
        small = y < 1.0
        if small.any():
            ss = s[small]
            acc1 = np.zeros_like(ss); acc2 = np.zeros_like(ss)
            for n in range(30, 1, -1):
                acc1 = acc1 + (-eta) ** n * (n - 1) * ss ** (n - 2) / math.factorial(n)
            for n in range(30, 2, -1):
                acc2 = acc2 + (-eta) ** n * (n - 1) * (n - 2) * ss ** (n - 3) / math.factorial(n)
            vp[small] = acc1; vpp[small] = acc2
        big = ~small
        if big.any():
            sb = s[big]; e = np.exp(-eta * sb)
            vp[big] = -eta * e / sb - (e - 1) / sb ** 2
            vpp[big] = eta ** 2 * e / sb + 2 * eta * e / sb ** 2 + 2 * (e - 1) / sb ** 3
        return vp, vpp

    def F(self, z, chunk=128, parts=False):
        z = np.atleast_2d(np.asarray(z, dtype=np.float64))
        out = np.empty(len(z))
        unc = np.empty(len(z))
        eta = self.eta
        for i0 in range(0, len(z), chunk):
            zz = z[i0:i0 + chunk]
            w = zz[:, None, :] - self.L[None, :, :]
            s = (w * w).sum(-1)
            e = np.exp(-eta * s)
            gp = -e * (eta / s + 1 / s ** 2)
            gpp = e * (eta ** 2 / s + 2 * eta / s ** 2 + 2 / s ** 3)
            wc0 = w[..., 0] + 1j * w[..., 1]
            wc1 = w[..., 2] + 1j * w[..., 3]
            h00 = (gp + gpp * (wc0.conj() * wc0).real).sum(1)
            h11 = (gp + gpp * (wc1.conj() * wc1).real).sum(1)
            h01 = (gpp * wc0.conj() * wc1).sum(1)
            y = eta * s
            s8 = (e * (1 + y + y * y / 2 + y ** 3 / 6) / s ** 4).sum(1)
            # l = 0 regular parts (S_reg: v(s) = (exp(-eta s) - 1)/s ; |x|^-8 : -P(4, eta s)/s^4)
            s0 = (zz * zz).sum(-1)
            zc0 = zz[:, 0] + 1j * zz[:, 1]
            zc1 = zz[:, 2] + 1j * zz[:, 3]
            vp, vpp = self.v_derivs(s0)
            h00 = h00 + vp + vpp * abs(zc0) ** 2
            h11 = h11 + vp + vpp * abs(zc1) ** 2
            h01 = h01 + vpp * zc0.conj() * zc1
            s8 = s8 - gammainc(4, eta * s0) / s0 ** 4
            # reciprocal
            th = 2 * math.pi * (zz @ self.k.T)
            cth = np.cos(th)
            f = -(self.coefS[None, :] * cth) * math.pi ** 2
            h00 = h00 + (f * (self.kc0.conj() * self.kc0).real[None, :]).sum(1)
            h11 = h11 + (f * (self.kc1.conj() * self.kc1).real[None, :]).sum(1)
            h01 = h01 + (f * (self.kc0.conj() * self.kc1)[None, :]).sum(1)
            s8 = s8 + (self.coef8[None, :] * cth).sum(1) + self.const8
            # K1 = ddbar(8 pi^2 s - S_reg/2)
            K00 = 8 * math.pi ** 2 - h00 / 2
            K11 = 8 * math.pi ** 2 - h11 / 2
            K01 = -h01 / 2
            K10 = K01.conj()
            # H0 = ddbar(-1/(2s)): (H0)_ij = delta_ij/(2 s^2) - conj(z_i) z_j / s^3
            H00 = 1 / (2 * s0 ** 2) - abs(zc0) ** 2 / s0 ** 3
            H11 = 1 / (2 * s0 ** 2) - abs(zc1) ** 2 / s0 ** 3
            H01 = -zc0.conj() * zc1 / s0 ** 3
            H10 = H01.conj()
            trH0K1 = (H00 * K00 + H01 * K10 + H10 * K01 + H11 * K11).real
            detK1 = (K00 * K11 - K01 * K10).real
            out[i0:i0 + chunk] = -trH0K1 + detK1 + (s8 - self.Z2) / 4
            # uncancelled for comparison: det(H0 + K1) + (s8 + 1/s^4 - Z2)/4
            A00 = H00 + K00; A11 = H11 + K11; A01 = H01 + K01; A10 = H10 + K10
            unc[i0:i0 + chunk] = (A00 * A11 - A01 * A10).real + (s8 + 1 / s0 ** 4 - self.Z2) / 4
        if parts:
            return out, unc
        return out


# ----------------------------------------------------------------------------------------------
def H6(x):
    x = np.atleast_2d(x)
    p2 = (x * x).sum(-1)
    m2 = (x ** 4).sum(-1)
    m3 = (x ** 6).sum(-1)
    return (6 * m3 - 7.5 * p2 * m2 + 1.875 * p2 ** 3) / p2 ** 3


def fmt(a, digits=None):
    return a.str(digits or 25, radius=True)


def c6_ewald(eta, T=None, K=None):
    """c6 from the Ewald form of S: P6[S_reg] = c6 H6 (spec L3(e) with psi1 = 8 pi^2|z|^2 - S_reg/2), so
    c6 = (8/3) (1/720) d^6/dx0^6 S_reg(0).  Real space: d^6/dy^6 exp(-t y^2) = t^3 H_6(sqrt t y) exp(-t y^2),
    H_6(u) = 64u^6 - 480u^4 + 720u^2 - 120, integrated over t in (eta, inf) in closed form; l = 0 regular part
    v(s) = (exp(-eta s) - 1)/s contributes 720 eta^4/24 = 30 eta^4; reciprocal: -(1/V) e^{-pi^2 k^2/eta}/k^2 (2 pi k0)^6.
    Double precision (a control, not an enclosure)."""
    if T is None:
        T = math.sqrt(60.0 / eta) + 0.5
    if K is None:
        K = math.sqrt(70.0 * eta) / math.pi
    lp = L_points_int(4, T).astype(np.float64) / 4
    lp = lp[(lp * lp).sum(1) > 0]
    y0 = -lp[:, 0]
    r2 = (lp * lp).sum(1)
    hc = {0: -120.0, 1: 720.0, 2: -480.0, 3: 64.0}
    tot = 0.0
    for m, h in hc.items():
        n = 3 + m
        y = eta * r2
        inc = math.factorial(n) * np.exp(-y) * sum(y ** i / math.factorial(i) for i in range(n + 1)) / r2 ** (n + 1)
        tot += (h * y0 ** (2 * m) * inc).sum()
    tot += 30 * eta ** 4
    kp = Ls_points_int(K).astype(np.float64)
    kk = (kp * kp).sum(1)
    tot += (-32.0 * np.exp(-math.pi ** 2 * kk / eta) / kk * (2 * math.pi * kp[:, 0]) ** 6).sum()
    return 8.0 / 3.0 * tot / 720.0


def main():
    res = {"seat": "R2", "object": "F8bulk of DEFINITIONS-AND-LEMMAS-v2 L4(e)"}
    # ---------------- controls on the lattice sums
    log("== C1: Z2 by generalised Ewald (two eta) vs closed form 24*4^4*(7/8)*zeta(4)*zeta(3) (spec 1.1)")
    Z2a, Ta, Ka = Z2_ewald(4)
    Z2b, Tb, Kb = Z2_ewald(9)
    Z2cf = 24 * arb(256) * arb(7) / 8 * arb(4).zeta() * arb(3).zeta()
    log("  Z2 eta=4 :", fmt(Z2a), "T,K =", Ta, Ka)
    log("  Z2 eta=9 :", fmt(Z2b), "T,K =", Tb, Kb)
    log("  closed   :", fmt(Z2cf))
    ok1 = Z2a.overlaps(Z2b) and Z2a.overlaps(Z2cf)
    log("  overlap:", ok1)
    Z2 = Z2a.intersection(Z2b)
    res["control_Z2"] = {"eta4": fmt(Z2a), "eta9": fmt(Z2b), "closed_form": fmt(Z2cf), "overlap": bool(ok1)}

    # ---------------- enclosure (1)
    z0 = [Fraction(1, 4), Fraction(0), Fraction(0), Fraction(0)]
    log("== (1) rigorous enclosure at z = (1/4,0,0,0)")
    encs = {}
    for eta in (4, 9):
        t = time.time()
        r = F8bulk_arb(z0, eta, Z2)
        encs[eta] = r
        log("  eta=%d: F8bulk = %s" % (eta, fmt(r["F"])))
        log("     detK = %s  tr2 = %s  trace(H_S) = %s (must be pi^2/V = 32 pi^2 = %s)" % (
            fmt(r["detK"], 20), fmt(r["tr2"], 20), fmt(r["trHS"], 20), (32 * PI * PI).str(20)))
        log("     T=%.2f K=%.2f  nL=%d nK=%d shells=%d  %.1fs" % (r["T"], r["K"], r["nL"], r["nK"], r["shells"], time.time() - t))
    ov = encs[4]["F"].overlaps(encs[9]["F"])
    Fint = encs[4]["F"].intersection(encs[9]["F"])
    trc = encs[4]["trHS"].overlaps(32 * PI * PI) and encs[9]["trHS"].overlaps(32 * PI * PI)
    log("  overlap eta4/eta9:", ov, " intersection:", fmt(Fint), " trace control:", trc)
    res["enclosure_z_quarter"] = {
        "z": "(1/4, 0, 0, 0)",
        "eta4": fmt(encs[4]["F"]), "eta9": fmt(encs[9]["F"]),
        "overlap": bool(ov), "intersection": fmt(Fint),
        "intersection_lo": Fint.lower().str(25, radius=False), "intersection_hi": Fint.upper().str(25, radius=False),
        "radii_eta4_T_K": [encs[4]["T"], encs[4]["K"]], "radii_eta9_T_K": [encs[9]["T"], encs[9]["K"]],
        "control_trace_ddbarS_equals_32pi2": bool(trc),
    }
    # ---------------- double evaluator vs arb
    Z2f = float(Z2.mid())
    fa = Fast(10.0, Z2f)
    fb = Fast(6.5, Z2f)
    log("== C2: double evaluator (eta=10, eta=6.5; cancelled and uncancelled) vs arb at (1/4,0,0,0)")
    zq = np.array([[0.25, 0, 0, 0]])
    va, ua = fa.F(zq, parts=True); vb = fb.F(zq)
    log("  eta10 canc %.15g  unc %.15g  eta6.5 %.15g  arb mid %s" % (va[0], ua[0], vb[0], Fint.mid().str(17)))
    res["control_double_vs_arb_at_quarter"] = {"eta10": repr(float(va[0])), "eta6.5": repr(float(vb[0])),
                                               "rel_diff_eta10_vs_arb": "%.2e" % abs(va[0] / float(Fint.mid()) - 1)}
    # a second rational check point (generic, near a bubble)
    zg = [Fraction(3, 32), Fraction(-1, 16), Fraction(5, 64), Fraction(1, 128)]
    rg = F8bulk_arb(zg, 9, Z2, target=1e-20)
    vg = fa.F(np.array([[float(q) for q in zg]]))[0]
    log("  z = (3/32,-1/16,5/64,1/128): arb %s  double %.15g  rel %.2e" % (fmt(rg["F"], 20), vg, abs(vg / float(rg["F"].mid()) - 1)))
    res["control_double_vs_arb_generic"] = {"z": "(3/32, -1/16, 5/64, 1/128)", "arb_eta9": fmt(rg["F"], 20),
                                            "double_eta10": repr(float(vg))}

    # ---------------- C3: lattice-point limit and the cell average (Remark V)
    log("== C3: near-lattice-point limit: F8bulk(t*u) for t -> 0 on the axis and on (1,1,0,0)/sqrt2")
    lim = {}
    for name, u in (("axis", np.array([1.0, 0, 0, 0])), ("diag", np.array([1, 1, 0, 0]) / math.sqrt(2))):
        ts = np.array([0.08, 0.04, 0.02, 0.01, 0.005])
        vals = fa.F(ts[:, None] * u[None, :])
        lim[name] = vals.tolist()
        log("  %s:" % name, " ".join("%.3f" % v for v in vals))
    res["control_limit"] = {"t": [0.08, 0.04, 0.02, 0.01, 0.005], **{k: [repr(x) for x in v] for k, v in lim.items()}}

    from scipy.stats import qmc
    log("== C4: cell average of F8bulk (QMC over the fundamental parallelepiped, 2^15 and 2^16 points)")
    avgs = []
    for m, seed in ((15, 11), (16, 12)):
        u = qmc.Sobol(d=4, scramble=True, seed=seed).random_base2(m)
        x = reduce_to_V0(u @ BASIS)
        vals = fa.F(x)
        avgs.append(vals.mean())
        log("  2^%d points: mean = %.2f  (std err naive %.2f)" % (m, vals.mean(), vals.std() / math.sqrt(len(vals))))
    pred = float((2 * PI ** 4 - Z2 / 128).mid()) * 32
    log("  Remark V prediction (2 pi^4 - Z2/128)/V = %.3f" % pred)
    res["control_cell_average"] = {"qmc_2^15": "%.2f" % avgs[0], "qmc_2^16": "%.2f" % avgs[1],
                                   "remark_V_prediction": "%.3f" % pred}

    log("== C5: c6 (sextic coefficient, spec L3(e)) from the Ewald form, two eta; sets the lattice-point limit -(5 c6/2) H6")
    c6a = c6_ewald(5.0); c6b = c6_ewald(11.0)
    c6v = 0.5 * (c6a + c6b)
    log("  c6 eta=5: %.10f  eta=11: %.10f  -> (15/16) c6 = %.6f" % (c6a, c6b, 15 / 16 * c6v))
    res["control_c6"] = {"eta5": "%.10f" % c6a, "eta11": "%.10f" % c6b, "limit_15c6_over_16": "%.6f" % (15 / 16 * c6v)}

    # ---------------- (2) sampled sup over V0
    log("== (2) sampled sup over V0: 2^15 scrambled Sobol points (seed 2026) reduced to V0, eta=10")
    u = qmc.Sobol(d=4, scramble=True, seed=2026).random_base2(15)
    X = reduce_to_V0(u @ BASIS)
    t = time.time()
    FX = fa.F(X)
    log("  evaluated %d points in %.1fs: max %.4f min %.4f sup|F| %.4f" % (len(X), time.time() - t, FX.max(), FX.min(), np.abs(FX).max()))
    res["sample_raw"] = {"n": int(len(X)), "max": "%.4f" % FX.max(), "min": "%.4f" % FX.min(),
                         "argmax_raw": X[FX.argmax()].tolist(), "argmin_raw": X[FX.argmin()].tolist()}
    from scipy.optimize import minimize

    def proj(q, rf):
        n = float(np.linalg.norm(q))
        if rf <= 0 or n >= rf:
            return q
        if n == 0:
            return np.array([rf, 0, 0, 0.0])
        return q * (rf / n)

    def refine(sign, rf):
        """Nelder-Mead from the best 20 sample points, then restarts with shrinking simplices until two
        successive maxima agree to 1e-3 relative.  rf > 0: points with |z| < rf are projected radially to
        |z| = rf (keeps the double evaluator away from its cancellation zone at the lattice point)."""
        order = np.argsort(-sign * FX)[:20]
        f1 = lambda q: -sign * fa.F(proj(q, rf)[None, :])[0]
        best = None
        for i in order:
            r = minimize(f1, X[i], method="Nelder-Mead",
                         options=dict(xatol=1e-9, fatol=1e-8, maxiter=1500, initial_simplex=X[i] + 0.01 * np.vstack([np.zeros(4), np.eye(4)])))
            if best is None or r.fun < best.fun:
                best = r
        hist = [-best.fun]
        p = best.x
        for it in range(6):
            r = minimize(f1, p, method="Nelder-Mead",
                         options=dict(xatol=1e-12, fatol=1e-10, maxiter=1500, initial_simplex=p + 10.0 ** (-3 - it) * np.vstack([np.zeros(4), np.eye(4)])))
            p = r.x
            hist.append(-r.fun)
            if abs(hist[-1] - hist[-2]) <= 1e-3 * abs(hist[-1]) and it >= 1:
                break
        p = reduce_to_V0(proj(p, rf)[None, :])[0]
        return p, sign * hist[-1], [sign * h for h in hist]

    pmin, vmin, hmin = refine(-1, 0.0)
    log("  refined MIN F = %.6f at %s |z| = %.6f (history %s)" % (vmin, np.array2string(pmin, precision=10), np.linalg.norm(pmin), ["%.6f" % h for h in hmin]))
    a1, u1 = fa.F(pmin[None, :], parts=True); b1 = fb.F(pmin[None, :])
    log("  at argmin: eta10 %.10f  eta6.5 %.10f  uncancelled %.10f" % (a1[0], b1[0], u1[0]))
    dy = [Fraction(int(round(c * 2 ** 16)), 2 ** 16) for c in pmin]
    emin = F8bulk_arb(dy, 9, Z2, target=1e-20)["F"]
    log("  arb at the 2^-16 dyadic rounding of argmin %s: %s" % ([str(q) for q in dy], fmt(emin, 20)))
    maxs = {}
    for rf in (0.01, 0.005):
        pmax, vmax, hmax = refine(+1, rf)
        maxs[rf] = (pmax, vmax, hmax)
        log("  refined MAX F with floor |z| >= %g: %.6f at %s |z| = %.6f H6(zhat) = %.6f (history %s)" % (
            rf, vmax, np.array2string(pmax, precision=8), np.linalg.norm(pmax), H6(pmax)[0], ["%.6f" % h for h in hmax]))
    # behaviour at the lattice point: the 24 directions (+-e_i +- e_j)/sqrt2 (H6 = -3/8) and Richardson in t^2
    dirs = []
    for i in range(4):
        for j in range(i + 1, 4):
            for si in (1, -1):
                for sj in (1, -1):
                    d = np.zeros(4); d[i] = si; d[j] = sj; dirs.append(d / math.sqrt(2))
    dirs = np.array(dirs)
    ts = [0.02, 0.01, 0.005]
    tab = np.array([fa.F(t * dirs) for t in ts])          # (3, 24)
    A = np.array([[1, t * t, t ** 4] for t in ts])
    lim24 = np.linalg.solve(A, tab)[0]
    log("  24 H6=-3/8 directions: F at t=0.02/0.01/0.005: max %.4f / %.4f / %.4f ; min %.4f / %.4f / %.4f" % (
        tab[0].max(), tab[1].max(), tab[2].max(), tab[0].min(), tab[1].min(), tab[2].min()))
    log("  Richardson (F0 + c t^2 + d t^4) limits over the 24 directions: min %.4f max %.4f ; (15/16) c6 = %.4f" % (
        lim24.min(), lim24.max(), 15 / 16 * c6v))
    zl = [Fraction(1, 512), Fraction(1, 512), Fraction(0), Fraction(0)]
    el = F8bulk_arb(zl, 9, Z2, target=1e-20)["F"]
    log("  arb at z = (2^-9, 2^-9, 0, 0): %s  (double: %.6f)" % (fmt(el, 20), fa.F(np.array([[2 ** -9, 2 ** -9, 0, 0]]))[0]))
    # special points: the 24 midpoints l/2 (l a minimal vector, |l| = 1/2) and the 24 deep holes (+-1/4, +-1/4, 0, 0)-type
    mins = []
    for i in range(4):
        for s in (0.5, -0.5):
            v = np.zeros(4); v[i] = s; mins.append(v)
    for sg in np.array(np.meshgrid(*[[0.25, -0.25]] * 4)).T.reshape(-1, 4):
        mins.append(sg)
    mins = np.array(mins) / 2
    holes = []
    for i in range(4):
        for j in range(i + 1, 4):
            for si in (0.25, -0.25):
                for sj in (0.25, -0.25):
                    v = np.zeros(4); v[i] = si; v[j] = sj; holes.append(v)
    holes = np.array(holes)
    fm = fa.F(mins); fh = fa.F(holes)
    log("  24 midpoints l/2: F in [%.6f, %.6f]; 24 deep holes: F in [%.6f, %.6f]" % (fm.min(), fm.max(), fh.min(), fh.max()))
    res["special_points"] = {"midpoints_l_over_2_min_max": ["%.6f" % fm.min(), "%.6f" % fm.max()],
                             "deep_holes_min_max": ["%.6f" % fh.min(), "%.6f" % fh.max()]}
    p1, v1, h1 = maxs[0.005]
    supv = max(abs(v1), abs(vmin), 15 / 16 * c6v)
    res["sampled_sup"] = {
        "n_points": 32768, "sobol_seed": 2026, "eta": 10.0,
        "min_F": "%.5f" % vmin, "argmin": [repr(c) for c in pmin.tolist()], "min_history": ["%.6f" % h for h in hmin],
        "arb_at_dyadic_argmin_2^-16": fmt(emin, 20), "argmin_dyadic": [str(q) for q in dy],
        "max_F_note": "sampled evidence, not a proof: no interior maximum found; F increases towards the lattice point along the H6 = -3/8 directions; sup F is the limit (15/16) c6, not attained on V0 minus {0}",
        "max_F_sup_limit_15c6_over_16": "%.4f" % (15 / 16 * c6v),
        "max_F_richardson_limit_24dirs": ["%.4f" % lim24.min(), "%.4f" % lim24.max()],
        "max_F_with_floor_0.01": "%.4f" % maxs[0.01][1], "argmax_floor_0.01": [repr(c) for c in maxs[0.01][0].tolist()],
        "max_F_with_floor_0.005": "%.4f" % maxs[0.005][1], "argmax_floor_0.005": [repr(c) for c in maxs[0.005][0].tolist()],
        "arb_at_(2^-9,2^-9,0,0)": fmt(el, 20),
        "sup_abs_F": "%.5f" % supv, "sup_abs_attained_at": "argmin (interior critical point, |z| = 1/4, midpoint between two neighbouring bubbles)" if abs(vmin) >= supv else "see max",
        "H6_dir_argmin": repr(H6(pmin)[0]),
    }
    # ---------------- (3) Q2
    q2 = q2_points()
    if q2 is not None:
        P, desc = q2
        d = np.sqrt((P * P).sum(1))
        sel = P[d >= 0.1]
        vals = fa.F(sel)
        valsb = fb.F(sel)
        smax = np.abs(vals).max()
        log("== (3) Q2: %s; %d points kept at distance >= 0.1; sup|F8bulk| = %.4f (eta6.5: %.4f), max %.4f, min %.4f" % (
            desc, len(sel), smax, np.abs(valsb).max(), vals.max(), vals.min()))
        res["Q2"] = {"point_set": desc, "n_total": int(len(P)), "n_kept": int(len(sel)),
                     "sup_abs_F8bulk": "%.4f" % smax, "sup_abs_eta6.5": "%.4f" % np.abs(valsb).max(),
                     "max": "%.4f" % vals.max(), "min": "%.4f" % vals.min(),
                     "rel_diff_vs_spec_14569.1": "%.2e" % abs(smax / 14569.1 - 1)}
    res["seconds"] = round(time.time() - T0, 1)
    with open(os.path.join(HERE, "rep2_result.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    with open(os.path.join(HERE, "rep2_run.log"), "w") as fh:
        fh.write("\n".join(LOG) + "\n")
    log("done %.1fs" % (time.time() - T0))


def q2_points():
    """The K1d 'cell' group (spec L4(e) check (1)): read off probe1/k1d.py main_sample and k1c.to_nearest
    (opened AFTER this file's evaluator was written; only the point generator was taken):
    rng = numpy default_rng(7); cell = rng.random((400, 4)) @ B (B = the same basis as BASIS, drawn first),
    z = x - nearest lattice point; keep |z| >= 0.1."""
    rng = np.random.default_rng(7)
    cell = rng.random((400, 4)) @ BASIS
    return reduce_to_V0(cell), "K1d cell group: default_rng(7).random((400,4)) @ B, reduced to nearest lattice point"


if __name__ == "__main__":
    main()
