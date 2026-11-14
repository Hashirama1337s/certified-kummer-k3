"""rep1.py - blind replication seat R1 (K3 Kummer window, Theorem 2 / Route 2).

Rigorous (arb ball arithmetic, python-flint) evaluation of the symmetric closeness
d2 = M2 e^eta - 1 of ANALYSIS-CHAIN.md Theorem 2, with
  * Lemma G evaluated two ways: (G) c0/Doeblin with a rigorous left Riemann sum on an
    exact dyadic time grid; (expl) the closed form Gamma_expl;
  * x chosen two ways: (rule) x = 1 + sqrt(delta/(4 C0 o)); (opt) float-optimised x,
    then evaluated rigorously at that exact x (any x > 0 is valid).
Floats are used ONLY to choose free parameters (s per time node, T, x); every reported
bound is the upper endpoint of an arb ball computed at those exact parameter values.

Run:  python rep1.py   -> prints, writes rep1_result.json next to this file.
"""
import json, os, sys, time
import numpy as np
from flint import arb, ctx

HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------- helpers
def A(v):
    """exact arb from str/int/float (floats are exact binary numbers)."""
    return arb(v)

def up(b):
    return b.upper()

def fl(b):
    return float(b.mid())

def s_up(b, n=10):
    """decimal string that is >= the upper endpoint (rounded up in the last digit)."""
    u = b.upper()
    m = float(u.mid())
    from decimal import Decimal, ROUND_CEILING, getcontext
    getcontext().prec = 40
    d = Decimal(u.mid().str(40, radius=False))  # mid of upper ball
    # pad by the ball radius, then round up to n significant digits
    d = d + Decimal(float(u.rad())) * 2
    q = Decimal(1).scaleb(d.adjusted() - n + 1)
    return str(d.quantize(q, rounding=ROUND_CEILING))

SQRT2_F = np.sqrt(2.0)

# ---------------------------------------------------------------- Lemma G, float side (parameter choice only)
def _f_obj(s, tau, k):
    """(sqrt P + sqrt Q)^2, P = 1/(4 tau (1-s)) - 2 log s, Q = 2 sqrt2 k tau (1-s); D scaled to 1."""
    P = 1.0 / (4.0 * tau * (1.0 - s)) - 2.0 * np.log(s)
    Q = 2.0 * SQRT2_F * k * tau * (1.0 - s)
    return (np.sqrt(P) + np.sqrt(Q)) ** 2

def best_s(taus, k, ngrid=800, iters=50):
    """vectorised: grid search over s in (0,1) then golden refinement of the bracket."""
    taus = np.asarray(taus, float)
    sg = np.linspace(0.0, 1.0, ngrid + 2)[1:-1]
    out = np.empty_like(taus)
    for c0 in range(0, len(taus), 512):
        tt = taus[c0:c0 + 512][:, None]
        F = _f_obj(sg[None, :], tt, k)
        j = np.argmin(F, axis=1)
        lo = sg[np.maximum(j - 1, 0)] * 0 + np.where(j > 0, sg[np.maximum(j - 1, 0)], 1e-9)
        hi = np.where(j < ngrid - 1, sg[np.minimum(j + 1, ngrid - 1)], 1 - 1e-9)
        t = tt[:, 0]
        g = (np.sqrt(5) - 1) / 2
        for _ in range(iters):
            x1 = hi - g * (hi - lo); x2 = lo + g * (hi - lo)
            f1 = _f_obj(x1, t, k); f2 = _f_obj(x2, t, k)
            m = f1 < f2
            hi = np.where(m, x2, hi); lo = np.where(m, lo, x1)
        out[c0:c0 + 512] = 0.5 * (lo + hi)
    return out

# ---------------------------------------------------------------- Lemma G, rigorous side
def c0_low(tau, s, k):
    """rigorous lower bound (arb, take .lower()) of exp(-(sqrt P + sqrt Q)^2) at exact tau, s; k an arb ball."""
    one = A(1)
    P = one / (4 * tau * (one - s)) - 2 * s.log()
    Q = 2 * A(2).sqrt() * k * tau * (one - s)
    return (-(P.sqrt() + Q.sqrt()) ** 2).exp()

_gh_cache = {}

def gammahat_G(k_ball, hexp=11, Tmax=3.0, refine=True):
    """Rigorous upper bound of Gammahat(k) = inf_T [int_0^T (1 - m) + T(1 - c0(T))/c0(T)], D = 1 scaling.
    Grid tau_i = i h, h = 2^-hexp exact.  Left Riemann sum is valid because m(t) = inf V H is
    nondecreasing (semigroup), so on [tau_i, tau_{i+1}] 1 - m <= 1 - m(tau_i) <= 1 - c0(tau_i).
    k_ball: arb enclosing K D^2 (any K >= true Ricci bound valid).  Returns (bound_ball, T, n_nodes)."""
    key = (k_ball.mid().str(30, radius=False), k_ball.rad().str(5, radius=False), hexp, Tmax, ctx.prec)
    if key in _gh_cache:
        return _gh_cache[key]
    h = 2.0 ** (-hexp)
    n = int(Tmax / h)
    kf = max(fl(k_ball) + float(k_ball.rad()), 0.0)
    taus = h * np.arange(1, n + 1)
    S = best_s(taus, kf)
    c0f = np.exp(-_f_obj(S, taus, kf))
    # float objective for each T_j = j h: h * sum_{i<j} (1 - c0_i) [c0_0 = 0] + T_j (1 - c0_j)/c0_j
    oneminus = np.concatenate([[1.0], 1.0 - c0f])          # index i = 0..n
    cums = h * np.cumsum(oneminus)[:-1]                    # sum_{i<j}, j = 1..n
    obj = cums + taus * (1.0 - c0f) / c0f
    j = int(np.argmin(obj)) + 1
    if j >= n - 2:
        raise RuntimeError("T optimum at grid edge; raise Tmax")
    # rigorous sum on the exact dyadic grid
    hb = A(h)
    tot = hb * 1  # i = 0 node: c0(0) := 0 -> contributes h * 1
    for i in range(1, j):
        tau = hb * i
        tot += hb * (1 - c0_low(tau, A(float(S[i - 1])), k_ball).lower())
    T = hb * j
    cT = c0_low(T, A(float(S[j - 1])), k_ball).lower()
    tot += T * (1 - cT) / cT
    res = (tot, float(T.mid()), j, float(cT.mid()))
    _gh_cache[key] = res
    return res

def gammahat_expl(k_ball):
    """Gamma_expl / D^2 = (3/4) exp((sqrt(1 + 2 ln(3/2)) + sqrt(k/sqrt2))^2)  (T = 3D^2/4, s = 2/3)."""
    p = (1 + 2 * (A(3) / 2).log()).sqrt()
    q = (k_ball / A(2).sqrt()).sqrt() if fl(k_ball) > 0 or k_ball.rad() > 0 else A(0)
    return A(3) / 4 * ((p + q) ** 2).exp()

def Gamma(K, Dsq, mode, hexp=11):
    """Gamma(K, D) = D^2 Gammahat(K D^2) (exact scaling t = D^2 tau)."""
    k = K * Dsq
    if mode == "G":
        return Dsq * gammahat_G(k, hexp=hexp)[0]
    return Dsq * gammahat_expl(k)

# ---------------------------------------------------------------- inputs
CERT = dict(
    Cq="17226.2012683894",
    k6="1731732.44",
    a2C0="2.0000352761032563",
    e="1.368159761580671e-4",
)

def D_L5(a, e):
    """L5 exact: D(a,e) = sqrt(2(1+e)) [4a + (17/16)^(1/4) (1/sqrt2 + 12(pi-2) a)]."""
    return (2 * (1 + e)).sqrt() * (4 * a + (A(17) / 16).root(4) * (1 / A(2).sqrt() + 12 * (arb.pi() - 2) * a))

def inputs_cert(a):
    e = A(CERT["e"])
    return dict(Cq=A(CERT["Cq"]), kappa=A(CERT["k6"]) * a ** 6, C0=A(CERT["a2C0"]) / a ** 2,
                D=D_L5(a, e))

def inputs_brief(a):
    """ANALYSIS-CHAIN sec. 2 brief inputs (control Q1)."""
    return dict(Cq=A("2.7e4"), kappa=A("2.7e6") * a ** 6, C0=A(2) / a ** 2,
                D=(1 + 140 * a ** 4).sqrt() * (A("1.0153") + A("25.33") * a))

# ---------------------------------------------------------------- Theorem 2 chain
def make_chain(a, inp, gmode, xmode, eta_both="spec", hexp=11):
    Cq, kappa, C0, D = inp["Cq"], inp["kappa"], inp["C0"], inp["D"]
    eta = -2 * (1 - Cq * a ** 8).log()
    sigma = 2 * kappa
    Dsq = D ** 2
    ee = eta.exp()
    hyp = dict(Cq_a8_lt_1=bool(Cq * a ** 8 < 1), e_eta_lt_3=bool(ee < 3), fourC0_gt_sigma=bool(4 * C0 > sigma))
    assert all(hyp.values()), hyp
    eta1 = 2 * (1 - (-eta / 2).exp())
    eta2 = 2 * ((eta / 2).exp() - 1)
    if eta_both == "eta2":
        eta1 = eta2
    G1 = Gamma(kappa, Dsq, gmode, hexp)          # Lemma G on g_a: K = kappa, diameter D
    delta = 4 * (1 - (-eta).exp()) + sigma / C0 * (-eta).exp()

    def osc(Lam, Kt):
        G2 = Gamma(Kt, Lam * Dsq, gmode, hexp)    # Lemma G on g_t: K_t, diameter sqrt(Lam) D
        return (2 * eta1 * G1 + 2 * eta2 * G2) / (1 - (ee - 1) / 2)

    def onepe2(o, x):
        return ee * (1 + x + ((x - 1) ** 2 + x * delta).sqrt()) / (2 * x) * (C0 * (1 + x) * o).exp()

    def M2(o):
        if xmode == "rule":
            x = 1 + (delta / (4 * C0 * o)).sqrt()
        else:
            x = A(opt_x(fl(C0), fl(o), fl(delta), fl(eta), a))
        q = onepe2(o, x)
        return q + (q ** 2 - (-eta).exp()).sqrt(), x

    return dict(eta=eta, sigma=sigma, D=D, delta=delta, G1=G1, osc=osc, M2=M2, ee=ee, hyp=hyp)

def opt_x(C0, o, delta, eta, a):
    """float choice of x minimising log(1+eps2) (mpmath, high precision); any x > 0 is admissible."""
    import mpmath as mp
    mp.mp.dps = 60
    C0, o, delta = mp.mpf(C0), mp.mpf(o), mp.mpf(delta)
    f = lambda u: mp.log((2 + mp.e ** u + mp.sqrt(mp.e ** (2 * u) + (1 + mp.e ** u) * delta)) / (2 * (1 + mp.e ** u))) \
        + C0 * (2 + mp.e ** u) * o
    u0 = mp.log(mp.sqrt(delta / (4 * C0 * o)))
    lo, hi = u0 - 6, u0 + 6
    g = (mp.sqrt(5) - 1) / 2
    for _ in range(200):
        x1 = hi - g * (hi - lo); x2 = lo + g * (hi - lo)
        if f(x1) < f(x2):
            hi = x2
        else:
            lo = x1
    return float(1 + mp.e ** ((lo + hi) / 2))

def theorem2(a, inp, gmode, xmode, eta_both="spec", hexp=11, iters=400):
    ch = make_chain(a, inp, gmode, xmode, eta_both, hexp)
    # smallness: M2(Lambda=2, kappa' = 2 kappa) e^eta < 2  (K_t = (1-t) 2 kappa, worst t = 0)
    o_s = ch["osc"](A(2), 2 * inp["kappa"])
    M2s, xs = ch["M2"](o_s)
    small = M2s * ch["ee"]
    small_ok = bool(small < 2)
    # t = 1 iteration: Lambda_{j+1} = max(U, 1/L) = M2(Lambda_j, K=0) e^eta, Lambda_0 = 2
    Lam = A(2)
    hist = []
    for j in range(iters):
        o1 = ch["osc"](Lam, A(0))
        M2v, x1 = ch["M2"](o1)
        new = up(M2v * ch["ee"])
        hist.append(fl(new))
        if j > 3 and fl(Lam - new) <= 1e-25 * fl(new - 1):
            Lam = new
            break
        Lam = new
    d2 = Lam - 1   # = M2(Lambda*, 0) e^eta - 1 with Lambda* the previous iterate
    return dict(d2=d2, small=small, small_ok=small_ok, Lambda_star=Lam, o1=o1, x1=x1, xs=xs,
                eta=ch["eta"], D=ch["D"], G1=ch["G1"], o_small=o_s, iters=len(hist), hyp=ch["hyp"])

# ---------------------------------------------------------------- limit a -> 0
def limit_coeff(gmode, hexp=11, e_val=None):
    """lim d2/a^3 = 2 sqrt( lim a^-6 C0 o ),  o ~ 2*(2 Cq a^8)*(Gamma(0,D0) + Gamma(0,D0)) [Lambda* -> 1],
    C0 a^2 -> c2 (certified constant bound), D0 = sqrt(1+e) (17/16)^(1/4).
    => lim = 2 sqrt(8 c2 Cq D0^2 Gammahat(0)) = 2 D0 sqrt(8 c2 Cq Gammahat(0))."""
    e = A(CERT["e"]) if e_val is None else A(e_val)
    D0sq = (1 + e) * (A(17) / 16).sqrt()
    gh0 = gammahat_G(A(0), hexp=hexp)[0] if gmode == "G" else gammahat_expl(A(0))
    return 2 * (8 * A(CERT["a2C0"]) * A(CERT["Cq"]) * D0sq * gh0).sqrt(), gh0

# ---------------------------------------------------------------- main
def main():
    t0 = time.time()
    res = {"seat": "R1", "spec": "ANALYSIS-CHAIN.md Thm 2 (sec 9), Prop C0 (sec 4), Lemma G (sec 3), closure (sec 7); DEFINITIONS-v2 L5",
           "arith": "python-flint arb ball arithmetic; strings are upper bounds rounded up"}
    out = {}
    for prec in (256, 512):
        ctx.prec = prec
        _gh_cache.clear()
        R = {}
        # Lemma G constants
        for hexp in (10, 11):
            g0, T0, n0, c0T = gammahat_G(A(0), hexp=hexp)
            R[f"gammahat0_G_h2^-{hexp}"] = (g0, T0, n0, c0T)
        R["gammahat0_expl"] = gammahat_expl(A(0))
        # (1) limits
        for gm in ("G", "expl"):
            L, _ = limit_coeff(gm)
            R[f"lim_{gm}"] = L
            Le0, _ = limit_coeff(gm, e_val=0)
            R[f"lim_{gm}_e0"] = Le0
        R["lim_G_h2^-10"] = limit_coeff("G", hexp=10)[0]
        # (2),(3)
        for aa in ("0.03", "0.01"):
            a = A(aa)
            for gm in ("G", "expl"):
                for xm in ("rule", "opt"):
                    R[f"T2_{aa}_{gm}_{xm}"] = theorem2(a, inputs_cert(a), gm, xm)
        # Gamma at a=0.03 with the coarser grid (convergence in h)
        a = A("0.03")
        R["T2_0.03_G_rule_h2^-10"] = theorem2(a, inputs_cert(a), "G", "rule", hexp=10)
        # convergence of d2/a^3 to the limit (small a), expl + G
        for aa in ("1e-3", "1e-4", "1e-5"):
            a = A(aa)
            for gm in ("G", "expl"):
                r = theorem2(a, inputs_cert(a), gm, "rule")
                R[f"ratio_{aa}_{gm}"] = r["d2"] / a ** 3
        # Q1 control
        for aa in ("1e-3", "1e-2"):
            a = A(aa)
            R[f"Q1_{aa}_spec"] = theorem2(a, inputs_brief(a), "expl", "rule")
            R[f"Q1_{aa}_eta2both"] = theorem2(a, inputs_brief(a), "expl", "rule", eta_both="eta2")
        out[prec] = R
        print(f"prec {prec} done, {time.time() - t0:.1f}s", flush=True)

    R = out[512]; Rl = out[256]
    def pr(k, b):
        print(f"{k:34s} {b.str(16)}")
    # precision agreement check
    def agree(k, f=lambda v: v):
        x, y = f(R[k]), f(Rl[k])
        return abs(fl(x) - fl(y)) <= 1e-12 * abs(fl(x)) + 1e-300

    J = {}
    J["controls"] = {}
    g10, g11 = R["gammahat0_G_h2^-10"], R["gammahat0_G_h2^-11"]
    J["controls"]["gammahat0_G_upper_h2^-11"] = s_up(g11[0])
    J["controls"]["gammahat0_G_upper_h2^-10"] = s_up(g10[0])
    J["controls"]["gammahat0_G_T_over_D2"] = f"{g11[1]:.6f}"
    J["controls"]["gammahat0_G_c0T"] = f"{g11[3]:.6f}"
    J["controls"]["gammahat0_expl_27e_over_16"] = s_up(R["gammahat0_expl"])
    for aa in ("1e-3", "1e-2"):
        J["controls"][f"Q1_d2_a{aa}_spec_eta1eta2"] = s_up(R[f"Q1_{aa}_spec"]["d2"], 6)
        J["controls"][f"Q1_d2_a{aa}_eta2both"] = s_up(R[f"Q1_{aa}_eta2both"]["d2"], 6)
    J["controls"]["Q1_table_targets_R2_DL5_expl"] = {"1e-3": "2.931e-6", "1e-2": "3.594e-3"}

    J["limit_d2_over_a3_upper"] = {
        "LemmaG": s_up(R["lim_G"], 8), "LemmaG_coarse_h2^-10": s_up(R["lim_G_h2^-10"], 8),
        "expl": s_up(R["lim_expl"], 8),
        "LemmaG_if_e_to_0": s_up(R["lim_G_e0"], 8), "expl_if_e_to_0": s_up(R["lim_expl_e0"], 8)}
    J["d2_over_a3_small_a_rule_x"] = {k[6:]: s_up(R[k], 8) for k in R if k.startswith("ratio_")}
    J["d2_upper"] = {}
    J["smallness_a0.03"] = {}
    J["details"] = {}
    for aa in ("0.03", "0.01"):
        for gm in ("G", "expl"):
            for xm in ("rule", "opt"):
                r = R[f"T2_{aa}_{gm}_{xm}"]
                J["d2_upper"][f"a{aa}_{gm}_{xm}"] = s_up(r["d2"], 8)
                J["details"][f"a{aa}_{gm}_{xm}"] = {
                    "eta": s_up(r["eta"], 8), "D_L5": s_up(r["D"], 10), "Gamma(kappa,D)": s_up(r["G1"], 8),
                    "o_at_t1": s_up(r["o1"], 8), "x_at_t1": f"{fl(r['x1']):.10g}",
                    "Lambda_star_iterations": r["iters"],
                    "smallness_M2_e^eta": s_up(r["small"], 8), "hypotheses": r["hyp"]}
                if aa == "0.03":
                    J["smallness_a0.03"][f"{gm}_{xm}"] = {"M2(2,2kappa)e^eta_upper": s_up(r["small"], 8),
                                                          "lt_2": r["small_ok"]}
    J["d2_upper"]["a0.03_G_rule_coarse_h2^-10"] = s_up(R["T2_0.03_G_rule_h2^-10"]["d2"], 8)
    # precision agreement
    keys = [k for k in R if k.startswith(("T2_", "Q1_"))]
    J["precision_256_vs_512_agree_1e-12"] = all(agree(k, lambda v: v["d2"]) for k in keys) and \
        all(agree(k) for k in R if k.startswith(("lim_", "ratio_")))
    J["runtime_s"] = round(time.time() - t0, 1)
    res.update(J)
    print(json.dumps(res, indent=1))
    with open(os.path.join(HERE, "rep1_result.json"), "w") as fh:
        json.dump(res, fh, indent=1)

if __name__ == "__main__":
    main()
