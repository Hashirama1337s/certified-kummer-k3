"""BAR-S2b runner: imports the sealed s2.py unchanged and replaces ONLY its double-precision comparator
theorem2_float by a cancellation-free version (log1p/expm1, same formulas, same explicit x rule), then runs s2.main.
Usage: memcap 8 py -3.11 s2b.py --main result_s2b.json   (refuses without GO-S2b)"""
import sys, os, math
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s2

def theorem2_float_cf(a, Cq, CR, c0s, D, gfun, Lam0=2.0, iters=6, Kfac=2.0):
    eta = -2 * math.log1p(-Cq * a ** 8)
    em = -math.expm1(-eta)                      # 1 - e^-eta, no cancellation
    eta1 = 2 * (-math.expm1(-eta / 2)); eta2 = 2 * math.expm1(eta / 2)
    kappa = CR * a ** 6; sigma = 2 * kappa; C0 = c0s / a ** 2
    delta = 4 * em + sigma / C0 * math.exp(-eta)
    den = 1 - math.expm1(eta) / 2
    def step(Lam, K2):
        G1 = gfun(kappa, D); G2 = gfun(K2, math.sqrt(Lam) * D)
        o = (2 * eta1 * G1 + 2 * eta2 * G2) / den
        xm1 = math.sqrt(delta / (4 * C0 * o)); x = 1 + xm1
        # log((1 + x + sqrt((x-1)^2 + x delta)) / (2x)) = log1p((sqrt(xm1^2 + x delta) - xm1) / (2x))
        q = math.sqrt(xm1 * xm1 + x * delta)
        L = eta + math.log1p((q - xm1) / (2 * x)) + C0 * (1 + x) * o
        eps2 = math.expm1(L)
        m1 = eps2 + math.sqrt(eps2 * (2 + eps2) + em)
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

s2.theorem2_float = theorem2_float_cf
s2.GO_FILE = os.path.join(HERE, 'GO-S2b')
s2.STATE_DIR = os.path.join(HERE, 'state_s2b')     # fresh state: no reuse of the S2 run
_orig_logger = s2.logger
s2.logger = lambda _path: _orig_logger(os.path.join(HERE, 'run_s2b.log'))   # separate log file
if __name__ == '__main__':
    if len(sys.argv) != 3 or sys.argv[1] != '--main': sys.exit('usage: s2b.py --main OUT.json')
    if not os.path.exists(os.path.join(HERE, 'GO-S2b')): sys.exit('REFUSED: GO-S2b missing')
    print('BAR-S2b: cancellation-free comparator installed', flush=True)
    s2.main(sys.argv[2], 520.0) if hasattr(s2, 'main') else sys.exit('s2.main not found')
