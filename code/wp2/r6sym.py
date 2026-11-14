# Symbolic proof attempt of the R6 angular reduction (v3 section 7 item 1), H = H6.
import sympy as sp
z1, z2, w1, w2 = sp.symbols('z1 z2 w1 w2')   # w = zbar, independent
P0, P1, P2, F = sp.symbols('P0 P1 P2 F')      # Pi, Pi', Pi'' at sigma; F = fh = sqrt(s^2+1)
z = [z1, z2]; w = [w1, w2]
s = z1*w1 + z2*w2
A1 = z1**5 - 5*z1*z2**4; A2 = z2**5 - 5*z1**4*z2
B1 = w1**5 - 5*w1*w2**4; B2 = w2**5 - 5*w1**4*w2
H = sp.Rational(3, 16)*((w1*A1 + w2*A2) + (z1*B1 + z2*B2))
# sanity: H real form equals 6 m3 - 15/2 p2 m2 + 15/8 p2^3 at a real point
x0, x1, x2, x3 = sp.Rational(3, 7), sp.Rational(-2, 5), sp.Rational(1, 3), sp.Rational(5, 11)
sub = {z1: x0+sp.I*x1, w1: x0-sp.I*x1, z2: x2+sp.I*x3, w2: x2-sp.I*x3}
p2 = x0**2+x1**2+x2**2+x3**2; m2 = x0**4+x1**4+x2**4+x3**4; m3 = x0**6+x1**6+x2**6+x3**6
print('H6 form check:', sp.simplify(sp.expand(H.subs(sub)) - (6*m3 - sp.Rational(15, 2)*p2*m2 + sp.Rational(15, 8)*p2**3)))
dH = [sp.diff(H, v) for v in z]; dbH = [sp.diff(H, v) for v in w]
K = [[sp.diff(H, z[i], w[j]) for j in range(2)] for i in range(2)]
# metric gh0_{i jbar} = (F/s) delta + phi'' w_i z_j, phi'' = (s/F - F/s)/s
php = (s/F - F/s)/s
G = sp.Matrix(2, 2, lambda i, j: (F/s)*(1 if i == j else 0) + php*w[i]*z[j])
M = sp.Matrix(2, 2, lambda i, j: (P1*(1 if i == j else 0) + P2*w[i]*z[j])*H + P1*(w[i]*dbH[j] + z[j]*dH[i]) + P0*K[i][j])
def red(e):
    e = sp.together(sp.expand(e))
    n, d = sp.fraction(e)
    n = sp.expand(n); d = sp.expand(d)
    n = sp.rem(sp.Poly(n, F), sp.Poly(F**2 - s**2 - 1, F)).as_expr()
    d = sp.rem(sp.Poly(d, F), sp.Poly(F**2 - s**2 - 1, F)).as_expr()
    return sp.expand(n), sp.expand(d)
detG = G.det()
n, d = red(detG - 1)
print('det G - 1 numerator zero:', sp.expand(n) == 0)
Ginv = sp.Matrix([[G[1, 1], -G[0, 1]], [-G[1, 0], G[0, 0]]])  # adj; det G = 1
b = (Ginv*M).trace()
detM = M.det()
# claimed invariants at unit x
Hx = H/s**3
rho = z1*dH[0] + z2*dH[1]
rhob = w1*dbH[0] + w2*dbH[1]
Imrho2 = -(rho - rhob)**2/4
U = (dH[0]*dbH[0] + dH[1]*dbH[1])/s**5 - Imrho2/s**6
k2 = sum(K[i][j]*K[j][i] for i in range(2) for j in range(2))/s**4
X = sum(z[i]*K[i][j]*dH[j] for i in range(2) for j in range(2))
Xb = sum(w[i]*K[j][i]*dbH[j] for i in range(2) for j in range(2))
wv = (X + Xb)/2/s**5
a0 = s**2*P0; a1 = s**3*P1; a2 = s**4*P2
beta = (s/F)*(8*a1 + a2) + (7*a1 + a2 + 5*a0)/(s*F)
Gam = beta**2/4 - 16*a1**2 - a1*a2 + 5*a0*a2
n1, d1 = red(b - Hx*beta)
print('b == H*beta :', n1 == 0)
disc = b**2/4 - detM
n2, d2 = red(disc - (Hx**2*Gam + U*a1**2 + k2*a0**2/2 + 2*a1*a0*wv))
print('disc identity:', n2 == 0)
if n2 != 0:
    print('residual terms (count):', len(sp.Add.make_args(n2)))
