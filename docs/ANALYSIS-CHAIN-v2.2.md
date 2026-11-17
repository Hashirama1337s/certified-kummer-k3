# K3 Kummer window: the analysis chain, v2.2 (editor seat, 2026-11-14; v2.1 + the text fixes of REFEREE-LEMMAS and REFEREE-A7, see CHANGES-v4.2.md; the v2 header follows)

v1 (ANALYSIS-CHAIN.md, WP4 proof seat) was refereed twice: REFEREE-WP4 (seats/REFEREE-WP4.md) and REFEREE-CHAIN (seats/REFEREE-CHAIN.md). Neither found an error in the mathematics of Route 2; both asked for text fixes. v2 makes every REFEREE-CHAIN fix (which include the open REFEREE-WP4 fixes) and replaces the brief inputs by the certified ones (WP2a, WP2b, WP3, WP5). **The v2 text itself is NOT re-refereed.** One line per change, with source: CHANGES-v4.md, part B. v1 is unchanged.

**Status words:** PROVED (argument complete here); OPENED (statement read in the source this session, locator given); STANDARD (textbook fact, source not opened); INPUT (named hypothesis; in v2 each input is CERTIFIED, i.e. enclosed by a sealed, scored run whose seal is given; exception: WP2a, whose c6, P* and R6* feed C_q and C_R, was scored FAILED-as-sealed on its rounded containment targets W1 and W2, see the rests-on list in section 0); GAP (not closed); CHECK (a computation compared, never a proof).

Scripts of v1 in <lab>/wp4/ (green_const.py, torus_ctrl.py, referee_gamma.py, chain.py, explicit.py, each with a .log); Li-Yau page images wp4/liyau/liyau1986_p14..17.png. The certified evaluation of Theorem 2 is <lab>/wp5/wp5.py (sealed 891fafa2d9), output wp5/result_wp5.json, verdict wp5/VERDICT-WP5.md (post-erratum).

Seals of the inputs: WP2a claim 2d19816000; WP2b claim cdf666b973; WP3 claim 4158aba4f8 (+ AMEND-WP3-1 eddd0d7e7b); WP5 claim 055009634d (+ AMEND-WP5-1 b4c1dabedd, instrument 891fafa2d9).

## 0. Result

**Theorem (certified instance of Theorem 2; computer-assisted; conditional only as listed at the end of this section).** Let omega_a be the approximant of DEFINITIONS-AND-LEMMAS-v4 1.5 (Kummer surface of the Hurwitz lattice, 16 Eguchi-Hanson bubbles of size a, C-infinity cut-off on the window [a/4, a]) and omega_CY the Ricci-flat Kahler form in its class. For every a in (0, 0.03] the continuity method of section 7 closes (smallness U = M2 e^eta <= 1.1694 < 2 on all 30 a-intervals of width 0.001), and

  (1 + d2(a))^-1 omega_a <= omega_CY <= (1 - d2(a))^-1 omega_a,

with the following certified upper bounds (VERDICT-WP5 post-erratum and result_wp5.json; every printed number is the certified value rounded UP):
- d2(a)/a^3 <= 2275.00 for every a in (0, 1e-4] (JSON 2274.997908).
- d2(a)/a^3 <= 2326.6 on (0, 0.001]; <= 2851.0 on [0.009, 0.01]; <= 4517.6 on [0.029, 0.03].
- d2(0.01) <= 2.851e-3; d2(0.03) <= 0.121974.
- limsup_{a -> 0} d2(a)/a^3 <= 2269.27. This is the a = 0 value of a continuous upper bound (JSON 2269.269445); it does not bound d2(a)/a^3 at any a > 0 (REFEREE-WP5 item 2: the exact bound is 2269.2701 at a = 1e-6). Its closed form is 8 D(0) sqrt(C_q Gammahat(0)) with D(0) <= 1.015342 (JSON 1.015341043), C_q <= 17226.21 (JSON 17226.20127) and Gammahat(0) <= 4.530841 (JSON 4.530840314).
- a0 = 0.03 (criterion in section 9). The binding limit is the range of the diameter lemma L5, not the smallness condition, which is far from failing at a = 0.03 (1.1694 < 2). a0 = 0.03 is where the certificate stops (the L5 range); it is not a property of the method.

Certified d2 per a-interval (each bound holds for every a in the interval; it is essentially the value at the upper end; BUILD-WP5 section 7 item 5). Rounded up from result_wp5.json `stageB`:

| a-interval | smallness U (continuity step) | d2 <= | d2/a^3 <= |
|---|---|---|---|
| [0, 0.001] | 1.000003 | 2.327e-6 | 2326.6 |
| [0.002, 0.003] | 1.000081 | 6.592e-5 | 2441.3 |
| [0.009, 0.01] | 1.003501 | 2.851e-3 | 2851.0 |
| [0.029, 0.03] | 1.1694 | 0.121974 | 4517.6 |

(a = 0.01 also lies in [0.01, 0.011], whose enclosure gives the weaker but valid 3.876e-3; VERDICT-WP5 V4.)

Leading order: d2 ~ 2 sqrt(C0 o) ~ 8 D sqrt(C_q gamma) a^3 with C0 = 2/a^2, o ~ 8 C_q gamma D^2 a^8, gamma = Gammahat(0).

**Route 1 (Theorem 1, Green + Chern-Lu) is not part of the certificate.** Its v1 numbers used the brief inputs and D = 1.1 and are withdrawn (historical note below).

**What the certified theorem rests on (not yet closed):**
- The v2 text of this chain was refereed in full by REFEREE-v4 (seats/REFEREE-v4.md: mathematics CORRECT; wording fixes applied in v2.1); the v2.2 changes (pointwise C0, REFEREE-LEMMAS and REFEREE-A7 fixes) were made by an editor pass from those reports and checked by the lead (consumer audit of C0, e, ||dR||, LOG 00:17); they have not had a separate referee pass.
- BUILD-WP5 lemmas S1-S12 and the scaled rewrite of Theorem 2 have had one referee pass (REFEREE-WP5: CORRECT / CORRECT WITH FIXES; fixes in DEFINITIONS-AND-LEMMAS-v4, Appendix).
- STANDARD facts, not opened: heat-kernel facts (positivity, symmetry, semigroup, mass 1), Kahler normal coordinates, Yau existence and Calabi uniqueness, Schauder/Fredholm theory for the implicit function step, |W(F4)| and its fundamental chamber (used for c0s through WP2b's W-reduction).
- Floating-point layer: IEEE round-to-nearest in parts of WP2b, WP3 and the WP5 sphere covers (S10).
- (a) The v4 write-ups of R6, F2, F3, the W-reduction, the S11 proof and the Z1 radius were first refereed in REFEREE-v4 (seats/REFEREE-v4.md), which found them correct with the fixes made in this version.
- (b) Lemma E, the binomial lemma and the multi-direction Cauchy estimate are PROVED; refereed: REFEREE-LEMMAS (seats/REFEREE-LEMMAS.md, CORRECT / CORRECT WITH FIXES, wording only, fixes made in DEFINITIONS v4.2); WP2b's e and c0s use them.
- (b2) Lemma A7 (chart transfer of jets) and zone Z1 are PROVED; refereed: REFEREE-A7 (seats/REFEREE-A7.md: CORRECT WITH FIXES; the C0/B+ definition is now the pointwise sup, A7 index, Wirtinger sentence and remark (2) fixed in DEFINITIONS v4.2); Z1 sets a^2 C0 = 2.0000353.
- (c) WP2a, whose c6, P* and R6* feed C_q and C_R, was scored FAILED-as-sealed on its rounded containment targets W1 and W2; the header's "each input is enclosed by a sealed, scored run" is to be read with that exception.
- (d) V5 (VERDICT-WP5) is POST HOC.
- Replication status (supersedes REFEREE-v4's "no blind second implementation yet"): blind replication now EXISTS for the assembly arithmetic and the bulk peak (rep/VERDICT-REP.md, verdict as sealed: REPLICATED), and the clean room reproduced the chain (cleanroom/CLEANROOM-REPORT.md, verdict: CLEAN ROOM REPRODUCED). The patch-interior (Stage A) blind replication is also REPLICATED (rep/VERDICT-REP3.md: blind core Ricci max 1,695,996.5 a^6 vs the certificate's own sample 1,695,994; independent rigorous Z2 bound 1.746877e6; results predate the instrument pin, disclosed).

### Historical note (v1 numbers, superseded; do not quote as results)
v1 headlined a0 = 0.0470 (Route 2) and 0.0283 (Route 1) and tabulated d at a = 1e-3, 3e-3, 1e-2, 3e-2 with the brief inputs C_q = 2.7e4, kappa = 2.7e6 a^6, B+ = C0 = 2/a^2 and D = 1.1 or D_L5 (with a placeholder e). These are withdrawn:
- D = 1.1 is a valid diameter bound only for a <= 0.0033 (section 10). The D = 1.1 rows at a = 1e-2 and 3e-2 used an invalid input, and a0 = 0.047 rests on them (REFEREE-WP4 item 6; REFEREE-CHAIN item 8 (f): even extrapolating D_L5 beyond its proved range, smallness fails between a = 0.04 (1.680) and 0.045 (2.736), so 0.047 is false on its own terms).
- The brief inputs are superseded by the certified ones. For the record, v1's "R2, D_L5, opt" row read 2.913e-6, 8.253e-5, 3.572e-3, 0.1581 at the four a, and WP5 control X1 reproduced the "R2, D_L5, expl" row (2.931e-6, 8.303e-5, 3.594e-3) from the same inputs.

## 1. Conventions (reconciled; unchanged from v1)
- n = 2 complex, m = 4 real.
- omega = i g_{jkbar} dz^j ^ dzbar^k, g_R = 2 Re g, so |X|^2 = g(X, Xbar) and g_R(X+Xbar, X+Xbar) = 2|X|^2. Volume dV_R = omega^2/2.
- Delta_c u = g^{jkbar} d_j d_kbar u = tr_omega(i d dbar u). Delta_R = div grad is Li-Yau's Laplacian (heat equation (Delta - d_t)u = 0, p.168). Delta_R = 2 Delta_c.
- R_{i jbar k lbar} = -d_k d_lbar g_{i jbar} + g^{p qbar} d_k g_{i qbar} d_lbar g_{p jbar}; Bisec(X,Y) = R(X,Xbar,Y,Ybar)/(|X|^2|Y|^2) (DEFINITIONS 1.6). Fubini-Study max Bisec 2; Eguchi-Hanson max 2/a^2 and min -2/a^2 (minimum on an orthogonal pair; DEFINITIONS-v4 Lemma EH-min).
- Ric_{i jbar} = g^{k lbar} R_{i jbar k lbar} = -d_i d_jbar log det g. Ric_R(X_R, X_R) = 2 Ric(X, Xbar), so "Ric >= -K" means the same with the same K in both conventions; Li-Yau's K is our kappa.
- F = log(omega_a^2 / omega_flat^2), omega_flat from dz1 ^ dz2 (invariant under -1, no zeros on X). F is global, Ric(omega_a) = -i d dbar F, omega_CY^2 = e^(-F+c) omega_a^2.
- Notation clash: DEFINITIONS uses eta for the Ewald parameter and B for f/s; here they mean the inputs below.
- REFEREE-CHAIN item 8 (a) checked that D, C0 and kappa use the same convention (g_R = 2 Re g); a mixed convention would have shown up as a factor sqrt 2 in d2, which is invariant under rescaling omega.

## 2. Inputs and path
All inputs below are CERTIFIED for every a in (0, 0.03] (seals in the header). Per-interval values are in result_wp5.json `stageB/<interval>/inputs` and `input_sources`.
- **(I1)** |F| <= C_q a^8 on X, with C_q a^8 < 1. Put eta_F := -log(1 - C_q a^8) and eta := 2 eta_F. Then eta_F >= C_q a^8 >= sup|F|, so osc F <= 2 sup|F| <= eta (valid, slightly conservative; REFEREE-CHAIN item 7). Certified: C_q = max(C_core_F, C_out) = 17226.21 (rounded up from 17226.20127; WP3 bulk C_out, set at the Voronoi facet centre), with the full core value C_core_F <= 14416.38 (WP5 Stage A). Both constants bound |F|/a^8 directly (via |log(1 + y)| <= |y|/(1 - |y|); BUILD-WP5 3.3, BUILD-WP3 F5). (v1 wrote I1 as |q| <= C_q a^8; REFEREE-CHAIN item 7.)
- **(I2)** |Ric(omega_a)| <= kappa := C_R a^6, sigma := sup|Delta_c^a F| <= 2 kappa (trace of a 2x2 form with eigenvalues in [-kappa, kappa]). Certified: C_R = max(C_core_R, R_out) <= 1,731,732.5 (WP5 Stage A core, JSON 1731732.443; the core sets the constant, R_out <= 65064.3 from WP3). The certificate uses the per-interval value (1695936.40 on [0, 0.001]).
- **(I3+)** Bisec(omega_a) <= B+ for all pairs, X = Y included. Used by Route 1 only. Certified anyway: a^2 B+ <= 2.0000353 (WP2b).
- **(I3-)** R^a(X, Xbar, Y, Ybar) >= -C0 |X|^2 |Y|^2 for g_a-orthogonal X, Y (Route 2). Source: DEFINITIONS-AND-LEMMAS-v4 L4 (g), Lemma M (the mirror lemma: -(K + D_p)/(1 - e_p)^2 <= Bisec_g <= (K + D_p)/(1 - e_p)^2) with K <= 2/a^2 by L1 and Lemma EH-min (the Eguchi-Hanson minimum -2/a^2 is attained by a g0-orthonormal pair on E). So C0 := the pointwise supremum C0 := sup_p (K(s_p) + D_p)/(1 - e_p)^2 over the points p of X, with K(s_p) <= 2/a^2 the Eguchi-Hanson curvature bound at the point and D_p, e_p the local ||dR|| and e bounds at p, enclosed zone by zone (on [0.025, 0.03], a^2 times the zone values: Z1 2.0000352761032563, Z2 1.9884300760979716, NB 0.7078465004722245, all upper bounds, from result_wp2b.json; the maximum is the Z1 value 2.0000352761), the (1 - e_p)^-2 factor kept (REFEREE-WP4 item 5; pointwise sup per REFEREE-A7 item 5). The global formula (2/a^2 + ||dR||)/(1 - e)^2 with the global e and ||dR|| is kept only as a stated cruder fallback (REFEREE-A7 item 5): it gives a^2 C0 = 2.0035992 on [0.025, 0.03], and any step reading e and ||dR|| separately must use the per-zone pairs. Lemma YA (DEFINITIONS 5.2) uses the curvature bound pointwise at the maximum point of log T - Au, so the pointwise sup is the right input; with the fallback 2.0036 the final d would grow by about 0.09% (d ~ 2 sqrt(C0 o) scales as sqrt(C0): sqrt(2.0035992/2.0000353) - 1 = 0.089%). Certified: **a^2 C0 <= 2.0000353 for every a in (0, 0.03]** (WP2b; JSON 2.0000352761 on [0.025, 0.03]; per-row values in DEFINITIONS-v4 L4 (g)). Lemma M, Lemma D, Lemma EH-min and Lemma M2 were refereed: REFEREE-v3 item 1, CORRECT. Lower bound at a = 0.03 on E: a^2 C0_true >= 1.99999862 (Lemma M2), so the input is tight. (v1: "B+ = C0 = 2/a^2 + small".)
- **(I4)** diam(g_R(omega_a)) <= D(a, e) from DEFINITIONS L5 (valid for a <= 0.03). The certificate uses L5's exact first expression with e as in (I5): D(0) <= 1.015342, D(0.03) <= 1.77518 (BUILD-WP5 section 8 item 4).
- **(I5)** omega_a > 0: e <= 1.37e-4 < 1 (the certificate's conservative max of WP2b 8.34e-5, WP3 1.37e-4 and the Stage A core value; DEFINITIONS-v4 4.4 reconciles the first two).
- **(I6)** Regularity: section 8 (route R1, everything C-infinity).

Path: h_t = -tF + c_t, c_t fixed by int e^{h_t} omega_a^2 = int omega_a^2. Then (omega_a + i d dbar u_t)^2 = e^{h_t} omega_a^2, int u_t omega_a^2 = 0.

**Lemma 2.1 (PROVED; REFEREE-CHAIN item 1: CORRECT).** (a) c_t in [t inf F, t sup F], so |h_t| <= t osc F <= eta. (b) c_0 = 0, u_0 = 0. (c) Ric(omega_t) = (1-t) Ric(omega_a); Ric(omega_1) = 0, so omega_1 = omega_CY (Calabi uniqueness). (d) Delta_c^a h_t >= -sigma. (e) If omega_a <= Lambda omega_t then Ric(omega_t) >= -(1-t) Lambda kappa omega_t. (f) If omega_t <= Lambda omega_a then diam_t <= sqrt(Lambda) D.

## 3. Lemma G (explicit oscillation bound)
**3.1 Opened.** Li-Yau, Acta Math 156 (1986), Thm 2.3, pp.168-169 (wp4/liyau/liyau1986_p15-16.png). Hypotheses: M compact, Ric >= -K, K >= 0, convex (here empty) boundary, u > 0 solving (Delta - d_t)u = 0. Conclusion: for alpha > 1, 0 < t1 < t2, u(x,t1) <= u(y,t2)(t2/t1)^{n alpha/2} exp((n alpha/sqrt2)(alpha-1)^-1 K(t2-t1) + alpha r^2/(4(t2-t1))). Thm 2.2(ii) (C6 not explicit) is not used; Cor 2.1 is the K = 0 case.
**n is the real dimension, n = dim_R M = 4 (v2: GAP 6 CLOSED).** OPENED, p.170 (wp4/liyau/liyau1986_p17.png), proof of Cor 2.1: the estimate follows "since H(x,x,t1) ~ (4 pi t1)^{-n/2} as t1 -> 0". This small-time on-diagonal asymptotic holds only for n = dim_R M (REFEREE-CHAIN item 2; REFEREE-WP4 item 1b reads the same from Cor 2.1 on p.169, whose lower bound (4 pi t)^{-n/2} exp(-r^2/4t) is the Euclidean kernel).

**3.2 Statement.** M^4 closed, Ric >= -K, diam <= D, Delta_R v >= -c. Then sup v - avg v <= c Gamma(K,D), where Gamma(K,D) := inf_T [int_0^T (1 - c0(t)) dt + T(1 - c0(T))/c0(T)] <= inf_T T/c0(T), c0(t) := sup_{0<s<1} exp(-(sqrt P + sqrt Q)^2), P = D^2/(4t(1-s)) - 2 log s, Q = 2 sqrt2 K t(1-s).
Closed form (T = 3D^2/4, s = 2/3): Gamma <= Gamma_expl := (3D^2/4) exp((sqrt(1 + 2 ln(3/2)) + sqrt(K D^2/sqrt2))^2); Gamma_expl(0, D) = (27e/16) D^2 = 4.5871 D^2.

**3.3 Proof** (heat kernel facts STANDARD; REFEREE-CHAIN item 2 and REFEREE-WP4 items 1a-1c: CORRECT).
(i) H(x,x,t) = int H(x,z,t/2)^2 dz >= 1/V (Cauchy-Schwarz).
(ii) Harnack with u = H(x,.,.), t1 = st, t2 = t: 1/V <= H(x,x,st) <= H(x,y,t) s^(-2 alpha) exp(...). With alpha = 1 + delta the exponent is -(1+delta)P - Q(1+delta)/delta; its sup over delta > 0 is -(sqrt P + sqrt Q)^2 at delta = sqrt(Q/P); K = 0 by the limit alpha -> 1. Hence V H >= c0(t).
(iii) m(t) := inf V H is nondecreasing (semigroup). Writing V H(x,z,kT) = m_k + r, r >= 0: m_{k+1} >= m_k + c0(1 - m_k), so 1 - m(t) <= (1 - c0)^k for t >= kT.
(iv) v(x) - int H(x,y,tau) v = -int_0^tau int (H - 1/V) Delta v. With Delta v = -c + w, w >= 0, int w = cV: -int (H - 1/V) Delta v <= c(1 - m(t)). Let tau -> inf (||H - 1/V||_L1 <= 2(1 - m) -> 0): v(x) - avg v <= c int_0^inf (1 - m).
(v) Sum the time intervals. QED. V cancels: the chain does not need Remark V.

**3.4 Numbers.** CHECK (double precision, v1): Gamma(K,D) = D^2 Gammahat(K D^2). Gammahat(0) = 4.5308 at T = 0.788 D^2, c0 = 0.1717. Gammahat(k): 1e-4: 4.638; 1e-3: 4.878; 1e-2: 5.735; 3e-2: 6.841; 0.1: 9.745; 0.5: 27.07 (reproduced by REFEREE-CHAIN with its own s-optimiser: 4.530828, 4.87803, 5.73453, 9.74512). These double-precision values are CHECK values, not certified bounds (the referee's own optimiser gives 4.530828, below the printed 4.5308); they would be a valid upper bound in principle because c0 is a max over a finite s-grid (any s is admissible), the time integral is a left Riemann sum (valid as m is nondecreasing), and Gammahat is nondecreasing in k.
CERTIFIED (WP5, BUILD-WP5 3.4, arb, N = 4000, s-grid 2001): **Gammahat(0) <= 4.530841** (JSON Gh0_upper = 4.530840314), +4e-5 above the double value, below the closed form 4.5871 (WP5 control X2). The certificate uses Gamma(K, D) := min(D^2 Gammahat_cert(upper(K D^2)), Gamma_expl(K, D)); REFEREE-WP5 item 4: CORRECT.

**3.5 Earlier referee's derivation.** G >= -(T/V)(1 + M/c0), M = 4 e^{D^2/4T} (on-diagonal bound from integrating Harnack over y): correct but lossy; reproduced 25.0739 at T = 1.069 D^2. M is unnecessary once Doeblin is applied to inf H.

**3.6 Positive control.** Flat torus R^4/Z^4, V = 1, D = 1: min H(0,.,t) >= c0(t) at every tested t; sharp -min G = 0.0702 <= 4.531 (about 65x slack, in the Harnack step; true int(1 - m) ~ 0.08). An earlier truncated image sum gave a false 20.36; fixed.

## 4. C^0 step
**Proposition C0 (PROVED; REFEREE-CHAIN item 3: CORRECT).** Assume Lambda^-1 omega_a <= omega_t <= Lambda omega_a. Set eta1 = 2(1 - e^{-eta/2}), eta2 = 2(e^{eta/2} - 1), K_t = (1-t) Lambda kappa, and assume e^eta < 3. Then osc u_t <= o := [2 eta1 Gamma(kappa, D) + 2 eta2 Gamma(K_t, sqrt(Lambda) D)] / (1 - (e^eta - 1)/2). The same holds with any K >= K_t in place of K_t, because Gamma is nondecreasing in K (used in section 7).
Proof. S1: Delta_c^a u = T - 2 >= 2e^{-eta/2} - 2 = -eta1; Delta_c^t u = 2 - S <= 2 - 2e^{-eta/2} = eta1 <= eta2 (AM-GM with det = e^h). S2: Riemannian factor 2; Lemma G on g_a (K = kappa, diameter D) for u and on g_t (K_t, diameter sqrt(Lambda) D) for -u, each with its own measure. Averages differ by |avg_t u - avg_a u| <= (osc/2)(e^eta - 1), since int (e^{h_t} - 1) dV_a = 0. The certificate uses eta2 (conservative). Leading order o ~ 2 eta gamma D^2 (1 + Lambda).

## 5. C^2 step
**5.1 Lemma CL (PROVED, constant 1; Route 1 only, not used by the certificate).** If Ric(g) >= -kappa' g and Bisec(ghat) <= B+, then Delta_g log tr_g ghat >= -kappa' - B+ tr_g ghat. g-normal coordinates at p with ghat diagonal (lambda_i): Delta S = sum Ric_{i ibar} lambda_i - sum Rhat_{i ibar k kbar} + sum |d_k ghat_{i pbar}|^2/lambda_p; Cauchy-Schwarz gives |dS|^2/S <= last sum; Rhat_{i ibar k kbar} <= B+ lambda_i lambda_k. Control: ghat = g gives Delta S = 0. Rubinstein Prop 7.2 has 2C3, not reconciled (GAP 4); using 2B raises d1 by 2^(1/4). REFEREE-CHAIN item 8 (e): the constant 1 matches the standard Schwarz-lemma form (STANDARD, not opened).
Maximum principle on log S - Au, A = B+(1+x): S <= 2(1+eps1), 1+eps1 = (1 + 1/x + kappa'/(2B+ x)) e^{B+(1+x)o}; optimum x ~ (B+ o)^(-1/2).

**5.2 Lemma YA (PROVED; constant verified by REFEREE-WP4 item 2 and REFEREE-CHAIN item 4).** Delta_t log T >= (Delta_a h + n^2 C0)/T - C0 S. g_a-normal coordinates with g_t diagonal (mu_i): curvature terms = (1/2) sum R^a_{i ibar k kbar}(mu_i - mu_k)^2/(mu_i mu_k) >= -C0(TS - n^2) (the i = k terms vanish, so only g_a-orthogonal pairs enter, and only the lower bound (I3-) is needed); gradient terms by Cauchy-Schwarz and Kahler symmetry d_k g_{i ibar} = d_i g_{k ibar}. Matches Kobayashi p.300 eq. (31) and Lye p.28 Prop 6.3 eq. (6.10) once R_{1 1bar 2 2bar} is replaced by its worst case -C0 (TS >= 4); Lye (6.10) OPENED in the lab's extraction wp4/lye_s6.txt (REFEREE-CHAIN item 4). Lye Lemma 6.5 and Kobayashi use A = 2R_a.
Maximum principle on log T - Au: in dimension 2, S = T e^{-h}, so at the maximum (A - C0) e^{-h} T^2 - 2AT + (4C0 - sigma) <= 0. With A = C0(1+x), delta = 4(1 - e^{-eta}) + (sigma/C0) e^{-eta}: T(p)/2 <= e^eta (1 + x + sqrt((x-1)^2 + x delta))/(2x), and T <= 2(1+eps2) with 1+eps2 = that x e^{C0(1+x)o}. The upper root decreases in E = e^{-h} for either sign of 4 - sigma/C0 (REFEREE-CHAIN item 5), so E = e^{-eta} is the worst case. When eta = sigma = 0 the root is exactly 2 for every x >= 1 (Kobayashi eq. (33)).
Choose x = 1 + sqrt(delta/(4 C0 o)): **eps2 = 2 C0 o + eta + O(sqrt(C0 o delta))**. The square root enters at the next order: sqrt(C0 o delta) = O(a^7) against 2 C0 o = O(a^6) (REFEREE-WP4 item 2; REFEREE-CHAIN item 4). (v1 said "no square root"; that is true only at leading order.)

**5.3 Why Route 2 is used.** Chern-Lu at the maximum has no +1/S term, so S(p) > 2 strictly, and balancing against e^{Ao} costs a square root at leading order. Yau's +n^2 C0/T term with det = e^h pins T(p) at 2.

## 6. Trace to eigenvalues (PROVED; REFEREE-CHAIN item 5: CORRECT)
Route 1: eigenvalues of omega_a relative to omega_t have sum S and product e^{-h}. With M1 = 1 + eps1 + sqrt((1+eps1)^2 - e^{-eta}): M1^-1 omega_a <= omega_t <= M1 e^eta omega_a; with d1 = M1 - 1 the symmetric form follows because e^eta (1 - d1^2) <= 1.
Route 2: eigenvalues of omega_t relative to omega_a have sum T and product e^h. mu_max <= M2 := (1+eps2) + sqrt((1+eps2)^2 - e^{-eta}), mu_min >= e^{-eta}/M2, so e^{-eta} M2^-1 omega_a <= omega_t <= M2 omega_a; with d2 = M2 e^eta - 1: M2(1 - d2) = 2M2 - M2^2 e^eta <= 1 - (M2-1)^2 <= 1 and the lower bound is exactly 1/(1 + d2), which gives the symmetric form.

## 7. Continuity closure (PROVED; REFEREE-CHAIN item 6: CORRECT WITH FIXES, the K_t fix made here)
1. Existence: Yau 1978 gives u_t in C-infinity for each t (STANDARD; Kobayashi pp.297, 302 the opened application). Calabi: uniqueness.
2. Continuity in t (IFT): Psi(u, c, t) = log((omega_a + i d dbar u)^2/omega_a^2) - h_t - c, from U in C^{4,beta}_0 x R x [0,1] to C^{2,beta}; linearisation D Psi(v, c') = Delta_{omega_u} v - c', an isomorphism (injective by integrating against omega_u^2, surjective by Schauder/Fredholm, STANDARD). c(t) = 0 by cohomology; the IFT branch is u_t by uniqueness. Hence t -> omega_t continuous in C^0.
3. Bootstrap set P = {t: Lambda0^-1 omega_a <= omega_t <= Lambda0 omega_a}, Lambda0 = 2. 0 in P; P closed under C^0 limits. On P: Ric(omega_t) >= -(1-t) Lambda0 kappa omega_t, diam <= sqrt(Lambda0) D, |h| <= eta, Delta_a h >= -sigma. The inputs are uniform in t; the worst case is t = 0, since Gammahat is nondecreasing in K. Smallness condition: U := max(upper, 1/lower) = M2 e^eta < Lambda0, evaluated with Lambda = Lambda0, diameter sqrt(Lambda0) D, and any K_t >= Lambda0 kappa = 2 kappa. Under it P is open (strict improvement), so P = [0,1].
   **K_t used by the certificate.** wp5.py evaluates the continuity step with K_t = 4 kappa (`step(Lam0, Lam0)`: Kfac x Lambda0 x kappa with Kfac = Lambda0 = 2; wp5.py l.797-800, l.820), not the 2 kappa stated in v1 and in BUILD-WP5 3.2. This is conservative and valid: a larger K only enlarges Gamma(K, sqrt(Lambda0) D), hence o, eps2, M2 and U (Proposition C0 holds with any K >= K_t). Effect at a = 0.03 (REFEREE-CHAIN item 6): faithful 2 kappa gives U = 1.161299, the certificate's 4 kappa gives 1.169325 (reported 1.1694).
4. Final bound at t = 1: Ric(omega_CY) = 0, so K_t = 0. Iterate Lambda_{j+1} = max(U, 1/L) at t = 1, each step using only the bound proven by the previous one (the certificate does six iterations per a-interval; result_wp5.json `Lambda_seq`).

## 8. Regularity (v2: rewritten; REFEREE-CHAIN item 6)
**Adopted: route R1.** The cut-off is the C-infinity step st(t) = h(t)/(h(t) + h(1-t)), h(t) = e^{-1/t}, on the window [a/4, a] (DEFINITIONS-v4 1.5; lead's decision, LOG 18:05). Then Phi_a and omega_a are C-infinity on X (DEFINITIONS L2 (iv)), and every step of the chain is classical:
- Yau existence and Calabi uniqueness apply to smooth data, so u_t and omega_t are C-infinity for every t.
- The IFT of section 7 needs omega_a in C^{2,beta} and F in C^{2,beta}: satisfied.
- Lemma G (Li-Yau uses the Bochner formula) needs smooth metrics g_a and g_t: satisfied.
- Lemmas CL and YA need Phi_a in C^4, u in C^4 and F in C^2: satisfied.
The constants that v1 marked "must be recomputed" for R1 (INPUT GAP) have been recomputed for this window and certified: C_q and the core and outer Ricci constants by WP2a (leading core profiles), WP3 (|z| >= a) and WP5 Stage A (full core); C0, B+ and e by WP2b (section 2).
(R2, the mollification of a C^{4,1} quintic-step approximant, is no longer needed and is not used. It is kept here only as a remark: it would remove a regularity gap for a non-smooth cut-off, modulo the STANDARD Hoelder d dbar-lemma and mollification.)

## 9. Theorems
Common data: o(Lambda, K) as in section 4, Gamma optimised (certified Gammahat_cert) or Gamma_expl.

**Theorem 1 (Route 1; NOT certified, not used).** For any x > 0: 1+eps1 = (1 + 1/x + kappa'/(2B+ x)) e^{B+(1+x)o}, M1 = 1 + eps1 + sqrt((1+eps1)^2 - e^{-eta}). Smallness: M1(Lambda=2, kappa'=2kappa, K=2kappa) e^eta < 2. Conclusion: (1+d1)^-1 omega_a <= omega_CY <= (1-d1)^-1 omega_a, d1 = M1(Lambda*, 0, 0) - 1, Lambda* from the t = 1 iteration. Explicit choice x = (B+ o)^(-1/2). Leading order d1 ~ 2(B+ o)^(1/4) ~ a^(3/2). v1's a0 values (0.0283 with D = 1.1, 0.0247 with D_L5) used the brief inputs and are withdrawn (section 0, historical note). GAP 4 (Rubinstein's 2C3) concerns this route only.

**Theorem 2 (Route 2).** Assume (I1)-(I6). The condition 4C0 > sigma = 2 kappa is **assumed for convenience**: it is not needed, because the upper root of section 5.2 decreases in E for either sign of 4 - sigma/C0 (REFEREE-CHAIN item 5); it holds for the certified inputs, and the certificate checks it on every a-interval. For any x > 0: 1+eps2 = e^eta (1 + x + sqrt((x-1)^2 + x delta))/(2x) e^{C0(1+x)o}; M2 = 1 + eps2 + sqrt((1+eps2)^2 - e^{-eta}). Smallness: M2 e^eta < 2, evaluated with Lambda = 2, diameter sqrt2 D and any K_t >= 2 kappa (the certificate: 4 kappa, section 7). Conclusion: e^{-eta} M2^-1 omega_a <= omega_CY <= M2 omega_a; symmetric form (1 + d2)^-1 omega_a <= omega_CY <= (1 - d2)^-1 omega_a with d2 = M2 e^eta - 1, M2 from the t = 1 iteration. Explicit choice x = 1 + sqrt(delta/(4 C0 o)) (any x > 0 is admissible; on an a-interval the certificate takes x = 1 + a xi with xi fixed per step, BUILD-WP5 section 8 item 7); x = 1 gives O(a^2).

**a0 criterion (v2: stated; REFEREE-WP4 item 6, REFEREE-CHAIN item 5).** a0 is the largest a such that, for every a' in (0, a]: (i) every input (I1)-(I6) is certified at a' (here: L5, WP2b, WP3 and WP5 cover a' <= 0.03); (ii) the continuity step closes, U < Lambda0 = 2; (iii) d2(a') < 1/2. For the certificate, (ii) and (iii) hold on all of (0, 0.03] (U <= 1.1694, d2 <= 0.121974), and (i) binds: **a0 = 0.03**.

**Certified instance:** the Theorem of section 0 (C_q, C_R, C0, D, e and Gamma as in sections 2 and 3.4; scaled exact rewrite of Theorem 2 in BUILD-WP5 3.2, refereed in REFEREE-WP5 item 2: CORRECT; at five points per interval including the upper end, where the minimum is attained, the exact unlinearised Theorem 2 stays below the certified bound, min ratio 1.0000000023, REFEREE-WP5 t2chk.py).

**Comparison with Kobayashi (38), (47) (wording per REFEREE-WP4 item 7).** In our units Kobayashi gets O(a), O(a^2) after Schauder, for an approximant with residual O(a^4). Our O(a^3) comes entirely from the approximant's a^8 residual: d ~ sqrt(a^-2 a^8). Explicit constants, O(a^3) from the sharper a^8 approximant; the qualitative rate of the same method on the same approximant is not beaten (applied to the same approximant, the qualitative method would give the same a^3 at the C^2 stage and better after a Schauder stage; REFEREE-WP4 item 7). The method is Kobayashi/Lye's with explicit constants and a Green/Harnack C^0 step instead of Moser iteration.

**Credit line (pasted verbatim from GATE-FINAL.md, "Honest wording we must use (credit lines)"; REFEREE-v4 item 6).**
"For the Kummer K3 surface of the Hurwitz (D4) lattice we give a computer-certified (interval arithmetic) two-sided bound (1+d)^-1 ω_a ≤ ω_CY ≤ (1−d)^-1 ω_a, valid for every bubble size a in (0, 0.03], with d(a) ≤ 2269.27 a^3 + o(a^3) and d(0.03) ≤ 0.122. To our knowledge (arXiv and Zenodo searched to 13 Nov 2026), this is the first explicit closeness bound with computed constants on a Kummer K3 surface. The qualitative estimate of this type is due to Kobayashi [1990, eq. (38), Thm 18], re-proved by Lye [2023, Thm 2.9, §6] and improved in weighted norms by Jiang [2025]. The convergent Green-operator representation and the leading Green-sum structure of the correction are due to Fu-Liu-Xiao-Yau [arXiv:2610.11077; Correction 2026-11-17: author list updated], who state that a certified value of a0 remained open. An explicit-interval bound of the same kind, proved by hand, exists for a special Schoen Calabi-Yau threefold (Fu-Liu-Xiao-Yau [arXiv:2610.21932], 0 < a ≤ 2^-2260). Certified Calabi-Yau metrics on the Fermat quartic, a different K3 at a fixed Kähler class, are due to Douglas-Gómez-Serrano-Lehmann-Platt-Qi [arXiv:2611.16908] and Liu-Zhang [arXiv:2611.15319]."

## 10. Input checks
- V_R = 1 confirmed; Remark V gives mu^2 - 16 pi^2 a^4 = 1 + O(a^8); V cancels in the chain.
- Diameter: the referee's D_R = 1.00 is the orbifold limit. L5 (v1 linear form, placeholder e; historical, not certified): D_L5(a) = sqrt(1+e)(1.0153 + 25.33a): 1.0406, 1.0913, 1.2686, 1.7753 at a = 1e-3, 3e-3, 1e-2, 3e-2. D = 1.1 is valid only for a <= 0.0033; the certificate never uses D = 1.1. The certificate uses L5's exact first expression with e = 1.37e-4: D(0.03) <= 1.77518 (BUILD-WP5 section 8 item 4; REFEREE-v3 item 8 (d): 1.775057 sqrt(1 + e)).
- eta: C_q = 17226.21 (WP3 bulk, the global maximum; the core is lower, 14416.38). v1's brief value 2.7e4 rounded the quintic-step 26873.4 and is superseded.
- kappa: C_R <= 1,731,732.5 (WP5 Stage A core). v1's brief value 2.7e6 is superseded.
- REFEREE-CHAIN item 7: WP5 control X4 compares the full F from the global potential at 616 core points per a (a = 0.03: |F|/a^8 = 14125.11 <= 14427.90 at control resolution); any constant offset between the core and bulk normalisations of F would have shown, and none is present.
- Independent float re-implementation of the chain with the WP5 inputs (REFEREE-CHAIN, scratch folder refchain.py): d2(0.01) = 2.851030e-3, d2(0.03) = 0.1219727, leading coefficient 2269.266, each at or below the certified value.

## 11. GAPs (v2 status)
1. Inputs C_q, kappa, B+, C0 (incl. ||dR||), e: **CLOSED** (CERTIFIED: WP2a, WP2b, WP3, WP5 Stage A). Residual risk: WP5 lemmas S1-S12 have had one referee pass (REFEREE-WP5); the float layer (IEEE round-to-nearest) in parts of WP2b, WP3, WP5.
2. C-infinity cut-off, C_q and kappa recomputed: **CLOSED** (section 8).
3. L5 holds only for a <= 0.03: **OPEN; it sets a0 = 0.03.**
4. Rubinstein's 2C3 not reconciled with constant 1: OPEN, Route 1 only (not used).
5. "opt" numbers in double precision: **CLOSED for Route 2**: the certificate evaluates Theorem 2 in arb with Gammahat_cert (WP5 Stage B).
6. Li-Yau's n: **CLOSED** (n = dim_R M, Li-Yau p.170; section 3.1).
STANDARD, not opened: heat kernel facts, Kahler normal coordinates, Yau/Calabi, Schauder/Fredholm; (R2 only, not used: Hoelder d dbar-lemma, mollification).
Closed in v1: Green constant; Chern-Lu constant; direction of the S4 sandwich; the path; continuity-method spaces; volume (cancels).

## 12. Possible improvements
- Sharper explicit heat-kernel lower bound (torus control shows 65x slack); d2 scales like sqrt(gamma).
- The remaining trace-to-eigenvalue square root needs a Schauder stage, not explicit today.
- L5 beyond a = 0.03 would raise a0 (smallness at 0.03 is 1.1694 < 2).
- The continuity step could use K_t = 2 kappa (U at a = 0.03: 1.161299 instead of 1.169325); this changes no final number, since d2 comes from the t = 1 iteration with K_t = 0.
