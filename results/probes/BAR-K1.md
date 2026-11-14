# BAR-K1: residual order of the second-order Kummer approximant (probe)

Sealed before any run of k1.py. Written 2026-11-13 by the lead. No number below was computed; all are predictions from the hand argument in "Reasoning".

## Object
Torus C^2/Lambda, bubbles at all points of L = Lambda/2, all of Eguchi-Hanson size a. Approximant potential
Phi = mu|z|^2 + sum_{l in L} R(|z-l|^2) - (a^4/2) S(z), as defined in the docstring of k1.py (R = Eguchi-Hanson minus its two leading far-field terms; S = Ewald lattice sum of 1/r^2; mu = 1 + a^4 pi^2/(4V)). No cut-off functions. Residual F = log det(d dbar Phi), target 0.
- D4: Lambda = Hurwitz quaternions (complex structure = left multiplication by i).
- SQ (null): Lambda = (Z+iZ)^2.
Sample (fixed seed 7): 3000 uniform points of the unit cell of L (points closer than 0.1a to a bubble centre dropped) plus 11 shells x 200 directions; "core shell" = 200 points at |z| = a. sup = max over the sample, not a proved supremum.

## Reasoning (why the predictions)
Near a bubble the others contribute a^4 x (regular part of S). Its Hermitian quadratic part is absorbed by mu. The next terms are lattice-invariant harmonic polynomials h_d of degree d; the cross term with the Eguchi-Hanson metric is about a^4 (f/s - f') p q h_d / s, of size a^6 s^(d/2-2) in the core and a^8 s^(d/2-3) outside. The point group of D4 is W(F4), whose invariants have degrees 2, 6, 8, 12, so d >= 6 and the residual is O(a^8) everywhere. For SQ a degree-4 invariant of type (2,2) exists (|z1|^4 + |z2|^4 - 4|z1|^2|z2|^2), giving a^6 on the core shell. Bulk terms (bubble-bubble products, mismatch eps^2 with eps = a^4 pi^2/(4V)) are O(a^8) for both.

## Predictions (PASS/FAIL each)
- P1 (order, D4): r(a) = sup|F|/a^8 at a = 0.1, 0.05, 0.025 satisfies max r / min r <= 2.
- P2 (order, D4, local): log2( sup|F|(0.05) / sup|F|(0.025) ) in [7.5, 8.5].
- P3 (size, D4): sup|F| at a = 0.05 is <= 1e-5.
- P4 (core, D4): core-shell sup at 0.05 over core-shell sup at 0.025 has log2 in [7.5, 8.5].
- P5 (null, SQ core): the same core-shell log2 ratio for SQ lies in [5.5, 6.5].
- P6 (null, SQ global): SQ global log2( sup(0.05)/sup(0.025) ) < 7.5, OR SQ core-shell sup at 0.025 is at least 10x the D4 core-shell sup at 0.025. (Stated as an OR because the SQ a^6 term may be hidden globally by a larger a^8 bulk constant; the core-shell comparison is the discriminating one.)

## Controls (all must pass or the run is VOID)
- C1 exact zero: one bubble alone, sup|F| <= 1e-9.
- C2 Ewald split independence: eta = 6 vs 9, max |dF| <= 1e-3 x sup|F| + 1e-13.
- C3 trace: tr(d dbar S) = pi^2/V, relative error <= 1e-8, both lattices.
- C4 independent implementation (k1_mp.py: potential in 50-digit arithmetic, finite-difference Hessian), 4 points, D4, a = 0.05: |F_mp - F_numpy| <= 1e-3 |F_mp| + 1e-12.
- No NaN anywhere.

## Verdict rule
ROUTE ALIVE if P1-P4 hold and controls pass. If P1/P2 fail with exponent < 7.5, the a^8 claim is FAILED and the Kummer route is dropped unless the failure is traced to a control. P5/P6 failing alone means the symmetry explanation is wrong or incomplete (report as FAILED; route decision then rests on P1-P4 with that caveat).
What this probe does NOT show: nothing about the linearised inverse bound, and sup is sampled.

## Seal fields
- files: BAR-K1.md, k1.py, k1_mp.py
- command: python <lab-tools>/memcap.py 8 python k1.py result_k1.json ; then k1_mp.py result_k1.json c4_k1.json
