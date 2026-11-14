# BAR-S2: Stage 2 probe 1+2 - sharper diameter (L5') and heat-kernel constant (Lemma G'), re-assembled on (0, 0.03]

Written 2026-11-14 by the lead before any S2 code exists. One-sided targets. Uses only the existing certified input files (wp2/result_wp2a.json, result_wp2b.json, wp3/result_wp3.json, wp5/result_wp5.json Stage A values); no input re-runs.

instruments: s2.py

## Object
stage2/STAGE2-DESIGN.md: Lemma L5' (D(a, e) = inf over 0 < rho <= 1/sqrt32 of sqrt(2(1+e)) [2 rho + (1 + a^4/rho^4)^(1/4) (1/sqrt2 + 2(pi-2) rho)], valid for a <= 0.05) and Lemma G' (Doeblin tail replaced by the spectral tail B(t0) e^{-lambda1 (t - t0)} with lambda1 >= pi^2/(D^2 cosh^3((D/2) sqrt(K/3))) from Andrews-Clutterbuck 2013 Thm 2). Outputs (arb, certified): Gammahat'(0) and Gamma'(K, D) as used by Prop C0; D'(a, e) on (0, 0.03]; the Route 2 assembly (ANALYSIS-CHAIN-v2.2 Theorem 2, WP5 Stage B engine) re-run on the same 30 a-intervals with L5' and Gamma' replacing L5 and Gamma.

## Expected (CHECK, double precision, design seat)
Gammahat'(0) ~ 0.5285; D(0.03) ~ 1.3220, D(0) -> 1; limsup d2/a^3 ~ 763; d2(0.01) ~ 8.7e-4; d2(0.03) ~ 0.029.

## Predictions
- S1: certified Gammahat'(0) <= 0.60.
- S2: certified D'(0.03, e) <= 1.33 and D'(a, e) <= 1.0001 as a -> 0 (on the first a-interval's lower end a = 0).
- S3: re-assembly on (0, 0.03]: smallness holds on all 30 intervals; certified limsup d2/a^3 <= 850; certified d2(0.03) <= 0.035; certified d2(0.01) <= 1.0e-3.
- S4: on every interval the certified d2 is >= the double-precision re-evaluation of Theorem 2 with the same certified inputs (same Gamma', same x, same D').

## Controls (any failure = VOID)
- Z1: brute-force count: over >= 1e5 random segments of length <= 1/sqrt2 in R^4, the number of open balls B(l, rho) (l in L, rho = 1/sqrt32) crossed is <= 2 (the L5' packing claim).
- Z2: flat-torus control: Lemma G' applied to the flat torus R^4/Z^4 (V = 1, D = 1, K = 0) gives an upper bound >= 0.0702 (the sharp -min G of ANALYSIS-CHAIN 3.6).
- Z3: known answer: with the spectral tail switched off (old Doeblin tail) the engine reproduces the certified Gammahat(0) <= 4.530841 to 1e-5 relative, and with L5' replaced by the old L5 the assembly reproduces wp5/result_wp5.json d2 on [0.029, 0.03] to 1e-6 relative.
- Z4: mutation: halving lambda1 in the spectral tail must raise Gammahat'(0).
- No NaN/inf; controls on the final sealed file (sha256 recorded; lab controls-check GO).

## Verdict rule
STAGE-2 PROBE CERTIFIED if S1-S4 hold and Z1-Z4 pass. Certified below a float = VOID. S1-S3 failing high: report the certified values; no retuning. If the Andrews-Clutterbuck statement cannot be confirmed (first NONZERO eigenvalue of a closed manifold with Ric >= -K), verdict GAP for Lemma G' and report L5' alone.

Noise floor/ceiling: float re-evaluations are the floor of every certified upper bound; the incumbent is the current certificate (WP5: limsup 2269.27, d2(0.03) <= 0.121974), which S3 must beat.

## Seal fields
- formula_ref: STAGE2-DESIGN.md:1-154, ../ANALYSIS-CHAIN-v2.2.md:140
- calibrated_on: 055009634d
- depends_on: 055009634d, cdf666b973, 4158aba4f8
- null_not_symmetry: Z3 switches the new ingredients off (known answers from WP5); Z4 halves lambda1 (changes the answer).
- target_form: S1 [0, 0.60]; S2 [0, 1.33] and [0, 1.0001]; S3 boolean, [0, 850], [0, 0.035], [0, 1.0e-3]; S4 boolean.
- tolerances: time-grid and rho-grid refinements must leave >= 5% headroom against S1-S3 (justified in BUILD-S2.md).
- parameters: the 30 WP5 a-intervals of width 0.001 covering (0, 0.03]; rho grid for L5' inf; fixed - inputs frozen to the certified files
- files: BAR-S2.md, s2.py
- command: memcap 8 python s2.py --main result_s2.json
