# BAR-K1d: remove only the odd-odd part, outside the core (probe)

Written 2026-11-13 by the lead before k1d.py exists. Informed by the settled K1b and K1c numbers (disclosed); no K1d quantity has been computed.

instruments: k1d.py

## Object and sample: SPEC-K1d.md sections 1, 2 and 5. sup = max over the sample, not a proved supremum.

## Reasoning
K1c showed the uncut odd-odd terms make Ricci grow like 1/s (exactly 100x from 0.1a to 0.01a) and that a cut-off inside the core costs about 4e7 a^6 in Ricci. Here only the odd-odd part (leading term about 7.7e3 a^4 |z|^6) is removed, at radius lambda a where the metric is close to flat. Hand estimate: the change in F is about the d d-bar of that term at |z| ~ lambda a, i.e. of order 1e4-1e6 lambda^4 a^8, and in Ricci one more 1/(lambda a)^2, i.e. of order 1e4-1e6 lambda^2 a^6. I do not trust my constants (K1c's were off by 50x), so the bars below are about ORDER, NO BLOW-UP and RELATIVE IMPROVEMENT, with generous size bars.

## Predictions (D4; each for lambda = 0.5 and lambda = 1 separately)
- D1 (Ricci order): log2( RIC(0.025)/RIC(0.0125) ) in [5.5, 6.5], and log2( RIC(0.05)/RIC(0.025) ) in [5.5, 6.5].
- D2 (improvement over K1c): RIC/a^6 at a = 0.025 <= 4.4e6 (10x below K1c's 4.4e7).
- D3 (no blow-up at the exceptional curve): at every a, RIC on the 0.01a shell <= 2 x RIC on the 0.1a shell.
- D4 (residual): sup|F| / a^8 <= 1e6 at every a.
- D5 (curvature): BIS x a^2 in [1.90, 2.05] at every a.
- D6 (positivity): POS >= 0.9 at every a.

## Controls (any failure = VOID): K-ZERO, K-Q, K-SMOOTH, K-MODEL as in SPEC-K1d.md section 4; model truncation at rA <= 1e-8 relative; no NaN.
Noise floor: Ricci jets have absolute noise of order 1e-12/(lambda a)^2 at most, below 1e-6 of any predicted RIC.
Baseline to beat: K1c (lambda = 0.1, cut of the whole regular part): RIC 4.4e7 a^6.

## Verdict rule
SHAPE FIXED for a given lambda if D1-D6 all hold for that lambda and the controls pass. If neither lambda passes, the cut-off line of attack goes back to the lead with no third retune without the owner's say.
Not shown: odd-odd terms of degree >= 12 are not removed (estimated relative size below 1e-10 at |z| <= 2a; a proof must remove them exactly); proved suprema; uniformity beyond three values of a.

## Seal fields
- formula_ref: SPEC-K1d.md:3-9, SPEC-K1c.md:12-16
- calibrated_on: 214e4b3168
- depends_on: 214e4b3168, ff0c52846b
- null_not_symmetry: the baseline K1c cuts the whole regular part inside the core; K1d cuts only the odd-odd part outside it. The two objects differ by smooth even-even terms and by the cut radius, which changes the Ricci constant, not a symmetry.
- target_form: D1 in [5.5, 6.5]; D2 RIC/a^6 in [0, 4.4e6]; D3 ratio in [0, 2]; D4 in [0, 1e6]; D5 in [1.90, 2.05]; D6 >= 0.9. Unrounded doubles.
- parameters: lambda = 0.5, 1; a = 0.05, 0.025, 0.0125
- files: BAR-K1d.md, SPEC-K1d.md, k1d.py
- command: memcap 8 python k1d.py --main result_k1d.json
