# VERDICT-WP2a (scored as sealed: claim 2d19816000, BAR file ba3dc16c79, instrument pinned + sealed 9f36f07d23; launch-check passed; run 213 s, peak 1.02 GB)
Certified enclosures (arb 128 bits; python-flint 0.9.0):
- c6 = [15476.9365415835862485180, 15476.9365415835862485185] (Ewald eta 6 and 8 and a direct radius-20 sum all overlap)
- quintic [0.5a, a]: P* in [9.26055627355716, 9.26055627356468]; R* x c6 in [2633433.2295, 2633433.2558]
- C-infinity [0.25a, a]: P* in [4.86705830003958, 4.86705830003983]; R* x c6 in [1695905.904, 1695905.922]
- consequence: core sup|F8| <= (c6/2)(3/8)P* = 14123.85 a^8 (C-infinity window)
| | sealed bar | measured | verdict |
|---|---|---|---|
| W1 quintic control | P* interval contains 9.260556, width <= 1e-6; R* bound in [2.6334e6, 2.66e6] | width 7.5e-12 OK; R* 2633433.26 OK; but 9.260556 (a ROUNDED constant) is not inside [9.2605562736, 9.2605562736] | FAILED (literal containment) |
| W2 c6 | contains 15476.9364, width <= 1e-3 | width 4e-19 OK; 15476.9364 not inside (true 15476.93654) | FAILED (literal containment) |
| W3 C-inf P* upper | in [4.867, 4.92] | 4.8670583 | HOLD |
| W4 C-inf R* upper | in [1.69e6, 2.0e6] | 1.6959e6 | HOLD |
| W5 sampled <= certified | boolean | all sampled values below the certified bounds (tightest: quintic R6 2.63343e6 vs 2633433.26) | HOLD |
Controls: C1 3.6e-14 PASS; C2 quintic axis values reproduced; C3 G_J invariance 4e-13 (W-quotient would be wrong: 323%); C4 halved st-derivative enclosures detected (50/65 vs mpmath; 2297/1747 R6 and 286/261 p violations in the pipeline; 0 unmutated); C5 NaN/inf guards PASS. Note: halving the Lemma C END-PIECE bound is undetectable by point tests because that bound is >= 4.63x loose (a property of the lemma, recorded).
Verdict as sealed: FAILED (W1, W2), because I sealed rounded constants as containment targets (doctrine #98-99: never a rounded constant; the lint had asked for exact forms). Substance: every certified enclosure agrees with every sampled value to the precision the sampled values had; the disclosed c6 = 15476.9364 was itself off by 1.4e-4 (finite-radius sums). The certified enclosures above are valid regardless of this BAR's predictions and are what the proof uses.
Weaknesses carried forward (BUILD-WP2a.md): the R6 invariant formula is checked numerically (7e-14), not proved symbolically; G_J quotient uses Lemma Q (numerically checked); covers the shell r1 <= |zeta| <= r2 only; DEFINITIONS-v2 needs the [0.25, 1] window edits and Lemma A6 extended to order 6.
