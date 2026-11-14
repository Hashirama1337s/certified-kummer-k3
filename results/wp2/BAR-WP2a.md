# BAR-WP2a: certified core profile constants (first certificate run)

Written 2026-11-13 by the lead before any WP2a code exists. Predictions use sampled double-precision values already on disk (disclosed below); the test is whether rigorous enclosures agree with them.

instruments: wp2a.py

## Object
DEFINITIONS-AND-LEMMAS-v2.md, section 4.1 (core zone, after Fix 1 normalisation) with the lead's decision: C-infinity step st(t) = h(t)/(h(t) + h(1-t)), h(t) = exp(-1/t), on the window [0.25a, a] (cut-off is 1 for |zeta| <= 0.25 and 0 for |zeta| >= 1; zeta = z/a). Quantities, all a-independent leading constants:
- P* = sup over sigma of |p(sigma)| (the one-variable profile in F8 = -(c6/2) H6(zeta/|zeta|) p(|zeta|^2));
- R* = sup over the core of the leading Ricci profile |Ric|_g / a^6 (R6 in the document);
- c6 enclosed as an interval from the lattice sum with a rigorous tail (Lemma A6).
Enclosures by interval/ball arithmetic (arb via python-flint, or mpmath.iv), with Lemma C for the cut-off derivatives near t = 0+ and 1-.

## Disclosed sampled values (double precision, not certified)
Window [0.25a, a], C-infinity step: sup|p| = 4.867, sampled sup R6 = 1.69e6. Window [0.5a, a] quintic (positive control): sup|p| = 9.260556, sup R6 = 2.6334e6 (refined 2.63343e6 at 0.75362). Window [0.5a, a] C-infinity: sup|p| = 8.714721, sup R6 = 9.3498e6. c6 = 15476.9364 (radius-7 sum 15476.936376).

## Predictions
- W1 (control, quintic [0.5a, a]): certified enclosure of sup|p| contains 9.260556 and has width <= 1e-6; certified upper bound of sup R6 lies in [2.6334e6, 2.66e6].
- W2 (c6): certified interval for c6 contains 15476.9364 and has width <= 1e-3.
- W3 (target, C-infinity [0.25a, a]): certified upper bound P* in [4.867, 4.92] (at most 1% above the sampled value).
- W4 (target): certified upper bound R* in [1.69e6, 2.0e6].
- W5 (consistency): sampled values never exceed the certified upper bounds (any sampled point above a certified bound = instrument bug = VOID).

## Controls (any failure = VOID)
- Interval code reproduces sup|p| for the quintic step in plain double precision to 1e-9 before switching to intervals.
- A mutated cut-off (derivative bound in Lemma C deliberately halved) must make the enclosure code report a violation of Lemma C on a test point.
- No NaN or infinite enclosure; every box of the cover reported.

Noise floor and ceiling: enclosure widths are set by the arithmetic precision (>= 128 bits) and the cover; the ceiling is the sampled value from below (true sup >= sampled), so a certified bound can only meet or exceed it.

## Verdict rule
CORE CONSTANTS CERTIFIED if W1-W5 hold and controls pass. If W3 or W4 fails high (certified bound far above sampled), the cover or the enclosure is too loose: report, do not retune silently. If a certified bound falls BELOW a sampled value, VOID (bug).
Not covered: ||dR|| (WP2b), bulk (WP3), neck remainders.

## Seal fields
- formula_ref: ../DEFINITIONS-AND-LEMMAS-v2.md:309-371, ../DEFINITIONS-AND-LEMMAS-v2.md:277
- calibrated_on: 4143022d85
- depends_on: 4143022d85
- null_not_symmetry: the quintic-window control uses a different cut-off profile, whose constants differ (9.26 vs 4.87); it is a known-answer control, not a symmetry.
- target_form: W1 interval containing 9.260556 with width <= 1e-6, R6 bound in [2.6334e6, 2.66e6]; W2 interval containing 15476.9364, width <= 1e-3; W3 [4.867, 4.92]; W4 [1.69e6, 2.0e6]; W5 boolean.
- parameters: window = [0.5a, a] quintic, [0.25a, a] C-infinity; Ewald eta = 6, 8 (c6 independence)
- files: BAR-WP2a.md, wp2a.py
- command: memcap 8 python wp2a.py --main result_wp2a.json
