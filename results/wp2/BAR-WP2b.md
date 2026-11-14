# BAR-WP2b: certified curvature bounds and metric perturbation of the approximant (uniform in a)

Written 2026-11-13 by the lead before any WP2b code exists. Lesson from WP2a: no rounded constants as containment targets; all targets below are one-sided bounds or ranges.

instruments: wp2b.py

## Object
DEFINITIONS-AND-LEMMAS-v2.md (L1, L4(g), Lemma A7, section 4.1 with Fix 1 and Fix 4), the C-infinity cut-off on the window [0.25a, a] (as certified in WP2a), lambda-window radii 0.25a and a. For EVERY a in (0, 0.03] (a enclosed as an interval, Fix 2), certify:
- e(a) := sup over X of the g0-operator norm of g0^-1 E (the approximant's metric relative to the exact Eguchi-Hanson metric of the nearest bubble);
- ||dR||(a): the sup of the curvature perturbation (L4(g)) in the norm needed there;
- hence B+(a) := sup Bisec(omega_a) and C0(a) := -inf Bisec(omega_a) over orthogonal pairs, via the mirror of L4(g) (Bisec >= -(2/a^2 + ||dR||)/(1 - e)^2).
Report every bound in dimensionless form: e, a^2 B+, a^2 C0, a^2 ||dR||.

## Predictions (all for every a in (0, 0.03])
- B1: e(a) <= 1e-3.
- B2: a^2 ||dR|| <= 1e-2.
- B3: a^2 B+(a) <= 2.05 and a^2 C0(a) <= 2.05.
- B4 (sharpness): a^2 B+(a) >= 1.99 at a = 0.03 (lower bound from a point evaluation on the exceptional curve; the Eguchi-Hanson value is 2).
- B5: the bounds are monotone-free: they are computed with a as an interval on at least 3 sub-intervals covering (0, 0.03], including one touching a = 0 (via the scaled remainder bounds), and every sub-interval passes B1-B3.

## Controls (any failure = VOID)
- K1: with the lattice terms switched off (pure Eguchi-Hanson), the code certifies a^2 B+ in an interval containing 2 and a^2 C0 in an interval containing 2, e = 0 exactly.
- K2: point-sample consistency: at 2000 random points (all regions) at a = 0.01 and 0.03, double-precision bisectional maxima/minima (k1d/k1c-style evaluation) lie inside the certified bounds.
- K3: mutation: scaling the lattice correction by 10 must raise the certified e and ||dR|| by a factor >= 5.
- No NaN/inf; every box reported.

## Verdict rule
CURVATURE INPUTS CERTIFIED if B1-B5 hold and controls pass. A certified bound below a sampled value = VOID (bug). If B1-B3 fail, report the certified values; the analysis chain is re-run with them (no retuning of this BAR).
Noise floor/ceiling: one-sided bounds; sampled values give the floor.

## Seal fields
- formula_ref: ../DEFINITIONS-AND-LEMMAS-v2.md:309-371
- calibrated_on: 2d19816000
- depends_on: 2d19816000
- null_not_symmetry: control K1 removes the lattice terms entirely (pure Eguchi-Hanson, known exact curvature 2/a^2); K3 scales them by 10; both change the answer.
- target_form: B1 e in [0, 1e-3]; B2 a^2||dR|| in [0, 1e-2]; B3 a^2 B+, a^2 C0 in [0, 2.05]; B4 a^2 B+ lower bound in [1.99, 2.05] at a = 0.03; B5 boolean.
- parameters: a sub-intervals >= 3 covering (0, 0.03]; sample a = 0.01, 0.03
- files: BAR-WP2b.md, wp2b.py
- command: memcap 8 python wp2b.py --main result_wp2b.json
