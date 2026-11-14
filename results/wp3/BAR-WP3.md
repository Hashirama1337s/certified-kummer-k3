# BAR-WP3: certified residual and Ricci bounds outside the cores (neck + bulk, uniform in a)

Written 2026-11-13 by the lead before any WP3 code exists. One-sided targets only (no rounded containment constants).

instruments: wp3.py

## Object
DEFINITIONS-AND-LEMMAS-v2.md: the exact neck/bulk formula Q(z; eps) (L4(d)), the cancelled bulk form F8bulk (L4(e)) with its remainder, the Ricci formulas (L4(f)), Lemma A6 (Ewald tails; extend to the derivative orders actually used and state the extension), section 4.2 (bulk interface) and the neck part of 4.1. Region: all points at distance >= a from every bubble (i.e. outside the WP2a core window [0, a] in |zeta|), for EVERY a in (0, 0.03] (a enclosed as an interval, Fix 2), split as neck a <= |z| <= r_c and bulk |z| >= r_c (r_c = 0.1 or the builder's justified choice), quotient by the holomorphic + antiholomorphic point symmetries actually valid for F (group G_J of order 192, as WP2a found).
Outputs: C_out := sup over the region and over a of |F| / a^8, and R_out := sup of |Ric|_g / a^6, both certified upper bounds; plus a certified lower bound for sup |F8bulk| from a point evaluation.

## Disclosed sampled values (not certified)
sup |F8bulk| over 3000 cell points = 16283 (min -16283, max 13019); sup |Q|/a^8 on shells |z| = a, 2a, 4a at a = 0.025: 9040, 12391, 13210; K1c/K1d cell Ricci about 4.9e4 a^6.

## Predictions
- U1: certified C_out <= 18000 (at most about 10% above the sampled 16283).
- U2: certified lower bound for sup |F8bulk| >= 16283 (the sampled value).
- U3: certified R_out <= 1.0e6 (below the certified core Ricci 1.696e6 a^6, so the core sets the global Ricci constant).
- U4: the a-dependence: C_out computed on >= 3 a-sub-intervals covering (0, 0.03] (one touching 0 through the scaled remainder) differs between sub-intervals by at most 5%.
- U5: consistency: no sampled value exceeds a certified bound.

## Controls (any failure = VOID)
- J1: the code's double-precision F8bulk reproduces the 396-point check of DEFINITIONS (sup 14569.1 at distance >= 0.1) to 1e-6 relative.
- J2: Ewald independence: certified enclosures at Ewald parameters 6 and 8 overlap at 100 points.
- J3: mutation: dropping the cancellation term in the cancelled bulk form must make the enclosure blow up near |z| = r_c (or the code must flag non-convergence).
- J4: the volume identity: the certified cell average of F8bulk overlaps (2 pi^4 - Z_2/128)/V (Remark V), or report the discrepancy (diagnostic, not VOID).
- No NaN/inf; every box reported.

## Verdict rule
OUTER INPUTS CERTIFIED if U1-U5 hold and controls J1-J3 pass. Certified bound below a sampled value = VOID. If U1 or U3 fail high, report the certified values (they feed the chain anyway); no retuning.
Noise floor/ceiling: sampled values are the floor of any upper bound.

## Seal fields
- formula_ref: ../DEFINITIONS-AND-LEMMAS-v2.md:325-371, ../DEFINITIONS-AND-LEMMAS-v2.md:268
- calibrated_on: 2d19816000
- depends_on: 2d19816000, 4143022d85
- null_not_symmetry: J3 removes a term of the formula (changes the answer near the bubbles); J4 tests an independent identity.
- target_form: U1 C_out in [0, 18000]; U2 lower bound in [16283, inf); U3 R_out in [0, 1.0e6]; U4 relative spread in [0, 0.05]; U5 boolean.
- parameters: a sub-intervals >= 3 covering (0, 0.03]; Ewald eta = 6, 8
- files: BAR-WP3.md, wp3.py
- command: memcap 8 python wp3.py --main result_wp3.json
