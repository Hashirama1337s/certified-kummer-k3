# BAR-K1c-AMEND1 (control bars only; written and sealed before any --main run; no prediction R1-R7 changes)

Amends BAR-K1c.md (seal 214e4b3168) and SPEC-K1c.md (file seal 59ed9e745f). Reason: the build seat (BUILD-K1c.md) ran the controls; four fail as specified. Lead's review of each:

1. K-EH expected values at |z| = 0.5a and a were WRONG in SPEC-K1c.md: 1.822966 and 0.696753 were copied from the referee's sampled maxima (REFEREE-SKELETON.md, eh_curv.py), which underestimate. Exact maximum bisectional curvature of Eguchi-Hanson (potential s phi' = f, f = sqrt(s^2 + a^4)) is 2a^4/f^3. Lead hand check: s = 0 gives 2/a^2; |z| = 0.5a gives 2/(1.0625)^(3/2) = 1.826151; |z| = a gives 2/2^(3/2) = 0.707107. Builder's analytic path and independent coordinate path agree with this to <= 1e-12. NEW BAR: BIS x a^2 at 0.1a, 0.5a, a equal 2/(1+t^4)^(3/2) for t = 0.1, 0.5, 1 to within 1e-6 relative.
2. K-K1b: k1b.py truncates lattice sums on a ball centred at the origin; k1c.py centres on the nearest bubble (periodic, consistent with delta). NEW BAR: on k1b's own ball, agreement <= 1e-12 absolute (builder: 1.0e-15); the nearest-centred difference is recorded, not barred.
3. K-MATCH F: per-point relative error is ill-defined where F nearly vanishes. NEW BAR: |F_A - F_B| <= 1e-6 x max|F| over the matching sample, per (lattice, a) (builder: <= 1.5e-7). g-matching bar unchanged (1e-6 relative; builder <= 7e-14).
4. Truncation of the polynomial model: unchanged bar (<= 1e-6 relative at r2; builder worst 2.8e-7, D4 a = 0.1 Hessian).

Prediction R4 (BIS x a^2 in [1.90, 2.05]) is unchanged; note for the record that the sample's innermost shell is 0.1a, where the exact Eguchi-Hanson value is 1.99970.
Process note: the SPEC's control values should have come from a formula, not from a sampled maximum (doctrine: target_form never a rounded/sampled constant).
