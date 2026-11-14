# BAR-WP5: assembly - certified explicit d(a) and a0 for the Kummer window (Route 2)

Written 2026-11-13 by the lead before any WP5 code exists. One-sided targets only (no rounded containment constants).

instruments: wp5.py

## Object
ANALYSIS-CHAIN.md Theorem 2 (Route 2) with Proposition C0, Lemma G (sec. 3) and the continuity closure (sec. 7), using the "expl" closed-form choices (x = 1 + sqrt(delta/(4 C0 o)), Gamma via Lemma G evaluated in interval arithmetic, not the double-precision "opt" grid), and L5 (DEFINITIONS-AND-LEMMAS-v2.md) for D(a, e). Inputs come ONLY from certified result files:
- wp2/result_wp2a.json (core leading constants, c6), wp2/result_wp2b.json (e, B+ = C0), wp3/result_wp3.json (outer C_out, R_out, e_out);
- plus a NEW core-remainder certificate inside wp5.py (stage A): the full |F|/a^8 and full |Ric|/a^6 in the core |zeta| <= 1 for every a in (0, 0.03] (WP2a certified only the a-independent leading profiles; the a^10-and-higher remainder must be bounded, or the gap reported).
e := max over all certified e values (WP2b and WP3 use different normalisations; take the max, the conservative choice). C_q := max(core full, outer) residual constant; kappa := max(core full, outer) Ricci constant.
Outputs: for every a in (0, 0.03] (a enclosed as an interval on sub-intervals): the smallness condition M2(2, 2 kappa) e^eta < 2 (certified true/false), d2(a) (certified upper bound), the leading coefficient lim_{a->0} d2(a)/a^3 (certified upper bound), and a0 := the largest a in (0, 0.03] up to which smallness is certified.

## Expected values (lead estimate, not certified)
Leading coefficient 8 D(0) sqrt(C_q Gammahat(0)) with D(0) = 1.0153, C_q = 17226.2, Gammahat(0) = 4.531: about 2270. Earlier double-precision chain (C_q = 2.7e4, kappa = 2.7e6): d2 = 3.57e-3 at a = 0.01, 0.158 at a = 0.03, a0 >= 0.03 (L5 range).

## Predictions
- V1: the core-remainder certificate gives core full |F|/a^8 <= 15000 and core full |Ric|/a^6 <= 1.9e6 for every a in (0, 0.03].
- V2: smallness is certified on every sub-interval of (0, 0.03], so a0 = 0.03 (L5's range is the binding limit).
- V3: certified upper bound for lim_{a->0} d2(a)/a^3 <= 2500.
- V4: certified d2(0.01) <= 3.2e-3 and d2(0.03) <= 0.15.
- V5: consistency: for every a the certified d2(a) is >= the double-precision "expl" formula re-evaluated with the same certified inputs (no certified value below its own float evaluation).

## Controls (any failure = VOID)
- X1: with the ANALYSIS-CHAIN brief inputs (C_q = 2.7e4, kappa = 2.7e6, B+ = C0 = 2/a^2, D_L5, Gamma from wp4), the interval code's float path reproduces the R2 "D_L5, expl" row of ANALYSIS-CHAIN.md sec. 0 at a = 1e-3, 3e-3, 1e-2 to 1e-3 relative.
- X2: Lemma G interval enclosure of Gammahat(0) contains the wp4 value 4.5308 (to its printed precision) or is below it (the interval bound may be sharper, never above the published closed form 4.5871).
- X3: mutation: doubling C_q must raise the certified leading coefficient by a factor in [1.35, 1.48] (about sqrt 2), and making e >= 1 must make the code refuse.
- X4: the core-remainder certificate dominates double-precision k1d/k1c-style samples of the full F and Ricci at >= 500 core points for a = 0.01 and 0.03.
- No NaN/inf; every input read from a result file is echoed with its source path and the seal id of its verdict.

## Verdict rule
ASSEMBLY CERTIFIED if V1-V5 hold and X1-X4 pass. A certified bound below a sampled or float value = VOID. If V1-V4 fail high, report the certified values (they are the result); no retuning of this BAR. If an input needed by the chain is not certified (e.g. the core remainder cannot be bounded), the verdict is GAP with the list, not HOLD.
Noise floor/ceiling: float evaluations are the floor of any certified upper bound.
Stopping tolerances (doctrine #116): any branch-and-bound or refinement tolerance in wp5.py must be justified against the headroom of V1-V4 (stated in BUILD-WP5.md).

## Seal fields
- formula_ref: ../ANALYSIS-CHAIN.md:99, ../ANALYSIS-CHAIN.md:69, ../DEFINITIONS-AND-LEMMAS-v2.md:296
- calibrated_on: cdf666b973, 4158aba4f8
- depends_on: 2d19816000, cdf666b973, 4158aba4f8
- null_not_symmetry: X3 doubles an input (changes the answer by about sqrt 2) and forces e >= 1 (must refuse); X1 is a known-answer control on the published float row.
- target_form: V1 core F in [0, 15000], core Ricci in [0, 1.9e6]; V2 boolean (a0 = 0.03); V3 coefficient in [0, 2500]; V4 d2(0.01) in [0, 3.2e-3], d2(0.03) in [0, 0.15]; V5 boolean.
- tolerances: any refinement tolerance must leave >= 2% headroom against V1-V4 (justified in BUILD-WP5.md).
- parameters: a sub-intervals >= 6 covering (0, 0.03]; sample a = 1e-3, 3e-3, 1e-2, 3e-2
- files: BAR-WP5.md, wp5.py
- command: memcap 8 python wp5.py --main result_wp5.json
