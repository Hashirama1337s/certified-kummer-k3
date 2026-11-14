# BAR-REP: blind replication of the two load-bearing numbers of the K3 certificate

Written 2026-11-13 by the lead before either blind seat starts. Two blind seats, each building from the spec only (fenced from wp2/, wp3/, wp5/ code and results, from DEFINITIONS v3, REFEREE-v3, LOG, WHERE-ELSE and memory).

instruments: rep1/rep1.py, rep2/rep2.py

## Object
R1 (assembly): from ANALYSIS-CHAIN.md (Theorem 2, Prop C0, Lemma G, continuity closure) and DEFINITIONS-AND-LEMMAS-v2.md L5, with the certified input constants given in the brief (C_q = 17226.2012683894, kappa = 1731732.44 a^6 for a <= 0.03 (core; a-interval values allowed if the seat re-derives nothing), a^2 C0 = a^2 B+ <= 2.0000352761032563, e <= 1.368159761580671e-4), an independent rigorous (interval) evaluation of: the leading coefficient lim_{a->0} d2(a)/a^3; d2 at a = 0.03; whether the smallness condition holds at a = 0.03.
R2 (bulk residual): from DEFINITIONS-AND-LEMMAS-v2.md (L4(e) cancelled bulk form F8bulk, Ewald sums, Lemma A6) only: a rigorous enclosure of F8bulk at the point z = (1/4, 0, 0, 0) and a double-precision sampled sup of |F8bulk| over the Voronoi cell V0 (>= 20000 points plus local refinement).

## Lead's certified values these must match (NOT shown to the seats)
R1: lim <= 2269.27 (certificate, Gammahat(0) <= 4.530840); d2(0.03) <= 0.12197; smallness 1.1694 < 2 at [0.029, 0.03].
R2: F8bulk(1/4,0,0,0) in [-16918.4629, -16918.2833] (WP3, eta 6 and 8); certified sup over V0 <= 17226.20; certified point lower bound 16918.28.

## Predictions
- P1: R1 blind leading coefficient (upper bound) lies in [2200, 2340] (an independent valid choice of Lemma G evaluation or of x may differ from 2269.27; the closed-form Gamma_expl would give about 2283).
- P2: R1 blind d2(0.03) (upper bound) lies in [0.10, 0.14], and smallness holds at a = 0.03.
- P3: R2 blind enclosure of F8bulk(1/4,0,0,0) overlaps [-16918.4629, -16918.2833].
- P4: R2 blind sampled sup |F8bulk| over V0 lies in [16900, 17226.20] (at least the certified point lower bound region, never above the certified upper bound).

## Controls (any failure = VOID)
- Q1 (R1): with the brief inputs of ANALYSIS-CHAIN sec. 2 (C_q = 2.7e4, kappa = 2.7e6 a^6, B+ = C0 = 2/a^2, D_L5 linear form, e = 140 a^4) the blind code reproduces the table row "R2, D_L5, expl" at a = 1e-3 and 1e-2 to 1e-2 relative.
- Q2 (R2): the blind F8bulk reproduces the DEFINITIONS-v2 396-point check value (sup 14569.1 at distance >= 0.1) to 1e-4 relative.
- Q3: each seat's sampled values never exceed its own certified bounds.

## Verdict rule
REPLICATED if P1-P4 hold and Q1-Q3 pass. P3 failing (disjoint enclosures) = a real discrepancy: stop, find the bug before any further claim. P4 above 17226.20 = the WP3 certificate is wrong (VOID the chain). P1/P2 outside their ranges: investigate the formulation difference; report, no retuning.
Noise floor/ceiling: R2 samples are a floor for the true sup; certified bounds are the ceiling.

## Seal fields
- formula_ref: ../ANALYSIS-CHAIN.md:99, ../DEFINITIONS-AND-LEMMAS-v2.md:296
- calibrated_on: 055009634d, 4158aba4f8
- depends_on: 055009634d, 4158aba4f8
- null_not_symmetry: Q1 changes the inputs to the published brief values (known answer); Q2 is the independent 396-point check.
- target_form: P1 [2200, 2340]; P2 [0.10, 0.14] and boolean; P3 interval overlap with [-16918.4629, -16918.2833]; P4 [16900, 17226.20].
- tolerances: R2 sampling refinement stops when successive local maxima agree to 1e-3 relative (headroom: P4 window is 1.9% wide); R1 needs no refinement.
- parameters: a = 0.03 and a -> 0 (R1); z = (1/4,0,0,0), >= 20000 cell samples (R2)
- files: BAR-REP.md, rep1/rep1.py, rep2/rep2.py
- command: memcap 8 python rep1/rep1.py; memcap 8 python rep2/rep2.py
