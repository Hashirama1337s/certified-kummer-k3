# BAR-REP3: blind replication of the core (patch-interior) certificate, WP5 Stage A

Written 2026-11-13 by the lead before the blind seat starts. Seat fenced from wp2/, wp3/, wp5/, DEFINITIONS v3, REFEREE-v3, LOG, WHERE-ELSE (fence 50ebed45 already active for lab:blind-replicator).

instruments: rep3/rep3.py

## Object
The full residual F = log(omega_a^2 / omega_flat^2) and the full Ricci form of the Kummer approximant inside one bubble core |z| <= a (zeta = z/a, |zeta| <= 1), with the C-infinity cut-off st(t) = h(t)/(h(t)+h(1-t)), h(t) = exp(-1/t), on the window [a/4, a] (the odd-odd part of the regular potential is removed with weight 1 - st((|z|/a - 1/4)/(3/4)) ... as specified in the brief), from DEFINITIONS-AND-LEMMAS-v2.md (sections 1-2, L2, L3, L4 (a)-(d), 4.1) only. For a in [0.029, 0.03] (a as an interval): (i) a double-precision sampled max of |F|/a^8 and |Ric|_g/a^6 over >= 20000 core points in the resolved chart plus refinement; (ii) an independent rigorous upper bound for both on |zeta| in [1/4, 1] (the window zone Z2) by any method the seat chooses (different from a leading-term-plus-remainder split if possible).

## Lead's certified values (NOT shown to the seat)
WP5 Stage A on [0.029, 0.03]: core |F|/a^8 <= 14416.38, core |Ric|/a^6 <= 1,731,732.5. WP2a leading profile values: sup|F8| <= 14123.85, R6* <= 1,695,905.922 (these are a-independent leading constants).

## Predictions
- P5: the blind sampled max of |Ric|_g/a^6 over the core at a = 0.03 is <= 1,731,732.5 and >= 1.69e6.
- P6: the blind sampled max of |F|/a^8 over the core at a = 0.03 is <= 14416.38 and >= 14100.
- P7: the blind rigorous upper bound for |Ric|_g/a^6 on Z2, a in [0.029, 0.03], lies in [1.69e6, 2.0e6].

## Controls (any failure = VOID)
- Q4: with the lattice terms switched off (pure Eguchi-Hanson) the blind code gives F = 0 and Ric = 0 to rounding.
- Q5: the blind sampled values never exceed the blind rigorous bound (P7 region).

## Verdict rule
REPLICATED if P5-P7 hold and Q4-Q5 pass. P5 or P6 with a blind sample ABOVE the WP5 certified bound = WP5 Stage A is wrong (VOID the chain, find the bug). P7 outside its window: report the formulation difference, no retuning.

## Seal fields
- formula_ref: ../DEFINITIONS-AND-LEMMAS-v2.md:309-371
- calibrated_on: 055009634d
- depends_on: 055009634d, 2d19816000
- null_not_symmetry: Q4 removes the lattice terms (answer must become 0).
- target_form: P5 [1.69e6, 1731732.5]; P6 [14100, 14416.38]; P7 [1.69e6, 2.0e6].
- tolerances: sampling refinement stops at 1e-4 relative agreement of successive maxima (P5 window lower slack 2%, upper is the certified bound itself).
- parameters: a = 0.03 (samples), a in [0.029, 0.03] (bound); fixed - the largest a is where the remainders are biggest
- files: BAR-REP3.md, rep3/rep3.py
- command: memcap 8 python rep3/rep3.py
