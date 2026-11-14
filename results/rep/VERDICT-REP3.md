# VERDICT-REP3 (scored as sealed: claim 93b705fabc, BAR file fa3e2bba0f; instrument rep3/rep3.py pinned after delivery - disclosed: results were produced by the seat before the pin; the seat reports a from-scratch rerun on the final instrument byte-identical to its cache; the lead did NOT rerun the 2008-s interval stage)

Fenced blind seat (fence 50ebed45). Disclosed exposure: DEFINITIONS v2 prints sampled values for this window (sup|p| 4.867, sup R6 1.69e6 "not refined", c6, F10(0)); not used as targets.
Independent method: exact Taylor polynomial of psi (degree 18; 14/12 rigorous) from a pairing formula for |z-l|^-2nu coefficients, jet arithmetic, exact 2x2 log det identity, resolved chart near E, polar Taylor-model boxes for the rigorous stage.

| | sealed bar | blind result | WP5 certificate | verdict |
|---|---|---|---|---|
| P5 | sampled max core |Ric|_g/a^6 at a = 0.03 in [1.69e6, 1731732.5] | 1.695996527e6 (axis directions, |zeta| = 0.883040; 9 digits trusted, 100-digit mpmath check) | <= 1,731,732.5; WP5's own X4 sample 1.695994e6 | HOLD |
| P6 | sampled max core |F|/a^8 at a = 0.03 in [14100, 14416.38] | 1.412511755e4 ((e_i +- e_j)/sqrt2, |zeta| = 0.743108) | <= 14416.38; X4 sample 14125.11 | HOLD |
| P7 | blind rigorous |Ric|_g/a^6 on 1/4 <= |zeta| <= 1, a in [0.029, 0.03], in [1.69e6, 2.0e6] | <= 1.746877e6 (507,037 boxes, 0 unresolved, queue bounded; also |F|/a^8 <= 1.483137e4) | <= 1,731,732.5 (tighter; both are upper bounds) | HOLD |
| Q4 | pure EH gives F = Ric = 0 | det g0 - 1 <= 2.06e-13; ddbar log det <= 3.6e-11; mpmath raw F -3.6e-42 | - | PASS |
| Q5 | samples <= own bound | 1.695996527e6 <= 1.746877e6; 14125.12 <= 14831.37 | - | PASS |

Verdict as sealed: REPLICATED.
Agreement: the blind sampled Ricci maximum 1,695,996.5 and WP5's independent X4 sample 1,695,994 agree to ~1.5e-6 relative; the blind sampled |F| maximum 14125.12 matches WP5's 14125.11. Both independent rigorous bounds sit above every sample.
Side finding: on E the refined Ricci maximum is 416,756.5 a^8 (a = 0.03), 0.15% above the spec's sampled "about 4.16e5" (irrelevant to the window).
