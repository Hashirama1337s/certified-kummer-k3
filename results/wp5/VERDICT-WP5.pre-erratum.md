# VERDICT-WP5 (scored as sealed: claim 055009634d, BAR file 488ee8329c, amendment AMEND-WP5-1 b4c1dabedd sealed before code; instrument wp5.py pinned + sealed 891fafa2d9; controls on the sealed file (sha256 b3919b17... in controls_wp5.json; lab controls-check GO); launch-check passed; main run 1 invocation, 63 s, peak 0.77 GB)

## Certified result (computer-assisted; conditional on the unrefereed written chain, see below)
For every a in (0, 0.03] the continuity method closes (smallness M2 e^eta <= 1.1694 < 2 on all 30 a-intervals of width 0.001), and
  (1 + d2)^-1 omega_a <= omega_CY <= (1 - d2)^-1 omega_a  (symmetric form of Theorem 2),
with certified upper bounds:
- lim_{a->0} d2(a)/a^3 <= 2269.27 (and d2/a^3 <= 2275.00 for every a in (0, 1e-4]); = 8 D(0) sqrt(C_q Gamma), C_q = 17226.20 (WP3 bulk), D(0) = 1.015341, Gammahat(0) <= 4.530840.
- d2(a)/a^3 <= 2326.6 on (0, 0.001], <= 2850.9 on [0.009, 0.01], <= 4517.5 on [0.029, 0.03].
- d2(0.01) <= 2.851e-3; d2(0.03) <= 0.12197.
- a0 = 0.03 (L5's range is the binding limit; smallness is far from failing there).
Stage A (new, closes the core-remainder GAP of DEFINITIONS v3 / REFEREE-v3): for every a in (0, 0.03], full core |F|/a^8 <= 14416.37 and full core |Ric|/a^6 <= 1,731,732.4 (leading values 14123.85 and 1,695,905.9 plus all remainders, on 30 a-intervals). So the global residual constant is C_q = 17226.20 (bulk) and the global Ricci constant is <= 1.7317e6 (core).

## Scoring
| | sealed bar | measured | verdict |
|---|---|---|---|
| V1 | core full F <= 15000, core full Ricci <= 1.9e6 | 14416.37; 1,731,732.4 | HOLD |
| V2 | smallness on every sub-interval, a0 = 0.03 | 30/30 close; a0 = 0.03 | HOLD |
| V3 | lim d2/a^3 <= 2500 | 2269.27 | HOLD |
| V4 | d2(0.01) <= 3.2e-3; d2(0.03) <= 0.15 | 2.851e-3 (interval [0.009, 0.01]; the [0.01, 0.011] enclosure also covers a = 0.01 and gives 3.875e-3, both valid); 0.12197 | HOLD |
| V5 | certified d2 >= float re-evaluation with the same certified inputs | 30/30 OK, min ratio 1.0000000023 (lead's v5check2.py: certificate's own G1, G2, x, inputs, wp4 r2 formulas in double) | HOLD |

V5 note (disclosed): the lead's FIRST check (v5check.py) used wp4/chain.py's double-precision Gamma grid (Gamh_up rounds K D^2 up to the next grid value) and the point-value x choice, i.e. different, looser valid inputs; there all 30 certified values came out 0.0006%-0.11% BELOW that float. Re-run with the certificate's own Gamma values and x (the sealed "same certified inputs"), all 30 are at or above. The certified Gammahat(0) = 4.530840 is above the double minimiser 4.5308 (as a rigorous upper bound must be) and below the closed form 4.5871. Both scripts: scratch folder, copied to wp5/v5check.py, v5check2.py with logs.

Controls (controls_wp5.json, ALL_CONTROLS_PASS = true, 36.5 s): X1 PASS (brief inputs reproduce the "R2, D_L5, expl" row: rel 1.2e-4, 2.5e-5, 8e-6; vs wp4/explicit.py <= 9e-12). X2 PASS per AMEND-WP5-1 (4.530840 in [0.0702, 4.5871]; +4.0e-5 above the double 4.5308). X3 PASS (2C_q: 2269.278 -> 3209.244, ratio sqrt 2; e = 1 refused). X4 PASS (616 core points per a; a = 0.03: |F|/a^8 14125.11 <= 14427.90, |Ric|/a^6 1.695994e6 <= 1.766609e6). NaN/inf 8/8.

Verdict as sealed: ASSEMBLY CERTIFIED (V1-V5 HOLD, X1-X4 pass).

## What this does NOT yet claim (blocks any public "theorem")
- The written chain is not fully refereed: ANALYSIS-CHAIN.md (WP4) had one referee round (REFEREE-WP4) whose fixes are not all re-checked; BUILD-WP5's new lemmas S1-S12 and the scaled rewrite of Theorem 2 are NOT REFEREED; DEFINITIONS v3 was refereed (REFEREE-v3: all CORRECT / CORRECT WITH FIXES), v4 fixes pending.
- STANDARD facts not opened (heat kernel facts, Yau/Calabi, Schauder); ANALYSIS-CHAIN GAP 6 (Li-Yau's n).
- Floating-point layer assumptions (IEEE round-to-nearest in parts of WP2b, WP3, WP5 sphere covers).
- e: the certificate used max(WP2b 8.34e-5, WP3 1.37e-4) (conservative); effect 3e-5 relative on D.
- No clean-room rebuild yet; no blind replication of any certificate.
