# VERDICT-REP (scored as sealed: claim 700278e11f, BAR file 6c8fed2f34; instruments rep1/rep1.py (49f3ac66...) and rep2/rep2.py (71714dff...) pinned, then re-run by the lead after the pin: rep1/lead_rerun.log, rep2/lead_rerun.log reproduce the seats' numbers)

Two fenced blind seats (fence 50ebed45: no wp2/, wp3/, wp5/, BAR-REP, DEFINITIONS v3, REFEREE-v3, LOG, WHERE-ELSE, memory). Both report no forbidden path opened. Disclosed exposures: R1's spec (ANALYSIS-CHAIN) shows the leading form 8 D sqrt(C_q gamma) a^3 and the old table; R2's spec (DEFINITIONS v2) shows c6 and the old sampled bulk sup 16283 (a floor, below the target). R1 cross-read wp4/explicit.py and one Li-Yau page AFTER its results; R2 read the probe1 point generator after its enclosure (for Q2 only).

| | sealed bar | blind result | lead's certified value | verdict |
|---|---|---|---|---|
| P1 | R1 leading coefficient in [2200, 2340] | 2269.2912 (Lemma G, own algorithm: direct k, dyadic time nodes); 2283.3351 with closed-form Gamma | 2269.27 | HOLD |
| P2 | R1 d2(0.03) in [0.10, 0.14], smallness holds | 0.12197295 (G) / 0.12271502 (expl); smallness 1.1613 < 2 | 0.121974; 1.1694 | HOLD |
| P3 | R2 enclosure of F8bulk(1/4,0,0,0) overlaps [-16918.4629, -16918.2833] | [-16918.37305985292092614861 +/- 4.88e-21] (Ewald eta 4 and 9, own |x|^-8 Ewald split) | WP3 [-16918.4629, -16918.2833] | HOLD |
| P4 | R2 sampled sup |F8bulk| over V0 in [16900, 17226.20] | 16918.37306 (attained at all 24 half-root points l/2; 32768 Sobol + refinement) | certified sup <= 17226.20 | HOLD |
| Q1 | R1 reproduces table row R2 D_L5 expl at 1e-3, 1e-2 to 1e-2 | 2.93066e-6, 3.59404e-3 vs 2.931e-6, 3.594e-3 | - | PASS |
| Q2 | R2 reproduces 396-point check 14569.1 to 1e-4 | 14569.1306 (rel 2.1e-6) | - | PASS |
| Q3 | samples never exceed own certified bounds | R1 n/a (no sampling); R2 samples <= its point enclosure region and below 17226.20 | - | PASS |

Verdict as sealed: REPLICATED.

Notes: the replicated quantities are the two load-bearing numbers of the chain (the assembly arithmetic d(a) from certified inputs, and the bulk residual peak that sets C_q). Not replicated by this BAR: WP2b curvature certificate, WP3's upper-bound cover (R2 checks the peak value and samples, not the cover), WP5 Stage A (separate BAR-REP3, seat running).
Side findings (waypoints, not claims): F8bulk at the 24 deep holes equals Z2 to 1.45e-31 (explained by D4 triality: coset sum = 5 Z for |x|^-8; WHERE-ELSE.md); sup of F8bulk itself is not attained in V0 (limit (15/16) c6 = 14509.628 at the lattice point), sup |F8bulk| = 16918.373 at l/2.
