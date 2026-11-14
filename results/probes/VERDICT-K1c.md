# VERDICT-K1c (scored as sealed: BAR-K1c 214e4b3168 + control amendment b1abf05f8a sealed before the run; instrument pinned 11fb62d874; launch-check passed)
Controls under the amended bars: K-EH exact (1.99970, 1.82615, 0.70711 = 2/(1+t^4)^1.5 to <= 1e-12); K-FS PASS; K-JET PASS (ddF 1.6e-10 rel, curvature 1e-14); K-MATCH F vs max <= 1.5e-8, g <= 7e-14; K-K1b on k1b's ball 1.0e-15; model truncation worst 2.8e-7; no NaN. All PASS.
| | sealed bar | measured (a = 0.1, 0.05, 0.025) | verdict |
|---|---|---|---|
| R1 D4 sup|F|/a^8 <= 2.5e4, spread <= 2 | | 1.378e4, 1.452e4, 1.457e4 (spread 1.06) | HOLD |
| R2 D4 Ricci <= 1e6 a^6 | | 7.1e7, 4.9e7, 4.4e7 a^6 (40-70x over), all inside the cut-off annulus | FAILED |
| R3 D4 Ricci order, log2 in [5.5, 6.5] | | 6.17 | HOLD |
| R4 D4 max bisectional x a^2 in [1.90, 2.05] | | 2.708, 2.0004, 1.9997 | FAILED (a = 0.1) |
| R5 D4 positivity >= 0.9 | | 0.992, 0.9994, 1.0000 | HOLD |
| R6 uncut blow-up: Ricci at 0.01a >= 10 x at 0.1a (D4, 0.05) | | 6.00 vs 0.0600 = 100x (like 1/s) | HOLD |
| R7 SQ null: sup|F| >= 10 x D4 at 0.025 | | 6.88e6 vs 1.46e4 a^8 = 472x | HOLD |
SHAPE FIXED condition (R1-R5) NOT met: the cut-off design (lambda = 0.1) is FAILED and goes back to the lead.
What held: the residual stays ~1.46e4 a^8 with the cut-off, the Ricci order is a^6 as predicted, positivity holds, and the referee's blow-up diagnosis is confirmed (uncut Ricci grows exactly like 1/s, i.e. 100x from 0.1a to 0.01a).
What failed: the Ricci constant. Cutting INSIDE the Eguchi-Hanson core (0.1a-0.2a) is expensive because the core metric is far from flat there; the BAR lint had warned that lambda appeared at one value only (lesson 108).
Post hoc (not tested): a cut-off OUTSIDE the core (lambda about 1-2, where the metric is nearly flat) should cost about c6 lambda^4 a^8 in F and about c6 lambda^2 a^6 in Ricci, with c6 about 1.5e4, i.e. a few 1e4 a^6 instead of 4e7. Needs a sealed K1d with at least two lambda values.
