# VERDICT-K1b (scored exactly as sealed, seal ff0c52846b; launch-check passed before the run)
Controls: CK worst relative difference 2.6e-10 (bar 1e-3) PASS; C2 1.6e-13 (bar 7.1e-10) PASS; C4 4/4 PASS (worst: 0.01a shell, 5.6e-13 against a bar of 2.65e-12); no NaN.
| | sealed bar | measured | verdict |
|---|---|---|---|
| Q1 D4 max r / min r (r = sup/a^8) | <= 2 | 1.045 (1.493e4, 1.555e4, 1.560e4) | HOLD |
| Q2 D4 global log2 ratio 0.05 -> 0.025 | [7.5, 8.5] | 7.996 | HOLD |
| Q3 D4 bounded core: deep <= rest, and 0.01a shell <= 0.1a shell, at all three a | inequalities | deep/rest = 0.027, 0.013, 0.010; shells decrease inward at all three a | HOLD |
| Q4 D4 sup at a = 0.05 | <= 2e-6 | 6.08e-7 | HOLD |
| Q5 SQ null: log2 ratio; SQ/D4 sup at 0.025 | [5.5, 6.5]; >= 10 | 6.000; 97.7 | HOLD |
ROUTE ALIVE by the sealed rule (Q1-Q4 hold, controls pass).
Measured, for the record: sup|F| ~ 1.56e4 a^8 on D4 (sampled), delta ~ 874 a^8 (D4), 328 a^8 (SQ).
Seals 8a8380e6ef and ad344efc9b are VOID: the BAR was re-sealed twice before any run to add required seal fields; no prediction or threshold changed between them.
Not shown: linearised inverse bound; proved supremum; any value of a0. Q4 was informed by K1 data (disclosed in the BAR).
