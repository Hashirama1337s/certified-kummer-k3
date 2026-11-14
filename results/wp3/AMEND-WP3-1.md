# AMEND-WP3-1: bulk refinement tolerances (instrument efficiency only)

Written by the lead 2026-11-13 after stopping the first main run at bulk iteration 700, BEFORE any WP3 result exists
(no bulk leaf bound had been reported; result_wp3.json does not exist). Parent claim 4158aba4f8; wp3.py (cb86c3f440) unchanged.

Why: with BULK_TOLF = 1e-3 (a box finishes only if its F bound is within 0.1% of the running max) and RCAP_FRAC = 0.8
(bulk Ricci boxes must sit below 0.8 x the neck's leading value), the waiting list grew without bound
(16k -> 1.25M boxes over 700 iterations, 95k evaluated). Neither tolerance is a prediction: the certified outputs
are the max over the FINAL boxes whatever the tolerance, and the sealed targets only need C_out <= 18000 and R_out <= 1e6.

Change (wp3_a1.py imports wp3.py and sets two module constants, nothing else):
- BULK_TOLF 1e-3 -> 3e-2 (the value wp3.py's own controls use, line 2738). Effect: C_out may be reported up to 3% above
  the certified point lower value (about 16775 -> at most about 17280), instead of 0.1%.
- RCAP_FRAC 0.8 -> 1.0. Effect: bulk boxes stop once their Ricci bound is below the neck's leading value L, which is
  itself a lower bound for R_out, so R_out is not loosened.
- The partial bulk state from the stopped run is moved to state_wp3/superseded/ (possibly mid-write when stopped);
  the bulk stage restarts from scratch. Lead and S2 neck stages (complete, same sealed code, same tolerances) are reused.
Predictions, controls and verdict rule of BAR-WP3 are unchanged.
