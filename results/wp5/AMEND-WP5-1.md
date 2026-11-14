# AMEND-WP5-1: control X2 restated one-sided (before any WP5 code or data)

Written by the lead 2026-11-13, minutes after sealing BAR-WP5 (488ee8329c); wp5.py does not exist yet. Parent claim 055009634d.
Why: X2 as sealed asks the interval enclosure of Gammahat(0) to "contain 4.5308 (to its printed precision) or be below it". 4.5308 is a rounded double-precision grid value (lint #98-99); a correct enclosure slightly above it would VOID the run, the WP2a failure mode.
Change: X2 now reads: "the certified upper bound for Gammahat(0) lies in [0.0702, 4.5871] (between the sharp flat-torus value 0.0702 of ANALYSIS-CHAIN 3.6 and the closed form (27e/16) = 4.5871 of 3.2); the BUILD file reports its distance from the double-precision 4.5308 (diagnostic)."
Thresholds: no threshold of V1-V5 changes; only X2's acceptance window changes as stated. What counts as FAILED / VOID is otherwise unchanged (an X2 value outside [0.0702, 4.5871] = VOID).
