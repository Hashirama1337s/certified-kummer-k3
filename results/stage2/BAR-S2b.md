# BAR-S2b: clean re-test of the Stage 2 probe with a cancellation-free comparator

Written 2026-11-14 by the lead BEFORE this run. Disclosure: this is a re-test after BAR-S2 (a4c39ac007) was VOID on its in-run comparator (log() cancellation in the double-precision float, REFEREE-S2CODE item 6). The S1-S3 thresholds are unchanged from BAR-S2 (their values were seen in the S2 run, so S1-S3 here are a regression check, not a fresh prediction); S4 is the test: a cancellation-free same-formula comparator, fixed in code before the run.

instruments: s2b.py, s2.py

## Object
As BAR-S2: L5' + Lemma G' re-assembly of Theorem 2 on the 30 WP5 a-intervals covering (0, 0.03] from the certified input files. s2b.py imports the sealed s2.py (9dca3e6cc1) unchanged and replaces only theorem2_float by theorem2_float_cf (log1p/expm1; same explicit x rule as s2.py's float); fresh state folder state_s2b, log run_s2b.log.

## Predictions
- S1: certified Gammahat'(0) <= 0.60.
- S2: certified D'(0.03, e) <= 1.33 and D'(a -> 0) <= 1.0001.
- S3: smallness on all 30 intervals; certified limsup d2/a^3 <= 850; d2(0.03) <= 0.035; d2(0.01) <= 1.0e-3.
- S4: on all 30 intervals the certified d2 and U are >= the cancellation-free float re-evaluation (s2b's theorem2_float_cf) times (1 - 1e-12) - i.e. result_s2b.json summary.float_recheck_all = true.

## Controls (any failure = VOID)
- Z5: the certified per-interval values (d2, U) reproduce result_s2.json bit for bit (same sealed s2.py engine, same inputs) - the runner changes only the comparator.
- Controls of s2.py (Z1-Z4) stand from controls_s2.json (same sealed s2.py; not re-run).
- No NaN/inf.

## Verdict rule
STAGE-2 PROBE CERTIFIED if S1-S4 hold and Z5 passes. Any interval with certified below the cancellation-free float = VOID and a bug hunt.
Noise floor/ceiling: the float evaluation is the floor of every certified upper bound; the comparator's own rounding noise is now O(1e-15) relative (log1p/expm1), well inside the 1e-12 tolerance.

## Seal fields
- formula_ref: STAGE2-DESIGN.md:1-154, ../ANALYSIS-CHAIN-v2.2.md:140
- calibrated_on: a4c39ac007
- depends_on: a4c39ac007, 055009634d
- null_not_symmetry: Z5 is a known-answer reproduction of the S2 engine output; S4's comparator is an independent float formula of the same theorem.
- target_form: S1 [0, 0.60]; S2 [0, 1.33] and [0, 1.0001]; S3 boolean, [0, 850], [0, 0.035], [0, 1.0e-3]; S4 boolean.
- tolerances: unchanged from BAR-S2 (time step 1e-3, justified in BUILD-S2.md); comparator tolerance 1e-12 relative.
- parameters: the 30 WP5 a-intervals of width 0.001 covering (0, 0.03]; fixed - re-test of a sealed engine on frozen certified inputs
- files: BAR-S2b.md, s2b.py
- command: memcap 8 python s2b.py --main result_s2b.json
