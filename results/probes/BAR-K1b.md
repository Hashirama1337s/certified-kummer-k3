# BAR-K1b: corrected second-order Kummer approximant (probe, follows BAR-K1 FAILED)

Sealed before any run of k1b.py. Written 2026-11-13 by the lead AFTER seeing the K1 results (result_k1.json, VERDICT-K1.md). This is a new approximant and a new test; it does not rescue K1, which stays FAILED.

instruments: k1b.py, k1b_mp.py, k1.py

## Change from K1
One change to the object: mu = 1 + a^4 pi^2/(4V) - delta, with delta = sum over the other bubbles of R' + |l|^2 R''/2 (the isotropic |z|^2 coefficient their far fields add at a bubble centre; order a^8). Reason: |z|^2 is not smooth on the resolved bubble, so any leftover multiple of it gives a residual growing like a^2/s at the bubble centre. One change to the evaluation: F is computed relative to the exact nearest bubble (log1p form), so the deep core can be sampled without rounding loss.

## Sample
Seed 7. "deep" = shells at 0.01a, 0.03a, 0.1a (200 directions each). "rest" = 3000 uniform cell points + shells at 0.2a, 0.5a, a, 2a, 4a, 0.05, 0.1, 0.15, 0.2, 0.25. sup = max over the sample (not a proved supremum).

## What I already know (disclosure)
From K1: away from the deep core F is about 3e-7 at a = 0.05 on D4 (three C4 points), core-shell exponent 7.99 (D4) and 6.01 (SQ). Q4 below is therefore informed, not blind. Q1-Q3 are the new content: K1 never measured the deep core with a bounded approximant.

## Predictions
- Q1 (order, D4): r(a) = sup|F|/a^8 at a = 0.1, 0.05, 0.025: max r / min r <= 2.
- Q2 (order, D4): log2( sup(0.05)/sup(0.025) ) in [7.5, 8.5].
- Q3 (bounded core, D4): at each of the three a, sup over "deep" <= sup over "rest"; and the 0.01a-shell sup <= the 0.1a-shell sup (residual decreases towards the centre).
- Q4 (size, D4): sup|F| at a = 0.05 <= 2e-6.
- Q5 (null, SQ): SQ global log2( sup(0.05)/sup(0.025) ) in [5.5, 6.5], and SQ sup at 0.025 >= 10 x D4 sup at 0.025.

## Controls (any failure = VOID)
- CK: with delta switched off, k1b.py reproduces the four 50-digit K1 values in c4_k1.json to relative 1e-3 (tests the new evaluation path against an independent implementation of the old object).
- C2: Ewald split eta = 6 vs 9, max |dF| <= 1e-3 x sup + 1e-13 (sample includes deep points).
- C4: k1b_mp.py (60-digit potential, finite-difference Hessian) at 4 points incl. 0.01a and 0.03a shells: |F_mp - F_numpy| <= 1e-3 |F_mp| + 1e-13.
- No NaN.
Noise floor: double rounding in the log1p form is about 1e-13 absolute; the smallest sealed quantity (D4 sup at 0.025) is expected near 1e-9 to 1e-8, well above it. Baseline to beat: K1's approximant (unbounded at the centre; sampled sup 3.4e-6 at a = 0.05).

## Verdict rule
ROUTE ALIVE if Q1-Q4 hold and controls pass. Otherwise FAILED, and no third version of this approximant is tried without the owner's say. Q5 failing alone is reported as FAILED for the symmetry explanation.
Not shown by this probe: the linearised inverse bound; a proved (not sampled) supremum; anything about a0.

## Seal fields
- formula_ref: k1b.py:1-9, k1b.py:36-41, k1.py:1-9, k1.py:52-62
- calibrated_on: 08b7d89f85
- depends_on: 08b7d89f85 (control CK compares against c4_k1.json)
- null_not_symmetry: the SQ lattice has a different point group (order-4 axes only), which admits a degree-4 invariant harmonic polynomial of type (2,2); D4 (point group W(F4)) admits none below degree 6. The null changes the predicted exponent from 8 to 6.
- target_form: Q1 max r/min r in [1, 2]; Q2 in [7.5, 8.5]; Q3 inequalities as written (<=); Q4 sup in [0, 2e-6]; Q5 in [5.5, 6.5] and ratio >= 10. All compared on unrounded doubles.
- parameters: a = 0.1, 0.05, 0.025; lattice = D4, SQ; Ewald eta = 6, 9 (control C2)
- files: BAR-K1b.md, k1b.py, k1b_mp.py, k1.py
- command: memcap 8 python k1b.py result_k1b.json ; memcap 4 python k1b_mp.py c4_k1b.json
