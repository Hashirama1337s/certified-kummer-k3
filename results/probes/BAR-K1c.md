# BAR-K1c: smooth approximant: residual, Ricci and curvature (probe)

Written 2026-11-13 by the lead before k1c.py exists. Predictions come from hand estimates and from the settled K1b numbers; no K1c quantity has been computed.

instruments: k1c.py

## Object and sample
Exactly SPEC-K1c.md sections 1, 2 and 5 (lambda = 0.1). sup = max over the sample, not a proved supremum.

## Reasoning
Outside r2 nothing changes from K1b (sup|F| about 1.56e4 a^8). In the annulus the cut-off regular part is about a^4 c6 |z|^6 with c6 about 1.5e4/2; its d d-bar measured against Eguchi-Hanson at |z| about lambda a costs about c6 lambda^2 a^8 x O(10), i.e. at most a few 1e3 a^8. Ricci = d d-bar F: core terms vary on scale a (so about C_F a^8 / a^2 = 1e4-1e5 a^6), annulus terms vary on scale lambda a (so up to about 2e3 a^8/(lambda a)^2 = 2e5 a^6). Curvature: Eguchi-Hanson's own (max bisectional 2/a^2 at the exceptional curve, 1.99/a^2 at 0.1a) plus O(a^6). Without the cut-off the degree-6 odd-odd terms make Ricci grow like 1/s towards the centre.

## Predictions
- R1 (residual): D4 r(a) = sup|F|/a^8 <= 2.5e4 at each of a = 0.1, 0.05, 0.025, and max r / min r <= 2.
- R2 (Ricci size): D4 RIC <= 1e6 a^6 at each a.
- R3 (Ricci order): D4 log2( RIC(0.05) / RIC(0.025) ) in [5.5, 6.5].
- R4 (curvature): D4 BIS x a^2 in [1.90, 2.05] at each a.
- R5 (positivity): D4 POS >= 0.9 at each a.
- R6 (blow-up control, uncut object, D4, a = 0.05): RIC on the 0.01a shell >= 10 x RIC on the 0.1a shell.
- R7 (null, SQ): SQ sup|F| at a = 0.025 >= 10 x D4 sup|F| at a = 0.025.

## Controls (any failure = VOID): K-EH, K-FS, K-JET, K-MATCH, K-K1b as in SPEC-K1c.md section 4, read from k1c_controls.json; polynomial-model truncation error at r2 <= 1e-6 relative; no NaN.
Noise floor: Ricci is computed from jets of the log1p form; expected absolute noise about 1e-12/(lambda a)^2, i.e. below 1e-6 x the predicted RIC at a = 0.025. Baseline: K1b's uncut approximant (R6).

## Verdict rule
SHAPE FIXED if R1-R5 hold and controls pass. R6 failing means the referee's blow-up diagnosis is not confirmed by this instrument (report FAILED for R6; R1-R5 still decide the shape). If R2/R3 fail, the cut-off design is FAILED and goes back to the lead; no silent retuning of lambda.
Not shown: proved suprema; uniformity in a beyond three values; anything about a0.

## Seal fields
- formula_ref: SPEC-K1c.md:3-9, SPEC-K1c.md:12-16
- calibrated_on: ff0c52846b
- depends_on: ff0c52846b
- null_not_symmetry: SQ has a different point group, which admits a degree-4 invariant harmonic polynomial of type (2,2); this changes the residual order from 8 to 6 (settled in ff0c52846b). The uncut object (R6) differs from the cut one exactly by the terms the cut-off removes.
- target_form: R1 r in [0, 2.5e4] and ratio in [1, 2]; R2 RIC/a^6 in [0, 1e6]; R3 in [5.5, 6.5]; R4 in [1.90, 2.05]; R5 >= 0.9; R6 ratio >= 10; R7 ratio >= 10. Unrounded doubles.
- parameters: a = 0.1, 0.05, 0.025; lattice = D4, SQ; object = cut, uncut
- files: BAR-K1c.md, SPEC-K1c.md, k1c.py
- command: memcap 8 python k1c.py --main result_k1c.json
