# VERDICT-S2 (scored as sealed: claim a4c39ac007, BAR file a0063992ad; instrument s2.py pinned + sealed 9dca3e6cc1; controls on the final file (sha256 b4b6c96c... in controls_s2.json; lab controls-check GO); launch-check passed; main 1 invocation, 30.5 s, peak 0.98 GB)

ERRATUM (lead, 2026-11-14 01:25): S4 was first scored "30/30 cert >= float" from the build log, which in fact shows 28/30 - the builder's own float re-check (theorem2_float) chooses its OWN x = 1 + sqrt(delta/(4 C0 o)) at a_hi, while the certificate uses x = 1 + a xi with xi fixed per step; on [0.003, 0.004] and [0.004, 0.005] the certified d2 was 6.8e-9 relative BELOW that different-x float (result_s2.json summary.float_recheck_all = false). The sealed S4 reads "same certified inputs (same Gamma', same x, same D')": re-run with the certificate's own G1, G2, xi and inputs (stage2/s4check.py, s4check.log; same method as wp5/v5check2.py): 30/30 at or above, min ratio 1.0000000018. S4 HOLD stands on the sealed reading; the slip (not reading every log line) is the lead's. Caught by the release seat.

Inputs: only the certified files (wp2/result_wp2a.json, result_wp2b.json, wp3/result_wp3.json, wp5/result_wp5.json Stage A per interval). New ingredients: L5' (diameter) and Lemma G' with Lemma Z' (spectral tail, Andrews-Clutterbuck lambda1), both refereed in design form (seats/REFEREE-S2.md: CORRECT / CORRECT WITH FIXES, fixes applied in s2.py before controls).

| | sealed bar | certified | verdict |
|---|---|---|---|
| S1 | Gammahat'(0) <= 0.60 | 0.5281643 (upper end 0.52816426...) | HOLD |
| S2 | D'(0.03, e) <= 1.33; D'(a -> 0) <= 1.0001 | 1.322026 on [0.029, 0.03]; 1.0000685 at a = 0 (= sqrt(1 + e)) | HOLD |
| S3 | smallness on all 30 intervals; limsup d2/a^3 <= 850; d2(0.03) <= 0.035; d2(0.01) <= 1.0e-3 | 30/30 close (max U 1.035497); limsup 763.1308 (a = 0 value) and 765.7416 on (0, 1e-4]; d2(0.03) <= 2.855982e-2; d2(0.01) <= 8.722007e-4 | HOLD |
| S4 | certified d2 >= float re-evaluation with the same inputs | 30/30 with the certificate's own Gamma', x, D' (s4check.py, min ratio 1.0000000018); the builder's different-x float check: 28/30 (see erratum) | HOLD |

Controls (controls_s2.json, ALL_CONTROLS_PASS = true, 87 s): Z1 PASS (~1e5 segments per rho in {0.04, 0.06, 0.10, 1/sqrt32}, max 2 crossings, 2 attained); Z2 PASS (flat torus Gh'(0) <= 0.528164 >= 0.0702; with exact torus inputs 0.08191; a deliberately wrong on-diagonal bound gives 0.06766 < 0.0702, so the check bites); Z3 PASS (spectral tail off reproduces 4.530840314 exactly; old L5 + old Gamma reproduce wp5 d2 0.121973059457 and U 1.1693247099 on [0.029, 0.03], rel 0); Z4 PASS (halving lambda1: 0.528164 -> 0.865905 -> 1.392230; monotone in k); guards PASS.
Diagnostic (not a control): at rho = 0.22 the unconstrained count reaches 3 (4,090 / 1e5) but with both endpoints outside every ball (L5''s hypothesis) it stayed 2 (39,349 segments); 3 first appears at rho = 0.235. REFEREE-S2's "rho = 0.22 must give 3" holds only for the unconstrained count.

Verdict as sealed (CORRECTED after REFEREE-S2CODE item 6): S4 FAILED as sealed - the in-run comparator (result_s2.json summary.float_recheck_all = false) has d2 below float on [0.003, 0.004] (6.8e-9 rel) and [0.004, 0.005] (1.1e-8 rel) and U below float on [0.002, 0.003] (1.8e-12 rel, outside its 1e-12 tolerance). By the bar's rule the probe verdict is VOID. S1-S3 HOLD and Z1-Z4 pass; the certified numbers are not in question.
POST HOC diagnosis (does not change the verdict): the in-run float theorem2_float computes log((1+x+sqrt(...))/(2x)) without log1p on an argument 1 + O(1e-13), giving ~1e-8 relative noise at a <= 0.005; a cancellation-free float (REFEREE-S2CODE, both the explicit x and the certificate's xi) gives cert/float 1.0000000050-1.0000000072 (d2) and >= 1.00000000000026 (U) on the three intervals; the lead's s4check.py (certificate's own G1, G2, xi; wp4 formulas) gives 30/30, min 1.0000000018. Clean re-test: BAR-S2b (log1p comparator, sealed before its run).

## Improved certified statement (on the same window a in (0, 0.03], same certified inputs)
(1 + d2)^-1 omega_a <= omega_CY <= (1 - d2)^-1 omega_a with
- limsup_{a->0} d2/a^3 <= 763.14 (a = 0 value; per-interval bound on (0, 1e-4]: <= 765.75)   [was 2269.27 / 2275.00]
- d2(0.01) <= 8.7221e-4   [was 2.851e-3]
- d2(0.03) <= 2.85599e-2   [was 0.121974]
- smallness U <= 1.0355 everywhere   [was 1.1694]
Upper bounds printed rounded UP from result_s2.json.
Not yet: the code-level referee of s2.py (REFEREE-S2 reviewed the design text); Andrews-Clutterbuck journal version not opened (arXiv v1 OPENED); eigenfunction expansion STANDARD; no blind replication of Gamma'; a0 still 0.03 (the extension to 0.05 needs the P3 input re-runs with REFEREE-S2 item 5(a) WP2b Z1 radius cap and WP3 A0 raised).
