# VERDICT-S2b (scored as sealed: claim 1e719251da, files BAR-S2b.md + s2b.py sealed 70c843f630, s2b.py + s2.py pinned; launch-check passed; GO-S2b; run 15.6 s)

Re-test of the VOID BAR-S2 (a4c39ac007; in-run comparator cancellation, REFEREE-S2CODE item 6) with the comparator fixed in code BEFORE the run. S1-S3 are a regression check (values seen in S2); S4 is the test.

| | sealed bar | result | verdict |
|---|---|---|---|
| S1 | Gammahat'(0) <= 0.60 | 0.5281643 | HOLD |
| S2 | D'(0.03) <= 1.33; D'(a -> 0) <= 1.0001 | 1.322026; 1.0000685 (a = 0 value, not a bound at any a > 0) | HOLD |
| S3 | smallness 30/30; limsup <= 850; d2(0.03) <= 0.035; d2(0.01) <= 1.0e-3 | 30/30 (max U 1.035497); limsup 763.1308 (a = 0 value; certified on (0, 1e-4]: 765.7416); 2.855982e-2; 8.722007e-4 | HOLD |
| S4 | certified >= cancellation-free float on all 30 intervals | 30/30 (run_s2b.log: 30 x "cert>=float True"; summary.float_recheck_all = true) | HOLD |
| Z5 | certified per-interval values identical to result_s2.json | 214/214 d2, U, d2/a^3 fields identical | PASS |

Verdict as sealed: STAGE-2 PROBE CERTIFIED.
Certified improvement on (0, 0.03] with the same certified inputs: d2/a^3 <= 765.75 on (0, 1e-4] (limsup value 763.14 at a = 0), d2(0.01) <= 8.7221e-4, d2(0.03) <= 2.85599e-2, U <= 1.0355 (was 2275.00 / 2269.27, 2.851e-3, 0.121974, 1.1694).
Code refereed (seats/REFEREE-S2CODE.md: L5', Gamma', wiring CORRECT; independent Gamma' and d2 agree on the safe side); design refereed (seats/REFEREE-S2.md). Remaining: Andrews-Clutterbuck journal version not opened (arXiv v1 opened); eigenfunction expansion STANDARD; no blind replication of Gamma' (referee G's independent re-implementation is the nearest).
