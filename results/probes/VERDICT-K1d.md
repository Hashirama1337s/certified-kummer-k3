# VERDICT-K1d (scored as sealed: BAR-K1d 4143022d85 + control amendment 172cbd900f sealed before the run; instrument pinned 174b6e7f29; launch-check passed)
Controls (k1d_controls.json, build seat, before the pin): K-Q 1.8e-14 PASS; K-SMOOTH ratio 1.00013 PASS (eta == 0 contrast 99.98x); K-MODEL g 3.1e-16, F 3.0e-16 PASS; K-ZERO F 7.9e-18 PASS, RIC band points 3.9e-6 <= 1e-5 (amended) PASS, other points <= 1.1e-13 PASS. No NaN.
| | lambda = 0.5 | lambda = 1 | bar | verdict |
|---|---|---|---|---|
| D1 Ricci order, log2 (0.05->0.025; 0.025->0.0125) | 6.0000; 6.0001 | 5.9998; 6.0000 | [5.5, 6.5] | HOLD both |
| D2 Ricci/a^6 at a = 0.025 | 1.398e6 (31x below K1c) | 3.491e6 (13x below K1c) | <= 4.4e6 | HOLD both |
| D3 Ricci 0.01a shell / 0.1a shell | 1.000 at every a | 1.000 at every a | <= 2 | HOLD both |
| D4 sup|F|/a^8 | 2.00e4 (all a) | 1.98e5 (all a) | <= 1e6 | HOLD both |
| D5 max bisectional x a^2 | 2.0000 (all a) | 2.0000 (all a) | [1.90, 2.05] | HOLD both |
| D6 positivity | >= 0.9994 | >= 0.9994 | >= 0.9 | HOLD both |
SHAPE FIXED for both lambda. Best: lambda = 0.5 (residual 2.0e4 a^8, Ricci 1.40e6 a^6, peak Ricci at 0.7a, inside the cut-off annulus).
Observed (not barred): every constant is the same at all three a to 4 digits, i.e. the approximant's defects scale exactly (a^8 residual, a^6 Ricci, 2/a^2 curvature), which is what a uniform-in-a certificate needs.
Not shown: odd-odd terms of degree >= 12 not removed (sampled effect below noise); sampled suprema; three values of a only.
