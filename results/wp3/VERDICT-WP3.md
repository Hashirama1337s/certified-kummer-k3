# VERDICT-WP3 (scored as sealed: claim 4158aba4f8, BAR file 61be266829, instrument pinned + sealed cb86c3f440, amendment AMEND-WP3-1 + wp3_a1.py sealed eddd0d7e7b BEFORE data; launch-check passed; GO 2026-11-14T01:37:07Z)

Run history: first main run stopped by the lead at bulk iteration 700, before any bulk result. The queue grew from 16k to 1.25M boxes with only 95k evaluated, because the bulk stopping tolerance (1e-3) was far tighter than the 7% headroom the target allows. AMEND-WP3-1 changed only BULK_TOLF (to 3e-2, the value the file's own controls use) and RCAP_FRAC (to 1.0), and restarted the bulk from scratch. Lead and S2 neck stages were reused, same sealed code. The amended bulk completed in 207 iterations (25,477 boxes evaluated, 22,601 leaves), 2 invocations, peak 3.23 GB.
Controls ran on the sealed file: controls_wp3.json (21:34) is newer than the last wp3.py edit (21:25). ALL_CONTROLS_PASS = true.

Certified (all a in (0, 0.03], outside the cores):
| a | C_out (|F| / a^8) | C_neck | R_out (|Ric| / a^6) | R_neck terms (L + a^2 S2 + sqrt(u) rest) | R_bulk |
|---|---|---|---|---|---|
| (0, 0.01] | 17226.20 | 15130.9 | 62359.3 | 61669.6 + 228.9 + 460.8 | 61663.1 |
| [0.01, 0.02] | 17226.20 | 15130.9 | 63239.1 | 61669.6 + 915.7 + 653.8 | 61663.1 |
| [0.02, 0.03] | 17226.20 | 15100.0 | 65064.3 | 61669.6 + 2060.3 + 1334.2 | 61663.1 |
Lower values: sup |F8bulk| >= 16918.28 (point (1/4, 0, 0, 0); Ewald eta 6 and 8 enclosures overlap); running box lower value 16761.13; R_out as a -> 0 >= 61663.76 (lead L interval [61663.76, 61669.56]).

| | sealed bar | measured | verdict |
|---|---|---|---|
| U1 | C_out <= 18000 | 17226.20 | HOLD |
| U2 | certified lower bound for sup |F8bulk| >= 16283 | 16918.28 | HOLD |
| U3 | R_out <= 1.0e6 | 65064.3 | HOLD |
| U4 | >= 3 a-sub-intervals, C_out spread <= 5% | 3 pieces, spread 0 (see note) | HOLD |
| U5 | no sampled value above a certified bound | samples 16283, 9040/12391/13210 (shells), Ricci 4.9e4 all below | HOLD |
Note on U4: the C_out bound is set by the bulk, whose bound is computed once with eps at its upper end (a = 0.03). So the three pieces are identical by construction, not by a measured flatness. The neck values (15131 / 15131 / 15100) do vary.

Controls: J1 PASS (double-precision F8bulk sup 14569.13068 vs 14569.136078, relative 3.7e-7; point enclosures contain all 396 values). J2 PASS (100/100 overlaps at eta 6 and 8). J3 PASS (dropping the cancellation partner blows up to 2.9e7 at |z| = 0.1; dropping N[psi1]/s^3 breaks 64/64). J3b (extra) PASS. J4 NOT DONE as a certified integral (diagnostic only per BAR; reason in BUILD-WP3.md). NaN/inf guards 8/8.

Verdict as sealed: OUTER INPUTS CERTIFIED (U1-U5 HOLD; J1-J3 pass; J4 diagnostic not run).

Consequences for the chain: the global residual constant is max(core 14123.85, outer 17226.20) = 17226.20 a^8, set by the bulk at the Voronoi facet centre. The doc's sampled 1.63e4 was low; the true sup >= 16918. The global Ricci constant is set by the core (1.696e6 a^6 vs outer 6.51e4 a^6).
Open items: e_upper 1.37e-4 here vs WP2b e <= 8.34e-5. They are likely different normalisations (to reconcile in WP5 before use). Hand-derived identities are checked numerically, not symbolically. A few libm ** calls are used with >= 1e-9 slack, and the near-point sums use a first-order running error bound (BUILD-WP3.md weaknesses). DEFINITIONS v3 must absorb BUILD-WP3's 5 doc problems.

ERRATUM 2 (lead, 2026-11-13 23:45, from the v4 editor pass CHANGES-v4 Part C): further printed values were rounded down; the authoritative certified values are the result JSON fields, rounded UP: C_out 17226.21 (JSON 17226.2013), R_out [0.01,0.02] 63239.14 (63239.138), C_neck [0.02,0.03] 15100.04 (15100.039), a^2 C0 rows per result_wp2b.json (e.g. 2.0000000281), Gammahat(0) <= 4.530841 (4.530840314), D(0) <= 1.015342 (1.015341043). BUILD-WP5 3.2 says K_t = 2 kappa in the continuity step; the code uses 4 kappa (conservative; final d2 unaffected). No certified number changed; the runs read full-precision JSON.
