# VERDICT-WP2b (scored as sealed: claim cdf666b973, BAR file 817e9a43b8, instrument pinned + sealed 6392aabef0; launch-check passed; GO 2026-11-14T00:12:16Z; main run 2 invocations, ~10.5 min, peak 0.70 GB)

Main resolution: tol 0.1, nT 256, nU = nW = 32; bulk cover 62,310 final boxes (70,685 evaluated), Ewald eta 16 (tails <= 4.9e-6); six a-sub-intervals of width 0.005 covering [0, 0.03].

| a | e | a^2 ||dR|| | a^2 B+ = a^2 C0 (upper) |
|---|---|---|---|
| [0, 0.005] | 6.43e-8 | 4.41e-8 | 2.0000000000004 |
| [0.005, 0.01] | 1.03e-6 | 1.64e-6 | 2.00000000046 |
| [0.01, 0.015] | 5.21e-6 | 1.47e-5 | 2.000000028 |
| [0.015, 0.02] | 1.65e-5 | 1.23e-4 | 2.00000053 |
| [0.02, 0.025] | 4.02e-5 | 7.46e-4 | 2.0000053 |
| [0.025, 0.03] | 8.34e-5 | 3.27e-3 | 2.0000353 |

Lower bound at a = 0.03 (point (u,w) = (0,0) on E): a^2 B+ >= 1.99999864, a^2 C0 >= 1.99999863.
Bulk |ddbar psi_1|_op: certified sup <= 102.86, box-centre lower value 93.93 (the doc's 95.9 is a sample, consistent).

| | sealed bar | measured | verdict |
|---|---|---|---|
| B1 | e <= 1e-3 | max 8.34e-5 | HOLD |
| B2 | a^2 ||dR|| <= 1e-2 | max 3.27e-3 (set by the crude Z2 bound) | HOLD |
| B3 | a^2 B+, a^2 C0 <= 2.05 | max 2.0000353 | HOLD |
| B4 | a^2 B+ lower bound at a = 0.03 in [1.99, 2.05] | 1.99999864 | HOLD |
| B5 | >= 3 sub-intervals covering (0, 0.03], one touching 0, each passes B1-B3 | 6 sub-intervals, all pass | HOLD |

Controls, re-run on the SEALED file (the build seat's controls_wp2b.json predates its last edit at 20:11, so it does not count): controls_wp2b_sealed.json/.log, ALL_CONTROLS_PASS = True, 182 s, peak 1.73 GB.
- K1 lattice off: e = 0 in all zones; a^2 B+ and a^2 C0 in [1.9999999999999991, 2.0000000000000018], containing 2. PASS.
- K2: 4000 double-precision samples (4 regions x 500 x a in {0.01, 0.03}) all inside the certified bounds (control resolution); also inside the main bounds (sampled max e 7.2e-5 vs main 8.34e-5; sampled max a^2 bisec 1.99999999999 vs 2.0000353). PASS.
- K3 lattice x10: e and ||dR|| ratios 10.000 (>= 5). PASS. Caveat from the builder: the E-bounds are linear in the lattice term by construction, so K3 tests wiring rather than sensitivity.
- NaN/inf guards: PASS.

Verdict as sealed: CURVATURE INPUTS CERTIFIED (all of B1-B5 HOLD, controls pass).

Carried to the analysis chain: for every a in (0, 0.03]: e(a) <= 8.34e-5, a^2 B+ = a^2 C0 <= 2.0000353 (mirror lemma, (1-e)^-2 kept).
Not refereed yet: the mirror lemma, Lemma D (||dR|| bound), Lemma M2 (BUILD-WP2b.md sections 1-2).
Weaknesses (BUILD-WP2b.md): the float layer relies on IEEE round-to-nearest plus Higham's dot-product bound; the outward rounding is swapped in memory into the imported wp2a module (the wp2a.py file is unchanged); the Weyl chamber and the 6 complex structures are checked by enumeration, not proved; bulk boxes cover a slight superset of V0 (conservative); Z2 dR bound is ~1e4 x the sampled value (it does not affect B+).

ERRATUM 2 (lead, 2026-11-13 23:45, from the v4 editor pass CHANGES-v4 Part C): further printed values were rounded down; the authoritative certified values are the result JSON fields, rounded UP: C_out 17226.21 (JSON 17226.2013), R_out [0.01,0.02] 63239.14 (63239.138), C_neck [0.02,0.03] 15100.04 (15100.039), a^2 C0 rows per result_wp2b.json (e.g. 2.0000000281), Gammahat(0) <= 4.530841 (4.530840314), D(0) <= 1.015342 (1.015341043). BUILD-WP5 3.2 says K_t = 2 kappa in the continuity step; the code uses 4 kappa (conservative; final d2 unaffected). No certified number changed; the runs read full-precision JSON.
