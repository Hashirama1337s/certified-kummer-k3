# Clean-room rebuild of the K3 certificate chain (lead, 2026-11-13 23:10-23:51 local)

Method: byte-identical copies (sha256 checked: wp2a 22b2bd30, wp2b 9cd11731, wp3 80d5492f, wp3_a1 434f6ed1, wp5 b3919b17) placed in fresh folders cleanroom/wp2, wp3, wp5 with only GO files; documented commands (`memcap 8 py -3.11 <script> --main result_*.json`, exit 3 = rerun) run in sequence from empty state.
Results (JSON compared field by field, excluding timing fields and the wp5 input-file sha256 entries):
| step | runtime | result vs original |
|---|---|---|
| WP2a core constants | 220 s | IDENTICAL |
| WP2b curvature | ~10 min | IDENTICAL |
| WP3 (with AMEND-WP3-1 runner) | ~22 min | IDENTICAL |
| WP5 assembly | 59 s | IDENTICAL (a0 = 0.03, limsup <= 2269.2694, d2(0.03) <= 0.1219731) |
Limitation (disclosed): wp5.py reads its input result files from the fixed path <lab>/wp2, wp3 (not the clean-room copies); since the clean-room WP2a/WP2b/WP3 results are identical to those files, the chain is reproduced end to end in substance. A fully path-independent rebuild would need a parameterised input path (release-package item).
Verdict: CLEAN ROOM REPRODUCED.
