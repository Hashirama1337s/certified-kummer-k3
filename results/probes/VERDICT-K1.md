# VERDICT-K1 (scored exactly as sealed, seal 08b7d89f85)
Controls: C1 4.8e-13 PASS; C2 1.8e-15 PASS; C3 1.0e-12 PASS (both lattices); C4 4/4 PASS (worst: argmax point, relative difference 7.9e-4 against a 1e-3 bar); no NaN.
| | sealed bar | measured | verdict |
|---|---|---|---|
| P1 D4 max r / min r (r = sup/a^8) | <= 2 | 1.48 (8.77e4, 8.75e4, 1.29e5) | HOLD |
| P2 D4 global log2 ratio 0.05 -> 0.025 | [7.5, 8.5] | 7.44 | FAILED |
| P3 D4 sup at a = 0.05 | <= 1e-5 | 3.42e-6 | HOLD |
| P4 D4 core-shell log2 ratio | [7.5, 8.5] | 7.99 | HOLD |
| P5 SQ core-shell log2 ratio (null) | [5.5, 6.5] | 6.01 | HOLD |
| P6 SQ global log2 ratio (null) | < 7.5 | 5.98 | HOLD |

By the sealed verdict rule the ROUTE-ALIVE condition (P1-P4 all hold) is NOT met: P2 FAILED.

PROCESS SLIP: the claim and BAR text were sealed before the run (08b7d89f85), but the file-seal step exited with warnings and launch-check printed DO NOT LAUNCH; a shell pipe swallowed the exit code and the run started anyway. BAR-K1.md, k1.py and k1_mp.py were not edited between sealing and the run.

POST HOC (not part of the sealed test): on D4 the sampled sup always sits on the innermost shell (0.1 a). Diagnosis by hand: the other bubbles' far-field a^8 terms add a small multiple delta ~ a^8 of |z|^2 at each bubble; |z|^2 is not a smooth function on the resolved bubble, so the residual delta x a^2/s grows without bound as s -> 0. The sampled sup therefore depends on where the innermost shell is, and the approximant as defined has an unbounded residual at the exceptional curve. Proposed fix (needs its own sealed probe): absorb delta by rescaling the bubble (lambda x Eguchi-Hanson of size b, lambda b^4 = a^4, target constant lambda^2).
