# A computer-certified closeness bound for the Kummer K3 surface of the Hurwitz lattice

**Moki & Julio**

Version 1.0.1 (2026-11-17). Every number below is traced to a result file of the accompanying code in
`CLAIMS.md`; certified upper bounds are printed rounded up, lower bounds rounded down.

---

## Abstract

Let X be the Kummer K3 surface of the torus C^2/Λ, where Λ is the Hurwitz lattice (the D4 lattice of Hurwitz
quaternion integers), and let ω_a be an explicit Kähler form on X obtained by gluing 16 Eguchi-Hanson bubbles of size a
into the flat orbifold, with a lattice Green-function correction and a C-infinity cut-off of the odd-odd part of the
regular potential on the window [a/4, a]. Let ω_CY be the Ricci-flat Kähler form in the same class. We prove, with
interval-arithmetic computer assistance, that for every a in (0, 0.03]

  (1 + d(a))^-1 ω_a ≤ ω_CY ≤ (1 − d(a))^-1 ω_a,

with certified bounds d(a)/a^3 ≤ 2275.00 for a in (0, 1e-4], d(0.01) ≤ 2.851e-3 and d(0.03) ≤ 0.121974;
the a = 0 value of a continuous upper bound for d(a)/a^3 is 2269.27, which bounds limsup d(a)/a^3. The proof follows the
continuity-method route of Kobayashi and Lye with every constant made explicit: a Li-Yau/Doeblin heat-kernel bound
for the C^0 step, Yau's second-order inequality for the C^2 step, and certified inputs for the residual, the Ricci
curvature, the bisectional curvature, the diameter and the positivity of ω_a. The value a0 = 0.03 is where the
certificate stops (the proved range of our diameter lemma); it is not a limit of the method. The certificate was
reproduced end to end from empty state, its three load-bearing pieces were replicated by blind independent code, and
the written chain went through nine hostile referee reports. Limitations are listed in Section 7.

---

> **In plain words**
>
> A K3 surface is a four-dimensional shape that carries a perfectly balanced geometry, a Ricci-flat metric. Yau proved
> in 1978 that this geometry exists, but nobody can write it down by a formula.
>
> The Kummer construction gets close to it by hand. Start with a flat four-dimensional torus built on a very symmetric
> grid (the Hurwitz lattice), fold it in half, and smooth each of the 16 sharp corner points with a tiny, known
> bubble of size a. Since 1990 it has been known that the perfect geometry is close to this glued one "if a is small
> enough". How small, and how close, was never computed.
>
> We computed both, and proved them with a computer that keeps a guaranteed error bar on every number. For every
> bubble size a up to 0.03, the perfect geometry stretches or shrinks the glued one by at most a factor 1/(1 − d), with
> d at most 2275.00 a^3 when a is below 1e-4, at most 0.002851 at a = 0.01 (under 0.3 percent), and at most 0.121974
> at a = 0.03.
>
> A picture: a tailor's suit and a perfect fit. Earlier work proved the suit fits well when the bubbles are small. We
> give the measurement, with a guarantee, for every bubble size up to 0.03. Where the picture breaks: the "fit" is
> checked at every point and in every direction at once, and the guarantee comes from a proof, not from measuring.
>
> The method and the qualitative statement are not new (Kobayashi 1990, Lye 2023); what is new is the explicit,
> certified size of the bound on a usable interval.

---

## 1. Result

**Theorem 1 (computer-assisted).** Let ω_a be the approximant of Section 2 (Kummer surface of the Hurwitz lattice,
16 Eguchi-Hanson bubbles of size a, C-infinity cut-off on the window [a/4, a]) and ω_CY the Ricci-flat Kähler form in
its class. For every a in (0, 0.03] the continuity method of Section 3 closes (smallness U = M2 e^η ≤ 1.1694 < 2 on all
30 a-intervals of width 0.001), and

  (1 + d2(a))^-1 ω_a ≤ ω_CY ≤ (1 − d2(a))^-1 ω_a,

with the following certified upper bounds (every printed number is the certified value rounded up):

- d2(a)/a^3 ≤ 2275.00 for every a in (0, 1e-4].
- d2(a)/a^3 ≤ 2326.6 on (0, 0.001]; ≤ 2851.0 on [0.009, 0.01]; ≤ 4517.6 on [0.029, 0.03].
- d2(0.01) ≤ 2.851e-3; d2(0.03) ≤ 0.121974.
- limsup_{a→0} d2(a)/a^3 ≤ 2269.27. This is the a = 0 value of a continuous upper bound; it does not bound d2(a)/a^3
  at any a > 0 (the exact bound at a = 1e-6 is 2269.2701). Its closed form is 8 D(0) sqrt(C_q Γ̂(0)) with
  D(0) ≤ 1.015342, C_q ≤ 17226.21 and Γ̂(0) ≤ 4.530841.
- a0 = 0.03. The binding limit is the proved range of the diameter lemma (L5), not the smallness condition, which is
  far from failing at a = 0.03 (1.1694 < 2). a0 = 0.03 is where the certificate stops; it is not a property of the
  method.

The theorem is conditional only as listed in Section 7 (floating-point layer assumptions, standard facts whose
sources we did not open, and the review status of individual lemmas).

Leading order: d2 ≈ 2 sqrt(C0 o) ≈ 8 D sqrt(C_q γ) a^3, with C0 = 2/a^2 the curvature input, o the oscillation bound of
Section 3.4, γ = Γ̂(0) the heat-kernel constant of Section 3.3, C_q the residual constant and D the diameter bound.

Certified d2 per a-interval (each bound holds for every a in the interval and is essentially the value at its upper
end). Appendix A lists all 30 intervals.

| a-interval | smallness U (continuity step) | d2 ≤ | d2/a^3 ≤ |
|---|---|---|---|
| [0, 0.001] | 1.000003 | 2.327e-6 | 2326.6 |
| [0.002, 0.003] | 1.000081 | 6.592e-5 | 2441.3 |
| [0.009, 0.01] | 1.003501 | 2.851e-3 | 2851.0 |
| [0.029, 0.03] | 1.1694 | 0.121974 | 4517.6 |

(a = 0.01 also lies in [0.01, 0.011], whose enclosure gives the weaker but valid 3.876e-3.)

**Comparison with the qualitative theorem.** In our units Kobayashi obtains closeness of order a, and a^2 after a
Schauder stage, for an approximant whose residual is of order a^4. Our order a^3 comes entirely from the approximant's
a^8 residual: d ~ sqrt(a^-2 a^8). We do not claim a better rate: applied to the same approximant, the qualitative method
would give the same a^3 at the C^2 stage and better after a Schauder stage. What is new is that the constants and the
interval are explicit and certified. The method is Kobayashi's and Lye's, with explicit constants and a Green/Harnack
C^0 step in place of Moser iteration.

**Sharpened bound (Stage 2; certified and refereed, not blind-replicated).** With the same certified inputs, a second
assembly replaces the diameter lemma L5 by a sharper lemma L5' and the heat-kernel constant by a spectral one built on
the eigenvalue bound of Andrews and Clutterbuck [AC]. On the same window (a0 = 0.03) it certifies
d2(a)/a^3 ≤ 765.75 for every a in (0, 1e-4] (limsup value 763.14 at a = 0), d2(0.01) ≤ 8.7221e-4 and
d2(0.03) ≤ 2.85599e-2. Its design and its code were refereed (L5', the spectral constant and the wiring correct; the
code referee's independent re-implementation agrees on the safe side); it has not been replicated blind. Its first
sealed test, BAR-S2, is VOID as sealed: the in-run floating-point comparator lost accuracy to cancellation in a
logarithm, and two d2 values fell below it by about 1e-8 relative and one U value by about 2e-12; the certified values are not in
question. The clean re-test, BAR-S2b, with the comparator fixed before the run, held on all 30 intervals, with
certified values bit-identical to the first run. Theorem 1 stays the WP5 certificate.

## 2. The construction

Full definitions, statements and proofs are in the companion document `docs/DEFINITIONS-AND-LEMMAS-v4.2.md`
(cited below as "Definitions"); section numbers refer to it.

- **Lattice and surface (Definitions 1.1).** Λ = Hurwitz integers = Z^4 ∪ (Z^4 + (1/2,1/2,1/2,1/2)) in C^2 = R^4,
  L = Λ/2 with covolume 1/32, T = C^2/Λ, and X the minimal resolution of T/{±1}, with 16 exceptional curves over the
  two-torsion points L/Λ. The point group W = Aut(L) = W(F4) has order 1152; its subgroup G_J of elements that commute
  or anticommute with the complex structure has order 192.
- **Green function (Definitions 1.2, Lemma S).** S, the L-periodic Green function of the flat Laplacian in Ewald form,
  mean zero over a cell; S − 1/|x − l|^2 is analytic across each lattice point l.
- **Eguchi-Hanson and the remainder kernel (Definitions 1.3).** With ε = a^4 and f(s) = sqrt(s^2 + ε), the
  Eguchi-Hanson potential φ_a(s) = f + a^2 log s − a^2 log(a^2 + f), s = |z|^2, and its far-field remainder R_a.
- **Regular part (Definitions 1.4).** μ = 1 + 8π^2 ε − δ(a), with δ(a) an explicit alternating lattice series, and
  ψ(z) = (μ − 1)|z|^2 + Σ_{l≠0} R_a(|z − l|^2) − (ε/2) S_reg(z) near a bubble.
- **Global potential and approximant (Definitions 1.5).** Φ^b = μ|z|^2 + Σ_{l∈L} R_a(|z − l|^2) − (ε/2) S. Near
  every lattice point Φ^b = φ_a + ψ + pluriharmonic (Lemma P). The odd-odd part ψ_oo of ψ (its Taylor degrees 2 mod 4)
  is not smooth on the resolution (Definitions L2 (b), (c)); it is removed inside each bubble with the C-infinity step
  st(t) = h(t)/(h(t) + h(1 − t)), h(t) = e^{-1/t}, on the window [a/4, a]:
  Φ_a = Φ^b − Σ_l χ_a(z − l) ψ_oo(z − l), ω_a = i∂∂̄Φ_a.
- **Symmetry and regularity (Lemma Q, L2 (iv)).** ω_a is L-periodic, invariant under G_J about every lattice point,
  descends to X, and is a C-infinity closed real (1,1)-form; where positive (certified below) it is a Kähler metric.
- **Expansion (L3).** ψ = ε ψ_1 + Σ_{k≥2} ε^k ψ_k, with ψ_1 harmonic, no quadratic term, no quartic term in ψ_1, and the
  sextic term −(c6/2) H6 with an explicit harmonic sextic H6 and the lattice constant c6 (certified in WP2a). In the
  scaled core the residual F = log(ω_a^2/ω_flat^2) is of order a^8 and the Ricci form of order a^6 (L4).
- **Diameter (L5).** For a ≤ 0.03 and positivity defect e < 1,
  diam(X, g_R) ≤ D(a, e) = sqrt(2(1 + e)) [4a + (17/16)^{1/4} (1/sqrt 2 + 12(π − 2)a)].

The construction uses the Green-sum structure that Fu, Liu, Xiao and Yau identified for the Kummer surface (Section 8);
what we add is a global closed-form potential with explicit cores, on the Hurwitz lattice, whose constants are
certified.

## 3. The analysis chain (proof outline)

Full text: companion `docs/ANALYSIS-CHAIN-v2.2.md` (cited as "Chain"). Conventions: complex dimension 2, real
dimension 4; g_R = 2 Re g; Δ_R = 2 Δ_c; F = log(ω_a^2/ω_flat^2), so Ric(ω_a) = −i∂∂̄F and ω_CY^2 = e^{−F+c} ω_a^2
(Chain 1).

### 3.1 Inputs (Chain 2)

Each input is certified for every a in (0, 0.03] by a sealed run (Section 4); per-interval values are in
`results/wp5/result_wp5.json`.

- (I1) |F| ≤ C_q a^8 with C_q a^8 < 1; η := −2 log(1 − C_q a^8) bounds osc F. Certified C_q = max(C_core_F, C_out)
  = 17226.21 (set by the bulk, WP3); the full core value is C_core_F ≤ 14416.38 (WP5 Stage A).
- (I2) |Ric(ω_a)| ≤ κ := C_R a^6 and σ := sup |Δ_c F| ≤ 2κ. Certified C_R ≤ 1,731,732.5 (the core sets it, WP5
  Stage A); the certificate uses per-interval values (1695936.40 on [0, 0.001]).
- (I3−) R(X, X̄, Y, Ȳ) ≥ −C0 |X|^2 |Y|^2 for g_a-orthogonal X, Y, from the mirror lemma (Definitions L4 (g), Lemma M).
  C0 is the pointwise bound C0 := sup_p (K(s_p) + D_p)/(1 − e_p)^2, where K(s_p) ≤ 2/a^2 is the Eguchi-Hanson
  curvature at p, D_p the local curvature perturbation and e_p the local positivity defect, enclosed zone by zone; the
  C^2 step uses the curvature bound pointwise, so this is the right input. Certified a^2 C0 ≤ 2.0000353 for every a in
  (0, 0.03] (WP2b; on [0.025, 0.03] the zone values are 2.0000352762 in the core zone Z1, 1.98843008 in Z2 and
  0.70784651 outside the cores). Combining the global e and the global ||dR|| instead would give the cruder
  2.0036. At a = 0.03 the true value is at least 1.99999862 (Lemma M2), so the input is tight.
- (I4) diam(g_R(ω_a)) ≤ D(a, e) from L5 (valid for a ≤ 0.03): D(0) ≤ 1.015342, D(0.03) ≤ 1.77518.
- (I5) Positivity: e ≤ 1.37e-4 < 1 (the conservative maximum of the WP2b, WP3 and Stage A values; WP2b alone gives
  e ≤ 8.34e-5).
- (I6) Regularity: everything is C-infinity (Chain 8).

### 3.2 Path (Chain 2, Lemma 2.1)

h_t = −tF + c_t with c_t fixed by ∫ e^{h_t} ω_a^2 = ∫ ω_a^2; solve (ω_a + i∂∂̄u_t)^2 = e^{h_t} ω_a^2. Lemma 2.1:
|h_t| ≤ η; Ric(ω_t) = (1 − t) Ric(ω_a), so ω_1 = ω_CY (Calabi uniqueness); Δ_c^a h_t ≥ −σ; the Ricci and diameter bounds
transfer to ω_t under a two-sided bound Λ^{-1} ω_a ≤ ω_t ≤ Λ ω_a.

### 3.3 Lemma G: explicit oscillation bound (Chain 3)

For a closed Riemannian 4-manifold with Ric ≥ −K, diam ≤ D and Δ_R v ≥ −c:
sup v − avg v ≤ c Γ(K, D). Proof: the Li-Yau parabolic Harnack inequality (Li-Yau 1986, Thm 2.3, pp. 168-169, with n
the real dimension, p. 170) gives a uniform lower bound V H(x, y, t) ≥ c0(t) for the heat kernel; a Doeblin argument
on the minimum of the heat kernel and the semigroup property give 1 − m(t) ≤ (1 − c0)^k; integrating the heat equation
in time bounds the oscillation (Chain 3.3 (i)-(v)). The volume cancels. Γ(K, D) = D^2 Γ̂(K D^2). Certified
(WP5, arb, N = 4000, s-grid 2001): Γ̂(0) ≤ 4.530841, below the closed form (27e/16) = 4.5871. The certificate uses
Γ = min(D^2 Γ̂_cert, Γ_expl).

### 3.4 C^0 step (Chain 4, Proposition C0)

Under Λ^{-1} ω_a ≤ ω_t ≤ Λ ω_a and e^η < 3, with η1 = 2(1 − e^{−η/2}), η2 = 2(e^{η/2} − 1):
osc u_t ≤ o := [2 η1 Γ(κ, D) + 2 η2 Γ(K_t, sqrt(Λ) D)] / (1 − (e^η − 1)/2), for any K_t ≥ (1 − t) Λ κ.

### 3.5 C^2 step (Chain 5.2, Lemma YA)

Yau's inequality Δ_t log T ≥ (Δ_a h + 4 C0)/T − C0 S (T = tr_{ω_a} ω_t, S = tr_{ω_t} ω_a; only g_a-orthogonal pairs
enter, so only (I3−) is needed), matching Kobayashi eq. (31) and Lye eq. (6.10). The maximum principle on log T − Au
with A = C0(1 + x), using det = e^h, gives T ≤ 2(1 + ε2) with
1 + ε2 = e^η (1 + x + sqrt((x − 1)^2 + x δ))/(2x) · e^{C0(1+x)o}, δ = 4(1 − e^{−η}) + (σ/C0) e^{−η}.

### 3.6 Trace to eigenvalues (Chain 6)

The eigenvalues of ω_t relative to ω_a have sum T and product e^h, so with M2 = (1 + ε2) + sqrt((1 + ε2)^2 − e^{−η}):
e^{−η} M2^{-1} ω_a ≤ ω_t ≤ M2 ω_a, and with d2 = M2 e^η − 1 the symmetric form (1 + d2)^{-1} ω_a ≤ ω_t ≤ (1 − d2)^{-1} ω_a.
These algebraic steps are checked in Lean (Section 6.4).

### 3.7 Continuity closure and the final bound (Chain 7, 9)

Yau's theorem gives smooth u_t for each t; the implicit function theorem gives continuity in t; the set of t where
2^{-1} ω_a ≤ ω_t ≤ 2 ω_a is closed, contains 0, and is open when the smallness condition U = M2 e^η < 2 holds with
Λ = 2, diameter sqrt 2 D and any K_t ≥ 2κ. (The certificate uses K_t = 4κ, which is conservative: Γ is nondecreasing in
K. The final d2 is unaffected because it comes from t = 1, where K_t = 0.) At t = 1 the bound is iterated six times.
Theorem 2 of the Chain states the result; Theorem 1 above is its certified instance on the 30 a-intervals. At five
points per interval including the upper end, the exact unlinearised Theorem 2 stays below the certified bound
(minimum ratio 1.0000000023, an independent referee check).

**a0 criterion (Chain 9).** a0 is the largest a such that for every a' ≤ a: (i) every input is certified at a',
(ii) U < 2, (iii) d2(a') < 1/2. Conditions (ii) and (iii) hold on all of (0, 0.03]; (i) binds at 0.03.

## 4. The certificates

All runs used interval (ball) arithmetic in arb via python-flint (precision 128 bits or more), with explicit
truncation tails; parts of WP2b, WP3 and the WP5 sphere covers use IEEE floating point with outward rounding
(Section 7). Each run was sealed (predictions, pass bars and the instrument file hashed) before its data existed.
Seal ids refer to the lab ledger; the sealed predictions and verdicts are in `results/`.

| Run | What it encloses, for every a in (0, 0.03] | Key certified outputs | Seal (claim; amendment) |
|---|---|---|---|
| WP2a | Leading core profiles on the window: the lattice constant c6, the profile maximum P*, the leading Ricci profile R6* | c6 in [15476.93654158358624852 ± 4e-19]; P* in [4.86705830003957, 4.86705830003983]; core sup of the leading residual ≤ 14123.85; R6* ≤ 1695905.922 | 2d19816000 |
| WP2b | Curvature inputs on all of X: positivity defect e, the curvature perturbation ||dR||, B+ and C0 | e ≤ 8.34e-5; a^2 ||dR|| ≤ 3.27e-3; a^2 C0 ≤ 2.0000353; lower bound 1.99999862 at a = 0.03 | cdf666b973 |
| WP3 | Residual and Ricci outside the cores (neck and bulk, |z| ≥ a) | C_out ≤ 17226.21; R_out ≤ 65064.3; certified lower value sup |F8bulk| ≥ 16918.28 | 4158aba4f8; eddd0d7e7b |
| WP5 Stage A | Full core residual and Ricci, all remainders included (30 a-intervals) | C_core_F ≤ 14416.38; C_core_R ≤ 1,731,732.5 | 055009634d; b4c1dabedd |
| WP5 Stage B | Theorem 2 in arb on 30 a-intervals; smallness; leading coefficient; a0 | Theorem 1 | 055009634d; b4c1dabedd |

The WP3 amendment changed only two refinement-stopping tolerances before any WP3 result existed. The WP5 amendment
restated one control one-sided before any WP5 code existed. The WP5 instrument was sealed as 891fafa2d9.

The global residual constant is set by the bulk next to the Voronoi facet centre (C_q = 17226.21); the global Ricci
constant is set by the core (C_R ≤ 1,731,732.5).

## 5. Sealed tests

Every test was sealed before its run. Settled: 8 HOLD, 3 FAILED, 4 VOID (the repository's `SCORECARD.md` and
`scorecard.svg`, generated from the verdict files). The FAILED tests are the probe K1 (one exponent outside its window),
the probe K1c (an inner cut-off that cost too much Ricci curvature, which led to the window design used here) and the
certificate run WP2a, which was scored FAILED as sealed because two of its predictions used rounded constants as
containment targets; its exact enclosures are valid and are what the proof uses. Three VOID seals were superseded
before any run or code; the fourth is the first Stage 2 test, BAR-S2 (Section 1). The prediction V5 of WP5 is labelled POST HOC: it passed only after the comparator was re-read.

## 6. Verification

### 6.1 Blind replication

Three fenced independent seats, without access to the certificate code, its results or its sealed predictions, rebuilt
the three load-bearing pieces from the written specification. Two seats (one sealed test, verdict REPLICATED, all four
predictions held, controls passed) covered:

- the assembly: leading coefficient 2269.2912 (against the certified 2269.27), d2(0.03) = 0.12197295 with smallness
  1.1612989 < 2 (their faithful 2κ continuity step; the certificate's conservative 4κ gives 1.1694);
- the bulk residual peak that sets C_q: the enclosure [−16918.37305985292092614861 ± 4.88e-21] of the leading bulk
  residual at (1/4, 0, 0, 0), inside the WP3 certified interval, and a sampled maximum of its absolute value over the
  Voronoi cell of 16918.37306, below the certified C_out ≤ 17226.21.

Disclosed exposures: one seat's specification showed the leading formula and an old table; the other's showed c6
and an old sampled bulk value. One seat read an earlier script and a source page after producing its results; the
other read a point generator after its enclosure.

A third seat (a second sealed test, verdict REPLICATED, all three predictions held, controls passed) rebuilt the full
core (WP5 Stage A) with its own method (an exact Taylor polynomial of the regular potential, jet arithmetic and polar
Taylor-model boxes). At a = 0.03 its refined sampled core maxima are 1.695996527e6 for |Ric|_g/a^6 and 14125.12 for
|F|/a^8, against the certificate's own samples 1,695,994 and 14125.11 and below the certified 1,731,732.5 and
14416.38; its independent rigorous bound on the window zone for a in [0.029, 0.03] is |Ric|_g/a^6 ≤ 1.746877e6, above
every sample. Disclosed: that seat's specification printed sampled values for this window (not used as targets), and
its results were produced before its instrument file was pinned; the seat reports a from-scratch rerun on the final
file that is byte-identical, which the lead did not repeat.

Not replicated: the WP2b curvature certificate and the WP3 upper-bound cover.

### 6.2 Clean room

Byte-identical copies of the five sealed scripts were run in fresh folders from empty state with only the launch
files present. All four results (WP2a, WP2b, WP3 with its amendment runner, WP5) were identical to the originals field
by field, excluding timing fields. Disclosed limitation: in that run WP5 read its inputs from the original folders,
whose content was identical. The release runner (`code/run_chain.py`) removes this limitation by pointing WP5 at its
own work folder.

### 6.3 Referee rounds

Nine hostile referee reports, each on a stated text; the referees re-derived the steps and ran their own checks:

1. Design skeleton: logic correct with fixes; the first approximant was wrong as built (unbounded curvature at the
   exceptional curve), which led to the cut-off.
2. Definitions v1: correct in substance; nine fixes.
3. Analysis chain v1 (first round): Route 2 mathematics correct; an earlier headline a0 = 0.047 used a diameter outside
   its valid range and was withdrawn; a rate comparison was wrong as framed and was withdrawn.
4. Analysis chain v1 (second round): every Route 2 step correct; an independent floating-point re-implementation
   reproduced d2(0.01) = 2.851030e-3, d2(0.03) = 0.1219727 and the leading coefficient 2269.266, each at or below the
   certified value.
5. Definitions v3: all eight items correct or correct with fixes.
6. The WP5 assembly: correct; six printed bounds had been rounded down and were corrected (erratum); V5 labelled POST HOC.
7. Analysis chain v2 and Definitions v4: mathematics correct; secondary prints rounded the wrong way, an incomplete
   rests-on list and a missing credit line, all fixed in v2.1 and v4.1 by mechanical, count-checked edits.
8. Lemma E, the binomial lemma and the multi-direction Cauchy estimate (used by WP2b): correct, with wording fixes;
   every constant the code uses re-derived in exact integer arithmetic (for example κ(1,4) = 696, κ(1,5) = 8400).
9. Lemma A7 and the core zone Z1 of WP2b: correct, with wording fixes; the Z1 numbers were rebuilt independently
   (a^2 C0 ≤ 2.0000352762 on [0.025, 0.03], matching the certified file). This report found that the documents defined
   C0 with the global e and ||dR||, which gives 2.0036; the certified value is the pointwise bound stated in
   Section 3.1, and the definition is corrected in the final text versions.

The final text versions of the chain and the definitions carry the fixes of rounds 8 and 9; their edits have not been
re-refereed.

### 6.4 Lean-checked algebraic steps

The real-number cores of two steps are checked by Lean 4 with Mathlib (`lean/K3Chain.lean`, 13 statements, and
`lean/K3Chain2.lean`, 5 statements; standard axioms only, no sorry): step S1 of Proposition C0, the Route 2 trace-to-
eigenvalue bounds and the symmetric sandwich, end to end. Two of the statements are refutations kept as controls (a
sign-flipped S1 and an exact form of a referee's square-root estimate are false). Lean checks the statements as
written; that they match the steps of the chain is our reading.

## 7. Limitations

- **Floating-point layer.** Parts of WP2b, WP3 and the WP5 sphere covers rely on IEEE round-to-nearest with outward
  rounding and Higham's dot-product bound; WP3 also uses libm powers in three places with at least 1e-9 relative slack
  and a first-order running error bound for near-term sums. Every certified line rests on this assumption.
- **Standard facts not opened.** Heat-kernel facts (positivity, symmetry, semigroup, mass 1), Kähler normal coordinates,
  Yau's existence theorem and Calabi's uniqueness, Schauder/Fredholm theory for the implicit-function step, |W(F4)| and
  its fundamental chamber (used through WP2b's W-reduction).
- **WP2a was FAILED as sealed** on its rounded containment targets W1 and W2. Its c6, P* and R6* feed C_q and C_R; the
  enclosures are valid regardless of those predictions, but the statement "every input is enclosed by a sealed run that
  held" has this exception.
- **V5 of WP5 is POST HOC.**
- **Review status.** The last text edits of the chain and the definitions (the fixes of the final referee rounds)
  are not re-refereed; the WP5 Stage A lemmas S1-S12 had one referee pass; Lemma E, the binomial lemma, the
  multi-direction Cauchy estimate and Lemma A7 had one referee pass each.
- **Blind replication of the full core** produced its results before its instrument file was pinned (disclosed above).
- **Diameter.** L5 is proved only for a ≤ 0.03; it sets a0.
- **One surface.** The result concerns one Kummer surface (Hurwitz lattice, 16 equal bubbles), not the K3 moduli.

## 8. Related work

For the Kummer K3 surface of the Hurwitz (D4) lattice we give a computer-certified (interval arithmetic) two-sided
bound (1+d)^-1 ω_a ≤ ω_CY ≤ (1−d)^-1 ω_a, valid for every bubble size a in (0, 0.03], with d(a) ≤ 2269.27 a^3 + o(a^3)
and d(0.03) ≤ 0.122. To our knowledge (arXiv and Zenodo searched to 13 Nov 2026), this is the first explicit closeness
bound with computed constants on a Kummer K3 surface. The qualitative estimate of this type is due to Kobayashi [1990,
eq. (38), Thm 18], re-proved by Lye [2023, Thm 2.9, §6] and improved in weighted norms by Jiang [2025]. The convergent
Green-operator representation and the leading Green-sum structure of the correction are due to Fu-Liu-Xiao-Yau
[arXiv:2610.11077], who state that a certified value of a0 remained open. An explicit-interval bound of the same kind,
proved by hand, exists for a special Schoen Calabi-Yau threefold (Fu-Liu-Xiao-Yau [arXiv:2610.21932], 0 < a ≤ 2^-2260).
Certified Calabi-Yau metrics on the Fermat quartic, a different K3 at a fixed Kähler class, are due to
Douglas-Gómez-Serrano-Lehmann-Platt-Qi [arXiv:2611.16908] and Liu-Zhang [arXiv:2611.15319].

Further context. Kobayashi and Lye work with an approximant of residual order a^4 (our units) and unspecified
constants; Fu, Liu, Xiao and Yau work on the square lattice and certify no number. Donaldson's and Shackleton's gluing
constructions are qualitative.

## 9. Reproducibility

The repository `certified-kummer-k3` contains the sealed scripts (path lines adapted, listed in `code/PATCHES.md`),
the result files, the sealed predictions and verdicts, the companion documents, the Lean files and a runner that
re-runs the chain in the documented order (WP2a, WP2b, WP3 through its amendment runner, WP5) and compares every
result with the shipped one. Requirements: Python 3.11, python-flint 0.9.0, numpy, mpmath. The full chain takes
about 35 to 70 minutes on one core (clean room about 36 minutes; release re-run about 67 minutes on a shared machine); the assembly alone, from the shipped inputs, about one minute.

## 10. Open doors

What this work leaves open, and what we noticed and did not follow up. "Open as far as our dated search found" reports our
own prior-art search; "UNTESTED" marks an observation of ours that we did not test.

- **a0 beyond 0.03.** a0 is set by L5 (Section 1; Section 7), not by the smallness condition (U ≤ 1.1694 < 2). A sharper
  diameter lemma, with the certified inputs (WP2b, WP3) recomputed on a larger range, could extend the window; not
  attempted. *Status:* UNTESTED observation from this work.
- **The aligned family.** Qualitative closeness for unequal sizes is due to Kobayashi [Kob] and Lye [Lye, Thm 2.9,
  under the bounded size ratio of Assumption 2.8], with constants not computed. We found no explicit constant, explicit
  parameter box or certificate uniform over lattice parameters or unequal sizes. In the parametrisation of [Jia, Thm 1.2]
  (10 + 3 × 16 = 58), aligned bubbles leave 10 + 2 + 16 = 28 parameters, the stratum a fixed-complex-structure
  (Monge-Ampère) chain can address. *Status:* open as far as our dated search found (2026-11-17; arXiv, Zenodo, web,
  citers of [Jia]).
- **Relative bubble rotations.** The full 58-parameter family is [Jia, Thm 1.2], qualitative; [Jia, Sec. 1] notes that
  fixing one complex structure does not give enough parameters. The remaining directions need a hyperkähler-triple chain
  with an explicit, a-uniform inverse bound for the linearised triple operator (Newton-Kantorovich in place of the maximum
  principle). Unchecked idea: on a hyperkähler background the triple linearisation may reduce, up to a perturbation, to
  three scalar Laplacians, whose gap Stage 2 bounds [AC]. *Status:* open as far as our dated search found (2026-11-17;
  arXiv, Zenodo, web, citers of [Jia]).
- **Other Kummer-type gluings.** The zone split of WP2b (bubble core, window, bulk) and the WP5 assembly should transfer to
  T^4/G gluings with an explicit single-bubble metric and a lattice Green term (for example T^4/Z_3 with A_2 bubbles). The
  bulk residual extremum of WP3/REP2 sits at the midpoints l/2 between neighbouring bubbles (code/rep/rep2/rep2.py);
  such symmetric-point extrema may recur in other gluings. *Status:* UNTESTED observation from this work.

## References

Only sources opened during this work are listed.

- [Kob] R. Kobayashi, Adv. Stud. Pure Math. 18-II (1990), 257-311. doi:10.2969/aspm/01820257 (eqs. (31), (38),
  (47); Thm 18).
- [LY] Li and Yau, Acta Math. 156 (1986), Thm 2.3, pp. 168-169, and the proof of Cor 2.1, p. 170.
- [Lye] Lye, Ann. Global Anal. Geom. 63:20 (2023). doi:10.1007/s10455-023-09898-w (Assumption 2.8; Thm 2.9; §6, eq. (6.10)).
- [Jia] Jiang, arXiv:2501.08491, version 2 (Thm 1.1; Remark 1.3; Thm 1.2; Sec. 1, pp. 3-4).
- [FLXY-K] J. Fu, J. Liu, Y. Xiao and S.-T. Yau, "Explicit Ricci-flat Metrics on Kummer K3 Surfaces", arXiv:2610.11077,
  version 2 (Thm 5.6; Prop 5.5; Sec 6.7). Version 1 (8 Oct 2026) lists three authors (Fu, Xiao, Yau). Version 2 (15 Oct
  2026, "Added an author") lists four on the PDF, although arXiv's metadata still shows three. *Correction 2026-11-17:
  author list updated.*
- [FLXY] Fu, Liu, Xiao and Yau, arXiv:2610.21932, version 1 (Thm 1.1; Thm 6.3).
- [DGS] Douglas, Gómez-Serrano, Lehmann, Platt and Qi, arXiv:2611.16908, version 1 (Thm 1.1).
- [LZ] Liu and Zhang, arXiv:2611.15319, version 1 (Thm A; Thm B).
- [Don] Donaldson, arXiv:1007.4218. [Sha] Shackleton, arXiv:2605.02046.
- [AC] Andrews and Clutterbuck, arXiv:1204.5079, version 1 (cited for the Stage 2 sharpened bound, Section 1).

## Appendix A. Certified bounds on all 30 a-intervals

Rounded up from `results/wp5/result_wp5.json` (`stageB`). Each bound holds for every a in the interval.

| a-interval | U (continuity step) ≤ | d2 ≤ | d2/a^3 ≤ |
|---|---|---|---|
| [0, 0.001] | 1.000003 | 2.327e-6 | 2326.6 |
| [0.001, 0.002] | 1.000024 | 1.908e-5 | 2383.9 |
| [0.002, 0.003] | 1.000081 | 6.592e-5 | 2441.3 |
| [0.003, 0.004] | 1.000196 | 1.600e-4 | 2498.9 |
| [0.004, 0.005] | 1.000392 | 3.196e-4 | 2556.6 |
| [0.005, 0.006] | 1.000692 | 5.648e-4 | 2614.6 |
| [0.006, 0.007] | 1.001124 | 9.169e-4 | 2673.0 |
| [0.007, 0.008] | 1.001715 | 1.399e-3 | 2731.7 |
| [0.008, 0.009] | 1.002497 | 2.035e-3 | 2791.0 |
| [0.009, 0.01] | 1.003501 | 2.851e-3 | 2851.0 |
| [0.01, 0.011] | 1.004763 | 3.876e-3 | 2911.7 |
| [0.011, 0.012] | 1.006323 | 5.138e-3 | 2973.2 |
| [0.012, 0.013] | 1.008221 | 6.670e-3 | 3035.9 |
| [0.013, 0.014] | 1.010503 | 8.506e-3 | 3099.7 |
| [0.014, 0.015] | 1.013218 | 1.069e-2 | 3164.9 |
| [0.015, 0.016] | 1.016423 | 1.324e-2 | 3231.7 |
| [0.016, 0.017] | 1.020178 | 1.622e-2 | 3300.3 |
| [0.017, 0.018] | 1.024551 | 1.966e-2 | 3371.0 |
| [0.018, 0.019] | 1.029619 | 2.363e-2 | 3443.9 |
| [0.019, 0.02] | 1.035468 | 2.816e-2 | 3519.5 |
| [0.02, 0.021] | 1.042197 | 3.333e-2 | 3598.1 |
| [0.021, 0.022] | 1.049920 | 3.919e-2 | 3679.9 |
| [0.022, 0.023] | 1.058767 | 4.582e-2 | 3765.5 |
| [0.023, 0.024] | 1.068891 | 5.330e-2 | 3855.3 |
| [0.024, 0.025] | 1.080471 | 6.172e-2 | 3949.9 |
| [0.025, 0.026] | 1.093720 | 7.118e-2 | 4049.8 |
| [0.026, 0.027] | 1.108890 | 8.180e-2 | 4155.6 |
| [0.027, 0.028] | 1.126286 | 9.370e-2 | 4268.3 |
| [0.028, 0.029] | 1.146279 | 0.1071 | 4388.6 |
| [0.029, 0.03] | 1.169325 | 0.1220 | 4517.6 |

## Appendix B. Companion documents and lemma map

- `docs/ANALYSIS-CHAIN-v2.2.md`: the analysis chain (Sections 0-12), with the status of every step.
- `docs/DEFINITIONS-AND-LEMMAS-v4.2.md`: definitions 1.1-1.6, Lemmas A1-A7, H, E, C, L1-L5, interfaces 4.1-4.4,
  status table (Section 6), identities used inside certified runs (Section 7) and the WP5 Stage A lemmas S1-S12
  (Appendix).
- `docs/CHANGES-v4.2.md`: the edits from the last referee rounds (pointwise C0; fixes from the reports on
  Lemma E, the binomial lemma, the multi-direction Cauchy estimate and Lemma A7).

These are the lab's working documents, reproduced verbatim except that absolute local paths were replaced by
`<lab>`. They refer to lab files (scripts, logs, page images, referee reports) of which only the scripts, results and
verdicts are included in the repository.

| Step of the proof | Paper | Chain | Definitions |
|---|---|---|---|
| Approximant, symmetry, smoothness | 2 | 8 | 1.5, Lemma P, Lemma Q, L2 |
| Residual and Ricci of order a^8, a^6 | 2, 3.1 | 2 (I1, I2) | L3, L4 (a)-(f), 4.1-4.3, Appendix S1-S12 |
| Curvature lower bound C0 | 3.1 | 2 (I3−) | L1, L4 (g), Lemma M, Lemma M2 |
| Diameter | 2, 3.1 | 2 (I4), 10 | L5 |
| Positivity e | 3.1 | 2 (I5) | 1.6, 4.4 |
| Path and Lemma 2.1 | 3.2 | 2 | - |
| Heat-kernel oscillation bound | 3.3 | 3 | - |
| C^0, C^2, eigenvalues | 3.4-3.6 | 4, 5.2, 6 | - |
| Continuity closure, a0 | 3.7 | 7, 9 | - |

---

*Prepared with AI assistance; what was checked, and by what, is stated above.*
