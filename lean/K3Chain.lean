import Mathlib

/-!
# K3Chain: the real-number core of two steps of the K3 analysis chain

Source: <lab>/ANALYSIS-CHAIN.md, section 4 (Proposition C0, step S1) and
section 6 (trace to eigenvalues, Route 2), and <lab>/seats/REFEREE-WP4.md item 4.

Setting (dimension 2): `μ₁, μ₂ > 0` are the eigenvalues of `ω_t` relative to `ω_a`,
so `T = μ₁ + μ₂` (trace), `S = 1/μ₁ + 1/μ₂`, and `μ₁ μ₂ = det = exp h` with `|h| ≤ η`.
-/

open Real

namespace K3Chain

/-! ## Helper: AM-GM with the product written as `exp h` -/

/-- Two positive numbers with product `exp h` have sum at least `2 exp (h/2)`. -/
lemma two_exp_half_le_sum {μ₁ μ₂ h : ℝ} (h1 : 0 < μ₁) (h2 : 0 < μ₂)
    (hp : μ₁ * μ₂ = exp h) : 2 * exp (h / 2) ≤ μ₁ + μ₂ := by
  have ha : exp (h / 2) * exp (h / 2) = exp h := by
    rw [← Real.exp_add]; congr 1; ring
  have hpos : 0 < exp (h / 2) := exp_pos _
  nlinarith [sq_nonneg (μ₁ - μ₂), sq_nonneg (μ₁ + μ₂ - 2 * exp (h / 2))]

/-! ## (a) Proposition C0, step S1 -/

/-- **Proposition C0, step S1 (real-number core).**
With `T = μ₁ + μ₂`, `S = 1/μ₁ + 1/μ₂`, `μ₁ μ₂ = exp h`, `|h| ≤ η`:
`T - 2 ≥ 2 e^{-η/2} - 2 = -η₁` and `2 - S ≤ 2 - 2 e^{-η/2} = η₁`,
where `η₁ = 2 (1 - e^{-η/2})`. -/
theorem propC0_S1 {μ₁ μ₂ h η : ℝ} (h1 : 0 < μ₁) (h2 : 0 < μ₂)
    (hp : μ₁ * μ₂ = exp h) (hη : |h| ≤ η) :
    (μ₁ + μ₂) - 2 ≥ 2 * exp (-η / 2) - 2 ∧
    2 - (1 / μ₁ + 1 / μ₂) ≤ 2 - 2 * exp (-η / 2) := by
  obtain ⟨hl, hu⟩ := abs_le.mp hη
  constructor
  · have hT := two_exp_half_le_sum h1 h2 hp
    have hexp : exp (-η / 2) ≤ exp (h / 2) := exp_le_exp.mpr (by linarith)
    linarith
  · have hinv : (1 / μ₁) * (1 / μ₂) = exp (-h) := by
      rw [div_mul_div_comm, one_mul, hp, Real.exp_neg, one_div]
    have hS := two_exp_half_le_sum (one_div_pos.mpr h1) (one_div_pos.mpr h2) hinv
    have hexp : exp (-η / 2) ≤ exp (-h / 2) := exp_le_exp.mpr (by linarith)
    linarith

/-- The same statement written with `η₁ = 2 (1 - e^{-η/2})` exactly as in the document. -/
theorem propC0_S1_eta1 {μ₁ μ₂ h η : ℝ} (h1 : 0 < μ₁) (h2 : 0 < μ₂)
    (hp : μ₁ * μ₂ = exp h) (hη : |h| ≤ η) :
    let η₁ := 2 * (1 - exp (-η / 2))
    (μ₁ + μ₂) - 2 ≥ -η₁ ∧ 2 - (1 / μ₁ + 1 / μ₂) ≤ η₁ := by
  intro η₁
  obtain ⟨hA, hB⟩ := propC0_S1 h1 h2 hp hη
  constructor <;> simp only [η₁] <;> linarith

/-- `η₁ ≤ η₂` for `η ≥ 0`, where `η₂ = 2 (e^{η/2} - 1)` (the conservative constant the code uses). -/
theorem eta1_le_eta2 {η : ℝ} (hη : 0 ≤ η) :
    2 * (1 - exp (-η / 2)) ≤ 2 * (exp (η / 2) - 1) := by
  have hprod : exp (η / 2) * exp (-η / 2) = 1 := by
    rw [← Real.exp_add]; ring_nf; simp
  have hge : 1 ≤ exp (η / 2) := Real.one_le_exp (by linarith)
  have hpos : 0 < exp (-η / 2) := exp_pos _
  nlinarith [sq_nonneg (exp (η / 2) - 1)]

/-! ## (b) Route 2: trace to eigenvalues -/

/-- One eigenvalue `μ` (with partner `ν`) of a pair with sum `≤ 2c` and product `p ≥ q > 0`
lies between the two roots of `x² - 2c x + q`. -/
lemma root_bounds {μ ν p q c : ℝ} (hμ : 0 < μ) (hs : μ + ν ≤ 2 * c)
    (hp : μ * ν = p) (hq : q ≤ p) (hq0 : 0 < q) :
    μ ≤ c + √(c ^ 2 - q) ∧ q / (c + √(c ^ 2 - q)) ≤ μ := by
  have key : (μ - c) ^ 2 ≤ c ^ 2 - q := by
    nlinarith [mul_le_mul_of_nonneg_left hs hμ.le]
  have hnn : 0 ≤ c ^ 2 - q := le_trans (sq_nonneg _) key
  have hr := Real.abs_le_sqrt key
  obtain ⟨hlo, hhi⟩ := abs_le.mp hr
  have hsq : √(c ^ 2 - q) ^ 2 = c ^ 2 - q := Real.sq_sqrt hnn
  have hν : 0 < ν := by
    by_contra hcon
    push Not at hcon
    nlinarith [mul_nonpos_of_nonneg_of_nonpos hμ.le hcon]
  have hc : 0 < c := by linarith
  have hM : 0 < c + √(c ^ 2 - q) := by positivity
  refine ⟨by linarith, ?_⟩
  rw [div_le_iff₀ hM]
  nlinarith [mul_le_mul_of_nonneg_right (show c - √(c ^ 2 - q) ≤ μ by linarith) hM.le]

/-- **Route 2, eigenvalue step.** If `μ₁ + μ₂ ≤ T = 2c` (here `c = 1 + ε₂`) and
`μ₁ μ₂ = e^h` with `|h| ≤ η`, then with `M = c + √(c² - e^{-η})` (the document's `M₂`):
`max μ ≤ M` and `min μ ≥ e^{-η}/M`. -/
theorem route2_eigen {μ₁ μ₂ h η c M : ℝ} (h1 : 0 < μ₁) (h2 : 0 < μ₂)
    (hs : μ₁ + μ₂ ≤ 2 * c) (hp : μ₁ * μ₂ = exp h) (hη : |h| ≤ η)
    (hM : M = c + √(c ^ 2 - exp (-η))) :
    max μ₁ μ₂ ≤ M ∧ exp (-η) / M ≤ min μ₁ μ₂ := by
  obtain ⟨hl, -⟩ := abs_le.mp hη
  have hq : exp (-η) ≤ exp h := exp_le_exp.mpr (by linarith)
  have hq0 : 0 < exp (-η) := exp_pos _
  obtain ⟨a1, b1⟩ := root_bounds h1 hs hp hq hq0
  obtain ⟨a2, b2⟩ := root_bounds h2 (by linarith : μ₂ + μ₁ ≤ 2 * c)
    (by rw [mul_comm]; exact hp) hq hq0
  subst hM
  exact ⟨max_le a1 a2, le_min b1 b2⟩

/-- **Referee's algebraic check (item 4).** For `η ≥ 0`:
`2M - M² e^η ≤ 1 - (M - 1)² ≤ 1`. -/
theorem referee_symmetric_algebra (M η : ℝ) (hη : 0 ≤ η) :
    2 * M - M ^ 2 * exp η ≤ 1 - (M - 1) ^ 2 ∧ 1 - (M - 1) ^ 2 ≤ 1 := by
  have he : 1 ≤ exp η := Real.one_le_exp hη
  constructor
  · nlinarith [sq_nonneg M]
  · nlinarith [sq_nonneg (M - 1)]

/-- **Route 2, symmetric form.** With `d = M e^η - 1 < 1` and `η ≥ 0`,
`e^{-η}/M = (1+d)⁻¹` and `M ≤ (1-d)⁻¹`; so any `μ` with `e^{-η}/M ≤ μ ≤ M`
satisfies `(1+d)⁻¹ ≤ μ ≤ (1-d)⁻¹`. -/
theorem route2_symmetric {μ M η d : ℝ} (hη : 0 ≤ η) (hMpos : 0 < M)
    (hd : d = M * exp η - 1) (hd1 : d < 1)
    (hlo : exp (-η) / M ≤ μ) (hup : μ ≤ M) :
    (1 + d)⁻¹ ≤ μ ∧ μ ≤ (1 - d)⁻¹ := by
  have heq : (1 + d)⁻¹ = exp (-η) / M := by
    rw [hd, Real.exp_neg]; field_simp; ring
  have hM1 : M ≤ (1 - d)⁻¹ := by
    have hpos : 0 < 1 - d := by linarith
    rw [inv_eq_one_div, le_div_iff₀ hpos]
    have := (referee_symmetric_algebra M η hη).1
    have := (referee_symmetric_algebra M η hη).2
    rw [hd]; nlinarith
  exact ⟨heq ▸ hlo, le_trans hup hM1⟩

/-- **Route 2, end to end**: both eigenvalues lie in `[(1+d)⁻¹, (1-d)⁻¹]`. -/
theorem route2_end_to_end {μ₁ μ₂ h η c M d : ℝ} (h1 : 0 < μ₁) (h2 : 0 < μ₂)
    (hs : μ₁ + μ₂ ≤ 2 * c) (hp : μ₁ * μ₂ = exp h) (hη : |h| ≤ η)
    (hM : M = c + √(c ^ 2 - exp (-η)))
    (hd : d = M * exp η - 1) (hd1 : d < 1) :
    ((1 + d)⁻¹ ≤ μ₁ ∧ μ₁ ≤ (1 - d)⁻¹) ∧ ((1 + d)⁻¹ ≤ μ₂ ∧ μ₂ ≤ (1 - d)⁻¹) := by
  have hη0 : 0 ≤ η := le_trans (abs_nonneg h) hη
  obtain ⟨hmax, hmin⟩ := route2_eigen h1 h2 hs hp hη hM
  have hMpos : 0 < M := lt_of_lt_of_le h1 (le_trans (le_max_left _ _) hmax)
  exact ⟨route2_symmetric hη0 hMpos hd hd1 (le_trans hmin (min_le_left _ _))
            (le_trans (le_max_left _ _) hmax),
         route2_symmetric hη0 hMpos hd hd1 (le_trans hmin (min_le_right _ _))
            (le_trans (le_max_right _ _) hmax)⟩

/-! ## Route 1 symmetric form (section 6): `e^η (1 - d₁²) ≤ 1` -/

/-- With `M₁ = 1 + ε + √((1+ε)² - e^{-η})`, `ε ≥ 0`, `d₁ = M₁ - 1`: `e^η (1 - d₁²) ≤ 1`. -/
theorem route1_symmetric {ε η M d : ℝ} (hε : 0 ≤ ε) (hη : 0 ≤ η)
    (hM : M = 1 + ε + √((1 + ε) ^ 2 - exp (-η))) (hd : d = M - 1) :
    exp η * (1 - d ^ 2) ≤ 1 := by
  have he : exp (-η) ≤ 1 := Real.exp_le_one_iff.mpr (by linarith)
  have hnn : 0 ≤ (1 + ε) ^ 2 - exp (-η) := by nlinarith
  have hsq := Real.sq_sqrt hnn
  have hr := Real.sqrt_nonneg ((1 + ε) ^ 2 - exp (-η))
  have hd' : d = ε + √((1 + ε) ^ 2 - exp (-η)) := by rw [hd, hM]; ring
  have h1 : 1 - d ^ 2 ≤ exp (-η) := by rw [hd']; nlinarith
  have hprod : exp η * exp (-η) = 1 := by rw [← Real.exp_add]; simp
  nlinarith [exp_pos η]

/-! ## Referee aside "(μ₁ - μ₂)² ≤ 8ε + 4η": true only to leading order -/

/-- Counterexample to the referee aside as an exact inequality:
`μ₁ = 4, μ₂ = 1/4, ε = 9/8, η = h = 0`. -/
theorem referee_sqrt_loss_exact_is_false :
    ¬ ∀ μ₁ μ₂ ε h η : ℝ, 0 < μ₁ → 0 < μ₂ → μ₁ + μ₂ ≤ 2 * (1 + ε) →
      μ₁ * μ₂ = exp h → |h| ≤ η → (μ₁ - μ₂) ^ 2 ≤ 8 * ε + 4 * η := by
  intro H
  have := H 4 (1 / 4) (9 / 8) 0 0 (by norm_num) (by norm_num) (by norm_num)
    (by simp) (by simp)
  norm_num at this

/-- Corrected exact form: `(μ₁ - μ₂)² ≤ 8ε + 4ε² + 4η`. -/
theorem referee_sqrt_loss_corrected {μ₁ μ₂ ε h η : ℝ} (h1 : 0 < μ₁) (h2 : 0 < μ₂)
    (hs : μ₁ + μ₂ ≤ 2 * (1 + ε)) (hp : μ₁ * μ₂ = exp h) (hη : |h| ≤ η) :
    (μ₁ - μ₂) ^ 2 ≤ 8 * ε + 4 * ε ^ 2 + 4 * η := by
  obtain ⟨hl, -⟩ := abs_le.mp hη
  have hq : exp (-η) ≤ exp h := exp_le_exp.mpr (by linarith)
  have hlin : 1 - η ≤ exp (-η) := by linarith [Real.add_one_le_exp (-η)]
  have hsq : (μ₁ + μ₂) ^ 2 ≤ (2 * (1 + ε)) ^ 2 := by
    apply pow_le_pow_left₀ (by linarith) hs
  nlinarith

/-! ## Positive control: a sign-flipped variant of S1 is FALSE, and Lean proves its negation -/

/-- The variant with `e^{+η/2}` in place of `e^{-η/2}` fails at `μ₁ = μ₂ = 1, h = 0, η = 1`. -/
theorem propC0_S1_signflip_false :
    ¬ ∀ μ₁ μ₂ h η : ℝ, 0 < μ₁ → 0 < μ₂ → μ₁ * μ₂ = exp h → |h| ≤ η →
      (μ₁ + μ₂) - 2 ≥ 2 * exp (η / 2) - 2 := by
  intro H
  have h0 := H 1 1 0 1 (by norm_num) (by norm_num) (by simp) (by simp)
  have : 1 < exp ((1 : ℝ) / 2) := Real.one_lt_exp_iff.mpr (by norm_num)
  linarith

end K3Chain

#print axioms K3Chain.propC0_S1
#print axioms K3Chain.propC0_S1_eta1
#print axioms K3Chain.eta1_le_eta2
#print axioms K3Chain.route2_eigen
#print axioms K3Chain.referee_symmetric_algebra
#print axioms K3Chain.route2_symmetric
#print axioms K3Chain.route2_end_to_end
#print axioms K3Chain.route1_symmetric
#print axioms K3Chain.referee_sqrt_loss_exact_is_false
#print axioms K3Chain.referee_sqrt_loss_corrected
#print axioms K3Chain.propC0_S1_signflip_false
