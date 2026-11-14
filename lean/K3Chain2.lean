import Mathlib

/-!
# K3 chain, part 2: algebraic steps of Theorem 2 and Proposition C0 (ANALYSIS-CHAIN, Route 2)

1. Symmetric sandwich (sec. 6 / Theorem 2): if the eigenvalues satisfy e^{-eta}/M <= mu <= M, then with
   d := M e^eta - 1 we have M <= 1/(1 - d) (when d < 1) via M (1 - d) <= 1, and e^{-eta}/M = 1/(1 + d).
2. Proposition C0, step S1: for positive eigenvalues with product e^h, the trace is at least 2 e^{h/2}
   and the inverse trace is at least 2 e^{-h/2}.
-/

namespace K3Chain2

open Real

/-- Theorem 2, symmetric form, upper side: M (1 - d) <= 1 where d = M e^eta - 1, for M >= 1, eta >= 0. -/
theorem sandwich_upper {M η : ℝ} (hM : 1 ≤ M) (hη : 0 ≤ η) :
    M * (1 - (M * exp η - 1)) ≤ 1 := by
  have he : 1 ≤ exp η := one_le_exp hη
  have hM0 : 0 ≤ M := le_trans zero_le_one hM
  have h1 : M * M ≤ M * M * exp η := by
    have := mul_le_mul_of_nonneg_left he (mul_nonneg hM0 hM0)
    linarith
  nlinarith [sq_nonneg (M - 1)]

/-- Theorem 2, symmetric form, lower side: e^{-eta}/M = 1/(1 + d) where d = M e^eta - 1. -/
theorem sandwich_lower {M η : ℝ} (hM : 0 < M) :
    exp (-η) / M = 1 / (1 + (M * exp η - 1)) := by
  have he : 0 < exp η := exp_pos η
  rw [exp_neg]
  field_simp
  ring

/-- Consequence: if d < 1 the upper eigenvalue bound M is at most 1/(1 - d). -/
theorem sandwich_upper_div {M η : ℝ} (hM : 1 ≤ M) (hη : 0 ≤ η)
    (hd : M * exp η - 1 < 1) : M ≤ 1 / (1 - (M * exp η - 1)) := by
  have hpos : 0 < 1 - (M * exp η - 1) := by linarith
  rw [le_div_iff₀ hpos]
  exact sandwich_upper hM hη

/-- Proposition C0, S1 (AM-GM with det = e^h): mu1 + mu2 >= 2 e^{h/2}. -/
theorem trace_lower {μ₁ μ₂ h : ℝ} (h1 : 0 < μ₁) (h2 : 0 < μ₂) (hp : μ₁ * μ₂ = exp h) :
    2 * exp (h / 2) ≤ μ₁ + μ₂ := by
  have hs : exp (h / 2) * exp (h / 2) = μ₁ * μ₂ := by
    rw [← exp_add, hp]; ring_nf
  nlinarith [sq_nonneg (μ₁ - μ₂), exp_pos (h / 2), sq_nonneg (μ₁ + μ₂ - 2 * exp (h / 2))]

/-- Proposition C0, S1 for the inverse metric: 1/mu1 + 1/mu2 >= 2 e^{-h/2}. -/
theorem inv_trace_lower {μ₁ μ₂ h : ℝ} (h1 : 0 < μ₁) (h2 : 0 < μ₂) (hp : μ₁ * μ₂ = exp h) :
    2 * exp (-h / 2) ≤ 1 / μ₁ + 1 / μ₂ := by
  have hp' : (1 / μ₁) * (1 / μ₂) = exp (-h) := by
    rw [exp_neg, ← hp]; field_simp
  have := trace_lower (one_div_pos.mpr h1) (one_div_pos.mpr h2) hp'
  have e : -h / 2 = (-h) / 2 := by ring
  rw [e]; exact this

end K3Chain2
