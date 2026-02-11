# Equivalent Formulations of the Riemann Hypothesis

## Overview

One of the remarkable features of the Riemann Hypothesis is the extraordinary number of seemingly unrelated statements that turn out to be equivalent to it. This document catalogs and analyzes these equivalences, organized by mathematical domain.

---

## Part I: Number-Theoretic Equivalences

### 1. Prime Counting Function

**RH ⟺** π(x) = Li(x) + O(√x log x)

where π(x) = #{p ≤ x : p prime} and Li(x) = ∫₂ˣ dt/ln t.

This is Riemann's original motivation: RH gives the best possible error term in the Prime Number Theorem. Without RH, the best known error term is:

π(x) = Li(x) + O(x · exp(-c(log x)^{3/5}(log log x)^{-1/5}))

The jump from this to O(√x log x) illustrates the enormous power of RH.

### 2. Chebyshev's ψ Function

**RH ⟺** ψ(x) = x + O(√x log²x)

where ψ(x) = Σ_{n≤x} Λ(n) = Σ_{p^k≤x} log p.

The explicit formula gives:

ψ(x) = x - Σ_ρ x^ρ/ρ - log(2π) - (1/2)log(1 - x^{-2})

The terms x^ρ/ρ have magnitude |x^ρ/ρ| = x^{Re(ρ)}/|ρ|. RH forces Re(ρ) = 1/2, giving each term size O(√x), which after summing gives the claimed error bound.

### 3. Möbius Function

**RH ⟺** M(x) = Σ_{n≤x} μ(n) = O(x^{1/2+ε}) for all ε > 0

More precisely, RH ⟺ the Dirichlet series 1/ζ(s) = Σ μ(n)/n^s converges for Re(s) > 1/2.

The Möbius function μ(n) encodes prime factorization:
- μ(n) = 0 if n has a squared prime factor
- μ(n) = (-1)^k if n is a product of k distinct primes

The cancellation in M(x) — roughly as much as in a random walk — is one way to "see" RH.

### 4. Mertens Function and Mertens Conjecture

The Mertens conjecture (now disproven): |M(x)| ≤ √x for all x ≥ 1.

This would imply RH, but it's actually false (Odlyzko-te Riele, 1985). However, the weaker statement M(x) = O(x^{1/2+ε}) is equivalent to RH.

### 5. Prime Gaps

**RH ⟹** p_{n+1} - p_n = O(√p_n log p_n)

where p_n is the n-th prime. The best unconditional result is p_{n+1} - p_n = O(p_n^{0.525}).

Under RH, Cramér's conjecture (p_{n+1} - p_n = O(log²p_n)) would be nearly optimal, but this remains unproven even under RH.

---

## Part II: Analytic Equivalences

### 6. Growth of 1/ζ(s)

**RH ⟺** 1/ζ(σ + it) = O(t^ε) for all ε > 0 and σ > 1/2

This says ζ has no zeros in the half-plane Re(s) > 1/2, and controls the growth of the reciprocal.

### 7. The Lindelöf Hypothesis (Implied by RH)

The Lindelöf Hypothesis (LH): ζ(1/2 + it) = O(t^ε) for all ε > 0.

**RH ⟹ LH**, but LH does not imply RH. The best unconditional bound is:

ζ(1/2 + it) = O(t^{13/84+ε})

(Bourgain, 2017). The exponent 13/84 ≈ 0.1548 is far from the conjectured 0.

### 8. Hardy's Z-Function

Define Z(t) = e^{iθ(t)}ζ(1/2 + it) where θ(t) is chosen so Z(t) ∈ ℝ for t ∈ ℝ.

**RH ⟺** All zeros of Z(t) are real.

The function Z(t) is the natural object for studying zeros on the critical line, as its sign changes correspond exactly to zeros of ζ on Re(s) = 1/2.

### 9. Li's Criterion

**RH ⟺** λ_n ≥ 0 for all n = 1, 2, 3, ...

where λ_n = (1/(n-1)!) (d/ds)^n [s^{n-1} log ξ(s)]|_{s=1}

= Σ_ρ [1 - (1 - 1/ρ)^n]

The first few values: λ₁ ≈ 0.0231, λ₂ ≈ 0.0923, λ₃ ≈ 0.2075, ...

These have been verified to be positive for n up to 10^8. The growth rate is λ_n ~ (n/2)log(n/2) under RH.

### 10. Generalized Li Coefficients

For τ ∈ ℝ, define:

λ_n(τ) = Σ_ρ [1 - ((ρ-τ)/(ρ-τ̄))^n]

**RH ⟺** λ_n(τ) ≥ 0 for all n ≥ 1 and all τ with Re(τ) = 1/2.

This family of criteria interpolates and provides additional structure.

---

## Part III: Function Space Equivalences

### 11. Nyman-Beurling Criterion

**RH ⟺** The indicator function χ_{(0,1)} is in the L²(0,∞) closure of the span of functions:

f_θ(x) = {θ/x} - θ{1/x}, 0 < θ ≤ 1

where {y} = y - ⌊y⌋ is the fractional part.

Equivalently, RH ⟺ the distance

d_N = inf_{c_k, θ_k} ‖χ_{(0,1)} - Σ_{k=1}^N c_k f_{θ_k}‖_{L²(0,∞)}

satisfies d_N → 0 as N → ∞.

### 12. Báez-Duarte's Refinement

**RH ⟺** d_N → 0 where d_N is the distance above, and moreover one can take θ_k = k/N.

This gives a concrete sequence of approximation problems whose solvability is equivalent to RH.

### 13. Volchkov's Criterion

**RH ⟺**

∫₀^∞ ∫₀^∞ (1 - 12xy) log|ζ(1/2 + it)| / ((1+4x²)(1+4y²)) dx dy = π(3 - γ)/32

where γ is the Euler-Mascheroni constant.

### 14. Riesz Criterion

**RH ⟺** Σ_{k=1}^∞ (-1)^{k+1} x^k / ((k-1)! ζ(2k)) = O(x^{1/4+ε}) as x → ∞

---

## Part IV: Arithmetic Equivalences

### 15. Robin's Inequality

**RH ⟺** For all n ≥ 5041:

σ(n) < e^γ · n · ln(ln(n))

where σ(n) = Σ_{d|n} d and γ ≈ 0.5772 is the Euler-Mascheroni constant.

The exceptions for n < 5041 are precisely characterized. The largest exception is n = 5040 = 7!.

### 16. Lagarias' Inequality

**RH ⟺** For all n ≥ 1:

σ(n) ≤ H_n + e^{H_n} · ln(H_n)

where H_n = 1 + 1/2 + ... + 1/n.

This is notable for being entirely elementary to state — no complex analysis or zeta function appears.

### 17. Nicolas' Criterion

**RH ⟺** For all k ≥ 1:

Π_{i=1}^k p_i / (Π_{i=1}^k (p_i - 1)) > e^γ · ln(ln(N_k))

where N_k = p₁p₂...p_k is the primorial and p_i is the i-th prime.

When this fails, the failure corresponds to a zero off the critical line.

---

## Part V: Operator/Spectral Equivalences

### 18. Hilbert-Pólya (Conjectural)

**RH ⟺** There exists a self-adjoint operator T on a Hilbert space H such that the eigenvalues of T are exactly {γ : ζ(1/2 + iγ) = 0}.

Self-adjointness forces the eigenvalues to be real, hence all zeros on Re(s) = 1/2.

### 19. de Branges Condition

**RH ⟺** The function E(z) = ξ(1/2 - iz) belongs to the de Branges class and generates a de Branges space with the appropriate positivity properties.

### 20. Weil's Positivity Condition

**RH ⟺** For all f ∈ C_c^∞(ℝ_{>0}):

Σ_ρ ĝ(ρ) ≥ 0

where g = f * f̃, f̃(x) = (1/x)f̄(1/x), and ĝ is the Mellin transform.

This is a positivity condition on a distribution, analogous to positive-definiteness.

---

## Part VI: The de Bruijn-Newman Constant

### 21. Λ = 0 Equivalence

Define the family of entire functions parametrized by t ∈ ℝ:

H_t(z) = ∫₀^∞ e^{tu²} Φ(u) cos(zu) du

where Φ(u) = Σ_{n=1}^∞ (2π²n⁴e^{9u} - 3πn²e^{5u}) exp(-πn²e^{4u}).

The de Bruijn-Newman constant Λ is defined by:
- H_t has only real zeros ⟺ t ≥ Λ

Key results:
- de Bruijn (1950): Λ ≤ 1/2
- Newman (1976): Λ ≥ 0 (conjectured, now proven)
- Rodgers-Tao (2020): Λ ≥ 0 (proved Newman's conjecture)
- Various: upper bounds have been reduced over time

**RH ⟺ Λ = 0 ⟺ Λ ≤ 0**

The current best upper bound is Λ ≤ 0.22 (Platt-Trudgian, 2021).

Proving RH via this route requires showing Λ = 0, i.e., the zeros of H_0 = Ξ are already all real — there is no "room to spare."

---

## Part VII: Probabilistic and Statistical Equivalences

### 22. Random Walk / Möbius Randomness

**RH ⟺** The partial sums of μ(n) behave like a random walk, i.e., M(x) = O(x^{1/2+ε}).

This is related to the "Möbius randomness law" — the idea that μ(n) is "orthogonal" to all structured sequences.

### 23. Denjoy's Probabilistic Interpretation

If the values μ(1), μ(2), μ(3), ... were independent random variables taking values -1, 0, +1 with the appropriate probabilities (probability 6/π² for ±1, and 1-6/π² for 0), then with probability 1:

Σ_{n≤x} μ(n) = O(x^{1/2+ε})

So RH asserts that μ(n) behaves "as randomly as possible."

---

## Part VIII: Structural Analysis

### Why So Many Equivalences?

The abundance of equivalent formulations reveals that RH sits at a nexus of mathematical structures:

1. **Primes ↔ Zeros**: The explicit formula creates a duality. RH constrains both sides.
2. **Analysis ↔ Arithmetic**: RH bridges continuous (zeta function) and discrete (primes) mathematics.
3. **Positivity**: Many equivalences are positivity conditions, suggesting a deep underlying positivity principle.
4. **Approximation**: Several criteria involve how well certain functions can be approximated, connecting to functional analysis.
5. **Spectral theory**: The zeros "want to be" eigenvalues of a self-adjoint operator.

### Which Equivalences Are Most Promising for a Proof?

**Tier 1 - Most structurally illuminating:**
- Weil's positivity (connects to geometry, has worked for function fields)
- Hilbert-Pólya/spectral (if the operator can be found)
- de Bruijn-Newman Λ = 0 (reduces to a single number)

**Tier 2 - Concrete and testable:**
- Li's criterion (concrete sequence, computable)
- Nyman-Beurling (functional analysis framework)
- Robin's inequality (elementary arithmetic)

**Tier 3 - Illuminating but difficult to leverage:**
- Prime distribution equivalences (essentially restate the problem)
- Möbius randomness (true but hard to prove from this angle)

### The Pattern of Barriers

Across all equivalences, the same fundamental obstacle appears in different guises:

| Equivalence | What must be proven | Why it's hard |
|-------------|-------------------|---------------|
| π(x) error term | Zeros have Re = 1/2 | Direct approach is circular |
| Li's criterion | λ_n ≥ 0 for all n | Requires uniform control over all zeros |
| Weil positivity | A distribution is positive | No geometric framework over ℤ |
| Nyman-Beurling | Density in L² | Requires precise arithmetic information |
| Λ = 0 | Heat equation zeros stay real | No monotonicity principle known |
| Robin | σ(n) bound | Would need deep understanding of divisor function |

Each equivalence provides a different lens, but the underlying hardness is invariant. A proof must engage with the deep arithmetic structure that none of these reformulations, by themselves, provide.

---

## Appendix: Table of Equivalences

| # | Name | Domain | Statement (brief) | Year |
|---|------|--------|-------------------|------|
| 1 | Prime counting | Number theory | π(x) = Li(x) + O(√x log x) | 1859 |
| 2 | Chebyshev | Number theory | ψ(x) = x + O(√x log²x) | 1859 |
| 3 | Möbius | Number theory | M(x) = O(x^{1/2+ε}) | ~1900 |
| 4 | 1/ζ growth | Analysis | 1/ζ(σ+it) = O(t^ε) for σ>1/2 | ~1900 |
| 5 | Li | Analysis | λ_n ≥ 0 for all n | 1997 |
| 6 | Nyman-Beurling | Functional analysis | Density in L² | 1950/1955 |
| 7 | Báez-Duarte | Functional analysis | Refined density | 2003 |
| 8 | Robin | Arithmetic | σ(n) < e^γ n ln ln n | 1984 |
| 9 | Lagarias | Arithmetic | σ(n) ≤ H_n + e^{H_n} ln H_n | 2002 |
| 10 | Nicolas | Arithmetic | Primorial inequality | 1983 |
| 11 | Weil | Analysis/Geometry | Positivity of distribution | 1952 |
| 12 | de Bruijn-Newman | Analysis | Λ = 0 | 1950/1976 |
| 13 | Riesz | Analysis | Growth bound on series | 1916 |
| 14 | Volchkov | Analysis | Integral identity | 1995 |
| 15 | Hilbert-Pólya | Spectral theory | Self-adjoint operator exists | ~1914 |
