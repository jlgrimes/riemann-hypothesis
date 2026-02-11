# The Structure of Prime-Prime Cross-Terms

## The Core Problem

The Arithmetic Positivity Theorem (APT) reduces to understanding the **cross-terms** in the arithmetic intersection pairing. This document formulates the problem precisely and identifies the key structural features.

---

## 1. The Exact Formulation

### 1.1 Weil Positivity Reformulation

APT is equivalent to Weil's positivity criterion:

**For all f ∈ C_c^∞(ℝ₊), the Weil functional satisfies W(f * f̃) ≥ 0.**

where f̃(x) = f̄(1/x)/x and:

W(g) = ĝ(1/2) + ĝ(-1/2) - Σ_p Σ_{m=1}^∞ (log p / p^{m/2}) · 2Re[g(m log p)] + Ω(g)

where Ω contains archimedean terms.

### 1.2 Expanding the Convolution

For g = f * f̃:

g(x) = ∫₀^∞ f(t) · f̄(te^{-x}) · e^{-x/2} dt

(in multiplicative coordinates on ℝ₊). Then:

g(m log p) = ∫₀^∞ f(t) · f̄(t · p^{-m}) · p^{-m/2} dt

### 1.3 The Bilinear Form

W(f * f̃) defines a bilinear form on test functions:

B(f, g) = W(f * g̃) = ⟨f, g⟩_{poles} - ⟨f, g⟩_{primes} + ⟨f, g⟩_{arch}

where:
- ⟨f,g⟩_{poles} = f̂(1/2)ĝ(1/2) + f̂(-1/2)ĝ(-1/2) (positive semi-definite, rank 2)
- ⟨f,g⟩_{primes} = Σ_p Σ_m (log p / p^{m/2}) · ∫ f(t)ḡ(t·p^{-m}) p^{-m/2} dt (THIS IS THE KEY)
- ⟨f,g⟩_{arch} = archimedean contribution

**APT = B(f,f) ≥ 0 for all f** = B is positive semi-definite.

### 1.4 Decomposing the Prime Term

The prime bilinear form decomposes as:

⟨f,f⟩_{primes} = Σ_p (DIAGONAL_p) + Σ_{p≠q} (CROSS_{p,q})

where:

DIAGONAL_p = Σ_m (log p / p^m) · |∫ f(t) · e^{-iγ m log p} dt|²
           = Σ_m (log p / p^m) · |f̂(γ)|² evaluated at certain points

Wait — this isn't quite right. Let me be more careful.

⟨f,f⟩_{primes} = Σ_p Σ_m (log p / p^{m/2}) · 2Re[(f*f̃)(m log p)]
                = Σ_p Σ_m (log p / p^{m/2}) · 2Re[∫ f(t) f̄(t·p^{-m}) p^{-m/2} dt]

The diagonal contribution from prime p involves only f evaluated at points related to powers of p:

DIAG_p = Σ_m (log p / p^m) · ∫ |f(t)|² dt + Σ_m (log p / p^m) · ∫ f(t)[f̄(t·p^{-m}) - f̄(t)] p^{-m/2} dt

The first part is purely diagonal (involves |f|² only). The second part involves f at t and at t·p^{-m} — this "reaches across" scales by factors of p^m.

### 1.5 Cross-Terms Between Primes

The cross-terms CROSS(p,q) for p ≠ q arise when we organize the sum differently. Consider:

⟨f,f⟩_{primes} involves evaluating f*f̃ at points {m log p : m ≥ 1, p prime}.

Two primes p, q contribute cross-terms when m log p ≈ n log q for some m, n, i.e., when p^m ≈ q^n.

**Key observation:** The cross-terms are controlled by how well powers of different primes can approximate each other. By the abc conjecture / Baker's theorem on linear forms in logarithms:

|m log p - n log q| ≥ C(ε) / (mn)^{1+ε}

This LOWER BOUND on the distance means the cross-terms involve f evaluated at DIFFERENT points, and the product f(m log p) · f(n log q) has controlled magnitude.

---

## 2. The Matrix Formulation

### 2.1 The Weil Matrix

Define the infinite matrix M with rows/columns indexed by pairs (p, m) (prime p, positive integer m):

M_{(p,m),(q,n)} = -(log p · log q)^{1/2} / (p^{m/2} · q^{n/2}) · K(m log p - n log q)

where K(x) is the kernel determined by the archimedean and pole contributions:

K(x) = (transform of pole and archimedean contributions)

**APT ⟺ M is negative semi-definite** (on the appropriate subspace).

### 2.2 Diagonal Dominance?

If M were **diagonally dominant** — meaning |M_{ii}| ≥ Σ_{j≠i} |M_{ij}| — then M would be negative semi-definite (by the Gershgorin circle theorem applied to -M).

The diagonal entries are:

M_{(p,m),(p,m)} = -(log p) / p^m · K(0)

The off-diagonal entries involve K(m log p - n log q) which is small when m log p ≠ n log q (because K decays).

**Question:** Does diagonal dominance hold?

For large primes, the diagonal entries are ~ -(log p)/p while the off-diagonal entries involve 1/(p^{1/2}q^{1/2}). The off-diagonal entries are SMALLER than the diagonal ones for large p, q. This suggests diagonal dominance holds "eventually."

For small primes (p = 2, 3, 5, ...), the analysis is more delicate and must be done explicitly.

### 2.3 The Spectral Gap

If M has a **spectral gap** (all eigenvalues bounded away from 0), then APT holds with room to spare. The spectral gap would measure "how far from failing" APT is.

Numerically, we expect the spectral gap to exist but to shrink as we include more primes. The rate of shrinking determines whether APT holds in the limit.

---

## 3. Three Key Observations

### 3.1 Observation 1: Log-Free Powers Never Match

For p ≠ q, the equation p^m = q^n has NO solutions (by unique prime factorization). Therefore m log p ≠ n log q always, and the cross-terms are always evaluated at DISTINCT points.

The minimum distance |m log p - n log q| is bounded below by Baker's theorem on linear forms in logarithms:

|m log p - n log q| ≥ exp(-C · log m · log n · log p · log q)

This is very small but POSITIVE. The cross-terms are non-zero but "spread out."

### 3.2 Observation 2: The Cross-Terms Oscillate

The kernel K(x) oscillates (it involves the function ψ(1/4 + ix/2) from the digamma function). This means the cross-terms M_{(p,m),(q,n)} can be positive or negative.

**Key:** If the oscillation causes cancellation in the sum Σ_{(q,n)} M_{(p,m),(q,n)}, then the row sums of M are small, and diagonal dominance (approximately) holds.

### 3.3 Observation 3: The Sum Is Convergent

The total contribution of all cross-terms:

Σ_{(p,m)≠(q,n)} |M_{(p,m),(q,n)}| ≤ Σ_{p,q} Σ_{m,n} (log p · log q) / (p^{m/2} q^{n/2}) < ∞

because Σ_p Σ_m (log p)/p^{m/2} = -ζ'(1/2)/ζ(1/2) (in a regularized sense). So the matrix M has absolutely convergent off-diagonal sums.

---

## 4. The Key Inequality

### 4.1 What We Need to Prove

For ALL vectors v = (v_{p,m}) with Σ v_{p,m} = 0 (primitivity):

Σ_{(p,m),(q,n)} v_{p,m} · M_{(p,m),(q,n)} · v_{q,n} ≤ 0

### 4.2 What We Know

- The diagonal is correct: Σ_{(p,m)} v_{p,m}² · M_{(p,m),(p,m)} ≤ 0 ✓
- The matrix is symmetric ✓
- The matrix has convergent row sums ✓
- The matrix reduces to a known negative-definite matrix in the function field case ✓

### 4.3 What We Don't Know

- Whether the off-diagonal terms can overwhelm the diagonal
- The exact magnitude of cancellation in the off-diagonal sums
- Whether the spectral gap survives in the infinite-dimensional limit

---

## 5. Approaches to Proving the Inequality

### Approach A: Diagonal Dominance

Show: |M_{(p,m),(p,m)}| ≥ Σ_{(q,n)≠(p,m)} |M_{(p,m),(q,n)}|

This requires: (log p)/p^m ≥ Σ_{q≠p} Σ_n (log p · log q)^{1/2} / (p^{m/2} q^{n/2}) · |K(m log p - n log q)|

The right side is bounded by:

(log p)^{1/2} / p^{m/2} · Σ_q Σ_n (log q)^{1/2} / q^{n/2} · sup|K|
= (log p)^{1/2} / p^{m/2} · C

where C = Σ_q Σ_n (log q)^{1/2}/q^{n/2} · sup|K| is a convergent constant.

For large p: (log p)/p^m >> (log p)^{1/2}/p^{m/2} (since p^{m/2} >> 1). ✓ Diagonal dominance holds.

For p = 2, m = 1: need (log 2)/2 ≥ C · (log 2)^{1/2} / √2. This depends on C.

**This is a computable check.** If C ≤ (log 2)^{1/2} / √2, diagonal dominance holds for all (p,m) and APT is proved.

### Approach B: Positive-Definite Kernel

If the function x ↦ K(x) is a **positive-definite function** (Fourier transform is non-negative), then the matrix M with entries K(x_i - x_j) is positive semi-definite. With the appropriate signs, this gives APT.

K is positive-definite iff K̂(ξ) ≥ 0 for all ξ.

By the explicit formula: K̂(ξ) is related to the pair correlation of zeros. If the pair correlation is (1 - (sin πξ / πξ)²) (GUE), then K̂(ξ) ≥ 0 for |ξ| ≤ 1 (Montgomery's result) but may not be non-negative everywhere.

**This approach works if we can prove K is positive-definite.** But K̂(ξ) ≥ 0 for all ξ is itself equivalent to RH (or very close to it).

### Approach C: Fourier Analytic / Tauberian

Use Tauberian methods to show that the slowly-varying part of the cross-terms is non-positive, and the oscillatory part cancels.

Decompose K(x) = K_slow(x) + K_osc(x) where K_slow is the "mean" and K_osc oscillates with mean zero.

If K_slow(x) ≤ 0 for x ≠ 0, then the cross-terms have negative mean contribution, and the oscillatory parts contribute zero on average. A large deviation bound would then give APT.

---

## 6. The Most Promising Path

Based on this analysis, **Approach A (diagonal dominance)** is the most concrete. It reduces APT to verifying a finite computation (checking that C ≤ (log 2)^{1/2}/√2) plus an asymptotic argument for large primes.

The computation of C requires evaluating the kernel K and summing over primes. This is the task for the computational agents.

If diagonal dominance fails marginally, a **perturbative approach** might work: show that the off-diagonal terms, while individually possibly positive, cancel in aggregate due to the oscillation of K. This would be a "cancellation" argument, requiring deeper analysis of the arithmetic of log p / log q ratios.
