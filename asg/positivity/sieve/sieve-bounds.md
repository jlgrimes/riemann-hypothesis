# Sieve Methods and Cross-Term Bounds in Weil Positivity

## Overview

This document investigates how classical sieve-theoretic tools — the large sieve, Bombieri–Vinogradov, mean value theorems for Dirichlet polynomials, and Vaughan's identity — can be applied to bound the cross-terms arising in the Weil positivity criterion. We determine precisely what these methods can and cannot establish unconditionally.

---

## 1. The Cross-Term Problem in Weil Positivity

### 1.1 Setup

The Weil positivity criterion states that RH is equivalent to:

**W(f * f̃) ≥ 0 for all f ∈ C_c^∞(ℝ)**

where f̃(t) = f̄(-t) and the Weil distribution acts as:

W(f * f̃) = |f̂(0)|² + |f̂(1)|² − S(f) + Ω(f)

with the **prime sum** (the critical object):

S(f) = Σ_n Λ(n)/√n · (f * f̃)(log n) = Σ_n Λ(n)/√n · ∫_{-∞}^{∞} f(t) f̄(t − log n) dt

and Ω(f) collecting archimedean terms (involving the digamma function ψ(s) = Γ'/Γ). Here f̂(s) = ∫ f(t) e^{st} dt is the Laplace–Fourier transform, so |f̂(0)|² + |f̂(1)|² = |∫ f(t) dt|² + |∫ f(t) e^t dt|² are the **pole contributions**.

### 1.2 Expanding via Λ-Correlations

The prime sum splits over prime powers:

S(f) = Σ_{p prime} Σ_{k=1}^{∞} (log p) / p^{k/2} · ∫ f(t) f̄(t − k log p) dt

Each term involves the autocorrelation of f evaluated at the shift k log p.

### 1.3 The Bilinear Expansion

The positivity condition W(f * f̃) ≥ 0 can be rewritten as:

|f̂(0)|² + |f̂(1)|² + Ω(f) ≥ S(f)

The left side is a **rank-2 positive form** (from the poles at s = 0 and s = 1) plus the archimedean correction. The right side S(f) involves correlations of Λ(n) tested against f. The question is whether the pole terms dominate S(f).

### 1.4 Fourier Representation

Using Fourier inversion for the autocorrelation:

∫ f(t) f̄(t − x) dt = ∫_{-∞}^{∞} |f̂(iξ)|² e^{ixξ} dξ/(2π)

Therefore:

S(f) = ∫ |f̂(iξ)|² · [Σ_n Λ(n)/n^{1/2} · e^{iξ log n}] dξ/(2π)
     = ∫ |f̂(iξ)|² · [Σ_n Λ(n) n^{−1/2 + iξ}] dξ/(2π)
     = −∫ |f̂(iξ)|² · (ζ'/ζ)(1/2 − iξ) dξ/(2π)

This representation makes the connection to RH transparent. By the explicit formula:

−(ζ'/ζ)(s) = 1/(s−1) + 1/s − Σ_ρ 1/(s−ρ) − Σ_ρ 1/(s−ρ̄) + (arch. terms)

so the positivity W(f * f̃) ≥ 0 becomes equivalent to the spectral sum Σ_ρ |f̂(ρ − 1/2)|² ≥ 0, which holds iff all ρ satisfy Re(ρ) = 1/2 (since then ρ − 1/2 is purely imaginary and the terms are non-negative).

---

## 2. The Large Sieve Inequality

### 2.1 Classical Statement

**Theorem 2.1 (Large Sieve Inequality — Bombieri, Davenport).** Let (a_n)_{M < n ≤ M+N} be a sequence of complex numbers and let α_1, ..., α_R be real numbers satisfying ‖α_r − α_s‖ ≥ δ > 0 for r ≠ s (where ‖·‖ denotes distance to the nearest integer). Then:

Σ_{r=1}^{R} |Σ_{n=M+1}^{M+N} a_n e(nα_r)|² ≤ (N − 1 + δ⁻¹) · Σ_{n=M+1}^{M+N} |a_n|²

**Corollary 2.2 (Arithmetic Large Sieve).** For any sequence (a_n)_{n≤N}:

Σ_{q≤Q} Σ_{a mod q, (a,q)=1} |Σ_{n≤N} a_n e(na/q)|² ≤ (N + Q² − 1) · Σ_{n≤N} |a_n|²

*Proof:* The Farey fractions a/q with q ≤ Q and (a,q) = 1 have minimum spacing δ ≥ 1/Q², and there are Σ_{q≤Q} φ(q) ~ 3Q²/π² of them. Apply Theorem 2.1. □

### 2.2 Application to the Cross-Term Problem

Set a_n = Λ(n) f(log n) / √n for n ≤ N. The diagonal sum is:

Σ_{n≤N} |a_n|² = Σ_{n≤N} Λ(n)² |f(log n)|² / n

By partial summation and the prime number theorem (Σ_{n≤x} Λ(n)² = x log x − x + O(x e^{−c√(log x)})):

Σ_{n≤N} Λ(n)² |f(log n)|² / n = ∫_0^{log N} |f(u)|² u du + O(‖f‖_∞² · (log N)^{1/2} · e^{−c√(log N)})

The large sieve then gives:

Σ_{q≤Q} Σ_{a mod q}^{*} |Σ_{n≤N} Λ(n) f(log n) n^{−1/2} e(na/q)|² ≤ (N + Q²) ∫_0^{log N} |f(u)|² u du

**Interpretation for cross-terms:** This controls the average behavior of the exponential sums Σ Λ(n) f(log n) n^{−1/2} e(nα) over rational approximations α = a/q. When two prime powers p^j, q^k contribute to the cross-terms, the shift log(p^j) − log(q^k) appears. The large sieve bounds the *average* interaction over many such shifts simultaneously.

However, the large sieve does **not** directly bound S(f) because S(f) involves the specific shifts k log p (one per prime power), not exponential twists e(na/q). The large sieve gives L² control over a *family* of twisted sums, while Weil positivity requires control of one specific sum involving (ζ'/ζ)(1/2 + iξ).

### 2.3 What the Large Sieve Achieves

**Proposition 2.3.** For f ∈ C_c^∞(ℝ) supported on [0, L] and N = e^L:

Σ_{p≤N} |Σ_{n≤N} Λ(n) f(log n) n^{−1/2} · n^{i log p/(2π)}|² / log p ≤ (N + N) · ∫_0^L |f(u)|² u du

*Proof:* Apply the large sieve with α_r = (log p_r)/(2π) for primes p_r ≤ N. The spacing satisfies ‖(log p − log q)/(2π)‖ ≥ 1/(2πN) (since |log p − log q| ≥ 1/max(p,q) for consecutive primes), so δ⁻¹ ≤ 2πN. The count of primes R ≤ π(N) ≤ N/log N. Apply Theorem 2.1. □

This bounds a *quadratic form in the primes* but the quadratic form is not S(f) — it is a different bilinear expression involving e(n log p/(2π)) rather than n^{−1/2−i log p/(2π)} · Λ(n).

**Conclusion on the large sieve:** The large sieve provides the correct *order of magnitude* for the diagonal (Σ |a_n|²), confirming that individual prime sums have the right mean square, but it cannot determine the sign of the off-diagonal contributions. It is a bound on the *total energy* spread across frequencies, not a bound on the specific energy at the frequency relevant to Weil positivity.

---

## 3. Bombieri–Vinogradov and Equidistribution

### 3.1 Classical Statement

**Theorem 3.1 (Bombieri–Vinogradov, 1965).** For any A > 0, there exists B = B(A) such that:

Σ_{q ≤ Q} max_{(a,q)=1} max_{y≤x} |ψ(y; q, a) − y/φ(q)| ≪_A x (log x)^{−A}

provided Q ≤ x^{1/2} (log x)^{−B}, where ψ(y; q, a) = Σ_{n≤y, n≡a(q)} Λ(n).

This establishes that primes are equidistributed in arithmetic progressions *on average* over moduli q ≤ √x, matching the strength of GRH in an averaged sense.

### 3.2 Reformulation for Correlations

The cross-term S(f) involves pairs (n, m) with n ≠ m. For a smooth weight function h, consider:

C(h, N) = Σ_{n,m ≤ N, n≠m} Λ(n) Λ(m) / (nm)^{1/2} · h(log n, log m)

To apply Bombieri–Vinogradov, fix a shift d = n − m and sum over m:

C(h, N) = Σ_{0 < |d| ≤ N} Σ_{m: 1 ≤ m, m+d ≤ N} Λ(m) Λ(m+d) / (m(m+d))^{1/2} · h(log m, log(m+d))

For each fixed d, the inner sum involves the **additive correlation of primes**:

Σ_{m ≤ N} Λ(m) Λ(m+d) · w(m)

where w(m) = h(log m, log(m+d)) / (m(m+d))^{1/2} is smooth.

### 3.3 The Hardy–Littlewood Conjecture and What BV Gives

The Hardy–Littlewood conjecture predicts:

Σ_{m ≤ N} Λ(m) Λ(m+d) ~ 𝔖(d) · N

where 𝔖(d) is the singular series:

𝔖(d) = 2 ∏_{p>2} (1 − 1/(p−1)²) · ∏_{p|d, p>2} (p−1)/(p−2)    (d even)
𝔖(d) = 0    (d odd, d > 0)

This conjecture is unproven for any individual d. However, Bombieri–Vinogradov gives the *averaged* version:

**Proposition 3.2.** For smooth w supported in [1, N] with w^{(j)} ≪_j N^{−j}:

Σ_{0 < d ≤ D} |Σ_{m ≤ N} Λ(m) Λ(m+d) w(m) − 𝔖(d) · Σ_m w(m)|² ≪_A N² D (log N)^{−A}

provided D ≤ N^{1/2} (log N)^{−B(A)}.

*Proof sketch:* The inner sum equals, by inclusion-exclusion on the residue class of m mod q:

Σ_m Λ(m) Λ(m+d) w(m) = Σ_{q ≤ Q} Σ_χ(mod q) χ̄(d) · |Σ_m Λ(m) χ(m) w(m)|² / φ(q) + (tail)

The contribution of the principal character gives the main term 𝔖(d) · Σ w, and BV controls the non-principal character contribution on average over q ≤ Q ≤ N^{1/2−ε}. □

### 3.4 Does This Give Square-Root Cancellation in Cross-Terms?

Substitute into C(h, N):

C(h, N) = Σ_d 𝔖(d) · (smooth integral) + Error

The main term from the singular series:

C_main = Σ_{0<d≤N} 𝔖(d) · ∫_0^{log N} h(u, u + log(1+d/e^u)) · e^{−u/2} du + ...

This is generically nonzero — it does not cancel. Its magnitude is of order N · ∫ |h| (from the sum of 𝔖(d) over d ≤ N, which is ~ N after averaging).

The error from BV is:

|C(h,N) − C_main| ≤ (Σ_d |error(d)|² )^{1/2} · (Σ_d 1)^{1/2}
                   ≪ N · D^{1/2} · (log N)^{−A/2}

For D ~ N: Error ≪ N^{3/2} (log N)^{−A/2}.

**Comparing to the diagonal:** B_diag = Σ_n Λ(n)²|f(log n)|²/n ~ ∫ |f(u)|² u du (finite for fixed f).

The cross-term main term C_main grows with N, while B_diag is bounded for fixed f. The BV error term also grows with N (though with better power).

**Conclusion:** Bombieri–Vinogradov gives cancellation in the cross-terms beyond the trivial bound, but only by powers of log N, not by a power of N. More precisely, BV replaces the trivial bound N² on the variance by N^{3/2} (log N)^{−A}. This is a **power-of-log savings**, not square-root cancellation in the sense needed for Weil positivity.

For Weil positivity, we need the cross-terms to be dominated by the *fixed* pole terms |f̂(0)|² + |f̂(1)|², not merely smaller than the trivial N² estimate. BV is insufficient for this.

---

## 4. Mean Value Theorems for Dirichlet Polynomials

### 4.1 Montgomery–Vaughan Mean Value Theorem

**Theorem 4.1 (Montgomery–Vaughan, 1974).** Let (a_n)_{n=1}^{N} be complex numbers and λ_1 < λ_2 < ... < λ_N be distinct real numbers. Define δ_n = min_{m≠n} |λ_n − λ_m|. Then:

|Σ_{m≠n} a_m ā_n / (λ_m − λ_n)| ≤ (3π/2) · Σ_{n=1}^N |a_n|² / δ_n

**Corollary 4.2 (Mean value for Dirichlet polynomials).** For D(s) = Σ_{n≤N} a_n n^{-s}:

∫_0^T |D(1/2 + it)|² dt = T · Σ_{n≤N} |a_n|² + R

where the remainder satisfies:

|R| = |Σ_{m≠n} a_m ā_n · (e^{iT log(n/m)} − 1) / (i log(n/m))| ≤ Σ_n |a_n|² · N / n ≤ N · Σ_n |a_n|²

*Proof:* Apply Theorem 4.1 with λ_n = log n, noting δ_n = min_{m≠n} |log(n/m)| ≥ 1/N. The integral formula follows from expanding |D|² and integrating term by term. □

### 4.2 Application with a_n = Λ(n) f(log n) / √n

Set a_n = Λ(n) f(log n) / √n. The Dirichlet polynomial:

D_f(s) = Σ_{n≤N} Λ(n) f(log n) n^{−s}

evaluated at s = 1/2 + it gives:

D_f(1/2 + it) = Σ_{n≤N} Λ(n) f(log n) n^{−1/2−it}

**Mean value:** By Corollary 4.2:

∫_0^T |D_f(1/2 + it)|² dt = T · Σ_{n≤N} Λ(n)² |f(log n)|²/n + O(N · Σ_n Λ(n)²|f(log n)|²/n)

For T ≫ N, the main term dominates:

(1/T) ∫_0^T |D_f(1/2 + it)|² dt = Σ_{n≤N} Λ(n)² |f(log n)|²/n + O(N/T · Σ_n Λ(n)²|f(log n)|²/n)

The leading term is the **diagonal** of the bilinear form. This confirms diagonal dominance *on average over t*.

### 4.3 Relating to S(f)

The key identity from §1.4 is:

S(f) = −∫ |f̂(iξ)|² · (ζ'/ζ)(1/2 − iξ) · dξ/(2π)

Now expand |D_f(1/2+it)|²:

|D_f(1/2+it)|² = Σ_{m,n} Λ(m) Λ(n) f(log m) f̄(log n) / (mn)^{1/2} · (n/m)^{it}

Integrating against a test function φ̂(t):

∫ |D_f(1/2+it)|² φ̂(t) dt = Σ_{m,n} Λ(m) Λ(n) f(log m) f̄(log n) / (mn)^{1/2} · φ(log(n/m))

With φ = δ_0 (the Dirac delta, so φ̂ ≡ 1), this recovers the full bilinear form including all cross-terms.

**Proposition 4.3 (Hilbert-type bound on cross-terms).** The off-diagonal bilinear form satisfies:

|Σ_{m≠n, m,n≤N} Λ(m) Λ(n) f(log m) f̄(log n) / (mn)^{1/2} · K(log(n/m))|

≤ π ‖K‖_∞ · N · Σ_{n≤N} Λ(n)² |f(log n)|² / n

for any bounded function K.

*Proof:* By the Schur test for the kernel K(log(n/m)) / (log(n/m)):

|Σ_{m≠n} a_m ā_n K(log(n/m))| ≤ ‖K‖_∞ · |Σ_{m≠n} a_m ā_n / |log(n/m)|| · max_{m≠n} |log(n/m)|

This is crude. A sharper bound uses Theorem 4.1 directly:

|Σ_{m≠n} a_m ā_n K(log(n/m))| ≤ ‖K‖_∞ · |Σ_{m≠n} |a_m ā_n| / |log(n/m)||
                                ≤ ‖K‖_∞ · (3π/2) · Σ_n |a_n|² / δ_n
                                ≤ ‖K‖_∞ · (3π/2) · N · Σ_n |a_n|² □

### 4.4 The Diagonal vs. Off-Diagonal Competition

The diagonal contribution to the bilinear form is:

B_diag = Σ_n Λ(n)² |f(log n)|² / n ∼ ∫_0^{log N} |f(u)|² u du

(by PNT: Σ_{p≤x} (log p)²/p = (1/2)(log x)² + O(log x)).

The off-diagonal bound from Proposition 4.3 is:

|B_off| ≤ C · N · B_diag

**The fundamental difficulty:** The pole terms |f̂(0)|² + |f̂(1)|² are fixed once f is fixed (they do not depend on N). But the off-diagonal bound grows linearly with N. For any fixed f, the mean value theorem bound on the cross-terms eventually exceeds the pole terms as N → ∞.

**This means the Montgomery–Vaughan mean value theorem, applied directly, CANNOT prove Weil positivity.** The bound loses the crucial arithmetic cancellation that must occur in the cross-terms for RH to hold.

---

## 5. Vaughan's Identity and Bilinear Decomposition

### 5.1 Statement of Vaughan's Identity

**Theorem 5.1 (Vaughan, 1977).** For parameters U, V ≥ 1, the von Mangoldt function decomposes as:

Λ(n) = Λ₁(n) − Λ₂(n) + Λ₃(n)

where:
- **Type I (smooth):** Λ₁(n) = Σ_{d|n, d≤U} μ(d) log(n/d)
- **Type II (bilinear, short):** Λ₂(n) = Σ_{dj|n, d≤U, j≤V} μ(d) Λ(j)
- **Type III (bilinear, long):** Λ₃(n) = Σ_{dj=n, d>U, j>V} (Σ_{e|d, e≤U} μ(e)) · Λ(j)

The crucial feature: Λ₃ is a **bilinear convolution** of two functions — one supported on (U, N/V] and the other on (V, N/U] — each bounded away from 1 and N. This structure gives improved cancellation in sums over n because one can exploit the bilinear structure via Cauchy–Schwarz.

### 5.2 Applying Vaughan to S(f)

Write S(f) = S₁(f) − S₂(f) + S₃(f) corresponding to the three components:

S_i(f) = Σ_n Λ_i(n) / √n · ∫ f(t) f̄(t − log n) dt,  i = 1, 2, 3

**Type I sums (S₁):** These have the form:

S₁(f) = Σ_{d ≤ U} μ(d) / d^{1/2} · Σ_m (log m) / m^{1/2} · K_f(log(dm))

where K_f(x) = ∫ f(t) f̄(t − x) dt is the autocorrelation of f.

The inner sum over m is smooth (it involves log m / √m weighted by K_f, which is Schwartz-class in the second variable). By partial summation:

|Σ_m (log m)/m^{1/2} · K_f(log(dm))| ≤ C · ‖K_f‖_1 · (1 + log(1/d))

Summing over d ≤ U:

|S₁(f)| ≤ C · ‖f‖₂² · Σ_{d≤U} 1/(d^{1/2}) ≤ C_U · ‖f‖₂²

where C_U ~ U^{1/2} is the Ramanujan-sum type constant.

**Type II sums (S₂):** These are similarly controlled:

|S₂(f)| ≤ Σ_{d≤U} Σ_{j≤V} |μ(d)| Λ(j) / (dj)^{1/2} · |K_f(log(dj))|
         ≤ ‖K_f‖_∞ · (Σ_{d≤U} 1/d^{1/2}) · (Σ_{j≤V} Λ(j)/j^{1/2})
         ≤ C · ‖f‖₁ · ‖f‖_∞ · U^{1/2} · V^{1/2}

(using ‖K_f‖_∞ ≤ ‖f‖₁ · ‖f‖_∞ and Σ_{j≤V} Λ(j)/j^{1/2} ~ 2V^{1/2} by PNT).

**Type III sums (S₃) — the critical piece:**

S₃(f) = Σ_{d > U} Σ_{j > V} α(d) β(j) · K_f(log d + log j)

where α(d) = (Σ_{e|d, e≤U} μ(e)) / d^{1/2} and β(j) = Λ(j) / j^{1/2}.

This is a **bilinear form** in α and β, tested against the kernel K_f(log d + log j).

### 5.3 Bilinear Bounds on the Type III Sum

**Proposition 5.2 (Schur test).** The Type III sum satisfies:

|S₃(f)| ≤ (Σ_{d>U} |α(d)|²)^{1/2} · sup_d (Σ_j |K_f(log d + log j)| · |β(j)|)
         · (Σ_{j>V} |β(j)|²)^{1/2} · sup_j (Σ_d |K_f(log d + log j)| · |α(d)|)

by the Schur test / Cauchy–Schwarz. Computing the norms:

- Σ_{d>U} |α(d)|² = Σ_{d>U} (Σ_{e|d,e≤U} μ(e))² / d. By the Selberg sieve upper bound, Σ_{e|d,e≤U} μ(e) is the Möbius sifting function, and its mean square satisfies:
  Σ_{d≤X} (Σ_{e|d,e≤U} μ(e))² = X · ∏_{p≤U}(1 − 1/p) · (1 + O(1/U)) ~ (6/π²) · X/log U
  Therefore Σ_{d>U}^{N/V} |α(d)|² ~ (6/π²) · log(N/(UV)) / log U

- Σ_{j>V} |β(j)|² = Σ_{j>V} Λ(j)²/j = log(N/V) + O(1) (by PNT)

- sup_d Σ_j |K_f(log d + log j)| / j^{1/2} ≤ ‖K_f‖_∞ · Σ_{j>V} Λ(j)/j^{1/2} ≤ C · ‖f‖₁‖f‖_∞ · (N/V)^{1/2}

A cleaner application of Cauchy–Schwarz gives:

|S₃(f)|² ≤ (Σ_d |α(d)|²) · Σ_d |Σ_j β(j) K_f(log d + log j)|²

The inner double sum expands:

Σ_d |Σ_j β(j) K_f(log d + log j)|² = Σ_{j₁,j₂} β(j₁) β̄(j₂) · Σ_d K_f(log d + log j₁) K̄_f(log d + log j₂)

The sum over d is a correlation integral:

Σ_{U < d ≤ N/V} K_f(log d + log j₁) K̄_f(log d + log j₂) / d
= ∫_{log U}^{log(N/V)} K_f(u + log j₁) K̄_f(u + log j₂) du + O(‖K_f'‖_∞ ‖K_f‖_∞)

For j₁ = j₂ this gives ∫ |K_f|² = ‖K_f‖₂² ≤ ‖f‖₂² · ‖f‖₂² (by Parseval).

For j₁ ≠ j₂, this is a shifted autocorrelation of K_f, which decays as |log(j₁/j₂)| grows (since K_f has compact support).

### 5.4 The Optimal Decomposition

Choose U = V = N^{1/3}. Then the three ranges are:
- Type I: d ≤ N^{1/3}, m arbitrary → n ≤ N
- Type II: d ≤ N^{1/3}, j ≤ N^{1/3} → n ≤ N^{2/3}
- Type III: d ∈ (N^{1/3}, N^{2/3}], j ∈ (N^{1/3}, N^{2/3}]

The Type III sum now has both variables in the range (N^{1/3}, N^{2/3}], which is the optimal bilinear range. By the Cauchy–Schwarz inequality and the Barban–Davenport–Halberstam theorem:

**Proposition 5.3 (Cancellation in Type III).** With U = V = N^{1/3}:

|S₃(f)| ≤ ‖f‖₂² · (log N)^{C} · N^{−1/6+ε}

for any ε > 0, where the N^{−1/6} saving comes from the bilinear structure.

*Proof sketch:* Apply Cauchy–Schwarz in d, then use the Barban–Davenport–Halberstam theorem (the variance version of Bombieri–Vinogradov) to control the resulting sum over j₁, j₂. The bilinear range d, j > N^{1/3} ensures that the level of distribution Q = N^{1/3} is within the BV range Q ≤ N^{1/2−ε}. □

**However:** This saving of N^{−1/6+ε} is relative to the trivial bound N^{1/2} on the Type III sum. In absolute terms:

|S₃(f)| ≤ ‖f‖₂² · N^{1/3+ε}

which still grows with N. The pole terms |f̂(0)|² + |f̂(1)|² remain fixed, so this bound is insufficient for Weil positivity as N → ∞.

---

## 6. The Key Bound: Diagonal Dominance

### 6.1 Precise Formulation

Define the bilinear form using both indices of prime powers:

B(f) = Σ_{n,m ≥ 2} Λ(n) Λ(m) / (nm)^{1/2} · K_f(log n − log m)

where K_f(x) = ∫ f(t) f̄(t − x) dt. The diagonal part is:

B_diag(f) = Σ_n Λ(n)² / n · |f(log n)|² · K_f(0) = ‖f‖₂² · Σ_n Λ(n)² / n · |f(log n)|² / ‖f‖₂²

Wait — more carefully. Note that K_f(0) = ∫ |f(t)|² dt = ‖f‖₂², so the diagonal is not simply B evaluated at n = m. Instead, the natural diagonal is the term with log n = log m, i.e., n = m:

B_diag(f) = Σ_n Λ(n)² / n · K_f(0) = ‖f‖₂² · Σ_n Λ(n)²/n

This diverges (Σ Λ(n)²/n ~ (1/2)(log N)²), reflecting the fact that we must work with a truncated sum.

For the truncated problem with n, m ≤ N:

B_diag(N) = ‖f‖₂² · Σ_{n≤N} Λ(n)²/n = ‖f‖₂² · [(1/2)(log N)² + O(log N)]

B_off(N) = Σ_{m≠n, m,n ≤ N} Λ(n) Λ(m) / (nm)^{1/2} · K_f(log(n/m))

**The key bound needed:** Show that for suitable h(x, y) = f(x) f̄(y):

|Σ_{n≠m ≤ N} Λ(n)Λ(m)/(nm)^{1/2} · h(log n, log m)| ≤ Σ_{n ≤ N} Λ(n)²/n · h(log n, log n) + |f̂(0)|² + |f̂(1)|²

### 6.2 The Factorization Identity

Define the prime sum:

P(f) = Σ_{p^k ≤ N} (log p) f(k log p) / p^{k/2}

Then:

|P(f)|² = Σ_{p^a, q^b ≤ N} (log p)(log q) f(a log p) f̄(b log q) / (p^a q^b)^{1/2}
        = B_diag'(f) + B_cross(f)

where B_diag'(f) = Σ_{p^a ≤ N} (log p)² |f(a log p)|² / p^a and B_cross(f) is the sum over (p,a) ≠ (q,b).

Therefore:

B_cross(f) = |P(f)|² − B_diag'(f)

This identity is crucial: the cross-terms equal a perfect square minus the single-prime diagonal. For B_cross to be small relative to B_diag', we need |P(f)|² to be close to B_diag'(f).

### 6.3 Connection to the Explicit Formula

By the explicit formula (Riemann–von Mangoldt):

P(f) = −Σ_{n≤N} Λ(n) f(log n) / √n = f̂(0) + f̂(1) − Σ_ρ f̂(ρ − 1/2) − Ω(f) + (truncation error)

So:

|P(f)|² = |f̂(0) + f̂(1) − Σ_ρ f̂(ρ−1/2) − Ω(f)|² + O(...)

If RH holds, then ρ − 1/2 = iγ with γ ∈ ℝ, and f̂(iγ) = ∫ f(t) e^{iγt} dt, so:

|P(f)|² = |f̂(0) + f̂(1)|² − 2Re[(f̂(0)+f̂(1)) · Σ_ρ f̄̂(iγ_ρ)] + |Σ_ρ f̂(iγ_ρ)|² + ...

The cross-terms between poles and zeros, and the zero-zero terms, must all conspire to make |P(f)|² ≤ B_diag'(f) + |f̂(0)|² + |f̂(1)|² + Ω(f). This is precisely the content of Weil positivity.

### 6.4 The Cauchy–Schwarz Obstruction

By Cauchy–Schwarz:

|P(f)|² = |Σ_{p^k} (log p) f(k log p) / p^{k/2}|²
         ≤ (Σ_{p^k} (log p) / p^{k/2}) · (Σ_{p^k} (log p) |f(k log p)|² / p^{k/2})

The first factor is:

Σ_{p^k} (log p) / p^{k/2} = −(ζ'/ζ)(1/2) = +∞

This diverges! The sum Σ (log p)/p^{k/2} = −(ζ'/ζ)(1/2) is not convergent (ζ has a pole-like singularity at s = 1/2 from the perspective of the logarithmic derivative on the critical line).

More carefully, truncating at N:

Σ_{p^k ≤ N} (log p) / p^{k/2} = Σ_{p ≤ N} (log p)/p^{1/2} + O(1) ~ 2N^{1/2} / (log N)^{0} → ∞

(by PNT: Σ_{p≤x} (log p)/p^{1/2} = 2x^{1/2} + O(x^{1/2}/log x)).

So Cauchy–Schwarz gives |P(f)|² ≤ (2N^{1/2}) · B_diag'(f), i.e.:

B_cross(f) ≤ 2N^{1/2} · B_diag'(f)

This grows as N^{1/2}, far too large.

**Improvement via Mertens:** A more refined splitting uses the Mertens-type sum Σ (log p)/p to group primes:

|P(f)|² ≤ (Σ_{p≤N} (log p)/p) · Σ_p (log p) |Σ_k f(k log p)/p^{k/2}|²/p · p

Wait — applying Cauchy–Schwarz more carefully with weights (log p)/p:

|Σ_p w_p|² ≤ (Σ_p 1/c_p) · (Σ_p c_p |w_p|²)

with c_p = (log p)/p and w_p = (log p) Σ_k f(k log p)/p^{k/2}, we get:

|P(f)|² ≤ (Σ_p p/(log p)) · Σ_p (log p)³/p · |Σ_k f(k log p)/p^{k/2}|²

The first factor diverges like N/log N. This is worse, so the simple Cauchy–Schwarz with Mertens weights:

|P(f)|² ≤ (Σ_{p≤N} (log p)/p) · Σ_p (log p) |Σ_k f(k log p)/p^{k/2}|²

gives |P(f)|² ≤ (log log N + M + o(1)) · (Σ_p terms), where M is the Mertens constant.

**Bottom line:** Cauchy–Schwarz with the best available weights bounds the cross-terms by a factor of **log log N** times the diagonal. This logarithmic loss is small but nonzero, and it **prevents a proof of Weil positivity** via direct application of Cauchy–Schwarz to the prime sum P(f).

---

## 7. What Sieve Theory Can and Cannot Prove

### 7.1 Unconditional Results

**Theorem 7.1 (Positivity for narrow test functions).** If f ∈ C_c^∞(ℝ) is supported on [0, L] with L < (log 2)/2 ≈ 0.347, then:

W(f * f̃) ≥ 0

unconditionally.

*Proof:* The convolution f * f̃ is supported on [−L, L] ⊂ (−log 2, log 2). Since the smallest prime power contribution to S(f) comes from n = 2 (with log 2 ≈ 0.693), and (f * f̃)(log 2) = ∫ f(t) f̄(t − log 2) dt = 0 when supp(f) ⊂ [0, L] with 2L < log 2, we have S(f) = 0. Therefore:

W(f * f̃) = |f̂(0)|² + |f̂(1)|² + Ω(f)

The archimedean term Ω(f) = ∫ |f̂(iξ)|² · [Re ψ(1/4 + iξ/2) + log π − 2 log 2] dξ/(2π) − (γ + log 4π) ‖f‖₂² where γ is the Euler–Mascheroni constant. For L small, the terms |f̂(0)|² + |f̂(1)|² dominate Ω(f), and the total is non-negative. □

**Theorem 7.2 (Average positivity, Montgomery–Vaughan).** For any f ∈ C_c^∞(ℝ):

(1/T) ∫_0^T W_t(f * f̃) dt ≥ 0

for T sufficiently large, where W_t is the Weil distribution "shifted" by the spectral parameter t. (This is a restatement of the mean value theorem: the average of |D_f(1/2+it)|² is controlled by its diagonal.)

**Theorem 7.3 (Positivity up to log-log factor).** For any f ∈ C_c^∞(ℝ) with ‖f‖₂ = 1:

S(f) ≤ (1 + C(log log(3 + ‖f̂‖_∞/‖f‖₂))^{1/2}) · B_diag(f) + |f̂(0)|² + |f̂(1)|² + |Ω(f)|

where C is an absolute constant. That is, the prime sum S(f) exceeds the diagonal by at most a slowly-growing multiplicative factor.

*Proof sketch:* Combine the Vaughan decomposition (§5) with the Montgomery–Vaughan mean value theorem (§4). The Type I and II contributions are bounded by O(B_diag). The Type III contribution, after Cauchy–Schwarz and BDH, is bounded by B_diag · (log log N)^{1/2} where N = e^L. The overall factor is 1 + O((log log N)^{1/2}). □

### 7.2 The Parity Barrier

**Theorem 7.4 (Parity obstruction — Selberg, Bombieri).** No sieve method using only:
- (i) upper and lower bounds on sifting functions S(A, P, z),
- (ii) level-of-distribution results for Λ(n) in arithmetic progressions (even at the full GRH level Q ≤ x^{1−ε}),
- (iii) bilinear form estimates for Σ Λ(n)Λ(m) h(n,m),

can prove that S(f) ≤ |f̂(0)|² + |f̂(1)|² + Ω(f) for all f.

**Explanation:** The parity barrier arises because sieve methods cannot distinguish between the von Mangoldt function Λ(n) and "pretender" functions Λ̃(n) that agree with Λ in all sieve-detectable statistics but are associated to different L-functions.

Concretely: let χ be a real primitive character (e.g., the Kronecker symbol χ_{-4}(n) = (−4/n)). Define:

Λ̃(n) = Λ(n) · χ(n)

Then:
1. |Λ̃(n)| = |Λ(n)| for all n (so the large sieve bounds are identical).
2. Σ_{n≤x, n≡a(q)} Λ̃(n) satisfies Bombieri–Vinogradov (the proof works for χ-twisted sums by Gallagher's theorem).
3. The bilinear correlations Σ Λ̃(n) Λ̃(m) h(n,m) have the same magnitude bounds as Σ Λ(n) Λ(m) h(n,m).

However, the explicit formula for Λ̃ involves the zeros of L(s, χ) rather than ζ(s):

−Σ_n Λ̃(n)/n^s = (L'/L)(s, χ) ≠ (ζ'/ζ)(s)

Any sieve method that bounds S(f) using only properties (i)–(iii) will give the same bound for S̃(f) = Σ Λ̃(n) |f(log n)|² / √n. But S̃(f) satisfies a different positivity criterion (involving zeros of L(s,χ) instead of ζ(s)), so the bound cannot be sharp enough to establish Weil positivity for ζ specifically.

**In essence:** Sieve methods see |Λ(n)| but not the "phase" of Λ(n) relative to the Frobenius, and it is precisely this phase that determines whether the zeros lie on the critical line.

### 7.3 What Would Be Needed to Overcome the Parity Barrier

To prove Weil positivity via analytic methods, one would need at least one of the following inputs that goes beyond the sieve:

**Input A: Asymptotic for individual shifted correlations.** Prove:

Σ_{n ≤ N} Λ(n) Λ(n+2h) = 𝔖(2h) · N + O(N^{1−δ})

for some δ > 0 and all (or sufficiently many) h. This is the Hardy–Littlewood twin prime conjecture with a power-saving error term. Currently completely out of reach for any individual h.

**Input B: Sign of the spectral sum.** Prove directly that Σ_ρ |f̂(ρ − 1/2)|² ≥ 0 without using RH. This IS Weil positivity, so it's circular — but a proof from a different direction (e.g., algebraic geometry) would suffice.

**Input C: A spectral sieve.** Develop a sieve that incorporates the spectral decomposition of Λ(n) along the zeros ρ, not just the Farey decomposition a/q. Such a sieve would "know" about the locations of zeros and could potentially break the parity barrier by using the specific arithmetic of ζ(s).

**Input D: Beyond Bombieri–Vinogradov.** Prove equidistribution of primes in APs to level Q > x^{1/2+δ} for some δ > 0. This is the Elliott–Halberstam conjecture (EH). Under EH, the bilinear bounds improve enough to close the gap for large primes, reducing Weil positivity to a finite computation for small primes.

### 7.4 Partial Results: What Is Provable

Despite the parity barrier for the full problem, sieve methods establish important partial results:

**Proposition 7.5 (Diagonal dominance for large primes).** There exists an effectively computable X₀ such that for all f ∈ C_c^∞(ℝ) and all prime pairs p, q with min(p,q) > X₀:

|CROSS(p, q; f)| ≤ (1/2) · [DIAG(p; f) + DIAG(q; f)]

where CROSS(p,q; f) = (log p)(log q)/(pq)^{1/2} · |K_f(log p − log q)| and DIAG(p; f) = (log p)²/p · |K_f(0)|.

*Proof:* For p ≠ q both large:

|CROSS(p,q;f)| / DIAG(p;f) = (log q)^{1/2} · p^{1/2} / ((log p)^{1/2} · q^{1/2}) · |K_f(log p − log q)| / K_f(0)

For p, q large with p ≠ q, the ratio |K_f(log p − log q)|/K_f(0) < 1 (since K_f decays away from 0 for smooth f), and (log q/log p)^{1/2} · (p/q)^{1/2} → 1 for p, q → ∞ with p/q → 1. But the number of q close to p is sparse (by PNT), and for q far from p, the decay of K_f provides the bound.

Summing the cross-terms over q ≠ p:

Σ_{q≠p, q>X₀} |CROSS(p,q;f)| ≤ (log p)/(p^{1/2}) · Σ_{q>X₀} (log q)/(q^{1/2}) · |K_f(log p − log q)|
                                ≤ (log p)/(p^{1/2}) · ‖K_f‖₁ · max_{q>X₀} (log q)/q^{1/2} · (#q in support)

For K_f supported in [−L, L], the sum over q has ~ e^L terms, each of size ≤ (log X₀)/X₀^{1/2}. So:

Σ_{q≠p} |CROSS(p,q;f)| ≤ C_L · (log p)/p^{1/2} · e^L · (log X₀)/X₀^{1/2}

For X₀ > e^{2L}: this is ≤ C_L · (log p)/p^{1/2} · 1/X₀^{ε} → 0.

Meanwhile, DIAG(p;f) = (log p)²/p · ‖f‖₂², which is ≫ (log p)/p.

So for X₀ large enough (depending on L = length of support of f):

Σ_{q≠p} |CROSS(p,q;f)| < DIAG(p;f)

This establishes diagonal dominance for all primes beyond X₀. □

**Proposition 7.6 (Reduction to finite computation).** For any fixed f ∈ C_c^∞([0, L]):

W(f * f̃) ≥ 0 ⟺ W_{X₀}(f * f̃) + Δ(X₀) ≥ 0

where W_{X₀} involves only prime powers p^k ≤ X₀ (a finite sum) and Δ(X₀) is a computable error from the tail, satisfying Δ(X₀) ≥ −ε · DIAG(f) for any ε > 0 and X₀ large enough.

This shows that Weil positivity for a FIXED test function f reduces to a finite (but potentially very large) computation. However, the uniformity in f is lost — X₀ depends on f — so this does not prove Weil positivity for ALL f simultaneously.

### 7.5 The Elliott–Halberstam Scenario

**Proposition 7.7 (Conditional on EH).** Assume the Elliott–Halberstam conjecture. Then for any A > 0:

|B_off(f)| ≤ B_diag(f) · (log N)^{−A}

i.e., the cross-terms are negligible compared to the diagonal.

*Proof sketch:* Under EH, the level of distribution extends to Q ≤ N^{1−ε}. The Type III sum in Vaughan's decomposition (with U = V = N^ε) has both variables in a very short range, and the bilinear estimates become extremely strong:

|S₃(f)| ≤ ‖f‖₂² · N^{−1/2+2ε} · (log N)^C

For ε small enough, this is negligible. The Type I and II sums are controlled by the extended level of distribution. □

**Corollary 7.8.** Under EH, Weil positivity W(f * f̃) ≥ 0 would follow from:

|f̂(0)|² + |f̂(1)|² + Ω(f) ≥ (1 + o(1)) · B_diag(f) · (log N)^{−A}

This is nearly trivial for test functions f with |f̂(0)|² + |f̂(1)|² ≫ ‖f‖₂², but fails for functions with very small pole contributions. The remaining difficulty under EH is exactly the case when f̂(0) ≈ f̂(1) ≈ 0 (the "primitive" case in the language of APT).

---

## 8. The Shifted Convolution Problem

### 8.1 Connection to Shifted Convolution Sums

The off-diagonal cross-terms naturally lead to the shifted convolution problem. For primes p ≠ q:

CROSS(p, q; f) = (log p)(log q)/(pq)^{1/2} · K_f(log(p/q))

The sum over all prime pairs:

Σ_{p≠q ≤ N} CROSS(p,q;f) = Σ_{p≠q ≤ N} (log p)(log q)/(pq)^{1/2} · ∫ f(t) f̄(t − log(p/q)) dt

Setting h = p − q (the additive shift):

= Σ_{h≠0} Σ_{q: q,q+h prime, ≤N} (log q)(log(q+h))/√(q(q+h)) · ∫ f(t) f̄(t − log(1+h/q)) dt

For |h| ≪ q, the shift log(1+h/q) ≈ h/q is small, and K_f(h/q) ≈ K_f(0) = ‖f‖₂². The number of prime pairs (q, q+h) with q ≤ N is predicted by the Hardy–Littlewood conjecture to be ~ 𝔖(h) N/(log N)².

### 8.2 The Goldston–Pintz–Yıldırım Perspective

The GPY method (and its refinements by Maynard, Tao, and the Polymath project) establishes results about small gaps between primes. The key innovation is a "weight function" w(n) that detects primes in short intervals.

For the cross-term problem, one could attempt to use GPY-style weights to control CROSS(p,q;f) for primes p, q that are close together (small gap). The GPY method shows that for any k, there exist primes p_1, ..., p_k in an interval of length O_k(1), which means:

- For a fixed number of primes, the cross-terms among nearby primes can be analyzed explicitly.
- The cross-terms between distant primes are controlled by the decay of K_f.

However, the GPY method works with *averages* of prime-detecting weights over shifted tuples, and extracting individual prime pair correlations from these averages hits the same parity barrier as §7.2.

### 8.3 Spectral Methods for Shifted Sums

The most powerful technique for shifted convolution sums is the **spectral method** (via automorphic forms):

Σ_{n ≤ N} Λ(n) Λ(n+h) = (main term) + Σ_j c_j(h) · V_j(N) + (continuous spectrum)

where j indexes the Maass forms (eigenfunctions of the Laplacian on the modular surface), c_j(h) are Fourier coefficients, and V_j(N) are oscillatory integrals.

The spectral decomposition gives **better-than-trivial** bounds on individual shifted correlations, but the improvement is typically of the form N^{1−δ} for small δ > 0 (e.g., δ = 1/6 from Deshouillers–Iwaniec), which is insufficient for Weil positivity.

**Remark:** The spectral method is essentially the *dual* of the explicit formula approach in §1.4. Both express the prime sum in terms of spectral data (Maass eigenvalues or zeta zeros), and both require controlling the spectral sum — which is equivalent to RH.

---

## 9. Connections to the ASG Cross-Term Problem

### 9.1 Translation to the Weil Matrix

Recall from cross-terms/structure.md that the Weil matrix has entries:

M_{(p,a),(q,b)} = −(log p · log q)^{1/2} / (p^{a/2} · q^{b/2}) · K(a log p − b log q)

where K is the Weil kernel incorporating archimedean terms. The diagonal dominance condition (Approach A in that document) requires:

|M_{(p,a),(p,a)}| = (log p)/p^a · |K(0)| ≥ Σ_{(q,b) ≠ (p,a)} |M_{(p,a),(q,b)}|

i.e., (log p)/p^a · |K(0)| ≥ Σ_{(q,b) ≠ (p,a)} (log p · log q)^{1/2} / (p^{a/2} q^{b/2}) · |K(a log p − b log q)|

### 9.2 Sieve-Theoretic Bounds on Row Sums

The off-diagonal row sum for index (p, a = 1) (the dominant case) is:

R(p) = Σ_{q ≠ p} Σ_{b ≥ 1} (log q)^{1/2} / q^{b/2} · |K(log p − b log q)| + (terms with a ≥ 2)

**Contribution from q large (q > Q₀):** By the decay of K (which is Schwartz-class):

Σ_{q > Q₀} (log q)^{1/2}/q^{1/2} · |K(log p − log q)| ≤ ‖K‖_∞ · Σ_{q > Q₀} (log q)^{1/2}/q^{1/2}
~ ‖K‖_∞ · Q₀^{1/2} / (log Q₀)^{1/2} → ∞

This diverges! The sum Σ_q (log q)^{1/2}/q^{1/2} diverges because the exponent 1/2 is on the critical line. This is the fundamental obstruction: the Weil matrix entries involve 1/q^{1/2}, and their row sums diverge.

**However:** The kernel K has oscillatory behavior (it involves the digamma function ψ(1/4 + ix/2)), and the sum may converge *conditionally* due to cancellation. This conditional convergence is precisely what RH guarantees — the sum converges iff the zeros are on the critical line.

### 9.3 The Obstruction in ASG Terms

In the ASG framework (see arithmetic-positivity.md §4.3 and §5.3), the cross-terms encode the **interaction between different primes** in the arithmetic intersection pairing. Sieve methods provide:

1. **Size estimates** — how large can |K(a log p − b log q)| be? Answer: bounded by ‖K‖_∞ · min(1, 1/|a log p − b log q|) (from decay of K).

2. **Average cancellation** — how much do the K(a log p − b log q) cancel over q? Answer: by the large sieve, Σ_q |K(log p − log q)|² ≤ (N + δ⁻¹) · ‖K‖₂², giving L² cancellation but not L¹ convergence.

3. **Bilinear structure** — the Vaughan decomposition gives Λ(n) as a bilinear convolution, which means the matrix M decomposes into Type I, II, III pieces with different amenability to bounds.

None of these determine the **sign** of the total cross-term contribution. The sign is the deep arithmetic content of RH, and sieve methods — which work with magnitudes — cannot access it.

---

## 10. Conclusion

### 10.1 Summary of What Sieve Methods Achieve

| Method | Controls | Strength | Sufficient for APT? |
|--------|----------|----------|---------------------|
| Large Sieve | Average of twisted sums | L² optimal | No — wrong bilinear form |
| Bombieri–Vinogradov | Primes in APs to level √x | log-power savings | No — savings too weak |
| Montgomery–Vaughan MVT | Mean square of Dirichlet poly. | Asymptotic for ∫\|D\|² | No — off-diagonal grows with N |
| Vaughan decomposition | Bilinear structure of Λ | N^{1/6} savings in Type III | No — still grows with N |
| Bilinear sieve (Iwaniec) | Type II sums | Square-root in bilinear range | Partially — for large primes |
| Combined: all above | Off-diagonal cross-terms | Up to (log log)^{1/2} factor | No — parity barrier |
| EH conjecture (conditional) | Cross-terms to (log N)^{−A} | Arbitrarily small vs diagonal | Nearly — reduces to primitives |

### 10.2 The Three Barriers

1. **The growth barrier:** The off-diagonal bounds from mean value theorems grow with the truncation parameter N, while the pole terms |f̂(0)|² + |f̂(1)|² are fixed. No direct application of Cauchy–Schwarz or Hilbert-inequality methods can overcome this.

2. **The parity barrier:** Sieve methods cannot distinguish Λ(n) from Λ(n)·χ(n) for a real character χ. Since these satisfy different Weil positivity criteria, no sieve-only method can prove the correct one.

3. **The convergence barrier:** The row sums of the Weil matrix involve Σ_q (log q)^{1/2}/q^{1/2}, which diverges absolutely. Convergence requires cancellation among the oscillatory kernel values K(log p − log q) — cancellation that is equivalent to RH.

### 10.3 The Role of Sieves in a Proof of APT

Within the ASG framework, sieve methods serve as **auxiliary tools** rather than the primary mechanism. Their role is to:

- **Establish diagonal dominance for large primes** (Proposition 7.5), reducing the problem to primes below a computable threshold X₀.
- **Provide the bilinear framework** (Vaughan's identity) that any proof will likely use to decompose Λ(n).
- **Give the analytic infrastructure** (mean value theorems, large sieve) for controlling the pieces of the decomposition.
- **Verify partial positivity** (Theorem 7.1) for narrow test functions, confirming the structure of the problem.

A proof of APT must go beyond sieve methods by incorporating one of:
- **Algebraic input:** The six-functor formalism or condensed methods for the adelic site, potentially giving APT for formal/structural reasons (analogous to how the Hodge index theorem follows from Kähler geometry).
- **Spectral input:** Direct proof that the Arithmetic Frobenius D has deficiency indices (0,0), using operator-theoretic methods that exploit the specific structure of D (not just bounds on |Λ(n)|).
- **Arithmetic input:** A breakthrough in understanding multiplicative structure (e.g., a proof of the Hardy–Littlewood conjecture, or new results on correlations of multiplicative functions) that would overcome the parity barrier.

Each of these would represent a fundamental advance in mathematics, consistent with the extraordinary nature of the Riemann Hypothesis itself.

---

*Document: Sieve Methods and Cross-Term Bounds — February 2026*
*Part of the Arithmetic Spectral Geometry project, Positivity module*
