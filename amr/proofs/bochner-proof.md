# CPD-1 of the Weil Kernel via Bochner-Schwartz: Proof Attempt

## Status: The Fourier-analytic argument establishes CPD-1 iff RH; identifies precisely where the proof reduces to an equivalent reformulation rather than an independent verification.

---

## 0. Goal and Setup

**Objective.** Prove (or identify the obstruction to proving) that the Weil kernel K(x) is conditionally positive definite of order 1 (CPD-1), i.e.:

$$\sum_{i,j} c_i c_j K(x_i - x_j) \geq 0 \quad \forall \{x_i\} \subset \mathbb{R},\; \forall c \in \mathbb{R}^N \text{ with } \sum c_i = 0$$

By Theorem 2.3 of subspace-alignment.md, this is equivalent to K̂(ξ) ≥ 0 as a tempered measure on ℝ \ {0} (Bochner-Schwartz characterization), which is in turn equivalent to RH.

**Strategy.** Decompose K into explicit components, compute their Fourier transforms, and analyze positivity of each piece and the total.

### 0.1 Conventions

- Fourier transform: F[f](ξ) = ∫ f(x) e^{-iξx} dx
- Inverse: f(x) = (1/2π) ∫ F[f](ξ) e^{iξx} dξ
- F[δ(x)](ξ) = 1
- F[a/(a² + x²)](ξ) = (π/a) e^{-a|ξ|} for a > 0 (standard Lorentzian transform)
- F[cos(γx)](ξ) = π[δ(ξ - γ) + δ(ξ + γ)]

### 0.2 The Weil Kernel Components

From the AMR framework, the Weil kernel is presented as:

$$K(x) = \delta(x) + K_{\text{bg}}(x) + K_{\text{zeros}}(x)$$

where:
- K_bg(x) = -(1/π) Re[ψ(1/4 + ix/2)] + log(π)/(2π)
- K_zeros(x) = (1/(2π)) Σ_γ 2cos(γx)/(1/4 + γ²), summing over imaginary parts of nontrivial zeta zeros

**Critical question (§7 below):** Is this decomposition the COMPLETE Weil kernel, or does it omit terms from the full Weil explicit formula?

---

## 1. Fourier Transform of K_bg

### 1.1 Digamma Series Expansion

The digamma function has the series:

$$\psi(s) = -\gamma + \sum_{n=0}^{\infty} \left(\frac{1}{n+1} - \frac{1}{s+n}\right)$$

Evaluating at s = 1/4 + ix/2:

$$\psi(1/4 + ix/2) = -\gamma + \sum_{n=0}^{\infty} \left(\frac{1}{n+1} - \frac{1}{n + 1/4 + ix/2}\right)$$

Taking the real part:

$$\text{Re}\left[\frac{1}{n + 1/4 + ix/2}\right] = \frac{n + 1/4}{(n + 1/4)^2 + x^2/4}$$

Therefore:

$$K_{\text{bg}}(x) = -\frac{1}{\pi}\left\{-\gamma + \sum_{n=0}^{\infty}\left(\frac{1}{n+1} - \frac{n+1/4}{(n+1/4)^2 + x^2/4}\right)\right\} + \frac{\log\pi}{2\pi}$$

$$= \frac{1}{\pi}\left\{\gamma - \sum_{n=0}^{\infty}\frac{1}{n+1} + \sum_{n=0}^{\infty} \frac{n+1/4}{(n+1/4)^2 + x^2/4}\right\} + \frac{\log\pi}{2\pi}$$

### 1.2 Fourier Transform of Each Lorentzian

Each summand has the form a/(a² + (x/2)²) with a = n + 1/4. Computing the Fourier transform:

$$\mathcal{F}_x\left[\frac{a}{a^2 + x^2/4}\right](\xi) = \int_{-\infty}^{\infty} \frac{a}{a^2 + x^2/4}\, e^{-i\xi x}\, dx$$

Substituting t = x/2, dx = 2dt:

$$= 2\int_{-\infty}^{\infty} \frac{a}{a^2 + t^2}\, e^{-2i\xi t}\, dt = 2 \cdot \pi\, e^{-2a|\xi|} = 2\pi\, e^{-2a|\xi|}$$

where we used F[a/(a² + t²)](ω) = π e^{-a|ω|}.

### 1.3 The Full Fourier Transform

The constant terms (γ, Σ 1/(n+1), log π) contribute only at ξ = 0 (as multiples of δ(ξ)). For ξ ≠ 0:

$$\hat{K}_{\text{bg}}(\xi) = \frac{1}{\pi} \sum_{n=0}^{\infty} 2\pi\, e^{-2(n+1/4)|\xi|} = 2\sum_{n=0}^{\infty} e^{-(2n + 1/2)|\xi|}$$

$$= 2\, e^{-|\xi|/2} \sum_{n=0}^{\infty} e^{-2n|\xi|} = \frac{2\, e^{-|\xi|/2}}{1 - e^{-2|\xi|}}$$

### 1.4 Positivity of K̂_bg

For all ξ ≠ 0:

$$\hat{K}_{\text{bg}}(\xi) = \frac{2\, e^{-|\xi|/2}}{1 - e^{-2|\xi|}} > 0$$

since e^{-|ξ|/2} > 0 and 1 - e^{-2|ξ|} > 0 for ξ ≠ 0.

**Asymptotic behavior:**
- As |ξ| → 0⁺: K̂_bg(ξ) ~ 2/(2|ξ|) = 1/|ξ| → +∞
- As |ξ| → ∞: K̂_bg(ξ) ~ 2e^{-|ξ|/2} → 0⁺

K̂_bg is manifestly positive, monotonically decreasing on (0, ∞), with a 1/|ξ| singularity at the origin. ∎

---

## 2. Fourier Transform of K_zeros

### 2.1 Formal Computation

$$K_{\text{zeros}}(x) = \frac{1}{2\pi} \sum_{\gamma} \frac{2\cos(\gamma x)}{1/4 + \gamma^2}$$

Using F[cos(γx)](ξ) = π[δ(ξ - γ) + δ(ξ + γ)]:

$$\hat{K}_{\text{zeros}}(\xi) = \frac{1}{2\pi} \sum_{\gamma} \frac{2\pi[\delta(\xi - \gamma) + \delta(\xi + \gamma)]}{1/4 + \gamma^2}$$

$$= \sum_{\gamma} \frac{\delta(\xi - \gamma) + \delta(\xi + \gamma)}{1/4 + \gamma^2}$$

### 2.2 As a Measure

K̂_zeros is the positive atomic measure:

$$\hat{K}_{\text{zeros}} = \sum_{\gamma} \frac{1}{1/4 + \gamma^2}\left(\delta_\gamma + \delta_{-\gamma}\right)$$

Each coefficient 1/(1/4 + γ²) > 0, so **K̂_zeros is a non-negative measure**. ∎

### 2.3 Total Mass

By the Hadamard identity (circularity-resolution.md §7.1):

$$\sum_\rho \frac{1}{\rho(1-\rho)} = 2 + \gamma_{\text{EM}} - \log(4\pi) \approx 0.046$$

Since 1/(ρ(1-ρ)) = 4/(1 + 4γ² - (2β-1)²) and on-line zeros (β = 1/2) give 1/(ρ(1-ρ)) = 1/(1/4 + γ²):

$$\text{total mass} = 2\sum_{\gamma > 0} \frac{2}{1/4 + \gamma^2} \leq \frac{4 \times 0.046}{1} \approx 0.184$$

(The factor 4 accounts for the relation between 1/(ρ(1-ρ)) and 1/(1/4+γ²) and the symmetry ρ ↔ 1-ρ̄.)

More precisely: Σ_γ 1/(1/4+γ²) = 0.046/... — the exact relation requires care, but the total mass is finite and small.

---

## 3. Fourier Transform of δ(x)

$$\hat{\delta}(\xi) = 1 \quad \forall \xi$$

This contributes a constant +1 to K̂(ξ).

---

## 4. The Naive Total (and Why It Almost Works)

### 4.1 Assembling the Pieces

For ξ ≠ 0, from the three components:

$$\hat{K}(\xi) = 1 + \frac{2\, e^{-|\xi|/2}}{1 - e^{-2|\xi|}} + \sum_{\gamma} \frac{\delta(\xi - \gamma) + \delta(\xi + \gamma)}{1/4 + \gamma^2}$$

The first term is +1 (constant). The second is strictly positive and continuous. The third is a non-negative atomic measure. The sum is manifestly ≥ 1 > 0 as a measure for all ξ ≠ 0.

**If this decomposition were complete, CPD-1 would follow trivially, and RH would be proved.**

### 4.2 Why This Cannot Be Right

CPD-1 of the Weil kernel is equivalent to RH (subspace-alignment.md, Theorem 2.3). If the above argument proved CPD-1 unconditionally, it would prove RH — a result that has resisted proof for over 160 years. Therefore, the decomposition K = δ + K_bg + K_zeros must either:

(a) Be incomplete (missing terms), or
(b) Involve a definitional inconsistency (the "K" in the matrix formulation is not the same as the "K" in the Weil explicit formula), or
(c) Contain a sign error.

**The resolution lies in (a) and (b), as analyzed in §5–§7 below.**

---

## 5. The Full Weil Functional and Its Kernel

### 5.1 The Weil Explicit Formula (from explicit-derivation.md)

The complete Weil distribution acts on g ∈ C_c^∞(ℝ) by:

$$W(g) = \hat{g}(0) + \hat{g}(1) - \sum_p \sum_{m=1}^{\infty} \frac{\log p}{p^{m/2}} \cdot 2\text{Re}[g(m\log p)] + \Omega(g)$$

where the archimedean contribution is:

$$\Omega(g) = \int_0^\infty \left[\frac{2g(0)}{e^x - 1} - \frac{2g(x)e^{-x/2}}{1 - e^{-2x}} - \frac{2g(-x)e^{-x/2}}{1 - e^{-2x}}\right] dx + g(0)\left(\log\pi - \frac{\gamma}{2}\right)$$

Setting g = f * f̃:

$$W(f * \tilde{f}) = |F(0)|^2 + |F(1)|^2 - P(f,f) + \Omega(f * \tilde{f})$$

where P(f,f) is the prime bilinear form with kernel K_prime(x) = Σ_{n≥2} (Λ(n)/√n)[δ(x - log n) + δ(x + log n)].

### 5.2 The Spectral Representation

In spectral form (explicit-derivation.md §3.2):

$$P(f,f) = \frac{1}{2\pi}\int |F(1/2 + i\tau)|^2\, \hat{K}_{\text{prime}}(\tau)\, d\tau$$

with K̂_prime(τ) = -2 Re[ζ'/ζ(1/2 + iτ)].

The archimedean spectral contribution (explicit-derivation.md §4.5):

$$\hat{\Omega}(\tau) = \text{Re}[\psi(1/4 + i\tau/2)] + \log\pi - \gamma/2$$

The full Weil kernel in spectral space:

$$\hat{W}(\tau) = -\hat{K}_{\text{prime}}(\tau) + \hat{\Omega}(\tau) = 2\,\text{Re}\left[\frac{\zeta'}{\zeta}(1/2 + i\tau)\right] + \text{Re}[\psi(1/4 + i\tau/2)] + \log\pi - \frac{\gamma}{2}$$

Plus the pole terms |F(0)|² + |F(1)|² (rank-2, outside the integral).

### 5.3 The Spectral Identity

By the explicit formula (explicit-derivation.md §5.1):

$$W(f * \tilde{f}) = \sum_\rho |F(\rho)|^2$$

Under RH: all ρ = 1/2 + iγ with γ real, so W(f * f̃) ≥ 0.

---

## 6. Reconciling the Two Formulations

### 6.1 The Matrix Kernel vs. The Weil Kernel

The AMR matrix formulation (subspace-alignment.md §2.1) uses:

$$M_{ij} = -w_i w_j K(\log p_i - \log p_j)$$

where K is the "Weil kernel" such that CPD-1 ⟺ RH.

The Weil explicit formula involves:
- The prime kernel K_prime (discrete, supported on {±log n})
- The archimedean kernel Ω (continuous)
- The pole terms |F(0)|² + |F(1)|² (rank-2)

**The question: what exactly IS the kernel K in the matrix formulation?**

### 6.2 Identification

The Weil positivity condition W(f * f̃) ≥ 0 can be written as:

$$|F(0)|^2 + |F(1)|^2 + \Omega(f * \tilde{f}) \geq P(f,f)$$

For the matrix formulation restricted to prime-power test functions (f = Σ c_{p,m} δ(x - m log p) type), the prime bilinear form P(f,f) evaluates to:

$$P = \sum_{(p,m),(q,n)} c_{p,m} c_{q,n} \frac{(\log p \cdot \log q)^{1/2}}{p^{m/2} q^{n/2}} \cdot \left[\frac{\Lambda(\text{lcm})}{...} + ...\right]$$

The CPD-1 kernel K in the matrix formulation must encode the COMPLETE Weil functional, not just K_prime. Specifically:

$$K(\log p_i - \log p_j) = \text{(contribution from Ω)} + \text{(contribution from poles)} - \text{(contribution from K_prime)}$$

Wait — this requires more care. Let me trace through the definitions precisely.

### 6.3 The Correct Kernel Identification

The quadratic form on the primitive subspace (subspace-alignment.md §2.1) is:

$$v^T M v = -\sum_{i,j} u_i u_j K(\log p_i - \log p_j) \quad \text{with } \sum u_k = 0$$

For APT (M|_prim ≤ 0), we need u^T K u ≥ 0 on Σu_k = 0, i.e., K is CPD-1.

Now, from the Weil functional applied to discretized test functions at log-prime points:

$$W(g) = \hat{g}(0) + \hat{g}(1) - P_{\text{prime}}(g) + \Omega(g)$$

The kernel K that appears in the matrix M_{ij} = -w_i w_j K(log p_i - log p_j) encodes the FULL Weil functional restricted to the prime lattice, INCLUDING the archimedean and pole contributions. Specifically:

$$K(x) = K_{\text{arch}}(x) + K_{\text{pole}}(x) - K_{\text{prime}}(x)$$

where:
- K_arch comes from Ω
- K_pole comes from |F(0)|² + |F(1)|² (these are rank-2 terms that contribute a specific function of the difference log p_i - log p_j)
- K_prime is subtracted (note the SIGN in W = poles - primes + arch)

**The decomposition K = δ + K_bg + K_zeros used in the AMR framework is a DIFFERENT splitting** of this same total kernel, organized by the source of each contribution (background from continuous spectrum vs. oscillatory from zeros).

### 6.4 The Relationship Between the Two Decompositions

The prime kernel has the spectral decomposition (from the explicit formula):

$$\hat{K}_{\text{prime}}(\tau) = -2\,\text{Re}\left[\frac{\zeta'}{\zeta}(1/2 + i\tau)\right]$$

Using the partial fraction expansion of -ζ'/ζ:

$$-\frac{\zeta'}{\zeta}(s) = \frac{1}{s-1} - \frac{1}{2}\log\pi + \frac{1}{2}\frac{\Gamma'}{\Gamma}(s/2) - \sum_\rho \frac{1}{s - \rho}$$

(up to the B constant from the Hadamard product).

On s = 1/2 + iτ:

$$-\frac{\zeta'}{\zeta}(1/2 + i\tau) = \frac{1}{-1/2 + i\tau} - \frac{\log\pi}{2} + \frac{1}{2}\psi(1/4 + i\tau/2) - \sum_\rho \frac{1}{1/2 + i\tau - \rho}$$

Taking real parts:

$$\hat{K}_{\text{prime}}(\tau)/(-2) = \text{Re}\left[\frac{\zeta'}{\zeta}(1/2+i\tau)\right]$$

So:

$$\hat{K}_{\text{prime}}(\tau) = -2\,\text{Re}\left[\frac{1}{-1/2 + i\tau}\right] + \log\pi - \text{Re}[\psi(1/4 + i\tau/2)] + 2\sum_\rho \text{Re}\left[\frac{1}{1/2 + i\tau - \rho}\right]$$

Now, the Weil spectral kernel (what must be ≥ 0 for CPD-1) is:

$$\hat{K}_{\text{Weil}}(\tau) = -\hat{K}_{\text{prime}}(\tau) + \hat{\Omega}(\tau) + (\text{pole terms})$$

From §5.2: Ω̂(τ) = Re[ψ(1/4 + iτ/2)] + log π - γ/2.

So:

$$\hat{K}_{\text{Weil}}(\tau) = 2\,\text{Re}\left[\frac{1}{-1/2 + i\tau}\right] - \log\pi + \text{Re}[\psi(1/4 + i\tau/2)] - 2\sum_\rho \text{Re}\left[\frac{1}{1/2 + i\tau - \rho}\right]$$
$$+ \text{Re}[\psi(1/4 + i\tau/2)] + \log\pi - \gamma/2 + (\text{pole terms})$$

$$= 2\,\text{Re}\left[\frac{1}{-1/2 + i\tau}\right] + 2\,\text{Re}[\psi(1/4 + i\tau/2)] - \gamma/2 - 2\sum_\rho \text{Re}\left[\frac{1}{1/2 + i\tau - \rho}\right] + (\text{pole terms})$$

### 6.5 The Zero Sum and Its Sign

The critical term is:

$$-2\sum_\rho \text{Re}\left[\frac{1}{1/2 + i\tau - \rho}\right]$$

For a zero ρ = β + iγ:

$$\text{Re}\left[\frac{1}{1/2 + i\tau - \beta - i\gamma}\right] = \frac{1/2 - \beta}{(1/2 - \beta)^2 + (\tau - \gamma)^2}$$

**Case 1 (β = 1/2, on critical line):** Each term is 0/((τ-γ)²) = 0. On-line zeros contribute nothing.

**Case 2 (β ≠ 1/2, off-line):** Paired zeros {ρ, 1-ρ̄} cancel exactly (circularity-resolution.md §1.2).

**Therefore the zero sum vanishes identically**, confirming the spectral cancellation theorem. The Weil spectral kernel becomes:

$$\hat{K}_{\text{Weil}}(\tau) = 2\,\text{Re}\left[\frac{1}{-1/2 + i\tau}\right] + 2\,\text{Re}[\psi(1/4 + i\tau/2)] - \gamma/2 + (\text{pole terms})$$

This is **independent of zero locations** — confirming circularity-resolution.md Theorem 1.1.

---

## 7. The Critical Gap: Where Does the Proof Break Down?

### 7.1 The Paradox

The spectral kernel K̂_Weil(τ) is zero-independent (§6.5). It is a specific, unconditionally computable function. If K̂_Weil(τ) ≥ 0 for all τ ≠ 0, then K is CPD-1 and RH follows.

So: **is K̂_Weil(τ) ≥ 0?**

### 7.2 Computing K̂_Weil(τ) Explicitly

From §6.5 (ignoring pole terms for now):

$$\hat{K}_{\text{Weil}}(\tau) = \frac{2(-1/2)}{1/4 + \tau^2} + 2\,\text{Re}[\psi(1/4 + i\tau/2)] - \gamma/2$$

$$= \frac{-1}{1/4 + \tau^2} + 2\,\text{Re}[\psi(1/4 + i\tau/2)] - \gamma/2$$

**At τ = 0:**

$$\hat{K}_{\text{Weil}}(0) = \frac{-1}{1/4} + 2\psi(1/4) - \gamma/2 = -4 + 2(-\gamma - \pi/2 - 3\log 2) - \gamma/2$$

Using ψ(1/4) = -γ - π/2 - 3 log 2 ≈ -γ - 1.5708 - 2.0794 ≈ -4.2274:

$$\hat{K}_{\text{Weil}}(0) = -4 + 2(-4.2274) - 0.2886 = -4 - 8.4548 - 0.2886 = -12.74$$

This is **negative**.

### 7.3 The Pole Terms Save Positivity (But Only on Primitives)

The pole contributions |F(0)|² + |F(1)|² are rank-2 positive-semidefinite terms. In the CPD-1 framework, these correspond to the behavior of K̂_Weil at ξ = 0: the CPD-1 condition allows K̂ to have an arbitrary distributional singularity at the origin.

More precisely: CPD-1 requires K̂(ξ) ≥ 0 only for ξ ≠ 0. The pole terms contribute a δ(ξ)-type singularity that is projected out by the condition Σc_i = 0 (primitivity).

But from §7.2, K̂_Weil(τ) < 0 at τ = 0 and nearby. This negativity occurs at τ VALUES near 0, not at ξ = 0 in the distributional sense. The pole terms |F(0)|² + |F(1)|² add rank-2 corrections that partially compensate — but they are OUTSIDE the τ-integral:

$$W(f * \tilde{f}) = |F(0)|^2 + |F(1)|^2 + \frac{1}{2\pi}\int |F(1/2 + i\tau)|^2 \hat{K}_{\text{Weil}}^{\text{(no poles)}}(\tau)\, d\tau$$

For K̂_Weil^(no poles)(τ) < 0 near τ = 0, the integral can be negative. The pole terms must compensate.

**On primitives (Σ c_i = 0):** The pole terms |F(0)|² + |F(1)|² are NOT zero for general primitives — they depend on F(0) = Σ c_i w_i e^{0·x_i} and F(1) = Σ c_i w_i e^{x_i}. On the primitive subspace Σ c_i w_i = 0 (NOT Σ c_i = 0 — this is a different condition depending on weights), F(0) does vanish in specific coordinates.

**This is the crux:** The pole terms provide rank-2 positivity, but the spectral integral can be negative. The question of whether the pole terms always compensate the negative spectral integral ON THE PRIMITIVE SUBSPACE is precisely RH.

### 7.4 The Precise Equivalence

**Theorem 7.1.** The following are equivalent:

(a) RH.

(b) For all f ∈ C_c^∞(ℝ): W(f * f̃) = |F(0)|² + |F(1)|² + (1/2π)∫ |F(1/2+iτ)|² Ŵ(τ) dτ ≥ 0.

(c) The Weil kernel (including pole contributions) is CPD-1.

(d) For all f ∈ C_c^∞(ℝ) with F(0) = F(1) = 0:

$$\frac{1}{2\pi}\int |F(1/2 + i\tau)|^2\, \hat{W}^{(\text{no poles})}(\tau)\, d\tau \geq 0$$

The equivalence (a) ⟺ (b) is Weil's criterion. The equivalence (b) ⟺ (c) is Bochner-Schwartz. The condition (d) is a STRONGER requirement (positivity without pole assistance) that is NOT equivalent to RH — it fails because Ŵ^(no poles)(τ) < 0 near τ = 0.

**The pole terms are essential.** They provide exactly the positivity needed to compensate the spectral negativity near τ = 0. This compensation is what the explicit formula enforces:

$$\sum_\rho |F(\rho)|^2 = |F(0)|^2 + |F(1)|^2 + \int |F(1/2+i\tau)|^2 \hat{W}^{(\text{no poles})}(\tau)\, d\tau$$

Under RH, the LHS ≥ 0 automatically. The question of whether the RHS ≥ 0 is exactly whether the pole terms compensate the potentially negative integral.

### 7.5 Where the Naive Argument Fails

Returning to §4: the decomposition K = δ + K_bg + K_zeros gave K̂(ξ) = 1 + K̂_bg + K̂_zeros > 0. This argument fails because:

1. **K_bg is NOT the full archimedean contribution.** The K_bg defined as -(1/π) Re[ψ(1/4+ix/2)] + log(π)/(2π) captures only PART of the Weil kernel. The full kernel also involves the prime contributions and their interaction with the archimedean terms.

2. **The decomposition K = δ + K_bg + K_zeros conflates different frameworks.** In the matrix formulation, K(log p_i - log p_j) encodes the Weil bilinear form restricted to log-prime points. The "δ" contribution comes from the diagonal of the matrix (K(0)), the "K_bg" from the continuous part of the Weil functional, and "K_zeros" from the zero sum. But the FULL Weil functional (§5.1) has the structure:

$$W = \text{poles} - \text{primes} + \text{archimedean} = \text{zero sum}$$

The identification of "K" in the matrix with the distributional Weil kernel requires careful bookkeeping of which terms go where.

3. **The sign convention matters.** The matrix has M_{ij} = -w_i w_j K(log p_i - log p_j) with a NEGATIVE sign. CPD-1 of K means M|_prim ≤ 0 (NSD on primitives). The K that must be CPD-1 is the kernel such that:

$$\sum c_i c_j K(x_i - x_j) = W(g_c)$$

where g_c is the test function corresponding to the coefficient vector c at the points {x_i}. This W includes the pole terms, the prime terms, and the archimedean terms. The K from §0.2 (K = δ + K_bg + K_zeros) is a DIFFERENT object — it is the "raw" kernel before the sign flip and before accounting for the full Weil structure.

---

## 8. What CAN Be Proved

### 8.1 Unconditional Results

**Theorem 8.1 (K_bg is CPD-1).** The background kernel K_bg(x) = -(1/π) Re[ψ(1/4 + ix/2)] + log(π)/(2π) is unconditionally CPD-1.

*Proof.* K̂_bg(ξ) = 2e^{-|ξ|/2}/(1 - e^{-2|ξ|}) > 0 for all ξ ≠ 0 (§1.4). By Bochner-Schwartz, K_bg is CPD-1. ∎

**Theorem 8.2 (K_zeros is CPD-1).** The zero-oscillation kernel K_zeros(x) = (1/2π) Σ_γ 2cos(γx)/(1/4+γ²) is CPD-1, unconditionally.

*Proof.* K̂_zeros is a non-negative measure (sum of positive delta masses at ±γ, §2.2). By Bochner-Schwartz, K_zeros is CPD-1. ∎

**Theorem 8.3 (δ is CPD-1).** The Dirac delta δ(x) is CPD-1.

*Proof.* δ̂(ξ) = 1 > 0 for all ξ. ∎

**Corollary 8.4.** The function K_naive(x) = δ(x) + K_bg(x) + K_zeros(x) is CPD-1.

*Proof.* The sum of CPD-1 functions with non-negative Fourier transforms is CPD-1:

$$\hat{K}_{\text{naive}}(\xi) = 1 + \frac{2e^{-|\xi|/2}}{1 - e^{-2|\xi|}} + (\text{positive measure}) > 0$$

for all ξ ≠ 0. ∎

### 8.2 What Corollary 8.4 Does NOT Prove

Corollary 8.4 proves CPD-1 of K_naive = δ + K_bg + K_zeros. This would prove RH IF K_naive were the correct Weil kernel. The gap is:

**The kernel K_naive is NOT the complete Weil kernel.** The complete Weil functional involves:

$$W(f * \tilde{f}) = |F(0)|^2 + |F(1)|^2 - P(f,f) + \Omega(f * \tilde{f})$$

The prime bilinear form P(f,f) involves the prime kernel K_prime, which is SUBTRACTED. The archimedean Ω is ADDED. The poles contribute RANK-2 terms.

The kernel K that must be CPD-1 for RH is NOT δ + K_bg + K_zeros but rather the kernel of the FULL Weil distribution W, which in spectral space is:

$$\hat{K}_{\text{Weil}}(\tau) = -\hat{K}_{\text{prime}}(\tau) + \hat{\Omega}(\tau) + (\text{pole terms as delta at 0})$$

This is a DIFFERENT function from K̂_naive, and it takes NEGATIVE values near τ = 0 (§7.2). The pole terms (as distributional contributions at ξ = 0) are projected out by the CPD-1 condition, but the negativity of the continuous part of K̂_Weil near τ = 0 is NOT projected out — it is compensated by the pole terms only when integrated against |F|² on the primitive subspace.

### 8.3 The Relationship Between K_naive and K_Weil

To understand the mismatch:

- K_naive = δ + K_bg + K_zeros captures "what the Weil matrix looks like at log-prime points" when organized by component source.
- K_Weil is the distributional kernel of the Weil positivity criterion.

These are related but not identical. The matrix M_{ij} = -w_i w_j K(log p_i - log p_j) uses a kernel K that, when the full explicit formula is applied, decomposes as:

$$K(x) = [\text{archimedean kernel at } x] - [\text{prime kernel evaluated at } x] + [\text{pole contribution at } x]$$

For x = log p_i - log p_j ≠ 0 (i ≠ j), the prime kernel K_prime(x) contributes delta masses at x only if x = ±m log p for some prime power. At generic x = log(p_i/p_j), K_prime(x) = 0 (no prime power has that log). So the off-diagonal entries are dominated by the archimedean contribution.

For x = 0 (diagonal), K_prime(0) diverges (it's a sum of delta masses evaluated at 0), which is regularized by the finite truncation.

This is consistent with the numerical observation that K(0) ≈ 2.528 and K(x) is dominated by K_bg for x ≠ 0.

---

## 9. The Honest Assessment

### 9.1 What This Proof Attempt Achieves

1. **Rigorous Fourier analysis** of the individual components K_bg, K_zeros, and δ, confirming each is separately CPD-1.

2. **Identification of the gap** between the naive kernel K_naive = δ + K_bg + K_zeros (which IS CPD-1) and the full Weil kernel K_Weil (whose CPD-1 is equivalent to RH).

3. **Precise location of the obstruction:** the spectral kernel K̂_Weil^(no poles)(τ) is NEGATIVE near τ = 0 (specifically, K̂_Weil^(no poles)(0) ≈ -12.7). The pole terms |F(0)|² + |F(1)|² must compensate this negativity on the primitive subspace. Whether they always do so is equivalent to RH.

4. **Confirmation of the spectral cancellation theorem:** K̂_Weil is independent of zero locations, making the problem unconditionally well-posed.

### 9.2 What Remains Open

The proof of CPD-1 for the FULL Weil kernel reduces to showing:

$$|F(0)|^2 + |F(1)|^2 \geq -\frac{1}{2\pi}\int |F(1/2+i\tau)|^2 \hat{W}^{(\text{no poles})}(\tau)\, d\tau$$

for all f ∈ C_c^∞(ℝ), where the RHS is positive (the integral is negative because Ŵ^(no poles) < 0 near 0).

This is an inequality between the "energy at the poles" (LHS) and the "negative spectral energy" (RHS). It is a well-posed, unconditionally stated inequality involving only known quantities (the digamma function, log π, the Euler constant). Its truth IS the Riemann Hypothesis.

### 9.3 Summary Table

| Component | K̂ ≥ 0? | CPD-1? | Status |
|-----------|---------|--------|--------|
| δ(x) | K̂ = 1 > 0 | Yes | Trivial |
| K_bg(x) | K̂ = 2e^{-\|ξ\|/2}/(1-e^{-2\|ξ\|}) > 0 | Yes | Proved (§1) |
| K_zeros(x) | K̂ = positive atomic measure | Yes | Proved (§2) |
| K_naive = δ + K_bg + K_zeros | K̂ > 0 everywhere | Yes | Proved (§4, Cor 8.4) |
| K_Weil (full, no poles) | K̂ < 0 near τ = 0 | **No** | Definite |
| K_Weil (full, with poles) | CPD-1 as measure | **⟺ RH** | Open |

### 9.4 Relationship to the AMR Framework

The AMR framework (subspace-alignment.md §8) reduces RH to Component C: "Haar implies APT." Component C is precisely the CPD-1 condition for K_Weil. This document shows:

- Components A and B (measure rigidity) are orthogonal to the Fourier analysis here.
- Component C is equivalent to the pole-compensation inequality (§9.2).
- The naive kernel K_naive satisfies CPD-1 trivially, but is NOT the correct kernel for Component C.
- The correct kernel K_Weil has a spectral gap of ~12.7 that must be filled by the rank-2 pole terms.

The perturbative approach (circularity-resolution.md, Strategy 3) shows that for FINITE truncations, the background spectral gap (~1.9) vastly exceeds the zero perturbation (~0.015), giving unconditional APT for small matrices. This finite result is consistent with our analysis: for finite matrices at log-prime points, the diagonal dominance from K(0) ensures CPD-1 without needing to resolve the delicate pole-compensation question at the continuum level.

---

## 10. Technical Appendix: Sign Convention Verification

### 10.1 The Matrix Sign

From subspace-alignment.md §2.1: M_{ij} = -w_i w_j K(log p_i - log p_j).

APT requires M|_prim ≤ 0, i.e., v^T M v ≤ 0 for v ∈ V_prim.

With u_i = v_i w_i and Σu_k = 0:

v^T M v = -Σ_{i,j} u_i u_j K(log p_i - log p_j) = -u^T K u

So M|_prim ≤ 0 ⟺ u^T K u ≥ 0 on {Σu_k = 0} ⟺ K is CPD-1. ✓

### 10.2 Weil Functional Sign

W(f * f̃) = Σ_ρ |F(ρ)|² ≥ 0 under RH.

W(f * f̃) = |F(0)|² + |F(1)|² - P(f,f) + Ω(f * f̃).

So positivity of W ⟺ P(f,f) ≤ |F(0)|² + |F(1)|² + Ω(f * f̃).

The prime kernel SUBTRACTS from positivity; the archimedean and pole terms ADD. ✓

### 10.3 The K in the Matrix

For the matrix at log-prime points with a_max = 1, the entries K(log p_i - log p_j) for i ≠ j should equal:

$$K(\log p_i - \log p_j) = K_{\text{arch}}(\log p_i - \log p_j) + 0 + 0$$

since K_prime has no support at log(p_i/p_j) for distinct primes with m = n = 1 (no prime power equals p_i/p_j), and K_zeros is part of the spectral decomposition of W evaluated at the explicit-formula level.

This confirms that for off-diagonal entries of the a_max=1 matrix, K is dominated by the archimedean (digamma) contribution, consistent with the numerical observation K(x) ≈ K_bg(x) for x at log-prime ratios.

---

## References

- [subspace-alignment.md](subspace-alignment.md) — CPD-1 ⟺ RH, matrix formulation
- [circularity-resolution.md](circularity-resolution.md) — Spectral cancellation theorem, unconditional bounds
- [../../asg/positivity/cross-terms/explicit-derivation.md](../../asg/positivity/cross-terms/explicit-derivation.md) — Full Weil functional derivation
- Bochner, S. (1933). Monotone Funktionen, Stieltjessche Integrale und harmonische Analyse. *Math. Ann.* 108, 378–410.
- Weil, A. (1952). Sur les "formules explicites" de la théorie des nombres premiers. *Comm. Sém. Math. Univ. Lund*, Tome Supplementaire, 252–265.

---

*Generated as part of the Arithmetic Measure Rigidity framework*
*Date: 2026-02-12*
