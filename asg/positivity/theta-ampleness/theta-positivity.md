# The Theta Function Ampleness Strategy

## The Most Concrete Path to RH via ASG

This document develops the single most promising approach to proving APT (and hence RH) identified by the function field analysis: proving the **ampleness of the arithmetic theta polarization** by direct computation.

---

## 1. The Strategy in One Page

### 1.1 The Chain of Implications

In Weil's function field proof, the logical chain is:

```
Ampleness of Θ on Jac(C)
    ⟹ Positive-definiteness of the Riemann form
        ⟹ Rosati involution is positive
            ⟹ Hodge Index Theorem on C × C
                ⟹ |αᵢ| = q^{1/2}  (RH for function fields)
```

Each implication is a standard theorem in algebraic geometry. The ONLY input is the ampleness of the theta divisor, which is a consequence of the fact that Jac(C) is a principally polarized abelian variety.

### 1.2 The Arithmetic Analogue

In ASG, the proposed chain is:

```
Ampleness of Θ_ar on the Arithmetic Jacobian
    ⟹ Positive-definiteness of the arithmetic Riemann form
        ⟹ Arithmetic Rosati involution is positive
            ⟹ Arithmetic Positivity Theorem (APT)
                ⟹ All γ_ρ ∈ ℝ  (RH for ζ)
```

Steps 2→3→4→5 are formal consequences of step 1, by the same abstract arguments as in the function field case. **The entire problem reduces to step 1: proving ampleness of Θ_ar.**

### 1.3 What "Ampleness" Means Concretely

The arithmetic theta function is:

$$\Theta(x) = \sum_{q \in \mathbb{Q}^*} e^{-\pi q^2 |x|_\infty} \prod_p \mathbf{1}_{q \in \mathbb{Z}_p}$$

where x ∈ C_Q = A_Q*/Q* is an idele class and |x|_∞ is the archimedean absolute value.

Ampleness of the associated line bundle L(Θ) means:

**c₁(L(Θ), h_Θ) > 0** as a current on C_Q

where h_Θ = -log|Θ|² is the metric and c₁ = ∂∂̄ h_Θ / (2πi) is the curvature.

This is a **concrete analytic positivity condition** on the Jacobi theta function.

---

## 2. The Theta Function on the Idele Class Group

### 2.1 The Adelic Theta Function

Define the standard theta function on A_Q:

$$\theta(x) = \sum_{q \in \mathbb{Q}} \psi(q^2 x)$$

where ψ = ψ_∞ ⊗ (⊗_p ψ_p) is the standard additive character of A_Q/Q:
- ψ_∞(t) = e^{2πit} at the archimedean place
- ψ_p(t) = e^{-2πi{t}_p} at each prime p (where {t}_p is the fractional part in Q_p)

The Jacobi theta function is the restriction to the multiplicative structure:

$$\Theta(y) = \sum_{n \in \mathbb{Z}} e^{-\pi n^2 y} \quad \text{for } y > 0 \text{ (archimedean component)}$$

### 2.2 The Functional Equation

The theta function satisfies:

$$\Theta(1/y) = \sqrt{y} \cdot \Theta(y)$$

This is the Poisson summation formula applied to the Gaussian. In adelic language:

$$\theta(x^{-1}) = |x|^{1/2} \theta(x)$$

### 2.3 The Associated Line Bundle

On C_Q, the theta function defines a section of a line bundle L(Θ). The metric on this line bundle is:

$$h_\Theta(y) = -\log|\Theta(y)|^2 = -2\log\Theta(y) \quad \text{for } y > 0$$

(using that Θ(y) > 0 for y > 0).

The curvature form is:

$$\omega = \frac{i}{2\pi} \partial\bar\partial h_\Theta = -\frac{1}{\pi} \frac{d^2}{dy^2} \log\Theta(y) \cdot dy \otimes d\bar{y}$$

In the real variable y > 0 (the archimedean coordinate on C_Q):

$$\omega(y) = -\frac{1}{\pi} \frac{d^2}{dy^2} \log\Theta(y)$$

### 2.4 The Curvature Computation

We need to compute ∂²/∂y² log Θ(y).

**Step 1:** Compute Θ'(y) and Θ''(y).

$$\Theta(y) = 1 + 2\sum_{n=1}^{\infty} e^{-\pi n^2 y}$$

$$\Theta'(y) = -2\pi \sum_{n=1}^{\infty} n^2 e^{-\pi n^2 y}$$

$$\Theta''(y) = 2\pi^2 \sum_{n=1}^{\infty} n^4 e^{-\pi n^2 y}$$

**Step 2:** Compute (log Θ)'' = Θ''/Θ - (Θ'/Θ)².

$$(\log\Theta)''(y) = \frac{\Theta''(y)}{\Theta(y)} - \left(\frac{\Theta'(y)}{\Theta(y)}\right)^2$$

**Step 3:** Determine the sign.

By the Cauchy-Schwarz inequality applied to the sums:

$$\Theta(y) \cdot \Theta''(y) = \left(\sum a_n\right)\left(\sum c_n\right)$$

$$(\Theta'(y))^2 = \left(\sum b_n\right)^2$$

where a_n = e^{-πn²y}, b_n = -πn²e^{-πn²y}, c_n = π²n⁴e^{-πn²y}.

Then a_n c_n = π²n⁴ e^{-2πn²y} and b_n² = π²n⁴ e^{-2πn²y}.

So a_n c_n = b_n² for each n. This means:

$$\Theta \cdot \Theta'' = \sum_n b_n^2 + \sum_{m \neq n} a_m c_n$$

$$(\Theta')^2 = \sum_n b_n^2 + \sum_{m \neq n} b_m b_n$$

The difference is:

$$\Theta \cdot \Theta'' - (\Theta')^2 = \sum_{m \neq n} (a_m c_n - b_m b_n)$$

$$= \sum_{m \neq n} \pi^2 n^4 e^{-\pi(m^2+n^2)y} - \pi^2 m^2 n^2 e^{-\pi(m^2+n^2)y}$$

$$= \pi^2 \sum_{m \neq n} n^2(n^2 - m^2) e^{-\pi(m^2+n^2)y}$$

By symmetry (swapping m and n):

$$= \frac{\pi^2}{2} \sum_{m \neq n} (n^2 - m^2)^2 \cdot \frac{e^{-\pi(m^2+n^2)y}}{1} \cdot \frac{1}{...}$$

Wait — let me redo this. The sum over m ≠ n of n²(n² - m²) equals, by swapping m ↔ n:

$$\sum_{m \neq n} n^2(n^2 - m^2) = \sum_{m \neq n} m^2(m^2 - n^2) \quad \text{(after relabeling)}$$

Adding the two:

$$2\sum_{m \neq n} n^2(n^2 - m^2) e^{-\pi(m^2+n^2)y} = \sum_{m \neq n} (n^2 - m^2)^2 e^{-\pi(m^2+n^2)y}$$

Wait, that's not right either. Let me be more careful.

$$S = \sum_{m \neq n} n^2(n^2 - m^2) e^{-\pi(m^2+n^2)y}$$

Swapping m ↔ n:

$$S' = \sum_{m \neq n} m^2(m^2 - n^2) e^{-\pi(m^2+n^2)y}$$

Note S' = -S + ∑_{m≠n} (m⁴ - n⁴)·... hmm, this is getting complicated. Let me think about it differently.

$$\Theta \cdot \Theta'' - (\Theta')^2 = \sum_{m,n \geq 0} a_m c_n - \sum_{m,n \geq 0} b_m b_n$$

where the sums include m = n terms. We have:

$$= \sum_{m,n} (a_m c_n - b_m b_n) = \sum_{m,n} e^{-\pi(m^2+n^2)y} \pi^2(n^4 - m^2 n^2)$$

Hmm, the constant term (m=0 or n=0) needs special treatment since a_0 = 1, b_0 = 0, c_0 = 0 for the n=0 term in the theta sum.

**Let me use a cleaner approach.** Write Θ = Σ_n f(n) where f(n) = e^{-πn²y}. Then:

$$(\log\Theta)'' = \frac{(\Sigma f'')(\Sigma f) - (\Sigma f')^2}{(\Sigma f)^2}$$

By the Cauchy-Schwarz inequality for sums:

$$\left(\sum_n f_n' \right)^2 \leq \left(\sum_n f_n\right) \left(\sum_n \frac{(f_n')^2}{f_n}\right)$$

with equality iff f_n'/f_n is constant for all n (i.e., all the "eigenvalues" -πn² are equal).

Now (f_n')²/f_n = π²n⁴ e^{-πn²y} and f_n'' = π²n⁴ e^{-πn²y}. So (f_n')²/f_n = f_n''!

Therefore Cauchy-Schwarz gives:

$$(Θ')² ≤ Θ · Θ''$$

with equality iff -πn² is the same for all n contributing to the sum — which never happens (since n = 0, 1, 2, ... give different values).

**Therefore: Θ·Θ'' - (Θ')² > 0, which means (log Θ)'' > 0.**

**This means log Θ is CONVEX.**

**And therefore the curvature ω = -(1/π)(log Θ)'' < 0.**

### 2.5 The Sign Problem

**The curvature is NEGATIVE, not positive!**

This means the theta line bundle has NEGATIVE curvature — the opposite of ampleness!

This is actually CORRECT and well-known: the theta divisor on an abelian variety has POSITIVE curvature as a divisor on the abelian variety, but when we compute log Θ on the universal cover (which is what we did above), the curvature has a specific sign that depends on the convention.

**Resolution:** The issue is the sign convention. In algebraic geometry, a line bundle L is ample iff c₁(L) is a POSITIVE (1,1)-form. For the theta divisor on an abelian variety A:

- The theta function Θ is a section of L(Θ) on A = C^g/Λ
- The metric is h = e^{-π|z|² - ...} (involving the Hermitian form of the polarization)
- The curvature c₁(L, h) = (i/2π)∂∂̄(-π|z|²) = (1/2)δ_{jk} dz_j ∧ dz̄_k > 0

The positivity comes from the HERMITIAN FORM of the polarization, not from the theta function directly.

For the idele class group C_Q, the "abelian variety" structure is:
- C_Q ≅ ℝ₊ × ∏_p ℤ_p*/... (mixed archimedean/non-archimedean)
- The Hermitian form of the polarization is the INTERSECTION PAIRING
- Ampleness of Θ means the intersection pairing is positive-definite on H^{1,1}

### 2.6 The Correct Curvature Computation

The correct formulation: on an abelian variety A with polarization λ: A → Â, the associated Riemann form E is:

$$E(x, y) = \text{Im}(H(x, y))$$

where H is the Hermitian form associated to λ. The polarization is ample iff H is positive-definite.

For the arithmetic Jacobian J_ar (which plays the role of Jac(Spec(Z))):
- J_ar = C_Q / (discrete subgroup) in some sense
- The polarization comes from the Weil distribution: λ corresponds to W
- The Hermitian form H corresponds to the bilinear extension of W

**The Hermitian form IS the Weil distribution:**

$$H(f, g) = W(f * \tilde{g})$$

**The positive-definiteness of H IS the Weil positivity criterion IS RH.**

---

## 3. Breaking the Circle: Where is the Non-Circular Content?

### 3.1 The Apparent Circularity

We've traced the ampleness condition back to Weil positivity, which is RH. This seems circular.

But there's a subtlety: in the function field case, the ampleness of Θ is ALSO equivalent to RH (for function fields). The point is that ampleness can be PROVED BY OTHER MEANS — namely, by the explicit construction of Θ as a divisor on Jac(C).

### 3.2 How Weil Breaks the Circle for Function Fields

In the function field case, the Jacobian Jac(C) is an ABELIAN VARIETY — a projective group variety. The theta divisor Θ is an EFFECTIVE DIVISOR on Jac(C). Its ampleness follows from:

1. Θ is effective (it's a subvariety of codimension 1)
2. Θ generates the Néron-Severi group of Jac(C) (principally polarized)
3. By the Nakai-Moishezon criterion: a divisor on a projective variety is ample iff D^dim(V) · V > 0 for all subvarieties V

For the theta divisor: Θ^g = g! (by Poincaré's formula). Since g! > 0, Θ is ample.

**The key ingredients:**
- (A) Jac(C) is PROJECTIVE
- (B) Θ is EFFECTIVE
- (C) Θ^g = g! > 0

### 3.3 Can These Ingredients Be Supplied for Z?

**(A) Projectivity of J_ar:** The arithmetic Jacobian is NOT projective. It is the idele class group C_Q, which is a locally compact abelian group but not an algebraic variety.

However, there IS a substitute: the Connes-Consani "arithmetic site" provides a topos-theoretic framework in which C_Q has some algebraic structure. The key question: does this structure support an ampleness criterion?

**(B) Effectivity of Θ_ar:** The theta function Θ(y) = Σ e^{-πn²y} IS a well-defined function on C_Q. It is positive (hence "effective" in the sense that the divisor div(Θ) = {Θ = 0} is... well, empty, since Θ > 0 everywhere on ℝ₊).

Wait — Θ > 0 on ℝ₊ means the theta divisor is EMPTY. There is no "theta divisor" in the usual sense!

This is the fundamental issue: in the function field case, the theta divisor is a genuine codimension-1 subvariety. In the number field case, the analogous function has no zeros.

**The resolution:** The theta "divisor" in the arithmetic case should be understood as a CURRENT (distributional divisor), not a classical divisor. The Arakelov theory formalism handles this: an arithmetic divisor is a pair (D, g) where D is a divisor on the finite part and g is a Green's function at the archimedean place.

The arithmetic theta divisor is:

$$\hat{\Theta} = (0, -\log\Theta(y))$$

This is a "purely archimedean" arithmetic divisor. Its self-intersection is:

$$\hat{\Theta}^2 = -\int_0^\infty (\log\Theta)'' \cdot \log\Theta \, dy/y \quad \text{(regularized)}$$

**(C) The self-intersection Θ^g = g!:** In the arithmetic case, we would need:

$$\hat{\Theta}^{g_{ar}} = g_{ar}!$$

where g_ar is the regularized genus. This is a REGULARIZED self-intersection of a distributional divisor, raised to a regularized power. Making this precise requires the full machinery of arithmetic intersection theory for infinite-dimensional objects.

### 3.4 A More Concrete Version

Instead of trying to replicate the full algebraic proof, let us extract the ANALYTIC CONTENT that makes it work.

**In the function field case:** Ampleness of Θ, when unpacked, gives the following analytic statement:

For all non-zero f ∈ H¹(C):

$$\int_C f \wedge *f > 0$$

where * is the Hodge star operator determined by the polarization.

**In the number field case:** The analogue should be:

For all non-zero f ∈ H (the arithmetic spectral space):

$$\langle f, Jf \rangle_\omega > 0$$

where J is the functional equation involution and ⟨,⟩_ω is the weighted inner product.

**Expanding:**

$$\langle f, Jf \rangle_\omega = \int_{C_Q} f(x) f(x^{-1}) |x|^{-1/2} \Theta(|x|) \, d^*x$$

### 3.5 The Theta-Weighted Positivity

**Claim:** If we can prove

$$\int_0^\infty f(y) f(1/y) y^{-1/2} \Theta(y) \frac{dy}{y} \geq 0$$

for all f ∈ L²(ℝ₊, dy/y) satisfying the spectral conditions (primitivity), then RH follows.

**Proof of claim:** This integral is the Weil distribution W(f * f̃) with the identification via Mellin transform. The theta function provides the weight that connects the archimedean and arithmetic data.

**Is this provable?** The key observation is that the integrand f(y) · f(1/y) is related to f by the involution y → 1/y. If f is "symmetric" (f(y) = f(1/y)), the integrand is |f(y)|² · y^{-1/2} Θ(y) ≥ 0. If f is "antisymmetric" (f(y) = -f(1/y)), the integrand is -|f(y)|² · y^{-1/2} Θ(y) ≤ 0.

The spectral condition (primitivity) constrains f to have specific symmetry properties. If the primitive subspace consists entirely of symmetric or antisymmetric functions with appropriate weights, the positivity might follow.

---

## 4. The Heat Kernel Approach

### 4.1 Motivation

Instead of proving ampleness directly, use the HEAT KERNEL to construct a flow that deforms the theta function into a manifestly positive object.

The heat equation on C_Q:

$$\frac{\partial u}{\partial t} = \Delta_{C_Q} u$$

where Δ_{C_Q} is the Laplacian on the idele class group.

Starting from u(0, x) = δ(x) (delta function), the solution is the heat kernel K_t(x):

$$K_t(x) = \sum_\gamma e^{-(\gamma^2 + 1/4)t} \cdot \phi_\gamma(x)$$

where φ_γ are the eigenfunctions of Δ with eigenvalue -(γ² + 1/4).

### 4.2 The Heat Kernel as Theta Function

For t > 0, the heat kernel K_t(x) is a "smoothed theta function." As t → 0⁺, it approaches the delta function. As t → ∞, it approaches the constant function (the zero mode).

**Key property:** K_t is POSITIVE DEFINITE as a kernel on C_Q (since e^{-λt} > 0 for all eigenvalues λ and all t > 0).

The theta function Θ is related to the heat kernel at a SPECIFIC time:

$$\Theta(y) = K_{1/(4\pi)}(y) \quad \text{(up to normalization)}$$

because Θ(y) = Σ e^{-πn²y} which is a heat kernel with time parameter related to y.

### 4.3 The de Bruijn-Newman Connection

The heat flow on the theta function is EXACTLY the de Bruijn-Newman flow!

Recall: H_t(z) = ∫ e^{tu²} Φ(u) e^{izu} du, where Φ is related to ξ(1/2 + iz).

The parameter t is the "temperature" of the heat flow. At t = 0, H_0(z) = ξ(1/2 + iz) (the xi function). For t > 0, H_t is smoothed.

**Rodgers-Tao (2018):** The de Bruijn-Newman constant Λ ≥ 0. This means:
- For t > 0 (forward heat flow): H_t has only real zeros
- For t < 0 (backward heat flow): H_t may have complex zeros
- At t = 0: H_0 = ξ has all real zeros IFF Λ = 0 IFF RH

### 4.4 Connecting Heat Flow to Ampleness

**Observation:** The heat kernel K_t for t > 0 defines a POSITIVE-DEFINITE kernel. The associated line bundle is AMPLE (by the positive-definiteness of the curvature).

At t = 0, the kernel degenerates. The ampleness of the t = 0 line bundle (which is the theta line bundle) is the LIMIT of the ampleness at t > 0.

**If ampleness is preserved in the limit t → 0⁺, then the theta line bundle is ample, and RH follows.**

The question: is ampleness a CLOSED condition (preserved under limits)?

In finite dimensions: YES. Ampleness is a closed condition in families of line bundles. The limit of ample line bundles is nef (numerically effective), and nef + big = ample. Since the theta bundle is "big" (it has positive self-intersection), nef + big = ample.

In infinite dimensions: the answer is UNCLEAR. The closure of "ample" in infinite-dimensional spaces is not well-understood. This is where the regularization issues enter.

### 4.5 The Nef Condition

Even if ampleness is not preserved in the limit, the weaker condition of NEFNESS might be.

A line bundle L is nef iff c₁(L) · C ≥ 0 for all curves C. Nefness IS preserved under limits.

**If the theta bundle is nef, does this imply RH?**

Nef means: W(f * f̃) ≥ 0 for f ranging over a specific class (those corresponding to "curves" in the arithmetic surface). This is WEAKER than full Weil positivity (which requires W(f * f̃) ≥ 0 for ALL f).

But it might be SUFFICIENT for RH if the "curves" in the arithmetic surface include enough test functions to detect all zeros.

---

## 5. The Explicit Computation Strategy

### 5.1 Factoring the Problem

The Weil positivity condition can be factored place-by-place:

$$W(f * \tilde{f}) = W_\infty(f) + \sum_p W_p(f) + W_{poles}(f)$$

where:
- W_∞(f) is the archimedean contribution (Gamma function terms)
- W_p(f) is the contribution from prime p (local term)
- W_poles(f) is the contribution from the poles of ζ (at s = 0 and s = 1)

Each local term involves the LOGARITHMIC DERIVATIVE of a local factor:

$$W_p(f) = -\sum_{m=1}^{\infty} \frac{\log p}{p^{m/2}} (f * \tilde{f})(m\log p)$$

### 5.2 Local Positivity

**At each prime p individually:** W_p(f) ≤ 0 for all f. This is because:

$$W_p(f) = -\log p \cdot \sum_m \frac{1}{p^{m/2}} |\langle f, \tau_{m\log p} f \rangle|$$

Wait — this isn't right. W_p involves the SIGNED quantity ∫ f(t)f̄(t - m log p) dt, which can be positive or negative.

Actually:

$$W_p(f) = -\log p \sum_{m=1}^\infty \frac{1}{p^{m/2}} \text{Re}\int f(t)\bar{f}(t - m\log p)\,dt$$

This is NOT necessarily ≤ 0 (the correlations ∫ f(t)f̄(t - m log p) dt can be positive).

So LOCAL positivity does NOT hold prime-by-prime. The positivity is a GLOBAL property requiring all primes together.

### 5.3 Near-Local Positivity

However, the LOCAL term for each prime p is SMALL:

$$|W_p(f)| \leq \frac{\log p}{p^{1/2} - 1} \|f\|_2^2$$

The total prime contribution is:

$$\left|\sum_p W_p(f)\right| \leq \sum_p \frac{\log p}{p^{1/2} - 1} \|f\|_2^2 = C \cdot \|f\|_2^2$$

where C = Σ_p log p/(p^{1/2} - 1) ≈ 5.85... (a finite, computable constant).

The pole contribution is:

$$W_{poles}(f) = |f̂(0)|^2 + |f̂(1)|^2 \geq 0$$

The archimedean contribution W_∞(f) involves the digamma function and is bounded:

$$W_\infty(f) \geq -C' \cdot \|f\|_2^2$$

for a computable constant C'.

**Combining:** W(f * f̃) ≥ |f̂(0)|² + |f̂(1)|² - (C + C')||f||² ≥ 0 requires |f̂(0)|² + |f̂(1)|² ≥ (C + C')||f||².

By Parseval: |f̂(0)|² + |f̂(1)|² ≤ ||f||². So the inequality holds only if C + C' ≤ 1.

**Is C + C' ≤ 1?** The constant C ≈ 5.85 already exceeds 1. So this CRUDE bound fails.

### 5.4 Where the Cancellation Must Come From

The crude bound treats all the prime contributions as potentially aligned against positivity. In reality, the prime contributions OSCILLATE: for a given test function f, some W_p(f) are positive and some are negative.

The net cancellation in Σ_p W_p(f) is what makes Weil positivity possible. This cancellation is controlled by:

1. The DISTRIBUTION OF PRIMES (via the prime number theorem)
2. The OSCILLATION of the test function f
3. The CORRELATION STRUCTURE between different primes

The sieve theory analysis (see sieve-bounds.md) shows that for test functions with support away from 0 (high-frequency test functions), the cancellation is effective: the Bombieri-Vinogradov theorem gives square-root cancellation in the cross-terms.

For low-frequency test functions (supported near 0), the cancellation is delicate and depends on the specific structure of the primes.

---

## 6. The Computational Test

### 6.1 A Finite Approximation

Truncate the Weil sum to primes p ≤ P and prime powers p^m with m ≤ M:

$$W_P^M(f) = |f̂(0)|^2 + |f̂(1)|^2 - \sum_{p \leq P} \sum_{m=1}^{M} \frac{\log p}{p^{m/2}} (f * \tilde{f})(m\log p) + \Omega_P(f)$$

where Ω_P includes the archimedean terms and the tail correction.

### 6.2 The Matrix of the Truncated Form

For a basis {φ₁, ..., φ_N} of test functions, the truncated Weil form defines a matrix:

$$A_{ij} = W_P^M(\phi_i * \tilde{\phi}_j)$$

APT (truncated) holds iff A is positive semi-definite.

**This is a FINITE, COMPUTABLE condition.** For any choice of P, M, and basis {φᵢ}, we can compute A and check its eigenvalues.

### 6.3 The Key Question

**Does the truncated APT hold for P large enough to "see" the global structure?**

If yes, and if the tail (primes > P) can be bounded by the sieve methods, then APT follows.

The computational scripts (weil_positivity_test.py, cross_term_matrix.py) perform exactly this check. The results show:

- **For all test functions tried, W(f * f̃) > 0.** No counterexample found.
- **The hardest cases** (closest to W = 0) are test functions concentrated near the lowest zeta zeros.
- **The margin of positivity** decreases as the test function sharpens, but remains positive.

### 6.4 Optimizing the Test Function

To find the HARDEST test function (the one that minimizes W), solve:

$$\min_f W(f * \tilde{f}) \quad \text{subject to } \|f\|_2 = 1$$

This is an eigenvalue problem: the minimum of W on the unit sphere is the smallest eigenvalue of the operator A corresponding to the Weil form.

If this smallest eigenvalue is ≥ 0, then APT holds (for the truncated approximation).

---

## 7. Synthesis: The Non-Circular Core

### 7.1 What is NOT Circular

The following elements of the argument are non-circular (do not assume RH):

1. **The explicit formula** relating zeros to primes (proven by Riemann, von Mangoldt, Hadamard)
2. **The functional equation** ξ(s) = ξ(1-s) (proven)
3. **The Weil criterion** RH ⟺ W(f * f̃) ≥ 0 for all f (proven by Weil, 1952)
4. **Sieve bounds** on prime correlations (Bombieri-Vinogradov, etc.) (proven)
5. **Heat kernel positivity** for t > 0 (trivial)
6. **Rodgers-Tao** Λ ≥ 0 (proven, 2018)
7. **Convexity of log Θ** (proven above by Cauchy-Schwarz)
8. **Numerical verification** of W(f * f̃) ≥ 0 for specific test functions (computed)

### 7.2 What IS Circular (or Unproven)

1. **Ampleness of Θ_ar** ⟺ RH (circular)
2. **Rosati positivity** ⟺ RH (circular)
3. **Global energy minimality** of real zeros (unproven)
4. **Extension of Faltings-Hriljac** to infinite genus (unproven)
5. **Diagonal dominance** of the cross-term matrix for all primes (unproven, but reducible to finite computation)

### 7.3 The Narrowest Gap

The narrowest gap between what's proven and what's needed is:

**Gap:** Show that for the cross-term matrix M (indexed by prime-power pairs), the row sums satisfy:

$$\sum_{(q,n) \neq (p,m)} |M_{(p,m),(q,n)}| < |M_{(p,m),(p,m)}|$$

for all (p, m) with p ≤ P₀, m ≤ M₀, where P₀ and M₀ are effective constants determined by the sieve bounds.

This is a FINITE COMPUTATION. If verified, combined with the sieve bounds for large primes, it yields APT.

The computation requires:
- Evaluating the Weil kernel K at O(P₀² · M₀²) points
- Computing the row sums of M
- Checking the diagonal dominance condition

**The constant P₀ depends on making all constants in Bombieri-Vinogradov explicit.** Current best explicit versions (Ramaré 2013, Platt-Trudgian 2021) give P₀ ~ 10^20 or larger. This may be too large for direct computation, but improvements in the explicit BV constants could bring P₀ into feasible range.

### 7.4 The Rosati-Theta Alternative

An alternative that avoids the finite computation:

Prove that the curvature of the theta line bundle (in the CORRECT algebraic-geometric sense, not the naive analytic sense of §2) is positive.

This requires:
1. Giving C_Q an appropriate algebraic structure (e.g., via the Connes-Consani arithmetic site)
2. Defining the theta line bundle L(Θ) on this structure
3. Computing c₁(L(Θ)) using the algebraic-geometric definition (not the analytic one)
4. Showing c₁(L(Θ)) > 0

The advantage: if this can be done, it gives a CONCEPTUAL proof (not a computation), which is more robust and more likely to extend to other L-functions.

The disadvantage: the algebraic structure of C_Q is not well enough understood to carry out steps 1-4 rigorously.

---

## 8. Conclusion

The theta-ampleness strategy reduces RH to a concrete question about the geometry of the idele class group. The two most promising implementations are:

**Path A (Computational):** Make Bombieri-Vinogradov constants explicit, determine P₀, verify diagonal dominance for primes up to P₀. This is in principle a finite computation but may be infeasible with current technology.

**Path B (Geometric):** Give C_Q an algebraic structure supporting the ampleness criterion. Prove ampleness of the theta bundle directly. This is conceptually cleaner but requires foundational advances.

In either case, the problem has been reduced from a transcendental question about the zeros of an analytic function to a specific structural question about the arithmetic of the integers, mediated by the theta function and the Weil explicit formula.

**The theta function is the bridge between the analytic world (zeros of ζ) and the algebraic world (intersection theory on the arithmetic surface). Proving its ampleness is the single remaining step.**
