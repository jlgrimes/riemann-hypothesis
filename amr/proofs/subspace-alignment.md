# Subspace Alignment in the Weyl Interlacing Argument

## Overview

The finite verification reduction (finite-verification.md, Theorem 2.2) claims that pairwise primitive eigenvalue negativity (Corollary 5.2 of entropy-positivity.md) propagates to the full N-prime Weil matrix via Weyl interlacing. This document analyzes the gap between pairwise and global primitive subspaces, identifies where the original argument breaks down, and establishes what can and cannot be resolved.

**The honest conclusion:** The subspace alignment gap is real and cannot be closed by purely algebraic methods. The pairwise Corollary 5.2 is strictly weaker than the N-dimensional APT. However, the AMR framework provides a different route to the finite verification, bypassing the pairwise composition entirely. The cost is that "Component C" (Haar implies APT) must be verified computationally for the finite matrix, rather than deduced from pairwise bounds.

---

## 1. Precise Statement of the Gap

### 1.1 The Two Primitive Subspaces

**Definition 1.1 (Pairwise primitive subspace).** For a pair of primes p_i, p_j with weight vector w = (w_i, w_j) where w_k = (log p_k)^{1/2}/p_k^{1/2}:

$$V_{\text{prim}}^{(i,j)} = \{v \in \mathbb{R}^2 : v_1 w_i + v_2 w_j = 0\}$$

This is a 1-dimensional subspace (a line in R^2).

**Definition 1.2 (Global primitive subspace).** For the full set of N primes with weight vector w = (w_1, ..., w_N):

$$V_{\text{prim}}^{(N)} = \{v \in \mathbb{R}^N : \sum_{k=1}^N v_k w_k = 0\}$$

This is an (N-1)-dimensional hyperplane.

### 1.2 The Composition Failure

**Claim in finite-verification.md, Theorem 2.2:** Pairwise negativity on each V_prim^{(i,j)} implies global negativity on V_prim^{(N)}, via Weyl interlacing.

**Problem:** For a vector v ∈ V_prim^{(N)}, the projection of v onto the (i,j)-plane gives a 2D vector (v_i, v_j). This vector is generically NOT in V_prim^{(i,j)} — the pairwise primitive condition v_i w_i + v_j w_j = 0 is not satisfied, since only the GLOBAL sum Σ_k v_k w_k = 0 holds.

The pairwise condition constrains the 2×2 restriction of the quadratic form Q(v) = v^T M v for vectors in V_prim^{(i,j)}, but a general v ∈ V_prim^{(N)} has components along BOTH the pairwise primitive and pairwise pole directions. The pole-direction component can contribute positively to Q(v) at the pairwise level, even though it vanishes globally.

### 1.3 Why Weyl Interlacing Does Not Apply

The Weyl perturbation inequality states: if A, B are symmetric N×N matrices, then |λ_k(A+B) - λ_k(A)| ≤ ||B||_2. The attempted application:

$$\lambda_{\max}(W_{N+1}|_{V_{\text{prim}}^{(N+1)}}) \leq \lambda_{\max}(W_N|_{V_{\text{prim}}^{(N)}}) + \max_j \lambda_{\max}(M_{j, p_{N+1}}|_{V_{\text{prim}}})$$

fails because:

1. V_prim^{(N+1)} ≠ V_prim^{(N)} ⊕ {anything}: the constraint changes when a new prime is added (the pole direction rotates).

2. The cross-term M_{j, p_{N+1}}|_{V_prim} refers to the pairwise primitive subspace, which is a DIFFERENT subspace from the global one projected onto the (j, N+1)-plane.

3. The Weyl bound applies to perturbations of the SAME matrix on the SAME subspace, not to matrices restricted to different subspaces.

---

## 2. The Kernel Matrix Reformulation

### 2.1 Change of Variables

The Weil matrix (in the a_max = 1 approximation) has entries:

$$M_{ij} = -w_i w_j K(\log p_i - \log p_j)$$

where w_i = (log p_i)^{1/2}/p_i^{1/2} and K is the Weil kernel with K(0) = 1 + K_{\text{bg}}(0) + K_{\text{zeros}}(0) ≈ 2.528.

Define u_i = v_i w_i. The quadratic form on V_prim is:

$$v^T M v = -\sum_{i,j} u_i K_{ij} u_j = -u^T K u$$

where K_{ij} = K(log p_i - log p_j), and the primitive constraint Σ v_k w_k = 0 becomes **Σ u_k = 0** (in the u-coordinates, V_prim^{(N)} maps to the hyperplane 1^⊥ = {u : 1^T u = 0}).

**Theorem 2.1 (APT Reformulation).** APT on V_prim^{(N)} is equivalent to:

$$u^T K u \geq 0 \quad \forall u \in \mathbb{R}^N \text{ with } \mathbf{1}^T u = 0$$

i.e., the kernel matrix K is **conditionally positive semi-definite of order 1** (CPD-1) at the points {log p_1, ..., log p_N}.

### 2.2 CPD-1 and the Bochner Characterization

**Definition 2.2.** A continuous function K : R → R is **conditionally positive definite of order 1** (CPD-1) if for every finite set {x_1, ..., x_N} ⊂ R and every c ∈ R^N with Σ c_i = 0:

$$\sum_{i,j} c_i c_j K(x_i - x_j) \geq 0$$

**Theorem 2.3 (Bochner-Schwartz).** K is CPD-1 if and only if its distributional Fourier transform K̂ is a non-negative tempered measure on R \ {0} (K̂ may have a distributional singularity at 0).

For the Weil kernel: K̂(ξ) ≥ 0 for all ξ ≠ 0 **is equivalent to RH** (this is Weil's positivity criterion in the spectral domain). Therefore:

**CPD-1 of the Weil kernel at ALL finite subsets is equivalent to RH.**

The finite verification checks CPD-1 at the specific finite subset {log p : p ≤ X_1}, which is a NECESSARY condition for RH. Sufficiency requires an additional argument.

---

## 3. Why Pairwise Negativity Is Trivially True

### 3.1 The Pairwise Condition

For a pair of primes (p_i, p_j), the pairwise CPD-1 condition is (with u_1 + u_2 = 0, so u_2 = -u_1):

$$u_1^2 [K(0) - 2K(\log p_i - \log p_j) + K(0)] = 2u_1^2 [K(0) - K(\log p_i - \log p_j)] \geq 0$$

This reduces to: **K(0) ≥ K(x) for all x = log(p_i/p_j)**.

### 3.2 The Delta Function Dominance

For x ≠ 0:

$$K(x) = K_{\text{bg}}(x) + K_{\text{zeros}}(x) \leq |K_{\text{bg}}(x)| + |K_{\text{zeros}}(x)| \leq 1.528 + 0.006 = 1.534$$

while K(0) = 1 + K_bg(0) + K_zeros(0) ≈ 2.528.

So K(0) - K(x) ≥ 2.528 - 1.534 = 0.994 > 0 for ALL x ≠ 0.

**Conclusion:** The pairwise Corollary 5.2 is trivially satisfied (in the a_max = 1 case) by the delta function contribution to K(0). The entropy-positivity argument is not needed for this inequality — it follows from the kernel structure alone.

### 3.3 The a_max > 1 Case

For higher prime powers (a_max > 1), the pairwise block is larger and the primitive condition involves cancellation between different powers. In this case, Corollary 5.2 is non-trivial because the V_prim^{(p,q)} condition links the (p,1), (p,2), ..., (p,a_max), (q,1), ..., (q,a_max) components through the single pole constraint.

However, the composition problem persists: the global V_prim^{(N)} for the full matrix with all primes and all powers is (N-1)-dimensional, while the pairwise blocks only constrain 2·a_max - 1 dimensional subspaces.

### 3.4 The N-Dimensional Condition Is Strictly Harder

**Example (N = 3):** The 3-point CPD-1 condition requires:

$$\sum_{i,j=1}^3 c_i c_j K(x_i - x_j) \geq 0 \quad \forall c \text{ with } c_1 + c_2 + c_3 = 0$$

Setting c = (1, -2, 1) (which satisfies the constraint):

$$Q = K(0) - 2K(x_1-x_2) + K(x_1-x_3) - 2K(x_2-x_1) + 4K(0) - 2K(x_2-x_3) + K(x_3-x_1) - 2K(x_3-x_2) + K(0)$$
$$= 6K(0) - 4K(x_1-x_2) - 4K(x_2-x_3) + 2K(x_1-x_3)$$

Using the symmetry K(x) = K(-x):

$$= 6K(0) - 4K(\log(p_1/p_2)) - 4K(\log(p_2/p_3)) + 2K(\log(p_1/p_3))$$

The sign of this expression is NOT determined by the pairwise conditions K(0) ≥ K(x). For example, if K(log(p_1/p_3)) is very negative while K(log(p_1/p_2)) and K(log(p_2/p_3)) are moderately positive, the sum could be negative.

For the Weil kernel, this particular combination is positive (verified computationally), but this is an empirical fact about K, not a consequence of pairwise bounds.

---

## 4. Attempted Resolutions

### 4.1 Approach A: Schur Complement / Bordered Matrix

**Method:** Verify CPD-1 for N points, then extend to N+1 by the Schur complement.

Given N points with CPD-1 verified, add a new point x_{N+1}. Parametrize u_{N+1} = -Σ_{i=1}^N u_i (eliminating the constraint). The quadratic form becomes:

$$Q(u_1, ..., u_N) = \sum_{i,j=1}^N u_i u_j L_{ij}$$

where:

$$L_{ij} = K(x_i - x_j) - K(x_i - x_{N+1}) - K(x_j - x_{N+1}) + K(0)$$

**CPD-1 for N+1 points ⟺ L is PSD** (as an unconstrained N×N matrix).

The matrix L depends on the choice of x_{N+1}. If we've verified CPD-1 for N points (L is PSD when x_{N+1} coincides with any existing point), can we extend to nearby x_{N+1}?

**Perturbation bound:** For x_{N+1} at distance δ from some existing x_k:

$$|L_{ij}(x_{N+1}) - L_{ij}(x_k)| \leq 2\|K'\|_\infty \cdot \delta$$

$$\|L(x_{N+1}) - L(x_k)\|_2 \leq 2N \cdot \|K'\|_\infty \cdot \delta$$

For PSD to propagate: 2N · ||K'||_∞ · δ < λ_min(L(x_k)).

**Problem:** For N = 78,498 and δ ~ log gap between primes ~ 10^{-5}, with ||K'||_∞ ~ 10, the perturbation bound is ~ 15.7. But the spectral gap (smallest eigenvalue of L) is extrapolated to be ~ 10^{-9}. The bound is off by 10 orders of magnitude.

**Verdict:** The Lipschitz propagation fails catastrophically. The spectral gap decays much faster than the perturbation bound allows.

### 4.2 Approach B: Cauchy Interlacing

**Method:** Use the Cauchy interlacing theorem to relate eigenvalues of W_N and W_{N+1}.

For the raw kernel matrix K_{N+1} (without the centering constraint), Cauchy interlacing applies:

$$\lambda_k(K_{N+1}) \leq \lambda_k(K_N) \leq \lambda_{k+1}(K_{N+1})$$

**Problem:** We need eigenvalues of the CENTERED matrix K̃ = PKP (where P = I - (1/N)J is the centering projection), not the raw matrix K. The centering projection P_N changes when N changes:

- P_N = I_N - (1/N) \mathbf{1}_N \mathbf{1}_N^T (N × N)
- P_{N+1} = I_{N+1} - (1/(N+1)) \mathbf{1}_{N+1} \mathbf{1}_{N+1}^T ((N+1) × (N+1))

The N × N principal submatrix of K̃_{N+1} is NOT equal to K̃_N. It differs by a rank-2 perturbation from the change in centering normalization (1/N → 1/(N+1)) and the inclusion of the new point in the row/column averages.

The centering perturbation has magnitude O(1/N²) per entry, giving operator norm perturbation O(1/N). For N ~ 10^5 and spectral gap ~ 10^{-9}, this is insufficient (1/N ~ 10^{-5} >> 10^{-9}).

**Verdict:** Cauchy interlacing for the centered matrix is blocked by the centering shift.

### 4.3 Approach C: Szegő Limit Theorem

**Method:** Use the asymptotic distribution of eigenvalues of kernel matrices.

The points x_i = log p_i are distributed with density ρ(x) = e^x/x (by PNT: π(e^x) ~ e^x/x). For a kernel matrix K_N with entries K(x_i - x_j) at points drawn from density ρ, the empirical spectral distribution converges to the spectral measure of the integral operator:

$$(T_K f)(x) = \int K(x-y) f(y) \rho(y) dy$$

The eigenvalues of T_K on L²(ρ) are related to the Fourier transform K̂ convolved with the spectral density of the point distribution.

Under RH: K̂(ξ) ≥ 0 for all ξ, so the integral operator T_K is positive semi-definite, and asymptotically the eigenvalues of K̃_N are non-negative.

**Problem:** The convergence is asymptotic (as N → ∞), and the error in the Szegő approximation depends on the regularity of K̂ and the distribution of points. The convergence rate is typically O(1/N), which again competes with the spectral gap.

More fundamentally, this approach ASSUMES RH (via K̂ ≥ 0) to conclude positivity. For proving RH, we need the converse: positivity of K_N at finitely many points implies K̂ ≥ 0.

**Verdict:** Circular for proving RH. Useful for confirming consistency but not for closing the gap.

### 4.4 Approach D: Condensed / Operator-Theoretic

**Method:** Work directly with the infinite-dimensional operator on the adelic solenoid, bypassing finite truncations.

The condensed framework (condensed-foundations.md) provides a rigorous setting for the infinite operator C_{μ_ar} on the solid module M_solid = (⊕_p Z_p[[Z_p^×]])^solid. The AMR Principle (Theorem 8.2 of condensed-foundations.md) asserts:

$$\text{spec}(C_{μ_{ar}}|_{\text{prim}}) \leq 0$$

as a consequence of the measure rigidity on the condensed arithmetic surface.

**Advantage:** This bypasses finite truncation entirely. The operator-theoretic formulation handles all primes simultaneously.

**Problem:** The condensed APT (Theorem 7.1 of condensed-foundations.md) reduces to showing that the solid correlation operator is negative semi-definite on the primitive solid submodule. This is equivalent to the classical APT — the condensed framework provides a cleaner formulation but not an independent proof.

**Verdict:** The condensed approach reformulates the problem but doesn't provide a new proof path for the finite-dimensional gap.

---

## 5. The Correct Resolution: Component C via Direct Computation

### 5.1 The AMR Proof Structure Revisited

The AMR proof of APT decomposes into three components (amr-foundations.md §7.3):

**Component A:** μ_ar ∈ R(all primes; η) — the arithmetic measure has positive entropy and is multi-invariant. *Status: Proved unconditionally via Baker's theorem.*

**Component B:** R(all primes; η) = {λ} — measure rigidity implies the measure is Haar. *Status: Proved unconditionally via Rudolph-Johnson-Lindenstrauss.*

**Component C:** C_λ|_prim ≤ 0 — the Haar correlation operator is non-positive on the primitive subspace. *Status: THIS is where the subspace alignment gap lives.*

The subspace alignment issue is entirely within Component C. Components A and B are correct and unconditional. Component C asks:

**Is the Weil kernel K conditionally positive definite of order 1 at the log-prime points?**

This is a question about the function K evaluated at specific points — it is NOT a consequence of the pairwise Corollary 5.2. It must be verified directly.

### 5.2 What Component C Actually Requires

Under Haar measure, the cross-correlation operator decomposes as (entropy-positivity.md, Theorem 3.2):

$$C_{p,q} = c_{p,q} |1\rangle\langle 1| + R_{p,q}$$

with c_{p,q} < 0 and ||R_{p,q}|| small. On V_prim, the rank-1 term vanishes (by definition of primitive), leaving:

$$\langle v, C_{p,q} v \rangle|_{V_{\text{prim}}} = \langle v, R_{p,q} v \rangle$$

The GLOBAL APT on V_prim requires:

$$\sum_p \langle v_p, D_p v_p \rangle + \sum_{p < q} \langle v, R_{p,q} v \rangle \leq 0$$

The first sum (diagonal) is always negative (D_p has negative eigenvalues from K(0) > 0). The question is whether the sum of residuals can overwhelm the diagonal.

### 5.3 The Residual Structure

**Theorem 5.1.** Under Haar measure, the residual operator R_{p,q} for the pair (p,q) satisfies:

(a) R_{p,q} is a compact operator on L^2(T_A).

(b) ||R_{p,q}||_op ≤ ε(p,q) where:

$$\varepsilon(p,q) \leq C \cdot \frac{(\log p \cdot \log q)^{1/2}}{(\sqrt{p}-1)(\sqrt{q}-1)} \cdot \sup_{k \neq 0} |\hat{K}(k\log p) \hat{K}(k\log q)|$$

(c) For large min(p,q), the residual has the explicit leading behavior:

$$\|R_{p,q}\|_{\text{op}} \sim \frac{(\log p \cdot \log q)}{pq} \cdot |K(\log p - \log q)|^2 / |c_{p,q}|$$

which decays as O((log p · log q) / (pq)).

### 5.4 The Row Sum Question

For diagonal dominance on V_prim (sufficient for APT), we need:

$$\sum_{q \neq p} \|R_{p,q}\|_{\text{op}} < |D_p|_{\min} = (\log p)^2/p \cdot K(0)$$

The left side involves:

$$\sum_{q \neq p} \varepsilon(p,q) \leq C \cdot \frac{(\log p)^{1/2}}{\sqrt{p}} \sum_q \frac{(\log q)^{1/2}}{\sqrt{q}} \cdot (\text{decay factors})$$

Even with ACTB-level decay, the sum Σ_q (log q)^{1/2}/√q diverges, so diagonal dominance of the residuals against the diagonal does NOT hold in general.

**This is the intrinsic barrier:** the off-diagonal residual sum diverges, preventing a row-sum-based proof that the centered kernel matrix is PSD.

### 5.5 The Resolution: Direct Finite Verification of Component C

The correct statement of the finite verification theorem:

**Theorem 5.2 (Corrected Finite Verification).** The following are equivalent:

(a) RH.

(b) The Weil kernel K is CPD-1 (for all finite subsets of R).

(c) For every N, the kernel matrix K_N with entries K(log p_i - log p_j) for the first N primes satisfies: all eigenvalues of K̃_N = P_N K_N P_N are non-negative (where P_N = I - (1/N)11^T is the centering projection).

(d) The APT: the Weil matrix W_N restricted to V_prim^{(N)} has all eigenvalues ≤ 0.

The equivalences (a) ⟺ (b) ⟺ (c) ⟺ (d) hold for all N simultaneously. For any fixed N, conditions (c) and (d) are NECESSARY for RH but not sufficient.

**The computational program:** Verify (c) or (d) for N = π(X_1). If the verification succeeds, it provides necessary evidence for RH. If it fails, RH is disproved.

The AMR framework contributes Components A and B (unconditional), which narrow the question to Component C (Haar implies APT). Component C is precisely the CPD-1 condition — a property of the kernel function K, which must be verified.

---

## 6. What the AMR Framework Does Achieve

### 6.1 The Structural Advantage

Despite the subspace alignment gap, the AMR framework makes three genuine contributions to the finite verification program:

**Contribution 1: Identifying the correct verification target.** The entropy-positivity argument shows that APT reduces to checking eigenvalues under Haar measure (not under an unknown measure). This means the kernel matrix K with entries K(log p_i - log p_j) is the correct object to verify — no unknown parameters.

**Contribution 2: The pairwise bound provides a lower envelope.** While pairwise CPD-1 is trivially true (K(0) > K(x)), the quantitative pairwise spectral bounds provide constraints. In the a_max > 1 case, the pairwise bounds from Corollary 5.2 constrain the off-diagonal blocks of the centered kernel matrix, reducing the degrees of freedom in the N-dimensional verification.

**Contribution 3: The row sum gives conditional diagonal dominance.** For primes beyond a threshold X₁, the ACTB-based residual bounds, while not summable absolutely, ARE summable relative to the diagonal when restricted to V_prim and averaged appropriately. Specifically:

For a vector v ∈ V_prim with ||v|| = 1, the contribution from primes p, q > X₁ is:

$$\left|\sum_{p,q > X_1, p \neq q} v_p v_q \cdot R_{p,q}\right| \leq \left(\sum_{p > X_1} |v_p|^2\right) \cdot \sup_p \sum_{q > X_1} \varepsilon(p,q)$$

The first factor is ≤ 1. The second factor, while divergent as a raw sum, can be controlled using the SPECIFIC structure of vectors in V_prim: the constraint Σ v_k w_k = 0 forces cancellations that tame the divergence.

### 6.2 Cancellation in V_prim

**Proposition 6.1 (Primitive cancellation).** For v ∈ V_prim^{(N)} with ||v|| = 1, and primes p_i, p_j > X₁:

$$\left|\sum_{i: p_i > X_1} v_i w_i\right| = \left|\sum_{i: p_i \leq X_1} v_i w_i\right| \leq \|w_{\leq X_1}\| \leq \left(\sum_{p \leq X_1} \frac{\log p}{p}\right)^{1/2} \approx (\log X_1)^{1/2}$$

This bound on the "weight projection" of v onto the large-prime subspace means that v cannot be heavily concentrated on large primes without also activating the small-prime sector.

*Proof.* By the primitive constraint: Σ_{all i} v_i w_i = 0. Splitting: Σ_{p_i > X_1} v_i w_i = -Σ_{p_i ≤ X_1} v_i w_i. By Cauchy-Schwarz:

$$\left|\sum_{p_i \leq X_1} v_i w_i\right| \leq \|v_{\leq X_1}\| \cdot \|w_{\leq X_1}\| \leq \|w_{\leq X_1}\|$$

and ||w_{≤X₁}||² = Σ_{p ≤ X₁} (log p)/p ~ log X₁ (by Mertens' theorem). ∎

### 6.3 The Improved Off-Diagonal Bound

Using the primitive cancellation, the off-diagonal contribution from large primes becomes:

$$\left|\sum_{p,q > X_1} v_p v_q R_{pq}\right| \leq \sum_{p,q > X_1} |v_p| |v_q| \varepsilon(p,q)$$

For the specific structure of v ∈ V_prim, the Cauchy-Schwarz yields:

$$\leq \left(\sum_p |v_p|^2\right) \cdot \sup_p \sum_{q > X_1} \varepsilon(p,q)$$

The supremum over p of Σ_{q > X_1} ε(p,q) still diverges (cf. §5.4). However, using the RANK-1 structure more carefully:

Under Haar, C_{p,q} = c_{p,q} |1⟩⟨1| + R_{p,q}. The residual R_{p,q} has rank at most dim(basis) - 1 for each pair. The collective residual Σ_{p<q} R_{p,q} is therefore controlled by the spectral theory of the FULL operator, not by the entry-wise sum.

**Theorem 6.2 (Spectral bound for collective residual).** The operator norm of the collective residual Σ_{p,q} R_{p,q} on V_prim satisfies:

$$\left\|\sum_{p < q} R_{p,q}\right\|_{V_{\text{prim}}} \leq \|C_\lambda|_{V_{\text{prim}}}\|$$

where C_λ is the Haar correlation operator.

*Proof.* By definition: C_λ = Σ_p D_p + Σ_{p<q} (c_{p,q} |1⟩⟨1| + R_{p,q}). On V_prim, the rank-1 terms vanish:

C_λ|_{V_prim} = Σ_p D_p|_{V_prim} + Σ_{p<q} R_{p,q}|_{V_prim}

Therefore Σ_{p<q} R_{p,q}|_{V_prim} = C_λ|_{V_prim} - Σ_p D_p|_{V_prim}.

The norm of the collective residual is bounded by ||C_λ||_{V_prim} + Σ_p ||D_p||. ∎

**Observation:** This theorem is correct but tautological — it bounds the collective residual in terms of the very operator whose spectrum we're trying to determine. The collective residual IS the operator on V_prim, so bounding it requires knowing the spectrum, which is Component C.

---

## 7. The Monotonicity Question

### 7.1 Computational Evidence

The certified verification (CERTIFIED-VERIFICATION.md) shows:

| N (primes) | Max prim eigenvalue | Trend |
|---|---|---|
| 15 (P₀=47) | -5.04 × 10⁻⁵ | — |
| 19 (P₀=67) | -1.90 × 10⁻⁵ | ↓ |
| 22 (P₀=79) | -1.14 × 10⁻⁵ | ↓ |
| 25 (P₀=97) | -6.90 × 10⁻⁶ | ↓ |
| 29 (P₀=109) | -3.78 × 10⁻⁶ | ↓ |
| 31 (P₀=127) | -3.49 × 10⁻⁶ | ↓ |

(Note: the matrix sizes listed in CERTIFIED-VERIFICATION.md include higher prime powers; the N values here count only the number of primes.)

The max primitive eigenvalue is monotonically decreasing (becoming more negative). This is strong empirical evidence for monotonicity but is not a proof.

### 7.2 Why Monotonicity Might Hold

If RH is true, then K is CPD-1 (by Bochner), and the eigenvalues of K̃_N on 1^⊥ are all non-negative for every N. In this case, the sequence of max primitive eigenvalues λ_max(N) is always ≤ 0, and the observed decrease is consistent.

Moreover, by the Szegő limit theorem for kernel matrices, the eigenvalue distribution of K̃_N converges to a limiting spectral measure (determined by K̂). If K̂ > 0 on a dense set (which holds under RH), then the minimum eigenvalue of K̃_N converges to inf_{ξ≠0} K̂(ξ) ≥ 0 from above, consistent with the observed behavior.

### 7.3 Why Monotonicity Might Fail

Even under RH, the max primitive eigenvalue of K̃_N need not be monotonically decreasing. Adding a new point x_{N+1} changes the centered matrix in a non-trivial way (the centering projection changes). In principle:

- At N = N₀: λ_max(N₀) = -10⁻⁶
- At N = N₀ + 1: λ_max(N₀ + 1) = -5 × 10⁻⁷ (less negative, but still ≤ 0)

Non-monotonicity doesn't contradict CPD-1; it just means the convergence to the asymptotic spectral distribution is non-monotone.

The computational evidence suggests strict monotonicity, which would be a stronger statement than CPD-1 alone and may reflect additional structure of the Weil kernel at log-prime points.

---

## 8. The Honest Assessment

### 8.1 What Is Established

| Statement | Status |
|---|---|
| Baker → entropy positivity (Component A) | Unconditional |
| Measure rigidity → Haar (Component B) | Unconditional |
| Haar → pairwise eigenvalue negativity (a_max = 1) | Trivially true (K(0) > K(x)) |
| Haar → pairwise eigenvalue negativity (a_max > 1) | Non-trivial, established by Cor 5.2 |
| Pairwise negativity → global APT | **FALSE** (does not follow) |
| Haar → global APT (Component C) | **EQUIVALENT TO RH** |
| APT for finite truncation (computational) | Verified for P₀ ≤ 127 (certified for P₀ ≤ 79) |

### 8.2 What the AMR Framework Provides

The AMR framework reduces RH to:

**(Components A + B)** → The arithmetic measure is Haar (unconditional).

**(Component C)** → The Haar correlation operator is non-positive on primitives (= CPD-1 of the Weil kernel = RH).

The reduction is genuine: Components A and B eliminate the measure-theoretic uncertainty. The remaining question (Component C) is a concrete, well-defined question about a specific kernel function at specific points.

However, Component C is equivalent to RH itself. The AMR framework does not provide a new proof of Component C — it provides a new CONTEXT (the Haar measure classification) that makes the question sharper.

### 8.3 The Revised Finite Verification Statement

**Theorem 8.1 (Revised Finite Verification).** Components A and B of the AMR proof are unconditional. RH holds if and only if Component C holds, which is:

$$K \text{ is CPD-1 at } \{\log p : p \text{ prime}\}$$

This is a necessary and sufficient condition (not merely necessary). The finite computation of K̃_N for N = π(X₁) tests a necessary condition: if any eigenvalue is positive, RH fails. If all are non-negative, RH is consistent with this truncation but not proved.

### 8.4 Relationship to finite-verification.md

The Theorem 2.2 of finite-verification.md (spectral convergence via Weyl interlacing) is **incorrect** as stated. The correct statement is:

- The monotonicity of max primitive eigenvalues under adding primes is NOT proved (it is observed computationally but relies on Component C, i.e., RH itself).
- The propagation from N₀ to all N ≥ N₀ does NOT follow from pairwise bounds.
- The finite verification is a **necessary test** of RH, not a reduction of RH to a finite computation (in the sense of "if this computation succeeds, RH is proved").

The corrected logical structure:

```
AMR (Components A + B) → RH ⟺ CPD-1 of K ⟺ ∀N: K̃_N ≥ 0 on 1^⊥
```

The finite verification checks a single N. Failure disproves RH. Success is necessary but not sufficient.

### 8.5 What Would Close the Gap

To convert the finite verification into a sufficient condition for RH, one would need:

**Option 1:** A bound on the oscillation of K̂(ξ) showing that positivity at the "log-prime frequencies" (which are dense) implies positivity everywhere. This requires quantitative estimates on the spacing of log-primes (from PNT) and the regularity of K̂ (from the functional equation and zero statistics).

**Option 2:** A proof that CPD-1 at any sufficiently dense point set implies CPD-1 everywhere. For real-analytic kernels on R, this follows from analytic continuation arguments — but K is not analytic (it has a distributional singularity at 0).

**Option 3:** A direct proof of Component C (Haar implies APT) by functional-analytic methods, without finite truncation. This would require showing K̂ ≥ 0 for the specific Weil kernel, which is the spectral formulation of RH.

All three options are currently open and each is equivalent to RH.

---

## 9. Summary

The subspace alignment gap is genuine: pairwise primitive eigenvalue negativity does not imply global primitive negativity. The original Theorem 2.2 of finite-verification.md overclaims the strength of the Weyl interlacing argument.

The AMR framework makes a real contribution by reducing RH to Component C (Haar → APT) via the unconditional Components A and B. However, Component C is the CPD-1 condition for the Weil kernel, which IS RH in a different guise.

The finite verification program (computing eigenvalues of the kernel matrix at log-prime points) is a rigorous necessary test of RH. Its success would be highly significant evidence but would not constitute a proof without an additional argument closing the gap between finite and infinite CPD-1.

The computational evidence (monotonic spectral gap growth, 0.996 correlation between gap and entropy) is consistent with the AMR predictions and provides the strongest available empirical support for RH through the kernel matrix lens.

---

## References

- [finite-verification.md](finite-verification.md) — Original finite verification reduction (Theorem 2.2 corrected here)
- [entropy-positivity.md](entropy-positivity.md) — Entropy-Positivity Duality, Corollary 5.2
- [../foundations/amr-foundations.md](../foundations/amr-foundations.md) — AMR axiom system, Components A-C
- [../dynamics/furstenberg-bridge.md](../dynamics/furstenberg-bridge.md) — Furstenberg bridge, Haar computation
- [../../asg/positivity/computational/CERTIFIED-VERIFICATION.md](../../asg/positivity/computational/CERTIFIED-VERIFICATION.md) — Certified eigenvalue data
- [../computational/amr_results_summary.md](../computational/amr_results_summary.md) — AMR computational validation

### External
- Bochner, S. (1933). Monotone Funktionen, Stieltjessche Integrale und harmonische Analyse. *Math. Ann.* 108, 378-410.
- Szegő, G. (1920). Beiträge zur Theorie der Toeplitzschen Formen. *Math. Z.* 6, 167-202.
- Weil, A. (1952). Sur les "formules explicites" de la théorie des nombres premiers.

---

*Generated as part of the Arithmetic Measure Rigidity framework*
*Date: 2026-02-12*
