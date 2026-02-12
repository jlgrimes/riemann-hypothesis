# Regularity Structures and Singular SPDEs for the Weil Kernel

## Status: Regularity structures provide a rigorous framework for classifying the distributional singularity of the Weil kernel and suggest deep algebraic connections (via Bruned-Hairer-Zambotti --> Connes-Kreimer --> motivic zeta values) but do NOT provide a path to proving CPD-1. The obstruction is categorical: regularity structures handle LOCAL singularity analysis, while CPD-1 is a GLOBAL positivity condition. The renormalization group acts on models (local descriptions), not on spectral measures (global objects).

---

## 0. Overview

The Weil kernel K(x) = delta(x) + K_bg(x) + K_zeros(x) has a distributional singularity at x = 0 (the Dirac delta). Hairer's regularity structures (2014 Fields Medal) provide the most powerful existing framework for handling distributional singularities in the context of stochastic PDEs. This document investigates whether that framework — or its algebraic extensions via Bruned-Hairer-Zambotti (2019) — can illuminate the CPD-1 property of K, which is equivalent to RH (bochner-proof.md, Theorem 7.1).

**Summary of findings:**

| Topic | Result | Section |
|-------|--------|---------|
| Regularity classification of K | alpha = -1 (from delta), K_bg is 0^- | 1 |
| K fits regularity structure framework | Yes, as a singular kernel | 1 |
| Renormalization of delta coefficient | Structurally suggestive but not constraining | 2 |
| BPHZ --> Connes-Kreimer --> zeta connection | Algebraically real, analytically incomplete | 3 |
| KPZ / stochastic quantization connection | Metaphorical only | 4 |
| Paracontrolled Littlewood-Paley for K | Natural decomposition, no new positivity | 5 |
| Local regularity --> global CPD-1 | **No** (categorical mismatch) | 6 |
| Overall viability for proving CPD-1 | **Not viable** | 7 |

---

## 1. The Weil Kernel's Regularity Classification

### 1.1 Holder Regularity of Distributions

In Hairer's framework, a distribution u on R^d is classified by its local Holder exponent alpha in R: informally, u in C^alpha means u "looks like |x|^alpha" at small scales. Precisely, u in C^alpha (Besov space B^alpha_{infty,infty}) iff for all Littlewood-Paley blocks:

||Delta_j u||_{L^infty} lesssim 2^{-j alpha}

For alpha >= 0: u is a genuine Holder-continuous function.
For alpha < 0: u is a distribution (not pointwise defined).

Key examples:
- Smooth function: alpha = +infty
- Holder-continuous f in C^alpha: exponent alpha in (0,1)
- Bounded measurable: alpha = 0
- Dirac delta on R: alpha = -1 (since <delta, phi^lambda_0> = phi(0) ~ lambda^0, but delta is in C^{-1-epsilon} for all epsilon > 0)
- White noise on R: alpha = -1/2 - epsilon
- Derivative of delta: alpha = -2

### 1.2 Components of the Weil Kernel

**The delta:** delta(x) has Holder regularity alpha = -1 on R. This is the worst singularity in K.

**The background K_bg:** Near x = 0, the digamma expansion gives:

K_bg(x) = K_bg(0) - C x^2 + O(x^4)

(since K_bg is even and smooth). So K_bg is C^infty away from any singularity. However, K_bg(x) ~ -(1/(2pi)) log|x| as |x| -> infty, so it is not in any global Holder space. Locally near any point, K_bg in C^infty. The "regularity" in the sense relevant to regularity structures is alpha = 0^- (the log singularity of the digamma function at s = 0 translates to a log divergence, which is just below Holder 0).

More precisely: Re[psi(1/4 + ix/2)] has a pole-type singularity at x = -i/2 (i.e., s = 0 in psi(s)), but for real x this pole is in the complex plane. For real x, K_bg is real-analytic. The regularity classification is alpha = +infty locally (smooth) but alpha = 0^- globally (log growth).

**The zero oscillation K_zeros:** Each term cos(gamma x)/(1/4 + gamma^2) is smooth and bounded. The sum converges absolutely (Sum 1/(1/4 + gamma^2) < infty by the Hadamard product). So K_zeros is in C^infty globally. Regularity alpha = +infty.

**The full kernel:**

K(x) = delta(x) + [smooth function]

has regularity alpha = -1, entirely due to the delta.

### 1.3 Regularity Structure for K

To encode K in a regularity structure T = (A, T, G):

**Index set:** A = {-1, 0, 1, 2, ...} (including -1 for the delta singularity).

**Model space:** T_{-1} = span{Xi} where Xi represents the delta distribution. T_0 = R (constants). T_k = span of k-th order polynomials for k >= 1.

**Model:** Pi_x(Xi) = delta(x - .) (the delta centered at x). Pi_x(1) = 1. Pi_x(X^k) = (. - x)^k.

The Weil kernel's local description at each point x_0:

K near x_0 = delta(x_0) . [Xi at x_0] + K_bg(x_0) . [1] + K_bg'(x_0) . [X] + ...

This is the Taylor-like expansion in the regularity structure, with the leading singular term Xi (at homogeneity -1) and smooth corrections in T_0, T_1, etc.

### 1.4 Significance

The regularity classification tells us:

1. K is a well-defined distribution (tempered, in fact) with singularity of order -1.
2. The singularity is "simple" — it is a single delta, not a derivative of delta or worse.
3. In the regularity structure hierarchy, K is at the same level as the noise in the 1D KPZ equation (where xi has regularity -3/2 but the solution h has regularity 1/2).
4. Products involving K (such as the bilinear form sum c_i c_j K(x_i - x_j)) are well-defined as long as the test functions are sufficiently regular — which discrete point evaluations ARE, since c_i delta_{x_i} is a measure.

**What this does NOT tell us:** Regularity classification describes the local singularity structure but says NOTHING about positivity, positive definiteness, or spectral properties. A distribution can have regularity alpha = -1 and be either positive definite or not. The delta itself is PD (its Fourier transform is 1 > 0), but -delta is equally regular and is negative definite.

---

## 2. Renormalization of the Delta Singularity

### 2.1 The Renormalization Question

In the Weil kernel K(x) = 1 . delta(x) + K_bg(x) + K_zeros(x), the coefficient of delta is exactly 1. In the matrix formulation:

K(0) = 1 + K_bg(0) + K_zeros(0) approx 1 + 1.528 + ... approx 2.528

**Question:** Is the coefficient "1" in front of delta determined by a renormalization condition? In regularity structures, when products of distributions require renormalization, the counterterms are constants (or lower-order terms) that absorb the divergences. Could the delta coefficient be such a counterterm?

### 2.2 The Analogy with KPZ Renormalization

In the KPZ equation d_t h = Delta h + (nabla h)^2 + xi, the ill-defined product (nabla h)^2 requires renormalization:

(nabla h)^2 --> (nabla h)^2 - C_epsilon

where C_epsilon = E[(nabla h_epsilon)^2] diverges as the regularization epsilon -> 0. The renormalized equation is:

d_t h = Delta h + :(nabla h)^2: + xi

with the Wick-ordered product :(nabla h)^2:. The "infinite constant" C_epsilon is a renormalization counterterm.

**For the Weil kernel:** The product structure in the CPD-1 condition is:

W(f * tilde{f}) = sum_rho |F(rho)|^2

The LHS involves the bilinear form sum c_i c_j K(x_i - x_j), where K(0) appears on the diagonal. If we "regularize" the delta by replacing delta(x) with delta_epsilon(x) = (1/epsilon) phi(x/epsilon) and take epsilon -> 0, the diagonal terms K(0) diverge as 1/epsilon. The "renormalization" is to subtract this divergence:

K^{ren}(x_i - x_j) = K(x_i - x_j) - (1/epsilon) phi((x_i - x_j)/epsilon) for i = j

This gives K^{ren}(0) = K_bg(0) + K_zeros(0), with the delta removed.

### 2.3 The Coefficient as a Normalization

The coefficient "1" in front of delta comes from the Fourier transform delta_hat(xi) = 1. In the Weil explicit formula:

W(f * tilde{f}) = sum_rho |F(rho)|^2 = |F(0)|^2 + |F(1)|^2 + (spectral integral)

The pole terms |F(0)|^2 + |F(1)|^2 contribute rank-2 pieces. The "1" from the delta in K_hat(xi) = 1 + K_bg_hat + K_zeros_hat is the contribution from the unit function in the Fourier domain, which relates to the residue of zeta at s = 1.

**Renormalization interpretation:** In the Connes-Kreimer framework, the residue of zeta at s = 1 is a "counterterm" — it appears in the Birkhoff factorization of the gamma-function of the underlying theory. The coefficient 1 in front of delta can be viewed as the simplest renormalization constant: it is the value that makes the Weil functional W well-defined on test functions that integrate to zero (i.e., on the primitive subspace).

### 2.4 Assessment

The renormalization interpretation is **structurally suggestive** but **analytically empty**. The coefficient 1 is determined by:
- The normalization of the Weil explicit formula (conventional)
- The residue of zeta at s = 1 (which is 1, unconditionally)
- The Fourier transform convention F[delta] = 1

None of these involve a choice or a divergent subtraction in the regularity-structures sense. The delta in K is not a "product of distributions requiring renormalization" — it is a single distribution that appears naturally in the kernel. The renormalization language is forced, not illuminating.

---

## 3. The Algebraic Chain: BHZ --> Connes-Kreimer --> Zeta

### 3.1 Bruned-Hairer-Zambotti (2019)

The BPHZ renormalization of regularity structures (Inventiones Math. 2019) is built on two Hopf algebras in cointeraction:

**H_+** (Butcher-Connes-Kreimer type): Coproduct Delta_+ implements recentering of models (translating Taylor-like expansions between basepoints). This is a direct descendant of the Connes-Kreimer coproduct on rooted trees.

**H_-** (extraction-contraction): Coproduct Delta_- implements extraction of divergent subforests, the analogue of the BPHZ forest formula. The counterterms are computed via the twisted antipode of H_-.

**Cointeraction:** Delta_- is a coalgebra morphism w.r.t. Delta_+, ensuring compatibility of renormalization and recentering.

The objects are **decorated colored forests**: trees/forests decorated with multi-indices on edges and vertices, encoding iterated stochastic integrals and their singular products.

### 3.2 Connection to Connes-Kreimer

The Connes-Kreimer Hopf algebra H for a renormalizable QFT is generated by 1PI Feynman graphs, with coproduct:

Delta(Gamma) = Gamma (x) 1 + 1 (x) Gamma + sum_{gamma subset Gamma, gamma divergent} gamma (x) (Gamma / gamma)

The Birkhoff factorization of characters of H gives the renormalized theory:

gamma(z) = gamma_-(z)^{-1} . gamma_+(z)

where gamma_- encodes counterterms and gamma_+ gives renormalized values.

**The BHZ extraction-contraction coproduct is the regularity-structures analogue of the Connes-Kreimer coproduct.** The Butcher-Connes-Kreimer Hopf algebra of rooted trees is a quotient of H_+. The full algebraic structure is:

Singular SPDEs --> (regularity structures) --> Decorated forest Hopf algebras (BHZ)
                                                            |
                                           (extraction-contraction = BPHZ)
                                                            |
                                                            v
                                                 Connes-Kreimer Hopf algebra
                                                            |
                                              (Birkhoff factorization = RH problem)
                                                            |
                                                            v
                                              Cosmic Galois group (motivic)
                                                            |
                                              (Tannakian = mixed Tate motives)
                                                            |
                                                            v
                                             Periods, MZVs, zeta values

### 3.3 Connes' Program and Zeta

Three pathways connect Connes-Kreimer renormalization to the Riemann zeta function:

**(a) Bost-Connes system.** The C*-dynamical system (C*(Q/Z) rtimes N^x, sigma_t) has partition function Z(beta) = zeta(beta). Phase transition at beta = 1: unique KMS state for beta <= 1, Galois-orbit of extremal KMS states for beta > 1. RH relates to positivity of certain functionals.

**(b) Adele class space trace formula.** On X = A_k / k^*, the Weil explicit formula becomes a trace formula:

sum_{rho: zeta(rho)=0} h(rho/i) = (principal value integral involving h-hat and arithmetic data)

RH is equivalent to: Tr(f * f^*) >= 0 on a suitable function space on the idele class group. This IS the Weil positivity criterion — the same condition as CPD-1 of the Weil kernel (bochner-proof.md, Theorem 7.1).

Recent work (Connes 2020, arXiv:2006.13771) identifies the "root of positivity" as the trace of the scaling action compressed onto the orthogonal complement of range-of-cutoff projections.

**(c) Cosmic Galois group.** Connes-Marcolli (2004): divergences in QFT are organized by a pro-unipotent group U (the cosmic Galois group), whose Lie algebra is freely generated by one generator e_n in each degree n >= 1. The category of equisingular flat vector bundles forms a Tannakian category equivalent to representations of U, which is reminiscent of the category of mixed Tate motives over Z.

### 3.4 Does the Chain Help?

The algebraic chain is:

BHZ (SPDE renormalization) --> Connes-Kreimer (QFT renormalization) --> Motivic Galois theory --> Zeta values

Each arrow is **mathematically rigorous**. But the composition does NOT yield a tool for proving CPD-1, because:

1. **The chain goes in the wrong direction.** It starts from SPDEs and arrives at zeta VALUES (special values zeta(n), multiple zeta values). The CPD-1 question concerns zeta ZEROS. Values and zeros are related by analytic continuation, but the algebraic framework does not cross this bridge.

2. **The Hopf algebra acts on counterterms, not on kernels.** The renormalization group of regularity structures acts on MODELS (local descriptions of solutions). The CPD-1 condition is a property of the KERNEL itself, not of any solution. The kernel K is an INPUT to the regularity structure, not an OUTPUT of the renormalization procedure.

3. **Connes' positivity criterion is the SAME as CPD-1.** The adele-class-space reformulation Tr(f * f^*) >= 0 is precisely the Weil positivity criterion, which IS CPD-1 (bochner-proof.md). Connes reformulates the problem in noncommutative geometry language but does not solve it. The regularity-structures algebraic tools connect to Connes' framework but inherit the same open problem.

---

## 4. KPZ, Stochastic Quantization, and Arithmetic

### 4.1 KPZ and the Prime Counting Function

The KPZ equation in 1+1 dimensions:

d_t h = Delta h + (nabla h)^2 + xi

models a growing interface with noise. The prime counting function pi(x) satisfies a "noisy step" equation with the von Mangoldt function Lambda(n) as "noise":

pi(x) = Li(x) + O(x^{1/2} log x) (under RH)

The error term pi(x) - Li(x) oscillates with amplitude ~ x^{1/2}, governed by the zeros of zeta.

**The analogy:**
- KPZ: smooth solution + nonlinear self-interaction + stochastic noise
- pi(x): smooth approximation Li(x) + nonlinear corrections + "arithmetic noise" from zeros

**Why the analogy fails:** The KPZ noise xi is STOCHASTIC (random, with known distribution). The "arithmetic noise" from zeta zeros is DETERMINISTIC (fixed, unknown). Regularity structures handle stochastic noise by averaging over realizations; there is no ensemble to average over for the zeta zeros.

### 4.2 Stochastic Quantization of Arithmetic QFT

The "arithmetic QFT" from the codebase constructs a QFT-like object from the Weil distribution. In stochastic quantization, a QFT measure is the invariant measure of a Langevin SPDE:

d_t phi = -delta S / delta phi + xi

**Question:** Is there an SPDE whose stationary measure is the Weil distribution?

If so, the positivity of the Weil distribution (= CPD-1 = RH) might follow from the SPDE structure — e.g., from reflection positivity of the Langevin dynamics.

**Assessment:** This is **speculative and problematic.** The Weil distribution W is a REAL distribution on test functions g: R -> R, not a probability measure on field space. The "stationary measure" interpretation requires:
- A configuration space (what are the "fields"?)
- A dynamics (what is the Langevin equation?)
- A Gibbs measure interpretation (W as e^{-S} for some action S)

None of these have natural definitions. The Weil functional W(g) = sum_rho hat{g}(rho) is a linear functional on test functions, not an integral over a measure space. Reinterpreting it as a stationary measure requires a non-trivial (and currently non-existent) construction.

### 4.3 The Phi^4_3 Model Analogy

The stochastic quantization equation for Phi^4 theory:

d_t phi = Delta phi - phi^3 + xi

requires BPHZ renormalization to be well-defined (Hairer 2014, Gubinelli-Hofmanova 2021). The renormalization constants absorb divergences from products of distributions.

**Parallel to Weil kernel:** The Weil kernel K = delta + K_bg + K_zeros involves a "product" in the CPD-1 condition: sum c_i c_j K(x_i - x_j). When the points x_i coincide (i = j), K(0) involves the delta's "value at 0" — formally infinite, regularized by the matrix truncation.

In Phi^4_3, the cubic term phi^3 is renormalized by subtracting an infinite mass: phi^3 --> :phi^3: - 3 C_epsilon phi. Could the Weil kernel's delta coefficient "1" be a similar renormalization constant?

**Answer: No.** The delta in K is not a divergent product of distributions — it is a single distribution that appears in the explicit formula. The matrix truncation (evaluating K at finitely many log-prime points) is a discretization, not a regularization of a product. There is no singular product in the CPD-1 bilinear form that requires renormalization.

---

## 5. Paracontrolled Distributions and Frequency Decomposition

### 5.1 Littlewood-Paley Decomposition of K

The Weil kernel has a natural frequency decomposition via Littlewood-Paley blocks Delta_j:

**Low frequencies (j small):** Dominated by K_bg, which has K_bg_hat(xi) ~ 1/|xi| as xi -> 0. The low-frequency content of K is entirely from the archimedean (Gamma factor) contribution.

**Medium frequencies:** K_bg_hat(xi) = 2e^{-|xi|/2}/(1 - e^{-2|xi|}) contributes a smooth, exponentially decaying spectrum. This is the "bulk" contribution.

**High frequencies at discrete locations:** K_zeros_hat is an atomic measure supported at {pm gamma : gamma imaginary part of zero}. Each zero contributes a point mass at frequency xi = gamma.

**Flat contribution from delta:** delta_hat = 1 contributes uniformly at all frequencies.

### 5.2 Bony Paraproduct Structure

For two distributions f, g, the Bony decomposition:

f . g = f prec g + f circ g + f succ g

decomposes the product into paraproduct (low x high), resonant (comparable frequencies), and adjoint paraproduct (high x low).

For the bilinear form B(c) = sum c_i c_j K(x_i - x_j), we can decompose:

B(c) = B_delta(c) + B_bg(c) + B_zeros(c)

where each term inherits the frequency structure of the corresponding kernel component.

**B_delta(c) = sum c_i^2:** The delta contribution is purely diagonal. On the primitive subspace (sum c_i = 0), this contributes sum c_i^2 >= 0. Always non-negative.

**B_bg(c):** The background contribution. Since K_bg_hat > 0 everywhere (bochner-proof.md, Theorem 8.1), B_bg(c) >= 0 on primitives. The low-frequency dominance of K_bg means this contribution is strongest when the "test function" c has slowly varying coefficients.

**B_zeros(c):** The zero-oscillation contribution. Each zero gamma contributes at frequency gamma. The sign of this contribution depends on whether zeros are on the critical line (bochner-proof.md, Section 5.12).

### 5.3 Paracontrolled Approach to CPD-1

**Idea:** If K were the solution to some PDE/SPDE, one could use the paracontrolled machinery to decompose K = K prec X + K^sharp, where X is a reference distribution and K^sharp is "smoother." The CPD-1 condition for K would then reduce to conditions on X and K^sharp separately.

**Problem:** K is not a solution — it is a given kernel. The paracontrolled framework applies to SOLUTIONS of equations, not to arbitrary distributions. One would need to identify an equation that K satisfies.

**Does K satisfy an equation?** The Weil explicit formula can be written as:

K_hat(xi) = 1 + H(xi) + sum_gamma a_gamma delta(xi - gamma)

where H(xi) = K_bg_hat(xi) is a known smooth function. This is a DECOMPOSITION, not a DIFFERENTIAL EQUATION. There is no known PDE or SPDE whose solution is the Weil kernel.

### 5.4 Assessment

The Littlewood-Paley decomposition provides a clean spectral view of K's frequency content, but this is EXACTLY the Fourier-analytic approach already pursued in bochner-proof.md. The paracontrolled calculus adds no new information because:
- K is not a solution to an equation (no PDE structure to exploit)
- The product involved in CPD-1 (bilinear form in c) is well-defined (no singular product to renormalize)
- The frequency decomposition K = K_bg + K_zeros + delta is already the natural Littlewood-Paley-type split

---

## 6. Local Regularity vs. Global CPD-1

### 6.1 The Categorical Mismatch

**Regularity structures are LOCAL.** They describe the behavior of a distribution near each point x_0:

Pi_{x_0}(tau)(y) = distribution in y, near x_0

The regularity alpha classifies the local singularity: how the distribution "looks" at scale lambda -> 0 near x_0.

**CPD-1 is GLOBAL.** It requires:

sum_{i,j} c_i c_j K(x_i - x_j) >= 0

for ALL finite subsets {x_i} subset R and ALL c with sum c_i = 0. This is a condition on the ENTIRE function K, evaluated at ARBITRARY separations x_i - x_j.

### 6.2 Can Local --> Global?

**For smooth kernels:** If K is smooth and PD locally (in a neighborhood of each point), and K is continuous, then K is globally PD. This is because PD is equivalent to K_hat >= 0 (Bochner), and the Fourier transform of a continuous function is determined globally.

**For distributions:** The local regularity of K (its Holder exponent at each point) tells us nothing about K_hat's sign. Examples:

- delta(x): regularity -1, CPD-1 (delta_hat = 1 > 0)
- -delta(x): regularity -1, NOT CPD-1 (-delta_hat = -1 < 0)
- delta'(x): regularity -2, neither PD nor CPD-1 (delta'_hat = i xi, purely imaginary for real xi)

All three have the same (or worse) local singularity structure. CPD-1 depends on the GLOBAL spectral content (sign of K_hat), which is invisible to local regularity.

### 6.3 What Would Be Needed

To bridge local -> global, one would need a theorem of the form:

"If K is a translation-invariant distribution on R with regularity alpha >= -1, AND K satisfies some additional local positivity condition at every point, THEN K is CPD-1."

**No such theorem exists or can exist in this generality.** The counterexamples above (delta vs. -delta) show that local conditions alone cannot determine global positivity.

One might add structural conditions: "K = delta + smooth, where the smooth part is CPD-1." Then K is automatically CPD-1 (sum of CPD-1 distributions). But this is the trivial observation that delta + CPD-1 = CPD-1, which requires ALREADY KNOWING the smooth part is CPD-1 — i.e., already knowing RH (bochner-proof.md, Section 9).

### 6.4 The Reconstruction Theorem Doesn't Help

Hairer's reconstruction theorem: given a modelled distribution f in D^gamma, there exists a unique distribution Rf consistent with the local descriptions Pi_x(f(x)). This is a LOCAL-to-GLOBAL result, but for EXISTENCE, not for POSITIVITY. The reconstructed distribution Rf is determined by its local Taylor-like descriptions, but its Fourier transform's sign is not constrained by the reconstruction.

---

## 7. Honest Assessment

### 7.1 What Regularity Structures Provide for the Weil Kernel

1. **A precise classification of K's singularity:** K in C^{-1} (regularity -1), with the singularity entirely from the delta at x = 0. The background and zero-oscillation parts are smooth.

2. **A framework for handling the delta rigorously:** Any computation involving K can be made rigorous using distribution theory (Schwartz) or regularity structures. The CPD-1 bilinear form is well-defined for any finite set of points.

3. **An algebraic connection to renormalization in QFT and number theory:** The chain BHZ --> Connes-Kreimer --> motivic Galois theory --> zeta values is mathematically rigorous. It shows that the algebraic structures underlying SPDE renormalization are closely related to those underlying the zeta function.

### 7.2 What Regularity Structures Do NOT Provide

1. **Any progress toward CPD-1:** The CPD-1 condition is a global spectral positivity requirement (K_hat(xi) >= 0 for xi != 0). Regularity structures analyze local singularities, not global spectral properties. The two are categorically different.

2. **A useful renormalization of the Weil kernel:** The delta in K is not a divergent product requiring renormalization — it is a datum of the explicit formula. The renormalization machinery of regularity structures has nothing to renormalize in the CPD-1 problem.

3. **A bridge from local to global positivity:** No theorem connects local regularity of a distribution to global positive definiteness. Such a theorem cannot exist in generality (delta vs. -delta counterexample).

4. **An SPDE whose solution or stationary measure is K:** No equation for K is known. Without an equation, the regularity structure / paracontrolled machinery has no dynamics to exploit.

### 7.3 The Obstruction Compared to Other Approaches

| Approach | Can handle K_bg? | Can handle K_zeros? | Obstruction |
|----------|-----------------|--------------------| ------------|
| Bochner (bochner-proof.md) | Yes (K_bg_hat > 0) | No (K_zeros_hat sign depends on RH) | Fourier sign of K_zeros |
| Schoenberg (schoenberg-attempt.md) | Yes (e^{tK_bg} PD) | No (e^{tK_zeros} PD only under RH) | Sine terms from off-line zeros |
| GMC (gmc-approach.md) | Yes (log-correlated field) | No (statistical vs. deterministic) | Averaged vs. pointwise |
| RMT (rmt-universality-approach.md) | Yes (background spectral gap) | No (no universality for deterministic matrices) | Statistical vs. deterministic |
| **Regularity structures** | Yes (smooth, CPD-1 by Fourier) | No (local analysis, no global positivity) | **Local vs. global** |

Every approach successfully handles K_bg and encounters a different manifestation of the same fundamental obstacle at K_zeros. For regularity structures, the obstacle is categorical: local singularity analysis cannot determine global spectral positivity.

### 7.4 Speculative Directions

Despite the negative assessment, two speculative directions merit mention:

**Direction A: The Connes trace formula as a "regularity structure."**

Connes' reformulation of the Weil explicit formula as Tr(f * f^*) >= 0 on the adele class space involves a trace on a space of operators. If this trace could be computed via a "regularity structure" on the adele class space — with the primes providing the model and the Gamma factor providing the kernel — then the BPHZ renormalization might give a canonical choice of counterterms that enforces positivity.

**Status:** Entirely speculative. The adele class space is not a manifold in the usual sense (it has a non-Hausdorff quotient topology), and regularity structures have not been defined on such spaces.

**Direction B: Noncommutative regularity structures.**

Recent work (arXiv:2509.07948, September 2025) extends regularity structures to noncommutative settings. Connes' program for RH is fundamentally noncommutative (it uses the adele class space as a noncommutative space). If noncommutative regularity structures could be combined with Connes' trace formula, the renormalization group might act on the space of "models" for the Weil distribution in a way that constrains its spectral properties.

**Status:** The noncommutative regularity structures paper does not mention number theory or zeta functions. The connection is structural (both involve noncommutative algebras + renormalization), not proven.

---

## 8. Summary

### 8.1 The Verdict

**Regularity structures are NOT a viable path to proving CPD-1 of the Weil kernel.** The framework is designed for LOCAL singularity analysis and STOCHASTIC dynamics, while CPD-1 is a GLOBAL positivity condition on a DETERMINISTIC distribution. The mismatch is categorical.

### 8.2 What Is Gained

The investigation clarifies the nature of the Weil kernel's singularity (regularity -1, from the delta) and maps out the algebraic chain connecting SPDE renormalization to number theory:

$$\text{Singular SPDEs} \xrightarrow{\text{BHZ}} \text{Hopf algebras} \xrightarrow{\text{CK}} \text{Motivic structure} \xrightarrow{\text{periods}} \text{Zeta values}$$

This chain is algebraically rigorous at each step but does not compose to a tool for CPD-1.

### 8.3 Comparison to the AMR Approach

The AMR approach (Baker -> entropy positivity -> Rudolph -> Haar -> CPD-1) avoids the local/global mismatch entirely. It proves CPD-1 at arithmetic points by showing that the arithmetic measure mu_ar must be Haar (a GLOBAL structure theorem), and that Haar measure produces the correct eigenvalue signs. The regularity-structures approach cannot replicate this because:
- It has no analogue of the measure-rigidity step (Baker's theorem, Rudolph's classification)
- It cannot constrain global spectral properties from local data
- The "renormalization" it performs is on models/solutions, not on spectral measures

---

## References

### Regularity Structures
- Hairer, M. (2014). A theory of regularity structures. *Invent. Math.* 198(2), 269-504. [arXiv:1303.5113](https://arxiv.org/abs/1303.5113)
- Hairer, M. (2014). Singular stochastic PDEs. *Proceedings of the ICM, Seoul*. [arXiv:1403.6353](https://arxiv.org/abs/1403.6353)
- Bruned, Y., Hairer, M., and Zambotti, L. (2019). Algebraic renormalisation of regularity structures. *Invent. Math.* 215(3), 1039-1156. [arXiv:1610.08468](https://arxiv.org/abs/1610.08468)

### Paracontrolled Distributions
- Gubinelli, M., Imkeller, P., and Perkowski, N. (2015). Paracontrolled distributions and singular PDEs. *Forum Math. Pi* 3, e6. [arXiv:1210.2684](https://arxiv.org/abs/1210.2684)

### Stochastic Quantization
- Gubinelli, M. and Hofmanova, M. (2021). A PDE Construction of the Euclidean Phi^4_3 Quantum Field Theory. *Comm. Math. Phys.* 384, 1-75. [arXiv:1810.01700](https://arxiv.org/abs/1810.01700)

### Connes-Kreimer and Zeta
- Connes, A. and Kreimer, D. (2000). Renormalization in quantum field theory and the Riemann-Hilbert problem I. *Comm. Math. Phys.* 210, 249-273. [arXiv:hep-th/9912092](https://arxiv.org/abs/hep-th/9912092)
- Connes, A. and Marcolli, M. (2004). Renormalization and motivic Galois theory. *Int. Math. Res. Not.* 2004(76), 4073-4091. [arXiv:math/0409306](https://arxiv.org/abs/math/0409306)
- Connes, A. (2020). Weil positivity and Trace formula the archimedean place. [arXiv:2006.13771](https://arxiv.org/abs/2006.13771)

### Noncommutative Extensions
- Noncommutative Regularity Structures (2025). [arXiv:2509.07948](https://arxiv.org/abs/2509.07948)

### Weil Kernel and CPD-1
- [bochner-proof.md](bochner-proof.md) — Fourier analysis, CPD-1 equivalence to RH
- [schoenberg-attempt.md](schoenberg-attempt.md) — Schoenberg representation, K_zeros obstruction
- [gmc-approach.md](gmc-approach.md) — Gaussian multiplicative chaos
- [rmt-universality-approach.md](rmt-universality-approach.md) — Random matrix universality

---

*Document: Regularity Structures Approach to the Weil Kernel*
*Part of the AMR (Arithmetic Measure Rigidity) proofs module*
*February 2026*
