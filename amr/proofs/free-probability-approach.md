# Free Probability and Free Convolution for the Weil Matrix Spectrum

## Status: Exploratory — identifies a viable theoretical framework but with significant gaps between the Weil matrix structure and standard free probability hypotheses

---

## 0. Goal

Investigate whether free probability (Voiculescu, 1980s) can prove that all eigenvalues of the Weil matrix M restricted to the primitive subspace are non-positive: spec(M|_prim) <= 0.

The Weil matrix decomposes as M = M_bg + M_zeros where:
- M_bg: digamma kernel (smooth, dominant, unconditionally NSD on primitives)
- M_zeros: zeta zero oscillations (small perturbation, ||M_zeros|| << spectral gap of M_bg)

Both are real symmetric. The question: can free probability give a *structural* proof that the sum preserves non-positive spectrum?

---

## 1. Free Probability: Core Tools

### 1.1 Free Independence (Voiculescu, 1983)

Two self-adjoint elements a, b in a non-commutative probability space (A, phi) are **freely independent** if their joint moments factor through the *free cumulant* formula: all mixed free cumulants vanish.

For N x N random matrices, A_N and B_N are **asymptotically freely independent** if:
- The empirical spectral measures mu_{A_N}, mu_{B_N} converge weakly to mu_A, mu_B
- The empirical spectral measure of A_N + B_N converges to mu_A boxplus mu_B (free additive convolution)

**Sufficient condition (Voiculescu, 1991):** A_N deterministic with converging spectral measure, B_N = U_N D_N U_N^* with U_N Haar-distributed on O(N) or U(N) and D_N deterministic with converging spectral measure => (A_N, B_N) are asymptotically free.

### 1.2 Free Additive Convolution (boxplus)

For probability measures mu, nu on R with compact support, the free convolution mu boxplus nu is defined by:

$$R_{\mu \boxplus \nu}(z) = R_\mu(z) + R_\nu(z)$$

where R_mu is the **R-transform** of mu, defined via the Cauchy transform:

$$G_\mu(z) = \int \frac{d\mu(t)}{z - t}, \quad K_\mu(G_\mu(z)) = z, \quad R_\mu(z) = K_\mu(z) - 1/z$$

**Key property (support):** If supp(mu) subset (-inf, a] and supp(nu) subset (-inf, b], then supp(mu boxplus nu) subset (-inf, a + b]. Free convolution *preserves* the non-positive support property.

### 1.3 R-transform

The R-transform linearizes free convolution (analogous to the log of the characteristic function for classical convolution). For a measure mu with moments m_1, m_2, ...:

$$R_\mu(z) = \sum_{n=0}^\infty \kappa_{n+1} z^n$$

where kappa_n are the free cumulants of mu. The first few:
- kappa_1 = m_1 (mean)
- kappa_2 = m_2 - m_1^2 (variance)
- kappa_3 = m_3 - 3m_2 m_1 + 2m_1^3

### 1.4 Free Entropy (Voiculescu, 1993)

The free entropy chi(mu) measures the "volume of matricial microstates" approximating mu:

$$\chi(\mu) = \lim_{N \to \infty} \left[\frac{1}{N^2} \log \text{vol}\{A \in M_N^{sa} : \mu_A \approx \mu\} + \frac{1}{2}\log N\right]$$

**Free Entropy Power Inequality (Szarek-Voiculescu, 1996):**

$$e^{2\chi(X+Y)} \geq e^{2\chi(X)} + e^{2\chi(Y)}$$

when X, Y are freely independent self-adjoint. This constrains the "spread" of the spectral distribution of the sum.

---

## 2. Are M_bg and M_zeros Asymptotically Freely Independent?

### 2.1 The Structure of Each Component

**M_bg** has entries:

$$M^{bg}_{(p,m),(q,n)} = -\frac{(\log p \cdot \log q)^{1/2}}{p^{m/2} q^{n/2}} K_{bg}(m\log p - n\log q)$$

where K_bg(x) = -(1/pi) Re[psi(1/4 + ix/2)] + log(pi)/(2pi). This is a smooth function of the log-prime difference. The eigenvectors of M_bg are determined by the **arithmetic geometry of the prime lattice** through the digamma function.

**M_zeros** has entries:

$$M^{zeros}_{(p,m),(q,n)} = -\frac{(\log p \cdot \log q)^{1/2}}{p^{m/2} q^{n/2}} K_{zeros}(m\log p - n\log q)$$

where K_zeros(x) = (1/2pi) sum_gamma 2cos(gamma x)/(1/4 + gamma^2). The eigenvectors of M_zeros are determined by the **oscillatory interference pattern** of zeta zeros at arithmetic points.

### 2.2 Arguments For Asymptotic Freeness

**Argument 1 (Eigenvector incoherence):** The eigenvectors of M_bg are determined by the smooth, monotone digamma function evaluated at log-prime differences. The eigenvectors of M_zeros are determined by the highly oscillatory cos(gamma_j * log(p/q)) terms where gamma_j are zeta zero ordinates. These involve fundamentally different structures:

- M_bg eigenvectors: smooth dependence on log p, slowly varying
- M_zeros eigenvectors: oscillatory dependence via sum of cos(gamma * log p) at incommensurate frequencies gamma_1, gamma_2, ...

The incommensurability of the gamma's (which follows from the simplicity of zeta zeros and the Hadamard product structure) means the eigenvectors of M_zeros rotate "quasi-randomly" relative to those of M_bg.

**Argument 2 (Scale separation):** As N -> infty (more primes), the matrix M_bg becomes increasingly dominated by its smooth structure (Toeplitz-like for equally spaced log-primes), while M_zeros involves an increasing number of incommensurate oscillatory components. The "direction" of M_zeros in matrix space becomes increasingly "generic" relative to M_bg.

**Argument 3 (GUE universality):** The local statistics of zeta zeros follow GUE (Montgomery-Odlyzko). If we model M_zeros as a function of GUE eigenvalues, then by universality results (Anderson-Farrell, 2014), the conjugation by the eigenvectors of M_zeros relative to M_bg is "asymptotically liberating."

### 2.3 Arguments Against Asymptotic Freeness

**Problem 1 (Both are deterministic):** Neither M_bg nor M_zeros is random. Standard asymptotic freeness results require at least one matrix to be conjugated by a Haar-random unitary (Voiculescu 1991) or an "asymptotically liberating" sequence (Anderson-Farrell 2014). Both our matrices are fully deterministic functions of primes and zeta zeros.

**Problem 2 (Common structure):** Both M_bg and M_zeros are functions of the SAME set of log-prime differences {m log p - n log q}. They share the same weight matrix W_{(p,m),(q,n)} = sqrt(log p log q)/(p^{m/2} q^{n/2}). This common algebraic structure creates correlations between their eigenvectors that could violate freeness.

**Problem 3 (No randomness model):** Free independence is a statement about joint moments. For deterministic matrices, the joint moments of M_bg and M_zeros are fixed numbers, not random variables. Asymptotic freeness would require showing these specific numbers converge to the free convolution prediction — a number-theoretic statement, not a probabilistic one.

### 2.4 Assessment

Asymptotic freeness of M_bg and M_zeros is **plausible but unproven**. The oscillatory nature of zeta zeros provides heuristic "randomness," but making this rigorous requires:

1. Identifying a random matrix model for M_zeros (e.g., via the CUE connection)
2. Proving the model captures the relevant spectral statistics
3. Applying a deterministic equivalent theorem (Male 2012, Pastur-Vasilchuk 2000)

**Difficulty level:** Very hard. This would likely require new results in the intersection of free probability and analytic number theory.

---

## 3. What Free Convolution Predicts

### 3.1 If Freeness Holds

Assume M_bg and M_zeros are asymptotically freely independent. Let:
- mu_bg = limiting spectral measure of M_bg|_prim
- mu_zeros = limiting spectral measure of M_zeros|_prim

Then:
- mu_M = mu_bg boxplus mu_zeros (free additive convolution)

**From known data:**
- supp(mu_bg) subset (-inf, -delta] where delta ~ 1.6-1.9 (spectral gap grows with N)
- mu_zeros is concentrated near 0 with ||M_zeros||_op <= 0.015

**Free convolution prediction:**
supp(mu_bg boxplus mu_zeros) subset (-inf, -delta + 0.015]

Since delta >> 0.015, the free convolution predicts **all eigenvalues remain strictly negative**. This is consistent with computational verification (all 159 primitive eigenvalues negative for 160 x 160 matrix).

### 3.2 R-transform Computation

**For M_bg:** The spectral measure mu_bg is approximately supported on [-2.3, -1.27e-8] (from 160x160 data). Its R-transform:

$$R_{bg}(z) = \kappa_1^{bg} + \kappa_2^{bg} z + \kappa_3^{bg} z^2 + \cdots$$

where kappa_1^{bg} = mean eigenvalue of M_bg|_prim (computable, approximately -0.5 for moderate N).

**For M_zeros:** The spectral measure mu_zeros is approximately supported on [-0.015, 0.015]. Its R-transform:

$$R_{zeros}(z) = \kappa_1^{zeros} + \kappa_2^{zeros} z + \cdots$$

where kappa_1^{zeros} ~ 0 (mean of zero oscillation is near zero by Hadamard identity) and kappa_2^{zeros} ~ 10^{-4} (tiny variance).

**Sum:** R_M(z) = R_bg(z) + R_zeros(z) ~ R_bg(z) + O(10^{-4})

The R-transform of M is dominated by the background R-transform, with negligible zero correction. The resulting spectral measure mu_M = mu_bg boxplus mu_zeros is a slight "smearing" of mu_bg — all eigenvalues remain negative with massive margin.

### 3.3 The Primitive Subspace Constraint

**Key subtlety:** Free convolution describes the FULL spectral measure, not the measure restricted to a subspace. The primitive subspace has codimension 1 (or 2 with pole terms).

**Resolution:** The primitive projection P = I - vv^T/N removes only the rank-1 direction v = (1,...,1)/sqrt(N). By interlacing (Cauchy):

$$\lambda_i(M|_{prim}) \leq \lambda_i(M) \leq \lambda_{i-1}(M|_{prim})$$

If all eigenvalues of M are <= 0 (which the free convolution predicts), then all eigenvalues of M|_prim are <= 0 except possibly the one interlaced with the removed direction. But the removed direction corresponds to the non-primitive eigenvector, which carries the pole contribution — its eigenvalue is the one that may be positive, and it is projected out.

Therefore, **free convolution on the full matrix implies the result on the primitive subspace**, provided the non-primitive eigenvalue is correctly identified with the pole contribution.

---

## 4. Connection to Keating-Snaith CUE Model

### 4.1 The CUE Framework

Keating and Snaith (2000) modeled zeta on the critical line by characteristic polynomials of the Circular Unitary Ensemble:

$$Z_N(\theta) = \det(I - e^{i\theta} U_N), \quad U_N \in U(N) \text{ Haar-distributed}$$

The key correspondence: N <-> log(T/(2pi)), where T is the height on the critical line.

**CUE statistics capture:**
- Pair correlation of zeros (Montgomery, 1973): R_2(u) = 1 - (sin pi u / pi u)^2
- Moments: E[|zeta(1/2+it)|^{2k}] matched by CUE predictions
- Spacing distribution: matches GUE (same as CUE for these statistics)

### 4.2 CUE Model for the Weil Matrix?

The Weil matrix M depends on zeta zeros through K_zeros. If we replace the "true" zeta zeros with eigenangles theta_1, ..., theta_N of a CUE matrix U_N, we get a **CUE model** for M_zeros:

$$K_{zeros}^{CUE}(x) = \frac{1}{2\pi} \sum_{j=1}^N \frac{2\cos(\theta_j x / (2\pi/\log T))}{1/4 + (\theta_j / (2\pi/\log T))^2}$$

(with appropriate scaling theta_j <-> gamma_j).

**In this model:**
- M_zeros^{CUE} is a random matrix (randomness from CUE)
- M_bg is deterministic
- By Voiculescu's theorem (1991): M_bg and U_N D U_N^* are asymptotically free if D is deterministic and U is Haar

**This would give:** M_bg and M_zeros^{CUE} are asymptotically free, and free convolution applies.

### 4.3 The Gap: CUE Model vs. Reality

The CUE model captures the *statistics* of zeta zeros but not their exact positions. The Weil matrix uses the ACTUAL zeros of zeta, not random CUE eigenvalues. The question is:

**Does asymptotic freeness hold for the actual zeros, not just for CUE-modeled zeros?**

This is a "universality" question: do the spectral statistics of M (computed from actual zeta zeros) match those computed from CUE zeros? By the Keating-Snaith philosophy, they should — but proving this rigorously is at the frontier of analytic number theory.

**Partial results:**
- Montgomery pair correlation (proven for restricted range under RH)
- Odlyzko numerical agreement (10^23 zeros match GUE to high precision)
- Rudnick-Sarnak n-point correlations (proven for restricted test functions)

None of these suffice to prove the full asymptotic freeness needed.

### 4.4 A CUE-Based Proof Strategy

If one could prove:

**Conjecture (CUE Universality for the Weil Matrix):** The empirical spectral distribution of M_zeros, constructed from zeta zeros up to height T, converges in distribution to the same limit as M_zeros^{CUE} constructed from CUE(N) eigenvalues with N = (1/2pi) log T.

Then:
1. M_bg and M_zeros^{CUE} are asymptotically free (by Voiculescu + Haar)
2. mu_M^{CUE} = mu_bg boxplus mu_zeros^{CUE} has support in (-inf, 0] (by Section 3)
3. By the universality conjecture, the actual mu_M has the same support
4. Therefore spec(M|_prim) <= 0

**Assessment:** This is a coherent proof strategy but relies on a conjecture (CUE universality) that is itself extremely difficult and currently unproven.

---

## 5. Deterministic Equivalents and Structured Matrices

### 5.1 Pastur-Vasilchuk Framework (2000)

For matrices of the form A_N + U_N B_N U_N^* where A_N, B_N are deterministic Hermitian and U_N is Haar unitary:

**Theorem (Pastur-Vasilchuk):** The empirical spectral measure of A_N + U_N B_N U_N^* converges a.s. to mu_A boxplus mu_B.

This is the cleanest result but requires Haar conjugation — not available for the Weil matrix.

### 5.2 Anderson-Farrell Framework (2014)

A sequence of unitary matrices {U_N} is **asymptotically liberating** if it "mixes" eigenbases sufficiently.

**Theorem (Anderson-Farrell):** If {U_N} is asymptotically liberating and A_N, B_N are deterministic with converging spectral measures, then:

$$\mu_{A_N + U_N B_N U_N^*} \to \mu_A \boxplus \mu_B$$

**Application to Weil matrix:** We need to find U_N such that M_zeros ~ U_N D_zeros U_N^* where D_zeros is diagonal and U_N is asymptotically liberating relative to the eigenbasis of M_bg.

The oscillatory structure cos(gamma * log(p/q)) suggests the eigenvectors of M_zeros are "rotated" by a matrix related to the discrete Fourier transform at zeta zero frequencies. If these frequencies are sufficiently "random" (in the Anderson-Farrell sense), liberating follows.

### 5.3 Male's Theorem (2012)

**Theorem (Male):** For polynomials P in GUE matrices and deterministic matrices with converging spectral distributions, the operator norm of P converges a.s. to the free probability prediction. In particular, **no eigenvalues escape** the limiting support.

**Relevance:** If M_zeros can be approximated by a GUE-type matrix, Male's theorem gives not just convergence of the spectral measure but **operator norm convergence** — meaning the largest eigenvalue of M|_prim converges to the edge of the free convolution support. This is exactly the "no positive eigenvalues" result we need.

**The catch:** M_zeros is not a polynomial in GUE matrices. It is a specific function of zeta zeros, which only *heuristically* resemble GUE eigenvalues.

---

## 6. Free Entropy Approach

### 6.1 Entropy of the Component Measures

**mu_bg (background spectral measure):**
- Supported on [-2.3, -delta] with delta ~ 10^{-8}
- Has density (continuous measure, not atomic for large N)
- Free entropy: chi(mu_bg) ~ -(1/2) int int log|s-t| dmu_bg(s) dmu_bg(t) + const

**mu_zeros (zero spectral measure):**
- Supported on [-epsilon, epsilon] with epsilon ~ 0.015
- Concentrated near 0 (most eigenvalues very small)
- Free entropy: chi(mu_zeros) is large (concentrated measures have high entropy in the free sense)

### 6.2 Free Entropy Power Inequality

If X ~ mu_bg and Y ~ mu_zeros are freely independent:

$$e^{2\chi(X+Y)} \geq e^{2\chi(X)} + e^{2\chi(Y)}$$

This gives a **lower bound on the spread** of the spectral distribution of M = M_bg + M_zeros. However, it constrains the *width* of the distribution, not its *support location*. For proving eigenvalue negativity, we need the support to stay in (-inf, 0], which the free EPI does not directly give.

### 6.3 What Free Entropy CAN Provide

The **free entropy-based concentration** (Guionnet, 2002):

$$P(\lambda_{max}(A_N + B_N) > t) \leq e^{-cN^2}$$

for t above the free convolution edge, when A_N and B_N satisfy freeness conditions. This exponential concentration means:

- The maximum eigenvalue of M is exponentially unlikely to exceed the free convolution prediction
- If the free convolution support is in (-inf, 0], the probability of a positive eigenvalue decays exponentially in N^2

For deterministic matrices, this "probability" statement needs reinterpretation (there's no randomness), but it suggests that "generic" perturbations preserve negativity with overwhelming margin.

---

## 7. Concrete Numerical Test

### 7.1 Protocol

To test whether free independence holds for the Weil matrix:

1. **Compute spectral distributions:** For N = 100 primes (m_max = 3, giving ~300-dim matrices):
   - Eigenvalues of M_bg|_prim -> empirical mu_bg
   - Eigenvalues of M_zeros|_prim -> empirical mu_zeros (expected: near-delta at 0)
   - Eigenvalues of M|_prim -> empirical mu_M

2. **Compute free convolution:** mu_bg boxplus mu_zeros numerically via:
   - R-transform: R_bg from mu_bg, R_zeros from mu_zeros
   - Sum: R_M = R_bg + R_zeros
   - Invert to get mu_bg boxplus mu_zeros

3. **Compare:** mu_M vs mu_bg boxplus mu_zeros
   - If Kolmogorov-Smirnov distance is small -> evidence for free independence
   - If the free convolution is supported on (-inf, 0] -> evidence for the approach

### 7.2 Expected Outcome

Given ||M_zeros|| << ||M_bg||, the free convolution mu_bg boxplus mu_zeros is approximately:

$$\mu_{bg} \boxplus \mu_{zeros} \approx \mu_{bg} \boxplus \delta_0 = \mu_{bg}$$

(free convolution with a point mass at 0 is the identity). So the prediction is:

$$\mu_M \approx \mu_{bg}$$

This is trivially true when the perturbation is small! The interesting regime would be when ||M_zeros|| is comparable to the spectral gap of M_bg — but for the Weil matrix, this never happens (the 110x margin from circularity-resolution.md ensures dominance).

### 7.3 Implication

The free probability prediction in this regime reduces to the **Weyl perturbation bound**:

$$|\lambda_i(M) - \lambda_i(M_{bg})| \leq ||M_{zeros}||_{op} \leq 0.015$$

which is already proven unconditionally (circularity-resolution.md, Strategy 3). Free probability doesn't add new information when the perturbation is this small.

---

## 8. Honest Assessment

### 8.1 Can Free Probability Prove Eigenvalue Negativity?

**In principle: Yes.** If one could establish:
1. Asymptotic freeness of M_bg and M_zeros (or a CUE model thereof)
2. That both component spectral measures have non-positive support
3. Apply Male's operator norm convergence to rule out escaping eigenvalues

Then spec(M|_prim) <= 0 would follow.

**In practice: Not with current tools.** The fundamental obstacles are:

| Obstacle | Nature | Difficulty |
|----------|--------|------------|
| Both matrices are deterministic | No Haar randomness to invoke Voiculescu | Very hard |
| CUE universality unproven | Connecting zeta zeros to CUE for spectral purposes | Frontier of ANT |
| Common algebraic structure | Both depend on same log-prime differences | May violate freeness |
| Small perturbation regime | ||M_zeros|| << spectral gap | Weyl bound already suffices |
| Infinite-dimensional limit | Need N -> infinity but matrices are finite | Structural issue |

### 8.2 What Free Probability DOES Provide

**Conceptual framework:** Free probability explains WHY the Weyl perturbation bound works so well for the Weil matrix — the oscillatory structure of M_zeros makes it "approximately free" from M_bg, meaning the eigenvalue perturbation is governed by the free convolution (which is the "generic" case for sums of matrices), not by worst-case alignment.

**Structural insight:** The fact that free convolution preserves non-positive support gives a "moral reason" for the negativity: two NSD matrices that are "generically oriented" remain NSD when added. The Weil matrix satisfies this because:
- M_bg is NSD on primitives (proven, actb-proof.md Theorem E)
- M_zeros has tiny spectral radius (proven, 0.015 bound)
- Their eigenbases are "incoherent" (heuristic, from oscillatory structure)

**Connection to random matrix theory:** Free probability is the algebraic backbone of the Keating-Snaith philosophy. The CUE model for zeta is, at its core, a free probability statement: the zeta function "behaves freely" relative to smooth test functions.

### 8.3 Comparison with Existing Approaches

| Approach | What it proves | Requires |
|----------|---------------|----------|
| **Weyl perturbation** (circularity-resolution.md) | spec(M|_prim) <= 0 for N x N truncations with 110x margin | Hadamard identity (unconditional) |
| **Measure rigidity** (actb-proof.md) | Full RH if mu_ar is ergodic | Ergodicity gap |
| **Free probability** (this document) | Conceptual explanation; would give full RH if CUE universality proven | CUE universality (open) |
| **BV + computation** (actb-proof.md) | Full RH if P_eff is feasible | Improved BV constants |

Free probability is **weaker than the Weyl perturbation approach** for finite truncations (it reduces to the same bound) but potentially **stronger for the infinite limit** (it could handle the N -> infinity regime where Weyl alone doesn't suffice).

### 8.4 The Real Value

The primary value of the free probability perspective is not as an independent proof technique but as:

1. **A bridge between the perturbative and structural approaches:** The Weyl bound handles finite truncations; measure rigidity handles the infinite limit. Free probability provides a unified framework where both are consequences of "freeness" at different scales.

2. **A potential route to CUE universality for the Weil matrix:** If one can prove that the Weil matrix "looks like" a CUE-derived matrix in a precise spectral sense, this would be a breakthrough connecting random matrix theory to the arithmetic of primes.

3. **A source of structural inequalities:** The free entropy power inequality, applied to the Weil matrix context, gives entropic bounds on the spectral distribution that complement the operator norm bounds from Weyl.

---

## 9. What Would Be Needed for a Rigorous Proof

### 9.1 Minimal Requirements

To make the free probability approach rigorous, one would need to prove ANY ONE of:

**(A) Deterministic asymptotic freeness for the Weil matrix.**
Show that for the specific sequences M_bg^{(N)} and M_zeros^{(N)} (truncated to N primes), the mixed moments satisfy the free moment formula asymptotically. This is a number-theoretic statement about products of digamma values and zeta zero oscillations.

**(B) CUE universality for K_zeros at arithmetic points.**
Show that the joint distribution of {K_zeros(m log p - n log q)} for primes p, q <= P is asymptotically the same as the corresponding distribution when zeta zeros are replaced by CUE eigenangles. This would follow from sufficiently strong n-correlation results for zeta zeros.

**(C) Anderson-Farrell liberating condition.**
Construct an explicit "liberating" sequence U_N such that M_zeros^{(N)} = U_N D_N U_N^* + o(1) where D_N captures the spectral measure and U_N is asymptotically liberating relative to M_bg^{(N)}'s eigenbasis.

### 9.2 Connections to Known Open Problems

Each of the above connects to major open problems:

- (A) relates to the **Mobius randomness law** and generalizations (correlations of arithmetic functions)
- (B) is a strengthening of the **Montgomery pair correlation conjecture** to matrix-level statistics
- (C) requires understanding the **eigenvector structure** of oscillatory kernel matrices, related to problems in harmonic analysis on the prime lattice

### 9.3 Most Promising Direction

**Male's theorem (2012) + CUE approximation** offers the most concrete path:

1. Approximate M_zeros by a polynomial in CUE matrices (using the spectral decomposition K_zeros = sum of cos(gamma * x) terms)
2. Apply Male's operator norm convergence to show no eigenvalues escape
3. The operator norm convergence gives spec(M|_prim) <= 0 directly

The main obstacle is Step 1: expressing K_zeros as a polynomial (or rational function) of CUE matrices requires a precise connection between the Weil matrix and the CUE characteristic polynomial.

---

## 10. Summary

| Question | Answer |
|----------|--------|
| Are M_bg and M_zeros asymptotically free? | Plausible but unproven; requires new number-theoretic results |
| Does free convolution predict negativity? | Yes — trivially, since ||M_zeros|| << spectral gap |
| Does this improve on Weyl perturbation? | Not for finite N; potentially for N -> infinity |
| Can CUE give a natural framework? | Yes, if CUE universality extends to matrix-level |
| Is this a viable proof route? | Only with major new results in ANT + free probability intersection |
| What's the main value? | Conceptual unification; bridge between perturbative and structural |

**Bottom line:** Free probability provides the correct *language* for understanding why the Weil matrix sum M = M_bg + M_zeros preserves negativity, but converting this understanding into a rigorous proof requires solving problems (CUE universality, deterministic freeness) that are themselves at the frontier of mathematics. For the current state of the AMR program, the Weyl perturbation bound (for finite N) and measure rigidity (for N -> infinity) remain the stronger approaches.

---

## References

### Free Probability
- Voiculescu, D. (1983). Symmetries of some reduced free product C*-algebras. *Operator Algebras and Their Connections with Topology and Ergodic Theory*, LNM 1132.
- Voiculescu, D. (1991). Limit laws for random matrices and free products. *Invent. Math.* 104, 201-220.
- Voiculescu, D. (1993). The analogues of entropy and of Fisher's information measure in free probability theory, I. *Comm. Math. Phys.* 155, 71-92.
- Szarek, S. & Voiculescu, D. (1996). Volumes of restricted Minkowski sums and the free analogue of the entropy power inequality. *Comm. Math. Phys.* 178, 563-570.

### Deterministic Equivalents and Asymptotic Freeness
- Pastur, L. & Vasilchuk, V. (2000). On the law of addition of random matrices. *Comm. Math. Phys.* 214, 249-286.
- Anderson, G. & Farrell, B. (2014). Asymptotically liberating sequences of random unitary matrices. *Adv. Math.* 252, 381-413.
- Male, C. (2012). The norm of polynomials in large random and deterministic matrices. *Probab. Theory Relat. Fields* 154, 477-532.

### Random Matrix Theory and Zeta
- Montgomery, H. (1973). The pair correlation of zeros of the zeta function. *Proc. Symp. Pure Math.* 24, 181-193.
- Keating, J. & Snaith, N. (2000). Random matrix theory and zeta(1/2+it). *Comm. Math. Phys.* 214, 57-89.
- Rudnick, Z. & Sarnak, P. (1996). Zeros of principal L-functions and random matrix theory. *Duke Math. J.* 81, 269-322.

### AMR Framework (Internal)
- [actb-proof.md](actb-proof.md) — ACTB via measure rigidity
- [circularity-resolution.md](circularity-resolution.md) — Spectral cancellation and unconditional bounds
- [bochner-proof.md](bochner-proof.md) — Fourier analysis of Weil kernel components
- [entropy-positivity.md](entropy-positivity.md) — Central duality theorem

---

*Generated as part of the Arithmetic Measure Rigidity framework*
*Task #11: Free probability exploration*
*Date: 2026-02-12*
