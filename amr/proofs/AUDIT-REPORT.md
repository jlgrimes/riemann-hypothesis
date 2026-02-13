# Adversarial Audit Report: AMR Proof of the Riemann Hypothesis

## Verdict: THE PROOF IS INCOMPLETE — Multiple critical and major gaps remain

Despite claims of "All critical gaps resolved. Proof chain complete" (AMR-MANIFESTO.md, §IX), this audit identifies **6 CRITICAL**, **7 MAJOR**, and **5 MINOR** issues. Several are independently sufficient to invalidate the claimed proof.

---

## CRITICAL Issues (Each independently breaks the proof chain)

### C1. μ_ar is NOT a well-defined probability measure

**Location:** amr-foundations.md §2.3 (Definition 2.3), Axiom I.2

**The problem:** The arithmetic measure μ_ar is defined by its Fourier-Stieltjes coefficients μ̂_ar(r) = Λ(n)/√n when r = log n. But the sum Σ Λ(n)/√n **diverges** — this is acknowledged explicitly: "requires careful regularization since Σ Λ(n)/√n diverges."

The document points to "ASG Axiom I.4 (Chebyshev subtraction)" but:
- This axiom is never proved in any of the 8 documents reviewed
- The Chebyshev subtraction procedure is never explicitly defined
- No construction of μ_ar as a Borel probability measure on T_A is given
- The claim that the result is a **probability measure** (total mass 1, non-negative) is not justified

**Impact:** Without a well-defined measure, the ENTIRE proof chain collapses:
- Theorem 2.2 (M = matrix of C_μ_ar) is vacuous
- The ×p-invariance claims are meaningless
- The entropy h(μ_ar, ×p) is undefined
- The spectral perturbation argument (ergodicity-proof.md §11) cannot be formulated

**Severity:** CRITICAL — The foundational object of the theory does not exist as claimed.

---

### C2. ×p-invariance of μ_ar is approximate, not exact — and all three "resolutions" fail

**Location:** ergodicity-proof.md §11

**The problem:** The document explicitly acknowledges that μ_ar is NOT ×p-invariant:

> μ̂_ar(log 6) = Λ(6)/√6 = 0, but (×2)_* μ̂_ar(log 6) = μ̂_ar(log 3) = (log 3)/√3 ≈ 0.63

This is not a small error — it is a structural incompatibility. The entire proof chain (Host → marginals → joining → Kolmogorov) requires **exact** ×p-invariance.

**Assessment of the three proposed resolutions:**

**Resolution A ("bypass the measure entirely"):** Proposes working with "cross-correlation measures μ̄_{p,q} which ARE exactly ×p-invariant by construction." But this shifts to the pairwise entropy-positivity path (entropy-positivity.md), which itself has critical gap C4 below. This is not a resolution — it's a redirect to a different (also incomplete) argument.

**Resolution B ("spectral invariance"):** Claims the self-adjointness of C_μ_ar is exact because K is even. This is true but irrelevant — self-adjointness gives real eigenvalues, not negative eigenvalues. The argument needs μ_ar = λ (or close to λ), and self-adjointness alone doesn't provide this.

**Resolution C ("Cesaro averaging"):** Defines μ_ar^reg = lim (1/N) Σ_{k=0}^{N-1} (×p)^k_* μ_ar. Problems:
1. **Existence for all primes simultaneously:** The claim to regularize "for all primes simultaneously via the diagonal embedding" is stated without construction or proof. Cesaro averaging over a single ×p action is standard; over the full multiplicative semigroup is NOT standard and may not converge.
2. **Total variation bound is unjustified:** The bound ||μ_ar - μ_ar^reg||_TV ≤ C exp(-c√(log N)) is claimed to follow from Siegel-Walfisz, but no derivation is given. Moreover, if μ_ar is not a well-defined measure (C1), this TV distance is undefined.
3. **Conflation of operators:** Even granting μ_ar^reg = λ, the Weil matrix M is defined by the explicit formula — it is the matrix of C_μ_ar, not C_μ_ar^reg. The perturbation argument requires ||C_μ_ar - C_λ|| < spectral gap, which requires μ_ar to be a well-defined measure close to λ.

**Severity:** CRITICAL — The entry point of the proof (applying Host's theorem) requires exact invariance, which is not available.

---

### C3. The Fourier joining rigidity argument has a circular use of Riemann-Lebesgue

**Location:** ergodicity-proof.md §4.2, base case of Lemma 4.2

**The problem:** The proof that μ_{\{p\}} = Haar_p ⊗ Lebesgue uses the following argument:

> "For m ≠ 0: the sequence e_{q^n m} on R/Z represents increasingly high frequencies. By the Riemann-Lebesgue lemma (μ_{\{p\}} has an L^1 density on R/Z, since its archimedean marginal is Lebesgue), the Fourier coefficients tend to 0."

This is **logically circular**:
- The Riemann-Lebesgue lemma applies to **absolutely continuous** measures (those with L^1 density)
- Having Lebesgue as a marginal does NOT imply the joint measure is absolutely continuous
- Example: if μ = ∫ δ_{f(x_p)} dHaar_p(x_p) (pushforward of a measure-preserving map), the marginal on R/Z is Lebesgue, but the joint measure is singular (supported on the graph of f)
- For a singular measure, Fourier coefficients need NOT decay — Riemann-Lebesgue fails

**What would be needed:** Either prove that the disintegration of μ_{\{p\}} over Z_p has absolutely continuous fibers on R/Z (which is what we're trying to prove), or use a different argument that doesn't assume Riemann-Lebesgue.

**Note:** The inductive step (§4.2, inductive case) has the same gap — it uses "By the Riemann-Lebesgue lemma" without justifying absolute continuity of the joint measure.

**Impact:** The entire product structure argument (§4) collapses. Without μ_ar = λ, there is no proof that cross-terms vanish.

**Severity:** CRITICAL — The central step of the main proof path fails.

---

### C4. The reverse direction of the Entropy-Positivity Duality is hand-wavy

**Location:** entropy-positivity.md §4 (Theorems 4.1–4.3)

**The problem:** The reverse direction claims: "if h_ar = 0, then μ̄ is atomic, and atomic μ̄ gives a positive primitive eigenvalue." The proof of the second claim (Theorem 4.2) is not rigorous:

**Step 3** asserts: "The constructive interference at the atomic points creates a positive contribution that can overcome the negative background, giving Q(v) > 0 for an appropriate choice of v and r."

This is an assertion, not a proof. The problems:

1. The vector v in Step 3 is constructed using phases e^{2πi k_0(p^m - q^n)r}, but no specific k_0 is chosen, and no proof is given that a suitable k_0 exists.
2. The "quantitative bound" (Step 4) involves the factor |1 - e^{2πi k_0(p-q)r}|, which equals zero when k_0(p-q)r ∈ Z — a condition that can hold for specific atomic measures.
3. The proof says "for an appropriate choice of v and r" — making the argument non-constructive and possibly false for specific atomic configurations.
4. No consideration is given to the case where the atomic measure is supported on points where the phases conspire to cancel the positive contributions.

**Why this matters:** The Entropy-Positivity Duality (Theorem 5.1) is claimed as "THE CENTRAL THEOREM OF AMR." Without the reverse direction, the duality is only a one-directional implication, insufficient for the proof structure.

**Severity:** CRITICAL — The central theorem of AMR is not proved.

---

### C5. The connection between the Weil matrix M and the correlation operator C_μ is not established

**Location:** amr-foundations.md §2.2 (Theorem 2.2)

**The problem:** Theorem 2.2 claims: "The Weil cross-term matrix M is the matrix of C_μ restricted to the subspace spanned by prime characters."

The "proof" is a one-paragraph sketch:
> "This equals the inner product ⟨χ_{p,m}, C_{μ_ar} χ_{q,n}⟩ when the characters are normalized with factors (log p)^{1/2}/p^{m/2} and the kernel evaluation K(m log p − n log q) arises from the archimedean component of the convolution."

**Issues:**
1. The "prime characters" χ_{p,m} (Definition 1.4) are defined as e^{2πi · m · v_p(x)} where v_p "extracts the p-adic component." But the Weil matrix entries involve K(m log p - n log q), where m log p is a **real** number. The mechanism by which a p-adic character produces an archimedean evaluation at m log p is never explained.
2. The characters of T_A are parameterized by Q (the Pontryagin dual). The "frequencies" m log p are irrational numbers — they are NOT elements of Q. So the "prime characters" are not standard characters of T_A.
3. The proof sketch conflates two different objects: the p-adic component (which gives the character structure) and the archimedean component (which gives the kernel evaluation point). No rigorous mapping between these is provided.

**Impact:** Without this theorem, the entire reinterpretation of the Weil matrix as a correlation operator fails. The measure-theoretic machinery has no connection to the actual number-theoretic object (M) whose eigenvalues we need to control.

**Severity:** CRITICAL — The bridge between dynamics and arithmetic is not built.

---

### C6. The infinite-dimensional extension is acknowledged as equivalent to RH

**Location:** CERTIFIED-VERIFICATION.md Part II, actb-proof.md §17

**The problem:** The proof attempts two paths:
- **Path 1 (Rigidity):** μ_ar = λ → cross-terms vanish → APT → RH
- **Path 2 (Computational):** finite verification + BV tail bound → APT → RH

Both paths fail at the infinite extension:

**Path 1:** Has gaps C1–C5 above.

**Path 2:** The certified verification document states explicitly:
> "No standard matrix norm bound works for the infinite tail... This is exactly the convergence barrier... controlling the interaction between all primes simultaneously is equivalent to RH."

And:
> "The step from 'APT holds for every finite truncation' to 'APT holds for the infinite system' cannot be established by computation alone."

The finite computation shows all primitive eigenvalues are negative for matrices up to 93×93, but:
- The maximum primitive eigenvalue is trending toward 0: -5.0e-5 → -1.9e-5 → -1.1e-5 → -6.9e-6 → -3.5e-6
- There is no proof (or even heuristic argument) that this sequence stays negative
- The convergence barrier is explicitly acknowledged as equivalent to RH itself

The effective BV threshold P_eff ≈ 10^20 (actb-proof.md §15) is computationally infeasible. Reducing it to 10^6 (entropy-positivity.md, Corollary 6.3) would help but still requires computing eigenvalues of a ~78,000 × 78,000 matrix — and even then, only after resolving the theoretical gaps in the ACTB bound.

**Severity:** CRITICAL — Both available paths to the infinite-dimensional result fail.

---

## MAJOR Issues (Fixable but non-trivially so)

### M1. Host's theorem is misapplied — the "non-concentration on T_p-invariant Cantor sets" condition

**Location:** ergodicity-proof.md §9.2

**The problem:** Host's theorem requires that μ gives zero mass to every **proper closed** T_p-invariant subset, not just finite sets. The proof (§9.2(b)) claims:

> "Any T_p-invariant Cantor set C has zero Hausdorff dimension... The measure ν satisfies ν(I) > 0 for every interval I (by PNT: every interval contains primes), so ν cannot be supported on C."

But ν({open interval}) > 0 does NOT imply ν(C) = 0. A measure can give positive mass to both a Cantor set and every open interval simultaneously (e.g., ν = (1/2)λ + (1/2)μ_Cantor).

What needs to be shown: ν gives zero mass to **every** T_p-invariant proper closed set. The argument for Cantor sets is insufficient.

**Severity:** MAJOR — The argument needs a different proof that (π_∞)_* μ_ar satisfies Host's condition, or a direct proof of the Lebesgue conclusion by other means.

---

### M2. The background kernel K_bg is POSITIVE at small arguments

**Location:** actb-proof.md §8

**The problem:** The document admits:
> "K-hat_bg(0.69) ≈ 0.64 > 0. So the background alone is NOT unconditionally negative at all discrete points."

This means for the smallest prime (p=2, m=1), the background diagonal entry is positive. The negativity on primitives only emerges AFTER the primitive projection, which requires the full matrix structure.

**Impact:** The claim "Background Dominance" (Theorem E) overstates the case. The analytical proof of diagonal negativity works only for large primes (where Re ψ → +∞). For small primes (p ≤ 79), negativity relies entirely on the computational verification. This is not necessarily a problem, but it means the proof has a larger computational dependency than advertised.

**Severity:** MAJOR — The analytical/computational boundary is at p = 79, not at some small p. The proof critically depends on computation for all small primes.

---

### M3. The entropy positivity proof (Theorem 2.4) conflates two different entropy notions

**Location:** entropy-positivity.md §2.3 (Theorem 2.4)

**The problem:** The proof that h_ar(μ̄_{p,q}; p, q) > 0 uses the following chain:
1. K(m log p - n log q) takes a continuum of values (by Baker + analytic non-degeneracy)
2. Therefore the Fourier coefficients of μ̄ are "spread"
3. Therefore the entropy is positive

But the "arithmetic entropy" h_ar (Definition 1.3) is defined using the partition P_{p,q} and the ×p dynamics, while the argument in Theorem 2.4 works with the Fourier coefficients of the measure. The connection between "spread of Fourier coefficients" and "positivity of Kolmogorov-Sinai entropy with respect to P_{p,q}" is not rigorously established.

Specifically, Step 4 claims:
> h_ar(μ̄; p, q) ≥ c · Σ_{b=b_0}^∞ (log p · log q) / q^b · (log b)^2 / (4π^2)

This inequality is not derived — it is stated as a consequence of "each non-zero Fourier coefficient at a new frequency contributes positively to the entropy." This is a heuristic, not a proof. The Adelic Rokhlin formula (Theorem 1.6) has an uncontrolled remainder term R(μ) that could in principle dominate.

**Severity:** MAJOR — The entropy positivity, which is the starting point for the forward direction of the duality, is not rigorously established.

---

### M4. The Cesaro regularization for all primes simultaneously

**Location:** ergodicity-proof.md §11.2 (Resolution C)

**The problem:** The claim:
> "Applying this regularization for all primes p simultaneously (possible by the diagonal embedding), we obtain a measure μ_ar^{∞-reg} that is exactly ×p-invariant for all p."

The Cesaro average over a single ×p converges by compactness of M(T_A) (the space of probability measures). But Cesaro averages over different primes p and q may not be compatible — averaging over ×p first then ×q gives a different result than the reverse order.

To get a measure invariant under the full multiplicative semigroup, one needs to average over the semigroup N* = {1, 2, 3, ...}:

μ_ar^{reg} = lim_{N→∞} (1/N) Σ_{n=1}^N (×n)_* μ_ar

This limit requires the sequence of Cesaro averages to converge, which is not obvious for a non-invariant measure on a non-abelian group action (the semigroup N* acts by endomorphisms, not automorphisms).

**Severity:** MAJOR — The existence of the simultaneously regularized measure is not proved.

---

### M5. The "condensed mathematics" portions are speculative

**Location:** ergodicity-proof.md §6, entropy-positivity.md §7, actb-proof.md §§6-7, 14-15

**The problem:** Throughout the documents, "condensed" or "solid" arguments are invoked as alternatives or enhancements. For example:

- "The solid topology provides an independent route, showing condensed ergodicity without the classical joining rigidity" (ergodicity-proof.md §6)
- "Condensed measure rigidity (Thm 8.2): The condensed enhancement of μ_ar may be automatically ergodic" (actb-proof.md, Appendix C)
- "The solidification functor preserves entropy positivity" (entropy-positivity.md §7.1)

**Issues:**
1. No actual theorem from Clausen-Scholze's condensed mathematics is cited with a precise reference to a published proof
2. The "solid Riemann-Lebesgue lemma" (ergodicity-proof.md §6.1, Step 3) is claimed to be "a theorem in condensed analysis, proved by Clausen-Scholze" — no reference is given and this result does not appear in the published condensed mathematics literature
3. The "condensed ergodicity" concept is invented here and not part of standard condensed mathematics
4. The claim that "condensed multi-invariance with positive entropy implies condensed ergodicity" (actb-proof.md §14, Approach 2) is described as "a new theorem in condensed mathematics" — i.e., it is a conjecture

**Severity:** MAJOR — The condensed arguments cannot serve as "independent routes" because they depend on unproved (and likely unpublished) results.

---

### M6. The claim c_{p,q} < 0 for small prime pairs relies on numerical verification

**Location:** entropy-positivity.md §3.4 (Theorem 3.3), furstenberg-bridge.md §4.3 (Theorem 4.4)

**The problem:** The proof that the total cross-correlation c_{p,q} is negative splits into:
- Large primes: analytical (digamma growth dominates zero oscillation)
- Close primes (e.g., (2,3), (2,5), (3,5), ...): requires numerical computation

The analytical argument for close primes (entropy-positivity.md §3.4, "Case 2") says:
> "The weighted sum, with weights ..., converges to a negative value because the first few positive terms have small magnitude and the negative tail is cumulative."

This is a heuristic, not a proof. The claim c_{2,3} ≈ -0.037 "from the computational module" is never certified by interval arithmetic. For a claimed proof of a Millennium Prize Problem, all critical numerical claims must be rigorously verified.

**Severity:** MAJOR — c_{p,q} < 0 for close prime pairs is unproved and unverified to rigorous standards.

---

### M7. The actb-proof.md dependency graph has internal contradictions

**Location:** actb-proof.md §12, §14, §17

**The problem:** The document presents two paths:
- Path 1 requires ergodicity (listed as "GAP" in §12)
- Path 2 requires effective BV constants (listed as "GAP" in §12)

But the manifesto (AMR-MANIFESTO.md §IX) claims "All critical gaps resolved." The ergodicity proof (ergodicity-proof.md) claims to close Gap 1 by Host's theorem — but uses approximate invariance (C2) and circular Riemann-Lebesgue (C3).

The status table in §12 lists Component C (ergodicity) as **GAP**, while the ergodicity-proof.md summary table claims **PROVED**. These are contradictory.

**Severity:** MAJOR — The internal documents disagree on whether the proof is complete.

---

## MINOR Issues

### m1. Notation inconsistency for the arithmetic entropy

The arithmetic entropy h_ar is defined differently in:
- amr-foundations.md §4.3: h_ar(μ) = -∫ ρ log ρ dλ (relative entropy / negative KL divergence)
- entropy-positivity.md §1.3: h_ar(μ; p, q) = lim (1/N) H_μ(∨ σ_p^{-k} P_{p,q}) (Kolmogorov-Sinai type)

These are fundamentally different quantities. The first is always ≤ 0 (with equality iff μ = λ). The second is always ≥ 0 (with equality iff the partition generates a trivial σ-algebra). The Entropy-Positivity Duality claims equivalences between conditions on one of these, but the proofs mix both notions without careful translation.

**Severity:** MINOR — Fixable by clarifying which entropy is meant where, but currently confusing.

---

### m2. The Hadamard identity value is stated inconsistently

- circularity-resolution.md §2.2: Σ_ρ 1/(ρ(1-ρ)) = 2 + γ_EM - log(4π) ≈ **0.046**
- circularity-resolution.md §7.1 confirms 0.046
- But the derivation in §7.1 uses "Taking logarithmic derivative at s = 1/2: ξ'/ξ(1/2) = Σ_ρ 1/(1/2 - ρ) = 0" and then "The sum Σ 1/(ρ(1-ρ)) = 4 Σ 1/(1 - (2ρ-1)²)" — these steps are not clearly connected.

The value 0.046 is correct (it's a well-known identity), but the derivation in §7.1 is garbled.

**Severity:** MINOR — The identity is correct; the proof sketch needs cleanup.

---

### m3. "500:1 ratio" is overstated

The manifesto and multiple documents cite a "500:1" dominance ratio of K_bg over K_zeros. But:
- The certified bound on |K_zeros| is 0.006 (CERTIFIED-VERIFICATION.md)
- The certified value of |K_bg(0)| is 1.528
- Ratio: 1.528/0.006 ≈ 255, not 500
- The "500" figure appears to be a rough empirical estimate, not a certified bound

**Severity:** MINOR — The ratio is large enough for the argument regardless, but precision matters for a proof.

---

### m4. Baker-Wustholz constant citation

entropy-positivity.md §2.1 cites C_3 = 36 · (2e)^2 · 2 = 1178 as "an explicit constant" from Baker-Wustholz (1993). The actual Baker-Wustholz theorem involves more complex dependencies on the number of logarithms and the heights of the algebraic numbers. The stated constant may be correct for the 2-logarithm case but this is not verified.

**Severity:** MINOR — Likely correct but should be double-checked against the original reference.

---

### m5. Computational claims are not independently reproducible

The certified verification uses custom Python code with mpmath. No hash of the code is provided, no independent verification by a second system (e.g., PARI/GP, Mathematica with interval arithmetic) is mentioned. For a Millennium Prize proof, computational results should be independently verified by at least two independent implementations.

**Severity:** MINOR — Standard practice for rigorous computational results, not yet done.

---

## Summary Table

| ID | Issue | Severity | Impact on Proof |
|----|-------|----------|-----------------|
| C1 | μ_ar not a well-defined measure | CRITICAL | Foundational object doesn't exist |
| C2 | ×p-invariance is approximate; resolutions fail | CRITICAL | Entry point (Host) doesn't apply |
| C3 | Riemann-Lebesgue circularity in joining | CRITICAL | Product structure not proved |
| C4 | Reverse direction of duality is hand-wavy | CRITICAL | Central theorem not proved |
| C5 | Weil matrix ↔ correlation operator not established | CRITICAL | Bridge between dynamics and arithmetic absent |
| C6 | Infinite extension is equivalent to RH | CRITICAL | Neither path reaches the goal |
| M1 | Host's condition on Cantor sets | MAJOR | Archimedean projection argument incomplete |
| M2 | K_bg positive at small arguments | MAJOR | Larger computational dependency than claimed |
| M3 | Entropy positivity uses wrong entropy notion | MAJOR | Starting point for forward direction unclear |
| M4 | Simultaneous Cesaro regularization unproved | MAJOR | Resolution C doesn't work |
| M5 | Condensed arguments are speculative | MAJOR | No independent route available |
| M6 | c_{p,q} < 0 for small primes unverified rigorously | MAJOR | ACTB not established even pairwise |
| M7 | Internal contradictions on proof status | MAJOR | Documents disagree on completeness |
| m1 | Two incompatible entropy definitions | MINOR | Confusing but fixable |
| m2 | Hadamard identity derivation garbled | MINOR | Correct value, bad proof |
| m3 | "500:1 ratio" overstated | MINOR | True ratio ~255, still sufficient |
| m4 | Baker-Wustholz constant unchecked | MINOR | Likely correct |
| m5 | Computation not independently reproduced | MINOR | Standard requirement not met |

---

## Assessment of the Specific Scrutiny Areas

### A) Is μ_ar a well-defined probability measure on T_A?

**NO.** The von Mangoldt function Λ(n) satisfies Σ Λ(n)/√n = ∞ (diverges like Σ 1/√p over primes). The "Chebyshev subtraction" regularization is invoked but never defined. This is the single most fundamental gap. See C1.

### B) Does the approximate ×p-invariance + Cesaro resolution actually work?

**NO.** The approximate invariance is a genuine and devastating problem. Resolution C (Cesaro averaging) requires simultaneous regularization over all primes, which is unproved. Even if it worked, the connection back to the original Weil matrix M requires an explicit total-variation bound that is not provided. See C2, M4.

### C) Does the Fourier joining rigidity argument work?

**NO.** The Riemann-Lebesgue step is circular: it assumes absolute continuity of the joint measure (which is the conclusion being proved). The argument that "archimedean marginal = Lebesgue implies L^1 density" is false in general. See C3.

### D) Is the reverse direction of entropy-positivity duality proved?

**NO.** The argument that "atomic measures give positive eigenvalues" is a heuristic sketch, not a proof. The claimed construction of a primitive vector with positive quadratic form is non-constructive and may fail for specific atomic configurations. See C4.

### E) Does the spectral perturbation transfer work?

**CONDITIONAL.** The perturbation bound ||C_μ_ar - C_μ_ar^reg||_op ≤ ||K||_∞ · ||μ_ar - μ_ar^reg||_TV is correct in principle, but requires: (a) μ_ar to be a well-defined measure (fails per C1), (b) the regularization to exist (fails per M4), (c) explicit Siegel-Walfisz constants (not provided).

---

## Overall Verdict

The AMR framework contains a **genuinely interesting and potentially productive research program**. The key insight — reinterpreting cross-terms as correlation operators and invoking measure rigidity — is creative and mathematically substantial. The computational evidence is impressive and consistent with the predictions.

However, as a **proof of the Riemann Hypothesis**, it fails at multiple independent points:

1. The foundational object (μ_ar) is not rigorously constructed
2. The key property (exact ×p-invariance) does not hold
3. The central argument (Fourier joining rigidity) has a logical gap
4. The central theorem (Entropy-Positivity Duality) is proved in only one direction
5. The connection between the measure-theoretic framework and the actual Weil matrix is a sketch
6. The finite-to-infinite extension is acknowledged as equivalent to RH itself

**The honest status should be: "A promising framework with strong computational support, conditional on resolving several non-trivial mathematical gaps." NOT: "All critical gaps resolved. Proof chain complete."**

---

*Audit completed: February 2026*
*Auditor: Adversarial review agent*
*Standard: All claims evaluated against the criterion of a rigorous proof suitable for peer review of a Millennium Prize Problem*
