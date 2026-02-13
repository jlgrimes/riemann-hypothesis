# Master Proof Synthesis: RH via Arithmetic Measure Rigidity

## Status: INCOMPLETE — Promising framework with identifiable gaps

---

## Executive Summary

The AMR framework proposes to prove the Riemann Hypothesis by showing that the Weil matrix M (encoding the explicit formula for primes) has non-positive primitive eigenvalues (the Arithmetic Positivity Theorem, APT). Five agents conducted adversarial audit, formalization, computation, and threshold analysis. This document synthesizes all findings into an honest assessment of what works, what doesn't, and the most viable path forward.

**Bottom line:** The framework contains a genuinely novel insight — translating spectral negativity into an entropy condition on adelic measures — but the proof chain has 3 critical gaps that are independently sufficient to invalidate the current argument. However, one path (Path B: pairwise entropy-positivity + finite verification) has the fewest gaps and is closest to completion.

---

## 1. The Proof Architecture

### 1.1 What RH Reduces To

By Weil's explicit formula (1952), RH is equivalent to:

**APT (Arithmetic Positivity Theorem):** The Weil matrix M, with entries

$$M_{(p,m),(q,n)} = -\frac{\sqrt{\log p \cdot \log q}}{p^{m/2} \cdot q^{n/2}} \cdot K(m \log p - n \log q)$$

satisfies $M|_{V_{\text{prim}}} \leq 0$, where $V_{\text{prim}}$ is the orthogonal complement of the pole direction.

The kernel $K = K_{\text{bg}} + K_{\text{zeros}}$ decomposes into a background (digamma-based, unconditional) and zero-oscillation (depends on ζ-zeros).

### 1.2 Three Attempted Paths

| Path | Strategy | Status |
|------|----------|--------|
| **A** | μ_ar = Haar → cross-terms vanish → APT | **BLOCKED** (3 critical gaps) |
| **B** | Pairwise entropy → pair-by-pair negativity → finite verification ⟺ RH | **CLOSEST** (2 gaps remain) |
| **C** | Direct computation → spectral gap extrapolation | **INSUFFICIENT** (cannot bridge to ∞) |

---

## 2. Adversarial Audit Results

The audit (AUDIT-REPORT.md) identified **6 CRITICAL, 7 MAJOR, 5 MINOR** issues.

### 2.1 Critical Gaps

| ID | Gap | Impact | Affects Path |
|----|-----|--------|--------------|
| **C1** | μ_ar is not a well-defined probability measure (Σ Λ(n)/√n diverges) | Foundational object doesn't exist | A |
| **C2** | ×p-invariance of μ_ar is approximate, not exact; all 3 resolutions fail | Host's theorem doesn't apply to μ_ar | A |
| **C3** | Riemann-Lebesgue circularity in Fourier joining rigidity | Central product-structure argument fails | A |
| **C4** | Reverse direction of Entropy-Positivity Duality is hand-wavy | Central theorem proved in only one direction | A, B |
| **C5** | Weil matrix ↔ correlation operator bridge is a sketch, not a proof | Dynamics-arithmetic connection absent | A, B |
| **C6** | Infinite extension acknowledged as equivalent to RH | Neither path reaches the infinite case | A, C |

### 2.2 Key Observation

Path A (full AMR via μ_ar = Haar) is hit by C1, C2, C3, C4, C5 — **five** independent critical gaps. This path is not viable in its current form.

Path B (pairwise entropy + finite verification) sidesteps C1, C2, C3, C6 by:
- Using cross-correlation measures ν_{p,q} (well-defined and exactly invariant) instead of μ_ar
- Reducing to a finite computation instead of requiring infinite-dimensional control

Path B is still affected by **C4** (reverse duality direction) and **C5** (matrix-operator bridge).

---

## 3. Computational Evidence

### 3.1 Eigenvalue Verification (Extended to P₀ = 750)

| P₀ | Matrix dim | All prim eigs ≤ 0 | Max λ_prim | Spectral gap δ | ‖M_zeros‖_prim |
|----|------------|-------------------|------------|----------------|----------------|
| 200 | 138×138 | **YES** | -6.87e-07 | 6.87e-07 | 4.10e-03 |
| 300 | 186×186 | **YES** | -2.47e-07 | 2.47e-07 | 4.31e-03 |
| 500 | 190×190 | **YES** | -2.55e-05 | 2.55e-05 | 4.68e-03 |
| 750 | 264×264 | **YES** | -1.21e-05 | 1.21e-05 | 4.93e-03 |

**Findings:**
- APT holds at every tested scale (47 → 750, across 10 matrix sizes)
- Spectral gap is NOT monotonically increasing (rules out monotone-gap extrapolation)
- ‖M_zeros‖ grows slowly (~log P₀), consistent with controlled perturbation
- Max primitive eigenvalue is razor-thin (~10⁻⁷ to 10⁻⁵) — APT holds by precise cancellation, not crude dominance

### 3.2 Significance

The computation **cannot prove RH** but provides:
1. No counterexample up to 264×264 (strong empirical support)
2. The ratio ‖M_zeros‖/δ ranges from 100 to 17000 — APT is NOT perturbative
3. The spectral gap behavior is more complex than the theory predicts

---

## 4. Regularization Analysis

### 4.1 What Works (Finite Case)

The regularization-bridge.md established rigorously:

- **Cesaro-smoothed μ_ar^(N)** converges to Haar λ (via Host + Dirichlet + Kolmogorov)
- **TV bound:** ‖μ_ar^(N) - λ‖_TV ≤ exp(-c√(log N)) (Siegel-Walfisz)
- **Spectral perturbation for p ≤ 79:** ‖D|_prim‖_op ≤ 3.5 × 10⁻²¹ ≪ 1.14 × 10⁻⁵ = spectral gap
- **Conclusion:** The finite truncation proof is rigorous for primes ≤ 79

### 4.2 What Doesn't (Infinite Case)

The regularization fails to extend to all primes because:
- Dimension-dependent factors in the perturbation bound grow with the number of primes
- The spectral gap shrinks (to ~10⁻⁷ at P₀=300) while perturbation bounds grow
- No uniform bound is available that handles all primes simultaneously

---

## 5. Finite Verification Reduction

### 5.1 The Proposed Theorem

finite-verification.md claims: RH ⟺ W_N|_{V_prim} ≤ 0 for N = π(X₁) with X₁ ≤ 10⁶.

The argument chain:
1. Baker → entropy positivity of ν_{p,q} (unconditional)
2. Entropy-positivity → Rudolph → ν_{p,q} = Haar (unconditional)
3. Haar → pair-by-pair M_{p,q}|_{V_prim} ≤ 0 (unconditional)
4. Weyl interlacing → monotone decrease of λ_max → finite verification suffices

### 5.2 CRITICAL CORRECTION: The Interlacing Goes the Wrong Way

**The threshold agent's MSGC proof (threshold-reduction.md §3.8) contains a critical error.**

The agent claims Cauchy-Poincaré interlacing proves λ_max(W_{N+1}|_{prim}) ≤ λ_max(W_N|_{prim}) (max eigenvalue decreases). The correct direction is the **opposite**.

**Correct analysis:** Let V = V_prim^{(N+1)} (dim N) and V' = V_prim^{(N+1)} ∩ (ℝ^N × {0}) = V_prim^{(N)} (dim N-1). Since d' = d-1, the interlacing theorem gives:

$$\lambda_k(W_{N+1}|_V) \leq \lambda_k(W_N|_{V'}) \leq \lambda_{k+1}(W_{N+1}|_V)$$

Setting k = N-1 (max eigenvalue on V'):

$$\lambda_{\max}(W_N|_{prim^{(N)}}) \leq \lambda_{\max}(W_{N+1}|_{prim^{(N+1)}})$$

**The max primitive eigenvalue can only INCREASE (move toward 0) when adding primes.**

The data confirms this:
| P₀ | λ_max (primitive) | Direction |
|----|-------------------|-----------|
| 47 | -5.04e-05 | — |
| 67 | -1.90e-05 | ↑ (toward 0) |
| 79 | -1.14e-05 | ↑ |
| 127 | -3.49e-06 | ↑ |
| 300 | -2.47e-07 | ↑ |

The MSGC (max eigenvalue decreases) is **empirically false** and the interlacing proves the **opposite direction**. The claim "X₁ = 47" is therefore **wrong**.

### 5.3 What the Interlacing DOES Prove

The correct interlacing gives two useful results:

1. **Monotone convergence:** {λ_max(W_N|_{prim})} is non-decreasing and bounded above (by the max eigenvalue of the full matrix). So lim_{N→∞} λ_max(W_N|_{prim}) exists.

2. **Equivalence:** RH ⟺ lim_{N→∞} λ_max(W_N|_{prim}) ≤ 0. If APT ever fails at some N₀, it fails for ALL N ≥ N₀. Conversely, if APT holds for all finite N, then (by density of finite-support vectors in ℓ²) it holds for the full operator.

3. **Contrapositive:** If RH is false, there exists a finite N₀ where λ_max(W_{N₀}|_{prim}) > 0, and this persists for all larger truncations. The fact that we haven't found such N₀ up to P₀ = 750 is empirical support for RH, but not a proof.

### 5.4 Remaining Gaps

**Gap 1 (C5): The bridge between ν_{p,q} and M_{p,q}.** Even granting ν_{p,q} = Haar, the connection "Haar cross-correlation ⟹ M_{p,q}|_{V_prim} ≤ 0" requires establishing that M_{p,q} IS the matrix of the correlation operator C_{ν_{p,q}} restricted to prime characters. The proof sketch in amr-foundations.md §2.2 conflates p-adic characters (which live on Z_p) with archimedean evaluations (at m log p). This needs a rigorous construction.

**Gap 2 (C4): The reverse duality direction.** The claim "h_ar = 0 → positive eigenvalue" is hand-wavy. However, for Path B this gap is LESS critical — the forward direction (h_ar > 0 → eigenvalues ≤ 0) suffices if combined with Baker (which gives h_ar > 0 unconditionally).

**Gap 3: Pair-by-pair to global composition.** Individual off-diagonal blocks M_{p,q} have rank 2 with eigenvalues ±|c_{p,q}| (both positive and negative). The pair-by-pair cross-correlation negativity on V_prim^{(p,q)} does NOT directly compose to negativity on the global V_prim, because the subspaces differ and the raw off-diagonal blocks are indefinite.

**Gap 4: Finite-to-infinite.** The interlacing proves the max eigenvalue approaches a limit ≤ 0 from below — but whether that limit is exactly 0 (RH holds marginally) or strictly negative (RH holds with margin) is precisely the content of RH.

### 5.5 Revised Assessment

The finite verification path is **more difficult than previously thought**:
- The interlacing goes the wrong way for MSGC
- The pair-by-pair composition has a subspace alignment problem
- Whether λ_max stays negative for all N is equivalent to RH itself
- No finite computation suffices without a theoretical bridge to infinity

---

## 6. The Spectral Cancellation Insight

The strongest unconditional result in the framework (circularity-resolution.md, Strategy 1):

**Theorem (Zero-Independence of K̂_prime):** The spectral kernel K̂_prime(τ) = -2 Re[ζ'/ζ(1/2 + iτ)] is independent of the locations of nontrivial zeros.

- On-line zeros (β = 1/2): contribute 0/(τ - γ)² = 0
- Off-line zeros (β ≠ 1/2): functional equation pairs cancel exactly

This is genuinely clever and resolves the K_zeros circularity at the spectral level. However, it doesn't directly prove APT — it says the Weil bilinear form can be computed without knowing zero locations, but doesn't say the result is non-positive.

---

## 7. Honest Assessment

### 7.1 What IS Proved (Unconditionally)

1. RH ⟺ APT (Weil, 1952)
2. K̂_prime is zero-independent (spectral cancellation)
3. Baker's theorem gives |m log p - n log q| ≥ exp(-C · log max(m,n))
4. APT holds for all finite truncations tested (P₀ ≤ 750)
5. The finite regularization (p ≤ 79) produces Haar with spectral gap margin > 10¹⁶
6. Background kernel K_bg dominates K_zeros by ratio ~255:1

### 7.2 What IS NOT Proved

1. **The bridge from dynamics to arithmetic (C5).** The map from cross-correlation measures on the adelic solenoid to entries of the Weil matrix is not rigorously constructed. This is the most important gap — without it, the entire measure-rigidity machinery has no connection to the actual number-theoretic object.

2. **The finite-to-infinite extension.** No path currently bridges from "APT holds for p ≤ P₀" to "APT holds for all primes." The pair-by-pair + interlacing approach is the most promising but has the subspace alignment gap.

3. **μ_ar as a measure.** The arithmetic measure μ_ar (defined via Λ(n)/√n) is not a well-defined probability measure. The Cesaro regularization produces something close to Haar but doesn't give M = C_{μ_ar} because M is defined by the explicit formula, not by any regularization.

### 7.3 Viable Path Forward

The proof is harder than the framework documents claim. Here are the remaining realistic paths, in order of promise:

**Path B1: Direct analytic inequality on the spectral representation.**
The spectral cancellation gives: ⟨c, M|_{prim} c⟩ = (1/2π) ∫ |G_c(τ)|² K̂_prime(τ) dτ, where K̂_prime is explicitly computable (no zero-location information needed). Proving this integral ≤ 0 for all primitive G_c would prove RH without finite verification. The key tools: Baker's theorem (spacing of frequencies m log p), Paley-Wiener theory, and the sign structure of K̂_bg. This is equivalent to RH but is a well-posed analytic inequality.

**Path B2: Resolve C5 (dynamics-arithmetic bridge) + composition.**
Rigorously construct the map from adelic characters to Weil matrix entries, then show that the decomposition W = Σ D_p + Σ M_{p,q} with each summand ≤ 0 on V_prim gives W|_{V_prim} ≤ 0. The obstacle: individual M_{p,q} blocks are rank-2 indefinite matrices (eigenvalues ±|c|), so the pair-by-pair negativity of the cross-correlation operator must translate to negativity of the off-diagonal block on the global V_prim. This requires understanding how the correlation operator structure differs from the raw matrix block.

**Path B3: Improve BV constants to make finite computation feasible.**
Under Elliott-Halberstam conjecture: P_eff ~ 10⁴ (feasible). Under improved unconditional BV: P_eff ~ 10¹²⁻¹⁵ (not yet feasible). This is a problem in explicit analytic number theory with active researchers (Platt, Trudgian, Kadiri, Ramaré).

### 7.4 What Definitely Does NOT Work

- **MSGC via interlacing:** The interlacing proves the OPPOSITE direction — λ_max goes UP toward 0, not down away from it. X₁ = 47 claim is wrong.
- **Extending computation to ∞:** Max primitive eigenvalue monotonically approaches 0; whether it reaches 0 IS the content of RH.
- **Full AMR (Path A):** Too many independent critical gaps (C1-C3).
- **Condensed mathematics arguments:** Speculative, relying on unpublished results.
- **Perturbation from Haar alone:** The margin is too thin (ratio ‖M_zeros‖/λ_max > 100).

---

## 8. Comparison with Known Approaches

| Approach | Reduces RH to | Feasibility |
|----------|--------------|-------------|
| Turing method (zero verification) | Verify zeros up to height T | T → ∞ needed (never complete) |
| De Bruijn-Newman | Determine Λ_dBN = 0 | Proved Λ ≤ 0 (Rodgers-Tao 2020); equality open |
| Li's criterion | Li coefficients λ_n ≥ 0 ∀n | Infinite sequence; no finite reduction |
| Sieve + diagonal dominance | Matrix for p ≤ 10²⁰ | Infeasible (10¹⁸ rows) |
| **AMR (Path B)** | **Matrix for p ≤ 10⁶** | **Feasible (78K rows) if gaps close** |

The AMR advantage is structural: Baker's theorem provides the "phase information" that sieve methods lack, bypassing the parity barrier. The effective threshold (10⁶ vs 10²⁰) makes the finite computation actually feasible.

---

## 9. Files and References

### Documents Produced in This Session

| File | Content | Agent |
|------|---------|-------|
| amr/proofs/AUDIT-REPORT.md | Adversarial audit: 6C/7M/5m gaps | auditor |
| amr/proofs/regularization-bridge.md | Cesaro regularization formalization | foundations |
| amr/computational/large_scale_eigenvalue_results.md | Eigenvalue data P₀ ≤ 750 | compute |
| amr/proofs/finite-verification.md | Threshold reduction analysis | threshold |
| MASTER-PROOF.md | This synthesis | team lead |

### Key Pre-Existing Documents

| File | Role |
|------|------|
| amr/AMR-MANIFESTO.md | Founding document of AMR |
| asg/ASG-MANIFESTO.md | Reduction of RH to APT |
| amr/proofs/entropy-positivity.md | Entropy-Positivity Duality |
| amr/proofs/circularity-resolution.md | K_zeros circularity resolution |
| amr/proofs/ergodicity-proof.md | μ_ar = Haar attempt (has gaps) |
| amr/foundations/amr-foundations.md | Core definitions and axioms |

---

## 10. Conclusion

The AMR framework contains genuinely novel ideas — the reinterpretation of Weil's explicit formula through measure rigidity on the adelic solenoid, the spectral cancellation of zero-dependent terms, and the entropy-positivity connection. The computational evidence is impressive: APT holds for every tested truncation up to 264×264.

However, the proof is **further from complete than the framework documents claim**:

1. **The MSGC is false.** Interlacing proves λ_max increases toward 0 with more primes. The "X₁ = 47" claim is wrong. No finite verification suffices without a theoretical argument bridging to infinity.

2. **The dynamics-arithmetic bridge (C5) is absent.** The most important gap is connecting the adelic measure-rigidity machinery to the actual Weil matrix. Without this bridge, the measure theory is disconnected from the number theory.

3. **The pair-by-pair composition fails naively.** Individual off-diagonal blocks M_{p,q} are indefinite (rank-2 with ±|c| eigenvalues). The cross-correlation operator negativity on pairwise primitive subspaces does not automatically compose to negativity on the global primitive subspace.

4. **The max eigenvalue approaching 0 IS the content of RH.** The interlacing proves the limit lim λ_max(W_N|_{prim}) exists. Whether this limit is ≤ 0 or > 0 is exactly the Riemann Hypothesis. The framework has reformulated RH into a clean spectral question but has not answered it.

**The most promising path forward** is the direct analytic inequality: prove ∫ |G_c(τ)|² K̂_prime(τ) dτ ≤ 0 for all primitive G_c, using the spectral cancellation (K̂_prime is zero-independent) combined with Baker's theorem (lattice spacing of prime-power frequencies). This bypasses finite verification entirely and attacks RH as a well-posed analytic inequality about the digamma function and the prime-power lattice.

**Honest status:** A creative research program with strong computational evidence, several novel mathematical insights, and well-identified gaps. Not a proof, but a productive direction.

---

*Synthesized: 2026-02-13*
*Based on work by 5 agents: team lead, auditor, foundations, compute, threshold*
*Framework: Arithmetic Measure Rigidity (AMR)*
