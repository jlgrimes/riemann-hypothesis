# Learnings: Failed & Circular Proof Approaches

These scripts represent proof strategies that were explored and found to be either circular (equivalent to restating RH) or incomplete (hitting irreducible gaps). They're preserved here as documentation of what was tried.

## Spectral Weight Approach (Circular)

**Files:** `weil_positivity_proof_attempt.py`, `spectral_weight_profile.py`, `eigenvector_anatomy.py`, `asymptotic_eigenvalue_fit.py`, `bg_vs_zeros_separation.py` + their `*_results.md`

**Strategy:** Decompose the Weil kernel K = K_bg + K_zeros, show K_hat_cont(tau) > 0 via Lorentzian/Bochner, show zero contributions are positive, conclude APT.

**What worked:**
- K_hat_cont(tau) > 1 for all tau > 0 (unconditionally proven)
- Background dominates zeros by ~400x at all tested scales
- Zero-independence of K_hat_prime (proven)
- Eigenvectors concentrate on large primes at highest m

**Why it failed:** The final assembly step -- proving K_hat_prime(tau) >= 0 as a distribution -- is equivalent to RH by the Weil criterion. The decomposition reformulates RH as RH.

## Diagonal Dominance (Insufficient)

**File:** `diagonal_dominance_proof.py`

**Strategy:** Show the Weil matrix is diagonally dominant, implying NSD.

**Why it failed:** Off-diagonal entries decay too slowly; diagonal dominance condition fails for large matrices.

## Structural Proof / Cauchy-Schwarz (Fatal Gap)

**Files:** `structural_proof.py`, `audit_cauchy_schwarz.py`

**Strategy:** Bound cross-terms via Cauchy-Schwarz inequality.

**Why it failed:** Missing N-factor in the Cauchy-Schwarz bound. The audit script (`audit_cauchy_schwarz.py`) explicitly identifies this fatal gap.

## Operator Norm Bound (Incomplete)

**File:** `operator_norm_proof.py`

**Strategy:** Show ||M_off||_2 < 1 to guarantee M is NSD.

**Why it failed:** Identifies the gap but cannot close it; the operator norm bound is not tight enough.

## AQFT Framework (Exploratory)

**File:** `aqft_proof.py`

**Strategy:** Use arithmetic QFT framework (Kallen-Lehmann representation, partition functions).

**Why it failed:** Novel formulation but no closure; exploratory without rigorous path to completion.

## Multiplicative Fourier (Incomplete)

**File:** `multiplicative_fourier_proof.py`

**Strategy:** Spectral domain approach with closed-form Phi_bg.

**Why it failed:** Hits spectral floor limit; cannot push past the critical threshold.

## Li Coefficient Analysis (Reformulation)

**File:** `li_coefficient_analysis.py`

**Strategy:** Analyze Li coefficients per-zero as equivalent RH criterion.

**Why it failed:** Valid equivalent formulation but doesn't advance main proof; just another way to state RH.

## Other Exploratory Scripts

- `spectral_screening.py` — "Arithmetic Kernel Screening" concept; interesting but incomplete
- `pair_prime_correlation.py` — Empirical pair correlation data; didn't connect to proof
- `jensen_rh.py` — Jensen polynomial hyperbolicity testing; valid formulation but orthogonal to ASG framework

## Key Lesson

Every approach that tries to prove the infinite-dimensional APT directly through spectral analysis of the Weil kernel eventually reduces to proving K_hat_prime >= 0, which IS RH. The finite certification (certified_weil_matrix.py) works because it sidesteps this by computing explicit matrices. The viable path forward is either (a) measure rigidity (AMR), (b) extending finite certification with explicit tail bounds, or (c) finding a monotonicity principle that bootstraps finite verification to infinity.
