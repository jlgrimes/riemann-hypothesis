# Summary: Structural Proof Attempt for the Riemann Hypothesis

## What Was Attempted

A proof of RH via the Arithmetic Positivity Theorem (APT):
- RH ⟺ All primitive eigenvalues of the Weil matrix are ≤ 0
- Decompose: Q(c) = ||f||² + Q_bg + Q_on + Q_off
- Show each component has the right sign

## What Is Rigorously Proven

1. **Weil explicit formula decomposition**: Q = ||f||² + Q_bg + Q_on + Q_off ✓
2. **||f||² > 0** for f = Dc, c ≠ 0 ✓
3. **Q_bg ≥ 0** on primitive (Lorentzian proof: K_bg is PSD) ✓
4. **Q_on ≥ 0** (Bochner's theorem: cos kernel is PSD) ✓
5. **Hadamard identity**: Σ_ρ 1/(ρ(1-ρ)) = 2 + γ_EM - log(4π) ≈ 0.046 ✓
6. **Contraction property**: D-weights tame cosh growth for σ < 1/2 ✓

## The Gap

The structural proof claimed: |F(ρ)·F(1-ρ)| ≤ ||u||·||v|| (Cauchy-Schwarz)

**This is wrong.** The correct bound is:
```
|F(ρ)·F(1-ρ)| ≤ N · ||u|| · ||v||
```
The missing factor of N (matrix dimension) makes the bound grow with matrix size.

**Verified numerically**: the bound is violated by factors up to 5x for random primitive vectors.

## What Computation Reveals

### Operator norm of M_off
- ||M_off||_2 < 0.016 for all N ≤ 285, σ ≤ 0.49
- Growth rate: approximately N^{0.17} for σ = 0.1
- 60x margin below the critical threshold of 1

### Critical σ* (the decisive finding)
For each matrix size N, there exists σ*(N) such that:
- σ < σ*(N): APT holds (all zeros could be off-line at σ)
- σ > σ*(N): APT fails (off-line zeros at σ break positivity)

| N | σ*(N) | σ*·log(N) |
|---|-------|-----------|
| 45 | 0.500 | 1.90 |
| 75 | 0.478 | 2.06 |
| 135 | 0.376 | 1.84 |

**Scaling: σ*(N) ∼ C / (log N)^{1.1}**

This means σ*(N) → 0 as N → ∞. APT cannot hold for fixed σ > 0 as N grows.

### Interpretation
The APT approach reproduces the **classical zero-free region** (σ < c/log T), but cannot prove RH. The off-line perturbation from zeros at σ > 0 eventually overcomes the positive floor ||f||² for large enough matrices.

## Framework for Future Attempts

The obstruction is fundamentally that:
1. For on-line zeros (σ = 0): cos(γx) is PSD → eigenvalues ≤ 0 ✓
2. For off-line zeros (σ > 0): cosh(σx)·cos(γx) has positive eigenvalues that grow with the kernel matrix size
3. The D-weights suppress this growth but not enough for σ > C/log(N)

To prove RH via this path, one would need to show that the ACTUAL zeros of ζ(s) (which are constrained by the functional equation, not arbitrary) cannot conspire to create positive eigenvalues. This requires number-theoretic input beyond what the structural/analytic framework provides alone.

## Files

- `structural_proof.py` — The proof attempt (gap at Step 3)
- `audit_cauchy_schwarz.py` — Demonstrates the N-factor gap
- `operator_norm_proof.py` — ||M_off||_2 growth analysis
- `kernel_eigenvalue_test.py` — Definitive APT test with off-line zeros
- `critical_sigma_test.py` — σ*(N) computation
- `lorentzian_proof.py` — K_bg PSD proof (valid)
- `bg_only_test.py` — Background matrix tests
- `multiplicative_fourier_proof.py` — Spectral domain framework
- `aqft_proof.py` — Arithmetic QFT framework
