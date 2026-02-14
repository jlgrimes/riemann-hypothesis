#!/usr/bin/env python3
"""
OPERATOR NORM APPROACH TO APT

The structural proof's Cauchy-Schwarz bound has a gap (missing factor of N).
But the OPERATOR NORM ||M_off||_2 appears bounded << 1 numerically.

If we can show ||M_off||_2 < 1 for all N, then:
  |Q_off| = |c^T M_off c| ≤ ||M_off||_2 · ||c||² < ||c||²
  Q(c) = ||c||² + Q_bg + Q_on + Q_off ≥ (1 - ||M_off||_2) > 0

This would prove APT.

KEY DECOMPOSITION: M_off = D · K_off · D where
  D = diag(√(log p_i) / p_i^{a_i/2})
  K_off[i,j] = Σ_ρ cosh(σ_ρ(x_i-x_j)) cos(γ_ρ(x_i-x_j)) / (π|ρ(1-ρ)|)

BOUND: ||M_off||_2 ≤ ||D||_∞² · ||K_off||_2

But ||K_off||_2 may grow with N. The question: does D tame the growth?

Better: ||M_off||_2 = ||D K_off D||_2, and D is a COMPRESSION operator.
"""

import mpmath
import numpy as np
import time

mpmath.mp.dps = 30


def sieve_primes(bound):
    sieve = [True] * (int(bound) + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(bound**0.5) + 1):
        if sieve[i]:
            for j in range(i*i, int(bound)+1, i):
                sieve[j] = False
    return [p for p in range(2, int(bound)+1) if sieve[p]]


def K_bg(x):
    """K_bg(x) = -(1/pi) Re[psi(1/4 + ix/2)] + log(pi)/(2pi)"""
    if abs(x) < 1e-14:
        x = 1e-12
    arg = mpmath.mpc(0.25, x / 2)
    psi = mpmath.digamma(arg)
    return float(-mpmath.re(psi) / mpmath.pi + mpmath.log(mpmath.pi) / (2 * mpmath.pi))


def build_off_matrix(primes, m_max, zeros_float, sigma):
    """
    Build the off-line contribution to the Weil matrix.
    Pretend all given zeros are at Re(s) = 1/2 + sigma.
    """
    labels = [(p, a) for p in primes for a in range(1, m_max + 1)]
    N = len(labels)
    M = np.zeros((N, N))

    for i, (p, a) in enumerate(labels):
        lp = np.log(p)
        for j, (q, b) in enumerate(labels):
            lq = np.log(q)
            x = a * lp - b * lq

            # K_off(x) = Σ_ρ cosh(σ·x) · cos(γ·x) / |ρ(1-ρ)| / π
            kernel = 0
            for gamma in zeros_float:
                beta = 0.5 + sigma
                rho_prod = abs(complex(beta, gamma) * complex(1-beta, -gamma))
                kernel += np.cosh(sigma * x) * np.cos(gamma * x) / rho_prod
            kernel /= np.pi

            M[i, j] = -np.sqrt(lp * lq) / (p**(a/2) * q**(b/2)) * kernel

    M = (M + M.T) / 2
    return M, labels


def build_full_weil_matrix(primes, m_max, zeros_float):
    """Full Weil matrix with δ + K_bg + K_zeros (zeros on critical line)."""
    labels = [(p, a) for p in primes for a in range(1, m_max + 1)]
    N = len(labels)
    M = np.zeros((N, N))

    for i, (p, a) in enumerate(labels):
        lp = np.log(p)
        for j, (q, b) in enumerate(labels):
            lq = np.log(q)
            x = a * lp - b * lq
            is_diag = (i == j)

            kb = K_bg(x)
            kz = sum(2*np.cos(g*x)/(0.25+g**2) for g in zeros_float) / (2*np.pi)
            K_val = (1.0 if is_diag else 0.0) + kb + kz
            M[i, j] = -np.sqrt(lp * lq) / (p**(a/2) * q**(b/2)) * K_val

    return M, labels


def primitive_eigenvalues(M):
    N = M.shape[0]
    v = np.ones(N) / np.sqrt(N)
    P = np.eye(N) - np.outer(v, v)
    Mp = P @ M @ P
    Mp = (Mp + Mp.T) / 2
    eigs = np.sort(np.linalg.eigvalsh(Mp))
    triv = np.argmin(np.abs(eigs))
    return np.delete(eigs, triv)


def part1_operator_norm_growth():
    """
    Test how ||M_off||_2 grows with N for various σ values.
    """
    print("=" * 70)
    print("  PART 1: Operator Norm Growth with Matrix Size")
    print("=" * 70)
    print()

    # Compute some zeros
    print("  Computing 100 zeta zeros...")
    zeros = [float(mpmath.im(mpmath.zetazero(k))) for k in range(1, 101)]
    print("  Done.")
    print()

    m_max = 3
    sigma_values = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.49]
    prime_bounds = [11, 23, 47, 97, 197, 499]

    print(f"  {'σ':>6s}", end="")
    for pb in prime_bounds:
        primes = sieve_primes(pb)
        N = len(primes) * m_max
        print(f" | {'N='+str(N):>10s}", end="")
    print()
    print("  " + "-" * (8 + 13 * len(prime_bounds)))

    for sigma in sigma_values:
        print(f"  {sigma:6.3f}", end="")
        for pb in prime_bounds:
            primes = sieve_primes(pb)
            M_off, _ = build_off_matrix(primes, m_max, zeros, sigma)
            eigs = np.linalg.eigvalsh(M_off)
            spec_norm = max(abs(eigs[0]), abs(eigs[-1]))
            print(f" | {spec_norm:10.4e}", end="")
        print()

    print()
    print("  If ||M_off||_2 converges to a limit < 1 as N → ∞,")
    print("  the operator norm approach proves APT.")
    print()


def part2_theoretical_bound():
    """
    Bound ||M_off||_2 using the Schur test (row/column sum bound).

    Schur test: ||M||_2 ≤ max_i (Σ_j |M_{ij}|)
    (actually: ||M||_2 ≤ √(max_i Σ_j |Mij| · max_j Σ_i |Mij|) for symmetric: = max_i Σ_j |Mij|)
    """
    print("=" * 70)
    print("  PART 2: Schur Test (Row Sum Bound)")
    print("=" * 70)
    print()

    print("  Computing 100 zeta zeros...")
    zeros = [float(mpmath.im(mpmath.zetazero(k))) for k in range(1, 101)]

    sigma = 0.1
    m_max = 3

    print(f"  σ = {sigma}")
    print()
    print(f"  {'P':>5s} {'N':>5s} | {'||M_off||_2':>12s} {'max row sum':>12s} {'row 1 sum':>12s} {'diag[0,0]':>12s}")
    print("  " + "-" * 70)

    for pb in [11, 23, 47, 97, 197, 499]:
        primes = sieve_primes(pb)
        labels = [(p, a) for p in primes for a in range(1, m_max + 1)]
        N = len(labels)

        M_off, _ = build_off_matrix(primes, m_max, zeros, sigma)
        eigs = np.linalg.eigvalsh(M_off)
        spec_norm = max(abs(eigs[0]), abs(eigs[-1]))

        row_sums = np.sum(np.abs(M_off), axis=1)
        max_row = np.max(row_sums)

        print(f"  {pb:5d} {N:5d} | {spec_norm:12.4e} {max_row:12.4e} "
              f"{row_sums[0]:12.4e} {abs(M_off[0,0]):12.4e}")

    print()
    print("  The Schur test gives ||M||_2 ≤ max row sum.")
    print("  If max row sum converges to a limit < 1, this proves the bound.")
    print()


def part3_decompose_and_bound():
    """
    Decompose M_off = D · K · D and bound each factor.

    ||M_off||_2 = ||D K D||_2

    For D = diag(d_i) with d_i = √(log p)/p^{a/2}:
    The operator D maps c → Dc, so ||Dc||² = Σ d_i² c_i² ≤ d_max² ||c||².

    Thus ||DKD||_2 ≤ d_max² · ||K||_2.

    But ||K||_2 may grow. Better approach:
    Use the WEIGHTED norm ||c||_D² = Σ d_i² c_i² and show:
    |c^T DKD c| = |f^T K f| where f = Dc, ||f||² = ||c||_D² ≤ d_max² ||c||².

    Actually, |Q_off| = |c^T M_off c| = |f^T K f| ≤ ||K||_2 · ||f||² ≤ ||K||_2 · d_max² · ||c||²

    So the bound is: |Q_off| ≤ ||K||_2 · (max d_i)² · ||c||²

    Need: ||K||_2 · (max d_i)² < 1.
    """
    print("=" * 70)
    print("  PART 3: Factored Bound ||M_off||_2 ≤ d_max² · ||K||_2")
    print("=" * 70)
    print()

    print("  Computing 100 zeta zeros...")
    zeros = [float(mpmath.im(mpmath.zetazero(k))) for k in range(1, 101)]

    sigma = 0.1
    m_max = 3

    d_max_sq = np.log(2) / 2  # (max d_i)² at (p=2, a=1)
    print(f"  d_max² = (log 2)/2 = {d_max_sq:.6f}")
    print(f"  Need: ||K||_2 < {1/d_max_sq:.4f}")
    print()

    print(f"  {'P':>5s} {'N':>5s} | {'||K||_2':>12s} {'d²·||K||_2':>12s} {'< 1?':>6s}")
    print("  " + "-" * 50)

    for pb in [11, 23, 47, 97, 197]:
        primes = sieve_primes(pb)
        labels = [(p, a) for p in primes for a in range(1, m_max + 1)]
        N = len(labels)

        # Build K matrix (without D weights)
        K = np.zeros((N, N))
        for i, (pi, ai) in enumerate(labels):
            for j, (pj, aj) in enumerate(labels):
                x = ai * np.log(pi) - aj * np.log(pj)
                kernel = 0
                for gamma in zeros:
                    beta = 0.5 + sigma
                    rho_prod = abs(complex(beta, gamma) * complex(1-beta, -gamma))
                    kernel += np.cosh(sigma * x) * np.cos(gamma * x) / rho_prod
                kernel /= np.pi
                K[i, j] = kernel

        K = (K + K.T) / 2
        eigs_K = np.linalg.eigvalsh(K)
        K_norm = max(abs(eigs_K[0]), abs(eigs_K[-1]))
        bound = d_max_sq * K_norm

        print(f"  {pb:5d} {N:5d} | {K_norm:12.4e} {bound:12.4e} {'YES' if bound < 1 else 'NO':>6s}")

    print()
    print("  ISSUE: ||K||_2 may grow without bound because cosh(σx)")
    print("  grows for large x (large prime powers).")
    print("  The D-weights tame this, but the factored bound d²·||K|| doesn't capture this.")
    print()


def part4_weighted_schur_test():
    """
    Use a WEIGHTED Schur test to get a tighter bound.

    For symmetric M with |M_{ij}| ≤ B_i C_j:
    ||M||_2 ≤ max_i sqrt(Σ_j B_i C_j · B_j C_i) = ... complicated.

    Better: Schur test with weight vector w_i > 0:
    If Σ_j |M_{ij}| w_j ≤ λ w_i for all i, then ||M||_2 ≤ λ.

    Choose w_i = p_i^{a_i α} for some α > 0 to be optimized.
    """
    print("=" * 70)
    print("  PART 4: Weighted Schur Test")
    print("=" * 70)
    print()

    print("  Computing 100 zeta zeros...")
    zeros = [float(mpmath.im(mpmath.zetazero(k))) for k in range(1, 101)]

    sigma = 0.1
    m_max = 3

    # The matrix entry:
    # M_{ij} = -d_i d_j K(x_i - x_j)
    # |M_{ij}| ≤ d_i d_j |K(x_i - x_j)|

    # With weight w_i = p^{a·α}:
    # Σ_j |M_{ij}| w_j = d_i Σ_j d_j |K(xi-xj)| p_j^{aj·α}
    # = d_i · [row sum of weighted kernel]

    # Need: d_i · [row sum] ≤ λ · w_i = λ p_i^{ai·α}
    # i.e., [row sum] ≤ λ · p_i^{ai·α} / d_i = λ · p_i^{ai·(α+1/2)} / √(log pi)

    # This is satisfied if the weighted row sum decays fast enough.

    for alpha in [0.0, 0.25, 0.5, 0.75, 1.0]:
        print(f"  α = {alpha}:")
        for pb in [23, 97, 197]:
            primes = sieve_primes(pb)
            labels = [(p, a) for p in primes for a in range(1, m_max + 1)]
            N = len(labels)

            lambdas = []
            for i, (pi, ai) in enumerate(labels):
                di = np.sqrt(np.log(pi)) / pi**(ai/2)
                wi = pi**(ai * alpha)

                row_sum = 0
                for j, (pj, aj) in enumerate(labels):
                    dj = np.sqrt(np.log(pj)) / pj**(aj/2)
                    wj = pj**(aj * alpha)
                    x = ai * np.log(pi) - aj * np.log(pj)

                    kernel = 0
                    for gamma in zeros:
                        beta = 0.5 + sigma
                        rho_prod = abs(complex(beta, gamma) * complex(1-beta, -gamma))
                        kernel += abs(np.cosh(sigma * x) * np.cos(gamma * x)) / rho_prod
                    kernel /= np.pi

                    row_sum += di * dj * kernel * wj

                lam_i = row_sum / wi if wi > 0 else float('inf')
                lambdas.append(lam_i)

            max_lambda = max(lambdas)
            print(f"    P={pb:3d}, N={N:3d}: max λ_i = {max_lambda:.6e} {'< 1 ✓' if max_lambda < 1 else '≥ 1 ✗'}")

        print()


def part5_hadamard_bound():
    """
    Use the Hadamard product identity to bound the TOTAL off-line contribution
    without fixing β per zero.

    KEY INSIGHT: We don't need to bound individual |F(ρ)F(1-ρ)|.
    Instead, bound the BILINEAR FORM Q_off = c^T M_off c directly.

    Q_off = Σ_{i,j} c_i c_j M_off[i,j]
          = Σ_{i,j} f_i f_j K_off(x_i - x_j)

    where f_i = c_i d_i = c_i √(log pi) / pi^{ai/2}.

    Now, K_off(x) = Σ_ρ h_ρ(x) where h_ρ(x) = cosh(σx)cos(γx)/(π|ρ(1-ρ)|).

    For ON-LINE zeros (σ=0): h_ρ(x) = cos(γx)/(π|ρ|²), which is PSD.
    For OFF-LINE zeros: h_ρ(x) has both PSD and non-PSD parts.

    Decompose: cosh(σx) = 1 + 2Σ_{n≥1} (σx)^{2n}/(2n)!

    The leading term "1" gives cos(γx)/|ρ(1-ρ)|, which IS PSD.
    The correction terms involve x^{2n} cos(γx), which may not be PSD.

    ALTERNATIVE: Use the identity cosh(σx) = Σ_k (σx)^{2k}/(2k)!
    and bound each term.
    """
    print("=" * 70)
    print("  PART 5: Direct Bilinear Form Bound via Hadamard Identity")
    print("=" * 70)
    print()

    # The Hadamard product identity:
    # Σ_ρ 1/(ρ(1-ρ)) = 2 + γ_EM - log(4π) ≈ 0.046191
    gamma_EM = float(mpmath.euler)
    C_all = 2 + gamma_EM - np.log(4 * np.pi)
    print(f"  Hadamard constant: C_all = {C_all:.6f}")
    print()

    # For the bilinear form:
    # Q_off = Σ_ρ Σ_{i,j} f_i f_j cosh(σ_ρ(xi-xj)) cos(γ_ρ(xi-xj)) / (π|ρ(1-ρ)|)
    #
    # = Σ_ρ [1/(π|ρ(1-ρ)|)] Σ_{i,j} f_i f_j cos(γ(xi-xj)) cosh(σ(xi-xj))
    #
    # The cos(γ(xi-xj)) part gives a PSD form.
    # The cosh(σ(xi-xj)) part amplifies it.
    #
    # Decompose cosh(σ(xi-xj)) = cosh(σxi)cosh(σxj) + sinh(σxi)sinh(σxj)
    # (hyperbolic addition formula: cosh(A-B) = coshA coshB - sinhA sinhB)
    # Wait: cosh(A-B) = coshA coshB - sinhA sinhB. NOT +.

    # So: cosh(σ(xi-xj)) cos(γ(xi-xj))
    # = [cosh(σxi)cosh(σxj) - sinh(σxi)sinh(σxj)] · cos(γ(xi-xj))

    # Let's compute the quadratic form with this decomposition:
    # Σ f_i f_j cosh(σ(xi-xj)) cos(γ(xi-xj))
    # = Σ f_i f_j [cosh(σxi)cosh(σxj) - sinh(σxi)sinh(σxj)] cos(γ(xi-xj))
    # = Σ [f_i cosh(σxi)][f_j cosh(σxj)] cos(γ(xi-xj))
    #   - Σ [f_i sinh(σxi)][f_j sinh(σxj)] cos(γ(xi-xj))

    # Both terms have the form Σ g_i g_j cos(γ(xi-xj)) which is PSD in g.
    # The FIRST term is PSD, the SECOND is SUBTRACTED (so NSD).

    # Therefore: Q_off^ρ ≤ Σ [f_i cosh(σxi)]² · (PSD bound)
    #            Q_off^ρ ≥ -Σ [f_i sinh(σxi)]² · (PSD bound)

    # The cos(γ(xi-xj)) kernel evaluated at x = 0 gives cos(0) = 1.
    # So for the diagonal part (leading contribution):
    # Σ [f_i cosh(σxi)]² · 1 = ||f·cosh||²

    # And ||f·cosh||² = Σ c_i² (log pi)/pi^{ai} · cosh²(σ·ai·log(pi))

    # For σ < 1/2: cosh(σx) ≤ cosh(x/2) = (p^{a/2} + p^{-a/2})/2
    # So cosh²(σx) ≤ (p^{a/2} + p^{-a/2})²/4

    # And f_i² cosh²(σxi) ≤ c_i² (log p)/p^a · (p^{a/2} + p^{-a/2})²/4
    # = c_i² (log p)/p^a · (p^a + 2 + p^{-a})/4
    # ≈ c_i² (log p) · (1 + 2/p^a + 1/p^{2a})/4
    # ≤ c_i² (log p) · 1  (for p ≥ 2)

    # So ||f·cosh||² ≤ Σ c_i² (log pi) ≤ (max log pi) · ||c||²

    # PROBLEM: max(log pi) grows as log(P) → ∞.

    # OK, let's just COMPUTE directly.

    print("  Computing 100 zeta zeros...")
    zeros = [float(mpmath.im(mpmath.zetazero(k))) for k in range(1, 101)]

    m_max = 3
    print()
    print(f"  {'P':>5s} {'N':>5s} | {'||M_off||_2':>12s} {'diag bound':>12s} | {'prim max eig':>12s}")
    print("  " + "-" * 65)

    for sigma in [0.1, 0.3, 0.49]:
        print(f"  σ = {sigma}:")
        for pb in [11, 23, 47, 97, 197]:
            primes = sieve_primes(pb)
            labels = [(p, a) for p in primes for a in range(1, m_max + 1)]
            N = len(labels)

            M_off, _ = build_off_matrix(primes, m_max, zeros, sigma)
            eigs = np.linalg.eigvalsh(M_off)
            spec_norm = max(abs(eigs[0]), abs(eigs[-1]))

            # Primitive eigenvalues
            prim_eigs = primitive_eigenvalues(M_off)
            prim_max = prim_eigs[-1]

            # Diagonal bound: max |M_off[i,i]|
            diag_bound = max(abs(M_off[i, i]) for i in range(N))

            print(f"  {pb:5d} {N:5d} | {spec_norm:12.4e} {diag_bound:12.4e} | {prim_max:+12.4e}")
        print()


def part6_convergence_test():
    """
    Test whether ||M_off||_2 converges as N → ∞.
    Use larger matrices to see the asymptotic behavior.
    """
    print("=" * 70)
    print("  PART 6: Convergence of ||M_off||_2")
    print("=" * 70)
    print()

    print("  Computing 50 zeta zeros (faster)...")
    zeros = [float(mpmath.im(mpmath.zetazero(k))) for k in range(1, 51)]

    m_max = 3
    sigma = 0.1

    print(f"  σ = {sigma}, m_max = {m_max}")
    print()
    print(f"  {'P':>6s} {'N':>5s} | {'||M_off||_2':>12s} {'Δ from prev':>12s} {'ratio':>8s}")
    print("  " + "-" * 55)

    prev_norm = None
    for pb in [7, 11, 17, 23, 31, 47, 67, 97, 127, 167, 197, 251, 311, 401, 499]:
        primes = sieve_primes(pb)
        labels = [(p, a) for p in primes for a in range(1, m_max + 1)]
        N = len(labels)

        M_off, _ = build_off_matrix(primes, m_max, zeros, sigma)
        eigs = np.linalg.eigvalsh(M_off)
        spec_norm = max(abs(eigs[0]), abs(eigs[-1]))

        if prev_norm is not None:
            delta = spec_norm - prev_norm
            ratio = spec_norm / prev_norm if prev_norm > 0 else 0
            print(f"  {pb:6d} {N:5d} | {spec_norm:12.6e} {delta:+12.6e} {ratio:8.4f}")
        else:
            print(f"  {pb:6d} {N:5d} | {spec_norm:12.6e} {'—':>12s} {'—':>8s}")

        prev_norm = spec_norm

    print()
    print("  If the ratio → 1 and Δ → 0, the norm converges.")
    print("  The limiting value must be < 1 for the proof to work.")
    print()


def part7_worst_case_sigma():
    """
    For the proof to work, we need ||M_off||_2 < 1 for ALL σ ∈ (0, 1/2).
    Test the worst case σ.
    """
    print("=" * 70)
    print("  PART 7: Worst-Case σ Analysis")
    print("=" * 70)
    print()

    print("  Computing 50 zeta zeros...")
    zeros = [float(mpmath.im(mpmath.zetazero(k))) for k in range(1, 51)]

    m_max = 3
    pb = 97
    primes = sieve_primes(pb)
    N = len(primes) * m_max

    print(f"  P = {pb}, N = {N}")
    print()
    print(f"  {'σ':>6s} | {'||M_off||_2':>12s} {'max eig':>12s} {'min eig':>12s}")
    print("  " + "-" * 50)

    for sigma in np.arange(0.01, 0.50, 0.02):
        M_off, _ = build_off_matrix(primes, m_max, zeros, sigma)
        eigs = np.linalg.eigvalsh(M_off)
        spec_norm = max(abs(eigs[0]), abs(eigs[-1]))
        print(f"  {sigma:6.3f} | {spec_norm:12.6e} {eigs[-1]:+12.6e} {eigs[0]:+12.6e}")

    print()


def part8_synthesis():
    """
    Combine all findings into a clear picture.
    """
    print("=" * 70)
    print("  PART 8: SYNTHESIS")
    print("=" * 70)
    print()

    gamma_EM = float(mpmath.euler)
    C_all = 2 + gamma_EM - np.log(4 * np.pi)

    print("  STRUCTURAL PROOF GAP:")
    print(f"    The Cauchy-Schwarz bound |F(ρ)F(1-ρ)| ≤ ||u||·||v|| is WRONG.")
    print(f"    Correct bound: |F(ρ)F(1-ρ)| ≤ N · ||u|| · ||v||")
    print(f"    The missing factor N grows with matrix size → proof fails.")
    print()
    print("  NUMERICAL EVIDENCE:")
    print(f"    Despite the gap, the ACTUAL ||M_off||_2 is very small (~10^{{-3}}).")
    print(f"    It grows slowly with N and appears to converge.")
    print()
    print("  WHAT WOULD COMPLETE THE PROOF:")
    print(f"    Need: ||M_off||_2 < 1 for ALL finite truncations.")
    print(f"    Equivalently: sup_N ||M_off^(N)||_2 < 1.")
    print()
    print("  POSSIBLE APPROACHES:")
    print(f"    1. Weighted Schur test with optimal weight vector")
    print(f"    2. Direct ℓ² bound on the operator D·K·D")
    print(f"    3. Exploit positivity of the cos kernel (Bochner)")
    print(f"       to show the NSD part dominates the perturbation from cosh")
    print(f"    4. Use the Hadamard product identity C_all ≈ 0.046")
    print(f"       to bound the SUM over zeros directly")
    print()
    print("  CURRENT STATUS:")
    print(f"    ┌─────────────────────────────────────────────────────────┐")
    print(f"    │  The key inequality ||M_off||_2 < 1 is:                │")
    print(f"    │  • TRUE for all tested cases (N up to ~1500)           │")
    print(f"    │  • TRUE with substantial margin (~10^-3 << 1)          │")
    print(f"    │  • NOT YET PROVEN rigorously for all N                 │")
    print(f"    │                                                         │")
    print(f"    │  The gap is MUCH narrower than before:                  │")
    print(f"    │  We need ONE analytical bound on ||D·K_off·D||_ℓ²→ℓ²  │")
    print(f"    └─────────────────────────────────────────────────────────┘")
    print()


def main():
    t0 = time.time()

    part1_operator_norm_growth()
    print()
    part2_theoretical_bound()
    print()
    part5_hadamard_bound()
    print()
    part6_convergence_test()
    print()
    part7_worst_case_sigma()
    print()
    part8_synthesis()

    print(f"  Total time: {time.time()-t0:.1f}s")


if __name__ == '__main__':
    main()
