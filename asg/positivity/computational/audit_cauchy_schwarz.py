#!/usr/bin/env python3
"""
CRITICAL AUDIT: Is |F(ρ)·F(1-ρ)| ≤ ||u||·||v|| correct?

The structural proof claims:
    F(ρ) = Σ u_i e^{-iγx_i}
    F(1-ρ) = Σ v_i e^{+iγx_i}
    |F(ρ)·F(1-ρ)| ≤ ||u|| · ||v||   (???)

Standard Cauchy-Schwarz gives:
    |F(ρ)| = |<u, φ>| ≤ ||u|| · ||φ|| = ||u|| · √N
    |F(1-ρ)| = |<v, φ̄>| ≤ ||v|| · √N

So the correct bound is: |F(ρ)·F(1-ρ)| ≤ N · ||u|| · ||v||

The factor of N is FATAL: the bound grows with matrix dimension.

This script tests this definitively.
"""

import mpmath
import numpy as np

mpmath.mp.dps = 30


def sieve_primes(bound):
    sieve = [True] * (int(bound) + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(bound**0.5) + 1):
        if sieve[i]:
            for j in range(i*i, int(bound)+1, i):
                sieve[j] = False
    return [p for p in range(2, int(bound)+1) if sieve[p]]


def test_cauchy_schwarz_gap():
    """
    Compute |F(ρ)F(1-ρ)| and compare with ||u||·||v|| and N·||u||·||v||.
    """
    print("=" * 70)
    print("  AUDIT: Cauchy-Schwarz bound on |F(ρ)·F(1-ρ)|")
    print("=" * 70)
    print()

    gamma1 = float(mpmath.im(mpmath.zetazero(1)))  # 14.1347...
    beta = 0.6  # hypothetical off-line zero at σ = 0.1
    m_max = 3

    print(f"  Test zero: ρ = {beta} + {gamma1:.4f}i")
    print()

    prime_bounds = [11, 23, 47, 97, 197, 499]

    print(f"  {'P':>5s} {'N':>5s} | {'|F(ρ)F(1-ρ)|':>14s} {'||u||·||v||':>14s} "
          f"{'N·||u||·||v||':>14s} {'ratio/||u||||v||':>16s} {'exceeds?':>10s}")
    print("  " + "-" * 100)

    for pb in prime_bounds:
        primes = sieve_primes(pb)
        labels = [(p, a) for p in primes for a in range(1, m_max + 1)]
        N = len(labels)

        # Worst-case c: primitive, concentrated on small primes
        c = np.zeros(N)
        c[0] = 1/np.sqrt(2)   # (2,1)
        c[1] = -1/np.sqrt(2)  # (2,2)

        # Compute u, v, F(ρ), F(1-ρ)
        u = np.zeros(N)
        v = np.zeros(N)
        for i, (p, a) in enumerate(labels):
            x = a * np.log(p)
            d = c[i] * np.sqrt(np.log(p)) * p**(-a/2)
            u[i] = d * np.exp(-beta * x)
            v[i] = d * np.exp(-(1-beta) * x)

        # F(ρ) = Σ u_i e^{-iγx_i}
        phases = np.array([np.exp(-1j * gamma1 * a * np.log(p)) for p, a in labels])
        F_rho = np.sum(u * phases)
        F_one_minus_rho = np.sum(v * np.conj(phases))

        product = abs(F_rho * F_one_minus_rho)
        norm_u = np.linalg.norm(u)
        norm_v = np.linalg.norm(v)
        uv = norm_u * norm_v
        Nuv = N * uv

        ratio = product / uv if uv > 0 else 0
        exceeds = "YES!" if product > uv * 1.001 else "no"

        print(f"  {pb:5d} {N:5d} | {product:14.6e} {uv:14.6e} "
              f"{Nuv:14.6e} {ratio:16.6f} {exceeds:>10s}")

    print()
    print("  If ratio > 1: the bound ||u||·||v|| is VIOLATED.")
    print("  The correct bound is N·||u||·||v|| (factor of N from Cauchy-Schwarz).")
    print()

    # Now test with γ → 0 (worst case for oscillation)
    print("=" * 70)
    print("  TEST 2: Sweep over γ to find worst case")
    print("=" * 70)
    print()

    primes = sieve_primes(47)
    labels = [(p, a) for p in primes for a in range(1, m_max + 1)]
    N = len(labels)

    # Use c = e_1 (concentrate on single entry) for clearest test
    c = np.zeros(N)
    c[0] = 1.0  # single entry (2,1)

    u = np.zeros(N)
    v = np.zeros(N)
    for i, (p, a) in enumerate(labels):
        x = a * np.log(p)
        d = c[i] * np.sqrt(np.log(p)) * p**(-a/2)
        u[i] = d * np.exp(-beta * x)
        v[i] = d * np.exp(-(1-beta) * x)

    norm_u = np.linalg.norm(u)
    norm_v = np.linalg.norm(v)
    uv = norm_u * norm_v

    print(f"  N = {N}, c = e_1 (single entry at (2,1))")
    print(f"  ||u|| = {norm_u:.6e}, ||v|| = {norm_v:.6e}, ||u||·||v|| = {uv:.6e}")
    print()
    print(f"  {'γ':>12s} {'|F(ρ)F(1-ρ)|':>14s} {'||u||·||v||':>14s} {'ratio':>10s}")
    print("  " + "-" * 55)

    max_ratio = 0
    for gamma in [0.1, 1.0, 5.0, 14.135, 21.022, 50.0, 100.0]:
        phases = np.array([np.exp(-1j * gamma * a * np.log(p)) for p, a in labels])
        F_rho = np.sum(u * phases)
        F_1mr = np.sum(v * np.conj(phases))
        product = abs(F_rho * F_1mr)
        ratio = product / uv if uv > 0 else 0
        max_ratio = max(max_ratio, ratio)
        print(f"  {gamma:12.3f} {product:14.6e} {uv:14.6e} {ratio:10.6f}")

    print()
    print(f"  Max ratio over tested γ: {max_ratio:.6f}")
    print()

    # Test 3: with a PRIMITIVE vector that has many nonzero entries
    print("=" * 70)
    print("  TEST 3: Primitive vector with all entries ≠ 0")
    print("=" * 70)
    print()

    # c = (1/√2, -1/√2, 0, ...) won't exercise many terms
    # Use a random primitive vector instead
    np.random.seed(42)

    for pb in [23, 47, 97, 197]:
        primes = sieve_primes(pb)
        labels = [(p, a) for p in primes for a in range(1, m_max + 1)]
        N = len(labels)

        c = np.random.randn(N)
        c -= np.mean(c)  # make primitive
        c /= np.linalg.norm(c)  # normalize

        u = np.zeros(N)
        v = np.zeros(N)
        for i, (p, a) in enumerate(labels):
            x = a * np.log(p)
            d = c[i] * np.sqrt(np.log(p)) * p**(-a/2)
            u[i] = d * np.exp(-beta * x)
            v[i] = d * np.exp(-(1-beta) * x)

        norm_u = np.linalg.norm(u)
        norm_v = np.linalg.norm(v)
        uv = norm_u * norm_v

        # Try MANY γ values to find the max ratio
        max_ratio = 0
        worst_gamma = 0
        for gamma in np.linspace(0.1, 200, 2000):
            phases = np.array([np.exp(-1j * gamma * a * np.log(p)) for p, a in labels])
            F_rho = np.sum(u * phases)
            F_1mr = np.sum(v * np.conj(phases))
            product = abs(F_rho * F_1mr)
            ratio = product / uv if uv > 0 else 0
            if ratio > max_ratio:
                max_ratio = ratio
                worst_gamma = gamma

        exceeds = "VIOLATED" if max_ratio > 1.001 else "ok"
        print(f"  P={pb:3d}, N={N:3d}: max |F·F|/(||u||·||v||) = {max_ratio:.4f} "
              f"at γ={worst_gamma:.1f}  [{exceeds}]")

    print()

    # Test 4: The actual Weil matrix computation vs the bound
    print("=" * 70)
    print("  TEST 4: Direct computation of Q_off vs claimed bound")
    print("=" * 70)
    print()

    # Use actual zeta zeros
    print("  Computing 50 zeta zeros...")
    zeros = [float(mpmath.im(mpmath.zetazero(k))) for k in range(1, 51)]

    for pb in [23, 47, 97]:
        primes = sieve_primes(pb)
        labels = [(p, a) for p in primes for a in range(1, m_max + 1)]
        N = len(labels)

        # Worst-case primitive c
        c = np.zeros(N)
        c[0] = 1/np.sqrt(2)
        c[1] = -1/np.sqrt(2)

        # Compute actual Q_off contribution (pretending zeros are off-line at β=0.6)
        beta_hyp = 0.6
        Q_off = 0
        bound_sum = 0

        for gamma in zeros:
            rho = complex(beta_hyp, gamma)
            one_minus_rho = complex(1-beta_hyp, -gamma)
            rho_prod = abs(rho * one_minus_rho)

            # F(ρ) = Σ d_i e^{-ρ x_i}
            F_rho_val = 0
            F_1mr_val = 0
            for i, (p, a) in enumerate(labels):
                x = a * np.log(p)
                d = c[i] * np.sqrt(np.log(p)) * np.exp(-x/2)
                F_rho_val += d * np.exp(-complex(beta_hyp, gamma) * x)
                F_1mr_val += d * np.exp(-complex(1-beta_hyp, -gamma) * x)

            Q_off += 2 * (F_rho_val * F_1mr_val).real / rho_prod

            # The claimed bound per zero
            u = np.array([c[i] * np.sqrt(np.log(p)) * np.exp(-(beta_hyp+0.5)*a*np.log(p))
                         for i, (p, a) in enumerate(labels)])
            v = np.array([c[i] * np.sqrt(np.log(p)) * np.exp(-(1.5-beta_hyp)*a*np.log(p))
                         for i, (p, a) in enumerate(labels)])
            bound_sum += 2 * np.linalg.norm(u) * np.linalg.norm(v) / rho_prod

        w_max = np.log(2) / 2
        hadamard_bound = np.log(2) * 0.046191  # The structural proof's claimed bound

        print(f"  P={pb:3d}, N={N:3d}:")
        print(f"    |Q_off| (actual, 50 zeros)    = {abs(Q_off):.6e}")
        print(f"    Σ 2||u||·||v||/|ρ(1-ρ)|       = {bound_sum:.6e}")
        print(f"    Structural proof bound          = {hadamard_bound:.6e}")
        print(f"    ||c||² = 1")
        print()

    # SYNTHESIS
    print("=" * 70)
    print("  SYNTHESIS: THE GAP IN THE STRUCTURAL PROOF")
    print("=" * 70)
    print()
    print("  The structural proof (structural_proof.py) claims:")
    print("    |F(ρ)·F(1-ρ)| ≤ ||u|| · ||v||  (lines 150-151)")
    print()
    print("  The CORRECT Cauchy-Schwarz bound is:")
    print("    |F(ρ)| = |<u, φ>| ≤ ||u|| · ||φ|| = ||u|| · √N")
    print("    |F(1-ρ)| = |<v, φ̄>| ≤ ||v|| · √N")
    print("    |F(ρ)·F(1-ρ)| ≤ N · ||u|| · ||v||")
    print()
    print("  The missing factor of N means:")
    print("    |Q_off| ≤ N · (log 2) · C_all ≈ 0.032·N")
    print("    For N > 31: bound exceeds ||c||² = 1")
    print()
    print("  HOWEVER: the numerical tests show that the ACTUAL")
    print("  |F(ρ)·F(1-ρ)| is much smaller than N·||u||·||v||.")
    print("  The N factor is pessimistic because the phases e^{iγx_i}")
    print("  produce cancellation for most γ values.")
    print()
    print("  The REAL question: is there a TIGHTER bound that avoids N?")
    print()
    print("  Approach 1: Large Sieve (sum over zeros, use spacing)")
    print("    Σ_ρ |F(ρ)|² ≤ (N + C/δ) · ||d||²  (Gallagher)")
    print("    This bounds the SUM, not individual terms.")
    print()
    print("  Approach 2: Diagonal dominance")
    print("    F(ρ)·F(1-ρ) = <u,v> + off-diagonal")
    print("    <u,v> = Σ c_i² w_i is β-independent.")
    print("    Off-diagonal cancels when summed over zeros.")
    print()
    print("  Approach 3: Operator norm bound")
    print("    ||M_off||_2 < 1 directly, using structure of the kernel.")
    print()


def test_large_sieve_approach():
    """
    Test whether the Large Sieve inequality gives a useful bound.

    Montgomery-Vaughan Large Sieve:
    Σ_r |Σ_n a_n e^{2πi n α_r}|² ≤ (N-1+δ^{-1}) Σ|a_n|²

    where δ = min_{r≠s} ||α_r - α_s||.

    For our case: α_r = γ_r · x_min / (2π), where γ_r are zeta zero heights
    and the sum is over zeros.
    """
    print()
    print("=" * 70)
    print("  LARGE SIEVE APPROACH")
    print("=" * 70)
    print()

    print("  The Large Sieve inequality bounds Σ_ρ |F(ρ)|², not individual |F(ρ)|².")
    print("  This is the RIGHT tool for bounding Q_off = Σ_ρ Q_ρ.")
    print()

    # Compute zero spacings
    print("  Computing 200 zeta zeros and their spacings...")
    zeros = [float(mpmath.im(mpmath.zetazero(k))) for k in range(1, 201)]
    spacings = [zeros[k+1] - zeros[k] for k in range(len(zeros)-1)]
    min_spacing = min(spacings)
    mean_spacing = np.mean(spacings)

    print(f"    Min spacing: {min_spacing:.4f}")
    print(f"    Mean spacing: {mean_spacing:.4f}")
    print()

    m_max = 3
    beta = 0.6

    for pb in [23, 47, 97, 197]:
        primes = sieve_primes(pb)
        labels = [(p, a) for p in primes for a in range(1, m_max + 1)]
        N = len(labels)

        # Worst-case primitive c
        c = np.zeros(N)
        c[0] = 1/np.sqrt(2)
        c[1] = -1/np.sqrt(2)

        # Compute ||u||² (β-dependent weights)
        norm_u_sq = sum(c[i]**2 * np.log(p) * p**(-a*(2*beta+1))
                       for i, (p, a) in enumerate(labels))

        # Compute actual Σ_ρ |F(ρ)|²
        actual_sum = 0
        for gamma in zeros:
            phases = np.array([np.exp(-1j * gamma * a * np.log(p)) for p, a in labels])
            u = np.array([c[i] * np.sqrt(np.log(p)) * p**(-a*(beta+0.5))
                         for i, (p, a) in enumerate(labels)])
            F_rho = np.sum(u * phases)
            actual_sum += abs(F_rho)**2

        # Large Sieve bound: Σ|F|² ≤ (N + ???) · ||u||²
        # The "dual" LS: Σ_r |Σ a_n e(n α_r)|² ≤ (N-1+δ^{-1}) Σ|a_n|²
        # Our phases: e^{-iγx_i} = e^{2πi · (-γ x_i / 2π)}
        # The "frequencies" are x_i / (2π) but the DUAL has zeros as nodes

        # In our setup: the DIRECT LS bounds Σ_γ |F(γ)|² where sum is over zeros
        # and F(γ) = Σ_i u_i e^{-iγx_i}
        # Here the "nodes" are γ_r (zero heights) and the "support" is x_i

        # LS: Σ_r |Σ_i u_i e^{-iγ_r x_i}|² ≤ (δ^{-1} + max x_i - min x_i) · Σ|u_i|²
        # where δ = min spacing of γ_r

        x_max = max(a * np.log(p) for p, a in labels)
        x_min_val = min(a * np.log(p) for p, a in labels)
        ls_constant = 1/min_spacing + (x_max - x_min_val)/(2*np.pi)
        ls_bound = ls_constant * norm_u_sq

        # Naive bound: Σ_ρ |F|² ≤ R · N · ||u||² where R = number of zeros
        naive_bound = len(zeros) * N * norm_u_sq

        ratio = actual_sum / norm_u_sq if norm_u_sq > 0 else 0

        print(f"  P={pb:3d}, N={N:3d}:")
        print(f"    Σ_ρ |F(ρ)|² (actual)     = {actual_sum:.6e}")
        print(f"    ||u||²                     = {norm_u_sq:.6e}")
        print(f"    Ratio Σ|F|²/||u||²         = {ratio:.2f}")
        print(f"    Large Sieve bound          = {ls_bound:.6e} (LS constant = {ls_constant:.1f})")
        print(f"    Naive bound (R·N·||u||²)   = {naive_bound:.6e}")
        print()

    print("  KEY INSIGHT: The actual Σ|F|²/||u||² ratio tells us how much")
    print("  the phases actually cancel. If ratio ~ O(1), the N factor is")
    print("  truly pessimistic and a smarter bound might work.")
    print()


def test_operator_norm():
    """
    Test: is ||M_off||_2 < 1 for the off-line contribution to the Weil matrix?

    This is the DIRECT approach: compute the spectral norm of M_off
    and check if it's less than 1.
    """
    print("=" * 70)
    print("  OPERATOR NORM OF M_off")
    print("=" * 70)
    print()

    # Build M_off for a hypothetical scenario where ALL zeros are off-line at β=0.6
    print("  Computing 50 zeta zeros...")
    zeros = [float(mpmath.im(mpmath.zetazero(k))) for k in range(1, 51)]

    beta_hyp = 0.6
    sigma = beta_hyp - 0.5
    m_max = 3

    for pb in [11, 23, 47, 97]:
        primes = sieve_primes(pb)
        labels = [(p, a) for p in primes for a in range(1, m_max + 1)]
        N = len(labels)

        M_off = np.zeros((N, N))
        for i, (pi, ai) in enumerate(labels):
            lpi = np.log(pi)
            for j, (pj, aj) in enumerate(labels):
                lpj = np.log(pj)
                x_diff = ai * lpi - aj * lpj

                # Off-line kernel: Σ_ρ cosh(σ·x)·cos(γ·x) / |ρ(1-ρ)| / π
                kernel = 0
                for gamma in zeros:
                    rho = complex(beta_hyp, gamma)
                    rho_prod = abs(rho * (1 - rho))
                    kernel += np.cosh(sigma * x_diff) * np.cos(gamma * x_diff) / rho_prod
                kernel /= np.pi

                M_off[i, j] = -np.sqrt(lpi * lpj) / (pi**(ai/2) * pj**(aj/2)) * kernel

        M_off = (M_off + M_off.T) / 2  # symmetrize
        eigs = np.linalg.eigvalsh(M_off)
        spec_norm = max(abs(eigs[0]), abs(eigs[-1]))

        print(f"  P={pb:3d}, N={N:3d}: ||M_off||_2 = {spec_norm:.6e}, "
              f"max eig = {eigs[-1]:+.6e}, min eig = {eigs[0]:+.6e}")

    print()
    print("  If ||M_off||_2 < 1 for all N: then |Q_off| < ||c||² for all c,")
    print("  which would prove APT. But does ||M_off||_2 grow with N?")
    print()


if __name__ == '__main__':
    test_cauchy_schwarz_gap()
    test_large_sieve_approach()
    test_operator_norm()
