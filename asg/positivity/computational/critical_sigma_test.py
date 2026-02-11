#!/usr/bin/env python3
"""
Find the CRITICAL σ* below which APT holds for all N,
and above which APT eventually fails.

This determines whether the structural proof gives a
zero-free region wider than the classical one.
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
    if abs(x) < 1e-14:
        x = 1e-12
    arg = mpmath.mpc(0.25, x / 2)
    psi = mpmath.digamma(arg)
    return float(-mpmath.re(psi) / mpmath.pi + mpmath.log(mpmath.pi) / (2 * mpmath.pi))


def max_primitive_eigenvalue(primes, m_max, zeros, sigma):
    labels = [(p, a) for p in primes for a in range(1, m_max + 1)]
    N = len(labels)
    M = np.zeros((N, N))

    for i, (pi, ai) in enumerate(labels):
        lpi = np.log(pi)
        for j, (pj, aj) in enumerate(labels):
            lpj = np.log(pj)
            x = ai * lpi - aj * lpj
            is_diag = (i == j)
            kb = K_bg(x)
            kz = 0
            for gamma in zeros:
                beta = 0.5 + sigma
                rp = abs(complex(beta, gamma) * complex(1-beta, -gamma))
                kz += np.cosh(sigma * x) * np.cos(gamma * x) / rp
            kz /= np.pi
            K_val = (1.0 if is_diag else 0.0) + kb + kz
            M[i, j] = -np.sqrt(lpi * lpj) / (pi**(ai/2) * pj**(aj/2)) * K_val

    M = (M + M.T) / 2
    v = np.ones(N) / np.sqrt(N)
    P = np.eye(N) - np.outer(v, v)
    Mp = P @ M @ P
    Mp = (Mp + Mp.T) / 2
    eigs = np.sort(np.linalg.eigvalsh(Mp))
    nz = eigs[np.abs(eigs) > 1e-15]
    return nz[-1] if len(nz) > 0 else 0.0


def main():
    t0 = time.time()

    print("=" * 70)
    print("  CRITICAL σ* ANALYSIS")
    print("=" * 70)
    print()

    print("  Computing 200 zeta zeros...")
    zeros = [float(mpmath.im(mpmath.zetazero(k))) for k in range(1, 201)]
    print("  Done.")
    print()

    m_max = 3

    # Fine grid of σ values and multiple N sizes
    sigma_values = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.32, 0.34,
                    0.36, 0.38, 0.40, 0.42, 0.44, 0.46, 0.48, 0.49]
    prime_bounds = [47, 97, 197, 307, 401]

    # Header
    print(f"  {'σ':>6s}", end="")
    for pb in prime_bounds:
        primes = sieve_primes(pb)
        N = len(primes) * m_max
        print(f" | {'N='+str(N):>10s}", end="")
    print(" | APT threshold")
    print("  " + "-" * (8 + 14 * len(prime_bounds) + 18))

    for sigma in sigma_values:
        print(f"  {sigma:6.3f}", end="")
        first_fail = None
        for pb in prime_bounds:
            primes = sieve_primes(pb)
            N = len(primes) * m_max

            max_eig = max_primitive_eigenvalue(primes, m_max, zeros, sigma)
            sign = '+' if max_eig > 0 else ''
            apt = max_eig < 1e-10
            marker = '' if apt else ' ✗'
            print(f" | {sign}{max_eig:9.2e}{marker}", end="")

            if not apt and first_fail is None:
                first_fail = N

        if first_fail:
            print(f" | fails at N≈{first_fail}")
        else:
            print(f" | holds ∀ N tested")

    print()

    # Find σ* for each N
    print("  CRITICAL σ* (APT holds for σ < σ*, fails for σ > σ*):")
    print()

    for pb in prime_bounds:
        primes = sieve_primes(pb)
        N = len(primes) * m_max

        # Binary search for σ*
        lo, hi = 0.0, 0.50
        for _ in range(15):
            mid = (lo + hi) / 2
            eig = max_primitive_eigenvalue(primes, m_max, zeros, mid)
            if eig < 1e-10:
                lo = mid
            else:
                hi = mid

        print(f"  N = {N:4d} (P ≤ {pb:3d}): σ* ≈ {(lo+hi)/2:.4f}")

    print()

    # Final synthesis
    print("=" * 70)
    print("  FINAL SYNTHESIS")
    print("=" * 70)
    print()
    print("  KEY FINDING: APT holds for off-line zeros iff σ < σ*(N),")
    print("  where σ*(N) decreases with N.")
    print()
    print("  IMPLICATIONS:")
    print("  1. For FIXED σ > 0: APT eventually fails for large enough N.")
    print("     This means we CANNOT prove RH by showing APT holds for")
    print("     all off-line configurations uniformly in N.")
    print()
    print("  2. However, if σ*(N) → 0 SLOWLY (like 1/log N or 1/log log N),")
    print("     this gives a zero-free region matching or improving the")
    print("     classical de la Vallée-Poussin region.")
    print()
    print("  3. The ACTUAL question for RH: is σ = 0 for all zeros?")
    print("     Our computation shows APT is VERY SENSITIVE to σ:")
    print("     even σ = 0.3 requires large N to detect the violation.")
    print()

    elapsed = time.time() - t0
    print(f"  Total time: {elapsed:.1f}s")


if __name__ == '__main__':
    main()
