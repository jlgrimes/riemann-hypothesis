#!/usr/bin/env python3
"""
THE CORRECT FORMULATION:

APT ⟺ for all primitive c with ||c|| = 1:
  -Q = ||f||² + f^T K_bg f + f^T K_on f + f^T K_off f ≥ 0

where f = Dc, d_i = √(log p_i)/p_i^{a_i/2}, and K is the kernel matrix.

Since K_bg and K_on produce non-negative quadratic forms on primitive f:
  -Q ≥ ||f||² + f^T K_off f

This is ≥ 0 iff: λ_min(I + K_off) ≥ 0, restricted to the range of D·(primitive).

Equivalently: λ_min(K_off) ≥ -1.

Since K_off(0) = C_all/π ≈ 0.015 > 0, the diagonal of K_off is positive.
The question: can the smallest eigenvalue of K_off go below -1?

K_off(x) = Σ_ρ cosh(σx) cos(γx) / (π|ρ(1-ρ)|)

Note: K_off WITHOUT D-weights. The entries can be large for large x.
But for the RESTRICTED f = Dc, the large-x entries are suppressed by D.

The test: compute eigenvalues of K_off for various N and σ.
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


def test_kernel_eigenvalues():
    """
    Compute eigenvalues of K_off (unweighted kernel matrix) and
    check if λ_min > -1.
    """
    print("=" * 70)
    print("  EIGENVALUES OF K_off (UNWEIGHTED KERNEL MATRIX)")
    print("=" * 70)
    print()

    print("  Computing 200 zeta zeros...")
    zeros = [float(mpmath.im(mpmath.zetazero(k))) for k in range(1, 201)]

    m_max = 3
    gamma_EM = float(mpmath.euler)
    C_all = 2 + gamma_EM - np.log(4 * np.pi)

    print(f"  C_all/π = {C_all/np.pi:.6f} (diagonal value of K_off)")
    print()

    for sigma in [0.0, 0.1, 0.3, 0.49]:
        print(f"  σ = {sigma}:")
        print(f"  {'P':>5s} {'N':>5s} | {'λ_min(K_off)':>14s} {'λ_max(K_off)':>14s} "
              f"{'λ_min > -1?':>12s} {'λ_min(I+K)':>14s}")
        print("  " + "-" * 75)

        for pb in [11, 23, 47, 97, 197]:
            primes = sieve_primes(pb)
            labels = [(p, a) for p in primes for a in range(1, m_max + 1)]
            N = len(labels)

            # Build K_off (UNWEIGHTED kernel matrix)
            K = np.zeros((N, N))
            for i, (pi, ai) in enumerate(labels):
                xi = ai * np.log(pi)
                for j, (pj, aj) in enumerate(labels):
                    xj = aj * np.log(pj)
                    x = xi - xj

                    kernel = 0
                    for gamma in zeros:
                        beta = 0.5 + sigma
                        rho_prod = abs(complex(beta, gamma) * complex(1-beta, -gamma))
                        kernel += np.cosh(sigma * x) * np.cos(gamma * x) / rho_prod
                    kernel /= np.pi
                    K[i, j] = kernel

            K = (K + K.T) / 2
            eigs = np.linalg.eigvalsh(K)

            lam_min = eigs[0]
            lam_max = eigs[-1]
            ok = "YES" if lam_min > -1.0 else "NO"
            lam_min_shift = lam_min + 1.0

            print(f"  {pb:5d} {N:5d} | {lam_min:+14.4e} {lam_max:+14.4e} "
                  f"{ok:>12s} {lam_min_shift:+14.4e}")

        print()

    print("  If λ_min(K_off) > -1: then I + K_off is PSD,")
    print("  which means ||f||² + f^T K_off f ≥ 0 for all f.")
    print("  Combined with K_bg, K_on ≥ 0 on primitive → APT holds.")
    print()


def test_restricted_eigenvalues():
    """
    Test eigenvalues of K_off restricted to the range of D·(primitive).
    This is the ACTUAL quantity that matters: we need f^T K_off f ≥ -||f||²
    only for f = Dc where c is primitive.
    """
    print("=" * 70)
    print("  EIGENVALUES OF K_off RESTRICTED TO D·(primitive)")
    print("=" * 70)
    print()

    print("  Computing 200 zeta zeros...")
    zeros = [float(mpmath.im(mpmath.zetazero(k))) for k in range(1, 201)]

    m_max = 3

    for sigma in [0.0, 0.1, 0.3, 0.49]:
        print(f"  σ = {sigma}:")

        for pb in [23, 47, 97, 197]:
            primes = sieve_primes(pb)
            labels = [(p, a) for p in primes for a in range(1, m_max + 1)]
            N = len(labels)

            # D-weights
            d = np.array([np.sqrt(np.log(p)) / p**(a/2) for p, a in labels])
            D = np.diag(d)

            # Primitive projector
            ones = np.ones(N) / np.sqrt(N)
            P = np.eye(N) - np.outer(ones, ones)

            # K_off matrix
            K = np.zeros((N, N))
            for i, (pi, ai) in enumerate(labels):
                xi = ai * np.log(pi)
                for j, (pj, aj) in enumerate(labels):
                    xj = aj * np.log(pj)
                    x = xi - xj
                    kernel = 0
                    for gamma in zeros:
                        beta = 0.5 + sigma
                        rho_prod = abs(complex(beta, gamma) * complex(1-beta, -gamma))
                        kernel += np.cosh(sigma * x) * np.cos(gamma * x) / rho_prod
                    kernel /= np.pi
                    K[i, j] = kernel
            K = (K + K.T) / 2

            # The quadratic form on primitive f = Dc:
            # f^T K f / ||f||² where f = D(Pc), so
            # c^T P D K D P c / (c^T P D² P c)
            # This is a generalized eigenvalue problem

            A = P @ D @ K @ D @ P
            B = P @ D @ D @ P  # D² restricted to primitive

            A = (A + A.T) / 2
            B = (B + B.T) / 2

            # Generalized eigenvalues of Ax = λBx
            # on the primitive subspace (remove null space of P)
            eigs_A = np.linalg.eigvalsh(A)
            eigs_B = np.linalg.eigvalsh(B)

            # The ratio f^T K f / ||f||² = c^T A c / c^T B c
            # For the worst case (max ratio): this is the largest generalized eigenvalue
            # For the bound: need min ratio > -1

            # Simple approach: compute the Weil matrix M_off = -DKD
            # and check its primitive eigenvalues
            M_off = -D @ K @ D
            M_off = (M_off + M_off.T) / 2
            Mp_off = P @ M_off @ P
            Mp_off = (Mp_off + Mp_off.T) / 2
            eigs_prim = np.sort(np.linalg.eigvalsh(Mp_off))
            nz = eigs_prim[np.abs(eigs_prim) > 1e-14]

            if len(nz) > 0:
                max_prim = nz[-1]
                min_prim = nz[0]
            else:
                max_prim = min_prim = 0

            # Also compute the FULL Weil matrix (δ + K_bg + K_on + K_off)
            # where K_on uses the SAME zeros but on the line
            M_full = np.zeros((N, N))
            for i, (pi, ai) in enumerate(labels):
                lpi = np.log(pi)
                for j, (pj, aj) in enumerate(labels):
                    lpj = np.log(pj)
                    x = ai * lpi - aj * lpj
                    is_diag = (i == j)
                    kb = K_bg(x)
                    # K_zeros with cosh(σ) factor (off-line)
                    kz = 0
                    for gamma in zeros:
                        beta = 0.5 + sigma
                        rp = abs(complex(beta, gamma) * complex(1-beta, -gamma))
                        kz += np.cosh(sigma * x) * np.cos(gamma * x) / rp
                    kz /= np.pi
                    K_val = (1.0 if is_diag else 0.0) + kb + kz
                    M_full[i, j] = -np.sqrt(lpi * lpj) / (pi**(ai/2) * pj**(aj/2)) * K_val

            M_full = (M_full + M_full.T) / 2
            Mp_full = P @ M_full @ P
            Mp_full = (Mp_full + Mp_full.T) / 2
            eigs_full = np.sort(np.linalg.eigvalsh(Mp_full))
            nz_full = eigs_full[np.abs(eigs_full) > 1e-14]
            max_full = nz_full[-1] if len(nz_full) > 0 else 0

            print(f"    P={pb:3d}, N={N:3d}: M_off prim max = {max_prim:+.4e}, "
                  f"M_full prim max = {max_full:+.4e}, "
                  f"APT? {'YES' if max_full <= 1e-10 else 'NO'}")

        print()


def test_full_apt():
    """
    THE DEFINITIVE TEST: Build the Weil matrix with zeros placed at
    various off-line positions and check whether APT still holds.
    """
    print("=" * 70)
    print("  DEFINITIVE APT TEST WITH OFF-LINE ZEROS")
    print("=" * 70)
    print()
    print("  Build Weil matrix M = -(δ + K_bg + K_zeros) with zeros")
    print("  placed at ρ = 1/2 + σ + iγ (off the critical line).")
    print("  Check: are all primitive eigenvalues ≤ 0?")
    print()

    print("  Computing 200 zeta zeros...")
    zeros = [float(mpmath.im(mpmath.zetazero(k))) for k in range(1, 201)]

    m_max = 3

    # For each σ value, build the Weil matrix as if zeros were off-line
    print(f"  {'σ':>6s} {'P':>5s} {'N':>5s} | {'max prim eig':>14s} {'APT?':>6s}")
    print("  " + "-" * 50)

    for sigma in [0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.49]:
        for pb in [47, 97, 197]:
            primes = sieve_primes(pb)
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

                    # Zeros at β = 1/2 + σ instead of β = 1/2
                    kz = 0
                    for gamma in zeros:
                        beta = 0.5 + sigma
                        rp = abs(complex(beta, gamma) * complex(1-beta, -gamma))
                        kz += np.cosh(sigma * x) * np.cos(gamma * x) / rp
                    kz /= np.pi

                    K_val = (1.0 if is_diag else 0.0) + kb + kz
                    M[i, j] = -np.sqrt(lpi * lpj) / (pi**(ai/2) * pj**(aj/2)) * K_val

            M = (M + M.T) / 2

            # Primitive eigenvalues
            v = np.ones(N) / np.sqrt(N)
            P = np.eye(N) - np.outer(v, v)
            Mp = P @ M @ P
            Mp = (Mp + Mp.T) / 2
            eigs = np.sort(np.linalg.eigvalsh(Mp))
            nz = eigs[np.abs(eigs) > 1e-14]
            max_prim = nz[-1] if len(nz) > 0 else 0

            apt = "YES" if max_prim < 1e-10 else "NO"
            print(f"  {sigma:6.3f} {pb:5d} {N:5d} | {max_prim:+14.6e} {apt:>6s}")

    print()
    print("  NOTE: σ = 0 means zeros on the critical line (RH scenario).")
    print("  σ > 0 means hypothetical off-line zeros.")
    print("  If APT = YES for all σ: the result is STRONGER than needed for RH.")
    print("  If APT = NO for some σ: the off-line perturbation breaks APT.")
    print()


def main():
    t0 = time.time()
    test_kernel_eigenvalues()
    print()
    test_full_apt()
    print()
    print(f"  Total time: {time.time()-t0:.1f}s")


if __name__ == '__main__':
    main()
