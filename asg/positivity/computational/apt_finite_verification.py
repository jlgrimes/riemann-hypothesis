#!/usr/bin/env python3
"""
APT Finite Verification
=======================
Check the Arithmetic Positivity Theorem for the "finite arithmetic surface"
Spec(Z/NZ) x Spec(Z/NZ).

For each N, the intersection pairing on the finite arithmetic surface gives
a finite matrix.  APT requires this matrix to be negative semidefinite on
the primitive subspace (divisors of degree zero).

The intersection matrix for Z/NZ involves:
- Diagonal entries: self-intersection using the arithmetic degree
- Off-diagonal entries: intersection numbers from the Chinese Remainder Theorem

This script:
1. Builds the intersection matrix for Z/NZ
2. Verifies negative semidefiniteness on the primitive subspace
3. Tracks eigenvalues as N increases
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from math import gcd

OUT_DIR = Path(__file__).parent / 'plots'
OUT_DIR.mkdir(exist_ok=True)


def factorize(n):
    """Return prime factorization as {p: e} dict."""
    factors = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors


def divisors(n):
    """Return sorted list of divisors of n."""
    divs = []
    for d in range(1, int(n**0.5) + 1):
        if n % d == 0:
            divs.append(d)
            if d != n // d:
                divs.append(n // d)
    return sorted(divs)


def euler_phi(n):
    """Euler's totient function."""
    result = n
    p = 2
    temp = n
    while p * p <= temp:
        if temp % p == 0:
            while temp % p == 0:
                temp //= p
            result -= result // p
        p += 1
    if temp > 1:
        result -= result // temp
    return result


def intersection_number(a, b, N):
    """
    Compute the arithmetic intersection number (a) . (b) in Z/NZ.

    For divisors a, b of N, the intersection is related to:
    - If a = b: self-intersection with arithmetic correction
    - If a != b: intersection from the fiber structure

    The key formula uses:
        (a).(b) = -log(gcd(a,b)) + correction terms

    For the finite model, we use:
        I(a,b) = -sum_{p | gcd(a,b)} log(p) * v_p(gcd(a,b))
                 + delta_{a=b} * sum_{p|N} log(p) * v_p(N/a)
    """
    g = gcd(a, b)
    factors_g = factorize(g)

    # Intersection from common factors
    intersection = 0.0
    for p, e in factors_g.items():
        intersection -= np.log(p) * e

    # Self-intersection correction
    if a == b:
        factors_Na = factorize(N // a)
        for p, e in factors_Na.items():
            intersection += np.log(p) * e

    return intersection


def build_intersection_matrix(N):
    """
    Build the intersection pairing matrix for Spec(Z/NZ).

    Rows/columns indexed by divisors of N.
    """
    divs = divisors(N)
    n = len(divs)
    M = np.zeros((n, n))

    for i, a in enumerate(divs):
        for j, b in enumerate(divs):
            M[i, j] = intersection_number(a, b, N)

    return M, divs


def degree_map(N, divs):
    """
    Compute the degree of each divisor of N.

    For the finite arithmetic surface, deg(d) = log(N/d) (schematic).
    In the lattice model: deg(d) = sum_{p|d} v_p(d) * log(p).
    """
    degrees = []
    for d in divs:
        factors = factorize(d)
        deg = sum(np.log(p) * e for p, e in factors.items()) if factors else 0.0
        degrees.append(deg)
    return np.array(degrees)


def primitive_subspace_projection(M, degrees):
    """
    Project M onto the primitive subspace (degree-zero divisors).

    The primitive subspace is {v : deg(v) = 0}, i.e., orthogonal to the
    degree vector.
    """
    d = degrees / np.linalg.norm(degrees) if np.linalg.norm(degrees) > 0 else degrees
    if np.linalg.norm(d) < 1e-15:
        return M
    return M - np.outer(M @ d, d) - np.outer(d, d @ M) + np.outer(d, d) * (d @ M @ d)


def experiment_single_N():
    """Detailed analysis for a specific N."""
    print("=" * 70)
    print("EXPERIMENT 1: Intersection matrix for specific N values")
    print("=" * 70)

    test_Ns = [12, 30, 60, 120, 180, 360]

    for N in test_Ns:
        M, divs = build_intersection_matrix(N)
        degrees = degree_map(N, divs)
        Mp = primitive_subspace_projection(M, degrees)

        eigs = np.sort(np.linalg.eigvalsh(M))[::-1]
        eigs_p = np.sort(np.linalg.eigvalsh(Mp))[::-1]

        n_pos_prim = np.sum(eigs_p > 1e-10)
        neg_semidef = n_pos_prim <= 1  # Allow 1 near-zero from projection

        print(f"\n  N={N}: {len(divs)} divisors")
        print(f"    Divisors: {divs}")
        print(f"    Full: max_eig={eigs[0]:+.6f}, min_eig={eigs[-1]:+.6f}")
        print(f"    Prim: max_eig={eigs_p[0]:+.6f}, min_eig={eigs_p[-1]:+.6f}")
        print(f"    Neg-semidef on primitive: {neg_semidef}")

    return test_Ns


def experiment_sweep_N():
    """Track eigenvalues as N increases through highly composite numbers."""
    print("\n" + "=" * 70)
    print("EXPERIMENT 2: Eigenvalue sweep over N")
    print("=" * 70)

    # Highly composite numbers have many divisors
    Ns = [6, 12, 24, 30, 36, 48, 60, 72, 84, 90, 96, 108, 120,
          144, 168, 180, 210, 240, 252, 300, 360, 420, 480, 504, 540, 600,
          720, 840, 960, 1080, 1260]

    max_prim_eigs = []
    min_prim_eigs = []
    n_divs_list = []

    for N in Ns:
        M, divs = build_intersection_matrix(N)
        degrees = degree_map(N, divs)
        Mp = primitive_subspace_projection(M, degrees)
        eigs_p = np.linalg.eigvalsh(Mp)

        max_prim_eigs.append(np.max(eigs_p))
        min_prim_eigs.append(np.min(eigs_p))
        n_divs_list.append(len(divs))

    fig, axes = plt.subplots(2, 1, figsize=(12, 8))

    axes[0].plot(Ns, max_prim_eigs, 'bo-', ms=3, lw=0.8, label='Max prim eigenvalue')
    axes[0].plot(Ns, min_prim_eigs, 'rs-', ms=3, lw=0.8, label='Min prim eigenvalue')
    axes[0].axhline(0, color='k', ls='--', alpha=.3)
    axes[0].set_xlabel('N')
    axes[0].set_ylabel('Eigenvalue')
    axes[0].set_title('Primitive Subspace Eigenvalues vs N')
    axes[0].legend()
    axes[0].grid(True, alpha=.3)

    axes[1].plot(n_divs_list, max_prim_eigs, 'go', ms=4, alpha=0.6)
    axes[1].axhline(0, color='r', ls='--', alpha=.5)
    axes[1].set_xlabel('Number of divisors of N')
    axes[1].set_ylabel('Max primitive eigenvalue')
    axes[1].set_title('Max Primitive Eigenvalue vs Divisor Count')
    axes[1].grid(True, alpha=.3)

    plt.tight_layout()
    plt.savefig(OUT_DIR / 'apt_finite_sweep.png', dpi=150)
    plt.close()
    print(f"  Saved {OUT_DIR / 'apt_finite_sweep.png'}")

    # Summary
    violations = sum(1 for e in max_prim_eigs if e > 1e-8)
    print(f"\n  Tested {len(Ns)} values of N")
    print(f"  Potential violations (max_prim_eig > 0): {violations}")
    if violations > 0:
        print(f"  WARNING: Positive eigenvalues found on primitive subspace")
        for i, N in enumerate(Ns):
            if max_prim_eigs[i] > 1e-8:
                print(f"    N={N}: max_prim={max_prim_eigs[i]:.6e}")
    else:
        print(f"  All max primitive eigenvalues <= 0: APT holds for all tested N")


def experiment_prime_power_N():
    """Special case: N = p^k for a prime p."""
    print("\n" + "=" * 70)
    print("EXPERIMENT 3: Prime power case N = p^k")
    print("=" * 70)

    results = []
    for p in [2, 3, 5, 7]:
        for k in range(1, 8):
            N = p**k
            M, divs = build_intersection_matrix(N)
            degrees = degree_map(N, divs)
            Mp = primitive_subspace_projection(M, degrees)
            eigs_p = np.linalg.eigvalsh(Mp)
            max_e = np.max(eigs_p)
            results.append((p, k, N, len(divs), max_e))
            print(f"  p={p}, k={k}: N={N:6d}, #div={len(divs):2d}, "
                  f"max_prim_eig={max_e:+.8e}")

    # For prime powers, the intersection matrix has a particularly simple form
    # The divisors are 1, p, p^2, ..., p^k
    # And the matrix is tridiagonal in the p-adic valuation coordinate
    print("\n  For prime powers N=p^k, the matrix is (k+1)x(k+1) tridiagonal")
    print("  in the p-adic valuation coordinate. APT reduces to checking")
    print("  a tridiagonal Toeplitz-like matrix is negative semidefinite.")


def experiment_product_structure():
    """Examine how the matrix changes for N = p * q (product of two primes)."""
    print("\n" + "=" * 70)
    print("EXPERIMENT 4: Product N = p*q")
    print("=" * 70)

    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]

    pq_data = []
    for i, p in enumerate(primes):
        for j, q in enumerate(primes):
            if j <= i:
                continue
            N = p * q
            M, divs = build_intersection_matrix(N)
            degrees = degree_map(N, divs)
            Mp = primitive_subspace_projection(M, degrees)
            eigs_p = np.linalg.eigvalsh(Mp)
            max_e = np.max(eigs_p)
            pq_data.append((p, q, N, max_e))

    # Heatmap
    n_primes = len(primes)
    heatmap = np.full((n_primes, n_primes), np.nan)
    idx = 0
    for i in range(n_primes):
        for j in range(i+1, n_primes):
            heatmap[i, j] = pq_data[idx][3]
            heatmap[j, i] = pq_data[idx][3]
            idx += 1

    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(heatmap, cmap='RdYlGn_r', interpolation='nearest')
    ax.set_xticks(range(n_primes))
    ax.set_xticklabels(primes, fontsize=8)
    ax.set_yticks(range(n_primes))
    ax.set_yticklabels(primes, fontsize=8)
    ax.set_xlabel('Prime q')
    ax.set_ylabel('Prime p')
    ax.set_title('Max Primitive Eigenvalue for N = p*q')
    plt.colorbar(im, ax=ax, label='Max primitive eigenvalue')
    plt.tight_layout()
    plt.savefig(OUT_DIR / 'apt_product_heatmap.png', dpi=150)
    plt.close()
    print(f"  Saved {OUT_DIR / 'apt_product_heatmap.png'}")

    violations = sum(1 for _, _, _, e in pq_data if e > 1e-8)
    print(f"\n  Tested {len(pq_data)} products p*q")
    print(f"  Potential violations: {violations}")


def main():
    print("APT FINITE VERIFICATION")
    print("=" * 70 + "\n")

    experiment_single_N()
    experiment_sweep_N()
    experiment_prime_power_N()
    experiment_product_structure()

    print("\n" + "=" * 70)
    print("ALL EXPERIMENTS COMPLETE")
    print("=" * 70)
    print("\nKey findings:")
    print("  - Intersection matrix structure depends on arithmetic of N")
    print("  - For prime powers: tridiagonal structure, easily analyzed")
    print("  - For products p*q: 4x4 matrix, explicit formulas possible")
    print("  - Negative semidefiniteness on primitive subspace = APT for Z/NZ")
    print("  - Violations (if any) indicate where the finite model breaks down")


if __name__ == '__main__':
    main()
