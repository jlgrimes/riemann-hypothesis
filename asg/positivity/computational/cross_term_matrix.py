#!/usr/bin/env python3
"""
Cross-Term Matrix Analysis
===========================
Build and analyze the matrix arising from the Weil explicit formula
when expanded in a prime-power basis.

For primes p, q define:
    M(p,q) = sum_{m,n>=1} (log p log q)/(p^{m/2} q^{n/2}) K(m log p, n log q)

where K(x,y) is the Weil kernel.  If the restriction of M to the
primitive subspace (orthogonal to the constant vector) is negative
semidefinite, APT holds for this finite approximation.

This script builds M, visualizes it as a heatmap, computes eigenvalues,
and tracks how the spectrum evolves as more primes are included.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mpmath
from pathlib import Path

OUT_DIR = Path(__file__).parent / 'plots'
OUT_DIR.mkdir(exist_ok=True)


def sieve_primes(bound):
    is_prime = [True] * (bound + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(bound**0.5) + 1):
        if is_prime[i]:
            for j in range(i*i, bound + 1, i):
                is_prime[j] = False
    return [p for p in range(2, bound + 1) if is_prime[p]]


def get_zeta_zeros(n_zeros=50):
    """First n non-trivial zeta zeros (imaginary parts)."""
    mpmath.mp.dps = 30
    return np.array([float(mpmath.zetazero(k).imag) for k in range(1, n_zeros + 1)])


def weil_kernel(x, y, zeros):
    """
    Discretized Weil kernel from the explicit formula.

    K(x,y) ~ delta_{x=y} - (1/2pi) sum_gamma 2 cos(gamma(x-y)) / (1/4+gamma^2)

    The zero sum encodes how zeta zeros constrain the correlation between
    prime-power points x = m log p and y = n log q.
    """
    diag = 1.0 if abs(x - y) < 1e-10 else 0.0
    diff = x - y
    zero_sum = np.sum(2 * np.cos(zeros * diff) / (0.25 + zeros**2))
    return diag - zero_sum / (2 * np.pi)


def build_matrix(primes, m_max=5, zeros=None):
    """Build cross-term matrix M(p_i, p_j)."""
    if zeros is None:
        zeros = get_zeta_zeros(50)
    N = len(primes)
    M = np.zeros((N, N))
    for i, p in enumerate(primes):
        lp = np.log(p)
        for j, q in enumerate(primes):
            lq = np.log(q)
            s = 0.0
            for m in range(1, m_max + 1):
                for n in range(1, m_max + 1):
                    s += (lp * lq) / (p**(m/2) * q**(n/2)) * weil_kernel(m*lp, n*lq, zeros)
            M[i, j] = s
    return M


def primitive_projection(M):
    """Project M onto the subspace orthogonal to the constant vector."""
    N = M.shape[0]
    v = np.ones(N) / np.sqrt(N)
    Mv = M @ v
    vM = v @ M
    return M - np.outer(Mv, v) - np.outer(v, vM) + np.outer(v, v) * (v @ M @ v)


# -- Experiments ---------------------------------------------------------------

def experiment_heatmap():
    """Build and visualize the cross-term matrix for primes up to 50."""
    print("=" * 70)
    print("EXPERIMENT 1: Cross-term matrix heatmap")
    print("=" * 70)

    primes = sieve_primes(50)
    print(f"  {len(primes)} primes up to 50")
    print("  Computing zeta zeros ...", end=" ", flush=True)
    zeros = get_zeta_zeros(50)
    print("done")
    print("  Building matrix ...", end=" ", flush=True)
    M = build_matrix(primes, m_max=4, zeros=zeros)
    print("done")

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    im0 = axes[0].imshow(M, cmap='RdBu_r', aspect='equal', interpolation='nearest')
    axes[0].set_title('M(p,q)')
    axes[0].set_xlabel('Prime q index')
    axes[0].set_ylabel('Prime p index')
    plt.colorbar(im0, ax=axes[0])

    M_abs = np.abs(M)
    M_abs[M_abs < 1e-15] = 1e-15
    im1 = axes[1].imshow(np.log10(M_abs), cmap='viridis', aspect='equal',
                         interpolation='nearest')
    axes[1].set_title('log10 |M(p,q)|')
    axes[1].set_xlabel('Prime q index')
    axes[1].set_ylabel('Prime p index')
    plt.colorbar(im1, ax=axes[1])

    plt.tight_layout()
    plt.savefig(OUT_DIR / 'cross_term_heatmap.png', dpi=150)
    plt.close()
    print(f"  Saved {OUT_DIR / 'cross_term_heatmap.png'}")
    return M, primes, zeros


def experiment_eigenvalues(M, primes):
    """Eigenvalue analysis of the cross-term matrix."""
    print("\n" + "=" * 70)
    print("EXPERIMENT 2: Eigenvalue analysis")
    print("=" * 70)

    eigs = np.sort(np.linalg.eigvalsh(M))[::-1]
    Mp = primitive_projection(M)
    eigs_p = np.sort(np.linalg.eigvalsh(Mp))[::-1]

    print(f"\n  Full spectrum ({M.shape[0]}x{M.shape[0]}):")
    print(f"    max eigenvalue:  {eigs[0]:+.8e}")
    print(f"    min eigenvalue:  {eigs[-1]:+.8e}")
    print(f"    trace:           {np.trace(M):+.8e}")
    print(f"    #positive:       {np.sum(eigs > 1e-12)}")
    print(f"    #negative:       {np.sum(eigs < -1e-12)}")

    print(f"\n  Primitive subspace:")
    print(f"    max eigenvalue:  {eigs_p[0]:+.8e}")
    print(f"    min eigenvalue:  {eigs_p[-1]:+.8e}")
    print(f"    #positive:       {np.sum(eigs_p > 1e-12)}")
    print(f"    #negative:       {np.sum(eigs_p < -1e-12)}")

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    axes[0].stem(range(len(eigs)), eigs, markerfmt='bo', linefmt='b-', basefmt='r-')
    axes[0].axhline(0, color='r', ls='--', alpha=.5)
    axes[0].set_xlabel('Index')
    axes[0].set_ylabel('Eigenvalue')
    axes[0].set_title('Full Spectrum')
    axes[0].grid(True, alpha=.3)

    axes[1].stem(range(len(eigs_p)), eigs_p, markerfmt='go', linefmt='g-', basefmt='r-')
    axes[1].axhline(0, color='r', ls='--', alpha=.5)
    axes[1].set_xlabel('Index')
    axes[1].set_ylabel('Eigenvalue')
    axes[1].set_title('Primitive Subspace')
    axes[1].grid(True, alpha=.3)

    plt.tight_layout()
    plt.savefig(OUT_DIR / 'cross_term_eigenvalues.png', dpi=150)
    plt.close()
    print(f"  Saved {OUT_DIR / 'cross_term_eigenvalues.png'}")
    return eigs, eigs_p


def experiment_eigenvalue_growth():
    """Track largest primitive eigenvalue as more primes are included."""
    print("\n" + "=" * 70)
    print("EXPERIMENT 3: Eigenvalue growth with prime count")
    print("=" * 70)

    all_primes = sieve_primes(100)
    zeros = get_zeta_zeros(40)
    counts = [5, 10, 15, 20, 25]
    max_full = []
    min_full = []
    max_prim = []

    for c in counts:
        primes = all_primes[:c]
        M = build_matrix(primes, m_max=3, zeros=zeros)
        eigs = np.linalg.eigvalsh(M)
        Mp = primitive_projection(M)
        eigs_p = np.linalg.eigvalsh(Mp)
        max_full.append(np.max(eigs))
        min_full.append(np.min(eigs))
        max_prim.append(np.max(eigs_p))
        print(f"  {c:2d} primes: max_full={max_full[-1]:+.6e}  max_prim={max_prim[-1]:+.6e}")

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(counts, max_full, 'bo-', lw=1.5, label='max eig (full)')
    ax.plot(counts, min_full, 'rs-', lw=1.5, label='min eig (full)')
    ax.plot(counts, max_prim, 'g^-', lw=1.5, label='max eig (primitive)')
    ax.axhline(0, color='k', ls='--', alpha=.3)
    ax.set_xlabel('Number of primes')
    ax.set_ylabel('Eigenvalue')
    ax.set_title('Eigenvalue Growth with Number of Primes')
    ax.legend()
    ax.grid(True, alpha=.3)
    plt.tight_layout()
    plt.savefig(OUT_DIR / 'eigenvalue_growth.png', dpi=150)
    plt.close()
    print(f"  Saved {OUT_DIR / 'eigenvalue_growth.png'}")


def main():
    print("CROSS-TERM MATRIX ANALYSIS")
    print("=" * 70 + "\n")

    M, primes, zeros = experiment_heatmap()
    experiment_eigenvalues(M, primes)
    experiment_eigenvalue_growth()

    print("\n" + "=" * 70)
    print("ALL EXPERIMENTS COMPLETE")
    print("=" * 70)
    print("\nKey findings:")
    print("  - Diagonal dominance reflects same-prime self-interaction")
    print("  - Off-diagonal structure encodes cross-prime correlations")
    print("  - Eigenvalue sign pattern on the primitive subspace constrains APT")


if __name__ == '__main__':
    main()
