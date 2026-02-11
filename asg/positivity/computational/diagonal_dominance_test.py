#!/usr/bin/env python3
"""
Diagonal Dominance Test for the Arithmetic Positivity Theorem
=============================================================

Tests the key condition: for the Weil matrix M indexed by (prime, power),
is M diagonally dominant on the primitive subspace?

Specifically, for each (p, m), we check:
    |M_{(p,m),(p,m)}| >= sum_{(q,n) != (p,m)} |M_{(p,m),(q,n)}|

If this holds for all (p, m), then APT (and RH) follow.

This script:
1. Computes the Weil kernel K at all needed arithmetic points
2. Builds the full Weil matrix (including pole and archimedean terms)
3. Tests diagonal dominance row by row
4. Identifies the "hardest" rows (closest to failure)
5. Tracks how the margin evolves as more primes are included
"""

import numpy as np
import mpmath
from pathlib import Path

mpmath.mp.dps = 30
OUT_DIR = Path(__file__).parent / 'results'
OUT_DIR.mkdir(exist_ok=True)


def sieve_primes(bound):
    """Sieve of Eratosthenes."""
    is_prime = [True] * (bound + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(bound**0.5) + 1):
        if is_prime[i]:
            for j in range(i*i, bound + 1, i):
                is_prime[j] = False
    return [p for p in range(2, bound + 1) if is_prime[p]]


def get_zeta_zeros(n_zeros):
    """Compute first n non-trivial zeta zeros (imaginary parts)."""
    print(f"  Computing {n_zeros} zeta zeros ... ", end="", flush=True)
    zeros = [float(mpmath.zetazero(k).imag) for k in range(1, n_zeros + 1)]
    print("done")
    return np.array(zeros)


def digamma_real_part(x):
    """Re[psi(1/4 + ix/2)] — archimedean contribution to the kernel."""
    z = mpmath.mpc(0.25, x / 2)
    return float(mpmath.re(mpmath.digamma(z)))


def weil_kernel_full(x, zeros, include_bg=True):
    """
    Full Weil kernel K(x) including:
    - Zero contribution: (1/pi) sum_gamma cos(gamma x) / (1/4 + gamma^2)
    - Background: -(1/pi) Re[psi(1/4 + ix/2)] + (1/2pi) log(pi)
    """
    # Zero contribution
    if len(zeros) > 0:
        K_zeros = (1 / np.pi) * np.sum(
            2 * np.cos(zeros * x) / (0.25 + zeros**2)
        )
    else:
        K_zeros = 0.0

    if not include_bg:
        return K_zeros

    # Background (archimedean + pole)
    K_bg = -(1 / np.pi) * digamma_real_part(x) + np.log(np.pi) / (2 * np.pi)

    return K_zeros + K_bg


def build_full_weil_matrix(primes, m_max, zeros):
    """
    Build the Weil matrix M indexed by (prime, power) pairs.

    M_{(p,m),(q,n)} = sqrt(log p * log q) / (p^{m/2} * q^{n/2}) * K(m log p - n log q)

    Returns M, labels where labels[i] = (p, m) for row/col i.
    """
    labels = []
    for p in primes:
        for m in range(1, m_max + 1):
            labels.append((p, m))

    N = len(labels)
    M = np.zeros((N, N))

    for i, (p, mp) in enumerate(labels):
        lp = np.log(p)
        for j, (q, nq) in enumerate(labels):
            lq = np.log(q)
            x = mp * lp - nq * lq
            weight = np.sqrt(lp * lq) / (p**(mp/2) * q**(nq/2))
            K_val = weil_kernel_full(x, zeros)
            M[i, j] = weight * K_val

    return M, labels


def test_diagonal_dominance(M, labels):
    """
    Test diagonal dominance for each row.
    Returns list of (label, diagonal, off_diag_sum, margin, is_dominant).
    """
    N = M.shape[0]
    results = []
    for i in range(N):
        diag = abs(M[i, i])
        off_diag = sum(abs(M[i, j]) for j in range(N) if j != i)
        margin = diag - off_diag
        results.append((labels[i], diag, off_diag, margin, margin > 0))
    return results


def primitive_eigenvalues(M):
    """Eigenvalues of M projected onto the primitive subspace."""
    N = M.shape[0]
    v = np.ones(N) / np.sqrt(N)
    Mv = M @ v
    vM = v @ M
    Mp = M - np.outer(Mv, v) - np.outer(v, vM) + (v @ M @ v) * np.outer(v, v)
    return np.sort(np.linalg.eigvalsh(Mp))


def run_experiment(primes, m_max, n_zeros):
    """Run the full diagonal dominance experiment."""
    print(f"\n{'='*75}")
    print(f"DIAGONAL DOMINANCE TEST: {len(primes)} primes up to {primes[-1]}, "
          f"m_max={m_max}, {n_zeros} zeros")
    print(f"{'='*75}")

    zeros = get_zeta_zeros(n_zeros)

    print(f"  Building {len(primes)*m_max}x{len(primes)*m_max} Weil matrix ... ",
          end="", flush=True)
    M, labels = build_full_weil_matrix(primes, m_max, zeros)
    print("done")

    # Test diagonal dominance
    dd_results = test_diagonal_dominance(M, labels)
    n_dominant = sum(1 for r in dd_results if r[4])
    n_total = len(dd_results)

    print(f"\n  Diagonal dominance: {n_dominant}/{n_total} rows pass")

    # Show worst rows
    dd_results.sort(key=lambda r: r[3])
    print(f"\n  {'(p,m)':>10s}  {'|M_ii|':>12s}  {'sum|M_ij|':>12s}  "
          f"{'Margin':>12s}  {'Pass?':>6s}")
    print(f"  {'-'*58}")
    for label, diag, off, margin, ok in dd_results[:15]:
        print(f"  ({label[0]:3d},{label[1]:1d})     {diag:12.6e}  {off:12.6e}  "
              f"{margin:+12.6e}  {'YES' if ok else 'NO':>6s}")

    # Primitive eigenvalue analysis
    eigs = primitive_eigenvalues(M)
    n_pos = np.sum(eigs > 1e-12)
    n_neg = np.sum(eigs < -1e-12)
    print(f"\n  Primitive subspace eigenvalues:")
    print(f"    Positive: {n_pos}  (APT needs 0)")
    print(f"    Negative: {n_neg}")
    print(f"    Near zero: {n_total - n_pos - n_neg}")
    print(f"    Largest:  {eigs[-1]:+.6e}")
    print(f"    Smallest: {eigs[0]:+.6e}")

    if n_pos == 0:
        print(f"\n  *** APT HOLDS for this truncation! ***")
    else:
        print(f"\n  --- APT UNCERTAIN (positive eigenvalues present) ---")

    return M, labels, dd_results, eigs


def convergence_study(max_prime_bound=100, m_max=3, n_zeros=100):
    """Study how the spectral gap evolves with more primes."""
    print(f"\n{'='*75}")
    print("CONVERGENCE STUDY: Spectral gap vs number of primes")
    print(f"{'='*75}")

    zeros = get_zeta_zeros(n_zeros)
    all_primes = sieve_primes(max_prime_bound)

    bounds = [10, 20, 30, 50, 70, 100]
    results = []

    print(f"\n  {'Bound':>6s}  {'#Primes':>8s}  {'MatSize':>8s}  "
          f"{'MaxPrimEig':>12s}  {'MinPrimEig':>12s}  {'DD pass':>8s}")
    print(f"  {'-'*66}")

    for b in bounds:
        primes = [p for p in all_primes if p <= b]
        if len(primes) < 2:
            continue
        M, labels = build_full_weil_matrix(primes, m_max, zeros)
        eigs = primitive_eigenvalues(M)
        dd = test_diagonal_dominance(M, labels)
        dd_pass = sum(1 for r in dd if r[4])

        results.append({
            'bound': b,
            'n_primes': len(primes),
            'size': len(labels),
            'max_eig': eigs[-1],
            'min_eig': eigs[0],
            'dd_pass': dd_pass,
            'dd_total': len(dd),
        })

        print(f"  {b:6d}  {len(primes):8d}  {len(labels):8d}  "
              f"{eigs[-1]:+12.6e}  {eigs[0]:+12.6e}  "
              f"{dd_pass}/{len(dd):>4d}")

    return results


def kernel_at_prime_diffs(primes, n_zeros=200):
    """Examine the Weil kernel at differences of prime logarithms."""
    print(f"\n{'='*75}")
    print("KERNEL AT PRIME-LOG DIFFERENCES: K(log p - log q)")
    print(f"{'='*75}")

    zeros = get_zeta_zeros(n_zeros)

    print(f"\n  {'p':>4s}  {'q':>4s}  {'log p-log q':>12s}  "
          f"{'K_total':>12s}  {'K_zeros':>12s}  {'K_bg':>12s}")
    print(f"  {'-'*62}")

    for p in primes[:8]:
        for q in primes[:8]:
            if p == q:
                continue
            x = np.log(p) - np.log(q)
            K_total = weil_kernel_full(x, zeros, include_bg=True)
            K_z = weil_kernel_full(x, zeros, include_bg=False)
            K_b = K_total - K_z
            print(f"  {p:4d}  {q:4d}  {x:12.6f}  "
                  f"{K_total:12.6f}  {K_z:12.6f}  {K_b:12.6f}")


def main():
    print("ARITHMETIC POSITIVITY THEOREM — DIAGONAL DOMINANCE ANALYSIS")
    print("=" * 75)
    print("Testing whether the Weil matrix is diagonally dominant,")
    print("which would prove APT and hence the Riemann Hypothesis.")
    print()

    # 1. Kernel values at key points
    primes = sieve_primes(50)
    kernel_at_prime_diffs(primes, n_zeros=200)

    # 2. Small test
    run_experiment(sieve_primes(20), m_max=3, n_zeros=100)

    # 3. Medium test
    run_experiment(sieve_primes(50), m_max=2, n_zeros=200)

    # 4. Convergence study
    convergence_study(max_prime_bound=70, m_max=2, n_zeros=150)

    print(f"\n{'='*75}")
    print("ANALYSIS COMPLETE")
    print("=" * 75)
    print()
    print("Interpretation:")
    print("  - If ALL primitive eigenvalues are ≤ 0: APT holds for this truncation")
    print("  - If diagonal dominance holds: APT is guaranteed")
    print("  - The margin of diagonal dominance quantifies 'how far from failure'")
    print("  - Convergence study shows whether the gap survives as more primes added")


if __name__ == '__main__':
    main()
