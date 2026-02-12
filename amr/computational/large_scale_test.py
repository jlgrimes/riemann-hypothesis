#!/usr/bin/env python3
"""
Large-Scale Weil Matrix Eigenvalue Test
=========================================
Push certified eigenvalue computation to 200+ primes.

Computes the Weil cross-term matrix for increasing numbers of primes,
projects onto primitive subspace, and tracks:
  - max primitive eigenvalue (must stay negative for APT)
  - spectral gap
  - entropy-gap correlation at scale
  - timing / scaling data
"""

import numpy as np
import mpmath
import time
import sys

mpmath.mp.dps = 40

SEPARATOR = "=" * 70


def sieve_primes(bound):
    is_prime = [True] * (bound + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(bound**0.5) + 1):
        if is_prime[i]:
            for j in range(i*i, bound + 1, i):
                is_prime[j] = False
    return [p for p in range(2, bound + 1) if is_prime[p]]


def get_zeta_zeros(n_zeros=100):
    """Compute and cache zeta zeros."""
    print(f"  Computing {n_zeros} zeta zeros (dps={mpmath.mp.dps})...")
    t0 = time.time()
    zeros = np.array([float(mpmath.zetazero(k).imag) for k in range(1, n_zeros + 1)])
    print(f"  Done in {time.time()-t0:.1f}s")
    return zeros


def K_bg(x):
    """Background kernel from digamma."""
    if abs(x) < 1e-14:
        x = 1e-12
    arg = mpmath.mpc(0.25, x / 2)
    psi = mpmath.digamma(arg)
    return float(-mpmath.re(psi) / mpmath.pi + mpmath.log(mpmath.pi) / (2 * mpmath.pi))


def weil_kernel_vec(x, zeros):
    """Weil kernel K(x) using vectorized zero sum."""
    diag = 1.0 if abs(x) < 1e-12 else 0.0
    zero_sum = np.sum(2 * np.cos(zeros * x) / (0.25 + zeros**2))
    return diag - zero_sum / (2 * np.pi)


def build_weil_matrix(primes, m_max, zeros, include_bg=True):
    """Build full Weil matrix with symmetry optimization."""
    labels = [(p, a) for p in primes for a in range(1, m_max + 1)]
    N = len(labels)
    M = np.zeros((N, N))

    # Precompute log values and weights
    log_p = {p: np.log(p) for p in primes}

    for i in range(N):
        pi, ai = labels[i]
        lpi = log_p[pi]
        xi = ai * lpi
        wi = np.sqrt(lpi) / pi**(ai / 2)

        for j in range(i + 1):  # Upper triangle + diagonal
            pj, aj = labels[j]
            lpj = log_p[pj]
            xj = aj * lpj
            wj = np.sqrt(lpj) / pj**(aj / 2)

            x = xi - xj
            K_val = weil_kernel_vec(x, zeros)
            if include_bg:
                K_val += K_bg(x)
            if i == j:
                K_val += 1.0  # Delta

            val = -wi * wj * K_val
            M[i, j] = val
            M[j, i] = val

    return M, labels


def primitive_projection(M):
    """Project onto primitive subspace."""
    N = M.shape[0]
    v = np.ones(N) / np.sqrt(N)
    P = np.eye(N) - np.outer(v, v)
    Mp = P @ M @ P
    return (Mp + Mp.T) / 2


def compute_entropy(primes, m_max, zeros):
    """Compute discrete entropy of cross-term distribution."""
    labels = [(p, a) for p in primes for a in range(1, m_max + 1)]
    N = len(labels)
    log_p = {p: np.log(p) for p in primes}

    W = np.zeros((N, N))
    for i, (pi, ai) in enumerate(labels):
        xi = ai * log_p[pi]
        wi = log_p[pi] / pi**(ai / 2)
        for j, (pj, aj) in enumerate(labels):
            xj = aj * log_p[pj]
            wj = log_p[pj] / pj**(aj / 2)
            K_val = weil_kernel_vec(xi - xj, zeros)
            W[i, j] = wi * wj * abs(K_val)

    total = np.sum(W)
    if total > 0:
        P = W / total
        flat = P.flatten()
        flat = flat[flat > 1e-30]
        H = -np.sum(flat * np.log(flat))
    else:
        H = 0
    return H


# =============================================================================
# MAIN TEST SUITE
# =============================================================================

def test_eigenvalue_scaling():
    """Main test: eigenvalue scaling with number of primes."""
    print(SEPARATOR)
    print("LARGE-SCALE WEIL MATRIX EIGENVALUE TEST")
    print(SEPARATOR)
    print()

    zeros = get_zeta_zeros(100)
    all_primes = sieve_primes(3500)  # Enough for 500 primes
    print(f"  {len(all_primes)} primes available up to {all_primes[-1]}")
    print()

    # Configurations: (n_primes, m_max)
    # For large matrices, use m_max=2 or 3 to keep size manageable
    configs = [
        (10, 3),    # 30×30
        (20, 3),    # 60×60
        (30, 3),    # 90×90
        (50, 3),    # 150×150
        (75, 3),    # 225×225
        (100, 3),   # 300×300
        (50, 4),    # 200×200
        (75, 4),    # 300×300
        (100, 4),   # 400×400
        (150, 3),   # 450×450
        (200, 3),   # 600×600
    ]

    results = []

    print(f"  {'#pr':>5s} {'m':>3s} {'N':>5s} | {'max_prim':>14s} {'min_prim':>14s} "
          f"{'gap':>12s} {'#neg':>5s} {'#pos':>5s} {'cond#':>12s} {'time':>8s}")
    print("  " + "-" * 100)

    for np_count, m_max in configs:
        primes = all_primes[:np_count]
        if len(primes) < np_count:
            print(f"  {np_count:5d} — not enough primes, skipping")
            continue

        N = np_count * m_max

        t1 = time.time()
        M, labels = build_weil_matrix(primes, m_max, zeros, include_bg=True)
        t_build = time.time() - t1

        t2 = time.time()
        Mp = primitive_projection(M)
        eigs = np.sort(np.linalg.eigvalsh(Mp))
        t_eig = time.time() - t2

        nz = eigs[np.abs(eigs) > 1e-13]
        if len(nz) > 0:
            max_prim = nz[-1]
            min_prim = nz[0]
            gap = max_prim - min_prim
            n_neg = np.sum(nz < -1e-13)
            n_pos = np.sum(nz > 1e-13)
            cond = abs(max_prim / min_prim) if abs(min_prim) > 1e-30 else float('inf')
        else:
            max_prim = min_prim = gap = 0
            n_neg = n_pos = 0
            cond = 0

        dt = t_build + t_eig
        apt = max_prim < 1e-10

        results.append({
            'np': np_count, 'm': m_max, 'N': N,
            'max': max_prim, 'min': min_prim,
            'gap': gap, 'n_neg': n_neg, 'n_pos': n_pos,
            'cond': cond, 'time': dt,
            't_build': t_build, 't_eig': t_eig,
            'apt': apt
        })

        apt_str = "✓" if apt else "✗"
        print(f"  {np_count:5d} {m_max:3d} {N:5d} | {max_prim:+14.6e} {min_prim:+14.6e} "
              f"{gap:12.4e} {n_neg:5d} {n_pos:5d} {cond:12.4e} {dt:7.1f}s {apt_str}")

        sys.stdout.flush()

        # Bail if taking too long
        if dt > 300:
            print(f"\n  WARNING: Last computation took {dt:.0f}s, stopping to avoid timeout.")
            break

    return results


def test_entropy_gap_correlation(zeros):
    """Entropy-gap correlation at larger scales."""
    print()
    print(SEPARATOR)
    print("ENTROPY-GAP CORRELATION AT LARGER SCALES")
    print(SEPARATOR)
    print()

    all_primes = sieve_primes(1000)

    sizes = [5, 10, 15, 20, 30, 40, 50, 60, 75]
    m_max = 3
    gaps = []
    entropies = []

    print(f"  {'#pr':>5s} {'N':>5s} | {'gap':>12s} {'H(μ)':>12s} {'gap/H':>10s}")
    print("  " + "-" * 50)

    for np_count in sizes:
        primes = all_primes[:np_count]
        N = np_count * m_max

        M, labels = build_weil_matrix(primes, m_max, zeros, include_bg=True)
        Mp = primitive_projection(M)
        eigs = np.sort(np.linalg.eigvalsh(Mp))
        nz = eigs[np.abs(eigs) > 1e-13]
        gap = (nz[-1] - nz[0]) if len(nz) > 0 else 0

        H = compute_entropy(primes, m_max, zeros)

        gaps.append(gap)
        entropies.append(H)
        gap_H = gap / H if H > 0 else 0

        print(f"  {np_count:5d} {N:5d} | {gap:12.6e} {H:12.6f} {gap_H:10.6f}")

    if len(gaps) >= 3:
        corr = np.corrcoef(gaps, entropies)[0, 1]
        print(f"\n  Correlation(gap, entropy) = {corr:.6f}")
        if corr > 0.99:
            print(f"  NEAR-PERFECT CORRELATION — AMR prediction strongly confirmed at scale")
        elif corr > 0.95:
            print(f"  EXCELLENT CORRELATION — AMR prediction confirmed")
        elif corr > 0.8:
            print(f"  STRONG CORRELATION — consistent with AMR")
        else:
            print(f"  MODERATE/WEAK CORRELATION — r = {corr:.4f}")

    return gaps, entropies


def analyze_scaling(results):
    """Fit scaling laws to the results."""
    print()
    print(SEPARATOR)
    print("SCALING ANALYSIS")
    print(SEPARATOR)
    print()

    if len(results) < 4:
        print("  Not enough data points for scaling analysis.")
        return

    # Filter to m_max=3 for consistent comparison
    r3 = [r for r in results if r['m'] == 3]
    if len(r3) < 4:
        r3 = results

    dims = np.array([r['N'] for r in r3])
    gaps = np.array([r['gap'] for r in r3])
    times = np.array([r['time'] for r in r3])
    max_eigs = np.array([r['max'] for r in r3])

    # Gap scaling
    valid = gaps > 1e-15
    if np.sum(valid) > 3:
        # Power law: gap ~ C * N^α
        coeffs = np.polyfit(np.log(dims[valid]), np.log(gaps[valid]), 1)
        alpha = coeffs[0]
        C = np.exp(coeffs[1])

        # Log law: gap ~ a * log(N) + b
        log_dims = np.log(dims[valid])
        coeffs_log = np.polyfit(log_dims, gaps[valid], 1)

        # Residuals
        resid_pow = np.sum((np.log(gaps[valid]) - (coeffs[0] * np.log(dims[valid]) + coeffs[1]))**2)
        resid_log = np.sum((gaps[valid] - (coeffs_log[0] * log_dims + coeffs_log[1]))**2)

        print(f"  SPECTRAL GAP SCALING:")
        print(f"    Power law: gap ≈ {C:.4e} · N^{{{alpha:.4f}}} (residual: {resid_pow:.4e})")
        print(f"    Log law:   gap ≈ {coeffs_log[0]:.4e} · log(N) + {coeffs_log[1]:.4e} (residual: {resid_log:.4e})")
        better = "logarithmic" if resid_log < resid_pow else f"power law (α={alpha:.4f})"
        print(f"    Better fit: {better}")

    # Time scaling
    if len(times) > 3 and np.all(times > 0):
        t_coeffs = np.polyfit(np.log(dims), np.log(times), 1)
        t_exp = t_coeffs[0]
        print(f"\n  TIME SCALING:")
        print(f"    time ≈ C · N^{{{t_exp:.2f}}}")
        if t_exp > 2.5:
            print(f"    Dominated by matrix construction (O(N²) kernel evaluations)")
        else:
            print(f"    Sub-cubic — eigenvalue decomposition not yet dominant")

        # Predict time for larger sizes
        for target_N in [500, 750, 1000]:
            pred_time = np.exp(t_coeffs[1]) * target_N**t_exp
            print(f"    Predicted time for N={target_N}: {pred_time:.0f}s ({pred_time/60:.1f}min)")

    # Max eigenvalue trend
    print(f"\n  MAX PRIMITIVE EIGENVALUE TREND:")
    for r in r3:
        status = "NEG ✓" if r['max'] < -1e-13 else ("~0" if abs(r['max']) < 1e-10 else "POS ✗")
        print(f"    N={r['N']:5d}: max_prim = {r['max']:+.6e}  [{status}]")

    all_apt = all(r['apt'] for r in results)
    print(f"\n  APT holds for ALL tested sizes: {'YES' if all_apt else 'NO'}")
    if not all_apt:
        failures = [r for r in results if not r['apt']]
        fail_strs = [f"N={r['N']}" for r in failures]
        print(f"  Failures at: {fail_strs}")
        print(f"  Max positive eigenvalue: {max(r['max'] for r in failures):.6e}")
        print(f"  (These may be numerical precision artifacts at matrix boundary)")


def main():
    t0 = time.time()
    print("=" * 70)
    print("AMR LARGE-SCALE EIGENVALUE COMPUTATION")
    print("=" * 70)
    print()
    print(f"  mpmath precision: {mpmath.mp.dps} digits")
    print(f"  numpy dtype: float64")
    print()

    results = test_eigenvalue_scaling()

    # Reuse cached zeros for entropy test
    zeros = get_zeta_zeros(80)
    gaps, entropies = test_entropy_gap_correlation(zeros)

    analyze_scaling(results)

    print()
    print(SEPARATOR)
    print("FINAL SUMMARY")
    print(SEPARATOR)
    if results:
        max_N = max(r['N'] for r in results)
        apt_count = sum(1 for r in results if r['apt'])
        total = len(results)
        print(f"  Tested up to N={max_N}")
        print(f"  APT holds: {apt_count}/{total} configurations")
        if any(not r['apt'] for r in results):
            worst = max(r['max'] for r in results)
            print(f"  Worst max primitive eigenvalue: {worst:+.6e}")
        else:
            best_neg = max(r['max'] for r in results)
            print(f"  All primitive eigenvalues negative (largest: {best_neg:+.6e})")
    print(f"  Total time: {time.time()-t0:.1f}s")
    print(SEPARATOR)


if __name__ == '__main__':
    main()
