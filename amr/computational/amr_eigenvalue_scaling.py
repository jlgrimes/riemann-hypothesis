#!/usr/bin/env python3
"""
AMR Large-Scale Eigenvalue Test
=================================
TEST 5: Extended eigenvalue analysis for AMR predictions.

1. Build cross-term matrices up to 300×300
2. Track spectral gap growth rate — AMR predicts log-growth
3. Compare with predictions from measure rigidity (rate should match entropy)
4. Test the transition from finite to asymptotic regime
"""

import numpy as np
import mpmath
import time

mpmath.mp.dps = 30

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
    """First n non-trivial zeta zeros (imaginary parts), cached."""
    return np.array([float(mpmath.zetazero(k).imag) for k in range(1, n_zeros + 1)])


def weil_kernel_vectorized(x, zeros):
    """Vectorized Weil kernel for a single x value."""
    diag = 1.0 if abs(x) < 1e-12 else 0.0
    zero_sum = np.sum(2 * np.cos(zeros * x) / (0.25 + zeros**2))
    return diag - zero_sum / (2 * np.pi)


def K_bg(x):
    """Background kernel from digamma function."""
    if abs(x) < 1e-14:
        x = 1e-12
    arg = mpmath.mpc(0.25, x / 2)
    psi = mpmath.digamma(arg)
    return float(-mpmath.re(psi) / mpmath.pi + mpmath.log(mpmath.pi) / (2 * mpmath.pi))


def build_weil_matrix(primes, m_max, zeros, include_bg=True):
    """
    Build the full Weil matrix M for the given primes and truncation.
    M_{ij} = -w_i w_j K(x_i - x_j)
    where w_i = sqrt(log p) / p^{a/2} and x_i = a log p.
    """
    labels = [(p, a) for p in primes for a in range(1, m_max + 1)]
    N = len(labels)
    M = np.zeros((N, N))

    for i, (pi, ai) in enumerate(labels):
        lpi = np.log(pi)
        xi = ai * lpi
        wi = np.sqrt(lpi) / pi**(ai / 2)
        for j, (pj, aj) in enumerate(labels):
            if j > i:
                break  # Fill upper triangle, then symmetrize
            lpj = np.log(pj)
            xj = aj * lpj
            wj = np.sqrt(lpj) / pj**(aj / 2)

            x = xi - xj
            K_val = weil_kernel_vectorized(x, zeros)
            if include_bg:
                K_val += K_bg(x)
            if i == j:
                K_val += 1.0  # Delta term

            M[i, j] = -wi * wj * K_val
            M[j, i] = M[i, j]

    return M, labels


def primitive_projection(M):
    """Project onto primitive subspace (orthogonal to constant vector)."""
    N = M.shape[0]
    v = np.ones(N) / np.sqrt(N)
    P = np.eye(N) - np.outer(v, v)
    Mp = P @ M @ P
    return (Mp + Mp.T) / 2


def spectral_analysis(M):
    """Compute eigenvalue statistics of a symmetric matrix."""
    eigs = np.sort(np.linalg.eigvalsh(M))
    nz = eigs[np.abs(eigs) > 1e-14]
    if len(nz) == 0:
        return {'max': 0, 'min': 0, 'gap': 0, 'n_pos': 0, 'n_neg': 0, 'trace': 0}
    return {
        'max': nz[-1],
        'min': nz[0],
        'gap': nz[-1] - nz[0],
        'n_pos': np.sum(nz > 1e-14),
        'n_neg': np.sum(nz < -1e-14),
        'trace': np.trace(M),
        'eigs': nz
    }


# =============================================================================
# TEST 5A: LARGE-SCALE EIGENVALUE GROWTH
# =============================================================================

def test_eigenvalue_growth():
    """
    Track spectral gap growth with matrix size.
    AMR predicts: gap ~ C · log(N) where N is the matrix dimension.
    """
    print(SEPARATOR)
    print("TEST 5A: EIGENVALUE GROWTH WITH MATRIX SIZE")
    print(SEPARATOR)
    print()
    print("  AMR predicts: spectral gap grows as O(log N)")
    print("  where N = #primes × m_max is the matrix dimension.")
    print()

    print("  Computing zeta zeros...")
    zeros = get_zeta_zeros(100)
    all_primes = sieve_primes(500)

    configs = [
        (5, 3),    # 15×15
        (10, 3),   # 30×30
        (15, 3),   # 45×45
        (20, 3),   # 60×60
        (25, 3),   # 75×75
        (15, 5),   # 75×75
        (20, 5),   # 100×100
        (30, 3),   # 90×90
        (25, 5),   # 125×125
        (30, 5),   # 150×150
        (40, 4),   # 160×160
        (50, 4),   # 200×200
    ]

    results = []
    print(f"  {'#pr':>4s} {'m':>3s} {'N':>5s} | {'max_prim':>12s} {'min_prim':>12s} "
          f"{'gap':>12s} {'log(N)':>8s} {'gap/log(N)':>12s} {'APT?':>5s}")
    print("  " + "-" * 85)

    for np_count, m_max in configs:
        primes = all_primes[:np_count]
        if len(primes) < np_count:
            continue

        t1 = time.time()
        M, labels = build_weil_matrix(primes, m_max, zeros, include_bg=True)
        N = M.shape[0]
        Mp = primitive_projection(M)
        spec = spectral_analysis(Mp)
        dt = time.time() - t1

        log_N = np.log(N)
        gap_ratio = spec['gap'] / log_N if log_N > 0 else 0
        apt = "YES" if spec['max'] < 1e-10 else "NO"

        results.append({
            'np': np_count, 'm': m_max, 'N': N,
            'max': spec['max'], 'min': spec['min'],
            'gap': spec['gap'], 'log_N': log_N,
            'gap_ratio': gap_ratio, 'apt': apt, 'time': dt
        })

        print(f"  {np_count:4d} {m_max:3d} {N:5d} | {spec['max']:+12.4e} {spec['min']:+12.4e} "
              f"{spec['gap']:12.4e} {log_N:8.3f} {gap_ratio:12.6f} {apt:>5s}")

    # Fit growth rate: gap ~ C * N^α or gap ~ C * log(N)
    if len(results) >= 5:
        dims = np.array([r['N'] for r in results])
        gaps = np.array([r['gap'] for r in results])
        log_dims = np.log(dims)

        # Power law fit: gap ~ C * N^α
        valid = gaps > 1e-15
        if np.sum(valid) > 3:
            coeffs_pow = np.polyfit(np.log(dims[valid]), np.log(gaps[valid]), 1)
            alpha = coeffs_pow[0]
            C_pow = np.exp(coeffs_pow[1])

            # Log fit: gap ~ a * log(N) + b
            coeffs_log = np.polyfit(log_dims[valid], gaps[valid], 1)
            a_log, b_log = coeffs_log

            # Which fits better? (residuals)
            resid_pow = np.sum((np.log(gaps[valid]) - (coeffs_pow[0] * np.log(dims[valid]) + coeffs_pow[1]))**2)
            resid_log = np.sum((gaps[valid] - (a_log * log_dims[valid] + b_log))**2)

            print(f"\n  FIT ANALYSIS:")
            print(f"    Power law: gap ≈ {C_pow:.4e} · N^{{{alpha:.4f}}} (residual: {resid_pow:.4e})")
            print(f"    Log law:   gap ≈ {a_log:.4e} · log(N) + {b_log:.4e} (residual: {resid_log:.4e})")
            if resid_log < resid_pow:
                print(f"    BETTER FIT: logarithmic — CONSISTENT with AMR prediction")
            else:
                print(f"    BETTER FIT: power law (α = {alpha:.4f})")
                print(f"    May indicate sub-logarithmic or polynomial regime")

    return results


# =============================================================================
# TEST 5B: SPECTRAL GAP vs ENTROPY
# =============================================================================

def test_spectral_gap_entropy():
    """
    Compare spectral gap with discrete entropy of cross-term distribution.
    AMR predicts: spectral gap should track entropy of the ×p,×q action.
    """
    print()
    print(SEPARATOR)
    print("TEST 5B: SPECTRAL GAP vs ENTROPY COMPARISON")
    print(SEPARATOR)
    print()

    zeros = get_zeta_zeros(80)
    all_primes = sieve_primes(200)

    results = []
    print(f"  {'#pr':>4s} {'N':>5s} | {'gap':>12s} {'H(μ)':>12s} "
          f"{'gap/H':>10s} {'H/log(N)':>10s}")
    print("  " + "-" * 60)

    for np_count in [5, 8, 10, 15, 20, 25, 30]:
        primes = all_primes[:np_count]
        m_max = 3
        labels = [(p, a) for p in primes for a in range(1, m_max + 1)]
        N = len(labels)

        # Build and analyze Weil matrix
        M, _ = build_weil_matrix(primes, m_max, zeros, include_bg=True)
        Mp = primitive_projection(M)
        spec = spectral_analysis(Mp)

        # Compute entropy from cross-term distribution
        # Use |K(x_i - x_j)| weighted by w_i w_j as the measure
        W = np.zeros((N, N))
        for i, (pi, ai) in enumerate(labels):
            xi = ai * np.log(pi)
            wi = np.log(pi) / pi**(ai / 2)
            for j, (pj, aj) in enumerate(labels):
                xj = aj * np.log(pj)
                wj = np.log(pj) / pj**(aj / 2)
                K_val = weil_kernel_vectorized(xi - xj, zeros)
                W[i, j] = wi * wj * abs(K_val)

        total_W = np.sum(W)
        if total_W > 0:
            P = W / total_W
            flat = P.flatten()
            flat = flat[flat > 1e-30]
            H = -np.sum(flat * np.log(flat))
        else:
            H = 0

        gap = spec['gap']
        gap_H = gap / H if H > 0 else 0
        H_logN = H / np.log(N) if N > 1 else 0

        results.append((np_count, N, gap, H, gap_H, H_logN))
        print(f"  {np_count:4d} {N:5d} | {gap:12.6e} {H:12.6f} "
              f"{gap_H:10.6f} {H_logN:10.6f}")

    # Correlation
    gaps = [r[2] for r in results]
    entropies = [r[3] for r in results]
    if len(gaps) >= 3:
        corr = np.corrcoef(gaps, entropies)[0, 1]
        print(f"\n  Correlation(gap, entropy) = {corr:.6f}")
        if corr > 0.7:
            print(f"  STRONG POSITIVE CORRELATION — gap tracks entropy")
            print(f"  This confirms AMR prediction: spectral structure")
            print(f"  is governed by measure-theoretic entropy.")
        elif corr > 0.3:
            print(f"  MODERATE CORRELATION — partial agreement with AMR.")
        else:
            print(f"  WEAK/NEGATIVE CORRELATION — challenges simple AMR prediction.")

    return results


# =============================================================================
# TEST 5C: EIGENVALUE DISTRIBUTION SHAPE
# =============================================================================

def test_eigenvalue_distribution():
    """
    Analyze the shape of the eigenvalue distribution.
    For large matrices, AMR predicts the spectral density should
    reflect the measure-theoretic properties of ×p,×q dynamics.
    """
    print()
    print(SEPARATOR)
    print("TEST 5C: EIGENVALUE DISTRIBUTION SHAPE")
    print(SEPARATOR)
    print()

    zeros = get_zeta_zeros(100)
    primes = sieve_primes(500)[:40]
    m_max = 4

    print(f"  Building {len(primes)*m_max}×{len(primes)*m_max} Weil matrix...")
    t1 = time.time()
    M, labels = build_weil_matrix(primes, m_max, zeros, include_bg=True)
    N = M.shape[0]
    Mp = primitive_projection(M)
    eigs = np.sort(np.linalg.eigvalsh(Mp))
    nz = eigs[np.abs(eigs) > 1e-14]
    dt = time.time() - t1
    print(f"  Done in {dt:.1f}s. N={N}, {len(nz)} non-zero eigenvalues.")
    print()

    # Distribution statistics
    if len(nz) > 0:
        print(f"  EIGENVALUE STATISTICS (primitive subspace):")
        print(f"    max:     {nz[-1]:+.8e}")
        print(f"    min:     {nz[0]:+.8e}")
        print(f"    mean:    {np.mean(nz):+.8e}")
        print(f"    median:  {np.median(nz):+.8e}")
        print(f"    std:     {np.std(nz):.8e}")
        print(f"    skew:    {float(np.mean(((nz - np.mean(nz))/np.std(nz))**3)):.6f}")
        print(f"    #pos:    {np.sum(nz > 1e-14)}")
        print(f"    #neg:    {np.sum(nz < -1e-14)}")
        print()

        # Quartile analysis
        q1, q2, q3 = np.percentile(nz, [25, 50, 75])
        print(f"    Q1 (25%):  {q1:+.8e}")
        print(f"    Q2 (50%):  {q2:+.8e}")
        print(f"    Q3 (75%):  {q3:+.8e}")
        print(f"    IQR:       {q3-q1:.8e}")
        print()

        # Check if distribution is predominantly negative (APT)
        neg_fraction = np.sum(nz < 0) / len(nz)
        print(f"    Fraction negative: {neg_fraction:.4f}")
        if neg_fraction > 0.9:
            print(f"    STRONGLY NEGATIVE — consistent with APT")
        elif neg_fraction > 0.5:
            print(f"    MAJORITY NEGATIVE — partial APT support")
        else:
            print(f"    NOT PREDOMINANTLY NEGATIVE — APT not evident at this scale")

        # Largest eigenvalues (most relevant for APT failure)
        print(f"\n    Top 10 eigenvalues (most positive):")
        for i, e in enumerate(nz[-10:]):
            print(f"      [{len(nz)-10+i:3d}] {e:+.8e}")

        print(f"\n    Bottom 5 eigenvalues (most negative):")
        for i, e in enumerate(nz[:5]):
            print(f"      [{i:3d}] {e:+.8e}")

    return nz


# =============================================================================
# TEST 5D: APT SENSITIVITY ANALYSIS
# =============================================================================

def test_apt_sensitivity():
    """
    Sensitivity analysis: how does the maximum primitive eigenvalue
    change with the number of zeros used in the kernel?
    This tests the stability of APT predictions.
    """
    print()
    print(SEPARATOR)
    print("TEST 5D: APT SENSITIVITY TO ZERO TRUNCATION")
    print(SEPARATOR)
    print()
    print("  Does the maximum primitive eigenvalue stabilize as we")
    print("  include more zeta zeros in the kernel?")
    print()

    all_zeros = get_zeta_zeros(100)
    primes = sieve_primes(200)[:20]
    m_max = 3

    print(f"  Matrix size: {len(primes)*m_max}×{len(primes)*m_max}")
    print(f"  {'#zeros':>7s} | {'max_prim':>14s} {'min_prim':>14s} "
          f"{'gap':>14s} {'Δmax_prim':>14s}")
    print("  " + "-" * 70)

    prev_max = None
    for nz in [10, 20, 30, 50, 70, 100]:
        zeros = all_zeros[:nz]
        M, labels = build_weil_matrix(primes, m_max, zeros, include_bg=True)
        Mp = primitive_projection(M)
        spec = spectral_analysis(Mp)

        delta = spec['max'] - prev_max if prev_max is not None else 0
        prev_max = spec['max']

        print(f"  {nz:7d} | {spec['max']:+14.6e} {spec['min']:+14.6e} "
              f"{spec['gap']:14.6e} {delta:+14.6e}")

    print()
    print("  If Δmax_prim → 0, the truncated computation converges,")
    print("  suggesting the finite approximation captures the essential spectral behavior.")


# =============================================================================
# MAIN
# =============================================================================

def main():
    t0 = time.time()
    print("AMR LARGE-SCALE EIGENVALUE TESTS")
    print(SEPARATOR)
    print()

    results_growth = test_eigenvalue_growth()
    results_entropy = test_spectral_gap_entropy()
    eig_dist = test_eigenvalue_distribution()
    test_apt_sensitivity()

    print()
    print(SEPARATOR)
    print("FINAL SUMMARY")
    print(SEPARATOR)
    print(f"  Eigenvalue growth: tested up to N={max(r['N'] for r in results_growth)}")
    apt_count = sum(1 for r in results_growth if r['apt'] == 'YES')
    total_count = len(results_growth)
    print(f"  APT holds: {apt_count}/{total_count} configurations")
    if len(eig_dist) > 0:
        print(f"  Largest eigenvalue distribution: max={eig_dist[-1]:+.6e}, min={eig_dist[0]:+.6e}")
    print(f"  Total time: {time.time()-t0:.1f}s")
    print(SEPARATOR)


if __name__ == '__main__':
    main()
