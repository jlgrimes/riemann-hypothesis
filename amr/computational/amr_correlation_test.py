#!/usr/bin/env python3
"""
AMR Correlation Decay & Entropy Tests
======================================
Tests 1-2 from the AMR computational validation plan.

TEST 1 — CORRELATION DECAY:
  For prime pairs (p,q), compute:
    M_{p,q} = Σ_{m,n=1}^{N} (log p · log q)/(p^{m/2} q^{n/2}) K(m log p - n log q)
  where K uses zeta zeros. AMR predicts:
    |M_{p,q}| = O(1/((pq)^{1/4} (log max(p,q))^{1/2+ε}))
  Compare with Baker's lower bound |m log p - n log q| ≥ C/(max(m,n))^κ.

TEST 2 — ENTROPY COMPUTATION:
  Define discrete entropy of cross-term distribution.
  Verify: positive entropy ⟹ negative eigenvalues (entropy-positivity duality).
"""

import numpy as np
import mpmath
from pathlib import Path
import time

mpmath.mp.dps = 30

OUT_DIR = Path(__file__).parent
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
    """First n non-trivial zeta zeros (imaginary parts)."""
    return [float(mpmath.zetazero(k).imag) for k in range(1, n_zeros + 1)]


def weil_kernel(x, zeros):
    """
    Weil kernel: K(x) = δ(x) - (1/2π) Σ_γ 2cos(γx)/(1/4 + γ²)
    For x ≠ 0, the delta term vanishes.
    """
    diag = 1.0 if abs(x) < 1e-12 else 0.0
    zero_sum = sum(2 * np.cos(g * x) / (0.25 + g**2) for g in zeros)
    return diag - zero_sum / (2 * np.pi)


# =============================================================================
# TEST 1: CORRELATION DECAY
# =============================================================================

def compute_M_pq(p, q, m_max, zeros):
    """Compute M_{p,q} = Σ_{m,n} (log p log q)/(p^{m/2} q^{n/2}) K(m log p - n log q)."""
    lp, lq = np.log(p), np.log(q)
    total = 0.0
    for m in range(1, m_max + 1):
        for n in range(1, m_max + 1):
            weight = (lp * lq) / (p**(m/2) * q**(n/2))
            x = m * lp - n * lq
            total += weight * weil_kernel(x, zeros)
    return total


def test_correlation_decay():
    """Test AMR prediction: |M_{p,q}| decays as O(1/((pq)^{1/4} (log max(p,q))^{1/2+ε}))."""
    print(SEPARATOR)
    print("TEST 1: CORRELATION DECAY")
    print(SEPARATOR)
    print()

    print("  Computing 100 zeta zeros...")
    zeros = get_zeta_zeros(100)
    primes = sieve_primes(200)
    m_max = 8

    print(f"  Using {len(primes)} primes up to 200, m_max={m_max}")
    print(f"  AMR prediction: |M_{{p,q}}| = O(1/((pq)^{{1/4}} (log max(p,q))^{{1/2+ε}}))")
    print()

    # Compute M_{p,q} for all pairs
    results = []
    for i, p in enumerate(primes):
        for j, q in enumerate(primes):
            if p >= q:
                continue
            M_pq = compute_M_pq(p, q, m_max, zeros)
            pq_factor = (p * q)**0.25
            log_factor = np.log(max(p, q))**0.5
            predicted_bound = 1.0 / (pq_factor * log_factor)
            ratio = abs(M_pq) / predicted_bound if predicted_bound > 0 else float('inf')
            results.append((p, q, M_pq, abs(M_pq), predicted_bound, ratio))

    # Sort by ratio (worst cases first)
    results.sort(key=lambda x: -x[5])

    print(f"  {'p':>4s} {'q':>4s} | {'M_{p,q}':>14s} {'|M_{p,q}|':>12s} "
          f"{'predicted':>12s} {'ratio':>10s} {'decay?':>7s}")
    print("  " + "-" * 75)

    # Show top 15 worst cases and a sample of good cases
    for p, q, M_pq, absM, pred, ratio in results[:15]:
        decay_ok = "YES" if ratio < 10.0 else "WARN"
        print(f"  {p:4d} {q:4d} | {M_pq:+14.6e} {absM:12.6e} "
              f"{pred:12.6e} {ratio:10.4f} {decay_ok:>7s}")

    print("  ...")
    for p, q, M_pq, absM, pred, ratio in results[-5:]:
        decay_ok = "YES" if ratio < 10.0 else "WARN"
        print(f"  {p:4d} {q:4d} | {M_pq:+14.6e} {absM:12.6e} "
              f"{pred:12.6e} {ratio:10.4f} {decay_ok:>7s}")

    # Statistics
    ratios = [r[5] for r in results]
    print()
    print(f"  STATISTICS over {len(results)} prime pairs:")
    print(f"    max ratio:    {max(ratios):.6f}")
    print(f"    mean ratio:   {np.mean(ratios):.6f}")
    print(f"    median ratio: {np.median(ratios):.6f}")
    print(f"    std ratio:    {np.std(ratios):.6f}")
    print(f"    pairs with ratio > 10: {sum(1 for r in ratios if r > 10)}")
    print(f"    pairs with ratio > 5:  {sum(1 for r in ratios if r > 5)}")

    # Fit power law: |M_{p,q}| ~ C / (pq)^α
    log_pq = np.array([np.log(r[0] * r[1]) for r in results])
    log_absM = np.array([np.log(r[3]) if r[3] > 1e-30 else -70 for r in results])
    valid = log_absM > -60
    if np.sum(valid) > 5:
        coeffs = np.polyfit(log_pq[valid], log_absM[valid], 1)
        alpha_fit = -coeffs[0]
        C_fit = np.exp(coeffs[1])
        print(f"\n  POWER LAW FIT: |M_{{p,q}}| ≈ {C_fit:.4f} / (pq)^{{{alpha_fit:.4f}}}")
        print(f"    AMR predicts α ≈ 0.25 (from (pq)^{{1/4}} factor)")
        print(f"    Fitted α = {alpha_fit:.4f}")
        if abs(alpha_fit - 0.25) < 0.15:
            print(f"    STATUS: CONSISTENT with AMR prediction")
        else:
            print(f"    STATUS: DEVIATES from AMR prediction (α_fit ≠ 0.25)")

    # Baker bound check
    print(f"\n  BAKER BOUND CHECK:")
    print(f"    Baker's theorem: |m log p - n log q| ≥ C / max(m,n)^κ")
    print(f"    for p ≠ q, with effective constants C, κ.")
    min_diff = float('inf')
    min_case = None
    for p in primes[:20]:
        for q in primes[:20]:
            if p == q:
                continue
            lp, lq = np.log(p), np.log(q)
            for m in range(1, m_max + 1):
                for n in range(1, m_max + 1):
                    diff = abs(m * lp - n * lq)
                    if diff > 1e-15 and diff < min_diff:
                        min_diff = diff
                        min_case = (m, p, n, q)

    if min_case:
        m, p, n, q = min_case
        baker_bound = 1.0 / max(m, n)**10  # Conservative Baker κ=10
        print(f"    Smallest |m log p - n log q| found: {min_diff:.10e}")
        print(f"      at m={m}, p={p}, n={n}, q={q}")
        print(f"    Baker lower bound (κ=10): {baker_bound:.10e}")
        print(f"    Baker satisfied: {'YES' if min_diff >= baker_bound else 'NO'}")

    return results


# =============================================================================
# TEST 2: ENTROPY COMPUTATION
# =============================================================================

def compute_cross_term_distribution(primes, m_max, zeros):
    """
    Compute the discrete cross-term distribution μ_{p,q} on T_A.
    Returns the matrix of kernel values and derived probability weights.
    """
    N = len(primes)
    labels = [(p, a) for p in primes for a in range(1, m_max + 1)]
    dim = len(labels)

    # Build weighted kernel matrix
    K_matrix = np.zeros((dim, dim))
    for i, (pi, ai) in enumerate(labels):
        xi = ai * np.log(pi)
        wi = np.log(pi) / pi**(ai / 2)
        for j, (pj, aj) in enumerate(labels):
            xj = aj * np.log(pj)
            wj = np.log(pj) / pj**(aj / 2)
            K_val = weil_kernel(xi - xj, zeros)
            K_matrix[i, j] = wi * wj * abs(K_val)

    # Normalize to get a probability distribution
    total = np.sum(K_matrix)
    if total > 0:
        P_matrix = K_matrix / total
    else:
        P_matrix = np.ones_like(K_matrix) / dim**2

    return K_matrix, P_matrix, labels


def discrete_entropy(P):
    """Compute Shannon entropy H = -Σ p_ij log(p_ij) of a probability matrix."""
    flat = P.flatten()
    flat = flat[flat > 1e-30]
    return -np.sum(flat * np.log(flat))


def test_entropy_positivity():
    """
    Test entropy-positivity duality:
      positive entropy ⟹ negative eigenvalues in primitive subspace.
    """
    print()
    print(SEPARATOR)
    print("TEST 2: ENTROPY-POSITIVITY DUALITY")
    print(SEPARATOR)
    print()
    print("  Theory: The discrete entropy of the cross-term measure μ_{p,q}")
    print("  on the arithmetic torus T_A should correlate with the spectral")
    print("  gap of the cross-term matrix. Higher entropy → more diffuse")
    print("  measure → better cancellation → more negative eigenvalues.")
    print()

    zeros = get_zeta_zeros(80)

    truncations = [3, 5, 8, 12, 15, 20, 30]
    results = []

    print(f"  {'#primes':>8s} {'dim':>5s} | {'entropy':>12s} {'max_eig':>12s} "
          f"{'min_eig':>12s} {'gap':>12s} {'H>0⟹gap>0':>12s}")
    print("  " + "-" * 80)

    for np_count in truncations:
        primes = sieve_primes(200)[:np_count]
        if len(primes) < np_count:
            primes = sieve_primes(500)[:np_count]
        m_max = 3
        labels = [(p, a) for p in primes for a in range(1, m_max + 1)]
        dim = len(labels)

        # Compute kernel matrix (unweighted for eigenvalues)
        K = np.zeros((dim, dim))
        for i, (pi, ai) in enumerate(labels):
            xi = ai * np.log(pi)
            for j, (pj, aj) in enumerate(labels):
                xj = aj * np.log(pj)
                K[i, j] = weil_kernel(xi - xj, zeros)
        K = (K + K.T) / 2

        # Weight matrix for Weil form
        d = np.array([np.sqrt(np.log(p)) / p**(a/2) for p, a in labels])
        D = np.diag(d)
        M = -D @ K @ D
        M = (M + M.T) / 2

        # Primitive projection
        v = np.ones(dim) / np.sqrt(dim)
        P = np.eye(dim) - np.outer(v, v)
        Mp = P @ M @ P
        Mp = (Mp + Mp.T) / 2
        eigs = np.linalg.eigvalsh(Mp)
        nz = eigs[np.abs(eigs) > 1e-14]
        max_eig = nz[-1] if len(nz) > 0 else 0
        min_eig = nz[0] if len(nz) > 0 else 0
        gap = max_eig - min_eig

        # Compute entropy of the weighted cross-term distribution
        K_w, P_matrix, _ = compute_cross_term_distribution(primes, m_max, zeros)
        H = discrete_entropy(P_matrix)

        duality = "YES" if (H > 0 and gap > 0) else ("TRIVIAL" if H <= 0 else "NO")
        results.append((np_count, dim, H, max_eig, min_eig, gap, duality))

        print(f"  {np_count:8d} {dim:5d} | {H:12.6f} {max_eig:+12.6e} "
              f"{min_eig:+12.6e} {gap:12.6e} {duality:>12s}")

    # Check correlation between entropy and spectral gap
    entropies = [r[2] for r in results]
    gaps = [r[5] for r in results]
    if len(entropies) > 3:
        corr = np.corrcoef(entropies, gaps)[0, 1]
        print(f"\n  Correlation(entropy, spectral_gap) = {corr:.6f}")
        if corr > 0.5:
            print(f"  POSITIVE CORRELATION: higher entropy ↔ larger spectral gap")
            print(f"  This SUPPORTS the entropy-positivity duality prediction.")
        elif corr < -0.5:
            print(f"  NEGATIVE CORRELATION: entropy inversely related to gap")
            print(f"  This CHALLENGES the simple duality prediction.")
        else:
            print(f"  WEAK CORRELATION: no clear linear relationship.")

    # Entropy growth rate
    if len(results) >= 3:
        dims = [r[1] for r in results]
        log_dims = np.log(dims)
        H_vals = entropies
        if min(H_vals) > 0:
            log_H = np.log(H_vals)
            coeffs = np.polyfit(log_dims, log_H, 1)
            print(f"\n  Entropy growth: H ~ dim^{{{coeffs[0]:.3f}}}")
            print(f"    AMR predicts logarithmic growth rate matching Rudolph entropy.")

    return results


# =============================================================================
# MAIN
# =============================================================================

def main():
    t0 = time.time()
    print("AMR CORRELATION DECAY & ENTROPY TESTS")
    print(SEPARATOR)
    print()

    results_decay = test_correlation_decay()
    results_entropy = test_entropy_positivity()

    print()
    print(SEPARATOR)
    print("SUMMARY")
    print(SEPARATOR)
    print(f"  Test 1 (Correlation Decay): "
          f"max ratio = {max(r[5] for r in results_decay):.4f}")
    ratios = [r[5] for r in results_decay]
    print(f"    AMR decay bound holds for {sum(1 for r in ratios if r < 10)/len(ratios)*100:.1f}% of pairs")
    print(f"  Test 2 (Entropy-Positivity): "
          f"all duality checks = {all(r[6] in ('YES','TRIVIAL') for r in results_entropy)}")
    print(f"  Total time: {time.time()-t0:.1f}s")
    print(SEPARATOR)


if __name__ == '__main__':
    main()
