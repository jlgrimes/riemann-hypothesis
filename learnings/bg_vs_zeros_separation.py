#!/usr/bin/env python3
"""
BG vs Zeros Separation at Scale
================================
Extends bg_only_test.py with the CORRECT primitive projection to larger scales.

For P_0 in {200, 300, 500, 750, 1000}, build THREE separate Weil matrices:
  1. M_bg:   uses only K_bg(x) + delta(x)       — no zeros
  2. M_zeros: uses only K_zeros(x)               — just the zero sum kernel
  3. M_full:  uses K_bg(x) + K_zeros(x) + delta(x) — everything

Key question: Does the M_bg spectral gap grow with P_0, swamping the zeros contribution?
"""

import mpmath
import numpy as np
import time
import os

mpmath.mp.dps = 30


# ────────────────────────── primes ──────────────────────────

def sieve_primes(bound):
    sieve = [True] * (int(bound) + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(bound**0.5) + 1):
        if sieve[i]:
            for j in range(i*i, int(bound)+1, i):
                sieve[j] = False
    return [p for p in range(2, int(bound)+1) if sieve[p]]


# ────────────────────────── kernels ──────────────────────────

def K_bg(x):
    """K_bg(x) = -(1/pi) Re[psi(1/4 + ix/2)] + log(pi)/(2pi)"""
    if abs(x) < 1e-14:
        x = 1e-12
    arg = mpmath.mpc(0.25, x / 2)
    psi = mpmath.digamma(arg)
    return float(-mpmath.re(psi) / mpmath.pi + mpmath.log(mpmath.pi) / (2 * mpmath.pi))


def K_zeros(x, zeros):
    """K_zeros(x) = (1/pi) sum_gamma cos(gamma*x) / (1/4 + gamma^2)"""
    return float(np.sum(np.cos(zeros * x) / (0.25 + zeros**2)) / np.pi)


# ────────────────────────── matrix builders ──────────────────────────

def _make_labels(primes, m_max):
    return [(p, a) for p in primes for a in range(1, m_max + 1)]


def build_bg_matrix(primes, m_max):
    """M_bg: uses only delta + K_bg (no zeros)."""
    labels = _make_labels(primes, m_max)
    N = len(labels)
    M = np.zeros((N, N))
    for i, (p, a) in enumerate(labels):
        lp = np.log(p)
        for j, (q, b) in enumerate(labels):
            lq = np.log(q)
            x = a * lp - b * lq
            is_diag = (i == j)
            kb = K_bg(x)
            K_val = (1.0 if is_diag else 0.0) + kb
            M[i, j] = -np.sqrt(lp * lq) / (p**(a/2) * q**(b/2)) * K_val
    return M, labels


def build_zeros_matrix(primes, m_max, zeros_float):
    """M_zeros: uses only K_zeros (no delta, no K_bg)."""
    labels = _make_labels(primes, m_max)
    N = len(labels)
    zeros_arr = np.array(zeros_float)
    M = np.zeros((N, N))
    for i, (p, a) in enumerate(labels):
        lp = np.log(p)
        for j, (q, b) in enumerate(labels):
            lq = np.log(q)
            x = a * lp - b * lq
            kz = K_zeros(x, zeros_arr)
            M[i, j] = -np.sqrt(lp * lq) / (p**(a/2) * q**(b/2)) * kz
    return M, labels


def build_full_matrix(primes, m_max, zeros_float):
    """M_full: uses delta + K_bg + K_zeros (everything)."""
    labels = _make_labels(primes, m_max)
    N = len(labels)
    zeros_arr = np.array(zeros_float)
    M = np.zeros((N, N))
    for i, (p, a) in enumerate(labels):
        lp = np.log(p)
        for j, (q, b) in enumerate(labels):
            lq = np.log(q)
            x = a * lp - b * lq
            is_diag = (i == j)
            kb = K_bg(x)
            kz = K_zeros(x, zeros_arr)
            K_val = (1.0 if is_diag else 0.0) + kb + kz
            M[i, j] = -np.sqrt(lp * lq) / (p**(a/2) * q**(b/2)) * K_val
    return M, labels


# ────────────────────────── primitive projection ──────────────────────────

def primitive_projection(M, labels):
    """Project M onto the primitive subspace using the CORRECT projection."""
    N = M.shape[0]
    v = np.array([np.sqrt(np.log(p)) / p**(a / 2.0) for p, a in labels])
    v = v / np.linalg.norm(v)
    P = np.eye(N) - np.outer(v, v)
    Mp = P @ M @ P
    return (Mp + Mp.T) / 2


def primitive_eigenvalues(Mp):
    """Compute eigenvalues on primitive subspace, removing the trivial zero."""
    eigs = np.sort(np.linalg.eigvalsh(Mp))
    triv = np.argmin(np.abs(eigs))
    return np.delete(eigs, triv)


# ────────────────────────── main ──────────────────────────

def main():
    t_start = time.time()

    print("=" * 80)
    print("  BG vs ZEROS SEPARATION — Scaling Analysis")
    print("  Does the M_bg spectral gap grow with P_0, swamping the zeros?")
    print("=" * 80)
    print()

    # ── Compute zeta zeros ──
    print("Computing 200 zeta zeros via mpmath.zetazero ...")
    t0 = time.time()
    zeros_float = []
    for k in range(1, 201):
        zeros_float.append(float(mpmath.im(mpmath.zetazero(k))))
        if k % 50 == 0:
            print(f"  ... {k}/200 zeros computed ({time.time()-t0:.1f}s)")
    print(f"  Done: 200 zeros in {time.time()-t0:.1f}s  "
          f"(gamma_1 = {zeros_float[0]:.6f}, gamma_200 = {zeros_float[-1]:.6f})")
    print()

    # ── Test configurations ──
    configs = [
        (200,  3),
        (300,  3),
        (500,  2),
        (750,  2),
        (1000, 2),
    ]

    results = []

    # ── Header ──
    hdr = (f"{'P0':>5s} {'m_max':>5s} {'|primes|':>8s} {'N':>5s} | "
           f"{'bg gap':>12s} {'zeros norm':>12s} {'full gap':>12s} | "
           f"{'bg/zeros':>10s} {'bg > zeros?':>12s}")
    print(hdr)
    print("-" * len(hdr))

    for P0, m_max in configs:
        primes = sieve_primes(P0)
        N = len(primes) * m_max
        print(f"\n>>> P0={P0}, m_max={m_max}, |primes|={len(primes)}, N={N}")

        # ── M_bg ──
        print(f"    Building M_bg ({N}x{N}) ...")
        t0 = time.time()
        M_bg, labels = build_bg_matrix(primes, m_max)
        print(f"    M_bg done in {time.time()-t0:.1f}s")

        Mp_bg = primitive_projection(M_bg, labels)
        eigs_bg = primitive_eigenvalues(Mp_bg)
        bg_gap = float(np.max(np.abs(eigs_bg)))  # spectral gap = |max eigenvalue|
        bg_max = float(eigs_bg[-1])
        bg_min = float(eigs_bg[0])

        # ── M_zeros ──
        print(f"    Building M_zeros ({N}x{N}) ...")
        t0 = time.time()
        M_zeros, _ = build_zeros_matrix(primes, m_max, zeros_float)
        print(f"    M_zeros done in {time.time()-t0:.1f}s")

        Mp_zeros = primitive_projection(M_zeros, labels)
        eigs_zeros = primitive_eigenvalues(Mp_zeros)
        zeros_norm = float(np.max(np.abs(eigs_zeros)))  # operator norm

        # ── M_full ──
        print(f"    Building M_full ({N}x{N}) ...")
        t0 = time.time()
        M_full, _ = build_full_matrix(primes, m_max, zeros_float)
        print(f"    M_full done in {time.time()-t0:.1f}s")

        Mp_full = primitive_projection(M_full, labels)
        eigs_full = primitive_eigenvalues(Mp_full)
        full_gap = float(np.max(np.abs(eigs_full)))
        full_max = float(eigs_full[-1])
        full_min = float(eigs_full[0])

        # ── Ratio and dominance ──
        ratio = bg_gap / zeros_norm if zeros_norm > 1e-30 else float('inf')
        bg_dominates = bg_gap > zeros_norm

        row = {
            'P0': P0, 'm_max': m_max, 'n_primes': len(primes), 'N': N,
            'bg_gap': bg_gap, 'bg_max': bg_max, 'bg_min': bg_min,
            'zeros_norm': zeros_norm,
            'full_gap': full_gap, 'full_max': full_max, 'full_min': full_min,
            'ratio': ratio, 'bg_dominates': bg_dominates,
        }
        results.append(row)

        tag = "YES" if bg_dominates else "NO"
        print(f"  {P0:5d} {m_max:5d} {len(primes):8d} {N:5d} | "
              f"{bg_gap:12.6e} {zeros_norm:12.6e} {full_gap:12.6e} | "
              f"{ratio:10.2f} {tag:>12s}")

    # ── Summary table ──
    print("\n")
    print("=" * 80)
    print("  SUMMARY TABLE")
    print("=" * 80)
    print()
    fmt_h = (f"{'P0':>5s} {'N':>5s} | {'bg gap':>12s} {'zeros norm':>12s} "
             f"{'full gap':>12s} | {'ratio':>8s} {'dom?':>5s} | "
             f"{'bg_max':>12s} {'full_max':>12s}")
    print(fmt_h)
    print("-" * len(fmt_h))
    for r in results:
        tag = "YES" if r['bg_dominates'] else "NO"
        print(f"{r['P0']:5d} {r['N']:5d} | {r['bg_gap']:12.6e} {r['zeros_norm']:12.6e} "
              f"{r['full_gap']:12.6e} | {r['ratio']:8.2f} {tag:>5s} | "
              f"{r['bg_max']:+12.6e} {r['full_max']:+12.6e}")

    # ── Growth analysis ──
    print("\n")
    print("=" * 80)
    print("  GROWTH ANALYSIS")
    print("=" * 80)
    print()

    if len(results) >= 2:
        bg_gaps = [r['bg_gap'] for r in results]
        zeros_norms = [r['zeros_norm'] for r in results]
        P0s = [r['P0'] for r in results]

        print("  P0 transitions:")
        for i in range(1, len(results)):
            bg_growth = bg_gaps[i] / bg_gaps[i-1] if bg_gaps[i-1] > 0 else float('inf')
            z_growth = zeros_norms[i] / zeros_norms[i-1] if zeros_norms[i-1] > 0 else float('inf')
            print(f"    P0 {P0s[i-1]} -> {P0s[i]}: "
                  f"bg_gap x{bg_growth:.3f}, zeros_norm x{z_growth:.3f}")

    # ── Verdict ──
    all_dominant = all(r['bg_dominates'] for r in results)
    ratios_growing = True
    for i in range(1, len(results)):
        if results[i]['ratio'] < results[i-1]['ratio']:
            ratios_growing = False

    print("\n  VERDICT:")
    if all_dominant and ratios_growing:
        print("    The M_bg spectral gap DOMINATES the M_zeros operator norm at all tested")
        print("    scales, and the dominance GROWS with P0.")
        print("    => The background kernel swamps the zeros contribution.")
    elif all_dominant:
        print("    The M_bg spectral gap dominates at all tested scales,")
        print("    but the dominance ratio is NOT monotonically growing.")
    else:
        dominant_count = sum(1 for r in results if r['bg_dominates'])
        print(f"    M_bg dominates at {dominant_count}/{len(results)} tested scales.")
        if not results[-1]['bg_dominates']:
            print("    At the LARGEST scale, zeros contribution is NOT swamped.")

    t_total = time.time() - t_start
    print(f"\n  Total runtime: {t_total:.1f}s")

    # ── Save markdown ──
    md_path = os.path.join(os.path.dirname(__file__),
                           "bg_vs_zeros_separation_results.md")
    save_markdown(results, md_path, all_dominant, ratios_growing)
    print(f"\n  Results saved to: {md_path}")


def save_markdown(results, path, all_dominant, ratios_growing):
    lines = []
    w = lines.append

    w("# BG vs Zeros Separation — Scaling Analysis")
    w("")
    w("**Question**: Does the M_bg spectral gap grow with P_0, swamping the zeros contribution?")
    w("")
    w("## Setup")
    w("")
    w("Three matrices built for each P_0:")
    w("- **M_bg**: delta(x) + K_bg(x) only (no zeros)")
    w("- **M_zeros**: K_zeros(x) only (just the zero-sum kernel, 200 zeros)")
    w("- **M_full**: delta(x) + K_bg(x) + K_zeros(x) (everything)")
    w("")
    w("All eigenvalues computed on the **primitive subspace** using the correct projection.")
    w("")
    w("## Results")
    w("")
    w("| P_0 | N | bg gap | zeros norm | full gap | ratio | bg dominates? | bg max prim | full max prim |")
    w("|-----|---|--------|------------|----------|-------|---------------|-------------|---------------|")
    for r in results:
        tag = "YES" if r['bg_dominates'] else "NO"
        w(f"| {r['P0']} | {r['N']} | {r['bg_gap']:.6e} | {r['zeros_norm']:.6e} | "
          f"{r['full_gap']:.6e} | {r['ratio']:.2f} | {tag} | "
          f"{r['bg_max']:+.6e} | {r['full_max']:+.6e} |")

    w("")
    w("## Growth Analysis")
    w("")
    if len(results) >= 2:
        w("| Transition | bg gap growth | zeros norm growth |")
        w("|------------|---------------|-------------------|")
        for i in range(1, len(results)):
            bg_g = results[i]['bg_gap'] / results[i-1]['bg_gap'] if results[i-1]['bg_gap'] > 0 else float('inf')
            z_g = results[i]['zeros_norm'] / results[i-1]['zeros_norm'] if results[i-1]['zeros_norm'] > 0 else float('inf')
            w(f"| P_0 {results[i-1]['P0']} -> {results[i]['P0']} | x{bg_g:.3f} | x{z_g:.3f} |")

    w("")
    w("## Conclusion")
    w("")
    if all_dominant and ratios_growing:
        w("**The M_bg spectral gap DOMINATES the M_zeros operator norm at all tested scales, "
          "and the dominance ratio GROWS with P_0.**")
        w("")
        w("This means the background (archimedean) contribution to the Weil matrix increasingly "
          "swamps the zeros contribution as the prime cutoff grows. The spectral structure is "
          "governed by K_bg, not by the explicit zero sum.")
    elif all_dominant:
        w("**M_bg dominates at all tested scales, but the dominance ratio does not grow monotonically.**")
        w("")
        w("Further investigation at larger scales may be needed to determine the asymptotic trend.")
    else:
        dominant_count = sum(1 for r in results if r['bg_dominates'])
        w(f"**M_bg dominates at {dominant_count}/{len(results)} tested scales.**")
        w("")
        w("The zeros contribution is significant and cannot be neglected in the spectral analysis.")

    with open(path, 'w') as f:
        f.write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
