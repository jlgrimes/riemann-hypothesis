#!/usr/bin/env python3
"""
Random vs Real Primes: Is APT a Prime-Specific Property?
=========================================================
The most important investigation. Tests whether the Arithmetic Positivity
Theorem (all primitive eigenvalues of the Weil matrix are <= 0) is specific
to actual primes, or whether it holds for arbitrary sets of numbers with
prime-like density.

If random numbers also satisfy APT -> the pattern is in the kernel, not primes.
If only real primes satisfy APT -> something is structurally special about primes.

Tests:
  1. Real primes (P_0=200, m_max=3)
  2. 50 random "pseudo-prime" sequences (same count, density ~ 1/log(n))
  3. Shifted primes {p+1}
  4. Composites only
"""

import numpy as np
from scipy.linalg import eigvalsh
import mpmath
import time
import sys
import os

mpmath.mp.dps = 30
SEPARATOR = "=" * 72


# ═══════════════════════════════════════════════════════════════════
#  Infrastructure
# ═══════════════════════════════════════════════════════════════════

def sieve_primes(bound):
    is_prime = [True] * (int(bound) + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(bound**0.5) + 1):
        if is_prime[i]:
            for j in range(i * i, int(bound) + 1, i):
                is_prime[j] = False
    return [p for p in range(2, int(bound) + 1) if is_prime[p]]


def compute_zeta_zeros(n_zeros):
    print(f"  Computing {n_zeros} zeta zeros (dps={mpmath.mp.dps})...")
    t0 = time.time()
    zeros = []
    for k in range(1, n_zeros + 1):
        g = float(mpmath.zetazero(k).imag)
        zeros.append(g)
        if k % 50 == 0:
            print(f"    {k}/{n_zeros} [{time.time()-t0:.1f}s]")
    zeros = np.array(zeros)
    print(f"  {n_zeros} zeros in {time.time()-t0:.1f}s")
    return zeros


def K_bg(x):
    arg = mpmath.mpc(0.25, float(x) / 2)
    psi = mpmath.digamma(arg)
    return float(-mpmath.re(psi) / mpmath.pi + mpmath.log(mpmath.pi) / (2 * mpmath.pi))


def K_zeros_val(x, zeros):
    return float(np.sum(np.cos(zeros * x) / (0.25 + zeros**2)) / np.pi)


def primitive_projection(M, labels):
    N = M.shape[0]
    v = np.array([np.sqrt(np.log(n)) / n**(a / 2.0) for n, a in labels])
    v = v / np.linalg.norm(v)
    P = np.eye(N) - np.outer(v, v)
    Mp = P @ M @ P
    return (Mp + Mp.T) / 2


def analyze_eigenvalues(Mp):
    eigs = eigvalsh(Mp)
    idx = np.argmin(np.abs(eigs))
    prim = np.delete(eigs, idx)
    return np.sort(prim)


# ═══════════════════════════════════════════════════════════════════
#  K_bg interpolation grid (for fast random trials)
# ═══════════════════════════════════════════════════════════════════

def precompute_K_bg_grid(x_min=-20.0, x_max=20.0, n_points=40001):
    """Precompute K_bg on a dense grid for fast interpolation."""
    print(f"  Precomputing K_bg on grid [{x_min}, {x_max}], {n_points} points...")
    t0 = time.time()
    grid_x = np.linspace(x_min, x_max, n_points)
    grid_y = np.empty(n_points)
    for i, x in enumerate(grid_x):
        grid_y[i] = K_bg(x)
        if (i + 1) % 10000 == 0:
            print(f"    {i+1}/{n_points} [{time.time()-t0:.1f}s]")
    print(f"  K_bg grid computed in {time.time()-t0:.1f}s")
    return grid_x, grid_y


# ═══════════════════════════════════════════════════════════════════
#  Matrix builders
# ═══════════════════════════════════════════════════════════════════

def build_weil_matrix_exact(numbers, m_max, zeros, verbose=True):
    """Build Weil matrix using exact K_bg (mpmath). For real primes."""
    labels = [(n, a) for n in numbers for a in range(1, m_max + 1)]
    N = len(labels)
    x_vals = np.array([a * np.log(n) for n, a in labels])
    weights = np.array([np.sqrt(np.log(n)) / n**(a / 2.0) for n, a in labels])

    M = np.zeros((N, N))
    t0 = time.time()
    for i in range(N):
        for j in range(i + 1):
            x = x_vals[i] - x_vals[j]
            diag = (i == j)
            k = K_bg(x) + K_zeros_val(x, zeros) + (1.0 if diag else 0.0)
            w = -weights[i] * weights[j]
            M[i, j] = w * k
            M[j, i] = M[i, j]
        if verbose and ((i + 1) % 30 == 0 or i + 1 == N):
            el = time.time() - t0
            frac = (i + 1) / N
            eta = el / frac * (1 - frac) if frac > 0 else 0
            print(f"    row {i+1}/{N}  [{el:.1f}s, ETA ~{eta:.0f}s]")
            sys.stdout.flush()
    if verbose:
        print(f"  {N}x{N} matrix built in {time.time()-t0:.1f}s")
    return M, labels


def build_weil_matrix_fast(numbers, m_max, zeros, grid_x, grid_y):
    """Build Weil matrix using interpolated K_bg. For random trials."""
    labels = [(n, a) for n in numbers for a in range(1, m_max + 1)]
    N = len(labels)
    x_vals = np.array([a * np.log(n) for n, a in labels])
    weights = np.array([np.sqrt(np.log(n)) / n**(a / 2.0) for n, a in labels])

    # Pairwise argument differences
    X = x_vals[:, None] - x_vals[None, :]

    # K_bg via interpolation
    Kb = np.interp(X.ravel(), grid_x, grid_y).reshape(N, N)

    # K_zeros vectorized: Kz[i,j] = (1/pi) sum_g cos(g * X[i,j]) / (1/4 + g^2)
    X_flat = X.ravel()
    cos_mat = np.cos(np.outer(X_flat, zeros))  # (N^2) x n_zeros
    w_zeros = 1.0 / (0.25 + zeros**2) / np.pi
    Kz = (cos_mat @ w_zeros).reshape(N, N)

    # Full kernel + delta
    K = Kb + Kz + np.eye(N)

    # Weight matrix
    W = -np.outer(weights, weights)
    M = W * K
    M = (M + M.T) / 2

    return M, labels


# ═══════════════════════════════════════════════════════════════════
#  Number generators
# ═══════════════════════════════════════════════════════════════════

def generate_pseudo_primes(count, max_val, rng):
    """Sample count numbers from [2, max_val] with density ~ 1/log(n)."""
    candidates = np.arange(2, max_val + 1)
    w = 1.0 / np.log(candidates.astype(float))
    w /= w.sum()
    return sorted(rng.choice(candidates, size=count, replace=False, p=w))


def generate_composites(count, max_val, primes_set, rng):
    """Sample count composite numbers from [4, max_val]."""
    composites = [n for n in range(4, max_val + 1) if n not in primes_set]
    return sorted(rng.choice(composites, size=min(count, len(composites)), replace=False))


# ═══════════════════════════════════════════════════════════════════
#  Main investigation
# ═══════════════════════════════════════════════════════════════════

def main():
    t_start = time.time()
    print(SEPARATOR)
    print("  RANDOM VS REAL PRIMES: IS APT A PRIME-SPECIFIC PROPERTY?")
    print(SEPARATOR)
    print()

    N_ZEROS = 200
    zeros = compute_zeta_zeros(N_ZEROS)
    print()

    # ── Precompute K_bg grid ──
    grid_x, grid_y = precompute_K_bg_grid()
    print()

    # ── Configuration ──
    P0 = 200
    m_max = 3
    real_primes = sieve_primes(P0)
    n_primes = len(real_primes)
    N_mat = n_primes * m_max
    n_trials = 50

    print(f"  P_0 = {P0}, {n_primes} primes, m_max = {m_max}, matrix {N_mat}x{N_mat}")
    print(f"  Random trials: {n_trials}")
    print()

    # ══════════════════════════════════════════════════════════
    #  Test 1: Real primes
    # ══════════════════════════════════════════════════════════
    print(f"{SEPARATOR}")
    print("  TEST 1: REAL PRIMES")
    print(f"{SEPARATOR}")

    M_real, labels_real = build_weil_matrix_exact(real_primes, m_max, zeros)
    Mp_real = primitive_projection(M_real, labels_real)
    eigs_real = analyze_eigenvalues(Mp_real)
    max_real = float(eigs_real[-1])
    min_real = float(eigs_real[0])
    apt_real = bool(np.all(eigs_real < 1e-12))

    print(f"\n  Real primes: max_prim = {max_real:+.8e}, min_prim = {min_real:+.8e}")
    print(f"  APT: {'YES' if apt_real else 'NO'}")
    print()

    # ══════════════════════════════════════════════════════════
    #  Test 2: Random pseudo-prime trials
    # ══════════════════════════════════════════════════════════
    print(f"{SEPARATOR}")
    print("  TEST 2: RANDOM PSEUDO-PRIME TRIALS")
    print(f"{SEPARATOR}")
    print(f"  Generating {n_trials} random sets of {n_primes} numbers (density ~ 1/log(n))...")
    print()

    rng = np.random.default_rng(42)
    random_max_eigs = []
    random_min_eigs = []
    random_apt = []

    t0 = time.time()
    for trial in range(n_trials):
        pseudo = generate_pseudo_primes(n_primes, P0, rng)
        M_rand, labels_rand = build_weil_matrix_fast(pseudo, m_max, zeros, grid_x, grid_y)
        Mp_rand = primitive_projection(M_rand, labels_rand)
        eigs_rand = analyze_eigenvalues(Mp_rand)
        max_rand = float(eigs_rand[-1])
        min_rand = float(eigs_rand[0])
        apt = bool(np.all(eigs_rand < 1e-12))

        random_max_eigs.append(max_rand)
        random_min_eigs.append(min_rand)
        random_apt.append(apt)

        if (trial + 1) % 10 == 0:
            elapsed = time.time() - t0
            print(f"    Trial {trial+1}/{n_trials}: max={max_rand:+.4e}, APT={'Y' if apt else 'N'} "
                  f"[{elapsed:.1f}s]")

    print(f"  {n_trials} trials completed in {time.time()-t0:.1f}s")

    random_max_arr = np.array(random_max_eigs)
    mean_random = float(np.mean(random_max_arr))
    std_random = float(np.std(random_max_arr))
    min_random_max = float(np.min(random_max_arr))
    max_random_max = float(np.max(random_max_arr))
    pass_rate = sum(random_apt) / n_trials

    if std_random > 1e-15:
        z_score = (max_real - mean_random) / std_random
    else:
        z_score = float('-inf') if max_real < mean_random else float('inf')

    print(f"\n  Random trial statistics (max primitive eigenvalue):")
    print(f"    Mean:   {mean_random:+.8e}")
    print(f"    Std:    {std_random:.8e}")
    print(f"    Min:    {min_random_max:+.8e}")
    print(f"    Max:    {max_random_max:+.8e}")
    print(f"    APT pass rate: {sum(random_apt)}/{n_trials} = {pass_rate:.0%}")
    print(f"\n  Real primes max:  {max_real:+.8e}")
    print(f"  Z-score (real vs random): {z_score:+.4f}")
    print()

    # ══════════════════════════════════════════════════════════
    #  Test 3: Shifted primes {p+1}
    # ══════════════════════════════════════════════════════════
    print(f"{SEPARATOR}")
    print("  TEST 3: SHIFTED PRIMES {p+1}")
    print(f"{SEPARATOR}")

    shifted = [p + 1 for p in real_primes]
    M_shift, labels_shift = build_weil_matrix_fast(shifted, m_max, zeros, grid_x, grid_y)
    Mp_shift = primitive_projection(M_shift, labels_shift)
    eigs_shift = analyze_eigenvalues(Mp_shift)
    max_shift = float(eigs_shift[-1])
    min_shift = float(eigs_shift[0])
    apt_shift = bool(np.all(eigs_shift < 1e-12))

    print(f"  Shifted primes: max_prim = {max_shift:+.8e}, min = {min_shift:+.8e}")
    print(f"  APT: {'YES' if apt_shift else 'NO'}")
    print()

    # ══════════════════════════════════════════════════════════
    #  Test 4: Composites only
    # ══════════════════════════════════════════════════════════
    print(f"{SEPARATOR}")
    print("  TEST 4: COMPOSITES ONLY")
    print(f"{SEPARATOR}")

    primes_set = set(real_primes)
    composites = generate_composites(n_primes, P0, primes_set, rng)
    M_comp, labels_comp = build_weil_matrix_fast(composites, m_max, zeros, grid_x, grid_y)
    Mp_comp = primitive_projection(M_comp, labels_comp)
    eigs_comp = analyze_eigenvalues(Mp_comp)
    max_comp = float(eigs_comp[-1])
    min_comp = float(eigs_comp[0])
    apt_comp = bool(np.all(eigs_comp < 1e-12))

    print(f"  Composites: max_prim = {max_comp:+.8e}, min = {min_comp:+.8e}")
    print(f"  APT: {'YES' if apt_comp else 'NO'}")
    print()

    # ══════════════════════════════════════════════════════════
    #  Summary
    # ══════════════════════════════════════════════════════════
    print(f"{SEPARATOR}")
    print("  SUMMARY")
    print(f"{SEPARATOR}")
    print()
    print(f"  {'Test':>20s} {'max_prim':>14s} {'APT':>5s} {'Z-score':>10s}")
    print("  " + "-" * 55)
    print(f"  {'Real primes':>20s} {max_real:+14.6e} {'YES' if apt_real else 'NO':>5s} {'---':>10s}")
    print(f"  {'Random (mean)':>20s} {mean_random:+14.6e} {f'{pass_rate:.0%}':>5s} {z_score:+10.4f}")
    print(f"  {'Shifted {p+1}':>20s} {max_shift:+14.6e} {'YES' if apt_shift else 'NO':>5s} {'---':>10s}")
    print(f"  {'Composites':>20s} {max_comp:+14.6e} {'YES' if apt_comp else 'NO':>5s} {'---':>10s}")
    print()

    if pass_rate > 0.9:
        print("  CONCLUSION: Random numbers ALSO satisfy APT at high rate.")
        print("  => APT is likely a property of the KERNEL, not specific to primes.")
    elif pass_rate < 0.1:
        print("  CONCLUSION: Random numbers RARELY satisfy APT.")
        print("  => APT appears to be PRIME-SPECIFIC.")
    else:
        print(f"  CONCLUSION: Random APT pass rate = {pass_rate:.0%} — intermediate.")
        print("  => Primes may have partial structural advantage.")

    total_time = time.time() - t_start
    print(f"\n  Total time: {total_time:.1f}s ({total_time/60:.1f} min)")
    print(SEPARATOR)

    # ══════════════════════════════════════════════════════════
    #  Write results markdown
    # ══════════════════════════════════════════════════════════
    write_results_md(
        real_primes, max_real, min_real, apt_real,
        random_max_eigs, random_min_eigs, random_apt,
        mean_random, std_random, z_score, pass_rate,
        shifted, max_shift, min_shift, apt_shift,
        composites, max_comp, min_comp, apt_comp,
        P0, m_max, n_primes, n_trials, total_time,
    )


def write_results_md(
    real_primes, max_real, min_real, apt_real,
    random_max_eigs, random_min_eigs, random_apt,
    mean_random, std_random, z_score, pass_rate,
    shifted, max_shift, min_shift, apt_shift,
    composites, max_comp, min_comp, apt_comp,
    P0, m_max, n_primes, n_trials, total_time,
):
    out_path = os.path.join(os.path.dirname(__file__), 'random_vs_real_primes_results.md')

    L = []
    w = L.append

    w("# Random vs Real Primes: Is APT a Prime-Specific Property?")
    w("")
    w("**The key question**: Does the Arithmetic Positivity Theorem (all primitive")
    w("eigenvalues of the Weil matrix <= 0) hold because of something special about")
    w("primes, or is it a property of the kernel K(x) = K_bg(x) + K_zeros(x)?")
    w("")
    w("## Configuration")
    w("")
    w(f"- P_0 = {P0}, m_max = {m_max}")
    w(f"- {n_primes} real primes, matrix size {n_primes * m_max}x{n_primes * m_max}")
    w(f"- {n_trials} random pseudo-prime trials (density ~ 1/log(n))")
    w(f"- 200 zeta zeros, mpmath dps = 30")
    w(f"- Total computation time: {total_time:.1f}s")
    w("")
    w("---")
    w("")

    w("## Results Summary")
    w("")
    w("| Test | max prim eig | min prim eig | APT? | Z-score |")
    w("|------|-------------|-------------|------|---------|")
    w(f"| **Real primes** | {max_real:+.6e} | {min_real:+.6e} | "
      f"{'**YES**' if apt_real else 'NO'} | — |")
    w(f"| Random (mean +/- std) | {mean_random:+.6e} +/- {std_random:.2e} | — | "
      f"{sum(random_apt)}/{n_trials} ({pass_rate:.0%}) | {z_score:+.2f} |")
    w(f"| Shifted primes (p+1) | {max_shift:+.6e} | {min_shift:+.6e} | "
      f"{'YES' if apt_shift else 'NO'} | — |")
    w(f"| Composites only | {max_comp:+.6e} | {min_comp:+.6e} | "
      f"{'YES' if apt_comp else 'NO'} | — |")
    w("")

    w("## Random Trial Distribution")
    w("")
    w("| Statistic | max prim eigenvalue |")
    w("|-----------|-------------------|")
    w(f"| Mean | {mean_random:+.8e} |")
    w(f"| Std | {std_random:.8e} |")
    w(f"| Min | {min(random_max_eigs):+.8e} |")
    w(f"| Max | {max(random_max_eigs):+.8e} |")
    w(f"| Median | {float(np.median(random_max_eigs)):+.8e} |")
    w(f"| APT pass rate | {sum(random_apt)}/{n_trials} ({pass_rate:.0%}) |")
    w("")
    w(f"**Real primes max eigenvalue**: {max_real:+.8e}")
    w(f"**Z-score**: {z_score:+.4f}")
    w("")

    w("### Individual Random Trial Results")
    w("")
    w("| Trial | max prim eig | APT? |")
    w("|-------|-------------|------|")
    for i in range(n_trials):
        apt_str = "YES" if random_apt[i] else "NO"
        w(f"| {i+1} | {random_max_eigs[i]:+.6e} | {apt_str} |")
    w("")

    w("---")
    w("")
    w("## Analysis")
    w("")

    if pass_rate > 0.9:
        w("### Finding: APT is a KERNEL property, not prime-specific")
        w("")
        w(f"Random numbers satisfy APT at a {pass_rate:.0%} rate, comparable to real primes.")
        w("This strongly suggests that the negative-definiteness of the Weil matrix on the")
        w("primitive subspace arises from the structure of the kernel K(x) = K_bg(x) + K_zeros(x),")
        w("not from any special property of the prime numbers themselves.")
        w("")
        w("**Implication**: The Weil matrix encodes the Riemann Hypothesis through its kernel,")
        w("which is built from the zeta zeros. The primes serve as a 'basis' for sampling")
        w("the kernel, but other bases work too.")
    elif pass_rate < 0.1:
        w("### Finding: APT is PRIME-SPECIFIC")
        w("")
        w(f"Random numbers satisfy APT at only a {pass_rate:.0%} rate.")
        w(f"Real primes achieve z-score = {z_score:+.2f} relative to the random distribution.")
        w("")
        w("This demonstrates that APT depends critically on the input being actual primes.")
        w("The explicit formula connects primes to zeros in a way that random numbers cannot")
        w("replicate. The Weil matrix 'knows' about primes.")
    else:
        w(f"### Finding: Intermediate result (pass rate = {pass_rate:.0%})")
        w("")
        w("The pass rate suggests partial prime-specificity. The kernel contributes")
        w("significant negative-definiteness, but primes may provide additional structure.")

    w("")

    if not apt_shift:
        w("**Shifted primes {p+1}**: Shifting primes by 1 BREAKS APT, showing that the")
        w("exact positions of primes matter, not just their approximate distribution.")
    else:
        w("**Shifted primes {p+1}**: Shifting primes by 1 still satisfies APT, suggesting")
        w("the property is robust to small perturbations of the input set.")
    w("")

    if not apt_comp:
        w("**Composites**: Using composites BREAKS APT, confirming that primality is essential.")
    else:
        w("**Composites**: Even composites satisfy APT, further suggesting the kernel dominates.")
    w("")

    w("---")
    w("")
    w("## Conclusion")
    w("")
    if pass_rate > 0.5:
        w("The Weil matrix's spectral negativity on the primitive subspace is primarily a")
        w("consequence of the kernel structure (K_bg + K_zeros) rather than a unique property")
        w("of prime numbers. The kernel, which encodes the zeta zeros via the explicit formula,")
        w("produces negative primitive eigenvalues regardless of which numbers are used as inputs.")
    else:
        w("The Weil matrix's spectral negativity on the primitive subspace is NOT merely a")
        w("consequence of the kernel. Real primes produce systematically more negative")
        w("eigenvalues than random numbers, revealing a structural signature of primes")
        w("encoded in the Weil matrix.")
    w("")
    w(f"*Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}*")

    with open(out_path, 'w') as f:
        f.write('\n'.join(L) + '\n')
    print(f"\n  Results written to {out_path}")


if __name__ == '__main__':
    main()
