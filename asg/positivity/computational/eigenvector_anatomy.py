#!/usr/bin/env python3
"""
Eigenvector Anatomy of Near-Zero Eigenvalues
=============================================
Analyzes eigenvectors corresponding to near-zero eigenvalues of the Weil matrix
restricted to the primitive subspace.

For each near-zero eigenvector c, we study:
  1. Which (p,m) indices carry the largest |c_i| coefficients
  2. Weighted version c_weighted_i = c_i * w_i where w_i = sqrt(log p)/p^{m/2}
  3. Concentration patterns: large primes vs high powers
  4. Spectral weight function |G_c(tau)|^2 where
     G_c(tau) = sum_i c_i * w_i * exp(i*tau*x_i),  x_i = m*log(p)
     Key question: does |G_c(tau)|^2 peak near tau ~ 2*pi?

Matrix: M_{(p,m),(q,n)} = -sqrt(log p * log q) / (p^{m/2} * q^{n/2}) * K(m log p - n log q)
Kernel: K(x) = K_bg(x) + K_zeros(x) + delta(x)
Primitive projection: P = I - v*v^T, v_{p,m} = sqrt(log p) / p^{m/2} (normalized)
"""

import numpy as np
import mpmath
import time
import sys
import os

mpmath.mp.dps = 30

SEPARATOR = "=" * 72
N_ZEROS = 200
N_NEAR_ZERO = 5
N_TOP_COEFFS = 10
N_TAU_POINTS = 1000
TAU_MAX = 20.0

CONFIGS = [
    (200, 3),
    (300, 3),
    (500, 2),
    (750, 2),
]


# ===========================================================================
# Utilities
# ===========================================================================

def sieve_primes(bound):
    is_prime = [True] * (bound + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(bound**0.5) + 1):
        if is_prime[i]:
            for j in range(i * i, bound + 1, i):
                is_prime[j] = False
    return [p for p in range(2, bound + 1) if is_prime[p]]


def compute_zeta_zeros(n_zeros):
    """Compute first n_zeros imaginary parts of non-trivial zeta zeros."""
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


# ===========================================================================
# Kernel functions
# ===========================================================================

def K_bg(x):
    """K_bg(x) = -(1/pi) Re psi(1/4 + ix/2) + (1/(2pi)) log pi"""
    arg = mpmath.mpc(0.25, float(x) / 2)
    psi = mpmath.digamma(arg)
    return float(-mpmath.re(psi) / mpmath.pi + mpmath.log(mpmath.pi) / (2 * mpmath.pi))


def K_zeros(x, zeros):
    """K_zeros(x) = (1/pi) sum_gamma cos(gamma*x)/(1/4 + gamma^2)"""
    return float(np.sum(np.cos(zeros * x) / (0.25 + zeros**2)) / np.pi)


# ===========================================================================
# Matrix construction
# ===========================================================================

def build_weil_matrix(primes, m_max, zeros):
    """
    Build Weil matrix M with entries:
    M_{(p,m),(q,n)} = -sqrt(log p * log q) / (p^{m/2} * q^{n/2}) * K(m log p - n log q)
    On diagonal: K_val = K_bg(0) + K_zeros(0) + 1
    """
    labels = [(p, a) for p in primes for a in range(1, m_max + 1)]
    N = len(labels)

    x_vals = np.array([a * np.log(p) for p, a in labels])
    weights = np.array([np.sqrt(np.log(p)) / p**(a / 2.0) for p, a in labels])

    M = np.zeros((N, N))
    t0 = time.time()

    for i in range(N):
        for j in range(i + 1):
            x = x_vals[i] - x_vals[j]
            diag = (i == j)

            k_val = K_bg(x) + K_zeros(x, zeros) + (1.0 if diag else 0.0)
            w = -weights[i] * weights[j]
            M[i, j] = w * k_val
            M[j, i] = M[i, j]

        if (i + 1) % 20 == 0 or i + 1 == N:
            el = time.time() - t0
            frac = (i + 1) / N
            eta = el / frac * (1 - frac) if frac > 0 else 0
            print(f"    row {i+1}/{N}  [{el:.1f}s, ETA ~{eta:.0f}s]")
            sys.stdout.flush()

    print(f"  {N}x{N} matrix built in {time.time()-t0:.1f}s")
    return M, labels


def primitive_projection(M, labels):
    """
    Project M onto the primitive subspace.
    Primitive subspace = orthogonal complement of v where v_{p,m} = (log p)^{1/2} / p^{m/2}.
    """
    N = M.shape[0]
    v = np.array([np.sqrt(np.log(p)) / p**(a / 2.0) for p, a in labels])
    v = v / np.linalg.norm(v)
    P = np.eye(N) - np.outer(v, v)
    Mp = P @ M @ P
    return (Mp + Mp.T) / 2


# ===========================================================================
# Eigenvector analysis
# ===========================================================================

def find_near_zero_eigenvectors(Mp, n_near=5):
    """
    Full eigendecomposition; return the n_near eigenvectors whose
    eigenvalues are closest to zero (excluding the trivial zero from projection).
    """
    eigenvalues, eigenvectors = np.linalg.eigh(Mp)

    # Remove trivial zero eigenvalue from the projection
    triv_idx = np.argmin(np.abs(eigenvalues))
    mask = np.ones(len(eigenvalues), dtype=bool)
    mask[triv_idx] = False
    eigenvalues = eigenvalues[mask]
    eigenvectors = eigenvectors[:, mask]

    # Find n_near eigenvalues closest to zero
    abs_eigs = np.abs(eigenvalues)
    near_zero_indices = np.argsort(abs_eigs)[:n_near]

    near_eigs = eigenvalues[near_zero_indices]
    near_vecs = eigenvectors[:, near_zero_indices]

    return near_eigs, near_vecs


def analyze_eigenvector(c, labels):
    """
    Analyze a single eigenvector c with given (p,m) labels.

    Returns:
        top_indices: indices of top N_TOP_COEFFS |c_i| entries
        c_weighted: c_i * w_i where w_i = sqrt(log p) / p^{m/2}
        concentration_info: dict with prime/power concentration metrics
    """
    N = len(c)
    abs_c = np.abs(c)

    # Top coefficients by |c_i|
    top_idx = np.argsort(abs_c)[::-1][:N_TOP_COEFFS]

    # Weighted version
    w = np.array([np.sqrt(np.log(p)) / p**(a / 2.0) for p, a in labels])
    c_weighted = c * w

    # Concentration analysis
    primes_arr = np.array([p for p, _ in labels])
    powers_arr = np.array([a for _, a in labels])

    # Weight of eigenvector on large primes (above median)
    median_p = np.median(primes_arr)
    large_prime_mask = primes_arr > median_p
    weight_large_primes = np.sum(c**2 * large_prime_mask)
    weight_small_primes = np.sum(c**2 * (~large_prime_mask))

    # Weight by power level
    weight_by_power = {}
    for m_val in sorted(set(powers_arr)):
        power_mask = powers_arr == m_val
        weight_by_power[int(m_val)] = float(np.sum(c**2 * power_mask))

    # Weighted version: concentration
    abs_cw = np.abs(c_weighted)
    total_cw = np.sum(abs_cw**2)
    weight_large_primes_cw = np.sum(c_weighted**2 * large_prime_mask) / total_cw if total_cw > 0 else 0.0
    weight_by_power_cw = {}
    for m_val in sorted(set(powers_arr)):
        power_mask = powers_arr == m_val
        weight_by_power_cw[int(m_val)] = float(np.sum(c_weighted**2 * power_mask) / total_cw) if total_cw > 0 else 0.0

    concentration_info = {
        'weight_large_primes': float(weight_large_primes),
        'weight_small_primes': float(weight_small_primes),
        'median_prime': float(median_p),
        'weight_by_power': weight_by_power,
        'weight_large_primes_weighted': float(weight_large_primes_cw),
        'weight_by_power_weighted': weight_by_power_cw,
    }

    return top_idx, c_weighted, concentration_info


def compute_spectral_weight(c, labels, tau_vals):
    """
    Compute |G_c(tau)|^2 where G_c(tau) = sum_i c_i * w_i * exp(i*tau*x_i).
    x_i = m * log(p), w_i = sqrt(log p) / p^{m/2}.
    """
    w = np.array([np.sqrt(np.log(p)) / p**(a / 2.0) for p, a in labels])
    x = np.array([a * np.log(p) for p, a in labels])
    cw = c * w  # c_i * w_i

    # G_c(tau) = sum_i cw_i * exp(i*tau*x_i)
    # |G_c(tau)|^2 = |sum_i cw_i * exp(i*tau*x_i)|^2
    G_sq = np.zeros(len(tau_vals))
    for k, tau in enumerate(tau_vals):
        phases = np.exp(1j * tau * x)
        G = np.sum(cw * phases)
        G_sq[k] = np.abs(G)**2

    return G_sq


# ===========================================================================
# Main driver
# ===========================================================================

def run_analysis(P0, m_max, zeros):
    """Run full eigenvector anatomy for a single (P0, m_max) configuration."""
    primes = sieve_primes(P0)
    n_primes = len(primes)
    N = n_primes * m_max

    print(f"\n{SEPARATOR}")
    print(f"  P_0 = {P0}  |  {n_primes} primes  |  m_max = {m_max}  |  {N}x{N} matrix")
    print(SEPARATOR)

    # Build matrix
    t0 = time.time()
    M, labels = build_weil_matrix(primes, m_max, zeros)
    t_build = time.time() - t0

    # Project onto primitive subspace
    Mp = primitive_projection(M, labels)

    # Full eigendecomposition
    print(f"  Computing full eigendecomposition ({N}x{N})...")
    t1 = time.time()
    near_eigs, near_vecs = find_near_zero_eigenvectors(Mp, N_NEAR_ZERO)
    t_eig = time.time() - t1
    print(f"  Eigendecomposition in {t_eig:.1f}s")

    # Tau grid for spectral weight
    tau_vals = np.linspace(0, TAU_MAX, N_TAU_POINTS)

    results = {
        'P0': P0,
        'm_max': m_max,
        'n_primes': n_primes,
        'N': N,
        't_build': t_build,
        't_eig': t_eig,
        'eigenvector_analyses': [],
    }

    for k in range(N_NEAR_ZERO):
        eig_val = near_eigs[k]
        c = near_vecs[:, k]

        print(f"\n  --- Near-zero eigenvector #{k+1}: lambda = {eig_val:+.10e} ---")

        top_idx, c_weighted, conc = analyze_eigenvector(c, labels)

        print(f"  Top {N_TOP_COEFFS} |c_i| coefficients:")
        for rank, idx in enumerate(top_idx):
            p, a = labels[idx]
            print(f"    #{rank+1}: (p={p}, m={a})  |c| = {abs(c[idx]):.6f}  "
                  f"c = {c[idx]:+.6f}  c_w = {c_weighted[idx]:+.6e}")

        print(f"  Concentration:")
        print(f"    Large primes (p > {conc['median_prime']:.0f}): "
              f"{conc['weight_large_primes']:.4f} of ||c||^2")
        print(f"    Small primes (p <= {conc['median_prime']:.0f}): "
              f"{conc['weight_small_primes']:.4f} of ||c||^2")
        for m_val, w_val in conc['weight_by_power'].items():
            print(f"    m={m_val}: {w_val:.4f} of ||c||^2")

        print(f"  Weighted concentration:")
        print(f"    Large primes: {conc['weight_large_primes_weighted']:.4f} of ||c_w||^2")
        for m_val, w_val in conc['weight_by_power_weighted'].items():
            print(f"    m={m_val}: {w_val:.4f} of ||c_w||^2")

        # Spectral weight
        print(f"  Computing |G_c(tau)|^2 ({N_TAU_POINTS} points, tau in [0, {TAU_MAX}])...")
        G_sq = compute_spectral_weight(c, labels, tau_vals)

        peak_idx = np.argmax(G_sq)
        peak_tau = tau_vals[peak_idx]
        peak_val = G_sq[peak_idx]
        val_at_2pi = G_sq[np.argmin(np.abs(tau_vals - 2 * np.pi))]

        print(f"  |G_c(tau)|^2 peak: tau = {peak_tau:.4f}, value = {peak_val:.6e}")
        print(f"  |G_c(2*pi)|^2 = {val_at_2pi:.6e}  (2*pi = {2*np.pi:.4f})")
        near_2pi = (abs(peak_tau - 2 * np.pi) < 1.0)
        print(f"  Peak near 2*pi? {'YES' if near_2pi else 'NO'} (peak at {peak_tau:.4f})")

        analysis = {
            'eigenvalue': float(eig_val),
            'top_indices': [(labels[i][0], labels[i][1], float(abs(c[i])),
                             float(c[i]), float(c_weighted[i])) for i in top_idx],
            'concentration': conc,
            'G_sq_peak_tau': float(peak_tau),
            'G_sq_peak_val': float(peak_val),
            'G_sq_at_2pi': float(val_at_2pi),
            'peak_near_2pi': bool(near_2pi),
            'G_sq': G_sq,
            'tau_vals': tau_vals,
        }
        results['eigenvector_analyses'].append(analysis)

    return results


def write_results_md(all_results, total_time):
    """Write results to markdown file."""
    out_path = os.path.join(os.path.dirname(__file__), 'eigenvector_anatomy_results.md')

    lines = []
    lines.append("# Eigenvector Anatomy of Near-Zero Eigenvalues")
    lines.append("")
    lines.append("## Configuration")
    lines.append("")
    lines.append(f"- **mpmath precision**: {mpmath.mp.dps} digits")
    lines.append(f"- **Zeta zeros used**: {N_ZEROS}")
    lines.append(f"- **Near-zero eigenvectors analyzed**: {N_NEAR_ZERO} per configuration")
    lines.append(f"- **Total computation time**: {total_time:.1f}s ({total_time/60:.1f}min)")
    lines.append("")
    lines.append("## Method")
    lines.append("")
    lines.append("For each P_0, we build the Weil matrix, project onto the primitive subspace,")
    lines.append("and compute a full eigendecomposition. We then extract the eigenvectors whose")
    lines.append("eigenvalues are closest to zero (excluding the trivial zero from projection).")
    lines.append("")
    lines.append("For each near-zero eigenvector c, we analyze:")
    lines.append("1. **Dominant indices**: which (p,m) pairs carry the largest |c_i|")
    lines.append("2. **Weighted coefficients**: c_weighted_i = c_i * sqrt(log p) / p^{m/2}")
    lines.append("3. **Concentration**: does weight concentrate on large primes or high powers?")
    lines.append("4. **Spectral weight**: |G_c(tau)|^2 where G_c(tau) = sum_i c_i * w_i * exp(i*tau*x_i)")
    lines.append("")
    lines.append("---")
    lines.append("")

    # Results for each P_0
    for res in all_results:
        P0 = res['P0']
        m_max = res['m_max']
        N = res['N']
        n_primes = res['n_primes']

        lines.append(f"## P_0 = {P0} ({n_primes} primes, m_max = {m_max}, {N}x{N} matrix)")
        lines.append("")
        lines.append(f"Build time: {res['t_build']:.1f}s | Eigendecomposition: {res['t_eig']:.1f}s")
        lines.append("")

        for k, analysis in enumerate(res['eigenvector_analyses']):
            eig = analysis['eigenvalue']
            conc = analysis['concentration']

            lines.append(f"### Near-zero eigenvector #{k+1}: lambda = {eig:+.10e}")
            lines.append("")

            # Top indices table
            lines.append("**Top 10 |c_i| coefficients:**")
            lines.append("")
            lines.append("| Rank | p | m | |c_i| | c_i | c_weighted_i |")
            lines.append("|------|---|---|-------|-----|--------------|")
            for rank, (p, m, abs_c, c_val, cw_val) in enumerate(analysis['top_indices']):
                lines.append(f"| {rank+1} | {p} | {m} | {abs_c:.6f} | {c_val:+.6f} | {cw_val:+.6e} |")
            lines.append("")

            # Concentration
            lines.append("**Concentration pattern:**")
            lines.append("")
            lines.append(f"- Large primes (p > {conc['median_prime']:.0f}): "
                         f"{conc['weight_large_primes']:.4f} of ||c||^2")
            lines.append(f"- Small primes (p <= {conc['median_prime']:.0f}): "
                         f"{conc['weight_small_primes']:.4f} of ||c||^2")
            for m_val, w_val in conc['weight_by_power'].items():
                lines.append(f"- m={m_val}: {w_val:.4f} of ||c||^2")
            lines.append("")
            lines.append("Weighted concentration (c_weighted):")
            lines.append("")
            lines.append(f"- Large primes: {conc['weight_large_primes_weighted']:.4f} of ||c_w||^2")
            for m_val, w_val in conc['weight_by_power_weighted'].items():
                lines.append(f"- m={m_val}: {w_val:.4f} of ||c_w||^2")
            lines.append("")

            # Spectral weight
            lines.append("**Spectral weight |G_c(tau)|^2:**")
            lines.append("")
            lines.append(f"- Peak location: tau = {analysis['G_sq_peak_tau']:.4f}")
            lines.append(f"- Peak value: {analysis['G_sq_peak_val']:.6e}")
            lines.append(f"- Value at tau = 2*pi (6.2832): {analysis['G_sq_at_2pi']:.6e}")
            lines.append(f"- Peak near 2*pi? **{'YES' if analysis['peak_near_2pi'] else 'NO'}**")
            lines.append("")

        lines.append("---")
        lines.append("")

    # Summary tables
    lines.append("## Summary Tables")
    lines.append("")

    # Table 1: eigenvalues and concentration
    lines.append("### Eigenvalues and Dominant Indices")
    lines.append("")
    lines.append("| P_0 | # | lambda | Top (p,m) | Large-prime frac | Dominant m |")
    lines.append("|-----|---|--------|-----------|-----------------|------------|")
    for res in all_results:
        P0 = res['P0']
        for k, analysis in enumerate(res['eigenvector_analyses']):
            eig = analysis['eigenvalue']
            top_pm = analysis['top_indices'][0]  # top 1
            large_frac = analysis['concentration']['weight_large_primes']
            # dominant power
            wbp = analysis['concentration']['weight_by_power']
            dom_m = max(wbp, key=wbp.get)
            dom_m_frac = wbp[dom_m]
            lines.append(f"| {P0} | {k+1} | {eig:+.4e} | ({top_pm[0]},{top_pm[1]}) "
                         f"| {large_frac:.3f} | m={dom_m} ({dom_m_frac:.3f}) |")
    lines.append("")

    # Table 2: G_c(tau) peak locations
    lines.append("### Spectral Weight G_c(tau) Peak Locations")
    lines.append("")
    lines.append("| P_0 | # | lambda | Peak tau | Peak value | |G(2pi)|^2 | Near 2pi? |")
    lines.append("|-----|---|--------|----------|------------|-----------|-----------|")
    for res in all_results:
        P0 = res['P0']
        for k, analysis in enumerate(res['eigenvector_analyses']):
            eig = analysis['eigenvalue']
            pt = analysis['G_sq_peak_tau']
            pv = analysis['G_sq_peak_val']
            g2pi = analysis['G_sq_at_2pi']
            near = "YES" if analysis['peak_near_2pi'] else "NO"
            lines.append(f"| {P0} | {k+1} | {eig:+.4e} | {pt:.4f} | {pv:.4e} "
                         f"| {g2pi:.4e} | {near} |")
    lines.append("")

    # Observations
    lines.append("## Observations")
    lines.append("")

    # Collect peak stats
    n_near_2pi = 0
    n_total = 0
    all_peaks = []
    for res in all_results:
        for analysis in res['eigenvector_analyses']:
            n_total += 1
            if analysis['peak_near_2pi']:
                n_near_2pi += 1
            all_peaks.append(analysis['G_sq_peak_tau'])

    lines.append(f"1. **G_c(tau) peak near 2*pi**: {n_near_2pi}/{n_total} eigenvectors have "
                 f"|G_c(tau)|^2 peaking within 1.0 of tau = 2*pi = 6.2832")
    lines.append(f"2. **Peak tau distribution**: min = {min(all_peaks):.4f}, "
                 f"max = {max(all_peaks):.4f}, "
                 f"mean = {np.mean(all_peaks):.4f}")
    lines.append("")

    # Concentration patterns
    all_large_frac = []
    all_m1_frac = []
    for res in all_results:
        for analysis in res['eigenvector_analyses']:
            all_large_frac.append(analysis['concentration']['weight_large_primes'])
            wbp = analysis['concentration']['weight_by_power']
            if 1 in wbp:
                all_m1_frac.append(wbp[1])

    if all_large_frac:
        lines.append(f"3. **Large-prime concentration**: mean = {np.mean(all_large_frac):.4f}, "
                     f"range = [{min(all_large_frac):.4f}, {max(all_large_frac):.4f}]")
    if all_m1_frac:
        lines.append(f"4. **m=1 concentration**: mean = {np.mean(all_m1_frac):.4f}, "
                     f"range = [{min(all_m1_frac):.4f}, {max(all_m1_frac):.4f}]")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append(f"*Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}*")

    with open(out_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    print(f"\n  Results written to {out_path}")


def main():
    t_start = time.time()
    print(SEPARATOR)
    print("  EIGENVECTOR ANATOMY OF NEAR-ZERO EIGENVALUES")
    print(SEPARATOR)
    print(f"  mpmath precision: {mpmath.mp.dps} digits")
    print(f"  Zeta zeros: {N_ZEROS}")
    print(f"  Near-zero eigenvectors per config: {N_NEAR_ZERO}")
    print(f"  Configs: {CONFIGS}")
    print()

    # Compute zeta zeros
    zeros = compute_zeta_zeros(N_ZEROS)
    print()

    all_results = []
    for P0, m_max in CONFIGS:
        res = run_analysis(P0, m_max, zeros)
        all_results.append(res)
        sys.stdout.flush()

    # Summary
    print(f"\n{SEPARATOR}")
    print("  SUMMARY")
    print(SEPARATOR)
    for res in all_results:
        P0 = res['P0']
        print(f"\n  P_0 = {P0} ({res['N']}x{res['N']}):")
        for k, analysis in enumerate(res['eigenvector_analyses']):
            eig = analysis['eigenvalue']
            peak = analysis['G_sq_peak_tau']
            near = "YES" if analysis['peak_near_2pi'] else "NO"
            top = analysis['top_indices'][0]
            print(f"    #{k+1}: lambda={eig:+.4e}  peak_tau={peak:.3f}  "
                  f"near_2pi={near}  top=({top[0]},{top[1]})")

    total_time = time.time() - t_start
    print(f"\n  Total time: {total_time:.1f}s ({total_time/60:.1f}min)")

    # Write markdown results
    write_results_md(all_results, total_time)

    print(SEPARATOR)
    return all_results


if __name__ == '__main__':
    main()
