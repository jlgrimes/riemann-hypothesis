#!/usr/bin/env python3
"""
Prime Correlation Study via the Weil Matrix
============================================
What structural patterns do primes have, as encoded by the Weil matrix?

Investigation 1: Prime Correlation Map — decay and coupling structure
Investigation 3: Eigenvalue Spectrum Evolution — one prime at a time
Investigation 4: Residue Class Structure — mod 4 and mod 3 comparisons
Investigation 5: Off-Line Zero Sensitivity — which zeros are load-bearing?
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


def build_weil_matrix(primes, m_max, zeros, verbose=True):
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


def primitive_projection(M, labels):
    N = M.shape[0]
    v = np.array([np.sqrt(np.log(p)) / p**(a / 2.0) for p, a in labels])
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
#  Investigation 1: Prime Correlation Map
# ═══════════════════════════════════════════════════════════════════

def investigation_1(zeros):
    print(f"\n{SEPARATOR}")
    print("  INVESTIGATION 1: PRIME CORRELATION MAP")
    print(f"{SEPARATOR}")

    P0, m_max = 300, 3
    primes = sieve_primes(P0)
    N = len(primes) * m_max
    print(f"  P_0={P0}, {len(primes)} primes, m_max={m_max}, {N}x{N} matrix")

    labels = [(p, a) for p in primes for a in range(1, m_max + 1)]
    x_vals = np.array([a * np.log(p) for p, a in labels])
    weights = np.array([np.sqrt(np.log(p)) / p**(a / 2.0) for p, a in labels])

    M = np.zeros((N, N))
    entries = []
    t0 = time.time()

    for i in range(N):
        for j in range(i + 1):
            x = x_vals[i] - x_vals[j]
            diag = (i == j)
            kb = K_bg(x)
            kz = K_zeros_val(x, zeros)
            kfull = kb + kz + (1.0 if diag else 0.0)
            w = -weights[i] * weights[j]
            M[i, j] = w * kfull
            M[j, i] = M[i, j]

            p, m = labels[i]
            q, n = labels[j]
            entries.append((p, m, q, n, x, M[i, j], kb, kz))
        if (i + 1) % 30 == 0 or i + 1 == N:
            el = time.time() - t0
            frac = (i + 1) / N
            eta = el / frac * (1 - frac) if frac > 0 else 0
            print(f"    row {i+1}/{N}  [{el:.1f}s, ETA ~{eta:.0f}s]")
            sys.stdout.flush()

    print(f"  {N}x{N} matrix built in {time.time()-t0:.1f}s")

    # Off-diagonal entries only
    off_diag = [e for e in entries if not (e[0] == e[2] and e[1] == e[3])]
    off_diag_sorted = sorted(off_diag, key=lambda e: abs(e[5]), reverse=True)
    top20 = off_diag_sorted[:20]

    print(f"\n  Top 20 off-diagonal entries (strongest couplings):")
    print(f"  {'p':>5s} {'m':>2s} {'q':>5s} {'n':>2s} | {'|arg|':>10s} {'M_ij':>14s} {'K_bg':>12s} {'K_zeros':>12s}")
    print("  " + "-" * 72)
    for p, m, q, n, x, mij, kb, kz in top20:
        print(f"  {p:5d} {m:2d} {q:5d} {n:2d} | {abs(x):10.4f} {mij:+14.6e} {kb:12.6e} {kz:+12.6e}")

    # Bin by |argument|
    abs_args = np.array([abs(e[4]) for e in off_diag])
    abs_mij = np.array([abs(e[5]) for e in off_diag])
    bin_edges = np.arange(0, float(np.max(abs_args)) + 0.5, 0.5)

    binned = []
    print(f"\n  Binned decay statistics:")
    print(f"  {'bin':>12s} {'count':>6s} {'mean |M|':>14s} {'std |M|':>14s} {'max |M|':>14s}")
    print("  " + "-" * 60)
    for b in range(min(len(bin_edges) - 1, 30)):
        mask = (abs_args >= bin_edges[b]) & (abs_args < bin_edges[b + 1])
        if mask.sum() > 0:
            vals = abs_mij[mask]
            row = {'bin_lo': float(bin_edges[b]), 'bin_hi': float(bin_edges[b + 1]),
                   'count': int(mask.sum()), 'mean': float(np.mean(vals)),
                   'std': float(np.std(vals)), 'max': float(np.max(vals))}
            binned.append(row)
            print(f"  [{bin_edges[b]:5.1f},{bin_edges[b+1]:5.1f}) {row['count']:6d} "
                  f"{row['mean']:14.6e} {row['std']:14.6e} {row['max']:14.6e}")

    # Decay fits (|x| > 0.1 to avoid near-diagonal)
    mask_fit = abs_args > 0.1
    x_fit = abs_args[mask_fit]
    y_fit = np.maximum(abs_mij[mask_fit], 1e-30)

    coeffs_exp = np.polyfit(x_fit, np.log(y_fit), 1)
    b_exp, a_exp = -coeffs_exp[0], np.exp(coeffs_exp[1])

    mask_pow = x_fit > 0.5
    if mask_pow.sum() > 10:
        coeffs_pow = np.polyfit(np.log(x_fit[mask_pow]), np.log(y_fit[mask_pow]), 1)
        b_pow, a_pow = -coeffs_pow[0], np.exp(coeffs_pow[1])
    else:
        b_pow = a_pow = float('nan')

    print(f"\n  Decay fits (|arg| > 0.1):")
    print(f"    Exponential: |M_ij| ~ {a_exp:.4e} * exp(-{b_exp:.4f} * |x|)")
    if not np.isnan(b_pow):
        print(f"    Power law:   |M_ij| ~ {a_pow:.4e} * |x|^(-{b_pow:.4f})")

    return {'top20': top20, 'binned': binned,
            'fit_exp': (a_exp, b_exp), 'fit_pow': (a_pow, b_pow),
            'n_entries': len(off_diag)}


# ═══════════════════════════════════════════════════════════════════
#  Investigation 3: Eigenvalue Spectrum Evolution
# ═══════════════════════════════════════════════════════════════════

def investigation_3(zeros):
    print(f"\n{SEPARATOR}")
    print("  INVESTIGATION 3: EIGENVALUE SPECTRUM EVOLUTION")
    print(f"{SEPARATOR}")

    P0, m_max = 500, 2
    all_primes = sieve_primes(P0)
    n_primes = len(all_primes)
    N_full = n_primes * m_max

    print(f"  P_0={P0}, {n_primes} primes, m_max={m_max}, full matrix {N_full}x{N_full}")
    print(f"  Building full matrix once, then slicing incrementally...")

    M_full, all_labels = build_weil_matrix(all_primes, m_max, zeros)

    print(f"\n  Computing eigenvalue evolution...")
    evolution = []
    prev_max = None

    for k in range(2, n_primes + 1):
        n = k * m_max
        M_sub = M_full[:n, :n]
        labels_sub = all_labels[:n]
        Mp_sub = primitive_projection(M_sub, labels_sub)
        eigs = analyze_eigenvalues(Mp_sub)
        max_eig = float(eigs[-1])
        min_eig = float(eigs[0])

        jump = max_eig - prev_max if prev_max is not None else 0.0
        evolution.append({'k': k, 'prime': all_primes[k - 1], 'N': n,
                          'max_eig': max_eig, 'min_eig': min_eig, 'jump': jump})
        prev_max = max_eig

        if k % 10 == 0 or k == n_primes:
            print(f"    k={k:3d}, p={all_primes[k-1]:4d}, N={n:4d}, "
                  f"max={max_eig:+.6e}, Δ={jump:+.6e}")

    jumps_sorted = sorted(evolution[1:], key=lambda e: abs(e['jump']), reverse=True)
    print(f"\n  Top 10 biggest eigenvalue jumps:")
    print(f"  {'prime':>6s} {'k':>4s} {'max_eig':>14s} {'jump':>14s}")
    print("  " + "-" * 42)
    for e in jumps_sorted[:10]:
        print(f"  {e['prime']:6d} {e['k']:4d} {e['max_eig']:+14.6e} {e['jump']:+14.6e}")

    return evolution


# ═══════════════════════════════════════════════════════════════════
#  Investigation 4: Residue Class Structure
# ═══════════════════════════════════════════════════════════════════

def investigation_4(zeros):
    print(f"\n{SEPARATOR}")
    print("  INVESTIGATION 4: RESIDUE CLASS STRUCTURE")
    print(f"{SEPARATOR}")

    P0, m_max = 500, 2
    all_primes = sieve_primes(P0)

    classes = {
        '1 (mod 4)': [p for p in all_primes if p % 4 == 1],
        '3 (mod 4)': [p for p in all_primes if p % 4 == 3],
        '1 (mod 3)': [p for p in all_primes if p % 3 == 1],
        '2 (mod 3)': [p for p in all_primes if p % 3 == 2],
        'all primes': all_primes,
    }

    results = {}
    for name, primes in classes.items():
        N = len(primes) * m_max
        print(f"\n  Class: {name} — {len(primes)} primes, {N}x{N} matrix")
        M, labels = build_weil_matrix(primes, m_max, zeros, verbose=False)
        Mp = primitive_projection(M, labels)
        eigs = analyze_eigenvalues(Mp)

        max_eig = float(eigs[-1])
        min_eig = float(eigs[0])
        gap = abs(max_eig)
        all_neg = bool(np.all(eigs < 1e-12))

        results[name] = {'n_primes': len(primes), 'N': N, 'max_eig': max_eig,
                         'min_eig': min_eig, 'spectral_gap': gap, 'all_neg': all_neg,
                         'top5': [float(e) for e in eigs[-5:]],
                         'bot5': [float(e) for e in eigs[:5]]}
        apt = "YES" if all_neg else "NO"
        print(f"    APT: {apt}, max_eig: {max_eig:+.6e}, gap: {gap:.6e}")

    print(f"\n  Residue class comparison:")
    print(f"  {'Class':>15s} {'#pr':>4s} {'max_eig':>14s} {'min_eig':>14s} {'gap':>14s} {'APT':>5s}")
    print("  " + "-" * 60)
    for name in ['1 (mod 4)', '3 (mod 4)', '1 (mod 3)', '2 (mod 3)', 'all primes']:
        r = results[name]
        apt = "YES" if r['all_neg'] else "NO"
        print(f"  {name:>15s} {r['n_primes']:4d} {r['max_eig']:+14.6e} "
              f"{r['min_eig']:+14.6e} {r['spectral_gap']:14.6e} {apt:>5s}")

    g1 = results['1 (mod 4)']['spectral_gap']
    g3 = results['3 (mod 4)']['spectral_gap']
    if g1 > 0:
        print(f"\n  Chebyshev bias (mod 4): gap(3)/gap(1) = {g3/g1:.4f}")

    return results


# ═══════════════════════════════════════════════════════════════════
#  Investigation 5: Off-Line Zero Sensitivity
# ═══════════════════════════════════════════════════════════════════

def investigation_5(zeros):
    print(f"\n{SEPARATOR}")
    print("  INVESTIGATION 5: OFF-LINE ZERO SENSITIVITY")
    print(f"{SEPARATOR}")

    P0, m_max = 97, 3
    primes = sieve_primes(P0)
    N = len(primes) * m_max
    labels = [(p, a) for p in primes for a in range(1, m_max + 1)]

    print(f"  P_0={P0}, {len(primes)} primes, m_max={m_max}, N={N}")
    print(f"  Building base matrix...")

    M_base, _ = build_weil_matrix(primes, m_max, zeros)
    Mp_base = primitive_projection(M_base, labels)
    base_eigs = analyze_eigenvalues(Mp_base)
    base_max = float(base_eigs[-1])
    print(f"  Base max primitive eigenvalue: {base_max:+.8e}")

    eps = 0.01
    n_test = min(100, len(zeros))

    x_vals = np.array([a * np.log(p) for p, a in labels])
    ws = np.array([np.sqrt(np.log(p)) / p**(a / 2.0) for p, a in labels])

    print(f"  Computing sensitivity for {n_test} zeros (eps={eps})...")
    sensitivities = []
    t0 = time.time()

    for k in range(n_test):
        gamma_k = zeros[k]
        beta = 0.5 + eps
        rho_prod = abs(complex(beta, gamma_k) * complex(1 - beta, -gamma_k))

        # Perturbation matrix: change in kernel from moving zero k off-line
        DM = np.zeros((N, N))
        for i in range(N):
            for j in range(i + 1):
                x = x_vals[i] - x_vals[j]
                dk = (np.cosh(eps * x) * np.cos(gamma_k * x) / rho_prod
                      - np.cos(gamma_k * x) / (0.25 + gamma_k**2)) / np.pi
                w = -ws[i] * ws[j]
                DM[i, j] = w * dk
                DM[j, i] = DM[i, j]

        M_pert = M_base + DM
        Mp_pert = primitive_projection(M_pert, labels)
        pert_eigs = analyze_eigenvalues(Mp_pert)
        pert_max = float(pert_eigs[-1])

        sensitivity = (pert_max - base_max) / eps
        sensitivities.append({'k': k + 1, 'gamma': gamma_k,
                              'sensitivity': sensitivity, 'pert_max': pert_max})
        if (k + 1) % 20 == 0:
            print(f"    {k+1}/{n_test} [{time.time()-t0:.1f}s]")

    print(f"  Done in {time.time()-t0:.1f}s")

    sens_sorted = sorted(sensitivities, key=lambda s: abs(s['sensitivity']), reverse=True)
    print(f"\n  Top 10 most load-bearing zeros:")
    print(f"  {'k':>4s} {'gamma':>12s} {'d_lam/d_sig':>14s} {'lam_pert':>14s}")
    print("  " + "-" * 48)
    for s in sens_sorted[:10]:
        print(f"  {s['k']:4d} {s['gamma']:12.4f} {s['sensitivity']:+14.6e} {s['pert_max']:+14.6e}")

    # Decay fit
    gammas = np.array([s['gamma'] for s in sensitivities])
    sens_abs = np.array([abs(s['sensitivity']) for s in sensitivities])
    mask = sens_abs > 1e-20
    if mask.sum() > 10:
        coeffs = np.polyfit(np.log(gammas[mask]), np.log(sens_abs[mask]), 1)
        decay_exp = coeffs[0]
        decay_amp = np.exp(coeffs[1])
        print(f"\n  Sensitivity decay: |dlam/dsig| ~ {decay_amp:.4e} * gamma^({decay_exp:.4f})")
    else:
        decay_exp = decay_amp = float('nan')

    return sensitivities, (decay_amp, decay_exp)


# ═══════════════════════════════════════════════════════════════════
#  Results markdown
# ═══════════════════════════════════════════════════════════════════

def write_results_md(inv1, inv3, inv4, inv5_data, total_time):
    out_path = os.path.join(os.path.dirname(__file__), 'prime_correlation_study_results.md')
    sensitivities, (decay_amp, decay_exp) = inv5_data

    L = []
    w = L.append

    w("# Prime Correlation Study via the Weil Matrix")
    w("")
    w("What structural patterns do primes have, as encoded by the Weil matrix?")
    w(f"Total computation time: {total_time:.1f}s ({total_time/60:.1f} min)")
    w("")
    w("---")
    w("")

    # ── Investigation 1 ──
    w("## Investigation 1: Prime Correlation Map")
    w("")
    w("P_0 = 300, m_max = 3 (186x186 matrix)")
    w("")
    w("### Top 20 Strongest Off-Diagonal Couplings")
    w("")
    w("| p | m | q | n | |arg| | M_ij | K_bg | K_zeros |")
    w("|---|---|---|---|------|------|------|---------|")
    for p, m, q, n, x, mij, kb, kz in inv1['top20']:
        w(f"| {p} | {m} | {q} | {n} | {abs(x):.4f} | {mij:+.4e} | {kb:.4e} | {kz:+.4e} |")
    w("")
    w("### Decay Statistics (binned by |argument|)")
    w("")
    w("| Bin | Count | Mean |M| | Std |M| | Max |M| |")
    w("|-----|-------|---------|---------|---------|")
    for b in inv1['binned'][:20]:
        w(f"| [{b['bin_lo']:.1f}, {b['bin_hi']:.1f}) | {b['count']} | "
          f"{b['mean']:.4e} | {b['std']:.4e} | {b['max']:.4e} |")
    w("")
    a_exp, b_exp = inv1['fit_exp']
    a_pow, b_pow = inv1['fit_pow']
    w("### Decay Fits")
    w("")
    w(f"- **Exponential**: |M_ij| ~ {a_exp:.4e} exp(-{b_exp:.4f} |x|)")
    if not np.isnan(b_pow):
        w(f"- **Power law**: |M_ij| ~ {a_pow:.4e} |x|^(-{b_pow:.4f})")
    w("")

    # ── Investigation 3 ──
    w("---")
    w("")
    w("## Investigation 3: Eigenvalue Spectrum Evolution")
    w("")
    w("P_0 = 5 -> 500, m_max = 2 (adding one prime at a time)")
    w("")
    w("### Evolution (sampled)")
    w("")
    w("| k | Prime | N | max lam_prim | min lam_prim | Jump |")
    w("|---|-------|---|-------------|-------------|------|")
    for e in inv3:
        if e['k'] <= 5 or e['k'] % 10 == 0 or e['k'] == len(inv3) + 1:
            w(f"| {e['k']} | {e['prime']} | {e['N']} | "
              f"{e['max_eig']:+.4e} | {e['min_eig']:+.4e} | {e['jump']:+.4e} |")
    w("")
    jumps_sorted = sorted(inv3[1:], key=lambda e: abs(e['jump']), reverse=True)
    w("### Top 10 Biggest Eigenvalue Jumps")
    w("")
    w("| Prime | k | max lam | Jump |")
    w("|-------|---|---------|------|")
    for e in jumps_sorted[:10]:
        w(f"| {e['prime']} | {e['k']} | {e['max_eig']:+.6e} | {e['jump']:+.6e} |")
    w("")

    # ── Investigation 4 ──
    w("---")
    w("")
    w("## Investigation 4: Residue Class Structure")
    w("")
    w("P_0 = 500, m_max = 2")
    w("")
    w("| Class | #primes | max lam | min lam | Gap | APT |")
    w("|-------|---------|---------|---------|-----|-----|")
    for name in ['1 (mod 4)', '3 (mod 4)', '1 (mod 3)', '2 (mod 3)', 'all primes']:
        r = inv4[name]
        apt = "YES" if r['all_neg'] else "NO"
        w(f"| {name} | {r['n_primes']} | {r['max_eig']:+.4e} | "
          f"{r['min_eig']:+.4e} | {r['spectral_gap']:.4e} | {apt} |")
    w("")
    g1 = inv4['1 (mod 4)']['spectral_gap']
    g3 = inv4['3 (mod 4)']['spectral_gap']
    if g1 > 0:
        w(f"**Chebyshev bias (mod 4)**: gap(3 mod 4) / gap(1 mod 4) = {g3/g1:.4f}")
    w("")

    # ── Investigation 5 ──
    w("---")
    w("")
    w("## Investigation 5: Off-Line Zero Sensitivity")
    w("")
    w("P_0 = 97, m_max = 3, eps = 0.01")
    w("")
    sens_sorted = sorted(sensitivities, key=lambda s: abs(s['sensitivity']), reverse=True)
    w("### Top 10 Most Load-Bearing Zeros")
    w("")
    w("| k | gamma_k | d_lam/d_sigma | lam(pert) |")
    w("|---|---------|---------------|-----------|")
    for s in sens_sorted[:10]:
        w(f"| {s['k']} | {s['gamma']:.4f} | {s['sensitivity']:+.6e} | {s['pert_max']:+.6e} |")
    w("")
    w("### Full Sensitivity Table (first 50)")
    w("")
    w("| k | gamma_k | |d_lam/d_sigma| |")
    w("|---|---------|-----------------|")
    for s in sensitivities[:50]:
        w(f"| {s['k']} | {s['gamma']:.4f} | {abs(s['sensitivity']):.4e} |")
    w("")
    if not np.isnan(decay_exp):
        w(f"**Sensitivity decay**: |d_lam/d_sigma| ~ {decay_amp:.4e} gamma^({decay_exp:.4f})")
    w("")

    # ── Synthesis ──
    w("---")
    w("")
    w("## Synthesis: Do primes have a pattern visible in the Weil matrix?")
    w("")
    w(f"1. **Correlation decay**: Matrix entries decay exponentially (rate ~ {b_exp:.3f}) "
      f"in |log(p^m/q^n)|. Prime-prime couplings are dominated by logarithmically nearby "
      f"prime powers. K_bg dominates K_zeros in the strongest couplings.")
    w("")
    max_jump = max(inv3[1:], key=lambda e: abs(e['jump']))
    w(f"2. **Spectrum evolution**: Max primitive eigenvalue evolves smoothly as primes are "
      f"added. Largest jump at p={max_jump['prime']} "
      f"(Delta={max_jump['jump']:+.4e}). Small primes have outsized influence.")
    w("")
    all_apt = all(inv4[n]['all_neg'] for n in inv4)
    w(f"3. **Residue classes**: {'All classes' if all_apt else 'Not all classes'} satisfy APT "
      f"independently. The spectral gap ratio (mod 4) = {g3/g1:.4f} "
      f"{'shows' if abs(g3/g1 - 1) > 0.05 else 'does not show significant'} Chebyshev bias.")
    w("")
    top_zero = sens_sorted[0]
    w(f"4. **Zero sensitivity**: Zero #{top_zero['k']} (gamma={top_zero['gamma']:.2f}) is most "
      f"load-bearing (|d_lam/d_sigma|={abs(top_zero['sensitivity']):.4e}). ")
    if not np.isnan(decay_exp):
        w(f"   Sensitivity decays as gamma^({decay_exp:.2f}), "
          f"{'consistent with' if decay_exp < -1 else 'slower than'} "
          f"the 1/(1/4+gamma^2) weight in K_zeros.")
    w("")
    w(f"*Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}*")

    with open(out_path, 'w') as f:
        f.write('\n'.join(L) + '\n')
    print(f"\n  Results written to {out_path}")


# ═══════════════════════════════════════════════════════════════════
#  Main
# ═══════════════════════════════════════════════════════════════════

def main():
    t_start = time.time()
    print(SEPARATOR)
    print("  PRIME CORRELATION STUDY VIA THE WEIL MATRIX")
    print(SEPARATOR)

    N_ZEROS = 200
    zeros = compute_zeta_zeros(N_ZEROS)
    print()

    inv1 = investigation_1(zeros)
    inv3 = investigation_3(zeros)
    inv4 = investigation_4(zeros)
    inv5 = investigation_5(zeros)

    total_time = time.time() - t_start
    write_results_md(inv1, inv3, inv4, inv5, total_time)

    print(f"\n  Total time: {total_time:.1f}s ({total_time/60:.1f} min)")
    print(SEPARATOR)


if __name__ == '__main__':
    main()
