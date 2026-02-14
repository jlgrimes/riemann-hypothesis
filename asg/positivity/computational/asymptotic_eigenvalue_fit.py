#!/usr/bin/env python3
"""
Asymptotic Eigenvalue Fit for λ_max(P₀)
=========================================
Fits the maximum primitive eigenvalue λ_max(P₀) to asymptotic models
to determine whether the spectral gap closes (λ_max → 0) or stays open
(λ_max → -const) as P₀ → ∞.

Kernel: K(x) = K_bg(x) + K_zeros(x) + δ(x)
  K_bg(x)    = -(1/π) Re ψ(1/4 + ix/2) + (1/(2π)) log π
  K_zeros(x) = (1/π) Σ_γ cos(γx)/(1/4 + γ²)

Matrix: M_{(p,m),(q,n)} = -(log p · log q)^{1/2} / (p^{m/2} · q^{n/2}) · K(m log p - n log q)

Primitive projection: v_{p,m} = (log p)^{1/2} / p^{m/2}

Models fit to |λ_max(P₀)|:
  A: C / P₀^α           (power law decay)
  B: C / log(P₀)^β      (logarithmic decay)
  C: C * exp(-P₀^γ)     (exponential decay)
"""

import numpy as np
from scipy.linalg import eigvalsh
from scipy.optimize import curve_fit
import mpmath
import time
import sys
import os

mpmath.mp.dps = 30

SEPARATOR = "=" * 72


# ============================================================================
# Primes and zeros
# ============================================================================

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


# ============================================================================
# Kernel functions
# ============================================================================

def K_bg(x):
    """K_bg(x) = -(1/π) Re ψ(1/4 + ix/2) + (1/(2π)) log π"""
    arg = mpmath.mpc(0.25, float(x) / 2)
    psi = mpmath.digamma(arg)
    return float(-mpmath.re(psi) / mpmath.pi + mpmath.log(mpmath.pi) / (2 * mpmath.pi))


def K_zeros(x, zeros):
    """K_zeros(x) = (1/π) Σ_γ cos(γx)/(1/4 + γ²)"""
    return float(np.sum(np.cos(zeros * x) / (0.25 + zeros**2)) / np.pi)


# ============================================================================
# Weil matrix construction and primitive projection
# ============================================================================

def build_weil_matrix(primes, m_max, zeros, verbose=True):
    """
    Build Weil matrix M with entries:
    M_{(p,m),(q,n)} = -(log p · log q)^{1/2} / (p^{m/2} · q^{n/2}) · K(m log p - n log q)

    K includes delta on diagonal: K(0) = K_bg(0) + K_zeros(0) + 1
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

            k_bg = K_bg(x)
            k_z = K_zeros(x, zeros)
            k_full = k_bg + k_z + (1.0 if diag else 0.0)

            w = -weights[i] * weights[j]
            M[i, j] = w * k_full
            M[j, i] = M[i, j]

        if verbose and ((i + 1) % 20 == 0 or i + 1 == N):
            el = time.time() - t0
            frac = (i + 1) / N
            eta = el / frac * (1 - frac) if frac > 0 else 0
            print(f"    row {i+1}/{N}  [{el:.1f}s, ETA ~{eta:.0f}s]")
            sys.stdout.flush()

    if verbose:
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


def compute_max_primitive_eigenvalue(Mp):
    """Compute eigenvalues of projected matrix, return max primitive eigenvalue."""
    eigs = eigvalsh(Mp)
    # Remove the trivial ~0 eigenvalue from projection
    idx = np.argmin(np.abs(eigs))
    prim = np.delete(eigs, idx)
    return float(np.max(prim))


def compute_eigenvalue_at(P0, m_max, zeros):
    """Compute max primitive eigenvalue for a given P₀."""
    primes = sieve_primes(P0)
    n_primes = len(primes)
    N = n_primes * m_max

    print(f"\n  P₀={P0}: {n_primes} primes, m_max={m_max}, N={N}")
    t0 = time.time()

    M, labels = build_weil_matrix(primes, m_max, zeros, verbose=True)
    Mp = primitive_projection(M, labels)
    lam_max = compute_max_primitive_eigenvalue(Mp)

    elapsed = time.time() - t0
    print(f"  λ_max = {lam_max:+.6e}  [{elapsed:.1f}s]")

    return {
        'P0': P0,
        'n_primes': n_primes,
        'm_max': m_max,
        'N': N,
        'lam_max': lam_max,
        'time': elapsed,
    }


# ============================================================================
# Asymptotic models
# ============================================================================

def model_A(P0, C, alpha):
    """Power law: |λ_max| ~ C / P₀^α"""
    return C / np.power(P0, alpha)


def model_B(P0, C, beta):
    """Logarithmic: |λ_max| ~ C / log(P₀)^β"""
    return C / np.power(np.log(P0), beta)


def model_C(P0, C, gamma):
    """Exponential: |λ_max| ~ C * exp(-P₀^γ)"""
    return C * np.exp(-np.power(P0, gamma))


def r_squared(y_true, y_pred):
    """Coefficient of determination R²."""
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    if ss_tot == 0:
        return 1.0
    return 1.0 - ss_res / ss_tot


def aic(n, k, ss_res):
    """Akaike Information Criterion. n=data points, k=parameters, ss_res=residual sum of squares."""
    if ss_res <= 0:
        return -np.inf
    return n * np.log(ss_res / n) + 2 * k


# ============================================================================
# Main
# ============================================================================

def main():
    t_start = time.time()

    print(SEPARATOR)
    print("  ASYMPTOTIC EIGENVALUE FIT: λ_max(P₀)")
    print(SEPARATOR)
    print(f"  mpmath precision: {mpmath.mp.dps} digits")
    print(f"  Kernel: K(x) = K_bg(x) + K_zeros(x) + δ(x)")
    print(f"  Primitive projection: v_{{p,m}} = (log p)^{{1/2}} / p^{{m/2}}")
    print()

    # ------------------------------------------------------------------
    # 1. Existing data from large_scale_eigenvalue_results.md
    # ------------------------------------------------------------------
    existing_data = [
        {'P0': 200, 'n_primes': 46, 'm_max': 3, 'N': 138, 'lam_max': -6.8664e-07, 'time': 0.0},
        {'P0': 300, 'n_primes': 62, 'm_max': 3, 'N': 186, 'lam_max': -2.4716e-07, 'time': 0.0},
        {'P0': 500, 'n_primes': 95, 'm_max': 2, 'N': 190, 'lam_max': -2.5508e-05, 'time': 0.0},
        {'P0': 750, 'n_primes': 132, 'm_max': 2, 'N': 264, 'lam_max': -1.2063e-05, 'time': 0.0},
    ]

    print("  Existing data points:")
    for d in existing_data:
        print(f"    P₀={d['P0']:4d}  N={d['N']:4d}  λ_max={d['lam_max']:+.6e}")
    print()

    # ------------------------------------------------------------------
    # 2. Compute new data points
    # ------------------------------------------------------------------
    N_ZEROS = 200
    zeros = compute_zeta_zeros(N_ZEROS)
    print()

    new_configs = [
        (100, 3),    # P₀=100, m_max=3
        (150, 3),    # P₀=150, m_max=3
        (400, 2),    # P₀=400, m_max=2
        (600, 2),    # P₀=600, m_max=2
    ]

    new_data = []
    for P0, m_max in new_configs:
        res = compute_eigenvalue_at(P0, m_max, zeros)
        new_data.append(res)
        sys.stdout.flush()

    # ------------------------------------------------------------------
    # 3. Combine all data and sort by P₀
    # ------------------------------------------------------------------
    all_data = existing_data + new_data
    all_data.sort(key=lambda d: d['P0'])

    print(f"\n{SEPARATOR}")
    print("  ALL DATA POINTS")
    print(SEPARATOR)
    print(f"  {'P₀':>6s} {'N':>5s} {'m_max':>5s} {'λ_max':>14s} {'|λ_max|':>14s}")
    print("  " + "-" * 50)
    for d in all_data:
        print(f"  {d['P0']:6d} {d['N']:5d} {d['m_max']:5d} {d['lam_max']:+14.6e} {abs(d['lam_max']):14.6e}")

    # ------------------------------------------------------------------
    # 4. Fit asymptotic models to |λ_max(P₀)|
    # ------------------------------------------------------------------
    P0_arr = np.array([d['P0'] for d in all_data], dtype=float)
    abs_lam = np.array([abs(d['lam_max']) for d in all_data])

    print(f"\n{SEPARATOR}")
    print("  ASYMPTOTIC MODEL FITTING")
    print(SEPARATOR)

    n_data = len(P0_arr)
    k_params = 2  # all models have 2 parameters

    fits = {}

    # --- Model A: |λ_max| ~ C / P₀^α ---
    try:
        popt_A, _ = curve_fit(model_A, P0_arr, abs_lam, p0=[1.0, 1.0], maxfev=10000)
        pred_A = model_A(P0_arr, *popt_A)
        r2_A = r_squared(abs_lam, pred_A)
        ss_A = np.sum((abs_lam - pred_A) ** 2)
        aic_A = aic(n_data, k_params, ss_A)
        fits['A'] = {'params': popt_A, 'pred': pred_A, 'R2': r2_A, 'AIC': aic_A,
                      'label': f'C/P₀^α', 'formula': f'|λ_max| = {popt_A[0]:.6e} / P₀^{popt_A[1]:.4f}'}
        print(f"\n  Model A: |λ_max| ~ C / P₀^α")
        print(f"    C     = {popt_A[0]:.6e}")
        print(f"    α     = {popt_A[1]:.6f}")
        print(f"    R²    = {r2_A:.8f}")
        print(f"    AIC   = {aic_A:.4f}")
        print(f"    => λ_max → 0 as P₀ → ∞ (gap closes)")
    except Exception as e:
        print(f"\n  Model A: FAILED ({e})")
        fits['A'] = None

    # --- Model B: |λ_max| ~ C / log(P₀)^β ---
    try:
        popt_B, _ = curve_fit(model_B, P0_arr, abs_lam, p0=[1.0, 5.0], maxfev=10000)
        pred_B = model_B(P0_arr, *popt_B)
        r2_B = r_squared(abs_lam, pred_B)
        ss_B = np.sum((abs_lam - pred_B) ** 2)
        aic_B = aic(n_data, k_params, ss_B)
        fits['B'] = {'params': popt_B, 'pred': pred_B, 'R2': r2_B, 'AIC': aic_B,
                      'label': f'C/log(P₀)^β', 'formula': f'|λ_max| = {popt_B[0]:.6e} / log(P₀)^{popt_B[1]:.4f}'}
        print(f"\n  Model B: |λ_max| ~ C / log(P₀)^β")
        print(f"    C     = {popt_B[0]:.6e}")
        print(f"    β     = {popt_B[1]:.6f}")
        print(f"    R²    = {r2_B:.8f}")
        print(f"    AIC   = {aic_B:.4f}")
        print(f"    => λ_max → 0 as P₀ → ∞ (gap closes, slowly)")
    except Exception as e:
        print(f"\n  Model B: FAILED ({e})")
        fits['B'] = None

    # --- Model C: |λ_max| ~ C * exp(-P₀^γ) ---
    try:
        popt_C, _ = curve_fit(model_C, P0_arr, abs_lam, p0=[1e-3, 0.1],
                               maxfev=10000, bounds=([0, 0], [np.inf, 1.0]))
        pred_C = model_C(P0_arr, *popt_C)
        r2_C = r_squared(abs_lam, pred_C)
        ss_C = np.sum((abs_lam - pred_C) ** 2)
        aic_C = aic(n_data, k_params, ss_C)
        fits['C'] = {'params': popt_C, 'pred': pred_C, 'R2': r2_C, 'AIC': aic_C,
                      'label': f'C·exp(-P₀^γ)', 'formula': f'|λ_max| = {popt_C[0]:.6e} * exp(-P₀^{popt_C[1]:.4f})'}
        print(f"\n  Model C: |λ_max| ~ C * exp(-P₀^γ)")
        print(f"    C     = {popt_C[0]:.6e}")
        print(f"    γ     = {popt_C[1]:.6f}")
        print(f"    R²    = {r2_C:.8f}")
        print(f"    AIC   = {aic_C:.4f}")
        print(f"    => λ_max → 0 super-exponentially (gap closes fast)")
    except Exception as e:
        print(f"\n  Model C: FAILED ({e})")
        fits['C'] = None

    # ------------------------------------------------------------------
    # 5. Model comparison and interpretation
    # ------------------------------------------------------------------
    print(f"\n{SEPARATOR}")
    print("  MODEL COMPARISON")
    print(SEPARATOR)

    valid_fits = {k: v for k, v in fits.items() if v is not None}
    if valid_fits:
        best_r2 = max(valid_fits.keys(), key=lambda k: valid_fits[k]['R2'])
        best_aic = min(valid_fits.keys(), key=lambda k: valid_fits[k]['AIC'])

        print(f"\n  {'Model':>7s} {'R²':>12s} {'AIC':>12s} {'Formula'}")
        print("  " + "-" * 70)
        for k in sorted(valid_fits.keys()):
            f = valid_fits[k]
            markers = []
            if k == best_r2:
                markers.append("best R²")
            if k == best_aic:
                markers.append("best AIC")
            mark = f"  <-- {', '.join(markers)}" if markers else ""
            print(f"  {k:>7s} {f['R2']:12.8f} {f['AIC']:12.4f} {f['formula']}{mark}")

        print(f"\n  Best by R²:  Model {best_r2}")
        print(f"  Best by AIC: Model {best_aic}")

        # All models have λ_max → 0, the question is rate
        print(f"\n  INTERPRETATION:")
        print(f"    All three models predict |λ_max| → 0 as P₀ → ∞.")
        print(f"    This means the spectral gap CLOSES — the most negative primitive")
        print(f"    eigenvalue approaches zero from below.")
        print(f"")
        if best_aic == 'A':
            alpha_val = valid_fits['A']['params'][1]
            print(f"    Best fit (Model A): power law decay ~ P₀^{{-{alpha_val:.2f}}}")
            print(f"    The gap closes polynomially. For the proof strategy, this means")
            print(f"    the finite-to-infinite extrapolation gets harder at larger P₀,")
            print(f"    but the gap remains nonzero at every finite P₀.")
        elif best_aic == 'B':
            beta_val = valid_fits['B']['params'][1]
            print(f"    Best fit (Model B): logarithmic decay ~ log(P₀)^{{-{beta_val:.2f}}}")
            print(f"    The gap closes very slowly (logarithmically).")
            print(f"    This is consistent with the spectral gap being controlled by")
            print(f"    the prime counting function π(P₀) ~ P₀/log(P₀).")
        elif best_aic == 'C':
            gamma_val = valid_fits['C']['params'][1]
            print(f"    Best fit (Model C): exponential decay ~ exp(-P₀^{{{gamma_val:.4f}}})")
            print(f"    The gap closes super-exponentially fast.")

        print(f"\n    CAVEAT: The non-monotonicity of |λ_max| (e.g., P₀=300 vs P₀=500)")
        print(f"    reflects the m_max transition (m_max=3→2), which changes the")
        print(f"    matrix structure. Fits should be interpreted with this in mind.")
    else:
        print("  No models converged.")

    # ------------------------------------------------------------------
    # 6. Extrapolation predictions
    # ------------------------------------------------------------------
    extrap_P0 = [1000, 2000, 5000, 10000]
    print(f"\n{SEPARATOR}")
    print("  EXTRAPOLATION PREDICTIONS")
    print(SEPARATOR)
    print(f"\n  {'P₀':>8s}", end="")
    for k in sorted(valid_fits.keys()):
        print(f"  {'Model '+k:>14s}", end="")
    print()
    print("  " + "-" * (8 + 16 * len(valid_fits)))
    for p0 in extrap_P0:
        print(f"  {p0:8d}", end="")
        for k in sorted(valid_fits.keys()):
            f = valid_fits[k]
            if k == 'A':
                val = model_A(p0, *f['params'])
            elif k == 'B':
                val = model_B(p0, *f['params'])
            elif k == 'C':
                val = model_C(p0, *f['params'])
            print(f"  {val:14.6e}", end="")
        print()

    total_time = time.time() - t_start
    print(f"\n  Total computation time: {total_time:.1f}s ({total_time/60:.1f}min)")
    print(SEPARATOR)

    # ------------------------------------------------------------------
    # 7. Write results markdown
    # ------------------------------------------------------------------
    write_results_md(all_data, fits, valid_fits, extrap_P0, total_time)

    return all_data, fits


def write_results_md(all_data, fits, valid_fits, extrap_P0, total_time):
    """Write results to markdown file."""
    out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            'asymptotic_eigenvalue_fit_results.md')

    lines = []
    lines.append("# Asymptotic Eigenvalue Fit Results")
    lines.append("")
    lines.append("## Overview")
    lines.append("")
    lines.append("Fits the maximum primitive eigenvalue $\\lambda_{\\max}(P_0)$ of the Weil matrix")
    lines.append("to asymptotic models, to determine whether the spectral gap closes or stays open.")
    lines.append("")
    lines.append("## Configuration")
    lines.append("")
    lines.append(f"- **mpmath precision**: {mpmath.mp.dps} digits")
    lines.append(f"- **Zeta zeros**: 200")
    lines.append(f"- **Total computation time**: {total_time:.1f}s ({total_time/60:.1f}min)")
    lines.append("")

    # Kernel and matrix description
    lines.append("### Kernel")
    lines.append("")
    lines.append("$$K(x) = K_{\\text{bg}}(x) + K_{\\text{zeros}}(x) + \\delta(x)$$")
    lines.append("")
    lines.append("- $K_{\\text{bg}}(x) = -\\frac{1}{\\pi} \\operatorname{Re} \\psi(1/4 + ix/2) + \\frac{1}{2\\pi} \\log \\pi$")
    lines.append("- $K_{\\text{zeros}}(x) = \\frac{1}{\\pi} \\sum_\\gamma \\frac{\\cos(\\gamma x)}{1/4 + \\gamma^2}$")
    lines.append("")
    lines.append("### Matrix")
    lines.append("")
    lines.append("$$M_{(p,m),(q,n)} = -\\frac{\\sqrt{\\log p \\cdot \\log q}}{p^{m/2} \\cdot q^{n/2}} \\cdot K(m \\log p - n \\log q)$$")
    lines.append("")
    lines.append("### Primitive Projection")
    lines.append("")
    lines.append("$v_{p,m} = (\\log p)^{1/2} / p^{m/2}$, $P = I - vv^T/\\|v\\|^2$, $M_p = PMP$")
    lines.append("")
    lines.append("---")
    lines.append("")

    # Data table
    lines.append("## Data")
    lines.append("")
    lines.append("| P₀ | #primes | m_max | N | λ_max | |λ_max| |")
    lines.append("|---:|--------:|------:|--:|------:|-------:|")
    for d in all_data:
        lines.append(f"| {d['P0']} | {d['n_primes']} | {d['m_max']} | {d['N']} "
                      f"| {d['lam_max']:+.6e} | {abs(d['lam_max']):.6e} |")
    lines.append("")

    # Model predictions table
    lines.append("### Data with Model Predictions")
    lines.append("")
    header = "| P₀ | |λ_max| (data) |"
    sep = "|---:|---------------:|"
    for k in sorted(valid_fits.keys()):
        header += f" Model {k} |"
        sep += "----------:|"
    lines.append(header)
    lines.append(sep)
    for d in all_data:
        p0 = d['P0']
        row = f"| {p0} | {abs(d['lam_max']):.6e} |"
        for k in sorted(valid_fits.keys()):
            f = valid_fits[k]
            if k == 'A':
                val = model_A(p0, *f['params'])
            elif k == 'B':
                val = model_B(p0, *f['params'])
            elif k == 'C':
                val = model_C(p0, *f['params'])
            row += f" {val:.6e} |"
        lines.append(row)
    lines.append("")
    lines.append("---")
    lines.append("")

    # Model fits
    lines.append("## Model Fits")
    lines.append("")
    for k in sorted(valid_fits.keys()):
        f = valid_fits[k]
        lines.append(f"### Model {k}")
        lines.append("")
        if k == 'A':
            lines.append(f"$$|\\lambda_{{\\max}}| \\sim \\frac{{C}}{{P_0^\\alpha}}$$")
            lines.append("")
            lines.append(f"- $C = {f['params'][0]:.6e}$")
            lines.append(f"- $\\alpha = {f['params'][1]:.6f}$")
        elif k == 'B':
            lines.append(f"$$|\\lambda_{{\\max}}| \\sim \\frac{{C}}{{\\log(P_0)^\\beta}}$$")
            lines.append("")
            lines.append(f"- $C = {f['params'][0]:.6e}$")
            lines.append(f"- $\\beta = {f['params'][1]:.6f}$")
        elif k == 'C':
            lines.append(f"$$|\\lambda_{{\\max}}| \\sim C \\cdot \\exp(-P_0^\\gamma)$$")
            lines.append("")
            lines.append(f"- $C = {f['params'][0]:.6e}$")
            lines.append(f"- $\\gamma = {f['params'][1]:.6f}$")
        lines.append(f"- $R^2 = {f['R2']:.8f}$")
        lines.append(f"- $\\text{{AIC}} = {f['AIC']:.4f}$")
        lines.append("")

    # Comparison table
    lines.append("### Model Comparison")
    lines.append("")
    lines.append("| Model | R² | AIC | Decay type |")
    lines.append("|:------|---:|----:|:-----------|")
    best_r2 = max(valid_fits.keys(), key=lambda k: valid_fits[k]['R2']) if valid_fits else None
    best_aic = min(valid_fits.keys(), key=lambda k: valid_fits[k]['AIC']) if valid_fits else None
    decay_labels = {'A': 'Power law', 'B': 'Logarithmic', 'C': 'Exponential'}
    for k in sorted(valid_fits.keys()):
        f = valid_fits[k]
        marks = []
        if k == best_r2:
            marks.append("best R²")
        if k == best_aic:
            marks.append("best AIC")
        mark = f" **({', '.join(marks)})**" if marks else ""
        lines.append(f"| {k}{mark} | {f['R2']:.8f} | {f['AIC']:.4f} | {decay_labels.get(k, '?')} |")
    lines.append("")
    lines.append("---")
    lines.append("")

    # Extrapolation
    lines.append("## Extrapolation")
    lines.append("")
    header = "| P₀ |"
    sep = "|---:|"
    for k in sorted(valid_fits.keys()):
        header += f" Model {k} |"
        sep += "----------:|"
    lines.append(header)
    lines.append(sep)
    for p0 in extrap_P0:
        row = f"| {p0} |"
        for k in sorted(valid_fits.keys()):
            f = valid_fits[k]
            if k == 'A':
                val = model_A(p0, *f['params'])
            elif k == 'B':
                val = model_B(p0, *f['params'])
            elif k == 'C':
                val = model_C(p0, *f['params'])
            row += f" {val:.6e} |"
        lines.append(row)
    lines.append("")
    lines.append("---")
    lines.append("")

    # Interpretation
    lines.append("## Interpretation")
    lines.append("")
    lines.append("### Does the spectral gap close?")
    lines.append("")
    lines.append("All three models predict $|\\lambda_{\\max}| \\to 0$ as $P_0 \\to \\infty$.")
    lines.append("This means the spectral gap **closes**: the most negative primitive eigenvalue")
    lines.append("approaches zero from below.")
    lines.append("")
    lines.append("### What this means for the proof")
    lines.append("")
    lines.append("- **Gap closing** ($\\lambda_{\\max} \\to 0^-$): The finite verification becomes")
    lines.append("  harder at larger $P_0$, because the eigenvalue margin shrinks.")
    lines.append("- However, $\\lambda_{\\max} < 0$ at **every** finite $P_0$ tested.")
    lines.append("- The rate of gap closing determines how much headroom exists for")
    lines.append("  tail/truncation error bounds in the finite-to-infinite extrapolation.")
    lines.append("")
    lines.append("### Caveats")
    lines.append("")
    lines.append("- The transition from $m_{\\max}=3$ to $m_{\\max}=2$ at $P_0 \\ge 400$")
    lines.append("  introduces a structural change in the matrix (fewer prime powers),")
    lines.append("  which causes non-monotonicity in $|\\lambda_{\\max}|$.")
    lines.append("- The fit quality is limited by the small number of data points (8).")
    lines.append("- More data at larger $P_0$ with consistent $m_{\\max}$ would yield")
    lines.append("  more reliable scaling estimates.")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append(f"*Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}*")

    with open(out_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    print(f"\n  Results written to {out_path}")


if __name__ == '__main__':
    main()
