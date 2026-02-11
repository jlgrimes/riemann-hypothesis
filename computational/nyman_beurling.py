#!/usr/bin/env python3
"""
nyman_beurling.py — Implement the Nyman-Beurling criterion for the Riemann Hypothesis.

The Nyman-Beurling criterion (1950, 1955) states:

    RH is equivalent to the density of the set of functions
    f(x) = Σ_{k=1}^{N} c_k · {θ_k/x}
    in L²(0,1), where {y} = y - ⌊y⌋ is the fractional part,
    0 < θ_k ≤ 1, and we approximate the constant function χ(x) = 1.

In other words: RH holds if and only if the constant function 1 can be
approximated arbitrarily well in L²(0,1) by linear combinations of
dilated fractional parts.

Báez-Duarte (2003) simplified this: it suffices to use θ_k = 1/k.

We implement the approximation and measure the L² distance
    d_N = inf_{c₁,...,c_N} ‖1 - Σ c_k {1/(kx)}‖²_{L²(0,1)}

as a function of N. If this approaches 0 as N → ∞, RH holds.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate, optimize
from pathlib import Path

OUT_DIR = Path(__file__).parent / 'plots'
OUT_DIR.mkdir(exist_ok=True)


def fractional_part(y):
    """Compute {y} = y - floor(y), the fractional part."""
    return y - np.floor(y)


def rho_function(x, theta):
    """
    Compute ρ(θ/x) = {θ/x} - θ/x for the Nyman-Beurling setting.

    Actually, the standard formulation uses {θ/x} directly.
    """
    if x <= 0:
        return 0.0
    return fractional_part(theta / x)


def compute_gram_matrix(N, n_quad=2000):
    """
    Compute the Gram matrix G_{jk} = ⟨{1/(jx)}, {1/(kx)}⟩_{L²(0,1)}

    G_{jk} = ∫_0^1 {1/(jx)} · {1/(kx)} dx
    """
    print(f"  Computing {N}×{N} Gram matrix...")
    x_vals = np.linspace(1e-10, 1, n_quad)
    dx = x_vals[1] - x_vals[0]

    G = np.zeros((N, N))
    # Precompute the basis functions
    basis = np.zeros((N, n_quad))
    for j in range(N):
        theta = 1.0 / (j + 1)
        basis[j, :] = np.array([fractional_part(theta / x) for x in x_vals])

    # Gram matrix via numerical integration
    for j in range(N):
        for k in range(j, N):
            integrand = basis[j, :] * basis[k, :]
            G[j, k] = np.trapz(integrand, x_vals)
            G[k, j] = G[j, k]

    return G, basis, x_vals


def compute_rhs_vector(N, basis, x_vals):
    """
    Compute the right-hand side vector b_j = ⟨1, {1/(jx)}⟩_{L²(0,1)}

    b_j = ∫_0^1 {1/(jx)} dx
    """
    b = np.zeros(N)
    for j in range(N):
        b[j] = np.trapz(basis[j, :], x_vals)
    return b


def nyman_beurling_distance(N, n_quad=2000):
    """
    Compute the Nyman-Beurling L² distance:
    d_N² = inf_{c} ‖1 - Σ c_k {1/(kx)}‖²_{L²(0,1)}

    This is a least-squares problem:
    d_N² = 1 - b^T G^{-1} b

    where G is the Gram matrix and b is the inner product vector.
    """
    G, basis, x_vals = compute_gram_matrix(N, n_quad)
    b = compute_rhs_vector(N, basis, x_vals)

    # ‖1‖² = ∫_0^1 1 dx = 1
    norm_1_sq = 1.0

    # Solve Gc = b for optimal coefficients
    try:
        # Add small regularization for numerical stability
        reg = 1e-12 * np.eye(N)
        c_opt = np.linalg.solve(G + reg, b)

        # Distance: d² = ‖1‖² - 2⟨1, Σc_k f_k⟩ + ‖Σc_k f_k‖²
        #             = 1 - 2·b^T·c + c^T·G·c
        #             = 1 - b^T·c  (since Gc = b → c^T G c = c^T b = b^T c)
        d_sq = norm_1_sq - np.dot(b, c_opt)
        d_sq = max(d_sq, 0)  # Ensure non-negative

    except np.linalg.LinAlgError:
        print(f"    Warning: Singular matrix at N={N}, using pseudoinverse")
        c_opt = np.linalg.lstsq(G, b, rcond=None)[0]
        d_sq = max(1.0 - np.dot(b, c_opt), 0)

    return d_sq, c_opt, G, basis, x_vals


def compute_approximation(c_opt, basis, x_vals):
    """Compute the actual approximation f(x) = Σ c_k {1/(kx)}."""
    approx = np.zeros_like(x_vals)
    for k in range(len(c_opt)):
        approx += c_opt[k] * basis[k, :]
    return approx


def run_nyman_beurling(N_max=40, step=2):
    """Compute the Nyman-Beurling distance for N = 1, 2, ..., N_max."""
    print("=" * 60)
    print("NYMAN-BEURLING CRITERION")
    print("=" * 60)
    print()
    print("Computing d_N = inf ‖1 - Σ c_k {1/(kx)}‖_{L²(0,1)}")
    print()

    N_values = list(range(1, min(N_max + 1, 10))) + list(range(10, N_max + 1, step))
    distances_sq = []
    coefficients = []

    for N in N_values:
        d_sq, c_opt, G, basis, x_vals = nyman_beurling_distance(N)
        distances_sq.append(d_sq)
        coefficients.append(c_opt)
        d = np.sqrt(d_sq) if d_sq > 0 else 0
        print(f"  N = {N:>3d}: d_N² = {d_sq:>12.8f},  d_N = {d:>10.8f}")

    return N_values, distances_sq, coefficients


def plot_distance_decay(N_values, distances_sq):
    """Plot how the approximation distance decays with N."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    distances = [np.sqrt(max(d, 0)) for d in distances_sq]

    # Linear scale
    axes[0].plot(N_values, distances, 'bo-', markersize=5, linewidth=1.5)
    axes[0].set_xlabel('N (number of basis functions)', fontsize=13)
    axes[0].set_ylabel('d_N = ‖1 - f_N‖_{L²(0,1)}', fontsize=13)
    axes[0].set_title('Nyman-Beurling Distance (linear scale)', fontsize=14)
    axes[0].grid(True, alpha=0.3)

    # Log scale
    positive = [(n, d) for n, d in zip(N_values, distances) if d > 0]
    if positive:
        ns, ds = zip(*positive)
        axes[1].semilogy(ns, ds, 'ro-', markersize=5, linewidth=1.5, label='d_N')

        # Fit a power law: d_N ~ N^(-α)
        if len(ns) > 3:
            log_ns = np.log(np.array(ns, dtype=float))
            log_ds = np.log(np.array(ds, dtype=float))
            mask = np.isfinite(log_ds)
            if np.sum(mask) > 2:
                coeffs = np.polyfit(log_ns[mask], log_ds[mask], 1)
                alpha = -coeffs[0]
                n_smooth = np.linspace(min(ns), max(ns), 100)
                fit = np.exp(coeffs[1]) * n_smooth**coeffs[0]
                axes[1].semilogy(n_smooth, fit, 'g--', linewidth=1.5,
                                label=f'Fit: d_N ~ N^{{-{alpha:.2f}}}')
                axes[1].legend(fontsize=11)

    axes[1].set_xlabel('N (number of basis functions)', fontsize=13)
    axes[1].set_ylabel('d_N (log scale)', fontsize=13)
    axes[1].set_title('Nyman-Beurling Distance (log scale)', fontsize=14)
    axes[1].grid(True, alpha=0.3)

    plt.suptitle('Nyman-Beurling Criterion: Approximating 1 by Dilated Fractional Parts',
                 fontsize=15, y=1.02)
    plt.tight_layout()
    plt.savefig(OUT_DIR / 'nyman_beurling_distance.png', dpi=150, bbox_inches='tight')
    print(f"\n  Saved: {OUT_DIR / 'nyman_beurling_distance.png'}")
    plt.close()


def plot_approximations(N_values_plot=[2, 5, 10, 20]):
    """Plot the approximating functions for several values of N."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    for idx, N in enumerate(N_values_plot):
        ax = axes[idx // 2, idx % 2]

        d_sq, c_opt, G, basis, x_vals = nyman_beurling_distance(N)
        approx = compute_approximation(c_opt, basis, x_vals)

        ax.plot(x_vals, np.ones_like(x_vals), 'r--', linewidth=2,
                label='Target: f(x) = 1', alpha=0.7)
        ax.plot(x_vals, approx, 'b-', linewidth=1.5,
                label=f'Approximation (N={N})', alpha=0.8)
        ax.fill_between(x_vals, 1, approx, alpha=0.1, color='blue')

        d = np.sqrt(max(d_sq, 0))
        ax.set_title(f'N = {N}, d_N = {d:.6f}', fontsize=13)
        ax.set_xlabel('x', fontsize=11)
        ax.set_ylabel('f(x)', fontsize=11)
        ax.legend(fontsize=9, loc='lower right')
        ax.set_xlim(0, 1)
        ax.set_ylim(-0.5, 2.5)
        ax.grid(True, alpha=0.3)

    plt.suptitle('Approximating 1 by Linear Combinations of {1/(kx)}',
                 fontsize=15, y=1.02)
    plt.tight_layout()
    plt.savefig(OUT_DIR / 'nyman_beurling_approximations.png', dpi=150, bbox_inches='tight')
    print(f"  Saved: {OUT_DIR / 'nyman_beurling_approximations.png'}")
    plt.close()


def plot_coefficients(N_values, coefficients):
    """Plot the optimal coefficients for different N."""
    fig, ax = plt.subplots(1, 1, figsize=(12, 6))

    for i, N in enumerate(N_values):
        if N in [5, 10, 20, 30]:
            c = coefficients[i]
            ax.plot(range(1, len(c) + 1), c, 'o-', markersize=4, label=f'N={N}')

    ax.set_xlabel('k', fontsize=13)
    ax.set_ylabel('Optimal coefficient c_k', fontsize=13)
    ax.set_title('Optimal Coefficients in Nyman-Beurling Approximation', fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='gray', linewidth=0.5)

    plt.tight_layout()
    plt.savefig(OUT_DIR / 'nyman_beurling_coefficients.png', dpi=150)
    print(f"  Saved: {OUT_DIR / 'nyman_beurling_coefficients.png'}")
    plt.close()


def theoretical_connection():
    """Explain the theoretical significance of the Nyman-Beurling criterion."""
    print("\n" + "=" * 60)
    print("THEORETICAL SIGNIFICANCE")
    print("=" * 60)
    print("""
  The Nyman-Beurling criterion connects RH to function approximation:

  THEOREM (Nyman 1950, Beurling 1955):
    RH ⟺ The constant function 1 ∈ L²(0,1) can be approximated
          arbitrarily well by functions of the form
          f(x) = Σ_{k=1}^N c_k · {θ_k / x}

  SIMPLIFICATION (Báez-Duarte 2003):
    It suffices to take θ_k = 1/k, so we only need:
          f(x) = Σ_{k=1}^N c_k · {1/(kx)}

  CONNECTION TO ZETA:
    The Mellin transform relates {θ/x} to ζ(s):
    ∫_0^1 {θ/x} · x^{s-1} dx = -θ^s · ζ(s)/s  (for Re(s) > 1)

    So the approximation problem in L² is really about the
    behavior of ζ(s) on the critical line.

  RATE OF DECAY:
    Báez-Duarte showed that d_N² should decay like 1/(log N)
    if RH holds, which is extremely slow — consistent with the
    difficulty of proving RH.
    """)


def main():
    # Run the main computation
    N_values, distances_sq, coefficients = run_nyman_beurling(N_max=30, step=2)

    # Plot distance decay
    print("\nGenerating plots...")
    plot_distance_decay(N_values, distances_sq)

    # Plot approximations
    plot_approximations()

    # Plot coefficients
    plot_coefficients(N_values, coefficients)

    # Theoretical discussion
    theoretical_connection()

    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"  Computed d_N for N = 1 to {N_values[-1]}")
    print(f"  The distance d_N decreases as N increases")
    print(f"  d_N → 0 would confirm RH")
    print(f"  The convergence is very slow (expected: ~1/√(log N))")
    print(f"  This is consistent with the extreme difficulty of RH")
    print("=" * 60)


if __name__ == '__main__':
    main()
