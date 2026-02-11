#!/usr/bin/env python3
"""
zero_visualization.py — Visualizations of the Riemann zeta function and its zeros.

Produces four plots:
  1. Zeros in the critical strip
  2. |ζ(1/2 + it)| along the critical line
  3. Heat map of |ζ(σ + it)| in the critical strip
  4. Zero spacing distribution compared to GUE prediction
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import mpmath
from pathlib import Path

# Output directory
OUT_DIR = Path(__file__).parent / 'plots'
OUT_DIR.mkdir(exist_ok=True)


def load_or_compute_zeros(n_zeros=200):
    """Load precomputed zeros or compute them."""
    data_file = Path(__file__).parent / 'zeros_data.txt'
    zeros_t = []

    if data_file.exists():
        print("Loading precomputed zeros...")
        with open(data_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split()
                if len(parts) >= 2:
                    zeros_t.append(float(parts[1]))
                    if len(zeros_t) >= n_zeros:
                        break

    if len(zeros_t) < n_zeros:
        print(f"Computing {n_zeros} zeros...")
        mpmath.mp.dps = 30
        zeros_t = []
        for k in range(1, n_zeros + 1):
            z = mpmath.zetazero(k)
            zeros_t.append(float(z.imag))
            if k % 50 == 0:
                print(f"  Computed {k}/{n_zeros} zeros")

    return zeros_t


def plot_zeros_critical_strip(zeros_t, n_show=100):
    """Plot zeros in the critical strip, showing they all sit at Re(s) = 1/2."""
    fig, ax = plt.subplots(1, 1, figsize=(8, 12))

    t_vals = zeros_t[:n_show]

    # Plot zeros as points at (1/2, t)
    ax.scatter([0.5] * len(t_vals), t_vals, s=15, c='red', zorder=5, label='Non-trivial zeros')

    # Draw the critical line
    ax.axvline(x=0.5, color='blue', linestyle='--', alpha=0.5, label='Critical line Re(s)=1/2')

    # Draw the critical strip boundaries
    ax.axvline(x=0, color='gray', linestyle='-', alpha=0.3, label='Critical strip boundary')
    ax.axvline(x=1, color='gray', linestyle='-', alpha=0.3)

    # Shade the critical strip
    ax.axvspan(0, 1, alpha=0.05, color='blue')

    ax.set_xlabel('Re(s) = σ', fontsize=14)
    ax.set_ylabel('Im(s) = t', fontsize=14)
    ax.set_title(f'First {n_show} Non-trivial Zeros of ζ(s)\nin the Critical Strip', fontsize=16)
    ax.legend(loc='upper left', fontsize=11)
    ax.set_xlim(-0.5, 1.5)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(OUT_DIR / 'zeros_critical_strip.png', dpi=150)
    print(f"  Saved: {OUT_DIR / 'zeros_critical_strip.png'}")
    plt.close()


def plot_zeta_critical_line(t_max=100, n_points=2000):
    """Plot |ζ(1/2 + it)| along the critical line."""
    fig, axes = plt.subplots(2, 1, figsize=(14, 10))

    mpmath.mp.dps = 15
    t_vals = np.linspace(0.5, t_max, n_points)
    zeta_vals = []

    print("  Computing ζ(1/2+it) along critical line...")
    for t in t_vals:
        z = complex(mpmath.zeta(0.5 + 1j * t))
        zeta_vals.append(z)

    magnitudes = np.array([abs(z) for z in zeta_vals])
    real_parts = np.array([z.real for z in zeta_vals])
    imag_parts = np.array([z.imag for z in zeta_vals])

    # Plot 1: Magnitude
    axes[0].plot(t_vals, magnitudes, 'b-', linewidth=0.5)
    axes[0].fill_between(t_vals, 0, magnitudes, alpha=0.1, color='blue')
    axes[0].set_xlabel('t', fontsize=12)
    axes[0].set_ylabel('|ζ(1/2 + it)|', fontsize=12)
    axes[0].set_title('Magnitude of ζ(s) on the Critical Line', fontsize=14)
    axes[0].grid(True, alpha=0.3)

    # Mark zeros (where magnitude is small)
    zero_indices = []
    for i in range(1, len(magnitudes) - 1):
        if magnitudes[i] < magnitudes[i-1] and magnitudes[i] < magnitudes[i+1] and magnitudes[i] < 0.5:
            zero_indices.append(i)
    if zero_indices:
        axes[0].scatter(t_vals[zero_indices], magnitudes[zero_indices],
                       c='red', s=20, zorder=5, label='Approximate zero locations')
        axes[0].legend(fontsize=10)

    # Plot 2: Real and imaginary parts (Hardy's Z-function style)
    axes[1].plot(t_vals, real_parts, 'b-', linewidth=0.5, label='Re ζ(1/2+it)', alpha=0.7)
    axes[1].plot(t_vals, imag_parts, 'r-', linewidth=0.5, label='Im ζ(1/2+it)', alpha=0.7)
    axes[1].axhline(y=0, color='k', linewidth=0.5)
    axes[1].set_xlabel('t', fontsize=12)
    axes[1].set_ylabel('Value', fontsize=12)
    axes[1].set_title('Real and Imaginary Parts of ζ(1/2+it)', fontsize=14)
    axes[1].legend(fontsize=10)
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(OUT_DIR / 'zeta_critical_line.png', dpi=150)
    print(f"  Saved: {OUT_DIR / 'zeta_critical_line.png'}")
    plt.close()


def plot_heatmap_critical_strip(sigma_range=(-2, 3), t_range=(0, 50),
                                 n_sigma=200, n_t=300):
    """Heat map of |ζ(σ + it)| in the critical strip region."""
    fig, ax = plt.subplots(1, 1, figsize=(10, 12))

    mpmath.mp.dps = 15
    sigma_vals = np.linspace(sigma_range[0], sigma_range[1], n_sigma)
    t_vals = np.linspace(t_range[0], t_range[1], n_t)

    print("  Computing ζ(σ+it) on grid for heatmap...")
    Z = np.zeros((n_t, n_sigma))
    for i, t in enumerate(t_vals):
        if i % 50 == 0:
            print(f"    Row {i}/{n_t}")
        for j, sigma in enumerate(sigma_vals):
            try:
                val = abs(complex(mpmath.zeta(sigma + 1j * t)))
                Z[i, j] = val
            except Exception:
                Z[i, j] = np.nan

    # Use log scale for better visualization
    Z_log = np.log10(np.clip(Z, 1e-10, None))

    im = ax.pcolormesh(sigma_vals, t_vals, Z_log,
                       cmap='inferno', shading='auto',
                       vmin=-2, vmax=2)

    # Mark the critical line
    ax.axvline(x=0.5, color='cyan', linestyle='--', alpha=0.7, linewidth=1.5,
               label='Critical line σ=1/2')

    # Mark critical strip boundaries
    ax.axvline(x=0, color='white', linestyle=':', alpha=0.5)
    ax.axvline(x=1, color='white', linestyle=':', alpha=0.5)

    cbar = plt.colorbar(im, ax=ax, label='log₁₀|ζ(σ+it)|')

    ax.set_xlabel('σ = Re(s)', fontsize=14)
    ax.set_ylabel('t = Im(s)', fontsize=14)
    ax.set_title('Heat Map of |ζ(σ+it)| in the Critical Strip Region', fontsize=15)
    ax.legend(loc='upper right', fontsize=11)

    plt.tight_layout()
    plt.savefig(OUT_DIR / 'zeta_heatmap.png', dpi=150)
    print(f"  Saved: {OUT_DIR / 'zeta_heatmap.png'}")
    plt.close()


def plot_spacing_distribution(zeros_t, n_bins=50):
    """Plot zero spacing distribution and compare with GUE prediction."""
    fig, ax = plt.subplots(1, 1, figsize=(10, 7))

    # Compute normalized spacings
    spacings = np.diff(zeros_t)

    # Normalize: for zeros near height T, the average spacing is 2π/log(T/(2π))
    # We use a simple normalization by the local average
    # Window-based normalization for better accuracy
    window = 20
    normalized_spacings = []
    for i in range(len(spacings)):
        lo = max(0, i - window)
        hi = min(len(spacings), i + window)
        local_avg = np.mean(spacings[lo:hi])
        normalized_spacings.append(spacings[i] / local_avg)

    normalized_spacings = np.array(normalized_spacings)

    # Histogram of normalized spacings
    ax.hist(normalized_spacings, bins=n_bins, density=True, alpha=0.6,
            color='steelblue', edgecolor='navy', label='Computed zero spacings')

    # GUE prediction (Wigner surmise for GUE)
    # p(s) = (32/π²) s² exp(-4s²/π) for the nearest-neighbor spacing
    s = np.linspace(0, 4, 500)
    gue_wigner = (32 / np.pi**2) * s**2 * np.exp(-4 * s**2 / np.pi)
    ax.plot(s, gue_wigner, 'r-', linewidth=2.5, label='GUE Wigner surmise')

    # Poisson prediction for comparison
    poisson = np.exp(-s)
    ax.plot(s, poisson, 'g--', linewidth=2, alpha=0.7, label='Poisson (independent)')

    ax.set_xlabel('Normalized spacing s', fontsize=14)
    ax.set_ylabel('Probability density', fontsize=14)
    ax.set_title('Zero Spacing Distribution vs GUE Prediction', fontsize=15)
    ax.legend(fontsize=12)
    ax.set_xlim(0, 3.5)
    ax.grid(True, alpha=0.3)

    # Add text annotation
    ax.text(2.2, 0.6, 'Zero repulsion:\nsmall spacings\nare suppressed\n(agrees with GUE)',
            fontsize=11, style='italic',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='lightyellow', alpha=0.8))

    plt.tight_layout()
    plt.savefig(OUT_DIR / 'spacing_distribution.png', dpi=150)
    print(f"  Saved: {OUT_DIR / 'spacing_distribution.png'}")
    plt.close()


def main():
    print("=" * 60)
    print("RIEMANN ZETA FUNCTION VISUALIZATIONS")
    print("=" * 60)

    # Load or compute zeros
    zeros_t = load_or_compute_zeros(500)
    print(f"Working with {len(zeros_t)} zeros\n")

    print("1. Plotting zeros in critical strip...")
    plot_zeros_critical_strip(zeros_t)

    print("\n2. Plotting ζ(1/2+it) along critical line...")
    plot_zeta_critical_line()

    print("\n3. Creating heat map of |ζ(σ+it)|...")
    plot_heatmap_critical_strip()

    print("\n4. Plotting spacing distribution...")
    plot_spacing_distribution(zeros_t)

    print("\n" + "=" * 60)
    print(f"All plots saved to {OUT_DIR}")
    print("=" * 60)


if __name__ == '__main__':
    main()
