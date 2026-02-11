#!/usr/bin/env python3
"""
newman_constant.py — Numerical exploration of the de Bruijn-Newman constant Λ.

The de Bruijn-Newman constant Λ is defined through the family of functions:

    H_t(z) = ∫_0^∞ Φ(u) · e^{tu²} · cos(zu) du

where Φ is a specific function derived from the Jacobi theta function:

    Φ(u) = Σ_{n=1}^∞ (2π²n⁴e^{9u} - 3πn²e^{5u}) · exp(-πn²e^{4u})

Key facts:
  - H_0(z) has the same zeros as ξ(1/2 + iz), the Riemann xi function
  - Λ is the infimum of t such that H_t has only real zeros
  - RH is equivalent to Λ ≤ 0
  - de Bruijn proved Λ ≤ 1/2
  - Newman conjectured Λ ≥ 0
  - Brad Rodgers and Terence Tao proved Λ ≥ 0 in 2018 (!)
  - Polymath 15 showed Λ ≤ 0.22
  - Best known: Λ ≤ 0.2 (Platt-Trudgian, 2021)
  - RH is equivalent to Λ = 0

We explore how the zeros of H_t move as t varies from negative to positive.
"""

import numpy as np
import matplotlib.pyplot as plt
import mpmath
from pathlib import Path

mpmath.mp.dps = 30

OUT_DIR = Path(__file__).parent / 'plots'
OUT_DIR.mkdir(exist_ok=True)


def phi_function(u, n_terms=20):
    """
    Compute Φ(u) = Σ_{n=1}^{n_terms} (2π²n⁴e^{9u} - 3πn²e^{5u}) · exp(-πn²e^{4u})

    This is the kernel function in the integral representation of the xi function.
    """
    result = mpmath.mpf(0)
    for n in range(1, n_terms + 1):
        n2 = mpmath.mpf(n)**2
        n4 = n2**2
        exp4u = mpmath.exp(4 * u)
        exponential = mpmath.exp(-mpmath.pi * n2 * exp4u)
        term = (2 * mpmath.pi**2 * n4 * mpmath.exp(9 * u) -
                3 * mpmath.pi * n2 * mpmath.exp(5 * u)) * exponential
        result += term
    return result


def H_t_numerical(z, t, n_quad=200, u_max=5.0):
    """
    Numerically compute H_t(z) = ∫_0^∞ Φ(u) · e^{tu²} · cos(zu) du

    Uses Gauss-Legendre quadrature on [0, u_max].
    """
    z = mpmath.mpf(z)
    t = mpmath.mpf(t)

    def integrand(u):
        if u == 0:
            u = mpmath.mpf('1e-30')
        return phi_function(u) * mpmath.exp(t * u**2) * mpmath.cos(z * u)

    result = mpmath.quad(integrand, [0, u_max], maxdegree=7)
    return result


def find_H_t_zeros_approximate(t, z_max=50, n_points=500):
    """
    Find approximate zeros of H_t(z) by looking for sign changes.

    Since H_t is an even function for real t, we only search z > 0.
    """
    z_vals = np.linspace(0.5, z_max, n_points)
    H_vals = []

    for z in z_vals:
        h = float(H_t_numerical(z, t))
        H_vals.append(h)

    H_vals = np.array(H_vals)

    # Find sign changes
    zeros = []
    for i in range(len(H_vals) - 1):
        if H_vals[i] * H_vals[i+1] < 0:
            # Linear interpolation to find approximate zero
            z0 = z_vals[i] - H_vals[i] * (z_vals[i+1] - z_vals[i]) / (H_vals[i+1] - H_vals[i])
            zeros.append(z0)

    return zeros, z_vals, H_vals


def approximate_Ht_from_xi(z, t, n_zeros=50):
    """
    Approximate H_t using the product formula and known zeta zeros.

    If ρ_k = 1/2 + iγ_k are zeta zeros, then H_0 has zeros at γ_k.
    Under the heat flow, the zeros evolve according to:
        dγ_k/dt = -Σ_{j≠k} 2/(γ_k - γ_j)

    For small t, the zeros move approximately as:
        γ_k(t) ≈ γ_k(0) + t · v_k

    where v_k = -Σ_{j≠k} 2/(γ_k(0) - γ_j(0))
    """
    mpmath.mp.dps = 30

    # Get zeta zeros
    gamma = []
    for k in range(1, n_zeros + 1):
        zero = mpmath.zetazero(k)
        gamma.append(float(zero.imag))

    gamma = np.array(gamma)

    # Compute velocity of each zero
    velocities = np.zeros(n_zeros)
    for k in range(n_zeros):
        v = 0.0
        for j in range(n_zeros):
            if j != k:
                diff = gamma[k] - gamma[j]
                if abs(diff) > 1e-10:
                    v -= 2.0 / diff
        velocities[k] = v

    # Evolved zeros (linear approximation for small t)
    gamma_t = gamma + t * velocities

    return gamma, gamma_t, velocities


def plot_zero_evolution():
    """Plot how zeros evolve under the heat flow as t varies."""
    print("Computing zero evolution under heat flow...")

    fig, axes = plt.subplots(1, 2, figsize=(16, 8))

    n_zeros = 30
    t_values = np.linspace(-0.1, 0.3, 50)

    # Get initial zeros and velocities
    gamma_0, _, velocities = approximate_Ht_from_xi(0, 0, n_zeros=n_zeros)

    # Track each zero as t varies
    for k in range(n_zeros):
        trajectory_real = gamma_0[k] + t_values * velocities[k]
        # For the imaginary part, zeros stay real when RH holds and t <= 0
        # They could become complex if t < Λ
        trajectory_imag = np.zeros_like(t_values)  # Zero imaginary part when on real axis

        axes[0].plot(t_values, trajectory_real, 'b-', linewidth=0.8, alpha=0.6)

    axes[0].axvline(x=0, color='red', linestyle='--', linewidth=1.5, label='t = 0 (RH)')
    axes[0].set_xlabel('t (deformation parameter)', fontsize=13)
    axes[0].set_ylabel('Zero location (real part)', fontsize=13)
    axes[0].set_title('Zeros of H_t(z) Under Heat Flow', fontsize=14)
    axes[0].legend(fontsize=11)
    axes[0].grid(True, alpha=0.3)

    # Plot the velocities
    axes[1].bar(range(1, n_zeros + 1), velocities, color='steelblue', alpha=0.7)
    axes[1].set_xlabel('Zero index k', fontsize=13)
    axes[1].set_ylabel('Velocity dγ_k/dt', fontsize=13)
    axes[1].set_title('Zero Velocities Under Heat Flow', fontsize=14)
    axes[1].grid(True, alpha=0.3)

    plt.suptitle('de Bruijn-Newman Constant: Zero Dynamics', fontsize=16, y=1.02)
    plt.tight_layout()
    plt.savefig(OUT_DIR / 'newman_zero_evolution.png', dpi=150, bbox_inches='tight')
    print(f"  Saved: {OUT_DIR / 'newman_zero_evolution.png'}")
    plt.close()


def plot_zero_spacing_under_flow():
    """Show how zero spacings change under the heat flow."""
    print("Analyzing spacing changes under heat flow...")

    fig, ax = plt.subplots(1, 1, figsize=(10, 7))

    n_zeros = 50
    gamma_0, _, velocities = approximate_Ht_from_xi(0, 0, n_zeros=n_zeros)

    t_values = np.linspace(-0.2, 0.3, 100)
    min_spacings = []

    for t in t_values:
        gamma_t = gamma_0 + t * velocities
        gamma_t_sorted = np.sort(gamma_t)
        spacings = np.diff(gamma_t_sorted)
        min_spacings.append(min(spacings) if len(spacings) > 0 else 0)

    ax.plot(t_values, min_spacings, 'b-', linewidth=2)
    ax.axvline(x=0, color='red', linestyle='--', label='t = 0 (RH ⟺ Λ = 0)')
    ax.axhline(y=0, color='gray', linestyle='-', alpha=0.5)
    ax.set_xlabel('t (deformation parameter)', fontsize=13)
    ax.set_ylabel('Minimum spacing between consecutive zeros', fontsize=13)
    ax.set_title('Minimum Zero Spacing Under Heat Flow\n'
                 '(Zeros collide when spacing → 0, then go off the real axis)', fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(OUT_DIR / 'newman_spacing_flow.png', dpi=150)
    print(f"  Saved: {OUT_DIR / 'newman_spacing_flow.png'}")
    plt.close()


def plot_Ht_function():
    """Plot H_t(z) for various values of t."""
    print("Computing H_t(z) for visualization (this may take a while)...")

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # We approximate H_t using the product over known zeros
    n_zeros = 30
    gamma_0, _, velocities = approximate_Ht_from_xi(0, 0, n_zeros=n_zeros)

    t_values_plot = [-0.1, 0.0, 0.1, 0.2]
    z_vals = np.linspace(0.5, 50, 500)

    for idx, t in enumerate(t_values_plot):
        ax = axes[idx // 2, idx % 2]

        # Approximate H_t(z) using evolved zeros
        gamma_t = gamma_0 + t * velocities
        H_approx = np.ones(len(z_vals))
        for gk in gamma_t:
            H_approx *= (1 - (z_vals / gk)**2)

        ax.plot(z_vals, H_approx, 'b-', linewidth=1)
        ax.axhline(y=0, color='gray', linewidth=0.5)
        ax.set_xlabel('z', fontsize=11)
        ax.set_ylabel('H_t(z) (approx.)', fontsize=11)
        ax.set_title(f't = {t:.1f}', fontsize=13)
        ax.grid(True, alpha=0.3)
        ax.set_ylim(-2, 2)

        # Mark approximate zero locations
        for gk in gamma_t:
            if 0 < gk < 50:
                ax.axvline(x=gk, color='red', alpha=0.3, linewidth=0.5)

    plt.suptitle('H_t(z) for Different Values of t\n'
                 '(Red lines mark zero locations)', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(OUT_DIR / 'newman_Ht_plots.png', dpi=150, bbox_inches='tight')
    print(f"  Saved: {OUT_DIR / 'newman_Ht_plots.png'}")
    plt.close()


def analyze_newman_bounds():
    """Summarize known bounds on Λ and what they mean for RH."""
    print("\n" + "=" * 60)
    print("KNOWN BOUNDS ON THE DE BRUIJN-NEWMAN CONSTANT Λ")
    print("=" * 60)

    bounds = [
        (1950, "de Bruijn", "Λ ≤ 1/2", "Upper bound from Fourier analysis"),
        (1976, "Newman", "Λ ≥ -∞ (conjecture: Λ ≥ 0)", "Newman's conjecture"),
        (1991, "Csordas-Norfolk-Varga", "Λ > -50", "First effective lower bound"),
        (1994, "Csordas-Odlyzko-Smith-Varga", "Λ > -5.895×10⁻⁹", "Improved lower bound"),
        (2009, "Saouter-Gourdon-Demichel", "Λ > -1.15×10⁻¹¹", "Further improvement"),
        (2018, "Rodgers-Tao", "Λ ≥ 0", "Newman's conjecture PROVED"),
        (2019, "Polymath 15", "Λ ≤ 0.22", "Collaborative upper bound"),
        (2021, "Platt-Trudgian", "Λ ≤ 0.2", "Best known upper bound"),
    ]

    for year, author, bound, note in bounds:
        print(f"  {year}: {bound:30s} ({author}) — {note}")

    print(f"\n  Current status: 0 ≤ Λ ≤ 0.2")
    print(f"  RH is equivalent to: Λ = 0")
    print(f"  If Λ > 0, then RH is FALSE")
    print(f"  If Λ = 0, then RH is TRUE")


def main():
    print("=" * 60)
    print("DE BRUIJN-NEWMAN CONSTANT Λ")
    print("=" * 60)
    print()
    print("The Riemann Hypothesis is equivalent to Λ = 0")
    print("where Λ is defined via the heat equation deformation")
    print("of the Riemann xi function.")
    print()

    # Known bounds
    analyze_newman_bounds()

    # Plot zero evolution
    print("\n" + "-" * 60)
    print("NUMERICAL EXPLORATION")
    print("-" * 60)

    plot_zero_evolution()
    plot_zero_spacing_under_flow()
    plot_Ht_function()

    # Estimate Λ from our computation
    print("\n" + "=" * 60)
    print("COMPUTATIONAL ESTIMATE")
    print("=" * 60)

    n_zeros = 30
    gamma_0, _, velocities = approximate_Ht_from_xi(0, 0, n_zeros=n_zeros)

    # Find the smallest t where two zeros collide (spacing → 0)
    # This gives an upper bound on Λ from our finite computation
    t_vals = np.linspace(0, 0.5, 10000)
    for t in t_vals:
        gamma_t = gamma_0 + t * velocities
        gamma_t_sorted = np.sort(gamma_t)
        min_spacing = min(np.diff(gamma_t_sorted))
        if min_spacing < 0.01:
            print(f"  First near-collision at t ≈ {t:.4f}")
            print(f"  (From {n_zeros} zeros — very rough upper bound on Λ)")
            break
    else:
        print(f"  No collision detected for t ≤ 0.5 with {n_zeros} zeros")

    print()
    print("  NOTE: Precise computation of Λ requires millions of zeros")
    print("  and sophisticated algorithms. Our rough estimate using 30 zeros")
    print("  should not be taken as a serious bound.")

    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print("  The de Bruijn-Newman constant Λ satisfies 0 ≤ Λ ≤ 0.2")
    print("  RH ⟺ Λ = 0")
    print("  Rodgers-Tao (2018) proved Λ ≥ 0 (Newman's conjecture)")
    print("  The gap between 0 and 0.2 remains to be closed")
    print("  Our visualizations show how zeros evolve under the heat flow")
    print("  and illustrate the mechanism by which violations could occur")
    print("=" * 60)


if __name__ == '__main__':
    main()
