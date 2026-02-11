#!/usr/bin/env python3
"""
explicit_formula.py — Implement the explicit formula for the prime counting function.

Riemann's explicit formula expresses the prime counting function π(x) in terms
of the non-trivial zeros of the zeta function:

    π(x) = li(x) - Σ_ρ li(x^ρ) - log(2) + ∫_x^∞ dt/(t(t²-1)log(t))

where the sum runs over all non-trivial zeros ρ of ζ(s), and
li(x) = ∫_0^x dt/log(t) is the logarithmic integral.

More precisely, we use the closely related formula for the Chebyshev function:

    ψ(x) = x - Σ_ρ x^ρ/ρ - log(2π) - (1/2)log(1-x^{-2})

and the von Mangoldt explicit formula. We show how including more zeros
improves the approximation to π(x).
"""

import numpy as np
import matplotlib.pyplot as plt
import mpmath
from pathlib import Path

mpmath.mp.dps = 30

OUT_DIR = Path(__file__).parent / 'plots'
OUT_DIR.mkdir(exist_ok=True)


def sieve_primes(n_max):
    """Generate all primes up to n_max using sieve of Eratosthenes."""
    sieve = np.ones(n_max + 1, dtype=bool)
    sieve[0] = sieve[1] = False
    for i in range(2, int(n_max**0.5) + 1):
        if sieve[i]:
            sieve[i*i::i] = False
    return np.where(sieve)[0]


def prime_counting(x_vals, primes):
    """Compute π(x) = #{p ≤ x : p prime}."""
    pi_x = np.zeros_like(x_vals, dtype=float)
    for i, x in enumerate(x_vals):
        pi_x[i] = np.searchsorted(primes, x, side='right')
    return pi_x


def logarithmic_integral(x):
    """Compute li(x) = ∫_0^x dt/log(t), the logarithmic integral."""
    if x <= 1:
        return 0.0
    return float(mpmath.li(x))


def chebyshev_psi(x_val, n_zeros=0, zeros_t=None):
    """
    Compute the Chebyshev ψ(x) function using the explicit formula:

    ψ(x) = x - Σ_ρ x^ρ/ρ - log(2π) - (1/2)log(1-x^{-2})

    where the sum is over the first n_zeros pairs of non-trivial zeros.
    """
    x = mpmath.mpf(x_val)

    # Main term
    result = x

    # Subtract contribution from zeros
    if zeros_t is not None and n_zeros > 0:
        for k in range(min(n_zeros, len(zeros_t))):
            gamma = mpmath.mpf(zeros_t[k])
            rho = mpmath.mpf('0.5') + 1j * gamma
            rho_conj = mpmath.mpf('0.5') - 1j * gamma

            # x^ρ/ρ + x^{ρ̄}/ρ̄ = 2·Re(x^ρ/ρ)
            x_rho = mpmath.power(x, rho)
            term = x_rho / rho
            result -= 2 * mpmath.re(term)

    # Constant term
    result -= mpmath.log(2 * mpmath.pi)

    # Trivial zeros contribution: -(1/2)log(1-x^{-2})
    if x > 1:
        result -= 0.5 * mpmath.log(1 - x**(-2))

    return float(result)


def explicit_formula_pi(x_val, n_zeros=0, zeros_t=None):
    """
    Compute π(x) using the Riemann explicit formula:

    π(x) ≈ R(x) - Σ_ρ R(x^ρ)

    where R(x) = Σ_{n=1}^∞ μ(n)/n · li(x^{1/n}) is the Riemann R function.
    We use a truncated version.
    """
    x = mpmath.mpf(x_val)
    if x <= 2:
        return 0.0

    # Riemann R function (truncated Gram series)
    def R(y):
        if y <= 1:
            return mpmath.mpf(0)
        result = mpmath.mpf(0)
        log_y = mpmath.log(y)
        term = mpmath.mpf(1)
        for n in range(1, 100):
            term *= log_y / n
            result += term / (n * mpmath.zeta(n + 1))
            if abs(term) < mpmath.mpf('1e-20'):
                break
        return 1 + result

    # Main term: R(x)
    result = float(R(x))

    # Subtract contribution from zeros
    if zeros_t is not None and n_zeros > 0:
        for k in range(min(n_zeros, len(zeros_t))):
            gamma = zeros_t[k]
            rho = 0.5 + 1j * gamma

            # R(x^ρ) — using the approximation R(y) ≈ li(y) for large y
            x_rho = complex(mpmath.power(x, mpmath.mpc(rho)))
            # li(x^ρ) approximation
            li_val = complex(mpmath.li(x_rho))

            # Add contribution from ρ and ρ̄
            result -= 2 * li_val.real

    return result


def load_or_compute_zeros(n_zeros=200):
    """Load precomputed zeros or compute them."""
    data_file = Path(__file__).parent / 'zeros_data.txt'
    zeros_t = []

    if data_file.exists():
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
        for k in range(len(zeros_t) + 1, n_zeros + 1):
            z = mpmath.zetazero(k)
            zeros_t.append(float(z.imag))
            if k % 50 == 0:
                print(f"  {k}/{n_zeros}")

    return zeros_t


def plot_explicit_formula_convergence(x_max=200, n_zeros_list=[0, 5, 20, 50, 100]):
    """Show how the explicit formula converges as more zeros are included."""
    print("Computing explicit formula with different numbers of zeros...")

    fig, axes = plt.subplots(len(n_zeros_list), 1, figsize=(14, 4 * len(n_zeros_list)),
                             sharex=True)

    # Compute actual π(x)
    primes = sieve_primes(x_max + 100)
    x_vals = np.linspace(2, x_max, 1000)
    pi_actual = prime_counting(x_vals, primes)

    # Load zeros
    max_zeros = max(n_zeros_list)
    zeros_t = load_or_compute_zeros(max_zeros)

    for idx, n_z in enumerate(n_zeros_list):
        ax = axes[idx]
        print(f"  Computing with {n_z} zeros...")

        # Compute the approximation using explicit formula
        pi_approx = np.zeros_like(x_vals)
        for i, x in enumerate(x_vals):
            if n_z == 0:
                # Just li(x)
                pi_approx[i] = logarithmic_integral(x)
            else:
                pi_approx[i] = explicit_formula_pi(x, n_z, zeros_t)

        # Plot
        ax.step(x_vals, pi_actual, 'k-', linewidth=1, label='π(x) actual', alpha=0.7)
        ax.plot(x_vals, pi_approx, 'r-', linewidth=1.5,
                label=f'Explicit formula ({n_z} zeros)' if n_z > 0 else 'li(x)')

        error = np.mean(np.abs(pi_approx - pi_actual))
        ax.set_ylabel('Count', fontsize=11)
        ax.set_title(f'{"li(x)" if n_z == 0 else f"{n_z} zeros"}: '
                     f'Mean |error| = {error:.2f}', fontsize=12)
        ax.legend(fontsize=10, loc='upper left')
        ax.grid(True, alpha=0.3)

    axes[-1].set_xlabel('x', fontsize=13)
    plt.suptitle("Riemann's Explicit Formula: Convergence with More Zeros",
                 fontsize=15, y=1.01)
    plt.tight_layout()
    plt.savefig(OUT_DIR / 'explicit_formula_convergence.png', dpi=150, bbox_inches='tight')
    print(f"  Saved: {OUT_DIR / 'explicit_formula_convergence.png'}")
    plt.close()


def plot_chebyshev_explicit(x_max=200, n_zeros_list=[0, 10, 50]):
    """Plot the Chebyshev ψ(x) using the explicit formula."""
    print("\nComputing Chebyshev ψ(x) explicit formula...")

    fig, axes = plt.subplots(1, len(n_zeros_list), figsize=(6 * len(n_zeros_list), 6),
                             sharey=True)

    zeros_t = load_or_compute_zeros(max(n_zeros_list))
    x_vals = np.linspace(2, x_max, 400)

    # Actual ψ(x) = Σ_{p^k ≤ x} log(p)
    primes = sieve_primes(x_max + 100)
    psi_actual = np.zeros_like(x_vals)
    for i, x in enumerate(x_vals):
        total = 0
        for p in primes:
            if p > x:
                break
            pk = p
            while pk <= x:
                total += np.log(p)
                pk *= p
        psi_actual[i] = total

    for idx, n_z in enumerate(n_zeros_list):
        ax = axes[idx]
        print(f"  ψ(x) with {n_z} zeros...")

        psi_approx = np.array([chebyshev_psi(x, n_z, zeros_t) for x in x_vals])

        ax.plot(x_vals, psi_actual, 'k-', linewidth=1, label='ψ(x) actual', alpha=0.7)
        ax.plot(x_vals, psi_approx, 'r-', linewidth=1.5,
                label=f'{n_z} zeros' if n_z > 0 else 'Main term only')

        ax.set_xlabel('x', fontsize=12)
        if idx == 0:
            ax.set_ylabel('ψ(x)', fontsize=12)
        ax.set_title(f'{"Main term x" if n_z == 0 else f"{n_z} zeros"}', fontsize=13)
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)

    plt.suptitle("Chebyshev ψ(x): Explicit Formula Approximation", fontsize=15, y=1.02)
    plt.tight_layout()
    plt.savefig(OUT_DIR / 'chebyshev_explicit.png', dpi=150, bbox_inches='tight')
    print(f"  Saved: {OUT_DIR / 'chebyshev_explicit.png'}")
    plt.close()


def plot_error_vs_zeros(x_test=100):
    """Plot how the approximation error decreases with number of zeros."""
    print(f"\nComputing error as function of number of zeros at x={x_test}...")

    zeros_t = load_or_compute_zeros(200)
    primes = sieve_primes(int(x_test) + 100)
    pi_actual = np.searchsorted(primes, x_test, side='right')

    n_zeros_range = list(range(0, 201, 5))
    errors = []

    for n_z in n_zeros_range:
        if n_z == 0:
            approx = logarithmic_integral(x_test)
        else:
            approx = explicit_formula_pi(x_test, n_z, zeros_t)
        errors.append(abs(approx - pi_actual))

    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    ax.plot(n_zeros_range, errors, 'bo-', markersize=3, linewidth=1)
    ax.set_xlabel('Number of zeros used', fontsize=13)
    ax.set_ylabel(f'|π({x_test}) - approximation|', fontsize=13)
    ax.set_title(f'Approximation Error vs Number of Zeros at x = {x_test}', fontsize=14)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(OUT_DIR / 'error_vs_zeros.png', dpi=150)
    print(f"  Saved: {OUT_DIR / 'error_vs_zeros.png'}")
    plt.close()


def plot_oscillatory_terms(x_max=200, n_zeros=20):
    """Visualize the oscillatory correction from each zero."""
    print(f"\nVisualizing oscillatory terms from first {n_zeros} zeros...")

    zeros_t = load_or_compute_zeros(n_zeros)
    x_vals = np.linspace(2, x_max, 500)

    fig, axes = plt.subplots(2, 1, figsize=(14, 10))

    # Individual oscillatory contributions
    for k in range(min(5, n_zeros)):
        gamma = zeros_t[k]
        contribution = np.zeros_like(x_vals)
        for i, x in enumerate(x_vals):
            rho = 0.5 + 1j * gamma
            x_rho = complex(x**rho)
            contribution[i] = -2 * (x_rho / rho).real

        axes[0].plot(x_vals, contribution, linewidth=1, alpha=0.7,
                    label=f'γ_{k+1} = {gamma:.2f}')

    axes[0].set_xlabel('x', fontsize=12)
    axes[0].set_ylabel('Oscillatory correction', fontsize=12)
    axes[0].set_title('Individual Zero Contributions to ψ(x)', fontsize=14)
    axes[0].legend(fontsize=9)
    axes[0].grid(True, alpha=0.3)

    # Cumulative oscillatory terms
    for n_z in [1, 5, 10, 20]:
        if n_z > n_zeros:
            break
        correction = np.zeros_like(x_vals)
        for k in range(n_z):
            gamma = zeros_t[k]
            for i, x in enumerate(x_vals):
                rho = 0.5 + 1j * gamma
                x_rho = complex(x**rho)
                correction[i] -= 2 * (x_rho / rho).real

        axes[1].plot(x_vals, correction, linewidth=1.5, alpha=0.8,
                    label=f'Sum of {n_z} zero contributions')

    axes[1].set_xlabel('x', fontsize=12)
    axes[1].set_ylabel('Cumulative correction', fontsize=12)
    axes[1].set_title('Cumulative Zero Corrections to ψ(x)', fontsize=14)
    axes[1].legend(fontsize=10)
    axes[1].grid(True, alpha=0.3)

    plt.suptitle('Oscillatory Terms in the Explicit Formula', fontsize=15, y=1.02)
    plt.tight_layout()
    plt.savefig(OUT_DIR / 'oscillatory_terms.png', dpi=150, bbox_inches='tight')
    print(f"  Saved: {OUT_DIR / 'oscillatory_terms.png'}")
    plt.close()


def main():
    print("=" * 60)
    print("RIEMANN'S EXPLICIT FORMULA FOR THE PRIME COUNTING FUNCTION")
    print("=" * 60)
    print()
    print("π(x) = R(x) - Σ_ρ R(x^ρ) + correction terms")
    print("where R(x) is the Riemann R-function and ρ are zeta zeros")
    print()

    # Plot convergence with more zeros
    plot_explicit_formula_convergence(x_max=200, n_zeros_list=[0, 5, 20, 50, 100])

    # Chebyshev ψ(x)
    plot_chebyshev_explicit(x_max=150, n_zeros_list=[0, 10, 50])

    # Error as a function of number of zeros
    plot_error_vs_zeros(x_test=100)

    # Visualize oscillatory terms
    plot_oscillatory_terms(x_max=200, n_zeros=20)

    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print("  The explicit formula demonstrates how prime distribution")
    print("  is encoded in the zeros of ζ(s).")
    print("  Each zero contributes an oscillatory correction.")
    print("  More zeros → better approximation to π(x).")
    print("  The zeros being on Re(s)=1/2 gives the tightest possible")
    print("  error bound: π(x) = li(x) + O(√x log x)")
    print("=" * 60)


if __name__ == '__main__':
    main()
