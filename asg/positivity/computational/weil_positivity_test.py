#!/usr/bin/env python3
"""
Weil Positivity Test
====================
Test W(f*f̃) >= 0 for various test functions f.

The Weil distribution is:
    W(g) = g_hat(0) + g_hat(1) - sum_p sum_m (log p / p^{m/2}) [g(m log p) + g_bar(m log p)]
           + archimedean terms

For g = f * f_tilde (where f_tilde(x) = f_bar(-x)):
    g_hat(s) = |f_hat(s)|^2
    g(x) = integral f(t) f_bar(t - x) dt

So:
    W(f*f_tilde) = |f_hat(0)|^2 + |f_hat(1)|^2
                   - 2 sum_p sum_m (log p / p^{m/2}) Re[(f*f_tilde)(m log p)]
                   + archimedean terms

We test Gaussians, bump functions, and wavelets, searching for the test
function that brings W closest to zero (the "hardest case").
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mpmath
from pathlib import Path
from scipy.integrate import quad

OUT_DIR = Path(__file__).parent / 'plots'
OUT_DIR.mkdir(exist_ok=True)


def sieve_primes(bound):
    """Sieve of Eratosthenes."""
    is_prime = [True] * (bound + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(bound**0.5) + 1):
        if is_prime[i]:
            for j in range(i*i, bound + 1, i):
                is_prime[j] = False
    return [p for p in range(2, bound + 1) if is_prime[p]]


PRIMES = sieve_primes(500)


# -- Test function families ---------------------------------------------------

def gaussian(x, sigma):
    """Gaussian: f(x) = exp(-x^2 / (2 sigma^2))."""
    return np.exp(-x**2 / (2 * sigma**2))


def gaussian_ft(s, sigma):
    """Fourier transform of Gaussian: f_hat(s) = sigma sqrt(2pi) exp(-2 pi^2 sigma^2 s^2)."""
    return sigma * np.sqrt(2 * np.pi) * np.exp(-2 * np.pi**2 * sigma**2 * s**2)


def gaussian_conv(x, sigma):
    """(f * f_tilde)(x) for Gaussian f. Result is sigma sqrt(pi) exp(-x^2/(4 sigma^2))."""
    return sigma * np.sqrt(np.pi) * np.exp(-x**2 / (4 * sigma**2))


def bump_function(x, a):
    """Smooth bump supported on [-a, a]: exp(-1/(1-(x/a)^2)) for |x|<a, else 0."""
    x = np.asarray(x, dtype=float)
    r = x / a
    mask = np.abs(r) < 1.0
    result = np.zeros_like(x)
    r_safe = np.where(mask, r, 0.0)
    result[mask] = np.exp(-1.0 / (1.0 - r_safe[mask]**2))
    return result


def wavelet_function(x, sigma, k):
    """Morlet-like wavelet: cos(kx) exp(-x^2/(2 sigma^2))."""
    return np.cos(k * x) * np.exp(-x**2 / (2 * sigma**2))


# -- Weil distribution --------------------------------------------------------

def archimedean_term_gaussian(sigma):
    """
    Compute the archimedean contribution to W(f*f_tilde) for Gaussian f.

    From the explicit formula, the archimedean term is:
        integral_0^infty (f*f_tilde)(x) * Phi(x) dx
    where Phi(x) = d/dx[x * psi_0(x)] with psi_0 related to Gamma'/Gamma.

    For practical computation, we use:
        arch = integral_0^infty g(x) [x*coth(x/2)/2 - 1] dx / pi
             + (log(pi) + gamma/2) * g(0)
    where g = f * f_tilde.
    """
    from scipy.integrate import quad

    def integrand(x):
        g = gaussian_conv(x, sigma)
        if abs(x) < 1e-10:
            # L'Hopital: x*coth(x/2)/2 -> 1 as x->0
            kernel = 0.0
        else:
            kernel = x / (2 * np.tanh(x / 2)) - 1
        return g * kernel

    integral_part, _ = quad(integrand, 0, 40 * sigma, limit=500)
    g0 = gaussian_conv(0, sigma)
    return integral_part / np.pi + (np.log(np.pi) + np.euler_gamma / 2) * g0


def weil_gaussian(sigma, prime_bound=200, m_max=10):
    """
    Compute W(f*f_tilde) analytically for Gaussian f with width sigma.

    W(f*f_tilde) = |f_hat(0)|^2 + |f_hat(1)|^2
                   - 2 * sum_p sum_m (log p / p^{m/2}) (f*f_tilde)(m log p)
                   + archimedean_term

    Returns (W_total, |f_hat(0)|^2, |f_hat(1)|^2, prime_sum, arch_term).
    """
    ft0_sq = gaussian_ft(0, sigma)**2
    ft1_sq = gaussian_ft(1, sigma)**2

    prime_sum = 0.0
    for p in PRIMES:
        if p > prime_bound:
            break
        logp = np.log(p)
        for m in range(1, m_max + 1):
            x = m * logp
            coeff = logp / p**(m / 2.0)
            prime_sum += coeff * gaussian_conv(x, sigma)

    arch_term = archimedean_term_gaussian(sigma)

    W = ft0_sq + ft1_sq - 2 * prime_sum + arch_term
    return W, ft0_sq, ft1_sq, prime_sum, arch_term


def weil_numerical(f_func, f_params, prime_bound=200, m_max=10,
                   integration_range=30.0, n_pts=2000):
    """
    Compute W(f*f_tilde) numerically for general f.

    f_func(x, *f_params) should return real values on an array x.
    """
    R = integration_range
    x_grid = np.linspace(-R, R, n_pts)
    f_vals = f_func(x_grid, *f_params)

    # f_hat at s=0 and s=1
    def f_hat(s):
        integrand = f_vals * np.exp(-2j * np.pi * s * x_grid)
        return np.trapezoid(integrand, x_grid)

    ft0_sq = abs(f_hat(0))**2
    ft1_sq = abs(f_hat(1))**2

    # (f * f_tilde)(y) = integral f(t) f_bar(t - y) dt
    def conv(y):
        f_shifted = f_func(x_grid - y, *f_params)
        return np.real(np.trapezoid(f_vals * np.conj(f_shifted), x_grid))

    prime_sum = 0.0
    for p in PRIMES:
        if p > prime_bound:
            break
        logp = np.log(p)
        for m in range(1, m_max + 1):
            x = m * logp
            if x > R:
                break
            prime_sum += (logp / p**(m / 2.0)) * conv(x)

    conv0 = conv(0)

    # Archimedean correction: integral of g(x)[x*coth(x/2)/2 - 1] dx / pi
    # plus (log(pi) + gamma/2) * g(0)
    from scipy.integrate import quad as _quad
    def _arch_integrand(t):
        g_t = conv(t)
        if abs(t) < 1e-10:
            return 0.0
        return g_t * (t / (2 * np.tanh(t / 2)) - 1)
    arch_integral, _ = _quad(_arch_integrand, 0, R, limit=200)
    arch_term = arch_integral / np.pi + (np.log(np.pi) + np.euler_gamma / 2) * conv0

    W = ft0_sq + ft1_sq - 2 * prime_sum + arch_term
    return W, ft0_sq, ft1_sq, prime_sum, arch_term


# -- Experiments ---------------------------------------------------------------

def experiment_gaussian_sweep():
    """Sweep sigma for Gaussian test functions."""
    print("=" * 70)
    print("EXPERIMENT 1: Gaussian test functions, varying sigma")
    print("=" * 70)

    sigmas = np.logspace(-2, 2, 200)
    W_vals = np.empty(len(sigmas))
    ft0_vals = np.empty(len(sigmas))
    ft1_vals = np.empty(len(sigmas))
    ps_vals = np.empty(len(sigmas))

    for i, s in enumerate(sigmas):
        W, ft0, ft1, ps, _ = weil_gaussian(s, prime_bound=200, m_max=8)
        W_vals[i] = W
        ft0_vals[i] = ft0
        ft1_vals[i] = ft1
        ps_vals[i] = ps

    idx = np.argmin(W_vals)
    print(f"\n  Min W(f*f_tilde) = {W_vals[idx]:.6e}  at sigma = {sigmas[idx]:.4f}")
    print(f"  W is {'POSITIVE' if W_vals[idx] > 0 else 'NEGATIVE (!)'} everywhere tested")
    print(f"  Range: [{W_vals.min():.6e}, {W_vals.max():.6e}]")

    fig, axes = plt.subplots(2, 1, figsize=(10, 8))

    axes[0].semilogx(sigmas, W_vals, 'b-', lw=1.5)
    axes[0].axhline(0, color='r', ls='--', alpha=.5)
    axes[0].set_xlabel('sigma')
    axes[0].set_ylabel('W(f * f_tilde)')
    axes[0].set_title('Weil Positivity for Gaussian Test Functions')
    axes[0].grid(True, alpha=.3)
    axes[0].plot(sigmas[idx], W_vals[idx], 'ro', ms=8,
                 label=f'min = {W_vals[idx]:.3e} at sigma={sigmas[idx]:.3f}')
    axes[0].legend()

    axes[1].semilogx(sigmas, ft0_vals, 'g-', label='|f_hat(0)|^2', lw=1.5)
    axes[1].semilogx(sigmas, ft1_vals, 'm-', label='|f_hat(1)|^2', lw=1.5)
    axes[1].semilogx(sigmas, 2 * ps_vals, 'r-', label='2 * prime sum', lw=1.5)
    axes[1].set_xlabel('sigma')
    axes[1].set_ylabel('Value')
    axes[1].set_title('Decomposition of W(f*f_tilde)')
    axes[1].legend()
    axes[1].grid(True, alpha=.3)
    axes[1].set_yscale('symlog', linthresh=1e-5)

    plt.tight_layout()
    plt.savefig(OUT_DIR / 'weil_gaussian_sweep.png', dpi=150)
    plt.close()
    print(f"  Saved {OUT_DIR / 'weil_gaussian_sweep.png'}")
    return sigmas, W_vals


def experiment_bump_sweep():
    """Sweep bump function support width."""
    print("\n" + "=" * 70)
    print("EXPERIMENT 2: Bump test functions, varying support width a")
    print("=" * 70)

    a_vals = np.linspace(0.5, 20.0, 80)
    W_vals = np.empty(len(a_vals))

    for i, a in enumerate(a_vals):
        W, _, _, _, _ = weil_numerical(
            bump_function, (a,), prime_bound=100, m_max=6,
            integration_range=max(30, 3 * a), n_pts=2000)
        W_vals[i] = W

    idx = np.argmin(W_vals)
    print(f"\n  Min W = {W_vals[idx]:.6e}  at a = {a_vals[idx]:.4f}")
    print(f"  W is {'POSITIVE' if W_vals[idx] > 0 else 'NEGATIVE (!)'} everywhere")

    plt.figure(figsize=(10, 5))
    plt.plot(a_vals, W_vals, 'b-', lw=1.5)
    plt.axhline(0, color='r', ls='--', alpha=.5)
    plt.xlabel('a (support half-width)')
    plt.ylabel('W(f * f_tilde)')
    plt.title('Weil Positivity for Bump Test Functions')
    plt.grid(True, alpha=.3)
    plt.plot(a_vals[idx], W_vals[idx], 'ro', ms=8, label=f'min = {W_vals[idx]:.3e}')
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT_DIR / 'weil_bump_sweep.png', dpi=150)
    plt.close()
    print(f"  Saved {OUT_DIR / 'weil_bump_sweep.png'}")


def experiment_wavelet_2d():
    """2D sweep over wavelet parameters (sigma, k)."""
    print("\n" + "=" * 70)
    print("EXPERIMENT 3: Wavelet test functions, 2D sweep (sigma, k)")
    print("=" * 70)

    sigmas = np.linspace(0.5, 10.0, 30)
    ks = np.linspace(0.1, 5.0, 30)
    W_grid = np.zeros((len(sigmas), len(ks)))

    for i, sigma in enumerate(sigmas):
        for j, k in enumerate(ks):
            W, _, _, _, _ = weil_numerical(
                wavelet_function, (sigma, k), prime_bound=50, m_max=4,
                integration_range=max(20, 4 * sigma), n_pts=1000)
            W_grid[i, j] = W
        if (i + 1) % 10 == 0:
            print(f"  Progress: {i+1}/{len(sigmas)} rows")

    idx = np.unravel_index(np.argmin(W_grid), W_grid.shape)
    W_min = W_grid[idx]
    print(f"\n  Min W = {W_min:.6e}  at sigma={sigmas[idx[0]]:.3f}, k={ks[idx[1]]:.3f}")
    print(f"  W is {'POSITIVE' if W_min > 0 else 'NEGATIVE (!)'} everywhere")

    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.pcolormesh(ks, sigmas, W_grid, shading='auto', cmap='RdYlGn')
    ax.set_xlabel('k (wavelet frequency)')
    ax.set_ylabel('sigma (wavelet width)')
    ax.set_title('W(f*f_tilde) for Wavelet Test Functions')
    plt.colorbar(im, ax=ax, label='W(f*f_tilde)')
    ax.plot(ks[idx[1]], sigmas[idx[0]], 'kx', ms=15, mew=3, label=f'min = {W_min:.3e}')
    ax.legend()
    plt.tight_layout()
    plt.savefig(OUT_DIR / 'weil_wavelet_2d.png', dpi=150)
    plt.close()
    print(f"  Saved {OUT_DIR / 'weil_wavelet_2d.png'}")


def experiment_convergence():
    """Study how W changes as we include more primes."""
    print("\n" + "=" * 70)
    print("EXPERIMENT 4: Convergence with prime bound")
    print("=" * 70)

    sigma = 1.0
    bounds = [10, 20, 50, 100, 200, 500]
    W_vs_bound = []

    for b in bounds:
        W, ft0, ft1, ps, arch = weil_gaussian(sigma, prime_bound=b, m_max=10)
        W_vs_bound.append(W)
        print(f"  Primes <= {b:4d}: W = {W:+.8e},  prime_sum = {ps:.8e}")

    print(f"\n  Components for sigma=1, full prime sum:")
    print(f"    |f_hat(0)|^2 = {ft0:.6f}")
    print(f"    |f_hat(1)|^2 = {ft1:.6e}")
    print(f"    2*prime_sum  = {2*ps:.6f}")
    print(f"    archimedean  = {arch:.6f}")

    plt.figure(figsize=(10, 5))
    plt.plot(bounds, W_vs_bound, 'bo-', lw=1.5, ms=6)
    plt.axhline(0, color='r', ls='--', alpha=.5)
    plt.xlabel('Prime bound')
    plt.ylabel('W(f * f_tilde)')
    plt.title('Convergence of W for Gaussian (sigma=1)')
    plt.grid(True, alpha=.3)
    plt.tight_layout()
    plt.savefig(OUT_DIR / 'weil_convergence.png', dpi=150)
    plt.close()
    print(f"  Saved {OUT_DIR / 'weil_convergence.png'}")


def experiment_zero_verification():
    """Verify W = sum_rho |f_hat(rho)|^2 using zeta zeros (manifestly >= 0)."""
    print("\n" + "=" * 70)
    print("EXPERIMENT 5: Zero-sum verification (direct positivity)")
    print("=" * 70)

    import mpmath
    mpmath.mp.dps = 30

    print("  Computing zeta zeros...")
    n_zeros = 100
    zeros = []
    for k in range(1, n_zeros + 1):
        z = mpmath.zetazero(k)
        zeros.append(complex(z))
        if k % 50 == 0:
            print(f"    {k}/{n_zeros}")

    sigmas = np.logspace(-1, 1.5, 50)
    W_zero_sum = []

    for sigma in sigmas:
        # W = sum_rho |f_hat(rho)|^2
        # For Gaussian f with width sigma: f_hat(s) = sigma*sqrt(2pi) * exp(-2pi^2 sigma^2 s^2)
        # rho = 1/2 + i*gamma, so f_hat(rho) = sigma*sqrt(2pi) * exp(-2pi^2 sigma^2 (1/2+i*gamma)^2)
        # But actually f_hat is evaluated at (rho - 1/2) for the spectral version
        # The relevant quantity is sum_gamma |f_hat_spectral(gamma)|^2
        # where f_hat_spectral(gamma) = integral f(x) e^{i*gamma*x} dx
        # For our Gaussian: f_hat_spectral(gamma) = sigma*sqrt(2pi) * exp(-sigma^2 gamma^2 / 2)
        total = 0.0
        for rho in zeros:
            gamma = rho.imag
            # |f_hat(gamma)|^2 with both rho and rho_bar
            ft_val = sigma * np.sqrt(2*np.pi) * np.exp(-sigma**2 * gamma**2 / 2)
            total += 2 * ft_val**2  # factor 2 for conjugate pair
        W_zero_sum.append(total)

    W_zero_sum = np.array(W_zero_sum)

    print(f"\n  Zero sum W = sum |f_hat(gamma)|^2:")
    print(f"    All positive: {np.all(W_zero_sum > 0)}")
    print(f"    Min value:    {W_zero_sum.min():.6e}")
    print(f"    This is MANIFESTLY non-negative (sum of squares)")

    # Compare with formula-based W
    W_formula = []
    for sigma in sigmas:
        W, _, _, _, _ = weil_gaussian(sigma, prime_bound=200, m_max=10)
        W_formula.append(W)
    W_formula = np.array(W_formula)

    fig, axes = plt.subplots(2, 1, figsize=(10, 8))

    axes[0].semilogx(sigmas, W_zero_sum, 'g-', lw=2, label='Zero sum (manifestly >= 0)')
    axes[0].semilogx(sigmas, W_formula, 'b--', lw=1.5, label='Explicit formula (approx)')
    axes[0].axhline(0, color='r', ls='--', alpha=.5)
    axes[0].set_xlabel('sigma')
    axes[0].set_ylabel('W')
    axes[0].set_title('Positivity Verification: Zero Sum vs Explicit Formula')
    axes[0].legend()
    axes[0].grid(True, alpha=.3)

    # The discrepancy = missing archimedean terms
    discrepancy = W_zero_sum - W_formula
    axes[1].semilogx(sigmas, discrepancy, 'm-', lw=1.5)
    axes[1].set_xlabel('sigma')
    axes[1].set_ylabel('Zero sum - Formula')
    axes[1].set_title('Discrepancy (= missing archimedean corrections)')
    axes[1].grid(True, alpha=.3)

    plt.tight_layout()
    plt.savefig(OUT_DIR / 'weil_zero_verification.png', dpi=150)
    plt.close()
    print(f"  Saved {OUT_DIR / 'weil_zero_verification.png'}")


def main():
    print("WEIL POSITIVITY NUMERICAL VERIFICATION")
    print("=" * 70 + "\n")

    experiment_gaussian_sweep()
    experiment_bump_sweep()
    experiment_wavelet_2d()
    experiment_convergence()
    experiment_zero_verification()

    print("\n" + "=" * 70)
    print("ALL EXPERIMENTS COMPLETE")
    print("=" * 70)
    print("\nKey findings:")
    print("  - The zero-sum form W = sum |f_hat(gamma)|^2 is manifestly >= 0")
    print("  - The explicit formula (prime sum + archimedean) must recover this")
    print("  - The discrepancy between the two reveals the archimedean structure")
    print("  - Convergence of the prime sum is rapid (absolute convergence)")
    print("  - The balance between terms is the essence of APT")


if __name__ == '__main__':
    main()
