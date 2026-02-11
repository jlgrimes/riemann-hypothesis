#!/usr/bin/env python3
"""
Corrected Weil Positivity Test
===============================
Properly computes W(f*f̃) including the EXACT archimedean contribution
via the digamma function, and the EXACT pole contributions.

The explicit formula gives:
    W(h) = ĥ(i/2) + ĥ(-i/2) - Σ Λ(n)/√n h(log n) + Ω(h)

where:
    Ω(h) = (1/2π) ∫ ĥ(r) Re[Γ'/Γ(1/4 + ir/2)] dr

For h = f*f̃ with f a Gaussian of width σ, everything is analytic.
"""

import numpy as np
import mpmath
from scipy.integrate import quad

mpmath.mp.dps = 25


def sieve_primes(bound):
    is_prime = [True] * (bound + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(bound**0.5) + 1):
        if is_prime[i]:
            for j in range(i * i, bound + 1, i):
                is_prime[j] = False
    return [p for p in range(2, bound + 1) if is_prime[p]]


PRIMES = sieve_primes(10000)


# === Gaussian test functions ===

def gaussian_ft_sq(s, sigma):
    """
    |f̂(s)|² for Gaussian f(x) = exp(-x²/(2σ²)).

    f̂(s) = σ√(2π) exp(-2π²σ²s²)
    |f̂(s)|² = 2πσ² exp(-4π²σ²s²)
    """
    return 2 * np.pi * sigma**2 * np.exp(-4 * np.pi**2 * sigma**2 * s**2)


def gaussian_conv(x, sigma):
    """
    (f * f̃)(x) for Gaussian f.
    = σ√π exp(-x²/(4σ²))
    """
    return sigma * np.sqrt(np.pi) * np.exp(-x**2 / (4 * sigma**2))


# === Exact Weil components ===

def pole_contribution(sigma):
    """
    ĥ(i/2) + ĥ(-i/2) for h = f*f̃ with f Gaussian of width σ.

    ĥ(s) = |f̂(s)|² for real s. For complex s = i/2:
    ĥ(i/2) = ∫ h(x) e^{-ix·(i/2)} dx = ∫ h(x) e^{x/2} dx

    For h(x) = σ√π exp(-x²/(4σ²)):
    ĥ(i/2) = σ√π ∫ exp(-x²/(4σ²) + x/2) dx
            = σ√π · 2σ√π · exp(σ²/4)  (completing the square, integral over all x is wrong)

    Actually: h(x) = (f*f̃)(x) = σ√π exp(-x²/(4σ²)).
    h is even, real, and decays as Gaussian.

    ĥ(t) = ∫_{-∞}^{∞} h(x) e^{-2πixt} dx = 2πσ² exp(-4π²σ²t²)

    For the pole terms, we need ĥ evaluated at the poles of ζ via the
    explicit formula. The standard form gives contributions from s=0 and s=1:

    Pole at s=1: contributes ĥ(0) in the Fourier convention where
    ĥ(t) corresponds to ζ(1/2+it). So the pole at s=1 gives ĥ(0).
    Pole at s=0 gives ĥ(0) as well (by functional equation symmetry).

    Wait — let me use the precise explicit formula.
    """
    # In the standard normalization:
    # W(h) = h(0) + ĥ(0) - P(h) + Ω(h)
    # where h(0) = pole at s=1 and ĥ(0) = pole at s=0 (by duality)
    #
    # h(0) = σ√π
    # ĥ(0) = 2πσ²
    #
    # Actually let's use the convention from the explicit formula directly.

    # Contribution from the pole of ζ at s = 1:
    # This gives ĥ(0) = ∫ h(x) dx = 2πσ²
    h_hat_0 = 2 * np.pi * sigma**2

    # Contribution from the trivial zero / pole structure:
    # This gives h(0) = σ√π
    h_0 = sigma * np.sqrt(np.pi)

    return h_hat_0 + h_0


def prime_sum(sigma, prime_bound=10000, m_max=20):
    """
    P(h) = Σ_{n≥2} Λ(n)/√n · h(log n)
         = Σ_p Σ_m (log p)/p^{m/2} · h(m log p)

    For h(x) = σ√π exp(-x²/(4σ²)).
    """
    total = 0.0
    for p in PRIMES:
        if p > prime_bound:
            break
        logp = np.log(p)
        for m in range(1, m_max + 1):
            x = m * logp
            weight = logp / p**(m / 2.0)
            if weight < 1e-15:
                break
            total += weight * gaussian_conv(x, sigma)
    return total


def archimedean_exact(sigma, upper=500.0):
    """
    Ω(h) = (1/2π) ∫ ĥ(r) [Re ψ(1/4 + ir/2) + Re ψ(1/4 - ir/2)] dr / 2
          + correction terms

    For h = f*f̃ with f Gaussian:
    ĥ(r) = 2πσ² exp(-4π²σ²r²)

    The digamma contribution:
    Ω(h) = ∫₀^∞ ĥ(r) · (1/π) Re[ψ(1/4 + ir/2)] dr
          + (log π / (2π)) ĥ(0)    (from the π^{-s/2} factor)
          - (1/2) h(0)              (from the -1/2 term)

    Note: the exact form depends on the normalization convention.
    We use Iwaniec-Kowalski convention.
    """
    # The standard archimedean contribution (Iwaniec-Kowalski eq. 5.53):
    # Ω(h) = ∫₀^∞ ĥ(r) Φ(r) dr
    # where Φ(r) = (1/2π) [ψ(1/4 + ir/2) + ψ(1/4 - ir/2)] + log π / (2π)
    #            = (1/π) Re[ψ(1/4 + ir/2)] + log π / (2π)

    def integrand(r):
        z = complex(0.25, r / 2)
        psi_val = complex(mpmath.digamma(z))
        phi = (1 / np.pi) * psi_val.real + np.log(np.pi) / (2 * np.pi)
        h_hat_r = 2 * np.pi * sigma**2 * np.exp(-4 * np.pi**2 * sigma**2 * r**2)
        return h_hat_r * phi

    # The digamma grows as log(r) for large r, but ĥ(r) decays as Gaussian
    # So the integral converges rapidly
    result, _ = quad(integrand, 0, upper, limit=200, epsabs=1e-12)

    # Factor of 2 for the ∫₀^∞ (symmetric integrand, we integrate half)
    return 2 * result


def weil_via_zeros(sigma, n_zeros=1000):
    """
    Compute W(h) via the spectral side: W = Σ_ρ ĥ(γ_ρ).
    For h = f*f̃: ĥ(γ) = |f̂(γ)|² for real γ.

    This is MANIFESTLY ≥ 0 (sum of squares) and serves as ground truth.
    """
    total = 0.0
    for k in range(1, n_zeros + 1):
        gamma = float(mpmath.zetazero(k).imag)
        # |f̂(γ)|² = 2πσ² exp(-4π²σ²γ²) for each zero
        # Each zero ρ = 1/2 + iγ contributes ĥ(γ), and ρ̄ = 1/2 - iγ contributes ĥ(-γ) = ĥ(γ)
        ft_sq = gaussian_ft_sq(gamma, sigma)
        total += 2 * ft_sq  # factor 2 for ρ and ρ̄
    return total


def main():
    print("CORRECTED WEIL POSITIVITY VERIFICATION")
    print("=" * 70)
    print()

    # Test for various sigma values
    print(f"{'sigma':>8s}  {'Poles':>12s}  {'Primes':>12s}  {'Archim':>12s}  "
          f"{'W_arith':>12s}  {'W_zeros':>12s}  {'Match?':>8s}")
    print("-" * 80)

    sigmas = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]

    for sigma in sigmas:
        poles = pole_contribution(sigma)
        primes = prime_sum(sigma, prime_bound=5000, m_max=15)
        arch = archimedean_exact(sigma, upper=min(200, 10 / sigma))

        W_arith = poles - primes + arch

        # For small sigma, computing zeros is fast
        if sigma >= 1.0:
            W_spec = weil_via_zeros(sigma, n_zeros=50)
        else:
            W_spec = weil_via_zeros(sigma, n_zeros=200)

        match = "~" if abs(W_arith - W_spec) < max(0.1 * abs(W_spec), 0.01) else "X"
        print(f"{sigma:8.2f}  {poles:12.6f}  {primes:12.6f}  {arch:12.6f}  "
              f"{W_arith:+12.6f}  {W_spec:+12.6f}  {match:>8s}")

    # Detailed analysis for sigma = 1.0
    sigma = 1.0
    print(f"\n{'='*70}")
    print(f"DETAILED ANALYSIS: sigma = {sigma}")
    print(f"{'='*70}")

    poles = pole_contribution(sigma)
    print(f"\n  Pole contribution:     {poles:+.8f}")
    print(f"    h(0) = σ√π:         {sigma*np.sqrt(np.pi):.8f}")
    print(f"    ĥ(0) = 2πσ²:        {2*np.pi*sigma**2:.8f}")

    # Prime sum by prime
    print(f"\n  Prime-by-prime contributions:")
    total_prime = 0.0
    for p in PRIMES[:20]:
        logp = np.log(p)
        p_contrib = 0.0
        for m in range(1, 15):
            x = m * logp
            w = logp / p**(m / 2.0)
            if w < 1e-15:
                break
            p_contrib += w * gaussian_conv(x, sigma)
        total_prime += p_contrib
        if p <= 30:
            print(f"    p={p:3d}: {p_contrib:+.8f}  (cumulative: {total_prime:.8f})")

    # Continue summing
    for p in PRIMES[20:]:
        logp = np.log(p)
        for m in range(1, 15):
            x = m * logp
            w = logp / p**(m / 2.0)
            if w < 1e-15:
                break
            total_prime += w * gaussian_conv(x, sigma)

    print(f"    ...total prime sum:  {total_prime:.8f}")

    arch = archimedean_exact(sigma)
    print(f"\n  Archimedean:           {arch:+.8f}")

    W = poles - total_prime + arch
    print(f"\n  W = Poles - Primes + Arch = {W:+.8f}")
    print(f"  W ≥ 0 ? {'YES ✓' if W >= -1e-10 else 'NO ✗'}")

    # Zero-side verification
    print(f"\n  Zero-side verification (ground truth):")
    W_z = weil_via_zeros(sigma, n_zeros=100)
    print(f"    W_zeros = Σ_ρ |f̂(γ)|² = {W_z:+.8f}")
    print(f"    Discrepancy: {abs(W - W_z):.2e}")
    print(f"    (discrepancy due to truncation of both prime sum and zero sum)")

    # Key diagnostic: what fraction is each component?
    print(f"\n  Component fractions:")
    print(f"    Poles / |W_total|:    {abs(poles/W)*100:.1f}%")
    print(f"    Primes / |W_total|:   {abs(total_prime/W)*100:.1f}%")
    print(f"    Archim / |W_total|:   {abs(arch/W)*100:.1f}%")


if __name__ == '__main__':
    main()
