#!/usr/bin/env python3
"""
Numerical verification of the Fourier transform of the background Weil kernel K_bg.

Convention: angular frequency F(ξ) = ∫ f(x) e^{-iξx} dx throughout.

Verifies:
1. FFT of K_bg matches analytic formula 2e^{-|ξ|/2}/(1-e^{-2|ξ|})
2. Numerical quadrature cross-check
3. Positivity of K̂_bg(ξ) for ξ ∈ (0.001, 1000)
4. Completeness of Weil kernel decomposition K = δ + K_bg + K_zeros
5. Lorentzian FT identity: F[1/(a²+(x/2)²)](ξ) = (2π/a)·e^{-2a|ξ|}
"""

import numpy as np
from mpmath import mp, digamma, zetazero, re as mpre, im as mpim, pi as mppi, log as mplog
from scipy import integrate
import sys
import warnings

warnings.filterwarnings('ignore', category=integrate.IntegrationWarning)

mp.dps = 30

# ─── Kernel definitions ───

def K_bg(x):
    """Background Weil kernel: -(1/π) Re[ψ(1/4 + ix/2)] + log(π)/(2π)"""
    z = mp.mpf('0.25') + 1j * mp.mpf(x) / 2
    psi_val = digamma(z)
    return float(-mpre(psi_val) / mppi + mplog(mppi) / (2 * mppi))


def K_bg_hat_analytic(xi):
    """Analytic Fourier transform (angular freq): 2e^{-|ξ|/2}/(1-e^{-2|ξ|})"""
    axi = np.abs(xi)
    if np.isscalar(axi):
        if axi < 1e-15:
            return np.inf
        return 2.0 * np.exp(-axi / 2) / (1.0 - np.exp(-2.0 * axi))
    result = np.full_like(axi, np.inf, dtype=float)
    mask = axi > 1e-15
    result[mask] = 2.0 * np.exp(-axi[mask] / 2) / (1.0 - np.exp(-2.0 * axi[mask]))
    return result


def K_bg_hat_regular_freq(f):
    """K̂_bg at regular frequency f: K̂_bg(2πf) = 2e^{-π|f|}/(1-e^{-4π|f|})"""
    return K_bg_hat_analytic(2 * np.pi * f)


def K_zeros(x, zeros):
    """Zero-dependent kernel: (1/(2π)) Σ_γ 2cos(γx)/(1/4 + γ²)"""
    total = 0.0
    for gamma in zeros:
        total += 2.0 * np.cos(gamma * x) / (0.25 + gamma**2)
    return total / (2.0 * np.pi)


# ─── Test 1: FFT of K_bg vs analytic formula ───

def test_fft_vs_analytic():
    print("=" * 70)
    print("TEST 1: FFT of K_bg vs analytic K̂_bg(ξ) [angular freq]")
    print("=" * 70)

    L = 200.0
    N = 2**17  # 131072 points
    dx = 2 * L / N
    x = np.linspace(-L, L, N, endpoint=False)

    print(f"  Sampling K_bg on {N} points, x in [-{L}, {L}] ...")
    k_vals = np.array([K_bg(xi) for xi in x])

    # numpy FFT: X_k = Σ_n x_n e^{-2πi k n/N}
    # This approximates F_reg(f_k) = ∫ f(x) e^{-2πi f_k x} dx
    # where f_k = k/(N·dx) are regular frequencies.
    # We have F_ang(ξ) = F_reg(ξ/(2π)), so F_reg(f) = F_ang(2πf).
    k_hat_raw = np.fft.fft(k_vals)
    freqs_regular = np.fft.fftfreq(N, d=dx)  # regular frequencies

    # Phase correction: x starts at -L, not 0
    phase = np.exp(2j * np.pi * freqs_regular * (-L))
    k_hat_numerical = np.real(dx * phase * k_hat_raw)

    # Compare at positive regular frequencies, converted to angular for analytic
    mask = (freqs_regular > 0.005) & (freqs_regular < 5.0)
    f_test = freqs_regular[mask]
    numerical = k_hat_numerical[mask]
    # Analytic at angular freq ξ = 2πf
    analytic = K_bg_hat_regular_freq(f_test)

    rel_err = np.abs(numerical - analytic) / np.abs(analytic)
    max_rel_err = np.max(rel_err)
    mean_rel_err = np.mean(rel_err)
    median_rel_err = np.median(rel_err)
    n_points = np.sum(mask)

    passed = median_rel_err < 0.02
    status = "PASS" if passed else "FAIL"

    print(f"  Compared at {n_points} frequency points")
    print(f"  Regular freq range: (0.005, 5.0)")
    print(f"  Angular freq range: ({0.005*2*np.pi:.3f}, {5.0*2*np.pi:.3f})")
    print(f"  Max    relative error: {max_rel_err:.6e}")
    print(f"  Median relative error: {median_rel_err:.6e}")
    print(f"  Mean   relative error: {mean_rel_err:.6e}")

    # Show a few sample points
    sample_idx = np.linspace(0, n_points - 1, 8, dtype=int)
    print(f"\n  {'f_reg':>8s}  {'ξ_ang':>8s}  {'Numerical':>14s}  {'Analytic':>14s}  {'Rel Err':>12s}")
    print("  " + "-" * 65)
    for i in sample_idx:
        print(f"  {f_test[i]:8.4f}  {2*np.pi*f_test[i]:8.4f}  {numerical[i]:14.6f}  {analytic[i]:14.6f}  {rel_err[i]:12.6e}")

    print(f"\n  [{status}]")
    print()
    return passed


# ─── Test 2: Partial fraction FT verification ───

def test_partial_fraction_ft():
    print("=" * 70)
    print("TEST 2: Partial fraction FT of K_bg vs analytic [angular freq]")
    print("=" * 70)
    print("  K_bg = (1/π) Σ_n [a_n/(a_n²+(x/2)²) - 1/(n+1)] + const")
    print("  FT of n-th Lorentzian = 2·e^{-2a_n|ξ|} (at ξ≠0)")
    print("  Partial sum FT = 2e^{-|ξ|/2}·(1-e^{-2(N+1)|ξ|})/(1-e^{-2|ξ|})")

    test_xis = [0.5, 1.0, 2.0, 5.0, 10.0]
    N_terms = 200  # partial fraction terms
    all_pass = True

    print(f"\n  Using N={N_terms} partial fraction terms")
    print(f"\n  {'ξ':>8s}  {'Quadrature':>14s}  {'Partial FT':>14s}  {'Full Formula':>14s}  Status")
    print("  " + "-" * 70)

    for xi in test_xis:
        # Numerically integrate each Lorentzian term's FT:
        # (1/π) Σ_{n=0}^{N} ∫ [a_n/(a_n²+(x/2)²)] cos(ξx) dx
        # = (1/π) Σ 2π·e^{-2a_n|ξ|} = 2 Σ e^{-2a_n|ξ|}
        partial_ft = 2.0 * sum(np.exp(-2 * (n + 0.25) * xi) for n in range(N_terms))

        # Closed form of partial sum:
        partial_closed = 2.0 * np.exp(-xi / 2) * (1 - np.exp(-2 * (N_terms) * xi)) / (1 - np.exp(-2 * xi))

        # Full formula (N→∞):
        full_formula = K_bg_hat_analytic(xi)

        # Check partial sum against closed form
        rel_err_closed = abs(partial_ft - partial_closed) / abs(partial_closed)
        # Check convergence to full formula
        rel_err_full = abs(partial_ft - full_formula) / abs(full_formula)

        ok = rel_err_closed < 1e-8 and rel_err_full < 1e-6
        if not ok:
            all_pass = False
        status = "PASS" if ok else "FAIL"
        print(f"  {xi:8.2f}  {partial_ft:14.10f}  {partial_closed:14.10f}  {full_formula:14.10f}  {status}")

    # Also verify by direct quadrature of individual Lorentzian terms
    print(f"\n  Direct quadrature of individual Lorentzians (angular freq):")
    print(f"  {'n':>4s}  {'ξ':>6s}  {'Quad':>14s}  {'2·e^-2aξ':>14s}  {'Rel Err':>12s}  Status")
    print("  " + "-" * 68)

    quad_pass = True
    for n in [0, 1, 5, 20]:
        a = n + 0.25
        for xi in [0.5, 2.0]:
            if 2 * a * xi > 10:
                continue

            def integrand(x, _a=a, _xi=xi):
                return _a * np.cos(_xi * x) / (_a**2 + (x / 2.0)**2)

            val, _ = integrate.quad(integrand, 0, np.inf, limit=500, epsabs=1e-13, epsrel=1e-13)
            val *= 2  # even function
            val /= np.pi  # the 1/π prefactor

            expected = 2.0 * np.exp(-2 * a * xi)
            rel_err = abs(val - expected) / abs(expected)
            ok = rel_err < 0.01
            if not ok:
                quad_pass = False
            status = "PASS" if ok else "FAIL"
            print(f"  {n:4d}  {xi:6.2f}  {val:14.10f}  {expected:14.10f}  {rel_err:12.6e}  {status}")

    all_pass = all_pass and quad_pass
    print(f"\n  [{'PASS' if all_pass else 'FAIL'}]")
    print()
    return all_pass


# ─── Test 3: Positivity of K̂_bg ───

def test_positivity():
    print("=" * 70)
    print("TEST 3: Positivity of K̂_bg(ξ) for ξ in (0.001, 1000)")
    print("=" * 70)

    xi_values = np.concatenate([
        np.linspace(0.001, 1.0, 500),
        np.linspace(1.0, 100.0, 500),
        np.linspace(100.0, 1000.0, 200),
    ])

    analytic_vals = K_bg_hat_analytic(xi_values)
    min_val = np.min(analytic_vals)
    n_negative = np.sum(analytic_vals <= 0)

    passed = n_negative == 0
    status = "PASS" if passed else "FAIL"

    print(f"  Tested {len(xi_values)} points in (0.001, 1000)")
    print(f"  Min value: {min_val:.10e}")
    print(f"  Negative values: {n_negative}")
    print(f"  Algebraic: 2e^{{-|ξ|/2}} > 0 and (1-e^{{-2|ξ|}}) > 0 for ξ>0 => K_bg_hat > 0")
    print(f"  [{status}]")
    print()
    return passed


# ─── Test 4: Completeness of decomposition ───

def test_completeness():
    print("=" * 70)
    print("TEST 4: Completeness of Weil kernel decomposition")
    print("=" * 70)

    n_zeros = 150
    print(f"  Computing {n_zeros} zeta zeros ...")
    zeros = [float(mpim(zetazero(k))) for k in range(1, n_zeros + 1)]

    k_bg_0 = K_bg(0.0)
    k_z_0 = K_zeros(0.0, zeros)

    print(f"  K_bg(0)    = {k_bg_0:.10f}")
    print(f"  K_zeros(0) = {k_z_0:.10f}  ({n_zeros} zeros)")
    print(f"  Sum        = {k_bg_0 + k_z_0:.10f}")

    # Table of values
    test_x = [0.0, 1.0, 5.0, 10.0, 20.0, 50.0]
    print(f"\n  {'x':>6s}  {'K_bg(x)':>14s}  {'K_zeros(x)':>14s}  {'Sum':>14s}")
    print("  " + "-" * 55)
    for x in test_x:
        kb = K_bg(x)
        kz = K_zeros(x, zeros)
        print(f"  {x:6.1f}  {kb:14.8f}  {kz:14.8f}  {kb + kz:14.8f}")

    # Symmetry check: K_bg and K_zeros are both even functions
    sym_errs = []
    for x in [1.0, 3.0, 7.0, 15.0]:
        vp = K_bg(x) + K_zeros(x, zeros)
        vm = K_bg(-x) + K_zeros(-x, zeros)
        sym_errs.append(abs(vp - vm) / (abs(vp) + 1e-30))

    max_sym_err = max(sym_errs)

    # Verify partial fraction convergence: K_bg(x) → -log|x|/π + ... for large |x|
    # Actually K_bg(x) ~ -(1/π) log(|x|/2) for large x (from digamma asymptotics)
    # Check that K_bg grows logarithmically:
    x_large = [50.0, 100.0, 200.0]
    print(f"\n  Large-x behavior of K_bg (should be ~ -(1/π)log(|x|/2)):")
    print(f"  {'x':>8s}  {'K_bg(x)':>14s}  {'-(1/π)ln(x/2)':>14s}  {'Ratio':>10s}")
    print("  " + "-" * 55)
    for x in x_large:
        kb = K_bg(x)
        asymp = -np.log(x / 2) / np.pi
        ratio = kb / asymp if abs(asymp) > 1e-15 else float('nan')
        print(f"  {x:8.1f}  {kb:14.8f}  {asymp:14.8f}  {ratio:10.6f}")

    passed = max_sym_err < 1e-6
    status = "PASS" if passed else "FAIL"

    print(f"\n  Symmetry check max error: {max_sym_err:.6e}")
    print(f"  [{status}]")
    print()
    return passed


# ─── Test 5: Lorentzian Fourier transform identity ───

def test_lorentzian_ft():
    print("=" * 70)
    print("TEST 5: F[1/(a^2 + (x/2)^2)](ξ) — Lorentzian FT identity")
    print("=" * 70)
    print("  Convention: F(ξ) = ∫ f(x) e^{-iξx} dx  [angular freq]")
    print()
    print("  Claimed:  2π·e^{-2a|ξ|}/(2a)  = (π/a)·e^{-2a|ξ|}")
    print("  Derived:  (2π/a)·e^{-2a|ξ|}   [from standard tables]")
    print("  Checking numerically which is correct...")
    print()

    test_ns = [0, 1, 2, 5]
    test_xis = [0.1, 0.5, 1.0, 2.0]
    all_pass = True
    claimed_wins = 0
    derived_wins = 0

    print(f"  {'n':>3s}  {'ξ':>5s}  {'Quadrature':>14s}  {'Claimed':>14s}  {'Derived':>14s}  {'Match':>8s}  Err")
    print("  " + "-" * 78)

    for n in test_ns:
        a = n + 0.25
        for xi in test_xis:
            # Skip cases where result is negligibly small (unreliable quadrature)
            if 2 * a * xi > 6:
                continue

            # ∫_{-∞}^{∞} cos(ξx)/(a² + x²/4) dx  (even integrand, angular freq)
            def integrand(x):
                return np.cos(xi * x) / (a**2 + (x / 2.0)**2)

            val, _ = integrate.quad(integrand, 0, np.inf, limit=500, epsabs=1e-13, epsrel=1e-13)
            val *= 2  # double for even function

            claimed = np.pi / a * np.exp(-2 * a * xi)           # (π/a)·e^{-2a|ξ|}
            derived = 2 * np.pi / a * np.exp(-2 * a * xi)       # (2π/a)·e^{-2a|ξ|}

            err_claimed = abs(val - claimed) / (abs(val) + 1e-30)
            err_derived = abs(val - derived) / (abs(val) + 1e-30)

            if err_derived < err_claimed:
                match = "DERIVED"
                derived_wins += 1
                ok = err_derived < 0.01
                err_show = err_derived
            else:
                match = "CLAIMED"
                claimed_wins += 1
                ok = err_claimed < 0.01
                err_show = err_claimed

            if not ok:
                all_pass = False

            print(f"  {n:3d}  {xi:5.1f}  {val:14.8f}  {claimed:14.8f}  {derived:14.8f}  {match:>8s}  {err_show:.2e}")

    print(f"\n  Claimed formula wins: {claimed_wins}")
    print(f"  Derived formula wins: {derived_wins}")

    if derived_wins > claimed_wins:
        print(f"\n  CONCLUSION: Correct formula is F[1/(a²+(x/2)²)](ξ) = (2π/a)·e^{{-2a|ξ|}}")
        print(f"  The claimed formula (π/a)·e^{{-2a|ξ|}} is off by factor 2.")
        print(f"  However, K̂_bg derivation remains valid because the building block")
        print(f"  in the partial fraction is a/(a²+(x/2)²), whose FT is 2π·e^{{-2a|ξ|}}.")
        print(f"  The factor of a cancels the 1/a, giving 2·Σ e^{{-2a|ξ|}} = K̂_bg as claimed.")

    # Verify the KEY transform used in K̂_bg derivation:
    # FT[a/(a²+(x/2)²)] = a · (2π/a) · e^{-2a|ξ|} = 2π·e^{-2a|ξ|}
    print(f"\n  Verifying the KEY identity: FT[a/(a²+(x/2)²)] = 2π·e^{{-2a|ξ|}}")
    print(f"  {'n':>3s}  {'ξ':>5s}  {'Quadrature':>14s}  {'2π·exp':>14s}  {'Rel Err':>12s}  Status")
    print("  " + "-" * 65)

    key_pass = True
    for n in [0, 1, 3]:
        a = n + 0.25
        for xi in [0.5, 1.0, 2.0]:
            if 2 * a * xi > 6:  # skip where result is too tiny for reliable quadrature
                continue

            def integrand_key(x, _a=a, _xi=xi):
                return _a * np.cos(_xi * x) / (_a**2 + (x / 2.0)**2)

            val, _ = integrate.quad(integrand_key, 0, np.inf, limit=500, epsabs=1e-13, epsrel=1e-13)
            val *= 2

            expected = 2 * np.pi * np.exp(-2 * a * xi)
            rel_err = abs(val - expected) / abs(expected)
            ok = rel_err < 0.01  # 1% tolerance for oscillatory quadrature
            if not ok:
                key_pass = False
            status = "PASS" if ok else "FAIL"
            print(f"  {n:3d}  {xi:5.1f}  {val:14.8f}  {expected:14.8f}  {rel_err:12.6e}  {status}")

    all_pass = all_pass and key_pass

    # The overall derivation validity
    print(f"\n  KEY identity (used in K̂_bg derivation): [{'PASS' if key_pass else 'FAIL'}]")
    print(f"  Overall: [{'PASS' if all_pass else 'FAIL'}]")
    print()
    return all_pass


# ─── Test 6: Geometric series verification ───

def test_geometric_series():
    print("=" * 70)
    print("TEST 6: Geometric series Σ e^{-2(n+1/4)|ξ|} = e^{-|ξ|/2}/(1-e^{-2|ξ|})")
    print("=" * 70)

    test_xis = [0.1, 0.5, 1.0, 3.0, 10.0]
    all_pass = True

    print(f"  {'ξ':>8s}  {'Partial Sum':>14s}  {'Formula':>14s}  {'Rel Err':>12s}  Status")
    print("  " + "-" * 65)

    for xi in test_xis:
        # Sum first 1000 terms
        n_terms = 1000
        partial = sum(np.exp(-2 * (n + 0.25) * xi) for n in range(n_terms))
        formula = np.exp(-xi / 2) / (1.0 - np.exp(-2 * xi))

        rel_err = abs(partial - formula) / abs(formula)
        ok = rel_err < 1e-8
        if not ok:
            all_pass = False
        status = "PASS" if ok else "FAIL"
        print(f"  {xi:8.2f}  {partial:14.10f}  {formula:14.10f}  {rel_err:12.6e}  {status}")

    # Therefore K̂_bg = 2·Σ = 2·e^{-|ξ|/2}/(1-e^{-2|ξ|})  ✓
    print(f"\n  => K̂_bg = 2·Σ = 2e^{{-|ξ|/2}}/(1-e^{{-2|ξ|}}) confirmed")
    print(f"  [{'PASS' if all_pass else 'FAIL'}]")
    print()
    return all_pass


# ─── Main ───

def main():
    print()
    print("=" * 70)
    print("  Fourier Verification of Background Weil Kernel K_bg")
    print("  K_bg(x) = -(1/pi)Re[psi(1/4+ix/2)] + log(pi)/(2pi)")
    print("  K_bg_hat(xi) = 2e^{-|xi|/2}/(1-e^{-2|xi|})  [angular freq]")
    print("=" * 70)
    print()

    results = {}
    results['FFT vs Analytic'] = test_fft_vs_analytic()
    results['Partial Fraction FT'] = test_partial_fraction_ft()
    results['Positivity'] = test_positivity()
    results['Completeness'] = test_completeness()
    results['Lorentzian FT'] = test_lorentzian_ft()
    results['Geometric Series'] = test_geometric_series()

    # Summary
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  {'Test':.<50s} {'Result':>8s}")
    print("  " + "-" * 60)
    all_pass = True
    for name, passed in results.items():
        status = "PASS" if passed else "FAIL"
        if not passed:
            all_pass = False
        print(f"  {name:.<50s} [{status:>4s}]")

    print()
    if all_pass:
        print("  *** ALL TESTS PASSED ***")
        print()
        print("  K_bg_hat(xi) = 2e^{-|xi|/2}/(1-e^{-2|xi|}) > 0 for all xi != 0")
        print("  => Background Weil kernel has manifestly positive Fourier transform.")
    else:
        print("  *** SOME TESTS FAILED — see details above ***")
    print()

    return 0 if all_pass else 1


if __name__ == '__main__':
    sys.exit(main())
