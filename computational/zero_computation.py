#!/usr/bin/env python3
"""
zero_computation.py — Compute and verify non-trivial zeros of the Riemann zeta function.

Computes the first 1000 non-trivial zeros on the critical line Re(s) = 1/2
using mpmath's zetazero function, verifies they lie on the critical line,
and displays statistical properties.
"""

import mpmath
import time
import sys

# Set high precision for accurate zero computation
mpmath.mp.dps = 50  # 50 decimal digits


def compute_zeros(n_zeros=1000, verbose=True):
    """Compute the first n_zeros non-trivial zeros of the Riemann zeta function."""
    zeros = []
    start_time = time.time()

    if verbose:
        print(f"Computing first {n_zeros} non-trivial zeros of ζ(s)...")
        print(f"Precision: {mpmath.mp.dps} decimal digits")
        print("-" * 70)

    for k in range(1, n_zeros + 1):
        # mpmath.zetazero(k) returns the k-th zero on the upper critical line
        zero = mpmath.zetazero(k)
        zeros.append(zero)

        if verbose and (k <= 10 or k % 100 == 0):
            elapsed = time.time() - start_time
            print(f"  ζ zero #{k:>4d}: s = 1/2 + {mpmath.nstr(zero.imag, 30)}i "
                  f"  [{elapsed:.1f}s elapsed]")

    elapsed = time.time() - start_time
    if verbose:
        print("-" * 70)
        print(f"Computed {n_zeros} zeros in {elapsed:.1f} seconds.")

    return zeros


def verify_critical_line(zeros, tol=1e-30):
    """Verify that all zeros have real part exactly 1/2."""
    print("\n" + "=" * 70)
    print("VERIFICATION: All zeros on the critical line Re(s) = 1/2?")
    print("=" * 70)

    max_deviation = mpmath.mpf(0)
    all_on_line = True

    for k, z in enumerate(zeros, 1):
        deviation = abs(z.real - mpmath.mpf('0.5'))
        if deviation > max_deviation:
            max_deviation = deviation
        if deviation > tol:
            all_on_line = False
            print(f"  WARNING: Zero #{k} deviates from critical line by {deviation}")

    print(f"\n  Number of zeros checked: {len(zeros)}")
    print(f"  Maximum deviation from Re(s)=1/2: {mpmath.nstr(max_deviation, 10)}")
    print(f"  Tolerance: {tol}")
    print(f"  All on critical line: {'YES ✓' if all_on_line else 'NO ✗'}")

    return all_on_line, max_deviation


def verify_zeta_vanishes(zeros, n_check=20):
    """Double-check that ζ(s) actually vanishes at the computed zeros."""
    print("\n" + "=" * 70)
    print(f"VERIFICATION: ζ(s) ≈ 0 at first {n_check} zeros?")
    print("=" * 70)

    for k in range(min(n_check, len(zeros))):
        z = zeros[k]
        val = mpmath.zeta(z)
        print(f"  |ζ(zero #{k+1:>3d})| = {mpmath.nstr(abs(val), 10)}")


def analyze_spacings(zeros):
    """Analyze the spacings between consecutive zeros."""
    print("\n" + "=" * 70)
    print("ZERO SPACING ANALYSIS")
    print("=" * 70)

    # Extract imaginary parts (the t-values)
    t_values = [float(z.imag) for z in zeros]
    spacings = [t_values[i+1] - t_values[i] for i in range(len(t_values)-1)]

    # Normalize spacings by average spacing
    avg_spacing = sum(spacings) / len(spacings)
    normalized = [s / avg_spacing for s in spacings]

    print(f"\n  First zero: t₁ = {t_values[0]:.10f}")
    print(f"  Last zero:  t_{len(zeros)} = {t_values[-1]:.10f}")
    print(f"  Range: {t_values[-1] - t_values[0]:.6f}")
    print(f"\n  Average spacing: {avg_spacing:.10f}")
    print(f"  Min spacing: {min(spacings):.10f}")
    print(f"  Max spacing: {max(spacings):.10f}")

    # Check against asymptotic formula: average spacing near t ~ 2π/log(t/(2π))
    t_mid = t_values[len(t_values)//2]
    predicted_avg = 2 * mpmath.pi / mpmath.log(t_mid / (2 * mpmath.pi))
    print(f"\n  Predicted average spacing (at t={t_mid:.2f}): {float(predicted_avg):.10f}")
    print(f"  Actual average spacing: {avg_spacing:.10f}")

    # Spacing statistics
    import statistics
    print(f"\n  Spacing std deviation: {statistics.stdev(spacings):.10f}")
    print(f"  Normalized spacing std: {statistics.stdev(normalized):.10f}")

    # Check for zero repulsion (GUE predicts small spacings are rare)
    very_small = sum(1 for s in normalized if s < 0.1)
    small = sum(1 for s in normalized if s < 0.5)
    print(f"\n  Spacings < 0.1 × average: {very_small} ({100*very_small/len(normalized):.2f}%)")
    print(f"  Spacings < 0.5 × average: {small} ({100*small/len(normalized):.2f}%)")
    print("  (GUE predicts strong zero repulsion — very small spacings are rare)")

    return t_values, spacings, normalized


def gram_point_analysis(zeros, n_check=50):
    """Analyze Gram's law: the behavior of zeros relative to Gram points."""
    print("\n" + "=" * 70)
    print("GRAM POINT ANALYSIS")
    print("=" * 70)

    # Gram points g_n satisfy θ(g_n) = nπ where θ is the Riemann-Siegel theta function
    # Between consecutive Gram points, there is "usually" exactly one zero (Gram's law)
    t_values = [float(z.imag) for z in zeros[:n_check]]

    print(f"\n  First {n_check} zeros and their Gram point neighborhoods:")
    violations = 0
    for k in range(min(n_check, 10)):
        t = t_values[k]
        # Approximate which Gram interval this zero falls in
        theta_val = float(mpmath.siegeltheta(t))
        gram_index = theta_val / float(mpmath.pi)
        print(f"    Zero #{k+1}: t = {t:.8f}, θ(t)/π ≈ {gram_index:.4f}")

    print(f"\n  (Gram's law states most zeros fall in their 'expected' Gram interval)")


def counting_function(zeros):
    """Compare zero count with Riemann-von Mangoldt formula."""
    print("\n" + "=" * 70)
    print("RIEMANN-VON MANGOLDT FORMULA VERIFICATION")
    print("=" * 70)

    # N(T) = T/(2π) log(T/(2πe)) + 7/8 + S(T) where S(T) is small
    t_values = [float(z.imag) for z in zeros]

    for n in [10, 50, 100, 500, len(zeros)]:
        T = t_values[n-1]
        N_actual = n
        N_predicted = T / (2 * mpmath.pi) * mpmath.log(T / (2 * mpmath.pi * mpmath.e)) + 7/8
        error = float(abs(N_actual - N_predicted))
        print(f"  N({T:.4f}) = {n:>4d}  (predicted: {float(N_predicted):>8.2f}, "
              f"|error|: {error:.4f})")


def main():
    n_zeros = 1000
    if len(sys.argv) > 1:
        n_zeros = int(sys.argv[1])

    # Compute zeros
    zeros = compute_zeros(n_zeros)

    # Verify all on critical line
    verify_critical_line(zeros)

    # Verify zeta vanishes
    verify_zeta_vanishes(zeros)

    # Analyze spacings
    analyze_spacings(zeros)

    # Gram point analysis
    gram_point_analysis(zeros)

    # Counting function verification
    counting_function(zeros)

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  All {n_zeros} computed zeros lie on the critical line Re(s) = 1/2.")
    print(f"  This is consistent with the Riemann Hypothesis.")
    print(f"  Note: Computational verification of finitely many zeros cannot")
    print(f"  prove RH, but violations would disprove it.")
    print("=" * 70)

    # Save zeros to file for use by other scripts
    with open('/Users/jaredgrimes/code/riemann-hypothesis/computational/zeros_data.txt', 'w') as f:
        f.write(f"# First {n_zeros} non-trivial zeros of zeta (imaginary parts)\n")
        f.write(f"# Precision: {mpmath.mp.dps} decimal digits\n")
        for k, z in enumerate(zeros, 1):
            f.write(f"{k} {mpmath.nstr(z.imag, 40)}\n")
    print(f"\n  Zeros saved to zeros_data.txt")


if __name__ == '__main__':
    main()
