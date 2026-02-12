#!/usr/bin/env python3
"""
CONSTRAINED OFF-LINE ZERO TEST

Previous tests moved ALL zeros off-line simultaneously → σ*(N) → 0.
But this is unphysical: if RH fails, only isolated zeros are off-line
while the vast majority remain on Re(s) = 1/2.

This script tests:
1. Keep all known zeros on-line EXCEPT one pair ρ_k, 1-ρ̄_k
2. Move that pair to 1/2 ± σ + iγ_k (off-line)
3. Check whether APT still holds

Key insight: the on-line zeros contribute K_on (PSD by Bochner) which
STRENGTHENS APT. Only the isolated off-line pair fights it.

If APT holds for all single-pair removals at any σ > 0, then
isolated off-line zeros are incompatible with the Weil positivity
framework — which would be strong evidence toward RH.

Additionally tests the EXPLICIT FORMULA CONSTRAINT:
Moving a zero off-line changes the explicit formula sum. The geometric
side (primes, poles, archimedean) is FIXED. So the remaining zeros
must compensate → the off-line configuration has a "cost."
"""

import mpmath
import numpy as np
import time
import sys

mpmath.mp.dps = 30


def sieve_primes(bound):
    sieve = [True] * (int(bound) + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(bound**0.5) + 1):
        if sieve[i]:
            for j in range(i*i, int(bound)+1, i):
                sieve[j] = False
    return [p for p in range(2, int(bound)+1) if sieve[p]]


def K_bg(x):
    """Archimedean background kernel."""
    if abs(x) < 1e-14:
        x = 1e-12
    arg = mpmath.mpc(0.25, x / 2)
    psi = mpmath.digamma(arg)
    return float(-mpmath.re(psi) / mpmath.pi + mpmath.log(mpmath.pi) / (2 * mpmath.pi))


def build_weil_matrix(primes, m_max, zeros_online, zeros_offline_pairs, sigma_off):
    """
    Build the Weil matrix with:
    - zeros_online: list of γ values on the critical line (Re = 1/2)
    - zeros_offline_pairs: list of γ values moved off-line to Re = 1/2 ± σ_off
    - sigma_off: how far off-line the pair is moved

    The kernel contribution from each zero type:
    - On-line (σ=0): cos(γx) / (1/4 + γ²)  (contributes to PSD part)
    - Off-line (σ>0): cosh(σx)·cos(γx) / |ρ(1-ρ)|  (can have either sign)
    """
    labels = [(p, a) for p in primes for a in range(1, m_max + 1)]
    N = len(labels)
    M = np.zeros((N, N))

    for i, (pi, ai) in enumerate(labels):
        lpi = np.log(pi)
        for j, (pj, aj) in enumerate(labels):
            lpj = np.log(pj)
            x = ai * lpi - aj * lpj
            is_diag = (i == j)

            # Background kernel (archimedean)
            kb = K_bg(x)

            # On-line zeros: standard contribution
            kz_on = 0
            for gamma in zeros_online:
                kz_on += np.cos(gamma * x) / (0.25 + gamma**2)
            kz_on /= np.pi

            # Off-line zeros: modified contribution with cosh factor
            kz_off = 0
            for gamma in zeros_offline_pairs:
                beta = 0.5 + sigma_off
                rho_prod = abs(complex(beta, gamma) * complex(1 - beta, -gamma))
                kz_off += np.cosh(sigma_off * x) * np.cos(gamma * x) / rho_prod
            kz_off /= np.pi

            K_val = (1.0 if is_diag else 0.0) + kb + kz_on + kz_off
            M[i, j] = -np.sqrt(lpi * lpj) / (pi**(ai/2) * pj**(aj/2)) * K_val

    M = (M + M.T) / 2
    return M, N


def max_primitive_eigenvalue(M, N):
    """Get max eigenvalue in the primitive subspace."""
    v = np.ones(N) / np.sqrt(N)
    P = np.eye(N) - np.outer(v, v)
    Mp = P @ M @ P
    Mp = (Mp + Mp.T) / 2
    eigs = np.sort(np.linalg.eigvalsh(Mp))
    nz = eigs[np.abs(eigs) > 1e-14]
    return nz[-1] if len(nz) > 0 else 0.0


def test_single_pair_offline():
    """
    Test 1: Move a SINGLE pair of zeros off-line.
    Keep all other zeros on the critical line.
    Check if APT holds.
    """
    print("=" * 70)
    print("  TEST 1: SINGLE ZERO PAIR OFF-LINE")
    print("=" * 70)
    print()
    print("  Move one pair (ρ_k, 1-ρ̄_k) off-line by σ,")
    print("  keep remaining zeros on Re(s) = 1/2.")
    print()

    n_zeros = 200
    print(f"  Computing {n_zeros} zeta zeros...")
    all_zeros = [float(mpmath.im(mpmath.zetazero(k))) for k in range(1, n_zeros + 1)]
    print("  Done.")
    print()

    pb = 97
    m_max = 3
    primes = sieve_primes(pb)

    # For each zero k, move it off-line and test APT
    sigma_values = [0.1, 0.2, 0.3, 0.4, 0.49]
    test_zeros = [0, 1, 2, 4, 9, 19, 49, 99]  # which zero to move (0-indexed)

    print(f"  Primes ≤ {pb}, m_max = {m_max}, N = {len(primes)*m_max}")
    print()
    print(f"  {'Zero#':>6s} {'γ':>10s}", end="")
    for sigma in sigma_values:
        print(f" | {'σ='+str(sigma):>10s}", end="")
    print()
    print("  " + "-" * (18 + 13 * len(sigma_values)))

    for k in test_zeros:
        gamma_k = all_zeros[k]
        print(f"  {k+1:6d} {gamma_k:10.3f}", end="")

        for sigma in sigma_values:
            # All zeros on-line except zero k
            online = [g for i, g in enumerate(all_zeros) if i != k]
            offline = [gamma_k]

            M, N = build_weil_matrix(primes, m_max, online, offline, sigma)
            max_eig = max_primitive_eigenvalue(M, N)
            apt = max_eig < 1e-10
            sign = '+' if max_eig > 0 else ''
            marker = ' ✓' if apt else ' ✗'
            print(f" | {sign}{max_eig:8.2e}{marker}", end="")

        print()

    print()
    print("  KEY: ✓ = APT holds (zero can't be off-line)")
    print("       ✗ = APT fails (off-line zero escapes detection)")
    print()


def test_multiple_pairs_offline():
    """
    Test 2: Move K pairs of zeros off-line simultaneously.
    How many simultaneous off-line zeros can APT detect?
    """
    print("=" * 70)
    print("  TEST 2: MULTIPLE ZERO PAIRS OFF-LINE")
    print("=" * 70)
    print()
    print("  Move the first K zeros off-line, keep rest on Re(s) = 1/2.")
    print()

    n_zeros = 200
    print(f"  Computing {n_zeros} zeta zeros...")
    all_zeros = [float(mpmath.im(mpmath.zetazero(k))) for k in range(1, n_zeros + 1)]
    print("  Done.")
    print()

    pb = 97
    m_max = 3
    primes = sieve_primes(pb)
    sigma = 0.3

    print(f"  σ = {sigma}, primes ≤ {pb}, N = {len(primes)*m_max}")
    print()
    print(f"  {'K pairs':>8s} {'max prim eig':>14s} {'APT?':>6s} {'ratio on/off':>14s}")
    print("  " + "-" * 50)

    for k_off in [1, 2, 3, 5, 10, 20, 50, 100, 150, 200]:
        online = all_zeros[k_off:]
        offline = all_zeros[:k_off]

        M, N = build_weil_matrix(primes, m_max, online, offline, sigma)
        max_eig = max_primitive_eigenvalue(M, N)
        apt = "YES" if max_eig < 1e-10 else "NO"
        ratio = len(online) / max(len(offline), 1)

        print(f"  {k_off:8d} {max_eig:+14.6e} {apt:>6s} {ratio:14.1f}")

    print()
    print("  This shows the TRANSITION: how many off-line zeros")
    print("  can accumulate before APT breaks down.")
    print()


def test_explicit_formula_cost():
    """
    Test 3: The explicit formula constraint.

    The explicit formula says:
      Σ_ρ h(γ_ρ) = geometric_side(h)

    for all test functions h. The geometric side is FIXED by primes.
    If we move a zero off-line, the sum changes, and other zeros must
    compensate. Quantify this "cost."

    For h(t) = e^{-αt²} (Gaussian):
      Σ_ρ e^{-α γ_ρ²} = geometric_side(α)

    Moving ρ_k from 1/2+iγ_k to 1/2+σ+iγ_k changes the LHS by:
      Δ = e^{-α(γ_k + iσ)²} - e^{-αγ_k²}  (approximately)

    This must be compensated by shifts in other zeros.
    """
    print("=" * 70)
    print("  TEST 3: EXPLICIT FORMULA COST OF OFF-LINE ZEROS")
    print("=" * 70)
    print()

    n_zeros = 200
    print(f"  Computing {n_zeros} zeta zeros...")
    all_zeros = [float(mpmath.im(mpmath.zetazero(k))) for k in range(1, n_zeros + 1)]
    print("  Done.")
    print()

    # For a Gaussian test function h(t) = e^{-αt²}
    # The spectral contribution from a zero at ρ = β + iγ is:
    #   ĥ(γ) when β = 1/2 (standard)
    # When β ≠ 1/2, the Weil functional contribution changes.
    #
    # The "cost" is the change in the Weil functional:
    #   ΔW = W(h; zeros with ρ_k off-line) - W(h; all zeros on-line)

    print("  For each zero ρ_k, compute the 'cost' of moving it off-line:")
    print("  ΔW = change in Weil positivity functional")
    print()

    pb = 47
    m_max = 3
    primes = sieve_primes(pb)
    labels = [(p, a) for p in primes for a in range(1, m_max + 1)]

    # Use specific test vectors that probe different scales
    N = len(labels)

    # Random primitive test vectors
    np.random.seed(42)
    n_tests = 50
    test_vecs = []
    for _ in range(n_tests):
        c = np.random.randn(N)
        c -= c.mean()  # make primitive
        c /= np.linalg.norm(c)
        test_vecs.append(c)

    sigmas = [0.1, 0.2, 0.3, 0.4, 0.49]

    print(f"  Primes ≤ {pb}, N = {N}, {n_tests} random primitive test vectors")
    print()
    print(f"  {'Zero#':>6s} {'γ':>10s}", end="")
    for sigma in sigmas:
        print(f" | {'σ='+str(sigma):>10s}", end="")
    print("    (avg ΔW / W_online)")
    print("  " + "-" * (18 + 13 * len(sigmas) + 20))

    for k in [0, 1, 2, 4, 9, 19, 49]:
        gamma_k = all_zeros[k]
        print(f"  {k+1:6d} {gamma_k:10.3f}", end="")

        # Baseline: all zeros on-line
        M_base, _ = build_weil_matrix(primes, m_max, all_zeros, [], 0.0)
        W_base_vals = [c @ M_base @ c for c in test_vecs]
        avg_W_base = np.mean(W_base_vals)

        for sigma in sigmas:
            online = [g for i, g in enumerate(all_zeros) if i != k]
            offline = [gamma_k]
            M_off, _ = build_weil_matrix(primes, m_max, online, offline, sigma)
            W_off_vals = [c @ M_off @ c for c in test_vecs]
            avg_W_off = np.mean(W_off_vals)

            delta = avg_W_off - avg_W_base
            ratio = delta / abs(avg_W_base) if abs(avg_W_base) > 1e-15 else 0

            print(f" | {ratio:+10.4f}", end="")

        print()

    print()
    print("  Positive ΔW/W means moving the zero off-line INCREASES")
    print("  the Weil functional (good for APT — makes it harder to violate).")
    print("  Negative means it DECREASES (bad — weakens APT).")
    print()


def test_worst_case_direction():
    """
    Test 4: Find the test vector that is MOST affected by a single
    off-line zero. This is the "worst case" direction.

    If even the worst-case direction maintains APT, then a single
    off-line zero is provably incompatible with the framework.
    """
    print("=" * 70)
    print("  TEST 4: WORST-CASE DIRECTION FOR SINGLE OFF-LINE ZERO")
    print("=" * 70)
    print()

    n_zeros = 200
    print(f"  Computing {n_zeros} zeta zeros...")
    all_zeros = [float(mpmath.im(mpmath.zetazero(k))) for k in range(1, n_zeros + 1)]
    print("  Done.")
    print()

    pb = 97
    m_max = 3
    primes = sieve_primes(pb)
    N_mat = len(primes) * m_max

    print(f"  Primes ≤ {pb}, N = {N_mat}")
    print()

    sigma = 0.3
    print(f"  σ = {sigma}")
    print()

    # For each zero, compute the PERTURBATION matrix ΔM = M_off - M_on
    # The worst-case direction is the eigenvector of ΔM with largest eigenvalue
    # (restricted to primitive subspace)

    print(f"  {'Zero#':>6s} {'γ':>10s} {'max Δλ':>12s} {'max λ(M_on)':>14s} "
          f"{'max λ(M_off)':>14s} {'APT?':>6s}")
    print("  " + "-" * 70)

    v = np.ones(N_mat) / np.sqrt(N_mat)
    P = np.eye(N_mat) - np.outer(v, v)

    # Baseline: all on-line
    M_base, _ = build_weil_matrix(primes, m_max, all_zeros, [], 0.0)
    Mp_base = P @ M_base @ P
    Mp_base = (Mp_base + Mp_base.T) / 2
    eigs_base = np.linalg.eigvalsh(Mp_base)
    nz_base = eigs_base[np.abs(eigs_base) > 1e-14]
    max_base = nz_base[-1] if len(nz_base) > 0 else 0

    for k in [0, 1, 2, 4, 9, 19, 49, 99]:
        gamma_k = all_zeros[k]

        online = [g for i, g in enumerate(all_zeros) if i != k]
        offline = [gamma_k]

        M_off, _ = build_weil_matrix(primes, m_max, online, offline, sigma)
        Mp_off = P @ M_off @ P
        Mp_off = (Mp_off + Mp_off.T) / 2

        # Perturbation
        delta_M = Mp_off - Mp_base
        delta_M = (delta_M + delta_M.T) / 2
        delta_eigs = np.linalg.eigvalsh(delta_M)
        max_delta = delta_eigs[-1]

        eigs_off = np.linalg.eigvalsh(Mp_off)
        nz_off = eigs_off[np.abs(eigs_off) > 1e-14]
        max_off = nz_off[-1] if len(nz_off) > 0 else 0

        apt = "YES" if max_off < 1e-10 else "NO"

        print(f"  {k+1:6d} {gamma_k:10.3f} {max_delta:+12.4e} {max_base:+14.6e} "
              f"{max_off:+14.6e} {apt:>6s}")

    print()

    # Also test: what if we move the zero that causes MAXIMUM perturbation?
    print("  SCALING TEST: same zero (#1), increasing matrix size:")
    print()
    print(f"  {'P':>5s} {'N':>5s} {'max λ(on)':>14s} {'max λ(off)':>14s} "
          f"{'Δλ':>14s} {'Δλ/|λ_on|':>12s} {'APT?':>6s}")
    print("  " + "-" * 75)

    for pb_test in [23, 47, 97, 197, 307]:
        primes_t = sieve_primes(pb_test)
        N_t = len(primes_t) * m_max

        v_t = np.ones(N_t) / np.sqrt(N_t)
        P_t = np.eye(N_t) - np.outer(v_t, v_t)

        # All on-line
        M_on, _ = build_weil_matrix(primes_t, m_max, all_zeros, [], 0.0)
        Mp_on = P_t @ M_on @ P_t
        Mp_on = (Mp_on + Mp_on.T) / 2
        eigs_on = np.linalg.eigvalsh(Mp_on)
        nz_on = eigs_on[np.abs(eigs_on) > 1e-14]
        max_on = nz_on[-1] if len(nz_on) > 0 else 0

        # Zero #1 off-line at σ = 0.3
        online = all_zeros[1:]
        offline = [all_zeros[0]]
        M_of, _ = build_weil_matrix(primes_t, m_max, online, offline, sigma)
        Mp_of = P_t @ M_of @ P_t
        Mp_of = (Mp_of + Mp_of.T) / 2
        eigs_of = np.linalg.eigvalsh(Mp_of)
        nz_of = eigs_of[np.abs(eigs_of) > 1e-14]
        max_of = nz_of[-1] if len(nz_of) > 0 else 0

        delta = max_of - max_on
        ratio = delta / abs(max_on) if abs(max_on) > 1e-15 else float('inf')
        apt = "YES" if max_of < 1e-10 else "NO"

        print(f"  {pb_test:5d} {N_t:5d} {max_on:+14.6e} {max_of:+14.6e} "
              f"{delta:+14.6e} {ratio:+12.4f} {apt:>6s}")

    print()


def test_density_constraint():
    """
    Test 5: Zero density estimates constrain off-line zeros.

    Ingham/Huxley: N(σ, T) ≤ C · T^{A(1-σ)} (log T)^B
    where N(σ,T) = number of zeros with Re(ρ) > σ, |Im(ρ)| < T.

    For σ close to 1/2: A ≈ 3 (Ingham), so N(σ,T) ∝ T^{3(1-2σ)·½} for σ near 1/2.
    The DENSITY of off-line zeros is bounded.

    Test: if off-line zeros are sparse (as density estimates demand),
    does APT hold?
    """
    print("=" * 70)
    print("  TEST 5: DENSITY-CONSTRAINED OFF-LINE ZEROS")
    print("=" * 70)
    print()

    n_zeros = 200
    print(f"  Computing {n_zeros} zeta zeros...")
    all_zeros = [float(mpmath.im(mpmath.zetazero(k))) for k in range(1, n_zeros + 1)]
    print("  Done.")
    print()

    pb = 97
    m_max = 3
    primes = sieve_primes(pb)

    # Ingham density estimate: at most C·T^{3(1-σ)} zeros with Re > σ up to height T
    # For our 200 zeros, T ≈ γ_200 ≈ 248.
    T_max = all_zeros[-1]

    print(f"  T_max ≈ {T_max:.1f} (height of 200th zero)")
    print()

    sigma = 0.3

    # Density estimate: N(σ, T) ≤ C · T^{3(1-σ)} log T
    # For σ = 0.3: N ≤ C · T^{2.1} log T ≈ huge → density doesn't help much
    # For σ = 0.4: N ≤ C · T^{1.8} log T → still large
    # For σ = 0.49: N ≤ C · T^{1.53} → still many possible

    # Better approach: use the PROPORTION. At σ > 0, at most a fraction
    # f(σ) of zeros can be off-line, where f depends on σ and the
    # density estimate. For σ close to 1/2, f → 0.

    # Selberg: at least 1/3 of zeros are on the line (Levinson: > 1/3, Conrey: > 2/5)
    # So at most 3/5 can be off-line

    proportions = [0.01, 0.05, 0.1, 0.2, 0.4, 0.6]

    print(f"  σ = {sigma}, primes ≤ {pb}")
    print()
    print(f"  {'Fraction off':>14s} {'K off':>6s} {'K on':>6s} "
          f"{'max prim eig':>14s} {'APT?':>6s}")
    print("  " + "-" * 55)

    for frac in proportions:
        k_off = max(1, int(frac * n_zeros))
        k_on = n_zeros - k_off

        # Move the FIRST k_off zeros off-line (worst case: low-lying zeros)
        online = all_zeros[k_off:]
        offline = all_zeros[:k_off]

        M, N = build_weil_matrix(primes, m_max, online, offline, sigma)
        max_eig = max_primitive_eigenvalue(M, N)
        apt = "YES" if max_eig < 1e-10 else "NO"

        print(f"  {frac:14.2f} {k_off:6d} {k_on:6d} "
              f"{max_eig:+14.6e} {apt:>6s}")

    print()

    # Now test with scattered off-line zeros (not just the first k)
    print("  SCATTERED OFF-LINE ZEROS (every Mth zero, not just the first K):")
    print()
    print(f"  {'Pattern':>20s} {'K off':>6s} {'max prim eig':>14s} {'APT?':>6s}")
    print("  " + "-" * 55)

    for spacing in [2, 5, 10, 20, 50]:
        offline_indices = list(range(0, n_zeros, spacing))
        k_off = len(offline_indices)
        online = [g for i, g in enumerate(all_zeros) if i not in offline_indices]
        offline = [all_zeros[i] for i in offline_indices]

        M, N = build_weil_matrix(primes, m_max, online, offline, sigma)
        max_eig = max_primitive_eigenvalue(M, N)
        apt = "YES" if max_eig < 1e-10 else "NO"

        print(f"  {'every '+str(spacing)+'th':>20s} {k_off:6d} "
              f"{max_eig:+14.6e} {apt:>6s}")

    print()


def main():
    t0 = time.time()

    print()
    print("╔" + "═" * 68 + "╗")
    print("║  CONSTRAINED OFF-LINE ZERO TEST                                    ║")
    print("║  Testing APT with realistic (sparse) off-line zero configurations  ║")
    print("╚" + "═" * 68 + "╝")
    print()

    test_single_pair_offline()
    print()
    test_worst_case_direction()
    print()
    test_multiple_pairs_offline()
    print()
    test_density_constraint()
    print()
    test_explicit_formula_cost()

    elapsed = time.time() - t0
    print()
    print(f"  Total time: {elapsed:.1f}s")
    print()


if __name__ == '__main__':
    main()
