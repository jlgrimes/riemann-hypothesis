# Random vs Real Primes: Is APT a Prime-Specific Property?

**The key question**: Does the Arithmetic Positivity Theorem (all primitive
eigenvalues of the Weil matrix <= 0) hold because of something special about
primes, or is it a property of the kernel K(x) = K_bg(x) + K_zeros(x)?

## Configuration

- P_0 = 200, m_max = 3
- 46 real primes, matrix size 138x138
- 50 random pseudo-prime trials (density ~ 1/log(n))
- 200 zeta zeros, mpmath dps = 30
- Total computation time: 35.2s

---

## Results Summary

| Test | max prim eig | min prim eig | APT? | Z-score |
|------|-------------|-------------|------|---------|
| **Real primes** | -6.866380e-07 | -1.724599e+00 | **YES** | — |
| Random (mean +/- std) | -7.540562e-07 +/- 7.94e-08 | — | 50/50 (100%) | +0.85 |
| Shifted primes (p+1) | -6.769503e-07 | -1.645057e+00 | YES | — |
| Composites only | -7.570363e-07 | -1.055619e+00 | YES | — |

## Random Trial Distribution

| Statistic | max prim eigenvalue |
|-----------|-------------------|
| Mean | -7.54056216e-07 |
| Std | 7.93539705e-08 |
| Min | -1.04090116e-06 |
| Max | -6.69181660e-07 |
| Median | -7.28233441e-07 |
| APT pass rate | 50/50 (100%) |

**Real primes max eigenvalue**: -6.86637954e-07
**Z-score**: +0.8496

### Individual Random Trial Results

| Trial | max prim eig | APT? |
|-------|-------------|------|
| 1 | -7.292799e-07 | YES |
| 2 | -8.752855e-07 | YES |
| 3 | -7.599101e-07 | YES |
| 4 | -7.147421e-07 | YES |
| 5 | -7.699325e-07 | YES |
| 6 | -8.841753e-07 | YES |
| 7 | -7.927716e-07 | YES |
| 8 | -6.692877e-07 | YES |
| 9 | -7.703326e-07 | YES |
| 10 | -1.040901e-06 | YES |
| 11 | -6.789136e-07 | YES |
| 12 | -7.256474e-07 | YES |
| 13 | -6.856840e-07 | YES |
| 14 | -7.664390e-07 | YES |
| 15 | -8.620551e-07 | YES |
| 16 | -6.789192e-07 | YES |
| 17 | -7.086015e-07 | YES |
| 18 | -7.599317e-07 | YES |
| 19 | -7.874471e-07 | YES |
| 20 | -6.788807e-07 | YES |
| 21 | -7.188781e-07 | YES |
| 22 | -6.770496e-07 | YES |
| 23 | -7.058832e-07 | YES |
| 24 | -8.297194e-07 | YES |
| 25 | -7.043616e-07 | YES |
| 26 | -7.472208e-07 | YES |
| 27 | -7.145882e-07 | YES |
| 28 | -7.293769e-07 | YES |
| 29 | -6.949944e-07 | YES |
| 30 | -7.831444e-07 | YES |
| 31 | -7.509491e-07 | YES |
| 32 | -6.951227e-07 | YES |
| 33 | -6.984932e-07 | YES |
| 34 | -7.271870e-07 | YES |
| 35 | -9.931300e-07 | YES |
| 36 | -6.691817e-07 | YES |
| 37 | -8.645536e-07 | YES |
| 38 | -7.293869e-07 | YES |
| 39 | -7.889012e-07 | YES |
| 40 | -7.059180e-07 | YES |
| 41 | -8.530955e-07 | YES |
| 42 | -6.789046e-07 | YES |
| 43 | -8.080422e-07 | YES |
| 44 | -7.159864e-07 | YES |
| 45 | -6.866222e-07 | YES |
| 46 | -6.887266e-07 | YES |
| 47 | -6.768628e-07 | YES |
| 48 | -8.073545e-07 | YES |
| 49 | -7.943223e-07 | YES |
| 50 | -7.257157e-07 | YES |

---

## Analysis

### Finding: APT is a KERNEL property, not prime-specific

Random numbers satisfy APT at a 100% rate, comparable to real primes.
This strongly suggests that the negative-definiteness of the Weil matrix on the
primitive subspace arises from the structure of the kernel K(x) = K_bg(x) + K_zeros(x),
not from any special property of the prime numbers themselves.

**Implication**: The Weil matrix encodes the Riemann Hypothesis through its kernel,
which is built from the zeta zeros. The primes serve as a 'basis' for sampling
the kernel, but other bases work too.

**Shifted primes {p+1}**: Shifting primes by 1 still satisfies APT, suggesting
the property is robust to small perturbations of the input set.

**Composites**: Even composites satisfy APT, further suggesting the kernel dominates.

---

## Conclusion

The Weil matrix's spectral negativity on the primitive subspace is primarily a
consequence of the kernel structure (K_bg + K_zeros) rather than a unique property
of prime numbers. The kernel, which encodes the zeta zeros via the explicit formula,
produces negative primitive eigenvalues regardless of which numbers are used as inputs.

*Generated: 2026-02-14 10:37:11*
