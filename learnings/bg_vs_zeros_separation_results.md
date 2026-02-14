# BG vs Zeros Separation — Scaling Analysis

**Question**: Does the M_bg spectral gap grow with P_0, swamping the zeros contribution?

## Setup

Three matrices built for each P_0:
- **M_bg**: delta(x) + K_bg(x) only (no zeros)
- **M_zeros**: K_zeros(x) only (just the zero-sum kernel, 200 zeros)
- **M_full**: delta(x) + K_bg(x) + K_zeros(x) (everything)

All eigenvalues computed on the **primitive subspace** using the correct projection.

## Results

| P_0 | N | bg gap | zeros norm | full gap | ratio | bg dominates? | bg max prim | full max prim |
|-----|---|--------|------------|----------|-------|---------------|-------------|---------------|
| 200 | 138 | 1.723670e+00 | 4.098190e-03 | 1.724599e+00 | 420.59 | YES | -6.835713e-07 | -6.866380e-07 |
| 300 | 186 | 1.817688e+00 | 4.306888e-03 | 1.818599e+00 | 422.04 | YES | -2.451340e-07 | -2.471602e-07 |
| 500 | 190 | 1.912594e+00 | 4.678131e-03 | 1.913649e+00 | 408.84 | YES | -2.536921e-05 | -2.550833e-05 |
| 750 | 264 | 2.007871e+00 | 4.929044e-03 | 2.008960e+00 | 407.35 | YES | -1.203407e-05 | -1.206270e-05 |
| 1000 | 336 | 2.076920e+00 | 5.226856e-03 | 2.078072e+00 | 397.36 | YES | -6.984299e-06 | -7.003253e-06 |

## Growth Analysis

| Transition | bg gap growth | zeros norm growth |
|------------|---------------|-------------------|
| P_0 200 -> 300 | x1.055 | x1.051 |
| P_0 300 -> 500 | x1.052 | x1.086 |
| P_0 500 -> 750 | x1.050 | x1.054 |
| P_0 750 -> 1000 | x1.034 | x1.060 |

## Conclusion

**M_bg dominates at all tested scales, but the dominance ratio does not grow monotonically.**

Further investigation at larger scales may be needed to determine the asymptotic trend.
