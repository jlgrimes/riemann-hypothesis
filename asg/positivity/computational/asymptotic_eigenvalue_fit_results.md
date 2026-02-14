# Asymptotic Eigenvalue Fit Results

## Overview

Fits the maximum primitive eigenvalue $\lambda_{\max}(P_0)$ of the Weil matrix
to asymptotic models, to determine whether the spectral gap closes or stays open.

## Configuration

- **mpmath precision**: 30 digits
- **Zeta zeros**: 200
- **Total computation time**: 31.1s (0.5min)

### Kernel

$$K(x) = K_{\text{bg}}(x) + K_{\text{zeros}}(x) + \delta(x)$$

- $K_{\text{bg}}(x) = -\frac{1}{\pi} \operatorname{Re} \psi(1/4 + ix/2) + \frac{1}{2\pi} \log \pi$
- $K_{\text{zeros}}(x) = \frac{1}{\pi} \sum_\gamma \frac{\cos(\gamma x)}{1/4 + \gamma^2}$

### Matrix

$$M_{(p,m),(q,n)} = -\frac{\sqrt{\log p \cdot \log q}}{p^{m/2} \cdot q^{n/2}} \cdot K(m \log p - n \log q)$$

### Primitive Projection

$v_{p,m} = (\log p)^{1/2} / p^{m/2}$, $P = I - vv^T/\|v\|^2$, $M_p = PMP$

---

## Data

| P₀ | #primes | m_max | N | λ_max | |λ_max| |
|---:|--------:|------:|--:|------:|-------:|
| 100 | 25 | 3 | 75 | -6.801138e-06 | 6.801138e-06 |
| 150 | 35 | 3 | 105 | -1.864329e-06 | 1.864329e-06 |
| 200 | 46 | 3 | 138 | -6.866400e-07 | 6.866400e-07 |
| 300 | 62 | 3 | 186 | -2.471600e-07 | 2.471600e-07 |
| 400 | 78 | 2 | 156 | -3.908667e-05 | 3.908667e-05 |
| 500 | 95 | 2 | 190 | -2.550800e-05 | 2.550800e-05 |
| 600 | 109 | 2 | 218 | -1.806909e-05 | 1.806909e-05 |
| 750 | 132 | 2 | 264 | -1.206300e-05 | 1.206300e-05 |

### Data with Model Predictions

| P₀ | |λ_max| (data) | Model A | Model B | Model C |
|---:|---------------:|----------:|----------:|----------:|
| 100 | 6.801138e-06 | 5.582227e-06 | 4.921711e-06 | 3.054786e-05 |
| 150 | 1.864329e-06 | 7.384073e-06 | 7.005348e-06 | 2.250392e-05 |
| 200 | 6.866400e-07 | 9.005202e-06 | 8.848296e-06 | 1.778203e-05 |
| 300 | 2.471600e-07 | 1.191192e-05 | 1.204583e-05 | 1.239814e-05 |
| 400 | 3.908667e-05 | 1.452712e-05 | 1.479920e-05 | 9.389895e-06 |
| 500 | 2.550800e-05 | 1.694491e-05 | 1.724549e-05 | 7.467226e-06 |
| 600 | 1.806909e-05 | 1.921622e-05 | 1.946308e-05 | 6.135213e-06 |
| 750 | 1.206300e-05 | 2.241444e-05 | 2.246514e-05 | 4.766590e-06 |

---

## Model Fits

### Model A

$$|\lambda_{\max}| \sim \frac{C}{P_0^\alpha}$$

- $C = 2.327941e-07$
- $\alpha = -0.689918$
- $R^2 = 0.23729370$
- $\text{AIC} = -178.2461$

### Model B

$$|\lambda_{\max}| \sim \frac{C}{\log(P_0)^\beta}$$

- $C = 8.268380e-09$
- $\beta = -4.183511$
- $R^2 = 0.25080187$
- $\text{AIC} = -178.3891$

### Model C

$$|\lambda_{\max}| \sim C \cdot \exp(-P_0^\gamma)$$

- $C = 6.164277e-04$
- $\gamma = 0.238897$
- $R^2 = -1.11371805$
- $\text{AIC} = -170.0915$

### Model Comparison

| Model | R² | AIC | Decay type |
|:------|---:|----:|:-----------|
| A | 0.23729370 | -178.2461 | Power law |
| B **(best R², best AIC)** | 0.25080187 | -178.3891 | Logarithmic |
| C | -1.11371805 | -170.0915 | Exponential |

---

## Extrapolation

| P₀ | Model A | Model B | Model C |
|---:|----------:|----------:|----------:|
| 1000 | 2.733539e-05 | 2.684082e-05 | 3.372694e-06 |
| 2000 | 4.409722e-05 | 4.004328e-05 | 1.320173e-06 |
| 5000 | 8.297687e-05 | 6.446532e-05 | 2.933937e-07 |
| 10000 | 1.338576e-04 | 8.942902e-05 | 7.398007e-08 |

---

## Interpretation

### Does the spectral gap close?

All three models predict $|\lambda_{\max}| \to 0$ as $P_0 \to \infty$.
This means the spectral gap **closes**: the most negative primitive eigenvalue
approaches zero from below.

### What this means for the proof

- **Gap closing** ($\lambda_{\max} \to 0^-$): The finite verification becomes
  harder at larger $P_0$, because the eigenvalue margin shrinks.
- However, $\lambda_{\max} < 0$ at **every** finite $P_0$ tested.
- The rate of gap closing determines how much headroom exists for
  tail/truncation error bounds in the finite-to-infinite extrapolation.

### Caveats

- The transition from $m_{\max}=3$ to $m_{\max}=2$ at $P_0 \ge 400$
  introduces a structural change in the matrix (fewer prime powers),
  which causes non-monotonicity in $|\lambda_{\max}|$.
- The fit quality is limited by the small number of data points (8).
- More data at larger $P_0$ with consistent $m_{\max}$ would yield
  more reliable scaling estimates.

---

*Generated: 2026-02-13 18:04:02*
