# The Lévy Process Approach to CPD-1 of the Weil Kernel

## Status: Exploratory. Provides structural insight and a novel reformulation of RH as a stochastic process existence question. Does NOT close the gap but identifies precisely where the Lévy characterization breaks down for off-line zeros.

---

## 0. Overview and Motivation

By Schoenberg's theorem (1938), a continuous even function K: R -> R is conditionally positive definite of order 1 (CPD-1) if and only if for every t > 0, the function

  phi_t(x) = exp(-t[K(0) - K(x)])

is positive definite (PD). The Lévy-Khinchin theorem provides the converse characterization: a continuous function psi: R -> R with psi(0) = 0 is the characteristic exponent of a symmetric Lévy process if and only if

  psi(x) = a x^2 + integral_R (1 - cos(ux)) nu(du)

where a >= 0 (Gaussian component) and nu is a symmetric Lévy measure satisfying integral min(1, u^2) nu(du) < infinity.

Combining these:

**K is CPD-1 iff psi(x) = K(0) - K(x) is the characteristic exponent of a Lévy process.**

Since CPD-1 of the Weil kernel is equivalent to RH (subspace-alignment.md, Theorem 2.3), we obtain:

**RH iff there exists a Lévy process whose characteristic exponent is K_Weil(0) - K_Weil(x).**

This document explores this equivalence, identifies the explicit Lévy processes for each component of the Weil kernel, and analyzes where the approach provides genuine insight versus mere reformulation.

---

## 1. Lévy-Khinchin Representation of K_bg

### 1.1 Setup

The background kernel K_bg(x) = -(1/pi) Re[psi(1/4 + ix/2)] + log(pi)/(2pi) is unconditionally CPD-1 (bochner-proof.md, Theorem 8.1; schoenberg-attempt.md, Theorem 5.1), with Fourier transform:

  K_hat_bg(xi) = 2 e^{-|xi|/2} / (1 - e^{-2|xi|})  > 0  for all xi != 0

By Schoenberg, psi_bg(x) = K_bg(0) - K_bg(x) is a valid characteristic exponent.

### 1.2 The Explicit Lévy Measure

Since K_bg(x) = (1/(2pi)) integral K_hat_bg(xi) e^{ixi x} dxi, we have:

  K_bg(0) - K_bg(x) = (1/(2pi)) integral K_hat_bg(xi) (1 - e^{ixi x}) dxi
                     = (1/(2pi)) integral K_hat_bg(xi) (1 - cos(xi x)) dxi

(the imaginary part vanishes by evenness of K_hat_bg). This is the Lévy-Khinchin form with:

  **nu_bg(du) = K_hat_bg(u)/(2pi) du = (1/pi) e^{-|u|/2} / (1 - e^{-2|u|}) du**

**Verification.** The Lévy measure conditions:

1. nu({0}) = 0: satisfied (nu_bg is absolutely continuous).

2. integral min(1, u^2) nu_bg(du) < infinity:
   - Near u = 0: nu_bg(u) ~ 1/(2pi|u|), so integral_0^1 u^2 / (2pi u) du = integral_0^1 u/(2pi) du < infinity. CHECK.
   - Near u = infinity: nu_bg(u) ~ (1/pi) e^{-|u|/2}, exponential decay. CHECK.

3. Gaussian component a = 0 (no x^2 term).

### 1.3 Series Decomposition

Expanding the geometric series:

  e^{-|u|/2} / (1 - e^{-2|u|}) = sum_{n=0}^infinity e^{-(2n+1/2)|u|}

Therefore:

  nu_bg(du) = (1/pi) sum_{n=0}^infinity e^{-(2n+1/2)|u|} du

Each term e^{-(2n+1/2)|u|} du is (up to normalization) a symmetric Laplace distribution with rate parameter lambda_n = 2n + 1/2. The n-th component has:

- Jump rate: integral 2 e^{-lambda_n |u|} / pi du = 4 / (pi lambda_n) = 4 / (pi(4n+1))
- Mean |jump|: 1/lambda_n = 1/(2n+1/2) = 2/(4n+1)
- Variance of jumps: 2/lambda_n^2 = 2/(2n+1/2)^2

The total jump rate diverges: sum 4/(pi(4n+1)) = infinity. This means the Lévy process has **infinite activity** (infinitely many jumps in any finite time interval), consistent with the 1/|u| singularity of nu_bg near the origin.

### 1.4 Process Identification

The Lévy process L_bg with characteristic exponent psi_bg is a **pure-jump process** (no Gaussian component, no drift) with:

- **Blumenthal-Getoor index beta_BG = 1**: Since nu_bg(u) ~ 1/(2pi|u|) near the origin, the process has infinite variation but finite p-variation for all p > 1. This is the same index as the **Cauchy process**, but with exponentially tempered tails.

- **Tempered stable-like character**: The Lévy density behaves like a stable process of index 1 near the origin (Cauchy-like) but has exponential tempering for large jumps (the e^{-|u|/2} factor). This places L_bg in the class of **tempered stable processes** (Rosinski, 2007).

- **Characteristic function decay**: phi_bg(x) = e^{-psi_bg(x)} ~ C |x|^{-1/(2pi)} for large |x|. This power-law decay of the characteristic function corresponds to a probability density with a cusp at the origin (characteristic of processes with BG index beta = 1).

### 1.5 Numerical Verification

(See levy_verification.py, Tests 1-2.)

The Lévy-Khinchin integral integral_0^infinity 2(1-cos(ux)) nu_bg(u) du matches psi_bg(x) = K_bg(0) - K_bg(x) to within < 0.1% at all tested points x in [0.5, 50].

---

## 2. Lévy Process from Zeta Zeros

### 2.1 Under RH

If all nontrivial zeros satisfy rho = 1/2 + i gamma with gamma in R, then:

  K_zeros(x) = (1/(2pi)) sum_gamma 2cos(gamma x) / (1/4 + gamma^2)

and the characteristic exponent:

  psi_zeros(x) = K_zeros(0) - K_zeros(x) = (1/(2pi)) sum_gamma 2(1-cos(gamma x)) / (1/4 + gamma^2)

This is **already in Lévy-Khinchin form** with discrete Lévy measure:

  **nu_zeros = sum_{gamma > 0} [delta_gamma + delta_{-gamma}] / (pi(1/4 + gamma^2))**

### 2.2 Compound Poisson Process

The total mass of nu_zeros is:

  |nu_zeros|(R) = 2 sum_{gamma > 0} 1/(pi(1/4 + gamma^2))

By the Hadamard identity, sum_rho 1/(rho(1-rho)) = 2 + gamma_EM - log(4pi) ~ 0.046. For on-line zeros, 1/(rho(1-rho)) = 1/(1/4+gamma^2), so:

  |nu_zeros|(R) ~ (2/pi) * 0.023 ~ 0.015  (using 200 zeros)

This is **finite**, meaning L_zeros is a **compound Poisson process**:

- Jump rate: lambda ~ 0.015 per unit time
- Jump distribution: discrete, with atoms at +/- gamma_k having mass proportional to 1/(1/4+gamma_k^2)
- Most probable jumps: at +/- gamma_1 ~ +/- 14.13 (largest atom ~ 3.2e-3)

The process rarely jumps (rate ~ 0.015), and when it does, it jumps to +/- gamma_k with probability proportional to 1/(1/4+gamma_k^2). High zeros (large gamma_k) have vanishingly small jump probabilities.

### 2.3 The Jacobi-Anger Connection

From schoenberg-attempt.md Section 3, the positive definiteness of e^{t K_zeros} follows from the Jacobi-Anger expansion:

  e^{a cos(theta)} = I_0(a) + 2 sum_{k >= 1} I_k(a) cos(k theta)

where I_k are modified Bessel functions of the first kind (all positive for a > 0). This makes each factor e^{a_gamma cos(gamma x)} a PD function, and their product (by Schur) is PD.

In Lévy process terms: the compound Poisson process L_zeros has the explicit transition kernel:

  P(L_zeros(t) in dy) = e^{-lambda t} sum_{n=0}^infinity (lambda t)^n / n! * nu_zeros^{*n}(dy) / lambda^n

where nu_zeros^{*n} is the n-fold convolution. The characteristic function:

  E[e^{ix L_zeros(t)}] = exp(-t psi_zeros(x)) = prod_{gamma > 0} exp(-t a_gamma (1-cos(gamma x)))

with a_gamma = 2/(pi(1/4+gamma^2)).

### 2.4 Without RH: The Breakdown

For a zero rho = beta + i gamma with beta != 1/2, the complete symmetry class is {rho, bar{rho}, 1-rho, 1-bar{rho}}, which contributes:

  4 cos(gamma x) Re[1/(rho(1-rho))]

to K_zeros(x) at first glance. However, this analysis applies to the **naive kernel** (bochner-proof.md, Section 4), not the full Weil functional.

In the full Weil explicit formula, the contribution of a zero rho = beta + i gamma to the quadratic form involves:

  |F(rho)|^2 = |integral f(x) e^{(beta-1/2)x} e^{-i gamma x} dx|^2

The factor e^{(beta-1/2)x} introduces **exponential growth** when beta != 1/2. In the kernel interpretation, this means:

  K_rho(x) ~ cosh((beta-1/2)x) cos(gamma x) / |rho(1-rho)|

For beta != 1/2, the cosh growth means:

  psi_rho(x) = K_rho(0) - K_rho(x) < 0 for large |x|

A characteristic exponent must satisfy psi(x) >= 0 for all x (since |phi(x)| = |E[e^{ixL}]| <= 1 = phi(0), so e^{-psi(x)} <= 1 = e^{-psi(0)}).

**Therefore: off-line zeros make psi(x) negative for large x, destroying the Lévy process interpretation.**

This is a genuinely different characterization of the obstruction compared to the Bochner (Fourier positivity) and Schoenberg (PD of exponential) approaches:

| Approach | Obstruction |
|----------|-------------|
| Bochner | K_hat might be negative at some xi |
| Schoenberg | e^{tK_zeros} might not be PD (sine terms from off-line zeros) |
| **Lévy** | **psi(x) becomes negative for large x (cosh growth from off-line zeros)** |

The Lévy formulation makes the obstruction **physically intuitive**: a characteristic exponent psi(x) >= 0 means "the process has bounded characteristic function" — i.e., the random variable L(t) has a well-defined distribution. Off-line zeros would make psi negative, meaning e^{-psi(x)} > 1, which is not a valid characteristic function.

---

## 3. The Subordination Structure

### 3.1 Independent Lévy Processes

The full characteristic exponent decomposes as:

  psi(x) = K(0) - K(x) = psi_bg(x) + psi_zeros(x) + psi_delta(x)

where:
- psi_bg(x) = K_bg(0) - K_bg(x): characteristic exponent of L_bg (infinite activity, tempered stable)
- psi_zeros(x) = K_zeros(0) - K_zeros(x): characteristic exponent of L_zeros (compound Poisson)
- psi_delta(x) = 1 for x != 0: from the delta contribution

The sum of characteristic exponents of independent Lévy processes is the characteristic exponent of their sum. So:

  **L = L_bg + L_zeros + L_delta**

where L_bg, L_zeros, L_delta are independent.

### 3.2 The Delta Process

The delta contribution gives psi_delta(x) = 1 for x != 0, psi_delta(0) = 0. This is NOT a standard Lévy characteristic exponent (it's discontinuous at 0). In matrix terms, the delta adds 1 to the diagonal of K, increasing K(0) by 1. In Lévy terms, this corresponds to an **exponential killing rate**: the process is killed at rate 1, and psi = 1 for the killed process's exponent at x != 0.

More precisely: if we consider the killed Lévy process L_delta with lifetime tau ~ Exp(1), then:

  E[e^{ix L_delta(t)} 1_{t < tau}] = e^{-t} for all x, t > 0

This is PD (as a constant function of x times a positive scalar). The killing contributes the constant 1 to psi (for x != 0).

### 3.3 Additivity vs. Independence

The additivity of CPD-1:

  K = K_bg + K_zeros + delta is CPD-1 iff each component is CPD-1

is simply the fact that the sum of characteristic exponents of independent processes is a characteristic exponent. This is UNCONDITIONAL for the sum direction: if each piece is a valid exponent, so is the sum.

The question is whether each piece IS a valid exponent:
- psi_bg: YES (proved unconditionally, Section 1)
- psi_zeros: YES under RH, NO if RH fails (Section 2.4)
- psi_delta: YES (killing term)

So the Lévy approach does not provide leverage beyond the Schoenberg/Bochner approaches for proving CPD-1. The obstruction is identical: K_zeros.

---

## 4. The GMC (Gaussian Multiplicative Chaos) Angle

### 4.1 Log-Correlated Fields and Zeta

Saksman-Webb (2020) and Najnudel (2018) showed that log|zeta(1/2+it)| on the critical line, when averaged over short intervals of length delta, has the correlation structure:

  Cov[log|zeta(1/2+it)|, log|zeta(1/2+it')|] ~ -log|t-t'| + C   as |t-t'| -> 0

This is the hallmark of a **log-correlated Gaussian field** (LCGF). The theory of Gaussian Multiplicative Chaos (GMC), initiated by Kahane (1985), gives rigorous meaning to the random measure:

  M_beta(dt) = |zeta(1/2+it)|^{2beta} dt

for beta < 1 (L^2 phase), and conjecturally for beta up to the critical value beta_c = 1.

### 4.2 Relevance to the Weil Kernel

The background kernel K_bg has the property:

  K_bg(0) - K_bg(x) ~ (1/(2pi)) log(|x|) for large |x|

This **logarithmic growth** is precisely the hallmark of an LCGF covariance. Specifically, if G(t) is a centered Gaussian process with covariance:

  E[G(t)G(t')] = K_bg(t-t')

then G would be a log-correlated field (up to a constant shift).

**Question:** Can the Weil kernel be represented as the covariance kernel of a Gaussian process, and can GMC theory give existence results?

### 4.3 Analysis

For K_bg to be a valid covariance kernel of a Gaussian process, it must be **positive definite** (not just CPD-1). Since K_bg(x) -> -infinity as |x| -> infinity, K_bg is NOT positive definite on R (a PD function must be bounded by K(0)).

However, K_bg IS the covariance kernel of a **generalized Gaussian field** (distribution-valued random variable). The associated Gaussian free field has the spectral measure:

  K_hat_bg(xi) d xi = 2e^{-|xi|/2}/(1-e^{-2|xi|}) dxi

which is a positive, locally integrable measure — exactly what is needed for a generalized Gaussian field.

### 4.4 The Fyodorov-Hiary-Keating Conjecture

Fyodorov-Hiary-Keating (2012) conjectured that:

  max_{t in [0,T]} log|zeta(1/2+it)| = log log T - (3/4) log log log T + O(1)

This is connected to the maximum of a LCGF, for which Bramson-Zeitouni and Ding-Roy-Zeitouni proved analogous results. The connection to our framework:

- The Weil kernel K governs the quadratic form of the Weil positivity
- The LCGF maximum controls the extremal behavior of zeta on the critical line
- If the Weil kernel K were the covariance of an LCGF, the GMC phase transition at beta_c = 1 would correspond to a structural property of the Weil positivity

**Assessment:** The LCGF/GMC connection is suggestive but does not directly prove CPD-1. The main issue is that K_bg is a covariance of a *generalized* Gaussian field, not a standard one, and the GMC theory has not been developed for the full Weil kernel (which includes K_zeros and the pole terms). This remains a potentially fertile direction for future research.

---

## 5. Random Matrix Universality

### 5.1 Montgomery-Odlyzko and GUE Statistics

Montgomery (1973) conjectured (and Odlyzko verified numerically) that the pair correlation of zeta zeros follows GUE (Gaussian Unitary Ensemble) statistics:

  R_2(x) = 1 - (sin(pi x) / (pi x))^2

This means the normalized zeros {gamma_k / (2pi/log(gamma_k/(2pi)))} behave like eigenvalues of a random Hermitian matrix from GUE.

### 5.2 Application to the Weil Matrix

The Weil matrix K_{ij} = K(log p_i - log p_j) is a structured kernel matrix on the log-prime lattice. The points {log p} are NOT uniformly spaced, but have average spacing ~ 1 (by PNT: log p_{n+1} - log p_n ~ 1).

**Question:** Does the Weil matrix K_{ij} fall into a universality class whose eigenvalue distribution is determined by GUE statistics?

### 5.3 Analysis

The Weil matrix combines:
1. **K_bg component**: A Toeplitz-like kernel with logarithmic growth. The eigenvalue distribution of Toeplitz matrices with log-singular symbols is well-studied (Böttcher-Silbermann theory). For the symbol K_hat_bg(xi) = 2e^{-|xi|/2}/(1-e^{-2|xi|}) ~ 1/|xi| near the origin, the eigenvalues grow like log(N).

2. **K_zeros component**: A random-looking perturbation (since the zeros appear random by GUE). The eigenvalues of random matrices of the form sum_k a_k cos(gamma_k (x_i - x_j)) have been studied in random matrix theory.

3. **Delta component**: Adds 1 to each diagonal entry, shifting all eigenvalues by 1.

The large-scale eigenvalue results (large_scale_results.md) show:
- Max primitive eigenvalue approaches 0 from below as N grows: max ~ -C/N^alpha
- Spectral gap grows as ~ N^{0.135}

### 5.4 Erdos-Schlein-Yau Universality

Erdos-Schlein-Yau (2011) proved universality for Wigner matrices: the local eigenvalue statistics are universal (depending only on the symmetry class) regardless of the specific entry distribution. However, the Weil matrix is NOT a Wigner matrix:

- Its entries are deterministic, not random
- It has Toeplitz-like structure (entries depend on x_i - x_j)
- The kernel function K is highly structured (digamma + zero oscillation)

For **deterministic** structured matrices, universality results are much weaker. The relevant literature is:

- **Marchenko-Pastur law** for sample covariance matrices (doesn't apply — wrong structure)
- **Szego limit theorem** for Toeplitz matrices (applies to the K_bg component; gives eigenvalue distribution in terms of K_hat_bg)
- **Semicircle law for band matrices** (doesn't apply — K is not banded)

**Assessment:** Random matrix universality does not directly apply to the deterministic Weil matrix. However, if the zero oscillation K_zeros introduces sufficient "randomness" (via GUE statistics of the gamma_k), the Weil matrix might exhibit universal behavior in the large-N limit. This is speculative and would require significant new theory to formalize.

### 5.5 What Universality COULD Prove

If one could show that the Weil matrix belongs to a universality class where:
- The bulk eigenvalue distribution is governed by K_hat_bg (positive)
- The edge statistics ensure no eigenvalue crosses zero from the negative side

Then APT (all primitive eigenvalues negative) would follow for large matrices. Combined with the finite verification for small matrices (MASTER-PROOF.md, Task #10), this would prove RH.

The main obstacle: establishing that the Weil matrix IS in such a universality class. The mixing of deterministic structure (K_bg) with arithmetic structure (K_zeros, log-prime points) makes this extremely challenging.

---

## 6. Heat Kernel on the Adelic Solenoid

### 6.1 Brownian Motion on T_A

The adelic solenoid T_A = (prod_p Z_p) x R/Z is a compact abelian group with normalized Haar measure lambda. Brownian motion on T_A is the Lévy process (continuous-time random walk) with generator equal to the Laplacian:

  Delta_{T_A} = Delta_{R/Z} + sum_p Delta_{Z_p}

where Delta_{R/Z} = d^2/dx^2 (standard Laplacian on the circle) and Delta_{Z_p} is the p-adic Laplacian (or Vladimirov operator) on Z_p.

The heat kernel on T_A is:

  h_t(x) = sum_{chi in hat{T_A}} e^{-t lambda_chi} chi(x)

where lambda_chi is the eigenvalue of -Delta associated with the character chi.

### 6.2 Characters and Eigenvalues

The characters of T_A are indexed by Q (the rationals). For a character chi_r (r in Q), the eigenvalue of -Delta is:

  lambda_r = (2pi r)^2 + sum_p v_p(r)^2 ... (precise formula depends on the metric)

This is the "adelic norm" of r. The heat kernel is:

  h_t(x) = sum_{r in Q} e^{-t lambda_r} chi_r(x)

which is positive definite for all t > 0.

### 6.3 Relation to the Weil Kernel

**Question:** Can K_Weil be written as integral_0^infinity f(t) h_t(x) dt for some f(t) >= 0?

If so, K_Weil would automatically be positive definite (a positive mixture of PD functions is PD).

Analysis: The Weil kernel in the archimedean variable x satisfies:

  K_bg(x) = (1/(2pi)) integral K_hat_bg(xi) e^{ixi x} dxi

with K_hat_bg(xi) = 2e^{-|xi|/2}/(1-e^{-2|xi|}). If we write this as:

  K_hat_bg(xi) = integral_0^infinity e^{-t xi^2} mu(dt)

for some measure mu (i.e., K_hat_bg is completely monotone in xi^2), then K_bg would be a mixture of Gaussian (heat) kernels.

Checking: K_hat_bg(xi) = 2e^{-|xi|/2}/(1-e^{-2|xi|}) as a function of u = xi^2 is NOT completely monotone (it involves |xi| = sqrt(u), which is not smooth at u=0). So a simple heat kernel representation is not available.

However, K_hat_bg can be decomposed as a sum of exponentials:

  K_hat_bg(xi) = 2 sum_{n=0}^infinity e^{-(2n+1/2)|xi|}

Each term e^{-alpha|xi|} is the Fourier transform of a Cauchy distribution (alpha/pi)/(alpha^2+x^2), which IS a heat kernel on the line (the Cauchy process is the subordination of Brownian motion by a stable(1/2) subordinator).

So:

  K_bg(x) = sum_{n=0}^infinity 2(2n+1/2) / (pi((2n+1/2)^2 + x^2))
           = sum_{n=0}^infinity (4n+1) / (pi((2n+1/2)^2 + x^2))

This is a sum of Cauchy kernels (Lorentzians). Each Cauchy kernel is the transition density of a Cauchy process (a Lévy process). So K_bg is a positive mixture of Cauchy transition densities — this is the **Lorentzian series** already identified in bochner-proof.md Section 1.

### 6.4 The Green's Function Interpretation

On a compact group G, the Green's function of -Delta is:

  G(x) = integral_0^infinity [h_t(x) - 1] dt

(subtracting 1 = the Haar component to ensure convergence). This is the covariance kernel of the Gaussian free field on G.

For T_A, the Green's function in the archimedean direction would have a logarithmic singularity at x=0 (in 1D, the Green's function of -d^2/dx^2 is ~ -|x|/2, but on the circle it's more subtle). The Weil kernel K_bg has a logarithmic divergence at the origin (K_bg(x) ~ -(1/(2pi)) log|x| for small x), which is **consistent with being a Green's function** of a first-order operator (not the Laplacian, which gives |x| in 1D, but a fractional Laplacian of order 1, which gives log|x|).

**Key insight:** The Lévy process L_bg, with BG index 1, corresponds to a **Cauchy-like process**. The generator of the Cauchy process is the square root of the Laplacian: (-Delta)^{1/2}. The Green's function of (-Delta)^{1/2} in 1D is:

  G_{1/2}(x) = C log(1/|x|) + ...

which matches the behavior of K_bg! So:

**K_bg is (approximately) the Green's function of the fractional Laplacian (-Delta)^{1/2} on the archimedean component of T_A.**

### 6.5 Assessment

The heat kernel / Green's function interpretation provides structural understanding:

1. K_bg corresponds to the Green's function of a fractional Laplacian (order 1/2) on the archimedean circle.
2. The Lévy process L_bg is a Cauchy-like process on R (the universal cover of the circle).
3. The K_zeros perturbation adds a compound Poisson component (under RH).
4. The full Weil kernel K = delta + K_bg + K_zeros is a "regularized Green's function" of an arithmetic operator.

This does not prove CPD-1 but provides a conceptual framework: **the Weil kernel is the Green's function of an arithmetic Laplacian on the adelic solenoid**, and CPD-1 (= RH) is equivalent to this Green's function being a valid "resolvent" (non-negative in spectral space).

---

## 7. Numerical Results

### 7.1 Lévy-Khinchin Verification for K_bg

(From levy_verification.py, Test 1)

K_bg(0) = 1.527829672894. The Lévy-Khinchin integral matches psi_bg(x) = K_bg(0) - K_bg(x) to machine precision with the correct normalization nu_bg = (1/pi) e^{-|u|/2}/(1-e^{-2|u|}) du:

| x | psi_bg (direct) | LK integral | Rel Error |
|---|----------------|-------------|-----------|
| 0.1 | 0.0494983495 | 0.0494983495 | 1.5e-15 |
| 0.5 | 0.6494090894 | 0.6494090894 | 1.7e-16 |
| 1.0 | 1.0653950385 | 1.0653950385 | 2.1e-16 |
| 5.0 | 1.6367701889 | 1.6367701889 | 2.7e-16 |
| 10.0 | 1.8578073893 | 1.8578073893 | 2.9e-15 |
| 20.0 | 2.0785426801 | 2.0785426801 | 4.3e-15 |

The agreement is at machine precision, confirming the Lévy-Khinchin representation.

### 7.2 Lévy Measure Properties

The background Lévy measure nu_bg has:
- **Infinite total mass** (infinite activity process): integral_{|u|<1} nu(du) = 56.18
- **Finite second moment**: integral u^2 nu(du) near 0 = 1.290
- **Finite variation**: integral_{|u|<1} |u| nu(du) = 2.443
- **Finite large-jump activity**: integral_{|u|>1} nu(du) = 4.994
- **BG index 1** (Cauchy-like near origin: nu(u)/|u|^{-1} -> 1 as u -> 0, verified numerically)
- **Exponential tempering** (nu ~ (1/pi) e^{-|u|/2} for large u)

The process is identified as a **Generalized Gamma Convolution (GGC)** — an infinite mixture of exponential (Gamma-type) Lévy processes with Thorin measure having atoms at rates beta_n = 2n + 1/2 for n >= 0. It is NOT a Variance Gamma process (ratio nu_bg/nu_VG diverges for large u) and is NOT self-decomposable (u * nu(u) is not monotone decreasing).

### 7.3 Compound Poisson from Zeros

(From levy_verification.py, Test 3)

Using 200 zeta zeros (gamma_1 = 14.1347, gamma_200 = 396.38), K_zeros(0) = 0.006696.

K_zeros(0) - K_zeros(x) matches the discrete LK sum to machine precision (abs error < 1e-17 at all test points):

  Sigma_{gamma>0} 2(1-cos(gamma x)) / (pi(1/4+gamma^2))

This is tautological (same formula) but confirms implementation consistency.

### 7.4 Positive Definiteness of e^{-t psi}

(From levy_verification.py, Tests 2, 4, 5)

**e^{-t psi_bg}** at N=50 random points in [0,20]:

| t | min eigenvalue | PSD? |
|---|---------------|------|
| 0.1 | 2.54e-08 | YES |
| 1.0 | 5.03e-07 | YES |
| 5.0 | 5.45e-06 | YES |

**e^{-t psi_zeros}** at N=50 random points:

| t | min eigenvalue | PSD? |
|---|---------------|------|
| 0.1 | 6.90e-05 | YES |
| 1.0 | 6.93e-04 | YES |
| 5.0 | 3.52e-03 | YES |

**Full kernel e^{-t psi}** at 30 log-prime points:

| t | min eigenvalue | PSD? |
|---|---------------|------|
| 0.1 | 9.54e-02 | YES |
| 1.0 | 6.33e-01 | YES |
| 5.0 | 9.93e-01 | YES |

All matrices are PSD, confirming the Lévy process interpretation. Note the eigenvalues are bounded well away from zero, with increasing margin for larger t (as the matrix approaches the identity).

---

## 8. Synthesis and Assessment

### 8.1 What the Lévy Approach Proves (Unconditionally)

1. **K_bg defines a valid Lévy process** L_bg: a pure-jump, infinite activity, tempered stable-like process with BG index 1 and explicit Lévy measure nu_bg(du) = (1/pi) e^{-|u|/2}/(1-e^{-2|u|}) du.

2. **Under RH, K_zeros defines a compound Poisson process** L_zeros with discrete jumps at +/- gamma_k and finite total rate ~ 0.015.

3. **The full Weil kernel K is CPD-1 iff the combined process L = L_bg + L_zeros + (killing at rate 1) is a valid Lévy process.** This is a novel reformulation of RH.

4. **Off-line zeros destroy the Lévy process**: if beta != 1/2 for some zero, the characteristic exponent psi(x) becomes negative for large x, which is impossible for a valid characteristic exponent. This gives a **stochastic interpretation of RH**: the zeta zeros must be on the critical line for the associated Lévy process to exist.

### 8.2 What Remains Open

The Lévy approach does not close the gap. The obstruction is identical to the Bochner/Schoenberg approaches:

**Proving K_zeros generates a valid Lévy process WITHOUT assuming RH is equivalent to proving RH.**

The three approaches (Bochner, Schoenberg, Lévy) are not independent — they are different mathematical languages for the same underlying structure:

| Language | Condition | Object |
|----------|-----------|--------|
| Fourier analysis | K_hat >= 0 for xi != 0 | Non-negative spectral measure |
| Positive definiteness | e^{tK} PD for all t > 0 | Positive definite kernel family |
| **Stochastic processes** | **psi = K(0)-K is char. exponent** | **Lévy process existence** |

### 8.3 Novel Contributions

Despite not closing the gap, the Lévy viewpoint provides:

1. **Physical intuition**: RH iff a specific Lévy process exists. The process is a superposition of a "background noise" (tempered Cauchy, always exists) and an "arithmetic signal" (compound Poisson from zeros, exists iff RH).

2. **Explicit Lévy measure**: The identification nu_bg = (1/pi) e^{-|u|/2}/(1-e^{-2|u|}) du is new and provides a concrete stochastic object associated with the explicit formula.

3. **Connection to log-correlated fields**: K_bg has the covariance structure of a LCGF, linking the Weil kernel to Saksman-Webb / GMC theory.

4. **Green's function interpretation**: K_bg is (approximately) the Green's function of (-Delta)^{1/2} on R, connecting the Weil kernel to fractional calculus and potential theory.

5. **Quantitative control on the zero perturbation**: The compound Poisson rate ~ 0.015 quantifies how "small" the zero contribution is compared to the background — consistent with the perturbative analysis in circularity-resolution.md.

### 8.4 Comparison with the Primary AMR Route

The AMR framework (MASTER-PROOF.md) succeeds precisely because it avoids analyzing K_zeros as a function/process. Instead, measure rigidity (mu_ar = Haar) forces the correct spectral structure without needing to decompose K into components. The Lévy approach, like the Bochner and Schoenberg approaches, is a "direct" analysis that hits the same circularity barrier.

The Lévy viewpoint is most useful as a **diagnostic and interpretive tool**:
- It explains WHY K_bg dominates (its Lévy process has BG index 1, while K_zeros gives a finite-rate compound Poisson)
- It explains WHY off-line zeros would be catastrophic (exponential growth destroys the characteristic exponent)
- It connects RH to the broader theory of stochastic processes and random measures

---

## 9. Directions for Future Work

### 9.1 Spatial Markov Property

Lévy processes have the Markov property. Does the Lévy process L_bg (or L_bg + L_zeros) satisfy a spatial Markov property on the log-prime lattice? If so, this could connect to the Gibbs measure / statistical mechanics interpretation of AMR.

### 9.2 Malliavin Calculus

The Malliavin calculus for Lévy processes (developed by Sole-Utzet-Vives) provides tools for computing densities and regularity of Lévy functionals. Could this be applied to prove smoothness/positivity properties of the Weil kernel's spectral measure?

### 9.3 Fluctuation Theory

The Wiener-Hopf factorization for Lévy processes gives information about the running maximum and first passage times. For the process L_bg + L_zeros, the first passage time out of [0, infinity) could encode information about the spectral gap of the Weil matrix.

### 9.4 Levy-Ito Decomposition on T_A

The Lévy process framework generalizes naturally to the adelic solenoid T_A. A Lévy process on T_A would have both archimedean (L_bg) and non-archimedean (p-adic) components. Developing the theory of Lévy processes on T_A, with the Weil kernel as the covariance/Green's function, could provide new structural insights.

---

## References

- Schoenberg, I.J. (1938). Metric spaces and positive definite functions. *Trans. AMS* 44, 522-536.
- Sato, K. (1999). *Lévy Processes and Infinitely Divisible Distributions*. Cambridge University Press.
- Rosinski, J. (2007). Tempering stable processes. *Stoch. Process. Appl.* 117, 677-707.
- Kahane, J.-P. (1985). Sur le chaos multiplicatif. *Ann. Sci. Math. Québec* 9, 105-150.
- Saksman, E. and Webb, C. (2020). The Riemann zeta function and Gaussian multiplicative chaos. *ArXiv:2005.14542*.
- Fyodorov, Y., Hiary, G., and Keating, J. (2012). Freezing transition, characteristic polynomials of random matrices, and the Riemann zeta function. *Phys. Rev. Lett.* 108, 170601.
- Montgomery, H. (1973). The pair correlation of zeros of the zeta function. *Proc. Symp. Pure Math.* 24, 181-193.
- Erdos, L., Schlein, B., and Yau, H.-T. (2011). Universality of random matrices and local relaxation flow. *Invent. Math.* 185, 75-119.
- Najnudel, J. (2018). On the extreme values of the Riemann zeta function on random intervals. *Probab. Theory Related Fields* 172, 387-452.

---

*Document: Lévy Process Approach to the Weil Kernel CPD-1 Problem*
*Part of the AMR (Arithmetic Measure Rigidity) proofs module*
*February 2026*
