# Explicit Derivation of Cross-Terms in the Weil Positivity Criterion

## Overview

This document derives the explicit form of the cross-terms arising in Weil's positivity criterion for the Riemann Hypothesis. The central object is a bilinear form on test functions whose positive-definiteness is equivalent to RH. We decompose this form, identify its integral kernel in terms of the von Mangoldt function, compute its spectral representation, analyze positive-definiteness, and connect to the pair correlation of zeros via the explicit formula.

---

## 1. The Full Bilinear Form as an Integral Operator

### 1.1 Setup and Conventions

We work in additive coordinates on R. Define:

- **Test functions:** f in C_c^infty(R) (smooth, compactly supported)
- **Involution:** f~(x) = f-bar(-x)
- **Convolution:** (f * f~)(x) = integral f(t) f-bar(t - x) dt
- **Mellin-adapted transform:** F(s) = integral f(t) e^{st} dt for s in C

On the critical line s = 1/2 + i*tau, F(1/2 + i*tau) is the Fourier transform of t -> f(t)e^{t/2} evaluated at tau. We have:

  (f * f~)^transform(s) = |F(s)|^2  when Re(s) = 1/2.

### 1.2 Weil's Explicit Formula as a Distribution

The Weil distribution W acts on test functions g in C_c^infty(R) by:

**W(g) = g-hat(0) + g-hat(1) - sum_p sum_{m=1}^infty (log p / p^{m/2}) * 2 Re[g(m log p)] + Omega(g)**

where:
- g-hat(s) = integral g(x) e^{sx} dx (Mellin-type transform)
- The sum is over all primes p and positive integers m
- Omega(g) is the archimedean contribution:

  Omega(g) = integral_0^infty [2g(0)/(e^x - 1) - 2g(x)e^{-x/2}/(1 - e^{-2x}) - 2g(-x)e^{-x/2}/(1 - e^{-2x})] dx
            + g(0)(log pi - gamma/2)

  (arising from the digamma/Gamma-factor terms in the functional equation of zeta).

**Weil's Criterion:** RH holds if and only if W(f * f~) >= 0 for all f in C_c^infty(R).

### 1.3 Substituting g = f * f~

Setting g = f * f~ and evaluating:

  g(m log p) = integral f(t) f-bar(t - m log p) dt

  g-hat(s) = |F(s)|^2

Therefore:

**W(f * f~) = |F(0)|^2 + |F(1)|^2 - P(f, f) + Omega(f * f~)**

where the **prime bilinear form** is:

  P(f, f) = sum_p sum_{m=1}^infty (log p / p^{m/2}) * 2 Re[integral f(t) f-bar(t - m log p) dt]

### 1.4 The Sesquilinear Extension

Polarizing to obtain the full sesquilinear (Hermitian) form:

  P(f, g) = sum_p sum_{m=1}^infty (log p / p^{m/2}) * [integral f(t) g-bar(t - m log p) dt + integral g(t) f-bar(t - m log p) dt]

Equivalently, incorporating both positive and negative shifts:

  P(f, g) = integral integral f(t) g-bar(s) K_prime(t - s) dt ds

where K_prime is the **prime kernel**. The bilinear form P is thus an integral operator with distributional kernel K_prime.

---

## 2. The Kernel in Terms of the Von Mangoldt Function

### 2.1 Definition of the Prime Kernel

**Definition.** The prime kernel is the distribution:

  K_prime(x) = sum_{n=2}^infty (Lambda(n) / sqrt(n)) * [delta(x - log n) + delta(x + log n)]

where Lambda(n) is the von Mangoldt function:

  Lambda(n) = log p   if n = p^m for some prime p and integer m >= 1
  Lambda(n) = 0       otherwise

Grouping by prime:

  K_prime(x) = sum_p sum_{m=1}^infty (log p / p^{m/2}) * [delta(x - m log p) + delta(x + m log p)]

**Verification.** Substituting into the integral operator form:

  integral integral f(t) f-bar(s) K_prime(t - s) dt ds

  = sum_{n>=2} (Lambda(n)/sqrt(n)) * [integral f(t) f-bar(t - log n) dt + integral f(t) f-bar(t + log n) dt]

  = sum_{n>=2} (Lambda(n)/sqrt(n)) * 2 Re[integral f(t) f-bar(t - log n) dt]

  = P(f, f)

as required. The symmetrization delta(x - log n) + delta(x + log n) ensures the kernel is even: K_prime(-x) = K_prime(x), which guarantees the sesquilinear form is Hermitian.

### 2.2 The Kernel as a Measure

As a Radon measure on R:

  K_prime = sum_{n=2}^infty (Lambda(n)/sqrt(n)) * (delta_{log n} + delta_{-log n})

**Support:** The support is the discrete set

  supp(K_prime) = {+/- log n : n >= 2, Lambda(n) > 0} = {+/- m log p : p prime, m >= 1}

These are the logarithms of prime powers, symmetrized about the origin. The set is countable, with accumulation points only at +/- infinity.

**Coefficient decay:** The weight assigned to the atom at x = m log p is log p / p^{m/2}, which decays exponentially in m (for fixed p) and polynomially in p (for fixed m). Specifically:

  For m = 1: weight at log p is (log p)/sqrt(p), which decreases as p -> infinity.
  For m >= 2: weight at m log p is (log p)/p^{m/2} << 1/p, negligible for all but the smallest primes.

### 2.3 Convergence of the Bilinear Form

The total mass of the kernel (summing absolute values of coefficients) is:

  2 * sum_{n=2}^infty Lambda(n)/sqrt(n)

This sum **diverges**: by the prime number theorem, sum_{n <= X} Lambda(n)/sqrt(n) ~ 2*sqrt(X) as X -> infinity.

However, P(f, f) converges for each f in C_c^infty(R). If supp(f) subset [-T, T], then f(t) * f-bar(t - log n) = 0 whenever log n > 2T, i.e., n > e^{2T}. The sum is therefore finite for each fixed f:

  P(f, f) = sum_{n=2}^{floor(e^{2T})} (Lambda(n)/sqrt(n)) * 2 Re[integral f(t) f-bar(t - log n) dt]

Moreover, for smooth f, the inner integrals decay rapidly as log n approaches the boundary of the support, providing additional effective convergence.

### 2.4 Decomposition by Prime

The kernel decomposes as a sum over individual primes:

  K_prime(x) = sum_p K_p(x)

where the **single-prime kernel** is:

  K_p(x) = sum_{m=1}^infty (log p / p^{m/2}) * [delta(x - m log p) + delta(x + m log p)]

K_p is supported on the lattice log p * Z \ {0}, with exponentially decaying coefficients.

The bilinear form decomposes correspondingly:

  P(f, f) = sum_p P_p(f, f)

where:

  P_p(f, f) = sum_{m=1}^infty (log p / p^{m/2}) * 2 Re[integral f(t) f-bar(t - m log p) dt]

**Critical observation:** Each term in the sum over n contributes independently to P(f, f). There are no "explicit" cross-terms P_{p,q} between distinct primes in this direct decomposition. The cross-interactions between primes are mediated entirely through the test function f.

---

## 3. Spectral Decomposition (Mellin/Fourier Transform) of the Kernel

### 3.1 Fourier Transform of K_prime

The distributional Fourier transform of K_prime is:

  K-hat_prime(tau) = integral K_prime(x) e^{-i*tau*x} dx

  = sum_{n=2}^infty (Lambda(n)/sqrt(n)) * [e^{-i*tau*log n} + e^{i*tau*log n}]

  = sum_{n=2}^infty (Lambda(n)/sqrt(n)) * 2 cos(tau * log n)

  = sum_{n=2}^infty Lambda(n) * [n^{-(1/2 + i*tau)} + n^{-(1/2 - i*tau)}]

Recognizing the Dirichlet series:

  sum_{n=1}^infty Lambda(n) * n^{-s} = -zeta'(s)/zeta(s)

we obtain:

**K-hat_prime(tau) = -zeta'(1/2 + i*tau)/zeta(1/2 + i*tau) - zeta'(1/2 - i*tau)/zeta(1/2 - i*tau)**

**= -2 Re[zeta'(1/2 + i*tau) / zeta(1/2 + i*tau)]**

This is valid in the sense of distributions (the Dirichlet series for -zeta'/zeta converges only for Re(s) > 1, but the distributional identity extends by analytic continuation and the theory of generalized functions).

### 3.2 The Bilinear Form in Spectral Space

By the Plancherel theorem, the bilinear form becomes:

  P(f, f) = (1/2*pi) * integral |F(1/2 + i*tau)|^2 * K-hat_prime(tau) d*tau

  = -(1/pi) * integral |F(1/2 + i*tau)|^2 * Re[zeta'(1/2 + i*tau)/zeta(1/2 + i*tau)] d*tau

where F(s) = integral f(t) e^{st} dt is the Mellin-adapted transform.

This spectral representation reveals P as a weighted integral of the **spectral density** |F(1/2 + i*tau)|^2 against the weight function -Re[zeta'/zeta(1/2 + i*tau)].

### 3.3 Singularity Structure of K-hat_prime

The function -zeta'/zeta(s) is meromorphic on C with:

- **Simple poles at non-trivial zeros** rho = beta + i*gamma of zeta(s), each with residue +1
- **Simple pole at s = 1** (from the pole of zeta), with residue -1
- **Simple poles at trivial zeros** s = -2, -4, -6, ..., each with residue +1

On the critical line s = 1/2 + i*tau, the behavior near a zero rho_0 = beta_0 + i*gamma_0 is:

  -zeta'(1/2 + i*tau)/zeta(1/2 + i*tau) ~ 1/((1/2 - beta_0) + i*(tau - gamma_0))

Taking the real part:

  Re[...] ~ (1/2 - beta_0) / [(1/2 - beta_0)^2 + (tau - gamma_0)^2]

**Case 1: Zero on the critical line** (beta_0 = 1/2): The real part contribution is 0/((tau - gamma_0)^2) = 0. Zeros on the critical line do NOT contribute to Re[-zeta'/zeta] on the critical line. They contribute only to the imaginary part.

**Case 2: Zero off the critical line** (beta_0 != 1/2): The contribution is a Lorentzian peak of height 1/(1/2 - beta_0) and width |1/2 - beta_0|. Since zeros come in pairs {rho, 1 - rho-bar}, the paired contributions are:

  (1/2 - beta_0)/[...] + (beta_0 - 1/2)/[...] = 0

These cancel exactly. Thus Re[-zeta'/zeta(1/2 + i*tau)] is independent of the locations of non-trivial zeros!

**Consequence:** K-hat_prime(tau) = -2 Re[zeta'/zeta(1/2 + i*tau)] is determined entirely by:
- The pole of zeta at s = 1
- The trivial zeros of zeta
- The Gamma-factor contributions

It does NOT directly encode individual zero locations. The positivity of W(f * f~) is a more subtle condition than pointwise non-negativity of K-hat_prime.

### 3.4 Spectral Transform of Individual Prime Kernels

For a single prime p:

  K-hat_p(tau) = sum_{m=1}^infty (log p / p^{m/2}) * 2 cos(m * tau * log p)

  = 2 * log p * Re[sum_{m=1}^infty (p^{-1/2} * e^{i*tau*log p})^m]

  = 2 * log p * Re[p^{-(1/2 - i*tau)} / (1 - p^{-(1/2 - i*tau)})]

  = 2 * log p * Re[1 / (p^{1/2 - i*tau} - 1)]

Setting z = p^{1/2 - i*tau} = p^{1/2} * e^{-i*tau*log p}:

  K-hat_p(tau) = 2 * log p * Re[1/(z - 1)]
               = 2 * log p * (Re(z) - 1) / |z - 1|^2
               = 2 * log p * (p^{1/2} * cos(tau * log p) - 1) / (p + 1 - 2*p^{1/2} * cos(tau * log p))

**Properties:**
- K-hat_p(0) = 2*log p / (p^{1/2} - 1) > 0
- K-hat_p(tau) oscillates with quasi-period 2*pi / log p
- K-hat_p(tau) < 0 when cos(tau * log p) < p^{-1/2}, i.e., when the numerator is negative

The full spectral kernel satisfies the Euler product identity:

  K-hat_prime(tau) = sum_p K-hat_p(tau) = -2 Re[zeta'/zeta(1/2 + i*tau)]

which is the additive decomposition of the logarithmic derivative over the Euler product.

### 3.5 Explicit Values

At tau = 0:

  K-hat_prime(0) = -2 * zeta'(1/2) / zeta(1/2)

Numerically: zeta(1/2) ~ -1.4604, zeta'(1/2) ~ -3.9228, so:

  K-hat_prime(0) ~ -2 * (-3.9228) / (-1.4604) ~ -5.373

This is **negative**, immediately showing that K_prime is NOT positive-definite.

For large tau, by the asymptotics of zeta'/zeta on the critical line:

  Re[-zeta'/zeta(1/2 + i*tau)] ~ (1/2) * log(|tau|/(2*pi)) + O(1/|tau|)

which is positive and growing. So K-hat_prime(tau) > 0 for |tau| sufficiently large.

---

## 4. Positive-Definiteness vs. Indefiniteness of the Kernel

### 4.1 The Criterion for Positive-Definiteness

A translation-invariant kernel K(x - y) on R defines a positive semi-definite bilinear form if and only if its Fourier transform satisfies:

  K-hat(tau) >= 0  for all tau in R

### 4.2 The Prime Kernel is Indefinite

**Theorem.** The prime kernel K_prime is NOT positive-definite.

**Proof.** By Section 3.5, K-hat_prime(0) ~ -5.373 < 0. Since the Fourier transform takes a negative value, the kernel fails to be positive-definite.

More precisely, there exist test functions f for which P(f, f) < 0. Choose f such that F(1/2 + i*tau) concentrates near tau = 0; then:

  P(f, f) ~ (1/2*pi) * |F(1/2)|^2 * K-hat_prime(0) < 0

### 4.3 The Prime Kernel is Also Not Negative-Definite

For large tau, K-hat_prime(tau) ~ log(|tau|/(2*pi)) > 0, so there exist f with P(f, f) > 0 (choose F concentrating at large tau).

**Conclusion:** The prime bilinear form P(f, f) is genuinely indefinite — it takes both positive and negative values depending on f.

### 4.4 The Complete Weil Kernel

The full Weil functional involves:

  W(f * f~) = |F(0)|^2 + |F(1)|^2 - P(f, f) + Omega(f * f~)

In spectral form:

  W(f * f~) = (1/2*pi) * integral |F(1/2 + i*tau)|^2 * W-hat(tau) d*tau  + |F(0)|^2 + |F(1)|^2

where:

  W-hat(tau) = -K-hat_prime(tau) + Omega-hat(tau)
             = 2 Re[zeta'/zeta(1/2 + i*tau)] + (archimedean spectral function)

By the explicit formula, the complete Weil kernel satisfies:

  W(f * f~) = sum_rho |F(rho)|^2

**Under RH:** All zeros have Re(rho) = 1/2, so F(rho) = F(1/2 + i*gamma) for real gamma, and sum_rho |F(rho)|^2 >= 0.

**Without RH:** If rho_0 = beta_0 + i*gamma_0 with beta_0 != 1/2, one can construct f with F concentrated near rho_0 and 1 - rho_0-bar to make the sum negative.

**Conclusion:** Positive-definiteness of the complete Weil kernel is equivalent to RH. The prime kernel alone is indefinite; it is the precise cancellation between P, the pole terms, and the archimedean contribution that yields (conditional) positivity.

### 4.5 The Positivity Budget

For the Weil functional to be non-negative, we need:

  P(f, f) <= |F(0)|^2 + |F(1)|^2 + Omega(f * f~)

The LHS is the prime contribution (indefinite). The RHS has three parts:

1. **Pole terms** |F(0)|^2 + |F(1)|^2: A rank-2 positive-semidefinite form. This provides positivity concentrated at two specific "frequencies" (s = 0 and s = 1).

2. **Archimedean term** Omega(f * f~): Controlled by the Gamma-factor. In spectral form:

    Omega-hat(tau) = Re[Psi(1/4 + i*tau/2)] + log pi - gamma/2

   where Psi = Gamma'/Gamma is the digamma function. For large tau: Omega-hat(tau) ~ (1/2)*log(|tau|/2) > 0. For small tau: Omega-hat(0) = Psi(1/4) + log pi - gamma/2 ~ -2.91.

3. **The interplay:** The pole terms compensate for the negativity of P near tau = 0, while the archimedean term compensates for large tau.

---

## 5. Connection to the Pair Correlation of Zeros via the Explicit Formula

### 5.1 The Fundamental Spectral Identity

The explicit formula for L-functions, applied to g = f * f~, yields the identity:

**sum_rho |F(rho)|^2 = |F(0)|^2 + |F(1)|^2 - P(f, f) + Omega(f * f~)**

This is the **spectral identity**: the sum over zeros (spectral side) equals the sum over primes plus poles plus archimedean (arithmetic side).

Rearranging:

**sum_rho |F(rho)|^2 = (pole terms) - (prime terms) + (archimedean terms)**

### 5.2 Derivation from the Explicit Formula

Starting from the Weil explicit formula for an even Schwartz function h:

  sum_rho h-hat(gamma_rho) = h-hat(i/2) + h-hat(-i/2) - sum_{n>=2} (Lambda(n)/sqrt(n)) * h(log n) + Omega_h

Set h(x) = (f * f~)(x). Then:
- h-hat(tau) = |F(1/2 + i*tau)|^2 (by the convolution theorem)
- h(log n) = integral f(t) f-bar(t - log n) dt
- h-hat(i/2) = |F(0)|^2 and h-hat(-i/2) = |F(1)|^2

The spectral side becomes:

  sum_rho h-hat(gamma_rho) = sum_rho |F(rho)|^2

(using h-hat(gamma) = |F(1/2 + i*gamma)|^2 and rho = 1/2 + i*gamma).

The arithmetic side becomes:

  |F(0)|^2 + |F(1)|^2 - P(f, f) + Omega(f * f~)

Thus the identity is established.

### 5.3 The Prime-Zero Duality

The spectral identity reveals a fundamental duality. The prime bilinear form can be expressed:

  P(f, f) = |F(0)|^2 + |F(1)|^2 + Omega(f * f~) - sum_rho |F(rho)|^2

This says: **the prime contribution is the deficit between the pole/archimedean terms and the spectral sum over zeros.**

Dually, the spectral sum is:

  sum_rho |F(rho)|^2 = |F(0)|^2 + |F(1)|^2 + Omega(f * f~) - P(f, f)

**The spectral sum measures how much the pole and archimedean terms exceed the prime contribution.**

The cross-interactions between distinct primes p and q in the bilinear form P are therefore dual to correlations between zeros rho and rho'. Specifically:

- The "diagonal" prime contributions P_p(f, f) (single prime, all powers) correspond to the "smooth" part of the spectral density.
- The "cross" contributions (how P_p and P_q interact through f) correspond to the pair correlation of zeros.

### 5.4 Montgomery's Pair Correlation

Montgomery's pair correlation conjecture states that for the normalized zero spacings:

  R_2(alpha) = 1 - (sin(pi*alpha)/(pi*alpha))^2 + delta(alpha)

(the GUE pair correlation, matching random matrix theory).

**Montgomery's Theorem (1973):** Assuming RH, for test functions whose Fourier transform is supported in [-1, 1]:

  (1/N(T)) * sum_{0 < gamma, gamma' <= T} f((gamma - gamma') * log T / (2*pi)) -> integral f(alpha) * R_2(alpha) d*alpha

This restricted result is proven using the explicit formula — specifically, by computing:

  sum_{gamma,gamma'} |F(rho) F(rho')|^2 = (explicit formula)^2

and comparing the "diagonal" rho = rho' terms with the "off-diagonal" rho != rho' terms.

**Connection to the prime kernel:** The pair correlation R_2 determines the two-point statistics of the spectral measure {gamma_rho}. Via the explicit formula, these statistics are dual to the two-point statistics of {log p^m} weighted by Lambda/sqrt(n). The restriction to Fourier support in [-1,1] corresponds precisely to the fact that the prime sum directly controls correlations only up to a certain scale.

### 5.5 The Fourier-Dual Picture

Taking the (inverse) Fourier transform of the spectral identity:

  sum_rho e^{i*gamma_rho * x} = delta(x) + e^{x/2} + Omega-kernel(x) - K_prime(x)

(in the distributional sense, with appropriate regularization).

The prime kernel K_prime(x) captures the "oscillatory corrections" to the smooth density of zeros:

  K_prime(x) = [smooth terms from poles and archimedean] - [oscillatory spectral sum]

At a point x = log(p^a / q^b) (the logarithm of a ratio of prime powers), the kernel K_prime has contributions from BOTH the delta functions at log(p^a) and log(q^b) AND from the spectral oscillations sum_rho e^{i*gamma*x}.

**The cross-interaction between primes p and q is encoded in the value of the spectral sum at the "difference point" x = a*log p - b*log q:**

  sum_rho e^{i*gamma_rho * (a*log p - b*log q)} = e^{i*gamma_rho * log(p^a/q^b)}

This oscillatory sum over zeros determines how the prime-power contributions at p^a and q^b interfere.

### 5.6 Baker's Theorem and Arithmetic Independence

By the fundamental theorem of arithmetic, p^a != q^b for distinct primes p, q (unless a = b = 0). More quantitatively, Baker's theorem on linear forms in logarithms gives:

  |a * log p - b * log q| >= exp(-C * log a * log b * log p * log q)

for an effective constant C. This lower bound ensures that the "difference points" a*log p - b*log q are always nonzero, and bounds how close they can get.

This arithmetic independence means:
- The delta functions in K_prime at log(p^a) and log(q^b) never coincide.
- The test function f is evaluated at genuinely distinct points when computing cross-interactions.
- The cross-correlation integral f(t) f-bar(t - (a*log p - b*log q)) dt involves f at two points separated by at least exp(-C * ...) > 0.

---

## 6. Test Functions That Come Closest to Violating Positivity

### 6.1 Formulation of the Extremal Problem

We seek to understand which f in C_c^infty(R) minimize the ratio:

  W(f * f~) / ||f||^2

where ||f||^2 = integral |f(t)|^2 dt.

By the spectral identity, W(f * f~) = sum_rho |F(rho)|^2. Under RH, this is always >= 0, and the infimum of the ratio is 0.

**Question:** How quickly can W(f * f~) / ||f||^2 approach 0? What functions achieve near-zero values?

### 6.2 Type 1: Concentrating at a Single Zero

Let rho_0 = 1/2 + i*gamma_0 be a zero of zeta on the critical line. Consider:

  f_T(t) = e^{(1/2 + i*gamma_0)*t} * phi(t/T)

where phi in C_c^infty(R) is a smooth bump with phi(0) = 1 and supp(phi) subset [-1, 1]. Then:

  F_T(s) = T * phi-hat(T*(s - 1/2 - i*gamma_0))

This concentrates F near s = rho_0. Since phi-hat is a Schwartz function:

  |F_T(rho)|^2 ~ T^2 * |phi-hat(T*(rho - rho_0))|^2

The dominant contributions to sum_rho |F_T(rho)|^2 come from zeros rho near rho_0 (within ~ 1/T).

  W(f_T * f_T~) / ||f_T||^2 ~ sum_rho |phi-hat(T*(rho - rho_0))|^2 / integral |phi(x)|^2 dx

For large T, only zeros within distance O(1/T) of rho_0 contribute significantly. If gamma_0 is a simple zero with nearest neighbor at distance d (the local gap), then for T >> 1/d:

  W / ||f||^2 ~ |phi-hat(0)|^2 / ||phi||^2 = const > 0

The ratio stabilizes at a positive constant. Concentrating at a single simple zero does NOT make W/||f||^2 approach 0.

### 6.3 Type 2: Exploiting Close Pairs of Zeros (Lehmer-Type Phenomena)

If two zeros gamma_1, gamma_2 have an exceptionally small gap delta = |gamma_1 - gamma_2|, one can choose f with F bridging both zeros. The minimum of W/||f||^2 for such f is approximately:

  W / ||f||^2 >= c * delta^2 (when delta << mean spacing ~ 2*pi/log T)

The Lehmer-type near-violations correspond to delta = o(1/log T), giving:

  W / ||f||^2 >= c / (log T)^2

which is small but positive.

**The closest approach to violating positivity** comes from the zero pair with the smallest normalized gap in [0, T]. Conjecturally (under GUE statistics), the minimum gap in [0, T] is ~ 1/(T * (log T)^2), giving:

  min W/||f||^2 ~ 1 / (T * (log T)^2)^2 * (log T)^2 -> 0 as T -> infinity

So the infimum IS 0, but it is never achieved for any fixed T.

### 6.4 Type 3: The Arithmetic Side of Near-Extremality

On the arithmetic side (prime sum), near-extremal functions are those for which P(f, f) nearly matches |F(0)|^2 + |F(1)|^2 + Omega(f * f~). This requires:

**Condition:** f must produce large, coherently positive auto-correlations integral f(t) f-bar(t - log n) dt for many values of n simultaneously.

This means f should be "quasi-periodic" with approximate periods in the set {log n : Lambda(n) > 0}. In other words, f must resonate simultaneously with many prime-power scales.

The function f(t) = sum_{j} a_j * e^{i*gamma_j * t} (a superposition of zero-frequencies) achieves this resonance because the explicit formula forces the prime sums and zero sums to balance. This is the duality of Section 5.3 in action: functions that resonate with primes automatically resonate with zeros, and vice versa.

### 6.5 Li's Criterion and the Natural Test Sequence

Li's criterion provides a canonical sequence of test quantities:

  lambda_n = sum_rho [1 - (1 - 1/rho)^n]

RH holds if and only if lambda_n >= 0 for all n >= 1.

The associated test functions probe the Weil functional at increasing "frequencies." Under RH:

  lambda_n ~ (n/2) * log n + c_1 * n + c_2 + O(1/n)

The positivity margin grows: lambda_n / n -> infinity as n -> infinity. This shows that for the Li test functions, positivity becomes EASIER to satisfy at high frequencies, not harder.

**The hard regime** is small n (n = 1, 2, 3, ...), where the numerical values are:

  lambda_1 ~ 0.0231, lambda_2 ~ 0.0923, lambda_3 ~ 0.2077, ...

All positive, but lambda_1 is closest to 0. This reflects the fact that the hardest test functions are those of "moderate" frequency — not too concentrated (missing the pole contributions) and not too spread out (averaging away the prime interactions).

### 6.6 The Extremal Principle: Why RH is Self-Consistent

The most dangerous test functions must simultaneously:

1. **Minimize spectral weight at the poles:** |F(0)|^2 + |F(1)|^2 should be small relative to ||f||^2, so the pole terms cannot rescue positivity.

2. **Maximize coherent prime correlations:** P(f, f) should be as large as possible, requiring f to correlate with many prime scales.

3. **Concentrate spectral weight near closely-spaced zeros:** To minimize sum_rho |F(rho)|^2 / ||f||^2.

The fundamental tension: requirement (2) forces f to have specific oscillation patterns dictated by the prime distribution, while requirement (3) forces f to match the zero distribution. By the explicit formula, these are dual constraints — optimizing for one automatically constrains the other.

**This is the arithmetic essence of why RH is self-consistent.** No test function can simultaneously exploit the prime distribution to make P(f, f) large AND exploit the zero distribution to make sum_rho |F(rho)|^2 small, because the two distributions are linked by the explicit formula. The explicit formula acts as a "conservation law" preventing the Weil functional from becoming negative.

---

## 7. Summary: The Structure of the Cross-Terms

### 7.1 The Complete Picture

The Weil positivity criterion decomposes into four terms:

  W(f * f~) = POLE(f) - PRIME(f) + ARCH(f) = SPECTRAL(f)

| Component | Expression | Sign | Nature |
|-----------|-----------|------|--------|
| POLE | \|F(0)\|^2 + \|F(1)\|^2 | >= 0 | Rank 2, positive semi-definite |
| PRIME | P(f,f) = integral integral f f-bar K_prime | Indefinite | Integral operator, kernel K_prime |
| ARCH | Omega(f * f~) | Variable | Digamma/Gamma contributions |
| SPECTRAL | sum_rho \|F(rho)\|^2 | >= 0 iff RH | Sum over non-trivial zeros |

### 7.2 The Prime Kernel

  K_prime(x) = sum_{n >= 2} (Lambda(n)/sqrt(n)) * [delta(x - log n) + delta(x + log n)]

  K-hat_prime(tau) = -2 Re[zeta'/zeta(1/2 + i*tau)]

Properties:
- Discrete support on {+/- log n : Lambda(n) > 0}
- Exponentially decaying coefficients in prime-power index m
- Spectrally indefinite: K-hat_prime(0) < 0, K-hat_prime(tau) > 0 for large tau
- Spectral transform is the log-derivative of zeta on the critical line
- Decomposes additively over primes: K_prime = sum_p K_p

### 7.3 The Nature of the Cross-Terms

The cross-terms between primes p and q do not appear as explicit terms in the decomposition P = sum_p P_p. Rather, they manifest as **implicit interactions mediated by the test function f**.

The interaction between p and q through f is controlled by:

  Auto-correlations: integral f(t) f-bar(t - x) dt  at  x = a*log p - b*log q

These "difference points" are always nonzero (by unique factorization) and bounded below by Baker's theorem. The arithmetic independence of log p and log q for distinct primes ensures no single test function can simultaneously achieve perfect resonance at all prime scales.

### 7.4 The Connection to the Matrix Formulation

Discretizing the bilinear form on a basis {e_j}:

  M_{ij} = P(e_i, e_j) = integral integral e_i(t) e_j-bar(s) K_prime(t-s) dt ds

APT requires: (pole matrix) - M + (archimedean matrix) is positive semi-definite.

The diagonal entries of M are controlled by single-prime contributions (sum over m of log p / p^{m/2} * auto-correlation at shift m*log p). The off-diagonal entries involve cross-correlations at "mixed" shifts, controlled by the arithmetic of log p / log q ratios.

Diagonal dominance of this matrix — |M_{ii}| > sum_{j != i} |M_{ij}| — is sufficient for the required semi-definiteness and holds for all large primes (as shown in structure.md). The verification for small primes reduces to a finite computation, which is the operational frontier of this approach.

### 7.5 The Irreducible Core

The cross-terms encode the same information as the pair correlations between zeros of zeta. Their control is equivalent to RH itself. A proof of APT must show that these cross-terms cannot overwhelm the diagonal contributions — equivalently, that the primes are "sufficiently independent" in the sense that no test function can coherently amplify their interactions beyond the budget provided by the pole and archimedean terms.

This is the arithmetic content of the Riemann Hypothesis, reformulated as a positivity statement about an integral operator with an arithmetically-defined kernel.
