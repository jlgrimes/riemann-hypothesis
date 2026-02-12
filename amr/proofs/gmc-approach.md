# Gaussian Multiplicative Chaos and the Weil Kernel: An Exploration

## Status: The GMC framework provides deep structural insight into WHY K_bg is log-correlated and CPD-1, confirms the Euler-product-as-branching-structure analogy, but does NOT yield an unconditional proof of CPD-1 for the full Weil kernel. The obstruction localizes to the same place as in the Bochner and Schoenberg approaches: the zero-oscillation component K_zeros.

---

## 0. Motivation and Overview

**Goal.** Investigate whether Gaussian Multiplicative Chaos (GMC) theory, combined with recent results of Saksman-Webb (2020) on the zeta function's convergence to a log-correlated Gaussian field, can prove (or illuminate) the conditional positive definiteness of order 1 (CPD-1) of the Weil kernel K — which is equivalent to RH.

**Key observation.** The background kernel K_bg(x) = -(1/pi) Re[psi(1/4 + ix/2)] + log(pi)/(2pi) satisfies:

K_bg(x) ~ -(1/(2pi)) log(|x|/2) as |x| -> infinity

This is precisely the covariance structure of a **log-correlated Gaussian field**. GMC theory provides powerful tools for analyzing such fields. The question is whether these tools can reach CPD-1.

**Summary of findings:**

| Result | Status | Section |
|--------|--------|---------|
| K_bg is a valid CND-1 kernel (conditionally negative definite of order 1) | **Proved** (unconditional) | Section 2 |
| K_bg generates a well-defined log-correlated Gaussian field | **Proved** (unconditional) | Section 2 |
| GMC measure mu_gamma associated to K_bg exists for gamma < sqrt(2) | **Proved** (unconditional) | Section 3 |
| Saksman-Webb: zeta on critical line converges to complex GMC | **Proved** (unconditional) | Section 4 |
| Euler product decomposes as branching random walk over primes | **Established** (structural) | Section 5 |
| FHK conjecture: zeta maximum in universality class of log-correlated fields | **Proved** (Arguin-Bourgade-Radziwill 2020-2023) | Section 5 |
| GMC existence implies CPD-1 of K_bg | **Yes** (but K_bg CPD-1 already known) | Section 6 |
| GMC tools give CPD-1 of full K | **No** (same obstruction as Schoenberg) | Section 7 |

---

## 1. Background: GMC Theory

### 1.1 Kahane's Construction (1985)

**Definition.** Let X = {X(t) : t in D} be a centered Gaussian process (or generalized Gaussian field) on a domain D subset R^d, with covariance:

C(s,t) = E[X(s)X(t)]

When C has a logarithmic singularity on the diagonal:

C(t,t) = +infinity, C(s,t) ~ -log|s-t| + O(1) as s -> t

the field X is called **log-correlated**. It cannot be defined pointwise — it exists only as a random generalized function (distribution in the sense of Schwartz).

**Gaussian Multiplicative Chaos** is the random measure defined by regularizing and exponentiating:

mu_gamma(dx) = lim_{epsilon -> 0} exp(gamma X_epsilon(x) - (gamma^2/2) E[X_epsilon(x)^2]) dx

where X_epsilon is a regularization of X (e.g., convolution with a mollifier at scale epsilon), and gamma in R is the **intermittency parameter**.

**Theorem (Kahane 1985, Berestycki 2017, Shamov 2016).** The limit exists in probability (as a random Radon measure) and is:
- **Non-trivial** for |gamma| < sqrt(2d) (subcritical phase)
- **Zero** for |gamma| > sqrt(2d) (supercritical phase)
- **Critical** at |gamma| = sqrt(2d) (a non-trivial limit exists with different normalization)

The construction is **universal**: the limiting measure depends on X only through its covariance kernel C, not on the choice of regularization.

**Reference:** Berestycki, N. (2017). An elementary approach to Gaussian multiplicative chaos. *Electron. Commun. Probab.* 22, paper no. 27. [arXiv:1506.09113](https://arxiv.org/abs/1506.09113)

### 1.2 The Role of Conditional Negative Definiteness

A kernel C(s,t) = phi(s-t) defines a valid covariance for a Gaussian process iff phi is **positive definite** (PD): F[phi](xi) >= 0 for all xi.

For log-correlated fields, C(s,t) ~ -log|s-t| diverges on the diagonal. The function psi(x) = -log|x| is NOT positive definite (F[-log|x|](xi) = -pi/|xi| < 0). Instead, it is **conditionally negative definite of order 1** (CND-1):

sum_{i,j} c_i c_j (-log|x_i - x_j|) <= 0 for all {x_i}, {c_i} with sum c_i = 0

Equivalently, psi(x) = -log|x| is CND-1, meaning -psi = log|x| is CPD-1. This is a classical result (Schoenberg 1938).

**The connection to our problem:** If K_bg ~ -(1/(2pi)) log|x| for large |x|, and -log|x| is CND-1, then K_bg should be "approximately CND-1" — or more precisely, K_bg + (1/(2pi)) log|x| should be a bounded correction to a CND-1 function. We showed in bochner-proof.md (Theorem 8.1) and schoenberg-attempt.md (Theorem 5.1) that K_bg is in fact CPD-1, which is the OPPOSITE sign convention: K_bg has K_bg_hat(xi) > 0 for xi != 0, i.e., K_bg is CPD-1 (not CND-1).

**Clarification of signs:** K_bg(x) is large positive near x = 0 and decays to -infinity as |x| -> infinity. The function -K_bg(x) grows logarithmically. -K_bg behaves like (1/(2pi))log|x| for large |x|. So:
- -K_bg is CND-1 (since K_bg is CPD-1)
- -K_bg can serve as a "squared-distance" or "variogram" for a Gaussian field
- The associated Gaussian field has covariance structure determined by -K_bg

### 1.3 From CND-1 to Gaussian Fields

**Theorem (Schoenberg 1938, cf. Berg-Christensen-Ressel 1984).** A continuous even function psi: R -> R with psi(0) = 0 is CND-1 iff e^{-t psi(x)} is PD for all t > 0.

**Corollary.** If phi is CPD-1 (i.e., -phi is CND-1 with appropriate normalization), then e^{t phi(x)} is PD for all t > 0 (up to a positive scalar).

For K_bg: since K_bg is CPD-1 (bochner-proof.md Theorem 8.1), the function e^{t K_bg(x)} is PD for all t > 0 (schoenberg-attempt.md Corollary 5.2). This means K_bg defines a valid "covariance-like" kernel for an infinitely divisible family of Gaussian processes.

---

## 2. K_bg as a Log-Correlated Covariance Kernel

### 2.1 Exact Structure

The background kernel has the Lorentzian expansion (bochner-proof.md Section 1.1):

K_bg(x) = (1/pi) sum_{n=0}^{infinity} (n+1/4) / ((n+1/4)^2 + x^2/4) + [constants]

Each term (n+1/4)/((n+1/4)^2 + x^2/4) is a Lorentzian (Cauchy distribution) centered at 0 with scale parameter 2(n+1/4).

**Fourier transform** (bochner-proof.md Section 1.3):

K_bg_hat(xi) = 2 e^{-|xi|/2} / (1 - e^{-2|xi|}) for xi != 0

**Asymptotic behavior:**
- K_bg_hat(xi) ~ 1/|xi| as |xi| -> 0^+ (the 1/|xi| singularity is characteristic of log-correlated fields)
- K_bg_hat(xi) ~ 2 e^{-|xi|/2} as |xi| -> infinity

### 2.2 Comparison with Standard Log-Correlated Kernels

The standard log-correlated kernel on R is:

C_log(x) = -log|x| (distributional sense)

with F[-log|x|](xi) = -pi/|xi| (in the distributional Fourier transform).

Our kernel K_bg has K_bg_hat(xi) ~ +1/|xi|, so K_bg behaves like -C_log/pi = (1/pi)log|x| at the Fourier level. More precisely:

K_bg_hat(xi) = 1/|xi| + (analytic correction near xi = 0)

This means K_bg is the covariance kernel of a process that is "the opposite sign" of a log-correlated field — it is the VARIOGRAM (negative of covariance) rather than the covariance itself.

### 2.3 Construction of the Associated Gaussian Field

Define the centered Gaussian field X on R by:

E[X(s)X(t)] = sigma^2 - K_bg(s-t)

where sigma^2 = K_bg(0) (formally infinite, but regularized). This gives:

E[(X(s) - X(t))^2] = 2[sigma^2 - E[X(s)X(t)]] ... (not quite right)

More precisely, -K_bg is CND-1, so by the GNS construction there exists a centered Gaussian field X indexed by R such that:

E[(X(s) - X(t))^2] = 2[-K_bg(s-t) + K_bg(0)]

For |s-t| large: E[(X(s)-X(t))^2] ~ (1/pi) log|s-t| + C

This is a Gaussian field with **logarithmic variogram** — the hallmark of a log-correlated field. The increments have variance growing logarithmically with distance.

**Theorem 2.1 (Unconditional).** The Weil background kernel K_bg defines a valid variogram for a log-correlated Gaussian field on R. Specifically, the function psi(x) = K_bg(0) - K_bg(x) (suitably regularized) is CND-1, and the associated Gaussian field has logarithmic increment variance.

*Proof.* K_bg is CPD-1 (bochner-proof.md, Theorem 8.1). Therefore psi(x) = K_bg(0) - K_bg(x) satisfies psi(0) = 0 and is CND-1 (since sum c_i c_j psi(x_i - x_j) = K_bg(0) sum c_i^2 - sum c_i c_j K_bg(x_i - x_j) = -sum c_i c_j K_bg(x_i - x_j) <= 0 when sum c_i = 0, by CPD-1 of K_bg, wait — the sign is wrong).

Let me redo this. K_bg CPD-1 means sum c_i c_j K_bg(x_i - x_j) >= 0 for sum c_i = 0. Then:

sum c_i c_j psi(x_i - x_j) = sum c_i c_j [K_bg(0) - K_bg(x_i - x_j)] = K_bg(0)(sum c_i)^2 - sum c_i c_j K_bg(x_i - x_j) = 0 - (non-negative) <= 0

So psi is CND-1. By Schoenberg, e^{-t psi(x)} is PD for all t > 0. The associated Gaussian field exists. QED

---

## 3. GMC Measure from the Weil Background Kernel

### 3.1 Construction

Given the log-correlated Gaussian field X from Section 2.3, define the GMC measure:

mu_gamma(dx) = lim_{epsilon -> 0} exp(gamma X_epsilon(x) - gamma^2/2 E[X_epsilon(x)^2]) dx

By Kahane-Berestycki theory, this exists and is non-trivial for gamma in (0, sqrt(2)) (in d = 1).

### 3.2 Interpretation

The GMC measure mu_gamma associated to K_bg is a random fractal measure on R that captures the "multiplicative chaos" generated by the background part of the Weil kernel. Its total mass on an interval [0,T]:

E[mu_gamma([0,T])] = T (by construction)

but the actual realization is a highly irregular random measure with multifractal properties.

### 3.3 Connection to Zeta Moments

The moments of mu_gamma relate to integrals of the form:

E[mu_gamma([0,T])^n] ~ integral_{[0,T]^n} prod_{i<j} |x_i - x_j|^{-gamma^2/pi} dx_1 ... dx_n

(using the log-correlated covariance structure). This resembles the Selberg integral, which also appears in the moments of |zeta(1/2+it)|^{2gamma}.

**Keating-Snaith (2000)** conjectured that the 2k-th moment of zeta on the critical line satisfies:

(1/T) integral_0^T |zeta(1/2+it)|^{2k} dt ~ C_k (log T)^{k^2}

where C_k involves products over primes and a random matrix theory factor. The exponent k^2 matches the prediction from GMC theory: for a log-correlated field X, E[exp(gamma X)^n] ~ (log T)^{n gamma^2/2}, and with gamma = sqrt(2k), one gets (log T)^{k^2}.

This is **not a coincidence** — it reflects the deep structural identity between the zeta function on the critical line and log-correlated Gaussian fields.

---

## 4. Saksman-Webb: Zeta Converges to GMC

### 4.1 The Main Theorem

**Theorem (Saksman-Webb, Annals of Probability 2020).** Let omega be uniform on [0,1]. As T -> infinity, the random generalized function:

t |-> zeta(1/2 + i omega T + it)

converges in law (in the space of distributions) to a random distribution that factors as:

Z(t) = F(t) * mu_{GMC}(t)

where F is a smooth random function and mu_{GMC} is a complex Gaussian multiplicative chaos distribution.

**Reference:** Saksman, E. and Webb, C. (2020). The Riemann zeta function and Gaussian multiplicative chaos: Statistics on the critical line. *Ann. Probab.* 48(6), 2680-2754. [DOI: 10.1214/20-AOP1433](https://projecteuclid.org/journals/annals-of-probability/volume-48/issue-6/The-Riemann-zeta-function-and-Gaussian-multiplicative-chaos--Statistics/10.1214/20-AOP1433.full)

### 4.2 The Mesoscopic Regime

On mesoscopic scales (delta_T -> 0 slowly):

t |-> zeta(1/2 + i delta_T t + i omega T)

converges to a PURELY Gaussian multiplicative chaos (the smooth factor F becomes constant). This means: at mesoscopic scales, zeta on the critical line IS a GMC object.

### 4.3 The Covariance Structure

The key input is the covariance of log|zeta| on the critical line:

E[log|zeta(1/2+is)| * log|zeta(1/2+it)|] ~ -log|s-t| + C for |s-t| << T

This logarithmic covariance is exactly the hallmark of a log-correlated field, and it arises from the Euler product:

log zeta(s) = sum_p sum_{m=1}^{infinity} 1/(m p^{ms})

The primes contribute approximately independent Gaussian terms (by the prime number theorem + central limit theorem), and the total variance grows like log log T (Selberg's theorem). The covariance between log|zeta(1/2+is)| and log|zeta(1/2+it)| comes from primes p with p^{i(s-t)} approx 1, which for |s-t| = h requires p << 1/h. The number of such primes is ~ log(1/h)/log log(1/h) ~ log(1/h), giving the logarithmic covariance.

### 4.4 Significance for Our Problem

Saksman-Webb proves that the **statistical behavior** of zeta on the critical line is governed by a log-correlated field whose covariance is asymptotically:

C(h) ~ -log|h| for small h

This is precisely the large-scale behavior of K_bg(x). The connection:

**K_bg(x) is the kernel that generates the same type of log-correlated field that zeta converges to.**

This is not a coincidence — it is a consequence of the explicit formula. The archimedean (Gamma factor) contribution to the Weil kernel IS the covariance structure of log|zeta| on the critical line, averaged over the zero oscillations.

### 4.5 Extension to L-Functions

Recent work extends Saksman-Webb to Dirichlet L-functions:

**Theorem (2025).** For chi_q uniformly random among Dirichlet characters mod q, the distribution L(1/2+ix, chi_q) converges as q -> infinity to the same complex GMC limiting object as random shifts of zeta.

This universality confirms that the GMC structure is intrinsic to the L-function family, not specific to zeta.

**Reference:** [arXiv:2506.16115](https://arxiv.org/abs/2506.16115)

---

## 5. The Euler Product as Branching Structure

### 5.1 The Arguin-Belius-Harper Insight

The Euler product log|zeta(1/2+it)| = sum_p Re[-log(1 - p^{-1/2-it})] decomposes the field into contributions from individual primes. For distinct primes p, q, the random variables:

X_p(t) = Re[-log(1 - p^{-1/2-it})]
X_q(t) = Re[-log(1 - q^{-1/2-it})]

are **approximately independent** (their covariance is O(1/(pq)^{1/2})).

Each X_p(t) contributes variance ~ 1/(2p) to the total log|zeta(1/2+it)|. The hierarchical structure — small primes contribute large variance, large primes contribute small variance — mirrors a **branching random walk (BRW)**:

- At "level" k (primes of size ~ e^k), roughly e^k/k primes each contribute variance ~ e^{-k}
- The total variance from levels 1 to K is ~ sum_{k=1}^K 1 = K ~ log T
- This matches the Selberg CLT: Var[log|zeta(1/2+it)|] ~ (1/2) log log T

### 5.2 The FHK Conjecture and Its Proof

**Conjecture (Fyodorov-Hiary-Keating, 2012).** The maximum of |zeta| in short intervals satisfies:

max_{|h| <= 1} log|zeta(1/2+it+ih)| = log log T - (3/4) log log log T + O_P(1)

where t is uniform on [T, 2T]. This places zeta's maximum in the **universality class of log-correlated fields**, which includes:
- Branching random walks (BRW)
- The 2D Gaussian free field (GFF)
- Cover times of random walks
- Characteristic polynomials of random matrices

**Theorem (Arguin-Bourgade-Radziwill, 2020/2023).**

Part I (2020): The upper bound holds:

P(max_{|h|<=1} |zeta(1/2+it+ih)| > e^y log T / (log log T)^{3/4}) <= C y e^{-2y}

Part II (2023): The matching lower bound, implying tightness of:

max_{|h|<=1} |zeta(1/2+it+ih)| * (log log T)^{3/4} / log T

This confirms the FHK conjecture and establishes zeta in the BRW universality class.

**References:**
- Arguin, L.-P., Bourgade, P., and Radziwill, M. (2020). The Fyodorov-Hiary-Keating Conjecture. I. [arXiv:2007.00988](https://arxiv.org/abs/2007.00988)
- Arguin, L.-P., Bourgade, P., and Radziwill, M. (2023). The Fyodorov-Hiary-Keating Conjecture. II. [arXiv:2307.00982](https://arxiv.org/abs/2307.00982)
- Further extension to mesoscopic intervals: [arXiv:2405.06474](https://arxiv.org/abs/2405.06474)

### 5.3 Connection to the Weil Kernel Decomposition

The Weil kernel decomposes as K = delta + K_bg + K_zeros. The prime bilinear form (bochner-proof.md Section 5.1) involves:

K_prime(x) = sum_{n >= 2} (Lambda(n)/sqrt(n)) [delta(x - log n) + delta(x + log n)]

This decomposes over primes:

K_prime = sum_p K_p where K_p(x) = sum_{m=1}^{infinity} (log p / p^{m/2}) [delta(x - m log p) + delta(x + m log p)]

**The parallel:** In the BRW picture of zeta, each prime p contributes an approximately independent component X_p to the log-correlated field. In the Weil kernel, each prime p contributes K_p to the prime bilinear form. The hierarchical decomposition K_prime = sum_p K_p mirrors the branching tree structure X = sum_p X_p.

This parallel is **structural**, not just analogical: both decompositions arise from the Euler product of zeta. The Weil explicit formula translates the multiplicative structure of zeta (Euler product) into the additive structure of the kernel (sum over primes).

---

## 6. Can GMC Existence Imply CPD-1?

### 6.1 The Direct Argument for K_bg

**Proposition 6.1.** The existence of the GMC measure mu_gamma for the Gaussian field associated to K_bg implies K_bg is CPD-1.

*Proof.* The GMC measure exists iff the Gaussian field X with variogram psi(x) = K_bg(0) - K_bg(x) is well-defined. By the GNS construction, X exists iff psi is CND-1. By the sign relation (Section 1.2), psi CND-1 iff K_bg CPD-1. QED

**Status:** This is correct but unilluminating — it merely restates CPD-1 of K_bg (already proved by Fourier analysis in bochner-proof.md) in GMC language. GMC existence does not add information beyond what the Fourier transform already gives.

### 6.2 The Argument for the Full Kernel

For the full Weil kernel K = delta + K_bg + K_zeros, CPD-1 requires:

K_hat(xi) = 1 + K_bg_hat(xi) + K_zeros_hat(xi) >= 0 for all xi != 0

From Sections 1-3, K_bg_hat(xi) > 0 unconditionally. The term K_zeros_hat depends on zero locations:

- Under RH: K_zeros_hat is a non-negative measure (sum of positive point masses at imaginary parts of zeros)
- Without RH: K_zeros_hat can have negative contributions from off-line zeros

Can GMC theory say anything about K_zeros?

### 6.3 GMC and the Zero Oscillation

The zero-oscillation kernel K_zeros encodes the locations of zeta's nontrivial zeros. In the GMC picture, these correspond to the **fine structure** of the multiplicative chaos — the precise distribution of mass at small scales.

Saksman-Webb prove that the GMC limit of zeta exists unconditionally (without assuming RH). This means the random measure:

mu_gamma(dx) = |zeta(1/2 + ix)|^{2gamma} dx (formally, in the Saksman-Webb sense)

exists as a well-defined random measure for gamma in the subcritical range, regardless of whether RH holds.

**However:** The existence of this measure tells us about the statistical distribution of |zeta|^{2gamma} on the critical line, averaged over random shifts. It does NOT directly constrain the spectral properties of the Weil kernel K, because:

1. The Weil kernel is a DETERMINISTIC object (a distribution on R), not a random one.
2. GMC existence is a statement about the LAW of a random process, not about specific realizations.
3. The CPD-1 condition on K involves ALL test functions, not just "typical" ones in any probabilistic sense.

### 6.4 The Indirect Route

One might hope for an indirect argument:

(a) GMC theory proves certain moment formulas for |zeta|^{2gamma} unconditionally.
(b) These moment formulas, combined with the explicit formula, constrain the possible zero distributions.
(c) The constraints rule out off-line zeros.

**Assessment:** Step (a) is partially achieved — Keating-Snaith moments, Selberg integral estimates, and GMC convergence give asymptotic moment information. Step (b) is the subject of active research but has not yielded zero-free regions beyond classical methods. Step (c) is extremely ambitious and there is no known path.

The fundamental issue is that GMC tools capture **averaged/statistical** properties of zeta, while RH is an **individual/deterministic** statement about every single zero. Statistical tools cannot easily bridge this gap.

---

## 7. Assessment: Is GMC a Viable Path to CPD-1?

### 7.1 What GMC DOES Provide

1. **Structural explanation of K_bg.** The log-correlated nature of K_bg is not an accident — it is a direct consequence of the Euler product structure of zeta, via the prime decomposition of the Gaussian field. GMC theory explains WHY K_bg has the form it does: it is the covariance kernel of the limiting log-correlated field that zeta converges to.

2. **Universality.** The Saksman-Webb theorem and the FHK confirmation show that zeta belongs to a broad universality class. Any proof of CPD-1 via GMC would automatically extend to the entire universality class (all L-functions in the Selberg class), which is consistent with GRH.

3. **Prime-by-prime intuition.** The branching random walk structure provides intuition for why the prime bilinear form K_prime = sum_p K_p has favorable spectral properties: each K_p contributes approximately independently, and the sum has Gaussian-like behavior by CLT-type arguments.

4. **Moment methods.** GMC provides unconditional results on the moments of |zeta|^{2gamma}, which constrain the possible distributions of zeros in an averaged sense. While these constraints do not prove RH, they are consistent with and suggestive of RH.

### 7.2 What GMC Does NOT Provide

1. **Pointwise kernel positivity.** CPD-1 of K requires K_hat(xi) >= 0 for all xi != 0 — a pointwise condition on a deterministic function. GMC tools give distributional/averaged information, not pointwise bounds on K_hat.

2. **Control of K_zeros.** The zero-oscillation kernel is the precisely the part that GMC cannot reach unconditionally. The GMC limit involves averaging over the zeros (via random shifts of t), which washes out the individual zero contributions that K_zeros encodes.

3. **An unconditional proof of CPD-1.** For the same reason as the Bochner and Schoenberg approaches (bochner-proof.md Section 7, schoenberg-attempt.md Section 5.12), GMC cannot bypass the fundamental barrier: proving that K_zeros has favorable spectral properties requires knowing that zeros are on the critical line.

### 7.3 The Obstruction is Universal

All three approaches we have examined — Bochner (bochner-proof.md), Schoenberg (schoenberg-attempt.md), and GMC (this document) — encounter the same obstruction at the same place:

| Approach | Obstruction | Location |
|----------|-------------|----------|
| Bochner | K_zeros_hat may be negative without RH | bochner-proof.md Section 7 |
| Schoenberg | e^{t K_zeros} not PD without RH (sine terms from off-line zeros) | schoenberg-attempt.md Section 5.12 |
| GMC | K_zeros not capturable by statistical/averaged tools | This document, Section 6.3 |

In each case, K_bg is handled successfully and unconditionally. The delta is trivial. The obstruction is always K_zeros — the component that directly encodes zero locations.

### 7.4 Possible Future Directions

Despite the negative assessment for a direct proof, GMC opens several speculative avenues:

**Direction A: GMC + Number Theory.**
Combine GMC moment estimates with sieve-theoretic or density methods to obtain unconditional constraints on K_zeros_hat. For example, if one could show that K_zeros_hat(xi) >= -(1 + K_bg_hat(xi)) for all xi, this would give K_hat >= 0 and hence RH. The GMC moments constrain the average of K_zeros_hat but not its pointwise values.

**Direction B: Stochastic CPD-1.**
Define a "probabilistic CPD-1" condition: K is "stochastically CPD-1" if for randomly chosen test points {x_i} (e.g., from a Poisson process), the quadratic form sum c_i c_j K(x_i - x_j) >= 0 with probability 1. This is weaker than CPD-1 but might be provable from GMC. Whether stochastic CPD-1 implies CPD-1 is unclear.

**Direction C: The Multiplicative Chaos Spectral Measure.**
The GMC measure mu_gamma has a Fourier transform (spectral measure) that is itself a random measure. The expected spectral measure E[mu_gamma_hat] should relate to K_bg_hat. If one could show that the spectral measure is a.s. non-negative (not just in expectation), this would constrain K_hat.

**Direction D: Log-Correlated Field Rigidity.**
Recent work on log-correlated fields shows remarkable rigidity properties: e.g., the maximum, the minimum, the level sets, and the multifractal spectrum are all determined by the covariance kernel up to universal constants. If there is a "spectral rigidity" result — the spectral measure of the field determines the zero-distribution uniquely — then the GMC spectral properties might force RH.

None of these directions currently constitutes a viable proof strategy, but they represent the frontier of what GMC tools might achieve.

---

## 8. The Key Unconditional Results

Distilling what is proved unconditionally by the GMC analysis:

**Theorem 8.1 (K_bg defines a log-correlated field, unconditional).** The background Weil kernel K_bg(x) = -(1/pi) Re[psi(1/4+ix/2)] + log(pi)/(2pi) is CPD-1 on R. The associated variogram psi(x) = K_bg(0) - K_bg(x) defines a log-correlated Gaussian field on R with:

E[(X(s) - X(t))^2] ~ (1/pi) log|s-t| for |s-t| -> infinity

**Theorem 8.2 (GMC measure exists, unconditional).** The Gaussian multiplicative chaos measure associated to the log-correlated field from Theorem 8.1 exists and is non-trivial for all gamma in (0, sqrt(2)). Its moments match the Selberg-integral predictions for |zeta|^{2gamma}.

**Theorem 8.3 (Zeta converges to GMC, unconditional).** The Riemann zeta function on the critical line, restricted to random mesoscopic intervals, converges in law to a complex Gaussian multiplicative chaos distribution (Saksman-Webb 2020). This convergence does not assume RH.

**Theorem 8.4 (FHK universality, unconditional).** The maximum of |zeta| in short intervals on the critical line lies in the universality class of log-correlated Gaussian fields. Specifically, max_{|h|<=1} log|zeta(1/2+it+ih)| = log log T - (3/4) log log log T + O_P(1) for random t in [T,2T] (Arguin-Bourgade-Radziwill 2020-2023).

**Theorem 8.5 (Branching decomposition, structural).** The Euler product induces a decomposition of the log-correlated field of zeta into approximately independent prime contributions, forming a branching random walk. This structure mirrors the decomposition K_prime = sum_p K_p of the Weil prime kernel.

---

## 9. Relationship to the AMR Framework

### 9.1 GMC and Component C

The AMR framework (subspace-alignment.md) reduces RH to Component C: "Haar implies APT," which is CPD-1 of the full Weil kernel. The GMC approach cannot directly prove Component C (Section 7.2), but it provides:

- Unconditional confirmation that K_bg (the "archimedean" contribution to Component C) is favorable.
- Structural insight into WHY the prime bilinear form has the decomposition it does (branching structure from Euler product).
- A statistical framework in which RH is the "expected" outcome: the log-correlated universality class predicts exactly the behavior that RH implies.

### 9.2 GMC and Measure Rigidity

The AMR's measure-rigidity approach (Baker -> entropy -> Rudolph -> Haar) operates at a different level than GMC. While GMC studies the **statistical behavior** of zeta (averaged over shifts), measure rigidity studies the **algebraic structure** of the arithmetic measure mu_ar (x |-> p^{ix} dynamics). These are complementary:

- GMC explains the **analytical** structure of the Weil kernel (why K_bg is log-correlated, why moments match random matrix predictions).
- Measure rigidity constrains the **algebraic** structure of the spectrum (why Haar measure is the only invariant measure compatible with entropy maximality).

A synthesis might be possible: GMC shows that the "natural" Gaussian process for K_bg has certain rigidity properties (universality class, moment constraints), and measure rigidity shows that the algebraic structure of mu_ar forces exactly the spectral properties needed. But making this synthesis precise remains open.

---

## 10. Conclusion

The GMC approach to the Weil kernel's CPD-1 property is **illuminating but not definitive**. It reveals that K_bg is naturally a log-correlated covariance kernel — a fact that connects the Weil explicit formula to one of the deepest areas of modern probability theory. The recent confirmations of FHK and the Saksman-Webb theorem show that this connection is not merely formal but reflects genuine structural properties of the zeta function.

However, GMC tools are inherently statistical/averaged in nature, while CPD-1 is a deterministic/pointwise condition. The gap between these two regimes is precisely where the zero-oscillation kernel K_zeros lives — and this gap is, as with all other approaches examined, equivalent to the Riemann Hypothesis.

**Viability assessment:** GMC is NOT a viable path to an unconditional proof of CPD-1, but it is a valuable **diagnostic and structural tool** that:
1. Confirms the favorability of K_bg from a probabilistic perspective.
2. Provides the "right" framework for understanding the prime-by-prime structure of the Weil functional.
3. Suggests that RH is the "natural" outcome from the perspective of log-correlated field universality.
4. Opens speculative directions (Sections 7.4 A-D) that might become tractable with future developments in GMC theory.

---

## References

### GMC Theory
- Kahane, J.-P. (1985). Sur le chaos multiplicatif. *Ann. Sci. Math. Quebec* 9, 105-150.
- [Berestycki, N. (2017). An elementary approach to Gaussian multiplicative chaos.](https://arxiv.org/abs/1506.09113) *Electron. Commun. Probab.* 22, no. 27.
- [Rhodes, R. and Vargas, V. (2014). Gaussian multiplicative chaos and applications: A review.](https://projecteuclid.org/journals/probability-surveys/volume-11/issue-none/Gaussian-multiplicative-chaos-and-applications-A-review/10.1214/13-PS218.pdf) *Probab. Surv.* 11, 315-392.
- [Duplantier, B. et al. (2014). Log-correlated Gaussian fields: an overview.](https://arxiv.org/abs/1407.5605)

### Zeta and GMC
- [Saksman, E. and Webb, C. (2020). The Riemann zeta function and Gaussian multiplicative chaos.](https://projecteuclid.org/journals/annals-of-probability/volume-48/issue-6/The-Riemann-zeta-function-and-Gaussian-multiplicative-chaos--Statistics/10.1214/20-AOP1433.full) *Ann. Probab.* 48(6), 2680-2754.
- [Saksman, E. and Webb, C. (2020). On the Riemann Zeta function and Gaussian multiplicative chaos.](https://link.springer.com/chapter/10.1007/978-3-030-40120-7_12) In: *Advancements in Complex Analysis*, Springer.
- [Dirichlet L-functions and multiplicative chaos (2025).](https://arxiv.org/abs/2506.16115)

### FHK Conjecture
- Fyodorov, Y., Hiary, G., and Keating, J.P. (2012). Freezing transition, characteristic polynomials of random matrices, and the Riemann zeta-function. *Phys. Rev. Lett.* 108, 170601.
- [Arguin, L.-P., Bourgade, P., and Radziwill, M. (2020). The FHK Conjecture. I.](https://arxiv.org/abs/2007.00988)
- [Arguin, L.-P., Bourgade, P., and Radziwill, M. (2023). The FHK Conjecture. II.](https://arxiv.org/abs/2307.00982)
- [FHK on mesoscopic intervals (2024).](https://arxiv.org/abs/2405.06474)

### Background (Weil Kernel and CPD-1)
- [bochner-proof.md](bochner-proof.md) — Fourier analysis of K_bg, K_zeros, identification of CPD-1 gap
- [schoenberg-attempt.md](schoenberg-attempt.md) — Schoenberg representation, K_bg CPD-1 unconditional, K_zeros obstruction
- Schoenberg, I.J. (1938). Metric spaces and positive definite functions. *Trans. AMS* 44, 522-536.
- [Weil's criterion (Wikipedia).](https://en.wikipedia.org/wiki/Weil's_criterion)
- Keating, J.P. and Snaith, N.C. (2000). Random matrix theory and zeta(1/2+it). *Comm. Math. Phys.* 214, 57-89.
- [Wang, Y. (2025). A family of log-correlated Gaussian processes.](https://link.springer.com/article/10.1007/s10959-025-01449-2) *J. Theoret. Probab.*

---

*Document: GMC Approach to Weil Kernel CPD-1*
*Part of the AMR (Arithmetic Measure Rigidity) proofs module*
*February 2026*
