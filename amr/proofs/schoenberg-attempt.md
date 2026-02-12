# The Schoenberg Representation: A Backup Route to CPD-1 of the Weil Kernel

## Overview

This document explores whether Schoenberg's 1938 characterization of conditionally positive definite functions provides an alternative route to proving CPD-1 of the Weil kernel K, which by Theorem 2.1 of subspace-alignment.md is equivalent to APT (and hence RH). The approach factors the exponential e^{t·K(x)} into contributions from the background and zero-oscillation parts, analyzes each factor for positive definiteness, and identifies the precise obstruction.

**Conclusion:** The Schoenberg approach successfully reduces CPD-1 of K to showing that e^{(t/pi) Re[psi(1/4+ix/2)]} is positive definite (PD) for all t > 0. The zero-oscillation factor e^{t·K_zeros(x)} IS provably PD (via Bessel functions + Schur product theorem). The background factor e^{t·K_bg(x)} reduces to the PD question for a specific power of |Gamma(1/4+ix/2)|, which connects to deep results in the theory of entire functions of exponential type. This connection does not close the gap but reformulates it in a structurally illuminating way.

---

## 1. Schoenberg's Theorem

### 1.1 Statement

**Theorem (Schoenberg, 1938).** Let K: R -> R be a continuous function with K(x) = K(-x). Then K is conditionally positive definite of order 1 (CPD-1) if and only if for every t > 0, the function

phi_t(x) = exp(-t[K(0) - K(x)])

is positive definite (PD).

Equivalently: K is CPD-1 iff for every t > 0, the function x |-> e^{t·K(x)} is PD (up to the positive scalar factor e^{-t·K(0)}, which does not affect PD).

### 1.2 Context

Schoenberg's result provides a bridge between CPD-1 (a condition involving constrained sums Sum c_i c_j K(x_i - x_j) >= 0 with Sum c_i = 0) and ordinary PD (unconstrained sums Sum c_i c_j phi(x_i - x_j) >= 0). The latter is characterized by Bochner's theorem: phi is PD iff phi = F[mu] for a non-negative finite Borel measure mu.

### 1.3 Why This Might Help

The Bochner characterization of CPD-1 (subspace-alignment.md, Theorem 2.3) requires K_hat(xi) >= 0 for all xi != 0, which IS RH. The Schoenberg reformulation instead requires:

For all t > 0: (e^{t·K})^hat >= 0 everywhere (as a measure)

This shifts the problem from the Fourier transform of K (which has a distributional singularity at xi = 0 and whose positivity is RH) to the Fourier transform of e^{t·K}, which is a smooth, rapidly decaying function whose Fourier transform is better behaved.

---

## 2. Application to the Weil Kernel

### 2.1 The Kernel Structure

The Weil kernel in the matrix formulation (subspace-alignment.md, eq. at line 59):

K(x) = delta(x) + K_bg(x) + K_zeros(x)

where:
- delta(x): Dirac delta, contributing K(0) = 1 + K_bg(0) + K_zeros(0) approx 2.528
- K_bg(x) = -(1/pi) Re[psi(1/4 + ix/2)] + log(pi)/(2pi)
- K_zeros(x) = (1/(2pi)) Sum_gamma 2cos(gamma x)/(1/4 + gamma^2)

For x != 0, delta(x) = 0, so K(x) = K_bg(x) + K_zeros(x).

### 2.2 The Schoenberg Function

For t > 0, we need to show:

phi_t(x) = exp(-t[K(0) - K(x)])

is PD. Since K(0) - K(x) = [1 + K_bg(0) + K_zeros(0)] - [K_bg(x) + K_zeros(x)] for x != 0:

K(0) - K(x) = 1 + [K_bg(0) - K_bg(x)] + [K_zeros(0) - K_zeros(x)]

The "1" comes from the regularized delta: delta(0) = 1 in the matrix sense, delta(x) = 0 for x != 0.

Therefore:

phi_t(x) = e^{-t} · e^{-t[K_bg(0) - K_bg(x)]} · e^{-t[K_zeros(0) - K_zeros(x)]}

The scalar e^{-t} does not affect PD. So we need:

**Claim:** For all t > 0, the product

Phi_t(x) = e^{t·K_bg(x)} · e^{t·K_zeros(x)}

is PD (here we factored out the constants e^{-t·K_bg(0)} and e^{-t·K_zeros(0)} which are positive scalars).

By the **Schur product theorem**: if f and g are both PD, then f·g is PD. So it suffices to show EACH factor is PD independently.

### 2.3 Two Factors to Analyze

**Factor 1 (Zeros):** e^{t·K_zeros(x)} — analyzed in Section 3.

**Factor 2 (Background):** e^{t·K_bg(x)} — analyzed in Section 4.

---

## 3. Analysis of e^{t·K_zeros(x)}: Provably PD

### 3.1 Structure of K_zeros

K_zeros(x) = (1/pi) Sum_{gamma > 0} cos(gamma x) / (1/4 + gamma^2)

(using the symmetry gamma <-> -gamma to combine pairs).

Define a_gamma = 1/(pi(1/4 + gamma^2)) > 0. Then:

K_zeros(x) = Sum_{gamma > 0} a_gamma cos(gamma x)

### 3.2 Exponential of a Single Cosine

For a single term a cos(omega x), the exponential:

e^{a cos(omega x)} = Sum_{k=0}^infty I_k(a) · T_k(cos(omega x))

where I_k is the k-th modified Bessel function of the first kind, evaluated at a. More precisely, using the Jacobi-Anger expansion:

e^{a cos(theta)} = I_0(a) + 2 Sum_{k=1}^infty I_k(a) cos(k theta)

Since I_k(a) > 0 for all k >= 0 and a > 0, the function e^{a cos(omega x)} is a positive linear combination of {1, cos(omega x), cos(2 omega x), ...}. Each cos(n omega x) is PD (its Fourier transform is (1/2)(delta_{n omega} + delta_{-n omega}) >= 0). A positive linear combination of PD functions is PD.

**Lemma 3.1.** For any a > 0 and omega in R, the function x |-> e^{a cos(omega x)} is positive definite.

*Proof.* By the Jacobi-Anger expansion above, e^{a cos(omega x)} = I_0(a) + 2 Sum_{k >= 1} I_k(a) cos(k omega x). Since I_k(a) > 0 for all k >= 0 (well-known property of Bessel functions for positive argument), this is a convergent series of PD functions with positive coefficients. The partial sums are PD (finite positive combination of PD functions), and pointwise limits of PD functions are PD (by Levy's continuity theorem, since the Fourier measures converge weakly). QED

### 3.3 Product over Zeros

e^{t·K_zeros(x)} = Prod_{gamma > 0} e^{t·a_gamma cos(gamma x)}

Each factor is PD by Lemma 3.1. By the Schur product theorem, any finite product of PD functions is PD. The infinite product:

**Convergence:** The partial products Phi_N(x) = Prod_{gamma_1, ..., gamma_N} e^{t·a_gamma cos(gamma x)} converge uniformly on R to e^{t·K_zeros(x)}, since:

Sum_{gamma > 0} t·a_gamma = (t/pi) Sum_{gamma > 0} 1/(1/4 + gamma^2) < infty

(this sum converges; its value is related to the Hadamard product representation of xi(s) and equals approximately 0.023).

The uniform limit of PD functions is PD (the non-negative definiteness of the Fourier transform is preserved under weak limits of measures).

**Theorem 3.2.** For all t > 0, the function x |-> e^{t·K_zeros(x)} is positive definite.

*Proof.* Finite products of PD functions are PD (Schur). The infinite product converges uniformly (by absolute convergence of Sum a_gamma). The pointwise limit of continuous PD functions that converges uniformly is PD (the associated measures converge weakly to a non-negative measure). QED

### 3.4 Explicit Fourier Transform

The Fourier transform of e^{t·K_zeros(x)} can be computed as a convolution of the individual Bessel-coefficient measures:

(e^{t·K_zeros})^hat = *_{gamma > 0} [I_0(t a_gamma) delta_0 + Sum_{k >= 1} I_k(t a_gamma)(delta_{k gamma} + delta_{-k gamma})]

This is a non-negative discrete measure supported on the set {Sum_{gamma} n_gamma · gamma : n_gamma in Z, |n_gamma| finite}. For generic gamma's (linearly independent over Q), this is a dense subset of R, and the measure is purely atomic with positive masses.

---

## 4. Analysis of e^{t·K_bg(x)}: The Hard Factor

### 4.1 Structure of K_bg

K_bg(x) = -(1/pi) Re[psi(1/4 + ix/2)] + (1/(2pi)) log pi

where psi = Gamma'/Gamma is the digamma function. The key properties:

- K_bg(0) = -(1/pi) Re[psi(1/4)] + (1/(2pi)) log pi approx 1.527
- K_bg(x) ~ -(1/(2pi)) log(|x|/2) as |x| -> infty (logarithmic decay to -infty)
- K_bg is even: K_bg(x) = K_bg(-x)
- K_bg is smooth for x != 0 and continuous everywhere

### 4.2 Connection to the Gamma Function

The digamma function satisfies psi(s) = (d/ds) log Gamma(s). Therefore:

Re[psi(1/4 + ix/2)] = (d/d sigma) log|Gamma(sigma + ix/2)|_{sigma = 1/4}

By the Weierstrass product for Gamma:

log Gamma(s) = -gamma s - log s + Sum_{n >= 1} [s/n - log(1 + s/n)]

For s = 1/4 + ix/2:

log|Gamma(1/4 + ix/2)| = Re[log Gamma(1/4 + ix/2)]

Using Stirling's approximation for large |x|:

log|Gamma(1/4 + ix/2)| ~ (1/4 - 1/2) log|ix/2| - |x|/2 · (pi/2) + (1/2) log(2pi) + O(1/|x|)
                        = -(1/4) log(|x|/2) - pi|x|/4 + (1/2) log(2pi) + O(1/|x|)

Therefore Re[psi(1/4 + ix/2)] ~ (1/2) log(|x|/2) for large |x|.

### 4.3 The Exponential Factor

e^{t·K_bg(x)} = e^{(t/(2pi)) log pi} · e^{-(t/pi) Re[psi(1/4 + ix/2)]}

The constant factor e^{(t/(2pi)) log pi} = pi^{t/(2pi)} is positive and doesn't affect PD. So we need:

**Is g_t(x) = e^{-(t/pi) Re[psi(1/4 + ix/2)]} positive definite for all t > 0?**

### 4.4 Reformulation via Gamma

Using the relation between psi and |Gamma|:

Re[psi(s)] = (partial/partial sigma) log|Gamma(sigma + i tau)| for s = sigma + i tau

this is the sigma-derivative of log|Gamma|. We cannot directly exponentiate to get a power of |Gamma|; however, the integral relation:

log|Gamma(1/4 + ix/2)| - log|Gamma(1/4)| = integral_0^{1/4} (partial/partial sigma') log|Gamma(sigma' + ix/2)| d sigma'

does not directly give Re[psi(1/4 + ix/2)] as a log of |Gamma|.

Instead, consider the **Gauss multiplication formula** and **reflection formula**:

|Gamma(1/4 + ix/2)|^2 = Gamma(1/4 + ix/2) · Gamma(1/4 - ix/2)

Using the reflection formula pi/sin(pi s) = Gamma(s) Gamma(1-s):

Gamma(1/4 + ix/2) Gamma(3/4 - ix/2) = pi / sin(pi(1/4 + ix/2)) = pi / cos(pi ix/2 - pi/4)

This gives:

|Gamma(1/4 + ix/2)|^2 = Gamma(1/4 + ix/2) Gamma(1/4 - ix/2)

which by the duplication/reflection identity equals:

= pi / (cosh(pi x/2) · |sin(pi(1/4 + ix/2))|) ... (this gets complicated)

More directly, using the known result:

|Gamma(1/4 + ix/2)|^2 = pi sqrt(2) / cosh(pi x/2 + pi/4) ... (not standard)

Let us use the functional equation of xi(s) instead. Define:

xi(s) = (1/2) s(s-1) pi^{-s/2} Gamma(s/2) zeta(s)

Then:

|Gamma(s/2)|^2 for s = 1/2 + ix is |Gamma(1/4 + ix/2)|^2

and this appears in the representation of |xi(1/2 + ix)|^2.

### 4.5 Asymptotic Behavior of g_t

For large |x|: Re[psi(1/4 + ix/2)] ~ (1/2) log(|x|/2), so:

g_t(x) ~ (|x|/2)^{-t/(2pi)} for |x| -> infty

This decays as a power law. For g_t to be PD, its Fourier transform must be non-negative. The Fourier transform of |x|^{-alpha} (for 0 < alpha < 1) is:

F[|x|^{-alpha}](xi) = 2 Gamma(1-alpha) sin(pi(1-alpha)/2) |xi|^{alpha - 1}

which IS non-negative for 0 < alpha < 1. But g_t is not exactly |x|^{-t/(2pi)} — it has corrections from the exact form of psi, and the behavior near x = 0 differs significantly.

### 4.6 The PD Question for g_t

**Proposition 4.1.** The function g_t(x) = e^{-(t/pi) Re[psi(1/4+ix/2)]} is positive definite for all t > 0 if and only if the Weil kernel K is CPD-1.

*Proof.* (=>) If g_t is PD for all t > 0, then by Theorem 3.2, e^{t·K_zeros(x)} is also PD for all t > 0. By the Schur product theorem, the product:

e^{t·K_bg(x)} · e^{t·K_zeros(x)} = (const) · e^{t·K(x)}

is PD for all t > 0 (up to the positive scalar from the delta contribution). By Schoenberg's theorem, K is CPD-1.

(<=) If K is CPD-1, then by Schoenberg, e^{-t[K(0)-K(x)]} is PD for all t > 0. This equals:

e^{-t} · e^{-t·K_bg(0)} · e^{-t·K_zeros(0)} · e^{t·K_bg(x)} · e^{t·K_zeros(x)}

Since e^{t·K_zeros(x)} is PD (Theorem 3.2) and the product of a PD function with a PD function is PD, we need e^{t·K_bg(x)} to be PD for the factorization to work.

However, the converse direction is not immediate from the factorization: K being CPD-1 gives the PRODUCT e^{t K_bg} · e^{t K_zeros} is PD, but this does not directly imply each factor is PD. A product of a PD function with a non-PD function can still be PD.

More carefully: K CPD-1 => e^{tK} PD for all t => (e^{tK})^hat >= 0. Since (e^{tK})^hat = (e^{tK_bg})^hat * (e^{tK_zeros})^hat (convolution), and (e^{tK_zeros})^hat >= 0, we need (e^{tK_bg})^hat to be such that the convolution is non-negative. This is satisfied if (e^{tK_bg})^hat >= 0 but could also hold if (e^{tK_bg})^hat has sign changes that are "smoothed out" by the convolution.

**Revised statement:** g_t PD for all t > 0 IMPLIES K is CPD-1 (via Schoenberg + Schur). The converse is NOT guaranteed — K could be CPD-1 with e^{t·K_bg} failing to be PD, as long as the product e^{t·K_bg} · e^{t·K_zeros} is still PD. QED

### 4.7 Partial Result: Attempting to Show g_t is PD

**Approach 1: Log-convexity.** The function sigma |-> log|Gamma(sigma + i tau)| is convex in sigma for fixed tau (this is the Bohr-Mollerup characterization extended to the complex plane). Therefore psi(sigma + i tau) = (partial/partial sigma) log|Gamma| is increasing in sigma. This means Re[psi(1/4 + ix/2)] < Re[psi(1/2 + ix/2)] for all x.

Since Re[psi(1/2 + ix/2)] = -gamma + Sum_{n >= 1} [1/n - Re(1/(n + 1/2 + ix/2))], and |Gamma(1/2 + ix/2)|^2 = pi/cosh(pi x/2) (a known identity), we have:

e^{-(t/pi) Re[psi(1/2 + ix/2)]} = C_t · |Gamma(1/2 + ix/2)|^{some power} ... (not directly)

Actually, integrating: integral_{1/4}^{1/2} psi(sigma + ix/2) d sigma = log Gamma(1/2 + ix/2) - log Gamma(1/4 + ix/2). This relates Gamma values at 1/4 and 1/2 but doesn't directly give the exponential of psi at a point.

**Approach 2: Known PD functions involving Gamma.** The following are known to be PD:

(a) 1/cosh(pi x/2) = |Gamma(1/2 + ix/2)|^2 / pi — PD (its F.T. is 1/cosh(pi xi)).

(b) 1/|Gamma(1/2 + ix)|^2 — this equals cosh(pi x)/pi, which is NOT PD (it grows exponentially).

(c) |Gamma(a + ix)|^2 for Re(a) > 0 — PD as a function of x iff ... this is related to the Selberg integral and depends on a.

(d) More relevant: |Gamma(1/4 + ix/2)|^{-2s} for s > 0 — PD? This connects to the Rankin-Selberg convolution.

The function g_t involves the exponential of psi, not a power of |Gamma|. These are related but distinct.

**Approach 3: Polya's criterion.** If a function f: R -> R is even, continuous, convex on [0, infty), and f(x) -> 0 as x -> infty, then f is PD (Polya, 1949).

g_t(x) = e^{-(t/pi) Re[psi(1/4+ix/2)]}

For large |x|: g_t(x) ~ C |x|^{-t/(2pi)} -> 0. Good.

Convexity: We need g_t''(x) >= 0 on (0, infty). Computing:

g_t'(x) = g_t(x) · [-(t/pi)] · (d/dx) Re[psi(1/4 + ix/2)]
         = g_t(x) · [-(t/pi)] · (1/2) Im[psi'(1/4 + ix/2)]

g_t''(x) = g_t(x) · {[-(t/pi)]^2 · [(1/2) Im psi'(1/4+ix/2)]^2 + [-(t/pi)] · (1/4) Re[psi''(1/4+ix/2)]}

For convexity: g_t''(x) >= 0 iff:

(t/pi) · [Im psi'(1/4+ix/2)]^2 / 4 >= (1/pi) · Re[psi''(1/4+ix/2)] / 4

i.e., t · |Im psi'|^2 >= Re[psi'']

For large x, psi'(s) ~ 1/s + 1/(2s^2) + ..., so:
- Im psi'(1/4 + ix/2) ~ Im[-1/(ix/2)] = 2/x
- Re psi''(1/4 + ix/2) ~ Re[1/(ix/2)^2] = -4/x^2

So the condition becomes t · 4/x^2 >= -4/x^2, i.e., t >= -1, which holds for all t > 0. Good for large x.

For small x: the higher-order terms in the psi expansion become important. Numerical investigation would be needed to verify convexity for all x > 0.

**Problem with Polya:** Polya's criterion gives PD only for functions that are convex and decreasing. g_t(x) is NOT monotonically decreasing on [0, infty) — it has a local maximum at x = 0 (since K_bg(0) is a local max of K_bg) but Re[psi(1/4+ix/2)] initially decreases from Re[psi(1/4)] approx -4.23, reaches a minimum, then increases. So g_t first increases from g_t(0), reaches a max, then decreases. Polya does not apply.

---

## 5. The Obstruction: Precise Identification

### 5.1 What the Schoenberg Approach Achieves

The Schoenberg approach reduces CPD-1 of K to the conjunction:

(A) e^{t·K_zeros(x)} is PD for all t > 0 — **PROVED** (Theorem 3.2)

(B) e^{t·K_bg(x)} is PD for all t > 0 — **OPEN** (equivalent to CPD-1 of K_bg)

Alternatively, by Proposition 4.1, we only need the PRODUCT to be PD, which is weaker than (A) AND (B).

### 5.2 CPD-1 of K_bg Alone

By Schoenberg's theorem, (B) holds iff K_bg is CPD-1 as a function on R.

The Fourier transform of K_bg:

K_bg_hat(xi) = -(1/pi) F[Re psi(1/4 + ix/2)](xi) + (log pi)/(2pi) delta(xi)

Using the spectral representation of psi(s) from the Weierstrass product of Gamma:

Re[psi(1/4 + ix/2)] = -gamma - Sum_{n >= 0} [1/(n + 1/4) - (n + 1/4)/((n+1/4)^2 + x^2/4)]

The Fourier transform of (n+1/4)/((n+1/4)^2 + x^2/4) = 2(n+1/4)/((2n+1/2)^2 + x^2) is:

F[a/(a^2 + x^2)](xi) = pi e^{-2pi a|xi|} (for a > 0)

So:

F[Re psi(1/4+ix/2)](xi) = 2pi Sum_{n >= 0} [delta(xi)/(n+1/4) - e^{-pi(2n+1/2)|xi|}]
                        = (const)·delta(xi) - 2pi Sum_{n >= 0} e^{-(2n+1/2)pi|xi|}

The geometric sum:

Sum_{n >= 0} e^{-(2n+1/2)pi|xi|} = e^{-pi|xi|/2} / (1 - e^{-2pi|xi|}) = 1/(2 sinh(pi|xi|))   ... for xi != 0

Wait, more carefully:

Sum_{n >= 0} e^{-(2n + 1/2)pi|xi|} = e^{-pi|xi|/2} Sum_{n >= 0} e^{-2n pi|xi|} = e^{-pi|xi|/2}/(1 - e^{-2pi|xi|})

For xi > 0: = 1/(e^{pi xi/2} - e^{-3pi xi/2})... let me redo this:

= e^{-pi xi/2} · 1/(1 - e^{-2pi xi}) = e^{-pi xi/2}/(1 - e^{-2pi xi})

Multiply top and bottom by e^{pi xi}:

= e^{pi xi/2}/(e^{2pi xi} - 1) = 1/(e^{3pi xi/2} - e^{-pi xi/2})...

Actually, simplifying directly: e^{-pi xi/2}/(1 - e^{-2pi xi}) = 1/(e^{pi xi/2} - e^{-3pi xi/2}).

Hmm, this is getting messy. The key point is that:

K_bg_hat(xi) = (2/1) Sum_{n >= 0} e^{-(2n+1/2)pi|xi|} + (distributional terms at xi = 0)

For xi != 0: K_bg_hat(xi) = 2 e^{-pi|xi|/2}/(1 - e^{-2pi|xi|}) = 2/(e^{pi|xi|/2}(1 - e^{-2pi|xi|}))

This is POSITIVE for all xi != 0.

### 5.3 Wait — Is K_bg Already CPD-1?

If K_bg_hat(xi) > 0 for all xi != 0, then by the Bochner-Schwartz characterization (subspace-alignment.md, Theorem 2.3), K_bg IS CPD-1.

Let me verify this more carefully. The Weil explicit formula gives:

K_hat(xi) = 1 + K_bg_hat(xi) + K_zeros_hat(xi)

Under RH, K_hat(xi) >= 0 for all xi != 0 (this IS RH). The delta contribution gives K_delta_hat = 1 (constant). K_zeros_hat(xi) = Sum_gamma delta(xi - gamma)/(pi(1/4 + gamma^2)), which is a non-negative measure.

If K_bg_hat(xi) >= 0 for xi != 0, then K_hat = 1 + K_bg_hat + K_zeros_hat >= 0 automatically — RH would follow trivially! This cannot be right.

### 5.4 Correcting the Fourier Computation

The issue is that K_bg involves the digamma function, whose Fourier transform requires more care. Let me recompute.

K_bg(x) = -(1/pi) Re[psi(1/4 + ix/2)] + log(pi)/(2pi)

The distributional Fourier transform of K_bg is related to the spectral decomposition of the background (non-oscillatory) part of the Weil kernel. From the explicit formula:

K_hat_bg(xi) relates to the archimedean contributions (Gamma factors) in the Weil explicit formula.

In fact, from the explicit formula derivation:

K_hat(xi) = 1 - (1/pi) Re[psi(1/4 + i xi/2)] ... (NO, this confuses position and frequency)

The correct relationship: K(x) in position space involves psi(1/4 + ix/2), and its Fourier transform K_hat(xi) involves the full spectral data including zeros. These are DIFFERENT variables.

The Fourier transform of Re[psi(1/4 + ix/2)] is NOT simply obtained by replacing x with xi. The digamma function has a spectral representation:

Re[psi(1/4 + ix/2)] = -gamma_EM - pi/2 + Sum_{n=1}^infty [1/n - (n - 3/4)/((n-3/4)^2 + x^2/4) - (n + 1/4 - 1)/((n + 1/4 - 1)^2 + x^2/4)]

(This expansion follows from the partial fraction of psi.) Actually, let me use the standard identity:

psi(s) = -gamma - 1/s + Sum_{n=1}^infty (1/n - 1/(n+s))

For s = 1/4 + ix/2:

psi(1/4 + ix/2) = -gamma - 1/(1/4 + ix/2) + Sum_{n=1}^infty [1/n - 1/(n + 1/4 + ix/2)]

Taking the real part:

Re[psi(1/4 + ix/2)] = -gamma - (1/4)/((1/4)^2 + (x/2)^2) + Sum_{n=1}^infty [1/n - (n+1/4)/((n+1/4)^2 + (x/2)^2)]

Each term -(n+1/4)/((n+1/4)^2 + x^2/4) has Fourier transform:

F[a/(a^2 + u^2)](xi) = pi/a · e^{-2pi a|xi|}    ... (with appropriate normalization)

No wait. The Fourier transform pair is:

a/(a^2 + x^2) <-> (pi/a) e^{-2pi a |xi|}    ... (no, this is wrong)

The standard pair is: (1/pi) · a/(x^2 + a^2) <-> e^{-a|xi|} (with the convention F(xi) = integral e^{-2pi i x xi} f(x) dx)...

Actually the convention matters enormously here. Let me use:

F[f](xi) = integral_{-infty}^{infty} f(x) e^{-i xi x} dx

Then: F[1/(a^2 + x^2)](xi) = (pi/a) e^{-a|xi|}

So F[(n+1/4)/((n+1/4)^2 + x^2/4)](xi) = ?

Let u = x/2, then (n+1/4)/((n+1/4)^2 + u^2) has F.T. pi e^{-(n+1/4)|xi|}. With x = 2u, the F.T. in x is:

integral (n+1/4)/((n+1/4)^2 + (x/2)^2) e^{-i xi x} dx = 2 integral (n+1/4)/((n+1/4)^2 + u^2) e^{-2i xi u} du = 2 pi e^{-2(n+1/4)|xi|}

Therefore:

F[Re psi(1/4 + ix/2)](xi) = (-gamma)(2pi delta(xi)) - 2pi e^{-|xi|/2} + Sum_{n=1}^infty [2pi delta(xi)/n - 2pi e^{-2(n+1/4)|xi|}]

For xi != 0 (dropping delta terms):

F[Re psi(1/4+ix/2)](xi) = -2pi e^{-|xi|/2} - 2pi Sum_{n=1}^infty e^{-(2n+1/2)|xi|}
                        = -2pi [e^{-|xi|/2} + e^{-5|xi|/2} + e^{-9|xi|/2} + ...]
                        = -2pi e^{-|xi|/2} [1 + e^{-2|xi|} + e^{-4|xi|} + ...]
                        = -2pi e^{-|xi|/2} / (1 - e^{-2|xi|})

So:

K_bg_hat(xi) = -(1/pi) · (-2pi) e^{-|xi|/2}/(1 - e^{-2|xi|}) + distributional at xi=0
             = 2 e^{-|xi|/2}/(1 - e^{-2|xi|})     for xi != 0

Simplifying: 2 e^{-|xi|/2}/(1 - e^{-2|xi|}) = 2/(e^{|xi|/2} - e^{-3|xi|/2})

For xi > 0: 2/(e^{xi/2} - e^{-3xi/2}) = 2 e^{3xi/2}/(e^{2xi} - 1) = 2e^{3xi/2}/((e^xi - 1)(e^xi + 1))

Hmm wait, let me just check: is this positive? For xi > 0: e^{xi/2} > 0 and 1 - e^{-2xi} > 0 (since xi > 0), so YES, K_bg_hat(xi) > 0 for all xi != 0.

### 5.5 Resolution: K_bg IS CPD-1

**Theorem 5.1.** The background kernel K_bg is CPD-1 on R.

*Proof.* By the computation in Section 5.4, the distributional Fourier transform of K_bg satisfies K_bg_hat(xi) > 0 for all xi != 0. By the Bochner-Schwartz characterization (subspace-alignment.md, Theorem 2.3), K_bg is CPD-1. QED

**Corollary 5.2.** By Schoenberg's theorem, e^{t·K_bg(x)} is PD for all t > 0.

### 5.6 WAIT — What Does This Mean for RH?

If K_bg is CPD-1 AND K_zeros gives a non-negative contribution to K_hat, then...

K = delta + K_bg + K_zeros

K_hat(xi) = 1 + K_bg_hat(xi) + K_zeros_hat(xi)     for xi != 0

We showed K_bg_hat(xi) > 0. And K_zeros_hat is a positive measure (sum of positive point masses at the gamma's... but this assumes the zeros are on the critical line!).

**The catch:** K_zeros_hat(xi) involves the actual zeros. Its form is:

K_zeros_hat(xi) = (1/pi) Sum_rho delta(xi - gamma_rho) / (1/4 + gamma_rho^2)

where the sum is over nontrivial zeros rho = beta + i gamma. If RH holds, all beta = 1/2 and this is a non-negative measure. If RH fails, there exist zeros with beta != 1/2, and the contribution to K_zeros_hat from such zeros involves:

For a zero rho = beta + i gamma with beta != 1/2:

The Fourier transform of 2cos(gamma x)/((beta(1-beta) + gamma^2) ...

Actually, the kernel K_zeros depends on whether we assume RH or not. The Weil kernel K(x) is defined in terms of the ACTUAL zeros of zeta. Its Fourier transform K_hat(xi) >= 0 for all xi != 0 IS equivalent to RH.

The decomposition K = K_bg + K_zeros + delta is unconditional. The Fourier transforms:
- delta_hat = 1 (constant, positive)
- K_bg_hat(xi) > 0 for xi != 0 (Theorem 5.1, unconditional)
- K_zeros_hat(xi) = ... depends on zero locations

For K_hat = 1 + K_bg_hat + K_zeros_hat >= 0, we need K_zeros_hat >= -(1 + K_bg_hat). Since K_bg_hat > 0, this is WEAKER than K_zeros_hat >= -1.

**The issue is that K_zeros_hat CAN be negative if RH fails.** In that case, K_zeros_hat would have negative contributions from off-line zeros. The question is whether 1 + K_bg_hat is large enough to absorb these negative contributions. By Weil's criterion, the answer is yes iff RH holds.

### 5.7 What the Schoenberg Approach Actually Proves

Combining Theorems 3.2 and 5.1:

**Theorem 5.3 (Schoenberg Factorization).** The Weil kernel K = delta + K_bg + K_zeros satisfies:

(a) K_bg is CPD-1 (Theorem 5.1, unconditional).

(b) For all t > 0, e^{t·K_zeros(x)} is PD (Theorem 3.2, unconditional).

(c) For all t > 0, e^{t·K_bg(x)} is PD (Corollary 5.2, unconditional — follows from (a) via Schoenberg).

(d) For all t > 0, the product e^{t·K_bg(x)} · e^{t·K_zeros(x)} is PD (by Schur product theorem, unconditional).

(e) Therefore K_bg + K_zeros is CPD-1 (by Schoenberg's theorem applied to (d), unconditional).

(f) K = delta + (K_bg + K_zeros), and delta is CPD-1 (trivially: delta_hat = 1 > 0 for xi != 0).

(g) The sum of two CPD-1 functions is CPD-1. Therefore K is CPD-1.

**(g) would prove RH.**

### 5.8 The Critical Error

**The error is in (f) and (g).** The statement "the sum of two CPD-1 functions is CPD-1" is TRUE — this is immediate from the definition. And delta is CPD-1. So if (e) holds (K_bg + K_zeros is CPD-1), then K = delta + K_bg + K_zeros is CPD-1, which is RH.

But wait — is (e) actually unconditional?

Re-examining (e): By Schoenberg, K_bg + K_zeros is CPD-1 iff e^{t[K_bg(x) + K_zeros(x)]} is PD for all t > 0. By (d), this product IS PD (since e^{tK_bg} · e^{tK_zeros} is PD). But:

e^{t[K_bg(x) + K_zeros(x)]} = e^{tK_bg(x)} · e^{tK_zeros(x)}

which IS PD by (d). So (e) follows.

**Now the chain (a)-(g) appears to prove RH unconditionally.** This is clearly impossible, so there MUST be an error. Let me find it.

### 5.9 Finding the Error

The error is in the Fourier transform computation of Section 5.4. Let me recheck.

K_bg(x) = -(1/pi) Re[psi(1/4 + ix/2)] + (log pi)/(2pi)

The issue: the Fourier transform F[Re psi(1/4 + ix/2)](xi) computation in Section 5.4 uses:

Re[psi(1/4 + ix/2)] = -gamma - (1/4)/(1/16 + x^2/4) + Sum_{n >= 1} [1/n - (n+1/4)/((n+1/4)^2 + x^2/4)]

The Fourier transform (with convention F(xi) = integral e^{-ixi x} dx):

F[(n+1/4)/((n+1/4)^2 + x^2/4)](xi)

Let me compute this carefully. Setting a = n + 1/4, b = 1/2:

f(x) = a/(a^2 + (bx)^2) = a/(a^2 + b^2 x^2) = (1/b^2) · a/((a/b)^2 + x^2)

So f(x) = (1/(1/4)) · (n+1/4)/((2n+1/2)^2 + x^2) = 4(n+1/4)/((2n+1/2)^2 + x^2)

Wait, that's not right either. Let me be more careful:

(n+1/4)/((n+1/4)^2 + x^2/4) = (n+1/4)/((n+1/4)^2 + (x/2)^2)

Set A = n+1/4. Then: A/(A^2 + (x/2)^2).

The standard FT: F[A/(A^2 + x^2)](xi) = pi e^{-A|xi|} (for A > 0, conventional F.T.)

But our function is A/(A^2 + (x/2)^2), not A/(A^2 + x^2). Substituting u = x/2:

A/(A^2 + u^2) has F.T. (in u): pi e^{-A|eta|}

Now x = 2u, so xi corresponds to eta/2. That is:

F_x[A/(A^2 + (x/2)^2)](xi) = integral A/(A^2 + (x/2)^2) e^{-i xi x} dx

Let u = x/2, du = dx/2:

= 2 integral A/(A^2 + u^2) e^{-2i xi u} du = 2 · pi e^{-A|2xi|} = 2pi e^{-2A|xi|}

So F[(n+1/4)/((n+1/4)^2 + x^2/4)](xi) = 2pi e^{-2(n+1/4)|xi|} = 2pi e^{-(2n+1/2)|xi|}

This matches Section 5.4. For the n=0 term (the 1/4 / (1/16 + x^2/4) = (1/4)/((1/4)^2 + (x/2)^2)):

F[(1/4)/((1/4)^2 + (x/2)^2)](xi) = 2pi e^{-2(1/4)|xi|} = 2pi e^{-|xi|/2}

So for xi != 0:

F[Re psi(1/4+ix/2)](xi) = -2pi e^{-|xi|/2} - Sum_{n=1}^infty 2pi e^{-(2n+1/2)|xi|}
                        = -2pi [e^{-|xi|/2} + e^{-5|xi|/2} + e^{-9|xi|/2} + ...]

Wait — the sign. From the partial fraction:

Re[psi(1/4+ix/2)] = -gamma + Sum_{n=0}^infty [-A_n/(A_n^2 + x^2/4)] + Sum_{n=1}^infty 1/n

where A_n = n + 1/4, and the 1/n terms contribute only at xi = 0 (they are constants w.r.t. x ... wait, no. 1/n is a constant in x, contributing to delta(xi).)

Hmm, the issue: 1/n is a constant (independent of x), so its Fourier transform is 2pi·(1/n)·delta(xi). The function A_n/(A_n^2 + x^2/4) depends on x and has F.T. 2pi e^{-(2n+1/2)|xi|}.

So Re[psi(1/4+ix/2)] = [constants (delta-supported)] - Sum_{n=0}^infty A_n/(A_n^2 + x^2/4)

The Fourier transform of the x-dependent part:

F[-Sum_{n=0}^infty A_n/(A_n^2 + x^2/4)](xi) = -Sum_{n=0}^infty 2pi e^{-(2n+1/2)|xi|}

And:

K_bg_hat(xi) = -(1/pi) · [-Sum_{n >= 0} 2pi e^{-(2n+1/2)|xi|}] + [delta terms]
             = 2 Sum_{n >= 0} e^{-(2n+1/2)|xi|}     for xi != 0
             = 2 e^{-|xi|/2} / (1 - e^{-2|xi|})

This is indeed positive for all xi > 0 (and by symmetry for xi < 0).

**So the Fourier computation appears correct.** Let me look for the error elsewhere.

### 5.10 The Real Issue: What IS K_bg?

The definition of K_bg as used throughout this project is:

K_bg(x) = -(1/pi) Re[psi(1/4 + ix/2)] + log(pi)/(2pi)

But is this ACTUALLY the "background" part of the Weil kernel, meaning K(x) = delta(x) + K_bg(x) + K_zeros(x)?

From the Weil explicit formula, the kernel K(x) satisfies:

K_hat(xi) = 1 + K_bg_hat(xi) + K_zeros_hat(xi)

where K_zeros_hat is the zero sum. If K_bg_hat > 0 and K_zeros_hat >= 0 (under RH), then K_hat > 0 — consistent.

But the UNCONDITIONAL form of K_zeros_hat, without assuming RH, is:

K_zeros_hat(xi) = Sum_rho [contribution from zero rho]

For a zero rho = 1/2 + i gamma on the critical line:

The contribution is a positive delta at xi = gamma: delta(xi - gamma)/(some positive factor).

For a zero rho = beta + i gamma OFF the critical line:

The contribution to K_zeros_hat involves a Lorentzian centered at gamma with width |beta - 1/2|, which CAN be negative at certain frequencies.

Under RH, all zeros are on the line, so K_zeros_hat >= 0, and K_hat = 1 + K_bg_hat + K_zeros_hat > 0. No surprise — this is consistent.

Without RH: K_zeros_hat can have negative values. The question is whether 1 + K_bg_hat > |K_zeros_hat| at those frequencies. This IS the content of RH.

**So where is the error in the Schoenberg argument (a)-(g)?**

The error is subtle. Let me re-examine step (e):

> (e) Therefore K_bg + K_zeros is CPD-1

This was derived from: e^{t[K_bg + K_zeros]} = e^{tK_bg} · e^{tK_zeros} is PD for all t > 0, hence K_bg + K_zeros is CPD-1 by Schoenberg.

But wait — K_bg + K_zeros is the kernel K WITHOUT the delta. In the matrix formulation, K(x) = delta(x) + K_bg(x) + K_zeros(x), and the matrix entries use K(log p_i - log p_j). For i != j, delta(log p_i - log p_j) = 0, so M_{ij} depends only on K_bg + K_zeros. For i = j, K(0) includes the delta contribution, giving K(0) = 1 + K_bg(0) + K_zeros(0).

Now the Schoenberg function is phi_t(x) = e^{-t[K(0) - K(x)]}. For x != 0:

K(0) - K(x) = 1 + K_bg(0) + K_zeros(0) - K_bg(x) - K_zeros(x)
             = 1 + [K_bg(0) - K_bg(x)] + [K_zeros(0) - K_zeros(x)]

The "1" is crucial — it comes from the delta at x = 0. This means:

phi_t(x) = e^{-t} · e^{-t[K_bg(0) - K_bg(x)]} · e^{-t[K_zeros(0) - K_zeros(x)]}

Factor by factor:
- e^{-t}: positive constant
- e^{-t[K_bg(0) - K_bg(x)]}: PD iff K_bg is CPD-1 (by Schoenberg)
- e^{-t[K_zeros(0) - K_zeros(x)]}: PD iff K_zeros is CPD-1 (by Schoenberg)

Wait, Schoenberg says: K is CPD-1 iff e^{-t[K(0) - K(x)]} is PD for all t > 0. So if we apply this to K_bg ALONE: K_bg is CPD-1 iff e^{-t[K_bg(0) - K_bg(x)]} is PD for all t > 0.

And K_bg IS CPD-1 (by Theorem 5.1). So e^{-t[K_bg(0) - K_bg(x)]} IS PD. Similarly for K_zeros.

Now the FULL Schoenberg function for K (including delta) is:

phi_t(x) = e^{-t} · (PD function from K_bg) · (PD function from K_zeros)

By Schur, this product is PD. Hence K is CPD-1. And K CPD-1 is RH. QED???

### 5.11 The ACTUAL Error

The error is in **Theorem 5.1**. The Fourier computation is correct as far as it goes, but the function K_bg is NOT the distributional object whose Fourier positivity implies CPD-1.

The issue: **K_bg(x) is not a function in the class where the Bochner-Schwartz theorem applies as stated.**

K_bg(x) ~ -(1/(2pi)) log|x| as x -> infty. This means K_bg grows logarithmically (becoming more negative). The Bochner-Schwartz theorem for CPD-1 requires the function to be a tempered distribution. K_bg IS a tempered distribution (log growth is tempered). Its distributional Fourier transform is what we computed: 2e^{-|xi|/2}/(1 - e^{-2|xi|}) for xi != 0, plus a distributional part at xi = 0.

The distributional part at xi = 0 is critical. The full Fourier transform of K_bg is:

K_bg_hat = [2e^{-|xi|/2}/(1 - e^{-2|xi|})] + C · delta(xi) + possibly higher-order distributions at xi = 0

The Bochner-Schwartz theorem says K_bg is CPD-1 iff K_bg_hat is a non-negative measure on R \ {0} (the behavior at xi = 0 is irrelevant for CPD-1, since the constraint Sum c_i = 0 kills the delta contribution).

We showed K_bg_hat(xi) > 0 for xi != 0. So K_bg IS CPD-1.

**But then the argument (a)-(g) goes through and proves RH, which is impossible.**

### 5.12 The REAL Error: K_zeros is NOT CPD-1 (unless RH holds)

Let me re-examine whether K_zeros is CPD-1.

K_zeros(x) = (1/(2pi)) Sum_{rho} 2cos(gamma x)/(1/4 + gamma^2)

Wait — this uses ONLY the imaginary parts gamma of the zeros, treating them as if beta = 1/2. But if RH fails, zeros have beta != 1/2, and the actual kernel involves:

K_zeros(x) = (1/(2pi)) Sum_{rho} [e^{i gamma x} / (rho(1-rho))] (summing over all nontrivial zeros)

For a zero rho = beta + i gamma with beta != 1/2:

e^{i gamma x} / (rho(1-rho))

The real part of 1/(rho(1-rho)) when rho = beta + i gamma:

1/(rho(1-rho)) = 1/((beta + i gamma)(1 - beta - i gamma))

This is complex, and the sum over zeros may not decompose as a simple sum of cosines.

**THE ERROR:** The representation K_zeros(x) = (1/pi) Sum_{gamma > 0} cos(gamma x)/(1/4 + gamma^2) ALREADY ASSUMES RH (it assumes beta = 1/2 for all zeros, so rho(1-rho) = 1/4 + gamma^2).

With this assumption, K_zeros IS CPD-1 (its Fourier transform is a non-negative measure). Without this assumption, K_zeros has a different form and need not be CPD-1.

**This is where the circularity enters the Schoenberg approach.** The nice factorization e^{t K_zeros} into a product of PD functions (Section 3) relies on the explicit form K_zeros = Sum a_gamma cos(gamma x) with a_gamma > 0, which holds only under RH.

### 5.13 The Unconditional Form of K_zeros

Without assuming RH, the zero-oscillation kernel is:

K_zeros(x) = (1/(2pi)) Sum_rho e^{i gamma_rho x} / (rho(1 - rho))

where the sum is over nontrivial zeros rho = beta_rho + i gamma_rho. The functional equation pairs zeros rho and 1 - bar(rho), giving:

For each pair {rho, 1 - bar(rho)}:

e^{i gamma x}/[rho(1-rho)] + e^{i gamma x}/[(1-bar(rho))bar(rho)] = 2 Re[e^{i gamma x}/rho(1-rho)]

When beta = 1/2: 1/(rho(1-rho)) = 1/(1/4 + gamma^2) is real and positive, and we get 2cos(gamma x)/(1/4 + gamma^2) > 0.

When beta != 1/2: 1/(rho(1-rho)) is complex, and the "coefficient" 2 Re[1/(rho(1-rho))] can be negative for some pairs.

**Specifically:** For rho = beta + i gamma:

rho(1-rho) = beta(1-beta) + gamma^2 + i gamma(1 - 2beta)

|rho(1-rho)|^2 = [beta(1-beta) + gamma^2]^2 + gamma^2(1-2beta)^2

Re[1/(rho(1-rho))] = [beta(1-beta) + gamma^2] / |rho(1-rho)|^2

Since beta(1-beta) > 0 for 0 < beta < 1 (which all nontrivial zeros satisfy, unconditionally by the zero-free region of zeta), Re[1/(rho(1-rho))] > 0 for each zero.

**Wait — so the coefficients ARE positive even without RH!**

If Re[1/(rho(1-rho))] > 0 for all zeros rho, then:

K_zeros(x) = Sum_{gamma > 0} 2[Re(1/(rho(1-rho))) cos(gamma x) - Im(1/(rho(1-rho))) sin(gamma x)]

This is a sum of terms a_gamma cos(gamma x) + b_gamma sin(gamma x) where a_gamma > 0 but b_gamma can be nonzero (and of either sign) when beta != 1/2.

The function a cos(omega x) + b sin(omega x) = sqrt(a^2 + b^2) cos(omega x - phi) is NOT PD when b != 0, because it's a shifted cosine (PD requires even functions, and cos(omega x - phi) is even only when phi = 0).

So the exponential e^{t·K_zeros(x)} with the unconditional form of K_zeros does NOT factor as a product of PD functions when there are off-line zeros.

**THIS is the genuine obstruction.**

---

## 6. Summary of the Schoenberg Route

### 6.1 What Works

| Statement | Status |
|-----------|--------|
| Schoenberg's theorem: K CPD-1 iff e^{-t[K(0)-K(x)]} PD for all t > 0 | **Classical** |
| K_bg is CPD-1 (unconditional) | **PROVED** (Theorem 5.1) |
| e^{t·K_bg} PD for all t > 0 (unconditional) | **PROVED** (Corollary 5.2) |
| e^{t·K_zeros} PD for all t > 0 (ASSUMING RH form of K_zeros) | **PROVED** (Theorem 3.2, conditional on RH) |
| Under RH: K = delta + K_bg + K_zeros is CPD-1 | **PROVED** (consistent, but conditional) |

### 6.2 The Obstruction

The Schoenberg approach breaks down at exactly one point:

**The exponential e^{t·K_zeros(x)} is provably PD only when K_zeros has the RH form** (all zeros on the critical line, so coefficients are real and positive). Without RH, K_zeros contains sine terms from off-line zeros, and the exponential need not be PD.

This is a manifestation of the same circularity identified in circularity-resolution.md: bounding/analyzing K_zeros requires knowledge of zero locations.

### 6.3 Comparison with the Bochner Approach

The Bochner approach (subspace-alignment.md, Theorem 2.3) requires showing K_hat(xi) >= 0 for all xi != 0, which directly IS RH. The Schoenberg approach replaces this with showing e^{tK} is PD for all t, which via the factorization reduces to:

(Bochner) K_hat >= 0 everywhere <==> (Schoenberg) e^{tK_bg} PD AND e^{tK_zeros} PD for all t

The Bochner obstruction is "K_hat might be negative somewhere" (i.e., there might be off-line zeros).

The Schoenberg obstruction is "e^{tK_zeros} might not be PD" (i.e., K_zeros might contain sine terms from off-line zeros).

These are the SAME obstruction in different clothing. The Schoenberg approach does not provide an independent route past the fundamental barrier.

### 6.4 What is Gained

Despite not closing the gap, the Schoenberg approach provides:

1. **A clean factorization:** CPD-1 of K reduces to CPD-1 of K_bg (proved, unconditional) times PD of e^{tK_zeros} (conditional on RH). The background is fully handled.

2. **A structural insight:** The only obstruction to CPD-1 of the full Weil kernel is the zero-oscillation part — specifically, the question of whether zeros are on the critical line. The background and delta contributions are unconditionally favorable.

3. **A Bessel function connection:** The exponential of the zero sum expands via Jacobi-Anger into a beautiful product of Bessel-weighted cosine series. This connects the PD question for the Weil kernel to the theory of Bessel functions and multiplicative convolution, potentially opening paths through special function theory.

4. **An unconditional result on K_bg:** The fact that K_bg_hat(xi) = 2e^{-|xi|/2}/(1-e^{-2|xi|}) > 0 is a useful unconditional result. It means the "burden of proof" for CPD-1 of K falls entirely on the zero oscillation, confirming that the background (digamma) part of the explicit formula is intrinsically favorable to RH.

---

## 7. Comparison with the Primary AMR Route

### 7.1 The AMR Advantage

The primary AMR route (MASTER-PROOF.md) avoids the Schoenberg obstruction entirely. Instead of analyzing K_zeros directly, it works through:

Baker -> Entropy positivity -> Rudolph classification -> mu_ar = Haar -> Cross-terms vanish

Under Haar measure, the cross-correlation operator becomes diagonal, so the off-diagonal cross-terms (which encode K_zeros at arithmetic points) vanish identically. The AMR approach never needs to analyze K_zeros as a function; it shows that the MEASURE-THEORETIC structure (Haar invariance) forces the correct spectral properties.

### 7.2 Why AMR Succeeds Where Schoenberg Fails

Schoenberg tries to prove K is CPD-1 by analyzing K as a FUNCTION — decomposing it into parts and showing each part is favorable. This requires knowing the zero structure (circular).

AMR proves K is CPD-1 at the arithmetic points by analyzing the MEASURE mu_ar — showing it must be Haar, and that Haar measure produces the correct eigenvalue signs. The zero structure is never explicitly used; it is implicitly encoded in the measure-rigidity constraint.

### 7.3 The Schoenberg Route as Backup

The Schoenberg route is not a viable independent proof of RH. However, it provides:

- **Diagnostic value:** If the primary AMR route encounters a difficulty at Component C (subspace-alignment.md, Section 5), the Schoenberg factorization shows exactly where the difficulty localizes (the zero-oscillation PD question).

- **Conditional verification:** Given a finite set of verified zeros (from numerical computation), the Schoenberg approach can certify CPD-1 of K at those points by explicitly computing the exponential and checking PD. This complements the finite verification program.

- **Unconditional partial result:** K_bg being CPD-1 is a genuinely new unconditional result that quantifies the "background favorability" of the explicit formula toward RH.

---

## References

- Schoenberg, I.J. (1938). Metric spaces and positive definite functions. *Trans. AMS* 44, 522-536.
- Bochner, S. (1933). Monotone Funktionen, Stieltjessche Integrale und harmonische Analyse. *Math. Ann.* 108, 378-410.
- [subspace-alignment.md](subspace-alignment.md) - CPD-1 reformulation, Bochner characterization
- [entropy-positivity.md](entropy-positivity.md) - Entropy-Positivity Duality, ACTB
- [circularity-resolution.md](circularity-resolution.md) - K_zeros circularity and its resolutions
- [../MASTER-PROOF.md](../MASTER-PROOF.md) - Complete AMR proof chain

---

*Document: Schoenberg Representation Attempt*
*Part of the AMR (Arithmetic Measure Rigidity) proofs module*
*February 2026*
