# Classical Approaches to the Riemann Hypothesis: A Comprehensive Survey

## Table of Contents

1. [Introduction](#1-introduction)
2. [Riemann's Original Work (1859)](#2-riemanns-original-work-1859)
3. [Hardy's 1914 Result](#3-hardys-1914-result)
4. [Hardy-Littlewood and Selberg: Positive Proportion](#4-hardy-littlewood-and-selberg-positive-proportion)
5. [Levinson's Method (1974)](#5-levinsons-method-1974)
6. [Conrey's Result (1989)](#6-conreys-result-1989)
7. [Zero-Free Regions](#7-zero-free-regions)
8. [The Lindelof Hypothesis](#8-the-lindelof-hypothesis)
9. [Li's Criterion](#9-lis-criterion)
10. [The Nyman-Beurling-Baez-Duarte Criterion](#10-the-nyman-beurling-baez-duarte-criterion)
11. [Weil's Explicit Formula](#11-weils-explicit-formula)
12. [The de Bruijn-Newman Constant](#12-the-de-bruijn-newman-constant)
13. [The Selberg Class](#13-the-selberg-class)
14. [Moment Conjectures](#14-moment-conjectures)
15. [Key Barriers and Fundamental Obstacles](#15-key-barriers-and-fundamental-obstacles)
16. [References](#16-references)

---

## 1. Introduction

The Riemann Hypothesis (RH) asserts that all non-trivial zeros of the Riemann zeta function

$$\zeta(s) = \sum_{n=1}^{\infty} \frac{1}{n^s}, \quad \Re(s) > 1,$$

lie on the critical line $\Re(s) = 1/2$. Equivalently, if $\zeta(\rho) = 0$ with $0 < \Re(\rho) < 1$, then $\Re(\rho) = 1/2$.

This document surveys all major classical approaches that have been attempted to prove or disprove the Riemann Hypothesis, from Riemann's 1859 memoir through modern developments. We aim for mathematical rigor and completeness, providing a roadmap for researchers seeking to understand the full landscape of attack strategies.

**Notation.** Throughout, $s = \sigma + it$ where $\sigma, t \in \mathbb{R}$. We write $\rho = \beta + i\gamma$ for a non-trivial zero of $\zeta(s)$. The number of zeros with $0 < \gamma \leq T$ is denoted $N(T)$, and those on the critical line $N_0(T)$.

---

## 2. Riemann's Original Work (1859)

### 2.1 The Memoir "Ueber die Anzahl der Primzahlen unter einer gegebenen Grosse"

Riemann's 1859 paper, only eight pages long, is arguably the most consequential short paper in the history of mathematics. In it, Riemann:

1. Extended the zeta function to a meromorphic function on all of $\mathbb{C}$.
2. Derived the functional equation.
3. Established the explicit formula connecting primes to zeros.
4. Stated the Riemann Hypothesis.
5. Sketched ideas toward a proof strategy.

### 2.2 Analytic Continuation and the Functional Equation

Riemann showed that $\zeta(s)$, initially defined for $\Re(s) > 1$, extends to a meromorphic function on $\mathbb{C}$ with a single simple pole at $s = 1$ with residue 1.

**The Functional Equation.** Define the completed zeta function (the xi function):

$$\xi(s) = \frac{1}{2} s(s-1) \pi^{-s/2} \Gamma(s/2) \zeta(s).$$

Then $\xi$ is an entire function of order 1, and satisfies:

$$\xi(s) = \xi(1-s).$$

Equivalently, one can write the functional equation in its asymmetric form:

$$\zeta(s) = 2^s \pi^{s-1} \sin\left(\frac{\pi s}{2}\right) \Gamma(1-s) \zeta(1-s).$$

**Riemann's derivation** proceeded via the integral representation. Starting from the Gamma function relation $\Gamma(s/2) \pi^{-s/2} n^{-s} = \int_0^\infty x^{s/2 - 1} e^{-n^2 \pi x} \, dx$, he obtained:

$$\pi^{-s/2} \Gamma(s/2) \zeta(s) = \int_0^\infty x^{s/2 - 1} \psi(x) \, dx,$$

where $\psi(x) = \sum_{n=1}^{\infty} e^{-n^2 \pi x}$ is related to the Jacobi theta function $\theta(x) = 1 + 2\psi(x)$. The Jacobi transformation formula $\theta(1/x) = \sqrt{x} \, \theta(x)$ then yields the functional equation.

### 2.3 The Explicit Formula

Riemann derived the remarkable connection between prime numbers and the zeros of $\zeta(s)$. In its modern form (due to von Mangoldt), the explicit formula reads:

$$\psi(x) = \sum_{n \leq x} \Lambda(n) = x - \sum_{\rho} \frac{x^{\rho}}{\rho} - \log(2\pi) - \frac{1}{2} \log(1 - x^{-2}),$$

where $\Lambda$ is the von Mangoldt function, $\psi$ is the Chebyshev function, and the sum runs over the non-trivial zeros $\rho$ of $\zeta(s)$, counted with multiplicity, with the sum understood in the symmetric sense $\lim_{T \to \infty} \sum_{|\Im(\rho)| < T}$.

Riemann's own version concerned the prime counting function $\pi(x)$:

$$\pi_0(x) = \operatorname{Li}(x) - \sum_{\rho} \operatorname{Li}(x^{\rho}) - \log 2 + \int_x^{\infty} \frac{dt}{t(t^2-1)\log t},$$

where $\pi_0(x) = \frac{1}{2}(\pi(x^+) + \pi(x^-))$ and $\operatorname{Li}(x) = \int_2^x \frac{dt}{\log t}$.

The key insight: **the error term in the Prime Number Theorem is controlled by the location of the zeros**. If all zeros satisfy $\Re(\rho) = 1/2$, then:

$$\psi(x) = x + O(x^{1/2} \log^2 x),$$

which is the best possible error term (up to the logarithmic factor). Conversely, any zero with $\Re(\rho) = \beta > 1/2$ would produce oscillatory terms of order $x^{\beta}$ in the error.

### 2.4 Riemann's Proof Strategy

Riemann's paper contains an often-overlooked passage where he discusses a potential approach to the hypothesis. He considered the function:

$$\Xi(t) = \xi(1/2 + it),$$

which is real for real $t$, and the Riemann Hypothesis is equivalent to $\Xi(t)$ having only real zeros.

Riemann indicated he believed one could prove this via a Hadamard product representation. He wrote that $\Xi(t)$ could be expressed as:

$$\Xi(t) = \Xi(0) \prod_{\alpha} \left(1 - \frac{t^2}{\alpha^2}\right),$$

where the product is over the zeros $\alpha$ of $\Xi$. If one could show that this product converges and represents $\Xi$, and if one could establish certain growth conditions, then the reality of the zeros would follow.

Riemann also mentioned a formula for $N(T)$ (later rigorously proved by von Mangoldt):

$$N(T) = \frac{T}{2\pi} \log \frac{T}{2\pi} - \frac{T}{2\pi} + \frac{7}{8} + S(T) + O(1/T),$$

where $S(T) = \frac{1}{\pi} \arg \zeta(1/2 + iT) = O(\log T)$.

### 2.5 The Product Formula and Counting Formula

The Hadamard product for $\xi(s)$ was later rigorously established by Hadamard (1893):

$$\xi(s) = e^{A+Bs} \prod_{\rho} \left(1 - \frac{s}{\rho}\right) e^{s/\rho},$$

where $A, B$ are constants and the product is over the non-trivial zeros. Hadamard showed $B = -\sum_\rho \Re(1/\rho) = \frac{1}{2}\log(4\pi) - 1 - \gamma/2$ where $\gamma$ is Euler's constant.

---

## 3. Hardy's 1914 Result

### 3.1 Statement

**Theorem (Hardy, 1914).** *The Riemann zeta function has infinitely many zeros on the critical line $\Re(s) = 1/2$.*

This was the first unconditional result showing that the critical line is special. Previously, it was conceivable (though not expected) that all zeros could lie off the critical line.

### 3.2 Method of Proof

Hardy's proof is a beautiful application of the theory of integral transforms and the properties of the $\Xi$ function.

**Step 1.** Define $\Xi(t) = \xi(1/2 + it)$, which is real for real $t$ and is an even entire function of $t$.

**Step 2.** Consider the integral:

$$I(T) = \int_0^T \Xi(t) \, dt.$$

Hardy showed that if $\Xi(t)$ has only finitely many real zeros, then one can derive contradictory growth estimates for this integral.

**Step 3.** The key identity Hardy used involves the function:

$$Z(t) = e^{i\theta(t)} \zeta(1/2 + it),$$

where $\theta(t) = \arg(\pi^{-it/2} \Gamma(1/4 + it/2))$ is the Riemann-Siegel theta function. The function $Z(t)$ is real for real $t$, and $|Z(t)| = |\zeta(1/2 + it)|$. A change of sign of $Z(t)$ therefore corresponds to a zero of $\zeta$ on the critical line.

**Step 4.** Hardy considered the integral

$$\int_1^\infty \Xi(t) \Phi(t) \, dt$$

for suitable test functions $\Phi$, and showed via a Mellin transform argument that this integral is nonzero, which forces $\Xi(t)$ to have sign changes, hence real zeros.

More specifically, Hardy used the relation:

$$\int_0^\infty \Xi(t) \cos(xt) \, dt = \frac{\pi}{2}\left(e^{x/2} - 2e^{-x/2}\psi(e^{-2x})\right) \cdot e^{-x/2},$$

where the right side can be analyzed directly to show it is not identically zero for large values, forcing sign changes.

### 3.3 Limitations

Hardy's proof shows infinitely many zeros on the critical line but gives no density result. It cannot tell us what proportion of zeros lie on the critical line. This deficiency motivated the work of Hardy-Littlewood and Selberg.

---

## 4. Hardy-Littlewood and Selberg: Positive Proportion

### 4.1 Hardy-Littlewood (1921)

Hardy and Littlewood refined Hardy's method to obtain a quantitative result:

**Theorem (Hardy-Littlewood, 1921).** *$N_0(T) \geq c T$ for some (unspecified, small) constant $c > 0$ and all sufficiently large $T$.*

Here $N_0(T)$ counts zeros on the critical line up to height $T$. Since $N(T) \sim \frac{T}{2\pi}\log T$, this shows that the critical line zeros have positive density relative to $T$ but not relative to $N(T)$.

Their method involved computing the mean value:

$$\int_0^T |Z(t)|^2 \, dt \sim T \log T,$$

and comparing it to:

$$\int_0^T Z(t)^2 \, dt,$$

using Cauchy-Schwarz to derive a lower bound on the number of sign changes of $Z(t)$.

### 4.2 Selberg's Breakthrough (1942)

Selberg achieved the definitive strengthening:

**Theorem (Selberg, 1942).** *$N_0(T) \geq c \cdot N(T)$ for some positive constant $c > 0$, i.e., a positive proportion of all non-trivial zeros lie on the critical line.*

Selberg's innovation was to consider **mollified** moments. Rather than looking at $Z(t)$ directly, he considered $Z(t) \cdot M(t)$ where $M(t)$ is a "mollifier" -- a short Dirichlet polynomial chosen to smooth out the erratic behavior of $Z(t)$.

**The mollifier idea.** Let:

$$M(s) = \sum_{n \leq y} \frac{\mu(n) a(n)}{n^s}$$

where $y = T^{\theta}$ for some small $\theta > 0$, $\mu$ is the Mobius function, and the $a(n)$ are suitable weights. The key is that $\zeta(s) \cdot M(s) \approx 1$ when $s$ is not near a zero, so $M$ "cancels" the large peaks of $1/\zeta$.

Selberg showed that by computing:

$$\int_0^T |Z(t) M(1/2 + it)|^2 \, dt \quad \text{and} \quad \int_0^T Z(t) M(1/2 + it) \, dt,$$

and applying the Cauchy-Schwarz inequality, one obtains a lower bound for $N_0(T)$ proportional to $N(T)$.

### 4.3 Selberg's Unpublished Improvement

Selberg reportedly obtained $N_0(T) \geq c \cdot N(T)$ with $c$ close to $1/2$ but never published the optimal constant from his method. His approach is the foundation upon which all subsequent improvements rest.

---

## 5. Levinson's Method (1974)

### 5.1 The Main Result

**Theorem (Levinson, 1974).** *At least one-third of the non-trivial zeros of $\zeta(s)$ lie on the critical line. More precisely, $N_0(T) \geq \frac{1}{3} N(T)$ for all sufficiently large $T$.*

This was a dramatic improvement over Selberg's unquantified constant.

### 5.2 Key Innovations

Levinson's breakthrough involved two main ideas:

**Innovation 1: Using the derivative.** Instead of mollifying $\zeta(s)$ directly, Levinson worked with a linear combination:

$$F(s) = \zeta(s) + \frac{1}{\log T} \zeta'(s).$$

The crucial observation is that $\zeta'(s)/\zeta(s) = -\sum_\rho \frac{1}{s - \rho} + \ldots$, so near a zero $\rho$, the derivative $\zeta'(\rho)$ captures information about the multiplicity and location of the zero in a way that $\zeta$ alone does not.

More precisely, Levinson considered:

$$A(s) = \zeta(s) M(s) + \frac{a}{\log T} (\zeta(s) M(s))',$$

where $a$ is a parameter to be optimized and $M(s)$ is a mollifier of length $T^{4/7 - \varepsilon}$.

**Innovation 2: Longer mollifiers.** Levinson's argument required mollifiers of length $y = T^{\theta}$ with $\theta$ close to $1/2$. The key estimate is:

$$\int_0^T |A(1/2 + it)|^2 \, dt \sim (1 + O(\delta)) T \log T,$$

where $\delta$ measures the deviation from 1 and can be made small.

### 5.3 Proof Sketch

**Step 1.** Show that $\zeta(s) M(s)$ has approximately the same zeros as $\zeta(s)$ in the critical strip (since $M(s)$ has no zeros there for appropriate choices).

**Step 2.** Using the argument principle (specifically, a variant of Littlewood's lemma applied to a rectangle), show that the number of zeros of $F(s) = A(s)$ near the critical line is related to mean-value integrals.

**Step 3.** Compute the mean values:

$$\frac{1}{T} \int_0^T |A(1/2 + it)|^2 \, dt$$

and

$$\frac{1}{T} \int_0^T |A(\sigma_0 + it)|^2 \, dt,$$

where $\sigma_0 = 1/2 + R/\log T$ for a suitable constant $R$.

**Step 4.** If the integral on $\sigma = \sigma_0$ is much smaller than the integral on $\sigma = 1/2$, then by a convexity-type argument, most zeros of $A(s)$ must be close to the critical line.

**Step 5.** By optimizing the parameter $a$ and the mollifier length, Levinson shows that the proportion of zeros forced onto (or very near) the critical line is at least $1/3$.

### 5.4 The Fraction 1/3

The fraction $1/3$ arises from the optimization. With mollifiers of length $y = T^{\theta}$, the method gives a proportion of approximately $1 - 6\theta/(3 + 2\theta)$, which is optimized at $\theta$ near $1/2$ (limited by the ability to estimate the mean values), yielding a proportion approaching $1/3$.

---

## 6. Conrey's Result (1989)

### 6.1 The Main Result

**Theorem (Conrey, 1989).** *At least two-fifths (40%) of the non-trivial zeros of $\zeta(s)$ lie on the critical line: $N_0(T) \geq 0.4 \cdot N(T)$.*

More precisely, Conrey proved $N_0(T) \geq 0.4088 \cdot N(T)$ for sufficiently large $T$.

### 6.2 The Refined Mollifier Technique

Conrey's improvement over Levinson came from a more sophisticated choice of mollifier and more powerful mean-value estimates.

**Longer mollifiers.** While Levinson used mollifiers of length $T^{1/2 - \varepsilon}$, Conrey extended this to $T^{4/7}$ by exploiting Deshouillers-Iwaniec's estimates for mean values of Dirichlet polynomials.

**Kloosterman sum estimates.** The key technical innovation was the use of bounds for Kloosterman sums. When computing the mean values of products of $\zeta$ with longer Dirichlet polynomials, one encounters sums of the form:

$$S(m,n;c) = \sum_{\substack{a \bmod c \\ (a,c)=1}} e\left(\frac{ma + n\bar{a}}{c}\right),$$

where $\bar{a}$ denotes the multiplicative inverse of $a$ modulo $c$. Weil's bound $|S(m,n;c)| \leq d(c) \sqrt{(m,n,c)} \sqrt{c}$ (a consequence of the Riemann Hypothesis for curves over finite fields) is essential.

### 6.3 The Asymptotic Large Sieve

Conrey used a result of Deshouillers and Iwaniec on the asymptotic behavior of mean values:

$$\int_0^T \left|\sum_{n \leq N} a_n n^{-it}\right|^2 dt = (T + O(N^{1+\varepsilon})) \sum_{n \leq N} |a_n|^2,$$

but with refined error terms that allow $N$ to be as large as $T^{4/7}$ (rather than $T^{1/2}$) while maintaining useful asymptotic control. This requires the Kloosterman sum technology.

### 6.4 Subsequent Improvements

The current best result (as of this survey) is:

**Theorem (Pratt, Robles, Zaharescu, Zeindler, 2020).** *At least $41.7\%$ of the non-trivial zeros lie on the critical line.*

This uses further refinements of the mollifier method, including "two-piece" mollifiers and optimized test functions.

---

## 7. Zero-Free Regions

### 7.1 The Prime Number Theorem

The Prime Number Theorem (PNT), independently proved by Hadamard and de la Vallee-Poussin in 1896, states:

$$\pi(x) \sim \frac{x}{\log x} \quad \text{as } x \to \infty.$$

The proof requires showing $\zeta(1 + it) \neq 0$ for all $t \in \mathbb{R}$, i.e., that $\zeta$ has no zeros on the line $\sigma = 1$.

### 7.2 Classical Zero-Free Region (de la Vallee-Poussin, 1899)

**Theorem.** *There exists an absolute constant $c > 0$ such that $\zeta(\sigma + it) \neq 0$ for*

$$\sigma \geq 1 - \frac{c}{\log(|t| + 2)}.$$

*Consequently:*

$$\psi(x) = x + O\left(x \exp\left(-c'\sqrt{\log x}\right)\right).$$

**Proof idea.** The starting point is the trigonometric inequality $3 + 4\cos\theta + \cos 2\theta \geq 0$, which implies:

$$\zeta(\sigma)^3 |\zeta(\sigma + it)|^4 |\zeta(\sigma + 2it)| \geq 1$$

for $\sigma > 1$. If $\zeta$ had a zero at $1 + it_0$, then $|\zeta(\sigma + it_0)|^4$ would vanish to fourth order as $\sigma \to 1^+$, but $\zeta(\sigma)^3$ has only a third-order pole, giving a contradiction.

### 7.3 Vinogradov-Korobov Zero-Free Region (1958)

**Theorem (Vinogradov, Korobov, independently, 1958).** *There exists a constant $c > 0$ such that $\zeta(\sigma + it) \neq 0$ for*

$$\sigma \geq 1 - \frac{c}{(\log |t|)^{2/3} (\log \log |t|)^{1/3}}, \quad |t| \geq 3.$$

This is a dramatic improvement over the classical region and remains **the best known zero-free region** as of today.

**Method.** The proof uses exponential sum estimates, specifically Vinogradov's method of estimating:

$$\sum_{n=M+1}^{M+N} n^{-it} = \sum_{n} e^{-it \log n}$$

via the method of Type I and Type II sums, combined with Weyl differencing and Poisson summation. The key estimate is:

$$\zeta(\sigma + it) \ll |t|^{c(1-\sigma)^{3/2}} (\log |t|)^{2/3} \quad \text{for } 1/2 \leq \sigma \leq 1.$$

### 7.4 Consequences for the PNT

The Vinogradov-Korobov zero-free region gives the error term:

$$\psi(x) = x + O\left(x \exp\left(-c (\log x)^{3/5} (\log \log x)^{-1/5}\right)\right).$$

Compare this with the **RH-conditional** bound:

$$\psi(x) = x + O\left(x^{1/2} \log^2 x\right).$$

The gap between these two is enormous and illustrates how far we are from proving RH.

### 7.5 The Barrier at $\sigma = 1$

A fundamental obstacle: all known zero-free regions have width approaching zero as $|t| \to \infty$. No one has been able to prove a zero-free region of the form $\sigma > 1 - \delta$ for any fixed $\delta > 0$. Such a result would be a quasi-Riemann Hypothesis and would already have spectacular arithmetic consequences. The reason for this barrier is deep: our estimates for $\zeta(s)$ near the critical line rely on mean-value theorems and exponential sum bounds that inherently lose information as we approach $\sigma = 1/2$.

---

## 8. The Lindelof Hypothesis

### 8.1 Statement

**The Lindelof Hypothesis (LH).** *For every $\varepsilon > 0$,*

$$\zeta(1/2 + it) = O(|t|^{\varepsilon}).$$

The convexity bound gives $\zeta(1/2 + it) \ll |t|^{1/4}$ (from the Phragmen-Lindelof principle applied to the strip $0 \leq \sigma \leq 1$). The best known subconvexity bound is:

$$\zeta(1/2 + it) \ll |t|^{13/84 + \varepsilon}$$

(Bourgain, 2017), where $13/84 \approx 0.1548$.

### 8.2 Relationship to RH

**RH implies LH.** If RH holds, then the Lindelof Hypothesis follows from the Riemann-von Mangoldt formula and partial summation.

**LH does NOT imply RH.** The Lindelof Hypothesis is strictly weaker than RH. However, LH implies:

$$N(\sigma, T) \ll T^{\varepsilon} \quad \text{for } \sigma > 1/2,$$

where $N(\sigma, T) = \#\{\rho : \zeta(\rho)=0, \Re(\rho) \geq \sigma, 0 < \Im(\rho) \leq T\}$. That is, LH implies that for any $\sigma > 1/2$, there are at most $O(T^\varepsilon)$ zeros with real part $\geq \sigma$ up to height $T$. This is an extremely strong zero-density estimate, but it doesn't rule out any zeros off the line.

### 8.3 Equivalent Formulations

LH is equivalent to each of the following:

1. **Moment estimates:** $\int_0^T |\zeta(1/2 + it)|^{2k} \, dt \ll T^{1+\varepsilon}$ for every positive integer $k$.
2. **Mean-value estimate:** $\sum_{n \leq N} d_k(n)^2 \cdot N^{-1} \ll N^{\varepsilon}$ in an appropriate sense relating to the Dirichlet divisor function.
3. **Exponential sum bounds:** Certain exponential sums over primes have square-root cancellation.

### 8.4 Current Status

LH is unproven. The Lindelof exponent $\mu = \inf\{\alpha : \zeta(1/2+it) \ll |t|^\alpha\}$ satisfies $0 \leq \mu \leq 13/84$. RH implies $\mu = 0$ (i.e., LH). Most progress has come from bounds on exponential sums (Weyl, van der Corput, Vinogradov, Bombieri-Iwaniec, Huxley, Bourgain).

---

## 9. Li's Criterion

### 9.1 Statement

**Theorem (Li, 1997).** *Define*

$$\lambda_n = \frac{1}{(n-1)!} \frac{d^n}{ds^n} \left[s^{n-1} \log \xi(s)\right]_{s=1}.$$

*Then the Riemann Hypothesis holds if and only if $\lambda_n \geq 0$ for all $n = 1, 2, 3, \ldots$*

Equivalently, using the zeros $\rho$:

$$\lambda_n = \sum_{\rho} \left[1 - \left(1 - \frac{1}{\rho}\right)^n\right],$$

where the sum is over all non-trivial zeros $\rho$ of $\zeta(s)$.

### 9.2 Derivation and Proof

The proof of Li's criterion is elegant and proceeds via the Hadamard product.

**Forward direction (RH $\Rightarrow$ $\lambda_n \geq 0$).** If RH holds, then every non-trivial zero has the form $\rho = 1/2 + i\gamma$. Write $w = 1 - 1/\rho$. Then $|w| = 1$ (since $\rho$ maps to the unit circle under this Mobius transformation when $\Re(\rho) = 1/2$). Thus $\lambda_n = \sum_\rho (1 - \Re(w^n))$ where each term $1 - \Re(w^n) \geq 0$.

**Reverse direction ($\lambda_n \geq 0$ for all $n$ $\Rightarrow$ RH).** If some $\rho$ has $\Re(\rho) \neq 1/2$, then $|1 - 1/\rho|$ is either greater than or less than 1. The terms $\left(1 - 1/\rho\right)^n$ grow exponentially for the zeros with $|1-1/\rho| > 1$, eventually dominating the sum and making some $\lambda_n < 0$.

### 9.3 Explicit Formulas

The $\lambda_n$ can be computed more explicitly. Bombieri and Lagarias (1999) showed:

$$\lambda_n = \sum_{j=1}^{n} \binom{n}{j} \frac{(-1)^{j-1}}{(j-1)!} \left.\frac{d^j}{ds^j}\left[(s-1)\zeta(s)\right]\right|_{s=1} \cdot [\text{Gamma-related terms}].$$

More practically:

$$\lambda_n = \frac{n}{2}(\log 4\pi + \gamma) - 1 - \sum_{k=2}^{n} \binom{n}{k} \frac{(-1)^k}{(k-1)\zeta(k)} \eta_k + [\text{correction terms}],$$

where $\eta_k$ involves Stieltjes constants.

### 9.4 Computational Verification

Maslanka (2004) and others have numerically computed $\lambda_n$ for $n$ up to several thousand, finding them all positive and growing roughly like $\frac{n}{2} \log n$. The asymptotic behavior is:

$$\lambda_n \sim \frac{n}{2}(\log n + \log 4\pi + \gamma - 1) \quad \text{(if RH holds)}.$$

### 9.5 Limitations

Li's criterion transforms the Riemann Hypothesis into a positivity statement about a sequence of real numbers. While elegant, it has not led to a proof because:

1. The $\lambda_n$ involve sums over *all* zeros, creating a circular dependence.
2. Proving positivity for all $n$ requires understanding the fine distribution of zeros, which is essentially equivalent to RH.
3. Partial verification (checking finitely many $\lambda_n$) does not suffice.

---

## 10. The Nyman-Beurling-Baez-Duarte Criterion

### 10.1 The Nyman-Beurling Approach

**Theorem (Nyman, 1950; Beurling, 1955).** *The Riemann Hypothesis is true if and only if the constant function 1 can be approximated in $L^2(0,1)$ by functions of the form*

$$f(x) = \sum_{k=1}^{N} c_k \left\lfloor \frac{\theta_k}{x} \right\rfloor - \sum_{k=1}^{N} c_k \frac{\theta_k}{x}$$

*... more precisely, by functions of the form*

$$f(x) = \sum_{k=1}^{N} c_k \rho\left(\frac{\theta_k}{x}\right),$$

*where $\rho(t) = t - \lfloor t \rfloor$ is the fractional part function, $0 < \theta_k \leq 1$, and the approximation is in the $L^2(0,1)$ norm.*

Equivalently: **RH holds if and only if the subspace**

$$\mathcal{B} = \operatorname{span}\left\{ x \mapsto \rho(\theta/x) : 0 < \theta \leq 1 \right\}$$

**is dense in $L^2(0,1)$.**

### 10.2 Connection to Zeta Zeros

The Mellin transform connects this function-space characterization to the zeta function. If $f(x) = \rho(\theta/x)$, then the Mellin transform is:

$$\hat{f}(s) = \int_0^1 f(x) x^{s-1} dx = -\frac{\theta^s}{s} \cdot \frac{\zeta(s)}{???}$$

More precisely, the connection arises because the Mellin transform of $\rho(\theta/x)$ on $(0,1)$ involves $\zeta(s)$, and the density of $\mathcal{B}$ in $L^2(0,1)$ is related (via Parseval/Plancherel on the Mellin transform) to whether $1/\zeta(s)$ has certain properties on the critical line.

### 10.3 The Baez-Duarte Formulation

**Theorem (Baez-Duarte, 2003).** *RH is equivalent to the statement:*

$$\lim_{N \to \infty} d_N = 0,$$

*where*

$$d_N^2 = \inf_{A_N} \int_0^\infty \left|1 - \zeta(1/2 + it) A_N(1/2+it)\right|^2 \frac{dt}{1/4 + t^2}$$

*and $A_N(s) = \sum_{n=1}^N a_n n^{-s}$ ranges over all Dirichlet polynomials of length $N$.*

Equivalently, Baez-Duarte showed RH is equivalent to:

$$\lim_{N \to \infty} \inf_{c_0, \ldots, c_N} \left\| 1 - \sum_{k=0}^N c_k e_k \right\|_{\mathcal{H}}^2 = 0,$$

where $e_k(s) = (1 - 1/2^s)(1 - 1/3^s) \cdots$ are specific basis functions in a Hardy space $\mathcal{H}$.

### 10.4 Computational Aspects

Baez-Duarte and others have computed $d_N$ for large $N$, finding slow convergence consistent with RH. The rate of decay of $d_N$ is related to how quickly the "best approximation" to $1/\zeta(s)$ by short Dirichlet polynomials converges. The expected rate, assuming RH, is:

$$d_N^2 \sim \frac{c}{\log N}.$$

### 10.5 Burnol's Refinement

Burnol (2002) connected the Nyman-Beurling space to the co-Poisson formula and de Branges spaces, providing a deeper functional-analytic framework. This connects to the operator-theoretic approaches discussed elsewhere.

---

## 11. Weil's Explicit Formula

### 11.1 The General Explicit Formula

**Theorem (Weil, 1952).** *For test functions $f$ in a suitable class (specifically, $f$ whose Fourier transform $\hat{f}$ is smooth and rapidly decreasing), define:*

$$W(f) = \hat{f}(0) + \hat{f}(1) - \sum_{\rho} \hat{f}(\rho) = \sum_p \sum_{m=1}^{\infty} \frac{\log p}{p^{m/2}} \left[f(m \log p) + f(-m \log p)\right] + \int_{-\infty}^{\infty} f(x) \, d\alpha(x),$$

*where the sum on the left is over non-trivial zeros $\rho$ of $\zeta(s)$, the sum on the right is over primes $p$, and $\alpha(x)$ involves the Gamma function and a logarithmic term.*

More precisely, Weil's formula takes the form:

$$\sum_\rho \hat{f}(\rho) = \hat{f}(0) + \hat{f}(1) - \sum_p \sum_m \frac{\log p}{p^{m/2}}(f(m\log p) + f(-m\log p)) - \int_{-\infty}^\infty \left(\frac{e^{x/2}}{e^x - 1} - \frac{e^{-|x|/2}}{|x|}\right) f(x) \, dx.$$

### 11.2 The Positivity Criterion

**Theorem (Weil).** *The Riemann Hypothesis is equivalent to the statement that the distribution*

$$W(f) \geq 0$$

*for all test functions $f$ of the form $f = g * \tilde{g}$ where $\tilde{g}(x) = \overline{g(-x)}$ and $*$ denotes convolution.*

That is, RH $\iff$ the distribution $W$ is **positive semi-definite**.

### 11.3 Interpretation

The Weil explicit formula can be interpreted as saying that the non-trivial zeros of $\zeta(s)$ form the "spectrum" of a certain distribution on the space of test functions. The positivity condition is then a statement about this "spectrum" being contained in a specific region (the critical line).

This is deeply related to the spectral interpretation of zeros: if there exists a self-adjoint operator $H$ on a Hilbert space such that the eigenvalues of $H$ are the imaginary parts of the non-trivial zeros, then RH follows because eigenvalues of self-adjoint operators are real.

### 11.4 Connection to Number Fields

Weil's explicit formula generalizes to Dedekind zeta functions and Hecke L-functions. For function fields over finite fields, Weil proved the analogue of RH using the theory of algebraic curves (the Weil conjectures for curves, proved in 1948). The positivity of the explicit formula in that case follows from the Castelnuovo-Severi inequality in algebraic geometry.

**The key question:** Can the positivity of Weil's functional in the number field case be established by purely arithmetic-geometric methods, analogous to the function field case?

### 11.5 Bombieri's Analysis

Bombieri (2000) analyzed Weil's positivity criterion in detail and showed that for certain restricted classes of test functions, the positivity can indeed be verified. However, extending this to all admissible test functions remains out of reach.

---

## 12. The de Bruijn-Newman Constant

### 12.1 Definition

Define the family of functions parameterized by $t \in \mathbb{R}$:

$$H_t(z) = \int_0^\infty e^{tu^2} \Phi(u) \cos(zu) \, du,$$

where $\Phi$ is the super-exponentially decaying function:

$$\Phi(u) = \sum_{n=1}^{\infty}(2\pi^2 n^4 e^{9u} - 3\pi n^2 e^{5u}) \exp(-\pi n^2 e^{4u}).$$

Note that $H_0(z) = \frac{1}{8}\Xi(z/2)$ (up to normalization), so the Riemann Hypothesis is equivalent to $H_0$ having only real zeros.

**Definition.** The **de Bruijn-Newman constant** $\Lambda$ is:

$$\Lambda = \inf\{t \in \mathbb{R} : H_t \text{ has only real zeros}\}.$$

### 12.2 Basic Properties

**Theorem (de Bruijn, 1950; Newman, 1976).**

1. *$H_t$ has only real zeros for all $t \geq 1/2$ (de Bruijn).*
2. *If $H_t$ has only real zeros, then $H_{t'}$ has only real zeros for all $t' > t$ (monotonicity).*
3. *$\Lambda \leq 1/2$.*

**Theorem (Newman, 1976).** $\Lambda \geq 0$ (Newman's conjecture), under the assumption that the zeros of $\Xi$ are "well-spaced" in a certain technical sense.

Newman conjectured $\Lambda \geq 0$, noting that it would imply that the Riemann Hypothesis, if true, is "just barely" true in a precise quantitative sense.

### 12.3 Equivalence

$$\textbf{RH} \iff \Lambda \leq 0.$$

Combined with $\Lambda \geq 0$, this means:

$$\textbf{RH} \iff \Lambda = 0.$$

### 12.4 Upper Bounds on $\Lambda$

Progress on upper bounds:

| Year | Authors | Upper bound on $\Lambda$ |
|------|---------|--------------------------|
| 1950 | de Bruijn | $\leq 1/2$ |
| 1986 | Csordas, Norfolk, Varga | $\leq 0.4999$ |
| 1994 | Norfolk, Ruttan, Varga | $\leq 0.0991$ |
| 2009 | Ki, Kim, Lee | $\leq 1/2$ (different method) |
| 2019 | Polymath 15 | $\leq 0.22$ |

### 12.5 The Rodgers-Tao Theorem (2020)

**Theorem (Rodgers-Tao, 2020, building on earlier work by Rodgers, 2017).** $\Lambda \geq 0$.

This confirms Newman's conjecture and establishes that $\Lambda = 0$ if and only if RH is true.

**Proof idea (very rough sketch).** The proof analyzes the behavior of the zeros of $H_t$ as $t$ varies ("zero dynamics"). The zeros satisfy a system of ODEs (the "electrostatic" or "Dyson" flow):

$$\frac{d\gamma_k}{dt} = -2 \sum_{j \neq k} \frac{1}{\gamma_k - \gamma_j},$$

where $\gamma_k(t)$ are the (real or complex) zeros of $H_t$. As $t$ decreases from $+\infty$ toward $-\infty$, real zeros can collide, form complex pairs, and move off the real axis.

Rodgers and Tao showed that if $\Lambda < 0$, then the zeros of $H_t$ (for $t$ slightly negative) would exhibit a statistical regularity incompatible with known results about the distribution of zeros of the zeta function (specifically, pair correlation results and GUE statistics from random matrix theory).

### 12.6 Significance

The result $\Lambda = 0 \iff$ RH is philosophically important: it says the Riemann Hypothesis is a "boundary" statement. The zeros of $H_t$ are all real for $t > 0$ (no matter how small), but at $t = 0$ (the actual zeta function) we are at the exact threshold. This suggests that any proof of RH must be "tight" -- there is no room for a soft argument.

---

## 13. The Selberg Class

### 13.1 Definition

The **Selberg class** $\mathcal{S}$ consists of Dirichlet series

$$F(s) = \sum_{n=1}^{\infty} \frac{a(n)}{n^s}$$

satisfying:

1. **(Ramanujan hypothesis)** $a(n) \ll n^{\varepsilon}$ for every $\varepsilon > 0$.
2. **(Analytic continuation)** $(s-1)^m F(s)$ extends to an entire function of finite order for some non-negative integer $m$.
3. **(Functional equation)** There exist $Q > 0$, $\alpha_j > 0$, $r_j \in \mathbb{C}$ with $\Re(r_j) \geq 0$, and $|w| = 1$ such that:

$$\gamma(s) F(s) = w \overline{\gamma(1-\bar{s})} \overline{F(1-\bar{s})},$$

where $\gamma(s) = Q^s \prod_{j=1}^r \Gamma(\alpha_j s + r_j)$.

4. **(Euler product)** $\log F(s) = \sum_{n=1}^{\infty} b(n)/n^s$ where $b(n) = 0$ unless $n = p^m$ for a prime $p$ and $m \geq 1$, and $b(n) \ll n^{\theta}$ for some $\theta < 1/2$.

### 13.2 The Grand Riemann Hypothesis

**Conjecture (GRH for the Selberg Class).** *Every function $F \in \mathcal{S}$ satisfies the Riemann Hypothesis: all non-trivial zeros of $F$ have real part $1/2$.*

### 13.3 Examples

- $\zeta(s)$ (the prototype)
- Dirichlet $L$-functions $L(s, \chi)$ for primitive characters $\chi$
- Dedekind zeta functions $\zeta_K(s)$ (conjectured)
- Hecke $L$-functions associated to modular forms
- Artin $L$-functions (conjectured to be in $\mathcal{S}$)
- Rankin-Selberg $L$-functions
- Symmetric power $L$-functions

### 13.4 Degree and Classification

The **degree** of $F \in \mathcal{S}$ is $d_F = 2 \sum_{j=1}^r \alpha_j$. Key results:

- $d_F = 0$: $F = 1$ (Conrey-Ghosh).
- $d_F = 1$: $F$ is the Riemann zeta function or a shifted Dirichlet $L$-function (Kaczorowski-Perelli).
- $0 < d_F < 1$: No such $F$ exists (Richert, Conrey-Ghosh).
- $1 < d_F < 2$: No such $F$ exists (Kaczorowski-Perelli).
- $d_F = 2$: Expected to consist of $L$-functions of GL(2) automorphic forms.

### 13.5 Selberg's Conjectures

Selberg made several deep conjectures about $\mathcal{S}$:

**Conjecture A (Normality).** For $F \in \mathcal{S}$, $\sum_{p \leq x} |a(p)|^2/p = d_F \log\log x + O(1)$.

**Conjecture B (Orthogonality).** For primitive $F, G \in \mathcal{S}$ with $F \neq G$:

$$\sum_{p \leq x} \frac{a_F(p) \overline{a_G(p)}}{p} = O(1).$$

These conjectures imply unique factorization in $\mathcal{S}$ and are deeply connected to the Langlands program.

### 13.6 Why the Selberg Class Matters for RH

The Selberg class framework suggests that RH is not an isolated statement about one function, but a universal phenomenon. Understanding the structure of $\mathcal{S}$ may reveal why zeros must lie on the critical line. The functional equation and Euler product together constrain the zeros, and the Selberg class axioms are an attempt to capture exactly the properties needed for RH.

---

## 14. Moment Conjectures

### 14.1 Moments of $\zeta(s)$

Define the $2k$-th moment:

$$I_k(T) = \int_0^T |\zeta(1/2 + it)|^{2k} \, dt.$$

### 14.2 Known Results

- $k = 1$ (second moment): $I_1(T) \sim T \log T$ (Hardy-Littlewood, 1918).
- $k = 2$ (fourth moment): $I_2(T) \sim \frac{1}{2\pi^2} T (\log T)^4$ (Ingham, 1926).
- $k \geq 3$: No asymptotic formula is known.

### 14.3 Conjectures

**Conjecture (Keating-Snaith, 2000, based on random matrix theory).** For integer $k \geq 1$:

$$I_k(T) \sim C_k T (\log T)^{k^2},$$

where $C_k = \frac{g_k \cdot a_k}{(k^2)!}$ with:

$$g_k = \prod_{j=0}^{k-1} j! = (G(k+1))^2 / G(2k+1)$$

(the Barnes $G$-function), and $a_k$ is an arithmetic factor:

$$a_k = \prod_p \left(1 - \frac{1}{p}\right)^{k^2} \sum_{m=0}^{\infty} \left(\frac{\Gamma(m+k)}{m! \, \Gamma(k)}\right)^2 p^{-m}.$$

### 14.4 Connection to RH

Moments encode deep information about the value distribution of $\zeta$ on the critical line. The Lindelof Hypothesis is equivalent to:

$$I_k(T) \ll_{\varepsilon, k} T^{1+\varepsilon} \quad \text{for all } k \geq 1.$$

Understanding moments is a prerequisite for understanding the distribution of values of $\zeta(1/2+it)$, which in turn constrains where zeros can lie.

### 14.5 The CFKRS Conjecture

Conrey, Farmer, Keating, Rubinstein, and Snaith (2005) formulated a much more precise conjecture that gives all lower-order terms in the moment asymptotics. Their "recipe" produces conjectures for moments of general $L$-functions. These conjectures have been verified in many special cases and in function field analogues (by Andrade-Keating, Florea, and others).

### 14.6 Shifted Moments and Ratios Conjectures

The **Ratios Conjecture** (Conrey, Farmer, Zirnbauer) predicts the asymptotic behavior of:

$$\int_0^T \frac{\zeta(1/2 + \alpha_1 + it) \cdots \zeta(1/2 + \alpha_K + it)}{\zeta(1/2 + \beta_1 + it) \cdots \zeta(1/2 + \beta_Q + it)} \, dt$$

for shifts $\alpha_j, \beta_j$ with positive real parts. These conjectures unify many results and conjectures about the zeta function and have applications to:

- Statistics of gaps between zeros
- Moments of $|\zeta'(\rho)|$
- Mollifier optimization (relevant to Sections 5-6)
- Pair correlation and higher correlation of zeros

---

## 15. Key Barriers and Fundamental Obstacles

### 15.1 Why Each Approach Has Stalled

#### 15.1.1 The Mollifier/Proportion Approach (Sections 4-6)

**Current state:** $\geq 41.7\%$ of zeros on the critical line.

**Barrier:** The mollifier method has an inherent ceiling. The fundamental constraint is the length of mollifiers for which mean-value theorems can be proved. With mollifiers of length $T^{\theta}$:

- $\theta < 1/2$: classical mean-value theorems suffice.
- $\theta = 1/2$: requires the asymptotic large sieve and Kloosterman sum estimates.
- $\theta > 1/2$: requires estimating off-diagonal terms that involve unproved hypotheses about exponential sums.

Even if one could take $\theta$ arbitrarily close to 1, the Levinson-Conrey method would yield $N_0(T) \geq (1 - \varepsilon) N(T)$ but **never** $N_0(T) = N(T)$. The method counts zeros near the critical line but cannot pin them exactly on it. It proves "almost all zeros" but not "all zeros."

#### 15.1.2 The Zero-Free Region Approach (Section 7)

**Current state:** $\sigma \geq 1 - c/(\log |t|)^{2/3}(\log\log |t|)^{1/3}$.

**Barrier:** The Vinogradov-Korobov bound comes from exponential sum estimates. Improving it would require fundamentally new bounds for exponential sums of the form:

$$\sum_{n \sim N} e^{2\pi i f(n)},$$

where $f$ is a smooth function. The current bounds lose a power of $\log$ compared to what square-root cancellation would give. Moreover, the zero-free region approach can never prove RH, since zero-free regions of width $c/(\log t)^A$ for any $A$ still allow zeros at $\sigma = 1/2 + \varepsilon$ for any fixed $\varepsilon > 0$.

#### 15.1.3 The Lindelof Hypothesis Approach (Section 8)

**Barrier:** Even LH is strictly weaker than RH. Moreover, proving LH itself seems to require understanding zeta on the critical line at a level comparable to RH.

#### 15.1.4 Li's Criterion (Section 9)

**Barrier:** The $\lambda_n$ are sums over all zeros, so proving their positivity requires global information about the zero distribution -- essentially equivalent to RH. No partial verification can suffice. As Bombieri and Lagarias observed, the coefficients $\lambda_n$ encode the same information as the zeros themselves, so Li's criterion is a reformulation rather than a simplification.

#### 15.1.5 Nyman-Beurling (Section 10)

**Barrier:** The approximation problem in $L^2(0,1)$ is inherently infinite-dimensional. While one can compute $d_N$ numerically and observe convergence to 0, proving $d_N \to 0$ requires understanding the fine structure of $\zeta$ on the critical line. The rate of convergence is closely connected to (and essentially equivalent to) the distribution of zeros.

#### 15.1.6 Weil Positivity (Section 11)

**Barrier:** The positivity of Weil's functional for all admissible test functions is an infinite-dimensional constraint. For the function field case, positivity follows from the Castelnuovo-Severi inequality in algebraic geometry. For number fields, the analogous geometric structure is absent (or at least unknown). This is one of the deepest barriers: **there is no known geometric object whose intersection theory would yield the positivity.**

#### 15.1.7 de Bruijn-Newman (Section 12)

**Barrier:** Rodgers-Tao proved $\Lambda \geq 0$, so RH $\iff \Lambda = 0$. Proving $\Lambda \leq 0$ (i.e., $H_t$ has only real zeros at $t=0$) requires showing that no complex zeros have formed as $t$ decreases to 0. The zero dynamics (Dyson flow) are chaotic and hard to control globally. Moreover, the result $\Lambda \geq 0$ means any proof must be "tight" -- soft perturbative arguments will not work.

#### 15.1.8 The Selberg Class (Section 13)

**Barrier:** While the Selberg class provides a beautiful framework, the axioms alone do not seem to force RH. There are "fake" Dirichlet series satisfying weakened versions of the axioms that have off-line zeros. The Euler product is essential but its role is not fully understood.

### 15.2 Overarching Fundamental Obstacles

#### 15.2.1 The Arithmetic-Analytic Divide

RH sits at the intersection of analysis (properties of an analytic function) and arithmetic (distribution of primes). Most approaches are predominantly analytic and fail to fully exploit the arithmetic structure. The function field proof (Weil) crucially uses arithmetic geometry -- a number field analogue is missing.

#### 15.2.2 Absence of a Spectral Interpretation

Many approaches would be resolved if one could construct a self-adjoint operator $H$ whose spectrum is $\{\gamma : \zeta(1/2 + i\gamma) = 0\}$. The Hilbert-Polya conjecture posits such an operator exists. Despite extensive efforts (by Berry, Connes, Deninger, and others), no such operator has been rigorously constructed. The random matrix theory evidence strongly suggests that such an operator exists and belongs to a specific universality class (GUE), but this has not been made rigorous.

#### 15.2.3 The Problem of Uniqueness

RH is a statement about one specific function ($\zeta(s)$). Many proof strategies attempt to work in generality (for a class of functions satisfying certain axioms) and then specialize. But RH seems to require exploiting **specific properties** of $\zeta(s)$ -- the Euler product over rational primes, the specific functional equation, the integrality of the coefficients, the connection to the arithmetic of $\mathbb{Z}$. General arguments may be insufficient.

#### 15.2.4 The Self-Referential Nature

Many equivalences of RH (Li, Nyman-Beurling, Weil positivity) ultimately encode the same information as the zeros themselves. These reformulations are valuable for understanding, but they do not simplify the problem -- they transform it into an equally hard problem in a different language. A proof likely requires a genuinely new idea that breaks this circularity.

#### 15.2.5 The Positivity Problem

Multiple formulations reduce RH to a positivity statement (Li's $\lambda_n \geq 0$, Weil positivity, Nyman-Beurling density). Proving positivity is generally hard and often requires geometric or algebraic input. The fact that RH reduces to positivity in so many different ways suggests that there is an underlying positivity phenomenon waiting to be discovered -- perhaps connected to an inner product, a metric, or a convexity structure that we have not yet identified.

### 15.3 What Would a Proof Require?

Based on this survey, a proof of RH would likely need to:

1. **Go beyond mean-value estimates** -- current analytic methods (mollifiers, moments) asymptotically approach RH but cannot reach it.

2. **Exploit the Euler product in a fundamentally new way** -- the multiplicative structure of $\zeta$ is not fully used by any current approach.

3. **Establish a spectral or geometric interpretation** -- either construct the Hilbert-Polya operator or find a geometric framework analogous to Weil's proof for function fields.

4. **Bridge the arithmetic-analytic divide** -- combine deep number-theoretic input (arithmetic of primes, structure of $\mathbb{Z}$) with analytic techniques in a new synthesis.

5. **Prove a positivity statement** -- in one of the known equivalent formulations or a new one, establish positivity through structural (not computational) means.

---

## 16. References

### Primary Sources

1. **Riemann, B.** (1859). "Ueber die Anzahl der Primzahlen unter einer gegebenen Grosse." *Monatsberichte der Berliner Akademie.*
2. **Hadamard, J.** (1896). "Sur la distribution des zeros de la fonction $\zeta(s)$ et ses consequences arithmetiques." *Bull. Soc. Math. France*, 24, 199--220.
3. **de la Vallee-Poussin, C.J.** (1896). "Recherches analytiques sur la theorie des nombres premiers." *Ann. Soc. Sci. Bruxelles*, 20, 183--256.
4. **Hardy, G.H.** (1914). "Sur les zeros de la fonction $\zeta(s)$ de Riemann." *C. R. Acad. Sci. Paris*, 158, 1012--1014.
5. **Hardy, G.H. and Littlewood, J.E.** (1921). "The zeros of Riemann's zeta-function on the critical line." *Math. Z.*, 10, 283--317.
6. **Selberg, A.** (1942). "On the zeros of Riemann's zeta-function." *Skr. Norske Vid. Akad. Oslo I*, No. 10.
7. **Weil, A.** (1952). "Sur les 'formules explicites' de la theorie des nombres premiers." *Comm. Sem. Math. Univ. Lund [Medd. Lunds Univ. Mat. Sem.]*, Tome Supplementaire, 252--265.
8. **Levinson, N.** (1974). "More than one third of zeros of Riemann's zeta-function are on $\sigma = 1/2$." *Adv. Math.*, 13, 383--436.
9. **Conrey, J.B.** (1989). "More than two fifths of the zeros of the Riemann zeta function are on the critical line." *J. Reine Angew. Math.*, 399, 1--26.

### Zero-Free Regions and Exponential Sums

10. **Vinogradov, I.M.** (1958). "A new estimate of the function $\zeta(1+it)$." *Izv. Akad. Nauk SSSR Ser. Mat.*, 22, 161--164.
11. **Korobov, N.M.** (1958). "Estimates of trigonometric sums and their applications." *Uspehi Mat. Nauk*, 13, 185--192.
12. **Bourgain, J.** (2017). "Decoupling, exponential sums and the Riemann zeta function." *J. Amer. Math. Soc.*, 30, 205--224.

### Criteria and Equivalences

13. **Li, X.-J.** (1997). "The positivity of a sequence of numbers and the Riemann hypothesis." *J. Number Theory*, 65, 325--333.
14. **Bombieri, E. and Lagarias, J.C.** (1999). "Complements to Li's criterion for the Riemann hypothesis." *J. Number Theory*, 77, 274--287.
15. **Nyman, B.** (1950). "On the one-dimensional translation group and semi-group in certain function spaces." Thesis, Uppsala.
16. **Beurling, A.** (1955). "A closure problem related to the Riemann zeta function." *Proc. Nat. Acad. Sci.*, 41, 312--314.
17. **Baez-Duarte, L.** (2003). "A strengthening of the Nyman-Beurling criterion for the Riemann hypothesis." *Atti Accad. Naz. Lincei Cl. Sci. Fis. Mat. Natur. Rend. Lincei (9) Mat. Appl.*, 14, 5--11.

### de Bruijn-Newman Constant

18. **de Bruijn, N.G.** (1950). "The roots of trigonometric integrals." *Duke Math. J.*, 17, 197--226.
19. **Newman, C.M.** (1976). "Fourier transforms with only real zeros." *Proc. Amer. Math. Soc.*, 61, 245--251.
20. **Rodgers, B. and Tao, T.** (2020). "The de Bruijn-Newman constant is non-negative." *Forum Math. Pi*, 8, e6.

### Moments and Random Matrix Theory

21. **Keating, J.P. and Snaith, N.C.** (2000). "Random matrix theory and $\zeta(1/2+it)$." *Comm. Math. Phys.*, 214, 57--89.
22. **Conrey, J.B., Farmer, D.W., Keating, J.P., Rubinstein, M.O., and Snaith, N.C.** (2005). "Integral moments of $L$-functions." *Proc. London Math. Soc.*, 91, 33--104.

### The Selberg Class

23. **Selberg, A.** (1992). "Old and new conjectures and results about a class of Dirichlet series." *Collected Papers*, Vol. II, Springer.
24. **Kaczorowski, J. and Perelli, A.** (1999). "On the structure of the Selberg class, I." *Acta Math.*, 182, 207--241.

### Surveys and Books

25. **Titchmarsh, E.C.** (1986). *The Theory of the Riemann Zeta-Function.* 2nd ed., revised by D.R. Heath-Brown, Oxford.
26. **Iwaniec, H. and Kowalski, E.** (2004). *Analytic Number Theory.* AMS Colloquium Publications.
27. **Conrey, J.B.** (2003). "The Riemann Hypothesis." *Notices of the AMS*, 50, 341--353.
28. **Bombieri, E.** (2000). "The Riemann Hypothesis." Official Problem Description, Clay Mathematics Institute.
29. **Sarnak, P.** (2005). "Problems of the Millennium: The Riemann Hypothesis." Princeton/IAS Lectures.

---

*This survey was prepared as part of a research project on the Riemann Hypothesis. It is intended as a reference document for understanding the landscape of classical approaches, their strengths, and their limitations.*
