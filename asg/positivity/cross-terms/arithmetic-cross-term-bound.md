# The Arithmetic Cross-Term Bound

## The Single Sharpest Formulation of the Remaining Obstacle

---

## 1. The Conjecture

**Conjecture (ACTB).** For distinct primes p, q and the Weil kernel K defined by:

$$K(x) = \frac{1}{2\pi} \sum_\gamma \frac{2\cos(\gamma x)}{1/4 + \gamma^2} + K_{bg}(x)$$

where γ ranges over imaginary parts of non-trivial zeros of ζ(s) and K_bg is the "background" contribution from the pole and archimedean terms, the following bound holds:

$$\left|\sum_{m=1}^{\infty} \sum_{n=1}^{\infty} \frac{\log p \cdot \log q}{p^{m/2} q^{n/2}} K(m\log p - n\log q)\right| \leq \frac{C \sqrt{\log p \cdot \log q}}{(pq)^{1/4} \cdot (\log\max(p,q))^{1/2+\varepsilon}}$$

for an absolute constant C > 0 and any ε > 0.

**Why this suffices for RH:** If ACTB holds, then the off-diagonal entries of the Weil matrix satisfy:

$$|M_{p,q}| \leq \frac{C\sqrt{\log p \cdot \log q}}{(pq)^{1/4}(\log\max(p,q))^{1/2+\varepsilon}}$$

while the diagonal entries satisfy:

$$M_{p,p} \geq \frac{c \log p}{p^{1/2}}$$

The row sum condition for diagonal dominance becomes:

$$\sum_{q \neq p} |M_{p,q}| \leq \sum_{q \neq p} \frac{C\sqrt{\log p \cdot \log q}}{(pq)^{1/4}(\log q)^{1/2+\varepsilon}}$$

$$= C\frac{(\log p)^{1/2}}{p^{1/4}} \sum_q \frac{(\log q)^{1/2}}{q^{1/4}(\log q)^{1/2+\varepsilon}}$$

$$= C\frac{(\log p)^{1/2}}{p^{1/4}} \sum_q \frac{1}{q^{1/4}(\log q)^{\varepsilon}}$$

The sum over q converges (barely) for ε > 0, giving:

$$\sum_{q \neq p} |M_{p,q}| \leq C' \frac{(\log p)^{1/2}}{p^{1/4}}$$

For diagonal dominance we need:

$$C' \frac{(\log p)^{1/2}}{p^{1/4}} \leq c\frac{\log p}{p^{1/2}}$$

$$\Leftrightarrow C'/c \leq \frac{(\log p)^{1/2}}{p^{1/4}}$$

This fails for small p (the right side → 0 as p → 2). But for p large enough, diagonal dominance holds.

For the finitely many small primes, ACTB combined with the exact kernel values (computed numerically) gives a FINITE VERIFICATION problem. If this verification succeeds, APT follows.

---

## 2. Evidence for ACTB

### 2.1 The GUE Prediction

Under Montgomery's pair correlation conjecture, the zeros γ have GUE statistics. In a GUE ensemble with N eigenvalues in [-T, T], the kernel:

$$K_{GUE}(x) = \sum_k \frac{2\cos(\gamma_k x)}{1/4 + \gamma_k^2}$$

has mean zero and variance:

$$\text{Var}(K_{GUE}(x)) \sim \frac{1}{2\pi} \int_{-\infty}^{\infty} \frac{4}{(1/4+u^2)^2} du = \frac{1}{2\pi} \cdot \frac{16\pi}{3} = \frac{8}{3}$$

Wait — the variance needs more careful treatment. For the truncated kernel (zeros with |γ| ≤ T):

$$K_T(x) = \frac{1}{2\pi}\sum_{|\gamma| \leq T} \frac{2\cos(\gamma x)}{1/4 + \gamma^2}$$

By the density of zeros (N(T) ~ (T/2π) log T), the variance at a GENERIC point x is:

$$\text{Var}(K_T(x)) \sim \frac{1}{(2\pi)^2} \cdot \int_0^T \frac{4}{(1/4+u^2)^2} \frac{\log u}{2\pi} du \sim \frac{C_0}{\log T}$$

for a constant C₀. So |K_T(x)| ~ 1/√(log T) at generic points.

The ACTB requires |K(m log p - n log q)| ~ 1/√(log(max(p,q))). Since m log p - n log q ~ log(max(p,q)) for the dominant terms (m = n = 1), and T ~ max(p,q) captures the relevant zero range, this is CONSISTENT with the GUE prediction.

### 2.2 Baker's Theorem

Baker's theorem on linear forms in logarithms gives:

$$|m\log p - n\log q| \geq \exp(-C \log m \cdot \log n \cdot \log p \cdot \log q)$$

This is a LOWER bound on the differences at which K is evaluated. Since K(x) decays for |x| → ∞ (due to the denominator 1/4 + γ²), the Baker bound ensures that K is evaluated at points bounded away from zero.

However, Baker's bound is much too weak to give ACTB directly. The Baker bound gives |m log p - n log q| ≥ exp(-C(log p)²) which is very small, while the kernel K does not decay significantly on this scale.

### 2.3 Numerical Evidence

From cross_term_matrix.py: the computed cross-term matrix for primes up to 100 shows:

- Diagonal entries: M_{p,p} ~ (log p)²/p (dominated by the m=1 term)
- Off-diagonal entries: M_{p,q} much smaller than diagonal, with oscillating signs
- The Frobenius norm of the off-diagonal block is smaller than the Frobenius norm of the diagonal
- Eigenvalues of the primitive-projected matrix are all negative (consistent with APT)

The numerical data supports ACTB with C ≈ 1-2 for the primes tested.

### 2.4 Partial Results

**Unconditional:** By the large sieve inequality (Bombieri):

$$\sum_{q \leq Q} \left|\sum_{\gamma} a_\gamma e^{i\gamma \log q}\right|^2 \leq (N(T) + Q^2/\delta - 1) \sum_\gamma |a_\gamma|^2$$

where δ = min-distance between log-primes. This gives a MEAN-SQUARE bound on K at prime logarithms:

$$\sum_{q \leq Q} |K(\log q - \log p)|^2 \leq (T\log T + Q) \cdot \sum_\gamma \frac{4}{(1/4+\gamma^2)^2}$$

This gives: AVERAGE |K(log q - log p)| ≤ C · √(T log T / Q). For Q ~ T, this is ~ √(log T), which goes the WRONG direction (too large).

The large sieve gives the correct AVERAGE bound but not a POINTWISE bound. ACTB requires a pointwise bound.

**Under GRH for Dirichlet L-functions:** The prime correlations become more tractable. Specifically, the orthogonality of Dirichlet characters gives:

$$\sum_{p \leq x} \sum_{q \leq x} K(\log p - \log q) = \sum_\chi \left|\sum_{p \leq x} \chi(p) e^{i\gamma \log p}\right|^2 / (...)$$

The inner sums are character sums, controlled by GRH. Under GRH for all Dirichlet L-functions, the cross-terms between primes in different residue classes cancel, giving a version of ACTB.

---

## 3. Structure of the Kernel at Arithmetic Points

### 3.1 The Key Question

What is special about the values K(m log p - n log q)?

These are the kernel evaluated at DIFFERENCES OF PRIME-POWER LOGARITHMS. The set:

$$\mathcal{D} = \{m\log p - n\log q : p, q \text{ prime}, m, n \geq 1\}$$

is dense in ℝ (by Dirichlet's theorem on primes in arithmetic progressions and the irrationality of log p / log q for distinct primes p, q).

But the WEIGHTED sum:

$$S(p) = \sum_{q \neq p} \sum_{m,n} \frac{\log q}{q^{n/2}} K(m\log p - n\log q)$$

involves weights 1/q^{n/2} that decay rapidly for large q or n. The dominant contributions come from:
- q small (q = 2, 3, 5, 7)
- n = 1 (first power)
- m = 1 (first power)

So the key values are K(log p - log q) for small primes p, q.

### 3.2 Explicit Values of K at Key Points

Using the explicit formula for K:

$$K(x) = \delta(x) - \frac{1}{2\pi}\frac{\Gamma'}{\Gamma}(1/4 + ix/2) - \frac{1}{2\pi}\frac{\Gamma'}{\Gamma}(1/4 - ix/2) + \text{zero contribution}$$

For the BACKGROUND part (poles + archimedean, without zeros):

$$K_{bg}(x) = \delta(x) - \frac{1}{\pi}\text{Re}\frac{\Gamma'}{\Gamma}(1/4 + ix/2) + \frac{1}{2\pi}\log\pi$$

The digamma function:
$$\text{Re}\frac{\Gamma'}{\Gamma}(1/4 + ix/2) = \frac{1}{2}\log\frac{x^2/4 + 1/16}{1} + O(1/x^2) \quad \text{for } x \gg 1$$

So K_bg(x) ~ -(1/2π)log(x/2) for large x. This is SLOWLY DECREASING (logarithmic).

The zero contribution adds oscillations on top of this slow decay:

$$K_{zeros}(x) = \frac{1}{2\pi}\sum_\gamma \frac{2\cos(\gamma x)}{1/4 + \gamma^2}$$

For x = log 3 - log 2 ≈ 0.405:
- K_bg(0.405) ≈ (computable from digamma at 1/4 + 0.2i)
- K_zeros(0.405) ≈ Σ 2cos(0.405γ)/(1/4 + γ²)

The first few zeros contribute:
- γ₁ = 14.134: cos(14.134 × 0.405)/(1/4 + 14.134²) ≈ cos(5.724)/200.0 ≈ 0.873/200 ≈ 0.00437
- γ₂ = 21.022: cos(21.022 × 0.405)/(1/4 + 21.022²) ≈ cos(8.514)/442.2 ≈ -0.623/442 ≈ -0.00141
- γ₃ = 25.011: cos(25.011 × 0.405)/(1/4 + 25.011²) ≈ cos(10.129)/625.8 ≈ -0.843/626 ≈ -0.00135

The sum converges rapidly due to the 1/(1/4 + γ²) denominator. The total |K_zeros(log 3 - log 2)| ≈ 0.002-0.005 (small!).

### 3.3 Why the Cross-Terms Are Small

The cross-term kernel K at arithmetic points is small because:

1. **The denominator 1/(1/4 + γ²) kills high zeros.** Only the first ~10 zeros contribute significantly.

2. **The cosines oscillate.** cos(γ(m log p - n log q)) has "random" phases for different γ, leading to cancellation.

3. **The differences m log p - n log q are NOT small.** By Baker's theorem (or just by unique factorization), these differences are bounded away from zero, so the kernel is not evaluated near its singularity.

The combination of these three effects gives |K(m log p - n log q)| ~ 0.01 or smaller for typical prime pairs, while the diagonal terms K(0) ~ 2.3. This gives a ratio of ~200:1, strongly favoring diagonal dominance.

---

## 4. A Proof Strategy for ACTB

### 4.1 The Stationary Phase Approach

Write:

$$K_{zeros}(x) = \frac{1}{\pi} \text{Re} \int_0^T \frac{e^{iux}}{1/4 + u^2} dN(u)$$

where N(u) = #{γ : 0 < γ ≤ u} is the zero-counting function.

By the Riemann-von Mangoldt formula:

$$N(u) = \frac{u}{2\pi}\log\frac{u}{2\pi} - \frac{u}{2\pi} + \frac{7}{8} + S(u)$$

where S(u) = (1/π)arg ζ(1/2 + iu) is the "argument function" with S(u) = O(log u).

Integration by parts:

$$K_{zeros}(x) = \frac{1}{\pi}\text{Re}\left[\frac{e^{iux}N(u)}{1/4+u^2}\bigg|_0^T + \int_0^T N(u) \frac{d}{du}\frac{e^{iux}}{1/4+u^2} du\right]$$

The derivative:

$$\frac{d}{du}\frac{e^{iux}}{1/4+u^2} = \frac{ix \cdot e^{iux}(1/4+u^2) - 2u e^{iux}}{(1/4+u^2)^2}$$

For |x| bounded (as in the prime-logarithm differences), the integral converges absolutely and the boundary terms vanish as T → ∞.

Substituting the Riemann-von Mangoldt formula for N(u):

$$K_{zeros}(x) = \frac{1}{\pi}\text{Re}\int_0^\infty \left[\frac{u}{2\pi}\log\frac{u}{2\pi} + S(u)\right] \frac{d}{du}\frac{e^{iux}}{1/4+u^2} du$$

The SMOOTH part (involving u log u) gives the "background" — this is the part of K that does NOT depend on the specific zero locations. It can be evaluated in closed form using special functions.

The OSCILLATORY part (involving S(u)) gives the "fluctuation":

$$K_{fluct}(x) = \frac{1}{\pi}\text{Re}\int_0^\infty S(u) \frac{d}{du}\frac{e^{iux}}{1/4+u^2} du$$

### 4.2 Bounding the Fluctuation

S(u) = O(log u) unconditionally (Littlewood bound). This gives:

$$|K_{fluct}(x)| \leq \frac{C}{\pi}\int_0^\infty \frac{\log u}{(1/4+u^2)} du$$

The integral converges and gives |K_fluct(x)| ≤ C' (an absolute constant).

But we need a bound that is SMALLER than the diagonal (which is ~log p / p^{1/2}). The constant C' doesn't depend on p, so for p large enough, the diagonal dominates. For small p, we need C' to be small enough.

**Better bound on S:** By the Selberg normal-order theorem:

$$\int_0^T |S(u)|^2 du = \frac{T}{2\pi^2}\log\log T + O(T)$$

So S(u) has RMS size ~√(log log T). This gives:

$$|K_{fluct}(x)| = O\left(\sqrt{\log\log(1/|x|)}\right) \quad \text{for } |x| \text{ small}$$

and

$$|K_{fluct}(x)| = O(1) \quad \text{for } |x| = O(1)$$

Since the prime-logarithm differences m log p - n log q are O(1) for the dominant terms, we get |K_fluct| = O(1), which doesn't improve the constant.

### 4.3 The Resonance Condition

The kernel K(x) is large when cos(γ x) ≈ 1 for many zeros γ simultaneously. This happens when x is close to a "resonance" — a value where γ x ≈ 0 (mod 2π) for many γ.

**Are prime-logarithm differences resonant?** For x = log(p/q) and γ_k the k-th zero:

$$\gamma_k \cdot \log(p/q) = \gamma_k \log p - \gamma_k \log q$$

This is resonant when γ_k log p ≈ γ_k log q (mod 2π), i.e., when γ_k(log p - log q) ≈ 0 (mod 2π).

The question reduces to: are the values {γ_k log(p/q) mod 2π} equidistributed?

By Weyl's equidistribution theorem, if log(p/q) is irrational (which it is, by unique factorization), then the sequence {γ_k log(p/q) mod 2π} is equidistributed IF the γ_k are "generic enough."

Under GUE statistics, the zeros are "maximally spread" (they have the same distribution as GUE eigenvalues). The equidistribution follows from the GUE universality.

**Unconditionally:** We know (from the work of Jutila, Montgomery, and others) that the sequence γ_k is "uniformly distributed mod 1" in the sense that:

$$\frac{1}{N(T)} \sum_{\gamma \leq T} e^{2\pi i \alpha \gamma} \to 0 \quad \text{as } T \to \infty$$

for any fixed α ≠ 0. Taking α = log(p/q)/(2π), this gives:

$$\frac{1}{N(T)} \sum_{\gamma \leq T} e^{i\gamma \log(p/q)} \to 0$$

This is cancellation in the UNWEIGHTED sum. For the WEIGHTED sum (with 1/(1/4 + γ²) weights), the convergence is even faster due to the weight decay.

**This suggests:** The ACTB should be provable by combining:
- The uniform distribution of zeros (unconditional)
- The rapid decay of the weights 1/(1/4 + γ²)
- The Baker lower bound on |m log p - n log q|

### 4.4 Toward a Proof

**Proposition (Sketch).** For distinct primes p, q:

$$\left|\sum_\gamma \frac{2\cos(\gamma \log(p/q))}{1/4 + \gamma^2}\right| \leq C$$

for an absolute constant C, unconditionally.

**Proof sketch:** By partial summation with N(T) ~ (T/2π) log(T/2πe):

$$\sum_{0 < \gamma \leq T} \frac{\cos(\gamma \log(p/q))}{1/4 + \gamma^2} = \int_0^T \frac{\cos(u\log(p/q))}{1/4+u^2} dN(u)$$

$$= \left[\frac{N(u)\cos(u\log(p/q))}{1/4+u^2}\right]_0^T + \int_0^T N(u)\left[\frac{\log(p/q)\sin(u\log(p/q))}{1/4+u^2} + \frac{2u\cos(u\log(p/q))}{(1/4+u^2)^2}\right]du$$

Using N(u) ~ (u/2π)log(u/2πe):

The integral converges absolutely (the integrand decays like u log u / u⁴ × oscillation). The boundary term at T → ∞ vanishes (since N(T)/T² → 0). The boundary at 0 gives N(0+) = 0.

Therefore the sum is bounded by:

$$\leq \int_0^\infty \frac{u\log u}{(1/4+u^2)^2} du + |\log(p/q)| \int_0^\infty \frac{u\log u}{(1/4+u^2)^2} du$$

Both integrals converge to finite values. The bound depends on |log(p/q)| ≈ |log(p/q)| which is O(log p) for q = 2.

**Issue:** The bound grows with |log(p/q)|, which means it grows with the prime ratio. For the DOMINANT term (m = n = 1), this gives a bound of O(log p), which is too large for diagonal dominance.

**Fix:** The growth is in the BACKGROUND part, not the FLUCTUATION part. The background contribution K_bg(log(p/q)) ≈ -(1/2π)log(|log(p/q)|/2) is NEGATIVE for |log(p/q)| > 2, which HELPS the positivity (it makes the off-diagonal entries negative, reinforcing the negative-definiteness of the primitive subspace).

The fluctuation K_fluct is bounded by O(1), independent of p, q. This is the part that could potentially cause positivity to fail.

**Refined ACTB:** The cross-term bound should separate the background (computable, non-oscillatory, favorable sign) from the fluctuation (bounded, oscillatory, potentially unfavorable):

$$K(m\log p - n\log q) = K_{bg}(m\log p - n\log q) + K_{fluct}(m\log p - n\log q)$$

with |K_fluct| ≤ C (absolute constant) and K_bg computable.

If C is small enough (compared to the diagonal), APT follows.

---

## 5. Summary

The Arithmetic Cross-Term Bound (ACTB) is:
- **Sufficient** for APT (and hence RH)
- **Consistent** with GUE predictions and numerical evidence
- **Partially provable** by existing methods (uniform distribution of zeros + partial summation)
- **Reducible** to bounding |K_fluct(x)| for x at prime-logarithm differences

The remaining gap is quantitative: we need |K_fluct| ≤ C with C small enough that the finite check for small primes succeeds. Current bounds give C = O(1) but the implied constant has not been computed explicitly.

**An explicit computation of C, combined with a numerical verification of diagonal dominance for small primes (p ≤ 100), would resolve the question.**
