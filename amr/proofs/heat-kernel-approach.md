# Heat Kernel on the Adelic Solenoid and the Weil Kernel

## Status: Exploratory. Identifies a precise structural parallel between K_bg and the Green's function of an archimedean Laplacian on T_A, but the full K_Weil requires arithmetic (multi-place) contributions that obstruct the naive heat-kernel representation. The Connes prolate-wave operator provides a more promising spectral framework.

---

## 0. Motivation

**Core question:** Can the Weil kernel K be represented as a positive mixture of heat kernels on the adelic solenoid T_A?

If K = ∫₀^∞ f(t) h_t dt for some f(t) ≥ 0 and heat kernel h_t on T_A, then K is automatically positive definite (PD), since each h_t is PD (being a transition density) and positive mixtures preserve PD. This would prove CPD-1 and hence RH.

**Weaker question:** Is K_bg (the background kernel) the Green's function of a natural Laplacian on T_A? This would give an independent proof that K_bg is CPD-1 (already known from Fourier analysis, bochner-proof.md §1.4) and illuminate the geometric meaning of K_bg.

---

## 1. The Adelic Solenoid

### 1.1 Definition and Structure

The adelic solenoid is:

$$\mathbb{T}_\mathbb{A} = \mathbb{R}/\mathbb{Z} \times \prod_p \mathbb{Z}_p$$

where the product is over all primes p, equipped with the product topology. This is a compact, metrizable, connected abelian group.

**Equivalent descriptions:**
- Inverse limit: T_A = lim_{←n} ℝ/nℤ (over positive integers ordered by divisibility)
- Adelic quotient: T_A ≅ A_Q / Q̂ where A_Q is the adele ring and Q̂ is the rational closure
- Pontryagin dual: T̂_A ≅ (ℚ, discrete), the rationals with discrete topology

### 1.2 Characters

Every continuous character χ: T_A → S¹ has the form:

$$\chi_r(x) = e^{2\pi i r \cdot \pi_\infty(x)} \qquad (r \in \mathbb{Q})$$

where π_∞: T_A → ℝ/ℤ is the archimedean projection. The characters are indexed by ℚ.

For r = a/b in lowest terms, χ_r factors through the finite quotient ℝ/bℤ and is well-defined on T_A because a/b · n ∈ ℤ for all n ∈ bℤ.

### 1.3 Haar Measure

The Haar probability measure λ on T_A is the product:

$$\lambda = \lambda_\infty \times \prod_p \mu_p$$

where λ_∞ is Lebesgue measure on ℝ/ℤ and μ_p is the Haar probability measure on ℤ_p.

### 1.4 Fourier Analysis

For f ∈ L²(T_A, λ), the Fourier expansion is:

$$f(x) = \sum_{r \in \mathbb{Q}} \hat{f}(r) \chi_r(x), \qquad \hat{f}(r) = \int_{\mathbb{T}_\mathbb{A}} f(x) \overline{\chi_r(x)}\, d\lambda(x)$$

The Plancherel theorem: ||f||² = Σ_{r ∈ ℚ} |f̂(r)|².

---

## 2. The Laplacian on T_A

### 2.1 Spectral Definition

On any compact abelian group G with dual Ĝ, a Laplacian Δ is specified by choosing eigenvalues: for each character χ ∈ Ĝ, set

$$\Delta \chi = -\lambda(\chi) \chi, \qquad \lambda(\chi) \geq 0, \quad \lambda(\chi) = 0 \iff \chi = 1$$

This defines a non-positive self-adjoint operator on L²(G) whose kernel is the constants.

For T_A with Ĝ = ℚ, a Laplacian is determined by a function λ: ℚ → [0,∞) with λ(0) = 0.

### 2.2 Natural Choices of Eigenvalues

**Choice A: Archimedean Laplacian.**

$$\lambda_\infty(r) = 4\pi^2 r^2$$

This is the pull-back of the standard Laplacian on ℝ/ℤ via the archimedean projection. The eigenvalues depend only on the archimedean absolute value |r|.

**Choice B: Full Adelic Laplacian.**

$$\lambda_{\text{ad}}(r) = 4\pi^2 r^2 + \sum_p (\log p)^2 \cdot v_p(r)_-^2$$

where v_p(r)_- = max(0, -v_p(r)) counts the "denominator valuation" at p. This incorporates contributions from all places: a rational r = a/b with p | b (i.e., v_p(r) < 0) incurs a p-adic penalty.

**Choice C: Connes' Arithmetic Operator.** (See §5 below.)

$$\lambda_D(r) = \left(\sum_p v_p(r) \log p\right)^2 + r^2$$

inspired by Connes' operator D = Σ_p (log p) D_p on adelic spaces.

### 2.3 The Heat Kernel

For any choice of eigenvalues λ(r), the heat kernel on T_A is:

$$h_t(x) = \sum_{r \in \mathbb{Q}} e^{-t\lambda(r)} \chi_r(x)$$

**Properties (standard for compact groups):**
1. h_t(x) > 0 for all t > 0, x ∈ T_A (by Weierstrass/compactness)
2. h_t is positive definite for each t > 0 (Fourier coefficients e^{-tλ(r)} ≥ 0)
3. h_t * h_s = h_{t+s} (semigroup property)
4. h_t → δ_0 as t → 0⁺ (approximation to identity)
5. ∫_{T_A} h_t dλ = 1 for all t > 0

---

## 3. Green's Function and the Digamma Connection

### 3.1 The Green's Function

The Green's function (resolvent kernel at spectral parameter 0) is formally:

$$G(x) = \int_0^\infty (h_t(x) - 1)\, dt = \sum_{r \in \mathbb{Q} \setminus \{0\}} \frac{\chi_r(x)}{\lambda(r)}$$

(The subtraction of 1 removes the zero mode; the integral converges iff Σ_{r≠0} 1/λ(r) < ∞.)

G satisfies: -ΔG = δ_0 - 1 (distributional identity on T_A).

### 3.2 Archimedean Green's Function

For the archimedean Laplacian (Choice A), λ_∞(r) = 4π²r²:

$$G_\infty(x) = \sum_{r \in \mathbb{Q} \setminus \{0\}} \frac{e^{2\pi i r \pi_\infty(x)}}{4\pi^2 r^2}$$

**Restricting to integers** (r ∈ ℤ \ {0}):

$$G_\infty^{(\mathbb{Z})}(\theta) = \sum_{n=1}^{\infty} \frac{2\cos(2\pi n \theta)}{4\pi^2 n^2} = \frac{1}{2\pi^2} \sum_{n=1}^{\infty} \frac{\cos(2\pi n \theta)}{n^2}$$

This equals (1/(2π²)) · [π²B_2(θ)] = B_2(θ)/2 where B_2(θ) = θ² - θ + 1/6 is the second Bernoulli polynomial (for θ ∈ [0,1]).

**Including all rationals** requires summing over ℚ \ {0}, which converges conditionally. Grouping by denominator b:

$$G_\infty(\theta) = \sum_{b=1}^{\infty} \sum_{\substack{a=1 \\ \gcd(a,b)=1}}^{b} \frac{e^{2\pi i (a/b) \theta}}{4\pi^2 (a/b)^2}$$

The inner sum involves Ramanujan sums and Euler's totient function.

### 3.3 Comparison with K_bg

Recall from bochner-proof.md §1.1:

$$K_{\text{bg}}(x) = -\frac{1}{\pi}\text{Re}[\psi(1/4 + ix/2)] + \frac{\log\pi}{2\pi}$$

Using the digamma series expansion:

$$K_{\text{bg}}(x) = \frac{1}{\pi}\sum_{n=0}^{\infty} \frac{n + 1/4}{(n+1/4)^2 + x^2/4} + (\text{constants contributing at } x = 0)$$

Each summand is a Lorentzian: (n + 1/4)/[(n + 1/4)² + x²/4].

**Key observation:** The Fourier transform of K_bg is (bochner-proof.md §1.3):

$$\hat{K}_{\text{bg}}(\xi) = \frac{2e^{-|\xi|/2}}{1 - e^{-2|\xi|}} = \frac{1}{\sinh(|\xi|)}$$

This can be written as:

$$\hat{K}_{\text{bg}}(\xi) = 2\sum_{n=0}^{\infty} e^{-(2n+1/2)|\xi|}$$

**Compare with a Green's function:** For a Laplacian with eigenvalues λ_n, the Green's function has Fourier coefficients 1/λ_n. If we set:

$$\frac{1}{\lambda_n} = 2 e^{-(2n+1/2)|\xi|}$$

This does NOT correspond to a standard eigenvalue problem, because the "eigenvalue" would depend on the Fourier variable ξ. The Green's function interpretation works only in the sense that:

$$K_{\text{bg}}(x) = \int_0^\infty w(t) h_t^\infty(x)\, dt$$

where h_t^∞ is the archimedean heat kernel on ℝ/ℤ restricted to T_A, and w(t) is an appropriate weight.

### 3.4 Heat Kernel Representation of K_bg

**Theorem 3.1.** *K_bg admits the representation:*

$$K_{\text{bg}}(x) = \int_0^\infty w(t) \cdot p_t(x)\, dt$$

*where p_t(x) = Σ_{n∈ℤ} (4πt)^{-1/2} exp(-(x-n)²/(4t)) is the heat kernel on ℝ (periodized to ℝ/ℤ, evaluated at the archimedean projection), and w(t) is a function determined by the relation:*

$$\hat{K}_{\text{bg}}(\xi) = \int_0^\infty w(t) e^{-4\pi^2 \xi^2 t}\, dt$$

*i.e., K̂_bg is the Laplace transform of w evaluated at 4π²ξ².*

**Proof sketch.** The Fourier coefficients of p_t at frequency n are e^{-4π²n²t}. Setting K̂_bg(n) = ∫₀^∞ w(t) e^{-4π²n²t} dt, we need:

$$\frac{2e^{-|n|/2}}{1 - e^{-2|n|}} = \int_0^\infty w(t) e^{-4\pi^2 n^2 t}\, dt$$

This is a Laplace inversion problem. Since the LHS is completely monotone in n² (as a function of n² > 0), the Bernstein-Widder theorem guarantees that w(t) ≥ 0 exists.

**Explicit computation:** Setting u = 4π²n²t, we need:

$$\frac{2e^{-|n|/2}}{1 - e^{-2|n|}} = \frac{1}{4\pi^2 n^2}\int_0^\infty w\left(\frac{u}{4\pi^2 n^2}\right) e^{-u}\, du$$

This is more delicate. The function K̂_bg(ξ) = 1/sinh(|ξ|) satisfies:

$$\frac{1}{\sinh s} = 2\int_0^\infty e^{-s\cosh\theta}\, d\theta \qquad (s > 0)$$

(This is a standard integral representation involving the modified Bessel function K_0.)

So:

$$\hat{K}_{\text{bg}}(\xi) = 2\int_0^\infty e^{-|\xi|\cosh\theta}\, d\theta$$

This shows K̂_bg is a Laplace transform of a POSITIVE measure in the variable |ξ|, but not in ξ². The heat kernel representation requires the Laplace transform in ξ², so an additional change of variables is needed, and the result may not preserve positivity of w in the t-variable directly.

**Result:** K_bg admits a heat kernel representation with a SIGNED weight w(t) (not necessarily non-negative), which means the naive route "K_bg = positive mixture of PD functions ⟹ K_bg is PD" does NOT directly apply from the heat kernel perspective. However, K_bg is still CPD-1 by the direct Fourier argument (bochner-proof.md §1.4). ∎

### 3.5 The Digamma as a Regularized Green's Function

A more illuminating connection: the digamma ψ(s) satisfies:

$$\psi(s) = -\gamma + \int_0^\infty \left(\frac{e^{-t}}{t} - \frac{e^{-st}}{1 - e^{-t}}\right) dt \qquad (\text{Re}(s) > 0)$$

Evaluating at s = 1/4 + ix/2:

$$\text{Re}[\psi(1/4 + ix/2)] = -\gamma + \int_0^\infty \left(\frac{e^{-t}}{t} - \frac{e^{-t/4}\cos(xt/2)}{1 - e^{-t}}\right) dt$$

So:

$$K_{\text{bg}}(x) = \frac{1}{\pi}\int_0^\infty \frac{e^{-t/4}\cos(xt/2)}{1 - e^{-t}}\, dt + (\text{constants})$$

The integrand e^{-t/4} cos(xt/2)/(1 - e^{-t}) can be interpreted as: at each "time" t, the function cos(xt/2) is a character evaluated at x (like a heat kernel mode), weighted by the positive measure e^{-t/4}/(1-e^{-t}) dt.

**This IS a positive mixture of cosine modes**, which is essentially Bochner's theorem applied directly. It is the real-space analogue of the Fourier positivity K̂_bg > 0.

---

## 4. The p-Adic Heat Kernel

### 4.1 The Vladimirov-Taibleson Operator

On ℚ_p, the Vladimirov operator of order α > 0 is:

$$D_p^\alpha f(x) = \frac{1 - p^\alpha}{1 - p^{-\alpha-1}} \int_{\mathbb{Q}_p} \frac{f(x) - f(y)}{|x - y|_p^{\alpha + 1}}\, dy$$

This is the p-adic analogue of the fractional Laplacian (-Δ)^{α/2} on ℝ.

### 4.2 Heat Kernel on ℤ_p

Restricting to ℤ_p (the p-adic integers), the natural Laplacian is the restriction of D_p^α. The characters of ℤ_p are indexed by ℚ_p/ℤ_p ≅ ℤ[1/p]/ℤ, represented as:

$$\chi_{a/p^k}(x) = e^{2\pi i \{ax/p^k\}} \qquad (a \in \{1, \ldots, p^k-1\}, \gcd(a,p) = 1)$$

The eigenvalues of D_p^α on these characters are:

$$\lambda_p(a/p^k) = p^{k\alpha}$$

(depending only on the "level" k = -v_p(a/p^k)).

The heat kernel on ℤ_p is:

$$h_t^{(p)}(x) = 1 + \sum_{k=1}^{\infty} e^{-t p^{k\alpha}} \sum_{\substack{a \bmod p^k \\ \gcd(a,p)=1}} \chi_{a/p^k}(x)$$

The inner sum equals the Ramanujan-type sum:

$$\sum_{\substack{a \bmod p^k \\ \gcd(a,p)=1}} e^{2\pi i a x / p^k} = \begin{cases} p^k - p^{k-1} & \text{if } x \in p^k \mathbb{Z}_p \\ -p^{k-1} & \text{if } x \in p^{k-1}\mathbb{Z}_p \setminus p^k\mathbb{Z}_p \\ 0 & \text{otherwise} \end{cases}$$

So the p-adic heat kernel is piecewise constant on cosets of p^k ℤ_p:

$$h_t^{(p)}(x) = \begin{cases} 1 + (p^k - p^{k-1}) e^{-tp^{k\alpha}} + \cdots & \text{if } x \in p^k\mathbb{Z}_p \\ \vdots & \end{cases}$$

**Key property:** h_t^{(p)}(x) is positive definite on ℤ_p for each t > 0 (being a heat kernel on a compact group).

### 4.3 Product Structure

On the full solenoid T_A = ℝ/ℤ × Π_p ℤ_p, if the Laplacian decomposes as:

$$\Delta = \Delta_\infty + \sum_p \Delta_p$$

(sum of independent Laplacians on each factor), then the heat kernel factorizes:

$$h_t(x) = h_t^\infty(\pi_\infty(x)) \cdot \prod_p h_t^{(p)}(\pi_p(x))$$

This is the **Euler product structure** of the solenoid heat kernel: the diffusion is independent on each local factor.

**The infinite product converges** because for almost all primes p and any fixed t > 0, the p-adic heat kernel satisfies h_t^{(p)} ≈ 1 + O(e^{-tp^α}) → 1 as p → ∞.

---

## 5. Connection to Connes' Arithmetic Operator

### 5.1 Connes' Spectral Realization

In Connes' noncommutative geometry approach (arXiv:math/9811068), the nontrivial zeros of ζ(s) arise as an absorption spectrum from an operator on the adele class space A_Q/Q*.

The key objects:
- The scaling action: σ_λ for λ ∈ C_Q = A_Q*/Q* (idele class group)
- The trace formula: Σ_ρ h(ρ) = (explicit formula terms)
- Weil positivity: W(f * f̃) = Σ_ρ |f̂(ρ)|² ≥ 0 under RH

### 5.2 The Prolate Wave Operator

Connes-Consani-Moscovici (arXiv:2310.18423, 2024) introduced the prolate wave operator W_Λ whose spectrum connects to zeta zeros:

- **Positive eigenvalues** of W_Λ realize the low-lying zeta zeros
- **Negative eigenvalues** correspond to the Sonin space, where Weil positivity at the archimedean place is proved
- The operator W_Λ = S² + Λ²T², where S is the scaling operator and T is related to the Fourier-type projection

### 5.3 Relation to Heat Kernels

The heat kernel e^{-tW_Λ} of the prolate operator would give a semigroup whose trace is:

$$\text{Tr}(e^{-tW_\Lambda}) = \sum_n e^{-t\mu_n}$$

where μ_n are eigenvalues of W_Λ. The zeta function of W_Λ is:

$$\zeta_{W_\Lambda}(s) = \text{Tr}(W_\Lambda^{-s}) = \sum_n \mu_n^{-s}$$

**Connection to the Riemann zeta:** The eigenvalues μ_n of the prolate operator reproduce (in the UV limit) the squares of the imaginary parts γ² of zeta zeros. So the "arithmetic heat kernel" e^{-tW_Λ} encodes, in its spectral decomposition, the same information as the zero-oscillation kernel K_zeros.

### 5.4 Why This Doesn't Directly Help

The Connes framework shows:
1. Zeta zeros = spectrum of an operator
2. RH ⟺ that operator is self-adjoint (spectrum real)
3. Self-adjointness ⟺ Weil positivity

But the heat kernel of the Connes operator does not directly give the Weil kernel K through a positive mixture. The relationship is:
- K_zeros involves the spectral decomposition of the operator (the zeros)
- The heat kernel involves e^{-tλ} for eigenvalues λ
- The Green's function involves 1/λ for eigenvalues λ

K_zeros has Fourier coefficients 1/(1/4 + γ²) at the zeros, which looks like a resolvent (Green's function) at spectral parameter 1/4, not like a heat kernel. More precisely:

$$K_{\text{zeros}}(x) \sim \sum_\gamma \frac{\cos(\gamma x)}{1/4 + \gamma^2}$$

Compare with the heat kernel: h_t(x) ~ Σ_γ e^{-tγ²} cos(γx).

And the resolvent: R_μ(x) ~ Σ_γ cos(γx)/(μ + γ²).

So **K_zeros is the resolvent kernel at μ = 1/4**, not a heat kernel. The resolvent CAN be written as a Laplace transform of the heat kernel:

$$\frac{1}{\mu + \gamma^2} = \int_0^\infty e^{-\mu t} e^{-\gamma^2 t}\, dt$$

Therefore:

$$K_{\text{zeros}}(x) = \int_0^\infty e^{-t/4} h_t^{(\text{zeros})}(x)\, dt$$

where h_t^{(zeros)}(x) = Σ_γ e^{-tγ²} cos(γx) is the "heat kernel on the zero spectrum." Since e^{-t/4} > 0 and h_t^{(zeros)} is PD (being a heat kernel), this shows K_zeros is PD — **but only under RH** (so that all γ are real and h_t^{(zeros)} is well-defined as a PD function).

---

## 6. The Euler Product and the Heat Kernel

### 6.1 Formal Euler Product of ζ as a Heat Trace

The Riemann zeta function has the Euler product:

$$\zeta(s) = \prod_p \frac{1}{1 - p^{-s}} \qquad (\text{Re}(s) > 1)$$

Taking logarithms:

$$\log\zeta(s) = -\sum_p \log(1 - p^{-s}) = \sum_p \sum_{m=1}^{\infty} \frac{p^{-ms}}{m}$$

The logarithmic derivative:

$$-\frac{\zeta'}{\zeta}(s) = \sum_p \frac{\log p}{p^s - 1} = \sum_p \sum_{m=1}^{\infty} \frac{\log p}{p^{ms}}$$

### 6.2 Zeta as a Spectral Zeta Function

If we define a Laplacian on T_A with eigenvalues λ(r) for r ∈ ℚ, the spectral zeta function is:

$$\zeta_\Delta(s) = \sum_{r \in \mathbb{Q} \setminus \{0\}} \lambda(r)^{-s}$$

For the archimedean Laplacian (λ_∞(n) = 4π²n² for n ∈ ℤ \ {0}):

$$\zeta_{\Delta_\infty}(s) = 2(4\pi^2)^{-s} \zeta_{\text{Riemann}}(2s)$$

This is the Epstein zeta function of the lattice ℤ. The factor of 2s (not s) means the Riemann zeta at s corresponds to the SQUARE ROOT of the Laplacian.

**For the full solenoid** with eigenvalues involving all rationals, the spectral zeta function over ℚ \ {0} does not directly equal the Riemann ζ(s), but is related through Dirichlet series involving arithmetic functions.

### 6.3 The Euler Product as Factorization of the Heat Trace

The heat trace on T_A (with product Laplacian) factorizes:

$$\Theta(t) := \text{Tr}(e^{t\Delta}) = \sum_{r \in \mathbb{Q}} e^{-t\lambda(r)} = \Theta_\infty(t) \cdot \prod_p \Theta_p(t)$$

where:
- Θ_∞(t) = Σ_{n∈ℤ} e^{-4π²n²t} = θ₃(0, e^{-4π²t}) (Jacobi theta)
- Θ_p(t) = 1 + Σ_{k=1}^∞ (p^k - p^{k-1}) e^{-tp^{kα}} (p-adic theta)

The Mellin transform of Θ(t) - 1 gives ζ_Δ(s):

$$\zeta_\Delta(s) = \frac{1}{\Gamma(s)} \int_0^\infty t^{s-1} (\Theta(t) - 1)\, dt$$

This connects the heat kernel on T_A to zeta-like functions, but the spectral zeta of the solenoid Laplacian is NOT the Riemann ζ(s) — it involves sums over ℚ with arithmetic weights, not sums over ℤ.

---

## 7. The Central Obstruction

### 7.1 Why K_Weil Cannot Be a Green's Function on T_A

The full Weil kernel K_Weil (bochner-proof.md §5–6) has the spectral form:

$$\hat{K}_{\text{Weil}}(\tau) = 2\,\text{Re}[\psi(1/4 + i\tau/2)] + 2\,\text{Re}[1/(-1/2+i\tau)] - \gamma/2 + (\text{pole terms})$$

From bochner-proof.md §7.2, K̂_Weil^(no poles)(0) ≈ -12.7. This is **negative**.

A Green's function G of any non-negative Laplacian Δ on T_A has Fourier coefficients Ĝ(r) = 1/λ(r) > 0 for r ≠ 0, making G automatically CPD-1.

**Therefore K_Weil (without pole terms) cannot be a Green's function on T_A**, since its Fourier transform takes negative values.

### 7.2 Why K_bg CAN Be (and Essentially IS) a Green's Function

K_bg has K̂_bg(ξ) = 2e^{-|ξ|/2}/(1 - e^{-2|ξ|}) > 0 for all ξ ≠ 0. This IS of the form 1/λ(ξ) for:

$$\lambda(\xi) = \frac{1 - e^{-2|\xi|}}{2e^{-|\xi|/2}} = \frac{e^{|\xi|/2}(1 - e^{-2|\xi|})}{2} = \frac{e^{|\xi|/2} - e^{-3|\xi|/2}}{2} = \sinh(|\xi|)$$

So K_bg is the Green's function of the operator whose eigenvalue at frequency ξ is sinh(|ξ|). On the solenoid with dual ℚ, this would be:

$$\Delta_{\text{bg}} \chi_r = -\sinh(|r|) \chi_r$$

This is a legitimate pseudodifferential operator on T_A (the eigenvalues sinh(|r|) grow as |r| → ∞ and vanish at r = 0), but it is NOT a differential operator — it has exponential growth rather than polynomial growth in the frequency.

**Interpretation:** K_bg is the resolvent of an exponential-type operator on T_A, related to but distinct from the standard Laplacian. This is consistent with K_bg involving the digamma function (which has logarithmic, not polynomial, growth).

### 7.3 The Role of the Pole Terms

The pole terms |F(0)|² + |F(1)|² in the Weil functional (bochner-proof.md §5.3) contribute distributional terms at ξ = 0 in the CPD-1 framework. They correspond to a rank-2 correction:

$$K_{\text{pole}}(x) = c_0 + c_1 e^x + c_1 e^{-x}$$

for appropriate constants c_0, c_1 (related to the poles of ζ at s = 0, 1).

In the heat-kernel framework, the pole terms would correspond to the zero modes of the Laplacian — the constant function (for s = 1 pole) and the character e^x (for the s = 0 pole). These are precisely the modes EXCLUDED from the Green's function (which inverts Δ only on the orthogonal complement of ker Δ).

**This is the geometric meaning of the CPD-1 condition:** it tests positivity of K after projecting out the zero modes, which corresponds to restricting to the primitive subspace (Σ c_i w_i = 0 annihilates the pole contributions).

---

## 8. Assessment

### 8.1 What Works

1. **K_bg as a Green's function.** K_bg is the Green's function (resolvent kernel) of the operator Δ_bg with eigenvalues sinh(|r|) on the solenoid T_A. This provides a geometric interpretation of K_bg as a potential function on T_A and an alternative proof of its CPD-1 property.

2. **K_zeros as a resolvent.** K_zeros = ∫₀^∞ e^{-t/4} h_t^{(zeros)} dt where h_t^{(zeros)} is the "heat kernel" on the zero spectrum. Under RH, this is a positive mixture of PD functions, hence PD. (But this is circular — RH is assumed.)

3. **Product structure.** The Euler product structure of ζ parallels the product decomposition of the heat kernel on T_A = ℝ/ℤ × Π_p ℤ_p. The local heat kernels are well-understood (circle diffusion + Vladimirov diffusion).

4. **Connes connection.** The prolate wave operator framework of Connes-Consani-Moscovici provides a rigorous spectral realization of zeta zeros. The Sonin space and Weil positivity at the archimedean place are proved within this framework.

### 8.2 What Doesn't Work

1. **K_Weil ≠ positive mixture of heat kernels.** The full Weil kernel K_Weil (without pole terms) has negative Fourier transform near ξ = 0, so it cannot be written as ∫ f(t) h_t dt with f ≥ 0 for ANY heat kernel h_t on T_A. The negative region near ξ = 0 is an essential feature (bochner-proof.md §7.2), not an artifact.

2. **The spectral zeta of T_A ≠ the Riemann ζ.** The heat trace on T_A gives zeta functions over ℚ, not the Riemann ζ(s) directly. The connection between the two involves the Mellin transform and change of variables that do not preserve the product structure cleanly.

3. **No unconditional PD proof.** Every route through the heat kernel to prove CPD-1 of K_Weil encounters the same obstruction: the pole terms are essential for compensating the spectral negativity near ξ = 0, and whether they always compensate is equivalent to RH.

### 8.3 The Most Promising Direction

The **Connes-Consani-Moscovici prolate framework** is the most promising "heat kernel-adjacent" approach, because:

- It provides a concrete self-adjoint operator W_Λ whose spectrum encodes zeta zeros
- Weil positivity at the archimedean place is PROVED (not conjectured)
- Extension to the semilocal case (finitely many primes) is underway
- The Sonin space provides the correct analogue of the "heat kernel compression" needed for positivity

The obstruction to completing the proof via this route is: extending from finitely many places to all places (the "infinite product" problem), which is analogous to the effective Bombieri-Vinogradov issue in the AMR framework.

### 8.4 Summary Table

| Approach | Target | Result | Obstruction |
|----------|--------|--------|-------------|
| K_bg = Green's fn on T_A | CPD-1 of K_bg | **Proved** (redundant with Fourier) | None |
| K_zeros as resolvent | PD of K_zeros | PD under RH only | Circular |
| K_Weil as heat mixture | CPD-1 of K_Weil | **Fails** | K̂_Weil < 0 near 0 |
| Euler product decomposition | Factor heat kernel | Structural parallel | Spectral ζ ≠ Riemann ζ |
| Connes prolate operator | Spectral realization | Archimedean case done | Extension to all places |

---

## 9. Technical Details

### 9.1 Convergence of the Green's Function on T_A

For the archimedean Laplacian, Σ_{r∈ℚ\{0}} 1/(4π²r²) diverges (the rationals are dense). The Green's function must be regularized:

$$G_\infty^{(\text{reg})}(x) = \lim_{N \to \infty} \left[\sum_{\substack{r \in \mathbb{Q} \\ 0 < |r| \leq N}} \frac{\chi_r(x)}{4\pi^2 r^2} - C(N)\right]$$

where C(N) is a divergent renormalization constant. This regularization is standard in potential theory on groups with dense duals.

### 9.2 The Bernstein-Widder Theorem Application

**Theorem (Bernstein-Widder).** A function f: (0,∞) → ℝ is completely monotone iff f(s) = ∫₀^∞ e^{-st} dμ(t) for a non-negative Borel measure μ.

K̂_bg(ξ) = 1/sinh(|ξ|) as a function of |ξ| ∈ (0,∞): Is it completely monotone?

$$\frac{d}{d\xi}\frac{1}{\sinh\xi} = -\frac{\cosh\xi}{\sinh^2\xi} < 0$$

$$\frac{d^2}{d\xi^2}\frac{1}{\sinh\xi} = \frac{2\cosh^2\xi - \sinh^2\xi}{\sinh^3\xi} = \frac{\cosh^2\xi + 1}{\sinh^3\xi} > 0$$

The alternating signs suggest complete monotonicity. Indeed, 1/sinh(ξ) = 2e^{-ξ}/(1-e^{-2ξ}) = 2Σ_{n=0}^∞ e^{-(2n+1)ξ} is a sum of exponentials with positive coefficients, confirming complete monotonicity. Therefore the Bernstein-Widder theorem applies, and K̂_bg IS a Laplace transform of a positive measure — but in the variable |ξ|, not ξ².

### 9.3 p-Adic Brownian Motion (Vladimirov Diffusion)

On ℤ_p with α = 2, the Vladimirov diffusion X_t is a Markov process that:
- Jumps between cosets of p^k ℤ_p at rate ~p^{2k}
- Has transition density h_t^{(p)}(x-y) given in §4.2
- Satisfies the recurrence relation for equilibration times: τ_k ~ p^{-2k} (time to mix within p^k ℤ_p)

The process is a "p-adic random walk" that makes larger jumps less frequently, in contrast to Euclidean Brownian motion which makes continuous small movements.

---

## 10. Connections to Other Approaches in the AMR Framework

### 10.1 Entropy-Positivity Duality

The heat kernel h_t on T_A defines a Markov semigroup. The entropy of this semigroup:

$$H(t) = -\int_{\mathbb{T}_\mathbb{A}} h_t \log h_t\, d\lambda$$

satisfies H'(t) ≤ 0 (entropy is non-increasing for the heat equation on a compact group). The long-time limit h_∞ = 1 (uniform/Haar measure) has maximal entropy.

This connects to the AMR entropy-positivity duality (entropy-positivity.md): the arithmetic entropy h_ar measures the mixing rate of the multiplication-by-p map, while the heat kernel entropy measures diffusion equilibration. Both converge to Haar measure, but via different dynamics (discrete multiplicative vs. continuous additive).

### 10.2 Lévy Processes

The heat kernel framework on T_A can be generalized to Lévy processes (processes with independent increments). On a compact abelian group, a Lévy process is characterized by its Lévy-Khintchine exponent:

$$\mathbb{E}[e^{i\langle r, X_t\rangle}] = e^{-t\psi(r)}$$

where ψ: ℚ → [0,∞) is a negative definite function. The Laplacian eigenvalues λ(r) correspond to ψ(r).

If K_Weil were the potential kernel of such a process (K̂_Weil = 1/ψ), then ψ = 1/K̂_Weil would need to be negative definite. But K̂_Weil takes negative values near ξ = 0 (§7.1), so 1/K̂_Weil has a pole, and ψ is not well-defined. This blocks the Lévy process interpretation of K_Weil.

---

## References

- Connes, A. (1999). Trace formula in noncommutative geometry and the zeros of the Riemann zeta function. [arXiv:math/9811068](https://arxiv.org/abs/math/9811068)
- Connes, A., Consani, C. (2020). Weil positivity and Trace formula, the archimedean place. [arXiv:2006.13771](https://arxiv.org/abs/2006.13771)
- Connes, A., Consani, C., Moscovici, H. (2024). Zeta zeros and prolate wave operators. [arXiv:2310.18423](https://arxiv.org/abs/2310.18423)
- Grigor'yan, A. (2023). Analysis on ultra-metric spaces via heat kernels. *p-Adic Numbers, Ultrametric Analysis and Applications*. [Springer](https://link.springer.com/article/10.1134/S2070046623030044)
- Kochubei, A.N. (2017). Linear and nonlinear heat equations on a p-adic ball. [arXiv:1708.03261](https://arxiv.org/abs/1708.03261)
- Zúñiga-Galindo, W.A. (2017). Parabolic-type equations and Markov processes on adeles. [Springer](https://link.springer.com/chapter/10.1007/978-3-319-46738-2_4)
- Chernoff, P.R. et al. Zeta functions, heat kernels, and spectral asymptotics on degenerating families of discrete tori. [Cambridge Core](https://www.cambridge.org/core/journals/nagoya-mathematical-journal/article/zeta-functions-heat-kernels-and-spectral-asymptotics-on-degenerating-families-of-discrete-tori/B184A6C302ABBAC52EDB57C136F93A1E)
- [bochner-proof.md](bochner-proof.md) — CPD-1 analysis of K_bg, K_zeros, K_Weil
- [schoenberg-attempt.md](schoenberg-attempt.md) — Schoenberg representation approach
- [entropy-positivity.md](../proofs/entropy-positivity.md) — Arithmetic entropy and positivity duality

---

*Generated as part of the AMR framework — Heat Kernel Exploration*
*Date: 2026-02-12*
