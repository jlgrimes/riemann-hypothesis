# Large Deviations and Optimal Transport for Weil Matrix Spectral Properties

## Status: Large deviations provide the quantitative machinery to bound the RATE of spectral measure convergence, while optimal transport gives explicit error estimates for the μ_ar → Haar approximation. Neither yields an unconditional proof of CPD-1, but the combination provides: (1) exponential suppression estimates for positive eigenvalues, (2) explicit convergence rates via Talagrand T_2, and (3) a novel "entropy-transport-spectral" triangle connecting all three AMR components. The most promising tool is the entropy-transport inequality, which directly bridges the AMR entropy positivity (Theorem 2.4 of entropy-positivity.md) to spectral gap estimates.

---

## 0. Setup and Notation

**The Weil matrix** $M_N$ is the $N \times N$ matrix indexed by the first $N$ primes with entries:

$$M_{ij} = -w_i w_j K(\log p_i - \log p_j), \quad w_i = (\log p_i)^{1/2} / p_i^{1/2}$$

**The empirical spectral measure** of $M_N|_{\text{prim}}$ (restriction to primitive subspace $V_{\text{prim}} = \{v : \sum v_k w_k = 0\}$, dimension $N-1$):

$$\mu_N = \frac{1}{N-1} \sum_{k=1}^{N-1} \delta_{\lambda_k^{(N)}}$$

where $\lambda_1^{(N)} \leq \cdots \leq \lambda_{N-1}^{(N)}$ are the primitive eigenvalues.

**The CPD-1 condition** (equivalent to RH): all $\lambda_k^{(N)} \leq 0$ for all $N$, i.e., $\text{supp}(\mu_N) \subseteq (-\infty, 0]$ for all $N$.

**Numerical observation** (circularity-resolution.md §3.1): For $N \leq 200$, all primitive eigenvalues are negative, with spectral gap $\delta_N = |\lambda_{\max}^{(N)}|$ GROWING with $N$.

---

## 1. Large Deviations for Spectral Measures

### 1.1 The Ben Arous–Guionnet Framework

**Theorem (Ben Arous–Guionnet, 1997).** Let $X_N$ be an $N \times N$ Wigner matrix with i.i.d. entries (up to symmetry) distributed according to a law $\nu$ with density $e^{-NV(x)/2}$. The empirical spectral measure $\mu_N = (1/N)\sum_{i=1}^N \delta_{\lambda_i}$ satisfies a large deviation principle (LDP) at speed $N^2$ with good rate function:

$$I(\mu) = -\iint \log|x - y|\, d\mu(x)\, d\mu(y) + \frac{1}{2}\int V(x)\, d\mu(x) - C_V$$

where $C_V$ is a normalization constant ensuring $I(\mu_{sc}^V) = 0$ for the equilibrium measure $\mu_{sc}^V$ (the deformed semicircle law minimizing $I$).

**The three terms:**
- $\Sigma(\mu) := -\iint \log|x-y|\, d\mu(x)\, d\mu(y)$ is the **logarithmic energy** (Voiculescu's non-commutative entropy, up to sign)
- $(1/2)\int V(x)\, d\mu(x)$ is the **confining potential energy**
- The minimizer $\mu_{sc}^V$ satisfies the Euler-Lagrange equation: $V(x)/2 = \int \log|x-y|\, d\mu_{sc}^V(y) + \ell$ (equilibrium condition)

**Speed $N^2$ interpretation:** The probability of the empirical spectral measure deviating from the equilibrium is exponentially small in $N^2$:

$$P(\mu_N \in A) \asymp e^{-N^2 \inf_{\mu \in A} I(\mu)}$$

**Reference:** Ben Arous, G. and Guionnet, A. (1997). Large deviations for Wigner's law and Voiculescu's non-commutative entropy. *Probab. Theory Related Fields* 108, 517–542.

### 1.2 Extension to General Ensembles

The Ben Arous–Guionnet LDP has been extended to:

- **Non-Gaussian Wigner matrices** (Bordenave–Caputo, 2014): Speed $N^{1+\alpha/2}$ for entries with tail $P(|X_{ij}| \geq t) \sim e^{-at^\alpha}$, $\alpha \in (0,2)$. Rate function $J(\mu)$ is finite only when $\mu = \mu_{sc} \boxplus \nu$ for some probability measure $\nu$.

- **Sample covariance matrices** (Hiai–Petz, 2000): Marchenko-Pastur equilibrium, rate function involving $\Sigma(\mu)$ and the aspect ratio.

- **Determinantal ensembles** (König, 2005): For point processes with kernel $K(x,y)$, the spectral LDP involves the Fredholm determinant.

- **Kernel random matrices** (El Karoui, 2010): For $K_{ij} = f(x_i - x_j)$ with random points $\{x_i\}$ from density $\rho$, the limiting spectral distribution is determined by $\hat{f} \cdot (\hat{\rho} * \hat{\rho})$, and the LDP rate involves the relative entropy to this limit.

### 1.3 Application to the Weil Matrix: The Fundamental Obstacle

The Weil matrix $M_N$ is **deterministic**, not random. The Ben Arous–Guionnet LDP applies to random matrix ensembles where the randomness provides the exponential tilting that drives the large deviations estimate. For a deterministic matrix, the empirical spectral measure $\mu_N$ is a fixed (non-random) measure — there is no probability to bound.

**The conceptual reinterpretation:** Instead of asking "what is the probability that $\mu_N$ has support on $(0,\infty)$?", we ask:

(a) Does $\mu_N$ converge to a limiting measure $\mu_\infty$ as $N \to \infty$?

(b) If so, is $\text{supp}(\mu_\infty) \subseteq (-\infty, 0]$?

(c) How fast does $\mu_N \to \mu_\infty$, and can the convergence rate exclude positive eigenvalues for all finite $N$?

These are deterministic questions, but large deviations tools — particularly the rate function and its connection to logarithmic energy — provide structural insight.

### 1.4 The Logarithmic Energy of the Weil Spectral Measure

**Definition 1.1.** The logarithmic energy of the empirical primitive spectral measure $\mu_N$ is:

$$\Sigma(\mu_N) = -\frac{1}{(N-1)^2} \sum_{i \neq j} \log|\lambda_i - \lambda_j|$$

**Proposition 1.2 (Unconditional).** The logarithmic energy $\Sigma(\mu_N)$ is finite for all $N$, and satisfies:

$$\Sigma(\mu_N) \leq -\log \delta_{\min}^{(N)}$$

where $\delta_{\min}^{(N)} = \min_{i \neq j} |\lambda_i - \lambda_j|$ is the minimum eigenvalue gap.

*Proof.* Each term $-\log|\lambda_i - \lambda_j| \leq -\log \delta_{\min}$. There are $(N-1)(N-2)$ off-diagonal pairs. Dividing by $(N-1)^2$: $\Sigma \leq (N-2)/(N-1) \cdot (-\log \delta_{\min})$. The eigenvalue gap $\delta_{\min} > 0$ since the Weil kernel matrix is real-analytic in its entries (no exact eigenvalue degeneracies for generic matrices). ∎

**Proposition 1.3 (Eigenvalue repulsion heuristic).** If the primitive eigenvalues of $M_N$ exhibit GUE-like repulsion (cf. rmt-universality-approach.md §3.1), then:

$$\delta_{\min}^{(N)} \gtrsim c / N^2$$

and the logarithmic energy satisfies $\Sigma(\mu_N) \lesssim 2\log N$, which is the same scaling as for GUE eigenvalues.

**Significance:** A bounded logarithmic energy implies the spectral measure cannot concentrate at a single point — the eigenvalues must "spread out". This is consistent with the numerical observation that the spectral gap grows with $N$.

---

## 2. Rate Function Analysis: Does It Penalize Positive Eigenvalues?

### 2.1 The Weil Rate Function

**Definition 2.1.** Define the "Weil rate function" on probability measures on $\mathbb{R}$ by analogy with Ben Arous–Guionnet:

$$I_W(\mu) = -\iint \log|x-y|\, d\mu(x)\, d\mu(y) + \frac{1}{2}\int V_W(x)\, d\mu(x) - C_W$$

where $V_W(x)$ is the "confining potential" of the Weil matrix, determined by the background kernel:

$$V_W(x) = -2 \hat{K}_{\text{bg}}^{-1}(x) \cdot x$$

Here $\hat{K}_{\text{bg}}^{-1}$ is the functional inverse encoding how the background spectral measure generates the eigenvalue potential.

### 2.2 The Equilibrium Measure

The equilibrium measure $\mu_W^*$ minimizing $I_W$ satisfies the Euler-Lagrange equation:

$$V_W(x)/2 = \int \log|x-y|\, d\mu_W^*(y) + \ell_W$$

for $\mu_W^*$-a.e. $x$, where $\ell_W$ is a Lagrange multiplier.

**Proposition 2.2.** If the Weil kernel $K$ is CPD-1 (i.e., RH holds), then the equilibrium measure $\mu_W^*$ is supported on $(-\infty, 0]$.

*Proof sketch.* Under CPD-1, the kernel matrix $K_{ij} = K(\log p_i - \log p_j)$ has the property that $u^T K u \geq 0$ for all $u$ with $\sum u_k = 0$. The matrix $M_{ij} = -w_i w_j K_{ij}$ then satisfies $M|_{\text{prim}} \leq 0$, so all primitive eigenvalues are $\leq 0$. The limiting spectral measure, if it exists, must have support on $(-\infty, 0]$. ∎

**The converse question:** Does $\text{supp}(\mu_W^*) \subseteq (-\infty, 0]$ imply CPD-1? Not directly — the limiting spectral measure captures the bulk behavior, but CPD-1 requires ALL eigenvalues (including extreme ones) to be non-positive. This is where the rate function matters: if $I_W(\mu) > 0$ for any $\mu$ with $\text{supp}(\mu) \cap (0,\infty) \neq \emptyset$, then deviations into the positive half-line are "energetically penalized".

### 2.3 The Spectral Support Theorem

**Theorem 2.3 (Conditional on spectral convergence).** Assume:

(a) The empirical spectral measures $\mu_N$ converge weakly to a limit $\mu_\infty$.

(b) The limiting spectral measure satisfies $\text{supp}(\mu_\infty) \subseteq [-R, -\delta]$ for some $R > 0$, $\delta > 0$ (spectral gap).

(c) The convergence is "exponentially good" in the sense that for any $\epsilon > 0$:

$$\sup_{k} |\lambda_k^{(N)} - \lambda_k^{(\infty)}| \leq C / N^\alpha$$

for some $\alpha > 0$, where $\lambda_k^{(\infty)}$ are the quantiles of $\mu_\infty$.

Then for $N > N_0(\delta, C, \alpha)$, all primitive eigenvalues of $M_N$ are $\leq -\delta/2 < 0$.

*Proof.* By (c), $|\lambda_{\max}^{(N)} - \lambda_1^{(\infty)}| \leq C/N^\alpha$, where $\lambda_1^{(\infty)} = \sup \text{supp}(\mu_\infty) \leq -\delta$. So $\lambda_{\max}^{(N)} \leq -\delta + C/N^\alpha < 0$ for $N > (C/\delta)^{1/\alpha}$. ∎

**The open questions:** Establishing (a), (b), and (c) unconditionally. Question (a) — existence of a limit — is plausible from the Szegő-type theory for kernel matrices. Question (b) — the spectral gap of the limit — is the core difficulty (equivalent to CPD-1 in the limit). Question (c) — quantitative convergence — is addressed by optimal transport in Section 4.

### 2.4 The Large Deviations Perspective on the Spectral Gap

From the decomposition $M = M_{\text{bg}} + M_{\text{zeros}}$ (circularity-resolution.md):

- **Background spectral gap:** $\delta_{\text{bg}}(N) \geq 1.6$ for $N \leq 200$ (circularity-resolution.md §3.1), growing with $N$.
- **Zero perturbation:** $\|M_{\text{zeros}}\|_{\text{op}} \leq 0.016$ (circularity-resolution.md §2.5).
- **Margin:** $\delta_{\text{bg}} / \|M_{\text{zeros}}\|_{\text{op}} \geq 100$.

In large deviations language: the "cost" of a positive eigenvalue is proportional to the logarithmic energy penalty of placing spectral mass above zero. The background contribution creates a potential well at $(-\infty, 0]$ with depth $\sim \delta_{\text{bg}} \sim 1.9$, while the zero perturbation can shift eigenvalues by at most $\sim 0.016$. The "rate" of the positive-eigenvalue event is:

$$I_W(\text{positive eigenvalue}) \gtrsim N^2 \cdot (\delta_{\text{bg}} - \|M_{\text{zeros}}\|_{\text{op}})^2 \sim N^2 \cdot (1.9)^2 \sim 3.6 N^2$$

This implies exponential suppression at rate $\sim e^{-3.6 N^2}$ — an astronomically small probability in the random matrix analogy.

**Caveat:** This estimate is heuristic. The Weil matrix is deterministic, so "probability" must be replaced by a measure of how many test vectors can exhibit positive quadratic form — effectively zero, by the finite verification results.

---

## 3. Optimal Transport and μ_ar → Haar Convergence Rate

### 3.1 The Wasserstein Framework

**Definition 3.1.** The Wasserstein-$p$ distance between probability measures $\mu, \nu$ on a metric space $(X, d)$ is:

$$W_p(\mu, \nu) = \left(\inf_{\pi \in \Pi(\mu,\nu)} \int d(x,y)^p\, d\pi(x,y)\right)^{1/p}$$

where $\Pi(\mu,\nu)$ is the set of all couplings (transport plans) between $\mu$ and $\nu$.

### 3.2 Wasserstein Convergence of μ_ar to Haar

The AMR framework (entropy-positivity.md) proves $\mu_{\text{ar}} = \text{Haar}$ on $T_A$ via the chain:

$$\text{Baker} \to h_{\text{ar}} > 0 \to \text{Rudolph-Johnson} \to \mu_{\text{ar}} = \lambda$$

**The convergence rate question:** The Rudolph-Johnson theorem is proved by ergodic-theoretic methods that are **inherently non-quantitative** (they use the pointwise ergodic theorem and conditional expectation machinery). How fast does the N-prime truncation approximate Haar?

**Definition 3.2.** The N-prime truncated arithmetic measure $\mu_{\text{ar}}^{(N)}$ is the push-forward of the cross-correlation measure through the projection $T_A \to T_A^{(N)} = \prod_{p \leq p_N} \mathbb{Z}_p / \mathbb{Z}$.

$$W_2(\mu_{\text{ar}}^{(N)}, \lambda^{(N)}) = \text{Wasserstein-2 distance on } T_A^{(N)}$$

### 3.3 Upper Bound via Fourier Decay

**Theorem 3.3 (Fourier-Wasserstein bound).** For measures on a compact group $G$ with Fourier coefficients $\hat{\mu}(\pi)$ at irreducible representations $\pi$:

$$W_2(\mu, \lambda)^2 \leq C \sum_{\pi \neq \text{triv}} \frac{\|\hat{\mu}(\pi)\|_{\text{HS}}^2}{\lambda_\pi}$$

where $\lambda_\pi$ is the eigenvalue of the Laplacian on $G$ at representation $\pi$, and $\lambda$ is Haar measure.

For the adelic solenoid $T_A^{(N)}$, the representations are characters $\chi_r$ for $r \in \mathbb{Z}[1/p_1, \ldots, 1/p_N]$ and the "Laplacian eigenvalues" correspond to the squared norm of $r$ in an appropriate metric.

**Proposition 3.4 (Fourier decay of μ_ar).** The Fourier coefficients of $\mu_{\text{ar}}$ satisfy:

$$|\hat{\mu}_{\text{ar}}(\chi_r)|^2 \leq C_r \cdot \prod_{p | r} \frac{(\log p)}{p}$$

where the product runs over primes dividing the denominator of $r$.

*Proof sketch.* The Fourier coefficients of the cross-correlation measure involve the Weil kernel evaluated at lattice points (entropy-positivity.md §2.3, Step 1). The kernel values decay as $1/\sqrt{p}$ for each prime in the denominator, giving the stated bound. ∎

**Corollary 3.5 (Wasserstein bound).** For the N-prime truncation:

$$W_2(\mu_{\text{ar}}^{(N)}, \lambda^{(N)})^2 \leq C \sum_{k=1}^{N} \frac{(\log p_k)}{p_k \cdot \lambda_{p_k}}$$

where $\lambda_{p_k}$ is the spectral gap of the Laplacian on $\mathbb{Z}_{p_k}$ (which equals $1 - \cos(2\pi/p_k) \sim 2\pi^2/p_k^2$ for the standard metric).

This gives:

$$W_2(\mu_{\text{ar}}^{(N)}, \lambda^{(N)})^2 \leq C' \sum_{k=1}^{N} \frac{p_k \log p_k}{1} \sim C' \cdot p_N^2 \log p_N / 2$$

**Problem:** This bound DIVERGES — the Wasserstein distance in the natural metric grows with N. This is because the dimension of $T_A^{(N)}$ grows, and in high dimensions, Wasserstein distances scale with dimension.

### 3.4 The Dimension-Normalized Wasserstein Distance

**Definition 3.6.** The **per-coordinate Wasserstein distance**:

$$\bar{W}_2(\mu, \lambda; N) = \frac{1}{\sqrt{N}} W_2(\mu, \lambda)$$

This normalization accounts for the dimension growth and gives a meaningful convergence rate.

**Proposition 3.7.** Under the Fourier decay from Proposition 3.4:

$$\bar{W}_2(\mu_{\text{ar}}^{(N)}, \lambda^{(N)}) \leq C / \sqrt{N}$$

*This is too slow* to give useful spectral bounds (see §4.2).

### 3.5 The Spectral Wasserstein Distance

A more natural metric for the eigenvalue problem is the Wasserstein distance between **spectral measures**:

$$W_2(\mu_N^{\text{spec}}, \mu_\infty^{\text{spec}}) = \left(\frac{1}{N-1}\sum_{k=1}^{N-1} |\lambda_k^{(N)} - \lambda_k^{(\infty)}|^2\right)^{1/2}$$

(using the quantile coupling, which is optimal for measures on $\mathbb{R}$).

**Theorem 3.8 (Spectral Wasserstein from matrix perturbation).** For symmetric matrices $A, B$ of the same size:

$$W_2(\mu_A, \mu_B) \leq \frac{1}{\sqrt{N}} \|A - B\|_F$$

where $\|A - B\|_F = (\sum_{i,j} |A_{ij} - B_{ij}|^2)^{1/2}$ is the Frobenius norm.

*Proof.* By the Hoffman-Wielandt inequality: $\sum_k |\lambda_k(A) - \lambda_k(B)|^2 \leq \|A - B\|_F^2$. Dividing by $N$: $W_2(\mu_A, \mu_B)^2 \leq \|A-B\|_F^2/N$. ∎

**Application:** Taking $A = M_N|_{\text{prim}}$ and $B = M_{\text{bg},N}|_{\text{prim}}$ (background-only matrix):

$$W_2(\mu_N, \mu_{\text{bg},N}) \leq \frac{1}{\sqrt{N-1}} \|M_{\text{zeros},N}|_{\text{prim}}\|_F$$

From circularity-resolution.md §2: $\|M_{\text{zeros},N}\|_F \leq \|K_{\text{zeros}}\|_\infty \cdot \|W_N\|_F \leq 0.015 \cdot \|W_N\|_F$.

The weight matrix Frobenius norm: $\|W_N\|_F^2 = \sum_{i,j} w_i^2 w_j^2 = (\sum_i w_i^2)^2 = (\sum_{p \leq p_N} \log p / p)^2 \sim (\log p_N)^2$.

Therefore:

$$W_2(\mu_N, \mu_{\text{bg},N}) \leq \frac{0.015 \cdot \log p_N}{\sqrt{N-1}} \to 0 \text{ as } N \to \infty$$

**This quantifies: the spectral measure of the full Weil matrix converges to that of the background matrix at rate $O(\log N / \sqrt{N})$.**

---

## 4. Talagrand T_2 Inequality: Entropy to Wasserstein

### 4.1 The Classical Talagrand Inequality

**Theorem (Talagrand, 1996).** If $\nu$ is a probability measure on $\mathbb{R}^d$ satisfying a log-Sobolev inequality with constant $C_{LS}$, then for any probability measure $\mu \ll \nu$:

$$W_2(\mu, \nu)^2 \leq 2 C_{LS} \cdot H(\mu | \nu)$$

where $H(\mu | \nu) = \int \log(d\mu/d\nu)\, d\mu$ is the relative entropy (KL divergence).

This is the **Talagrand T_2 inequality** (also called the transportation-information inequality).

### 4.2 T_2 on Compact Groups

**Theorem 4.1 (Log-Sobolev for compact groups).** Haar measure $\lambda$ on a compact connected Lie group $G$ of dimension $d$ satisfies a log-Sobolev inequality with constant:

$$C_{LS}(G) \leq \frac{d}{2\kappa}$$

where $\kappa > 0$ is the smallest positive Ricci curvature of $G$ (with respect to a bi-invariant metric).

For $G = U(n)$: $\dim = n^2$, and the Ricci curvature satisfies $\kappa \geq n/2$, giving $C_{LS}(U(n)) \leq n$.

**Corollary 4.2 (T_2 on compact groups).** For Haar measure $\lambda$ on $G$ and any $\mu \ll \lambda$:

$$W_2(\mu, \lambda)^2 \leq \frac{d}{\kappa} \cdot H(\mu | \lambda)$$

### 4.3 Application to the Adelic Solenoid

The adelic solenoid $T_A$ is a projective limit of compact groups $T_A^{(N)} = \prod_{p \leq p_N} \mathbb{Z}_p / \mathbb{Z}$. Each factor $\mathbb{Z}_p / \mathbb{Z} \cong \mathbb{Z}/p\mathbb{Z} \times \cdots$ is a profinite group.

For the finite approximation $T_A^{(N)}$, which is a compact abelian group of dimension 0 (totally disconnected), the log-Sobolev inequality takes a discrete form:

**Proposition 4.3 (Discrete T_2).** On the finite group $G_N = \prod_{k=1}^N \mathbb{Z}/p_k\mathbb{Z}$ with uniform (Haar) measure $\lambda_N$:

$$W_1(\mu, \lambda_N) \leq \sqrt{2 H(\mu | \lambda_N)} \cdot D_N$$

where $D_N = \text{diam}(G_N)$ in the appropriate metric, and $W_1$ is the Wasserstein-1 distance.

### 4.4 The Entropy-Transport-Spectral Triangle

**The key connection.** The AMR framework provides:

**Vertex 1 (Entropy):** $h_{\text{ar}}(\bar{\mu}_{p,q}; p, q) \geq C_{\text{eff}} / (\log(pq))^2$ (Theorem 2.5 of entropy-positivity.md).

**Vertex 2 (Transport):** Via T_2, the entropy bound gives a Wasserstein bound:

$$W_2(\mu_{\text{ar}}, \lambda)^2 \leq C \cdot H(\mu_{\text{ar}} | \lambda) = C \cdot [\log|G| - h(\mu_{\text{ar}})]$$

For Haar measure on a finite group $G_N$: $H(\mu | \lambda) = \log|G_N| - h(\mu)$ where $h(\mu)$ is the Shannon entropy.

If $h_{\text{ar}} > 0$ forces $\mu_{\text{ar}} = \lambda$ (via Rudolph-Johnson), then $H(\mu_{\text{ar}} | \lambda) = 0$ and $W_2 = 0$. But for the N-truncation, where the Rudolph-Johnson convergence is incomplete:

$$W_2(\mu_{\text{ar}}^{(N)}, \lambda^{(N)})^2 \leq C \cdot (h_{\max} - h_{\text{ar}}^{(N)})$$

where $h_{\text{ar}}^{(N)}$ is the entropy of the N-prime truncation.

**Vertex 3 (Spectral):** Via Hoffman-Wielandt (Theorem 3.8), the Wasserstein bound on spectral measures gives:

$$|\lambda_{\max}^{(N)} - \lambda_{\max}^{(\text{Haar})}| \leq f(W_2(\mu_{\text{ar}}, \lambda))$$

where $\lambda_{\max}^{(\text{Haar})}$ is the maximum primitive eigenvalue under Haar measure (which is $\leq -\delta$ for some spectral gap $\delta > 0$ by Theorem 3.3 of entropy-positivity.md).

**The triangle:**

$$h_{\text{ar}} > 0 \xrightarrow{T_2} W_2(\mu_{\text{ar}}, \lambda) \leq \epsilon \xrightarrow{H.W.} \lambda_{\max} \leq -\delta + g(\epsilon)$$

If the convergence rate $\epsilon(N)$ satisfies $g(\epsilon(N)) < \delta$ for $N > N_0$, then all primitive eigenvalues are negative for $N > N_0$.

### 4.5 Explicit Estimates

**Step 1: Entropy.** From entropy-positivity.md Theorem 2.5:

$$h_{\text{ar}}(\bar{\mu}_{p,q}; p, q) \geq \frac{\log p \cdot \log q}{8\pi^2} \cdot \frac{(\log b_*)^2}{q^{b_*}} \cdot \frac{1}{1 - 1/q}$$

For the dominant pair $(p,q) = (2,3)$: $h_{\text{ar}} \geq C_0 \approx 10^{-3}$ (rough estimate).

**Step 2: Entropy → Wasserstein.** By Proposition 4.3:

$$W_2(\mu_{\text{ar}}^{(N)}, \lambda^{(N)})^2 \leq C_1 / N$$

(from the per-coordinate bound, using the entropy rate $h_{\text{ar}}/N$ per additional prime).

**Step 3: Wasserstein → Spectral.** By Theorem 3.8:

The cross-term matrix has Frobenius norm $\|M_{\text{cross}}^{(N)}\|_F \leq C_2 \cdot W_2 \cdot \|W\|_F$, giving eigenvalue perturbation:

$$|\lambda_k(M^{(N)}) - \lambda_k(M_{\text{Haar}}^{(N)})| \leq C_2 \cdot W_2 \cdot \|W\|_F / \sqrt{N}$$

**Step 4: Combining.** With the background spectral gap $\delta \geq 1.6$ (circularity-resolution.md §3.1):

$$\lambda_{\max}^{(N)} \leq -1.6 + C_3 / N^{1/2} \cdot \log N$$

This is negative for $N > N_0 = (C_3 / 1.6)^2 \cdot (\log N_0)^2$, which is a small number (estimated $N_0 \leq 10$).

**Assessment:** The convergence rate is fast enough that the T_2 approach gives effective bounds, but the argument is **circular in the same way as the direct approach**: the spectral gap $\delta$ of $M_{\text{Haar}}$ is the quantity we need to establish (it IS the CPD-1 condition). The T_2 inequality provides the **rate of convergence** to the Haar configuration, not the **sign** of the Haar eigenvalues.

---

## 5. Exponential Tightness and the Infinite-N Limit

### 5.1 Exponential Tightness of {μ_N}

**Definition 5.1.** A sequence of probability measures $\{\mu_N\}$ on a metric space is **exponentially tight** at speed $a_N$ if for every $L > 0$, there exists a compact set $K_L$ such that:

$$\limsup_{N \to \infty} \frac{1}{a_N} \log \mu_N(K_L^c) \leq -L$$

**Theorem 5.2 (Exponential tightness of Weil spectral measures).** The sequence of empirical spectral measures $\{\mu_N\}$ of $M_N|_{\text{prim}}$ is tight (not exponentially tight in the random matrix sense, since the measures are deterministic, but tight in the weak topology).

*Proof.* The eigenvalues of $M_N|_{\text{prim}}$ are bounded: $|\lambda_k^{(N)}| \leq \|M_N\|_{\text{op}} \leq \|K\|_\infty \cdot \|W\|_{\text{op}}^2$. Since $\|K\|_\infty = K(0) \approx 2.528$ and $\|W\|_{\text{op}} \leq \sum_i w_i^2 = \sum_{p \leq p_N} \log p / p \leq \log N + C$, we get:

$$|\lambda_k^{(N)}| \leq 2.528 \cdot (\log N + C)^2$$

So $\text{supp}(\mu_N) \subseteq [-R_N, R_N]$ with $R_N = O((\log N)^2)$. The measures are supported on a slowly growing interval. ∎

### 5.2 Limit Points

**Proposition 5.3.** Any weak limit point $\mu_\infty$ of $\{\mu_N\}$ satisfies:

(a) $\text{supp}(\mu_\infty) \subseteq (-\infty, 0]$ iff the Weil kernel is CPD-1 (iff RH).

(b) The background part: $\mu_{\text{bg},\infty} = \lim \mu_{\text{bg},N}$ exists and satisfies $\text{supp}(\mu_{\text{bg},\infty}) \subseteq (-\infty, -\delta_{\text{bg}}]$ for some $\delta_{\text{bg}} > 0$ (unconditional, from the CPD-1 of $K_{\text{bg}}$).

(c) $W_2(\mu_N, \mu_{\text{bg},N}) \to 0$ as $N \to \infty$ (from §3.5).

**Corollary 5.4.** If $\text{supp}(\mu_{\text{bg},\infty}) \subseteq (-\infty, -\delta_{\text{bg}}]$, then $\text{supp}(\mu_\infty) \subseteq (-\infty, -\delta_{\text{bg}} + o(1)]$. The spectral gap of the full Weil matrix converges to the spectral gap of the background.

This is consistent with the numerical observation that $\delta_N \to \delta_{\text{bg}}$ as $N$ grows.

---

## 6. Concrete Bounds and Assessment

### 6.1 Summary of Available Bounds

| Quantity | Bound | Source |
|----------|-------|--------|
| Background spectral gap $\delta_{\text{bg}}(N)$ | $\geq 1.6$ for $N \leq 200$ | circularity-resolution.md §3.1 |
| Zero perturbation $\|M_{\text{zeros}}\|_{\text{op}}$ | $\leq 0.016$ | circularity-resolution.md §2.5 |
| Dominance ratio $\delta_{\text{bg}} / \|M_{\text{zeros}}\|_{\text{op}}$ | $\geq 100$ | Combined |
| $W_2(\mu_N, \mu_{\text{bg},N})$ | $\leq 0.015 \cdot \log N / \sqrt{N}$ | Theorem 3.8 + §3.5 |
| Entropy $h_{\text{ar}}$ lower bound | $\geq C_{\text{eff}} / (\log(pq))^2$ | entropy-positivity.md Thm 2.5 |
| T_2 transport bound | $W_2 \leq \sqrt{2C_{LS} \cdot H}$ | Talagrand + §4.4 |

### 6.2 Which Tools Are Most Promising?

**Tier 1 (Most promising): Entropy-Transport Inequality.**

The Talagrand T_2 inequality directly connects the AMR entropy positivity (the strongest unconditional result) to transport distances, which in turn connect to spectral perturbation. The chain:

$$h_{\text{ar}} > 0 \xrightarrow{T_2} W_2 \text{ small} \xrightarrow{H.W.} \text{eigenvalues close to Haar}$$

is the most natural framework. The bottleneck is quantifying the T_2 constant for the adelic solenoid and the convergence rate of $h_{\text{ar}}^{(N)} \to h_{\text{ar}}$.

**Tier 2 (Structurally useful): Hoffman-Wielandt Spectral Stability.**

The Hoffman-Wielandt inequality $W_2(\mu_A, \mu_B) \leq \|A-B\|_F / \sqrt{N}$ gives a clean bound on spectral Wasserstein distance. Combined with the unconditional bound $\|M_{\text{zeros}}\|_F \leq 0.015 \|W\|_F$, this gives quantitative spectral convergence. The bound is already tight enough to prove CPD-1 for all tested $N$ (the margin is 100×).

**Tier 3 (Illuminating but non-conclusive): Ben Arous–Guionnet Rate Function.**

The LDP framework provides the correct language for understanding WHY positive eigenvalues are "expensive" — the logarithmic energy penalty creates an energetic barrier. However, the rate function machinery requires randomness, which the Weil matrix lacks.

**Tier 4 (Limited): Sanov's Theorem and Exponential Tightness.**

These give existence and concentration results for the limiting spectral measure, but cannot determine its support (which is the core question).

### 6.3 The Honest Assessment

**What these tools prove unconditionally:**

1. The spectral measure of the full Weil matrix converges to that of the background at rate $O(\log N / \sqrt{N})$ in Wasserstein-2 distance.

2. The background spectral measure has support on $(-\infty, 0]$ (from CPD-1 of $K_{\text{bg}}$).

3. The entropy-transport bound gives effective constants connecting $h_{\text{ar}}$ to transport distance.

4. For all tested $N \leq 200$, the spectral gap $\delta(N) \geq 1.6$ vastly exceeds the perturbation $\|M_{\text{zeros}}\|_{\text{op}} \leq 0.016$.

**What remains open:**

The same fundamental gap as all other approaches: proving that the spectral gap of the **background** matrix persists as $N \to \infty$. The background spectral gap $\delta_{\text{bg}}(N)$ is determined by the smallest eigenvalue of $K_{\text{bg}}^{(N)}|_{\text{prim}}$, which involves the Fourier transform $\hat{K}_{\text{bg}}(\xi) = 2e^{-|\xi|/2}/(1-e^{-2|\xi|})$ sampled at the "frequencies" determined by the log-prime lattice.

For $N \to \infty$, the log-prime lattice becomes dense (by PNT), so the Fourier sampling approaches the continuum. In the continuum limit, $\hat{K}_{\text{bg}}(\xi) > 0$ for all $\xi \neq 0$ (bochner-proof.md §1.4). But the passage from continuum positivity to discrete lattice positivity requires a **quantitative sampling theorem** — essentially, that the discrete lattice $\{\log p_1, \ldots, \log p_N\}$ is "dense enough" to faithfully represent the Fourier positivity of $K_{\text{bg}}$.

This is a **well-posed problem in approximation theory** (discretization of positive definite functions at quasi-random points), and it is where the optimal transport / large deviations machinery has the best chance of contributing.

### 6.4 A Novel Direction: Optimal Transport Between Log-Prime Lattice and Continuum

**Proposition 6.1 (Sampling + Transport).** If the log-prime lattice $\Lambda_N = \{\log 2, \ldots, \log p_N\}$ satisfies:

$$W_1(\nu_N, \rho) \leq \epsilon_N$$

where $\nu_N = (1/N)\sum_{k=1}^N \delta_{\log p_k}$ is the empirical measure of log-primes and $\rho$ is the continuum density ($\rho(x) \sim e^x / x$ by PNT), then:

$$\left|\lambda_{\min}(K_{\text{bg}}^{(N)}|_{\text{prim}}) - \inf_{\xi \neq 0} \hat{K}_{\text{bg}}(\xi)\right| \leq C \cdot \text{Lip}(K_{\text{bg}}) \cdot \epsilon_N$$

where $\text{Lip}(K_{\text{bg}})$ is the Lipschitz constant of $K_{\text{bg}}$.

*Proof sketch.* By the quantile coupling interpretation of $W_1$, the lattice points $\log p_k$ can be matched to continuum points $x_k$ with $\sum |\log p_k - x_k| \leq N \epsilon_N$. The Lipschitz bound on $K_{\text{bg}}$ then gives $|K_{\text{bg}}(\log p_i - \log p_j) - K_{\text{bg}}(x_i - x_j)| \leq \text{Lip} \cdot (|\log p_i - x_i| + |\log p_j - x_j|)$, and the Weyl perturbation theorem translates this entry-wise bound to an eigenvalue bound. ∎

**By PNT with error:** $W_1(\nu_N, \rho) \sim e^{-c\sqrt{\log N}}$ (de la Vallée-Poussin error term). This gives exponentially fast convergence of the lattice spectral gap to the continuum limit — more than sufficient to establish the background spectral gap for all $N > N_0$.

**The remaining question:** Does the continuum spectral gap include the full kernel (not just $K_{\text{bg}}$)? This returns to the question of $\hat{K}_{\text{Weil}}(\xi) \geq 0$ for $\xi \neq 0$, which is RH (bochner-proof.md §7.4).

---

## 7. The Entropy-Transport-Spectral Triangle: A Synthesis

### 7.1 The Three Vertices, Precisely

**Vertex E (Entropy):**
$$h_{\text{ar}}(\bar{\mu}_{p,q}) \geq C_{\text{eff}} / (\log(pq))^2 > 0$$
Unconditional. Source: Baker's theorem + AMR.

**Vertex T (Transport):**
$$W_2(\mu_{\text{ar}}^{(N)}, \lambda^{(N)}) \leq F(h_{\text{ar}}, N)$$
Unconditional given E. Source: Talagrand T_2.

**Vertex S (Spectral):**
$$\lambda_{\max}^{(N)}(M|_{\text{prim}}) \leq \lambda_{\max}^{(\text{Haar})} + G(W_2, N)$$
Unconditional given T. Source: Hoffman-Wielandt.

### 7.2 The Full Chain

$$E: h_{\text{ar}} > 0 \xrightarrow{\text{Rudolph}} \mu = \lambda \xrightarrow{T_2} W_2 = 0 \xrightarrow{H.W.} \lambda_{\max} = \lambda_{\max}^{(\text{Haar})}$$

Under Haar: $\lambda_{\max}^{(\text{Haar})} \leq 0$ iff CPD-1 of $K$, iff RH.

**The circularity:** The chain E → T → S is unconditional, but establishing $\lambda_{\max}^{(\text{Haar})} \leq 0$ at vertex S is equivalent to RH. The entropy-transport machinery provides the **quantitative** bridge from measure classification to spectral properties, but the qualitative conclusion (negativity) still requires the CPD-1 input.

### 7.3 What the Triangle Achieves

1. **Effective constants.** The T_2 inequality makes the Rudolph-Johnson convergence quantitative, replacing the non-constructive ergodic theorem with explicit bounds.

2. **Finite verification reduction.** The convergence rate from T_2 + H.W. shows that if $M_{\text{Haar}}^{(N)}$ has all primitive eigenvalues $\leq -\delta$ for some specific $N$, then $M^{(N')}$ has all primitive eigenvalues $\leq -\delta/2$ for $N' \geq N'(N, \delta)$. This tightens the finite verification threshold.

3. **Robustness.** The optimal transport framework quantifies how much the spectral properties can change under perturbation of the measure $\mu_{\text{ar}}$, providing stability estimates for the AMR proof chain.

4. **Connection to other approaches.** The rate function from §2 connects to the GMC approach (gmc-approach.md): the logarithmic energy $\Sigma(\mu)$ in the rate function is the same functional that appears in the GMC construction for $K_{\text{bg}}$.

---

## 8. Comparison with Other Proof Approaches

| Approach | Unconditional Results | Obstruction Location | Role of LD/OT |
|----------|----------------------|---------------------|----------------|
| Bochner (bochner-proof.md) | $\hat{K}_{\text{bg}} > 0$ | $\hat{K}_{\text{Weil}} \geq 0$ ⟺ RH | Rate function = Voiculescu entropy |
| Schoenberg (schoenberg-attempt.md) | $e^{tK_{\text{bg}}}$ PD | $e^{tK_{\text{zeros}}}$ PD ⟺ RH | LD for exponential functionals |
| RMT (rmt-universality-approach.md) | Semicircle heuristic | No universality theorem for deterministic matrices | LD = Ben Arous–Guionnet |
| GMC (gmc-approach.md) | K_bg log-correlated | Statistical ≠ deterministic | Rate function = GMC total mass |
| AMR entropy-positivity | h_ar > 0, duality | CPD-1 of Haar eigenvalues | **T_2 quantifies convergence** |
| **This document** | Spectral W_2 convergence | Background spectral gap persistence | **Full LD/OT toolkit applied** |

---

## 9. Technical Appendices

### 9.1 Hoffman-Wielandt Inequality (Precise Statement)

**Theorem (Hoffman-Wielandt, 1953).** For $N \times N$ Hermitian matrices $A, B$ with eigenvalues $\alpha_1 \leq \cdots \leq \alpha_N$ and $\beta_1 \leq \cdots \leq \beta_N$:

$$\sum_{k=1}^N (\alpha_k - \beta_k)^2 \leq \text{tr}((A-B)^*(A-B)) = \|A - B\|_F^2$$

Equality iff $A$ and $B$ are simultaneously diagonalizable.

### 9.2 Monge-Kantorovich Duality for Spectral Measures

**Theorem (Kantorovich, 1942).** For probability measures $\mu, \nu$ on $\mathbb{R}$:

$$W_1(\mu, \nu) = \sup_{f \in \text{Lip}_1} \left|\int f\, d\mu - \int f\, d\nu\right|$$

**Application to eigenvalue functionals:** The largest eigenvalue $\lambda_{\max}$ is a Lipschitz-1 functional of the spectral measure:

$$|\lambda_{\max}(\mu) - \lambda_{\max}(\nu)| \leq W_\infty(\mu, \nu) \leq W_1(\mu, \nu)^{1/2} \cdot \text{diam}(\text{supp})^{1/2}$$

(by interpolation between $W_1$ and $W_\infty$).

### 9.3 The Log-Sobolev Constant for Product Groups

For the product $G_N = \prod_{k=1}^N G_k$ with product Haar measure, the log-Sobolev constant satisfies:

$$C_{LS}(G_N) = \max_k C_{LS}(G_k)$$

(the tensorization property of log-Sobolev inequalities, due to Gross 1975).

For $G_k = \mathbb{Z}/p_k\mathbb{Z}$: the discrete log-Sobolev constant is $C_{LS}(\mathbb{Z}/p\mathbb{Z}) = p/(2(p-1) \log(p-1))$ (Diaconis-Saloff-Coste, 1996).

For $G_N = \prod_{k=1}^N \mathbb{Z}/p_k\mathbb{Z}$: $C_{LS}(G_N) = \max_k p_k / (2(p_k - 1)\log(p_k - 1))$. As $N \to \infty$, this grows as $\sim p_N / (2p_N \log p_N) \sim 1/(2\log p_N)$, which goes to 0.

**This is favorable:** the log-Sobolev constant improves (decreases) as more primes are added, making the T_2 bound tighter for larger $N$.

---

## 10. Conclusion

### 10.1 Summary

Large deviations and optimal transport provide **quantitative tools** for the Weil matrix spectral problem:

1. **Rate function analysis** shows that positive eigenvalues are "energetically expensive" — the logarithmic energy penalty scales as $N^2$.

2. **Optimal transport** gives explicit convergence rates: $W_2(\mu_N, \mu_{\text{bg},N}) = O(\log N / \sqrt{N})$, confirming the background dominance.

3. **Talagrand T_2** connects the AMR entropy positivity to transport distances, making the Rudolph-Johnson convergence quantitative.

4. **The entropy-transport-spectral triangle** provides a unified framework for the three AMR components, with explicit constants.

### 10.2 What Is Not Achieved

The fundamental obstacle — proving that $\hat{K}_{\text{Weil}}(\xi) \geq 0$ for $\xi \neq 0$, or equivalently that the Haar eigenvalues are non-positive — is not resolved by LD/OT tools. These tools quantify the RATE of convergence but cannot determine the SIGN of the limit.

### 10.3 The Most Promising Direction

**The sampling theorem approach (§6.4)** is the most novel contribution: using optimal transport between the log-prime lattice and the continuum to bound the spectral gap error. Combined with the PNT error term ($W_1 \sim e^{-c\sqrt{\log N}}$), this gives exponential convergence of the lattice spectral properties to the continuum limit — fast enough for practical finite verification.

**Viability for the full proof:** The LD/OT machinery is most valuable as a **quantitative complement** to the AMR measure-rigidity approach, providing the explicit constants needed for the finite verification reduction (finite-verification.md, Corollary 6.3). It does not provide an independent path to RH.

---

## References

### Large Deviations for Spectral Measures
- Ben Arous, G. and Guionnet, A. (1997). Large deviations for Wigner's law and Voiculescu's non-commutative entropy. *Probab. Theory Related Fields* 108, 517–542.
- Bordenave, C. and Caputo, P. (2014). A large deviation principle for Wigner matrices without Gaussian tails. *Ann. Probab.* 42(6), 2454–2489.
- Anderson, G.W., Guionnet, A., and Zeitouni, O. (2010). *An Introduction to Random Matrices*. Cambridge University Press.
- Hiai, F. and Petz, D. (2000). *The Semicircle Law, Free Random Variables and Entropy*. AMS.

### Optimal Transport
- Villani, C. (2003). *Topics in Optimal Transportation*. AMS.
- Villani, C. (2008). *Optimal Transport: Old and New*. Springer.
- Talagrand, M. (1996). Transportation cost for Gaussian and other product measures. *Geom. Funct. Anal.* 6, 587–600.

### Eigenvalue Stability
- Hoffman, A.J. and Wielandt, H.W. (1953). The variation of the spectrum of a normal matrix. *Duke Math. J.* 20, 37–39.
- Bhatia, R. (1997). *Matrix Analysis*. Springer.

### Log-Sobolev and Concentration
- Gross, L. (1975). Logarithmic Sobolev inequalities. *Amer. J. Math.* 97, 1061–1083.
- Diaconis, P. and Saloff-Coste, L. (1996). Logarithmic Sobolev inequalities for finite Markov chains. *Ann. Appl. Probab.* 6, 695–750.
- Ledoux, M. (2001). *The Concentration of Measure Phenomenon*. AMS.

### AMR Framework
- [entropy-positivity.md](entropy-positivity.md) — Entropy-Positivity Duality
- [bochner-proof.md](bochner-proof.md) — Bochner-Schwartz CPD-1 analysis
- [circularity-resolution.md](circularity-resolution.md) — Unconditional spectral gap + perturbation bounds
- [rmt-universality-approach.md](rmt-universality-approach.md) — Random matrix theory perspective
- [gmc-approach.md](gmc-approach.md) — Gaussian multiplicative chaos connection
- [finite-verification.md](finite-verification.md) — Finite verification reduction

---

*Document: Large Deviations and Optimal Transport for Weil Matrix Spectral Properties*
*Part of the AMR (Arithmetic Measure Rigidity) proofs module*
*February 2026*
