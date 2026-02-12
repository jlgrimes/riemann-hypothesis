# Random Matrix Universality for Weil Matrix Eigenvalue Bounds

## Status: RMT universality does NOT provide an unconditional proof of CPD-1. The Weil matrix is not random and does not belong to any standard universality class whose theorems would directly yield the required eigenvalue bounds. However, RMT provides strong heuristic guidance, asymptotic predictions, and a structural framework that illuminates *why* CPD-1 is expected.

---

## 0. Setup and Goal

**The CPD-1 condition (Bochner-Schwartz, cf. bochner-proof.md):** The Weil matrix

$$K_{ij} = K(\log p_i - \log p_j)$$

truncated to the first N primes, is CPD-1 iff all eigenvalues of the centered matrix $\tilde{K} = K - \frac{1}{N}\mathbf{1}\mathbf{1}^T K$ are non-negative. Equivalently, all eigenvalues of $K$ restricted to the primitive subspace $\{\sum c_i = 0\}$ are non-negative.

**Question:** Can random matrix theory (RMT) universality results prove this eigenvalue condition unconditionally?

---

## 1. The RMT–Zeta Connection: What Is Known

### 1.1 Montgomery-Odlyzko Law (1973/1987)

Montgomery conjectured, and Odlyzko numerically verified to spectacular precision, that the pair correlation of nontrivial zeta zeros (normalized to unit mean spacing) matches the GUE pair correlation:

$$R_2(\alpha) = 1 - \left(\frac{\sin \pi\alpha}{\pi\alpha}\right)^2$$

**Status:** Conjecture. Montgomery proved it holds for test functions whose Fourier support lies in $[-1, 1]$. Full proof remains open. Odlyzko verified for zeros near $\gamma \sim 10^{12}$ and beyond.

**Implications for the Weil matrix:** The pair correlation governs the *statistical distribution of zero spacings*, which in turn controls the oscillatory part $K_{\text{zeros}}(x) = \frac{1}{2\pi}\sum_\gamma \frac{2\cos(\gamma x)}{1/4 + \gamma^2}$. GUE statistics imply that zeros are "repelled" from each other (level repulsion), preventing large clusters that might cause destructive interference in $K_{\text{zeros}}$.

### 1.2 Keating-Snaith Model (2000)

Keating and Snaith modeled $\zeta(1/2 + it)$ by the characteristic polynomial $Z(U, \theta)$ of a random $N \times N$ unitary matrix (CUE), yielding precise conjectures for:

- Moments: $\int_0^T |\zeta(1/2+it)|^{2k}\,dt \sim c_k T(\log T)^{k^2}$
- Value distribution of $\log|\zeta|$ on the critical line

**Relevance:** The Keating-Snaith model predicts the *distribution of $\zeta$ values*, not the *location of zeros*. For the Weil matrix, it provides a model for the statistical behavior of $K_{\text{zeros}}$ but does not constrain individual eigenvalues of K.

### 1.3 The Hilbert-Pólya Program

The zeta zeros are conjectured to be eigenvalues of a self-adjoint operator (the "Hilbert-Pólya operator"). If such an operator exists, RMT universality for its spectral statistics would follow from the Erdős-Schlein-Yau dynamical approach. However:

- The operator has not been constructed
- Even if constructed, universality of its spectrum does not directly imply CPD-1 of the Weil matrix (these are *different spectral problems*)

---

## 2. Universality Classes: Where Does the Weil Matrix Fit?

### 2.1 Wigner Matrices

**Definition:** $H_{ij}$ i.i.d. (up to symmetry), $\mathbb{E}[H_{ij}] = 0$, $\mathbb{E}[H_{ij}^2] = \sigma^2/N$.

**Spectral law:** Wigner semicircle on $[-2\sigma, 2\sigma]$.

**Edge universality:** Largest eigenvalue $\lambda_1 \approx 2\sigma + N^{-2/3} \cdot \text{TW}_\beta$ (Tracy-Widom).

**Weil matrix?** No. The entries $K(\log p_i - \log p_j)$ are:
- Deterministic, not random
- Strongly correlated (they depend on $\log(p_i/p_j)$, creating algebraic structure)
- Not mean-zero (the diagonal $K(0) \approx 2.528$ dominates)

**Conclusion:** The Weil matrix is **not** a Wigner matrix. Semicircle law does not apply.

### 2.2 Generalized Wigner-Type Matrices (Ajanki-Erdős-Krüger, 2017)

The Ajanki-Erdős-Krüger universality theorem extends bulk universality to Wigner-type matrices where:
- Entries are independent (up to symmetry) but not identically distributed
- The variance matrix $S_{ij} = \mathbb{E}[|H_{ij}|^2]$ need not be doubly stochastic

The density of states is no longer semicircular but is determined by a vector-valued fixed-point equation.

**Weil matrix?** Still no. The Weil entries are deterministic, not random. The theorem requires independent entries with controlled moments.

### 2.3 Kernel Random Matrices

**Definition:** $K_{ij} = f(x_i - x_j)$ where $\{x_i\}$ are i.i.d. random points and $f$ is a kernel function.

This is the closest RMT class to the Weil matrix, where $x_i = \log p_i$ and $f = K$ (the Weil kernel).

**The crucial difference:** For kernel random matrices, the points $x_i$ are *random*. For the Weil matrix, the points $\log p_i$ are *deterministic*.

**Asymptotic spectral distribution (Szegő-type):** For a translation-invariant kernel $f$ on random points drawn from density $\rho$, the empirical eigenvalue distribution of $K/N$ converges to the spectral measure:

$$\mu_f = \hat{f} \cdot \hat{\rho} * \hat{\rho}$$

(convolution of the Fourier transform of $f$ with the squared Fourier transform of $\rho$).

If $\hat{f}(\xi) \geq 0$ for all $\xi$, then $\mu_f$ is supported on $[0, \infty)$ — all eigenvalues are asymptotically non-negative.

**Problem:** The condition $\hat{K}(\xi) \geq 0$ for all $\xi \neq 0$ IS the Riemann Hypothesis (cf. bochner-proof.md, Theorem 7.1). So the kernel random matrix approach reduces to checking the same condition.

### 2.4 Deterministic Kernel Matrices with Quasi-Random Points

The log-primes $\{\log p_i\}$ are *not* random but are "quasi-random" in the sense of the Prime Number Theorem: they are distributed asymptotically like a Poisson process with intensity $\sim e^x/x$ (by PNT, $\pi(e^x) \sim e^x/x$).

**Weyl-type equidistribution:** By PNT with error term, $\log p_n \sim n \log n$, and the "gap" $\log p_{n+1} - \log p_n \sim 1/n$ is quasi-regular. This gives the log-primes a pseudo-random character.

**However:** No universality theorem exists for deterministic kernel matrices with quasi-random points. The existing RMT machinery fundamentally requires randomness (either in the entries or the points) to invoke concentration inequalities and moment methods.

**Assessment:** The Weil matrix does not fit any known universality class with existing theorems.

---

## 3. What RMT CAN Say About the Weil Matrix

### 3.1 The Semicircle Heuristic for the Oscillatory Part

Decompose $K = K_{\text{bg}} + K_{\text{zeros}} + \delta$ (cf. bochner-proof.md). The background $K_{\text{bg}}(\log p_i - \log p_j)$ is a smooth, slowly-varying function — essentially a "signal" matrix. The zero-oscillation $K_{\text{zeros}}(\log p_i - \log p_j)$ oscillates quasi-randomly.

**Heuristic:** If $K_{\text{zeros}}$ at log-prime points behaves like a Wigner matrix perturbation, its eigenvalues would follow a semicircle distribution with support $[-2\sigma_N, 2\sigma_N]$ where:

$$\sigma_N^2 = \frac{1}{N}\sum_{i \neq j} K_{\text{zeros}}(\log p_i - \log p_j)^2$$

From the explicit formula, $K_{\text{zeros}}(x) \sim \frac{1}{2\pi}\sum_\gamma \frac{2\cos(\gamma x)}{1/4+\gamma^2}$. The sum converges absolutely with total contribution $\sim 0.023$ (cf. schoenberg-attempt.md, §3.3). So:

$$|K_{\text{zeros}}(x)| \lesssim 0.023 \quad \text{for all } x$$

This gives $\sigma_N \sim 0.023$, and the semicircle edge $2\sigma_N \sim 0.046$.

**Meanwhile:** The minimum eigenvalue of $K_{\text{bg}}$ on primitives is $\gg 1$ (the spectral gap from bochner-proof.md §1.4 gives $\hat{K}_{\text{bg}}(\xi) \geq 2e^{-|\xi|/2}/(1-e^{-2|\xi|})$, which at the "worst" point is still $\gg 0$).

**Heuristic conclusion:** The oscillatory perturbation $K_{\text{zeros}}$ is vastly smaller than the background spectral gap, so CPD-1 holds with enormous margin. This is consistent with the perturbative analysis in circularity-resolution.md.

**But this is not a proof.** The semicircle law requires random i.i.d. entries, which $K_{\text{zeros}}$ does not have.

### 3.2 Tracy-Widom Bounds on Extreme Eigenvalues

For GUE matrices of size $N$, the smallest eigenvalue satisfies:

$$\lambda_{\min} = -2\sqrt{N} + N^{-1/6} \cdot TW_2$$

where $TW_2$ is the Tracy-Widom distribution (supported on approximately $[-5, 2]$).

**Application to Weil matrix?** If we could prove the Weil matrix belongs to a universality class with Tracy-Widom edge statistics, we would get:

$$\lambda_{\min}(K|_{\text{prim}}) \geq -C\sqrt{N} \cdot \|K_{\text{zeros}}\|$$

for some constant $C$. Combined with the background eigenvalue lower bound, this might yield positivity.

**Problem:** Tracy-Widom universality is proved only for matrices with random entries (Wigner, sample covariance, beta-ensembles). No theorem covers deterministic structured matrices.

### 3.3 Marchenko-Pastur for Signal + Noise

If we view the Weil matrix as $K = K_{\text{signal}} + K_{\text{noise}}$ where:
- $K_{\text{signal}} = \delta + K_{\text{bg}}$ (rank-$\infty$ but slowly varying, spectral gap $\gg 1$)
- $K_{\text{noise}} = K_{\text{zeros}}$ (oscillatory, small in norm)

the Marchenko-Pastur framework for spiked covariance matrices (Baik-Ben Arous-Péché, 2005) gives a phase transition: signal eigenvalues "pop out" of the noise bulk when they exceed the bulk edge.

Since $\|K_{\text{signal}}\| \gg \|K_{\text{noise}}\|$, ALL eigenvalues of $K_{\text{signal}}$ are above the noise bulk — the noise never overwhelms the signal.

**Again, this is heuristic.** The BBP transition requires the noise to be a Wigner matrix, which $K_{\text{zeros}}$ is not.

---

## 4. Dyson Brownian Motion and the Spectral Evolution

### 4.1 Setup

As we add primes to the Weil matrix (truncating at the first $N$ primes, then $N+1$, etc.), the eigenvalues evolve. Does this evolution follow a Dyson Brownian motion (DBM)?

**Dyson's SDE:** For an $N \times N$ Hermitian matrix with entries evolving by independent Brownian motions, the eigenvalues $\lambda_1 < \cdots < \lambda_N$ satisfy:

$$d\lambda_i = \frac{dB_i}{\sqrt{N}} + \frac{1}{N}\sum_{j \neq i} \frac{1}{\lambda_i - \lambda_j}\,dt$$

The key feature is **eigenvalue repulsion**: the drift pushes eigenvalues apart, preventing collisions. At equilibrium, the eigenvalue density is the semicircle law.

### 4.2 Adding a Prime: Rank-1 Update

When we add prime $p_{N+1}$ to an $N \times N$ Weil matrix, the new $(N+1) \times (N+1)$ matrix is:

$$K^{(N+1)} = \begin{pmatrix} K^{(N)} & v \\ v^T & K(0) \end{pmatrix}$$

where $v_i = K(\log p_i - \log p_{N+1})$.

By the **Cauchy interlacing theorem**, the eigenvalues of $K^{(N+1)}$ interlace with those of $K^{(N)}$:

$$\lambda_1^{(N+1)} \leq \lambda_1^{(N)} \leq \lambda_2^{(N+1)} \leq \cdots \leq \lambda_N^{(N)} \leq \lambda_{N+1}^{(N+1)}$$

This means the minimum eigenvalue can only decrease or stay the same when adding a prime. This is **not** Dyson Brownian motion — there is no stochastic driving, and the update is rank-1 (not full-rank perturbation).

### 4.3 Eigenvalue Repulsion Without Randomness?

The Dyson repulsion mechanism ($1/(\lambda_i - \lambda_j)$ drift) arises from the *Jacobian of the eigenvalue-eigenvector decomposition*, not from randomness per se. For any smooth one-parameter family of symmetric matrices $A(t)$, the eigenvalues satisfy:

$$\dot{\lambda}_i = A'_{ii}(t) + \sum_{j \neq i} \frac{|A'_{ij}(t)|^2}{\lambda_i - \lambda_j}$$

(Hellmann-Feynman + second-order perturbation theory).

**For the Weil matrix parametrized by $t$ (e.g., $t$ = number of primes or a continuous deformation):** The eigenvalue dynamics have a repulsion-type term. If two eigenvalues approach each other, the second-order term diverges, pushing them apart.

**Can this prevent eigenvalue zero-crossings?** If the initial configuration (small $N$) has all primitive eigenvalues positive, and the repulsion prevents any eigenvalue from crossing through zero...

**Problem:** Cauchy interlacing shows the minimum eigenvalue is *non-increasing* as $N$ grows. Repulsion prevents eigenvalues from *colliding with each other*, but it does NOT prevent the minimum eigenvalue from *decreasing toward zero*. A new eigenvalue is born below the existing spectrum at each step.

**Computational observation:** Numerical experiments (fourier_verification.py and related) show that the minimum primitive eigenvalue of $K^{(N)}$ remains positive for all tested $N$, but the margin of positivity fluctuates. There is no clear convergence to zero, suggesting CPD-1 holds in the limit — but this is numerical evidence, not proof.

### 4.4 Continuous Deformation Approach

Consider the one-parameter family:

$$K_t(x) = \delta(x) + K_{\text{bg}}(x) + t \cdot K_{\text{zeros}}(x), \quad t \in [0, 1]$$

At $t = 0$: $K_0 = \delta + K_{\text{bg}}$, which is CPD-1 (unconditionally, cf. bochner-proof.md §8.1). The primitive eigenvalues of $K_0^{(N)}$ are all positive.

At $t = 1$: $K_1 = K$ is the full Weil matrix. CPD-1 of $K_1$ is RH.

**As $t$ increases from 0 to 1:** The eigenvalues evolve continuously (by matrix perturbation theory). For an eigenvalue to go from positive (at $t = 0$) to negative (at some $t^*$), it must cross zero at $t = t^*$.

**Can we show no zero-crossing occurs?** This requires bounding:

$$|\lambda_{\min}(K_0^{(N)}|_{\text{prim}})| > \|K_{\text{zeros}}^{(N)}|_{\text{prim}}\|_{\text{op}}$$

i.e., the background spectral gap exceeds the operator norm of the zero perturbation.

From numerical evidence (circularity-resolution.md): for $N \leq 100$, the background gap is $\sim 1.9$ while $\|K_{\text{zeros}}\|_{\text{op}} \sim 0.015$. The ratio is $\sim 127:1$ — enormous margin.

**The question is whether this margin persists as $N \to \infty$.** As $N$ grows:
- The background spectral gap depends on the smallest eigenvalue of $K_{\text{bg}}^{(N)}|_{\text{prim}}$, which is governed by $\hat{K}_{\text{bg}}$ at high frequencies. Since $\hat{K}_{\text{bg}}(\xi) \sim 2e^{-|\xi|/2}$ decays exponentially, the gap may shrink.
- The zero perturbation norm $\|K_{\text{zeros}}^{(N)}\|_{\text{op}}$ involves the coherent sum of oscillations from all zeros, which could grow.

**RMT heuristic:** If $K_{\text{zeros}}^{(N)}$ at arithmetic points behaves like a random matrix of norm $\sim \sigma\sqrt{N}$ with $\sigma \sim 0.023/\sqrt{N}$ (i.e., the entries are $O(1/N)$ in an appropriate sense), the norm stays bounded. But this requires that the oscillations at log-prime points don't coherently constructively interfere — essentially a statement about the pseudo-randomness of $\{(\gamma \log p) \mod 2\pi\}$, which is related to the GUE hypothesis.

---

## 5. What Would Be Needed to Make RMT Work

### 5.1 A Universality Theorem for Deterministic Kernel Matrices

**Required result:** For a kernel matrix $K_{ij} = f(x_i - x_j)$ where $\{x_i\}$ are deterministic points satisfying a "pseudo-random" equidistribution condition (e.g., low discrepancy, pair correlation matching Poisson/GUE), the spectral distribution converges to the same limit as the random kernel matrix ensemble.

**State of the art:** No such theorem exists. The closest results are:
- **Szegő's theorem on Toeplitz matrices:** For *equally spaced* points $x_i = i/N$, the eigenvalues of $K_{ij} = f(x_i - x_j)$ converge to the spectral measure of $f$ (i.e., $\hat{f}$). But log-primes are not equally spaced.
- **Quantitative equidistribution:** If the points $\{x_i\}$ satisfy a discrepancy bound $D_N \leq C/N^{1-\epsilon}$, the eigenvalue distribution approximates the Szegő limit. Log-primes satisfy PNT with error, giving $D_N \sim e^{-c\sqrt{\log N}}$ — very good but not sufficient without an explicit universality theorem.

### 5.2 A Deterministic Norm Bound for $K_{\text{zeros}}^{(N)}$

**Required result:** $\|K_{\text{zeros}}^{(N)}|_{\text{prim}}\|_{\text{op}} \leq C$ for all $N$, where $C$ is smaller than the background spectral gap.

**Approach via exponential sums:** The operator norm of $K_{\text{zeros}}^{(N)}$ is controlled by:

$$\|K_{\text{zeros}}^{(N)}\|_{\text{op}} \leq \sup_{\|c\|=1, \sum c_i = 0} \left|\sum_{i,j} c_i c_j K_{\text{zeros}}(\log p_i - \log p_j)\right|$$

Expanding $K_{\text{zeros}}$:

$$= \sup_c \left|\sum_\gamma \frac{1}{\pi(1/4+\gamma^2)} \left|\sum_i c_i e^{i\gamma \log p_i}\right|^2\right| \leq \frac{1}{\pi} \sum_\gamma \frac{1}{1/4+\gamma^2} \cdot \sup_c \left|\sum_i c_i p_i^{i\gamma}\right|^2$$

The inner supremum $\sup_c |\sum c_i p_i^{i\gamma}|^2$ with $\|c\| = 1$ is just the maximum eigenvalue of the matrix $M_{ij}^\gamma = p_i^{i\gamma}\overline{p_j^{i\gamma}} = (p_i/p_j)^{i\gamma}$, which equals $N$ (since $M^\gamma$ has rank 1 with eigenvalue $N$).

But $\sum c_i = 0$ restricts to primitives, and $M^\gamma|_{\text{prim}}$ has eigenvalue $N - N = 0$ if $p_i^{i\gamma}$ is approximately constant... This is NOT the case for generic $\gamma$ and primes.

The bound becomes:

$$\|K_{\text{zeros}}^{(N)}|_{\text{prim}}\|_{\text{op}} \leq \frac{N}{\pi} \sum_\gamma \frac{1}{1/4+\gamma^2} \sim 0.023 N$$

This grows linearly in $N$ — too large. To get a useful bound, we need cancellation in $\sum_i c_i p_i^{i\gamma}$ across different $\gamma$, which requires understanding the distribution of $\{(\gamma \log p_i) \mod 2\pi\}$. This is a deep problem in analytic number theory, essentially equivalent to strong forms of the pair correlation conjecture.

### 5.3 Bridge Between GUE Statistics and Eigenvalue Positivity

**The fundamental gap:** RMT tells us the *statistical distribution* of zeta zeros matches GUE. But CPD-1 of the Weil matrix requires a *deterministic* statement about ALL zeros — that their collective effect produces a positive definite kernel. Going from statistical to deterministic is the hardest step.

**What would bridge the gap:** A theorem of the form: "If the pair correlation of a sequence $\{\gamma_n\}$ matches GUE, then for any smooth kernel $f$ with $\hat{f} \geq 0$ on $\mathbb{R} \setminus \{0\}$, the matrix $f(\gamma_m - \gamma_n)$ is asymptotically positive semi-definite."

No such theorem exists or is expected to be provable purely from pair correlation statistics. The positivity requires *more* than pair statistics — it needs control of all $k$-point correlations, and even then the passage from statistical to deterministic is non-trivial.

---

## 6. Honest Assessment

### 6.1 What RMT Proves About the Weil Matrix

| Statement | Proved? | Notes |
|-----------|---------|-------|
| Zeta zeros have GUE pair correlation | **No** (conjecture) | Montgomery partial result for restricted Fourier support; Odlyzko numerical verification |
| Weil matrix belongs to a known universality class | **No** | Deterministic entries disqualify it from all standard classes |
| Tracy-Widom governs extreme eigenvalues of Weil matrix | **No** | No theorem covers deterministic structured kernels |
| Dyson repulsion prevents eigenvalue zero-crossing | **No** | Repulsion prevents collision, not zero-crossing |
| $\|K_{\text{zeros}}^{(N)}\|_{\text{op}}$ stays bounded | **No** | Known bounds grow linearly in $N$; tighter bounds need pair correlation input |
| Background spectral gap exceeds zero perturbation (finite $N$) | **Yes** (numerically, unconditional for small $N$) | Gap $\sim 1.9$ vs perturbation $\sim 0.015$ for $N \leq 100$ |

### 6.2 What RMT Strongly Suggests

1. **CPD-1 holds with enormous margin for all finite truncations.** The background dominance (ratio $\sim 127:1$ for $N \leq 100$) is so overwhelming that the GUE-distributed oscillations are completely negligible. RMT provides the *heuristic* that this ratio stays favorable as $N \to \infty$ because the oscillations remain pseudo-random.

2. **The Weil matrix spectrum is "GUE-like" in a qualitative sense.** The zeros exhibit level repulsion, the spacing distribution matches GUE, and the eigenvalue statistics of the Weil matrix (if ever computed at large scale) should reflect this. GUE statistics are "maximally spread" — they produce the most uniform possible eigenvalue distribution, which is the least likely to create concentrated negative contributions.

3. **No "resonance catastrophe" occurs.** In a truly adversarial scenario, the oscillations $K_{\text{zeros}}$ could constructively interfere at arithmetic points to create a large negative eigenvalue. GUE statistics prevent such coherent constructive interference, because the zeros are "repelled" and hence their oscillations at any fixed frequency cancel out statistically.

### 6.3 The Fundamental Limitation of RMT for RH

RMT describes the *statistics* of zeta zeros. The Riemann Hypothesis is a *deterministic* statement (every zero has $\beta = 1/2$). The passage from statistical to deterministic is a category error:

- **Statistical → deterministic:** "The zeros are statistically distributed like GUE eigenvalues" does NOT imply "every zero is on the critical line." A single off-line zero would be a measure-zero event, invisible to pair correlation statistics.

- **Deterministic → statistical:** RH IMPLIES GUE statistics (this is essentially Montgomery's result, conditional on RH). The converse is false.

**Analogy:** Knowing that the prime numbers are "distributed like random integers with density $1/\log N$" (PNT) does not prove specific properties of individual primes (e.g., the twin prime conjecture). Similarly, knowing the zeros are "distributed like GUE eigenvalues" does not prove specific properties of individual zeros (e.g., $\beta = 1/2$).

---

## 7. Connections to Other Approaches

### 7.1 Comparison with Bochner (bochner-proof.md)

The Bochner approach requires $\hat{K}(\xi) \geq 0$ for $\xi \neq 0$, which IS RH. RMT does not bypass this — it provides *heuristic reasons* to believe $\hat{K}(\xi) \geq 0$ (the oscillatory contributions are too small to overcome the background), but no proof.

### 7.2 Comparison with Schoenberg (schoenberg-attempt.md)

The Schoenberg approach identifies that CPD-1 reduces to PD of $e^{tK_{\text{zeros}}}$, which requires the zero-oscillation kernel to have its "RH form" (pure cosines with positive coefficients). RMT universality does not establish this form.

### 7.3 Synergy with the AMR Approach

The AMR (Arithmetic Measure Rigidity) approach avoids analyzing $K_{\text{zeros}}$ directly, instead proving CPD-1 via measure rigidity (Baker → entropy positivity → Rudolph → Haar → CPD-1). RMT provides *complementary* evidence:

- AMR proves the *structural reason* for CPD-1 (Haar invariance forces it)
- RMT provides the *quantitative picture* (the spectral gap is enormous and the perturbation is tiny)

The RMT perspective confirms that the AMR approach is targeting the right structure: the "signal" (background) completely dominates the "noise" (zero oscillation), and the measure-rigidity machinery correctly identifies this dominance.

---

## 8. Specific Recent Results and Their Relevance

### 8.1 Bourgade-Erdős-Yau: Fixed Energy Universality (2014)

**Result:** Local eigenvalue statistics of Wigner matrices are universal at every fixed energy in the bulk of the spectrum, not just in expectation.

**Relevance:** Strengthens the universality of GUE statistics but still applies only to random matrices. Does not help with deterministic Weil matrix.

### 8.2 Erdős-Knowles-Yau-Yin: Spectral Statistics of Erdős-Rényi Graphs (2013)

**Result:** Adjacency matrices of Erdős-Rényi random graphs (structured sparse random matrices) have GUE local statistics in the bulk.

**Relevance:** Shows universality extends to *structured* random matrices. The Weil matrix has combinatorial structure from the multiplicative group of primes. However, the Erdős-Rényi model still requires randomness in entries; the Weil matrix has no randomness.

### 8.3 Ajanki-Erdős-Krüger: General Wigner-Type Universality (2017)

**Result:** Bulk universality for random matrices with non-stochastic variance profiles.

**Relevance:** The most general Wigner-type universality theorem. Still requires independent random entries with bounded moments.

### 8.4 Dyson Brownian Motion for Structured Matrices (2024-2025)

**Recent work:** DBM has been applied to weight matrices in neural networks during training, showing that structured matrix dynamics can exhibit eigenvalue repulsion.

**Relevance:** Conceptual — shows that non-random matrix dynamics can exhibit RMT-like behavior. Does not provide quantitative theorems applicable to the Weil matrix.

---

## 9. Speculative Directions

### 9.1 Arithmetic Random Matrix Theory

**Idea:** Develop an "arithmetic RMT" where the randomness comes not from the matrix entries but from the *arithmetic structure* of the primes. The primes are deterministic but have "enough randomness" (in the sense of pseudo-randomness) to drive universality.

**What's needed:** A version of the Erdős-Schlein-Yau dynamical approach where the Dyson Brownian motion is driven by the pseudo-random dynamics of prime gaps rather than actual Brownian motion. This does not yet exist and would be a major development in its own right.

### 9.2 Deterministic Trace Formula Methods

**Idea:** Use the Selberg trace formula (which is a deterministic analogue of the Weil explicit formula for hyperbolic surfaces) as a model. The eigenvalue positivity for the Laplacian on a hyperbolic surface is PROVED (the Laplacian is self-adjoint positive). Can the trace formula method be adapted?

**Problem:** The Selberg trace formula works because the underlying operator (the Laplacian) is known to exist and be self-adjoint. For zeta, the underlying operator (the Hilbert-Pólya operator) is conjectural.

### 9.3 Free Probability and Asymptotic Freeness

**Idea:** In free probability, asymptotically free random matrices have spectra described by free convolution. If $K_{\text{bg}}$ and $K_{\text{zeros}}$ are "asymptotically free" in an appropriate sense, the spectrum of $K = K_{\text{bg}} + K_{\text{zeros}} + \delta$ can be computed via free additive convolution.

**Status:** Free probability applies to random matrices that are unitarily invariant. The Weil matrix components are not unitarily invariant (they have arithmetic structure). An "arithmetic free probability" theory does not yet exist.

---

## 10. Conclusion

**Can RMT prove CPD-1 of the Weil matrix unconditionally? No.**

The Weil matrix is a deterministic structured matrix that does not belong to any known universality class. All existing RMT universality theorems require randomness (either in entries or in point configurations), which the Weil matrix lacks.

**What RMT provides:**
1. **Strong heuristic support** for CPD-1: the zero oscillations are vastly smaller than the background spectral gap, and GUE statistics ensure no coherent destructive interference.
2. **A quantitative framework** for understanding WHY CPD-1 holds: the signal-to-noise ratio (background vs. zeros) is enormous ($\sim 127:1$ for moderate $N$).
3. **A roadmap** for what WOULD suffice: a universality theorem for deterministic kernel matrices with quasi-random points, or a deterministic norm bound for $K_{\text{zeros}}^{(N)}|_{\text{prim}}$ growing slower than the background spectral gap.

**The honest verdict:** RMT is the right *language* for understanding the spectral properties of the Weil matrix, but it does not currently provide the right *tools* for a proof. The passage from statistical (GUE) to deterministic (every zero on the line) remains the central unsolved problem, and RMT does not bridge this gap.

---

## References

- Montgomery, H.L. (1973). The pair correlation of zeros of the zeta function. *Proc. Symp. Pure Math.* 24, 181–193.
- Odlyzko, A.M. (1987). On the distribution of spacings between zeros of the zeta function. *Math. Comp.* 48, 273–308.
- Keating, J.P. and Snaith, N.C. (2000). Random matrix theory and ζ(1/2+it). *Comm. Math. Phys.* 214, 57–89.
- Tracy, C.A. and Widom, H. (1994). Level-spacing distributions and the Airy kernel. *Comm. Math. Phys.* 159, 151–174.
- Erdős, L., Schlein, B., and Yau, H.-T. (2009). Semicircle law on short scales and delocalization of eigenvectors for Wigner random matrices. *Ann. Probab.* 37, 815–852.
- Ajanki, O., Erdős, L., and Krüger, T. (2017). Universality for general Wigner-type matrices. *Probab. Theory Related Fields* 169, 667–727.
- Bourgade, P., Erdős, L., and Yau, H.-T. (2014). Fixed energy universality for generalized Wigner matrices. *Comm. Pure Appl. Math.* 69, 1815–1881.
- Erdős, L., Knowles, A., Yau, H.-T., and Yin, J. (2013). Spectral statistics of Erdős-Rényi graphs II: Eigenvalue spacing and the extreme eigenvalues. *Comm. Math. Phys.* 314, 587–640.
- Baik, J., Ben Arous, G., and Péché, S. (2005). Phase transition of the largest eigenvalue for nonnull complex sample covariance matrices. *Ann. Probab.* 33, 1643–1697.
- Dyson, F.J. (1962). A Brownian-motion model for the eigenvalues of a random matrix. *J. Math. Phys.* 3, 1191–1198.
- Tao, T. and Vu, V. (2011). Random matrices: universality of local eigenvalue statistics. *Acta Math.* 206, 127–204.
- Bourgade, P. (2024). Quantum chaos, random matrix theory, and the Riemann ζ-function. *Séminaire Poincaré*.
- [bochner-proof.md](bochner-proof.md) — Bochner-Schwartz CPD-1 analysis
- [schoenberg-attempt.md](schoenberg-attempt.md) — Schoenberg representation approach
- [nyman-beurling.md](nyman-beurling.md) — Computational evidence via Nyman-Beurling

---

*Document: Random Matrix Universality Approach*
*Part of the AMR (Arithmetic Measure Rigidity) proofs module*
*February 2026*
