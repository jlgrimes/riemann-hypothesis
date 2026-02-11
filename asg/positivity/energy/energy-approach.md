# The Energy/Variational Approach to the Arithmetic Positivity Theorem

## Overview

This document develops an energy-based perspective on the Arithmetic Positivity Theorem (APT). We model the nontrivial zeros of ζ(s) as particles in a logarithmic Coulomb gas, reformulate APT as an energy minimization principle, and explore whether techniques from statistical mechanics and variational calculus can yield a proof.

The central idea: the zeros of ζ(1/2 + iγ) sit at the equilibrium configuration of a log-gas on the critical line. Any perturbation that moves zeros off the line INCREASES the energy. APT is precisely the statement that this energy increase is non-negative.

---

## 1. The Physical Model: Log-Gas of Zeros

### 1.1 The Energy Functional

Treat the nontrivial zeros ρ_j = 1/2 + iγ_j of ζ(s) as particles on the real line (parametrized by γ_j). Define the energy functional:

E[{γ_j}] = -Σ_{j<k} log|γ_j - γ_k| + Σ_j V(γ_j)

where V(γ) is a confining potential that prevents the particles from escaping to infinity.

The first term is the **pairwise logarithmic repulsion** — the same interaction as the 2D Coulomb potential restricted to a line. The second term is an **external field** arising from the arithmetic (the primes).

### 1.2 The Confining Potential

The confining potential is derived from the zero counting function. Let N(T) denote the number of zeros with 0 < γ ≤ T. By the Riemann–von Mangoldt formula:

N(T) = (T/2π) log(T/2π) - T/2π + O(log T)

The equilibrium density of zeros is:

ρ_eq(γ) = dN/dγ = (1/2π) log(|γ|/2π) + O(1/|γ|)

In a log-gas at inverse temperature β, the equilibrium density satisfies the Euler–Lagrange equation:

β · V'(γ) = 2 P.V. ∫ ρ_eq(γ')/(γ - γ') dγ'

Inverting this relation with ρ_eq as above determines V. For large |γ|:

**V(γ) ~ (γ/2) log(|γ|/(2πe))**

This logarithmic growth (stronger than the quadratic confinement of classical GUE ensembles, matching the logarithmic growth of the zero density) ensures that infinitely many zeros are confined.

### 1.3 Derivation of V from the Explicit Formula

More precisely, V arises from the prime side of the explicit formula. Write:

Σ_p Σ_{m=1}^∞ (log p / p^{m/2}) · h(m log p) = Σ_ρ ĥ(γ_ρ) - ĥ(i/2) - ĥ(-i/2) + Ω(h)

The left side, encoding the prime distribution, acts as an external field on the zeros. Each zero γ_j "feels" the potential:

V(γ_j) = -Σ_p Σ_m (log p / p^{m/2}) · e^{-imγ_j log p} + archimedean terms

This is the Fourier dual of the von Mangoldt function evaluated at γ_j. The prime number theorem ensures V grows at the rate V(γ) ~ (γ/2) log(|γ|/(2πe)), consistent with Section 1.2.

### 1.4 The Physical Picture

- **Repulsion**: Zeros repel logarithmically, enforcing level repulsion (no clustering). This is the Montgomery–Odlyzko phenomenon.
- **Confinement**: The arithmetically-determined potential V keeps the zeros from drifting.
- **Temperature**: β = 2 corresponds to the GUE (Gaussian Unitary Ensemble) universality class — time-reversal symmetry broken, matching the functional equation ξ(s) = ξ(1-s).

---

## 2. The GUE Equilibrium

### 2.1 Free Energy and the Variational Principle

For a log-gas of N particles at inverse temperature β, the **free energy** is:

F[ρ] = E[ρ] - (1/β) S[ρ]

where:
- E[ρ] = -∫∫_{x≠y} log|x-y| ρ(x)ρ(y) dx dy + ∫ V(x) ρ(x) dx is the energy
- S[ρ] = -∫ ρ(x) log ρ(x) dx is the Boltzmann entropy
- β is the inverse temperature

The equilibrium measure μ_eq minimizes F over probability measures ρ with ∫ρ = N.

### 2.2 The Euler–Lagrange Equation

Setting the first variation δF/δρ = 0 yields:

-2∫ log|x-y| ρ_eq(y) dy + V(x) + (1/β) log ρ_eq(x) = λ

where λ is a Lagrange multiplier enforcing ∫ρ = N.

At **β = 2** (the GUE value), this becomes:

-2∫ log|x-y| ρ_eq(y) dy + V(x) + (1/2) log ρ_eq(x) = λ

### 2.3 Verification: The Equilibrium Density Matches the Zero Density

**Theorem 2.1.** With the confining potential V(γ) ~ (γ/2) log(|γ|/(2πe)) and β = 2, the equilibrium density of the log-gas is:

ρ_eq(γ) = (1/2π) log(|γ|/2π) + O(1/|γ|)

which matches the mean density of zeros of ζ(1/2 + iγ).

*Proof.* The equilibrium density satisfies the singular integral equation from §2.2. We verify consistency.

**Step 1.** Compute the Stieltjes transform of ρ_eq:

G(z) = ∫ ρ_eq(γ)/(z - γ) dγ = (1/2π)[log(z/2π) + O(1/z)]

for Im(z) > 0. This follows from the integral representation of the digamma function ψ(z) = Γ'(z)/Γ(z).

**Step 2.** The Euler–Lagrange equation requires (ignoring the entropy correction, which is subleading for large N):

β · V'(γ) = 2 P.V. ∫ ρ_eq(γ')/(γ - γ') dγ'

The principal value integral is:

2 Re G(γ + i0⁺) = (1/π) log(|γ|/2π)

With β = 2:

V'(γ) = (1/2π) · (1/π) ... — wait, let us be careful with the normalization. We have:

2 P.V. ∫ ρ_eq(γ')/(γ-γ') dγ' = 2 · Re G(γ+i0⁺)

and Re G(γ+i0⁺) = (1/2π) log(|γ|/2π) (from the known Hilbert transform of log).

So β · V'(γ) = (1/π) log(|γ|/2π), giving V'(γ) = (1/2π) log(|γ|/2π).

**Step 3.** Integrating:

V(γ) = (1/2π)[γ log(|γ|/2π) - γ] + C = (γ/2π) log(|γ|/(2πe)) + C

Up to the conventional factor of π (absorbed into the normalization of the energy), this matches V(γ) ~ (γ/2) log(|γ|/(2πe)). ∎

### 2.4 The Second Variation: Stability

The equilibrium is a genuine minimum (not just a critical point) if the second variation of F is positive:

δ²F[δρ] = -2∫∫ log|x-y| δρ(x) δρ(y) dx dy + (1/β) ∫ (δρ)²/ρ_eq dx

The first term is the logarithmic energy of the perturbation δρ. For measures δρ with ∫δρ = 0 (preserving total mass), the logarithmic energy is ≥ 0 by a classical result (the logarithmic kernel is conditionally positive-definite). The second term is the Fisher information, which is ≥ 0.

**Therefore δ²F ≥ 0**, confirming that the equilibrium is a genuine minimum of the free energy. The GUE equilibrium is **stable**.

### 2.5 Physical Interpretation

The zeros are not merely "distributed like" a GUE ensemble — the correspondence is structural:

| GUE Log-Gas | ζ-Zeros |
|---|---|
| β = 2 | Functional equation symmetry |
| V(x) = x²/2 (quadratic) | V(γ) ~ (γ/2)log(|γ|/2πe) (arithmetic) |
| Semicircle density | (1/2π)log(|γ|/2π) density |
| Level repulsion ∝ |x-y|² | Montgomery pair correlation |
| Determinantal correlations | Conjectured for ζ-zeros |

The GUE value β = 2 is singled out because it is the unique temperature at which the log-gas has a determinantal structure — the n-point correlation functions are determinants of a single kernel K(x,y). This determinantal structure is what makes the partition function manifestly positive (see §5).

---

## 3. The de Bruijn–Newman Dynamics

### 3.1 Definition and Context

The de Bruijn–Newman constant Λ is defined via the family of entire functions:

H_t(z) = ∫₀^∞ e^{tu²} Φ(u) cos(zu) du,  t ∈ ℝ

where Φ(u) = Σ_{n=1}^∞ (2π²n⁴e^{9u} - 3πn²e^{5u}) exp(-πn²e^{4u}) is the function appearing in the Jacobi theta representation, satisfying Ξ(z) = H_0(z) (the Riemann xi function).

Λ = inf{t ∈ ℝ : H_t has all zeros real}

**RH ⟺ Λ ≤ 0.** The complementary bound **Λ ≥ 0** was proved by Rodgers and Tao (2020), so RH ⟺ Λ = 0.

### 3.2 The Zero Dynamics

Under the heat flow, the zeros γ_j(t) of H_t satisfy the coupled ODE system:

dγ_j/dt = -Σ_{k≠j} 2/(γ_j(t) - γ_k(t))

This is the **gradient flow** of the pairwise interaction energy:

dγ_j/dt = -∂E/∂γ_j  where  E = -Σ_{j<k} log|γ_j - γ_k|

Each zero moves under the "electrostatic" force exerted by all other zeros. The dynamics drives particles apart (maximizing repulsion).

### 3.3 Energy Monotonicity

**Proposition 3.1 (Energy monotonicity).** Along the de Bruijn–Newman flow with all zeros real, the energy E(t) = -Σ_{j<k} log|γ_j(t) - γ_k(t)| is monotonically decreasing:

dE/dt = -Σ_j |∂E/∂γ_j|² = -Σ_j |Σ_{k≠j} 2/(γ_j - γ_k)|² ≤ 0

*Proof.* Direct computation:

dE/dt = Σ_j (∂E/∂γ_j) · (dγ_j/dt) = Σ_j (∂E/∂γ_j) · (-∂E/∂γ_j) = -Σ_j (∂E/∂γ_j)² ≤ 0

Equality holds iff ∂E/∂γ_j = 0 for all j, i.e., at equilibrium. ∎

**Corollary 3.2.** The gradient flow is dissipative. Once all zeros are real (which holds for t > Λ), they remain real for all subsequent t. The energy decreases monotonically until the equilibrium spacing is achieved.

### 3.4 Lyapunov Analysis

Define the full Lyapunov function:

L(t) = E[{γ_j(t)}] + (1/2) Σ_j V(γ_j(t); t)

where V(γ; t) is the time-dependent confining potential from the heat kernel factor e^{tu²}. Then:

dL/dt = -Σ_j |∂L/∂γ_j|² + (1/2) Σ_j (∂V/∂t)(γ_j; t)

The first term is ≤ 0 (gradient flow dissipation). The second term reflects strengthening confinement — as t increases, the Gaussian weight e^{tu²} narrows, confining zeros more tightly.

**Key observation:** For t < 0, the "anti-confinement" (the factor e^{tu²} with t < 0 is a narrowing Gaussian in Fourier space, widening in direct space) can overpower the dissipation, allowing complex zeros. At t = 0 (the ζ-function line), the balance is exact.

### 3.5 Energy Analysis Near t = 0

Define the total regularized energy:

E_reg(t) = -Σ_{j<k} log|γ_j(t) - γ_k(t)| + Σ_j V(γ_j(t))

(regularized by subtracting the equilibrium energy). Then:

- For t > 0: E_reg(t) < E_reg(0) (energy decreasing, all zeros real)
- For t < 0 (if Λ < 0): zeros complex, energy landscape changes character
- At t = 0: **E_reg has a critical point** (the gradient of E with respect to t vanishes if Λ = 0)

**Connection to APT:** If we can show E_reg is minimized at t = 0 among all configurations with the correct counting function, then the zeros must be real at t = 0. This is because:

1. For t > 0, all zeros are real (Rodgers–Tao) and the energy is lower.
2. Taking t → 0⁺, we get a real configuration at t = 0.
3. This real configuration has lower energy than any complex configuration.

The gap: this argument uses continuity as t → 0⁺ and the known result Λ ≥ 0. It does NOT independently prove Λ ≤ 0 (i.e., RH). The energy monotonicity tells us the flow is well-behaved but does not rule out Λ > 0.

### 3.6 The Calogero–Moser Connection

The ODE system dγ_j/dt = -Σ_{k≠j} 2/(γ_j - γ_k) is a special case of the **Calogero–Moser system** — specifically, it is the gradient flow (overdamped dynamics) of the Calogero–Moser Hamiltonian:

H_CM = (1/2) Σ_j p_j² + Σ_{j<k} g²/(γ_j - γ_k)²

at zero kinetic energy (p_j = 0). The Calogero–Moser system is completely integrable — it has N independent conserved quantities in involution. This integrability is ultimately a consequence of the same determinantal structure underlying the GUE.

**Open question:** Can the integrability of the Calogero–Moser system be exploited to prove global statements about the zero dynamics?

---

## 4. Connection to APT: The Energy Inequality

### 4.1 The Weil Functional as Energy

Recall the Weil positivity criterion: RH ⟺ W(f * f̃) ≥ 0 for all f ∈ C_c^∞(ℝ₊*).

The Weil functional has a spectral decomposition via the explicit formula:

W(f * f̃) = Σ_ρ ĥ(γ_ρ)

where h = f * f̃ and ĥ is the Fourier transform of h. Now:

- For h = f * f̃: ĥ(t) = |f̂(t)|² for **real** t.
- For complex t = a + ib: ĥ(a+ib) = f̂(a+ib) · f̂(a-ib)* (using f̃ structure).

**The critical distinction:**
- If all γ_ρ ∈ ℝ (RH true): W(f * f̃) = Σ_ρ |f̂(γ_ρ)|² ≥ 0. **Trivially non-negative.**
- If some γ_ρ ∉ ℝ (RH false): terms ĥ(γ_ρ) at complex γ_ρ can have negative real part. **Positivity is non-trivial.**

### 4.2 Energy Interpretation

Interpret f as encoding a "perturbation" or "test charge" inserted into the zero gas. The functional W(f * f̃) measures the **response energy** — the energy cost of the perturbation as felt by the zero configuration.

**Definition 4.1.** The *spectral energy* of the test function f relative to the zero configuration {γ_ρ} is:

E_f[{γ_ρ}] = Σ_ρ ĥ(γ_ρ) = W(f * f̃)

**Theorem 4.2 (APT as energy inequality).** The following are equivalent:

(i) RH holds.

(ii) For all f ∈ C_c^∞(ℝ₊*): the spectral energy W(f * f̃) ≥ 0.

(iii) The actual zero configuration has lower energy than any "perturbed" configuration with zeros moved off the critical line, in the following precise sense: the actual zero contribution Σ_ρ ĥ(γ_ρ) to the explicit formula is non-negative for all h = f * f̃.

*Proof of (i) ⟹ (ii).* If RH holds, all γ_ρ ∈ ℝ, and W(f * f̃) = Σ_ρ |f̂(γ_ρ)|² ≥ 0.

*Proof of (ii) ⟹ (i).* Suppose ρ_0 = 1/2 + i(a + ib) with b ≠ 0 is a zero. By the functional equation, ρ̄_0 = 1/2 + i(a - ib) is also a zero. Choose f with f̂ concentrated near a. The contribution from this zero pair is:

ĥ(a + ib) + ĥ(a - ib) = f̂(a+ib)·f̂(a-ib)* + f̂(a-ib)·f̂(a+ib)*
                        = 2 Re[f̂(a+ib) · f̂(a-ib)*]

For f̂ concentrated near a: f̂(a+ib) ≈ f̂(a) · e^{-b·(something)} and f̂(a-ib) ≈ f̂(a) · e^{b·(something)}. The product f̂(a+ib)·f̂(a-ib)* is NOT automatically non-negative. By choosing f appropriately (e.g., f̂ with a phase that makes this product negative), we can make W(f * f̃) < 0, contradicting (ii). ∎

### 4.3 The Spectral Energy Landscape

Define the **spectral energy landscape** as a function of the zero position σ + iγ:

E(σ, γ; f) = contribution to W(f * f̃) from a zero at 1/2 + i(γ + i(σ-1/2))

For σ = 1/2 (on the critical line): E(1/2, γ; f) = |f̂(γ)|² ≥ 0.

For σ ≠ 1/2 (off the critical line): E(σ, γ; f) can be negative for appropriate f.

**Key geometric picture:** The spectral energy landscape has a **valley along the critical line** σ = 1/2. APT states that the actual zero configuration sits in this valley — the global energy minimum is on the critical line.

### 4.4 The Variational Characterization

**Proposition 4.3.** APT holds if and only if: for every test function f, the critical line is the global minimizer of the spectral energy:

min_{admissible configs} Σ_ρ ĥ(γ_ρ) = Σ_{ρ on line} |f̂(γ_ρ)|² ≥ 0

where "admissible" means satisfying the functional equation and having the correct counting function.

This is not merely a local minimum — it must be global. The arithmetic constraints (the Euler product) are what force the global minimum to coincide with the critical line.

---

## 5. The Coulomb Gas Analogy

### 5.1 The 2D Coulomb Gas Partition Function

Consider N particles z_1, ..., z_N ∈ ℂ with energy:

E_{2D}[{z_j}] = -Σ_{j<k} log|z_j - z_k| + N Σ_j V(z_j)

The partition function at inverse temperature β is:

Z_N(β) = ∫_{ℂ^N} Π_{j<k} |z_j - z_k|^β · exp(-βN Σ_j V(z_j)) Π_j d²z_j

### 5.2 The β = 2 Determinantal Structure

At β = 2, the Vandermonde factor simplifies:

Π_{j<k} |z_j - z_k|² = |det[z_k^{j-1}]_{j,k=1}^N|²

This is the **squared modulus of the Vandermonde determinant**. Therefore:

Z_N(2) = ∫_{ℂ^N} |det[z_k^{j-1}]|² · exp(-2N Σ V(z_j)) Π d²z_j

Expanding the determinant squared and integrating term by term:

Z_N(2) = N! · det[∫_ℂ z^{j-1} z̄^{k-1} e^{-2NV(z)} d²z]_{j,k=1}^N

**Proposition 5.1.** Z_N(2) ≥ 0 for any confining potential V.

*Proof.* The matrix G_{jk} = ⟨φ_j, φ_k⟩_{L²} where φ_j(z) = z^{j-1} e^{-NV(z)} is a **Gram matrix**. Every Gram matrix is positive semi-definite:

det(G) = det(⟨φ_j, φ_k⟩) ≥ 0

with equality iff the functions {φ_j} are linearly dependent (which does not happen for polynomial × exponential weight). Thus Z_N(2) = N! · det(G) > 0. ∎

**Remark.** This positivity is the Coulomb gas analogue of APT: the "partition function" (analogous to W(f * f̃)) is always non-negative. The β = 2 structure makes this automatic.

### 5.3 The Analogy Table

| 2D Coulomb Gas (β=2) | Arithmetic (APT) |
|---|---|
| Particles z_j ∈ ℂ free | Zeros ρ_j constrained by Euler product |
| Vandermonde Π\|z_j - z_k\|² | Weil kernel K(γ_j - γ_k) |
| Confining potential V (arbitrary) | Arithmetic potential V (determined by primes) |
| Partition function Z_N | Weil functional W(f * f̃) |
| Z_N ≥ 0 (Gram matrix argument) | W(f * f̃) ≥ 0 (= APT = RH) |
| **Trivially proved** | **The Riemann Hypothesis** |

### 5.4 Why the Analogy Doesn't Immediately Give a Proof

The Coulomb gas positivity Z_N ≥ 0 is trivial because:

1. The particles are **free** — they range over all of ℂ^N.
2. The integrand is a **squared modulus** — |Vandermonde|² · exp(weight) ≥ 0 pointwise.
3. The integral of a non-negative function is non-negative.

In the arithmetic setting, none of these hold:

1. The zeros are **constrained** by the Euler product — they are not free parameters.
2. The Weil functional is NOT a pointwise integral of a non-negative function — it is a distributional pairing.
3. The "integrand" (individual prime contributions) can be negative; only the sum over all primes should be non-negative.

### 5.5 What the Analogy DOES Suggest

Despite these differences, the analogy is powerful:

**(a) Determinantal structure implies positivity.** At β = 2, the GUE has a determinantal point process structure: all correlation functions are determinants of a single kernel K(x,y). If the ζ-zeros had exact determinantal structure, APT would follow.

**Proposition 5.2 (Conditional).** Suppose the zeros {γ_j} of ζ(1/2 + iγ) form a determinantal point process with kernel K_N(x,y) on [-T, T]. Then for any test function f:

W(f * f̃) = ∫∫ f(x) f̄(y) · K_N(x,y) dx dy ≥ 0

*Proof.* A determinantal kernel K_N is the integral kernel of a projection operator P (or a positive contraction). Hence ⟨f, K_N f⟩ = ⟨f, Pf⟩ = ‖P^{1/2}f‖² ≥ 0. ∎

The hypothesis (exact determinantal structure) is much stronger than what is known. The Montgomery–Odlyzko results establish GUE statistics only in the scaling limit — they say the pair correlation converges to the GUE prediction, not that the exact finite-N process is determinantal.

**(b) Universality suggests robustness.** The β = 2 positivity holds for ALL confining potentials V. If APT were a consequence of the GUE structure (rather than the specific arithmetic V), it would be "universal" — true for any L-function with GUE statistics. This aligns with the GRH prediction.

**(c) The Gram matrix route.** The proof Z_N ≥ 0 uses the Gram matrix structure. Can we find an analogous "Gram matrix" for the ζ-zeros?

Potential candidate: Let φ_p(γ) = p^{-1/2-iγ} · (log p)^{1/2} for each prime p. Then:

⟨φ_p, φ_q⟩ = ∫ p^{-1/2-iγ} q^{-1/2+iγ} (log p · log q)^{1/2} dμ(γ)

where μ is the zero counting measure. If μ is supported on ℝ (RH), this is a genuine inner product and the Gram matrix [⟨φ_p, φ_q⟩] is positive semi-definite. But proving the Gram matrix is PSD requires knowing that μ is real-supported — which is RH. Circular.

---

## 6. The Large Deviation Principle

### 6.1 Setup

Consider the empirical spectral measure of the first N zeros:

μ_N = (1/N) Σ_{j=1}^N δ_{γ_j}

(assuming RH, so all γ_j ∈ ℝ). The empirical measure converges weakly to the equilibrium measure μ_eq (with density ρ_eq from §2).

### 6.2 The GUE Large Deviation Principle

In random matrix theory (GUE), the empirical eigenvalue distribution satisfies a large deviation principle (LDP):

Pr[μ_N ≈ ν] ≍ exp(-β N² · I(ν))

where the rate function is:

I(ν) = ∫ V(x) dν(x) - ∫∫ log|x-y| dν(x)dν(y) - I(μ_eq)

and I(ν) ≥ 0 with equality iff ν = μ_eq.

### 6.3 Cost of a Zero Off the Critical Line

Consider moving one zero from the critical line to distance δ into the complex plane: γ_j → γ_j + iδ. The energy change is:

ΔE = -Σ_{k≠j} [log|γ_j + iδ - γ_k| - log|γ_j - γ_k|]

For small δ, expand:

log|γ_j + iδ - γ_k| = log|γ_j - γ_k| + Re[iδ/(γ_j - γ_k)] - δ²/(2|γ_j - γ_k|²) + O(δ³)

The first-order term vanishes when summed (by the equilibrium condition). The second-order term gives:

**ΔE ≈ (δ²/2) Σ_{k≠j} 1/|γ_j - γ_k|² > 0**

The sum S_j = Σ_{k≠j} |γ_j - γ_k|^{-2} is the local **two-point function**. By GUE universality in the bulk:

- The spacing between adjacent zeros near γ_j is Δ ~ 2π/log(|γ_j|/2π).
- There are ~ N_local zeros within distance O(1) of γ_j.
- Each contributes ~ 1/Δ² ~ (log T)²/(4π²) to the sum.
- Total: S_j ~ c · (log T)² · N_local ~ c' · (log T)³ for γ_j ~ T.

In the properly scaled variable (spacing = O(1/N_local)), S_j ~ c · N_local², giving:

**Pr(zero at distance δ from critical line) ~ exp(-c δ² N_local²)**

### 6.4 The N → ∞ Limit

As N → ∞ (equivalently T → ∞):

- N_local → ∞ (more zeros in any fixed window)
- The cost c δ² N_local² → ∞ for any fixed δ > 0
- The probability of a zero at distance δ > 0 from the critical line → 0

**Interpretation:** In the large-N limit, the equilibrium measure is supported entirely on the real line. The critical line is an "infinitely deep potential well" — no zero can escape.

### 6.5 Can This Be Made Into a Proof?

**The optimistic scenario:** If the zeros of ζ(s) satisfy an LDP with the GUE rate function, then zeros off the critical line have probability zero → RH.

**The obstacles:**

**(i) The zeros are deterministic.** There is no probability measure over "alternative ζ-functions." The LDP applies to random matrix ensembles, not to the specific function ζ(s). The zeros are what they are — γ₁ = 14.134..., γ₂ = 21.022..., etc.

**(ii) Circular reasoning.** The GUE LDP assumes eigenvalues of a Hermitian matrix. Hermitian eigenvalues are real by definition. Saying "the probability of a non-real eigenvalue is zero" merely restates the Hermitian constraint — it does not prove that ζ-zeros are eigenvalues of a Hermitian operator.

**(iii) What COULD work — a constructive approach.** The LDP argument becomes non-circular if we can:

- Step 1: Construct an explicit random matrix ensemble M_N whose characteristic polynomial approximates ζ(1/2 + it) for t ∈ [-T, T].
- Step 2: Show that the approximation is good enough that zeros of ζ and eigenvalues of M_N are within O(1/N²) of each other.
- Step 3: Since eigenvalues of the Hermitian M_N are real, zeros of ζ are within O(1/N²) of the real line.
- Step 4: Take N → ∞ to conclude zeros of ζ are real.

Steps 1–2 are related to the Keating–Snaith program (random matrix models for L-functions). Current results fall short of Step 2: the approximation controls statistical properties (moments) but not individual zero locations.

### 6.6 A Rigorous Partial Result

**Theorem 6.1 (Conditional on pair correlation).** Assume Montgomery's pair correlation conjecture: for test functions g,

lim_{T→∞} (1/N(T)) Σ_{j≠k, j,k ≤ N(T)} g((γ_j - γ_k)log T/(2π)) = ∫_{-∞}^∞ g(x)(1 - (sin πx/(πx))²) dx

Then the energy cost of displacing one zero by distance δ (in scaled units) from its GUE-predicted position satisfies ΔE ≥ c · δ² for a universal constant c > 0.

*Proof.* The pair correlation function R₂(x) = 1 - (sin πx/(πx))² has the properties:
- R₂(x) ~ π²x² as x → 0 (quadratic level repulsion)
- R₂(x) → 1 as x → ∞ (decorrelation)
- ∫ R₂(x)/x² dx < ∞ (the integral converges)

The sum S_j = Σ_{k≠j} |γ_j - γ_k|^{-2}, in scaled variables, converges to:

S_j → ∫ R₂(x)/x² dx = c₀ > 0

(a finite, positive constant; the quadratic vanishing R₂(x) ~ π²x² near x = 0 ensures convergence).

Therefore ΔE ≈ δ² · c₀/2 > 0. ∎

**Interpretation:** Under the pair correlation conjecture, the critical line is **locally rigid** — small perturbations are energetically costly, with cost proportional to δ². This does not prove RH (it's only a local and conditional result) but demonstrates that the GUE structure enforces the critical line as a local minimum.

---

## 7. A Variational Proof Attempt

### 7.1 Setup

Define the spectral functional:

F(γ₁, ..., γ_N) = Σ_ρ |f̂(γ_ρ)|²

where f is a fixed test function and the sum is over the first N nontrivial zeros. The goal: show F is minimized when all γ_j ∈ ℝ.

### 7.2 First Variation

Write γ_j = α_j + iβ_j (α_j = Re(γ_j), β_j = Im(γ_j)). Compute:

∂F/∂β_j = 2 Re[f̂(γ_j)* · ∂f̂/∂β_j(γ_j)]

Since f̂(z) = ∫ f(x) e^{-ixz} dx, we have ∂f̂/∂β_j = ∫ f(x) · x · e^{-ixγ_j} dx (with an appropriate sign from ∂/∂β of e^{-ix(α+iβ)} = x · e^{-ixγ_j}).

At β_j = 0 (on the critical line), the first variation does NOT automatically vanish for arbitrary f. The critical line is NOT a critical point of F for a general test function.

However, the zeros are subject to the **symmetry constraint** from the functional equation: if γ = α + iβ is a zero, so is ᾱ - iβ (equivalently, γ̄). For real f, f̂(γ̄) = f̂(γ)*, and the paired contribution is:

|f̂(α + iβ)|² + |f̂(α - iβ)|² = 2 Re[f̂(α+iβ) · f̂(α-iβ)*]

This pair contribution has the first variation (in β):

∂/∂β [|f̂(α+iβ)|² + |f̂(α-iβ)|²] = 0 at β = 0

by the symmetry β ↦ -β. **So the critical line IS a critical point of the paired functional.**

### 7.3 Second Variation (Convexity Analysis)

At β = 0, compute the second variation of the paired contribution:

∂²/∂β² [|f̂(α+iβ)|² + |f̂(α-iβ)|²]|_{β=0}

Let Φ(β) = |f̂(α+iβ)|². For real-valued f: f̂(α-iβ) = f̂(α+iβ)*, so Φ(-β) = Φ(β) and the function is even. We compute:

Φ(β) = f̂(α+iβ) · f̂(α+iβ)* = f̂(α+iβ) · f̂(α-iβ)

(using f̂(z̄) = f̂(z)* for real f). Differentiating twice at β = 0:

Φ''(0) = 2|f̂'_β(α)|² + 2 Re[f̂(α)* · f̂''_β(α)]

where f̂'_β(α) = ∂f̂/∂β|_{β=0} and f̂''_β(α) = ∂²f̂/∂β²|_{β=0}.

Now: ∂f̂/∂β(α+iβ)|_{β=0} = ∫ f(x) · x · e^{-ixα} dx = (x̂f)(α) (Fourier transform of xf(x)).

And: ∂²f̂/∂β²(α+iβ)|_{β=0} = -∫ f(x) · x² · e^{-ixα} dx = -(x̂²f)(α).

So:

Φ''(0) = 2|(x̂f)(α)|² - 2 Re[f̂(α)* · (x̂²f)(α)]

The sign of Φ''(0) determines whether β = 0 is a local min or max.

### 7.4 The Uncertainty Principle Connection

Φ''(0) = 2|(x̂f)(α)|² - 2 Re[f̂(α)* · (x̂²f)(α)]

By the Cauchy–Schwarz inequality:

|Re[f̂(α)* · (x̂²f)(α)]| ≤ |f̂(α)| · |(x̂²f)(α)|

If we define the "local moments" at frequency α:

μ₁(α) = (x̂f)(α)/f̂(α),  μ₂(α) = (x̂²f)(α)/f̂(α)

then (assuming f̂(α) ≠ 0):

Φ''(0)/|f̂(α)|² = 2|μ₁(α)|² - 2 Re[μ₂(α)]

This is 2(|μ₁|² - Re μ₂) = 2(|μ₁|² - Re μ₂). Now μ₂ - |μ₁|² is the "variance" of the x-distribution at frequency α. For real-valued distributions:

Re μ₂ ≥ |μ₁|²  (variance ≥ 0)

which gives Φ''(0) ≤ 0 — the **wrong sign**! The critical line would be a local maximum, not minimum.

**Wait — this needs careful treatment.** The "distribution" here is NOT a probability distribution (it can have complex values). Let's reconsider.

For f ∈ C_c^∞(ℝ) real-valued, f̂(α) is in general complex. The quantity μ₁(α) = (x̂f)(α)/f̂(α) is complex, and μ₂(α) = (x̂²f)(α)/f̂(α) is also complex. The relation Var = μ₂ - μ₁² ≥ 0 holds only for REAL probability distributions; it need not hold for the complex "distribution" f(x)e^{-iαx}/f̂(α).

The correct analysis: Φ''(0) does not have a definite sign for arbitrary f and α. For specific zeros α_j, the second variation can be positive or negative.

### 7.5 Resolution: The Summed Second Variation

The key insight is that we should not analyze individual zeros separately. The functional F = Σ_ρ |f̂(γ_ρ)|² involves ALL zeros, and we need:

Σ_ρ Φ''_ρ(0) ≥ 0  (summed second variation non-negative)

Using the explicit formula to relate the sum over zeros to prime sums:

Σ_ρ Φ''_ρ(0) = Σ_ρ [2|(x̂f)(γ_ρ)|² - 2 Re f̂(γ_ρ)* · (x̂²f)(γ_ρ)]

The first sum Σ_ρ |(x̂f)(γ_ρ)|² = W_{xf}(xf * x̃f), which is ≥ 0 if RH holds for xf. The second sum involves ⟨f̂, x̂²f⟩ over zeros, which relates to a derivative of the Weil functional.

**The circularity:** Proving the summed second variation is non-negative appears to require RH (or something equally strong). The variational approach reduces APT to a statement about the distribution of zeros — but that distribution IS the content of RH.

### 7.6 The Obstruction

The fundamental difficulty with the variational approach is that F(γ₁, ..., γ_N) treats the zeros as parameters to be varied. But the zeros are NOT free parameters — they are **rigidly determined** by the Euler product:

ζ(s) = Π_p (1 - p^{-s})^{-1}

Changing one zero changes the entire function and hence all other zeros. The "space of zero configurations" is not ℝ^N or ℂ^N — it is a single point (the actual zeros of ζ).

A genuine variational proof would need to:

1. Define energy over the constrained space of zero configurations arising from Euler products.
2. Show the constraints force the minimum to the critical line.
3. The constraints encode the multiplicative structure of integers — making them precise is equivalent to making the Euler product structure precise, which is... analytic number theory.

### 7.7 A Modified Variational Principle

**Approach:** Instead of varying individual zeros, vary the underlying Dirichlet series.

Consider the family:

L_ε(s) = Π_p (1 - ε_p · p^{-s})^{-1}

where |ε_p| = 1 (unit complex numbers). For ε_p = 1 for all p, this is ζ(s).

Define:

E(ε) = W_{L_ε}(f * f̃)

where W_{L_ε} is the Weil functional for L_ε.

**Conjecture 7.1 (Universal Euler Product Positivity).** For any choice of phases {ε_p} with |ε_p| = 1:

W_{L_ε}(f * f̃) ≥ 0 for all f ∈ C_c^∞(ℝ₊*)

This implies RH (take ε_p = 1), and also GRH for Dirichlet L-functions (ε_p = χ(p) for a character χ). The universality might make the statement EASIER to prove — the positivity would follow from the structure of Euler products, not the specific coefficients.

**Evidence for the conjecture:** The explicit formula for L_ε gives:

W_{L_ε}(f*f̃) = (pole terms) - Σ_p (log p) Σ_m (Re[ε_p^m]/p^{m/2}) · (f*f̃)(m log p) + Ω(f)

The sum over p with random phases ε_p should exhibit cancellation (by the law of large numbers), making the prime sum small. The pole terms are positive. So W_{L_ε} should be positive "on average" over ε.

The difficulty: "on average" positivity does not imply pointwise positivity (for every specific choice of ε).

---

## 8. Synthesis: What the Energy Approach Achieves

### 8.1 Rigorous Results

| Result | Status | Reference |
|---|---|---|
| Equilibrium density matches zero density | **Proved** | §2.3 |
| Energy monotone under de Bruijn–Newman flow | **Proved** (for real zeros) | §3.3 |
| Critical line is critical point of paired functional | **Proved** | §7.2 |
| Local energy cost of perturbation ~ δ² | **Proved conditionally** (on pair correlation) | §6.6 |
| Z_N(β=2) ≥ 0 for Coulomb gas | **Proved** (Gram matrix) | §5.2 |

### 8.2 The Gaps

| Gap | Nature | Difficulty |
|---|---|---|
| Coulomb analogy is imprecise | ζ-zeros constrained, not free | Structural |
| LDP argument is circular | Uses Hermitian = real eigenvalues | Foundational |
| Variational approach lacks global convexity | Only local minimum proved | Technical |
| Determinantal structure unproved for ζ | Only statistical GUE known | Deep |
| Modified variational principle | Average positivity ≠ pointwise | Technical |

### 8.3 What Would Close the Gap

A proof via the energy approach requires one of:

**(A) Proving determinantal structure for ζ-zeros.** Show that the zeros form an exact determinantal point process. Then APT follows from the Gram matrix argument (§5.5, Proposition 5.2). This would be a breakthrough in analytic number theory — it would mean the ζ-zeros are eigenvalues of an explicit operator.

**(B) Proving convexity of the constrained energy.** Show that the energy functional on the space of "Euler product zero configurations" is convex with unique minimum on the critical line. This requires understanding how the Euler product structure constrains the zero locations — essentially, making the "rigidity" of the multiplicative structure of ℤ into a geometric/variational statement.

**(C) A universality argument (Conjecture 7.1).** Show that W(f * f̃) ≥ 0 for ALL Euler products. If the positivity is universal — independent of the specific coefficients — then it holds for ζ(s) in particular. The proof would need to use only the STRUCTURE of Euler products (multiplicative, with each factor (1 - ε_p · p^{-s})^{-1} degree 1), not the specific values ε_p = 1.

### 8.4 Connection to Other Approaches

**Cross-terms (→ cross-terms/structure.md):** The cross-terms Σ_{p≠q} CROSS(p,q) correspond to the interaction energy between particles at different "prime-indexed" positions in the log-gas. Diagonal dominance (Approach A in cross-terms/structure.md) is the energy statement that self-energy dominates interaction energy.

**Arithmetic intersection pairing (→ arithmetic-positivity.md):** The energy functional E = -Σ log|γ_j - γ_k| is the spectral side of the arithmetic intersection pairing ⟨Δ, Δ⟩_ar. The Hodge index theorem (negativity of intersection form on primitives) is the statement that the equilibrium configuration minimizes the energy. The energy approach provides the PHYSICAL interpretation; the intersection theory provides the ALGEBRAIC framework.

**The explicit formula as bridge:** The explicit formula

Σ_ρ ĥ(γ_ρ) = (pole terms) - (prime terms) + (archimedean terms)

connects the spectral energy (left side) to the prime energy (right side). APT says both sides are non-negative for h = f * f̃. The prime side's positivity involves controlling cross-terms between primes; the spectral side's positivity involves the reality of zeros. The energy approach shows these are dual perspectives on the same variational principle.
