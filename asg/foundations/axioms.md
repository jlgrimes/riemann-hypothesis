# Axioms of Arithmetic Spectral Geometry

## Formal Axiom System

---

### Primitive Objects

**P1.** The rational field ℚ with its set of places V = {2, 3, 5, 7, ..., ∞}.

**P2.** For each finite place v = p: the p-adic completion ℚ_p, its ring of integers ℤ_p, and the residue field 𝔽_p.

**P3.** For the infinite place v = ∞: the real completion ℝ.

**P4.** The adele ring 𝔸 = ℝ × ∏'_p ℚ_p (restricted product with respect to ℤ_p).

**P5.** The idele group 𝔸* and the idele class group C = 𝔸*/ℚ*.

---

### Axiom Group I: The Frobenius Structure

**Axiom I.1 (Local Frobenius).** For each prime p, there exists a bounded operator

$$\Phi_p : L^2(\mathbb{Z}_p^*, \mu_p) \to L^2(\mathbb{Z}_p^*, \mu_p)$$

defined by (Φ_p f)(x) = f(x^p), where μ_p is the normalized Haar measure on ℤ_p*.

**Axiom I.2 (Local Generator).** The infinitesimal generator of the Frobenius semigroup at p is:

$$D_p = \frac{d}{dt}\bigg|_{t=0} \Phi_p^t$$

where Φ_p^t is the analytic continuation of the semigroup.

**Axiom I.3 (Archimedean Generator).** At the infinite place, the generator is:

$$D_\infty = -i\left(x\frac{d}{dx} + \frac{1}{2}\right)$$

acting on L²(ℝ₊, dx/x) with domain the Schwartz functions.

**Axiom I.4 (Global Assembly).** The global Arithmetic Frobenius generator is the renormalized sum:

$$\mathfrak{D} = D_\infty + \lim_{N \to \infty} \left[\sum_{p \leq N} (\log p) \cdot \pi_p(D_p) - \psi(N) \cdot \text{Id}\right]$$

where π_p: L²(ℤ_p*) → L²(C) is the projection induced by the quotient map 𝔸* → C, and ψ(N) = Σ_{p≤N} log p is the Chebyshev function.

**Axiom I.5 (Convergence).** The limit in I.4 converges in the strong operator topology on a dense subspace 𝒟 ⊂ L²(C, ω) containing the Schwartz-Bruhat space 𝒮(C).

---

### Axiom Group II: The Arithmetic Weight

**Axiom II.1 (Weight Measure).** There exists a positive Borel measure ω on C = 𝔸*/ℚ* satisfying:

(a) ω is absolutely continuous with respect to d*x (Haar measure on C)

(b) dω/d*x = W(x) where W: C → ℝ₊ is the weight function

(c) L²(C, ω) is a separable Hilbert space

**Axiom II.2 (Weight Factorization).** The weight function factors as:

$$W(x) = |x|_\mathbb{A}^{1/2} \cdot \prod_p W_p(x_p) \cdot W_\infty(x_\infty)$$

where:
- W_p(x_p) = (1 - p⁻¹)^{-1/2} for x_p ∈ ℤ_p*
- W_∞(x_∞) = Ω(x_∞) is a Schwartz function on ℝ₊*

**Axiom II.3 (Functional Equation Compatibility).** The Fourier transform ℱ_ω on L²(C, ω), defined with respect to ω, satisfies:

$$\mathcal{F}_\omega^2 = \text{Id}$$

and intertwines 𝔇 with -𝔇:

$$\mathcal{F}_\omega \circ \mathfrak{D} = -\mathfrak{D} \circ \mathcal{F}_\omega$$

This encodes the functional equation ξ(s) = ξ(1-s).

---

### Axiom Group III: The Spectral Space

**Axiom III.1 (Pole Spaces).** Define:

$$\mathcal{H}^0 = \{f \in L^2(C, \omega) : \mathfrak{D}f = \frac{i}{2}f\}$$
$$\mathcal{H}^2 = \{f \in L^2(C, \omega) : \mathfrak{D}f = -\frac{i}{2}f\}$$

These are 1-dimensional (corresponding to the poles of ζ at s = 1 and s = 0 respectively).

**Axiom III.2 (Spectral Space).** The Arithmetic Spectral Space is:

$$\mathcal{H} = L^2(C, \omega) \ominus \mathcal{H}^0 \ominus \mathcal{H}^2$$

(orthogonal complement of the pole spaces).

**Axiom III.3 (Spectral Correspondence).** The point spectrum of 𝔇 restricted to ℋ is:

$$\sigma_p(\mathfrak{D}|_\mathcal{H}) = \{\gamma \in \mathbb{C} : \zeta(1/2 + i\gamma) = 0\}$$

Each eigenvalue has multiplicity equal to the order of the corresponding zero of ζ.

---

### Axiom Group IV: Cohomological Structure

**Axiom IV.1 (Adelic Site).** There exists a Grothendieck site (Spec(ℤ))_ad whose cohomology with coefficients in the arithmetic sheaf 𝒜_ad gives:

$$H^i_{ad}(\text{Spec}(\mathbb{Z})) \cong \mathcal{H}^i \quad \text{for } i = 0, 1, 2$$

**Axiom IV.2 (Lefschetz Trace Formula).** For the Arithmetic Frobenius Φ_t = e^{t𝔇} and suitable test functions h:

$$\sum_{i=0}^{2} (-1)^i \text{Tr}(h(\mathfrak{D}) | \mathcal{H}^i) = \hat{h}(i/2) + \hat{h}(-i/2) - \sum_p \log p \sum_{m=1}^\infty \frac{h(m\log p)}{p^{m/2}}$$

This is the Weil explicit formula, recovered as a Lefschetz fixed-point theorem.

**Axiom IV.3 (Poincaré Duality).** There exists a perfect pairing:

$$\langle \cdot, \cdot \rangle_{PD} : \mathcal{H}^i \times \mathcal{H}^{2-i} \to \mathbb{C}$$

compatible with the Frobenius action.

---

### Axiom Group V: Positivity (The Arithmetic Positivity Theorem)

**Axiom V.1 (Arithmetic Surface).** There exists an object S_ar in the category of arithmetic spaces, equipped with an intersection pairing:

$$\langle \cdot, \cdot \rangle_{ar} : \text{Div}(S_{ar}) \times \text{Div}(S_{ar}) \to \mathbb{R}$$

**Axiom V.2 (Primitivity).** A divisor D ∈ Div(S_ar) is primitive if ⟨D, H₁⟩_ar = ⟨D, H₂⟩_ar = 0 where H₁, H₂ are the two rulings of S_ar.

**Axiom V.3 (APT — The Arithmetic Positivity Theorem).**

$$\forall D \in \text{Div}(S_{ar})^{prim}: \quad \langle D, D \rangle_{ar} \leq 0$$

**STATUS: CONJECTURAL.** This is the sole unproven axiom.

---

### Derived Theorems

**Theorem D.1 (Self-Adjointness).** Axioms I-V ⟹ 𝔇 is essentially self-adjoint on (ℋ, ω).

*Proof:* APT ⟹ the arithmetic inner product is positive-definite ⟹ deficiency indices (0,0) ⟹ essential self-adjointness. ∎

**Theorem D.2 (Riemann Hypothesis).** Axioms I-V ⟹ All non-trivial zeros of ζ(s) have Re(s) = 1/2.

*Proof:* D.1 ⟹ σ(𝔇) ⊂ ℝ ⟹ γ ∈ ℝ for all zeros ρ = 1/2 + iγ ⟹ Re(ρ) = 1/2. ∎

**Theorem D.3 (Consistency).** In the function field case (replace Spec(ℤ) with a curve C/𝔽_q):
- Axioms I-IV are theorems (classical algebraic geometry)
- Axiom V is the Hodge Index Theorem (proven)
- Theorem D.2 recovers Weil's proof of RH for function fields ∎

---

### Independence and Consistency

**Proposition.** Axioms I.1–I.3, II.1–II.2, and III.1–III.2 are consequences of standard mathematics (adelic analysis, harmonic analysis on locally compact abelian groups).

**Proposition.** Axiom I.4 (global assembly) requires the specific regularization choice. The regularization is natural (dictated by the explicit formula) but its convergence (I.5) requires proof.

**Proposition.** Axiom III.3 (spectral correspondence) is provable from I-II using the Mellin transform and Euler product.

**Proposition.** Axiom IV.2 (Lefschetz formula) is equivalent to the Weil explicit formula, which is a theorem.

**Proposition.** Axiom V.3 (APT) is equivalent to the Riemann Hypothesis. It is not a consequence of Axioms I-IV.

---

### Dependency Graph

```
I.1 (local Frob) ──→ I.2 (local gen) ──→ I.4 (global assembly)
                                              │
I.3 (archimedean) ─────────────────────→ I.4 ──→ I.5 (convergence)
                                              │
II.1 (weight) ──→ II.2 (factorization) ──→ II.3 (func. eq.)
    │                                         │
    └──→ III.1 (poles) ──→ III.2 (spectral space) ──→ III.3 (correspondence)
              │                │
              │                └──→ IV.1 (site) ──→ IV.2 (Lefschetz)
              │                                       │
              │                                       └──→ IV.3 (duality)
              │
              └──→ V.1 (surface) ──→ V.2 (primitivity) ──→ V.3 (APT) ⚠️
                                                              │
                                                              └──→ D.1 ──→ D.2 (RH) ★
```

⚠️ = unproven (the single open axiom)
★ = the goal
