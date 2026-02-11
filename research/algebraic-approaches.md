# Arithmetic Geometry and Algebraic Approaches to the Riemann Hypothesis

## A Comprehensive Research Survey

---

## Table of Contents

1. [The Function Field Analogue](#1-the-function-field-analogue)
2. [Weil's Proof (1948)](#2-weils-proof-1948)
3. [Deligne's Proof (1974)](#3-delignes-proof-1974)
4. [Why the Proof Doesn't Transfer](#4-why-the-proof-doesnt-transfer)
5. [Deninger's Program](#5-deningers-program)
6. [Connes' Approach via the Adele Class Space](#6-connes-approach-via-the-adele-class-space)
7. [The Field with One Element F₁](#7-the-field-with-one-element-f₁)
8. [Absolute Geometry](#8-absolute-geometry)
9. [Mochizuki's Inter-Universal Teichmüller Theory](#9-mochizukis-inter-universal-teichmüller-theory)
10. [The Arithmetic Riemann-Roch Theorem](#10-the-arithmetic-riemann-roch-theorem)
11. [What's Missing: Towards a Proof](#11-whats-missing-towards-a-proof)

---

## 1. The Function Field Analogue

### 1.1 Setup and Motivation

The single most important paradigm for approaching the Riemann Hypothesis algebraically is the **analogy between number fields and function fields**. This analogy, elevated to a guiding principle by Weil, proceeds as follows:

| Number Field Side | Function Field Side |
|---|---|
| ℤ | 𝔽_q[t] |
| ℚ | 𝔽_q(t) |
| Spec(ℤ) | Affine line 𝔸¹ over 𝔽_q |
| primes p | monic irreducible polynomials |
| ζ(s) (Riemann) | Z(C/𝔽_q, s) (Hasse-Weil) |
| |p| = p | |f| = q^(deg f) |

The function field analogue of RH has been **proved** — first by Hasse for elliptic curves (genus 1, 1933), then by Weil for all curves (1948), and finally by Deligne for arbitrary varieties (1974). Understanding *why* these proofs work, and *why* they don't generalize, is the central question of algebraic approaches to RH.

### 1.2 Zeta Functions of Curves over Finite Fields

Let C be a smooth, projective, geometrically irreducible curve of genus g over 𝔽_q. The **zeta function** of C is defined as:

$$Z(C/\mathbb{F}_q, T) = \exp\left(\sum_{n=1}^{\infty} \frac{|C(\mathbb{F}_{q^n})|}{n} T^n\right)$$

where T = q^(-s) and |C(𝔽_{q^n})| denotes the number of 𝔽_{q^n}-rational points.

Equivalently, there is an Euler product:

$$Z(C/\mathbb{F}_q, T) = \prod_{x \in |C|} \frac{1}{1 - T^{\deg(x)}}$$

where |C| denotes the set of closed points (i.e., orbits of geometric points under Frobenius), and deg(x) = [κ(x) : 𝔽_q].

### 1.3 The Weil Conjectures for Curves

**Theorem (Weil, 1948).** Let C/𝔽_q be a smooth projective curve of genus g. Then:

**(i) Rationality.** Z(C/𝔽_q, T) is a rational function of T:

$$Z(C/\mathbb{F}_q, T) = \frac{P(T)}{(1-T)(1-qT)}$$

where P(T) ∈ ℤ[T] has degree 2g.

**(ii) Functional Equation.** Z satisfies:

$$Z(C/\mathbb{F}_q, 1/(qT)) = q^{1-g} T^{2-2g} Z(C/\mathbb{F}_q, T)$$

Equivalently, if P(T) = ∏(1 - α_i T), then α_i ↦ q/α_i is a bijection on the roots.

**(iii) Riemann Hypothesis.** The roots α₁, ..., α_{2g} of P(T) satisfy:

$$|α_i| = q^{1/2} \quad \text{for all } i = 1, \ldots, 2g.$$

Equivalently, writing P(T) = ∏(1 - α_i T), all zeros of Z(C, q^{-s}) in the critical strip have Re(s) = 1/2.

### 1.4 Consequences: Point Counting

The Riemann Hypothesis for curves gives the **Hasse-Weil bound**:

$$\left| |C(\mathbb{F}_{q^n})| - (q^n + 1) \right| \leq 2g \cdot q^{n/2}$$

This is the function field analogue of the prime number theorem with optimal error term. For n = 1:

$$|C(\mathbb{F}_q)| = q + 1 - \sum_{i=1}^{2g} \alpha_i, \quad |\alpha_i| = \sqrt{q}$$

### 1.5 The Cohomological Framework

The modern understanding, via Grothendieck, interprets these results through **ℓ-adic étale cohomology**. For a prime ℓ �174 not dividing q:

$$Z(C/\mathbb{F}_q, T) = \prod_{i=0}^{2} \det(1 - T \cdot \mathrm{Frob}_q \mid H^i_{\text{ét}}(C_{\overline{\mathbb{F}_q}}, \mathbb{Q}_\ell))^{(-1)^{i+1}}$$

where:
- H⁰ = ℚ_ℓ contributes (1-T)
- H² = ℚ_ℓ(-1) contributes (1-qT)
- H¹ has dimension 2g and contributes P(T)

The Riemann Hypothesis becomes a statement about the **eigenvalues of Frobenius on H¹**: they all have absolute value q^{1/2}.

---

## 2. Weil's Proof (1948)

### 2.1 Overview of Strategy

Weil's proof is a tour de force of algebraic geometry over finite fields, relying on:
1. The theory of correspondences on C × C
2. The Frobenius endomorphism and its properties
3. Positivity results from intersection theory (the Castelnuovo-Severi inequality)
4. The Hodge index theorem in its algebraic incarnation

The key insight is that the Riemann Hypothesis is equivalent to a **positivity statement** in the intersection theory of the surface C × C.

### 2.2 The Frobenius Endomorphism

Let C/𝔽_q be our smooth projective curve. The **q-Frobenius endomorphism** is the morphism:

$$\mathrm{Fr}_q: C \to C, \quad (x_0 : x_1 : \ldots : x_n) \mapsto (x_0^q : x_1^q : \ldots : x_n^q)$$

This is a purely inseparable morphism of degree q. The fixed points of Fr_q^n are precisely the 𝔽_{q^n}-rational points:

$$C(\mathbb{F}_{q^n}) = \mathrm{Fix}(\mathrm{Fr}_q^n)$$

### 2.3 Correspondences and the Algebra of Divisors on C × C

A **correspondence** on C is a divisor on C × C (up to algebraic equivalence). The key correspondences are:

- **Δ**: the diagonal {(x,x) : x ∈ C}, the identity correspondence
- **Γ_n**: the graph of Fr_q^n, i.e., {(x, Fr_q^n(x)) : x ∈ C}
- **Γ_n^t**: the transpose of Γ_n

The intersection number Γ_n · Δ (on the surface C × C) counts the fixed points of Fr_q^n with multiplicity:

$$\Gamma_n \cdot \Delta = |C(\mathbb{F}_{q^n})|$$

More generally, for correspondences D₁, D₂ on C × C, we have the trace formula:

$$D_1 \cdot D_2 = d_1 d_2' + d_1' d_2 - \sum_{i=1}^{2g} \alpha_i(D_1)\overline{\alpha_i(D_2)}$$

where d, d' are the degrees of the correspondence projected to each factor.

### 2.4 The Key Computation

For the Frobenius correspondence Γ = Γ₁ (graph of Fr_q):

$$\Gamma \cdot \Delta = q + 1 - \sum_{i=1}^{2g} \alpha_i$$

This is the Lefschetz trace formula: the number of fixed points equals the alternating sum of traces on cohomology. Here:
- Trace on H⁰ = 1
- Trace on H² = q
- Trace on H¹ = ∑ α_i

So |C(𝔽_q)| = q + 1 - ∑ α_i.

### 2.5 The Castelnuovo-Severi Inequality

The crucial positivity result is the **Castelnuovo-Severi inequality** (also called the Hodge index theorem for surfaces, or the algebraic analogue of the Cauchy-Schwarz inequality for the intersection pairing).

**Theorem (Castelnuovo-Severi Inequality).** Let C be a curve of genus g over an algebraically closed field. For any correspondence D on C × C of degree d (resp. d') on the first (resp. second) factor:

$$(D \cdot D) \leq 2dd'$$

More precisely, if D is algebraically equivalent to d·(C × pt) + d'·(pt × C) + D₀ where D₀ has degree 0 on both factors, then:

$$(D_0 \cdot D_0) \leq 0$$

This is equivalent to the **negative semi-definiteness** of the intersection pairing on the trace-zero part of the Néron-Severi group of C × C (restricted to correspondences of degree zero on both factors).

### 2.6 Deriving the Riemann Hypothesis

Consider the correspondence D_n = nΔ - Γ for integer n. Then:

- Degree on first factor: n - q
- Degree on second factor: n - 1

The self-intersection:

$$D_n \cdot D_n = (n\Delta - \Gamma) \cdot (n\Delta - \Gamma) = n^2(\Delta \cdot \Delta) - 2n(\Delta \cdot \Gamma) + (\Gamma \cdot \Gamma)$$

Using the known values:
- Δ · Δ = 2 - 2g (by adjunction on C × C)
- Γ · Δ = |C(𝔽_q)| = q + 1 - ∑α_i
- Γ · Γ = q(2 - 2g) (since Γ^t · Γ = degree of Γ^t ∘ Γ)

The Castelnuovo-Severi inequality applied to the degree-zero-on-both-factors part gives:

$$D_n \cdot D_n \leq 2(n-q)(n-1) \quad \text{for all } n$$

Expanding and simplifying, this yields:

$$n^2(2-2g) - 2n(q+1-\sum \alpha_i) + q(2-2g) \leq 2(n-q)(n-1)$$

After simplification:

$$2n\sum \alpha_i \leq 2n(q+1) + 2g(n^2 + q) - 2(n^2 - n(q+1) + q)$$

By applying this for all integers n (or more carefully, for the correspondence nΔ - Γ_m for all m and n), one obtains:

$$|\alpha_i| \leq q^{1/2} \quad \text{for all } i$$

Combined with the functional equation (which gives α_i · α_{2g+1-i} = q), this yields |α_i| = q^{1/2}.

### 2.7 Why Positivity is the Heart

The entire proof reduces to one key fact: **the intersection pairing on the surface C × C has the right signature**. Specifically, the Hodge index theorem says that on a surface, the intersection pairing is negative definite on the orthogonal complement of an ample class. This positivity, when applied to carefully chosen correspondences involving Frobenius, yields the bound on eigenvalues.

This is the template that all algebraic approaches to the classical RH attempt to replicate.

---

## 3. Deligne's Proof (1974)

### 3.1 The Full Weil Conjectures

For a smooth projective variety X/𝔽_q of dimension d, the Weil conjectures (proved by Deligne in 1974, building on Grothendieck's framework) state:

**(i) Rationality** (Dwork, 1960):

$$Z(X/\mathbb{F}_q, T) = \frac{P_1(T) P_3(T) \cdots P_{2d-1}(T)}{P_0(T) P_2(T) \cdots P_{2d}(T)}$$

where P_i(T) = det(1 - T·Frob_q | H^i_ét(X_{ℤar{F}_q}, ℚ_ℓ)).

**(ii) Functional equation:**

$$Z(X, 1/(q^d T)) = \pm q^{d\chi/2} T^{\chi} Z(X, T)$$

where χ = ∑(-1)^i dim H^i is the Euler characteristic.

**(iii) Riemann Hypothesis** (Deligne, 1974): The eigenvalues α of Frob_q on H^i_ét have:

$$|\alpha| = q^{i/2}$$

for every archimedean absolute value. We say α has **weight i**.

**(iv) Betti numbers:** dim H^i_ét(X_{ℤar{F}_q}, ℚ_ℓ) equals the i-th Betti number of the corresponding complex variety (when X lifts to characteristic 0).

### 3.2 Grothendieck's Cohomological Framework

Grothendieck's program reduced the Weil conjectures to properties of **ℓ-adic étale cohomology** H^i_ét(X, ℚ_ℓ). The key tools:

**Étale cohomology** is defined as the inverse limit:

$$H^i_{\text{ét}}(X, \mathbb{Q}_\ell) = \left(\varprojlim_n H^i_{\text{ét}}(X, \mathbb{Z}/\ell^n\mathbb{Z})\right) \otimes_{\mathbb{Z}_\ell} \mathbb{Q}_\ell$$

Key properties that make the proof work:
- **Lefschetz trace formula**: |X(𝔽_{q^n})| = ∑_i (-1)^i Tr(Frob_q^n | H^i)
- **Poincaré duality**: H^i × H^{2d-i} → H^{2d} ≅ ℚ_ℓ(-d)
- **Künneth formula**: H^*(X × Y) ≅ H^*(X) ⊗ H^*(Y)
- **Lefschetz hyperplane theorem**: controls cohomology under hyperplane sections
- **Hard Lefschetz theorem**: the operator L^{d-i}: H^i → H^{2d-i} is an isomorphism

Grothendieck showed that the Weil conjectures would follow from a **standard conjectures** program. However, the standard conjectures remain unproved in general.

### 3.3 Deligne's Strategy

Deligne's proof bypassed the standard conjectures by using a brilliant inductive argument:

**Step 1.** Reduce to showing that eigenvalues of Frobenius on H^i have absolute value ≤ q^{i/2} (the upper bound suffices, since the lower bound follows from Poincaré duality).

**Step 2.** Use a Lefschetz pencil: embed X in projective space, take hyperplane sections X_t = X ∩ H_t for t varying in ℙ¹. This gives a fibration π: X̃ → ℙ¹ (after blowing up the base locus).

**Step 3.** Analyze the higher direct images R^i π_* ℚ_ℓ as ℓ-adic sheaves on ℙ¹. Apply the Leray spectral sequence:

$$E_2^{p,q} = H^p(\mathbb{P}^1, R^q \pi_* \mathbb{Q}_\ell) \Rightarrow H^{p+q}(X, \mathbb{Q}_\ell)$$

**Step 4.** The key innovation: use the **Rankin-Selberg method** (a convolution trick). For a smooth sheaf ℱ on an open U ⊆ ℙ¹, consider the tensor product ℱ ⊗ ℱ^∨. The L-function L(ℱ ⊗ ℱ^∨, s) has a pole at s = 1 (from the trivial component), and this pole controls the absolute values of eigenvalues.

**Step 5.** Use Grothendieck's theory of L-functions for ℓ-adic sheaves and the Euler-Poincaré formula to relate the order of the pole to the sum of squares of absolute values of eigenvalues. The rationality and integrality constraints then force the eigenvalues to have the correct absolute value.

### 3.4 The Concept of Weights

Deligne introduced the notion of **weight** for ℓ-adic sheaves:

**Definition.** A smooth ℚ_ℓ-sheaf ℱ on a scheme X₀ of finite type over 𝔽_q is **pure of weight w** if for every closed point x ∈ X₀ with residue field 𝔽_{q^{deg(x)}}, all eigenvalues of Frob_x on the stalk ℱ_x̄ have absolute value q^{w·deg(x)/2}.

**Definition.** ℱ is **mixed of weight ≤ w** if it has a filtration whose graded pieces are pure of weight ≤ w.

**Theorem (Deligne, Weil II, 1980).** Let f: X → Y be a morphism of schemes of finite type over 𝔽_q. If ℱ on X is mixed of weight ≤ w, then R^i f_* ℱ is mixed of weight ≤ w + i.

This vastly generalizes the original Weil conjectures and provides a powerful machine for computing weights in geometric situations.

### 3.5 Significance for RH

Deligne's proof shows that the Riemann Hypothesis over finite fields is fundamentally a consequence of:
1. The existence of a good cohomology theory (étale cohomology)
2. The Lefschetz trace formula
3. Hard Lefschetz / Poincaré duality
4. A clever analytic argument (the Rankin-Selberg trick) that converts positivity/integrality into the weight bound

The question for classical RH is: **can any of this be transported to Spec(ℤ)?**

---

## 4. Why the Proof Doesn't Transfer

### 4.1 The Fundamental Obstructions

The analogy between function fields and number fields, while extremely powerful as a source of conjectures and intuition, breaks down at precisely the points needed to prove RH. Here we catalogue the obstructions.

### 4.2 Obstruction 1: The Absence of Frobenius

Over 𝔽_q, the Frobenius endomorphism Fr_q: x ↦ x^q is:
- An endomorphism of every 𝔽_q-scheme
- Generates Gal(𝔽̄_q/𝔽_q) ≅ ℤ̂
- Has a well-defined action on cohomology whose eigenvalues encode arithmetic data

Over ℤ, there is **no single endomorphism** playing this role. The "Frobenius at p" exists for each prime p (as an element of the decomposition group D_p ⊂ Gal(ℚ̄/ℚ)), but:
- Different primes give different Frobenius elements
- They don't cohere into a single global endomorphism
- Gal(ℚ̄/ℚ) is not procyclic — it's enormously complicated

**The lack of a single Frobenius means there is no single operator whose eigenvalues we need to bound.** Instead, we would need to simultaneously control the eigenvalues at every prime.

### 4.3 Obstruction 2: The Missing Cohomology Theory

For curves over 𝔽_q, ℓ-adic étale cohomology provides:
- Finite-dimensional vector spaces H^i
- A Frobenius action with the right traces (Lefschetz formula)
- Poincaré duality, Künneth, functoriality

For Spec(ℤ), we would need a cohomology theory H^i(Spec(ℤ)) such that:

$$\xi(s) = \pi^{-s/2}\Gamma(s/2)\zeta(s) = \frac{\det(s \cdot \mathrm{Id} - \Theta \mid H^1)}{\det(s \cdot \mathrm{Id} - \Theta \mid H^0) \cdot \det(s \cdot \mathrm{Id} - \Theta \mid H^2)}$$

where Θ is some operator playing the role of "log Frobenius." This would require:
- H⁰ = one-dimensional, eigenvalue of Θ is 0
- H² = one-dimensional, eigenvalue of Θ is 1
- H¹ = **infinite-dimensional**, with eigenvalues at the zeros ρ of ζ(s)

No such cohomology theory is known to exist. This is the central goal of the programs of Deninger and Connes.

### 4.4 Obstruction 3: The Missing Geometric Structure

Spec(ℤ) is a one-dimensional scheme with:
- Generic point η = Spec(ℚ)
- Closed points = primes p (with residue fields 𝔽_p of varying characteristic)
- No uniform characteristic (unlike 𝔽_q[t], which lives in characteristic p)

For a curve C/𝔽_q, the "absolute" surface is C × C, and the diagonal Δ ⊂ C × C is where intersection theory happens. For Spec(ℤ), we would need:
- A "base" Spec(𝔽₁) over which Spec(ℤ) is a curve
- A product Spec(ℤ) ×_{𝔽₁} Spec(ℤ) with good intersection theory
- An analogue of the Hodge index theorem on this "arithmetic surface"

None of these objects have satisfactory definitions, though many candidates have been proposed (see Sections 7 and 8).

### 4.5 Obstruction 4: Positivity

In Weil's proof, the Castelnuovo-Severi inequality provides the crucial positivity. On a surface S, the Hodge index theorem says the intersection pairing on NS(S) ⊗ ℝ has signature (1, ρ-1) where ρ = rank NS(S).

For an "arithmetic surface" (in the sense of Arakelov geometry), the analogous positivity is related to the **Generalized Riemann Hypothesis** (GRH), not merely to RH! In Arakelov geometry, the arithmetic intersection pairing does exist, but proving the analogue of the Hodge index theorem for arithmetic surfaces **is equivalent to GRH** (by work of Szpiro, Ullmo, and Zhang).

This is a stunning circularity: the positivity needed for the proof IS the statement we're trying to prove.

### 4.6 Obstruction 5: The Archimedean Place

Over 𝔽_q, all places are non-archimedean. The completed zeta function Z(C, T) has no gamma factors. Over ℚ, the archimedean place ℝ (or ℂ) contributes:
- The gamma factor Γ(s/2) in the functional equation
- An entire theory of Hodge structures, periods, and regulators
- The difficulty of handling "analysis at infinity"

Incorporating the archimedean place into a purely algebraic framework is one of the deepest challenges. Arakelov geometry (Section 10) is the most developed attempt.

### 4.7 Summary Table of Obstructions

| What's needed | Function field status | Number field status |
|---|---|---|
| Frobenius endomorphism | Fr_q exists globally | Only local Frob_p |
| Finite-dimensional cohomology | ℓ-adic étale cohomology | No known theory |
| Base field 𝔽₁ | 𝔽_q exists | Hypothetical |
| Product C ×_{base} C | C ×_{𝔽_q} C | Spec(ℤ) ×_{𝔽₁} Spec(ℤ)? |
| Positivity (Hodge index) | Castelnuovo-Severi | Equivalent to GRH! |
| Trace formula | Lefschetz trace formula | Weil explicit formula (partial) |
| All places non-archimedean | Yes | No (ℝ is archimedean) |

---

## 5. Deninger's Program

### 5.1 Motivation and Vision

Christopher Deninger, beginning in 1991, proposed a conjectural cohomological framework that would explain the "shape" of the zeta function and potentially lead to a proof of RH. His program seeks to construct a cohomology theory that makes the analogy with function fields precise.

### 5.2 The Expected Properties

Deninger postulated the existence of:

**(A) Infinite-dimensional cohomology groups** H^i(X̄, ℝ) for arithmetic schemes X (where X̄ is a compactification including the archimedean places), equipped with:

**(B) A "flow" Φ_t** (a one-parameter group of automorphisms, the analogue of Frobenius), whose infinitesimal generator Θ plays the role of "log Frobenius."

**(C) A regularized determinant** formula:

$$\xi(s) = \frac{\det_\infty(s \cdot \mathrm{Id} - \Theta \mid H^1(\overline{\mathrm{Spec}(\mathbb{Z})}))}{det_\infty(s \cdot \mathrm{Id} - \Theta \mid H^0(\overline{\mathrm{Spec}(\mathbb{Z})})) \cdot \det_\infty(s \cdot \mathrm{Id} - \Theta \mid H^2(\overline{\mathrm{Spec}(\mathbb{Z})}))}$$

where det_∞ is a regularized (zeta-regularized) determinant in the sense of Ray-Singer.

### 5.3 Regularized Determinants

For an operator Θ on an infinite-dimensional space with discrete spectrum {λ_n}, the **zeta-regularized determinant** is defined via:

$$\zeta_\Theta(w) = \sum_n \lambda_n^{-w}, \quad \det_\infty(s - \Theta) = \exp\left(-\frac{d}{dw}\Big|_{w=0} \sum_n (s - \lambda_n)^{-w}\right)$$

Deninger showed that if Θ has eigenvalues at the nontrivial zeros ρ of ζ(s), then:

$$\det_\infty(s - \Theta \mid H^1) \sim \xi(s)$$

and the simple eigenvalues 0 (on H⁰) and 1 (on H²) produce the trivial factors.

### 5.4 The Explicit Formulas as Trace Formulas

A key insight of Deninger is that the **Weil explicit formula** should be interpreted as a Lefschetz trace formula for the flow Φ_t:

$$\sum_\rho h(\rho) = h(0) + h(1) - \sum_p \sum_{m=1}^{\infty} \frac{\log p}{p^{m/2}} \left(h\left(\frac{1}{2} + \frac{im\log p}{\log q}\right) + \cdots \right)$$

The left side is the spectral side (sum over eigenvalues of Θ = zeros of ζ), and the right side is the geometric side (sum over "closed orbits" = prime powers). This is directly analogous to the Selberg trace formula and to the Lefschetz fixed-point formula for flows.

### 5.5 Deninger's Explicit Constructions

In more recent work (2000s-2020s), Deninger has attempted to construct parts of this framework:

**Foliated dynamical systems.** Deninger showed that laminated spaces (solenoids) with transverse measures can produce zeta functions via dynamical Lefschetz trace formulas. Specifically:

- Consider a foliation (S, ℱ) with a transverse measure μ
- The "foliation zeta function" ζ(S, ℱ, s) can be defined via regularized determinants of the leafwise Laplacian
- For specific foliations, this reproduces arithmetic zeta functions

**The arithmetic site.** In joint work (and influenced by Connes-Consani), Deninger has considered:
- The space of adele classes as a candidate geometric space
- Foliated structures where leaves correspond to primes
- The "arithmetic topos" approach where Spec(ℤ) is viewed through the lens of topos theory

### 5.6 What the Program Has Achieved

- Provides a **conceptual framework** explaining WHY ζ(s) should have the properties it has
- Shows that RH is equivalent to a **positivity statement** for the Θ operator: namely that the eigenvalues of Θ on H¹ lie on Re(s) = 1/2, which would follow from Θ being self-adjoint with appropriate shifts
- Correctly predicts the shape of gamma factors and functional equations
- Illuminates the role of the archimedean place

### 5.7 What Remains

- The cohomology groups H^i have not been constructed as concrete mathematical objects
- The "flow" Φ_t has not been realized in a way that yields arithmetic information
- The positivity/self-adjointness of Θ remains conjectural
- No mechanism to deduce RH from the framework has been found, even conditionally

---

## 6. Connes' Approach via the Adele Class Space

### 6.1 Overview

Alain Connes has developed the most ambitious and detailed program connecting noncommutative geometry to the Riemann Hypothesis. The central idea is to construct a geometric space whose "geometry" encodes the zeros of ζ(s).

### 6.2 The Adele Class Space

**Definition.** The **adele class space** of ℚ is:

$$C_\mathbb{Q} = \mathbb{A}_\mathbb{Q}^* / \mathbb{Q}^*$$

where 𝔸_ℚ = ℝ × ∏'_p ℚ_p is the ring of adeles of ℚ and 𝔸_ℚ^* is the group of ideles.

More precisely, Connes considers the quotient space:

$$X_\mathbb{Q} = \mathbb{A}_\mathbb{Q} / \mathbb{Q}^*$$

where ℚ* acts by multiplication. This is a "bad" quotient in classical topology (not Hausdorff), but is well-suited to noncommutative geometry.

### 6.3 The Spectral Realization

**Theorem (Connes, 1999).** There exists a natural representation of the idele class group C_ℚ = 𝔸_ℚ*/ℚ* on a Hilbert space ℋ, and a self-adjoint operator D (playing the role of "Dirac operator" or "1/2 + iΘ"), such that:

The **absorption spectrum** of D (the values λ where the representation degenerates) is precisely the set:

$$\{\rho : \zeta(\rho) = 0, \; 0 < \mathrm{Re}(\rho) < 1\}$$

More precisely, the spectral realization works as follows:

**Step 1.** Consider the space L²(C_ℚ) of square-integrable functions on the idele class group.

**Step 2.** Define the "Weil operator" W as the action by scaling: (W_λ f)(x) = f(λ^{-1}x) for λ ∈ C_ℚ.

**Step 3.** The scaling action by ℝ₊* ⊂ C_ℚ gives a one-parameter group whose generator H has the property:

$$\mathrm{Tr}(f(H)) = \sum_\rho \hat{f}(\rho) + \text{(terms from poles and archimedean factors)}$$

for suitable test functions f, where the sum is over zeros ρ of ζ(s).

### 6.4 Connection to Weil's Explicit Formula

The Weil explicit formula in its distributional form states:

$$\sum_\rho h(\rho) = h(0) + h(1) - \sum_v \int_{\mathbb{Q}_v^*} \frac{h(|x|_v)}{|1-x|_v} d^*x$$

Connes showed that this formula is exactly the **trace formula** for the action of C_ℚ on a suitable space. The left side is the spectral side, the right side is the geometric side, and the sum over places v (both finite primes p and the archimedean place ∞) gives the orbital integrals.

This is the arithmetic analogue of the Selberg trace formula, and it confirms Deninger's intuition that the explicit formula should be viewed cohomologically.

### 6.5 The Semi-Local Trace Formula and the Weil Positivity Criterion

**Theorem (Weil, 1952; reformulated by Connes).** The Riemann Hypothesis is equivalent to the following positivity condition:

For all f ∈ S(C_ℚ) (Schwartz-Bruhat functions on the idele class group) such that f = f * f^* (i.e., f is a "positive type" function), we have:

$$\sum_\rho \hat{f}(\rho) \geq 0$$

where the sum is over all nontrivial zeros ρ of ζ(s), counted with multiplicity.

Connes reformulated this as: RH ⟺ a certain operator is nonnegative. The challenge is to prove this positivity.

### 6.6 The Connes-Consani Arithmetic Site

Starting around 2014, Connes and Consani introduced the **Arithmetic Site** 𝒜:

**Definition.** The arithmetic site is the topos of sheaves on the small category with:
- Objects: the natural numbers ℕ× (under multiplication)
- A single object ★
- Morphisms: multiplication by n ∈ ℕ×

The structure sheaf is the semiring ℤ_max = (ℤ ∪ {-∞}, max, +), the tropical integers.

**Key properties:**
1. The "points" of 𝒜 over ℝ_max^+ recover the adele class space C_ℚ
2. The scaling action by ℝ₊* corresponds to the "Frobenius flow"
3. The arithmetic site provides a framework where Spec(ℤ) is viewed as a "curve" over 𝔽₁ (in the tropical/characteristic one sense)

### 6.7 The Scaling Site and Weil's Proof Analogy

Connes and Consani further refined this to the **Scaling Site**:

- Objects are given by the action of ℕ× on ℝ₊*
- The topos structure provides sheaf cohomology
- The goal is to implement Weil's proof strategy:
  1. Work on the "surface" 𝒜 ×_{𝔽₁} 𝒜
  2. Identify the diagonal and the "Frobenius correspondence"
  3. Use a form of Riemann-Roch and positivity to bound eigenvalues

### 6.8 Current Status

As of the mid-2020s:
- The spectral realization of zeros is well-established
- The trace formula interpretation is rigorous
- The arithmetic site provides a genuine topos-theoretic framework
- **The crucial positivity step remains unproved**
- The scaling site cohomology is still being developed
- Connes-Consani have published extensive work on the tropical/characteristic one geometry aspects
- A complete proof remains out of reach, but the framework continues to be refined

### 6.9 The Role of Tropical Geometry

A significant insight from Connes-Consani is that the "geometry over 𝔽₁" naturally involves **tropical geometry** and **semirings**:

- The semifield ℝ_max = (ℝ ∪ {-∞}, max, +) replaces ordinary fields
- Tropical curves and their Jacobians provide analogues of algebraic curves
- The Riemann-Roch theorem has tropical analogues (Baker-Norine theorem)
- This tropical/semiring geometry may provide the right categorical framework for 𝔽₁

---

## 7. The Field with One Element 𝔽₁

### 7.1 Motivation

The idea of a "field with one element" 𝔽₁ arose from multiple independent sources:

**Tits' observation (1957).** Many counting formulas for structures over 𝔽_q have meaningful limits as q → 1. For example:
- |GL_n(𝔽_q)| = (q^n - 1)(q^n - q)···(q^n - q^{n-1}) → n! as q → 1
- |Grass(k,n)(𝔽_q)| = [n choose k]_q → (n choose k) as q → 1
- The Weyl group W appears as "GL_n(𝔽₁)"

This suggests that 𝔽₁ is a "field" with one element (just 0), where linear algebra reduces to combinatorics.

**Manin's observation (1995).** If Spec(ℤ) is to be a "curve over 𝔽₁," then:
- Spec(ℤ) ×_{Spec(𝔽₁)} Spec(ℤ) should be a "surface"
- Weil's proof strategy could potentially be applied
- The "genus" of this curve should be related to the zeros of ζ(s)

### 7.2 The Counting Formula

For a smooth projective curve C of genus g over 𝔽_q:

$$|C(\mathbb{F}_q)| = q + 1 - \sum_{i=1}^{2g} \alpha_i$$

If we formally set q = 1, this gives:

$$|C(\mathbb{F}_1)| = 2 - 2g$$

For Spec(ℤ) viewed as a curve over 𝔽₁, comparing with the explicit formula for ζ(s) (via Weil's explicit formula), Manin computed a formal "genus":

$$2 - 2g = 2 + \frac{1}{2\pi} \int_{-\infty}^{\infty} \frac{\Gamma'}{\Gamma}\left(\frac{1}{4} + \frac{it}{2}\right) dt - \sum_p \frac{\log p}{p - 1}$$

This diverges, suggesting that Spec(ℤ) has **infinite genus** as a curve over 𝔽₁ — consistent with the zeta function having infinitely many zeros.

### 7.3 Soulé's Approach (2004)

Soulé proposed a definition of varieties over 𝔽₁ based on the idea that 𝔽₁-structures should be seen in the "underlying combinatorial data" of algebraic varieties:

**Definition (Soulé).** An 𝔽₁-variety is a covariant functor X: FinAb → Sets (from finite abelian groups to sets) together with a complex variety X_ℂ and a natural transformation X(A) → X_ℂ(ℂ) compatible with base change.

The category of 𝔽₁-varieties includes:
- Toric varieties (the prototypical examples)
- Flag varieties, Grassmannians
- Split reductive groups

The base change to ℤ is: X_ℤ = X ⊗_{𝔽₁} ℤ, where the construction glues the combinatorial data with the ℤ-scheme structure.

### 7.4 Connes-Consani's Approach

Connes and Consani proposed several frameworks:

**(a) 𝔽₁-schemes via monoids (2010).** An 𝔽₁-scheme is glued from spectra of monoids (replacing commutative rings with commutative monoids). The base change to ℤ is:
$$\mathrm{Spec}(M) \otimes_{\mathbb{F}_1} \mathbb{Z} = \mathrm{Spec}(\mathbb{Z}[M])$$
where ℤ[M] is the monoid ring.

**(b) Hyperrings and the Krasner hyperfield (2011).** 𝔽₁ is related to the **Krasner hyperfield** 𝕂 = {0, 1} with 1 + 1 = {0, 1}. Hyperring theory provides a framework where:
- Spec(ℤ) has a natural structure over 𝕂
- The projective line ℙ¹(𝕂) classifies valuations

**(c) Λ-rings and Borger's approach (2009).** James Borger proposed that the "descent data from ℤ to 𝔽₁" should be a **Λ-ring structure** (i.e., a collection of commuting Frobenius lifts ψ_p for all primes p):

**Definition (Borger).** An 𝔽₁-algebra is a Λ-ring, i.e., a commutative ring R equipped with commuting ring endomorphisms ψ_p: R → R for all primes p such that ψ_p(x) ≡ x^p (mod p) for all x ∈ R.

This approach has the advantage that:
- It recovers the "Frobenius at every prime" simultaneously
- The category of Λ-rings has excellent formal properties
- Borger's 𝔽₁ base change is: X ↦ X with its canonical Λ-structure
- The Witt vector construction W(𝔽_p) = ℤ_p is naturally a Λ-ring

### 7.5 Borger's Approach in Detail

Borger's framework is perhaps the most algebraically satisfying:

**Key insight.** Over 𝔽_q, the Frobenius endomorphism x ↦ x^q is what makes the whole theory work. Over ℤ, we have "Frobenius at p" for each prime p, but no single Frobenius. Borger's idea: **use ALL the Frobenius lifts simultaneously**.

A Λ-structure on a ring R is:
- Adams operations ψ^n: R → R for all n ≥ 1
- ψ^1 = id, ψ^m ∘ ψ^n = ψ^{mn}
- ψ^p(x) ≡ x^p (mod p) for all primes p
- Equivalently, a coaction of the big Witt vector functor W

**Theorem (Borger).** The category of Λ-rings is equivalent to the category of "𝔽₁-algebras" in a precise sense, and the forgetful functor to commutative rings gives the base change:

$$- \otimes_{\mathbb{F}_1} \mathbb{Z}: \text{Λ-Rings} \to \text{CRing}$$

### 7.6 How Spec(ℤ) Should Be a Curve over 𝔽₁

In the 𝔽₁-picture, the structural analogy is:

| Function field | Number field (𝔽₁ picture) |
|---|---|
| 𝔽_q | 𝔽₁ |
| 𝔽_q[t] | ℤ |
| C/𝔽_q (projective curve) | Spec(ℤ) ∪ {∞} (compactified) |
| 𝔽̄_q = lim 𝔽_{q^n} | "𝔽̄₁" = roots of unity μ_∞ |
| Gal(𝔽̄_q/𝔽_q) ≅ ℤ̂ | Gal(ℚ^{ab}/ℚ) ≅ ℤ̂* (class field theory) |
| Frobenius Fr_q | ??? |

The compactification of Spec(ℤ) by adding the "archimedean place" ∞ (corresponding to the absolute value | · |_∞ on ℚ) gives "Spec(ℤ̄)" — the complete "curve over 𝔽₁."

### 7.7 Current Status of 𝔽₁ Geometry

Despite decades of work:
- Multiple inequivalent definitions of 𝔽₁ exist, each capturing different aspects
- No single framework has proved capable of attacking RH
- The relationship between different approaches is partially understood (through work of López Peña, Lorscheid, and others)
- The "correct" definition likely requires a synthesis of several approaches
- Borger's Λ-ring approach and Connes-Consani's arithmetic site approach appear most promising but remain far from proving RH

---

## 8. Absolute Geometry

### 8.1 The Quest for Foundations

"Absolute geometry" refers to the broader program of developing algebraic geometry in a setting where the base is not a ring or field, but something more primitive. The goal is to find the right generalization of schemes that accommodates 𝔽₁-geometry and ultimately Spec(ℤ) as a curve.

### 8.2 Durov's Generalized Rings

Nikolai Durov (2007) introduced **generalized rings** (also called "generalized commutative algebraic monads"):

**Definition.** A generalized ring is a commutative algebraic monad on the category of sets, i.e., a functor Σ: Sets → Sets equipped with a monoid structure (unit η: Id → Σ, multiplication μ: Σ ∘ Σ → Σ) satisfying commutativity and certain algebraic conditions.

**Key examples:**
- Ordinary commutative rings R give generalized rings via Σ_R(S) = R^{(S)} (free R-module on S)
- 𝔽₁ corresponds to the identity monad Σ(S) = S
- The "field of signs" 𝔽_{±1} corresponds to Σ(S) = {formal sums ±s : s ∈ S}
- ℤ_∞ (integers with archimedean absolute value ≤ 1) gives a generalized ring encoding the archimedean place

**Relevance:** Durov showed that generalized algebraic geometry naturally incorporates both the archimedean and non-archimedean places of ℚ, providing a unified framework for Arakelov geometry.

### 8.3 Haran's F-Rings

Shai Haran proposed **F-rings** (2007, 2009) as foundations for absolute geometry:

**Definition.** An F-ring is a triple (A, A⁺, A×) where:
- A is a set with two operations (addition and multiplication)
- A⁺ and A× encode "non-archimedean" and "archimedean" data
- Axioms unifying the theory of rings with absolute values

The key insight is that the passage from 𝔽_q to 𝔽₁ should be controlled by the interplay between addition and multiplication — and in the 𝔽₁ limit, addition "disappears," leaving only multiplicative structure.

### 8.4 Lorscheid's Blueprints

Oliver Lorscheid (2012) introduced **blueprints** as a compromise between monoids and rings:

**Definition.** A blueprint is a pair B = (A, R) where:
- A is a commutative monoid (with absorbing element 0)
- R ⊆ ℕ[A] × ℕ[A] is a relation (the "blueprint relation") on formal sums of elements of A
- R is an equivalence relation compatible with the monoid structure

**Key property:** Blueprints form a category containing both:
- Commutative monoids (when R is trivial)
- Commutative rings (when R identifies all equal sums)
- Various intermediate structures

**Blue schemes** are obtained by gluing spectra of blueprints, just as ordinary schemes are obtained from rings.

Lorscheid showed that blue schemes provide a unified framework encompassing most other 𝔽₁-geometry proposals, and that many classical constructions (Grassmannians, flag varieties, toric varieties) have natural blue scheme structures.

### 8.5 Toen-Vaquié's Approach

Toen and Vaquié (2009) proposed using **schemes relative to symmetric monoidal categories** (C, ⊗):

**Definition.** A scheme relative to (C, ⊗) is constructed by:
1. Defining commutative monoid objects in C as "generalized rings"
2. Defining their spectra and structure sheaves
3. Gluing to form schemes

Taking C = (Sets, ×) gives "𝔽₁-schemes" (monoid schemes). Taking C = (Ab, ⊗) gives ordinary schemes.

### 8.6 Comparison and Synthesis

The various approaches to absolute geometry are related:

```
                    Durov (generalized rings)
                          ↕
    Soulé (gadgets) ← Lorscheid (blueprints) → Toen-Vaquié (C-schemes)
                          ↕
                    Connes-Consani (Λ-rings, hyperrings, arithmetic site)
                          ↕
                    Borger (Λ-rings)
```

Lorscheid's blueprints have emerged as a partial unifying framework, but significant gaps remain:
- The "right" cohomology theory in absolute geometry is unknown
- Intersection theory for blue schemes is underdeveloped
- The analogue of the Hodge index theorem in this setting has not been formulated

---

## 9. Mochizuki's Inter-Universal Teichmüller Theory

### 9.1 Overview

Shinichi Mochizuki's Inter-Universal Teichmüller Theory (IUT), developed over 2000-2012 and published in a series of four papers (2012, revised through 2021), was initially claimed to prove the ABC conjecture. The theory introduces radically new mathematical structures:

- **Hodge theaters**: elaborate combinatorial/categorical structures encoding the arithmetic of number fields
- **Frobenioids**: categorical frameworks abstracting divisor theory and Frobenius-like operations
- **The Θ-link**: a key construction connecting different "universes" of mathematical structures
- **Mono-anabelian transport**: the process of moving between different copies of arithmetic data

### 9.2 Connection to RH

Mochizuki's IUT theory does **not** directly address the Riemann Hypothesis. The ABC conjecture, if proved, would have some implications for the distribution of primes, but these fall short of RH.

However, IUT is relevant to the broader algebraic approach to RH in the following indirect ways:

**(a) Frobenius-like operations.** The "Frobenioids" in IUT are designed to capture Frobenius-like structures over number fields. This resonates with the need for a "global Frobenius" in algebraic approaches to RH.

**(b) Anabelian geometry.** Mochizuki's work builds on Grothendieck's anabelian program, which seeks to reconstruct arithmetic objects from their fundamental groups. If number fields were fully determined by their étale fundamental groups, this might provide a new handle on their zeta functions.

**(c) Deconstruction and reconstruction.** The IUT philosophy of "dismantling" ring structures (separating addition from multiplication) and reconstructing them through categorical links is philosophically aligned with 𝔽₁-geometry and Haran's F-rings.

### 9.3 Controversy and Status

The mathematical status of IUT is deeply controversial:

- **Scholze-Stix objection (2018):** Peter Scholze and Jakob Stix identified what they consider a fundamental gap in Mochizuki's proof, specifically in Corollary 3.12 of the third IUT paper. They argue that the claimed inequality becomes trivial (giving only 0 ≤ 0) when the indeterminacies are properly tracked.

- **Mochizuki's response:** Mochizuki maintains that Scholze-Stix misunderstand the theory and that their objection does not apply. He has published detailed rebuttals.

- **Community status:** As of the mid-2020s, no mathematical consensus has emerged. The majority of number theorists who have engaged with the material (outside Mochizuki's immediate circle) remain unconvinced. The work was published in PRIMS (a journal where Mochizuki serves as editor-in-chief), which has not resolved the controversy.

- **Impact on RH:** Even if IUT is correct, it does not directly address RH. The techniques, if valid, might eventually contribute to understanding the multiplicative-additive interplay relevant to zeta functions, but this remains speculative.

### 9.4 Relevance to This Survey

IUT is included here primarily for completeness. Its relevance to algebraic approaches to RH is **indirect at best**. The core ideas that might connect to RH — abstracting Frobenius, handling the addition-multiplication interplay, and moving between different arithmetic "universes" — are more directly addressed by the programs of Connes, Deninger, and Borger.

---

## 10. The Arithmetic Riemann-Roch Theorem

### 10.1 Arakelov Geometry: Foundations

Suren Arakelov (1974) introduced a framework for doing intersection theory on arithmetic surfaces by incorporating the archimedean places. The key idea: to compactify Spec(ℤ) (and more generally, arithmetic varieties), we must include data at ∞.

**Definition.** An **arithmetic variety** is a pair (X, {h_σ}) where:
- X is a projective variety over Spec(ℤ) (or more generally, over the ring of integers of a number field)
- {h_σ} is a collection of Hermitian metrics on the complex points X_σ(ℂ) for each archimedean embedding σ

**Definition.** An **Arakelov divisor** on an arithmetic curve C (i.e., a regular model of a curve over ℚ) is a formal sum:

$$\hat{D} = D + \sum_{\sigma} \lambda_\sigma \cdot [\sigma]$$

where D is a Weil divisor on C and λ_σ ∈ ℝ are real numbers associated to the archimedean places σ.

### 10.2 The Arithmetic Intersection Pairing

For an arithmetic surface X → Spec(ℤ_K) (where ℤ_K is the ring of integers of a number field K), Arakelov defined an intersection pairing on arithmetic divisors:

$$\hat{D}_1 \cdot \hat{D}_2 = (D_1 \cdot D_2)_{\text{fin}} + (D_1 \cdot D_2)_{\text{inf}}$$

where:
- The finite part (D₁ · D₂)_fin is the usual intersection number on the generic and special fibers
- The infinite part (D₁ · D₂)_inf involves Green's functions on the Riemann surfaces X_σ(ℂ)

**The Green's function.** For a Riemann surface Σ with a metric μ, the Arakelov-Green's function g_μ(P, Q) satisfies:

$$\partial_P \bar{\partial}_P g_\mu(P, Q) = \pi i (\delta_Q - \mu)$$

and the intersection at infinity is:

$$(D_1 \cdot D_2)_\infty = -\sum_\sigma \sum_{P, Q} n_P m_Q \cdot g_\sigma(P, Q)$$

### 10.3 Faltings' Contributions

Gerd Faltings extended Arakelov's work significantly:

**Theorem (Faltings, 1984).** The arithmetic intersection pairing on an arithmetic surface satisfies an analogue of the Hodge index theorem: on the subspace of arithmetic divisors of degree 0 on the generic fiber, the pairing is **negative definite** (after accounting for the appropriate normalizations).

This was used in Faltings' proof of the Mordell conjecture (every curve of genus ≥ 2 over ℚ has finitely many rational points).

**Faltings' Riemann-Roch.** Faltings proved an arithmetic Riemann-Roch formula relating:
- The arithmetic degree deg(R π_* L̄) of the direct image of a metrized line bundle
- Arithmetic intersection numbers
- Analytic invariants (spectral zeta functions of Laplacians)

### 10.4 The Arithmetic Riemann-Roch Theorem of Gillet-Soulé

Henri Gillet and Christophe Soulé (1992) proved the definitive version:

**Theorem (Gillet-Soulé Arithmetic Riemann-Roch).** Let π: X → Spec(ℤ) be a smooth projective arithmetic variety of dimension d+1, and let L̄ = (L, h) be a Hermitian line bundle on X. Then:

$$\widehat{\mathrm{ch}}(\pi_! \bar{L}) = \pi_*\left(\widehat{\mathrm{ch}}(\bar{L}) \cdot \widehat{\mathrm{Td}}(\bar{T}_{X/\mathbb{Z}})\right) - a(\widetilde{T}(\bar{T}_{X/\mathbb{Z}(\mathbb{C})}))$$

where:
- ĉh is the arithmetic Chern character (in the arithmetic Chow ring ĈH*(X))
- T̂d is the arithmetic Todd class
- a is the map from differential forms to the arithmetic Chow ring
- T̃ is the R-genus (higher analytic torsion form)

### 10.5 Arithmetic Chow Groups

**Definition.** The **arithmetic Chow group** ĈH^p(X) of an arithmetic variety X consists of pairs (Z, g_Z) where:
- Z is a codimension-p algebraic cycle on X
- g_Z is a **Green current** for Z on X(ℂ): a (p-1, p-1)-current satisfying dd^c g_Z + δ_Z = ω for some smooth form ω

The arithmetic intersection product is:

$$(Z_1, g_1) \cdot (Z_2, g_2) = (Z_1 \cdot Z_2, g_1 * g_2)$$

where the Green current product g₁ * g₂ involves the star product of currents.

### 10.6 Connection to Zeta Functions

The arithmetic Riemann-Roch theorem connects to zeta functions through:

**(a) Heights.** The arithmetic degree map deg: ĈH^{d+1}(X) → ℝ computes **heights** of arithmetic cycles. Heights of points on varieties are intimately connected to Diophantine properties and L-functions.

**(b) Arithmetic ampleness.** An arithmetic line bundle L̄ is **arithmetically ample** if it has positive degree on every arithmetic curve. The arithmetic Nakai-Moishezon criterion (Zhang, 1995) characterizes this in terms of the positivity of the arithmetic intersection pairing.

**(c) The zeta function connection.** For a number field K with ring of integers 𝒪_K:
- The arithmetic degree of the metrized canonical bundle on Spec(𝒪_K) involves the discriminant Δ_K
- The Dedekind zeta function ζ_K(s) appears in the arithmetic Riemann-Roch theorem through its special values
- The leading term of ζ_K(s) at s = 0 is given by the analytic class number formula, which has an Arakelov-theoretic interpretation

**(d) The Hodge index theorem and RH.** As noted in Section 4.5, the arithmetic Hodge index theorem for the "curve" Spec(ℤ̄) = Spec(ℤ) ∪ {∞} would be directly equivalent to RH. Specifically:

The statement that the arithmetic intersection pairing on ĈH¹(Spec(ℤ̄)) restricted to divisors of degree 0 is negative definite **is equivalent to the non-vanishing of ζ(s) on Re(s) > 1/2** (and hence equivalent to RH by the functional equation).

### 10.7 Yuan-Zhang and Beyond

Xinyi Yuan and Shou-Wu Zhang have further developed arithmetic intersection theory:

- **Arithmetic Hodge index theorem for higher dimensions** (2017): Extending Faltings' result
- **Adelic line bundles**: A framework using adelic metrics that naturally incorporates all places
- **Height formulas**: Connecting heights of cycles to special values of L-functions (Gross-Zagier type formulas)

These developments show that arithmetic intersection theory is a powerful and active area, but the central challenge — proving the positivity that would imply RH — remains as daunting as ever.

---

## 11. What's Missing: Towards a Proof

### 11.1 The Rosetta Stone

The "Rosetta Stone" connecting number fields to function fields for RH would consist of:

| Ingredient | Function Field (KNOWN) | Number Field (NEEDED) |
|---|---|---|
| Base field | 𝔽_q | 𝔽₁ (unknown) |
| Curve | C/𝔽_q | Spec(ℤ̄) |
| Frobenius | Fr_q: x ↦ x^q | ??? (no known analogue) |
| Cohomology | H^i_ét(C, ℚ_ℓ) | H^i(Spec(ℤ̄), ???) |
| Trace formula | Lefschetz fixed-point | Weil explicit formula ✓ |
| Product surface | C ×_{𝔽_q} C | Spec(ℤ̄) ×_{𝔽₁} Spec(ℤ̄) |
| Intersection theory | Classical on surfaces | Arithmetic (Arakelov) ≈ |
| Positivity | Hodge index / C-S inequality | ≡ RH (circular!) |

### 11.2 What Mathematical Objects Would Need to Exist

For an algebraic/geometric proof of RH, one would likely need at least the following:

**(1) A cohomology theory for Spec(ℤ)** producing:
- H⁰ ≅ ℝ (or ℂ), one-dimensional
- H¹ = infinite-dimensional, with a "spectral" decomposition indexed by zeros of ζ(s)
- H² ≅ ℝ (or ℂ), one-dimensional
- A trace formula: ∑ eigenvalues on H¹ = ∑ contributions from primes
- Poincaré duality: H⁰ ↔ H²
- Functorial with respect to morphisms of arithmetic schemes

Deninger's program outlines the properties this should have. Connes' spectral realization achieves a form of this but not in a way that has yet proved RH.

**(2) A "Frobenius" operator or flow** Θ on H¹ such that:
- The spectrum of Θ (eigenvalues or spectral measure) encodes the zeros of ζ(s)
- Θ is "self-adjoint" in a suitable sense (which would force eigenvalues to be real, i.e., zeros on the critical line)
- The Weil explicit formula IS the trace formula for Θ

The self-adjointness of Θ is essentially equivalent to RH (the Hilbert-Pólya conjecture in a geometric guise).

**(3) A product formula / intersection theory** on Spec(ℤ̄) ×_{𝔽₁} Spec(ℤ̄):
- An arithmetic surface where we can do intersection theory
- A Hodge index theorem that isn't circularly equivalent to RH
- This requires breaking the circularity identified in Section 4.5

**(4) A positivity principle** that is:
- Provable from the structure of the space, not equivalent to RH
- Analogous to the Castelnuovo-Severi inequality
- Applicable in the infinite-dimensional setting

### 11.3 The Circularity Problem

The deepest challenge is the **circularity** of positivity:

In Weil's proof, positivity (Hodge index theorem) is a consequence of the underlying geometry of the surface C × C. It is provable from the axioms of algebraic geometry.

For Spec(ℤ), the analogous positivity **IS** the Riemann Hypothesis. To break this circularity, one would need to:
- Construct the "arithmetic surface" in a way where positivity follows from its definition, not from the properties of ζ(s)
- Find a NEW source of positivity not currently known in arithmetic geometry
- Possibly use non-commutative geometry, where positivity can arise from different mechanisms (e.g., the positivity of states in C*-algebras)

Connes' program specifically attempts to address this through:
- The trace formula on the adele class space
- The positivity of the Weil distribution
- C*-algebraic positivity principles from noncommutative geometry

But even Connes has not succeeded in proving the required positivity.

### 11.4 Potential Breakthrough Directions

Based on the current state of research, the most promising directions for an algebraic proof include:

**(A) Complete the Connes-Consani arithmetic site program:**
- Develop the sheaf cohomology of the scaling site to the point where a Riemann-Roch theorem can be proved
- Establish the correct analogue of the Castelnuovo-Severi inequality in this tropical/semiring setting
- The tropical Riemann-Roch theorem (Baker-Norine) provides a template

**(B) Find the right notion of "absolute cohomology":**
- A synthesis of Deninger's infinite-dimensional cohomology with Connes' noncommutative geometry
- Should simultaneously be:
  - Concrete enough to compute with
  - Abstract enough to apply to Spec(ℤ)
  - Equipped with a self-adjoint "Frobenius" operator

**(C) Develop Borger's Λ-ring approach to full maturity:**
- Use the Λ-structure (simultaneous Frobenius lifts) as the replacement for Frobenius
- Develop intersection theory for Λ-schemes
- The descent from ℤ to 𝔽₁ via Λ-structures should produce the right cohomology

**(D) Connect to the Langlands program:**
- The Langlands program provides automorphic representations that encode L-functions
- Automorphic forms on GL(1) give ζ(s)
- The geometric Langlands program (over function fields) uses the full machinery of algebraic geometry
- Transporting geometric Langlands to the number field setting could provide the missing link

**(E) Exploit the analogy with physics:**
- Statistical mechanics on the Bost-Connes system
- The partition function of the Bost-Connes system exhibits a phase transition at β = 1, related to the pole of ζ(s)
- Quantum statistical mechanics may provide new positivity principles

### 11.5 Assessment of Each Program

| Program | Maturity | Main Achievement | Main Gap |
|---|---|---|---|
| Deninger | Conceptual | Correct framework shape | No concrete construction |
| Connes | High | Spectral realization, trace formula | Positivity unproved |
| 𝔽₁ (various) | Medium | Multiple frameworks | No unified theory, no cohomology |
| Borger Λ-rings | Medium | Clean algebraic framework | No intersection theory |
| Arakelov/Arithmetic RR | High | Full intersection theory | Positivity ≡ RH (circular) |
| Arithmetic site (C-C) | Growing | Topos-theoretic framework | Riemann-Roch incomplete |

### 11.6 The Meta-Question

Perhaps the deepest lesson from the algebraic approaches is this: **the Riemann Hypothesis may require mathematics that does not yet exist**. Every existing framework, when pushed far enough, either:
1. Becomes circular (positivity ≡ RH), or
2. Encounters an impassable technical barrier, or
3. Is too underdeveloped to apply

The history of the function field case is instructive: Weil's proof required developing entirely new foundations for algebraic geometry (algebraic geometry over arbitrary fields, intersection theory on surfaces without assuming the ground field is algebraically closed). The classical RH may similarly require a foundational revolution — perhaps in the direction of absolute geometry, perhaps in a direction not yet imagined.

### 11.7 Concrete Desiderata

For a complete algebraic proof, we would need to verify the following chain:

1. **Construct** a cohomology theory H^i(Spec(ℤ̄)) with the right dimensions and functorial properties
2. **Identify** an operator Θ on H¹ whose spectrum is the set of zeros of ζ(s)
3. **Prove** a trace formula relating Tr(f(Θ)) to a sum over primes (= Weil explicit formula)
4. **Establish** a positivity principle for Θ (self-adjointness, or a Hodge-index-type inequality) from STRUCTURAL properties of the cohomology, not from properties of ζ(s)
5. **Deduce** that all eigenvalues of Θ have Re(ρ) = 1/2

Steps 1-3 are partially achieved by Connes and Deninger (in different frameworks). Step 4 is the core unsolved problem. Step 5 follows formally from Step 4.

The search for Step 4 — a source of positivity that isn't circular — is the holy grail of algebraic approaches to the Riemann Hypothesis.

---

## References and Further Reading

### Primary Sources

- A. Weil, *Sur les courbes algébriques et les variétés qui s'en déduisent* (1948)
- A. Weil, "Sur les 'formules explicites' de la théorie des nombres premiers" (1952)
- P. Deligne, "La conjecture de Weil. I" (1974), *Publ. Math. IHÉS* 43
- P. Deligne, "La conjecture de Weil. II" (1980), *Publ. Math. IHÉS* 52
- A. Grothendieck, *SGA 4, 4½, 5* (étale cohomology)

### Deninger's Program
- C. Deninger, "On the Γ-factors attached to motives" (1991), *Invent. Math.*
- C. Deninger, "Motivic L-functions and regularized determinants" (1994)
- C. Deninger, "Some analogies between number theory and dynamical systems on foliated spaces" (1996)
- C. Deninger, "A note on arithmetic topology and dynamical systems" (2002)

### Connes' Approach
- A. Connes, "Trace formula in noncommutative geometry and the zeros of the Riemann zeta function" (1999)
- A. Connes, *Noncommutative Geometry* (Academic Press, 1994)
- A. Connes and C. Consani, "The Arithmetic Site" (2014), *C. R. Math.*
- A. Connes and C. Consani, "Geometry of the Scaling Site" (2017)
- A. Connes and C. Consani, "The scaling site" (2018), *Advances in Math.*

### 𝔽₁ Geometry
- J. Tits, "Sur les analogues algébriques des groupes semi-simples complexes" (1957)
- Y. Manin, "Lectures on zeta functions and motives" (1995)
- C. Soulé, "Les variétés sur le corps à un élément" (2004)
- J. Borger, "Λ-rings and the field with one element" (2009)
- A. Connes and C. Consani, "Schemes over 𝔽₁ and zeta functions" (2010)
- O. Lorscheid, "Blueprints — towards absolute arithmetic?" (2012)

### Arakelov Geometry
- S. Arakelov, "Intersection theory of divisors on an arithmetic surface" (1974)
- G. Faltings, "Calculus on arithmetic surfaces" (1984)
- H. Gillet and C. Soulé, "An arithmetic Riemann-Roch theorem" (1992)
- S.-W. Zhang, "Small points and adelic metrics" (1995)
- X. Yuan and S.-W. Zhang, "The arithmetic Hodge index theorem for adelic line bundles" (2017)

### Absolute Geometry
- N. Durov, "New approach to Arakelov geometry" (2007)
- S. Haran, *Non-Additive Geometry* (2007)
- B. Toen and M. Vaquié, "Au-dessous de Spec ℤ" (2009)

### IUT Theory
- S. Mochizuki, "Inter-universal Teichmüller Theory I-IV" (2012-2021)
- P. Scholze and J. Stix, "Why abc is still a conjecture" (2018)

### Surveys
- J.-P. Serre, "Analogues kählériens de certaines conjectures de Weil" (1960)
- N. Katz, "An overview of Deligne's proof of the Riemann hypothesis for varieties over finite fields" (1976)
- E. Bombieri, "The Rosetta Stone" in *The Millennium Prize Problems* (2006)
- A. Connes, "An essay on the Riemann Hypothesis" (2015), in *Open Problems in Mathematics*
- J. López Peña and O. Lorscheid, "Mapping 𝔽₁-land" (2011)

---

*Document prepared as part of the Riemann Hypothesis research project.*
*This survey focuses on the algebraic and arithmetic-geometric approaches; see companion documents for analytic, spectral-theoretic, and computational approaches.*
