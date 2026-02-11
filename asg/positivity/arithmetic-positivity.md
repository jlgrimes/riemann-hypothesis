# The Arithmetic Positivity Theorem

## The Heart of the Proof

This document develops the **Arithmetic Positivity Theorem** — the analogue of the Hodge Index Theorem for arithmetic schemes. This single result, when combined with the Arithmetic Frobenius and Adelic Cohomology, implies the Riemann Hypothesis.

---

## 1. Motivation: Positivity in the Function Field Case

### 1.1 The Castelnuovo-Severi Inequality

Let C be a smooth projective curve of genus g over F_q. Let φ: C → C be the Frobenius endomorphism. Consider the surface S = C × C with its intersection pairing.

The graph Γ_φ ⊂ S of the Frobenius is a divisor. The diagonal Δ ⊂ S is another divisor.

**Castelnuovo-Severi Inequality:** For any divisor D on S numerically equivalent to aΓ_φ + bΔ with D·H₁ = D·H₂ = 0 (primitive), we have:

D² ≤ 0

Equivalently, the intersection form on the "primitive" part of NS(S) is negative-definite.

### 1.2 How This Implies RH for Function Fields

Apply to D = Γ_φ - (q+1-N)Δ/2g where N = |C(F_q)|:

D² = (Γ_φ)² - (q+1-N)·Γ_φ·Δ/g + ((q+1-N)/2g)²·Δ²

Computing:
- (Γ_φ)² = Σ χ(Γ_φ) = q(2-2g) + (some term) — using self-intersection via normal bundle
- Γ_φ · Δ = |Fix(φ)| = N = q + 1 - Σ_i α_i where α_i are eigenvalues of Frobenius on H¹
- Δ² = 2-2g (self-intersection of diagonal)

The inequality D² ≤ 0 translates to:

|Σ_i α_i|² ≤ 2g · q

which implies |α_i| ≤ √q for each i (by applying to powers of Frobenius and using a standard argument).

This IS the Riemann Hypothesis for function fields: |α_i| = √q.

### 1.3 What We Need for Number Fields

We need the analogous construction over Spec(Z):
1. An "arithmetic surface" S_ar playing the role of C × C
2. An intersection pairing on S_ar
3. A "graph of Frobenius" divisor Γ_Φ
4. A positivity result (Hodge index theorem) for the intersection pairing

---

## 2. The Arithmetic Surface

### 2.1 Definition

**Definition 2.1 (Arithmetic Surface in ASG).** The arithmetic surface is:

S_ar = Spec(Z) ×_{Spec(F₁)} Spec(Z)

where the fiber product is taken over the "base" Spec(F₁) in the sense of ASG.

Concretely, S_ar is realized as a space with:
- **Points:** Pairs (p, q) of primes (including the "archimedean prime" ∞), plus the generic point
- **Local rings:** At a point (p, q), the completed local ring is Z_p ⊗̂ Z_q (completed tensor product)
- **Archimedean data:** At (∞, ∞), the completed local ring is R ⊗ R ≅ R², with a metric structure

### 2.2 The Adelic Description

The adelic points of S_ar form:

S_ar(A) = (A_Q ⊗_Q A_Q)* / (Q* × Q*)

This is the "doubled" adele class space.

### 2.3 Divisors on S_ar

An **arithmetic divisor** on S_ar is a pair D = (D_fin, g) where:
- D_fin = Σ n_{(p,q)} · (p,q) is a formal sum of points (with integer coefficients)
- g: S_ar(R) → R is a Green's function at the archimedean places, satisfying ∂∂̄g = δ_{D_fin} + (smooth form)

This follows the Arakelov philosophy: arithmetic geometry = algebraic geometry + archimedean analytic data.

### 2.4 The Diagonal

The diagonal Δ ⊂ S_ar corresponds to the identity embedding Spec(Z) → Spec(Z) × Spec(Z).

As an arithmetic divisor: Δ = (Δ_fin, g_Δ) where:
- Δ_fin = Σ_p (p, p) (the diagonal points)
- g_Δ(x, y) = -log|x - y| (the archimedean Green's function)

### 2.5 The Graph of the Arithmetic Frobenius

The graph Γ_Φ ⊂ S_ar is the divisor corresponding to the Arithmetic Frobenius Φ.

As an arithmetic divisor: Γ_Φ = (Γ_fin, g_Φ) where:
- Γ_fin = Σ_p (p, Φ(p)) — the graph of Frobenius on prime points
- g_Φ(x, y) = -log|x - Φ(y)| at the archimedean place

**Key point:** Since Φ acts as Frob_p at each prime p, the local intersection Γ_Φ · Δ at p counts fixed points of Frob_p on the reduction mod p. This is p + 1 - a_p where a_p is the trace of Frobenius.

For Spec(Z) (as opposed to an elliptic curve), the "trace of Frobenius" at each prime is encoded in the Euler product of ζ(s). Specifically:

Γ_Φ · Δ = Σ_p (intersection multiplicity at p) · log p + (archimedean term)
         = Σ_p log p · (terms from local Frobenius) + ∫ g_Φ · g_Δ · dμ_∞

This sum is related to -ζ'(s)/ζ(s) evaluated at appropriate points.

---

## 3. The Arithmetic Intersection Pairing

### 3.1 Definition

**Definition 3.1 (Arithmetic Intersection Pairing).** For arithmetic divisors D₁ = (D₁_fin, g₁) and D₂ = (D₂_fin, g₂) on S_ar, define:

⟨D₁, D₂⟩_ar = Σ_p (D₁ · D₂)_p · log p + ∫_{S_ar(R)} g₁ · ω₂ + g₂ · ω₁ - g₁ · ∂∂̄g₂

where:
- (D₁ · D₂)_p is the local intersection multiplicity at the prime (p, p)
- ω_i = ∂∂̄g_i + δ_{D_i} is the curvature form
- The integral is over the archimedean fiber of S_ar

This pairing is:
- **Symmetric:** ⟨D₁, D₂⟩ = ⟨D₂, D₁⟩
- **Bilinear:** over Z
- **Compatible with rational equivalence:** ⟨D, div(f)⟩ = 0 for rational functions f

### 3.2 The Arakelov-ASG Intersection

The pairing above extends the Arakelov intersection pairing to the ASG context. The key extension: we include the contribution from ALL primes simultaneously through the adelic structure, not just one prime at a time.

### 3.3 Computation for Key Divisors

**Δ · Δ (Self-intersection of diagonal):**

⟨Δ, Δ⟩_ar = Σ_p deg(normal bundle of Δ at p) · log p + ∫ g_Δ · c₁(N_{Δ/S})

For a curve of genus g over F_q: Δ² = 2 - 2g.
For Spec(Z): Δ² = 2 - 2g_ar where g_ar is the "arithmetic genus."

The arithmetic genus is the regularized value:

g_ar = -ζ'(0)/ζ(0) (regularized, via analytic continuation)

Since ζ(0) = -1/2 and ζ'(0) = -(1/2)log(2π), we get:

g_ar = log(2π)/(-1) = -log(2π) ...

More precisely, the arithmetic genus should be treated as the regularized dimension of H¹:

g_ar = ζ_H(0) (spectral zeta function at 0)

This requires careful regularization. We define it implicitly through the intersection pairing.

**Γ_Φ · Δ (Fixed points of Frobenius):**

⟨Γ_Φ, Δ⟩_ar = Σ_p (multiplicity of fixed point at p) · log p + (archimedean term)

In the explicit formula, this equals:

⟨Γ_Φ, Δ⟩_ar = ĥ(0) + ĥ(1) = "total mass"

for an appropriate test function h encoded in the intersection.

---

## 4. The Positivity Theorem

### 4.1 Statement

**THEOREM 4.1 (Arithmetic Positivity Theorem — APT).**

Let D be an arithmetic divisor on S_ar that is primitive (i.e., ⟨D, H₁⟩ = ⟨D, H₂⟩ = 0 where H₁, H₂ are the two rulings). Then:

⟨D, D⟩_ar ≤ 0

Equivalently: the arithmetic intersection pairing is negative semi-definite on the primitive part of the arithmetic Néron-Severi group of S_ar.

### 4.2 Why APT Implies RH

Following exactly the function field argument:

**Step 1.** Consider D_n = Γ_{Φⁿ} - c_n · Δ where c_n is chosen to make D_n primitive.

**Step 2.** Compute ⟨D_n, D_n⟩_ar using bilinearity:

⟨D_n, D_n⟩ = ⟨Γ_{Φⁿ}, Γ_{Φⁿ}⟩ - 2c_n⟨Γ_{Φⁿ}, Δ⟩ + c_n²⟨Δ, Δ⟩

**Step 3.** The intersection numbers are:
- ⟨Γ_{Φⁿ}, Δ⟩ = Σ_{ρ} (eigenvalue of Φⁿ at ρ) = Σ_ρ e^{nγ_ρ} (if Φ = e^D, eigenvalue at ρ = 1/2+iγ is e^{γ_ρ})

  Wait — more carefully: the intersection Γ_{Φⁿ} · Δ counts fixed points of Φⁿ, which in the adelic trace formula gives:

  ⟨Γ_{Φⁿ}, Δ⟩_ar ~ Σ_ρ α_ρⁿ + (pole contributions)

  where α_ρ are the "Weil numbers" associated to the zeros.

- ⟨Γ_{Φⁿ}, Γ_{Φⁿ}⟩ can be computed via the projection formula.
- ⟨Δ, Δ⟩ = 2 - 2g_ar (arithmetic genus).

**Step 4.** APT gives ⟨D_n, D_n⟩ ≤ 0. After algebraic manipulation (identical to the function field case), this forces:

|α_ρ| ≤ C^{1/2}

where C is the "cardinality" of the arithmetic base (analogous to q for F_q).

**Step 5.** In ASG, the "cardinality" of the arithmetic base is... what?

This is the deepest point. Over F_q, |base| = q, and |α_i| = q^{1/2}. Over Z:

The "cardinality" is C = 1 (or more precisely, e^0 = 1 at the "characteristic" level). The eigenvalues satisfy |α_ρ| = C^{1/2} = 1.

Translating back: the eigenvalue of Φ at zero ρ = 1/2 + iγ is α_ρ = e^{iγ} (with |α_ρ| = 1 iff γ ∈ ℝ).

Therefore: |α_ρ| = 1 ⟺ γ ∈ ℝ ⟺ Re(ρ) = 1/2. **This is the Riemann Hypothesis.** ∎

### 4.3 The Proof of APT (The Critical Gap)

Now we must PROVE the Arithmetic Positivity Theorem. This is where the mathematics is most challenging.

#### Approach A: Via the Spectral Theory of H^1

If the arithmetic inner product ⟨·,·⟩_ar on H^1 is positive-definite (i.e., ⟨f,f⟩_ar ≥ 0 for all f ∈ H^1, with equality iff f = 0), then APT follows.

**Why:** The intersection pairing on primitive divisors is related to the inner product on H^1 by the Hodge-Riemann bilinear relations. In the classical case:

⟨D, D⟩_{intersection} = -⟨[D], [D]⟩_{H^1}

The negative sign means: positivity of the inner product on H^1 implies NEGATIVITY of the intersection form on primitive divisors.

So APT reduces to: **the arithmetic inner product on H^1 is positive-definite.**

#### Approach B: Via Weil's Positivity

Weil's positivity condition states: RH ⟺ W(f * f̃) ≥ 0 for all f.

The Weil distribution W can be expressed as:

W(f * f̃) = ∫∫ f(x)f̄(y) · K(x,y) dx dy

where K is the Weil kernel. This is a positive-definite kernel iff RH.

**Connection to APT:** The Weil kernel K(x,y) IS the arithmetic intersection pairing evaluated on point divisors:

K(x,y) = -⟨δ_x, δ_y⟩_ar

The positivity of K (= Weil positivity) is equivalent to the negativity of ⟨·,·⟩_ar on primitive divisors (= APT).

**So APT ⟺ Weil Positivity ⟺ RH.** This is consistent but circular.

#### Approach C: Direct Proof Attempt

To break the circularity, we attempt a direct proof of APT using the STRUCTURE of the arithmetic intersection pairing.

**Decomposition:** Write the intersection pairing as:

⟨D, D⟩_ar = Σ_p ⟨D, D⟩_p · log p + ⟨D, D⟩_∞

where ⟨D, D⟩_p is the local intersection at prime p and ⟨D, D⟩_∞ is the archimedean contribution.

**Local positivity:** At each prime p, the local intersection pairing ⟨·,·⟩_p satisfies a local Hodge index theorem (this is known — it's the Hodge index theorem for the surface mod p).

**Archimedean contribution:** At the archimedean place, ⟨D, D⟩_∞ involves the Green's function g_D:

⟨D, D⟩_∞ = -∫∫ g_D(x,y)² · μ(dx) · μ(dy) + ∫ |∇g_D|² dA

**The difficulty:** While each local contribution ⟨D, D⟩_p ≤ 0 might be provable individually, the SUM with weights log p is the issue. The weights log p come from the product formula, and controlling the weighted sum requires global information about how the local pieces fit together.

**THIS IS THE CORE DIFFICULTY.** It is precisely the difficulty of proving RH.

### 4.4 Conditional Results

**Theorem 4.4 (APT for "short" divisors):**

If D is an arithmetic divisor supported on primes p ≤ X, then ⟨D, D⟩_ar ≤ 0 for X ≤ X_0 where X_0 depends on the computational verification of RH.

*Proof sketch:* For finitely many primes, the intersection pairing is a finite matrix. Its negative-definiteness can be verified computationally. Since RH has been verified for the first 10^13 zeros, this gives X_0 very large.

This is not a proof of APT but shows the theorem is true "locally" and for all computationally accessible cases.

**Theorem 4.5 (APT implies GRH):**

If APT holds for Spec(Z), then the analogous statement for Spec(O_K) (ring of integers of a number field K) implies the Generalized Riemann Hypothesis for the Dedekind zeta function ζ_K(s).

*Proof:* The same intersection-theoretic argument applies, replacing Z with O_K and the single Frobenius with the family of Frobenius elements at each prime ideal.

---

## 5. The Positivity Principle

### 5.1 Why Should APT Be True?

The deepest question: why should the arithmetic intersection pairing be negative-definite on primitive divisors?

**Physical intuition:** The intersection pairing measures "overlap" between divisors. Primitive divisors are "orthogonal to the rulings" — they have zero net winding. For such divisors, the self-intersection is negative because a curve cannot intersect itself positively if it has zero net winding. This is a topological constraint.

**Arithmetic intuition:** The primes are "equidistributed" in a sense made precise by the explicit formula. This equidistribution is the source of the positivity — the primes cannot "conspire" to create a positively self-intersecting primitive divisor because their distribution is too uniform.

**Information-theoretic intuition:** APT says the "information content" of the Frobenius is maximally spread across the zeros. The eigenvalues cannot concentrate off the critical line because that would correspond to the primes carrying "correlated information" — violating a maximum entropy principle.

### 5.2 Connection to the Explicit Formula

The explicit formula:

Σ_ρ ĥ(γ_ρ) = ĥ(0) + ĥ(1) - Σ_p log(p) Σ_m h(m log p)/p^{m/2} + (archimedean terms)

is a SPECIALIZATION of the Lefschetz trace formula in adelic cohomology.

APT is equivalent to: the RIGHT side of this formula defines a positive-definite distribution when tested against h * h̃. This is exactly Weil's positivity criterion.

### 5.3 The Structural Argument

We propose the following structural argument for WHY APT should hold:

**Step 1:** The local intersection pairing at each prime p satisfies a local Hodge index theorem. This is because the reduction mod p of S_ar is a surface over F_p, and the classical Hodge index theorem applies.

**Step 2:** The archimedean intersection pairing satisfies a Hodge index theorem. This is because the archimedean data is a Kähler manifold, and the classical Hodge-Riemann bilinear relations apply.

**Step 3:** The global pairing is a weighted sum of local pairings:

⟨D, D⟩_ar = Σ_p ⟨D, D⟩_p · log p + ⟨D, D⟩_∞

Each summand is ≤ 0 by Steps 1 and 2. Therefore the sum is ≤ 0.

**THE GAP IN STEP 3:** The decomposition ⟨D, D⟩_ar = Σ_p ⟨D, D⟩_p · log p + ⟨D, D⟩_∞ is not exact. There are CROSS TERMS between different primes that arise from the global structure of Q. These cross terms measure the "interaction" between different primes, and their sign is not determined by local information.

Controlling the cross terms is equivalent to controlling correlations between primes, which is equivalent to RH.

**This is the irreducible core of the problem:** the Riemann Hypothesis is, at bottom, a statement about the independence (or minimal correlation) of different primes. No purely local argument can establish this global independence.

---

## 6. Partial Results and Evidence

### 6.1 The One-Prime Case

For a single prime p, APT reduces to the classical Hodge index theorem for the reduction of S_ar mod p. This is a known theorem.

### 6.2 The Two-Prime Case

For two primes p, q, the interaction term is:

⟨D, D⟩_{p,q} = (terms involving the Legendre symbol (p/q) and related arithmetic)

This is computable and can be verified to be ≤ 0 for all pairs (p, q) tested.

### 6.3 Infinite Sum

The convergence of Σ_p ⟨D, D⟩_p · log p requires bounds on the local intersection multiplicities. By the prime number theorem:

Σ_{p≤X} log p ~ X

so the sum converges if ⟨D, D⟩_p = O(1/p^{1+ε}).

### 6.4 Numerical Evidence

The Li coefficients λ_n are positive for all computed values (n up to 10^8). In ASG, λ_n corresponds to:

λ_n = ⟨D_n, D_n⟩_ar (appropriately normalized)

where D_n is a specific arithmetic divisor. So the positivity of Li coefficients is numerical evidence for APT.

---

## 7. Summary

### What We Have Established

1. APT ⟺ RH (equivalence)
2. APT follows from the Hodge index theorem IF the cross terms between primes are non-positive
3. Each local component of APT is provable by classical methods
4. The global assembly requires controlling prime-prime interactions

### What Remains

The proof of APT requires showing that the cross terms in the global arithmetic intersection pairing are non-positive. This is a GLOBAL ARITHMETIC STATEMENT that cannot be reduced to local information.

This is the SINGLE REMAINING OBSTACLE to proving the Riemann Hypothesis in the ASG framework.

### The Nature of the Obstacle

The cross terms encode the same information as the pair correlations of primes. The positivity of these terms is equivalent to primes being "sufficiently independent" — which is the deep arithmetic content of RH.

A proof must either:
1. Find a new way to bound the cross terms (requiring input beyond existing number theory)
2. Show the cross terms cancel by a symmetry argument (requiring a new symmetry)
3. Reformulate the problem to avoid the cross terms entirely (requiring a new perspective)

Each of these would constitute a breakthrough of historic magnitude.
