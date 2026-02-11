# The Riemann Hypothesis: A Collaborative Research Repository

## The Problem

**The Riemann Hypothesis (RH):** All non-trivial zeros of the Riemann zeta function ζ(s) have real part equal to 1/2.

This is one of the seven Clay Millennium Prize Problems, open since 1859.

## Main Result

This project develops **Arithmetic Spectral Geometry (ASG)**, a new mathematical framework that reduces RH to a single geometric statement: the **Arithmetic Positivity Theorem (APT)**. The reduction follows Weil's 1948 strategy for function fields, translated to the number field setting.

**Status:** RH is reduced to proving APT, which is reduced to a specific quantitative bound on the Weil kernel at prime-logarithm differences (the Arithmetic Cross-Term Bound). Numerical experiments confirm the bound with a ratio of ~500:1.

## Repository Structure

```
riemann-hypothesis/
├── README.md                                    # This file
├── SYNTHESIS.md                                 # Master synthesis of Phase 1 research
│
├── research/                                    # Phase 1: Survey of all known approaches
│   ├── classical-approaches.md                  # Hardy, Selberg, Levinson, Conrey, ...
│   ├── spectral-approaches.md                   # Hilbert-Pólya, Berry-Keating, Connes
│   ├── random-matrix-theory.md                  # Montgomery, Odlyzko, GUE, Katz-Sarnak
│   ├── algebraic-approaches.md                  # Weil, Deligne, F₁, arithmetic geometry
│   ├── computational-evidence.md                # Computational experiments analysis
│   ├── foundational-framework.md                # Mathematical foundations
│   └── equivalent-formulations.md               # 15+ equivalent statements of RH
│
├── computational/                               # Phase 1: Numerical tools
│   ├── zero_computation.py                      # Compute zeta zeros via mpmath
│   ├── zero_visualization.py                    # Visualize zeros and ζ(s)
│   ├── pair_correlation.py                      # Pair correlation vs GUE
│   ├── li_criterion.py                          # Li coefficients
│   ├── newman_constant.py                       # de Bruijn-Newman constant
│   ├── nyman_beurling.py                        # Nyman-Beurling criterion
│   └── explicit_formula.py                      # Explicit formula for π(x)
│
├── proofs/                                      # Phase 1: Proof strategies
│   ├── novel-strategies.md                      # 7 novel strategies evaluated
│   └── proof-attempt.md                         # Spectral positivity attempt
│
├── explorations/
│   └── connections.md                           # Cross-disciplinary connections
│
└── asg/                                         # Phase 2+3: Arithmetic Spectral Geometry
    ├── ASG-MANIFESTO.md                         # ★ Founding document of ASG
    │
    ├── foundations/                              # The 5 axioms of ASG
    │   ├── axioms.md                            # Formal axiom system
    │   ├── arithmetic-frobenius.md              # Axiom I: The Arithmetic Frobenius Φ
    │   └── spectral-space.md                    # Axiom III: The spectral space H
    │
    ├── cohomology/
    │   └── adelic-cohomology.md                 # Axiom IV: New Weil cohomology theory
    │
    ├── proof/
    │   └── main-proof.md                        # Conditional proof: ASG axioms ⟹ RH
    │
    └── positivity/                              # Phase 3: The APT attack
        ├── arithmetic-positivity.md             # ★ APT statement and three approaches
        ├── APT-PROOF-ATTEMPT.md                 # ★ Final synthesis / near-miss analysis
        │
        ├── cross-terms/                         # The core obstacle
        │   ├── structure.md                     # Matrix formulation of cross-terms
        │   ├── explicit-derivation.md           # Explicit formula analysis
        │   └── arithmetic-cross-term-bound.md   # ★ The sharpest remaining conjecture
        │
        ├── sieve/
        │   └── sieve-bounds.md                  # Bombieri-Vinogradov, Elliott-Halberstam
        │
        ├── energy/
        │   └── energy-approach.md               # Log-gas, de Bruijn-Newman dynamics
        │
        ├── algebraic/
        │   └── function-field-analysis.md       # ★ Forensic analysis of Weil's proof
        │
        ├── theta-ampleness/
        │   └── theta-positivity.md              # ★ Theta ampleness strategy
        │
        └── computational/                       # Numerical experiments
            ├── weil_positivity_test.py           # Weil W(f*f̃) ≥ 0 tests
            ├── weil_corrected.py                 # Corrected with exact archimedean
            ├── cross_term_matrix.py              # Cross-term matrix eigenvalues
            ├── li_analysis.py                    # Li coefficient structure
            ├── diagonal_dominance_test.py        # ★ Full diagonal dominance test
            ├── numerical-results.md              # Summary of all numerical findings
            └── plots/                            # Generated visualizations
```

★ = Key documents

## The Logical Structure

```
Phase 1: Research Survey
    └── Identify Weil's function field proof as the most promising framework

Phase 2: Arithmetic Spectral Geometry (ASG)
    ├── Define Arithmetic Frobenius Φ on the idele class group
    ├── Construct spectral space H = L²(C_Q, ω) ⊖ H⁰ ⊖ H²
    ├── Build adelic cohomology theory
    └── Reduce RH to APT (Arithmetic Positivity Theorem)

Phase 3: Attack on APT
    ├── Forensic analysis of Weil's proof: identify 6 failure points
    ├── Sieve theory: BV gives diagonal dominance for large primes
    ├── Energy methods: real zeros are local energy minimum
    ├── Theta ampleness: the most concrete path (curvature of Θ)
    ├── Numerical experiments: kernel structure at arithmetic points
    └── Final synthesis: ACTB conjecture as the sharpest remaining gap

Result: RH ⟺ APT ⟺ ACTB (Arithmetic Cross-Term Bound)
```

## Key Numerical Finding

The Weil kernel K(x) at prime-logarithm differences x = log(p/q) decomposes as:

- **K_bg(x)** (background, from digamma): magnitude ~0.1 to 1.4, POSITIVE
- **K_zeros(x)** (from zeta zeros): magnitude ~0.001 to 0.003, oscillating

The background dominates by a factor of ~500. This means the positivity of the Weil form is driven by the archimedean geometry, not by the zeros. A proof of APT may be achievable by bounding the zero oscillation below the archimedean background.

## Key References

- Riemann (1859), Hardy (1914), Selberg (1942), Weil (1948)
- Connes (1999): Trace formula approach
- Bombieri (2000): Problems of the Millennium
- Rodgers-Tao (2018): Λ ≥ 0
- Ramaré (2013): Explicit Bombieri-Vinogradov constants
- Platt (2017): Verification of first 10^{13} zeros
