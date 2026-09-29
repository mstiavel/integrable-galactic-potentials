# Non-classical integrable potentials for ellipticals: session research note

*Working note, 26–27 September 2026. Covers the analytic and numerical work done in this session on whether a smooth axisymmetric potential with a spherical (cusped or cored) centre and an outward-growing flattening can carry a polynomial third integral.*

---

## 1. The problem and the framing

**Goal.** Find integrable potentials suitable for real elliptical galaxies: centrally concentrated, spherical (or nearly so) at the centre, flattened or triaxial in the envelope, with a smooth transition. Stäckel potentials are excluded from the outset (they cannot host an isolated central cusp), and the aim is a genuinely new family, not a repaired Stäckel one.

**Framing established early.**

- *Why Stäckel fails, geometrically.* For a natural Hamiltonian, a second integral quadratic in momenta is equivalent to orthogonal separability (Eisenhart/Benenti). Singularities of a Stäckel potential live on the singular locus of the coordinate system: focal ellipse/hyperbola in ellipsoidal coordinates, focal ring (oblate) or focal segment (prolate) in spheroidal ones. An isolated point singularity at the origin is available only in the spherical family and in parabolic coordinates (Sridhar & Touma's cuspy but lopsided models). So "spherical cusp + non-spherical envelope" is excluded by the coordinate systems, not by any subtlety of the density functions. An oblate cusped centre is excluded for the same reason (focal ring); only a spherical centre is available from Stäckel pieces.
- *The two proposed ideas are one idea.* Non-quadratic integrals ⇔ non-separable integrability. The atlas picture (action–angle charts per orbit family) is the ordinary Liouville–Arnold situation; the sharper form of it is an integral that is smooth but not analytic, glued family by family. Bolsinov & Taimanov's geodesic flows with a C^∞ integral and positive topological entropy show this is mathematically legitimate.
- *Local abundance, global rarity.* For polynomial integrals of degree n ≥ 3 the direct-method PDEs have functional freedom locally (Bialy–Mironov, hydrodynamic-type reformulation); the obstruction is global. On a ball with a cusp this global problem is essentially unstudied.
- *The centre fixes the seed.* With a spherical centre, any polynomial integral must reduce there to a function of E, L² and L_z (axisymmetric), which is reducible; the irreducible part of the integral is born from the departure from sphericity. This motivates a perturbative expansion in the shape function at fixed polynomial degree.

**Standing caveats for everything below.** All tests use the two invariant rectilinear orbits of an axisymmetric reflection-symmetric potential (the symmetry axis and the equatorial radial orbit with L_z = 0). Even-in-transverse integrals only. Cusp results are local at the centre (hence independent of the envelope); core results are perturbative in the flattening.

---

## 2. Tool 1: Morales–Ramis via numerical monodromy

**NVE along the invariant orbits.** For Φ = A(r) + B(r)P₂(cos θ), with z as independent variable along the orbit:

    2(E − U) ξ'' − U' ξ' + (U'/z − 3σB/z²) ξ = 0

with (U, σ) = (A + B, +1) on the axis and (A − B/2, −1) in the plane. Verified against the Cartesian Hessian. The shape enters only through B/r² (the ellipticity function); at B = 0 the exact solution ξ = z (tilted radial orbit) exists, hence the spherical case has a unipotent Galois group.

**Method.** For rational A, B the equation is Fuchsian, so the Galois group is the Zariski closure of the monodromy group (Schlesinger). The fundamental matrix is propagated numerically around every singular point in the complex z-plane. Two facts together force G⁰ = SL(2,ℂ): an element of infinite order (|tr| > 2 in SL₂ normalization) and a non-abelian derived subgroup (kills both the dihedral and the unipotent-Borel cases). Controls: spherical case gives commutators at 10⁻¹¹–10⁻¹³; the product of all finite loops reproduces the unipotent monodromy at infinity predicted by the exponents (0, −1) to six decimals in every run.

**Results.**
- Hernquist cusp A = −1/(1+r), B = εr²/(1+r)⁴: G⁰ = SL(2,ℂ) on both orbits for ε ∈ {0.3, 1}, E ∈ {−0.5, −0.2}. The obstruction is O(1) already at ε = 10⁻³: the Galois group jumps for any ε ≠ 0 (as in Hénon–Heiles).
- Cusp-to-envelope profiles realizing the endpoint-admissible flattenings (family F: Hernquist inside, Φ ∝ r outside, B = εr²/(1+r); family G: cored inside): G⁰ = SL(2,ℂ) at the candidate values ε = 0.25, −0.2 and at five non-candidate values, both orbits.

**Interpretation.** No additional integral meromorphic on a finite covering exists for these profiles. Morales–Ramis says nothing against smooth non-analytic integrals, which is exactly the atlas class.

---

## 3. Tool 2: the endpoint (scale-free envelope) condition

At a potential-dominated endpoint A ~ a rᵘ, B ~ b rᵘ (r → ∞ for μ > 0), the NVE exponent difference is δ = Δ/2 with

    Δ² = (μ+2)² − 24 b_eff/u_eff,   (u_eff, b_eff) = (a+b, b) on the axis, (a − b/2, −b) in the plane.

Writing Δ² = (μ−2)² + 8μλ with λ = 1 − 3b_eff/(μu_eff) shows this is Yoshida's integrability coefficient, so the Morales–Ramis table applies: for general μ, δ = μ(p+½) or δ = |μ(p−½)+1| (specific rational μ can add exceptional values, not included). Requiring both orbits admissible gives a discrete set of (μ, b/a) — a couple of dozen with |p|, |q| ≤ 3, positive density, |b/a| ≤ 0.7, for γ_out ∈ (0.5, 1.7). Cleanest: μ = 1, b/a = +1/2 (oblate, λ_axis = 0, zero density on the axis) and −2/5 (prolate). The formula was validated end to end by the monodromy runs (loop-at-infinity eigenvalues match, including hyperbolic cases).

The endpoint condition is necessary but far from sufficient: profiles satisfying it still gave SL(2,ℂ) (Section 2).

---

## 4. Tool 3: the axis-restricted direct method

**Construction.** Expanding {I₃, H} = 0 to second order in the transverse variables at an invariant rectilinear orbit, with I₃ = f(H₀) + I₂ + O(4), gives dI₂/dt = f'(E)·½Ẇξ², so

    Q = I₂ − f'(E)H₂ = aξ² + 2bξξ̇ + cξ̇²

is an exactly conserved quadratic form of the NVE with coefficients polynomial in ż of degrees (n, n−1, n−2) for an integral of degree n. A conserved quadratic form forces G ⊂ O(2,ℂ), the same virtually abelian structure Morales–Ramis demands, seen constructively. With D = ż∂_z − U'∂_ż the conditions Da = 2bW, Db = cW − a, Dc = −2b expand into a finite hierarchy in powers of ż: the free data are the ξ̇²-coefficients c_k(z), the top one is forced quadratic in z, and the leftover equations are constraints. (Odd-in-transverse integrals give a linear NVE solution polynomial in ż instead — a parallel hierarchy, not run.)

**n = 2 (validation).** c₀ = κ₀ + κ₂z² and

    W(κ₀+κ₂z²)² = κ₂z(κ₀+κ₂z²)U' + 2κ₀κ₂U + C.

Contains the spherical case (κ₀ = 0, C = 0), the spherical-Stäckel term A + g(θ)/r² (κ₀ = 0, C ≠ 0), and the anisotropic harmonic case. Eliminating B between axis and equator inside the P₂ family gives a linear second-order ODE on A whose exponents at z = 0 are (0, 1, 5) with the z² term forced by the inhomogeneity: a spherical centre implies a harmonic core. The Stäckel no-go is reproduced from the axis+equator shadow alone. (The α ≠ 0 branch has B ≈ −αz: a non-spherical cusp Φ − Φ₀ ∝ R²/r, consistent with cusps living on focal curves.)

**n = 3 (exact).** The ż⁰ equation is algebraic in W, W = (c₁U'' + 3U'c₁')/(4c₁), and the residual ODE for G = U' with c₁ = κ₀ + κ₂z² has the two-parameter solution

    U' = [a₄(κ₀² − 6κ₀κ₂z² + κ₂²z⁴) + a₃κ₂z(κ₂z² − κ₀)] / (κ₂²(κ₀+κ₂z²)³).

*Correction made during the session:* this family is analytic in z but allows U'(0) ≠ 0, i.e. a linear (γ = 1) cusp in r; what actually excludes it is physics: for κ₀κ₂ > 0 the numerator is negative at z² = κ₀/κ₂ for any a₃ (force reverses sign); for κ₀κ₂ < 0 the denominator vanishes at a real radius. So an even cubic integral admits no everywhere-attractive regular axis potential except Kepler and 1/r². Jaffe (U' ∝ 1/z) and cored-log profiles are outside the family (verified).

**n = 4, first order in the flattening, Hernquist cusp.** Zeroth-order seed (α+βE)(zξ̇−żξ)²; first-order unknowns B₁ and the deformation γ of c₀; c₂ allowed to shift by a quadratic (constants k₀, k₁, k₂). The elimination is degenerate in a structured way (the compatibility condition does not involve γ), and yields one third-order linear ODE for B₁. Axis and equator give the *same* ODE up to sign because at linear order both see only the −3B/z² non-spherical curvature (spherical perturbations are annihilated). Near the cusp:

    2(α−β)[z²B''' + 6zB'' + 6B'] + 6(4α−β)B + (4k₀/3)/z = 0,

indicial polynomial s(s+1)(s+2) for α ≠ β (and (s+2)(2s+3)(s+3) for α = β); the k₀ source resonates with s = 0 into a log. No solution with B₁ = o(z): every flattening compatible with an even quartic integral is singular at the centre (constant quadrupole, point quadrupole P₂/r, or the 1/r² Stäckel term). Control: a cored A gives exponents (2, 1, −2), i.e. the anisotropic-harmonic and g(θ)/r² deformations appear where they must.

**n = 4, P₂ + P₄.** Axis condition = the P₂-only operator applied to (3B₂ + 10B₄)/3; equator = minus the same operator on (3B₂ − 7.5B₄)/3 (verified symbolically). Both combinations must therefore vanish for a regular centre: B₂ = B₄ = 0. General structure: at first order each orbit sees Σ l(l+1)B_l (axis) and Σ P_l''(0)B_l (equator); with three or more multipoles these can be made to vanish identically, e.g. g(θ) = P₆ − 0.45P₄ − 5.5P₂. Such perturbations have g'' = 0 at pole and equator, so both rectilinear NVEs are exactly those of a spherical potential *to all orders in ε*. This is the blind class of every rectilinear-orbit test (Morales–Ramis, endpoint condition, axis direct method).

**n = 4 and n = 6, Jaffe cusp (Φ ∝ ln r inside).** n = 3 excluded. n = 4 first order: exponents (−1, −2, −2) for the L² seed; for β ≠ 0 the seed's ln(z/(1+z)) dominates as z → 0 and leaves the same indicial polynomial up to O(1/ln z) corrections; the source resonates with the double root. Stronger than Hernquist (no constant quadrupole allowed). n = 6: (s+1)²(s+2)³, all negative.

**n = 6, Hernquist cusp** (jet-space elimination, validated by reproducing s(s+1)(s+2) at n = 4): B₁ ODE of order 5, indicial s²(s−1)(s+1)(s+2), exponents (−2, −1, 0, 0, 1); degenerate seed gives (0, −1, −2, −3, −3/2). The s = 1 root means B₁ ∝ z: a *flattened* scale-free cusp Φ ≈ Φ₀ + r(1 + εP₂), not a spherical one. So no spherical cusp with outward flattening at degree 6 either; but this is the first degree at which a slightly flattened μ = 1 cusp passes the local test (see Section 5).

**Cored centre (cored isothermal A = ½ln(1+r²)), second order in the flattening.** Exact power series at the core. First order gives a family of B₁ regular at the centre (as it must: cored spheres have Stäckel-type flattened neighbours in general). At second order the axis and equator conditions separate (different quadratic sources); solvability conditions are polynomial in the first-order parameters. Tracking the surviving B₁-family as the truncation order M grows:
- n = 2 restriction (c₂ ≡ 0): dimension 0 at M = 6 and 9. No Stäckel-type P₂ flattening of the cored isothermal sphere.
- n = 4, β = 0 seed: dimension 3 (M = 6) → 1 (M = 9) → inconsistent at M = 12 (least-squares residual 3×10⁻² vs 10⁻⁵³ at M = 9 with b₂ normalized to 1). The M = 9 survivor was a truncation leftover.
- n = 4, β = 1/3 seed: inconsistent at M = 12 (residual 1.3×10⁻²).

---

## 5. The scale-free μ = 1 candidates: exact treatment

For Φ = r(a + bP₂), t = b/a:

- *Morales–Ramis, exact.* The NVE is hypergeometric with exponent differences (1, ½, Δ/2); Kimura's criterion reduces to Δ an odd integer. Axis: t = (9−Δ²)/(15+Δ²) ∈ {1/2, 0, −2/5, −5/8, −3/4, …}; equator: t ∈ {−2/5, 0, 1/2, 10/11, 6/5, …}. Intersection: t ∈ {0, 1/2, −2/5}. A generic small flattening does **not** survive.
- *Direct method, exact.* For a homogeneous potential the integral is quasi-homogeneous, so c_k = C_k z^{w−1−k/2} and the hierarchy is algebraic. Per orbit the admissible t at degree n are exactly the odd-Δ values with Δ ≤ n+1; the two-orbit intersection is {0, 1/2, −2/5} at n = 4, 6, 8.
- *Why these two pass.* t = 1/2 makes the axis transverse curvature (1−2t)/z vanish identically (zero on-axis density, ρ ∝ sin²θ/r); t = −2/5 makes the equatorial one (1+5t/2)/x vanish (prolate, ρ ∝ (1.6 + 1.2cos²θ)/r). Each trivializes one rectilinear NVE and lands the other on Δ = 5. They are members of the blind class; passing is not evidence of integrability.
- *Surfaces of section* (L_z = 0 meridional problem, E = 1, quick integrator): t = 0.25 (non-admissible control) shows island chains; t = 1/2 shows a strikingly regular four-lobed structure; t = −2/5 shows clean inner curves. Caveat: the quick integrator loses accuracy for near-radial orbits at the cusp (visible as scatter in the t = 0 control), so the outer regions are not diagnostic, and the single-shot Lyapunov indicator at T = 600 does not separate any case.

---

## 6. State of the polynomial route

| centre | n = 2 | n = 3 | n = 4 | n = 6 |
|---|---|---|---|---|
| Hernquist cusp, spherical, any envelope | excluded (Stäckel no-go) | excluded (exact) | excluded (local exponents; P₂ and P₂+P₄) | excluded (local exponents) |
| Jaffe cusp, spherical, any envelope | excluded | excluded | excluded (stronger) | excluded |
| cored (harmonic centre, log middle) | excluded at 2nd order | excluded (exact) | excluded at 2nd order, two seeds | not run |
| flattened scale-free μ = 1 cusp | — | — | only t ∈ {1/2, −2/5} | same |

Within the class the two rectilinear orbits can see, the polynomial route for a smooth spherical-centre-to-flattened elliptical is closed through degree 6 (cusps) and degree 4 (cores), and the closure is not fine-tuned: the obstruction appeared at O(1) in every test. The consistent message is that a cusp carrying a polynomial integral wants to be non-spherical, and a spherical cusp wants a non-polynomial integral — in line with the mirror-project structure and the atlas intuition.

---

## 7. Open items and recommended next steps

1. **t = −2/5 prolate cusp.** Well-posed candidate that nothing so far excludes. Numerical: regularized integrator + NAFF at L_z = 0 (scale-free, one energy, a few hundred orbits to ~10³ periods); a single resonance with a measurable chaotic layer kills it. Analytic: second-order variational equation (Morales–Ramis–Simó) along the axis, the orbit whose NVE is non-trivial. t = 1/2 is less interesting (zero on-axis density).
2. **The blind class.** Potentials with g'' = 0 at pole and equator need either the full meridional-plane 2D direct method (Bertrand–Darboux-type PDE system, perturbative around a spherical seed) or variational equations along non-rectilinear periodic orbits (numerical).
3. **Odd-in-transverse integrals** (linear NVE solutions polynomial in ż): parallel hierarchy, not yet run.
4. **Other cusp slopes** (0 < γ < 2 with non-integer powers) and **n = 6 second order for cores**: same code, longer runs.
5. **Transcendental / non-analytic construction.** The results above point here. A smooth non-analytic integral requires the chaotic set to have measure zero; for triaxial cusps the realistic target is a quasi-integral on the regular set plus a structural reason for exponentially thin layers (the coboundary mechanism restated).

---

## 8. Corrections made during the session

- The claim that numerical surveys "prove" cuspy triaxial potentials are chaotic was withdrawn; only generic non-integrability is a theorem.
- The homogeneous-limit argument at the centre was noted to be avoidable by a spherical centre; the analysis then moved to the transition region.
- The n = 3 statement "U must be analytic ⇒ no cusp" was corrected: the exact family allows U'(0) ≠ 0; the exclusion comes from force reversal / real singularities instead.
- The first numerical surviving n = 4 family at second order (cored case, M = 9) was identified as a truncation leftover by pushing to M = 12.

---

## 9. Code index (all attached in the conversation)

- `nve_derive.py` — symbolic derivation of both NVEs and endpoint exponents.
- `monodromy.py`, `check.py` — numerical monodromy / Galois test for the Hernquist family, robustness checks.
- `rational_monodromy.py` — general rational A, B (sympy reduction, singular points, monodromy), families F and G.
- `endpoint.py` — endpoint admissibility search (both orbits).
- `direct_axis.py` — the ż-hierarchy for any degree n and parity.
- `n4_pert.py`, `n4_elim2.py` — first-order n = 4 problem and elimination to the B₁ ODE.
- `multipole_n4.py` — P₂ + P₄ structure.
- `n6_local2.py`, `n6_jet.py`, `n6_more.py` — jet-space elimination, n = 6 exponents, Jaffe local form.
- `second_order.py`, `so_num.py` — cored second-order test, exact and numerical.
- `scalefree_exact.py` — exact quasi-homogeneous direct method for r(1 + tP₂).
- `sos.py`, `sos2.py` — surfaces of section and Lyapunov indicator; images `sos_0.25_0.5.png`, `sos_0.0_-0.4.png`.

---

## 10. Addendum: exact scale-free direct method (Hietarinta) — degrees 4 and 6

*Added after the note above. The effort was redirected to analytic methods; orbit calculations are used only to discard cases.*

**Reduction.** For V = rᵘf(θ) in the meridional plane, a polynomial integral of degree n can be taken quasi-homogeneous; its top part is a rank-n Killing tensor L^d Q_{n−d}(p) (d even by reflection symmetry) and every lower coefficient factorizes as r^{a}×(function of θ). {I,H} = 0 becomes a system of linear ODEs in θ for the lower coefficient functions, with f in the coefficients and sources linear in the Killing parameters. Regular solutions are trigonometric polynomials (integrating factors are pure powers of sin θ, cos θ; regularity kills the logarithms). Eliminating the unknown functions by jet algebra leaves conditions on f alone. Painlevé analysis is vacuous for 0 < μ < 2 (no movable singularities), so the direct method is the tool.

**Degree 4, μ = 1 (all weights).** Weight 4 (L²(a p_x² + b p_z²)): the elimination is singular only on the sphere (determinant ∝ f′ sin θ) and reduces to the master ODE

    (a sin²θ + b cos²θ)(f″ − μ(μ+1)f) = (2μ+1)(a − b) sin θ cos θ f′,

with residual conditions ∝ a·b·f′³, so a·b = 0; with b = 0 the general solution is f = A(1 + cos θ)^{μ+1} + B(1 − cos θ)^{μ+1}, the degree-μ parabolic-Stäckel family (Sridhar–Touma), whose quartic integral is the square of the quadratic one. Established explicitly at μ = 1, 1/2, 3/2 (identical structure). Weight 6 (L⁴): only the sphere. Weight 2 (constant Q₄): Groebner basis of the residual conditions gives only (1, 2, 1) = H², apart from the degenerate locus f + f″ = 0 (uniform field). The earlier scale-free "candidates" t = 1/2 and −2/5 are exactly the two parabolic members (separable in parabolic coordinates about z, resp. about x); t = −2/5 is separable only in each meridional plane and has no 3D third integral.

**Degree 6, μ = 1 (all Killing types).** d = 6: sphere. d = 4: f″ + f = K√(a₀sin²θ + a₁cos²θ) plus a second-order condition; only the sphere. d = 2: f″ + f = K·Q₄(sin θ, cos θ), so f = c₀ + c₂cos 2θ + c₄cos 4θ with the Killing tensor fixed by (c₀, c₂, c₄) (invertible map); exact linear algebra on the full θ-system over the (c₂, c₄) plane finds only the sphere (E²L²) and the parabolic profiles c₂ = ±1/3, c₄ = 0 ((L·p)²H); nothing with c₄ ≠ 0. d = 0: the level-5 condition factorizes as (f + f″)·P(sin θ) with P ≡ 0 iff Q₆ ∝ (p²)³; only H³ and the uniform field.

**Conclusion.** Through degree 6, for every angular profile, the only scale-free degree-1 planar potentials with a polynomial integral are the sphere, the uniform field and the parabolic-Stäckel family; all higher-degree integrals are products of lower ones. The mechanism (first compatibility condition fixes f + f″ from the Killing form; the rest is finite algebra) is expected to persist at degree 8. The polynomial route in the scale-free meridional plane is closed to this order; 3D (L_z ≠ 0) structure only adds constraints.

**Code (this addendum).** `sf_direct.py`, `sf_linalg.py`, `sf_scan.py` (system, fixed-t solve, Legendre scans); `sf_classify*.py`, `sf_master*.py`, `sf_w6.py`, `sf_w2c.py`, `sf_mu.py` (degree-4 classification); `sf_general.py`, `sf_fast.py`, `sf_first.py` (general degree, fast elimination, level-5 first condition); `deg6_d4_solve.py`, `deg6_close.py`, `deg6_close2.py` (degree-6 closure).

---

## 11. Addendum: first-order 2D direct method around a spherical seed (transition region)

Layer-3 tests along the rectilinear orbits (Morales–Ramis–Simó) are structurally blind to the class with g″ = 0 at pole and equator through VE₃ (VE₂ involves only Φ_xxz, fixed by the axis data; VE₃ is abelian for any Φ_xxxx once W = U′/z — verified numerically with a 17-dimensional monomial system, controls: sphere abelian to 10⁻¹³, pure P₂ non-abelian at 0.33). The efficient tool for the transition region is instead the first-order 2D direct method: seed H₀ spherical (Hernquist), seed integral (α + βH₀)L² + δL⁴, perturbation Σ_l B_l(r)P_l; {I₁, H₀} = −{I₀, H₁} decouples into Fourier modes in θ, mode m sourced by l ≥ m; per mode nine unknown functions of r, twelve equations, jet elimination over ℚ(i)(r), conditions linear in the B-jets.

Results (quartic integrals, Hernquist seed, any seed parameters): modes 6 and 4 force B_l ∝ 1/r² (the spherical-Stäckel g(θ)/r² term; consistency identically satisfied); mode 2 gives a seed-independent second-order ODE, r²(1+r)(1+3r)B″ + 2r(9r²+8r+2)B′ + 2(9r²+4r+1)B = 0, with exponents (−1, −2) at the centre (exact 1/r² solution plus a 1/r-type solution); mode 0 gives nothing. Hence no first-order quartic deformation with a regular spherical centre exists for any finite multipole content; the blind class (regular P₆ component) is excluded, and so is the constant central quadrupole left open by the axis test. Code: `fo2d.py` (build and eliminate per mode), `fo2d_m6*.py`, `fo2d_gen.py`, `fo2d_m2.py`, `fo2d_m0.py`; VE₃: `ve3*.py`.

---

## 12. Addendum: transition region — other seeds, degree 6, and second order for cores

**Other seeds, degree 4, first order (2D direct method, `fo2d2.py`).** Jaffe cusp A = ln(r/(1+r)): mode 2 gives r²(1+r)(2+3r)B″ + 2r(9r²+13r+5)B′ + 2(9r²+11r+4)B = 0 with exponents (−2, −2) at the centre (1/r² and (ln r)/r²); mode 4 forces B₄ ∝ 1/r². No regular flattening. Cored seed A = −1/(1+r²): mode 2 exponents (2, −2), regular solution exactly B₂ = r²/(1+r²)² (the linear P₂ part of −1/(1 + R² + z²/q²)); first order is permissive for cores, as expected.

**Degree 6, first order, Hernquist (sextic seed (α+βH₀+γH₀²)L² + (δ+ηH₀)L⁴ + κL⁶, sextic I₁; 100×96 jet systems per mode, monomial-decomposed multi-column solve).** Mode 6: B₆ ∝ 1/r² only, despite the rank-6 Killing tensors reaching angular mode 6. Mode 2: the same ODE and exponents (−1, −2) as at degree 4. Nothing new at degree 6.

**Second order for the cored seed (`so2d_seq.py`, `so2d_series.py`).** The first-order integral I₁ was solved explicitly by sequential quadrature: its mode-2 part is unique and rational (K_c40 = 3bβ/16 removes the logarithms; the remaining constants are fixed by the degree-3 and degree-1 compatibilities), its mode-0 part is elementary with free constants equal to polynomials in H₀, L². The O(ε²) condition {I₂, H₀} + {I₁, H₁} + {I₀, H₂} = 0, with H₂ = C₀ + C₂P₂ + C₄P₄ regular, was tested by exact power series at the core: solvable in modes 0, 2, 4 (mode 4 to N = 14), already with a quadratic seed; C₄ is unique, C₄ = b²(−(18/35)r⁴ + (54/35)r⁶ + …), which is exactly the P₄ part of −1/(1 + R² + z²/q²) expanded in the ellipticity (b = 2ε/3). The quartic seed adds no freedom (identical rank and free-direction counts in modes 4 and 2). Conclusion: the second-order continuation of the flattened core is the Kuzmin/Stäckel-type deformation with a quadratic integral, not a non-classical quartic one.

**State of the transition-region program.** Cusped spherical centres (Hernquist, Jaffe): no regular polynomial deformation at first order for any finite multipole content, at degree 4 and (Hernquist) degree 6. Cored centres: only the classical Stäckel-type family through second order at degree 4. Everything found anywhere in this program is separable.
