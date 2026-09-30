# 03 — Exact scale-free direct method (Hietarinta) for V = r^mu f(theta)

## Idea
For a homogeneous potential the integral is quasi-homogeneous; its top part is a rank-n Killing
tensor L^d Q_{n-d}(p) and every lower coefficient factorizes as r^a x (function of theta).
{I, H} = 0 becomes linear ODEs in theta for the coefficient functions with f in the coefficients.
Regular solutions are trigonometric polynomials (the integrating factors are pure powers of sin,
cos and regularity kills the logarithms), so the problem is exact algebra.

## Scripts
* `sf_direct.py` — builds the theta-system for weights 2, 4, 6 (degree 4) at slope mu.
* `sf_linalg.py` — exact solve at fixed rational t for f = 1 + t P2 (finds the parabolic squares).
* `sf_scan.py` — residual scan over f = 1 + t2 P2 + t4 P4 for given weight and mu.
* `scalefree_exact.py` — quasi-homogeneous exact conditions for r(1 + t P2), degrees 4, 6, 8.
* Classification with f free (degree 4): `sf_classify.py` -> `sf_classify2.py` -> `sf_classify3.py`
  -> `sf_struct.py` -> `sf_master.py`, `sf_master2.py`, `sf_degen.py` (weight 4, mu = 1);
  `sf_w6.py <w>` (weights 6 and 2), `sf_w2c.py` (weight 2 Groebner step); `sf_mu.py <mu>`
  (self-contained weight-4 pipeline at any slope).
* Degree 6: `sf_general.py` (any degree n, Killing type d; validated on degree 4),
  `sf_fast.py` (polynomial-arithmetic elimination, rank-aware, DomainMatrix solve),
  `sf_first.py` (first compatibility condition from the top level only — the cheap decisive step),
  `deg6_d4_solve.py`, `deg6_close.py`, `deg6_close2.py` (closure at d = 4 and d = 2).

## Results
* Degree 4, weight 4, any mu: master ODE
  (a sin^2 + b cos^2)(f'' - mu(mu+1) f) = (2mu+1)(a-b) sin cos f', residuals ~ a b f'^3 => ab = 0,
  general solution f = A(1+cos)^(mu+1) + B(1-cos)^(mu+1): the parabolic-Stäckel family only.
  Weight 6: sphere only. Weight 2: H^2 only (plus the uniform-field degenerate locus).
* Degree 6 (mu = 1): d = 6 sphere; d = 4 sphere (f'' + f = K sqrt(Q2(sin,cos)) plus a second-order
  condition); d = 2 sphere and parabolic only (f'' + f = K Q4(sin,cos) confines f to degree-4
  trigonometric polynomials; exact linear algebra over (c2,c4) finds only the classical points);
  d = 0 H^3 and uniform field ((f + f'') P(sin) = 0).
* The earlier Morales–Ramis "candidates" t = 1/2, -2/5 are separable in parabolic coordinates
  about z, resp. x; t = -2/5 is separable only in each meridional plane and has no 3D third integral.

## Long runs
`sf_fast.py 6 2 1` takes ~20 min and `sf_fast.py 6 0 1` more than 3 h (measured on an Apple-silicon
laptop, sympy 1.14); `sf_first.py` gives the decisive first
condition in seconds and is preferred.
