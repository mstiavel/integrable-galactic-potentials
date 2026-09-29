# CLAUDE.md — context for continuing this project with Claude Code

## What this project is
An analytic search for non-Stäckel integrable potentials for elliptical galaxies (spherical centre,
flattened envelope). Read `docs/research_note.md` first: it is the complete record of methods,
results, corrections and open items. The README in each numbered directory explains that
directory's scripts and the order to run them.

## State of the problem (September 2026)
* Polynomial integrals: excluded through degree 6 in every setting reached (see README table).
* Meromorphic integrals: excluded by Morales–Ramis wherever the rectilinear-orbit NVE is non-trivial.
* Open: the "blind class" (angular profiles with g'' = 0 at pole and equator) for non-polynomial
  analytic integrals; the Birkhoff normal-form test was non-discriminating at order ≤ 20.
* Candidates that survived a test were all classical in disguise: t = 1/2, -2/5 scale-free
  (parabolic-Stäckel); the flattened core -1/(1+R^2+z^2/q^2) (oblate Stäckel, Delta^2 = 1-q^2).

## Conventions used throughout the code
* Meridional plane (L_z = 0 unless stated), coordinates (x, z) or (r, theta) with theta from the
  symmetry axis; momenta (p_x, p_z) or (p_r, p_theta); L = x p_z - z p_x.
* Along a rectilinear orbit the independent variable is z (or x); NVE:
  `2(E-U) xi'' - U' xi' + (U'/z - 3 sigma B/z^2) xi = 0`, U = A+B on the axis (sigma=+1),
  U = A-B/2 in the plane (sigma=-1), for Phi = A(r) + B(r) P2(cos theta).
* Trigonometric canonical form: S = sin theta, C = cos theta with C^2 -> 1 - S^2.
* Jet elimination: unknown functions are replaced by jet symbols `name_k` (k-th derivative), with
  derivative rules implemented by a `D`/`Dth`/`Dz` function in each script.
* Exact linear algebra over fraction fields via `sympy.polys.matrices.DomainMatrix`; keep the field
  generators minimal (only symbols that actually occur) — enlarging the field was the main cause
  of timeouts.
* Numerical monodromy: fundamental matrix propagated with `solve_ivp(DOP853)` around each singular
  point; Galois test = derived subgroup non-abelian + an element of infinite order.

## Practical lessons (read before running long jobs)
* Never `pkill -f` with a pattern that also matches your own shell command.
* Background jobs: launch from a script with `setsid nohup ... &`; two parallel FFT-heavy jobs
  (Birkhoff at order ≥ 24) exceed 4 GB.
* Sympy `series` of nested square roots is slow; use the truncated power-series class in
  `05_birkhoff_normal_form/bnf2.py` instead.
* The rescaling K(z) = lam^2 H(z/lam) is required for the Birkhoff series beyond order ~16 in
  double precision.

## Suggested next steps
1. VE5 along the axis for the blind class (monomial linearization as in `ve3.py`, dimension ~ 60).
2. Birkhoff normal form to order 30+ with a sparse representation and mpmath, with an
   optimal-truncation diagnostic; calibrate on `flat` vs `kk` before using on any candidate.
3. Second-order 2D direct method for other cored seeds (change `A` in `so2d_seq.py`, `so2d_series.py`).
4. Scale-free classification at degree 8 with the angular variable u = cos 2 theta.
