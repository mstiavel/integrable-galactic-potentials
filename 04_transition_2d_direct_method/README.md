# 04 — First- and second-order 2D direct method around a spherical seed; VE3

## Idea
Seed H0 = 1/2(p_r^2 + p_theta^2/r^2) + A(r), seed integral (alpha + beta H0) L^2 + delta L^4
(or the sextic analogue), perturbation H1 = sum_l B_l(r) P_l(cos theta). The O(eps) condition
{I1, H0} = -{I0, H1} decouples into Fourier modes in theta; mode m is sourced only by l >= m.
Per mode: nine (quartic) or sixteen (sextic) unknown functions of r, twelve or twenty equations,
jet elimination over Q(i)(r); the conditions are linear ODEs on the B_l.

## Scripts
* `fo2d.py m l1,l2,...` — quartic, Hernquist seed, one mode (original version).
* `fo2d2.py seed n m l1,l2,...` — general: seeds `hernquist`, `jaffe` (log adjoined), `cored`;
  n = 4 or 6; rank-aware, monomial-decomposed multi-column solve (fast).
* `fo2d_m6.py`, `fo2d_m6b.py`, `fo2d_gen.py`, `fo2d_m2.py`, `fo2d_m0.py`, `sext_an.py` — analysis
  of the per-mode conditions (jet nullspace, B'/B, second-order ODE and exponents).
* Second order for the cored seed -1/(1+r^2): `so2d_seq.py <m>` (explicit first-order integral by
  sequential quadrature, modes 2 and 0), `so2d_step2.py` (second-order equation per mode, slow
  exact route), `so2d_series.py m N [quad]` (exact power-series solvability at the core; the
  practical route). `so2d_step1.py` is the failed rational-ansatz attempt kept for the record.
* `ve3_setup.py`, `ve3.py`, `ve3b.py`, `ve3c.py`, `ve3d.py` — Morales–Ramis–Simó third-order
  variational equation along the axis (17-dimensional monomial linearization); shows VE2 and VE3
  are structurally blind to the blind class.

## Results
* Hernquist, degree 4: modes 6 and 4 force B_l ~ 1/r^2 (the spherical-Stäckel g(theta)/r^2 term);
  mode 2 obeys the seed-independent ODE
  r^2(1+r)(1+3r) B'' + 2r(9r^2+8r+2) B' + 2(9r^2+4r+1) B = 0 with exponents (-1,-2); mode 0 free.
  => no regular first-order flattening for any finite multipole content; the blind class is excluded.
  Degree 6: identical (mode 6 -> 1/r^2; mode 2 the same ODE).
* Jaffe: mode 2 exponents (-2,-2); mode 4 -> 1/r^2. Cored: mode 2 exponents (2,-2), regular
  solution B2 = r^2/(1+r^2)^2 exactly.
* Cored, second order: solvable in modes 0, 2, 4 already with a quadratic seed; C4 unique,
  C4 = b^2(-(18/35) r^4 + (54/35) r^6 + ...), i.e. the P4 part of -1/(1+R^2+z^2/q^2), which is
  Stäckel in oblate spheroidal coordinates with Delta^2 = 1 - q^2 (verified exactly in 01/mr_core.py
  and by the Stäckel condition). The quartic seed adds nothing.
* VE3 along the axis: abelian for the blind class and indeed for any Phi_xxxx once W = U'/z;
  controls: sphere 1e-13, pure P2 0.33.
