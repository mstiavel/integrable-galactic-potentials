# 02 — Direct method restricted to the rectilinear orbits

## Idea
Expanding {I, H} = 0 at a rectilinear orbit shows that Q = I2 - f'(E) H2 is an exactly conserved
quadratic form of the NVE with coefficients polynomial in the velocity zdot (degrees n, n-1, n-2
for an integral of degree n). Writing out conservation gives a finite hierarchy in powers of zdot
whose leftover equations are conditions on the on-axis data (U, W) and hence on (A, B).

## Scripts (run order)
1. `direct_axis.py` — builds the hierarchy for any degree n and parity (`hierarchy(n, parity)`).
2. `n4_pert.py` -> `n4_elim2.py` — first order in the flattening at degree 4 around a spherical
   seed (Hernquist by default): elimination to a single third-order ODE for B1; exponents at the
   cusp. `multipole_n4.py` — the same with P2 + P4.
3. `n6_local2.py`, `n6_jet.py`, `n6_more.py` — degree 6 by jet-space elimination (Hernquist and
   Jaffe local forms; degenerate seeds). `n6_jet.py` is the general elimination tool.
4. `second_order.py`, `so_num.py`, `so_run.py` — second order in the flattening for the cored
   isothermal seed 1/2 ln(1+r^2), by exact power series at the core (`so_num.py` solves the
   resulting quadratic compatibility system numerically). Long runs (10–30 min).

## Results
* n = 2 reproduces the Stäckel no-go (spherical centre => harmonic core) from the two-orbit shadow.
* n = 3: exact family of on-axis potentials U' = [a4(k0^2 - 6k0k2 z^2 + k2^2 z^4) + a3 k2 z(k2 z^2 - k0)]
  /(k2^2 (k0 + k2 z^2)^3); no everywhere-attractive regular member except Kepler and 1/r^2.
* n = 4 (Hernquist): master ODE near the cusp
  2(alpha-beta)[z^2 B''' + 6 z B'' + 6 B'] + 6(4 alpha - beta) B + (4k0/3)/z = 0,
  exponents s(s+1)(s+2): no regular flattening; P2+P4 gives the same operator on
  (3B2+10B4) and (3B2-7.5B4); Jaffe exponents (-1,-2,-2).
* n = 6: exponents (-2,-1,0,0,1) Hernquist, (-1,-1,-2,-2,-2) Jaffe: no spherical cusp.
* Cored isothermal, second order: n = 2 and n = 4 excluded (the P2-only truncation cannot
  represent the P4 component that the true Stäckel continuation needs; see 04).

## Notes
The "blind class" (angular profiles with sum l(l+1) B_l = 0 on the axis and sum P_l''(0) B_l = 0 on
the equator) is invisible to every test in this directory; it is closed in `04`.
