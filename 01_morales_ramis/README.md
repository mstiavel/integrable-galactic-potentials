# 01 — Morales–Ramis via numerical monodromy; endpoint condition

## Purpose
Test whether a given profile Phi = A(r) + B(r) P2 admits *any* meromorphic third integral, using the
Morales–Ramis theorem: the identity component of the differential Galois group of the normal
variational equation (NVE) along a particular solution must be abelian. Along the rectilinear orbits
the NVE is Fuchsian, so the Galois group is the Zariski closure of the monodromy group (Schlesinger),
which is computed numerically.

## Scripts (run order)
1. `nve_derive.py` — symbolic derivation of the NVE along the axis and the equatorial orbit; verifies
   the transverse curvature formulas against the Cartesian Hessian; local exponents at a
   power-law endpoint (the Yoshida coefficient).
2. `monodromy.py`, `check.py` — monodromy/Galois test for the Hernquist family
   A = -1/(1+r), B = eps r^2/(1+r)^4; robustness checks (base point, tolerance, loop radius,
   loop-at-infinity identity, scale-invariant commutator measures).
3. `rational_monodromy.py` — general rational A(z), B(z) as sympy expressions (exact reduction,
   singular points, monodromy); families F (Hernquist cusp -> r^1 envelope) and G (cored -> r^1).
4. `endpoint.py` — endpoint (scale-free envelope) admissibility: delta = Delta/2 with
   Delta^2 = (mu+2)^2 - 24 b_eff/u_eff must be a Morales–Ramis-admissible value on both orbits.
5. `mr_core.py` — the same test for the flattened core -1/(1+R^2+z^2/q^2) (passes: dihedral group,
   consistent with it being Stäckel).

## Results
* Hernquist + P2 flattening: G0 = SL(2,C) for every eps != 0 (already at eps = 1e-3), both orbits,
  two energies. No meromorphic integral. Spherical control: unipotent group, commutators ~1e-13.
* Families F and G: G0 = SL(2,C) at the endpoint-admissible flattenings and at generic ones.
* Endpoint condition: discrete admissible (mu, b/a); for mu = 1 both orbits admit only
  b/a in {0, 1/2, -2/5} — later identified as the parabolic-Stäckel members.
