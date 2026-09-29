# Integrable galactic potentials beyond Stäckel: an analytic search

**Summary (M. Stiavelli):**

This was an attempt to find non-classical potentials suitable to describe an elliptical galaxy and possessing three integrals of the motion. The tests were developed through interaction with Claude Fable 5.1 that developed the codes and ran the tests. My idea was to do this to be able to exclude quickly cases not leading to integrability. If we had a positive identification of a candidate I planned to take over the analysis. No matter how it is found an integral of the motion can be easily verified. Unfortunately, no non-classical integral was found and I decided to share the null result and how it was obtained to potentially save time to other researchers.

## What this repository contains

An analytic (and, where needed, numerically assisted) investigation of whether a smooth
axisymmetric galactic potential with a spherical centre — cusped (Hernquist, Jaffe) or cored —
and an outward-growing flattening can carry a third integral of motion that is not of the
classical Stäckel type. The work was carried out interactively with Claude (Anthropic) in
September 2026; the full narrative, results and caveats are in `docs/research_note.md`.

The investigation proceeds through six method families, one directory each, in the order in which
they were developed. Each directory has its own README with the purpose, the scripts, the order in
which to run them, and the results they establish.

| directory | method | what it decides |
|---|---|---|
| `01_morales_ramis` | Morales–Ramis via numerical monodromy of the normal variational equation along the rectilinear orbits; endpoint (scale-free) admissibility | excludes *all meromorphic* integrals for a given profile |
| `02_axis_direct_method` | direct method restricted to the symmetry axis / equatorial orbit (conserved quadratic forms of the NVE polynomial in the velocity) | excludes polynomial integrals of degree 2–6 at a cusp; second order in the flattening for cores |
| `03_scalefree_direct_method` | exact Hietarinta direct method for scale-free potentials `r^mu f(theta)` (ODEs in theta) | complete classification of quartic and sextic integrals for arbitrary angular profile |
| `04_transition_2d_direct_method` | first-order (and second-order) 2D direct method around a spherical seed, Fourier-mode decomposition; third-order variational equations | closes the multipole "blind class"; identifies the cored continuation as Stäckel |
| `05_birkhoff_normal_form` | Birkhoff–Lie normalization around the equatorial circular orbit, convergence diagnostic | probe of analytic non-polynomial integrability (non-discriminating at the orders reached) |
| `06_numerical_checks` | surfaces of section and Lyapunov indicators | used only to discard cases |

## Headline results

* Cusped spherical centres (Hernquist, Jaffe) admit no polynomial third integral of degree ≤ 6 for
  any outward flattening, for any finite multipole content (first order in the flattening; degree 4
  and 6 in 2D, degree 2–6 on the axis).
* Scale-free potentials `r f(theta)`: every quartic integral (all weights) and every sextic integral
  (all Killing types) belongs to the sphere, the uniform field, or the parabolic-Stäckel family;
  higher-degree integrals are products of lower ones. The weight-4 classification is governed by
  one linear ODE, `(a sin^2 + b cos^2)(f'' - mu(mu+1) f) = (2mu+1)(a-b) sin cos f'`.
* Cored centres: the unique regular flattening continues to second order and is exactly the
  Stäckel potential `-1/(1+R^2+z^2/q^2)` (oblate spheroidal, `Delta^2 = 1-q^2`).
* Morales–Ramis excludes all meromorphic integrals wherever it was applicable; the one class it
  cannot see (angular profiles flat to second order at pole and equator) is excluded only at the
  polynomial level and remains open for non-polynomial analytic integrals.

## Requirements

Python ≥ 3.10 with `sympy`, `numpy`, `scipy`, `matplotlib` (see `requirements.txt`). Most scripts run in
seconds to minutes; the few long runs are flagged in the section READMEs. No network access is needed.

## Reproducibility notes

The scripts were written incrementally during the investigation and are kept as they were run
(with the small fixes noted in each README). Several scripts write intermediate pickles that later
scripts read; the run order in each README respects those dependencies.

## Citation

See `CITATION.cff`. Please cite the Zenodo DOI of the release you used.

## License

Code: MIT (see `LICENSE`). Research note, READMEs and figures: CC BY 4.0 (see `LICENSE-docs`).
