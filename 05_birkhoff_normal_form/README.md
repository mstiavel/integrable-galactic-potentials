# 05 — Birkhoff normal form around the equatorial circular orbit

## Idea
By Ito's theorem the Birkhoff normalizing transformation converges for an analytically integrable
system with nonresonant frequencies and diverges (Gevrey-1) otherwise. Diagnostic:
(max|chi_n|)^(1/n)/n -> 0 (convergent) or -> const > 0 (divergent).

## Scripts
* `bnf.py N [cases]` — reference version with symbolic Taylor expansion (slow beyond N ~ 12).
* `bnf2.py N lam case` — fast version: truncated bivariate power-series arithmetic for the Taylor
  expansion, dense 4-index arrays for the Lie series, rescaling K(z) = lam^2 H(z/lam) (lam = 3)
  to keep coefficients O(1) across orders (required beyond order ~16 in double precision).
  Cases: `sph` (spherical Hernquist), `flat` (Hernquist + 0.3 r^2/(1+r)^4 P2), `kk`
  (Kuzmin–Kutuzov a=1, c=0.5), `core` (-1/(1+R^2+z^2/q^2), q=0.8).
* `bnf_report.py` — prints the diagnostic from the logs `bnfs_<case>.log`.

## Results (order 20)
All four cases show a falling diagnostic (0.09 -> 0.05); the provably non-integrable `flat` is
indistinguishable from the integrable controls at this order. The test is therefore
non-discriminating at order <= 20 for perturbations of this size; order 30+ with extended
precision and an optimal-truncation diagnostic would be needed. Memory: one run at a time
beyond order 22 (padded FFT arrays).
