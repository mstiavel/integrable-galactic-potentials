# 06 — Numerical checks (used only to discard cases)

* `sos.py`, `sos2.py` — surfaces of section (z = 0, p_z > 0) and a crude finite-time Lyapunov
  indicator for the scale-free meridional problem V = r(1 + t P2) at E = 1.
* `sos_0.25_0.5.png`, `sos_0.0_-0.4.png` — sections for t = 0.25 (island chains, non-integrable),
  t = 0.5 and t = -0.4 (the two parabolic-Stäckel members), t = 0 (integrable control).

Caveat: the quick integrator loses accuracy for near-radial orbits at the cusp (visible as scatter
in the t = 0 control), so only the inner regions of the sections are diagnostic. These plots played
no role in any exclusion; every exclusion in this repository is analytic.
