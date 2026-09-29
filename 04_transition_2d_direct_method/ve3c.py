import numpy as np, itertools, sympy as sp, sys
from ve3 import *
from ve3b import analyse2
A = -1/(1+r); B = r**2/(1+r)**4
blind = {2: -sp.Rational(11, 2), 4: -sp.Rational(9, 20), 6: 1}
for E in [-0.5, -0.2]:
    for label, gc, ev in [('blind eps=0.2', blind, sp.Rational(1,5))]:
        U, W, Y = axis_data(A, B, gc, ev); sing = singular_points(U, W, Y, E); Am = build_system(U, W, Y, E)
        for tol, nseg in [(1e-11, 128), (1e-13, 256)]:
            Ms = monodromies(Am, sing, nseg=nseg, rtol=tol, atol=tol*1e-2)
            Pinf = np.eye(17, dtype=complex)
            for k in np.argsort([np.angle(s - (0.6+0.9j)) for s in sing]): Pinf = Ms[k] @ Pinf
            print(f"E={E} tol={tol}: |tr(inf)-17|={abs(np.trace(Pinf)-17):.1e}", end='; ')
            analyse2(Ms, label); sys.stdout.flush()
