import numpy as np, itertools, sympy as sp, sys
from ve3 import *
def analyse2(Ms, label):
    I = np.eye(17)
    cand = []
    for M in Ms:
        cand += [M, M @ M, np.linalg.matrix_power(M, 4)]
    for M in Ms:
        Minv = np.linalg.inv(M)
        for N in list(cand): cand.append(M @ N @ Minv)
    unip = [U for U in cand if np.allclose(np.linalg.eigvals(U), 1, atol=1e-4) and np.linalg.norm(U - I) > 1e-8]
    def worst(sub):
        w = 0.0
        for P, Q in itertools.combinations(unip, 2):
            Ps, Qs = P[np.ix_(sub, sub)], Q[np.ix_(sub, sub)]
            w = max(w, np.linalg.norm(Ps @ Qs - Qs @ Ps)/(np.linalg.norm(Ps)*np.linalg.norm(Qs)))
        return w
    full = list(range(17)); low = list(range(15))
    print(f"  {label}: {len(unip)} unipotent G0-elements; max commutator full {worst(full):.2e}, VE1+VE2 block {worst(low):.2e}")
A = -1/(1+r); B = r**2/(1+r)**4
blind = {2: -sp.Rational(11, 2), 4: -sp.Rational(9, 20), 6: 1}
E = -0.5
for label, gc, ev in [('spherical', blind, 0), ('blind eps=0.05', blind, sp.Rational(1,20)), ('blind eps=0.2', blind, sp.Rational(1,5)), ('pure P2 eps=0.2', {2: 1}, sp.Rational(1,5))]:
    U, W, Y = axis_data(A, B, gc, ev); sing = singular_points(U, W, Y, E)
    Ms = monodromies(build_system(U, W, Y, E), sing)
    analyse2(Ms, label); sys.stdout.flush()
