import numpy as np, sympy as sp, itertools, sys
from rational_monodromy import singular_points, integrate_path, sl2, comm
z = sp.symbols('z'); q = sp.Rational(4, 5)
def analyse(Us, Ws, E, base=0.6+0.9j, frac=0.4, nseg=96):
    dUs = sp.diff(Us, z)
    U, dU, W = [sp.lambdify(z, e, 'numpy') for e in (Us, dUs, Ws)]
    sing = singular_points(Us, Ws, E)
    order = np.argsort([np.angle(s - base) for s in sing])
    Ms = []
    for k in order:
        s = sing[k]; others = [t for t in sing if abs(t - s) > 1e-9]
        rad = frac*min([abs(s - t) for t in others] + [abs(s - base)])
        d = (base - s)/abs(base - s); start = s + rad*d; th0 = np.angle(d)
        circ = [s + rad*np.exp(1j*(th0 + 2*np.pi*j/nseg)) for j in range(nseg + 1)]
        Ms.append(integrate_path(U, dU, W, E, [base, start] + circ[1:] + [start, base]))
    S = [sl2(M) for M in Ms]
    comms = [comm(S[i], S[j]) for i, j in itertools.combinations(range(len(S)), 2)]
    relC = max(np.linalg.norm(C1@C2 - C2@C1)/(np.linalg.norm(C1)*np.linalg.norm(C2)) for C1, C2 in itertools.combinations(comms, 2))
    hyper = max(abs(np.trace(S[i]@S[j])) for i, j in itertools.permutations(range(len(S)), 2))
    return sing, relC, hyper
for label, Us, Ws in [
    ('axis  q=0.8', -1/(1 + z**2/q**2), 2/(1 + z**2/q**2)**2),
    ('plane q=0.8', -1/(1 + z**2), 2/(q**2*(1 + z**2)**2)),
    ('axis  q=1 (spherical control)', -1/(1 + z**2), 2/(1 + z**2)**2)]:
    for E in [-0.5, -0.2]:
        sing, relC, hyper = analyse(Us, Ws, E)
        print(f"{label:32s} E={E}: singular points {len(sing)}; derived-subgroup noncommutativity {relC:.2e}; max|tr| {hyper:.2f}")
