import numpy as np, itertools
from monodromy import *

def analyse(eps, n, E, orbit, base, rtol, atol, frac=0.4):
    U, dU, W = make(eps, n, orbit)
    tps = turning_points(eps, n, E, orbit)
    sing = [0.0+0j, -1.0+0j] + list(tps)
    # order singular points by argument around base (for the product identity)
    order = np.argsort([np.angle(s-base) for s in sing])
    Ms = {}
    for k in order:
        s = sing[k]
        others = [t for t in sing if abs(t-s) > 1e-9]
        rad = frac*min([abs(s-t) for t in others] + [abs(s-base)])
        d = (base-s)/abs(base-s); start = s + rad*d; th0 = np.angle(d)
        circ = [s + rad*np.exp(1j*(th0 + 2*np.pi*j/96)) for j in range(97)]
        pts = [base, start] + circ[1:] + [start, base]
        Ms[k] = integrate_path(U, dU, W, E, pts, rtol, atol)
    # product in angular order (one orientation) ~ loop around infinity
    P = np.eye(2, dtype=complex)
    for k in order: P = Ms[k] @ P
    S = {k: sl2(M) for k, M in Ms.items()}
    rel = 0.0
    for i, j in itertools.combinations(S, 2):
        C = S[i]@S[j] - S[j]@S[i]
        rel = max(rel, np.linalg.norm(C)/(np.linalg.norm(S[i])*np.linalg.norm(S[j])))
    # derived subgroup abelian?
    comms = [comm(S[i], S[j]) for i, j in itertools.combinations(S, 2)]
    relC = 0.0
    for C1, C2 in itertools.combinations(comms, 2):
        X = C1@C2 - C2@C1
        relC = max(relC, np.linalg.norm(X)/(np.linalg.norm(C1)*np.linalg.norm(C2)))
    maxnorm = max(np.linalg.norm(M) for M in Ms.values())
    return dict(P_eig=np.linalg.eigvals(P), P_tr=np.trace(P), rel_noncomm=rel, rel_derived_noncomm=relC, maxnorm=maxnorm, Ms=Ms)

for orbit in ['axis','eq']:
  for eps, E in [(0.0,-0.5),(0.3,-0.5),(1.0,-0.5)]:
    print(f"\n--- orbit={orbit} eps={eps} E={E}")
    for base, rtol, atol, frac in [(0.6+0.9j,1e-11,1e-13,0.4),(0.3+1.4j,1e-12,1e-14,0.3),(-0.2+0.6j,1e-12,1e-14,0.5)]:
        r = analyse(eps, 4, E, orbit, base, rtol, atol, frac)
        print(f" base={base} rtol={rtol}: loop@inf eig={np.round(r['P_eig'],6)} tr={r['P_tr']:.6f} | max|M|={r['maxnorm']:.2e} | rel noncomm gens={r['rel_noncomm']:.2e} | rel noncomm of derived subgroup={r['rel_derived_noncomm']:.2e}")
